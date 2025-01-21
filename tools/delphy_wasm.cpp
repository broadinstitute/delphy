#include <cstdint>
#include <endian.h>
#include <cstddef>
#include <memory>
#include <iostream>
#include <ranges>
#include <sstream>
#include <cstdlib>

#include <emscripten/emscripten.h>
#include <boost/iostreams/stream.hpp>

#include "version.h"
#include "run.h"
#include "sequence.h"
#include "io.h"
#include "cmdline.h"
#include "dates.h"
#include "api.h"
#include "phylo_tree_calc.h"
#include "mcc_tree.h"
#include "beasty_input.h"
#include "beasty_output.h"

// We use the delphy namespace instead of entering it so that the symbols below are
// patently placed in the global namespace, but the symbols in the delphy namespace are
// easy to use.  This is the only `using namespace` declaration of the whole project.
using namespace delphy;

// The "*_async" methods in this file take a "callback_id" that is provided from the JS
// side by glue code in delphy_api.js.  Given that callback_id, they can complete a
// calculation on the JS main thread with a single result with the following magic
// snippet:
//
//  MAIN_THREAD_ASYNC_EM_ASM({
//      delphyRunCallback($0, {result1: $1, result2: $2});
//    }, callback_id, result1, result2);
//
// The callback value can be any JS value.  We leverage some dark Emscripten magic to
// ferry the variables inside a call intact to the main thread, where they are assembled
// into the result object just before the promise completes.

// Similarly, some functions accept hooks that can be called during the execution, e.g.,
// to report progress and warnings while loading a file.  On the C++ side, the functions
// take in a `hook_id` that is provided from the JS side by glue code in delphy_api.js.
// The C++ side calls the hook with the following magic snippet:
//
//  MAIN_THREAD_ASYNC_EM_ASM({
//      delphyRunHook($0, {result1: $1, result2: $2});
//    }, hook_id, result1, result2);
//
// The hook parameters are entirely unconstrained.


// A Delphy context is a set of global variables visible across all API invocations.  Separate
// contexts are completely independent
struct Delphy_context {
  Delphy_context(uint32_t prng_seed, int num_threads)
      : bitgen_{prng_seed}, thread_pool_{num_threads} {}
  
  std::mt19937 bitgen_;
  ctpl::thread_pool thread_pool_;
  Local_arena_scope scope_;
};

// Initialize the main thread pool with N worker threads
// * prng_seed: if 0, used std::random_device to initialize prng
// * num_threads: if 0, default to largest of thread::hardware_concurrency() and 1
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_create_context(uint32_t prng_seed, int num_threads) -> Delphy_context* {
  auto real_prng_seed =
      prng_seed != 0 ? prng_seed : std::random_device{}();
  auto real_num_threads =
      num_threads != 0 ? num_threads : std::max(1, estd::as_signed(std::thread::hardware_concurrency()));
  
  return new Delphy_context(real_prng_seed, real_num_threads);
}

// Delete a context
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_delete_context(Delphy_context* ctx) -> void {
  delete ctx;
}

// Versioning
// ----------

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_get_version_string(Delphy_context& /*ctx*/) -> const char* {
  return k_delphy_version_string.c_str();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_get_build_number(Delphy_context& /*ctx*/) -> int {
  return k_delphy_build_number;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_get_commit_string(Delphy_context& /*ctx*/) -> const char* {
  return k_delphy_commit_string.c_str();;
}

// In principle, the JS side can take any usable input format like a FASTA file and set up
// an initial tree.  But it's easier to expose the C++ implementations for input formats
// that we use for the command-line version.  We don't stream file contents: if a file is
// too large to hold in memory, too bad.
//
// The parsing proceeds through three stages:
// 1. Reading FASTA file
// 2. Analyzing sequences (diffing them from consensus, mapping gaps, simplifying ambiguities)
// 3. Building the initial tree, one tip at a time
//
// The `stage_progress_hook` is called a number {1,2,3} when entering each stage
// (completion is notified by calling `callback_id`, which will usually cause an
// associated JS Promise to complete).
//
// Every time a sequence is read from the original input file, the hook
// `read_progress_hook` is called with the number of sequences read so far.
//
// Every time a sequence is analyzed, the hook `analysis_progress_hook` is called with the
// number of sequences analyzed so far and the total number of sequences.
//
// Every time a tip is added to the initial tree, the hook `initial_build_progress_hook`
// is called with the number of tips added so far and the total number of tips.
//
// Warnings throughout result in a callback to the `warning_hook`, with a sequence ID (or
// an empty string if not a sequence-specific warning), a detail code
// and a JS object with details, if any (see below).

struct Sequence_warning_hook_bridge {
  int warning_hook_id;
  auto operator()(const std::string& seq_id, const Sequence_warning& warning) -> void {
    // JS calls are synchronous!  seq_id may be destroyed as soon as this hook call ends
    auto whid = warning_hook_id;
    std::visit(estd::overloaded{
        [whid, &seq_id](const Sequence_warnings::No_valid_date&) {
          MAIN_THREAD_EM_ASM(  
              {delphyRunHook($0, UTF8ToString($1), 1, {});},
              whid, seq_id.c_str());
        },
        [whid, &seq_id](const Sequence_warnings::Ambiguity_precision_loss& w) {
          MAIN_THREAD_EM_ASM(  
              {delphyRunHook($0, UTF8ToString($1), 2, {originalState: String.fromCharCode($2), site: $3});},
              whid, seq_id.c_str(), to_char(w.original_state), w.site);
        },
        [whid, &seq_id](const Sequence_warnings::Invalid_state& w) {
          MAIN_THREAD_EM_ASM(  
              {delphyRunHook($0, UTF8ToString($1), 3, {state: String.fromCharCode($2)});},
              whid, seq_id.c_str(), w.state_letter);
        },
        [whid, &seq_id](const Sequence_warnings::Invalid_gap& w) {
          MAIN_THREAD_EM_ASM(  
              {delphyRunHook($0, UTF8ToString($1), 4, {startSite: $2, endSite: $3});},
              whid, seq_id.c_str(), w.start_site, w.end_site);
        },
        [whid, &seq_id](const Sequence_warnings::Invalid_mutation& w) {
          MAIN_THREAD_EM_ASM(  
              {delphyRunHook($0, UTF8ToString($1), 5, {from: String.fromCharCode($2), site: $3, to: String.fromCharCode($4)});},
              whid, seq_id.c_str(), to_char(w.from), w.site, to_char(w.to));
        }
      }, warning);
  }
};

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_parse_fasta_into_initial_tree_async(
    Delphy_context& ctx,
    const char* fasta_bytes,
    size_t num_fasta_bytes,
    int stage_progress_hook_id,         // (stage: int) -> void
    int read_progress_hook_id,          // (seqs_so_far: int, bytes_so_far: size_t, total_bytes: size_t) -> void
    int analysis_progress_hook_id,      // (seqs_so_far: int, total_seqs: int) -> void
    int initial_build_progress_hook_id, // (tips_so_far: int, total_tips: int) -> void
    int warning_hook_id,                // (seq_id: const char*, warning_code: int, detail: object) -> void
    //                                         (called synchronously)
    int callback_id)   // Signature: (PhyloTree*) -> void (success), or (msg: const char*) -> void (failure)
    -> void {

  auto subbitgen_seed = ctx.bitgen_();
  ctx.thread_pool_.push([fasta_bytes, num_fasta_bytes,
                         stage_progress_hook_id,
                         read_progress_hook_id,
                         analysis_progress_hook_id,
                         initial_build_progress_hook_id,
                         warning_hook_id,
                         callback_id, subbitgen_seed](int /*thread_id*/) {
    try {
      // Stage 1: Reading FASTA file
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 1);}, stage_progress_hook_id);
      auto in_fasta_is = boost::iostreams::stream<boost::iostreams::array_source>{fasta_bytes, num_fasta_bytes};
      auto in_fasta = read_fasta(
          in_fasta_is,
          [read_progress_hook_id, num_fasta_bytes](int seqs_so_far, size_t bytes_so_far) {
            MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2, $3);},
                                     read_progress_hook_id, seqs_so_far, bytes_so_far, num_fasta_bytes);
          },
          Sequence_warning_hook_bridge{warning_hook_id});

      // Stage 2: Analyse it (= transform it into a MAPLE file)
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 2);}, stage_progress_hook_id);
      auto maple_file = fasta_to_maple(
          in_fasta, 
          [analysis_progress_hook_id](int seqs_so_far, int total_seqs) {
            MAIN_THREAD_ASYNC_EM_ASM(
                {delphyRunHook($0, $1, $2);},
                analysis_progress_hook_id, seqs_so_far, total_seqs);
          },
          Sequence_warning_hook_bridge{warning_hook_id});
      
      // Stage 3: Initial tree build file
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 3);}, stage_progress_hook_id);
      auto init_random = false;  // true = random, false = UShER-like
      auto bitgen = std::mt19937{subbitgen_seed};
      auto tree = new Phylo_tree{
        build_rough_initial_tree_from_maple(
            std::move(maple_file), init_random, bitgen,
            [initial_build_progress_hook_id](int tips_so_far, int total_tips) {
              MAIN_THREAD_ASYNC_EM_ASM(
                  {delphyRunHook($0, $1, $2);},
                  initial_build_progress_hook_id, tips_so_far, total_tips);
            })
      };
      
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunCallback($0, $1);}, callback_id, tree);
    } catch (std::exception& e) {
      // Sync!  e.what() may be destroyed as soon as this hook call ends 
      MAIN_THREAD_EM_ASM({delphyFailCallback($0, UTF8ToString($1));}, callback_id, e.what());
    }
  });
}

// For MAPLE files, the parsing proceeds through two stages (the analysis stage of FASTA
// file loading is essentially a FASTA -> MAPLE conversion):
// 1. Reading MAPLE file
// 2. Building the initial tree, one tip at a time
//
// The `stage_progress_hook` is called a number {1,2,3} when entering each stage
// (completion is notified by calling `callback_id`, which will usually cause an
// associated JS Promise to complete).
//
// Every time a sequence is read from the original input file, the hook
// `read_progress_hook` is called with the number of sequences read so far.
//
// Every time a tip is added to the initial tree, the hook `initial_build_progress_hook`
// is called with the number of tips added so far and the total number of tips.
//
// Warnings throughout result in a callback to the `warning_hook`, with a sequence ID (or
// an empty string if not a sequence-specific warning), a detail code
// and a JS object with details, if any (see below).

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_parse_maple_into_initial_tree_async(
    Delphy_context& ctx,
    const char* maple_bytes,
    size_t num_maple_bytes,
    int stage_progress_hook_id,         // (stage: int) -> void
    int read_progress_hook_id,          // (seqs_so_far: int, bytes_so_far: size_t, total_bytes: size_t) -> void
    int initial_build_progress_hook_id, // (tips_so_far: int, total_tips: int) -> void
    int warning_hook_id,                // (seq_id: const char*, warning_code: int, detail: object) -> void
    //                                         (called synchronously)
    int callback_id)   // Signature: (PhyloTree*) -> void (success), or (msg: const char*) -> void (failure)
    -> void {

  auto subbitgen_seed = ctx.bitgen_();
  ctx.thread_pool_.push([maple_bytes, num_maple_bytes,
                         stage_progress_hook_id,
                         read_progress_hook_id,
                         initial_build_progress_hook_id,
                         warning_hook_id,
                         callback_id, subbitgen_seed](int /*thread_id*/) {
    try {
      // Stage 1: Reading MAPLE file
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 1);}, stage_progress_hook_id);
      auto in_maple_is = boost::iostreams::stream<boost::iostreams::array_source>{maple_bytes, num_maple_bytes};
      auto in_maple = read_maple(
          in_maple_is,
          [read_progress_hook_id, num_maple_bytes](int seqs_so_far, size_t bytes_so_far) {
            MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, $1, $2, $3);},
                                     read_progress_hook_id, seqs_so_far, bytes_so_far, num_maple_bytes);
          },
          Sequence_warning_hook_bridge{warning_hook_id});

      // Stage 2: Initial tree build file
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunHook($0, 2);}, stage_progress_hook_id);
      auto init_random = false;  // true = random, false = UShER-like
      auto bitgen = std::mt19937{subbitgen_seed};
      auto tree = new Phylo_tree{
        build_rough_initial_tree_from_maple(
            std::move(in_maple), init_random, bitgen,
            [initial_build_progress_hook_id](int tips_so_far, int total_tips) {
              MAIN_THREAD_ASYNC_EM_ASM(
                  {delphyRunHook($0, $1, $2);},
                  initial_build_progress_hook_id, tips_so_far, total_tips);
            })
      };
      
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunCallback($0, $1);}, callback_id, tree);
    } catch (std::exception& e) {
      // Sync!  e.what() may be destroyed as soon as this hook call ends 
      MAIN_THREAD_EM_ASM({delphyFailCallback($0, UTF8ToString($1));}, callback_id, e.what());
    }
  });
}


// Fb_holder
// ---------

// WASM's ability to return anything but a primitive data type is basically non-existent.  We instead have to
// pass in pointers to pre-allocated structs in linear memory where results will go.  It's up to the *caller*
// to allocate these (with `Module._malloc`) and free them (with `Module._free`)

// All flatbuffers are size-prefixed, i.e., the first 4 bytes pointed to by the buffer are a little-ending uint32
// saying how many subsequent bytes make up the actual flatbuffer.

using Fb_holder = flatbuffers::DetachedBuffer;

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_fb_holder_sizeof(Delphy_context& /*ctx*/) -> int32_t {
  return sizeof(Fb_holder);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_fb_holder_construct(Delphy_context& /*ctx*/, Fb_holder* out_fb_holder) -> void {
  new (out_fb_holder) Fb_holder{};
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_fb_holder_delete(Delphy_context& /*ctx*/, Fb_holder* fb_holder) -> void {
  delete fb_holder;
}

// Return the size of the (size-prefixed) flatbuffer held by the Fb_holder (including the 4-byte size prefix!)
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_fb_holder_get_size(Delphy_context& /*ctx*/, Fb_holder* fb_holder) -> int32_t {
  return flatbuffers::GetPrefixedSize(fb_holder->data()) + 4;  // +4 => account for the size prefix!
}

// Return a pointer to the start of the (size-prefixed) flatbuffer held by the Fb_holder
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_fb_holder_get_fb(Delphy_context& /*ctx*/, Fb_holder* fb_holder) -> uint8_t* {
  return fb_holder->data();
}


// String
// ======
//
// We sometimes want to return C++ string objects, and these methods allow us to work with them

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_string_delete(Delphy_context& /*ctx*/, std::string* str) -> void {
  delete str;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_string_data(Delphy_context& /*ctx*/, const std::string* str) -> const char* {
  return str->data();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_string_size(Delphy_context& /*ctx*/, const std::string* str) -> int32_t {
  return str->size();
}


// Phylo_tree
// ==========

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_copy(Delphy_context& /*ctx*/, const Phylo_tree& tree) -> Phylo_tree* {
  return new Phylo_tree{tree};
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_copy_from(Delphy_context& /*ctx*/, Phylo_tree& dst, const Phylo_tree& src) -> void {
  dst = src;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_from_flatbuffers(
    Delphy_context& /*ctx*/,
    const uint8_t* tree_fb,
    const uint8_t* tree_info_fb)
    -> Phylo_tree* {
  try {
    return new Phylo_tree{api_tree_and_tree_info_to_phylo_tree(tree_fb, tree_info_fb)};
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return nullptr;
  }
}

// Caller is responsible for freeing the flatbuffer (using fb_holder_delete).  Root object type: `Tree`
// On failure, clears the fb_holder
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_to_flatbuffer(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Fb_holder* fb_holder)
    -> void {
  try {
    *fb_holder = phylo_tree_to_api_tree(tree);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    *fb_holder = {};
  }
}

// Caller is responsible for freeing the flatbuffer (using fb_holder_delete).  Root object type: `TreeInfo`
// On failure, clears the fb_holder
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_to_info_flatbuffer(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Fb_holder* fb_holder)
    -> void {
  try {
    *fb_holder = phylo_tree_to_api_tree_info(tree);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    *fb_holder = {};
  }
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_delete(Delphy_context& /*ctx*/, Phylo_tree* tree) -> void {
  delete tree;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_size(Delphy_context& /*ctx*/, const Phylo_tree& tree) -> Node_index {
  return tree.size();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_root(Delphy_context& /*ctx*/, const Phylo_tree& tree) -> Node_index {
  return tree.root;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_root(Delphy_context& /*ctx*/, Phylo_tree& tree, Node_index new_root) -> void {
  tree.root = new_root;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_num_sites(Delphy_context& /*ctx*/, const Phylo_tree& tree) -> Site_index {
  return tree.num_sites();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_root_sequence(Delphy_context& /*ctx*/, const Phylo_tree& tree) -> const Real_seq_letter* {
  // TODO: Either explicitly rename to "get_ref_sequence" or ensure that mutations above root are applied to ref seq.
  // In current usage, this method is only ever called at points in time where the tree has been referenced to the
  // root sequence, so there's no difference.  But someday in the future, there will be trouble...
  return tree.ref_sequence.data();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_root_sequence(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    const Real_seq_letter* new_root_sequence_data,
    Site_index new_root_sequence_length)
    -> void {
  // TODO: Either explicitly rename to "set_ref_sequence" or ensure that mutations above root are applied to ref seq.
  // In current usage, this method is only ever called at points in time where the tree has been referenced to the
  // root sequence, so there's no difference.  But someday in the future, there will be trouble...
  tree.ref_sequence.assign(new_root_sequence_data, new_root_sequence_data + new_root_sequence_length);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_parent_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).parent;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_parent_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    Node_index new_parent)
    -> void {
  tree.at(node).parent = new_parent;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_num_children_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> int {
  return std::ssize(tree.at(node).children);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_left_child_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).left_child();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_right_child_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).right_child();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_clear_children_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node)
    -> void {
  tree.at(node).children = {};
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_children_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    Node_index new_left_child,
    Node_index new_right_child)
    -> void {
  tree.at(node).children = {new_left_child, new_right_child};
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_name_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> const char* {
  return tree.at(node).name.c_str();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_name_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    const char* new_name)  // null-terminated, please!
    -> void {
  tree.at(node).name = new_name;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_time_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> double {
  return tree.at(node).t;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_time_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    double new_time)
    -> void {
  tree.at(node).t = new_time;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_min_time_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> double {
  return tree.at(node).t_min;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_min_time_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    double new_time)
    -> void {
  tree.at(node).t_min = new_time;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_max_time_of(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    Node_index node)
    -> double {
  return tree.at(node).t_max;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_set_max_time_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    double new_time)
    -> void {
  tree.at(node).t_max = new_time;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_get_mutation_list_iterators_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node,
    Mutation_list<>::iterator* out_begin,
    Mutation_list<>::iterator* out_end)
    -> void {
  *out_begin = tree.at(node).mutations.begin();
  *out_end = tree.at(node).mutations.end();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_phylo_tree_clear_mutations_of(
    Delphy_context& /*ctx*/,
    Phylo_tree& tree,
    Node_index node)
    -> void {
  tree.at(node).mutations = {};
}

// Mutation_list
// =============

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_get_iterator_size(Delphy_context& /*ctx*/) -> size_t {
  return sizeof(Mutation_list<>::iterator);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterators_are_equal(
    Delphy_context& /*ctx*/,
    const Mutation_list<>::iterator& lhs,
    const Mutation_list<>::iterator& rhs)
    -> bool {
  return lhs == rhs;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterator_advance(
    Delphy_context& /*ctx*/,
    Mutation_list<>::iterator& it)
    -> void {
  ++it;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterator_get_from(
    Delphy_context& /*ctx*/,
    Mutation_list<>::iterator& it)
    -> Real_seq_letter {
  return it->from;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterator_get_site(
    Delphy_context& /*ctx*/,
    Mutation_list<>::iterator& it)
    -> Site_index {
  return it->site;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterator_get_to(
    Delphy_context& /*ctx*/,
    Mutation_list<>::iterator& it)
    -> Real_seq_letter {
  return it->to;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mutation_list_iterator_get_time(
    Delphy_context& /*ctx*/,
    Mutation_list<>::iterator& it)
    -> double {
  return it->t;
}


// TODO: Continue with missations & expose in JS land


// Run
// ===

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_create_run(
    Delphy_context& ctx,
    Phylo_tree& phylo_tree,
    uint32_t prng_seed)
    -> Run* {
  
  try {
    auto prng = std::mt19937{prng_seed == 0 ? absl::Uniform<std::mt19937::result_type>(ctx.bitgen_) : prng_seed};
    return new Run(ctx.thread_pool_, prng, std::move(phylo_tree));  // TODO: don't move!
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return nullptr;
  }
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_delete_run(
    Delphy_context& /*ctx*/,
    Run* run)
    -> void {
  
  delete run;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_steps_async(
    Delphy_context& ctx,
    Run& run,
    uint32_t num_steps,
    int callback_id)    // Signature () -> void (success), or (msg: const char*) -> void (failure)
    -> void {
  
  ctx.thread_pool_.push([&run, num_steps, callback_id](int /*thread_id*/) {
    try {
      run.do_mcmc_steps(num_steps);
      
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunCallback($0);}, callback_id);
    } catch (std::exception& e) {
      // Sync!  e.what() may be destroyed as soon as this hook call ends 
      MAIN_THREAD_EM_ASM({delphyFailCallback($0, UTF8ToString($1));}, callback_id, e.what());
    }
  });
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_tree(Delphy_context& /*ctx*/, Run& run) -> Phylo_tree* {
  return &run.tree();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_step(Delphy_context& /*ctx*/, const Run& run) -> int64_t {
  return run.step();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_step(Delphy_context& /*ctx*/, Run& run, int64_t num_iters) -> void {
  return run.set_step(num_iters);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_local_moves_per_global_move(Delphy_context& /*ctx*/, const Run& run) -> int32_t {
  return run.local_moves_per_global_move();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_local_moves_per_global_move(Delphy_context& /*ctx*/, Run& run, int32_t r) -> void {
  return run.set_local_moves_per_global_move(r);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_num_parts(Delphy_context& /*ctx*/, const Run& run) -> int32_t {
  return run.num_parts();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_num_parts(Delphy_context& /*ctx*/, Run& run, int num_parts) -> void {
  run.set_num_parts(num_parts);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_mu(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.mu();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_mu(Delphy_context& /*ctx*/, Run& run, double mu) -> void {
  run.set_mu(mu);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_alpha(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.alpha();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_alpha(Delphy_context& /*ctx*/, Run& run, double alpha) -> void {
  run.set_alpha(alpha);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_hky_kappa(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.hky_kappa();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_hky_kappa(Delphy_context& /*ctx*/, Run& run, double hky_kappa) -> void {
  run.set_hky_kappa(hky_kappa);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_hky_pi_A(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.hky_pi()[Real_seq_letter::A];
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_hky_pi_A(Delphy_context& /*ctx*/, Run& run, double hky_pi_A) -> void {
  auto hky_pi = run.hky_pi();
  hky_pi[Real_seq_letter::A] = hky_pi_A;
  run.set_hky_pi(hky_pi);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_hky_pi_C(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.hky_pi()[Real_seq_letter::C];
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_hky_pi_C(Delphy_context& /*ctx*/, Run& run, double hky_pi_C) -> void {
  auto hky_pi = run.hky_pi();
  hky_pi[Real_seq_letter::C] = hky_pi_C;
  run.set_hky_pi(hky_pi);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_hky_pi_G(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.hky_pi()[Real_seq_letter::G];
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_hky_pi_G(Delphy_context& /*ctx*/, Run& run, double hky_pi_G) -> void {
  auto hky_pi = run.hky_pi();
  hky_pi[Real_seq_letter::G] = hky_pi_G;
  run.set_hky_pi(hky_pi);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_hky_pi_T(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.hky_pi()[Real_seq_letter::T];
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_hky_pi_T(Delphy_context& /*ctx*/, Run& run, double hky_pi_T) -> void {
  auto hky_pi = run.hky_pi();
  hky_pi[Real_seq_letter::T] = hky_pi_T;
  run.set_hky_pi(hky_pi);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_skygrid_tau(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.skygrid_tau();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_skygrid_tau(Delphy_context& /*ctx*/, Run& run, double skygrid_tau) -> void {
  run.set_skygrid_tau(skygrid_tau);
}

// Delphy retains ownership of the returned Pop_model object
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_pop_model(Delphy_context& /*ctx*/, const Run& run) -> const Pop_model* {
  return &run.pop_model();
}

// Delphy takes ownership of the Pop_model object
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_pop_model(
    Delphy_context& /*ctx*/,
    Run& run,
    const Pop_model* pop_model)
    -> void {
  
  run.set_pop_model(std::shared_ptr<const Pop_model>{pop_model});
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_log_G(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.log_G();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_log_posterior(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.log_G() + run.log_coalescent_prior() + run.log_other_priors();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_num_mutations(Delphy_context& /*ctx*/, const Run& run) -> int {
  return run.num_muts();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_log_coalescent_prior(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.log_coalescent_prior();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_log_other_priors(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.log_other_priors();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_total_branch_length(Delphy_context& /*ctx*/, const Run& run) -> double {
  return calc_T(run.tree());
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_mpox_hack_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.mpox_hack_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_mpox_hack_enabled(Delphy_context& /*ctx*/, Run& run, bool mpox_hack_enabled) -> void {
  run.set_mpox_hack_enabled(mpox_hack_enabled);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_mpox_mu(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.mpox_mu();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_mpox_mu(Delphy_context& /*ctx*/, Run& run, double mpox_mu) -> void {
  run.set_mpox_mu(mpox_mu);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_mpox_mu_star(Delphy_context& /*ctx*/, const Run& run) -> double {
  return run.mpox_mu_star();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_mpox_mu_star(Delphy_context& /*ctx*/, Run& run, double mpox_mu_star) -> void {
  run.set_mpox_mu_star(mpox_mu_star);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_only_displacing_inner_nodes(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.only_displacing_inner_nodes();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_only_displacing_inner_nodes(Delphy_context& /*ctx*/, Run& run, bool b) -> void {
  run.set_only_displacing_inner_nodes(b);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_are_topology_moves_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.topology_moves_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_topology_moves_enabled(Delphy_context& /*ctx*/, Run& run, bool b) -> void {
  run.set_topology_moves_enabled(b);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_repartitioning_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.repartitioning_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_repartitioning_enabled(Delphy_context& /*ctx*/, Run& run, bool b) -> void {
  run.set_repartitioning_enabled(b);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_alpha_move_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.alpha_move_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_alpha_move_enabled(Delphy_context& /*ctx*/, Run& run, bool alpha_move_enabled) -> void {
  run.set_alpha_move_enabled(alpha_move_enabled);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_mu_move_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.mu_move_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_mu_move_enabled(Delphy_context& /*ctx*/, Run& run, bool mu_move_enabled) -> void {
  run.set_mu_move_enabled(mu_move_enabled);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_final_pop_size_move_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.final_pop_size_move_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_final_pop_size_move_enabled(Delphy_context& /*ctx*/, Run& run, bool final_pop_size_move_enabled) -> void {
  run.set_final_pop_size_move_enabled(final_pop_size_move_enabled);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_is_pop_growth_rate_move_enabled(Delphy_context& /*ctx*/, const Run& run) -> bool {
  return run.pop_growth_rate_move_enabled();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_pop_growth_rate_move_enabled(Delphy_context& /*ctx*/, Run& run, bool pop_growth_rate_move_enabled) -> void {
  run.set_pop_growth_rate_move_enabled(pop_growth_rate_move_enabled);
}

// Caller is responsible for freeing the flatbuffer (using fb_holder_delete).  Root object type: `Params`
// On failure, clears the fb_holder
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_get_params_to_flatbuffer(
    Delphy_context& /*ctx*/,
    const Run& run,
    Fb_holder* fb_holder)
    -> void {
  try {
    *fb_holder = run_to_api_params(run);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    *fb_holder = {};
  }
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_set_params_from_flatbuffer(Delphy_context& /*ctx*/, Run& run, const uint8_t* params_fb) -> void {
  try {
    apply_api_params_to_run(params_fb, run);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_run_export_beast_input(
    Delphy_context /*ctx*/,
    const Run& run)
    -> std::string* {
  auto ss = std::ostringstream{};
  export_beast_input(run, ss);
  return new std::string{std::move(ss).str()};  // Idiom to move string out of ss instead of copying it
}


// Mcc_tree
// ========

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_derive_mcc_tree_async(
    Delphy_context& ctx,
    Phylo_tree** trees,
    size_t num_trees,
    int callback_id)   // Signature: (MccTree*) -> void (success), or (msg: const char*) -> void (failure)
    -> void {

  auto seed = absl::Uniform<std::mt19937::result_type>(ctx.bitgen_);
  ctx.thread_pool_.push([seed, trees, num_trees, callback_id](int /*thread_id*/) {
    try {
      auto bitgen = std::mt19937{seed};
      auto tree = new Mcc_tree{derive_mcc_tree(std::vector<Phylo_tree*>{trees, trees + num_trees}, bitgen)};
      
      MAIN_THREAD_ASYNC_EM_ASM({delphyRunCallback($0, $1);}, callback_id, tree);
    } catch (std::exception& e) {
      // Sync!  e.what() may be destroyed as soon as this hook call ends 
      MAIN_THREAD_EM_ASM({delphyFailCallback($0, UTF8ToString($1));}, callback_id, e.what());
    }
  });
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_delete(Delphy_context& /*ctx*/, Mcc_tree* tree) -> void {
  delete tree;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_num_base_trees(Delphy_context& /*ctx*/, const Mcc_tree& tree) -> Base_tree_index {
  return std::ssize(tree.base_trees());
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_master_base_tree_index(Delphy_context& /*ctx*/, const Mcc_tree& tree) -> Base_tree_index {
  return tree.master_base_tree_index();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_base_tree(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Base_tree_index base_tree_index)
    -> const Phylo_tree* {
  return tree.base_trees().at(base_tree_index);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_size(Delphy_context& /*ctx*/, const Mcc_tree& tree) -> Node_index {
  return tree.size();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_root(Delphy_context& /*ctx*/, const Mcc_tree& tree) -> Node_index {
  return tree.root;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_parent_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).parent;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_num_children_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> int {
  return std::ssize(tree.at(node).children);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_left_child_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).left_child();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_right_child_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> Node_index {
  return tree.at(node).right_child();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_time_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> double {
  return tree.at(node).t();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_mrca_time_of(
    Delphy_context& /*ctx*/,
    const Mcc_tree& tree,
    Node_index node)
    -> double {
  return tree.at(node).t_mrca();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_get_corresponding_node_in_base_tree(
    Delphy_context& /*ctx*/,
    const Mcc_tree& mcc_tree,
    Node_index mcc_node,
    Base_tree_index base_tree_index)
    -> Node_index {
  return corresponding_node_to(mcc_tree, mcc_node, base_tree_index);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_is_exact_match_in_base_tree(  // i.e., are tips below MCC node monophyletic in base tree?
    Delphy_context& /*ctx*/,
    const Mcc_tree& mcc_tree,
    Node_index mcc_node,
    Base_tree_index base_tree_index)
    -> Node_index {
  return mcc_tree.at(mcc_node).corresponding_node_infos().at(base_tree_index).is_monophyletic_in_base_tree;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_mcc_tree_export(
    Delphy_context /*ctx*/,
    const Mcc_tree& tree,
    bool inner_nodes_defined_as_mrcas_of_tips)
    -> std::string* {
  return new std::string{api_export_mcc_tree(tree, inner_nodes_defined_as_mrcas_of_tips)};
}


// Population model
// ================
//
// On the C++ side, `Pop_model` is an abstract base class with different concrete implementations.
// This creates some unavoidable awkwardness in the WASM interface, which we try to limit to
// the following section

// Hard-coded IDs used only at the WASM/JS boundary (see InternalPopModelClass in delphy_api.ts)
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_pop_model_get_class(Delphy_context& /*ctx*/, Pop_model* pop_model) -> int32_t {
  if (typeid(*pop_model) == typeid(Exp_pop_model)) {
    return 1;
  } else if (typeid(*pop_model) == typeid(Skygrid_pop_model)) {
    return 2;
  } else {
    std::cerr << "Unknown Pop_model class: " << typeid(*pop_model).name() << std::endl;
    return -1;
  }
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_pop_model_delete(Delphy_context& /*ctx*/, Pop_model* pop_model) -> void {
  delete pop_model;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_new(
    Delphy_context& /*ctx*/, 
    double t0,
    double n0,
    double g)
    -> Pop_model* {
  return new Exp_pop_model{t0, n0, g};
}

// Assumes that exp_pop_model really is an `Exp_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_get_t0(Delphy_context& /*ctx*/, const Exp_pop_model& exp_pop_model) -> double {
  return exp_pop_model.t0();
}

// Assumes that exp_pop_model really is an `Exp_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_get_n0(Delphy_context& /*ctx*/, const Exp_pop_model& exp_pop_model) -> double {
  return exp_pop_model.pop_at_t0();
}

// Assumes that exp_pop_model really is an `Exp_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_exp_pop_model_get_g(Delphy_context& /*ctx*/, const Exp_pop_model& exp_pop_model) -> double {
  return exp_pop_model.growth_rate();
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_skygrid_pop_model_new(
    Delphy_context& /*ctx*/,
    int32_t num_knots,
    double* in_x,      // array of `num_knots` values
    double* in_gamma,  // array of `num_knots` values
    Skygrid_pop_model::Type type)
    -> Pop_model* {
  return new Skygrid_pop_model{
    std::vector<double>{in_x, in_x + num_knots},
    std::vector<double>{in_gamma, in_gamma + num_knots},
    type};
}

// Assumes that skygrid_pop_model really is an `Skygrid_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_skygrid_pop_model_get_type(
    Delphy_context& /*ctx*/,
    const Skygrid_pop_model& skygrid_pop_model)
    -> Skygrid_pop_model::Type {
  
  return skygrid_pop_model.type();
}

// Assumes that skygrid_pop_model really is an `Skygrid_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_skygrid_pop_model_get_num_knots(
    Delphy_context& /*ctx*/,
    const Skygrid_pop_model& skygrid_pop_model)
    -> int32_t {
  
  return std::ssize(skygrid_pop_model.x());
}

// Assumes that skygrid_pop_model really is an `Skygrid_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_skygrid_pop_model_get_x(
    Delphy_context& /*ctx*/,
    const Skygrid_pop_model& skygrid_pop_model,
    int32_t num_knots,
    double* out_x)
    -> void {

  if (std::ssize(skygrid_pop_model.x()) != num_knots) {
    std::cerr
        << absl::StreamFormat(
            "ERROR in delphy_skygrid_pop_model_get_x! "
            "num_knots = %d differs from actual number of knots (%d)",
            num_knots, std::ssize(skygrid_pop_model.x()))
        << std::endl;
    return;
  }

  for (auto k = 0; k != num_knots; ++k) {
    out_x[k] = skygrid_pop_model.x(k);
  }
}

// Assumes that skygrid_pop_model really is an `Skygrid_pop_model` !
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_skygrid_pop_model_get_gamma(
    Delphy_context& /*ctx*/,
    const Skygrid_pop_model& skygrid_pop_model,
    int32_t num_knots,
    double* out_gamma)
    -> void {

  if (std::ssize(skygrid_pop_model.gamma()) != num_knots) {
    std::cerr
        << absl::StreamFormat(
            "ERROR in delphy_skygrid_pop_model_get_gamma! "
            "num_knots = %d differs from actual number of knots (%d)",
            num_knots, std::ssize(skygrid_pop_model.gamma()))
        << std::endl;
    return;
  }

  for (auto k = 0; k != num_knots; ++k) {
    out_gamma[k] = skygrid_pop_model.gamma(k);
  }
}



// For the given population model parameters, renders the population curve at discrete intervals.
//
// out_values should point to an array of num_t_cells doubles.
// The entry out_values[i] is the average size of the population over the period [t_{i-1}, t_i], where
//
//   t_i = t_start + (i+1) * t_step,
//
// and t_step := (t_end - t_start) / num_t_cells.
//
// IMPORTANT: the size of out_values should be sizeof(double) * num_t_cells.
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_pop_model_render_population_curve(
    Delphy_context& /*ctx*/,
    const Pop_model& pop_model,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  using namespace delphy;
  try {
    api_render_population_curve(pop_model, t_start, t_end, num_t_cells, out_values);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}

// For a given base tree and associated population model parameters, simulates where a random
// probe sample at time t would coalesce into the tree and what its state at a given site would thus be.
// The time t is sampled uniformly at num_t_cells values covering (t_start, t_end] (note that the first
// sample is some time after t_start).
//
// There's a bug in the calculation for now that requires t_start to be earlier than the time of the root node.
//
// out_values should point to an array of 4 * num_t_cells doubles.  The entry out_values[s * num_t_cells + i]
// corresponds to the probability of observing state s (A = 0, C = 1, G = 2, T = 3) at time t = t_i, where
//
//   t_i = t_start + (i+1) * t_step,
//
// and t_step := (t_end - t_start) / num_t_cells.
//
// As everywhere else in the WASM interface, site is *0-based* here.
//
// IMPORTANT: the size of out_values should be sizeof(double) * 4 * num_t_cells.
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_pop_model_probe_site_states_on_tree(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    int32_t site,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  using namespace delphy;
  try {
    api_probe_site_states_on_tree(tree, pop_model, site, t_start, t_end, num_t_cells, out_values);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}

// For a given base tree and associated population model parameters, simulates where a random
// probe sample at time t would coalesce into the tree and which among a set of marked nodes
// would thus be its closest ancestor.
// The time t is sampled uniformly at num_t_cells values covering (t_start, t_end] (note that the first
// sample is some time after t_start).
//
// There's a bug in the calculation for now that requires t_start to be earlier than the time of the root node.
//
// marked_ancestors should point to an array of num_marked_ancestors int32_t's.  Each element is
// a node in the given tree that should be considered "marked".  An node can be k_no_node (-1)
// to indicate a node that doesn't exist on this tree (e.g., if iterating over all base trees in an MCC,
// and some node in the MCC isn't linked to any node in some base tree).
//
// out_values should point to an array of (num_marked_ancestors + 1) * num_t_cells doubles.
// The entry out_values[j * num_t_cells + i] corresponds to the probability of observing ancestor j
// as the closest marked ancestor of the probe at time t = t_i, where
//
//   t_i = t_start + (i+1) * t_step,
//
// and t_step := (t_end - t_start) / num_t_cells.  Note that j == num_marked_ancestors means "none".
//
// IMPORTANT: the size of out_values should be sizeof(double) * num_t_cells * (num_marked_ancestors + 1).
EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_pop_model_probe_ancestors_on_tree(
    Delphy_context& /*ctx*/,
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    const int32_t* marked_ancestors,
    int32_t num_marked_ancestors,
    double t_start,
    double t_end,
    int32_t num_t_cells,
    double* out_values)
    -> void {
  
  using namespace delphy;
  try {
    api_probe_ancestors_on_tree(tree, pop_model, marked_ancestors, num_marked_ancestors,
                                t_start, t_end, num_t_cells, out_values);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}


// Beasty_output
// =============

struct Beasty_output {
  std::ostringstream* log_os;
  std::ostringstream* tree_os;
  Beasty_log_output log_output;
  Beasty_trees_output trees_output;
};

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_create_beasty_output(Delphy_context& /*ctx*/, Run& run) -> Beasty_output* {
  auto log_os = new std::ostringstream{};
  auto tree_os = new std::ostringstream{};
  auto bout = new Beasty_output{
    .log_os = log_os,
    .tree_os = tree_os,
    .log_output = {log_os, true},
    .trees_output = {tree_os, true}
  };
  bout->log_output.output_headers(run);
  bout->trees_output.output_headers(run);
  return bout;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_delete_beasty_output(Delphy_context& /*ctx*/, Beasty_output* bout) -> void {
  delete bout;
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_beasty_output_snapshot(Delphy_context& /*ctx*/, Beasty_output* bout, Run& run) -> void {
  bout->log_output.output_log(run);
  bout->trees_output.output_tree(run);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_beasty_output_finalize(Delphy_context& /*ctx*/, Beasty_output* bout, Run& run) -> void {
  bout->log_output.output_footers(run);
  bout->trees_output.output_footers(run);
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_beasty_output_extract_log(Delphy_context& /*ctx*/, Beasty_output* bout) -> std::string* {
  return new std::string(std::move(*bout->log_os).str()); // Idiom to move-construct result instead of copying
}

EMSCRIPTEN_KEEPALIVE
extern "C"
auto delphy_beasty_output_extract_trees(Delphy_context& /*ctx*/, Beasty_output* bout) -> std::string* {
  return new std::string(std::move(*bout->tree_os).str()); // Idiom to move-construct result instead of copying
}
