#ifndef DELPHY_TREE_EDITING_H
#define DELPHY_TREE_EDITING_H

#include "evo_model.h"
#include "phylo_tree.h"
#include "phylo_tree_calc.h"
#include "site_deltas.h"

namespace delphy {

// When we start a tree-editing sessions to move a node X somewhere else, we
// remove all mutations above X and only keep track of the current delta between X's parent P and X.
// When the editing session ends, we synthesize mutations in the middle of the P-X branch that
// are compatible with that delta.
//
struct Tree_editing_session {
  Tree_editing_session(
      Phylo_tree& tree,
      Node_index X,
      const Global_evo_model& evo,
      Node_vector<double>& lambda_i,
      const std::vector<double>& ref_cum_Q_l,
      Node_vector<int>& num_sites_missing_at_every_node,
      Scratch_space& scratch);

  bool ongoing;
  
  Phylo_tree *tree;
  Node_index X;
  const Global_evo_model* evo;
  Node_vector<double>* lambda_i;
  const std::vector<double>* ref_cum_Q_l;
  Node_vector<int>* num_sites_missing_at_every_node;
  Scratch_space* scratch;

  Site_deltas deltas_nexus_to_X;

  // Slides the parent P of a node X to another point along the branch between P's parent G
  // and X's sibling S.
  // Only the time & sequence of P change in this elementary move (but not its missing sites):
  // the time, sequences and missing sites of G, X and S remain the same.
  //
  // Example:
  //
  //                                          +C3N --(C2G)--- X
  //                                          |
  //    G -- A4C ----- A1G ---- T2C --- A3C --+G1N -- C4T --- S
  //                                          P
  //                                     
  // becomes
  //
  //                +A3N ------------(A1G,T2G)--------------- X
  //                |
  //    G -- A4C ---+A1N------- T2C --- A3C --------- C4T --- S
  //                P
  //
  auto slide_P_along_branch(double new_t_P) -> void;

  // Given that P is immediately below its own parent G (no intervening mutations),
  // place P immediately *above* G, keeping X as a child and relinquishing the other to G.
  // If G is the root, then P becomes the new root.
  // Only the missing sites of P change in this elementary move:
  // the time, sequences and missing sites of G, U, X and S remain unchanged.
  //
  // Example:
  //
  //                +A1N---------------- X
  //                |
  //            +G3N+C2N-------- A4C --- S
  //            |   P
  //   GG ------+ G
  //            |
  //            +C2N---- U
  //
  // becomes
  //
  //        +A1N,G3N-------------------- X
  //        |
  //        |   +G3N------------ A4C --- S
  //        |   |
  //   GG --+C2N+ G
  //        P   |
  //            +------- U
  //
  auto hop_up() -> void;

  // Assume that X's parent P is immediately below its own parent G (no intervening mutations).
  // Conceptually, the P--X is attached to a branch G-S.
  // flip() flips it to the G-U branch.
  // Only the missing sites of P change in this elementary move:
  // the time, sequences and missing sites of G, U, S and X remain unchanged.
  //
  // Example:
  //
  //          P
  //      +G3N+A1N---------------- S
  //      |   |
  //      |   |
  //    --+ G +C2N-------- A4C --- X
  //      |
  //      |
  //      +C2N---- U
  //
  // becomes
  //
  //      +A1N,G3N---------------- S
  //      |
  //      |
  //    --+ G +G3N-------- A4C --- X
  //      |   |
  //      |   |
  //      +C2N+------- U
  //          P
  //
  auto flip() -> void;

  // Inverse operation of hop_up: see before/after picture above.
  // Note that hop_down(S') === do_hop_up(U'):
  //
  // Example:
  //
  //             +A1N---------------- U'
  //         P   |
  //    G ---+G3N+ U
  //         |   |
  //         |   +C2N-------- A4C --- S'
  //         |
  //         +C2N-------- X
  //
  // becomes
  //
  //             +A1N,G3N------------ U'
  //             |
  //    G -------+ U
  //             |   P
  //             +C2N+G3N---- A4C --- S'
  //                 |
  //                 +--- X
  //
  auto hop_down(Node_index SS) -> void;

  auto end() -> void;

 private:
  auto slide_root(double new_t_root) -> void;
  auto do_hop_up(Node_index X) -> void;  // Hop for an arbitrary X, not just this->X; common core of hop_up and hop_down
};

}  // namespace delphy

#endif // DELPHY_TREE_EDITING_H
