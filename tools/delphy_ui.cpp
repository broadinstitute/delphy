#include <chrono>
#include <iostream>
#include <GL/glut.h>
#include <GL/freeglut.h>  // For glutMainLoopEvent()

#include "absl/log/initialize.h"

#include "mcc_tree.h"
#include "incremental_mcc_tree.h"
#include "io.h"
#include "phylo_tree_calc.h"
#include "beasty_input.h"
#include "beasty_output.h"
#include "dates.h"
#include "cmdline.h"
#include "api.h"
#include "run.h"

namespace delphy {

static Run *ui_run = nullptr;
static Node_vector<Node_index> leaf_index{};  // 0 for inner nodes, 1..N for original numbering of tree
static int64_t steps_per_refresh = 0;
static int mcc_refresh_freq = 1;

struct Timestamp {
  int64_t step;
  std::chrono::time_point<std::chrono::high_resolution_clock> time;
};
static std::deque<Timestamp> timestamps;
inline constexpr size_t k_max_timestamps = 10;

//static std::unique_ptr<Incremental_mcc_tree> incremental_mcc_tree;
static std::deque<std::shared_ptr<Phylo_tree>> tree_snapshots;
static std::vector<std::shared_ptr<Phylo_tree>> trees_in_mcc;
static size_t max_tree_snapshots = 20;
static Mcc_tree old_mcc_tree;
static Mcc_tree mcc_tree;
static bool mcc_nodes_are_mrcas = false;
static bool ladderize_tree = true;

static Processed_cmd_line *cmd = nullptr;
static Beasty_log_output *ui_log_output = nullptr;
static Beasty_trees_output *ui_trees_output = nullptr;

static double t_min_x_axis = std::numeric_limits<double>::max();
static double t_max_x_axis = -std::numeric_limits<double>::max();
static double t_min_scale_bar, t_max_scale_bar;
static auto rescale_axes() {
  const auto& tree = ui_run->tree();
  t_min_scale_bar = std::numeric_limits<double>::max();
  t_max_scale_bar = -std::numeric_limits<double>::max();
  
  for (const auto& node : index_order_traversal(tree)) {
    auto t = tree.at(node).t;
    t_min_x_axis = std::min(t, t_min_x_axis);
    t_max_x_axis = std::max(t, t_max_x_axis);
    if (tree.at(node).is_tip()) {
      t_min_scale_bar = std::min(t, t_min_scale_bar);
      t_max_scale_bar = std::max(t, t_max_scale_bar);
    }
  }
}

template<typename T>
static auto canonicalize_tree(T& tree) -> void {
  struct Accum {
    double min = 0.0;
    double sum = 0.0;
    double count = 0.0;
    double avg() const { return sum / count; }
  };

  // First, calculate for each node the average index of all its tips
  auto canonicalize_result = Node_vector<Accum>(std::ssize(tree), Accum{});
  for (const auto& node : post_order_traversal(tree)) { // postorder to accumulate from tips to root
    if (tree.at(node).is_tip()) {
      canonicalize_result[node] =
          Accum{static_cast<double>(node), static_cast<double>(node), 1.0};
    } else {
      auto& this_result = canonicalize_result[node];
      this_result = {std::numeric_limits<double>::max(), 0.0, 0.0};
      for (const auto& child : tree.at(node).children) {
        this_result.min = std::min(this_result.min, canonicalize_result[child].min);
        this_result.sum += canonicalize_result[child].sum;
        this_result.count += canonicalize_result[child].count;
      }
    }
  }

  // Now sort children of each inner node to appear in increasing avg index
  for (const auto& node : pre_order_traversal(tree)) {  // pre-order because we mutate child indices!
    if (tree.at(node).is_inner_node()) {
      auto swap_cond = false;
      if (ladderize_tree) {
        swap_cond = canonicalize_result[tree.at(node).left_child()].count >
            canonicalize_result[tree.at(node).right_child()].count;
      } else {
        swap_cond = canonicalize_result[tree.at(node).left_child()].avg() >
            canonicalize_result[tree.at(node).right_child()].avg();
      }
      if (swap_cond) {
        tree.at(node).children = {tree.at(node).right_child(), tree.at(node).left_child()};
      }
    }
  }
}

template<typename T>
static auto position_tree_nodes(const T& tree) -> Node_vector<double> {
  // We want the height of the tree to be 0.9 of the height of the window,
  // and each tip is separated from the next by `space_between_tips`
  auto num_tips = (std::ssize(tree) + 1) / 2;
  auto space_between_tips = 0.9 / (num_tips - 1);

  auto result = Node_vector<double>(std::ssize(tree), 0.0);
  auto next_y = 0.05;
  for (const auto& node : post_order_traversal(tree)) {  // postorder to accumulate from tips to root
    if (tree.at(node).is_tip()) {
      result[node] = next_y;
      next_y += space_between_tips;
    } else {
      result[node] =
          0.5 * (result[tree.at(node).left_child()] + result[tree.at(node).right_child()]);
    }
  }
  return result;
}

static auto x_for(double t) -> double {
  auto x = 0.05 + 0.90 * (t - t_min_x_axis) / (t_max_x_axis - t_min_x_axis);
  return x;
}

static auto render_tree() -> void {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // white background
  glClear(GL_COLOR_BUFFER_BIT);

  auto& tree = ui_run->tree();
  
  // Define coordinates so (0,0) is bottom-left and (1,1) is top-right
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 2.0, 0.0, 1.0);

  glColor3d(0.0, 0.0, 0.0);  // black lines

  rescale_axes();
  
  glBegin(GL_LINES);
  glVertex2d(x_for(t_min_scale_bar), 0.03);
  glVertex2d(x_for(t_max_scale_bar), 0.03);
  glEnd();

  canonicalize_tree(tree);
  auto node_ys = position_tree_nodes(tree);

  // Draw current tree
  static auto partition_colors = std::vector<std::array<double, 3>>{
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}
  };
  for (const auto& node : index_order_traversal(tree)) {

    auto partition_index = ui_run->partition().tree_index_to_partition_index()[node];
    auto color_index = partition_index % std::ssize(partition_colors);
    auto [r,g,b] = partition_colors[color_index];
    glColor3d(r, g, b);

    glBegin(GL_LINES);
    if (node != tree.root) {
      glVertex2d(x_for(tree.at(node).t), node_ys[node]);
      glVertex2d(x_for(tree.at_parent_of(node).t), node_ys[node]);
    }
    if (tree.at(node).is_inner_node()) {
      glVertex2d(x_for(tree.at(node).t), node_ys[tree.at(node).left_child()]);
      glVertex2d(x_for(tree.at(node).t), node_ys[tree.at(node).right_child()]);
    }
    glEnd();
  }

  // Draw tips
  auto max_leaf_index = std::ranges::max(leaf_index);
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      auto half_side = 0.0025;
      auto cx = x_for(tree.at(node).t);
      auto cy = node_ys[node];

      auto N = max_leaf_index;
      auto c = static_cast<double>(leaf_index[node]) / N;
      glColor3d(0.0, 1.0 * (1 - c), 1.0 * c);

      glBegin(GL_QUADS);
      glVertex2d(cx - half_side, cy - half_side);
      glVertex2d(cx - half_side, cy + half_side);
      glVertex2d(cx + half_side, cy + half_side);
      glVertex2d(cx + half_side, cy - half_side);
      glEnd();
    }
  }

  // Draw mutations
  const auto& evo = ui_run->evo();
  for (const auto& node : index_order_traversal(tree)) {
    for (const auto& m : tree.at(node).mutations) {
      auto cx = x_for(m.t);
      auto cy = node_ys[node];
      auto beta = evo.partition_for_site[m.site];
      auto mut_radius = 0.001;
      if (beta == 0) {
        glColor3d(1.0, 0.0, 0.0);  // mutations in red
      } else {
        glColor3d(0.0, 0.0, 1.0);  // mutations in partitions other than 0 in blue
        mut_radius *= 5;
      }
      glBegin(GL_TRIANGLE_FAN);
      glVertex2d(cx, cy);
      for (auto i = 0; i <= 10; ++i) {
        glVertex2d(
            cx + mut_radius * std::cos(2 * M_PI * i / 10.0),
            cy + mut_radius * std::sin(2 * M_PI * i / 10.0)
        );
      }
      glEnd();
    }
  }

  // Draw missations
  auto missation_radius = 0.002;
  glColor3d(0.0, 1.0, 0.0);  // missations in green
  for (const auto& node : index_order_traversal(tree)) {
    if (not tree.at(node).missations.empty()) {
      auto t = (node == tree.root ? tree.at(node).t : tree.at_parent_of(node).t);
      auto cx = x_for(t);
      auto cy = node_ys[node];
      glBegin(GL_TRIANGLE_FAN);
      glVertex2d(cx, cy);
      for (auto i = 0; i <= 10; ++i) {
        glVertex2d(
            cx + missation_radius * std::cos(2 * M_PI * i / 10.0),
            cy + missation_radius * std::sin(2 * M_PI * i / 10.0)
        );
      }
      glEnd();
    }
  }

  // // Draw num_active_parts
  // auto y_min = 0.02;
  // auto y_scale = 0.02 / std::ssize(ui_run->partition().parts());

  // glBegin(GL_LINES);
  // glColor3d(0.0, 0.0, 0.0);
  // glVertex2d(x_for(ui_run->coalescent_prior_.t_ref_), y_min);
  // glVertex2d(x_for(ui_run->coalescent_prior_.t_ref_), 0.05);
  // glEnd();
  
  // for (auto p = 0; p != std::ssize(ui_run->partition().parts()); ++p) {
  //   const auto& vscp = ui_run->coalescent_prior_parts_.at(p);
  //   auto num_cells = std::ssize(vscp.num_active_parts_);
    
  //   auto color_index = p % std::ssize(partition_colors);
  //   auto [r,g,b] = partition_colors[color_index];
  //   glColor3d(r, g, b);

  //   glBegin(GL_LINE_STRIP);
  //   glVertex2d(x_for(vscp.t_ref_), y_min);
    
  //   for (auto cell = 0; cell != num_cells; ++cell) {
  //     auto t_max = vscp.t_ref_ - cell*vscp.t_step_;
  //     auto t_min = t_max - vscp.t_step_;

  //     auto y = y_min + y_scale*vscp.num_active_parts_[cell];
  //     glVertex2d(x_for(t_max), y);
  //     glVertex2d(x_for(t_min), y);
  //   }

  //   glVertex2d(x_for(vscp.t_ref_ - num_cells*vscp.t_step_), y_min);
  //   glEnd();
  // }
}


static auto render_mcc() -> void {
  // Define coordinates so (0,0) is bottom-left and (1,1) is top-right
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1.0, 1.0, 0.0, 1.0);

  canonicalize_tree(mcc_tree);
  auto node_ys = position_tree_nodes(mcc_tree);

  if (std::ssize(mcc_tree) > 0) {
    // Calculate tMRCA distribution of root
    auto root = mcc_tree.root;
    auto tmrcas = std::vector<double>{};
    for (auto i = 0; i != mcc_tree.num_base_trees(); ++i) {
      const auto& base_tree = *mcc_tree.base_trees()[i];
      auto corresponding_node = corresponding_node_to(mcc_tree, root, i);
      tmrcas.push_back(base_tree.at(corresponding_node).t);
    }
    auto sum_t = 0.0;
    auto sum_t_2 = 0.0;
    for (auto t : tmrcas) {
      sum_t += t;
      sum_t_2 += t * t;
    }
    auto mean_t = sum_t / tmrcas.size();
    auto mean_t2 = sum_t_2 / tmrcas.size();
    auto var_t = mean_t2 - mean_t * mean_t;
    auto std_t = std::sqrt(var_t);

    auto bin_width = std::max(0.01, 3 * std_t / std::sqrt(tmrcas.size()));

    auto num_bins = 40;
    auto min_bin = mean_t - (num_bins / 2) * bin_width;
    auto bin_for = [min_bin, bin_width](double t) {
      return std::floor((t - min_bin) / bin_width);
    };

    auto bin_freqs = std::vector<double>(num_bins, 0.0);
    for (auto t : tmrcas) {
      auto bin = bin_for(t);
      if (0 <= bin && bin < num_bins) {
        bin_freqs[bin] += 1.0 / tmrcas.size();
      }
    }

    for (auto bin = 0; bin < num_bins; ++bin) {
      auto bin_left = x_for(min_bin + bin * bin_width);
      auto bin_right = x_for(min_bin + (bin+1) * bin_width);

      glColor3d(0.5, 0.5, 0.5);
      glBegin(GL_QUADS);
      glVertex2d(bin_left, 0.04);
      glVertex2d(bin_left, 0.04 + 3 * bin_freqs[bin] / bin_width);
      glVertex2d(bin_right, 0.04 + 3 * bin_freqs[bin] / bin_width);
      glVertex2d(bin_right, 0.04);
      glEnd();
    }
  }

  // Axis
  glColor3d(0.0, 0.0, 0.0);  // black lines

  glBegin(GL_LINES);
  glVertex2d(x_for(t_min_scale_bar), 0.03);
  glVertex2d(x_for(t_max_scale_bar), 0.03);
  glEnd();

  // Pick time to display
  auto t_for = [&](const auto& node) -> double {
    return mcc_nodes_are_mrcas ? mcc_tree.at(node).t_mrca() : mcc_tree.at(node).t();
  };
  
  // Draw current tree
  for (const auto& node : index_order_traversal(mcc_tree)) {

    glBegin(GL_LINES);
    if (node != mcc_tree.root) {
      glVertex2d(x_for(t_for(node)), node_ys[node]);
      glVertex2d(x_for(t_for(mcc_tree.at(node).parent)), node_ys[node]);
    }
    if (mcc_tree.at(node).is_inner_node()) {
      glVertex2d(x_for(t_for(node)), node_ys[mcc_tree.at(node).left_child()]);
      glVertex2d(x_for(t_for(node)), node_ys[mcc_tree.at(node).right_child()]);
    }
    glEnd();
  }

  // Draw tips
  auto max_leaf_index = std::ranges::max(leaf_index);
  for (const auto& node : index_order_traversal(mcc_tree)) {
    if (mcc_tree.at(node).is_tip()) {
      auto half_side = 0.0025;
      auto cx = x_for(t_for(node));
      auto cy = node_ys[node];

      auto N = max_leaf_index;
      auto c = static_cast<double>(leaf_index[node]) / N;
      glColor3d(0.0, 1.0 * (1 - c), 1.0 * c);

      glBegin(GL_QUADS);
      glVertex2d(cx - half_side, cy - half_side);
      glVertex2d(cx - half_side, cy + half_side);
      glVertex2d(cx + half_side, cy + half_side);
      glVertex2d(cx + half_side, cy - half_side);
      glEnd();
    }
  }
}

static auto render_all() -> void {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // white background
  glClear(GL_COLOR_BUFFER_BIT);

  render_tree();
  render_mcc();

  glutSwapBuffers();
}

static auto print_stats_line() -> void {
  using enum Real_seq_letter;

  const auto& tree = ui_run->tree();
  
  auto M_ts = [&M_ab = ui_run->num_muts_ab()]() {
    return M_ab[A][G] + M_ab[G][A] + M_ab[C][T] + M_ab[T][C];
  }();

  auto sum_nu = 0.0;
  auto sum_nu2 = 0.0;
  for (const auto& nu : ui_run->nu()) {
    sum_nu += nu;
    sum_nu2 += nu * nu;
  }
  auto mean_nu = sum_nu / tree.num_sites();
  auto mean_nu2 = sum_nu2 / tree.num_sites();
  auto sigma_nu = std::sqrt(mean_nu2 - mean_nu*mean_nu);
  
  auto steps_per_s = 0.0;
  if (timestamps.size() > 1) {
    const auto& earliest_timestamp = timestamps.front();
    const auto& latest_timestamp = timestamps.back();
    steps_per_s = static_cast<double>(latest_timestamp.step - earliest_timestamp.step)
        / (1e-9 * static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(
            latest_timestamp.time - earliest_timestamp.time).count()));
  }
  
  auto state_frequencies_at_root = ui_run->state_frequencies_of_ref_sequence();
  for (const auto& m : tree.at_root().mutations) {
    --state_frequencies_at_root[m.from];
    ++state_frequencies_at_root[m.to];
  }
  for (const auto& [mi_site, mi_from] : tree.at_root().missations.slow_elements(tree.ref_sequence)) {
    --state_frequencies_at_root[mi_from];
  }
  
  std::cerr << absl::StreamFormat("Step %12d, ", ui_run->step())
            << absl::StreamFormat("%.1f M steps/s, ", steps_per_s / 1e6)
            << absl::StreamFormat("%2d trees in MCC, ", tree_snapshots.size())
      //<< absl::StreamFormat("%2d trees in iMCC, ",
      //                           incremental_mcc_tree == nullptr ? 0 : incremental_mcc_tree->num_trees())
            << absl::StreamFormat("%7d steps/refresh, ", steps_per_refresh)
            << absl::StreamFormat("num_muts = %d, ", ui_run->num_muts())
            << absl::StreamFormat("T = %.3g, ", calc_T(tree))
            << absl::StreamFormat("t_MRCA = %s, ", to_iso_date(tree.at_root().t))
            << absl::StreamFormat("log(posterior) = %.3g, ", ui_run->log_posterior())
            << absl::StreamFormat("log(G) = %.3g, ", ui_run->log_G())
            << absl::StreamFormat("log_coal = %.3g, ", ui_run->log_coalescent_prior())
            << absl::StreamFormat("log_other_priors = %.3g, ", ui_run->log_other_priors())
            << absl::StreamFormat("n0 = %.4g, ", ui_run->pop_model().pop_at_t0())
            << absl::StreamFormat("g = %.4g, ", ui_run->pop_model().growth_rate());
  if (not ui_run->mpox_hack_enabled()) {
    std::cerr << absl::StreamFormat("mu = %.2g * 10^-3 subst / site / year, ", ui_run->mu() * 365 * 1000)
              << absl::StreamFormat("M_ts = %d, M_tv = %d, ", M_ts, ui_run->num_muts() - M_ts)
              << absl::StreamFormat("kappa = %.3g, ", ui_run->hky_kappa())
              << absl::StreamFormat("pi = [%.2g, %.2g, %.2g, %.2g], ",
                                    ui_run->hky_pi()[A],
                                    ui_run->hky_pi()[C],
                                    ui_run->hky_pi()[G],
                                    ui_run->hky_pi()[T]);
  } else {
    std::cerr << absl::StreamFormat("mpox_mu = %.2g * 10^-6 subst / site / year, ", ui_run->mpox_mu()*365*1e6)
              << absl::StreamFormat("mpox_mu* = %.2g * 10^-4 subst / site / year, ", ui_run->mpox_mu_star()*365*1e4);

    auto M_apobec = 0;
    for (const auto& node : index_order_traversal(tree)) {
      if (node == tree.root) { continue; }
      for (const auto& m : tree.at(node).mutations) {
        if (ui_run->evo().partition_for_site[m.site] == 1) {
          if ((m.from == Real_seq_letter::C && m.to == Real_seq_letter::T) ||
              (m.from == Real_seq_letter::G && m.to == Real_seq_letter::A)) {
            ++M_apobec;
          }
        }
      }
    }

    std::cerr << absl::StreamFormat("M_apobec = %d of %d, ", M_apobec, ui_run->num_muts());
  }
  
  std::cerr << (ui_run->alpha_move_enabled()
                ? absl::StreamFormat("alpha = %.3g, ", ui_run->alpha())
                : absl::StreamFormat("alpha OFF, "))
            << (ui_run->alpha_move_enabled()
                ? absl::StreamFormat("nu ~ %.3g +/- %.3g, ", mean_nu, sigma_nu)
                : absl::StreamFormat("nu = 1.000, "))
            << absl::StreamFormat("Root Seq counts: [%d, %d, %d, %d], ",
                                  state_frequencies_at_root[A],
                                  state_frequencies_at_root[C],
                                  state_frequencies_at_root[G],
                                  state_frequencies_at_root[T])
            << absl::StreamFormat("num_muts_ab counts: [%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d], ",
                                  ui_run->num_muts_ab()[A][C],
                                  ui_run->num_muts_ab()[A][G],
                                  ui_run->num_muts_ab()[A][T],
                                  ui_run->num_muts_ab()[C][A],
                                  ui_run->num_muts_ab()[C][G],
                                  ui_run->num_muts_ab()[C][T],
                                  ui_run->num_muts_ab()[G][A],
                                  ui_run->num_muts_ab()[G][C],
                                  ui_run->num_muts_ab()[G][T],
                                  ui_run->num_muts_ab()[T][A],
                                  ui_run->num_muts_ab()[T][C],
                                  ui_run->num_muts_ab()[T][G])
            << "\n";
}

//auto reset_imcc() -> void {
//  incremental_mcc_tree = std::make_unique<Incremental_mcc_tree>(max_tree_snapshots, 5, ui_run->tree().size(), ui_run->prng());
//}

static auto idle_func() -> void {

  if (steps_per_refresh <= 0) { return; }

  auto step_granularity = steps_per_refresh;
  if (ui_log_output) {
    step_granularity = std::gcd(step_granularity, cmd->log_every);
  }
  if (ui_trees_output) {
    step_granularity = std::gcd(step_granularity, cmd->tree_every);
  }
  auto target_step = ((ui_run->step()/steps_per_refresh)+1)*steps_per_refresh;
  
  while (ui_run->step() < target_step) {
    auto steps_to_next_grain = ((ui_run->step()/step_granularity)+1)*step_granularity - ui_run->step();
    ui_run->do_mcmc_steps(steps_to_next_grain);

    if (ui_log_output) {
      auto log_now = (ui_run->step() % cmd->log_every) == 0;
      if (log_now) {
        ui_log_output->output_log(*ui_run);
        ui_log_output->flush();
      }
    }

    if (ui_trees_output) {
      auto tree_now = (ui_run->step() % cmd->tree_every) == 0;
      if (tree_now) {
        ui_trees_output->output_tree(*ui_run);
        ui_trees_output->flush();
      }
    }
  }
  
  timestamps.push_back(Timestamp{ui_run->step(), std::chrono::high_resolution_clock::now()});
  if (timestamps.size() > k_max_timestamps) {
    timestamps.pop_front();
  }
  
  tree_snapshots.push_back(std::make_shared<Phylo_tree>(ui_run->tree()));
  while (tree_snapshots.size() > max_tree_snapshots) {
    tree_snapshots.pop_front();
  }
  //incremental_mcc_tree->add_base_tree(ui_run->tree());
  
  static auto i = 0;
  if (++i % mcc_refresh_freq == 0) {
    trees_in_mcc.assign(begin(tree_snapshots), end(tree_snapshots));  // Destroys any dangling snapshots
    auto trees_in_mcc_ptrs = std::vector<Phylo_tree*>{};
    for (const auto& tree_in_mcc : trees_in_mcc) {
      trees_in_mcc_ptrs.push_back(tree_in_mcc.get());  // Removed std::shared_ptr wrapping
    }
    mcc_tree = derive_mcc_tree(std::move(trees_in_mcc_ptrs), ui_run->bitgen());
    //mcc_tree = incremental_mcc_tree->derive_mcc_tree();
  }
  
  glutPostRedisplay();
  print_stats_line();
}

auto randomize_tree_times(Phylo_tree& tree) -> void {
  for (const auto& node : post_order_traversal(tree)) {  // Reset child times before parents
    if (tree.at(node).is_inner_node()) {
      auto earliest_child_t = std::min(tree.at(tree.at(node).left_child()).t, tree.at(tree.at(node).right_child()).t);
      tree.at(node).t = earliest_child_t - absl::Uniform(absl::IntervalOpenClosed, ui_run->bitgen(), 0.0, 10.0);
    }
  }

  auto scratch = Scratch_space{};
  randomize_mutation_times(tree, ui_run->bitgen(), scratch);
}

auto start_log_output() -> void {
  if (ui_log_output == nullptr) {
    auto os_ptr = std::make_unique<std::ofstream>(cmd->log_filename.value_or("output.log"));
    ui_log_output = new Beasty_log_output{os_ptr.get(), true};
    os_ptr.release();
    ui_log_output->output_headers(*ui_run);
  }
}

auto stop_log_output() -> void {
  if (ui_log_output != nullptr) {
    ui_log_output->output_footers(*ui_run);
    ui_log_output->flush();
    delete ui_log_output;
    ui_log_output = nullptr;
  }
}

auto start_trees_output() -> void {
  if (ui_trees_output == nullptr) {
    auto os_ptr = std::make_unique<std::ofstream>(cmd->trees_filename.value_or("output.trees"));
    ui_trees_output = new Beasty_trees_output{os_ptr.get(), true};
    os_ptr.release();
    ui_trees_output->output_headers(*ui_run);
  }
}

auto stop_trees_output() -> void {
  if (ui_trees_output != nullptr) {
    ui_trees_output->output_footers(*ui_run);
    ui_trees_output->flush();
    delete ui_trees_output;
    ui_trees_output = nullptr;
  }
}

auto keyboard_func(unsigned char key, int, int) -> void {
  switch (key) {
    case 's':
    case 'S': {
      t_min_x_axis = std::numeric_limits<double>::max();
      t_max_x_axis = -std::numeric_limits<double>::max();
      rescale_axes();
      std::cerr << "*** AXES RESCALED ***" << std::endl;
      break;
    }
    case 'l':
    case 'L': {
      ladderize_tree = not ladderize_tree;
      if (ladderize_tree) {
        std::cerr << "*** TREE LADDERIZED FOR DISPLAY ***" << std::endl;
      } else {
        std::cerr << "*** TREE ROTATED TO MAXIMIZE OVERLAP WITH INPUT TIP ORDER ***" << std::endl;
      }
      break;
    }
    case 'r': {
      randomize_tree(ui_run->tree(), ui_run->bitgen());
      ui_run->tree_modified();
      std::cerr << "*** TREE RANDOMIZED ***" << std::endl;
      ui_run->set_alpha(1.0);  // Force reset of site relative rates, which may be so low as to cause underflows
      ui_run->repartition();
      tree_snapshots.clear();
      mcc_tree = {};
      //reset_imcc();
      break;
    }
    case 't': {
      randomize_tree_times(ui_run->tree());
      ui_run->tree_modified();
      std::cerr << "*** TREE TIMES RANDOMIZED ***" << std::endl;
      ui_run->invalidate_derived_quantities();
      ui_run->repartition();
      break;
    }
    case 'T':
    case 'B': {
      if (key == 'B') {
        ui_run->set_only_displacing_inner_nodes(not ui_run->only_displacing_inner_nodes());
      } else if (key == 'T') {
        ui_run->set_topology_moves_enabled(not ui_run->topology_moves_enabled());
      }
      if (ui_run->only_displacing_inner_nodes()) {
        std::cerr << "*** ONLY INNER NODE DISPLACEMENT MOVES ENABLED ***" << std::endl;
      } else if (not ui_run->topology_moves_enabled()) {
        std::cerr << "*** ONLY TOPOLOGY-PRESERVING MOVES ENABLED ***" << std::endl;
      } else {
        std::cerr << "*** ALL MCMC MOVES ENABLED ***" << std::endl;
      }
      break;
    }
    case 'C': {
      ui_run->set_repartitioning_enabled(not ui_run->repartitioning_enabled());
      if (ui_run->repartitioning_enabled()) {
        std::cerr << "*** REPARTITIONING ENABLED ***" << std::endl;
      } else {
        std::cerr << "*** REPARTITIONING DISABLED ***" << std::endl;
      }
      break;
    }
    case 'A': {
      if (ui_run->alpha_move_enabled()) {
        ui_run->set_alpha_move_enabled(false);
        ui_run->set_alpha(10.0);
        std::cerr << "*** SITE RATE HETEROGENEITY OFF ***" << std::endl;
      } else {
        ui_run->set_alpha_move_enabled(true);
        std::cerr << "*** SITE RATE HETEROGENEITY ON ***" << std::endl;
      }
      break;
    }
    case 'a': {
      ui_run->set_alpha(10.0);
      std::cerr << "*** ALPHA RESET TO 10.0, ALL SITE RELATIVE RATES (nu) RESET TO 1.0 ***" << std::endl;
      break;
    }
    case 'U': {
      if (ui_run->mu_move_enabled()) {
        ui_run->set_mu_move_enabled(false);
        std::cerr << "*** MUTATION RATE FIXED ***" << std::endl;
      } else {
        ui_run->set_mu_move_enabled(true);
        std::cerr << "*** MUTATION RATE WILL BE INFERRED ***" << std::endl;
      }
      break;
    }
    case 'u': {
      ui_run->set_mu(1.0e-3 / 365.0);
      std::cerr << "*** MU RESET TO 1.0x10^-3 / site / year ***" << std::endl;
      break;
    }
    case '=':
    case '+':  // Shift-=
    case '-':
    case '0': {
      switch (key) {
        case '=': steps_per_refresh = (steps_per_refresh < 1) ? 1 : steps_per_refresh * 10; break;
        case '+': steps_per_refresh = 100'000; break;
        case '-': steps_per_refresh = std::max(0L, steps_per_refresh / 10); break;
        case '0': steps_per_refresh = 0; break;
      }
      std::cerr << "*** " << steps_per_refresh << " steps per refresh ***" << std::endl;
      if (steps_per_refresh > 0) {
         tree_snapshots.clear();
         mcc_tree = {};
         //reset_imcc();
      }
      break;
    }
    case '!': {  // Shift-1
      std::cerr << "*** SINGLE STEP ***" << std::endl;
      steps_per_refresh = 1;
      idle_func();
      steps_per_refresh = 0;
      return;  // idle_func has already refreshed the display and printed a stats line
    }
    case 'q': {
      glutLeaveMainLoop();
      return;
    }
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9': {
      auto new_num_parts = static_cast<int>(key - '0');
      if (new_num_parts == 9) { new_num_parts = 96; }
      ui_run->set_num_parts(new_num_parts);
      std::cerr << "*** USING " << ui_run->num_parts() << " SUBTREES IN PARALLEL ***" << std::endl;
      timestamps.clear();
      mcc_tree = {};
      break;
    }
    case 'm':
    case 'M': {
      if (key == 'm' && max_tree_snapshots > 1) {
        max_tree_snapshots /= 2;
      } else if (key == 'M') {
        max_tree_snapshots *= 2;
      }
      std::cerr << "*** USING UP TO " << max_tree_snapshots << " TREES TO BUILD MCC ***" << std::endl;
      //reset_imcc();
      break;
    }
    case 'p':
    case 'P': {
      if (mcc_nodes_are_mrcas) {
        mcc_nodes_are_mrcas = false;
        std::cerr << "*** MCC NODES ARE DEFINED TRADITIONALLY ***" << std::endl;
      } else {
        mcc_nodes_are_mrcas = true;
        std::cerr << "*** MCC NODES ARE MRCAS ***" << std::endl;
      }
      break;
    }
    case 'n':
    case 'N': {
      if (key == 'n' && mcc_refresh_freq > 1) {
        mcc_refresh_freq /= 2;
      } else if (key == 'N') {
        mcc_refresh_freq *= 2;
      }
      std::cerr << "*** RECALCUlATING MCC EVERY " << mcc_refresh_freq << " TREES ***" << std::endl;
      break;
    }
    default: {
      return;  // Avoid redisplay and stats line
    }
    case 'i':
    case 'I': {
      std::cerr << "*** EXPORTING ANALOGOUS BEAST INPUT FILE ('beast2.xml') ***" << std::endl;
      auto os = std::ofstream("beast2.xml");
      if (os) {
        export_beast_input(*ui_run, os);
      }
      std::cerr << "*** EXPORT COMPLETE ***" << std::endl;
      break;
    }
    case 'o':
    case 'O': {
      if (ui_log_output == nullptr && ui_trees_output == nullptr) {
        start_log_output();
        start_trees_output();
        std::cerr << "*** BEAST-COMPATIBLE OUTPUT ENABLED ***" << std::endl;
      } else {
        stop_log_output();
        stop_trees_output();
        std::cerr << "*** BEAST-COMPATIBLE OUTPUT DISABLED ***" << std::endl;
      }
      break;
    }
    case 'x':
    case 'X': {
      ui_run->set_mpox_hack_enabled(not ui_run->mpox_hack_enabled());
      if (ui_run->mpox_hack_enabled()) {
        std::cerr << "*** MPOX HACK ENABLED ***" << std::endl;
      } else {
        std::cerr << "*** MPOX HACK DISABLED ***" << std::endl;
      }
      break;
    }
  }
  glutPostRedisplay();
  print_stats_line();
}

auto ui_init(int *argc, char **argv) -> void {
  glutInit(argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowPosition(100, 20);
  glutInitWindowSize(1800, 1800);
  glutCreateWindow("Delphy");
}

auto ui_main_loop(Processed_cmd_line& c) -> int {
  cmd = &c;
  ui_run = c.run.get();

  //reset_imcc();

  if (cmd->log_filename.has_value()) {
    start_log_output();
  }
  if (cmd->trees_filename.has_value()) {
    start_trees_output();
  }

  // Set up leaf_index Node_vector
  leaf_index.assign(std::ssize(ui_run->tree()), 0);
  auto next_leaf_index = Node_index{0};
  const auto& tree = ui_run->tree();
  for (const auto& node : index_order_traversal(tree)) {  // index-order => leaf indices match input order in FASTA/MAPLE
    if (tree.at(node).is_tip()) {
      leaf_index[node] = next_leaf_index;
      ++next_leaf_index;
    }
  }

  print_stats_line();

  glutDisplayFunc(render_all);
  glutIdleFunc(idle_func);
  glutKeyboardFunc(keyboard_func);
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

  glutMainLoop();

  stop_log_output();
  stop_trees_output();

  ui_run = nullptr;
  
  return 1;
}

}  // namespace delphy

auto main(int argc, char** argv) -> int {
  using namespace delphy;

  absl::InitializeLog();
  ui_init(&argc, argv);

  auto c = process_args(argc, argv);

  return ui_main_loop(c);
}
