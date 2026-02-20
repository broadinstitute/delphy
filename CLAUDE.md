# Claude Code Instructions

**STOP. You MUST read [AGENTS.md](AGENTS.md) in full before doing anything.**

AGENTS.md contains critical coding conventions, build instructions, and codebase architecture that you must understand before writing any code.

---

## Conceptual Guide to Delphy's Core

Delphy performs Bayesian phylogenetic inference via MCMC. The core loop proposes changes to
a timed phylogenetic tree and accepts or rejects them using the Metropolis-Hastings criterion.
The most complex move is the SPR (Subtree Prune and Regraft) move.

### Missing Data (Missations)

Missing data is tracked per-branch via `Missation_map`, which stores half-open site intervals
`[start, end)` plus `from_states` recording the sequence state entering each missing region.

**Missation factoring invariant:** At each internal node P with children X and S, the
missation intervals of X, S, and P are pairwise disjoint. `P.missations` holds sites missing
below *both* children; `X.missations` and `S.missations` hold sites missing below only that
child. This invariant is maintained by bottom-up factoring (pushing the intersection of
children's missations to the parent). A key consequence: `S.missations ∩ miss_X_total = ∅`,
where `miss_X_total` is the union of all missation intervals from X up to the root.

### Lambda and Site Counting

`lambda_i[node]` is the cumulative mutation rate at a node across all sites present at that
node. It is computed cumulatively from the root:
`lambda_i[child] = lambda_i[parent] + calc_delta_lambda_across_branch(child)`.
Note that `lambda_i[root]` includes the delta for `root.missations`.

`num_sites_missing_at_every_node[node]` (`nsmn`) counts the total number of sites missing
at a node (summing missation intervals from the node to the root).

Both are preserved through tree editing operations.

### SPR Move Architecture

An SPR move prunes the subtree at X and regrafts it elsewhere. The process:

1. **Analyze** (`analyze_graft`): catalogue all mutations on the "hot path" upstream of X
   that would need to change if X were pruned.
2. **Peel** (`peel_graft`): remove hot-path mutations from the tree, placing temporary fake
   mutations on the P-X branch to maintain tree consistency.
3. **Study** (`Spr_study`): explore candidate reattachment regions, weighted by how many
   mutations each would require.
4. **Move** (`move`): physically relocate P to the new graft point via `Tree_editing_session`.
5. **Propose** (`propose_new_graft`): re-analyze at the new position and sample new
   mutational histories for the hot paths.
6. **Accept/Reject**: compute MH ratio from `delta_log_G`, `log_alpha_mut`, and study
   probabilities. Either `apply_graft(new_graft)` or undo via `move` + `apply_graft(old_graft)`.

### Inner vs Rooty Grafts

Dispatched by whether `parent(X) == root` at the time of the call:

- **Inner graft** (P != root): hot paths go from P-X upward through successive ancestors.
  Each level handles a disjoint set of "hot sites" determined by sibling missations.
  The final path may be "open" (unconstrained start state) if it reaches the root.

- **Rooty graft** (P == root): three paths handle three disjoint site categories:
  P-X (sites in `miss_S`, open), P-S (sites in `miss_X`, open), and S-P-X (sites present
  everywhere along the path, closed). Only S-P-X mutations contribute to `min_muts`.

### Hot and Warm Sites

In inner grafts, "warm sites" at a given level are sites whose missations from S are still
sliding upward — their mutational state is not yet fully determined. "Hot sites" are the
subset of warm sites that become determined at that level (because the sibling at that level
is NOT missing those sites). Hot sites are disjoint across branch_infos.

### SPR Study and Candidate Regions

`Spr_study_builder` performs a restricted depth-first traversal of the tree, tracking
`cur_to_X_deltas` (the sequence delta from the current position to X). It skips mutations at
sites missing at X. Each candidate region records its `min_muts = ssize(cur_to_X_deltas)`.

`Spr_study` weights candidate regions by an annealed likelihood factor (controlled by
`annealing_factor < 1`), enabling the MH criterion to add corrections.

### Tree Editing Sessions

`Tree_editing_session` provides a structured way to relocate X's parent P through the tree
via `slide_P_along_branch`, `hop_up`, `flip`, and `hop_down`. The constructor removes
X.mutations and stores them as `deltas_nexus_to_X`; `end()` places them back. During slides,
these deltas and the missation `from_states` are updated to maintain consistency. These
operations preserve `lambda_i`, `nsmn`, and the missation factoring invariant.

### Arena Allocator (Scratch_vector)

`Scratch_vector`, `Scratch_flat_hash_map`, etc. use a thread-local arena (bump) allocator.
Allocation is fast; `free()` is a no-op (memory persists until the arena is reset). This
means old memory remains readable after reallocation, but relying on this is undefined
behavior.

### Mutational History Sampling

`sample_mutational_history` (closed paths): rejection sampling per delta site, plus a
geometric-skip algorithm for null sites. Uses `L = tree->num_sites()` for sampling, then
filters to hot sites. Probability computation uses effective L (hot sites only); correctness
follows from site independence.

`sample_unconstrained_mutational_history` (open paths): Gillespie algorithm with no endpoint
constraint. `adjust_mutational_history` rotates states of non-delta-site mutations to match
the actual sequence state at the path endpoint.
