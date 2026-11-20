# Show the blob tree of the current posterior tree as a third panel in `delphy_ui`

## Summary

`delphy_ui` currently draws two panels side by side: the live MCMC tree (left half, with the
population curve overlaid at the top) and the MCC tree (right half).  This plan adds a third
panel that draws the POPLAR blob tree (`core/poplar.h`) derived from the current posterior tree,
so that the blob decomposition can be watched evolving as the MCMC runs.

Most of the change is in `tools/delphy_ui.cpp`.  The only change in `core/` is adding
`Blob_vector` and `Blob_map` aliases to `core/poplar.h`.  Along the way, `delphy_ui` stops
re-rendering a frame on every pass through its main loop and renders only when something has
changed.

## Background

### How `delphy_ui` is organized today

Everything lives in `tools/delphy_ui.cpp` (about 1,240 lines, SDL2 + SDL2_ttf, software-ish 2D
drawing through a thin shim: `color_3d`, `draw_line`, `draw_rect`, `fill_rect`, `draw_text`).

* **Coordinates.**  `render_all()` sets a unit-square coordinate system with (0,0) at the
  bottom-left (`ox`, `oy`, `dx`, `dy`), and `px()`/`py()` convert to pixels.
* **Panels.**  There is no panel abstraction.  `render_tree()` and `render_pop_curve()` each do
  `dx *= 0.5; defer { dx *= 2.0; };` to confine themselves to the left half, and `render_mcc()`
  does `ox += 0.5*dx; dx *= 0.5; defer { dx *= 2.0; ox -= 0.5*dx; };` for the right half.
* **Shared time axis.**  `x_for(t)` maps time to x using the global `t_min_x_axis` /
  `t_max_x_axis`, updated by `rescale_axes()` from the live tree.  All panels use it, so they
  are time-aligned.
* **Main loop.**  `ui_sdl_loop()` calls `render_all()` on *every* iteration (even when nothing
  has changed and the MCMC is paused, where `idle_func()` just sleeps 10 ms), followed by
  `idle_func()`, which advances the MCMC by `steps_per_refresh` steps (0 means paused), pushes a
  tree snapshot, rebuilds the MCC tree and prints the stats line.
* **Input.**  Keyboard only (`SDL_TEXTINPUT` -> `keyboard_func`).  There is no mouse handling,
  and window events other than `SDL_QUIT` are ignored.  Free letter keys include `v`, `V`,
  `w`, `W`, `e`, `E`, `d`, `D`, `f`, `F`, `j`, `J`, `y`, `Y`.  Upper and lower case are
  sometimes synonyms (`s`/`S`) and sometimes distinct commands (`a` vs `A`, `u` vs `U`).
* **Tree-modifying keys.**  `r` (randomize tree) and `t` (randomize times) change the live tree
  without running any MCMC steps, so cached views of the tree must be invalidated there too.

### What the blob tree API gives us

`derive_blob_tree(const Phylo_tree&) -> Blob_tree`.  A `Blob_tree` is a `Tree<Blob_node>`
with n-ary nodes (`children` is a `std::vector`), and per blob:

* `fingerprint` (64-bit hash of the blob's sequence, identical for identical sequences),
* `nodes`, `branches` (member `Phylo_tree` nodes and branches), `seed_node` / `seed_branch`,
* `active_mutation_extents` (non-empty means the blob has ambiguous sites),
* `sequence_delta` vs. its parent blob.

`Blob_tree` also has `node_2_blob` / `branch_2_blob`.  The generic traversals in `core/tree.h`
(`pre_order_traversal`, `post_order_traversal`, `index_order_traversal`) work on it, since they
are templated on `Tree_like`.

### Things that get in the way

1. **The existing tree helpers are binary-only.**  `canonicalize_tree` and
   `position_tree_nodes` use `left_child()` / `right_child()`, so they cannot be reused for
   the n-ary blob tree.
2. **`render_tree()` calls `canonicalize_tree(ui_run->tree())`, which reorders the live tree's
   children every frame.**  That is harmless for a `Blob_tree` (it is keyed by node index, not
   child order), but it means the blob tree must not be derived lazily from inside a render
   function at an arbitrary time; derive it at well-defined points instead (see below).
3. **Rendering happens on every loop iteration**, which is wasteful by itself (see section 1)
   and would be worse if rendering derived a blob tree each time.  Even with rendering only on
   change, derivation belongs at the points where the tree changes, not in the render path, so
   the derived data must be cached.
4. **The blob tree has no time axis of its own.**  Times have to be recovered from the member
   nodes and branches in the `Phylo_tree`.
5. **Blob layout must not depend on which tips a blob contains.**  Tips can migrate from one
   blob to another between consecutive blob trees, so anything derived from contained tips
   (ordering, positions) would make the picture jump around for no good reason.  Layout may only
   depend on blob tree structure and on blob fingerprints.

## Design

### 1. Render only when something has changed (`ui_dirty`)

Now is as good a time as any to stop rendering a frame on every pass through
`ui_sdl_loop()`.  Add a flag:

```cpp
static bool ui_dirty = true;  // Set by anything that could change what is on screen
```

and change the loop to:

```cpp
while (sdl_running) {
  // ... poll events; set ui_dirty on any that could change the picture (see below) ...

  if (ui_dirty) {
    ui_dirty = false;
    render_all();
  }
  idle_func();
}
```

`ui_dirty` is set by:

* **Any key press**: at the top of `keyboard_func()`, unconditionally (almost every key changes
  something, and the cost of a spurious redraw is nil).
* **`idle_func()` having advanced the MCMC**: at the end of the path that runs steps (this covers
  the single-step key `!` as well, which calls `idle_func()` directly).
* **Window events that invalidate the picture**: handle `SDL_WINDOWEVENT` and set the flag for
  `SDL_WINDOWEVENT_SIZE_CHANGED`, `SDL_WINDOWEVENT_EXPOSED`, `SDL_WINDOWEVENT_RESTORED` and
  `SDL_WINDOWEVENT_SHOWN`.  Today these events are dropped on the floor (the window just gets
  repainted on the next loop iteration, so nothing special is needed), but without continuous
  rendering the window would stay stale or blank after a resize or after being uncovered.
* **Startup**: initialized to `true`, so the first frame is drawn.

Notes:

* When paused (`steps_per_refresh <= 0`), the loop now only polls events and sleeps: no drawing.
* `render_all()` has side effects (`render_tree()` canonicalizes the live tree's child order
  and `rescale_axes()` widens the time axis), so these now happen only on dirty frames.  Both
  are idempotent for an unchanged tree, so skipping them on clean frames changes nothing visible.
* Anything added later that changes the picture outside of these paths must set `ui_dirty`; a
  comment next to the declaration says so.

### 2. Panel layout: a small `Panel_scope` helper

Replace the three ad-hoc `dx`/`ox` manipulations with one RAII helper, which also avoids the
floating-point drift of the `dx *= 0.5` / `dx *= 2.0` round trip:

```cpp
// Temporarily restrict drawing to panel `index` (0-based) out of `num_panels` equal-width
// vertical strips of the window
struct Panel_scope {
  float saved_ox;
  float saved_dx;

  Panel_scope(int index, int num_panels) : saved_ox{ox}, saved_dx{dx} {
    dx /= num_panels;
    ox += index * dx;
  }
  ~Panel_scope() {
    ox = saved_ox;
    dx = saved_dx;
  }
};
```

`render_all()` computes `num_panels = show_blob_tree ? 3 : 2` and wraps each group of render
calls in a `Panel_scope`, so the render functions themselves no longer know anything about
the window layout:

| Panel | Rendered by                              | Notes                                   |
|-------|------------------------------------------|-----------------------------------------|
| 0     | `render_tree()`, `render_pop_curve()`    | population curve still overlaid         |
| 1     | `render_mcc()`                           | unchanged otherwise                     |
| 2     | `render_blob_tree()`                     | new; only when `show_blob_tree`         |

With the blob panel hidden the layout is exactly today's two halves.

### 3. Derive once per tree change, and cache what drawing needs

#### `Blob_vector` and `Blob_map` in `core/poplar.h`

Blobs are indexed by `Blob_index` (an alias for `Node_index`), and per-blob information deserves
the same decoration helpers that nodes and branches have in `core/tree.h` (`Node_vector`,
`Node_map`, `Branch_vector`, ...).  Add, next to `Blob_index` / `k_no_blob`:

```cpp
// Associate information to blobs in the same way as we do to nodes
template<typename T, typename Alloc = std::allocator<T>>
using Blob_vector = Node_vector<T, Alloc>;

template<typename T>
using Blob_map = Node_map<T>;
```

There is currently no container in `poplar.h` / `poplar.cpp` that is indexed *by blob*
(`node_2_blob` and `branch_2_blob` are indexed by node and branch; blob data lives in the
`Blob_node` fields themselves), so the first user is the UI cache below.  They are also what
the TODOs at the end of `poplar.h` (matching blobs across many posterior trees) will want.

#### The cache

All per-blob rendering data is gathered in one struct, held in a `Blob_vector`, instead of a
bunch of parallel vectors:

```cpp
struct Blob_render_info {
  double t_lo = 0.0;              // Time at which the blob starts (see below)
  int depth = 0;                  // Depth in the blob tree, relative to the blob with the reference tip
  int num_descendant_blobs = 0;   // Blobs strictly below this one in the blob tree (for ladderizing)
  double y = 0.0;                 // Vertical position (see 4.)
  std::array<double, 3> color{};  // From the fingerprint (see 5.)
};

struct Blob_view {
  Blob_tree blob_tree;            // Children sorted for display (see 4.)
  Blob_vector<Blob_render_info> info;
  int num_leaf_blobs = 0;
};
static Blob_view blob_view;
static bool show_blob_tree = true;

// Lowest and highest blob depths seen so far, used to scale the depth axis (see 4.)
static int blob_depth_low_watermark = 0;
static int blob_depth_high_watermark = 0;
```

Per-blob data deliberately does *not* include the number of tips behind each blob: nothing
about the picture depends on which tips a blob contains.

`refresh_blob_view()` does the following, using `ui_run->tree()`:

1. `blob_tree = derive_blob_tree(tree)`, followed by `assert_tree_integrity(blob_tree)`
   (self-gated to debug builds).
2. For each blob, compute the **time at which it starts** from the underlying `Phylo_tree`:
   `t_lo` = min of `tree.at(n).t` over member nodes `n`, and `tree.at_parent_of(b).t` over
   member branches `b`.  A child blob starts where its parent's region ends, so it never starts
   before its parent: for a blob seeded at node `n`, `t_lo = t(n)`, which is also the end time
   of the member branch above `n` in the parent blob.  Also derive the blob's color from its
   fingerprint (section 5).
3. Compute `depth` and `num_descendant_blobs` over the blob tree.  Depth is first computed
   from the root blob (pre-order, setting the root blob's depth to 0 explicitly), and then
   shifted so that the blob containing a fixed **reference tip** (the first tip in index order,
   i.e., the first tip in the input) has depth 0, which sets all other depths relative to it.
   The root blob therefore gets a negative depth.  Update `blob_depth_low_watermark` and
   `blob_depth_high_watermark` with the lowest and highest depths found.
   `num_descendant_blobs` is computed in post-order, again setting it explicitly for leaf
   blobs.
4. Lay the blob tree out (section 4): ladderize the children of each blob, and compute the
   vertical positions.  This is done inside `refresh_blob_view()` as well, rather than in a
   separate function, since it only needs the data computed above.

Caching means that all the data needed for drawing lives in `blob_view`, so `render_blob_tree()`
never touches the live `Phylo_tree` (and so never sees a tree that has changed since the blob
tree was derived).

**When to refresh** (only when `show_blob_tree`; otherwise skip the cost entirely):

* at the end of `idle_func()`, after the MCMC steps and the snapshot bookkeeping, before
  `print_stats_line()`;
* in `keyboard_func()` after any key that modifies the live tree without stepping: `r`
  (randomize tree) and `t` (randomize times);
* once in `ui_main_loop()` before entering the SDL loop (so the panel is populated before the
  first step), and whenever the panel is toggled on.

### 4. Layout of an n-ary blob tree

The layout follows the style of `position_tree_nodes()`, adapted to n-ary trees.  Per the
constraint above, **nothing about the layout depends on the tips contained in the blobs**:

* **Child order** (ladderizing): sort each blob's children in the cached `blob_tree` (it is our
  private copy) by ascending `num_descendant_blobs`, so the smaller subtrees come first, as the
  existing ladderization does for trees.  Break ties by ascending blob `fingerprint`, which is
  arbitrary but stable across different blob trees.  Two sibling blobs can have the same
  fingerprint (blobs with equal sequences that recur in the tree are kept distinct), so fall
  back on blob index as a last resort (arbitrary, but deterministic within a tree).  The
  existing `l` toggle (which can rotate trees to match tip order) has no counterpart here and
  does not affect this panel.
* **y**: in post-order, leaf blobs get consecutive slots spaced `0.9 / (num_leaf_blobs - 1)`
  apart starting at 0.05 (same margins as the tree panels); an inner blob sits at the mean
  of the y values of its children (a blob with a single child therefore lines up with it).
  Guard the `num_leaf_blobs == 1` case, where the spacing formula divides by zero.
* **x**: a blob is a point, whose horizontal position depends on a display mode (toggled with
  `F`, see section 6):
  * **time mode** (default): `x = x_for(t_lo)`, so the blob panel stays time-aligned with the
    live and MCC panels;
  * **depth mode**: `x = 0.05 + 0.90 * (depth - low) / (high - low + 1)`, where `low` and
    `high` are the low and high depth watermarks, so blob generations are equally spaced.  This
    is the classic cladogram view, which ignores the (possibly very uneven) time structure and
    shows the topology of the blob tree on its own.

    Two things about the blob tree fluctuate a little from one posterior tree to the next: the
    root blob, and the maximum depth.  A depth axis measured from the root blob and scaled by
    the current maximum depth would therefore make every blob "breathe" distractingly.
    Instead, depth is measured relative to the blob containing a fixed reference tip (see
    section 3), which makes that blob's depth a constant 0.  Only the root blob (at a negative
    depth) and the deepest blobs then move, and the axis is scaled by **low and high
    watermarks**: the lowest and highest blob depths observed so far, which only ever widen.
    Both watermarks are reset to 0 whenever the time scale is rescaled (the `s` key, see
    section 6), at which point the blob view is refreshed immediately so that they restart
    from the current blob tree rather than from 0.

### 5. Drawing

`render_blob_tree()` is called inside `Panel_scope{2, 3}`; `rescale_axes()` is *not* needed
again here (the live tree panel already did it this frame; the cached times use the same
scale).  Everything that is drawn is a function of the blob tree structure and the
fingerprints only:

* **Blobs**: filled circles of one common size at `(x, y)`, drawn in index order.  The radius
  is fixed (4 pixels), shrinking only if leaf blobs are so crowded that neighboring circles
  would overlap (down to a minimum of 1.5 pixels).  SDL has no circle primitive, so add a
  `fill_circle(x, y, radius_in_px)` to the drawing shim next to `fill_rect`, which draws
  horizontal scanlines; the radius is in pixels so the circle stays round even though the
  horizontal and vertical scales differ.
* **Color**: derived from the fingerprint, e.g. three bytes of `fingerprint` mapped into
  `[0.15, 0.85]` per channel so that no color vanishes against the white background.  Equal
  sequences get equal colors (within one tree, and across successive refreshes), so blobs that
  persist while the MCMC runs keep their color, which is the main thing worth seeing in this
  panel.  Other coloring options are left for follow-up work.
* **Connectors**: a thin black straight line from each non-root blob to its parent, drawn
  before the circles.
* **Axis**: in time mode, the same scale-bar line as the other panels at y = 0.03 (omitted in
  depth mode, where there is no time axis).
* **Text**: add `Blobs: N` below the existing `Step:` and `Muts:` lines in `render_all()`.

Blob size is not used to encode anything (no tip counts, no ambiguity marking): both were tried
and did not help, see "Alternatives considered".

### 6. Controls

* `f` / `F` are used for the blob panel, in honor of Frank Liu, who did most of the POPLAR
  work preceding this implementation in Delphy:
  * `f` toggles `show_blob_tree`, with the usual stderr message in the existing style
    (`*** BLOB TREE PANEL SHOWN ***` / `*** BLOB TREE PANEL HIDDEN ***`).  Toggling on calls
    `refresh_blob_view()` immediately.
  * `F` toggles the blob panel's x axis between time mode and depth mode (message:
    `*** BLOB TREE X AXIS IS NOW TIME ***` / `*** BLOB TREE X AXIS IS NOW DEPTH ***`).  This
    only changes how the cache is drawn; no re-derivation is needed.  Upper and lower case are
    distinct commands here, as for `a`/`A` and `u`/`U`.
* The existing `s` / `S` key (rescale the time axis) also resets both blob depth watermarks to
  0 and refreshes the blob view, so that the depth axis restarts from the current blob tree.
* Default: panel shown, x axis in time mode.
* `print_stats_line()` gains a `%d blobs` item.

### 7. Performance

`derive_blob_tree` is `O(M * depth)` in the number of mutations `M`: it maps every mutation to
its extent in the most straightforward way (a documented, deliberate choice in `poplar.cpp`).
In practice this is cheap: in `delphy-web` it takes single-digit milliseconds on the 3.3k-tip
H5N1 tree.  That is why the panel is shown by default and refreshed on every refresh, without
any throttling.  To be confirmed in `delphy_ui`:

* Measure `refresh_blob_view()` on the largest datasets we use and compare with the time per
  refresh implied by the displayed steps/s (add the derivation time in ms to the stats line
  while testing).
* Only if that turns out to matter, mitigations in increasing order of invasiveness: (a) refresh
  only every k-th `idle_func()` call, as is done for the MCC tree via `mcc_refresh_freq`; (b) speed
  up `map_tree_mutation_extents` in core (a separate piece of work, not part of this plan).  The
  panel can always be hidden with `f`.

## Files changed

| File                  | Action | What |
|-----------------------|--------|------|
| `core/poplar.h`       | Modify | Add `Blob_vector` and `Blob_map` aliases next to `Blob_index` |
| `tools/delphy_ui.cpp` | Modify | Add `#include "poplar.h"`; add `ui_dirty` and render only when it is set (key presses, MCMC progress, window events); add `Panel_scope`, remove the ad-hoc half-screen code from `render_tree`, `render_pop_curve` and `render_mcc`, and apply the scopes in `render_all()`; add `fill_circle` to the drawing shim; add `Blob_render_info`, `Blob_view`, `show_blob_tree`, the blob depth low and high watermarks, the x-axis mode, `refresh_blob_view()` (derivation, layout and colors), `blob_x()` and `render_blob_tree()`; update `render_all()`, `idle_func()`, `ui_main_loop()`, `ui_sdl_loop()`, `keyboard_func()` (new `f` and `F` keys, refresh after `r` and `t`, watermark reset and refresh on `s`) and `print_stats_line()` |

No changes to `tests/`, the WASM interface, or the CMake files: `delphy_ui` already links
`delphy_core`, which provides `poplar.h`.

## Alternatives considered

* **Aligning each blob's y with the tips it contains in the live tree panel.**  Besides the
  constraint that layout must not depend on contained tips, nested blobs would overlap
  vertically (a parent blob covers the tips of all its descendants), so a separate dendrogram
  layout is clearer.  Highlighting the tips behind a blob would be a natural follow-up, but
  needs mouse support, which `delphy_ui` does not have.
* **Putting the layout code in `core/` with unit tests.**  The layout is a pure function of the
  blob tree and the `Phylo_tree`, so it could live there, but `delphy_ui` is a development tool
  whose drawing code (including `position_tree_nodes`) all lives in one file.  Keeping the
  new code there matches the existing style; it can be lifted into `core/` later if the web
  UI wants to share it.
* **Ladderizing by number of tips, or breaking ties by tip indices.**  Rejected because tips
  migrate between blobs from one blob tree to the next, which would make the layout jitter.
* **Drawing blobs as horizontal bars spanning their time range, with thickness proportional to
  the number of tips (and an outline on blobs with ambiguous sites).**  This was the first
  implementation, and it did not work: the thickness varies wildly and obscures the structure,
  and the outline did not tell anything useful.  Replaced by same-size circles at the start
  time of each blob, with no size encoding at all (not even a size scaling with the number
  of tips, which would also change from one blob tree to the next as tips migrate).
* **Elbow connectors**, in the style of the tree panels.  Straight lines between parent and
  child blobs look better now that blobs are points.
* **Scaling the depth axis by the maximum depth of the current blob tree.**  Rejected because
  the maximum depth fluctuates from one posterior tree to the next, which makes the axis breathe
  distractingly; replaced by the watermarks described in section 4.
* **Measuring depth from the root blob.**  Rejected because the root blob itself fluctuates
  between posterior trees, which shifts the depth of every other blob at the same time.  Depth
  is measured from the blob containing a fixed reference tip instead.

## Verification

`delphy_ui` has no automated tests, so verification is by running it:

1. Build the `delphy_ui` target in both debug and release configurations.  Confirm no new
   warnings (`-Wall -Wextra -Werror`).
2. Run `delphy_ui` with the same command line as `delphy` (see `delphy --help`; for example
   `--v0-init mp-plus-timing` against one of the paper datasets under
   `~/now/delphy-2026-paper-data/<dataset>/delphy_inputs/*.maple`) and check, with the MCMC
   paused (`0`) and then running (`=`):
   * three equal-width panels appear (live tree, MCC tree, blob tree), and `f` switches
     between that and today's two halves with no distortion or leftover offsets;
   * the number of circles matches `Blobs: N` and the stats line;
   * in time mode, circles are time-aligned with the live tree panel (e.g. the root blob sits
     at the root's time); `F` switches to depth mode, where blob generations are equally
     spaced, and back; all circles have the same size and are joined to their parents by
     straight lines;
   * with the MCMC paused, colors and layout are stable from frame to frame; with it running,
     blobs that persist across refreshes keep their color, and the order of blobs does not
     jitter as tips migrate between blobs;
   * in depth mode while the MCMC runs, the depth axis does not breathe: the blob containing
     the first tip stays put, the root blob (at a negative depth) may move left as long as it
     reaches new lows, and the axis only ever stretches (as a taller or deeper blob tree is
     observed); `s` resets it to the current extent;
   * `r`, `t`, and single step (`!`) all update the blob panel.
3. `ui_dirty` behavior, with the MCMC paused: CPU use of `delphy_ui` is near idle (no
   continuous rendering); every key press, resizing the window, and uncovering or restoring it
   redraw correctly; running the MCMC (`=`) redraws after every refresh; and there is no
   visual regression in the existing two panels.
4. Debug builds: confirm `assert_tree_integrity(blob_tree)` and the blob tree `CHECK`s do not
   fire over a few minutes of running on several datasets (this exercises `derive_blob_tree` on
   thousands of real posterior trees, which our unit tests do not).
5. Measure the cost of `refresh_blob_view()` as described in section 7 and record the numbers
   in this plan, in particular to confirm the single-digit-millisecond figure seen in
   `delphy-web` on the large H5N1 tree.
6. Run the full unit test suite and build the WASM target to confirm nothing else was affected
   (the only change in `core/` is the pair of aliases in `poplar.h`).

## Decisions made during review

* The blob panel is shown by default.
* `f` toggles the panel and `F` toggles the panel's x axis between time and depth.
* Fingerprint coloring is the initial scheme; other options (by tips, by ambiguities, by depth)
  are for follow-up work.

## Changes made after trying the first implementation

After running the first implementation, these changes were made (and are reflected above):

* Blobs are drawn as same-size filled circles at their start time (or at their depth), joined
  to their parents by straight lines, instead of bars with thickness proportional to the number
  of tips.  The number of tips behind a blob is no longer computed at all.
* Blobs with ambiguous sites are no longer highlighted.
* The depth axis is scaled by low and high watermarks of the depths observed, reset when the
  time scale is rescaled, instead of by the maximum depth of the current blob tree.
* Depth is measured relative to the blob containing a fixed reference tip (the first tip in the
  input), because the root blob oscillates and would make the depths of all blobs oscillate
  with it.  The root blob consequently has a negative, fluctuating depth, which is why there is
  a low watermark as well as a high one.
