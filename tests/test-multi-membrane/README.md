# Multi-membrane steric interaction

Reusable demos of `MultiMembraneSystem` — several closed membranes
("organelles") that are fully independent except for short-range face-face
excluded-volume repulsion.

All cases run the same binary and differ in how the load is applied:

| case | setup | what it stresses |
|---|---|---|
| `approach` | two organelles, constant body force pushing their centres of mass together | can the repulsion **hold a load**? |
| `inflate` | two organelles starting close; one's preferred volume is raised so it swells | can it **transmit a quasi-static push** through a growing contact patch? |
| `cluster` | N ≥ 3 on a ring, each pushed toward the origin until they jam | the **N-body path** — `MultiMembraneSystem::step()` and its C(N,2) pair loop |
| `bvhcheck` | evaluates the BVH and exhaustive narrow phases from the same state | that the BVH acceleration changes **no physics** |

## Build

Needs the C++20 API directly — `MultiMembraneSystem` is a C++ class with no
plain-C entry point, unlike the Cytosim coupling.

```bash
cmake -S tests/test-multi-membrane -B build-multi -DCMAKE_BUILD_TYPE=Release -DMEM3DG_ROOT="$MEM3DG" -DMEM3DG_BUILD="$BUILD/mem3dg" -DEIGEN_SRC="$BUILD/mem3dg/_deps/eigen-src" && cmake --build build-multi -j8
```

`MEM3DG_WITH_NETCDF` is read out of the library's own `CMakeCache.txt` and
matched automatically. Do not hardcode it — see the note at the bottom.

## Case 1 — pushing centres of mass together

```bash
multi_membrane_driver approach tests/icosphere.obj frames-approach 6000 1e-4 0.10 500 500 1.0
```

Arguments: `<scenario> <mesh> <outDir> [nSteps] [dt] [d] [mu] [saveEvery] [push] [grow] [Kv] [Ksg]`

Two equal organelles (R = 0.8, 1280 faces), contact distance `d = 0.10`,
stiffness `mu = 500`, body force 1.0 per vertex on each, inward.

| step | surface gap | centre distance | volume | \|F_rep\| | Newton residual |
|---|---|---|---|---|---|
| 0 | 0.4800 | 2.0800 | 2.126 | 0 | 0 |
| 1600 | 0.1600 | 1.7600 | 2.126 | 0 | 0 |
| 2000 | 0.0998 | 1.6800 | 2.126 | 6.98 | 0 |
| 3600 | 0.0998 | 1.3778 | 2.060 | 71.2 | 6.8e-13 |
| 6000 | 0.0998 | 0.9854 | 1.779 | 172.1 | 5.7e-13 |

Reading it: the gap closes at a constant 2e-4/step with **zero** repulsion
force — the cutoff is respected, nothing acts at range. At step 2000 the gap
reaches `d` and then **locks at 0.0998 for the remaining 4000 steps** while the
push keeps loading it. The centre distance keeps falling (2.08 → 0.985) because
the organelles flatten against each other rather than merging, and their volume
drops 2.126 → 1.779 — this scenario applies no volume constraint, so squashing
them is allowed. `|F_rep|` climbs monotonically to 172 as the load grows: a
harmonic wall pushing back harder, at 0.2% compression of `d`.

## Case 2 — inflating one organelle against the other

```bash
multi_membrane_driver inflate tests/icosphere.obj frames-inflate 12000 1e-4 0.10 500 750 0 2.0 200 0.01
```

Organelle A is full size, B is scaled to 0.8. They start `1.5*d` apart —
just clear of contact, because the thing under test is the swelling, not the
travel. A's preferred volume is set to 2.0x its initial value.

| step | surface gap | centre distance | V_A | V_B | \|F_rep\| |
|---|---|---|---|---|---|
| 0 | 0.1500 | 1.5900 | 2.126 | 1.0886 | 0 |
| 750 | 0.0997 | 1.5902 | 2.879 | 1.0881 | 3.72 |
| 1500 | 0.0998 | 1.5914 | 3.506 | 1.0849 | 6.06 |
| 3000 | 0.1000 | 1.5933 | 4.099 | 1.0803 | 1.70 |
| 12000 | 0.1000 | 1.5938 | 4.252 | 1.0791 | 0.015 |

A reaches 4.2520 against a 2.0x target of 4.2524 — 0.01%. The gap closes to
exactly `d` and holds. Three things confirm the push is really being
transmitted rather than absorbed: the centre distance grows monotonically
(B is displaced), **B's own volume drops** 1.0886 → 1.0791 (B is being
compressed, not just translated), and `|F_rep|` peaks at 6.06 during the fast
growth phase then decays to 0.015 as A stops growing and the system relaxes —
the load history, not a constant.

## Case 3 — three or more organelles (`cluster`)

The other two scenarios hand-roll their step loop and call
`computeInterMembraneRepulsionForce` on the single pair directly, so
**`MultiMembraneSystem::step()` never ran** and neither did
`applyPairwiseRepulsion()`'s C(N,2) loop. With two membranes that loop has one
iteration, so it cannot expose the failure mode that matters at N ≥ 3: a
membrane must *accumulate* force from N−1 partners into one buffer, and a pair
that is skipped, double-counted, or written to the wrong `System` is invisible
in a 2-body test.

`cluster` therefore calls `multi.step()` deliberately, and drives the push
through `parameters.external.form` — the only hook that fires *inside* `step()`,
via `addNonconservativeForcing()`. Injecting the force afterwards would mean
bypassing the code under test again.

N organelles sit equally spaced on a ring and are each pushed toward the origin
until they jam.

```bash
multi_membrane_driver cluster tests/icosphere.obj frames-cluster5 3000 1e-4 0.10 500 1000 1.0 0 0 0 5
```

The last argument is N. Results at step 3000:

| N | pairs | contacts | ring radius | min gap | \|F\|max | momentum residual | runtime |
|---|---|---|---|---|---|---|---|
| 3 | 3 | 3/3 | 0.7385 | 9.942e-02 | 134.1 | 1.6e-15 | 45.6 s |
| 4 | 6 | 4/6 | 0.9529 | 9.942e-02 | 76.4 | 1.6e-15 | 56.3 s |
| 5 | 10 | 5/10 | 1.1970 | 9.967e-02 | 42.3 | 1.7e-15 | 64.8 s |
| 6 | 15 | 6/15 | 1.4549 | 9.965e-02 | 25.4 | 1.9e-15 | 71.9 s |

**The new check is global momentum.** After `step()`, each member's
`selfAvoidanceForceVec` holds only its inter-membrane contribution
(`selfAvoidance.mu` is 0 and `computeConservativeForcing()` zeroed the buffer),
so summing over all N totals the entire N-body repulsion, which must vanish. It
stays at ~1e-15 for every N. A pairwise check cannot make this statement.

`contacts` is N, not C(N,2), and that is correct: on a ring only adjacent
neighbours touch at this level of collapse — diagonals are farther apart
(2ρ vs 2ρ·sin(π/N)). For N = 4 a diagonal does come into contact if the run is
pushed further, and the count rises to 5/6.

**Symmetry is reported, not asserted.** The organelles are *translated* copies of
one mesh, not rotated ones, so the icosphere's triangles do not line up with the
ring's N-fold axis and each contact patch samples different faces. A few percent
spread is discretization. The aggregate is also restricted to steps where the
contact set is adjacent-only, since once diagonals touch the configuration
genuinely is not symmetric and a large spread there would be right physics.

## Performance

Two optimizations, in the solver rather than the test.

**1. Organelle-level broad phase** (`applyPairwiseRepulsion`). One bounding
sphere per membrane, computed once per step, rejects whole pairs before any face
work: two closed surfaces can only come within `d` if
`‖c_a − c_b‖ − r_a − r_b < d`, which is a true lower bound on their
surface-surface distance. O(N²) cheap tests replace O(N²·nF²).

**2. BVH over faces with refit** (`face_bvh.h`). One tree per membrane, built
once and refitted every step — refit is O(nF) with no allocation, whereas a
rebuild is O(nF log nF) and is only triggered when the face count changes, which
is how remeshing is detected. A BVH rather than a uniform grid because Mem3DG's
mesh mutator deliberately produces widely varying triangle sizes, which is
exactly where a grid's fixed cell size degrades.

Measured, 2000 steps:

| | N = 3 | N = 5 |
|---|---|---|
| broad phase only | 62.1 s | 103.9 s |
| broad phase + BVH | 25.4 s | 37.0 s |

Before either change, N = 5 did not finish 6000 steps in 10 minutes (it reached
about step 3000), so the combined speedup at N = 5 is roughly an order of
magnitude. Scaling also improved: from N = 3 to N = 6 the pair count grows 5×
(3 → 15) but runtime only 1.6×.

### Verifying the BVH did not change the physics

Comparing two trajectories cannot answer this. Contact dynamics amplify
round-off, so merely reordering the summation into a vertex shows up as a
visibly different `|F|max` after a few thousand steps — 85.88 vs 88.92 in the
first comparison here. That says the runs diverged, not that the physics differs.

The `bvhcheck` scenario instead evaluates both paths from the **same**
configuration, one step at a time, and reports the per-vertex difference:

```bash
multi_membrane_driver bvhcheck tests/icosphere.obj /tmp/bc 1200 1e-4 0.10 500 0 1.0 0 0 0 3
```

```
pair evaluations compared : 7206
evaluations that differed : 5031
worst per-vertex |dF|     : 2.664535e-15
peak per-vertex |F|       : 3.119635e+00
worst relative difference : 8.541e-16
verdict: IDENTICAL to round-off
```

So the BVH finds the identical contact set; 5031 evaluations differ in the last
bit because the tree returns candidates in traversal order rather than index
order. That is what amplified into the `|F|max` gap above.

The reason the sets match is that the BVH prune is conservative: if two
triangles are within `t`, their boxes padded by `t` must overlap on every axis.
Pruned pairs are farther apart than `d`, where the harmonic force is identically
zero.

## Physics verification

**Newton's third law.** The repulsion is isolated by snapshotting
`mechanicalForceVec`, calling `computeInterMembraneRepulsionForce()`, and
differencing. `ΣF_A + ΣF_B` must vanish, which holds only if the DDG barycentric
weights sum to 1 on *both* faces.

| case | worst absolute | worst relative |
|---|---|---|
| `approach` | 2.6e-12 | 2.7e-12 |
| `inflate` | 6.5e-14 | 1.7e-13 |

Round-off. This is the strongest single check on the force decomposition.

**Non-interpenetration, against a control.** Each case was also run with
`mu = 0`. That comparison is the point — "they did not interpenetrate" means
nothing without a run that did:

| | `approach` | `approach` control | `inflate` | `inflate` control |
|---|---|---|---|---|
| min gap ever | 9.98e-02 | 0.0 | 9.97e-02 | 0.0 |
| verdict | NONE | **DETECTED** | NONE | **DETECTED** |
| neighbour volume | 1.779 (squashed) | 2.126 (untouched) | 1.0791 | 1.0886 (untouched) |
| centre distance | 0.985 | 0.880 | 1.5938 | 1.5900 (never moved) |

With `mu = 0` the meshes pass straight through each other as ghosts: the
neighbour's volume never changes and, in `inflate`, its centre never moves at
all — no force is transmitted.

**Independent cross-check of the metric.** The driver's brute-force
minimum-surface-distance was re-derived from the exported OBJs with a separate
vectorized point-to-triangle routine:

| frame | vertex-vertex | A_vert→B_tri | B_vert→A_tri | driver `sep` |
|---|---|---|---|---|
| 0000 | 0.15000 | 0.15000 | 0.15000 | 0.15000 |
| 0008 | 0.10023 | 0.10023 | 0.10000 | 0.10000 |
| 0016 | 0.10056 | 0.10022 | 0.10000 | 0.10000 |

Note that comparing x-extents of the two meshes is *not* a valid proxy here and
will suggest a smaller gap: once the caps flatten, their extreme-x vertices sit
on the rim of the contact disc and are laterally offset from each other, so they
are not the closest pair.

## Viewing

```
frames-<case>/membrane0_%04d.obj
frames-<case>/membrane1_%04d.obj
frames-<case>/metrics.csv
```

Point `blender/import_trajectory.py` at the directory. It picks up every
`membrane<K>_` sequence, imports them into one frame, and gives each organelle
its own colour so neighbours are tellable apart.

## Bugs this exercise found

1. **`MEM3DG_WITH_NETCDF` struct-layout mismatch** — the worst one. That define
   adds members to `System` and `Integrator`; compiling the test without it while
   the library has it gave the two different field offsets for the same class.
   Not a link error — *corruption*. The symptom was a crash in
   `NcFile::close()` at exit, on a pointer whose bit pattern was a double,
   **after** the run had printed plausible-looking physics. The first
   `approach` run reported the organelles closing 160x faster than they really
   do. The CMakeLists now reads the flag out of the library's `CMakeCache.txt`
   so it cannot be got wrong.
2. **`triTriDist` missed the piercing case** — it covered vertex-triangle and
   edge-edge, but two triangles can cross with an edge through the other's
   interior, every vertex far away and no edge-edge pair touching. Those pairs
   reported a *positive* distance, so the `inflate` control was scored
   "interpenetration: NONE" while A had grown straight through B. Fixed with a
   Möller–Trumbore segment-triangle test, matching the sub-case
   `minimumEdgeFaceContact()` already handles in the solver.
3. **Separation printed as `%.6f`** — a gap of 1e-9 displayed as `0.000000`,
   and the verdict tested `minSep > 0.0` exactly. A failing run could read as a
   passing one. Now scientific notation with a 1e-12 tolerance.
4. **`inflate` started too far apart** (0.64) to ever reach contact, and
   `push`, `grow`, `Kv`, `Ksg` were hardcoded at values that produced no visible
   motion. All four are now CLI arguments; `inflate` derives its initial
   clearance from `d`.
5. **The face-pair reject was not a conservative bound.** It tested
   `‖cA − cB‖ > cutoffFactor · d0`, which ignores triangle size entirely. On the
   1280-face icosphere the circumradius is ≈0.070, so a correct reject needs
   `‖cA − cB‖ ≤ d0 + 0.139`; with `d0 = 0.1` the old `cutoffFactor = 3` gave
   0.3 and passed with 1.25× margin — by luck. At the `d0 = 0.0035` the Cytosim
   coupling uses it gave **0.0105**, far inside `rA + rB`, silently discarding
   real contacts. Now `‖cA − cB‖ − rA − rB > d0` with `cutoffFactor` demoted to
   optional extra slack that cannot cause a miss.
6. **BVH right child is not at `left + 1`.** The left call emits its entire
   subtree first, so `left + 1` is the left child's own left child. Assuming
   adjacency searches the wrong subtree and loses contacts silently; both child
   indices are now stored.
7. **BVH query dropped nodes on stack overflow**, which would also have meant
   missed contacts. The fallback now emits the node's whole face range — a
   superset the narrow phase filters, i.e. slower rather than wrong.
8. **The Blender importer never matched multi-membrane files** — it globbed
   `membrane_*.obj`, which does not match `membrane0_0000.obj` (the next
   character is a digit, not an underscore). Multi-membrane runs wrote frames
   that the viewer silently ignored.

## Caveat on the non-interpenetration check

The geometry check is sampled at `saveEvery`, not every step, because the
brute-force verifier is O(nF²) per pair. Interpenetration that appeared and
healed between two probes would be missed. The momentum check in `cluster` does
run every step. If you need the geometric guarantee at full resolution, lower
`saveEvery` and accept the cost.

## Tuning note

`Kv` and `Ksg` compete: the osmotic term drives volume up, surface tension
resists the area increase that requires. Laplace balance (`p ≈ 2σ/R`) means
nothing inflates unless `Kv` is a few times `Ksg`. The shipped defaults
(`Kv = Ksg = 1e-2`) sit right at that boundary and produce 0.02% growth in 8000
steps — which is why the case above passes `Kv = 200, Ksg = 0.01`.
