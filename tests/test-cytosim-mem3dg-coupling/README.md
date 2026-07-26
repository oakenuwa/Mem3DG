# Coupled Cytosim + Mem3DG test dataset

A short coupled run, kept as a fixture for testing the Blender import path and
as a regression check that the coupling still transfers force in both
directions.

## Contents

Two scenarios, differing in membrane topology.

| Path | What |
|---|---|
| `config.cym`, `run.sh`, `frames/` | **OPEN** membrane: `../frame9.ply`, a cap with 78 boundary edges (Euler 1) |
| `config-icosphere.cym`, `run-icosphere.sh`, `frames-icosphere/` | **CLOSED** membrane: `../icosphere.obj`, R=0.8, 642 v / 1280 f (Euler 2, no boundary) |

Each `frames*/` holds `membrane_%04d.obj` + `fibers_%04d.txt`.

### Which to use

**Prefer the icosphere.** A closed surface is the only one for which "inside"
is defined, so it is the prerequisite for any future `SpaceMem3DG`
(`confine = inside, ...`) and the only one where enclosed volume is a
meaningful diagnostic. It is also smaller (1280 vs 3324 faces, so faster) and
couples far more strongly — filaments genuinely straddle the surface:

```
                   contacts   |F_fib|
open cap           58 - 151   0.01 - 0.04
closed icosphere   23 -  46   2.0  - 10.6
```

The open cap is retained because it was the original test asset and exercises
the boundary-edge path.

### Two steric modes (driver arg 9)

| mode | function | behaviour |
|---|---|---|
| `0` | `mem3dg_system_filament_repulsion` | **symmetric**: pushes a segment off the surface on whichever side it is, within `cutoff`. Works on open or closed meshes. **Does not confine.** |
| `1` | `mem3dg_system_filament_confinement` | **one-sided wall** at `R − shell`: keeps filaments *inside* and pushes the membrane *outward*. Requires a **closed** mesh (returns −1 otherwise). |

Mode 0 cannot confine, and this is measurable rather than theoretical. In an
earlier icosphere run under mode 0 the fraction of filament points outside the
membrane rose monotonically **80.2% → 89.6%**: a point that drifts outside is
pushed *further* out, and once it is more than `cutoff` away it feels nothing
and never returns.

Mode 1 uses the outward face normal of the *nearest* face to get a signed
distance `s`, and for `s > −shell` applies `∓stiffness·(s + shell)·n̂_out`
— inward on the filament, outward on the membrane. There is no cutoff, so an
escaped filament is always pulled back.

### The closed-membrane run (mode 1, growing filaments)

**2 filaments**, initial length 0.5, `growing_speed = 2.0`, seeded inside
`r < 0.3`; membrane initially at `r = 0.8`, `shell = 0.05`.
**10 000 steps → 1001 exported frames** (~25 s, 63 MB).

```
frame  totLen   max_r   membrane_r   escaped    volume    dVol%   closed
    0    1.00   0.3353    0.8000        0      2.12620  +0.000%    yes
  200    1.80   0.7005    0.8000        0      2.12620  +0.000%    yes
  400    2.57   0.7697    0.8093        0      2.12659  +0.018%    yes
  600    3.27   0.7990    0.8183        0      2.12832  +0.100%    yes
  800    4.00   0.7986    0.8180        0      2.12987  +0.172%    yes
 1000    4.76   0.7751    0.8177        0      2.13151  +0.249%    yes
```

- filaments grow **4.8×** in total length (1.00 → 4.76)
- **zero escapes** at every sampled frame
- membrane radius grows **0.8000 → 0.8183** and volume rises **monotonically to
  +0.249%** — the filaments visibly inflate it
- mesh verified **closed** throughout

Note the escape test is **membrane-relative**: each filament point is compared
against the membrane radius *in its own direction* at that frame. A fixed
`r > 0.8` threshold would be wrong here precisely because the membrane inflates
past 0.8 — filaments legitimately follow it outward while still being inside.

Both sequences are written by `coupled_driver` **in the same loop**, so
`membrane_0003.obj` and `fibers_0003.txt` are the same simulated instant by
construction. Coordinates are raw micrometres in a shared frame — see
`../../blender/COORDINATES.md`.

## Viewing in Blender

`../../blender/import_trajectory.py` already points at this directory. Open it
in Blender's Scripting workspace and run it, then scrub the timeline. Or:

```bash
blender --python /Users/hakenuwa/Desktop/Mem3DG/blender/import_trajectory.py
```

## Regenerating

```bash
./run.sh
```

## What this run demonstrates

The filaments and membrane genuinely interact. Contacts are counted as
(filament segment, membrane triangle) pairs by Mem3DG's edge-face kernel,
reached through `mem3dg_system_filament_repulsion()`:

```
step    0  contacts    82  |F_fib| 0.01435  |F_mem| 0.02514
step   90  contacts   147  |F_fib| 0.03801  |F_mem| 0.05047
step  190  contacts   100  |F_fib| 0.01633  |F_mem| 0.02962
```

An earlier version computed contacts in the driver, point-to-vertex, and saw
only 2-7. That test is blind to a filament crossing the middle of a triangle
away from its corners -- which is exactly how filaments tunnel through a
membrane. The kernel now measures true segment-to-triangle distance, including
the piercing case.

This matters because it uses a **physical cutoff (0.05)**. Earlier smoke tests
used `cym/actin.cym`, whose filaments occupy roughly ±0.13 while the membrane
spans ±0.78 — they never touched, and the only way to get contacts was to force
the cutoff to 1.0, which is not a physical interaction. The config here places
filaments in the membrane's own region instead.

Measured extents at frame 0 (they overlap, which is the registration check):

```
membrane  x[-0.78 0.68]  y[-0.74 0.68]  z[ 0.00 0.71]
fibers    x[-0.67 0.87]  y[-0.84 0.80]  z[-0.44 0.87]
```

## Caveats

- **Coupling stiffness has to be large to matter.** This run uses `20.0`. At the
  earlier default of `0.02` the membrane moved *identically to zero coupling*
  (max |dz| = 2.4265e-04 either way) — the filament force was ~1e-2 against the
  membrane's own bending force of ~1e4, six orders of magnitude smaller. If you
  change parameters, re-check against a `stiffness = 0` run: if the deformation
  is unchanged, the coupling is not doing anything physically.
- Deformation over this run rises smoothly to max |dz| ≈ 1.0e-3 and plateaus —
  stable, no blow-up. Explicit splitting means raising stiffness further will
  eventually need a smaller `dt`.
- Cost is O(nSegments × nFaces) per contact evaluation, with a centroid
  quick-reject inside Mem3DG. A BVH is the next optimisation.
- Frames are `.obj`, not `.ply`: `mem3dg_system_save_ply()` forwards to
  geometry-central's writer, which only implements OBJ. Requesting `.ply`
  silently produces a 0-byte file.
