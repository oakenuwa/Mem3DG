# Placing singles / couples on a Mem3DG membrane

This test checks that Cytosim's ordinary placement and confinement keywords work
against a Mem3DG membrane, exactly as they do against a built-in `sphere` or
`capsule`.

## How it works

The driver owns the membrane; Cytosim borrows it.

```
coupled_driver                          libcytosim
  mem3dg_system_create_from_file  ->  SpaceMem3DG::setHandle(mem)   (before Parser)
                                            |
                                      shape = mem3dg
                                            |
                                      inside() / project() / boundaries()
                                            |
                                      mem3dg_system_inside/_project/_bounds
```

`SpaceMem3DG` holds only a `static void*` handle -- it never owns or frees the
membrane, so there is exactly one membrane and no lifetime question. The handle
must be set *before* `Parser::readConfig`, because `position = surface` is
resolved during parsing.

Cytosim's base `Space` derives every placement/confinement operation from three
primitives, so forwarding those three is enough to get all the keywords:

| Cytosim call     | built from              | config keyword       |
|------------------|-------------------------|----------------------|
| `place()`        | `boundaries` + `inside` | `position = inside`  |
| `placeOnEdge(R)` | `boundaries` + `project`| `position = surface` |
| `bounce()`       | `project` + `inside`    | `confine = inside`   |

Singles and couples accept only `off | on | inside` -- there is no `outside`
key for them (only fibers have one). See `single_prop.cc`.

## Speed

`project` and `inside` are the hot path (once per object per step), so both go
through a KD-tree rather than scanning all faces:

- `NearestNeighborFinder` (nanoflann, via `geometrycentral/utilities/knn.h`)
  finds the k=8 nearest *vertices*, and only their incident faces are tested.
  Faces are de-duplicated with a generation stamp instead of a `set`.
- The 1-ring is expanded only when the nearest vertex is not clearly closer than
  the best face found so far -- the case where the true closest point could lie
  on a face not incident to any of the k vertices.
- `inside` uses the **angle-weighted pseudonormal** (Bærentzen & Aanæs) of
  whatever feature the closest point lands on, so the sign is exact at edges and
  vertices, where a plain face normal is ambiguous. This replaces ray casting;
  geometry-central's `MeshRayTracer` is not compiled in this build.

Measured on the remeshed 3806-face mesh: **0/3000 `project` results differ from
brute force** (worst 2.2e-16), **0/2886 `inside` results wrong**, and it is
**10x faster** (5.27 vs 53.76 us/query).

## Running it

```bash
coupled_driver config.cym ../icosphere.obj 200 0.05 100.0 1 100 frames 1 8.22e-6 100.0 100.0 0.98 0.0 0
```

`verify_placement.cc` is the check: it parses the config, measures every object's
distance to the membrane surface, steps the simulation, and measures again.

## Result

Distance to the membrane surface, 642-vertex icosphere, before and after 3000 steps:

| object          | placement          | confine  | dist after 3000 steps | inside |
|-----------------|--------------------|----------|-----------------------|--------|
| surfaceAnchor   | `surface 0.01`     | (fixed)  | max 1.6e-16           | 200/200|
| surfaceDiffuser | `surface 0.01`     | `on`     | max 1.1e-16           | 27/100 |
| innerDiffuser   | `inside`           | `inside` | mean 0.235            | 100/100|
| surfaceCross    | `surface 0.01`     | `on`     | max 1.6e-16           | 34/150 |
| innerCross      | `inside`           | `inside` | mean 0.177            | 50/50  |

Surface objects sit on the membrane to machine epsilon and stay there while
diffusing. Inside objects stay strictly inside -- 150/150 across both types,
with none escaping.

The low `inside` counts for the *surface* rows are not failures: those objects
are at distance 0, where the inside/outside sign is a coin flip. The distance
column is the meaningful one for `confine = on`.

## Surface nucleation run

`config-nucleate.cym` is the dynamic counterpart: singles and couples diffuse on
the membrane and nucleate filaments where they sit, so every filament is born AT
the surface rather than in the lumen. `nucleate_in_plane = 1` orients each new
filament tangent to the membrane -- the nucleator asks the Space for
`normalToEdge(pos)`, which `SpaceMem3DG` answers from the same pseudonormal used
by `inside`.

```bash
coupled_driver config-nucleate.cym ../icosphere.obj 4000 0.0035 10 1 100 frames-nucleate 1 8.22e-6 100.0 100.0 0.98 0.0 0
```

Confinement range 0.0035, strength 10 (driver args 4 and 5).

Filament count grows 0 -> 13 -> 20 -> 30 -> 38 over 41 exported frames, and
`|F_fib|` equals `|F_mem|` at every step, so the steric exchange is
action-reaction exact.

Distance to the **deformed** membrane at frame 40 (exact point-to-triangle):

| class                             |  n  | mean d   | note |
|-----------------------------------|-----|----------|------|
| `single surfaceNucleator free`    | 36  | 3.1e-07  | on the surface |
| `single surfaceDiffuser free`     | 99  | 3.0e-07  | on the surface |
| `couple surfaceNucleatorPair free`| 25  | 2.8e-07  | on the surface |
| `single surfaceNucleator bound`   | 24  | 1.2e-01  | carried by its filament |
| `couple innerCross free`          | 49  | 2.3e-01  | inside, as configured |

The 3e-07 floor is the `%g` precision of the exported coordinates, not drift --
free surface objects track the membrane as it deforms.

`bound` objects leaving the surface is correct, not a confinement failure: once a
hand attaches, the object's position is the hand's position along the filament,
and the filament runs into the lumen. Confinement governs the diffusing state.

## Blender

`points_%04d.txt` is written next to the membrane and fiber frames and is picked
up by `blender/import_trajectory.py` with no extra setup -- `_TEST` points at
`frames-nucleate`. Each block becomes one vertex-instanced mesh: hundreds of
real sphere objects per frame would make scrubbing unusable, so one object per
class carries all its points and instances a shared sphere mesh on every vertex.

Blocks are keyed by `<kind> <name> <state>`, so binding state gets its own
material: grey translucent = `free`, orange = `bound`, red = `bridge` (a couple
holding two filaments). The state comes free from Cytosim's list structure --
`SingleSet` already keeps free and attached apart, `CoupleSet` keeps all four
combinations -- so no per-object queries are needed.

## Files

- `config.cym` -- static placement check: 200 pinned + 100 diffusing surface
  singles, 100 inside singles, 150 surface couples, 50 inside couples
- `config-nucleate.cym` -- surface diffusion + nucleation (the run above).
  `set fiber actin` must precede `set hand seeder`: the hand's `nucleate` key is
  resolved to a fiber class at parse time
- `config-bound.cym` -- variant with pre-seeded filaments, to exercise the
  `bound` point classes
- `verify_placement.cc` -- the placement/confinement measurement
