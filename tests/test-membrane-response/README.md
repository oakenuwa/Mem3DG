# Membrane response to local deformation

Diagnoses why filament protrusions came out as **needles** — a few vertices
dragged far out with hugely stretched triangles — instead of smooth bulges, and
fixes it.

## Contents

| Path | What |
|---|---|
| `membrane_response.cc` | pulls one vertex out, reports the force on each topological ring around it |
| `response.txt` | its output |
| `stretch_metric.py` | edge stretch vs the initial mesh, per frame |
| `spread.py` | how displacement is distributed with distance from the filaments |
| `config-icosphere.cym` | the Cytosim side |
| `frames/` | the fixed run: 1001 frames |
| `final_run.txt` | verification of that run |

## The question

Does a local deformation propagate, i.e. do neighbouring vertices respond?

`membrane_response.cc` answers it directly: displace vertex 0 outward by 0.15,
then read the membrane's own internal forcing and average |force| over each
topological ring around that vertex.

```
                                       ring0     ring1     ring2     ring3     ring4
bending only                         7.997e-04 3.446e-04 4.266e-05 6.267e-08 1.988e-07
+ tension Ksg=100, At=A0             1.055e-01 1.415e-02 1.111e-02 1.109e-02 1.196e-02
+ tension Ksg=100, PRE-TENSED .98    6.836e-01 9.260e-02 7.219e-02 7.232e-02 7.797e-02
+ Ksg=100 pre + Kv=100               6.843e-01 9.252e-02 7.300e-02 7.314e-02 7.885e-02
+ ... + EDGE SPRINGS Kse=0.1         3.423e+00 7.560e-01 7.300e-02 7.314e-02 7.885e-02
+ ... + EDGE SPRINGS Kse=1.0         2.807e+01 6.772e+00 7.300e-02 7.314e-02 7.885e-02
```

Reading it:

- **Bending alone is local.** Force falls four orders of magnitude by ring 3.
  Nothing beyond the immediate neighbourhood knows the deformation happened.
- **Tension is global but weak, and only when pre-tensed.** With `At = A0` the
  membrane is tensionless at rest and ring 0 sees only 0.105. Setting
  `At = 0.98·A0` raises that 6.5x. Note rings 2-4 all sit at ~0.073: tension is
  a *uniform* far-field response, not a decaying one.
- **Edge springs act where it matters.** `Kse = 0.1` raises ring 0 to 3.4 and
  ring 1 to 0.76 — an 8x jump in the FIRST ring specifically. That is the term
  that makes the immediate neighbours share the load.

## Why tension was not enough

Tension penalises **total** area. The mesh can therefore concentrate all its
stretch into a handful of edges at essentially no cost. Measured on the run that
produced the needles:

```
max edge stretch 4.49x     mean edge stretch 0.999
```

A few edges stretched 4.5x while the mean was unchanged — so tension registered
nothing wrong. `Kse` penalises each edge **individually** against its reference
length, which is what forces the deformation to be shared.

## Choosing Kse

Stiffer springs suppress the needle but also suppress the deformation:

| Kse | max edge stretch | max displacement | r_max |
|---|---|---|---|
| 0 | **4.08** | 0.377 | 1.105 |
| 0.03 | 2.00 | 0.207 | 0.967 |
| **0.1** | **1.44** | **0.135** | 0.923 |
| 0.3 | 1.18 | 0.088 | 0.886 |
| 1.0 | 1.04 | 0.052 | 0.852 |

`Kse = 0.1` keeps stretch under ~1.5 while still permitting a visible bulge.
Above ~0.3 the mesh is effectively rigid.

## The fixed run

```
frame  maxStretch  maxDisp   r_max   escaped  closed
    0      1.00    0.0000   0.800    0/12    yes
  500      1.36    0.1313   0.931    0/38    yes
 1000      1.60    0.1923   0.987    0/60    yes

total escaped: 0
```

Deformation is now smooth (stretch 1.6 rather than 4.5), confinement still
holds, and the mesh stays closed.

Regenerate with:

```
coupled_driver config-icosphere.cym ../icosphere.obj 40000 0.05 100.0 1 40 \
    frames 1 8.22e-6 100.0 100.0 0.98 0.1
                          Kb      Ksg   Kv  areaFrac Kse
```

## Two ways to fix it — pick one

`frames/` uses **edge springs**, `frames-remesh/` uses **adaptive remeshing**.
They are alternatives, not complements: `Kse` measures every edge against a
*reference* mesh, and remeshing invalidates that reference for every edge it
creates. The driver warns if both are on.

### (a) Edge springs — `Kse` (frames/)

Resists the stretch. Cheap, no topology change, vertex count fixed. But it
resists the *tube*, not just the artefact: `r_max` only reaches 0.99 because the
mesh is stiffer everywhere.

### (b) Adaptive remeshing — `remesh=1` (frames-remesh/)

Resolves the stretch by adding material where it is needed:

```
frame   nVerts  nFaces   maxEdge  meanEdge   r_max   closed
    0     642    1280   0.1317   0.1206   0.800    yes
  250    1097    2190   0.1410   0.0907   0.806    yes
  500    1320    2636   0.1414   0.0795   0.866    yes
 1000    1905    3806   0.1408   0.0625   0.914    yes
```

The mesh nearly **triples** (642 -> 1905 vertices) while `maxEdge` stays *capped*
near 0.141 instead of drifting to 0.22. Refinement is local: the fraction of
vertices within 0.15 of a filament grows 0% -> 53%, i.e. new triangles appear
exactly where the filaments push. Confinement holds (0 escapes) and the mesh
stays closed.

### A trap in Mem3DG's mutator

`splitLong` does **not** use `maximumEdgeLength`, despite the name. Its test is
*relative*:

```cpp
edgeLength(e) > edgeLength(next) + edgeLength(twin.next)
```

so an edge must exceed the sum of two neighbours -- roughly a doubling against
the local mesh, which a smoothly stretching tube never reaches. Setting
`maximumEdgeLength` alone produces **no splits at all** (verified: vertex count
stayed at 642 while edges reached 0.217 against a 0.1275 threshold).
`maximumEdgeLength` only guards against *collapsing* long edges.

The absolute criterion is `splitLarge`, on face area:

```cpp
faceArea(he) + faceArea(twin) > 4 * maximumFaceArea
```

`mem3dg_system_set_remeshing()` therefore derives
`maximumFaceArea = 0.2165 * max_edge^2`, which places the threshold at edges of
length `max_edge` (two equilateral triangles of side L sum to 0.866 L^2).

### Usage

```
                                                    Kse  remesh maxEdge every
coupled_driver ... frames-remesh 1 8.22e-6 100 100 0.98 0.0    1    0     50
```

`maxEdge = 0` auto-derives ~1.5x the initial mean edge length.

## Caveat

Remeshing changes the vertex and face count, so any buffer sized from
`num_vertices()` must be re-queried after a step that mutates. The driver does
this; the internal external-force buffer is resized and **cleared** automatically
(the old per-vertex forces no longer refer to anything).
