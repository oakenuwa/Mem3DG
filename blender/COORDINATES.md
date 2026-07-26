# Shared coordinate convention

The whole point of a shared convention is that the membrane and the filaments
land in the **same place** in Blender without anyone hand-tuning offsets. This
file is the single source of truth for that. If the two ever fail to line up,
the bug is here.

## The rule

**Both simulators export raw world coordinates, in micrometres, about a common
origin, with no transform applied on export.** All registration happens by
sharing one origin and one scale — never by shifting one object to match the
other.

| | Mem3DG (membrane) | Cytosim (filaments) |
|---|---|---|
| Exporter | `coupled_driver` → `saveMembraneFrame()` | `export_fibers` tool |
| Format | Wavefront **OBJ** (`membrane_%04d.obj`) | text polylines (`fibers_%04d.txt`) |
| Units | micrometres | micrometres |
| Origin | mesh origin as loaded / evolved | Cytosim world origin (space centre) |
| Transform on export | **none** | **none** |
| Z in 2D | n/a (Mem3DG is 3D) | written as `0.0` |

The Blender importer applies exactly **one** transform: a single
`WORLD_SCALE` multiplier, identical for both objects. That rescales the scene
for convenient viewing; it does not move one relative to the other.

## Why raw, why no per-export offset

Every membrane force in the coupling already goes through minimum-image /
projection math that assumes both bodies live in one coordinate frame — the
filament forces and membrane forces are computed from `pos_fiber - pos_vertex`
directly (`contactForces()` in `coupled_driver.cc`). So the two are *already*
in a common frame at run time. Exporting raw preserves that frame. Any offset
applied at export would silently desynchronise the picture from the physics.

## Frame index = time

`membrane_000K.obj` and `fibers_000K.txt` are the **same simulated instant**.
This only holds if the two exporters use the same cadence:

- membrane: `coupled_driver ... <saveEvery>` (7th CLI arg) — a frame every
  `saveEvery` coupled steps.
- fibers: `export_fibers period=N` — every `N`th trajectory frame.

Choose these so the emitted counts match. The importer pairs them by export
index, not by wall-clock time, so a mismatch shows up as the membrane and
filaments drifting out of sync as the animation plays.

## Registration checklist (when they don't line up)

1. **Same origin?** Cytosim's origin is the centre of its `space`. Mem3DG's is
   wherever the input mesh sits. If the membrane mesh was authored off-origin,
   translate the *input mesh* once, before the run — not on export.
2. **Same scale?** Both must be in the same length unit. If one config used nm
   and the other µm, fix the config; do not paper over it with `WORLD_SCALE`.
3. **Same handedness?** Both are right-handed XYZ. No axis swap is applied. If
   Blender looks mirrored, check you didn't add one.
4. **2D vs 3D Cytosim.** A `DIM=2` build writes `z=0`; the membrane is fully 3D.
   They still share X and Y, so they overlay in the Z=0 plane — expected.

## What is deliberately NOT done here

- No per-frame recentering (would make the membrane appear to jitter).
- No bounding-box normalisation (would break the shared scale).
- No unit auto-detection (silent and error-prone; units are a config decision).
