# Nucleation from a membrane patch: filament-filament steric on vs off

100 nucleators pinned to a 1 µm patch of a Mem3DG icosphere, no starting
filaments. Each nucleator fires once and holds the filament it makes, so the run
develops exactly 100 filaments anchored at the membrane and growing into the
lumen, confined by filament-membrane steric.

Two configs, differing in **one line**:

| | `config-steric.cym` | `config-no-steric.cym` |
|---|---|---|
| `set simul { steric }` | `1, 500` | `0` |

Everything else is byte-identical — verify with

```bash
diff <(grep -vE "^\s*%|^\s*$" config-steric.cym) <(grep -vE "^\s*%|^\s*$" config-no-steric.cym)
```

Filament-**membrane** steric is active in both: the driver applies it through
the Mem3DG C ABI (one-sided wall, `mode = 1`), independent of Cytosim's
`simul:steric`, which governs only fiber-fiber.

## Running

```bash
coupled_driver config-steric.cym ../icosphere.obj 12000 0.0035 10 1 400 frames-steric 1 8.22e-6 100.0 100.0 0.98 0.0 0
```

Confinement range 0.0035 µm, strength 10. Swap the config and output directory
for the other arm.

## Setup, and four things that had to be fixed to get it right

**Patch placement.** The nucleators sit on `{z > 0.6245}`, which is the 1 µm
diameter cap: rim radius 0.5 on R = 0.8 gives z = √(0.8² − 0.5²) = 0.62450.
Imposed with

```
position  = surface
placement = surface, membrane, (Z > 0.6245)
```

The obvious alternative, `position = cap 0.8 0.1755`, **does not work**. It
generates points at radius exactly 0.8, but the mesh is a polyhedron *inscribed*
in that sphere — measured inradius 0.796377 — so every generated point lies
outside the membrane and Cytosim rejects all 1024 placement trials. Sampling
`surface` and filtering by z puts them on the real triangulated surface to
machine precision.

**`single:stiffness` defaults to 0.** This was the significant one. With the
default, a nucleator's hand is attached in Cytosim's bookkeeping but **no force
holds the filament to its anchor**. Measured: filament minus ends drifted up to
**0.735 µm** from their nucleators — further than the filament's own length — so
the render showed a cluster of anchors with filaments scattered anywhere. It also
suppressed the physics: membrane contacts ~15 and `|F_mem|` ~0.6. Setting
`stiffness = 2000` brought the link stretch to ≤ 3.9 nm and contacts to ~110
with `|F_mem|` ~2.9.

**An attached single is two points, not one.** `posFoot()` is the anchor,
`posHand()` is the grip on the filament; they coincide only while free. The
exporter now emits both, as classes `anchor` and `hand`, so the viewer shows the
membrane anchor and the filament grip separately and the link stretch is visible.

**A single membrane got no material in Blender.** `_membrane_material` was only
applied when more than one membrane sequence was present, so a one-membrane run
imported an untextured opaque mesh that simply hid every filament inside it —
which is why the membrane appeared "missing". It is now always assigned a
material, translucent by default (`MEMBRANE_COLOR`, alpha 0.25, `blend_method =
BLEND`, backface culling off).

## Verification

Measured on the exported frames with an independent point-to-triangle routine:

| frame | anchor → membrane (max) | anchor → hand (max) | hand → minus end (max) | filaments |
|---|---|---|---|---|
| 0000 | 1.95e-06 | 0.0 | 0.0 | 2 |
| 0001 | 1.62e-03 | 3.51e-03 | 0.0 | 99 |
| 0005 | 4.25e-03 | 3.72e-03 | 0.0 | 100 |
| 0010 | 4.97e-03 | 3.94e-03 | 0.0 | 100 |

- **Filaments are located with the nucleators.** `hand → minus end` is exactly 0
  — the hand grips the filament's minus end by construction — and the anchor is
  within 3.9 nm of it.
- **Filaments stay bound.** All 100 remain attached for the whole run; the link
  never stretches past 4 nm.
- **Nucleators stay on the membrane**, to within 5 nm at frame 10.

That last number is small but **grows**, and the reason is worth knowing:
`activity = fixed` pins a nucleator in the **lab frame**, not to the membrane. As
the membrane deforms under the filaments it moves out from under them. At 5 nm on
an 800 nm sphere (0.6%) this is negligible, but it is drift, not noise, so a much
longer or more strongly deforming run would need a membrane-anchored single
instead.

## Known limitation: filaments are not normal to the membrane

`nucleate_normal = 1` sets the **launch** direction and nothing maintains it. A
newly nucleated filament is ~0.09 µm long and anchored at one end, giving
D_rot ≈ 193 rad²/s at the default `kT = 0.0042` — about 159° of rotation per
0.02 s output interval. Measured:

| | launch (frame 0) | frame 1 | frame 10 |
|---|---|---|---|
| angle to inward normal, mean | 1.7° | 55.4° | 47.2° |

(90° would be fully isotropic.) The stiff anchor clamp does **not** fix this: it
pins the minus-end position, not the direction.

This run deliberately keeps thermal noise, so the array is normal only at
nucleation. `kT = 0` in `set simul` would make the filaments stay normal and
would also make the two arms a controlled comparison rather than two noise
realizations — worth doing if the normal geometry matters for the question being
asked.

## Viewing

```
frames-<case>/membrane_%04d.obj
frames-<case>/fibers_%04d.txt
frames-<case>/points_%04d.txt
```

Point `blender/import_trajectory.py` at either directory. The membrane renders
translucent so the filaments inside are visible; anchors are blue, filament grips
orange.
