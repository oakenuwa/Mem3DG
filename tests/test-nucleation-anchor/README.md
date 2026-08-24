# Plus-end-held nucleation from a 0.7 µm patch — filament-filament steric on/off

100 non-diffusing anchors in a single 0.7 µm patch of a Mem3DG icosphere. Each
holds and **tracks the growing plus end** of one filament while the filament
body extends into the lumen — the formin/filopodium geometry. Two runs,
identical except for Cytosim's fiber-fiber steric engine.

| | `config-steric.cym` | `config-no-steric.cym` |
|---|---|---|
| `set simul { steric }` | `1, 500` | `0` |
| everything else | identical | identical |

The functional diff is exactly those lines. Worth re-checking with `diff` after
any edit — an earlier attempt at this pair silently failed to apply the change
and produced two identical runs that looked like a finished experiment.

## Running

```bash
coupled_driver config-steric.cym    ../icosphere.obj 12000 0.0035 10 1 400 frames-steric    1 8.22e-6 100.0 100.0 0.98 0.0 0
coupled_driver config-no-steric.cym ../icosphere.obj 12000 0.0035 10 1 400 frames-no-steric 1 8.22e-6 100.0 100.0 0.98 0.0 0
```

Filament-membrane steric — range 0.0035, strength 10, one-sided wall (mode 1) —
is applied by the driver through the Mem3DG C ABI and is **on in both runs**;
`simul:steric` does not reach it.

## The plus-end configuration

```
track_end = plus_end
hold_end  = plus_end
nucleate  = 200, actin, (length=0.05; plus_end=grow;)
```

`track_end` is what makes the hand *follow* the growing tip. With `hold_end`
alone the hand stays at the abscissa it first gripped, which stops being the plus
end as soon as the fiber elongates — the anchor would end up holding an interior
point. `NucleatorProp::read` sets `hold_end = track_end` when `track_end` is
given and `complete()` rejects the pair if they disagree; both are written out to
keep the intent explicit.

### A Cytosim fix this required

`nucleate_normal` had to be made `hold_end`-aware. In `Nucleator::makeFiber` the
fiber is rotated so its minus→plus direction aligns with `dir`, then translated
so the **held** end lands on the nucleator. So:

| `hold_end` | body extends along | `dir` must be |
|---|---|---|
| `minus_end` | `+dir` | inward, `-out` |
| `plus_end` | `-dir` | outward, `+out` |

The generic `if (hold_end == PLUS_END) dir.negate()` earlier in the function does
this, but the `nucleate_normal` block runs *after* it and assigned `dir = -out`
unconditionally, discarding the flip. Fixed to
`dir = (hold_end == PLUS_END) ? out : -out`.

This fails silently, which is why it is worth stating: with the wrong sign the
filament is simply built on the outside of the membrane and immediately fights
the confinement, rather than erroring.

## Verification

`verify_plusend.py`. Measured at the **nucleation frame** (0000):

| | steric | no-steric |
|---|---|---|
| body angle to inward normal at launch | 0.1° | 1.75°, 0.78° |
| \|plus end\| at launch | 0.8001 | 0.7978 |

The launch is normal to the surface to within ~2°.

**Check frame 0000, not 0001.** `saveEvery` is 400 steps, so by frame 0001 the
filaments have had 0.02 s to reorient under load and read ~57° even when the
launch was exact. That is a measurement artifact of the export cadence, not a
setup error.

The mechanism holds for the whole run — plus end pinned at the membrane, minus
end travelling into the lumen:

| frame | mean length | \|plus end\| (held) | \|minus end\| (free) | body angle |
|---|---|---|---|---|
| 0001 | 0.0728 | 0.7985 | 0.7649 | 57.4° |
| 0010 | 0.3281 | 0.7985 | 0.6203 | 43.8° |
| 0030 | 0.8738 | 0.7984 (min 0.7953, max 0.8027) | 0.5730 (min 0.1312) | 39.7° |

(steric run; the no-steric run is the same to within a few percent). The plus end
never leaves the 0.795–0.803 shell, so the anchor really is holding the growing
tip; the minus end reaches 0.13 from the centre.

The 40–57° tilt is real reorientation, not a launch error: a filament pinned at
its tip with a free body, in a converging bundle of 100, tilts under load.

> **Resolution caveat, added after the deformation runs below.** The numbers in
> this section were computed on the un-remeshed 642-vertex mesh, which
> **under-resolves the filament-membrane contact by ~20x in force**. The
> steric on/off comparison is still valid (both sides used the same mesh), but
> the absolute forces and deformation magnitudes here are understated. See
> "Resolution: the coarse mesh was wrong" below.

## Result — does fiber-fiber steric change anything?

At the last frame, over all 4950 filament pairs. Steric radius 0.01 µm, so two
filaments repel once their axes come within 0.020 µm:

| | steric ON | steric OFF |
|---|---|---|
| pairs closer than 0.020 (the steric range) | **41** (0.8%) | **145** (2.9%) |
| pairs closer than 0.010 | 6 (0.1%) | 33 (0.7%) |
| 1st percentile of pair distance | **0.02108** | 0.01214 |
| 5th percentile | 0.03674 | 0.02508 |
| median pair distance | 0.22955 | 0.18410 |
| minimum pair distance | 0.00139 | 0.00091 |

**3.5× fewer overlapping pairs**, 5.5× fewer deep overlaps, and the 1st
percentile now sits *above* the 0.020 interaction range — 99% of pairs are
properly excluded. Median spacing is 25% larger, so the converging bundle is
genuinely held more open.

Steric also keeps the bundle more radial: body angle 43.8° vs 55.4° at frame 10
and 39.7° vs 42.1° at frame 30. Mutual repulsion resists filaments crossing, so
they stay closer to the normal they launched along.

**The minimum distance is not the right statistic.** It is 0.0014 with steric ON,
far inside the 0.020 range, which reads like the steric is not working. It is one
extreme of 4950 pairs, and a finite-stiffness harmonic wall under load always
leaves a penetrating tail. The percentile rows are the honest measure.

Independent confirmation the toggle took effect: the steric run's startup prints
a `LocusGrid` line (Cytosim's fiber-fiber steric grid); the no-steric run does
not.

## Comparison with minus-end holding

`frames-steric-minus-end/` and `frames-no-steric-minus-end/` are the same
experiment with `hold_end = minus_end` (Arp2/3-like: the capped end sits at the
membrane, the growing tip travels inward). Kept for contrast:

| pairs closer than 0.020 | steric ON | steric OFF | ratio |
|---|---|---|---|
| plus end held (this run) | 41 | 145 | 3.5× |
| minus end held | 61 | 144 | 2.4× |

The steric effect is cleaner in the plus-end geometry. That is consistent with
where the filaments actually crowd: with the tip pinned at the membrane the
bodies fan into the lumen and stay roughly radial, whereas minus-end holding
sends the growing tips through the convergence focus and then along the far wall,
where they pile up regardless of the steric setting.

## Files

- `config-steric.cym` / `config-no-steric.cym` — the pair
- `verify_plusend.py` — launch angle, held/free end tracking, steric contrast
- `verify_patch.py` — the earlier verifier; its angle column assumes the
  **minus** end is held, so use `verify_plusend.py` for these runs

`verify_patch.py` also had to learn that an attached anchor exports **two**
points — `anchor` (the pinned foot) and `hand` (where the hand sits on the
filament). It originally summed both and reported 101/103 nucleators for 100.

## Viewing

`blender/import_trajectory.py` picks up `membrane_%04d.obj`, `fibers_%04d.txt`
and `points_%04d.txt` from either directory unchanged.


---

# High filament-membrane repulsion, with remeshing

Same two configs, driven hard enough to deform the membrane, with adaptive
remeshing on:

```bash
coupled_driver config-steric.cym    ../icosphere.obj 12000 0.05 100 1 400 frames-steric-deform    1 8.22e-6 100 100 0.98 0.0 1 0.12 100
coupled_driver config-no-steric.cym ../icosphere.obj 12000 0.05 100 1 400 frames-no-steric-deform 1 8.22e-6 100 100 0.98 0.0 1 0.12 100
```

Shell 0.05 (14x the earlier 0.0035), stiffness 100, remeshing every 100 steps at
`maxEdge = 0.12`.

## Calibration: stiffness is the wrong lever

Raising the steric *stiffness* **reduces** the force it delivers, because a
stiffer wall ejects the filament tip out of the interaction shell and it then
feels nothing. Contact counts over a 1500-step run, shell fixed at 0.0035:

| stiffness | contacts, start -> end | \|F_mem\| end |
|---|---|---|
| 10 | 118 -> 87 | 1.4 |
| 200 | 36 -> 5 | 3.3 |
| 1000 | 8 -> **0** | **0** |

At fixed shell the stiffness barely matters (shell 0.02 gives bulge +0.031 at
both k=50 and k=200). The lever is the **shell width**, which is what keeps the
point in range while it is being pushed:

| shell | k | bulge | rms deviation | contacts |
|---|---|---|---|---|
| 0.0035 | 10 | +0.0043 | 0.0013 | 87 |
| 0.02 | 50 | +0.0312 | 0.0061 | 23 |
| 0.05 | 100 | +0.0631 | 0.0136 | 24 |
| 0.08 | 150 | +0.0942 | 0.0209 | 27 |
| 0.10 | 100 | +0.1140 | 0.0260 | 36 |

Note the shell is a modelling choice, not a molecular radius: it is the
clearance the one-sided wall enforces between filament and membrane.

## Resolution: the coarse mesh was wrong

Enabling remeshing changed `|F_mem|` from 12.4 to 327.8 at the same step with the
same parameters and an essentially unchanged filament configuration (397 vs 415
filament points). That is not remeshing corrupting the physics — the coarse mesh
was under-resolving it. Holding physics fixed (shell 0.05, k 100, 800 steps) and
varying only the mesh:

| vertices | contacts | \|F_mem\| |
|---|---|---|
| 642 | 27 | 8.8 |
| 1869 | 124 | 179.4 |
| 2059 | 129 | 182.6 |
| 5661 | 128 | 175.2 |
| 8156 | 125 | 172.1 |

Converged from 1869 onward — a further 4.4x refinement moves it 4%. The
642-vertex mesh is the outlier, and the reason is countable:

| mesh | faces | faces inside the 0.7 um patch | faces per filament tip |
|---|---|---|---|
| initial icosphere | 1280 | 62 | **0.62** |
| remeshed, maxEdge 0.12 | 3734 | 2109 | 21.1 |
| remeshed, maxEdge 0.06 | 11318 | 2585 | 25.9 |

**The original mesh has fewer triangles in the patch than there are filaments**,
so per-tip contact is not representable at all. Remeshing puts 56% of its faces
into the 5% of surface area where the deformation is.

Consequence: steric parameters calibrated at 1280 faces do not transfer to a
remeshed run, and every earlier run in this directory had the filament-membrane
force about 20x too low.

## Result: the membrane deforms, and the deformation is localized

642 -> 2923 (steric) / 2859 (no-steric) vertices. Radius by region at the last
frame, initial radius 0.80:

| region | vertices | mean r | max r |
|---|---|---|---|
| patch, z > 0.719 | 1299 | **0.8361** | **0.8520** |
| equator | 1553 | 0.7915 | 0.8187 |
| antipode, z < −0.719 | 71 | 0.7898 | 0.7929 |

The patch bulges out ~0.036 µm on average and 0.052 µm at the peak, while the
equator and the far pole pull in ~0.009 µm. That is volume-conserving
redistribution under `Kv = 100`, and the protrusion is exactly where the
nucleators are. Deformation rms is 0.0235 against 0.0034 for the earlier
low-repulsion run — a 7x increase.

The load relaxes as the membrane accommodates the tips: `|F_mem|` 46.8 at step
4000, 19.9 at 8000, 12.2 at 12000.

## Result: filament-filament steric does NOT change the membrane shape here

The filament bodies are clearly separated by the steric engine:

| | steric ON | steric OFF |
|---|---|---|
| pairs closer than 0.020 | **37** (0.7%) | **150** (3.0%) |
| pairs closer than 0.010 | 7 | 34 |
| 1st percentile of pair distance | 0.02152 | 0.01185 |
| median pair distance | 0.22616 | 0.15695 |

4.1x fewer overlaps, median spacing 44% larger. And yet the membrane is
identical between the two runs:

| | steric ON | steric OFF |
|---|---|---|
| patch mean r | 0.8361 | 0.8364 |
| patch max r | 0.8520 | 0.8520 |
| rms deviation | 0.0235 | 0.0234 |

**Why, and why it is not a bug.** `activity = anchor` pins each foot in the lab
frame with `stiffness = 2000`, and the plus end is held there. The membrane is
loaded by those pinned tips pushing outward through the steric wall, so the tip
positions — and therefore the load — are set by the *anchor placement*, not by
what the filament bodies do to each other. Fiber-fiber steric rearranges the
bodies in the lumen and does not move the tips.

So this is a real negative result, not a failed run: in a plus-end-anchored
geometry, filament-filament steric cannot influence membrane shape. To make it
influence the membrane you need the loading points free to be repositioned by
neighbour interactions — hold the **minus** end instead (growing tips free to
press the far membrane), soften the anchor, or drop the anchors entirely.

## Files

- `frames-steric-deform/`, `frames-no-steric-deform/` — the runs above
- `frames-steric/`, `frames-no-steric/` — plus-end held, low repulsion, no remeshing
- `frames-steric-minus-end/`, `frames-no-steric-minus-end/` — minus-end held
