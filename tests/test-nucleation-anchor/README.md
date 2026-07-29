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
