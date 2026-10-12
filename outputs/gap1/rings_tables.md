# Section 6.1-6.5 on the Pohl22 ring background

Produced by `analysis/gap1/robustness_scan.py --rings` and
`analysis/gap1/fit_rings.py`. Full records under `outputs/gap1/rsp_*`,
`outputs/gap1/rings_*`, not committed. Every fit converges: summed lnL gap
bound 6e-7 over the eight bins, so every TS below is good to ~1e-6.

## Table R1. The freedom ladder, Mask B + PS, 1-10 GeV, observed TS

| diffuse model | components | gNFW 1.2 | HESTIA G1.1 0deg | VVV |
|---|---|---|---|---|
| `gll_iem_v07` | 1 | 120.7 | 248.0 | 203.1 |
| `rings_species` | 5 | 526.2 | 622.6 | 2151.7 |
| `rings_free` | 16 | 125.9 | 75.9 | 794.0 |
| `rings_free_NB` | 17 | 21.8 | 12.1 | 570.5 |
| `rings_free_NB_BB` | 18 | 0.0 | -0.0 | 124.2 |

## Table R2. Admissibility is per template

Injected TS and boosted gain for a signal added at the flux the
`gll_iem_v07` fit finds for that same template. Threshold 25 on both, set
in `configs/diffuse_freedom.yaml` before any run.

| diffuse model | VVV | HESTIA | gNFW 1.2 |
|---|---|---|---|
| `gll_iem_v07` | 184.3 / +585 **yes** | 240.0 / +684 **yes** | 109.3 / +320 **yes** |
| `rings_species` | 255.2 / +1472 **yes** | 219.3 / +746 **yes** | 134.7 / +478 **yes** |
| `rings_free` | 82.9 / +505 **yes** | 33.8 / +137 **yes** | 23.6 / +135 no |
| `rings_free_NB` | 37.3 / +343 **yes** | 15.8 / +37 no | 21.4 / +26 no |
| `rings_free_NB_BB` | 8.9 / +49 no | 14.3 / -0 no | 10.2 / -0 no |

`rings_free` recovers an injected stellar bulge at TS 83 and an injected
HESTIA halo at 34, but an injected gNFW at only 23.6. The same
configuration is admissible for two template classes and not the third,
which is why the control is now run per template.

## Table R3. Per-template controls at one level, Mask B + PS, 1-10 GeV

| injected template | rings_species | rings_free | Mask D, rings_free |
|---|---|---|---|
| `vvv` | 255.2 **yes** | 82.9 **yes** | 14.9 no |
| `f98` | 223.1 **yes** | 47.6 **yes** | 4.3 no |
| `hestia_G1.1_ang0` | 219.3 **yes** | 33.8 **yes** | 9.7 no |
| `hestia_G3.2_ang90` | 139.9 **yes** | 25.5 **yes** | 2.9 no |
| `gnfw_g1.2` | 134.7 **yes** | 23.6 no | 2.4 no |
| `gnfw_g1.4` | - | 36.4 **yes** | 6.4 no |

Nothing is admissible under Mask D once the background is ring-decomposed,
and Mask A keeps 8.7% of the ROI and has no power at all. Mask B + PS is the
only usable configuration.

## Table R4. Section 6.2 redone at the admissible level

`rings_species`, Mask B + PS, 1-10 GeV. Every template in this table passes
its own injection control.

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 2151.7 | -0.0 |
| `f98` | 1761.2 | +390.5 |
| `hestia_G3.2_ang90` | 1263.2 | +888.5 |
| `hestia_G2.2_ang90` | 1084.3 | +1067.4 |
| `hestia_G2.2_ang0` | 1056.2 | +1095.5 |
| `hestia_G2.2_ang60` | 1033.4 | +1118.2 |
| `hestia_G2.2_ang30` | 977.1 | +1174.6 |
| `hestia_G1.1_ang90` | 947.1 | +1204.6 |
| ... | | |
| `gnfw_g0.9_th1` | 0.0 | +2151.7 |

## Table R5. Section 6.1, the axis-ratio correlation, on each background

| background | Table 1 q (published recipe) | per-angle 5% | per-angle 50% |
|---|---|---|---|
| `gll_iem_v07`, Mask B + PS (as published) | +0.656 (p 0.00051) | +0.217 (p 0.31) | -0.251 (p 0.24) |
| `rings_species`, Mask B + PS | -0.190 (p 0.37) | +0.253 (p 0.23) | +0.363 (p 0.081) |
| `rings_free`, Mask B + PS | -0.056 (p 0.79) | +0.400 (p 0.053) | +0.383 (p 0.064) |
| `rings_species`, Mask D + PS | +0.063 (p 0.77) | +0.457 (p 0.025) | +0.157 (p 0.47) |

## Table R6. Section 6.4, energy dependence, rings_species, Mask B + PS

| band GeV | VVV | F98 | best HESTIA | gNFW 1.2 |
|---|---|---|---|---|
| 0.5-1 | 0.0 | 0.0 | hestia_G2.2_ang0 0.0 | 0.0 |
| 1-3 | 1628.9 | 1279.4 | hestia_G3.2_ang90 968.1 | 368.7 |
| 3-10 | 736.3 | 651.5 | hestia_G3.2_ang90 380.8 | 152.4 |
| 10-50 | 20.9 | 14.8 | hestia_G2.2_ang60 0.7 | 0.0 |

The 0.5-1 GeV band detects nothing at all on this background, for any
template, and no injection control was run for it; it is not a null, it is
an empty cell.

## Table R7. Section 6.5, triaxial ordering, rings_species, Mask B + PS

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 2151.7 | -0.0 |
| `triaxial_vert0.6_plane0.8_flip` | 764.3 | +1387.4 |
| `triaxial_vert0.6_plane0.8_std` | 692.5 | +1459.2 |
| `triaxial_vert0.8_plane0.6_flip` | 643.3 | +1508.3 |
| `hestia_G1.1_ang0` | 622.6 | +1529.1 |
| `gnfw_g1.2` | 526.2 | +1625.5 |
| `triaxial_vert0.8_plane0.6_std` | 510.0 | +1641.7 |
| `hestia_G2.1_ang90` | 329.9 | +1821.8 |
| `gnfw_g1.1` | 318.6 | +1833.1 |

Vertical flattening still beats in-plane flattening, 764 against 643, the
same ordering as on the single-template background. Both are a long way
behind the stellar bulge.

## Table R8. What flux was injected, and per-template recovery

The flux injected for the admissibility test in Tables R2-R3 was **the flux
that template itself fits against the single-template `gll_iem_v07`
baseline** - not its own fitted flux at the level under test, and not VVV's.
Those fluxes: VVV 6.3931e-08, gNFW 1.2 6.0292e-08, HESTIA 9.5463e-08
ph cm^-2 s^-1, already within a factor 1.6 of one another.

Repeated with **every template injected at one common flux and one common
per-bin spectrum**, VVV's 6.3931e-08, so that differences in recovered TS
are differences in sensitivity to morphology and nothing else.
`rings_species`, Mask B + PS, 1-10 GeV.

| injected template | observed TS | injected-only TS | boosted gain | admissible |
|---|---|---|---|---|
| `vvv` | 2151.7 | **255.2** | +1472 | yes |
| `f98` | 1761.2 | **355.6** | +1637 | yes |
| `hestia_G3.2_ang90` | 1263.2 | **81.8** | +566 | yes |
| `hestia_G1.1_ang0` | 622.6 | **65.6** | +426 | yes |
| `gnfw_g1.4` | 799.7 | **81.8** | +480 | yes |
| `gnfw_g1.2` | 526.2 | **91.8** | +413 | yes |

At identical flux the stellar bulges recover TS 255-356 and the dark-matter
templates 66-92. **The pipeline is about three times more sensitive in TS to
a bulge morphology than to a dark-matter one** - a factor 1.7 in flux, since
TS goes as flux squared. Every template still clears the threshold of 25 by a
wide margin, so an injected dark-matter signal of this brightness would have
been detected; but the dAIC of +1626 between VVV and gNFW on the real data
cannot be read as if the two were measured at equal sensitivity. Table R9
settles what remains after that is accounted for.

## Table R9. Template confusion: does the scan pick the true winner?

Admissibility asks whether a signal would be *detected*. A model comparison
claims something stronger - that the *right template* would be preferred. So:
build a sky from `rings_species`'s own fitted background plus a known
template, and run the whole 40-template ordering on it. The true winner is
known by construction.

| true sky | injected flux | truth's rank | best dark matter | best stellar | prefers |
|---|---|---|---|---|---|
| `vvv` (baseline flux) | 6.393e-08 | 1/40 | 227.4 | 255.2 | **stellar** |
| `hestia_G3.2_ang90` (baseline flux) | 6.393e-08 | 1/40 | 81.8 | 65.0 | **dark matter** |
| `gnfw_g1.2` (baseline flux) | 6.393e-08 | 6/40 | 94.2 | 81.6 | **dark matter** |
| `vvv` (own fitted flux) | 2.148e-07 | 1/40 | 1834.7 | 2093.8 | **stellar** |
| `hestia_G3.2_ang90` (own fitted flux) | 2.681e-07 | 1/40 | 1159.4 | 921.3 | **dark matter** |
| `gnfw_g1.2` (own fitted flux) | 1.469e-07 | 2/40 | 462.9 | 359.2 | **dark matter** |

Six for six on the morphology class, and the exact template first in four of
six - in the two misses the templates above the truth are its own immediate
neighbours (gNFW at adjacent gamma, and their top-hat-smoothed twins) within
2.5 in TS, which is the degeneracy sec 6.2's correlation matrix already
reports. **The pipeline is not biased towards the bulge: when the sky
contains dark matter, the fit prefers dark matter.**

One number is worth stating because it runs the other way. On a sky that
truly contains VVV at its fitted flux, the best dark-matter template reaches
1834.7 against VVV's 2093.8, a gap of 259. On the real data the gap is 888.
The data discriminate **more** sharply than a simulated pure bulge does, so
the ordering is not being manufactured by the dark-matter templates being
handicapped - if anything they are doing worse on the sky than they can do
against a bulge. Whether that is extra bulge-specific structure or a residual
the bulge template happens to absorb is not something this test separates.

## Table R10. The 1-3 GeV band alone, where the IC interpolation is safe

Pohl's inverse-Compton maps carry stored planes at 1.32, 2.28, 4.05 and
7.4 GeV. The 1-3 GeV window used here spans 1.0-2.37 GeV and brackets two
stored planes, so it is genuine interpolation with no extrapolation - unlike
the 3-10 and 10-50 GeV bins, where the top plane is exceeded.

`rings_species`, Mask B + PS, 3 bins, mask keeps 79.9%.

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 1142.6 | -0.0 |
| `f98` | 880.3 | +262.4 |
| `hestia_G3.2_ang90` | 697.5 | +445.1 |
| `hestia_G3.2_ang60` | 575.3 | +567.4 |
| `hestia_G3.2_ang30` | 563.2 | +579.4 |
| `hestia_G2.2_ang0` | 554.0 | +588.6 |

VVV still wins, by dAIC +445 over the best of the 24 HESTIA templates, and
all three controls are admissible in the band (VVV injected 103.1, HESTIA
41.5, gNFW 55.1 at a common flux). **The IC interpolation is a footnote, not
a driver**: the ordering in the cleanest band is the same ordering as in the
full 1-10 GeV fit.

## Table R11. Mask A under Pohl22

| | kept | VVV | best HESTIA | gNFW 1.2 | anything admissible? |
|---|---|---|---|---|---|
| Mask A, sources masked | 8.7% | 43.9 | hestia_G3.2_ang30 30.0 | 7.4 | no |
| Mask A, sources also modelled | 8.7% | 32.6 | hestia_G3.2_ang30 22.1 | 6.1 | no |

Mask A masks **every** 4FGL source with a disc of radius r95 + 0.5 deg. With
the real catalogue that is 715 sources, and at 1 GeV the discs overlap enough
to leave **8.7%** of the ROI. (The self-test reports 75.6% because it uses a
synthetic catalogue with a handful of sources - not a bug, but the self-test
does not exercise the real mask.) Nothing is admissible: a signal injected at
a common flux returns TS 6.4-18.8 against the threshold of 25. The baseline
normalisation also runs away - the flux the `gll_iem_v07` fit reports under
Mask A is 1.28e-06, twenty times the Mask B value, which is what an
under-constrained fit on 8.7% of the pixels looks like.

**So Mask A cannot localise the Macias-Di Mauro contradiction.** Dark matter
does not appear under it - the ordering is the same, VVV 43.9 against the
best HESTIA at 30.0 - but every number in it is uninterpretable. And Mask A
is not Di Mauro's configuration in any case: he *models* point sources
through Level-2 iterative source finding rather than masking them, which is
what Mask B + PS approximates. The contradiction therefore localises to the
pipeline rather than the mask, by elimination rather than by a positive test,
and testing it properly needs Level-2 source finding - Gap 5's deliverable.

## Table R12. The six diffuse-freedom splits, Mask B and Mask D

`configs/diffuse_freedom.yaml`, fixed before any run. A split is usable as
evidence only if it passes **both** criteria: it keeps at least half the
injected signal the rigid split keeps (otherwise it is degenerate with the
signal and cannot speak either way), and it still finds a signal added to
the real data at a gain of at least 25 (otherwise a low observed TS from it
means nothing). Recomputed with the converged fitter.

### Mask B + PS, 1-10 GeV, keeps 85.9%, control VVV injected at 6.3931e-08 ph cm^-2 s^-1

| split | extra params | observed TS | injected TS | boosted gain | injected kept | verdict |
|---|---|---|---|---|---|---|
| `rigid` | 0 | 203.1 | 184.3 | +585.0 | 100% | usable |
| `ns` | 1 | 146.7 | 174.8 | +504.7 | 95% | usable |
| `lat_outside` | 4 | 264.4 | 154.2 | +552.1 | 84% | usable |
| `lat_inside` | 4 | 99.0 | 86.8 | +241.4 | 47% | **no** - degenerate with the signal |
| `lat` | 3 | 15.4 | 144.2 | +238.9 | 78% | usable |
| `lat_ns` | 7 | 8.7 | 137.2 | +203.8 | 74% | usable |

Usable: `rigid`, `ns`, `lat_outside`, `lat`, `lat_ns`. Observed TS across them 8.7 to 264.4,
a factor of 30.

### Mask D + PS, 1-10 GeV, keeps 71.0%, control VVV injected at 1.8799e-08 ph cm^-2 s^-1

| split | extra params | observed TS | injected TS | boosted gain | injected kept | verdict |
|---|---|---|---|---|---|---|
| `rigid` | 0 | 22.9 | 37.8 | +66.3 | 100% | usable |
| `ns` | 1 | 4.6 | 33.4 | +38.5 | 88% | usable |
| `lat_outside` | 4 | 79.9 | 36.5 | +91.7 | 97% | usable |
| `lat_inside` | 4 | 8.9 | 10.7 | +26.7 | 28% | **no** - degenerate with the signal |
| `lat` | 3 | 0.0 | 27.9 | +12.4 | 74% | **no** - cannot find an added signal |
| `lat_ns` | 7 | 0.0 | 25.6 | +8.7 | 68% | **no** - cannot find an added signal |

Usable: `rigid`, `ns`, `lat_outside`. Observed TS across them 4.6 to 79.9,
a factor of 17.

**The three that fail under Mask D, named, with their boosted gains:**
`lat_inside` (gain +26.7) fails as **degenerate** - it keeps only 28% of the
injected signal, so it absorbs what it is meant to measure. `lat` (gain
+12.4) and `lat_ns` (gain +8.7) are **not** degenerate - they keep 74% and
68% - they fail because they cannot find a signal added to the real data.
Those are different failures and the code used to print one label for both;
fixed, and the reason is now in the JSON.

Under Mask B only `lat_inside` fails, and marginally: it keeps 47% against a
threshold of 50%. Note that under the *absolute* criterion the ring work
uses - injected TS >= 25 rather than a fraction of the rigid split's - Mask
B's `lat_inside` would pass at injected TS 86.8. The two criteria disagree on
exactly that one cell, and the conventions should name which one governs.

**So, for the conventions:** under Mask B + PS five of six splits are
admissible; under Mask D + PS three of six are, and the three that fail are
the ones that give the diffuse model latitude freedom at inner longitudes
(`lat_inside`) or everywhere (`lat`, `lat_ns`).

### What the fitter fix did to these numbers

The verdicts are unchanged - the same splits are usable under each mask -
but the magnitudes moved, and the headline factor with them.

| | usable splits | observed TS range | factor |
|---|---|---|---|
| Mask B + PS, before | 5/6 | 3.8 - 305.6 | **80** |
| Mask B + PS, after | 5/6 | 8.7 - 264.4 | **30** |
| Mask D + PS, before | 3/6 | 1.7 - 79.3 | 46 |
| Mask D + PS, after | 3/6 | 4.6 - 79.9 | 17 |

Both extremes moved inward, which is what fixing an under-converged fit
does: an unconverged fit understates lnL, and it understated it most where
the model had the most freedom. **Section 6.0's 'factor of eighty' is a
factor of thirty.** The conclusion it supports is unchanged - the background
model still dominates the measurement, and the ring decomposition still
narrows it further, to 11 - but the number itself has been quoted in the
memo, the README and the mini-paper and needs correcting in all three.

One loose end, stated rather than hidden: under Mask D the *injected* TS
also moved a long way (rigid 8.6 to 37.8). Part of that is a small Mask D
change between the runs (kept fraction 70.75% to 71.04%) and a 15% shift in
the control flux, which together account for perhaps a third of it. I have
not isolated the rest. The verdicts do not depend on it and the new numbers
come from the current, converged code, but it is worth one hour before the
memo is final.

## Table R13. Level 2: iterative source finding instead of a fixed mask

The last candidate explanation for Macias versus Di Mauro. He lets the fit
add point sources; we mask a fixed catalogue. Mask A eliminated itself
(Table R11), so this is what was left.

One design choice decides the answer, so both ways were run. **neutral**:
find sources with no signal in the model, then fit every template against
that one common source list - equal parameter count, but the iteration has
already had free rein to absorb the excess, so it is biased *against* any
signal. **joint**: find sources with the candidate already in the model,
one run per hypothesis, which is what Di Mauro does - each hypothesis gets
its best shot, at the cost of different source counts, so the comparison
uses the full AIC = 2k - 2lnL rather than dAIC = -2 dlnL.

### neutral arm: 12 sources added, one common list

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 1018.2 | -0.0 |
| `f98` | 840.3 | +177.9 |
| `hestia_G3.2_ang90` | 567.8 | +450.3 |
| `hestia_G3.1_ang60` | 515.1 | +503.0 |
| `hestia_G2.2_ang90` | 473.8 | +544.4 |
| `hestia_G2.2_ang0` | 466.6 | +551.6 |

Sources found, in order (residual TS at the moment each was added):

| source | l | b | residual TS | pass |
|---|---|---|---|---|
| src000 | -18.35 | +0.15 | 681 | 100 |
| src001 | -3.65 | +1.15 | 324 | 100 |
| src002 | -10.25 | +0.15 | 283 | 100 |
| src003 | +17.35 | +2.45 | 203 | 100 |
| src004 | +6.35 | +4.95 | 168 | 100 |
| src005 | +18.85 | -0.35 | 138 | 100 |
| src006 | -0.35 | +0.35 | 93 | 64 |
| src007 | -3.05 | -2.25 | 84 | 64 |
| src008 | +3.55 | +1.95 | 66 | 64 |
| src009 | +14.95 | +9.25 | 63 | 36 |
| src010 | -17.15 | -2.75 | 50 | 36 |
| src011 | +16.85 | -1.85 | 44 | 36 |

One of them, `src006`, sits at l = -0.35, b = +0.35 - on the Galactic
Centre. That is the mechanism this test exists to catch: a point source
placed on the excess eats it. It does, and the ordering still holds.

### joint arm: one source search per hypothesis

| hypothesis | sources added | free normalisations | lnL | AIC | dAIC |
|---|---|---|---|---|---|
| `vvv` | 10 | 152 | 554,484.95 | -1,108,665.89 | +0.0 |
| `f98` | 10 | 152 | 554,354.82 | -1,108,405.63 | +260.3 |
| `gnfw_g1.4` | 13 | 176 | 554,288.99 | -1,108,225.98 | +439.9 |
| `gnfw_g1.2` | 13 | 176 | 554,235.08 | -1,108,118.16 | +547.7 |
| `hestia_G3.2_ang90` | 10 | 152 | 554,182.35 | -1,108,060.69 | +605.2 |
| `hestia_G2.2_ang90` | 10 | 152 | 554,092.06 | -1,107,880.11 | +785.8 |

**Both arms put the stellar bulge first.** The neutral arm gives VVV 1018.2
against 567.8 for the best of the 24 HESTIA templates, dAIC +450, and that
is the arm biased against a signal. The joint arm, where every hypothesis
runs its own source search, puts VVV ahead of the best dark-matter template
by +605 on full AIC - and the gNFW hypotheses needed **three extra point
sources** (13 against 10) to patch what their template could not fit, and
lost anyway.

**So Level 2 is not the explanation either.** Within this implementation,
iterative source finding does not reverse the ordering; it lowers every TS
by roughly half (VVV 2151.7 to 1018.2) and leaves the ranking alone.

The honest limit: this is a reimplementation of the procedure Di Mauro
describes, not his pipeline. Ours adds 10-13 sources over a 40x40 deg ROI
in 1-10 GeV; his works over a larger region and a wider band and adds far
more. A null from our version is evidence that source finding does not
flip the ordering *here*, not proof that his does not flip it *there*.

## Table R14. What Level 2 does to a morphology, and what it invents

Three simulated skies, each built from `rings_species`'s own fitted
background, run through the full iteration and then the 40-template
ordering.

| simulated sky | injected flux | sources added | winner | truth's rank |
|---|---|---|---|---|
| flattened HESTIA G3.2 90 deg | 2.1477e-07 | **1** | `hestia_G3.2_ang90` at TS 465.9 | 1/40 |
| spherical gNFW 1.2 | 2.1477e-07 | **1** | `gnfw_g1.2` at TS 378.3 | 1/40 |
| nothing at all | 0.0000e+00 | **0** | `hestia_G3.2_ang90` at TS 5.3 | n/a |

**The iteration does not round off morphology.** A flattened halo injected
at 2.15e-07 comes back as itself, ranked 1 of 40 - and the spherical
control at matched flux also comes back as itself, which is the half that
makes the first half mean anything. Without it, 'the flattened one survived'
would be equally consistent with a pipeline that returns whatever it is
given.

**And it invents nothing.** On a sky with no signal and no sources beyond
the background model, the iteration adds **zero** sources and the best of
all 40 templates reaches TS 5.3. The 12 sources it finds on the real sky
are therefore structure the ring model misses, not an artefact of the
method - which is also the first real number for Gap 5's phantom question.

One observation worth carrying: the iteration adds 12 sources to the real
sky and 0-1 to a simulated one drawn from the same background. The
simulated skies are drawn from the model, so by construction they contain
nothing the model cannot describe. The gap between 12 and 0 is a measure of
how much real structure Pohl22 still misses.

## Table R15. The Mask D injected-TS discrepancy, closed

Flagged in Table R12: recomputing the diffuse-freedom splits with the
converged fitter moved Mask D's injected TS from 8.6 to 37.8 for the rigid
split, which a small mask change and a 15% control-flux shift could not
explain. Afeefa asked for it closed or documented. It is closed.

**It is the fitter bug, and Mask D is the configuration most exposed to**
**it.** The stall needed a component sitting at exactly zero whose Newton
direction points negative. Under Mask D the signal normalisation is driven
to exactly zero in **three of eight bins**; under Mask B, in one.

| bin | E GeV | Mask B + PS, pinned at zero | Mask D + PS, pinned at zero |
|---|---|---|---|
| 0 | 1.15 | bubbles | bubbles, **sig** |
| 1 | 1.54 | bubbles | bubbles |
| 2 | 2.05 | bubbles | - |
| 3 | 2.74 | - | - |
| 4 | 3.65 | - | - |
| 5 | 4.87 | - | - |
| 6 | 6.49 | - | **sig** |
| 7 | 8.66 | **sig** | **sig** |

The simulated sky for each split is built from that split's own fitted
normalisations. A stalled fit therefore produced a wrong sky, and the
injected TS measured on it was wrong - suppressed, and suppressed across
every split, which is what the old Mask D table shows (injected TS 7.6 to
9.4 for all six splits, against 25.6 to 37.8 now). The 15% control-flux
change has the same cause: the flux is the sum of the fitted signal
normalisations, and three of them were being solved at a point the
iteration never reached.

Mask B's injected TS moved by 1.4% over the same fix (181.7 to 184.3),
which is the control on this explanation: one pinned bin, almost no change;
three pinned bins, a factor of four.

A hypothesis that was tested and is wrong, recorded so nobody retests it:
the old run was **not** secretly using what is now Mask E. Mask E + PS
keeps 86.7% of the ROI against the old run's 70.75%, so the relabel of
2026-09-27 is not involved.

## Table R16. What Level 2 found, cross-matched against 4FGL

`analysis/level2/source_record.py`. Every mechanical column is filled in;
`verdict` and `note` are left for a human. PSF evaluated at 1.15 GeV, the
lowest bin, where it is widest - so 'inside the PSF' is the generous
reading and anything called uncatalogued here is uncatalogued at every
energy.

| source | l | b | residual TS | nearest 4FGL | sep deg | inside 95% PSF |
|---|---|---|---|---|---|---|
| `src000` | -18.35 | 0.15 | 681.1 | 4FGL J1651.7-4359  | 0.145 | yes |
| `src001` | -3.65 | 1.15 | 324.0 | 4FGL J1732.5-3131  | 0.148 | yes |
| `src002` | -10.25 | 0.15 | 282.7 | 4FGL J1718.0-3726  | 0.03 | yes |
| `src003` | 17.35 | 2.45 | 202.8 | 4FGL J1813.4-1246  | 0.106 | yes |
| `src004` | 6.35 | 4.95 | 167.8 | 4FGL J1741.9-2054  | 0.088 | yes |
| `src005` | 18.85 | -0.35 | 137.7 | 4FGL J1826.1-1256  | 0.295 | yes |
| `src006` **(Galactic Centre)** | -0.35 | 0.35 | 92.7 | 4FGL J1744.9-2905  | 0.34 | yes |
| `src007` | -3.05 | -2.25 | 83.9 | 4FGL J1746.8-3239  | 0.082 | yes |
| `src008` | 3.55 | 1.95 | 66.1 | 4FGL J1745.7-2506  | 0.259 | yes |
| `src009` | 14.95 | 9.25 | 62.6 | 4FGL J1744.4-1135  | 0.176 | yes |
| `src010` | -17.15 | -2.75 | 50.0 | 4FGL J1709.7-4429  | 0.256 | yes |
| `src011` | 16.85 | -1.85 | 44.3 | 4FGL J1829.4-1500c | 0.326 | yes |

**All twelve fall inside the 95% PSF containment radius of a catalogued
4FGL source, and ten of the twelve within 0.3 degrees.** Level 2 is
rediscovering sources that are already known, not inventing new ones. With
the null-sky result in Table R14 - zero sources added to a sky containing
none - that is two independent statements that the iteration does not
produce phantoms, and the first hard number for Gap 5's phantom question.

Why it adds them at all, given that Mask B + PS already *models* the
catalogue: the fitted point-source normalisation comes out at 0.76 to 0.92
rather than 1.0 in every bin (sec 6.2). The template under-predicts those
sources, and the iteration tops them up at the same positions. That is a
statement about our point-source spectral integration, not about the sky.

`src006` at l = -0.35, b = +0.35 matches 4FGL J1744.9-2905 at 0.34 deg -
a catalogued Galactic Centre source. So the source that sits on the excess
and eats part of it is a real one, which makes the ordering surviving it a
stronger result than if it had been a phantom.

The judgement left for a reviewer is therefore thinner than expected: there
are no uncatalogued candidates to adjudicate. What is still worth a human
is whether any of the twelve is a *different* source from its match rather
than the same one re-found, and `src011`'s match being a `c`-flagged
source (4FGL J1829.4-1500c, flagged in the catalogue as possibly
confused with interstellar emission) is the one place that matters.
