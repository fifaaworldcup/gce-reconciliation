# Our HESTIA pipeline against Muru's, convention by convention

What this answers: whether our axis ratios can be compared to the published
HESTIA numbers, and whether the J-factor result we want to publish survives
being redone Muru's way. Source read: `gitlab.aip.de/muru/gce_in_hestia` at
`b50114d9` — `src/GalacticCoordinates.jl`, `src/BulgeDensityProjection.jl`,
`src/calc_methods_densityprojections.jl`, `notebooks/azimuthal_analysis_erebos.jl`.
Code: `analysis/gap1/muru_match.py` and the three drivers beside it. Results:
`outputs/gap1/muru_match.json`.

Short version. Of the nine conventions that matter, seven are fixed by the
public code and we now match all seven. The two that are not fixed are the two
that set the axis ratio — the disc normal and the in-plane azimuth reference —
and both come from `get_galvectors_from_profile` in `HestiaUtils.jl`, which the
repository calls but does not contain. Matching the other seven does not close
the gap: our reconstruction comes out 0.08–0.14 rounder than his maps and
*anti*-correlated with them across the six galaxies. The J-factor result is
unaffected, because it compares two maps built in whatever frame we chose, and
the comparison is the same under two normals 47° apart.

## 1. The nine conventions

| | Muru, from the source | ours before | ours now |
|---|---|---|---|
| observer distance | `8.0*0.677` kpc/h = 8.0 kpc | 8.2 kpc | 8.0 kpc |
| line-of-sight cut | `d <= 15*0.677` kpc/h, 15 kpc **from the observer** | ±10 kpc box about the centre | 15 kpc from the observer |
| sight-line kernel | top-hat of fixed **angular** radius, KDTree on (l,b) in degrees | Gaussian on a 3D CIC grid | top-hat, 1.0° and 3.0° |
| grid | `-20:0.1:20` in l and b | ±20°, 0.1° | same |
| particles | `ptype == :dm`, everything in the halo file | DM within 10 kpc | all DM in the field |
| squaring | `dens .^ 2` applied to the finished column map | same | same |
| axis ratio | min/max singular value of the mean-centred (l,b) of every pixel ≥ frac × peak, unweighted | eigenvalues of the same covariance — algebraically identical | identical |
| isophote | frac ∈ 0.7, 0.5, 0.3 of peak | 0.05 of peak | both, and a scan |
| disc normal, azimuth zero | `get_galvectors_from_profile`, **not in the repository** | stellar inertia tensor, stars < 10 kpc | AHF profile row at r ≈ 10 kpc/h, and the stellar one for comparison |

Three of those are worth spelling out.

**His map is not a column density.** The cone has a fixed *angular* radius, so
its solid angle is the same at every distance and its cross-section grows as
d². Summing the masses inside it gives ∫ρ s² ds per unit solid angle, not
∫ρ ds. The far half of the 15 kpc sight line is weighted an order of magnitude
above the near half. Nothing in his paper depends on the absolute value, and
our fits float a normalisation per bin, so this changes no published number —
but it means "projected column density" in the paper is a d²-weighted one, and
anyone deriving a flux from these maps needs to know.

**The axis ratio is orientation-blind, and that settles the reciprocal
question.** `minimum(Σ)/maximum(Σ)` is the minor over the major axis of the
super-level region whichever way round they lie, so it is ≤ 1 by construction.
In all 24 of the maps he sent, the major axis is in longitude, so his number is
the vertical over the in-plane extent — **our q, not its reciprocal.** Di
Mauro's r still needs inverting; Muru's does not.

**Squaring a column map cannot change its shape.** `{m² ≥ f·max(m²)}` is
`{m ≥ √f·max(m)}`, the same set of pixels. So under his own estimator
q(m², f) = q(m, √f) exactly — squaring does not produce a new morphology, it
relabels the level sequence. `muru_match.selftest()` checks this to the last
bit at three levels. This is the Abazajian, Kumar & Macias point proved inside
the code it is about, and it is why the only informative comparison is against
a genuine ∫ρ² ds, computed before projection.

## 2. The J-factor result, redone Muru's way

J-factor from the particles directly: ρ at each particle from its 32nd nearest
neighbour, floored at the 220 pc softening, then Σ m_i ρ_i in the same cone
that gives Σ m_i for the column. Both weightings therefore carry the same d²,
and the comparison isolates the squaring alone.

At the 50% isophote with his 3° cone, over all six galaxies and four viewing
angles, the J-factor map is **flatter than the squared column by
q(J) − q(col²) = −0.032 ± 0.020**, the same sign in 22 of 24 cases. Under the
stellar disc normal instead — a frame 17° away for the G1 pair and 47–48° away
for the other four — it is **−0.041 ± 0.016, the same sign in 24 of 24.**

That is a factor of four tighter than the ±0.088 we reported from the earlier
reconstruction, and the sign is now consistent, where before it was not. The
difference is real and it is small: squaring the column over-rounds the halo by
about 0.03–0.04 in q, which matters for a measurement quoting q to ±0.04 and
not otherwise.

## 3. Why the per-galaxy comparison to his numbers cannot be made

`get_galvectors_from_profile(id, sim)` returns the disc normal and an in-plane
reference vector. It is in `HestiaUtils.jl`, which the README calls "a
non-listed package". The AHF profile files we have offer the pieces it must be
built from, and they are not unique:

- Its shape-tensor minor axis `Ec` at the r = 10.28 kpc/h row (15.2 kpc) sits
  **0.1°–3.1° from the angular momentum direction** at the same radius, so
  those two candidates are the same choice.
- That axis sits **17° from the stellar disc normal for G1.1 and G1.2, and
  47–48° away for G2.1, G2.2, G3.1 and G3.2.** AHF's profile is over all
  particles and is reporting the dark halo, which in four of these six is not
  aligned with the stars.
- Switching between the two moves the angle-marginalised q by up to 0.19 for a
  single galaxy.

The in-plane reference is worse. The azimuth zero point has to be a direction
in the disc plane, and the AHF profile's in-plane axes are degenerate: b/a is
**0.977 to 0.994** at that row, so the major axis within the plane is
numerically undefined. Our "angle 0" is at an unknown offset from his.

Both problems were tested rather than assumed.

- **Same-angle match.** Shape-only correlation (log map minus its own
  azimuthal average) between our reconstruction and his template at the same
  quoted angle: the right halo wins 7 of 24 at the 1° cone. Not a match.
- **Azimuth free.** Scanning our azimuth over 0–355° does not fix it. Our
  G1.1 is the best match to *twenty-two of his twenty-four maps*, at r = 0.84
  to 0.93, and the best-matching azimuth does not step with his quoted angle.
- **The control that explains it.** His own 24 maps correlate with each other
  at mean r = 0.81 in the shape residual and 0.966 raw — and the
  same-galaxy-different-angle correlation (0.81–0.90) is no higher than the
  different-galaxy-same-angle correlation (0.80–0.82). His maps cannot be told
  apart from each other by shape, so no map-level comparison can establish
  identity in either direction, ours or his.
- **Angle-marginalised q**, which needs no zero point: ours is +0.135 rounder
  than his at the 3° cone (+0.084 at 1°, +0.079 with the stellar normal) and
  the per-galaxy ordering is anti-correlated, Spearman −0.66 to −0.83.

So the comparison is not available from the public material, and the reason is
specific: one missing orientation function. That is the only thing worth
asking for if we ever ask — not a discussion, a ten-line function, or the AHF
halo catalogue columns it reads.

## 4. What this does to §6.1

The §6.1 correlation is ρ(q, ΔAIC vs VVV) over 24 halo × angle fits. Its q
comes from Table 1: our estimator, at 5% of peak, on Muru's **uncropped** ±40°
maps, averaged over the four viewing angles. Table 1 reproduces from the
cropped templates to ±0.007 for five of the six galaxies, so the number itself
is sound. The correlation built on it is not stable under any of the three
choices nobody fixed in advance.

| what is varied | ρ | p | n |
|---|---|---|---|
| as published (Table 1 q, repeated per angle) | **+0.656** | 5×10⁻⁴ | 24 |
| the same, six galaxies only | +0.771 | 0.072 | 6 |
| q measured on the **cropped** ±20° map instead | +0.268 | 0.21 | 24 |
| q varying **per angle** rather than per galaxy | +0.217 | 0.31 | 24 |
| isophote at 10% of peak instead of 5% | +0.463 | 0.023 | 24 |
| isophote at 50% of peak, where Muru measures | −0.251 | 0.24 | 24 |
| isophote at 2% of peak | −0.653 | 5×10⁻⁴ | 24 |

Eleven of the 24 level sets at 5% of peak touch the ±20° crop edge, which is
why the cropped and uncropped numbers differ; at 2% all 24 do and the shape is
the ROI's, not the halo's (Table 2). Within a single galaxy, across its four
viewing angles, ρ is +0.4, +0.4, +0.4, −0.4, −0.2 and 0.0 — a coin flip. The
whole correlation is the between-galaxy spread of six numbers, and at n = 6 it
is p = 0.07.

Afeefa's four-way run, all at 5% of peak against the same 24 fits:

| | ρ | p |
|---|---|---|
| squared column, current reconstruction | +0.428 | 0.037 |
| J-factor, current reconstruction | +0.047 | 0.83 |
| squared column, matched reconstruction | +0.223 | 0.30 |
| J-factor, matched reconstruction | −0.062 | 0.77 |

None of them reproduces +0.73. §6.1 was already withdrawn in §6.0 because the
background model moves the detection by a factor of eighty. This is a second
and independent reason: the shape statistic it correlates is not robust to the
isophote, the crop, or whether q is allowed to vary with viewing angle.

## 5. Labels

The particle file names are wrong and the map file names cannot be checked.

The six Arrow files are named `..._G11_...` through `..._G32_...`, and within
each pair the two are swapped: the file named `G11` holds the halo at G1.2's
AHF centre, and G31/G32 likewise, 850+ kpc apart. `reduce.py` ignores the
names and identifies each file by matching its dark-matter centroid to the AHF
centres in the data readme, which is an independent handle. Fixed there, and
`whichcentre.py` prints the match for anyone who wants to re-check it.

The map files carry the simulation ID and halo ID inside their names, and
`hestia_to_template.GALAXIES` maps those to G-labels in the same grouping
Muru's own `scripts/plot_6_galaxy_projections.jl` uses. So §6.1's per-galaxy
split — G2.1 and G1.1 beating VVV, G3.1 and G3.2 losing — does not come from
the renamed particle files and the swap does not touch it.

But it cannot be verified either. The check would be to rebuild each map from
the correctly identified particles and see which published map it matches, and
§3 above shows that test has no power: his maps are mutually indistinguishable
by shape. The §6.1 ordering rests on his filenames being right, and we have no
independent handle on them, which is worth stating given that the one place a
handle existed is the one place the labels turned out to be wrong.

## 6. Reproduce

```bash
python3 analysis/gap1/muru_match.py                 # self-tests, no data needed
python3 analysis/gap1/muru_match_run.py \
  --reduced <dir> --profiles <dir> \
  --current outputs/gap1/hestia_jfactor.json
python3 analysis/gap1/muru_match_check.py <reduced> <profiles>   # label test
python3 analysis/gap1/muru_match_azimuth.py <reduced> <profiles> # azimuth free
python3 analysis/gap1/muru_match_noise.py <reduced> <profiles>   # shot noise
```

`<reduced>` is the output of `reduce.py` on the six Arrow files: DM positions
in physical kpc about the AHF centre and masses in 10¹⁰ M☉. `<profiles>` holds
the six `profile_G**_AHF.txt`.
