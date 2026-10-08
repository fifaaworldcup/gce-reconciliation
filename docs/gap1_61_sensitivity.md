# Section 6.1, rewritten as a sensitivity study

Section 6.1 reported that the projected axis ratio of a HESTIA template
predicts how well it fits the GCE relative to the VVV bulge, Spearman
rho = +0.73 at p = 5e-5 over 24 halo x angle combinations, with the crossover
at q ~ 0.6. It was withdrawn in sec 6.0 because the background model moves the
detection by a factor of eighty. It stays withdrawn, and this is the
replacement: not a result, a measurement of how much of that correlation was
the sky and how much was our analysis choices.

The honest summary is that rho is not a number this dataset has. It ranges
from -0.65 to +0.46 across choices that were never fixed in advance, and the
published +0.73 sits at one end of that range.

## 1. What varies, holding the data and the fits fixed

Same 24 templates, same 24 fits, same maps. Only the recipe for q changes.

| q measured as | rho | p |
|---|---|---|
| ellipse fit at 5% of peak, uncropped map, averaged over viewing angle — **as published** | **+0.656** | 5e-4 |
| the same, 6 galaxies rather than 24 rows | +0.771 | 0.072 |
| the same estimator and level, on the cropped +-20 deg ROI | +0.268 | 0.21 |
| per-angle rather than per-galaxy, 5% of peak | +0.217 | 0.31 |
| super-level SVD at 10% of peak | +0.463 | 0.023 |
| super-level SVD at 20% of peak | -0.130 | 0.54 |
| super-level SVD at 50% of peak — the level Muru uses | -0.251 | 0.24 |
| super-level SVD at 2% of peak | -0.653 | 5e-4 |

Three separate choices move it, and none was pre-registered.

**The isophote.** q is not a property of a halo, it is a property of a halo at
a radius. These halos get rounder outward until the level set reaches the ROI
boundary and then abruptly rounder still. Measured in the core (50% of peak,
about 1.6-4.2 deg) the correlation is negative; measured in the outskirts (10%,
about 9 deg) it is positive and nominally significant; measured at 2%, where
all 24 level sets have hit the edge of the map, it is strongly negative and is
an artefact of the crop, not of the sky (this is the same effect as the fake
boxiness in Table 2).

**The crop.** At 5% of peak, 11 of the 24 level sets touch the +-20 deg ROI
edge. Table 1 was measured on Muru's uncropped +-40 deg maps and reproduces
from the cropped templates to +-0.007 for five of six galaxies — but for G2.1
the crop moves q from 0.590 to 0.677, which is enough to reorder the sample.
So the published correlation compares the fit quality of a template inside the
ROI against the shape of that template half of which the fit never saw.
Defensible, but a choice.

**Per-angle versus per-galaxy.** Table 1 averages q over the four viewing
angles, so the 24 rows carry only six distinct values of q. The scatter within
a galaxy is small (+-0.014 to +-0.024) against the spread between galaxies
(0.55 to 0.69), so the correlation is entirely between-galaxy and the n = 24
p-value of 5e-5 is not a valid test of anything: the effective sample is six,
where p = 0.07. Letting q take its actual per-row value drops rho to +0.217.
Within a single galaxy, across its four angles, rho is +0.4, +0.4, +0.4, -0.4,
-0.2 and 0.0.

## 2. What varies when the quantity changes

Afeefa's four-way run, all at 5% of peak against the same 24 fits.

| | rho | p |
|---|---|---|
| squared column, our particle reconstruction | +0.428 | 0.037 |
| annihilation integral, our particle reconstruction | +0.047 | 0.83 |
| squared column, reconstruction matched to Muru's conventions | +0.223 | 0.30 |
| annihilation integral, matched reconstruction | -0.062 | 0.77 |

All four differ, so by her rubric both the quantity and the reconstruction are
in play. But neither is the dominant term: the spread across isophotes on a
single map, +0.46 to -0.65, is wider than the spread across quantities and
reconstructions at a fixed isophote, +0.43 to -0.06.

## 3. What happens on a background built from gas astronomy

Everything above holds the background fixed at a single all-sky interstellar
template, which sec 6.0 showed is the dominant systematic. Replacing it with
Pohl22's sixteen-component ring decomposition removes the correlation
altogether.

| q used, Pohl22 rings_free, Mask B + PS, 1-10 GeV | rho | p |
|---|---|---|
| Table 1 q — the published recipe | **-0.056** | 0.79 |
| per-angle q, 5% of peak | +0.400 | 0.053 |
| per-angle q, 10% of peak | +0.028 | 0.90 |
| per-angle q, 50% of peak | +0.383 | 0.064 |

The published recipe gives nothing. Two other recipes give a marginal positive
at p ~ 0.05-0.06, which is what one expects from four looks at a null.

The reason is visible in the fit itself rather than the correlation. On the
ring background VVV reaches TS 794 and the best of the 24 HESTIA templates
reaches 260, a gap of dAIC +534; the second stellar bulge, F98, is +133. On
the single-template background the same data had HESTIA G2.1 ahead of VVV by
dAIC -58. So the ordering the correlation was built to explain is reversed by
the background model, and the quantity the correlation claimed the data
constrain — q, not "dark matter versus stars" — is no longer the quantity
that separates the templates.

## 4. What to carry forward

- The axis-ratio correlation is withdrawn and is not being re-reported. Any
  paper that quotes a q for the GCE from template fitting should state the
  isophote, whether the shape was measured inside or outside the fitted
  region, and whether the degrees of freedom are templates or halos. Ours did
  not, and the number moved by 1.1 in rho across those choices.
- The useful residue is the method, not the measurement: hold the data and the
  fits fixed and vary only the shape statistic. On this dataset that single
  test would have caught the result before it was written down.
- The incompleteness calibration of sec 6.6 is untouched by all of this. It is
  a differential measurement on simulated skies where the background is exact
  by construction, and it still says that past a certain unresolved
  point-source flux any injected morphology is driven to q ~ 0.59. That is a
  reason to distrust a q ~ 0.6 measurement independently of everything above.
