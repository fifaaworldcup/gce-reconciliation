# The diffuse-freedom test

**Oren · 2026-09-29 · `analysis/diffuse_freedom/` · `configs/diffuse_freedom.yaml`**

> Macias and Di Mauro use the same bulge templates and the same gas maps and reach opposite conclusions: NFW² at 17.6σ dropping to 2.4σ once the X-bulge and nuclear bulge are included, against TS 8×10²–1.5×10⁴ for dark matter with NB and BB included. The difference between them is residuals, 30–40% against ≤10%.
>
> This test reproduces that disagreement inside a single pipeline, on one dataset, with one set of templates. **A defensible change in how much the interstellar model is allowed to bend moves a detection by a factor of eighty.** No difference in data, templates, mask or event selection is required.

---

## 1. The question

A likelihood fit reports a TS for a signal template against a background. If the background model is too rigid to describe the sky, the fit will use whatever component it has to make up the difference — and a signal template that happens to resemble the error is the cheapest correction available.

Gap 1's background is one interstellar template, `gll_iem_v07`, with one free normalisation per energy bin. Four numbers across 1–10 GeV. Published analyses decompose the diffuse emission into Galactocentric rings and gas phases; Di Mauro (2026) fits about twenty ring templates per IEM. The question is what that difference is worth in TS.

## 2. The method

Split the interstellar plane into components that each carry their own free normalisation, so the model has freedom in shape as well as scale, and refit. Six splits, fixed in `configs/diffuse_freedom.yaml` before any run and applied identically to every template:

| split | what it allows | extra free parameters |
|---|---|---|
| `rigid` | one normalisation per energy bin | 0 |
| `ns` | north and south independently | 1 |
| `lat_outside` | four \|b\| bands, but only at \|l\| > 10° | 4 |
| `lat_inside` | four \|b\| bands, but only at \|l\| ≤ 10° | 4 |
| `lat` | four \|b\| bands across the ROI | 3 |
| `lat_ns` | those bands, north and south separately | 7 |

Bands are \|b\| = 0–2°, 2–5°, 5–10°, 10–20°. The two gated splits exist to localise the effect: they ask whether the correction is being applied where the signal lives or away from it.

### The controls are the test

A split that removes a signal proves nothing on its own — a split flexible enough to absorb a real signal would do the same. So every split runs against three skies:

**observed** — the data.

**injected** — a Poisson realisation of the fitted background plus a signal of known strength. Built from the `rigid` fit's own normalisations, with the control template at the flux that fit attributed to it. This sky contains no mismodelling by construction, so it asks one thing: can the split see a signal at all?

**boosted** — the observed data with one more signal of that same known strength added on top. Every real residual is retained, so it asks the harder question: can the split see a signal *through* the mismodelling that is actually there?

A split is admissible as evidence only if `injected` survives it. If `injected` and `boosted` both survive and `observed` does not, the observed excess was the background model's error.

## 3. Gap 1's result

Mask B + PS, 1–10 GeV, 4,933,089 photons, mask keeps 86.0%. IEM family 1 (`gll_iem_v07`). Control template VVV (Coleman et al. 2019), injected at **7.371×10⁻⁸ ph cm⁻² s⁻¹**, the flux the rigid fit attributed to it. Seed 42.

### 3.1 All three skies, all three templates

TS of each signal template. `boosted` is the absolute TS, not the gain.

| split | sky | gNFW γ=1.2 | HESTIA G1.1 | VVV |
|---|---|---|---|---|
| `rigid` | observed | 113.5 | 270.7 | **226.4** |
| | injected | 129.6 | 132.2 | 181.7 |
| | boosted | 559.2 | 880.7 | 904.5 |
| `ns` | observed | 82.9 | 181.7 | 159.4 |
| | injected | 124.3 | 125.9 | 175.5 |
| | boosted | 467.2 | 685.9 | 740.9 |
| `lat_outside` | observed | 163.5 | 206.6 | **305.6** |
| | injected | 104.5 | 109.0 | 148.1 |
| | boosted | 586.8 | 682.7 | 944.1 |
| `lat_inside` | observed | 19.4 | 111.1 | 113.7 |
| | injected | 52.4 | 48.7 | **64.8** |
| | boosted | 147.4 | 336.8 | 397.1 |
| `lat` | observed | 2.0 | 4.5 | 8.2 |
| | injected | 105.0 | 106.0 | 145.6 |
| | boosted | 132.4 | 152.5 | 253.4 |
| `lat_ns` | observed | 0.3 | 0.0 | **3.8** |
| | injected | 101.1 | 104.4 | 143.7 |
| | boosted | 101.7 | 89.1 | 195.6 |

### 3.2 Which splits are admissible

VVV, the control template. "Injected kept" is relative to the rigid split. A split is admissible if it keeps at least half the injected signal and still gains at least 25 in TS on the boosted sky — both thresholds are in the config and were fixed before the run.

| split | extra params | observed TS | injected kept | boosted gain | admissible |
|---|---|---|---|---|---|
| `rigid` | 0 | 226.4 | 100% | +678.1 | yes |
| `ns` | 1 | 159.4 | 97% | +581.5 | yes |
| `lat_outside` | 4 | **305.6** | 81% | +638.5 | yes |
| `lat_inside` | 4 | 113.7 | **36%** | +283.4 | **no — degenerate** |
| `lat` | 3 | 8.2 | 80% | +245.2 | yes |
| `lat_ns` | 7 | **3.8** | 79% | +191.8 | yes |

**Across the five admissible splits the observed TS runs from 3.8 to 305.6, a factor of 80.** Figure: `outputs/gap1/fig8_diffuse_freedom.png`.

### 3.3 What the controls establish

- The splits are not too flexible. Four of the five admissible ones keep 79–100% of a signal injected at the observed flux.
- They are not blind to a signal buried in real residuals either. Adding one more real-strength signal to the actual data is detected at +192 to +678 in TS under every admissible split.
- `lat_inside` is genuinely degenerate — it keeps only 36% of the injected signal — and is therefore reported but excluded. That is the split doing the most obvious thing wrong, and the criterion catches it without anyone deciding after the fact.
- One single extra parameter, the north/south ratio, moves the observed TS by 30% while moving the injected TS by 3%. The asymmetry is the whole finding in miniature.

### 3.4 Where the ambiguity lives

Letting the model bend **away** from the Galactic Centre raises the detection, 226 → 306: the rigid template was mis-fitting the outer plane and dragging the inner fit with it. Letting it bend **across the whole ROI**, through the longitudes the excess occupies, destroys it, 226 → 8.2.

The difference between those two is a single question — may the inner-longitude latitude profile of the diffuse emission differ from what `gll_iem_v07` says? — and **a single all-sky template cannot answer it.** Both answers are defensible before looking at the data, and they disagree by eighty.

### 3.5 The axis ratio says the same thing

Fitting a ladder of gNFW templates that differ only in vertical axis ratio to the real counts (`analysis/gap1/measure_q.py`) drives q to the bottom rung, q ≤ 0.30, under both Mask B + PS and Mask D + PS. Di Mauro reports r = 1.10 ± 0.05, which in our convention is q = 0.909 ± 0.041. A likelihood that wants the flattest template available whatever the floor is buying a correction to the diffuse model, not measuring a shape.

### 3.6 Other masks

| mask | rigid observed TS | `lat_ns` | admissible splits | spread across them |
|---|---|---|---|---|
| B + PS (keeps 86.0%) | 226.4 | 3.8 | 5 of 6 | **×80** (3.8 to 305.6) |
| D = Cholis T3 + PS (70.8%) | 18.6 | 0.0 | 2 of 6 | ×46 (18.6 to 79.3) |
| A (8.5%) | 4.7 | 1.3 | **0 of 6** | cannot be measured |

Under Mask D, three splits including `lat` and `lat_ns` fail the admissibility test: with the model freed they cannot recover a signal added to the real data (boosted gain +0.3 to +1.7 against a threshold of +25). Under Mask A no split is admissible, not even the rigid one — the mask keeps 8.5% of the ROI and there is nothing left to measure with. Both belong with a ring-decomposed IEM or not at all, which is the constraint that goes in the Conventions.

## 4. What follows

**For Gap 1.** §6.1–6.5 are withdrawn. Every ΔAIC, the axis-ratio correlation, the mask comparison, the energy dependence and the triaxial ordering are computed at one arbitrary point on an eighty-fold range. Untouched: the PSF correction (§6.7), H₀ᴮ (§6.6) and the catalogue-incompleteness calibration (§6.6), all of which are differential measurements on simulated skies where the background is correct by construction.

**For every gap.** A TS quoted against a background model whose freedom has not been varied is not yet a measurement. The test needs a counts cube, a diffuse model and the templates under test; it does not need anything specific to Gap 1. See `analysis/diffuse_freedom/README.md`.

**For a gap whose IEM is already ring-decomposed** (Gap 5 by construction), the informative comparison is not adding latitude bands on top — it is running the same data through a single-template IEM and through the rings, and asking whether the spread closes. That is the measurement that says whether ring models remove this ambiguity or only hide it.

**For Gap 7.** The diffuse-freedom axis belongs in the covariance as a sixth dimension. Within Gap 1 alone it is the largest single source of variance, larger than mask, larger than point-source treatment, larger than energy binning.

## 5. Limits

- Latitude bands are our geometry, not gas astronomy. They are enough to show that the background choice dominates and not enough to say which choice is right. A ring decomposition is the fix, not a better choice of bands.
- One IEM family. The spread is measured within `gll_iem_v07`; the spread *between* families is a separate and probably larger number, and is what Families 2 and 3 are for.
- One ROI, one event selection, 1–10 GeV, Level 1 only.
- The admissibility thresholds (50% of the injected signal, +25 on the boosted sky) are judgement calls. They were fixed before the run and they are in the config, but a reader who prefers different ones should re-run rather than re-read the table.

## 6. Reproducing it

```bash
python analysis/diffuse_freedom/run_test.py \
  --ccube gap1_ccube_full.fits --expcube gap1_expcube_full.fits \
  --iem gap1_iem_roi.fits --catalog gll_psc_catalog.fit \
  --mask B --ps --out outputs/gap1/diffuse_freedom_maskB_ps

python analysis/gap1/make_diffuse_figure.py \
  --json outputs/gap1/diffuse_freedom_maskB_ps/diffuse_freedom.json \
  --labels "Mask B + PS, 1-10 GeV, VVV template" \
  --out outputs/gap1/fig8_diffuse_freedom.png
```

Full records including every template and every sky: `outputs/gap1/diffuse_freedom_*/diffuse_freedom.json`.
