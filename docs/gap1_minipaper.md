# Gap 1: Can a Baryonically Contracted Boxy Dark Matter Halo Mimic the VVV Stellar Bulge?

**Team:** Oren (Gap 1 lead), + project lead support. Open to Jason / Quentin.
**Status:** Sections 1–5 complete. Section 6 pending first results.
**Repo path:** `docs/gap1_minipaper.md` | **Analysis path:** `analysis/gap1/`
**Version:** draft v0.1 — Week 1
**Git commit:** _(insert hash at first results run)_

---

## 1. Specific Research Question

**Primary question.** Under a common high-precision analysis pipeline, can a boxy dark matter template derived from the baryonically contracted inner halos of the HESTIA constrained cosmological simulations reproduce the Galactic Center Excess as well as, or better than, the VVV red-clump stellar bulge template?

**Secondary question (the pipeline question).** Does iterative source-finding of the kind used by Di Mauro (2026) — which resolves new point sources until residuals fall below ~10% — itself bias the recovered GCE morphology toward sphericity by absorbing boxy asymmetry into newly "resolved" sources? Equivalently: **is the reported sphericity of the GCE a property of the sky, or a property of the pipeline?**

These are deliberately posed as one gap rather than two, because a Level-2 preference for VVV is only interpretable once we know whether the Level-2 pipeline can recover a boxy signal at all. Answering (1) without (2) is what the field has already done.

---

## 2. Null Hypothesis

We state two nulls, both of which must be tested; the second is the control on the first.

**H₀ᴬ — Morphological indistinguishability.**
The HESTIA boxy DM template and the VVV bulge template describe the GCE equally well. Operationally: the AIC difference between the two single-template fits satisfies |ΔAIC| ≤ 10 for the majority of the 3 masks × 3 IEM families × 3 HESTIA realizations grid, at both Level 1 and Level 2, with no consistent sign.

**H₀ᴮ — Pipeline fidelity (no morphological bias).**
The Level-2 iterative pipeline is morphologically unbiased. Operationally: for boxy HESTIA-shaped signals injected at realistic flux, the recovered boxiness statistic satisfies B_rec / B_inj = 1 within statistical uncertainty, and the number and total flux of newly resolved point sources in boxy-injection runs is statistically consistent with that in spherical-injection runs at matched injected flux.

**Rejection thresholds (proposed, for team sign-off).** ΔAIC ≥ 10 in a consistent direction across ≥ 7 of 9 mask × IEM combinations is treated as a real preference. ΔAIC between 4 and 10, or preference that flips sign with mask or IEM, is reported as "systematics-limited" and handed to Gap 7 rather than claimed. Fixing these thresholds *before* unblinding is deliberate: the 3 × 3 × 3 grid is large enough that post-hoc selection of the most favorable configuration would be trivially easy.

---

## 3. Data and Methods

**Dataset.** Fermi-LAT Pass 8, 16.5 yr (MJD 54682–60710), P8R3_SOURCE_V3, evclass=128, evtype=3, zenith < 90°, rocking angle < 52°, DATA_QUAL>0 && LAT_CONFIG==1, FL16Y transients excluded. ROI |l| < 20°, |b| < 20°, CAR projection in Galactic coordinates, 0.1° pixels, 400 × 400 grid. (Project Conventions §1–2.)

**Templates used**

| Template | Source | Role |
|---|---|---|
| HESTIA boxy DM (ρ²) | Muru et al. 2025; Hussein et al. 2025 — 3 realizations, standard + flipped orientation | Signal hypothesis under test |
| VVV stellar bulge | Coleman et al. 2019 — nuclear and boxy/peanut components separated where possible | Competing hypothesis |
| Spherical gNFW | r_s = 20 kpc, R_⊙ = 8.2 kpc, γ free in [0.8, 1.4] | Benchmark / reference morphology |
| Triaxial DM (Hu et al. 2026) | q₁ = 0.6, q₂ = 0.8, θ = 25°, standard + flipped | Cross-check; interface to Gap 4 |
| Galactic bar | Macias et al. 2019 gas-correlated | Confounder — must float, not be assumed absent |
| Fermi Bubbles | sharp- and fuzzy-edge variants from `templates/bubbles/` | Background, free normalization |
| Isotropic + IEM + FL16Y point sources | per IEM family | Background |

All templates normalized to unit integral over the ROI and convolved with the energy-dependent LAT PSF using the shared `utils/psf.py` implementation — no gap-local PSF code.

**Masks.** A (aggressive: all FL16Y sources at 95% containment + 0.5°), B (minimal: TS > 25 and flux > 5 × 10⁻⁹ ph cm⁻² s⁻¹), C (Mask A + latitude cut, |b| < 2° for |l| < 15° tapering to |b| < 0.5° at |l| = 20°). All three, all fits.

**IEMs.** All three families: (1) `gll_iem_v07.fits` + `iso_P8R3_SOURCE_V3_V7.txt`; (2) the 12 extreme McDermott et al. 2023 models selected by Gap 5; (3) `galp21.fits` + Pohl22 isotropic (exact filename pending Gap 5 confirmation). Each with isotropic + free-normalization Bubbles in both edge variants.

**Energy bins.** 0.5–1, 1–3, 3–10, 10–50 GeV. Fits performed per-bin and jointly; morphology conclusions must hold in the 1–3 and 3–10 GeV bins where the GCE dominates.

**Statistical method.** Poissonian binned likelihood, Fermipy 1.2 with MINUIT2, max 5000 iterations. Report TS and AIC for every template and every mask × IEM combination. Two analysis levels run on identical data:

- **Level 1 — single-pass template fitting.** Fixed FL16Y source list, no source addition. This is the field's standard comparison and our baseline.
- **Level 2 — iterative source-finding.** Replicates the Di Mauro (2026) procedure: fit → residual TS map → add point sources above threshold → refit, iterating until fractional residuals fall below the stopping criterion. The stopping rule and TS threshold are fixed in `configs/` before running and are *not* tuned per template — using different stopping behavior for the boxy and VVV fits would manufacture the very bias we are testing for.

**Software.** Python 3.12, Fermitools 2.5.3, numpy, scipy, astropy. (NPTFit / PyTorch not required for Gap 1 unless we extend to photon-count cross-checks.)

**Special analysis steps**

1. **HESTIA map preparation.** Obtain ρ² maps; perform line-of-sight integration from the solar position; project to the Galactic CAR grid at 0.1°; produce one FITS template per realization × orientation. Validate against any published J-factor or angular profile in the source papers before use.
2. **Orientation scan.** Beyond the standard/flipped pair, scan bar-angle orientation in steps across the plausible range, since the boxy signature is orientation-sensitive by construction.
3. **Injection–recovery (the core test).** Inject synthetic boxy HESTIA-shaped signals into simulated skies (seed = 42) at several flux levels bracketing the observed GCE, run them through the full Level-2 pipeline, and measure the recovered morphology.
4. **Spherical control injection.** Inject a spherical gNFW signal of matched flux through the identical pipeline. Bias is only demonstrated if boxy-in-spherical-out occurs *while* spherical-in-spherical-out holds — otherwise we are measuring a generic reconstruction failure, not a morphological one.
5. **Source-absorption bookkeeping.** For every Level-2 run, log the number, position, and flux of newly resolved sources. If boxy injections systematically produce more new sources along the major axis of the injected boxiness, that is the direct mechanism of the bias, not merely a symptom.
6. **Template correlation matrix.** Compute pixel-space correlations among HESTIA, VVV, gNFW, and the Macias bar template over each mask, and report the condition number. Two templates that are ~95% correlated within the ROI cannot be separated at any statistical significance, and we should say so rather than report a spurious ΔAIC.
7. **Results handoff.** Emit the full ΔAIC / TS grid in the conventions' file format to `outputs/gap1/` for Gap 7's systematic covariance computation.

---

## 4. Expected Results and Interpretation

| Outcome | Level 1 | Level 2 | Injection test | Interpretation |
|---|---|---|---|---|
| **A. Boxy DM competitive** | HESTIA ≈ or > VVV | HESTIA ≈ or > VVV | unbiased | The morphological argument for MSPs collapses. Boxiness ceases to be evidence for a stellar origin, and the DM interpretation is restored to parity. Strongest possible result for DM. |
| **B. Bulge robustly preferred** | VVV > HESTIA | VVV > HESTIA | unbiased | The MSP/stellar interpretation is substantially strengthened — it survives the most serious morphological challenge available to it. A genuinely decisive null for DM morphology. |
| **C. Pipeline-induced sphericity** | HESTIA ≈ VVV | VVV > HESTIA | **biased** (boxy in → spherical out) | Di Mauro (2026)'s spherical GCE is an artifact of iterative source-finding, not a measurement. This would be the most consequential outcome of the whole project: it invalidates a headline result and means high-precision pipelines must be validated on injected morphologies before their morphological conclusions are trusted. |
| **D. Degenerate** | \|ΔAIC\| small | \|ΔAIC\| small | unbiased | Morphology alone cannot discriminate DM from the stellar bulge at Fermi-LAT resolution. This retires a decade of morphology-based arguments on *both* sides and pushes the field toward spectra, photon-count statistics, and multiwavelength discriminants. A publishable and clarifying result, not a failure. |
| **E. Systematics-limited** | preference flips with mask/IEM | ditto | either | The apparent preference is a background-modeling artifact. Result is handed to Gap 7 for covariance treatment; no morphological claim is made. |

Note that outcomes C and D are, in our judgment, at least as likely as A or B, and both are more interesting than the binary "DM vs. pulsars" framing the field usually reports. We are not designing this to produce a winner.

---

## 5. Potential Pitfalls and Mitigations

| Pitfall | Why it matters | Mitigation |
|---|---|---|
| **HESTIA maps may not be publicly available** in projectable form | Blocks the entire gap at step 1; this is the critical-path risk for Week 1 | Contact Muru et al. / Hussein et al. immediately. Fallback: construct an approximate boxy template by fitting a parametric triaxial + boxiness profile to the published HESTIA figures, clearly labeled as approximate and never used for a headline claim. |
| Realization / orientation / baryonic-physics dependence | The boxy template is not unique; a single realization could over- or under-state the effect | Three independent realizations × orientation scan; report the full envelope of ΔAIC, not the best member. |
| ρ² normalization and line-of-sight integration errors | A normalization error is invisible in a fit with free amplitude but corrupts any cross-comparison and all J-factor statements | Independent recomputation by a second person; validate against published profiles; unit-integral normalization over ROI per conventions. |
| Simulation resolution near the Galactic Center | HESTIA softening length may exceed the angular scale we are fitting at 0.1° pixels | Quantify the effective angular resolution of each map; smooth all templates to a common floor; exclude the innermost region if it lies below the resolution limit, and state the exclusion radius. |
| Template degeneracy (HESTIA vs. VVV vs. bar vs. Bubbles) | Highly correlated templates produce ΔAIC values with no discriminating power | Correlation matrix and condition number reported alongside every ΔAIC; declare degeneracy explicitly where present. |
| Level-2 stopping criterion is somewhat arbitrary | The result could depend on where you stop iterating | Fix the criterion in `configs/` in advance; additionally report ΔAIC *as a function of iteration number* so the reader sees the trajectory rather than one endpoint. |
| Tuning the pipeline per template | Would silently create or erase the effect under test | Identical pipeline configuration for all templates; enforced by shared config files and reviewed in code review. |
| Single injected flux level | Bias may appear only at particular signal strengths | Inject at ≥ 3 flux levels bracketing the observed GCE. |
| Multiple comparisons across the 3 × 3 × 3 grid | Cherry-picking a favorable configuration is easy and tempting | Pre-registered thresholds (§2); require consistency across ≥ 7 of 9 mask × IEM combinations; publish the full grid including unfavorable entries. |
| AIC applied to non-nested models | HESTIA and VVV are not nested; naive likelihood-ratio tests are invalid | Use AIC and TS as specified in the conventions; reserve likelihood-ratio tests for genuinely nested cases (e.g. γ within the gNFW family). |
| Mask B leaves faint sources unmasked | Could inject spurious small-scale structure that the iterative pipeline then chases | Compare Mask B against A and C; treat any A/B/C disagreement as a systematics flag for Gap 7. |
| Shared PSF implementation could be wrong for everyone at once | A common-mode error would not show up as inter-team disagreement | Validate `utils/psf.py` against Fermitools `gtpsf` outputs before Week 2 and record the check in `utils/README.md`. |

---

## 6. Preliminary Findings

**Pending.** To be filled as results come in. Planned contents: the ΔAIC / TS grid over masks × IEM families × HESTIA realizations at both analysis levels; injection–recovery curves of recovered vs. injected boxiness; the source-absorption diagnostic (new-source count and flux vs. injected morphology); template correlation matrix; and residual maps for the best- and worst-fitting configurations.

---

## 7. References

**Boxy dark matter / simulations**
- Muru et al. (2025), *DM morphology in Milky Way simulations* — arXiv:2508.06314
- Hussein et al. (2025), HESTIA inner-halo shapes
- Hu, Cholis & Zhong (2026), *Generic triaxial DM halo* — arXiv:2602.20252

**Stellar bulge / morphology**
- Coleman et al. (2019), *Galactic bulge and nuclear bulge templates* — arXiv:1911.04714
- Macias et al. (2016), *Galactic Bulge Preferred Over Dark Matter* — arXiv:1611.06644
- Bartels et al. (2018), *GCE as a tracer of stellar mass in the bulge* — arXiv:1711.04756
- Ramirez et al. (2024), *Gaussian-process morphology* — arXiv:2410.21367
- Song et al. (2024), *Robust inference of GCE spatial properties* — arXiv:2402.05449

**Pipeline / high-precision analyses**
- Di Mauro (2026), *Precise Fermi-LAT GCE morphology and spectrum* — arXiv:2605.22913
- Di Mauro (2021), *11 Years of Fermi-LAT* — arXiv:2101.04694
- Zhong & Cholis (2024), *Robustness Against Masking* — arXiv:2401.02481
- McDermott, Zhong & Cholis (2023), *A Phantom Menace* — arXiv:2209.00006

**Systematics and mismodeling**
- Calore, Cholis & Weniger (2015), *Background systematics* — arXiv:1502.02805
- Leane & Slatyer (2020), *Spurious Point Source Signals* — arXiv:2002.12370
- Leane & Slatyer (2020), *The Enigmatic Galactic Center Excess* — arXiv:2002.12371
- Leane & Slatyer (2019), *Revival of the Dark Matter Hypothesis* — arXiv:1904.08436

**Context**
- Hooper & Goodenough (2011) — arXiv:1010.2752
- Daylan et al. (2016), *Characterization of the central Milky Way* — arXiv:1402.6703
- Fermi-LAT Collaboration (2017) — arXiv:1704.03910
- Fermi-LAT 16-year Source List (FL16Y) — arXiv:2602.22148

_Reference metadata (years, exact titles) to be verified against the shared library before merge into the final paper._
