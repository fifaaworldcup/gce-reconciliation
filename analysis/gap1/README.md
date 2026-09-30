# Gap 1: Flattened Dark Matter Morphology vs. the VVV Bulge

**Team:** Oren · **Updated:** 2026-09-27 · **Commit:** see git log for this file
**Status:** the pipeline is built, validated and running. **The Level-1 morphology results in §6.1–6.5 are withdrawn pending a ring-decomposed diffuse model — see §6.0, where the same data give TS 3.8 to 306 depending on how much the background may bend.** What stands: the method work (§6.6, §6.7), the catalogue-incompleteness calibration (§6.6), and the diffuse-freedom test itself (§6.0).

> **Framing — why this is no longer titled "boxy DM".** Muru's `dm2` maps are the square of a projected column density, (∫ρ ds)², not the J-factor ∫ρ² ds. Those are not proportional, so **until Muru confirms otherwise these are not annihilation templates and nothing here should be read as a statement about an annihilation signal.** What §6 tests is whether a *flattened dark-matter-shaped* morphology fits the GCE as well as the stellar bulge, which is well posed either way. `rho2_definition_test.py` shows boxiness is immune to the difference (|Δc₄| ≤ 0.001) but the radial profile differs by 65–70%, so every HESTIA ΔAIC here inherits it. Revert the title and this note if Muru confirms ∫ρ²ds.

*Masks, IEM families, energy bins, ROI and software versions follow Project Conventions §1–4 and §7 — not repeated here. Sections 1 and 2 as already in the doc.*

> **For Afeefa / the team — read §6.0 first, it supersedes points 1 and 3 below.**
> 0. **Our detection moves by a factor of eighty depending on how much the diffuse model is allowed to bend (§6.0).** Mask B + PS, VVV, 1–10 GeV: TS 306 if the background may bend away from the Centre, 226 with one normalisation as in §6.2, 8.2 if it may bend across the whole ROI, 3.8 with north and south separate. Controls rule out the test being at fault — a signal we inject keeps 80% under every non-degenerate split, and one added on top of the real data is still found at ΔTS +192 to +639. It is not that the GCE is mismodelling; it is that **with one all-sky interstellar template nothing in our setup picks between background models that disagree by eighty**, so no TS and no morphology can be quoted. §6.1–6.5 are withdrawn until a ring-decomposed IEM exists. **Pohl22 alone would fix it** — that archive has four H I rings, four H₂ rings and two dust residual maps, and we only kept the Bubbles file out of it.
> 1. ~~Mask A is a null result~~ — **still true but now for a different reason (§6.0):** Mask A leaves TS 4.7 before any diffuse freedom, so there was nothing there to discriminate with. Original text: **Mask A is a null result, reported as one (§6.3).** The aggressive convention of Di Mauro 2026 and Song et al. 2024 is insensitive to GCE morphology in our hands — every |ΔAIC| < 10 across all templates. It keeps 8.7% of the ROI at 4FGL-DR4 density. **Mask C inherits this.** Mask D (flux-scaled radius) added as a mask that *can* distinguish.
> 2. **`utils/psf.py` was up to 30% too narrow above 10 GeV** and is now driven by a real `gtpsf` table (§6.7). It inflated every TS by about 2× and invented a 0.5–1 GeV detection. Anyone using it for shape work must pull the new version.
> 3. ~~Withdrawn, see §6.0.~~ **HESTIA vs VVV is systematics-limited, not a preference.** It flips sign with point-source treatment (§6.2). Per §2 it goes to Gap 7 rather than being claimed. The headline is the axis-ratio correlation (§6.1), which survives every configuration.
> 4. **For Gap 4:** Hu's q₁/q₂ axis assignment decides whether triaxial beats spherical at all (§6.5).
> 5. **Level 2 built and run; H₀ᴮ answered (§6.6).** No morphology-specific rounding at either source count tested (~1 and ~29 added sources): the flattened-minus-spherical differential is −0.003 ± 0.020. **But the starting catalogue matters far more than the iteration does** — withholding 90% of 4FGL biases the recovered flattening by −0.12 to −0.32, in the same direction as §6.1's headline, so q ≈ 0.6 is now an upper bound. Sizing the real bias is Gap 5's NPTF question. **Also for Gap 5:** the iteration resolves the peak of a diffuse excess into point sources that are not in the simulated sky. Module is `analysis/level2/`, offered as shared.
> 6. **Mask names changed to match yours.** *Mask D* is now Cholis 2022 Table III as you specified — TS-split radii (θₛ for TS ≤ 49, θ_l = (10/3)θₛ above it, both energy-dependent) plus the |b| < 2° disc cut. Our old flux-scaled mask is now **Mask E**. Every Mask-D number in §6.2–6.4 predates the rename and refers to Mask E; they are labelled.
> 7. **Fermitools 2.5.3 vs Di Mauro's 2.2.0:** documented as a systematic, not matched. Rationale in §6.10.
> 8. Conventions §1 (16.7 → 16.50 yr) and §7 (Fermitools 2.2) — Afeefa is fixing. One more for the Conventions: **Daylan et al. 2016's axis ratio is the reciprocal of Di Mauro's and Cholis's.** Any table that pools them without inverting Daylan manufactures a trend. Worth a line in §11.

---

## 3. Data and Methods

### Templates

| Template | Source | Role |
|---|---|---|
| HESTIA DM (ρ²) | Muru et al. 2025, private communication — 6 galaxies × 4 viewing angles | Hypothesis under test |
| VVV stellar bulge | Coleman et al. 2019; nuclear and boxy/peanut split where possible | Competing hypothesis |
| Spherical gNFW | r_s = 20 kpc, R_⊙ = 8.2 kpc, γ free in [0.8, 1.4] | Benchmark; ties to Di Mauro's γ = 1.15 ± 0.02 |
| Triaxial DM | Hu et al. 2026: q₁ = 0.6, q₂ = 0.8, θ = 25° | Cross-check; interface to Gap 4 |
| Galactic bar | Macias et al. 2019, gas-correlated | Confounder — floated, never omitted |

All normalised to unit integral over the ROI and PSF-convolved via the shared `utils/psf.py`.

### The HESTIA inputs, as delivered

Files are `density_projection_{sim}_8192_halo{id}_angle{A}_dm2_{S}deg.{fits,tsv}`, on an 800 × 800 grid of pixel centres from −39.95° to +39.95° in 0.1° steps. Galaxy naming per Muru: G1.1 = 09_18/…003, G2.1 = 17_11/…003, G3.1 = 37_11/…002, G1.2 = 09_18/…002, G2.2 = 17_11/…002, G3.2 = 37_11/…001.

Four properties of these files, each of which corrupts the template silently if missed:

- **The maps are already squared.** `dm2` *is* the annihilation map. Squaring again gives ρ⁴.
- **Their `CDELT1` is +0.1**, ours is −0.1. Longitude must be flipped or the template comes out mirrored east–west.
- **The `.tsv` files are the transpose of the `.fits` files** — verified directly, `np.loadtxt(tsv) == fits.T`. Mixing them swaps l and b.
- Their pixel centres lie on the same 0.1° lattice as ours, so our ROI is an exact **crop** (indices 200:600) — no interpolation, no half-pixel offset.

Code: `analysis/gap1/hestia_to_template.py` (reader, regridder, shape statistics) and `analysis/gap1/project_template.py` (analytic gNFW builder and geometry self-tests). Both self-test with no data present.

### Analysis steps specific to this gap

1. **Two analysis levels on identical data.** *Level 1* — single pass, source list fixed. *Level 2* — iterative source-finding replicating Di Mauro (2026): fit → residual TS map → add sources above threshold → refit, with thresholds and stopping rule fixed in `configs/` in advance and applied identically to every template.
2. **Common resolution floor.** The HESTIA maps carry their own smoothing. Every competing template is convolved with the same top-hat before comparison.
3. **Shape measured on the uncropped ±40° source maps**, not inside the ROI (pitfall 1).
4. **Boxiness statistic** — a c₄ Fourier/superellipse measure rather than an axis ratio, calibrated against superellipses of known shape.
5. **Orientation sampling** across the four delivered viewing angles (0°, 30°, 60°, 90°).
6. **Injection–recovery** of HESTIA-shaped signals into simulated skies (seed = 42) at ≥ 3 flux levels, through the full Level-2 pipeline.
7. **Spherical control injection** at matched flux through the identical pipeline.
8. **Source-absorption bookkeeping** — position and flux of every source the iteration adds, to test whether new sources cluster along the major axis of the injected shape.
9. **Template correlation matrix** and condition number, published alongside every ΔAIC.
10. **Results handoff** to `outputs/gap1/` in the Conventions §11 format for Gap 7.

## 4. Expected Results and Interpretation

If HESTIA DM fits as well as or better than the VVV bulge, the MSP argument is weakened. If VVV is better, the MSP interpretation is strengthened. Four ways this can land:

- **DM competitive** — HESTIA fits as well as or better than VVV at both levels, and the injection test comes back clean. The morphological argument for MSPs collapses and the DM interpretation is back at parity.
- **Bulge robustly preferred** — VVV wins at both levels with a clean injection test. The MSP interpretation gets stronger, because it survives the best morphological challenge available to it.
- **Pipeline-induced sphericity** — the two templates tie at Level 1, VVV wins at Level 2, and the injection test shows the pipeline flattens shape. Then Di Mauro's spherical GCE is an artifact of method rather than a measurement, and high-precision pipelines need validating against injected morphologies before anyone trusts their shape conclusions. This is the most consequential outcome available to us.
- **Degenerate** — ΔAIC stays small everywhere. Morphology just can't discriminate at Fermi-LAT resolution, which retires morphology-based arguments on both sides and pushes the field toward spectra and photon-count statistics. Worth publishing, not a failure.

Our §6 results tilt this. The HESTIA halos are flattened but not boxy, while the VVV boxy/peanut bulge is genuinely boxy — so boxiness may still separate them even where axis ratio doesn't. That favours the second or fourth outcome, but by a route neither Di Mauro's *r* nor Muru's SVD axis ratio can test, since an axis ratio only ever describes an ellipse.

## 5. Potential Pitfalls and Mitigations

The first three aren't hypothetical — we hit them while building the templates.

- **A square ROI manufactures boxiness.** Once a contour reaches the ROI edge it gets clipped along the edges but not the corners, which injects a 4-fold signal — exactly what c₄ measures. Measured inside the ROI we got c₄ = +0.013 to +0.021 for all six galaxies, which reads as a confident boxy detection; the same contours on the uncropped map give −0.006 to +0.002. We nearly reported a false positive. Fix: measure every shape statistic on the uncropped source maps and use the crop only for the template that gets fitted. This applies to any gap measuring morphology, not just us.
- **Smoothing changes the answer.** G1.2's axis ratio runs 0.59 at 0.1° smoothing up to 0.79 at 3° — blurring makes everything rounder. So an unmatched resolution difference between templates doesn't add noise, it manufactures or erases the effect under test. Fix: one common top-hat floor on every template, and report how things move when that floor changes. Also the 0.1° and 0.3° maps are shot-noise dominated (25% and 18% pixel scatter) and not usable.
- **Don't measure below the simulation's resolution.** HESTIA's 220 pc softening is about 1.6° at 8 kpc, and the 50%-of-peak contour sits at roughly that radius. Fix: measure at 2–10% of peak, which probes 5–20°, and quote the contour radius with every number.
- **The ρ² definition is still open.** `dm2` is the square of a projected column density, (∫ρ ds)², not the J-factor ∫ρ² ds — confirmed from the upstream source. Those aren't proportional, so fitting the former and calling it an annihilation template mislabels the whole result. Question is out to Muru; meanwhile the convention is stamped into every FITS header and the reader refuses to assume one. Can also be cross-checked using the repo's own mock-NFW generator against an analytic J-factor.
- **Boxiness needs the right statistic.** Di Mauro scans an elliptical axis ratio and the upstream HESTIA analysis uses an SVD axis ratio; neither can express boxiness, since a shape can be strongly boxy at axis ratio 1. Using one would make this gap unable to detect the thing it exists to detect. Fix: the c₄ measure, frozen before unblinding, calibrated on superellipses where the answer is known (ellipse → 0.000, boxy n=4 → +0.019, disky n=1.3 → −0.025), with the axis ratio reported alongside for comparability.
- **The templates may be too similar to separate.** HESTIA, VVV, the Macias bar and a mildly elongated gNFW are all centrally concentrated blobs of similar size in the same box. Two templates correlated at ~95% inside the mask can't be separated at any significance, and a ΔAIC between them measures nothing. Fix: publish the correlation matrix and condition number beside every ΔAIC, and where they're near-degenerate report that honestly instead of a spurious preference.
- **Halo and orientation dependence.** One galaxy or one viewing angle could over- or under-state the effect, so we use all six galaxies across all four delivered angles and report the envelope rather than the best member.
- **Don't tune the Level-2 pipeline per template.** Different stopping behaviour for HESTIA and VVV would manufacture the exact bias we're testing for, undetectably. Thresholds and stopping rule go in `configs/` before any run, identical for every template, checked in review. We also report ΔAIC as a function of iteration number, not just at convergence.
- **Injection without a control proves nothing.** Non-spherical → spherical on its own is equally consistent with a pipeline that mangles every morphology. It only demonstrates directional bias alongside spherical → spherical at matched flux, and we inject at three or more flux levels since bias can show up only near the detection threshold.
- **Our Level 2 is a reimplementation, not a replication.** Di Mauro uses SOURCEVETO events, 0.08° pixels and 1–10 GeV, and our conventions differ. So any discrepancy could get blamed on event selection rather than on the injection test. Fix: run one Level-2 configuration matched to his selection as a validation anchor.

## 6. Preliminary Findings

### 6.0 Our background model choice moves the detection by a factor of eighty

**2026-09-28. This section supersedes the morphology conclusions in §6.1 to §6.5.**

**The worry.** Our background is one interstellar template — `gll_iem_v07` — with one free normalisation per energy bin. Four numbers across 1–10 GeV. Every published analysis of this excess decomposes the diffuse emission into Galactocentric rings and gas phases with independent normalisations; Di Mauro (2026) fits about twenty ring templates per IEM. §6.10 has said since the first round that ours leaves a north/south plane residual in every variant, and that was filed as a caveat rather than treated as a threat. If a rigid template mispredicts the plane, the flattest component in the fit gets recruited to fix it — and §6.1's headline was precisely that flatter templates fit better.

**The test** (`analysis/diffuse_freedom/`, shared module; full write-up in `docs/diffuse_freedom_test.md`). Refit with the diffuse model given freedom in *shape* as well as scale, in six ways, and watch the TS. Three skies each, because a test that destroys the signal proves nothing on its own: (a) the real sky, (b) a simulated sky built from our own fitted background plus a VVV-shaped signal at exactly the flux the rigid fit attributed to it, and (c) **the real data with one more real-strength signal added on top**, which keeps every real residual and asks whether a GCE-shaped excess is still findable in their presence.

**Table 11.** TS of the VVV bulge template, Mask B + PS, 1–10 GeV — the only configuration with any sensitivity (see below). "Added" is the *gain* from injecting one more real-strength signal into the real data.

| how the diffuse model may bend | extra params | observed | injected-only sky | added (ΔTS) |
|---|---|---|---|---|
| one normalisation per bin (as in §6.2) | 0 | 226.4 | 181.7 | +678.1 |
| north/south | 1 | 159.4 | 175.5 | +581.5 |
| four \|b\| bands, but only at \|l\| > 10° — away from the GCE | 4 | **305.6** | 148.1 | +638.5 |
| four \|b\| bands, but only at \|l\| ≤ 10° — on top of the GCE | 4 | 113.7 | **64.8** | +283.4 |
| four \|b\| bands across the whole ROI | 3 | 8.2 | 145.6 | +245.2 |
| those bands north and south separately | 7 | **3.8** | 143.7 | +191.8 |

**What it says.** The observed TS runs from **3.8 to 305.6** — a factor of eighty — across six background models, none of which could be ruled out before looking at the data. That is the headline, and it is worse than a single number going away.

The controls say which rows can be trusted. The injected-only sky keeps 80% of its signal under every split except the one restricted to \|l\| ≤ 10°, which keeps 36% — so *that* row is genuinely degenerate with the signal and cannot be used as evidence either way. The other splits are not degenerate: they keep four fifths of a signal we put there ourselves, and they still find a real-strength signal added on top of the actual data at ΔTS +192 to +639. So when the observed excess collapses to 3.8 under bands-everywhere, that is not the test failing to see a GCE.

**Where the ambiguity actually lives.** Letting the diffuse model bend *away* from the Galactic Centre makes the detection **stronger**, 226 → 306: the rigid model was mis-fitting the outer plane and dragging the inner fit with it. Letting it bend *across the whole ROI*, including through the longitudes the excess occupies, destroys it. The difference between those two is one question — is the inner-longitude latitude profile of the diffuse emission allowed to differ from what `gll_iem_v07` says? — and **our data cannot answer it with a single all-sky template.** One extra parameter, just north over south, already moves the observed TS by 30% while moving the injected one by 3%.

**So the honest position is not "the GCE was all mismodelling."** A real signal survives every non-degenerate split, and Di Mauro reports TS of order 10⁴ with ring-decomposed backgrounds. It is that **no TS and no morphology can be stated from this configuration**, because the background model choice dominates the measurement by a factor of eighty and nothing in our setup picks between the choices. Every ΔAIC in §6.2, the axis-ratio correlation in §6.1, the mask comparison in §6.3, the energy dependence in §6.4 and the triaxial ordering in §6.5 are computed at one arbitrary point on that range. Figure 8.

The direct axis-ratio fit makes the same point from the other side. Running the §6.6 q-ladder on the real data with the rigid background drives q to the bottom rung, q ≤ 0.30 under both Mask B + PS and Mask D + PS, against Di Mauro's r = 1.10 ± 0.05, i.e. q = 0.909 ± 0.041. A likelihood that wants the flattest template available whatever the floor is buying a correction to the diffuse model, not measuring a shape. `analysis/gap1/measure_q.py`.

**The other two masks have nothing to measure.** Mask A gives TS 4.7 before any freedom is added, so §6.3's null result is real but its cause is simpler than argued there: there was no signal in it to discriminate with. Mask D (Cholis) gives 18.6, and once the diffuse model is freed that configuration cannot detect even a signal we put there ourselves (ΔTS +0.3). A \|b\| < 2° cut plus a single-template background leaves no morphology information in our setup at all.

**What this means for the plan.** The fix is not more clever masking, and it is not a better choice of latitude bands — those are our arbitrary geometry standing in for gas astronomy. It is a diffuse model whose freedom is set by the gas distribution: a ring decomposition. **IEM families 2 and 3 were a robustness check and are now the critical path.** Pohl22 alone would do it, since that archive carries four H I rings, four H₂ rings and two dust residual maps (§6.9).

**Done, 2026-09-30 — and it changes the picture twice.** Pohl22 is in and regridded (`analysis/gap1/fit_rings.py`, §7 of `docs/diffuse_freedom_test.md`). Swapping the single template for sixteen ring components moves VVV from TS 226 to 867–2365 and rescues Mask D, which had no sensitivity at all before and now reaches 480–1048. Adding Pohl's nuclear bulge then collapses the dark-matter templates, gNFW 129 → 13 and HESTIA 79 → 4, while the stellar bulge holds at 620 — apparently Macias's result. **But injecting a gNFW signal at a fixed reference flux shows that configuration can no longer recover dark matter that is definitely there** (TS 16.5, below threshold), so the collapse is a loss of discriminating power rather than evidence against dark matter. A null obtained by adding background components means nothing without that control. The 18-component fit returns a negative injected TS and is not converging; reported, not interpreted.

**What survives untouched.** Everything that is about method rather than about the sky: the PSF correction (§6.7), the injection–recovery machinery and the H₀ᴮ answer (§6.6), and the catalogue-incompleteness calibration (§6.6, Table 10). Those are differential measurements on simulated skies where the background is exactly right by construction, so a background error cannot bias them. The Level-2 module, the mask builders, the fitter and its self-tests are unaffected.

**And the test is itself a deliverable.** §2 frames this project around the Macias–Di Mauro contradiction: the same templates and gas maps giving opposite conclusions, with residuals of 30–40% in one case and ≤10% in the other. Table 11 is that disagreement reproduced inside one pipeline, and it says the disagreement does not need different data or different templates — a defensible change in how much the background may bend is enough to move a detection by eighty. Any team quoting a TS should run it. It takes a counts cube, a diffuse model and the templates under test, and needs nothing specific to Gap 1. Splits live in `configs/diffuse_freedom.yaml`, fixed before any run and identical for every template, on the same footing as `configs/level2.yaml`.

### 6.1 Headline result

> **Withdrawn 2026-09-28 (§6.0).** The correlation below is computed at one arbitrary point on a background-model range that moves the detection by a factor of eighty, and a flattened template is exactly what a rigid diffuse model recruits to fix its plane residual. Kept as written so the correction is visible. Re-run when a ring-decomposed IEM lands.

**The Fermi data select a flattening, and do not care where it comes from.**

Each HESTIA map's projected axis ratio q, measured on the raw template before any fit (Table 1), predicts how well it fits relative to the VVV bulge. Spearman ρ across all 24 halo × angle combinations:

| cell | ρ(q, ΔAIC vs VVV) | p | n |
|---|---|---|---|
| Mask B + PS, 1–10 GeV | **+0.73** | 5×10⁻⁵ | 24 |
| Mask B + PS, 1–3 GeV | +0.68 | 2×10⁻⁴ | 24 |
| Mask B + PS, 3–10 GeV | +0.67 | 4×10⁻⁴ | 24 |
| Mask D + PS, 1–10 GeV | +0.51 | 1×10⁻² | 24 |
| Mask D, 1–10 GeV | +0.59 | 3×10⁻³ | 24 |

Positive and significant in every cell, in both energy bins where the GCE is detected, under two different masks and with point sources either masked or modelled. Flatter halos fit better, and the crossover sits at q ≈ 0.6 — VVV's own axis ratio. Boxiness (c₄) shows nothing consistent: its sign flips between contour levels.

This is the result, not the head-count of how many HESTIA halos beat VVV. It also **explains** that head-count: the split is by galaxy because the galaxies have different axis ratios, and G3.1 and G3.2 lose everywhere because they are the two roundest (q = 0.69 and 0.63 against 0.55–0.61 for the rest).

Two things follow. A dark matter halo flattened like the bulge fits the GCE as well as the bulge does, so morphology alone does not separate them — H₀ᴬ is not rejected. And the quantity the data actually constrain is q, not "DM vs stars", which is a cleaner thing to report than a template ranking.

Caveats that belong with it: the four viewing angles per galaxy are not independent, so the honest sample is ~6 galaxies, where ρ ≈ 0.77 at p ≈ 0.07 — suggestive, not conclusive. This is one IEM of three. And **§6.6 found that an incomplete source catalogue biases the recovered flattening in exactly this direction** — worse, that past a certain point it drives *any* injected morphology to q ≈ 0.59, which is the number measured here. That is now calibrated rather than open-ended (Table 10): the correction depends on one quantity, the unresolved point-source flux relative to the GCE, which Gap 5 measures. Below about a third of the GCE flux the correction is smaller than the scatter and q ≈ 0.6 stands as a measurement; above about the GCE flux it does not stand at all. **Until Gap 5 supplies that number, quote q ≈ 0.6 with the conditional attached, not bare.**

### 6.2 The fit

> **Withdrawn 2026-09-28 (§6.0).** Every TS and ΔAIC here is measured against a single-template background. The same data give TS 3.8 to 306 across six defensible background models, and nothing in our setup picks between them. Numbers kept for the record and for comparison after the rerun.

Level-1 template fit, 16.5 yr of P8R3_SOURCE_V3, 1–10 GeV, 4.93 M photons. Background: `gll_iem_v07` + isotropic + fuzzy Bubbles, free normalisation per component per bin; signal likewise. Binned Poisson likelihood, EM then Newton to a provable lnL gap < 10⁻⁶ in every bin. Same parameter count for every model, so ΔAIC = −2ΔlnL. HESTIA below is G1.1 at 0°; the other 23 are in §6.1 and Figure 4.

**All numbers use the real LAT PSF** (`gtpsf` on our own livetime cube, §6.7). Nothing here comes from the earlier analytic approximation.

**Table 3.** Detection and normalisation. Flux is the fitted signal over 1–10 GeV and the whole ROI.

| mask (sky kept) | template | TS | σ | flux 1–10 GeV (10⁻⁷ ph cm⁻² s⁻¹) | residual |
|---|---|---|---|---|---|
| A (8.7%) | gNFW γ=1.2 | 2.9 | 1.7 | 0.32 | 1.013σ |
| A | HESTIA | 5.5 | 2.3 | 0.39 | 1.013σ |
| A | VVV | 9.8 | 3.1 | (degenerate) | 1.013σ |
| B (85.9%) | gNFW γ=1.2 | 1,402 | 37.4 | 2.70 | 1.224σ |
| B | HESTIA | 2,166 | 46.5 | 3.30 | 1.221σ |
| B | VVV | 1,986 | 44.6 | 2.29 | 1.221σ |
| B + PS | gNFW γ=1.2 | 121 | 11.0 | 0.60 | 1.052σ |
| B + PS | HESTIA | 248 | 15.8 | 0.95 | 1.051σ |
| B + PS | VVV | 203 | 14.2 | 0.64 | 1.051σ |
| D (86.4%) | gNFW γ=1.2 | 1,354 | 36.8 | 1.51 | 1.098σ |
| D | HESTIA | 1,609 | 40.1 | 2.39 | 1.097σ |
| D | VVV | 1,794 | 42.4 | 1.62 | 1.096σ |
| **D + PS** | gNFW γ=1.2 | 9.2 | 3.0 | 0.10 | 1.060σ |
| **D + PS** | HESTIA | 230 | 15.2 | 0.89 | 1.059σ |
| **D + PS** | VVV | 143 | 12.0 | 0.46 | 1.060σ |

**Table 4.** Model comparison. Negative ΔAIC favours the first-named template.

| mask | comparison | ΔlnL | ΔAIC | decisive? |
|---|---|---|---|---|
| A | HESTIA vs gNFW | +1.3 | −2.6 | no |
| A | VVV vs gNFW | +3.5 | −7.0 | no |
| A | HESTIA vs VVV | −2.2 | +4.4 | no |
| B | HESTIA vs gNFW | +381.7 | −763.5 | yes |
| B | VVV vs gNFW | +292.0 | −584.0 | yes |
| B | HESTIA vs VVV | +89.7 | −179.4 | yes |
| B + PS | HESTIA vs gNFW | +63.7 | −127.3 | yes |
| B + PS | VVV vs gNFW | +41.2 | −82.5 | yes |
| B + PS | HESTIA vs VVV | +22.4 | −44.9 | yes |
| D | HESTIA vs gNFW | +127.5 | −254.9 | yes |
| D | VVV vs gNFW | +219.8 | −439.7 | yes |
| D | HESTIA vs VVV | **−92.4** | **+184.7** | yes (favours VVV) |
| D + PS | HESTIA vs gNFW | +110.3 | −220.7 | yes |
| D + PS | VVV vs gNFW | +67.1 | −134.2 | yes |
| D + PS | HESTIA vs VVV | **+43.3** | **−86.5** | yes (favours HESTIA) |

**What is robust and what is not.**

- **Robust — the spherical gNFW loses.** Every mask with a detection, every γ in [0.8, 1.4], all 24 HESTIA halos, both energy bins. ΔAIC against VVV is +584 (B), +113 (B + PS), +440 (D), +134 (D + PS). Best γ is 1.1–1.2 wherever the GCE is detected, matching Di Mauro's 1.15 ± 0.02. Passing the gNFW through Muru's own 1° top-hat makes it *worse* at every γ, so this is not a smoothing artifact.
- **Not robust — HESTIA vs VVV.** The preference *flips with how point sources are treated*: under Mask D, masking them favours VVV by ΔAIC +185 (0 of 24 HESTIA halos beat VVV), modelling them favours HESTIA by −87 (18 of 24 beat VVV). Under Mask B + PS it is a near-even split, 9 of 24. By §2's own rule a preference that flips sign with the analysis configuration is **systematics-limited and is not claimed**; it goes to Gap 7. Figure 6.
- **Robust — VVV beats the older F98 boxy bulge** in every cell (ΔAIC +9 to +89), consistent with Coleman et al. 2019.

**Table 5.** Shape, cropped and uncropped (Mask B + PS). *Fitted* is the model signal map the likelihood saw, after PSF and exposure. A cropped contour reaching the ROI edge is refused, not reported, because that is where a square crop fakes boxiness (Table 2).

| template | contour | q fitted | c₄ fitted | q cropped | c₄ cropped | q uncropped | c₄ uncropped |
|---|---|---|---|---|---|---|---|
| HESTIA G1.1 0° | 30% | 0.693 | +0.0009 | 0.698 | −0.0012 | 0.698 | −0.0012 |
| | 20% | 0.692 | +0.0094 | 0.674 | +0.0107 | 0.674 | +0.0107 |
| | 10% | 0.621 | −0.0016 | 0.630 | −0.0010 | 0.630 | −0.0010 |
| | 5% | refused | | refused | | 0.532 | −0.0064 |
| VVV | 30% | 0.601 | +0.0062 | 0.592 | +0.0056 | — | — |
| | 20% | 0.599 | +0.0094 | 0.590 | +0.0086 | — | — |
| | 10% | 0.601 | +0.0123 | 0.597 | +0.0122 | — | — |

Where a contour sits inside the ROI, cropped and uncropped agree exactly — the radius guard working. PSF and exposure weighting move q by ≤ 0.02.

### 6.3 Masks: A is a null result, E is the replacement

> **Renamed 2026-09-28.** Afeefa specified *Mask D* = Cholis 2022 Table III (TS-split radii plus a |b| < 2° disc cut), so what this section calls Mask D is now **Mask E** and every Mask-D number in §6.2 to §6.4 is a Mask-E number. Mask D as now implemented is in `masks/make_masks.py` and its results are in §6.0's Table 11 only.
> **Partly superseded (§6.0).** Mask A's insensitivity is confirmed but its cause is simpler than argued here: the mask leaves TS 4.7 before any diffuse freedom is added, so there was no signal in it to discriminate with.

**Mask A, reported as a null result.** The aggressive masking convention — all catalogue sources at the 95% containment radius + 0.5°, as used by Di Mauro (2026) and Song et al. (2024) — **is insensitive to GCE morphology in our hands: every |ΔAIC| < 10 across all templates** (Table 4). With 4FGL-DR4 that is 761 discs over the ROI; it keeps 8.7% of the sky averaged over 1–10 GeV and 0.0% in the lowest bin, and the inner degrees where the morphology is decided are gone entirely. This is a finding about the convention, not a failure of the data: it says a fixed-radius aggressive mask cannot be used to make morphological statements at 4FGL-DR4 source density. It got worse, not better, once the true PSF replaced the too-narrow approximation, because the real 95% containment is wider. **Mask C inherits the limitation**, since it is Mask A plus a latitude cut.

**Mask D (new).** A mask whose radius is set by each source rather than by one number for all of them. Mask out to where the source's PSF-spread surface brightness falls to a fraction *f* of the diffuse background **at that source's position**:

    F_s(E) · PSF_E(r) = f · B_E(l_s, b_s)

with F_s the source's flux in that bin from its own 4FGL spectrum and B the diffuse model there. Both sides are ph cm⁻² s⁻¹ sr⁻¹, so the criterion is physical rather than tuned. Bright sources get wide discs, faint ones essentially none, the radius shrinks with energy through the PSF, and a source on the bright plane must be brighter to earn the same disc as one at high latitude. Code in `masks/make_masks.py`, self-tested against the PSF inversion.

**f = 0.1 was fixed before any Mask D fit**, on the contamination argument "mask where a source contributes ≥10% of the local background", not on an outcome. Sensitivity, for transparency:

| f | 1.0 | 0.3 | **0.1** | 0.03 | 0.01 |
|---|---|---|---|---|---|
| ROI kept | 98% | 95% | **87%** | 70% | 50% |

At f = 0.1 Mask D keeps 86.4%, comparable to Mask B's 85.9%, and unlike Mask A it is decisive: every comparison in Table 4 exceeds |ΔAIC| = 10. It is also what exposed the point-source-treatment dependence in §6.2, which a single mask would have hidden.

### 6.4 Energy dependence (Conventions bins)

Same fit (Mask B + PS, 40 templates) per Conventions bin, on a 0.5–50 GeV cube built from the same events; its 1–10 GeV subset reproduces the original cube's 4,933,089 photons exactly. Livetime cube reused, exposure sampled at 21 energies.

**Table 6.**

| bin | VVV TS | HESTIA TS (min/med/max) | beat / lose VVV | best gNFW ΔAIC (γ) | ρ(q, ΔAIC) | PS norm |
|---|---|---|---|---|---|---|
| 0.5–1 GeV | 0.0 | 0 / 0 / 0 | — | — | (n/a) | 0.80, 0.86 |
| **1–3 GeV** | **178.3** | 128 / 166 / 222 | 4 / 13 | **+79.9 (1.2)** | **+0.68** | 0.86, 0.85 |
| **3–10 GeV** | **48.1** | 24 / 52 / 69 | 8 / 5 | **+29.0 (0.9)** | **+0.67** | 0.88, 0.81 |
| 10–50 GeV | 0.0 | 0 / 0 / 0 | — | — | (n/a) | 0.56, 0.52 |
| 0.5–50 GeV | 226.4 | 165 / 216 / 286 | 9 / 12 | +112.6 (1.1) | +0.73 | — |

- The §6.1 correlation holds in both bins where the GCE is detected, which §3 requires of any morphology claim (Figure 5).
- **With the real PSF the GCE is detected only between 1 and 10 GeV.** Outside that every template gets TS 0, and 0.5–50 GeV reproduces 1–10 GeV exactly because the outer bins carry no signal. Under the old too-narrow PSF the 0.5–1 GeV bin appeared to hold a TS-19 VVV signal; it does not. Plausible — the GCE spectrum peaks near 2 GeV and at 0.5 GeV the PSF is r68 = 1.25° — but **a non-detection in a bin the Conventions ask for should be reproduced by an independent pipeline before it goes in the paper.**
- **Exposure sampling is a ~10% systematic:** the same photons with exposure at 9 vs 21 energies move VVV's TS from 203 to 226. Below every effect reported here, not zero.

### 6.5 Triaxial DM cross-check (interface to Gap 4)

Conventions give Hu et al. (2026) as q₁ = 0.6, q₂ = 0.8, θ = 25°, standard and flipped, without saying which axis each q is. Since the data select *vertical* flattening (§6.1), that assignment decides the outcome, so both readings were built (`build_triaxial.py`: gNFW γ = 1.2, ∫ρ²ds, same sampler and 15 kpc cut as every other template). Self-tests: q = 1 reproduces the spherical template to 3×10⁻¹¹; in-plane rotation leaves an axisymmetric halo unchanged; q_vert = 0.6 projects to q = 0.609.

**Table 7.** 1–10 GeV, Mask B + PS. Negative ΔAIC beats VVV.

| template | projected q (5%) | TS | ΔAIC vs VVV |
|---|---|---|---|
| HESTIA G2.1 90° | — | 286.5 | −60.1 |
| HESTIA G1.1 0° | 0.532 | 270.7 | −44.3 |
| VVV | 0.60 | 226.4 | 0 |
| Triaxial, vertical 0.6 / in-plane 0.8, flipped | 0.712 | 167.7 | +58.7 |
| Triaxial, vertical 0.6 / in-plane 0.8, standard | 0.712 | 155.2 | +71.2 |
| Spherical gNFW, γ = 1.1 | 1.000 | 113.7 | +112.6 |
| Triaxial, vertical 0.8 / in-plane 0.6, flipped | 0.839 | 103.0 | +123.4 |
| Triaxial, vertical 0.8 / in-plane 0.6, standard | 0.839 | 87.6 | +138.8 |

**The axis assignment decides whether triaxial beats spherical at all.** With 0.6 vertical it beats spherical by ~54 in ΔAIC; with 0.6 in the plane it loses to spherical. That is §6.1 again, from Gap 4's own template, ordered exactly by projected axis ratio. "Flipped" beats "standard" by 13–15 in both readings (our sign convention, to be matched to Hu's). No triaxial reading reaches VVV or the flatter HESTIA halos. Same ordering over 0.5–50 GeV.

### 6.6 H₀ᴮ: does iterative source-finding round off the morphology?

**Answer, at the level this test can reach: no.** Δq = **+0.001 ± 0.002** for a flattened injection, against **−0.001 ± 0.001** for a spherical control. H₀ᴮ is not rejected.

**The test** (`analysis/gap1/injection_recovery.py`). Simulate a sky from our own fitted background (diffuse + isotropic + Bubbles + the 4FGL point-source template, at the normalisations the real fit returned) plus a signal of known shape. Run the full Level-2 iteration on it. Measure the recovered flattening *before* and *after* the iteration on one fixed ruler: a ladder of gNFW-ρ² templates with vertical axis ratio q = 0.4 … 1.0, fitted one at a time, best q by lnL with a parabolic refinement between rungs. Five noise realizations per configuration; the pipeline is given the 4FGL catalogue, so the iteration only adds sources *beyond* it, as Di Mauro's does. Stopping rule from `configs/level2.yaml`, identical for every injection.

The spherical control carries the argument. "Flattened comes back rounder" on its own is equally consistent with a pipeline that mangles every morphology, so a spherical injection at matched flux goes through the identical pipeline and must come back unshifted.

**Table 8.** Recovered flattening before and after iteration. Errors are the standard error over 5 realizations.

| injected | flux | q level 1 | q level 2 | shift | sources added | signal flux lost |
|---|---|---|---|---|---|---|
| flattened (HESTIA G1.1 0°) | 1× | 0.642 | 0.640 | −0.002 ± 0.002 | 0.4 | 0.5% |
| flattened | 2× | 0.659 | 0.663 | +0.003 ± 0.003 | 1.2 | 2.3% |
| spherical control (gNFW γ=1.2) | 1× | 0.893 | 0.891 | −0.003 ± 0.003 | 0.4 | 4.7% |
| spherical control | 2× | 0.982 | 0.982 | −0.000 ± 0.000 | 1.0 | 4.1% |
| **pooled** | | | | **flattened +0.001 ± 0.002, control −0.001 ± 0.001** | | |

Reading the ladder: it recovers a flattened injection at q ≈ 0.64 rather than HESTIA's own 0.53, and a spherical one at 0.89–0.98 rather than 1.00, because the ruler's radial profile is gNFW's and the signal-to-noise is finite. That is why the quantity reported is the *shift* between levels on the same ruler, not the absolute value.

**What the iteration actually does.** It does not reshape the excess, but it does bite a piece out of the middle of it. Of 15 sources added across 20 runs, **14 landed within 2° of the Galactic Centre** (median 0.35°) — the iteration converts the peak of a diffuse excess into a "resolved" point source. It takes 0.5–4.7% of the signal flux with it, and more from the cuspy spherical injection (4–5%) than the extended flattened one (0.5–2%), which is what you would expect if what is being resolved is the cusp. **This is Gap 5's phantom-source question showing up inside Gap 1's test**, and it is worth their attention: these sources are not in the simulated sky, and the pipeline finds them at TS 38–154.

**No major-axis clustering detectable yet.** Of the 8 sources added under flattened injections, 7 sit within 0.6° of the centre and one at l = +7.4°. Mean |l| vs mean |b| is 1.15 vs 0.40 *only because of that one outlier*; excluding it, 0.26 vs 0.32 — no preference. The clustering test promised in §3 needs many more added sources than this configuration produces.

**At realistic source counts: still no morphology bias, but catalogue incompleteness is a large one.** The run above hands the pipeline the whole 4FGL catalogue, so its iteration only adds ~1 source. Di Mauro's resolves far more. Repeating the test with only the **brightest 10% of the catalogue given to the pipeline** — 643 sources present in the simulated sky but withheld from it — brings the iteration to 29 added sources per run, and changes the picture:

**Table 9.** Sparse-catalogue run, 1× flux, 3 realizations, **paired**. The simulated sky does not depend on how much of the catalogue the pipeline is shown, so each sparse run is differenced against the full-catalogue run on the *identical* sky and the Poisson noise cancels rather than being averaged over.

| injected | reference q | level-1 q | level-1 bias | level-2 q | shift | sources added |
|---|---|---|---|---|---|---|
| flattened | 0.642 | 0.521 | **−0.092 ± 0.029** | 0.567 | +0.046 ± 0.012 | 29.7 |
| spherical control | 0.893 | 0.574 | **−0.292 ± 0.042** | 0.623 | +0.049 ± 0.016 | 28.3 |

*Corrected 2026-09-26.* The first version of this table differenced the mean of 3 sparse realizations against the mean of 5 full-catalogue ones and quoted −0.120 and −0.320 with no error. Level-1 q scatters by ±0.05 (flattened) and ±0.08 (spherical) between realizations, so an unpaired difference of two small samples is largely noise. Pairing the identical skies gives a smaller bias and, for the first time, an uncertainty on it. Nothing downstream changes.

Two things, and the second is the important one.

- **The iteration still shows no morphology-specific bias.** Both injections shift by the same amount: the difference between them is **−0.003 ± 0.020**, consistent with zero. And the shift is *toward* the reference, not away — iterative source-finding is partially undoing a bias, not creating one. H₀ᴮ is not rejected at 29 added sources per run, which is the regime the question was about.
- **An incomplete source catalogue makes the excess look flattened**, by −0.12 for a flattened injection and −0.32 for a spherical one. Unmodelled point sources cluster along the Galactic plane, so the emission left over is plane-hugging and the fit absorbs it into the signal template as flattening. The iteration recovers only about a third of that for the flattened case and a sixth for the spherical one, and the sources it adds take 15% of the signal flux with them.

**This one cuts at our own headline.** §6.1 reports that the data select a flattening of q ≈ 0.6. Catalogue incompleteness biases exactly that quantity, in exactly that direction, and 4FGL-DR4 is not complete in the inner Galaxy. The 10%-catalogue test is a deliberate extreme and is not a estimate of the real bias — but it establishes the sign, and the sign is unfavourable: **the true excess may be rounder than we measure it.** Bounding this properly needs the real source-count distribution below the 4FGL threshold, which is Gap 5's NPTF territory, not something this test can supply. Until then §6.1's q ≈ 0.6 should be read as an upper bound on flattening, not a measurement of it.

**The calibration curve (new, 2026-09-26).** Two points are a bracket, not a bound, and "the catalogue is 10% complete" is not a quantity anyone measures. So the sweep was extended to thirteen catalogue fractions and re-expressed against something an NPTF analysis does report: the flux in point sources the fit was not given, as a fraction of the GCE flux. Call it *f*. `analysis/gap1/incompleteness_curve.py`, 5 realizations, level 1 only, 130 fits.

Two things make this a clean differential measurement. The simulated sky does not depend on how much of the catalogue the fit is shown, so one sky per realization is reused across every point on the curve and the noise is common to all of them. And the bias is a level-1 effect — the fit absorbing unmodelled sources into the signal template — so no source-finding is run and each point costs ten seconds instead of four minutes.

**Table 10.** Bias in recovered q against hidden flux. Errors are the standard error of the paired difference over 5 realizations. Reference q (f = 0) is 0.642 flattened, 0.893 spherical — reproducing Table 8 exactly, which is the check that this run and that one are the same experiment.

| f = hidden flux / GCE flux | sources hidden | bias, flattened | bias, spherical |
|---|---|---|---|
| 0.013 | 15 | −0.001 ± 0.001 | +0.002 ± 0.001 |
| 0.042 | 36 | −0.018 ± 0.002 | −0.006 ± 0.002 |
| 0.108 | 72 | −0.017 ± 0.003 | −0.017 ± 0.005 |
| 0.285 | 143 | −0.008 ± 0.005 | −0.033 ± 0.010 |
| 0.518 | 215 | −0.052 ± 0.012 | −0.100 ± 0.027 |
| 0.804 | 286 | −0.036 ± 0.013 | −0.124 ± 0.013 |
| 1.157 | 357 | −0.094 ± 0.019 | −0.233 ± 0.020 |
| 1.845 | 465 | −0.056 ± 0.018 | −0.218 ± 0.023 |
| 2.458 | 536 | −0.069 ± 0.020 | −0.255 ± 0.029 |
| 3.344 | 607 | −0.054 ± 0.020 | −0.246 ± 0.030 |
| 4.045 | 643 | −0.115 ± 0.022 | −0.311 ± 0.030 |
| 5.135 | 679 | −0.093 ± 0.023 | −0.286 ± 0.031 |

**The curve saturates, and that matters more than the slope.** Fitted over f ≤ 1.2 — the regime where an NPTF measurement of the inner Galaxy lives — the bias is **−0.079 ± 0.030 × f** for a flattened injection and **−0.168 ± 0.014 × f** for a spherical one. Fitted over the whole sweep instead, one straight line gives −0.028 and −0.086, a factor of 2–3 shallower, because the far points have already bottomed out. A global linear fit would therefore have *under-corrected* §6.1 by a factor of three, in the direction that flatters our own result. The saturating form −A f^p gives p = 0.66 and 0.83. The table, not either fit, is what should be used; `correct()` in the script interpolates it. The f ≤ 1.2 cut was chosen after seeing the curve bend, which is why both numbers are reported.

**The thing that should worry us: q ≈ 0.6 is an attractor.** A flattened and a spherical injection start 0.25 apart on this ruler. As hidden flux grows the gap closes monotonically — 0.25, 0.20, 0.16, 0.11, 0.09, 0.07, 0.06 — and both converge near **q = 0.59**. An incomplete catalogue does not merely bias the answer; past f ≈ 2 it erases the difference between a flattened and a round excess and returns 0.6 either way. Section 6.1 measures 0.60.

| f | q, flattened truth | q, spherical truth | gap |
|---|---|---|---|
| 0 | 0.642 | 0.893 | 0.252 |
| 0.29 | 0.634 | 0.860 | 0.227 |
| 0.80 | 0.606 | 0.769 | 0.163 |
| 1.16 | 0.547 | 0.660 | 0.113 |
| 2.46 | 0.572 | 0.639 | 0.066 |
| 5.14 | 0.549 | 0.607 | 0.059 |

Figure 7.

**So what does this do to §6.1?** It converts an open-ended worry into a one-number dependency, and the news is better than the bracket suggested:

| Gap 5's unresolved fraction f | correction to add to measured q |
|---|---|
| 0.05 | +0.02 |
| 0.10 | +0.02 |
| 0.20 | +0.01 to +0.03 |
| 0.30 | +0.01 to +0.04 |
| 0.50 | +0.05 to +0.10 |
| 1.00 | +0.07 to +0.19 |

Below f ≈ 0.3 the correction is smaller than the ±0.05 realization scatter on q itself, so **if Gap 5's unresolved point-source flux is below about a third of the GCE flux, §6.1's q ≈ 0.6 is a measurement and not an upper bound.** Between 0.3 and 1 it is a measurement with a stated correction. Above about 1 the measurement degenerates toward the attractor and should not be quoted at all. Published NPTF estimates put unresolved point sources at or below the GCE flux in this region, so the middle row is the expected case — but that is someone else's measurement to make on our ROI, mask and energy range, not a literature value to borrow.

**Does the literature show it? One paper does, and it is the only controlled test available.**

If an incomplete catalogue fakes flattening, then analyses that resolve fewer point sources should report a more plane-hugging excess. Across papers this is untestable and the raw ordering actually runs the other way: only four report an axis ratio at all, and going Daylan 2016 (2FGL) → Di Mauro 2021 (4FGL) → Cholis 2022 (4FGL-DR2) → Di Mauro 2026 (4FGL-DR4 plus ~30–50 iteratively found sources) the reported elongation along the plane *increases* with depth, 0.77 → 0.8–1.2 → 1.0–1.4 → 1.10 ± 0.05. The confound is obvious and fatal: those papers mask the plane differently (|b| < 1°, none, |b| < 2°, |b| < 2°) and Cholis 2022's own Appendix A shows that changing the longitude range of the mask moves ε by as much as the whole spread being compared.

The one place where source depth is varied with everything else held fixed is **Zhong & Cholis 2024** (arXiv:2401.02481), and it agrees with us. Same data, same mask geometry, same pipeline: with standard 4FGL point-source masks they find ε ≈ 0.9–1.5; with wavelet-based masks that remove *additional* source-like structure beyond the catalogue, ε drops to **0.4–0.8**. Removing more unresolved source-like emission makes the excess markedly less plane-hugging, which is §6.6's effect measured by someone else in a different code. What Table 10 adds is a calibration: how much shift per unit hidden flux, rather than a single before-and-after.

Two cautions. Their wavelet masks preferentially cut along the plane, because that is where the sources are, so part of that shift is geometric rather than astrophysical — the same confound as everywhere else in this subsection. And a *convention hazard worth putting in the Conventions doc*: Daylan et al. define the axis ratio as the reciprocal of Di Mauro's and Cholis's. Pooling them without inverting Daylan manufactures a trend in exactly the direction one might want.

**The clean test, for whoever gets there first.** Fix the mask, then vary only the source depth. Di Mauro's own pipeline already does three iterative passes at TS > 100, 64 and 36; reading the axis ratio off each pass would give an ε-versus-resolved-sources curve at fixed geometry, which is the measurement the literature is missing. Our Level-2 module uses the same three thresholds, so we can run it ourselves — after §6.0's diffuse-model problem is fixed, not before.

**What this curve is not.** The hidden sources are the faintest members of 4FGL-DR4, which are still brighter than the genuinely unresolved population below the detection threshold. A real unresolved population of the same total flux is made of more, fainter, more smoothly distributed sources, and would plausibly bias *less* per unit flux than this test does, because it looks more like diffuse emission and less like lumps. The curve is therefore likely conservative, but that is an argument, not a measurement, and it is the reason the x-axis is flux rather than source count.

**The honest limit on this result.** The two runs bracket the source count: ~1 added source with the full catalogue, ~29 with 10% of it. Di Mauro's pipeline sits somewhere in between and resolves at 0.08° pixels rather than our 0.1°. Neither run finds a morphology-specific bias, and the 95% limits on one are |Δq| ≲ 0.005 (full) and ≲ 0.04 (sparse). What the sparse run does establish is that the *starting catalogue* matters far more than the iteration does.

Also: one IEM, one mask, two morphologies, 1–10 GeV, and a ruler whose radial profile is not the injected one.

### 6.7 PSF correction — 2026-09-23

**If you use `utils/psf.py` for shape work, pull the current version.** Everything before commit `fc411c0` used an analytic King fit to the published containment curve, flagged in the file as needing validation. `gtpsf` on our own livetime cube (`scripts/run_gtpsf.sh`) says it was fine below 10 GeV and badly wrong above:

| E (GeV) | r68 gtpsf | r68 analytic | off by |
|---|---|---|---|
| 1.78 | 0.508° | 0.508° | 0.2% |
| 8.72 | 0.167° | 0.154° | 7.9% |
| 10.22 | 0.154° | 0.138° | 10.2% |
| 50.0 | 0.100° | 0.069° | **30.3%** |

`psf.py` now interpolates the gtpsf table directly — the real radial profile, not a King fit — with the analytic form kept only as a fallback. Agreement is 0.0% by construction; the table is committed at `utils/psf_data/gtpsf_P8R3_SOURCE_V3_evtype3.fits`.

**Scale of the effect, so nobody has to guess.** A too-narrow model PSF over-predicts every point source's peak, so the fit pushes the 4FGL normalisation down and the signal templates absorb the difference:

- Mask B GCE TS **4,236 → 1,986**. Every TS was inflated by roughly a factor of two.
- The **0.5–1 GeV detection disappeared entirely** (VVV TS 19 → 0). A whole energy bin's result was an artifact of a PSF error that only exceeded 10% above 10 GeV — the low-energy bin was affected through the mask and the point-source template, not directly.
- Mask A's kept fraction **22% → 8.7%**, because the real 95% containment is wider.
- Point-source normalisation **1.00 → 0.85**. The apparent perfect agreement had been two errors cancelling.

All 16 fits were rerun. This is why the validation step is worth doing before, not after.

### 6.8 Handoff to Gap 7

`outputs/gap1/gap1_results_for_gap7.csv`: 408 rows, one per (mask, point-source treatment, energy bin, template), `result_id` in the Conventions §11 form `gap1_{mask}_{iem}_{template}_{ebin}`. Columns: TS, lnL, ΔAIC vs VVV and vs the best gNFW, convergence, kept sky fraction, source run, git commit. Each cell comes from one run, so every row in a cell shares one counts cube and one exposure. ΔAIC compares within a cell, not across cells. Ten cells, all IEM family 1, all with the real PSF. Regenerate: `python analysis/gap1/export_for_gap7.py`.

**The covariance can be computed on these now, but conclusions should wait for Gap 4 and Gap 5** — with one IEM family the mask-vs-IEM split cannot be separated, and the HESTIA-vs-VVV entries are already known to be configuration-dependent (§6.2).

### 6.9 Still blocked

| what | on whom | why it matters |
|---|---|---|
| **A ring-decomposed diffuse model — now the critical path, not a robustness check** | Jack (download) → then us | §6.0: with one template and one normalisation per bin, 98% of our TS is background error. Nothing in §6.1–6.5 can be re-stated without this. Status of each family below. |
| **Family 2 — Cholis et al. 2022, 80 models. AVAILABLE.** | Jack, by hand | Afeefa was right: Zenodo record **6423495**, one file `GCE_TEMPLATES_ZENODO_FILES_v3.zip`, **2.0 GB**, CC-BY, paper arXiv:2112.09706. `zenodo.org` is blocked by the egress proxy in both the cloud sandbox and the desktop workspace, so it has to be downloaded in a browser. Once it is in the `fermi data` folder the ingestion is ours. |
| ~~Family 3, Pohl22~~ **DONE 2026-09-30** | — | Zenodo **6276721**, the archive the Bubbles template came from. Searched both machines: only `bubbles_fuzzy` and `bubbles_sharp` survived, the rest of the archive is gone. It holds four H I rings (`HI_pohl_Ts_bestTexc_ring_1..4`), four H₂ rings (`H2_pohl_ring_1..4`), two dust residual maps and a boxy bulge. **That is a ring decomposition, and it is the one thing that would let §6.1–6.5 be redone.** Downloaded and regridded; results in §7 of `docs/diffuse_freedom_test.md`. Family 2 is next. |
| **Family 3, Galp21 — NOT publicly available** | the authors | Checked the GALPROP site, its code and WebRun pages, the supplementary-material index, SourceForge, and every Zenodo/Dataverse deposit tied to arXiv:2112.12745. There is no precomputed SA0/SA50/SA100 product anywhere; what is public is GALPROP v57 plus the GALDEF configuration files, i.e. you regenerate the maps yourself. Di Mauro's "distributed with the public GALPROP release" means the configuration, not the maps. **This one needs an email to Porter or Jóhannesson, or a GALPROP run.** |
| **One number**: unresolved point-source flux ÷ GCE flux, in our ROI, mask and energy range | Gap 5 (NPTF) | No longer open-ended. §6.6's Table 10 calibrates the bias against exactly this quantity, so one number closes it. Below ≈0.3 §6.1 is a measurement; above ≈1 it is not quotable. |
| Gap 5 review of `analysis/level2/` | Gap 5 | Offered as a shared module. Worth their review before any run that gets compared between teams. |
| ρ² definition | Muru | Decides whether these are annihilation templates at all (see the framing note under the title). Ask for the unnormalised maps in the same email — §6.10 shows the ones we have carry no absolute J-factor. |
| G3.1 / G3.2: MW or M31 analogues? | Muru | They are the two roundest halos and lose in every cell. If they are the M31 of their pairs they leave the sample and §6.1's spread narrows. **Two concrete things to put in the email.** The halo IDs are not numbered consistently: simulations 09_18 and 17_11 both supply halos …0002 and …0003, but 37_11 supplies …0001 and …0002 — so a rule like "the lower ID is the MW" cannot be right for all three, and 37_11 is exactly the pair in question. And the maps are normalised to a common total (§6.10), so we cannot resolve it ourselves by mass. |
| Independent check of the 1–10 GeV detection window | **runnable now — needs a laptop** | §6.4: TS 0 outside 1–10 GeV rests on our own fitter. Config, driver and pass/fail criteria are written (`analysis/gap1/fermipy_check/`, `scripts/run_fermipy_check.sh`); fermitools and fermipy are conda-forge only and the sandbox's egress proxy blocks conda-forge, so this is a laptop job. |
| FL16Y catalogue | project | Using 4FGL-DR4 meanwhile. |

### 6.10 Caveats

Ordered by size. The first one is larger than everything under it put together and is why §6.1–6.5 are withdrawn.

**1. The background model choice dominates the measurement.** §6.0: the same data give VVV TS 3.8 to 305.6 across five splits of the diffuse model that a control shows are not degenerate with the signal. One all-sky interstellar template with one normalisation per bin cannot decide between them. Full write-up, including every injection control number, in `docs/diffuse_freedom_test.md`.

**2. One IEM family.** The eighty-fold spread above is *within* `gll_iem_v07`. The spread between families is a separate number and probably a larger one. A large-scale north/south residual remains in every variant (Figure 3): over-predicted above the plane, under-predicted below, strongest at l < −10°. That residual is the thing §6.0 shows the signal templates were absorbing.

**3. Templates are close to degenerate.** Inside Mask B at 1.15 GeV, PSF-convolved and exposure-weighted: HESTIA–VVV 0.943, gNFW–HESTIA 0.948, gNFW–VVV 0.918, VVV–F98 0.984; condition number 1.2×10⁴. With 5 M photons they separate statistically, but at 94% overlap a ΔAIC between them is the quantity most exposed to background mismodelling — which is exactly what §6.0 found.

**4. Mask D plus a single-template background has almost no sensitivity.** TS 18.6 rigid, 0.0 once the model is freed across the ROI, and three of its six splits fail the admissibility test outright because they cannot recover a signal injected at the observed flux. The two that survive still span a factor of 46 (18.6 to 79.3). **Gap 1 cannot use Mask D until Pohl22 or Family 2 is loaded.** Gap 5 can, because its IEMs are ring-decomposed by construction; Gaps 4 and 7 inherit the constraint. **Mask A is worse: no split at all is admissible**, including the rigid one, so that configuration cannot measure a TS in either direction. Belongs in the Conventions.

**5. Catalogue incompleteness biases recovered flattening,** by −0.079 ± 0.030 per unit unresolved-flux fraction for a flattened injection and −0.168 ± 0.014 for a spherical one, over f ≤ 1.2 (§6.6, Table 10). Below f ≈ 0.3 the correction is subdominant to the realization scatter. Independently supported by Zhong & Cholis 2024, the only published case where source depth is varied at fixed mask. Documented as a systematic, with the calibration curve, pending Gap 5's NPTF number.

**6. Point-source normalisation is not 1, and its sign depends on the mask:** 0.81–0.88 under Mask B, 1.12–1.84 under Mask E. Flagged, not chased. Not a clean exposure validation, and the mask dependence argues against a single global exposure error — it tracks which sources each mask leaves in. Re-check across IEM families.

**7. Fermitools 2.5.3, not Di Mauro's 2.2.0 — documented as a systematic rather than matched.** 2.2 is not installable on osx-arm64 through the current conda-forge channel, containerising it would cost days, and the parts of the chain we use (gtselect, gtmktime, gtbin, gtltcube, gtexpcube2, gtpsf) are stable across that range. Revisit only if a discrepancy with Di Mauro survives the diffuse-model fix.

**8. Latitude bands are not ring templates.** §6.0's splits give the diffuse model freedom in latitude, with longitude used only to separate inside the GCE from outside it. Our geometry standing in for gas astronomy: enough to show the background choice dominates, not enough to settle which choice is right. Do not use any of them as a background model for a result.

**9. Axis-ratio conventions are reciprocal between papers.** Daylan et al. 2016 define it as the inverse of Di Mauro's and Cholis's. Pooling without inverting Daylan manufactures a trend. Conventions §11.

**10. HESTIA's ρ².** See the framing note under the title.

**11. The maps carry shape only — there is no absolute J-factor in them.** All fifteen of Muru's `dm2` maps sum to 328280.6 over the 80°×80° field, agreeing to 2×10⁻⁷ across every halo, simulation and smoothing scale, while their peaks differ by a factor of 5.5. Nothing here depends on it, since every fit floats an independent normalisation per bin, but amplitudes cannot be compared between halos, Gap 4 and Gap 7 cannot get an absolute annihilation flux from these files, and the MW-vs-M31 question cannot be settled by mass ranking. `analysis/gap1/check_hestia_normalisation.py`.

**12. Level 1 only, 4FGL-DR4 rather than FL16Y, one ROI, one event selection.**

### 6.11 Template shapes before any fit

*(From the first round, unchanged.)* 72 shape measurements (six galaxies × four angles × three contour levels) on Muru's uncropped maps, in `outputs/gap1/shape_measurements.csv`.

- **The flattening reproduces.** Axis ratios of 0.52–0.71 across all six galaxies (Table 1), stable from ~5° to ~20° and across viewing angles, matching Muru et al. §6.2 now shows this is the property the data actually select on.
- **The halos are flattened, not boxy.** Across all 72 measurements c₄ = −0.003 ± 0.007, where a boxy n = 4 superellipse gives +0.019 on the same estimator.
- **The square ROI fakes boxiness.** Measured inside the ROI, all six gave c₄ = +0.011 to +0.021, entirely a crop artifact (Table 2). The radius guard in every shape function now refuses those contours.

### Figures

**Figure 1** — `outputs/gap1/hestia_shape_summary.png`. HESTIA maps with contours; axis ratio vs viewing angle; the ROI crop artifact.

**Figure 2** — `outputs/gap1/fig2_flattening_and_calibration.png`. Axis ratio vs contour radius; c₄ calibration on superellipses.

**Figure 3** — `outputs/gap1/fig3_fit_residuals_{maskA,maskB,maskB_ps}.png`. For each variant: fitted signal map and residual significance for gNFW, HESTIA, VVV.

**Figure 4** — `outputs/gap1/fig4_hestia_envelope.png`. ΔAIC vs VVV for all 24 HESTIA configurations in all three variants, with the best-γ gNFW marked.

**Figure 5** — `outputs/gap1/fig5_envelope_by_energy.png`. The same, per Conventions energy bin (1–3 GeV, 3–10 GeV) and for 0.5–50 GeV, Mask B + PS.

**Figure 6** — `outputs/gap1/fig6_maskD.png`. Mask D with sources masked, Mask D with them modelled, and Mask B with them modelled. The gNFW loses under all three; HESTIA vs VVV crosses from one side to the other.

**Figure 8** — `outputs/gap1/fig8_diffuse_freedom.png`. TS of the VVV template against how the diffuse model is allowed to bend, Mask B + PS, with both controls. The observed excess runs from 306 to 3.8; a signal we injected and a signal added to the real sky both survive throughout, so the range is a property of the background choice rather than of the test.

**Figure 7** — `outputs/gap1/fig7_incompleteness.png`. Left: bias in recovered axis ratio against the flux in sources the fit was not given, for a flattened and a spherical injection, with the f ≤ 1.2 calibration and the saturating fit. Right: the same runs as recovered q, showing two truths 0.25 apart converging on q ≈ 0.59.

### Tables

**Table 1.** Shape of the HESTIA J-factor maps at the 5%-of-peak contour (~7–13°), averaged over the four viewing angles. Errors are the spread across angles. A pure ellipse gives c₄ = 0.000, a boxy n = 4 superellipse +0.019.

| Galaxy | Simulation / halo | axis ratio q | boxiness c₄ |
|---|---|---|---|
| G1.1 | 09_18 / …003 | 0.548 ± 0.017 | -0.009 ± 0.008 |
| G2.1 | 17_11 / …003 | 0.590 ± 0.004 | +0.004 ± 0.004 |
| G3.1 | 37_11 / …002 | 0.693 ± 0.013 | -0.002 ± 0.003 |
| G1.2 | 09_18 / …002 | 0.579 ± 0.017 | -0.001 ± 0.003 |
| G2.2 | 17_11 / …002 | 0.610 ± 0.021 | -0.010 ± 0.003 |
| G3.2 | 37_11 / …001 | 0.633 ± 0.016 | -0.013 ± 0.004 |

**Table 2.** The ROI artifact: c₄ of the 2%-of-peak contour inside the ROI vs on Muru's uncropped ±40° maps.

| Galaxy | c₄ inside the ROI | c₄ uncropped |
|---|---|---|
| G1.1 | **+0.020** | -0.003 |
| G2.1 | **+0.013** | -0.004 |
| G3.1 | **+0.021** | -0.003 |
| G1.2 | **+0.015** | -0.004 |
| G2.2 | **+0.011** | -0.006 |
| G3.2 | **+0.015** | +0.002 |

**Tables 3–5** are in §6.2, **6** in §6.4, **7** in §6.5, **8–9** in §6.6. Per-configuration numbers for all 40 templates in every cell: `outputs/gap1/{scan_mask*,eb_scan_*,tri_scan_*}/scan.csv`, or the single combined table in §6.7. Full fit records: `outputs/gap1/fit_*/fit_results.json`.

### Reproduce

```bash
# laptop: scripts/run_chain.sh   -> gap1_ccube / expcube / iem_roi
python analysis/gap1/fit_cubes.py --ccube gap1_ccube.fits --expcube gap1_expcube.fits \
  --iem gap1_iem_roi.fits --catalog gll_psc_v35.fit --mask B --bubbles --ps \
  --signals gnfw,hestia,vvv --hestia-wide <Muru G1.1 0deg 1deg map> --out outputs/gap1/fit_maskB_ps
python analysis/gap1/robustness_scan.py ... --mask B --ps --out outputs/gap1/scan_maskB_ps
python analysis/gap1/make_fit_figures.py --fitdir outputs/gap1/fit_maskB_ps
python analysis/gap1/make_envelope_figure.py masks
# conventions energy bins: laptop scripts/run_ebins.sh, then per bin
python analysis/gap1/robustness_scan.py ... --mask B --ps --emin 1 --emax 3 --out outputs/gap1/eb_scan_b1-3
python analysis/gap1/build_triaxial.py      # Hu et al. readings
python analysis/gap1/robustness_scan.py ... --mask B --ps --library triaxial --emin 1 --emax 10 --out outputs/gap1/tri_scan_b1-10
python analysis/gap1/make_envelope_figure.py energy
python analysis/gap1/export_for_gap7.py
# mask D (flux-scaled), and the level-2 engine's self-test
python analysis/gap1/fit_cubes.py ... --mask D --bubbles --ps --emin 1 --emax 10 --out outputs/gap1/fit_maskD_ps
python analysis/level2/iterate.py
# catalogue-incompleteness calibration (sec 6.6, Table 10, Figure 7)
python analysis/gap1/incompleteness_curve.py --ccube gap1_ccube_full.fits \
  --expcube gap1_expcube_full.fits --iem gap1_iem_roi.fits \
  --catalog gll_psc_catalog.fit --fit outputs/gap1/eb_fit_all/fit_results.json \
  --realizations 5 --out outputs/gap1/incompleteness
python analysis/gap1/make_incompleteness_figure.py \
  --calib outputs/gap1/incompleteness/calibration.json --out outputs/gap1/fig7_incompleteness.png
# sec 6.0: how much of the detection is the diffuse model's error
python analysis/diffuse_freedom/run_test.py --ccube gap1_ccube_full.fits \
  --expcube gap1_expcube_full.fits --iem gap1_iem_roi.fits \
  --catalog gll_psc_catalog.fit --mask B --ps --out outputs/gap1/diffuse_freedom_maskB_ps
python analysis/gap1/measure_q.py --ccube ... --mask B --ps --out outputs/gap1/measure_q_maskB_ps
python analysis/gap1/make_diffuse_figure.py \
  --json outputs/gap1/diffuse_freedom_maskB_ps/diffuse_freedom.json \
  --labels "Mask B + PS, 1-10 GeV, VVV template" --out outputs/gap1/fig8_diffuse_freedom.png
# laptop only: independent fitter check of the energy window (sec 6.4)
bash scripts/run_fermipy_check.sh
```

## 7. References

**Flattened / boxy dark matter and simulations**
- Muru, Silk, Libeskind, Gottlöber & Hoffman (2025), *Fermi-LAT Galactic Center Excess Morphology of Dark Matter in Simulations of the Milky Way Galaxy*, Phys. Rev. Lett. — arXiv:2508.06314
- Hu, Cholis & Zhong (2026), *Generic triaxial dark matter halo* — arXiv:2602.20252

**Stellar bulge and morphology**
- Coleman et al. (2019), *Galactic bulge and nuclear bulge templates* — arXiv:1911.04714
- Macias et al. (2016), *Galactic Bulge Preferred Over Dark Matter for the Galactic Centre Gamma-Ray Excess* — arXiv:1611.06644
- Bartels et al. (2018), *The Fermi-LAT GeV excess as a tracer of stellar mass in the Galactic bulge* — arXiv:1711.04756
- Song et al. (2024), *Robust inference of the Galactic Centre gamma-ray excess spatial properties* — arXiv:2402.05449
- Ramirez et al. (2024), *Inferring the morphology of the GCE with Gaussian processes* — arXiv:2410.21367

**Pipeline and high-precision analyses**
- Di Mauro (2026), *A Precise Measurement of the Fermi-LAT Galactic Center Excess Morphology and Spectrum* — arXiv:2605.22913
- Di Mauro (2021), *The characteristics of the GCE measured with 11 years of Fermi-LAT data* — arXiv:2101.04694
- Zhong & Cholis (2024), *Robustness Against Masking* — arXiv:2401.02481

**Systematics and mismodelling**
- Calore, Cholis & Weniger (2015), *Background model systematics for the Fermi GeV excess* — arXiv:1502.02805
- Leane & Slatyer (2020), *Spurious point source signals in the GCE* — arXiv:2002.12370
- Leane & Slatyer (2020), *The enigmatic Galactic Center excess* — arXiv:2002.12371
- Leane & Slatyer (2019), *Revival of the dark matter hypothesis* — arXiv:1904.08436

**Context**
- Hooper & Goodenough (2011) — arXiv:1010.2752
- Daylan et al. (2016), *Characterization of the central Milky Way* — arXiv:1402.6703
- Fermi-LAT Collaboration (2017) — arXiv:1704.03910
