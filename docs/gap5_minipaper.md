# Gap 5: Phantom Menace in the Iterative Pipeline

*20 Oct 2026. Builds on Hussain's draft and IEM selection. Complete for Mask B, 1–10 GeV: 2 reference backgrounds and 23 Cholis models, with ordering controls. Results in §6.4–6.5; changes and amendments at the end.*

---

## 1. Specific Research Question

Gap 1 found that, inside the Level 2 iterative pipeline, the VVV stellar bulge beats both spherical gNFW and HESTIA flattened DM. That result used a ring-decomposed background (Pohl22). Di Mauro (2026) runs a similar iterative pipeline, with Galp21 SA100 as his best-fitting background, and finds a spherical, DM-like excess. McDermott, Zhong & Cholis (2023) found that GALPROP-template backgrounds fit the inner Galaxy far better than ring backgrounds (Δχ² ≈ 6500), and that under them the best-fitting excess has a DM-like morphology.

**Question:** When the ring background is replaced by the extreme GALPROP models of Cholis et al. (2022), run through the same Level 2 pipeline, does the morphology verdict change? And does the excess survive at all?

This is the last of the four candidate explanations for the Macias vs Di Mauro contradiction. Gap 1 has ruled out the mask, the template quantity (squared column vs J-factor), and the pipeline. The background model is the one still open.

## 2. Null Hypothesis

**H₀:** The GCE morphology and significance do not depend on the extreme IEM variations.

**How H₀ is tested.** For each Family 2 model that passes both checks below, Level 2 gives the same answer as Gap 1: VVV is preferred over gNFW and HESTIA (ΔAIC(gNFW − VVV) > 0), and the best template stays detected at TS ≥ 25.

**H₀ is rejected** if any model that passes both checks reverses the ordering or loses the excess.

The two checks, both fixed before any run:

- **Admissible.** This is the Gap 1 criterion, unchanged. An injected signal at the observed flux must be recovered at ≥ 50% of its rigid-split TS, and at TS ≥ 25. A model that fails is reported as "cannot measure", not as a null.
- **Competitive.** A model only counts against H₀ if the data don't already reject it more strongly than it moves the verdict. The rule: count a flip only when the model's background-fit penalty against the best Family 2 model, in full AIC, is smaller than the |ΔAIC| of the flip it produces. Without this rule, a model the data reject by thousands could "kill" the excess. That would say nothing about the GCE.

All models are reported either way. The rule decides what counts as evidence; it does not hide rows.

## 3. Data and Methods

**Data.** Project conventions: 16.5 yr, P8R3_SOURCE_V3, evtype 3, |l|, |b| < 20°, 0.1° CAR.

**Templates.** These are the same four Gap 1 Level 2 used: VVV (Coleman 2019), F98 boxy bulge, gNFW γ = 1.2, and HESTIA G3.2 at 90°. Keeping the set identical is what makes Gap 5's answer comparable to Gap 1's.

**Masks.** Mask **B** is primary and Mask **D** is secondary. Mask D is the TS/flux-scaled mask from Table III of the same Cholis et al. (2022) paper the IEMs come from, so it is the natural pairing.
*Masks A and C are dropped.* Main paper §2.2 shows Mask A has no admissible split, not even the rigid one, under any background. Mask C equals Mask A at GeV energies and is marked diagnostic only.

**IEMs.** There are twelve Family 2 models from Cholis et al. (2022), selected as described in §5.1. Two reference backgrounds go through the identical code path: gll_iem_v07 (Family 1) and Gap 1's Pohl22 rings_species.

**Background freedom.** Gap 1 showed that background freedom moves the morphology verdict by a factor of 30. Each Cholis model is therefore fitted as separate GALPROP components (π⁰, bremsstrahlung, inverse Compton), each with a free normalisation, plus isotropic and Bubbles. If the release only provides total maps, we fit them with one normalisation and flag the freedom mismatch against Pohl22. The parameter count k is reported for every row, and models with different k are compared on full AIC = 2k − 2lnL.

**Energy.**
- Primary: **1–10 GeV.** This matches Di Mauro's morphology range and Gap 1's Level 2, so all three results sit on the same grid.
- Secondary: **0.5–3 GeV and 3–50 GeV.** Masking effects are energy-dependent, and Cholis et al. find the gas assumptions matter most at low energies, which is exactly where these IEMs differ most.

**Statistics and software.** Binned Poisson likelihood. Level 2 iterative source finding runs through `analysis/level2/run_level2.py` (the Gap 1 module, with a new `--components` option), in both arms:
- **Neutral** (run first): one common source list, then all 24 HESTIA maps, 7 gNFW slopes, VVV and F98 fitted at equal parameter count.
- **Joint:** each hypothesis gets its own source search, as Di Mauro does. The list is Gap 1's six templates **plus the best HESTIA and best gNFW from that background's own neutral arm** (amendment 1, below).

Fermipy is used as a cross-check on one model only. Its gammapy/regions dependency clash is still unresolved, and the Gap 1 fitter is the validated path.

**Per-model quantities.** These are the "Phantom Menace in the pipeline" measurement itself:
- the background-only fit quality (lnL, k, AIC)
- the TS of each template, and ΔAIC versus VVV
- the number of point sources Level 2 adds
- whether any added source lands within 0.5° of the GC, and how much excess flux it absorbs

The worry this tests: a worse background leaves larger residuals, so Level 2 adds more sources, which can absorb more of the excess. Gap 1 already saw one added source at (l, b) = (−0.35°, +0.35°) absorb part of the excess under Pohl22.

**Controls** (the Gap 1 standard):
1. Admissibility injection, run per model × mask.
2. Direct-injection ordering on 4 models spanning the range. Inject gNFW and require gNFW to win; inject VVV and require VVV to win. This is the only control that licenses a "morphology X is preferred" claim.
3. A null sky per control model. With no excess injected, no template may be detected, and no sources may be invented at the GC.

**Pre-registration.** All rules are in `configs/gap5.yaml`. Two amendments were made on 19 Oct 2026, both after running only the two reference backgrounds and before any Family 2 model:
1. Gap 1's joint list holds the two halos that were best *under Pohl22*. Under gll_iem_v07 a different halo (G2.1) is best, and the fixed list never tested it. The joint arm now also gets the best HESTIA and best gNFW of each background's own neutral arm.
2. The q ladder was extended from 0.30 down to 0.10, because under Pohl22 the q profile was still rising at 0.30.

**Run order** (first answer after the first 12 runs):
1. Mask B, 1–10 GeV, joint arm, all 12 models
2. Neutral arm
3. Mask D
4. Energy split
5. Controls

## 4. Expected Results and Interpretation

| Outcome | What it means |
|---|---|
| VVV first under every admissible, competitive model | Gap 1's bulge result is not a ring-background artifact. Phantom Menace's DM-like preference does not survive Level 2. The IEM is eliminated as the source of the contradiction, which closes all four explanations. |
| gNFW first under some admissible, competitive models | The morphology verdict depends on the background family. The Macias vs Di Mauro contradiction is, at least in part, an IEM contradiction. |
| Excess lost (TS < 25) under some competitive model | The excess itself is background-dependent inside an iterative pipeline. This is the strongest Phantom-Menace-type result. |
| Most models inadmissible | Family 2 cannot measure morphology in this pipeline. Reported as such, not as a null. |

**One cross-cutting result.** Comparing background-fit quality between Family 2 and Pohl22 in our pipeline tests McDermott et al.'s claim that templates beat rings. If the Family 2 models fit much better, their verdict should carry more weight than Gap 1's.

## 5. Potential Pitfalls and Mitigations

### 5.1 The IEM selection

There are four problems with the current 12. Points 1–3 are checked by `analysis/gap5/iem_selection_audit.py` (in the repo) and the paper's Table II.

1. **The library's extremes were dropped.** The six Table II models seeded the selection but are not in the final 12. Those six hold the library's extremes.
   - The library spans vA = 0–50 km/s and dv_c/dz = 0–200 km/s/kpc (Cholis et al. 2022). The 12 stop at **32.2** and **100**.
   - The 50 km/s models are II and XV; the 200 km/s/kpc model is VI.
   - Several picks are diluted copies of those seeds. **LI and LII share VI's halo height, diffusion, gas, ISRF and B-field, but have convection 60 and 100 instead of 200.**
2. **Near-duplicates.** LI and LII differ in exactly one parameter, convection. The 12 fall into only **4 distinct ISRF/gas/B-field configurations**, and **6 of the 12** share one of them. Nearest-neighbour distances vary 6× across the set (0.26 to 1.55, normalised).
   This is the failure the draft's own pitfall note warned about. It comes from the objective: "largest *sum* of distances" rewards being far from the set's average, and two outliers next to each other both score highly. "Largest *minimum* distance" (farthest-point selection) is the objective that guarantees spread.
3. **Confounded electron profile.** Nuclei use the pulsar profile in all 12. The Pul/SNR electron split (6/6) lines up exactly with the configuration split: all six Pul/Pul models are the ISRF 0.7 / gas 9,4,1 family. Any difference between the halves cannot be pinned on the electron profile.
4. **Normalisation and encoding.** "Divide by the largest value" maps a column to [min/max, 1], not [0, 1], so narrow-range parameters like δ (0.33–0.5) are under-weighted. The B-field entry is a packed code (B₀, r_c, z_c), and gas and ISRF are category labels. None of these can enter a Euclidean distance as a single raw number.

**Mitigation (recommended).** Re-select 12 as follows:
- Keep the **6 Table II models (I, II, III, VI, X, XV)** as fixed anchors. They hold the extremes, two of the paper's five best fits (X, XV), its reference model (I) and one worst fit (II).
- Add **6 more by farthest-point selection** over all 80 models.
- Use min-max scaling on the continuous columns, unpack the B-field into B₀, r_c and z_c, and score the categorical columns as match/mismatch.
- Require every gas configuration in the library to appear at least once.

This keeps the count at 12, as the conventions specify. If the full Table VIII is not available in time, the fallback is to run the current 12 plus the 6 anchors (18 models).

The anchors also guard against the "competitive" rule leaving no test at all. At present only one of the 12 (LIII) is among the paper's five best fits. If most of the set fits poorly, the rule could exclude nearly everything.

### 5.2 Other pitfalls

- **Format and energy grid.** GALPROP outputs may be HEALPix, at different resolution, with a different energy grid and units. Mitigation: regrid with `utils/regrid.py` and interpolate log-log onto our bins. Gap 1's inverse-Compton ring bug (4 energy planes against 8 bins) was exactly this failure, so every model is checked for plane count and units before fitting.
- **Freedom mismatch.** This is handled by the component split and full AIC (§3). Without it, a Family 2 vs Pohl22 difference could be a freedom difference rather than a physics difference.
- **One flip read as a verdict.** The result is the pattern across all 12. Every row is reported, with its admissible and competitive flags.
- **Level 2 is a reimplementation.** One configuration is matched to Di Mauro's selection (SOURCEVETO, 0.08°, 1–10 GeV) as a validation anchor, as in Gap 1.

## 6. Preliminary Findings

### 6.1 Current IEM selection

The paper labels models with Roman numerals. The numbers below map as 5 = V, 7 = VII, … 56 = LVI.

| # | zL | D0 | δ | vA | dv_c/dz | αp1, αp2 | αe1, αe2 | Np, Ne | ISRF | H₂, HI, HII | S_N/S_e | B-field |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V | 10 | 10 | 0.33 | 32.2 | 0 | 1.7, 2.39 | 1.6, 2.33 | 1, 2.66 | 1.4, 1.4 | 0, 5, 2 | Pul/SNR | 200040050 |
| VII | 10 | 8 | 0.33 | 0 | 0 | 1.4, 1.8 | 1.4, 2.3 | 1.3, 3.33 | 1.4, 1.4 | 0, 5, 2 | Pul/SNR | 200040050 |
| XVII | 5.6 | 4.85 | 0.4 | 24 | 1 | 2, 2.25 | 1.6, 2.3 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| XIX | 6 | 6.5 | 0.33 | 30 | 0 | 2.04, 2.41 | 1.6, 2.43 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| XXII | 5.5 | 5.5 | 0.37 | 30 | 2.5 | 2, 2.38 | 1.6, 2.43 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| XXV | 5.7 | 3.9 | 0.45 | 25.7 | 6 | 1.99, 2.36 | 1.6, 2.43 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| XXVII | 6 | 3.1 | 0.5 | 23 | 9 | 2.02, 2.38 | 1.6, 2.43 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| XXX | 3 | 2.67 | 0.4 | 22 | 3 | 2.08, 2.41 | 1.6, 2.43 | 5.8, 2 | 0.7, 0.7 | 9, 4, 1 | Pul/Pul | 025200010 |
| LI | 6 | 2 | 0.33 | 0 | 60 | 1.6, 2.3 | 1.6, 2.26 | 1.5, 5.9 | 1.4, 1.4 | 9, 5, 1 | Pul/SNR | 200030050 |
| LII | 6 | 2 | 0.33 | 0 | 100 | 1.6, 2.3 | 1.6, 2.26 | 1.5, 5.9 | 1.4, 1.4 | 9, 5, 1 | Pul/SNR | 200030050 |
| LIII | 10 | 8 | 0.33 | 32.2 | 100 | 1.4, 1.8 | 1.4, 2.3 | 1.3, 3.33 | 1.4, 1.4 | 0, 5, 2 | Pul/SNR | 200040050 |
| LVI | 6 | 40 | 0.33 | 0 | 0 | 1.89, 2.1 | 1.4, 2.1 | 2, 1.1 | 1, 1 | 0, 5, 2 | Pul/SNR | 050100020 |

Units: zL in kpc; D0 in 10²⁸ cm²/s; vA in km/s; dv_c/dz in km/s/kpc. The B-field code packs B₀ × 10 (µG), r_c × 10 (kpc) and z_c × 10 (kpc).

### 6.2 What the 12 cover

This replaces the draft's verification paragraph, which overstated the coverage.

**Covered:**
- The 12 span the library's halo height (zL = 3–10 kpc) and diffusion normalisation (D0 = 2–40 × 10²⁸ cm²/s), including LVI, the only D0 = 40 model.
- They include the highest electron normalisation (Ne = 5.9).
- They include LIII, one of the paper's five best-fitting models.

**Not covered:**
- **Reacceleration:** the 12 stop at 32.2 km/s against the library's 50.
- **Convection:** the 12 stop at 100 km/s/kpc against the library's 200.

**Structure of the set:**
- The electron source profile splits 6/6 (pulsar/SNR), but the split coincides exactly with the ISRF/gas/B-field configuration split, and nuclei are pulsar-like in all 12.
- The set contains 4 distinct configurations, with 6 models in one.

See §5.1 for the fix.

### 6.3 Results so far: pipeline validation and the two reference backgrounds

**Validation.** Gap 1's Pohl22 background, put through the new component path, gives back Gap 1's Level 2 joint arm exactly: all six templates have the same lnL to four decimal places, the same number of added sources at the same positions, and the same parameter count. At Level 1 the two code paths agree to 2×10⁻⁵ in lnL. gll_iem_v07 admissibility through the Gap 5 driver reproduces Gap 1's ring table (VVV 184.3 / +585, gNFW 109.3 / +320).

**Mask B + PS, 1–10 GeV, Level 2.** ΔAIC is against VVV; negative means better than VVV.

| background | free diffuse comps | admissible | joint winner (margin to next class) | best HESTIA, neutral | best gNFW, neutral | sources added | at GC |
|---|---|---|---|---|---|---|---|
| Pohl22 rings_species (Gap 1) | 5 | yes | **VVV** (+439.9) | G3.2 90°: +450.3 | γ 1.4: +708.5 | 10–13 | 1 |
| gll_iem_v07 | 1 | yes | **HESTIA G2.1 0°** (+41.9) | G2.1 0°: **−48.6** | γ 1.1: +85.0 | 4–5 | 0 |

Under gll_iem_v07, 7 of the 24 HESTIA maps beat VVV by more than 10. In the joint arm, HESTIA G2.1 wins with the same 4 sources and the same parameter count as VVV, so the win isn't bought with extra point sources.

**Fit quality.** Same data, mask and bins: the best fit under gll_iem_v07 has an AIC **24,358 lower** than the best fit under Pohl22, with 72 free normalisations against 152.

**Caveat.** gll_iem_v07 is the Fermi team's point-source model. Its documented predecessor (Acero et al. 2016) was fit to LAT data and carries a data-derived "extended excess" component, so part of that advantage is built in. It is not independent evidence that gll is the better physical background.

**Mask D.** Neither background passes admissibility. An injected gNFW at the reference flux comes back at TS 10.6 under Pohl22 and 6.6 under gll, both below 25. Mask D therefore cannot measure the spherical template under either background, which matches Gap 1's finding that admissibility is per template. Under Pohl22, VVV still leads in both arms: by 379 over the best HESTIA in the neutral arm, and by 329 over HESTIA and 893 over gNFW in the joint arm.

**Other energy ranges, Mask B.**

| range | background | admissible | neutral: best HESTIA vs VVV | neutral: best gNFW vs VVV | joint winner (margin) | sources added |
|---|---|---|---|---|---|---|
| 0.5–3 GeV | gll_iem_v07 | yes | −46.6 | +93.2 | F98 (+2.4, not decisive) | 3 |
| 0.5–3 GeV | Pohl22 | yes | *running* | *running* | *running* | **68** (neutral arm) |
| 3–50 GeV | gll_iem_v07 | no (gNFW TS 16.6) | −23.8 | +25.8 | HESTIA G2.1 (+23.8) | 5 |
| 3–50 GeV | Pohl22 | no (gNFW TS 10.0) | *running* | *running* | *running* | *running* |

At 3–50 GeV neither background can see an injected gNFW, so that range cannot test the spherical template.

At 0.5–3 GeV, Pohl22's neutral arm adds 68 point sources, against 3 under gll. This is the "phantom menace in the pipeline" mechanism itself. Where the background fits worst (Cholis et al. find gas assumptions matter most at low energy), Level 2 patches the gap with dozens of point sources.

**Shape of the excess (q ladder, Level 1, Mask B).** q is the vertical / in-plane axis ratio; Di Mauro's r = 1.10 ± 0.05 corresponds to q = 0.909 ± 0.041.

| background | 0.5–3 GeV | 1–10 GeV | 3–50 GeV |
|---|---|---|---|
| gll_iem_v07 | 0.369 | 0.382 [0.369, 0.425] | 0.322 |
| Pohl22 | 0.272 | 0.285 [0.280, 0.297] | 0.348 |

- **Every measurement is strongly flattened, 9–16σ from Di Mauro's q.** That is consistent with spherical gNFW losing everywhere.
- **Caveat.** q is lowest under the worse-fitting background. That is what you would expect if mismodelled Galactic-plane emission is part of what a flattened template absorbs, so q should not be read as a clean shape measurement of the excess alone.
- Where the profile is sharper than one ladder step, the interval comes from its curvature. It agrees with the crossing method wherever both exist.

**What this says so far.**
1. **Spherical gNFW loses under every background, mask, energy range and arm run so far.** That is the robust part of the Gap 1 result.
2. **Bulge vs flattened DM depends on the background.** It is VVV by 440 under Pohl22 and HESTIA by 42–49 under gll_iem_v07. Gap 1's "VVV wins at Level 2" is a Pohl22 result, and Pohl22 is the worse-fitting of the two backgrounds by 24,000 in AIC.
3. The Cholis models therefore decide two things: where GALPROP-type backgrounds put the verdict, and how well they fit compared with both references. McDermott et al. (2023) found GALPROP templates fit far better than rings, so this is now the deciding test for the paper's headline.

*Still running at the time of writing: Pohl22 at 0.5–3 and 3–50 GeV. Family 2 is pending the Zenodo download.*

### 6.4 Family 2 results: 23 Cholis models, Mask B, 1–10 GeV

**Fit quality.** Every Cholis model fits the data better than Gap 1's Pohl22, by 3,900 to 16,500 in AIC (same data, mask and bins). The standard Fermi background, gll_iem_v07, fits best of all. The best Cholis model is LIII, one of Cholis et al.'s own five best fits, and it sits 7,867 behind gll_iem_v07. This is McDermott et al.'s templates-beat-rings result, reproduced inside our pipeline.

**Results, ordered by fit quality.** Fit gap is AIC above the best background. Margins are ΔAIC to the next class.

| background | set | fit gap | admissible | neutral winner (margin) | joint winner (margin) | Level 2 sources |
|---|---|---|---|---|---|---|
| gll_iem_v07 | ref | 0 | yes | flattened (48.6) | flattened (41.9) | 5 |
| LIII | original | 7,867 | yes | spherical (58.7) | spherical (38.7) | 19 |
| X | re-selected | 7,947 | yes | spherical (49.1) | flattened (24.7) | 19 |
| V | original | 7,958 | yes | spherical (45.5) | spherical (15.9) | 19 |
| VII | original | 8,002 | yes | spherical (44.6) | spherical (0.4) | 19 |
| XV | re-selected | 8,307 | yes | bulge (98.8) | bulge (22.7) | 20 |
| LVI | original | 8,410 | yes | bulge (127.3) | bulge (256.3) | 19 |
| LI | original | 9,748 | yes | spherical (162.8) | spherical (69.3) | 17 |
| LII | original | 9,751 | yes | spherical (166.1) | flattened (93.5) | 17 |
| VI | re-selected | 9,803 | yes | spherical (173.6) | spherical (36.3) | 17 |
| I | re-selected | 10,006 | yes | spherical (135.6) | spherical (110.3) | 18 |
| XXXIX | re-selected | 10,278 | yes | spherical (147.0) | spherical (378.0) | 20 |
| XLV | re-selected | 10,592 | yes | spherical (31.5) | spherical (31.6) | 15 |
| IX | re-selected | 10,615 | yes | flattened (20.1) | flattened (20.3) | 37 |
| LXVII | re-selected | 11,065 | yes | bulge (143.4) | bulge (63.5) | 25 |
| III | re-selected | 11,439 | yes | spherical (103.2) | flattened (272.8) | 17 |
| XXX | original | 12,243 | yes | spherical (53.5) | spherical (117.1) | 17 |
| XXVI | re-selected | 12,243 | no | spherical (88.7) | spherical (59.1) | 34 |
| XVII | original | 12,420 | yes | spherical (74.4) | spherical (69.6) | 23 |
| XXVII | re-selected+original | 12,435 | yes | spherical (53.5) | spherical (93.1) | 21 |
| XXII | original | 12,452 | yes | spherical (71.9) | spherical (69.1) | 23 |
| XXV | original | 12,625 | yes | spherical (63.5) | flattened (17.7) | 21 |
| XIX | original | 12,793 | yes | spherical (57.6) | bulge (224.8) | 24 |
| II | re-selected | 20,454 | yes | flattened (119.9) | flattened (182.6) | 14 |
| pohl22_rings_species | ref | 24,358 | yes | bulge (450.3) | bulge (439.9) | 12 |

**Counts.**

| arm | spherical | flattened | bulge |
|---|---|---|---|
| neutral (equal sources, equal k) | 18 | 2 | 3 |
| joint (own source search, as Di Mauro) | 13 | 6 | 4 |

**The arms disagree in 5 models** (III, LII, X, XIX, XXV). In every one, the joint-arm winner is the template whose own source search found more point sources than the gNFW run under the same background:

| model | joint winner | its sources | gNFW sources |
|---|---|---|---|
| III | HESTIA | 19 | 17 |
| LII | HESTIA | 18 | 16 |
| X | HESTIA | 19 | 18 |
| XIX | F98 | 24 | 22 |
| XXV | HESTIA | 21 | 19 |

This is the "phantom menace in the pipeline" mechanism itself: the morphology verdict follows whichever hypothesis gets to convert more residual into point sources. LI and LII show it most cleanly. They differ only in convection (60 vs 100 km/s/kpc) and agree in the neutral arm (spherical by 163–166), yet their joint verdicts are spherical and flattened.

**The best-fitting backgrounds.** The six best-fitting Cholis models are within 543 AIC of each other.
- Neutral arm: spherical in 4 of them (LIII, X, V, VII), bulge in 2 (XV, LVI).
- Joint arm: spherical 2 plus a tie (VII, 0.4), flattened 1 (X), bulge 2.
- Joint-arm margins among them are 0.4–38.7, except LVI (bulge by 256). Neutral-arm margins are 45–127.

**Ordering controls** (pre-registered: V, XVII, LI, LVI; each template injected at the flux it fits on the real sky):
- Injected gNFW is recovered as spherical in 3 of 4 (by 16.5–21.5). Under XVII it comes back flattened, with spherical third, 2.9 behind.
- Injected VVV is recovered as bulge in 4 of 4, though by only +2.2 under LI.
- With nothing injected, no model invents a source and no template wins by more than 2.3.

So the pipeline can tell the classes apart, but only by margins of tens. That is the same size as most real-sky margins under the best-fitting backgrounds. One limitation: the simulated skies contain no unresolved sources, so Level 2 adds none and the controls do not exercise the source-absorption mechanism above.

**The excess survives.** The winning template is detected under every background. Under two models (XXVI, XLV) the excess vanishes at Level 1 (TS 0–1.5); the background absorbs it. It reappears at Level 2 (best gNFW TS 307 and 74) once the iteration has added 34 and 15 sources. XXVI fails admissibility on the boost test and is reported as "cannot measure". XLV passes.

**H₀, by the pre-registered rule.**

| set | H₀ | rejecting models |
|---|---|---|
| all 23 | rejected | LIII, V |
| original 12 | rejected | LIII, V |
| re-selected 12 | holds (11 of 12 measured; XXVI cannot measure) | none |

The rule rejects H₀ through LIII (spherical by 38.7, penalty 0) and V (spherical by 15.9, penalty 91, below its swing of 96 against the bulge). XIX also trips the "excess lost" clause, but only on a technicality: its joint winner (F98) has neutral TS 16, while its best gNFW has TS 99. It is not needed for the verdict.

The re-selected set's own non-bulge verdicts all come from backgrounds that fit worse than LIII by more than the size of their flip, so the rule does not count them. The verdict therefore depends on whether LIII and V are in the set, and that should be stated whenever H₀ is quoted.

**q is unusable under Family 2.** The profile pins to the thinnest rung: under X, q = 0.10 is preferred over 0.15 by 451 in lnL. The flattened ladder is absorbing Galactic-plane residuals, not measuring the excess.

**Not run.** Family 2 was not run on Mask D or the secondary energy ranges. Both reference backgrounds already fail gNFW admissibility there (Mask D: TS 10.6 and 6.6; 3–50 GeV: 16.6 and 10.0), so those configurations cannot test the spherical template.

### 6.5 What Gap 5 says

1. **Gap 1's bulge result does not survive.** VVV wins by 440 only under Pohl22, which fits the data worst of all 25 backgrounds run here. Under the Cholis backgrounds the bulge wins 3 of 23 (neutral) or 4 of 23 (joint).
2. **The excess itself is robust.** It is detected under every background, including the two where Level 1 alone loses it.
3. **Spherical is the most common verdict, not a settled one.** It wins 18 of 23 at equal point sources and 13 of 23 in Di Mauro's procedure. But among the best-fitting backgrounds all three classes appear, and the margins are as small as what the injection controls can resolve. The standard Fermi background, which fits best overall, prefers flattened HESTIA.
4. **The iterative pipeline is a systematic in its own right.** Letting each hypothesis run its own source search changes the class in 5 of 23 backgrounds, every time toward the template that absorbed more residual into point sources.

**Answer to the Gap 5 question.** Di Mauro's spherical GCE survives most extreme IEMs inside the iterative pipeline, but not all of them. Spherical against flattened against bulge is decided by the choice of background and by how the source search is run, not by the data alone. The Macias vs Di Mauro contradiction is, in large part, an IEM contradiction.

## 7. References

- Cholis, Zhong, McDermott & Surdutovich (2022), *The Return of the Templates: Revisiting the Galactic Center Excess with Multi-Messenger Observations*, arXiv:2112.09706. Model maps: Zenodo record 6423495.
- McDermott, Zhong & Cholis (2023), *A Phantom Menace: On the Morphology of the Galactic Center Excess*, arXiv:2209.00006.
- Di Mauro (2026). Iterative source-finding, Galp21/Pohl22.
- Pohl et al. (2022). Ring-decomposed IEM.
- Acero et al. (2016), *Development of the Model of Galactic Interstellar Emission for Standard Point-source Analysis of Fermi LAT Data*, ApJS 223, 26 (arXiv:1602.07246).
- Gap 1 mini-paper and `analysis/gap1/README.md` §6.0c (Level 2 result).

---

## What changed from the draft, and why

| Draft | Revised | Reason |
|---|---|---|
| Masks A, B, C | Masks B, D | A has no admissible split; C = A at GeV; D is Cholis's own mask |
| "Does the spherical GCE survive" | "Does the verdict change" | Gap 1 already found VVV wins at Level 2 under Pohl22, so the open question is whether that depends on the background |
| H₀ in words | H₀ with an operational test | Admissible + competitive rules fixed before running |
| 12 models from max-sum selection | 6 anchors + 6 farthest-point | Restores the dropped extremes; removes near-duplicates |
| Verification paragraph | Corrected | vA and convection ranges were not fully covered |
| 0.5–50 GeV only | 1–10 primary, plus 0.5–3 and 3–50 | Comparable to Gap 1 and Di Mauro; low energy is where the IEMs differ |
| Fermipy | Gap 1 Level 2 module; Fermipy as cross-check | Fermipy dependency clash unresolved |
| — | Per-model source count and GC absorption | This is the "in the pipeline" half of the question |
| Joint list fixed to Gap 1's six | + best HESTIA and gNFW of each background's neutral arm | Gap 1's list held the halos best under Pohl22 only (amendment 1) |
| q ladder 0.30–1.40 | 0.10–1.40 | under Pohl22 the profile still rose at 0.30 (amendment 2) |
| re-selection: every gas configuration | every value of each gas setting, both source profiles | 27 gas configurations cannot fit in 12 (amendment 3, before any Family 2 fit) |
| run the 12 or the 6 + 6 | run both (23 models) | the re-selection lost the D0 = 40 and zL = 3 extremes; the original 12 hold them |
| joint arm, 8 templates | VVV, F98, best gNFW and best HESTIA of the neutral arm (Family 2 only; I and V kept 8) | ~100 min a model under these backgrounds; the verdict only needs the best of each class (amendment 4) |
