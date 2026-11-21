# A Quantitative Reconciliation of Conflicting Galactic Center Excess Inferences

**Repository:** `fifaaworldcup/gce-reconciliation`  
**Manuscript:** *A Quantitative Reconciliation of Conflicting Galactic Center Excess Inferences*  
**Authors:** Afeefa Anwar, Oren Pernica, Prithiga Ravichandran, Hussain Dalal, and Tianbo Shen *(equal contribution: Afeefa Anwar and Oren Pernica)*

This repository contains the analysis code, configuration files, spatial templates, documentation, and derived products for a systematic investigation into how modeling and data-analysis choices affect inferences about the Galactic Center Gamma-Ray Excess (GCE) in *Fermi*-LAT data.

> **Main Conclusion:** The GCE is robustly detected across all tested interstellar emission models (IEMs), but its inferred spatial morphology is **not** robust. The preferred signal template shifts with the diffuse-background model and, in several cases, with the iterative source-finding procedure. Current morphological evidence does not decisively establish either a dark-matter or a stellar-bulge origin.
>
> **Manuscript Precedence:** All numerical results and interpretations below follow the submitted manuscript. Where older component notes or intermediate outputs report earlier revision numbers, treat the manuscript values as authoritative.

---

## Table of Contents

1. [Research Overview](#1-research-overview)
2. [Key Results](#2-key-results)
3. [Data, Masks, Templates, and IEMs](#3-data-masks-templates-and-iems)
4. [Repository Structure](#4-repository-structure)
5. [Reproducing the Analyses](#5-reproducing-the-analyses)
6. [External Data and Provenance](#6-external-data-and-provenance)
7. [Software Stack and Implemented Modules](#7-software-stack-and-implemented-modules)
8. [Scope, Caveats, and Limitations](#8-scope-caveats-and-limitations)
9. [Citation and Acknowledgements](#9-citation-and-acknowledgements)

---

## 1. Research Overview

The GCE is an extended excess of GeV gamma rays toward the inner Milky Way whose physical origin remains unresolved. Leading hypotheses include a smooth dark-matter distribution, unresolved millisecond pulsars associated with the Galactic bulge, and residuals from imperfect modeling of diffuse Galactic emission.

Because the excess is subdominant to the Galactic diffuse foreground, choices across foreground modeling, point-source/plane masking, source finding, and signal morphology can produce conflicting conclusions on the same photon dataset. This project evaluates those dependencies under a unified framework using a single 16.5-year *Fermi*-LAT selection across six methodological axes:

1. **Morphology:** Spherical generalized-NFW (gNFW) profiles, flattened HESTIA dark-matter projections, and stellar-bulge templates.
2. **Masking:** Alternative point-source and Galactic-plane masks audited against 4FGL-DR4.
3. **Interstellar Emission Model (IEM):** Standard *Fermi* diffuse model (`gll_iem_v07`), Cholis et al. (2022) models, and ring-decomposed Pohl22 components.
4. **Source Finding:** Single-pass fitting versus staged iterative source finding.
5. **Energy:** Full-band and energy-resolved comparisons.
6. **Diffuse Freedom:** Alternative parameterizations controlling spatial flexibility of the diffuse background.

Every apparent detection or null result is paired with injection/admissibility controls to verify whether the pipeline can recover a known signal under that configuration.

---

## 2. Key Results

### 2.1 Excess Recovery Across All 25 Tested IEMs
In the common **Mask B**, **1–10 GeV** Level 2 comparison, the GCE is detected across **all 25 IEMs** evaluated: 23 models from the Cholis et al. (2022) library, the ring-decomposed Pohl22 model, and the standard *Fermi* diffuse model `gll_iem_v07`. Detection is consistent across these backgrounds, whereas inferred morphology varies markedly.

### 2.2 Sensitivity to Diffuse-Background Freedom
Under **Mask B** using `gll_iem_v07` and the VVV stellar-bulge template, six diffuse-freedom configurations were evaluated against pre-registered admissibility rules: a configuration is admissible only if it recovers an injected signal at **$\ge 50\%$ of the rigid-split TS** and reaches **$\text{TS} \ge 25$**.

| Diffuse-Freedom Configuration | Extra Nuisance Params | Observed VVV TS | Injected Signal Retained | Admissible? |
| :--- | :---: | :---: | :---: | :---: |
| `rigid` | 0 | 264.4 | 100% | Yes |
| `ns` | 1 | 159.4 | 97% | Yes |
| `lat_outside` | 4 | 305.6 | 81% | Yes |
| `lat_inside` | 4 | 113.7 | 47% | No (fails retention criterion) |
| `lat` | 3 | 8.7 | 80% | Yes |
| `lat_ns` | 7 | 8.7 | 79% | Yes |

Across the **five admissible** configurations, the observed VVV TS spans **8.7–305.6** (an approximate 30-fold spread reported in the manuscript). Granting spatial freedom outside central longitudes increases fitted TS, whereas full-region latitude-dependent freedom reduces it.

* Under **Mask D**, three usable splits span TS **4.6–79.9** (a ~17-fold range).
* Under **Mask A**, no split is admissible; a low-energy mask audit shows Mask A retains only **0.4%** of the region of interest (ROI) at **0.59 GeV**. **Mask C** similarly has a restrictive low-energy footprint and serves as a diagnostic mask only.

### 2.3 Pohl22 Bulge Preference vs. 25-IEM Scan
* **Level 1 (Pohl22 ring-decomposed):** The VVV stellar-bulge template outperforms every tested dark-matter template; the best HESTIA template trails VVV by **$\Delta\text{AIC} = +888$** (where $\Delta\text{AIC} = \text{AIC}_{\text{template}} - \text{AIC}_{\text{ref}}$, so positive values indicate disfavor relative to the reference).
* **Level 2 (25-IEM common scan):** Pohl22 is the **worst-fitting** of the 25 backgrounds, with every Cholis model preferred over Pohl22 by **$\Delta\text{AIC}$ of 3,904 to 16,491**.
  * Across the **23 Cholis models**, stellar-bulge templates win in **3/23** cases with a common source list and **4/23** with per-hypothesis source finding.
  * Spherical **gNFW** wins in **18/23** cases with a common source list and **13/23** with per-hypothesis source finding; flattened **HESTIA** templates win the remaining cases.
  * Among the six best-fitting backgrounds (**LIII, X, V, VII, XV, LVI**), morphological winners do not converge; four sit near the resolution limit set by injection controls, and only **LIII** and **LVI** are decisive in the reported comparison.

### 2.4 Source Finding as a Second Axis of Variation
Applying Level 2 iterative source finding substantially reduces individual template TS values while preserving the Pohl22 template ordering (for example, shifting VVV TS from **2151.7** to **1018.2**).

* Across the 23 Cholis models, running source finding **per hypothesis** rather than using a **common source list** flips the winning morphology in **5/23** cases. Every flip favors the hypothesis whose search adds more point sources ("source absorption").
* On real data, the Level 2 run adds **12 sources**: **9** match a 4FGL-DR4 source within $0.3^\circ$, **3** are uncertain, and **0** are unmatched. Removing the one uncertain candidate near the Galactic Center leaves the template ordering unchanged. On a null simulation, **0** sources are added.

### 2.5 Variance Decomposition Across Template Classes
A descriptive sum-of-squares decomposition of $\Delta\text{AIC}$ across **960 observations** (24 diffuse backgrounds $\times$ 40 templates across six HESTIA halos) partitions variation into background main effect, template main effect, and interaction remainder:

| Template Class | Background Share | Template Share | Interaction Remainder |
| :--- | :---: | :---: | :---: |
| **Stellar (VVV, F98)** | 96.89% | 1.87% | 1.24% |
| **HESTIA simulated dark matter** | 62.72% | 14.43% | 22.85% |
| **Spherical gNFW** | 20.94% | 57.40% | 21.66% |

The pooled background main effect across all classes is **33.48%**. Because HESTIA and spherical gNFW exhibit non-negligible interaction remainders, these shares are descriptive summaries of observed sums of squares rather than clean, isolated causal parameters. A six-halo bootstrap on HESTIA yields a median background share of **64.79%** (95% interval **54.85%–78.86%**), which is indicative given the six-halo sample size.

### 2.6 Signal-Injection Controls
* **Direct morphology injection:** Across six synthetic skies (VVV, HESTIA G3.2, or spherical gNFW at two flux levels), the pipeline recovers the true morphology class in **6/6** cases and ranks the exact injected template first in **4/6** (the two misses land on neighboring templates within $\Delta\text{TS} \le 2.5$). Dark-matter injections yield dark-matter winners in every case.
* **Flux sensitivity asymmetry:** The pipeline is roughly **1.7$\times$ less sensitive in flux** to the tested dark-matter morphologies than to bulge morphologies. This asymmetry is reported directly rather than rescaled.

---

## 3. Data, Masks, Templates, and IEMs

### 3.1 *Fermi*-LAT Data Selection
* **Observation window:** 16.5 years, 2008 August 4 to February 2025 (MJD 54682–60710).
* **Event selection & ROI:** `P8R3_SOURCE_V3`, `evtype=3`, within $\vert{}l\vert{} < 20^\circ$, $\vert{}b\vert{} < 20^\circ$.
* **Energy binning:** 0.5–50 GeV split into **0.5–1**, **1–3**, **3–10**, and **10–50 GeV** bins.
* **Quality & transient cuts:** Zenith angle $< 90^\circ$, rocking angle $< 52^\circ$, `DATA_QUAL > 0`, `LAT_CONFIG == 1`, and exclusion of bright transients via the `FL16Y` interval list.
* **Cross-checks:** `P8R3_CLEAN` above 10 GeV for cosmic-ray residual checks; a Level 2 comparison against Di Mauro (2026) uses `P8R3_SOURCEVETO`, $0.08^\circ$ spatial pixels, 1–10 GeV for morphology, and 0.5–1000 GeV for the spectral-energy distribution.

### 3.2 Masks (Audited Against 4FGL-DR4)
* **Mask A:** Aggressive 95%-containment source masks plus $0.5^\circ$; structurally inadmissible at low GCE energies (0.4% ROI retained at 0.59 GeV).
* **Mask B:** Minimal source mask; primary baseline comparison mask.
* **Mask C:** Galactic-plane robustness mask; diagnostic only due to low-energy area loss.
* **Mask D:** TS- and flux-scaled mask following Cholis et al. (2022) with energy-dependent radii and inner-Galaxy disc cut.
* **Mask E:** Flux-scaled cross-check superseded by Mask D.

### 3.3 Spatial Templates
* **Spherical gNFW:** Scale radius $r_s = 20\text{ kpc}$, Solar radius $R_\odot = 8.2\text{ kpc}$, inner slope $\gamma \in [0.8, 1.4]$.
* **Flattened HESTIA dark matter:** Density-squared projections from six constrained simulations across four viewing angles, standardized to $1^\circ$ smoothing. Direct J-factor checks show axis-ratio differences of $-0.032 \pm 0.020$ relative to squared-column estimates (sign-consistent in 22/24 maps) and $-0.041 \pm 0.016$ in the disc frame (24/24).
* **VVV stellar bulge:** Maximum-entropy stellar-mass map from VVV red-clump data (boxiness $c_4 = +0.009$).
* **F98 boxy bulge:** Parametric boxy-bulge reference ($c_4 = +0.024$).
* **Galactic bar:** Included in the broader template set, outside the paper's covariance analysis.

### 3.4 Interstellar Emission Models & Likelihood Fitting
* **25-IEM scan:** `gll_iem_v07`, 23 Cholis et al. (2022) GALDEF models, and ring-decomposed Pohl22.
* **Galp21 (`SA0`, `SA50`, `SA100`):** Referenced in literature comparisons; precomputed sky maps are not bundled with public GALDEF configurations.
* **Fitting & convergence:** Binned Poisson likelihood via Fermipy and `MINUIT2`, reporting TS and $\Delta\text{AIC}$ with explicit reference definitions, plus a per-bin log-likelihood gap bound to catch boundary-frozen active sets beyond boolean convergence flags.

---

## 4. Repository Structure

```text
.
├── analysis/
│   ├── diffuse_freedom/   # diffuse-background flexibility and admissibility tests
│   ├── level2/            # iterative source-finding module
│   ├── gap1/              # GCE template fits, data preparation, and related analyses
│   └── gap7/              # mask–IEM covariance and variance decomposition
├── configs/               # fixed analysis settings and run configurations
├── data/                  # local Fermi-LAT inputs (large/raw inputs untracked)
├── docs/                  # methods, data provenance, and analysis notes
├── iem/                   # external IEM libraries (untracked)
├── masks/                 # source and ROI masks
├── outputs/               # derived tables, figures, and compact result products
├── scripts/               # data acquisition and preparation utilities
├── templates/             # spatial templates used by the analyses
└── utils/                 # shared utilities
```

### Key Paths and Entry Points

| Path | Role |
| :--- | :--- |
| `docs/data_provenance.md` | Sources, sizes, and availability of external datasets. |
| `docs/diffuse_freedom_test.md` | Detailed write-up and reproduction commands for the diffuse-freedom test (check Section 2.2 & Section 8 for final manuscript numbers after projected-Newton fixes). |
| `analysis/diffuse_freedom/` | Implementation of the diffuse-freedom and admissibility tests. |
| `analysis/level2/` | Iterative source-finding implementation. |
| `analysis/gap7/` | Mask–IEM covariance and variance-decomposition scripts. |
| `configs/level2.yaml` | Locked Level 2 source-finding parameters. |
| `outputs/` | Derived tables, figures, and result products supporting the manuscript. |

---

## 5. Reproducing the Analyses

A fresh clone will not reproduce every scan out of the box because raw *Fermi*-LAT photon/spacecraft files and large external diffuse/simulation archives are not stored in Git. Download the required inputs (Section 6), verify paths in `docs/data_provenance.md`, and run the module scripts.

### Diffuse-Freedom Test Entry Point

```bash
python analysis/diffuse_freedom/run_test.py \
  --ccube gap1_ccube_full.fits \
  --expcube gap1_expcube_full.fits \
  --iem gap1_iem_roi.fits \
  --catalog gll_psc_catalog.fit \
  --mask B \
  --ps \
  --out outputs/gap1/diffuse_freedom_maskB_ps
```

See `docs/diffuse_freedom_test.md` for figure-generation commands and JSON output paths. Run Level 2 and `gap7` covariance/variance scripts from their respective `analysis/` subdirectories using the committed configuration files.

### Fixed Level 2 Configuration (`configs/level2.yaml`)

To avoid hypothesis-dependent tuning, all templates, masks, and IEMs share fixed stopping rules:
* **Passes and TS thresholds:** 3 passes at $\text{TS} = 100$, $64$, and $36$.
* **Per-pass cap:** Up to 40 added sources per pass.
* **Spatial constraints:** $0.30^\circ$ minimum source separation; $1.0^\circ$ ROI-edge margin.
* **Refitting & RNG:** Full refit after each added source; random seed `42`.

---

## 6. External Data and Provenance

| Input | Source | Notes / Availability |
| :--- | :--- | :--- |
| *Fermi*-LAT photon, spacecraft, counts, exposure, and livetime cubes | *Fermi* Science Support Center (FSSC) | Public; prepared via repository data-prep scripts/docs. |
| `gll_iem_v07.fits` and `iso_P8R3_SOURCE_V3_v1.txt` | *Fermi* diffuse model products | Public; helper script at `scripts/download_aux.sh`. |
| 4FGL-DR4 catalog | *Fermi* Science Support Center (FSSC) | Public point-source reference catalog. |
| Cholis et al. (2022) IEM library | Zenodo record `6423495` | 80-model archive; download separately. |
| Pohl22 gas-ring and related maps | Zenodo record `6276721` | Download separately. |
| Galp21 (`SA0`, `SA50`, `SA100`) maps | GALPROP v57 / Galp21 authors | Precomputed sky maps not bundled in public GALDEF files. |
| HESTIA J-factor projections & particle data | Muru et al. / HESTIA collaboration | Obtain per `docs/data_provenance.md` or authors. |

> **HESTIA File-Naming Warning:** In the supplied HESTIA particle files, filenames are swapped within the `G1` and `G3` pairs. Do **not** identify halos by filename alone—match each dark-matter particle centroid to its corresponding AHF halo center (`scripts/hestia/whichcentre.py`, applied during reduction in `scripts/hestia/reduce.py`).

`.gitignore` excludes raw data, large IEM archives, and local paths while selectively tracking key figures, tables, and compact outputs needed for the paper.

---

## 7. Software Stack and Implemented Modules

* **Core likelihood & analysis stack:** *Fermi* Science Tools, **Fermipy**, **MINUIT2**, **healpy**, and **NPTFit**.
* **Recorded validation environment:** Manuscript Level 2 validation notes Fermipy `v1.2.0` and Fermitools `v2.2.0`, with nearest documented Fermitools `2.5.3` on Python `3.12`.
* **Unreported framework modules:** The broader codebase also contains modules for NPTF, joint spatial–spectral neural simulation-based inference (SBI), and catalog evolution. Those modules are **not** used to generate the results reported in the manuscript.

---

## 8. Scope, Caveats, and Limitations

* **25-IEM scan scope:** Reported for **Mask B** and **1–10 GeV**. Mask D and multi-band extensions were omitted for Cholis models because the gNFW injection control failed under those settings.
* **Diffuse-freedom splits:** Latitude-band divisions are diagnostic parameterizations rather than claims of physical gas structure.
* **Descriptive variance decomposition:** Reported percentages reflect observed $\Delta\text{AIC}$ sum-of-squares shares, not random-effects variance components.
* **HESTIA reconstruction differences:** Reconstructed HESTIA maps differ from published maps by an average axis-ratio offset of **0.135** with a per-galaxy ordering Spearman correlation of **$-0.77$**, and axis-ratio correlations vary across isophote, crop, and per-angle fitting choices. Treat comparisons as internal differential results of this pipeline.
* **Open questions & unresolved disputes:** This study does not resolve the high-energy disagreement between Manconi et al. and Di Mauro, nor does it establish the physical origin of the GCE.
* **Interpreting low TS:** A low TS in a configuration that fails signal-injection recovery does **not** demonstrate the absence of a signal.

---

## 9. Citation and Acknowledgements

If you use this code, configurations, or derived products, please cite the associated manuscript and the specific repository commit used (updating this section with the final journal reference, DOI, or preprint ID once assigned):

* **Title:** *A Quantitative Reconciliation of Conflicting Galactic Center Excess Inferences*
* **Authors:** Afeefa Anwar, Oren Pernica, Prithiga Ravichandran, Hussain Dalal, and Tianbo Shen

We thank Moorits Mihkel Muru for providing HESTIA density maps and particle data and for correspondence on simulation conventions, and Troy Porter and Guðlaugur Jóhannesson for guidance on Galp21 configurations. This work relies on the *Fermi* Science Tools, Fermipy, NPTFit, and healpy.
```
