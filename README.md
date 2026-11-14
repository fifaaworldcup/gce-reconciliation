# GCE Reconciliation

**Code and analysis accompanying the manuscript:**

> *A Quantitative Reconciliation of Conflicting Galactic Center Excess Inferences*
> A. Anwar, O. Pernica, P. Ravichandran, H. Dalal, T. Shen

---

## ⚠️ Status

- **The manuscript is not yet on arXiv.**
- **It has not been published or peer reviewed.**

Nothing here should be treated as a citable result yet. When the
manuscript is posted to arXiv, this section will be updated with the arXiv
identifier.

---

## What this repository is

This is the analysis code for a paper that studies the **Galactic Center
Gamma-Ray Excess (GCE)** — an unexplained excess of GeV gamma rays from the
inner Galaxy — and asks a narrow question:

> **How much of the disagreement in the literature is caused by analysis
> choices, rather than by the data?**

The paper does not argue for dark matter, millisecond pulsars, or background
mismodeling. It measures how the inferred result moves when the analysis
choices are varied.

The two conclusions reported in the manuscript are:

1. **The excess is robustly detected.** It is recovered under all 25
   interstellar emission models (IEMs) tested.
2. **Its morphology is not robust.** The preferred spatial template changes
   with the IEM. Most of the variation in the stellar-bulge templates is
   attributable to the choice of IEM.

---

## What this repository is **not**

- It is **not** a published result. There is no arXiv ID or DOI yet.
- It is **not** a general-purpose Fermi-LAT pipeline. It is the code that
  produced one paper.
- It is **not** a data repository. Fermi-LAT data and the large external IEM
  libraries are **not** included here; see *Data* below.
- It is **not** complete. Folders are being filled in as the paper is
  finalized.

---

## Repository layout

The structure below is the intended layout. **Folders marked "(pending)" are
not yet fully populated.**

| Folder | Purpose | Status |
|---|---|---|
| `analysis/` | Gap 1, Gap 5, Gap 7 pipelines; Level 2 source finding | pending |
| `configs/` | Analysis configuration files and pre-registered thresholds | pending |
| `iem/` | Interstellar emission model interfaces | pending |
| `masks/` | Mask definitions and mask audit tool | pending |
| `templates/` | Spatial templates (gNFW, HESTIA, VVV, F98, bar) | pending |
| `docs/` | Notes on data provenance and each analysis component | pending |
| `outputs/` | Derived products (not tracked by Git) | not tracked |
| `utils/` | Shared helper code | pending |

---

## Data and external dependencies

This repository does **not** contain raw data. To run the code you will need
to obtain, separately:

- **Fermi-LAT data** — public, from the
  [Fermi Science Support Center](https://fermi.gsfc.nasa.gov/ssc/data/access/).
- **Interstellar emission models** — e.g. the Cholis et al. GALDEF library
  (Zenodo record 6423495) and the Fermi diffuse model `gll_iem_v07`.
- **HESTIA dark matter maps and particle data** — provided by the HESTIA
  authors; see `docs/` once populated.
- **Template maps** — VVV bulge, F98, Galactic bar; see `docs/`.

Local paths are not committed. Once `configs/` is populated, copy the
template path file and edit it for your machine:

```bash
cp configs/paths_template.yaml configs/paths.yaml
# then edit configs/paths.yaml
