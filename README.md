## GCE Reconciliation Project

Yo guys, welcome to the central repository for the GCE (Galactic Center Excess) Reconciliation Project.

**Important:** This repository houses the code skeleton and shared architecture. For the complete, agreed-upon parameters, statistical methods, and dataset cuts, **the canonical source is the shared Google Doc (Project Conventions)**. For convenience, the full content of those conventions is also replicated in this README.

---

## Start here: the diffuse-freedom test

**`docs/diffuse_freedom_test.md`.** Macias and Di Mauro use the same bulge
templates and the same gas maps and reach opposite conclusions about dark
matter. That disagreement is reproduced here inside one pipeline on one
dataset: a defensible change in how much the interstellar model is allowed to
bend moves a detection from TS 3.8 to 306, with injection controls showing the
test is not simply eating the signal.

It is a shared module, `analysis/diffuse_freedom/`, meant to be run by every
gap. Splits are fixed in `configs/diffuse_freedom.yaml` before any run and
applied identically to every template, on the same footing as
`configs/level2.yaml`.

Two constraints from it that belong in the Conventions:

* Mask D with a single-template IEM has almost no sensitivity and Mask A has
  none at all, so Gap 1 cannot use either until a ring-decomposed IEM is
  loaded. Gap 5 can, its IEMs are ring-decomposed by construction. Gaps 4 and
  7 inherit it.
* Daylan et al. 2016 define the axis ratio as the reciprocal of Di Mauro's and
  Cholis's. Pooling them without inverting Daylan manufactures a trend.

## Where the large data comes from

Raw Fermi products and the external IEM libraries are not committed.
`docs/data_provenance.md` says what each one is, where it comes from and
whether it is reproducible, downloadable or unavailable.

---

## Directory Structure & What Goes Where

Right now, folders are empty (containing only hidden `.gitkeep` files to preserve the structure on GitHub). As work progresses, use the following layout:

* **`analysis/`**: Contains subdirectories for each gap team (`gap1`, `gap4`, `gap5`, `gap7`). Place all team-specific Python scripts, Jupyter notebooks, and local `README.md` files (following the Mini-Paper template) inside your designated folder.
* **`configs/`**: Stores YAML configuration files for pipeline runs. Your local `paths.yaml` configuration file will live here.
* **`data/`**: The designated location for your raw Fermi-LAT data. This folder is ignored by Git, ensuring heavy datasets stay strictly on your local machine or cluster.
* **`docs/`**: Shared documentation, reference papers, and notes.
* **`iem/`**: Interstellar Emission Model files (GALPROP, McDermott extreme library, Galp21/Pohl22) and associated helper scripts.
* **`masks/`**: Shared mask definitions (Mask A, B, and C) and code for generating them.
* **`outputs/`**: Local destination for all generated results (FITS, HDF5, CSV) and logs. Ignored by Git.
* **`templates/`**: Shared spatial templates (gNFW, HESTIA, VVV, etc.). Fermi Bubbles templates (sharp and fuzzy variants) are stored inside `templates/bubbles/`.
* **`utils/`**: Shared utility scripts. (`psf.py`) is for shape measurement, SBI inputs, and plots. It is not used for fitting, since Fermipy handles PSF convolution automatically via gtsrcmaps.



---

## Getting Started: Local Setup

Because raw datasets and heavy outputs are not tracked by Git, every team member must configure their local environment after cloning.

### 1. Clone the Repository

```bash
git clone https://github.com/YOUR-USERNAME/gce-reconciliation.git
cd gce-reconciliation
```

### 2. Configure Local Paths

Do not hard-code absolute file paths (such as `/Users/name/Desktop/...`) in your analysis scripts.

* Navigate to the `configs/` folder.
* Duplicate `paths_template.yaml` and rename the copy to `paths.yaml`.
* Update `paths.yaml` with the absolute paths pointing to your local Fermi data, IEM models, and output directories.
* *Note: `paths.yaml` is intentionally ignored by Git to prevent overwriting teammates' local configurations.*

### 3. Link Your Data

Place your Fermi-LAT data directly in the `data/` folder, or create a symlink pointing to where the data is stored on your machine.

### 4. Reference the Google Doc
The full Project Conventions and Mini-Paper Template are in the Google Doc. Always check the Google Doc for any updates.

All the best for all of us! 

