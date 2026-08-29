## GCE Reconciliation Project

Yo guys, welcome to the central repository for the GCE (Galactic Center Excess) Reconciliation Project.

**Important:** This repository houses the code skeleton and shared architecture. For the complete, agreed-upon parameters, statistical methods, and dataset cuts, **the canonical source is the shared Google Doc (Project Conventions)**. For convenience, the full content of those conventions is also replicated in this README.

---

## Directory Structure & What Goes Where

Right now, most of these folders are empty (containing only hidden `.gitkeep` files to preserve the structure on GitHub). As work progresses, use the following layout:

* **`analysis/`**: Contains subdirectories for each gap team (`gap1`, `gap4`, `gap5`, `gap7`). Place all team-specific Python scripts, Jupyter notebooks, and local `README.md` files (following the Mini-Paper template) inside your designated folder.
* **`configs/`**: Stores YAML configuration files for pipeline runs. Your local `paths.yaml` configuration file will live here.
* **`data/`**: The designated location for your raw Fermi-LAT data. This folder is ignored by Git, ensuring heavy datasets stay strictly on your local machine or cluster.
* **`docs/`**: Shared documentation, reference papers, and notes.
* **`iem/`**: Interstellar Emission Model files (GALPROP, McDermott extreme library, Galp21/Pohl22) and associated helper scripts.
* **`masks/`**: Shared mask definitions (Mask A, B, and C) and code for generating them.
* **`outputs/`**: Local destination for all generated results (FITS, HDF5, CSV) and logs. Ignored by Git.
* **`templates/`**: Shared spatial templates (gNFW, HESTIA, VVV, etc.). Fermi Bubbles templates (sharp and fuzzy variants) are stored inside `templates/bubbles/`.
* **`utils/`**: Shared utility scripts. The standard PSF convolution function (`psf.py`) lives here so all gap teams use the exact same implementation.

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

The full Project Conventions and Mini-Paper Template are provided below for convenience, but the Google Doc remains the authoritative source. Always check the Google Doc for any updates.
```
