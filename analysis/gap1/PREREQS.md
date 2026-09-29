# Gap 1 baseline fit — what has to exist before step 1

Nothing in Afeefa's five steps can run yet. Four things are missing, and only one
of them is ours to fix. This is the list, in dependency order.

## 1. A machine that can actually run Fermitools

`fermipy` is on PyPI, but it is only a Python wrapper — the actual likelihood
engine (`gtselect`, `gtbin`, `gtsrcmaps`, `gtlike`, `gtpsf`) ships **only** as
the conda `fermitools` package. There is no pip install for it.

```bash
# needs conda/mamba and open internet
conda create -n fermi -c conda-forge -c fermi python=3.11 fermitools
conda activate fermi
pip install fermipy==1.2
python -c "import gt_apps, fermipy; print(fermipy.__version__)"
```

Checked in our cloud sandbox: no conda, and the egress proxy blocks
`conda.anaconda.org`, `fermi.gsfc.nasa.gov` and `heasarc.gsfc.nasa.gov`
(only pypi.org resolves). **So the fit cannot be run there.** It needs a laptop
on normal internet, or a university cluster. Everything else in this repo runs
anywhere.

## 2. The Fermi-LAT photon data

Nobody has downloaded it. From the FSSC data server
(`https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi`), request:

| Field | Value |
|---|---|
| Coordinates | `0, 0` galactic |
| Search radius | `30` deg (a 40x40 deg square ROI needs the corners, at 28.3 deg) |
| Time range (MET) | `239557417` to `760492800` — check these against MJD 54682–60710 |
| Energy range | `500, 50000` MeV |
| Spacecraft data | yes — the `SC00` file is required for exposure |

That returns `*_PH??.fits` photon files plus one `*_SC00.fits`. Expect several GB.
Then the standard cuts (Conventions §1):

```bash
gtselect infile=@events.txt outfile=gap1_filtered.fits \
  ra=266.4050 dec=-28.9362 rad=30 \
  tmin=INDEF tmax=INDEF emin=500 emax=50000 \
  zmax=90 evclass=128 evtype=3

gtmktime scfile=SC00.fits filter="DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
  roicut=no evfile=gap1_filtered.fits outfile=gap1_gti.fits
```

Note `ra=266.4050 dec=-28.9362` is the galactic centre in J2000 — gtselect wants
celestial coordinates even though our ROI is defined in galactic.

## 3. Background models and the catalog

```bash
# diffuse models (FSSC)
wget https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits
wget https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/iso_P8R3_SOURCE_V3_v1.txt
```

Two problems to resolve before running:

- Conventions §4 names the isotropic file `iso_P8R3_SOURCE_V3_V7.txt`. The FSSC
  filename is `iso_P8R3_SOURCE_V3_v1.txt`. One of the two is wrong — confirm
  which before anyone hard-codes it.
- **FL16Y** is the 16-year source list (arXiv:2602.22148). We need the catalog
  FITS, and Fermipy's `catalogs` setting has to accept it — it ships knowing
  4FGL-DR3/DR4, so FL16Y may need to be passed as a custom catalog file.

Fermi Bubbles templates (sharp and fuzzy) are supposed to live in
`templates/bubbles/` per Conventions §4. They do not exist yet and nobody owns
making them.

## 4. The VVV bulge template — SOLVED, no email needed

It is public: **github.com/chrisgordon1/galactic_bulge_templates**, the official
template release for Coleman et al. 2019 (arXiv:1911.04714). Cloned and in
`templates/bulge/`. Two files, both from Table 5 of that paper:

- `Bulge_modulated_Coleman_etal_2019_Normalized.fits` — the **non-parametric
  (VVV) bulge**. This is the one Afeefa's step 1 wants.
- `F98_BoxyBulge_arxv1611.06644_Normalized.fits` — Freudenreich 1998 S-bulge,
  the parametric boxy model used by Macias et al. 2016. Free bonus: an
  independent boxy template, useful as a positive control for the c4 statistic.

Both are 200 x 200 at 0.2 deg on a GLON/GLAT-CAR grid covering exactly our
40 x 40 deg ROI, and both are normalised to unit integral with a flat pixel
area — the same convention Muru used. Their `CDELT1` is already negative, so
unlike the HESTIA files they need **no longitude flip**. Their pixel centres
sit half of *our* pixel off ours, so the regrid to 0.1 deg is an interpolation,
not a 2x block upsample. Handled in `analysis/gap1/bulge_to_template.py`.

Caveat worth knowing: because these maps are exactly ROI-sized, there is no
"uncropped" view of them, so the ROI crop artifact cannot be checked against a
wider map the way it can for HESTIA. Stay inside the radius guard.

## 5. What is ready

- 24 HESTIA templates on the conventions grid, in `templates/`
- Coleman (VVV) and F98 bulge templates on the same grid
- gNFW builder, bulge loader, HESTIA loader, shape statistics — all self-testing
- `utils/psf.py` — shared PSF, shape-measurement use only (see below)
- `masks/make_masks.py` — masks A, B and C as gtlike weights maps
- Config and fit driver, written and reviewed but never executed

### The run, once Fermitools exists

```bash
python analysis/gap1/run_fit.py \
  --config configs/config_gap1_maskA_iem1.yaml \
  --gce hestia=templates/gap1_hestia_G1.1_ang0_1deg.fits \
  --gce coleman=templates/gap1_bulge_coleman.fits \
  --gce gnfw=templates/gnfw_gamma1.2.fits \
  --wide hestia=<path to Muru's 800x800 density_projection_..._1deg.fits> \
  --out outputs/gap1/baseline_fit.json
```

Each template is fitted against the *same* stage-1 background, which is what
makes the delta-lnL comparable. `--wide` is optional and only feeds the
uncropped shape column.

## The PSF convolution problem — raise this with the team

Conventions §5 said every template must be PSF-convolved with a shared
`utils/psf.py` before use. **For anything fed to Fermipy that is wrong**, and
Afeefa has now updated the conventions: raw templates go into Fermipy directly,
and `utils/psf.py` is used only for shape measurement, SBI inputs and plots.
`gtsrcmaps` convolves a `SpatialMap` template with the PSF and exposure itself
when it builds the source maps. Handing it a pre-convolved template blurs the
signal twice, which makes every template rounder — and "rounder" is exactly the
thing Gap 1 and Gap 4 are measuring.

So the rule should be:

- **Templates going into a fit:** intrinsic sky maps, NOT pre-convolved. Let
  gtsrcmaps do it.
- **Templates used for shape measurement** (comparing q and c4 at LAT
  resolution): convolve with a shared PSF, because there gtsrcmaps is not in
  the loop.

Our template FITS headers already record `PSFCONV = False` so this cannot be
lost track of. Worth confirming against the Fermipy docs and then fixing the
conventions, since it affects all four gap teams identically.
