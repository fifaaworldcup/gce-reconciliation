# Gap 5 — the Cholis et al. (2022) models through Level 2

Question: Gap 1 found the VVV bulge beats spherical gNFW and HESTIA at Level 2
on Pohl22's ring background. Does that verdict change when the background is
swapped for the extreme GALPROP models of Cholis et al. (2022)?

All rules are fixed in `configs/gap5.yaml` before any Family 2 model is fitted.
The mini-paper is `docs/gap5_minipaper.md`.

## Files

| file | what it does |
|---|---|
| `convert_iem.py` | GALPROP maps → one cube per component (π⁰, bremss, ICS) on our grid. Reads CAR, HEALPix and the Cholis release layout. `--cholis <ZENODO_FILES> --models ...`, `--check-against` for the orientation test, `--selftest` |
| `select_models.py` | Table VIII rebuilt from the 80 GALDEF files, and the re-selection (6 Table II anchors + coverage + farthest point) |
| `export_reference_components.py` | writes gll_iem_v07 and Pohl22 rings_species in the same format. `--check` |
| `run_gap5.py` | driver: level1, admit, q, l2, controls per model/mask/energy |
| `summarise_gap5.py` | table + verdict from the pre-registered rule |
| `iem_selection_audit.py` | audit of the 12-model selection |

`analysis/level2/run_level2.py` gained `--components DIR` (one free
normalisation per cube in DIR), `--label` and `--checkpoint`. `iterate.run`
takes an optional checkpoint file. With checkpointing off, Level 2 is unchanged
(the Gap 1 reproduction below was run that way). With it on, a killed run
resumes exactly: same sources at the same positions, same lnL and TS to the
last digit (tested by killing a gll_iem_v07 neutral run after 3 sources).

## The Cholis release (Zenodo 6423495)

The maps carry no WCS. They are 240×240 over |l|, |b| < 30° at 0.25°, with 38
planes log-spaced from 50 MeV to 814 GeV. Longitude rises with column,
latitude with row, not transposed. This layout was found by testing all 8
orientations against gll_iem_v07 at 0.9 GeV: corr 0.991, next best 0.907.
Every converted model is re-checked (0.979–0.991 for the 23 run). The model
↔ file code mapping is `NAMING_CONVENTION_OF_DIFFUSE_EMISSION_MODELS.dat`.

## Reproduce

```bash
D=<folder with gap1_ccube.fits, gap1_expcube.fits, gll_psc_catalog.fit, gap1_iem_roi.fits>

# 1. reference backgrounds, plus the check that the component path = fit_rings
python analysis/gap5/export_reference_components.py --ccube $D/gap1_ccube.fits \
    --iem $D/gap1_iem_roi.fits --rings <pohl22 folder> --check \
    --expcube $D/gap1_expcube.fits --catalog $D/gll_psc_catalog.fit

# 2. each Cholis model (one model's files at a time)
python analysis/gap5/select_models.py --zenodo <ZENODO_FILES> --out outputs/gap5/selection
python analysis/gap5/convert_iem.py --cholis <ZENODO_FILES> --models I,II,... \
    --out outputs/gap5/components --check-against $D/gap1_iem_roi.fits

# 3. everything
python analysis/gap5/run_gap5.py --ccube $D/gap1_ccube.fits \
    --expcube $D/gap1_expcube.fits --catalog $D/gll_psc_catalog.fit \
    --models refs,family2 --steps level1,admit,q,l2
python analysis/gap5/run_gap5.py ... --models family2 --masks B --steps controls
python analysis/gap5/summarise_gap5.py
```

Outputs go to `outputs/gap5/<model>/<mask>/<emin>-<emax>/` (not committed).
A Level 2 run is about 9 minutes per arm on 2 cores.

## Validation

Pohl22 rings_species through `--components` against Gap 1's Level 2 joint arm
(mask B + PS, 1–10 GeV, 8 bins):

| template | Gap 1 lnL / sources | Gap 5 path lnL / sources |
|---|---|---|
| vvv | 554,484.945 / 10 | 554,484.95 / 10 |
| f98 | 554,354.816 / 10 | 554,354.82 / 10 |
| gnfw_g1.2 | 554,235.082 / 13 | 554,235.08 / 13 |
| gnfw_g1.4 | 554,288.988 / 13 | 554,288.99 / 13 |
| hestia_G3.2_ang90 | 554,182.346 / 10 | 554,182.35 / 10 |
| hestia_G2.2_ang90 | 554,092.057 / 10 | 554,092.06 / 10 |

Same source positions in every row.

Level 1 for the same background, the two code paths: lnL differs by 2e-5.
gll_iem_v07 admissibility through `run_gap5.py` gives the same numbers as Gap 1's
ring table (vvv 184.3 / +585, gnfw_g1.2 109.3 / +320).

q at Level 1 uses the 8-bin cube Level 2 runs on. Gap 1's `measure_q` used 4
bins over the same range, so the two q values are not expected to agree
exactly (gll, mask B: 0.382 here, 0.359 in Gap 1).

## Status

- [x] pipeline built and validated against Gap 1
- [x] both reference backgrounds run (masks B and D, 1–10 GeV; mask B at 0.5–3 and 3–50) — see `outputs/gap5/gap5_table.md` and the mini-paper §6.3
- [x] Cholis maps converted, 23 models (re-selected 12 + original 12) run on mask B, 1–10 GeV, both arms
- [x] ordering controls on V, XVII, LI, LVI (gNFW, VVV, nothing injected)
- [x] Table VIII rebuilt from the GALDEF files; re-selection in `outputs/gap5/selection/`
- [ ] Family 2 on mask D and the secondary energy ranges (not run: both
      references already fail gNFW admissibility under mask D and at 3–50 GeV)
