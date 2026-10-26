# covariance_input.csv

Gap 1 rows (`outputs/gap1/gap1_results_for_gap7.csv`) plus the Gap 5
Level 2 rows. Same columns as the Gap 1 file, with three added at the end:

- `arm`: `neutral` (one common source list, all 40 templates) or `joint`
  (each template runs its own source search). Empty for Gap 1 rows.
- `iem_set`: `reference`, `reselected`, `original`, or `reselected+original` (XXVII).
- `fit_gap_vs_best_iem`: best AIC under that background minus best AIC under
  any background (mask B, 1-10 GeV). gll_iem_v07 is 0.

`iem_family` for Gap 5 rows: 1 = gll_iem_v07, 2 = Pohl22 rings_species,
3 = Cholis et al. (2022). The Gap 5 Pohl22 rows reproduce the Gap 1 Level 2
rows (same lnL to 0.01), so filter on `gap` to avoid counting them twice.

Use rows with `status` starting `OK` and `admissible` true, and always filter
on `energy_bin` as well as `mask`: the same background appears at several
energy ranges.

## What `admissible` means

Gap 1 is per template. Each row's configuration (background, mask, energy)
was given that template injected at the flux it fits under gll_iem_v07, and
the template had to come back with injected TS >= 25 and boosted gain >= 25.
`control_template` says which control decided the row. For mask B and D at
1-10 GeV every template has its own control. Elsewhere a control of the same
class is used (gNFW for gNFW, HESTIA for HESTIA, bulge for bulge); a control
of another class never counts. If there is none the field is empty.

Gap 5 is per background, as pre-registered in `configs/gap5.yaml`: the
background is admissible only if both the gNFW g=1.2 and the VVV controls
pass, and then every template under it is.

Gap 1 Level 2 rows take the Level 1 controls of the same background, mask and
energy.

## Changes 2026-10-26

The Gap 1 export dropped rows whose `result_id` collided. `robustness_scan.py`
never wrote its energy range, so every run defaulted to 1-10 GeV and the
mask B Pohl22 scans (rings_species and rings_free) were dropped in favour of
the 0.5-1 GeV runs, which were labelled 1-10 GeV. It also applied each run's
VVV control to every template, so mask D gNFW rows were marked admissible
when the gNFW control fails there (injected TS 11). Both fixed; the Gap 1
part went from 612 to 1436 rows.

Rebuilt by `analysis/gap1/export_for_gap7.py`, then `analysis/gap5/export_tables.py`.
