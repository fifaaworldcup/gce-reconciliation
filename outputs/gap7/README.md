# covariance_input.csv

Gap 1 rows (`outputs/gap1/gap1_results_for_gap7.csv`, unchanged) plus the Gap 5
Level 2 rows. Same columns as the Gap 1 file, with three added at the end:

- `arm`: `neutral` (one common source list, all 40 templates) or `joint`
  (each template runs its own source search). Empty for Gap 1 rows.
- `iem_set`: `reference`, `reselected`, `original`, or `reselected+original` (XXVII).
- `fit_gap_vs_best_iem`: best AIC under that background minus best AIC under
  any background (mask B, 1-10 GeV). gll_iem_v07 is 0.

`iem_family` for Gap 5 rows: 1 = gll_iem_v07, 2 = Pohl22 rings_species,
3 = Cholis et al. (2022). The Gap 5 Pohl22 rows reproduce the Gap 1 Level 2
rows (same lnL to 0.01), so filter on `gap` to avoid counting them twice.

Use rows with `status` starting `OK`. Rebuilt by `analysis/gap5/export_tables.py`.
