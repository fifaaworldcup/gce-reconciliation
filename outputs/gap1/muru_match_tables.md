# Matched-convention HESTIA shape measurements

Produced by `analysis/gap1/muru_match_run.py`. Conventions and what
they mean: `docs/muru_conventions.md`. Full records in
`outputs/gap1/muru_match.json`, not committed.

## Axis ratio at the 50% isophote, Muru's SVD estimator, 3 deg cone

| galaxy | angle | column | column squared | J-factor | q(J) - q(col^2) |
|---|---|---|---|---|---|
| G1.1 | 0 | 0.730 | 0.793 | 0.743 | -0.050 |
| G1.1 | 30 | 0.684 | 0.717 | 0.694 | -0.023 |
| G1.1 | 60 | 0.617 | 0.612 | 0.604 | -0.008 |
| G1.1 | 90 | 0.568 | 0.587 | 0.567 | -0.020 |
| G1.2 | 0 | 0.687 | 0.755 | 0.717 | -0.037 |
| G1.2 | 30 | 0.727 | 0.787 | 0.758 | -0.029 |
| G1.2 | 60 | 0.767 | 0.831 | 0.786 | -0.045 |
| G1.2 | 90 | 0.812 | 0.878 | 0.828 | -0.050 |
| G2.1 | 0 | 0.682 | 0.738 | 0.692 | -0.046 |
| G2.1 | 30 | 0.821 | 0.831 | 0.823 | -0.008 |
| G2.1 | 60 | 0.873 | 0.869 | 0.873 | +0.004 |
| G2.1 | 90 | 0.749 | 0.749 | 0.735 | -0.015 |
| G2.2 | 0 | 0.651 | 0.735 | 0.683 | -0.052 |
| G2.2 | 30 | 0.709 | 0.800 | 0.751 | -0.049 |
| G2.2 | 60 | 0.797 | 0.913 | 0.854 | -0.059 |
| G2.2 | 90 | 0.935 | 0.974 | 0.979 | +0.005 |
| G3.1 | 0 | 0.655 | 0.720 | 0.667 | -0.054 |
| G3.1 | 30 | 0.736 | 0.798 | 0.750 | -0.048 |
| G3.1 | 60 | 0.805 | 0.846 | 0.820 | -0.026 |
| G3.1 | 90 | 0.839 | 0.873 | 0.863 | -0.010 |
| G3.2 | 0 | 0.863 | 0.879 | 0.846 | -0.032 |
| G3.2 | 30 | 0.924 | 0.933 | 0.923 | -0.011 |
| G3.2 | 60 | 0.780 | 0.819 | 0.765 | -0.054 |
| G3.2 | 90 | 0.641 | 0.700 | 0.649 | -0.052 |

Mean -0.032, sd 0.020, range -0.059 to +0.005, same sign in 22 of 24.

## Disc frame from the AHF profile

| galaxy | profile row | particles | b/a | c/a | minor axis vs L |
|---|---|---|---|---|---|
| G1.1 | 10.28 kpc/h = 15.19 kpc | 4,120,545 | 0.977 | 0.677 | 1.4 deg |
| G1.2 | 10.45 kpc/h = 15.43 kpc | 5,246,155 | 0.990 | 0.789 | 0.7 deg |
| G2.1 | 10.68 kpc/h = 15.78 kpc | 3,452,716 | 0.994 | 0.750 | 0.1 deg |
| G2.2 | 9.03 kpc/h = 13.34 kpc | 3,621,096 | 0.994 | 0.727 | 0.5 deg |
| G3.1 | 10.62 kpc/h = 15.69 kpc | 3,191,090 | 0.994 | 0.717 | 3.1 deg |
| G3.2 | 10.29 kpc/h = 15.19 kpc | 2,222,072 | 0.988 | 0.708 | 0.3 deg |

b/a of 0.977-0.994 is why the in-plane azimuth reference is
undefined: there is no distinguishable major axis within the plane.

## Axis ratio against isophote, on the maps Muru sent

| level, fraction of peak | mean q | min | max | equivalent radius | level set at ROI edge |
|---|---|---|---|---|---|
| 0.70 | 0.682 | 0.497 | 0.971 | 1.6 deg | 0/24 |
| 0.50 | 0.663 | 0.423 | 0.850 | 2.6 deg | 0/24 |
| 0.30 | 0.634 | 0.467 | 0.760 | 4.2 deg | 0/24 |
| 0.20 | 0.613 | 0.518 | 0.745 | 5.7 deg | 0/24 |
| 0.10 | 0.614 | 0.534 | 0.700 | 9.1 deg | 0/24 |
| 0.05 | 0.625 | 0.542 | 0.707 | 13.5 deg | 11/24 |
| 0.02 | 0.796 | 0.667 | 0.984 | 19.1 deg | 24/24 |

## Spearman rho against dAIC vs VVV, Mask B + PS, 1-10 GeV, n = 24

| q used | rho | p |
|---|---|---|
| Muru's own maps, squared column, q10 | +0.463 | 0.023 |
| current squared column, 5% contour | +0.428 | 0.037 |
| matched J-factor, cone 3, q10 | +0.252 | 0.23 |
| matched squared column, cone 1, q50 | +0.252 | 0.23 |
| Muru's own maps, squared column, q50 | -0.251 | 0.24 |
| matched J-factor, cone 1, q10 | +0.248 | 0.24 |
| Muru's own maps, squared column, q30 | -0.227 | 0.29 |
| matched squared column, cone 3, q5 | +0.223 | 0.3 |
| Muru's own maps, squared column, q5 | +0.217 | 0.31 |
| Muru's own maps, squared column, q70 | -0.202 | 0.34 |
| matched J-factor, cone 1, q70 | +0.195 | 0.36 |
| matched J-factor, cone 3, q20 | +0.187 | 0.38 |
| matched squared column, cone 1, q5 | +0.164 | 0.44 |
| matched squared column, cone 1, q70 | +0.140 | 0.51 |
| matched J-factor, cone 1, q20 | +0.139 | 0.52 |
| matched squared column, cone 3, q10 | +0.135 | 0.53 |
| Muru's own maps, squared column, q20 | -0.130 | 0.54 |
| matched squared column, cone 1, q10 | +0.129 | 0.55 |
| matched J-factor, cone 3, q70 | +0.120 | 0.58 |
| matched J-factor, cone 3, q30 | +0.119 | 0.58 |
| matched squared column, cone 3, q70 | +0.118 | 0.58 |
| matched squared column, cone 3, q50 | +0.095 | 0.66 |
| matched J-factor, cone 1, q5 | +0.087 | 0.69 |
| matched J-factor, cone 1, q30 | +0.084 | 0.7 |
| matched J-factor, cone 3, q5 | -0.062 | 0.77 |
| current J-factor, 5% contour | +0.047 | 0.83 |
| matched squared column, cone 3, q20 | +0.046 | 0.83 |
| matched J-factor, cone 3, q50 | +0.039 | 0.86 |
| matched squared column, cone 3, q30 | +0.029 | 0.89 |
| matched squared column, cone 1, q30 | +0.024 | 0.91 |
| matched J-factor, cone 1, q50 | +0.022 | 0.92 |
| matched squared column, cone 1, q20 | -0.003 | 0.99 |
