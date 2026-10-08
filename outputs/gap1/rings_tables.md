# Section 6.1-6.5 on the Pohl22 ring background

Produced by `analysis/gap1/robustness_scan.py --rings` and
`analysis/gap1/fit_rings.py`. Full records under `outputs/gap1/rsp_*`,
`outputs/gap1/rings_*`, not committed. Every fit converges: summed lnL gap
bound 6e-7 over the eight bins, so every TS below is good to ~1e-6.

## Table R1. The freedom ladder, Mask B + PS, 1-10 GeV, observed TS

| diffuse model | components | gNFW 1.2 | HESTIA G1.1 0deg | VVV |
|---|---|---|---|---|
| `gll_iem_v07` | 1 | 120.7 | 248.0 | 203.1 |
| `rings_species` | 5 | 526.2 | 622.6 | 2151.7 |
| `rings_free` | 16 | 125.9 | 75.9 | 794.0 |
| `rings_free_NB` | 17 | 21.8 | 12.1 | 570.5 |
| `rings_free_NB_BB` | 18 | 0.0 | -0.0 | 124.2 |

## Table R2. Admissibility is per template

Injected TS and boosted gain for a signal added at the flux the
`gll_iem_v07` fit finds for that same template. Threshold 25 on both, set
in `configs/diffuse_freedom.yaml` before any run.

| diffuse model | VVV | HESTIA | gNFW 1.2 |
|---|---|---|---|
| `gll_iem_v07` | 184.3 / +585 **yes** | 240.0 / +684 **yes** | 109.3 / +320 **yes** |
| `rings_species` | 255.2 / +1472 **yes** | 219.3 / +746 **yes** | 134.7 / +478 **yes** |
| `rings_free` | 82.9 / +505 **yes** | 33.8 / +137 **yes** | 23.6 / +135 no |
| `rings_free_NB` | 37.3 / +343 **yes** | 15.8 / +37 no | 21.4 / +26 no |
| `rings_free_NB_BB` | 8.9 / +49 no | 14.3 / -0 no | 10.2 / -0 no |

`rings_free` recovers an injected stellar bulge at TS 83 and an injected
HESTIA halo at 34, but an injected gNFW at only 23.6. The same
configuration is admissible for two template classes and not the third,
which is why the control is now run per template.

## Table R3. Per-template controls at one level, Mask B + PS, 1-10 GeV

| injected template | rings_species | rings_free | Mask D, rings_free |
|---|---|---|---|
| `vvv` | 255.2 **yes** | 82.9 **yes** | 14.9 no |
| `f98` | 223.1 **yes** | 47.6 **yes** | 4.3 no |
| `hestia_G1.1_ang0` | 219.3 **yes** | 33.8 **yes** | 9.7 no |
| `hestia_G3.2_ang90` | 139.9 **yes** | 25.5 **yes** | 2.9 no |
| `gnfw_g1.2` | 134.7 **yes** | 23.6 no | 2.4 no |
| `gnfw_g1.4` | - | 36.4 **yes** | 6.4 no |

Nothing is admissible under Mask D once the background is ring-decomposed,
and Mask A keeps 8.7% of the ROI and has no power at all. Mask B + PS is the
only usable configuration.

## Table R4. Section 6.2 redone at the admissible level

`rings_species`, Mask B + PS, 1-10 GeV. Every template in this table passes
its own injection control.

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 2151.7 | -0.0 |
| `f98` | 1761.2 | +390.5 |
| `hestia_G3.2_ang90` | 1263.2 | +888.5 |
| `hestia_G2.2_ang90` | 1084.3 | +1067.4 |
| `hestia_G2.2_ang0` | 1056.2 | +1095.5 |
| `hestia_G2.2_ang60` | 1033.4 | +1118.2 |
| `hestia_G2.2_ang30` | 977.1 | +1174.6 |
| `hestia_G1.1_ang90` | 947.1 | +1204.6 |
| ... | | |
| `gnfw_g0.9_th1` | 0.0 | +2151.7 |

## Table R5. Section 6.1, the axis-ratio correlation, on each background

| background | Table 1 q (published recipe) | per-angle 5% | per-angle 50% |
|---|---|---|---|
| `gll_iem_v07`, Mask B + PS (as published) | +0.656 (p 0.00051) | +0.217 (p 0.31) | -0.251 (p 0.24) |
| `rings_species`, Mask B + PS | -0.190 (p 0.37) | +0.253 (p 0.23) | +0.363 (p 0.081) |
| `rings_free`, Mask B + PS | -0.056 (p 0.79) | +0.400 (p 0.053) | +0.383 (p 0.064) |
| `rings_species`, Mask D + PS | +0.063 (p 0.77) | +0.457 (p 0.025) | +0.157 (p 0.47) |

## Table R6. Section 6.4, energy dependence, rings_species, Mask B + PS

| band GeV | VVV | F98 | best HESTIA | gNFW 1.2 |
|---|---|---|---|---|
| 0.5-1 | 0.0 | 0.0 | hestia_G2.2_ang0 0.0 | 0.0 |
| 1-3 | 1628.9 | 1279.4 | hestia_G3.2_ang90 968.1 | 368.7 |
| 3-10 | 736.3 | 651.5 | hestia_G3.2_ang90 380.8 | 152.4 |
| 10-50 | 20.9 | 14.8 | hestia_G2.2_ang60 0.7 | 0.0 |

The 0.5-1 GeV band detects nothing at all on this background, for any
template, and no injection control was run for it; it is not a null, it is
an empty cell.

## Table R7. Section 6.5, triaxial ordering, rings_species, Mask B + PS

| template | TS | dAIC vs VVV |
|---|---|---|
| `vvv` | 2151.7 | -0.0 |
| `triaxial_vert0.6_plane0.8_flip` | 764.3 | +1387.4 |
| `triaxial_vert0.6_plane0.8_std` | 692.5 | +1459.2 |
| `triaxial_vert0.8_plane0.6_flip` | 643.3 | +1508.3 |
| `hestia_G1.1_ang0` | 622.6 | +1529.1 |
| `gnfw_g1.2` | 526.2 | +1625.5 |
| `triaxial_vert0.8_plane0.6_std` | 510.0 | +1641.7 |
| `hestia_G2.1_ang90` | 329.9 | +1821.8 |
| `gnfw_g1.1` | 318.6 | +1833.1 |

Vertical flattening still beats in-plane flattening, 764 against 643, the
same ordering as on the single-template background. Both are a long way
behind the stellar bulge.
