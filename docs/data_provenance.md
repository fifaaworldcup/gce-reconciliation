# Data provenance

What the large inputs are, where they come from, and whether they can be
regenerated. None of them are committed.

Put them where the scripts expect them: a working folder for the chain
outputs, `iem/` for the libraries. Both are gitignored.

| what | where it comes from | size | status |
|---|---|---|---|
| Photon and spacecraft files, counts / exposure / livetime cubes | FSSC query, then `scripts/run_chain.sh` and `scripts/run_ebins.sh`. Selection in `analysis/gap1/GET_DATA.md`. | ~6 GB | reproducible |
| `gll_iem_v07.fits`, `iso_P8R3_SOURCE_V3_v1.txt` | Fermi diffuse models, `scripts/download_aux.sh` | ~4 GB | reproducible |
| 4FGL-DR4 catalogue `gll_psc_*.fit` | FSSC | 8 MB | reproducible |
| **IEM family 2** — Cholis, Zhong, McDermott & Surdutovich 2022, 80 models | Zenodo record 6423495, `GCE_TEMPLATES_ZENODO_FILES_v3.zip`. Paper arXiv:2112.09706. CC-BY. | 2.0 GB | download by hand |
| **IEM family 3, Pohl22** — gas rings and dust residuals | Zenodo record 6276721. Contains `HI_pohl_Ts_bestTexc_ring_1..4`, `H2_pohl_ring_1..4`, dust residual maps, boxy bulge, Fermi Bubbles. The committed `templates/bubbles/*` came from here. | — | download by hand |
| **IEM family 3, Galp21** — SA0 / SA50 / SA100 | Porter, Jóhannesson & Moskalenko 2022, arXiv:2112.12745. **No precomputed maps are published anywhere**: GALPROP ships the GALDEF configuration files, not the skymaps, so these must be regenerated with GALPROP v57 or obtained from the authors. Checked 2026-09-28. | — | not public |
| HESTIA J-factor projections | Muru et al. 2025, arXiv:2508.06314, from the authors. Note these are normalised to a common total and carry no absolute J-factor — see Gap 1 §6.10. | — | on request |

Zenodo is blocked by some institutional egress proxies; both archives above
had to be fetched from an ordinary browser.

## HESTIA halo particle files, 2026-10-05

Six Arrow files from `cloud.aip.de/index.php/s/3k766wSLtbzxgi9`, ~1 GB each,
written by Julia `Arrow.jl`. Columns: particle ID, `ptype` (dm / bh / gas /
stars, a Julia Symbol), coordinates and velocities, masses, stellar formation
time. Coordinates are absolute box coordinates in Mpc/h with h = 0.677; AHF
halo centres, the softening (220 pc) and the DM particle mass (1.5 × 10⁵ M☉)
are in the readme shipped with them.

**The file names are swapped within each pair.** The file named `G11` holds the
halo at G1.2's AHF centre, `G12` holds G1.1's, and G31/G32 likewise — 850+ kpc
apart, so it is not a near miss. Identify each file by matching its dark-matter
centroid to the AHF centres, never by its name: `scripts/hestia/whichcentre.py`
prints the match, and `scripts/hestia/reduce.py` does it before writing
anything. G21/G22 happen to be correct, which is the trap.

`reduce.py` writes DM and star positions in physical kpc about the AHF centre
and masses in 10¹⁰ M☉, streaming one column at a time as float32 — the full
17 M particles as float64 exhausts memory.

Conventions for everything downstream of these files, read off Muru's Julia
rather than assumed: `docs/muru_conventions.md`.
