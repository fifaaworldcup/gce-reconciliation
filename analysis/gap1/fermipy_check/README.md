# Independent check of the 1–10 GeV window

Section 6.4 reports that the GCE is detected only between 1 and 10 GeV: TS ≈ 178
(1–3 GeV) and ≈ 48 (3–10 GeV), and **TS = 0** in both 0.5–1 GeV and 10–50 GeV.
The two outer bins are bins the Conventions ask for, and the 0.5–1 GeV
non-detection replaced an apparent TS-19 signal that turned out to be an
artefact of the too-narrow PSF approximation (§6.7). A null result that only
one fitter has ever produced should not go in a paper.

This directory runs the same data and the same templates through Fermipy —
the Science Tools' own binned likelihood — and prints the two side by side.

## What would count as agreement

The two fitters are not doing the same thing. Fermipy applies energy
dispersion, fits the spectral shapes of the diffuse components and of nearby
catalogue sources, and uses its own exposure calculation; ours fits
normalisations only, per energy bin, with no dispersion. The TS values will
differ, possibly by tens of percent in the detected bins. What must agree is
the **pattern**:

| | pass | fail |
|---|---|---|
| 0.5–1 GeV | TS < 9 | TS > 9 |
| 1–3 GeV | TS > 25 | TS < 25 |
| 3–10 GeV | TS > 25 | TS < 25 |
| 10–50 GeV | TS < 9 | TS > 9 |

A fail in either outer bin means §6.4 is wrong and the energy window in §6.1
and §6.4 has to be rebuilt. `run_check.py` exits non-zero in that case, on
purpose.

## Running it

Requires the `fermi` conda env from `PREREQS.md`, plus fermipy:

    conda activate fermi
    conda install -c conda-forge fermipy

Then, from the `fermi data` folder:

    cp analysis/gap1/fermipy_check/gap1_fermipy.yaml .
    cp analysis/gap1/fermipy_check/run_check.py .
    mkdir -p fermipy_templates
    cp templates/gap1_bulge_coleman.fits templates/gap1_hestia_G1.1_ang0_1deg.fits fermipy_templates/
    bash scripts/run_fermipy_check.sh

`run_fermipy_check.sh` does the env check, the input check and the install for
you, and logs to `fermipy_check.log`. It reuses `gap1_ltcube.fits`, so it does
not repeat the multi-hour livetime step; it does build its own counts and
exposure cubes, which takes 30–60 minutes on a laptop.

## Why fermipy cannot run in the cloud sandbox

Fermitools and fermipy are conda-forge only, and `conda.anaconda.org` is not
reachable from the sandbox (403 at the egress proxy), so this is a laptop job.
Checked 2026-09-26.

## Notes on the config

`gap1_fermipy.yaml` copies the selection out of `scripts/run_chain.sh` and
`scripts/run_ebins.sh` verbatim — same event class and type, same zenith and
rocking cuts, same time range, same 40°×40° CAR grid at 0.1°, same two log
sub-bins per Conventions bin. Those values exist so a disagreement is a
difference of method rather than of data. Changing one of them to make the
check pass would defeat the point of running it.
