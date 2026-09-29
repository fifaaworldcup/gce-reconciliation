#!/bin/bash
# gap1: the real LAT PSF for our exact data, from the livetime cube.
# ~1 minute. utils/psf.py gets checked against this, and if it's off the fits
# get rerun with the real one.
cd "$(dirname "$0")"
exec > >(tee -a gtpsf.log) 2>&1
source "$HOME/miniforge3/etc/profile.d/conda.sh" || exit 1
conda activate fermi || exit 1
[ -s gap1_ltcube.fits ] || { echo "no gap1_ltcube.fits"; exit 1; }
gtpsf expcube=gap1_ltcube.fits outfile=gap1_psf.fits irfs=P8R3_SOURCE_V3 evtype=3 \
  ra=266.404988 dec=-28.936178 emin=500 emax=50000 nenergies=30 \
  thetamax=10 ntheta=500 || { echo "!!! gtpsf failed"; exit 1; }
ls -la gap1_psf.fits && echo "ALL DONE"
