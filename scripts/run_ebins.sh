#!/bin/bash
# gap1: counts + exposure over the full conventions range, 0.5-50 GeV.
# reuses the livetime cube (energy independent) so no 40 min wait.
# each conventions bin (0.5-1, 1-3, 3-10, 10-50) is split in two, so the
# fit can report per conventions bin AND the combined range.
cd "$(dirname "$0")"
exec > >(tee -a ebins.log) 2>&1
echo "================ started $(date) ================"

source "$HOME/miniforge3/etc/profile.d/conda.sh" || exit 1
conda activate fermi || exit 1
command -v gtbin >/dev/null || { echo "env not active"; exit 1; }

for f in gap1_gti.fits gap1_ltcube.fits; do
  [ -s "$f" ] && [ -f ".done_$f" ] || { echo "missing or unfinished $f - rerun run_chain.sh first"; exit 1; }
done
echo "reusing gap1_gti.fits and gap1_ltcube.fits"

# ---- 1. four gtbin runs, one per conventions bin, 2 log sub-bins each.
# doing it this way means the conventions edges are exact without needing a
# custom bin-definition file.
EDGES="500:1000 1000:3000 3000:10000 10000:50000"
i=0
for pair in $EDGES; do
  lo=${pair%:*}; hi=${pair#*:}
  out="gap1_cc_part$i.fits"
  if [ -s "$out" ] && [ -f ".done_$out" ]; then
    echo "--- gtbin $lo-$hi MeV: done already"
  else
    echo "--- gtbin $lo-$hi MeV -> $out  ($(date +%H:%M:%S))"
    gtbin evfile=gap1_gti.fits scfile=NONE outfile="$out" \
      algorithm=CCUBE nxpix=400 nypix=400 binsz=0.1 \
      coordsys=GAL xref=0 yref=0 axisrot=0 proj=CAR \
      ebinalg=LOG emin=$lo emax=$hi enumbins=2 || { echo "!!! gtbin failed"; exit 1; }
    touch ".done_$out"
  fi
  i=$((i+1))
done

# ---- 2. stitch them into one cube
python - <<'PY' || exit 1
import numpy as np
from astropy.io import fits
parts = [f"gap1_cc_part{i}.fits" for i in range(4)]
data, lo, hi, first = [], [], [], None
for p in parts:
    hd = fits.open(p)
    if first is None:
        first = hd
    data.append(np.asarray(hd[0].data))
    lo += list(hd["EBOUNDS"].data["E_MIN"]); hi += list(hd["EBOUNDS"].data["E_MAX"])
cube = np.concatenate(data, axis=0)
pri = fits.PrimaryHDU(cube, first[0].header)
pri.header["NAXIS3"] = cube.shape[0]
eb = fits.BinTableHDU.from_columns([
    fits.Column(name="CHANNEL", format="I", array=np.arange(1, len(lo) + 1)),
    fits.Column(name="E_MIN", format="1E", unit="keV", array=np.array(lo, np.float32)),
    fits.Column(name="E_MAX", format="1E", unit="keV", array=np.array(hi, np.float32)),
], name="EBOUNDS", header=first["EBOUNDS"].header)
hdus = [pri, eb] + [e.copy() for e in first[2:]]
fits.HDUList(hdus).writeto("gap1_ccube_full.fits", overwrite=True)
print(f"stitched gap1_ccube_full.fits: {cube.shape}")
for a, b, pl in zip(lo, hi, cube):
    print(f"   {a/1e6:7.3f} - {b/1e6:7.3f} GeV   {pl.sum():>11,.0f}")
if (cube.reshape(len(lo), -1).sum(1) == 0).any():
    print("EMPTY BIN - STOP, DO NOT USE THIS CUBE"); raise SystemExit(1)
PY

# ---- 3. exposure, finely sampled; the fit interpolates to bin centres itself
if [ -s gap1_expcube_full.fits ] && [ -f .done_gap1_expcube_full.fits ]; then
  echo "--- gtexpcube2: done already"
else
  echo "--- gtexpcube2 0.5-50 GeV  ($(date +%H:%M:%S)), a few minutes"
  gtexpcube2 infile=gap1_ltcube.fits cmap=none outfile=gap1_expcube_full.fits \
    irfs=P8R3_SOURCE_V3 evtype=3 \
    nxpix=440 nypix=440 binsz=0.1 coordsys=GAL \
    xref=0 yref=0 axisrot=0 proj=CAR \
    ebinalg=LOG emin=500 emax=50000 enumbins=20 || { echo "!!! gtexpcube2 failed"; exit 1; }
  touch .done_gap1_expcube_full.fits
fi

echo
echo "================ files to send over ================"
for f in gap1_ccube_full.fits gap1_expcube_full.fits; do
  printf '%-26s %s\n' "$f" "$([ -s "$f" ] && du -h "$f" | cut -f1 || echo MISSING)"
done
echo "ALL DONE"
echo "================ finished $(date) ================"
