#!/bin/bash
# builds the livetime + exposure cubes from the data already on disk.
#
# the livetime cube depends only on the spacecraft file and the good-time
# intervals, neither of which depends on energy, so this is very likely the
# same file the new download will need. claude checks the GTIs match before
# reusing it. worst case it gets rebuilt and we lost nothing but idle time.
#
# does NOT touch the photon files or rerun gtselect/gtmktime/gtbin.
cd "$(dirname "$0")"
exec > >(tee -a dry_run.log) 2>&1

echo "================ started $(date) ================"
source "$HOME/miniforge3/etc/profile.d/conda.sh" || exit 1
conda activate fermi || exit 1
command -v gtltcube >/dev/null || { echo "env not active"; exit 1; }

SC=$(ls *SC00.fits 2>/dev/null | head -1)
for f in gap1_gti.fits gap1_ccube.fits "$SC"; do
  [ -s "$f" ] || { echo "missing $f - run the main chain first"; exit 1; }
done
echo "using existing gap1_gti.fits and gap1_ccube.fits (5-10 GeV data)"
echo "spacecraft: $SC"
echo

step() {
  local out="$1" label="$2"; shift 2
  if [ -f ".done_$out" ] && [ -s "$out" ]; then
    echo "--- $label: already done, skipping"; return 0
  fi
  rm -f "$out"
  echo "--- $label -> $out   ($(date +%H:%M:%S))"
  "$@"
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$out" ]; then echo "!!! $label FAILED ($rc)"; return 1; fi
  touch ".done_$out"
  echo "    done $(date +%H:%M:%S), $(du -h "$out" | cut -f1)"
}

echo "--- gtltcube: 20-40 min. prints a dot occasionally. this is the file"
echo "    we're hoping to reuse for the new data."
step gap1_ltcube.fits "1/3 gtltcube" \
  gtltcube evfile=gap1_gti.fits scfile="$SC" \
    outfile=gap1_ltcube.fits dcostheta=0.025 binsz=1 zmax=90 || exit 1

step gap1_expcube.fits "2/3 gtexpcube2" \
  gtexpcube2 infile=gap1_ltcube.fits cmap=none outfile=gap1_expcube.fits \
    irfs=P8R3_SOURCE_V3 evtype=3 \
    nxpix=440 nypix=440 binsz=0.1 coordsys=GAL \
    xref=0 yref=0 axisrot=0 proj=CAR \
    emin=1000 emax=10000 enumbins=8 || exit 1

step gap1_iem_roi.fits "3/3 reduce_iem" \
  python reduce_iem.py gll_iem_v07.fits gap1_iem_roi.fits || exit 1

echo
echo "================ done - tell claude ================"
for f in gap1_ccube.fits gap1_expcube.fits gap1_iem_roi.fits gap1_ltcube.fits; do
  printf '%-26s %s\n' "$f" "$([ -s "$f" ] && du -h "$f" | cut -f1 || echo MISSING)"
done
echo "================ finished $(date) ================"
