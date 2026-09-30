#!/bin/bash
# gap1 fermitools chain. resumable, and it checks the photon files actually
# cover the band we're asking for BEFORE spending 40 minutes on a livetime cube.
cd "$(dirname "$0")"
exec > >(tee -a chain.log) 2>&1

echo "================ started $(date) ================"

source "$HOME/miniforge3/etc/profile.d/conda.sh" || exit 1
conda activate fermi || exit 1
command -v gtselect >/dev/null || { echo "gtselect missing, env not active"; exit 1; }
echo "env: $CONDA_PREFIX"

SC=$(ls *SC00.fits 2>/dev/null | head -1)
[ -n "$SC" ] || { echo "no spacecraft file found"; exit 1; }
echo "spacecraft: $SC"

# ---- quarantine any photon file from the old 5 GeV download ---------------
mkdir -p old_5gev
python - <<'PY'
import glob, os, shutil
from astropy.io import fits
moved = []
for f in sorted(glob.glob("*_PH*.fits")):
    try:
        with fits.open(f) as hd:
            e = hd[1].data["ENERGY"]
            lo = float(e.min())
    except Exception as ex:
        print(f"  {f}: unreadable ({ex}), leaving it alone")
        continue
    if lo >= 4000.0:
        shutil.move(f, os.path.join("old_5gev", f))
        moved.append(f)
print(f"  quarantined {len(moved)} file(s) from the 5 GeV download")
PY

ls *_PH*.fits > events.txt 2>/dev/null
N=$(wc -l < events.txt | tr -d ' ')
echo "event files: $N"
[ "$N" -gt 0 ] || { echo "no photon files left - download the new data first"; exit 1; }

# ---- fail fast if the new data still doesn't reach low enough -------------
echo
echo "--- checking energy coverage before doing any real work ---"
python - <<'PY' || exit 1
import sys
import numpy as np
from astropy.io import fits

NEED_LO = 600.0     # MeV. we fit from 1 GeV, and want margin for dispersion
lo_all, hi_all, n = [], [], 0
with open("events.txt") as f:
    files = [x.strip() for x in f if x.strip()]
for fn in files:
    with fits.open(fn) as hd:
        e = np.asarray(hd[1].data["ENERGY"], float)
        lo_all.append(e.min()); hi_all.append(e.max()); n += len(e)
        h = hd[1].header
        cut = None
        for k in h:
            if k.startswith("DSTYP") and h[k].strip() == "ENERGY":
                cut = h.get("DSVAL" + k[5:])
    print(f"  {fn}: {len(e):>10,} events, {e.min():9.1f} - {e.max():9.1f} MeV"
          f"   cut={cut}")
lo, hi = min(lo_all), max(hi_all)
print(f"\n  TOTAL {n:,} events, {lo:.1f} - {hi:.1f} MeV")
if lo > NEED_LO:
    print(f"\n  STOP: lowest event is {lo:.0f} MeV but the fit starts at "
          f"1000 MeV and needs data from about {NEED_LO:.0f} MeV.")
    print("  The download's energy range is still too high. Re-query FSSC")
    print("  with energy 500 to 500000 MeV and rerun this script.")
    sys.exit(1)
print("  energy coverage OK")
PY

# ---- if the photon set changed, derived products are stale ----------------
NEWSUM=$(md5 -q events.txt 2>/dev/null || md5sum events.txt | cut -d' ' -f1)
OLDSUM=$(cat .events.md5 2>/dev/null)
if [ "$NEWSUM" != "$OLDSUM" ]; then
  echo
  echo "--- photon file set changed, clearing products built from the old set"
  for f in gap1_sel.fits gap1_gti.fits gap1_ccube.fits gap1_ltcube.fits \
           gap1_expcube.fits; do
    [ -e "$f" ] && { mv -f "$f" "old_5gev/$f" 2>/dev/null && echo "    moved $f"; }
    rm -f ".done_$f"
  done
  echo "$NEWSUM" > .events.md5
fi
echo

# a step is complete only when its marker exists. checking the output file
# itself is not enough: a run killed mid-write leaves a truncated file that
# looks finished and gets silently skipped. that nearly happened already.
step() {
  local out="$1" label="$2"; shift 2
  if [ -f ".done_$out" ] && [ -s "$out" ]; then
    echo "--- $label: $out already done, skipping"
    return 0
  fi
  rm -f "$out"
  echo "--- $label -> $out   ($(date +%H:%M:%S))"
  "$@"
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$out" ]; then
    echo "!!! $label FAILED (exit $rc)"
    return 1
  fi
  touch ".done_$out"
  echo "    done $(date +%H:%M:%S), $(du -h "$out" | cut -f1)"
  return 0
}

step gap1_sel.fits "1/6 gtselect" \
  gtselect infile=@events.txt outfile=gap1_sel.fits \
    ra=266.404988 dec=-28.936178 rad=30 \
    tmin=239557417 tmax=760320000 emin=500 emax=50000 \
    zmax=90 evclass=128 evtype=3 || exit 1

step gap1_gti.fits "2/6 gtmktime" \
  gtmktime scfile="$SC" \
    filter="DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
    roicut=no evfile=gap1_sel.fits outfile=gap1_gti.fits || exit 1

step gap1_ccube.fits "3/6 gtbin CCUBE" \
  gtbin evfile=gap1_gti.fits scfile=NONE outfile=gap1_ccube.fits \
    algorithm=CCUBE nxpix=400 nypix=400 binsz=0.1 \
    coordsys=GAL xref=0 yref=0 axisrot=0 proj=CAR \
    ebinalg=LOG emin=1000 emax=10000 enumbins=8 || exit 1

echo "--- 4/6 gtltcube: the slow one, 20-40 min. it prints a dot now and then."
step gap1_ltcube.fits "4/6 gtltcube" \
  gtltcube evfile=gap1_gti.fits scfile="$SC" \
    outfile=gap1_ltcube.fits dcostheta=0.025 binsz=1 zmax=90 || exit 1

step gap1_expcube.fits "5/6 gtexpcube2" \
  gtexpcube2 infile=gap1_ltcube.fits cmap=none outfile=gap1_expcube.fits \
    irfs=P8R3_SOURCE_V3 evtype=3 \
    nxpix=440 nypix=440 binsz=0.1 coordsys=GAL \
    xref=0 yref=0 axisrot=0 proj=CAR \
    emin=1000 emax=10000 enumbins=8 || exit 1

step gap1_iem_roi.fits "6/6 reduce_iem" \
  python reduce_iem.py gll_iem_v07.fits gap1_iem_roi.fits || exit 1

echo
echo "================ sanity check ================"
python - <<'PY'
from astropy.io import fits
import numpy as np
c = fits.getdata("gap1_ccube.fits")
eb = fits.getdata("gap1_ccube.fits", "EBOUNDS")
print(f"counts cube {c.shape}, {np.sum(c):,.0f} photons")
bad = 0
for i in range(c.shape[0]):
    lo, hi = eb["E_MIN"][i] / 1e6, eb["E_MAX"][i] / 1e6
    s = float(np.sum(c[i]))
    flag = "  <-- EMPTY, SOMETHING IS WRONG" if s == 0 else ""
    if s == 0:
        bad += 1
    print(f"  bin {i}: {lo:6.3f} - {hi:6.3f} GeV   {s:>12,.0f}{flag}")
print("ALL BINS POPULATED" if bad == 0 else f"{bad} EMPTY BINS - STOP, DO NOT USE THIS CUBE")
for f in ("gap1_expcube.fits", "gap1_iem_roi.fits"):
    d = fits.getdata(f)
    print(f"{f}: {d.shape} sum={np.nansum(d):.6g}")
PY

echo
echo "================ files to send over ================"
for f in gap1_ccube.fits gap1_expcube.fits gap1_iem_roi.fits gll_psc_catalog.fit; do
  if [ -s "$f" ]; then printf '%-28s %s\n' "$f" "$(du -h "$f" | cut -f1)"
  else printf '%-28s %s\n' "$f" "MISSING"; fi
done
echo "================ finished $(date) ================"
