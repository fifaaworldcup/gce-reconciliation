#!/bin/bash
# Fermipy cross-check of sec 6.4's 1-10 GeV window. Run from the fermi data
# folder. Reuses gap1_ltcube.fits.
cd "$(dirname "$0")"
exec > >(tee -a fermipy_check.log) 2>&1
echo "================ started $(date) ================"

source "$HOME/miniforge3/etc/profile.d/conda.sh" || exit 1
conda activate fermi || exit 1

# ---- 1. environment -------------------------------------------------------
if ! python -c "from fermipy.gtanalysis import GTAnalysis" 2>/dev/null; then
  echo "fermipy.gtanalysis will not load. Installing/repairing fermipy and"
  echo "fermitools together so their versions match (one time, several"
  echo "minutes)."
  conda install -y -c conda-forge fermipy fermitools || {
    echo "!!! fermipy install failed - stop and tell claude what it printed"
    exit 1; }
fi
python - <<'EOF' || exit 1
import fermipy
print("fermipy", fermipy.__version__)
from fermipy.gtanalysis import GTAnalysis          # the import that matters
import pyLikelihood
print("Science Tools import OK")
EOF

# ---- 2. inputs it needs ---------------------------------------------------
SC=$(ls *_SC00.fits 2>/dev/null | head -1)
[ -n "$SC" ] || { echo "no spacecraft file here"; exit 1; }
ln -sf "$SC" SC.fits
ls *_PH*.fits > events.txt 2>/dev/null
[ -s events.txt ] || { echo "no photon files here"; exit 1; }
echo "events.txt: $(wc -l < events.txt) file(s); SC = $SC"

for f in gap1_ltcube.fits gll_iem_v07.fits iso_P8R3_SOURCE_V3_v1.txt; do
  [ -s "$f" ] || { echo "missing $f - it should be in this folder"; exit 1; }
done

# fermipy ships its own 4FGL-DR4 copy, so nothing to do here

# ---- 3. the templates -----------------------------------------------------
mkdir -p fermipy_templates
for t in gap1_bulge_coleman.fits gap1_hestia_G1.1_ang0_1deg.fits; do
  if [ ! -s "fermipy_templates/$t" ]; then
    echo "!!! missing fermipy_templates/$t"
    echo "    copy it from the repo's templates/ folder into that directory"
    exit 1
  fi
done

# ---- 4. run ---------------------------------------------------------------
cp -f gap1_fermipy.yaml fermipy_run.yaml 2>/dev/null || {
  echo "!!! gap1_fermipy.yaml not here - copy it from"
  echo "    analysis/gap1/fermipy_check/ in the repo"; exit 1; }

echo "--- fermipy setup + fit. First run builds its own counts and exposure"
echo "    cubes, so expect 30-60 minutes. It reuses gap1_ltcube.fits."
python run_check.py --config fermipy_run.yaml \
  --templates fermipy_templates \
  --out fermipy_check_results.json
rc=$?
echo "================ finished $(date), exit $rc ================"
echo
if [ $rc -eq 0 ]; then
  echo "Section 6.4 survives an independent fitter. Send claude"
  echo "fermipy_check_results.json."
else
  echo "Either the run failed or fermipy disagrees. Either way send claude"
  echo "fermipy_check.log - do not edit section 6.4 yourself."
fi
exit $rc
