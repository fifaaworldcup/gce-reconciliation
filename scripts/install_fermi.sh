#!/bin/bash
# creates the fermi conda env with fermitools. logs to install.log
cd "$(dirname "$0")"
exec > >(tee install.log) 2>&1

echo "=== started $(date) ==="
source "$HOME/miniforge3/etc/profile.d/conda.sh" || { echo "cant source conda"; exit 1; }

echo
echo "arch: $(uname -m)"
echo "rosetta check:"
/usr/bin/pgrep -q oahd && echo "  rosetta daemon running" || echo "  rosetta daemon not detected (may still be fine)"

# already done?
if conda env list | grep -qE '^fermi\s'; then
  echo
  echo "env fermi already exists, checking it"
  conda activate fermi
  if command -v gtselect >/dev/null; then
    echo "gtselect already present at $(command -v gtselect)"
    echo "=== nothing to do ==="
    exit 0
  fi
  echo "env exists but no gtselect, removing and rebuilding"
  conda deactivate
  conda env remove -n fermi -y
fi

echo
echo "=== attempt 1: native $(uname -m) ==="
conda create -n fermi -c conda-forge -c fermi fermitools -y
A1=$?
echo "attempt 1 exit: $A1"

if [ $A1 -ne 0 ]; then
  echo
  echo "=== attempt 1 failed, cleaning up ==="
  conda env remove -n fermi -y 2>/dev/null
  echo
  echo "=== attempt 2: forcing osx-64 (runs under rosetta) ==="
  CONDA_SUBDIR=osx-64 conda create -n fermi -c conda-forge -c fermi fermitools -y
  A2=$?
  echo "attempt 2 exit: $A2"
  if [ $A2 -ne 0 ]; then
    echo
    echo "=== BOTH ATTEMPTS FAILED - stopping, see errors above ==="
    exit 1
  fi
  conda activate fermi
  conda config --env --set subdir osx-64
  echo "pinned env subdir to osx-64"
else
  conda activate fermi
fi

echo
echo "=== verify ==="
echo "CONDA_PREFIX=$CONDA_PREFIX"
for t in gtselect gtmktime gtbin gtltcube gtexpcube2 python; do
  printf '%-12s %s\n' "$t" "$(command -v $t || echo MISSING)"
done

echo
echo "=== versions ==="
conda list -n fermi 2>/dev/null | grep -iE 'fermitools|^python ' 

echo
echo "=== FERMI_DIR / refdata ==="
echo "FERMI_DIR=$FERMI_DIR"
ls "$CONDA_PREFIX/share/fermitools/refdata/fermi" 2>/dev/null | head || echo "refdata dir not where expected"

echo
if command -v gtselect >/dev/null; then
  echo "=== SUCCESS - fermitools installed ==="
else
  echo "=== FAILED - env built but gt tools missing ==="
fi
echo "=== finished $(date) ==="
