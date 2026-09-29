#!/bin/bash
# downloads the 4FGL catalog, galactic diffuse model and isotropic spectrum
# from FSSC. resumable - safe to rerun if it dies partway.
cd "$(dirname "$0")"
exec > >(tee -a download_aux.log) 2>&1

echo "=== started $(date) ==="

get() {
  local out="$1"; shift
  if [ -s "$out" ]; then
    echo "[$out] already present ($(du -h "$out" | cut -f1)), skipping"
    return 0
  fi
  for url in "$@"; do
    echo "[$out] trying $url"
    if curl -fL -C - --retry 3 --retry-delay 5 -o "$out" "$url"; then
      echo "[$out] OK  $(du -h "$out" | cut -f1)  from $url"
      return 0
    fi
    echo "[$out] failed on that url"
    rm -f "$out"
  done
  echo "[$out] ALL URLS FAILED"
  return 1
}

B=https://fermi.gsfc.nasa.gov/ssc/data

echo
echo "--- 1/3 isotropic spectrum (tiny) ---"
get iso_P8R3_SOURCE_V3_v1.txt \
  "$B/analysis/software/aux/4fgl/iso_P8R3_SOURCE_V3_v1.txt" \
  "$B/analysis/software/aux/iso_P8R3_SOURCE_V3_v1.txt" \
  "$B/analysis/software/aux/4fgl/iso_P8R3_SOURCE_V3_v01.txt"

echo
echo "--- 2/3 4FGL source catalog (~100 MB) ---"
get gll_psc_catalog.fit \
  "$B/access/lat/14yr_catalog/gll_psc_v35.fit" \
  "$B/access/lat/14yr_catalog/gll_psc_v34.fit" \
  "$B/access/lat/14yr_catalog/gll_psc_v32.fit" \
  "$B/access/lat/12yr_catalog/gll_psc_v31.fit" \
  "$B/access/lat/12yr_catalog/gll_psc_v27.fit"

echo
echo "--- 3/3 galactic diffuse model (~3.7 GB, this is the slow one) ---"
get gll_iem_v07.fits \
  "$B/analysis/software/aux/4fgl/gll_iem_v07.fits" \
  "$B/analysis/software/aux/gll_iem_v07.fits"

echo
echo "=== results ==="
for f in iso_P8R3_SOURCE_V3_v1.txt gll_psc_catalog.fit gll_iem_v07.fits; do
  if [ -s "$f" ]; then printf '%-32s %s\n' "$f" "$(du -h "$f" | cut -f1)"
  else printf '%-32s %s\n' "$f" "MISSING"; fi
done
echo "=== finished $(date) ==="
