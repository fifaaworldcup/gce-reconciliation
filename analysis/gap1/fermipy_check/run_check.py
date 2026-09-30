"""Cross-check sec 6.4's energy window with Fermipy instead of our fitter.

Sec 6.4 says the GCE is detected only between 1 and 10 GeV, TS 0 in the
0.5-1 and 10-50 GeV bins. Those are Conventions bins and the 0.5-1
non-detection replaced an artefact of the old PSF approximation, so it should
not rest on one fitter.

Same data, same templates, through the Science Tools' binned likelihood.
The TS values will differ - fermipy applies energy dispersion and fits
spectral shapes, we fit normalisations only. What must agree is the pattern.

  PASS  outer bins TS < 9, inner bins TS > 25
  FAIL  either outer bin significant, in which case sec 6.4 is wrong

Run via scripts/run_fermipy_check.sh.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

# what section 6.4 claims, for the comparison table at the end
OURS = {
    "0.5-1":  {"vvv_ts": 0.0,   "hestia_ts_med": 0.0},
    "1-3":    {"vvv_ts": 178.3, "hestia_ts_med": 166.0},
    "3-10":   {"vvv_ts": 48.1,  "hestia_ts_med": 52.0},
    "10-50":  {"vvv_ts": 0.0,   "hestia_ts_med": 0.0},
}
BINS = [(0.5, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 50.0)]
TEMPLATES = {
    "vvv":    "gap1_bulge_coleman.fits",
    "hestia": "gap1_hestia_G1.1_ang0_1deg.fits",
}


def as_spatialmap(src_fits, dst_fits):
    """Strip our bookkeeping keywords and peak-normalise. Shape unchanged."""
    from astropy.io import fits
    d = fits.getdata(src_fits).astype(float)
    h = fits.getheader(src_fits)
    d = np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)
    d[d < 0] = 0.0
    if d.max() <= 0:
        raise SystemExit(f"{src_fits} is empty after cleaning")
    d /= d.max()
    for k in ("NORMED", "PSFCONV", "RHO2CONV", "LONFLIP", "SMOOTHIN"):
        h.pop(k, None)
    fits.PrimaryHDU(d.astype("float32"), h).writeto(dst_fits, overwrite=True)
    return dst_fits


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="gap1_fermipy.yaml")
    ap.add_argument("--templates", default=".",
                    help="directory holding the template FITS files")
    ap.add_argument("--out", default="fermipy_check_results.json")
    ap.add_argument("--free-radius", type=float, default=3.0,
                    help="free the spectral parameters of catalogue sources "
                         "within this many degrees of the centre")
    a = ap.parse_args()

    try:
        from fermipy.gtanalysis import GTAnalysis
    except Exception as exc:
        # fermipy.gtanalysis pulls in the Science Tools (pyLikelihood, GtApp),
        # so this fails for reasons that have nothing to do with fermipy being
        # absent. Print what actually went wrong instead of guessing.
        import traceback
        traceback.print_exc()
        try:
            import fermipy
            have = f"fermipy {fermipy.__version__} imports fine"
        except Exception as e2:
            have = f"fermipy itself will not import: {e2}"
        try:
            import pyLikelihood          # noqa: F401
            st = "pyLikelihood imports fine"
        except Exception as e3:
            st = f"pyLikelihood (Science Tools) will not import: {e3}"
        raise SystemExit(
            f"\ncould not load fermipy.gtanalysis: {type(exc).__name__}: {exc}\n"
            f"  {have}\n  {st}\n"
            "Send the traceback above with the log. The usual cause is fermipy "
            "and fermitools being different versions in the same env, not a "
            "missing package.")

    gta = GTAnalysis(a.config, logging={"verbosity": 3})
    gta.setup()
    # free the diffuse and nearby sources, as a standard analysis would
    gta.free_sources(free=False)
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_sources(distance=a.free_radius, pars="norm")
    gta.optimize()
    gta.fit()
    base_ll = -gta.like()
    print(f"baseline lnL (no GCE component) = {base_ll:.2f}", flush=True)

    out = {"baseline_lnL": float(base_ll), "bins": {}, "templates": {}}
    for name, fn in TEMPLATES.items():
        src = os.path.join(a.templates, fn)
        if not os.path.exists(src):
            print(f"  skipping {name}: {src} not found")
            continue
        mapfile = as_spatialmap(src, f"spatialmap_{name}.fits")
        gta.add_source("gce", {
            "SpatialModel": "SpatialMap", "Spatial_Filename": mapfile,
            "SpectrumType": "PowerLaw", "Index": 2.0,
            "Scale": 3000.0, "Prefactor": 1e-11,
        }, free=True)
        gta.fit()
        full = {"ts_total": float(gta.roi["gce"]["ts"])}
        # per-bin TS as sec 6.4 does it: refit each bin alone
        per_bin = {}
        for lo, hi in BINS:
            key = f"{lo:g}-{hi:g}"
            gta.set_energy_range(np.log10(lo * 1e3), np.log10(hi * 1e3))
            gta.fit()
            per_bin[key] = float(gta.roi["gce"]["ts"])
            print(f"  {name} {key} GeV: TS = {per_bin[key]:.1f}", flush=True)
        gta.set_energy_range(None, None)
        full["ts_per_bin"] = per_bin
        out["templates"][name] = full
        gta.delete_source("gce")

    print("\n================ section 6.4 cross-check ================")
    print(f"{'bin (GeV)':>10} {'ours VVV':>10} {'fermipy VVV':>12} "
          f"{'ours HESTIA':>12} {'fermipy HESTIA':>15}")
    verdict = True
    for lo, hi in BINS:
        key = f"{lo:g}-{hi:g}"
        fv = out["templates"].get("vvv", {}).get("ts_per_bin", {}).get(key)
        fh = out["templates"].get("hestia", {}).get("ts_per_bin", {}).get(key)
        print(f"{key:>10} {OURS[key]['vvv_ts']:10.1f} "
              f"{(fv if fv is not None else float('nan')):12.1f} "
              f"{OURS[key]['hestia_ts_med']:12.1f} "
              f"{(fh if fh is not None else float('nan')):15.1f}")
        if fv is None:
            continue
        inner = key in ("1-3", "3-10")
        if inner and fv < 25:
            verdict = False
        if not inner and fv > 9:
            verdict = False
    out["agrees_with_section_6_4"] = bool(verdict)
    print("\nVERDICT: " + ("fermipy reproduces the 1-10 GeV window - section "
                           "6.4 stands"
                           if verdict else
                           "fermipy DISAGREES - section 6.4 must be rewritten "
                           "before the paper goes anywhere"))
    with open(a.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"written {a.out}")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
