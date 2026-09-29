"""
Gap 1 baseline fit - Afeefa's steps 4 and 5, once with HESTIA G1.1 and once
with a second GCE template, on the same data and background model.

!!! NEVER EXECUTED. Fermitools cannot be installed in our sandbox (conda only,
and the proxy blocks conda.anaconda.org and fermi.gsfc.nasa.gov), so this has
been written and reviewed but not run. Expect to debug it on first contact.

Usage:
    python run_fit.py --config ../../configs/config_gap1_maskA_iem1.yaml \
        --gce hestia=../../templates/gap1_hestia_G1.1_ang0_1deg.fits \
        --gce gnfw=../../templates/gnfw_gamma1.2.fits

Templates are passed as name=path. Each is fitted as the GCE in turn against
the SAME stage-1 background fit, which is what makes the delta-lnL meaningful.

IMPORTANT: the templates handed to fermipy must NOT be PSF-convolved.
gtsrcmaps convolves a SpatialMap with the PSF and exposure itself; feeding it a
pre-convolved map blurs twice and makes every template artificially round. Our
FITS headers carry PSFCONV=False for exactly this reason - the script refuses to
run if that keyword says otherwise.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hestia_to_template as H   # for the shape statistics


# ---------------------------------------------------------------- guards
def check_not_preconvolved(path):
    hdr = fits.getheader(path)
    if bool(hdr.get("PSFCONV", False)):
        raise SystemExit(
            f"{path} is already PSF-convolved. gtsrcmaps will convolve again, "
            "which double-blurs the template and biases it round. Use the "
            "intrinsic map instead.")
    if not bool(hdr.get("NORMED", False)):
        print(f"  warning: {os.path.basename(path)} has no NORMED keyword - "
              "the fitted normalisation will not be comparable across templates")


# ------------------------------------------- shape of an INPUT template
def shape_of_input(path, wide_path=None, max_radius_deg=15.0):
    """Shape of the template as handed to the fit, measured BOTH ways.

    Afeefa asked for the ROI-cropped q and c4 alongside the uncropped ones, so
    the size of the crop artifact is visible instead of assumed. The difference
    between the two columns IS the artifact: on the HESTIA maps it was the whole
    apparent boxiness (+0.02 cropped, ~0.00 uncropped).

    'uncropped' needs a map wider than the ROI. Our saved templates are already
    cropped to 400x400, so pass wide_path = the ORIGINAL wide file (e.g. Muru's
    800x800 density_projection_*.fits) via --wide on the command line. The
    Coleman/F98 bulge maps are natively +-20 deg, i.e. the ROI itself, so no
    wider view exists for them and it is reported as None rather than silently
    reusing the cropped number.
    """
    l = (H.L_HALF - H.PIX / 2.0) - H.PIX * np.arange(H.NPIX)
    b = (-H.B_HALF + H.PIX / 2.0) + H.PIX * np.arange(H.NPIX)
    ll, bb = np.meshgrid(l, b)
    theta = np.hypot(ll, bb)

    src = np.asarray(fits.getdata(path), float)
    if src.shape == (H.NPIX, H.NPIX):
        cropped = src
    elif src.shape == (H.SRC_NPIX, H.SRC_NPIX):
        cropped = H.to_conventions_grid(src)
    else:
        return {"error": f"unrecognised template grid {src.shape}"}

    wide = None
    if wide_path:
        w = np.asarray(fits.getdata(wide_path), float)
        if w.shape != (H.SRC_NPIX, H.SRC_NPIX):
            return {"error": f"wide map {wide_path} is {w.shape}, "
                             f"expected {H.SRC_NPIX}^2"}
        wide = w

    out = {}
    for lv in (0.30, 0.20, 0.10, 0.05):
        m = cropped >= lv * cropped.max()
        if not m.any():
            continue
        rmax = float(theta[m].max())
        entry = {"mean_radius_deg": round(float(theta[m].mean()), 2),
                 "max_radius_deg": round(rmax, 2)}
        if rmax > max_radius_deg:
            entry["cropped"] = None
            entry["note"] = "contour reaches the ROI edge; cropped value would be an artifact"
        else:
            q, c4 = H._shape(ll[m], bb[m])
            entry["cropped"] = {"q": round(q, 4), "c4": round(c4, 4)}
        if wide is not None:
            qw, c4w, rmw, _ = H.shape_on_source(wide, lv)
            entry["uncropped"] = {"q": round(qw, 4), "c4": round(c4w, 4),
                                  "mean_radius_deg": round(rmw, 2)}
            if entry["cropped"]:
                entry["roi_artifact_c4"] = round(entry["cropped"]["c4"] - c4w, 4)
        else:
            entry["uncropped"] = None
            entry["note_uncropped"] = "source map is only +-20 deg; no wider view exists"
        out[f"{lv:.2f}"] = entry
    return out


# ------------------------------------------------- shape of a fitted map
def shape_of_model_map(model_counts, max_radius_deg=15.0):
    """Axis ratio and c4 of a fitted GCE model map, on the ROI grid.

    Afeefa asked for shape measured on the FITTED template rather than the raw
    one - i.e. after PSF convolution and exposure weighting, which is what the
    likelihood actually saw.

    The ROI is square, and a square crop fakes boxiness once a contour nears the
    boundary (that is Gap 1 pitfall 1 - it produced a false c4 = +0.02 in all
    six halos). Here the ROI *is* the analysis region so we cannot dodge it by
    using a wider map. Instead we refuse any contour level whose extent reaches
    past max_radius_deg and say so, rather than quietly returning a number.
    """
    l = (H.L_HALF - H.PIX / 2.0) - H.PIX * np.arange(H.NPIX)
    b = (-H.B_HALF + H.PIX / 2.0) + H.PIX * np.arange(H.NPIX)
    ll, bb = np.meshgrid(l, b)
    theta = np.hypot(ll, bb)

    out = {}
    for lv in (0.30, 0.20, 0.10, 0.05):
        m = model_counts >= lv * model_counts.max()
        if not m.any():
            continue
        rmax = float(theta[m].max())
        if rmax > max_radius_deg:
            out[f"{lv:.2f}"] = {
                "q": None, "c4": None, "max_radius_deg": round(rmax, 2),
                "skipped": "contour reaches the ROI edge; c4 would be a crop artifact"}
            continue
        q, c4 = H._shape(ll[m], bb[m])
        out[f"{lv:.2f}"] = {"q": round(q, 4), "c4": round(c4, 4),
                            "mean_radius_deg": round(float(theta[m].mean()), 2),
                            "max_radius_deg": round(rmax, 2)}
    return out


# ---------------------------------------------------------------- the fit
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--gce", action="append", required=True,
                    metavar="NAME=PATH", help="repeatable; e.g. hestia=/path.fits")
    ap.add_argument("--wide", action="append", default=[], metavar="NAME=PATH",
                    help="optional wide (+-40 deg) source map per template, used "
                         "for the uncropped shape column")
    ap.add_argument("--out", default="fit_results.json")
    args = ap.parse_args()

    templates = {}
    for spec in args.gce:
        name, _, path = spec.partition("=")
        if not path:
            raise SystemExit(f"--gce needs NAME=PATH, got {spec!r}")
        check_not_preconvolved(path)
        templates[name] = os.path.abspath(path)

    wides = {}
    for spec in args.wide:
        name, _, path = spec.partition("=")
        if name not in templates:
            raise SystemExit(f"--wide {name} has no matching --gce")
        wides[name] = os.path.abspath(path)

    from fermipy.gta import GTAnalysis

    results = {"config": args.config, "templates": templates, "stages": {}}

    # ---- setup ---------------------------------------------------------
    gta = GTAnalysis(args.config, logging={"verbosity": 3})
    gta.setup()

    # ---- STAGE 1: background + catalog only ----------------------------
    # Free the diffuse normalisations and the brightest inner sources only.
    # Freeing everything in this field does not converge in 5000 iterations.
    print("\n=== stage 1: background and catalog sources ===")
    gta.free_sources(free=False)
    gta.free_source("galdiff", free=True)
    gta.free_source("isodiff", free=True)
    gta.free_source("bubbles", free=True)
    gta.free_sources(minmax_ts=[100.0, None], distance=10.0, free=True,
                     pars="norm")
    fit1 = gta.fit()
    # fermipy's "loglike" is the log-likelihood itself, so a HIGHER value is a
    # better fit and differences are taken as (with) - (without).
    lnL_bkg = float(fit1["loglike"])
    results["stages"]["background_only"] = {
        "fit_quality": fit1.get("fit_quality"),
        "fit_status": fit1.get("fit_status"),
        "loglike": lnL_bkg,
    }
    print(f"  stage 1 loglike = {lnL_bkg:.3f}  quality={fit1.get('fit_quality')}")
    gta.write_roi("stage1_background")

    # ---- STAGE 2: add each GCE template in turn ------------------------
    for name, path in templates.items():
        print(f"\n=== stage 2: GCE = {name} ===")
        gta.load_roi("stage1_background")     # identical starting point every time

        srcname = f"gce_{name}"
        gta.add_source(srcname, {
            "SpatialModel": "SpatialMap",
            "Spatial_Filename": path,
            "SpectrumType": "PowerLaw",
            "Index": 2.4,                     # GCE is soft; refit below
            "Scale": 3000.0,
            "Prefactor": 1e-11,
        })
        gta.free_source(srcname, free=True)
        fit2 = gta.fit()

        ts = float(gta.roi[srcname]["ts"])
        norm = float(gta.roi[srcname]["param_values"][0])
        norm_err = float(gta.roi[srcname]["param_errors"][0])
        lnL = float(fit2["loglike"])

        # residual map with the GCE in the model
        resid = gta.residmap(f"resid_{name}", model={}, make_plots=True)

        # shape of what the likelihood actually fitted
        mm = gta.model_counts_map(srcname)
        model_2d = np.asarray(mm.data, float).sum(axis=0)   # sum over energy
        shape = shape_of_model_map(model_2d)

        results["stages"][name] = {
            "input_shape": shape_of_input(path, wides.get(name)),
            "ts": ts,
            "sqrt_ts_sigma": round(float(np.sqrt(max(ts, 0.0))), 2),
            "prefactor": norm,
            "prefactor_err": norm_err,
            "loglike": lnL,
            "delta_loglike_vs_background": lnL - lnL_bkg,
            "fit_quality": fit2.get("fit_quality"),
            "fit_status": fit2.get("fit_status"),
            "fitted_shape": shape,
            "residmap": f"resid_{name}",
        }
        print(f"  TS = {ts:.1f}  ({np.sqrt(max(ts,0)):.1f} sigma)")
        print(f"  prefactor = {norm:.3e} +/- {norm_err:.1e}")
        print(f"  d(loglike) vs background-only = {lnL - lnL_bkg:.2f}")
        inp = results["stages"][name]["input_shape"]
        print("  input template shape (cropped vs uncropped):")
        for lv, e in sorted(inp.items()):
            if isinstance(e, str) or "cropped" not in e:
                continue
            cr = e["cropped"]; un = e["uncropped"]
            cs = f"q={cr['q']:.3f} c4={cr['c4']:+.3f}" if cr else "skipped (ROI edge)"
            us = f"q={un['q']:.3f} c4={un['c4']:+.3f}" if un else "n/a (map is ROI-sized)"
            art = e.get("roi_artifact_c4")
            astr = f"   crop artifact in c4: {art:+.3f}" if art is not None else ""
            print(f"    {lv}: cropped {cs} | uncropped {us}{astr}")
        for lv, s in shape.items():
            if s.get("q") is None:
                print(f"  shape @ {lv}: skipped - {s['skipped']}")
            else:
                print(f"  shape @ {lv}: q={s['q']:.3f} c4={s['c4']:+.3f} "
                      f"(r~{s['mean_radius_deg']:.0f} deg)")

        gta.write_roi(f"stage2_{name}")
        gta.delete_source(srcname)

    # ---- head to head ---------------------------------------------------
    names = list(templates)
    if len(names) >= 2:
        a, b = names[0], names[1]
        dlnL = results["stages"][a]["loglike"] - results["stages"][b]["loglike"]
        # Both GCE models add the same number of free parameters here, so AIC
        # differences reduce to 2*dlnL. If that ever stops being true, compute
        # AIC properly instead of reusing this shortcut.
        results["comparison"] = {
            "a": a, "b": b,
            "delta_loglike_a_minus_b": dlnL,
            "delta_aic_a_minus_b": -2.0 * dlnL,
            "note": "equal parameter counts, so dAIC = -2*dlnL; recheck if that changes",
        }
        print(f"\n=== {a} vs {b} ===")
        print(f"  d(lnL) = {dlnL:+.2f}   dAIC = {-2*dlnL:+.2f} "
              f"({'favours '+a if dlnL>0 else 'favours '+b})")

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
