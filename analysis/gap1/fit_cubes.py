"""
Gap 1 Level-1 fit: run the template comparison on the real counts cube.

This is the step Afeefa asked for - gNFW vs HESTIA vs VVV, Mask A, 1-10 GeV,
recording TS, best-fit normalisation, delta ln L, residual maps, and both the
cropped and uncropped q and c4.

Division of labour, because Fermitools is conda-only and the sandbox can't
reach fermi.gsfc.nasa.gov:

  laptop   gtselect/gtmktime/gtbin/gtltcube/gtexpcube2 -> counts + exposure cube
           reduce_iem.py                              -> diffuse model, ROI cut
  here     everything else

Why not gtlike. gtlike models every point source from an XML catalogue; the
project design masks them instead (Conventions section 3, Mask A). Once they're
masked there is nothing left for gtlike to do that the EM fitter in
fit_likelihood.py doesn't already do, and that one is closed-loop tested. It
also means no make4FGLxml step, which is a large and fragile dependency.

Each energy bin gets its own free normalisation per component, so no spectral
model is assumed anywhere and only the per-bin SPATIAL shape matters. That is
what makes this a morphology test rather than a spectral one.

Usage:
    python fit_cubes.py --ccube gap1_ccube.fits --expcube gap1_expcube.fits \
        --iem gap1_iem_roi.fits --catalog gll_psc_catalog.fit --out outputs/gap1/fit
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_likelihood as FL          # noqa: E402
import regrid as RG                  # noqa: E402
import run_fit as RF                 # noqa: E402
import make_masks as MM              # noqa: E402
import ps_template as PS             # noqa: E402

NPIX = 400
PIX = 0.1


# --------------------------------------------------------------- readers
def read_ccube(path):
    """Counts cube from gtbin CCUBE. Returns (counts, e_lo, e_hi) in GeV."""
    with fits.open(path) as hd:
        counts = np.asarray(hd[0].data, float)
        hdr = hd[0].header
        eb = None
        for e in hd[1:]:
            if e.name.upper() == "EBOUNDS":
                eb = e.data
        if eb is None:
            raise SystemExit(f"{path} has no EBOUNDS extension - is it a CCUBE?")
        # EBOUNDS is keV by the FITS standard, whatever the gtbin inputs were
        e_lo = np.asarray(eb["E_MIN"], float) / 1.0e6
        e_hi = np.asarray(eb["E_MAX"], float) / 1.0e6

    if counts.ndim != 3:
        raise SystemExit(f"counts cube is {counts.ndim}D, expected 3D")
    if counts.shape[1:] != (NPIX, NPIX):
        raise SystemExit(f"counts cube is {counts.shape[1:]}, expected "
                         f"{(NPIX, NPIX)} - was gtbin run with the project grid?")
    if counts.shape[0] != len(e_lo):
        raise SystemExit(f"{counts.shape[0]} image planes but {len(e_lo)} "
                         f"EBOUNDS rows")
    return counts, e_lo, e_hi, hdr


def _log_interp_cube(values, grid_e, target_e):
    """Interpolate a (nP, ny, nx) stack in log-log to one energy."""
    lg = np.log(np.asarray(grid_e, float))
    lt = float(np.log(target_e))
    if len(lg) == 1:
        return np.asarray(values[0], float)
    j = int(np.clip(np.searchsorted(lg, lt), 1, len(lg) - 1))
    f = (lt - lg[j - 1]) / (lg[j] - lg[j - 1])
    v0 = np.log(np.maximum(values[j - 1], 1e-300))
    v1 = np.log(np.maximum(values[j], 1e-300))
    return np.exp(v0 + f * (v1 - v0))


def read_expcube(path, energies_gev, ccube_hdr):
    """Exposure on the counts-cube grid at the counts-cube bin energies.

    gtexpcube2 is deliberately run wider than the ROI (440 vs 400 pixels) so
    gtsrcmaps has margin for the PSF wings, and it writes exposure at the bin
    BOUNDARIES, not the centres. Both have to be undone here. The crop offset
    is taken from the two WCSs rather than assumed to be (440-400)/2, because
    if either geometry ever changes a hardcoded 20 would shift the exposure
    silently and bias every normalisation.
    """
    with fits.open(path) as hd:
        exp = np.asarray(hd[0].data, float)
        ewcs = WCS(hd[0].header).celestial
        en = None
        for e in hd[1:]:
            if e.name.upper() in ("ENERGIES", "EBOUNDS"):
                col = e.data.names[0]
                en = np.asarray(e.data[col], float) / 1.0e3   # MeV -> GeV
                break
    if en is None:
        raise SystemExit(f"{path} has no ENERGIES extension")
    if exp.ndim != 3:
        raise SystemExit(f"exposure cube is {exp.ndim}D, expected 3D")

    cwcs = WCS(ccube_hdr).celestial
    # where does the counts grid's first pixel centre land on the exposure grid
    l0, b0 = cwcs.wcs_pix2world(0.0, 0.0, 0)
    x0, y0 = ewcs.wcs_world2pix(l0, b0, 0)
    ox, oy = int(round(float(x0))), int(round(float(y0)))
    if (abs(float(x0) - ox) > 1e-3) or (abs(float(y0) - oy) > 1e-3):
        raise SystemExit(f"exposure grid is offset from the counts grid by a "
                         f"non-integer number of pixels ({x0:.3f}, {y0:.3f}); "
                         f"resampling would be needed and isn't implemented")
    if ox < 0 or oy < 0 or ox + NPIX > exp.shape[2] or oy + NPIX > exp.shape[1]:
        raise SystemExit(f"exposure cube {exp.shape[1:]} does not cover the "
                         f"counts grid at offset ({ox}, {oy})")
    exp = exp[:, oy:oy + NPIX, ox:ox + NPIX]

    out = np.stack([_log_interp_cube(exp, en, E) for E in energies_gev])
    info = {"crop_offset_px": [ox, oy],
            "planes_in_file": int(len(en)),
            "energies_gev": [round(float(x), 4) for x in en],
            "bin_energies_gev": [round(float(x), 4) for x in energies_gev]}
    return out, info


def read_iem(path, energies_gev):
    """Galactic diffuse model, regridded and interpolated to our bins."""
    with fits.open(path) as hd:
        cube = np.asarray(hd[0].data, float)
        w = WCS(hd[0].header).celestial
        en = None
        for e in hd[1:]:
            if e.name.upper() in ("ENERGIES", "EBOUNDS"):
                col = e.data.names[0]
                en = np.asarray(e.data[col], float) / 1.0e3   # MeV -> GeV
                break
    if en is None:
        raise SystemExit(f"{path} has no ENERGIES extension")
    if cube.ndim != 3 or cube.shape[0] != len(en):
        raise SystemExit(f"diffuse cube {cube.shape} vs {len(en)} energies")

    planes, cover = [], []
    for k in range(cube.shape[0]):
        m, c = RG.regrid_array(cube[k], w)
        planes.append(m)
        cover.append(c)
    planes = np.stack(planes)

    out = np.stack([_log_interp_cube(planes, en, E) for E in energies_gev])
    return out, {"coverage_min": round(float(np.min(cover)), 4),
                 "planes": int(cube.shape[0]),
                 "energy_range_gev": [round(float(en[0]), 4),
                                      round(float(en[-1]), 4)]}


# ------------------------------------------------------------ components
def load_template(path):
    m, cover = RG.load(path)
    if not np.isfinite(m).all():
        raise SystemExit(f"{path} has non-finite pixels")
    return RG.normalise(m), cover


def build_masks(catalog, energies_gev, which="A", e_lo=None, e_hi=None, iem=None):
    if catalog and os.path.exists(catalog):
        kw = {}
        if which == "E":
            if iem is None:
                raise SystemExit("mask E needs the diffuse model: pass --iem")
            kw = {"e_lo": e_lo, "e_hi": e_hi, "iem": iem}
        cube = MM.build_at_energies(catalog, which, energies_gev, **kw)
        info = {"mask": which, "catalog": os.path.basename(catalog),
                "kept_fraction": round(float(cube.mean()), 4)}
        if which == "E":
            info["mask_e_fraction"] = MM.MASK_E_FRACTION
        return cube, info
    return (np.ones((len(energies_gev), NPIX, NPIX), np.float32),
            {"mask": "NONE", "catalog": None, "kept_fraction": 1.0,
             "warning": "no catalogue supplied - point sources are NOT masked, "
                        "so every normalisation is biased high and the TS "
                        "values are meaningless as absolute numbers"})


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem")
    ap.add_argument("--catalog")
    ap.add_argument("--mask", default="A", choices=list("ABCDE"))
    ap.add_argument("--templates", default=os.path.join(ROOT, "templates"))
    ap.add_argument("--signals", default="gnfw,hestia,vvv",
                    help="comma separated: gnfw, hestia, vvv, f98")
    ap.add_argument("--hestia", default="gap1_hestia_G1.1_ang0_1deg.fits")
    ap.add_argument("--hestia-wide",
                    help="Muru's original 800x800 map for the same galaxy and "
                         "angle. Needed for the UNCROPPED q and c4 - our saved "
                         "templates are already cut to the ROI.")
    ap.add_argument("--bubbles", action="store_true",
                    help="add the fuzzy bubbles template as a 4th component")
    ap.add_argument("--ps", action="store_true",
                    help="add the 4FGL point-source template (physical units, "
                         "fitted norm should come out near 1)")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "gap1", "fit"))
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rec = {"inputs": {}, "components": {}, "models": {}, "comparisons": {},
           "caveats": []}

    # ---- data
    counts, e_lo, e_hi, chdr = read_ccube(args.ccube)
    emid = np.sqrt(e_lo * e_hi)
    keep = (e_lo >= args.emin - 1e-9) & (e_hi <= args.emax + 1e-9)
    if not keep.any():
        raise SystemExit(f"no counts-cube bin lies inside "
                         f"{args.emin}-{args.emax} GeV; bins are "
                         f"{list(zip(e_lo.round(3), e_hi.round(3)))}")
    if not keep.all():
        print(f"using {keep.sum()}/{len(keep)} bins inside "
              f"{args.emin}-{args.emax} GeV")
    counts, e_lo, e_hi, emid = counts[keep], e_lo[keep], e_hi[keep], emid[keep]
    print(f"counts cube: {counts.shape}, {counts.sum():,.0f} photons, "
          f"{e_lo[0]:.2f}-{e_hi[-1]:.2f} GeV")

    # An empty bin is never physical here - the GC is bright across this whole
    # range - so it means the events file doesn't cover the band being asked
    # for. That happened once already: the FSSC download was made with a 5 GeV
    # floor instead of 500 MeV, and gtselect can only intersect with what it is
    # given, so five of eight bins came through as exact zeros. Fitting that
    # returns normalisations of zero and a TS of zero, which reads like a null
    # result rather than a missing file. Refuse instead.
    per_bin = counts.reshape(counts.shape[0], -1).sum(1)
    empty = np.flatnonzero(per_bin == 0)
    if empty.size:
        rng = ", ".join(f"{e_lo[i]:.2f}-{e_hi[i]:.2f} GeV" for i in empty)
        raise SystemExit(
            f"{empty.size} of {len(per_bin)} energy bins contain zero counts "
            f"({rng}).\nThe counts cube does not cover the requested band. "
            f"Check the ENERGY selection recorded in the events file "
            f"(DSVAL for DSTYP=ENERGY) against --emin/--emax.")

    # a bin far below its neighbours is the partial version of the same fault
    nz = per_bin[per_bin > 0]
    if nz.size >= 3 and nz.min() < 0.02 * np.median(nz):
        print(f"  WARNING one bin holds {nz.min():,.0f} counts against a "
              f"median of {np.median(nz):,.0f} - check the events file covers "
              f"the full band, this is what a clipped energy range looks like")
        rec["caveats"].append("an energy bin is far below its neighbours; the "
                              "events file may not span the full fit range")

    exposure, expinfo = read_expcube(args.expcube, emid, chdr)
    print(f"exposure: crop offset {expinfo['crop_offset_px']}, "
          f"median {np.median(exposure):.3e} cm2 s")
    rec["inputs"] = {"ccube": os.path.basename(args.ccube),
                     "expcube": os.path.basename(args.expcube),
                     "photons": int(counts.sum()),
                     "bins_gev": [[round(float(a), 3), round(float(b), 3)]
                                  for a, b in zip(e_lo, e_hi)],
                     "exposure": expinfo}

    # ---- the diffuse model comes first: mask D is defined against it
    iem_cube = None
    if args.iem and os.path.exists(args.iem):
        iem_cube, iemi = read_iem(args.iem, emid)
        print(f"diffuse model: {iemi['planes']} planes, coverage "
              f"{iemi['coverage_min']:.1%}")
        rec["components"]["iem"] = iemi

    # ---- masks
    mask, maskinfo = build_masks(args.catalog, emid, args.mask,
                                 e_lo=e_lo, e_hi=e_hi, iem=iem_cube)
    print(f"mask {maskinfo['mask']}: keeps {maskinfo['kept_fraction']:.1%}")
    if "warning" in maskinfo:
        print("  WARNING " + maskinfo["warning"])
        rec["caveats"].append(maskinfo["warning"])
    rec["inputs"]["mask"] = maskinfo

    # ---- background components
    background = {}
    if iem_cube is not None:
        iem_per_bin = [RG.normalise(iem_cube[e]) for e in range(len(emid))]
    else:
        msg = ("no galactic diffuse model supplied - the dominant background "
               "is missing, so the normalisations, TS values and residual "
               "maps here describe the pipeline, not the sky")
        print("WARNING " + msg)
        rec["caveats"].append(msg)
        iem_per_bin = None

    flat = RG.normalise(np.ones((NPIX, NPIX)))
    ps_cube = None
    if args.ps:
        if not (args.catalog and os.path.exists(args.catalog)):
            raise SystemExit("--ps needs --catalog")
        ps_cube, psi = PS.build(args.catalog, e_lo, e_hi)
        rec["components"]["ps"] = psi
        print(f"point-source template: {psi['sources_placed']} sources")
    bub = None
    if args.bubbles:
        bub, _ = load_template(os.path.join(args.templates,
                                            "gap1_bubbles_fuzzy.fits"))

    # ---- signal templates
    paths = {"gnfw": "gnfw_gamma1.2.fits",
             "hestia": args.hestia,
             "vvv": "gap1_bulge_coleman.fits",
             "f98": "gap1_bulge_f98.fits"}
    signals = {}
    for name in [s.strip() for s in args.signals.split(",") if s.strip()]:
        fp = os.path.join(args.templates, paths[name])
        if not os.path.exists(fp):
            raise SystemExit(f"missing template {fp}")
        t, cover = load_template(fp)
        signals[name] = t
        rec["components"][name] = {"file": paths[name],
                                   "coverage": round(cover, 4)}
        print(f"signal template {name:7} <- {paths[name]}  "
              f"coverage {cover:.1%}")

    # ---- fit each signal model against the same background
    for name, tmpl in signals.items():
        print(f"\n=== fitting {name} ===")
        comps = []
        for e in range(len(emid)):
            c = {"iso": flat, name: tmpl}
            if iem_per_bin is not None:
                c["iem"] = iem_per_bin[e]
            if bub is not None:
                c["bubbles"] = bub
            if ps_cube is not None:
                c["ps"] = ps_cube[e]
            comps.append(c)

        ts, with_it, without = FL.ts_of(counts, exposure, comps, emid, name,
                                       mask_cube=mask)
        resid = FL.residual_map(counts, exposure, comps, emid, with_it,
                                mask_cube=mask)

        # the fitted signal map: what the likelihood actually saw, i.e. after
        # PSF convolution and exposure weighting, summed over the fitted bins
        model = np.zeros((NPIX, NPIX))
        for e in range(len(emid)):
            conv = FL.convolve_templates({name: comps[e][name]}, emid[e])[name]
            model += with_it["bins"][e]["norms"][name] * exposure[e] * conv

        shape_fitted = RF.shape_of_model_map(model)
        # uncropped needs a view wider than the ROI. only HESTIA has one:
        # gNFW is circular so a square crop can't distort it until r > 20 deg,
        # and the Coleman map is natively exactly the ROI, so no wider view of
        # it exists. shape_of_input reports those as None rather than reusing
        # the cropped number.
        wide = args.hestia_wide if name == "hestia" else None
        shape_input = RF.shape_of_input(os.path.join(args.templates,
                                                     paths[name]),
                                        wide_path=wide)

        counts_in_model = float(model.sum())
        print(f"  TS {ts:,.1f}   ({np.sqrt(max(ts, 0)):.1f} sigma)")
        print(f"  lnL {with_it['lnL_total']:,.1f}  "
              f"(without {name}: {without['lnL_total']:,.1f})")
        print(f"  model photons attributed to {name}: {counts_in_model:,.0f}")
        print(f"  residual: mean {np.nanmean(resid):+.3f} sigma, "
              f"std {np.nanstd(resid):.3f}")

        fits.PrimaryHDU(resid.astype(np.float32)).writeto(
            os.path.join(args.out, f"residual_{name}.fits"), overwrite=True)
        fits.PrimaryHDU(model.astype(np.float32)).writeto(
            os.path.join(args.out, f"model_{name}.fits"), overwrite=True)

        rec["models"][name] = {
            "TS": round(float(ts), 2),
            "sigma": round(float(np.sqrt(max(ts, 0.0))), 2),
            "lnL_total": round(float(with_it["lnL_total"]), 3),
            "lnL_without": round(float(without["lnL_total"]), 3),
            "delta_lnL_vs_no_signal": round(float(
                with_it["lnL_total"] - without["lnL_total"]), 3),
            "model_photons": round(counts_in_model, 1),
            "norms_per_bin": [{"energy_gev": round(float(emid[e]), 3),
                               **{k: float(v) for k, v in
                                  with_it["bins"][e]["norms"].items()}}
                              for e in range(len(emid))],
            "all_bins_converged": all(b["converged"]
                                      for b in with_it["bins"]),
            "residual_mean_sigma": round(float(np.nanmean(resid)), 4),
            "residual_std_sigma": round(float(np.nanstd(resid)), 4),
            "shape_of_fitted_map": shape_fitted,
            "shape_of_input_template": shape_input,
        }

    # ---- model comparison. same number of free parameters either way, so
    # AIC reduces to -2 dlnL and the pre-registered threshold applies directly
    names = list(signals)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = (rec["models"][b]["lnL_total"] - rec["models"][a]["lnL_total"])
            rec["comparisons"][f"{b}_minus_{a}"] = {
                "delta_lnL": round(float(d), 3),
                "delta_AIC": round(float(-2.0 * d), 3),
                "favoured": b if d > 0 else a,
                "decisive_at_dAIC10": bool(abs(2.0 * d) >= 10.0),
            }
            print(f"\n{b} vs {a}: dlnL = {d:+.2f}, dAIC = {-2 * d:+.2f}  "
                  f"-> {'favours ' + (b if d > 0 else a)}"
                  f"{'  (decisive)' if abs(2 * d) >= 10 else '  (not decisive)'}")

    with open(os.path.join(args.out, "fit_results.json"), "w") as f:
        json.dump(rec, f, indent=2)
    print(f"\nwrote {os.path.join(args.out, 'fit_results.json')}")
    if rec["caveats"]:
        print("\nCAVEATS")
        for c in rec["caveats"]:
            print("  - " + c)
    return rec


if __name__ == "__main__":
    main()
