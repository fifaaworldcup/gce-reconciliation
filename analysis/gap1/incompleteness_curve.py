"""Calibrate the catalogue-incompleteness bias on recovered flattening.

Sec 6.6 had two points: full catalogue, no bias; brightest 10%, a round
injection reading q = 0.57. This sweeps the catalogue fraction and plots the
bias against a quantity NPTF reports - unresolved point-source flux over GCE
flux - so Gap 5 supplies one number and the correction is read off.

The sky does not depend on how much of the catalogue the fit sees, so one sky
per realization is reused across every point and the Poisson noise cancels in
the paired difference. The bias is a level-1 effect, so no source finding runs
and each point costs seconds rather than minutes.

Same q-ladder as injection_recovery.py, imported not redeclared.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks"),
          os.path.join(ROOT, "analysis", "level2")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import iterate as L2              # noqa: E402
from injection_recovery import ladder, best_q, LADDER_Q   # noqa: E402

PIX = 0.1
NPIX = 400
HALF = 20.0
DEFAULT_FRACS = "1.0,0.9,0.75,0.6,0.5,0.4,0.3,0.2,0.15,0.1,0.05"


def catalogue_split(catalog_path, e_lo, e_hi, frac):
    """ps_seen cube for the brightest `frac` of the catalogue by count. The
    withheld sources are the faintest, as in a flux-limited catalogue."""
    from astropy.io import fits as _fits
    t = _fits.getdata(catalog_path, 1)
    glon = np.asarray(t["GLON"], float)
    glon = np.where(glon > 180, glon - 360, glon)
    glat = np.asarray(t["GLAT"], float)
    per_bin = np.array([PS.band_flux(t, lo * 1e3, hi * 1e3)
                        for lo, hi in zip(e_lo, e_hi)])       # (nE, nsrc)
    tot = per_bin.sum(axis=0)

    j = np.floor((HALF - glon) / PIX).astype(int)
    i = np.floor((glat + HALF) / PIX).astype(int)
    inside = (i >= 0) & (i < NPIX) & (j >= 0) & (j < NPIX)
    inner10 = inside & (np.abs(glon) < 10) & (np.abs(glat) < 10)

    if frac >= 1.0:
        keep = np.ones(len(tot), bool)
    else:
        cut = np.quantile(tot[inside], 1.0 - frac)
        keep = tot >= cut

    seen = np.zeros((len(e_lo), NPIX, NPIX))
    sel = inside & keep
    for k in range(len(e_lo)):
        np.add.at(seen[k], (i[sel], j[sel]), per_bin[k][sel])

    withheld = inside & ~keep
    book = {
        "catalog_frac_requested": float(frac),
        "n_in_roi": int(inside.sum()),
        "n_shown": int(sel.sum()),
        "n_withheld": int(withheld.sum()),
        "flux_roi_total": float(tot[inside].sum()),
        "flux_withheld": float(tot[withheld].sum()),
        "flux_withheld_inner10": float(tot[inner10 & ~keep].sum()),
        "flux_roi_inner10": float(tot[inner10].sum()),
        "faintest_shown": float(tot[sel].min()) if sel.any() else 0.0,
    }
    return seen, book


LINEAR_RANGE = 1.2      # see fit_local_slope

def fit_through_origin(x, y, yerr):
    """b = -k x, weighted, through the origin. sigma_k is inflated by
    sqrt(chi2/dof): adjacent points share most of their hidden sources and the
    ladder quantises, so the per-realization error understates the scatter."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    w = 1.0 / np.maximum(np.asarray(yerr, float), 1e-9) ** 2
    k = -np.sum(w * x * y) / np.sum(w * x * x)
    sk = 1.0 / np.sqrt(np.sum(w * x * x))
    resid = y + k * x
    dof = max(len(x) - 1, 1)
    chi2 = float(np.sum(w * resid ** 2) / dof)
    return float(k), float(sk * max(np.sqrt(chi2), 1.0)), chi2


def fit_local_slope(pts, xmax=LINEAR_RANGE):
    """Slope over f <= xmax, the regime NPTF reports.

    The curve saturates, so one line through the whole sweep is dragged flat
    by the far points and understates the slope near the origin by a factor of
    a few - in the direction that under-corrects sec 6.1. The cut was chosen
    after seeing the bend, so both fits are reported and the lookup table is
    the primary product."""
    use = [p for p in pts if 0.0 < p["f_unres_gce"] <= xmax]
    if len(use) < 3:
        return None
    k, sk, chi2 = fit_through_origin([p["f_unres_gce"] for p in use],
                                     [p["bias"] for p in use],
                                     [max(p["bias_err"], 0.002) for p in use])
    return {"slope_k": k, "slope_err": sk, "chi2_per_dof": chi2,
            "f_max": xmax, "n_points": len(use)}


def fit_power(x, y):
    """b = -A x^p over the whole sweep. p < 1 is the saturation. A
    description of the curve, not something to extrapolate with."""
    m = (np.asarray(y) < 0) & (np.asarray(x) > 0)
    if m.sum() < 3:
        return None
    lx = np.log(np.asarray(x)[m]); ly = np.log(-np.asarray(y)[m])
    p, lA = np.polyfit(lx, ly, 1)
    return {"A": float(np.exp(lA)), "p": float(p), "n_points": int(m.sum())}


def correct(pts, f_unres):
    """Correction interpolated from the table, no functional form assumed.
    Returns the bias (negative): q_true = q_measured - bias."""
    x = np.array([0.0] + [p["f_unres_gce"] for p in pts if p["f_unres_gce"] > 0])
    y = np.array([0.0] + [p["bias"] for p in pts if p["f_unres_gce"] > 0])
    o = np.argsort(x)
    return float(np.interp(f_unres, x[o], y[o]))


def summarise(rows, fracs, shape_names, out_dir, flux_level=1.0,
              realizations=None):
    """Paired bias and calibration from rows already computed."""
    # paired bias and calibration
    summary = {}
    print("\n================ catalogue-incompleteness calibration ==========")
    for shape_name in shape_names:
        R = [r for r in rows if r["shape"] == shape_name]
        ref = {r["realization"]: r["q"] for r in R if r["catalog_frac"] == 1.0}
        pts = []
        print(f"\n{shape_name}:  reference q = "
              f"{np.mean(list(ref.values())):.3f} "
              f"(spread {np.std(list(ref.values())):.3f} over "
              f"{len(ref)} realizations)")
        print(f"  {'f_cat':>6} {'hidden':>7} {'unres/GCE':>10} {'q':>7} "
              f"{'bias':>16}")
        for f in fracs:
            sub = [r for r in R if r["catalog_frac"] == f]
            d = np.array([r["q"] - ref[r["realization"]] for r in sub])
            err = (d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1
                   else float("nan"))
            x = float(np.mean([r["f_unres_gce"] for r in sub]))
            pts.append({"catalog_frac": f, "f_unres_gce": x,
                        "f_unres_ps": float(np.mean([r["f_unres_ps"]
                                                     for r in sub])),
                        "n_withheld": sub[0]["n_withheld"],
                        "q_mean": float(np.mean([r["q"] for r in sub])),
                        "bias": float(d.mean()),
                        "bias_err": float(err)})
            print(f"  {f:6.2f} {sub[0]['n_withheld']:7d} {x:10.3f} "
                  f"{np.mean([r['q'] for r in sub]):7.3f} "
                  f"{d.mean():+8.3f} +- {err:5.3f}")
        use = [p for p in pts if p["catalog_frac"] < 1.0]
        gk, gsk, gchi = fit_through_origin(
            [p["f_unres_gce"] for p in use], [p["bias"] for p in use],
            [max(p["bias_err"], 0.002) for p in use])
        loc = fit_local_slope(use)
        pw = fit_power([p["f_unres_gce"] for p in use],
                       [p["bias"] for p in use])
        summary[shape_name] = {
            "reference_q": float(np.mean(list(ref.values()))),
            "reference_q_spread": float(np.std(list(ref.values()))),
            "points": pts,
            "calibration": loc,
            "global_line": {"slope_k": gk, "slope_err": gsk,
                            "chi2_per_dof": gchi},
            "power_law": pw}
        if loc:
            print(f"  CALIBRATION (f <= {loc['f_max']}, {loc['n_points']} pts): "
                  f"bias = -{loc['slope_k']:.3f} +- {loc['slope_err']:.3f} x f"
                  f"   chi2/dof = {loc['chi2_per_dof']:.1f}")
        print(f"  whole sweep as one line: -{gk:.3f} +- {gsk:.3f} x f "
              f"(chi2/dof {gchi:.1f}) - flattened by saturation, do not use")
        if pw:
            print(f"  saturating form:    bias = -{pw['A']:.3f} f^{pw['p']:.2f}"
                  f"   (p < 1 is the saturation)")
    # attractor: the two injections start 0.25 apart. If hiding sources pulls
    # every morphology toward one value, the gap closes - and that value is
    # what an incomplete catalogue makes any excess look like.
    attractor = []
    print("\n---- do the two injections converge on one q? ----")
    print(f"  {'f':>7} {'q flattened':>12} {'q spherical':>12} {'gap':>7}")
    for f in fracs:
        qs = {}
        for shape_name in shape_names:
            sub = [r for r in rows if r["shape"] == shape_name
                   and r["catalog_frac"] == f]
            if sub:
                qs[shape_name] = float(np.mean([r["q"] for r in sub]))
                x = float(np.mean([r["f_unres_gce"] for r in sub]))
        if len(qs) == 2:
            gap = abs(qs["spherical"] - qs["hestia"])
            attractor.append({"f_unres_gce": x, "q_hestia": qs["hestia"],
                              "q_spherical": qs["spherical"], "gap": gap,
                              "midpoint": 0.5 * (qs["hestia"]
                                                 + qs["spherical"])})
            print(f"  {x:7.3f} {qs['hestia']:12.3f} "
                  f"{qs['spherical']:12.3f} {gap:7.3f}")
    if attractor:
        far = [p for p in attractor if p["f_unres_gce"] > 2.0]
        if far:
            print(f"  gap closes from {attractor[0]['gap']:.3f} at f=0 to "
                  f"{np.mean([p['gap'] for p in far]):.3f} beyond f=2, "
                  f"converging near q = "
                  f"{np.mean([p['midpoint'] for p in far]):.2f}")

    with open(os.path.join(out_dir, "calibration.json"), "w") as fh:
        json.dump({"flux_level": flux_level,
                   "realizations": (realizations if realizations is not None
                                    else len({r["realization"]
                                              for r in rows})),
                   "summary": summary, "attractor": attractor}, fh, indent=2)

    print("\nTo use: q_true = q_measured - bias(f), with f = unresolved")
    print("point-source flux / GCE flux from Gap 5's NPTF, in the same ROI and")
    print("energy range. Read off the measured curve, no functional form:")
    print(f"  {'f':>6} " + " ".join(f"{k:>22}" for k in summary))
    for fu in (0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00):
        cells = []
        for k, sv in summary.items():
            b = correct(sv["points"], fu)
            cells.append(f"q_true = q_meas {-b:+.3f}")
        print(f"  {fu:6.2f} " + " ".join(f"{c:>22}" for c in cells))
    return summary



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--fit", default=os.path.join(ROOT, "outputs", "gap1",
                                                  "fit_maskB_ps",
                                                  "fit_results.json"))
    ap.add_argument("--fracs", default=DEFAULT_FRACS)
    ap.add_argument("--flux-level", type=float, default=1.0)
    ap.add_argument("--realizations", type=int, default=5)
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cfg = L2.load_config()
    rng_seed = int(cfg["seed"])
    fracs = [float(x) for x in a.fracs.split(",")]
    if 1.0 not in fracs:
        raise SystemExit("fracs must include 1.0 - it is the reference every "
                         "other point is differenced against")

    counts0, e_lo, e_hi, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo >= a.emin - 1e-6) & (e_hi <= a.emax + 1e-6)
    e_lo, e_hi = e_lo[keep], e_hi[keep]
    emid = np.sqrt(e_lo * e_hi)
    exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
    iem, _ = FC.read_iem(a.iem, emid)
    mask, minfo = FC.build_masks(a.catalog, emid, "B")
    ps_cube, _ = PS.build(a.catalog, e_lo, e_hi)
    flat = RG.normalise(np.ones((NPIX, NPIX)))
    bub, _ = FC.load_template(os.path.join(ROOT, "templates",
                                           "gap1_bubbles_fuzzy.fits"))
    T = os.path.join(ROOT, "templates")
    shapes = {
        "hestia": FC.load_template(os.path.join(T,
                    "gap1_hestia_G1.1_ang0_1deg.fits"))[0],
        "spherical": FC.load_template(os.path.join(T, "gnfw_gamma1.2.fits"))[0],
    }

    allbins = json.load(open(a.fit))["models"]["hestia"]["norms_per_bin"]
    ref_e = np.array([r["energy_gev"] for r in allbins])
    real = []
    for e in emid:
        k = int(np.argmin(np.abs(np.log(ref_e) - np.log(e))))
        if abs(np.log(ref_e[k] / e)) > np.log(1.6):
            raise SystemExit(f"no reference bin near {e:.2f} GeV in {a.fit}")
        real.append(allbins[k])
    bg_truth = [{"iso": r["iso"], "iem": r["iem"], "bubbles": r["bubbles"],
                 "ps": r["ps"]} for r in real]
    sig_truth = [r["hestia"] for r in real]

    # catalogue ladder, built once
    seen_cubes, books = {}, {}
    for f in fracs:
        seen_cubes[f], books[f] = catalogue_split(a.catalog, e_lo, e_hi, f)
    # f=1 must reproduce the sky's own PS cube, or the withheld effect is
    # contaminated by a reconstruction mismatch
    dev = float(np.max(np.abs(seen_cubes[1.0] - ps_cube)))
    scale = float(np.max(np.abs(ps_cube)))
    print(f"catalogue reconstruction check: max |seen(f=1) - truth| = {dev:.3e} "
          f"({dev / max(scale, 1e-300):.2e} of peak)", flush=True)
    if dev > 1e-12 * max(scale, 1e-300):
        raise SystemExit("f=1 reconstruction does not match the sky's own "
                         "point-source cube - fix that before reading any "
                         "bias off this curve")

    bg_truth_maps = [{"iso": flat, "iem": RG.normalise(iem[e]), "bubbles": bub,
                      "ps": ps_cube[e]} for e in range(len(emid))]
    rung = ladder()
    print(f"ruler: q = {LADDER_Q};  mask B keeps {minfo['kept_fraction']:.1%};  "
          f"{len(emid)} bins {e_lo[0]:.2f}-{e_hi[-1]:.2f} GeV", flush=True)
    print(f"catalogue ladder: " + ", ".join(
        f"{f:.2f}->{books[f]['n_withheld']}src" for f in fracs), flush=True)

    rows = []
    for shape_name, shape in shapes.items():
        inj_flux = a.flux_level * sum(sig_truth) * float(shape.sum())
        for rz in range(a.realizations):
            # one sky, reused for every catalogue fraction
            rng = np.random.default_rng(rng_seed + 1000 * rz)
            counts = np.zeros((len(emid), NPIX, NPIX))
            for e in range(len(emid)):
                comp = {**bg_truth_maps[e], "sig": shape}
                conv = FL.convolve_templates(comp, emid[e])
                mu = exposure[e] * (
                    sum(bg_truth[e][k] * conv[k] for k in bg_truth_maps[e])
                    + a.flux_level * sig_truth[e] * conv["sig"])
                counts[e] = rng.poisson(np.maximum(mu, 0))
            print(f"\n=== {shape_name}, realization {rz}: "
                  f"{counts.sum():,.0f} photons", flush=True)

            for f in fracs:
                t0 = time.time()
                bg_maps = [{"iso": flat, "iem": RG.normalise(iem[e]),
                            "bubbles": bub, "ps": seen_cubes[f][e]}
                           for e in range(len(emid))]
                q, lnl, sf = best_q(counts, exposure, bg_maps, emid, mask, rung)
                b = books[f]
                rows.append({
                    "shape": shape_name, "realization": rz,
                    "catalog_frac": f,
                    "q": q, "signal_flux": sf, "lnL_ladder": lnl,
                    "injected_signal_flux": inj_flux,
                    "f_unres_gce": b["flux_withheld"] / inj_flux,
                    "f_unres_ps": b["flux_withheld"] / b["flux_roi_total"],
                    "f_unres_gce_inner10": (b["flux_withheld_inner10"]
                                            / inj_flux),
                    **{k: v for k, v in b.items() if k != "catalog_frac_requested"},
                })
                print(f"  f={f:5.2f}  {b['n_withheld']:4d} hidden  "
                      f"unres/GCE={rows[-1]['f_unres_gce']:6.3f}  "
                      f"q={q:.3f}  ({time.time() - t0:.0f}s)", flush=True)
            with open(os.path.join(a.out, "incompleteness_curve.json"), "w") as fh:
                json.dump({"config": cfg, "ladder_q": LADDER_Q,
                           "flux_level": a.flux_level, "rows": rows}, fh,
                          indent=2)

    summary = summarise(rows, fracs, list(shapes), a.out,
                        flux_level=a.flux_level,
                        realizations=a.realizations)
    return summary


if __name__ == "__main__":
    main()
