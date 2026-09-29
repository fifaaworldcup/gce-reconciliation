"""
H0_B: does iterative source-finding bias the recovered GCE morphology toward
sphericity?

The test. Build a simulated sky from our own fitted background plus a signal of
known shape, run the full Level-2 iteration on it, and measure the recovered
flattening BEFORE and AFTER the iteration on one fixed ruler. If the pipeline
rounds things off, the recovered q rises after iteration.

The ruler is a ladder of gNFW-rho^2 templates with vertical axis ratio
q = 0.4 ... 1.0, fitted one at a time; the best-fit q is the ladder member with
the highest lnL, refined by a parabola through its two neighbours. Both
injections are measured on the same ladder, so the quantity that matters is the
SHIFT between levels, not the absolute value - a HESTIA-shaped injection need
not land exactly on its own q when measured with a gNFW-profile ruler.

The spherical control is not optional. "Flattened comes back rounder" on its own
is equally consistent with a pipeline that mangles every morphology, so a
spherical injection at matched flux goes through the identical pipeline and must
come back unshifted.

The simulated sky contains the real 4FGL point sources, and the pipeline is
given them, exactly as Di Mauro's is. The iteration therefore only adds sources
BEYOND the catalogue, which is the thing under test.

Everything about the iteration - thresholds, stopping, separation - comes from
configs/level2.yaml and is identical for every injection.
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
import build_triaxial as TRI      # noqa: E402
import iterate as L2              # noqa: E402

LADDER_Q = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def ladder(cache=os.path.join(ROOT, "templates", "qladder")):
    """gNFW rho^2 templates flattened vertically by q. One fixed ruler."""
    os.makedirs(cache, exist_ok=True)
    out = {}
    for q in LADDER_Q:
        f = os.path.join(cache, f"qladder_{q:.1f}.npy")
        if os.path.exists(f):
            t = np.load(f)
        else:
            t = TRI.triaxial(q_plane=1.0, q_vert=q, theta_deg=0.0)
            np.save(f, t)
        out[q] = RG.normalise(t)
    return out


def best_q(counts, exposure, comps_base, energies, mask, rung):
    """Fit each ladder member; return (q_hat, lnL per rung, signal flux at the
    best rung). The flux is what the absorption bookkeeping compares."""
    lnl, flux = {}, {}
    for q, t in rung.items():
        comps = [{**b, "sig": t} for b in comps_base]
        r = FL.fit(counts, exposure, comps, energies, mask)
        lnl[q] = r["lnL_total"]
        flux[q] = sum(b["norms"]["sig"] for b in r["bins"]) * float(t.sum())
    qs = np.array(sorted(lnl)); ls = np.array([lnl[q] for q in qs])
    i = int(np.argmax(ls))
    if 0 < i < len(qs) - 1:                     # parabolic refinement
        y0, y1, y2 = ls[i - 1], ls[i], ls[i + 1]
        denom = (y0 - 2 * y1 + y2)
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        q_hat = float(qs[i] + np.clip(shift, -1, 1) * (qs[1] - qs[0]))
    else:
        q_hat = float(qs[i])
    q_near = min(rung, key=lambda q: abs(q - q_hat))
    return (q_hat, {float(k): float(v) for k, v in lnl.items()},
            float(flux[q_near]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--fit", default=os.path.join(ROOT, "outputs", "gap1",
                                                  "fit_maskB_ps",
                                                  "fit_results.json"),
                    help="real fit whose background normalisations define the "
                         "simulated sky")
    ap.add_argument("--flux-levels", default="0.5,1.0,2.0")
    ap.add_argument("--seed-catalog-frac", type=float, default=1.0,
                    help="fraction of catalogue sources, brightest first, that "
                         "the PIPELINE is given. The simulated sky always "
                         "contains all of them. 1.0 reproduces sec 6.6; a small "
                         "value forces the iteration to resolve the rest, which "
                         "is what drives its source count toward Di Mauro's.")
    ap.add_argument("--realizations", type=int, default=1,
                    help="noise realizations per configuration. One tells you "
                         "nothing about whether a shift is real.")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cfg = L2.load_config()
    rng_seed = int(cfg["seed"])

    counts0, e_lo, e_hi, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo >= a.emin - 1e-6) & (e_hi <= a.emax + 1e-6)
    e_lo, e_hi = e_lo[keep], e_hi[keep]
    emid = np.sqrt(e_lo * e_hi)
    exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
    iem, _ = FC.read_iem(a.iem, emid)
    mask, minfo = FC.build_masks(a.catalog, emid, "B")
    ps_cube, _ = PS.build(a.catalog, e_lo, e_hi)
    flat = RG.normalise(np.ones((400, 400)))
    bub, _ = FC.load_template(os.path.join(ROOT, "templates",
                                           "gap1_bubbles_fuzzy.fits"))
    T = os.path.join(ROOT, "templates")
    shapes = {
        "hestia": FC.load_template(os.path.join(T, "gap1_hestia_G1.1_ang0_1deg.fits"))[0],
        "spherical": FC.load_template(os.path.join(T, "gnfw_gamma1.2.fits"))[0],
    }

    # background truth: our own fitted normalisations, per bin
    allbins = json.load(open(a.fit))["models"]["hestia"]["norms_per_bin"]
    # match by ENERGY, never by position: the reference fit may be binned
    # differently from the cube we are simulating, and slicing silently pairs a
    # 1.2 GeV normalisation with a 4 GeV bin.
    ref_e = np.array([r["energy_gev"] for r in allbins])
    real = []
    for e in emid:
        k = int(np.argmin(np.abs(np.log(ref_e) - np.log(e))))
        if abs(np.log(ref_e[k] / e)) > np.log(1.6):
            raise SystemExit(f"no reference bin near {e:.2f} GeV in {a.fit} "
                             f"(closest {ref_e[k]:.2f}) - use a fit with "
                             f"matching binning")
        real.append(allbins[k])
    print("background truth taken from " + os.path.basename(os.path.dirname(a.fit))
          + ": bins " + ", ".join(f"{r['energy_gev']:.2f}" for r in real)
          + f" GeV for cube bins " + ", ".join(f"{e:.2f}" for e in emid),
          flush=True)
    bg_truth = [{"iso": r["iso"], "iem": r["iem"], "bubbles": r["bubbles"],
                 "ps": r["ps"]} for r in real]
    sig_truth = [r["hestia"] for r in real]      # the observed GCE level

    # The sky always contains every catalogue source. What the PIPELINE gets
    # may be only the brightest fraction of them, so the iteration has to
    # resolve the remainder - the regime Di Mauro's pipeline actually runs in.
    ps_truth = ps_cube
    if a.seed_catalog_frac >= 1.0:
        ps_seen = ps_cube
        n_withheld = 0
    else:
        from astropy.io import fits as _fits
        t = _fits.getdata(a.catalog, 1)
        tot = sum(PS.band_flux(t, lo * 1e3, hi * 1e3)
                  for lo, hi in zip(e_lo, e_hi))
        glon = np.asarray(t["GLON"], float)
        glon = np.where(glon > 180, glon - 360, glon)
        inroi = (np.abs(glon) < 20) & (np.abs(np.asarray(t["GLAT"], float)) < 20)
        cut = np.quantile(tot[inroi], 1.0 - a.seed_catalog_frac)
        keep_src = tot >= cut
        n_withheld = int((inroi & ~keep_src).sum())
        ps_seen = np.zeros_like(ps_cube)
        j = np.floor((20.0 - glon) / 0.1).astype(int)
        i = np.floor((np.asarray(t["GLAT"], float) + 20.0) / 0.1).astype(int)
        ins = (i >= 0) & (i < 400) & (j >= 0) & (j < 400) & keep_src
        for k, (lo, hi) in enumerate(zip(e_lo, e_hi)):
            fl = PS.band_flux(t, lo * 1e3, hi * 1e3)
            np.add.at(ps_seen[k], (i[ins], j[ins]), fl[ins])
        print(f"pipeline seeded with the brightest {a.seed_catalog_frac:.0%} of "
              f"the catalogue; {n_withheld} sources inside the ROI are in the "
              f"sky but NOT given to it", flush=True)

    bg_truth_maps = [{"iso": flat, "iem": RG.normalise(iem[e]), "bubbles": bub,
                      "ps": ps_truth[e]} for e in range(len(emid))]
    bg_maps = [{"iso": flat, "iem": RG.normalise(iem[e]), "bubbles": bub,
                "ps": ps_seen[e]} for e in range(len(emid))]
    rung = ladder()
    print(f"ruler: q = {LADDER_Q}", flush=True)
    print(f"mask B keeps {minfo['kept_fraction']:.1%}; "
          f"{len(emid)} bins {e_lo[0]:.2f}-{e_hi[-1]:.2f} GeV", flush=True)

    results = []
    for shape_name, shape in shapes.items():
        for lvl in [float(x) for x in a.flux_levels.split(",")]:
          for rz in range(a.realizations):
            t0 = time.time()
            rng = np.random.default_rng(rng_seed + 1000 * rz)
            counts = np.zeros((len(emid), 400, 400))
            for e in range(len(emid)):
                comp = {**bg_truth_maps[e], "sig": shape}
                conv = FL.convolve_templates(comp, emid[e])
                mu = exposure[e] * (
                    sum(bg_truth[e][k] * conv[k] for k in bg_truth_maps[e])
                    + lvl * sig_truth[e] * conv["sig"])
                counts[e] = rng.poisson(np.maximum(mu, 0))
            print(f"\n=== inject {shape_name} at {lvl}x the observed GCE flux, "
                  f"realization {rz}, {counts.sum():,.0f} photons", flush=True)

            # LEVEL 1: background as given, no source finding
            q1, lnl1, f1 = best_q(counts, exposure, bg_maps, emid, mask, rung)
            print(f"  level 1 recovered q = {q1:.3f}", flush=True)

            # LEVEL 2: same background, then iterate
            comps2, fit2, found = L2.run(counts, exposure, bg_maps, emid,
                                         mask=mask, config=cfg, verbose=True)
            base2 = [{k: v for k, v in c.items()} for c in comps2]
            q2, lnl2, f2 = best_q(counts, exposure, base2, emid, mask, rung)
            print(f"  level 2 recovered q = {q2:.3f}  "
                  f"after adding {len(found)} sources  "
                  f"({time.time() - t0:.0f}s)", flush=True)

            # source-absorption bookkeeping: how much flux did the added
            # sources take, against the signal that was injected?
            # absorption = the signal flux the iteration took away from the
            # GCE template, as a fraction of what level 1 attributed to it.
            # (Flux landing IN the new sources is not the right number: a new
            # source also eats background, so it is not bounded by the signal.)
            inj_flux = lvl * sum(sig_truth) * float(shape.sum())
            added_flux = sum(fit2["bins"][e]["norms"].get(f[0], 0.0)
                             for e in range(len(emid)) for f in found)
            absorbed = (f1 - f2) / f1 if f1 > 0 else None
            results.append({
                "injected_shape": shape_name, "flux_level": lvl,
                "realization": rz,
                "injected_signal_flux": inj_flux,
                "signal_flux_level1": f1, "signal_flux_level2": f2,
                "flux_in_added_sources": added_flux,
                "absorbed_fraction": absorbed,
                "q_level1": q1, "q_level2": q2, "delta_q": q2 - q1,
                "n_sources_added": len(found),
                "n_withheld": n_withheld,
                "seed_catalog_frac": a.seed_catalog_frac,
                "sources": [{"name": f[0], "ts": f[1], "l": f[2], "b": f[3],
                             "pass": f[4]} for f in found],
                "lnL_ladder_level1": lnl1, "lnL_ladder_level2": lnl2,
            })
            with open(os.path.join(a.out, "injection_recovery.json"), "w") as f:
                json.dump({"config": cfg, "ladder_q": LADDER_Q,
                           "results": results}, f, indent=2)

    print("\n================ H0_B ================")
    print(f"{'injected':10} {'flux':>5} {'n':>3} {'q lvl1':>8} {'q lvl2':>8} "
          f"{'shift':>16} {'src':>4} {'sig lost':>9}")
    import collections
    grp = collections.defaultdict(list)
    for r in results:
        grp[(r["injected_shape"], r["flux_level"])].append(r)
    for (sh, lvl), rs in sorted(grp.items()):
        d = np.array([r["delta_q"] for r in rs])
        q1 = np.mean([r["q_level1"] for r in rs])
        q2 = np.mean([r["q_level2"] for r in rs])
        ab = np.mean([(r["absorbed_fraction"] or 0.0) for r in rs])
        err = d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else float("nan")
        print(f"{sh:10} {lvl:5.1f} {len(rs):3d} {q1:8.3f} {q2:8.3f} "
              f"{d.mean():+8.3f} +- {err:5.3f} "
              f"{np.mean([r['n_sources_added'] for r in rs]):4.1f} {ab:8.1%}")
    def pooled(shape):
        d = np.array([r["delta_q"] for r in results
                      if r["injected_shape"] == shape])
        return d.mean(), (d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1
                          else float("nan"))
    fm, fe = pooled("hestia"); sm, se = pooled("spherical")
    print(f"\npooled shift: flattened {fm:+.3f} +- {fe:.3f}, "
          f"spherical control {sm:+.3f} +- {se:.3f}")
    print("a positive shift for the flattened injection AND ~0 for the control "
          "is what H0_B predicts against")
    return results


if __name__ == "__main__":
    main()
