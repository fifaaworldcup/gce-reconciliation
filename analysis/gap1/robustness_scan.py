"""
Robustness checks on the first Gap 1 fit. Three questions a reader should ask
before believing "HESTIA beats VVV beats gNFW":

  1. Is it one lucky halo?  Fit all six HESTIA galaxies x four viewing angles
     and report the envelope, not the best member (mini-paper pitfall 5).
  2. Is gNFW a straw man?  The first fit froze gamma = 1.2. Conventions sec 5
     say gamma is free in [0.8, 1.4]. Scan it.
  3. Is it smoothing, not shape?  Muru's maps carry a 1 deg top-hat on top of
     the PSF; the gNFW does not. So every gNFW is fitted raw AND through the
     identical top-hat, which isolates core smoothing from morphology.

Plus the template correlation matrix the mini-paper promised beside every
delta-AIC: two templates correlated at ~0.95 inside the mask can't be told
apart, and a delta-AIC between them would measure the background, not them.

Every signal is fitted against the same background (iso + diffuse + bubbles),
one free normalisation per component per bin, so every model has the same
number of parameters and delta-AIC = -2 delta-lnL throughout.
"""

import argparse
import csv
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks")):
    sys.path.insert(0, p)

import fit_cubes as FC               # noqa: E402
import fit_likelihood as FL          # noqa: E402
import regrid as RG                  # noqa: E402
import project_template as PT        # noqa: E402
import hestia_to_template as H       # noqa: E402
import ps_template as PS             # noqa: E402

GAMMAS = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]


def gnfw_templates(cache):
    os.makedirs(cache, exist_ok=True)
    out = {}
    for g in GAMMAS:
        path = os.path.join(cache, f"gnfw_g{g:.1f}.npy")
        if os.path.exists(path):
            raw = np.load(path)
        else:
            raw = PT.normalise(PT.analytic(PT.gnfw(g)))
            np.save(path, raw)
        out[f"gnfw_g{g:.1f}"] = RG.normalise(raw)
        out[f"gnfw_g{g:.1f}_th1"] = RG.normalise(H.tophat_smooth(raw, 1.0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--mask", default="B", choices=list("ABCDE"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--library", default="full", choices=["full", "triaxial"],
                    help="triaxial: the four Hu et al. readings against VVV, "
                         "gNFW and the fiducial HESTIA, for the Gap 4 cross-check")
    ap.add_argument("--emin", type=float, default=0.0)
    ap.add_argument("--emax", type=float, default=1e9)
    ap.add_argument("--ps", action="store_true",
                    help="add the 4FGL point-source template to the background")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = os.path.join(ROOT, "templates")

    counts, e_lo, e_hi, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo >= a.emin - 1e-6) & (e_hi <= a.emax + 1e-6)
    if not keep.any():
        raise SystemExit("no bins in the requested energy range")
    counts, e_lo, e_hi = counts[keep], e_lo[keep], e_hi[keep]
    emid = np.sqrt(e_lo * e_hi)
    print(f"energy {e_lo[0]:.2f}-{e_hi[-1]:.2f} GeV, {keep.sum()} bins", flush=True)
    exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
    iem, _ = FC.read_iem(a.iem, emid)
    mask, minfo = FC.build_masks(a.catalog, emid, a.mask,
                                 e_lo=e_lo, e_hi=e_hi, iem=iem)
    flat = RG.normalise(np.ones((400, 400)))
    bub, _ = FC.load_template(os.path.join(T, "gap1_bubbles_fuzzy.fits"))
    bg = [{"iso": flat, "iem": RG.normalise(iem[e]), "bubbles": bub}
          for e in range(len(emid))]
    if a.ps:
        ps, psinfo = PS.build(a.catalog, e_lo, e_hi)
        for e in range(len(emid)):
            bg[e]["ps"] = ps[e]          # physical units: expect norm ~ 1
        print(f"point-source template: {psinfo['sources_placed']} sources",
              flush=True)
    print(f"mask {a.mask}: keeps {minfo['kept_fraction']:.1%}, "
          f"{counts.sum():,.0f} photons, {len(emid)} bins", flush=True)

    # ---- the signal library
    sig = {}
    pat = re.compile(r"gap1_hestia_(G\d\.\d)_ang(\d+)_1deg\.fits$")
    for f in sorted(glob.glob(os.path.join(T, "gap1_hestia_G*_ang*_1deg.fits"))):
        m = pat.search(f)
        sig[f"hestia_{m.group(1)}_ang{m.group(2)}"] = FC.load_template(f)[0]
    sig["vvv"] = FC.load_template(os.path.join(T, "gap1_bulge_coleman.fits"))[0]
    sig["f98"] = FC.load_template(os.path.join(T, "gap1_bulge_f98.fits"))[0]
    sig.update(gnfw_templates(os.path.join(T, "gnfw_scan")))
    if a.library == "triaxial":
        keepk = ["vvv", "hestia_G1.1_ang0", "hestia_G2.1_ang90", "gnfw_g1.2",
                 "gnfw_g1.1"]
        sig = {k: sig[k] for k in keepk}
        for f in sorted(glob.glob(os.path.join(T, "triaxial", "triaxial_*.fits"))):
            name = os.path.basename(f)[:-5]
            sig[name] = FC.load_template(f)[0]
    print(f"{len(sig)} signal templates", flush=True)

    t0 = time.time()
    base = FL.fit(counts, exposure, bg, emid, mask)
    lnl0 = base["lnL_total"]
    if a.ps:
        print("  fitted point-source normalisation per bin (want ~1): " +
              " ".join(f"{b['norms']['ps']:.3f}" for b in base["bins"]),
              flush=True)
    print(f"background only: lnL {lnl0:,.2f}  "
          f"converged={all(b['converged'] for b in base['bins'])}  "
          f"({time.time() - t0:.0f}s)", flush=True)

    rows = []
    for i, (name, t) in enumerate(sig.items(), 1):
        t1 = time.time()
        comps = [{**b, "sig": t} for b in bg]
        r = FL.fit(counts, exposure, comps, emid, mask)
        lnl = r["lnL_total"]
        # photons the fit hands to the signal, summed over bins
        ph = 0.0
        for e in range(len(emid)):
            conv = FL.convolve_templates({"sig": t}, emid[e])["sig"]
            ph += r["bins"][e]["norms"]["sig"] * float(
                (exposure[e] * conv)[mask[e].astype(bool)].sum())
        rows.append({"template": name, "lnL": lnl,
                     "TS": max(2.0 * (lnl - lnl0), 0.0),
                     "photons": ph,
                     "converged": all(b["converged"] for b in r["bins"])})
        print(f"  [{i:2d}/{len(sig)}] {name:22s} TS {rows[-1]['TS']:9.1f}  "
              f"conv={rows[-1]['converged']}  ({time.time() - t1:.0f}s)",
              flush=True)

    ref = next(r for r in rows if r["template"] == "vvv")["lnL"]
    for r in rows:
        r["dAIC_vs_vvv"] = -2.0 * (r["lnL"] - ref)     # negative = beats VVV

    with open(os.path.join(a.out, "scan.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    # ---- correlation matrix, what the likelihood actually sees: PSF-convolved,
    # exposure-weighted, inside the mask, at the bin with the most photons
    e = int(np.argmax(counts.reshape(len(emid), -1).sum(1)))
    m = mask[e].astype(bool)
    best_g = min((r for r in rows if r["template"].startswith("gnfw_")
                  and not r["template"].endswith("_th1")), key=lambda r: -r["lnL"])
    keys = ["iem", "iso", "bubbles", "gnfw_g1.2", best_g["template"],
            "hestia_G1.1_ang0", "vvv", "f98"]
    keys = [k for k in dict.fromkeys(keys) if k in {**bg[0], **sig}]
    lib = {**bg[e], **sig}
    conv = FL.convolve_templates({k: lib[k] for k in keys}, emid[e])
    V = np.stack([(exposure[e] * conv[k])[m] for k in keys])
    C = np.corrcoef(V)
    G = V @ V.T
    d = np.sqrt(np.diag(G))
    cond = float(np.linalg.cond(G / np.outer(d, d)))

    summary = {"mask": a.mask, "kept_fraction": minfo["kept_fraction"],
               "point_source_template": bool(a.ps),
               "ps_norm_per_bin": ([b["norms"]["ps"] for b in base["bins"]]
                                   if a.ps else None),
               "lnL_background_only": lnl0,
               "correlation_energy_gev": float(emid[e]),
               "correlation_keys": keys,
               "correlation": C.round(4).tolist(),
               "condition_number_normalised_gram": cond,
               "rows": rows}
    with open(os.path.join(a.out, "scan.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\ncorrelation at {emid[e]:.2f} GeV inside mask {a.mask} "
          f"(condition number {cond:,.0f}):")
    print("            " + " ".join(f"{k[:10]:>10}" for k in keys))
    for k, row in zip(keys, C):
        print(f"  {k[:10]:>10}" + " ".join(f"{v:10.3f}" for v in row))
    print(f"\ntotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
