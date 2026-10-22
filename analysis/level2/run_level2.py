"""Level 2 on the real sky: iterative source finding instead of a fixed mask.

Gap 1 sec 6.0b left one question open. On a gas-decomposed background the
stellar bulge beats every dark-matter template, and Di Mauro gets the opposite
with ring backgrounds too. The difference we have not controlled is that he
lets the fit ADD point sources and we mask a fixed catalogue. Mask A
eliminated itself (it keeps 0.4% of the ROI in the lowest bin), so this is the
last candidate explanation.

There is one design choice here and it decides the answer, so it is made in
the open and both ways are run.

  neutral   Find sources with NO signal template in the model, then fit every
            template against that one common source list. Every template then
            has the same parameter count, so dAIC = -2 dlnL as everywhere
            else. But the iteration has had free rein to absorb the excess
            into point sources first, so this arm is biased AGAINST any
            signal. It is the conservative bound, not the answer.

  joint     Find sources with the candidate template already in the model, one
            run per hypothesis, which is what Di Mauro does. Each hypothesis
            gets its best shot. The cost is that two hypotheses can end up
            with different numbers of sources, so dAIC = -2 dlnL no longer
            holds and the comparison uses the full AIC = 2k - 2lnL with k
            counted per bin.

Reporting only `neutral` would understate every signal; reporting only `joint`
would let a template buy likelihood with extra sources. Both are printed.

  python3 analysis/level2/run_level2.py --ccube ... --expcube ... --iem ...
      --catalog ... --rings <dir> --level rings_species --mask B --ps
      --arm neutral --out outputs/gap1/l2_neutral
"""

from __future__ import annotations

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
for p in (HERE, os.path.join(ROOT, "analysis", "gap1"),
          os.path.join(ROOT, "analysis", "diffuse_freedom"),
          os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import fit_rings as FR            # noqa: E402
import freedom as DF              # noqa: E402
import iterate as L2              # noqa: E402
import project_template as PT     # noqa: E402
import hestia_to_template as H    # noqa: E402

GAMMAS = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]


def n_params(comps, nbins):
    """Free normalisations: one per component per energy bin."""
    return int(sum(len(comps[e]) for e in range(nbins)))


def signal_library(tdir, cache):
    sig = {}
    pat = re.compile(r"gap1_hestia_(G\d\.\d)_ang(\d+)_1deg\.fits$")
    for f in sorted(glob.glob(os.path.join(tdir, "gap1_hestia_G*_ang*_1deg.fits"))):
        m = pat.search(f)
        sig[f"hestia_{m.group(1)}_ang{m.group(2)}"] = FC.load_template(f)[0]
    sig["vvv"] = FC.load_template(os.path.join(tdir, "gap1_bulge_coleman.fits"))[0]
    sig["f98"] = FC.load_template(os.path.join(tdir, "gap1_bulge_f98.fits"))[0]
    os.makedirs(cache, exist_ok=True)
    for g in GAMMAS:
        path = os.path.join(cache, f"gnfw_g{g:.1f}.npy")
        raw = np.load(path) if os.path.exists(path) else PT.normalise(
            PT.analytic(PT.gnfw(g)))
        if not os.path.exists(path):
            np.save(path, raw)
        sig[f"gnfw_g{g:.1f}"] = RG.normalise(raw)
        sig[f"gnfw_g{g:.1f}_th1"] = RG.normalise(H.tophat_smooth(raw, 1.0))
    return sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", default="",
                    help="single diffuse cube; needed unless --components")
    ap.add_argument("--components", default="",
                    help="folder of diffuse cubes, one free norm each (gap 5)")
    ap.add_argument("--label", default="",
                    help="diffuse model name for the output")
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--rings", default="")
    ap.add_argument("--level", default="rings_species")
    ap.add_argument("--mask", default="B", choices=list("ABCDE"))
    ap.add_argument("--ps", action="store_true")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--arm", default="neutral", choices=("neutral", "joint"))
    ap.add_argument("--joint-templates", default="vvv,f98,gnfw_g1.2,"
                                                 "hestia_G3.2_ang90")
    ap.add_argument("--simulate-sky", default="",
                    help="run on a simulated sky containing this template, at "
                         "the flux --flux-template fits, so the recovered "
                         "morphology can be compared with what went in")
    ap.add_argument("--flux-template", default="")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--checkpoint", action="store_true",
                    help="save progress in --out so a killed run resumes")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = os.path.join(ROOT, "templates")

    counts, e_lo, e_hi, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo >= a.emin - 1e-6) & (e_hi <= a.emax + 1e-6)
    counts, e_lo, e_hi = counts[keep], e_lo[keep], e_hi[keep]
    emid = np.sqrt(e_lo * e_hi)
    nb = len(emid)
    exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
    comp_cubes = {}
    if a.components:
        files = sorted(glob.glob(os.path.join(a.components, "*.fits")))
        if not files:
            raise SystemExit(f"no .fits cubes in {a.components}")
        for f in files:
            comp_cubes[os.path.basename(f)[:-5]] = FC.read_iem(f, emid)[0]
    if a.iem:
        iem, _ = FC.read_iem(a.iem, emid)
    elif comp_cubes:
        # only mask E needs this
        iem = sum(comp_cubes.values())
    else:
        raise SystemExit("need --iem or --components")
    mask, minfo = FC.build_masks(a.catalog, emid, a.mask,
                                 e_lo=e_lo, e_hi=e_hi, iem=iem)
    flat = RG.normalise(np.ones(counts.shape[1:]))
    bub, _ = FC.load_template(os.path.join(T, "gap1_bubbles_fuzzy.fits"))

    diffuse = [{"iem": RG.normalise(iem[e])} for e in range(nb)]
    if comp_cubes:
        if a.rings:
            raise SystemExit("--components and --rings are alternatives")
        diffuse = [{f"iem_{k}": RG.normalise(np.maximum(c[e], 0.0))
                    for k, c in comp_cubes.items()} for e in range(nb)]
    if a.rings:
        ladder = FR.levels(FR.load_rings(a.rings, nb, emid), iem, nb)
        diffuse = ladder[a.level]
    bg = [{"iso": flat, "bubbles": bub, **diffuse[e]} for e in range(nb)]
    if a.ps:
        ps, psinfo = PS.build(a.catalog, e_lo, e_hi)
        for e in range(nb):
            bg[e]["ps"] = ps[e]

    model_name = a.label or (os.path.basename(os.path.normpath(a.components))
                             if comp_cubes else
                             (a.level if a.rings else "gll_iem_v07"))
    sig = signal_library(T, os.path.join(T, "gnfw_scan"))
    print(f"mask {a.mask}{'+PS' if a.ps else ''} keeps "
          f"{minfo['kept_fraction']:.1%}; {nb} bins {a.emin:g}-{a.emax:g} GeV; "
          f"background {model_name} "
          f"({len(bg[0])} components); {len(sig)} templates", flush=True)

    # optional simulated sky with a known morphology in it
    truth_info = None
    if a.simulate_sky == "null":
        # A sky with the fitted background and NO signal and no extra point
        # sources. Whatever the iteration adds here it invented: this is the
        # phantom rate, and it is the number Gap 5 wants.
        here = FL.fit(counts, exposure, bg, emid, mask)
        rng = np.random.default_rng(a.seed)
        counts = DF.simulate(counts, exposure, bg, emid,
                             [here["bins"][e]["norms"] for e in range(nb)], rng)
        truth_info = {"template": None, "flux": 0.0}
        print("simulated sky contains NO signal and no sources beyond the "
              "background model\n", flush=True)
    elif a.simulate_sky:
        truth = sig[a.simulate_sky]
        fluxt = sig[a.flux_template] if a.flux_template else truth
        ref = FL.fit(counts, exposure, [{**c, "sig": fluxt} for c in bg],
                     emid, mask)
        scale = float(fluxt.sum()) / float(truth.sum())
        tn = [ref["bins"][e]["norms"]["sig"] * scale for e in range(nb)]
        here = FL.fit(counts, exposure, [{**c, "sig": truth} for c in bg],
                      emid, mask)
        norms = [{**here["bins"][e]["norms"], "sig": tn[e]} for e in range(nb)]
        rng = np.random.default_rng(a.seed)
        counts = DF.simulate(counts, exposure, bg, emid, norms, rng,
                             extra=(truth,))
        truth_info = {"template": a.simulate_sky,
                      "flux": float(sum(tn) * float(truth.sum()))}
        print(f"simulated sky contains {a.simulate_sky} at "
              f"{truth_info['flux']:.4e} ph cm^-2 s^-1\n", flush=True)

    cfg = L2.load_config()
    print(f"stopping rule: passes {cfg['passes']}, "
          f"max {cfg['max_sources_per_pass']}/pass, "
          f"min separation {cfg['min_separation_deg']} deg\n", flush=True)

    t0 = time.time()
    rows, record = [], {"arm": a.arm, "mask": a.mask, "ps": bool(a.ps),
                        "emin": a.emin, "emax": a.emax,
                        "diffuse_model": model_name,
                        "n_diffuse_components": len(diffuse[0]),
                        "kept_fraction": minfo["kept_fraction"],
                        "stopping_rule": cfg, "simulated_sky": truth_info,
                        "seed": a.seed}

    if a.arm == "neutral":
        print("finding sources with NO signal in the model", flush=True)
        ck = os.path.join(a.out, "ckpt_neutral_sources.json") if a.checkpoint else None
        comps, fit0, found = L2.run(counts, exposure, bg, emid, mask, cfg, checkpoint=ck)
        rows_ck = os.path.join(a.out, "ckpt_neutral_rows.json")
        done = json.load(open(rows_ck)) if (a.checkpoint and os.path.exists(rows_ck)) else {}
        k0 = n_params(comps, nb)
        print(f"\n{len(found)} sources added; background lnL "
              f"{fit0['lnL_total']:,.2f}; {k0} free normalisations "
              f"({time.time() - t0:.0f}s)\n", flush=True)
        record["sources"] = [dict(zip(
            ("name", "residual_ts", "l", "b", "pass", "index"), f))
            for f in found]
        for i, (name, t) in enumerate(sig.items(), 1):
            if name in done:
                rows.append(done[name])
                continue
            c = [{**comps[e], "sig": t} for e in range(nb)]
            r = FL.fit(counts, exposure, c, emid, mask)
            rows.append({"template": name, "lnL": r["lnL_total"],
                         "TS": max(2.0 * (r["lnL_total"] - fit0["lnL_total"]), 0.0),
                         "n_sources": len(found),
                         "n_params": n_params(c, nb),
                         "converged": all(b["converged"] for b in r["bins"]),
                         "max_lnL_gap": max(b["lnL_gap_bound"] for b in r["bins"])})
            print(f"  [{i:2d}/{len(sig)}] {name:22s} TS {rows[-1]['TS']:9.1f}",
                  flush=True)
            if a.checkpoint:
                done[name] = rows[-1]
                json.dump(done, open(rows_ck + ".tmp", "w"))
                os.replace(rows_ck + ".tmp", rows_ck)
    else:
        want = [s for s in a.joint_templates.split(",") if s]
        rows_ck = os.path.join(a.out, "ckpt_joint_rows.json")
        done = json.load(open(rows_ck)) if (a.checkpoint and os.path.exists(rows_ck)) else {}
        for i, name in enumerate(want, 1):
            if name in done:
                rows.append(done[name])
                continue
            t1 = time.time()
            print(f"\n[{i}/{len(want)}] finding sources WITH {name} in the model",
                  flush=True)
            bgs = [{**bg[e], "sig": sig[name]} for e in range(nb)]
            ck = (os.path.join(a.out, f"ckpt_joint_{name}_sources.json")
                  if a.checkpoint else None)
            comps, fit, found = L2.run(counts, exposure, bgs, emid, mask, cfg, checkpoint=ck)
            k = n_params(comps, nb)
            rows.append({"template": name, "lnL": fit["lnL_total"],
                         "n_sources": len(found), "n_params": k,
                         "AIC": 2.0 * k - 2.0 * fit["lnL_total"],
                         "converged": all(b["converged"] for b in fit["bins"]),
                         "max_lnL_gap": max(b["lnL_gap_bound"] for b in fit["bins"]),
                         "sources": [dict(zip(
                             ("name", "residual_ts", "l", "b", "pass", "index"), f))
                             for f in found]})
            print(f"  {name}: {len(found)} sources, lnL {fit['lnL_total']:,.2f}, "
                  f"k {k}, AIC {rows[-1]['AIC']:,.2f} "
                  f"({time.time() - t1:.0f}s)", flush=True)
            if a.checkpoint:
                done[name] = rows[-1]
                json.dump(done, open(rows_ck + ".tmp", "w"))
                os.replace(rows_ck + ".tmp", rows_ck)
        best = min(r["AIC"] for r in rows)
        for r in rows:
            r["dAIC_vs_best"] = r["AIC"] - best
        print(f"\n  {'template':22s} {'sources':>8} {'k':>6} {'lnL':>14} "
              f"{'AIC':>14} {'dAIC':>9}")
        for r in sorted(rows, key=lambda r: r["AIC"]):
            print(f"  {r['template']:22s} {r['n_sources']:8d} {r['n_params']:6d} "
                  f"{r['lnL']:14,.2f} {r['AIC']:14,.2f} {r['dAIC_vs_best']:9.1f}")

    if a.arm == "neutral" and rows:
        ref = next(r for r in rows if r["template"] == "vvv")["lnL"]
        for r in rows:
            r["dAIC_vs_vvv"] = -2.0 * (r["lnL"] - ref)
        print(f"\n  {'template':22s} {'TS':>10} {'dAIC vs vvv':>13}")
        for r in sorted(rows, key=lambda r: r["dAIC_vs_vvv"])[:8]:
            print(f"  {r['template']:22s} {r['TS']:10.1f} "
                  f"{r['dAIC_vs_vvv']:13.1f}")

    record["rows"] = rows
    record["wall_seconds"] = time.time() - t0
    with open(os.path.join(a.out, "level2.json"), "w") as f:
        json.dump(record, f, indent=1)
    if a.arm == "neutral":
        # scan.json in the same shape robustness_scan writes, so the Gap 7
        # exporter reads a Level-2 run without knowing it is one
        with open(os.path.join(a.out, "scan.json"), "w") as f:
            json.dump({"mask": a.mask,
                       "point_source_template": bool(a.ps),
                       "kept_fraction": minfo["kept_fraction"],
                       "emin": a.emin, "emax": a.emax,
                       "diffuse_model": model_name,
                       "n_diffuse_components": len(diffuse[0]),
                       "level": 2, "n_sources": len(found),
                       "simulated_sky": truth_info,
                       "admissibility": None, "rows": rows}, f, indent=1)
        cols = ["template", "TS", "lnL", "dAIC_vs_vvv", "n_sources",
                "n_params", "converged", "max_lnL_gap"]
        with open(os.path.join(a.out, "scan.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
    print(f"\nwrote {a.out}  ({time.time() - t0:.0f}s total)")


if __name__ == "__main__":
    main()
