"""CLI for the shared diffuse-freedom test, Gap 1's inputs.

Other gaps should call freedom.run() with their own components rather than
extend this driver - see analysis/diffuse_freedom/README.md.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks"),
          os.path.join(ROOT, "analysis", "gap1")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import freedom as DF              # noqa: E402

TEMPLATES = {
    "gnfw1.2": "gnfw_gamma1.2.fits",
    "hestia": "gap1_hestia_G1.1_ang0_1deg.fits",
    "vvv": "gap1_bulge_coleman.fits",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--mask", default="B", choices=list("ABCDE"))
    ap.add_argument("--ps", action="store_true")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cfg = DF.load_config(a.config)

    counts0, e_lo0, e_hi0, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo0 >= a.emin - 1e-6) & (e_hi0 <= a.emax + 1e-6)
    e_lo, e_hi = e_lo0[keep], e_hi0[keep]
    if not len(e_lo):
        raise SystemExit(f"no bins inside {a.emin}-{a.emax} GeV")
    emid = np.sqrt(e_lo * e_hi)
    counts = counts0[keep]
    exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
    iem, _ = FC.read_iem(a.iem, emid)
    mask, minfo = FC.build_masks(a.catalog, emid, a.mask,
                                 e_lo=e_lo, e_hi=e_hi, iem=iem)
    flat = RG.normalise(np.ones(counts.shape[1:]))
    bub, _ = FC.load_template(os.path.join(ROOT, "templates",
                                           "gap1_bubbles_fuzzy.fits"))
    T = os.path.join(ROOT, "templates")
    signals = {k: FC.load_template(os.path.join(T, v))[0]
               for k, v in TEMPLATES.items()}
    comps_base = [{"iso": flat, "bubbles": bub} for _ in emid]
    if a.ps:
        ps_cube, _ = PS.build(a.catalog, e_lo, e_hi)
        for e in range(len(emid)):
            comps_base[e]["ps"] = ps_cube[e]

    print(f"mask {a.mask}{'+PS' if a.ps else ''} keeps "
          f"{minfo['kept_fraction']:.1%}; {len(emid)} bins "
          f"{a.emin:g}-{a.emax:g} GeV; {counts.sum():,.0f} photons\n"
          f"observed TS per split:")
    res = DF.run(counts, exposure, comps_base, emid, mask, signals,
                 [iem[e] for e in range(len(emid))], config=cfg)
    v = DF.verdict(res,
                   usable_injected_frac=cfg["control"]["usable_injected_fraction"],
                   boosted_min_gain=cfg["control"]["boosted_min_gain"])

    print(f"\n{'split':12} {'extra':>6} {'observed':>9} {'injected':>9} "
          f"{'boosted+':>9} {'inj kept':>9}  usable")
    for r in v["rows"]:
        print(f"{r['split']:12} {r['extra_params']:6d} {r['observed_ts']:9.1f} "
              f"{r['injected_ts']:9.1f} {r['boosted_gain']:9.1f} "
              f"{r['injected_fraction_kept']:8.0%}  "
              f"{'yes' if r['usable'] else 'no - degenerate'}")
    print(f"\n{v['note']}")
    out = {"mask": a.mask, "ps": bool(a.ps), "emin": a.emin, "emax": a.emax,
           "kept_fraction": float(minfo["kept_fraction"]),
           "photons": float(counts.sum()), **res, "verdict": v}
    with open(os.path.join(a.out, "diffuse_freedom.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"written {os.path.join(a.out, 'diffuse_freedom.json')}")


if __name__ == "__main__":
    main()
