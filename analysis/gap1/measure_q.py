"""Fit the GCE axis ratio directly, for comparison with Di Mauro (2026).

Section 6.1 reports a crossover between template families, not a shape,
so it cannot be compared with a published number. The q-ladder from the
H0_B test can be: gNFW-rho^2 templates differing only in vertical axis
ratio, fitted to the real counts, giving a profile likelihood in q.

Conventions are reciprocal. Ours is q = vertical / in-plane, so q < 1 is
squashed toward the plane. Di Mauro's r is along-plane / perpendicular.
  q = 1/r,  so r = 1.10 +- 0.05  ->  q = 0.909 +- 0.041

Errors from the likelihood profile at Delta lnL = 0.5, not ladder spacing.
Level 1 only: 4FGL as given, no iterative source finding.
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
          os.path.join(ROOT, "analysis", "level2")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import build_triaxial as TRI      # noqa: E402
import hestia_to_template as HT   # noqa: E402

DIMAURO_R = 1.10
DIMAURO_R_ERR = 0.05


def ladder(qs, cache=os.path.join(ROOT, "templates", "qladder")):
    os.makedirs(cache, exist_ok=True)
    out = {}
    for q in qs:
        f = os.path.join(cache, f"qladder_{q:.2f}.npy")
        if os.path.exists(f):
            t = np.load(f)
        else:
            t = TRI.triaxial(q_plane=1.0, q_vert=float(q))
            np.save(f, t)
        out[float(q)] = RG.normalise(t)
    return out


def profile(lnl):
    """Peak and Delta lnL = 0.5 interval, interpolated on the profile."""
    qs = np.array(sorted(lnl))
    ls = np.array([lnl[q] for q in qs])
    i = int(np.argmax(ls))
    if 0 < i < len(qs) - 1:
        y0, y1, y2 = ls[i - 1], ls[i], ls[i + 1]
        den = y0 - 2 * y1 + y2
        sh = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        q_hat = float(qs[i] + np.clip(sh, -1, 1) * (qs[i + 1] - qs[i]))
        peak = float(y1 - 0.125 * (y0 - y2) ** 2 / den) if den != 0 else float(y1)
    else:
        q_hat, peak = float(qs[i]), float(ls[i])

    def cross(side):
        target = peak - 0.5
        if side < 0:
            idx = range(i, 0, -1)
        else:
            idx = range(i, len(qs) - 1)
        for k in idx:
            k2 = k - 1 if side < 0 else k + 1
            if ls[k] >= target >= ls[k2]:
                w = (ls[k] - target) / (ls[k] - ls[k2])
                return float(qs[k] + w * (qs[k2] - qs[k]))
        return float("nan")

    return q_hat, peak, cross(-1), cross(+1)



def split_iem(iem_plane, mode):
    """Split the IEM into components that each carry a free normalisation,
    so the diffuse model has freedom in shape and not only in scale. A stand-in
    for the ring decomposition published analyses use. Components sum to the
    input."""
    b = (-20.0 + 0.1 / 2.0) + 0.1 * np.arange(400)
    _, BB = np.meshgrid(np.zeros(400), b)
    BB = np.repeat(b[:, None], 400, axis=1)
    if mode == "none":
        return {"iem": iem_plane}
    if mode == "ns":
        # one extra parameter. Sec 6.10's residual is north/south.
        return {"iem_N": iem_plane * (BB >= 0), "iem_S": iem_plane * (BB < 0)}
    LL = np.repeat(((20.0 - 0.1 / 2.0) - 0.1 * np.arange(400))[None, :], 400,
                   axis=0)
    inner = np.abs(LL) <= 10.0
    edges = [0.0, 2.0, 5.0, 10.0, 20.001]
    if mode in ("latout", "latin"):
        # localise the correction. latout frees everywhere except the GCE's
        # longitudes, latin frees only those - the degenerate case, included
        # to bound how much of the loss is degeneracy rather than mismodelling.
        free = inner if mode == "latin" else ~inner
        out = {"iem_rest": iem_plane * (~free)}
        for k in range(len(edges) - 1):
            sel = (np.abs(BB) >= edges[k]) & (np.abs(BB) < edges[k + 1]) & free
            m = iem_plane * sel
            if m.sum() > 0:
                out[f"iem_b{edges[k]:g}-{edges[k+1]:g}"] = m
        return out
    out = {}
    for k in range(len(edges) - 1):
        sel = (np.abs(BB) >= edges[k]) & (np.abs(BB) < edges[k + 1])
        if mode == "lat":
            m = iem_plane * sel
            if m.sum() > 0:
                out[f"iem_b{edges[k]:g}-{edges[k+1]:g}"] = m
        elif mode == "latns":
            for side, ssel in (("N", BB >= 0), ("S", BB < 0)):
                m = iem_plane * (sel & ssel)
                if m.sum() > 0:
                    out[f"iem_b{edges[k]:g}-{edges[k+1]:g}{side}"] = m
        else:
            raise ValueError(mode)
    return out


def run_one(counts, exposure, comps_base, energies, mask, rung):
    lnl, norms = {}, {}
    for q, t in rung.items():
        comps = [{**b, "sig": t} for b in comps_base]
        r = FL.fit(counts, exposure, comps, energies, mask)
        lnl[q] = r["lnL_total"]
        norms[q] = [b["norms"]["sig"] for b in r["bins"]]
    # TS of the signal at the best rung, against no signal at all
    r0 = FL.fit(counts, exposure, [dict(b) for b in comps_base], energies, mask)
    q_hat, peak, lo, hi = profile(lnl)
    return {"lnL": {float(k): float(v) for k, v in lnl.items()},
            "q_hat": q_hat, "q_lo": lo, "q_hi": hi,
            "lnL_peak": peak, "lnL_nosignal": float(r0["lnL_total"]),
            "ts_signal": float(2 * (peak - r0["lnL_total"])),
            "norms_at_best": norms[min(rung, key=lambda q: abs(q - q_hat))]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--mask", default="B")
    ap.add_argument("--ps", action="store_true")
    ap.add_argument("--qmin", type=float, default=0.30)
    ap.add_argument("--qmax", type=float, default=1.40)
    ap.add_argument("--qstep", type=float, default=0.05)
    ap.add_argument("--iem-split", default="none",
                    choices=["none", "ns", "lat", "latns", "latout", "latin"],
                    help="give the diffuse model freedom in shape as well as "
                         "scale: independent normalisations per latitude band "
                         "(lat) or per band and hemisphere (latns)")
    ap.add_argument("--bins", default="1-10,1-3,3-10",
                    help="energy ranges in GeV, comma separated")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    qs = np.round(np.arange(a.qmin, a.qmax + 1e-9, a.qstep), 2)
    rung = ladder(qs)
    print(f"ladder: {len(qs)} rungs, q = {qs[0]:.2f} .. {qs[-1]:.2f}", flush=True)

    # projected axis ratio of each rung, same estimator as sec 6.11
    proj = {}
    for q, t in rung.items():
        proj[q] = float(HT.axis_ratio(t, level_frac=0.05))

    counts0, e_lo0, e_hi0, chdr = FC.read_ccube(a.ccube)
    results = {}
    for spec in a.bins.split(","):
        lo, hi = (float(x) for x in spec.split("-"))
        keep = (e_lo0 >= lo - 1e-6) & (e_hi0 <= hi + 1e-6)
        if not keep.any():
            print(f"  {spec} GeV: no bins, skipped")
            continue
        e_lo, e_hi = e_lo0[keep], e_hi0[keep]
        emid = np.sqrt(e_lo * e_hi)
        counts = counts0[keep]
        exposure, _ = FC.read_expcube(a.expcube, emid, chdr)
        iem, _ = FC.read_iem(a.iem, emid)
        mask, minfo = FC.build_masks(a.catalog, emid, a.mask,
                                     e_lo=e_lo, e_hi=e_hi, iem=iem)
        flat = RG.normalise(np.ones((400, 400)))
        bub, _ = FC.load_template(os.path.join(ROOT, "templates",
                                               "gap1_bubbles_fuzzy.fits"))
        comps = []
        for e in range(len(emid)):
            parts = split_iem(np.asarray(iem[e], float), a.iem_split)
            c = {"iso": flat, "bubbles": bub}
            for k, v in parts.items():
                c[k] = RG.normalise(v)
            comps.append(c)
        if a.ps:
            ps_cube, _ = PS.build(a.catalog, e_lo, e_hi)
            for e in range(len(emid)):
                comps[e]["ps"] = ps_cube[e]
        r = run_one(counts, exposure, comps, emid, mask, rung)
        r["kept_fraction"] = float(minfo["kept_fraction"])
        r["n_bins"] = int(len(emid))
        r["projected_q_at_best"] = float(np.interp(
            r["q_hat"], qs, [proj[q] for q in qs]))
        results[spec] = r
        print(f"\n  {spec} GeV  ({len(emid)} bins, "
              f"{len(comps[0])} background components, mask {a.mask}"
              f"{'+PS' if a.ps else ''}, keeps {minfo['kept_fraction']:.1%})")
        print(f"    q = {r['q_hat']:.3f}  [{r['q_lo']:.3f}, {r['q_hi']:.3f}]"
              f"   (projected {r['projected_q_at_best']:.3f})")
        print(f"    TS of the signal at that q = {r['ts_signal']:.0f}")

    q_dm = 1.0 / DIMAURO_R
    q_dm_err = DIMAURO_R_ERR / DIMAURO_R ** 2
    print("\n================ against Di Mauro (2026) ================")
    print(f"Di Mauro r = {DIMAURO_R} +- {DIMAURO_R_ERR}  (along plane / "
          f"perpendicular)  ->  q = {q_dm:.3f} +- {q_dm_err:.3f}")
    for spec, r in results.items():
        qh = r["projected_q_at_best"]
        halfwidth = 0.5 * abs(r["q_hi"] - r["q_lo"])
        sig = ((q_dm - qh) / np.hypot(max(halfwidth, 1e-3), q_dm_err)
               if np.isfinite(halfwidth) else float("nan"))
        print(f"  {spec:>6} GeV: ours q = {qh:.3f} +- {halfwidth:.3f}, "
              f"Di Mauro {q_dm:.3f} +- {q_dm_err:.3f}  -> {sig:+.1f} sigma apart")
    out = {"convention": "q = vertical / in-plane; Di Mauro r = 1/q",
           "dimauro_r": DIMAURO_R, "dimauro_r_err": DIMAURO_R_ERR,
           "dimauro_q": q_dm, "dimauro_q_err": q_dm_err,
           "mask": a.mask, "ps": bool(a.ps),
           "iem_split": a.iem_split,
           "ladder_q": [float(q) for q in qs],
           "projected_q": {str(q): proj[q] for q in qs},
           "results": results}
    with open(os.path.join(a.out, "measure_q.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwritten {os.path.join(a.out, 'measure_q.json')}")
    return out


if __name__ == "__main__":
    main()
