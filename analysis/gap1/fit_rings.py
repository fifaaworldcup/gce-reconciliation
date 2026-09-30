"""Redo the sec 6.0 test with Pohl22's ring decomposition instead of one
all-sky template.

Sec 6.0 showed the detection runs from TS 3.8 to 306 depending on how much a
single interstellar template is allowed to bend, and said the fix is a diffuse
model whose freedom comes from gas astronomy rather than from our choice of
latitude bands. Pohl22 (Zenodo 6276721) supplies that: four HI rings, four H2
rings, six inverse-Compton rings and a dust residual map, all on the sky.

The freedom ladder here is physical rather than geometric:

  gll_iem_v07   our old background, one template, one normalisation per bin
  rings_tied    the same rings summed to one template - Pohl's shape, our
                rigidity. Isolates shape from freedom.
  rings_species HI, H2, IC and dust each free
  rings_free    all sixteen components free

Same three skies as the shared module: observed, a simulation containing a
known signal, and the real data with that signal added on top.

Dust_Negative_Residuals is not used. The fitter constrains normalisations to
be non-negative, so a map meant to be subtracted cannot be given its proper
sign; including it as a positive component would let it soak up emission it
should be removing.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks"),
          os.path.join(ROOT, "analysis", "diffuse_freedom")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import freedom as DF              # noqa: E402

TEMPLATES = {"gnfw1.2": "gnfw_gamma1.2.fits",
             "hestia": "gap1_hestia_G1.1_ang0_1deg.fits",
             "vvv": "gap1_bulge_coleman.fits"}
SPECIES = {"hi": "pohl22_hi_ring", "h2": "pohl22_h2_ring",
           "ic": "pohl22_ic_ring"}
EXTRA = {"dust": "pohl22_dust_pos.fits", "loopI": "pohl22_loopI.fits"}
# not part of the diffuse model - added only in the NB / NB+BB levels
BULGES = {"nb": "pohl22_nb.fits", "bb": "pohl22_bb.fits"}


def load_rings(d, nbins):
    """{name: [plane per energy bin]} for every Pohl22 component."""
    out = {}
    for tag, stem in SPECIES.items():
        for f in sorted(glob.glob(os.path.join(d, stem + "*.fits"))):
            a = np.asarray(fits.getdata(f), float)
            name = os.path.basename(f)[len("pohl22_"):-len(".fits")]
            if a.ndim == 3:
                if a.shape[0] != nbins:
                    raise SystemExit(f"{name} has {a.shape[0]} planes, "
                                     f"cube has {nbins} bins")
                out[name] = [a[e] for e in range(nbins)]
            else:
                out[name] = [a] * nbins
    for tag, fn in list(EXTRA.items()) + list(BULGES.items()):
        f = os.path.join(d, fn)
        if os.path.exists(f):
            out[tag] = [np.asarray(fits.getdata(f), float)] * nbins
    if not out:
        raise SystemExit(f"no pohl22_*.fits under {d}")
    return out


def levels(rings, iem, nbins):
    """The freedom ladder: {level: [per-bin dict of diffuse components]}."""
    def norm(m):
        return RG.normalise(np.maximum(m, 0.0))

    lv = {}
    lv["gll_iem_v07"] = [{"iem": norm(iem[e])} for e in range(nbins)]
    diffuse_only = {k: v for k, v in rings.items() if k not in BULGES}
    lv["rings_species"] = []
    for e in range(nbins):
        d = {}
        for tag in list(SPECIES) + list(EXTRA):
            parts = [v[e] for k, v in diffuse_only.items()
                     if k.startswith(tag)]
            if parts:
                d[tag] = norm(sum(parts))
        lv["rings_species"].append(d)
    lv["rings_free"] = [{k: norm(v[e]) for k, v in diffuse_only.items()}
                        for e in range(nbins)]

    # Macias's configuration, with Pohl's own bulge templates. Macias reports
    # NFW^2 falling from 17.6 sigma to 2.4 once the nuclear bulge and the
    # X/boxy bulge are in the background; Di Mauro reports dark matter holding
    # at TS 8e2 to 1.5e4 with the same two included. Same templates, same gas
    # maps, opposite conclusions - so run it and see which way ours falls.
    nb = [k for k in ("nb",) if k in rings]
    bb = [k for k in ("bb",) if k in rings]
    if nb:
        lv["rings_free_NB"] = [{**lv["rings_free"][e],
                                **{k: norm(rings[k][e]) for k in nb}}
                               for e in range(nbins)]
    if nb and bb:
        lv["rings_free_NB_BB"] = [{**lv["rings_free"][e],
                                   **{k: norm(rings[k][e]) for k in nb + bb}}
                                  for e in range(nbins)]
    return lv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--rings", required=True, help="dir of pohl22_*.fits")
    ap.add_argument("--mask", default="B", choices=list("ABCDE"))
    ap.add_argument("--ps", action="store_true")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--control", default="vvv")
    ap.add_argument("--inject-flux-from", default="gll_iem_v07",
                    help="level whose fitted signal normalisation sets the "
                         "injected flux. Injecting at each level's OWN fitted "
                         "flux is vacuous where that flux is near zero: it "
                         "only shows the fit is self-consistent. Holding the "
                         "flux fixed asks the question that matters - could "
                         "this background have seen a signal of the size the "
                         "baseline found?")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    counts0, e_lo0, e_hi0, chdr = FC.read_ccube(a.ccube)
    keep = (e_lo0 >= a.emin - 1e-6) & (e_hi0 <= a.emax + 1e-6)
    e_lo, e_hi = e_lo0[keep], e_hi0[keep]
    emid = np.sqrt(e_lo * e_hi)
    counts = counts0[keep]
    nb = len(emid)
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
    base = [{"iso": flat, "bubbles": bub} for _ in range(nb)]
    if a.ps:
        ps_cube, _ = PS.build(a.catalog, e_lo, e_hi)
        for e in range(nb):
            base[e]["ps"] = ps_cube[e]

    rings = load_rings(a.rings, nb)
    print(f"mask {a.mask}{'+PS' if a.ps else ''} keeps "
          f"{minfo['kept_fraction']:.1%}; {nb} bins {a.emin:g}-{a.emax:g} GeV; "
          f"{counts.sum():,.0f} photons")
    print(f"pohl22 components: {', '.join(sorted(rings))}\n")

    lv = levels(rings, iem, nb)

    # reference signal normalisation, common to every level's injection
    ref_lv = a.inject_flux_from
    if ref_lv not in lv:
        raise SystemExit(f"--inject-flux-from {ref_lv} is not a level")
    ref_comps = [{**base[e], **lv[ref_lv][e]} for e in range(nb)]
    ref_fit = FL.fit(counts, exposure,
                     [{**c, "sig": signals[a.control]} for c in ref_comps],
                     emid, mask)
    sig_norms = [ref_fit["bins"][e]["norms"]["sig"] for e in range(nb)]
    inj_flux = float(sum(sig_norms) * float(signals[a.control].sum()))
    print(f"injecting {a.control} at {inj_flux:.4e} ph cm^-2 s^-1, "
          f"the flux the {ref_lv} fit found\n")

    res = {}
    for name, diff in lv.items():
        comps = [{**base[e], **diff[e]} for e in range(nb)]
        # Each level gets its OWN simulated skies, built from its own fit.
        # A sky simulated from one background and fitted with another carries
        # the mismatch between them, which the signal template then absorbs -
        # the control has to be internal to the model under test.
        ref = FL.fit(counts, exposure,
                     [{**c, "sig": signals[a.control]} for c in comps],
                     emid, mask)
        # this level's own background, but the signal at the common reference
        # flux rather than whatever this level happened to fit
        norms = [{**ref["bins"][e]["norms"], "sig": sig_norms[e]}
                 for e in range(nb)]
        rng = np.random.default_rng(a.seed)
        sim = DF.simulate(counts, exposure, comps, emid, norms, rng,
                          extra=(signals[a.control],))
        add = DF.simulate(counts, exposure, [{} for _ in range(nb)], emid,
                          [{"sig": sig_norms[e]} for e in range(nb)], rng,
                          extra=(signals[a.control],))
        skies = {"observed": counts, "injected": sim, "boosted": counts + add}
        res[name] = {"n_components": len(comps[0]),
                     "n_diffuse": len(diff[0]),
                     "control_flux": inj_flux,
                     "own_fitted_flux": float(
                         sum(ref["bins"][e]["norms"]["sig"] for e in range(nb))
                         * float(signals[a.control].sum()))}
        for sky_name, sky in skies.items():
            res[name][sky_name] = {}
            for sname, sig in signals.items():
                ts, flux = DF._ts(sky, exposure, comps, emid, mask, sig)
                res[name][sky_name][sname] = {"ts": ts, "flux": flux}
        r = res[name]
        print(f"  {name:14} ({r['n_diffuse']:2d} diffuse comps): " +
              "  ".join(f"{n} {r['observed'][n]['ts']:8.1f}" for n in signals),
              flush=True)

    s = a.control
    print(f"\n{'level':14} {'diffuse':>8} {'observed':>9} {'injected':>9} "
          f"{'boosted+':>9} {'inj/obs':>9}  can it see a signal?")
    rows = []
    for name, r in res.items():
        obs = r["observed"][s]["ts"]
        inj = r["injected"][s]["ts"]
        gain = r["boosted"][s]["ts"] - obs
        # each level's injected sky contains the same flux it itself fitted,
        # so the ratio to that level's own observed TS is the like-for-like
        # check on whether the model can see a signal it is given
        # the injection is at a fixed flux for every level, so the question
        # is simply whether the level can see it
        kept = inj / max(obs, 1e-9)
        ok = inj >= 25.0 and gain >= 25.0
        rows.append({"level": name, "n_diffuse": r["n_diffuse"],
                     "observed_ts": obs, "injected_ts": inj,
                     "boosted_gain": gain, "injected_fraction_kept": kept,
                     "usable": ok})
        print(f"{name:14} {r['n_diffuse']:8d} {obs:9.1f} {inj:9.1f} "
              f"{gain:9.1f} {kept:8.0%}  "
              f"{'yes' if ok else 'NO'}")
    use = [r["observed_ts"] for r in rows if r["usable"]]
    spread = max(use) / max(min(use), 1e-9) if len(use) > 1 else float("nan")
    print(f"\nacross admissible levels the observed TS spans a factor of "
          f"{spread:.0f}")
    print("sec 6.0, with one all-sky template and geometric splits, gave 80")
    with open(os.path.join(a.out, "fit_rings.json"), "w") as f:
        json.dump({"mask": a.mask, "ps": bool(a.ps), "control": s,
                   "seed": a.seed, "kept_fraction": float(minfo["kept_fraction"]),
           "inject_flux_from": ref_lv, "injected_flux": inj_flux,
                   "levels": res, "summary": rows,
                   "spread_factor": spread}, f, indent=2)
    print(f"written {os.path.join(a.out, 'fit_rings.json')}")


if __name__ == "__main__":
    main()
