"""Shared diffuse-freedom test.

How much of a reported TS is the diffuse model's error rather than a signal?

Refit with the interstellar model split into components that each carry a free
normalisation, so it has freedom in shape and not only in scale, and watch the
TS. On its own that proves nothing: a split flexible enough to eat a real
signal would give the same answer. So every split is run on three skies.

  observed    the data
  injected    a simulation from the fitted background plus a known signal
  boosted     the data with that known signal added on top

A split is usable as evidence only if `injected` survives it. If `boosted` also
survives and `observed` does not, the observed excess is the background model's
error.

Gap-agnostic: callers pass their own component maps, energies, mask and signal
templates. Splits come from configs/diffuse_freedom.yaml so every gap runs the
same ones, fixed before any run.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (os.path.join(ROOT, "utils"), os.path.join(ROOT, "analysis", "gap1")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402

DEFAULT_CONFIG = os.path.join(ROOT, "configs", "diffuse_freedom.yaml")


def load_config(path=None):
    import yaml
    with open(path or DEFAULT_CONFIG) as f:
        return yaml.safe_load(f)


def grid(nl, nb, half_l, half_b):
    """(l, b) per pixel for a CAR grid with l descending, as Conventions §3."""
    l = (half_l - (2 * half_l / nl) / 2.0) - (2 * half_l / nl) * np.arange(nl)
    b = (-half_b + (2 * half_b / nb) / 2.0) + (2 * half_b / nb) * np.arange(nb)
    return np.repeat(l[None, :], nb, axis=0), np.repeat(b[:, None], nl, axis=1)


def split_iem(iem_plane, spec, LL, BB):
    """Split one diffuse plane into components that sum to it.

    spec is a dict from the config:
      bands      latitude edges in |b|, e.g. [0, 2, 5, 10, 20.001]
      hemisphere true to separate north from south
      lon_gate   null, "inside" or "outside" - restrict the freedom to
                 |l| <= lon_deg or to the rest, leaving the other part on a
                 single normalisation. Use it to ask where the correction is
                 being applied.
      lon_deg    the gate longitude
    """
    bands = spec.get("bands")
    if not bands:
        return {"iem": iem_plane}
    gate = spec.get("lon_gate")
    if gate is None:
        free = np.ones_like(LL, bool)
        out = {}
    else:
        inside = np.abs(LL) <= float(spec.get("lon_deg", 10.0))
        free = inside if gate == "inside" else ~inside
        out = {"iem_fixed": iem_plane * (~free)}
    sides = ((("N", BB >= 0), ("S", BB < 0)) if spec.get("hemisphere")
             else (("", np.ones_like(BB, bool)),))
    for k in range(len(bands) - 1):
        sel = (np.abs(BB) >= bands[k]) & (np.abs(BB) < bands[k + 1]) & free
        for tag, ssel in sides:
            m = iem_plane * (sel & ssel)
            if m.sum() > 0:
                out[f"iem_b{bands[k]:g}-{bands[k+1]:g}{tag}"] = m
    return out


def _ts(counts, exposure, comps, energies, mask, sig):
    r0 = FL.fit(counts, exposure, [dict(c) for c in comps], energies, mask)
    r1 = FL.fit(counts, exposure, [{**c, "sig": sig} for c in comps],
                energies, mask)
    flux = sum(b["norms"]["sig"] for b in r1["bins"]) * float(sig.sum())
    return float(2 * (r1["lnL_total"] - r0["lnL_total"])), float(flux)


def simulate(counts_like, exposure, comps, energies, norms, rng, extra=None):
    """Poisson sky from a fitted model. `extra` adds one more component."""
    out = np.zeros_like(counts_like, dtype=float)
    for e in range(len(energies)):
        comp = dict(comps[e])
        if extra is not None:
            comp["sig"] = extra[0]
        conv = FL.convolve_templates(comp, energies[e])
        mu = sum(norms[e].get(k, 0.0) * conv[k] for k in comp)
        out[e] = rng.poisson(np.maximum(exposure[e] * mu, 0))
    return out


def run(counts, exposure, comps_base, energies, mask, signals, iem_planes,
        config=None, control_signal=None, seed=None, verbose=True):
    """Run the test.

    comps_base  per-bin dicts of background components EXCLUDING the diffuse
                plane, which is passed separately as iem_planes[e] so it can
                be split
    signals     {name: map} templates to test
    returns     {"splits": {name: {sky: {signal: {ts, flux}}}}, ...}
    """
    cfg = config or load_config()
    seed = int(cfg.get("seed", 42)) if seed is None else seed
    nb, nl = iem_planes[0].shape
    LL, BB = grid(nl, nb, cfg["roi"]["half_l_deg"], cfg["roi"]["half_b_deg"])
    ctrl = control_signal or cfg["control"]["signal"]
    if ctrl not in signals:
        raise SystemExit(f"control signal {ctrl!r} is not in {list(signals)}")

    def build(spec):
        out = []
        for e in range(len(energies)):
            c = dict(comps_base[e])
            for k, v in split_iem(np.asarray(iem_planes[e], float),
                                  spec, LL, BB).items():
                c[k] = RG.normalise(v)
            out.append(c)
        return out

    specs = cfg["splits"]
    rigid = build(specs[0])
    ref = FL.fit(counts, exposure,
                 [{**c, "sig": signals[ctrl]} for c in rigid], energies, mask)
    rng = np.random.default_rng(seed)
    norms = [ref["bins"][e]["norms"] for e in range(len(energies))]
    sim = simulate(counts, exposure, rigid, energies, norms, rng,
                   extra=(signals[ctrl],))
    add = simulate(counts, exposure, [{} for _ in energies], energies,
                   [{"sig": n["sig"]} for n in norms], rng,
                   extra=(signals[ctrl],))
    skies = {"observed": counts, "injected": sim, "boosted": counts + add}

    res = {}
    for spec in specs:
        comps = build(spec)
        res[spec["name"]] = {"n_components": len(comps[0]),
                             "extra_params": len(comps[0]) - len(rigid[0])}
        for sky_name, sky in skies.items():
            res[spec["name"]][sky_name] = {}
            for sname, sig in signals.items():
                ts, flux = _ts(sky, exposure, comps, energies, mask, sig)
                res[spec["name"]][sky_name][sname] = {"ts": ts, "flux": flux}
        if verbose:
            r = res[spec["name"]]
            print(f"  {spec['name']:10} ({r['n_components']:2d} comps): " +
                  "  ".join(f"{n} {r['observed'][n]['ts']:7.1f}"
                            for n in signals), flush=True)
    return {"control_signal": ctrl, "seed": seed,
            "control_flux": float(sum(n["sig"] for n in norms)
                                  * float(signals[ctrl].sum())),
            "splits": res}


def verdict(res, signal=None, usable_injected_frac=0.5, boosted_min_gain=25.0):
    """Summarise. A split counts as usable only if it keeps
    `usable_injected_frac` of the injected signal."""
    s = signal or res["control_signal"]
    sp = res["splits"]
    rigid = list(sp)[0]
    base_obs = sp[rigid]["observed"][s]["ts"]
    base_inj = sp[rigid]["injected"][s]["ts"]
    rows, usable = [], []
    for name, r in sp.items():
        obs = r["observed"][s]["ts"]
        inj = r["injected"][s]["ts"]
        gain = r["boosted"][s]["ts"] - obs
        keep_inj = inj / base_inj if base_inj > 0 else float("nan")
        ok = keep_inj >= usable_injected_frac and gain >= boosted_min_gain
        rows.append({"split": name, "extra_params": r["extra_params"],
                     "observed_ts": obs, "injected_ts": inj,
                     "boosted_gain": gain,
                     "injected_fraction_kept": keep_inj, "usable": ok})
        if ok:
            usable.append(obs)
    if len(usable) > 1:
        spread = max(usable) / max(min(usable), 1e-9)
        note = ("TS spans a factor of %.0f across splits that a control shows "
                "are not degenerate with the signal" % spread)
    elif len(usable) == 1:
        spread = 1.0
        note = ("only one split is admissible, so the spread is untested here "
                "rather than small")
    else:
        spread = float("nan")
        note = ("no split is admissible: this configuration cannot recover a "
                "signal injected at the observed flux, so it cannot measure "
                "a TS either way")
    return {"signal": s, "rigid_split": rigid, "rigid_observed_ts": base_obs,
            "rows": rows, "usable_observed_ts": usable,
            "spread_factor": spread, "note": note}
