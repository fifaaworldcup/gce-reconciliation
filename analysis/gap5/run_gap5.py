"""
gap 5 driver. per model / mask / energy range, under outputs/gap5/:
  level1    background fit + each template's TS
  admit     injection test from configs/gap5.yaml
  q         q ladder, compared to di mauro
  l2        run_level2.py, joint and neutral
  controls  run_level2.py on simulated skies with a known template
each step skips itself if its output exists, so it can be restarted.

  python analysis/gap5/run_gap5.py --ccube ... --expcube ... --catalog ... \
      --models refs,family2 --steps level1,admit,q,l2
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "analysis", "gap1"),
          os.path.join(ROOT, "analysis", "level2"),
          os.path.join(ROOT, "analysis", "diffuse_freedom"),
          os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC            # noqa: E402
import fit_likelihood as FL       # noqa: E402
import regrid as RG               # noqa: E402
import ps_template as PS          # noqa: E402
import freedom as DF              # noqa: E402
import measure_q as MQ            # noqa: E402
import run_level2 as RL           # noqa: E402

CONFIG = os.path.join(ROOT, "configs", "gap5.yaml")
OUT = os.path.join(ROOT, "outputs", "gap5")
STEPS = ("level1", "admit", "q", "l2", "controls")


def load_config(path=CONFIG):
    with open(path) as f:
        return yaml.safe_load(f)


def model_dirs(spec, cfg, root):
    out = []
    roman = cfg["selection"]["roman"]
    for tok in [t.strip() for t in spec.split(",") if t.strip()]:
        if tok == "refs":
            out += list(cfg["references"])
        elif tok == "family2":
            out += [f"cholis_{roman[i]}" for i in cfg["selection"]["ids"]]
        elif tok in ("reselected", "original"):
            out += [f"cholis_{r}" for r in cfg["selection"][tok]]
        elif tok == "all":
            seen = cfg["selection"]["reselected"] + cfg["selection"]["original"]
            out += [f"cholis_{r}" for r in dict.fromkeys(seen)]
        elif tok == "anchors":
            to_r = {1: "I", 2: "II", 3: "III", 6: "VI", 10: "X", 15: "XV"}
            out += [f"cholis_{to_r[i]}" for i in cfg["selection"]["anchors"]]
        else:
            out.append(tok)
    missing = [m for m in out if not os.path.isdir(os.path.join(root, m))]
    if missing:
        print(f"no component folder yet for: {', '.join(missing)} (skipped)")
    return [m for m in out if m not in missing]


class Setup:

    def __init__(self, a, comp_root, model, mask_name):
        counts, e_lo, e_hi, chdr = FC.read_ccube(a.ccube)
        keep = (e_lo >= a.emin - 1e-6) & (e_hi <= a.emax + 1e-6)
        if not keep.any():
            raise SystemExit(f"no bins of {a.ccube} inside {a.emin}-{a.emax} GeV")
        self.counts, self.e_lo, self.e_hi = counts[keep], e_lo[keep], e_hi[keep]
        self.emid = np.sqrt(self.e_lo * self.e_hi)
        nb = self.nb = len(self.emid)
        self.exposure, _ = FC.read_expcube(a.expcube, self.emid, chdr)
        self.mask, minfo = FC.build_masks(a.catalog, self.emid, mask_name)
        self.kept = minfo["kept_fraction"]
        d = os.path.join(comp_root, model)
        files = sorted(f for f in os.listdir(d) if f.endswith(".fits"))
        cubes = {f[:-5]: FC.read_iem(os.path.join(d, f), self.emid)[0] for f in files}
        flat = RG.normalise(np.ones((400, 400)))
        bub, _ = FC.load_template(os.path.join(ROOT, "templates", "gap1_bubbles_fuzzy.fits"))
        ps, _ = PS.build(a.catalog, self.e_lo, self.e_hi)
        self.bg = [{"iso": flat, "bubbles": bub, "ps": ps[e],
                    **{f"iem_{k}": RG.normalise(np.maximum(c[e], 0.0))
                       for k, c in cubes.items()}} for e in range(nb)]
        self.n_components = len(cubes)
        self.sig = RL.signal_library(os.path.join(ROOT, "templates"),
                                     os.path.join(ROOT, "templates", "gnfw_scan"))


def step_level1(S, templates):
    t0 = time.time()
    r0 = FL.fit(S.counts, S.exposure, S.bg, S.emid, S.mask)
    k0 = RL.n_params(S.bg, S.nb)
    rows = []
    for name in templates:
        t = S.sig[name]
        r = FL.fit(S.counts, S.exposure, [{**c, "sig": t} for c in S.bg], S.emid, S.mask)
        flux = sum(b["norms"]["sig"] for b in r["bins"]) * float(t.sum())
        rows.append({"template": name, "lnL": r["lnL_total"],
                     "TS": max(2 * (r["lnL_total"] - r0["lnL_total"]), 0.0),
                     "flux": flux, "n_params": k0 + S.nb,
                     "converged": all(b["converged"] for b in r["bins"])})
    ref = next((r["lnL"] for r in rows if r["template"] == "vvv"), None)
    for r in rows:
        r["dAIC_vs_vvv"] = None if ref is None else -2 * (r["lnL"] - ref)
    return {"lnL_background": r0["lnL_total"], "n_params_background": k0,
            "AIC_background": 2 * k0 - 2 * r0["lnL_total"],
            "background_converged": all(b["converged"] for b in r0["bins"]),
            "rows": rows, "seconds": time.time() - t0}


def step_admit(S, cfg, ref_flux, seed):
    # same test as gap 1's ring work
    A = cfg["admissibility"]
    out = {"flux_from": A["flux_from"], "controls": {}}
    for ctrl_name in A["controls"]:
        ctrl = S.sig[ctrl_name]
        flux_total = ref_flux[ctrl_name]["flux"]
        per_bin = np.asarray(ref_flux[ctrl_name]["norm_per_bin"], float)
        here = FL.fit(S.counts, S.exposure, [{**c, "sig": ctrl} for c in S.bg],
                      S.emid, S.mask)
        norms = [{**here["bins"][e]["norms"], "sig": per_bin[e]} for e in range(S.nb)]
        rng = np.random.default_rng(seed)
        simulated = DF.simulate(S.counts, S.exposure, S.bg, S.emid, norms, rng,
                                extra=(ctrl,))
        added = DF.simulate(S.counts, S.exposure, [{} for _ in range(S.nb)], S.emid,
                            [{"sig": per_bin[e]} for e in range(S.nb)], rng,
                            extra=(ctrl,))
        ts_obs, _ = DF._ts(S.counts, S.exposure, S.bg, S.emid, S.mask, ctrl)
        ts_inj, _ = DF._ts(simulated, S.exposure, S.bg, S.emid, S.mask, ctrl)
        ts_boost, _ = DF._ts(S.counts + added, S.exposure, S.bg, S.emid, S.mask, ctrl)
        ok = ts_inj >= A["ts_injected_min"] and (ts_boost - ts_obs) >= A["boosted_gain_min"]
        out["controls"][ctrl_name] = {
            "injected_flux": flux_total, "ts_observed": ts_obs,
            "ts_injected": ts_inj, "ts_boosted": ts_boost,
            "boosted_gain": ts_boost - ts_obs, "admissible": bool(ok)}
        print(f"    admissibility {ctrl_name:10s} injected TS {ts_inj:8.1f}  "
              f"gain {ts_boost - ts_obs:+8.1f}  -> {'ok' if ok else 'FAILS'}", flush=True)
    out["admissible"] = all(v["admissible"] for v in out["controls"].values())
    return out


def reference_flux(a, comp_root, cfg, mask_name, outdir):
    path = os.path.join(outdir, "admit_reference_flux.json")
    if os.path.exists(path):
        return json.load(open(path))
    S = Setup(a, comp_root, cfg["admissibility"]["flux_from"], mask_name)
    res = {}
    for name in cfg["admissibility"]["controls"]:
        t = S.sig[name]
        r = FL.fit(S.counts, S.exposure, [{**c, "sig": t} for c in S.bg], S.emid, S.mask)
        nb_ = [b["norms"]["sig"] for b in r["bins"]]
        res[name] = {"norm_per_bin": nb_, "flux": float(sum(nb_) * float(t.sum()))}
    os.makedirs(outdir, exist_ok=True)
    json.dump(res, open(path, "w"), indent=1)
    return res


def q_interval(lnl):
    # measure_q.profile looks for the dlnL = 0.5 crossing between rungs. when
    # the profile is sharper than one rung it finds none, so fall back to the
    # curvature of the parabola through the peak and its two neighbours
    qs = np.array(sorted(lnl, key=float), float)
    ls = np.array([lnl[k] for k in sorted(lnl, key=float)], float)
    i = int(np.argmax(ls))
    if not 0 < i < len(qs) - 1:
        return float("nan"), float("nan")
    h = qs[i + 1] - qs[i]
    curv = (ls[i - 1] - 2 * ls[i] + ls[i + 1]) / h ** 2
    if curv >= 0:
        return float("nan"), float("nan")
    q0 = qs[i] + h * 0.5 * (ls[i - 1] - ls[i + 1]) / (ls[i - 1] - 2 * ls[i] + ls[i + 1])
    w = float(np.sqrt(1.0 / -curv))
    return float(q0 - w), float(q0 + w)


def finish_q(r, qs, proj):
    if not (np.isfinite(r["q_lo"]) and np.isfinite(r["q_hi"])):
        lo, hi = q_interval(r["lnL"])
        r["q_lo"], r["q_hi"], r["interval_from"] = lo, hi, "curvature"
    else:
        r["interval_from"] = "crossing"
    r["projected_q_at_best"] = float(np.interp(r["q_hat"], qs, proj))
    q_dm = 1.0 / MQ.DIMAURO_R
    q_dm_err = MQ.DIMAURO_R_ERR / MQ.DIMAURO_R ** 2
    half = 0.5 * abs(r["q_hi"] - r["q_lo"])
    r["dimauro_q"], r["dimauro_q_err"] = q_dm, q_dm_err
    r["sigma_from_dimauro"] = (float((q_dm - r["projected_q_at_best"])
                                     / np.hypot(max(half, 1e-3), q_dm_err))
                               if np.isfinite(half) else None)
    r["at_ladder_edge"] = bool(r["q_hat"] <= qs[0] + 1e-9 or r["q_hat"] >= qs[-1] - 1e-9)
    return r


def step_q(S, cfg):
    L = cfg["ladder"]
    qs = np.round(np.arange(L["q_min"], L["q_max"] + 1e-9, L["q_step"]), 2)
    rung = MQ.ladder(qs)
    r = MQ.run_one(S.counts, S.exposure, S.bg, S.emid, S.mask, rung)
    import hestia_to_template as HT
    proj = [float(HT.axis_ratio(rung[q], level_frac=0.05)) for q in qs]
    r["ladder_q"] = [float(q) for q in qs]
    r["projected_q"] = proj
    return finish_q(r, qs, proj)


def joint_list(cfg, neutral_json, model=""):
    # gap 1's list (or the short one for cholis models) plus the best of each
    # pattern in this model's neutral arm
    import fnmatch
    short = cfg["analysis"].get("joint_templates_family2")
    tl = list(short if (short and model.startswith("cholis_"))
              else cfg["analysis"]["joint_templates"])
    pats = cfg["analysis"].get("joint_add_best_of_neutral", [])
    if pats and not os.path.exists(neutral_json):
        raise SystemExit(f"joint arm needs the neutral arm first: {neutral_json}")
    rows = json.load(open(neutral_json))["rows"] if pats else []
    for pat in pats:
        m = [r for r in rows if fnmatch.fnmatch(r["template"], pat)
             and not r["template"].endswith("_th1")]
        if m:
            best = max(m, key=lambda r: r["lnL"])["template"]
            if best not in tl:
                tl.append(best)
    return tl


def run_l2(a, comp_root, model, mask_name, arm, outdir, cfg, extra=()):
    cmd = [sys.executable, os.path.join(ROOT, "analysis", "level2", "run_level2.py"),
           "--ccube", a.ccube, "--expcube", a.expcube, "--catalog", a.catalog,
           "--components", os.path.join(comp_root, model), "--label", model,
           "--mask", mask_name, "--ps", "--arm", arm,
           "--emin", str(a.emin), "--emax", str(a.emax),
           "--seed", str(cfg["seed"]), "--checkpoint", *extra]
    if "--joint-templates" not in extra:
        cmd += ["--joint-templates", ",".join(cfg["analysis"]["joint_templates"])]
    cmd += ["--out", outdir]
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "run.log"), "a") as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        print(f"    run_level2 failed ({rc}); see {outdir}/run.log", flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--expcube", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--models", default="refs,family2")
    ap.add_argument("--masks", default="", help="default: analysis.masks in the config")
    ap.add_argument("--steps", default="level1,admit,q,l2")
    ap.add_argument("--emin", type=float, default=1.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--components", default=os.path.join(OUT, "components"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--config", default=CONFIG)
    a = ap.parse_args()
    cfg = load_config(a.config)
    steps = [s for s in a.steps.split(",") if s]
    bad = set(steps) - set(STEPS)
    if bad:
        ap.error(f"unknown steps {bad}")
    masks = a.masks.split(",") if a.masks else cfg["analysis"]["masks"]
    models = model_dirs(a.models, cfg, a.components)
    tag = f"{a.emin:g}-{a.emax:g}"
    templates = cfg["analysis"]["joint_templates"]
    print(f"{len(models)} models x masks {masks}, {tag} GeV, steps {steps}", flush=True)

    for mask_name in masks:
        refdir = os.path.join(a.out, "_reference", mask_name, tag)
        ref_flux = (reference_flux(a, a.components, cfg, mask_name, refdir)
                    if "admit" in steps else None)
        for model in models:
            od = os.path.join(a.out, model, mask_name, tag)
            os.makedirs(od, exist_ok=True)
            print(f"\n=== {model}  mask {mask_name}  {tag} GeV ===", flush=True)
            S = None

            def setup():
                nonlocal S
                if S is None:
                    S = Setup(a, a.components, model, mask_name)
                    print(f"    {S.n_components} diffuse components, mask keeps "
                          f"{S.kept:.1%}, {S.nb} bins", flush=True)
                return S

            jobs = {"level1": lambda: step_level1(setup(), templates),
                    "admit": lambda: step_admit(setup(), cfg, ref_flux, cfg["seed"]),
                    "q": lambda: step_q(setup(), cfg)}
            for st in ("level1", "admit", "q"):
                if st not in steps:
                    continue
                path = os.path.join(od, f"{st}.json")
                if os.path.exists(path):
                    print(f"    {st}: done already", flush=True)
                    continue
                t0 = time.time()
                res = jobs[st]()
                res["model"], res["mask"], res["energy"] = model, mask_name, [a.emin, a.emax]
                json.dump(res, open(path, "w"), indent=1)
                print(f"    {st}: written ({time.time() - t0:.0f}s)", flush=True)

            if "l2" in steps:
                for arm in cfg["analysis"]["arms"]:
                    d = os.path.join(od, f"l2_{arm}")
                    if os.path.exists(os.path.join(d, "level2.json")):
                        print(f"    l2 {arm}: done already", flush=True)
                        continue
                    extra = ()
                    if arm == "joint":
                        tl = joint_list(cfg, os.path.join(od, "l2_neutral", "level2.json"), model)
                        extra = ("--joint-templates", ",".join(tl))
                        print(f"    joint templates: {', '.join(tl)}", flush=True)
                    t0 = time.time()
                    run_l2(a, a.components, model, mask_name, arm, d, cfg, extra)
                    print(f"    l2 {arm}: {time.time() - t0:.0f}s", flush=True)

            if "controls" in steps:
                C = cfg["controls"]
                ids = {f"cholis_{cfg['selection']['roman'][i]}" for i in C["models"]}
                if model in ids and mask_name == C["mask"]:
                    for inj in C["inject"]:
                        d = os.path.join(od, f"control_{inj}")
                        if os.path.exists(os.path.join(d, "level2.json")):
                            continue
                        tl = joint_list(cfg, os.path.join(od, "l2_neutral", "level2.json"), model)
                        extra = ["--simulate-sky", inj, "--joint-templates", ",".join(tl)]
                        t0 = time.time()
                        run_l2(a, a.components, model, mask_name, C["arm"], d, cfg, extra)
                        print(f"    control {inj}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
