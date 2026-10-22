"""
one table across all models + the verdict, using decision_rule in
configs/gap5.yaml. missing files are blank, not zero. models that fail
admissibility are "cannot measure". counts are also given per family.
writes gap5_summary.json / .csv and gap5_table.md in outputs/gap5.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import run_gap5 as RG5            # noqa: E402

GC_RADIUS = 0.5


def _load(path):
    return json.load(open(path)) if os.path.exists(path) else None


def class_of(name, classes):
    import fnmatch
    for c, pats in classes.items():
        if any(fnmatch.fnmatch(name, p) for p in pats):
            return c
    return None


_TABLE = None


def family_of(model, cfg):
    # family = ISRF / gas / source profiles / B-field, from table viii
    global _TABLE
    if model in cfg["references"]:
        return "reference"
    if _TABLE is None:
        p = os.path.join(RG5.OUT, "selection", "table_viii.csv")
        _TABLE = {r["roman"]: r for r in csv.DictReader(open(p))} if os.path.exists(p) else {}
    t = _TABLE.get(model.replace("cholis_", ""))
    if not t:
        return "?"
    return (f"ISRF {t['ISRF'].split(',')[0]} | gas {t['xco']},{t['hi']},{t['hii']} | "
            f"{t['src_nuc']}/{t['src_el']} | B {t['B_field']}")


def gc_sources(srcs):
    return sum(1 for s in srcs or [] if np.hypot(s["l"], s["b"]) < GC_RADIUS)


def joint_verdict(j, classes, decisive):
    rows = [r for r in j["rows"] if class_of(r["template"], classes)]
    if not rows:
        return None
    best = min(rows, key=lambda r: r["AIC"])
    win = class_of(best["template"], classes)
    others = [r for r in rows if class_of(r["template"], classes) != win]
    gap = (min(r["AIC"] for r in others) - best["AIC"]) if others else None
    bulge = [r for r in rows if class_of(r["template"], classes) == "bulge"]
    sph = [r for r in rows if class_of(r["template"], classes) == "spherical"]
    b_best = min(bulge, key=lambda r: r["AIC"]) if bulge else None
    s_best = min(sph, key=lambda r: r["AIC"]) if sph else None
    return {
        "winner": best["template"], "winner_class": win,
        "winner_n_sources": best["n_sources"],
        "winner_gc_sources": gc_sources(best.get("sources")),
        "margin_to_next_class": gap,
        "decisive": gap is not None and gap >= decisive,
        "AIC_best": best["AIC"],
        "dAIC_winner_vs_best_bulge": (best["AIC"] - b_best["AIC"]) if b_best else None,
        "dAIC_spherical_minus_bulge": ((s_best["AIC"] - b_best["AIC"])
                                       if (s_best and b_best) else None),
        "spherical_n_sources": s_best["n_sources"] if s_best else None,
        "bulge_n_sources": b_best["n_sources"] if b_best else None,
        "all_converged": all(r.get("converged", True) for r in j["rows"]),
    }


def load_rows(outdir, cfg):
    classes = cfg["decision_rule"]["classes"]
    decisive = cfg["decision_rule"]["decisive_dAIC"]
    rows = []
    for model in sorted(os.listdir(outdir)):
        mdir = os.path.join(outdir, model)
        if model.startswith("_") or model.startswith("components") \
                or model in ("logs", "validate", "selection") or not os.path.isdir(mdir):
            continue
        for mask in sorted(os.listdir(mdir)):
            if not os.path.isdir(os.path.join(mdir, mask)):
                continue
            for tag in sorted(os.listdir(os.path.join(mdir, mask))):
                d = os.path.join(mdir, mask, tag)
                if not os.path.isdir(d):
                    continue
                l1 = _load(os.path.join(d, "level1.json"))
                ad = _load(os.path.join(d, "admit.json"))
                q = _load(os.path.join(d, "q.json"))
                jn = _load(os.path.join(d, "l2_joint", "level2.json"))
                nt = _load(os.path.join(d, "l2_neutral", "level2.json"))
                r = {"model": model, "family": family_of(model, cfg), "mask": mask,
                     "energy": tag}
                if l1:
                    r["lnL_background_L1"] = l1["lnL_background"]
                    r["AIC_background_L1"] = l1["AIC_background"]
                    r["n_params_background_L1"] = l1["n_params_background"]
                    for x in l1["rows"]:
                        r[f"L1_TS_{x['template']}"] = x["TS"]
                if ad:
                    r["admissible"] = ad["admissible"]
                    for c, v in ad["controls"].items():
                        r[f"admit_ts_injected_{c}"] = v["ts_injected"]
                        r[f"admit_gain_{c}"] = v["boosted_gain"]
                if q:
                    r["q_hat"] = q["projected_q_at_best"]
                    r["q_interval"] = [q["q_lo"], q["q_hi"]]
                    r["q_sigma_from_dimauro"] = q.get("sigma_from_dimauro")
                    r["q_at_ladder_edge"] = q.get("at_ladder_edge")
                if jn:
                    r["n_diffuse_components"] = jn.get("n_diffuse_components")
                    r.update(joint_verdict(jn, classes, decisive) or {})
                if nt:
                    by = {x["template"]: x for x in nt["rows"]}
                    r["neutral_n_sources"] = nt["rows"][0]["n_sources"] if nt["rows"] else None
                    r["neutral_gc_sources"] = gc_sources(nt.get("sources"))
                    for t in cfg["analysis"]["joint_templates"]:
                        if t in by:
                            r[f"neutral_TS_{t}"] = by[t]["TS"]
                    for cname in ("flattened", "spherical"):
                        m = [x for t, x in by.items() if class_of(t, {cname: classes[cname]})
                             and not t.endswith("_th1")]
                        if m and "vvv" in by:
                            b = max(m, key=lambda x: x["lnL"])
                            r[f"neutral_best_{cname}"] = b["template"]
                            r[f"neutral_dAIC_best_{cname}_minus_vvv"] = -2 * (b["lnL"] - by["vvv"]["lnL"])
                    if "winner" in r and r["winner"] in by:
                        r["winner_TS_neutral"] = by[r["winner"]]["TS"]
                    sph = [x for t, x in by.items() if class_of(t, {"s": classes["spherical"]})
                           and not t.endswith("_th1")]
                    if sph and "vvv" in by:
                        s = max(sph, key=lambda x: x["lnL"])
                        r["neutral_dAIC_spherical_minus_vvv"] = -2 * (s["lnL"] - by["vvv"]["lnL"])
                rows.append(r)
    return rows


def controls(outdir, cfg):
    classes = cfg["decision_rule"]["classes"]
    C = cfg["controls"]
    out = []
    for i in C["models"]:
        model = f"cholis_{cfg['selection']['roman'][i]}"
        for inj in C["inject"]:
            d = os.path.join(outdir, model, C["mask"],
                             f"{cfg['decision_rule']['primary']['energy'][0]:g}-"
                             f"{cfg['decision_rule']['primary']['energy'][1]:g}",
                             f"control_{inj}", "level2.json")
            j = _load(d)
            if not j:
                out.append({"model": model, "injected": inj, "status": "not run"})
                continue
            v = joint_verdict(j, classes, cfg["decision_rule"]["decisive_dAIC"])
            if inj == "null":
                gcs = sum(gc_sources(r.get("sources")) for r in j["rows"])
                ok = gcs == 0
                out.append({"model": model, "injected": "null",
                            "gc_sources_invented": gcs,
                            "n_sources": [r["n_sources"] for r in j["rows"]],
                            "pass": ok})
            else:
                want = class_of(inj, classes)
                out.append({"model": model, "injected": inj,
                            "recovered": v["winner"], "recovered_class": v["winner_class"],
                            "pass": v["winner_class"] == want})
    return out


def apply_decision_rule(rows, cfg):
    D = cfg["decision_rule"]
    P = D["primary"]
    tag = f"{P['energy'][0]:g}-{P['energy'][1]:g}"
    prim = [r for r in rows if r["mask"] == P["mask"] and r["energy"] == tag]
    fam2 = [r for r in prim if r["family"] not in ("reference",) and "AIC_best" in r]
    best_f2 = min((r["AIC_best"] for r in fam2), default=None)
    for r in prim:
        if "AIC_best" not in r:
            r["status"] = "incomplete"
            continue
        r["penalty_vs_best_family2"] = (r["AIC_best"] - best_f2) if best_f2 is not None else None
        flips = r["winner_class"] != D["gap1_verdict"] and r["decisive"]
        r["flip"] = bool(flips)
        swing = abs(r["dAIC_winner_vs_best_bulge"] or 0.0)
        r["competitive"] = (r["penalty_vs_best_family2"] is not None
                            and r["penalty_vs_best_family2"] < swing) if flips else None
        ts = r.get("winner_TS_neutral")
        r["excess_detected"] = None if ts is None else ts >= D["detected_ts_min"]
        if r.get("admissible") is False:
            r["status"] = "cannot measure"
        elif r.get("admissible") is None:
            r["status"] = "admissibility not run"
        else:
            r["status"] = "measured"

    allbest = min((r["AIC_best"] for r in prim if "AIC_best" in r), default=None)
    for r in prim:
        if "AIC_best" in r and allbest is not None:
            r["fit_gap_vs_best_model"] = r["AIC_best"] - allbest
    sets = {"reselected": {f"cholis_{m}" for m in cfg["selection"]["reselected"]},
            "original": {f"cholis_{m}" for m in cfg["selection"]["original"]}}
    sets["all"] = sets["reselected"] | sets["original"]
    by_set = {}
    for name, members in sets.items():
        rs = [r for r in prim if r["model"] in members]
        meas = [r for r in rs if r.get("status") == "measured"]
        rej = [r["model"] for r in meas
               if (r.get("flip") and r.get("competitive")) or r.get("excess_detected") is False]
        by_set[name] = {"in_table": len(rs), "measured": len(meas),
                        "winner_classes": {c: sum(r.get("winner_class") == c for r in rs)
                                           for c in ("bulge", "spherical", "flattened")},
                        "rejecting": rej,
                        "H0": ("not tested yet" if not meas else
                               ("rejected" if rej else "holds on the models measured"))}
    counted = [r for r in prim if r["family"] != "reference" and r.get("status") == "measured"]
    rejecting = [r["model"] for r in counted
                 if (r["flip"] and r["competitive"]) or r.get("excess_detected") is False]
    per_family = {}
    for f in sorted({r["family"] for r in prim}):
        m = [r for r in prim if r["family"] == f and "winner_class" in r]
        per_family[f] = {"models": len(m),
                         "bulge": sum(r["winner_class"] == "bulge" for r in m),
                         "spherical": sum(r["winner_class"] == "spherical" for r in m),
                         "flattened": sum(r["winner_class"] == "flattened" for r in m),
                         "admissible": sum(bool(r.get("admissible")) for r in m)}
    n_f2 = sum(1 for r in prim if r["family"] not in ("reference",))
    verdict = {
        "primary": P, "family2_models_in_table": n_f2,
        "family2_models_measured": len(counted),
        "models_that_reject_H0": rejecting,
        "H0": ("not tested yet" if not counted else
               ("rejected" if rejecting else "holds on the models measured")),
        "per_family": per_family,
        "by_set": by_set,
    }
    return verdict


def fmt(x, nd=1):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        return f"{x:,.{nd}f}"
    return str(x)


def print_table(rows, verdict, ctrl, path=None):
    P = verdict["primary"]
    tag = f"{P['energy'][0]:g}-{P['energy'][1]:g}"
    lines = [f"# Gap 5 summary — mask {P['mask']}, {tag} GeV, joint arm", ""]
    lines += ["| model | family | k | admissible | winner | margin | ΔAIC sph−bulge | "
              "src (bulge/sph) | GC src | winner TS | fit gap | penalty | q | status |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted([r for r in rows if r["mask"] == P["mask"] and r["energy"] == tag],
                    key=lambda r: (r["family"] != "reference", r["family"], r["model"])):
        lines.append(
            f"| {r['model']} | {r['family']} | {fmt(r.get('n_diffuse_components'))} | "
            f"{fmt(r.get('admissible'))} | {fmt(r.get('winner'))} | "
            f"{fmt(r.get('margin_to_next_class'))} | {fmt(r.get('dAIC_spherical_minus_bulge'))} | "
            f"{fmt(r.get('bulge_n_sources'))}/{fmt(r.get('spherical_n_sources'))} | "
            f"{fmt(r.get('winner_gc_sources'))} | {fmt(r.get('winner_TS_neutral'))} | "
            f"{fmt(r.get('fit_gap_vs_best_model'))} | "
            f"{fmt(r.get('penalty_vs_best_family2'))} | "
            f"{('≤' if r.get('q_at_ladder_edge') and (r.get('q_hat') or 1) < 0.5 else '')}"
            f"{fmt(r.get('q_hat'), 3)} | "
            f"{fmt(r.get('status'))} |")
    lines += ["", "fit gap: AIC of this model's best fit minus the best fit under any model "
              "(same data, mask and bins). penalty: the same against Family 2 only.",
              "", f"**H₀: {verdict['H0']}.** "
              f"{verdict['family2_models_measured']} of {verdict['family2_models_in_table']} "
              f"Family 2 models measured; rejecting: "
              f"{', '.join(verdict['models_that_reject_H0']) or 'none'}.", "",
              "| family | models | bulge | spherical | flattened | admissible |",
              "|---|---|---|---|---|---|"]
    for f, v in verdict["per_family"].items():
        lines.append(f"| {f} | {v['models']} | {v['bulge']} | {v['spherical']} | "
                     f"{v['flattened']} | {v['admissible']} |")
    if ctrl:
        lines += ["", "| control model | injected | recovered | pass |", "|---|---|---|---|"]
        for c in ctrl:
            rec = c.get("recovered", c.get("gc_sources_invented", c.get("status")))
            lines.append(f"| {c['model']} | {c['injected']} | {rec} | {fmt(c.get('pass'))} |")
    bs = verdict.get("by_set", {})
    if bs:
        lines += ["", "| selection set | measured | bulge | spherical | flattened | H₀ | rejecting |",
                  "|---|---|---|---|---|---|---|"]
        for name, v in bs.items():
            w = v["winner_classes"]
            lines.append(f"| {name} | {v['measured']}/{v['in_table']} | {w['bulge']} | "
                         f"{w['spherical']} | {w['flattened']} | {v['H0']} | "
                         f"{', '.join(v['rejecting']) or 'none'} |")
    txt = "\n".join(lines)
    print(txt)
    if path:
        open(path, "w").write(txt + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=RG5.OUT)
    ap.add_argument("--config", default=RG5.CONFIG)
    a = ap.parse_args()
    cfg = RG5.load_config(a.config)
    rows = load_rows(a.outdir, cfg)
    verdict = apply_decision_rule(rows, cfg)
    ctrl = controls(a.outdir, cfg)
    json.dump({"verdict": verdict, "rows": rows, "controls": ctrl},
              open(os.path.join(a.outdir, "gap5_summary.json"), "w"), indent=1,
              default=float)
    keys = sorted({k for r in rows for k in r if not isinstance(r[k], (list, dict))})
    with open(os.path.join(a.outdir, "gap5_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print_table(rows, verdict, ctrl, os.path.join(a.outdir, "gap5_table.md"))


if __name__ == "__main__":
    main()
