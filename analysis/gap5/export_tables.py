"""Flat CSVs from the Gap 5 runs.

    outputs/gap5/gap5_models.csv     one row per background: winners, margins, flips
    outputs/gap5/gap5_fits.csv       every Level 2 fit (background x arm x template)
    outputs/gap5/gap5_controls.csv   the injection controls
    outputs/gap5/level2_sources.csv  every source Level 2 added, every run
    outputs/gap7/covariance_input.csv  Gap 1 rows + Gap 5 rows, Gap 1 column layout

Mask B, 1-10 GeV unless the run says otherwise. AIC = 2k - 2 lnL throughout.
"""
import argparse
import csv
import datetime
import json
import subprocess
from fnmatch import fnmatch
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CLASSES = {"bulge": ["vvv", "f98"], "spherical": ["gnfw_*"], "flattened": ["hestia_*"]}
GAP7_CLASS = {"bulge": "stellar", "spherical": "dark_matter_spherical",
              "flattened": "dark_matter_simulated"}
REFS = {"gll_iem_v07": 1, "pohl22_rings_species": 2}


def cls(t):
    for c, pats in CLASSES.items():
        if any(fnmatch(t, p) for p in pats):
            return c
    return "other"


def aic(r):
    return r["AIC"] if "AIC" in r else 2 * r["n_params"] - 2 * r["lnL"]


def best_by_class(rows):
    out = {}
    for r in rows:
        c, a = cls(r["template"]), aic(r)
        if c not in out or a < out[c][1]:
            out[c] = (r["template"], a)
    return out


def winner(rows):
    b = best_by_class(rows)
    order = sorted(b.items(), key=lambda kv: kv[1][1])
    (wc, (wt, wa)), margin = order[0], (order[1][1][1] - order[0][1][1] if len(order) > 1 else None)
    return wc, wt, wa, margin, b


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def roman(model):
    return model.replace("cholis_", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap5", default=str(ROOT / "outputs/gap5"))
    ap.add_argument("--gap1-csv", default=str(ROOT / "outputs/gap1/gap1_results_for_gap7.csv"))
    ap.add_argument("--out7", default=str(ROOT / "outputs/gap7"))
    a = ap.parse_args()
    g5 = Path(a.gap5)
    cfg = yaml.safe_load((ROOT / "configs/gap5.yaml").read_text())
    resel = set(cfg["selection"]["reselected"])
    orig = set(cfg["selection"]["original"])
    summary = json.loads((g5 / "gap5_summary.json").read_text())
    srows = {r["model"]: r for r in summary["rows"] if r["mask"] == "B" and r["energy"] == "1-10"}

    try:
        commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                                         text=True).strip()
    except Exception:
        commit = ""
    today = datetime.date.today().isoformat()

    models, fits, sources, gap7 = [], [], [], []
    ranked = sorted((m for m in srows if "AIC_best" in srows[m]),
                    key=lambda m: srows[m]["fit_gap_vs_best_model"])

    for run_dir in sorted(g5.glob("*/*/*/l2_*")):
        model, mask, erange, arm = run_dir.parts[-4], run_dir.parts[-3], run_dir.parts[-2], run_dir.name
        d = load(run_dir / "level2.json")
        if not d or model == "_reference":
            continue
        arm = arm.replace("l2_", "")
        best = min(aic(r) for r in d["rows"])
        vvv = next((aic(r) for r in d["rows"] if r["template"] == "vvv"), None)
        gn = [aic(r) for r in d["rows"] if cls(r["template"]) == "spherical"]
        for r in d["rows"]:
            fits.append({"background": model, "mask": mask, "energy_gev": erange, "arm": arm,
                         "template": r["template"], "class": cls(r["template"]),
                         "lnL": round(r["lnL"], 3), "n_params": r["n_params"],
                         "AIC": round(aic(r), 3), "dAIC_vs_best_in_run": round(aic(r) - best, 3),
                         "n_sources": r["n_sources"], "converged": r.get("converged")})
            if arm in ("neutral", "joint") and model not in ("_reference",):
                ref = REFS.get(model)
                emin, emax = d["emin"], d["emax"]
                gap7.append({
                    "result_id": f"gap5_mask{mask}ps_{model}_L2{arm}_{r['template']}_{erange}GeV",
                    "gap": 5, "mask": mask, "point_source_template": d.get("ps", True),
                    "iem_family": ref if ref else 3, "diffuse_model": model,
                    "n_diffuse_components": d.get("n_diffuse_components"), "level": 2,
                    "simulated_sky": "", "n_sources_added": r["n_sources"],
                    "energy_bin": f"{erange}GeV", "emin_gev": emin, "emax_gev": emax,
                    "template": r["template"], "template_class": GAP7_CLASS.get(cls(r["template"]), ""),
                    "halo": r["template"].split("_")[1] if r["template"].startswith("hestia_") else "",
                    "view_angle": r["template"].split("ang")[1] if "_ang" in r["template"] else "",
                    "TS": round(r["TS"], 3) if "TS" in r else "", "lnL": round(r["lnL"], 3),
                    "dAIC_vs_vvv": round(aic(r) - vvv, 3) if vvv is not None else "",
                    "dAIC_vs_best_gnfw": round(aic(r) - min(gn), 3) if gn else "",
                    "signal_photons": "", "converged": r.get("converged"),
                    "max_lnL_gap_bound": r.get("max_lnL_gap", ""),
                    "mask_kept_fraction": d.get("kept_fraction"),
                    "admissible": srows.get(model, {}).get("admissible", "") if (mask, erange) == ("B", "1-10") else "",
                    "control_injected_ts": "", "control_template": "",
                    "status": "OK" if arm == "joint" else "OK_common_source_list",
                    "source_run": f"gap5_l2_{arm}", "git_commit": commit, "export_date": today,
                    "arm": arm, "iem_set": ("reference" if ref else
                                            "+".join(s for s, ok in (("reselected", roman(model) in resel),
                                                                     ("original", roman(model) in orig)) if ok)),
                    "fit_gap_vs_best_iem": (round(srows[model]["fit_gap_vs_best_model"], 1)
                                            if (mask, erange) == ("B", "1-10") and model in srows
                                            and "fit_gap_vs_best_model" in srows[model] else ""),
                })
        if arm == "neutral":
            for s in d.get("sources", []):
                sources.append({"background": model, "mask": mask, "energy_gev": erange, "arm": arm,
                                "template": "(common list)", "source": s["name"], "l_deg": round(s["l"], 3),
                                "b_deg": round(s["b"], 3), "residual_TS": round(s["residual_ts"], 2),
                                "pass_threshold": s["pass"], "order_added": s["index"]})
        else:
            for r in d["rows"]:
                for s in r.get("sources", []):
                    sources.append({"background": model, "mask": mask, "energy_gev": erange, "arm": arm,
                                    "template": r["template"], "source": s["name"], "l_deg": round(s["l"], 3),
                                    "b_deg": round(s["b"], 3), "residual_TS": round(s["residual_ts"], 2),
                                    "pass_threshold": s["pass"], "order_added": s["index"]})

    for rank, m in enumerate(ranked, 1):
        base = g5 / m / "B" / "1-10"
        neu, joi = load(base / "l2_neutral/level2.json"), load(base / "l2_joint/level2.json")
        row = {"background": roman(m), "family": srows[m]["family"],
               "set": ("reference" if m in REFS else
                       "+".join(s for s, ok in (("reselected", roman(m) in resel),
                                                ("original", roman(m) in orig)) if ok)),
               "fit_rank": rank, "fit_gap_vs_best": round(srows[m]["fit_gap_vs_best_model"], 1),
               "AIC_best": round(srows[m]["AIC_best"], 1), "admissible": srows[m].get("admissible")}
        for name, d in (("neutral", neu), ("joint", joi)):
            if not d:
                continue
            wc, wt, wa, mg, b = winner(d["rows"])
            row[f"{name}_winner_class"], row[f"{name}_winner"] = wc, wt
            row[f"{name}_margin_to_next_class"] = round(mg, 1) if mg is not None else ""
            for c in ("bulge", "spherical", "flattened"):
                if c in b:
                    row[f"{name}_best_{c}"] = b[c][0]
                    row[f"{name}_AIC_best_{c}"] = round(b[c][1], 1)
            nsrc = {r["template"]: r["n_sources"] for r in d["rows"]}
            row[f"{name}_n_sources_winner"] = nsrc[wt]
            if name == "joint" and "spherical" in b:
                row["joint_n_sources_best_spherical"] = nsrc[b["spherical"][0]]
        row["neutral_joint_disagree"] = row.get("neutral_winner_class") != row.get("joint_winner_class")
        row["winner_TS_neutral"] = round(srows[m].get("winner_TS_neutral") or 0, 1)
        row["q_hat"] = srows[m].get("q_hat")
        row["status"] = srows[m].get("status")
        models.append(row)

    controls = []
    for c in summary["controls"]:
        d = load(g5 / c["model"] / "B" / "1-10" / f"control_{c['injected']}" / "level2.json")
        out = {"background": roman(c["model"]), "injected": c["injected"],
               "recovered": c.get("recovered", "nothing"), "recovered_class": c.get("recovered_class", ""),
               "pass": c["pass"]}
        if d:
            wc, wt, wa, mg, b = winner(d["rows"])
            out["margin_to_next_class"] = round(mg, 1) if mg is not None else ""
            inj = cls(c["injected"]) if c["injected"] != "null" else None
            if inj and inj != wc and inj in b:
                out["true_class_behind_by"] = round(b[inj][1] - wa, 1)
            out["n_sources"] = max(r["n_sources"] for r in d["rows"])
        controls.append(out)

    def write(path, rows, fields=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        fields = fields or list(dict.fromkeys(k for r in rows for k in r))
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(path.relative_to(ROOT), len(rows))

    write(g5 / "gap5_models.csv", models)
    write(g5 / "gap5_fits.csv", fits)
    write(g5 / "gap5_controls.csv", controls)
    write(g5 / "level2_sources.csv", sources)

    with open(a.gap1_csv) as f:
        g1 = list(csv.DictReader(f))
    fields = list(g1[0].keys()) + ["arm", "iem_set", "fit_gap_vs_best_iem"]
    write(Path(a.out7) / "covariance_input.csv", g1 + gap7, fields)


if __name__ == "__main__":
    main()
