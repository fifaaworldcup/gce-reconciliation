"""Gap 1 results in one long table for Gap 7's covariance work.

One row per (background model, mask, point-source treatment, analysis level,
energy range, template). Rewritten 2026-10-11: the previous version globbed
every directory with "scan" in the name, hardcoded `iem_family: 1`,
`level: 1`, and `status: WITHDRAWN`, and had no column for which diffuse
model a row came from. After the Pohl22 re-run none of that is true any more
and the withdrawn flag was on every row, including the ones that replaced the
withdrawn ones.

Four things this table now carries that the covariance needs and the old one
could not express.

`diffuse_model` and `n_diffuse_components`. The background is the dominant
systematic in this gap, and "which interstellar model" and "how much freedom
inside it" are two different axes. gll_iem_v07 with one normalisation and
Pohl22 with sixteen are not two points on one scale.

`admissible`. Whether that configuration passed its own injection control -
inject this template at a stated flux and check it comes back above the
pre-registered threshold. Rows marked false are not weak evidence, they are
no evidence: the configuration could not have seen the signal if it were
there. **They must be excluded before any covariance, not down-weighted.**

`halo` and `view_angle`. The 24 HESTIA rows are six galaxies seen from four
angles each, not 24 independent measurements, and every row in a cell is
fitted to the same photons besides. Treating rows as independent samples is
what took a p-value from 0.07 to 5e-5 in the first version of sec 6.1. Group
on `halo` before counting degrees of freedom.

`status`. Per run, not hardcoded: OK for the ring-decomposed results,
WITHDRAWN_* for the single-template ones kept for comparison, and
CONTROL_NOT_A_RESULT for runs on simulated skies, which are in the file
because leaving them out silently is worse, and must never enter a fit.
"""
import csv
import datetime
import glob
import json
import os
import re
import subprocess

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "gap1", "gap1_results_for_gap7.csv")

# Every run directory we export, with what it is. Explicit rather than
# globbed: a glob cannot tell a result from a control.
RUNS = {
    # Pohl22 ring decomposition - the current results
    "rsp_scan_maskB_ps":  dict(iem=2, status="OK"),
    "rsp_scan_maskD_ps":  dict(iem=2, status="OK"),
    "rsp_tri_maskB_ps":   dict(iem=2, status="OK"),
    "rsp_eb_1-3":         dict(iem=2, status="OK"),
    "rsp_eb_3-10":        dict(iem=2, status="OK"),
    "rsp_eb_0_5-1":       dict(iem=2, status="OK_EMPTY_no_detection_any_template"),
    "rsp_eb_10-50":       dict(iem=2, status="OK"),
    "rings_scan_maskB_ps": dict(iem=2, status="OK"),
    "rings_scan_maskD_ps": dict(iem=2, status="NOT_ADMISSIBLE_see_table_R3"),
    "rings_scan_maskA":    dict(iem=2, status="NOT_ADMISSIBLE_mask_keeps_8.7pct"),
    "rsp_scan_maskA":      dict(iem=2, status="NOT_ADMISSIBLE_mask_keeps_8.7pct"),
    "rsp_scan_maskA_ps":   dict(iem=2, status="NOT_ADMISSIBLE_mask_keeps_8.7pct"),
    # --emin 1 --emax 3 on the 1-10 GeV cube keeps whole bins only: 1.00-2.37
    "rsp_scan_1-3":        dict(iem=2, status="OK", energy=(1.0, 2.371)),
    "rings_eb_0_5-1":     dict(iem=2, status="OK_EMPTY_no_detection_any_template"),
    "rings_eb_1-3":       dict(iem=2, status="OK"),
    "rings_eb_3-10":      dict(iem=2, status="OK"),
    "rings_eb_10-50":     dict(iem=2, status="OK"),
    "rings_tri_maskB_ps": dict(iem=2, status="OK"),
    # Level 2 iterative source finding, one common source list. Its
    # admissibility is the Level 1 control on the same background, mask and
    # energy range (the Gap 5 rule)
    "l2_neutral_maskB_ps": dict(iem=2, status="OK", level=2,
                                admit_from="rsp_admit_maskB_ps"),
    # single all-sky template - superseded, kept so the change is visible
    "scan_maskB_ps":      dict(iem=1, status="WITHDRAWN_see_section_6.0_rigid_IEM"),
    "scan_maskA":         dict(iem=1, status="WITHDRAWN_see_section_6.0_rigid_IEM"),
    "eb_scan_all":        dict(iem=1, status="WITHDRAWN_see_section_6.0_rigid_IEM"),
    "tri_scan_b1-10":     dict(iem=1, status="WITHDRAWN_see_section_6.0_rigid_IEM"),
}
CONTROL_PREFIXES = ("rsp_sim_", "rsp_simfit_", "l2_sim_")
# Every pre-Pohl22 run is named scan_*, eb_* or tri_*, and all of them are
# single-template-IEM fits superseded by sec 6.0b. A rule rather than a list,
# so a legacy run that turns up later is classified correctly instead of
# landing in the file unlabelled.
LEGACY_PREFIXES = ("scan_", "eb_", "tri_")
_HESTIA = re.compile(r"^hestia_(G\d\.\d)_ang(\d+)$")


def commit():
    try:
        return subprocess.check_output(
            ["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def ebin_label(lo, hi):
    f = lambda x: f"{x:g}".replace(".", "p")
    return f"{f(lo)}-{f(hi)}GeV"


def template_class(name):
    if name in ("vvv", "f98"):
        return "stellar"
    if name.startswith("gnfw_"):
        return "dark_matter_spherical"
    if name.startswith("hestia_"):
        return "dark_matter_simulated"
    if name.startswith("triaxial"):
        return "dark_matter_triaxial"
    return "other"


def energy_range(name, S):
    """robustness_scan.py takes --emin/--emax but never wrote them to
    scan.json, so read them from the run name. Every unnamed run is 1-10 GeV
    (gap1_ccube.fits); *_all is the whole 0.5-50 GeV cube."""
    if S.get("emin") is not None and S.get("emax") is not None:
        return float(S["emin"]), float(S["emax"])
    m = re.search(r"_b?(\d+(?:_\d+)?)-(\d+(?:_\d+)?)(?:_|$)", name)
    if m:
        return (float(m.group(1).replace("_", ".")),
                float(m.group(2).replace("_", ".")))
    if name.endswith("_all"):
        return 0.5, 50.0
    return 1.0, 10.0


def _read_admit(d, template):
    f = os.path.join(d, f"admissibility_{template}.json")
    if os.path.exists(f):
        a = json.load(open(f))
        return a["admissible"], a["ts_injected"], a["control"]
    return None


def admissibility(outdir, template, admit_from=None):
    """Did this configuration pass an injection control for this template?

    1. the template's own control, in the run or its *_admit_* directory
    2. a control of the same class (gnfw, hestia, bulge) in those directories;
       if there are several, all of them must pass
    3. the run's own control, again only if it is the same class

    A control of another class says nothing: a background that recovers an
    injected bulge has not been shown to recover a halo. Before 2026-10-26 the
    run's VVV control was applied to every template, which marked mask D's
    gNFW rows admissible when the gNFW control fails there.
    """
    dirs = [outdir, outdir.replace("_scan", "_admit").replace("scan_", "admit_")]
    if admit_from:
        dirs.insert(0, os.path.join(os.path.dirname(outdir), admit_from))
    for d in dirs:
        a = _read_admit(d, template)
        if a:
            return a
    cls = template_class(template)
    same = []
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "admissibility_*.json"))):
            t = os.path.basename(f)[len("admissibility_"):-len(".json")]
            if template_class(t) == cls:
                same.append(_read_admit(d, t))
        if same:
            break
    if same:
        bad = [x for x in same if not x[0]]
        pick = min(bad or same, key=lambda x: x[1])
        return (not bad), pick[1], pick[2]
    js = os.path.join(outdir, "scan.json")
    if os.path.exists(js):
        a = (json.load(open(js)) or {}).get("admissibility")
        if a and template_class(a["control"]) == cls:
            return a["admissible"], a["ts_injected"], a["control"]
    return "", "", ""


def main():
    rows, c, today = [], commit(), datetime.date.today().isoformat()
    found_dirs = set()
    for path in sorted(glob.glob(os.path.join(ROOT, "outputs", "gap1", "*"))):
        name = os.path.basename(path)
        js = os.path.join(path, "scan.json")
        if not os.path.exists(js):
            continue
        found_dirs.add(name)
        if name in RUNS:
            meta = RUNS[name]
        elif name.startswith(CONTROL_PREFIXES):
            meta = dict(iem=2, status="CONTROL_NOT_A_RESULT")
        elif name.startswith(LEGACY_PREFIXES):
            meta = dict(iem=1, status="WITHDRAWN_see_section_6.0_rigid_IEM")
        else:
            meta = dict(iem=0, status="UNCLASSIFIED_run_not_in_export_manifest")
        S = json.load(open(js))
        lo, hi = meta.get("energy") or energy_range(name, S)
        mask = S["mask"] + ("ps" if S.get("point_source_template") else "")
        # robustness_scan writes simulated_sky as a template name or null;
        # run_level2 writes a dict with the name and the injected flux
        raw = S.get("simulated_sky")
        sim = raw.get("template") if isinstance(raw, dict) else raw
        simtag = f"_sim{sim}" if sim else ("_simnull" if raw else "")
        g = {r["template"]: r for r in S["rows"]}
        gn = [r for k, r in g.items()
              if k.startswith("gnfw_") and not k.endswith("_th1")]
        best_gnfw = max(gn, key=lambda r: r["lnL"])["lnL"] if gn else None
        for k, r in g.items():
            m = _HESTIA.match(k)
            adm, inj, ctl = admissibility(path, k, meta.get("admit_from"))
            rows.append({
                "result_id": f"gap1_mask{mask}_iem{meta['iem']}_"
                             f"{S.get('diffuse_model', 'gll_iem_v07')}_"
                             f"L{meta.get('level', 1)}_{k}_"
                             f"{ebin_label(lo, hi)}{simtag}",
                "gap": 1, "mask": S["mask"],
                "point_source_template": bool(S.get("point_source_template")),
                "iem_family": meta["iem"],
                "diffuse_model": S.get("diffuse_model", "gll_iem_v07"),
                "n_diffuse_components": S.get("n_diffuse_components", 1),
                "level": meta.get("level", 1),
                "simulated_sky": sim or ("null" if raw else ""),
                "n_sources_added": r.get("n_sources", 0),
                "energy_bin": ebin_label(lo, hi), "emin_gev": lo, "emax_gev": hi,
                "template": k, "template_class": template_class(k),
                "halo": m.group(1) if m else "",
                "view_angle": int(m.group(2)) if m else "",
                "TS": round(r["TS"], 3), "lnL": round(r["lnL"], 3),
                "dAIC_vs_vvv": round(r["dAIC_vs_vvv"], 3),
                "dAIC_vs_best_gnfw": (round(-2 * (r["lnL"] - best_gnfw), 3)
                                      if best_gnfw is not None else ""),
                "signal_photons": round(r["photons"], 1) if "photons" in r else "",
                "converged": r["converged"],
                "max_lnL_gap_bound": r.get("max_lnL_gap", ""),
                "mask_kept_fraction": round(S["kept_fraction"], 4),
                "admissible": adm, "control_injected_ts": inj,
                "control_template": ctl,
                "status": meta["status"],
                "source_run": name, "git_commit": c, "export_date": today,
            })
    missing = [k for k in RUNS if k not in found_dirs]
    # two runs can share a configuration (eb_scan_b1-10 and scan_maskB_ps
    # are both gll, mask B + PS, 1-10 GeV, from different cubes). Keep both
    # and say which is which; dropping one silently is how the mask B Pohl22
    # scan went missing from this file before.
    count = {}
    for r in rows:
        count[r["result_id"]] = count.get(r["result_id"], 0) + 1
    uniq = []
    for r in rows:
        if count[r["result_id"]] > 1:
            r["result_id"] += f"_{r['source_run']}"
        uniq.append(r)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(uniq[0]))
        w.writeheader(); w.writerows(uniq)
    by = {}
    for r in uniq:
        by[r["status"]] = by.get(r["status"], 0) + 1
    print(f"wrote {OUT}: {len(uniq)} rows from {len(found_dirs)} runs")
    for k, v in sorted(by.items(), key=lambda kv: -kv[1]):
        print(f"  {v:5d}  {k}")
    if missing:
        print(f"  manifest lists {len(missing)} runs with no scan.json yet: "
              f"{', '.join(sorted(missing))}")
    usable = [r for r in uniq if r["status"] == "OK" and r["admissible"] is True]
    halos = {r["halo"] for r in usable if r["halo"]}
    print(f"\n  rows safe for a covariance (status OK and admissible true): "
          f"{len(usable)}")
    print(f"  distinct HESTIA halos among them: {len(halos)} - group on `halo`, "
          f"do not count rows")


if __name__ == "__main__":
    main()
