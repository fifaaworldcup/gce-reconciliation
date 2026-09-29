"""
Gap 1 results in one long table for Gap 7's covariance work.

One row per (mask, IEM family, point-source treatment, energy range,
template). result_id follows Conventions sec 11:
    gap{number}_{mask}_{iem_family}_{template}_{energy_bin}
with the point-source template recorded as a mask suffix ("maskBps") since
the Conventions have no slot for it.

dAIC columns: every model in a row group shares the same background and the
same parameter count, so dAIC = -2 dlnL and compares directly within a group.
It does NOT compare across groups (different masks keep different pixels).
"""
import csv
import datetime
import glob
import json
import os
import subprocess

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "outputs", "gap1", "gap1_results_for_gap7.csv")


def commit():
    try:
        return subprocess.check_output(["git", "-C", ROOT, "rev-parse", "--short",
                                        "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def ebin_label(lo, hi):
    f = lambda x: (f"{x:g}").replace(".", "p")
    return f"{f(lo)}-{f(hi)}GeV"


def main():
    rows = []
    c = commit()
    today = datetime.date.today().isoformat()
    # full-range cube with finely sampled exposure (eb_*, tri_*) supersedes
    # the original 1-10 GeV cube (scan_mask*) wherever both cover a cell, so
    # every row in a cell shares one cube and one exposure
    def prio(d):
        b = os.path.basename(d)
        return (0 if b.startswith("eb_") else 1 if b.startswith("tri_") else 2, b)
    cell_source = {}
    for d in sorted(glob.glob(os.path.join(ROOT, "outputs", "gap1", "*scan*")), key=prio):
        js = os.path.join(d, "scan.json")
        if not os.path.exists(js):
            continue
        S = json.load(open(js))
        name = os.path.basename(d)
        # energy range: eb_scan_* carry it in the log; older scans were 1-10
        lo, hi = 1.0, 10.0
        log = os.path.join(d, "run.log")
        if os.path.exists(log):
            for line in open(log):
                if line.startswith("energy "):
                    a, b = line.split()[1].split("-")
                    lo, hi = float(a), float(b.rstrip(",GeV"))
                    break
        # 2026-09-27: what these runs called mask D is now mask E. Mask D is
        # Cholis 2022 Table III (see masks/make_masks.py). Relabel on the way
        # out so nobody downstream compares an E row to a real D row.
        mask_letter = "E" if S["mask"] == "D" else S["mask"]
        S = {**S, "mask": mask_letter}
        mask = mask_letter + ("ps" if S.get("point_source_template") else "")
        cell = (mask, ebin_label(lo, hi))
        cube = "full" if not name.startswith("scan_") else "orig"
        if cell_source.setdefault(cell, cube) != cube:
            continue
        g = {r["template"]: r for r in S["rows"]}
        gn = [r for k, r in g.items() if k.startswith("gnfw_") and not k.endswith("_th1")]
        best_gnfw = max(gn, key=lambda r: r["lnL"])["lnL"] if gn else None
        for k, r in g.items():
            rows.append({
                "result_id": f"gap1_mask{mask}_iem1_{k}_{ebin_label(lo, hi)}",
                "gap": 1, "mask": S["mask"],
                "point_source_template": bool(S.get("point_source_template")),
                "iem_family": 1, "level": 1,
                "energy_bin": ebin_label(lo, hi), "emin_gev": lo, "emax_gev": hi,
                "template": k,
                "TS": round(r["TS"], 3), "lnL": round(r["lnL"], 3),
                "dAIC_vs_vvv": round(r["dAIC_vs_vvv"], 3),
                "dAIC_vs_best_gnfw": (round(-2 * (r["lnL"] - best_gnfw), 3)
                                      if best_gnfw is not None else ""),
                "signal_photons": round(r["photons"], 1),
                "converged": r["converged"],
                "mask_kept_fraction": round(S["kept_fraction"], 4),
                # sec 6.0: these are all single-template-IEM fits, and 98% of
                # the TS in the only sensitive configuration is absorbed by
                # giving that template latitude freedom. Do not draw a
                # conclusion from these rows until a ring-decomposed IEM
                # exists. The column is here so the CSV cannot be used
                # without seeing it.
                "status": "WITHDRAWN_see_section_6.0_rigid_IEM",
                "source_run": name, "git_commit": c, "export_date": today,
            })
    # a triaxial library run can repeat templates from the full scan in the same
    # cell; keep the first so every result_id is unique
    seen, uniq = set(), []
    for r in rows:
        if r["result_id"] in seen:
            continue
        seen.add(r["result_id"]); uniq.append(r)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(uniq[0]))
        w.writeheader(); w.writerows(uniq)
    cells = sorted({(r["mask"], r["point_source_template"], r["energy_bin"]) for r in uniq})
    print(f"wrote {OUT}: {len(uniq)} rows, {len(cells)} mask x PS x energy cells")
    for m, ps, e in cells:
        print(f"   mask {m}{' +PS' if ps else '    '}  {e}")


if __name__ == "__main__":
    main()
