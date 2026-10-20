"""
table viii rebuilt from the 80 galdef files in the zenodo release, and the
re-selection from configs/gap5.yaml: the six table II models as anchors, then
farthest-point (largest distance to the nearest chosen model) for the rest.

distance: continuous columns min-max scaled to [0, 1], B-field unpacked into
B0, r_c, z_c, categorical columns 0 if equal, 1 if not. every value of each
categorical column has to be in the set before the free picks start.

  python analysis/gap5/select_models.py --zenodo <ZENODO_FILES folder>
"""

from __future__ import annotations

import argparse
import csv
import json
import os

import numpy as np

HI_CODE = {"rbands_hi12_v2_qdeg_zmax1_Ts125.fits": 1,
           "rbands_hi12_v2_qdeg_zmax1_Ts100000_EBV_mag2_limit.fits.gz": 2,
           "rbands_hi12_v2_qdeg_zmax1_Ts100000_EBV_mag5_limit.fits.gz": 3,
           "rbands_hi12_v2_qdeg_zmax1_Ts150_EBV_mag2_limit.fits.gz": 4,
           "rbands_hi12_v2_qdeg_zmax1_Ts150_EBV_mag5_limit.fits.gz": 5}
SRC = {"3": "Pul", "2": "SNR"}
CONT = ["zL", "D0", "delta", "vA", "dvc", "ap1", "ap2", "ae1", "ae2", "Ne",
        "isrf", "B0", "rc", "zc"]
CAT = ["xco", "hi", "hii", "src_nuc", "src_el"]


def read_galdef(path):
    d = {}
    for line in open(path, errors="replace"):
        if line.lstrip().startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        v = v.split()
        d[k.strip()] = v[0] if v else ""
    return d


def table(zen):
    names = []
    for line in open(os.path.join(zen, "NAMING_CONVENTION_OF_DIFFUSE_EMISSION_MODELS.dat")):
        if line.strip() and not line.startswith("#"):
            names.append(line.split())
    rows = []
    for roman, code, gfile in names:
        d = read_galdef(os.path.join(zen, "GALDEF_FILES", f"galdef_54_{gfile}"))
        f = lambda k: float(d[k].replace("+", ""))
        b = d["B_field_model"]
        rows.append({
            "roman": roman, "code": code, "galdef": gfile,
            "zL": f("z_max"), "D0": f("D0_xx") / 1e28, "delta": f("D_g_2"),
            "vA": f("v_Alfven"), "dvc": f("dvdz_conv") if d["convection"] == "1" else 0.0,
            "ap1": f("nuc_g_1"), "ap2": f("nuc_g_2"),
            "ae1": f("electron_g_0"), "ae2": f("electron_g_1"),
            "Ne": f("electron_norm_flux") * 1e9,
            "isrf": float(d["ISRF_factors"].split(",")[0]), "ISRF": d["ISRF_factors"],
            "B_field": b, "B0": int(b[0:3]) / 10, "rc": int(b[3:6]) / 10, "zc": int(b[6:9]) / 10,
            "xco": d["n_X_CO"], "hi": HI_CODE[d["HIR_filename"]], "hii": d["nHII_model"],
            "src_nuc": SRC[d["source_model"]], "src_el": SRC[d["source_model_elec"]],
        })
    return rows


def features(rows):
    X = np.array([[r[c] for c in CONT] for r in rows], float)
    lo, hi = X.min(0), X.max(0)
    Xs = (X - lo) / np.where(hi > lo, hi - lo, 1.0)
    C = np.array([[str(r[c]) for c in CAT] for r in rows])
    return Xs, C


def dist(Xs, C, i, j):
    return float(np.sqrt(((Xs[i] - Xs[j]) ** 2).sum() + (C[i] != C[j]).sum()))


def select(rows, anchors, n_total):
    Xs, C = features(rows)
    idx = {r["roman"]: k for k, r in enumerate(rows)}
    chosen = [idx[a] for a in anchors]
    log = [(r, "anchor") for r in anchors]

    def mind(k):
        return min(dist(Xs, C, k, j) for j in chosen)

    need = {(c, v) for c in range(len(CAT)) for v in set(C[:, c])}
    need -= {(c, C[j, c]) for j in chosen for c in range(len(CAT))}
    while need and len(chosen) < n_total:
        cand = [k for k in range(len(rows)) if k not in chosen]
        gain = {k: sum((c, C[k, c]) in need for c in range(len(CAT))) for k in cand}
        top = max(gain.values())
        if top == 0:
            break
        k = max((k for k in cand if gain[k] == top), key=mind)
        covers = [f"{CAT[c]}={C[k, c]}" for c in range(len(CAT)) if (c, C[k, c]) in need]
        need -= {(c, C[k, c]) for c in range(len(CAT))}
        chosen.append(k)
        log.append((rows[k]["roman"], "covers " + ", ".join(covers)))
    while len(chosen) < n_total:
        cand = [k for k in range(len(rows)) if k not in chosen]
        k = max(cand, key=mind)
        log.append((rows[k]["roman"], f"farthest point, nearest chosen at {mind(k):.2f}"))
        chosen.append(k)
    return [rows[k]["roman"] for k in chosen], log, need


def spread(rows, picks):
    Xs, C = features(rows)
    idx = {r["roman"]: k for k, r in enumerate(rows)}
    ks = [idx[p] for p in picks]
    nn = [min(dist(Xs, C, a, b) for b in ks if b != a) for a in ks]
    return min(nn), float(np.median(nn))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zenodo", required=True)
    ap.add_argument("--anchors", default="I,II,III,VI,X,XV")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--compare", default="V,VII,XVII,XIX,XXII,XXV,XXVII,XXX,LI,LII,LIII,LVI")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    rows = table(a.zenodo)
    picks, log, left = select(rows, a.anchors.split(","), a.n)
    for r, why in log:
        print(f"  {r:7s} {why}")
    if left:
        print(f"  not covered within {a.n}: {sorted(left)}")
    old = a.compare.split(",")
    for name, s in (("re-selected", picks), ("original 12", old)):
        sub = [r for r in rows if r["roman"] in s]
        mn, md = spread(rows, s)
        print(f"\n{name}: nearest-neighbour distance min {mn:.2f}, median {md:.2f}")
        for c in ("vA", "dvc", "zL", "D0"):
            print(f"  {c:5s} {min(r[c] for r in sub):g} - {max(r[c] for r in sub):g}"
                  f"   (library {min(r[c] for r in rows):g} - {max(r[c] for r in rows):g})")
        for c in CAT:
            print(f"  {c:7s} {sorted({str(r[c]) for r in sub})}   (library {sorted({str(r[c]) for r in rows})})")
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "table_viii.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        json.dump({"reselected": picks, "original": old, "log": log},
                  open(os.path.join(a.out, "selection.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
