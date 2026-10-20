"""
writes gll_iem_v07 and gap 1's pohl22 rings_species in the same component
format as the cholis models, so everything runs through one path.
--check confirms the component version gives the same lnL as fit_rings.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "analysis", "gap1"), os.path.join(ROOT, "utils"),
          os.path.join(ROOT, "masks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_cubes as FC          # noqa: E402
import fit_rings as FR          # noqa: E402
import convert_iem as CI        # noqa: E402


def export(ccube, iem_path, rings, out):
    _, e_lo, e_hi, _ = FC.read_ccube(ccube)
    emid = np.sqrt(e_lo * e_hi)
    nb = len(emid)

    gdir = os.path.join(out, "gll_iem_v07")
    os.makedirs(gdir, exist_ok=True)
    shutil.copyfile(iem_path, os.path.join(gdir, "total.fits"))

    iem, _ = FC.read_iem(iem_path, emid)
    lv = FR.levels(FR.load_rings(rings, nb, emid), iem, nb)["rings_species"]
    pdir = os.path.join(out, "pohl22_rings_species")
    for name in lv[0]:
        cube = np.stack([lv[e][name] for e in range(nb)])
        CI.write_cube(cube, emid * 1e3, os.path.join(pdir, f"{name}.fits"),
                      {"MODEL": "pohl22_rings_species", "COMP": name,
                       "EBINSRC": (os.path.basename(ccube)[:60],
                                   "planes are this cube's bin centres")})
    print(f"wrote {gdir} (1 component) and {pdir} ({len(lv[0])} components: "
          f"{', '.join(lv[0])})")
    return lv


def check_reference(ccube, expcube, catalog, iem_path, rings, out, mask="B"):
    import fit_likelihood as FL
    import ps_template as PS
    import regrid as RG
    counts, e_lo, e_hi, chdr = FC.read_ccube(ccube)
    emid = np.sqrt(e_lo * e_hi)
    nb = len(emid)
    exposure, _ = FC.read_expcube(expcube, emid, chdr)
    iem, _ = FC.read_iem(iem_path, emid)
    m, _ = FC.build_masks(catalog, emid, mask)
    flat = RG.normalise(np.ones((400, 400)))
    bub, _ = FC.load_template(os.path.join(ROOT, "templates", "gap1_bubbles_fuzzy.fits"))
    ps, _ = PS.build(catalog, e_lo, e_hi)
    base = [{"iso": flat, "bubbles": bub, "ps": ps[e]} for e in range(nb)]

    lv = FR.levels(FR.load_rings(rings, nb, emid), iem, nb)["rings_species"]
    a = FL.fit(counts, exposure, [{**base[e], **lv[e]} for e in range(nb)], emid, m)

    pdir = os.path.join(out, "pohl22_rings_species")
    cubes = {f[:-5]: FC.read_iem(os.path.join(pdir, f), emid)[0]
             for f in sorted(os.listdir(pdir)) if f.endswith(".fits")}
    comp = [{f"iem_{k}": RG.normalise(np.maximum(c[e], 0.0)) for k, c in cubes.items()}
            for e in range(nb)]
    b = FL.fit(counts, exposure, [{**base[e], **comp[e]} for e in range(nb)], emid, m)
    d = abs(a["lnL_total"] - b["lnL_total"])
    print(f"Pohl22 rings_species at Level 1, mask {mask}: lnL via rings "
          f"{a['lnL_total']:,.4f}, via components {b['lnL_total']:,.4f}, "
          f"difference {d:.2e}")
    ok = d < 1e-3
    print("REFERENCE CHECK PASS" if ok else "REFERENCE CHECK FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccube", required=True)
    ap.add_argument("--iem", required=True)
    ap.add_argument("--rings", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "gap5", "components"))
    ap.add_argument("--check", action="store_true",
                    help="also confirm the two paths agree (needs --expcube, --catalog)")
    ap.add_argument("--expcube", default="")
    ap.add_argument("--catalog", default="")
    a = ap.parse_args()
    export(a.ccube, a.iem, a.rings, a.out)
    if a.check:
        sys.exit(0 if check_reference(a.ccube, a.expcube, a.catalog, a.iem,
                                      a.rings, a.out) else 1)


if __name__ == "__main__":
    main()
