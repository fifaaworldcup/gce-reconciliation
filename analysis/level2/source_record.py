"""The record of what Level 2 added, with every mechanical column filled in.

Afeefa asked for a form someone else can complete. Everything a machine can
decide is decided here - position, residual TS, pass, the nearest 4FGL source
and whether it sits inside the PSF, whether that source is extended, and
whether the candidate is near the mask boundary or the ROI edge. What is left
is the judgement: catalogued, plausibly real, or likely an artefact, and one
line saying why.

Leaving the mechanical columns blank would have handed a reviewer four hours
of arithmetic and no judgement. This way the four hours are all judgement.

  python3 analysis/level2/source_record.py \
      --level2 outputs/gap1/l2_neutral_maskB_ps/level2.json \
      --catalog <4FGL file> --out outputs/gap1/l2_source_record.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (os.path.join(ROOT, "utils"), os.path.join(ROOT, "masks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import psf as PSF              # noqa: E402

HALF = 20.0
EXTENDED_FLAG_COLS = ("Extended_Source_Name", "Source_Name")


def read_4fgl(path):
    with fits.open(path) as hd:
        t = hd[1].data
        names = [c.name for c in hd[1].columns]
        glon = np.asarray(t["GLON"], float)
        glon = np.where(glon > 180, glon - 360, glon)
        out = {"glon": glon, "glat": np.asarray(t["GLAT"], float),
               "name": np.asarray(t["Source_Name"], str),
               "ts": np.asarray(t["Signif_Avg"], float) ** 2
               if "Signif_Avg" in names else np.full(len(glon), np.nan)}
        out["ext"] = (np.asarray(t["Extended_Source_Name"], str)
                      if "Extended_Source_Name" in names
                      else np.full(len(glon), "", dtype=object))
        out["class"] = (np.asarray(t["CLASS1"], str) if "CLASS1" in names
                        else np.full(len(glon), "", dtype=object))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level2", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--energy-gev", type=float, default=1.15,
                    help="energy at which the PSF containment radius is "
                         "evaluated. The default is the lowest bin, where the "
                         "PSF is widest, so 'inside the PSF' is the generous "
                         "reading and a source called uncatalogued here is "
                         "uncatalogued at every energy.")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    j = json.load(open(a.level2))
    cat = read_4fgl(a.catalog)
    r95 = PSF.containment_radius(0.95, a.energy_gev)
    r68 = PSF.containment_radius(0.68, a.energy_gev)

    rows = []
    for s in j["sources"]:
        l, b = float(s["l"]), float(s["b"])
        # angular separation to every catalogue source, small-angle in l
        d = np.hypot((cat["glon"] - l) * np.cos(np.radians(b)),
                     cat["glat"] - b)
        i = int(np.argmin(d))
        sep = float(d[i])
        edge = min(HALF - abs(l), HALF - abs(b))
        rows.append({
            "source": s["name"],
            "l_deg": round(l, 2), "b_deg": round(b, 2),
            "residual_TS": round(float(s["residual_ts"]), 1),
            "pass_threshold": s["pass"],
            "order_added": s["index"],
            "nearest_4fgl": cat["name"][i],
            "nearest_4fgl_sep_deg": round(sep, 3),
            "nearest_4fgl_class": str(cat["class"][i]).strip(),
            "nearest_4fgl_ts": ("" if not np.isfinite(cat["ts"][i])
                                else round(float(cat["ts"][i]), 1)),
            "inside_psf68": bool(sep <= r68),
            "inside_psf95": bool(sep <= r95),
            "nearest_is_extended": bool(str(cat["ext"][i]).strip() not in ("", "--")),
            "deg_from_roi_edge": round(float(edge), 2),
            "near_roi_edge_1deg": bool(edge < 1.0),
            "on_galactic_centre_1deg": bool(np.hypot(l, b) < 1.0),
            # ---- for a human to fill in ----
            "verdict": "",          # catalogued | plausibly_real | likely_artefact
            "note": "",
        })
    rows.sort(key=lambda r: -r["residual_TS"])

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    print(f"wrote {a.out}: {len(rows)} sources\n")
    print(f"PSF at {a.energy_gev:g} GeV: 68% containment {r68:.2f} deg, "
          f"95% {r95:.2f} deg\n")
    print(f"  {'source':8s} {'l':>7} {'b':>7} {'TS':>7} {'nearest 4FGL':>20} "
          f"{'sep':>7} {'in 95%':>7}  flags")
    for r in rows:
        flags = []
        if r["on_galactic_centre_1deg"]:
            flags.append("GALACTIC CENTRE")
        if r["near_roi_edge_1deg"]:
            flags.append("near ROI edge")
        if r["nearest_is_extended"]:
            flags.append("nearest is extended")
        print(f"  {r['source']:8s} {r['l_deg']:7.2f} {r['b_deg']:7.2f} "
              f"{r['residual_TS']:7.1f} {r['nearest_4fgl']:>20} "
              f"{r['nearest_4fgl_sep_deg']:7.3f} "
              f"{'yes' if r['inside_psf95'] else 'NO':>7}  "
              f"{', '.join(flags)}")
    n_in = sum(r["inside_psf95"] for r in rows)
    print(f"\n  {n_in} of {len(rows)} fall inside the 95% PSF of a catalogued "
          f"source; {len(rows) - n_in} do not and need a judgement.")
    print("  Fill `verdict` with catalogued | plausibly_real | likely_artefact "
          "and one line in `note`.")


if __name__ == "__main__":
    main()
