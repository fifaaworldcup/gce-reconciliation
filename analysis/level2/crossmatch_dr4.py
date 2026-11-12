"""Cross-match the Level 2 sources against 4FGL-DR4 (gll_psc_v35.fit).

Rule: every DR4 source within 0.3 deg is a candidate. One candidate is a
match. Two or more, or none inside 0.3 deg but one within 0.4 deg, is
uncertain, with the reason written out. Nothing within 0.4 deg is no match.

    python analysis/level2/crossmatch_dr4.py --catalog gll_psc_v35.fit \
        --sources outputs/gap1/l2_source_record.csv \
        --out outputs/gap1/l2_source_crossmatch_dr4.csv
"""
import argparse
import csv

import numpy as np
from astropy.io import fits

RADIUS = 0.3
NEAR_MISS = 0.4


def sep_deg(l1, b1, l2, b2):
    l1, b1, l2, b2 = map(np.radians, (l1, b1, l2, b2))
    c = np.sin(b1) * np.sin(b2) + np.cos(b1) * np.cos(b2) * np.cos(l1 - l2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def label(cls):
    cls = cls.strip()
    if not cls:
        return "unassociated"
    return {"psr": "pulsar", "msp": "millisecond pulsar", "snr": "supernova remnant",
            "unk": "unknown-type association", "spp": "SNR/PWN-type", "pwn": "pulsar wind nebula",
            "glc": "globular cluster", "bcu": "blazar", "fsrq": "blazar", "bll": "blazar",
            "hmb": "high-mass binary", "bin": "binary"}.get(cls.lower(), cls)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--sources", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    h = fits.open(a.catalog)
    if h[1].header.get("VERSION") != "v35":
        raise SystemExit(f"expected 4FGL-DR4 v35, got {h[1].header.get('VERSION')}")
    cat = h[1].data
    glon = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"])
    glat = cat["GLAT"]

    rows = []
    for s in csv.DictReader(open(a.sources)):
        l, b = float(s["l_deg"]), float(s["b_deg"])
        d = sep_deg(l, b, glon, glat)
        order = np.argsort(d)
        inside = [i for i in order if d[i] <= RADIUS]
        near = order[0]

        def desc(i):
            return (f"{cat['Source_Name'][i].strip()} ({d[i]:.3f} deg, "
                    f"{label(cat['CLASS1'][i])}, F1000={cat['Flux1000'][i]:.2e})")

        if len(inside) == 1:
            status, why = "match", ""
        elif len(inside) > 1:
            status = "uncertain"
            why = f"{len(inside)} DR4 sources within {RADIUS} deg: " + "; ".join(desc(i) for i in inside)
        elif d[near] <= NEAR_MISS:
            status = "uncertain"
            why = f"nearest DR4 source is just outside {RADIUS} deg: {desc(near)}"
        else:
            status, why = "no match", f"nearest DR4 source {d[near]:.2f} deg away"

        best = inside[0] if inside else near
        rows.append({
            "source": s["source"], "l_deg": l, "b_deg": b,
            "residual_TS": s["residual_TS"], "order_added": s["order_added"],
            "status": status,
            "n_dr4_within_0p3deg": len(inside),
            "dr4_name": cat["Source_Name"][best].strip(),
            "sep_deg": round(float(d[best]), 3),
            "dr4_class": cat["CLASS1"][best].strip(),
            "dr4_type": label(cat["CLASS1"][best]),
            "dr4_assoc": cat["ASSOC1"][best].strip(),
            "dr4_flux1000_ph_cm2_s": f"{cat['Flux1000'][best]:.3e}",
            "dr4_signif": round(float(cat["Signif_Avg"][best]), 1),
            "dr4_extended": bool(cat["Extended_Source_Name"][best].strip()),
            "on_galactic_centre_1deg": s.get("on_galactic_centre_1deg", ""),
            "near_roi_edge_1deg": s.get("near_roi_edge_1deg", ""),
            "reason": why,
        })

    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['source']}  {r['status']:9s} {r['dr4_name']:20s} {r['sep_deg']:.3f}  "
              f"{r['dr4_type']:24s} {r['dr4_assoc']:18s} {r['reason']}")


if __name__ == "__main__":
    main()
