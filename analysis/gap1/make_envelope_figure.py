"""
Figure 4: the HESTIA envelope. dAIC of every HESTIA galaxy x viewing angle
against the VVV bulge, in each fit variant, with the best-gamma gNFW marked for
scale. The point of the figure is the SPREAD - one halo alone told the wrong
story in the first fit.

Colours validated with the dataviz palette checker (blue/orange, light
surface: CVD dE 24.7, contrast >= 3:1). Angle is carried by marker shape, not
colour, so nothing depends on hue alone.
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
GALS = ["G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2"]
ANG = {"0": "o", "30": "s", "60": "^", "90": "D"}
SETS = {
    "masks": ([("scan_maskA", "Mask A (22% of sky kept)"),
               ("scan_maskB", "Mask B"),
               ("scan_maskB_ps", "Mask B + 4FGL point-source template")],
              "outputs/gap1/fig4_hestia_envelope.png",
              "HESTIA vs VVV depends on which simulated galaxy; every HESTIA halo "
              "beats the best spherical gNFW"),
    "masks2": ([("scan_maskD", "Mask D, sources masked only"),
                ("scan_maskD_ps", "Mask D + point-source template"),
                ("eb_scan_b1-10", "Mask B + point-source template")],
               "outputs/gap1/fig6_maskD.png",
               "gNFW loses under every treatment; HESTIA vs VVV does not survive "
               "the choice of how point sources are handled"),
    "energy": ([("eb_scan_b1-3", "1–3 GeV  (Mask B + PS)"),
                ("eb_scan_b3-10", "3–10 GeV  (Mask B + PS)"),
                ("eb_scan_all", "0.5–50 GeV  (Mask B + PS)")],
               "outputs/gap1/fig5_envelope_by_energy.png",
               "Same pattern in both bins where the GCE is detected"),
}
VARIANTS, OUT, TITLE = SETS[sys.argv[1] if len(sys.argv) > 1 else "masks"]


def load(d):
    return {r["template"]: float(r["dAIC_vs_vvv"])
            for r in csv.DictReader(open(f"outputs/gap1/{d}/scan.csv"))}


fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), facecolor=SURF)
for ax, (d, title) in zip(axes, VARIANTS):
    R = load(d)
    ax.set_facecolor(SURF)
    g_best = min((k for k in R if k.startswith("gnfw_") and not k.endswith("_th1")),
                 key=lambda k: R[k])
    xs = [R[k] for k in R if k.startswith("hestia")] + [R[g_best], 0.0]
    lo, hi = min(xs), max(xs)
    pad = 0.08 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad * 2.2)
    ax.axvspan(-10, 10, color=GRID, zorder=0, lw=0)
    ax.axvline(0, color=INK2, lw=1, zorder=1)
    ax.axvline(R[g_best], color=ORANGE, lw=2, ls=(0, (4, 3)), zorder=2)
    ax.text(R[g_best], -0.55, f"best gNFW γ={g_best.split('_g')[1]} ",
            color=INK2, fontsize=8, va="bottom", ha="right")
    for y, g in enumerate(GALS):
        for a, m in ANG.items():
            k = f"hestia_{g}_ang{a}"
            ax.plot(R[k], y, m, ms=8, mfc=BLUE, mec=SURF, mew=1.5, zorder=3)
    ax.set_yticks(range(len(GALS)))
    ax.set_yticklabels(GALS if ax is axes[0] else [])
    ax.set_ylim(len(GALS) - 0.4, -0.95)
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_xlabel("ΔAIC vs VVV bulge", fontsize=9, color=INK2)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(axis="x", color=GRID, lw=0.6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.01, -0.13, "← favours HESTIA", transform=ax.transAxes,
            fontsize=8, color=INK2)
    ax.text(0.99, -0.13, "favours VVV →", transform=ax.transAxes,
            fontsize=8, color=INK2, ha="right")

axes[0].set_ylabel("HESTIA galaxy", fontsize=9, color=INK2)
handles = [plt.Line2D([], [], ls="", marker=m, ms=8, mfc=BLUE, mec=SURF,
                      label=f"{a}° view") for a, m in ANG.items()]
handles.append(plt.Rectangle((0, 0), 1, 1, color=GRID, label="|ΔAIC| < 10, not decisive"))
fig.legend(handles=handles, loc="upper left", ncol=5, fontsize=8, frameon=False,
           bbox_to_anchor=(0.005, 0.93))
fig.suptitle(TITLE, fontsize=11, color=INK, x=0.01, ha="left", y=0.99)
fig.tight_layout(rect=[0, 0.02, 1, 0.86])
out = OUT
fig.savefig(out, dpi=160, facecolor=SURF)
print("wrote", out)
