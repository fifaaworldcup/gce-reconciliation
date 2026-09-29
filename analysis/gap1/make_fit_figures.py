"""
Figures for mini-paper section 6, built from fit_cubes.py output.

Three panels per signal model: the fitted signal map, the residual
significance map, and the residual with the signal component removed - i.e.
what the GCE looks like to the fit. Plus a summary bar of delta ln L.

Residual maps use a diverging colormap centred exactly on zero with symmetric
limits, because the eye reads the colour midpoint as "no excess". Clipping at
+-5 sigma rather than autoscaling keeps the panels comparable between models,
which is the entire point of putting them side by side.
"""

import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

HALF = 20.0
EXTENT = [HALF, -HALF, -HALF, HALF]      # l runs the other way, as on the sky
SIGMA_CLIP = 5.0


def _sky(ax, im, title, cmap, vmin, vmax, cbar_label):
    h = ax.imshow(im, origin="lower", extent=EXTENT, cmap=cmap,
                  vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("$\\ell$ [deg]", fontsize=9)
    ax.set_ylabel("$b$ [deg]", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.axhline(0, color="0.4", lw=0.4, alpha=0.6)
    ax.axvline(0, color="0.4", lw=0.4, alpha=0.6)
    cb = plt.colorbar(h, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitdir", default="outputs/gap1/fit")
    ap.add_argument("--out", default="outputs/gap1/fig3_fit_residuals.png")
    args = ap.parse_args()

    with open(os.path.join(args.fitdir, "fit_results.json")) as f:
        rec = json.load(f)

    names = list(rec["models"])
    n = len(names)
    fig, axes = plt.subplots(n, 2, figsize=(9.0, 3.7 * n), squeeze=False)

    for i, name in enumerate(names):
        m = rec["models"][name]
        model = fits.getdata(os.path.join(args.fitdir, f"model_{name}.fits"))
        resid = np.asarray(fits.getdata(os.path.join(args.fitdir,
                                                     f"residual_{name}.fits")),
                           float)

        mx = float(np.nanmax(model))
        _sky(axes[i][0], model / mx if mx > 0 else model,
             f"{name}: fitted signal map", "magma", 0.0, 1.0,
             "fraction of peak")

        _sky(axes[i][1], resid,
             f"{name}: residual  (TS {m['TS']:,.0f}, "
             f"{m['residual_std_sigma']:.2f}$\\sigma$ scatter)",
             "RdBu_r", -SIGMA_CLIP, SIGMA_CLIP,
             "residual [$\\sigma$]")

    fig.suptitle("Gap 1 first fit: signal templates and residuals, "
                 f"{rec['inputs']['bins_gev'][0][0]:.0f}"
                 f"-{rec['inputs']['bins_gev'][-1][1]:.0f} GeV, "
                 f"Mask {rec['inputs']['mask']['mask']}"
                 + (" + 4FGL point-source template"
                    if "ps" in rec.get("components", {}) else ""),
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

    # ---- comparison summary table for section 6
    lines = ["| model | TS | sigma | ln L | photons | residual sigma |",
             "|---|---|---|---|---|---|"]
    for name in names:
        m = rec["models"][name]
        lines.append(f"| {name} | {m['TS']:,.1f} | {m['sigma']:.1f} | "
                     f"{m['lnL_total']:,.1f} | {m['model_photons']:,.0f} | "
                     f"{m['residual_mean_sigma']:+.3f} +- "
                     f"{m['residual_std_sigma']:.3f} |")
    lines += ["", "| comparison | dlnL | dAIC | favoured | decisive? |",
              "|---|---|---|---|---|"]
    for k, v in rec["comparisons"].items():
        lines.append(f"| {k.replace('_minus_', ' - ')} | {v['delta_lnL']:+.2f} | "
                     f"{v['delta_AIC']:+.2f} | {v['favoured']} | "
                     f"{'yes' if v['decisive_at_dAIC10'] else 'no'} |")
    if rec.get("caveats"):
        lines += ["", "**Caveats**", ""]
        lines += [f"- {c}" for c in rec["caveats"]]

    tbl = os.path.join(os.path.dirname(args.out), "section6_fit_tables.md")
    with open(tbl, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {tbl}")
    print()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
