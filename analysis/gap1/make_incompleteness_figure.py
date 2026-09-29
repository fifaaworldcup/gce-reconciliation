"""Figure 7: catalogue-incompleteness calibration curve."""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

C = {"hestia": "#1f4e79", "spherical": "#a6402a"}
LABEL = {"hestia": "flattened injection (HESTIA G1.1)",
         "spherical": "spherical injection (gNFW $\\gamma$=1.2)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", required=True)
    ap.add_argument("--q-measured", type=float, default=0.60,
                    help="the flattening section 6.1 reports")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S = json.load(open(a.calib))["summary"]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))

    ax = axes[0]
    xmax = 0
    for name, s in S.items():
        p = [q for q in s["points"] if q["catalog_frac"] < 1.0]
        x = np.array([q["f_unres_gce"] for q in p])
        y = np.array([q["bias"] for q in p])
        e = np.array([q["bias_err"] for q in p])
        xmax = max(xmax, x.max())
        ax.errorbar(x, y, yerr=e, fmt="o", ms=4.5, lw=0, elinewidth=1.1,
                    capsize=2, color=C[name], label=LABEL[name])
        cal = s.get("calibration") or {}
        if cal:
            xc = np.linspace(0, cal["f_max"], 50)
            ax.plot(xc, -cal["slope_k"] * xc, "-", lw=1.4, color=C[name],
                    alpha=0.85)
        if s.get("power_law"):
            pw = s["power_law"]
            xs = np.linspace(0.01, x.max() * 1.05, 120)
            ax.plot(xs, -pw["A"] * xs ** pw["p"], "--", lw=1.0,
                    color=C[name], alpha=0.5)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlabel("flux in sources the fit was not given  /  GCE flux")
    ax.set_ylabel("bias in recovered axis ratio  $q - q_{\\rm ref}$")
    ax.set_title("Hiding faint sources makes the excess look flattened",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower left", frameon=False)
    ax.set_xlim(0, xmax * 1.05)
    ax.text(0.98, 0.95, "solid: calibration, $f\\leq1.2$   "
            "dashed: $-Af^{\\,p}$, whole sweep",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            color="0.35")

    ax = axes[1]
    A = json.load(open(a.calib)).get("attractor", [])
    if A:
        A = sorted(A, key=lambda p: p["f_unres_gce"])
        x = np.array([max(p["f_unres_gce"], 0.008) for p in A])
        qh = [p["q_hestia"] for p in A]
        qs_ = [p["q_spherical"] for p in A]
        ax.plot(x, qh, "o-", ms=4, lw=1.4, color=C["hestia"],
                label="flattened injection")
        ax.plot(x, qs_, "s-", ms=4, lw=1.4, color=C["spherical"],
                label="spherical injection")
        ax.fill_between(x, qh, qs_, color="0.6", alpha=0.18, lw=0)
        ax.annotate("", xy=(x[0], qs_[0]), xytext=(x[0], qh[0]),
                    arrowprops=dict(arrowstyle="<->", color="0.35", lw=0.9))
        ax.text(x[0] * 1.15, 0.5 * (qh[0] + qs_[0]),
                f"truth differs\nby {qs_[0] - qh[0]:.2f}", fontsize=7.5,
                color="0.3", va="center")
        far = [p for p in A if p["f_unres_gce"] > 2.0]
        if far:
            mid = float(np.mean([p["midpoint"] for p in far]))
            ax.axhline(mid, color="0.25", lw=0.9, ls="--")
            ax.text(x.max(), mid - 0.004, f"both read q = {mid:.2f} ",
                    ha="right", va="top", fontsize=8, color="0.2")
    ax.axhline(a.q_measured, color="#2e7d32", lw=1.2, ls=":")
    ax.text(0.99, 0.055, f"section 6.1 measures q = {a.q_measured:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="#2e7d32")
    ax.set_xlabel("flux in sources the fit was not given  /  GCE flux")
    ax.set_ylabel("recovered axis ratio $q$")
    ax.set_title("Two very different truths converge on one answer",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.set_xscale("log")
    ax.set_xlim(0.007, 7.0)
    ax.set_ylim(0.50, 0.95)

    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print("wrote " + a.out)


if __name__ == "__main__":
    main()
