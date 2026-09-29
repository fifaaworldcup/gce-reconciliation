"""Figure 8: TS against how much the diffuse model is allowed to bend."""
from __future__ import annotations
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ORDER = ["rigid", "ns", "lat_outside", "lat_inside", "lat", "lat_ns"]
XLAB = ["1 norm\n(as in §6.2)", "+ N/S\n(1 extra)",
        "|b| bands, but only\nwhere the GCE isn't",
        "|b| bands, only\nwhere the GCE is",
        "|b| bands\neverywhere",
        "bands × N/S\neverywhere"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    n = len(a.json)
    fig, axes = plt.subplots(1, n, figsize=(max(4.0 * n, 7.2), 5.0), squeeze=False)
    for k, (path, lab) in enumerate(zip(a.json, a.labels)):
        d = json.load(open(path))
        ax = axes[0][k]
        x = np.arange(len(ORDER))
        sp = d["splits"]
        real = [sp[s]["observed"]["vvv"]["ts"] for s in ORDER]
        injd = [sp[s]["injected"]["vvv"]["ts"] for s in ORDER]
        gain = [sp[s]["boosted"]["vvv"]["ts"] - sp[s]["observed"]["vvv"]["ts"]
                for s in ORDER]
        bad = {r["split"] for r in d["verdict"]["rows"] if not r["usable"]}
        ax.plot(x, real, "o-", lw=2.0, ms=6, color="#a6402a",
                label="observed excess")
        ax.plot(x, injd, "s--", lw=1.6, ms=5, color="#1f4e79",
                label="signal we injected ourselves")
        ax.plot(x, gain, "^:", lw=1.4, ms=5, color="#2e7d32",
                label="extra signal added to the real sky")
        for j, nm in enumerate(ORDER):
            if nm in bad:
                ax.axvspan(j - 0.35, j + 0.35, color="0.85", zorder=0)
                ax.text(j, 1.3, "degenerate\nwith the signal", ha="center",
                        va="bottom", fontsize=6.5, color="0.35")
        ax.axhline(25, color="0.5", lw=0.8, ls="--")
        ax.text(len(ORDER) - 1, 26, "TS 25 ", ha="right", va="bottom",
                fontsize=7.5, color="0.4")
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xticks(x)
        ax.set_xticklabels(XLAB, fontsize=6.5, rotation=30,
                           ha="right")
        ax.set_title(lab, fontsize=10)
        if k == 0:
            ax.set_ylabel("TS of the VVV bulge template")
            ax.legend(fontsize=8, frameon=False, loc="lower left")
        ax.set_xlabel("how the diffuse model is allowed to bend")
        ax.set_ylim(-1, 1200)
    fig.suptitle("The detection spans a factor of eighty across\nbackground "
                 "models we cannot choose between", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(a.out, dpi=150)
    print("wrote " + a.out)


if __name__ == "__main__":
    main()
