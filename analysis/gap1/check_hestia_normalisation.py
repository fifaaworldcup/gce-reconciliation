"""Do Muru's maps carry absolute J-factors, or only shapes?

Decides whether the halos can be ranked by mass (and so whether MW and M31
analogues can be told apart without asking), and whether Gap 4 or Gap 7 can
get an absolute annihilation flux from these files.

Sum each map. If the sums agree across halos, simulations and smoothing, the
normalisation was imposed and the absolute information is gone.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits


def main(root="~/muru_data/jf"):
    root = os.path.expanduser(root)
    files = sorted(glob.glob(os.path.join(root, "density_projection_*_angle0.0_dm2_*deg.fits")))
    if not files:
        raise SystemExit(f"no maps under {root}")
    rows = []
    for f in files:
        m = re.search(r"projection_(\d+_\d+)_8192_halo(\d+)_angle[\d.]+_dm2_([\d.]+)deg", f)
        if not m:
            continue
        d = np.nan_to_num(fits.getdata(f).astype(float))
        rows.append((m.group(1), m.group(2)[-4:], m.group(3), d.sum(), d.max()))

    print(f"{'sim':8} {'halo':>6} {'smooth':>7} {'total':>16} {'peak':>10}")
    for r in sorted(rows):
        print(f"{r[0]:8} {r[1]:>6} {r[2]:>7} {r[3]:16.6f} {r[4]:10.4g}")

    tot = np.array([r[3] for r in rows])
    pk = np.array([r[4] for r in rows])
    spread = (tot.max() - tot.min()) / tot.mean()
    print(f"\n{len(rows)} maps")
    print(f"total spread across all of them : {spread:.3e}")
    print(f"peak  spread across all of them : "
          f"{(pk.max() - pk.min()) / pk.mean():.3f}  "
          f"(factor {pk.max() / pk.min():.2f} between the most and least "
          f"concentrated)")
    if spread < 1e-5:
        print("\nVERDICT: the maps are normalised to a common total. They carry")
        print("SHAPE ONLY - no absolute J-factor, and no way to rank the halos")
        print("by mass, so MW vs M31 cannot be settled from these files.")
        print("Anyone downstream who needs absolute flux must ask Muru for the")
        print("unnormalised maps.")
        return 0
    print("\nVERDICT: totals differ - absolute information is present.")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
