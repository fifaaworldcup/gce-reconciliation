"""Does the matched reconstruction reproduce the maps Muru actually sent?

A straight Pearson r between two of these maps is 0.88-0.99 no matter which
two, because every map is a bright central cusp and the cusp carries no
morphology. So the test here strips the cusp: take log of the map, subtract
its own azimuthal average as a function of angular radius from the centre,
and correlate what is left. That residual IS the shape -- the elongation and
the boxiness -- and it is what the templates are being fitted for.

Two labellings have to agree. The particle files were identified by matching
the dark-matter centroid to the AHF centres in the data readme, because their
file names are swapped within each pair. The map files were identified by the
simulation and halo ID inside their names. The diagonal of the matrix below
has to win, and the l-mirrored arm says whether Muru's backwards longitude
(e_y = e_x cross e_z, the opposite handedness to galactic) survived into the
templates.
"""

import os
import sys

import numpy as np
from astropy.io import fits

import muru_match as mm
from muru_match_run import GALAXIES, ANGLES, load


def published(template_dir, galaxy, angle):
    path = os.path.join(template_dir, f"gap1_hestia_{galaxy}_ang{angle}_1deg.fits")
    with fits.open(path) as hd:
        return np.asarray(hd[0].data, float).T          # [b,l] -> [l,b]


def regrid(m):
    return 0.25 * (m[:-1, :-1] + m[1:, :-1] + m[:-1, 1:] + m[1:, 1:])


def shape(m, nring=40, rmax=15.0):
    """log map minus its own azimuthal average, inside rmax degrees."""
    n = m.shape[0]
    ax = (np.arange(n) - (n - 1) / 2.0) * (40.0 / n)
    r = np.hypot(ax[:, None], ax[None, :])
    v = np.log(np.maximum(m, np.nanmax(m) * 1e-8))
    out = np.full_like(v, np.nan)
    edges = np.linspace(0.0, rmax, nring + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (r >= lo) & (r < hi)
        if sel.sum() > 3:
            out[sel] = v[sel] - v[sel].mean()
    return out


def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    x, y = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    return float(x @ y / np.sqrt((x @ x) * (y @ y)))


def main():
    reduced, profiles = sys.argv[1], sys.argv[2]
    tdir = sys.argv[3] if len(sys.argv) > 3 else "templates"

    mine, mine_mirror, qmine = {}, {}, {}
    for g in GALAXIES:
        pos, mass, frame = load(reduced, profiles, g)
        rho = mm.knn_density(pos, mass)
        normal, disk = frame["Ec"], frame["Ea"]
        disk = disk - normal * float(disk @ normal)
        fv = {"normal": normal, "disk": disk / np.linalg.norm(disk)}
        for ang in ANGLES:
            mp = mm.maps_for(pos, mass, rho, fv, ang, cone=1.0)
            sq = regrid(mp["col_cone"] ** 2)
            mine[(g, ang)] = shape(sq)
            mine_mirror[(g, ang)] = shape(sq[::-1])
            la, ba = mp["_axes"]
            qmine[(g, ang)] = mm.axis_ratio_svd(mp["col_cone"] ** 2, 0.5, la, ba)["q"]
        print(f"built {g}")

    pub, qpub = {}, {}
    for g in GALAXIES:
        for ang in ANGLES:
            p = published(tdir, g, ang)
            pub[(g, ang)] = shape(p)
            n = p.shape[0]
            ax = (np.arange(n) - (n - 1) / 2.0) * 0.1
            qpub[(g, ang)] = mm.axis_ratio_svd(p, 0.5, ax, ax)["q"]

    for name, src in (("as built", mine), ("mirrored in l", mine_mirror)):
        print(f"\nshape-only r, mine {name} (rows) vs Muru's template (columns)")
        total = 0
        for ang in ANGLES:
            hits = 0
            print(f"  angle {ang:3d}      " + "".join(f"{h:>8s}" for h in GALAXIES))
            for g in GALAXIES:
                rs = [corr(src[(g, ang)], pub[(h, ang)]) for h in GALAXIES]
                best = GALAXIES[int(np.argmax(rs))]
                hits += best == g
                tail = "  ok" if best == g else f"  best {best}"
                print(f"           {g:>7s}" + "".join(f"{r:8.3f}" for r in rs) + tail)
            print(f"           diagonal {hits}/6")
            total += hits
        print(f"  total diagonal {total}/24")

    print("\naxis ratio q50 of the squared column map, mine vs his")
    print("  galaxy    " + "".join(f"  ang{a:<3d}        " for a in ANGLES))
    for g in GALAXIES:
        row = "".join(f"  {qmine[(g, a)]:.3f}/{qpub[(g, a)]:.3f}  " for a in ANGLES)
        print(f"  {g:>7s} " + row)
    dm = np.array([qmine[(g, a)] for g in GALAXIES for a in ANGLES])
    dp = np.array([qpub[(g, a)] for g in GALAXIES for a in ANGLES])
    print(f"  mine mean {dm.mean():.3f}  his mean {dp.mean():.3f}"
          f"  offset {dm.mean() - dp.mean():+.3f}"
          f"  rms difference {np.sqrt(((dm - dp) ** 2).mean()):.3f}")
    from scipy.stats import spearmanr
    print(f"  Spearman between them rho {spearmanr(dm, dp).statistic:+.3f}")


if __name__ == "__main__":
    main()
