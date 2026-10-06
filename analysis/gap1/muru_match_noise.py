"""How much of our reconstruction is shape and how much is Poisson noise?

The 1 deg cone holds only a few hundred dark-matter particles per sight line
in these files, so the shape residual -- the map minus its own azimuthal
average -- can be mostly shot noise. This measures that directly: the residual
power on the pixel scale, which is noise, against the residual power on the
few-degree scale, which is shape. It does the same to Muru's maps, which is
how we find out whether the 1.0deg in the file names he sent is the cone
radius or a smoothing applied afterwards.
"""

import sys

import numpy as np
from scipy.ndimage import gaussian_filter

import muru_match as mm
from muru_match_run import GALAXIES, ANGLES, load
from muru_match_check import published, regrid, shape, corr


def power(res, sigma_pix=5.0):
    """rms of the residual, and of the part of it that survives smoothing."""
    v = np.nan_to_num(res)
    m = np.isfinite(res).astype(float)
    sm = gaussian_filter(v, sigma_pix) / np.maximum(gaussian_filter(m, sigma_pix), 1e-6)
    ok = np.isfinite(res)
    return float(np.sqrt((v[ok] ** 2).mean())), float(np.sqrt((sm[ok] ** 2).mean()))


def main():
    reduced, profiles = sys.argv[1], sys.argv[2]
    tdir = sys.argv[3] if len(sys.argv) > 3 else "templates"

    print("residual rms, total vs smoothed on 0.5 deg; the ratio is how much")
    print("of the shape survives smoothing, so low means noise-dominated\n")
    print(f"  {'map':24s} {'total':>8s} {'smooth':>8s} {'kept':>7s}")
    for g in GALAXIES:
        p = published(tdir, g, 0)
        t, s = power(shape(p))
        print(f"  {'his ' + g + ' ang0':24s} {t:8.4f} {s:8.4f} {s / t:7.2f}")

    for g in GALAXIES:
        pos, mass, frame = load(reduced, profiles, g)
        rho = mm.knn_density(pos, mass)
        normal, disk = frame["Ec"], frame["Ea"]
        disk = disk - normal * float(disk @ normal)
        fv = {"normal": normal, "disk": disk / np.linalg.norm(disk)}
        for cone in (1.0, 3.0):
            mp = mm.maps_for(pos, mass, rho, fv, 0, cone=cone)
            t, s = power(shape(regrid(mp["col_cone"] ** 2)))
            print(f"  {'mine ' + g + f' cone{cone:g}':24s} "
                  f"{t:8.4f} {s:8.4f} {s / t:7.2f}"
                  f"   {mp['_nkept']:,} particles in view")

    print("\nfree-azimuth shape match at the 3 deg cone")
    SCAN = range(0, 360, 10)
    pub = {(g, a): shape(published(tdir, g, a)) for g in GALAXIES for a in ANGLES}
    mine = {}
    for g in GALAXIES:
        pos, mass, frame = load(reduced, profiles, g)
        rho = mm.knn_density(pos, mass)
        normal, disk = frame["Ec"], frame["Ea"]
        disk = disk - normal * float(disk @ normal)
        fv = {"normal": normal, "disk": disk / np.linalg.norm(disk)}
        for az in SCAN:
            mp = mm.maps_for(pos, mass, rho, fv, az, cone=3.0)
            mine[(g, az)] = shape(regrid(mp["col_cone"] ** 2)).astype(np.float32)
        print(f"  scanned {g}")
    print(f"\n  {'his map':14s} {'best halo':>10s} {'az':>5s} {'r':>7s} "
          f"{'r diagonal':>11s}")
    for g in GALAXIES:
        for a in ANGLES:
            cands = []
            for h in GALAXIES:
                rs = [corr(mine[(h, az)], pub[(g, a)]) for az in SCAN]
                k = int(np.argmax(rs))
                cands.append((rs[k], h, list(SCAN)[k]))
            cands.sort(reverse=True)
            diag = max(corr(mine[(g, az)], pub[(g, a)]) for az in SCAN)
            print(f"  {g} ang{a:<3d}   {cands[0][1]:>10s} {cands[0][2]:5d} "
                  f"{cands[0][0]:7.3f} {diag:11.3f}")


if __name__ == "__main__":
    main()
