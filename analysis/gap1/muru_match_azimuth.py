"""Find the azimuth zero point and the halo identity at the same time.

Matching my reconstruction to Muru's map at the SAME quoted viewing angle
assumes his in-plane reference vector v_disk is the AHF major axis Ea. Nothing
says it is: v_disk comes from HestiaUtils, which is not in the repository. If
it is some other in-plane direction then my angle 0 is his angle 0 + delta,
and a same-angle comparison fails even when everything else is right.

So the azimuth is left free. For each of his maps, scan mine over 0-355 deg and
keep the best shape correlation. A real match shows up as one delta per galaxy
that is the same at all four of his quoted angles, with his angle steps
reproduced by equal steps in mine.
"""

import os
import sys

import numpy as np
from astropy.io import fits

import muru_match as mm
from muru_match_run import GALAXIES, ANGLES, load
from muru_match_check import published, regrid, shape, corr

STEP = 5
SCAN = range(0, 360, STEP)


def main():
    reduced, profiles = sys.argv[1], sys.argv[2]
    tdir = sys.argv[3] if len(sys.argv) > 3 else "templates"

    pub = {(g, a): shape(published(tdir, g, a)) for g in GALAXIES for a in ANGLES}

    mine = {}
    for g in GALAXIES:
        pos, mass, frame = load(reduced, profiles, g)
        rho = mm.knn_density(pos, mass)
        normal, disk = frame["Ec"], frame["Ea"]
        disk = disk - normal * float(disk @ normal)
        fv = {"normal": normal, "disk": disk / np.linalg.norm(disk)}
        for az in SCAN:
            mp = mm.maps_for(pos, mass, rho, fv, az, cone=1.0)
            mine[(g, az)] = shape(regrid(mp["col_cone"] ** 2)).astype(np.float32)
        print(f"scanned {g}")

    print("\nbest match for each of Muru's maps, azimuth free")
    print("  his map        my halo   my azimuth    r      r on the"
          "   second best")
    print("                                               diagonal")
    for g in GALAXIES:
        for a in ANGLES:
            best, second = None, None
            for h in GALAXIES:
                rs = [corr(mine[(h, az)], pub[(g, a)]) for az in SCAN]
                k = int(np.argmax(rs))
                cand = (rs[k], h, list(SCAN)[k])
                if best is None or cand[0] > best[0]:
                    best, second = cand, best
                elif second is None or cand[0] > second[0]:
                    second = cand
            diag = max(corr(mine[(g, az)], pub[(g, a)]) for az in SCAN)
            print(f"  {g} ang{a:<3d}      {best[1]}     {best[2]:3d} deg "
                  f"  {best[0]:.3f}    {diag:.3f}      "
                  f"{second[1]} {second[2]:3d} ({second[0]:.3f})")

    print("\nhow well the best azimuth tracks his angle steps")
    for g in GALAXIES:
        for h in GALAXIES:
            azs = []
            for a in ANGLES:
                rs = [corr(mine[(h, az)], pub[(g, a)]) for az in SCAN]
                azs.append(list(SCAN)[int(np.argmax(rs))])
            d = np.diff(np.unwrap(np.radians(azs)) * 180 / np.pi)
            if np.all(np.abs(np.abs(d) - 30) < 15):
                print(f"  his {g} <- my {h}: azimuths {azs}, steps "
                      f"{np.round(d).astype(int).tolist()}  CONSISTENT")


if __name__ == "__main__":
    main()
