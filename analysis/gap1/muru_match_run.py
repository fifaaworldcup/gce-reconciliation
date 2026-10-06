"""Run the matched-convention HESTIA reconstruction and the four correlations.

Three things come out of this:

  * q for every galaxy and viewing angle under Muru's conventions, for the
    column map, its square, and the two J-factor weightings. His axis-ratio
    estimator is used throughout, so his numbers and ours are the same
    quantity measured the same way.
  * the same estimator applied to the maps he actually sent us, which is the
    only honest way to ask whether our reconstruction reproduces his.
  * Spearman rho against the sec 6.1 dAIC table, four ways: squared column
    and J-factor, each under the current and the matched reconstruction.

  python3 analysis/gap1/muru_match_run.py --reduced DIR --profiles DIR
"""

import argparse
import glob
import json
import os
import re

import numpy as np
from scipy.stats import spearmanr

import muru_match as mm

GALAXIES = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
ANGLES = (0, 30, 60, 90)
CONES = (1.0, 3.0)
TAG = {g: g.replace(".", "") for g in GALAXIES}


def load(reduced, profiles, galaxy, r_keep=26.0):
    t = TAG[galaxy]
    pos = np.load(os.path.join(reduced, f"{t}_dm_pos.npy")).astype(np.float64)
    mass = np.load(os.path.join(reduced, f"{t}_dm_mass.npy")).astype(np.float64)
    near = (pos ** 2).sum(axis=1) < r_keep ** 2
    pos, mass = pos[near], mass[near]
    frame = mm.read_ahf_frame(os.path.join(profiles, f"profile_{t}_AHF.txt"))
    return pos, mass, frame


def published_q(template_dir):
    """Muru's own maps, measured with his own estimator.

    The files he sent are the 1 deg cone maps, already squared (dm2), so the
    unsquared column is the square root and both arms are available from one
    file.
    """
    from astropy.io import fits
    out = {}
    for path in sorted(glob.glob(os.path.join(template_dir,
                                              "gap1_hestia_*_ang*_1deg.fits"))):
        m = re.search(r"gap1_hestia_(G\d\.\d)_ang(\d+)_1deg\.fits$", path)
        if not m:
            continue
        with fits.open(path) as hd:
            arr = np.asarray(hd[0].data, float)
            cd = abs(float(hd[0].header["CDELT1"]))
            n = arr.shape[0]
        axis = (np.arange(n) - (n - 1) / 2.0) * cd
        arr = arr.T                       # file is [b, l]; estimator wants [l, b]
        row = {}
        for f in (0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02):
            row[f"sq_q{int(f * 100)}"] = mm.axis_ratio_svd(arr, f, axis, axis)
            row[f"col_q{int(f * 100)}"] = mm.axis_ratio_svd(
                np.sqrt(np.maximum(arr, 0)), f, axis, axis)
            row[f"sq_q{int(f * 100)}"]["touches_roi_edge"] = _edge(arr, f)
        out[f"{m.group(1)}_ang{int(m.group(2))}"] = row
    return out


def _edge(m, frac):
    """True if the super-level set reaches the map boundary, where the crop
    rather than the halo sets its shape (sec 6.1, Table 2)."""
    s = m >= float(np.nanmax(m)) * frac
    return bool(s[0].any() or s[-1].any() or s[:, 0].any() or s[:, -1].any())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reduced", required=True)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--templates", default="templates")
    ap.add_argument("--scan", default="outputs/gap1/scan_maskB_ps/scan.csv")
    ap.add_argument("--current", default="")
    ap.add_argument("--out", default="outputs/gap1/muru_match.json")
    ap.add_argument("--normal", default="Ec", choices=("Ec", "L"))
    a = ap.parse_args()

    res = {"conventions": {
        "r_sun_kpc": mm.R_SUN, "d_max_from_observer_kpc": mm.D_MAX,
        "grid_half_deg": mm.HALF, "grid_step_deg": mm.STEP,
        "cones_deg": list(CONES), "ahf_profile_row_kpc_over_h": mm.AHF_R_HINV,
        "disc_normal": a.normal, "knn": mm.KNN,
        "softening_kpc": mm.SOFTENING,
        "axis_ratio": "min/max singular value of super-level pixels, unweighted",
    }, "frames": {}, "rows": {}}

    for g in GALAXIES:
        pos, mass, frame = load(a.reduced, a.profiles, g)
        rho = mm.knn_density(pos, mass)
        normal = frame[a.normal]
        disk = frame["Ea"]
        # Ea need not be exactly in the disc plane; project it there, which is
        # what a disc azimuth means.
        disk = disk - normal * float(disk @ normal)
        frame_v = {"normal": normal, "disk": disk / np.linalg.norm(disk)}
        res["frames"][g] = {
            "ahf_r_kpc": frame["r_kpc"], "ahf_r_hinv": frame["r_hinv"],
            "ahf_npart": frame["npart"],
            "ahf_b_over_a": frame["b_over_a"], "ahf_c_over_a": frame["c_over_a"],
            "normal_vs_angular_momentum_deg": frame["normal_vs_L_deg"],
            "n_dm_within_26kpc": int(pos.shape[0]),
            "rho_median": float(np.median(rho)),
        }
        print(f"{g}: AHF row r = {frame['r_kpc']:.2f} kpc "
              f"({frame['npart']:,} particles), c/a {frame['c_over_a']:.3f}, "
              f"minor axis {frame['normal_vs_L_deg']:.1f} deg from L, "
              f"{pos.shape[0]:,} DM within 26 kpc")

        for cone in CONES:
            for ang in ANGLES:
                maps = mm.maps_for(pos, mass, rho, frame_v, ang, cone=cone)
                key = f"{g}_ang{ang}_cone{cone:g}"
                res["rows"][key] = {"n_particles_in_view": maps["_nkept"],
                                    **mm.summarise(maps)}
                s = res["rows"][key]
                print(f"   cone {cone:.0f} deg  angle {ang:3d}:  "
                      f"col {s['col_cone']['q50']['q']:.3f}  "
                      f"col^2 {s['col_cone_sq']['q50']['q']:.3f}  "
                      f"J {s['j_cone']['q50']['q']:.3f}  "
                      f"colLOS {s['col_los']['q50']['q']:.3f}  "
                      f"JLOS {s['j_los']['q50']['q']:.3f}")

    res["published"] = published_q(a.templates)

    # ---- the four correlations -------------------------------------------
    daic = {}
    with open(a.scan) as fh:
        head = fh.readline().strip().split(",")
        ic = head.index("dAIC_vs_vvv")
        for line in fh:
            p = line.strip().split(",")
            m = re.match(r"hestia_(G\d\.\d)_ang(\d+)$", p[0])
            if m:
                daic[f"{m.group(1)}_ang{int(m.group(2))}"] = float(p[ic])

    current = {}
    if a.current and os.path.exists(a.current):
        with open(a.current) as fh:
            cur = json.load(fh)
        for row in cur["rows"]:
            g = row["galaxy"]
            g = g if "." in g else f"{g[:2]}.{g[2:]}"
            current[f"{g}_ang{int(row['azimuth'])}"] = row

    arms = {}
    for f in (70, 50, 30, 20, 10, 5):
        for cone in CONES:
            tail = f"_cone{cone:g}"
            arms[f"matched squared column, cone {cone:g}, q{f}"] = {
                k[: -len(tail)]: v["col_cone_sq"][f"q{f}"]["q"]
                for k, v in res["rows"].items() if k.endswith(tail)}
            arms[f"matched J-factor, cone {cone:g}, q{f}"] = {
                k[: -len(tail)]: v["j_cone"][f"q{f}"]["q"]
                for k, v in res["rows"].items() if k.endswith(tail)}
        arms[f"Muru's own maps, squared column, q{f}"] = {
            k: v[f"sq_q{f}"]["q"] for k, v in res["published"].items()}
    if current:
        for name, field in (("current squared column, 5% contour", "q_dm2"),
                            ("current J-factor, 5% contour", "q_jfactor")):
            vals = {k: v.get(field) for k, v in current.items()
                    if isinstance(v, dict) and v.get(field) is not None}
            if vals:
                arms[name] = vals

    res["correlations"] = {}
    print("\nSpearman rho against dAIC vs VVV (Mask B + PS, 1-10 GeV)")
    for name, vals in arms.items():
        keys = [k for k in vals if k in daic and np.isfinite(vals[k])]
        if len(keys) < 5:
            continue
        rho, p = spearmanr([vals[k] for k in keys], [daic[k] for k in keys])
        res["correlations"][name] = {"rho": float(rho), "p": float(p),
                                     "n": len(keys)}
        print(f"  {name:44s} rho {rho:+.3f}  p {p:.2g}  n {len(keys)}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
