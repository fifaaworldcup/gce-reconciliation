"""J-factor maps for the six HESTIA galaxies, from Muru's particle data.

Abazajian, Kumar & Macias (arXiv:2609.34155) point out that Muru's released
dm2 maps are (int rho ds)^2, not the annihilation integral int rho^2 ds, and
that squaring a non-negative map leaves its isophote shapes unchanged - so a
morphology conclusion drawn from dm2 describes projected mass. With the
particle data this can be computed rather than argued.

For each galaxy both quantities are built from the same particles, the same
observer and the same smoothing, so the only difference between them is the
one under dispute.

The observer sits in the stellar disc plane at R_SUN, at a chosen azimuth from
the in-plane major axis. The disc frame comes from the reduced-inertia tensor
of the stars (computed on the particle files, see disc_frames.json).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project_template as PT          # noqa: E402
import hestia_to_template as H         # noqa: E402
import jfactor_from_particles as JP     # noqa: E402

R_SUN = PT.R_SUN
S_MAX = PT.DEFAULT_LOS_CUT_KPC


def observer_frame(disc_R, azimuth_deg):
    """Rotation taking simulation coords into the frame project_template
    assumes: observer on +x at R_SUN, galactic north along +z.

    disc_R rows are (in-plane major, in-plane minor, disc normal). The
    observer is placed in the disc plane at `azimuth` from the major axis, so
    azimuth 0 looks along it.
    """
    maj, mino, nrm = (np.asarray(disc_R, float)[i] for i in range(3))
    a = np.radians(azimuth_deg)
    obs = np.cos(a) * maj + np.sin(a) * mino        # Sun direction
    yhat = np.cross(nrm, obs)
    yhat /= np.linalg.norm(yhat)
    return np.stack([obs / np.linalg.norm(obs), yhat, nrm])


def maps_for(pos, mass, disc_R, azimuth_deg, half, cells, smooth,
             s_max=S_MAX, fwhm=1.0):
    """(J map, dm2 map) on the Conventions grid, both 1 deg smoothed."""
    R = observer_frame(disc_R, azimuth_deg)
    p = np.asarray(pos, float) @ R.T
    rho, h = JP.density_grid(p, mass, half, cells, smooth_cells=smooth)
    j = JP.jmap(rho, h, half, s_max=s_max, square=True)
    col = JP.jmap(rho, h, half, s_max=s_max, square=False)
    return (PT.normalise(PT.smooth_to(j, fwhm)),
            PT.normalise(PT.smooth_to(col ** 2, fwhm)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reduced", required=True,
                    help="dir with <TAG>_dm_pos.npy and disc_frames.json")
    ap.add_argument("--azimuths", default="0,30,60,90")
    ap.add_argument("--half-kpc", type=float, default=JP.DEFAULT_HALF_KPC)
    ap.add_argument("--cells", type=int, default=JP.DEFAULT_CELLS)
    ap.add_argument("--smooth-cells", type=float, default=1.0)
    ap.add_argument("--level", type=float, default=0.05)
    ap.add_argument("--out")
    a = ap.parse_args()

    frames = json.load(open(os.path.join(a.reduced, "disc_frames.json")))
    tags = sorted(t for t in frames
                  if os.path.exists(os.path.join(a.reduced, t + "_dm_pos.npy")))
    az = [float(x) for x in a.azimuths.split(",")]
    print(f"{len(tags)} galaxies, azimuths {az}, "
          f"grid {a.cells}^3 over +-{a.half_kpc:g} kpc, "
          f"contour at {a.level:.0%} of peak\n")
    print(f"{'galaxy':7} {'azim':>5} {'q(J)':>7} {'q(dm2)':>8} {'dq':>7} "
          f"{'c4(J)':>8} {'c4(dm2)':>9}")
    rows = []
    for t in tags:
        pos = np.load(os.path.join(a.reduced, t + "_dm_pos.npy"))
        mass = np.load(os.path.join(a.reduced, t + "_dm_mass.npy"))
        for z in az:
            j, d = maps_for(pos, mass, frames[t]["frame"], z,
                            a.half_kpc, a.cells, a.smooth_cells)
            qj, c4j = H._shape(*_contour(j, a.level))
            qd, c4d = H._shape(*_contour(d, a.level))
            rows.append({"galaxy": t, "azimuth": z, "q_jfactor": float(qj),
                         "q_dm2": float(qd), "dq": float(qj - qd),
                         "c4_jfactor": float(c4j), "c4_dm2": float(c4d)})
            print(f"{t:7} {z:5.0f} {qj:7.3f} {qd:8.3f} {qj-qd:+7.3f} "
                  f"{c4j:+8.4f} {c4d:+9.4f}")
    dq = np.array([r["dq"] for r in rows])
    print(f"\nacross all {len(rows)} galaxy-angle combinations:")
    print(f"  mean q(J) - q(dm2) = {dq.mean():+.3f}  "
          f"spread {dq.std():.3f}  worst {dq[np.argmax(abs(dq))]:+.3f}")
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        json.dump({"config": vars(a), "rows": rows}, open(a.out, "w"), indent=1)
        print(f"  wrote {a.out}")


def _contour(t, level):
    ls = (20.0 - 0.1 / 2.0) - 0.1 * np.arange(400)
    bs = (-20.0 + 0.1 / 2.0) + 0.1 * np.arange(400)
    LL, BB = np.meshgrid(ls, bs)
    m = t >= level * t.max()
    return LL[m], BB[m]


if __name__ == "__main__":
    main()
