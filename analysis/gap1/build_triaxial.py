"""
Triaxial gNFW J-factor templates, for the Gap 4 cross-check.

Conventions sec 5 gives Hu et al. (2026) as q1 = 0.6, q2 = 0.8, theta = 25 deg,
"standard and flipped", but not which axis each q belongs to. That is not a
detail here: Gap 1's fits select a VERTICAL flattening of q ~ 0.6, so reading
q1 as vertical predicts the triaxial halo does well and reading it as in-plane
predicts it doesn't. So build both readings rather than guess, and let Gap 4
say which one is Hu's.

Density:  rho(m) = gNFW(m),  m^2 = x'^2 + (y'/q_plane)^2 + (z/q_vert)^2
with (x', y') the galactocentric plane rotated by theta about the pole; x
points from the Sun to the centre. theta > 0 = "standard", < 0 = "flipped"
(sign convention to confirm with Gap 4). Template = integral of rho^2 ds,
same sight-line sampler and 15 kpc cut as every other template here.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project_template as PT          # noqa: E402
import hestia_to_template as H         # noqa: E402

GAMMA = 1.2


def triaxial(gamma=GAMMA, q_plane=1.0, q_vert=1.0, theta_deg=0.0,
             s_max=PT.DEFAULT_LOS_CUT_KPC, n_steps=400):
    ll, bb = PT.sky_grid()
    rho = PT.gnfw(gamma)
    th = np.radians(theta_deg)
    c, s_ = np.cos(th), np.sin(th)
    out = np.zeros_like(ll)
    for si, w in PT._los_samples(ll, bb, s_max, n_steps):
        x, y, z = PT._los_xyz(ll, bb, si)
        xp = c * x + s_ * y
        yp = -s_ * x + c * y
        m = np.sqrt(xp ** 2 + (yp / q_plane) ** 2 + (z / q_vert) ** 2)
        out += rho(m) ** 2 * w
    return PT.normalise(out)


VARIANTS = {
    # name: (q_plane, q_vert, theta)
    "vert0.6_plane0.8_std":  (0.8, 0.6, +25.0),
    "vert0.6_plane0.8_flip": (0.8, 0.6, -25.0),
    "vert0.8_plane0.6_std":  (0.6, 0.8, +25.0),
    "vert0.8_plane0.6_flip": (0.6, 0.8, -25.0),
}


def selftest():
    ok = True
    sph = PT.normalise(PT.analytic(PT.gnfw(GAMMA)))
    tri1 = triaxial(q_plane=1.0, q_vert=1.0, theta_deg=33.0)
    d = np.max(np.abs(tri1 - sph)) / sph.max()
    print(f"  q=1 triaxial vs spherical builder: max diff {d:.2e} of peak")
    if d > 1e-6:
        print("FAIL q=1 should reproduce the spherical template"); ok = False
    # vertical squash must show up as a projected axis ratio < 1, and rotating
    # in the plane must not change a vertically-only squashed halo
    v = triaxial(q_plane=1.0, q_vert=0.6, theta_deg=0.0)
    v_rot = triaxial(q_plane=1.0, q_vert=0.6, theta_deg=25.0)
    dr = np.max(np.abs(v - v_rot)) / v.max()
    print(f"  vertical-only squash, rotated vs not: {dr:.2e} (want ~0)")
    if dr > 1e-6:
        print("FAIL in-plane rotation changed an axisymmetric halo"); ok = False
    ll, bb = PT.sky_grid()
    m = v >= 0.05 * v.max()
    q, _ = H._shape(ll[m], bb[m])
    print(f"  q_vert=0.6 projects to axis ratio {q:.3f} at the 5% contour")
    if not q < 0.9:
        print("FAIL vertical squash did not flatten the projection"); ok = False
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    from astropy.io import fits  # noqa: F401
    if not selftest():
        sys.exit(1)
    outdir = os.path.join(HERE, "..", "..", "templates", "triaxial")
    os.makedirs(outdir, exist_ok=True)
    ll, bb = PT.sky_grid()
    for name, (qp, qv, th) in VARIANTS.items():
        t = triaxial(q_plane=qp, q_vert=qv, theta_deg=th)
        path = os.path.join(outdir, f"triaxial_{name}.fits")
        PT.write_fits(t, path, f"triaxial_{name}",
                      {"GAMMA": GAMMA, "QPLANE": qp, "QVERT": qv,
                       "THETA": th, "RSCALE": PT.R_S, "RSUN": PT.R_SUN})
        m = t >= 0.05 * t.max()
        q, c4 = H._shape(ll[m], bb[m])
        print(f"  wrote {name:24s} projected q={q:.3f} c4={c4:+.4f}")
