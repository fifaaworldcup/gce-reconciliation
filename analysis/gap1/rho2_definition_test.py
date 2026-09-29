# how much does the rho^2 definition actually matter for OUR conclusion?
#
# muru's dm2 maps are (integral rho ds)^2 - the projected column density,
# squared. the annihilation signal is the J-factor, integral rho^2 ds. those
# are different quantities.
#
# for a SPHERICAL halo both are circular, so shape statistics can't tell them
# apart and the question looks harmless. the thing gap 1 actually measures is
# flattening and boxiness, so the real question is:
#
#   does squaring the column instead of integrating the square change the
#   apparent q and c4 of a NON-spherical halo?
#
# if not, muru's maps are still usable for our morphology claim even though the
# label is wrong. if yes, finding 2 is built on the wrong quantity.
#
# tested here on analytic halos where we know the truth: a flattened one and a
# boxy one, both integrated to 15 kpc like muru does.

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hestia_to_template as H
import project_template as P

R_SUN = 8.2
S_MAX = 15.0          # muru truncates the line of sight here
N_STEPS = 600


def density_fn(gamma=1.2, r_s=20.0, q3d=1.0, boxy_n=2.0):
    """gNFW in a deformed radius.

    q3d < 1 flattens along z. boxy_n > 2 makes the isodensity surfaces boxy
    in the x-y plane (superellipsoid), which is the shape HESTIA is claimed
    to have.
    """
    def rho(x, y, z):
        rxy = (np.abs(x) ** boxy_n + np.abs(y) ** boxy_n) ** (1.0 / boxy_n)
        m = np.sqrt(rxy ** 2 + (z / q3d) ** 2)
        m = np.maximum(m, 1e-3)
        return (m / r_s) ** (-gamma) * (1.0 + m / r_s) ** (-(3.0 - gamma))
    return rho


def project(rho, power):
    """Integrate rho**power along every sight line on the conventions grid."""
    ll, bb = P.sky_grid()
    out = np.zeros_like(ll)
    for s, w in P._los_samples(ll, bb, S_MAX, N_STEPS):
        x, y, z = P._los_xyz(ll, bb, s)
        out += rho(x, y, z) ** power * w
    return out


SMOOTH_DEG = 1.0        # match the HESTIA maps we actually measure

def _grid():
    l = (H.L_HALF - H.PIX / 2.0) - H.PIX * np.arange(H.NPIX)
    b = (-H.B_HALF + H.PIX / 2.0) + H.PIX * np.arange(H.NPIX)
    ll, bb = np.meshgrid(l, b)
    return ll, bb, np.hypot(ll, bb)


def _level_for_radius(t, th, target_deg):
    """Find the contour level whose enclosed region has this mean radius.

    A gNFW with r_s = 20 kpc is far more peaked than the 1-deg-smoothed HESTIA
    maps, so fixed fractions of peak (30%, 5%...) all land inside half a degree
    and c4 has too few pixels to fit. Selecting by radius instead puts the
    comparison at the 5-12 deg scales where gap 1 actually measures.
    """
    lo, hi = 1e-9, 1.0
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        m = t >= mid * t.max()
        if not m.any():
            hi = mid; continue
        r = float(th[m].mean())
        if r > target_deg:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def shape(t, targets=(4.0, 6.0, 8.0), guard=15.0):
    from scipy.signal import fftconvolve
    ll, bb, th = _grid()
    r = int(round(SMOOTH_DEG / H.PIX))
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    k = ((x ** 2 + y ** 2) <= r ** 2).astype(float)
    t = fftconvolve(t, k / k.sum(), mode="same")

    rows = []
    for tgt in targets:
        lv = _level_for_radius(t, th, tgt)
        m = t >= lv * t.max()
        if not m.any() or float(th[m].max()) > guard:
            # a contour touching the ROI edge gets a spurious 4-fold signal.
            # relaxing this guard produced c4 = +0.022 for an axisymmetric
            # flattened halo that must be a pure ellipse - the crop artifact,
            # found again independently.
            continue
        q, c4 = H._shape(ll[m], bb[m])
        rows.append((tgt, q, c4, float(th[m].mean())))
    return rows


def compare(label, **kw):
    rho = density_fn(**kw)
    j = project(rho, 2.0)                 # integral rho^2 ds  (correct)
    col = project(rho, 1.0) ** 2          # (integral rho ds)^2  (what muru has)
    print(f"\n{label}")
    print(f"  {'target':>7} {'<r>':>6} | {'q (J)':>7} {'q (col^2)':>10} {'dq':>7} "
          f"| {'c4 (J)':>8} {'c4 (col^2)':>11} {'dc4':>8}")
    sj, sc = shape(j), shape(col)
    worst_q = worst_c = 0.0
    for (lv, qj, cj, r), (_, qc, cc, _) in zip(sj, sc):
        dq, dc = qc - qj, cc - cj
        worst_q = max(worst_q, abs(dq)); worst_c = max(worst_c, abs(dc))
        print(f"  {lv:6.0f}d {r:5.1f}d | {qj:7.3f} {qc:10.3f} {dq:+7.3f} "
              f"| {cj:+8.3f} {cc:+11.3f} {dc:+8.3f}")
    # how different are the radial profiles, which is what a fit would see
    jn, cn = j / j.sum(), col / col.sum()
    m = jn > 0.02 * jn.max()
    print(f"  worst |dq| {worst_q:.3f}   worst |dc4| {worst_c:.3f}   "
          f"median profile difference {np.median(np.abs(jn-cn)[m]/jn[m]):.1%}")
    return worst_q, worst_c


if __name__ == "__main__":
    print("rho^2 definition test: integral(rho^2)ds  vs  (integral rho ds)^2")
    print(f"gNFW gamma=1.2, r_s=20 kpc, observer {R_SUN} kpc, LOS cut {S_MAX} kpc")
    r = []
    r.append(compare("spherical (q3d=1.0, n=2)", q3d=1.0, boxy_n=2.0))
    r.append(compare("flattened (q3d=0.6, n=2)", q3d=0.6, boxy_n=2.0))
    r.append(compare("flattened + boxy (q3d=0.6, n=4)", q3d=0.6, boxy_n=4.0))
    wq = max(x[0] for x in r); wc = max(x[1] for x in r)
    print(f"\nacross all three: worst |dq| = {wq:.3f}, worst |dc4| = {wc:.3f}")
    print("for scale: HESTIA-vs-Coleman c4 separation is 0.006, and the")
    print("c4 noise floor from contour-level scatter is about 0.007")
