"""HESTIA maps rebuilt with gce_in_hestia's own conventions.

Everything here is read off Muru's Julia source, not inferred:

  src/GalacticCoordinates.jl
    d_galcenter = 8.0*0.677 kpc/h, i.e. the observer sits at 8.0 physical kpc
    R = [e_centre, e_centre x e_normal, e_normal]; l = atan2(y,x), b = asin(z/d)
    p_sun = centre + d*(cos a * v_disk + sin a * (v_normal x v_disk))
    v_normal and v_disk come from the AHF profile, not from a fit to the
    particles, so the disc frame is AHF's inertia tensor over all particles.

  src/calc_methods_densityprojections.jl
    filter d <= 15*0.677 kpc/h   -- 15 physical kpc from the SUN, not the centre
    filter ptype == :dm
    grid -20:0.1:20 in l and in b
    calc_density_on_grid(..., 3.0)  -- cone radius in degrees

  src/BulgeDensityProjection.jl
    map value = sum of the masses whose (l,b) fall within the cone radius of
    the sight line. A KDTree on (l,b) in degrees under the Euclidean metric,
    so the cone is a disc in the l-b plane and its solid angle is the same at
    every distance. The far end of the sight line therefore carries a d^2
    weight: the map is int rho s^2 ds, not int rho ds.

  notebooks/azimuthal_analysis_erebos.jl
    dm2 = map .^ 2
    q   = min/max singular value of the mean-centred (l,b) of every pixel at
          or above frac * peak, unweighted, frac in (0.7, 0.5, 0.3).

That last line settles the squaring question by itself. The set
{m^2 >= f*max(m^2)} is the set {m >= sqrt(f)*max(m)}, so under his own
estimator q(m^2, f) = q(m, sqrt(f)) exactly: squaring a column map does not
make a new shape, it steps along the level sequence. selftest() checks it.

Positions are expected in physical kpc about the AHF centre, masses in
10^10 Msol -- the output of reduce.py. The 220 pc softening is the floor on
the kNN density.
"""

import json
import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree

H = 0.677
R_SUN = 8.0                 # kpc, 8.0*0.677 kpc/h in the Julia
D_MAX = 15.0                # kpc from the observer, 15*0.677 kpc/h
HALF = 20.0                 # grid half-width, deg
STEP = 0.1                  # grid step, deg
CONE = 3.0                  # cone radius, deg
AHF_R_HINV = 10.0           # profile row to read the disc frame off, kpc/h
SOFTENING = 0.220           # kpc
D_MIN = 1.0                 # kpc, floor on the 1/d^2 weighting
KNN = 32

# AHF profile columns, 1-based as in the AHF manual (pp. 170-171)
_COL = {"r": 1, "npart": 2, "Lx": 9, "Ly": 10, "Lz": 11, "b": 12, "c": 13,
        "Eax": 14, "Eay": 15, "Eaz": 16, "Ebx": 17, "Eby": 18, "Ebz": 19,
        "Ecx": 20, "Ecy": 21, "Ecz": 22}


def read_ahf_frame(path, r_hinv=AHF_R_HINV):
    """Disc vectors from the AHF profile row nearest r_hinv kpc/h.

    Returns the minor-axis eigenvector Ec (the disc normal), the major axis
    Ea, the angular momentum direction, and the row actually used. AHF writes
    the radius negative for bins it does not consider converged; those rows
    are skipped.
    """
    rows = np.loadtxt(path, comments="#")
    r = rows[:, _COL["r"] - 1]
    ok = r > 0
    if not ok.any():
        raise ValueError(f"{path}: no converged profile rows")
    idx = np.where(ok)[0][np.argmin(np.abs(r[ok] - r_hinv))]
    row = rows[idx]
    get = lambda *k: np.array([row[_COL[x] - 1] for x in k], float)
    out = {
        "r_hinv": float(row[_COL["r"] - 1]),
        "r_kpc": float(row[_COL["r"] - 1] / H),
        "npart": int(row[_COL["npart"] - 1]),
        "b_over_a": float(row[_COL["b"] - 1]),
        "c_over_a": float(row[_COL["c"] - 1]),
        "Ea": _unit(get("Eax", "Eay", "Eaz")),
        "Eb": _unit(get("Ebx", "Eby", "Ebz")),
        "Ec": _unit(get("Ecx", "Ecy", "Ecz")),
        "L": _unit(get("Lx", "Ly", "Lz")),
    }
    out["normal_vs_L_deg"] = float(
        np.degrees(np.arccos(min(1.0, abs(float(out["Ec"] @ out["L"]))))))
    return out


def _unit(v):
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("zero vector")
    return v / n


def sun_position(normal, disk, azimuth_deg, r_sun=R_SUN):
    """generate_sun_position with the galactic centre at the origin."""
    v1 = _unit(np.asarray(normal, float))
    v2 = _unit(np.asarray(disk, float))
    a = np.radians(azimuth_deg)
    return r_sun * _unit(np.cos(a) * v2 + np.sin(a) * np.cross(v1, v2))


def rotation(normal, observer):
    """calc_rotmat_for_lb_conversion with p_galcenter at the origin.

    e_x points from the observer to the centre, e_z is the disc normal,
    e_y = e_x x e_z. That is the opposite handedness to standard galactic
    coordinates, so l runs backwards; it mirrors the map in l and leaves any
    axis ratio alone.
    """
    e_n = _unit(np.asarray(normal, float))
    e_c = _unit(-np.asarray(observer, float))
    e_y = _unit(np.cross(e_c, e_n))
    return np.stack([e_c, e_y, e_n])


def galactic(pos, observer, R):
    """l, b in degrees and d in kpc for positions about the galactic centre."""
    p = (np.asarray(pos, np.float64) - observer) @ R.T
    d = np.linalg.norm(p, axis=1)
    l = np.degrees(np.arctan2(p[:, 1], p[:, 0]))
    b = np.degrees(np.arcsin(np.clip(p[:, 2] / np.maximum(d, 1e-12), -1, 1)))
    return l, b, d


def knn_density(pos, mass, k=KNN, softening=SOFTENING):
    """Local density per particle from the k-th neighbour distance.

    rho_i = (mass inside the k-neighbour ball) / (volume of that ball), with
    the radius floored at the simulation softening so the cusp cannot invent
    a density the simulation never resolved.
    """
    tree = cKDTree(np.asarray(pos, np.float64))
    dist, idx = tree.query(pos, k=k + 1, workers=-1)
    r_k = np.maximum(dist[:, -1], softening)
    m_in = np.asarray(mass, np.float64)[idx].sum(axis=1)
    return m_in / (4.0 / 3.0 * np.pi * r_k ** 3)


def cone_map(l, b, weights, half=HALF, step=STEP, cone=CONE):
    """Sum of weights within `cone` degrees of each grid point.

    Muru walks the KDTree once per grid point; the same thing is a top-hat
    convolution of the binned particles, which is exact up to binning the
    particle positions to the 0.1 deg cell they already land in.

    Returns (map, l_axis, b_axis) with the map indexed [il, ib].
    """
    pad = cone + 2 * step
    edge = half + pad
    n = int(round(2 * edge / step)) + 1
    ctr = -edge + step * np.arange(n)
    bins = np.concatenate([ctr - step / 2, [ctr[-1] + step / 2]])
    counts, _, _ = np.histogram2d(l, b, bins=[bins, bins], weights=weights)

    rpix = int(round(cone / step))
    g = step * np.arange(-rpix, rpix + 1)
    kern = (g[:, None] ** 2 + g[None, :] ** 2) <= cone ** 2 + 1e-9
    full = fftconvolve(counts, kern.astype(np.float64), mode="same")

    k0 = int(round(pad / step))
    k1 = k0 + int(round(2 * half / step)) + 1
    axis = ctr[k0:k1]
    return full[k0:k1, k0:k1], axis, axis


def axis_ratio_svd(m, frac, l_axis, b_axis):
    """calc_axis_ratio_svd: min/max singular value of the super-level pixels.

    Also returns which grid axis the long one is, since min/max is blind to
    orientation and only equals (vertical / in-plane) while the structure is
    flattened towards the plane.
    """
    lvl = float(np.nanmax(m)) * frac
    il, ib = np.nonzero(m >= lvl)
    if il.size < 3:
        return {"q": float("nan"), "npix": int(il.size)}
    xy = np.stack([l_axis[il], b_axis[ib]], axis=1)
    xy = xy - xy.mean(axis=0)
    s = np.linalg.svd(xy, compute_uv=False)
    u = np.linalg.svd(xy, full_matrices=False)[2]
    major = u[0]
    return {
        "q": float(s.min() / s.max()),
        "npix": int(il.size),
        "major_axis_pa_deg": float(np.degrees(np.arctan2(major[1], major[0]))),
        "major_is_longitude": bool(abs(major[0]) >= abs(major[1])),
        "sigma_l": float(xy[:, 0].std()),
        "sigma_b": float(xy[:, 1].std()),
    }


def maps_for(pos, mass, rho, frame, azimuth_deg, cone=CONE,
             d_max=D_MAX, half=HALF, step=STEP):
    """The four weightings of one sight-line map for one viewing angle."""
    obs = sun_position(frame["normal"], frame["disk"], azimuth_deg)
    R = rotation(frame["normal"], obs)
    l, b, d = galactic(pos, obs, R)

    keep = (d <= d_max) & (np.abs(l) <= half + cone + 1.0) \
        & (np.abs(b) <= half + cone + 1.0)
    l, b, d = l[keep], b[keep], d[keep]
    m, r = np.asarray(mass, np.float64)[keep], np.asarray(rho, np.float64)[keep]
    # The 1/d^2 arms are a diagnostic, not Muru's convention: they turn his
    # cone sum into a real surface density. Without a floor a particle a few
    # hundred pc from the observer outweighs the whole Galactic Centre and the
    # map's peak lands in a random direction, so the level sets become
    # meaningless. D_MIN is that floor.
    inv2 = 1.0 / np.maximum(d, D_MIN) ** 2

    out = {}
    for name, w in (("col_cone", m), ("j_cone", m * r),
                    ("col_los", m * inv2), ("j_los", m * r * inv2)):
        out[name], la, ba = cone_map(l, b, w, half=half, step=step, cone=cone)
    out["_axes"] = (la, ba)
    out["_nkept"] = int(l.size)
    return out


def summarise(maps, fracs=(0.7, 0.5, 0.3, 0.2, 0.1, 0.05)):
    la, ba = maps["_axes"]
    rows = {}
    for name in ("col_cone", "j_cone", "col_los", "j_los"):
        m = maps[name]
        rows[name] = {f"q{int(f * 100)}": axis_ratio_svd(m, f, la, ba)
                      for f in fracs}
        rows[name + "_sq"] = {f"q{int(f * 100)}": axis_ratio_svd(m ** 2, f, la, ba)
                              for f in fracs}
    return rows


def selftest():
    """Three checks that do not need the simulation.

    1. Squaring a map cannot change the axis ratio sequence, only shift it:
       q(m^2, f) must equal q(m, sqrt(f)) to the last bit.
    2. The cone convolution must reproduce a direct KDTree neighbour sum at
       the production grid, where the only error is binning particles to the
       0.1 deg cell they already sit in.
    3. That error must not move the axis ratio.
    """
    rng = np.random.default_rng(0)
    l = rng.normal(0, 4.0, 400_000)
    b = rng.normal(0, 2.0, 400_000)
    w = np.ones_like(l)
    m, la, ba = cone_map(l, b, w)

    ok1 = True
    for f in (0.7, 0.5, 0.3):
        a = axis_ratio_svd(m ** 2, f, la, ba)["q"]
        c = axis_ratio_svd(m, np.sqrt(f), la, ba)["q"]
        ok1 &= a == c
        print(f"  frac {f:.2f}: q(m^2) {a:.6f}  q(m, sqrt f) {c:.6f}"
              f"  {'same' if a == c else 'DIFFER'}")

    tree = cKDTree(np.stack([l, b], axis=1))
    sub = rng.choice(la.size, 40, replace=False)
    pts = np.stack([la[sub], ba[rng.choice(ba.size, 40, replace=False)]], axis=1)
    direct = np.array([len(tree.query_ball_point(p, CONE)) for p in pts], float)
    conv = np.array([m[np.argmin(abs(la - p[0])), np.argmin(abs(ba - p[1]))]
                     for p in pts])
    rel = np.abs(conv - direct) / np.maximum(direct, 1)
    tot = np.abs(conv - direct).sum() / max(direct.sum(), 1)
    busy = direct > 100
    print(f"  cone vs direct, 40 sight lines: median {np.median(rel):.5f}"
          f"  worst well-populated {rel[busy].max() if busy.any() else 0:.5f}"
          f"  summed {tot:.5f}")
    ok2 = tot < 0.005 and (not busy.any() or rel[busy].max() < 0.01)

    # the same map from the tree, to see whether the binning moves q at all
    gl, gb = np.meshgrid(la, ba, indexing="ij")
    flat = np.stack([gl.ravel(), gb.ravel()], axis=1)
    hits = tree.query_ball_point(flat, CONE, workers=-1, return_length=True)
    exact = hits.reshape(gl.shape).astype(float)
    ok3 = True
    for f in (0.7, 0.5, 0.3):
        a = axis_ratio_svd(m, f, la, ba)["q"]
        e = axis_ratio_svd(exact, f, la, ba)["q"]
        ok3 &= abs(a - e) < 2e-3
        print(f"  frac {f:.2f}: q binned {a:.5f}  q exact {e:.5f}"
              f"  diff {a - e:+.5f}")

    print(f"  {'PASS' if ok1 and ok2 and ok3 else 'FAIL'}")
    return ok1 and ok2 and ok3


if __name__ == "__main__":
    selftest()
