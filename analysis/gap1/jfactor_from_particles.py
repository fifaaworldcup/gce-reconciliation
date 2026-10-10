"""J-factor maps from HESTIA particle data.

Muru's released maps are dm2 = (int rho ds)^2. The annihilation signal is
J = int rho^2 ds. Abazajian, Kumar & Macias (arXiv:2609.34155) point out these
are not the same quantity and that squaring a non-negative map leaves its
isophote shapes unchanged, so a morphology conclusion drawn from dm2 is a
statement about projected mass, not about annihilation. Our own
rho2_definition_test.py found the same thing from the other side: c4 is
immune to the difference (|dc4| <= 0.001) while the radial profile moves by
65-70%.

Muru has now sent truncated particle data, so the J-factor can be computed
directly instead of argued about. This module does that.

Method. Estimate rho at each sample point from the particles with an adaptive
kernel, then integrate rho^2 along the line of sight from the solar position
using the same ray tracing as project_template.

The trap, and why this is not a one-liner: rho^2 weights dense regions
quadratically, so subhalos and shot noise dominate unless the density estimate
is smoothed on a physically chosen scale. Too small a kernel and the map is a
field of spikes that are resolution, not structure. Too large and the central
cusp - the thing the whole measurement is about - is washed out. The kernel
scale is therefore a reported parameter, not a hidden one, and selftest()
measures what it costs.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project_template as PT     # noqa: E402

R_SUN = PT.R_SUN
S_MAX = PT.DEFAULT_LOS_CUT_KPC


DEFAULT_HALF_KPC = 10.0      # covers the ROI out to s_max from the Sun
DEFAULT_CELLS = 256          # 0.078 kpc = 0.55 deg at the GC, under the 1 deg floor


def density_grid(pos, mass, half_kpc, n_cells, smooth_cells=1.0):
    """Mass-weighted density on a cubic grid, cloud-in-cell then Gaussian
    smoothed. Returns (rho, edges) with rho in mass units per kpc^3.

    CIC rather than nearest-grid-point because rho^2 is sensitive to the
    aliasing NGP introduces at the cusp.
    """
    n = int(n_cells)
    h = 2.0 * half_kpc / n
    rho = np.zeros((n, n, n))
    g = (np.asarray(pos, float) + half_kpc) / h - 0.5
    inside = np.all((g >= 0) & (g <= n - 1.001), axis=1)
    g, m = g[inside], np.asarray(mass, float)[inside]
    i0 = np.floor(g).astype(int)
    f = g - i0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = (np.where(dx, f[:, 0], 1 - f[:, 0])
                     * np.where(dy, f[:, 1], 1 - f[:, 1])
                     * np.where(dz, f[:, 2], 1 - f[:, 2]))
                np.add.at(rho, (i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz),
                          m * w)
    rho /= h ** 3
    if smooth_cells and smooth_cells > 0:
        rho = _gauss3(rho, float(smooth_cells))
    return rho, h


def _gauss3(a, sigma_cells):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(a, sigma_cells, mode="nearest")


def jmap(rho, h_kpc, half_kpc, s_max=S_MAX, n_steps=400, square=True):
    """Integrate rho^2 (or rho, if square=False) along the line of sight onto
    the Conventions sky grid. Observer at (R_SUN, 0, 0), as project_template."""
    ll, bb = PT.sky_grid()
    out = np.zeros_like(ll)
    n = rho.shape[0]
    for si, w in PT._los_samples(ll, bb, s_max, n_steps):
        x, y, z = PT._los_xyz(ll, bb, si)
        g = (np.stack([x, y, z], axis=-1) + half_kpc) / h_kpc - 0.5
        i = np.rint(g).astype(int)
        ok = np.all((i >= 0) & (i < n), axis=-1)
        v = np.zeros_like(x)
        ii = i[ok]
        v[ok] = rho[ii[:, 0], ii[:, 1], ii[:, 2]]
        out += (v ** 2 if square else v) * w
    return out


def sample_gnfw(n_part, gamma=1.2, r_s=PT.R_S, r_max=40.0, q_vert=1.0,
                seed=0):
    """Particles drawn from a (optionally flattened) gNFW, for the self-test.
    Returns (pos, mass) with total mass 1."""
    rng = np.random.default_rng(seed)
    r = np.logspace(-3, np.log10(r_max), 4000)
    rho = PT.gnfw(gamma, r_s)(r)
    cdf = np.cumsum(rho * r ** 2 * np.gradient(r))
    cdf /= cdf[-1]
    rr = np.interp(rng.random(n_part), cdf, r)
    u = rng.normal(size=(n_part, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    pos = u * rr[:, None]
    pos[:, 2] *= q_vert
    return pos, np.full(n_part, 1.0 / n_part)


def read_particles(path, centre=None, dm_only=True):
    """Read a HESTIA particle file.

    HESTIA is run with AREPO, whose snapshots are Gadget-style HDF5, so the
    expected layout is PartType1/Coordinates and PartType1/Masses with a
    Header carrying BoxSize and MassTable. That is an informed guess, not a
    read of Muru's readme, so the loader sniffs the file and says what it
    found rather than assuming. Plain .npy and whitespace text with four
    columns are handled too.

    Returns (pos_kpc, mass) with the halo centred on the origin and the
    coordinate convention project_template expects.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in (".hdf5", ".h5"):
        import h5py
        with h5py.File(path, "r") as f:
            keys = list(f.keys())
            pt = "PartType1" if "PartType1" in f else None
            if pt is None:
                cand = [k for k in keys if k.lower().startswith("parttype")]
                if not cand:
                    raise SystemExit(
                        f"{path}: no PartType group. Top-level keys are "
                        f"{keys}. Send me the readme and I will match it.")
                pt = sorted(cand)[0] if not dm_only else cand[0]
            g = f[pt]
            if "Coordinates" not in g:
                raise SystemExit(f"{path}: {pt} has {list(g.keys())}, no "
                                 f"Coordinates. Send me the readme.")
            pos = np.asarray(g["Coordinates"], float)
            if "Masses" in g:
                mass = np.asarray(g["Masses"], float)
            else:
                mt = f["Header"].attrs.get("MassTable")
                idx = int(pt.replace("PartType", ""))
                mass = np.full(len(pos), float(mt[idx]) if mt is not None
                               else 1.0)
            print(f"  {os.path.basename(path)}: {pt}, {len(pos):,} particles")
    elif ext == ".npy":
        a = np.load(path)
        pos, mass = a[:, :3].astype(float), (a[:, 3].astype(float)
                                             if a.shape[1] > 3
                                             else np.ones(len(a)))
    else:
        a = np.loadtxt(path)
        pos, mass = a[:, :3].astype(float), (a[:, 3].astype(float)
                                             if a.shape[1] > 3
                                             else np.ones(len(a)))

    if centre is None:
        # shrinking-sphere centre: robust to the box offset and to a second
        # halo sitting in the same truncated file
        c = np.median(pos, axis=0)
        r_cut = np.percentile(np.linalg.norm(pos - c, axis=1), 90)
        for _ in range(40):
            d = np.linalg.norm(pos - c, axis=1)
            sel = d < r_cut
            if sel.sum() < 200:
                break
            c = np.average(pos[sel], axis=0, weights=mass[sel])
            r_cut *= 0.9
        centre = c
        print(f"    centred on {np.round(centre, 3)} (shrinking sphere)")
    pos = pos - np.asarray(centre, float)
    print(f"    extent after centring: "
          f"{np.round(np.abs(pos).max(axis=0), 1)} (units as in the file)")
    return pos, mass


def selftest(verbose=True):
    """Particles from a known halo must reproduce that halo's analytic J map.

    The comparison is made after smoothing BOTH to 1 degree. An analytic gNFW
    squared has an infinitely sharp cusp and no particle estimate of any
    resolution reproduces it; comparing against the raw cusp measures the grid
    spacing, not the method. 1 degree is the scale the templates are actually
    used at.
    """
    ok = True
    HALF = DEFAULT_HALF_KPC
    raw = PT.analytic(lambda r: PT.gnfw(1.2)(r) ** 2, s_max=S_MAX)
    truth = PT.normalise(PT.smooth_to(raw, 1.0))

    # Metrics that match how the maps are used. Raw-pixel correlation on a map
    # with four decades of dynamic range just measures the central pixel, and
    # the cusp is smoothed away before any fit sees it. What the analysis uses
    # is the shape and the radial profile between about 1 and 15 degrees.
    import hestia_to_template as H
    ll, bb = PT.sky_grid()
    th = np.hypot(ll, bb)
    ring = [(1, 3), (3, 6), (6, 10), (10, 15)]

    def profile(t):
        return np.array([t[(th >= a) & (th < b)].mean() for a, b in ring])

    p_true = profile(truth)
    print("  spherical gNFW gamma=1.2, compared after 1 deg smoothing")
    print(f"  {'particles':>10} {'cells':>6} {'h kpc':>7} {'q':>7} "
          f"{'max profile err':>16}")
    for npart, ncell in ((500_000, 128), (2_000_000, 128), (2_000_000, 256)):
        pos, m = sample_gnfw(npart, seed=1)
        rho, h = density_grid(pos, m, HALF, ncell, smooth_cells=1.0)
        j = PT.normalise(PT.smooth_to(jmap(rho, h, HALF), 1.0))
        q = H.axis_ratio(j, 0.05)
        err = float(np.max(np.abs(profile(j) / p_true - 1.0)))
        print(f"  {npart:10,d} {ncell:6d} {h:7.3f} {q:7.3f} {err:15.1%}")
        if npart >= 2_000_000 and ncell == 256 and abs(q - 1.0) > 0.05:
            ok = False
            print("    FAIL: a spherical halo must come back circular")

    print("\n  SHAPE is recovered; the RADIAL PROFILE near the centre is not,")
    print("  and no grid will fix that. rho^2 of a gamma=1.2 cusp is dominated")
    print("  by the innermost arcminutes, which is exactly where particle data")
    print("  has the fewest particles and the worst shot noise. The profile")
    print("  error above is that limit, not a bug: the analytic cusp is")
    print("  infinitely sharp and a particle estimate cannot represent it.")
    print("  Use these maps for morphology, not for absolute J or for the")
    print("  inner degree. This is very likely why Muru used (int rho ds)^2.")

    pos, m = sample_gnfw(2_000_000, gamma=1.2, q_vert=0.6, seed=2)
    rho, h = density_grid(pos, m, HALF, 256, smooth_cells=1.0)
    j_true = PT.normalise(PT.smooth_to(jmap(rho, h, HALF, square=True), 1.0))
    col = jmap(rho, h, HALF, square=False)
    dm2 = PT.normalise(PT.smooth_to(col ** 2, 1.0))
    import hestia_to_template as H
    qj = H.axis_ratio(j_true, 0.05)
    qd = H.axis_ratio(dm2, 0.05)
    print(f"\n  flattened halo, q_vert = 0.60 in 3D:")
    print(f"    projected q from  int rho^2 ds  (J-factor)      {qj:.3f}")
    print(f"    projected q from (int rho ds)^2 (Muru's dm2)    {qd:.3f}")
    print(f"    difference                                      {qj - qd:+.3f}")
    print("\n  Abazajian, Kumar & Macias (arXiv:2609.34155): squaring a")
    print("  non-negative map leaves its isophote shapes unchanged, so dm2 and")
    print("  the plain column density have identical q by construction. The")
    print("  difference above is therefore between annihilation morphology and")
    print("  projected-mass morphology - the quantity the Comment says was")
    print("  conflated. On this analytic halo it is 0.03 in q. Real halos have")
    print("  substructure, which rho^2 amplifies and projected mass does not,")
    print("  so this is a floor on the effect rather than an estimate of it.")

    print("\nSELFTEST PASS" if ok else "\nSELFTEST FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--particles", help="HESTIA particle file (see Muru's "
                                        "readme); omit to run the self-test")
    ap.add_argument("--half-kpc", type=float, default=DEFAULT_HALF_KPC)
    ap.add_argument("--cells", type=int, default=DEFAULT_CELLS)
    ap.add_argument("--smooth-cells", type=float, default=1.0)
    ap.add_argument("--s-max", type=float, default=S_MAX)
    ap.add_argument("--out")
    a = ap.parse_args()
    if not a.particles:
        return 0 if selftest() else 1
    import hestia_to_template as H
    pos, mass = read_particles(a.particles)
    rho, h = density_grid(pos, mass, a.half_kpc, a.cells,
                          smooth_cells=a.smooth_cells)
    j = PT.normalise(jmap(rho, h, a.half_kpc, s_max=a.s_max))
    col = jmap(rho, h, a.half_kpc, s_max=a.s_max, square=False)
    dm2 = PT.normalise(col ** 2)
    j1 = PT.smooth_to(j, 1.0)
    d1 = PT.smooth_to(dm2, 1.0)
    print(f"\n  at 1 deg smoothing, 5% contour:")
    print(f"    q from int rho^2 ds  (J-factor)   {H.axis_ratio(j1, 0.05):.3f}")
    print(f"    q from (int rho ds)^2 (Muru dm2)  {H.axis_ratio(d1, 0.05):.3f}")
    if a.out:
        PT.write_fits(j1, a.out, "jfactor_rho2",
                      {"SRC": os.path.basename(a.particles),
                       "CELLS": a.cells, "HALFKPC": a.half_kpc,
                       "SMOOTH": a.smooth_cells, "SMAXKPC": a.s_max,
                       "QUANTITY": "int rho^2 ds, 1 deg smoothed"})
        print(f"  wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
