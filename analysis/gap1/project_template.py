"""
Gap 1 - project a 3-D dark matter density field into a rho^2 sky template
on the project's standard grid.

Follows Project Conventions sections 2 (ROI/grid), 5 (template library) and 9
(no hard-coded paths).  PSF convolution is deliberately NOT done here - it is
shared infrastructure and belongs in utils/psf.py so every gap team uses the
identical implementation (Conventions section 5).

Two input paths:

  1. analytic(density_fn)   - any callable rho(r_kpc) -> density.  Used to
                              validate the machinery against a case where the
                              answer is known before any simulation data is
                              involved.
  2. from_grid(cube, ...)   - a regularly-spaced 3-D density cube in galactic
                              cartesian kpc, trilinearly sampled along each
                              sight line.  This is the HESTIA path.

Both return the line-of-sight integral of rho^2, i.e. the J-factor map, in
whatever density units squared x kpc the caller supplied.  Absolute units do
not matter: the template is normalised to unit integral over the ROI, so the
fitted normalisation carries the physical scale.

Run directly to build and self-test a spherical gNFW template:
    python project_template.py --out ../../templates/gnfw_gamma1.2.fits
"""

from __future__ import annotations

import argparse
import datetime as _dt
import subprocess

import numpy as np
from astropy.io import fits

# --------------------------------------------------------------------------
# Project Conventions section 2 - the one true grid.
# --------------------------------------------------------------------------
L_HALF = 20.0          # |l| < 20 deg
B_HALF = 20.0          # |b| < 20 deg
PIX = 0.1              # deg per pixel
NPIX_L = int(round(2 * L_HALF / PIX))   # 400
NPIX_B = int(round(2 * B_HALF / PIX))   # 400

# Conventions section 5
R_SUN = 8.2            # kpc
R_S = 20.0             # kpc, gNFW scale radius

# Muru et al. 2025: 220 pc softening at 8 kpc -> ~1.6 deg.  Anything finer than
# this in a HESTIA template is interpolation, not physics.  Applied as a common
# smoothing floor to EVERY competing template, or the comparison is unfair.
# Upstream truncates every sight line at 15 kpc from the observer
# (`filter!(:d => <=(15*0.677))`).  Our default s_max matches that, so a gNFW
# built here is directly comparable to a HESTIA map.  Raise it only if you also
# know the HESTIA side was not truncated.
DEFAULT_LOS_CUT_KPC = 15.0

HESTIA_SOFTENING_KPC = 0.220
HESTIA_RES_DEG = np.degrees(HESTIA_SOFTENING_KPC / 8.0)   # ~1.575 deg


def sky_grid():
    """Pixel-centre (l, b) in degrees for the standard CAR grid.

    Galactic longitude increases to the left on the sky, so column 0 is
    l = +20 deg.  Returns arrays of shape (NPIX_B, NPIX_L) = (row, col).
    """
    l_edges_first = L_HALF - PIX / 2.0
    l = l_edges_first - PIX * np.arange(NPIX_L)          # +20 -> -20
    b = -B_HALF + PIX / 2.0 + PIX * np.arange(NPIX_B)    # -20 -> +20
    return np.meshgrid(l, b)


def solid_angle():
    """Per-pixel solid angle in steradians for CAR: dl*db*cos(b)."""
    _, bb = sky_grid()
    return np.radians(PIX) ** 2 * np.cos(np.radians(bb))


def _los_samples(l_deg, b_deg, s_max, n_steps):
    """Yield (s, weight) pairs sampling each sight line adaptively.

    A uniform step in s badly undersamples sight lines that pass close to the
    galactic centre, where rho^2 spikes - the first version of this file did
    that and the circular-symmetry self-test caught it (63% scatter at fixed
    angular radius, which would have looked like real morphology).

    Substituting  s = s_peak + r_min*sinh(u)  and stepping uniformly in u puts
    the samples where the integrand actually is: s_peak is the point of closest
    approach to the centre and r_min is that closest distance, both per pixel.
    """
    cos_psi = np.clip(np.cos(np.radians(b_deg)) * np.cos(np.radians(l_deg)), -1, 1)
    sin_psi = np.sqrt(np.maximum(0.0, 1.0 - cos_psi**2))
    s_peak = R_SUN * cos_psi
    r_min = np.maximum(R_SUN * sin_psi, 1e-3)      # floor avoids u -> inf at (0,0)
    u_lo = np.arcsinh((0.0 - s_peak) / r_min)
    u_hi = np.arcsinh((s_max - s_peak) / r_min)
    du = (u_hi - u_lo) / n_steps
    for k in range(n_steps):
        u = u_lo + du * (k + 0.5)
        yield s_peak + r_min * np.sinh(u), r_min * np.cosh(u) * du


def _los_radii(l_deg, b_deg, s):
    """Galactocentric radius r(kpc) at distance s along sight line (l, b)."""
    cos_psi = np.cos(np.radians(b_deg)) * np.cos(np.radians(l_deg))
    return np.sqrt(np.maximum(R_SUN**2 + s**2 - 2.0 * R_SUN * s * cos_psi, 0.0))


def _los_xyz(l_deg, b_deg, s):
    """Galactocentric cartesian (x, y, z) in kpc along the sight line.

    Observer at (R_SUN, 0, 0); x runs from the Sun toward the centre, z is the
    galactic north pole.
    """
    lr, br = np.radians(l_deg), np.radians(b_deg)
    x = R_SUN - s * np.cos(br) * np.cos(lr)
    y = -s * np.cos(br) * np.sin(lr)
    z = s * np.sin(br)
    return x, y, z


def gnfw(gamma=1.2, r_s=R_S):
    """Generalised NFW density profile, arbitrary normalisation."""
    def rho(r):
        r = np.maximum(r, 1e-4)          # avoid the r=0 divergence
        return (r / r_s) ** (-gamma) * (1.0 + r / r_s) ** (-(3.0 - gamma))
    return rho


def analytic(density_fn, s_max=15.0, n_steps=400):
    """Integrate rho^2 along every sight line for a spherical rho(r)."""
    ll, bb = sky_grid()
    out = np.zeros_like(ll)
    for si, w in _los_samples(ll, bb, s_max, n_steps):
        out += density_fn(_los_radii(ll, bb, si)) ** 2 * w
    return out


def from_grid(cube, origin_kpc, spacing_kpc, s_max=15.0, n_steps=400):
    """Integrate rho^2 through a regular 3-D density cube (the HESTIA path).

    cube        : ndarray (nx, ny, nz), density on a regular galactocentric grid
    origin_kpc  : (x0, y0, z0) coordinates of cube[0, 0, 0]
    spacing_kpc : scalar, or (dx, dy, dz)

    Anything outside the cube contributes zero, so make sure the cube covers
    the sight lines you care about before trusting the outer ROI.
    """
    from scipy.ndimage import map_coordinates

    ll, bb = sky_grid()
    spacing = np.broadcast_to(np.asarray(spacing_kpc, float), (3,)).astype(float)
    origin = np.asarray(origin_kpc, float)
    out = np.zeros_like(ll)
    for si, w in _los_samples(ll, bb, s_max, n_steps):
        x, y, z = _los_xyz(ll, bb, si)
        idx = np.stack([(x - origin[0]) / spacing[0],
                        (y - origin[1]) / spacing[1],
                        (z - origin[2]) / spacing[2]])
        rho = map_coordinates(cube, idx, order=1, mode="constant", cval=0.0)
        out += rho ** 2 * w
    return out


def smooth_to(template, fwhm_deg):
    """Gaussian-smooth to a common angular resolution floor.

    Apply the SAME fwhm_deg to every competing template before comparing, or
    a resolution difference will masquerade as a morphology difference.
    """
    from scipy.ndimage import gaussian_filter
    sigma_pix = fwhm_deg / 2.3548 / PIX
    return gaussian_filter(template, sigma_pix, mode="nearest")


def normalise(template):
    """Unit integral over the ROI (Conventions section 5)."""
    total = np.sum(template * solid_angle())
    if total <= 0:
        raise ValueError("template integrates to zero - check inputs")
    return template / total


def _git_hash():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "unknown"


def write_fits(template, path, name, extra=None):
    """Write with WCS keywords so the file is self-describing."""
    h = fits.Header()
    h["CTYPE1"] = "GLON-CAR"
    h["CTYPE2"] = "GLAT-CAR"
    h["CRPIX1"] = (NPIX_L + 1) / 2.0
    h["CRPIX2"] = (NPIX_B + 1) / 2.0
    h["CRVAL1"] = 0.0
    h["CRVAL2"] = 0.0
    h["CDELT1"] = -PIX          # negative: l increases to the left
    h["CDELT2"] = PIX
    h["EQUINOX"] = 2000.0
    h["BUNIT"] = "sr-1"
    h["TEMPNAME"] = (name, "template identifier")
    h["NORMED"] = (True, "unit integral over ROI")
    h["PSFCONV"] = (False, "PSF convolution NOT applied - use utils/psf.py")
    h["GITHASH"] = (_git_hash(), "code version")
    h["DATE"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    for k, v in (extra or {}).items():
        h[k] = v
    fits.PrimaryHDU(data=template.astype(np.float32), header=h).writeto(
        path, overwrite=True)
    return path


# --------------------------------------------------------------------------
# Self-tests.  A spherical input MUST produce a template that depends only on
# angular distance from the centre - if it does not, the geometry is wrong.
# --------------------------------------------------------------------------
def selftest(verbose=True):
    ll, bb = sky_grid()
    ok = True

    if (NPIX_L, NPIX_B) != (400, 400):
        print("FAIL grid is not 400x400"); ok = False

    t = analytic(gnfw(gamma=1.2))

    # 1. peak sits at the centre pixel
    pk = np.unravel_index(np.argmax(t), t.shape)
    if abs(ll[pk]) > PIX or abs(bb[pk]) > PIX:
        print(f"FAIL peak at l={ll[pk]:.2f} b={bb[pk]:.2f}, expected (0,0)"); ok = False

    # 2. spherical in -> circularly symmetric out.
    #    Comparing pixels inside a 1-degree-wide annulus does NOT test this:
    #    the radial falloff across the annulus dominates and you get ~60%
    #    "scatter" from a perfectly symmetric map.  Compare the four cardinal
    #    rays at MATCHED angular radius instead.
    row0 = int(np.argmin(np.abs(bb[:, 0])))
    col0 = int(np.argmin(np.abs(ll[0, :])))
    th = np.arange(1.0, 15.01, 0.25)
    lpos = ll[row0, :]
    bpos = bb[:, col0]
    rays = [
        np.interp(th, lpos[lpos > 0][::-1], t[row0, lpos > 0][::-1]),   # +l
        np.interp(th, -lpos[lpos < 0], t[row0, lpos < 0]),              # -l
        np.interp(th, bpos[bpos > 0], t[bpos > 0, col0]),               # +b
        np.interp(th, -bpos[bpos < 0][::-1], t[bpos < 0, col0][::-1]),  # -b
    ]
    rays = np.array(rays)
    worst = float(np.max((rays.max(0) - rays.min(0)) / rays.mean(0)))
    if worst > 0.01:
        print(f"FAIL circular symmetry broken, {worst:.3%} between cardinal rays")
        ok = False

    theta = np.degrees(np.arccos(np.clip(
        np.cos(np.radians(bb)) * np.cos(np.radians(ll)), -1, 1)))

    # 3. steeper inner slope -> more centrally concentrated
    def frac_inside(tt, deg):
        w = solid_angle()
        return np.sum(tt[theta < deg] * w[theta < deg]) / np.sum(tt * w)
    f12 = frac_inside(analytic(gnfw(1.2)), 5.0)
    f08 = frac_inside(analytic(gnfw(0.8)), 5.0)
    if not f12 > f08:
        print(f"FAIL gamma=1.2 not more concentrated than 0.8 ({f12:.3f} vs {f08:.3f})")
        ok = False

    # 4. normalisation
    n = normalise(t)
    integral = np.sum(n * solid_angle())
    if abs(integral - 1.0) > 1e-9:
        print(f"FAIL normalised integral = {integral}"); ok = False

    if verbose:
        print(f"  grid                {NPIX_L} x {NPIX_B} @ {PIX} deg")
        print(f"  peak pixel          l={ll[pk]:+.2f} b={bb[pk]:+.2f}")
        print(f"  circular symmetry   {worst:.4%} across cardinal rays (want < 1%)")
        print(f"  flux within 5 deg   gamma=1.2 {f12:.3f} | gamma=0.8 {f08:.3f}")
        print(f"  normalised integral {integral:.12f}")
        print(f"  HESTIA res floor    {HESTIA_RES_DEG:.2f} deg = {HESTIA_RES_DEG/PIX:.0f} pixels")
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="gnfw_gamma1.2.fits")
    ap.add_argument("--gamma", type=float, default=1.2)
    ap.add_argument("--smooth", type=float, default=0.0,
                    help="FWHM deg; pass %.2f to match HESTIA resolution"
                         % HESTIA_RES_DEG)
    a = ap.parse_args()
    selftest()
    t = analytic(gnfw(a.gamma))
    if a.smooth > 0:
        t = smooth_to(t, a.smooth)
    t = normalise(t)
    print("wrote", write_fits(t, a.out, f"gnfw_gamma{a.gamma}",
                              {"GAMMA": a.gamma, "RSCALE": R_S, "RSUN": R_SUN,
                               "SMOOTH": a.smooth}))
