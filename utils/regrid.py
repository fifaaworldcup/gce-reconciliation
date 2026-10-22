# WCS-based regridder onto the project grid
#
# every template we've been handed sits on a slightly different CAR grid:
#   Muru HESTIA  800x800 @0.1 deg, CRPIX 400.5, CDELT1 = +0.1  (needs l flip)
#   Coleman/F98  200x200 @0.2 deg, CRPIX 100.5, CDELT1 = -0.2
#   bubbles       80x80  @0.5 deg, CRPIX 40    , CDELT1 = -0.5  (half-pixel offset)
#
# hand-coding the index maths per file is how you get a silent 0.05 deg shift.
# this reads the WCS and lets astropy do it, so a new file with any CRPIX,
# CDELT sign or size lands correctly without us re-deriving anything.

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

PIX = 0.1
NPIX = 400
HALF = 20.0


def conventions_wcs():
    """The project grid as a WCS: 400x400, 0.1 deg, CAR, CDELT1 negative."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
    w.wcs.crpix = [(NPIX + 1) / 2.0, (NPIX + 1) / 2.0]
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cdelt = [-PIX, PIX]
    return w


def load(path, hdu=0, order=1):
    """Read any CAR template and return it on the conventions grid as [b, l].

    Bilinear by default. Returns (map, coverage_fraction).

    Edge policy matters and is easy to get wrong. A source pixel CENTRE sits
    half a pixel inside that pixel's footprint, so sky just beyond the outermost
    centre is still inside the map. Zeroing it puts a ring of zeros around the
    template at the ROI rim. So we accept anything within half a source pixel of
    the outermost centre and fill it by nearest-neighbour; only sky genuinely
    off the map comes back 0. Caught by cross-checking against the hand-coded
    Coleman loader, which disagreed by 100% at l = +19.95 and nowhere else.
    """
    with fits.open(path) as hd:
        src = np.asarray(hd[hdu].data, float)
        swcs = WCS(hd[hdu].header)
    return regrid_array(src, swcs, order=order)


def regrid_array(src, src_wcs, order=1):
    """Same as load() but for an array already in memory.

    Split out so the diffuse model cube can be regridded plane by plane without
    writing 30 intermediate files. load() goes through here, so both paths share
    the edge policy and there is only one place to get it wrong.
    """
    from scipy.ndimage import map_coordinates

    src = np.asarray(src, float)
    swcs = src_wcs.celestial if hasattr(src_wcs, "celestial") else src_wcs

    l = (HALF - PIX / 2.0) - PIX * np.arange(NPIX)
    b = (-HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    ll, bb = np.meshgrid(l, b)

    # our sky coords -> their pixel coords, whatever their WCS happens to be
    x, y = swcs.wcs_world2pix(ll, bb, 0)
    tol = 0.5                                   # half a source pixel
    inside = ((x >= -tol) & (x <= src.shape[1] - 1 + tol) &
              (y >= -tol) & (y <= src.shape[0] - 1 + tol))
    out = map_coordinates(src, np.stack([y, x]), order=order, mode="nearest")
    out[~inside] = 0.0
    return out, float(inside.mean())


def solid_angle():
    b = (-HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    return np.radians(PIX) ** 2 * np.cos(np.radians(b))[:, None] * np.ones((1, NPIX))


def normalise(t):
    tot = float(np.sum(t * solid_angle()))
    if not np.isfinite(tot) or tot <= 0:
        raise ValueError("template integrates to zero")
    return t / tot


def selftest():
    """Cross-check against the two hand-coded loaders.

    If the WCS path and the hand-coded path agree, both are probably right.
    If they disagree, at least one has an offset and we want to know now.
    """
    import os, sys
    ok = True
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "..", "analysis", "gap1"))
    import hestia_to_template as H
    import bulge_to_template as B

    checks = []

    hest_src = os.path.expanduser(
        "~/muru_data/jf/density_projection_09_18_8192"
        "_halo127000000000003_angle0.0_dm2_1deg.fits")
    if os.path.exists(hest_src):
        hand = H.normalise(H.to_conventions_grid(H.read_projection(hest_src)[0]))
        wcsv, cov = load(hest_src)
        checks.append(("HESTIA 800x800 CDELT1>0", hand, normalise(wcsv), cov))

    col = os.path.join(here, "..", "templates", "bulge",
                       B.FILES["coleman"])
    if os.path.exists(col):
        hand = H.normalise(B.to_conventions_grid(B.read(col)))
        wcsv, cov = load(col)
        checks.append(("Coleman 200x200 CDELT1<0", hand, normalise(wcsv), cov))

    for name, a, b, cov in checks:
        m = a > 0.02 * a.max()
        worst = float(np.max(np.abs(a - b)[m] / a[m])) if m.any() else float("nan")
        flag = "" if worst < 0.02 else "   <-- DISAGREES"
        print(f"  {name:26} coverage {cov:5.1%}  worst diff {worst:6.2%}{flag}")
        if not (worst < 0.02):
            ok = False

    # bubbles: just report coverage and the half-pixel offset we found
    for which in ("sharp", "fuzzy"):
        p = os.path.join(here, "..", "templates", "bubbles", f"bubbles_{which}.fits")
        if os.path.exists(p):
            t, cov = load(p)
            print(f"  bubbles_{which:<6} coverage {cov:5.1%}  "
                  f"integral after renorm {np.sum(normalise(t)*solid_angle()):.4f}")

    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    selftest()
