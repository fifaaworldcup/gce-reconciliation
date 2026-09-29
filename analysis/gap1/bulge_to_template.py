# loads the Coleman/Macias bulge templates onto our grid
# source: github.com/chrisgordon1/galactic_bulge_templates (Coleman et al 2019,
# arXiv:1911.04714). two files:
#   Bulge_modulated_Coleman_etal_2019_Normalized.fits - non-parametric (VVV) bulge
#   F98_BoxyBulge_arxv1611.06644_Normalized.fits      - Freudenreich 98 S-bulge
#
# their grid: 200x200, CDELT 0.2 deg, CRPIX 100.5, CRVAL 0 -> centres at +-19.9
# our grid:   400x400, CDELT 0.1 deg                      -> centres at +-19.95
# so their centres sit half of OUR pixel off ours. that is not a 2x block
# upsample, it needs interpolation, or the map shifts by 0.05 deg.
#
# CDELT1 is negative in both, same as ours, so no longitude flip here (unlike
# the HESTIA files, which are +0.1 and do need flipping).
#
# both files sum to 82070.2, i.e. unit integral with a FLAT pixel area and no
# cos(b) factor - same convention Muru used. we renormalise with cos(b) so
# every template in this project integrates the same way.

import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hestia_to_template as H

SRC_NPIX = 200
SRC_PIX = 0.2
SRC_HALF_CENTRE = 19.9          # |l| of the outermost pixel centre

FILES = {
    "coleman": "Bulge_modulated_Coleman_etal_2019_Normalized.fits",
    "f98": "F98_BoxyBulge_arxv1611.06644_Normalized.fits",
}


def read(path):
    with fits.open(path) as hd:
        arr = np.asarray(hd[0].data, float)
        hdr = hd[0].header
    if arr.shape != (SRC_NPIX, SRC_NPIX):
        raise ValueError(f"{path}: expected {SRC_NPIX}^2, got {arr.shape}")
    for k, want in (("CDELT1", -SRC_PIX), ("CDELT2", SRC_PIX),
                    ("CRPIX1", 100.5), ("CRPIX2", 100.5),
                    ("CRVAL1", 0.0), ("CRVAL2", 0.0)):
        if k in hdr and not np.isclose(float(hdr[k]), want):
            raise ValueError(f"{path}: {k}={hdr[k]} not {want}; grid changed")
    return arr


def to_conventions_grid(arr_bl):
    """Interpolate the 0.2 deg map onto our 0.1 deg grid. No flip needed."""
    from scipy.ndimage import map_coordinates
    l = (H.L_HALF - H.PIX / 2.0) - H.PIX * np.arange(H.NPIX)    # +19.95 -> -19.95
    b = (-H.B_HALF + H.PIX / 2.0) + H.PIX * np.arange(H.NPIX)
    ll, bb = np.meshgrid(l, b)
    # fractional index into their grid: index 0 is l=+19.9, l decreases
    li = (SRC_HALF_CENTRE - ll) / SRC_PIX
    bi = (bb + SRC_HALF_CENTRE) / SRC_PIX
    return map_coordinates(arr_bl, np.stack([bi, li]), order=1, mode="nearest")


def build(which, root="../../templates/bulge"):
    path = os.path.join(root, FILES[which])
    t = H.normalise(to_conventions_grid(read(path)))
    return t, {
        "TEMPNAME": f"bulge_{which}",
        "BULGEMOD": (which, "coleman=non-parametric VVV, f98=Freudenreich S-bulge"),
        "SRCFILE": FILES[which][:60],
        "SRCREF": "arXiv:1911.04714 via github.com/chrisgordon1",
        "SRCPIX": (SRC_PIX, "deg, interpolated up to 0.1"),
        "LONFLIP": (False, "source CDELT1 already negative"),
    }


def selftest():
    ok = True
    # an asymmetric synthetic map must survive the regrid in the right place
    l = SRC_HALF_CENTRE - SRC_PIX * np.arange(SRC_NPIX)
    b = -SRC_HALF_CENTRE + SRC_PIX * np.arange(SRC_NPIX)
    L, B = np.meshgrid(l, b)
    fake = np.exp(-((L / 8.0) ** 2 + ((B - 4.0) / 3.0) ** 2))
    out = to_conventions_grid(fake)
    lc = (H.L_HALF - H.PIX / 2.0) - H.PIX * np.arange(H.NPIX)
    bc = (-H.B_HALF + H.PIX / 2.0) + H.PIX * np.arange(H.NPIX)
    r, c = np.unravel_index(np.argmax(out), out.shape)
    if abs(lc[c]) > 0.11 or abs(bc[r] - 4.0) > 0.11:
        print(f"FAIL regrid peak at l={lc[c]:+.2f} b={bc[r]:+.2f}, want (0, +4)")
        ok = False
    wl = out[r, :].sum() / out[r, :].max()
    wb = out[:, c].sum() / out[:, c].max()
    if not wl > 1.5 * wb:
        print(f"FAIL elongation lost or transposed (l {wl:.1f} vs b {wb:.1f})")
        ok = False
    print(f"  regrid peak l={lc[c]:+.2f} b={bc[r]:+.2f}; widths l {wl:.0f} b {wb:.0f}")
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    selftest()
    for which in FILES:
        t, hdr = build(which)
        out = f"../../templates/gap1_bulge_{which}.fits"
        print(f"  {which:9} -> {H.write_fits(t, out, hdr)}")
