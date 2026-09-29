"""
Cuts the galactic diffuse model down to something small enough to move.

gll_iem_v07.fits is all-sky and ~3.7 GB, well past the file transfer cap. We
only ever look at |l|,|b| < 20, so this takes a slightly padded box around the
GC at the model's own 0.125 deg resolution and writes it out. Regridding onto
the project 0.1 deg grid happens on the other side with utils/regrid.py, which
is already cross-checked against the two hand-written loaders.

numpy + astropy only, no scipy - keeps it runnable in the bare fermi env.
"""

import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

SRC = "gll_iem_v07.fits"
OUT = "gap1_iem_roi.fits"
PAD = 21.0          # deg half-width, 1 deg past the ROI so interpolation at
                    # the edge has something to chew on


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else SRC
    out = sys.argv[2] if len(sys.argv) > 2 else OUT

    with fits.open(src, memmap=True) as hd:
        hdu = hd[0]
        data = hdu.data
        hdr = hdu.header
        print(f"source shape {data.shape} dtype {data.dtype}")
        print(f"CRVAL1={hdr.get('CRVAL1')} CRPIX1={hdr.get('CRPIX1')} "
              f"CDELT1={hdr.get('CDELT1')}")
        print(f"CRVAL2={hdr.get('CRVAL2')} CRPIX2={hdr.get('CRPIX2')} "
              f"CDELT2={hdr.get('CDELT2')}")

        if data.ndim != 3:
            sys.exit(f"expected a 3D cube, got {data.ndim}D")
        nE, ny, nx = data.shape

        w = WCS(hdr).celestial
        cd1 = float(hdr["CDELT1"])
        cd2 = float(hdr["CDELT2"])

        # latitude: straightforward, no wrap
        yb = []
        for b in (-PAD, PAD):
            _, y = w.wcs_world2pix(0.0, b, 0)
            yb.append(float(y))
        y0 = int(np.floor(min(yb))); y1 = int(np.ceil(max(yb)))
        y0 = max(y0, 0); y1 = min(y1, ny - 1)
        print(f"lat pixel range {y0}..{y1}  ({y1 - y0 + 1} rows)")

        # longitude wraps, so work out the centre pixel and step outward by a
        # fixed count rather than trusting world2pix either side of the seam
        xc, _ = w.wcs_world2pix(0.0, 0.0, 0)
        xc = float(xc)
        halfpix = int(np.ceil(PAD / abs(cd1)))
        xidx = (np.round(xc).astype(int) + np.arange(-halfpix, halfpix + 1)) % nx
        print(f"lon centre pixel {xc:.2f}, taking +-{halfpix} px "
              f"({xidx.size} cols), wrapped")

        sub = np.asarray(data[:, y0:y1 + 1, :], dtype=np.float32)
        sub = sub[:, :, xidx]
        print(f"cut shape {sub.shape}  "
              f"{sub.nbytes / 1e6:.1f} MB  "
              f"min={np.nanmin(sub):.3e} max={np.nanmax(sub):.3e}")

        # new header: same projection and pixel size, recentred on the GC
        h = fits.Header()
        h["CTYPE1"] = hdr.get("CTYPE1", "GLON-CAR")
        h["CTYPE2"] = hdr.get("CTYPE2", "GLAT-CAR")
        h["CDELT1"] = cd1
        h["CDELT2"] = cd2
        h["CRVAL1"] = 0.0
        h["CRVAL2"] = 0.0
        # centre pixel of the cut, 1-based
        h["CRPIX1"] = float(halfpix + 1)
        _, ygc = w.wcs_world2pix(0.0, 0.0, 0)
        h["CRPIX2"] = float(ygc) - y0 + 1.0
        for k in ("CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3", "BUNIT"):
            if k in hdr:
                h[k] = hdr[k]
        h["HISTORY"] = f"ROI cut from {src} by reduce_iem.py, |l|,|b|<{PAD}"
        h["SRCNAXIS"] = (f"{nx}x{ny}x{nE}", "original all-sky dimensions")

        hdus = [fits.PrimaryHDU(sub, h)]
        # carry the energy table across verbatim - the fit needs it
        for e in hd[1:]:
            if e.name.upper() in ("ENERGIES", "EBOUNDS"):
                print(f"copying extension {e.name}")
                hdus.append(e.copy())
        fits.HDUList(hdus).writeto(out, overwrite=True)

    import os
    print(f"wrote {out}  {os.path.getsize(out) / 1e6:.1f} MB")
    print("OK")


if __name__ == "__main__":
    main()
