# turns muru's hestia maps into template fits files on our grid
# their files: 800x800, -39.95 to 39.95 deg, 0.1 steps, already squared (dm2)

import datetime, os, re
import numpy as np
from astropy.io import fits

PIX = 0.1
NPIX = 400
L_HALF = B_HALF = 20.0
SRC_NPIX = 800
SRC_HALF = 40.0
CROP = slice(200, 600)   # their grid lines up with ours so its just a crop

GALAXIES = {
    ("09_18", "127000000000003"): "G1.1",
    ("17_11", "127000000000003"): "G2.1",
    ("37_11", "127000000000002"): "G3.1",
    ("09_18", "127000000000002"): "G1.2",
    ("17_11", "127000000000002"): "G2.2",
    ("37_11", "127000000000001"): "G3.2",
}

_PAT = re.compile(r"density_projection_(?P<sim>\d+_\d+)_(?P<npart>\d+)_halo(?P<halo>\d+)"
                  r"_angle(?P<angle>[\d.]+)_dm2_(?P<smooth>[\d.]+)deg\.(?P<ext>fits|tsv)$")


def parse_name(path):
    m = _PAT.match(os.path.basename(path))
    if not m:
        raise ValueError(f"{path}: bad filename")
    g = m.groupdict()
    g["galaxy"] = GALAXIES.get((g["sim"], g["halo"]), "unknown")
    g["angle"] = float(g["angle"])
    g["smooth"] = float(g["smooth"])
    return g


def read_projection(path):
    """returns the map as [b, l]"""
    meta = parse_name(path)
    if meta["ext"] == "fits":
        with fits.open(path) as hd:
            arr, hdr = np.asarray(hd[0].data, float), hd[0].header
        # bail out if a future file isnt on the same grid
        for k, want in (("CRPIX1", 400.5), ("CRPIX2", 400.5), ("CDELT1", 0.1),
                        ("CDELT2", 0.1), ("CRVAL1", 0.0), ("CRVAL2", 0.0)):
            if k in hdr and not np.isclose(float(hdr[k]), want):
                raise ValueError(f"{path}: {k}={hdr[k]} not {want}, grid changed")
    else:
        arr = np.loadtxt(path).T   # the tsv is the transpose of the fits, checked
    if arr.shape != (SRC_NPIX, SRC_NPIX):
        raise ValueError(f"{path}: got {arr.shape}")
    return arr, meta


def to_conventions_grid(arr_bl):
    # their CDELT1 is +0.1, ours is -0.1, so flip l or the map ends up mirrored
    return arr_bl[CROP, CROP][:, ::-1]


def solid_angle():
    b = (-B_HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    return np.radians(PIX) ** 2 * np.cos(np.radians(b))[:, None] * np.ones((1, NPIX))


def normalise(t):
    total = float(np.sum(t * solid_angle()))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("empty template")
    return t / total


def tophat_smooth(t, radius_deg):
    # same top hat muru uses. put every template through this before comparing
    from scipy.signal import fftconvolve
    r = int(round(radius_deg / PIX))
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    k = ((x ** 2 + y ** 2) <= r ** 2).astype(float)
    return fftconvolve(t, k / k.sum(), mode="same")


def build(path):
    # no squaring here, dm2 means the files are already squared
    arr, meta = read_projection(path)
    t = normalise(to_conventions_grid(arr))
    return t, {
        "TEMPNAME": f"hestia_{meta['galaxy']}",
        "GALAXY": meta["galaxy"],
        "SIMID": meta["sim"],
        "HALOID": meta["halo"],
        "VIEWANG": (meta["angle"], "deg"),
        "SMOOTHIN": (meta["smooth"], "deg, already in the source map"),
        "RHO2CONV": ("square_of_column", "dm2 is (column)^2 not int rho^2 ds"),
        "SRCFILE": os.path.basename(path)[:60],
        "LONFLIP": (True, "l flipped for CDELT1 < 0"),
    }


def write_fits(t, path, extra=None):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (NPIX + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = -PIX, PIX
    h["EQUINOX"] = 2000.0
    h["NORMED"] = (True, "unit integral over roi")
    h["PSFCONV"] = (False, "not psf convolved, use utils/psf.py")
    h["DATE"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    for k, v in (extra or {}).items():
        h[k] = v
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fits.PrimaryHDU(np.asarray(t, np.float32), h).writeto(path, overwrite=True)
    return path


def _shape(x, y):
    x, y = x - x.mean(), y - y.mean()
    w, v = np.linalg.eigh(np.cov(np.stack([x, y])))
    q = float(np.sqrt(w.min() / w.max()))
    xy = (np.stack([x, y]).T @ v) / np.sqrt(w)
    r, phi = np.hypot(xy[:, 0], xy[:, 1]), np.arctan2(xy[:, 1], xy[:, 0])
    sel = r >= np.percentile(r, 97) * 0.9
    if sel.sum() < 50:
        return q, float("nan")
    c4 = float(-2 * np.mean((r[sel] / r[sel].mean() - 1) * np.cos(4 * phi[sel])))
    return q, c4


def shape_on_source(arr_bl, level_frac):
    """q and c4 on the full +-40 map. returns (q, c4, mean_r, max_r)

    use this, not the cropped roi. a square crop cuts the contour at the edges
    but not the corners which fakes a 4 fold signal, ie fake boxiness. on these
    files the crop gives c4 = +0.013..+0.021 (looks boxy) and the same contour
    uncropped gives -0.006..+0.002 (plain ellipse).

    c4 > 0 is boxy, < 0 disky, 0 is an ellipse. checked on superellipses:
    n=2 -> 0.000, n=4 -> +0.019, n=1.3 -> -0.025
    """
    ls = (-SRC_HALF + PIX / 2.0) + PIX * np.arange(SRC_NPIX)
    LL, BB = np.meshgrid(ls, ls)
    m = arr_bl >= level_frac * arr_bl.max()
    q, c4 = _shape(LL[m], BB[m])
    th = np.hypot(LL, BB)[m]
    return q, c4, float(th.mean()), float(th.max())


def axis_ratio(t, level_frac=0.5):
    # muru's svd axis ratio, for comparing to the paper. cant see boxiness
    l = (L_HALF - PIX / 2.0) - PIX * np.arange(NPIX)
    b = (-B_HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    ll, bb = np.meshgrid(l, b)
    m = t >= level_frac * t.max()
    return _shape(ll[m], bb[m])[0]


def boxiness(t, level_frac=0.5):
    # cropped version, only here to show the roi artifact. use shape_on_source
    l = (L_HALF - PIX / 2.0) - PIX * np.arange(NPIX)
    b = (-B_HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    ll, bb = np.meshgrid(l, b)
    m = t >= level_frac * t.max()
    return _shape(ll[m], bb[m])[1]


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        t, hdr = build(p)
        arr, _ = read_projection(p)
        q, c4, r, _ = shape_on_source(arr, 0.05)
        out = (f"../../templates/gap1_hestia_{hdr['GALAXY']}"
               f"_ang{int(hdr['VIEWANG'][0])}_{hdr['SMOOTHIN'][0]}deg.fits")
        print(f"{hdr['GALAXY']:5} ang{hdr['VIEWANG'][0]:5.1f}  q={q:.3f} c4={c4:+.3f} "
              f"(r~{r:.0f}deg)  -> {write_fits(t, out, hdr)}")
