# shared PSF convolution. one implementation so every gap team blurs identically
#
# READ THIS FIRST - WHAT THIS IS AND ISN'T FOR
#
# do NOT pre-convolve templates you hand to fermipy. gtsrcmaps convolves a
# SpatialMap with the PSF and exposure itself when it builds source maps, so
# pre-convolving blurs twice. double blurring makes every template rounder,
# and rounder is exactly what gaps 1 and 4 are measuring.
#
# use this module for:
#   - comparing template SHAPES at LAT resolution (q, c4), where gtsrcmaps
#     is not in the loop
#   - putting competing templates on a common resolution floor
#
# THE PSF USED HERE IS NOW THE REAL ONE. psf_data/gtpsf_*.fits is a gtpsf run
# on our own livetime cube (P8R3_SOURCE_V3, evtype=3), and every function below
# interpolates that table when it is present. The old analytic King fit is kept
# only as a fallback, and validate_against_gtpsf() still checks whichever path
# is live.
#
# Why this mattered: the analytic fit agreed with gtpsf to <5% below 10 GeV but
# was 30% TOO NARROW at 50 GeV. A too-narrow model PSF over-predicts the peak
# of every point source, so the 10-50 GeV fit pushed the 4FGL normalisation
# down to ~0.5 and swallowed the excess with it - the GCE came out at TS 0
# there. Nothing errored; the numbers were just wrong.

import os

import numpy as np

PIX = 0.1          # deg, Conventions section 2
KING_GAMMA = 2.5   # shape param of the analytic FALLBACK only

# gtpsf table: psf(energy, theta) in sr^-1, plus its two axes. Loaded lazily
# from psf_data/ so importing this module stays cheap and no caller has to know
# about it.
_TABLE = None
_DEFAULT_TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "psf_data",
                              "gtpsf_P8R3_SOURCE_V3_evtype3.fits")


def load_gtpsf(path=None, required=False):
    """Load a gtpsf table. Returns the table dict, or None if there isn't one."""
    global _TABLE
    path = path or _DEFAULT_TABLE
    if _TABLE is not None and _TABLE["path"] == path:
        return _TABLE
    if not os.path.exists(path):
        if required:
            raise FileNotFoundError(path)
        return None
    from astropy.io import fits
    with fits.open(path) as hd:
        th = np.asarray(hd["THETA"].data["Theta"], float)          # deg
        e = np.asarray(hd["PSF"].data["ENERGY"], float) / 1000.0   # MeV->GeV
        psf = np.asarray(hd["PSF"].data["PSF"], float)             # sr^-1
    # cumulative containment per energy, on the theta grid
    w = 2.0 * np.pi * np.radians(th) * np.radians(np.gradient(th))
    cum = np.cumsum(psf * w, axis=1)
    cum /= cum[:, -1:]
    _TABLE = {"path": path, "theta": th, "energy": e, "psf": psf, "cum": cum}
    return _TABLE


def _interp_rows(tab, energy_gev):
    """Linear-in-log-energy weights for the two bracketing table rows."""
    le = np.log(tab["energy"])
    lt = float(np.log(np.clip(energy_gev, tab["energy"][0], tab["energy"][-1])))
    j = int(np.clip(np.searchsorted(le, lt), 1, len(le) - 1))
    f = (lt - le[j - 1]) / (le[j] - le[j - 1])
    return j - 1, j, f


def r68_deg(energy_gev):
    """68% containment radius. gtpsf table if available, else the old fit."""
    e = np.atleast_1d(np.asarray(energy_gev, float))
    tab = load_gtpsf()
    if tab is not None:
        return np.array([containment_radius(0.68, float(x)) for x in e])
    return np.sqrt((0.8 * e ** -0.8) ** 2 + 0.06 ** 2)


def _r68_analytic(energy_gev):
    e = np.atleast_1d(np.asarray(energy_gev, float))
    return np.sqrt((0.8 * e ** -0.8) ** 2 + 0.06 ** 2)


def _sigma_from_r68(r68, gamma=KING_GAMMA):
    # King containment: F(x) = 1 - (1 + x^2/(2*gamma*sigma^2))^(1-gamma)
    # solve F(r68) = 0.68 for sigma
    k = (1.0 - 0.68) ** (1.0 / (1.0 - gamma)) - 1.0
    return r68 / np.sqrt(2.0 * gamma * k)


def king(x_deg, energy_gev, gamma=KING_GAMMA):
    """PSF value at angular offset x, normalised so the 2D integral is 1.

    With a gtpsf table loaded this interpolates the real radial profile in
    log-log, which is a better shape than a single King at high energy where
    the LAT PSF has a wide tail the analytic form misses.
    """
    tab = load_gtpsf()
    if tab is not None:
        i, j, f = _interp_rows(tab, energy_gev)
        lp = ((1 - f) * np.log(np.maximum(tab["psf"][i], 1e-300))
              + f * np.log(np.maximum(tab["psf"][j], 1e-300)))
        x = np.asarray(x_deg, float)
        lth = np.log(np.maximum(tab["theta"], 1e-6))
        out = np.exp(np.interp(np.log(np.maximum(x, 1e-6)), lth, lp))
        out[x > tab["theta"][-1]] = 0.0
        return out
    s = _sigma_from_r68(_r68_analytic(energy_gev), gamma)[0]
    x = np.asarray(x_deg, float)
    norm = 1.0 / (2.0 * np.pi * s ** 2) * (1.0 - 1.0 / gamma)
    return norm * (1.0 + x ** 2 / (2.0 * gamma * s ** 2)) ** (-gamma)


def containment_radius(frac, energy_gev, gamma=KING_GAMMA):
    """Radius containing `frac` of the flux. Used by masks/make_masks.py."""
    tab = load_gtpsf()
    if tab is not None:
        i, j, f = _interp_rows(tab, energy_gev)
        cum = (1 - f) * tab["cum"][i] + f * tab["cum"][j]
        return float(np.interp(frac, cum, tab["theta"]))
    s = _sigma_from_r68(_r68_analytic(energy_gev), gamma)[0]
    k = (1.0 - frac) ** (1.0 / (1.0 - gamma)) - 1.0
    return float(s * np.sqrt(2.0 * gamma * k))


def kernel(energy_gev, pix=PIX, truncate_at_frac=0.995):
    """PSF kernel on the analysis grid, sized to hold `truncate_at_frac`."""
    r = containment_radius(truncate_at_frac, energy_gev)
    n = max(int(np.ceil(r / pix)), 1)
    y, x = np.mgrid[-n:n + 1, -n:n + 1] * pix
    k = king(np.hypot(x, y), energy_gev)
    return k / k.sum()


def convolve(template, energy_gev, pix=PIX):
    """Blur a template to LAT resolution at one energy. Flux conserving."""
    from scipy.signal import fftconvolve
    return fftconvolve(template, kernel(energy_gev, pix), mode="same")


def validate_against_gtpsf(gtpsf_fits):
    """Compare this module's containment radii to a real gtpsf run.

    Run this before trusting any shape number. Writes nothing - prints a table
    and returns the worst fractional disagreement so a caller can assert on it.

        gtpsf expcube=... outfile=psf.fits irfs=P8R3_SOURCE_V3 \
              ra=266.405 dec=-28.936 emin=500 emax=50000 nenergies=20 \
              thetamax=5 ntheta=300
    """
    from astropy.io import fits
    with fits.open(gtpsf_fits) as hd:
        th = np.asarray(hd["THETA"].data["Theta"], float)         # deg
        e = np.asarray(hd["PSF"].data["ENERGY"], float) / 1000.0  # MeV -> GeV
        psf = np.asarray(hd["PSF"].data["PSF"], float)            # sr^-1

    worst = 0.0
    print(f"{'E [GeV]':>9} {'r68 gtpsf':>10} {'r68 ours':>9} {'diff':>7}")
    for i, en in enumerate(e):
        w = 2.0 * np.pi * np.radians(th) * np.radians(np.gradient(th))
        cum = np.cumsum(psf[i] * w)
        cum /= cum[-1]
        r68_true = float(np.interp(0.68, cum, th))
        r68_ours = float(r68_deg(en)[0])
        d = abs(r68_ours - r68_true) / r68_true
        worst = max(worst, d)
        print(f"{en:9.2f} {r68_true:10.3f} {r68_ours:9.3f} {d:6.1%}")
    print(f"\nworst disagreement: {worst:.1%}")
    if worst > 0.10:
        print("MORE THAN 10% OFF - replace r68_deg() with the gtpsf values "
              "before using this for anything that goes in a paper")
    return worst


def selftest():
    ok = True

    # 1. kernel normalised
    for e in (0.75, 2.0, 6.0, 22.0):
        k = kernel(e)
        if abs(k.sum() - 1.0) > 1e-9:
            print(f"FAIL kernel at {e} GeV sums to {k.sum()}"); ok = False

    # 2. containment monotonic in fraction and decreasing with energy
    for e in (0.75, 2.0, 6.0, 22.0):
        r = [containment_radius(f, e) for f in (0.5, 0.68, 0.95, 0.99)]
        if not all(r[i] < r[i + 1] for i in range(len(r) - 1)):
            print(f"FAIL containment not monotonic at {e} GeV: {r}"); ok = False
    r95 = [containment_radius(0.95, e) for e in (0.75, 2.0, 6.0, 22.0)]
    if not all(r95[i] > r95[i + 1] for i in range(len(r95) - 1)):
        print(f"FAIL r95 not decreasing with energy: {r95}"); ok = False

    # 3. convolution conserves flux and blurs a point into the right width
    img = np.zeros((400, 400)); img[200, 200] = 1.0
    for e in (2.0, 6.0):
        out = convolve(img, e)
        if abs(out.sum() - 1.0) > 2e-3:
            print(f"FAIL convolution lost flux at {e} GeV: {out.sum():.4f}")
            ok = False
        # measure r68 back off the convolved point source
        yy, xx = np.mgrid[0:400, 0:400]
        rr = np.hypot(xx - 200, yy - 200) * PIX
        order = np.argsort(rr.ravel())
        cum = np.cumsum(out.ravel()[order]); cum /= cum[-1]
        r68_meas = float(rr.ravel()[order][np.searchsorted(cum, 0.68)])
        r68_want = float(r68_deg(e)[0])
        if abs(r68_meas - r68_want) > 0.12:
            print(f"FAIL recovered r68 {r68_meas:.3f} vs wanted {r68_want:.3f} "
                  f"at {e} GeV"); ok = False

    print(f"  {'bin':>12} {'r68':>7} {'r95':>7} {'kernel':>9}")
    for (lo, hi), e in zip([(0.5,1),(1,3),(3,10),(10,50)],
                           [0.707, 1.732, 5.477, 22.36]):
        print(f"  {str(lo)+'-'+str(hi)+' GeV':>12} {r68_deg(e)[0]:7.3f} "
              f"{containment_radius(0.95,e):7.3f} {kernel(e).shape[0]:5d} px")
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    selftest()
