# builds masks A, B, C, D from a catalog, as likelihood WEIGHTS MAPS for gtlike
# gtlike has no pixel-mask option so 0 = excluded, 1 = kept
#
# masks are energy dependent because the 95% containment radius shrinks with
# energy, so the output is a 3D cube matching the analysis energy bins

import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import psf as PSF   # single source of truth for containment radii

PIX = 0.1
NPIX = 400
L_HALF = B_HALF = 20.0

# Conventions section 1 analysis bins, GeV
EBINS = [(0.5, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 50.0)]

# 95% containment per bin, taken from utils/psf.py at the geometric mean energy
# so the masks and the shape measurements can never drift apart.
# still approximate until psf.validate_against_gtpsf() has been run.
EMID = [np.sqrt(lo * hi) for lo, hi in EBINS]
PSF95_DEG = {i: PSF.containment_radius(0.95, EMID[i]) for i in range(len(EBINS))}

MASK_B_TS_MIN = 25.0
MASK_B_FLUX_MIN = 5e-9        # ph cm^-2 s^-1

# ---------------------------------------------------------------- mask D
# Masks A and C cut a fixed 95%-containment disc around every catalogue source.
# With 4FGL-DR4 that is 761 discs over our ROI and it keeps 8.7% of the sky,
# the inner degrees included, and no template comparison inside it is decisive
# (mini-paper section 6.1). A fixed radius is the problem: it treats the
# faintest source in the catalogue like the brightest.
#
# Mask D scales the radius by the source itself. Mask out to the radius where
# the source's PSF-spread surface brightness drops to MASK_E_FRACTION times the
# local diffuse background:
#
#     F_s(E) * PSF_E(r)  =  MASK_E_FRACTION * B_E
#
# with F_s the source's own flux in that energy bin from its catalogue spectrum
# and B_E the diffuse model's intensity AT THAT SOURCE'S POSITION. Both sides
# are ph cm^-2 s^-1 sr^-1, so the comparison is physical rather than tuned.
# Bright sources get wide discs, faint ones get almost none, the radius shrinks
# with energy through the PSF, and a source sitting on the bright Galactic
# plane has to be brighter to earn the same disc as one at high latitude. Using
# a single ROI-average background instead would hold those two to the same
# threshold, which is most of what is wrong with a fixed radius.
# FRACTION FIXED BEFORE ANY MASK-D FIT, on a contamination argument, not on an
# outcome: 0.1 means "mask wherever a catalogue source contributes 10% or more
# of the local diffuse background". That keeps 87% of the ROI, comparable to
# mask B's 84%, but with radii set by each source instead of one number for
# all of them. Sensitivity to this choice is reported in the mini-paper
# (1.0 -> 98% kept, 0.3 -> 95%, 0.1 -> 87%, 0.03 -> 70%, 0.01 -> 50%).
MASK_E_FRACTION = 0.1
MASK_E_RMIN = 0.10            # deg, one pixel
MASK_E_RMAX = 5.0             # deg, don't let one bright source eat the ROI


def grid():
    l = (L_HALF - PIX / 2.0) - PIX * np.arange(NPIX)     # +19.95 -> -19.95
    b = (-B_HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    return np.meshgrid(l, b)


def read_catalog(path):
    """Pull glon, glat, ts, flux out of a 4FGL/FL16Y style catalog FITS.

    Column names differ between releases so we try the usual spellings and
    fail loudly rather than guessing.
    """
    with fits.open(path) as hd:
        t = hd[1].data
        cols = {c.upper(): c for c in hd[1].columns.names}

    def pick(*names):
        for n in names:
            if n.upper() in cols:
                return np.asarray(t[cols[n.upper()]], float)
        raise KeyError(f"catalog has none of {names}; columns are {sorted(cols)}")

    glon = pick("GLON")
    glat = pick("GLAT")
    ts = pick("Signif_Avg", "TS", "Test_Statistic")
    # Signif_Avg is sigma, TS is sigma^2 - normalise to TS
    if "SIGNIF_AVG" in cols and "TS" not in cols:
        ts = ts ** 2
    flux = pick("Flux1000", "Energy_Flux100", "Flux_Density")
    glon = np.where(glon > 180.0, glon - 360.0, glon)
    return glon, glat, ts, flux


def _disc_mask(ll, bb, glon, glat, radii):
    """0 inside any source disc, 1 outside. Angular distance, not euclidean -
    at |b| up to 20 deg a flat sqrt(dl^2+db^2) is wrong by several percent.

    The first version tested every source against every pixel. Fine for the
    three-source self-test, hopeless for 4FGL-DR4 (7195 sources x 160k pixels
    x 8 bins, timed out). Now sources that cannot reach the ROI are dropped
    and each remaining disc only touches its own bounding box of pixels.
    Same answer, a few hundred times less work.
    """
    keep = np.ones_like(ll, dtype=bool)
    glon = np.asarray(glon, float); glat = np.asarray(glat, float)
    radii = np.asarray(radii, float)
    if glon.size == 0:
        return keep

    l_axis = ll[0, :]                  # descending in l
    b_axis = bb[:, 0]                  # ascending in b
    pix = abs(b_axis[1] - b_axis[0])
    lmax, lmin = l_axis.max(), l_axis.min()
    bmax, bmin = b_axis.max(), b_axis.min()

    dl0 = (glon + 180.0) % 360.0 - 180.0
    for lo, la, r in zip(dl0, glat, radii):
        # the l half-width of a disc grows as 1/cos(b); take the widest row
        cmin = np.cos(np.radians(min(abs(la) + r, 89.0)))
        rl = r / max(cmin, 1e-3)
        if (lo - rl > lmax + pix or lo + rl < lmin - pix or
                la - r > bmax + pix or la + r < bmin - pix):
            continue
        i0 = max(int(np.floor((la - r - bmin) / pix)) - 1, 0)
        i1 = min(int(np.ceil((la + r - bmin) / pix)) + 2, bb.shape[0])
        j0 = max(int(np.floor((lmax - (lo + rl)) / pix)) - 1, 0)
        j1 = min(int(np.ceil((lmax - (lo - rl)) / pix)) + 2, ll.shape[1])
        if i0 >= i1 or j0 >= j1:
            continue
        sl = (slice(i0, i1), slice(j0, j1))
        cb = np.cos(np.radians(bb[sl]))
        dl = (ll[sl] - lo + 180.0) % 360.0 - 180.0
        sep = np.hypot(dl * cb, bb[sl] - la)      # small-angle, exact enough here
        keep[sl] &= sep > r
    return keep


def _disc_mask_bruteforce(ll, bb, glon, glat, radii):
    """The original all-pixels version, kept only so the self-test can prove
    the fast one gives the identical mask."""
    keep = np.ones_like(ll, dtype=bool)
    cb = np.cos(np.radians(bb))
    for lo, la, r in zip(glon, glat, radii):
        dl = (ll - lo + 180.0) % 360.0 - 180.0
        sep = np.hypot(dl * cb, bb - la)
        keep &= sep > r
    return keep


def mask_a(ll, bb, cat, r95, **_):
    glon, glat, ts, flux = cat
    r = np.full(glon.shape, r95 + 0.5)
    return _disc_mask(ll, bb, glon, glat, r)


def mask_b(ll, bb, cat, r95, **_):
    glon, glat, ts, flux = cat
    sel = (ts > MASK_B_TS_MIN) & (flux > MASK_B_FLUX_MIN)
    r = np.full(int(sel.sum()), r95 + 0.5)
    return _disc_mask(ll, bb, glon[sel], glat[sel], r)


def mask_c(ll, bb, cat, r95, **_):
    """Mask A plus the latitude cut: |b| < 2 deg for |l| < 15, tapering
    linearly to |b| < 0.5 deg at |l| = 20."""
    keep = mask_a(ll, bb, cat, r95)
    al = np.abs(ll)
    bcut = np.where(al <= 15.0, 2.0,
                    np.clip(2.0 + (al - 15.0) * (0.5 - 2.0) / 5.0, 0.5, 2.0))
    return keep & (np.abs(bb) >= bcut)


def _bkg_at_sources(iem_plane, cat_table):
    """Diffuse intensity at each catalogue source's pixel. Sources outside the
    ROI fall back to the ROI median, which is all they need - their disc only
    matters where it clips our edge."""
    glon = np.asarray(cat_table["GLON"], float)
    glon = np.where(glon > 180, glon - 360, glon)
    glat = np.asarray(cat_table["GLAT"], float)
    j = np.floor((L_HALF - glon) / PIX).astype(int)
    i = np.floor((glat + B_HALF) / PIX).astype(int)
    inside = (i >= 0) & (i < NPIX) & (j >= 0) & (j < NPIX)
    out = np.full(glon.shape, float(np.median(iem_plane)))
    out[inside] = iem_plane[i[inside], j[inside]]
    return out


def mask_e_radii(flux_bin, energy_gev, bkg_intensity,
                 fraction=MASK_E_FRACTION, rmin=MASK_E_RMIN, rmax=MASK_E_RMAX):
    """Per-source mask radius for mask D. flux_bin in ph cm^-2 s^-1 for the
    bin, bkg_intensity in ph cm^-2 s^-1 sr^-1. Returns radii in degrees.

    PSF_E(r) is monotonically decreasing in r, so invert it by interpolating on
    a log grid rather than solving per source.
    """
    r = np.logspace(np.log10(rmin), np.log10(rmax), 400)
    k = np.asarray(PSF.king(r, float(energy_gev)), float)      # sr^-1
    flux_bin = np.asarray(flux_bin, float)
    b = np.broadcast_to(np.asarray(bkg_intensity, float), flux_bin.shape)
    target = np.where(flux_bin > 0,
                      fraction * b / np.maximum(flux_bin, 1e-300),
                      np.inf)
    # k is decreasing; np.interp needs increasing x, so reverse both
    out = np.interp(target, k[::-1], r[::-1], left=rmax, right=rmin)
    return np.clip(out, rmin, rmax)


def mask_e(ll, bb, cat, r95, flux_bin=None, energy_gev=None, bkg=None):
    """Flux-scaled point-source mask. Needs the per-bin source fluxes and a
    background level, so it is only reachable through build_at_energies().
    """
    if flux_bin is None or energy_gev is None or bkg is None:
        raise ValueError("mask D needs flux_bin, energy_gev and bkg - "
                         "use build_at_energies(which='D', iem=...)")
    glon, glat, ts, flux = cat
    r = mask_e_radii(flux_bin, energy_gev, bkg)
    return _disc_mask(ll, bb, glon, glat, r)



# Mask D: Cholis et al 2022 (arXiv:2112.09706) Table III, the mask the
# Family 2 IEM library was built against.
#   TS <= 49 -> theta_s,  TS > 49 -> theta_l = (10/3) theta_s,  both
#   energy-dependent, plus |b| < 2 deg at all longitudes and energies.
# The disc cut matters: masks A, B and E have none, and without one a
# flattened signal template can absorb plane residual freely.
# Verified against two independent fetches; theta_l/theta_s = 10/3 in all
# fourteen rows.
CHOLIS_T3 = [
    # e_lo GeV, e_hi GeV, theta_s deg, theta_l deg
    (0.275, 0.357, 1.125, 3.750),
    (0.357, 0.464, 0.975, 3.250),
    (0.464, 0.603, 0.788, 2.630),
    (0.603, 0.784, 0.600, 2.000),
    (0.784, 1.020, 0.450, 1.500),
    (1.020, 1.320, 0.375, 1.250),
    (1.320, 1.720, 0.300, 1.000),
    (1.720, 2.240, 0.225, 0.750),
    (2.240, 2.910, 0.188, 0.625),
    (2.910, 3.780, 0.162, 0.540),
    (3.780, 4.910, 0.125, 0.417),
    (4.910, 10.800, 0.100, 0.333),
    (10.800, 23.700, 0.060, 0.200),
    (23.700, 51.900, 0.053, 0.175),
]
MASK_D_TS_SPLIT = 49.0
MASK_D_BCUT_DEG = 2.0


def cholis_radii(energy_gev):
    """(theta_s, theta_l) at any bin energy, log-log interpolated across
    Table III so a binning unlike theirs does not jump between rows."""
    e = np.array([np.sqrt(r[0] * r[1]) for r in CHOLIS_T3])
    ts_ = np.array([r[2] for r in CHOLIS_T3])
    tl_ = np.array([r[3] for r in CHOLIS_T3])
    x = np.log(np.clip(float(energy_gev), e[0], e[-1]))
    return (float(np.exp(np.interp(x, np.log(e), np.log(ts_)))),
            float(np.exp(np.interp(x, np.log(e), np.log(tl_)))))


def mask_d(ll, bb, cat, r95, energy_gev=None, **_):
    """Cholis 2022 Table III: TS-split radii plus the |b| < 2 deg disc cut."""
    glon, glat, ts, flux = cat
    if energy_gev is None:
        raise ValueError("mask D needs the bin energy")
    th_s, th_l = cholis_radii(energy_gev)
    r = np.where(np.asarray(ts, float) > MASK_D_TS_SPLIT, th_l, th_s)
    keep = _disc_mask(ll, bb, glon, glat, r)
    return keep & (np.abs(bb) >= MASK_D_BCUT_DEG)


BUILDERS = {"A": mask_a, "B": mask_b, "C": mask_c, "D": mask_d,
            "E": mask_e}


def build(catalog_path, which):
    ll, bb = grid()
    cat = read_catalog(catalog_path)
    cube = np.stack([
        BUILDERS[which](ll, bb, cat, PSF95_DEG[i],
                        **({"energy_gev": float(np.sqrt(EBINS[i][0] * EBINS[i][1]))}
                           if which == "D" else {})).astype(np.float32)
        for i in range(len(EBINS))])
    return cube


def build_at_energies(catalog_path, which, energies_gev, e_lo=None, e_hi=None,
                      iem=None):
    """same as build() but for whatever binning the counts cube actually has.
    the conventions 4 bins are a reporting choice, the fit doesnt have to match
    them, and the mask radius has to follow the real bin energy either way.

    which="E" also needs the bin edges and a diffuse cube on our grid
    (iem[e] in ph cm^-2 s^-1 MeV^-1 sr^-1 per bin, as fit_cubes.read_iem
    returns), because its radius comes from flux against background.
    """
    ll, bb = grid()
    cat = read_catalog(catalog_path)
    kw_per_bin = [{} for _ in energies_gev]
    if which == "E":
        if e_lo is None or e_hi is None or iem is None:
            raise ValueError("mask E needs e_lo, e_hi and iem")
        sys.path.insert(0, os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "gap1"))
        import ps_template as PS      # the validated 4FGL spectral integration
        with fits.open(catalog_path) as hd:
            t = hd[1].data
            for i, (lo, hi) in enumerate(zip(e_lo, e_hi)):
                # bin-integrated source flux, and the median diffuse intensity
                # over the same bin, both ph cm^-2 s^-1 (sr^-1 for the diffuse)
                fb = PS.band_flux(t, lo * 1e3, hi * 1e3)
                dE = (hi - lo) * 1e3                      # MeV
                bkg = _bkg_at_sources(iem[i], t) * dE
                kw_per_bin[i] = {"flux_bin": fb,
                                 "energy_gev": float(energies_gev[i]),
                                 "bkg": bkg}
    cube = []
    for i, e in enumerate(energies_gev):
        r95 = PSF.containment_radius(0.95, float(e))
        kw = dict(kw_per_bin[i])
        if which == "D":
            kw["energy_gev"] = float(e)
        cube.append(BUILDERS[which](ll, bb, cat, r95,
                                    **kw).astype(np.float32))
    return np.stack(cube)


def write(cube, path, which, catalog_path):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (NPIX + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = -PIX, PIX
    h["CTYPE3"] = "ENERGY"
    h["MASK"] = (which, "Conventions section 3")
    h["CATALOG"] = (catalog_path.split("/")[-1][:40], "source list used")
    h["PSF95APX"] = (True, "PSF95 radii approximate - validate vs gtpsf")
    fits.PrimaryHDU(cube, h).writeto(path, overwrite=True)
    return path


def selftest():
    ok = True
    ll, bb = grid()

    # synthetic catalog: one bright source dead centre, one faint at (10, 5),
    # one just outside the ROI that should still clip the corner
    glon = np.array([0.0, 10.0, 21.0])
    glat = np.array([0.0, 5.0, 0.0])
    ts = np.array([10000.0, 30.0, 500.0])
    flux = np.array([1e-7, 1e-9, 8e-9])
    cat = (glon, glat, ts, flux)

    a = mask_a(ll, bb, cat, PSF95_DEG[1])
    b = mask_b(ll, bb, cat, PSF95_DEG[1])
    c = mask_c(ll, bb, cat, PSF95_DEG[1])

    r = PSF95_DEG[1] + 0.5
    # centre pixel must be masked in A
    ci = NPIX // 2
    if a[ci, ci]:
        print("FAIL centre not masked in A"); ok = False
    # a point 1.5*r away from every source must survive
    if not a[ci, ci + int(2 * r / PIX)]:
        print("FAIL A masks too wide"); ok = False
    # B keeps the faint source (flux 1e-9 < 5e-9) so (10,5) stays unmasked
    j = int(round((L_HALF - PIX / 2.0 - 10.0) / PIX))
    i = int(round((5.0 + B_HALF - PIX / 2.0) / PIX))
    if a[i, j]:
        pass
    if not b[i, j]:
        print("FAIL B masked a source below its flux threshold"); ok = False
    # C must remove the plane
    if c[ci, ci + 100]:      # (l=-10, b~0)
        print("FAIL C did not cut the plane at |l|<15"); ok = False
    # C taper, tested on its own so nearby sources can't confound it.
    # |b| cut should be 2.0 deg inside |l|=15 and 0.5 deg at |l|=20.
    empty = (np.array([]), np.array([]), np.array([]), np.array([]))
    c_only = mask_c(ll, bb, empty, PSF95_DEG[1])
    def kept(l_deg, b_deg):
        j = int(round((L_HALF - PIX / 2.0 - l_deg) / PIX))
        i = int(round((b_deg + B_HALF - PIX / 2.0) / PIX))
        return bool(c_only[i, j])
    for l_deg, b_deg, want, why in [
            (5.0, 1.0, False, "|b|=1 inside |l|=15 must be cut"),
            (5.0, 3.0, True,  "|b|=3 inside |l|=15 must survive"),
            (19.9, 1.0, True, "|b|=1 at |l|=20 must survive the taper"),
            (19.9, 0.2, False, "|b|=0.2 at |l|=20 must still be cut")]:
        if kept(l_deg, b_deg) != want:
            print(f"FAIL C taper: {why}"); ok = False
    # ordering: C is a subset of A is a subset of the full sky
    if not (c.sum() <= a.sum() <= b.sum()):
        print(f"FAIL expected C <= A <= B, got {c.sum()} {a.sum()} {b.sum()}")
        ok = False

    # fast disc mask must equal the brute-force one exactly, including sources
    # outside the ROI whose discs clip it, high-|b| sources, and l near +-180
    rng = np.random.default_rng(3)
    gl = np.concatenate([rng.uniform(-25, 25, 300), [21.5, -21.5, 179.0, -179.5, 0.0]])
    gb = np.concatenate([rng.uniform(-25, 25, 300), [0.0, 19.0, 0.0, 5.0, 21.8]])
    rr = rng.uniform(0.3, 3.5, gl.size)
    fast = _disc_mask(ll, bb, gl, gb, rr)
    slow = _disc_mask_bruteforce(ll, bb, gl, gb, rr)
    ndiff = int((fast != slow).sum())
    print(f"  fast vs brute-force disc mask: {ndiff} differing pixels")
    if ndiff:
        print("FAIL fast mask disagrees with brute force"); ok = False

    # mask D: radius must rise with source flux, invert the PSF exactly, and
    # keep more sky than the fixed-radius mask A
    E, B = PSF95_DEG[1] and 2.0, 1e-5
    fl = np.array([1e-7, 1e-8, 1e-9, 1e-11])
    rd = mask_e_radii(fl, E, B)
    print("  mask D radii for flux 1e-7..1e-11: " + ", ".join(f"{x:.2f}" for x in rd))
    if not np.all(np.diff(rd) <= 0):
        print("FAIL mask D radius should not grow as the source gets fainter"); ok = False
    for f_, r_ in zip(fl[:2], rd[:2]):
        lhs = f_ * float(PSF.king(np.array([r_]), E)[0])
        if abs(lhs / (MASK_E_FRACTION * B) - 1) > 1e-3:
            print(f"FAIL mask D radius does not invert the PSF ({lhs:.3e} vs {B:.3e})")
            ok = False
    d_test = _disc_mask(ll, bb, glon, glat, mask_e_radii(flux, E, B))
    print(f"  mask D keeps {d_test.mean():6.1%} of this synthetic ROI "
          f"(A keeps {a.mean():.1%})")
    if d_test.mean() < a.mean():
        print("FAIL mask D should keep at least as much sky as mask A"); ok = False

    print(f"  mask A keeps {a.mean():6.1%} of the ROI")
    print(f"  mask B keeps {b.mean():6.1%} of the ROI")
    print(f"  mask C keeps {c.mean():6.1%} of the ROI")
    print(f"  radii by energy bin: " +
          ", ".join(f"{EBINS[i][0]}-{EBINS[i][1]}GeV {PSF95_DEG[i]+0.5:.2f}deg"
                    for i in range(4)))
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        selftest()
    else:
        cat = sys.argv[1]
        for w in "ABC":
            print(write(build(cat, w), f"mask_{w.lower()}_wmap.fits", w, cat))
