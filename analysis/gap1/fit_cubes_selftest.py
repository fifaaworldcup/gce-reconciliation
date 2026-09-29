"""
End-to-end test of fit_cubes.py against mock cubes in the exact FITS layout
gtbin / gtexpcube2 / reduce_iem.py emit.

The point is not to check the maths - fit_likelihood.py has its own closed-loop
test for that. The point is the plumbing: keV vs MeV vs GeV, bin edges vs bin
centres, the 440 -> 400 exposure crop, a diffuse model on a different pixel
size with the opposite CDELT1 sign, and the catalogue column spellings. Every
one of those is a silent wrong-answer bug rather than a crash, so they need a
test where the right answer is known.

Injects a known signal, then checks the fit finds it, recovers its
normalisation, and prefers the template that was actually injected.
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (HERE, os.path.join(ROOT, "utils")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_likelihood as FL      # noqa: E402
import regrid as RG              # noqa: E402

NPIX, PIX = 400, 0.1
NE = 8
EMIN_GEV, EMAX_GEV = 1.0, 10.0
EXP_NPIX = 440               # gtexpcube2 is run wider than the ROI on purpose
IEM_PIX = 0.125              # the diffuse model's native pixel size
IEM_HALF = 21.0


def _edges():
    return np.logspace(np.log10(EMIN_GEV), np.log10(EMAX_GEV), NE + 1)


def write_ccube(path, counts):
    """gtbin CCUBE layout: 3D primary + EBOUNDS in keV."""
    ed = _edges()
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (NPIX + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = -PIX, PIX
    h["CTYPE3"] = "ENERGY"
    pri = fits.PrimaryHDU(counts.astype(np.float32), h)
    eb = fits.BinTableHDU.from_columns([
        fits.Column(name="CHANNEL", format="J", array=np.arange(1, NE + 1)),
        fits.Column(name="E_MIN", format="E", unit="keV",
                    array=(ed[:-1] * 1e6).astype(np.float32)),
        fits.Column(name="E_MAX", format="E", unit="keV",
                    array=(ed[1:] * 1e6).astype(np.float32)),
    ], name="EBOUNDS")
    fits.HDUList([pri, eb]).writeto(path, overwrite=True)


def write_expcube(path, exp_planes, plane_energies_gev):
    """gtexpcube2 layout: wider grid, energies at bin BOUNDARIES, in MeV."""
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (EXP_NPIX + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = -PIX, PIX
    pri = fits.PrimaryHDU(exp_planes.astype(np.float32), h)
    en = fits.BinTableHDU.from_columns([
        fits.Column(name="Energy", format="D", unit="MeV",
                    array=np.asarray(plane_energies_gev) * 1e3),
    ], name="ENERGIES")
    fits.HDUList([pri, en]).writeto(path, overwrite=True)


def write_iem(path, cube, energies_gev):
    """reduce_iem.py layout: coarser pixels, CDELT1 POSITIVE (the all-sky
    diffuse model runs l increasing with x, opposite to our grid). If the
    regridder silently ignored that the map would come out mirrored."""
    n = cube.shape[-1]
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (n + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = +IEM_PIX, IEM_PIX
    pri = fits.PrimaryHDU(cube.astype(np.float32), h)
    en = fits.BinTableHDU.from_columns([
        fits.Column(name="Energy", format="D", unit="MeV",
                    array=np.asarray(energies_gev) * 1e3),
    ], name="ENERGIES")
    fits.HDUList([pri, en]).writeto(path, overwrite=True)


def write_catalog(path, glon, glat, sigma, flux):
    """4FGL column spellings: GLON / GLAT / Signif_Avg / Flux1000."""
    t = fits.BinTableHDU.from_columns([
        fits.Column(name="GLON", format="E", array=np.asarray(glon, np.float32)),
        fits.Column(name="GLAT", format="E", array=np.asarray(glat, np.float32)),
        fits.Column(name="Signif_Avg", format="E", array=np.asarray(sigma, np.float32)),
        fits.Column(name="Flux1000", format="E", array=np.asarray(flux, np.float32)),
    ], name="LAT_Point_Source_Catalog")
    fits.HDUList([fits.PrimaryHDU(), t]).writeto(path, overwrite=True)


def build_mock(tmp, signal_template, rng):
    ed = _edges()
    emid = np.sqrt(ed[:-1] * ed[1:])

    l = (20 - PIX / 2) - PIX * np.arange(NPIX)
    b = (-20 + PIX / 2) + PIX * np.arange(NPIX)
    ll, bb = np.meshgrid(l, b)

    # ---- diffuse model on its own coarse grid, l increasing with x
    n = int(round(2 * IEM_HALF / IEM_PIX)) + 1
    li = (-IEM_HALF) + IEM_PIX * np.arange(n)      # matches CDELT1 > 0
    bi = (-IEM_HALF) + IEM_PIX * np.arange(n)
    lli, bbi = np.meshgrid(li, bi)
    iem_e = np.logspace(np.log10(0.5), np.log10(50.0), 12)
    # a plane that narrows with energy, so an energy mix-up shows up
    iem_cube = np.stack([
        (np.exp(-(bbi / (4.0 * (E / 3.0) ** -0.25)) ** 2) * (1 + 0.3 * np.cos(np.radians(lli)))
         ) * E ** -2.7
        for E in iem_e])
    write_iem(os.path.join(tmp, "iem.fits"), iem_cube, iem_e)

    # the truth on OUR grid: same function, evaluated directly
    iem_true = np.stack([
        RG.normalise(np.exp(-(bb / (4.0 * (E / 3.0) ** -0.25)) ** 2)
                     * (1 + 0.3 * np.cos(np.radians(ll))))
        for E in emid])

    # ---- exposure: wider grid, mild gradient, falling with energy
    le = (22 - PIX / 2) - PIX * np.arange(EXP_NPIX)
    be = (-22 + PIX / 2) + PIX * np.arange(EXP_NPIX)
    lle, bbe = np.meshgrid(le, be)
    the = np.hypot(lle, bbe)
    exp_planes = np.stack([5.0e10 * (1 - 0.12 * the / 31.0) * (E / 3.0) ** -0.1
                           for E in ed])          # at BOUNDARIES, NE+1 planes
    write_expcube(os.path.join(tmp, "exp.fits"), exp_planes, ed)

    th = np.hypot(ll, bb)
    exp_true = np.stack([5.0e10 * (1 - 0.12 * th / 31.0) * (E / 3.0) ** -0.1
                         for E in emid])

    # ---- components and the injected truth
    flat = RG.normalise(np.ones((NPIX, NPIX)))
    sig = RG.normalise(signal_template)
    want = {"iem": 6.0e5, "iso": 1.2e5, "sig": 8.0e4}       # counts per bin

    truth = {}
    for k, tm in (("iem", iem_true[0]), ("iso", flat), ("sig", sig)):
        conv = FL.convolve_templates({k: tm}, emid[0])[k]
        truth[k] = want[k] / float((exp_true[0] * conv).sum())

    counts = np.zeros((NE, NPIX, NPIX))
    for e in range(NE):
        comp = {"iem": iem_true[e], "iso": flat, "sig": sig}
        conv = FL.convolve_templates(comp, emid[e])
        mu = exp_true[e] * sum(truth[k] * conv[k] for k in comp)
        counts[e] = rng.poisson(mu)
    write_ccube(os.path.join(tmp, "ccube.fits"), counts)

    # ---- catalogue: a handful of sources to mask
    ns = 40
    write_catalog(os.path.join(tmp, "cat.fits"),
                  rng.uniform(-19, 19, ns), rng.uniform(-19, 19, ns),
                  rng.uniform(4, 40, ns), rng.uniform(1e-10, 5e-8, ns))

    return counts, truth, exp_true, emid


def main():
    ok = True
    rng = np.random.default_rng(7)
    tmp = tempfile.mkdtemp(prefix="gap1fit_")
    try:
        tdir = os.path.join(ROOT, "templates")
        # inject the HESTIA template, then ask the fit to choose between it
        # and gNFW. it has to pick the one that is actually there.
        inj_path = os.path.join(tdir, "gap1_hestia_G1.1_ang0_1deg.fits")
        inj, _ = RG.load(inj_path)

        counts, truth, exp_true, emid = build_mock(tmp, inj, rng)
        print(f"  mock sky: {counts.sum():,.0f} photons in {NE} bins")

        sys.argv = ["fit_cubes.py",
                    "--ccube", os.path.join(tmp, "ccube.fits"),
                    "--expcube", os.path.join(tmp, "exp.fits"),
                    "--iem", os.path.join(tmp, "iem.fits"),
                    "--catalog", os.path.join(tmp, "cat.fits"),
                    "--signals", "hestia,gnfw",
                    "--out", os.path.join(tmp, "out")]
        import fit_cubes
        rec = fit_cubes.main()

        print()
        # 1. the exposure crop must have been found, not guessed
        off = rec["inputs"]["exposure"]["crop_offset_px"]
        want_off = (EXP_NPIX - NPIX) // 2
        if off != [want_off, want_off]:
            print(f"  FAIL exposure crop offset {off}, expected "
                  f"[{want_off}, {want_off}]"); ok = False
        else:
            print(f"  exposure crop offset {off} ok")

        # 2. energies read back as GeV, not MeV or keV
        bins = rec["inputs"]["bins_gev"]
        if not (abs(bins[0][0] - EMIN_GEV) < 1e-3 and
                abs(bins[-1][1] - EMAX_GEV) < 1e-2):
            print(f"  FAIL energy bins came back as {bins[0]}..{bins[-1]}, "
                  f"expected {EMIN_GEV}..{EMAX_GEV} GeV"); ok = False
        else:
            print(f"  energy bins {bins[0][0]:.2f}-{bins[-1][1]:.2f} GeV ok")

        # 3. the injected template must win
        cmp_key = [k for k in rec["comparisons"]][0]
        fav = rec["comparisons"][cmp_key]["favoured"]
        d = rec["comparisons"][cmp_key]["delta_lnL"]
        print(f"  {cmp_key}: dlnL={d:+.1f}, favours {fav}")
        if fav != "hestia":
            print(f"  FAIL injected hestia but the fit preferred {fav}")
            ok = False

        # 4. its normalisation must come back right
        rec_norm = np.mean([b["hestia"] for b in
                            rec["models"]["hestia"]["norms_per_bin"]])
        ratio = rec_norm / truth["sig"]
        print(f"  hestia normalisation: injected {truth['sig']:.3e}, "
              f"recovered {rec_norm:.3e}, ratio {ratio:.3f}")
        if not 0.90 < ratio < 1.10:
            print("  FAIL normalisation off by more than 10%"); ok = False

        # 5. TS must be large, and the residuals must look like noise
        ts = rec["models"]["hestia"]["TS"]
        rs = rec["models"]["hestia"]["residual_std_sigma"]
        rm = rec["models"]["hestia"]["residual_mean_sigma"]
        print(f"  hestia TS {ts:,.0f}, residual {rm:+.3f} +- {rs:.3f} sigma")
        if ts < 100:
            print("  FAIL injected signal should be a strong detection"); ok = False
        if not 0.8 < rs < 1.25:
            print("  FAIL residual scatter is not ~1 sigma"); ok = False

        # 6. the mask was actually applied
        kf = rec["inputs"]["mask"]["kept_fraction"]
        print(f"  mask A kept {kf:.1%} of the ROI")
        if not 0.5 < kf < 0.999:
            print("  FAIL mask kept a suspicious fraction"); ok = False

        if not rec["models"]["hestia"]["all_bins_converged"]:
            print("  FAIL a bin did not converge"); ok = False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
