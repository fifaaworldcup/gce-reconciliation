"""
Point-source template from the 4FGL catalogue.

Why: the Mask B residual map is covered in unmasked faint sources (~720 of
them in the ROI). Nothing in the model accounts for them, and an extended GCE
template can soak up the ones near the Galactic Centre in a way a narrow cusp
cannot - so a template comparison without them can reward extension for the
wrong reason. That is the classic GCE confusion and exactly what Gap 1 has to
be careful about.

So: every catalogue source dropped at its position, weighted by its own
spectral model integrated over each of our energy bins, PSF-convolved by the
fitter like any other template, with one free normalisation per bin.

Kept in physical units (ph cm^-2 s^-1 per pixel), NOT normalised, so the
fitted normalisation should come out near 1 if the catalogue, our exposure
cube and our PSF agree. That is a free end-to-end check on the exposure.

Spectral forms are the 4FGL-DR4 ones:
  PowerLaw           K (E/E0)^-G
  LogParabola        K (E/E0)^-(a + b ln(E/E0))
  PLSuperExpCutoff4  K (E/E0)^(d/b - Gs) exp(d/b^2 (1 - (E/E0)^b)),  b != 0
                     K (E/E0)^(-Gs - d/2 ln(E/E0)),                  b == 0
"""

import numpy as np
from astropy.io import fits

NPIX, PIX, HALF = 400, 0.1, 20.0


def dnde(t, E_mev):
    """Differential flux of every source at energies E (MeV). (nsrc, nE)."""
    E = np.atleast_1d(np.asarray(E_mev, float))[None, :]
    typ = np.char.strip(np.asarray(t["SpectrumType"]).astype(str))
    E0 = np.asarray(t["Pivot_Energy"], float)[:, None]
    x = E / E0
    lx = np.log(x)
    out = np.zeros((len(typ), E.shape[1]))

    m = typ == "PowerLaw"
    out[m] = (np.asarray(t["PL_Flux_Density"], float)[m, None]
              * x[m] ** -np.asarray(t["PL_Index"], float)[m, None])

    m = typ == "LogParabola"
    a = np.asarray(t["LP_Index"], float)[m, None]
    b = np.asarray(t["LP_beta"], float)[m, None]
    out[m] = (np.asarray(t["LP_Flux_Density"], float)[m, None]
              * x[m] ** -(a + b * lx[m]))

    m = np.char.startswith(typ, "PLSuperExpCutoff")
    K = np.asarray(t["PLEC_Flux_Density"], float)[m, None]
    Gs = np.asarray(t["PLEC_IndexS"], float)[m, None]
    d = np.asarray(t["PLEC_ExpfactorS"], float)[m, None]
    bb = np.asarray(t["PLEC_Exp_Index"], float)[m, None]
    xm, lxm = x[m], lx[m]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        nz = np.abs(bb) > 1e-6
        f_nz = K * xm ** (d / bb - Gs) * np.exp(d / bb ** 2 * (1 - xm ** bb))
        f_z = K * xm ** (-Gs - d / 2.0 * lxm)
        out[m] = np.where(nz, f_nz, f_z)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def band_flux(t, lo_mev, hi_mev, n=33):
    """Integrated photon flux in [lo, hi] for every source, log-trapezoid."""
    E = np.logspace(np.log10(lo_mev), np.log10(hi_mev), n)
    f = dnde(t, E)
    return np.trapezoid(f * E[None, :], np.log(E), axis=1)


def build(catalog_path, e_lo_gev, e_hi_gev):
    """(nE, 400, 400) point-source cube on the conventions grid."""
    t = fits.getdata(catalog_path, 1)
    glon = np.asarray(t["GLON"], float)
    glon = np.where(glon > 180, glon - 360, glon)
    glat = np.asarray(t["GLAT"], float)
    j = np.floor((HALF - glon) / PIX).astype(int)          # l descends
    i = np.floor((glat + HALF) / PIX).astype(int)
    inside = (i >= 0) & (i < NPIX) & (j >= 0) & (j < NPIX)
    cube = np.zeros((len(e_lo_gev), NPIX, NPIX))
    for k, (lo, hi) in enumerate(zip(e_lo_gev, e_hi_gev)):
        fl = band_flux(t, lo * 1e3, hi * 1e3)
        np.add.at(cube[k], (i[inside], j[inside]), fl[inside])
    return cube, {"sources_placed": int(inside.sum()),
                  "note": "extended 4FGL sources placed as points"}


def selftest(catalog_path):
    """Our spectral integration from 1-100 GeV must reproduce the catalogue's
    own Flux1000 column - same models, same integral, so any disagreement
    is a bug in the spectral forms above."""
    t = fits.getdata(catalog_path, 1)
    ours = band_flux(t, 1e3, 1e5, n=129)
    cat = np.asarray(t["Flux1000"], float)
    ok_rows = (cat > 0) & np.isfinite(ours)
    ratio = ours[ok_rows] / cat[ok_rows]
    typ = np.char.strip(np.asarray(t["SpectrumType"]).astype(str))[ok_rows]
    good = True
    for ty in sorted(set(typ)):
        r = ratio[typ == ty]
        print(f"  {ty:18s} n={r.size:5d}  median ratio {np.median(r):.4f}  "
              f"worst {np.max(np.abs(r - 1)):.2%}")
        if abs(np.median(r) - 1) > 0.01:
            good = False
    print("SELFTEST PASS" if good else "SELFTEST FAIL")
    return good


if __name__ == "__main__":
    import sys
    selftest(sys.argv[1])
