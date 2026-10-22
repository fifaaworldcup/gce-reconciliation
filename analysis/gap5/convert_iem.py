"""
galprop maps -> one cube per component (pi0, bremss, ics) on our grid.
reads CAR mapcubes (incl. 0..360 maps with the GC on the seam) and HEALPix
tables. output cubes have an ENERGIES table in MeV, same as read_iem expects.
if a total map and per-ring maps are both there, the total is used.

  python analysis/gap5/convert_iem.py --files <one model's files> --out <dir>
  python analysis/gap5/convert_iem.py --selftest
"""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

NPIX, PIX, HALF = 400, 0.1, 20.0

# first match wins
COMPONENTS = [
    ("pi0", re.compile(r"^pion[_-]?decay|^pi0_", re.I)),
    ("bremss", re.compile(r"^brem", re.I)),
    ("ics", re.compile(r"^ics|^inverse[_-]?compton", re.I)),
    ("synch", re.compile(r"^synch", re.I)),
]
# ring / gas phase parts, only used if there's no total
RING = re.compile(r"(HIR|H2R|HII|_ring_?\d|_r\d+_|isotropic_comp|anisotropic_comp"
                  r"|_optical|_infrared|_cmb)", re.I)


def our_grid(oversample=1):
    k = oversample
    step = PIX / k
    l = (HALF - step / 2.0) - step * np.arange(NPIX * k)
    b = (-HALF + step / 2.0) + step * np.arange(NPIX * k)
    return np.meshgrid(l, b)


def _block_mean(a, k):
    if k == 1:
        return a
    return a.reshape(NPIX, k, NPIX, k).mean(axis=(1, 3))


# ------------------------------------------------------------------ readers
def _energies_from_table(hdul):
    for h in hdul[1:]:
        if h.name.upper() in ("ENERGIES", "EBOUNDS", "ENERGY"):
            d = h.data
            col = d.names[0]
            e = np.asarray(d[col], float)
            unit = (h.columns[col].unit or "MeV").strip().lower()
            if unit == "gev":
                e = e * 1e3
            elif unit == "kev":
                e = e / 1e3
            return e
    return None


def _energies_from_axis(hdr, n):
    # old mapcubes keep log10(E/MeV) on axis 3
    if "CRVAL3" not in hdr:
        return None
    c0, d, p = float(hdr["CRVAL3"]), float(hdr.get("CDELT3", 0)), float(hdr.get("CRPIX3", 1))
    v = c0 + (np.arange(n) + 1 - p) * d
    ctype = str(hdr.get("CTYPE3", "")).lower()
    if "log" in ctype or (0 < c0 < 8 and 0 < d < 1):
        return 10.0 ** v
    return v


def _cholis_wcs(hdr, shape):
    # the zenodo 6423495 maps carry no WCS, only EMIN/EMAX. 240x240 over
    # |l|, |b| < 30 at 0.25 deg; l rises with column, b with row, not
    # transposed. checked against gll_iem_v07 at 0.9 GeV: corr 0.991 for this
    # layout, 0.907 with l flipped, < 0.1 transposed (see check_orientation)
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"], h["CRVAL1"], h["CDELT1"] = (shape[1] + 1) / 2.0, 0.0, 60.0 / shape[1]
    h["CRPIX2"], h["CRVAL2"], h["CDELT2"] = (shape[0] + 1) / 2.0, 0.0, 60.0 / shape[0]
    return WCS(h)


def read_galprop(path):
    with fits.open(path) as hd:
        prim = hd[0]
        if prim.data is not None and np.asarray(prim.data).ndim >= 3:
            a = np.asarray(prim.data, float)
            while a.ndim > 3:
                if a.shape[0] != 1:
                    raise SystemExit(f"{path}: {a.ndim}D cube with leading axis "
                                     f"{a.shape[0]}, not a single model")
                a = a[0]
            hdr = prim.header
            e = _energies_from_table(hd)
            if e is None:
                e = _energies_from_axis(hdr, a.shape[0])
            if e is None and "EMIN" in hdr and "EMAX" in hdr:
                # cholis maps: planes log-spaced from EMIN to EMAX, MeV
                e = float(hdr["EMIN"]) * (float(hdr["EMAX"]) / float(hdr["EMIN"])) ** (
                    np.arange(a.shape[0]) / (a.shape[0] - 1))
            if e is None or len(e) != a.shape[0]:
                raise SystemExit(f"{path}: cannot place {a.shape[0]} planes in energy")
            if "CTYPE1" in hdr:
                w = WCS(hdr).celestial
            elif "InnerGalaxy_60x60" in os.path.basename(path):
                w = _cholis_wcs(hdr, a.shape[1:])
            else:
                raise SystemExit(f"{path}: no WCS and not a known layout")
            return a, e, _car_sampler(w, a.shape[1:])
        for h in hd[1:]:
            hdr = h.header
            if str(hdr.get("PIXTYPE", "")).upper() == "HEALPIX" or "NSIDE" in hdr:
                nside = int(hdr["NSIDE"])
                nest = str(hdr.get("ORDERING", "RING")).upper().startswith("NEST")
                cols = h.data.names
                if len(cols) == 1 and np.asarray(h.data[cols[0]]).ndim == 2:
                    a = np.asarray(h.data[cols[0]], float).T          # (nE, npix)
                else:
                    bins = [c for c in cols if re.match(r"^bin\d+$", c, re.I)]
                    use = sorted(bins, key=lambda c: int(c[3:])) if bins else cols
                    a = np.stack([np.asarray(h.data[c], float).ravel() for c in use])
                e = _energies_from_table(hd)
                if e is None or len(e) != a.shape[0]:
                    raise SystemExit(f"{path}: HEALPix with {a.shape[0]} planes "
                                     f"and no matching ENERGIES table")
                return a, e, _healpix_sampler(nside, nest)
    raise SystemExit(f"{path}: neither a CAR cube nor a HEALPix table")


def _car_sampler(w, shape):
    from scipy.ndimage import map_coordinates
    ny, nx = shape

    # astropy folds longitude itself, so wrap 0..360 maps in pixels
    per360 = 360.0 / abs(float(np.diag(w.pixel_scale_matrix)[0]))

    def sample(plane, l, b):
        best = np.full(l.shape, np.nan)
        x0, y = w.wcs_world2pix(l, b, 0)
        for shift in (0.0, per360, -per360):
            x = x0 + shift
            ok = (x >= -0.5) & (x <= nx - 0.5) & (y >= -0.5) & (y <= ny - 0.5)
            ok &= np.isnan(best)
            if ok.any():
                xv = np.clip(x[ok], 0, nx - 1)
                yv = np.clip(y[ok], 0, ny - 1)
                best[ok] = map_coordinates(plane, np.stack([yv, xv]), order=1,
                                           mode="nearest")
        return best
    return sample


def _healpix_sampler(nside, nest):
    import healpy as hp

    def sample(plane, l, b):
        return hp.get_interp_val(plane, np.mod(l, 360.0), b, nest=nest,
                                 lonlat=True)
    return sample


def to_grid(path, oversample=3):
    planes, e, sampler = read_galprop(path)
    L, B = our_grid(oversample)
    out = np.empty((len(e), NPIX, NPIX))
    for k in range(len(e)):
        v = sampler(planes[k], L, B)
        if np.isnan(v).any():
            frac = float(np.isnan(v).mean())
            raise SystemExit(f"{path}: plane {k} does not cover the ROI "
                             f"({frac:.1%} of samples outside the map)")
        out[k] = _block_mean(v, oversample)
    return out, e


# ------------------------------------------------------------- components
def component_of(path):
    name = os.path.basename(path)
    for comp, pat in COMPONENTS:
        if pat.search(name):
            return comp, bool(RING.search(name))
    return None, False


def plan(files):
    by = {}
    for f in files:
        c, is_part = component_of(f)
        if c is None:
            print(f"  skipped (not an emission component): {os.path.basename(f)}")
            continue
        by.setdefault(c, {"total": [], "parts": []})
        by[c]["total" if not is_part else "parts"].append(f)
    out = {}
    for c, d in by.items():
        if len(d["total"]) > 1:
            raise SystemExit(f"{c}: {len(d['total'])} total maps "
                             f"{[os.path.basename(x) for x in d['total']]} - "
                             "pass the files for ONE model")
        out[c] = d["total"] if d["total"] else d["parts"]
    return out


def write_cube(cube, e_mev, path, meta):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRPIX1"] = h["CRPIX2"] = (NPIX + 1) / 2.0
    h["CRVAL1"] = h["CRVAL2"] = 0.0
    h["CDELT1"], h["CDELT2"] = -PIX, PIX
    h["EQUINOX"] = 2000.0
    h["DATE"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    for k, v in meta.items():
        h[k] = v
    en = fits.BinTableHDU.from_columns(
        [fits.Column(name="Energy", format="D", unit="MeV", array=np.asarray(e_mev))],
        name="ENERGIES")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(cube.astype(np.float32), h), en]).writeto(
        path, overwrite=True)


def check_orientation(path, ref_cube, e_mev=1000.0, min_corr=0.95):
    # correlate one plane with a reference diffuse model on our grid, in log,
    # away from the plane. a flipped or transposed map fails this
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gap1"))
    import fit_cubes as FC
    planes, e, sampler = read_galprop(path)
    k = int(np.argmin(np.abs(np.log(e / e_mev))))
    L, B = our_grid(1)
    v = sampler(planes[k], L, B)
    ref, _ = FC.read_iem(ref_cube, np.array([e[k] / 1e3]))
    m = (np.abs(B) > 0.5) & (v > 0) & (ref[0] > 0)
    c = float(np.corrcoef(np.log(v[m]), np.log(ref[0][m]))[0, 1])
    return c, c >= min_corr


def convert(files, out, model="", oversample=3):
    todo = plan(files)
    if not todo:
        raise SystemExit("no pion_decay / bremss / ics maps among the inputs")
    report = {}
    for comp, fs in sorted(todo.items()):
        cube, e = None, None
        for f in fs:
            c, ef = to_grid(f, oversample)
            if e is None:
                cube, e = c, ef
            else:
                if len(ef) != len(e) or not np.allclose(ef, e, rtol=1e-4):
                    raise SystemExit(f"{comp}: parts have different energy grids")
                cube = cube + c
        if (cube < 0).any():
            neg = float((cube < 0).mean())
            print(f"  {comp}: {neg:.2%} negative pixels clipped to zero")
            cube = np.maximum(cube, 0.0)
        path = os.path.join(out, f"{comp}.fits")
        write_cube(cube, e, path, {
            "MODEL": (model[:60], "Cholis et al. 2022 model id"),
            "COMP": comp,
            "NSRC": (len(fs), "source files summed"),
            "SRCFILE": os.path.basename(fs[0])[:60],
            "EMINMEV": float(e[0]), "EMAXMEV": float(e[-1])})
        report[comp] = {"files": [os.path.basename(f) for f in fs],
                        "planes": int(len(e)),
                        "energy_mev": [float(e[0]), float(e[-1])]}
        print(f"  {comp:7s} <- {len(fs)} file(s), {len(e)} planes "
              f"{e[0]:.3g}-{e[-1]:.3g} MeV -> {path}")
    return report


# --------------------------------------------------------------- self-test
def _truth(l, b, e_mev):
    lw = np.where(l > 180, l - 360, l)
    shape = np.exp(-(lw - 2.0) ** 2 / 60.0 - (b + 0.5) ** 2 / 9.0) + 0.05
    return shape * (e_mev / 1e3) ** (-2.0 - 0.02 * lw / 20.0)


def selftest():
    import tempfile
    import healpy as hp
    ok = True
    e = np.logspace(2, 5, 13)                         # 100 MeV - 100 GeV
    L, B = our_grid(3)
    want = np.stack([_block_mean(_truth(L, B, x), 3) for x in e])
    core = want[0] > 0.2 * want[0].max()

    def check(tag, got, eg):
        nonlocal ok
        if len(eg) != len(e) or not np.allclose(eg, e, rtol=1e-6):
            print(f"  FAIL {tag}: energies read back wrong {eg[:3]}"); ok = False
            return
        err = np.max(np.abs(got / want - 1)[:, core])
        print(f"  {tag:34s} max error in the bright region {err:.3%}")
        if err > 0.02:
            print(f"  FAIL {tag}"); ok = False

    with tempfile.TemporaryDirectory() as tmp:
        # healpix, one column per energy
        nside = 512
        th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
        cols = [fits.Column(name=f"Bin{k}", format="D", array=_truth(th, ph, x))
                for k, x in enumerate(e)]
        t = fits.BinTableHDU.from_columns(cols)
        t.header["PIXTYPE"], t.header["NSIDE"], t.header["ORDERING"] = "HEALPIX", nside, "RING"
        en = fits.BinTableHDU.from_columns(
            [fits.Column(name="Energy", format="D", unit="MeV", array=e)], name="ENERGIES")
        p1 = os.path.join(tmp, "pion_decay_healpix_test.fits")
        fits.HDUList([fits.PrimaryHDU(), t, en]).writeto(p1)
        check("HEALPix RING, Bin columns", *to_grid(p1))

        # CAR 0..360, energy on a log axis
        dl = 0.125
        nl, nb = int(360 / dl), int(180 / dl)
        lc = (np.arange(nl)) * dl                     # pixel 1 at l = 0
        bc = -90 + dl / 2 + np.arange(nb) * dl
        LL, BB = np.meshgrid(lc, bc)
        cube = np.stack([_truth(LL, BB, x) for x in e])
        h = fits.Header()
        h["CTYPE1"], h["CTYPE2"], h["CTYPE3"] = "GLON-CAR", "GLAT-CAR", "log_Energy"
        h["CRPIX1"], h["CRVAL1"], h["CDELT1"] = 1.0, 0.0, dl
        h["CRPIX2"], h["CRVAL2"], h["CDELT2"] = (nb + 1) / 2.0, 0.0, dl
        h["CRPIX3"], h["CRVAL3"], h["CDELT3"] = 1.0, 2.0, 0.25
        p2 = os.path.join(tmp, "bremss_mapcube_test.fits")
        fits.PrimaryHDU(cube, h).writeto(p2)
        check("CAR 0..360, log10 energy axis", *to_grid(p2))

        # totals win over ring parts
        fs = ["pion_decay_mapcube_54_X.gz", "pion_decay_HIR_ring_1_mapcube_54_X.gz",
              "bremss_mapcube_54_X.gz", "ics_isotropic_mapcube_54_X.gz",
              "ics_isotropic_comp_1_mapcube_54_X.gz", "nuclei_full_54_X.gz"]
        pl = plan(fs)
        want_plan = {"pi0": [fs[0]], "bremss": [fs[2]], "ics": [fs[3]]}
        print(f"  component plan: { {k: [os.path.basename(x) for x in v] for k, v in pl.items()} }")
        if pl != want_plan:
            print(f"  FAIL plan, wanted {want_plan}"); ok = False

        # read back the way the fit reads it
        here = os.path.dirname(os.path.abspath(__file__))
        root = os.path.normpath(os.path.join(here, "..", ".."))
        for p in (os.path.join(root, "analysis", "gap1"), os.path.join(root, "utils"),
                  os.path.join(root, "masks")):
            if p not in sys.path:
                sys.path.insert(0, p)
        import fit_cubes as FC
        rep = convert([p1, p2], os.path.join(tmp, "out"), model="selftest")
        emid = np.array([1.155, 3.16, 7.5])
        back, info = FC.read_iem(os.path.join(tmp, "out", "pi0.fits"), emid)
        direct = np.stack([_block_mean(_truth(L, B, x * 1e3), 3) for x in emid])
        err = np.max(np.abs(back / direct - 1)[:, core])
        print(f"  through fit_cubes.read_iem at bin centres: max error {err:.3%}, "
              f"coverage {info['coverage_min']:.1%}")
        if err > 0.03 or info["coverage_min"] < 0.999:
            print("  FAIL round trip"); ok = False
        if sorted(rep) != ["bremss", "pi0"]:
            print(f"  FAIL wrote {sorted(rep)}"); ok = False

    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="*", default=[])
    ap.add_argument("--glob", default="", help="alternative to --files")
    ap.add_argument("--out", default="")
    ap.add_argument("--model", default="")
    ap.add_argument("--oversample", type=int, default=3)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--cholis", default="", help="ZENODO_FILES folder; with --models")
    ap.add_argument("--models", default="", help="roman numerals, comma separated")
    ap.add_argument("--check-against", default="", help="reference cube for the orientation check")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.cholis:
        names = {}
        for line in open(os.path.join(a.cholis, "NAMING_CONVENTION_OF_DIFFUSE_EMISSION_MODELS.dat")):
            if line.strip() and not line.startswith("#"):
                r, code, _ = line.split()
                names[r] = code
        mdir = os.path.join(a.cholis, "GALACTIC_DIFFUSE_EMISSION_MAPS_0p25deg")
        for r in [m for m in a.models.split(",") if m]:
            fs = [os.path.join(mdir, f"{c}_{names[r]}_Map_flux_E_50-814008_MeV_InnerGalaxy_60x60.fits")
                  for c in ("pi0", "bremss", "ICS")]
            miss = [os.path.basename(f) for f in fs if not os.path.exists(f)]
            if miss:
                print(f"{r}: missing {', '.join(miss)} - skipped")
                continue
            if a.check_against:
                c, ok = check_orientation(fs[0], a.check_against)
                print(f"{r}: orientation check corr {c:.3f}")
                if not ok:
                    raise SystemExit(f"{r}: pi0 map does not line up with the reference")
            print(f"{r} (code {names[r]}):")
            convert(fs, os.path.join(a.out or ".", f"cholis_{r}"), model=r, oversample=a.oversample)
        return
    files = list(a.files) + (sorted(glob.glob(a.glob)) if a.glob else [])
    if not files or not a.out:
        ap.error("need --files/--glob and --out")
    convert(files, a.out, a.model, a.oversample)


if __name__ == "__main__":
    main()
