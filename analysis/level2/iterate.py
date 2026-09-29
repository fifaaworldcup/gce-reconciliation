"""
Iterative source-finding (Level 2), shared by Gap 1 and Gap 5.

Fit -> residual TS map -> add the best new point source above threshold ->
refit -> repeat, in successive passes at decreasing threshold, exactly as
Di Mauro (2026) describes. The stopping rule is read from configs/level2.yaml
and nothing here lets a caller override it per template: that is the point.

Residual TS map. For a point source at pixel p, with A(p) the PSF-convolved
exposure-weighted unit-flux map of that source, the score test for adding it
to the current fit is

    TS(p) = S(p)^2 / F(p)  for S > 0,
    S(p)  = sum_i A_i(p) (n_i/mu_i - 1),
    F(p)  = sum_i A_i(p)^2 n_i / mu_i^2

Because A_i(p) = exposure_i * PSF(|i - p|), both sums are cross-correlations
with the PSF kernel, so a whole map costs two convolutions per energy bin
instead of one fit per pixel.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for p in (os.path.join(ROOT, "analysis", "gap1"), os.path.join(ROOT, "utils")):
    if p not in sys.path:
        sys.path.insert(0, p)

import fit_likelihood as FL      # noqa: E402
import psf as PSF                # noqa: E402

NPIX, PIX, HALF = 400, 0.1, 20.0
CONFIG = os.path.join(ROOT, "configs", "level2.yaml")


def load_config(path=CONFIG):
    """Read the pre-registered stopping rule. No defaults in code: if the file
    is missing we stop, rather than quietly inventing thresholds."""
    import yaml
    with open(path) as f:
        c = yaml.safe_load(f)
    need = ("passes", "max_sources_per_pass", "min_separation_deg",
            "edge_margin_deg", "seed")
    missing = [k for k in need if k not in c]
    if missing:
        raise SystemExit(f"{path} is missing {missing}")
    return c


def _grid():
    l = (HALF - PIX / 2.0) - PIX * np.arange(NPIX)
    b = (-HALF + PIX / 2.0) + PIX * np.arange(NPIX)
    return np.meshgrid(l, b)


def point_source_map(i, j):
    """Unit-flux point source at pixel (i, j), before PSF convolution."""
    t = np.zeros((NPIX, NPIX))
    t[i, j] = 1.0
    return t


def residual_ts_map(counts, exposure, comps, energies, result, mask=None):
    """Score-test TS for adding one point source at every pixel."""
    from scipy.signal import fftconvolve
    S = np.zeros((NPIX, NPIX))
    F = np.zeros((NPIX, NPIX))
    for e in range(counts.shape[0]):
        tb = FL.bin_templates(comps, e)
        conv = FL.convolve_templates(tb, energies[e])
        norms = result["bins"][e]["norms"]
        mu = exposure[e] * sum(norms[k] * conv[k] for k in tb)
        mu = np.maximum(mu, 1e-12)
        m = (np.ones(mu.shape, bool) if mask is None
             else mask[e].astype(bool)) & np.isfinite(exposure[e])
        k = PSF.kernel(energies[e])
        # S: correlate (n/mu - 1) * exposure with the PSF
        a = np.where(m, (counts[e] / mu - 1.0) * exposure[e], 0.0)
        S += fftconvolve(a, k[::-1, ::-1], mode="same")
        # F: correlate n * exposure^2 / mu^2 with the PSF squared
        b = np.where(m, counts[e] * (exposure[e] / mu) ** 2, 0.0)
        F += fftconvolve(b, (k ** 2)[::-1, ::-1], mode="same")
    ts = np.where(S > 0, S ** 2 / np.maximum(F, 1e-300), 0.0)
    return ts


def _blocked(ll, bb, found, min_sep, edge):
    """Pixels too close to an existing source or to the ROI edge."""
    block = (np.abs(ll) > HALF - edge) | (np.abs(bb) > HALF - edge)
    for (_, _, l0, b0, *_rest) in found:
        block |= np.hypot((ll - l0) * np.cos(np.radians(bb)), bb - b0) < min_sep
    return block


def run(counts, exposure, background, energies, mask=None, config=None,
        verbose=True):
    """Iterate. `background` is the component set BEFORE any new sources:
    one dict, or one dict per energy bin. Returns (components, fit, found)."""
    cfg = config or load_config()
    ll, bb = _grid()
    comps = ([dict(background) for _ in energies]
             if isinstance(background, dict) else [dict(b) for b in background])

    fit = FL.fit(counts, exposure, comps, energies, mask)
    found = []
    for thresh in cfg["passes"]:
        for _ in range(int(cfg["max_sources_per_pass"])):
            ts = residual_ts_map(counts, exposure, comps, energies, fit, mask)
            ts = np.where(_blocked(ll, bb, found, cfg["min_separation_deg"],
                                   cfg["edge_margin_deg"]), 0.0, ts)
            i, j = np.unravel_index(int(np.argmax(ts)), ts.shape)
            if ts[i, j] < thresh:
                break
            name = f"src{len(found):03d}"
            src = point_source_map(i, j)
            for e in range(len(energies)):
                comps[e][name] = src
            # Warm start: only one component changed, so the previous solution
            # is already almost right. Refitting cold re-derives a hundred
            # normalisations that converged a moment ago, and once the model
            # holds ~35 sources that is minutes per iteration.
            x0 = []
            for e in range(len(energies)):
                prev = dict(fit["bins"][e]["norms"])
                conv = FL.convolve_templates({name: src}, energies[e])[name]
                scale = float((exposure[e] * conv).sum())
                tot = float(counts[e].sum())
                prev[name] = 0.1 * tot / max(len(comps[e]), 1) / max(scale, 1e-30)
                x0.append(prev)
            fit = FL.fit(counts, exposure, comps, energies, mask, x0=x0)
            found.append((name, float(ts[i, j]), float(ll[i, j]), float(bb[i, j]),
                          int(thresh), len(found) + 1))
            if verbose:
                print(f"    + {name} at l={ll[i, j]:+.2f} b={bb[i, j]:+.2f} "
                      f"residual TS {ts[i, j]:.0f} (pass {thresh})", flush=True)
    return comps, fit, found


# ------------------------------------------------------------- self-test
def selftest():
    """Inject point sources into a simulated sky and check the iteration finds
    them, in the right places, and then stops."""
    rng = np.random.default_rng(42)
    ok = True
    ll, bb = _grid()
    th = np.hypot(ll, bb)
    energies = np.array([1.5, 4.0])
    exposure = np.broadcast_to(5.0e11, (2, NPIX, NPIX)).copy()

    sa = np.radians(PIX) ** 2 * np.cos(np.radians(bb))
    flat = np.ones((NPIX, NPIX)); flat /= float((flat * sa).sum())
    plane = np.exp(-(bb / 4.0) ** 2); plane /= float((plane * sa).sum())
    truth_bg = {"iso": 1.0e-11, "diffuse": 6.0e-11}

    # three sources, well separated, away from the edge
    inject = [(6.0, 5.0, 4.0e-9), (-8.0, -3.0, 3.0e-9), (2.0, -9.0, 2.5e-9)]
    counts = np.zeros((2, NPIX, NPIX))
    for e, en in enumerate(energies):
        conv = FL.convolve_templates({"iso": flat, "diffuse": plane}, en)
        mu = exposure[e] * sum(truth_bg[k] * conv[k] for k in conv)
        for (l0, b0, f0) in inject:
            j = int(round((HALF - PIX / 2.0 - l0) / PIX))
            i = int(round((b0 + HALF - PIX / 2.0) / PIX))
            mu = mu + exposure[e] * f0 * FL.convolve_templates(
                {"s": point_source_map(i, j)}, en)["s"]
        counts[e] = rng.poisson(mu)

    cfg = load_config()
    print(f"  stopping rule from configs/level2.yaml: passes {cfg['passes']}, "
          f"min separation {cfg['min_separation_deg']} deg")
    comps, fit, found = run(counts, exposure,
                            {"iso": flat, "diffuse": plane}, energies,
                            config=cfg)

    print(f"  injected {len(inject)} sources, recovered {len(found)}")
    for (l0, b0, f0) in inject:
        d = min(np.hypot((f[2] - l0) * np.cos(np.radians(b0)), f[3] - b0)
                for f in found) if found else 99
        print(f"    injected (l={l0:+.1f}, b={b0:+.1f}): nearest recovered "
              f"{d:.2f} deg away")
        if d > 0.25:
            print("    FAIL not recovered within a quarter degree"); ok = False
    # and it must stop: no spurious extras
    if len(found) > len(inject):
        print(f"  FAIL added {len(found) - len(inject)} spurious source(s)")
        ok = False
    # the background normalisations must survive the procedure
    for k, want in truth_bg.items():
        got = np.mean([b["norms"][k] for b in fit["bins"]])
        print(f"  {k:8} injected {want:.3e} recovered {got:.3e} "
              f"ratio {got/want:.3f}")
        if not 0.9 < got / want < 1.1:
            print(f"  FAIL {k} normalisation moved by more than 10%"); ok = False
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)
