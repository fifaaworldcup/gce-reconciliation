"""
Binned Poisson likelihood template fit, in pure Python.

Why this exists: Fermitools is conda-only and can't be installed in our cloud
sandbox, so the Fermipy path has to run on a laptop. But the actual maths of a
Level-1 template comparison is small, and we only need Fermitools to MAKE the
counts cube and exposure cube. Given those two files this module does the rest.

Model, per energy bin E:

    mu_i = exposure_i * sum_k  N_kE * (PSF_E convolved with T_k)_i

Each component gets its own free normalisation in each energy bin. No spectral
model is assumed - that is the bin-by-bin approach used for SED extraction, and
it means a wrong spectral assumption can't leak into a morphology comparison,
which is the whole point of gap 1.

Poisson log-likelihood, masked pixels dropped:

    lnL = sum_i  n_i * ln(mu_i) - mu_i

Templates go in RAW. This module does the PSF convolution (utils/psf.py),
because unlike gtsrcmaps nothing else here will.

Not a replacement for Fermipy, and not the Di Mauro Level-2 iterative pipeline.
It is an independent Level-1 implementation - useful now, and useful later as a
cross-check on the Fermipy numbers.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "utils"))
import psf as PSF

NPIX = 400
PIX = 0.1


# ----------------------------------------------------------------- model
# Convolutions dominate the cost and most of them are recomputed unchanged.
# Level 2 refits after every source it adds, and the sources already in the
# model never change, so without a cache the iteration costs O(N^2)
# convolutions of identical maps - 30 minutes instead of 3 for one run.
# Keyed on the IDENTITY of the template array, and the array itself is held in
# the cache entry so the id cannot be recycled under us. Templates are never
# mutated in place anywhere in this codebase; if that changes, this breaks.
_CONV_CACHE = {}
_CONV_CACHE_MAX = 400          # x 400x400 float64 = about 500 MB ceiling


def clear_conv_cache():
    _CONV_CACHE.clear()


def convolve_templates(templates, energy_gev):
    """PSF-convolve each template at one energy. Returns dict name -> 2D."""
    out = {}
    for k, v in templates.items():
        key = (id(v), float(energy_gev))
        hit = _CONV_CACHE.get(key)
        if hit is not None and hit[0] is v:
            out[k] = hit[1]
            continue
        c = PSF.convolve(v, energy_gev)
        if len(_CONV_CACHE) >= _CONV_CACHE_MAX:
            _CONV_CACHE.pop(next(iter(_CONV_CACHE)))
        _CONV_CACHE[key] = (v, c)
        out[k] = c
    return out


def predict(norms, tmpl_conv, exposure):
    """mu = exposure * sum_k N_k T_k, with N_k >= 0 enforced by the caller."""
    mu = np.zeros_like(exposure)
    for n, t in zip(norms, tmpl_conv):
        mu += n * t
    return exposure * mu


def poisson_lnL(norms, tmpl_conv, exposure, counts, mask):
    mu = np.maximum(predict(norms, tmpl_conv, exposure)[mask], 1e-12)
    n = counts[mask]
    return float(np.sum(n * np.log(mu) - mu))


def fit_one_bin(counts, exposure, templates, energy_gev, mask=None, x0=None,
                max_iter=150, tol=1e-10):
    """Maximise the Poisson likelihood in one energy bin, by ML-EM.

    The first version of this used L-BFGS-B on log(N) and it was wrong: with
    the sharp likelihood surface you get from a million photons it stopped
    short of the maximum, and the resulting TS of an ABSENT component came out
    NEGATIVE - impossible, since the bigger model can always set that
    normalisation to zero and match the smaller one.

    The Poisson likelihood is LINEAR in the normalisations, so the right tool
    is the multiplicative EM update

        N_k  <-  N_k * (sum_i n_i T_ki / mu_i) / (sum_i E_i T_ki)

    which is monotonic in the likelihood by construction and can never go
    negative. No step size, no line search, nothing to tune.

    Returns (norms dict, lnL, info).
    """
    if mask is None:
        mask = np.ones(counts.shape, bool)
    mask = mask & np.isfinite(exposure) & (exposure > 0)

    names = list(templates)
    conv_all = convolve_templates(templates, energy_gev)
    T = np.stack([conv_all[k][mask] for k in names])       # (k, npix)
    E = exposure[mask]
    n = counts[mask]

    A = E[None, :] * T                                     # A_ki = E_i T_ki
    denom = A.sum(1)                                       # sum_i A_ki
    denom = np.where(denom > 0, denom, np.inf)

    if x0 is not None:
        N = np.array([max(float(x0[k]), 1e-30) for k in names])
    else:
        N = (max(n.sum(), 1.0) / len(names)) / denom

    lnL_prev = -np.inf
    N_prev = np.zeros_like(N)
    it = 0
    for it in range(1, max_iter + 1):
        mu = np.maximum((N[:, None] * A).sum(0), 1e-12)
        # numerator is sum_i n_i A_ki / mu_i -- the A, not T. getting this
        # wrong (using T) rescales every normalisation by the exposure and
        # the fit silently returns values ~1e10 too small.
        N = N * ((n / mu)[None, :] * A).sum(1) / denom
        N = np.maximum(N, 0.0)
        if it % 25 == 0 or it == max_iter:
            mu_now = np.maximum((N[:, None] * A).sum(0), 1e-12)
            lnL = float(np.sum(n * np.log(mu_now) - mu_now))
            # require BOTH the likelihood and the parameters to be stable. a
            # component seeded near zero grows multiplicatively and barely
            # moves lnL for many iterations, so an lnL-only test declares
            # convergence while it is still climbing - which returned TS = 0
            # for a component that was plainly in the data.
            moved = np.max(np.abs(N - N_prev) / np.maximum(N, 1e-300))
            if (np.isfinite(lnL_prev)
                    and abs(lnL - lnL_prev) < tol * max(abs(lnL), 1.0)
                    and moved < 1e-4):
                break
            lnL_prev, N_prev = lnL, N.copy()

    # EM gets into the right basin fast but then crawls when components are
    # broad and overlapping (iso, bubbles, diffuse, signal all are), and on the
    # real sky it hit max_iter in every bin. An unconverged fit understates
    # lnL by an amount that differs between templates, which is exactly what
    # biases a delta-lnL. So finish with Newton, which converges quadratically
    # and comes with an honest stopping rule.
    N, pol = _newton_polish(N, A, n)
    mu = np.maximum((N[:, None] * A).sum(0), 1e-12)
    lnL = float(np.sum(n * np.log(mu) - mu))
    return dict(zip(names, N)), lnL, {"iterations": it,
                                      "em_converged": it < max_iter,
                                      "newton_iterations": pol["iterations"],
                                      "lnL_gap_bound": pol["gap"],
                                      "converged": pol["converged"]}


def _newton_polish(N, A, n, max_iter=200, gap_tol=1e-6):
    """Projected Newton on the Poisson likelihood, from the EM solution.

    lnL(N) = sum_i n_i ln(mu_i) - mu_i with mu = N @ A is concave in N, so a
    Newton step on the free parameters, with a backtracking line search and
    non-negativity kept by clamping, reaches the maximum in a handful of steps.

    Parameters are rescaled to x_k = N_k * sum_i A_ki, the expected counts of
    each component. The raw N_k are ~1e-12 and badly scaled; the counts are
    O(1e5) and well conditioned. That scaling was the reason L-BFGS-B failed
    in the very first version of this file.

    Stopping: half the Newton decrement, g.d/2, estimates how far lnL still is
    below its maximum. Stop when that is below gap_tol, i.e. when the fit is
    provably within 1e-6 in lnL of the optimum on the free set, and every
    component pinned at zero has a gradient pointing further negative (KKT).
    """
    s = A.sum(1)
    s = np.where(s > 0, s, 1.0)
    B = A / s[:, None]
    x = np.maximum(N * s, 0.0)

    def lnl(x):
        mu = np.maximum(x @ B, 1e-300)
        return float(np.sum(n * np.log(mu) - mu)), mu

    f, mu = lnl(x)
    gap = np.inf
    it = 0
    ok = False
    for it in range(1, max_iter + 1):
        g = B @ (n / mu - 1.0)
        free = (x > 0) | (g > 0)
        if not free.any():
            ok = True; gap = 0.0; break
        W = n / (mu * mu)
        Bf = B[free]
        H = (Bf * W[None, :]) @ Bf.T               # minus the Hessian, PSD
        gf = g[free]
        try:
            d = np.linalg.solve(H + 1e-12 * np.trace(H) * np.eye(len(gf)), gf)
        except np.linalg.LinAlgError:
            d = np.linalg.lstsq(H, gf, rcond=None)[0]
        gap = 0.5 * float(gf @ d)
        kkt = np.all(g[~free] <= 1e-8 * max(1.0, float(np.abs(g).max())))
        if gap < gap_tol and kkt:
            ok = True
            break
        step = np.zeros_like(x); step[free] = d
        # largest step that keeps every free parameter >= 0
        neg = step < 0
        tmax = float(np.min(-x[neg] / step[neg])) if neg.any() else np.inf
        t = min(1.0, tmax)
        improved = False
        for _ in range(60):
            xn = np.maximum(x + t * step, 0.0)
            if t == tmax:
                # whatever hit the boundary sits exactly on it
                xn[neg & (np.abs(x + t * step) < 1e-9 * np.maximum(x, 1.0))] = 0.0
            fn, mun = lnl(xn)
            if fn >= f + 1e-4 * t * float(g @ (xn - x)) - 1e-12 * abs(f):
                improved = True
                break
            t *= 0.5
        if not improved:
            # no ascent possible along the Newton direction: at the optimum to
            # within floating point
            ok = gap < 1e3 * gap_tol
            break
        x, f, mu = xn, fn, mun

    return x / s, {"iterations": it, "gap": float(gap), "converged": bool(ok)}


def bin_templates(templates, e):
    """templates is either one dict reused in every bin, or a list of dicts,
    one per bin. the galactic diffuse model changes shape with energy, so it
    has to be the list form - a single dict would freeze it at one energy.
    """
    if isinstance(templates, (list, tuple)):
        return templates[e]
    return templates


def component_names(templates):
    return list(bin_templates(templates, 0))


def drop_component(templates, target):
    """Same structure back, minus one component."""
    if isinstance(templates, (list, tuple)):
        return [{k: v for k, v in t.items() if k != target} for t in templates]
    return {k: v for k, v in templates.items() if k != target}


def fit(counts_cube, exposure_cube, templates, energies_gev, mask_cube=None,
        x0=None):
    """Fit every energy bin. counts/exposure are (nE, 400, 400)."""
    out = {"bins": [], "lnL_total": 0.0}
    for e in range(counts_cube.shape[0]):
        m = None if mask_cube is None else mask_cube[e].astype(bool)
        norms, lnL, info = fit_one_bin(counts_cube[e], exposure_cube[e],
                                       bin_templates(templates, e),
                                       energies_gev[e], m,
                                       x0=None if x0 is None else x0[e])
        out["bins"].append({"energy_gev": float(energies_gev[e]),
                            "norms": {k: float(v) for k, v in norms.items()},
                            "lnL": lnL, "converged": bool(info["converged"]),
                            "iterations": info["iterations"]})
        out["lnL_total"] += lnL
    return out


def ts_of(counts_cube, exposure_cube, templates, energies_gev, target,
          mask_cube=None):
    """TS of one component: 2 * (lnL with it - lnL without it)."""
    without = drop_component(templates, target)
    if not component_names(without):
        raise ValueError("need at least one other component")
    wo = fit(counts_cube, exposure_cube, without, energies_gev, mask_cube)

    # Run the bigger fit two ways and keep the better one.
    #   (a) seeded from the smaller solution, extra component small but NOT
    #       negligible - seeding it at ~0 leaves EM climbing for thousands of
    #       iterations and it reports a spurious TS = 0.
    #   (b) from the default cold start.
    # taking the higher likelihood of the two keeps TS >= 0 in practice, which
    # it must be: the bigger model can always reproduce the smaller one.
    warm = []
    for e, bn in enumerate(wo["bins"]):
        tot = float(counts_cube[e].sum())
        share = 0.1 * tot / max(len(component_names(templates)), 1)
        tb = bin_templates(templates, e)
        conv = convolve_templates({target: tb[target]},
                                  energies_gev[e])[target]
        scale = float((exposure_cube[e] * conv).sum())
        warm.append({**bn["norms"], target: share / max(scale, 1e-30)})

    a = fit(counts_cube, exposure_cube, templates, energies_gev, mask_cube,
            x0=warm)
    b = fit(counts_cube, exposure_cube, templates, energies_gev, mask_cube)
    with_it = a if a["lnL_total"] >= b["lnL_total"] else b

    ts = 2.0 * (with_it["lnL_total"] - wo["lnL_total"])
    if ts < -1e-6:
        print(f"    warning: TS = {ts:.3g} < 0, the fit did not converge")
    return max(ts, 0.0), with_it, wo


def residual_map(counts_cube, exposure_cube, templates, energies_gev, result,
                 mask_cube=None):
    """(data - model) summed over energy, in sigma units per pixel."""
    resid = np.zeros((NPIX, NPIX))
    var = np.zeros((NPIX, NPIX))
    nused = np.zeros((NPIX, NPIX))
    for e in range(counts_cube.shape[0]):
        tb = bin_templates(templates, e)
        conv = convolve_templates(tb, energies_gev[e])
        norms = result["bins"][e]["norms"]
        mu = exposure_cube[e] * sum(norms[k] * conv[k] for k in tb)
        m = np.ones(mu.shape, bool) if mask_cube is None else mask_cube[e].astype(bool)
        resid += np.where(m, counts_cube[e] - mu, 0.0)
        var += np.where(m, np.maximum(mu, 1e-12), 0.0)
        nused += m
    # a pixel masked in every bin contributed nothing, so it is not a zero
    # residual, it is no measurement. say so with nan instead of drawing a hole
    # in the middle of the significance map and letting someone read it as fit
    # quality.
    sig = resid / np.sqrt(np.maximum(var, 1e-12))
    sig[nused == 0] = np.nan
    return sig


# ------------------------------------------------------------- self-test
def selftest():
    """Closed loop: build a fake sky from known normalisations, fit, check
    we get them back. If injected and recovered disagree, the fitter is wrong
    and nothing downstream can be trusted."""
    rng = np.random.default_rng(42)
    ok = True

    l = (20 - PIX / 2) - PIX * np.arange(NPIX)
    b = (-20 + PIX / 2) + PIX * np.arange(NPIX)
    ll, bb = np.meshgrid(l, b)
    th = np.hypot(ll, bb)

    # three components with quite different shapes
    flat = np.ones((NPIX, NPIX))
    plane = np.exp(-(bb / 4.0) ** 2)
    gce = np.exp(-(th / 5.0) ** 2)
    sa = np.radians(PIX) ** 2 * np.cos(np.radians(bb))
    templates = {k: v / float((v * sa).sum())
                 for k, v in (("iso", flat), ("diffuse", plane), ("gce", gce))}

    energies = np.array([0.7, 1.7, 5.5, 22.0])
    exposure = np.broadcast_to(5.0e10 * (1 - 0.15 * (th / 28.0)),
                               (len(energies), NPIX, NPIX)).copy()

    # pick normalisations that give a REALISTIC number of photons. the first
    # version of this test simulated 1.4e12 counts, which is a million times
    # what Fermi actually collects here, and the absurdly sharp likelihood
    # hid the convergence bug rather than exposing it.
    want = {"iso": 2.0e5, "diffuse": 9.0e5, "gce": 1.5e5}    # counts per bin
    truth = {}
    for k in ("iso", "diffuse", "gce"):
        conv = convolve_templates({k: templates[k]}, energies[0])[k]
        truth[k] = want[k] / float((exposure[0] * conv).sum())

    counts = np.zeros((len(energies), NPIX, NPIX))
    for e, en in enumerate(energies):
        conv = convolve_templates(templates, en)
        mu = exposure[e] * sum(truth[k] * conv[k] for k in templates)
        counts[e] = rng.poisson(mu)

    print(f"  simulated {counts.sum():,.0f} photons over {len(energies)} bins")

    res = fit(counts, exposure, templates, energies)
    print(f"  {'component':10} {'injected':>10} {'recovered':>10} {'ratio':>7}")
    for k in templates:
        rec = np.mean([bn["norms"][k] for bn in res["bins"]])
        ratio = rec / truth[k]
        print(f"  {k:10} {truth[k]:10.3e} {rec:10.3e} {ratio:7.3f}")
        if not 0.95 < ratio < 1.05:
            print(f"    FAIL {k} recovered {ratio:.1%} of injected"); ok = False

    if not all(bn["converged"] for bn in res["bins"]):
        print("  FAIL a bin did not converge"); ok = False

    ts, _, _ = ts_of(counts, exposure, templates, energies, "gce")
    print(f"  TS of the injected GCE: {ts:,.0f}  ({np.sqrt(max(ts,0)):.0f} sigma)")
    if ts < 100:
        print("  FAIL injected GCE should be detected at high TS"); ok = False

    # null test: a GCE that is NOT in the data must come back with low TS
    counts0 = np.zeros_like(counts)
    for e, en in enumerate(energies):
        conv = convolve_templates(templates, en)
        mu = exposure[e] * sum(truth[k] * conv[k] for k in ("iso", "diffuse"))
        counts0[e] = rng.poisson(mu)
    ts0, _, _ = ts_of(counts0, exposure, templates, energies, "gce")
    print(f"  TS of an absent GCE:    {ts0:,.1f}  (should be small)")

    # the per-bin template form is a separate code path, so check it agrees.
    # same templates, just handed over as one dict per bin - the answer must be
    # identical or the energy-dependent diffuse model can't be trusted.
    per_bin = [dict(templates) for _ in energies]
    res_pb = fit(counts, exposure, per_bin, energies)
    dlnL = abs(res_pb["lnL_total"] - res["lnL_total"])
    print(f"  per-bin template form: dlnL vs shared form = {dlnL:.2e}")
    if dlnL > 1e-6 * max(abs(res["lnL_total"]), 1.0):
        print("  FAIL per-bin templates gave a different likelihood"); ok = False
    ts_pb, _, _ = ts_of(counts, exposure, per_bin, energies, "gce")
    if abs(ts_pb - ts) > max(1e-3 * abs(ts), 1.0):
        print(f"  FAIL per-bin TS {ts_pb:.1f} != shared TS {ts:.1f}"); ok = False

    # and check it actually USES the per-bin maps: swap the gce shape out for
    # noise in half the bins and the TS must drop
    bad = [dict(templates) for _ in energies]
    for e in range(len(energies) // 2):
        bad[e] = {**templates, "gce": templates["iso"]}
    ts_bad, _, _ = ts_of(counts, exposure, bad, energies, "gce")
    print(f"  per-bin form is live: TS {ts:,.0f} -> {ts_bad:,.0f} "
          f"when half the bins get the wrong shape")
    if ts_bad >= ts:
        print("  FAIL per-bin templates are being ignored"); ok = False
    if ts0 > 30:
        print("  FAIL absent GCE detected - the fitter invents signal"); ok = False

    r = residual_map(counts, exposure, templates, energies, res)
    print(f"  residual map: mean {r.mean():+.3f} sigma, std {r.std():.3f} "
          f"(want ~0 and ~1)")
    if abs(r.mean()) > 0.05 or not 0.8 < r.std() < 1.25:
        print("  FAIL residuals are not unit-normal"); ok = False

    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return ok


if __name__ == "__main__":
    selftest()
