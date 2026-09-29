# Level 2: iterative source-finding

Shared module. Gap 1 needs it for H₀ᴮ ("is the reported sphericity of the GCE a
property of the sky or of the pipeline?") and Gap 5's phantom-source question
is the same machinery, so it is written once here rather than twice. Offered to
Gap 5 — not yet agreed.

## What it does

Replicates the Di Mauro (2026) procedure: fit, make a residual TS map, add the
most significant new point source above threshold, refit, repeat. Successive
passes at TS > 100, 64, 36.

The stopping rule lives in `configs/level2.yaml` and is fixed before any run.
`iterate.py` refuses to start if it is asked to use anything else, because the
whole test is worthless if the pipeline can be tuned per template.

## The residual TS map

For a candidate point source at pixel p with PSF-convolved, exposure-weighted
template A(p), the one-parameter score test for adding it is

    TS(p) = S(p)^2 / F(p),    S = sum_i A_i (n_i/mu_i - 1),
                              F = sum_i A_i^2 n_i / mu_i^2

evaluated at the current best fit, for S > 0. Both sums are correlations of a
per-pixel image with the PSF, so the whole map is two convolutions per energy
bin rather than a fit per pixel. Same quantity `gttsmap` computes the slow way.

## What the test needs it for

- **Gap 1 H₀ᴮ:** inject a HESTIA-shaped signal, run the full iteration, and
  measure whether the recovered morphology is rounder than what went in.
  A spherical injection at matched flux is the control — without it, "boxy
  becomes round" is equally consistent with a pipeline that mangles every
  morphology.
- **Gap 5:** whether the added sources are real or phantoms.
- **Both:** bookkeeping of every source the iteration adds — position, TS,
  flux, pass number, iteration — so it can be asked whether new sources cluster
  along the major axis of the injected shape.
