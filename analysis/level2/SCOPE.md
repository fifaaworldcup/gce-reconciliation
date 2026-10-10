# Level 2: what exists, what is missing, and who decides what

Written so whoever owns this can start on Monday. Gap 1 needs Level 2 for
H₀ᴮ; Gap 5 needs the same machinery for phantom sources; and after the Pohl22
re-run it is also **the only remaining route to deciding Macias versus
Di Mauro**, because Mask A — the other candidate explanation — turned out to
have no sensitivity at all (Gap 1 §6.0b, Table R11).

Offered to Gap 5. Not yet agreed. Nothing below assumes an owner.

## Why this is now on the critical path

Gap 1's Pohl22 result says: on a ring-decomposed background, the stellar
bulge beats every dark-matter template by ΔAIC +888, the fit demonstrably
picks the correct morphology class on six simulated skies, and the ordering
survives in the one energy band free of IC interpolation. Macias's answer, on
this background.

Di Mauro reports the opposite with ring-decomposed backgrounds too. The
difference between the two analyses that we have **not** controlled is that he
runs iterative source finding — he lets the fit add point sources to the
model, Level 2 — and we mask a fixed catalogue. Two candidate explanations
existed. Mask A eliminated itself: with the real 4FGL it keeps 8.7% of the ROI
over 1–10 GeV and 0.4% in the lowest bin, and nothing injected into it is
recoverable. That leaves the pipeline. Level 2 is how we test it.

**The experiment, stated before it is run:** take the Pohl22 `rings_species`
background and the same 40 templates, run Level 2 iterative source finding
instead of catalogue masking, and re-run the ordering. If the bulge still
wins, Di Mauro's dark-matter preference is not a source-finding effect and the
disagreement is somewhere we have not looked. If dark matter wins, iterative
source finding is the hinge, and the next question is whether the sources it
adds are real — which is Gap 5's question, reached from the other direction.

## What already exists and is tested

- `configs/level2.yaml` — the stopping rule, frozen: passes at residual TS
  100, 64, 36; at most 40 sources per pass; no new source within 0.30° of an
  existing one; candidates within 1.0° of the ROI edge ignored; refit after
  every addition; seed 42. `iterate.py` refuses to run with anything else,
  which is the whole point — a pipeline tunable per template would
  manufacture exactly the bias Gap 1 exists to measure.
- `iterate.py::residual_ts_map` — the score test for adding one point source
  at every pixel, as two PSF cross-correlations per energy bin rather than a
  fit per pixel. Same quantity `gttsmap` computes the slow way.
- `iterate.py::run` — the loop, with a warm start so the hundred-odd
  normalisations are not re-derived cold after every addition.
- `iterate.py::selftest` — injects three point sources into a simulated sky
  and checks the iteration finds them, in the right places, and stops.

## What is missing

**1. It has never been run on real data.** The self-test is synthetic. Before
any science it needs one run on the real counts cube with the Pohl22
background, with wall-clock and source count recorded. At ~35 sources and a
refit per addition this is the expensive piece of the project.

**2. No recovery control.** Finding sources is not the test. The test is
whether the *morphology* survives the iteration. That needs, at minimum:
  - inject a HESTIA-shaped signal at a known flux, run the iteration, and
    measure the axis ratio of what is recovered;
  - inject a **spherical** signal at matched flux as the control. Without it,
    "boxy comes out round" is equally consistent with a pipeline that rounds
    off everything, and we would be making the §6.1 mistake a fourth time;
  - inject **nothing** as the null, to see how many sources the iteration adds
    to a sky that contains none.

**3. No phantom bookkeeping yet usable by Gap 5.** `run()` returns
`found` as a list of tuples `(name, ts, l, b, pass, index)`. Gap 5 needs that
written out with, per source, the fitted flux per energy bin, whether a 4FGL
source lies within the PSF, and the iteration at which it appeared. Small
change, but it is the interface — worth agreeing before either of us builds
on it.

**4. Two design questions that need a decision, not a guess.**
  - *Does the signal template stay in the model during source finding?*
    Di Mauro fits the GCE and the sources together. If we instead find sources
    against a background with no signal component, the iteration will absorb
    the excess into point sources and the answer is predetermined. This is the
    single most consequential choice in the module and it should be made
    explicitly, written into `level2.yaml`, and run both ways once.
  - *Energy-dependent or energy-integrated source finding?* A residual TS map
    summed over bins finds different sources than per-bin maps. Di Mauro does
    the former; ours currently does the former. Worth confirming against his
    text rather than our memory.

**5. Admissibility, per template, as everywhere else.** Whatever Level 2
concludes, it needs the same control the ring work needed: inject each
template at a stated flux, check it is recovered above the pre-registered
threshold. Gap 1 §6.0b found that a background admissible for the stellar
bulge was *not* admissible for gNFW at the same flux. Assume nothing carries
over.

## Interface

```python
from analysis.level2 import iterate

comps, fit, found = iterate.run(
    counts,        # (nE, 400, 400) counts cube
    exposure,      # (nE, 400, 400) cm^2 s, same grid
    background,    # one dict of named templates, or one per energy bin
    energies,      # (nE,) bin centre energies in GeV
    mask=None,     # (nE, 400, 400) 0/1, or None
    config=None,   # defaults to configs/level2.yaml; do not override per run
)
```

`comps` is the background plus every added source, ready to pass to
`fit_likelihood.fit`. `fit` is the final fit record. `found` lists the added
sources. Templates go in **raw**: `fit_likelihood` does the PSF convolution.

## Suggested division, if Gap 5 takes it

The module is one piece of machinery answering two questions, so splitting it
by file would be worse than splitting it by question.

- **Shared, build once, agree first:** items 3 and 4 above — the source
  record and the two design decisions. An afternoon, and everything else
  depends on it.
- **Gap 1 side:** the morphology-recovery controls (item 2) and the
  ordering re-run on the Pohl22 background.
- **Gap 5 side:** whether the added sources are real — cross-match against
  4FGL, the distribution of their fluxes against the known source count, and
  whether they cluster along the major axis of an injected shape.

The last of those is the interesting one, and it is where the two questions
meet: a pipeline that invents sources along the major axis of whatever it is
given is both Gap 5's phantom problem and Gap 1's rounding mechanism.
