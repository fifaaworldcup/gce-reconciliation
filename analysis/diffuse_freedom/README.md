# Diffuse-freedom test (shared)

One question: how much of a reported TS is the diffuse model's error rather
than a signal?

Refit with the interstellar model split into components that each carry a free
normalisation, so it has freedom in shape and not only in scale, and watch the
TS. Splits come from `configs/diffuse_freedom.yaml`, fixed before any run and
identical for every template.

A split on its own proves nothing, so each one runs on three skies:

| sky | what it is | what it tells you |
|---|---|---|
| observed | the data | the number under test |
| injected | a simulation from the fitted background plus a known signal | whether the split can see a signal at all |
| boosted | the data with that known signal added on top | whether it can see one through real residuals |

A split is evidence only if `injected` survives it. `verdict()` drops the rest
and reports the spread of `observed` across the survivors.

## Using it from another gap

`freedom.run()` takes your own components, energies, mask and templates. It
does not know about Gap 1. Pass the diffuse plane separately from everything
else, because that is the thing being split.

```python
import sys; sys.path.insert(0, "analysis/diffuse_freedom")
import freedom as DF

res = DF.run(counts, exposure,
             comps_base=[{"iso": flat, "ps": ps[e]} for e in range(nE)],
             energies=emid, mask=mask,
             signals={"vvv": vvv, "gnfw1.2": gnfw},
             iem_planes=[iem[e] for e in range(nE)])
print(DF.verdict(res)["note"])
```

`comps_base` excludes the diffuse plane. `iem_planes[e]` is the plane for bin
`e` on your grid; the ROI half-widths in the config must match that grid.

A gap whose IEM is already ring-decomposed should split the rings instead of
adding latitude bands: pass the rings as separate entries of `comps_base` and
run the `rigid` split alone as a baseline, then compare against the same gap's
single-template run. That is the comparison that says whether a ring model
removes the ambiguity or only hides it.

## Gap 1's result

`docs/diffuse_freedom_test.md`. TS 3.8 to 306 across splits a control shows
are not degenerate, Mask B + PS, 1–10 GeV, one all-sky IEM.

## Running Gap 1's configuration

```bash
python analysis/diffuse_freedom/run_test.py \
  --ccube gap1_ccube_full.fits --expcube gap1_expcube_full.fits \
  --iem gap1_iem_roi.fits --catalog gll_psc_catalog.fit \
  --mask B --ps --out outputs/gap1/diffuse_freedom_maskB_ps
```
