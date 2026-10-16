# Section 3.1 — The squared column, the annihilation integral, and what
# separates them

Paste-ready. Three claims in order of strength: an exact identity, a measured
difference, and a geometric correction. The first needs no data; the second is
measured on the HESTIA particle data in Muru's own frame and estimator; the
third follows from reading the projection code. Supporting numbers:
`docs/muru_conventions.md`, `outputs/gap1/muru_match_tables.md`. Code:
`analysis/gap1/muru_match.py`, self-tested without data.

---

## 3.1.1 Squaring a projected map is a relabeling of its isophotes, not a
## transformation of its shape

Let `m(l,b) >= 0` be any map and let the shape of its level set at a fraction
`f` of the peak be summarised by any functional `Q` of the point set

    S(m, f) = { (l,b) : m(l,b) >= f * max m } .

Because `t -> t^2` is strictly increasing on the non-negative reals,

    m(l,b)^2 >= f * max(m^2)   <=>   m(l,b) >= sqrt(f) * max(m) ,

so

    S(m^2, f) = S(m, sqrt(f))      exactly, for every f in (0, 1].          (1)

Any shape statistic computed on a super-level set therefore satisfies

    Q(m^2, f) = Q(m, sqrt(f)) .                                             (2)

The squared map's one-parameter family of isophote shapes is the unsquared
map's family with the level axis reparameterised. It contains no shape the
unsquared map does not, and omits none that it has.

This is not an approximation and does not depend on the map, the smoothing,
the grid, or the choice of `Q`. It applies in particular to the axis ratio
used in Muru et al. — the ratio of the smaller to the larger singular value
of the mean-centred coordinates of the pixels above `f * max`
(`notebooks/azimuthal_analysis_erebos.jl`) — and to the ellipse-fitted axis
ratio and boxiness used in this work. We verify (2) numerically to the last
representable bit at `f` = 0.7, 0.5 and 0.3 in `muru_match.selftest()`.

**Consequence.** The `dm2` maps are the `dm` maps as far as morphology is
concerned. A measurement that squares a projected column and then compares
shapes is reporting the shape of the column at a shifted isophote. In
particular, the agreement between a squared-column template and the Galactic
Centre Excess cannot be read as evidence about annihilation morphology
specifically, because the same agreement is available from the column itself.

The Comment of Abazajian, Kumar & Macias (arXiv:2609.34155) observes that
`(int rho ds)^2` is not `int rho^2 ds`. Equation (1) sharpens the objection:
the squaring step is not merely the wrong operation, it is a shape-neutral
one, so the published morphology is the morphology of the projected column,
whatever was intended.

## 3.1.2 The annihilation integral is flatter than the squared column by
## Dq = -0.032 +- 0.020

Equation (1) says what squaring *cannot* do; it says nothing about how far
`int rho^2 ds` is from `(int rho ds)^2`. That is an empirical question about
the halos, and the HESTIA particle data answer it.

Method. For each of the six HESTIA galaxies we estimate the local dark-matter
density at every particle from its 32nd nearest neighbour, with the radius
floored at the simulation's 220 pc softening so the estimator cannot invent
structure the simulation does not resolve. Along each sight line we then form

    column:       C(l,b) = sum_i m_i
    annihilation: J(l,b) = sum_i m_i rho_i

over the same particle selection and the same sight-line kernel, so the two
maps differ only by the factor `rho_i`. Observer at 8.0 kpc, particles within
15 kpc of the observer, kernel a top-hat of 3 deg angular radius, grid
`-20..20 deg` at 0.1 deg — Muru's conventions throughout (§1 of
`docs/muru_conventions.md`). Axis ratios use his estimator.

Result, at the 50% isophote over all six galaxies and four viewing angles:

    q(J) - q(C^2) = -0.032 +- 0.020,     same sign in 22 of 24

and, computed in the stellar disc frame instead of the AHF halo frame — an
orientation 17 deg away for two galaxies and 47-48 deg away for the other
four —

    q(J) - q(C^2) = -0.041 +- 0.016,     same sign in 24 of 24 .

The measurement is therefore frame-independent, which matters because the
orientation convention of the published maps cannot be reconstructed from the
public material (§3.1.3).

Interpretation. The squared column over-rounds the halo by about 0.03 to 0.04
in projected axis ratio. This is far below the systematic that dominates every
current GCE morphology measurement — our own interstellar-emission test moves
the detection significance of the same data by a factor of eighty — so no
conclusion in Muru et al. turns on it. It is not below the precision at which
axis ratios are quoted. Templates built from `dm2` should therefore not be
used to support a flattening claim at better than roughly 0.05 in `q`, and
analyses that need the annihilation morphology should build it from the
particle data, which is now straightforward.

## 3.1.3 Two geometric caveats that belong with the above

**The published maps are not column densities.** The sight-line kernel in
`src/BulgeDensityProjection.jl` is a top-hat of fixed *angular* radius, built
as a KD-tree query in `(l,b)` degrees. Its solid angle is the same at every
distance, so its physical cross-section grows as `d^2`, and the sum of
particle masses inside it is

    sum_i m_i  =  int rho(s) s^2 ds  per unit solid angle,

not `int rho ds`. Over a 15 kpc sight line with the observer 8 kpc from the
centre, material behind the Galactic Centre is weighted roughly an order of
magnitude above material in front of it. Shapes in the published figures are
therefore shapes of a `d^2`-weighted column. The normalisation is divided out
and every fit here floats a normalisation per energy bin, so no published
number changes; but a flux derived from these maps would be wrong, and the
`J` and `C` maps in §3.1.2 are deliberately built with the same weighting so
that the comparison isolates the squaring alone.

**The axis ratio is orientation-blind, so it needs no inversion.** Muru's
statistic is `min/max` of the singular values and is `<= 1` by construction
whichever way the structure lies. In all 24 published maps the major axis is
in longitude, so the published number is the vertical over the in-plane
extent — the same convention as `q` in this work, not its reciprocal. Di
Mauro's `r` is the reciprocal and still requires inversion; Muru's does not.
This resolves a convention ambiguity that would otherwise have inverted every
cross-comparison in §6.

**Halo identity in the published maps cannot be verified.** The particle files
as distributed have their names swapped within each pair (the file named
`G11` holds the halo at G1.2's centre; §3 of the Gap 1 notes), caught only
because the AHF centres give an independent handle. The map files have no such
handle, and the maps themselves cannot supply one: their shape residuals
correlate with each other at mean `r = 0.81`, with the same-galaxy
different-angle correlation no higher than the different-galaxy same-angle
correlation. Any per-halo ordering taken from these maps rests on the file
names.
