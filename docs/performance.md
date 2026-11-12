# Performance of the Feynman-trick recursion

How a run spends its time, which options matter, and how the largest family
is computed. Measurements are on a 4-core x86-64 Linux container (GCC 13,
FLINT 3.4.0, Release build); wall times are in [timings.md](timings.md).

## Where the time goes

A run has an exact preparation (IBP reductions per level, cached) and a
numerical phase. In the numerical phase about three quarters is the
transport of endpoint functionals along each level's merge parameter; the
rest is exact endpoint work (series, projections, matching). The JSON
`timings` list every phase, and the phases add up to `numerical_seconds`:
`ordinary_seconds` (transport), `endpoint_gauge_seconds` (the gauge of each
endpoint to Fuchsian form), `endpoint_series_seconds`,
`endpoint_frame_seconds`, `endpoint_observable_seconds` (projection and
regulated integration of every requested row), `endpoint_operator_seconds`,
`basis_selection_seconds` and `leaf_seconds`. For one sample of the double
pentagon (one thread): transport 31 s, observables 10 s, series 6.5 s, basis
scan 3 s, frames and operators 1.2 s each. In a sample of a small family the
scalar leaf is the largest item (four-loop banana: 2.0 of 4.5 s), because it
is still summed as a series in epsilon.

A multiply-add of the Taylor recurrence is one exact coefficient times one
working-precision value, and costs one pass over the value per limb of the
coefficient. Expanded runs have coefficients of one or two limbs (2.4 to 5
passes per product). In a sample the leg is complex and epsilon is a
number, so the shifted coefficients are integers of up to seven limbs: 8 to
18 passes per product. The working precision itself costs little (four-loop
banana, expanded: 27, 30 and 33 s at 256, 384 and 512 bits).

One transport leg costs about

    rows x connection size x epsilon orders x charts x Taylor order.

Evaluated by rows (`--accumulation rows`), the rows of a level are the
integrals its parent needs, so the one-row top level is cheap and the
middle levels dominate. The default, `--accumulation weighted`, evaluates
the same product of local operators from the top: a level receives a weight
for each of its requested integrals from the level above and transports the
one weighted combination of its rows, with the weighted combination of
their forcings, so every level has one row per requested integral of the
family. Such a row costs the whole connection plus all forcing rows of its
level. Counted in multiplications per Taylor order, the inner levels become
cheaper by 6.2 (double pentagon level 2, 69 rows on 33 masters), 3.9 (its
level 3), 4.5 (Higgs+jet family F, level 2), 5.4 (four-loop banana, level
2) and 3.9 (kite, level 2); the top level does not change. These are
counts, not measured times. The epsilon orders carried are the same: a
level's weights start at the poles of the levels above it, and its window
ends that much lower.

Each leg reports its size on stderr. For the planar double box by rows (one
thread):

| Level | Rows x masters | Epsilon orders | Charts (lower / upper) | Seconds (lower / upper) |
| ---: | --- | ---: | --- | --- |
| 5 | 3 x 2 | 11 | 13 / 12 | 0.01 / 0.03 |
| 4 | 7 x 3 | 11 | 15 / 34 | 0.10 / 0.44 |
| 3 | 11 x 7 | 11-12 | 20 / 29 | 1.6 / 2.2 |
| 2 | 12 x 11 | 9-10 | 35 / 26 | 1.7 / 3.3 |
| 1 | 1 x 12 | 8-9 | 34 / 24 | 0.4 / 0.2 |

Transport time is linear in the number of epsilon orders carried.

## What the defaults do

- **Taylor charts on compact numbers.** A chart's recurrence runs on
  contiguous fixed-size limb vectors with scaled-double radii; the ball
  recurrence (`expand_balls`) is the reference and fallback. The Taylor tail
  admitted per chart is ten digits below the requested leaf digits (2^-127
  by default; `--transport-bits`).
- **One weighted row per level** (`--accumulation weighted`). The levels
  are combined from the top (see above). The weights are Laurent series in
  epsilon with ball coefficients. A row's forcing is carried by one constant
  component of the transported state per forcing row, which holds that
  row's weight; the weights' uncertainty travels with the row, and a
  transport checkpoint is keyed by them. Windows are sized from the same
  structural losses as by rows, and the leaf is asked for the same first
  order. An endpoint constraint of a level is carried down as one more row
  and checked on the leaf values, through the requested epsilon order (by
  rows it is checked through the order its level is asked for). A level
  with an endpoint boundary multiplies its weights by the boundary map.
  Evaluated by rows instead: `--method factored` and `auto`, and more than
  four requested integrals (`--henn-basis`), where weighted rows would not
  be fewer. `--accumulation rows` selects that everywhere; the library
  default (`NumericalOptions::accumulation`) is `rows`. In weighted mode
  only the top level has a stored linear expression.
- **Paths.** A leg stays on the real axis and passes a real singular point
  over a triangle on the causal side, which is derived for each point. An
  epsilon sample uses a path lifted to the contour height instead, crossing
  the axis between points that ask for opposite sides.
- **Uncertainty by probes** (`--uncertainty probe`). Midpoints are
  transported and the size of discarded radii is followed along two probe
  rows (one, with the second direction in its imaginary part, beside a
  single transported row). This is an estimate, not an enclosure;
  `--uncertainty map` gives the rigorous chart map at several times the
  cost. Steps are not shortened by inherited or inflated interval radii:
  once a leg transports midpoints, the radii of a chart are the rounding
  that chart itself commits (the integrated defect of its recurrence, not a
  bound propagated through the absolute values of the connection), and a
  step is halved only if that loses half the working precision. Near an
  endpoint where the master basis has higher-order poles the interval bound
  used to force each chart down to 1/128 of its geometric step.
- **Endpoint series.** The exact series at each endpoint has order 22 and is
  summed at most 1/64 of the way to the nearest other singular point, a
  2^-128 tail (`--endpoint-order`, `--overlap-ratio`; samples use 32 and
  1/16). The exact work grows like the order to the power 2.4, while a
  smaller overlap costs the transport about 1.25 charts per halving. The
  overlap is chosen per endpoint, from the distance of that endpoint to its
  own nearest singular point: a singular point next to x = 0 does not move
  the overlap point at x = 1 (the progress line `endpoint overlap=` shows
  both).
- **Moving poles removed** (`--moving-poles remove`). Denominator factors
  that mix the merge parameter with epsilon are removed by exchanging a
  master. This narrows no epsilon window but shrinks the connection where
  such a factor exists; see [moving-poles.md](moving-poles.md). A factor
  whose removal would leave the connection more than four times larger is
  kept (a factor of degree seven in x of a three-loop level: 706235 to
  26971415 characters, after which one column of the level's endpoint series
  took 110 s). The outcome of the search is recorded beside the exact cache
  (`moving-poles/`), so that later runs and the samples of a fit apply the
  recorded steps instead of searching again.
- **Working precision** (`--precision-bits`, 384 by default). Every level
  loses precision between its overlap points and its matching point: 30 to
  140 bits a level in the samples of the large families (the transported
  row is a small combination of solutions that grow towards the endpoint).
  A seven-level sample can lose all 384. When rounding leaves a coefficient
  fewer digits than `--leaf-digits` asks for, the numbers are computed
  again with that many bits more; the exact endpoint plans are kept. A
  precision stated on the command line or in the family file is not
  changed, and nothing is raised beyond 4096 bits.
- **Merge selection** (`--merge-selection auto`). Which two propagators a
  level merges decides how large the levels below become. Candidates are
  the pairs on one loop momentum, whose merged line is again an ordinary
  propagator. For each, the number of master integrals of the merged graph
  is counted from its Symanzik polynomials: the critical points of the
  Lee-Pomeransky polynomial U + F, found with Singular over a finite field
  in seconds (`critical_points.hpp`; it agrees with the IBP counts: 16 for
  the kite, 12 for the double box, 75 and 49 for two merges of the
  Higgs+jet family F, 69 for the double pentagon). The pair with the fewest
  is merged, the first pair on a tie. Without Singular the masters are
  counted by IBP eliminations instead, for levels of 24 or more masters. An explicit `merges` or `anchors` entry, or stated causal signs,
  keep the ladder as written. For family F (seven propagators, four
  massive) the choice turns levels of 75, 47, 21 masters into 49, 19, 13:
  preparation about a minute instead of eleven, one epsilon sample under two
  minutes instead of 27. The kite's second level has 7 masters instead of
  10.
- **Endpoint boundary** (`--boundary auto`). A level needs its masters at
  one value of its merge parameter x, which the level below computes. At an
  interior anchor (1/3, moved if singular) that level is the merged family
  with the kinematics of the anchor. At x = 0 or x = 1 the merged line is
  one of the two merged propagators, so the level below is the family with
  the other one removed, at its original kinematics: fewer masters and
  fewer singular points. A level matches at an endpoint if all of the
  following hold; otherwise it keeps its anchor.
  - The family's U and F have positive monomials, the propagators that make
    up the two lines share one loop momentum, and the line that remains is
    massive. The level's integrals are then analytic at the endpoint: U
    depends on the two parameters through their sum, and F(x) - F(0) is
    bounded by a multiple of |x| F(0) on the whole simplex. With a massless
    remaining line the endpoint has collinear regions (exponents -eps,
    -2 eps) that the smaller family does not provide, and no endpoint is
    taken.
  - No master of the level has a numerator that depends on the propagators
    at the endpoint (an unused numerator slot that does is replaced).
  - The level's equation, expanded exactly at the endpoint, confirms that
    the masters' values there fix the analytic solution: the combinations of
    local solutions with integer exponents that have no negative power and
    no logarithm are determined by their constant terms. Solutions with any
    other exponent (x^(-2 eps), x^(-1/2)) have no part in an analytic
    function. A level whose analytic solutions would need Taylor
    coefficients beyond the values keeps its anchor.

  The masters at the overlap point next to the endpoint are the exact series
  of that solution applied to the masters at the endpoint
  (`endpoint_boundary.hpp`); the functionals of the other endpoint are
  transported to that point in one arm. Such a level uses the adjoint
  method. Explicit `anchors`, stated causal sides and `--henn-basis` keep
  the anchors; `--boundary anchor` keeps them everywhere, and
  `--boundary endpoint` fails if no level qualifies. The JSON report lists
  each level's masters and boundary under `ladder`.
  The kite matches its first two levels at x = 0 (13 of 16 and 4 of 7 local
  solutions are analytic) and its third level has 2 masters instead of 3;
  `pentagon_massive` matches all three levels at x = 0. For the Higgs+jet
  family F the first and third levels match at x = 0 (22 of the 49 local
  solutions of the first are analytic) and the levels are 49, 12, 8, 7, 3
  instead of 49, 19, 13, 8, 3: the transport below the first level takes
  1.5 s instead of 12 s per epsilon sample, while the first level, which
  dominates the sample, is unchanged. Results agree with the interior
  anchors to the accuracy of the run (41 digits for the kite, 34 for a
  sample of family F). Families without a massive pair on one loop momentum
  (the bananas, the massless boxes and pentagon) are computed exactly as
  with `--boundary anchor`.
- **Solver limits.** A reduction of the built-in IBP solver that meets its
  degree or sample bound is repeated with the bounds doubled; the time and
  memory budgets (a day and 4 GiB by default) are what ends a preparation.
- **Threads.** With OpenMP the streams and rows of a Taylor order run in
  parallel once an order costs 5000 multiply-adds; endpoint series solve
  their columns on separate threads; the built-in IBP solver works ahead on
  free threads ([ibp-solver.md](ibp-solver.md)). Results do not depend on
  the thread count. The CLI leaves cores that are busy at start-up alone;
  `OMP_NUM_THREADS` and `--ibp-threads` override.

Options that do not help: a lower working precision (the four-loop banana
loses about 160 bits expanded and 180 in a sample), a different contour
height, and `--method factored`, which transports the solution map instead
of the rows and is slower and less accurate on every family measured.

## Epsilon sampling

`diffexp ft --epsilon-strategy sample` prepares the recursion once,
evaluates the integral at the shifted dimensions epsilon = +-k*step in
parallel processes and fits the Laurent coefficients
(`epsilon_sampling.hpp`). At a generic shift no level has epsilon poles,
every level needs only order zero of its child, and the levels above the
leaf are specialised to the numeric dimension.

- `--samples N` (default 2*loops + epsilon order + 5, at least 8),
  `--sample-step Q` (default 1/400) and `--sample-jobs N` (default one per
  idle hardware thread) shape the run. The family file takes
  `epsilon_strategy`, `samples` and `sample_step` under `numerical`.
- The order of the pole is detected from the growth of the samples towards
  zero (at most 2*loops); `--sample-poles N` fixes it.
- The fit is the polynomial through epsilon^poles times the samples, in
  ball arithmetic. Each coefficient's radius is widened by the change it
  makes to drop the two outermost samples: an estimate of the interpolation
  error, not a bound. JSON also reports it as
  `fit_change_without_outermost_samples`.
- The scalar leaf is re-expanded from its Laurent series at zero, which
  costs digits in deep recursions: samples ask the leaf for at most 16
  digits, 12 from six levels on, unless `--leaf-digits` is given.
- One sample alone is `--epsilon-shift Q --epsilon-order 0`.
- A sample that had to raise its working precision (see above) reports it,
  and the samples started after it begin at that precision
  (`--precision-floor`). For a three-loop family with nine propagators and
  seven levels a sample has 8 digits at 384 bits and 18 at 448, where the
  Taylor tail of 16 leaf digits takes over.

Accuracy of the finite part of the four-loop equal-mass banana from samples
at +-k/400 (reference 39.655526834297652...):

| Samples | Error |
| ---: | ---: |
| 2 | 2.0e-2 |
| 4 | 1.5e-5 |
| 6 | 1.8e-8 |
| 8 | 3.3e-11 |

Each pair of samples buys about 2.7 digits there. Sampling is slower than
the expanded recursion for the small families and is the only practical
route for the double pentagon, whose inner levels carry epsilon through 12
to 16 orders when expanded. Samples are independent, so with one core per
sample the wall time is that of one sample.

## The double pentagon

`henn_double_pentagon_x0` (two loops, eight propagators, seven levels, 69
masters at the top) prepares in 158 s on four threads. One sample at
`--leaf-digits 12` takes about 150 s on four threads. Twelve samples at
epsilon = +-k/400 give epsilon^4 times the scalar integral:

| Order | Value | Change from the last two samples |
| --- | --- | ---: |
| epsilon^-4 | -0.77777777777777779 + 0 i | 1e-14 |
| epsilon^-3 | -0.71244912607989874 - 17.453292519943297 i | 3e-14 |
| epsilon^-2 | 436.218113787 - 5.202032870 i | 2e-9 |
| epsilon^-1 | -88.950617318 + 4207.353609604 i | 6e-9 |
| epsilon^0 | -29747.6631 - 3507.3533 i | 1e-4 |

There is no independent reference for this integral, so it is not an
accepted result. Its leading terms have the structure they must have: the
epsilon^-4 coefficient is -7/9 to 14 digits, and the epsilon^-3 coefficient
matches -1 - (5/9) log 3 + (14/9) gamma_E - (50/9) pi i to 15 digits. The
sample values are in
[validation/henn-double-pentagon-epsilon-samples.json](validation/henn-double-pentagon-epsilon-samples.json).

## Limits

- A deep recursion with many scales loses far more precision expanded than
  sampled. The three-loop nine-propagator family above (levels of 24, 28,
  29, 19, 12, 6, 2 masters, epsilon windows 17 to 23 orders wide) completes
  its expanded run in about ten minutes and returns radii of 1e+1301 at 384
  bits; no precision below 4096 bits repairs that, and none is tried. Its
  samples take 26 to 37 s each.

- The expanded recursion of the double pentagon does not finish in hours:
  levels 2 and 1 have 64 x 33 and 69-master systems with epsilon through 12.
- The exact endpoint series cannot be continued numerically without
  redoing the exact pole-cancelling step on the full series; a timing-only
  run suggests about 1.5 times for the four-loop banana if it could.
- In a sample the endpoint observables and series are still exact
  arithmetic, about 20 s of 150 s for the double pentagon.
