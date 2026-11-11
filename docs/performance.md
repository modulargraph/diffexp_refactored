# Performance of the Feynman-trick recursion

How a run spends its time, which options matter, and how the largest family
is computed. Measurements are on a 4-core x86-64 Linux container (GCC 13,
FLINT 3.4.0, Release build); wall times are in [timings.md](timings.md).

## Where the time goes

A run has an exact preparation (IBP reductions per level, cached) and a
numerical phase. In the numerical phase about three quarters is the
transport of endpoint functionals along each level's merge parameter; the
rest is exact endpoint work (series, projections, matching).

One transport leg costs about

    rows x connection size x epsilon orders x charts x Taylor order.

The rows are the integrals the parent level needs, so the one-row top level
is cheap and the middle levels dominate. Each leg reports its size on
stderr. For the planar double box (one thread):

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
- **Paths.** A leg stays on the real axis and passes a real singular point
  over a triangle on the causal side, which is derived for each point. An
  epsilon sample uses a path lifted to the contour height instead, crossing
  the axis between points that ask for opposite sides.
- **Uncertainty by probes** (`--uncertainty probe`). Midpoints are
  transported and the size of discarded radii is followed along two probe
  rows. This is an estimate, not an enclosure; `--uncertainty map` gives the
  rigorous chart map at several times the cost.
- **Endpoint series.** The exact series at each endpoint has order 22 and is
  summed at most 1/64 of the way to the nearest other singular point, a
  2^-128 tail (`--endpoint-order`, `--overlap-ratio`; samples use 32 and
  1/16). The exact work grows like the order to the power 2.4, while a
  smaller overlap costs the transport about 1.25 charts per halving.
- **Moving poles removed** (`--moving-poles remove`). Denominator factors
  that mix the merge parameter with epsilon are removed by exchanging a
  master. This narrows no epsilon window but shrinks the connection where
  such a factor exists; see [moving-poles.md](moving-poles.md).
- **Merge selection** (`--merge-selection auto`). Which two propagators a
  level merges decides how large the levels below become. Candidates are
  the pairs on one loop momentum, whose merged line is again an ordinary
  propagator. For each, the number of master integrals of the merged graph
  is counted from its Symanzik polynomials: the critical points of the
  Lee-Pomeransky polynomial U + F, found with Singular over a finite field
  in seconds (`critical_points.hpp`; it agrees with the IBP counts: 16 for
  the kite, 12 for the double box, 75 and 49 for two merges of the
  Higgs+jet family F, 69 for the double pentagon). A pair with fewer masters
  than the first pair is taken once a quick IBP count confirms that its
  level is smaller, since a smaller graph is not always a cheaper reduction.
  Without Singular only IBP counts are used, for levels of 24 or more
  masters. An explicit `merges` or `anchors` entry, or stated causal signs,
  keep the ladder as written. For family F (seven propagators, four
  massive) the choice turns levels of 75, 47, 21 masters into 49, 19, 13:
  preparation about a minute instead of eleven, one epsilon sample under two
  minutes instead of 27. The kite's second level has 7 masters instead of
  10.
- **Solver limits.** A reduction of the built-in IBP solver that meets its
  degree or sample bound is repeated with the bounds doubled; the time and
  memory budgets (a day and 4 GiB by default) are what ends a preparation.
- **Threads.** With OpenMP the streams and rows of a Taylor order run in
  parallel once an order costs 5000 multiply-adds; endpoint series solve
  their columns on separate threads; the built-in IBP solver works ahead on
  free threads ([ibp-solver.md](ibp-solver.md)). Results do not depend on
  the thread count. The CLI leaves cores that are busy at start-up alone;
  `OMP_NUM_THREADS` and `--ibp-threads` override.

Options that do not help: a different working precision (384 bits is
needed by the inner levels and more buys nothing), a different contour
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

- The expanded recursion of the double pentagon does not finish in hours:
  levels 2 and 1 have 64 x 33 and 69-master systems with epsilon through 12.
- The exact endpoint series cannot be continued numerically without
  redoing the exact pole-cancelling step on the full series; a timing-only
  run suggests about 1.5 times for the four-loop banana if it could.
- In a sample the endpoint observables and series are still exact
  arithmetic, about 20 s of 150 s for the double pentagon.
