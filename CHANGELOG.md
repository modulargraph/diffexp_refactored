# Changelog

## 2.2.0

- `diffexp ft --moving-poles remove|exchange` (experimental, default `keep`)
  removes the denominator factors of a level's connection that mix the merge
  parameter with epsilon, by exchanging a master for another integral of the
  cached reduction or for an exact combination (`moving_poles.hpp`).
  `diffexp_master_basis_report` lists the denominator factors of every level.
  On one core the planar double box takes 13.7 s of numerical time instead
  of 19.6; a double-pentagon sample and all epsilon windows are unchanged
  ([docs/moving-poles.md](docs/moving-poles.md)).
- Mathematica `WorkingPrecision` defaults to `Automatic` (twice the
  `AccuracyGoal` plus 30 digits, 100 digits without a goal) instead of 500.
- `--fire-seconds` also bounds one elimination of the built-in IBP solver
  (it was fixed at 20 s).
- Threads. Endpoint series solve their independent columns on one thread per
  hardware thread. With OpenMP (optional, `DIFFEXP_OPENMP`, on by default)
  the adjoint chart expansion runs the independent streams and rows of each
  Taylor order in parallel once an order costs 5000 multiply-adds or more;
  the CLI leaves the cores that are busy at start-up to their work
  (`OMP_NUM_THREADS` overrides). Results are
  bit-identical for any thread count.
- `diffexp --help` is sectioned and lists every default; `diffexp_ft_examples`
  uses the built-in IBP solver when no FIRE path is given.
- The Taylor recurrence of a chart runs on compact numbers (contiguous
  fixed-size limb vectors, limb-aligned exponents, scaled-double radii) with
  the semantics of the ball recurrence, which stays as `expand_balls`,
  reference and fallback. Transport legs of the expanded recursion follow
  the real axis and detour only around real singular points. The Taylor tail
  admitted per chart follows the requested leaf digits. On one core the
  kite takes 5.1 s of numerical time instead of 12.5, the planar double box
  18.7 instead of 40.8, the four-loop banana 49.4 instead of 240.8.
- The order of the exact endpoint series and the overlap point where the
  transport takes over are chosen together: overlap at most 1/64 of the way
  to the nearest other singular point with order 22 (1/16 and 32 before, and
  still for `--epsilon-shift` samples), the same 2^-128 tail for less than
  half the exact endpoint work. `--overlap-ratio R` (configuration key
  `overlap_ratio`) and `--endpoint-order N` set the pair; the example
  configurations no longer fix the order. One core: kite 4.2 s, planar
  double box 16.3 s, four-loop banana 43.7 s.
- An `--epsilon-shift` sample inverts the leading coefficient of each
  endpoint frame by a certified ball inverse instead of an exact elimination
  (`affine_matching::Options::numeric_leading_gauge`): 12 s less per double
  pentagon sample.
- Uncertainty propagation in the recursion's transports. A chart that would
  wrap its inherited uncertainty used to be evaluated through the chart map,
  one expansion per master. The recursion now transports midpoints from the
  first such chart of a leg and follows the size of the discarded radii along
  two probe rows (`--uncertainty probe`, default): an estimate, not an
  enclosure. `--uncertainty map` restores the rigorous map, which remains the
  default of `transport_adjoint_rows` itself. The four-loop banana drops from
  168 s to 54 s of numerical time with an unchanged result; JSON reports the
  mode as `ft_transport.uncertainty`.
- `--epsilon-shift Q` with `--epsilon-order 0` specialises every level above
  the leaf to the dimension d0-2Q: one sample of the four-loop banana takes
  14 s instead of 64 s, one of the double pentagon under four minutes instead
  of 23.
- `--transport-bits N` sets the Taylor tail admitted per chart.
- Each adjoint transport leg reports its rows, masters, epsilon window, chart
  counts and its compile, prepare and recurrence time in the progress output.
- Experimental `diffexp ft --epsilon-shift Q` evaluates the recursion expanded
  around epsilon = Q, and `tools/experiments/epsilon_sampling.py` fits Laurent
  coefficients from such samples (complex values, `--jobs` samples in
  parallel). It reproduces the expanded recursion; it is slower for small
  families and the practical route for the double pentagon; see
  `docs/recursion-cost-experiments.md`, which also records the
  precision experiments and the cost model of a transport leg.
- A requested recursive anchor (`preparation.anchors`) that is singular in the
  master basis the reduction selected moves to the simplest regular rational
  within 3/20 of it that has no other real singular point of the level in
  between, so a supplied causal prescription keeps its meaning; the progress
  output says so. It used to abort. The family-name anchor table is gone:
  ladders live in the family file, and the double-pentagon file carries the
  preparation budgets it needs.
- The built-in IBP solver reconstructs rational coefficients with up to 32
  primes (was 8). A preparation that runs out of time names the three budgets
  in force and how to raise them.
- Removed the inherited prepared-input compatibility runtime (`src/kernel`,
  `tests/kernel`, the `kernel` command and `DIFFEXP_BUILD_KERNEL_RUNTIME`).
  The exact `series` command now uses a direct recurrence; it also corrects
  the highest retained epsilon coefficient of a Laurent window when the
  cleared denominator depends on epsilon.
- Feynman-trick pipeline cleanup. Removed: the FT spectral and ultraspherical
  transport backends (`ft_spectral*.hpp`, `ft_ultraspherical*.hpp`,
  `diagonal_log_gauge.hpp`; Taylor charts were as fast or faster on every
  measured leg), the legacy "values" method, the on-disk endpoint-series cache
  (`cached_affine.hpp`; recomputing is faster than verifying a stored series),
  the endpoint principal-part plumbing and the tools-only epsilon-basis
  optimiser of the evaluator, and four no-op endpoint options
  (`univariate_epsilon_recurrence`, `univariate_epsilon_projection`,
  `cleared_epsilon_projection`, `finite_lag_projection`).
- CLI and configuration: `--ft-transport`, `--ft-transport-digits` and
  `--method values` (config: `numerical.transport` other than
  `"taylor"`/`"auto"`, `numerical.transport_digits`, `numerical.method:
  "values"`) are rejected with a "removed" message. JSON: `ft_transport` is now
  `{method, direct_endpoints}`; `timings.spectral_seconds`,
  `basis_selection.scans_skipped` and `basis_selection.endpoint_series_reused`
  are gone. The cache directory holds `physical-bases` (was
  `affine-series/physical-bases-v1`) and `ordinary-transport`;
  `--no-numerical-cache` disables these two.
- The physical-basis scan runs at every level (the two skip rules are gone: a
  smaller connection pays for the scan there as well), no longer tries pure
  epsilon rescalings, and builds candidates by rank-one updates.
- Removed tools: `profile_accepted_chart`, `import_centered_checkpoint`,
  `benchmark_circuit_source`, `benchmark_adjoint_dot`, `benchmark_affine_domain`,
  `verify_affine_cache`, `diagnose_endpoint_sectors`, and the
  `tools/performance` experiments built on removed kernels. The series
  verifier of the former cache is kept as a test oracle
  (`tests/affine_series_verifier.hpp`).
- Feynman-trick charts (`rational_line`, the adjoint executors) estimate their
  Taylor tail and halve a step whose polynomial has not converged at the
  designed rate. Adjoint checkpoint
  identities advance, so earlier numerical checkpoints are not reused.
- The physical-basis scan ranks candidates by the size of their exact
  connection (the child window may not widen), a chosen basis is installed by
  exact pullback instead of rebuilding endpoint series, endpoint operators are
  prepared once per level and extended in epsilon, and the adjoint path
  requests the per-column structural child window. `epsilon_low` can rise
  where a leading order was an exact structural zero. Basis cache revision v5.
- Endpoint constants and divergent-mode constraints fail with a clear error
  when the endpoint x order does not reach the fixed sector, instead of being
  silently truncated.
- Built-in IBP preparation plans the derivative closure first, validates
  held-out points on the pilot's dependency subset of original equations with
  one complete elimination per attempt (`--ibp-validation strict` restores a
  complete elimination at every point), skips and records exceptional sample
  points, and continues primal-only when a trace exceeds its budget.
- Generic transport accumulates omitted-tail estimates beside the balls for
  systems with epsilon-zero feedback. Mathematica `ToPiecewise` evaluates saved
  Taylor segments on undeformed paths instead of replaying the transport.
- CLI: numerical options are validated before IBP preparation, ignored
  `--fire` paths and reused level closures are reported, and the JSON report gains per-phase timings, `ibp_cache` and
  `basis_selection`.

- Generic transport reuses a chart's Taylor coefficients when an accuracy
  refinement shortens its step, predicts the shorter step from the truncation
  estimate, and no longer charges a missed truncation target to the finite-lag
  recurrence.
- Rational physical-basis entries use an O(N deg) product recurrence, and the
  uncertainty map of a centered chart uses a short mantissa.
- A pole on the transport path or at its endpoint is reported at once, with
  its path parameter.
- Fix the kite example with the built-in IBP solver: per-variable degree boxes
  have their own cap and are no longer cut by a total degree; the CLI enables
  coordinate slices; budget errors name the exhausted limit.
- Fix "recursive local operator lookahead invariant failed" after exact-zero
  padding was trimmed from a child value but not from its stored map.
- Reject scalar-leaf power divergences that epsilon does not regulate, fix the
  spectral coupled-block node budget, and repair the stale batched expectation
  in the conditioned-adjoint test.

- Retry deficient reduced-precision uncertainty maps with the compact recurrence
  at full precision before invoking legacy reference columns.
- Avoid unnecessary epsilon refinements caused by exact-zero padding while
  preserving unknown upper coefficients and uncertain intervals.
- Include coefficient values and radii in endpoint constraint failures, and
  report homogeneous fallback and precision-retry counters on failures.

- Schedule factored Feynman-trick transport by per-component epsilon demand and
  compose exact adjacent maps before truncation.
- Rank physical bases by coefficient demand before expression size, reuse
  endpoint maps during screening, and persist verified choices across restarts.
- Add exact whole-functional counterterm and rectangular-realization certificates.
- Compile compact rational recurrences with shared chart preparation, bounded
  auxiliary histories, and batched FLINT dot products.
- Store homogeneous uncertainty maps sparsely, omit exact-input columns, and
  shorten map windows using the support of inherited uncertainty.
- Add algebraic Lean contracts and tests for recurrence dependencies, finite
  coefficient windows, and enclosure selection.

See [the transport backend design](docs/backend-review-implementation.md) and
[epsilon scheduling](docs/epsilon-demand-design.md) for the algorithms and scope.

## 2.1.1

- Add generic JSON-defined Feynman families, including editable configurations
  for every supplied example and arbitrary requested integral powers.
- Add native rational/canonical/algebraic transport and partial Frobenius
  boundary matching, exposed through the original Mathematica configuration,
  boundary, transport and plotting workflow using RunProcess.
- Preserve boundary precision, report preparation/numerical/wrapper timings,
  and retain explicit error estimates and omitted-tail status.
- Accept explicit algebraic basis prefactors for principal endpoint
  normalization after root continuation, including through the Mathematica wrapper.
- Use guarded grouped-dot rational charts with the scalar fallback, and share
  canonical letter expansions and source convolutions across matrix entries.
- Retain carried boundary uncertainty through local homogeneous maps where
  feedback makes direct interval recurrence poorly conditioned.
- Reuse nearby physical samples for single-coordinate Mathematica plotting.
- Include the original equal-banana reference data needed by a fresh checkout.


## 2.1.0

- Publish the standalone C++20 mathematical backend.
- Retain series recurrences, regular-singular matching, Feynman-trick recursion,
  native FIRE/FIRE7p integration and verified durable checkpoints.
- Preserve shared-boundary correlations across recursion levels; both four-loop
  banana examples pass independent finite-part comparisons.
- Add a small Mathematica `RunProcess` wrapper and JSON CLI output. Remove the
  former Mathematica/LibraryLink application from the current repository tree.
- Install as CMake package `DiffExp`, target `DiffExp::core`, namespace `diffexp`,
  executable `diffexp`, version 2.1.0.
- Keep all original example groups as opt-in data comparisons. Full Henn
  108-component FT reconstruction remains experimental and frozen.

This is an API break from 2.0. Existing mathematical cache formats remain readable. The old application remains available in Git history.
