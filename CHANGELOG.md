# Changelog

## 2.2.1

- Cancellation between the levels. The exact row gauge of an endpoint
  frame takes a pivot within 2^-20 of the largest admissible one instead of
  the first entry that is not zero: the expanded run of a three-loop family
  with nine propagators returned radii of 1e+1301 at any working precision
  and now has 30 digits at 384 bits. Weighted rows drop the epsilon orders
  that cancelled to rounding (`NumericalOptions::clear_cancelled_orders`,
  `cleared_epsilon_orders` in the JSON report): unchanged values with 17 to
  31 per cent fewer products in the expanded runs of the bananas and of
  that family. A fit no longer fails on a sample that lowered its admitted
  Taylor tail, and later samples start at that tail (`--tail-floor`).
  [performance.md](docs/performance.md) has the measurements of where the
  cancellation comes from.
- Numerical phase, counted work. A moving pole is kept when its removal
  would leave the level's connection more than four times larger
  (`moving_poles::Options::max_growth`), and the outcome of the moving-pole
  search is recorded beside the exact cache (`moving-poles/`) and applied
  again instead of being searched for by every run and every sample. The
  gauge of an endpoint to Fuchsian form applies a shear entry by entry
  instead of through an exact inverse and three dense products. A transport
  leg reads its poles from the compiled circuit instead of isolating the
  roots of every denominator again. An endpoint series converts its
  coefficients for projection once. `timings` lists
  `endpoint_gauge_seconds`, `endpoint_frame_seconds` and
  `endpoint_observable_seconds`; the phases now add up to
  `numerical_seconds`. Results of the default runs are unchanged.
- Working precision. Without `--precision-bits` (or `working_bits` in the
  family file) a run starts at 384 bits and computes its numbers again with
  more when rounding leaves a coefficient fewer digits than `--leaf-digits`
  asks for (`precision` in the JSON report; at most 4096 bits). The samples
  of a fit start at the precision an earlier sample needed
  (`--precision-floor`).
- Built-in IBP solver, preparation. Seeds start small and grow on demand:
  the default numerator degree is 1 (was 2; `--ibp-numerators`,
  `preparation.ibp_numerators`), and a requested integral no longer raises
  the seeds of every sector in advance (`--ibp-seeds demand|sector|
  subsectors|all`, `ibp_solver::Options::target_seeds`; `all` with
  `--ibp-numerators 2` is the previous behaviour). An integral left
  unreduced outside the seeds of its sector widens that sector alone, as
  before, now up to 8 times per batch. A complete elimination continues
  from the pivots of the previous one at the same point instead of being
  repeated: closure-planning passes share their two points, the same
  equations are not eliminated again and wider seeds insert only the new
  equations (`continued_solves`, `continued_equations`, `reused_solves` in
  `ibp_statistics`). Samples pass from the solver to the reconstruction as
  machine words and reconstructed coefficients are built from integer
  terms; nothing is written to text and parsed back inside one process.
  The level closures of kite, banana, banana4, double box, family F and
  the double pentagon are unchanged (same content identifiers). Caches
  written by earlier versions are not reused for the solver's samples and
  closure plans; verified level closures are.
- Fixed: a chart centre of few bits deep below the binary point (the first
  charts next to an endpoint) had the coefficient polynomials of degree
  above ten rounded while they were shifted to it, at a precision that did
  not follow the working precision. Their relative radius of about 2^-300
  was harmless at 384 bits and exhausted the conditioning halvings from
  about 640 bits on ("adjoint arithmetic-conditioning step budget
  exhausted", four-loop banana with `--precision-bits 640
  --transport-bits 180`). The shift is exact again.

- The scalar leaf of a recursion around a shifted epsilon is expanded at
  the shifted dimension itself (`CertifiedDeepestBetaOptions::epsilon_shift`,
  `feynman::tadpole_branch` with a shift), through the requested order,
  instead of summing some twenty orders of its series at epsilon=0 and
  re-expanding them. The leaf of an epsilon sample costs 15 to 50 times
  less (double box 1.07 s -> 0.04 s, four-loop banana 1.8 s -> 0.05 s),
  half of the numerical time of a sample of a small family. A shift on a
  pole of the gamma function or resonant with an endpoint power falls back
  to the series (`NumericalOptions::direct_shifted_leaf`).

- Endpoint operators of a recursion at a numeric dimension: the leading
  coefficients of a functional composed with the ball inverse that vanish
  exactly are exact zeros again. They are proved like the zeros of the
  inverse, by evaluating the product modulo two primes at the exact point
  (`affine_operator::impose_exact_zeros`, `Options::operator_exact_zeros`);
  no size is compared with anything. Rows built on such an operator no
  longer carry columns that are not theirs: the double pentagon sample by
  rows has 1201 such entries and 40% fewer products, Higgs+jet family F 62.
  The report has `basis_selection.operator_exact_zeros`.
- Adjoint transport finds the singular points of a connection once and maps
  them to every leg of the path, instead of once per leg from the pulled-back
  coefficients.

- Reported radii contain an estimate of the Taylor tails the transport
  admits. What a chart cuts off is estimated from its last two retained
  terms (a geometric tail, with its sign) and added to the probe rows, which
  the flow carries to the end of the path. The probes are returned as rows
  (`AdjointOptions::error`) and every further linear map is applied to them
  as to the result (`LaurentBoundary::error`, eight times of which
  `diffexp ft` adds to the printed radii): carried as a radius the estimate
  would not cancel where the result does, and a recursion of several levels
  would report it larger by the product of their condition numbers.
  Weighted accumulation carries the error rows from level to level. An
  evaluation by rows estimates its error by accumulating one weighted
  combination of all requested integrals more
  (`NumericalOptions::rows_error_estimate`, selected by `diffexp ft`).
  The probes are a batch of their own, at half the Taylor order and a third
  of the working precision; the progress line of a transport counts their
  products apart. The factored method is unchanged and has no estimate.
  The fixed widening of `diffexp ft` is one tail relative to the largest
  coefficient of an integral (it was 2^6 tails).

- Transport: the endpoint overlap is chosen per endpoint (a singular point
  next to one end no longer costs charts at the other); in probe mode a
  chart's radii are its own rounding (defect radii of the midpoint
  recurrence), so steps are no longer halved by inflated interval bounds;
  a single transported row is followed by one probe row instead of two.
- One weighted row per level (`diffexp ft --accumulation weighted|rows`,
  default `weighted`; `NumericalOptions::accumulation`, whose library
  default is `rows`). The product of the levels' local operators is formed
  from the top: a level transports the combination of its rows weighted by
  what the level above passes down, one row per requested integral of the
  family, and hands the weights on its masters to its child; the leaf
  values enter last. Counted in multiplications per Taylor order the inner
  levels cost 4 to 6 times less (docs/performance.md); not yet timed.
  Covered: the adjoint method with at most four requested integrals
  (`weighted_max_targets`), including endpoint-boundary levels, internal
  bases and endpoint constraints, which are carried as further rows and
  checked through the requested epsilon order. Other requests are
  evaluated by rows as before. The report has `ft_transport.accumulation`.

- The Feynman i0 is derived (`causal.hpp`): `diffexp ft` evaluates families
  in physical regions from their definition, with no `causal` object. Every
  real singular point of a level is passed on the side
  `Im x = -sign(dF/dx)` of its Landau pinch of the level's Symanzik `F`,
  point by point (two cuts of one level may need opposite sides); points
  without a pinch are free. The scalar leaf passes simple zeros of `F(x)`
  inside `(0,1)` off the axis on the side `-sign(F')`. Decisions are made in
  ball arithmetic: monotonicity of `F` in the level parameter for any number
  of loops, the face-by-face Landau analysis for one-loop levels. A point
  that is not decided (a multi-loop level above a threshold where `F` is not
  monotone, a one-loop pinch on the boundary of a face) stops the run with a
  message naming the level and the point; a repeated zero of the leaf's `F`
  (a family exactly at threshold) stops it as a pinched contour. `causal`
  remains as an override: with `level_signs` it states one side per level,
  with `points` the side of single points, and without either it only sets
  `f_rim`. When points are treated one by one, no detour encloses a complex
  singular point of its level. Families with positive `U` and `F` keep
  their contour. The report has `causal.mode` and the oriented and free
  points.
  Examples `bubble_timelike`, `bubble_above_threshold`, `triangle_physical`
  and `box_physical`; test `cli.causal`.
- `diffexp ft --epsilon-strategy sample` (family key
  `numerical.epsilon_strategy`): the Laurent coefficients from samples at
  shifted dimensions. The recursion is prepared once, the samples run as
  parallel processes (`--samples`, `--sample-step`, `--sample-jobs`), the
  pole order is detected (`--sample-poles` fixes it), and the fit reports
  each coefficient with an estimated error. `epsilon_sampling.hpp` holds the
  fit; `tools/experiments/epsilon_sampling.py` is removed.
- Endpoint boundary (`--boundary anchor|endpoint|auto`, default `auto`;
  `recursion::Options::boundary`, whose library default is `anchor`): a
  level may take its boundary data at x = 0 or x = 1 of its merge parameter
  instead of at an interior anchor. The level below is then the family with
  one of the two merged propagators removed, at its original kinematics. An
  endpoint qualifies if the family's U and F have positive monomials, the
  two lines share one loop momentum, the line that remains is massive, no
  master's numerator degenerates there, and the level's exact endpoint
  series shows that the masters' values fix its analytic solution
  (`endpoint_boundary.hpp`: `kinematics`, `examine`, `determine`,
  `taylor_map`). Every other level, and every level under explicit
  `anchors`, stated causal sides or `--henn-basis`, keeps its anchor;
  families without a qualifying level are computed as before. The kite's
  third level has 2 masters instead of 3; the Higgs+jet family F has levels
  of 49, 12, 8, 7, 3 masters instead of 49, 19, 13, 8, 3. The report has a
  `ladder` array with the masters and the boundary of each level. Tests
  `endpoint_boundary` and `cli.endpoint_boundary`.
- Reported radii of an uncertified result include the Taylor tail admitted
  per chart (2^6 charts' worth, relative to the largest coefficient), so a
  sample asked for few leaf digits no longer reports more digits than it has.
- Merge selection (`--merge-selection auto|first`, default `auto`): every
  level merges the pair of propagators whose merged graph has the fewest
  master integrals. The count is the number of critical points of U + F,
  computed with Singular in seconds for the pairs on one loop momentum
  (`critical_points.hpp`); without Singular, IBP counts
  (`ibp_solver::Session::count_masters`) decide for large levels.
  Higgs+jet family F: levels of 75, 47, 21 masters become 49, 19, 13; the
  kite's second level 7 instead of 10.
- The executables are compiled for the host CPU by default
  (`DIFFEXP_NATIVE`, with floating-point contraction off so results do not
  depend on the CPU); `tools/build-tuned-deps.sh PREFIX` builds GMP and FLINT
  for the host CPU.
- The built-in IBP solver no longer stops at fixed sizes: degree and sample
  bounds double automatically when met, seed and term counts are bounded by
  memory only, and the default time budgets are an hour per elimination and
  a day per level and in total.
- The built-in IBP solver takes sectors proved scaleless by the scaling
  criterion on U + F as zero sectors (`scaleless_sectors.hpp`), not only
  those with a free loop momentum, and no longer accepts as a master an
  integral left unreduced outside the seeds of its sector: the seeds of that
  sector are widened and the elimination repeated
  (`seed_widenings` in `ibp_statistics`). A level whose requests carry two
  dots on one line used to close on ever higher dots. Planar double box with
  its second level merging the two outer lines of the second loop: 9
  masters in 2 s, where the plan did not settle in 15 minutes; that merge is
  now the one selected. First level of the double box: 18910 equations
  instead of 55010. Solver caches written before this change are not reused.
- With `--transport-bits` below 87 the endpoint-constraint check uses half
  the bits of the admitted Taylor tail as its tolerance.
- `--ibp-max-degree N` raises the degree of reduction coefficients the
  built-in solver reconstructs (default 16), with its per-variable degree and
  sample budgets.
- Singularity polynomials are factored into irreducibles before their roots
  are isolated; a level with a denominator of degree in the hundreds no
  longer stalls there.
- Documentation describes the current state only; timings are near the top
  of the README.

## 2.2.0

- `diffexp ft` removes the denominator factors of a level's connection that
  mix the merge parameter with epsilon, by exchanging a master for another
  integral of the cached reduction or for an exact combination
  (`moving_poles.hpp`; `--moving-poles keep|remove|exchange`, default
  `remove`, `keep` with `--epsilon-shift`). `diffexp_master_basis_report`
  lists the denominator factors of every level. The gain is a smaller
  connection, not a narrower epsilon window: on one core the planar double
  box takes 12.3 s of numerical time instead of 16.3; families without such
  factors, a double-pentagon sample and all epsilon windows are unchanged
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
- The built-in IBP solver works ahead on several threads (`--ibp-threads N`,
  default the hardware threads idle at start-up; `ibp_solver::Options::threads`).
  Trace replays of announced sample points run in batches, and the complete
  eliminations of closure planning and validation, whose points are fixed in
  advance, start early on free threads. Reductions, sample points, counts and
  cached records are identical for every thread count; `ibp_statistics`
  reports `threads`, `batched_replays`, `early_solves` and
  `early_solves_used`. Cold preparation on four threads against one: planar
  double box 3.4 s against 5.2, four-loop banana 15.3 s against 33.7, double
  pentagon 158 s against 298 (peak memory 1.24 GB against 0.85). See
  `docs/ibp-solver.md`.
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
  (`affine_matching::Options::numeric_leading_gauge`), with the exact zeros
  of the inverse restored by a modular check: 13 s less per double pentagon
  sample.
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
  `docs/performance.md`, which also records the
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
