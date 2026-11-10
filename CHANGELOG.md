# Changelog

## Unreleased

- Feynman-trick charts (`rational_line`, the adjoint executors) estimate their
  Taylor tail and halve a step whose polynomial has not converged at the
  designed rate; `max_truncation_halvings` bounds it. Adjoint checkpoint
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
- A stored endpoint-series column rejected by the verifier is discarded and
  recomputed by default; fail-closed remains selectable.
- Built-in IBP preparation plans the derivative closure first, validates
  held-out points on the pilot's dependency subset of original equations with
  one complete elimination per attempt (`--ibp-validation strict` restores a
  complete elimination at every point), skips and records exceptional sample
  points, and continues primal-only when a trace exceeds its budget.
- Spectral FT transport no longer subdivides on inherited input uncertainty,
  retries an exhausted coupled-block level by subdivision, honours
  `matrix_bits`, and pivots in the ultraspherical factorization. The two FT
  spectral backends share one driver.
- Generic transport accumulates omitted-tail estimates beside the balls for
  systems with epsilon-zero feedback. Mathematica `ToPiecewise` evaluates saved
  Taylor segments on undeformed paths instead of replaying the transport.
- CLI: numerical options are validated before IBP preparation, `--ft-transport
  ultraspherical` is accepted, ignored `--fire` paths and reused level closures
  are reported, and the JSON report gains per-phase timings, `ibp_cache` and
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
