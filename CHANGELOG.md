# Changelog

## Unreleased

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
