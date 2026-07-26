# Changelog

## Unreleased

- Reproduced every example notebook in the pinned public DiffExp revision:
  the Henn and one-mass canonical pentagons, multiple polylogarithms, and
  both equal- and unequal-mass Banana endpoints, with hash-pinned inputs,
  saved-reference checks, and opt-in timing gates.
- Parallelized explicitly bounded terminal endpoint batches while preserving
  atomic, request-ordered publication. The 15-component unequal-mass Banana
  terminal projection fell from about 494 to 208 seconds on the development
  machine.
- Bound deferred native value-handoff accuracy to the observable run's
  matching target, preventing low-digit transports from silently performing
  high-precision tail-certificate retries.
- Upgraded Feynman-trick preparation to FIRE 7.1 with finite-field Zippel
  sampling, multiprime MPI execution, exact rational reconstruction, and an
  explicitly selectable Classical parity backend.
- Added source/runtime-bound modular caches, atomic completion markers and
  locks, exact-artifact validation, process-tree timeouts, independent basis
  probes, and strict table-coverage checks.
- Expanded basis discovery with first-shell dotted and numerator seeds, fixing
  request-dependent sunrise masters; real bubble and sunrise reductions now
  agree exactly between Classical and modular FIRE7.

## 2.0.0 — 2026-07-11

- Replaced the strategy/fallback stack with one strict finite-width recurrence
  solver for ordinary, regular-singular, resonant, logarithmic, fractional,
  and inhomogeneous sectors.
- Added the C++20/FLINT recurrence backend and made it the release default,
  while retaining an explicit Wolfram reference backend for parity checks.
- Preserved exact local `x^(a+b eps) Log[x]^p` sectors through transport,
  endpoint limits, and analytically regularized line integration.
- Restored the classic coupled segmentation geometry, including complex-root
  projections and matching at `+1/DivisionOrder` and
  `-1/DivisionOrder` in adjacent regular charts.
- Added a stable `DiffExp2`` umbrella API, named result schemas, closed-form
  regular boundary preparation, inspectable segments, and piecewise
  evaluation.
- Added a typed Feynman-trick facade with prepared FIRE caches, atomic ladder
  checkpoints, strict subprocess settings, and safe invocation-local IBP
  batching.
- Added unequal-mass banana comparisons, independent Bessel oracles, and
  focused Wolfram/C++ parity, analytic-regulator, and edge-case tests.

This is a source-incompatible major release. See
[Docs/Migration.md](Docs/Migration.md) for the DiffExp 1 mapping and current
limitations.
