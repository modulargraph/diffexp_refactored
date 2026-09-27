# Transport backend design

DiffExp 2.1 uses compiled rational recurrences, sparse uncertainty propagation
and explicit coefficient-demand contracts for ordinary Feynman-trick transport.
FLINT provides exact algebra and Arb/ACB provides numerical ball arithmetic.

## Compiled recurrence

`rational_circuit.hpp` separates three lifetimes:

- A compiled plan stores typed polynomial blocks, normalized Gaussian-rational
  denominators and shared streams identified by source column and denominator.
  Each stream represents `W = Y/Q`; numerator blocks consume that stream.
- A prepared chart fixes the center, precision, Taylor order and epsilon window.
  Shifted coefficients, denominator pivots, inverse pivots and ring lengths can
  be shared by independent right-hand sides and homogeneous-map columns.
- An executor holds independent ball payloads and bounded auxiliary histories.
  Unit denominators alias source history; reachability removes inactive rows
  and streams. Work and live-cell counters describe the actual execution.

Observable rows are batched against a complete workspace budget, including
auxiliary storage. The backend does not lower precision or order to fit the
budget. Pooled FLINT `acb_dot` calls accumulate polynomial convolutions before
rounding, with scratch storage included in the workspace accounting.

## Uncertainty propagation

When inherited uncertainty already exceeds the arithmetic reserve, the solver
can propagate a guarded midpoint plus a homogeneous uncertainty map directly.
This avoids a redundant uncertain-input recurrence. Existing domain,
conditioning, precision-retry and numerical acceptance checks remain active.

Homogeneous maps are stored sparsely. Structurally zero connection columns have
exact identity action. Exact-input columns need no uncertainty-map solve, and
reachability eliminates components that cannot affect a retained output.
Independent columns share prepared coefficients and respect the worker budget.

The incoming uncertainty support is computed across all observable batches
before any batch advances. If all incoming coefficients below index `L` vanish
exactly, an output window of width `W` needs only `W-L` map coefficients for that
column. Application checks every omitted uncertainty coefficient. Unknown data
and small zero-containing balls are never treated as exact zeros.

[Lean proofs](proofs/UncertaintyWindowVerification.md) establish the relevant
finite-convolution and affine-decomposition identities. They do not formally
verify the C++ implementation, ACB rounding or omitted analytic tails.

## Epsilon demand and basis selection

[Component scheduling](epsilon-demand-design.md) propagates required orders
backward through the differential system and its consumers. Exact adjacent maps
are composed before truncation to retain cancellations. Known upper orders
remain attached to unequal-width expressions until composition is complete.

Physical-basis candidates are scored by coefficient demand before expression
size. Screening reuses prepared endpoint maps; only a bounded improving
shortlist receives full reconstruction and contour checks. Accepted coordinate
systems are persisted and revalidated, rather than silently changed on restart.
Exact whole-functional counterterms are accepted only with the required
endpoint and differential identities.

## Checkpoints and numerical policy

Checkpoint identities include the exact equations, initial balls, path, retained
orders, precision and numerical policy. A changed representation cannot silently
reuse an incompatible numerical continuation. The explicit checkpoint importer
validates immutable input snapshots and records transformations between
supported numerical policies.

The default compact policy enables `rational_circuit_recurrence`,
`compact_centered_recovery`, `compact_centered_only` and `circuit_grouped_dot`.
Their corresponding options allow controlled comparisons. Reference recurrences
remain available for verification and conditioning recovery.

Exact rational extraction, constant construction, integer coefficients and
valuations use FLINT objects directly. Rational coefficients enter Arb through
`arb_set_fmpq`, avoiding decimal formatting and reparsing in that conversion.
Serialization and interning still use strings where appropriate.

## Verification

Tests cover rational recurrence agreement, exact support and Laurent-window
handling, multiple observable batches, incoming uncertainty, checkpoint identity,
physical-basis persistence and complete small-family recursion. Reproduction
requires only the packaged tests and their documented optional dependencies;
no local run archives are required.

General recursive outputs retain their explicit omitted-tail status. Algebraic
identities, retained-coefficient agreement and analytic error bounds are
separate guarantees.
