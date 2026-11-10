# Epsilon scheduling and physical-basis search

DiffExp 2.1 schedules Feynman-trick transport by component demand, in both
linear methods (adjoint and factored). The aim is to reduce unnecessary
coefficients while preserving the requested physical outputs and their
uncertainty. This does not guarantee that every recursion level can use the same
epsilon order as the final requested output.

## Demand propagation

For `Y' = A Y`, let `h_i` be the highest coefficient needed from component `i`.
A nonzero entry `A_ij` imposes

```
h_j >= h_i - valuation_epsilon(A_ij).
```

`epsilon_demand.hpp` closes these inequalities, starting from the actual
endpoint and integral consumers. Unreachable components need no coefficients.
A reachable cycle that demands indefinitely higher grades is reported explicitly.
Exact adjacent coordinate maps are multiplied before Laurent expansion so that
cancellation can reduce their demand. This multiplication does not move a map
through an intervening transport.

Each retained row carries its own known upper order. Composition checks this
metadata and rejects reads beyond it; omitted coefficients are never assumed
zero. Numerical zero-containing balls remain uncertain. Only certified exact
zeros can shorten a lower prefix. The scheduler compiles the unequal row ranges
into epsilon-scaled coordinates for the existing transport backend, which retains
its uncertainty propagation and chart checks. Counters report the actual packed
storage, including accumulator tails, rather than just logical output widths.

The recursive child evaluator accepts a common maximum order: the largest of
the per-component windows. The adjoint method uses the same bound. Its local
operator is first assembled with one common low order (a column shift declares
it for all columns); coefficients below a column's structural window are then
replaced by exact zeros, provided the computed ball contains zero, and the
empty leading orders are dropped. A ball that excludes zero is kept, and the
refinement loop then requests the wider window. A fully unequal-width interface
through all recursion levels is not implemented.

## Basis selection

The structural child windows are invariant under every diagonal epsilon
rescaling of the internal basis, and in the examples no candidate basis ever
lowered them. What a different basis does change is the size of the exact
connection, and with it the cost of transport. Selection therefore looks for
the smallest connection that does not widen the child window.

Candidates are cached dotted/numerator reductions, sparse row combinations
(with epsilon weights) and bounded compositions of two changes; pure epsilon
rescalings are not tried. Each candidate is a rank-one update of its base and
is screened only for new parameter poles; it is ranked by the size of its
connection, with no endpoint data. The search is a bounded beam (256 trials,
`basis_scan_seconds`, default 30 s per level), not a proof of optimality.

The best candidate is verified in full (both inverse products and the exact
gauge identity) and installed by pullback: if `Y = D_old U_old = D_new U_new`,
the endpoint functionals at the ordinary cuts are `E_old V` with the exact
`V = D_old^{-1} D_new`, at `V(h)` for the lower and `V(1-h)` for the upper cut.
No endpoint series is rebuilt. It is accepted only if the child window of the
resulting plan is not wider and its gauged connection is smaller.

The chosen coordinates are stored under `basis_cache_directory` (CLI:
`physical-bases` in the cache) and pinned on restart, so checkpoints match. A
search cut by its time budget, or ended by an error, stores nothing.
Selection revision v5 separates these choices from earlier revisions.

## Whole-functional cancellation

For an ordinary-cut functional

```
Phi = Ea Y(a) + integral_a^b B Y dx + Eb Y(b),
```

an exact rational counterterm `C` preserves it under

```
(Ea, B, Eb) -> (Ea-C(a), B-C'-C A, Eb+C(b)).
```

`ft_functional_reduction.hpp` implements the identity, exact certificates,
and a bounded principal-part search that enforces cancellation in the bulk
and at both endpoints together. It also verifies supplied rectangular closed
realizations, including transition, initialization and final-output identities.
It does not discover general invariant submodules automatically.

The pipeline only tries the counterterm search when both endpoint primitives
satisfy their **full exact differential residual**, in addition to the existing
endpoint prescription checks. A spatially truncated primitive or retained
numerical Frobenius operator is insufficient. The default pipeline search uses
polynomial spatial counterterms of degree at most two, at most 64 unknowns and
a one-second cooperative allowance. Failure means no certificate was obtained
within this ansatz, not that a better realization is impossible.

The component scheduler and exact-functional reduction are enabled in
`recursion::NumericalOptions` by default; set `component_epsilon_demands=false`
or `exact_functional_reduction=false` for controlled comparisons.

## Validation and scope

Tests cover weighted dependency cycles, negative Laurent orders, unequal
per-row demands, unknown-tail rejection, uncertain input preservation, and
agreement with the existing complete sunrise evaluator. Exact synthetic tests
check counterterm and realization identities, rejection of merely truncated
endpoint primitives, persisted-basis pinning and the upper-cut pullback.

A smaller matrix or narrower coefficient window is not by itself a speedup
measurement. The implementation reports actual packed slots and recurrence
work so applications can compare complete costs for their own problems.
Automatic invariant-module discovery, joint basis optimization across levels
and general singular-endpoint counterterm discovery are not implemented.
