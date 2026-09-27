# Epsilon scheduling and physical-basis search

DiffExp 2.1 now schedules factored Feynman-trick transport by component demand.
The aim is to reduce unnecessary coefficients while preserving the requested
physical outputs and their uncertainty. This does not guarantee that every recursion level can use the same
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

The recursive child evaluator still accepts a common maximum order. Its
per-component ledger is available, but a fully unequal-width interface through
all recursion levels is not yet implemented. This is one remaining source of work.

## Basis search

Candidates include cached dotted/numerator reductions, epsilon prefactors,
sparse combinations and bounded compositions of changes. Every proposed
coordinate system must pass exact inverse and differential-gauge identities.
The search is a bounded beam, not an exhaustive proof of optimality.

The production score prioritizes maximum child lookahead, then the sum of
positive component lookaheads, then connection expression size. A larger matrix
can therefore win if it needs fewer coefficients. This is a per-level proxy;
there is not yet a joint optimizer of the full recursion or adjacent-level bases.

Most candidates reuse the original prepared endpoint maps. If `Y = D_old U_old`
and `Y = D_new U_new`, their ordinary-cut maps are pulled back by
`V = D_old^{-1} D_new`. The lower cut uses `V(h)` and the upper cut uses
`V(1-h)`. Required epsilon lookahead is supplied before composition. At most
eight improving candidates receive full endpoint reconstruction and contour
validation, stopping at the first certified improvement. Final validation has
a separate cooperative allowance capped at 30 s, so the next expensive rebuild
is not started after that allowance. No trial numerical transport is needed
to rank a basis.

When automatic scanning is enabled, its default discovery allowance is 30 s.
This is a cooperative allowance: an in-progress exact operation can exceed it,
and final candidate validation is additional work. Selected coordinates remain
pinned on restart and are revalidated. Selection revision v3 separates these
choices from earlier expression-first and full-rebuild search decisions.

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
endpoint primitives, persisted-basis pinning and upper-cut coordinate screening.

A smaller matrix or narrower coefficient window is not by itself a speedup
measurement. The implementation reports actual packed slots and recurrence
work so applications can compare complete costs for their own problems.
Automatic invariant-module discovery, joint basis optimization across levels
and general singular-endpoint counterterm discovery are not implemented.
