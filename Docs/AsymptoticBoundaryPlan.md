# DiffExp 2 Singular/Asymptotic Boundary Plan

## Status and objective

DiffExp 2 already transports exact tagged sectors with the requested
`+/- i0` convention, and `PrepareLaurentBoundary` preserves Laurent-valued
point data at a regular anchor. `SegmentLine` still rejects a singular
starting point. This document specifies the missing DiffExp 1 asymptotic
boundary mode without weakening the existing branch or epsilon-window
contracts.

The objective is a typed constructor and a singular-start transport path for
exact local data of the form

```text
sum c[k,n,p] eps^k x^(a+b eps+n) Log[x]^p,
```

with all exponents, log powers, Laurent bounds, and prescription data retained
structurally.  It is not a generic replacement-rule facility.

## Public object and constructor

Add

```text
PrepareAsymptoticBoundary[expressions, localParameter, options]
```

returning a validated `DiffExp2.AsymptoticBoundary/v1` record containing:

- the exact `LocalSolution` at the singular chart;
- its honest epsilon and Taylor windows;
- the singular physical point and affine chart map;
- pulled-back prescription factors and their provenance;
- the derived chart imaginary sign;
- a construction certificate recording the parsed source terms and the
  exact reconstruction residual.

The constructor accepts either prescriptions already expressed in the local
parameter or a physical line map from which every physical factor can be
pulled back.  Rewriting a bare `Log[x]` globally is forbidden: simultaneous
physical factors can require different orientations, and only the pulled-back
factor/multiplicity/leading-coefficient records determine the local sign.

## Exact construction algorithm

1. Canonicalize the epsilon symbol and pin variables through `Config`.
2. Pull each physical prescription through the line map.
3. At the singular start, factor its leading local monomial
   `c x^m (1+O(x))`; record `m`, `Sign[c]`, and provenance.
4. Reuse `SectorSeries` `ChartImSign` as the sole owner of the local `+/- i0`
   decision.  Odd multiplicities constrain the sign; even multiplicities do
   not.  Conflicts are loud.
5. Parse each boundary summand into exact `(a,b,p)` tags plus a finite
   epsilon-by-Taylor coefficient tensor.  Accept rational/algebraic exact
   `a`, exact `b`, and nonnegative integer `p`; reject inexact tags.
6. Convert ordinary `Log[x]^p` to the package normalization
   `(eps Log[x])^p/p!` by the exact epsilon shift and factorial weight.  This
   can deepen the lower Laurent window and must never be hidden by padding.
7. Merge only sectors allowed by `CanonicalizeLocalSolution`, retaining
   honest lower and complete-upper windows.
8. Reconstruct the accepted expression from the tagged object on a symbolic
   positive arm and certify an exact zero residual through the stored order.
9. Validate the final `LocalSolution`; no numerical exponent fitting or
   `SeriesData`-dependent inference is permitted.

## Branch convention

For a negative local coordinate `x=-r`, the constructor and marcher use

```text
+i0: Log[x] = Log[r] + i pi,  x^alpha = r^alpha exp(+i pi alpha)
-i0: Log[x] = Log[r] - i pi,  x^alpha = r^alpha exp(-i pi alpha).
```

Thus the lower-rim value is Mathematica's principal value with
`Log[x] -> Log[x]-2 pi i` and
`x^alpha -> exp(-2 pi i alpha) x^alpha`.  The rule applies to the complete
`alpha=a+b eps` and to the triangular mixing of every log-chain sector.

## Singular-start planning and marching

`SegmentLine` may accept a singular `from` only when paired with a validated
`AsymptoticBoundary` whose center, variable, prescriptions, and system
dimension agree with the planned first chart.

The marcher must then:

1. install the supplied local object directly as the first chart solution;
2. never evaluate it at local coordinate zero;
3. choose the first handoff at the ordinary certified nonzero overlap point;
4. evaluate that handoff on the correct outgoing arm;
5. continue with the existing matching, crossing, residual, and error-probe
   machinery;
6. record that the first chart was boundary-seeded in checkpoints and cache
   keys.

Reversing a line is not implemented by conjugating numbers.  It constructs a
new plan, re-derives the local side from the same physical prescriptions, and
uses the exact crossing operator where required.

## Error contract

The new path is loud when:

- the singular start has no validated typed boundary;
- a material multivalued term lacks a derivable sign;
- pulled-back factors conflict or have an indeterminate real leading sign;
- an expression is not a finite sum of supported tagged monomials;
- a requested coefficient lies outside the honest epsilon/Taylor window;
- the first overlap point is outside either convergence disk;
- reconstruction or ODE residual certification fails.

`AbortOnAnalyticContinuationFail -> False` must not choose a sheet.  If a
future soft mode is retained, it must return an explicitly flagged,
one-sided-only result that public numerical accessors refuse to present as a
complete transport.

## Acceptance gates

Minimum exact tests:

1. `Log[x]`, half-integer power, general rational power, and `x^(b eps)`
   starts on both `+i0` and `-i0` sheets.
2. Odd multiplicities 1 and 3, an even multiplicity, reversed polynomial
   orientation, consistent simultaneous factors, and a loud conflict.
3. Singular-start transport away from the point and a round trip back to a
   certified nonzero overlap point in both directions.
4. A resonant log chain with the exact lower Laurent shifts from
   `(eps Log[x])^p/p!` normalization.
5. A boundary containing several sectors and components with unequal honest
   lower windows.
6. Missing-prescription, malformed-expression, insufficient-window, and
   incompatible-plan failures.
7. Parity against a small DiffExp 1 boundary example on both sheets.
8. Wolfram/C++ parity after the boundary has been converted to the common
   rational recurrence representation.

## Non-goals for the first slice

- arbitrary Mathematica branch-function rewriting;
- starting from an irregular singularity unsupported by the recurrence;
- inferring prescriptions from numerical complex samples;
- square-root matrices themselves (covered by `AlgebraicSquareRootPlan.md`);
- changing the established real PV/Hadamard convention for integer-power
  meromorphic interior integration.
