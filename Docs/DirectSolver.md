# Direct Solver

The direct solver transports a boundary vector for

```text
d f(x,eps) / dx = A(x,eps) f(x,eps)
```

where `A` is exact and rational in the line parameter and epsilon.  Local
solutions are built by finite-width coefficient recurrences and matched along
a certified chain of affine charts.

## 1. Load the package and configuration

```mathematica
Get["/absolute/path/to/DiffExp2/DiffExp2/DiffExp2.m"];

DiffExp2`Config`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 200,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 4,
  "DivisionOrder" -> 3,
  "RadiusOfConvergence" -> 1,
  "DeltaPrescriptions" -> {}
}];
```

`LoadConfiguration` resets every key to its schema default and then applies
the supplied rules.  `UpdateConfiguration` changes only the supplied keys.
Both validate the whole update before modifying global solver state.

## 2. Load an exact system

From a Wolfram expression:

```mathematica
x = Global`x;
eps = Global`eps;

sys = DiffExp2`API`LoadSystem[<|
  "Matrix" -> {
    {0, 1/x},
    {eps/(1 - x), 0}
  },
  "Variable" -> x
|>];
```

Or from an exact full-matrix export:

```mathematica
sys = DiffExp2`API`LoadSystem[<|
  "FullMatrixFile" -> "/absolute/path/to/dx_full.m",
  "Variable" -> x
|>];
```

The file must evaluate to a matrix.  Truncated `dx_0.m`, `dx_1.m`, … slice
directories are not accepted by this API because exact indicial and recurrence
decisions require the full epsilon-rational matrix.

`LoadSystem` returns:

```text
<|
  "Matrix" -> A,
  "Variable" -> x,
  "SingularFactors" -> {...}
|>
```

and clears solver caches belonging to the previous system.

## 3. Supply boundary values

For a finite boundary at `x0`, use a rectangular list:

```text
{
  {f1[eps^0], f1[eps^1], ..., f1[eps^K]},
  {f2[eps^0], f2[eps^1], ..., f2[eps^K]},
  ...
}
```

For example:

```mathematica
boundary = {
  {1, 0, 0, 0, 0},
  {2, -Log[2], Log[2]^2/2, -Log[2]^3/6, Log[2]^4/24}
};
```

The plain-array interface starts at epsilon order zero.  Internally, DiffExp 2
uses `EpsSeries` records with explicit lower and complete upper window bounds;
it never pads an unknown Laurent coefficient with zero.  Expert callers can
pass a validated `LocalSolution` to `TransportLine`, but the release examples
use finite arrays at regular anchors.

## 4. Plan and transport a line

The short form is:

```mathematica
result = DiffExp2`API`TransportEndpoint[
  sys, boundary, x0, x1,
  "ExtraSingularFactors" -> extraFactors
];
```

`ExtraSingularFactors` adds factors introduced by an observable or IBP
coefficient even when they are absent from the differential matrix.  This is
important: a rational coefficient with a pole inside a chart must be included
in segmentation before it is multiplied into a local solution.

To inspect or reuse the segmentation:

```mathematica
transportSystem = Join[sys, <|
  "ExtraSingularFactors" -> extraFactors
|>];

singularities = DiffExp2`Transport`FindSingularities[transportSystem];
plan = DiffExp2`Transport`SegmentLine[transportSystem, {x0, x1}];
DiffExp2`Transport`ValidatePlan[plan];

result = DiffExp2`Transport`TransportLine[
  transportSystem, boundary, plan
];
```

The planner uses exact roots to decide singular points and true complex-plane
distances for convergence radii.  The real projections `Re[z]` and
`Re[z] +/- Im[z]` of complex roots are retained as regular waypoints when
needed for stable classic matching; they are not fake singular charts.

## Transport result

The implemented result record contains:

| Key | Meaning |
| --- | --- |
| `"Value"` | `EpsSeries` of the vector at a regular endpoint; `None` at a singular endpoint |
| `"Final"` | final `LocalSolution`, including exact sectors |
| `"Charts"` | ordered list of `<|"Chart", "LocalSolution"|>` records |
| `"SegmentCount"` | number of retained chart solutions |
| `"EndpointIsSingular"` | whether the endpoint is a true singular point |
| `"ErrorEstimate"` | accumulated full-vs-reduced-order estimate by epsilon order |

Read a coefficient without discarding window metadata:

```mathematica
eps0Vector = DiffExp2`EpsSeries`ESCoefficient[result["Value"], 0];
window = DiffExp2`EpsSeries`ESWindow[result["Value"]];
```

Convert to a plain expression only at an API or presentation boundary:

```mathematica
plain = DiffExp2`EpsSeries`ESToExpression[result["Value"], eps];
```

## 5. Inspect segments and evaluate a chart

Each chart records its physical center, scale, convergence radius, incoming
match point, singular flag, and prescriptions.  The local coordinate is

```text
t = (x - chartCenter) / chartScale.
```

Evaluate a retained local solution inside its radius:

```mathematica
entry = First[result["Charts"]];
chart = entry["Chart"];
local = entry["LocalSolution"];
t = (probeX - chart["Center"])/chart["Scale"];

evaluated = DiffExp2`SectorSeries`EvaluateLocalSolution[local, t];
probeVector = DiffExp2`EpsSeries`ESCoefficient[
  evaluated["Value"], 0
];
```

`EvaluateLocalSolution` refuses points outside the chart radius.  It also
refuses a negative point on a multivalued chart unless an unambiguous
prescription supplies the imaginary side.

There is not yet a public `PlotSolution` convenience function.  The release
example builds a plot from these records without private symbols:
[SingularEndpointAndSegments.wl](../Examples/Direct/SingularEndpointAndSegments.wl).

## 6. Exact local behavior and singular endpoints

At a singular chart, use:

```mathematica
decomposition = DiffExp2`SectorSeries`SectorDecomposition[result["Final"]];
tags = KeyTake[#, {"a", "b", "p"}] & /@ decomposition["Sectors"];
```

A sector represents

```text
t^(a + b eps) (eps Log[t])^p / p!
```

times a regular vector Taylor series.  The `a`, `b`, and `p` values are exact
indicial/log-chain data.

For an epsilon-independent scalar observable at a singular endpoint:

```mathematica
limit = DiffExp2`API`EndpointLimitValues[result, coefficientVector];
```

The components are combined before the dimensional-regularization drop rule
and divergence gate.  `coefficientVector` must already be evaluated at the
endpoint and must be epsilon-free.  More general rational endpoint
coefficients should first be multiplied into component `LocalSolution`
objects; the current public API does not yet wrap that expert operation.

## 7. Integrate over a line

```mathematica
integral = DiffExp2`API`LineIntegral[
  sys,
  boundary,
  anchor,
  {lower, upper},
  coefficientVector,
  "ExtraSingularFactors" -> extraFactors
];
```

This computes

```text
Integral[ coefficientVector(x,eps) . f(x,eps), {x, lower, upper} ].
```

The coefficient vector must have one exact rational entry per master.  DiffExp
2 transports from the anchor in both directions when necessary, assigns each
interval tile to a certified half-radius chart, multiplies all active master
components, combines them, and only then applies endpoint regularization.  The
order is deliberate: divergent component integrals may cancel in the requested
observable.

If transports were already computed, `"PrecomputedCharts" -> charts` reuses
their chart records.  This option is implemented for the Feynman-trick ladder;
ordinary callers should normally let `LineIntegral` plan the interval.

## Failure semantics

DiffExp 2 errors are structured `Failure` objects thrown on the tag
`"DiffExp2Error"`.  Tests and advanced applications can catch them with:

```mathematica
Catch[expression, "DiffExp2Error"]
```

Unsupported indicial structure, missing branch prescriptions, incomplete
epsilon windows, failed recurrence backends, and invalid matching residuals
are loud failures.  None authorizes a fallback to a different solver.
