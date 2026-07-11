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
Get["/absolute/path/to/diffexp2/DiffExp2.m"];

DiffExp2`LoadConfiguration[
  "WorkingPrecision" -> 200,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 4,
  "DivisionOrder" -> 3,
  "RadiusOfConvergence" -> 1,
  "DeltaPrescriptions" -> {}
];
```

`LoadConfiguration` resets every key to its schema default and then applies
the supplied rules.  `UpdateConfiguration` changes only the supplied keys.
Both validate the whole update before modifying global solver state.

## 2. Load an exact system

From a Wolfram expression:

```mathematica
x = Global`x;
eps = Global`eps;

sys = DiffExp2`LoadSystem[<|
  "Matrix" -> {
    {0, 1/x},
    {eps/(1 - x), 0}
  },
  "Variable" -> x
|>];
```

Or from an exact full-matrix export:

```mathematica
sys = DiffExp2`LoadSystem[
  "/absolute/path/to/dx_full.m", "Variable" -> x
];
```

The file must evaluate to a matrix.  Truncated `dx_0.m`, `dx_1.m`, … slice
directories are not accepted by this API because exact indicial and recurrence
decisions require the full epsilon-rational matrix.

`LoadSystem` returns:

```text
<|
  "Schema" -> "DiffExp2.System/v1",
  "Matrix" -> A,
  "Variable" -> x,
  "SingularFactors" -> {...},       (* epsilon-zero planner alphabet *)
  "SingularFactorsExact" -> {...},  (* full epsilon-dependent factors *)
  "Dimension" -> Length[A],
  "Source" -> ...
|>
```

and clears solver caches belonging to the previous system.  `Matrix` is the
unchanged exact epsilon-rational matrix.  For segmentation only, each exact
singular factor is projected to its first nonzero epsilon coefficient after
removing any overall epsilon valuation.  Thus both `x + eps` and
`eps (x + eps)` place the limiting singularity at `x = 0`, while
`1 + eps x` contributes no finite planner root.  The unprojected factors stay
available in `SingularFactorsExact` for local and indicial calculations.
If a projected moving matrix pole itself becomes a required chart center, the
current finite-width local ansatz rejects the degenerating cleared denominator
with `E3`; the projection prevents incorrect line geometry but does not invent
unsupported `x/eps` asymptotics.

## 3. Supply boundary values

For a finite boundary at `x0`, use a rectangular list:

```text
{
  {f1[eps^0], f1[eps^1], ..., f1[eps^K]},
  {f2[eps^0], f2[eps^1], ..., f2[eps^K]},
  ...
}
```

For example, closed expressions can be expanded without manually transposing
coefficient tables:

```mathematica
boundary = DiffExp2`PrepareBoundary[
  {1, 2^(1-eps)}, "EpsilonOrder" -> 4
];
```

The equivalent explicit array is:

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
result = DiffExp2`TransportEndpoint[
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
plan = DiffExp2`PlanLine[
  sys, {x0, x1}, "ExtraSingularFactors" -> extraFactors
];
result = DiffExp2`TransportLine[
  sys, boundary, plan
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
| `"Schema"` | `"DiffExp2.TransportResult/v1"` |
| `"Plan"` | the validated line plan used for this transport |
| `"Value"` | `EpsSeries` of the vector at a regular endpoint; `None` at a singular endpoint |
| `"Final"` | final `LocalSolution`, including exact sectors |
| `"Charts"` | ordered list of `<|"Chart", "LocalSolution"|>` records |
| `"SegmentCount"` | number of retained chart solutions |
| `"EndpointIsSingular"` | whether the endpoint is a true singular point |
| `"ErrorEstimate"` | accumulated full-vs-reduced-order estimate by epsilon order |

Read a coefficient without discarding window metadata:

```mathematica
eps0Vector = DiffExp2`EpsilonCoefficient[result, 0];
window = DiffExp2`EpsilonWindow[result];
```

Convert to a plain expression only at an API or presentation boundary:

```mathematica
plain = DiffExp2`EpsilonExpression[result["Value"], eps];
```

## 5. Inspect segments and evaluate a chart

Each chart records its physical center, scale, convergence radius, incoming
match point, singular flag, and prescriptions.  The local coordinate is

```text
t = (x - chartCenter) / chartScale.
```

Evaluate a retained local solution inside its radius:

```mathematica
segment = First[DiffExp2`LineSegments[result]];
local = segment["LocalSolution"];
t = (probeX - segment["Center"])/segment["Scale"];

evaluated = DiffExp2`EvaluateLocal[local, t];
probeVector = DiffExp2`EpsilonCoefficient[evaluated, 0];
```

`EvaluateLocal` refuses points outside the chart radius. It also
refuses a negative point on a multivalued chart unless an unambiguous
prescription supplies the imaginary side.

`DiffExp2`PiecewiseSolution[result]["Function"]` is suitable for sampling in
`Plot` or for building custom data. There is no dedicated `PlotSolution`
styling wrapper. The release example builds a singular-endpoint plot from the
same public records:
[SingularEndpointAndSegments.wl](../Examples/Direct/SingularEndpointAndSegments.wl).

## 6. Exact local behavior and singular endpoints

At a singular chart, use:

```mathematica
sectors = DiffExp2`ExactSectors[result];
tags = KeyTake[#, {"a", "b", "p", "Exponent", "LogPower"}] & /@ sectors;
```

A sector represents

```text
t^(a + b eps) (eps Log[t])^p / p!
```

times a regular vector Taylor series.  The `a`, `b`, and `p` values are exact
indicial/log-chain data.

For an epsilon-independent scalar observable at a singular endpoint:

```mathematica
limit = DiffExp2`EndpointLimit[result, coefficientVector];
```

The components are combined before the dimensional-regularization drop rule
and divergence gate.  `coefficientVector` must already be evaluated at the
endpoint and must be epsilon-free.  More general rational endpoint
coefficients should first be multiplied into component `LocalSolution`
objects; the current public API does not yet wrap that expert operation.

## 7. Integrate over a line

```mathematica
integral = DiffExp2`IntegrateLine[
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
their chart records. This is primarily an expert/Feynman-trick optimization;
ordinary callers should normally let `IntegrateLine` plan the interval.

## Failure semantics

DiffExp 2 errors are structured `Failure` objects thrown on the tag
`"DiffExp2Error"`.  Tests and advanced applications can catch them with:

```mathematica
Catch[expression, "DiffExp2Error"]
```

Unsupported indicial structure, missing branch prescriptions, incomplete
epsilon windows, failed recurrence backends, and invalid matching residuals
are loud failures.  None authorizes a fallback to a different solver.
