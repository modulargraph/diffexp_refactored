# Quick Start

This page gives one direct transport and one Feynman-trick run.  The direct
workflow is the stable core API.  The Feynman-trick workflow currently uses the
repository driver rather than a polished one-call Wolfram Language function.

## Direct transport

Run the complete example:

```sh
wolframscript -file Examples/Direct/MinimalTransport.wl
```

Its essential steps are:

```mathematica
Get["/absolute/path/to/DiffExp2/DiffExp2/DiffExp2.m"];

DiffExp2`Config`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 100,
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 2,
  "DivisionOrder" -> 3
}];

x = Global`x;
eps = Global`eps;

sys = DiffExp2`API`LoadSystem[
  <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x|>
];

transport = DiffExp2`API`TransportEndpoint[
  sys, {{1/2, 0, 0}}, 0, 1
];

value = transport["Value"];
DiffExp2`EpsSeries`ESToExpression[value, eps]
(* {1/4} *)
```

The differential equation is `f'(x)=f(x)/(x-2)`, with `f(0)=1/2`.  The
boundary array has one row per master and one column per epsilon coefficient,
starting at epsilon order zero.

## Integrate a transported solution

`LineIntegral` transports both sides of the anchor as needed, tiles the
interval with certified charts, combines master components before endpoint
regularization, and returns an `EpsSeries`:

```mathematica
sys2 = DiffExp2`API`LoadSystem[
  <|"Matrix" -> {{eps/x}}, "Variable" -> x|>
];

bvals = Transpose@Table[
  {SeriesCoefficient[(11/23)^eps, {eps, 0, k}]},
  {k, 0, 3}
];

integral = DiffExp2`API`LineIntegral[
  sys2, bvals, 11/23, {0, 1}, {1}
];

DiffExp2`EpsSeries`ESToExpression[integral, eps]
(* 1 - eps + eps^2 - eps^3 + ... *)
```

Here the exact solution is `x^eps`, and dimensional regularization gives
`Integral[0,1] x^eps dx = 1/(1+eps)`.

## Inspect exact endpoint behavior

Transport to the singular endpoint and read its exact sector tags:

```mathematica
endpoint = DiffExp2`API`TransportEndpoint[
  sys2, bvals, 11/23, 0
];

decomposition = DiffExp2`SectorSeries`SectorDecomposition[
  endpoint["Final"]
];

KeyTake[#, {"a", "b", "p"}] & /@ decomposition["Sectors"]
```

The tags describe `x^(a+b eps) (eps Log[x])^p/p!` exactly.  They are solver
data, not exponents fitted from numerical epsilon coefficients.  See the
[singular endpoint example](../Examples/Direct/SingularEndpointAndSegments.wl)
for segment tables and a plot.

## Feynman-trick bubble

After installing FIRE and building the C++ backend:

```sh
sh Examples/FeynmanTrick/Bubble.sh
```

The runner:

1. constructs each Feynman-trick level;
2. asks FIRE for master bases, reductions, and differential matrices;
3. obtains the deepest analytic boundary condition;
4. transports to the lower and upper endpoints with DiffExp 2;
5. performs exact sector-aware limits or line integrals;
6. repeats until level zero and prints `STEPWISE` and `FINAL` JSON records.

The first run can be FIRE-dominated.  Later runs reuse `FT_PREP_CACHE_DIR` and
resume transport from `FT_LADDER_CHECKPOINT_DIR`.  Larger examples and the
meaning of every numerical setting are in [Feynman Trick](FeynmanTrick.md).

## Where to go next

- [Direct Solver](DirectSolver.md) explains systems, boundaries, plans,
  transport results, and integration.
- [Analytic Continuation](AnalyticContinuation.md) explains singularities,
  prescriptions, and exact local sectors.
- [API Reference](API.md) lists the implemented public symbols and current
  facade gaps.
- [Migration](Migration.md) maps the original DiffExp workflow onto DiffExp 2.
