# Quick Start

This page gives a direct transport, a regularized line integral, and a
Feynman-trick run through the release-facing Wolfram Language APIs.

Set the repository path once and load the root package:

```mathematica
repo = "/absolute/path/to/diffexp2";
Get[FileNameJoin[{repo, "DiffExp2.m"}]];
```

The root loader selects the compiled recurrence backend and installs a
validated default configuration.

## Direct transport

```mathematica
DiffExp2`LoadConfiguration[
  "WorkingPrecision" -> 100,
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 3,
  "DivisionOrder" -> 3
];

x = Global`x;
eps = Global`eps;

sys = DiffExp2`LoadSystem[
  <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x|>
];

boundary = DiffExp2`PrepareBoundary[{1/2}];

transport = DiffExp2`TransportEndpoint[
  sys,
  boundary,
  0,
  1
];

DiffExp2`EpsilonCoefficient[transport, 0]
(* {1/4} *)

DiffExp2`EpsilonExpression[transport, eps]
(* {1/4} *)
```

The equation is `f'(x) = f(x)/(x - 2)`, with `f(0) = 1/2`. The boundary
helper expands one finite regular-anchor expression per master into the
masters-by-epsilon-coefficients array expected by transport.

Transport records retain every local chart. Inspect or evaluate them through
the umbrella API:

```mathematica
segments = DiffExp2`LineSegments[transport];
atThreeQuarters = DiffExp2`EvaluateLine[transport, 3/4];
piecewise = DiffExp2`PiecewiseSolution[transport];
```

## Integrate a transported solution

`IntegrateLine` transports both sides of the anchor as needed, tiles the
interval with certified charts, combines master components before endpoint
regularization, and returns an honest epsilon series:

```mathematica
sys2 = DiffExp2`LoadSystem[
  <|"Matrix" -> {{eps/x}}, "Variable" -> x|>
];

bvals = DiffExp2`PrepareBoundary[
  {(11/23)^eps},
  "EpsilonOrder" -> 3
];

integral = DiffExp2`IntegrateLine[
  sys2,
  bvals,
  11/23,
  {0, 1},
  {1}
];

DiffExp2`EpsilonCoefficientList[integral, 0, 3]
(* {1, -1, 1, -1} *)

DiffExp2`EpsilonExpression[integral, eps]
(* 1 - eps + eps^2 - eps^3 *)
```

The exact solution is `x^eps`, and dimensional regularization gives
`Integral[0,1] x^eps dx = 1/(1 + eps)`.

## Inspect exact endpoint behavior

Transport to the singular endpoint and read its exact sector tags:

```mathematica
endpoint = DiffExp2`TransportEndpoint[
  sys2,
  bvals,
  11/23,
  0
];

sectors = DiffExp2`ExactSectors[endpoint];
KeyTake[#, {"a", "b", "p"}] & /@ sectors
```

The tags describe `x^(a + b eps) Log[x]^p` exactly. They are solver data,
not exponents inferred from numerical epsilon samples. Use `LocalBehavior`
for the complete canonical decomposition or `EndpointLimit` for the
dimensionally regularized componentwise limit.

## Feynman-trick bubble

After installing FIRE and building the C++ backend, load the root facade:

```mathematica
Get[FileNameJoin[{repo, "FeynmanTrick.m"}]];

plan = FeynmanTrick`PipelinePlan[
  "bubble",
  "WorkingPrecision" -> 300,
  "ExpansionOrder" -> 40,
  "CppThreads" -> 4
];

run = FeynmanTrick`RunIntegrationPipeline[plan];

run["Status"]
run["Final"]
```

`PipelinePlan` is inspectable and does not execute the ladder. Pass it to
`RunIntegrationPipeline` to execute exactly those settings, or call
`RunIntegrationPipeline["bubble", opts]` directly. The facade launches a
clean `wolframscript` subprocess, so calling it from an existing kernel may
occupy a second Wolfram license seat while the run is active. Its current
`/usr/bin/env` launcher is supported on macOS and Linux.

The runner:

1. constructs each Feynman-trick level;
2. asks FIRE for master bases, reductions, and differential matrices;
3. obtains the deepest analytic boundary condition;
4. transports to the lower and upper endpoints with DiffExp 2;
5. performs exact sector-aware limits or line integrals;
6. repeats to level zero and returns parsed `STEPWISE` and `FINAL` records.

Prepared FIRE data and transport checkpoints are independent, so a solver
rerun need not repeat the reduction. The command-line driver remains
available for advanced or multi-example runs:

```sh
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
FT_EXAMPLES=bubble \
FT_WORKING_PRECISION=300 \
FT_EXPANSION_ORDER=40 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

See [Feynman Trick](FeynmanTrick.md) for registry names, cache controls,
checkpoint resumption, and every numerical setting.

## Where to go next

- [API Reference](API.md) defines the stable umbrella symbols and named
  records.
- [Direct Solver](DirectSolver.md) explains systems, boundaries, plans,
  transport, and integration.
- [Analytic Continuation](AnalyticContinuation.md) explains singularities,
  prescriptions, and exact local sectors.
- [Migration](Migration.md) maps the original DiffExp workflow onto DiffExp 2.
