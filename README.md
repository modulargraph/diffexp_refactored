# DiffExp 2

DiffExp 2 transports local recurrence solutions of differential systems for
dimensionally regularized Feynman integrals. It keeps exact local behavior

```text
x^(a + b eps) Log[x]^p
```

through transport, matching, endpoint limits, and analytic regularization.
The C++20/FLINT recurrence backend is the release default; Wolfram Language
retains exact singularity analysis, branch control, matching, integration,
and validation. Selecting either backend is strict: failures never trigger a
silent fallback to the other implementation.

## Quick start

Build the compiled backend as described in
[Installation](Docs/Installation.md), then load the root package:

```mathematica
repo = "/absolute/path/to/diffexp2";
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

DiffExp2`LoadConfiguration[
  "WorkingPrecision" -> 100,
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 2
];

x = Global`x;
sys = DiffExp2`LoadSystem[
  <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x|>
];

boundary = DiffExp2`PrepareBoundary[{1/2}];

result = DiffExp2`TransportEndpoint[
  sys,
  boundary,
  0,
  1
];

DiffExp2`EpsilonCoefficient[result, 0]
(* {1/4} *)

DiffExp2`EpsilonExpression[result, Global`eps]
(* {1/4} *)
```

The root `DiffExp2.m` loader installs the stable umbrella API and
selects the compiled backend by default.

## Feynman trick

Load the root Feynman-trick package to plan or run one validated example:

```mathematica
Get[FileNameJoin[{repo, "FeynmanTrick.m"}]];

plan = FeynmanTrick`PipelinePlan["bubble"];
result = FeynmanTrick`RunIntegrationPipeline[plan];
```

The facade runs the tested ladder in a clean Wolfram subprocess, reuses
prepared FIRE data and atomic ladder checkpoints, and returns named plan,
process, or result records. The current process facade requires
`/usr/bin/env`, so it is supported on macOS and Linux. Direct DiffExp 2
transport is not subject to that platform restriction.

The underlying driver remains available for advanced and batch use:

```sh
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
FT_EXAMPLES=bubble \
FT_WORKING_PRECISION=300 \
FT_EXPANSION_ORDER=40 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

See [Feynman Trick](Docs/FeynmanTrick.md) for registry names, numerical
settings, cache controls, and checkpoint resumption.

## Capabilities

- exact epsilon-rational systems at ordinary and regular-singular points;
- finite-width recurrences, resonances, log sectors, fractional powers, and
  inhomogeneous sources;
- exact real-line segmentation with complex-root projections and explicit
  `+i0`/`-i0` prescriptions;
- inspectable line segments and piecewise evaluation;
- exact sector-aware endpoint limits and analytically regularized line
  integrals;
- FIRE preparation and recursive Feynman-trick integration;
- strict Wolfram/C++ recurrence parity tests and parallel native solves.

## Documentation

- [Installation](Docs/Installation.md)
- [Quick Start](Docs/QuickStart.md)
- [API Reference](Docs/API.md)
- [Direct Solver](Docs/DirectSolver.md)
- [Feynman Trick](Docs/FeynmanTrick.md)
- [Analytic Continuation](Docs/AnalyticContinuation.md)
- [Verified Results](Docs/Results.md)
- [Migration from DiffExp 1](Docs/Migration.md)
- [C++ Backend](Docs/CppBackend.md)
- [Release Manifest](Docs/ReleaseManifest.md)
- [Citation](Docs/Citation.md)
- [Changelog](CHANGELOG.md)

## Current limitations

- The differential matrix must be exact and rational in the line variable
  and epsilon. `LoadSystem` rejects non-integer powers with a
  variable-dependent base, such as `Sqrt[x]`. Exact algebraic constants,
  singular locations, and local exponents are a different case and are not
  rejected merely for being algebraic.
- Indicial eigenvalues must be affine in epsilon, `a + b eps`.
- There is no dedicated plotting command. Use `LineSegments` or
  `PiecewiseSolution` with the plotting style appropriate to the application.
- FIRE is required only to prepare new Feynman-trick families. Transporting
  an existing exact differential matrix does not require FIRE.

DiffExp 2 is the recurrence-based successor to
[DiffExp](https://gitlab.com/hiddingm/diffexp).

## License

Copyright 2024–2026 Martijn Hidding. DiffExp 2 is free software distributed
under the [GNU General Public License, version 3 or later](LICENSE).
