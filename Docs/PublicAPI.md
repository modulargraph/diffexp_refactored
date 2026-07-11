# DiffExp 2 public API prototype

This document fixes the first release-facing slice. The recurrence,
matching, analytic-continuation, and regularized-integration algorithms stay
in the existing `DiffExp2/*` implementation modules. Users load the root
`DiffExp2.m` file and work only with the small `DiffExp2`` context.

The compiled C++ recurrence backend is the default. Selecting
`"RecurrenceBackend" -> "Wolfram"` remains possible for diagnostics, but
there is no silent fallback when the selected backend fails.

## First workflow

```wl
Get["/path/to/DiffExp2.m"];

DiffExp2`LoadConfiguration[
  "WorkingPrecision" -> 200,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 4
];

sys = DiffExp2`LoadSystem[
  <|"Matrix" -> A, "Variable" -> x|>
];

plan = DiffExp2`PlanLine[sys, {11/23, 1/3}];
result = DiffExp2`TransportLine[sys, boundary, plan];

DiffExp2`EpsilonCoefficient[result, 0]
DiffExp2`LineSegments[result]

endpoint = DiffExp2`TransportEndpoint[sys, boundary, 11/23, 0];
DiffExp2`LocalBehavior[endpoint]
DiffExp2`ExactSectors[endpoint]
DiffExp2`EndpointLimit[endpoint]
```

`LoadSystem` also accepts one exact `d<variable>_full.m` file, or a directory
containing exactly one such file. Multiple kinematic variables and legacy
epsilon-slice directories are intentionally not guessed by this first
slice; they require a later explicit line-parameterization API.

## Stable umbrella exports

Configuration:

- `LoadConfiguration`, `UpdateConfiguration`, `CurrentConfiguration`
- `BackendInformation`

Systems and lines:

- `LoadSystem`
- `PlanLine`, `TransportEndpoint`, `TransportLine`
- `LineSegments`, `LineSegment`, `EvaluateLine`, `PiecewiseSolution`
- `IntegrateLine`, `EndpointLimit`

Local and epsilon data:

- `EvaluateLocal`, `LocalBehavior`, `ExactSectors`
- `EpsilonWindow`, `EpsilonCoefficient`, `EpsilonCoefficientList`

Implementation symbols remain available in named subcontexts for package
development, but they are not part of the release compatibility promise.

## Named records

The umbrella adds a `"Schema"` key while retaining the proven internal data:

| Record | Schema | Important keys |
|---|---|---|
| system | `DiffExp2.System/v1` | `Matrix`, `Variable`, `SingularFactors`, `Dimension`, `Source` |
| plan | `DiffExp2.LinePlan/v1` | `From`, `To`, `Charts`, `Direction`, `ExtraSingularFactors` |
| transport | `DiffExp2.TransportResult/v1` | `Plan`, `Charts`, `Final`, `Value`, `EndpointIsSingular`, `ErrorEstimate` |
| segment | `DiffExp2.LineSegment/v1` | `Index`, `Domain`, `Center`, `Scale`, `Radius`, `LocalSolution` |
| piecewise | `DiffExp2.PiecewiseSolution/v1` | `Domain`, `Segments`, `Function` |

`PiecewiseSolution[result]["Function"]` evaluates through the certified
half-radius segment coverage. A singular endpoint is not sampled directly;
use `EndpointLimit` or inspect `LocalBehavior`.

`ExactSectors` exposes each exact local term as

```wl
<|
  "a" -> a,
  "b" -> b,
  "p" -> p,
  "Exponent" -> a + b eps,
  "LogPower" -> p,
  "Coefficients" -> tensor
|>
```

Thus the public API never has to fit a collapsed logarithmic series to guess
an `x^(a+b eps)` exponent.

Current limitation: matrices whose basis itself contains algebraic roots of
the line variable (for example square-root letters represented directly in
the differential-equation matrix) are not yet supported. `LoadSystem`
rejects fractional powers of variable-dependent expressions loudly. Exact
algebraic singularity locations of an otherwise rational matrix are
supported by the line planner; these are a different case.

## Feynman-trick facade

Loading `FeynmanTrick.m` also loads the DiffExp 2 umbrella and exports:

```wl
plan = FeynmanTrick`PipelinePlan["banana4_unequal"];
run = FeynmanTrick`RunIntegrationPipeline["banana4_unequal"];
resumed = FeynmanTrick`ResumeIntegrationPipeline[
  "banana4_unequal", checkpointFile
];
```

The minimal facade executes the already-tested
`Scripts/run_ft_stepwise2.m` ladder in a clean Wolfram subprocess. It uses
an argv list, not a shell string, and makes all environment settings visible
in `PipelinePlan`. This preserves the prepared FIRE cache and atomic
lower/upper-arm ladder checkpoints while the numerical loop is progressively
moved into package modules.

Pipeline schemas:

- `FeynmanTrick.PipelinePlan/v1`: command, working directory, exact
  environment, settings, cache/checkpoint directories, resume source.
- `FeynmanTrick.PipelineResult/v1`: status, exit code, parsed `FINAL` and
  `STEPWISE` records, stdout/stderr, and the originating plan.
- `FeynmanTrick.PipelineProcess/v1`: asynchronous process object and plan.

The facade defaults to C++, value-vector transport, endpoint-arm batching,
and up to ten native workers (bounded by `$ProcessorCount`). Custom topology
objects are the next extraction slice; the prototype accepts the validated
example-registry names used by the runner.
