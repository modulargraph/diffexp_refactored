# API Reference

This page defines the release-facing Wolfram Language interfaces. Load the
root files; implementation modules under `DiffExp2/` and `FeynmanTrick/` are
not compatibility entry points.

```mathematica
repo = "/absolute/path/to/diffexp2";
Get[FileNameJoin[{repo, "DiffExp2.m"}]];
```

The loader installs the stable root umbrella, validates a default
configuration, and selects the compiled C++ recurrence backend. There is no
silent backend fallback.

## Umbrella exports

Configuration and backend:

- `$DiffExp2Version`
- `LoadConfiguration`, `UpdateConfiguration`, `CurrentConfiguration`
- `BackendInformation`

Systems, planning, and transport:

- `LoadSystem`, `PrepareBoundary`
- `PlanLine`, `TransportEndpoint`, `TransportLine`
- `LineSegments`, `LineSegment`, `EvaluateLine`, `PiecewiseSolution`

Local solutions, endpoint operations, and epsilon data:

- `EvaluateLocal`, `LocalBehavior`, `ExactSectors`
- `EndpointLimit`, `IntegrateLine`
- `EpsilonWindow`, `EpsilonCoefficient`, `EpsilonCoefficientList`
- `EpsilonExpression`

## Configuration

```mathematica
DiffExp2`LoadConfiguration[]
DiffExp2`LoadConfiguration[rules...]
DiffExp2`LoadConfiguration[{rules...}]
DiffExp2`LoadConfiguration[association]

DiffExp2`UpdateConfiguration[rules...]
DiffExp2`UpdateConfiguration[{rules...}]
DiffExp2`UpdateConfiguration[association]

DiffExp2`CurrentConfiguration[]
```

`LoadConfiguration` resets every key before applying its arguments.
`UpdateConfiguration` validates and merges only the supplied keys. Both are
atomic and return the resolved, string-keyed configuration. When no backend
is supplied, the release umbrella chooses `"Cpp"`; selecting `"Wolfram"`
explicitly enables the reference recurrence.

Common defaults at the release boundary are:

| Key | Default | Meaning |
| --- | ---: | --- |
| `"RecurrenceBackend"` | `"Cpp"` | strict recurrence implementation |
| `"WorkingPrecision"` | `500` | recurrence and evaluation precision |
| `"ExpansionOrder"` | `50` | complete local Taylor order |
| `"EpsilonOrder"` | `4` | requested complete epsilon top |
| `"DivisionOrder"` | `3` | coupled chart placement and matching divisor |
| `"RadiusOfConvergence"` | `1` | affine chart-coordinate normalization |
| `"DeltaPrescriptions"` | `{}` | signed physical-side prescriptions |
| `"UsePade"` | `False` | default local-series evaluation mode |
| `"Variables"` | `{}` | additional declared kinematic symbols |
| `"Verbosity"` | `1` | user-facing progress level |

Use `CurrentConfiguration[]` to inspect every resolved key and automatic
tolerance.

```mathematica
DiffExp2`BackendInformation[]
```

`BackendInformation` reports compiled-library discovery, loading, version,
and runtime details through the stable umbrella.

## Systems

```mathematica
DiffExp2`LoadSystem[association]
DiffExp2`LoadSystem[file, "Variable" -> Automatic]
DiffExp2`LoadSystem[directory, "Variable" -> Automatic]
```

An in-memory system has the form:

```mathematica
<|
  "Matrix" -> matrix,
  "Variable" -> x
|>
```

An association may instead provide an exact matrix file:

```mathematica
<|
  "FullMatrixFile" -> "/path/to/dx_full.m",
  "Variable" -> x
|>
```

A string may name one exact `d<variable>_full.m` file or a directory
containing exactly one such file. The filename determines the variable unless
the `"Variable"` option overrides it. The matrix must be exact and rational
in the line variable and epsilon.

The result has schema `DiffExp2.System/v1` and includes:

| Key | Meaning |
| --- | --- |
| `"Matrix"` | exact differential matrix |
| `"Variable"` | line variable |
| `"SingularFactors"` | exact matrix-denominator factors |
| `"Dimension"` | number of masters |
| `"Source"` | in-memory marker or expanded matrix path |

## Regular-anchor boundaries

```mathematica
DiffExp2`PrepareBoundary[expressions, opts]
```

`expressions` contains one closed-form expression per master. The function
returns the masters-by-epsilon-coefficients array accepted by
`TransportEndpoint`, `TransportLine`, and `IntegrateLine`.

| Option | Default | Meaning |
| --- | ---: | --- |
| `"EpsilonSymbol"` | `Automatic` | expansion symbol; automatic uses the canonical global `eps` symbol |
| `"EpsilonOrder"` | `Automatic` | nonnegative expansion top; automatic uses the current configuration |

`PrepareBoundary` is intentionally limited to regular-anchor expressions
finite at epsilon zero. A pole would lose a nonzero Laurent lower window, so
the function rejects it; pass an explicit exact local solution when such a
window must be preserved.

## Planning and transport

```mathematica
DiffExp2`PlanLine[sys, {from, to}, opts]

DiffExp2`TransportEndpoint[sys, boundary, from, to, opts]
DiffExp2`TransportLine[sys, boundary, {from, to}, opts]
DiffExp2`TransportLine[sys, boundary, plan]
```

The forms that build a plan accept:

| Option | Default | Meaning |
| --- | ---: | --- |
| `"ExtraSingularFactors"` | `{}` | observable or IBP denominator factors that also constrain segmentation |

`boundary` is a masters-by-epsilon-coefficients array. `PlanLine` performs
exact singularity analysis, constructs the affine chart sequence, and
validates its coverage without solving. Passing the resulting plan to
`TransportLine` reuses it.

A plan has schema `DiffExp2.LinePlan/v1`. Its important keys are `"From"`,
`"To"`, `"Direction"`, `"Charts"`, `"EndpointIsSingular"`,
`"Singularities"`, and `"ExtraSingularFactors"`.

A completed transport has schema `DiffExp2.TransportResult/v1` and retains
`"Plan"`, `"Charts"`, `"Final"`, `"Value"`,
`"EndpointIsSingular"`, and `"ErrorEstimate"`. At a regular endpoint,
`"Value"` is an honest epsilon series. At a singular endpoint, `"Value"` is
`None` and `"Final"` is the exact local solution; use `EndpointLimit` or the
local-sector functions rather than sampling the endpoint.

## Segments and piecewise evaluation

```mathematica
DiffExp2`LineSegments[result]
DiffExp2`LineSegment[result, i]

DiffExp2`EvaluateLine[result, point, opts]
DiffExp2`PiecewiseSolution[result]
```

`LineSegments` returns records with schema `DiffExp2.LineSegment/v1`. Each
contains `"Index"`, `"Domain"`, `"Center"`, `"Scale"`, `"Radius"`,
`"Chart"`, and `"LocalSolution"`. `LineSegment` uses one-based indexing.

`EvaluateLine` chooses a certified covering segment and returns its local
evaluation record. Its options are the same as `EvaluateLocal`:

| Option | Default | Meaning |
| --- | ---: | --- |
| `"UsePade"` | `Automatic` | inherit the configured Pade choice or override it |
| `"TOrderReduction"` | `0` | omit this many highest Taylor columns |
| `"ImSign"` | `Automatic` | derive or explicitly select the branch side |
| `"ComputeTailEstimates"` | `True` | compute advisory Taylor-tail diagnostics |

`PiecewiseSolution` returns schema `DiffExp2.PiecewiseSolution/v1` with
`"Domain"`, `"Segments"`, and a callable `"Function"`. Its function uses
the certified half-radius coverage. It deliberately refuses to sample a
singular endpoint directly.

## Exact local behavior

```mathematica
DiffExp2`EvaluateLocal[localOrRecord, t, opts]
DiffExp2`LocalBehavior[localOrResult]
DiffExp2`ExactSectors[localOrResult]
```

`EvaluateLocal` accepts an exact local solution, a segment record, or a
transport record and uses the options listed above. `LocalBehavior` returns
the canonical sector decomposition. `ExactSectors` exposes each term as:

```mathematica
<|
  "a" -> a,
  "b" -> b,
  "p" -> p,
  "Exponent" -> a + b eps,
  "LogPower" -> p,
  "Coefficients" -> tensor
|>
```

These tags are exact recurrence data. The API never fits sampled epsilon
coefficients to guess an exponent of the form `x^(a + b eps)`.

## Endpoint limits and line integration

```mathematica
DiffExp2`EndpointLimit[singularTransport]
DiffExp2`EndpointLimit[singularTransport, weights]

DiffExp2`IntegrateLine[
  sys, boundary, from, {lo, hi}, coefficients, opts
]
```

`EndpointLimit[result]` returns the dimensionally regularized limit of every
component. The weighted form combines components before the cancellation
gate; its weights must be epsilon-free values already evaluated at the
endpoint.

`IntegrateLine` transports both directions from the anchor, combines the
rational master coefficients before regularization, and returns an honest
epsilon series. Its options are:

| Option | Default | Meaning |
| --- | ---: | --- |
| `"ExtraSingularFactors"` | `{}` | include coefficient-denominator factors in segmentation |
| `"PrecomputedCharts"` | `None` | reuse a complete chart list, primarily for the Feynman ladder |

## Epsilon windows

```mathematica
DiffExp2`EpsilonWindow[value]
DiffExp2`EpsilonCoefficient[value, k]
DiffExp2`EpsilonCoefficientList[value, k1, k2]
DiffExp2`EpsilonExpression[value]
DiffExp2`EpsilonExpression[value, eps]
```

An honest epsilon series records a Laurent window:

```text
<|"Min" -> kmin, "CompleteMax" -> kmax|>
```

`EpsilonWindow` accepts an epsilon series, an evaluation or transport record,
a local solution, or a list of them. `EpsilonCoefficient` and
`EpsilonCoefficientList` accept epsilon series, regular evaluation or
transport records whose `"Value"` is an epsilon series, and lists. Reads
above `"CompleteMax"`, or reads that would discard nonzero lower-window
content, fail loudly. A singular transport has no regular `"Value"`; inspect
its exact sectors or take its endpoint limit instead.

`EpsilonExpression` accepts one complete epsilon series, a regular evaluation
or transport record, or a list. It converts the certified Laurent window to
an ordinary expression in the requested symbol for presentation. The result
deliberately no longer carries honest-window metadata, so retain the original
record for further solver operations.

## Named record schemas

| Record | Schema |
| --- | --- |
| system | `DiffExp2.System/v1` |
| line plan | `DiffExp2.LinePlan/v1` |
| transport | `DiffExp2.TransportResult/v1` |
| line segment | `DiffExp2.LineSegment/v1` |
| piecewise solution | `DiffExp2.PiecewiseSolution/v1` |

Schema tags are intended for validation and serialization boundaries. The
associations retain the underlying exact mathematical data rather than
collapsing it into sampled values.

## Feynman-trick facade

Load its root file:

```mathematica
Get[FileNameJoin[{repo, "FeynmanTrick.m"}]];
```

The release-facing pipeline functions are:

```mathematica
FeynmanTrick`$FeynmanTrickVersion
FeynmanTrick`SupportedExamples[]
FeynmanTrick`PipelinePlan[example, opts]
FeynmanTrick`RunIntegrationPipeline[example, opts]
FeynmanTrick`RunIntegrationPipeline[plan]
FeynmanTrick`ResumeIntegrationPipeline[example, checkpoint, opts]
```

`example` must be one registry name, not a comma-separated list or a custom
topology object. The current registry is:

```text
bubble, sunrise, banana, banana_unequal, banana4, banana4_unequal,
kite, box, pentagon, box_bubble, box_triangle, double_box_planar
```

`PipelinePlan` validates and exposes the exact command, environment, cache
locations, and settings without running. `RunIntegrationPipeline` either
builds a plan from an example name and options or executes an existing
`FeynmanTrick.PipelinePlan/v1` record. `ResumeIntegrationPipeline` sets the
resume checkpoint and runs.

All three pipeline functions share these options:

| Option | Default | Meaning |
| --- | ---: | --- |
| `"WorkingPrecision"` | `500` | Wolfram numerical precision |
| `"ExpansionOrder"` | `50` | local Taylor order |
| `"EpsilonOrder"` | `0` | requested final epsilon top |
| `"BoundaryExtraOrder"` | `4` | extra level-boundary epsilon orders |
| `"LevelEpsilonHalos"` | `{0}` | per-level internal epsilon lookahead |
| `"DivisionOrder"` | `3` | coupled placement and `+/-1/k` matching divisor |
| `"RadiusOfConvergence"` | `1` | affine chart-coordinate normalization |
| `"RecurrenceBackend"` | `"Cpp"` | strict recurrence backend |
| `"CppThreads"` | `Automatic` | up to ten native workers, bounded by `$ProcessorCount` |
| `"ValueTransport"` | `True` | use regular-chart value transport |
| `"BatchEndpointArms"` | `True` | request paired native endpoint-arm prewarming when applicable |
| `"SingularMatchPrecondition"` | `False` | enable singular-match preconditioning |
| `"PreparedCacheDirectory"` | `Automatic` | FIRE preparation cache root |
| `"FIREPath"` | `Automatic` | FIRE6 installation; defaults to `Dependencies/fire/FIRE6` |
| `"CheckpointDirectory"` | `Automatic` | ladder checkpoint directory |
| `"ResumeFrom"` | `None` | existing checkpoint file |
| `"RebuildPreparation"` | `False` | ignore and rebuild prepared FIRE data |
| `"AllowStaleCheckpoint"` | `False` | permit explicitly stale ladder metadata |
| `"StopAfterBoundaryLevel"` | `None` | stop after a nonnegative level |
| `"FIRETimeoutSeconds"` | `1800` | watchdog for one FIRE invocation |
| `"Runner"` | `Automatic` | path to the tested runner script |
| `"WolframScript"` | `Automatic` | `wolframscript` executable |
| `"WorkingDirectory"` | `Automatic` | subprocess working directory |
| `"ExtraEnvironment"` | `<||>` | additional string-valued environment entries |
| `"Asynchronous"` | `False` | return a running process record instead of waiting |

The facade returns one of:

| Record | Schema | Important keys |
| --- | --- | --- |
| plan | `FeynmanTrick.PipelinePlan/v1` | `Example`, `Command`, `Environment`, `Settings`, cache and checkpoint paths |
| result | `FeynmanTrick.PipelineResult/v1` | `Status`, `ExitCode`, parsed `Final` and `Stepwise`, output streams, `Plan` |
| asynchronous process | `FeynmanTrick.PipelineProcess/v1` | `Status`, `Process`, `Plan` |

The current facade launches the runner through `/usr/bin/env`, so process
execution is supported on macOS and Linux. It starts a clean
`wolframscript` kernel; when called from another kernel, the two processes may
occupy two Wolfram license seats. Native C++ worker threads do not each start
a Wolfram kernel.

The command-line driver remains available for advanced multi-example runs:

```sh
FT_EXAMPLES=bubble,banana \
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

## Compatibility boundary

Named implementation contexts remain accessible for package development and
focused tests, but they are unstable implementation details. Release code
should load only the root files and call the umbrella symbols documented on
this page. In particular, do not depend on low-level configuration defaults,
private record shapes, or backend serialization functions.
