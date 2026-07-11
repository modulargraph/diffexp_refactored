# API Reference

The current loader is:

```mathematica
Get["/absolute/path/to/DiffExp2/DiffExp2/DiffExp2.m"];
```

Symbols are currently organized in subcontexts.  This page distinguishes the
small user-facing surface from lower-level expert functions.  A top-level
top-level DiffExp 2 convenience facade is a release usability gap and is not invented
here.

## Core user API

### ``DiffExp2`Config`LoadConfiguration[rules]``

Reset every configuration key to its default, validate `rules`, install the
new tolerance state, and return the resolved configuration association.

```mathematica
DiffExp2`Config`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 300,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 4
}]
```

### ``DiffExp2`Config`UpdateConfiguration[rules]``

Validate and atomically merge only the supplied keys.

### ``DiffExp2`Config`CurrentConfiguration[]``

Return every known key and its resolved value.  Automatic chop and
rationalization thresholds are resolved in the returned association.

### ``DiffExp2`API`LoadSystem[spec]``

Load an exact epsilon-rational system.

```mathematica
LoadSystem[<|"Matrix" -> matrix, "Variable" -> x|>]
LoadSystem[<|"FullMatrixFile" -> path, "Variable" -> x|>]
```

Returns `<|"Matrix", "Variable", "SingularFactors"|>`.

### ``DiffExp2`API`TransportEndpoint[sys,boundary,from,to,opts]``

Plan and transport from a regular anchor to one endpoint.

Option:

- `"ExtraSingularFactors" -> {}` adds observable/IBP denominator factors to
  segmentation.

At a regular endpoint, `result["Value"]` is an `EpsSeries`.  At a singular
endpoint it is `None`, and `result["Final"]` is the exact `LocalSolution`.

### ``DiffExp2`API`LineIntegral[sys,boundary,from,{lo,hi},cvec,opts]``

Return the `EpsSeries` for
`Integrate[cvec(x,eps).f(x,eps),{x,lo,hi}]`.

Options:

- `"ExtraSingularFactors" -> {}`;
- `"PrecomputedCharts" -> None`, used by the Feynman-trick ladder to reuse
  already completed endpoint arms.

`cvec` must have one rational function per system component.  Components are
combined before regularization and divergence checks.

### ``DiffExp2`API`EndpointLimitValues[transportResult,cvec]``

Return the dimensional-regularization limit of an epsilon-free scalar
combination at a singular endpoint.  `cvec` entries must already be evaluated
at that endpoint.

## Epsilon-series accessors

An `EpsSeries` has the structural form:

```text
<|
  "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "Coeffs" -> {c[kmin], ..., c[kmax]}
|>
```

The most useful public functions are:

| Function | Purpose |
| --- | --- |
| `ESNew[kmin,coeffs]` | construct a complete Laurent window |
| `ESZero[kmax]` | exact zero known complete through `kmax` |
| `ESWindow[s]` | read `{Min,CompleteMax}` metadata |
| `ESCoefficient[s,k]` | read one coefficient; loud above `CompleteMax` |
| `ESCoefficientList[s,k1,k2]` | read a certified range without dropping material lower content |
| `ESToExpression[s,eps]` | convert to a plain expression and lose window metadata |
| `ESFromExpression[expr,eps,kmax]` | exact Laurent conversion at an API boundary |
| `ESAdd`, `ESTimes`, `ESDivide` | honest-window Laurent arithmetic |
| `ESTrim`, `ESTruncate` | remove negligible leading content or lower the complete top |

These symbols live in the DiffExp 2 EpsSeries module.

## Local-solution and sector API

The user-relevant functions in the DiffExp 2 SectorSeries module are:

| Function | Purpose |
| --- | --- |
| `EvaluateLocalSolution[ls,t,opts]` | evaluate inside one chart |
| `SectorDecomposition[ls]` | expose exact `{a,b,p}` sectors and windows |
| `MultiplyRational[ls,c,t]` | multiply by a rational chart-coordinate factor |
| `CombineLocalSolutions[weights,lss]` | exact linear combination on identical charts |
| `DifferentiateLocalSolution[ls]` | chart-coordinate derivative |
| `ReexpandLocalSolution[ls,newCenter,targetOrder]` | re-expand inside a regular part of the chart |
| `ParseTaggedPower[expr,var,eps]` | parse one `c var^(a+b eps) Log[var]^p` monomial |

`EvaluateLocalSolution` options are:

- `"UsePade" -> Automatic`;
- `"TOrderReduction" -> 0`;
- `"ImSign" -> Automatic`.

The internal option `"ComputeTailEstimates"` is exposed by the implementation
but is not a release-stable user control.

## Planning and transport API

The expert functions in the DiffExp 2 Transport module are:

| Function | Purpose |
| --- | --- |
| `FindSingularities[sys]` | exact roots, real roots, factors, and projection waypoints |
| `ChartRadius[center,allSingularities]` | distance to the nearest other true singularity |
| `SegmentLine[sys,{from,to}]` | build a complete `SegmentPlan` |
| `ValidatePlan[plan]` | audit overlap, half-radius probes, and singular approach sides |
| `TransportLine[sys,boundary,plan]` | execute an already built plan |

`MatchWeights`, `ApplyCrossing`, and `SegmentErrorProbe` are exported from the
module for focused testing and expert diagnostics, but they are not ordinary
entry points.

## Integration primitives

The expert functions in the DiffExp 2 Integrate module are:

| Function | Purpose |
| --- | --- |
| `IntegrateLocalSolution[ls,{t1,t2}]` | integrate one chart, including endpoint/PV cases |
| `EndpointSectorLimit[ls]` | componentwise sector limit at `t=0` |
| `SectorMonomialIntegral[...]` | exact monomial integral used by the chart integrator |

Most applications should use `LineIntegral` or `EndpointLimitValues`, which
combine master components before the cancellation gate.

## C++ backend inspection

```mathematica
DiffExp2`CppBackend`BackendAvailableQ[]
DiffExp2`CppBackend`BackendInformation[]
DiffExp2`CppBackend`ResetBackend[]
```

The serialization functions are implementation interfaces and should not be
called by ordinary users.

## Common configuration keys

| Key | Implemented default in prototype | Meaning |
| --- | ---: | --- |
| `"RecurrenceBackend"` | `"Wolfram"` | `"Cpp"` is the intended release default; mismatch still to patch |
| `"WorkingPrecision"` | `500` | recurrence and evaluation precision |
| `"ChopPrecision"` | `Automatic` | resolved from working precision |
| `"LinearSolveChopPrecision"` | `Automatic` | matching/pivot threshold digits |
| `"ExpansionOrder"` | `50` | complete chart Taylor order |
| `"EpsilonOrder"` | `4` | requested complete epsilon top |
| `"DivisionOrder"` | `3` | coupled placement/matching divisor |
| `"RadiusOfConvergence"` | `1` | affine local coordinate normalization |
| `"DeltaPrescriptions"` | `{}` | physical sides of singular factors |
| `"UsePade"` | `False` | Pade evaluation of local Taylor slices |
| `"EstimateError"` | `"Fast"` | retained compatibility/error-estimate mode |
| `"Variables"` | `{}` | declared kinematic symbols other than line parameter |
| `"Verbosity"` | `1` | user-facing progress output |

Use `ConfigSchema[]` and `CurrentConfiguration[]` rather than depending on
this summary for every advanced key.

## Feynman-trick API status

The FeynmanTrick package exports configuration and lower-level topology/iteration
functions, but the sector-native DiffExp 2 ladder is currently implemented by
`Scripts/run_ft_stepwise2.m`.  The supported release workflow is therefore the
documented CLI.  Promoting that driver to a stable Wolfram Language function
is an explicit release-interface task, not an undocumented existing feature.
