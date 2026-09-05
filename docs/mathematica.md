# Mathematica interface

Load `DiffExp.m` from the repository, or the installed
`share/diffexp/Mathematica/DiffExp.wl`. The interface uses `RunProcess` with an
argument list and JSON over stdin/stdout. There is no LibraryLink, WSTP, compiled
Mathematica extension or persistent kernel session in the C++ backend.

The wrapper loads configuration and matrix files and converts legacy boundary
expressions to epsilon coefficients. C++ performs differential-equation
pullback, series recurrence, Frobenius matching, analytic continuation and
Feynman-trick recursion. The C++ package itself does not require Mathematica.

```wolfram
Get["/path/to/DiffExp2.1/DiffExp.m"];
$DiffExpExecutable = "/path/to/diffexp"; (* optional if built beside the package *)
LoadConfiguration[{
  MatrixDirectory -> "/path/to/matrices",
  EpsilonOrder -> 4, ExpansionOrder -> 60,
  WorkingPrecision -> 150, AccuracyGoal -> 25
}];
boundary = PrepareBoundaryConditions[initialValues, {t -> -1}];
result = TransportTo[boundary, {t -> 20}];
DiffExpLastTimings[]
```

## Original API

| Function | Behavior |
| --- | --- |
| `LoadConfiguration[rulesOrAssociation]` | Reset defaults and load matrices |
| `UpdateConfiguration[rulesOrAssociation]` | Update the current configuration |
| `CurrentConfiguration[]` | Current options as an Association |
| `PrepareBoundaryConditions[values, pointOrLine]` | Prepare epsilon coefficients and partial asymptotic input |
| `TransportTo[boundary, pointOrLine, to:1, save:False]` | Native transport; `{point, values, estimatedErrors}` |
| `ToPiecewise[saved, pade:False]` | Matrix of numerical functions callable as `functions[[i,k]][x]` |
| `IntegrateSystem[boundary,line]` | Numerical functions using the native transport interface |
| `DiffExpLastTimings[]` | Total wrapper time, individual native calls and chart count |

Ordinary files `dVARIABLE_ORDER.m` and canonical logarithmic `d_1.m` matrices
are supported, including `SparseArray`. Missing epsilon matrices are zero.
The wrapper accepts point/line Associations or rule lists, explicit nonlinear
paths, chained transport results, lists of epsilon coefficients, closed boundary
expressions in `eps`, unknown components `"?"`, and power/log asymptotic input
with `SeriesData` cutoffs. The original exported epsilon symbol and `eps` alias
are available. Causal coordinate deformations follow `DeltaPrescriptions`.

The native command is generic: it receives the actual matrix entries and path,
not an example identifier. Algebraic transport continues square-root sheets
along the path. The current algebraic implementation requires independent
polynomial radicands; nested algebraic towers are rejected explicitly.

Saved output has the original outer shape `{{point,values,errors}, metadata}`.
For a physical plotting point, `ToPiecewise` uses native transport and caches
the result. For a line varying one coordinate affinely, it starts from the
nearest previously computed physical sample, including the saved endpoint and
initial boundary. Each new point is evaluated on its physical path; local
polynomials on a complex contour are not extrapolated onto the real plotting
axis. A failed continuation from a cached sample falls back once to the original
boundary. More general curves replay from their original boundary. Saved
metadata preserves the original configuration and boundary so later
configuration changes do not alter the saved result. Numerical
plot coordinates are interpreted at their nominal values; excessive plotting
precision is capped at the configured integration precision before exact path
construction.

`WorkingPrecision`, `ExpansionOrder`, `EpsilonOrder`, `DivisionOrder` and
`AccuracyGoal` control the native calculation. Error estimates include retained
arithmetic and local truncation estimates; they are not full tail certificates.
Legacy strategy, Padé, Möbius and chopping options are retained for call
compatibility, while the native engine chooses its continuation method.
`ToPiecewise`'s Padé argument does not invoke Wolfram Padé integration.
Symbolic general solutions with unfixed constants and symbolic `SeriesData`
output from `IntegrateSystem` are not currently provided. Closed dimension-form
`dVARIABLE_d.m` loading and the fifth sampling-list argument are not yet part of
the compatibility surface. These limits are separate from running the original
notebook example workflows.

## Direct native operations

| Function | Result |
| --- | --- |
| `DiffExpBackendInfo[]` | Backend information |
| `DiffExpSeries[request]` | Exact regular-series coefficient strings |
| `DiffExpFamilyTemplate[name]` | Complete editable Feynman-family Association |
| `DiffExpFeynmanTrick[configuration,arguments:{}]` | Generic recursive integral calculation |
| `DiffExpRun[arguments,input:""]` | Process result, or `Failure` |

See [family configurations](feynman-families.md) for arbitrary FT families.
These direct operations accept `"Executable" -> "/path/to/diffexp"`.
Otherwise executable selection uses `$DiffExpExecutable`, the
`DIFFEXP_EXECUTABLE` environment variable, a nearby build/installation, then
`diffexp` on PATH. Compatibility transport uses the same global selection.

Process failures preserve stderr/stdout inside `Failure`. Response numbers are
parsed with a decimal grammar and never evaluated as Wolfram code. Precision
annotations on input boundary numbers are preserved. Exact-series coefficients
and FT Arb intervals remain strings to avoid machine-precision JSON rounding.

The optional smoke and compatibility tests run with:

```sh
cmake -S . -B build -DBUILD_TESTING=ON \
  -DDIFFEXP_WOLFRAMSCRIPT=/path/to/wolframscript \
  -DDIFFEXP_WOLFRAM_KERNEL=/path/to/WolframKernel
ctest --test-dir build -R mathematica --output-on-failure
```

Fetch original ancillary data as described in [examples.md](examples.md), then
run the notebook computational workflows sequentially:

```sh
python3 scripts/check_original_mathematica.py --kernel /path/to/WolframKernel
# The original 128-digit planar demonstration is a separate longer check:
python3 scripts/check_original_mathematica.py --kernel /path/to/WolframKernel --case zzz-high --timeout 3600
```

Use `--case banana --case banana-routes`, for example, to select a workflow.
The runner records per-case logs and elapsed time under
`build-reference/mathematica-workflows`, with a finite per-case timeout. It
checks original polylogarithms, all four planar families, supplied-boundary
Henn transport, partial banana initialization, chained equal/unequal banana
routes and saved-output numerical evaluation. The high-precision demonstration
is opt-in. No full Henn FT reconstruction is started.

The [native JSON protocol](transport-protocol.md) is available for other language
wrappers. Mathematica tests are optional and are not required to build or test
the standalone native core.
