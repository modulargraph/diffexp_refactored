# DiffExp 2

DiffExp 2 solves systems of differential equations for dimensionally
regularized Feynman integrals by transporting local recurrence solutions along
a line in kinematic space.  It keeps the local behavior

```text
x^(a + b eps) Log[x]^p
```

explicit, so singular endpoint limits and analytically regularized integrals
do not have to reconstruct those sectors from sampled epsilon expansions.

It is the recurrence-based successor to
[DiffExp](https://gitlab.com/hiddingm/diffexp), with exact local sectors,
analytic regularization, Feynman-trick recursion, and a compiled C++/FLINT
recurrence backend.

The release design uses the compiled C++/FLINT recurrence backend by default.
Wolfram Language still owns exact singularity analysis, resonance and log-chain
decisions, analytic continuation, matching, integration, and validation.  The
compiled code accelerates the finite-width coefficient recurrence; it is not a
second mathematical implementation.

The root `DiffExp2.m` loader selects the C++ backend by default. A missing or
unsupported compiled backend fails loudly and never falls back silently; the
Wolfram recurrence remains available as an explicit reference mode.

## What it supports

- exact epsilon-rational differential systems at ordinary and regular-singular
  points;
- denominator-cleared finite-width recurrences, resonances, log sectors,
  fractional powers, inhomogeneous sources, and rank reduction of supported
  higher-order poles;
- real-line transport with exact singularity geometry, complex-root projection
  waypoints, and explicit `+i0`/`-i0` prescriptions;
- exact sector-aware endpoint limits and line integrals compatible with
  analytic regularization;
- a Feynman-trick ladder which generates differential systems with FIRE and
  solves every level with DiffExp 2;
- a C++20/FLINT recurrence engine with strict Wolfram/C++ parity checks and
  parallel homogeneous-column solves.

## Quick start

Build the recurrence backend first; see [Installation](Docs/Installation.md).
Then run a complete one-component transport:

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

result = DiffExp2`TransportEndpoint[
  sys,
  {{1/2, 0, 0}},              (* one master, eps^0 through eps^2 *)
  0,
  1
];

DiffExp2`EpsilonCoefficient[result, 0]
(* {1/4} *)
```

The runnable version is
[Examples/Direct/MinimalTransport.wl](Examples/Direct/MinimalTransport.wl).
For individual chart segments, exact local powers, and a plot near a singular
endpoint, see
[Examples/Direct/SingularEndpointAndSegments.wl](Examples/Direct/SingularEndpointAndSegments.wl).

## Feynman trick

The supported end-to-end driver currently remains a command-line workflow:

```sh
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
FT_EXAMPLES=bubble \
FT_WORKING_PRECISION=300 \
FT_EXPANSION_ORDER=40 \
FT_BOUNDARY_EXTRA_ORDER=10 \
wolframscript -file Scripts/run_ft_stepwise2.m
```

The runner selects C++ by default. `DE2_CPP_THREADS=4` changes its native
worker budget; `DE2_RECURRENCE_BACKEND=Wolfram` explicitly selects the
reference recurrence for diagnostics. Neither backend silently falls back to
the other.

`FT_EXAMPLES` is a comma-separated list (default `bubble`); available names
are `bubble`, `sunrise`, `banana`, `box`, `pentagon`, `box_bubble`,
`box_triangle`, and `double_box_planar`.  The other variables set the working
precision, local expansion order, requested epsilon order, extra epsilon
orders retained at level boundaries, optional comma-separated per-level
epsilon halos (levels listed from 1 upward), the classic coupled
predivision/matching divisor, and the affine chart-coordinate radius,
respectively. Adjacent regular charts meet at
`+1/DivisionOrder` and `-1/DivisionOrder`; the runner overrides a differing
legacy `FT_STEP_DIVISION_ORDER` value. A halo computes extra internal
coefficients without raising the downstream requested order; it compensates
only certified finite-window losses such as analytic log/Laurent basis
width. The package facade enables the regular-chart value-transport path by
default; the low-level script retains the explicit `DE2_VALUE_TRANSPORT`
switch.

It caches FIRE preparation separately from transport-ladder checkpoints, so a
solver rerun does not have to repeat the reduction.  Start with the
[bubble example](Examples/FeynmanTrick/Bubble.sh) and read the
[Feynman-trick guide](Docs/FeynmanTrick.md) before running the larger examples.

The same workflow is available through a typed Wolfram Language facade:

```mathematica
Get[FileNameJoin[{repo, "FeynmanTrick", "FeynmanTrick.m"}]];
plan = FeynmanTrick`PipelinePlan["bubble"];
result = FeynmanTrick`RunIntegrationPipeline[plan];
```

## Documentation

- [Installation](Docs/Installation.md)
- [Quick start](Docs/QuickStart.md)
- [Direct solver](Docs/DirectSolver.md)
- [Feynman trick](Docs/FeynmanTrick.md)
- [Analytic continuation and exact local sectors](Docs/AnalyticContinuation.md)
- [Verified results](Docs/Results.md)
- [API reference](Docs/API.md)
- [Migration from DiffExp 1](Docs/Migration.md)
- [C++ backend details](Docs/CppBackend.md)
- [Release file-disposition manifest](Docs/ReleaseManifest.md)

## Important current limitations

- The input differential matrix must be exact and rational in the line
  variable and epsilon.  `LoadSystem` rejects non-integer powers with a
  line-variable-dependent base, such as `Sqrt[x]`.  This is the precise
  source-verified limitation sometimes summarized as “roots in the basis are
  not supported”; algebraic constants, algebraic singular locations, and
  algebraic local exponents are not rejected merely for being algebraic.
- Indicial eigenvalues must be affine in epsilon, `a + b eps`.  Unsupported
  structures fail rather than being sampled or fitted.
- Plotting is assembled from `LineSegments` or `PiecewiseSolution`; the package
  deliberately returns inspectable data rather than owning a plotting style.
- FIRE is required only for preparing new Feynman-trick families.  Direct
  transport of an existing exact matrix does not require FIRE.

DiffExp 2 is the successor to the original
[DiffExp](https://gitlab.com/hiddingm/diffexp).  The numerical conventions are
being preserved where they remain compatible with explicit sector data, but
DiffExp 2 deliberately does not silently fall back to the old
Wronskian/Frobenius strategy stack.
