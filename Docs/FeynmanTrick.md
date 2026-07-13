# Feynman Trick

The Feynman-trick workflow recursively combines propagators,

```text
1/(D_i^a D_j^b)
  = Gamma[a+b]/(Gamma[a] Gamma[b])
    Integral_0^1 dx x^(a-1) (1-x)^(b-1)
    / (x D_i + (1-x) D_j)^(a+b),
```

until the deepest level has analytic tadpole-type boundary data.  FIRE builds
the master bases, IBP reductions, and exact differential matrix at each level.
DiffExp 2 transports that data to both endpoints and applies a sector-aware
limit or integral to obtain the next level's boundary.

## Public pipeline facade

Load the root package and inspect the exact registry names:

```mathematica
Get["/absolute/path/to/diffexp2/FeynmanTrick.m"];
FeynmanTrick`$FeynmanTrickVersion
FeynmanTrick`SupportedExamples[]
```

Build a reproducible plan without starting FIRE or another kernel, then run
that same plan:

```mathematica
plan = FeynmanTrick`PipelinePlan["bubble",
  "WorkingPrecision" -> 300,
  "ExpansionOrder" -> 40,
  "BoundaryExtraOrder" -> 10,
  "CppThreads" -> 4
];

result = FeynmanTrick`RunIntegrationPipeline[plan];
```

The result has schema `FeynmanTrick.PipelineResult/v1` and records status,
exit code, parsed `FINAL`/`STEPWISE` rows, stdout/stderr, and the complete plan.
Resume an atomic checkpoint with:

```mathematica
result = FeynmanTrick`ResumeIntegrationPipeline[
  "bubble", "/absolute/path/to/bubble_level0_boundary.mx"
];
```

`"Asynchronous" -> True` returns a
`FeynmanTrick.PipelineProcess/v1` record containing the process object and
plan. A run accepts either one validated registry name or one exact
user-defined family. The underlying registry CLI also accepts comma-separated
names.

### User-defined families and output selection

Pass an exact family together with one index vector, an ordered list of index
vectors, or `All`:

```mathematica
family = <|
  "Name" -> "massive_bubble_custom",
  "LoopMomenta" -> {l},
  "ExternalMomenta" -> {p},
  "Propagators" -> {1 - l^2, 3/2 - (l - p)^2},
  "Replacements" -> {p^2 -> s},
  "NumericalPoint" -> {s -> -1},
  "Dimension" -> 2 - 2 FeynmanTrick`FTeps
|>;

selected = FeynmanTrick`RunIntegrationPipeline[
  family, {{1, 1}}, "EpsilonOrder" -> 2
];

allMasters = FeynmanTrick`RunIntegrationPipeline[
  family, All, "EpsilonOrder" -> 2
];
```

`CreateFamily` validates and canonicalizes the mathematical input without
starting FIRE. `PipelinePlan[family, targets]` then creates an exact,
content-addressed WXF request without writing it; execution publishes that
request atomically before starting the child. Explicit targets retain their
order and multiplicity. `NumericalPoint` is applied both to scalar-product
replacement rules and to symbols occurring directly in propagators, so exact
symbolic masses may be fixed there. `All` performs deterministic level-zero
FIRE master discovery at execution and caches a source- and runtime-bound
resolution.
The parent independently validates every returned target identity and, for
`All`, the discovery manifest.

The facade runs `Scripts/run_ft_stepwise2.m` in a clean `wolframscript`
subprocess. Its argv and environment are explicit in `PipelinePlan`; no shell
string is constructed. The current process launcher uses `/usr/bin/env`, so
the facade is presently supported on macOS/Linux/Unix. The direct DiffExp 2
solver itself has no such platform restriction. Calling the facade from an
already running kernel can occupy a second Wolfram license seat; its native
C++ worker threads do not consume additional seats.

## Command-line alternative

```sh
DE2_RECURRENCE_BACKEND=Cpp \
DE2_CPP_THREADS=4 \
FT_EXAMPLES=bubble \
FT_WORKING_PRECISION=300 \
FT_EXPANSION_ORDER=40 \
FT_EPS_ORDER=0 \
FT_BOUNDARY_EXTRA_ORDER=10 \
FT_DIVISION_ORDER=3 \
FT_RADIUS_OF_CONVERGENCE=1 \
FT_PREP_CACHE_DIR="$HOME/.cache/diffexp2/fire" \
FT_LADDER_CHECKPOINT_DIR="$HOME/.cache/diffexp2/ladder/bubble" \
wolframscript -file Scripts/run_ft_stepwise2.m
```

This is the command generated conceptually by the facade and remains useful
for terminals and batch systems. The first run prepares FIRE data and writes a
reusable preparation snapshot.
Later runs print `FTPREP CACHE HIT` and skip that work.

## Built-in Euclidean examples

The current fixtures in `Scripts/FTExamples.m` are:

| Name | Integral and point | Dimension |
| --- | --- | --- |
| `bubble` | equal-mass one-loop bubble, `p^2=-1` | `2-2 eps` |
| `sunrise` | equal-mass two-loop sunrise, `p^2=-1` | `2-2 eps` |
| `banana` | equal-mass three-loop banana, `p^2=-1` | `2-2 eps` |
| `banana_unequal` | squared masses `{2,3/2,4/3,1}`, `p^2=-1` | `2-2 eps` |
| `banana4` | equal-mass four-loop banana, `p^2=-1` | `2-2 eps` |
| `banana4_unequal` | squared masses `{2,3/2,4/3,5/4,1}`, `p^2=-1` | `2-2 eps` |
| `kite` | fully massive equal-mass two-loop kite, `p^2=-1` | `2-2 eps` |
| `box` | massless on-shell box, `s=-1`, `t=-1/3` | `4-2 eps` |
| `box_bubble` | massless two-loop three-point subfamily, same Euclidean invariants | `4-2 eps` |
| `box_triangle` | massless two-loop four-point subfamily | `4-2 eps` |
| `double_box_planar` | massless planar double box | `4-2 eps` |
| `pentagon` | massless pentagon at fixed negative adjacent invariants | `4-2 eps` |
| `pentagon_massive` | fully massive pentagon at a symmetric Euclidean point | `4-2 eps` |

Release examples are in [Examples/FeynmanTrick](../Examples/FeynmanTrick/README.md).
`bubble`, `sunrise`, and `banana_unequal` are the recommended progression.
The four-loop unequal banana is intentionally marked experimental because a
complete, source-controlled DiffExp 2 ladder result is not yet present.

## Facade defaults

The typed facade defaults to C++, working precision 500, expansion order 50,
epsilon order 0, boundary lookahead 4, halos `{0}`, division order 3, radius
1, value transport enabled, and endpoint-arm batching requested. Native thread
count is `Min[10,$ProcessorCount]` with a floor of one. Every setting can be
overridden by a named `PipelinePlan` option and is serialized into its
`"Environment"` record.

## Numerical controls

| Environment variable | Meaning |
| --- | --- |
| `DE2_RECURRENCE_BACKEND` | `Cpp` for the release path, or explicit `Wolfram` reference mode |
| `DE2_CPP_THREADS` | native recurrence worker budget; no extra Wolfram kernels or licenses |
| `FT_WORKING_PRECISION` | working precision for transport and level handoffs |
| `FT_EXPANSION_ORDER` | Taylor order retained in every local chart |
| `FT_EPS_ORDER` | highest final epsilon coefficient requested |
| `FT_BOUNDARY_EXTRA_ORDER` | extra internal epsilon lookahead at each level |
| `FT_LEVEL_EPS_HALOS` | comma-separated extra lookahead by level, listed from level 1 upward |
| `FT_DIVISION_ORDER` | coupled chart placement/matching divisor; adjacent regular charts meet at `+1/k` and `-1/k` |
| `FT_RADIUS_OF_CONVERGENCE` | affine local-coordinate radius normalization |
| `FT_FIRE_TIMEOUT_SECONDS` | watchdog timeout for one FIRE run |
| `FT_FIRE_PATH` | FIRE6 installation path; facade option `"FIREPath"` |
| `FT_CPP_BATCH_ENDPOINT_ARMS` | request paired C++ prewarming for small lower/upper chart bases; CLI default on |
| `DE2_VALUE_TRANSPORT` | regular-chart value transport; facade default on, direct CLI default off |

More precision does not replace sufficient expansion order.  Conversely, a
large epsilon halo does not increase the requested final order; it only
supplies internal coefficients consumed by Laurent shifts, resonances, and
endpoint integrations.  If an honest window is too short, the runner stops
with `FTLADDER INCOMPLETE` instead of padding it.

At each level the runner also checks the exact FIRE differential matrix for
epsilon poles.  When needed it transports the diagonal basis
`J_i = eps^k_i I_i`, using `FindEpsPrefactors` and `ApplyEpsPrefactors`, and
converts the incoming boundary arrays and IBP coefficients by the same
per-master powers.  A common shift keeps every boundary conversion finite;
the stored complete epsilon ceiling is not extended.  Consequently a relative
basis shift can consume upper physical orders.  Supply those orders with
`FT_BOUNDARY_EXTRA_ORDER` or `FT_LEVEL_EPS_HALOS`; the runner reports an
incomplete handoff rather than assuming the missing coefficients vanish.

## Preparation cache and ladder checkpoints

`FT_PREP_CACHE_DIR` stores completed Feynman-trick/FIRE preparation and the
reduction cache.  Its v3 contract records the exact topology, combination
sequence, dimension, preparation-affecting configuration, Wolfram/FIRE
runtime, and repository-relative hashes of the three preparation modules:
`PropagatorAlgebra.m`, `FIREInterface.m`, and `FeynmanTrickIteration.m`.
Custom-family contracts additionally hash `FamilySpec.m` and
`PipelineRequest.m` and bind the exact output selection.
`FeynmanTrick.m` configuration semantics are represented by their evaluated
contract values; `MatrixExport.m` writes artifacts only; and
`LevelReduction.m` derives transport-time requests whose exact hardened
reduction keys are revalidated on every load.  Consequently facade, runner,
boundary, LevelReduction, DiffExp2 solver, and transport-only changes do not
force FIRE to rerun.  A snapshot stores and exactly rechecks the full contract
and all required reduction keys, rather than trusting its digest alone.

Legacy v1/v2 snapshots are rejected by default because they do not contain
that source contract.  A one-time migration can be requested explicitly with
`FT_MIGRATE_LEGACY_PREP=1` (or the facade option
`"MigrateLegacyPreparation" -> True`).  Migration succeeds only when exactly
one candidate passes the current topology, sequence, configuration,
prepared-level, and reduction-key checks; otherwise it rejects the snapshot
and preparation runs normally.  V2 must additionally carry FIRE setup records
matching the current runtime.  V1 predates those records and hardened keys: its
old tuple keys must identify exactly one retained level topology, are
normalized and re-keyed with the full current fallback topology/configuration
record, and must cover every current boundary request without a nonidentical
collision.  The migrated v3 snapshot records that its historical FIRE setup
and source provenance remain unverified.  Exact re-keyed cache hits never
launch FIRE; any future missing reduction is still rejected because no
verified setup record was fabricated.

`FT_LADDER_CHECKPOINT_DIR` stores two kinds of checkpoint:

- a transport checkpoint after each completed lower or upper endpoint arm;
- a boundary checkpoint after the next level's boundary vector is assembled.

Resume explicitly:

```sh
FT_RESUME_LADDER_CHECKPOINT=/absolute/path/to/example_level2_boundary.mx \
DE2_RECURRENCE_BACKEND=Cpp \
FT_EXAMPLES=example \
wolframscript -file Scripts/run_ft_stepwise2.m
```

The checkpoint records the example, backend, precision, expansion and epsilon
settings, prepared-data key, level metadata, the exact per-master epsilon
basis and normalized-matrix hash, and a source fingerprint. Custom-family
checkpoints additionally bind the family, pipeline request, and any resolved
`All` selection identities.
Stale/unversioned checkpoints are rejected unless
`FT_ALLOW_STALE_LADDER_CHECKPOINT=1` is set; descendants of an explicitly
accepted stale checkpoint remain marked tainted.

## Lower and upper endpoint scheduling

A single Wolfram kernel cannot execute two Wolfram marching loops at once.
With `FT_CPP_BATCH_ENDPOINT_ARMS=1`, however, small boundary-independent
homogeneous bases from the two arms can be submitted together to the native
worker pool.  Planning, matching, analytic continuation, integration, and
checkpoint writes remain sequential. The runner skips this schedule when a
single chart already fills the thread budget or when value transport makes the
recurrence boundary-dependent. Therefore the facade requesting both value
transport and arm batching does not imply that homogeneous prewarming runs.

## Output

The runner prints machine-readable records:

```text
STEPWISE {"Example":...,"Level":...,"Master":...,"RawMinPower":...,"Coefficients":...,"Certification":...}
FINAL {"Example":...,"Finite":...,"RawMinPower":...,"Certification":...}
```

Custom-family runs emit one ordered `FINAL` row per requested integral, with
its family, pipeline-request, physical-integral, and ordinal identities. An
`All` run first emits exactly one `OUTPUT_RESOLUTION` manifest; the parent
requires the final rows to match that manifest. The typed result exposes all
accepted rows through `result["Outputs"]`.

Use the `STEPWISE` rows for Laurent coefficients and intermediate-level
audits.  Historical logs used `"Finite"` inconsistently when a result had
poles; [Results](Results.md) reports coefficients explicitly and does not rely
on that label alone.  Native integral rows carry the exported `Scope`,
`ErrorGuarantee`, and `ErrorEnvelope`; a `stored_truncation` scope is reported
honestly with no guarantee and is not a run failure.  Endpoint limits, direct
handoffs, proved-zero observables, and the initial deepest boundary instead
carry an explicit `"Applicability":"not-applicable"` record.  The terminal
`FINAL` row repeats the certification of its corresponding level-zero master.

## Custom-family limits

The public family path requires exact propagators, replacement rules,
kinematic points, dimensions, and a complete merge sequence ending in one
active propagator. `All` discovery fails before transport when FIRE introduces
additional numerator slots. Any requested or discovered numerator at a merge
position is rejected, while auxiliary slots introduced during explicit-target
preparation remain a lower-level implementation detail and are not advertised
as `All` outputs. Custom analytic-prescription and kinematic-assumption fields
are also rejected until their branch/assumption semantics are wired into the
production request path.

The lower-level `DefineTopology`, `DefineFTIteration`, and `RunFullIteration`
objects remain available for development work, but copying runner-private
helpers into notebooks is not recommended.

## Failure and reproducibility rules

- FIRE failure, timeout, or incomplete master data aborts preparation.
- A requested C++ backend never falls back silently.
- Missing analytic continuation prescriptions for material multivalued
  sectors are errors.
- Incomplete epsilon windows are errors.
- Exact propagator definitions, squared-mass conventions, dimensions, and
  Euclidean points must be recorded with every result.
- Warm transport timings must be distinguished from cold runs that include
  FIRE preparation.

Verified values and timing scopes are summarized in [Results](Results.md).
