# Persistent Native Solver

Status: implementation contract for the migration from the stateless
recurrence bridge to one persistent C++ solver session.

## Outcome

A loaded differential system owns one long-lived native session.  After the
exact system and transport plan have crossed the LibraryLink boundary, chart
marching no longer serializes complete recurrence problems or materializes
intermediate `LocalSolution` tensors in Wolfram.  The native session owns:

1. prepared regular and singular chart operators and the exact SCC graph;
2. regular, regular-singular, resonant, logarithmic, fractional-power, and
   inhomogeneous recurrences;
3. all off-diagonal SCC source propagation;
4. matching, iterative refinement, and chart-to-chart state;
5. local evaluation and residual/error certificates;
6. endpoint limits and line/tile integration;
7. versioned checkpoint state sufficient to resume without reconstructing
   completed native work.

Wolfram remains the public orchestration and exact-input layer.  It may build
or verify exact algebraic metadata before session creation, but it must not be
the per-chart numeric execution engine.

## Non-negotiable contracts

- Selecting the C++ backend is strict.  Unsupported input or a failed native
  certificate is loud; there is no automatic Wolfram-solver fallback.
- Resonance, SCC edges, epsilon valuations, and branch choices are never
  inferred from floating-point midpoints.
- Every epsilon tensor carries both its Laurent frame and its completeness
  metadata.  A shift that discards known lower-frame content is an error.
- Sector tags `(a,b,p)` remain distinct and exact.  Logarithmic sectors are
  not collapsed or numerically fitted.
- `+i0` and `-i0` are explicit chart data.  For `t=-r`, the native branch is
  `Log(t)=Log(r)+i*pi*sigma`; the same `sigma` controls the power phase and
  log-chain mixing.  Missing or conflicting prescriptions are loud.
- Analytic regularization keeps the existing endpoint rule: `b != 0`
  endpoint sectors are continued in epsilon and dropped by the limit;
  genuinely divergent `b == 0` content must cancel before the endpoint gate.
- Numeric production uses complex Arb balls.  Exact rational and exact
  analytic-regulator coefficient fields remain available where structural
  decisions require them.
- A checkpoint is accepted only when its schema, solver build, exact-system
  identity, configuration identity, scalar domain, and branch-plan identity
  agree with the requesting run.

## Session API

LibraryLink exposes a single command entry point in addition to backend
information.  Commands and responses use a versioned envelope.  Session and
object handles are opaque process-local strings:

```text
request  = {schema: 2, op: string, session?: string, ...operation fields}
response = {status: "ok"|"error", session?: string, result?: object,
            id?: string, detail?: string}
```

The implemented recurrence-session commands are:

| Command | Purpose |
| --- | --- |
| `session.create` | Allocate a session with immutable precision/domain/configuration policy. |
| `chart.prepare` | Parse and retain one chart/frame recurrence operator and its SCC certificate. |
| `chart.solve` | Run one recurrence against a retained operator. |
| `chart.solve_batch` | Run an ordered bounded batch of columns/sectors sharing one operator. |
| `session.solve_many` | Run one ordered worker pool across several retained charts in the same session. |
| `chart.release` | Release one retained operator. |
| `local.solve` | Run recurrence plus assembly and retain the typed local sector family without coefficient JSON. |
| `local.evaluate` | Evaluate retained value and theta-value with explicit prescription semantics. |
| `local.stats` | Inspect one retained local solution and its evaluation counters. |
| `local.release` | Release one retained local solution. |
| `scc.prepare` | Validate and retain a typed composite SCC chart; execution is not enabled yet. |
| `scc.stats` | Inspect one retained composite manifest and signed-shift bounds. |
| `scc.release` | Release one retained composite SCC chart. |
| `session.stats` | Return preparation/run timings and retained-object counters. |
| `session.close` | Release the session deterministically. |

The current Wolfram preparation seam is
``Solve`PrepareNativeSCCComposite[sccEnvelope, req]``.  It captures each
diagonal block's ordinary grouped homogeneous request without executing it,
requires one shared coefficient field and work rectangle, builds the exact
parent/block/coupling manifest, and calls
``CppBackend`PreparePersistentSCC[groups, manifest]``.  This first slice is
strictly limited to regular collision-free blocks with identity gauge and
spectral transforms.  It always prepares `BlockSequentialStrict`; it never
constructs a monolithic alternative and is not yet used by production
`SolveHomogeneous` or transport.  A bounded collision-checked Solve cache
reuses live handles, checks native stats before an early hit, and fails loudly
at capacity rather than evicting a public SCC object.

Subsequent milestones extend the same session with complete multi-sector
composition plus `transport.*`, `endpoint.*`, `integration.*`, and
`checkpoint.*` operations; they do not introduce a second stateless
execution path.

All commands are serialized per session.  Independent sessions may run in
parallel.  A session may execute independent lower/upper endpoint arms in its
own bounded worker pool because no Wolfram evaluation is needed between
native chart steps.

## Native state

```text
SolverSession
  Identity
    schema/build/system/config/plan hashes and collision-proof source records
  Domain
    precision, output digits, analytic-regulator symbols, chop policy
  System
    exact SCC graph, topological order, components, coupling matrices
  Charts[id]
    geometry and branch metadata
    cleared theta operator and denominator groups
    regular/singular indicial and Jordan data
    T/P/R schedules and epsilon valuations
    prepared V, V^-1, and gauge operators
    prepared multiplier and integration kernels
  TransportStates[id]
    current chart, local sector slabs, incoming value, error certificate
  Results
    completed arms, endpoint limits, line/tile integrals
  Counters
    preparation, recurrence, matching, residual, integration, allocation
```

Chart IDs are content addressed by the complete exact chart input, not only a
digest.  Hashes index the registry; the stored source record proves equality
and makes collisions loud.

## Data model

### Framed epsilon tensors

The existing C++ `Frame`, `FrameBlock`, validity tensor, and strict lower-frame
guards become the common storage for all numeric stages.  Operations consume
explicit `(frame_base, frame_width, complete_max)` records.  Intermediate SCC
sources never round-trip through JSON.

### Local sectors

```text
LocalSector
  ExactTag a, b
  uint32 log_power
  Dense coefficient slab [epsilon][Taylor][component]
  per-row/per-component completeness

LocalSolution
  chart_id, radius, prescriptions
  epsilon and Taylor windows
  vector<LocalSector>
  tail/error certificate
```

Regular charts use the same representation with the single `(0,0,0)` sector.
This avoids a second transport representation and permits regular-to-singular
and singular-to-regular handoffs without conversion through Wolfram.

### Analytic-regulator numeric coefficients

The current exact `SymbolicRational` domain is sufficient for recurrence
preparation but not for all matching/evaluation operations: evaluating
`t^(a+b eps)` introduces Arb constants such as `Log(t)` and `pi` multiplying
rational functions of extra regulators.  The persistent engine therefore
needs a coefficient domain whose elements are rational functions in declared
regulators with complex-ball coefficients, or an equivalent expression DAG
closed under the required finite arithmetic.  Regulator sampling is not an
acceptable substitute.  Until this domain exists, native matching of formal
analytic-regulator runs is incomplete even if numeric runs work.

The current `ComplexBall` precision and `SymbolicRational` FLINT variable
context are process-global.  A persistent multi-session implementation must
replace them with immutable per-session arithmetic contexts, or serialize all
incompatible domains behind a verified context guard.  Merely storing two
handles while allowing one session to reconfigure live values from another is
not safe.

## Execution pipeline

### 1. Creation and preparation

Wolfram sends immutable configuration and an exact system manifest once.
C++ validates dimensions and exact SCC metadata, stores the full source record,
and prepares every chart operator once.  Repeated segments reference IDs.

### 2. Recurrence and SCC propagation

The session runs diagonal SCC blocks in topological order.  For a target block
it forms all incoming coupling sources directly from retained predecessor
slabs, applies prepared rational matrices, combines sources in their honest
window intersection, and invokes the regular/singular recurrence.  No JSON
request is built per block and no source tensor is decoded in Wolfram.

### 3. Matching and refinement

The producer and receiver are evaluated at the common geometric match point
inside C++.  Laurent-graded elimination preserves formal valuations.  At most
the configured number of deterministic refinement steps solves the residual
correction using the same factorization.  The accepted state stores the
residual certificate and the error probe; only compact progress metadata is
returned to Wolfram.

### 4. Certification

Local evaluation, theta differentiation, and matrix application share retained
operator data.  A completed SCC solve receives one full original-system ODE
certificate.  Implementation-intermediate diagonal blocks are not redundantly
certified.  Arb enclosures and completeness windows are both part of the
certificate.

The certificate separates mechanisms instead of collapsing them into one
dimensionless list:

```text
ErrorCertificate
  work_window, delivered_window
  acb_roundoff_upper[]
  matching_residual[]
  ode_residual[]
  local_tail_absolute[], local_tail_relative[]
  propagated_tail_upper[], reference_scale[]
  status = certified | heuristic | failed
```

The existing full-vs-reduced-order difference at `0.51 Radius` is recorded as
`heuristic`, keyed by its actual epsilon powers.  It must not be list-position
added across differing Laurent minima or trigger an unconditional
`absolute > 1` abort.  A rigorous truncation certificate ultimately needs a
recurrence majorant and chart-transfer sensitivity:

```text
E_out <= abs(transfer_operator) * E_in + E_local.
```

### 5. Endpoint and integration

Endpoint limits operate on already-combined native sectors.  Line integration
owns the certified half-radius tiling, prepared rational multipliers, and
closed-form monomial primitives.  Identical primitive towers and multiplier
kernels are retained across masters and components.  Lower and upper arms may
execute concurrently under one global native worker budget.

The level-facing operation accepts all lower-master observables together.  It
tiles each arm once, prepares each distinct rational multiplier once, builds
each primitive tower once, and contracts every requested observable before
returning the next boundary vector.  Calling a single-master line integral in
a Wolfram loop is a compatibility wrapper, not the production native path.

### 6. Checkpoints

The checkpoint is a binary, endian-stable container with a small JSON header
and checksummed typed sections.  At minimum it contains:

- schema/build/domain/precision/configuration identity;
- exact system and branch-plan identity records;
- prepared SCC and chart operator tables, or content-addressed references to
  an independently validated preparation cache;
- completed transport-arm states and their certificates;
- endpoint and line/tile results already completed;
- counters needed to distinguish resumed work from recomputation.

Writes use `temporary -> fsync -> atomic rename`.  Restore validates every
identity before exposing a session handle.  Unknown mandatory sections fail;
unknown explicitly optional sections may be skipped.

## Migration milestones

1. Persistent session registry and retained prepared recurrence operators.
2. Retained chart/SCC graph plus native target-block source propagation.
3. Native local-sector evaluation, matching, refinement, and ODE certificate.
4. Native endpoint limits and line/tile integration.
5. Native checkpoint/restore and concurrent endpoint arms.
6. Remove the stateless full-problem recurrence path from production; retain
   it only as a focused reference/parity seam.

Each milestone must make the production call graph more native.  A wrapper
that stores JSON but reconstructs the same complete Wolfram-side source and
result tensors on every chart does not satisfy a milestone.

## Validation and performance gates

Validation is deliberately staged rather than repeated after every edit:

1. one tiny scalar regular/singular C++ smoke;
2. sunrise;
3. box;
4. box with bubble;
5. one saved difficult L3 singular segment with before/after phase timings;
6. L2/double-box only after the previous gates pass;
7. one final focused regression suite.

For each case compare sector tags, epsilon/Taylor windows, prescriptions,
endpoint or integral coefficients, and certificate status.  Performance
reports separate one-time preparation, recurrence, SCC propagation, matching,
certification, integration, serialization, and checkpoint I/O.  The decisive
metric is end-to-end warm native wall time, not only recurrence-kernel time.

## Completion evidence

The migration is complete only when all seven ownership items in **Outcome**
have a production C++ implementation, the normal C++ run performs no
per-chart Wolfram recurrence/matching/evaluation/integration, analytic
regularization and both delta branches pass parity cases, checkpoint restore
does not recompute completed work, and the staged examples complete with
recorded timings and precision.
