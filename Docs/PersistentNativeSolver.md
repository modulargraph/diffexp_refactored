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
| `local.evaluate` | Evaluate retained value and theta-value with explicit prescription semantics; an exact `certified_tail_radius_exact` may request the attached regular-tail proof. |
| `local.certify_residual` | Certify the retained local's stored-truncation theta residual with native Arb magnitudes. |
| `local.match` | Exact-match retained regular locals and retain the Laurent lattice transformation and weights without coefficient JSON. |
| `local.match_acb` | Evaluate retained Acb locals at an exact common physical point, apply an exact-Rational saturation transformation, and retain bounded-refinement residual diagnostics without coefficient JSON. |
| `local.endpoint_limit` | Apply the native sector endpoint gate to one retained local and retain the specialized endpoint vector without coefficient JSON. |
| `local.stats` | Inspect one retained local solution and its evaluation counters. |
| `local.release` | Release one retained local solution. |
| `endpoint.stats` | Inspect one retained endpoint result and its exact analytic/branch provenance. |
| `endpoint.export` | Explicitly export a final specialized Acb epsilon vector for compatibility. |
| `endpoint.release` | Release one retained endpoint result. |
| `tile.plan` | Build and retain immutable exact lower/upper arm plans from retained chart geometry, branch topology, and `DivisionOrder`. |
| `tile.match_interval` | Read one exact physical handoff and its opposite-sign producing/receiving local coordinates from a retained plan. |
| `tile.integration_interval` | Read one exact physical/local tile interval and its affine Jacobian from a retained plan. |
| `tile.match_advance` | Derive one Rational/Acb match entirely from a retained plan and retain the receiving-basis weights with strong plan/local ownership. |
| `tile.stats` / `tile.release` | Inspect or release one retained independent-arm plan. |
| `integration.line` | Integrate one retained local over a plan-selected tile, apply the exact affine Jacobian, and retain the physical epsilon vector; optional `certify_tail:true` promotes only after proof. |
| `integration.stats` / `integration.export` / `integration.release` | Inspect, explicitly export, or release one retained tile integral. |
| `match.stats` | Inspect one retained native match and its exact provenance. |
| `match.materialize_local` | Apply a plan-driven match's retained Laurent weights to its receiving basis and publish the next retained local without coefficient JSON. |
| `match.release` | Release one retained native match state. |
| `scc.prepare` | Validate and retain a typed composite SCC chart. |
| `scc.solve_column` | Execute one retained regular exact-Rational or Acb SCC basis column, or an exact-Rational regular-singular Jordan column, and retain the parent local without coefficient JSON. |
| `scc.solve_columns` | Execute an ordered SCC basis-column batch in one bounded native worker pool and publish every retained local atomically, or none on failure. |
| `scc.stats` | Inspect one retained composite manifest and signed-shift bounds. |
| `scc.release` | Release one retained composite SCC chart. |
| `session.stats` | Return preparation/run timings and retained-object counters. |
| `session.close` | Release the session deterministically. |

The current Wolfram preparation seam is
``Solve`PrepareNativeSCCComposite[sccEnvelope, req]``.  It captures each
diagonal block's ordinary grouped homogeneous request without executing it,
requires one shared coefficient field, takes the exact union of each
diagonal block's independently proved work rectangle, builds the exact
parent/block/coupling manifest, and calls
``CppBackend`PreparePersistentSCC[groups, manifest]``.  This first slice is
limited to regular or exact affine-Jordan blocks with identity
gauge and spectral transforms.  C++ reconstructs and retains the complete
Rational Jordan indicial certificate; `regular` is a classification rather
than an admission predicate. It always prepares `BlockSequentialStrict`; it
never
constructs a monolithic alternative and is not yet used by production
`SolveHomogeneous` or transport.  A bounded collision-checked Solve cache
reuses live handles, checks native stats before an early hit, and fails loudly
at capacity rather than evicting a public SCC object.

The explicit
``Solve`SolveNativeSCCBasisColumn[sccEnvelope, req, seedBlock,
seedLocalComponent:1]`` seam consumes the compact run-only data in that cache
for regular exact-Rational or Acb scalar/multidimensional scopes and the
exact-Rational regular-singular scalar/Jordan scopes. Every capture retains
exact producer-side `(a,b,P)` task metadata; roots and log ceilings are never
reconstructed from Acb values. The seam derives every component from its
flattened eps^0 unit tensor, checked against independently encoded exact
zero/one values in the selected coefficient field, selects the requested seed
without trusting capture order, orders reachable targets by the exact SCC
certificate, and propagates the exact log ceiling through the DAG.
Each target's zero-initial particular request and complete T/P/R schedule are
built once from that target's exact Jordan data, then independently
re-certified by C++. Regular Acb runs are separately restricted to exact
`a=b=0`, `p=0`, exact integer Taylor shifts, and the retained singleton
schedule; numerical enclosures never decide a resonance or pseudo case.
Native coupling and
recurrence operate on retained `LocalSolution` values and return only an
opaque local handle plus exact column provenance.  Exact parent geometry and
analytic prescriptions are preserved.  The three-argument scalar-v1 wire
request and returned record are unchanged, and neither scope is connected to
production `SolveHomogeneous` or transport dispatch.

``Solve`SolveNativeSCCBasis[sccEnvelope, req, threads]`` builds the same
certified request for every physical basis index, submits them through one
ordered `scc.solve_columns` worker pool, and returns opaque local handles
sorted by that physical index. Native retention is atomic: a failed column
publishes none of the batch, while successful batches never return coefficient
slabs. This complete-basis seam is the receiving state used by native
plan-driven matching; automatic production transport dispatch is still being
migrated.

For the incoming side, ``Solve`SolveNativeValueRegular`` turns one honest
center-value vector directly into a retained full-system local. When its chart
was prepared as an SCC skeleton, it materializes the exact identity physical
frame for that one value recursion while preserving the parent's native
session identity. Thus a boundary seed and later retained SCC bases coexist in
one C++ session without constructing or exporting a Wolfram fundamental
matrix.

``Solve`SolveNativeRegularBasis`` is the chart-generic receiving-basis entry
point. It selects the atomic SCC column batch for a multi-block envelope and
uses exact eps^0 unit seeds in one retained full-system chart when the regular
dependency graph is a single strongly connected block. Both forms expose the
same ordered opaque-column contract to transport orchestration.

The native protocol advertises both
`exact-rational-regular-singular-scalar-block-dag-column-v1` and
`exact-rational-regular-singular-jordan-block-dag-column-v2`. At least one
diagonal block is singular; every block has exact identity gauge/assembly, and
every active coupling has center pole order zero so it preserves the exact
`(a,b,p)` sector family. `no_pseudo` is retained as truthful producer
provenance, not used as an admission decision. Couplings may
have a nonzero constant Taylor coefficient, allowing an exact source resonance
to create a new log member, and signed epsilon shifts retain the ordinary
strict lower-frame guard.

This capability does not trust a caller-consistent schedule. For each
retained Rational chart, C++ reconstructs the complete affine Jordan operator
from exact `Nhat_0 d_0^-1`, accepting only affine epsilon support and proving
the declared block partition, equal roots, unit superdiagonal, and zero
off-block entries. Submitted captured steps must then satisfy exactly
`da=a_target+n-a_i` and `db=b_target-b_i`; the seed's resonant `n=0` step
therefore binds its tag to the retained root, and every descendant tag is
also checked against its in-memory predecessor source. Exact CASE-P hits are
compensated by certificate-derived homogeneous Jordan columns; their negative
Laurent weights are proved over the stored Taylor overlap, and the resulting
exact `(a,b)` sectors are split, propagated, and recombined through every
downstream SCC. Incomplete polar windows, cyclic compensation dependencies,
nonidentity spectral frames, and center-pole tag shifts remain loud errors.
The Wolfram producer executes both scalar and multidimensional Jordan
scopes through the explicit seam; production dispatch is unchanged.

The exact-rational `local.match` migration seam consumes retained regular
locals from one session. It proves that the two chart coordinates name the
same rational physical point, performs finite Laurent lattice saturation and
matching in C++, verifies an exact stored-truncation residual through the
requested CompleteMax, and retains the transformation and weights behind an
opaque match handle.  The separate `local.match_acb` operation leaves that v1
wire path unchanged.  It receives an exact-Rational evaluated-lattice witness,
derives and retains its saturation transformation `T` natively, then evaluates
the retained Acb basis and incoming locals at coordinates proved to name one
exact rational physical point.  Acb never decides lattice support, zero
structure, or resonance: every solve pivot must exclude zero as an enclosure.
One certified factorization is reused for a bounded number of residual
corrections, and the opaque match records pass/fail/inconclusive diagnostics
with its honest epsilon CompleteMax.  Exact chart/checkpoint identities,
prescriptions, requested/effective rim signs, and the canonical witness
binding are part of the retained provenance.  The low-level Wolfram entry is
`RunPersistentAcbLocalMatch`; `match.stats` and `match.release` are shared by
both match kinds.

`local.endpoint_limit` is a separate retained-result lifecycle.  Admission
strongly owns its source local, so concurrent `local.release` cannot invalidate
an in-flight limit, and session close prevents a completed in-flight result
from being published.  The request must bind both checkpoint identities, the
exact approach direction, an optional `+1`/`-1` rim, and one of two explicit
cancellation policies: `exact-coefficient-field` or
`exact-or-acb-singleton`.  The latter accepts an Acb divergent sum only when
it is the exact singleton zero; no smallness tolerance is supported.  Any
regulator slope certified exactly nonzero is dropped by dimensional analytic
continuation, including an exact symbolic regulator tag without sampling.
Unresolved symbolic coefficient fields, regulator slopes without an exact
zero decision, and non-rational unregulated powers fail loudly.  The opaque
result records every sector tag, prescription, effective rim, and source
checkpoint; coefficients cross the bridge only through the explicit
final-compatibility `endpoint.export` call.

`tile.plan` is the first retained native path/integration lifecycle.  The
request names retained prepared charts only; C++ reads their authoritative
exact center, scale, radius, and full prescription records, then runs the
exact path planner.  It retains lower and upper arms as separate immutable
snapshots sharing one strongly owned anchor.  The usual handoff is exactly
`+1/DivisionOrder` in the producing local coordinate and
`-1/DivisionOrder` in the receiver (with signs reversed on a reversed arm);
unsafe or forbidden geometry is resolved or rejected by the native exact
planner.  Complex `Re(z)` and `Re(z) +/- |Im(z)|` projections and every branch
sheet remain explicit topology rather than midpoint guesses.

The plan itself feeds both `tile.match_interval` and
`tile.integration_interval`; Wolfram does not recompute these coordinates.
`integration.line` accepts no endpoint coordinates.  It selects the retained
tile, verifies the local's chart and checkpoint identity plus its complete
prescription record, integrates the stored Taylor truncation, and applies the
exact `dx = scale dt` Jacobian with path orientation.  Without
`certify_tail:true`, the opaque line result remains `StoredTruncation` with
`ErrorGuarantee -> None` as before.

Eligible direct Rational and Acb `local.solve` calls now attach the regular
homogeneous Gronwall/Cauchy model described in
[Native local-solution core](NativeLocalSolutionCore.md).  A
`local.evaluate` request can provide
`options.certified_tail_radius_exact` as an exact positive rational.  Its
value and theta frames carry `Certified` error envelopes only if that disk is
proved; otherwise the top-level `tail_certificate` is explicitly
`inconclusive` or `unsupported` and no envelope is published.  `local.stats`
records the attached-model status and request outcome counters.

For `integration.line` with `certify_tail:true`, C++ derives an exact witness
radius strictly between the tile's outer local endpoint and its exact chart
radius.  A successful proof promotes the result to
`FullLocalWithCertifiedTail`, then multiplies both coefficients and error
bounds by the exact physical Jacobian.  Failed/unsupported attempts remain
`StoredTruncation` and retain the non-promotion reason.  The integrator still
rejects a nonempty input error envelope, incomplete epsilon halo,
nonintegrable unseen center monomial, and merely-small Acb cancellation.
Coefficients cross JSON only through `integration.export`.

Tail-model checkpoint serialization is deliberately not part of schema v2.
Local summaries expose `checkpoint_serialized:false`; a restored local reports
the model unavailable and must be re-solved before certified-tail promotion.
This is explicit non-promotion, never reconstruction from presentation
digits or fallback to the advisory last-column estimate.
Successful results strongly own their local and plan snapshots, so releasing
the public source handles does not invalidate inspection or export.

Subsequent milestones extend the same session with complete multi-sector
composition plus `transport.*`, `endpoint.*`, `integration.*`, and
`checkpoint.*` operations; they do not introduce a second stateless
execution path.

Registry admission and publication are serialized per session.  Immutable
native work executes outside that lock, so lower and upper tile integrations
in one session may run concurrently without Wolfram evaluation between chart
steps.  Independent sessions may also run in parallel.

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

The implemented schema-2 checkpoint covers quiescent prepared chart operators
and composite SCC graphs, retained Rational/Acb local solutions, ordinary and
plan-driven exact-rational or exact-lattice-guided Acb matches,
match-materialized receiving locals, endpoint results, exact tile plans, and
completed line results. Local tensors, exact
tag specializations, pseudo-hit state, error envelopes, match
weights/residuals, and every coefficient-level refinement diagnostic use exact
FLINT/Arb dump encodings; presentation decimals are never used for restart
state. The serialized ownership closure includes released locals,
planned-match hops, plans, and prepared charts still strongly owned by a
retained match, materialized local, or line, while a separate visibility
manifest prevents those dependency-only handles from becoming public again
after restore. Plan-driven records retain their typed exact/Acb match payload,
exact handoff, branch prescription, plan/checkpoint provenance, and
materialization owner lineage; an unknown embedded match kind fails instead
of degrading to another solver path. Restore reproduces stable scoped handles,
runtime and next-handle counters, exact branch sheets/rims, and exact
source/result/checkpoint provenance before exposing the new session. Symbolic
locals remain the sole deferred scalar kind; pending tile/line operations still
prevent a quiescent checkpoint.

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
