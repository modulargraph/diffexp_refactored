# Persistent native solver protocol

## Scope of schema 2

Schema 2 is the lifecycle boundary for the persistent C++ engine.  Its first
milestone retains, as typed FLINT/Arb objects, the immutable recurrence data
for one chart and epsilon frame:

- cleared `d` and `Nhat` lags;
- rational denominators and exact epsilon valuations;
- inverse `d0`, Jordan blocks, and the prepared assembly matrix;
- the exact SCC structural graph, partition, condensation DAG, and depth;
- analytic-regulator and branch-prescription identity.

The recurrence seed, source frames, exact T/P/R schedule, and target
exponents remain run data.  A run holds a read-only view of the retained
operator; no high-precision static tensor is copied into a run.

Schema 1 remains supported unchanged.

## Lifecycle

All commands are sent through the existing `de2RunRecurrence` LibraryLink
entry point.

### Create a session

```json
{
  "schema": 2,
  "op": "session.create",
  "domain": "acb",
  "precision_bits": 1694,
  "output_digits": 520,
  "chart_capacity": 256,
  "analytic": {
    "regulators": [],
    "branch_policy": "prescription-specialized"
  }
}
```

The result contains an opaque process-local handle such as `s:1`.  Acb
precision is installed per worker through a guarded thread-local lease, so
sessions at the same precision may run concurrently and cached sessions at
different precisions may coexist.  Operations requiring incompatible live
precision state are serialized by the lease; changing Wolfram
`WorkingPrecision` does not require manual cache destruction.

For `domain: "symbolic"`, `symbols` is mandatory and names the exact
rational-function coefficient field.  A different field cannot be selected
while any symbolic values in the previous field remain alive.

### Prepare a chart/window operator

```json
{
  "schema": 2,
  "op": "chart.prepare",
  "session": "s:1",
  "key": "system/chart/frame cache key",
  "identity": "collision-certified canonical identity",
  "analytic": {
    "delta_sign": -1,
    "branch_signature": "...",
    "geometry": {
      "center_exact": "0",
      "scale_exact": "1",
      "radius_exact": "2",
      "infinite_radius": false,
      "prescriptions": []
    }
  },
  "scc": {
    "components": [[0, 2], [1]],
    "structural_edges": [[0, 2], [2, 0], [0, 1]],
    "condensation_edges": [[0, 1]],
    "topological_order": [0, 1],
    "coupling_depth": 1
  },
  "problem": {
    "domain": "acb",
    "precision_bits": 1694,
    "d": 3,
    "fb": -7,
    "w": 70,
    "d_lags": [],
    "denominators": [],
    "nhat_lags": [],
    "d0_inverse": null,
    "blocks": [],
    "assembly": null,
    "chop_digits": 25
  }
}
```

Indices in the SCC record are zero-based.  `structural_edges` uses the
source-to-target convention (matrix column to matrix row).  C++ recomputes
the SCC partition and condensation edges from that exact graph and rejects a
certificate that is merely plausible but not identical.

The frame base and width are part of the prepared-operator identity.  A
different work rectangle requires a different prepared chart handle.

`analytic.geometry` is optional for compatibility with recurrence-only
charts, but mandatory when a chart is referenced by a native composite SCC
object.  Composite preparation compares the complete canonical JSON record;
it never infers chart equality from Acb midpoint overlap or from a principal
operator hash.

`key` is only an index.  Reuse is allowed only if the complete static
operator JSON, `identity`, SCC certificate, session analytic identity, and
chart analytic identity are byte-identical after canonical protocol
assembly.  An unequal duplicate is a loud collision error.

### Run a retained operator

```json
{
  "schema": 2,
  "op": "chart.solve",
  "session": "s:1",
  "chart": "c:1",
  "run": {
    "nmax": 50,
    "p": 0,
    "has_initial": true,
    "adaptive_probe": false,
    "a_target": ["0", "0"],
    "b_target": ["0", "0"],
    "a_shift_min": 0,
    "a_shifts": [],
    "schedule": [],
    "initial": [],
    "initial_validity": [],
    "source": null,
    "return_u": false
  }
}
```

Every dynamic field is required.  This prevents accidental inheritance of a
previous column's seed, source, or resonance schedule.  The ordinary schema
1 result is returned with a `persistent` diagnostic containing preparation
time, run-parse time, run number, SCC size/depth, and
`static_tensor_copies: 0`.

Several independent columns or source sectors sharing the same prepared
operator may be submitted in one ordered batch:

```json
{
  "schema": 2,
  "op": "chart.solve_batch",
  "session": "s:1",
  "chart": "c:1",
  "output_digits": 520,
  "threads": 4,
  "runs": [{"nmax": 50}, {"nmax": 50}]
}
```

Each entry of `runs` has the complete dynamic shape of `chart.solve`'s
`run` object; the abbreviated objects above are illustrative only.  Results
remain in input order and each slot carries its own success or error record.
The worker count is bounded by the request size and the native limit.
Symbolic-rational batches are deliberately serialized until the shared FLINT
coefficient context has a separately verified parallel ownership model.

Independent runs spanning several already-retained charts in one session use
`session.solve_many`:

```json
{
  "schema": 2,
  "op": "session.solve_many",
  "session": "s:1",
  "threads": 10,
  "jobs": [
    {"chart": "c:1", "run": {}, "output_digits": 520},
    {"chart": "c:2", "run": {}, "output_digits": 520}
  ]
}
```

Every `run` is complete as above.  All chart handles are resolved under the
session lock before workers start, so an invalid job fails before partial
execution and a concurrent release cannot invalidate an accepted job.
Results preserve job order and carry per-slot typed errors.  This is the
native primitive for concurrent lower/upper-arm prewarming; it must not be
replaced by the old schema-1 batch, which reparses every static operator.

### Retain and evaluate a local solution

`local.solve` runs recurrence plus the retained assembly matrix and moves the
result directly into a session-owned typed `LocalSolution`.  It returns an
opaque `l:N` handle and summary metadata, never the coefficient slab.  Exact
chart geometry, exact `(a,b)` descriptors, normalized prescription records,
and checkpoint identity are mandatory and are bound to the recurrence target.

```json
{
  "schema": 2,
  "op": "local.solve",
  "session": "s:1",
  "chart": "c:1",
  "run": {"nmax": 50},
  "metadata": {
    "chart": {
      "center_exact": "0",
      "scale_exact": "1",
      "radius": "2",
      "infinite_radius": false
    },
    "tag": {
      "a": {"domain": "rational", "canonical": "1/2"},
      "b": {"domain": "rational", "canonical": "0"}
    },
    "prescriptions": [{
      "factor_exact": "t",
      "sign": 1,
      "multiplicity": 1,
      "leading_coefficient_sign": 1
    }],
    "checkpoint_identity": "exact-local-identity"
  }
}
```

The abbreviated `run` above is illustrative; it has the same complete
dynamic shape as `chart.solve`.  `local.evaluate` accepts an exact rational
real point and returns epsilon-framed value and theta-value vectors.  Its
optional explicit imaginary sign selects the requested rim; otherwise the
stored normalized prescriptions derive it.  Unresolved symbolic coefficient
fields are rejected rather than numerically sampled.

`local.certify_residual` evaluates the retained local internally and applies
the native Acb theta-residual certifier to a caller-supplied epsilon-framed
theta operator and optional inhomogeneous value at that point. The response is
`pass`, `fail`, or `inconclusive`; it is explicitly scoped to the stored Taylor
truncation unless a future certified tail majorant permits a full-local
conclusion. Operator, source, and checkpoint identities remain in the response
provenance. By default no solution or residual coefficient slab crosses the
bridge; the reported double bounds are diagnostics, while the verdict is
decided from native Arb magnitudes.

```json
{"schema":2,"op":"local.evaluate","session":"s:1","local":"l:1",
 "point":{"exact":"-1/2"},
 "options":{"imaginary_sign":-1,"t_order_reduction":0,
            "tail_estimate":false}}
{"schema":2,"op":"local.certify_residual","session":"s:1","local":"l:1",
 "point":{"exact":"-1/2"},
 "options":{"imaginary_sign":-1,"tail_estimate":false},
 "theta_operator":{"min":0,"max":0,"dimension":1,
                   "coefficients":["1/2"]},
 "source":null,"relative_tolerance":"1e-60",
 "scope":"stored_truncation","include_residual":false,
 "operator_identity":"theta-at-minus-half-v1",
 "checkpoint_identity":"residual-check-v1"}
{"schema":2,"op":"local.stats","session":"s:1","local":"l:1"}
{"schema":2,"op":"local.release","session":"s:1","local":"l:1"}
```

A retained local solution survives `chart.release`; `session.close` releases
both charts and locals.

#### Wolfram production seam

``DiffExp2`CppBackend` `` exposes the schema-2 lifecycle without exposing the
prepared-chart cache internals:

```wl
created = RunPersistentLocalSolve[request, persistentMetadata, localMetadata];
value = EvaluatePersistentLocal[created, <|"exact" -> "-1/2"|>,
  <|"imaginary_sign" -> -1, "tail_estimate" -> False|>, 80];
certificate = CertifyPersistentLocalResidual[created, <|
  "point" -> <|"exact" -> "-1/2"|>,
  "theta_operator" -> <|"min" -> 0, "max" -> 0, "dimension" -> 1,
    "coefficients" -> {"1/2"}|>,
  "relative_tolerance" -> "1e-60",
  "operator_identity" -> "theta-at-minus-half-v1",
  "checkpoint_identity" -> "residual-check-v1"|>];
stats = PersistentLocalStatistics[created];
ReleasePersistentLocal[created];
```

The solve wrapper goes through the same full-signature/session/chart
collision certificates as `RunPersistentRequest`; it does not construct a
second chart registry.  The three lifecycle wrappers accept either the raw
`local.solve` response or an association containing exact `Session` and
`Local` tokens. Residual certification additionally requires explicit operator
and checkpoint identities, and a source identity whenever a source is present.

``DiffExp2`Solve`SolveNativeLocalFamily[cs, req,
<|"a"->a,"b"->b,"p"->p|>, init]`` is the first narrow solver-level entry.
It builds the ordinary framed recurrence request, retains the assembled
result in C++, and returns an opaque handle/summary association with no
coefficient slab.  `p` remains the exact nonnegative integer `run.p`; `a`
and `b` are separately serialized as exact descriptors with an Acb
specialization (or as exact rationals in rational parity mode).

This entry is deliberately not a replacement for `SolveHomogeneous` yet. It
requires one non-SCC ChartSystem, the grouped spectral transform, an exact
identity gauge, a prepared frame with no pseudo-resonant compensation, and
a fully specialized numeric coefficient field.  In particular it rejects an
unresolved analytic regulator instead of sampling it.  Rank-reduced charts
must wait for a native *sequential* `V`-then-gauge assembly-chain operator:
composing the two matrices is not equivalent to the existing epsilon/Taylor
completeness accounting when exact or near cancellations occur.

### Prepare a typed composite SCC chart

`scc.prepare` is the schema-2 preparation boundary for native block-sequential
source propagation. It retains typed diagonal-chart pointers and typed
coupling multipliers; type erasure occurs only in the session registry.
General SCC execution remains disabled, but the deliberately narrow
`scc.solve_column` operation is enabled when stats certify
`exact-rational-regular-scalar-block-dag-column-v1` or the multidimensional
`exact-rational-regular-block-dag-column-v2`. The latter requires a regular
singleton recurrence partition inside every retained diagonal block. The
explicit Wolfram `SolveNativeSCCBasisColumn` migration seam accepts both
scopes, validates the corresponding statistics/capability string, and fails
loudly for any other retained composite.

The parent manifest carries two machine-checkable `D` by `D` row-major exact
matrix records. Every cell has the shape
`{"exact": string, "proven_zero": bool}`. Original and theta zero facts must
agree and reproduce the supplied source-to-target structural graph exactly.
No zero decision is made from an Acb coefficient slab.

```json
{
  "schema": 2, "op": "scc.prepare", "session": "s:1",
  "key": "collision-certified composite key",
  "identity": "exact parent identity",
  "parent": {
    "dimension": 2,
    "exact_system_record": [
      [{"exact":"0","proven_zero":true},
       {"exact":"0","proven_zero":true}],
      [{"exact":"g","proven_zero":false},
       {"exact":"0","proven_zero":true}]],
    "exact_theta_record": [
      [{"exact":"0","proven_zero":true},
       {"exact":"0","proven_zero":true}],
      [{"exact":"theta-g","proven_zero":false},
       {"exact":"0","proven_zero":true}]],
    "chart": {
      "center_exact": "0", "scale_exact": "1",
      "radius_exact": "2", "infinite_radius": false,
      "prescriptions": []
    },
    "scc": {
      "components": [[0],[1]], "structural_edges": [[0,1]],
      "condensation_edges": [[0,1]], "topological_order": [0,1],
      "coupling_depth": 1
    },
    "execution": {"mode":"BlockSequentialStrict", "work_t_order":50},
    "work_contract": {
      "work_min": -4, "requested_min": -2, "requested_max": 6,
      "work_complete_max": 20, "public_t_order": 44,
      "wolfram_coupling_depth": 2
    }
  },
  "blocks": [
    {"block": 0, "vertices": [0], "chart": "c:1",
     "principal_identity": "exact diagonal identity 0",
     "regular": true, "identity_gauge": true,
     "identity_v": true, "no_pseudo": true},
    {"block": 1, "vertices": [1], "chart": "c:2",
     "principal_identity": "exact diagonal identity 1",
     "regular": true, "identity_gauge": true,
     "identity_v": true, "no_pseudo": true}],
  "couplings": [{
    "source_block": 0, "target_block": 1,
    "source_vertices": [0], "target_vertices": [1],
    "rows": 1, "columns": 1,
    "exact_identity": "exact coupling matrix identity",
    "domain": "rational", "symbols": [],
    "entries": [{
      "row": 0, "column": 0,
      "source_vertex": 0, "target_vertex": 1,
      "exact_original_entry": "g", "exact_theta_entry": "theta-g",
      "multiplier": {
        "epsilon_shift": -1, "center_pole_order": 0,
        "kernels": [], "exact_identity": "theta-g",
        "proven_zero": false
      }
    }]
  }]
}
```

The matrices and kernels above are abbreviated: both parent matrices are
square of size `dimension`, active multiplier kernels cover the complete
epsilon/Taylor work rectangle, and the SCC certificate has its ordinary full
schema. Blocks, coupling groups, and entries are in deterministic
block/edge/row-column order. Each entry's global and local indices, exact
parent cells, theta identity, coefficient domain, and structural-zero fact
are checked one-to-one.
The coupling group's caller `exact_identity` remains part of the collision
signature only. C++ derives the authoritative retained matrix identity from
the validated source/target bases and indexed parent original/theta cells.

Composite-capable `chart.prepare` requests additionally bind these optional
analytic records into the ordinary chart collision signature:

- `analytic.geometry`, byte-equal to the normalized parent chart record;
- `analytic.principal_matrix`, the indexed local original-matrix subblock;
- `analytic.native_scc_capabilities` with `regular`, `identity_gauge`,
  `identity_v`, and `no_pseudo` facts.

The first slice requires all four capabilities true. Only `identity_v` is
native-proved here, against the retained identity assembly operator;
`regular`, `identity_gauge`, and `no_pseudo` remain collision-bound producer
certificates which execution revalidates or refuses.  Stats label that
distinction explicitly. All diagonal charts must share the scalar domain,
geometry, work frame, and exact SCC principal graph. This first slice accepts
only exact rational center/scale and finite radius (the scale is nonzero and
radius strictly positive); unsupported algebraic geometry is loud. Signed
epsilon shifts are retained, and execution proves that `work_min` supplies
enough lower halo.

The first Wolfram producer is
``Solve`PrepareNativeSCCComposite[sccEnvelope, req]``.  It captures each
diagonal block through the existing grouped homogeneous request builder under
isolated Solve caches and stops before execution.  The producer accepts only
regular, collision-free blocks with exact identity gauge and spectral
transforms.  Every captured group must share one scalar field and identical
`fb`, `w`, and guarded `nmax`; no alternate preparation is used to force a
common rectangle.  It enriches each block's ordinary `chart.prepare`
analytic record with `geometry`, the nested exact `principal_matrix`, and
`native_scc_capabilities`, then supplies the handle-free parent/block/coupling
manifest to ``CppBackend`PreparePersistentSCC``.  The bridge alone injects
opaque chart handles and exact prepared-chart identities.

This explicit API always describes `BlockSequentialStrict`; it neither
consults nor executes the older Wolfram monolithic-coarsening performance
heuristic, and no production solve dispatch uses the returned handle yet.
Its bounded Solve-level cache compares the full parent chart/request/config
signature after a hash lookup, validates a hit through `scc.stats`, rebuilds a
directly released stale handle, and refuses capacity exhaustion rather than
evicting a live public SCC object.  `ClearSolveCaches[]` clears that cache and
closes the native sessions.  Its compact private execution descriptor is
accepted only within that full-signature cache entry and binds the same public
handle.  It retains run-only block records, block dimensions, graph,
field/frame contract, identity, precision, and last checked stats; large
captured static operators and block systems are released after preparation.
The public preparation result remains opaque.

``Solve`SolveNativeSCCBasisColumn[sccEnvelope, req, seedBlock,
seedLocalComponent:1]`` is the first consumer.  It requires one certified
canonical rational homogeneous run per local block component (`p=0`, `a=b=0`,
one exact eps^0 unit seed, fixed frame, and one exact regular `R/T` schedule
step per retained singleton).  Wolfram derives each run's component from the
flattened initial tensor, proves that every block supplies a complete
permutation of its local unit columns, and sends the selected seed's captured
dynamic run unchanged.  For every reachable descendant in the certified
topological order it derives a correctly dimensioned exact-zero/no-initial run
from that descendant's captured component-one run.  C++ multiplies and
combines retained predecessor locals, injects the vector source directly into
the target recurrence, embeds each block result into the parent dimension,
caps to the public window, and retains the final local.  The response contains
an opaque local handle and exact SCC/seed/global-basis provenance with
`json_coefficients: 0`; for v2 Wolfram also verifies the encoded local
component in the exact column identity.  There is no coefficient slab
roundtrip.  Parent chart geometry and ordered analytic prescriptions are
rechecked on every local metadata record.  The three-argument scalar-v1 call
keeps its original wire request and result shape.  This explicit column API is
not a `SolveHomogeneous` or transport dispatch path.

```json
{
  "schema":2, "op":"scc.solve_column", "session":"s:1", "scc":"scc:1",
  "checkpoint_identity":"exact column checkpoint",
  "seed":{"block":0,"run":{"...":"canonical captured run"},
          "metadata":{"...":"exact parent local metadata"}},
  "targets":[
    {"block":1,"run":{"...":"derived zero-initial run"},
     "metadata":{"...":"exact parent local metadata"}}
  ]
}
```

```json
{"schema":2,"op":"scc.stats","session":"s:1","scc":"scc:1"}
{"schema":2,"op":"scc.release","session":"s:1","scc":"scc:1"}
```

### Inspect and destroy

```json
{"schema":2,"op":"session.stats","session":"s:1"}
{"schema":2,"op":"chart.release","session":"s:1","chart":"c:1"}
{"schema":2,"op":"session.close","session":"s:1"}
```

LibraryLink unload also destroys every retained chart before its FLINT
coefficient context is torn down.

## Analytic regularization contract

The native recurrence never guesses a `+i delta` or `-i delta` prescription
from an Acb enclosure.  Wolfram continues to make the exact structural and
branch decision.  The resulting prescription identity is retained with the
chart and participates in cache collision certification.  Exact regulator
coefficients may use the symbolic rational domain; prescribed numerical
coefficients may use Acb.  A later native matching/evaluation milestone must
consume the same stored analytic identity rather than re-derive it on the
real axis. Native evaluation and stored-truncation residual certification now
do so through the retained local; matching must preserve the same rule.

## Wolfram integration sequence

1. Split `cppRunRecursionCore`'s existing schema 1 association into a static
   operator record and a dynamic run record.
2. Create one session per loaded differential system/configuration and close
   it from `ClearSolveCaches`/`LoadSystem`.
3. Cache chart handles on the existing collision-certified chart/frame key.
4. Submit homogeneous columns and particular sectors as `chart.solve`
   commands.  Upper and lower arms may share handles and run concurrently.
5. Move SCC coupling multiplication/source propagation into the session so
   the large source tensor no longer returns to Wolfram between blocks.

The next native milestones extend the retained session rather than adding
new stateless calls: matching/refinement, promotion from stored-truncation to
certified full-local residuals once rigorous tail majorants exist, endpoint
limits and line/tile integration, and finally native checkpoint state.
