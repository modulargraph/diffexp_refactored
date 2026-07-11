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
  "analytic": {"delta_sign": -1, "branch_signature": "..."},
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

```json
{"schema":2,"op":"local.evaluate","session":"s:1","local":"l:1",
 "point":{"exact":"-1/2"},
 "options":{"imaginary_sign":-1,"t_order_reduction":0,
            "tail_estimate":false}}
{"schema":2,"op":"local.stats","session":"s:1","local":"l:1"}
{"schema":2,"op":"local.release","session":"s:1","local":"l:1"}
```

A retained local solution survives `chart.release`; `session.close` releases
both charts and locals.

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
real axis.

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
new stateless calls: matching/refinement and residual certification, then
endpoint limits and line/tile integration, and finally native checkpoint
state.
