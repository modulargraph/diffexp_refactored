# Persistent FIRE Storage

> **Status: NO-GO (2026-07-11).** Real FIRE 6.5.2 cold/warm parity fails.
> FIRE's `#storage` databases are not a safe extensible cache across separate
> invocations with new integral requests. Do not enable or integrate this
> prototype. The code and benchmark remain only as negative research evidence.

Successive `CloseLevelMasterBasis` iterations often ask FIRE for derivatives
of masters discovered by the preceding reduction.  The exact integral cache
cannot predict those masters, while FIRE's default temporary database is
discarded after every invocation.  Rebuilding that database dominated the
observed unequal-mass banana preparation time.

The prototype added an opt-in configuration switch intended to reuse FIRE's
internal sector database across those successive requests:

```mathematica
FeynmanTrick`SetFTOption["PersistentFIREStorage", True];
FeynmanTrick`SetFTOption["FIREStorageDirectory", "/fast/local/cache"];
```

The research runner exposes equivalent environment settings, but they must not
be used for production calculations:

```sh
FT_PERSISTENT_FIRE_STORAGE=1 \
FT_FIRE_STORAGE_DIR=/fast/local/cache \
wolframscript -file Scripts/run_ft_stepwise2.m
```

When the directory is `Automatic` (the default), storage lives below the
configured FeynmanTrick work directory.  Persistence is disabled by default.

## Safety model

The cache key is a SHA-256 digest over the exact `.start` file, hashes of the
FIRE/FLAME/Fermat binaries and relevant libraries, variable order, problem and
topology metadata, FIRE thread/restriction options, and platform identifiers.
Changing any of those inputs selects a new database rather than risking an
incompatible reuse.

Each content key has one atomic `CURRENT` pointer to an immutable generation.
The committed generation is never opened by FIRE.  A reduction instead:

1. acquires a per-key directory lock;
2. clones `CURRENT` into an attempt-local directory;
3. runs FIRE with `#storage <attempt>` and `#keepall`;
4. requires a zero FIRE exit code and successfully parsed rules and masters;
5. moves the attempt into a new immutable generation and atomically replaces
   `CURRENT`.

Failed executions, retries, missing tables, parse failures, and interrupted
publication discard the attempt and preserve `CURRENT`.  Retries are reset
from the committed generation, so a partially mutated database is never used
as their starting point.  Concurrent writers for the same key fail loudly.
This lock does not broaden the existing FIRE work-directory concurrency
contract: independent processes must still use independent work directories.
An uncatchable process kill or power loss can leave `LOCK` behind deliberately;
inspect `LOCK/owner.json` and remove that lock only after confirming the owner
process is gone.  The last `CURRENT` generation remains untouched.

The `!` modifier is intentionally not used.  FIRE copies reusable forward
reduction databases to ordinary `#storage` before substitution.  With `!`, it
overwrites those files after substitution with state specialized to the current
request set; a later disjoint request can then inherit invalid virtual points.
The real-FIRE parity benchmark covers this distinction.

The returned FIRE master list is also inserted into the in-memory exact
reduction cache as identity reductions.  This avoids needlessly asking FIRE
for a master it just reported, while the fixed-point basis-closure loop remains
in place for correctness.

## Validation

`Tests/test_fire_persistent_storage.m` uses a fake FIRE executable and stubbed
table readers.  It checks generation cloning, key invalidation, writer
serialization, successful publication, master-identity caching, and rollback
after both FIRE-exit and table-parse failures without starting FIRE itself.

The repository includes a bounded real-FIRE microbenchmark for this purpose.
It uses two disjoint request batches for a fixed Euclidean one-loop bubble and
disables the in-memory reduction cache:

```sh
FT_FIRE_PATH=/path/to/FIRE6 \
wolframscript -file Scripts/bench_fire_persistent_storage.m
```

Set `FT_FIRE_STORAGE_BENCH_KEEP=1` to preserve its logs and generated storage.
The benchmark intentionally exits nonzero when it reproduces the no-go.

### Real FIRE 6.5.2 result

The fixed Euclidean one-loop bubble used disjoint closure-like batches:

```text
round 1 = {{1,1}, {2,1}, {1,2}}
round 2 = {{3,1}, {2,2}, {1,3}}
```

Persistence-off and persistence-on cold round 1 agreed exactly, with master
set `{{0,1},{1,0},{1,1}}`. The warm round loaded four database files (11 copy
events total) and FIRE reported every sector as `nothing to do`, proving that
storage was consumed. Nevertheless, its table contained invalid virtual point
identifiers such as `80-10-1`. The strict unresolved-master guard rejected the
table, discarded the attempt, and left `CURRENT` on the cold generation.

| configuration | off R1 | off R2 | on cold R1 | on warm R2 | result |
|---|---:|---:|---:|---:|---|
| wall seconds, plain storage + keepall | 2.098626 | 2.093978 | 2.663787 | 2.168900 | warm rejected |
| FIRE seconds, plain storage + keepall | 0.134010 | 0.125736 | 0.196977 | 0.182933 | warm rejected |
| FIRE seconds, cumulative R1 union R2 | 0.132988 | 0.121922 | 0.184110 | 0.179930 | same invalid virtual points |
| FIRE seconds, plain storage without keepall | 0.134768 | 0.127994 | 0.205186 | 0.200028 | silently wrong |

Removing `#keepall` is worse: warm round 2 completes and publishes, but returns
the three requested integrals as identities and reports them as masters:

```text
G[1,{3,1}], G[1,{2,2}], G[1,{1,3}]
```

The correct persistence-off reductions instead use the master set
`{{0,1},{1,0},{1,1}}`. Explicitly supplying `round 1 union round 2` on the warm
invocation does not repair the virtual points. The `!` storage modifier also
fails. Thus neither request accumulation nor storage-mode selection makes the
internal databases safe across FIRE processes, and there is no measured speed
benefit even in the failing cases.

The exact correct and invalid warm FIRE tables are frozen in
`Tests/refs/fire_persistent_storage/`. The recommended optimization direction
is one-invocation closure batching, not cross-invocation database reuse.
