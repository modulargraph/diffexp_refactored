# Persistent FIRE Storage

Successive `CloseLevelMasterBasis` iterations often ask FIRE for derivatives
of masters discovered by the preceding reduction.  The exact integral cache
cannot predict those masters, while FIRE's default temporary database is
discarded after every invocation.  Rebuilding that database dominated the
observed unequal-mass banana preparation time.

`PersistentFIREStorage` is an experimental, opt-in configuration switch that
reuses FIRE's internal sector database across those successive requests:

```mathematica
FeynmanTrick`SetFTOption["PersistentFIREStorage", True];
FeynmanTrick`SetFTOption["FIREStorageDirectory", "/fast/local/cache"];
```

For `Scripts/run_ft_stepwise2.m`, the equivalent environment settings are:

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

Before making this mode the default, run an actual FIRE parity/performance
trial from a cold cache, then repeat from the populated cache.  The reductions
and master set must agree exactly with persistence disabled.

The repository includes a bounded real-FIRE microbenchmark for this purpose.
It uses two disjoint request batches for a fixed Euclidean one-loop bubble and
disables the in-memory reduction cache:

```sh
FT_FIRE_PATH=/path/to/FIRE6 \
wolframscript -file Scripts/bench_fire_persistent_storage.m
```

Set `FT_FIRE_STORAGE_BENCH_KEEP=1` to preserve its logs and generated storage.
