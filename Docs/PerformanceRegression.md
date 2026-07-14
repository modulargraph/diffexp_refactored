# Performance regression gates

The full banana timing regression is an opt-in release gate. It is not part
of the ordinary test suite because it occupies a Wolfram license and takes
several minutes.

Run it directly:

```bash
FT_FIRE_PATH=/path/to/FIRE7 \
FT_PREP_CACHE_DIR=/path/to/prepared-cache \
Scripts/run_banana_timing_regression.sh
```

Or include it after the normal release suite:

```bash
DE2_RUN_BANANA_TIMING=1 \
FT_FIRE_PATH=/path/to/FIRE7 \
FT_PREP_CACHE_DIR=/path/to/prepared-cache \
Scripts/run_release_tests.sh
```

On a machine without a prepared banana cache, set
`BANANA_TIMING_WARM_CACHE=1` for the first run. FIRE preparation happens
before the timer. The measured run must report `FTPREP CACHE HIT`; this keeps
external FIRE performance from hiding a DiffExp2 transport regression.

The pinned production configuration is:

- working precision 500;
- matching contract 20 decimal digits;
- Taylor expansion order 50;
- division order 3;
- the C++ recurrence backend with 10 threads by default.

The gate exercises the ordinary production runner, including automatic
private epsilon-reservoir retries. It checks the final banana value against
`8.2681045358689687315430153454799888687` with tolerance `1e-20` and fails
if wall time exceeds 300 seconds. The ceiling can be overridden explicitly
with `BANANA_TIMING_MAX_SECONDS` when qualifying different hardware; such an
override should be recorded with the release results.

Before the full timer starts, the harness also solves both singular level-1
endpoint bases at order 50 and requires a combined time below 60 seconds.
This diagnostic sub-gate makes endpoint/Fuchsian recurrence regressions
immediately visible. Set `BANANA_SINGULAR_MAX_SECONDS` to qualify different
hardware, or `BANANA_TIMING_SKIP_SINGULAR_GATE=1` when only reproducing the
historical end-to-end number.

The July 2026 recovery baseline on the development Mac was 264.5 seconds for
the ordinary production path and 162.2 seconds with the already-known private
epsilon halos. The production-path number is the regression baseline because
it includes the bounded discovery retries users actually encounter.
