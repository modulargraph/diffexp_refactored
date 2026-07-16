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

## Pentagon gate

The massless-pentagon timing regression is also opt-in. It exercises the
algebraic-chart transport that was recovered in July 2026 and verifies the
finite value as well as the elapsed time:

```bash
FT_FIRE_PATH=/path/to/FIRE7 \
FT_PREP_CACHE_DIR=/path/to/prepared-cache \
Scripts/run_pentagon_timing_regression.sh
```

Or append it to the ordinary release suite with
`DE2_RUN_PENTAGON_TIMING=1`. On a new machine, set
`PENTAGON_TIMING_WARM_CACHE=1` once; as for the banana gate, FIRE preparation
is outside the measured region and the timed run must report a cache hit.

The pinned configuration is working precision 300, matching contract 8
decimal digits, Taylor order 25, boundary extra order 16, division order 3,
and 10 C++ threads by default. The runner starts with fresh transport
checkpoints and supplies the exact private level-1 matching halo of two orders;
this leaves the public epsilon request unchanged while avoiding a discovery
replay. It checks the finite value against
`-0.025411779885306218280920810082412842534323133089462` with tolerance
`1e-12` and fails above 300 seconds. Override the hardware qualification
ceiling explicitly with `PENTAGON_TIMING_MAX_SECONDS`.

The first integrated recovery run took 74.54 seconds while restoring the
level-1 boundary; profiling showed that rebuilding the short upstream ladder
adds roughly 5--6 seconds. The five-minute gate is intentionally conservative
enough for normal machine variability while still catching a return to the
former repeated-owner/retry behavior.

## Planar double-box gate

The planar double-box timing regression is an opt-in end-to-end gate for the
multi-level rational-chart transport and its complete Laurent result:

```bash
FT_FIRE_PATH=/path/to/FIRE7 \
FT_PREP_CACHE_DIR=/path/to/prepared-cache \
Scripts/run_double_box_timing_regression.sh
```

Or append it to the ordinary release suite with
`DE2_RUN_DOUBLE_BOX_TIMING=1`. On a new machine, set
`DOUBLE_BOX_TIMING_WARM_CACHE=1` once. FIRE preparation remains outside the
timer, and the measured run must report a preparation-cache hit.

The pinned configuration is working precision 150, matching contract 8
decimal digits, Taylor order 25, boundary extra order 16, division order 3,
epsilon order zero, and 10 C++ threads by default. The run starts with fresh
transport checkpoints and calls the internal example runner with the learned
private matching profile `<|1 -> 2, 2 -> 5, 3 -> 1|>`; ambient learned-profile
loading is disabled, the public epsilon request stays at zero, and any
matching-halo discovery retry fails the gate.

After transport, `Scripts/verify_double_box_planar.m` checks every coefficient
from epsilon^-4 through epsilon^0 at 8 digits. The default wall-time ceiling is
1500 seconds. The July 2026 recovery replay took about 1189 seconds after the
final matching profile was known; this initial ceiling leaves normal machine
variance while still catching a return to discovery replay. Override it
explicitly with `DOUBLE_BOX_TIMING_MAX_SECONDS` when
qualifying different hardware and record that override with the release
results.
