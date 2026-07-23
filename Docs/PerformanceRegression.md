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
Every measured command runs in its own process group under a hard wall-clock
deadline. Exceeding the configured ceiling terminates the Wolfram process and
its children with status 124, rather than waiting indefinitely and merely
reporting a late result after it eventually finishes.

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

## Four-loop banana gates

The four-loop banana has an opt-in boundary gate that crosses both of its first
two singular handoffs and stops after the level-2 boundary. It is a
production-chain correctness regression:
manifest forwarding, exact Rational-shadow retention, singular physical-ODE
construction, real-ray terminal composed-adjoint certification, and the
near-endpoint rational-pole case are all exercised together.

```bash
FT_FIRE_PATH=/path/to/FIRE7 \
FT_PREP_CACHE_DIR=/path/to/prepared-cache \
Scripts/run_banana4_boundary_regression.sh
```

Or append it to the release suite with `DE2_RUN_BANANA4_BOUNDARY=1`. The
pinned configuration is working precision 500, matching contract 25 digits,
Taylor order 50, public epsilon order 4, private matching halos
`20,14,8,0`, private composed-adjoint order 100, and 10 C++ threads by
default. A fresh transport checkpoint is mandatory. The gate requires a
center-ending composed-adjoint tail which also satisfies its downstream
accuracy budget in every executed ladder level. It rejects epsilon-reservoir,
stalled/exhausted Taylor, and unsupported markers; a bounded Taylor retry is
allowed only while its coefficient verdict strictly improves. A noncenter
terminal tile is accepted only when it is explicitly classified as outside
the one-contraction theorem.

For a quicker diagnostic that only produces the level-3 boundary, set
`BANANA4_BOUNDARY_TARGET_LEVEL=3`. This is not the release default because it
does not exercise the near-endpoint pole that previously failed at the next
handoff. `BANANA4_BOUNDARY_ADJOINT_ORDER` may be raised for experiments but
may not be smaller than `BANANA4_BOUNDARY_EXPANSION_ORDER`.

The default ceiling is 420 seconds and can be changed with
`BANANA4_BOUNDARY_MAX_SECONDS` when qualifying hardware. As with the other
gates, the measured run must be a preparation-cache hit. Set
`BANANA4_BOUNDARY_WARM_CACHE=1` once when populating a new cache; preparation
is not charged to the timer.

The stronger release gate executes level 1, requires final publication, and
compares the finite value with the independent five-propagator Bessel oracle
to eight decimal digits:

```bash
BANANA4_BOUNDARY_TARGET_LEVEL=final \
Scripts/run_banana4_boundary_regression.sh
```

Or append it to the release suite with `DE2_RUN_BANANA4_FULL=1`. Full mode
uses a 15-digit internal matching contract and a 2400-second ceiling by
default, leaving seven guard digits beyond the public eight-digit oracle while
covering the regular near-endpoint overlaps and final pole audit that the
level-2 boundary gate cannot observe.

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
level-1 boundary.  After the frame-independent regular-owner transport was
wired through production, a warm-preparation run with fresh transport
checkpoints completed in 173 seconds on July 16, 2026.  Its finite coefficient
differed from the pinned oracle by about `2.10e-14`, and it entered no matching
discovery retry.  The five-minute gate is intentionally conservative enough
for normal machine variability while still catching a return to the former
repeated-owner/retry behavior.

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
