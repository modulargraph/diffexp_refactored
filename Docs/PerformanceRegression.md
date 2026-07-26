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

The full four-loop investigation is intentionally paused as of 2026-07-26
while other families are developed. No final banana4 result has passed the
independent Bessel oracle. The boundary gate remains useful as an opt-in
regression, but the full gate is diagnostic rather than release evidence until
the correlated clearance-seed terminal path completes in production. See
[`Banana4Status.md`](Banana4Status.md) for the precise stopping point, durable
artifacts, closed experiments, and resume protocol.

The four-loop banana has an opt-in boundary gate that crosses both of its first
two singular handoffs and stops after the level-2 boundary. It is a
production-chain correctness regression:
manifest forwarding, exact Rational-shadow retention, singular physical-ODE
construction, real-ray terminal composed-adjoint certification, and the
near-endpoint rational-pole case are all exercised together.
The ordinary native suite complements this slow gate with two deterministic
producer/consumer contract tests: one distinguishes an inconclusive propagated
enclosure from a genuine midpoint residual (so Taylor order is not retried
blindly), and one carries unequal rigorous specializations of the same exact
algebraic chart radius through matching and checkpoint replay. These fixtures
are intentionally end-to-end at the native ownership boundary; local SCC
solver tests alone cannot detect either failure.

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
the one-contraction theorem. Native stage records are enabled as part of the
gate: singular-chart basis crossings remain allowed, but any basis fallback
whose receiver is regular fails the gate. This protects the private epsilon
reservoir and prevents one under-certified value hop from silently becoming
a long fallback cascade.

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
currently pins the B4-130 diagnostic producer profile
`<|2->19,3->25,4->25|>` and a 2400-second ceiling. B4-130 is a recorded
failure and must not be used as release evidence; the next focused replay
must follow the B4-131 C++ value/tail audit before this full-gate profile is
advanced. B4-131 disproved internal-certification clipping as the cause of
the regular-hop regression, so that experimental policy change was removed.
The published result will still require the independent eight-digit oracle
plus the regular near-endpoint and final-pole audits.

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
checkpoints and supplies the learned private matching halos
`<|1 -> 3, 2 -> 1|>`. Producer accuracy is pinned separately at levels 2
through 4 to 12, 14, and 16 digits. These are private transport resources:
the public epsilon request remains order zero and the gate performs no
discovery replay. It checks the finite value against
`-0.025411779885306218280920810082412842534323133089462` with tolerance
`5e-10`, appropriate to the deliberately low eight-digit matching profile and
compact FINAL record, and fails above 300 seconds. Override the hardware qualification
ceiling explicitly with `PENTAGON_TIMING_MAX_SECONDS`.

The first integrated recovery run took 74.54 seconds while restoring the
level-1 boundary.  After the frame-independent regular-owner transport was
wired through production, a warm-preparation run with fresh transport
checkpoints completed in 173 seconds on July 16, 2026.  Its finite coefficient
differed from the pinned oracle by about `2.10e-14`, and it entered no matching
discovery retry.  The five-minute gate is intentionally conservative enough
for normal machine variability while still catching a return to the former
repeated-owner/retry behavior.

The July 23 recovery added source-specific enclosure diagnostics and a
certified midpoint-candidate retry for regular matching weights. A fresh
transport run then completed in 92.67 seconds with finite value
`-0.0254117799`. The midpoint solve is proposal-only: the untouched Acb basis
and incoming boundary must still pass the full-ball forward residual before
the candidate is published.

## Original Henn canonical pentagon gate

The independent 108-master canonical example from the original DiffExp
repository is opt-in:

```bash
Scripts/fetch_henn_nonplanar_data.sh
Scripts/run_henn_nonplanar_timing_regression.sh
```

It can be appended to the release suite with
`DE2_RUN_HENN_NONPLANAR_TIMING=1`. The gate checks all published endpoint
coefficients through epsilon order 4 with maximum absolute tolerance `1e-9`.
Its default settings are working precision 50, expansion order 25, and Padé
chart matching.

The default hard deadline is 120 seconds. The transport-only ceiling is 30
seconds, compared with the original notebook's saved 49.765221-second
order-25 transport. Set `HENN_NONPLANAR_DEADLINE_SECONDS` or
`HENN_NONPLANAR_MAX_TRANSPORT_SECONDS` only for an explicit hardware
qualification. Ancillary data are hash-pinned and must be fetched before the
gate; ordinary tests never access the network.

## Original planar one-mass canonical pentagon gate

The four canonical systems from the original DiffExp
`5pPlanar1Mass.nb` example are opt-in:

```bash
Scripts/fetch_planar_one_mass_data.sh
Scripts/run_planar_one_mass_timing_regression.sh
```

They can be appended to the release suite with
`DE2_RUN_PLANAR_ONE_MASS_TIMING=1`. The gate checks the PH1-to-PH6 endpoint
for every master and every coefficient from epsilon order 0 through 4, with
maximum absolute tolerance `1e-8`. Default numerical settings are working
precision 50, expansion order 25, no Padé matching, and clearance-certified
algebraic charts.

The hard deadline is 900 seconds. Default comparable-time ceilings are 60
seconds for `1loop`, 300 seconds for `zmz` and `mzz`, and 350 seconds for
`zzz`. They can be explicitly qualified on other hardware with
`PLANAR_ONE_MASS_MAX_1LOOP_SECONDS`,
`PLANAR_ONE_MASS_MAX_ZMZ_SECONDS`,
`PLANAR_ONE_MASS_MAX_MZZ_SECONDS`, and
`PLANAR_ONE_MASS_MAX_ZZZ_SECONDS`. Ancillary data are hash-pinned and
ordinary tests never access the network.

## Original multiple-polylogarithm gate

The seven evaluations from the original DiffExp
`MultiplePolylogarithms.nb` example are opt-in:

```bash
Scripts/run_mpl_timing_regression.sh
```

They can be appended to the release suite with `DE2_RUN_MPL_TIMING=1`.
The gate checks the three weight-3/4 examples at expansion orders 50 and 75,
then the weight-20 example at expansion order 100. The default hard deadline
is 180 seconds and the weight-20 ceiling is 120 seconds. Override them with
`MPL_DEADLINE_SECONDS` and `MPL_MAX_WEIGHT20_SECONDS`.

The short examples require errors below `1e-23` or `1e-34`; the weight-20
stress case requires error below `1e-40`. No external data or network access
is needed.

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
