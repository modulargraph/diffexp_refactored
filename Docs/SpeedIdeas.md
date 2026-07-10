# SpeedIdeas: measured per-chart profile and ranked levers

## 2026-07-10 algebraic regular-chart numeric handoff

Classic balanced segmentation can produce exact quadratic-algebraic regular
centers and scales.  Small transformed operator coefficients stayed exact
under the ByteCount gate, so a scalar EO50 level-3 recurrence accumulated
rational/algebraic expression swell for more than two minutes.  Preparation
itself took only 0.006 s; the cost was entirely inside `runRecursion`.

For regular charts whose center or scale is an exact non-rational numeric
quantity, `prepareCleared` now grounds every exact numeric operator
coefficient once at 2x WorkingPrecision.  Planning geometry and structural
data remain exact, as do rational charts and symbolic analytic regulators.
The predicate must be `NumericQ`, not `NumberQ`: Wolfram returns `False` from
`NumberQ` for expressions such as `a + b Sqrt[n]` even though they are exact
numeric quantities.

On the exact fresh banana level-3 upper chart at WP1000/EO50, total
`SolveChart` time fell from more than 120 s to 0.090 s: recurrence 0.038 s,
assembly 0.003 s, and residual certification 0.027 s.  The focused grouped,
solve, and transport suites pass 17/17, 40/40, and 62/62 respectively.

The end-to-end gate exposed the same algebraic cancellation again at the
integration seam.  Numerical grounding now consistently uses
`RootReduce -> N[2 WP] -> SetPrecision[2 WP] -> Chop[ChopFloor]` inside a
raised `$MaxExtraPrecision` block, but only after exact planning and pole
classification.  This policy applies to prepared solve coefficients,
rational-multiplier Taylor coefficients, and final local integration bounds;
symbolic analytic regulators remain exact.  `ESTrim` also uses a unit
decision scale only for an all-centered-zero frame, preventing resolved
`0``A` rows from manufacturing dozens of fake Laurent pole orders.

At EO50/division-order 3, the level-2 boundary data is truncation-limited to
roughly 25 digits regardless of WP1000, so the production validation run uses
WP500; the larger precision paid bignum cost without improving that boundary.
After the remaining ladder integrations, the completed banana finite part
`8.268104535868963836...` agrees with the existing oracle for roughly 14
digits and retains residual epsilon poles of order `10^-16`.  Raising WP
alone cannot recover the final digits lost to the fixed Taylor order.

The WP500 gate also exposed a distinct singular-endpoint window failure:
the adaptive recurrence correctly found physical basis starts at epsilon
orders -2/-1, but algebraic cancellations stored as high-accuracy centered
zeros (`0``A`) made the structural scan retain scratch starts at -52/-51.
The subsequent honest Cauchy product then had no overlapping complete
window and returned only zeros.  Framed solve assembly now applies an
uncertainty-certified absolute chop before its first-nonzero scan: `Chop`
is only a candidate filter, and the full magnitude upper bound must be below
`ChopFloor`.  Thus unresolved low-accuracy zeros remain material and
symbolic regulators remain exact.  On the exact level-2 upper banana replay,
the matched endpoint window is restored from `-52..-39` to `-2..11` without
changing the incoming values or matching weights.  The resumed WP500/EO50
ladder then completed through level 0; its level-1 lower and upper singular
endpoint solves took about 163 s and 626 s respectively, establishing the
upper endpoint recurrence/assembly as the dominant remaining hotspot.

## 2026-07-10 adaptive singular lower epsilon frames

The scalar `recurrencePoleDepth` bound remains the authoritative upper
completeness halo, but singular homogeneous solves no longer have to mirror
that pessimistic depth below the physical epsilon window.  They first run the
same strict recurrence with the deepest exact single-use pole across all
prepared Taylor lags below the window.
Negative matrix shifts apply their complete coefficient matrix before testing
underflow, so rank-one/square-zero annihilation is visible.  Genuine content
at the lower boundary triggers monotone geometric widening; the terminal
retry is exactly the former scalar rectangle.  There is no alternate solver,
tolerance trimming, or assumption that rank-one maps are nilpotent.

Focused controls cover scalar repeated poles, rank-one idempotents,
square-zero maps, a deep high-Taylor-lag pole hidden behind a pole-free first
lag, symbolic analytic regulators, and a rational upper tail that later
crosses into the delivered window.  The latter is why only the
lower frame is adaptive in this first production slice: the proven full upper
halo and `UValid`/`CompleteMax` propagation are unchanged.

Measured on the real banana L1 endpoint charts at WP100, epsilon order 2;
every adaptive result was coefficient-identical to a forced terminal-frame
run:

| endpoint/order | adaptive width | terminal width | time adaptive | time terminal | speedup |
|---|---:|---:|---:|---:|---:|
| x=0, N=10 | 29 | 47 | 1.025 s | 1.218 s | 1.19x |
| x=0, N=20 | 49 | 87 | 5.236 s | 6.459 s | 1.23x |
| x=1, N=10 | 29 | 47 | 2.491 s | 2.812 s | 1.13x |
| x=1, N=20 | 49 | 87 | 14.986 s | 15.077 s | 1.01x |

All seven banana columns completed on their first adaptive rectangle.  The
width ratio approaches two as Taylor order grows because this slice removes
the false lower `2 N` depth while deliberately retaining the upper `2 N`
certificate.  At x=1 the grouped `Q=1+5 eps` causal operator work currently
dominates and scales less strongly with this lower-only reduction.  This is a
conservative first-stage gain, not an order-100 measured speedup claim.

The reproducible harness is `Scripts/bench_adaptive_lower_frame.m`; the
strict parity/unit gate is `Tests/test_adaptive_lower_frame.m`.

## 2026-07-10 integration coefficient-frame fast path

Two algebraic changes remove the current integration bottleneck without
changing WorkingPrecision, ExpansionOrder, epsilon windows, or tolerances:

1. `IntegrateLocalSolution` computes each primitive term's honest product
   window, takes the sum intersection once, and accumulates all finite Cauchy
   products directly into coefficient slabs.  This replaces one
   `ESNew`/`ESTimes`/ever-growing `ESAdd` chain per
   `(sector,Taylor-order,component)`.
2. `antiderivativeAtLog` now emits the coefficients of
   `Exp[b eps Log T] (m+1+b eps)^(-j-1)` directly from the finite binomial
   convolution.  It is algebraically identical to the former repeated
   `ESInvert`/`ESTimes` tower.  Same-sign `m=-1,b!=0` bounds use the manifestly
   combined interval formula, as required by `Integrate.md` 2.2.1, so no
   artificial Laurent row or CompleteMax loss is introduced.

Measured with the committed `Scripts/bench_chart.m` fixtures at WP120,
EO40, epsilon order 5, three repetitions/minimum:

| fixture/phase | pre-change s | direct slabs s | direct primitive s | total speedup |
|---|---:|---:|---:|---:|
| sunrise cold integrate | 0.435663 | 0.313533 | 0.029059 | 14.99x |
| sunrise warm integrate | 0.148134 | 0.024243 | 0.024288 | 6.10x |
| banana cold integrate | 0.607841 | 0.324169 | 0.032313 | 18.81x |
| banana warm integrate | 0.322845 | 0.028576 | 0.027656 | 11.67x |
| banana sliced scalar unit | 0.071462 | 0.029918 | 0.029073 | 2.46x |

On the saved production L2 checkpoint (WP1000, EO100, unchanged halo and
tolerances), banana master 1 fell from 213.6 s to 57.7 s: **3.70x
end-to-end** including all rational multiplications and 36-tile API overhead.
The exact differential coverage is in `Tests/test_integrate_fastpath.m`:
40 randomized honest-window contractions, 42 direct-vs-legacy primitive
cases, symbolic/arbitrary-precision inputs, exact-zero window constraints,
and same-sign analytically regularized pole sectors on both arms through log
depth 3.

## 2026-07-10 regular match-frame normalization

The restored `+1/DivisionOrder` / `-1/DivisionOrder` geometry removed the
excess segment count, but the WP1000 banana trace exposed a second,
independent conditioning cost: each fundamental basis is center-normalized,
while its weights are solved at the incoming edge.  The raw regular-chart
frame lost about 60–90 digits at every match and reached segment 9 with only
about seven digits, even though its epsilon lattice and truncation tail were
certified.

After lattice saturation, regular charts now use the constant right action
`P = F_eps0^-1` at the match point.  The transformed LocalSolutions are
re-evaluated and must pass identity, window, and repeat-saturation proofs.
This consumes no epsilon order and does not change tags or prescriptions.
On the same saved level-1 boundary, the observed per-match loss fell to about
28–35 digits; segment 9 entered with 391.8 digits and its weights retained
361.0 digits, clearing the former failure.  Singular/resonant frames retain
the existing saturated Laurent path to avoid forcing their unequal tagged
windows to a common top.

Date: 2026-06-12.  Baseline: master 1ce6038 (M5g-8).  Companion to
Docs/PerfGapAnalysis.md (structural old-vs-new analysis; its lever 1 —
value-vector transport — is being landed separately and is NOT
re-litigated here).  Everything below is grounded in MEASURED numbers
from `Scripts/bench_chart.m` on the committed FT level-1 fixtures
(`Tests/refs/bench/{sunrise_L1,banana_L1}.m`); the benchmark never loads
the FT layer.

## 0. Method

- Fixtures: the REAL level-1 systems (sunrise d = 3, banana d = 7, exact
  eps-rational in xx1, IBP factors included), frozen once by
  `Scripts/gen_bench_fixtures.m` from the same DefineFTIteration +
  RunFullIteration chain the runner uses.
- Settings: the FT-runner level-1 reality — WP 120, ExpansionOrder 40,
  EpsilonOrder 5, DivisionOrder 4, StepDivisionOrder 4.
- Per phase: AbsoluteTiming, 3 reps, min.  Full default run ≈ 3 min.
- The bench chart is the first interior REGULAR chart of the
  `SegmentLine[sys, {11/23, 0}]` plan — the marching steady state.
- Cross-RUN comparisons carry ~10-20% ambient load noise (the second
  licensed kernel often runs a marathon); within-run ratios are the
  trustworthy ones.

## 1. Measured baseline (library at 1ce6038) vs this branch

banana_L1 (d = 7; plan: 13 charts, 1 singular).  "opt" = after the §2
changes, same harness, same settings:

| phase | base s | opt s | note |
|---|---|---|---|
| Plan (SegmentLine, whole line) | 3.31 | 3.21 | |
| PrepareChartCold | 0.118 | 0.122 | 0.9% of chart |
| PrepareChartMemo | 8e-6 | 7e-6 | free |
| SolveChartCold | 9.90 | 9.97 | recursion+residual bound |
| — PrepCleared | 0.23 | 0.23 | |
| — Recursion1Col (x7 = 5.7) | 0.81 | 0.82 | 43% of chart — untouched |
| — Assemble1Col (x7 = 0.6) | 0.089 | 0.093 | |
| — ResidualCheck | 3.08 | 2.94 | see §2.3/R4 |
| SolveChartMemo | 3.1e-5 | 3.6e-5 | lo/hi anchor replay is free |
| EvalBasis (7 cols) | 0.206 | 0.214 | |
| MatchWeights | 2.59 | **0.037** | **x71** (§2.2) |
| Combine | 0.268 | 0.253 | basis arrays still exact-mixed |
| EvalCombined | 0.033 | 0.0099 | combined object now numeric |
| Validate (1 read) | 0.0066 | 0.0063 | |
| ErrorProbe | 0.066 | 0.021 | |
| **PerChartBasis (regular marching cost)** | **13.18** | **10.63** | x1.24 |
| LineEstimate (13 x regular) | 171 | 138 | |
| PrepareSingCold | 0.26 | 0.26 | |
| **SolveSingCold (endpoint x = 0)** | **68.6** | **67.7** | recursion/ladder bound |
| MultiplyRational (tile, full 7-vector) | 0.96 | **0.145** | x6.6 |
| IntegrateTile, cold | 6.05 | 4.75 | |
| IntegrateTile, warm | — | **3.54** | §2.1 memo |
| TileSliced1Comp (warm, 1 component) | — | **0.535** | the R1 prize |

sunrise_L1 (d = 3; 8 charts, 1 singular), base -> opt: PerChartBasis
1.30 -> 1.24 (SolveChartCold 1.14 -> 1.11 with ResidualCheck
0.294 -> 0.280); MatchWeights 0.025 -> 0.0037; MultiplyRational
0.30 -> 0.070; IntegrateTile 3.88 -> 1.47 warm (2.76 cold);
TileSliced1Comp 0.51; SolveSingCold 1.56 -> 1.49.

### What the profile says

1. **The tile phase, not transport, dominates the banana run.**  At
   baseline one tile integral on the 7-component object costs 7.0 s
   (0.96 multiply + 6.05 integrate).  `API.m LineIntegral` pays this PER
   nonzero cvec component PER lower master PER tile, then keeps ONE
   component of the 7-vector result (`ri2["Values"][[ci]]`).  Real-run
   scale: ~25 tiles (lo+hi kept charts) x ~5 masters x ~3-4 nonzero
   components ≈ 400-500 units x 7 s ≈ **45-60 min** — the core of the
   2 h banana run (L2/L3 tiles and limit cases add more of the same).
2. **Transports ≈ 8 min/level at L1**: honest line = 12 x 13.18 +
   singular ≈ 230 s, two lines minus the shared anchor solve.
3. Within a regular chart: d-column recursion 43%, residual self-check
   23%, MatchWeights 20%.  The per-center exact rational prep
   (PrepareChart + PrepCleared = 0.34 s) is only **2.6%** — the
   PerfGapAnalysis lever-2 projection (x1.5-2.5 from line-global prep
   amortization) is dead on this profile; coarse-rational centers
   already killed that grind.
4. **The pipeline was exact-rational bound, not bignum bound.**
   WP 80 vs WP 120 (adjacent runs): SolveChartCold x1.03, ResidualCheck
   x1.00, MatchWeights x1.00, PerChartBasis **x1.02**.  Harness meta on
   the MatchWeights inputs at baseline: **430 of 560** nonzero evaluated
   basis coefficients were EXACT rationals (max 528 B, compounding under
   field ops at ~50 µs/op — 25x a WP-120 bignum op).  The recursion's
   ByteCount <= 500 gate keeps small exact entries exact; evaluation at
   EXACT rational match points contracts them against exact t-powers,
   and the consumers ground exact-giant arithmetic.
5. **One singular chart ≈ 5-6 regular charts** (68 s).  After value
   transport lands (regular charts ~7x cheaper), the singular endpoint
   chart becomes THE per-line floor.
6. Memo hits are free (7-45 µs): the lo/hi shared-anchor replay and
   repeated-center preps cost nothing.

## 2. Landed in this branch (every number from the harness, same run)

### 2.1 Antiderivative memo in Integrate.m

`antiderivativeAt[m, b, p, T, kMaxOut]` is a pure function of exact
arguments (+ WP, which fixes the numericized log) yet was rebuilt per
call; the per-n towers are identical across the tile's components AND
across the ~20 (master, cvec-component) units that integrate the SAME
tile.  Added a bounded memo (flush at 4096 entries, the $shCache
pattern; entries are system-independent because the function is pure).
Measured: IntegrateTile 4.75 cold -> 3.54 warm (banana), 2.76 -> 1.47
(sunrise); the real runner is warm for every unit after a tile's first.

### 2.2 Numeric handoff at the matching/residual seams (the exact-giant killer)

`EvaluateLocalSolution` keeps its exact-in/exact-out module contract
(pinned by test_sectorseries t06-t16).  The PIPELINE, whose policy is
already "coefficients are the only numerics" (the runner numericizes
boundary values and level handoffs at WP+20), now applies that policy at
its own seams:
- `Transport.m TransportLine`: the evaluated incoming values (vvals) and
  the evaluated basis matrix (F) numericize at WP+20 before MatchWeights
  (`numHandoff`); weights, the combined object, and every downstream
  consumer then run on bignums;
- `Solve.m ODEResidualCheck`: the evaluated f / theta-f values
  numericize before the d^2 esTimes grid (`numV`).
Values only; tags, windows, structural decisions untouched; exact zeros
stay exact (N[0] = 0); symbols pass through.
Measured (banana, same-window): MatchWeights 2.59 -> **0.037** (meta:
exact nonzero coefficients 430 -> 0), EvalCombined 0.033 -> 0.0099,
ErrorProbe 0.066 -> 0.021, MultiplyRational 0.96 -> **0.145**,
IntegrateTile 6.05 -> 3.54 warm.  PerChartBasis 13.18 -> **10.63**
(x1.24); the tile unit 7.0 -> **3.69** (x1.9).  Side effect: the known
N::meprec storms at the d = 7 singular chart (campaign backlog, 4
warnings per bench run) are GONE — they originated in numeric
comparisons over the exact-giant residual/limit content.

### 2.3 ODEResidualCheck B(t0) hoist (hygiene)

The d^2 `Together[entry /. t -> t0]` + eps-expansion build ran once per
column (7x at banana L1).  Hoisted + memoized per eps-window; identical
values.  Measured saving small (~0.1-0.4 s, partly masked by load
noise): the check's real remaining cost is the 7 x
(DifferentiateLocalSolution + 2 evaluations) on exact-mixed coefficient
arrays — see R4.

Validation: test_eps_series, test_sectorseries, test_solve,
test_transport, test_integrate, test_api, test_pins green (see commit).

## 3. Ranked levers (factors against the post-§2 profile)

### R1 — Tile component slicing in LineIntegral / limitCombined (x6-7 on the dominant remaining phase; ~30 lines; zero numerics risk)

**Measured basis (same run, post-§2).**  Warm full-vector tile unit =
0.145 (multiply) + 3.54 (integrate) = 3.69 s; warm 1-component
projection of the SAME work (`TileSliced1Comp`) = **0.535 s** — x6.9.
`sumSectorIntegrals*` spends the difference on the per-(n, component)
EpsSeries build/multiply over components the caller throws away:
LineIntegral keeps `Values[[ci]]`, the runner's limitCombined keeps
`[[j]]`.

**Sketch.**  In the per-(tile, ci) loop, project the LocalSolution to
component ci before `MultiplyRational` (per sector:
`Coeffs[[All, All, {ci}]]` — the dims invariant holds, metadata
unchanged), integrate, read `Values[[1]]`.  Arithmetic on the kept row
is IDENTICAL.  Same surgery in `limitCombined` / `EndpointLimitValues`.

**Arithmetic.**  Banana L1 tile phase ≈ 450 units x 3.69 ≈ 28 min
(post-§2) -> 450 x 0.54 ≈ **4 min**.

**Risk.**  Low.  PV pairing and the divergence-gate verdicts run per
component already; two flagged deltas: (i) the gate's relative-zero
SCALE currently spans all components — slicing narrows it to the kept
one (stricter, hence still honest); (ii) a divergent-but-UNUSED
component no longer aborts a unit that never consumes it (the drop-rule
design only guards consumed values).  A/B: the harness rows + one
`FT_EXAMPLES=sunrise` STEPWISE diff.

### R2 — Two-seat parallelism: split each level's lo/hi-line work across the 2 licensed kernels (x1.8-1.9 on everything after FIRE; zero numerics risk)

The license allows two concurrent kernels (campaign-verified).  Both
halves of a level split by line: trLo/trHi transports are independent
pure functions of (sys, bvals, anchor, extraFacs); every tile belongs
to exactly one kept chart (lo or hi set); per-master tile sums are
associative.  FIRE6 is a C++ subprocess and holds no kernel seat, so
this is purely about the WL phases.

**Sketch.**  `Scripts/ft2_worker.m` reads a Put spec {matrix, var,
extraFacs, bvals, from, to, config, cvecs}, computes its line's
transport + its tiles' per-master partial sums, Puts the result; the
parent runs the lo half in-process and adds partials.  Serialization ≈
2-4 MB exact text ≈ 1-3 s — noise.  The worker runs only while the
parent waits (never a third kernel; acquire the same /tmp lock
convention before spawning).

**Arithmetic.**  Post-R1 banana level ≈ transports 6.5 min + tiles 4 min
-> **~5-6 min wall**.  Multiplies with every other lever and with value
transport.

**Risk.**  Orchestration only: seat contention (lock), worker crash ->
serial fallback, Put/Get bignum fidelity (exact InputForm round-trip).
Numerics identical.  A/B: per-level wall + STEPWISE diff under
FT2_PARALLEL=0/1.

### R3 — Singular-chart cost program (the emerging per-line floor)

SolveSingCold 67.7 s vs 10.6 s regular — and §2 did NOT move it: the
cost is the recursion/ladder machinery, not the seams.  After value
transport lands (regular ~7x cheaper), a banana line ≈ 12 x ~1.5 + ~68 —
the singular chart is then ~75-80% of the line.  Next profiling pass
(point the harness's Private-seam phases at the singular chart):
- frame width W = K + 3 Pmax + 2 cdMax + 7 vs K + 7 regular: every
  frConv/frInv is O(W^2); audit Pmax/cdMax against the chart's actual
  log ceiling for over-conservatism;
- CASE R/P ladders multiply per-n work by (P+1) l-levels per column;
- applyGauge (FuchsianReduce charts) full-array convolutions;
- its residual check evaluates d columns x multiple sectors (R4 helps).
Even x2 here is x1.6 on a post-value-transport line.  Measurement-first;
no invariant changes proposed.

### R4 — Recursion/array-content numericization (extend §2.2 upstream; x1.3-1.6 on the chart)

Post-§2 the regular chart is again solver-bound: 7 x 0.82 = 5.8 s
recursion + 2.9 s residual of 10.6 total.  Both residues are the SAME
disease one level deeper: the recursion's U-content and the basis
arrays keep sub-500-byte EXACT entries (the ratEpsList gate); frConv
(ListConvolve) on mixed exact lists does exact work that grows with n,
and the residual check pays 7 x (DifferentiateLocalSolution + 2
evaluations) on those exact-mixed arrays (why §2.3 bought so little).
WP 80 vs 120 = x1.02 proves digits are not yet the cost.
Options, honesty-ordered: (a) numericize recursion INIT and U-content at
WP+20 above a much smaller ByteCount inside runRecursion (pure value
policy; tags/frames untouched); (b) with (a) in place the frame lists
become homogeneous bignum -> packed-array slabs for the Dot/shift
kernels (the M6 target in the campaign notes).  A/B: Recursion1Col +
ResidualCheck rows; pins for value drift within MatchTol.

### R5 — Line-global prep reuse (Taylor-shift the line once): DEFER

Measured PrepareChartCold + PrepCleared = 0.35 s = 2.6% of the baseline
chart (3.3% post-§2) ⇒ x1.03 at best today.  Revisit only after value
transport lands (chart ≈ 1.5 s makes it ~20% ⇒ x1.2): numericize the
line's cleared num/den coefficient polynomials once per line and
Taylor-shift packed arrays per center (O(deg^2) numeric vs exact
Cancel@Together + PolynomialLCM + ratEpsList per center), with the
regular/singular DECISION still made exactly against the factored root
list (I-6 intact).

### R6 — eps-first cascade: REJECT (measured confirmation)

Measured FrameW = 12 at K = 5 ⇒ the coupled symbolic-eps frame costs
<= W/(K+1) = 2.0x per delivered eps order vs per-order numeric systems —
before the cascade pays its own (K+1)-fold source assembly and loses
the solve-once-per-segment basis economy.  Value transport already
removes the d-column multiplicity on regular charts; a second per-order
solver would violate the ONE-solver design (Config $droppedReasons,
RewritePlan I2) for <= x1.5.  Keep only the micro-variant: CASE-T steps
with a single eps^0 delta entry can replace frInv + one frConv with a
scalar scale (O(W^2) -> O(W)); bounded by the Recursion1Col share ⇒
<= x1.1, fold into R3/R4 frame work.

### R7 — WP staging: DEAD on measurement (x1.02) — re-test after R4

WP 80 vs 120 changed PerChartBasis by 2% (§1.4): the cost was exact
content, not digits.  §2.2 made the SEAM values numeric; R4 would make
the recursion numeric too — only then can staged WP (regular interior
charts at lower WP, full WP near singular charts; ErrorEstimate
accumulation and the DigitBudget E3 gate already exist) pay its
theoretical x1.3-1.5.  Not before the content levers.

### R8 — Small fry (bundle, don't standalone)

- SegmentLine 3.2-3.3 s/line (banana): RootReduce-heavy dedup per line;
  FindSingularities depends only on (sys, extras) — cache per level:
  -6 s/banana level.
- Validate-on-every-read: 6.3 ms x ~24 reads/chart ≈ 0.15 s/chart
  (1.4%); the construction-time choke point (PerfGap L4) is right but
  small here.  Bundle with any SectorSeries rework.
- Combine still 0.25 s (exact-mixed basis arrays x numeric weights) —
  falls out of R4 automatically.

## 4. Where this leaves banana L1 (per level, measured-derived)

| stage | tiles (~450 units) | transports (2 lines) | L1 wall |
|---|---|---|---|
| baseline 1ce6038 | ~50 min | ~8 min | ~58 min |
| **+ §2 (landed here)** | **~28 min** | **~7 min** | **~35 min** |
| + R1 slice | ~4 min | ~7 min | ~11 min |
| + R2 two-seat | ~2 min | ~3.5 min | **~6 min** |
| + value transport (other agent) | ~2 min | ~12 x 1.5 s + 68 s/line ≈ 3 min | ~5 min |
| + R3 singular + R4 content | ~1.5 min | ~1.5 min | **~3 min** |

L2/L3 levels, the limit cases and the FT layer add minutes more; the
remaining gap to the old engine's ~30 s is per-op machinery
(PerfGapAnalysis §4) — R4(b) packed slabs is the structural answer.

## 5. How to A/B anything here

1. `Scripts/bench_chart.m`: one kernel, ~3 min, JSON rows; scope with
   BENCH_FIXTURES / BENCH_REPS / BENCH_WP / BENCH_EPS_ORDER; trust
   within-run ratios over cross-run (ambient load noise 10-20%).
2. Runner-level claims (R1, R2): `FT_EXAMPLES=sunrise` stepwise run,
   STEPWISE rows must match to MatchTol; banana for wall-clock.
3. Anything touching values (R4, R7): the 7-suite battery
   (test_eps_series/sectorseries/solve/transport/integrate/api/pins) +
   a bubble/sunrise oracle diff before promoting — the §2.2 change went
   through exactly that gate.
