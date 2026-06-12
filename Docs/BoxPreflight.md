# Box readiness report (DiffExp2 M5 ladder preflight)

Generated 2026-06-12 by the box preflight task; zero full runs were spent.
Method: vendor the static per-level artifacts (FIRE + DiffMatrix + IBP
factors) once, then audit geometry/indicial structure offline — the lesson
from banana losing three 20-minute runs to plan-geometry errors.

Artifacts (committed):

- `Tests/refs/bench/box_L{1,2,3}.m` — per-level fixtures: seam-normalized
  DiffMatrix, CollectLevelIBPSingularFactors, masters, boundary requests,
  normalized FIRE reductions, FT settings.  Generator:
  `Scripts/box_preflight_fixtures.m` (kernel + FIRE, ~18 s).
- `Scripts/box_preflight_audit.m` — plans + ValidatePlan + indicial audit
  from the fixtures (kernel, no FIRE, ~5 s).  Raw log: AUDIT JSON rows.
- `Scripts/box_preflight_repro.m` — transport reproducer for the findings
  below (kernel, no FIRE).
- `Scripts/box_verdict.py` — the one-command GO/NO-GO judge for a completed
  box run (no kernel).

## 1. Verdict summary

| level | line geometry (ValidatePlan) | structure | verdict |
|-------|------------------------------|-----------|---------|
| L3 (d=1, xx3) | both plans OK, 4+4 charts, 0 E8 | endpoint 1 REGULAR (limitUpper takes the `Value` path); x^(1-eps) sector at 0 | GO |
| L2 (d=2, xx2) | both plans OK, 6+10 charts, 0 E8 | APPARENT singularity at xx2 = 7/11 crossed by the hi line; its local solution carries t^(-1+eps) | **NO-GO — blocker B1** (runtime crossing gate, not a geometry error) |
| L1 (d=2, xx1) | both plans OK, 13+4 charts, 0 E8 | interior IBP pole at 1/4 (regular-indicial chart); first interior-PV integration exercise | GO with watch item W1 |

Box CANNOT complete until B1 is fixed: the runner needs the L2 hi transport
for BOTH of its boundary cases (integrate and limitUpper), so the run dies
at level 2 after level 3 completes.  Everything else is green: all FOUR
SegmentLine plans validate statically with zero E8s, no FuchsianReduce
charts exist anywhere (all pole orders <= 1 — no t05 class), the eps-depth
budgets (7/6/5) follow the post-c0b24f3 full-depth rule, and the reference
wiring is verified end-to-end against the vendored oracle log.

## 2. Fixtures and plans (FT settings: WP 120, ExpansionOrder 40, DivisionOrder 4, classic stride 4, anchor 11/23)

Level/case structure (from the fixtures; matches the old-core
`Tests/refs/oracle_logs/l2_box.log` exactly):

| level | var | d | masters | cases below | extra (IBP) factors | matrix denominator factors |
|---|---|---|---|---|---|---|
| 3 | xx3 | 1 | {1,0,0,0} | integrate {2,0,0,0}; limitUpper {1,0,0,0} | xx3, 23+10·xx3 | xx3, 23+10·xx3 |
| 2 | xx2 | 2 | {1,0,0,1},{1,0,0,0} | integrate {2,0,0,0}; limitUpper {1,0,0,1} | xx2, xx2−1 | xx2, **11·xx2−7**, xx2−1 |
| 1 | xx1 | 2 | {1,0,1,0},{1,0,0,1} | integrate {2,0,1,1} | xx1, **4·xx1−1**, xx1−1 | xx1, xx1−1 |

ValidatePlan: **OK on all four plans** ([11/23 → 0] and [11/23 → 1] at every
level).  Every incoming match point lies strictly inside the producing
chart's disk; the two singular-handoff (FixWithin) clips that exist — into
7/11 on the L2 hi line and into 1/4 on the L1 lo line — are geometrically
possible.  DigitsNeeded = 60 everywhere (AccuracyGoal "?", WP 120: ample).

Chart chains (centers; S = singular chart; full exact radii + match-point
chains vendored in `Tests/refs/bench/box_plan_audit.log` AUDIT JSON rows,
reproducible in ~5 s via `Scripts/box_preflight_audit.m`):

- L3→0: 0.478, 0.375, 0.300, **0(S)** — 4 charts.
- L3→1: 0.478, 0.600, 0.714, 0.833 — 4 charts, endpoint REGULAR
  (`EndpointIsSingular` False: 23+10·xx3 has its root at −2.3; nothing at 1).
- L2→0: 0.478, 0.375, 0.300, 0.235, 0.1875, **0(S)** — 6 charts.
- L2→1: 0.478, 0.514, 0.538, **7/11(S)**, 0.714, 0.762, 0.8125, 0.85,
  0.88, **1(S)** — 10 charts; the 7/11 chart is CROSSED (radius 4/11,
  far-side handoff at 0.695 — well inside the disk).
- L1→0: 0.478, 0.429, 0.393, 0.364, 0.341, 0.321, **1/4(S)**, 0.200,
  0.167, 0.133, 0.105, 0.083, **0(S)** — 13 charts; 1/4 is crossed.
- L1→1: 0.478, 0.571, 0.667, **1(S)** — 4 charts.

Total 41 charts across all 6 lines (transports are computed once per level
and reused by every master via the runner's trLo/trHi + PrecomputedCharts
cache).

## 3. Structure audit

Indicial data at every singular chart (exact, from `ChartIndicial` via
`PrepareChart`):

| chart | pole order | spectrum (a, b, mult) | family class | notes |
|---|---|---|---|---|
| L3 @ 0 | 1 | (1, −1, 1) | Single | sector x^(1−eps); no resonance |
| L2 @ 0 | 1 | (0, −1), (1, −1) | TrueResonant, LogMax 1, 1 collision | integer-spaced same-b: CASE R log ladder (exercised class) |
| L2 @ 7/11 | 1 | (−1, 1), (0, 0) | Pseudo, JointSolve, 1 LaurentShift collision | THE apparent singularity (campaign pin: old-core eigenvalue −1 at eps=0) |
| L2 @ 1 | 1 | (0, 0), (1, −1) | Pseudo, JointSolve | the exact L2 endpoint-tower chart (closed form Γ(1+eps)B(1−eps,1−eps)(4/23)^(−eps), PINS §3.1) |
| L1 @ 1/4 | **0 (REGULAR)** | (0, 0) ×2 | Single | IBP-only singular chart: matrix regular, crossing is a trivial phase |
| L1 @ 0 | 1 | (0, −1), (0, 0) | Pseudo, JointSolve | classic dimreg x^(−eps) ⊕ 1 mix (bubble/sunrise class) |
| L1 @ 1 | 1 | (0, −1), (0, 0) | Pseudo, JointSolve | same |

- FuchsianReduce/gauge: **never triggers** (no pole order ≥ 2 anywhere) —
  box has no banana-L1-style nilpotent double poles, no t05-class charts.
- EpsZeroDegeneracy = 0 at every singular chart (RecombineBasis stays
  inert).  No EpsDegenerate families.
- Apparent-singularity check (matrix denominator roots NOT among IBP
  singular factors): L2 root 7/11 of 11·xx2−7 is the ONLY one (L3: none —
  matrix and IBP factor sets coincide; L1: matrix roots {0,1} are a SUBSET
  of the IBP roots {0, 1/4, 1}, the 1/4 chart exists only through
  ExtraSingularFactors).  7/11 ≈ 0.636 sits INSIDE (anchor, 1), 0.158 from
  the anchor — comfortably far (banana's x = 1/2 sat 0.022 away).
- Singularity layout relative to [0,1] and anchor 11/23 ≈ 0.478:
  L3 real roots {0, −23/10}: only 0 in range; upper endpoint regular.
  L2 real roots {0, 7/11, 1}: one interior point (hi side).
  L1 real roots {0, 1/4, 1}: one interior point (lo side).
  No roots are NEAR the anchor (min distance 0.158) and none are complex.
- Boundary-case inventory (BoundaryRequestRecords): L3 integrate +
  limitUpper, L2 integrate + limitUpper, L1 integrate only.  No
  limitLower, no direct cases.  The L3 limitUpper is the first
  REGULAR-endpoint limit in the DiffExp2 ladder — it takes the
  `tr["Value"]` branch (TransportLine evaluates the final chart at the
  endpoint when `EndpointIsSingular` is False; verified present in
  Transport.m's return record).

## 4. Blockers and watch items

### B1 (BLOCKER): crossing the apparent chart at xx2 = 7/11 hits the multivalued-crossing E8

Statically derived from Transport.m + SectorSeries.m and confirmed by
runtime reproducer (see §6):

1. The L2 hi line must CROSS the singular chart at 7/11 (plan above);
   TransportLine's crossing branch runs
   `sigma = ChartImSign[current]` and, when sigma is not ±1, raises
   `E8 "crossing a multivalued singular chart without a derivable Im-sign
   (missing DeltaPrescriptions)"` if ANY sector has non-integer a, b ≠ 0,
   or p > 0 (Transport.m crossing gate).
2. `ChartImSign` derives ONLY from `ls["Prescriptions"]` (SectorSeries.m),
   SegmentLine builds every chart with `"Prescriptions" -> {}`, and the box
   runner configures no DeltaPrescriptions → sigma = None.
3. The 7/11 chart's local solution carries the (a,b) = (−1, 1) sector —
   i.e. t^(−1+eps), b ≠ 0 → the gate fires.  The PHYSICAL weight of that
   sector is ~0 (the transported function is analytic there — the campaign
   pin), but the gate is SYNTACTIC on sector tags:
   `CanonicalizeLocalSolution` drops exact-zero sectors only (F10), and the
   matched weight is numerically tiny, never the exact integer 0.
4. No validated example ever exercised this: bubble/sunrise crossed no
   b ≠ 0 singular charts (their interior lines are clean), banana's L3/L2
   passes did not either.  Box L2 is the first.

Cross-prediction worth checking on seat 1: the in-flight banana run's L1
line crosses the x = 1/2 resonant point; if that chart's sectors carry
p > 0 or b ≠ 0, the SAME gate fires there (grep its log for
`E8` / "crossing a multivalued singular chart").

Suggested fix (one site): make the crossing gate magnitude-aware — filter
`current["Sectors"]` through the established numerical-negligibility test
(the `ESCoeffZeroQ`/LaurentLeadTol relative gate that MatchWeights' mwTrim
already uses, scale = the object's coefficient scale) BEFORE the AnyTrue
multivalued test, and let ApplyCrossing skip negligible sectors.  A sector
that carries no weight needs no branch choice.  Genuine branch points with
REAL weight still error exactly as designed (the pentagon/DEC-21 gate is
untouched).  Alternative unblock (no code change, uglier): pass a
DeltaPrescription for 11·xx2−7 with either sign — for a zero-weight sector
the sign cannot matter — but that pollutes the FT runner with a per-example
config and normalizes prescribing APPARENT factors; prefer the gate fix.

Minimal reproducer (~40 s kernel, no FIRE, no boundary data needed):

    env WolframKernel='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel' \
      wolframscript -file Scripts/box_preflight_repro.m
    # expected before the fix: "TRANSPORT L2 -> 1: FAIL E8 ... crossing a
    # multivalued singular chart"; all five other transports OK.

### W1 (watch): L1 integrate is DiffExp2's first interior-PV exercise

The L1→L0 integrand {2,0,1,1}(t) has IBP-coefficient poles at t = 0, 1/4, 1
with the 1/4 pole INSIDE the integration range (campaign first; old core
needed the complex-log-leakage + Hadamard finite-part fixes here, commit
9aeb300).  In DiffExp2 the 1/4 tile becomes an interior-PV
IntegrateLocalSolution call (PV parity (−1)^m, Log 2 pin) on a
regular-indicial chart.  Geometry and indicial data are clean; the risk is
purely in the PV/divergentGate value path, which unit tests cover but no
example has driven end-to-end.  Triage refs if L0 eps^0 comes out wrong
while poles are right: the pointwise transport pins straddling 1/4
(`Tests/refs/pysecdec/results_box_t*.json`, PINS §4) — eps^-1 = −1200/19
exactly at t = 1/20 AND t = 19/20, −25/2 at t = 2/5.  Also watch for
imaginary-part leakage at L0 eps^0 (gate: |Im| < 1e-10 in the verdict).

### Non-blockers verified

- Budgets: EpsilonOrder 7/6/5 at L3/L2/L1 (= FT_EPS_ORDER 0 + level + 4)
  with full incoming depth carried per level (the c0b24f3 rule the runner
  already implements; box_bubble's starvation class cannot recur here).
- All six transports are needed and planned (integrate at every level
  forces lo+hi); no limitLower/direct paths are exercised at all.
- No FuchsianReduce, no EpsDegenerate recombination, no resonant
  TrueResonant family except L2@0 (same-b integer ladder, the standard
  CASE R class covered by su08-class units and sunrise/banana endpoints).

## 5. Reference wiring (verified, no kernel)

- Old-core oracle: `Tests/refs/oracle_logs/l2_box.log` — full 50-digit
  STEPWISE rows for L3/L2/L1/L0 plus FINAL; values verified exact against
  independent sympy/mpmath quadrature per PINS §3.2.  NOTE its FINAL row
  prints `Finite -> 12.0` — the OLD runner's "Finite" is the leading
  (min-power) coefficient, while run_ft_stepwise2.m prints the true eps^0
  coefficient.  Judge FINAL only with the new semantics.
- pySecDec: `Tests/refs/pysecdec/box1l_pin_result.log` (THE L0 pin
  {12, −0.334914246809752, −41.28416739757452, −70.1654426808142});
  `results_box1l.json` (family route, levels 0–3 at t = 11/23);
  `results_box_t{1_20,2_5,19_20}.json` (pointwise {2,0,1,1} triage pins).
- `compare_stepwise_log.py --refs results_box1l.json` verified on the
  vendored oracle log: it checks ONLY eps^-1 per row (the pySecDec summary
  constant term has no `*eps^0` suffix, so TERM_RE skips eps^0
  everywhere), and it reports EXACTLY one expected failure — the
  known-garbage box_L0_1x1x1x1 entry (ref −175.59; MANIFEST trust notes).
  Its exit code is therefore permanently 1 on a perfect run: NOT usable as
  a verdict.
- PINS→runner field mapping at FT_EPS_ORDER=0: STEPWISE `Coefficients`
  pairs [k, c_k] for k = RawMinPower..0 (the printer caps display at
  eps^0, so the eps^1 pin −70.165… is out of scope for this runner —
  PINS §6.6); L0 row carries eps^-2/-1/0; FINAL `Finite` = the eps^0
  coefficient = the pin value −41.28416739757452.
- `Scripts/box_verdict.py` (new, validated against the oracle log: 13/13
  STEPWISE pins PASS; the single FAIL it reports there is the old-FINAL
  semantic noted above, which is correct behavior): pins every level
  against the 50-digit oracle values at 1e-9 relative, L0 against the
  pySecDec pin (eps^0 at 5e-9 absolute — the pin's own Cuhre floor), and
  enforces real-valuedness.

THE VERDICT COMMAND (run + judge) once B1 is fixed:

    env WolframKernel='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel' \
      FT_EXAMPLES=box FT_EPS_ORDER=0 FT_WORKING_PRECISION=120 \
      wolframscript -file Scripts/run_ft_stepwise2.m > /tmp/ft2_box.log 2>&1
    python3 Scripts/box_verdict.py /tmp/ft2_box.log
    # exit 0 = GO (all 14 pins), exit 1 = NO-GO with per-order diffs

Optional cross-check (informational; expect exit 1 with EXACTLY one FAIL,
the L0 garbage ref): `python3 Scripts/compare_stepwise_log.py
/tmp/ft2_box.log --refs Tests/refs/pysecdec/results_box1l.json`.

## 6. Runtime confirmation (reproducer, generic boundary values) and estimate

`Scripts/box_preflight_repro.m` ran all six transports at FT settings
(17.6 s kernel total, log /tmp/box_preflight_repro.log, 2026-06-12):

| transport | result | charts | wall |
|---|---|---|---|
| L2 → 1 | **FAIL E8 at seg5@0.714286** — `"crossing a multivalued singular chart without a derivable Im-sign (missing DeltaPrescriptions)"` — B1 confirmed at exactly the predicted chart | failed at 5/10 | 2.0 s |
| L1 → 0 | OK (the 1/4 regular-indicial crossing PASSES, as predicted) | 13 | 4.3 s |
| L1 → 1 | OK | 4 | 1.3 s |
| L2 → 0 | OK | 6 | 3.0 s |
| L3 → 0 | OK | 4 | 0.6 s |
| L3 → 1 | OK (regular endpoint, `EndpointIsSingular` False) | 4 | 0.6 s |

Per-chart cost at d ≤ 2, WP 120, eo 40: ~0.15–0.7 s (vs sunrise's 1–4 s at
d = 3).

### W2 (watch, cosmetic): meprec/ztest1 warning storms at L2 charts

During BOTH L2 transports the log carries `N::meprec` (4) and
`PossibleZeroQ::ztest1` (4, then General::stop) warnings on GIANT exact
rational + Log[...] combinations — the known exact-giant zeroQ class
(campaign ledger; banana L1 tolerates 4 of the same).  L1/L3 are clean.
Transports succeed regardless; this is a perf/cleanliness tax tied to the
L2 resonant/pseudo charts, not a correctness gate.  If the box run grinds
at L2, this is the first suspect (`grep -c meprec <log>`).

### Full-run estimate (after the B1 fix)

FIRE + iteration ~18 s (measured in phase A) + transports ~15 s (measured
above; one two-way transport per level serves all masters via the runner's
chart cache) + LineIntegral tiles: 8×1 (L3) + 16×1 (L2) + 17×2 (L1) = 58
component-tile integrations at ~0.3–1.5 s each (sunrise's d=3 tiles ran
~2 s) ≈ 20–90 s + two endpoint limits + handoffs (cheap).
**Total ≈ 1–3 minutes** at FT_EPS_ORDER=0, WP 120 — between bubble (13 s)
and sunrise (~1.5 min) per level-count scaling, far below banana (~12 min).

## 7. Regeneration recipes

    # fixtures (kernel + FIRE, ~20 s):
    wolframscript -file Scripts/box_preflight_fixtures.m
    # static audit: plans, ValidatePlan, indicial (kernel, no FIRE, ~5 s):
    wolframscript -file Scripts/box_preflight_audit.m
    # transport reproducer for B1/W1 (kernel, no FIRE, ~40 s):
    wolframscript -file Scripts/box_preflight_repro.m

(all with the Wolfram Player kernel env prefix as in §5.)
