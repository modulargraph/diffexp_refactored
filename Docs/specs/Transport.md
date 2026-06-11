# Module spec: DiffExp2/Transport.m

Status: M0 deliverable, implements RewritePlan.md section 3.2 (Transport.m,
~600 lines).  This document is a contract.  An implementation agent must be
able to build the module from THIS file, Docs/RewritePlan.md, and the cited
old code alone.  Old-code citations are against the frozen oracle tree at
the repo root (DiffExp.m + DiffExp/, FeynmanTrick/), current master.

---

## 1. PURPOSE

Transport.m owns everything between "a loaded system plus boundary data at
one point on a line" and "a LocalSolution (or numeric values) at another
point on that line": it locates ALL singularities of the system on and near
the line exactly (matrix-denominator factors plus configured
ExtraSingularFactors, complex roots included, with chart radii equal to true
complex-plane distances), plans the chain of expansion charts (predivision,
two passes over ONE shared geometry function, with the segment-count digit
budget), and runs the marching loop that at each chart solves (via Solve.m),
matches the incoming data onto the chart's fundamental basis with an
eps-GRADED Laurent weight solve (recombining the basis unimodularly at eps=0
where sector tags collide, asserting ord_eps det F = 0 afterwards), applies
the explicit crossing operator (phase times unipotent log-chain mixing,
convention pinned to the old principal-branch-far-side rules) when the path
passes a singular chart center, runs the two-point full-vs-reduced-order
error probe per segment, and maintains honest EpsWindow/TWindow/
ErrorEstimate metadata end-to-end with no silent fallback anywhere.

---

## 2. PUBLIC SYMBOLS

All symbols live in `` DiffExp2`Transport` `` and are exported with ::usage.
"EpsSeries" below means the truncated eps-Laurent array object of
Docs/specs/EpsSeries.md (window-carrying); "exact" means an exact Wolfram
expression (Rational, algebraic number, `Root`), never a float.

### 2.1 `FindSingularities[sys_Association] -> Association`

Input: the LoadSystem record (API.md contract); the fields consumed are
`sys["SingularFactors"]` (the irreducible x-dependent factors of all matrix
denominators, the analogue of `DiffExp`State`MatricesIrreducibleFactors`,
old MatrixLoading.m:217-231) and `sys["ExtraSingularFactors"]` (config; the
FT layer's IBP-coefficient denominators, old
FeynmanTrick/DiffExpIntegration.m:204-222 `appendMatrixFactors` and
:1326-1385 `CollectLevelIBPSingularFactors`).

Output:
```
<|
  "All"    -> { exact roots, complex included, deduplicated EXACTLY },
  "Real"   -> { the subset with exactly real value },
  "Factors"-> { factor -> {its roots} associations, for error messages }
|>
```

Behavior contract:
- Roots are obtained by exact `Solve`/`Roots` on each factor in the line
  parameter.  NO `Quiet` (old code Quiets, LineSegmentation.m:76).  A factor
  that cannot be solved exactly (non-algebraic in x) is a loud error E1.
- Deduplication is exact (`SameQ`/`RootReduce`), NOT numeric chop
  (old code dedups by `N[#, ChopPrecisionVal]`, LineSegmentation.m:77 — this
  can merge distinct nearby roots; box_triangle's interior pole at
  12167/12651 ~ 0.9617 sits 0.038 from the endpoint 1 and must never merge).
- Complex roots are RETAINED as complex numbers.  The old ghost projection
  (each complex root projected to real points Re and Re±Im unless another
  singularity already lies in that interval, LineSegmentation.m:81-99) is
  NOT reimplemented (forbidden fallback F4); complex roots enter only
  through `ChartRadius`.

### 2.2 `ChartRadius[center_, allSingularities_List] -> radius`

True complex-plane distance: `Min[Abs[center - s]]` over all s in
"All" with s != center, evaluated at >= WorkingPrecision with the exact
inputs retained alongside.  Returns the exact expression when cheap
(rational centers, rational/quadratic roots), else a >= WP-precision number.
This replaces the old real-line-only `minDistance = Min[at - leftBound,
rightBound - at]` over PROJECTED singularities (Mobius.m:50-61).
RewritePlan 3.1: "Radius -> distance to nearest singularity IN THE COMPLEX
PLANE".

### 2.3 `SegmentLine[sys_Association, {from_, to_}] -> SegmentPlan`

Pass 1 of the predivision strategy: pure geometry, no matrix expansion, no
solving.  Replaces the old dry-run block (old Transport.m:666-750, which is
a hand-maintained DUPLICATE of the marching code — comment at old
Transport.m:666; the duplication itself is forbidden, see F9).

Behavior contract:
- Direction = Sign[to - from]; bidirectional (to < from supported, as old
  Transport.m handles throughout, e.g. :691-694, :1112-1116).
- Singular charts are placed at every member of `FindSingularities["Real"]`
  lying strictly between from and to, plus at `to` itself when `to` is a
  singularity (EndpointIsSingular mode, old Transport.m:605-613).  Endpoint
  membership is tested EXACTLY (old code uses chopped numeric MemberQ,
  Transport.m:606-609 — forbidden, F13).
- Regular chart centers are produced by `NextCenter` (2.4) so that every
  match point sits at radius/DivisionOrder of BOTH adjacent charts.
- Match points: each chart's outgoing match point is at chart-coordinate
  distance Radius/DivisionOrder from its center on the marching side
  (old Transport.m:679, 1041-1046: interval ±RoC/DivisionOrder in the
  rescaled chart coordinate); when the next chart is singular the match
  point is additionally clipped into the intersection of both charts'
  validity intervals and centered there (old Transport.m:1124-1136
  `FixWithin` midpoint rule) — preserve the midpoint rule.
- `DigitsNeeded` is computed from the resulting segment count (2.5) and
  stored in the plan.
- The plan is the SINGLE source of truth: `TransportLine` executes it
  verbatim and asserts it never deviates (invariant I1).

### 2.4 `NextCenter[prevMatch_, singularities_, direction_] -> center`

The GetCPL/GetCPR geometry, ported from old Mobius.m:98-142 (affine charts)
and Mobius.m:72-90 (Mobius charts): given the previous match point x_b and
its two neighbouring singularities {z_min, z_max}, solve for the next
center x_new such that x_b sits at exactly (new chart radius)/k of the NEW
chart, k = DivisionOrder.  Closed form (affine, finite bounds, right-moving,
old Mobius.m:110-119): with s the new radius,
`s = k(z_max - x_b)/(1+k)` when the new chart is capped on the right (and
the mirror branch when capped on the left), `x_new = x_b + s/k`.
This both-adjacent-charts guarantee is the actual ill-conditioning defense
for the matching solve (RewritePlan R3; legacy review finding 12).
DivisionOrder: classic default 3 (State.m:112), FT pins 4
(FeynmanTrick/DiffExpIntegration.m:272).

NEW vs old: the bounding singularities for the s-solve are the REAL
on-path cap points, but the resulting radius must additionally be capped by
`ChartRadius` (complex distance).  If the complex cap is the binding one,
re-solve the 1/k placement against the capped radius (the geometry formula
holds with z-cap replaced by center ± capped radius).

### 2.5 `DigitBudget[accuracyGoal_, segmentCount_] -> Integer`

```
DigitsNeeded = accuracyGoal + Ceiling[Log10[segmentCount]] + safetyDigits
```
with safetyDigits from Tolerances.m (old value 2, `ISafetyDigits`,
State.m:215).  Old code: Transport.m:758-760.  Additional NEW contract:
if `DigitsNeeded + chopReserve > WorkingPrecision` (chopReserve from
Tolerances.m; old code keeps WP - ChopPrecision = 250 of reserve,
State.m:109,135), LOUD ERROR E3 — the old code never checks and silently
loses digits on long lines.

### 2.6 `TransportLine[sys_Association, boundary_, plan_SegmentPlan] -> TransportResult`

The marching loop (pass 2).  `boundary` is either
(a) a LocalSolution anchored at `plan["From"]` (general case: includes
    asymptotic/singular boundary data — the classic "line bcs" of old
    LineSegmentation.m:148-236 arrive here already normalized by API.m), or
(b) a plain values array `vals[[integral, epsorder]]` at a regular `from`
    (wrapped internally into the one-sector regular LocalSolution).
Symbolic indeterminate coefficients in `boundary` are allowed and tracked
per-indeterminate by the error probe (old Transport.m:565-569, 947-962);
they force UsePade off in SectorSeries evaluation (old Transport.m:566-569).

Loop body per chart (in plan order):
1. `basis, particular, meta = DiffExp2`Solve`SolveChart[sys, chart, ...]`
   — fundamental LocalSolution basis + particular + metadata, including
   `meta["CouplingDepth"]` (longest chain in the coupled-block dependency
   DAG; old proxy `MaxCouplingOrder` = largest coupled block,
   MatrixLoading.m:357-385, esp. :383).
2. v = incoming data evaluated at chart's MatchIn point (SectorSeries
   evaluate on the previous chart's LocalSolution; if MatchIn lies on the
   far side of a singular previous chart, evaluation goes through
   `ApplyCrossing` — exactly once, invariant I6).
3. `w = MatchWeights[...]` (2.7) after `RecombineBasis` (2.8).
4. Assemble the chart's LocalSolution = basis.w + particular (SectorSeries
   algebra), set EpsWindow by honest window arithmetic, set
   `TWindow["CompleteMax"] = expansionOrder - (meta["CouplingDepth"] - 1)`
   (legacy finding 6; old GetLargestTerm discount, LineSegmentation.m:109-113).
5. `SegmentErrorProbe` (2.10); accumulate additively into the running
   ErrorEstimate array; LOUD ERROR E10 when any entry > 1
   (old Transport.m:991-993).
6. Advance.  At a singular chart the marching continues on the far side;
   the side flip is recorded so step 2 of the NEXT chart applies the
   crossing.

Termination: when the plan's last chart is reached, either evaluate at `to`
(regular endpoint; values + error estimates) or, when
`plan["EndpointIsSingular"]`, return the final chart's LocalSolution itself
WITHOUT evaluating (old Transport.m:605-613, 1053-1096; consumed by the FT
layer and by EndpointLimit in API.m).

### 2.7 `MatchWeights[basisValues_, incomingValues_, chartLabel_] -> EpsSeries vector`

The eps-graded Laurent weight solve.  `basisValues` = N x N matrix of
EpsSeries (basis solution i, component c, evaluated at the match point;
columns pre-normalized by eps^(-kmin_i) so every column's window starts at
eps^0 — note homogeneous true-resonance sectors have kmin = -p BY
CONSTRUCTION, RewritePlan 3.1 invariant, and this normalization is exempt
window arithmetic, not an error).  `incomingValues` = EpsSeries vector.

Contract:
- Graded solve: with F~ = Sum_k F[k] eps^k after recombination (2.8),
  assert |det F[0]| > rankTol * scale (Tolerances.m), then solve order by
  order.  ord_eps det F~ != 0 after recombination is LOUD ERROR E5 — never
  a silent window shift (RewritePlan 3.4: matchingShift "target 0,
  asserted"; math review finding 4).
- Result window: honest min over the windows of incomingValues and
  basisValues; no padding.
- Residual assert: per eps order and component,
  |F~.w - v| <= matchTol * scale (Tolerances.m); violation = LOUD ERROR E6
  carrying the residual.  This replaces BOTH the old LinearSolve-with-
  ZeroTest (old Transport.m:341) and the old checked least-squares rescue
  (old Transport.m:348-374) — the rescue is FORBIDDEN (F6).
- Numerically-zero leading rows/pivots: a pivot with
  |pivot| < snapTol * scale is snapped to exact 0; snapTol*scale <= |pivot|
  <= rankTol*scale is the gray zone = LOUD ERROR E7 (no silent rounding;
  RewritePlan section 5: "numerical-zero leading-coefficient skipping
  (generalizes to matching solves)").
- Underdetermined systems: FORBIDDEN inside the marching loop (F7).  The
  old NullSpace free-parameter path (old Transport.m:392-410) survives ONLY
  as the API-level `"?"` wildcard contract (API.md): wildcards enter
  `TransportLine` as symbolic indeterminates already present in `boundary`,
  which keeps the solve square.

### 2.8 `RecombineBasis[basis_List, chart_] -> <|"Basis" -> {LocalSolution..}, "Record" -> ...|>`

Removes the eps=0 degeneracy of the fundamental matrix that regular
incoming data otherwise pays for with 1/eps weights — math review finding 4;
canonical class: log x = (x^(2eps) - 1)/(2eps), i.e. sectors with equal a
(mod Z), different b, whose eps->0 indicial eigenvectors collide.

Contract:
- Candidate families are determined from EXACT data: sector tags with
  a_i - a_j in Z and b_i != b_j, restricted to families flagged eps=0-
  degenerate by Indicial.m metadata (`chart["EpsDegenerateFamilies"]`:
  the eps->0 limit of the residue has a Jordan block that nonzero eps
  resolves).  See OPEN QUESTIONS Q2 for the Indicial.md cross-contract.
- Replacement rule, per degenerate family ordered by b:
  `B_1 = S_1`, `B_m = (S_m - S_1)/((b_m - b_1) eps)` for m >= 2; applied
  recursively to {B_2, ...} while degeneracy persists (depth = size of the
  eps=0 Jordan block).  Each division is an exact EpsSeries window SHIFT.
  Every B_m is a genuine solution (exact linear combination of solutions
  divided by an x-independent scalar).
- Assert after each division: the eps^(-1) (post-shift leading) coefficient
  of B_m's value data is zero after snapping (snapTol); gray zone = E7.
- det effect: each division multiplies det F by 1/((b_m - b_1) eps),
  cancelling exactly one eps order of det degeneracy; "unimodular at eps=0"
  means the FINAL recombined fundamental matrix has det = O(1), assert E5.
- The recombined basis REPLACES the original for all downstream use on this
  chart (assembly, evaluation, probe); the Record field stores the exact
  transformation for debugging/parity dumps.

### 2.9 `CrossingOperator[sector_, sigma_] -> operator` and `ApplyCrossing[ls_LocalSolution, sigma_] -> LocalSolution`

The analytic-continuation rule for evaluating a singular chart's
LocalSolution on the far side (chart coordinate t < 0, u := -t > 0).
sigma in {+1, -1} is the chart's derived Im-sign from its Prescriptions
list (RewritePlan 3.1; consistency rules in 5/E8 below).

Definition (math review finding 5; THE formula, all cases uniform):
```
t^(a + b eps) (eps Log t)^p / p!
  -> e^(sigma I Pi (a + b eps)) *
     Sum[ ((sigma I Pi eps)^j / j!) * u^(a + b eps) (eps Log u)^(p-j) / (p-j)!,
          {j, 0, p} ]
```
i.e. scalar phase e^(sigma I Pi (a + b eps)) TIMES the unipotent
lower-triangular log-chain mixing matrix `M_{p -> p-j} = (sigma I Pi eps)^j / j!`
within the (a,b) family.  Equivalent substitution form:
`Log t -> Log u + sigma I Pi`, `t^(a+b eps) -> e^(sigma I Pi (a+b eps)) u^(a+b eps)`.
Non-integer a matters: the phase carries `a`, not just `b eps`
(e.g. a = 1/2 gives a factor sigma*I).

Convention pinning (MUST match the frozen oracle,
DiffExp/AnalyticContinuation.m:18-90, replacement rules at :70-79):
- sigma = +1 corresponds to the old code applying NO replacement
  (`AnalyticContinuationReplacements = {}` unless uniqueSigns === {-1},
  AnalyticContinuation.m:70-79) and then evaluating `Logx -> Log[at]` and
  powers at NEGATIVE numeric arguments under Wolfram's PRINCIPAL branch
  (Pade.m:70-77 SEval2).  Principal branch at negative real = the
  Im > 0 side: `Log(-u) = Log u + I Pi`, `(-u)^q = e^(I Pi q) u^q`.
  The formula above with sigma = +1 reproduces this exactly.
- sigma = -1 corresponds to the old full-rotation correction on the theta_m
  (far) side: `Logx -> Logx - 2 Pi I theta_m` (a SHIFT that mixes log-chain
  members binomially — not a phase), `x^b -> e^(-2 Pi I b) x^b` applied only
  for Denominator[b] > 2, and the half-integer case routed through
  `x^(b-1/2) (theta_p - theta_m) Sqrt[x]` (AnalyticContinuation.m:72-77).
  Composed with the principal branch these all equal the single formula
  with sigma = -1; unit test T8 pins the equality case by case.
- The old theta_p/theta_m two-sided representation and
  `Project\[Theta]s` (AnalyticContinuation.m:93-103) are SUBSUMED: DiffExp2
  keeps one-sided data plus the explicit operator.
- 9aeb300 interior-split lesson (commit message + current
  DiffExp/RegularizedIntegration.m:523, :1842): imaginary parts produced by
  a crossing must be applied EXACTLY ONCE, by this operator, from the
  prescription — never by accidental principal-branch evaluation of Log at
  a negative number, and PV-paired half-segments (Integrate.m) use real
  logs of the distance with the phases accounted at the object level.
  Consequences for Transport.m: (i) SectorSeries evaluation is only ever
  called with POSITIVE chart-coordinate arguments; `ApplyCrossing` happens
  first (invariant I7); (ii) `CrossingOperator` is the single exported
  authority that Integrate.m reuses for its interior-pole pairing.

Family closure: applying the operator to a p > 0 sector POPULATES ALL
lower-p members of the same (a,b) family.  Target sectors that do not yet
exist in the LocalSolution are CREATED; dropping any contribution is
forbidden (F8).  Window arithmetic: the eps^j factor shifts that
contribution's window up by j; the merged sector window is the honest min
(EpsSeries contract).

### 2.10 `SegmentErrorProbe[ls_LocalSolution, {tIn_, tOut_}, indeterminates_List] -> errs[[component, epsorder]]`

The two-point full-vs-reduced-order probe, ported from old
Transport.m:905-993:
- Evaluate ls at tOut and (when tIn != 0) at tIn, twice: at full order and
  with the t-series truncated down by
  `probeDecrement = Ceiling[0.7 * couplingDepth] + 2`
  (old `ICurrEvalErrorSeriesDecrease`, State.m:222; reduction applied as in
  old Transport.m:913 via DecreaseSeriesOrderBy).  Both evaluations go
  through SectorSeries (Pade applied there when configured, with ITS loud
  fallback — Transport must not Quiet it, F5).
- tIn == 0 (chart center) is SKIPPED — at the center every log-bearing term
  vanishes spuriously and the error would read 0 (old comment and hack,
  Transport.m:915-923); the tOut value is used for both points in that case
  (parity with old :920-922).
- Per-indeterminate errors: for each symbolic indeterminate, the error of
  its coefficient is tracked separately, plus the constant part
  (old ComputeErrorsPerIndeterminate, Transport.m:947-962).
- The per-(component, epsorder) error = max over the two probe points of
  |full - reduced| (old Transport.m:970-976).
- Accumulation across segments is ADDITIVE into
  `TransportResult["ErrorEstimate"]` (old Transport.m:980-984); abort > 1 is
  E10.  Seeding: when `boundary` carries no error data the seed is exact 0
  with `"ErrorSeeded" -> False` metadata — the old Accuracy[]-based seeding
  and PadRight zero-padding (Transport.m:583-591) are FORBIDDEN (F11).
- Per-segment digit check: when AccuracyGoal is numeric, the per-segment
  error must satisfy err <= 10^-DigitsNeeded; violation = LOUD ERROR E11
  naming the segment and the order increase that would be needed.  (This
  replaces the old adaptive AccuracyGoalValidate Before/After machinery —
  see section 9 cuts and Q4.)

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution / Sector / EpsWindow (RewritePlan 3.1, VERBATIM)

```
LocalSolution = <|
  "Center" -> exact x0, "ChartMap" -> affine (Mobius optional, see 3.2),
  "Radius" -> distance to nearest singularity IN THE COMPLEX PLANE
              (complex singularities are real: pentagon/unequal-mass
              lines have them; old code projects ghosts Re, Re±Im —
              ledger item; new code uses true complex distance),
  "Sectors" -> { Sector.. },
  "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "TWindow"   -> <|"CompleteMax" -> nmax|>,   (* t-order honesty: top
              (couplingDepth−1) orders of chained particular solutions
              are degraded (old MaxCouplingOrder/ISafetyExpansionSubtract
              lesson); tracked, not guessed *)
  "ErrorEstimate" -> per (eps-order) accumulated error (two-point
              full-vs-reduced-order probe, additive across segments,
              abort > 1 — ported old machinery, user-facing),
  "Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
              "LeadingCoeffSign"|> with the derived chart Im-sign;
              consistency-checked at construction (even multiplicity = no
              constraint; conflict or missing prescription at a chart with
              b!=0/p>0 sectors = LOUD ERROR; sqrt factors auto-prescribed
              as in old State.m DEqnSquareRoots)
|>

Sector = <|
  "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
  "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT:
              solutions are vectors; matching and ODE residuals need all
              components *)
|>
Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                Σ_k eps^k Σ_n c[k,n,comp] t^n
```

Window exemption (RewritePlan 3.1 invariants): homogeneous true-resonance
sectors have kmin = -p BY CONSTRUCTION; Transport's column normalization
(2.7) and window arithmetic must account for this, not error.

### 3.2 Chart (produced by SegmentLine, consumed by TransportLine/Solve.m)

```
Chart = <|
  "Index"      -> k (1-based position in the plan),
  "Center"     -> exact x0,
  "Map"        -> <|"Type" -> "Affine" | "Mobius",
                    (* Affine: x = Center + t * Scale, with
                       Scale = trueRadius / RadiusOfConvergence so the
                       chart-coordinate radius is RoC and high-order
                       coefficients stay O(1) — old Mobius.m:47,65;
                       banana REQUIRES RoC = 10 (ledger).
                       Mobius: the old GetMobius map, Mobius.m:22-35,
                       composed with t -> t/RoC; optional, see Q3. *)
                    "Scale" -> ..., "Params" -> ...|>,
  "Radius"     -> chart-coordinate radius (= RoC for affine charts; for
                  Mobius charts the min |chart image| over ALL
                  singularities, complex included),
  "IsSingular" -> True | False,
  "MatchIn"    -> chart-coordinate t of the incoming match point (None for
                  the first chart; |MatchIn| = Radius/DivisionOrder up to
                  the FixWithin midpoint clipping, 2.3),
  "MatchOut"   -> chart-coordinate t of the outgoing match point, or
                  "Endpoint",
  "CrossSign"  -> +1 | -1 | None  (derived chart Im-sign; None iff not
                  crossed or no multivalued content),
  "Prescriptions" -> as in 3.1,
  "EpsDegenerateFamilies" -> from Indicial.m (see 2.8, Q2),
  "Label"      -> string "chart #k at x0=<N[center]> (singular|regular)"
                  used verbatim in every error message
|>
```

### 3.3 SegmentPlan

```
SegmentPlan = <|
  "From" -> exact, "To" -> exact, "Direction" -> +1 | -1,
  "Charts" -> {Chart..}, "SegmentCount" -> n,
  "DigitsNeeded" -> Integer, "EndpointIsSingular" -> True | False
|>
```

### 3.4 TransportResult

```
TransportResult = <|
  "Values"   -> vals[[integral, epsorder]]      (* regular endpoint *)
   | "Solution" -> LocalSolution,               (* singular endpoint *)
  "ErrorEstimate" -> errs[[integral, epsorder]] (* additive, per 2.10 *),
  "ErrorSeeded" -> True | False,
  "EpsWindow" -> as 3.1, "TWindow" -> as 3.1,
  "SegmentCount" -> n,
  "SavedCharts" -> {(* per-chart LocalSolutions, only when the
                      SaveExpansions-equivalent config is on *)} | Missing
|>
```
API.m wraps this into the user-facing TransportTo association (kinematic
point, NumIntegrals, orders — API.md contract; old return shape at old
Transport.m:1234-1246).

### 3.5 Weights

`MatchWeights` returns an EpsSeries vector with an explicit window
`<|"Min" -> 0 (asserted for regular incoming data), "CompleteMax" -> ...|>`;
any Min < 0 for data whose own window starts at 0 is E5 territory (it means
recombination failed).

---

## 4. INVARIANTS (always-on cheap asserts)

I1  Plan fidelity: TransportLine executes exactly plan["SegmentCount"]
    segments with exactly the planned centers/match points; any deviation
    aborts (kills the old duplicated-dry-run drift, old Transport.m:666).
I2  Geometry: for every interior match point m between charts i and i+1,
    |t_i(m)| / Radius_i == 1/DivisionOrder == |t_{i+1}(m)| / Radius_{i+1}
    within snapTol (the GetCPL/GetCPR guarantee, Mobius.m:98-142); every
    evaluation point satisfies |t| <= Radius/DivisionOrder < Radius.
I3  ord_eps det F~ == 0 after recombination, and |det F~[0]| > rankTol*scale
    (per 2.7/2.8).
I4  Matching residual <= matchTol*scale per eps order and component.
I5  Windows are monotone: no operation increases EpsWindow["CompleteMax"]
    or TWindow["CompleteMax"]; no padding ever (kmin = -p exemption per
    3.1).
I6  Crossing is applied exactly once per traversed singular chart (a
    per-chart "Crossed" flag checked at evaluation time).
I7  SectorSeries evaluation is only invoked with strictly positive chart
    coordinates; far-side requests must already be crossed (9aeb300
    lesson).
I8  Precision: after input raising (section 6, L1), every numeric input has
    Accuracy >= 2*WP at entry, and all evaluation blocks run under
    $MinPrecision = WP (old Pade.m:34,70,80) with $MaxExtraPrecision
    headroom (old Mobius.m:39).
I9  ErrorEstimate entries are monotone nondecreasing across segments.
I10 Regular chart <=> exactly one sector (a=0, b=0, p=0) among the
    HOMOGENEOUS solutions (RewritePlan 3.1, scope per math finding 8:
    ODE solutions only).
I11 ODE residual spot-check at one random interior point per chart
    (RewritePlan 3.1 invariant; executed on the assembled chart
    LocalSolution, residual <= matchTol*scale).

---

## 5. ERROR CONTRACT

Every error below is LOUD (a thrown, aborting error in the DiffExp2 error
idiom) and its message MUST carry the Chart "Label" (chart index + center +
singular/regular) where applicable, plus the named fields listed.  There is
NO error class in this module that downgrades to a warning, with the single
structured exception noted in E8.

E1  UNSOLVABLE SINGULAR FACTOR.  A matrix-denominator or extra singular
    factor has no exact root solution in the line parameter.
    Carries: the factor, the variable.  (Old code Quiets Solve,
    LineSegmentation.m:76.)
E2  LINE ON SINGULARITY.  A singular factor vanishes identically along the
    line.  Carries: the vanishing factors.  (Port of old
    Transport.m:571-574.)
E3  DIGIT BUDGET EXCEEDS PRECISION.  DigitsNeeded + chopReserve > WP.
    Carries: DigitsNeeded, segment count, WP, the minimal WP that would
    suffice.  (NEW — old code has no check.)
E4  SOLVE/INDICIAL FAILURE PROPAGATION.  Errors raised by Indicial.m or
    Solve.m (indicial contract violation, non-Fuchsian beyond rank
    reduction, resonance machinery) pass through UNCAUGHT.  Transport never
    wraps Solve calls in Check/Quiet (F12).
E5  SINGULAR FUNDAMENTAL MATRIX AT eps=0.  ord_eps det F~ > 0 after
    recombination, or |det F~[0]| <= rankTol*scale.
    Carries: chart label, the full sector tag list {a_i, b_i, p_i},
    the measured ord_eps det F, the match point, which recombinations were
    already applied.  NO silent EpsWindow shift (RewritePlan 3.4
    matchingShift assert).
E6  MATCHING RESIDUAL.  |F~.w - v| > matchTol*scale at some (eps order k,
    component c).  Carries: chart label, k, c, residual, scale.
    (Replaces old solveFailed -> least-squares -> ReportError chain,
    Transport.m:337-387.)
E7  TOLERANCE GRAY ZONE.  A pivot/leading coefficient lies between
    snapTol*scale and rankTol*scale (matching solve or recombination
    post-shift leading coefficient).  Carries: chart label, the value, both
    thresholds, the operation.
E8  PRESCRIPTION CONFLICT OR MISSING.  At a singular chart whose solutions
    contain multivalued sectors (b != 0, p > 0, or non-integer a), the
    Prescriptions list yields no consistent derived Im-sign (conflicting
    signs across simultaneously vanishing factors after the
    multiplicity-parity and leading-coefficient-sign reduction —
    derivation rules per old AnalyticContinuation.m:45-68: even vanishing
    multiplicity = no constraint; required sign = prescribed sign divided
    by leading-coefficient sign) or no prescription covers the chart at
    all.  Carries: chart label, every vanishing factor with its prescribed
    sign/multiplicity/leading-coefficient sign, the per-factor derived
    Im-signs, and the hint to add/fix DeltaPrescriptions.
    The old behavior of flagging AnalyticContinuationFailed and CONTINUING
    with sign +1 (AnalyticContinuation.m:66-67) is FORBIDDEN (F1).  The old
    one-sided-validity warning path (Transport.m:494-496) and the
    MultivaluedFail return-last-segment path (Transport.m:899-943,
    1185-1188, 1225-1227) are FORBIDDEN (F2).  Structured exception: iff
    config AbortOnAnalyticContinuationFail is explicitly False (the FT
    setting, FeynmanTrick/DiffExpIntegration.m:411-414) AND the chart is
    the FINAL one in the plan, TransportLine may return a TransportResult
    whose completeness metadata marks every order INCOMPLETE from that
    chart on — an honest partial result, never a silently one-sided value.
    Anywhere else: error.  (See Q5.)
E9  CROSSING WITHOUT SIGN.  ApplyCrossing invoked with CrossSign None on a
    chart with multivalued content (should be unreachable given E8; kept as
    a cheap assert).  Carries: chart label, sector tags.
E10 ERROR ESTIMATE > 1.  Any accumulated ErrorEstimate entry exceeds 1.
    Carries: chart label of the segment that tipped it, the (integral,
    eps-order) entry, the per-segment and accumulated values.  (Port of old
    Transport.m:991-993.)
E11 SEGMENT ACCURACY MISS.  Numeric AccuracyGoal configured and a segment's
    probe error > 10^-DigitsNeeded.  Carries: chart label, measured error,
    DigitsNeeded, current ExpansionOrder, suggested ExpansionOrder.
E12 WINDOW UNDERFLOW.  A consumer (matching, crossing, evaluation,
    assembly) requests an eps order > EpsWindow["CompleteMax"] or a t order
    > TWindow["CompleteMax"] of any operand.  Carries: chart label, sector
    tags, requested order, available CompleteMax.  (RewritePlan 3.1: "fail
    loudly naming (chart, sector, order)".)
E13 NO MARCHING PROGRESS.  A planned step does not strictly advance toward
    `to`, or executed segments exceed plan["SegmentCount"] (I1).  Carries:
    chart label, match points.  (Old code can loop forever; NEW.)
E14 AMBIGUOUS ENDPOINT.  `to` is within snapTol of a singularity but not
    exactly equal to it (exact arithmetic distinguishes; chop-based
    detection of the old Transport.m:606-609 is forbidden, F13).  Carries:
    `to`, the nearby root, their exact difference.

### Forbidden fallbacks (exhaustive; each is a place a fallback is tempting)

F1  Prescription conflict -> proceed with sign +1 and a flag
    (old AnalyticContinuation.m:55-68).  FORBIDDEN -> E8.
F2  Multivalued-without-prescription -> warn, abort early, return last
    segment as the result (old Transport.m:470-478, 899-943, 1185-1188).
    FORBIDDEN -> E8 (structured partial result only under the narrow E8
    exception).
F3  Numeric chop-deduplication of singularity roots and chop-membership
    endpoint tests (old LineSegmentation.m:77, Transport.m:606-609).
    FORBIDDEN -> exact arithmetic + E14.
F4  Ghost projection of complex roots to Re, Re±Im real points
    (old LineSegmentation.m:81-99, suppression rule included).  FORBIDDEN
    -> true complex distance in ChartRadius (the suppression lesson — do
    not create extra charts for complex roots — is preserved by
    construction: complex roots only cap radii).
F5  Pade failure -> silently evaluate the plain partial sum
    (old Pade.m:44-46 warns at verbosity and proceeds).  Pade lives in
    SectorSeries (RewritePlan 3.2); Transport MUST NOT Quiet/Check around
    SectorSeries evaluation.
F6  Singular matching solve -> residual-checked LeastSquares rescue
    (old Transport.m:348-374).  FORBIDDEN -> E5/E6.  The rescue compensated
    representation loss that the sector-native design removes; if the solve
    is singular, the basis or the windows are wrong and the user must see
    it.
F7  Underdetermined matching -> introduce free parameters via NullSpace
    (old Transport.m:392-410).  FORBIDDEN inside marching -> E6; wildcard
    `"?"` boundary data is an API.m-level contract that arrives as
    symbolic indeterminates.
F8  Crossing a p > 0 sector -> apply only the scalar phase and drop the
    lower-p contributions (what a naive port of the old per-eps-order rules
    would do).  FORBIDDEN: the unipotent mixing is mandatory; absent target
    sectors are created.
F9  Maintaining two copies of the segmentation logic (dry run vs marching,
    old Transport.m:666-750 vs :999-1220).  FORBIDDEN: one SegmentLine, one
    plan object, I1 asserts.
F10 Silent EpsWindow shift to absorb ord_eps det F > 0 ("account the shift
    and move on").  FORBIDDEN -> E5.  RewritePlan 3.4 allows a matchingShift
    term ONLY as an asserted-zero quantity.
F11 Seeding/padding error arrays from input Accuracy or PadRight
    (old Transport.m:583-591).  FORBIDDEN -> exact 0 seed +
    "ErrorSeeded" -> False metadata.
F12 Catching errors from Solve.m/Indicial.m/SectorSeries (Check/Quiet) to
    keep marching.  FORBIDDEN -> E4.
F13 Chop-based "endpoint is a singularity" detection (old
    Transport.m:606-609).  FORBIDDEN -> exact membership + E14.
F14 Skipping the error probe when an evaluation fails or when in a hurry
    (EstimateError must be explicitly disabled in config to skip; when
    disabled, ErrorEstimate -> Missing["NotComputed"], never zeros — old
    Transport.m:1238-1239 parity).
F15 Falling back from the exact `FindSingularities` root set to a numeric
    root finder when Solve is slow.  FORBIDDEN; performance issues are
    R2-class (flagged, not worked around).

---

## 6. ABSORBED OLD CODE

### 6.1 Replaced functionality (file:line)

| Old code | What it did | Disposition |
|---|---|---|
| DiffExp/Transport.m:514-1247 `TransportTo` | marching loop, segmentation driving, error accumulation, return assoc | core of `TransportLine`/`SegmentLine` |
| DiffExp/Transport.m:666-750 | predivision dry-run segment counting (duplicated code) | `SegmentLine` (single implementation, F9) |
| DiffExp/Transport.m:758-760 | DigitsNeeded = AccuracyGoal + Ceiling[Log10[#segments]] + ISafetyDigits | `DigitBudget` (2.5) + E3 |
| DiffExp/Transport.m:541-562 (also :38-41, :170-176) | input precision raising to 2x WP with warning | kept; I8/L1 |
| DiffExp/Transport.m:571-574 | line-on-singularity error | E2 |
| DiffExp/Transport.m:603, :615-626 | FindMatrixSingularities call + per-singularity chart construction | `FindSingularities` + `SegmentLine` |
| DiffExp/Transport.m:605-613, 1053-1096 | singular-endpoint mode (integrate final segment, return series) | `TransportLine` termination, 2.6 |
| DiffExp/Transport.m:635-655, 679, 1041-1046 | predivision pole intervals at ±RoC/DivisionOrder | match-point placement, 2.3 |
| DiffExp/Transport.m:1124-1136 | FixWithin midpoint clipping before a singular chart | preserved, 2.3 |
| DiffExp/Transport.m:1147-1168 | regular-step next center via FindNextCenterPointL/R | `NextCenter` (2.4) |
| DiffExp/Transport.m:776-841, 1191-1213 | AccuracyGoalValidate "Before" adaptive order search / "After" redo | NOT ported; replaced by E11 (see 9, Q4) |
| DiffExp/Transport.m:287-441 (in IntegrateSystem) | boundary fixing: equations, LinearSolve, least-squares rescue, NullSpace free params | `MatchWeights` (2.7) + E5/E6/E7; rescue/NullSpace forbidden F6/F7 |
| DiffExp/Transport.m:317-327 | FixAt != 0 coefficient extraction by differentiation through nested SeriesData | obsolete: exact tags + SectorSeries evaluation |
| DiffExp/Transport.m:905-993 | SEval1/SEval2 evaluation + two-point error probe + per-indeterminate errors + abort > 1 | `SegmentErrorProbe` (2.10) + E10 |
| DiffExp/Transport.m:947-962 | ComputeErrorsPerIndeterminate | inside 2.10 |
| DiffExp/Transport.m:1234-1246 | return association | TransportResult (3.4) via API.m |
| DiffExp/LineSegmentation.m:65-106 `FindMatrixSingularities` | factor collection, Quiet Solve, chop dedup, ghost projection | `FindSingularities` (2.1); F3/F4 |
| DiffExp/LineSegmentation.m:109-113 `GetLargestTerm` | matrix-truncation magnitude at coupling-discounted order | subsumed by TWindow discount + probe (2.6 step 4, 2.10) |
| DiffExp/LineSegmentation.m:116-145 `GetMatricesPrecisionDistance` | Dynamic segmentation intervals | DROPPED (SegmentationStrategy = Predivision only in v1, RewritePlan Config.m) |
| DiffExp/LineSegmentation.m:20-62, 239-258 `RelateLines`/`GetMatchingPoint` | Solve-based relation between segment parametrizations | obsolete: explicit invertible ChartMaps composed exactly |
| DiffExp/LineSegmentation.m:148-236 `CheckBoundaryConditionsAndReparametrize` | bcs validation/reparametrization | API.m's contract; Transport asserts normalized input |
| DiffExp/Mobius.m:22-35 `GetMobius` | Mobius chart maps | Chart "Map" Type "Mobius" (optional, Q3) |
| DiffExp/Mobius.m:38-69 `GetLineRescaled` | chart construction + RoC rescaling + 2x WP SetPrecision | Chart construction in SegmentLine; L1/L4 |
| DiffExp/Mobius.m:72-142 GetMobiusCPL/R, GetCPL/GetCPR | next-center geometry (1/k of both charts) | `NextCenter` (2.4) |
| DiffExp/Mobius.m:145-155 FindNextCenterPointL/R | bound selection wrapper | inside `NextCenter` |
| DiffExp/AnalyticContinuation.m:18-90 `PrepareAnalyticContinuation` | per-chart vanishing-factor collection, sign derivation, replacement rules, conflict flag | sign derivation -> Prescriptions construction (3.1, shared contract with SectorSeries/API); rules -> `CrossingOperator` (2.9); conflict -> E8 |
| DiffExp/AnalyticContinuation.m:93-103 `Project\[Theta]s` | theta-projection of two-sided data | subsumed by one-sided data + explicit crossing |
| DiffExp/Pade.m:55-91 SEval1/SEval2/SEval | evaluation with continuation + theta resolution | evaluation -> SectorSeries (with Pade, RewritePlan 3.2); theta resolution -> ApplyCrossing |
| DiffExp/MatrixLoading.m:357-385 `InitializeIntegrationSequence` | coupled-block ordering + MaxCouplingOrder | block ordering -> Solve.m; depth metadata consumed here for TWindow/probe |
| DiffExp/MatrixLoading.m:217-231 MatricesIrreducibleFactors | singular-factor inventory | provided by LoadSystem (API.md); consumed by 2.1 |
| FeynmanTrick/DiffExpIntegration.m:204-238 appendMatrixFactors / deltaPrescriptionsForFactors | merging IBP singular factors + endpoint prescriptions into DiffExp state | first-class config inputs: sys["ExtraSingularFactors"], Prescriptions |
| DiffExp/Transport.m:277-282 segmentCaches/blockCache | per-segment block caching across eps orders | largely moot under symbolic-eps solve; per-chart Solve results are computed once per chart by construction (L9) |

### 6.2 Numerical lessons that MUST be preserved (read from the old code)

L1  2x WP input raising + precision floors.  Inputs (line, endpoint, bcs)
    raised to 2*WorkingPrecision when Accuracy < WP, with a visible warning
    (old Transport.m:541-562, :38-41); chart maps built under
    SetPrecision[..., 2 WP] (Mobius.m:48,66) with $MaxExtraPrecision = 1000
    (Mobius.m:39); evaluations under $MinPrecision = WP (Pade.m:34,70,80).
    The 2x headroom is what survives chart-map composition (legacy review
    finding 18).
L2  Segment-count digit budgeting.  DigitsNeeded = AccuracyGoal +
    Ceiling[Log10[#segments]] + 2, computed AFTER a geometry-only counting
    pass (old Transport.m:666-760; ISafetyDigits = 2, State.m:215).
L3  Both-adjacent-charts match-point geometry.  The match point must sit at
    radius/DivisionOrder of BOTH the chart that produced the values and the
    chart that consumes them (GetCPL/GetCPR solve for this, Mobius.m:98-142)
    — this, not the solver, is the conditioning defense for matching (R3).
    Classic default k=3 (State.m:112); FT pins k=4
    (DiffExpIntegration.m:272).
L4  RoC chart rescaling.  Chart coordinate scaled by trueRadius/RoC so that
    the series' convergence radius maps to RoC and high-order coefficients
    stay O(1) (Mobius.m:47,65); banana classic line REQUIRES RoC = 10
    (Reference/Examples Banana_example.m; RewritePlan section 5).
L5  Complex singularities are real and cap radii.  Old ghost projection
    with suppression (LineSegmentation.m:81-99) exists BECAUSE complex
    roots limit convergence; the new true-complex-distance radius preserves
    the cap while deleting the ghosts (legacy finding 4; pentagon and
    unequal-mass banana have complex roots on generic lines).
L6  Two-point error probe.  Full vs reduced order at BOTH the incoming
    match point and the evaluation point, skipping the chart center where
    log terms vanish spuriously (old Transport.m:905-925, hack at
    :916-923); per-indeterminate error split (:947-962); additive
    accumulation (:980-984); hard abort when > 1 (:991-993).
L7  Coupling-depth t-order degradation.  The top (couplingDepth - 1)
    t-orders of chained particular solutions are garbage; the old code
    discounts them in the matrix-error read-off
    (ExpansionOrder - ISafetyExpansionSubtract - (MaxCouplingOrder - 1),
    LineSegmentation.m:109-113, ISafetyExpansionSubtract = 5 State.m:216)
    and in the probe decrement (Ceiling[0.7 MaxCouplingOrder] + 2,
    State.m:222).  New: TWindow["CompleteMax"] = N - (couplingDepth - 1),
    probe decrement formula preserved verbatim.
L8  Crossing convention.  sign +1 = NO replacement, relying on
    principal-branch Log/Power on the far side; sign -1 = Logx ->
    Logx - 2 Pi I theta_m SHIFT (binomial log-chain mixing), e^(-2 Pi I b)
    only for Denominator[b] > 2, half-integers through explicit Sqrt with
    (theta_p - theta_m) (AnalyticContinuation.m:70-79).  The single formula
    of 2.9 must reproduce all of these (unit test T8).  9aeb300 lesson:
    crossing phases come from the prescription exactly once; PV-paired
    interior splits use real logs ($InteriorSplitRealLog,
    RegularizedIntegration.m:523, :1842).
L9  Per-chart work is computed once.  Old code caches solver blocks across
    eps orders within a segment (Transport.m:277-282) and drops matrix
    expansions per segment (ClearMatrices, Transport.m:1215-1218,
    MatrixLoading.m:235-249); new code solves each chart exactly once
    (symbolic eps) and must release per-chart data after its outgoing match
    (memory discipline on long lines).
L10 Gray-zone honesty in numerical-zero decisions.  Leading-coefficient /
    pivot zero tests must be relative and three-valued (zero / nonzero /
    loud gray zone) — generalization of the FT LaurentTrim lesson
    (RewritePlan section 5; Tolerances.md snapTol/rankTol).

---

## 7. DEPENDENCIES

Order (acyclic, RewritePlan):
`Tolerances < Config < EpsSeries < SectorSeries < Indicial < Solve <
Transport/Integrate < API`.

Transport.m MAY call:
- Tolerances.m: matchTol, snapTol, rankTol, safetyDigits, chopReserve,
  probe constants.
- Config.m: WorkingPrecision, AccuracyGoal(+Validate), DivisionOrder,
  RadiusOfConvergence, ExpansionOrder, EpsilonOrder, EstimateError,
  SegmentationStrategy (must be "Predivision"; anything else is Config's
  loud error), UseMobius, UsePade (read-only; enforcement in SectorSeries),
  DeltaPrescriptions, AbortOnAnalyticContinuationFail, Verbosity,
  SaveExpansions-equivalents.
- EpsSeries.m: all window-carrying eps-Laurent arithmetic, the graded/
  Laurent linear solve used by MatchWeights, exact eps-division (window
  shift) used by RecombineBasis.
- SectorSeries.m: evaluate (with Pade and ITS loud fallback), multiply,
  re-expand, assemble/add LocalSolutions, differentiate.
- Indicial.m: per-chart sector specs, Prescriptions derivation inputs,
  EpsDegenerateFamilies metadata (Q2).
- Solve.m: SolveChart (fundamental basis + particular + CouplingDepth and
  related metadata).

Transport.m MUST NOT call: Integrate.m (sibling; Integrate may call
Transport's exported CrossingOperator), API.m, anything in FeynmanTrick/,
anything in the frozen DiffExp/ tree.

---

## 8. UNIT TESTS

Closed-form/pure tests run without the old library; tests marked (kernel,
oracle) consume M0-pregenerated artifacts and belong to the M4 gate but are
enumerated here because they pin THIS module's behavior.

T1  `test_segmentation_exact_roots`
    Factors {1 - 4 x, x, 1 - x}; SegmentLine from 11/23 to 0.  Asserts: the
    singular roots are EXACTLY {0, 1/4, 1} (SameQ); the plan contains a
    singular chart with Center === 1/4; plan["EndpointIsSingular"] ===
    True (endpoint 0).  (Box L1 geometry; campaign interior pole at 1/4.)
T2  `test_segmentation_complex_radius`
    Factors {x^2 + 1, x - 3}.  Asserts: ChartRadius[0, all] == 1 — the
    COMPLEX roots ±I are binding (real root 3 is farther), and there are
    no ghost projections at ±1; ChartRadius[2, all] == Min[|2 - 3|,
    |2 - I|] == Min[1, Sqrt[5]] == 1 — the real root is binding; no chart
    in any plan is CENTERED at a non-real point or at a ghost projection.
T3  `test_segmentation_exact_complex_roots`
    Factor x^2 - 2 x + 2 (roots 1 ± I).  Asserts ChartRadius[1, all] == 1
    with the exact algebraic roots retained in "All".
T4  `test_nextcenter_both_charts`
    Singularities {0, 1}, previous match point 1/2, k = 4, right-moving.
    Asserts NextCenter == 3/5 exactly, new radius == 2/5, and
    |1/2 - 3/5| == (2/5)/4 == 1/10 (the GetCPL closed form,
    Mobius.m:110-119).  Property sweep: for 50 random rational
    {z_min, x_b, z_max, k in 2..6}, the 1/k relation holds for BOTH
    adjacent charts within snapTol (I2).
T5  `test_digit_budget`
    DigitBudget[20, 12] == 20 + Ceiling[Log10[12]] + 2 == 24.  And: E3
    fires when DigitsNeeded + chopReserve > WP (construct WP = 30 case);
    message names the minimal sufficient WP.
T6  `test_plan_fidelity`
    A 5-segment toy plan; TransportLine result has SegmentCount == 5 and
    the executed centers equal the plan's (I1).  Mutating the plan
    mid-flight (test hook) triggers E13.
T7  `test_precision_raise`
    Endpoint given at MachinePrecision; asserts the visible warning fires
    and internal chart maps carry Accuracy >= 2*WP (I8, L1).
T8  `test_crossing_convention_parity`
    For sigma in {+1, -1} and (a, b, p) in
    {(1/2, 0, 0), (1/3, 0, 0), (1, 0, 0), (0, 1/3, 0), (0, 2, 0)}:
    CrossingOperator's scalar equals the OLD-rule composition —
    sigma=+1: principal branch, e.g. a=1/2 -> +I (since
    (-u)^(1/2) = I Sqrt[u]); sigma=-1: e^(-2 Pi I a) * principal, e.g.
    a=1/2 -> -I, a=1/3 -> e^(-I Pi/3); b-only sectors: e^(sigma I Pi b eps)
    expanded in eps matches substituting Log t -> Log u + sigma I Pi in the
    eps-expansion of t^(b eps).  Pins AnalyticContinuation.m:70-79 (L8).
T9  `test_crossing_log_chain_mixing`
    p = 2 family, sigma = -1: ApplyCrossing of (eps Log t)^2/2! yields
    contributions { (eps Log u)^2/2!, (-I Pi eps)(eps Log u),
    (-I Pi eps)^2/2! } to p = {2, 1, 0}, each times e^(-I Pi(a + b eps)).
    Verified against direct substitution Log t = Log u - I Pi.
T10 `test_crossing_populates_lower_p`
    LocalSolution with a single (0, 0, 2) sector; after ApplyCrossing the
    (0,0,1) and (0,0,0) sectors EXIST with the T9 coefficients; a check
    that no contribution was dropped: value parity at u = 1/7 between
    direct substitution and the crossed object to 10^-(WP/2).
T11 `test_matching_regular`
    One-sector regular chart, 2x2 toy basis with F[0] = IdentityMatrix;
    MatchWeights returns w == v exactly; residual 0; weight window Min 0.
T12 `test_matching_recombination_logx`
    System A = {{0, 1/x}, {0, 2 eps/x}} (the log x = (x^(2eps)-1)/(2eps)
    class).  Naive basis S1 = (1, 0), S2 = (x^(2eps), 2 eps x^(2eps)):
    ord_eps det F == 1 at any match point.  Asserts: RecombineBasis
    produces B2 = (S2 - S1)/(2 eps); det F~[0] == 1 (unimodular); for
    incoming data = the regular solution ((x^(2eps)-1)/(2 eps),
    x^(2eps)-1) evaluated at t_m = 1/3, MatchWeights == (0, 1) to matchTol,
    with weight window starting at eps^0 (NOT eps^-1).
T13 `test_matching_det_assert_loud`
    Two identical basis columns with IDENTICAL tags (not recombinable):
    E5 fires; message contains the chart label string, all sector tags,
    and ord_eps det F.
T14 `test_matching_residual_loud`
    Perturb T11's v by 10^(-AG/2) in one component at eps order 1: E6
    fires naming (eps order 1, that component, the residual).
T15 `test_gray_zone_loud`
    Pivot constructed at geometric mean of snapTol and rankTol scales: E7
    fires naming both thresholds.
T16 `test_error_probe_closed_form`
    ls = single regular sector with c[0, n] = 1 for n <= 10 (f = sum t^n),
    probeDecrement = 3, probe points {0, 1/3}.  Asserts: probe SKIPS t = 0
    and uses the t = 1/3 value for both points (L6 parity); error ==
    Sum[(1/3)^n, {n, 8, 10}] == 13/59049 exactly (relative to evaluation
    precision).  With boundary data v + c*w (symbolic c): separate error
    entries for the c-coefficient and the constant part.
T17 `test_error_abort`
    Inject a segment error of 2 (test hook on the probe): E10 fires naming
    the segment and the (integral, eps-order) entry.
T18 `test_twindow_coupling`
    Solve metadata CouplingDepth = 3 at ExpansionOrder 20:
    TWindow["CompleteMax"] == 18; probe decrement == Ceiling[0.7*3] + 2
    == 5 (L7; State.m:222 formula).
T19 `test_singular_endpoint_series`
    1x1 system A = b eps / x with exact b = 2; transport from 1/2 (boundary
    value = the eps-expansion of (1/2)^(2 eps) to EpsilonOrder) to the
    singular endpoint 0.  Asserts: TransportResult has "Solution" (no
    "Values"); the LocalSolution has exactly one sector (a, b, p) ==
    (0, 2, 0) with coefficient series == 1 (eps^0) and 0 at all higher
    orders, complete to the requested EpsWindow.
T20 `test_bidirectional_roundtrip`
    Same 1x1 system; transport 1/2 -> 1/5 -> 1/2 (regular points only):
    returned values match the start values to 10^-AG; the two plans have
    Direction -1 and +1 respectively.
T21 `test_no_ghost_charts_near_endpoint`
    Singularities {0, 12167/12651, 1} (box_triangle L3 geometry, the
    0.9617 interior pole): plan from 11/23 to 1 keeps 12167/12651 and 1 as
    DISTINCT exact centers (F3), all match points satisfy I2, and the plan
    terminates (E13 never fires).
T22 (kernel, oracle) `test_marching_parity_bubble`
    Classic bubble line vs M0 checkpoint dumps
    (Scripts/dump_transport_checkpoints.m output): values at every segment
    boundary agree to 1e-25 under the example-pinned old config
    (RewritePlan M4).
T23 (kernel, oracle) `test_interior_crossing_box_quarter`
    Box L1 line crossing the pole at exactly 1/4 with prescription sign +1:
    far-side values match the campaign pin (RewritePlan M4 "interior-pole
    crossing (box L1 pole at 1/4, campaign pin)"); imaginary parts appear
    iff the prescription says so (9aeb300 class).

---

## 9. LINE BUDGET

RewritePlan 3.2: Transport.m ~600 lines (module map; core total ceiling
3.5k).  Indicative allocation: segmentation + geometry (FindSingularities,
ChartRadius, NextCenter, SegmentLine) ~180; DigitBudget + precision raising
~40; matching (MatchWeights + RecombineBasis) ~140; crossing ~80; marching
loop + TransportResult ~110; error probe ~70.  Total ~620 — already tight.

If over budget, cut in this order (function first, sugar last):
1. AccuracyGoalValidate adaptive machinery: ALREADY cut by this spec
   (old Transport.m:776-841 "Before" search and :1191-1213 "After" redo are
   replaced by the single post-hoc check E11 with a suggested
   ExpansionOrder).  Do not reintroduce it to "save" a failing example —
   raise ExpansionOrder in the example config instead.  (Ledger waiver
   required, RewritePlan section 5 seed "expansion-order adaptive search".)
2. Mobius chart support (Q3): if UseMobius is deferred entirely to a
   follow-up, charts are affine-only and ~40 lines of map handling go away;
   requires re-baselining the M4 banana parity config (the oracle run must
   then also use UseMobius -> False) — orchestrator decision, not silent.
3. SavedCharts/SaveExpansions plumbing (~30 lines): move entirely into
   API.m by having TransportLine call a caller-supplied per-chart callback.
4. Benchmark/timing instrumentation: keep a single per-segment timing hook;
   the old BenchmarkData jungle (old Transport.m:161-164, 233-237 etc.) is
   not ported.
NEVER cut: complex-radius segmentation, two-pass digit budget, eps-graded
matching + recombination + asserts, crossing operator mixing, two-point
probe, TWindow tracking, any error of section 5.

---

## 10. OPEN QUESTIONS

Q1  Match-point clipping near singular charts: the old FixWithin midpoint
    rule (Transport.m:1124-1136) intersects THREE intervals (current chart
    validity, pole interval, center-to-pole span) and takes the midpoint.
    This spec preserves the midpoint rule, but with complex-capped radii
    the intersection can be smaller than the old code's; confirm on the
    banana line that segment counts stay comparable (cost, not
    correctness).
Q2  Indicial.md cross-contract: this spec requires per-chart
    "EpsDegenerateFamilies" (sector families whose eps->0 indicial data is
    degenerate — same a mod Z, distinct b, colliding eps=0 eigenvectors) so
    that RecombineBasis is tag-driven rather than numerically detected.
    Indicial.md must export it; if Indicial cannot certify degeneracy
    exactly in some case, the fallback is NOT numeric detection — it is E5
    firing at the match point with a message pointing at the family.
    (Spec-review agents: reconcile with Docs/specs/Indicial.md.)
Q3  Mobius charts in v1 Transport: the M4 parity gate pins banana classic =
    UnequalMassConfiguration WITH UseMobius (RewritePlan M4), but Integrate.m
    rejects Mobius charts (RewritePlan 3.2) and FT requires UseMobius ->
    False (DiffExpIntegration.m:352).  This spec includes Mobius as a Chart
    Map type for classic-transport parity only; if the orchestrator
    re-baselines the banana oracle without Mobius (and the RoC=10 rescaling
    proves sufficient), cut per section 9 item 2.
Q4  AccuracyGoal semantics: this spec demotes AccuracyGoalValidate from
    adaptive search to validation (E11).  Config.md owns the final
    kept/waived table entry for "AccuracyGoal[+Validate]" (RewritePlan 3.2
    Config.m); align before M1.
Q5  AbortOnAnalyticContinuationFail = False: the narrow E8 exception
    (honest INCOMPLETE result at the final chart only) is this spec's
    proposal for the FT use case (DiffExpIntegration.m:411-414 sets False;
    the pentagon "not recognized as a branch point" warning is the disease,
    old Transport.m:470-478).  M0 task 16 (pentagon DeltaPrescriptions
    triage) may conclude the flag can be dropped entirely; Config.md
    decides.
Q6  Probe decrement formula: Ceiling[0.7 * couplingDepth] + 2 is inherited
    verbatim (State.m:222).  With TWindow now tracked honestly, a simpler
    "probe at TWindow CompleteMax minus fixed 3" may be equivalent; keep
    the old formula until M4 parity passes, then revisit (ledger note).
Q7  Exactness ceiling in geometry: NextCenter on Root-object singularities
    can produce nested Root algebra on long lines.  Policy proposal: keep
    centers exact when the expression depth is bounded (LeafCount below a
    Tolerances constant), else rationalize the CENTER (not the
    singularities) at 2x DigitsNeeded with an exactness note in the chart
    record.  Needs a reviewer sign-off because it is a tolerance-gated
    decision in otherwise-exact code.
Q8  Direction-dependent prescription signs: the derived chart Im-sign is
    defined for the chart's negative-t side with orientation-preserving
    chart maps (this spec).  Old code resolves theta_p/theta_m by the SIGN
    of the evaluation point in chart coordinates (Pade.m:70-77) and never
    flips maps; verify during M4 that backward transport (to < from)
    reproduces old left-moving runs (old Transport.m:1026-1031 restricts
    some first-segment cases — that restriction is NOT ported; if a
    genuine asymmetry surfaces, it becomes a loud error, not a silent
    restriction).
