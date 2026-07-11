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
complex-plane distances and old-style real projection waypoints), plans the
chain of expansion charts (predivision,
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
Docs/specs/EpsSeries.md (window-carrying); per DEC-13 its shape is the
EpsSeries.md canonical `<|"EpsWindow" -> <|"Min","CompleteMax"|>,
"Coeffs" -> ...|>` — no alternative eps-Laurent shape appears in this
module.  "exact" means an exact Wolfram expression (Rational, algebraic
number, `Root`), never a float.

### 2.1 `FindSingularities[sys_Association] -> Association`

Input: the LoadSystem record (API.md contract); the fields consumed are
`sys["SingularFactorsExact"]` (the full epsilon-dependent, irreducible
x-dependent factors of all matrix denominators; legacy records may instead
supply `sys["SingularFactors"]`) and `sys["ExtraSingularFactors"]` (threaded
PER CALL from the "ExtraSingularFactors" option of TransportTo and
IntegrateOverLine — DEC-10, REVIEW-math D20; NOT system state and NOT a
config key: the FT layer's IBP-coefficient denominators, old
FeynmanTrick/DiffExpIntegration.m:204-222 `appendMatrixFactors` and
:1326-1385 `CollectLevelIBPSingularFactors`, cannot place charts without
it).

Before root finding, every matrix and extra factor is mapped to its exact
epsilon-zero zero locus: remove the overall epsilon valuation of the exact
rational factor, take its first nonzero epsilon coefficient, and use the
numerator of that coefficient.  Epsilon-independent factors are preserved
verbatim; a constant leading coefficient is dropped (the corresponding root
moves to infinity).  This projection affects only segmentation.  It never
changes `sys["Matrix"]` or `sys["SingularFactorsExact"]`, which remain the
source of truth for chart preparation and indicial analysis.

`DeltaPrescriptions` are projected by the identical rule before chart
attachment.  Otherwise a prescribed factor such as `x + eps` would plan a
limiting chart at `x = 0` but fail the literal vanishing test there, losing its
branch sign.  This projection does not override Solve's `E3` contract for a
moving matrix denominator that degenerates onto the chart origin.

Output:
```
<|
  "All"    -> { exact roots, complex included, deduplicated EXACTLY },
  "Real"   -> { the subset with exactly real value },
  "Projected" -> { old suppressed real projections Re(z), Re(z)+/-|Im(z)| },
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
  A numeric SnapTol pre-filter MAY be used to find merge CANDIDATES on
  large root sets, but two roots are merged only after exact RootReduce
  confirmation (Tolerances.md Q3); an unconfirmed near-pair is kept
  distinct.
- Complex roots are RETAINED as complex numbers and remain the sole source
  for `ChartRadius`.  In addition, the old suppressed real projections
  (Re and Re±|Im|, LineSegmentation.m:81-99) are restored as REGULAR
  predivision waypoints.  They are never treated as matrix poles.  The
  planner simplifies non-real algebraic waypoints to nearby small rationals
  so exact chart preparation does not inherit a large algebraic field.

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
  ordinary, unclipped match point sits at projected geometry
  radius/DivisionOrder of BOTH adjacent charts.  Target-clipped and singular
  transitions use the bounded exceptions below.
- Match points: each ordinary chart's outgoing match point is at physical
  distance `Min[MatchRadius/DivisionOrder, Radius/2]` from its center on the
  marching side.  The first term is the old projected-geometry construction;
  the true half-radius cap aligns every handoff with the error probe
  (old Transport.m:679, 1041-1046: interval ±RoC/DivisionOrder in the
  rescaled chart coordinate); when the next chart is singular the match
  point is additionally clipped into the intersection of both charts'
  validity intervals and centered there (old Transport.m:1124-1136
  `FixWithin` midpoint rule) — preserve the midpoint rule.

  AMENDMENT (recorded, M5i): the FixWithin clip is implemented as the
  BALANCED point of the validity intersection rather than the literal
  interval midpoint.  (The old midpoint operated on ±RoC/DivisionOrder
  design intervals, old Transport.m:635-655; a midpoint of full-radius
  validity intervals would sit at 1/2 of the singular radius whenever the
  singular interval is interior to the producing disk — a 0.5 evaluation
  ratio.)  With gap = dir·(z − prevCenter) > 0, margin = 9/10: the
  intersection of {|m − prevCenter| ≤ margin·prevRad} with the approach
  interval (z − radius_sing, z) is nonempty iff
  gap < den := margin·prevRad + radius_sing (note prevRad ≤ gap always,
  since z is a singularity), and the clip point
    m* = prevCenter + dir·gap·margin·prevRad/den
  equalizes — and thereby minimizes the max of — the two normalized
  evaluation ratios |m − prevCenter|/(margin·prevRad) = |m − z|/radius_sing
  = gap/den.  SegmentLine's cover test admits the singular chart only once
  |m* − prevCenter| ≤ prevRad/k, which bounds the EXECUTED ratios to the
  design classes ≤ 1/k (producing side) and ≤ ~(9/8)/(margin·k) (singular
  side); the un-clipped point z − dir·radius_sing/k is geometry-blind on
  the producing side and lands outside small producing disks (banana
  level 1: anchor radius 1/46 vs singular radius ~0.45 — the naive point
  sits ~4 anchor radii behind the anchor).  Empty intersection = loud E8.
  ONE shared helper (`chartMatchPoint`/`singularMatchPoint`) produces the
  point for SegmentLine's cover target, TransportLine's evaluation, and
  the static audit `ValidatePlan[plan]` (run unconditionally at the top of
  TransportLine): every chart's incoming match point must lie strictly
  inside the producing chart's disk (first chart: the anchor inside its
  own disk; singular charts: additionally on the approach side strictly
  inside their own disk) — violations raise E8 with the full chain
  geometry before any solve is attempted.
- `DigitsNeeded` is computed from the resulting segment count (2.5) and
  stored in the plan.
- The plan is the SINGLE source of truth: `TransportLine` executes it
  verbatim and asserts it never deviates (invariant I1).

### 2.4 `NextCenter[prevMatch_, singularities_, direction_] -> center`

The GetCPL/GetCPR geometry, ported from old Mobius.m:98-142 (affine
charts only; the Mobius variants Mobius.m:72-90 are NOT ported — Mobius
charts are dropped from the new core, DEC-18): given the previous match
point x_b and
its two neighbouring singularities {z_min, z_max}, solve for the next
center x_new such that x_b sits at exactly (new chart radius)/k of the NEW
chart, k = DivisionOrder.  Closed form (affine, finite bounds, right-moving,
old Mobius.m:110-119): with s the new radius,
`s = k(z_max - x_b)/(1+k)` when the new chart is capped on the right (and
the mirror branch when capped on the left), `x_new = x_b + s/k`.
This both-adjacent-charts guarantee is the actual ill-conditioning defense
for the matching solve (RewritePlan R3; legacy review finding 12).
DivisionOrder: classic default 3 (State.m:112); the FT stepwise runner now
uses that same default.

The bounding points for this solve are the restored real projection
waypoints.  Their distance is a conservative geometry scale; the separately
stored true complex-plane `Radius` remains the validity bound.  Every
ordinary adjacent pair therefore shares one exact match point at +1/k in
the producing chart's geometry scale and -1/k in the receiving chart's
scale.  This is the old GetCPL/GetCPR conditioning trick.  The attempted
`StepDivisionOrder` decoupling is retired for this path: a value such as 16
created 28/29 banana-L1 charts and destroyed significance through needless
matching.  The FT runner couples both divisors and defaults them to 3.
Each plan still begins with an anchor chart centered exactly at `from`.
When the mandatory regular chart after a singular crossing is clipped to a
nearby target, the natural pre-clip point MUST be discarded.  The shared
point is recomputed from the actual receiver's
`-Min[MatchRadius/k, Radius/2]` side and accepted
only inside the producer's design and physical disks; otherwise an
intermediate center is inserted.  `ValidatePlan` independently requires all
actual handoffs to lie within one half of both true physical radii, matching
the conservative `SegmentErrorProbe` envelope.

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

With AccuracyGoal == "?": DigitsNeeded := ChopDigits (the chop reserve is
then the only budget), E3 still applies, and E11 is skipped — recorded in
the plan as `"BudgetMode" -> "Unvalidated"`, never silently (REVIEW-math
D25; the old code silently skipped budgeting, Transport.m:776-777 class).

### 2.6 `TransportLine[sys_Association, boundary_, plan_SegmentPlan] -> TransportResult`

The marching loop (pass 2).  `boundary` is either
(a) a LocalSolution anchored at `plan["From"]` (general case: includes
    asymptotic/singular boundary data — the classic "line bcs" of old
    LineSegmentation.m:148-236 arrive here already normalized by API.m), or
(b) a plain values array `vals[[integral, epsorder]]` at a regular `from`
    (wrapped internally into the one-sector regular LocalSolution).
Symbolic indeterminate coefficients in `boundary` are allowed and tracked
per-indeterminate by the error probe (old Transport.m:565-569, 947-962);
they force UsePade off in SectorSeries evaluation via an explicit per-call
option override, NEVER by mutating the configuration store
(REVIEW-minimalism 20; the old global mutation at Transport.m:566-569 is
the forbidden pattern).

Loop body per chart (in plan order):
1. `cs = DiffExp2`Solve`PrepareChart[sys, chart]` — Solve.m assembles the
   lazy exact-SCC envelope (chart map applied to the loaded exact matrix and
   original theta form; diagonal indicial frames are prepared on demand;
   DEC-7, Solve.md §2).  Then `sol = DiffExp2`Solve`SolveChart[cs, req]`
   -> `<|"Basis" -> fundamental LocalSolution basis, "Particular" ->
   LocalSolution | None, "CouplingDepth" -> _Integer|>`, solved
   from the exact original-theta SCC plan.  Depth-zero downstream recurrence
   targets use `SCCExecutionMode -> "BlockSequentialStrict"`.  If a target
   has positive recurrence-pole depth at the guarded work Taylor order, Solve
   lazily full-prepares the same original system and uses
   `"MonolithicStrict"` to avoid repeatedly propagating future epsilon
   orders.  Both modes use the configured strict recurrence solver; this is
   execution coarsening, not an alternate-solver fallback (DEC-9; Solve.md
   §2.4).  `CouplingDepth` remains the longest chain in the exact block
   dependency DAG (old proxy `MaxCouplingOrder` = largest coupled block,
   MatrixLoading.m:357-385, esp. :383).
2. v = incoming data evaluated at chart's MatchIn point (SectorSeries
   evaluate on the previous chart's LocalSolution).  After a singular
   previous chart, `ApplyCrossing` is used exactly when that point has
   negative previous-chart coordinate `t`; a rightward far-side point has
   `t > 0` and must not receive a second branch phase (invariant I6).
3. `w = MatchWeights[...]` (2.7) after `RecombineBasis` (2.8).
4. Assemble the chart's LocalSolution = basis.w + particular via
   `SectorSeries`CombineLocalSolutions` (SectorSeries.md 2.10), set
   EpsWindow by honest window arithmetic, set
   `TWindow["CompleteMax"] = expansionOrder` — DEC-9: there is NO
   (couplingDepth - 1) degradation in the new core (recursion matrices
   are exact polynomials, so a particular built from a lower-block
   solution complete to TOrder is itself complete to TOrder; the old
   GetLargestTerm/MaxCouplingOrder discount, LineSegmentation.m:109-113,
   is waived in the LessonsLedger as "subsumed: exact polynomial
   recursion + ErrorEstimate covers it").
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
EpsSeries (basis solution i, component c: columns are the EpsSeries VALUES
of the basis solutions at the match point). Direct callers may supply an
honest Laurent frame.  The marching path first applies tag-driven
recombination and the epsilon-adic lattice saturation in 2.8.  On an
ODE-regular chart it then applies the constant right-frame normalization
`P = F_eps0^-1` at that match point and re-evaluates the transformed
LocalSolutions, so the leading value matrix is the identity.  It finally
requires the resulting weights to start at or above eps^0.
`incomingValues` = EpsSeries vector.

Contract:
- Graded solve: with F~ = Sum_k F[k] eps^k after recombination/saturation
  (2.8),
  assert det F[0] is nonzero via ``Tolerances`NumericallyZeroQ[detF0,
  scale, Tol["RankTol"], <chart label>]`` (False required; True -> E5,
  band hit -> E7; DEC-3), then solve order by order. In the marching path,
  `ord_eps det F~ != 0` after saturation is LOUD ERROR E5—never a silent
  window shift (RewritePlan 3.4: matchingShift "target 0, asserted"; math
  review finding 4).
- Constant regular-chart frame normalization is epsilon-independent and
  therefore preserves sector tags `(a,b,p)`, prescriptions, valuations, and
  each operand's epsilon window; it is local to the march and MUST NOT be
  installed in the cached chart basis.  The implementation MUST NOT
  threshold small entries of `P`.  It certifies `F_eps0.P == I` with a
  contribution-aware uncertainty scale, recombines whole LocalSolutions,
  re-evaluates them, certifies `F'_eps0 == I`, reruns the saturation audit
  (zero shifts and zero steps required), and asserts that the shared
  `CompleteMax` did not fall or cross the requested top.  Any failed proof is
  E5.  Singular/resonant frames stay on the saturated Laurent path because
  dense mixing of unequal tagged windows can impose an unnecessary common
  top, even though the constant GL action is algebraically valid there.
- Result window: honest min over the windows of incomingValues and
  basisValues; no padding.
- Residual assert: per eps order and component, the ORIGINAL, untrimmed
  system must satisfy
  |F~.w - v| <= Max[matchTol, laurentLeadTol] * scale (Tolerances.m).
  The structural floor is the same relative threshold used to trim Laurent
  pivots; it is dimension-independent.  Resolved nonzero inexact residuals
  also carry their own accuracy uncertainty.  An inexact residual stored
  exactly at zero is accepted in production; `StrictMatchingUncertainty ->
  True` is the opt-in diagnostic mode that additionally requires its
  uncertainty ball to satisfy the same bound.  Violation = LOUD ERROR E6
  carrying the residual, uncertainty, scale and order/component.
- A primary violation does not loosen that contract.  At most two
  deterministic iterative-refinement steps solve
  `F deltaW = v - F.w` using the same Laurent-field elimination, add the
  correction only when it preserves every weight CompleteMax, and then
  rerun the unchanged original-system assertion.  This recovers coherent
  sums of individually sub-floor row-operation terms without a dimension
  factor, tolerance inflation, least-squares rescue, or window shrink.
  The old LinearSolve-with-ZeroTest (old Transport.m:341) and checked
  least-squares rescue (old Transport.m:348-374) remain FORBIDDEN (F6).
- Pivot/leading-coefficient zero classification: EXACTLY
  ``Tolerances`NumericallyZeroQ[pivot, scale, Tol["RankTol"], <chart
  label + operation>]``.  True -> the pivot is exact 0 (snap); False ->
  nonzero; the ambiguity band IS the gray zone and aborts loudly (E7 is
  the band error, not a separate mechanism — REVIEW-minimalism 10,
  DEC-3/DEC-14).  SnapTol plays NO role in matching solves; it is
  reserved for coordinate/endpoint snapping.  (No silent rounding;
  RewritePlan section 5: "numerical-zero leading-coefficient skipping
  (generalizes to matching solves)".)
- After input/rank classification, production cancellation trimming advances
  past exact zeros or inexact coefficients stored exactly at zero.  A
  resolved nonzero Schur coefficient is retained however small.
  `StrictMatchingUncertainty -> True` additionally requires a centered
  inexact zero's full uncertainty ball to lie below the matching residual
  contract; overlap is then E5, never a formal pivot.
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
  a_i - a_j in Z and b_i != b_j, restricted to the families returned by
  ``Indicial`EpsDegenerateFamilies`` (DEC-7: per-chart families with the
  same a mod Z, distinct b, whose eps=0 eigenvectors collide; detected
  EXACTLY by Indicial — chain-top vectors normalized by their
  eps-valuation, eps -> 0 substituted exactly, rank over Q(alpha); the
  record carries the rank deficiency r0 = column count − eps=0 rank,
  REVIEW-math D8).  Division of labor: Indicial proposes the exact common
  families and Transport performs the deterministic divided differences.
- Replacement rule, per degenerate family ordered by b:
  `B_1 = S_1`, `B_m = (S_m - S_1)/((b_m - b_1) eps)` for m >= 2; applied
  recursively to {B_2, ...} while degeneracy persists (recursion depth =
  the family's r0 from the EpsDegenerateFamilies record, REVIEW-math D8).  Each division is an exact EpsSeries window SHIFT.
  Every B_m is a genuine solution (exact linear combination of solutions
  divided by an x-independent scalar).
- Assert after each division: the eps^(-1) (post-shift leading) coefficient
  of B_m's value data is zero via NumericallyZeroQ at Tol["RankTol"]; band
  hit = loud E7 (REVIEW-minimalism 10).
- det effect: each division multiplies det F by 1/((b_m - b_1) eps),
  cancelling exactly one eps order of det degeneracy; "unimodular at eps=0"
  means the FINAL recombined fundamental matrix has det = O(1), assert E5.
- The recombined basis REPLACES the original for all downstream use on this
  chart (assembly, evaluation, probe); the Record field stores the exact
  transformation for debugging/parity dumps.
- After the tag-driven pass, a rank-triggered epsilon-adic closure handles
  any remaining lattice defect. Columns are first normalized by their
  certified epsilon valuations. A division-free determinant series fixes
  `delta = ord_eps det F`; importantly this distinguishes a small constant
  unit from a zero leading determinant followed by material higher orders.
  While the eps^0 matrix is rank deficient, full-pivot elimination produces
  a certified right-null relation `a`, and one target column is replaced by
  `(Sum_i a_i B_i)/eps`. The same operation is applied to the actual
  LocalSolution basis and its match-point frame.
- Each closure step must certify the eps^0 numerator row using same-order
  scales, consumes exactly one complete top coefficient, and lowers the
  determinant valuation by one. Exactly `delta` steps must produce a
  full-rank leading matrix. Gray rank, an unresolved determinant window,
  failed divisibility, early/late termination, or Laurent matching weights
  after closure are loud E5/E7 failures.
- Constant column magnitudes are units, not epsilon valuations. Rank
  elimination normalizes those magnitudes, and a determinant series of
  valuation zero prevents an ill-conditioned but ordinary constant matrix
  from triggering an epsilon division.

### 2.9 `CrossingOperator[sector_, sigma_] -> operator` and `ApplyCrossing[ls_LocalSolution, sigma_] -> LocalSolution`

The analytic-continuation rule for representing a singular chart's negative
arm (chart coordinate t < 0) in the positive reflected coordinate
u := -t > 0.
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
  Consequences for Transport.m: (i) once `ApplyCrossing` has reflected an
  object, that object is evaluated only at POSITIVE u.  Direct negative-t
  evaluations used to match onto the approach arm instead pass the same
  sigma to SectorSeries and do not also apply the operator.  In particular,
  leftward far-side continuation reflects once, while rightward far-side
  continuation already has positive t and does not reflect (invariant I7);
  (ii) the sector-level mixing operator is
  Transport-internal and is DEFINED as the SectorSeries sigma rule
  (its 2.4.1) applied tagwise; Integrate.m's interior pairing and
  negative-arm boundary terms use the SectorSeries rule directly (its
  2.2.5), guaranteeing the conventions coincide (REVIEW-math D28 —
  SectorSeries owns the scalar branch rule; Integrate does NOT call
  Transport).  Sign derivation is NOT reimplemented here:
  Chart["CrossSign"] := ``SectorSeries`ChartImSign`` applied to the
  chart's Prescriptions record; conflict/missing surfaces as ChartImSign's
  `::branchconflict`/`::branchmissing` wrapped into E8's payload
  (REVIEW-minimalism 18).

Family closure: applying the operator to a p > 0 sector POPULATES ALL
lower-p members of the same (a,b) family.  Target sectors that do not yet
exist in the LocalSolution are CREATED; dropping any contribution is
forbidden (F8).  Window arithmetic: the eps^j factor shifts that
contribution's window up by j; the merged sector window is the honest min
(EpsSeries contract).

### 2.10 `SegmentErrorProbe[ls_LocalSolution, tProbe_, couplingDepth_] -> errs[[epsorder]]`

The full-vs-reduced-order probe, ported from old
Transport.m:905-993:
- Evaluate ls at `tProbe`, twice: at full order and
  with the t-series truncated down by
  `probeDecrement = Ceiling[0.7 * couplingDepth] + 2`
  (old `ICurrEvalErrorSeriesDecrease`, State.m:222; reduction applied as in
  old Transport.m:913 via DecreaseSeriesOrderBy).  Both evaluations go
  through SectorSeries (Pade applied there when configured, with ITS loud
  fallback — Transport must not Quiet it, F5).
- `TransportLine` calls the probe at both `-0.51 Radius` and `+0.51 Radius`
  for every regular chart and takes the per-order maximum.
  A regular LocalSolution can be consumed on either side by the next match
  or by `LineIntegral`, so an incoming-side-only probe is nonconservative.
  Interior singular charts are also probed on both signs because
  `LineIntegral` consumes their prescription-aware two-arm tile directly.
  Only a singular endpoint is probed on its actual incoming side.
  The 0.51 radius is deliberately outside the planner/API's accepted
  half-radius envelope and never the chart center, where log-bearing terms
  vanish spuriously.
- The per-epsorder error is `|full - reduced|`; the marching call takes the
  maximum over the two regular-chart signs before accumulation.
- Accumulation across segments is ADDITIVE into
  `TransportResult["ErrorEstimate"]` (old Transport.m:980-984); abort > 1 is
  E10.  Seeding: when `boundary` carries no error data the seed is exact 0
  with `"ErrorSeeded" -> False` metadata — the old Accuracy[]-based seeding
  and PadRight zero-padding (Transport.m:583-591) are FORBIDDEN (F11).
- Per-segment digit check: when AccuracyGoalValidate is enabled (any
  non-False value; default False per DEC-5) and AccuracyGoal is numeric,
  the per-segment error must satisfy err <= 10^-DigitsNeeded; violation =
  LOUD ERROR E11 naming the segment and the order increase that would be
  needed.  (DEC-11: AccuracyGoalValidate is demoted to validation-only;
  the old adaptive Before/After machinery is NOT ported — Config's table
  records the waiver "replaced by exact recursion + ErrorEstimate gate".
  See section 9 cuts.)

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution / Sector / EpsWindow (RewritePlan 3.1; the ChartMap
and TWindow comments are superseded by DEC-18 and DEC-9 as noted)

```
LocalSolution = <|
  "Center" -> exact x0, "ChartMap" -> affine (Mobius DROPPED, DEC-18),
  "Radius" -> true complex-plane distance expressed in the local affine
              coordinate t, i.e. PhysicalRadius/ChartMap["Scale"],
  "Sectors" -> { Sector.. },
  "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "TWindow"   -> <|"CompleteMax" -> nmax|>,   (* t-order truncation
              record ONLY (DEC-9): no coupling-depth degradation in the
              new core — recursion matrices are exact polynomials; the
              old MaxCouplingOrder/ISafetyExpansionSubtract discount is
              waived in the LessonsLedger *)
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
  "Scale"      -> projectedGeometryRadius / RadiusOfConvergence,
                  (* physical x = Center + Scale t; this is the affine part
                     of old GetLineRescaled.  The projected radius preserves
                     the classic Re/Re±Im geometry, while true roots still
                     own the validity bound. *)
  "Radius"     -> physical true complex-plane radius (used by plan
                  validation and API physical tiling),
  "MatchRadius" -> physical projected-geometry radius,
  "LocalRadius" -> Radius/Scale (copied to LocalSolution["Radius"]),
  "IsSingular" -> True | False,
  "IncomingMatchPoint" -> physical x of the incoming match point (absent
                  for the first chart; converted through Scale only at the
                  solve/evaluation boundary),
  "CrossSign"  -> +1 | -1 | None  (derived chart Im-sign; None iff not
                  crossed or no multivalued content),
  "Prescriptions" -> as in 3.1,
  "EpsDegenerateFamilies" -> from `` Indicial`EpsDegenerateFamilies ``
                  (DEC-7); populated when the chart is prepared/solved
                  (2.6 step 1), NOT by SegmentLine geometry (see 2.8),
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
  "SegmentCount" -> n
|>
```
SavedCharts/SaveExpansions is NOT a TransportResult field: per-chart
LocalSolution saving lives in API.m, which passes TransportLine an
optional caller-supplied per-chart callback (pre-authorized cut 3,
section 9; REVIEW-minimalism defect 4).  API.m wraps this into the
user-facing TransportTo association (kinematic point, NumIntegrals,
orders — API.md contract; old return shape at old Transport.m:1234-1246).

### 3.5 Weights

`MatchWeights` returns an EpsSeries vector with an explicit window
`<|"Min" -> 0, "CompleteMax" -> ...|>`; the eps^0 window start is asserted
for ALL columns post-trim (2.7; REVIEW-math D7 — not only "regular
incoming data"); any Min < 0 for data whose own window starts at 0 is E5
territory (it means recombination failed).

---

## 4. INVARIANTS (always-on cheap asserts)

I1  Plan fidelity: TransportLine executes exactly plan["SegmentCount"]
    segments with exactly the planned centers/match points; any deviation
    aborts (kills the old duplicated-dry-run drift, old Transport.m:666).
I2  Geometry: every ordinary, unclipped match point m between charts i and
    i+1 obeys +1/DivisionOrder and -1/DivisionOrder in the two projected
    geometry scales within snapTol (the GetCPL/GetCPR guarantee,
    Mobius.m:98-142).  Singular/target-clipped transitions may be asymmetric,
    but the stored ACTUAL shared point must lie within one half of both true
    physical radii; the receiver-side point is recomputed after every clip.
    This half-radius bound is the `SegmentErrorProbe` contract, not merely a
    convergence check.
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
I10 Regular chart <=> exactly one sector (a=0, b=0, p=0) PER HOMOGENEOUS
    basis column (RewritePlan 3.1, scope per math finding 8: ODE solutions
    only).  DEC-6 alignment: at a regular chart Indicial returns ONE
    (0,0,0) family with d-dimensional coefficient space (d sector specs),
    so the fundamental basis is d columns, each with exactly one (0,0,0)
    sector; the assert checks every column.
I11 ODE residual spot-check at one random interior point per chart
    (RewritePlan 3.1 invariant; executed on the assembled chart
    LocalSolution, residual <= matchTol*scale).

---

## 5. ERROR CONTRACT

All errors are raised via ``Tolerances`DE2Error[id, payload]`` (DEC-1: THE
library-wide loud-error primitive — prints a one-line summary and
`Throw[Failure["DiffExp2", payload], "DiffExp2Error"]`, caught only at
API.m entry points); this module defines no other error mechanism.  The
payload always carries "ID", "Module" -> "Transport", the Chart "Label"
(chart index + center + singular/regular) where applicable, plus the named
fields listed.  There is NO error class in this module that downgrades to
a warning, with the single structured exception noted in E8.

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
E7  TOLERANCE GRAY ZONE.  A pivot/leading coefficient falls inside the
    NumericallyZeroQ ambiguity band around Tol["RankTol"]*scale (matching
    solve or recombination post-shift leading coefficient; E7 IS the band
    error of DEC-3, not a separate mechanism — REVIEW-minimalism 10).
    Carries: chart label, the value, both band edges, the operation.
E8  PRESCRIPTION CONFLICT OR MISSING.  Fires when the chart is multivalued
    AT ALL: b != 0 OR p > 0 OR Denominator[a] > 1 (DEC-16; the
    pentagon-triage Kallen charts are fractional-a with b = 0, p = 0 and
    MUST trigger this) and the Prescriptions list yields no consistent
    derived Im-sign or no prescription covers the chart at all.  The sign
    derivation is ``SectorSeries`ChartImSign`` (REVIEW-minimalism 18 —
    not reimplemented here; rules per old AnalyticContinuation.m:45-68:
    even vanishing multiplicity = no constraint; required sign =
    prescribed sign divided by leading-coefficient sign); its
    `::branchconflict`/`::branchmissing` are wrapped into this error's
    payload.  Prescription-factor dedup is SIGN-AWARE (DEC-16:
    {-1+x, +1} vs {1-x, +1} flip the implied i-delta side; conflicting
    normalizations are an error here — the old sign-blind
    deltaPrescriptionsForFactors dedup is forbidden).
    Carries: chart label, every vanishing factor with its prescribed
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
    Anywhere else: error.  (DEC-15: this final-chart-only rule is
    NORMATIVE across Transport/API/Config — REVIEW-math D14; pentagon
    triage does not override it.)
E9  CROSSING WITHOUT SIGN.  ApplyCrossing invoked with CrossSign None on a
    chart with multivalued content (should be unreachable given E8; kept as
    a cheap assert).  Carries: chart label, sector tags.
E10 ERROR ESTIMATE > 1.  Any accumulated ErrorEstimate entry exceeds 1.
    Carries: chart label of the segment that tipped it, the (integral,
    eps-order) entry, the per-segment and accumulated values.  (Port of old
    Transport.m:991-993.)
E11 SEGMENT ACCURACY MISS.  AccuracyGoalValidate enabled (DEC-11
    validation-only semantics; DEC-5 default False), AccuracyGoal numeric,
    and a segment's probe error > 10^-DigitsNeeded.  Carries: chart label,
    measured error, DigitsNeeded, current ExpansionOrder, suggested
    ExpansionOrder.
E12 WINDOW UNDERFLOW.  A consumer (matching, crossing, evaluation,
    assembly) requests an eps order > EpsWindow["CompleteMax"] or a t order
    > TWindow["CompleteMax"] of any operand.  Carries: chart label, sector
    tags, requested order, available CompleteMax.  (RewritePlan 3.1: "fail
    loudly naming (chart, sector, order)".)
E13 NO MARCHING PROGRESS.  A planned step does not strictly advance toward
    `to`, or executed segments exceed plan["SegmentCount"] (I1).  Carries:
    chart label, match points.  (Old code can loop forever; NEW.)
E14 SUSPICIOUSLY NEAR-SINGULAR ENDPOINT.  Input-accuracy guard rule
    (Tolerances.md $NearSingularityGuardDecades entry; DEC-14: user-input
    target snapping is INPUT-scaled, never wp-scaled — a 1e-16-off user
    target must snap, not create a degenerate chart).  An INEXACT `to`
    with |to - root| < 10^(-Floor[Accuracy[to]/2]) for some exact
    singularity root is SNAPPED to the root with a visible warning (the
    input's own accuracy sets the snap scale); an inexact `to` with
    10^(-Floor[Accuracy[to]/2]) <= |to - root| <
    10^(-$NearSingularityGuardDecades) is THIS loud error ("target is
    suspiciously near singularity <root>; pass the exact value, or an
    exact offset, to confirm intent" — REVIEW-minimalism 11).  Exact `to`
    is never snapped and never guarded — exact membership decides
    (chop-based detection of the old Transport.m:606-609 is forbidden,
    F13).  Carries: `to`, the nearby root, their exact difference, both
    thresholds.

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
F4  Treating projected Re, Re±Im waypoints as actual matrix singularities.
    FORBIDDEN.  The waypoints are required for classic predivision geometry,
    but only true roots in "All" cap `ChartRadius` or produce singular charts.
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
| DiffExp/Transport.m:776-841, 1191-1213 | AccuracyGoalValidate "Before" adaptive order search / "After" redo | NOT ported; replaced by E11 (DEC-11 validation-only; Config's table records the waiver "replaced by exact recursion + ErrorEstimate gate") |
| DiffExp/Transport.m:287-441 (in IntegrateSystem) | boundary fixing: equations, LinearSolve, least-squares rescue, NullSpace free params | `MatchWeights` (2.7) + E5/E6/E7; rescue/NullSpace forbidden F6/F7 |
| DiffExp/Transport.m:317-327 | FixAt != 0 coefficient extraction by differentiation through nested SeriesData | obsolete: exact tags + SectorSeries evaluation |
| DiffExp/Transport.m:905-993 | SEval1/SEval2 evaluation + two-point error probe + per-indeterminate errors + abort > 1 | `SegmentErrorProbe` (2.10) + E10 |
| DiffExp/Transport.m:947-962 | ComputeErrorsPerIndeterminate | inside 2.10 |
| DiffExp/Transport.m:1234-1246 | return association | TransportResult (3.4) via API.m |
| DiffExp/LineSegmentation.m:65-106 `FindMatrixSingularities` | factor collection, Quiet Solve, chop dedup, ghost projection | `FindSingularities` (2.1); F3/F4 |
| DiffExp/LineSegmentation.m:109-113 `GetLargestTerm` | matrix-truncation magnitude at coupling-discounted order | WAIVED (DEC-9: the degradation premise — truncated matrix expansions — is removed; error control = probe (2.10) + E11; ledger entry "subsumed: exact polynomial recursion + ErrorEstimate covers it") |
| DiffExp/LineSegmentation.m:116-145 `GetMatricesPrecisionDistance` | Dynamic segmentation intervals | DROPPED (SegmentationStrategy = Predivision only in v1, RewritePlan Config.m) |
| DiffExp/LineSegmentation.m:20-62, 239-258 `RelateLines`/`GetMatchingPoint` | Solve-based relation between segment parametrizations | obsolete: explicit invertible ChartMaps composed exactly |
| DiffExp/LineSegmentation.m:148-236 `CheckBoundaryConditionsAndReparametrize` | bcs validation/reparametrization | API.m's contract; Transport asserts normalized input |
| DiffExp/Mobius.m:22-35 `GetMobius` | Mobius chart maps | DROPPED (DEC-18: Mobius is out of the new core entirely; the RoC rescaling — an affine rescaling, not a Mobius map — is KEPT) |
| DiffExp/Mobius.m:38-69 `GetLineRescaled` | chart construction + RoC rescaling + 2x WP SetPrecision | Chart construction in SegmentLine; L1/L4 |
| DiffExp/Mobius.m:72-142 GetMobiusCPL/R, GetCPL/GetCPR | next-center geometry (1/k of both charts) | `NextCenter` (2.4; affine GetCPL/GetCPR only — Mobius variants dropped, DEC-18) |
| DiffExp/Mobius.m:145-155 FindNextCenterPointL/R | bound selection wrapper | inside `NextCenter` |
| DiffExp/AnalyticContinuation.m:18-90 `PrepareAnalyticContinuation` | per-chart vanishing-factor collection, sign derivation, replacement rules, conflict flag | sign derivation -> Prescriptions construction (3.1, shared contract with SectorSeries/API); rules -> `CrossingOperator` (2.9); conflict -> E8 |
| DiffExp/AnalyticContinuation.m:93-103 `Project\[Theta]s` | theta-projection of two-sided data | subsumed by one-sided data + explicit crossing |
| DiffExp/Pade.m:55-91 SEval1/SEval2/SEval | evaluation with continuation + theta resolution | evaluation -> SectorSeries (with Pade, RewritePlan 3.2); theta resolution -> ApplyCrossing |
| DiffExp/MatrixLoading.m:357-385 `InitializeIntegrationSequence` | coupled-block ordering + MaxCouplingOrder | block ordering -> Solve.m (SolveChart, DEC-9); CouplingDepth consumed here for the probe decrement ONLY (not TWindow — DEC-9) |
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
    projected geometry radius/DivisionOrder of BOTH charts for every natural
    unclipped join (GetCPL/GetCPR, Mobius.m:98-142).  Target/singular clips
    recompute the actual receiver point and remain within both true
    half-radius envelopes.  This, not the solver alone, is the conditioning
    defense for matching (R3).  Classic and FT default k=3.
L4  RoC chart rescaling.  Chart coordinate scaled by
    projectedGeometryRadius/RoC so high-order coefficients stay O(1), while
    the separately stored true complex radius divided by that scale remains
    the LocalSolution validity bound.  Banana classic line uses RoC = 10
    (Reference/Examples Banana_example.m; RewritePlan section 5).
L5  Complex singularities have two distinct roles: exact complex roots cap
    convergence radii, while the old suppressed Re, Re±Im projections are
    regular real waypoints for stable symmetric predivision matching.
L6  Two-point error probe.  Full vs reduced order at BOTH the incoming
    and outgoing 0.51-radius points of regular and interior-singular charts
    (incoming side only for a singular endpoint), skipping the chart center
    where log terms vanish spuriously; additive accumulation and hard abort
    when > 1.
L7  Coupling-depth t-order degradation — premise REMOVED (DEC-9).  The
    old top-(couplingDepth - 1) t-orders were garbage because chained
    particulars were built from TRUNCATED numeric matrix expansions; the
    old code discounts them in the matrix-error read-off
    (ExpansionOrder - ISafetyExpansionSubtract - (MaxCouplingOrder - 1),
    LineSegmentation.m:109-113, ISafetyExpansionSubtract = 5 State.m:216)
    and in the probe decrement (Ceiling[0.7 MaxCouplingOrder] + 2,
    State.m:222).  New core: recursion matrices are exact polynomials, so
    TWindow["CompleteMax"] = ExpansionOrder with NO discount (LessonsLedger
    entry "subsumed: exact polynomial recursion + ErrorEstimate covers
    it"); the probe decrement formula is preserved verbatim as a heuristic
    input (couplingDepth from Solve metadata).
L8  Crossing convention.  sign +1 = NO replacement, relying on
    principal-branch Log/Power on the far side; sign -1 = Logx ->
    Logx - 2 Pi I theta_m SHIFT (binomial log-chain mixing), e^(-2 Pi I b)
    only for Denominator[b] > 2, half-integers through explicit Sqrt with
    (theta_p - theta_m) (AnalyticContinuation.m:70-79).  The single formula
    of 2.9 must reproduce all of these (unit test T8).  9aeb300 lesson:
    crossing phases come from the prescription exactly once.  Reflection
    `t -> -u` is required only for a negative target coordinate; applying it
    to an already-positive rightward far-side point adds a spurious full
    monodromy.  PV-paired
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
    (RewritePlan section 5; DEC-3: the one library-wide predicate is
    ``Tolerances`NumericallyZeroQ`` at Tol["RankTol"], whose ambiguity
    band IS the gray zone).

---

## 7. DEPENDENCIES

Order (acyclic, RewritePlan):
`Tolerances < Config < EpsSeries < SectorSeries < Indicial < Solve <
Transport/Integrate < API`.

Transport.m MAY call:
- Tolerances.m: DE2Error (DEC-1), MatchTol, RankTol (only via
  NumericallyZeroQ, DEC-3), ResidTol, SnapTol (computed-value/geometry
  snapping only — no role in matching solves), ChopReserve, $SafetyDigits,
  $InputPrecisionFactor, EvalErrorSeriesDecrease (probe constants),
  $MinExpansionOrder, $NearSingularityGuardDecades (E14 input-accuracy
  guard).  (Names per the Tolerances §3 consumer table as amended by
  REVIEW-math D16 / REVIEW-minimalism 17.)
- Config.m: WorkingPrecision, AccuracyGoal, AccuracyGoalValidate
  (validation-only gate for E11 — DEC-11; default False — DEC-5),
  DivisionOrder, RadiusOfConvergence, ExpansionOrder, EpsilonOrder,
  EstimateError, SegmentationStrategy (must be "Predivision"; anything
  else is Config's loud error), UsePade (read-only; enforcement in
  SectorSeries), DeltaPrescriptions, AbortOnAnalyticContinuationFail,
  Verbosity via PrintInfo/PrintWarning (the only print helpers,
  REVIEW-minimalism 5).  UseMobius is NOT read (Mobius dropped, DEC-18);
  SaveExpansions is API.m's business (caller callback, 3.4).
- EpsSeries.m: all window-carrying eps-Laurent arithmetic, the graded/
  Laurent linear solve used by MatchWeights, exact eps-division (window
  shift) used by RecombineBasis, ESTrim/ESCoeffZeroQ (2.7 column trim).
- SectorSeries.m: evaluate (with Pade and ITS loud fallback), multiply,
  re-expand, CombineLocalSolutions (SectorSeries.md 2.10), differentiate,
  ChartImSign (the single sign-derivation authority, REVIEW-minimalism
  18).
- Indicial.m: per-chart sector specs, Prescriptions derivation inputs,
  EpsDegenerateFamilies (DEC-7).
- Solve.m: PrepareChart (lazy exact-SCC envelope plus on-demand diagonal or
  full indicial preparation — DEC-7) and SolveChart (strict
  `BlockSequentialStrict` or cost-coarsened `MonolithicStrict` fundamental
  basis + particular + exact CouplingDepth — DEC-9).

Transport.m MUST NOT call: Integrate.m (sibling — no calls in EITHER
direction; the shared branch convention is SectorSeries' sigma rule, not
a sibling export — REVIEW-math D28), API.m, anything in FeynmanTrick/,
anything in the frozen DiffExp/ tree.

---

## 8. UNIT TESTS

Closed-form/pure tests run without the old library; tests marked (kernel,
oracle) consume M0-pregenerated artifacts and belong to the M4 gate but are
enumerated here because they pin THIS module's behavior.

Error tests assert via `Catch[..., "DiffExp2Error"]` returning the
`Failure["DiffExp2", ...]` with the named "ID"/"Module" payload fields
(DEC-1).  COVERAGE RULE (M0 review): every error ID and warning ID of
section 5 has at least one unit test that triggers it and asserts its
required payload fields, OR a one-line waiver in the test file naming the
ID and why it is untestable in isolation (e.g. "unreachable by
construction, guarded by assert X").  The IDs currently missing tests are
enumerated in Docs/specs/REVIEW-minimalism.md defect 25; the
implementation agent closes the list before the module's milestone gate.

T1  `test_segmentation_exact_roots`
    Factors {1 - 4 x, x, 1 - x}; SegmentLine from 11/23 to 0.  Asserts: the
    singular roots are EXACTLY {0, 1/4, 1} (SameQ); the plan contains a
    singular chart with Center === 1/4; plan["EndpointIsSingular"] ===
    True (endpoint 0).  (Box L1 geometry; campaign interior pole at 1/4.)
T2  `test_segmentation_complex_radius`
    Factors {x^2 + 1, x - 3}.  Asserts: ChartRadius[0, all] == 1 — the
    COMPLEX roots ±I are binding (real root 3 is farther), Projected contains
    {-1,0,1,3}, and the -1,0,1 charts are regular.  ChartRadius[2, all] ==
    Min[|2 - 3|, |2 - I|] == 1: projections never cap true radius.
T3  `test_segmentation_exact_complex_roots`
    Factor x^2 - 2 x + 2 (roots 1 ± I).  Asserts ChartRadius[1, all] == 1
    with the exact algebraic roots retained in "All".
T4  `test_nextcenter_both_charts`
    Singularities {0, 1}, previous match point 1/2, k = 4, right-moving.
    Asserts NextCenter == 3/5 exactly, new radius == 2/5, and
    |1/2 - 3/5| == (2/5)/4 == 1/10 (the GetCPL closed form,
    Mobius.m:110-119).  Property sweep: for 50 random rational
    {z_min, x_b, z_max, k in 2..6}, the natural unclipped 1/k relation holds
    for BOTH adjacent projected geometry scales within snapTol (I2).  Separate
    target/singular-clip properties assert the actual point stays within both
    true half-radius envelopes.
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
    incoming data = the eps-regular solution B2 itself, i.e.
    ((x^(2eps)-1)/(2 eps), x^(2eps)) evaluated at t_m = 1/3 (eps^0 value
    (Log[1/3], 1)), MatchWeights == (0, 1) to matchTol, with weight window
    starting at eps^0 (NOT eps^-1).  (Second component corrected at M0
    review: x^(2eps)-1 is not a solution of this system — REVIEW-math D9 /
    REVIEW-minimalism 22.)
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
