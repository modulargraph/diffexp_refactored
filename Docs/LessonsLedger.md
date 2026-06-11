# Lessons Ledger (M0 deliverable)

Status: v1, 2026-06-11.  Execution-contract companion to `Docs/RewritePlan.md`
(section 5 seeds; all seeds appear below with dispositions).  Sources: full
read of `DiffExp.m`, every file under `DiffExp/` (incl.
`DiffExp/IntegrationStrategies/`), `FeynmanTrick/DiffExpIntegration.m`
(eps-window/trim/prefactor logic), the campaign docs
(`Docs/FeynmanTrickBoxFamilyStatus.md`, `Docs/FeynmanTrickBananaStatus.md`),
`Scripts/run_ft_stepwise.m`, and the three-lens review
(`Docs/reviews/rewrite_plan_review_3lens.json`; all 18 "legacy" findings are
folded in, plus the operationally relevant "math" findings).

This ledger exists so the rewrite cannot silently lose a decade of
numerical-operational wisdom.  Implementation agents: every DISPOSITION line
is binding.  An item marked `implement-in-<module>` must appear in that
module's spec/tests; `subsumed-by-<feature>` means the mechanism dies but the
*requirement it served* is discharged by the named design feature (verify the
discharge, do not just delete); `waived-because-<reason>` requires no code.

Line numbers refer to the working tree at the M0 baseline (current master,
old library frozen as parity oracle).  "WP" = WorkingPrecision, "CP" =
ChopPrecision, "FEC" = the `DiffExp`State`` config association.

---

## A. Precision management

### L1 Raise input precision to 2x WP
WHAT: Every external numeric input (line, endpoint, boundary values) with
`Accuracy < FEWorkingPrecision` is re-fixed with `SetPrecision[..., 2*WP]`:
Transport.m:38-41 (PrepareBoundaryConditions), Transport.m:172-175
(IntegrateSystem), Transport.m:542-544 and 559-562 (TransportTo line and
endpoint); chart rescaling re-fixes every produced segment line at 2x WP
(Mobius.m:48, Mobius.m:66).  Evaluations run under
`$MinPrecision = FEWorkingPrecision` blocks (Pade.m:34, 70, 80).
WHY: The 2x headroom is what survives chart-map composition and repeated
`Together/Collect` over hundreds of segments; without the `$MinPrecision`
floor, significance arithmetic silently bleeds digits at every substitution
and Mathematica downgrades to machine precision mid-chain.
DISPOSITION: implement-in-Transport.m (input normalization) and
implement-in-Tolerances.m (the 2x factor and `$MinPrecision` policy are named
constants, not folklore).

### L2 $MaxExtraPrecision = 1000 in chart construction
WHAT: `GetLineRescaled` runs inside `Block[{$MaxExtraPrecision = 1000}, ...]`
(Mobius.m:39).
WHY: Exact-rational Mobius/affine maps composed with high-precision numerics
need deep extra-precision headroom for `Simplify/Together` zero tests;
the default 50 makes exact cancellation tests fail spuriously near chart
edges.
DISPOSITION: implement-in-Transport.m (chart-map construction block).

### L3 ChopPrecision < WorkingPrecision invariant; LinearSolveChopPrecision auto-sync
WHAT: Any update of WP or CP errors if `CP >= WP` (DiffExp/DiffExp.m:121-125)
and resets `"LinearSolveChopPrecision" = FEC[ChopPrecision]` on EVERY config
update, overridable only explicitly afterwards (DiffExp/DiffExp.m:126-130).
Defaults: WP 500, CP 250 (State.m:109, 119, 135).  The FT layer pins
`ChopPrecision -> precision - 50` (FeynmanTrick/DiffExpIntegration.m:349),
i.e. a 50-digit chop margin rather than WP/2.
WHY: Chopping above working precision deletes real data; a stale
LinearSolveChopPrecision after a WP change makes `ZeroTest` in every
`LinearSolve` wrong (silent rank misjudgement).  The FT margin choice (WP-50)
exists because FT towers carry ~30 digits of legitimate cancellation residue
(see L9) that WP/2 chopping would misclassify.
DISPOSITION: implement-in-Tolerances.m (all chop scales derived from one WP,
with the relationships asserted at config time) and implement-in-Config.m
(validation, single update path).

### L4 Re-fix precision before linear solves (zero-accuracy zeros)
WHAT: `SafeNumericLinearSolve` re-applies `SetPrecision[..., wp]` whenever
`Precision[input] < wp` before solving (ResonantRecurrence.m:536-544, comment
"the banana x = 1/2 eps^-1 impurity").
WHY: Catastrophic cancellations upstream leave zero-accuracy zeros (`0``a`)
that drag the whole array's Precision to 0; LinearSolve/LeastSquares then
fail or return a zero "solution" for a nonzero rhs, silently dropping a
recurrence source term.  The data is exact at WP by construction, so
re-fixing is correct; letting significance arithmetic veto the solve is not.
DISPOSITION: implement-in-Solve.m (linear-solve wrapper) — DiffExp2's
eps-Laurent coefficient solves hit exactly the same hazard.

### L5 Compare residuals as "accept unless PROVEN above tolerance"
WHAT: The unified-particular residual gate is
`If[!TrueQ[maxResidual > tolerance], accept]` — not
`If[maxResidual < tolerance]` (ResonantRecurrence.m:1230-1239).
WHY: Residuals of steps fed by accuracy-limited data are significance zeros
whose comparison against a tolerance is *undetermined*; the naive comparison
rejects exact solves.  General rule: undetermined numeric comparisons must
default to the non-destructive branch and only PROVEN violations may trigger
rejection/retry.
DISPOSITION: implement-in-Tolerances.m (a named comparison helper with these
semantics, used by all residual gates).

---

## B. Tolerances and zero tests

### L6 The tolerance map (disease D3 catalogue)
WHAT: The old code uses, scattered: `10^-CP` (PChop, Utilities.m:103),
`10^-LinearSolveChopPrecision` (LSPChop, Utilities.m:104; all LinearSolve
ZeroTests/NullSpace Tolerances, e.g. Transport.m:341, Wronskian.m:62-63),
`10^-30` fixed crosscheck chop (State.m:149), `10^-5` multivalued check
(State.m:212), `10^(-CP/2)` eigenvalue/exponent snapping
(Recurrence.m:136, 150, 618; ResonantRecurrence.m:371, 1047),
`10^(-CP/3)` singular-step/Jordan tolerances (ResonantRecurrence.m:602, 657,
893, 1164), `Max[tol^(1/3), 10^-12]` fit consistency (RegularizedIntegration.m:941),
`Max[10^(-CP/2), 10^-24]` leading-skip floor (Recurrence.m:618),
`10^-2` root clustering (RegularizedIntegration.m:862),
`RationalizationTolerance = 10^-10` structural rationalization (State.m:123),
`10^-40` absolute FT Laurent-zero test (DiffExpIntegration.m:1177-1178) and
snap fallback (DiffExpIntegration.m:154-155, 182-183), `10^-Min[80,
Max[20, Floor[LSCP/4]]]` least-squares acceptance (Transport.m:362),
`100*tolerance` solve-residual warning (ResonantRecurrence.m:577),
`10^(-activePrecision/2)` Hankel rank gate (RegularizedIntegration.m:1047),
`10^100` evaluation blow-up gate (RegularizedIntegration.m:1845).
WHY: Each constant encodes a discovered failure scale (documented in the
items below); their SCATTER is the disease.  Several pairs interact (e.g. a
chop at `10^-CP` upstream makes a `10^(-CP/2)` test downstream see exact
zeros) and have caused real bugs when changed independently.
DISPOSITION: implement-in-Tolerances.m — ONE module defining named,
WP-derived tolerances (chopFloor, matchTol, snapTol, rankTol,
laurentLeadTol, residualTol, blowupGuard) with comments mapping each to the
ledger item that motivates it.

### L7 Chop-as-you-go series policy
WHAT: `SExpand` applies `PChop@*Expand` to every series coefficient after
every operation (SeriesOps.m:64-67); all strategy outputs pass through
`PChop` (e.g. Recurrence.m:272, 789; ResonantRecurrence.m:853, 1042).
WHY: Without continuous chopping, cancellation residues at `10^-CP` scale
accumulate, masquerade as nonzero leading coefficients (L9), inflate Logx
power scans, and poison structural decisions.  Chopping is part of the
numerical contract, not cosmetics.
DISPOSITION: implement-in-EpsSeries.m and implement-in-SectorSeries.m
(coefficient normalization after every arithmetic op, at the Tolerances.m
chopFloor).

### L8 RationalizationTolerance is for STRUCTURE, not values
WHAT: `RationalizationTolerance` (default `10^-10`, State.m:123) is used
only to recognize exact structure in numerics: integer Logx powers
(SeriesOps.m:95-106), rational x-powers and extracted exponents a, b
(SingularityDecomposition.m:52-64, 336-364), zero/negative tests on
exponents (RegularizedIntegration.m:258-284), endpoint snapping.  It never
chops solution coefficients.
WHY: Structural snapping needs a LOOSE tolerance (these quantities are small
exact rationals polluted by fit noise far above chopFloor), while value
chopping needs a TIGHT one; conflating them either destroys data or fails to
recognize structure.  This is the seed of the Tolerances.m "named semantics"
requirement.
DISPOSITION: implement-in-Tolerances.m (snapTol distinct from chopFloor);
mostly subsumed-by-exact-sector-tags (DiffExp2 carries a, b, p exactly so
recognition largely disappears; snapping survives only at LoadSystem
validation and root-finding edges).

### L9 Numerical-zero leading-coefficient skipping is RELATIVE with a 1e-24 floor
WHAT: `ComputeSingularParticular` determines the source's leading power by
skipping coefficients that are small RELATIVE to the series' own magnitude,
probing Logx at an irrational value and BOTH theta branches, with threshold
`tolRel = Max[10^(-CP/2), 10^-24]` (Recurrence.m:587-639, esp. 594-618; the
comment calls the 1e-24 floor load-bearing).
WHY: Box L2 apparent-singularity sources carry cancellation residues at
~1e-29 RELATIVE at WP 300 — far above PChop and above `10^(-CP/2)` —
below the true leading power; reading nmin naively shifts the Frobenius
ansatz exponent s onto an eigenvalue (s - lambda = 0) and silently corrupts
the particular from that eps order onward (the -6.69 box eps^0 deficit,
Docs/FeynmanTrickBoxFamilyStatus.md:104-119).  Genuine leading content is
either O(1) relative or exactly zero in exact arithmetic, so a loose
relative floor is safe AND necessary.
DISPOSITION: implement-in-Solve.m (source leading-power detection) AND
implement-in-Transport.m (the same skip generalizes to leading-weight
detection in matching solves — seed text: "generalizes to matching solves").

### L10 Two distinct snapping scales: spectrum vs step
WHAT: Eigenvalues/exponents are snapped to rationals at `10^(-CP/2)`
(Recurrence.m:135-136; ResonantRecurrence.m:380-384, 424-425,
NormalizeSourcePower ResonantRecurrence.m:1045-1048) while per-step
singularity decisions (is L_n singular, divisor zero) use the looser
`10^(-CP/3)` (ResonantRecurrence.m:602-607, 657, 716, 736).
WHY: Snapping must be tight enough not to merge genuinely distinct
eigenvalues but loose enough to absorb float noise; step-singularity
detection must be looser because noise has been amplified by the recursion
by the time a near-singular L_n appears.  Getting these backwards causes
either missed resonances (division by software zero) or spurious ones.
DISPOSITION: implement-in-Indicial.m and implement-in-Tolerances.m — note
that exact indicial algebra removes MOST of this (resonance becomes exact
integer arithmetic on a, b), so the surviving need is only LoadSystem-side
validation; record the constants for the legacy parity path.

### L11 Zero tests must evaluate BOTH theta branches and probe Logx
WHAT: `EffectiveZeroExprQ`/`EffectivelyZero` substitute both
`{thetap->1,thetam->0}` and `{thetap->0,thetam->1}` before judging zero
(RegularizedIntegration.m:286-308; SingularityDecomposition.m:257-299);
`CheckParticularResidual` evaluates residuals on both branches instead of
letting `FiniteAbsMax` treat symbolic theta content as 0
(ResonantRecurrence.m:1003-1015); leading-coefficient scans probe
`Logx -> N[Pi]` (Recurrence.m:600-607).
WHY: Coefficients carry analytic-continuation theta symbols and Logx towers;
a zero test that sees an unevaluated symbol and answers "not numeric, treat
as 0" (or "nonzero") makes drop/keep decisions on the wrong branch — this
blinded the residual check during the banana campaign.
DISPOSITION: subsumed-by-exact-sector-tags for Logx (log content is the
explicit p tag, never an embedded symbol); implement-in-Tolerances.m for the
branch-aware zero test where theta-like branch data still appears
(Transport.m crossing data).

### L12 Never IntegerQ (or unguarded Round) on floats for resonance tests
WHAT: The non-resonance guard tests integer distance numerically:
`rounded = Round[Re[diff]]; rounded >= 0 &&
PossibleZeroQ[PChop[diff - rounded]]` (Recurrence.m:650-663), explicitly
replacing an earlier `IntegerQ` test that is always False on approximate
numbers (Docs/FeynmanTrickBoxFamilyStatus.md:109-111, disease D1).
WHY: `IntegerQ[2.0000000...]` is False; the resonant case then slips through
the guard into a division by a software zero, poisoning silently.
DISPOSITION: subsumed-by-exact-indicial-arithmetic (a, b are exact; resonance
is exact integer arithmetic in DiffExp2).  Keep as a review checklist item
for any residual numeric comparison.

### L13 FiniteAbsMax maps non-finite/non-numeric to 0 — a deliberate but dangerous default
WHAT: `FiniteAbsMax[expr]` returns 0 when expr contains Overflow,
Indeterminate, infinities, or unevaluated LinearSolve/LeastSquares heads
(Utilities.m:139-151, pattern at 129-130).
WHY: It exists so error ESTIMATION never crashes — but it means a completely
broken segment can report error 0.  Anywhere it guards a correctness
decision (not just a printout) it is a D1 silent-degradation site.
DISPOSITION: implement-in-Transport.m with inverted semantics: non-finite
content in an error probe is itself a loud error; the 0-default survives
only for cosmetic printing.

### L14 FT LaurentTrim's ABSOLUTE 1e-40 zero test
WHAT: `zeroCoeffQ[c] := PossibleZeroQ[c] || (NumericQ[c] && Abs[c] < 10^-40)`
advances a Laurent window's MinPower past "numerically zero" leading
coefficients (FeynmanTrick/DiffExpIntegration.m:1177-1190).
WHY: An absolute test misclassifies when the tower's natural scale is tiny
or huge; a wrong MinPower shifts the entire Laurent series AND the uniform
per-level prefactor shift built on it (L62), i.e. one bad trim re-labels
every eps order of every master below that level.
DISPOSITION: implement-in-Tolerances.m (laurentLeadTol, RELATIVE to the
tower's own magnitude) consumed by implement-in-FT-layer-hardening (Stage A1
edit of LaurentTrim); EpsSeries.m windows make the trim explicit-metadata
rather than inference.

---

## C. Order budgets, windows, error estimation

### L15 Coupling-depth t-order degradation (MaxCouplingOrder)
WHAT: `InitializeIntegrationSequence` builds the coupling graph of the eps^0
system, orders blocks, and sets `MaxCouplingOrder` = longest coupled chain
(MatrixLoading.m:358-385).  Downstream: the matrix-truncation error term is
read at order `ExpansionOrder - ISafetyExpansionSubtract -
(MaxCouplingOrder - 1)` (LineSegmentation.m:109-113, with
ISafetyExpansionSubtract = 5, State.m:216); the evaluation error probe
re-evaluates at order reduced by `ICurrEvalErrorSeriesDecrease :=
Ceiling[0.7*MaxCouplingOrder] + 2` (State.m:222, used Transport.m:913);
`IDecreaseOrderByErrorPrecise = MaxCouplingOrder` (State.m:223).
WHY: Each inhomogeneous solve in a chained block consumes top t-orders: the
top (couplingDepth-1) orders of a chained particular solution are degraded
garbage.  Error estimates and matching solves that trust them are wrong; the
0.7 factor + 2 is the empirically calibrated probe depth.
DISPOSITION: implement-in-Transport.m / implement-in-Solve.m as the TWindow
"CompleteMax" metadata of plan 3.1 (tracked, not guessed) — the formulas
above become the *defaults* for how much to discount, and the probe (L16)
verifies them.

### L16 Two-point full-vs-reduced-order error probe; avoid x = 0; per-indeterminate; abort > 1
WHAT: EstimateError = "Fast" evaluates each segment at full and
order-reduced truncation at BOTH the eval point and the matching point
(Transport.m:905-924), explicitly avoiding the point 0 "because there might
be logarithms, and the error should be manifestly zero" (Transport.m:915-923
— when BoundaryFixPoint == 0 it reuses the endpoint values as a small hack);
errors are computed per (integral, eps-order) and per symbolic indeterminate
(`ComputeErrorsPerIndeterminate`, Transport.m:947-962), accumulated
additively across segments (Transport.m:980-984), and the run aborts when a
segment error exceeds 1 (Transport.m:991-993).  `ErrorEstimates` is a
user-facing output consumed by tests (Transport.m:1238-1239).
WHY: The reduced-order re-evaluation is the only cheap truncation-error
estimator that sees the ACTUAL local solution; probing at x=0 is blind
because every Logx and positive power vanishes there; per-indeterminate
errors are required for symbolic-BC transport; the abort>1 guard catches
numerical instability before it propagates further.
DISPOSITION: implement-in-Transport.m ("ErrorEstimate" field of
LocalSolution + per-segment probe per plan 3.1/3.2) and
implement-in-API.m (ErrorEstimates output).

### L17 Incoming-error seeding from Accuracy, and its padding hazard
WHAT: TransportTo seeds accumulated errors from a third BC argument if
given, else from `10^-Accuracy` of the incoming boundary values
(Transport.m:578-587), then pads the per-order error array with
`10^-WP` entries when the eps dimension is short (Transport.m:589-591).
WHY: Chained transports (FT levels) must inherit upstream error; BUT the
`PadRight` with the *best possible* error (10^-WP) for unknown orders is a
zero-padding fallback of class D1 — unknown orders should inherit the WORST
known error or be flagged, not be assumed perfect.
DISPOSITION: implement-in-Transport.m (error chaining as explicit metadata);
the padding direction is fixed as part of A1 (unknown = pessimistic + named).

### L18 Predivision digit budgeting: + log10(#segments)
WHAT: Predivision dry-runs the entire segmentation without solving to count
segments (Transport.m:666-750), then sets
`DigitsNeeded = AccuracyGoal + Ceiling[Log10[SegmentsToIntegrate]] +
ISafetyDigits` with ISafetyDigits = 2 (Transport.m:758-760, State.m:215;
non-predivision baseline Transport.m:527, LineSegmentation.m:121).
WHY: Per-segment errors add; N segments need log10(N) extra digits per
segment to hit the global accuracy goal.  Long lines (banana to t = 10
through three thresholds) silently lose exactly these digits without it.
DISPOSITION: implement-in-Transport.m (two-pass predivision with digit
budget; plan 3.2 names it).

### L19 AccuracyGoalValidate "Before": adaptive expansion-order search
WHAT: Per segment, starting from the rolling mean of the last
`IExpansionOrdersAveraging = 3` used orders (Transport.m:780, State.m:217),
the matrix-truncation error term (L15's GetLargestTerm at the matching and
eval points, Transport.m:762-772) is compared to DigitsNeeded; the order is
raised in steps of `IExpansionOrderIncrease = 10` until sufficient, and on
first try lowered in steps of `IExpansionOrderDecrease = 10` while the
surplus exceeds `IDigitsSurplusDecreaseExpansionOrder = 3` digits, never
below `IMinExpansionOrder = 10` (Transport.m:776-841; State.m:218-224).
WHY: Expansion order is the dominant cost knob; fixed global orders either
waste large factors of runtime on easy segments or under-resolve hard ones.
The averaging window damps oscillation between adjacent segments; the
3-digit surplus hysteresis prevents thrashing.
DISPOSITION: implement-in-Transport.m (same control loop on the new chart
objects; constants live in Tolerances.m/Config.m).

### L20 AccuracyGoalValidate "After": redo segment at +25 with full state restore
WHAT: After integrating a segment, if the measured error exceeds
`10^-AccuracyGoal`, the expansion order is raised by
`IExpansionOrderIncrease2 = 25` (State.m:220) and the SEGMENT IS REDONE:
all loop state is restored from `CurrStatusBackup` (bcs, line, center,
accumulated errors), saved expansions are popped, the cached matrix
expansion is dropped, and the counter decremented (Transport.m:1191-1213,
backup at Transport.m:1011).
WHY: "Before" only bounds the matrix-truncation error; the solve itself can
lose more.  The redo path is the only self-healing mechanism, and the state
restore list documents exactly what constitutes transport state — a
non-obvious inventory that must survive in DiffExp2's marching loop.
DISPOSITION: implement-in-Transport.m (retry policy with immutable segment
inputs so "restore" is trivial by construction).

### L21 Single-level eps budgets underestimate deep chains — carry full depth
WHAT: The runner deliberately carries the FULL incoming boundary depth into
each level's transport instead of cutting to
`RequiredTransportEpsilonOrder`; the single-level estimate is kept only as a
lower-bound diagnostic (Scripts/run_ft_stepwise.m:155-167; pipeline-side
trim+warn at FeynmanTrick/DiffExpIntegration.m:1759-1777).
WHY: Each level's IBP eps-poles and prefactor shifts consume orders
CUMULATIVELY; cutting per level starved deep chains — box_bubble L2 was
capped at eps order 2 while needing 9 (and L1 needed 11) despite a generous
deepest-level budget (Docs/FeynmanTrickBoxFamilyStatus.md:19-26).  The plan
3.4 static budget MUST reproduce these measured numbers (M5 unit).
DISPOSITION: implement-in-FT-layer-hardening (static budget per plan 3.4
with cumulative prefactorShift convention) backstopped by
subsumed-by-EpsWindow-propagation (any miscount becomes a named error, never
a silent cut).

### L22 IntegrationPoleAllowance = 4 = 1 pole + 2 sector-solve + 1 truncation guard
WHAT: Integration steps add `IntegrationPoleAllowance` (config default 4,
env-overridable via FT_INTEGRATION_POLE_ALLOWANCE) extra eps orders to the
transport requirement (FeynmanTrick/DiffExpIntegration.m:1421-1436,
1455-1465); the decomposition of the 4 is documented in
Docs/FeynmanTrickBananaStatus.md:154-158: each endpoint pole deepens the
Laurent window by one order, solving sector weights at offset q consumes
offsets q..q+sectorCount-1, and the top offsets are the truncation boundary.
The deepest-level boundary is computed at `epsOrder + nLevels`
(FeynmanTrick/DiffExpIntegration.m:1717).
WHY: This is the empirically calibrated lookahead that made the banana L1
sector fits converge; with allowance 1 the fitter consumed partially
combined top orders as data (complex eps^-1 garbage).
DISPOSITION: subsumed-by-static-order-budget (plan 3.4: the "enhancement"
and "+1 safety" terms replace the fit-driven parts; sector-solve lookahead
dies with the fitter) — but the budget validation must dominate the recorded
allowance-4 behavior on banana.

### L23 Endpoint pole enhancement counts DEPTH max(p+1), not sector count
WHAT: Review-derived correction (legacy finding 5, math finding 7):
`Integrate[t^(-1+b eps) (eps Log t)^p]` produces `1/(b eps)^(p+1)` — a
log-chain endpoint pole enhances by p+1; several plain p=0 pole sectors
each contribute one inverse power that ADDS (max), not compounds.  The old
code never computed this number explicitly — it paid for it implicitly via
allowance-4 plus trims (L21/L22).
WHY: Both overcounting (several p=0 sectors) and undercounting (one p=2
sector = depth 3) corrupt the static budget; banana endpoints have p>0
in-scope from day one.
DISPOSITION: implement-in-Integrate.m (case table produces the actual depth)
feeding the plan 3.4 budget term `enhancement = max over endpoint pole
sectors of (p+1)`.

### L24 True-resonance homogeneous sectors have kmin = -p by construction
WHAT: Review-derived (math finding 14) + campaign data: under the
`(eps Logx)^p / p!` normalization, log-chain solutions of a genuine Jordan
block (banana L1 x = 1/2 chart, eigenvalues {0,0,0,0,1}, size-3 block —
Docs/FeynmanTrickBananaStatus.md "level-0" section) carry eps-INDEPENDENT
weights, i.e. coefficient windows starting at eps^(-p).
WHY: Window arithmetic that treats negative kmin in homogeneous solutions
as an error would reject correct solutions; this is the one systematic
EXEMPTION to "negative window = enhancement happened".
DISPOSITION: implement-in-EpsSeries.m (window arithmetic accounts for the
exemption; plan 3.1 invariant text).

### L25 Pseudo-resonance is in-scope from day one
WHAT: Campaign fact: the banana segment-1 endpoint mixes x^0, x^(-1+eps),
x^(2 eps) sectors (effective b values {2, 1} over a = -1;
Docs/FeynmanTrickBananaStatus.md:96-104) — integer-spaced a with DIFFERENT
b.  In the old per-eps-order representation this appears only as residual
Logx towers; symbolically the Frobenius denominator at the collision is
exactly `(b_i - b_j) eps`, a pure Laurent shift.
WHY: Treating these collisions as window losses silently degrades every
downstream order; the campaign's fitter existed largely to undo this.  The
plan's I2 joint-solve spec (solve integer-spaced-a families jointly,
choosing the eps-regular combination) is the discharge.
DISPOSITION: subsumed-by-I2-joint-solve (implement-in-Solve.m); M3 carries
the banana-derived pseudo-resonance closed-form unit.

---

## D. Segmentation and chart geometry

### L26 Complex singularities: ghost projection Re, Re±Im with suppression
WHAT: `FindMatrixSingularities` solves ALL irreducible factors (complex
roots included), sorts by real part, and projects each complex root onto the
real line as ghost points `Re - Im`, `Re`, `Re + Im` — where the side ghosts
are SUPPRESSED if another singularity already lies in the corresponding real
interval (LineSegmentation.m:65-106, projection at 81-99).  TransportTo
keeps the unprojected complex list as `imaginarySingularities` and treats
ghost-centered charts as regular (no branch handling), only announcing
"centered at singularity" for genuine ones (Transport.m:603, 1014-1016).
WHY: A complex singularity at distance d from the real line caps the
convergence radius of every nearby real chart; ignoring it produces
divergent partial sums with no real singularity in sight.  Pentagon and the
unequal-mass banana Kallen denominators put complex roots on generic lines.
The projection is the old code's approximation of "radius = distance in the
complex plane".
DISPOSITION: subsumed-by-complex-radius (plan 3.1: Radius = true complex
distance; segmentation solves exact roots incl. complex; ghost CHARTS
disappear, but the suppression lesson survives as "complex roots influence
radius only, never indicial structure") — implement-in-Transport.m.

### L27 DivisionOrder controls THREE things; the match point must satisfy BOTH adjacent charts
WHAT: k = DivisionOrder governs (a) the evaluation/matching point at chart
radius/k (Transport.m:646, 679, 1041); (b) next-center placement: GetCPL/
GetCPR solve for the next center such that the previous eval point ALSO sits
at the 1/k fraction of the NEW chart's radius (Mobius.m:98-142, Mobius
variants 72-90, used Transport.m:734-736, 1161-1163); (c) through the
segment count, the digit budget (L18).  Classic default k = 3 (State.m:112);
FT pins k = 4 (FeynmanTrick/DiffExpIntegration.m:272).
WHY: (b) is the actual defense against ill-conditioned weight matching: a
matching point deep inside BOTH charts keeps both local series
well-converged at the solve point (risk R3).  FT's k = 4 was needed because
integration consumes series accuracy across the whole segment, not just at
the match point.
DISPOSITION: implement-in-Transport.m (chart sizing + matching-point
geometry exactly as GetCPL/GetCPR; k in Config.m).

### L28 RadiusOfConvergence chart rescaling keeps high-order coefficients O(1)
WHAT: Every segment line is rescaled `x -> x/RadiusOfConvergence` at
construction (Mobius.m:47, 65); banana classic transport REQUIRES
RoC -> 10: "Without setting this option, the intermediate expansions blow
up" (Reference/Examples/Banana_example.m:75).
WHY: With chart radius r != 1, coefficient c_n ~ r^-n; at order 50+ this
overflows/underflows the working precision budget and poisons Pade and
matching.  Rescaling normalizes the geometric growth so coefficients stay
O(1).
DISPOSITION: implement-in-Transport.m (chart coordinate normalized so the
local radius is O(1) ALWAYS, automatically — not a user knob; keep the
config key accepted for compatibility, mapped or ignored with a note in
Config.md).

### L29 Mobius charts: needed for classic lines, forbidden for integration
WHAT: UseMobius -> True maps (prev sing, center, next sing) -> (-1, 0, +1)
(Mobius.m:22-35) doubling usable chart range; classic examples require it
(Banana_example.m:19, 76; FivePointNonPlanar_example.m:150;
MultiplePolylogarithms_example.m:56).  The FT integration layer hard-pins
`UseMobius -> False  (* Required for integration! *)`
(FeynmanTrick/DiffExpIntegration.m:352; RegularizedIntegration.m:8-10) —
the integration code assumes affine local-to-main maps (constant Jacobian,
L78).
WHY: Mobius is a 2x-or-better segment-count saver on long classic lines, but
the regularized-integration coordinate bookkeeping (bounds interpolation,
Jacobians at midpoints, prefactor power mapping) is only correct for affine
maps.  Silence here recreates exactly a seam-class bug.
DISPOSITION: implement-in-Transport.m (Mobius charts as optional chart map
for classic transport, M4 parity needs them for banana) AND
implement-in-Integrate.m (LOUD rejection of non-affine charts, plan 3.2).

### L30 Interior-pole approach: evaluation at the midpoint of the admissible window
WHAT: When the next singularity's pole interval overlaps the current chart,
the matching point is placed at the MIDPOINT of
`(current interval) ∩ (pole interval) ∩ (current center, pole center)`
(Transport.m:700-710, 1121-1139).
WHY: The midpoint balances convergence of the outgoing chart against the
incoming singular chart; biasing toward either edge makes one of the two
series marginal at the match point.
DISPOSITION: implement-in-Transport.m (with the L27 geometry; assert both
radii cover the point with margin).

### L31 Lines lying ON a singularity are rejected up front
WHAT: TransportTo aborts if any irreducible factor vanishes identically
along the line, naming the factors (Transport.m:572-574); boundary
points on singularities only warn when the asymptotic limit can be finite
(Transport.m:51-53).
WHY: Identically-singular lines produce 0/0 matrices after substitution;
detecting this late yields inscrutable solver failures.
DISPOSITION: implement-in-Transport.m (LoadSystem/segmentation precondition,
loud error naming factors).

### L32 Only affine line segments; RelateLines crosschecks numerically
WHAT: Non-linear segments are rejected (LineSegmentation.m:153-155, also
RelateLines' multiple-solution error 34-38); `RelateLines` solves one
component for the reparametrization, then verifies ALL components agree
numerically under PChop before accepting (LineSegmentation.m:20-59).
WHY: All downstream geometry (intervals, jacobians, GetCP*) assumes affine
maps; the full-component crosscheck catches lines that agree in one variable
but differ in another (a real failure mode with multi-variable kinematics).
DISPOSITION: implement-in-Transport.m (line algebra with the same
verify-all-components contract).

### L33 Dynamic segmentation strategy
WHAT: SegmentationStrategy -> "Dynamic" sizes each segment from the matrix
expansion error (`GetMatricesPrecisionDistance`, LineSegmentation.m:116-145,
root-solve for where the largest truncated term hits 10^-DigitsNeeded;
even/odd order root selection at 132-141).
WHY: Predecessor of predivision; survives for boundary cases but no example
in the pin suite uses it and the digit budgeting (L18) only exists for
Predivision.
DISPOSITION: waived-because-Predivision-only-v1 (plan Config table:
SegmentationStrategy accepts only "Predivision"; "Dynamic" = loud
not-supported error).

---

## E. Branch data and analytic continuation

### L34 Per-chart idelta sign is DERIVED from all vanishing prescription factors
WHAT: `PrepareAnalyticContinuation` collects every DeltaPrescription factor
whose leading series coefficient on the chart vanishes at x = 0
(AnalyticContinuation.m:21-27), and for each derives a required
`Sign[Im x]` as (prescribed sign) / (sign of the factor's leading
coefficient, with `Sign[x] -> 1`) PROVIDED the vanishing exponent in x is
exactly 1; any factor vanishing with exponent k != 1 yields "?"
(AnalyticContinuation.m:45-50).  Multiple distinct signs, or any "?", set
`AnalyticContinuationFailed` and default to +1 (AnalyticContinuation.m:55-68).
Results are cached per line (AnalyticContinuation.m:81-89).
WHY: The Im-sign of the chart variable is not free data — it is forced
jointly by every physical prescription that vanishes there; storing "one
factor, one sign" (as an early plan draft did) cannot represent conflicts or
multiplicity.  The reviewer-amended design upgrades "?": EVEN multiplicity
genuinely imposes no constraint (the factor's phase winds fully), while the
old code conservatively flags it failed.
DISPOSITION: implement-in-Transport.m as the plan 3.1 "Prescriptions" record:
LIST of (factor, sign, multiplicity, leading-coeff sign) -> consistency-
checked single Im-sign at chart construction; even multiplicity = no
constraint; conflict or missing prescription at a chart with b != 0 / p > 0
sectors = LOUD ERROR.

### L35 The crossing convention: +1 is a NO-OP; -1 is a Logx SHIFT, not a phase
WHAT: For derived sign +1 NO replacement is applied at all — the code relies
on principal-branch `Log`/`Power` at negative argument (this PINS which of
±i pi appears).  For sign -1 the replacements are:
`Logx -> (thetap + thetam) Logx - 2 pi I thetam` (a SHIFT of the log);
half-integer powers via `x^b -> x^(b-1/2) (thetap - thetam) Sqrt[x]`;
and `x^b -> (thetap + Exp[-2 pi I b] thetam) x^b` ONLY for
`Denominator[b] > 2` (AnalyticContinuation.m:70-79).  Theta symbols are
resolved per evaluation point sign (Pade.m:70-77, SEval/SEval2).
CONSEQUENCE (math review finding 5): because the -1 rule SHIFTS Logx, a
sector with (eps Logx)^p mixes BINOMIALLY into all lower-p members of its
confluent family under crossing — the crossing operator is
phase x unipotent matrix `M_{p->p-j} = (±i pi eps)^j / j!`, with the phase
`e^(±i pi (a + b eps))` (non-integer a occurs; the sqrt rule above is its
half-integer special case).
WHY: Getting the convention wrong flips imaginary parts of every result past
a branch point; the old form is the principal-branch-on-the-far-side
convention, NOT the |t| convention.  Both M4 parity and the 9aeb300
interior-split fix (L37) depend on knowing exactly which convention the
coefficients carry.
DISPOSITION: implement-in-Transport.m (explicit crossing operator =
phase x log-chain mixing matrix, convention pinned to the old
principal-branch-far-side form; plan 3.2).

### L36 Square-root factors are auto-prescribed +idelta and flipped in the matrices
WHAT: At load, irreducible square-root arguments found in the matrices
(`DEqnSquareRoots`, MatrixLoading.m:167-179) that the user did not prescribe
are added as {root, +1} prescriptions (`SquareRootPrescriptionsAdded`,
State.m:227-230; MatrixLoading.m:181-190), and all matrix sqrt's are
rewritten so the branch matches the prescription:
`Sqrt[-sigma] -> -I Sqrt[sigma]` per sign (MatrixLoading.m:193-212).
`PrepareAnalyticContinuation` separately tracks whether the current chart's
prescription came ONLY from an auto-added sqrt
(`CurrentSingularityWasAddedFromSquareRoot`, AnalyticContinuation.m:29-43),
which tightens the multivalued check (L38).
WHY: Sqrt branch points are singularities of the SOLUTIONS even where the
matrix is finite; without auto-prescription, transport across them silently
picks an arbitrary branch.  The matrix flip keeps the DE coefficients
consistent with the prescribed sheet.
DISPOSITION: waived-because-sqrt-x-dependence-out-of-scope-v1 (plan
non-goal: irrational x-dependence = loud LoadSystem error) — BUT the
auto-prescription RULE is recorded for v1.1 and the
"prescription provenance" distinction is kept in the Prescriptions record
(L34) for algebraic-exponent charts.

### L37 Interior-split real-log convention (commit 9aeb300)
WHAT: When an integration interval straddles the local singular point, it is
split at zero and BOTH half-segments are integrated under
`Block[{$InteriorSplitRealLog = True}, ...]`
(RegularizedIntegration.m:504-530); in that mode the negative arm evaluates
`Logx -> Log[Abs[pt]]` (REAL log of the distance) instead of complex
`Log[pt]` (RegularizedIntegration.m:1836-1843).  A loud warning fires if the
straddled crossing has b != 0 on either branch (the real-log split assumes a
meromorphic crossing; RegularizedIntegration.m:510-521).  The Hadamard
finite-part path (L72) makes the formal divergences cancel between the two
arms (PV for the residue).
WHY: Both arms hold real series of a meromorphic crossing; complex
Log[negative] leaks i*pi at one order and -pi^2 at the next into a real
result.  The endpoint-regularized paths, by contrast, RELY on complex Log
pairing with the theta-resolved branch phases — so the convention is
conditional, not global.  First exercised by the box L1 IBP pole at
xx1 = 1/4 (Docs/FeynmanTrickBoxFamilyStatus.md:48-55).
WHY-PAIRING: the cancellation only happens if the two half-segments are
actually paired; integrating one arm alone is wrong.
DISPOSITION: implement-in-Integrate.m (interior-pole PV/idelta with the
pairing ENFORCED by assertion, plan 3.2; real-log mode as the meromorphic
crossing rule; loud error for b != 0 straddles until specced).

### L38 SingularityCheck (default-on) and the multivalued-result guard
WHAT: With flag "SingularityCheck" active (CrosscheckFlags level 0 =>
ALWAYS on by default, State.m:196-207 with CrosscheckLevel default 0,
DiffExp/DiffExp.m:115-120), IntegrateSystem inspects the fixed solutions:
if Logx or fractional powers survive (chop at 10^-ICheckMultivaluedChop =
10^-5, State.m:212) at a chart with NO idelta prescription — or with only an
auto-sqrt prescription but Logx present — it raises the "current point is
not recognized as a branch point ... add DeltaPrescriptions" error
(Transport.m:470-492).  `AbortOnAnalyticContinuationFail` chooses warn vs
abort (Transport.m:471-477); on warn, `MultivaluedFail` aborts the transport
loop and returns the last line's data (Transport.m:899-943, 1185-1188).
WHY: This is the configuration-gap detector: multivalued local content with
no prescription means the user's DeltaPrescriptions are incomplete and ALL
downstream continuation is undefined.  The pentagon's standing warning is
exactly this (config gap, not solver bug — plan section 2 / M0 task 16).
DISPOSITION: implement-in-Transport.m (chart construction: b != 0 or p > 0
sectors at a chart with no consistent prescription = LOUD ERROR by default;
the FT layer's `AbortOnAnalyticContinuationFail -> False` accommodation
becomes an explicit per-chart override recorded in the output, not a global
silencer).

### L39 One-sided-validity semantics when continuation fails at the start
WHAT: If sign derivation fails on the FIRST segment of a transport (or at a
backward transport start), results are declared valid only on the side
(x > 0 or x < 0) where the BCs were given (Transport.m:494-496); on any
LATER segment it is a hard error telling the user to separate singularities
(Transport.m:1026-1032).
WHY: A one-sided expansion is still useful at the starting chart (the other
side never gets evaluated), but mid-line it means crossing an unprescribed
cut — silent wrongness.
DISPOSITION: implement-in-Transport.m (same two-tier rule; the "valid side"
becomes chart metadata instead of a printed warning).

### L40 Theta resolution happens per evaluation point; Logx stays symbolic
WHAT: FT comment block: "Do NOT substitute Logx -> Log[x]! ... Logx is
resolved only at the final boundary evaluation step"
(FeynmanTrick/DiffExpIntegration.m:52-62); theta rules are chosen from the
sign of the actual local point (or an explicit direction from the segment's
other bound), never once per segment
(FeynmanTrick/DiffExpIntegration.m:107-148, 125-138;
RegularizedIntegration.m:1772-1787).
WHY: A segment can straddle local x = 0 (interior pole) so the two arms need
opposite theta branches; substituting Log[x] early evaluates `0*Log[0] =
Indeterminate` at endpoints and destroys the structural Logx-power scans.
DISPOSITION: subsumed-by-exact-sector-tags (log content and branch data are
explicit tags; evaluation applies the branch rule at the point) —
implement-in-SectorSeries.m (evaluate with branch rule).

---

## F. Solver mechanics (regular and singular recurrences)

### L41 Denominator-clearing polynomial recurrence is the performance backbone
WHAT: Both the regular and singular recurrences first try "rational mode":
Together the factored block matrix, take the polynomial LCM denominator D(x)
(factoring x out first in the singular case), and run the recursion on
polynomial coefficient data so each step costs O(deg) matrix-vector products
instead of O(n) (RationalizeAMatrixCore, Recurrence.m:13-77; fundamental
matrix Recurrence.m:241-273, 475-545; particular Recurrence.m:697-748).
Guards: the attempt is wrapped in `TimeConstrained[..., 5.0, {False}]`
(Recurrence.m:25-74) and is rejected as not beneficial when the cleared
numerator degree exceeds ExpansionOrder/2 (Recurrence.m:66-68); fallback is
direct series-coefficient convolution (Recurrence.m:234-238, 302-307).
D(0) != 0 is required (Recurrence.m:57, 61-63).
WHY: For long classic lines this is the constant-factor difference between
minutes and hours (the R6/M4 2x benchmark has no other mechanism); the
degree gate and timeout exist because pathological rational matrices make
clearing slower than convolving.
DISPOSITION: implement-in-Solve.m (plan 3.2: x-denominators of A cleared up
front so recursion coefficients are polynomial; keep BOTH the cost gate and
the loud fallback to convolution mode).

### L42 Per-segment block caching across eps orders
WHAT: IntegrateSystem keeps `segmentCaches[intind]` across the eps-order
loop and threads it through DispatchStrategy (Transport.m:276-285); each
strategy computes its homogeneous data once at epsord 0 and reuses it
("FMat"/"VOP" Default.m:29-175, VOP.m:22-74, 185-282; "RR"
Recurrence.m:292-315; "SingRR" Recurrence.m:804-843; "GenSingRR"
ResonantRecurrence.m:1307-1337; "FuchsianRR" with nested Subcache
ResonantRecurrence.m:1463-1496).  Dispatch also short-circuits applicability
re-checks when the cache key exists (Dispatch.m:30-41).
WHY: The homogeneous fundamental system is eps-order independent in the
per-order formulation; recomputing it per order multiplies segment cost by
(EpsilonOrder+1).
DISPOSITION: subsumed-by-symbolic-eps-solving (DiffExp2 solves all eps
orders at once per chart, so the cache dimension disappears) — but
implement-in-Transport.m for the cross-CHART reuse aspect: chart-level
solver artifacts (indicial data, cleared-denominator data) are computed once
per chart and shared by transport and integration consumers.

### L43 Logx powers are solved TOP-DOWN in the particular recursion
WHAT: The particular recurrence iterates log power k from highest to 0
inside each x-order m, because d/dx Logx^(k+1) contributes to the Logx^k
equation at the same power of x (Recurrence.m:344-393, esp. comment 345-347
and the `logDerivativeCoeff` term 368-378; same ordering
ResonantRecurrence.m:969-976).
WHY: Solving bottom-up uses not-yet-computed higher-log data (silently zero)
— a wrong-result bug, not a crash.
DISPOSITION: implement-in-Solve.m (log-chain descent ordering; in the
symbolic-eps formulation this becomes the explicit lower-triangular
log-chain coupling, same ordering requirement).

### L44 Singular recurrence: diagonalize in the residue eigenbasis, reject degenerate snapped spectra
WHAT: `PrepareSingularRecurrence` requires min order exactly -1, nonzero
residue, eigenvalues snapped at 10^(-CP/2), then REJECTS the block if the
snapped spectrum has any multiplicity (Recurrence.m:103-172, degenerate
rejection 137-146 with the comment: float noise splits degenerate
eigenvalues by ~sqrt(noise), so the eigenvector matrix looks full-rank while
being catastrophically ill-conditioned), checks eigenvector rank and
non-resonance (no positive-integer differences; 148-160), and routes
rejected blocks to the general solver.
WHY: Banana L1 blocks were once "solved" by plain diagonalization with a
1e111-amplified basis (Docs/FeynmanTrickBananaStatus.md, fix 3) — the rank
test cannot see near-parallel eigenvectors.  sqrt(noise) splitting is the
signature failure of numeric spectra of structurally degenerate matrices.
DISPOSITION: subsumed-by-exact-indicial-algebra (Indicial.m factors the
characteristic polynomial exactly; degeneracy is exact data) — record the
sqrt(noise) phenomenon in Indicial.m's spec for the numeric-validation path.

### L45 Normalize source TYPES before structural decisions (inert SeriesCoefficient tails)
WHAT: `ComputeSingularParticular` rewrites source entries whose head is not
SeriesData: inert `SeriesCoefficient[ser, k]` requests AT OR BEYOND the
series window are replaced by 0 and the entry re-SExpand-ed; a non-series
survivor triggers a loud warning (Recurrence.m:558-585).
WHY: Upstream series arithmetic leaves one-past-the-end SeriesCoefficient
tail requests (numerically ~1e-25); their SYMBOLIC presence turns the
entry's head into Plus, and the old `Head === SeriesData` probe then
silently fell back to leading power 0 — dropping genuine x^-1 source content
(box eps^0 deficit chain, Docs/FeynmanTrickBoxFamilyStatus.md:90-103: "their
TYPE was the poison, not their size").  Disease D1: type poisoning.
DISPOSITION: implement-in-EpsSeries.m / implement-in-SectorSeries.m by
construction (typed coefficient arrays cannot carry inert symbolic heads;
constructors validate types loudly).

### L46 Never divide by a software zero: zero-divisor guards defer to the general solver
WHAT: Both particular recursions guard every divisor with
`PossibleZeroQ[PChop[divisor]]` and `Throw` a resonance signal instead of
dividing (Recurrence.m:714-723, 734-741, caught at 695-747 -> `$Failed` ->
SolveSingularRecurrence falls through to SolveGeneralSingularRecurrence,
Recurrence.m:855-861).
WHY: A divisor that should be nonzero by the resonance guard can still be a
software zero after noise; dividing poisons silently and the corruption
surfaces segments later.
DISPOSITION: implement-in-Solve.m — in symbolic eps the analogue is: a
denominator `(n + a_i - a_j + (b_i - b_j) eps)` that vanishes IDENTICALLY
must route to the log-chain/joint-solve path, and a numerically-zero
eps^0 part must route to the Laurent-shift path; no silent division ever.

### L47 Normalize the particular back to SeriesData on return
WHAT: After multiplying by x^s, each particular component is re-SExpand-ed
and warned about if its head is not SeriesData (Recurrence.m:770-787).
WHY: Downstream boundary fixing and per-order source assembly only reliably
consume SeriesData; head-Plus particulars were silently dropped (box status
doc item 5, line 112-114).
DISPOSITION: subsumed-by-typed-representation (same as L45).

### L48 Boundary constants are extracted by DIFFERENTIATION, not Coefficient
WHAT: When fixing constants at FixAt != 0, fGeneral's dependence on each c_i
is extracted as `D[term, c_i] /. (all c -> 0)` and evaluated with SEval,
plus the c-free part separately (Transport.m:311-327, comment 311-316).
WHY: Recurrence-generated singular solutions are compound expressions with
SeriesData whose COEFFICIENTS contain the c_i; outer-level `Coefficient`
misses those nested constants, and sending symbolic constants into
PadeApproximant breaks it.  fGeneral is linear in c_i so differentiation is
exact.
DISPOSITION: implement-in-Transport.m (matching solve assembles the
fundamental-matrix columns explicitly from the solver output — by
construction in DiffExp2, but the parity-path code that evaluates legacy
expressions must keep this rule).

### L49 Boundary matching: LinearSolve with chop ZeroTest, then CHECKED least-squares
WHAT: The boundary solve uses
`LinearSolve[..., ZeroTest -> (N[LSPChop@..., LSCP] == 0 &)]`
(Transport.m:341); on failure (or non-finite result) it tries LeastSquares
and accepts ONLY if the relative residual passes
`lsTol = 10^-Min[80, Max[20, Floor[LSCP/4]]]`
(Transport.m:348-374); on final failure the full context is stored in
`LastErrorContext` and a loud error names the integrals (Transport.m:376-388).
WHY: High-precision rank decisions need the chop-based ZeroTest (default
exact ZeroTest stalls or misjudges); raw LeastSquares would silently project
inconsistent systems — the residual gate makes it a checked fallback, not a
silent one.
DISPOSITION: implement-in-Transport.m (matching = eps-graded Laurent solve
per plan 3.2, with the same checked-fallback discipline; tolerances from
Tolerances.m).

### L50 Underdetermined boundary data introduces NAMED free parameters and disables Pade
WHAT: If the boundary matrix has a nullspace, the solution gains explicit
symbols `Subscript[c, epsord, intind, i]` along the null vectors, a warning
names them, the full matching context is stored, and Pade is turned off
(`TurnOffPade`, since approximants cannot carry symbols)
(Transport.m:392-434, 178-181); same when no boundary terms exist at all
(Transport.m:439-441).  Indeterminate coefficients in INPUT BCs likewise
disable Pade with a notice (Transport.m:565-569).
WHY: Partially-known BCs are a first-class workflow ("?" wildcards, symbolic
transport); silently zeroing the free directions would be D1.  The
per-indeterminate error machinery (L16) exists to support exactly this.
DISPOSITION: implement-in-API.m / implement-in-Transport.m (symbolic
indeterminate support is in the plan's API list; Pade exclusion rule kept).

### L51 Logx-bearing Wronskians need the Frobenius-method inverse
WHAT: InvWronskSolver "Auto" picks: Logx present in the Wronskian =>
"Frobenius" (solve the ADJOINT-derived nth-order equation for the inverse
basis, Default.m:70-131), else plain "Inverse" (Default.m:133-135); the
Frobenius route normalizes `WronskInvPrime . Wronsk` by its CONSTANT part
(`/. x^k -> 0 /. Logx -> 0`, Default.m:128) before inverting, with the
identity crosscheck behind flag "WronskInv" (Default.m:121-126).  A direct
series `Inverse` on Logx data is wrong (Logx is not algebraically
independent of x under d/dx); `MatrixLogxInverse` (Wronskian.m:17-32) exists
as the order-by-order alternative.
WHY: Singular charts produce log solutions; inverting their Wronskian
naively was an early wrong-result source.  The constant-normalization step
encodes that the product is I + (nilpotent series), invertible by its
constant term.
DISPOSITION: subsumed-by-one-solver (DiffExp2 never inverts log-bearing
fundamental matrices — matching solves use the eps-graded linear solve on
coefficient data, and variation-of-parameters disappears with VOP.m).

### L52 Pivot selection with fallback (Mtilde singularity)
WHAT: `CombineDifferentialEquationsWithPivotSelection` tries every pivot
integral until the derived nth-order ODE has a 1-dimensional nullspace /
invertible Mtilde, falling back to the VOPAlt strategy when all pivots fail
(Wronskian.m:48-108; Default.m:36-50); Mtilde inversion itself has a
DivisionFreeRowReduction-with-Check fallback (Default.m:143-152).
WHY: Reducing a first-order system to one nth-order scalar ODE fails when
the chosen component decouples; pivoting is the cheap fix, VOPAlt the
expensive one.  Encodes: NEVER assume component 1 is cyclic.
DISPOSITION: subsumed-by-one-solver (no scalar-ODE reduction in DiffExp2's
direct vector Frobenius recursion).

### L53 HomogeneousSolve "Expand" vs "DontExpand"; FT pins "DontExpand"
WHAT: Two modes thread through the strategy stack: numeric series arithmetic
("Expand") vs exact `Together` on the factored matrix ("DontExpand")
(Default.m:38-42, Helpers.m:9-18, Wronskian.m:53-67); NullSpace failures on
expanded data print "Try setting HomogeneousSolve -> DontExpand"
(Wronskian.m:35-45).  FT TransportLevel defaults to "DontExpand"
(FeynmanTrick/DiffExpIntegration.m:277).
WHY: Numeric nullspaces of series matrices can produce empty/garbage results
(negative-order empty SeriesData detected at Wronskian.m:38); exact
rational manipulation is slower but robust — the FT pipeline chose
robustness.
DISPOSITION: subsumed-by-exact-eps-rational-input (DiffExp2's recursion
coefficients come from the exact matrix, cleared of denominators; no
expanded-vs-exact mode split) — record in Config.md as accepted+ignored.

### L54 Crosscheck flag system: leveled, ODE-residual based
WHAT: `CrosscheckFlags` assigns each internal verification a level
(FrobeniusSolutions 1, MatrixDelta 1, Wronskians 1, WronskInv 0,
PeriodMatrix 1, GeneralSolutionMatrix 2, GeneralSolution 1,
VariationOfParameters 1, SingularityCheck 0; State.m:196-207); active set =
flags with level <= CrosscheckLevel plus any explicitly named
(DiffExp/DiffExp.m:115-120).  The checks substitute the candidate back into
the ODE and compare at CPChop = 10^-30 over the first
ICrossCheckVerifyResultOrder = 5 orders (Default.m:98-126, 157-167, 182-192;
Transport.m:443-452; Frobenius.m:110-134).
WHY: Cheap residual spot-checks at low order catch most numerical
catastrophes for ~zero cost; the leveling lets heavy checks stay off by
default.  Plan A2 makes a subset always-on.
DISPOSITION: implement-in-Solve.m / implement-in-Transport.m (plan 3.1
invariant: ODE residual spot-check at a random interior point per chart,
ALWAYS on; the leveled system collapses to always-on-cheap + debug-only-
expensive).

### L55 Indicial roots of degree > 2 are suspect
WHAT: Frobenius1 warns when the largest indicial root has
`Denominator[rMax] > 2` unless `"IgnoreIndicialCheck" -> True`
(Frobenius.m:42-45; config State.m:116); matrix loading separately REJECTS
higher-than-square roots in the input matrices (MatrixLoading.m:159-164).
WHY: Cube-and-higher roots were never exercised by the validated examples —
the warning marks the tested boundary of the algebra (half-integer powers
have dedicated handling everywhere: crossing rules, sqrt machinery).
DISPOSITION: implement-in-Indicial.m (the I1 char-poly factorization
contract subsumes the check: any exact algebraic a, b is now REPRESENTABLE,
but a denominator > 2 still deserves an INFO note since no pin covers it;
R2 flags performance not correctness).

### L56 Frobenius instability path: revalidate, then continue or abort
WHAT: If the Frobenius linear solve emits messages, the code distinguishes
hard failure (unevaluated LinearSolve head -> error) from possible
instability: it substitutes the candidate back into the ODE and continues
only if the residual is exactly zero, else aborts (Frobenius.m:58-79).
WHY: At high WP, LinearSolve warnings are often spurious; recomputing the
residual converts "maybe wrong" into a definite verdict instead of either
ignoring messages (D1) or failing good runs.
DISPOSITION: subsumed-by-always-on-residual-checks (L54 discharge); the
verify-don't-trust-messages pattern is recorded for any remaining
Mathematica-solver call sites.

### L57 Log-depth auto-extension and the NormalizeLogPower export trap (commit ebc4724)
WHAT: `DiffExpIntegrate` scans its input for the max Logx power and extends
the integration rule table when it exceeds `IMaxLogOrder`
(Integration.m:23-41, UpdateIntReps 62-74; also extended preemptively by the
resonant solvers, ResonantRecurrence.m:1256-1259, 1324-1329, and the
rational recurrence, Recurrence.m:330-332).  The scan rationalizes numeric
powers via `NormalizeLogPower` (SeriesOps.m:95-106).  THE BUG:
NormalizeLogPower had no ::usage export, so Integration.m's fully qualified
cross-package call silently returned UNEVALUATED, the integer filter saw no
powers, the auto-extension never fired, and `x^-1 Logx^k` (k >= 2) terms
were integrated AS CONSTANTS by the `b*Const` fallback — corruption starting
exactly at eps order 3 (Docs/FeynmanTrickBoxFamilyStatus.md:57-73; memory
note "Wolfram package context traps").
WHY: Two lessons: (1) capacity that grows on demand via pattern detection is
D1 — when detection fails, work is silently mangled rather than refused;
(2) Wolfram cross-package calls to unexported symbols no-op silently.
DISPOSITION: subsumed-by-closed-form-integration (Integrate.m's case table
has no capacity to extend — every (a, b, p, n) is handled exactly or errors)
+ implement-in-process: DiffExp2 package layout keeps ALL cross-module
symbols exported and a load-time assertion checks for `Symbol::shdw`-class
no-ops.

### L58 The IntReps pattern-table integrator and its b*Const fallback
WHAT: Series integration is done by replacement rules mapping
`Log[x]^n x^m -> antiderivative` built up to IMaxLogOrder
(Integration.m:44-74); any term matched by NO rule falls through to the
`b*Const` constant-integration path (Integration.m:49-52: `Out - Const +
b Const`), i.e. unmatched structures are integrated as constants.
WHY: This is the D1 mechanism behind L57's corruption; a closed-form
dispatcher must REFUSE structures it does not recognize.
DISPOSITION: subsumed-by-Integrate.m-case-table (exhaustive {b}x{a+n+1}x{p}
dispatch, loud error on anything else).

---

## G. Higher-order poles and resonance machinery

### L59 Local Fuchsianization must run in EXACT arithmetic with a degeneracy guard
WHAT: `BuildFuchsianizedRecurrenceData` rationalizes the local theta-matrix
(`Rationalize[..., 10^(-Max[20, WP-50])]`) BEFORE calling FuchsianizeLocal /
TrimFuchsianLattice (ResonantRecurrence.m:144-167), and rejects the result
if any lattice entry's x-degree exceeds `4*systemSize + 8`
(ResonantRecurrence.m:179-205).  The reduction itself: lattice saturation
adjoining minimum-valuation columns (LocalSeries.m:133-216, MaxSteps 200
with a loud non-termination Failure) and a trim pass removing unnecessary
column poles (LocalSeries.m:218-253).
WHY: Exact pole-cancellation tests NEVER fire on floats, so the trim pass
kept "succeeding" on noise and ratcheted lattice columns by x once per pass
up to x^50-scale entries while precision collapsed to ~15 digits (banana L0
fix 1, Docs/FeynmanTrickBananaStatus.md).  A local balance only needs power
shifts of order systemSize + pole depth — anything larger means the
reduction went numerically astray.
DISPOSITION: implement-in-Indicial.m (ported rank reduction per plan 3.2:
exact eps-rational input by construction, MaxSteps loud error, degree guard
kept as an assertion; banana L1 nilpotent double pole is the M2/M3 pin).

### L60 Degenerate spectra need toleranced generalized-eigenvector chains (NumericJordanData)
WHAT: When the snapped spectrum has multiplicities,
`ComputeResonanceStructure` builds Jordan data via SVD-based toleranced
nested null spaces of (M - lambda)^k with explicit chain assembly
(NumericJordanData, ResonantRecurrence.m:219-351), because (a) exact
JordanDecomposition on rationalized noisy residues sees noise-split
"distinct" eigenvalues with ill-scaled eigenvectors, and (b)
`NullSpace[..., Tolerance]` is unreliable for high-precision input (it can
underflow to machine arithmetic and return nothing — comment at 237-240).
Fallback ladder: numeric chains -> exact JordanDecomposition on rationalized
residue -> numeric JordanDecomposition -> loud error
(ResonantRecurrence.m:390-419).
WHY: This recovered the genuine size-3 Jordan block (log^2 solutions) at the
banana x = 1/2 chart that the old path missed entirely (banana doc fix 2).
DISPOSITION: subsumed-by-exact-indicial-algebra (Jordan structure of
A_{-1}(eps) computed exactly in Indicial.m under the I1 contract); the
NullSpace-Tolerance unreliability note is kept for any numeric validation
code.

### L61 Max log power = (position in block - 1) + sum of block sizes of LARGER eigenvalues in the resonance class
WHAT: For each solution, `maxLogPowers` = Jordan-chain position plus the
total size of all Jordan blocks of eigenvalues strictly greater (by integer
spacing) in its resonance class (ResonantRecurrence.m:440-508, esp.
489-494); resonance classes group eigenvalues differing by integers,
sorted smallest-first (440-460).
WHY: This is the exact combinatorial bound on log depth for resonant
Frobenius solutions — too small loses solutions, too large wastes a
multiplicative cost factor.
DISPOSITION: implement-in-Indicial.m (the {a, b, p} spec generator must
reproduce this rule; with exact algebra it is exact).

### L62 n = 0 Jordan-chain descent: resolve kernel freedom from the NEXT equation's solvability
WHAT: At n = 0 (L_0 singular by construction) the chain f_{0,k} for
k = initK-1 .. 0 is solved with the same kernel-freedom resolution as the
n >= 1 resonant steps: when a step is singular and the previous step left
nullspace freedom, the overlap system `leftNull . prevNull` determines the
free coefficients so the next equation becomes consistent
(ResonantRecurrence.m:684-733; the n >= 1 analogue 761-830).
WHY: The minimal-norm pseudo-inverse choice strips kernel components that
the true generalized-eigenvector chain NEEDS, making the next chain equation
inconsistent — which SolveRecurrenceStep would silently project away (banana
doc fix 4).  General lesson: pseudo-inverse "solutions" of singular chains
are wrong unless freedom is resolved forward.
DISPOSITION: subsumed-by-explicit-log-chains (true-resonance log-chain
construction in Solve.m is closed-form, no kernel guessing); the lesson
binds any remaining singular linear solves (assert consistency, never
project silently).

### L63 Singular steps must be solved as ONE coupled block across all log powers
WHAT: At a singular step n the k-by-k descent with PseudoInverse FAILS for
Jordan blocks of size > 1 (the null/left-null overlap is zero), so
`SolveSingularNBlock` assembles the full (kMax+1)*size system with L_n on
the diagonal and (k+1) I on the superdiagonal and solves it simultaneously
(ResonantRecurrence.m:884-935, dispatch at 951-977).
WHY: The log-chain coupling IS the system; solving levels independently and
patching with pseudo-inverses cannot represent it.
DISPOSITION: implement-in-Solve.m (the true-resonance construction solves
the coupled family jointly by design; this item is its numerical-era proof
of necessity).

### L64 Residual-driven log-depth growth, capped, with EVERY-step residual checks
WHAT: `ComputeUnifiedParticular` starts with
`kMax = bMaxLogK + Max[MaxLogPowers] + 1` (ResonantRecurrence.m:1191-1195),
runs the recurrence, checks the residual at EVERY step — not only singular
ones — because a non-singular step whose solve silently degraded is
otherwise invisible (CheckParticularResidual, ResonantRecurrence.m:985-1023,
comment 992-996), and grows kMax by 1 up to 3 times if the residual is
PROVEN above tolerance (1216-1254, L5 semantics), else returns $Failed
loudly (refusing fallback, 1356-1368).
WHY: The needed log depth of a resonant particular is not knowable a priori
in the numeric formulation; growth-with-verification was the safe protocol.
The "check every step" rule comes from the banana x = 1/2 eps^-1 impurity
(a dropped source at a NON-singular step).
DISPOSITION: subsumed-by-resonant-source-log-bump-rule (L65: in exact
algebra the needed p-bump is computable, source tag hits homogeneous tag =>
p -> p+1) — keep the full-residual verification as the M3 closed-form test
oracle.

### L65 The resonant-source empty-particular hole
WHAT: The general singular solver returns an EMPTY particular (with complex
leaks at higher orders) for resonant-with-source inputs — currently
bypassed, never properly fixed (Docs/FeynmanTrickBoxFamilyStatus.md:139-143;
the L64 machinery papers over the common cases).  No old-code parity oracle
exists for this path.
WHY: When the source tag satisfies a_s + n = a_i with b_s = b_i the
denominator is identically zero in eps and the ansatz NEEDS a log bump
(p -> p+1) — the box L1 endpoint sources "gain one Logx power per epsilon
order".  This is the one known correctness hole that survives in the frozen
oracle.
DISPOSITION: implement-in-Solve.m (explicit resonant-source log-augmentation
rule, plan 3.2) with an M3 closed-form pin (Integral of x^-1 log^k x class)
since no oracle exists.

### L66 Trailing non-finite series tails are trimmed; interior non-finites are fatal
WHAT: `TrimNonFiniteSeriesTails` drops a CONTIGUOUS trailing run of
non-finite coefficients from solver output, recording the dropped count and
first dropped power in a warning; non-finite coefficients NOT forming a pure
tail are a hard error with full context (ResonantRecurrence.m:62-120,
applied 1369-1378).
WHY: Overflow at the truncation boundary is benign (the orders were garbage
anyway, and the TWindow shrinks); overflow mid-series means the solve is
broken.  The distinction prevents both silent corruption and false alarms.
DISPOSITION: implement-in-Solve.m (same tail-vs-interior rule, feeding
TWindow CompleteMax instead of a warning).

---

## H. Endpoint sectors, drop rules, regularized integration

### L67 DecomposeSingularity's b-extraction guards (keep collapsed towers explicit)
WHAT: b is inferred from the tower as `b = [Logx^1 coeff at eps^(n+1)] /
[Logx^0 coeff at eps^n]` at the leading power, then: tiny imaginary part
stripped, rationalized at RationalizationTolerance, and REJECTED back to
b = 0 (towers kept explicit) if it does not rationalize cleanly —
denominator > 16, |b| > 100, or complex (SingularityDecomposition.m:329-365,
guard comment 351-357: a non-clean b is the collapsed AVERAGE of several
sectors and extracting it "would only reshuffle the epsilon tower against a
meaningless exponent").  `MultiplyByXMinusBEps` implements the x^(-b eps)
removal as an explicit eps-order convolution with (-b Logx)^k/k!
(SingularityDecomposition.m:225-243).  maxIter = 20 safety
(SingularityDecomposition.m:311).
WHY: Multi-sector endpoints (several x^(a + b_i eps) sharing integer a)
CANNOT be described by one exponent; an averaged or complex b poisons every
downstream drop/keep decision.  Keeping b = 0 with explicit towers at least
hands the recovery problem to the fitter intact (D2 made survivable).
DISPOSITION: subsumed-by-exact-sector-tags (DiffExp2 never reconstructs b
from towers; the indicial data IS the b-spectrum) — the FT call site
(FeynmanTrick/DiffExpIntegration.m:822-class) gets the plan's named
SectorSeries replacement API.

### L68 The collapsed-exponent disease and the residual-sector fitter (D2 core)
WHAT: `FitResidualEndpointSectors` reconstructs residual x^(r eps)
(eps Logx)^p sectors from the eps-tower of a fixed local power via shifted
Prony/Hankel moment diagonals `m_k = k! T(q0+k, k)`, scanning the reference
order q0 because leading weight slices can cancel (s0 = 0, the banana
{1,0,2,1} mechanism), with confluent (repeated-root) sectors as polynomial-
in-k weights (RegularizedIntegration.m:829-1331; header derivation 829-859;
N-root Hankel block 1000-1086).  Its consumers re-apply the drop rule per
recovered sector (L71) and integrate per sector (L73).
WHY: This ~1k-line machine exists ONLY because the per-eps-order
representation discarded exact indicial data (disease D2).  Pentagon proves
even the fixed fitter cannot always recover it (the drop rules upstream of
the fit already lost content).
DISPOSITION: subsumed-by-sector-native-representation (the entire fitter,
its salvage, and its drop rules die; plan section 2, I1/I2).  The M4 gate
ports its closed-form test cases (test_multisector_fit's non-fitter
assertions) as exact-sector tests.

### L69 Fit epistemology: falsifiability gate, strict dominance, retreat-one-order
WHAT: An N >= 2 sector candidate is accepted only if its validated run
contains a nonzero slice at log power k >= N (the k < N slices are consumed
by the weight solve, so without such a slice the extra roots were never
falsifiable; RegularizedIntegration.m:1249-1260); a candidate with more
roots must explain at least one extra validated offset PER extra root
(dominance rule 1261-1276); when validation fails mid-tower the usable
window retreats ONE EXTRA order because the by-construction rows touching
the failing region may have absorbed pollution invisible until one offset
later (1219-1234); failed weight solves invalidate the offset rather than
writing silent zeros (1156-1169).
WHY: These rules are what made the fitter stop hallucinating sectors from
truncation noise; the underlying epistemic lesson — never trust a model
component that the data could not have falsified — applies verbatim to any
DiffExp2 validation/assertion that compares reconstructions against data.
DISPOSITION: subsumed-by-exact-sector-tags (no fitting); record the
falsifiability principle in the testing guidelines for M3-M5 closed-form
gates (a gate must include data that would fail if the feature were wrong).

### L70 Root snapping: the values 0 and -branchB are DOWNSTREAM-SPECIAL
WHAT: Cluster representatives of fitted roots snap to nearby integers at the
loose cluster tolerance (1e-2) EXCEPT to 0 and to -branchB, which are
allowed only within the fit noise floor relTol — because the absolute
exponent branchB + root vanishing flips consumers onto the b = 0
resonant/limit-surviving paths (clusterRootSpecs,
RegularizedIntegration.m:855-907).
WHY: Snapping to a regular exponent trades one regular value for another
(harmless); snapping ONTO a special value silently changes which drop rule
fires — a discrete wrong-branch error.
DISPOSITION: subsumed-by-exact-sector-tags; the general rule (tolerance for
snapping onto SEMANTICALLY SPECIAL values must be the noise floor, not the
convenience tolerance) goes to Tolerances.m's documentation.

### L71 Salvage convention and trust warnings naming the first affected order
WHAT: Offsets beyond the fit's validated run keep their NON-log content on
the plain branch exponent (truncation-boundary convention) and drop residual
Logx content, counting both; the warnings distinguish "omitted data affects
the reported Laurent window starting at eps order k — results at and above k
are NOT trustworthy; rerun with more lookahead" from "only affects orders
beyond the requested window" (RegularizedIntegration.m:1560-1596,
1668-1689; limit-path analogue 2066-2137).  On stale dumps the flagged order
was exactly where pySecDec disagreement started (banana doc).
WHY: Honest completeness metadata in USER-VISIBLE form: the single most
useful diagnostic of the campaign.  This is the prototype of the plan's
EpsWindow CompleteMax semantics.
DISPOSITION: subsumed-by-EpsWindow/TWindow-metadata (plan 3.1: consumers
check need <= CompleteMax and fail loudly naming chart/sector/order; the
warning text's two-tier semantics becomes the error/ok boundary).

### L72 The zero-regulator drop rule (the pentagon killer)
WHAT: `addMonomialWithB` DROPS a contribution (with a warning) when the
sector's absolute exponent basisB is numerically zero AND the local power is
resonant (a + lp + 1 = 0): "encountered a resonant endpoint coefficient with
zero epsilon regulator; dropping contribution"
(RegularizedIntegration.m:1467-1477).
WHY: On COLLAPSED towers, a genuine x^(b eps)/x^(-b eps) pair that averages
to zero exponent hits this rule and real content is discarded — this is
where the pentagon loses its leading pole, UPSTREAM of any fitting (plan
section 2 diagnosis).  In a correct representation the case is simply the
b = 0, a + n + 1 = 0 divergence cell of the case table, which must either
cancel in the assembled combination or be a loud divergence error.
DISPOSITION: subsumed-by-Integrate.m-case-table + object-level cancellation
check (plan 3.2: b = 0 divergence is an error UNLESS cancelled in the
assembled combination, checked at the object level — never a per-term drop).

### L73 Endpoint limits: drop b != 0 per sector on the ABSOLUTE exponent; dimreg zero convention
WHAT: `EvaluateEndpointLimitSectors` resolves each integer power's tower
into sectors and applies the paper's prescription per sector: only the
absolute-b = 0 sector at absolute power x^0 contributes the limit; b != 0
sectors are set to zero EVEN at negative powers (analytic regularization /
dimreg `t^(a + b eps)|_{t=0} := 0`); a b = 0 sector at negative power, or a
confluent b = 0 sector with surviving logs, is a LOUD divergence warning
(RegularizedIntegration.m:1954-2141, esp. 2030-2062;
collapsed-exponent predecessor EvaluateLimitAtSingularity 1921-1952; at
integration bounds `FiniteEndpointConstant` keeps only the Logx-free,
x-free part, 1694-1699, used at 1828-1832).
WHY: This IS the dimensional-regularization boundary convention the whole FT
recursion rests on (limitUpper/limitLower cases); applying it to a collapsed
exponent instead of per sector discards or keeps entire towers wholesale —
the banana {1,0,0,1} limitUpper failure.
DISPOSITION: implement-in-Integrate.m + implement-in-API.m (EndpointLimit:
constant of the (0,0,0) sector; b != 0 dropped EXACTLY; divergence loud —
plan 3.3), trivially exact under sector-native tags.

### L74 Meromorphic (b = 0) regular powers take Hadamard finite parts — never the unit regulator
WHAT: For b = 0, integer a <= -1, Logx-free integer-step series at an
endpoint, the singular powers are integrated with the Hadamard finite-part
antiderivative directly (formal divergence at x = 0 dropped, pairing with
the partner arm / neighboring segment giving PV), and the REGULAR remainder
is peeled off and integrated exactly (RegularizedIntegration.m:564-652).
WHY: Routing regular powers through the unit-regulated subtraction machinery
leaks boundary terms of higher-eps residues into LOWER eps orders (e.g. an
eps^-1 coefficient c1 x leaks c1 (B log B - A log|A| - (B - A)) into eps^0)
— the second 9aeb300 bug ("unit-regulator over-regularization",
Docs/FeynmanTrickBoxFamilyStatus.md:53-55).  NO epsilon dependence may be
introduced for meromorphic content.
DISPOSITION: implement-in-Integrate.m (the b = 0 column of the case table is
exact finite-part/PV with no eps mixing, by construction).

### L75 The formal unit regulator only applies when no Logx towers survive
WHAT: When the EXTRACTED exponent vanishes but a divergent power is present,
a formal regulator b = 1 is substituted ONLY if the coefficients carry no
residual Logx towers; towers always take the subtraction path because they
carry the TRUE sector exponents — shifting the basis exponent by one would
displace every recovered sector (RegularizedIntegration.m:560-661, comments
560-567).
WHY: Banana doc fix 3: a vanishing averaged exponent with perfectly
cancelling weights otherwise had every recovered sector shifted by one.
DISPOSITION: subsumed-by-exact-sector-tags (no extracted-average exponent
exists; each sector integrates against its own exact b).

### L76 Endpoint integration = subtract divergent powers + closed-form monomial basis with pole depth p+1
WHAT: The subtraction integrator removes from g every local power with
a + lp + 1 <= 0, integrates the remainder with the ordinary antiderivative,
and adds back each subtracted monomial's exact closed form: for nonzero
total exponent s = p' + b eps, the basis is
`Exp[s logValue] Sum_j Binomial[logPower, j] logValue^(logPower-j) (-1)^j
j! / s^(j+1)` expanded in eps (denominators (a + n + 1 + b eps)^(j+1),
j = 0..p — pole DEPTH p+1); for zero exponent the leading term is
`(-1)^p p! b^(-p-1)` at eps^(-p-1) (monomialBasis,
RegularizedIntegration.m:1395-1465; assembly 1467-1657).
WHY: This is the exact integral counterpart of the sector representation and
the source of truth for Integrate.m's case table: one closed form per
(sign of a + n + 1, b zero/nonzero, p) cell, with the eps-Laurent window
shifting by the pole depth.
DISPOSITION: implement-in-Integrate.m (the full case table, plan 3.2; the
old closed forms are the M3/M4 oracle).

### L77 1/(b eps) is a Laurent SHIFT, not an inversion (ApplyRegularizationStep)
WHAT: The IBP-style regularization step
`Integrate[x^(a+b eps) g] = Integrate[x^(a+1+b eps)/(1+a+b eps) ((2+a+b eps)/c g - (1-x/c) g')]`
handles a = -1 by shifting epsMinPower down by 1 and scaling by 1/b —
explicitly a window shift — while a != -1 expands the prefactor as a finite
eps convolution (RegularizedIntegration.m:347-466; pole-shift counting for
window pre-extension 791-827).
WHY: Two representational lessons: (1) eps-pole creation must be tracked as
window metadata (epsMinPower), never by materializing 1/eps symbols;
(2) the step needs the gList pre-extended by the number of expected pole
shifts or top orders are silently lost (extraRegOrders, 808-812).
DISPOSITION: implement-in-EpsSeries.m (Laurent division/shift semantics:
denominator with vanishing eps^0 part shifts kmin AND CompleteMax — plan
3.2); the IBP stepping itself is subsumed-by-Integrate.m closed forms.

### L78 Top-eps-order Logx content is a truncation boundary, not data
WHAT: Explicit Logx terms in the HIGHEST available eps coefficient of an
endpoint subtraction are dropped and counted
(`droppedTopOrderLogTerms`, RegularizedIntegration.m:1643-1651, warning
1659-1667): "treated as truncation-boundary terms and not allowed to
generate lower Laurent poles".
WHY: A Logx at the top order would, through the monomial basis, deposit
content at LOWER eps orders — but the top order is incomplete by
construction (its partner terms lie beyond the window), so trusting it
manufactures spurious poles.
DISPOSITION: subsumed-by-EpsWindow-honesty (CompleteMax excludes the
incomplete top by construction; no content beyond CompleteMax may influence
lower orders — assert this in Integrate.m).

### L79 Evaluation blow-up rescue via Pade
WHAT: `EvaluateIntegralAtPoint` detects |value| > 10^100 at a regular
evaluation and retries through a Pade approximant of the integrand series,
keeping the smaller result (RegularizedIntegration.m:1845-1870).
WHY: Truncated antiderivatives evaluated near the chart edge can blow up
polynomially while the true function is finite; Pade is the cheap
resummation.  The acceptance rule (keep the SMALLER) is a heuristic, flagged
as such.
DISPOSITION: waived-because-symptom-of-bad-chart-geometry — DiffExp2's
chart-sizing (L27) plus the re-expansion truncation contract (plan 3.2,
SectorSeries) make evaluation points interior by construction; if the gate
ever fires in the new code it is a LOUD error, not a rescue.

### L80 Prefactor integration mechanics: endpoint snapping, side phases, power absorption, affine Jacobians
WHAT: `IntegrateSegmentWithPrefactor[Laurent]` (RegularizedIntegration.m:
2628-3010): (a) segment maps are SNAPPED so chart endpoints hit the exact
integration bounds 0/1 and IBP-factor roots
(snapMainExpression/snapValuesFromFactors,
FeynmanTrick/DiffExpIntegration.m:151-202; snapLocalCoordinate,
RegularizedIntegration.m:2388-2398) — tolerance currently the 10^-40
fallback (L14); (b) boundary prefactor powers (x - lower)^alpha,
(upper - x)^beta become singular-power shifts of the local exponent a with
`|jacobian|^alpha` scaling and an explicit branch phase `Exp[-I Pi power]`
when the local coordinate runs negative (localSidePhase,
RegularizedIntegration.m:1798-1806, used 2698-2723); (c) smooth prefactor
parts are series-multiplied onto g; (d) any negative powers the rational
factor introduces are absorbed into the exponent a so g starts at x^0
(2784-2822, 2952-2989); (e) Jacobians are evaluated AT THE SEGMENT MIDPOINT
— valid only for affine maps (2663-2664, 2869-2870; cf. L29).
WHY: Each sub-mechanism encodes a discovered failure: unsnapped endpoints
put the singular point epsilon-off the chart center (wrong drop rules);
missing side phases flip the imaginary part of (1-x)^(beta) integrals;
unabsorbed rational poles feed g-series with negative powers that the
integrator's case analysis does not expect.
DISPOSITION: implement-in-Integrate.m (prefactor spec handling) and
implement-in-SectorSeries.m (rational multiply with pole absorption =
partial-fraction split across charts, plan 3.2); affine-only enforced per
L29.

---

## I. FT pipeline: eps windows, trims, prefactors

### L81 Combine the IBP-weighted integrand BEFORE integrating
WHAT: `IntegrateCombinedMasters` builds `Sum_j c_j(x) f_j(x)` per segment
per eps order FIRST and integrates the combined series; comment marked
CRITICAL: individual c_j have poles at x = 0, 1 that CANCEL in the sum —
term-by-term integration would diverge
(FeynmanTrick/DiffExpIntegration.m:533-560, 681-794).  The same
combine-before-limit requirement applies to EvaluateLimitFromTransport
(survey/A1 list: it currently weights AFTER the per-master limit,
DiffExpIntegration.m:1000-1109 — a latent ordering hazard when coefficients
diverge, see L88).
WHY: Pole cancellation between IBP coefficients is structural, not
numerical; any path that evaluates masters separately against divergent
weights destroys it.
DISPOSITION: implement-in-FT-layer-hardening (combine-before-limit
everywhere) + implement-in-Integrate.m (rational-multiply + exact
cancellation at object level, plan 3.3).

### L82 Trim incomplete combined top orders; DELETE the "keep at least one order" floor
WHAT: A combined output order n is complete only if EVERY active master
still has transport data at shifted index n + k_j - k for every nonzero IBP
coefficient order k; the combine trims the incomplete top orders and warns
when the trim cuts into the requested window
(FeynmanTrick/DiffExpIntegration.m:633-667).  BUT line 658-661 keeps "at
least one output order so callers receive a well-formed result" even when
NOTHING is complete — masking an empty window.
WHY: Partially combined top orders are incomplete-but-present data; the old
fitter consumed them as truth (banana root cause 1).  The floor converts
"no trustworthy data" into "one untrustworthy order" — D1.
DISPOSITION: implement-in-FT-layer-hardening (trim kept, floor DELETED per
plan A1) — long-term subsumed-by-EpsWindow-arithmetic (window intersection
across masters is the same computation with honest metadata).

### L83 Cumulative eps-prefactor threading and ShiftRawBoundariesToFinite
WHAT: Masters are transported in the prefactored basis J_j = eps^(k_j) I_j;
the combine convolves with shift k_j
(FeynmanTrick/DiffExpIntegration.m:669-679), the per-level requirement adds
prefactors CUMULATIVELY (`epsOrder + prefacs[j] - coeff MinPower +
allowance`, 1455-1465), and after each level
`ShiftRawBoundariesToFinite` applies ONE UNIFORM shift = max pole depth
across all masters and ZERO-PADS every master up to the level max
(1258-1280).  Prefactor metadata is preferred over re-inference; the
`CheckEpsPoles` inference guard tests the d-form and can NEVER fire
(1495-1510; FeynmanTrickIteration.m:455, box doc latent list).
`ExpandIBPCoeffLaurent` silently Floors a non-integer MinPower (1224).
WHY: The uniform shift is what makes "boundary values" a plain non-negative
list for DiffExp — but the zero-padding manufactures fake known-zero orders
(D1), and the dead inference guard means wrong bases would not be caught.
The cumulative convention is the one the plan 3.4 budget formula must match.
DISPOSITION: implement-in-FT-layer-hardening (zero-padding replaced by
CompleteMaxPower metadata per A1; dead guard removed; Floor made an error);
the basis bookkeeping is subsumed-by-exact-Laurent-bookkeeping in eps (plan
goal (b)).

### L84 IBP coefficient poles extend the SEGMENTATION alphabet
WHAT: Before each level's transport, the denominators of the level's IBP
reduction coefficients are factored and injected into
`DiffExp`State`MatricesIrreducibleFactors` (appendMatrixFactors,
FeynmanTrick/DiffExpIntegration.m:204-222; CollectLevelIBPSingularFactors
1326-1391; threading 1788-1807) and their real roots in [0,1] become snap
targets (180-202, 309-315).
WHY: The integrand's singularities are NOT only the DE matrix's: IBP
coefficients multiply the solution, and their poles (box L1's t = 1/4,
box_triangle L3's 0.9617) must cap chart radii and trigger interior-pole
handling.  Without this, segments straddle invisible poles.
DISPOSITION: implement-in-Transport.m (segmentation accepts an
ExtraSingularFactors input as first-class; plan 3.2 names
CollectLevelIBPSingularFactors in the multiplication-closure contract).

### L85 Variable-context pinning (FIRE Global` vs package contexts)
WHAT: IBP coefficients (from FIRE, in Global`) and DiffExp segments
(variables auto-detected from matrix filenames into whatever context is
current) can hold SAME-NAMED symbols in different contexts; the code remaps
by `SymbolName` matching before any substitution — in the combine
(FeynmanTrick/DiffExpIntegration.m:617-631), in the limits (977-991), and in
factor injection (374-380).  The FT layer reaches into
`DiffExp`State`FEC[Variables]` for this (the x27 State reach-in class).
WHY: A context mismatch makes `coeff /. var -> expr` a silent no-op — the
coefficient stays symbolic and downstream numeric checks treat it as zero or
garbage.  This bit repeatedly (memory note "Wolfram package context traps").
DISPOSITION: implement-in-Config.m (variable-context pinning at LoadSystem:
ONE canonical symbol per variable name, all inputs normalized to it on
entry, mismatches loud).

### L86 Limit segments are selected by FEYNMAN-PARAMETER value, not path bounds
WHAT: `EvaluateLimitFromTransport` picks the segment whose mapped Feynman
parameter actually reaches the boundary (minimum distance of mapped local
bounds to the target), explicitly NOT segment element 3 (path-parameter
bounds), "using it can confuse the upper endpoint with the reversed lower
transport" (FeynmanTrick/DiffExpIntegration.m:930-958); the lower transport's
segments are stored REVERSED when combined (499-501).
WHY: The two-direction transport (fixed point -> 0 and -> 1) makes path
order and parameter order disagree on one side; selecting by path bounds
evaluated the wrong segment's series at the wrong endpoint.
DISPOSITION: implement-in-API.m (EndpointLimit takes the chart at the
endpoint by CHART CENTER identity, not list position — bidirectional
transport is first-class in the plan's API).

### L87 The Quiet[Check[..., 0]] silent-zero chains (A1 sites)
WHAT: Endpoint evaluation falls back to 0 on ANY error:
series/symbolic endpoint values (FeynmanTrick/DiffExpIntegration.m:1037-1059)
and IBP coefficients at boundaries via `Limit -> substitution -> 0`
(1076-1098); the limit path also silently zeroes endpoint-divergent IBP
coefficients, and NEGATIVE-INDEX boundary requests fall through to "direct"
which can silently delete numerator powers
(Docs/FeynmanTrickBoxFamilyStatus.md:143-148).
WHY: Every one of these converts an exceptional condition into a silent 0 in
a BOUNDARY VALUE — the worst possible D1 placement, since boundary values
seed all lower levels.
DISPOSITION: implement-in-FT-layer-hardening (Stage A1: each site becomes
either a handled case with a closed-form rule or a loud error; the
divergent-coefficient case must combine-before-limit per L81).

### L88 Hardcoded expOrd = 30 and phantom-order laundering
WHAT: The symbolic-transport-value branch of the combine series-expands at a
hardcoded order 30 (FeynmanTrick/DiffExpIntegration.m:753-764); afterwards
the combined SeriesData is Normal-ized and RE-Series-ed at
`Ceiling[nmax/den]` (774-784), which can mint one phantom top order beyond
the data (the same laundering exists in `ApplyAnalyticContinuation`,
SeriesOps.m:81-84: re-Series of Normal at Floor[nmax/den]).
WHY: 30 is unrelated to ExpansionOrder (silently truncates if the transport
used more); the re-Series laundering re-labels a truncated window as if it
extended one order further — both are honest-window violations.
DISPOSITION: implement-in-FT-layer-hardening (order from the transport
metadata; window-preserving normalization) / subsumed-by-EpsSeries+TWindow
arithmetic for the new core.

### L89 Config reset trap: prescriptions must be re-applied after LoadConfiguration
WHAT: `LoadConfiguration` resets to defaults (DiffExp/DiffExp.m:95-99), so
TransportLevel re-applies DeltaPrescriptions and
AbortOnAnalyticContinuationFail after EVERY reload — including the reload
between the lower and upper transports (FeynmanTrick/DiffExpIntegration.m:
409-414, comment "CRITICAL" at 465-469); injected ExtraSingularFactors are
likewise re-appended (461-463).
WHY: A missed re-application silently transports without prescriptions —
producing exactly the pentagon-class "not recognized as branch point"
warnings on configs that LOOK correct.
DISPOSITION: implement-in-Config.m (immutable config snapshot passed to each
transport call; no global mutable config to desynchronize).

### L90 FT operating-point conventions
WHAT: Fixed Feynman-parameter anchor 11/23 (rational, away from 0, 1, 1/2
and observed singular-factor roots; FeynmanTrick/DiffExpIntegration.m:270);
both transports save expansions (SaveExpansions = True, 440-446, 472-480);
prescriptions {x, +1} and {1-x, +1} always included
(deltaPrescriptionsForFactors, 224-238); UseMobius False (L29), UsePade
False (approximants would have to be re-derived per combined integrand),
HomogeneousSolve DontExpand (L53), DivisionOrder 4 (L27), Predivision.
WHY: This tuple is the configuration under which every campaign validation
holds; M5's cutover must reproduce results under an equivalent declared
configuration, not rediscover it.
DISPOSITION: implement-in-Config.m (documented FT profile) +
implement-in-Tests (the per-example pinned configs, see L96).

---

## J. Loading, configuration, API conventions

### L91 Per-eps-order slice exports are TRUNCATED; exact work needs the full export
WHAT: The classic matrix dir format is `d<var>_<ord>.m` slices loaded up to
EpsilonOrder (MatrixLoading.m:29-55, 91-111); FT's exporter writes slices by
`SeriesCoefficient` at a fixed order (MatrixExport.m:42-83) — a
rational-in-eps entry is silently Taylor-truncated.  The exact format
already exists: `d<var>_d.m` closed form (MatrixLoading.m:42-46, 114-124)
and `ExportGeneralMatrix`/`d<var>_full.m` (MatrixExport.m:91-114), but
nothing consumed it.
WHY: The I1 indicial contract (exact eigenvalue factorization) is
UNCERTIFIABLE from truncated slices — banana ships slices 0..4, so A_{-1}(eps)
looks quartic; FT systems are rational in eps.  "Nothing ever needs fitting"
degrades to "fitting replaced by truncation error" without the full format
(legacy findings 1-2).
DISPOSITION: implement-in-API.m (LoadSystem consumes the exact full export;
slice dirs accepted ONLY for legacy parity transport; M5 switches
MatrixExport to the full format — plan I1).

### L92 Missing slice files are silently zero-filled
WHAT: A missing `d<var>_<ord>.m` becomes a zero matrix with only a
Verbosity-1 note "Assuming M[var][ord] is zero" (MatrixLoading.m:104-109,
141-144).
WHY: Legitimate (sparse eps-dependence is common) but indistinguishable from
a typo'd filename or a failed export — a D1 trap that has burned real runs.
DISPOSITION: implement-in-API.m (LoadSystem: zero slices must be DECLARED
(manifest or explicit zero file) or the load errors listing assumed-zero
entries; full-format loading makes this moot for FT).

### L93 Input-matrix algebra constraints checked at load
WHAT: LoadMatrices rejects unsupported function heads
(MatrixLoading.m:148-157), higher-than-square roots (159-164), and
reducible square-root arguments (176-179); DeltaPrescription polynomials
must be irreducible (DiffExp/DiffExp.m:164-168); the closed-form path also
admits eps in the variable check (MatrixLoading.m:118).
WHY: Every constraint marks the validated algebra boundary; violations
beyond it failed in undefined ways downstream.
DISPOSITION: implement-in-API.m (LoadSystem validation, incl. the v1 loud
error for irrational x-dependence per plan non-goals).

### L94 The singularity alphabet is built at eps -> 0
WHAT: `MatricesIrreducibleFactors` collects irreducible denominator factors
and sqrt arguments of all matrices WITH `eps -> 0` applied, deduplicated up
to sign (MatrixLoading.m:214-231); segmentation derives all singularities
from it (L26) plus FT's injected factors (L84).
WHY: eps-dependent denominator factors (possible in raw FT exports) would
otherwise generate spurious eps-dependent "singularities"; but the eps -> 0
projection also means a factor vanishing ONLY at eps = 0 is invisible —
fine for Fuchsian-in-x systems, an assumption worth asserting.
DISPOSITION: implement-in-Transport.m (segmentation from the exact matrix's
x-denominators at eps = 0, with an assertion that clearing eps from
denominators is exact — ties into I1's eps-rational input).

### L95 The kept configuration surface and the validated accessor
WHAT: Live keys verified in code/examples: AccuracyGoal +
"AccuracyGoalValidate" (L19/L20), ChopPrecision (+ sync, L3),
DeltaPrescriptions (parse `poly ± I delta` -> {poly, sign},
DiffExp/DiffExp.m:140-148, user list tracked separately 146-148),
DivisionOrder, EpsilonOrder, "EstimateError" (True/"Fast"/False parse,
DiffExp/DiffExp.m:150-158), ExpansionOrder, LineParameter (collision check
vs Variables, DiffExp/DiffExp.m:131-135), LogFile (appends to $Output AND
$Messages, 107-114), MatrixDirectory, RadiusOfConvergence,
RationalizationTolerance, SegmentationStrategy, UseMobius, UsePade,
UseRationalRecurrence, Variables, Verbosity (+VerbosityDebug),
WorkingPrecision, CrosscheckLevel/CrosscheckFlags, SaveExpansionsCompress
(+Directory/+Order, Transport.m:852-866), "AbortOnAnalyticContinuationFail",
"IgnoreIndicialCheck", "KeepMatrixExpansions", "HomogeneousSolve",
"InvWronskSolver", "LinearSolveChopPrecision", IntegrationStrategy,
"Parallel" (dead).  `CurrentConfiguration[]` exposes a curated subset
(DiffExp/DiffExp.m:86-92).  FEC is a plain association: a misspelled key
read returns Missing and flows on silently (the FEC-class bug).
WHY: Each key is load-bearing for some pin (Docs/specs/Config.md per M0
carries the kept/waived table); the silent-miss accessor is the original D1
config disease.
DISPOSITION: implement-in-Config.m (schema-validated accessor — unknown key
read/write is an error by construction; every kept key mapped or waived in
Docs/specs/Config.md).

### L96 Parity oracles must pin the EXACT old config per example
WHAT: Checkpoint locations and segment structure depend on
UseMobius/UsePade/DivisionOrder/RoC (e.g. banana classic =
UnequalMassConfiguration with UsePade, UseMobius, RoC -> 10,
Banana_example.m:73-78; FT profile per L90); legacy review finding 15.
WHY: Comparing new-code checkpoints against an oracle run under a different
geometry compares apples to oranges and wastes kernel time.
DISPOSITION: implement-in-Tests (Tests/PINS.md records the full old-config
tuple per parity example; M4 gate text references it).

### L97 Both `eps` and `\[Epsilon]` are accepted
WHAT: `eps := \[Epsilon]` as a delayed alias (Symbols.m:22); matrix loading
admits both in closed-form files (MatrixLoading.m:118).
WHY: User inputs and exported matrices mix the two freely; rejecting either
breaks every existing example.
DISPOSITION: implement-in-API.m (both accepted, normalized at entry — plan
3.2 API list).

### L98 PrepareBoundaryConditions conventions
WHAT: (a) per-integral "?" wildcard = ignore this integral's BCs
(Transport.m:63-65, consumed via IgnorePositions at 291-294, used in
Banana_example.m); (b) closed-form eps expressions are auto-expanded with a
sanity check that the eps series starts at order >= 0 (Transport.m:67-74);
(c) for LINE (asymptotic) BCs, every order below the leading eps coefficient
is set to an exact zero SERIES with a printed assumption (Transport.m:76-84)
— a documented convention, not a guess; (d) leading-asymptotics extraction
takes the leading series coefficient of each entry along the line
(Transport.m:91-106 via LeadingCoefficientSeries, SeriesOps.m:174-189
including its fractional-power denominator merge); (e)
`Log[a x] -> Log[a] + Logx`, `Log[x] -> Logx` rewriting at the end
(Transport.m:120).
WHY: These define the exact input contract of every example; (c) in
particular is load-bearing for asymptotic BCs (the eps-leading behavior
defines all lower orders to be exactly zero in the limit).
DISPOSITION: implement-in-API.m (PrepareBoundaryConditions compatibility per
plan API list, semantics identical).

### L99 Asymptotic-BC reparametrization: orientation check and log rescaling
WHAT: When BCs live at an asymptotic limit of a different line
parametrization, `CheckBoundaryConditionsAndReparametrize` requires the
inverse relation's derivative at 0 to be POSITIVE ("Asymptotic boundary
conditions should be oriented in the same direction as the integration
line", LineSegmentation.m:196-201), substitutes a SERIES of the relation at
2x ExpansionOrder (ISeriesChangeCoefficient = 2, State.m:210;
LineSegmentation.m:204-207), and rewrites the emergent
`Log[a yy] -> Log[±a] + Logx` with sign handling, erroring if any Log
survives (LineSegmentation.m:209-222).
WHY: A reversed orientation silently flips which branch the x^(b eps)
asymptotics live on; surviving Log heads mean the rescaling failed and
constants would be wrong.
DISPOSITION: implement-in-API.m (same checks; in sector representation the
rescale acts exactly on tags: x -> c x shifts coefficients by c^(a + b eps),
log shift Logx -> Logx + Log[c]).

### L100 Singular-endpoint mode returns the SERIES
WHAT: TransportTo detects a target endpoint that is a singularity (chop
comparison at ChopPrecision, Transport.m:605-613) and, on the final segment,
integrates but does NOT evaluate: it stores the analytic-continued local
series as the result with `"EndpointIsSingularity" -> True`
(Transport.m:1053-1106, output 1234-1246).
WHY: The FT pipeline's entire endpoint-limit/integration machinery consumes
exactly this mode (TransportLevel transports TO 0 and TO 1); evaluating at
the singularity would be Indeterminate.
DISPOSITION: implement-in-API.m (first-class singular-endpoint mode
returning the LocalSolution at the endpoint chart — plan API list; in
DiffExp2 this is just "return the chart object").

### L101 SaveExpansions / ToPiecewise / SampleAtList output surface
WHAT: SaveExpansions stores per segment: {line, lineRelation, main-bounds,
local-bounds, continued series}, optionally Compress-ed to disk keyed by
content hash (Transport.m:845-868, 1061-1081; SaveExpansionsCompress
Directory/Order config); ToPiecewise rebuilds piecewise functions with
optional Pade (Transport.m:1254-1313); SampleAtList exports values at
requested points during transport (Transport.m:870-897).  This 5-tuple is
the FT pipeline's transport interchange format (RegularizedIntegration and
DiffExpIntegration consume it).
WHY: The interchange format IS the FT<->DiffExp contract; the rewrite must
either keep it as a shim or replace it with the LocalSolution chain plus an
adapter (R4 shim audit).
DISPOSITION: implement-in-API.m (SaveExpansions/ToPiecewise equivalents per
plan; the 5-tuple documented in the shim contract deliverable).

### L102 Pade evaluation: per-Logx-coefficient, half-window diagonal, currently a SILENT fallback
WHAT: `GetPade` applies PadeApproximant per Logx coefficient at diagonal
order `Floor[(window+1)/2]` under `$MinPrecision = WP`, chopping inputs at
CP (Pade.m:30-53); on PadeApproximant failure it WARNS and silently
evaluates the plain truncated sum instead (Pade.m:43-48).  It is
load-bearing in the marching loop: segment-boundary evaluation and the
error probe both flow through SEval/SEval1/SEval2 (Transport.m:905-924,
317-323); classic examples require UsePade -> True (Banana, Pentagon, 2F1
resonant tests).
WHY: Matching points sit at radius/k where plain partial sums converge
slowly — Pade buys accuracy-per-order, not luxury; the silent fallback means
a failed approximant quietly degrades the segment's accuracy with only a
print (the error probe partially catches it, but the probe uses the SAME
evaluator).  The R6 benchmark is against old-with-Pade.
DISPOSITION: implement-in-SectorSeries.m (Pade as the evaluation accelerator
per plan 3.2 decision, ~+100 lines, with the fallback made LOUD and recorded
in the segment's error estimate).

### L103 d_1.m canonical-form path
WHAT: A canonical `d_1.m` (dlog-form matrix, entries `Sum c_i Log[l_i]`)
is detected, validated, and its alphabet logs differentiated into the eps^1
matrix at expansion time (MatrixLoading.m:59-85, 343-353).
WHY: Supported the original canonical-form workflow
(MultiplePolylogarithms example); no FT example and no M4 parity line
(bubble/sunrise/2F1/banana) uses it.
DISPOSITION: waived-because-no-pin-suite-example-consumes-it (loud
not-supported error at LoadSystem if d_1.m present; revisit if a classic
parity target is added — record in ExportDisposition).

### L104 Memory policy: matrix expansions are dropped per segment
WHAT: After each segment the previous chart's cached expansions are dropped
unless `"KeepMatrixExpansions" -> True` (ClearMatrices,
MatrixLoading.m:235-249; called Transport.m:1215-1218, 216, 503).
WHY: Per-chart expanded matrices at order 50+/WP 500 are large; hundreds of
segments would exhaust memory.  The config escape exists for repeated
same-line experiments.
DISPOSITION: implement-in-Transport.m (chart artifacts dropped when the
marching frontier passes, kept behind a debug flag).

### L105 Diagnostic environment hooks are part of the operational surface
WHAT: Env-gated, deliberately retained: DEBUG_DUMP_DISPATCH_DIR (per-call
ctx/bVec dumps — the M3 replay harness, Dispatch.m:10-16);
DIFFEXP_DUMP_LAURENT_DIR + DIFFEXP_DUMP_LAURENT_ABORT_AFTER (definite-
integral input dumps + replay via Scripts/eval_dump_generic.m,
RegularizedIntegration.m:225-256, 2585-2600); DEBUG_SING_PART,
DEBUG_POWER_COEFF, DEBUG_SECTOR_FIT (Recurrence.m:643-648, 673-685;
ResonantRecurrence.m:1130-1135; RegularizedIntegration.m:1013-1097);
$LogStrategyDispatch trace (Dispatch.m:61-85); $DebugFuchsianizedCheck /
$DebugBlockResidualSeries self-checks (ResonantRecurrence.m:1387-1457);
FT_INTEGRATION_POLE_ALLOWANCE, FT_SAVE_TRANSPORT_DIR,
FT_STOP_AFTER_BOUNDARY_LEVEL (runner).
WHY: The dump-replay loop (capture once, iterate in seconds without the
kernel-expensive pipeline) is what made the campaign debuggable under a
single license; M0's oracle pre-generation uses these hooks directly.
DISPOSITION: implement-in-API.m / implement-in-Transport.m (equivalent
dump/replay hooks in DiffExp2 from day one; the old hooks stay in the frozen
oracle for M3/M4 dump generation).

### L106 LastErrorContext: every loud error ships its context
WHAT: Nearly every ReportError site first stores the relevant local state in
`DiffExp`State`LastErrorContext` (declared State.m:79; e.g.
Transport.m:376-381, 903; Frobenius.m:64; ResonantRecurrence.m:75-83,
413, 567, 1360; LineSegmentation.m:26, 41, 52, 157, 180; Wronskian.m:40),
and `LastBoundaryMatchingContext` captures underdetermined matches
(Transport.m:394-405, 419-426).
WHY: With a single kernel and multi-hour runs, a post-mortem inspectable
context is the difference between one debug cycle and five.
DISPOSITION: implement-in-API.m (structured error objects carrying context —
the plan's "fail loudly naming (chart, sector, order)" plus payload).

---

## Cross-reference: RewritePlan section 5 seeds -> ledger items

| Seed (plan section 5) | Ledger |
|---|---|
| input precision 2x WP + $MinPrecision/$MaxExtraPrecision floors | L1, L2 |
| ChopPrecision/LinearSolveChopPrecision sync | L3 |
| expansion-order adaptive search (AccuracyGoalValidate Before/After) | L19, L20 |
| DivisionOrder both-charts geometry (GetCPL/GetCPR), FT k=4 | L27 |
| predivision segment-count digit budgeting | L18 |
| coupling-depth t-order degradation (MaxCouplingOrder, ISafety, ICurrEvalErrorSeriesDecrease) | L15 |
| two-point error probe avoiding x=0, per-indeterminate, abort>1 | L16 |
| RoC chart rescaling O(1) coefficients; banana REQUIRES 10 | L28 |
| Pade evaluation rationale + silent-fallback fix | L102 |
| Mobius charts (who needs them; banana classic does) | L29 |
| complex-singularity ghost projection (Re, Re±Im, suppression) | L26 |
| delta-prescription derivation (parity, leading sign, conflicts, sqrt auto, SingularityCheck) | L34, L36, L38 |
| interior-split real-log convention (9aeb300) | L37 |
| principal-branch crossing convention (+1 no-op; -1 Logx shift; denom>2) | L35 |
| denominator-clearing recursion fast mode | L41 |
| segment block caching | L42 |
| eps-window trimming semantics incl. uniform per-level shift | L82, L83, L14 |
| FIRE variable-context pinning | L85 |
| numerical-zero leading-coefficient skipping (generalizes to matching) | L9 |
| resonant-source empty-particular hole | L65 |

Legacy-review findings 1-18 map to: 1->L91 (+L10, L44 context), 2->L91,
3->L34, 4->L26/L27(radius), 5->L23, 6->L15, 7->L16, 8->L18, 9->L98-L101
(+L86, L97), 10->L95, 11->L102+L28, 12->L27, 13->L14+L82+L83+L87, 14->L35,
15->L96, 16->L29, 17->(plan-text fix; LocalSeries disposition in L59),
18->L1.
