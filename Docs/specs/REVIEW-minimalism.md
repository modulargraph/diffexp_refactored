# Adversarial spec review — minimalism and error contracts (M0 task 13)

Scope: Docs/specs/{Tolerances,Config,EpsSeries,SectorSeries,Indicial,Solve,
Transport,Integrate,API}.md, cross-checked against Docs/RewritePlan.md (v2)
and the frozen old tree.  FTShimContract.md, LessonsLedger.md, and
ExportDisposition.md did not exist at review time; defects that depend on
them say so.  Old-code claims load-bearing for a defect were re-verified
against the working tree (citations below).  Defects only; no praise.
Severity: blocker = implementation from the spec as written produces wrong
results, an unimplementable interface, or fails its own milestone gate;
major = contradiction/gap that will surface as a silent fallback, dead
surface, or cross-module breakage; minor = consistency/coverage debt.

---

## BLOCKERS

### 1. [blocker] Tolerances.md §2 (`LaurentLeadTol`) vs EpsSeries.md §6 N-1/N-2 + §8 test 11 vs Integrate.md §2.2.3 — the central zero-classification tolerance is mis-derived; the specs' own tests fail under it

Problem.  Tolerances.md defines `LaurentLeadTol[chopDigits] ->
10^(-Floor[chopDigits/2])` (= 10^-125 at the default chopDigits 250, 10^-75
at WP 300) and bans absolute floors (F9, "there is deliberately NO ...").
But the only campaign-validated calibration of this exact test is
`tolRel = Max[10^(-ChopPrecision/2), 10^(-24)]` — verified at
DiffExp/IntegrationStrategies/Recurrence.m:618 with the in-code comment
"Cancellation residues at apparent singularities sit far above
10^(-ChopPrecision/2) (e.g. ~1e-29 relative at WP 300 ...)".  EpsSeries.md
itself pins this requirement (N-2: "laurentLeadTol must be loose enough to
absorb ~30 digits of cancellation at WP 300") and then writes a unit test
that is arithmetically impossible under Tolerances' derivation: EpsSeries
test 11(a) requires `10^-60` at scale 2 to classify ZERO, but
`10^-60 > LaurentLeadTol[250]*2 = 2*10^-125` by 65 decades — the test fails,
and worse, in production every genuine ~1e-29-relative cancellation residue
classifies NONZERO, so Integrate.md's cancellation gate (§2.2.3 step 3) fires
false E2 divergence errors on correct data (the gate cites the campaign class
`Max[10^(−ChopPrecision/2), 10^−24]` while naming `chopFloor`, which
Tolerances defines as the ABSOLUTE 10^-chopDigits noise floor — a second,
independent mis-binding).  Solve.md §3.6 step 3 has the same mis-binding
("below chopFloor (RELATIVE to the column's own coefficient scale)").
Additionally, with any threshold near 10^-24, Tolerances'
`$AmbiguityBandDecades = 10` puts the documented residue population
(~1e-29 relative) INSIDE the band [tol/10^10, tol*10^10] = [10^-34, 10^-14],
so NumericallyZeroQ aborts on correct data.

Amendment (exact text):

In Tolerances.md §2, REPLACE the `LaurentLeadTol` derivation line with:

> - `LaurentLeadTol[chopDigits_Integer] -> Max[10^(-Floor[chopDigits/2]), 10^-24]`
>   RELATIVE leading-coefficient zero test for eps-Laurent window trimming
>   and for object-level cancellation gates.  The 10^-24 floor is the
>   campaign-validated calibration (DiffExp/IntegrationStrategies/
>   Recurrence.m:613-618; Docs/FeynmanTrickBoxFamilyStatus.md:116-119 "the
>   1e-24 floor is load-bearing"): imperfect cancellations at apparent
>   singularities leave residues up to ~1e-29 RELATIVE at WP 300 that MUST
>   classify zero, while genuine leading content is O(1) relative or exactly
>   zero.  This floor is a RELATIVE quantity and does not violate F9 (which
>   bans absolute tests on raw coefficients).  Calls through NumericallyZeroQ
>   with this tolerance use a band of `$LaurentLeadBandDecades = 4` (NOT the
>   default 10): the populations to separate are <= 1e-29 (zero side) and
>   >= ~1e-20 (signal side) — band edges 10^-28 / 10^-20.  M1 unit must
>   verify both campaign populations land outside the band.

In Tolerances.md §2, ADD next to `$AmbiguityBandDecades`:
`$LaurentLeadBandDecades = 4` (band override for LaurentLeadTol-based calls;
NumericallyZeroQ takes the band as an optional fifth argument defaulting to
`$AmbiguityBandDecades`).  Recompute Tolerances T4 (`LaurentLeadTol[250] ===
10^-24`), T9 (band edges become tol/10^4 and tol*10^4 for laurentLeadTol
calls), and EpsSeries test 11 expectations (which then pass: 10^-60 vs scale
2 is firmly below 10^-24*2/10^4).  In Integrate.md §2.2.3 step 3 and §7,
REPLACE `chopFloor` with `LaurentLeadTol` (chopFloor remains the absolute
coefficient noise floor and is not a classification tolerance).  In Solve.md
§3.6 step 3 and §6 L7, REPLACE `chopFloor (RELATIVE ...)` with
`LaurentLeadTol (relative, Tolerances.md)`.

### 2. [blocker] Tolerances.md §2/§3 (`NumericallyZeroQ`) vs EpsSeries.md §2 (`ESCoeffZeroQ`) — two incompatible zero predicates, and the window-moving one is the wrong one

Problem.  Tolerances.md mandates: "EpsWindow ... movement is decided in
EpsSeries.m, but ONLY through `NumericallyZeroQ` with `LaurentLeadTol` — no
other predicate may move a window boundary" (§3), where NumericallyZeroQ is
TERNARY (True / False / loud error in the ambiguity band, context string
mandatory).  EpsSeries.md instead defines its own BINARY predicate
`ESCoeffZeroQ[c, scale]` ("True iff ... `Abs[N[c, wp]] < laurentLeadTol *
scale`", no band, no error, no context) and uses it for ESTrim/ESLeading/
ESInvert — i.e. every window move.  Transport.md L10 sides with Tolerances
("relative and three-valued").  As written, EpsSeries reintroduces exactly
the threshold-straddling silent misclassification (disease D1) that the
banded predicate exists to kill, and the two specs cannot both be
implemented.

Amendment (exact text).  In EpsSeries.md §2, REPLACE the `ESCoeffZeroQ`
definition paragraph with:

> - `ESCoeffZeroQ[c_, scale_] -> True | False | LOUD ERROR`
>   Thin wrapper, definitionally:
>   `ESCoeffZeroQ[c_, scale_] := Tolerances`NumericallyZeroQ[c, scale,
>   Tol["LaurentLeadTol"], $ESErrorContext, $LaurentLeadBandDecades]`.
>   All semantics (PossibleZeroQ short-circuit, symbolic -> False, zero
>   scale -> False, relative band, loud ambiguity error embedding the
>   caller's context) are Tolerances.md's; EpsSeries adds nothing.  The
>   ternary outcome is deliberate: a coefficient inside the ambiguity band
>   ABORTS the trim/lead/inversion that asked, naming the window and
>   context — it is never silently classified.

And in EpsSeries.md §8 ADD to test 11: "(c) a coefficient constructed inside
the ambiguity band (e.g. `10^-26` at scale 1 with tol `10^-24`, band 4
decades) makes `ESTrim` abort quoting `$ESErrorContext` (E5-class)".  ESSameQ
keeps binary behavior explicitly: ADD to ESSameQ's paragraph "ESSameQ uses
NumericallyZeroQ with matchTol and treats a band abort as False is FORBIDDEN
— a band hit inside ESSameQ is the same loud error; equality testing of
borderline data must be resolved by the caller, not guessed."

### 3. [blocker] Transport.md §2.6 / API.md §6 vs Solve.md §2 — the Solve interface Transport consumes does not exist; block structure / CouplingDepth / the marching particular have no owner

Problem.  Transport.md §2.6 step 1 calls
``DiffExp2`Solve`SolveChart[sys, chart, ...]`` returning "fundamental
LocalSolution basis + particular + metadata, including
`meta["CouplingDepth"]`", and step 4 sets `TWindow["CompleteMax"] =
expansionOrder - (meta["CouplingDepth"] - 1)`; §2.10 derives the probe
decrement from couplingDepth.  Solve.md exports "Exactly three" symbols —
SolveHomogeneous, SolveParticular, ODEResidualCheck — no SolveChart, no
CouplingDepth, no block/coupling concept anywhere.  API.md §6 separately
assigns "InitializeIntegrationSequence (:357-385) -> Solve.m (block ordering
+ MaxCouplingOrder for TWindow)" — also absent from Solve.md.  Deeper
unresolved design seam: Solve.md's recursion solves one d x d chart system
jointly; the old code's particulars (and the entire coupling-depth t-order
degradation, RewritePlan 3.1 TWindow comment) arise from BLOCK-SEQUENTIAL
solving where solved lower blocks become sources for upper blocks.  Solve.md
never says which policy holds, so (a) nobody produces the `particular` that
Transport's step 1/4 consumes, and (b) TWindow degradation has a consumer
(Transport) and no producer.  Note: with exact polynomial matrices
(denominators cleared, §3.6 step 1) the OLD degradation mechanism (truncated
matrix expansions) does not exist, so the t-orders of chained particulars
are NOT degraded in the new core — the RewritePlan 3.1 comment is old-code
wisdom whose premise the rewrite removes; this needs an explicit decision,
not silence.

Amendment (exact text).  ADD to Solve.md §2 a fourth export:

> ### 2.4 SolveChart
> ```
> SolveChart[chart_Association, req_Association]
>   -> <|"Basis" -> FundamentalSystem, "Particular" -> LocalSolution | None,
>        "CouplingDepth" -> _Integer|>
> ```
> The Transport-facing wrapper.  POLICY (binding): the chart system is
> solved BLOCK-SEQUENTIALLY along the exact block-triangular structure of
> ThetaMatrix (blocks = strongly connected components of the entry-sparsity
> graph, topologically ordered — the old InitializeIntegrationSequence role,
> DiffExp/MatrixLoading.m:357-385).  Per block: SolveHomogeneous on the
> diagonal block; the coupling submatrix times already-solved lower-block
> LocalSolutions (SectorSeries rational-multiply, exact polynomial
> coefficients) forms the SourceData for SolveParticular.  "CouplingDepth" =
> longest chain in the block dependency DAG.  "Particular" is None for a
> chart whose system is a single block (then transport is purely
> homogeneous).  T-ORDER NOTE (binding on Transport.md): because the
> coupling coefficients are exact polynomials, a particular built from a
> lower-block solution complete to TOrder is itself complete to TOrder —
> there is NO (couplingDepth-1) degradation in the new core.  Transport.md
> 2.6 step 4 is amended to `TWindow["CompleteMax"] = expansionOrder`, with
> the old MaxCouplingOrder discount recorded as a LessonsLedger waiver
> ("premise removed: matrices are exact polynomials, not truncated
> expansions"); the probe decrement formula of Transport 2.10 keeps
> couplingDepth as input unchanged (heuristic only).  If the orchestrator
> rejects the no-degradation analysis, the conservative fallback is to keep
> Transport's discount verbatim — but the decision must be recorded in the
> ledger either way; silence is not an option.

And in Transport.md §2.6 step 1 REPLACE "basis, particular, meta =
DiffExp2`Solve`SolveChart[sys, chart, ...]" with a cross-reference to
Solve.md §2.4 and the amended step 4.  Budget: +~35 lines in Solve.m, paid
by Solve.md §9 cut 3 (memoization plumbing) if needed.

### 4. [blocker] All specs §9 — the per-module line estimates sum to the hard ceiling BEFORE the functionality this review found missing; the budget does not add up

Problem.  RewritePlan 3.2 budgets sum to 3300 (Tolerances 100, Config 150,
EpsSeries 300, SectorSeries 400, Indicial 300, Solve 700, Transport 600,
Integrate 500, API 250; "Core total target: ~3.3k (ceiling 3.5k)" — i.e. the
listed figures ARE the 3.3k, so SectorSeries' "+100 Pade" reading is double
counting).  The specs' own §9 estimates: Tolerances 100, Config 150,
EpsSeries 305, SectorSeries 530 (60+25+90+100+110+70+35+15+25 — 130 over the
plan figure even on the generous 500 reading), Indicial 330 (self-declared
over), Solve 700, Transport 620 ("already tight"), Integrate 505, API 250.
Sum: 3490 of a 3500 ceiling — 10 lines of slack before implementation drift.
On top of that, defects 3, 8, 9, 10, 13 of this review add unbudgeted
functionality (~150 lines: SolveChart/blocks ~35, CombineLocalSolutions ~25,
tagged-power/analytic prefactor multiply ~40, bc tag parser ~30, error/print
primitives ~20).  As specified, the project busts its own hard ceiling.

Amendment (exact text).  ADD to each spec's §9 the following pre-authorized
cuts, to be taken at implementation time WITHOUT further sign-off, freeing
~200 lines to fund the missing operations:

> PRE-AUTHORIZED at M0 review (REVIEW-minimalism.md defect 4):
> - SectorSeries.md §9 cuts 1, 3, 5 are taken up front (Main-coordinates
>   option, multi-point evaluation form, public PadeEvaluate): -55; the
>   SectorSeries target is restated as 400 + 75 Pade = 475 TOTAL.
> - Indicial.md §9 cuts C-2 and C-3 are taken up front: -18 (C-1 trim-pass
>   stays IN until M2 evidence, because T-12/N-1b regression-pins it).
> - Transport.md §9 cut 3 (SavedCharts via caller callback): -30.
> - Tolerances.md: the dead adaptive-search constants are DELETED, not
>   folded (defect 8): -10.
> - EpsSeries.md §9 cuts 1-2 (binary ESAdd, drop ESZero): -13.
> - Integrate.md §9 cut 2 (b != 0 interior phase-paired path demoted to a
>   loud not-implemented error + ledger waiver) is pre-authorized but NOT
>   taken by default; it is the designated reserve (-50) if the post-fund
>   total exceeds 3450.
> Restated totals after cuts and additions: ~3455 of 3500.  Any module
> landing >10% over its restated figure stops and reports to the
> orchestrator BEFORE writing more code.

---

## MAJORS

### 5. [major] Tolerances.md §5 / EpsSeries.md §5 O-1 / Indicial.md §5.0 OQ-2 / SectorSeries.md §5.1 / Solve.md §5 / Transport.md §5 / Integrate.md §5 / API.md §5 — six different error mechanisms; the only nominated owner (API.m) is a dependency inversion

Problem.  Tolerances.md: "the shared DiffExp2 loud-error primitive (Abort
with printed context; the API.m spec owns its exact form)" — but Tolerances
is the BOTTOM module and cannot call API.m (its own §7: "May call: NOTHING").
EpsSeries invents `esError` (print + Abort).  SectorSeries uses
`Message[SectorSeries::tag] + Throw`.  Indicial assumes
`Throw[Failure["DiffExp2", payload], "DiffExp2"]` "or the shared loud-error
primitive if Config.md defines one" — Config.md defines none.  Solve uses "a
single error head, e.g. DiffExp2::solve".  Integrate uses "one helper that
prints/returns a structured Failure".  Every spec hedges with "align later";
M0 spec review is "later".  Same gap for verbosity printing: Config owns the
Verbosity/VerbosityDebug KEYS but no spec defines PrintInfo/PrintWarning
helpers, while Integrate §2.2.3 ("info-level Verbosity>=2 print"), API §5
(W1-W5, mandatory I1 print) and SectorSeries verbosity output all need them.

Amendment (exact text).  ADD to Tolerances.md §2:

> - `LoudError[module_String, tag_String, payload_Association] -> (never returns)`
>   THE library-wide error primitive.  Behavior:
>   `Throw[Failure["DiffExp2", Join[<|"Module" -> module, "Error" -> tag|>,
>   payload]], "DiffExp2"]` after printing one line
>   `"DiffExp2 <module> error <tag>: <payload as key=value>"`.  API.m's
>   entry points are the ONLY Catch["DiffExp2"] sites and convert the
>   Failure to a user-facing abort.  Mandatory payload keys where
>   applicable: chart center/label, sector tag (a,b,p), eps order, t order.
>   (+~12 lines, funded by defect 8's deleted constants.)

ADD to Config.md §2:

> - `PrintInfo[level_Integer, args__]` / `PrintWarning[args__]`: the only
>   verbosity-gated print helpers in DiffExp2 (gates: CFG["Verbosity"],
>   CFG["VerbosityDebug"]).  PrintWarning is never gated below visibility of
>   level 1.  No other module defines print helpers.  (+~8 lines.)

REPLACE the mechanism paragraphs of EpsSeries §5 (esError body), SectorSeries
§5.1 (Message+Throw — keep the per-tag NAMES as the `tag` argument),
Indicial §5.0, Solve §5 head, Transport §5 head, Integrate §5 head, and
Tolerances §5 head with: "All errors are raised via
``Tolerances`LoudError["<Module>", tag, payload]``; this module defines no
other error mechanism."  DELETE Tolerances.md's "the API.m spec owns its
exact form" clause.  Open questions O-1 (EpsSeries), OQ-2 (Indicial), and
SectorSeries' never-Quiet note are CLOSED by this amendment.

### 6. [major] Solve.md §3.7/§7 and Transport.md §2.6/§7 vs SectorSeries.md §2 — the LocalSolution linear-combination operations both consumers require are not in SectorSeries' public surface or budget

Problem.  Solve.md §3.7 step 2 needs "an EpsSeries-scalar × LocalSolution
operation; SectorSeries.m provides it" and §7 lists "EpsSeries-scalar ×
LocalSolution (compensation terms)".  Transport.md §2.6 step 4 assembles
"the chart's LocalSolution = basis.w + particular (SectorSeries algebra)"
and §7 lists "assemble/add LocalSolutions".  SectorSeries.md's public
symbols (2.1-2.9) contain NO addition of LocalSolutions and NO
multiplication by an EpsSeries scalar (MultiplyRational takes a closed-form
rational expression, not a window-carrying EpsSeries object).  Two modules
consume an operation no module provides or budgets.

Amendment (exact text).  ADD to SectorSeries.md §2:

> ### 2.10 `CombineLocalSolutions[weights_List, lss_List] -> LocalSolution`
> Exact linear combination Σ_i w_i · ls_i where each `w_i` is an EpsSeries
> (EpsSeries.md object; a plain exact number is auto-wrapped as a width-1
> series with window [0, +inf-equivalent: the partner's CompleteMax]).  All
> ls_i must share Center/ChartMap/Radius/Prescriptions exactly (else
> `::dims`).  Per sector: coefficient rows are EpsSeries-multiplied by w_i
> (window per ESTimes; Laurent weights shift kmin honestly — this is how
> Solve's compensation terms with polar γ(eps) and Transport's matched
> weights enter), same-tag sectors merged by ESAdd, result canonicalized
> (2.2).  EpsWindow = honest min over contributions; TWindow = min;
> ErrorEstimate = entrywise sum (I-10).  Budget +~25 lines (funded per
> defect 4).  Unit test t34: 2-sector combine with weights {1, 1/(2 eps)}
> reproduces a hand-computed window shift and merged coefficients exactly.

UPDATE SectorSeries §9 allocation accordingly; Solve.md §3.7 step 2 and
Transport.md §2.6 step 4 reference 2.10 by name.

### 7. [major] API.md §2.6 + §8 test 17 vs SectorSeries.md §2.5 — the IntegrateOverLine "Prefactor" contract is unimplementable: eps-dependent endpoint powers are not rational and no module can multiply them in

Problem.  API.md types the prefactor as `"PowerAtLower" -> v1(eps)` etc.,
and test 17 passes `PowerAtLower -> -1+eps` (the Beta-function pin).
Integrate.md's disposition table routes the prefactor to "SectorSeries.m's
rational-multiply (a-shift at the endpoint charts is exact tag arithmetic)".
But `(x - x0)^(-1+eps)` is not rational in (x, eps): SectorSeries
MultiplyRational raises `::nonrational` on it by contract (§2.5: "Anything
non-rational ... -> LOUD error"), and at INTERIOR charts the factor
`x^(b·eps) = e^(b·eps·Log x)` is an analytic non-rational function needing
an eps-series-of-t-series multiply that no spec defines.  The in-repo FT
pipeline only ever passes INTEGER powers (v1, v2 are propagator exponents;
FeynmanTrick/DiffExpIntegration.m:868-877, 1127-1128 — verified), for which
plain MultiplyRational suffices; the eps-dependent surface is creep that
happens to also be the only thing test 17 exercises.

Amendment (exact text).  In API.md §2.6 REPLACE the Prefactor option
description with:

> "Prefactor" -> <|"PowerAtLower" -> v1, "PowerAtUpper" -> v2,
>     "RationalFactor" -> r(x, eps)|> with v1, v2 EXACT RATIONALS (in-repo
> FT usage is integer powers; grep-verified the only call site,
> FeynmanTrick/DiffExpIntegration.m:868-877, passes propagator exponents).
> eps-DEPENDENT powers are out of scope v1 and raise a loud error naming the
> power and this section ("eps-dependent prefactor exponents are not
> supported in v1; fold the eps-power into the transported system or the
> combination").  With rational powers the whole prefactor is rational in
> (x, eps) and is applied per chart by SectorSeries MultiplyRational (center
> pole -> exact a-shift; other-chart poles -> analytic fold), per
> Integrate.md's disposition table.

And REPLACE API test 17 with:

> 17. test_api_integrate_over_line_beta: 1x1 EXACT system
>     dz_full.m = {{(-1+eps)/z - (-1+eps)/(1-z)... }} i.e.
>     A = (-1+eps)*(1/z + 1/(z-1)) so that f = z^(-1+eps) (1-z)^(-1+eps) is
>     the solution; transport with SaveExpansions across [0,1] anchored by
>     the interior value at z = 1/2; IntegrateOverLine over {0,1} with NO
>     prefactor; assert the LaurentValue equals
>     Series[Gamma[eps]^2/Gamma[2 eps], {eps, 0, CompleteMax}]
>     coefficientwise to 10^-25 (leading 2/eps).  Exercises: pole cells at
>     both endpoint charts, regular interior charts, assembly — without the
>     out-of-scope prefactor surface.

If the orchestrator instead wants eps-dependent powers in v1, the required
amendment is a new SectorSeries export `MultiplyTaggedPower[ls, alpha, beta]`
(per-sector exact tag shift a->a+alpha, b->b+beta, legal ONLY when the
power's base point is the chart center; ~12 lines) PLUS an interior-chart
analytic-factor multiply (~30 lines) — both must then be specified, budgeted
and tested in SectorSeries.md, not assumed.

### 8. [major] Config.md §3.2 + Tolerances.md §2 + API.md §2.3 vs Transport.md §6.1/§9 — AccuracyGoalValidate: three specs ship the adaptive-search machinery that the fourth deletes; six exported constants are dead on arrival

Problem.  Transport.md does NOT port the adaptive expansion-order search
(§6.1 table: "NOT ported; replaced by E11"; §9 cut 1: "ALREADY cut by this
spec").  Config.md nevertheless specifies "AccuracyGoalValidate" as a kept
key with full "Before"-search/"After"-redo semantics (schema row), a
cross-field error E8, and test C14.  Tolerances.md exports six constants
that exist only for that machinery ($ExpansionOrdersAveraging,
$ExpansionOrderIncrease, $ExpansionOrderDecrease,
$ExpansionOrderIncreaseValidate, $DigitsSurplusDecrease) plus
$SafetyExpansionSubtract (whose consumer, the matrix-truncation error
read-off LineSegmentation.m:109-113, also does not survive — no spec uses
it), and pins them in test T12.  API.md §2.3 asserts the machinery exists
("AccuracyGoalValidate 'Before' adaptive expansion-order search ... and
'After' redo ... implemented in Transport.m").  A kept config key whose
specified behavior nothing implements is precisely the silent-no-op disease
(D1/FEC class) this rewrite exists to kill.  Grep evidence: ZERO in-repo
consumers of AccuracyGoalValidate outside DiffExp/ (Reference/Examples,
Tests, Scripts, FeynmanTrick all clean — verified).

Amendment (exact text).  In Config.md §3.2, MOVE "AccuracyGoalValidate" to
the dropped-keys table with the row:

> | "AccuracyGoalValidate" | the adaptive expansion-order search (old
> Transport.m:776-841, 1191-1213) is not ported; validation is ALWAYS ON
> when "AccuracyGoal" is numeric (Transport.md E11: per-segment probe error
> <= 10^-DigitsNeeded, error naming the segment and the suggested
> ExpansionOrder).  LessonsLedger waiver "expansion-order adaptive search"
> required. | State.m:108 |

DELETE Config E8 and retarget C14 to assert the dedicated dropped-key error.
In Tolerances.md §2 DELETE $ExpansionOrdersAveraging, $ExpansionOrderIncrease,
$ExpansionOrderDecrease, $ExpansionOrderIncreaseValidate,
$DigitsSurplusDecrease, $SafetyExpansionSubtract (KEEP $MinExpansionOrder —
Config's ExpansionOrder validator references it); prune T12 accordingly and
DELETE §9 cut option (1) (nothing left to fold).  In API.md §2.3 DELETE the
sentence claiming the Before/After machinery and substitute "AccuracyGoal
numeric => always-on per-segment validation (Transport.md E11)".  ALSO (same
disease, same fix pattern): Config.md must state the meaning of
`"EstimateError" -> True` vs `"Fast"` — no spec distinguishes them; ADD to
the schema row: "True is accepted as an alias of \"Fast\" (one probe
implementation exists in DiffExp2; the old slow/fast split is waived in the
ledger)."

### 9. [major] API.md §2.1/§2.8 vs Config.md §1/§3.2 — "Config.m forwards to LoadSystem" is an upward call that violates the module order Config.md itself declares

Problem.  API.md: "LoadConfiguration with MatrixDirectory still works
(Config.m forwards to LoadSystem)" (§2.1) and "LoadConfiguration with
MatrixDirectory triggers LoadSystem" (§2.8).  Config.md: "Config performs NO
loading" (MatrixDirectory schema row), the auto-load is "WAIVED", and
Config.m may call ONLY Tolerances (§7).  Config calling API.m's LoadSystem
is a dependency inversion; as specified, the classic call pattern either
breaks (Config doesn't forward) or the order breaks (it does).

Amendment (exact text).  In API.md §2.8 REPLACE the forwarding sentence
with:

> API.m re-exports LoadConfiguration/UpdateConfiguration as thin wrappers:
> `DiffExp2`LoadConfiguration[rules]` calls
> ``Config`LoadConfiguration[rules]`` and THEN, iff the resulting
> configuration has a non-empty "MatrixDirectory", calls
> `LoadSystem[CFG["MatrixDirectory"]]`.  The forwarding lives in API.m
> (which may call both Config.m and its own LoadSystem); Config.m itself
> performs no loading (Config.md §1, MatrixDirectory row).  Classic scripts
> that call LoadConfiguration via the DiffExp2` context therefore keep the
> old auto-load behavior; code that calls ``DiffExp2`Config`
> LoadConfiguration`` directly gets configuration only.

And in API.md §2.1 REPLACE "(Config.m forwards to LoadSystem)" with
"(the API.m re-export wrapper forwards; see 2.8)".

### 10. [major] Transport.md §2.7/§8 T15 vs Tolerances.md §2 — SnapTol === RankTol makes Transport's pivot gray zone empty and T15 ill-posed; the snap/rank pair duplicates NumericallyZeroQ

Problem.  Tolerances defines `SnapTol[wp] = RankTol[wp] =
10^(-Floor[wp/2])` — identical values.  Transport 2.7 specifies "a pivot
with |pivot| < snapTol*scale is snapped to exact 0; snapTol*scale <= |pivot|
<= rankTol*scale is the gray zone = LOUD ERROR E7" — with equal tolerances
the gray zone is a single point and E7 is unreachable; T15 ("pivot
constructed at geometric mean of snapTol and rankTol scales") constructs the
boundary point itself.  Moreover this hand-rolled three-zone logic is a
second implementation of exactly what Tolerances' NumericallyZeroQ band
already provides, violating Tolerances §2's own rule that rank decisions
"MUST go through NumericallyZeroQ ... never a raw comparison".

Amendment (exact text).  In Transport.md §2.7 REPLACE the
"Numerically-zero leading rows/pivots" bullet with:

> - Pivot/leading-coefficient zero classification: EXACTLY
>   ``Tolerances`NumericallyZeroQ[pivot, scale, Tol["RankTol"], <chart
>   label + operation>]``.  True -> the pivot is exact 0 (snap); False ->
>   nonzero; the ambiguity band IS the gray zone and aborts loudly (E7 is
>   the band error, not a separate mechanism).  SnapTol plays NO role in
>   matching solves; it is reserved for coordinate/endpoint snapping.

REPLACE T15 with: "T15 `test_gray_zone_loud`: pivot constructed inside the
NumericallyZeroQ band around Tol[\"RankTol\"]*scale (e.g. exactly
rankTol*scale): the band error fires naming both band edges and the chart
label."  In Transport.md §2.8 the post-shift leading-coefficient check gets
the same one-line replacement ("assert via NumericallyZeroQ at RankTol; band
hit = loud").  In Tolerances.md §3's consumer table, ADD RankTol-via-
NumericallyZeroQ to Transport's row and REMOVE the matching-solve role from
SnapTol's description in §2.

### 11. [major] Tolerances.md §2 (SnapTol derivation) + Transport.md E14 + API.md W5 — the default snap tolerance 10^(-wp/2) cannot perform its stated user-facing job; numeric endpoints near singularities silently become regular points

Problem.  SnapTol[500] = 10^-250.  API W5 snaps "a numeric target ... within
Tolerances.snapTol of an exact root"; Transport E14 errors when `to` is
"within snapTol of a singularity but not exactly equal".  A user-supplied
`to = N[1/4, 30]` (distance ~10^-31 from the exact singularity 1/4) is
NEITHER snapped NOR flagged: it silently becomes a regular endpoint sitting
10^-31 inside a chart of radius 10^-31 — segment explosion or garbage, the
exact silent-degradation class the plan bans.  The old absolute 10^-10
(State.m:123) existed precisely for user-numeric inputs; deriving the
replacement from wp instead of from the INPUT's accuracy misses the point of
the lesson.

Amendment (exact text).  ADD to Tolerances.md §2:

> - `$NearSingularityGuardDecades = 6` (new).  Used by Transport/API
>   endpoint classification: an INEXACT target/bound `to` with
>   `|to - root| < 10^(-Floor[Accuracy[to]/2])` for some exact singularity
>   root is SNAPPED to the root with warning W5 (the input's own accuracy,
>   not wp, sets the snap scale — a 30-digit 0.25 snaps to 1/4); an inexact
>   `to` with `10^(-Floor[Accuracy[to]/2]) <= |to - root| <
>   10^(-$NearSingularityGuardDecades)` is a LOUD ERROR ("target is
>   suspiciously near singularity <root>; pass the exact value, or an
>   exact offset, to confirm intent").  Exact inputs are never snapped and
>   never guarded (exact arithmetic decides).  SnapTol (10^(-Floor[wp/2]))
>   remains for snapping COMPUTED full-precision quantities only
>   (chart-map endpoint values, factor-root dedup pre-filter).

REPLACE Transport.md E14's condition with the guard rule above (the
"within snapTol" clause is deleted), and API.md W5's clause "within
Tolerances.snapTol" with "per the input-accuracy snap rule
(Tolerances.md $NearSingularityGuardDecades entry)".

### 12. [major] API.md §7 vs Integrate.md §6 vs SectorSeries.md §2.4.1 — the endpoint-limit drop rule (RewritePlan 3.3's first boundary case) has no owner

Problem.  API.md §7: "EndpointLimit -> SectorSeries (combination) +
Integrate (drop rule/divergence checks)".  Integrate.md §6 export table:
"`EvaluateLimitAtSingularity` ... NOT Integrate.m: endpoint limits are
SectorSeries evaluate-with-branch-rule / API.m EndpointLimit ...
Cross-reference only."  SectorSeries `::originlimit` (§2.4.1) explicitly
refuses t=0 evaluation and points "to the EndpointLimit semantics in
API.m/Integrate.m (the b≠0 drop rule and the divergence checks live THERE,
not here)".  API.m's budget (§9: EndpointLimit ~20 lines) rests on "API.m
contains no mathematics".  Three pointers, zero implementations: the
(0,0,0)-sector readout, the exact b≠0 drop, and the b=0/a<=0/p>0 divergence
gate of RewritePlan 3.3 are specified nowhere as owned code.

Amendment (exact text).  ADD to Integrate.md §2 a fourth public symbol (and
count it in §9, +~20 lines, inside the existing chart-level block):

> ### 2.4 `EndpointSectorLimit[ls_Association, direction_:1] -> {EpsLaurent per component}`
> The RewritePlan 3.3 limitUpper/limitLower primitive, owned HERE because it
> shares the 2.2.3 gate: lim_{t->0} of the assembled LocalSolution.  Per
> sector: b != 0 -> dropped EXACTLY (dimreg convention of 2.1.2, applied to
> the value: t^(b eps)|_{t=0} := 0); b = 0 with a < 0 (any p), or a = 0 with
> p > 0, present at any n with a + n <= -1 resp. (a+n = 0, p > 0) -> run the
> 2.2.3 merged-coefficient cancellation gate per eps order; survivor ->
> LOUD ERROR (API E26 payload: component, tag, eps order); b = 0, a = 0,
> p = 0 -> the n = 0 column IS the limit.  `direction` resolves nothing here
> (the limit value is branch-independent for the surviving terms) but is
> validated against the chart Prescriptions for consistency (E3 if a b != 0
> sector would have needed a sign AND no prescription exists — the check
> stays so configuration gaps surface even on the drop path).
> API.m's EndpointLimit = combination (SectorSeries.CombineLocalSolutions /
> MultiplyRational) + this function.  Unit test 23: the three-sector mix
> x^0 (constant 7/3), x^(-1+eps), x^(2 eps) returns exactly 7/3 at eps^0;
> flipping the middle tag to (-1, 0, 0) with nonzero merged coefficient
> raises the divergence error.

UPDATE API.md §7's EndpointLimit line to "SectorSeries (combination) +
Integrate.EndpointSectorLimit (drop rule + divergence gate, Integrate.md
2.4)" and Integrate.md's §6 EvaluateLimitAtSingularity row to "ABSORBED ->
EndpointSectorLimit (2.4)".

### 13. [major] API.md §2.2/§9 vs SectorSeries.md §2 — the boundary-condition closed-form/asymptotic -> sector-tag parser is assigned to SectorSeries by API's budget but absent from SectorSeries' spec

Problem.  API.md §9 counting guide: "PrepareBoundaryConditions ~60
(entry-type dispatch; tag parsing helpers live in SectorSeries)"; §2.2
requires parsing `(-1/t)^(1+3eps) /. t -> -1/x` into exact tags {a,b,p},
plus the x-SeriesData form `0 + O[x]^(1/2)` (test 6 asserts tags {1,0},
{1,1}, {1,2}, {1,3}).  SectorSeries.md §2 has no such symbol, §9 no such
lines.  ~30 lines of real parsing (Power/Times/Log pattern walk + eps-affine
exponent extraction + loud failure for unparseable forms) are budgeted
nowhere and error-contracted nowhere (what happens on `x^Sqrt[eps]`?).

Amendment (exact text).  ADD to SectorSeries.md §2:

> ### 2.11 `ParseTaggedPower[expr_, var_Symbol, epsSym_Symbol] -> <|"a", "b", "p", "Coefficient"|> | $FailedParse`
> Boundary-ingestion helper (sole caller: API.m PrepareBoundaryConditions).
> Recognizes c * var^(a + b*epsSym) * Log[var]^p products with exact a, b
> (affine-in-eps exponent REQUIRED: a non-affine exponent is a LOUD error
> naming the exponent — same I1 discipline as Indicial, never a numeric
> fit), p integer >= 0; rewrites Log[k*var] -> Log[k] + Log[var] first (old
> Transport.m:120 rule).  Returns the inert marker $FailedParse for
> expressions that are not of this product form — the CALLER (API.m)
> decides between the documented minimum contract (old eps-expansion into
> Logx polynomials, old Transport.m:91-120) and an error; $FailedParse is
> data for that documented branch, not a silent fallback.  Budget +~25
> lines (funded per defect 4).  Unit t35: the banana entry-3 forms of
> API.md test 6 parse to tags {1,0},{1,1},{1,2},{1,3} exactly;
> `x^Sqrt[eps]` errors naming the exponent.

UPDATE SectorSeries §9 and API.md §9 item 4 (which currently offers to
"move the parser into SectorSeries" as a CUT — it is the baseline, not a
cut; delete that item).

### 14. [major] Transport.md §2.8/Q2 vs Indicial.md §2.5/§3.6 — Transport requires per-chart "EpsDegenerateFamilies" that Indicial explicitly declines to compute

Problem.  Transport 2.8: recombination candidates are "restricted to
families flagged eps=0-degenerate by Indicial.m metadata
(`chart["EpsDegenerateFamilies"]`)".  Indicial 2.5: the eps->0 eigenvector
collision "is NOT Jordan structure at generic eps and is NOT detected here
... Division of labor: Indicial = generic-eps structure; Transport =
eps -> 0 limits."  IndicialData (§3.6) has no such field.  Transport's own
Q2 fallback ("E5 firing at the match point") converts every legitimate
recombination case into an abort.  Note JointSolve == True (>= 2 distinct
b) is necessary but NOT sufficient for eps=0 degeneracy (diag(eps, 2 eps)
is a pseudo family with independent eps=0 eigenvectors), so Transport
cannot derive the flag from the existing Family record alone.

Amendment (exact text).  ADD to Indicial.md §2.6 step 6 and §3.4:

> - `"EpsZeroDegenerate" -> True | False` per Family: True iff the Jordan
>   type of the eps -> 0 residue R(0) restricted to the family's root
>   subspace differs from the eps-generic type — computed EXACTLY by
>   running the §2.5 nullity sequence once more on R /. eps -> 0 for the
>   family's eps->0 eigenvalue(s) and comparing block multisets (~10
>   lines; all-exact, no tolerance).  This is the tag-driven input to
>   Transport.m's RecombineBasis (Transport.md 2.8); Transport performs NO
>   numeric degeneracy detection.  Unit T-23: diag(eps, 2 eps) ->
>   EpsZeroDegenerate False; {{0,1},{0,2 eps}} (the log x class, Transport
>   T12 matrix residue) -> True.

UPDATE Transport.md 2.8 first bullet to cite the field name
"EpsZeroDegenerate" (singular per family) and DELETE Transport Q2.

### 15. [major] Solve.md §3.2/OQ1 vs Indicial.md §3.3-3.6 — the producer/consumer data contracts diverge in content, not just key names: V/VInv are consumed but never produced, and ChartSystem has no builder

Problem.  Solve's ChartSystem requires per family exact `V(eps)`,
`V^(-1)(eps)`, `J`, collision keys {From, To, Offset, DeltaB},
CollisionDepth, plus chart fields (ChartVar, ChartMap, Radius,
Prescriptions, PreparedCacheKey) and a certified "det V(eps) vanishes at
eps = 0 only through the collision factors" property.  Indicial produces
EigRecords with "Chains" (vectors, not matrices), no inverse, collisions
keyed {n, LowerIdx, UpperIdx, DeltaB, Type}, no CollisionDepth, and
ChartIndicial's input ChartRef carries only Name/Center/Variable — so no
module ever assembles the record Solve consumes.  Both specs file this
under "open question"; M0 spec review is where it must close, and "key
names to be reconciled" understates a missing inverse computation and a
missing certification.

Amendment (exact text).  ADD to Solve.md §3.2 (replacing the OQ1 hedge):

> RESOLUTION (M0 review): (i) Solve.m consumes Indicial's Family/EigRecord
> shapes VERBATIM (Docs/specs/Indicial.md 3.3/3.4 are normative; collision
> keys n/LowerIdx/UpperIdx/DeltaB/Type — Solve's From/To/Offset names in
> this spec are renamed to those).  (ii) Solve.m BUILDS the per-family
> spectral frame itself, once per chart: V = the matrix whose columns are
> the family's chain vectors in spec order, VInv = exact Inverse[V]
> (Q(alpha)(eps) arithmetic; +~15 lines, inside the §3.6 step-2 budget),
> and certifies at that point that det V's eps -> 0 vanishing order equals
> the count implied by the recorded LaurentShift collisions (the former E2
> third check; failure stays E2).  (iii) CollisionDepth = longest directed
> chain in the family's collision DAG, computed here (~5 lines).
> (iv) ChartSystem is ASSEMBLED BY TRANSPORT.M (the only module holding
> both IndicialData and chart geometry): ChartSystem = chart geometry
> fields (Transport.md 3.2) + IndicialData fields, with
> "ThetaMatrix"/"Gauge"/"GaugeInverse"/"Residue" lifted from
> IndicialData["Reduction"].  Transport.md 2.6 step 1 gains the sentence
> "TransportLine builds the ChartSystem record from SegmentLine's Chart and
> ChartIndicial's IndicialData before calling SolveChart."

DELETE Solve OQ1 and Indicial OQ-5 (subsumed); Indicial keeps emitting
exactly what its §3 already specifies.

### 16. [major] Config.md §2 / EpsSeries.md §4 I-3, §7 / SectorSeries.md §3.2, §7 — the canonical eps symbol and the declared-indeterminates registry are consumed from Config, which defines neither

Problem.  EpsSeries I-3 validates eps-freeness against "the symbol set
[that] comes from Config.m's validated accessor"; SectorSeries §7 reads
"the pinned variable and eps symbols, declared indeterminates (I-5)" from
Config.  Config.md's schema and public symbols contain NO eps-symbol entry
and NO indeterminates registry (API.md's BoundaryConditions carries
"Indeterminates" per call, and SectorSeries cannot call API).  As written,
SectorSeries I-5 ("linear combinations of DECLARED indeterminate symbols")
is unvalidatable and EpsSeries I-3 reads a nonexistent key.

Amendment (exact text).  ADD to Config.md §2:

> - `EpsSymbols[] -> {Global`eps, \[Epsilon]}` and
>   `CanonicalEps[] -> Global`eps`: constants (not schema keys — not user
>   configurable), THE library-wide answer to "which symbols are the
>   regulator" (RewritePlan 3.2 API.m: both accepted).  EpsSeries I-3 and
>   SectorSeries validation read these; API.m normalizes inputs to
>   CanonicalEps[] at ingestion.

ADD to SectorSeries.md §3.1/3.2: "LocalSolution carries an OPTIONAL
`"Indeterminates" -> {sym..}` key (default {}; populated by API.m from the
BoundaryConditions record).  I-5's 'declared indeterminate symbols' means
EXACTLY this list: a symbolic coefficient entry whose symbols are not all
in the object's own Indeterminates list is `::badcoeff`.  Config holds no
indeterminate state."  UPDATE EpsSeries §7's Config bullet to cite
EpsSymbols[]/CanonicalEps[] by name, and DELETE its claim that Indicial
consumes EpsSeries (§7 "Consumed by" list — Indicial.md §7 says EpsSeries
is NOT USED; Indicial is right, the consumer list is wrong).

### 17. [major] Tolerances.md §3 consumer table vs every consumer spec — the binding table is wrong in five places and two consumed names don't exist

Problem.  The table is declared binding ("anything else is a spec violation
found in review").  Violations found: (i) EpsSeries consumes matchTol
(ESSameQ; EpsSeries §7) — table allows only LaurentLeadTol + ChopFloor;
(ii) the table assigns EpsSeries "ChopFloor (coefficient chop after
add/mul/div)" — EpsSeries I-6 FORBIDS chopping stored coefficients
entirely (the EpsSeries position is the principled one); (iii) the table
assigns Indicial "ChopFloor for numeric matrix entries" — Indicial.md §7
declares Tolerances "NOT USED. Deliberate and load-bearing" and I-8/E1
reject numeric entries outright (Indicial is right); (iv) Transport
consumes `chopReserve` (§2.5, E3) — no such Tolerances export exists;
(v) SectorSeries' `::geomambiguous` (§2.5 step 3) uses "the Tolerances
guard" — no named export exists, violating Tolerances' own F4 ("a
different semantics ... gets a new NAMED export").

Amendment (exact text).  In Tolerances.md §3 REPLACE the consumer table
rows for EpsSeries, Indicial, Transport, Integrate with:

> - EpsSeries.m: `LaurentLeadTol` (window/lead classification via
>   NumericallyZeroQ only), `MatchTol` (ESSameQ).  NO chopping: EpsSeries
>   never alters stored coefficients (its I-6).
> - Indicial.m: NOTHING.  (Exact algebra only; its §7 statement governs.)
> - Solve.m: `MatchTol`, `RankTol`, `LaurentLeadTol` (assembly-time
>   log-member trim, defect-1 amendment), `ResidTol`.
> - Transport.m: `MatchTol`, `RankTol` (via NumericallyZeroQ),
>   `ResidTol`, `SnapTol` (computed-value snapping only), `ChopReserve`,
>   `$SafetyDigits`, `$InputPrecisionFactor`, `EvalErrorSeriesDecrease`,
>   `$MinExpansionOrder`, `$NearSingularityGuardDecades`.
> - Integrate.m: `LaurentLeadTol` (cancellation gate + result trim),
>   `SnapTol` (tiling/endpoint snapping), `MatchTol` (additivity
>   spot-check I6).

ADD to Tolerances.md §2: "`ChopReserve[wp_Integer, chopDigits_Integer] ->
wp - chopDigits` — the digit reserve between working precision and the
chop floor; Transport.md E3 requires DigitsNeeded + ChopReserve <= wp (old
implicit reserve: WP - ChopPrecision = 250, State.m:109,135)."  ADD:
"`GeomGuardTol[wp_Integer] -> 10^(-Floor[wp/2])` — the
radius-vs-pole-modulus comparison guard consumed by SectorSeries
`::geomambiguous` and Transport geometry asserts; comparisons whose
difference magnitude falls below it are LOUD errors, not decisions."
UPDATE SectorSeries §2.5 step 3 and §7 to name GeomGuardTol.

### 18. [major] SectorSeries.md §2.3 vs Transport.md §2.9/E8 vs Integrate.md §2.2.4/E3 — the prescription-sign derivation is specified three times, and the Transport/Integrate dependency statements contradict each other

Problem.  SectorSeries.ChartImSign (§2.3) is the ported sign-derivation
(AnalyticContinuation.m:18-68).  Transport E8 re-states the same derivation
rules inline ("derivation rules per old AnalyticContinuation.m:45-68 ...")
and Chart["CrossSign"] is "derived" without naming the function; Integrate
E3 derives the sign a third time ("from ls['Prescriptions']; missing/
ambiguous => E3").  Three implementations of one rule is how conventions
drift (the exact disease the 9aeb300 lesson documents).  Additionally:
Transport.md §7 says "Integrate may call Transport's exported
CrossingOperator" while Integrate.md §7 says Integrate "MUST NOT call ...
Transport.m (same level)" — and Integrate indeed reimplements the pairing
phases internally (§2.2.4).  Finally Config.md §3.3 still specifies the
sqrt auto-prescription union ("LoadSystem can union in sqrt-derived
auto-prescriptions, old State.m:227-230") while API.md E6 rejects every
sqrt-bearing matrix in v1 — the union source is empty by construction.

Amendment (exact text).  (i) In Transport.md §2.9 and E8 ADD: "Sign
derivation is NOT reimplemented here: Chart['CrossSign'] :=
SectorSeries`ChartImSign applied to the chart's Prescriptions record;
conflict/missing surfaces as ChartImSign's `::branchconflict` /
`::branchmissing` wrapped into E8's payload."  (ii) In Integrate.md §2.2.4
b != 0 bullet and §2.2.5 ADD: "σ := SectorSeries`ChartImSign[ls]; E3 fires
iff it returns None while a phase is needed, or throws conflict."
(iii) In Transport.md §7 DELETE "(sibling; Integrate may call Transport's
exported CrossingOperator)" — Integrate.md's prohibition stands; the shared
convention is pinned by SectorSeries' sigma rule (2.4.1) plus the
cross-module parity tests (Transport T8, Integrate test 14), not by a
sibling call.  (iv) In Config.md §3.3 REPLACE the sqrt-union sentence with:
"The old sqrt auto-prescription union (State.m:227-230,
MatrixLoading.m:181-190) has no v1 source — LoadSystem rejects sqrt
matrices (API.md E6).  The user/effective list separation is kept (it is
~2 lines and the v1.1 sqrt ledger item lands on it), with a comment that
the effective list equals the user list in v1."

### 19. [major] API.md §2.3 (singular-endpoint "SeriesValues") + §10 Q3 — an unowned tag->Logx-series converter is kept alive for a consumer already scheduled for retargeting

Problem.  At a singular endpoint API keeps "SeriesValues" as
"per-(integral, eps-order) x-series evaluation forms ... (Normal-izable,
Logx-bearing)".  Producing those from sector-native data requires
collapsing exact tags into per-eps-order Logx towers — a converter no spec
owns, in the one direction (exact -> collapsed) the whole design calls
disease D2; its sole consumer, test_singular_endpoint.m, is ALREADY
designated for M6 retargeting to LocalSolution (API compat table row 18),
and API Q3 admits the forms can then be dropped.  Surface creep with an
unbudgeted implementation cost and an anti-pattern shape.

Amendment (exact text).  In API.md §2.3 REPLACE the compatibility sentence
with:

> At a singular endpoint, "SeriesValues" -> Missing["SingularEndpoint"]
> and the typed result is "LocalSolution" (RewritePlan 3.1).  The old
> Logx-bearing SeriesData tables (old Transport.m:1084-1085) are NOT
> reproduced: collapsing exact tags into per-eps-order Logx towers is the
> D2 representation this library exists to delete, and the only in-repo
> consumer (Tests/test_singular_endpoint.m:119-176) is retargeted to
> LocalSolution at M6 (compat table row 18).  Anyone needing display forms
> evaluates the LocalSolution or formats the Sectors list.

DELETE API Q3.  (If the orchestrator insists on the compat forms for the
M4 window, the converter must be added to SectorSeries.md with budget and
tests — it may not exist as unowned API glue.)

### 20. [major] API.md §2.3 (symbolic indeterminates) — keeping the old global UsePade mutation "for compat" reintroduces the exact forbidden fallback SectorSeries F9 bans

Problem.  API.md: "UsePade is force-disabled with a LOUD warning ... (old
Transport.m:565-569 prints Info and mutates global config — the global
mutation is kept for compat ...)".  SectorSeries F9 forbids precisely this
("silently mutate global config (old: TurnOffPade ...)"), and Config.md's
design makes UpdateConfiguration the only sanctioned mutation path.  A
transport call that flips a config key as a side effect leaves the NEXT
unrelated call running without Pade — hidden cross-call state.

Amendment (exact text).  In API.md §2.3 REPLACE "(... the global mutation
is kept for compat but the result metadata is new)" with:

> The disable is PER CALL: API.m passes an explicit "UsePade" -> False
> override down the Transport/SectorSeries call chain for this transport
> only (SectorSeries' evaluate already takes the option, SectorSeries.md
> 2.4); the configuration store is NEVER mutated as a side effect
> (SectorSeries.md F9; Config.md F-g discipline).  The warning W3 and
> "PadeDisabled" -> True result metadata are the visible record.

### 21. [major] SectorSeries.md §2.4.2/§2.9 vs §8 t14 — the Pade order formula and the load-bearing accuracy test contradict each other (and the old code's formula has an off-by-one the spec silently inherits)

Problem.  §2.4.2: `m = Floor[(numberOfCoefficients + 1)/2]`; §2.9:
`m = Floor[(N + 2)/2]` for coefficients {c0..cN} — both give m = 6 for the
11 coefficients of t14, i.e. [6/6].  t14 asserts "[5/5] reproduces 4/7 to
better than 10^-(wp-10)".  These cannot both hold; worse, [6/6] BUILT FROM
AN 11-COEFFICIENT POLYNOMIAL (the old code applies PadeApproximant to the
Normal-ized polynomial — verified DiffExp/Pade.m:35,42 with maxPadeOrder
computed from the one-past-end SeriesData slot, which overcounts the
coefficient count by one) treats the unknown orders 11-12 as exact zeros
and does NOT reproduce 1/(1+t): its error at t = 3/4 is ~the truncation
tail, i.e. no better than the plain sum — t14's accuracy claim fails, and
with it the documented rationale "diagonal Pade of a rational function is
exact" (L8, the R6 benchmark backbone).

Amendment (exact text).  In SectorSeries.md §2.4.2 and §2.9 REPLACE the
order formula with:

> diagonal order [m/m] with `m = Floor[numberOfCoefficients/2]` (so 2m+1 <=
> numberOfCoefficients: the approximant is determined ENTIRELY by known
> coefficients; for {c0..cN}, m = Floor[(N+1)/2]).  DELIBERATE one-slot
> deviation from old Pade.m:35, whose SeriesData arithmetic
> (`(a[[5]]-a[[4]])/a[[6]] + 1`) overcounts by one and silently consumes a
> phantom top order — LessonsLedger entry required ("Pade order off-by-one:
> old code built [m/m] needing 2m+1 coefficients from 2m known ones,
> treating the truncation boundary as exact zeros").  t14's [5/5] from 11
> coefficients is then consistent, and the rational-function exactness
> claim holds.

### 22. [major] Transport.md §8 T12 — the recombination test's incoming data is not a solution of its own system and the asserted weights are unreachable

Problem.  For A = {{0, 1/x}, {0, 2 eps/x}} the test's basis S2 =
(x^(2eps), 2 eps x^(2eps)) is a solution (verified), and B2 = (S2 - S1)/
(2 eps) = ((x^(2eps)-1)/(2 eps), x^(2eps)).  The stated incoming data
"((x^(2eps)-1)/(2 eps), x^(2eps)-1)" is NOT a solution of A (check:
f2 = x^(2eps)-1 gives f2' = 2 eps x^(2eps)/x != 2 eps f2/x), and no
constant weights reproduce it in the basis span — the asserted result
"MatchWeights == (0, 1) to matchTol" is unsatisfiable, so the test would
fail (or worse, pass only against a buggy implementation).

Amendment (exact text).  In Transport.md T12 REPLACE the incoming-data
clause with:

> for incoming data = the eps-regular solution B2 itself, i.e.
> ((x^(2eps)-1)/(2 eps), x^(2eps)) evaluated at t_m = 1/3 (eps^0 value
> (Log[1/3], 1)), MatchWeights == (0, 1) to matchTol, with weight window
> starting at eps^0 (NOT eps^-1).

### 23. [major] SectorSeries.md §3.4 / Integrate.md §3.4 / API.md §3.8 vs EpsSeries.md §3.1 — four different eps-Laurent value shapes; the reconciliation every spec defers is due now

Problem.  EpsSeries (normative owner):
`<|"EpsWindow" -> <|"Min","CompleteMax"|>, "Coeffs" -> {...}|>`.
SectorSeries EpsLaurentValue: `<|"MinPower","CompleteMax","Coefficients",
"PadeFallbacks","TailEstimates"|>`.  Integrate EpsLaurent:
`<|"Min","CompleteMax","Coefficients"|>`.  API LaurentValue:
`<|"MinPower","Coefficients"|>` (+"CompleteMax" on output).  Each spec
files a Q1-class open question pointing at the others; implementation
agents reading one spec each will ship four shapes.

Amendment (exact text).  ADD to EpsSeries.md §3.1 (and cite from the other
three, deleting SectorSeries Q1, Integrate Q1, and the Integrate §3.4
"assumed shape" block):

> CANONICAL SHAPE RULING (M0 review): every INTERNAL eps-Laurent value in
> DiffExp2 is this EpsSeries association, verbatim — SectorSeries'
> EvaluateLocalSolution returns `<|"Value" -> EpsSeries-per-component (a
> length-Ncomp list of EpsSeries, or one EpsSeries of vector coefficients
> per its 3.2 choice), "PadeFallbacks" -> ..., "TailEstimates" -> ...|>`;
> Integrate consumes/produces EpsSeries objects directly.  The ONLY other
> shape in the library is API.md's user-facing LaurentValue
> (`<|"MinPower","Coefficients","CompleteMax"|>`, FT-compatible per
> Scripts/run_ft_stepwise.m:47-53), produced exclusively by an API.m
> boundary converter (ESWindow/ESCoefficientList -> LaurentValue, ~4
> lines) — the ESToExpression usage rule applies to it analogously: no
> core module constructs or consumes LaurentValue.

### 24. [major] Transport.md E8 (structured exception) vs API.md E19 — the AbortOnAnalyticContinuationFail = False semantics differ in scope between the two specs

Problem.  Transport E8: with the flag False, a degraded (honest-partial)
result is permitted ONLY "iff ... the chart is the FINAL one in the plan";
anywhere else it errors.  API E19: with the flag False, "degrade to
WARNING + 'MultivaluedFail' -> True in the result" — no final-chart
restriction.  The FT use case (FeynmanTrick/DiffExpIntegration.m:411-414)
gets different behavior depending on which spec the implementer read; the
pentagon triage (M0 task 16, RewritePlan R8) depends on which one is real.

Amendment (exact text).  In API.md E19 REPLACE the False-branch sentence
with:

> With config "AbortOnAnalyticContinuationFail" -> False the degraded path
> exists ONLY when the offending chart is the FINAL chart of the plan
> (Transport.md E8's structured exception — the FT singular-endpoint use
> case, where the unprescribed chart IS the target): warning +
> "MultivaluedFail" -> True + completeness metadata marking every order
> incomplete from that chart on.  An unprescribed branch point in the
> INTERIOR of the path errors regardless of the flag (a crossing cannot be
> half-performed honestly).  Chaining from a MultivaluedFail result is E25.

### 25. [major] Missing unit tests for declared error contracts — every "loud error" without a test is a future silent regression

Problem.  Spec-by-spec, the following declared error/invariant IDs have NO
test exercising them (test lists §8 of each spec, cross-checked): Config —
E7 (delta sign/nonlinear δ), E13 (LogFile; see defect 27).  EpsSeries —
ESMap (any behavior incl. its eps-introduction error), the ESSameQ band
behavior (after defect 2).  SectorSeries — `::geomambiguous`,
`::inexactdenominator`, `::complexpoint`, invariant I-10 (ErrorEstimate
monotonicity).  Solve — E2 (spectral mismatch), E4 (work-window overflow),
E5 (ladder inconsistency/I-5 failure path), E6 (empty delivered window),
E8 (non-finite coefficient).  Transport — E2 (line on singularity), E8
(prescription conflict/missing — the PENTAGON-CLASS error, the single most
campaign-relevant error in the module), E9, E11 (accuracy miss), E12
(window underflow), E14, invariant I11 (chart-level ODE residual).
Integrate — E7 (duplicate tags), E8 (estimate abort), E9a, E9b, E10
(malformed tag), E11, I9.  API — E12, E13, E14, E16, E17, E21, E22, E23,
E27, E28, E31, E33, E34; W2 and W4 have no positive-fire test.

Amendment (exact text).  ADD to each spec's §8 header line:

> COVERAGE RULE (M0 review): every error ID and warning ID of section 5
> has at least one unit test that triggers it and asserts its required
> payload fields, OR a one-line waiver in the test file naming the ID and
> why it is untestable in isolation (e.g. "unreachable by construction,
> guarded by assert X").  The IDs currently missing tests are enumerated
> in Docs/specs/REVIEW-minimalism.md defect 25; the implementation agent
> closes the list before the module's milestone gate.

Transport E8 specifically: ADD test "T24 `test_prescription_conflict_loud`:
two odd-multiplicity prescriptions deriving opposite sigmas at one chart
on a crossing path -> E8 with the per-factor table in the payload; same
chart as FINAL chart with AbortOnAnalyticContinuationFail -> False ->
honest-partial result with MultivaluedFail metadata (defect 24 semantics)."

---

## MINORS

### 26. [minor] Config.md §3.2 + §8 — Crosscheck keys are surface with zero specified consumers

Problem.  "CrosscheckLevel"/"CrosscheckFlags" are kept keys, but no module
spec defines a single optional check, a registry name, or a gated code
path; Config Q2's load-time registration mechanism is itself unresolved.
Grep: one in-repo consumer total, passing `{}` (Tests/
test_singular_recurrence.m:278 — an old-core test that dies in M6).  Kept
keys with no behavior = silent no-ops.

Amendment.  In Config.md §3.2 MOVE both keys to the dropped-keys table:
"| 'CrosscheckLevel' / 'CrosscheckFlags' | DiffExp2's invariants are
always-on (RewritePlan A2) and not configurable; no optional-check
registry exists in v1.  Re-add (fixed name list, no runtime registration)
only when a module specifies a concrete expensive check.  Ledger waiver
required. | State.m:110, 196-207 |"  DELETE C20 and Q2; retarget C20 to
the dropped-key error.  Alternative (if the orchestrator wants the hook):
keep the keys but then Solve.md/Transport.md MUST each name one registered
check in their §2 — currently none does, and that, not the keys, is the
defect.

### 27. [minor] Config.md §3.2 "LogFile" — no consumer anywhere in the repo

Problem.  Grep over Reference/Examples, Tests, Scripts, FeynmanTrick: zero
hits.  The key buys an error path (E13), schema surface, and the only file
I/O in Config (its own §7 carve-out exists solely for this key).

Amendment.  MOVE "LogFile" to the dropped-keys table ("no in-repo
consumer; session logging is the shell's job; E9 names the removal"),
DELETE E13 and the §7 I/O carve-out, and DELETE §9 cut option (3) (moot).

### 28. [minor] Integrate.md §2.2.2 — the t-tail estimator divides by zero when the top kept power is the pole power

Problem.  The tail estimate `|c[k,nmax,comp]| · T^(a+nmax+1) / |a+nmax+1|`
is undefined at `a + nmax + 1 = 0` (top kept term IS the t^-1 cell); the
parenthetical "for the pole-cell sector use the same bound at n = nmax"
does not fix the division.

Amendment.  REPLACE the estimator sentence with: "tail estimate per eps
order: `|c[k,nmax,comp]| * T^(a+nmax+1)/|a+nmax+1|` when
`a+nmax+1 != 0`, and `|c[k,nmax,comp]| * |Log T|` when `a+nmax+1 == 0`
(the next term's antiderivative magnitude in the log cell)."

### 29. [minor] SectorSeries.md §2.4/§3.4 — EvaluateLocalSolution returns "TailEstimates" but only re-expansion (§2.6) defines how they are computed

Problem.  §3.4 lists `"TailEstimates" -> per-eps-order list or Missing[]`
on the evaluation result; §2.4.1 never says when it is computed or what
formula evaluation uses (the 2.6 formula is specific to re-expansion's new
chart).  An implementer must guess — guessed error metadata is worse than
none.

Amendment.  ADD to §2.4.1 after the value-window paragraph: "TailEstimates
on EVALUATION results: per eps order, `Max over comp of
|c[k, ntop, comp] * tval^ntop| * q/(1-q)` with `q = |tval|/Radius` and
`ntop` the last used t-column — the geometric tail at the actual
evaluation point.  Computed always (cheap); Missing[] never occurs on the
evaluation path (the §3.4 Missing[] option applies only to objects that
bypassed evaluation).  Transport's probe (its 2.10) remains the
authoritative error feed; TailEstimates is advisory metadata."

### 30. [minor] API.md §2.1(c)/§5 E1 vs RewritePlan M4 — verify that no M4 parity line crosses a singular chart with a LegacySlices directory, or the gate is unrunnable

Problem.  E1 fires for TransportTo "when the path contains or targets a
singular chart" in LegacySlices mode.  Classic parity lines (M4: bubble/
sunrise/2F1/banana) cross on-path singular points (banana crosses its
thresholds; that is what the crossing machinery is for).  If any M4 oracle
directory is slices-only, M4 cannot run under DiffExp2 as specified.  The
classic examples ship closed-form matrices (mode (b)), so this is likely
fine — but nothing in the specs pins it, and the M0 oracle-generation task
could innocently dump slices.

Amendment.  ADD to API.md §2.1(c): "M0/M4 NOTE: every parity-oracle matrix
directory vendored for M4 (RewritePlan M0 KERNEL item) MUST be Exact-mode
(closed-form d<var>_d.m or full export) — the M4 lines cross singular
charts, which LegacySlices mode refuses by design.  The M0 vendoring task
asserts this at generation time."

### 31. [minor] Tolerances.md §8 T3 — the unit suite of the bottom module silently depends on Config

Problem.  T3 installs state "via a minimal Config round-trip", making
test_tolerances.m unrunnable before Config.m exists, inverting the build
order for tests (M1 builds Tolerances+Config together so it happens to
work, but the file-level dependency contradicts §7's "depends on NOTHING"
and will bite anyone bisecting M1).

Amendment.  REPLACE T3 with: "T3 `test_match_tol_sync_default`: install
directly via `InstallToleranceState` with the record Config WOULD produce
at wp=500 defaults (the §3.4 Config.md table is the source of the
expected values); assert `Tol["MatchTol"] === Tol["ChopFloor"] ===
10^-250`.  The Config-side half of the sync behavior is tested in
test_config.m C6, not here."

### 32. [minor] Solve.md §5 E9 — the "source in d/dt form detected" clause is mathematically incoherent

Problem.  E9's middle clause ("a source sector with a <= -1 integer and
b = 0 colliding with nothing is flagged when its particular would
integrate t^{-1}: that specific divergence is a loud error here") describes
a divergence that does not exist in formal Frobenius solving: a theta-form
source with tag a = -1, b = 0 at a regular chart solves through CASE T
(offset -1+n vs root 0; no division by zero at any n... the n = 1 step has
offset 0 -> CASE R log bump, which the spec already handles).  There is no
"divergence" for Solve to detect; divergences are Integrate.m's domain.
An implementer trying to honor this clause will invent a spurious error
path (or worse, a heuristic d/dt-form detector).

Amendment.  REPLACE E9 with: "E9 MALFORMED REQUEST/SOURCE: req windows
inverted (Min > CompleteMax); source sector tags non-exact (inexact a/b
heads, negative or non-integer p); source TWindow/EpsWindow inconsistent
with its Coeffs dimensions.  The theta-form normalization contract (3.3)
is NOT mechanically detectable and is not guessed at: a d/dt-form source
produces wrong (shifted-tag) results that the always-on ODE residual check
(E7, which includes the source term) catches — that is the designed
backstop, and a code comment at the SourceData validator must say so."

### 33. [minor] Indicial.md §7 vs defect-5 amendment — Indicial's Config dependency exists only for facilities Config doesn't provide

Problem.  Indicial §7: "May call Config.m: ONLY for the shared loud-error
primitive / chart-name formatting and verbosity-gated info logging (if
Config.md provides them; OQ-2)."  After defect 5's resolution the error
primitive is Tolerances.m's and the print gates are Config.m's; Indicial's
hedged sentence should bind to the real names so the module's dependency
list is checkable.

Amendment.  REPLACE the Config bullet of Indicial §7 with: "May call:
``Tolerances`LoudError`` (sole Tolerances use — the I-8 'no numeric
tolerances' invariant is about thresholds, not the error helper) and
``Config`PrintInfo`` for verbosity-gated logging.  No configuration key
changes Indicial's mathematics."  Update I-8's grep test T-21 to allow the
LoudError symbol.

### 34. [minor] Transport.md §2.10 vs Tolerances.md Q3 — the promised exact-confirmation rule for snap-prefiltered dedup never landed in Transport

Problem.  Tolerances Q3: "numeric dedup at SnapTol is allowed only as a
PRE-filter, with exact RootReduce-level confirmation required before two
roots are merged; Transport.m spec must repeat this."  Transport.md §2.1
specifies exact-only dedup and never mentions the pre-filter rule — fine
in itself (stricter), but the dangling cross-reference invites a future
implementer to add the pre-filter WITHOUT the confirmation clause.

Amendment.  ADD one sentence to Transport.md §2.1's dedup bullet: "A
numeric SnapTol pre-filter MAY be used to find merge CANDIDATES on large
root sets, but two roots are merged only after exact RootReduce
confirmation (Tolerances.md Q3); an unconfirmed near-pair is kept
distinct."  Mark Tolerances Q3 resolved.

### 35. [minor] Indicial.md §9 vs RewritePlan 3.2 — budget self-declared over with the most valuable cut (C-4) requiring a Solve.md change that Solve.md doesn't anticipate

Problem.  Indicial lands at ~330 vs ~300 and its C-4 ("move chain
construction to Solve.m ... adds ~40 there") would push Solve over ITS
budget (700, fully allocated) — the cut is an inter-module budget
transfer presented as a local option.  After defect 15's amendment (Solve
builds V/VInv from Chains), C-4 becomes outright incompatible (Solve needs
the chains).

Amendment.  In Indicial.md §9 DELETE cut C-4 and ADD: "C-4 is withdrawn
(REVIEW-minimalism defect 35): Solve.m consumes Chains verbatim (defect-15
resolution), so chain construction cannot leave Indicial.  Remaining
headroom comes from C-2/C-3 (pre-authorized, defect 4) and, on M2
evidence only, C-1."

### 36. [minor] EpsSeries.md §2 (`ESCoefficientList` ERR-DROP-BELOW) vs Solve.md §7 ("polar-part extraction")

Problem.  Solve lists "polar-part extraction" among the EpsSeries services
it consumes; EpsSeries exports no such function, and the natural spelling
(`ESCoefficientList[s, ESMinPower[s], -1]`) is legal but nowhere blessed —
an implementer may instead slice with a positive k1 and trip
ERR-DROP-BELOW in a context where dropping is intended (compensation
construction needs polar AND regular parts separately).

Amendment.  ADD to EpsSeries.md §2 under ESCoefficientList: "POLAR-PART
IDIOM (blessed): the polar part of `s` is `ESCoefficientList[s,
ESMinPower[s], -1]` (empty list when Min >= 0) and the regular complement
is `ESCoefficientList[ESShift[ESTrim[...]]...]`-free: simply
`ESCoefficientList[s, 0, ESCompleteMax[s]]` AFTER the caller has handled
the polar part — this pair never triggers ERR-DROP-BELOW and is the form
Solve.m's §3.7 compensation uses.  No new export is added."

---

## Tally

Blockers: 4 (defects 1-4).  Majors: 21 (defects 5-25).  Minors: 11
(defects 26-36).  Recurring failure modes worth the orchestrator's
attention when merging amendments: (a) every module hedged shared
decisions behind "open question — align later" and the aggregate of those
hedges is itself the defect (error primitive, eps-Laurent shape, spectral
data contract, eps symbol ownership — defects 5, 15, 16, 23); (b) three
specs ship machinery a fourth deletes (defect 8) — the kept-key table of
Config.md must be regenerated AFTER Transport/Integrate amendments land,
not before; (c) the budget has no slack (defect 4): treat every "+~N
lines" in the amendments above as already spent against the pre-authorized
cuts.  LessonsLedger.md and ExportDisposition.md did not exist at review
time; defects 8, 18, 21, 26, 27 each create a ledger-waiver obligation
that must be checked into the ledger when it lands.
