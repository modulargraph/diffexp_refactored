# Module spec: DiffExp2/Tolerances.m

Status: M0 deliverable per Docs/RewritePlan.md section 3.2 / section 6 (M0 task 4-11).
Context name: ``DiffExp2`Tolerances` ``.  Line budget: ~100 (see section 9).
Implementable from: this file + Docs/RewritePlan.md + the old-code citations below.
Module order (acyclic, hard rule): Tolerances < Config < EpsSeries < SectorSeries <
Indicial < Solve < Transport/Integrate < API.  Tolerances is the BOTTOM module:
it depends on NOTHING (not even Config).  Every threshold is a pure function of
explicitly passed precision integers; Config.m computes/validates those integers
and installs a shared tolerance state here (one-way call, Config -> Tolerances).

---

## 1. PURPOSE

Tolerances.m is the single home for every numeric threshold in DiffExp2, replacing
the old library's disease-D3 scatter of ad-hoc constants (absolute `10^-40` in
FeynmanTrick/DiffExpIntegration.m:1177-1178, absolute `10^-10`
RationalizationTolerance in DiffExp/State.m:123, `10^-2`/`10^-12` cluster and
moment tolerances in DiffExp/RegularizedIntegration.m:862,941, hardwired option
defaults `"Tolerance" -> 10^-10` in DiffExp/LocalSeries.m:480,492, chop constants
in DiffExp/State.m:109,119,149,212).  Each threshold gets a NAME with stated
semantics (chopFloor, matchTol, snapTol, inputSnapTol, rankTol, geomGuardTol,
laurentLeadTol, residTol, chopReserve), is DERIVED from WorkingPrecision or
from the input's own accuracy (never a free literal at a call site), and zero
tests that classify structure (leading-coefficient trimming, rank decisions) are
RELATIVE with a loud-error ambiguity band — except at the universal 10^-24
floor, where the campaign calibration makes the test binary (DEC-2) — so a
borderline value can never be silently classified either way.  The module also
centralizes the SURVIVING operational safety constants of the old library (the
`I*` constants of DiffExp/State.m:209-224 minus the dropped adaptive-search
set, DEC-11) so that no tuning number lives anywhere else in DiffExp2.  As the
BOTTOM module, Tolerances.m additionally exports the ONE library-wide
loud-error primitive `DE2Error` (DEC-1) that every other module raises through.

## 2. PUBLIC SYMBOLS

All exported with `::usage` before `Begin["`Private`"]` (the unexported-symbol
silent no-op trap — see memory/wolfram-package-context-traps.md and the
NormalizeLogPower incident, commit ebc4724 — must be impossible: an in-repo test
asserts every symbol listed here evaluates from a foreign context, see section 8).

Derivation functions (pure; `wp`, `chopDigits`, `matchDigits` are positive
machine Integers; all return exact rationals `10^-k` unless stated):

- `ChopDigits[wp_Integer] -> Integer`
  Default chop-digit count: `Floor[wp/2]`.  (Old default pair: ChopPrecision 250
  at WorkingPrecision 500, DiffExp/State.m:109,135 — same ratio.  The FT campaign
  pinned `precision - 50` at DiffExpIntegration.m:349; Config may override, this
  function only supplies the default.)
- `ChopFloor[chopDigits_Integer] -> 10^-chopDigits`
  ABSOLUTE noise floor for series/array coefficients.  Replaces
  `PChop = Chop[#, 10^-FEC[ChopPrecision]]&` (DiffExp/Utilities.m:103) and every
  PChop call site that survives into DiffExp2.
- `MatchTol[matchDigits_Integer] -> 10^-matchDigits`
  ABSOLUTE zero test inside boundary-matching / linear-solve pivoting.  Replaces
  `LSPChop` (DiffExp/Utilities.m:104), the `LinearSolve ZeroTest` lambdas
  (DiffExp/Frobenius.m:61, DiffExp/Transport.m:341,
  DiffExp/IntegrationStrategies/Helpers.m:82, VOP.m:118,246) and
  `LinearSolveChopPrecisionVal` (DiffExp/State.m:148).  `matchDigits` defaults to
  `chopDigits` (the auto-sync lesson, DiffExp/DiffExp.m:126-130).
- `SnapTol[wp_Integer] -> 10^(-Floor[wp/2])`
  Tolerance for snapping a computed numeric value to a nearby EXACT target
  (chart-map endpoint values to 0/1, dedup of numeric singularity positions).
  Replaces the absolute `RationalizationTolerance = 10^-10` (DiffExp/State.m:123)
  in its two SURVIVING roles (FeynmanTrick/DiffExpIntegration.m:151-178
  snapMainExpression, 180-202 snapValuesFromFactors).  Precedent for the wp/2
  scaling is in-tree: `Chop[Uval - 1, 10^(-precision/2)]`
  (FeynmanTrick/BoundaryConditions.m:191,195).  SnapTol applies to COMPUTED
  full-precision quantities ONLY (chart-map endpoint values, factor-root dedup
  pre-filter); it plays NO role in matching/rank solves (REVIEW-minimalism 10)
  and EXTERNAL numeric inputs are snapped/guarded by InputSnapTol and
  `$NearSingularityGuardDecades` below (DEC-14).
- `InputSnapTol[v_] -> 10^(-Floor[3*Min[Precision[v], Accuracy[v]]/4])`
  (argument: an INEXACT numeric input; exact input uses exact membership only.)
  Input-scaled snap tolerance for EXTERNAL targets/bounds (DEC-14;
  REVIEW-math D19): external numeric inputs carry input-resolution error;
  snapping them at the internal SnapTol (wp/2 decades) is a dead test — the
  lesson behind the old absolute 10^-10 (State.m:123), made input-aware instead
  of absolute.  Endpoint classification rule (consumed by Transport.md E14/W5
  and API.md W5): an inexact target `to` with `|to - root| < InputSnapTol[to]`
  for some exact singularity root is SNAPPED to the root with warning W5 (a
  1e-16-off machine-precision target and a 30-digit 0.25 both snap — the
  input's own accuracy, not wp, sets the snap scale); an inexact `to` with
  `InputSnapTol[to] <= |to - root| < 10^-$NearSingularityGuardDecades` is a
  LOUD ERROR ("target is suspiciously near singularity <root>; pass the exact
  value, or an exact offset, to confirm intent").  Exact inputs are never
  snapped and never guarded (exact arithmetic decides).
- `$NearSingularityGuardDecades = 6` (new; REVIEW-minimalism 11): outer edge
  of the loud-guard zone in the InputSnapTol rule above.
- `RankTol[wp_Integer] -> 10^(-Floor[wp/4])`
  RELATIVE rank/nullspace threshold: a singular value (or pivot) counts as zero
  iff `sv < RankTol[wp] * svMax`.  Rank/pivot cuts must sit a band above the
  snap floor (DEC-14; REVIEW-math D15): the gray zone
  `[SnapTol*scale, RankTol*scale]` (width wp/4 decades) is the loud-error
  region of Transport.md E7, and the NumericallyZeroQ ambiguity band at RankTol
  IS that gray zone (REVIEW-minimalism 10) — not a separate mechanism.
  Replaces the absolute
  `NullSpace[..., Tolerance -> 10^-LinearSolveChopPrecisionVal]` sites
  (DiffExp/Wronskian.m:63, IntegrationStrategies/Helpers.m:32,84, VOP.m:120,248)
  and the fitter's `Min[svals] >= 10^(-activeNumericPrecision[]/2) * ...` test
  (DiffExp/RegularizedIntegration.m:1047).  Rank decisions create/destroy
  structure (disease D1) and therefore MUST be relative and MUST go through
  `NumericallyZeroQ` (ambiguity band) — never a raw comparison.
- `GeomGuardTol[wp_Integer] -> 10^(-Floor[wp/2])`
  Radius-vs-pole-modulus comparison guard consumed by SectorSeries
  `::geomambiguous` and Transport geometry asserts: comparisons whose
  difference magnitude falls below `GeomGuardTol[wp] * Max[|t_i|, Radius]` are
  LOUD errors, not decisions (REVIEW-minimalism 17; supersedes REVIEW-math
  D26's RankTol-based guard — the named export keeps the wp/2 width that D26
  intended after RankTol moved to wp/4, honoring F4's named-export rule).
- `ChopReserve[wp_Integer, chopDigits_Integer] -> wp - chopDigits`
  The digit reserve between working precision and the chop floor; Transport.md
  E3 requires `DigitsNeeded + ChopReserve <= wp` (old implicit reserve
  WP - ChopPrecision = 250, State.m:109,135).
- `LaurentLeadTol[chopDigits_Integer] -> Max[10^(-Floor[chopDigits/2]), 10^-24]`
  RELATIVE leading-coefficient zero test for eps-Laurent window trimming and
  for object-level cancellation gates (DEC-2).  This is the replacement
  demanded by RewritePlan 3.2 for the FT `zeroCoeffQ` ABSOLUTE `10^-40` test
  (FeynmanTrick/DiffExpIntegration.m:1177-1178, consumed by LaurentTrim at
  1180-1190).  The 10^-24 floor is the campaign-verified calibration
  (DiffExp/IntegrationStrategies/Recurrence.m:613-618;
  Docs/FeynmanTrickBoxFamilyStatus.md "the 1e-24 floor is load-bearing"):
  incomplete cancellations at apparent singularities leave residues up to
  ~1e-29 RELATIVE at WP 300 which MUST classify zero, while genuine leading
  content is O(1) relative or exactly zero.  The floor is a HARD CONSTANT
  (not configurable, not separately exported — DEC-2) and is a RELATIVE
  quantity, so it does not reproduce forbidden fallback F1 (which bans
  absolute tests on raw coefficients).  Semantics: leading coefficient `c` of
  a Laurent array with window scale `s = Max[Abs /@ allWindowCoefficients]` is
  zero iff `NumericallyZeroQ[c, s, LaurentLeadTol[chopDigits], context]`;
  `s == 0` => only exact zeros trim.  When the Max picks the floor
  (chopDigits >= 48 — every campaign setting) the classification is BINARY:
  the ambiguity band does NOT apply at the floor (DEC-2; a band abort there
  would fire on the correct ~1e-29 residue population).  For the Integrate.m
  cancellation gate the scale is the max |coefficient| of the merged
  combination at that eps order (DEC-4).
- `ResidTol[wp_Integer] -> 10^(-Floor[wp/10])`
  RELATIVE threshold for the always-on ODE-residual spot-check (RewritePlan 3.1
  invariants) and for construction-time cross-checks.  Replaces the old absolute
  30-digit `CrosscheckChopPrecision`/CPChop (DiffExp/State.m:149,
  DiffExp/Utilities.m:105): correct solutions give residuals at the noise floor
  (`~10^-chopDigits` relative), wrong ones give O(1); wp/10 (default 50 decades)
  sits far from both.
- `NumericMagnitude[c_?NumericQ, digits_Integer] -> nonnegative numeric scalar`
  Stable mathematical modulus for numeric coefficient and residual decisions.
  For a complex value, inspect the real and imaginary stored centers
  independently and evaluate a scaled Euclidean norm.  This avoids Wolfram's
  `Abs` precision poisoning when (for example) a resolved real component is
  paired with a centered-zero imaginary component carrying poor accuracy,
  while retaining the true modulus rather than the smaller infinity norm.
  The helper reports the central magnitude only; residual proof sites account
  for the original value's precision uncertainty separately.  Geometry keeps
  using exact/numeric `Abs` directly.
- `NumericMagnitudeBounds[c_?NumericQ, digits_Integer] -> {lower, upper}`
  Componentwise uncertainty enclosure for the same modulus.  Exact inputs
  remain exact; inexact real/imaginary centers receive independent
  `10^-Accuracy` radii evaluated at arbitrary precision.  Residual proofs and
  error bounds use `upper`; rank classification uses both endpoints.

The shared classification predicate:

- `NumericallyZeroQ[c_, scale_, tol_, context_String,
   bandDecades_:$AmbiguityBandDecades, binaryFloor_:False]
   -> True | False | LOUD ERROR`
  THE one zero-classification predicate for all relative tests in DiffExp2
  (DEC-3): exactly one library-wide; EpsSeries' `ESCoeffZeroQ` is a thin
  wrapper over it that adds window context to the error payload — same
  semantics, never a second predicate.  Contract: (1) for exact data only,
  `PossibleZeroQ[c]` exactly True -> True; inexact approximate zeros proceed
  through the uncertainty proof.  (2) `c` not NumericQ -> False (symbolic
  content is never "numerically zero"; callers handle symbols structurally).
  (3) `scale == 0` (no reference magnitude) -> False unless branch (1) fired.
  (4) Otherwise derive componentwise lower/upper magnitude bounds `(lo, hi)`
  from the stored centers and an arbitrary-precision
  `10^-Accuracy[component]` uncertainty (so it cannot machine-underflow).  A
  poor centered-zero imaginary component
  therefore cannot erase a separately resolved real lower bound.
  FLOOR EXEMPTION (DEC-2/DEC-3): only a Laurent-semantic caller for which
  `LaurentLeadTol` actually selected the hard `10^-24` floor passes
  `binaryFloor -> True`; then `hi <= tol*scale -> True`, else False.  The mode
  is explicit because `RankTol[96]` also equals `10^-24` numerically and must
  remain ternary.  (5) Otherwise ternary: `hi < tol*scale/10^bandDecades ->
  True`; `lo > tol*scale*10^bandDecades -> False`; ELSE LOUD ERROR (section 5, E5)
  quoting `context` verbatim — `context` MUST already contain the chart /
  sector tag (a,b,p) / eps-order / t-order identifiers the caller has.
  `bandDecades` defaults to the exported constant `$AmbiguityBandDecades = 4`
  (recalibrated from 10 by DEC-2/DEC-3: the two documented populations are
  ~1e-29-relative cancellation residues vs O(1)-relative signal; every
  in-scope zero/nonzero population pair is separated by far more than
  2*4 decades, and the narrower band aborts on strictly less).
  This predicate is what makes "IntegerQ-on-floats" and threshold-straddling
  misclassification (disease D1; the LaurentTrim MinPower-shift hazard of
  legacy-review finding 13a, Docs/reviews/rewrite_plan_review_3lens.json) impossible.

The shared loud-error primitive (DEC-1; REVIEW-math D17, REVIEW-minimalism 5 —
exported from the bottom of the dependency order so every module can raise
through it without inverting the module order):

- `DE2Error[id_String, payload_Association] -> (never returns)`
  THE library-wide error primitive.  Behavior: prints one line
  `"DiffExp2 error <id>: <payload as key=value summary>"`, then
  `Throw[Failure["DiffExp2", Join[<|"ID" -> id|>, payload]], "DiffExp2Error"]`.
  `payload` always carries `"Module"` (and `"ID"` after the join) plus
  whatever chart / sector / eps-order / t-order context the caller's
  `$ESErrorContext`-style scoping provides.  The catch sits at every API.m
  entry point — API.m's entry points are the ONLY `Catch[..., "DiffExp2Error"]`
  sites and convert the Failure to a user-facing abort.  No module defines its
  own Abort/Message idiom.  DE2Error has no state dependency (it must work
  before any tolerance install).  (~12 lines, funded by the deleted
  adaptive-search constants — REVIEW-minimalism 5/8.)

Installed-state accessors (state is write-once-per-configuration, installed by
Config.m ONLY):

- `InstallToleranceState[assoc_Association] -> Null`
  `assoc` must contain EXACTLY the keys
  `{"WorkingPrecision","ChopDigits","MatchDigits","ChopFloor","MatchTol",
    "SnapTol","RankTol","LaurentLeadTol","ResidTol"}`
  (String keys; Integer for the first three, exact positive rationals for the
  rest).  Missing key, extra key, or type mismatch -> LOUD ERROR.  Re-install
  (configuration update) is allowed and replaces the whole record atomically.
- `Tol[name_String] -> value`
  Validated read of the installed state.  Unknown name -> LOUD ERROR listing the
  nine valid names.  Called before any install -> LOUD ERROR ("no configuration
  loaded; call DiffExp2 LoadConfiguration/LoadSystem first").  NEVER returns
  Missing[...] or a default.
- `ToleranceStateInstalledQ[] -> True|False` (for API.m error paths and tests).

Exported operational constants (values are the old library's, with citations;
these are data, not functions):

- `$SafetyDigits = 2` (ISafetyDigits, DiffExp/State.m:215; digit budget
  `DigitsNeeded = AccuracyGoal + Ceiling[Log10[#segments]] + $SafetyDigits`,
  DiffExp/Transport.m:527,759 — consumed by Transport.m).
- `$InputPrecisionFactor = 2` (inputs raised to 2x WorkingPrecision:
  DiffExp/Transport.m:543,560; DiffExp/Mobius.m:48,66 — consumed by
  Transport.m/API.m; ledger seed "2x headroom survives chart-map composition").
- `$MinPrecisionFloor -> wp` semantics: all evaluation happens inside
  `Block[{$MinPrecision = wp}]` (DiffExp/Pade.m:34,70,80) with
  `$MaxExtraPrecision = $MaxExtraPrecisionValue = 1000` (DiffExp/Mobius.m:39).
  Tolerances exports the two constants; the Blocks live at the consumer.
- `$MinExpansionOrder = 10` (DiffExp/State.m:223; floor consumed by Config's
  ExpansionOrder validator — the ONE survivor of the State.m:216-224 set).
  The other six adaptive-search constants ($SafetyExpansionSubtract,
  $ExpansionOrdersAveraging, $ExpansionOrderIncrease, $ExpansionOrderDecrease,
  $ExpansionOrderIncreaseValidate, $DigitsSurplusDecrease) are DELETED, not
  ported and not folded: their only consumer, the AccuracyGoalValidate
  adaptive expansion-order search (DiffExp/Transport.m:776-841,1191-1213) and
  the LineSegmentation.m:109-113 truncation read-off, is not ported (DEC-11;
  REVIEW-minimalism 8 — LessonsLedger waiver "expansion-order adaptive
  search"; the -10 lines fund DE2Error per REVIEW-minimalism 4/5).
- `EvalErrorSeriesDecrease[couplingDepth_Integer] := Ceiling[0.7*couplingDepth] + 2`
  (ICurrEvalErrorSeriesDecrease, DiffExp/State.m:222; the two-point error probe
  order reduction, DiffExp/Transport.m:913).
- `$AmbiguityBandDecades = 4` (new; see NumericallyZeroQ — recalibrated from
  10 by DEC-2/DEC-3, and exempt at the 10^-24 floor).

## 3. DATA CONTRACTS

- Tolerance state record (produced by Config.m, consumed here, read by every
  higher module through `Tol`):
  `<|"WorkingPrecision"->Integer, "ChopDigits"->Integer, "MatchDigits"->Integer,
    "ChopFloor"->Rational, "MatchTol"->Rational, "SnapTol"->Rational,
    "RankTol"->Rational, "LaurentLeadTol"->Rational, "ResidTol"->Rational|>`
  All tolerance values exact `10^-k` rationals (never machine reals — they are
  compared against arbitrary-precision numbers).
- Consumers and which name they may use (binding; anything else is a spec
  violation found in review — table per REVIEW-minimalism 17 / REVIEW-math
  D16; DE2Error is an error helper, not a tolerance, and is implicitly
  available to every module):
  - EpsSeries.m: `LaurentLeadTol` (window/lead classification via
    NumericallyZeroQ only — ESCoeffZeroQ is the thin wrapper, DEC-3),
    `MatchTol` (ESSameQ).  NO chopping: EpsSeries never alters stored
    coefficients (its I-6).
  - SectorSeries.m: `ChopFloor`, `ResidTol` (evaluate/re-expand checks),
    `GeomGuardTol` (geometry ambiguity guard, `::geomambiguous`).
  - Indicial.m: NOTHING — indicial data is EXACT (RewritePlan I1; its §7
    "NOT USED" statement governs, and I-8/E1 reject numeric entries
    outright).  Exponent rationalization tolerances are DELETED (section 6).
  - Solve.m: `MatchTol`, `RankTol`, `LaurentLeadTol` (assembly-time
    log-member trim, REVIEW-minimalism 1), `ResidTol`.
  - Transport.m: `MatchTol`, `RankTol` (via NumericallyZeroQ — the band IS
    the E7 gray zone), `ResidTol`, `SnapTol` (computed-value snapping only),
    `InputSnapTol` + `$NearSingularityGuardDecades` (E14/W5 endpoint rule),
    `ChopReserve` (E3 digit budget), `GeomGuardTol` (geometry asserts),
    `$SafetyDigits`, `$InputPrecisionFactor`, `$MinExpansionOrder`,
    `EvalErrorSeriesDecrease`.
  - Integrate.m: `LaurentLeadTol` (cancellation gate at the DEC-4
    merged-combination per-eps-order scale, + result trim), `SnapTol`
    (tiling/endpoint snapping), `MatchTol` (additivity spot-check I6).
  - API.m / FT layer: `SnapTol` (chart-map endpoint snapping, singularity
    dedup pre-filter), `InputSnapTol` (W5 user-target snapping), `ChopFloor`.
- LocalSolution / Sector / EpsWindow (RewritePlan 3.1, verbatim definitions are
  normative there): Tolerances.m never constructs or stores these; it only
  classifies their numeric coefficients.  In particular EpsWindow
  `<|"Min"->kmin,"CompleteMax"->kmax|>` movement is decided in EpsSeries.m, but
  ONLY through `NumericallyZeroQ` with `LaurentLeadTol` — no other predicate may
  move a window boundary.

## 4. INVARIANTS (always on, cheap)

- Every derivation function: argument positive Integer, result positive exact
  rational < 1; violated -> loud error (not silent coercion).
- `chopDigits < wp` and `matchDigits <= chopDigits` checked at
  `InstallToleranceState` (the ChopPrecision < WorkingPrecision invariant,
  DiffExp/DiffExp.m:121-125, hoisted here so it can never be bypassed).
- Installed record: exactly the nine keys, types as in section 3 (asserted at
  install, so `Tol` reads need no per-call checks).
- `NumericallyZeroQ` never returns for ambiguous input (it throws); asserted by
  unit T9.
- All exported symbols carry `::usage` (export-completeness test T13).

## 5. ERROR CONTRACT (no silent fallbacks anywhere)

All errors are raised via `DE2Error` (section 2; DEC-1); this module defines
no other error mechanism, and no other spec owns the form (the old "API.m spec
owns its exact form" clause is deleted — REVIEW-math D17 / REVIEW-minimalism 5;
API.m's entry points only own the `Catch[..., "DiffExp2Error"]`).  Each
condition must carry the strings shown.  Enumerated conditions:

E1 `Tol` before install: "Tolerances: no tolerance state installed (key <name>
   requested). Load a configuration first."  MUST NOT return a default.
E2 `Tol` unknown name: "Tolerances: unknown tolerance name <name>; valid names:
   <list>".  This is the same key-discipline as Config.m: a read of a
   nonexistent entry is an ERROR, never Missing[].
E3 `InstallToleranceState` schema violation: names the missing/extra keys and
   offending types.
E4 `InstallToleranceState` with `chopDigits >= wp`: message text mirrors the old
   one — "The value of ChopPrecision should be smaller than the value of
   WorkingPrecision" (DiffExp/DiffExp.m:123) — plus the two numbers.
E5 `NumericallyZeroQ` ambiguity band: "Tolerances: cannot classify coefficient
   <N[c,6]> against scale <N[scale,6]> at tolerance <tol> (within
   10^±$AmbiguityBandDecades). Context: <context>."  The caller-supplied
   `context` string MUST contain chart center, sector tag (a,b,p) and eps/t
   order when the call originates from window or matching logic; EpsSeries.m,
   Solve.m and Transport.m specs repeat this obligation on their side.
E6 Derivation function with non-Integer / non-positive argument: names the
   function and the argument.

Forbidden fallbacks (each was a real site; reproducing any of them is a defect):

F1 Absolute zero test on a Laurent leading coefficient
   (`Abs[c] < 10^-40`, FeynmanTrick/DiffExpIntegration.m:1177-1178).  Forbidden:
   only `NumericallyZeroQ` with `LaurentLeadTol` and a window scale.
F2 Hard-coded tolerance default when a config read fails
   (`Quiet[Check[FEC[RationalizationTolerance], 10^-40]]`,
   FeynmanTrick/DiffExpIntegration.m:154-155,182-183).  Forbidden: `Tol` either
   answers or aborts; there is nothing to Check.
F3 Returning a permissive result when the state is missing (old `FEC[key]` ->
   `Missing["KeyAbsent",...]` flowing into `TrueQ[x < Missing[...]] == False`,
   live today at FeynmanTrick/DiffExpIntegration.m:98,109,127,924).  Forbidden
   by E1/E2.
F4 Silently widening/narrowing a tolerance at a call site ("Max[tol^(1/3),
   10^-12]", DiffExp/RegularizedIntegration.m:941; "clusterTol = 10^-2",
   RegularizedIntegration.m:862).  Forbidden: call sites may not transform
   tolerance values; if a different semantics is needed it gets a new NAMED
   export here.
F5 Rationalizing/snapping exponents or eigenvalues numerically
   (`Rationalize[eigenvalues, 10^(-ChopPrecision/2)]`,
   DiffExp/IntegrationStrategies/Recurrence.m:136; `Rationalize[#,
   10^-ChopPrecisionVal]`, DiffExp/Frobenius.m:36; NormalizeXPower / b-snapping,
   DiffExp/SingularityDecomposition.m:52-64,336-362; NormalizeLogPower,
   DiffExp/SeriesOps.m:95-106).  Forbidden everywhere in DiffExp2 core: sector
   tags are exact algebra (RewritePlan I1/I2); there is deliberately NO
   "exponentTol" export, so any attempt to reintroduce one is visible in review.
F6 Chopping with a literal (`Chop[x, 10^-20]`-style debug chops are fine in
   PRINT statements only, as in RegularizedIntegration.m:1016; any literal
   tolerance in non-print logic is forbidden).

## 6. ABSORBED OLD CODE (and the numerical lessons preserved)

- DiffExp/Utilities.m:103-105 (PChop/LSPChop/CPChop) -> ChopFloor / MatchTol /
  ResidTol.  Lesson: three DISTINCT semantics existed but were all absolute
  chops; keep the distinction, fix the relativity where classification happens.
- DiffExp/State.m:109,119,135 + DiffExp/DiffExp.m:121-130 (ChopPrecision,
  LinearSolveChopPrecision, auto-sync, `< WP` invariant) -> ChopDigits /
  MatchTol defaults + install-time invariant E4.  Lesson: LinearSolveChop MUST
  track ChopPrecision unless explicitly overridden, or matching solves use a
  stale zero test after a precision change.
- DiffExp/State.m:123 RationalizationTolerance (default 10^-10) and its readers:
  - exponent/structure uses (DiffExp/SeriesOps.m:96; DiffExp/
    SingularityDecomposition.m:54,106,259,284,336; DiffExp/
    RegularizedIntegration.m:97,259,266,274,287,368,496) — DELETED WITH THE
    REPRESENTATION: exact tags mean nothing is reconstructed numerically (D2).
  - value-snapping uses (FeynmanTrick/DiffExpIntegration.m:151-202) -> SnapTol.
  Lesson: 10^-10 absolute at 500-digit precision was simultaneously far too
  loose for snapping (can snap distinct nearby roots together) and meaningless
  as an exponent test; wp-derived SnapTol + exact-confirmation policy (a snap
  must be CONFIRMED exactly when the snapped object is later used as exact —
  the consumer spec for LoadSystem/Transport carries that rule).
- FeynmanTrick/DiffExpIntegration.m:1177-1190 (zeroCoeffQ + LaurentTrim) ->
  LaurentLeadTol + NumericallyZeroQ.  Lesson (legacy review finding 13a):
  misclassifying ONE leading coefficient shifts the whole Laurent window AND the
  uniform per-level prefactor shift (ShiftRawBoundariesToFinite,
  DiffExpIntegration.m:1258-1280) — this single absolute test silently
  reshaped entire levels; hence relative test + ambiguity abort.
- DiffExp/Wronskian.m:63, IntegrationStrategies/Helpers.m:32,84, VOP.m:120,248,
  RegularizedIntegration.m:1047 (rank/nullspace tolerances) -> RankTol.
  Lesson: the fitter's relative `10^(-prec/2)` SVD cut was the only relative
  test in the old code and the only one that never misfired in the campaign;
  tightened to wp/4 to open the [SnapTol, RankTol] gray zone (REVIEW-math
  D15) — the fitter itself is deleted.
- DiffExp/State.m:209-224 (`I*` constants), DiffExp/LineSegmentation.m:109-113,
  121,132 (DigitsNeeded, safety-subtract geometry), DiffExp/Transport.m:527,
  759-760,776-841,913,1191-1213 -> the SURVIVING operational constants
  ($SafetyDigits, $MinExpansionOrder, EvalErrorSeriesDecrease); the six
  adaptive-search constants are DELETED with their machinery (DEC-11;
  REVIEW-minimalism 8; ledger waiver "expansion-order adaptive search").
  Lessons kept: per-segment digit budget needs `Ceiling[Log10[#segments]]`
  (predivision two-pass, Transport.m:666-760); error-probe order reduction
  grows with coupling depth (top `(couplingDepth-1)` t-orders of chained
  particular solutions are degraded — MaxCouplingOrder itself is subsumed by
  exact polynomial recursion + ErrorEstimate, DEC-9, but the probe-order
  decrement survives as EvalErrorSeriesDecrease).
- DiffExp/State.m:212 ICheckMultivaluedChop (= 5, used at Transport.m:481 to
  chop before the SingularityCheck) — DROPPED WITH REASON: in DiffExp2 leftover
  multivaluedness is a STRUCTURAL question on exact tags (is there a sector
  with b != 0 or p > 0 at a chart with no prescription?), not a numeric one; no
  tolerance needed (see Config.m spec, Prescriptions).
- DiffExp/State.m:149 CrosscheckChopPrecision (= 30) -> ResidTol (relative).
- DiffExp/LocalSeries.m:480,492 (`"Tolerance" -> 10^-10` option defaults of the
  finite-width solver) — the solver is subsumed by Solve.m (RewritePlan I2);
  any residual check it needs uses ResidTol; no per-function tolerance options
  exist in DiffExp2.
- DiffExp/Pade.m:34,70,80 + DiffExp/Mobius.m:39,48,66 + DiffExp/Transport.m:
  541-562 (precision floors and the 2x input raise) -> $InputPrecisionFactor,
  $MaxExtraPrecisionValue, $MinPrecisionFloor semantics.  Lesson: the 2x
  headroom is what survives chart-map composition; $MinPrecision floors stop
  Mathematica's significance arithmetic from bleeding digits in long
  evaluation chains.

## 7. DEPENDENCIES

May call: NOTHING (System` only).  May be called by: every other DiffExp2
module and (post-M5) the FT layer.  In particular Tolerances.m MUST NOT read
configuration, files, or global mutable state other than its own installed
record.  This is what makes the acyclic order sound.

## 8. UNIT TESTS (Tests/test_tolerances.m; names binding)

Error assertions ("aborts"/"ABORTS" below) are `Catch[..., "DiffExp2Error"]`
returning the `Failure["DiffExp2", ...]` with the named "ID"/"Module" payload
fields (DEC-1; REVIEW-math D17) — not CheckAbort, not message assertions.

T1  `test_chop_digits_default`: `ChopDigits[500] === 250`, `ChopDigits[100] === 50`.
T2  `test_chop_floor_exact`: `ChopFloor[250] === 10^-250` and is an exact
    Rational (Head === Rational), not a Real.
T3  `test_match_tol_sync_default`: install directly via
    `InstallToleranceState` with the record Config WOULD produce at wp=500
    defaults (the Config.md §3 table is the source of the expected values; NO
    Config call — this suite depends on NOTHING, section 7, and must run
    before Config.m exists; REVIEW-minimalism 31); assert
    `Tol["MatchTol"] === Tol["ChopFloor"] === 10^-250`
    (DiffExp/DiffExp.m:126 sync behavior).  The Config-side half of the sync
    is tested in test_config.m, not here.
T4  `test_snap_rank_derivation`: `SnapTol[500] === 10^-250`;
    `RankTol[500] === 10^-125` (strictly looser — the E7 gray zone is
    non-empty, REVIEW-math D15); `GeomGuardTol[500] === 10^-250`;
    `ChopReserve[500, 250] === 250`; `LaurentLeadTol[250] === 10^-24` (the
    floor); `LaurentLeadTol[40] === 10^-20` (the Max picks the derived value
    below chopDigits 48); `ResidTol[500] === 10^-50`.
T5  `test_install_schema_loud`: InstallToleranceState with (a) a missing key,
    (b) an extra key `"Foo"`, (c) `"ChopDigits" -> 250.0` (Real) each abort with
    messages naming the offender (E3).
T6  `test_install_chop_ge_wp_loud`: `"ChopDigits" -> 500` at
    `"WorkingPrecision" -> 500` aborts with the E4 text.
T7  `test_tol_unknown_key_loud`: `Tol["ChopPrecision"]` (old name!) aborts
    listing the nine valid names (E2) — guards against silent old-name reads.
T8  `test_tol_before_install_loud`: fresh kernel state, `Tol["ChopFloor"]`
    aborts (E1).
T9  `test_numerically_zero_band`: with tol = 10^-125, scale = 1:
    `NumericallyZeroQ[10^-300, 1, 10^-125, "t"] === True`;
    `NumericallyZeroQ[10^-40, 1, 10^-125, "t"] === False`;
    `NumericallyZeroQ[10^-125, 1, 10^-125, "t"]` ABORTS quoting "t" (E5);
    band edges: `10^-130` -> True, `10^-120` -> False (closed-form: the
    ERROR band is the closed interval `[tol/10^4, tol*10^4] =
    [10^-129, 10^-121]`; the True/False regions are the open exteriors —
    REVIEW-math D29 wording).  Floor exemption (DEC-2/DEC-3):
    `NumericallyZeroQ[10^-29, 1, 10^-24, "t"] === True` and
    `NumericallyZeroQ[10^-22, 1, 10^-24, "t"] === False`, with NO error in
    either case — at the 10^-24 floor the test is binary, so the two
    campaign populations (~1e-29-relative residues, O(1)-relative signal)
    both classify cleanly.
T10 `test_numerically_zero_exact_and_symbolic`:
    `NumericallyZeroQ[0, 0, tol, "t"] === True` (exact zero, zero scale);
    `NumericallyZeroQ[2 - 2, anything...] === True`;
    `NumericallyZeroQ[Sym, 1, tol, "t"] === False` for an undefined symbol;
    `NumericallyZeroQ[10^-300, 0, tol, "t"] === False` (zero scale blocks
    numeric trimming).
T11 `test_laurent_trim_lesson_relative`: the FT regression in miniature —
    coefficient `10^-45` with window scale `10^-44` (tiny but REAL leading
    term, relative size 10^-1) classifies False under
    `LaurentLeadTol[250] = 10^-24` (binary at the floor: 10^-1 >> 10^-24),
    where the old absolute 10^-40 test (DiffExpIntegration.m:1178) would
    have trimmed it and shifted MinPower.
T12 `test_operational_constants_pinned`: `$SafetyDigits === 2`,
    `$MinExpansionOrder === 10`, `$InputPrecisionFactor === 2`,
    `$MaxExtraPrecisionValue === 1000`, `$AmbiguityBandDecades === 4`,
    `$NearSingularityGuardDecades === 6`,
    `EvalErrorSeriesDecrease[1] === 3`, `EvalErrorSeriesDecrease[5] === 6`
    (Ceiling[0.7*5]+2 = 6; DiffExp/State.m:222 formula); additionally the
    six deleted adaptive-search symbols (section 2) do NOT exist in the
    `DiffExp2`Tolerances` context (DEC-11; REVIEW-minimalism 8).
T13 `test_exports_visible_cross_context`: from a scratch context with an empty
    $ContextPath, every symbol in section 2 evaluates under its fully
    qualified name (e.g. ``DiffExp2`Tolerances`ChopDigits[500]`` returns 250,
    not unevaluated) — the wolfram-package-context-traps regression.

## 9. LINE BUDGET

~100 lines (RewritePlan 3.2).  PRE-AUTHORIZED at M0 review
(REVIEW-minimalism 4): the six dead adaptive-search constants are DELETED,
not folded (-10), funding DE2Error (+~12, REVIEW-minimalism 5).  Estimated:
derivations ~30 (incl. InputSnapTol/GeomGuardTol/ChopReserve),
NumericallyZeroQ ~20, DE2Error ~12, install/Tol state ~25, constants ~7,
usage strings ~15 — total ~109 against the restated ceiling.  If over
budget, cut in this order: (1) drop `ToleranceStateInstalledQ` (API.m can
probe via a Check on `Tol`); (2) shorten error messages to the mandatory
fields only.  (The former cut "fold the expansion-order constants into
$TransportTuning" is gone — nothing left to fold, REVIEW-minimalism 8.)
DO NOT cut: the ambiguity band, the floor exemption, the install-time schema
check, or any usage declaration (export visibility is a correctness property
here, not documentation).

## 10. OPEN QUESTIONS

Q1 RESOLVED (DEC-4): for combined/merged objects the cancellation-gate scale
   is the max |coefficient| of the merged combination at that eps order —
   global per merged object per eps order, not per-master.  For a single
   Laurent object the scale stays "max |coeff| over the SAME object".  No
   open question remains; M5 ladder data is observational only.
Q2 `ResidTol = 10^(-wp/10)` is a judgment call (old absolute 30 digits at
   wp=500 corresponds to wp/16.7).  Any value in [wp/16, wp/8] separates the
   observed populations; revisit only if M3 closed-form units show residuals
   above 10^-50 for correct solutions at wp=500.
Q3 Should `SnapTol` also govern the dedup of EXACT algebraic singularity
   positions in Transport.m segmentation (old: DeleteDuplicatesBy
   N[#, ChopPrecisionVal], DiffExp/LineSegmentation.m:77,97)?  Current answer:
   numeric dedup at SnapTol is allowed only as a PRE-filter, with exact
   `RootReduce`-level confirmation required before two roots are merged;
   Transport.m spec must repeat this.  Flagged here because it is a place a
   silent absolute test could creep back in.
Q4 The FT layer keeps its own `$FTConfig` until M5/M6; its
   "IntegrationPoleAllowance" (FeynmanTrick/FeynmanTrick.m:50) is an ORDER
   budget, not a tolerance, and stays out of this module.  Confirm at M5 that
   no FT numeric literal other than the ones absorbed above survives
   (grep gate: `10\^-` under FeynmanTrick/ should hit only print statements).
