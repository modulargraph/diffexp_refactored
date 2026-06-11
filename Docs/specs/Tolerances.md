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
semantics (chopFloor, matchTol, snapTol, rankTol, laurentLeadTol, residTol), is
DERIVED from WorkingPrecision (never a free literal at a call site), and zero
tests that classify structure (leading-coefficient trimming, rank decisions) are
RELATIVE with a loud-error ambiguity band, so a borderline value can never be
silently classified either way.  The module also centralizes the old library's
operational safety constants (the `I*` constants of DiffExp/State.m:209-224) so
that no tuning number lives anywhere else in DiffExp2.

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
  (FeynmanTrick/BoundaryConditions.m:191,195).
- `RankTol[wp_Integer] -> 10^(-Floor[wp/2])`
  RELATIVE rank/nullspace threshold: a singular value (or pivot) counts as zero
  iff `sv < RankTol[wp] * svMax`.  Replaces the absolute
  `NullSpace[..., Tolerance -> 10^-LinearSolveChopPrecisionVal]` sites
  (DiffExp/Wronskian.m:63, IntegrationStrategies/Helpers.m:32,84, VOP.m:120,248)
  and the fitter's `Min[svals] >= 10^(-activeNumericPrecision[]/2) * ...` test
  (DiffExp/RegularizedIntegration.m:1047).  Rank decisions create/destroy
  structure (disease D1) and therefore MUST be relative and MUST go through
  `NumericallyZeroQ` (ambiguity band) — never a raw comparison.
- `LaurentLeadTol[chopDigits_Integer] -> 10^(-Floor[chopDigits/2])`
  RELATIVE leading-coefficient zero test for eps-Laurent window trimming.  This
  is the replacement demanded by RewritePlan 3.2 for the FT `zeroCoeffQ`
  ABSOLUTE `10^-40` test (FeynmanTrick/DiffExpIntegration.m:1177-1178, consumed
  by LaurentTrim at 1180-1190).  Semantics: leading coefficient `c` of a Laurent
  array with window scale `s = Max[Abs /@ allWindowCoefficients]` is zero iff
  `NumericallyZeroQ[c, s, LaurentLeadTol[chopDigits], context]`.  `s == 0` =>
  only exact zeros trim.  Rationale for the midpoint exponent: post-chop noise
  sits at relative `<= 10^-chopDigits`, real signal at relative `>= ~10^-40`
  empirically (the campaign's constant); `chopDigits/2` (default 125 decades)
  separates both by >= 80 decades, and the ambiguity band (below) turns any
  violation of that separation into a loud error instead of a wrong window.
- `ResidTol[wp_Integer] -> 10^(-Floor[wp/10])`
  RELATIVE threshold for the always-on ODE-residual spot-check (RewritePlan 3.1
  invariants) and for construction-time cross-checks.  Replaces the old absolute
  30-digit `CrosscheckChopPrecision`/CPChop (DiffExp/State.m:149,
  DiffExp/Utilities.m:105): correct solutions give residuals at the noise floor
  (`~10^-chopDigits` relative), wrong ones give O(1); wp/10 (default 50 decades)
  sits far from both.

The shared classification predicate:

- `NumericallyZeroQ[c_, scale_, tol_, context_String] -> True | False | LOUD ERROR`
  THE one zero-classification routine for all relative tests in DiffExp2.
  Contract: (1) `PossibleZeroQ[c]` exactly True -> True.  (2) `c` not NumericQ
  -> False (symbolic content is never "numerically zero"; callers handle
  symbols structurally).  (3) `scale == 0` (no reference magnitude) -> False
  unless branch (1) fired.  (4) Otherwise let `r = Abs[N[c, dig]]` with `dig`
  the installed ChopDigits: `r < tol*scale/10^bandDecades -> True`;
  `r > tol*scale*10^bandDecades -> False`; ELSE LOUD ERROR (section 5) quoting
  `context` verbatim — `context` MUST already contain the chart / sector tag
  (a,b,p) / eps-order / t-order identifiers the caller has.  `bandDecades`
  is the exported constant `$AmbiguityBandDecades = 10`.
  This predicate is what makes "IntegerQ-on-floats" and threshold-straddling
  misclassification (disease D1; the LaurentTrim MinPower-shift hazard of
  legacy-review finding 13a, Docs/reviews/rewrite_plan_review_3lens.json) impossible.

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
- `$SafetyExpansionSubtract = 5` (ISafetyExpansionSubtract, DiffExp/State.m:216;
  matrix-truncation error read at order
  `ExpansionOrder - $SafetyExpansionSubtract - (couplingDepth - 1)`,
  DiffExp/LineSegmentation.m:109-113).
- `$ExpansionOrdersAveraging = 3`, `$ExpansionOrderIncrease = 10`,
  `$ExpansionOrderDecrease = 10`, `$ExpansionOrderIncreaseValidate = 25`,
  `$DigitsSurplusDecrease = 3`, `$MinExpansionOrder = 10`
  (DiffExp/State.m:217-224; the AccuracyGoalValidate "Before" adaptive search,
  DiffExp/Transport.m:776-841, and "After" redo, Transport.m:1191-1213).
- `EvalErrorSeriesDecrease[couplingDepth_Integer] := Ceiling[0.7*couplingDepth] + 2`
  (ICurrEvalErrorSeriesDecrease, DiffExp/State.m:222; the two-point error probe
  order reduction, DiffExp/Transport.m:913).
- `$AmbiguityBandDecades = 10` (new; see NumericallyZeroQ).

## 3. DATA CONTRACTS

- Tolerance state record (produced by Config.m, consumed here, read by every
  higher module through `Tol`):
  `<|"WorkingPrecision"->Integer, "ChopDigits"->Integer, "MatchDigits"->Integer,
    "ChopFloor"->Rational, "MatchTol"->Rational, "SnapTol"->Rational,
    "RankTol"->Rational, "LaurentLeadTol"->Rational, "ResidTol"->Rational|>`
  All tolerance values exact `10^-k` rationals (never machine reals — they are
  compared against arbitrary-precision numbers).
- Consumers and which name they may use (binding; anything else is a spec
  violation found in review):
  - EpsSeries.m: `LaurentLeadTol` (window trimming via NumericallyZeroQ only),
    `ChopFloor` (coefficient chop after add/mul/div).
  - SectorSeries.m: `ChopFloor`, `ResidTol` (evaluate/re-expand checks),
    `SnapTol` (none expected; flag in review if used).
  - Indicial.m: NONE of the numeric tolerances for exponents — indicial data is
    EXACT (RewritePlan I1); only `ChopFloor` for numeric matrix entries.
    Exponent rationalization tolerances are DELETED (section 6).
  - Solve.m: `MatchTol`, `RankTol`, `ChopFloor`, `ResidTol`.
  - Transport.m: `MatchTol`, `RankTol`, `ResidTol`, `$SafetyDigits`,
    `$InputPrecisionFactor`, the expansion-order constants,
    `EvalErrorSeriesDecrease`.
  - Integrate.m: `ChopFloor`, `LaurentLeadTol` (assembled-combination
    cancellation checks at the object level), `SnapTol` (endpoint snapping).
  - API.m / FT layer: `SnapTol` (chart-map endpoint snapping, singularity
    dedup), `ChopFloor`.
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

Every error is raised via the shared DiffExp2 loud-error primitive (Abort with
printed context; the API.m spec owns its exact form) and must carry the strings
shown.  Enumerated conditions:

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
  test in the old code and the only one that never misfired in the campaign.
- DiffExp/State.m:209-224 (`I*` constants), DiffExp/LineSegmentation.m:109-113,
  121,132 (DigitsNeeded, safety-subtract geometry), DiffExp/Transport.m:527,
  759-760,776-841,913,1191-1213 -> the exported operational constants.  Lessons:
  per-segment digit budget needs `Ceiling[Log10[#segments]]` (predivision
  two-pass, Transport.m:666-760); error-probe order reduction grows with
  coupling depth (top `(couplingDepth-1)` t-orders of chained particular
  solutions are degraded — MaxCouplingOrder, DiffExp/MatrixLoading.m:383).
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

T1  `test_chop_digits_default`: `ChopDigits[500] === 250`, `ChopDigits[100] === 50`.
T2  `test_chop_floor_exact`: `ChopFloor[250] === 10^-250` and is an exact
    Rational (Head === Rational), not a Real.
T3  `test_match_tol_sync_default`: with a state installed via a minimal Config
    round-trip at wp=500 and no explicit overrides,
    `Tol["MatchTol"] === Tol["ChopFloor"] === 10^-250`
    (DiffExp/DiffExp.m:126 sync behavior).
T4  `test_snap_rank_derivation`: `SnapTol[500] === RankTol[500] === 10^-250`;
    `LaurentLeadTol[250] === 10^-125`; `ResidTol[500] === 10^-50`.
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
    band edges: `10^-136` -> True, `10^-114` -> False (closed-form:
    band is `[tol/10^10, tol*10^10] = [10^-135, 10^-115]`, exclusive).
T10 `test_numerically_zero_exact_and_symbolic`:
    `NumericallyZeroQ[0, 0, tol, "t"] === True` (exact zero, zero scale);
    `NumericallyZeroQ[2 - 2, anything...] === True`;
    `NumericallyZeroQ[Sym, 1, tol, "t"] === False` for an undefined symbol;
    `NumericallyZeroQ[10^-300, 0, tol, "t"] === False` (zero scale blocks
    numeric trimming).
T11 `test_laurent_trim_lesson_relative`: the FT regression in miniature —
    coefficient `10^-45` with window scale `10^-44` (tiny but REAL leading
    term, relative size 10^-1) classifies False under
    `LaurentLeadTol[250] = 10^-125`, where the old absolute 10^-40 test
    (DiffExpIntegration.m:1178) would have trimmed it and shifted MinPower.
T12 `test_operational_constants_pinned`: `$SafetyDigits === 2`,
    `$SafetyExpansionSubtract === 5`, `$ExpansionOrdersAveraging === 3`,
    `$ExpansionOrderIncrease === 10`, `$ExpansionOrderDecrease === 10`,
    `$ExpansionOrderIncreaseValidate === 25`, `$DigitsSurplusDecrease === 3`,
    `$MinExpansionOrder === 10`, `$InputPrecisionFactor === 2`,
    `$MaxExtraPrecisionValue === 1000`, `$AmbiguityBandDecades === 10`,
    `EvalErrorSeriesDecrease[1] === 3`, `EvalErrorSeriesDecrease[5] === 6`
    (Ceiling[0.7*5]+2 = 6; DiffExp/State.m:222 formula).
T13 `test_exports_visible_cross_context`: from a scratch context with an empty
    $ContextPath, every symbol in section 2 evaluates under its fully
    qualified name (e.g. ``DiffExp2`Tolerances`ChopDigits[500]`` returns 250,
    not unevaluated) — the wolfram-package-context-traps regression.

## 9. LINE BUDGET

~100 lines (RewritePlan 3.2).  Estimated: derivations ~25, NumericallyZeroQ
~20, install/Tol state ~25, constants ~15, usage strings ~15.  If over budget,
cut in this order: (1) fold the six expansion-order constants into one exported
association `$TransportTuning` (saves ~8 lines of usage text); (2) drop
`ToleranceStateInstalledQ` (API.m can probe via a Check on `Tol`); (3) shorten
error messages to the mandatory fields only.  DO NOT cut: the ambiguity band,
the install-time schema check, or any usage declaration (export visibility is a
correctness property here, not documentation).

## 10. OPEN QUESTIONS

Q1 `LaurentLeadTol` window scale: spec says "max |coeff| over the SAME Laurent
   object".  For multi-master combined objects, should the scale be per-master
   or global per level?  Per-object is the default; M5 ladder data (box_bubble
   level boundaries) should confirm no cross-master scale disparity > 10^80
   (which would put real coefficients in the band).  Owner decision if it does.
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
