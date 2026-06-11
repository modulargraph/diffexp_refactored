# Module spec: DiffExp2/Config.m

Status: M0 deliverable per Docs/RewritePlan.md section 3.2 / section 6 (M0 task 4-11).
Context name: ``DiffExp2`Config` ``.  Line budget: ~150 (see section 9).
Implementable from: this file + Docs/RewritePlan.md + the old-code citations below.
Module order (acyclic, hard rule): Tolerances < Config < EpsSeries < SectorSeries <
Indicial < Solve < Transport/Integrate < API.  Config may call ONLY Tolerances.

---

## 1. PURPOSE

Config.m is the single validated option store for DiffExp2: a fixed,
string-keyed schema with typed validators, a read accessor that makes an
unknown or unset key a LOUD ERROR at call time, and canonicalization of all
symbol-valued inputs (kinematic variables, line parameter) into one pinned
context.  It exists to make the FEC-class silent-miss bug impossible BY
CONSTRUCTION: in the old library the configuration was a raw Association with
mixed Symbol/String keys (DiffExp/State.m:106-136) read through the
unvalidated alias `FEC := DiffExpConfiguration` (DiffExp/State.m:144), so a
bare-symbol key written from a package whose context path lacks
``DiffExp`State` `` parses into the WRONG context and the lookup silently
returns `Missing["KeyAbsent", ...]` — live today at
FeynmanTrick/DiffExpIntegration.m:98,109,127,924, where bare
`RationalizationTolerance` resolves to
``FeynmanTrick`DiffExpIntegration`Private`RationalizationTolerance`` and every
read misses (the guarded variants at 154-155,182-183 then silently substitute
a hard-coded 10^-40, ignoring whatever the user configured; the file's own
author comment at 339-341 — "All config keys must be fully qualified ... since
this package's BeginPackage restricts the context path" — documents the trap
for writes but the reads were still missed).  This is the same disease class
as the unexported-symbol no-op of memory/wolfram-package-context-traps.md.
Config.m also owns cross-field invariants (ChopPrecision/WorkingPrecision,
LinearSolveChopPrecision sync), DeltaPrescriptions parsing/validation, and the
installation of the tolerance state into Tolerances.m; it does NOT load
matrices or touch solver state (the old UpdateConfiguration's hidden
LoadMatrices call and cache wipes, DiffExp/DiffExp.m:171-190, move to
API.m/LoadSystem).

## 2. PUBLIC SYMBOLS

All exported with `::usage` before `Begin["`Private`"]` (cross-context
visibility test C18 mirrors Tolerances T13).

- `LoadConfiguration[rules : {___Rule} | __Rule | _Association] -> Association`
  Reset every key to its schema default, then apply `rules` as
  UpdateConfiguration.  (Old semantics DiffExp/DiffExp.m:95-99, minus the
  DiffExpExtensions hook list, State.m:138/DiffExp.m:97,192 — no extension
  mechanism in v1, dropped with reason.)  Returns `CurrentConfiguration[]`.
- `UpdateConfiguration[rules : {___Rule} | __Rule | _Association] -> Association`
  Validate-and-merge: every key canonicalized (below), checked against the
  schema, value validated, cross-field invariants re-checked, tolerance state
  recomputed and installed via ``Tolerances`InstallToleranceState``.  Any
  failure aborts BEFORE any mutation (all-or-nothing update; the old code
  merged first and validated after, DiffExp/DiffExp.m:105, so a failed update
  left half-applied state).
- `CurrentConfiguration[] -> Association`
  Full copy of all known keys and current values, String-keyed.  (Old version
  returned a filtered subset, DiffExp/DiffExp.m:86-92, which hid live keys like
  "EstimateError" and "AccuracyGoalValidate" from the user.)
- `CFG[key_] -> value`
  THE validated read accessor (replaces every `FEC[...]` read).  `key` may be a
  String or a Symbol; Symbols are canonicalized by `SymbolName` (so
  ``DiffExp2`Config`WorkingPrecision``, `System`WorkingPrecision`, and a
  stray `SomePkg`Private`WorkingPrecision` all reach "WorkingPrecision" — and
  if the NAME is unknown the error prints the symbol's full context, which is
  exactly the diagnostic the old silent miss withheld).  Unknown key -> LOUD
  ERROR.  Never returns `Missing[...]`, never returns unevaluated.
- `ConfigSchema[] -> Association` Read-only copy of the schema (for tests,
  the FT shim audit, and ExportDisposition checks).
- `CanonicalKey[key_String | key_Symbol] -> String` The canonicalization used
  by CFG/Update (exported for the FT shim and tests).
- `PinnedVariable[sym_Symbol] -> Symbol`
  Variable-context pinning: returns `Symbol[SymbolName[sym]]` evaluated in the
  pinned context `Global`` (the FIRE/Fermat rule: FIRE requires plain Global`
  symbols, FeynmanTrick/FIREInterface.m:203-226; the FT layer re-derives "the
  detected variable" by SymbolName for the same reason,
  FeynmanTrick/DiffExpIntegration.m:371-378).  All Variables/LineParameter
  values are passed through this at update time, so two symbols that agree in
  name but differ in context can never coexist in DiffExp2 state.
- `ConfiguredQ[] -> True|False` True after the first successful
  Load/UpdateConfiguration (used by API.m for its "no configuration loaded"
  error path).

Note: there is deliberately NO exported mutable association and NO `FEC`-style
alias.  Higher modules read through `CFG` only; Tests assert (C19) that the
string "DiffExpConfiguration" does not occur outside Config.m.

## 3. DATA CONTRACTS

### 3.1 Schema entry (internal shape, fixed at load; one per key)

`<|"Type" -> validator, "Default" -> value, "Normalize" -> f | None|>` where
`validator` is a predicate over the NORMALIZED value.  The schema is a
constant; nothing may add keys at runtime (no plugin keys in v1).

### 3.2 The configuration schema (all keys; String-keyed)

Kept keys, their DiffExp2 meaning, type, default, and the old-code anchor:

| Key | Type / values | Default | DiffExp2 meaning (consumer) | Old anchor |
|---|---|---|---|---|
| "WorkingPrecision" | Integer >= 20 | 500 | master precision; drives ALL tolerances | State.m:135 |
| "ChopPrecision" | Integer or Automatic | Automatic -> `Tolerances`ChopDigits[wp]` | chopDigits for Tolerances; must be < WP | State.m:109; DiffExp.m:121-125 |
| "LinearSolveChopPrecision" | Integer or Automatic | Automatic -> ChopPrecision | matchDigits; auto-syncs to ChopPrecision unless explicitly set IN THE SAME update | State.m:119; DiffExp.m:126-130 |
| "AccuracyGoal" | Integer or "?" | "?" | requested digits at evaluation points; digit budgeting + adaptive order search (Transport.m module) | State.m:107; old Transport.m:527,759 |
| "AccuracyGoalValidate" | "Before" \| "After" \| None | "Before" | "Before": adaptive expansion-order search (averaging window $ExpansionOrdersAveraging, +/-10 steps, 3-digit surplus, min order 10; old Transport.m:776-841); "After": redo segment at +$ExpansionOrderIncreaseValidate (old Transport.m:1191-1213); None: fixed order | State.m:108 |
| "DeltaPrescriptions" | list, see 3.3 | {} | branch-point prescriptions feeding LocalSolution "Prescriptions" | DiffExp.m:140-148,164-168 |
| "DivisionOrder" | Integer >= 2 | 3 | chart geometry k: match/eval point at radius/k of BOTH adjacent charts (GetCPL/GetCPR lesson, old Mobius.m:98-142; old Transport.m:679); FT pins k=4 (DiffExpIntegration.m:272) | State.m:112 |
| "EpsilonOrder" | Integer >= 0 | 4 | requested top eps order; EpsWindow tracks honest completeness against it | State.m:113 |
| "EstimateError" | False \| True \| "Fast" | "Fast" | two-point error probe per segment, additive accumulation, abort > 1 (old Transport.m:907-925,980-993); feeds "ErrorEstimate" of LocalSolution and "ErrorEstimates" output | State.m:114 |
| "ExpansionOrder" | Integer >= `Tolerances`$MinExpansionOrder` | 50 | t-truncation order of chart series | State.m:115; old Transport.m:526 |
| "LineParameter" | Symbol | Global`x | user-facing line-parameter NAME (pinned via PinnedVariable); must not be in Variables | State.m:120; DiffExp.m:131-138 |
| "LogFile" | String or None | None | append session output to file | DiffExp.m:107-114 |
| "MatrixDirectory" | String | "" | data only; consumed by API.m LoadSystem (full-format `d<var>_full.m` primary, legacy slices for parity transport only, per RewritePlan I1) — Config performs NO loading | State.m:121; old auto-load DiffExp.m:171-179 WAIVED |
| "RadiusOfConvergence" | positive rational | 1 | chart-coordinate rescaling keeping high-order coefficients O(1); banana REQUIRES 10 (ledger seed; old Mobius.m:47,65) | State.m:122 |
| "RationalizationTolerance" | positive rational or Automatic | Automatic -> `Tolerances`SnapTol[wp]` | REDEFINED: snap/dedup tolerance at the FT seam ONLY (chart-map endpoint snapping, singular-factor root dedup; old FeynmanTrick/DiffExpIntegration.m:151-202). Exponent-rationalization uses are DELETED (exact tags, RewritePlan I1/I2) | State.m:123 |
| "SegmentationStrategy" | "Predivision" | "Predivision" | v1 supports Predivision ONLY; "Dynamic" -> loud error E10 (its machinery, old LineSegmentation.m:116-145, is not ported; RewritePlan 3.2) | State.m:128 |
| "UseMobius" | Boolean | False | Mobius chart maps for classic transport (banana classic line); Integrate.m REJECTS Mobius charts loudly (old FT hard-required False: "Required for integration!", DiffExpIntegration.m:352) | State.m:130 |
| "UsePade" | Boolean | False | Pade evaluation accelerator in SectorSeries.m (ported per RewritePlan 3.2; old GetPade applied per-Logx-coefficient, Pade.m:32-53; its silent fallback Pade.m:44-46 becomes loud there) | State.m:131 |
| "Variables" | list of Symbols | {} | kinematic invariants/masses; pinned via PinnedVariable; auto-detected by LoadSystem from filenames when {} (old MatrixLoading.m:27-55,88) — see E12 for the overwrite rule | State.m:132 |
| "Verbosity" | Integer >= 0 | 1 | PrintInfo gate (old Utilities.m:60) | State.m:133 |
| "VerbosityDebug" | Integer >= 0 | 0 | PrintDebug gate (old Utilities.m:56) | State.m:134 |
| "CrosscheckLevel" | Integer >= 0 | 0 | enables OPTIONAL extra cross-checks by level threshold (old DiffExp.m:115-117, State.m:196-206). DiffExp2's always-on invariants (RewritePlan 3.1) are NOT gated by this and cannot be disabled | State.m:110 |
| "CrosscheckFlags" | list of Strings from the published check registry | {} | enable named optional checks regardless of level (old DiffExp.m:118-120). Unknown flag name -> loud error. Old flag "SingularityCheck" (default-OFF at level 0, State.m:205) is SUPERSEDED: prescription consistency is an always-on construction check in DiffExp2 (RewritePlan 3.1 "Prescriptions"; legacy review finding 3) and the flag name is rejected with a pointer | State.m:196-207 |
| "SaveExpansionsCompress" | Boolean | False | compress saved segment series (old Transport.m:852-865,1278-1287) | State.m:124 |
| "SaveExpansionsCompressDirectory" | String or None | None | spill compressed segments to files (old Transport.m:854-855,1070-1071; old code probed this key with `!FEC[...] === "?"` against an ABSENT key — declared default kills that pattern, see F-h) | Transport.m:854 |
| "SaveExpansionsOrder" | Integer or None | None | truncation order applied when saving (old Transport.m:853,1069: KeyExistsQ-gated) | Transport.m:853 |
| "AbortOnAnalyticContinuationFail" | Boolean | True | True: missing/conflicting prescription at a chart that needs one is a hard abort; False: downgrade to warning + flagged, FT pipeline mode (old Transport.m:470-478; FT sets False at DiffExpIntegration.m:413,468). Old default was IMPLICIT (key absent -> `Missing =!= False` -> abort-ish branch); now declared | Transport.m:471 |

Dropped keys — `UpdateConfiguration`/`CFG` reject each with a DEDICATED error
(E9) naming the replacement, never the generic unknown-key error:

| Key | Reason |
|---|---|
| "IntegrationStrategy" ("Default"/"VOP") | strategy stack deleted (RewritePlan I2); ONE solver | State.m:129, IntegrationStrategies/* |
| "UseRationalRecurrence" | ditto; the denominator-cleared recursion is the only path in Solve.m | State.m:127 |
| "InvWronskSolver" | Wronskian machinery deleted | State.m:117 |
| "HomogeneousSolve" ("Expand"/"DontExpand") | strategy-stack switch; symbolic-eps Frobenius has one homogeneous path (FT passes "DontExpand" today, DiffExpIntegration.m:276 — shim audit removes it at M5) | State.m:125 |
| "KeepMatrixExpansions" | cache policy is internal to DiffExp2 (old MatrixLoading.m:236-249) | State.m:118 |
| "Parallel" | no consumer in the old core (grep-verified: defined only at State.m:126) | State.m:126 |
| "IgnoreIndicialCheck" | suppressed the indicial-root warning (old Frobenius.m:42); the DiffExp2 I1 contract violation (char poly not factoring as Π(λ−a−b eps)) is a HARD error and must not be suppressible (RewritePlan R1) | State.m:116, Frobenius.m:42 |

### 3.3 DeltaPrescriptions normalized form

Accepted input elements (old parser DiffExp/DiffExp.m:140-148): either
`poly ± I*Global`δ` (the δ-coefficient/I gives the sign) or `{poly, sign}`.
Normalized stored form: list of `{poly, sign}` with `sign ∈ {+1, -1}` and
`poly` an irreducible polynomial in the PINNED kinematic variables
(irreducibility check via FactorList as in DiffExp/DiffExp.m:164-168: more
than one nontrivial factor -> E6; message text "Physical singularities should
be irreducible polynomials!").  The USER list is stored separately from the
effective list (old UserDeltaPrescriptions, DiffExp.m:147, State.m:191) so
that LoadSystem can union in sqrt-derived auto-prescriptions
(old State.m:227-230 SquareRootPrescriptionsAdded; MatrixLoading.m:181-190
assigns +iδ to new sqrt factors with an info print) without clobbering user
data on re-parse.  Downstream contract: Transport.m derives from this list the
per-chart record required by RewritePlan 3.1 (quoted verbatim):

```
"Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
            "LeadingCoeffSign"|> with the derived chart Im-sign;
            consistency-checked at construction (even multiplicity = no
            constraint; conflict or missing prescription at a chart with
            b!=0/p>0 sectors = LOUD ERROR; sqrt factors auto-prescribed
            as in old State.m DEqnSquareRoots)
```

Config.m's responsibility ends at the validated `{poly, sign}` list +
separate user list; multiplicity/leading-coefficient derivation is
Transport.m's (old PrepareAnalyticContinuation, AnalyticContinuation.m:18-68).

### 3.4 Tolerance state handed to Tolerances.m

On every successful update, Config computes and installs (see the Tolerances.md
data contract; keys verbatim):
`<|"WorkingPrecision" -> wp, "ChopDigits" -> cd, "MatchDigits" -> md,
  "ChopFloor" -> 10^-cd, "MatchTol" -> 10^-md,
  "SnapTol" -> (RationalizationTolerance if numeric, else Tolerances`SnapTol[wp]),
  "RankTol" -> Tolerances`RankTol[wp],
  "LaurentLeadTol" -> Tolerances`LaurentLeadTol[cd],
  "ResidTol" -> Tolerances`ResidTol[wp]|>`
with `cd` = ChopPrecision (resolved from Automatic), `md` =
LinearSolveChopPrecision (resolved per the sync rule).

### 3.5 LocalSolution / Sector / EpsWindow

Config.m neither produces nor consumes these objects (RewritePlan 3.1 is
normative for their shape); it only supplies the scalars (EpsilonOrder,
ExpansionOrder, DivisionOrder, RadiusOfConvergence, prescriptions) from which
other modules build them.

## 4. INVARIANTS (always on, cheap)

- The store always contains EXACTLY the schema's keys (no more, no fewer);
  checked after every update.
- `ChopPrecision < WorkingPrecision` (resolved values; old invariant
  DiffExp/DiffExp.m:121-125) and `LinearSolveChopPrecision <= ChopPrecision`.
- `LineParameter ∉ Variables` by SymbolName (old DiffExp.m:131-135 and
  MatrixLoading.m:52-54 both enforced the two directions; here it is one
  check on canonicalized names).
- Every Symbol stored under "Variables"/"LineParameter" has
  `Context[s] === "Global`"` (post-pinning).
- If "AccuracyGoalValidate" =!= None then "AccuracyGoal" is numeric (the old
  code SILENTLY SKIPPED the order search when AccuracyGoal was "?",
  Transport.m:776-777 `NumericQ` gate — now a cross-field error E8).
- Tolerance state in Tolerances.m is never stale: re-installed atomically on
  every successful update (so `Tol` and `CFG` can never disagree about wp).

## 5. ERROR CONTRACT (no silent fallbacks anywhere)

All errors abort with messages carrying the named fields.  Conditions:

E1  READ of unknown key: "Config: unknown configuration key <name>. Valid keys:
    <schema list>."  When the key was a Symbol, append "(symbol context:
    <Context[sym]> — a non-Global context here usually means a package-path
    symbol-resolution trap; see Docs/specs/Config.md)".  THIS is the
    constructive fix for the FEC-class bug: the failure mode that silently
    returned Missing (FeynmanTrick/DiffExpIntegration.m:98,109,127,924) now
    names itself and its context.
E2  WRITE of unknown key: same message + "no keys were updated" (all-or-nothing).
E3  Type/value violation: "Config: invalid value for <key>: got <value>
    (expected <type/values>)."  Applies to EVERY key; in particular
    "EstimateError" -> anything outside {False, True, "Fast"} is an ERROR —
    the old parser silently coerced unknown values (and the STRING "False")
    to False (DiffExp/DiffExp.m:150-158).
E4  `ChopPrecision >= WorkingPrecision`: old text preserved — "The value of
    ChopPrecision should be smaller than the value of WorkingPrecision."
    (DiffExp/DiffExp.m:123) + both numbers.
E5  `LineParameter` ∈ Variables (by name): old text preserved — "The symbol
    for the line parameter can't be equal to one of the kinematic invariants
    or masses." (DiffExp/DiffExp.m:133).
E6  Reducible DeltaPrescriptions polynomial (3.3).
E7  DeltaPrescriptions sign not ±1 after normalization, or δ appearing
    nonlinearly: error naming the offending element (old parser would compute
    a garbage Coefficient silently, DiffExp.m:142).
E8  "AccuracyGoalValidate" set with non-numeric "AccuracyGoal" (section 4).
E9  Dropped-key write/read: "Config: <key> was removed in DiffExp2: <reason
    from 3.2 table>."  (Migration aid; distinct from E1 so old scripts fail
    with instructions, not confusion.)
E10 "SegmentationStrategy" -> "Dynamic": "Dynamic segmentation is not ported
    in DiffExp2 v1 (RewritePlan 3.2); use Predivision."
E11 `CFG` before any LoadConfiguration: "Config: no configuration loaded."
    (`ConfiguredQ[]` exists so API.m can give friendlier context.)
E12 Variables auto-detect conflict (enforced here, raised by LoadSystem
    through a Config helper): if the user supplied non-empty "Variables" and
    the matrix directory's detected names differ as a SET (by SymbolName), the
    old code silently OVERWROTE the configured value
    (DiffExp/MatrixLoading.m:88 writes FEC[Variables] unconditionally); in
    DiffExp2 this is an error naming both lists.  Empty user list -> detected
    list is adopted (and pinned) as before.
E13 "LogFile" unopenable: error (old code would propagate a failed
    OpenAppend stream into $Output silently, DiffExp.m:108-114).

Forbidden fallbacks (every place a fallback is tempting; each is banned):

F-a Returning `Missing[...]`/unevaluated on a key miss (raw Association
    semantics, old FEC State.m:144).  CFG either answers or aborts.
F-b Accepting and merging unknown keys (old `Merge[{config, assoc}, Last]`,
    DiffExp/DiffExp.m:105, accepted ANY key including mis-contexted symbol
    twins of real keys — the write-side half of the FEC trap).
F-c Coercing invalid option values to a default (the "EstimateError" parser,
    DiffExp/DiffExp.m:150-158).  E3 instead.
F-d Hard-coded tolerance default on read failure
    (`Quiet[Check[FEC[RationalizationTolerance], 10^-40]]`,
    FeynmanTrick/DiffExpIntegration.m:154-155,182-183 — note `10^-40` is not
    even the documented default `10^-10`, so the fallback silently changed
    semantics).  After M5 the FT layer reads through CFG and the Quiet/Check
    wrapper is deleted by the shim audit.
F-e Implicit behavior from ABSENT keys: `FEC["AbortOnAnalyticContinuationFail"]
    === False` (Transport.m:471), `KeyExistsQ[..., "SaveExpansionsOrder"]`
    (Transport.m:853), `!FEC["SaveExpansionsCompressDirectory"] === "?"`
    (Transport.m:854) — every key now has a declared default; KeyExistsQ-style
    probing of the store is forbidden (and impossible through CFG).
F-f Auto-detected variables silently replacing user config
    (MatrixLoading.m:88).  E12 instead.
F-g Config writes mutating global symbol bindings: the old store assigned
    ``DiffExp`Symbols`x = FEC[LineParameter]`` at load and on update
    (State.m:154, DiffExp.m:138), rebinding a shared symbol as a side effect.
    DiffExp2 chart coordinates are module-local; Config stores the
    LineParameter symbol as DATA only.
F-h Validating after merging (partial updates surviving an abort; old order
    DiffExp.m:105 vs 121ff).  Validate-then-commit only.
F-i Filtering the user-visible config view (old CurrentConfiguration hid live
    keys, DiffExp.m:86-92).  CurrentConfiguration returns everything.
F-j Quietly accepting both-symbol-and-string twins of one key as DIFFERENT
    entries (possible in the old store: `System`AccuracyGoal` vs
    "AccuracyGoalValidate" mixed key types, State.m:107-108).  Canonicalization
    makes the twin collapse, and invariant 1 (exact key set) makes a stray
    entry impossible.

## 6. ABSORBED OLD CODE (and the numerical/operational lessons preserved)

- DiffExp/State.m:104-144 (DefaultConfiguration, DiffExpConfiguration, FEC)
  -> the schema + store + CFG.  Lesson: one mutable Association with mixed key
  types and no read validation is the root enabler of disease D1 in config
  space.
- DiffExp/State.m:146-163 (the `*Val` accessor layer: ChopPrecisionVal,
  LinearSolveChopPrecisionVal, FEWorkingPrecision, FEAccuracyGoal,
  DeltaPrescriptionsVal, UseMobiusVal, RadiusOfConvergenceVal,
  DivisionOrderVal, MatrixDirectoryVal, EpsilonOrderVal, LineParameterVal,
  ExternalScalesVal) -> all become `CFG["..."]` reads; the accessor layer is
  NOT reproduced (it added a second, also-unvalidated spelling for every key).
  ExportDisposition.md must map each name to its CFG key for the FT shim
  (State reach-ins x27 are mostly FEC[RationalizationTolerance]/FEC[Variables],
  per the execution review finding 3).
- DiffExp/DiffExp.m:86-197 (CurrentConfiguration/LoadConfiguration/
  UpdateConfiguration) -> sections 2/5 above.  Specific lessons preserved:
  ChopPrecision sync (121-130) = invariant 2 + Automatic resolution; LogFile
  (107-114) = E13; LineParameter checks (131-138) = E5/F-g; DeltaPrescriptions
  parsing + irreducibility (140-168) = 3.3/E6/E7; EstimateError parsing
  (150-158) = E3; the auto-LoadMatrices + cache-wipe block (171-190) = MOVED
  to API.m LoadSystem (Config is side-effect-free beyond the tolerance
  install); CrosscheckLevel/Flags resolution (115-120) = kept keys with the
  always-on carve-out.
- DiffExp/State.m:196-207 (CrosscheckFlags registry incl. "SingularityCheck"
  default-off) -> "CrosscheckLevel"/"CrosscheckFlags" keys; the
  SingularityCheck OFF-by-default footgun (an unprescribed branch point could
  pass silently unless the flag was enabled; Transport.m:480-492) is closed by
  making prescription checking structural and always-on (RewritePlan 3.1).
- FeynmanTrick/DiffExpIntegration.m:339-361,409-414,465-469 (the FT layer's
  fully-qualified config block, its UseMobius->False pin, ChopPrecision =
  precision-50, DivisionOrder 4, AbortOnAnalyticContinuationFail->False,
  re-adding DeltaPrescriptions after every LoadConfiguration because reset
  wiped them) -> after M5 the FT layer issues ONE
  `LoadConfiguration[...]`-equivalent through the DiffExp2 API with string
  keys; the "CRITICAL: Re-add delta prescriptions after LoadConfiguration
  reset" dance (465) disappears because LoadSystem, not LoadConfiguration,
  owns matrix-coupled state.
- FeynmanTrick/FIREInterface.m:203-226 (extractVariables/buildFIRESubstitution
  Global`-pinning) -> PinnedVariable.  Lesson: every external tool boundary
  (FIRE/Fermat, matrix files, user notebooks) speaks bare names; the ONLY safe
  internal representation is one pinned context plus SymbolName comparison.
- DiffExp/MatrixLoading.m:27-55,88,118 (variable auto-detection from `d*_*.m`
  filenames; the `eps`/`\[Epsilon]` acceptance at 118) -> detection itself is
  API.m/LoadSystem work, but the adopt-or-error rule (E12) and pinning are
  Config's.  Both `eps` and `\[Epsilon]` remain accepted (RewritePlan 3.2
  API.m bullet).

## 7. DEPENDENCIES

May call: ``DiffExp2`Tolerances` `` ONLY (derivation functions +
InstallToleranceState).  May be called by: EpsSeries, SectorSeries, Indicial,
Solve, Transport, Integrate, API, and (post-M5) the FT layer through the API
shim.  Config performs NO file I/O except "LogFile" (E13) and NO kernel-state
mutation outside its own store + the Tolerances install (F-g).

## 8. UNIT TESTS (Tests/test_config.m; names binding)

C1  `test_defaults_complete`: after `LoadConfiguration[{}]`,
    `CurrentConfiguration[]` has exactly the 3.2 kept-key set; spot values:
    "WorkingPrecision" -> 500, "ChopPrecision" resolves to 250,
    "LinearSolveChopPrecision" resolves to 250, "DivisionOrder" -> 3,
    "EpsilonOrder" -> 4, "EstimateError" -> "Fast", "ExpansionOrder" -> 50,
    "RadiusOfConvergence" -> 1, "SegmentationStrategy" -> "Predivision",
    "UseMobius" -> False, "UsePade" -> False, "Verbosity" -> 1,
    "AbortOnAnalyticContinuationFail" -> True.
C2  `test_read_unknown_key_loud`: `CFG["RationalizationTollerance"]` (typo)
    aborts with E1 listing valid keys.
C3  `test_symbol_key_canonicalization`: `CFG[System`WorkingPrecision] === 500`
    and, after creating a fresh symbol `Foo`Bar`WorkingPrecision`,
    `CFG[Foo`Bar`WorkingPrecision] === 500` — the context trap is neutralized
    by name canonicalization (PURPOSE paragraph scenario).
C4  `test_symbol_key_unknown_context_diagnostic`: `CFG[Foo`Bar`NoSuchKey]`
    aborts and the message contains "Foo`Bar`" (E1 context diagnostic).
C5  `test_write_unknown_key_loud_atomic`:
    `UpdateConfiguration[{"WorkingPrecision" -> 300, "NoSuchKey" -> 1}]` aborts
    (E2) AND `CFG["WorkingPrecision"]` still returns the prior value
    (all-or-nothing).
C6  `test_chop_sync`: `UpdateConfiguration[{"WorkingPrecision" -> 300}]` =>
    ChopPrecision resolves 150, LinearSolveChopPrecision 150, and
    `Tolerances`Tol["MatchTol"] === 10^-150` (install propagated).  Explicit
    `{"ChopPrecision" -> 200, "LinearSolveChopPrecision" -> 180}` => Tol pair
    (10^-200, 10^-180).  (DiffExp.m:121-130 semantics.)
C7  `test_chop_ge_wp_loud`: `{"WorkingPrecision" -> 100, "ChopPrecision" -> 100}`
    aborts with the E4 text.
C8  `test_estimate_error_strict`: value `"Slow"` aborts (E3); old-style string
    `"False"` ALSO aborts (the old parser accepted it as a coerce-to-False,
    DiffExp.m:151; DiffExp2 requires the Boolean False).
C9  `test_delta_prescriptions_parse`: input
    `{m^2 - s - I*Global`δ, {s - 4 m^2, 1}}` normalizes to
    `{{m^2 - s, -1}, {s - 4 m^2, 1}}` (δ-coefficient/I rule, DiffExp.m:142)
    with all symbols pinned to Global`; the user list is retrievable and equal
    to the effective list before any LoadSystem.
C10 `test_delta_prescriptions_reducible_loud`: `{(s-1)^2 (s-2) + I*Global`δ}`
    aborts with the E6 text ("irreducible polynomials").
C11 `test_line_parameter_clash_loud`:
    `{"Variables" -> {Global`s}, "LineParameter" -> Global`s}` aborts (E5).
C12 `test_variable_pinning`: `UpdateConfiguration[{"Variables" ->
    {Foo`Bar`s, Global`m}}]` stores symbols whose Context is "Global`" and
    whose names are {"s","m"} (PinnedVariable applied).
C13 `test_dropped_keys_dedicated_errors`: each of "IntegrationStrategy",
    "UseRationalRecurrence", "InvWronskSolver", "HomogeneousSolve",
    "KeepMatrixExpansions", "Parallel", "IgnoreIndicialCheck" aborts on write
    with E9 text containing the word "removed" and the per-key reason; and
    "SegmentationStrategy" -> "Dynamic" aborts with E10.
C14 `test_accuracy_goal_cross_field`: `{"AccuracyGoalValidate" -> "Before"}`
    with "AccuracyGoal" left at "?" aborts (E8);
    `{"AccuracyGoal" -> 30, "AccuracyGoalValidate" -> "Before"}` succeeds.
C15 `test_load_resets_then_applies`: after `UpdateConfiguration[
    {"Verbosity" -> 3}]`, `LoadConfiguration[{"WorkingPrecision" -> 300}]`
    yields "Verbosity" -> 1 (reset) and "WorkingPrecision" -> 300.
C16 `test_rationalization_tolerance_automatic`: default resolves
    `CFG["RationalizationTolerance"] === Tolerances`SnapTol[500] === 10^-250`;
    explicit `{"RationalizationTolerance" -> 10^-30}` => Tol["SnapTol"] ===
    10^-30.  (Kills the F-d hard-coded 10^-40 class: the configured value is
    what Tolerances serves.)
C17 `test_cfg_before_load_loud`: fresh state, `CFG["WorkingPrecision"]`
    aborts (E11); `ConfiguredQ[] === False` before, True after
    `LoadConfiguration[{}]`.
C18 `test_exports_visible_cross_context`: every section-2 symbol evaluates
    fully qualified from a scratch context (wolfram-package-context-traps
    regression; mirrors Tolerances T13).
C19 `test_no_raw_store_reads_in_tree`: source grep over DiffExp2/*.m — the
    private store symbol name appears only in Config.m, and the regexes
    `FEC\[` / `DiffExpConfiguration` have zero hits outside Config.m and
    Legacy/.  (Static test, runs without a kernel solve.)
C20 `test_crosscheck_flag_unknown_loud`: `{"CrosscheckFlags" ->
    {"FrobeniusSolutions"}}` aborts naming the flag as not in the DiffExp2
    registry (old registry State.m:196-206 does not carry over names whose
    machinery died); `{"CrosscheckFlags" -> {"SingularityCheck"}}` aborts with
    the always-on supersession pointer (3.2 table).

## 9. LINE BUDGET

~150 lines (RewritePlan 3.2; the execution review finding 3 sized Config at
~150 precisely so MatrixLoading's 389 lines do NOT creep in here — matrix
loading is API.m's).  Estimate: schema table ~45 (one line per key), CFG +
canonicalization ~20, Update/Load (validate-then-commit + sync + tolerance
install) ~40, DeltaPrescriptions parsing ~15, errors/usage ~30.  If over
budget, cut in this order: (1) collapse E9 per-key reasons into one lookup
table line each; (2) drop the E1 nearest-key suggestion text (keep the key
list); (3) move "LogFile" handling to API.m (it is the only I/O here).  DO NOT
cut: canonicalization, all-or-nothing commit, cross-field checks, dedicated
dropped-key errors, or any ::usage declaration.

## 10. OPEN QUESTIONS

Q1 Pinned context for variables is `Global`` (matches FIRE and user
   notebooks).  If DiffExp2 is ever loaded in an environment that sandboxes
   Global` (cloud kernels), PinnedVariable needs a configurable target
   context; out of scope v1 — flag if M5 hits it.
Q2 Should "CrosscheckLevel"/"CrosscheckFlags" survive at all, given always-on
   invariants?  Kept as the on-switch for EXPENSIVE optional checks (full
   ODE-residual sweeps vs spot checks).  The DiffExp2 check registry (names ->
   level) is owned by the modules that implement the checks; Config validates
   names against a registry constant that those modules register at load.  If
   that registration mechanism proves > ~10 lines, simplify to a fixed name
   list in Config.m and accept the coupling.
Q3 "AccuracyGoal" accepts "?" for backward compatibility.  Should DiffExp2
   instead require None?  Kept "?" because Reference/Examples and Tests pass
   it (old default State.m:107); revisit at M6 docs rewrite.
Q4 The FT layer's own `$FTConfig` (FeynmanTrick/FeynmanTrick.m:26-56,
   `SetFTOption` is unvalidated assignment) is OUT OF SCOPE for this module,
   but it has the same disease (string keys, no read validation,
   `Lookup[..., default]` fallbacks).  M6 hardening should either route it
   through a second Config schema instance or document why not — tracked for
   the M6 cleanup gate, not for M1.
Q5 LoadConfiguration resets DeltaPrescriptions (by design, full reset).  The
   old FT workaround re-added them after each reset (DiffExpIntegration.m:
   465-469).  After M5, does any caller still need a "preserve prescriptions
   across reset" mode?  Expected no (single configuration per LoadSystem);
   if the shim audit finds one, add an explicit `LoadConfiguration[...,
   "Keep" -> {"DeltaPrescriptions"}]` rather than reintroducing the dance.
