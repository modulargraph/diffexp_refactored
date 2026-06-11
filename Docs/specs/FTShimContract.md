# FT-layer reach-in shim contract (M0 task 3)

Status: M0 deliverable per Docs/RewritePlan.md section 6 (M0 task 3) and risk
R4.  This document is the binding contract for the M5 FT cutover: it
enumerates EVERY reference from the retained FT layer (FeynmanTrick/ and
Scripts/) into the old library's contexts, states what each reach-in
semantically NEEDS (as opposed to the mechanism it uses today), and names the
DiffExp2 replacement — a public API call per RewritePlan 3.2 /
Docs/specs/API.md, a config key per Docs/specs/Config.md, or
"eliminate: <design reason>".  M5 may not close while any row here is
unresolved (RewritePlan M5 "shim audit per disposition table"; R4).

Method: `grep -rn "DiffExp\`" FeynmanTrick/ Scripts/` plus targeted greps for
`Get[`/`Needs[` and bare top-level calls, every hit read in context, every
old-library target symbol located in DiffExp/ source.  Cross-checked against
Docs/reviews/survey_feynmantrick_refactor.md (items 1, 15, 16, 17, 24) and
the execution-review finding 3 / legacy-review finding 13 in
Docs/reviews/rewrite_plan_review_3lens.json.

---

## 0. Inventory totals (exact, vs RewritePlan's estimates)

| Old context | RewritePlan estimate | Found in FeynmanTrick/*.m | Found in Scripts/*.m | Total |
|---|---|---|---|---|
| ``DiffExp`State` `` | ~27 | 26 (all in DiffExpIntegration.m) | 2 (run_ft_stepwise.m:43,45) | 28 |
| ``DiffExp`Symbols` `` | ~10 | 10 (DiffExpIntegration.m) | 7 occurrences on 5 lines (check_transport_ode_residual.m:29,56,58,85,86) | 17 |
| ``DiffExp`Utilities` `` | ~5 | 5 (DiffExpIntegration.m) | 1 (check_transport_ode_residual.m:35) | 6 |
| ``DiffExp`Transport` `` | ~2 | 2 (DiffExpIntegration.m:440,474) | 0 | 2 |
| ``DiffExp`RegularizedIntegration` `` | ~2 | 2 (DiffExpIntegration.m:875,1015) | 6, of which 5 into `Private`` (check_transport_ode_residual.m:36,46; eval_dump_generic.m:14,21,34,37) | 8 |
| ``DiffExp`SingularityDecomposition` `` | ~1 | 1 (DiffExpIntegration.m:822) | 0 | 1 |
| ``DiffExp`Pade` `` | — | 0 | 0 | 0 |
| bare ``DiffExp` `` top level | "LoadConfiguration, UpdateConfiguration, TransportTo, IntegrateSystem etc." | 4 (LoadConfiguration x2: 367,460; UpdateConfiguration x2: 411,466) | 0 | 4 |
| `Get`/`Needs` of DiffExp files | — | 1 (DiffExpIntegration.m:331, lazy) | 3 (run_ft_stepwise.m:8, check_transport_ode_residual.m:13, eval_dump_generic.m:7) | 4 |

Notes on the estimates: the "x27" for State undercounts by one because line
DiffExpIntegration.m:211-213 is a single statement containing TWO
`MatricesIrreducibleFactors` occurrences, and it predates the two
script-side debug-flag pokes.  There are NO calls to ``DiffExp`TransportTo``
at top level (both transport calls are fully qualified
``DiffExp`Transport`TransportTo``) and NO calls to `IntegrateSystem` anywhere
in FeynmanTrick/ or Scripts/ (grep-verified) — the RewritePlan's "etc." list
item is empty for IntegrateSystem; its public-name disposition is handled in
Docs/specs/API.md 2.3 ("IntegrateSystem as a public name is DROPPED").
``DiffExp`Pade`` has zero reach-ins; the only Pade-related FT line is the
config write ``DiffExp`State`UsePade -> False`` (DiffExpIntegration.m:353),
counted under State.  Non-code mention:
FeynmanTrick/IMPLEMENTATION_STATUS.md:131 cites
``DiffExp`State`ChopPrecision`` in prose; the file is stale (survey item 25),
no shim action — fold/mark superseded at M6.

Abbreviations below: DEI = FeynmanTrick/DiffExpIntegration.m,
REG = DiffExp/RegularizedIntegration.m.  "CFG[...]"/key names refer to
Docs/specs/Config.md section 2/3.2; API entry points refer to
Docs/specs/API.md section 2.  All DiffExp2 public symbols live in
``DiffExp2` `` (RewritePlan section 1).

---

## 1. THE KNOWN CAMPAIGN TRAP, VERIFIED: bare-symbol FEC keys resolving in the wrong context

This is the FEC-key-context latent bug found by the FT survey
(Docs/reviews/survey_feynmantrick_refactor.md item 1) — same disease class as
the NormalizeLogPower unexported-symbol no-op
(memory/wolfram-package-context-traps.md; commit ebc4724).

Mechanism (verified in source):

- DEI declares `BeginPackage["FeynmanTrick`DiffExpIntegration`",
  {"FeynmanTrick`"}]` (DEI:7) and `Begin["`Private`"]` (DEI:44).  Its context
  path during parsing therefore does NOT include ``DiffExp`State` ``.
- The old config store is a plain Association whose tolerance entry is keyed
  by the SYMBOL ``DiffExp`State`RationalizationTolerance``
  (DiffExp/State.m:123, default 10^-10; ::usage at State.m:22), read through
  the unvalidated alias `FEC := DiffExpConfiguration` (DiffExp/State.m:144).
- Every DEI lookup spelled ``DiffExp`State`FEC[RationalizationTolerance]``
  qualifies FEC but NOT the key, so the key parses as the fresh symbol
  ``FeynmanTrick`DiffExpIntegration`Private`RationalizationTolerance`` and
  the lookup silently returns `Missing["KeyAbsent", ...]`.
- DEI's own author comment (DEI:339-341, "All config keys must be fully
  qualified to match DiffExp's internal symbols, since this package's
  BeginPackage restricts the context path") documents the trap for the WRITE
  side (DEI:346-361 writes ARE fully qualified) — the READ sites were missed.

The six sites and the dormant behavior each one silently disables:

| Site | What the Missing tol does today |
|---|---|
| DEI:98 (`realNumericAtActivePrecision`) | `TrueQ[Abs[Im[val]] < Missing[...]]` is False, so the imaginary-part chop NEVER fires: nominally-real endpoint values keep spurious tiny imaginary parts. |
| DEI:109 (`thetaRulesAtLocalPoint`) | the numeric `Abs[...] < tol` and `... < 0` comparisons go unevaluated → `TrueQ` False → `sign` falls through to the default `True, 1` branch: theta-branch (iδ side) detection for NEGATIVE local endpoints is effectively disabled; everything gets the θ+ rules. |
| DEI:127 (`localEndpointDirection`) | the numeric near-equality filter on segment bounds is dead; only exact `PossibleZeroQ` dedup runs, so numerically coincident bounds are not recognized and the inferred approach direction can be wrong/`Automatic`. |
| DEI:154 (`snapMainExpression`) | `Quiet[Check[...]]` catches nothing (Missing is not an error) but the `!NumericQ[tol]` guard at DEI:155 self-heals to the HARDCODED 10^-40 — which is NOT the documented default 10^-10 (State.m:123), so the user-configured value is ignored either way. |
| DEI:182 (`snapValuesFromFactors`) | same self-heal to hardcoded 10^-40 (DEI:183); configured tolerance ignored. |
| DEI:924 (`EvaluateLimitFromTransport`, `zeroQ`) | `zeroQ`'s numeric branch (DEI:925-926) is dead; endpoint-at-zero detection degrades to exact `PossibleZeroQ` only, so a numerically-zero-but-inexact local endpoint takes the "nonsingular endpoint inside the segment" branch (DEI:1028-1066) and is EVALUATED AT a singular endpoint instead of sector-limited. |

A seventh, related-but-different site: DEI:78 reads
``DiffExp`State`FEWorkingPrecision`` — that one is a correctly-qualified
EXPORTED accessor (`FEWorkingPrecision := FEC[System`WorkingPrecision]`,
DiffExp/State.m:159, ::usage State.m:43) and works; its defect is only the
`Quiet[Check[..., $FTConfig fallback, ... 500]]` chain (DEI:75-86), an A1
silent-fallback item.

DiffExp2 closure (constructive, per Docs/specs/Config.md): the validated
accessor `CFG[key]` canonicalizes Symbol keys by `SymbolName` — a
mis-contexted symbol twin reaches the right entry, and an unknown NAME is a
loud error printing the symbol's full context (Config.md E1, tests C2-C4).
The class is impossible by construction; the shim work in section 6 below
decides, per site, whether the read survives (→ `CFG`) or the surrounding
code is deleted outright.

---

## 2. Reach-in catalog: ``DiffExp`State` `` (28 sites)

### 2.1 Configuration writes — the TransportLevel config block

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:347 `MatrixDirectory -> matrixDir` | points the old loader at the level's exported slice dir | load THIS level's DE system | `LoadSystem[matrixDir]` (API.md 2.1); the directory now contains the exact full format `d<var>_full.m` (RewritePlan I1; M5 switches MatrixExport — see checklist C1).  `LoadConfiguration[..., "MatrixDirectory" -> dir]` also still triggers LoadSystem (API.md 2.8) if the one-call shape is preferred |
| DEI:348 `System`WorkingPrecision -> precision` | master precision | master precision | `"WorkingPrecision"` key (Config.md 3.2) |
| DEI:349 `ChopPrecision -> precision - 50` | chop threshold synced to WP (IMPLEMENTATION_STATUS.md:131 records the bug this fixed: stale absolute chop killing real signal at high WP) | a chop threshold that scales with WP | `"ChopPrecision"` key; either keep the explicit `precision - 50` (bit-compatible with the campaign) or adopt `Automatic` (resolves to wp/2 per Config.md 3.2).  DECISION: keep explicit `precision - 50` at cutover, revisit at M6 (risk note in checklist C3) |
| DEI:350 `ExpansionOrder -> expOrder` | t-truncation order | same | `"ExpansionOrder"` key |
| DEI:351 `EpsilonOrder -> epsOrder` | requested eps depth | same, with honest completeness | `"EpsilonOrder"` key; EpsWindow honesty replaces the FT-side trimming/warning plumbing (DEI:642-667, 1759-1777) |
| DEI:352 `UseMobius -> False` ("Required for integration!") | affine charts so the integrator's closed forms apply | guarantee no Mobius charts reach Integrate | `"UseMobius"` key, default ALREADY False (Config.md 3.2); belt-and-braces: Integrate.m/IntegrateOverLine REJECTS Mobius charts loudly (RewritePlan 3.2 Integrate.m; API.md 2.6), so the pin is enforced even if a user flips the key |
| DEI:353 `UsePade -> False` | avoid Pade interfering with saved raw series | same | `"UsePade"` key, default False.  This is the only Pade-adjacent FT line; no ``DiffExp`Pade`` reach-in exists |
| DEI:354 `DivisionOrder -> divisionOrder` (FT default 4, DEI:272) | chart sizing safety (GetCPL/GetCPR lesson) | same | `"DivisionOrder"` key (Config default 3; FT keeps passing 4 — ledger seed, RewritePlan section 5) |
| DEI:355 `Verbosity -> verbosity` | print gate | same | `"Verbosity"` key |
| DEI:356 `SegmentationStrategy -> "Predivision"` | predivision segmentation | same | `"SegmentationStrategy"` key — "Predivision" is the ONLY v1 value (Config.md E10), so the write becomes a no-op default; keep or drop |
| DEI:357 `UseRationalRecurrence -> useRationalRecurrence` | select the recursive solver fast path | a solver that is fast AND never silently falls back | ELIMINATE: dropped key (Config.md 3.2 dropped-key table, E9).  The denominator-cleared recursion is the only Solve.m path (RewritePlan I2).  Remove the TransportLevel option plumbing (DEI:278,302,357) and the run_ft_stepwise.m:191 pass-through |
| DEI:358 `IntegrationStrategy -> integrationStrategy` | select strategy stack member | nothing — strategy stack deleted | ELIMINATE: dropped key (Config.md E9); remove option plumbing DEI:279,303,358 |
| DEI:359 `"EstimateError" -> estimateError` (string key, correct by accident) | two-point error probe mode | same | `"EstimateError"` key (kept, Config.md 3.2; values strictly False/True/"Fast", E3) |
| DEI:360 `"HomogeneousSolve" -> homogeneousSolve` ("DontExpand", DEI:276) | strategy-stack switch | nothing — one homogeneous path | ELIMINATE: dropped key (Config.md E9 names this exact FT write); remove option plumbing DEI:276,301,360 |

### 2.2 Configuration reads

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:78 `FEWorkingPrecision` (in `activeNumericPrecision`, DEI:75-86) | active numeric precision for N[]/SetPrecision, with Quiet/Check fallback to $FTConfig then 500 | the loaded system's working precision, loudly | `CFG["WorkingPrecision"]` (Config.md 2); DELETE the Quiet/Check/500 chain (A1; Config.md F-d class).  Most callers of `numericAtActivePrecision` die with the manual combine/limit code (sections 4-5); the helper survives only for the FT seam (boundary printing, snap dedup) |
| DEI:98 `FEC[RationalizationTolerance]` (BUG, sec. 1) | im-part chop tolerance for "this value is really real" | a snap tolerance derived from WP | `CFG["RationalizationTolerance"]` (Config.md 3.2: REDEFINED as the FT-seam snap/dedup tolerance, Automatic -> `Tolerances`SnapTol[wp]`).  NOTE: the consuming helper survives only where section 6 keeps it; the EndpointLimit path that needed it most is replaced wholesale (sec. 5.2) |
| DEI:109 `FEC[RationalizationTolerance]` (BUG) | endpoint-sign tolerance for θ-branch choice | which SIDE of the expansion point is being approached | ELIMINATE with the theta machinery: DiffExp2 evaluation takes `"Direction" -> +1\|-1` explicitly (API.md 2.4/2.5) and the branch is resolved from the chart's "Prescriptions" record (RewritePlan 3.1); no tolerance-based sign inference remains |
| DEI:127 `FEC[RationalizationTolerance]` (BUG) | dedup of segment bounds to infer approach direction | the approach direction at a chart endpoint | ELIMINATE: SegmentRecord carries exact "LocalInterval"/"MainInterval" (API.md 3.7) and EndpointLimit takes "Direction"; no numeric inference |
| DEI:154 `Quiet[Check[FEC[RationalizationTolerance], 10^-40]]` (BUG + hardcoded fallback) | snap tolerance for chart-map endpoint snapping | exact chart maps | ELIMINATE `snapMainExpression` (DEI:151-178) entirely: DiffExp2 charts have EXACT centers and affine maps from the exact singularity solve (RewritePlan 3.2 Transport.m; LocalSolution "Center" -> exact x0), so there is nothing to snap.  Forbidden-fallback citation: Config.md F-d names this exact line |
| DEI:182 `Quiet[Check[FEC[RationalizationTolerance], 10^-40]]` (BUG + hardcoded fallback) | dedup tolerance for numeric roots of IBP singular factors in [0,1] (`snapValuesFromFactors`, DEI:180-202) | dedup of candidate singular-factor roots | KEEP the seam, fix the read: `CFG["RationalizationTolerance"]`, delete Quiet/Check (Config.md C16 pins this).  The function itself feeds the ExtraSingularFactors mechanism (sec. 2.3); under DiffExp2 root solving is exact so the dedup shrinks to exact `Union` plus one numeric in-[0,1] filter |
| DEI:328 `If[!ValueQ[DiffExp`State`FEC], ...]` | "is the old library loaded yet?" lazy-load guard | DiffExp2 loaded exactly once | ELIMINATE: FeynmanTrick.m loads DiffExp2 at package load (checklist C2); `ConfiguredQ[]`/the exactly-one-system invariant (API.md sec. 4 inv. 1) covers the "not configured yet" error path loudly |
| DEI:372 `First[DiffExp`State`FEC[System`Variables]]` | recover the variable symbol the old loader CREATED from filenames, to align FT's points/factors with the loader's context (`detectedVar`, plus the SymbolName remap dance DEI:374-393, 417) | one agreed-on variable symbol across FIRE, matrix files, and the solver | `LoadSystem[...]` returns `SystemInfo["Variables"]` (API.md 3.4) AND Config pins every variable to ``Global` `` via `PinnedVariable` (Config.md 2, sec. 6 "FIREInterface Global`-pinning" lesson) — two symbols agreeing in name but differing in context can no longer coexist, so the entire detectedVar/SymbolName-remap dance (DEI:371-393, 605-631, 977-991) is deleted, not ported |
| DEI:418 `DiffExp`State`ExternalScalesVal` (verbosity print; accessor at DiffExp/State.m:150) | print the loaded variables | same | `SystemInfo["Variables"]` from the LoadSystem return (or `CurrentConfiguration[]["Variables"]`) |
| DEI:924 `FEC[RationalizationTolerance]` (BUG) | `zeroQ` numeric-zero test for "is this local endpoint the chart center?" | reliable endpoint-at-center classification | ELIMINATE with `EvaluateLimitFromTransport` (sec. 5.2): EndpointLimit consumes the typed LocalSolution/TransportResult — whether the endpoint is the chart center is a property of the SegmentRecord, not a numeric guess |

### 2.3 Mutable solver-state writes — the segmentation-alphabet extension

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:211-213 (`appendMatrixFactors`, 2 occurrences) | appends IBP-coefficient denominator factors to ``DiffExp`State`MatricesIrreducibleFactors`` (declared State.m:74; written by the loader at MatrixLoading.m:217; consumed by segmentation/transport at LineSegmentation.m:71, Transport.m:51,572-573) so transport segments at the roots of c_j(x) | the segmentation alphabet must include singular points of the IBP coefficients, not just of the matrices, so that (a) chart boundaries land on them and (b) prescriptions can be attached | NEW NAMED INPUT, not a state poke: `TransportTo`/`IntegrateOverLine` option `"ExtraSingularFactors" -> {poly..}` (exact polynomials in the pinned variables) unioned EXACTLY into the segmentation alphabet for that call.  ** API GAP: Docs/specs/API.md 2.3 does not yet list this option — it must be added; flagged as an open question of this contract.**  Raw writes to solver state are forbidden in DiffExp2 (no exported mutable state) |
| DEI:402 (read of `MatricesIrreducibleFactors` inside the `PrescribeMatrixFactors` branch) | unions matrix factors into the delta-prescription list | prescriptions for every factor the transport may cross | the matrix factors are `SystemInfo["IrreducibleFactors"]`/`"IrreducibleFactorsExact"` (API.md 3.4); the FT builds its `{poly, sign}` list from that return value + its extra factors and passes it as the `"DeltaPrescriptions"` config key (Config.md 3.3) |

### 2.4 Debug-flag pokes (Scripts)

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| Scripts/run_ft_stepwise.m:43 ``DiffExp`State`$DebugFuchsianizedCheck = True`` (env DEBUG_FUCHS_CHECK) | enables per-piece solve-vs-block-ODE verification inside the old resonant-recurrence strategy (consumed at DiffExp/IntegrationStrategies/ResonantRecurrence.m:1226,1380,1502,1522; never declared in State.m — created by assignment) | confidence that every local solve satisfies its ODE | ELIMINATE: the strategy stack is deleted (RewritePlan I2) and DiffExp2 has an ALWAYS-ON ODE residual spot-check per chart (RewritePlan 3.1 invariants); expensive full sweeps hang off `"CrosscheckLevel"`/`"CrosscheckFlags"` (Config.md 3.2).  Keep the env hook working against Legacy/ for oracle replays until M6 |
| Scripts/run_ft_stepwise.m:45 ``DiffExp`State`$DebugBlockResidualSeries = True`` (env DEBUG_BLOCK_RESID; consumed ResonantRecurrence.m:956,1339,1438) | prints escaping residual series | same | same as above |

---

## 3. Reach-in catalog: bare ``DiffExp` `` top-level calls (4 sites)

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:367 `DiffExp`LoadConfiguration[diffExpConfig]` (old def DiffExp/DiffExp.m:80,95-99; side effect: triggers LoadMatrices, DiffExp.m:171-179) | reset + apply config AND load matrices for the level | configure, then load the system | `DiffExp2`LoadConfiguration[stringKeyedRules]` (Config.m semantics: reset-then-apply, validate-then-commit) + `LoadSystem[matrixDir]`.  Loading is an explicit verb (Config.md sec. 1); the compat path `LoadConfiguration[..., "MatrixDirectory" -> dir]` forwarding to LoadSystem also exists (API.md 2.8) |
| DEI:411 `DiffExp`UpdateConfiguration[{DeltaPrescriptions -> ..., "AbortOnAnalyticContinuationFail" -> False}]` (old def DiffExp.m:81,102-104) | attach prescriptions AFTER load (the loader reset wiped them) + soften continuation failures | prescriptions attached once; FT-pipeline failure mode | `UpdateConfiguration[{"DeltaPrescriptions" -> ..., "AbortOnAnalyticContinuationFail" -> False}]` — both keys kept (Config.md 3.2).  The ORDER dance disappears: LoadSystem, not LoadConfiguration, owns matrix-coupled state (Config.md sec. 6, Q5) |
| DEI:460 second `DiffExp`LoadConfiguration[diffExpConfig]` ("Reload config for transport in other direction") | reset solver state between the two transport directions | two independent transports from one anchor | ELIMINATE: DiffExp2 TransportTo is stateless with respect to other calls (bidirectional, RewritePlan 3.2 API.m; API.md 2.3); call TransportTo twice, no reload.  The re-`appendMatrixFactors` at DEI:461-463 dies with sec. 2.3 |
| DEI:466 second `DiffExp`UpdateConfiguration[...]` ("CRITICAL: Re-add delta prescriptions after LoadConfiguration reset") | re-attach prescriptions wiped by the reload | prescriptions to persist for the session | ELIMINATE with the DEI:460 reload (Config.md sec. 6 names this exact dance; Q5 records the contingency if a keep-across-reset mode is ever truly needed) |

---

## 4. Reach-in catalog: ``DiffExp`Symbols` `` (17 sites)

The old library publishes shared mutable symbols: `x` (the internal line
parameter, DiffExp/Symbols.m:17, REBOUND by config at State.m:154), `Logx`
(Symbols.m:12), and the analytic-continuation markers `θp`/`θm`
(Symbols.m:13-14, introduced into series by the crossing rules at
DiffExp/AnalyticContinuation.m:73-101).  DiffExp2 has NO shared symbol
surface: chart coordinates are module-local data (Config.md F-g), log powers
are the exact sector tag `p`, and branch choice is the explicit "Direction"
option resolved through the chart's "Prescriptions" record (crossing applied
at transport time as phase x unipotent log-chain mixing, RewritePlan 3.2
Transport.m).  Consequently EVERY Symbols reach-in is "eliminate"; the table
records what each site was for so the M5 rewrite can verify nothing semantic
is lost.

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:64-65, 69-70 (`$thetaPlusRules`/`$thetaMinusRules`: `θp -> 1, θm -> 0` and reverse) | resolve saved-series branch markers to a side of the cut | evaluate a local solution on a chosen side | ELIMINATE: "Direction" option of EndpointLimit/SolveAtPoint (API.md 2.4/2.5); branch baked in at evaluation by SectorSeries' branch rule.  (`$thetaRules` at DEI:73 is already dead — never read; survey item 11) |
| DEI:143, 153, 700, 928 (`xLocal = DiffExp`Symbols`x`) | name the local chart coordinate of saved positional SegmentData series for substitution/series ops | manipulate transported local series | ELIMINATE: FT stops touching raw series — the manual combine loop (DEI:681-794), endpoint evaluation (DEI:140-149, 1028-1066), and segment-map handling consume the typed SegmentRecord/LocalSolution through IntegrateOverLine / EndpointLimit / ToPiecewise (API.md 2.5-2.7, compat row 18 "BREAK -> named SegmentRecord") |
| DEI:146 (`Logx -> Log[xLocal]` substitution at final evaluation) | turn the symbolic log marker into a number away from x=0 | evaluate log sectors at a point | ELIMINATE: SectorSeries evaluation handles (eps Log t)^p/p! tags exactly (RewritePlan 3.1) |
| DEI:831 (`SeriesCoefficient[orders[[n]], {DiffExp`Symbols`x, 0, minPow}]`, Verbosity>=3 debug block) | inspect leading coefficients of combined series | diagnose negative leading powers | ELIMINATE with the debug block (DEI:807-842): LocalSolution carries exact sector tags and honest windows; the diagnostic becomes a one-line print of `"Sectors"` |
| check_transport_ode_residual.m:29 (`xLocal`), 56,58 (`Logx -> LL`), 85-86 (θ branch rules x4 occurrences) | rebuild Normal forms of saved series to check the ODE residual of a saved transport segment | independent verification that transported series solve the ODE | ELIMINATE/RETARGET: DiffExp2 runs an ODE residual spot-check per chart always-on (RewritePlan 3.1); the standalone checker is retired or rewritten against SegmentRecord + a public chart-evaluation call at M6.  Until then the script remains a Legacy/-oracle tool (RewritePlan R8) |

---

## 5. Reach-in catalog: ``DiffExp`Utilities` ``, ``DiffExp`Transport` ``, ``DiffExp`RegularizedIntegration` ``, ``DiffExp`SingularityDecomposition` ``

### 5.1 Utilities (6 sites — all `PChop`)

`PChop` chops below ChopPrecision (def DiffExp/Utilities.m:103, ::usage :28;
note its body reads `FEC[ChopPrecision]` with the key resolving correctly
only because Utilities.m compiles inside the DiffExp package — the same trap
one context away).

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:177 (in `snapMainExpression`) | chop the snapped chart map | nothing — exact charts | ELIMINATE with snapMainExpression (sec. 2.2, DEI:154 row) |
| DEI:705 (chop `xMainExpr`) | clean the segment's local→main map | exact "LineMap" | ELIMINATE: SegmentRecord "LineMap" is exact (API.md 3.7) |
| DEI:733, 777, 783 (chop IBP-coefficient series and combined totals inside the manual combine loop) | numerical hygiene of the hand-rolled combination | exact object-level combination | ELIMINATE with the combine loop: `IntegrateOverLine[..., "Combination" -> ibpRow, ...]` forms Σ c_j f_j at object level with exact cancellation (API.md 2.6; RewritePlan 3.3 "integrate").  Chopping inside DiffExp2 is owned by Tolerances.m (`chopFloor`), never by the caller |
| check_transport_ode_residual.m:35 | chop the reconstructed map in the checker script | — | ELIMINATE/RETARGET with the script (sec. 4 last row) |

### 5.2 Transport (2 sites) and RegularizedIntegration (8 sites)

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:440-446, 474-480 ``DiffExp`Transport`TransportTo[{startPoint, bcs}, <\|var -> endpoint\|>, 1, True]`` (old def Transport.m:514, ::usage :25) | transport the level's masters from the anchor 11/23 to each endpoint of [0,1], saving expansions; endpoints are SINGULAR (series returned, not evaluated) | bidirectional transport to singular endpoints with SaveExpansions | `DiffExp2`TransportTo` — SAME positional call shape (API.md 2.3 lists this exact FT pattern as compat row 9, KEEP); singular-endpoint mode now returns the typed LocalSolution + named SegmentRecords (API.md 2.3 SINGULAR-ENDPOINT MODE).  The two directions need no config reload between them (sec. 3) |
| DEI:875 ``DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[combinedData, {0,1}, prefactorSpec, combinedMinPower]`` (old def REG:2549, ::usage REG:76) | integrate the hand-combined single-master series over [0,1] with x^(v1-1)(1-x)^(v2-1) prefactor, Laurent-aware | ∫₀¹ x^(v1-1)(1-x)^(v2-1) Σ_j c_j(x,eps) f_j(x,eps) dx with exact Laurent bookkeeping across interior singular points | `IntegrateOverLine[transportResult, {0,1}, "Prefactor" -> <\|"PowerAtLower" -> v1-1, "PowerAtUpper" -> v2-1, "RationalFactor" -> 1\|>, "Combination" -> ibpCoeffRow]` (API.md 2.6).  This single call replaces ALL of IntegrateCombinedMasters' manual machinery: ExpandIBPCoeffLaurent windows (EpsSeries.m owns Laurent windows), the per-segment combine loop DEI:681-794, the `expOrd = 30` truncation cliff (DEI:757; survey item 4), the incomplete-top-order trim DEI:642-667 and its "keep at least one output order" floor (A1: deleted), and the epsMinPower side channel (inside LaurentValue, API.md 2.6) |
| DEI:1015 ``DiffExp`RegularizedIntegration`EvaluateEndpointLimitSectors[seriesAtMaster, direction]`` (old def REG:1965, ::usage REG:59) | per-sector endpoint limit with the b≠0 drop rule, per master, then hand-weighted by IBP coefficients with Quiet[Check[...,0]] chains (DEI:1037-1059, 1080-1093) | lim_{x->0 or 1} Σ_j c_j(x,eps) f_j(x,eps) with the dimreg drop rule applied per EXACT sector | `EndpointLimit[src, "Direction" -> ±1, "Combination" -> ibpCoeffRow]` (API.md 2.5): combination formed at object level BEFORE the limit; b≠0 sectors dropped exactly; divergence is a loud error.  The silent-zero Check fallbacks are FORBIDDEN in the replacement (API.md 2.5 names these lines) |
| eval_dump_generic.m:21 `DefiniteIntegralWithPrefactorLaurent[...]` (public) | replay a saved integral dump | reproduce one integral call in isolation (the dump-replay fast-iteration workflow) | RETARGET at M5/M6: Integrate.m's dump-on-failure artifact (RewritePlan section 6 KERNEL discipline, section 7) replaces the old `laurentIntegralDump`; the replay script calls `IntegrateOverLine` on the dumped TransportResult.  Until then the script replays Legacy dumps only |
| eval_dump_generic.m:14 ``REG`Private`laurentIntegralDump`` (assigned REG:247) | read the dump payload symbol | a STABLE dump file format | same as above — the new dump is a documented association, not a Private symbol (the reload-no-op hazard of Private reach-ins is the survey's item 24) |
| eval_dump_generic.m:34 ``REG`Private`segmentActualBounds`` (def REG:2366) | recover a segment's main-variable bounds | per-segment bounds | SegmentRecord "MainInterval" (API.md 3.7) |
| eval_dump_generic.m:37 ``REG`Private`IntegrateSegmentWithPrefactorLaurent`` (def REG:2840) | re-run one segment's integral | per-segment diagnostic integration | `IntegrateOverLine[result, segMainInterval, ...]` on the single segment's range (IntegrateOverLine accepts any [x0,x1] covered by SegmentData) |
| check_transport_ode_residual.m:36 ``REG`Private`segmentMainExpression`` (def REG:2355) | extract the local→main map from the positional 5-tuple | the chart map | SegmentRecord "LineMap"; script retargeted/retired (sec. 4) |
| check_transport_ode_residual.m:46 ``REG`Private`uncompressSeriesData`` (def REG:2413; DEI also carries its own copy at DEI:240-248) | uncompress saved series element 5 | readable saved chart data | ELIMINATE: compression is transparent inside DiffExp2 (ToPiecewise/IntegrateOverLine uncompress internally, API.md 2.7); the FT-side copy `uncompressSeriesData`/`prepareTransportSegments` (DEI:240-254, 643, 690, 995) is deleted with the manual series handling |

### 5.3 SingularityDecomposition (1 site)

| Site | What it does | What it NEEDS | DiffExp2 replacement |
|---|---|---|---|
| DEI:822 ``DiffExp`SingularityDecomposition`DecomposeSingularity[orders]`` (def DiffExp/SingularityDecomposition.m:305, ::usage :25 — note its usage text says it "uses RationalizationTolerance from configuration": it numerically FITS a and b from collapsed towers) | Verbosity>=3 diagnostic: print the {a,b} sectors of the combined segment series | knowledge of the sector spectrum at a chart | ELIMINATE: the LocalSolution carries `"Sectors" -> {<\|"a","b","p",...\|>..}` as EXACT data (RewritePlan I1/I2; API.md 2.4 explicitly: "DecomposeSingularity becomes unnecessary").  The execution review (finding 4) requires this call site in the R4 shim audit: it is hereby dispositioned — debug print of exact tags, no fit, no tolerance.  SectorSeries.m absorbs the partial-fraction role of SingularityDecomposition.m (RewritePlan 3.2) |

---

## 6. Get/Needs of DiffExp files, and FT-internal Private reach-ins (adjacent hygiene)

| Site | What it does | DiffExp2 replacement |
|---|---|---|
| DEI:322-337 (path build + `Quiet[Get[diffExpPath], {General::shdw, Symbol::shdw}]` inside `Block[{$ContextPath}]`, guarded by `ValueQ[DiffExp`State`FEC]`) | lazy-load the old library on first TransportLevel call, suppressing shadowing warnings | FeynmanTrick.m loads DiffExp2 ONCE at package load (a `Get`/`Needs` next to its own component loads, FeynmanTrick/FeynmanTrick.m:15-21); the shdw-Quiet and $ContextPath block die — DiffExp2 exports a clean, small symbol set (API.md sec. 2) with no Global` shadow collisions |
| Scripts/run_ft_stepwise.m:7-8, check_transport_ode_residual.m:12-13, eval_dump_generic.m:7 (`Get["DiffExp.m"]` etc.) | script bootstrap | swap to the DiffExp2 entry file at M5 cutover (run_ft_stepwise), at retarget time (other two); Legacy/ runs keep the old Get against the frozen tree (RewritePlan M6/R10 path sweep) |

FT-internal `Private`` reach-ins (NOT DiffExp reach-ins; listed for
completeness because the M5 edit touches the same lines — full treatment is
survey item 24): run_ft_stepwise.m:163
``FeynmanTrick`DiffExpIntegration`Private`RequiredTransportEpsilonOrder``;
check_transport_ode_residual.m:37 ``...Private`snapMainExpression`` (dies
with snapMainExpression); export_pysecdec_family_specs.m:171
``...Private`BoundaryRequestRecords`` and :104
``FeynmanTrick`Private`DimensionExpression``.  Export-with-::usage or
eliminate per the survey; `RequiredTransportEpsilonOrder` is additionally
superseded by the RewritePlan 3.4 static order budget at M5 (kept as a
lower-bound diagnostic only, per the run_ft_stepwise.m:154-161 comment and
the box_bubble 9/11 lesson).

---

## 7. API gaps this audit exposes (feed back into Docs/specs/API.md before M5)

G1 `"ExtraSingularFactors"` option (sec. 2.3): TransportTo/IntegrateOverLine
   must accept exact extra alphabet polynomials.  Without it the FT pipeline
   cannot place charts at IBP-coefficient poles and the "Combination"
   rational multiply is not closed on charts that contain such a pole in
   their interior.  Owner: API.md 2.3/2.6 + Transport.m spec.  (The
   alternative — silently widening the alphabet from the Combination's
   denominators — is rejected: implicit state derived from a later call's
   arguments is the D1 pattern again.  Explicit option, loud if a
   Combination denominator root is not in the alphabet.)
G2 Prescription signs for extra factors: old `deltaPrescriptionsForFactors`
   (DEI:224-238) auto-prescribed {var, sign}, {1-var, sign} and every extra
   factor with one global sign (default +1, DEI:280).  DiffExp2 keeps this an
   FT-layer responsibility (build the `{poly, sign}` list, pass via
   "DeltaPrescriptions"); the contract is that Transport.m's
   construction-time consistency check (RewritePlan 3.1 "Prescriptions")
   makes a WRONG sign loud, which the old pipeline never checked.  Pentagon's
   "unrecognized branch point" warning is exactly this class (M0 task 16).
G3 LaurentValue boundary ingestion: PrepareBoundaryConditions form (iii)
   (API.md 2.2) replaces the EpsPrefactors side channel
   (`ShiftRawBoundariesToFinite`, DEI:1258-1280; "EpsPrefactors" plumbing
   DEI:316-319, 514, 576-577, 918-922, 1495-1510; run_ft_stepwise.m:142-153).
   Already specced — listed here because five DEI functions and the stepwise
   driver change shape when it lands.

---

## 8. M5 cutover checklist (ordered; each step ends with battery + `FT_EXAMPLES=bubble` stepwise gate per the survey's batch discipline)

Pre-flight (blocking, from M0): Docs/ExportDisposition.md rows exist for
every symbol named in sections 2-5; G1 resolved in API.md/Transport spec;
pentagon prescription triage (M0 task 16) has produced the corrected
FTExamples DeltaPrescriptions so step 4 has a truth to validate against.

1. **MatrixExport -> full format; LoadSystem consumes it.**
   FeynmanTrickIteration `ExportLevel[..., "diffexp", order]` call sites
   (DEI:1779-1781, 1858-1861; run_ft_stepwise.m:173-175) switch to writing
   `d<var>_full.m` via ExportGeneralMatrix (FeynmanTrick/MatrixExport.m:91-114)
   instead of eps-sliced ExportDiffExpMatrix (MatrixExport.m:42-83); the
   per-level `transportEpsOrder` re-export dance (DEI:1731-1738) becomes
   unnecessary for matrices (exact files are order-independent) — KEEP the
   order bookkeeping for boundary windows only.
   RISK: low-medium.  Exact rational matrices are bigger and Together/Factor
   at export can be slow for pentagon/dbox-class systems; benchmark at first
   ladder rung.  Slice export retained only for Legacy parity runs.
2. **Load mechanics.**  FeynmanTrick.m loads DiffExp2; delete DEI:321-337
   (lazy-load guard + shdw Quiet).  RISK: low.  Watch
   Tests/test_package_loading.m, test_symbols_namespacing.m.
3. **Config block.**  Replace DEI:346-367 with
   `DiffExp2`LoadConfiguration[{"WorkingPrecision" -> precision,
   "ChopPrecision" -> precision - 50, "ExpansionOrder" -> expOrder,
   "EpsilonOrder" -> epsOrder, "UseMobius" -> False, "UsePade" -> False,
   "DivisionOrder" -> divisionOrder, "Verbosity" -> verbosity,
   "EstimateError" -> estimateError}]` + `LoadSystem[matrixDir]`; delete the
   `UseRationalRecurrence`/`IntegrationStrategy`/`HomogeneousSolve`/
   `SegmentationStrategy` writes and their TransportLevel option plumbing
   (DEI:276-279, 301-303, 356-360) and the run_ft_stepwise.m:191 pass-through.
   RISK: low-medium.  Config E9 errors will LOUDLY flush any forgotten
   dropped-key write — that is the design working, not a regression.  Keeping
   explicit `precision - 50` (not Automatic) preserves campaign numerics;
   note the choice in the ledger.
4. **Prescriptions + extra singular factors.**  Replace `appendMatrixFactors`
   (DEI:204-222, calls at 379-381, 461-463) and the
   MatricesIrreducibleFactors read (DEI:402) with: factors :=
   union(SystemInfo["IrreducibleFactors"], CollectLevelIBPSingularFactors
   output) passed as `"ExtraSingularFactors"` (G1) and
   `deltaPrescriptionsForFactors`-built `{poly, sign}` list passed through
   `UpdateConfiguration[{"DeltaPrescriptions" -> ...,
   "AbortOnAnalyticContinuationFail" -> False}]` ONCE (delete the second
   config round-trip DEI:459-469 entirely).  The detectedVar remap
   (DEI:371-393) collapses to nothing under PinnedVariable; `snapValues`
   bookkeeping survives only as the exact-root in-[0,1] filter feeding
   nothing (chart maps exact) — delete `snapMainExpression` (DEI:151-178)
   and the "SnapValues" result key (DEI:515, consumed at 707).
   RISK: medium.  This is where pentagon's configuration gap lives;
   Transport.m's prescription consistency check turns formerly-silent wrong
   signs into loud errors.  Triage protocol: prescription/config-class per
   the M5 gate taxonomy (RewritePlan section 6 M5).
5. **TransportLevel transport calls.**  DEI:440-446/474-480 →
   `DiffExp2`TransportTo` (same positional shape, compat row 9); result
   handling switches from positional 5-tuples to SegmentRecords; the
   Reverse-and-join of lower/upper halves (DEI:488-516) keeps its logic but
   operates on named records; drop the write-only "LowerResult"/"UpperResult"
   keys (survey item 11).  RISK: low-medium (shape change ripples into steps
   6-7; do not gate until 7 lands — steps 5-7 are one battery unit).
6. **IntegrateCombinedMasters -> IntegrateOverLine.**  Replace the body
   (DEI:562-892) with: build `ibpCoeffRow` (exact rational c_j(x, eps) —
   NO pre-expansion in eps; eps-window arithmetic is EpsSeries.m's), call
   `IntegrateOverLine[transportResult, {0,1}, "Prefactor" -> <\|...v1-1,
   v2-1, 1\|>, "Combination" -> ibpCoeffRow]`, scale by Gamma prefactor via
   LaurentValue arithmetic.  Deletes: ExpandIBPCoeffLaurent's truncation
   windows for this path, the combine loop (DEI:681-794) with its
   `expOrd = 30` cliff and PChop hygiene, the incomplete-top-order trim +
   "keep at least one output order" floor (DEI:642-667), the all-zero
   shortcut (DEI:847-864 — IntegrateOverLine of a zero combination is
   exactly zero), the DecomposeSingularity debug block (DEI:807-842).
   RISK: HIGH (the load-bearing change of the rewrite).  Gates: closed-form
   multisector pins (M4), box_bubble 9/11 budget unit, then the ladder.
7. **EvaluateLimitFromTransport -> EndpointLimit.**  Replace DEI:908-1115
   with `EndpointLimit[src, "Direction" -> dir, "Combination" ->
   ibpCoeffRow]` where src is the TransportResult of the matching direction
   (singular endpoint carries the LocalSolution).  Deletes: theta rule sets
   (DEI:63-73) and `thetaRulesAtLocalPoint`/`localEndpointDirection`
   (DEI:107-138), segment selection by numeric distance (DEI:934-967),
   `zeroQ` (DEI:924-926), every `Quiet[Check[..., 0]]` silent zero
   (DEI:1037-1059, 1080-1093).  RISK: medium-high.  The combine-before-limit
   semantics can change values where the old per-master limit silently
   zeroed a divergent coefficient — any difference is a FINDING (D1 bug
   surfacing), to be triaged against pySecDec references, not papered over.
8. **Boundary plumbing on LaurentValue (G3).**  ComputeLevelBoundary returns
   LaurentValues end-to-end; `ShiftRawBoundariesToFinite` (DEI:1258-1280)
   survives only as the adapter producing the legacy "BoundaryValues"+
   "EpsPrefactors" pair for printing; PrepareBoundaryConditions form (iii)
   consumes the raw LaurentValues for the next level's transport.  The five
   wrong-length-prefactor zero-fill sites (survey item 2) die with the side
   channel.  RISK: medium (touches every level boundary; STEPWISE JSON rows
   must stay bit-comparable via Scripts/compare_stepwise_log.py).
9. **Tolerance/precision read survivors.**  Whatever remains of
   `activeNumericPrecision`/`realNumericAtActivePrecision` after steps 4-8
   (expected: only output cleaning in run_ft_stepwise and root dedup in
   `snapValuesFromFactors`) reads `CFG["WorkingPrecision"]` /
   `CFG["RationalizationTolerance"]`; delete every Quiet/Check fallback.
   RISK: low AFTER steps 4-8 (the dormant-logic activation risk the survey
   rated medium applies to fixing the reads UNDER the old code paths; here
   those paths are already gone — that is why this step is sequenced late).
10. **Scripts.**  run_ft_stepwise.m: drop the two State debug-flag pokes
    (sec. 2.4) or fence them behind a Legacy-run branch; level loop adapts to
    steps 5-8 shapes; keep STEPWISE/FINAL JSON format EXACTLY (the
    comparator contract).  check_transport_ode_residual.m and
    eval_dump_generic.m: keep pointed at Legacy/ as oracle tools; retarget or
    retire at M6 (battery disposition).  export_pysecdec_family_specs.m: no
    DiffExp reach-ins; only the FT-internal Private exports (sec. 6).
    RISK: low.
11. **Audit close-out.**  Re-run the section-0 greps: zero hits for
    ``DiffExp`` ` in FeynmanTrick/*.m; Scripts hits only in the declared
    Legacy-oracle tools.  Tests C18/C19-style static checks (Config.md sec. 8)
    extended with this grep become the permanent regression gate.  Then the
    M5 ladder per RewritePlan section 6 (bubble → ... → pentagon →
    double_box_planar, >= 6/8 with mandatory triage).

Sequencing rationale: 1-3 are independent of the math and flush config-class
breakage first; 4 isolates the prescription/alphabet semantics (pentagon's
known gap) before any numerics change; 5-7 are the representation cutover and
gate as one unit; 8 changes shapes between levels; 9 is deliberately AFTER
the code it would have activated is deleted; 10-11 close.

## 9. Open questions

Q1 (= G1) "ExtraSingularFactors" option is not yet in API.md 2.3/2.6 — must
   be added before M5 step 4.  Decision needed from the API/Transport spec
   owners on whether it lives on TransportTo, IntegrateOverLine, or both.
Q2 Step 3 keeps explicit `ChopPrecision -> precision - 50` rather than
   Automatic (wp/2).  At wp=500 these are 450 vs 250 — materially different.
   Campaign numerics were validated at precision-50; revisit at M6 with a
   one-off A/B on the ladder.
Q3 `RequiredTransportEpsilonOrder` (FT Private, run_ft_stepwise.m:163) vs
   the RewritePlan 3.4 static budget: the contract assumes the 3.4 formula
   becomes the authority at M5 (validated by the box_bubble 9/11 unit) and
   the old function stays as a cross-check print only.  If the formula unit
   fails, the cutover order of step 8 must hold the old carry-all-depth
   behavior (run_ft_stepwise.m:154-161 comment) until triaged.
