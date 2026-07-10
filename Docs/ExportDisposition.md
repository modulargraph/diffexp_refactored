# Export disposition table (RewritePlan M0 deliverable, agent task 2)

Status: complete, 2026-06-11.  Scope frozen at the working tree corresponding to
RewritePlan baseline (old library frozen as parity oracle).

## Method and scope

ENUMERATION.  Every symbol with a `::usage` definition before
`Begin["`Private`"]` in `DiffExp/DiffExp.m` (the umbrella package) and every
file under `DiffExp/` was extracted.  Total: **229 exported symbols**.
Verification notes:

- The repo-root `DiffExp.m` is a 7-line path loader (no exports).
- The six files under `DiffExp/IntegrationStrategies/` (`Helpers.m`,
  `Default.m`, `VOP.m`, `Recurrence.m`, `ResonantRecurrence.m`, `Dispatch.m`)
  are loaded INSIDE `Begin["`Private`"]` of `DiffExp/IntegrationStrategies.m`
  (IntegrationStrategies.m:32-39) and declare no `::usage` of their own; they
  provide the definitions for the 13 symbols exported by
  `DiffExp/IntegrationStrategies.m`.  They contribute no additional exports.
- A pre-`Private` sweep of all files found NO exported symbols lacking
  `::usage` (everything before `Begin["`Private`"]` is usage text, comments,
  or `Get`/`BeginPackage` plumbing).
- `DiffExp/State.m` exports **79** symbols (the plan's "~90" estimate
  overcounts), plus it uses three `System`` symbols as config keys WITHOUT
  declaring usage (State.m:9-11): `System`AccuracyGoal`, `System`Variables`,
  `System`WorkingPrecision`.  These and the STRING config keys are covered in
  the "Non-symbol API surface" section at the end — they are part of the
  compatibility contract even though they are not DiffExp-owned exports.
- `DiffExp/RegularizedIntegration.m` exports exactly **15** symbols
  (RegularizedIntegration.m:45-76), as the plan states.

CROSS-REFERENCE SCOPES (word-boundary grep, context-qualified references
included):

- "FT layer" = `FeynmanTrick.m` (root loader), `FeynmanTrick/*.m`,
  `Scripts/*` .  FT references DiffExp ONLY via fully qualified context paths
  (its `BeginPackage["FeynmanTrick`"]` imports no DiffExp contexts;
  DiffExpIntegration.m:339-341 states this explicitly), so FT hits are exact.
  Hits on FT-internal homonyms (FT's own `"Verbosity"`/`"ExpansionOrder"`/
  `"DivisionOrder"`/... STRING option keys, FT-local `eps`/`x` variables) were
  inspected line-by-line and excluded; only genuine DiffExp-symbol references
  are listed.
- "Tests/Examples" = `Tests/*.m` + `Reference/Examples/*.m`.
  `Reference/DiffExp_original.m` is EXCLUDED (it is the upstream monolith
  source — every symbol "matches" there as a definition, not a usage).
  Counts are given as yes(n)/no; n = matching lines in scope.

DISPOSITION CLASSES (per RewritePlan §6 M0 item (2), module map of §3.2):

- **kept-in-API (Name)** — DiffExp2 provides a public equivalent in API.m /
  Config.m; old name or the stated new name.  Final signatures live in
  `Docs/specs/API.md` / `Docs/specs/Config.md`.
- **absorbed → Module.m** — the functionality is covered inside the named
  DiffExp2 module; no public symbol unless the row says a named replacement
  is required (FT shim).
- **dropped (reason)** — functionality intentionally not carried.

HARD RULE (RewritePlan §6 M5, R4): **M5 cutover may not begin while any
FT-referenced symbol is unclassified.**  All 28 FT-referenced symbols are
classified below and consolidated in the shim list at the end.  Exact call
sites are given for every FT-referenced symbol.

DiffExp2 module key: Tol = Tolerances.m, Cfg = Config.m, ES = EpsSeries.m,
SS = SectorSeries.m, Ind = Indicial.m, Sol = Solve.m, Tr = Transport.m,
Int = Integrate.m, API = API.m.

---

## 1. DiffExp/DiffExp.m (umbrella package) — 3 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| CurrentConfiguration | DiffExp/DiffExp.m:79 | no | yes(2) | kept-in-API (CurrentConfiguration; Cfg-backed) |
| LoadConfiguration | DiffExp/DiffExp.m:80 | **YES** — FeynmanTrick/DiffExpIntegration.m:367, 460 | yes(40) | kept-in-API (LoadConfiguration; resets to validated defaults then applies — same reset semantics, FT relies on it for the lower/upper re-transport at DiffExpIntegration.m:460) |
| UpdateConfiguration | DiffExp/DiffExp.m:81 | **YES** — FeynmanTrick/DiffExpIntegration.m:411, 466 | yes(7) | kept-in-API (UpdateConfiguration; FT uses it to re-add DeltaPrescriptions after the LoadConfiguration reset — that re-add dance must keep working or be obsoleted by a documented "prescriptions survive reload" rule in Cfg) |

## 2. DiffExp/AnalyticContinuation.m — 2 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| PrepareAnalyticContinuation | DiffExp/AnalyticContinuation.m:9 | no | yes(2) | absorbed → Tr (per-chart prescription derivation: all vanishing factors, multiplicity parity, leading-coeff sign, consistency check → LocalSolution "Prescriptions" record per §3.1; review legacy-finding 3) |
| Project\[Theta]s | DiffExp/AnalyticContinuation.m:10 | no | no | absorbed → Tr (theta-deduplication internalized in the crossing operator; the new crossing data is structural, not expression-level theta scrubbing) |

## 3. DiffExp/Frobenius.m — 2 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| Frobenius1 | DiffExp/Frobenius.m:9 | no | no | dropped (single-solution scalar Frobenius helper of the strategy stack; dies with it per §2 I2) |
| FrobeniusSolutions | DiffExp/Frobenius.m:10 | no | yes(1) | absorbed → Sol (homogeneous local solutions are the core output of the ONE symbolic-eps solver) |

## 4. DiffExp/Integration.m — 4 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| DiffExpIntegrate | DiffExp/Integration.m:9 | no | yes(14) | absorbed → Int (series antiderivative primitive becomes the exact per-sector case table {b=0/b≠0}×{a+n+1<0/=0/>0}×{p=0/p>0}) |
| DiffExpIntegrate1 | DiffExp/Integration.m:10 | no | no | dropped (inner worker of DiffExpIntegrate) |
| UpdateIntReps | DiffExp/Integration.m:11 | no | yes(5) | dropped (precomputed ∫x^n Logx^k replacement-rule cache keyed by IMaxLogOrder; log powers are exact `p` tags now, closed forms in Int replace the rule industry) |
| IntReps | DiffExp/Integration.m:12 | no | yes(3) | dropped (the rule table itself; same reason) |

## 5. DiffExp/IntegrationStrategies.m — 13 exports
(definitions live in DiffExp/IntegrationStrategies/{Helpers,Default,VOP,Recurrence,ResonantRecurrence,Dispatch}.m; whole stack deleted per §2 I2)

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| SolveSimple | DiffExp/IntegrationStrategies.m:18 | no | no | dropped (no-homogeneous-part special case; trivial branch of the one solver) |
| SolveDefault | DiffExp/IntegrationStrategies.m:19 | no | yes(4) | dropped (per-eps-order Frobenius/Wronskian method — the D2 disease carrier; capability covered by Sol) |
| SolveVOP | DiffExp/IntegrationStrategies.m:20 | no | no | dropped (variation-of-parameters needs Wronskian inversion of collapsed towers; obsolete under symbolic-eps Frobenius) |
| SolveVOPAlt | DiffExp/IntegrationStrategies.m:21 | no | no | dropped (same) |
| SolveRationalRecurrence | DiffExp/IntegrationStrategies.m:22 | no | no | absorbed → Sol (the denominator-cleared polynomial recursion is ported as Sol's performance backbone, §3.2 Solve.m; the strategy wrapper dies) |
| RationalRecurrenceApplicableQ | DiffExp/IntegrationStrategies.m:23 | no | yes(4) | dropped (dispatch predicate; no dispatch in a one-solver design) |
| SolveSingularRecurrence | DiffExp/IntegrationStrategies.m:24 | no | no | absorbed → Sol (simple-pole non-resonant recursion = Sol's base case at symbolic eps) |
| SingularRecurrenceApplicableQ | DiffExp/IntegrationStrategies.m:25 | no | yes(4) | dropped (dispatch predicate) |
| DispatchStrategy | DiffExp/IntegrationStrategies.m:26 | no | yes(1) | absorbed → Sol (the single solver entry point Transport calls; the ctx/bVec/cache interface dies — note: DEBUG_DUMP_DISPATCH_DIR dumps (Dispatch.m:10) are the M3 replay-parity oracle, generated in M0 before this code freezes) |
| SolveGeneralSingularRecurrence | DiffExp/IntegrationStrategies.m:27 | no | no | absorbed → Sol (resonant / non-diagonalizable residue handling → explicit true-resonance log-chains + joint pseudo-resonant solve per §2 I2) |
| GeneralSingularRecurrenceApplicableQ | DiffExp/IntegrationStrategies.m:28 | no | no | dropped (dispatch predicate; Ind decides the structure exactly) |
| SolveFuchsianizedSingularRecurrence | DiffExp/IntegrationStrategies.m:29 | no | no | absorbed → Ind (the rank-reduction core (FuchsianizeLocal) is ported as Indicial preprocessing; solving then proceeds in Sol — banana L1 nilpotent double pole is the in-scope day-one case) |
| FuchsianizedSingularRecurrenceApplicableQ | DiffExp/IntegrationStrategies.m:30 | no | no | dropped (predicate; Ind reads the pole order from the exact matrix) |

## 6. DiffExp/LineSegmentation.m — 7 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| RelateLines | DiffExp/LineSegmentation.m:9 | no | yes(1) | absorbed → Tr (chart-map relation between adjacent segments; chart maps are recorded on LocalSolution per §3.1) |
| RelateLinesPoint | DiffExp/LineSegmentation.m:10 | no | no | absorbed → Tr (point variant of the same) |
| FindMatrixSingularities | DiffExp/LineSegmentation.m:11 | no | no | absorbed → Tr (exact complex roots determine true radii; suppressed Re/Re±Im projections remain regular predivision waypoints for symmetric matching — ledger L26) |
| GetLargestTerm | DiffExp/LineSegmentation.m:12 | no | no | absorbed → Tr (coupling-depth-discounted magnitude probe feeding the error estimate; reads order ExpansionOrder−ISafetyExpansionSubtract−(MaxCouplingOrder−1) at LineSegmentation.m:109-113 — the TWindow lesson) |
| GetMatricesPrecisionDistance | DiffExp/LineSegmentation.m:13 | no | no | dropped (Dynamic segmentation strategy only; v1 is Predivision-only per §3.2 Config.m) |
| CheckBoundaryConditionsAndReparametrize | DiffExp/LineSegmentation.m:14 | no | no | absorbed → API (TransportTo/PrepareBoundaryConditions input validation and bcs-line reparametrization) |
| GetMatchingPoint | DiffExp/LineSegmentation.m:15 | no | no | absorbed → Tr (match-point geometry; radius/DivisionOrder of BOTH adjacent charts — GetCPL/GetCPR lesson, R3) |

## 7. DiffExp/LocalSeries.m — 15 exports
(per §2 I2: "the finite-width machinery is subsumed; FuchsianizeLocal is ported")

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| LocalZValuation | DiffExp/LocalSeries.m:12 | no | yes(2) | absorbed → Ind (z-adic valuation support for the ported rank reduction) |
| LocalLeadingCoefficient | DiffExp/LocalSeries.m:13 | no | no | absorbed → Ind (leading coefficient at the chart center) |
| LocalConnectionMatrix | DiffExp/LocalSeries.m:14 | no | yes(1) | absorbed → Ind (gauge-transform bookkeeping of the reduction) |
| FuchsianizeLocal | DiffExp/LocalSeries.m:15 | no | yes(3) | absorbed → Ind (EXPLICIT port per §3.2 Indicial.m: Moser/shearing rank reduction with loud non-termination error, R11) |
| TrimFuchsianLattice | DiffExp/LocalSeries.m:16 | no | no | absorbed → Ind (lattice cleanup of the ported reduction; port-or-waive decided in Docs/specs/Indicial.md) |
| LaurentCoefficientsRational | DiffExp/LocalSeries.m:17 | no | yes(2) | absorbed → ES (rational expression → Laurent coefficient arrays without asymptotic objects; shared primitive) |
| ClearZDenominators | DiffExp/LocalSeries.m:18 | no | no | absorbed → Sol (denominator clearing before the recursion — the fast-path backbone, §3.2 Solve.m) |
| PrepareFiniteWidthData | DiffExp/LocalSeries.m:19 | no | yes(1) | absorbed → Sol (recursion coefficient dictionaries) |
| RecursiveFiniteWidthSolve | DiffExp/LocalSeries.m:20 | no | yes(13) | absorbed → Sol (one-sector recursion → symbolic-eps sector-family recursion) |
| SolveLocalFuchsianSeries | DiffExp/LocalSeries.m:21 | no | yes(1) | absorbed → Sol (the end-to-end local solve = Ind preprocessing + Sol recursion pipeline) |
| FiniteWidthCoefficient | DiffExp/LocalSeries.m:22 | no | yes(15) | absorbed → SS (c[k,n,comp] accessor on Sector "Coeffs", §3.1) |
| FiniteWidthEvaluate | DiffExp/LocalSeries.m:23 | no | no | absorbed → SS (truncated local-solution evaluation) |
| RationalMatrixZEpsLaurentAssoc | DiffExp/LocalSeries.m:24 | no | yes(1) | absorbed → ES (matrix → coefficient-array expansion in (z, eps)) |
| ApplyGaugeToSolution | DiffExp/LocalSeries.m:25 | no | yes(1) | absorbed → Ind (undo the rank-reduction gauge on solutions) |
| ValidateFiniteWidthSolution | DiffExp/LocalSeries.m:26 | no | yes(3) | absorbed → Sol (ODE-residual check becomes the always-on per-chart invariant of §3.1) |

## 8. DiffExp/MatrixLoading.m — 8 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| LoadMatrices | DiffExp/MatrixLoading.m:9 | no (FT loads via MatrixDirectory config → LoadConfiguration; never calls LoadMatrices directly) | yes(1) | kept-in-API (LoadSystem: full-format `d<var>_full.m` primary, closed form, legacy slice dirs for parity transport only; accepts both `eps` and `\[Epsilon]` as MatrixLoading.m:118 does) |
| PrepareMatrices | DiffExp/MatrixLoading.m:10 | no | yes(1) | absorbed → Tr (line-restricted system preparation) |
| PrepareMatricesFrom1 | DiffExp/MatrixLoading.m:11 | no | no | absorbed → Tr (reuse-previous-segment variant; superseded by exact-matrix + per-chart expansion) |
| PrepareMatricesFrom | DiffExp/MatrixLoading.m:12 | no | no | absorbed → Tr (same) |
| PrepareMatricesFactored | DiffExp/MatrixLoading.m:13 | no | no | absorbed → Tr (the Factored cache layer dies; exact eps-rational matrix carried, charts expand on demand with segment block caching — ledger) |
| PrepareMatricesExpanded | DiffExp/MatrixLoading.m:14 | no | no | absorbed → Tr (same for the Expanded layer) |
| ClearMatrices | DiffExp/MatrixLoading.m:15 | no | no | dropped (invalidation hook for cache layers that no longer exist; Cfg reload semantics cover state reset) |
| InitializeIntegrationSequence | DiffExp/MatrixLoading.m:16 | no | no | absorbed → Sol (block dependency ordering; the chain lengths feed couplingDepth → TWindow per §3.1) |

## 9. DiffExp/Mobius.m — 9 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| GetMobius | DiffExp/Mobius.m:9 | no | yes(1) | absorbed → Tr (Mobius chart maps; "ChartMap" affine with Mobius optional per §3.1; Int REJECTS Mobius charts loudly v1 per §3.2) |
| GetLineRescaled | DiffExp/Mobius.m:10 | no | no | absorbed → Tr (chart rescaling incl. RadiusOfConvergence coefficient-magnitude control and the SetPrecision[...,2·WP] input raise at Mobius.m:48,66 — ledger) |
| GetMobiusCPL | DiffExp/Mobius.m:11 | no | no | absorbed → Tr (chart sizing, Mobius variant) |
| GetMobiusCPR | DiffExp/Mobius.m:12 | no | no | absorbed → Tr |
| GetCPLRep | DiffExp/Mobius.m:13 | no | no | absorbed → Tr (center-point replacement plumbing) |
| GetCPL | DiffExp/Mobius.m:14 | no | no | absorbed → Tr (next-center solved so the eval point sits at radius/k of BOTH adjacent charts (Mobius.m:98-142) — the R3 ill-conditioning defense, ledger) |
| GetCPR | DiffExp/Mobius.m:15 | no | no | absorbed → Tr (same, rightward) |
| FindNextCenterPointL | DiffExp/Mobius.m:16 | no | no | absorbed → Tr (next-center placement) |
| FindNextCenterPointR | DiffExp/Mobius.m:17 | no | no | absorbed → Tr |

## 10. DiffExp/Pade.m — 4 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| GetPade | DiffExp/Pade.m:9 | no | yes(1) | absorbed → SS (Pade evaluation ported as the evaluation accelerator per §3.2 SectorSeries.m DECISION; the silent fallback at Pade.m:44-46 becomes LOUD) |
| SEval1 | DiffExp/Pade.m:10 | no | no | absorbed → SS (evaluate with optional Pade) |
| SEval2 | DiffExp/Pade.m:11 | no | no | absorbed → SS (evaluate at a point with analytic continuation → SS "evaluate (with branch rule)") |
| SEval | DiffExp/Pade.m:12 | no | no | absorbed → SS (same family; $MinPrecision-floored evaluation blocks at Pade.m:34,70,80 — ledger) |

## 11. DiffExp/RegularizedIntegration.m — 15 exports
(3015 lines; this table is the per-export delete-vs-reimplement record the execution review demanded)

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| ApplyRegularizationStep | DiffExp/RegularizedIntegration.m:45 | no | yes(3) | dropped (the power-raising integration-by-parts ladder exists only because exponents were collapsed; Int's closed case table with (a+n+1+b·eps)^(j+1) denominators integrates sectors exactly) |
| RegularizeIntegrand | DiffExp/RegularizedIntegration.m:47 | no | yes(5) | dropped (driver of the same ladder) |
| IntegrateSingularTerm | DiffExp/RegularizedIntegration.m:49 | no | no | absorbed → Int (per-sector endpoint integration) |
| IntegrateSingularTermLaurent | DiffExp/RegularizedIntegration.m:51 | no | yes(15) | absorbed → Int (Laurent-aware variant; its MinPower/Coefficients honest-window output becomes EpsWindow arithmetic in ES) |
| IntegrateDecomposition | DiffExp/RegularizedIntegration.m:53 | no | yes(1) | absorbed → Int (integrate a sector decomposition over an interval = Int acting on LocalSolution sectors) |
| IntegrateDecompositionLaurent | DiffExp/RegularizedIntegration.m:55 | no | yes(9) | absorbed → Int (same, Laurent windows) |
| EvaluateLimitAtSingularity | DiffExp/RegularizedIntegration.m:57 | no | yes(2) | absorbed → API (EndpointLimit: constant of the (0,0,0)-sector, b≠0 dropped exactly, divergence loud — §3.3) |
| EvaluateEndpointLimitSectors | DiffExp/RegularizedIntegration.m:59 | **YES** — FeynmanTrick/DiffExpIntegration.m:1015 | yes(4) | kept-in-API (**EndpointLimit**; the per-sector resolution this function fakes by fitting is native — FT SHIM required) |
| FitResidualEndpointSectors | DiffExp/RegularizedIntegration.m:61 | no | yes(4) | dropped (**dies with the Prony/N-root fitter**: sector structure is exact indicial data in DiffExp2, nothing to fit; pentagon proved the drops precede fitting anyway — §2 D2, §10) |
| IntegrateSegmentData | DiffExp/RegularizedIntegration.m:63 | no | no | absorbed → Int (per-segment integration of saved expansions) |
| IntegratePiecewiseSaved | DiffExp/RegularizedIntegration.m:65 | no | no | absorbed → API (IntegrateOverLine assembly across segments) |
| DefiniteIntegral | DiffExp/RegularizedIntegration.m:67 | no | yes(12) | kept-in-API (**IntegrateOverLine**, plain form; battery tests pin its values) |
| IndefiniteIntegral | DiffExp/RegularizedIntegration.m:69 | no | yes(6) | absorbed → Int (per-sector exact antiderivatives exist as objects; whether API exposes a piecewise indefinite form is an API.md decision — open question 2) |
| DefiniteIntegralWithPrefactor | DiffExp/RegularizedIntegration.m:71 | no (DiffExpIntegration.m:846 is a comment only) | no | absorbed → Int (prefactor x^α(c−x)^β·r(x) = rational multiply at object level + closed forms, §3.3 "integrate") |
| DefiniteIntegralWithPrefactorLaurent | DiffExp/RegularizedIntegration.m:76 | **YES** — FeynmanTrick/DiffExpIntegration.m:875; Scripts/eval_dump_generic.m:21 | yes(3) | kept-in-API (**IntegrateOverLine**, Laurent-aware prefactor form — THE main FT integration entry point; FT SHIM required; honest EpsWindow replaces the epsMinPower offset convention) |

## 12. DiffExp/SeriesOps.m — 23 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| SApply | DiffExp/SeriesOps.m:9 | no | no | absorbed → SS (coefficient-array map) |
| SExpand | DiffExp/SeriesOps.m:10 | no | yes(11) | absorbed → SS (normalize+chop; chop threshold from Tol chopFloor) |
| SN | DiffExp/SeriesOps.m:11 | no | no | absorbed → SS (numericize coefficients) |
| SSN | DiffExp/SeriesOps.m:12 | no | no | absorbed → SS (same) |
| SMultiply | DiffExp/SeriesOps.m:13 | no | no | absorbed → SS (series product with honest truncation windows) |
| SeriesCoefficientMinus | DiffExp/SeriesOps.m:14 | no | no | absorbed → SS (end-indexed coefficient access; trivial on arrays) |
| ApplyAnalyticContinuation | DiffExp/SeriesOps.m:15 | no | no | absorbed → Tr (application of crossing data to series → the crossing operator) |
| SafeReplaceSeries11 | DiffExp/SeriesOps.m:16 | no | no | dropped (SeriesData-replacement fragility workaround; array representation removes the failure mode) |
| NormalizeLogPower | DiffExp/SeriesOps.m:19 | no | yes(6) | dropped (numeric rationalization of Logx powers; `p` is exact — D2 cure.  Also kills the unexported-symbol cross-package no-op bug class this symbol caused, see memory/wolfram-package-context-traps.md) |
| MaxLogxPower | DiffExp/SeriesOps.m:20 | no | yes(1) | dropped (log-power scanning; exact `p` tags) |
| LogxPowerRange | DiffExp/SeriesOps.m:21 | no | no | dropped (same) |
| LogxCoeff | DiffExp/SeriesOps.m:22 | no | no | absorbed → SS (coefficient of Logx^p = the p-graded sector slice, a direct structural read) |
| LogxCoeffNS | DiffExp/SeriesOps.m:23 | no | yes(1) | absorbed → SS (non-series variant of the same) |
| LogxCoeffList | DiffExp/SeriesOps.m:24 | no | no | absorbed → SS (all p-slices) |
| MatrixMultiplySExpand | DiffExp/SeriesOps.m:27 | no | no | absorbed → SS (matrix×series with truncation) |
| MatrixPowerSExpand | DiffExp/SeriesOps.m:28 | no | no | absorbed → SS |
| DiffExpSeries | DiffExp/SeriesOps.m:31 | no | no | absorbed → SS (series constructor) |
| SeriesAlways | DiffExp/SeriesOps.m:32 | no | no | absorbed → SS (constant-term-safe constructor; trivial on arrays) |
| LeadingCoefficientSeries | DiffExp/SeriesOps.m:33 | no | no | absorbed → SS (leading coefficient; numerical-zero-skipping lesson generalizes to matching solves — ledger) |
| SeriesMinPower | DiffExp/SeriesOps.m:34 | no | no | absorbed → SS (array bounds / TWindow read) |
| SeriesMaxPower | DiffExp/SeriesOps.m:35 | no | no | absorbed → SS |
| DecreaseSeriesOrderBy | DiffExp/SeriesOps.m:36 | no | no | absorbed → SS (truncation-order management; used by the two-point error probe) |
| SD | DiffExp/SeriesOps.m:39 | no | no | absorbed → SS (differentiate; log-tag-aware derivative is structural) |

## 13. DiffExp/SingularityDecomposition.m — 3 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| DecomposeSingularity | DiffExp/SingularityDecomposition.m:25 | **YES** — FeynmanTrick/DiffExpIntegration.m:822 | yes(9) | absorbed → SS (sector-native output makes the decomposition a direct read of LocalSolution "Sectors"; §3.2 SectorSeries.m: "absorbs old SingularityDecomposition.m's role; the FT call site gets a named replacement API" — FT SHIM required) |
| DecomposeSingularityAll | DiffExp/SingularityDecomposition.m:27 | no | yes(1) | absorbed → SS (map over integrals; same replacement API) |
| PrintDecomposition | DiffExp/SingularityDecomposition.m:29 | no | yes(1) | dropped (pretty-printer for a defunct intermediate format) |

## 14. DiffExp/State.m — 79 exports

### 14a. Configuration option keys (State.m:12-27) — 16

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| ChopPrecision | DiffExp/State.m:12 | **YES** — DiffExpIntegration.m:349 | yes(23) | kept-in-API (Cfg key; ChopPrecision<WorkingPrecision invariant + LinearSolveChopPrecision auto-sync, DiffExp.m:121-130; thresholds derived in Tol) |
| DeltaPrescriptions | DiffExp/State.m:13 | **YES** — DiffExpIntegration.m:412, 467 | yes(17) | kept-in-API (Cfg key; ±Iδ parsing + irreducibility validation, DiffExp.m:140-168; feeds LocalSolution "Prescriptions") |
| DivisionOrder | DiffExp/State.m:14 | **YES** — DiffExpIntegration.m:354 | yes(19) | kept-in-API (Cfg key; FT pins k=4 at DiffExpIntegration.m:272, classic default 3) |
| EpsilonOrder | DiffExp/State.m:15 | **YES** — DiffExpIntegration.m:351 | yes(14) | kept-in-API (Cfg key; in DiffExp2 the requested order feeds the static budget of §3.4 and EpsWindow honesty) |
| ExpansionOrder | DiffExp/State.m:16 | **YES** — DiffExpIntegration.m:350 | yes(32) | kept-in-API (Cfg key; t-truncation order, with TWindow honesty) |
| IntegrationStrategy | DiffExp/State.m:17 | **YES** — DiffExpIntegration.m:358 | yes(10) | dropped (ONE solver, no strategy selection; **FT call site must be deleted at M5** — see shim list) |
| LineParameter | DiffExp/State.m:18 | no | no | kept-in-API (Cfg key; binds the working variable, replaces the `DiffExp`Symbols`x` mutation at DiffExp.m:136-139) |
| LogFile | DiffExp/State.m:19 | no | no | dropped (session logging via $Output/$Messages append, DiffExp.m:107-114; NOT in the §3.2 kept-config surface — waiver to be recorded in Docs/specs/Config.md; see open question 4) |
| MatrixDirectory | DiffExp/State.m:20 | **YES** — DiffExpIntegration.m:347 | yes(30) | kept-in-API (Cfg key; LoadSystem source; M5 switches FT export to full format) |
| RadiusOfConvergence | DiffExp/State.m:21 | no | yes(5) | kept-in-API (Cfg key; chart rescaling keeps high-order coefficients O(1) — banana REQUIRES RoC=10, ledger) |
| RationalizationTolerance | DiffExp/State.m:22 | **YES** — as FEC key at DiffExpIntegration.m:98, 109, 127, 154, 182, 924 | yes(4) | kept-in-API (Cfg key with semantics moved to Tol: named relative tolerance replacing FT's 10^-40 absolute fallbacks at the cited sites) |
| SegmentationStrategy | DiffExp/State.m:23 | **YES** — DiffExpIntegration.m:356 (passes "Predivision") | yes(1) | kept-in-API (Cfg key; Predivision only v1, other values = loud error) |
| UseMobius | DiffExp/State.m:24 | **YES** — DiffExpIntegration.m:352 (False, "Required for integration!") | yes(30) | kept-in-API (Cfg key; Int rejects Mobius charts loudly v1 — same contract, now enforced instead of commented) |
| UsePade | DiffExp/State.m:25 | **YES** — DiffExpIntegration.m:353 | yes(30) | kept-in-API (Cfg key; SS Pade accelerator toggle) |
| UseRationalRecurrence | DiffExp/State.m:26 | **YES** — DiffExpIntegration.m:357 | yes(33) | dropped (the recursive-strategy opt-in dies with the strategy stack; denominator-cleared recursion is Sol's unconditional fast path; **FT call site deleted at M5**, also Scripts/run_ft_stepwise.m:191 passes the FT-level option string) |
| Verbosity | DiffExp/State.m:27 | **YES** — DiffExpIntegration.m:355 (the other 64 FT grep hits are FT's own "Verbosity" string option) | yes(42) | kept-in-API (Cfg key) |

### 14b. State accessors (State.m:30-32) — 3

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| FEC | DiffExp/State.m:30 | **YES** — DiffExpIntegration.m:98, 109, 127, 154, 182, 328 (ValueQ load-guard), 372 (FEC[System`Variables]), 924 | yes(2) | absorbed → Cfg (the validated accessor of §3.2 — silent-miss-by-construction-impossible; **FT SHIM required**: all 8 sites move to the new accessor at M5) |
| DiffExpConfiguration | DiffExp/State.m:31 | no | yes(21) | absorbed → Cfg (internal storage; tests that poke it retarget at M6 battery disposition) |
| DefaultConfiguration | DiffExp/State.m:32 | no | no | absorbed → Cfg (option schema defaults) |

### 14c. Value accessors (State.m:35-49) — 15
(cached `...Val` mirrors of config; the mirror pattern itself is a D1 source and dies — the validated accessor is the only read path)

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| ChopPrecisionVal | DiffExp/State.m:35 | no | yes(1) | absorbed → Tol (chopFloor derivation) |
| LinearSolveChopPrecisionVal | DiffExp/State.m:36 | no | yes(1) | absorbed → Tol (ChopPrecision/LinearSolveChopPrecision sync — ledger) |
| CrosscheckChopPrecision | DiffExp/State.m:37 | no | no | absorbed → Tol (crosscheck threshold, named) |
| ExternalScalesVal | DiffExp/State.m:38 | **YES** — DiffExpIntegration.m:418 (verbosity≥1 debug Print only) | yes(3) | absorbed → Cfg (kinematic-scale list accessor; FT site is a debug print — shim trivially or delete the print at M5) |
| LineParameterVal | DiffExp/State.m:39 | no | no | absorbed → Cfg |
| MatrixDirectoryVal | DiffExp/State.m:40 | no | no | absorbed → Cfg |
| EpsilonOrderVal | DiffExp/State.m:41 | no | no | absorbed → Cfg |
| FEAccuracyGoal | DiffExp/State.m:42 | no | no | absorbed → Cfg |
| FEWorkingPrecision | DiffExp/State.m:43 | **YES** — DiffExpIntegration.m:78 (activeNumericPrecision fallback chain) | yes(1) | absorbed → Cfg (**FT SHIM required**: validated WorkingPrecision accessor; FT's Quiet[Check[...]] fallback to its own config dies with A1) |
| DeltaPrescriptionsVal | DiffExp/State.m:44 | no | no | absorbed → Cfg |
| UseMobiusVal | DiffExp/State.m:45 | no | no | absorbed → Cfg |
| RadiusOfConvergenceVal | DiffExp/State.m:46 | no | no | absorbed → Cfg |
| DivisionOrderVal | DiffExp/State.m:47 | no | no | absorbed → Cfg |
| ExpansionOrderVal | DiffExp/State.m:48 | no | yes(1) | absorbed → Cfg |
| MaxCouplingOrder | DiffExp/State.m:49 | no | no | absorbed → Sol (longest coupled chain, MatrixLoading.m:383; becomes couplingDepth feeding TWindow per §3.1 — legacy finding 6) |

### 14d. State variables (State.m:52-79) — 28

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| AnalyticContinuationFailed | DiffExp/State.m:52 | no | no | absorbed → Tr (prescription conflict/indeterminacy → consistency-checked Prescriptions record with loud error; the soft one-sided-validity mode survives only behind the kept "AbortOnAnalyticContinuationFail" key) |
| AnalyticContinuationReplacements | DiffExp/State.m:53 | no | no | absorbed → Tr (crossing-operator data, now structural: phase × unipotent log-chain mixing, §3.2 Transport.m) |
| AnalyticContinuationReplacementsAssociation | DiffExp/State.m:54 | no | no | absorbed → Tr (per-line cache of the same) |
| BenchmarkData | DiffExp/State.m:55 | no | no | dropped (ad-hoc timing store, Transport.m:157-237; M4 benchmark gate uses external harness) |
| CurrentSingularityWasAddedFromSquareRoot | DiffExp/State.m:56 | no | no | dropped (sqrt-derived-prescription flag; sqrt x-dependence OUT OF SCOPE v1 per §1 NON-GOALS — ledger v1.1; see open question 5) |
| CurrentSingularityHasIDeltaPrescription | DiffExp/State.m:57 | no | no | absorbed → Tr (per-chart prescription presence is a field of the LocalSolution Prescriptions record, not a mutable global) |
| DEqnSquareRoots | DiffExp/State.m:58 | no | no | dropped (square roots in DE matrices; OUT OF SCOPE v1, loud error at LoadSystem; machinery ledger-documented for v1.1.  NOTE §3.1 keeps the auto-prescription RULE ("sqrt factors auto-prescribed as in old State.m DEqnSquareRoots") — the rule moves to Tr prescription derivation, the matrix-sqrt support does not — open question 5) |
| MultivaluedFail | DiffExp/State.m:59 | no | no | dropped (silent failure flag — D1; loud error replaces) |
| UserDeltaPrescriptions | DiffExp/State.m:60 | no | no | absorbed → Cfg (user-vs-auto-derived prescription separation, kept for reload semantics) |
| UsingClosedFormMatrix | DiffExp/State.m:61 | no | no | absorbed → API (LoadSystem input-format handling) |
| DEqnMatricesFactored | DiffExp/State.m:62 | no | yes(2) | absorbed → Tr (per-line factored matrix cache → exact matrix + per-chart expansion + segment block caching) |
| DEqnMatricesFactoredClosedForm | DiffExp/State.m:63 | no | no | absorbed → Tr (closed-form variant of the same cache) |
| DEqnMatricesExpanded | DiffExp/State.m:64 | no | yes(1) | absorbed → Tr (expanded variant) |
| NumIntegrals | DiffExp/State.m:65 | no (FT's `transportResult["NumIntegrals"]` reads at DiffExpIntegration.m:917 / check_transport_ode_residual.m:28 are the STRING result key, not this symbol) | yes(9) | absorbed → API (system dimension on the loaded-system object; the TransportTo RESULT key "NumIntegrals" (old Transport.m:506,1243) is preserved — see Non-symbol API surface) |
| IntegrationSequence | DiffExp/State.m:66 | no | no | absorbed → Sol (block solve ordering, MatrixLoading.m:358-375) |
| ExpansionMatrices | DiffExp/State.m:67 | no | no | absorbed → API (loaded partial-derivative matrices = LoadSystem state) |
| ExpansionMatricesCanonical1 | DiffExp/State.m:68 | no | no | absorbed → API (canonical dlog `d_1.m` input path, MatrixLoading.m:62-82; LIVE via Reference/Examples/FivePointNonPlanar_example.m:117 and test_five_point_nonplanar.m — open question 1: §3.2 LoadSystem does not name this format) |
| ExpansionMatricesClosedForm | DiffExp/State.m:69 | no | no | absorbed → API (closed-form input matrices) |
| AlphabetLogs | DiffExp/State.m:70 | no | no | absorbed → API (canonical-form letter extraction, MatrixLoading.m:71-76; tied to open question 1) |
| AlphabetLogRules | DiffExp/State.m:71 | no | no | absorbed → API (Log[letter] placeholder rules for the canonical path) |
| AlphabetLogRulesFactored | DiffExp/State.m:72 | no | no | dropped (dead state: written empty at DiffExp.m:185, never read anywhere in DiffExp/) |
| AlphabetLogRulesExpanded | DiffExp/State.m:73 | no | no | dropped (dead state: DiffExp.m:186 only) |
| MatricesIrreducibleFactors | DiffExp/State.m:74 | **YES** — DiffExpIntegration.m:211, 213 (FT WRITES it: appendMatrixFactors extends the alphabet with IBP factors), 402 (reads it for prescription derivation) | no | absorbed → API (the singular-factor alphabet of the loaded system; **FT SHIM required and it is the hardest one**: FT mutates frozen-library state — DiffExp2 must expose a sanctioned read-AND-EXTEND API (e.g. LoadSystem option or AddSingularFactors) so segmentation and prescriptions see the IBP-induced factors) |
| CurrCrosscheckFlags | DiffExp/State.m:75 | no | no | absorbed → Cfg (active crosscheck set from CrosscheckLevel/CrosscheckFlags keys, DiffExp.m:115-120) |
| CrosscheckFlags | DiffExp/State.m:76 | no | yes(1) | absorbed → Cfg (flag registry; Crosscheck* keys kept per §3.2; SingularityCheck default-on — ledger) |
| DiffExpExtensions | DiffExp/State.m:77 | no | no | dropped (extension-hook registry, DiffExp.m:97,192; zero consumers in repo) |
| LogStream | DiffExp/State.m:78 | no | no | dropped (LogFile machinery; see LogFile row) |
| LastErrorContext | DiffExp/State.m:79 | no | no | dropped (post-hoc error-context global; replaced by loud errors that NAME (chart, sector, order) in the payload per §3.1) |

### 14e. Internal constants (State.m:82-97) — 16

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| ISeriesChangeCoefficient | DiffExp/State.m:82 | no | no | absorbed → Tr (series-inversion order multiplier in RelateLines, LineSegmentation.m:204; subsumed by SS's re-expand truncation contract (target order + (Δ/R)^N tail bound)) |
| IMaxLogOrder | DiffExp/State.m:83 | no | yes(3) | dropped (log-order cap; exact `p` tags make caps meaningless) |
| IMaxLogOrderDefault | DiffExp/State.m:84 | no | no | dropped (same) |
| ICheckMultivaluedChop | DiffExp/State.m:85 | no | no | absorbed → Tol (SingularityCheck chop at Transport.m:481 becomes a named WP-derived threshold) |
| ICrossCheckPrintResultOrder | DiffExp/State.m:86 | no | no | absorbed → Cfg (crosscheck reporting orders) |
| ICrossCheckVerifyResultOrder | DiffExp/State.m:87 | no | no | absorbed → Cfg |
| ISafetyDigits | DiffExp/State.m:88 | no | no | absorbed → Tr (predivision digit budget: DigitsNeeded = AccuracyGoal + Ceiling[Log10[#segments]] + safety, Transport.m:758-760 — ledger) |
| ISafetyExpansionSubtract | DiffExp/State.m:89 | no | no | absorbed → Tr (coupling-depth t-order discount — TWindow lesson, LineSegmentation.m:109-113) |
| IExpansionOrdersAveraging | DiffExp/State.m:90 | no | no | absorbed → Tr (adaptive expansion-order search, AccuracyGoalValidate "Before": averaging window 3, Transport.m:776-841) |
| IExpansionOrderIncrease | DiffExp/State.m:91 | no | no | absorbed → Tr (search step sizes) |
| IExpansionOrderDecrease | DiffExp/State.m:92 | no | no | absorbed → Tr |
| IExpansionOrderIncrease2 | DiffExp/State.m:93 | no | no | absorbed → Tr |
| IDigitsSurplusDecreaseExpansionOrder | DiffExp/State.m:94 | no | no | absorbed → Tr (3-digit surplus threshold of the search) |
| ICurrEvalErrorSeriesDecrease | DiffExp/State.m:95 | no | no | absorbed → Tr (error-probe order reduction Ceiling[0.7·MaxCouplingOrder]+2, State.m:222 / Transport.m:913 — ledger formula) |
| IDecreaseOrderByErrorPrecise | DiffExp/State.m:96 | no | no | absorbed → Tr (precise-mode variant) |
| IMinExpansionOrder | DiffExp/State.m:97 | no | no | absorbed → Tr (min order 10 floor of the adaptive search) |

(These become NAMED constants in the Tr/Tol specs; any numeric threshold
among them is derived in Tol per the one-tolerance-module rule.)

### 14f. Helper (State.m:100) — 1

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| SquareRootPrescriptionsAdded | DiffExp/State.m:100 | no | no | dropped (sqrt machinery v1.1 ledger item; the auto-prescription rule itself is ported into Tr prescription derivation per §3.1 — open question 5) |

## 15. DiffExp/Symbols.m — 6 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| \[Epsilon] | DiffExp/Symbols.m:10 | no (FT uses its own FTeps; bridge is via exported matrix files) | yes(61) | kept-in-API (the regulator symbol; both `eps` and `\[Epsilon]` accepted everywhere per §3.2 API.m) |
| eps | DiffExp/Symbols.m:11 | no direct code reference, but a HARD IMPLICIT dependency: FT's MatrixExport strips contexts so matrix files contain the bare name `eps` (MatrixExport.m:43-72, 92-94), which only works because `eps := \[Epsilon]` (Symbols.m:22) | yes(63) | kept-in-API (alias; LoadSystem and the full-format reader MUST resolve bare `eps` in matrix files — this is load-bearing for the entire FT pipeline and for M5's ExportGeneralMatrix cutover) |
| Logx | DiffExp/Symbols.m:12 | **YES** — DiffExpIntegration.m:146 (Logx → Log[x] at final evaluation); Scripts/check_transport_ode_residual.m:56, 58 | yes(43) | absorbed → SS (log tags are exact (eps·Logx)^p/p! structure; a symbolic log carrier appears only in rendered output (ToPiecewise-equivalent); **FT SHIM**: structured sector access replaces Logx scraping — FT's own comment at DiffExpIntegration.m:57-61 documents how fragile the current contract is) |
| \[Theta]p | DiffExp/Symbols.m:13 | **YES** — DiffExpIntegration.m:64, 69 ($thetaPlusRules/$thetaMinusRules); Scripts/check_transport_ode_residual.m:85, 86 | yes(1) | absorbed → SS evaluate branch rule + Tr crossing operator (side-of-cut becomes an explicit branch argument; **FT SHIM**: theta substitution rules replaced by passing the branch/direction) |
| \[Theta]m | DiffExp/Symbols.m:14 | **YES** — DiffExpIntegration.m:65, 70; Scripts/check_transport_ode_residual.m:85, 86 | yes(1) | absorbed → SS/Tr (same; note the principal-branch convention: sign −1 applies Logx → Logx − 2πi·θm as a SHIFT — ledger) |
| x | DiffExp/Symbols.m:17 | **YES** — DiffExpIntegration.m:143, 153, 700, 831, 928 (FT reads the bound line-parameter symbol); Scripts/check_transport_ode_residual.m:29 | yes(220 — almost all are test-LOCAL variables named x; genuine coupling is via the LineParameter binding at DiffExp.m:138) | absorbed → Cfg (variable-context pinning per §3.2 Config.m; the mutable global alias dies; **FT SHIM**: a Cfg accessor for the working variable replaces `DiffExp`Symbols`x` reads — this also fixes the FIRE variable-context pinning lesson) |

## 16. DiffExp/Transport.m — 4 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| PrepareBoundaryConditions | DiffExp/Transport.m:23 | no | yes(21) | kept-in-API (PrepareBoundaryConditions: closed-form eps expressions with auto-expansion, leading-asymptotics extraction incl. x^(a+b·eps) scaling, Log→Logx rewriting (Transport.m:31-120), per-integral "?" wildcards (Transport.m:63-65) — full contract in Docs/specs/API.md per legacy finding 9) |
| IntegrateSystem | DiffExp/Transport.m:24 | no | yes(4) | kept-in-API (general-series-solutions-along-a-line mode; covered by TransportTo's symbolic-indeterminate mode / SolveAtPoint — exact surface name fixed in Docs/specs/API.md, open question 3) |
| TransportTo | DiffExp/Transport.m:25 | **YES** — DiffExpIntegration.m:440, 474 (4-arg form: bcs-with-point, target assoc, endpoint 1, SaveExpansions True) | yes(68) | kept-in-API (TransportTo: assoc chaining, BIDIRECTIONAL, singular-endpoint mode returning the series (Transport.m:605-613, 1053-1096), SaveExpansions arg, result keys "ErrorEstimates"/"SegmentData"/"NumIntegrals"/"EpsilonOrder" (Transport.m:1238-1243) — FT consumes the result keys, see Non-symbol API surface) |
| ToPiecewise | DiffExp/Transport.m:26 | no | yes(13) | kept-in-API (ToPiecewise equivalent over saved segment data; note DiffExp.m:75's claim it lives in Pade` is stale — it is Transport.m:26) |

## 17. DiffExp/Utilities.m — 28 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| PrintDebug | DiffExp/Utilities.m:9 | no | no | absorbed → Cfg (Verbosity-gated reporting primitive shared by all modules) |
| PrintInfo | DiffExp/Utilities.m:10 | no | yes(1) | absorbed → Cfg |
| PrintWarning | DiffExp/Utilities.m:11 | no | yes(3) | absorbed → Cfg (warnings only where the result is still complete; A1 forbids warn-and-continue degradation) |
| ReportError | DiffExp/Utilities.m:12 | no | no | absorbed → Cfg (the loud-error primitive; payload must name (chart, sector, order) per §3.1) |
| PrintMobiusNormalized | DiffExp/Utilities.m:13 | no | no | dropped (cosmetic Mobius-term printing) |
| AllSameQ | DiffExp/Utilities.m:16 | no | no | dropped (generic one-liner) |
| CA | DiffExp/Utilities.m:17 | no | no | dropped (ConstantArray alias) |
| GetCases | DiffExp/Utilities.m:18 | no | no | dropped (generic Cases wrapper) |
| DependsQ | DiffExp/Utilities.m:19 | no | yes(8) | dropped (generic dependence test; reimplement locally if needed) |
| ZeroQ | DiffExp/Utilities.m:20 | no | no | dropped (exact-zero test inline) |
| R | DiffExp/Utilities.m:21 | no | no | dropped (ReplaceAll alias) |
| FirstOrNull | DiffExp/Utilities.m:22 | no | no | dropped (generic) |
| FindPivots | DiffExp/Utilities.m:23 | no | no | dropped (pivot scan for Wronskian pivot selection; dies with Wronskian.m — matching solves use rankTol from Tol) |
| SplitTimes | DiffExp/Utilities.m:24 | no | no | dropped (generic) |
| SplitSum | DiffExp/Utilities.m:25 | no | no | dropped (generic) |
| PChop | DiffExp/Utilities.m:28 | **YES** — DiffExpIntegration.m:177, 705, 733, 777, 783; Scripts/check_transport_ode_residual.m:35 | no | absorbed → Tol (chopFloor with named, WP-derived semantics; **FT SHIM required**: the 6 sites call the Tol chop) |
| LSPChop | DiffExp/Utilities.m:29 | no | no | absorbed → Tol (LinearSolveChopPrecision chop) |
| CPChop | DiffExp/Utilities.m:30 | no | no | absorbed → Tol (crosscheck chop) |
| IsPoint | DiffExp/Utilities.m:33 | no | no | absorbed → API (line/point input classification in LoadSystem/TransportTo parsing) |
| IsLine | DiffExp/Utilities.m:34 | no | no | absorbed → API |
| IntervalOverlapQ | DiffExp/Utilities.m:37 | no | no | absorbed → Tr (segment-interval bookkeeping in the marching loop — Transport.m:696, 1117; NOT dynamic-segmentation-only) |
| IntervalIntersec | DiffExp/Utilities.m:38 | no | no | absorbed → Tr (Transport.m:704-708, 1126-1132) |
| IntervalContainsQ | DiffExp/Utilities.m:39 | no | no | absorbed → Tr (Transport.m:682, 1049) |
| ExactLineQ | DiffExp/Utilities.m:42 | no | no | absorbed → Tr (exact-vs-inexact line normalization; Factor vs Together choice in matrix prep, MatrixLoading.m:276-315) |
| FactorOrTogether | DiffExp/Utilities.m:43 | no | no | absorbed → Tr (same) |
| NonFiniteExpressionQ | DiffExp/Utilities.m:46 | no | no | absorbed → Tr (matching-solve sanity guards, Transport.m:340-364; becomes part of the A2 structural asserts) |
| FiniteAbsMax | DiffExp/Utilities.m:47 | no | no | absorbed → Tr (error-probe magnitude extraction, Transport.m:580, 768-770, 961) |
| ReplaceSparseArrays | DiffExp/Utilities.m:48 | no | no | absorbed → Tr (normalization helper of the same guards) |

## 18. DiffExp/Wronskian.m — 4 exports

| Symbol | Defined at | FT layer? | Tests/Examples? | DISPOSITION |
|---|---|---|---|---|
| MatrixLogxInverse | DiffExp/Wronskian.m:9 | no | yes(1) | dropped (inversion of Logx-laden matrices; dies with the Wronskian/VOP industry per §2 I2 — sector recursion never inverts a fundamental matrix in this representation) |
| NullSpaceTryAgainOnFail | DiffExp/Wronskian.m:10 | no | no | dropped (precision-retry nullspace — D3-class ad-hoc tolerance; matching solves use Tol rankTol semantics) |
| CombineDifferentialEquationsHomogeneous | DiffExp/Wronskian.m:11 | no | no | dropped (reduction of first-order systems to higher-order scalar ODEs for Frobenius; Sol's matrix-level sector recursion replaces the scalar detour) |
| CombineDifferentialEquationsWithPivotSelection | DiffExp/Wronskian.m:12 | no | no | dropped (pivot-retry wrapper of the same; its consumer "HomogeneousSolve" config key dies too — see Non-symbol API surface) |

---

## Summary counts

| Disposition | Count |
|---|---|
| kept-in-API | **26** |
| absorbed | **151** |
| dropped | **52** |
| **Total** | **229** |

kept-in-API (26): CurrentConfiguration, LoadConfiguration, UpdateConfiguration,
LoadMatrices (→LoadSystem), EvaluateEndpointLimitSectors (→EndpointLimit),
DefiniteIntegral (→IntegrateOverLine), DefiniteIntegralWithPrefactorLaurent
(→IntegrateOverLine, Laurent prefactor form), ChopPrecision,
DeltaPrescriptions, DivisionOrder, EpsilonOrder, ExpansionOrder, LineParameter,
MatrixDirectory, RadiusOfConvergence, RationalizationTolerance,
SegmentationStrategy, UseMobius, UsePade, Verbosity, \[Epsilon], eps,
PrepareBoundaryConditions, IntegrateSystem, TransportTo, ToPiecewise.

Absorbed, by receiving module (sums to 151):

| Module | Count | Notable absorptions |
|---|---|---|
| Transport.m | 48 | prescriptions/AC machinery incl. ApplyAnalyticContinuation (7), segmentation+match geometry (5), matrix prep (5), Mobius/chart maps (9), matrix caches (3), expansion-order/error-probe/digit-budget constants (11), interval+line+finiteness helpers (8) |
| SectorSeries.m | 29 | SeriesOps algebra (18), Pade evaluation (4), DecomposeSingularity(+All) (2), FiniteWidth accessors (2), Logx/θp/θm (3) |
| Config.m | 24 | FEC + storage (3), value accessors (11), prescriptions/crosscheck state (5), print/error primitives (4), Symbols`x (1) |
| Solve.m | 13 | FrobeniusSolutions, 4 strategy capabilities, 5 finite-width solver pieces, InitializeIntegrationSequence, MaxCouplingOrder, IntegrationSequence |
| API.m | 13 | bcs validation, EndpointLimit + piecewise integration (2), LoadSystem state (8: incl. canonical-d_1 path), IsPoint/IsLine |
| Integrate.m | 8 | DiffExpIntegrate, 7 RegularizedIntegration integration entry points |
| Tolerances.m | 7 | PChop/LSPChop/CPChop, 3 chop-precision accessors, ICheckMultivaluedChop |
| Indicial.m | 7 | FuchsianizeLocal + 5 rank-reduction helpers, SolveFuchsianizedSingularRecurrence |
| EpsSeries.m | 2 | LaurentCoefficientsRational, RationalMatrixZEpsLaurentAssoc |

Dropped (52), by reason class:
- Strategy/Wronskian/Frobenius machinery obsolete under symbolic-eps solving
  (D2 cure): Frobenius1, SolveSimple, SolveDefault, SolveVOP, SolveVOPAlt,
  4×ApplicableQ predicates, MatrixLogxInverse, NullSpaceTryAgainOnFail,
  CombineDifferentialEquationsHomogeneous,
  CombineDifferentialEquationsWithPivotSelection, IntegrationStrategy,
  UseRationalRecurrence — 15.
- Logx/log-order reconstruction industry (exact p tags): NormalizeLogPower,
  MaxLogxPower, LogxPowerRange, SafeReplaceSeries11, UpdateIntReps, IntReps,
  DiffExpIntegrate1, IMaxLogOrder, IMaxLogOrderDefault — 9.
- Dies with the Prony/N-root fitter or the collapsed-exponent representation:
  FitResidualEndpointSectors, ApplyRegularizationStep, RegularizeIntegrand,
  PrintDecomposition — 4.
- sqrt-matrix machinery out of scope v1 (ledger v1.1): DEqnSquareRoots,
  CurrentSingularityWasAddedFromSquareRoot, SquareRootPrescriptionsAdded — 3.
- Silent-failure flags / post-hoc diagnostics (D1): MultivaluedFail,
  LastErrorContext, BenchmarkData — 3.
- Dynamic segmentation (v1 Predivision-only): GetMatricesPrecisionDistance — 1.
- Dead state / no consumers: AlphabetLogRulesFactored, AlphabetLogRulesExpanded,
  DiffExpExtensions, ClearMatrices — 4.
- Logging (waived config surface): LogFile, LogStream — 2.
- Generic one-line helpers with no API value: PrintMobiusNormalized, AllSameQ,
  CA, GetCases, DependsQ, ZeroQ, R, FirstOrNull, FindPivots, SplitTimes,
  SplitSum — 11.

---

## FT-referenced symbols: the M5 shim ledger

28 exported symbols are referenced by the FT layer (FeynmanTrick/ + Scripts/).
Site counts match the execution review's grep (State×27, Symbols×10,
Utilities×5, Transport×2, RegularizedIntegration×2,
SingularityDecomposition×1, plus DiffExp` top-level config calls×4).
Per the plan, M5 may not start until each has one of: a kept API name, a named
replacement accessor, or a scheduled call-site deletion.  Cross-reference:
M0 agent task (3) writes the full shim CONTRACT; this is the authoritative
symbol list for it.

A. Kept API — FT keeps calling (possibly under the new name):
1.  LoadConfiguration — DiffExpIntegration.m:367, 460
2.  UpdateConfiguration — DiffExpIntegration.m:411, 466
3.  TransportTo — DiffExpIntegration.m:440, 474
4.  EvaluateEndpointLimitSectors → EndpointLimit — DiffExpIntegration.m:1015
5.  DefiniteIntegralWithPrefactorLaurent → IntegrateOverLine —
    DiffExpIntegration.m:875; Scripts/eval_dump_generic.m:21
6.  MatrixDirectory — DiffExpIntegration.m:347
7.  ChopPrecision — DiffExpIntegration.m:349
8.  ExpansionOrder — DiffExpIntegration.m:350
9.  EpsilonOrder — DiffExpIntegration.m:351
10. UseMobius — DiffExpIntegration.m:352
11. UsePade — DiffExpIntegration.m:353
12. DivisionOrder — DiffExpIntegration.m:354
13. Verbosity — DiffExpIntegration.m:355
14. SegmentationStrategy — DiffExpIntegration.m:356
15. DeltaPrescriptions — DiffExpIntegration.m:412, 467
16. RationalizationTolerance — DiffExpIntegration.m:98, 109, 127, 154, 182,
    924 (semantics now Tol-derived; FT's hardcoded 10^-40 fallbacks at :154,
    :182 are deleted under A1)

B. Absorbed — FT gets a NAMED replacement accessor/API (shim work):
17. FEC → Cfg validated accessor — DiffExpIntegration.m:98, 109, 127, 154,
    182, 328 (load-guard `ValueQ[DiffExp`State`FEC]` needs a "is DiffExp2
    loaded" predicate), 372 (FEC[System`Variables])
18. FEWorkingPrecision → Cfg accessor — DiffExpIntegration.m:78
19. ExternalScalesVal → Cfg accessor (or delete the debug print) —
    DiffExpIntegration.m:418
20. MatricesIrreducibleFactors → API read-and-EXTEND of the singular-factor
    alphabet — DiffExpIntegration.m:211, 213 (WRITE), 402 (read).
    HIGHEST-RISK SHIM: FT currently mutates library state to make IBP-induced
    poles segmentable and prescribable.
21. DecomposeSingularity → SS named replacement (sector read of
    LocalSolution) — DiffExpIntegration.m:822
22. x (DiffExp`Symbols`x) → Cfg working-variable accessor —
    DiffExpIntegration.m:143, 153, 700, 831, 928;
    Scripts/check_transport_ode_residual.m:29
23. Logx → structured sector/evaluation access (no symbol scraping) —
    DiffExpIntegration.m:146; Scripts/check_transport_ode_residual.m:56, 58
24. \[Theta]p → branch argument to evaluate — DiffExpIntegration.m:64, 69;
    Scripts/check_transport_ode_residual.m:85, 86
25. \[Theta]m → branch argument — DiffExpIntegration.m:65, 70;
    Scripts/check_transport_ode_residual.m:85, 86
26. PChop → Tol chop — DiffExpIntegration.m:177, 705, 733, 777, 783;
    Scripts/check_transport_ode_residual.m:35

C. Dropped — FT call sites DELETED at M5:
27. IntegrationStrategy — DiffExpIntegration.m:358 (+ FT-level option
    plumbing :278, :303, :1670, :1687, :1805)
28. UseRationalRecurrence — DiffExpIntegration.m:357 (+ FT-level option
    plumbing :277, :302, :1669, :1686, :1804;
    Scripts/run_ft_stepwise.m:191)

## Non-symbol API surface FT also depends on (for the shim contract)

Not `::usage` exports, but breakages here are M5 blockers all the same:

- STRING config keys passed by FT: `"EstimateError"` (DiffExpIntegration.m:359),
  `"HomogeneousSolve"` (:360 — consumer is Wronskian.m:53,74 which is DROPPED;
  FT must stop passing it at M5), `"AbortOnAnalyticContinuationFail"`
  (:413, :468 — kept per §3.2 Config.m), `System`WorkingPrecision` (:348),
  `System`Variables` (read back at :372).  Old default key inventory:
  State.m:106-233.
- TransportTo RESULT keys consumed by FT/Scripts: `"SegmentData"`,
  `"NumIntegrals"`, `"EpsilonOrder"`, `"ErrorEstimates"` (produced at old
  Transport.m:506-507, 1238-1243; consumed e.g. DiffExpIntegration.m:917 and
  Scripts/check_transport_ode_residual.m:28).  DiffExp2 TransportTo must keep
  these names (or the shim maps them).
- Private-context reach-ins from Scripts (replace with public API or retarget
  the scripts at M6): `DiffExp`RegularizedIntegration`Private`laurentIntegralDump`
  (Scripts/eval_dump_generic.m:14), `...`Private`segmentActualBounds` (:34),
  `...`Private`IntegrateSegmentWithPrefactorLaurent` (:37),
  `...`Private`segmentMainExpression` (Scripts/check_transport_ode_residual.m:36),
  `...`Private`uncompressSeriesData` (:46); debug globals
  `DiffExp`State`$DebugFuchsianizedCheck` / `$DebugBlockResidualSeries`
  (Scripts/run_ft_stepwise.m:43, 45).
- Matrix-file symbol bridge: exported matrices spell the regulator as bare
  `eps` (MatrixExport.m:66 stripContexts); see the `eps` row.

## Open questions flagged for the orchestrator / spec agents

1. **Canonical dlog `d_1.m` input format** (ExpansionMatricesCanonical1,
   AlphabetLogs, AlphabetLogRules): live in
   Reference/Examples/FivePointNonPlanar_example.m:117 and
   Tests/test_five_point_nonplanar.m, but §3.2 LoadSystem names only
   full-format / closed-form / legacy slices.  API.md must either add the
   canonical format or declare the classic-pentagon example a
   legacy-parity-only configuration.
2. **IndefiniteIntegral**: tests use it (6 sites); §3.2 API.m lists only
   IntegrateOverLine.  Decide whether API exposes a piecewise indefinite form
   or M6 retargets those tests.
3. **IntegrateSystem surface name**: kept-in-API here; whether it maps to a
   TransportTo mode, SolveAtPoint, or keeps its own name is an API.md call.
4. **LogFile/LogStream waiver**: dropped here because §3.2's kept-key list
   omits LogFile, but the legacy review (finding 10) lists it among live keys.
   Config.md must record the waiver explicitly.
5. **sqrt tension**: §1 NON-GOALS drops sqrt x-dependence (loud error at
   LoadSystem) while §3.1 keeps "sqrt factors auto-prescribed as in old
   State.m DEqnSquareRoots".  Resolution assumed here: the auto-prescription
   RULE survives in Tr for prescription-supplied factors; the matrix-sqrt
   machinery (DEqnSquareRoots + flags) is dropped/v1.1.  Indicial.md/
   Transport.md specs must state this the same way.
6. **State.m export count**: 79 actual vs "~90" in the plan — no missing
   symbols; the plan's estimate included the System`-key and string-key
   config surface, which is covered by the Non-symbol section above.
