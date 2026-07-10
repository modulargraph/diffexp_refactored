# Module spec: SectorSeries.m (DiffExp2)

Status: M0 deliverable, v1, amended per Docs/specs/DECISIONS-M0.md and the
M0 reviews (REVIEW-math.md, REVIEW-minimalism.md).  Contract for the
implementation agent.  Read together with Docs/RewritePlan.md (sections 3.1,
3.2, 3.3) and the old-code citations below.  Budget: 400 + 75 Pade = 475
restated TOTAL, +50 review-added scope (2.10, 2.11) — section 9.

Module context: ``DiffExp2`SectorSeries` ``.  Position in the acyclic module
order: Tolerances < Config < EpsSeries < **SectorSeries** < Indicial < Solve <
Transport/Integrate < API.

---

## 1. PURPOSE

SectorSeries.m owns the fundamental data object of DiffExp2 — the
`LocalSolution` of RewritePlan 3.1, a vector-valued local series carrying exact
sector tags x^(a + b·eps)·(eps·Log t)^p/p! — and the closed algebra on it:
evaluate at a point (with the branch rule derived from the chart's
Prescriptions, accelerated by ported Pade approximants), multiply by a rational
function c(x,eps) (closed under partial fractions, with poles split across
charts), re-expand around a new regular center under an explicit truncation
contract, differentiate, linearly combine LocalSolutions with EpsSeries
weights (2.10, for Solve's compensation terms and Transport's matched
assembly), and parse boundary-condition tagged powers for API.m (2.11).  It
replaces the role of the old
DiffExp/SingularityDecomposition.m (435 lines) — which numerically *fitted*
sector exponents back out of collapsed per-eps-order towers — with exact tag
bookkeeping: tags are constructed once (by Indicial/Solve), never inferred, and
the named replacement API `SectorDecomposition` hands the exact sector list to
the retained FeynmanTrick layer.  It also absorbs the evaluation half of
DiffExp/Pade.m (GetPade/SEval/SEval1/SEval2) with that file's silent fallback
made loud.  This module makes NO tolerance-based structural decisions: numeric
smallness never changes a tag, a window, or a sector count.

---

## 2. PUBLIC SYMBOLS

All exported symbols MUST carry a `::usage` declaration before
`Begin["`Private`"]` (lesson L10, section 6: an unexported symbol called
fully-qualified from another package returns unevaluated and conditions fall
through silently — this killed the box family for days).

### 2.1 `ValidateLocalSolution[ls_Association] -> ls`

Checks every structural invariant of section 4 on `ls` and returns it
unchanged; LOUD error otherwise (section 5).  Called at the head of every
other public function in this module.  Cheap: shape checks, exactness checks
on tags, window-arithmetic consistency; no numerics.

### 2.2 `CanonicalizeLocalSolution[ls_Association] -> ls`

Returns the canonical form:
- sectors with identical exact tag `(a, b, p)` merged by adding `Coeffs`;
- sectors with equal `(b, p)` and integer-spaced `a` are considered greedily
  in ascending `a`.  To merge a higher-`a` slab into a lower-`a` slab, let
  `s = aHigh - aLow >= 0`: its `n` columns are shifted right by `s` inside
  the existing fixed Taylor width.  The merge is performed only when every
  coefficient that truncation would discard (the final `s` columns, or the
  whole slab when `s >= ncols`) is syntactically zero.  Otherwise the towers
  remain separate.  Thus canonicalization never erases known finite-`TOrder`
  data merely to eliminate an integer-spaced tag pair;
- sectors whose `Coeffs` are SYNTACTICALLY all zero are dropped (exact zeros
  only — tolerance-based dropping is forbidden, F10 in section 5.3), except
  that an all-zero `b = 0`, `a < 0` tower is retained when
  `Floor[-a] >= ncols`: its stored zero prefix does not certify the unseen
  endpoint-relevant tail.  If every sector would otherwise be dropped, one
  zero representative is retained (in particular the regular `(0,0,0)`
  sector is never lost);
- sectors sorted lexicographically by `(a, b, p)` (canonical order shared with
  Integrate.m and the FT shim).

### 2.3 `ChartImSign[ls_Association] -> +1 | -1 | None`

Derives the chart's Im-sign from `ls["Prescriptions"]` (the LIST of
`<|"Factor","Sign","Multiplicity","LeadingCoeffSign"|>` records, RewritePlan
3.1).  Rule, port of old DiffExp/AnalyticContinuation.m:18-68:

- entries with EVEN `"Multiplicity"` contribute no constraint
  (AnalyticContinuation.m:46-50 marked these "?"; the old code then treated
  any "?" as failure — the new rule per RewritePlan 3.1 is: even multiplicity
  = no constraint);
- each ODD-multiplicity entry requires `sigma = Sign * LeadingCoeffSign`
  (old form: prescribed sign divided by the sign of the leading coefficient,
  AnalyticContinuation.m:45-53; division and multiplication coincide on ±1);
- all required sigmas must agree -> return that sigma;
- no constraining entries -> return `None`;
- conflict -> LOUD error `SectorSeries::branchconflict` (section 5.1).  The
  old behavior — set `AnalyticContinuationFailed = True` and silently continue
  with `uniqueSigns = {1}` (AnalyticContinuation.m:55-68) — is FORBIDDEN (F8).

Deterministic and cheap; called by `EvaluateLocalSolution` and
`ReexpandLocalSolution` whenever a branch is needed.

SINGLE OWNER (math review D28; minimalism review 18): `ChartImSign` is the
ONE sign-derivation in DiffExp2.  Transport's `Chart["CrossSign"]` and
Integrate's pairing/negative-arm sigma are DEFINED as this function applied
to the chart's Prescriptions record; neither module reimplements the rule.
SIGN-AWARENESS (DEC-16): the derivation is sign-aware through
`LeadingCoeffSign` — two prescription entries whose factors differ only by
leading-coefficient normalization (`{-1+x, +1}` vs `{1-x, +1}` flip the
implied i-delta side) derive OPPOSITE sigmas and MUST surface as
`::branchconflict`; a sign-blind upstream dedup (the old
deltaPrescriptionsForFactors class) that would hide the conflict is
forbidden.

### 2.4 `EvaluateLocalSolution[ls_Association, tval, opts] -> EvaluationResult`

Evaluates the LocalSolution at the chart-coordinate point `tval`.

The multi-point form is CUT (taken up front per REVIEW-minimalism defect 4,
cut 3): callers needing several points (the old pattern: `SEval1` once,
`SEval2` per point — Transport.m:905-918 evaluates the same segment at the
eval point, the boundary fix point, and two reduced-order probe points) loop
over single-point calls; the cost of rebuilding Pade per point is ACCEPTED
and recorded in the R6 benchmark notes.

Arguments:
- `tval`: a real number (exact rational, or arbitrary-precision Real).
  Complex or symbolic points -> LOUD error (v1; open question
  Q5).  `0 < |tval| < ls["Radius"]` required for singular charts;
  `|tval| < Radius` (including 0) for regular charts.
- Options:
  - `"UsePade" -> Automatic | True | False`.  Automatic reads Config key
    `UsePade`.  See 2.4.2.
  - `"TOrderReduction" -> d_Integer >= 0` (default 0): evaluate with the top
    `d` t-orders of every sector's array dropped (exact truncation, no
    re-validation of TWindow).  This is the seam for Transport's two-point
    full-vs-reduced error probe (old `DecreaseSeriesOrderBy`,
    DiffExp/SeriesOps.m:196-203, used at Transport.m:913).
  - `"ImSign" -> Automatic | 1 | -1` (default Automatic = `ChartImSign[ls]`):
    test/Transport override of the branch sign.

Return value `EvaluationResult` (section 3.4): the eps-Laurent vector value
as a canonical EpsSeries.m object (DEC-13) with sibling evaluation metadata
(`"PadeFallbacks"`, `"TailEstimates"`).

#### 2.4.1 Evaluation formula and branch rule

Let `wp = Config[WorkingPrecision]`.  All numerics run inside
`Block[{$MinPrecision = wp}, ...]` and the result is `SetPrecision[..., wp]`
(port of Pade.m:34, 70-77, 80-88 — lesson L5).

For `tval > 0`: `L = Log[tval]` (real), `t^a = tval^a` (real positive root for
fractional a).

For `tval < 0`: requires `sigma = ChartImSign` (or the `"ImSign"` override)
whenever any sector has non-integer `a`, `b != 0`, or `p > 0`; if no such
sector exists the branch is irrelevant and `sigma` is not consulted.  Then

    L            = Log[-tval] + I Pi sigma
    tval^(a+b eps) = (-tval)^(a + b eps) * Exp[I Pi sigma (a + b eps)]

This single rule subsumes ALL the old per-case replacements
(AnalyticContinuation.m:70-79): old sign +1 applied no replacement at all
(principal branch at negative argument — `Log[t] = Log[-t] + I Pi`,
`t^a = |t|^a e^{I Pi a}` — exactly sigma = +1 above); old sign −1 applied
`Logx -> Logx - 2 Pi I θm` plus `x^b -> Exp[-2 Pi I b] x^b` (denominator > 2)
and the `Sqrt` `(θp - θm)` half-integer case — all equal to principal times
`e^{-2 Pi I (a + b eps)}`, i.e. sigma = −1 above.  State in code comments that
sigma = +1 is the principal branch (no replacement), pinning the convention
per the 9aeb300 interior-split lesson and RewritePlan §5.  Coefficients are
theta-free by invariant I-9; no `Projectθs`-style projection is ever needed
at evaluation.

Per sector `(a, b, p)` with rows `c[k, n, comp]`, the contribution to the
value's eps-order `K` coefficient (component vector) is

    tval^a * L^p / p! * Sum[ alpha[k] * (b L)^j / j! ,
                             {j, 0, K - p - kmin}, k = K - p - j ]

where `alpha[k][comp]` is the evaluation of the inner t-series
`Sum[c[k,n,comp] tval^n, {n, 0, nmax}]` (direct sum, or Pade per 2.4.2), and
`tval^a`, `L` are branch-resolved as above.  This is the exact-algebra version
of the old `x^(b eps) = Exp[b eps Logx]` convolution
(SingularityDecomposition.m:222-243, `MultiplyByXMinusBEps`) — same math, now
applied forward instead of being fitted backward.

Value window (the returned EpsSeries' `EpsWindow`, section 3.4):
`"Min" = Min over sectors of (first present row + p)`;
`"CompleteMax" = ls["EpsWindow"]["CompleteMax"]` (every value order
`K <= CompleteMax` needs only rows `k = K - p - j <= K <= CompleteMax`, all
present).  Requesting orders above CompleteMax is impossible by construction
(the return object carries the window; consumers check — RewritePlan 3.1
invariants).

TailEstimates on EVALUATION results (minimalism review 29): per eps order,
`Max over comp of |c[k, ntop, comp] * tval^ntop| * q/(1-q)` with
`q = |tval|/Radius` and `ntop` the last used t-column — the geometric tail
at the actual evaluation point.  Computed always (cheap); Missing[] never
occurs on the evaluation path (the §3.4 Missing[] option applies only to
objects that bypassed evaluation).  Transport's probe (its 2.10) remains the
authoritative error feed; TailEstimates is advisory metadata.

Evaluation at `tval == 0`:
- regular chart (single `(0,0,0)` sector): exact, returns the `n = 0` column.
- any sector with `a < 0`, fractional `a`, `b != 0`, or `p > 0` present:
  LOUD error `SectorSeries::originlimit` directing the caller to
  `Integrate`EndpointSectorLimit` (Integrate.md 2.4, per REVIEW-minimalism
  defect 12) under API.m's EndpointLimit (the b≠0 drop rule and the
  divergence gate live THERE, not here; RewritePlan 3.3).  Silently
  returning the `(0,0,0)` constant would re-create the old wholesale-drop
  disease (banana doc, Docs/FeynmanTrickBananaStatus.md:314-329).

#### 2.4.2 Pade acceleration (ported, loud)

Why it exists (legacy review finding 11): matching and evaluation points sit
at `radius/DivisionOrder` of the chart (FT pins DivisionOrder = 4,
FeynmanTrick/DiffExpIntegration.m:272; classic default 3; the GetCPL/GetCPR
geometry, DiffExp/Mobius.m:93-142, places the point at the 1/k fraction of
BOTH adjacent charts), where plain partial sums converge slowly.  Pade is
accuracy per order, not luxury; the R6/M4 benchmark compares old-with-Pade vs
new-with-Pade.  Classic examples REQUIRE `UsePade -> True`
(Reference/Examples/Banana_example.m, FivePointNonPlanar_example.m).

Port of old GetPade (DiffExp/Pade.m:30-53), restructured: the old code applied
PadeApproximant per-Logx-coefficient of each per-eps-order SeriesData
(Pade.m:40-52, using `MaxLogxPower`/`LogxCoeff` from SeriesOps.m:108-136).  In
the new representation the log structure is the exact tag `p`, so Pade applies
to the inner t-polynomial per `(sector, eps-row k, component)`:

- diagonal order `[m/m]` with `m = Floor[numberOfCoefficients/2]` (so
  2m+1 <= numberOfCoefficients: the approximant is determined ENTIRELY by
  known coefficients; for `{c0..cN}`, `m = Floor[(N+1)/2]`).  DELIBERATE
  one-slot deviation from old Pade.m:35, whose SeriesData arithmetic
  (`(a[[5]]-a[[4]])/a[[6]] + 1`) overcounts by one and silently consumes a
  phantom top order — LessonsLedger entry required ("Pade order off-by-one:
  old code built [m/m] needing 2m+1 coefficients from 2m known ones,
  treating the truncation boundary as exact zeros").  Normative per DEC-12;
  unit t14 is the gate;
- coefficients chopped before construction at the Tolerances-derived
  `chopFloor` (old: `Chop[#, 10^-ChopPrecisionVal]`, Pade.m:43 — the constant
  moves to Tolerances.m);
- constructed once per (ls, options) call (the multi-point form is cut per
  REVIEW-minimalism defect 4 — callers evaluating several points rebuild;
  cost accepted, see 2.4);
- FAILURE (Mathematica's `PadeApproximant` returns unevaluated — detected
  exactly as the old code does, by an embedded `PadeApproximant[__]` head,
  Pade.m:44): the old code printed a level-1 warning and silently evaluated
  normally (Pade.m:44-48).  NEW CONTRACT: emit the loud, never-Quiet-ed
  message `SectorSeries::padefail` naming chart center, sector tag `(a,b,p)`,
  eps row `k`, and component index; fall back to the direct partial sum FOR
  THAT SERIES ONLY; and record the event in the returned value's
  `"PadeFallbacks"` list so Transport's error probe and the user-facing
  ErrorEstimate see it.  This is the ONLY permitted fallback in the module
  (section 5.3, F1), permitted because it is loud AND recorded AND
  mathematically safe (the partial sum is exact data, merely slower to
  converge).
- coefficients containing symbolic indeterminates with Pade requested ->
  LOUD error `SectorSeries::indetpade`.  The old code silently mutated the
  global config (`TurnOffPade`, Transport.m:178-180, 566-568); the
  disable-Pade-for-indeterminate-runs DECISION belongs to Transport/API and
  must be made there, loudly, once.

Out-of-window embedded `SeriesCoefficient` tails — the old
`SafeSeriesCoefficient`/`NormalizeEmbeddedSeries` machinery (Pade.m:16-27,
60-66) existed because upstream series arithmetic left inert one-past-the-end
`SeriesCoefficient[...]` heads in coefficients (box L0 fix, item 1 of
Docs/FeynmanTrickBoxFamilyStatus.md "RESOLVED: the eps^0 deficit").  In
DiffExp2 coefficients are plain array entries (numbers or linear forms in
declared indeterminates, invariant I-8); a symbolic head in `Coeffs` is a
construction-time LOUD error, so this failure class is excluded by
representation and NO normalization shim is ported.

### 2.5 `MultiplyRational[ls_Association, c_, opts] -> ls'`

Multiplies the LocalSolution by a rational function `c` of (variable, eps).
The algebra is CLOSED: the result is again a LocalSolution on the same chart.

Arguments:
- `c`: an expression rational in the CHART coordinate symbol and the eps
  symbol (eps symbol per Config's `CanonicalEps[]`).  Anything non-rational
  (Sqrt of the variable, Log, special functions) -> LOUD error
  `SectorSeries::nonrational` (RewritePlan non-goal: algebraic x-dependence
  is out of scope v1).
- The `"Coordinates" -> "Main"` conversion option is CUT (taken up front per
  REVIEW-minimalism defect 4, cut 1): the FT shim converts main-line-variable
  IBP coefficients through `ls["ChartMap"]` BEFORE calling (old conversion
  site FeynmanTrick/DiffExpIntegration.m:731-736 maps to the shim); t23
  retargets the shim's test file.

Algorithm (normative):

1. `c = Together[c]`, numerator/denominator exact.  If the denominator is
   divisible by `eps^r`, factor `eps^r` out first: this is a pure eps-Laurent
   SHIFT of the result windows (`MinPower` and `CompleteMax` both drop by
   `r`), per the EpsSeries.m Laurent-division contract (RewritePlan 3.2,
   EpsSeries: "the 1/(b eps) enhancement is a shift, not an inversion").
   This is NOT an error.
2. Expand `c` as an eps-Laurent series with rational-in-t coefficients
   `c_j(t)`, `j` from the shifted minimum up to what the target
   `CompleteMax` requires (EpsSeries.m arithmetic; honest window
   intersection).
3. Solve the t-denominator's roots EXACTLY at eps = 0 (the poles of every
   `c_j(t)` lie at roots of `Denominator[c](t, 0)`; eps-dependent root motion
   only raises multiplicities order by order in eps).  Exact algebraic roots
   (Root objects allowed).  Classify each root `t_i` (this is the
   "partial fractions split poles ACROSS charts" contract, RewritePlan 3.2,
   and math-review finding 8):
   - `t_i == 0` (the chart center): the pole shifts sector exponents.
     Implemented by writing `c_j(t) = t^(-M) q_j(t)` with `q_j` analytic at 0
     and `M = ` max center-pole order over the needed eps range; every sector
     `(a, b, p)` maps to `(a - M, b, p)`.
   - `|t_i| >= Radius`: the pole belongs to ANOTHER chart; its partial
     fraction is analytic on this chart's disk and is folded into the regular
     Taylor part of `q_j` (coefficients exact rationals/algebraics, computed
     in closed form, truncated at the t-order needed — the truncation of an
     EXACT geometric tail, complete to the kept order by construction).
   - `0 < |t_i| < Radius`: LOUD error `SectorSeries::interiorpole`
     (section 5.1).  NO radius-shrink fallback (F6): Transport.m's
     segmentation (matrix denominators plus ExtraSingularFactors, cf.
     CollectLevelIBPSingularFactors) is REQUIRED to make this branch
     unreachable by placing a chart at every pole; reaching it is a
     segmentation bug, not a representable state.
   - `|t_i|` vs `Radius` comparison numerically ambiguous at WorkingPrecision
     — difference magnitude below `Tol["GeomGuardTol"] * Max[|t_i|, Radius]`
     (the named Tolerances.m export `GeomGuardTol[wp] = 10^(-Floor[wp/2])`,
     per REVIEW-minimalism defect 17 / math review D26): LOUD error
     `SectorSeries::geomambiguous` — never pick a side silently.  The same
     guard governs the `|D|` vs `Radius` comparison of 2.6.
4. Convolve: `c'[k, n, comp] = Sum[ Q_jrow[n1] * c[krow, n2, comp] ]` over
   `jrow + krow = k`, `n1 + n2 = n`, where `Q` is the assembled analytic part.
   eps windows combine per EpsSeries.m multiplication (honest min);
   t windows: the array index range `0..nmax` is preserved (Q is exact to any
   order), TWindow metadata carried over unchanged.
5. Canonicalize (2.2) and return.  Exact cancellations may produce syntactic
   zero rows — droppable; NUMERICALLY small leading entries are NOT droppable
   here (F10) — object-level cancellation checks belong to Integrate.m
   (RewritePlan 3.2, Integrate: "checked at the object level").

Note the post-multiply object at an ODE-regular chart legitimately has
sectors with `a != 0`; the regular-chart one-sector invariant is scoped to
ODE solutions and is NOT re-asserted here (math finding 8; invariant I-7).

### 2.6 `ReexpandLocalSolution[ls_Association, newCenter_, targetOrder_Integer, opts] -> ls'`

Re-expands the LocalSolution around a new center inside the chart, producing
a REGULAR LocalSolution (exactly one `(0,0,0)` sector).

Arguments:
- `newCenter`: real, in CHART coordinates; `0 < |newCenter| < Radius` for a
  singular chart, `|newCenter| < Radius` and `newCenter != 0` for regular
  ones (re-expanding around the same center is the identity; allowed but
  must round-trip exactly).  Re-expansion AROUND A SINGULAR POINT (i.e. with
  the new center placed at a singularity, or any request to produce tagged
  sectors at the new center) is FORBIDDEN: tags at singular charts come ONLY
  from fresh Frobenius (Indicial.m + Solve.m).  Violations ->
  `SectorSeries::resingular`.
- `targetOrder`: the t'-Taylor order N of the new expansion (t' = t −
  newCenter).
- Options: `"ImSign"` as in 2.4 (needed when `newCenter < 0`).

Mathematics (normative), per sector `(a, b, p)` with `D = newCenter`,
`w = t'/D`:

    t^(a + b eps)     = D^a * Exp[b eps Log[D]] * (1 + w)^(a + b eps)
    (1 + w)^(a+b eps) = Sum[Binomial[a + b eps, m] w^m, {m, 0, N}]   (exact,
                        Binomial expands to an eps-polynomial of degree m)
    (eps Log t)^p/p!  = eps^p (Log[D] + Log[1 + w])^p / p!           (Log[1+w]
                        as the exact alternating series to order N)
    t^n               = D^n (1 + w)^n

with `Log[D]`, `D^a` branch-resolved by the sigma rule of 2.4.1 when `D < 0`.
Collect into eps-orders and t'-powers; the output rows of the single
`(0,0,0)` sector ARE value eps-orders.

TRUNCATION CONTRACT (math-review finding 12 — this is the load-bearing part):
- target order: the caller's `targetOrder` N; the output array has exactly
  N+1 columns;
- tail bound: the new chart's honest radius is
  `rho = Min[R_old - |D|, If[chart singular, |D|, Infinity]]`
  (triangle inequality; the old center is itself a singularity at distance
  |D| when the chart is singular).  For each output eps order k the module
  records (math review D11 — the formula must be rho-aware; the unit-radius
  form under- or over-estimates by rho^(±N) on small/large charts):

      TailEstimate[k] = Max over comp of Abs[ c'[k, N, comp] ] * r^N * q/(1-q),
      with r = rho/DivisionOrder (the design evaluation radius) and
      q = 1/DivisionOrder (so q/(1-q) = 1/(DivisionOrder - 1)).

  Equivalently, in the radius-normalized coordinate w = t'/rho the familiar
  form |ĉ_N| q^(N+1)/(1-q) holds with ĉ_N = c'_N rho^N; the implementation
  must use whichever form matches its stored coordinate.  This is the
  geometric tail at the DESIGN evaluation radius rho/DivisionOrder
  (DivisionOrder from Config; FT defaults to 3).  The geometric model (|c'_n| ~
  rho^-n, top kept coefficient representative) is valid because RoC chart
  rescaling keeps high-order coefficients O(1) (ledger lesson; banana
  requires RoC = 10).  These estimates are ADDED per eps order into
  `ls'["ErrorEstimate"]` (never replacing — ErrorEstimate is additive across
  operations and segments, RewritePlan 3.1), and returned in the
  `"TailEstimates"` key.  This is the CurrEvalError-equivalent feed named in
  RewritePlan §5.
- eps windows: `CompleteMax' = CompleteMax` (every output row K assembles
  from input rows `K - p - j <= K`, all present);
  `MinPower' = Min over sectors (first present row + p)`.

STATED EXPLICITLY (required by the plan): re-expansion COLLAPSES the exact
b-structure into per-eps-order logs — `Exp[b eps Log D]` becomes
`Log[D]^j/j!` constants inside eps-order coefficients, and the tags are gone
from the new object.  This is ACCEPTABLE for the pipeline because tags are
never propagated through regular regions: at the next singular chart the
sector spectrum is re-derived by fresh Frobenius from the exact indicial data
(Indicial.m contract I1), not reconstructed from these coefficients.  Any
attempt to reconstruct tags from a re-expanded object is the old D2 disease
and is forbidden.

New object fields: `Center = old Center + chart-map image of D` (ChartMap
composed with the affine shift `t = D + t'`), `Radius = rho`,
`Prescriptions = {}` (regular center — no vanishing prescription factors),
`Sectors = {single (0,0,0) sector}`, windows and ErrorEstimate per above.

### 2.7 `DifferentiateLocalSolution[ls_Association] -> ls'`

Exact termwise derivative, in CHART coordinates.

Chart-coordinate rule (normative).  Each sector `(a, b, p)` with rows
`c[k, n, comp]` contributes:

- to sector `(a - 1, b, p)`:  rows `(a + n) c[k, n, comp] + b c[k-1, n, comp]`
  (the eps-linear factor `(a + n + b eps)` acts as an exact window-respecting
  shift-and-add: output row k consumes input rows k and k−1, so CompleteMax
  is preserved);
- to sector `(a - 1, b, p - 1)` when `p >= 1`:  rows `c[k-1, n, comp]`
  (from `d/dt (eps Log t)^p/p! = (eps Log t)^(p-1)/(p-1)! * eps/t`; the row
  shift +1 maps the resonant exemption window `kmin = -p` exactly onto the
  target sector's `kmin = -(p-1)` — no window violation by construction).

Output sectors are canonicalized (2.2).  Value-level t-completeness drops by
one absolute power (the derivative of a function known mod `t^(a+N+1)` is
known mod `t^(a+N)`); since the output exponent is `a - 1` the INDEX range
`0..nmax` and the TWindow metadata are unchanged.  This replaces the old
per-Logx-power derivative `SD` (DiffExp/SeriesOps.m:206-214) with exact tag
algebra.

The `"Coordinates" -> "Main"` option is CUT (taken up front per
REVIEW-minimalism defect 4, cut 1): the caller (FT shim / Transport) applies
the exact affine factor `1/beta` of `x = x0 + beta*t` itself; t29 retargets
the shim's test file.  ChartMap is AFFINE ONLY in the new core (DEC-18:
Mobius is dropped entirely; RoC chart rescaling is an affine rescaling and
is kept) — no Mobius chain rule exists.

### 2.8 `SectorDecomposition[ls_Association] -> Association`

THE NAMED REPLACEMENT API for the old
``DiffExp`SingularityDecomposition`DecomposeSingularity`` (live FT call site:
FeynmanTrick/DiffExpIntegration.m:822; shim contract doc lists it as the one
SingularityDecomposition reach-in).  Returns

    <| "Sectors"   -> { <|"a" -> a, "b" -> b, "p" -> p, "Coeffs" -> c|> .. },
       "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
       "TWindow"   -> <|"CompleteMax" -> nmax|>,
       "Center"    -> center, "ChartMap" -> map, "Radius" -> radius |>

with `Sectors` the canonicalized (2.2) exact sector list — i.e. the
LocalSolution's own data, validated and sorted; NO computation, NO fitting,
NO tolerance.  Mapping to the old output shape `{<|"a","b","g"|>..}`
(SingularityDecomposition.m:25, Docs/SingularityDecomposition.md): old `g` =
per-eps-order series of one collapsed term; new `Coeffs` = the exact
per-sector array with `p` explicit.  Differences the FT layer must absorb
(documented here because the shim agent's contract consumes this section):
several sectors instead of one collapsed term; `p` explicit instead of
residual `Logx` towers; vector-valued (component index); exact tags instead
of `RationalizationTolerance`-rationalized fits.  The drop rule consumers
(limitUpper/limitLower: "constant of the (0,0,0)-sector; b != 0 dropped
exactly", RewritePlan 3.3) read this output through
`Integrate`EndpointSectorLimit` (Integrate.md 2.4) under API.m's
EndpointLimit; the FT
verbosity printout at DiffExpIntegration.m:822-826 maps to
`{#a, #b, #p} & /@ result["Sectors"]`.

### 2.9 `padeEvaluate` (Private — NOT exported)

Internal workhorse, demoted from the public surface (cut 5 taken up front
per REVIEW-minimalism defect 4; t14/t15 exercise it through
EvaluateLocalSolution): given one plain coefficient vector
`{c0, c1, ..., cN}` (numbers, arbitrary precision), build the `[m/m]`
diagonal Pade approximant with `m = Floor[numberOfCoefficients/2]
= Floor[(N + 1)/2]` (the 2.4.2 formula, normative per DEC-12 — NOT old
Pade.m:35's off-by-one) of `Sum[c_n t^n]` and evaluate at the point.  Same
chop, precision-block, failure-detection, and loud-fallback contract as
2.4.2 (failure here raises `SectorSeries::padefail` with sector/order
context supplied by the caller via an internal hook).

### 2.10 `CombineLocalSolutions[weights_List, lss_List] -> LocalSolution`

(Added per REVIEW-minimalism defect 6 — the operation Solve.md §3.7 step 2
and Transport.md §2.6 step 4 consume.)  Exact linear combination
Σ_i w_i · ls_i where each `w_i` is an EpsSeries (EpsSeries.md object; a
plain exact number is auto-wrapped as a width-1 series with window
[0, +inf-equivalent: the partner's CompleteMax]).  All ls_i must share
Center/ChartMap/Radius/Prescriptions exactly (else `::dims`).  Per sector:
coefficient rows are EpsSeries-multiplied by w_i (window per ESTimes;
Laurent weights shift kmin honestly — this is how Solve's compensation terms
with polar γ(eps) and Transport's matched weights enter), same-tag sectors
merged by ESAdd, result canonicalized (2.2).  EpsWindow = honest min over
contributions; TWindow = min; ErrorEstimate = entrywise sum (I-10).  Budget
+~25 lines (funded per REVIEW-minimalism defect 4).  Unit test t34.

### 2.11 `ParseTaggedPower[expr_, var_Symbol, epsSym_Symbol] -> <|"a", "b", "p", "Coefficient"|> | $FailedParse`

(Added per REVIEW-minimalism defect 13.)  Boundary-ingestion helper (sole
caller: API.m PrepareBoundaryConditions).  Recognizes
c * var^(a + b*epsSym) * Log[var]^p products with exact a, b (affine-in-eps
exponent REQUIRED: a non-affine exponent is a LOUD error naming the exponent
— same I1 discipline as Indicial, never a numeric fit), p integer >= 0;
rewrites Log[k*var] -> Log[k] + Log[var] first (old Transport.m:120 rule).
Returns the inert marker $FailedParse for expressions that are not of this
product form — the CALLER (API.m) decides between the documented minimum
contract (old eps-expansion into Logx polynomials, old Transport.m:91-120)
and an error; $FailedParse is data for that documented branch, not a silent
fallback.  Budget +~25 lines (funded per REVIEW-minimalism defect 4).  Unit
test t35.

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution (RewritePlan 3.1 as amended by DECISIONS-M0
DEC-9/DEC-16/DEC-18 — normative)

    LocalSolution = <|
      "Center" -> exact x0, "ChartMap" -> AFFINE ONLY (DEC-18: Mobius is
                  dropped from the new core entirely; RoC chart rescaling is
                  an affine rescaling and is kept),
      "Radius" -> distance to nearest true singularity expressed in the
                  local affine chart coordinate t (physical distance divided
                  by ChartMap["Scale"]; restored Re, Re±Im points determine
                  the affine scale but are regular planner waypoints),
      "Sectors" -> { Sector.. },
      "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
      "TWindow"   -> <|"CompleteMax" -> nmax|>,   (* the truncation-order
                  record ONLY (DEC-9): recursion matrices are exact
                  polynomials, so there is NO coupling-depth degradation in
                  the new core — the old MaxCouplingOrder discount was a
                  numeric-matrix artifact, subsumed by exact recursion +
                  ErrorEstimate *)
      "ErrorEstimate" -> per (eps-order) accumulated error (two-point
                  full-vs-reduced-order probe, additive across segments,
                  abort > 1 — ported old machinery, user-facing),
      "Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
                  "LeadingCoeffSign"|> with the derived chart Im-sign;
                  consistency-checked at construction (even multiplicity = no
                  constraint; conflict or missing prescription at a chart that
                  is multivalued AT ALL — b != 0 OR p > 0 OR Denominator[a] >
                  1, DEC-16: the Kallen charts are fractional-a with b=0, p=0
                  — = LOUD ERROR; dedup of prescription factors is SIGN-AWARE
                  per DEC-16 and 2.3; the old sqrt auto-prescription union has
                  no v1 source — LoadSystem rejects sqrt matrices, API.md E6)
    |>

    Sector = <|
      "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
      "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT:
                  solutions are vectors; matching and ODE residuals need all
                  components *)
    |>
    Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                    Σ_k eps^k Σ_n c[k,n,comp] t^n

### 3.2 Index conventions (this spec's binding refinement of 3.1)

- `Coeffs` is a full (dense) 3-D array.  First index: eps COEFFICIENT rows of
  the inner sum, row `i` ↔ eps power `EpsWindow["Min"] + i - 1`; ALL sectors
  share the LocalSolution-level row range `[Min, CompleteMax]`.  Rows below a
  sector's intrinsic support are EXACT zeros — this is not padding in the
  forbidden sense (padding = claiming completeness for absent data; an exact
  zero is data).  Rows above `CompleteMax` are ABSENT, never zero-filled.
- Resonant exemption (RewritePlan 3.1 invariants, math finding 14):
  homogeneous true-resonance sectors have row support starting at `-p` BY
  CONSTRUCTION (log-chain weights are eps-independent under the
  `(eps Logx)^p/p!` normalization, so the VALUE starts at eps^0).  Window
  arithmetic must account for `Min <= -p`, not error.
- Second index: t-powers `n = 0, 1, ..` (column `i` ↔ `t^(i-1)`).  `n` is a
  NON-NEGATIVE INTEGER within every sector; all fractional and negative
  leading behavior lives in the exact tag `a`.  Sectors of one Frobenius
  family step in integers; distinct fractional offsets are distinct sectors.
- Third index: component, `1..Ncomp`, identical across sectors (math finding
  10: solutions are VECTORS; matching solves and ODE residuals need all
  components).
- `TWindow["CompleteMax"]` is the truncation-order record ONLY (DEC-9): the
  exact polynomial recursion delivers every stored column complete, and no
  coupling-depth degradation exists in the new core.  Columns beyond
  `TWindow["CompleteMax"]` do not arise from the solver; consumers treat
  CompleteMax as the completeness bound.  Evaluation honesty is supplied by
  the two-point error probe via `"TOrderReduction"`.
- Coefficient entries: exact numbers, arbitrary-precision numbers, or linear
  combinations of DECLARED indeterminate symbols (the API.m symbolic-BC
  feature).  `SeriesData`, `SeriesCoefficient`, theta symbols, or any other
  compound head -> construction-time LOUD error.
- LocalSolution carries an OPTIONAL `"Indeterminates" -> {sym..}` key
  (default {}; populated by API.m from the BoundaryConditions record).
  I-5's "declared indeterminate symbols" means EXACTLY this list: a symbolic
  coefficient entry whose symbols are not all in the object's own
  Indeterminates list is `::badcoeff`.  Config holds no indeterminate state
  (REVIEW-minimalism defect 16).

### 3.3 EpsWindow / TWindow (verbatim semantics)

`EpsWindow = <|"Min" -> kmin, "CompleteMax" -> kmax|>`: rows outside are
absent; consumers check `need <= CompleteMax` and fail loudly naming (chart,
sector, order).  `TWindow = <|"CompleteMax" -> nmax|>`: t-order honesty as in
3.2.  No SectorSeries operation ever widens a window without exact
justification (the enumerated cases: eps-Laurent shift in 2.5 step 1 moves
both ends; everything else is min-arithmetic).

### 3.4 EvaluationResult (returned by EvaluateLocalSolution)

(Rewritten per math review D18 / REVIEW-minimalism defect 23 / DEC-13.)

    <| "Value"         -> the canonical EpsSeries.m object
                          <|"EpsWindow" -> <|"Min" -> m, "CompleteMax" -> M|>,
                            "Coeffs" -> { v_m, ..., v_M }|>
                          (key names VERBATIM from Docs/specs/EpsSeries.md
                          §3.1 — the ONLY in-core eps-Laurent shape, DEC-13;
                          each coefficient v is a length-Ncomp numeric List),
       "PadeFallbacks" -> { <|"Sector" -> {a,b,p}, "EpsRow" -> k,
                             "Component" -> comp|> .. },
       "TailEstimates" -> per-eps-order list (always computed on the
                          evaluation path, 2.4.1; Missing[] only on objects
                          that bypassed evaluation) |>

Evaluation metadata is a SIBLING of the series value, never carried inside
the EpsSeries object.  API.md's user-facing LaurentValue
(`"MinPower"/"Coefficients"/"CompleteMax"`) is the FT-boundary shape only,
produced exclusively by an API.m converter
(`ESMinPower`/`ESCoefficientList`/`ESCompleteMax`); no core module
constructs or consumes it.

### 3.5 MultiplyRational input

`c`: Mathematica expression, rational in (coordinate symbol, eps symbol) with
exact (rational/algebraic) coefficients; arbitrary-precision numeric
coefficients accepted (FT reality) but the DENOMINATOR structure must be
exact — an inexact denominator root makes the interior-pole classification
tolerance-dependent and raises `SectorSeries::inexactdenominator`.

---

## 4. INVARIANTS (always-on, cheap; checked by ValidateLocalSolution and re-asserted by every op on its output)

- I-1 Tags exact: `a`, `b` carry no inexact numbers (`Precision === Infinity`
  on every numeric part); `p` a non-negative machine integer.
- I-2 Shape: every sector's `Coeffs` has dimensions
  `{kmax - kmin + 1, ncols, Ncomp}` with the SAME row range and `Ncomp`
  across sectors; `ncols >= TWindow["CompleteMax"] + 1`.
- I-3 Windows: `kmin <= kmax`, integers; `kmin` may be negative (resonant
  exemption `kmin = -p`); `TWindow["CompleteMax"] >= 0`.
- I-4 `n` columns are integer-indexed from 0; no fractional steps (fractional
  structure in `a` only).
- I-5 Coefficient entries per 3.2 (numbers or declared-indeterminate linear
  forms; no compound heads, no theta symbols).
- I-6 Lossless canonical tag set: after `CanonicalizeLocalSolution`, exact
  `(a,b,p)` tags are pairwise distinct.  A same-`(b,p)` integer-spaced pair
  may remain when the greedy fixed-width shift was rejected because it would
  discard a syntactically nonzero known tail; such a pair is canonical and
  downstream consumers must accept it.  Under-expanded all-zero negative
  `b=0` towers are structural completeness records, not droppable zero data.
- I-7 Regular chart ⇔ exactly one `(0,0,0)` sector — asserted FOR ODE-SOLUTION
  OBJECTS at construction (by Solve.m/Transport.m), NOT re-asserted after
  `MultiplyRational` (math finding 8: products with IBP coefficients
  legitimately violate it).  SectorSeries itself only asserts the implication
  "single (0,0,0) sector ⇒ branch-free evaluation".
- I-8 `Radius > 0`, real; every evaluation/re-expansion point strictly inside.
- I-9 Coefficients are branch-resolved: free of theta-like symbols; the
  branch enters ONLY through the sigma rule at evaluation/re-expansion and
  through Transport's crossing operator.
- I-10 ErrorEstimate monotone: no SectorSeries operation decreases any entry.
- I-11 No tolerance-based structural decisions in this module: tags, windows,
  sector counts, and pole classifications never depend on numeric smallness
  (the single permitted numeric comparison — pole modulus vs Radius in 2.5 —
  errors out when ambiguous instead of deciding).

(The ODE-residual spot-check invariant of RewritePlan 3.1 is owned by Solve.m
/ Transport.m — it needs the matrix; SectorSeries provides the evaluation
primitive it uses.)

---

## 5. ERROR CONTRACT

### 5.1 Loud-error catalogue

Every message is a named `Message[SectorSeries::tag, ...]` followed by a hard
abort of the operation (Throw to the module's error handler; never a returned
`$Failed` that a caller could ignore silently).  Required payload for EVERY
message: chart `Center` (and its main-variable image via `ChartMap`).
Additional payload per message:

| message | fires when | must additionally name |
|---|---|---|
| `::badtag` | tag inexact / p negative / wrong head (I-1) | the offending sector tag |
| `::dims` | shape violation (I-2/I-4/F4/F11) | sector tag, expected vs found dimensions |
| `::badcoeff` | compound/symbolic head in Coeffs (I-5) | sector tag, eps row, t column, component, the head |
| `::window` | consumer requests eps order > CompleteMax or < Min | sector tag (or "value"), requested order, window |
| `::radius` | evaluation/re-expansion point with \|t\| >= Radius | the point, Radius |
| `::originlimit` | evaluation at t = 0 of a chart with a<0 / fractional a / b≠0 / p>0 sectors | the offending sector tags; pointer to API EndpointLimit |
| `::branchmissing` | t < 0 (or D < 0) needed branch, Prescriptions give None, and a branch-sensitive sector exists | the sector tags needing a branch; the (empty/unconstraining) prescription list |
| `::branchconflict` | odd-multiplicity prescriptions derive conflicting sigmas | every factor with its derived sigma |
| `::interiorpole` | MultiplyRational pole with 0 < \|t_i\| < Radius | the exact pole, its modulus, Radius; the instruction that segmentation must split the line here |
| `::geomambiguous` | \|t_i\| vs Radius (or \|D\| vs Radius) not decidable at WorkingPrecision | both quantities and their precision |
| `::nonrational` | c not rational in (var, eps) | the non-rational subexpression |
| `::inexactdenominator` | denominator of c has inexact coefficients | the denominator |
| `::resingular` | re-expansion target at a singularity / request for tagged output | the target point |
| `::indetpade` | Pade requested with indeterminate-bearing coefficients | the indeterminate symbols found |
| `::complexpoint` | tval not real | the point |
| `::padefail` | PadeApproximant construction fails (LOUD RECORDED FALLBACK, not abort — see 2.4.2) | sector tag, eps row, component |

All messages fire even under user-level `Quiet` of generic messages (use a
dedicated message group; never wrap module internals in `Quiet[Check[...]]`).

### 5.2 The single permitted fallback

`::padefail` -> direct partial-sum evaluation, loud + recorded in
`"PadeFallbacks"` + visible to ErrorEstimate.  Nothing else in this module
falls back to anything, ever.

### 5.3 Forbidden fallbacks (enumerated temptations — each one exists in the old code; all are BANNED)

- F1 Pade failure -> silent normal evaluation with a swallowable warning
  (old: Pade.m:44-48 inside `Quiet[...]`).  Replaced by 5.2.
- F2 Out-of-range coefficient request -> return 0 (old: SafeSeriesCoefficient
  Pade.m:16-27; GetCoefficientAtPower SingularityDecomposition.m:90-99
  returns 0 for `idx < 1 || idx > Length`).  New: `::window`.
- F3 Unrationalizable/complex fitted exponent -> set b = 0 and keep towers
  explicit (old: SingularityDecomposition.m:351-362; Docs/
  FeynmanTrickBananaStatus.md:256-259).  Cannot arise — tags are never
  fitted; ANY numeric inference of a/b/p anywhere in DiffExp2 is forbidden.
- F4 Series shape mismatch -> generic `SExpand` re-expansion fallback (old:
  AddCompatibleSeries SingularityDecomposition.m:155-180).  New: `::dims`.
- F5 Iteration cap -> print-and-return-partial (old: `maxIter = 20` +
  "Warning: ... reached maximum iterations",
  SingularityDecomposition.m:311, 397-399).  No iterative extraction exists
  in the new module; returning partial results is forbidden everywhere.
- F6 Interior pole -> silently shrink Radius, or Taylor-expand a divergent
  partial fraction.  New: `::interiorpole` (math finding 8 decision).
- F7 Point outside Radius -> extrapolate the series/Pade anyway.  New:
  `::radius`.
- F8 Missing/conflicting branch -> proceed with sign +1 and a flag (old:
  AnalyticContinuation.m:55-68 sets `AnalyticContinuationFailed = True`,
  forces `uniqueSigns = {1}`, keeps going).  New: `::branchmissing` /
  `::branchconflict`.  (The old `AbortOnAnalyticContinuationFail` soft mode
  is an API/Transport policy decision; SectorSeries never proceeds.)
- F9 Indeterminates + Pade -> silently mutate global config
  (old: TurnOffPade Transport.m:178-180, 566-568).  New: `::indetpade`.
- F10 Tolerance-based structural decisions: dropping "effectively zero"
  sectors/coefficient rows (old: EffectivelyZero/EffectivelyZeroExpr
  SingularityDecomposition.m:251-299), snapping numeric exponents
  (NormalizeXPower SingularityDecomposition.m:52-64; NormalizeLogPower
  SeriesOps.m:95-106), rationalizing fitted b
  (SingularityDecomposition.m:336-365).  Only SYNTACTIC zeros may be dropped;
  the only numeric chop is the Tolerances-derived pre-Pade chop of
  COEFFICIENT VALUES (2.4.2), which never changes structure.
- F11 Zero-padding or clamping arrays to make shapes meet (the FT
  ShiftRawBoundaries zero-padding class, plan A1).  New: `::dims`/`::window`.
- F12 Evaluating t = 0 limits by substituting and Quiet-ing Power::infy (the
  old endpoint chains, e.g. FeynmanTrick/DiffExpIntegration.m:744-747's
  `Quiet[Series[...], {Power::infy, ...}]` and the Quiet[Check[..., 0]]
  endpoint fallbacks at DiffExpIntegration.m:1038-1093).  New:
  `::originlimit`; limits are an API.m/Integrate.m operation with the
  per-sector drop rule.

---

## 6. ABSORBED OLD CODE

### 6.1 DiffExp/SingularityDecomposition.m (435 lines) — role replaced entirely

What it did: reconstructed `f = Σ_i x^(a_i + b_i eps) g_i(x, eps)` from
per-eps-order towers AFTER the representation had collapsed it (D2 disease).
Functions and their dispositions:

| old | lines | disposition |
|---|---|---|
| `DecomposeSingularity` | SingularityDecomposition.m:305-402 | replaced by `SectorDecomposition` (2.8): exact data readout, no fitting |
| leading-power extraction `GetLeadingPower`/`TermXPower`/`NormalizeXPower` | 49-86, 52-64 | deleted; structure is never inferred from floats (F10) |
| `GetCoefficientAtPower` | 90-118 | deleted; out-of-range reads are `::window`, not 0 (F2) |
| `ShiftSeriesByPower` | 122-136 | becomes the exact `a`-shift inside MultiplyRational (2.5 step 3) |
| `AddCompatibleSeries` | 138-220 | deleted; shape mismatch is `::dims` (F4), aligned-array addition is trivial on dense arrays |
| `MultiplyByXMinusBEps` | 225-243 | the eps-convolution MATH is kept (correct) and reappears forward-direction in 2.4.1 and 2.6 |
| b-determination from Logx ratios + rationalization + keep-b=0 guard | 332-365 | deleted; lesson L2 below |
| `EffectivelyZero*` | 251-299 | deleted (F10) |
| `DecomposeSingularityAll` | 404-408 | per-integral mapping moves to the FT shim over `SectorDecomposition` |
| `PrintDecomposition` | 414-431 | dropped (one-line Map at the call site; see budget) |

Old-core internal consumers that die with the old representation (for the
record; their NEEDS move to Integrate.m/API.m over this module's objects):
`EvaluateEndpointLimitSectors` + `FitResidualEndpointSectors` (the Prony
fitter; DiffExp/RegularizedIntegration.m:59-61, 1954-2130),
`EvaluateLimitAtSingularity` (1921-1952), `IntegrateDecomposition[Laurent]`
(1879-1919), and the calls at RegularizedIntegration.m:1984, 2201, 2734,
2912.  Live FT call sites today: FeynmanTrick/DiffExpIntegration.m:822
(direct, verbosity-gated) and indirectly via
`EvaluateEndpointLimitSectors`/`DefiniteIntegralWithPrefactorLaurent`
(DiffExpIntegration.m:1015, 875).

Numerical lessons that MUST survive:

- L1 MULTI-SECTOR COLLAPSE.  `DecomposeSingularity` extracts ONE exponent
  with `a` = most negative power; multi-sector endpoints (banana level-2
  upper endpoint mixes `x^0`, `x^(-1+eps)`, `x^(2 eps)`,
  Docs/FeynmanTrickBananaStatus.md:78-84, 341-345) collapse to the wrong
  single tag (there: `(a,b) = (-1,1)`), and the b≠0 drop rule then discards
  or keeps whole towers (FeynmanTrickBananaStatus.md:314-329; the level-1
  master `{1,0,0,1}` came out IDENTICALLY ZERO).  Preservation: tags are
  exact and per-sector from birth; `SectorDecomposition` returns ALL sectors;
  the drop rule is applied per sector downstream.  Unit T15 pins exactly the
  banana mix.
- L2 NEVER AVERAGE EXPONENTS.  The old guard
  (SingularityDecomposition.m:351-362) kept `b = 0` when the fitted b did not
  rationalize cleanly (denominator > 16, |b| > 100, complex) because a
  collapsed average is meaningless.  In DiffExp2 the situation is
  unrepresentable; any code path that would "estimate" a tag numerically is
  forbidden (F3/F10).
- L3 THE eps-CONVOLUTION IS THE TRANSPORT BETWEEN TAG AND TOWER.
  `x^(b eps) = Exp[b eps Logx] = Σ (b Logx)^j eps^j / j!`
  (SingularityDecomposition.m:222-243) — correct mathematics, reused
  forward in evaluation (2.4.1) and re-expansion (2.6), exact.
- L4 RATIONALIZATION TOLERANCE CONTAGION.  The old module's every structural
  step ran through `FEC[RationalizationTolerance]` (lines 54, 106, 259, 284,
  336-349).  DiffExp2: zero structural tolerances in this module (I-11).

### 6.2 DiffExp/Pade.m (95 lines) — evaluation half ported, made loud

| old | lines | disposition |
|---|---|---|
| `GetPade` | Pade.m:30-53 | ported as the internal Pade builder of 2.4.2/2.9: diagonal `[m/m]` (line 35), pre-chop at `10^-ChopPrecisionVal` (line 43, constant -> Tolerances.m), per-Logx-coefficient application (lines 40-52) becomes per-(sector,row,component) |
| silent fallback | Pade.m:44-48 | made LOUD + RECORDED (`::padefail`, 5.2) — the defining change |
| `SEval1` | Pade.m:56-66 | absorbed: "build evaluable form once" = the multi-point evaluation path of 2.4 |
| `SEval2`/`SEval` | Pade.m:69-91 | absorbed into 2.4.1: branch by sign of the point (old `at >= 0` theta choice, lines 72-75/82-86) generalized to the sigma rule; `$MinPrecision = FEWorkingPrecision` blocks (lines 34, 70, 80) and `SetPrecision[..., WP]` of results (lines 76, 87) ported verbatim as lesson L5 |
| `SafeSeriesCoefficient`/`NormalizeEmbeddedSeries` | Pade.m:16-27, 60-66 | NOT ported; failure class excluded by representation (see 2.4.2 and box-doc lesson L7) |

Numerical lessons:

- L5 PRECISION FLOORS: evaluate under `$MinPrecision = WP`, SetPrecision
  results to WP (Pade.m:34, 70-77, 80-88).  The 2x-WP input raise is
  Transport.m's job (ledger), not this module's.
- L6 BRANCH AT EVALUATION: the old θp/θm projection at `at >= 0` / `at < 0`
  with `AnalyticContinuationReplacements` (Pade.m:71-75, 82-86 +
  AnalyticContinuation.m:70-79).  Ported as the sigma rule (2.4.1) with the
  principal-branch pinning stated (sign +1 = NO replacement; the old σ=−1
  rules `Logx -> Logx - 2 Pi I θm`, `x^b -> Exp[-2 Pi I b] x^b` for
  denominator > 2, and the Sqrt `(θp - θm)` case are all the single factor
  `Exp[-I Pi (a + b eps)]` relative to principal).
- L7 INERT SERIES-COEFFICIENT TAILS: upstream arithmetic used to leave
  symbolic `SeriesCoefficient[sd, {x,0,nmax}]` heads (~1e-25) in coefficients
  (box L0 root cause item 1, Docs/FeynmanTrickBoxFamilyStatus.md "RESOLVED:
  the eps^0 deficit"); Pade.m's SafeSeriesCoefficient existed to neutralize
  them.  Dense arrays + I-5 make the class unconstructible; `::badcoeff`
  guards the door.
- L8 PADE IS LOAD-BEARING (legacy review finding 11): matching/eval points at
  radius/DivisionOrder (FT pins k = 4, FeynmanTrick/DiffExpIntegration.m:272;
  classic default 3; GetCPL/GetCPR two-sided geometry Mobius.m:93-142) make
  plain partial sums the accuracy bottleneck; classic configs set
  `UsePade -> True` (Banana_example.m, FivePointNonPlanar_example.m).  The
  R6/M4 benchmark is old-with-Pade vs new-with-Pade.  Old evaluation call
  sites being replaced: Transport.m:318 (boundary-condition evaluation),
  888 (SampleAt export), 905-918 (segment evaluation + two-point error
  probe).

### 6.3 Parts of DiffExp/SeriesOps.m

- `SD` (SeriesOps.m:206-214, per-Logx-power derivative with
  `D[Log[x]^k] -> Logx` bookkeeping) -> replaced by the exact tag algebra of
  2.7.
- `LogxCoeff`/`LogxCoeffNS`/`MaxLogxPower`/`LogxPowerRange`
  (SeriesOps.m:108-136) -> unnecessary: log powers are the integer tag `p`.
- `NormalizeLogPower` (SeriesOps.m:95-106) -> deleted; lesson L10: this
  symbol's MISSING `::usage` export made the fully-qualified cross-package
  call return unevaluated, the integer filter saw no powers, the log-depth
  auto-extension never fired, and `x^-1 Logx^k` (k >= 2) sources integrated
  as constants — endpoint ODE residuals wrong from eps order 3
  (Docs/FeynmanTrickBoxFamilyStatus.md "RESOLVED: the endpoint-series
  corruption"; memory note wolfram-package-context-traps).  Preserved as:
  (a) integer `p` tags from birth — no numeric log-power recognition exists
  to silently no-op; (b) the module hygiene rule of section 2 (every export
  has ::usage) plus unit T23.
- `DecreaseSeriesOrderBy` (SeriesOps.m:196-203) -> the `"TOrderReduction"`
  option (2.4), serving the error probe.

### 6.4 What this module does NOT absorb (for boundary clarity)

The crossing operator (phase × unipotent log-chain mixing) and the matching
solve: Transport.m.  The endpoint case table, dimreg drop rule, and PV/i-delta
pairing: Integrate.m.  Tag derivation (eigenvalues, Jordan/confluence, rank
reduction): Indicial.m.  Frobenius recursion: Solve.m.  EndpointLimit user
semantics: API.m.  RegularizedIntegration.m's fitter/salvage stack: deleted
with the old core (RewritePlan I2/D2), needs replaced by exactness here.

---

## 7. DEPENDENCIES

May call (acyclic order Tolerances < Config < EpsSeries < SectorSeries):

- `Tolerances.m`: `chopFloor` (pre-Pade coefficient chop, 2.4.2), the
  geometry-ambiguity guard used by `::geomambiguous` (2.5 step 3), precision
  constants.
- `Config.m`: validated accessor for `WorkingPrecision`, `UsePade`,
  `DivisionOrder` (tail-estimate design radius, 2.6), the pinned variable and
  eps symbols, declared indeterminates (I-5).
- `EpsSeries.m`: eps-Laurent array arithmetic — expansion of rational
  c(·,eps) with Laurent-division window shifts (2.5 steps 1-2), window
  min-arithmetic for products, the canonical eps-Laurent value object (3.4).

MUST NOT call: Indicial.m, Solve.m, Transport.m, Integrate.m, API.m, anything
under ``DiffExp` `` (old core) or ``FeynmanTrick` ``.  Mathematica built-ins
used: `PadeApproximant`, exact `Roots`/`Factor`/`Root` for denominators,
standard arithmetic.  No global mutable state except the message definitions;
no writes to Config.

---

## 8. UNIT TESTS

File: Tests/test_sectorseries.m (M2 gate, RewritePlan §6).  All expected
values exact unless noted; `wp = 100` unless noted.  Every loud-error test
asserts BOTH that the named message fired AND that the payload contains the
required names (chart/sector/order).

1. `t01_validate_roundtrip` — a hand-built two-sector LocalSolution passes
   `ValidateLocalSolution` unchanged (SameQ).
2. `t02_validate_inexact_tag` — sector with `a = 0.5` (machine float) ->
   `::badtag` naming the sector.
3. `t03_validate_shape` — sectors with mismatched Ncomp -> `::dims` with
   expected/found dimensions.
4. `t04_validate_badcoeff` — a `SeriesCoefficient[...]` head planted in
   Coeffs -> `::badcoeff` naming (sector, row, column, component).
5. `t05_eval_regular_exp` — single `(0,0,0)` sector, `Ncomp = 1`, eps row 0
   only, `c[0, n] = 1/n!` for `n = 0..20`; evaluate at `t = 1/10`, UsePade
   False.  Assert `|value - Exp[1/10]| < 10^-20` (tail `(1/10)^21/21!`) and
   `Precision >= wp` (lesson L5).
6. `t06_eval_branch_halfinteger` — sector `(1/2, 0, 0)`, `g = 1`; at
   `t = -1/4`: sigma +1 -> exactly `I/2`; sigma −1 -> `-I/2` (via the
   `"ImSign"` override).  Pins the principal-branch convention (L6).
7. `t07_eval_branch_logchain` — sector `(0, 0, 1)` with row support `{0}`,
   `c[0,0] = 1` (value `eps Log t`); at `t = -1/4`, sigma +1: eps-order-1
   coefficient exactly `-Log[4] + I Pi`.
8. `t08_eval_beps_tower` — sector `(0, 2, 0)`, `g = 1`, CompleteMax = 5; at
   `t = 1/3`: value eps-order k coefficient exactly `(2 Log[1/3])^k / k!`
   for `k = 0..5`; result CompleteMax = 5 (window honesty for the entire
   exp-factor).
9. `t09_eval_resonant_kmin` — sector `(0, 0, 2)` with rows from `kmin = -2`,
   `c[-2, 0] = 1` (true log-chain: value `Log[t]^2/2`); at `t = 1/2`: value
   eps-order 0 exactly `Log[1/2]^2/2`; `MinPower` of the result = 0; no
   window error (the RewritePlan 3.1 exemption).
10. `t10_eval_radius_error` — point at `t = Radius` -> `::radius` naming the
    point and Radius.
11. `t11_eval_origin_error` — singular chart (any b≠0 sector) evaluated at
    `t = 0` -> `::originlimit` listing the offending tags.
12. `t12_eval_branch_missing` — `Prescriptions -> {}` chart with a
    `(1/2,0,0)` sector at `t = -1/10` -> `::branchmissing`; same chart with
    ONLY a `(0,0,0)` sector at `t = -1/10` evaluates fine (branch
    irrelevant).
13. `t13_chart_imsign` — Prescriptions list with one odd-multiplicity entry
    (`Sign = -1`, `LeadingCoeffSign = -1`) -> sigma = +1; adding an even-
    multiplicity entry with opposite sign changes nothing; adding a second
    odd entry deriving −1 -> `::branchconflict` naming both factors.
14. `t14_pade_exact_rational` — coefficients `(-1)^n`, `n = 0..10` (truncated
    `1/(1+t)`); evaluate at `t = 3/4`: plain sum error `>= 10^-3`
    (`(3/4)^11/(7/4) ≈ 0.024`), Pade `[5/5]` reproduces `4/7` to better than
    `10^-(wp-10)` (diagonal Pade of a rational function is exact).  This is
    the "Pade is load-bearing" pin (L8).
15. `t15_pade_loud_fallback` — via the test hook (internal option forcing
    PadeApproximant failure): `::padefail` fires naming (sector, row,
    component); returned value equals the plain-sum value exactly;
    `"PadeFallbacks"` records exactly one entry.  Assert the message is NOT
    suppressed by `Quiet[..., General::xxx]`-style generic quieting.
16. `t16_pade_indeterminate` — coefficients containing a declared
    indeterminate symbol with `"UsePade" -> True` -> `::indetpade`; with
    `"UsePade" -> False` evaluation succeeds and returns linear forms in the
    indeterminate.
17. `t17_multipoint_consistency` — `EvaluateLocalSolution[ls, {t1, t2}]`
    equals `{EvaluateLocalSolution[ls, t1], EvaluateLocalSolution[ls, t2]}`
    entrywise (with UsePade True).
18. `t18_torder_reduction` — `"TOrderReduction" -> 3` on t05's object equals
    evaluation of the manually 3-column-truncated object (the error-probe
    seam).
19. `t19_mul_center_pole` — `f ≡ 1` (regular one-sector object, Radius 1)
    times `c = 1/(t (1 - t))`: result is the SINGLE canonical sector
    `(-1, 0, 0)` with `c[0, n] = 1` for all kept n (partial fractions
    `1/t + 1/(1-t)`, the second pole at `|1| >= Radius` folded as the
    geometric series, then integer-spaced-merge).  Exact coefficient check.
20. `t20_mul_eps_laurent_shift` — `c = 1/(2 eps)`: result windows shift down
    by exactly 1 (`Min` and `CompleteMax`), coefficients halved, NO error
    (the shift-not-inversion contract).  `c = 1/(1 + eps)`: windows
    unchanged, rows are the alternating geometric convolution — check rows
    0..3 exactly.
21. `t21_mul_interior_pole` — Radius 1, `c = 1/(t - 1/2)` ->
    `::interiorpole` naming pole `1/2` and Radius 1.
22. `t22_mul_nonrational` — `c = Sqrt[t]` -> `::nonrational`.
23. `t23_mul_main_coordinates` — affine ChartMap `x = 3 + t`, `c = 1/(x - 3)`
    with `"Coordinates" -> "Main"`: equals `t19`-style center-pole shift
    (`a -> a - 1`).  Exact.
24. `t24_reexpand_singular_source` — sector `(-1, 0, 0)`, `g = 1` (i.e.
    `1/t`), Radius 1; re-expand at `D = 1/4`, N = 10: expected single
    `(0,0,0)` sector with `c[0, m] = (-1)^m 4^(m+1)` (from
    `1/(D + u) = Σ (-1)^m D^(-m-1) u^m`), new `Radius = Min[1/4, 3/4] = 1/4`,
    `Prescriptions = {}`; recorded TailEstimate matches the formula in 2.6
    with the test's DivisionOrder.  Exact coefficient check for `m = 0..10`.
25. `t25_reexpand_beps_collapse` — `f = t^eps` (sector `(0,1,0)`, `g = 1`,
    CompleteMax 3); re-expand at `D = 1/2`, N = 4: single `(0,0,0)` sector;
    eps row k, column 0 = `Log[1/2]^k / k!`; eps row 1 columns 1..4 =
    the series of `Log[1 + 2 u]`: `{2, -2, 8/3, -4}`.  Asserts the documented
    b-collapse (and that the output carries NO `b != 0` sector).
26. `t26_reexpand_errors` — `D = 0` on a singular chart -> `::resingular`;
    `|D| >= Radius` -> `::radius`.
27. `t27_diff_tagged` — sector `(a, b, 1) = (1/2, 3, 1)` with row support
    `{0}`, `c[0, n] = δ_{n,0}` (value `t^(1/2 + 3 eps) eps Log t`):
    derivative = sectors `(-1/2, 3, 1)` with rows
    `{k=0: 1/2, k=1: 3}` at `n = 0`, and `(-1/2, 3, 0)` with row
    `{k=1: 1}` — matches
    `d/dt [t^(1/2+3eps) eps Log t] = (1/2 + 3 eps) t^(-1/2+3eps) eps Log t
    + t^(-1/2+3eps) eps` exactly.  Also: differentiating t05's exp object
    reproduces the factorial shift.
28. `t28_diff_resonant_window` — differentiate t09's `(0,0,2)`/`kmin = -2`
    object: the produced `(−1,0,1)` sector has rows from `-1` (= `-(p-1)`),
    no window violation.
29. `t29_diff_main_affine` — ChartMap `x = 3 + 2 t`: main-variable derivative
    = chart derivative times exactly `1/2`.
30. `t30_sector_decomposition_banana_mix` — build the banana level-2 endpoint
    mix exactly (FeynmanTrickBananaStatus.md:344): sectors
    `(0,0,0)`, `(-1,1,0)`, `(0,2,0)` with distinct marker coefficients.
    `SectorDecomposition` returns exactly these three tags, sorted, with the
    marker arrays untouched and windows attached.  THE FT-shim replacement
    pin: the old code collapsed this very input to the single fitted tag
    `(-1, 1)` (FeynmanTrickBananaStatus.md:341-349).
31. `t31_canonicalize_merge` — sectors `(-1, 1, 0)` and `(0, 1, 0)` with
    known arrays merge to one `(-1, 1, 0)` sector with the second's rows
    shifted one column right and added; duplicate `(0, 2, 0)` sectors add;
    an all-syntactic-zero sector is dropped; a numerically tiny
    (`10^-200`) but nonzero sector is NOT dropped (F10 pin).  The companion
    `t31_canonicalize_preserves_overflowing_integer_shift` pin uses a shift
    wider than the stored slab and asserts that its nonzero coefficient
    remains at the higher-`a` tag instead of being truncated away.
32. `t32_window_request_error` — a consumer helper requesting value order
    `CompleteMax + 1` from t08's evaluation -> `::window` naming the order
    and window.
33. `t33_export_hygiene` — every symbol referenced as
    ``DiffExp2`SectorSeries`...`` from the module's own test file resolves to
    a defined symbol with a `::usage` string (the L10/NormalizeLogPower trap
    pin; implement by checking `Names["DiffExp2`SectorSeries`*"]` against the
    spec's public list and `::usage =!= ""`).

---

## 9. LINE BUDGET

RewritePlan 3.2: SectorSeries.m ~400 lines, plus the Pade port decision
~+100 -> ~500 total inside the 3.3k core target (3.5k hard ceiling).

Indicative allocation: validation/canonicalization 60; ChartImSign 25;
evaluation incl. branch rule and windows 90; Pade port 100; MultiplyRational
(partial fractions, pole classification, convolution) 110; Reexpand incl.
truncation contract 70; Differentiate 35; SectorDecomposition 15; messages 25.

If over budget, cut IN THIS ORDER (never cut loud errors, window honesty, or
Pade loudness):
1. `"Coordinates" -> "Main"` option on MultiplyRational/Differentiate (~-30):
   the FT shim converts coordinates instead; t23/t29 retarget the shim.
2. Mobius chain rule in DifferentiateLocalSolution (~-15): loud
   not-implemented error naming the chart (Integrate.m rejects Mobius anyway;
   classic Mobius lines never differentiate through this path).
3. Multi-point evaluation form (~-20): callers loop; ACCEPT the cost of
   rebuilding Pade per point and SAY SO in the R6 benchmark notes.
4. `CanonicalizeLocalSolution` integer-spaced-`a` merging (~-15): keep
   duplicate-tag merging only; Solve.m already emits family-merged sectors;
   I-6 weakens accordingly (documented).
5. `PadeEvaluate` as a public symbol (~-5): fold into Private with t14/t15
   going through EvaluateLocalSolution.

What may NEVER be cut: the truncation contract of 2.6 (math finding 12), the
interior-pole error (math finding 8), the vector component index (math
finding 10), the recorded-loud Pade fallback, the branch-conflict error, the
resonant `kmin = -p` window exemption.

---

## 10. OPEN QUESTIONS

- Q1 EpsSeries object naming: this spec assumes
  `<|"MinPower","CompleteMax","Coefficients"|>` for the eps-Laurent array
  (3.4).  Docs/specs/EpsSeries.md is normative for those key names; reconcile
  at M1 before implementation.
- Q2 Should `::padefail` escalate to a hard abort under a strict Config flag
  (e.g. `"StrictPade" -> True`) for parity-gate runs, so a fallback can never
  hide inside a green benchmark?  (Default behavior is fixed by 5.2 either
  way.)
- Q3 Evaluation guard band: is `|t| < Radius` strict-inequality enough, or
  should Tolerances.m define a fractional guard (e.g. reject
  `|t| > Radius (1 - 10^-3)`) given Transport places points at Radius/k?
  Decision belongs to Tolerances.m/Transport.m; SectorSeries takes the
  threshold as given.
- Q4 `SectorDecomposition` ordering: lexicographic `(a, b, p)` is proposed
  here; Integrate.m and the FT shim doc must adopt the same order — confirm
  in their specs.
- Q5 Complex evaluation points: v1 rejects them (`::complexpoint`); the
  complex-singularity ledger item (ghost projection) may eventually want
  evaluation slightly off-axis.  Defer; revisit with the pentagon line work.
- Q6 Re-expansion consumers: in the v1 pipeline only SaveExpansions/
  ToPiecewise-equivalents and off-grid evaluation need Reexpand (Transport
  matches by pointwise evaluation).  If M4 shows it unused on the parity
  lines, consider demoting to the cut list — but the truncation contract
  stays specified, since API.m's ToPiecewise needs it.
- Q7 Indeterminate-bearing coefficients and the tail estimate (2.6): the
  Max-over-comp of an indeterminate linear form is ill-defined; proposal:
  per-indeterminate tail estimates mirroring the old
  ComputeErrorsPerIndeterminate (Transport.m:947-962).  Decide with API.m's
  symbolic-BC spec.
