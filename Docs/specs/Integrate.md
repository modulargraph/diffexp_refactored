# Module spec: DiffExp2/Integrate.m  (~500 lines)

Status: M0 deliverable per Docs/RewritePlan.md section 3.2 and 6 (task 4-11).
Contract document — an implementation agent must be able to build this module
from THIS file + Docs/RewritePlan.md + the cited old code alone.
Old-code citations refer to the frozen oracle tree (DiffExp/, FeynmanTrick/)
at the plan baseline.

---

## 1. PURPOSE

Integrate.m computes definite integrals of LocalSolutions over chart windows
on a real line — the FT pipeline's `∫_0^1 dx` of a combined master series
across a chain of charts — using exact per-sector antiderivatives in the
sector-native representation: every integrand monomial
`t^(a+n+b·eps) (eps·Log t)^p / p!` has a closed-form antiderivative whose
eps-Laurent structure (denominators `(a+n+1+b·eps)^(j+1)`, j = 0..p) is
computed exactly by EpsSeries window arithmetic, so the endpoint
regularization that old `DiffExp/RegularizedIntegration.m` (3015 lines)
reached through iterated IBP regularization steps, endpoint-sector Prony
fitting, salvage paths, and averaged-exponent fallbacks becomes a twelve-cell
case table with one load-bearing boundary rule (the dimreg convention
`t^(a+n+1+b·eps)|_{t=0} := 0` for `b ≠ 0`) and one loud-error class (`b = 0`
divergences that fail the object-level cancellation check). The module also
handles interior singular points (PV / finite-part for `b = 0`, exact
phase-paired iδ crossing for `b ≠ 0`, both with the two-half-segment pairing
enforced by assert), owns the endpoint-limit drop rule and divergence gate
(`EndpointSectorLimit`, 2.4 — it shares the 2.2.3 gate), rejects non-affine
charts loudly, and assembles the segment chain into one final eps-Laurent per
component with honest EpsWindow and TWindow metadata — no padding, no silent
drops, no fitting.

---

## 2. PUBLIC SYMBOLS

Context: `DiffExp2`Integrate``. Exactly four exports. Everything else is
`Private``.

### 2.1 `SectorMonomialIntegral[m, b, p, T, kMaxOut]`

The case-table primitive. Computes, exactly,

    ∫_0^T  t^(m + b·eps) · (eps·Log t)^p / p!  dt          (T > 0)

with the lower boundary handled by the rules of section 2.1.2.

Arguments:
- `m`     — exact rational (or exact algebraic) local power `a + n`.
            NEVER a machine/arbitrary-precision float (assert; see E10).
- `b`     — exact rational/algebraic eps-slope of the sector tag. The test
            `b == 0` is EXACT structural comparison on the tag, never a
            numerical zero test (disease D1: the old IntegerQ-on-floats and
            `NumericZeroQ[effectiveB]` classes,
            DiffExp/RegularizedIntegration.m:497-498, 571, die here).
- `p`     — integer ≥ 0 (assert).
- `T`     — positive number at working precision (raised to 2×WP on entry,
            ledger lesson), or exact rational. `T > 0` is asserted; negative
            endpoints never reach this primitive (the branch-resolved wrappers
            in 2.2 handle the negative arm).
- `kMaxOut` — integer: highest eps order requested in the result. Needed
            because the closed form is exact to all orders; the result is
            truncated, with `CompleteMax = kMaxOut` honestly.

Returns: an EpsSeries value (section 3.4; DEC-13 shape) with
`EpsWindow["Min"] = p − (p+1) = −1` in the pole cell, `= p` otherwise, and
`EpsWindow["CompleteMax"] = kMaxOut`.

#### 2.1.1 The antiderivative (exact closed forms)

Let `α := m + 1 + b·eps`. For `α ≢ 0` (i.e. NOT(`b = 0` and `m + 1 = 0`)):

    F(t) = t^α · Σ_{j=0}^{p} (−1)^j · eps^p · Log^(p−j) t / ((p−j)! · α^(j+1))

(Verification: d/dt F = t^(m+b·eps) (eps Log t)^p / p!; the j-sum telescopes.
This is old `IntegratePowerLogMonomial`,
DiffExp/RegularizedIntegration.m:1703-1724, and old `monomialBasis`'s
`basisExpr`, DiffExp/RegularizedIntegration.m:1428-1436, unified — same
formula, `Binomial[p,j]·j! = p!/(p−j)!`.)

For `α ≡ 0` (`b = 0` AND `m + 1 = 0`):

    F(t) = eps^p · Log^(p+1) t / ((p+1) · p!)

(old DiffExp/RegularizedIntegration.m:1718-1719.)

`1/α^(j+1)` is expanded in eps by EpsSeries: when `m+1 ≠ 0` it is the
eps-regular geometric series in `(−b/(m+1))·eps`; when `m+1 = 0`, `b ≠ 0`
it is the exact Laurent shift `1/(b·eps)^(j+1)` (EpsSeries Laurent-division
contract, math review finding 6(iv)). `T^(b·eps) = Σ_r (b·Log T)^r eps^r/r!`
is assembled to the needed length with `Log T` evaluated once per chart at
≥ 2×WP under `$MinPrecision = WP` (ledger: precision floors). No symbolic
`Series[]` calls; everything is finite EpsSeries window arithmetic.

#### 2.1.2 Boundary evaluation at t = 0 — THE RULES

- **b ≠ 0 (any m, any p): `F(0) := 0`.** This is the dimreg convention
  `t^(m+1+b·eps) Log^q t |_{t=0} := 0`, by analytic continuation in eps
  (continue from the region Re(m+1+b·eps) > 0 where the limit is genuinely
  zero). It is load-bearing: it IS the integral counterpart of the FT drop
  rule "sectors with b ≠ 0 are set to zero at the endpoint, even when a < 0"
  (old usage text DiffExp/RegularizedIntegration.m:57-59; paper rule quoted
  at FeynmanTrick/DiffExpIntegration.m:899-907). State it in code as a named
  function `dimregLowerBoundary[]` returning exact 0 with a comment citing
  this spec section. Every endpoint regularization the old code performed by
  iterated IBP steps (`ApplyRegularizationStep`/`RegularizeIntegrand`,
  DiffExp/RegularizedIntegration.m:363-466) or by the subtraction machinery
  (`IntegrateAnalyticRegularizedBySubtractionLaurent`,
  DiffExp/RegularizedIntegration.m:1333-1692) is THIS rule applied to the
  closed form — the new module needs no iteration and no subtraction.
- **b = 0, m + 1 > 0: `F(0) = 0`** (true limit, `t^(m+1) Log^q t → 0`).
- **b = 0, m + 1 ≤ 0 (any p): divergent.** LOUD ERROR E2 — see the
  cancellation contract in 2.2.3. Never reached inside
  `SectorMonomialIntegral` itself: the wrapper must have removed (cancelled)
  or errored on such monomials BEFORE dispatching here; if one arrives
  anyway, raise E2 (defense in depth).

#### 2.1.3 THE FULL CASE TABLE

`{b = 0 / b ≠ 0} × {m+1 < 0 / = 0 / > 0} × {p = 0 / p > 0}`, twelve cells.
"shift" = eps-window shift of the result relative to the regular-cell
baseline `Min = p` (the explicit `eps^p` of the tag normalization).

| # | b | m+1 | p | result over [0,T] | eps structure / window |
|---|---|-----|---|-------------------|------------------------|
| 1 | =0 | >0 | =0 | `T^(m+1)/(m+1)` | exact, eps-independent; contributes at coefficient order k; no shift |
| 2 | =0 | >0 | >0 | `T^(m+1) · Σ_{j=0}^p (−1)^j eps^p Log^(p−j)T / ((p−j)!(m+1)^(j+1))` | finite eps-polynomial of degree p; orders k+p; no shift |
| 3 | =0 | =0 | =0 | DIVERGENT (`Log t → −∞`) | LOUD ERROR E2 at a true endpoint, unless the merged coefficient is zero (2.2.3); PV pairing if interior (2.2.4) |
| 4 | =0 | =0 | >0 | DIVERGENT (`Log^(p+1) t`) | LOUD ERROR E2 — note: ANY p, not only p = 0 (math review finding 6(ii)); interior: PV pairing gives `eps^p (Log^(p+1) t2 − Log^(p+1)\|t1\|)/((p+1)·p!)` |
| 5 | =0 | <0 | =0 | DIVERGENT (power) | LOUD ERROR E2 unless cancelled; interior: symmetric-cutoff finite part (2.2.4) |
| 6 | =0 | <0 | >0 | DIVERGENT (power × log) | LOUD ERROR E2 unless cancelled; interior: finite part (2.2.4) |
| 7 | ≠0 | >0 | =0 | `T^α/α` | `1/α` eps-regular (geometric); orders ≥ k; no shift |
| 8 | ≠0 | >0 | >0 | full j-sum of 2.1.1 at t = T | all `α^(j+1)` eps-regular; orders ≥ k+p; no shift |
| 9 | ≠0 | =0 | =0 | `T^(b·eps)/(b·eps)` | exact Laurent SHIFT by 1: orders ≥ k−1; `Min` drops by 1 AND `CompleteMax` drops by 1 (output order K needs coefficient order K+1) |
| 10 | ≠0 | =0 | >0 | `Σ_{j=0}^p (−1)^j eps^p Log^(p−j)T · T^(b·eps) / ((p−j)!(b·eps)^(j+1))` | denominators of depth 1..p+1: an eps-POLE OF DEPTH p+1 below the tag's natural `eps^p` scaling (j = p term: `(−1)^p T^(b·eps)/(b^(p+1) eps)` after the `eps^p` cancellation) — NOT a single `1/(b·eps)` (math review finding 6(i)). Orders ≥ k−1; window shift p+1 relative to the regular cells' k+p |
| 11 | ≠0 | <0 | =0 | `T^α/α`, `F(0) = 0` by the dimreg rule | `1/α` eps-regular since `m+1 ≠ 0`; NO error, NO shift — this cell is the analytic regularization of non-integrable powers; the old IBP-step machinery (RegularizedIntegration.m:363-466) collapses to it |
| 12 | ≠0 | <0 | >0 | full j-sum, `F(0) = 0` | eps-regular; no shift |

Cell selection is EXACT on the tags: `m + 1 = 0` requires the rational
`a + n + 1` to be exactly zero (only possible for integer `a ≤ −1`); no
snapping, no tolerance (kills the `IntegerPower`/`NumericZeroQ` snap-and-drop
class, DiffExp/RegularizedIntegration.m:273-284, 705-741).

Pole-depth bookkeeping, stated precisely (reconciles finding 6(i) with the
RewritePlan 3.4 budget term `enhancement = max(p+1)`):

- A single coefficient order k of a pole-cell sector (`b ≠ 0`, `a+n+1 = 0`
  reachable) contributes to ALL output orders ≥ k−1; completeness of output
  order K therefore needs coefficient orders up to K+1. With the EpsWindow
  `[κmin, κmax]` of the Coeffs array, the sector's integral has window
  `[κmin−1, κmax−1]` versus `[κmin+p, κmax+p]` for non-pole sectors — a
  relative loss of p+1 orders.
- Homogeneous true-resonance log-chains have `κmin = −p` by construction
  (RewritePlan 3.1 EXEMPTION), so the pole cell reaches eps order
  `−p−1`: the integral of an O(eps^0) log-chain sector has a genuine
  eps-pole of depth p+1. This is the precise content of "pole depth p+1".
- The static budget's `enhancement(Lk) = max over pole-hit endpoint sectors
  of (p+1)` (RewritePlan 3.4) is exactly the window-shift difference above;
  Integrate.m does NOT consume the budget — EpsWindow propagation is the
  enforcement (any miscount upstream surfaces as error E4 here).

### 2.2 `IntegrateLocalSolution[ls, {t1, t2}]`

Definite integral of a LocalSolution over a chart-coordinate interval.

Arguments:
- `ls` — LocalSolution (section 3.1), affine chart (else E1).
- `{t1, t2}` — chart-coordinate bounds, `t1 < t2` after normalization
  (callers pass main-line-ordered bounds; orientation is handled in 2.3).
  Both must satisfy `|t| < ls["Radius"]` (assert E9a). Endpoint-equality
  with the chart center is decided exactly or within `snapTol`
  (Tolerances.m); a bound that is "almost" 0 but not snapped is E5.

Returns: `<| "Values" -> {EpsSeries value per component (3.4)},
"EpsWindow" -> <|"Min"->.., "CompleteMax"->..|> (of the Values),
"ErrorEstimate" -> <|k -> est..|> (t-truncation tail per eps order),
"TWindowUsed" -> nmax |>` — the metadata keys are SIBLINGS of the
series values, never fields inside them (3.4).

#### 2.2.1 Endpoint classification (exact, three cases)

1. **Endpoint at the chart center** (`t1 = 0` or `t2 = 0`): per-sector
   endpoint integrals via `SectorMonomialIntegral` (boundary rules 2.1.2),
   after the divergence/cancellation gate 2.2.3. `t2 = 0` with `t1 < 0` is
   the mirrored case: substitute `t → −t` is NOT performed; instead the
   negative arm uses the branch-resolved evaluation 2.2.5 with the dimreg
   rule at 0 unchanged (it is branch-independent: the limit is 0 from either
   side under continuation).
2. **Center strictly inside** (`t1 < 0 < t2`): the interior-crossing path
   2.2.4. BOTH half-windows are processed inside this single call — the
   two-half-segment pairing is enforced here by construction and by assert
   I4; callers cannot obtain one half of a crossing (passing `{t1, 0}` alone
   is legal ONLY when the caller's main-line integration bound genuinely
   ends at the chart center, which the segment-assembly layer 2.3 verifies
   against `{xLo, xHi}`; see E6).
3. **Center outside** `[t1, t2]` (both bounds same sign, neither zero):
   plain antiderivative difference `F(t2) − F(t1)` per monomial, both
   boundary terms evaluated per 2.2.5. No drop rule, no divergence possible
   (all `Log` and power values finite). Window shift only from pole-cell
   denominators if `b ≠ 0`, `m+1 = 0` monomials exist — note these are NOT
   singular here, but `1/(b·eps)` still appears in `F` at each endpoint
   separately; implement this case with the COMBINED closed form
   `(t2^β − t1^β)/β = Σ_{r≥1} β^(r−1)(Log^r t2 − Log^r t1)/r!` (β = b·eps,
   branch-resolved logs) so the spurious per-endpoint shift never enters
   the window arithmetic. Same for the `Log^(p+1)` α ≡ 0 form (difference
   of logs, finite).

#### 2.2.2 Per-sector assembly

For each sector `s = <|a, b, p, Coeffs c[k,n,comp]|>` and component:

    contribution = Σ_{k=κmin}^{κmax} Σ_{n=0}^{nmax} c[k,n,comp] · eps^k ·
                   SectorMonomialIntegral[a+n, b, p, T, K_s − k]

where `kMaxOut = K_s − k`, with `K_s = κmax + p` for sectors that cannot hit
the pole cell and `K_s = κmax − 1` for pole-hit sectors (`b ≠ 0` with
`a + n + 1 = 0` reachable), so that after the `eps^k` shift every
contribution shares the sector's uniform CompleteMax `K_s` of 2.1.3
(REVIEW-math D22). Implemented sector-at-a-time (factor `T^(b·eps)` and the
`Log^q T` powers once per sector; per-n only the geometric `1/α^(j+1)`
series differ), with `nmax = ls["TWindow"]["CompleteMax"]` — t-orders beyond
CompleteMax do not exist in the object (TWindow is the truncation-order
record ONLY, DEC-9: recursion matrices are exact polynomials, so the old
coupling-depth degradation has no premise in the new core). The t-truncation
tail per eps order is estimated as
`|c[k,nmax,comp]| · T^(a+nmax+1) / |a+nmax+1|` when `a+nmax+1 ≠ 0`, and
`|c[k,nmax,comp]| · |Log T|` when `a+nmax+1 == 0` (the next term's
antiderivative magnitude in the log cell; REVIEW-minimalism 28), accumulated into
`ErrorEstimate[k']` for every output order k' the term feeds, and ADDED to
the LocalSolution's incoming `ErrorEstimate`. Estimate > 1 (relative to the
result's leading coefficient at that order) ⇒ abort, error E8 (the ported
abort>1 rule, RewritePlan 3.1 "ErrorEstimate").

Output EpsWindow: computed BY EpsSeries operations (add/scale/Laurent-divide
carry windows); additionally the static per-sector account
(`[κmin+p, κmax+p]` regular, `[κmin−1, κmax−1]` pole-hit) is computed
independently and asserted equal (invariant I5).

#### 2.2.3 The b = 0 divergence gate and the object-level cancellation check

Precondition: the LocalSolution handed to Integrate is the ASSEMBLED
combination — in FT use, `Σ_j c_j(x,eps)·f_j(x,eps)` has already been formed
by SectorSeries rational-multiply + add (RewritePlan 3.3 "integrate:
rational-multiply + exact cancellation at object level"), and same-tag
sectors are merged (one sector per exact `(a,b,p)`; invariant I2), and
integer-spaced same-(b,p) sectors merged to minimal a (SectorSeries.md 2.2;
assert at entry: no two sectors share b, p with integer a-difference —
violation is E7-class; REVIEW-math D12). The gate's per-monomial
coefficients are well-defined only on this canonical form. The old
warning "Combines the integrand before integration to handle cancellation of
poles in IBP coefficients" (FeynmanTrick/DiffExpIntegration.m:16-21,
532-545) is here an input invariant, not a hope.

The gate, run BEFORE any monomial dispatch when an integration bound sits at
the chart center:

For every `b = 0` sector and every divergent local power
(`m = a + n ≤ −1`, any p), and every eps order k in the window:

1. The relevant coefficient is `c[k, n, comp]` of the MERGED sector — i.e.
   the sum of the offending coefficients across the master combination has
   already been taken (merging same-tag sectors IS the sum; do not error on
   pre-merge data). Distinct p values are distinct monomials and may NOT be
   summed against each other.
2. If the coefficient is exactly zero: drop, no event.
3. If it is numerically zero — `ESCoeffZeroQ[c, scale_k]` at
   `Tol["LaurentLeadTol"]` (= `Max[10^(-Floor[chopDigits/2]), 10^-24]` per
   Tolerances.md as amended by REVIEW-math D1; the campaign class
   Recurrence.m:613-618; DEC-2's hard 10^-24 floor and 4-decade band), with
   `scale_k = max_n' |c[k,n',comp]|` over the merged combination at that eps
   order (DEC-4: RELATIVE, per-eps-order; the divergence error fires only if
   the offending coefficient exceeds laurentLeadTol relative to that scale)
   — drop it and count the event (result metadata
   `"CancelledDivergent" -> count`; an info-level `Config`PrintInfo` at
   Verbosity≥2, never a warning).
4. Otherwise: LOUD ERROR E2, carrying every field listed in section 5.
   There is NO step 5. In particular the old "resonant endpoint coefficient
   with zero epsilon regulator; dropping contribution" silent drop
   (DiffExp/RegularizedIntegration.m:1470-1477) — the proximate killer of
   the pentagon leading pole (RewritePlan sections 2, 10) — is forbidden
   (F1), as is substituting a formal unit regulator `b := 1`
   (DiffExp/RegularizedIntegration.m:654-660; forbidden F2 — the old code's
   own comment documents why: regulating regular powers "injects artificial
   higher-epsilon terms", RegularizedIntegration.m:562-570).

No deferral case exists: coefficients above `EpsWindow["CompleteMax"]` do
not exist in the object, so "divergence only affects orders beyond the
requested window" (the old warning split,
DiffExp/RegularizedIntegration.m:1668-1689) cannot arise. Do not implement a
deferral branch.

#### 2.2.4 Interior crossing (t1 < 0 < t2): PV / iδ with enforced pairing

All monomials are integrated as a SINGLE paired closed form over both arms —
never as two independent half-integrals later added (that is what leaked
complex logs and spurious eps-shifts in the old code; commit 9aeb300,
Docs/FeynmanTrickBoxFamilyStatus.md "Fixed this campaign").

- **b = 0, m+1 > 0 (and all p):** regular; plain `F(t2) − F(t1)` with REAL
  values (integer m+1; real logs per the next bullet for p > 0).
- **b = 0, divergent monomials (m+1 ≤ 0) surviving the gate 2.2.3 chop:**
  symmetric-cutoff finite part with the REAL-log convention on the negative
  arm (`Log t → Log|t|`; old `$InteriorSplitRealLog`,
  DiffExp/RegularizedIntegration.m:523-529, 1833-1843): compute
  `FP = [F_real(t)]_{t1}^{t2}` where `F_real` uses `Log|t|` and real powers,
  and the formal boundary terms at `0^±` are dropped ONCE, jointly (they are
  computed as a single symmetric-cutoff limit, not per-arm). Closed forms:
  - `m = −1, p = 0`: `c · Log(t2/|t1|)` — the principal value, REAL
    (campaign pin: Tests/test_interior_singular_integration.m:41-51,
    `Log(B/|A|) + (B−A)` for `g = 1+x` over `{−1/5, 1/4}` — the box-family
    interior IBP pole at xx1 = 1/4 class).
  - `m = −1, p > 0`: `eps^p (Log^(p+1) t2 − Log^(p+1)|t1|)/((p+1)·p!)`
    (the cutoff terms cancel identically between arms).
  - `m < −1`: `[t^(m+1)/(m+1)]` evaluated at the outer endpoints with real
    arithmetic (`t1^(m+1)` real for integer m), cutoff terms dropped — the
    Hadamard finite part (old fpAt, DiffExp/RegularizedIntegration.m:606-650;
    pin: `(−1/B + 1/A) + (B − A)` for `x^(−2)(1+x^2)`,
    Tests/test_interior_singular_integration.m:53-64). Note: for even m the
    cutoff terms do NOT cancel pointwise between arms; dropping them is the
    finite-part PRESCRIPTION, stated here as the defined convention (it is
    deterministic and exact, not a fallback). Count such drops in result
    metadata `"FinitePartDrops" -> count`.
  - Reality assert: for real input coefficients these paths must produce
    real results (Im ≤ tol; invariant I7). This is the 9aeb300 regression
    pin.
- **b ≠ 0 monomials:** exact iδ crossing using the chart prescription sign
  σ ∈ {+1, −1}: `σ := SectorSeries`ChartImSign[ls]`; E3 fires iff it returns
  None while a phase is needed, or throws conflict (the sign derivation is
  NOT reimplemented here — REVIEW-minimalism 18). With
  the principal-branch-far-side convention (negative arm:
  `t^α := e^(iπσα)|t|^α`, `Log t := Log|t| + iπσ`; pinned to SectorSeries'
  sigma rule (its 2.4.1), the single owner of the scalar branch rule —
  REVIEW-math D28; agreement with Transport.m's sector-level mixing operator
  is guaranteed by the cross-module parity tests (Transport T8, test 14
  here), never by a sibling call), the paired
  integral of the pole cell is computed as ONE regular series:

      ∫_{t1}^{t2} t^(−1+b·eps) dt
        = (t2^β − e^(iπσβ)|t1|^β)/β                          (β = b·eps)
        = Σ_{r≥1} β^(r−1) (Log^r t2 − (Log|t1| + iπσ)^r) / r!
        = Log(t2/|t1|) − iπσ + O(eps).

  The `1/(b·eps)` of each arm cancels in the pair; implement the combined
  series directly (never difference-of-two-shifted-Laurents), so the
  crossing contributes NO eps-window shift — assert I4 verifies the
  assembled window equals the no-pole static account (math review
  finding 7(v): `(e^(iπb·eps) − 1)/(b·eps)` is regular). General p > 0:
  same pairing applied to every `1/(b·eps)^(j+1)` term; the combined series
  is `Σ_j (−1)^j eps^p/(p−j)! · d^j/dβ^j`-style regular forms — implement
  via the generating identity
  `∫_{t1}^{t2} t^(−1+β)(eps Log t)^p/p! dt = (eps^p/p!) ∂_β^p[(t2^β − e^(iπσβ)|t1|^β)/β]`
  expanded as a finite double sum with branch-resolved logs.
  The old code WARNED and produced wrong results here
  (DiffExp/RegularizedIntegration.m:507-530 "branch phases on the negative
  side are NOT applied"); the new module computes it exactly.

#### 2.2.5 Boundary terms at nonzero endpoints

`F(T)` at `T > 0`: principal real `Log T`, real powers, evaluated at ≥ 2×WP
under `$MinPrecision = WP`. `F(t)` at `t < 0` (case-3 windows on the
negative side, and the b ≠ 0 interior arms): the chart's resolved branch
sign σ — `σ := SectorSeries`ChartImSign[ls]`; E3 fires iff it returns None
while a phase is needed, or throws conflict (REVIEW-minimalism 18) — gives
`Log t = Log|t| + iπσ`,
`t^(m+b·eps) = e^(iπσ(m+b·eps))|t|^(m+b·eps)` — note the phase carries the
FULL exponent `a + n + b·eps` including non-integer a (math review
finding 5(i); old half-integer special-casing,
DiffExp/AnalyticContinuation.m:74-77). A `b = 0` sector with non-integer `a`
on a negative arm also requires σ (E3 if missing). If the evaluated boundary
term is non-numeric or exceeds `10^(WP/2)` in magnitude, raise E9b — the old
silent Pade rescue (DiffExp/RegularizedIntegration.m:1845-1870) is forbidden
(F5).

### 2.3 `IntegrateSegmentedLine[segments, {xLo, xHi}]`

Assembly of a chart chain into the final eps-Laurent.

Arguments:
- `segments` — list of IntegrationSegment (section 3.3), in main-line order.
- `{xLo, xHi}` — main-line bounds, `xLo < xHi`, each either exactly a
  segment boundary / chart-center image or strictly inside a segment.

Algorithm:
1. Validate every chart is affine (E1) and every segment's LocalSolution has
   unique sector tags (E7).
2. Tiling check: the segments' MainWindows must exactly tile `[xLo, xHi]`
   after clipping — successive boundaries equal within `snapTol`, no gaps,
   no overlaps, full coverage (E5). Snap-or-error: a boundary within
   `snapTol` of a neighbor's is snapped TO it; outside `snapTol` it is E5
   (no silent passthrough — forbidden F14; old silent snap fallback,
   DiffExp/RegularizedIntegration.m:2388-2398, and FT-side
   `snapMainExpression`, FeynmanTrick/DiffExpIntegration.m:151-178).
3. Per segment: map the clipped main window to chart coordinates through the
   affine map `x = x0 + s·t`: `t-window = sorted[{(x−x0)/s}]`, jacobian
   contribution `∫ f dx = |s| · ∫_{tlo}^{thi} f dt` — EXACT constant `|s|`,
   never evaluated at a midpoint (the old midpoint jacobian,
   DiffExp/RegularizedIntegration.m:2186-2194, is correct only because maps
   are linear; Mobius rejection makes the exact rule total). Call
   `IntegrateLocalSolution`.
4. Interior-crossing bookkeeping: if a chart center's main-line image lies
   strictly inside `(xLo, xHi)`, the segment(s) covering it must present
   Integrate with a t-window straddling 0 in ONE IntegrateLocalSolution call
   (2.2.1 case 2). If the transport layer split the chart into two
   half-segments (lower/upper transport concatenation,
   FeynmanTrick/DiffExpIntegration.m:489-501), the assembler MERGES the two
   half-windows of the same chart (same `Center`, same LocalSolution
   identity — assert) before integrating; an unmatched half whose inner
   boundary is a chart center NOT equal to `xLo`/`xHi` is E6. This is the
   "two-half-segment pairing ENFORCED by assert".
5. Sum per component via EpsSeries add. Final window:
   `Min = min_seg Min_seg`, `CompleteMax = min_seg CompleteMax_seg` — honest
   min, no padding (F10). ErrorEstimates add across segments.
6. Result-level leading-zero trim: leading coefficients are trimmed only by
   the RELATIVE `laurentLeadTol` test from Tolerances.m (the FT absolute
   `10^-40` test, FeynmanTrick/DiffExpIntegration.m:1177-1190, is forbidden
   F12). Interior zeros are never trimmed.

Returns `LineIntegralResult` (section 3.5).

Env-gated diagnostics (campaign-proven, keep): if
`DIFFEXP2_DUMP_INTEGRATE_DIR` is set, dump the full input
(segments + bounds) to a numbered file before integrating, replayable by a
one-shot script (old pattern: `maybeDumpLaurentDefiniteIntegral`,
DiffExp/RegularizedIntegration.m:225-256 + `Scripts/eval_dump_generic.m`;
memory brief "dump-replay loop"). ~25 lines; prints one line per dump.

### 2.4 `EndpointSectorLimit[ls, direction]`

    EndpointSectorLimit[ls_Association, direction_:1]
      -> {EpsSeries value per component (3.4)}

The RewritePlan 3.3 limitUpper/limitLower primitive, owned HERE because it
shares the 2.2.3 gate (REVIEW-minimalism 12): `lim_{t→0}` of the assembled
LocalSolution. Per sector:

- `b ≠ 0` (any a, any p): dropped EXACTLY — the dimreg convention of 2.1.2
  applied to the value, `t^(b·eps) Log^q t |_{t=0} := 0`.
- `b = 0` content that diverges at t = 0 — any monomial `t^m` with
  `m = a + n ≤ −1` (any p), or `t^0` content (`a + n = 0`) with `p > 0` —
  runs the 2.2.3 merged-coefficient cancellation gate per eps order; a
  surviving coefficient is a LOUD ERROR (E2-class; API.md E26 payload:
  component, tag, eps order — REVIEW-math D13).
- The limit readout is the coefficient of the monomial `t^0`: Σ over sectors
  with `b == 0`, `p == 0` and integer `a ≤ 0` of `c[k, n = −a, comp]` (exact
  tag selection; on the canonical merged input — I2 — at most one such
  sector exists per component). NOT "the (0,0,0) sector's constant": after
  the integer-spaced merge the former (0,0,0) content may sit at column
  n = −a of a merged a < 0 sector (REVIEW-math D13).

`direction` resolves nothing here (the limit value is branch-independent for
the surviving terms) but is validated against the chart Prescriptions for
consistency: E3 fires if the chart is multivalued AT ALL (`b ≠ 0` OR `p > 0`
OR `Denominator[a] > 1` — DEC-16) and no prescription exists — the check
stays so configuration gaps surface even on the drop path. API.m's
EndpointLimit = combination (SectorSeries CombineLocalSolutions /
MultiplyRational) + this function (API.md §7).

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution (RewritePlan 3.1, verbatim)

    LocalSolution = <|
      "Center" -> exact x0, "ChartMap" -> affine (DEC-18: Mobius is dropped
                  from the new core entirely; RoC rescaling is an affine
                  rescaling and folds into Scale — see 3.2),
      "Radius" -> distance to nearest singularity IN THE COMPLEX PLANE
                  (complex singularities are real: pentagon/unequal-mass
                  lines have them; old code projects ghosts Re, Re±Im —
                  ledger item; new code uses true complex distance),
      "Sectors" -> { Sector.. },
      "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
      "TWindow"   -> <|"CompleteMax" -> nmax|>,   (* truncation-order record
                  ONLY (DEC-9): recursion matrices are exact polynomials, so
                  the old coupling-depth degradation has no premise in the
                  new core; the MaxCouplingOrder lesson is subsumed by exact
                  polynomial recursion + ErrorEstimate *)
      "ErrorEstimate" -> per (eps-order) accumulated error (two-point
                  full-vs-reduced-order probe, additive across segments,
                  abort > 1 — ported old machinery, user-facing),
      "Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
                  "LeadingCoeffSign"|> with the derived chart Im-sign;
                  consistency-checked at construction (even multiplicity = no
                  constraint; conflict or missing prescription at a chart
                  that is multivalued AT ALL — b != 0 OR p > 0 OR
                  Denominator[a] > 1, DEC-16 — = LOUD ERROR;
                  prescription-factor dedup is SIGN-AWARE per DEC-16; sqrt
                  factors auto-prescribed as in old State.m DEqnSquareRoots)
    |>

### 3.2 Sector (RewritePlan 3.1, verbatim)

    Sector = <|
      "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
      "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT:
                  solutions are vectors; matching and ODE residuals need all
                  components *)
    |>
    Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                    Σ_k eps^k Σ_n c[k,n,comp] t^n

Integrate.m's reading of the windows: `EpsWindow` indexes the eps-order k of
`Coeffs` (k ∈ [Min, CompleteMax]; below Min coefficients are exactly zero;
above CompleteMax unknown — RewritePlan 3.1 EXEMPTION: homogeneous
true-resonance sectors have `Min = −p`). `TWindow["CompleteMax"]` is the
highest trustworthy t-power n. Integrate consumes ONLY k ≤ CompleteMax and
n ≤ TWindow CompleteMax.

`ChartMap` canonical affine form (affine-only is now library-wide — DEC-18
drops Mobius from the new core; field-name reconciliation with
Transport.md remains, Q2):

    "ChartMap" -> <|"Type" -> "Affine", "Center" -> x0, "Scale" -> s|>
    meaning x = x0 + s·t,  s ≠ 0 real (includes any RoC rescaling).

Any other `"Type"` (e.g. "Mobius") ⇒ E1.

### 3.3 IntegrationSegment (input to 2.3; produced by Transport.m, threaded
by API.m's IntegrateOverLine)

    IntegrationSegment = <|
      "LocalSolution" -> LocalSolution,
      "MainWindow" -> {xL, xR}      (* xL < xR, main-line coordinates;
                                       the sub-interval this chart covers *)
    |>

### 3.4 Eps-Laurent values (owned by Docs/specs/EpsSeries.md; DEC-13)

Every eps-Laurent value this module produces or consumes is the canonical
EpsSeries object, verbatim (DEC-13; REVIEW-math D18 / REVIEW-minimalism 23):

    <|"EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
      "Coeffs" -> {c_kmin, ..., c_kmax}|>

with add / scale / multiply / Laurent-divide carrying windows per the
EpsSeries window-shift semantics (denominator with vanishing eps^0 part
shifts Min AND CompleteMax). Metadata ("ErrorEstimate",
"CancelledDivergent", tail estimates) lives in SIBLING association keys
(2.2 return shape, 3.5), never inside the series object. No other in-core
shape exists; API.md's user-facing LaurentValue is produced exclusively by
an API.m boundary converter — this module never constructs or consumes it.

### 3.5 LineIntegralResult (output of 2.3)

    LineIntegralResult = <|
      "Values" -> {EpsSeries value per component (3.4)},
      "EpsWindow" -> <|"Min"->.., "CompleteMax"->..|>,
      "ErrorEstimate" -> <|k -> est ..|>,
      "CancelledDivergent" -> n1, "FinitePartDrops" -> n2,
      "SegmentCount" -> n
    |>

Consumers (API.m / FT layer) must check `need ≤ CompleteMax` and fail loudly
naming (chart, sector, order) — RewritePlan 3.1 invariant. Integrate's
output makes the FT-side `ShiftRawBoundariesToFinite` zero-padding
(FeynmanTrick/DiffExpIntegration.m:1258-1280) and the "keep at least one
output order" floor (FeynmanTrick/DiffExpIntegration.m:658-661) both
unnecessary and forbidden (F8, F10).

---

## 4. INVARIANTS (always on, cheap)

- **I1 Affine-only:** every consumed chart has `ChartMap["Type"] === "Affine"`
  with numeric nonzero Scale. Violation ⇒ E1 (assert at entry of 2.2 and 2.3).
- **I2 Unique tags, merged canonical form:** within each LocalSolution,
  sector tags `(a,b,p)` are pairwise distinct under exact comparison, AND no
  two sectors share (b, p) with integer a-difference (the integer-spaced
  merge of SectorSeries.md 2.2 has been applied — REVIEW-math D12).
  Violation ⇒ E7. (Cancellation correctness depends on this: merging is
  summation, and the 2.2.3 gate's per-monomial coefficients are well-defined
  only on the merged canonical form.)
- **I3 Exact cell selection:** `b == 0` and `a+n+1 == 0` decided on exact
  tags only; assert tags are exact (no inexact numbers in a, b, p).
- **I4 Pairing/window assert:** after any interior crossing, the assembled
  EpsWindow equals the static no-enhancement account (the paired closed
  forms contribute no shift); after any chart-level integral, the
  EpsSeries-computed window equals the static per-sector account of 2.1.3.
  Mismatch ⇒ internal error (E12), never a quiet window adjustment.
- **I5 Window honesty:** output windows are produced exclusively by EpsSeries
  ops + the min-rule of 2.3 step 5; no literal window construction anywhere
  else in the module (code-review checkable).
- **I6 Additivity spot-check:** per IntegrateSegmentedLine call, pick one
  random segment and one random interior point v of its t-window; assert
  `∫[u,v] + ∫[v,w] = ∫[u,w]` to relative `matchTol` on the in-window orders.
  Catches branch-rule and pairing inconsistencies for ~1 extra chart
  integral per line.
- **I7 Euclidean reality:** if all input coefficients are real and no
  crossing phase was applied (no negative-arm b≠0/non-integer-a content),
  assert `|Im| ≤ tol·(1+|Re|)` on every output coefficient (the 9aeb300
  complex-leak regression, Tests/test_interior_singular_integration.m:39,51,64).
- **I8 Bounds-in-chart:** `|t| < Radius` for every integration bound
  (E9a); record `max |t|/Radius` in verbosity output.
- **I9 t-tail recorded:** every chart integral adds its t-truncation
  estimate to ErrorEstimate (2.2.2); estimate > 1 aborts (E8).

---

## 5. ERROR CONTRACT

All errors are raised via ``Tolerances`DE2Error[id, payload]`` (DEC-1:
prints a one-line summary and
`Throw[Failure["DiffExp2", payload], "DiffExp2Error"]`; the catch sits at
every API.m entry point); this module defines no other error mechanism
(REVIEW-math D17 / REVIEW-minimalism 5). The payload carries, ALWAYS: "ID",
"Module" -> "Integrate", the entry point name, the chart Center x0
(main-line), the integration bounds (main-line and chart coordinates), and
the component index; plus the per-error fields below. No error in this
module is catchable-and-defaulted internally (no `Quiet[Check[..., 0]]`
anywhere — the FT-layer instances at FeynmanTrick/DiffExpIntegration.m:
1037-1059 and 1080-1093 are the named anti-pattern, F9).

- **E1 MobiusChartRejected.** Fires when any consumed `ChartMap` is not
  affine. Message: chart Center, the map's Type, and the remediation "Mobius
  charts are dropped from the new core entirely (DEC-18); this segment came
  from the legacy pipeline or a foreign producer — re-transport with the
  DiffExp2 core (affine charts only; RoC rescaling is affine and stays)".
  Same contract as the old hard requirement
  (FeynmanTrick/DiffExpIntegration.m:352
  `UseMobius -> False, (* Required for integration! *)`;
  DiffExp/RegularizedIntegration.m:8-10), now loud instead of conventional.
- **E2 NonIntegrableB0Divergence.** Fires per 2.2.3 step 4: a `b = 0` sector
  with `a+n+1 ≤ 0` (ANY p, including the `a+n+1 = 0`, p ≥ 0 log classes)
  has a merged coefficient failing the RELATIVE `laurentLeadTol`
  cancellation test of 2.2.3 step 3 (DEC-4; REVIEW-math D3 — never the
  absolute chop floor) at an integration endpoint
  at the chart center. Message: chart Center, sector tag (a, b, p), local
  power n, absolute power m = a+n, eps order k, the coefficient value, the
  relative scale and threshold used, which bound (Lower/Upper) hit the
  center, and the hint "divergent meromorphic content must cancel in the
  assembled master combination; check IBP combination and upstream windows".
- **E3 MissingOrConflictingPrescription.** Fires when a negative-arm
  boundary term or an interior b ≠ 0 crossing needs the chart Im-sign and
  `Prescriptions` yields none or a conflict (cf. RewritePlan 3.1
  Prescriptions; old AnalyticContinuationFailed class,
  DiffExp/AnalyticContinuation.m:45-68). Message: chart Center, the factors
  and signs found, the sector tags requiring the phase, eps order range
  affected. (The pentagon "current point is not recognized as a branch
  point... add DeltaPrescriptions" warning class becomes THIS error —
  configuration gaps fail loudly, RewritePlan section 2.)
- **E4 EpsWindowExceeded.** Fires when a consumer-facing request (from
  API.m) needs order K > output CompleteMax, or when an input
  LocalSolution's window is too small to produce ANY requested order.
  Message: chart, sector tag responsible for the binding shift (the
  arg-min of the CompleteMax account), requested vs available order.
- **E5 SegmentTilingError.** Gaps, overlaps, unsorted segments, bounds
  outside coverage, or snap failures (boundary mismatch > snapTol).
  Message: the two offending boundary values, their difference, snapTol,
  and both charts' Centers.
- **E6 UnpairedHalfSegment.** A chart-center main-line image lies strictly
  inside (xLo, xHi) but the assembler could not merge both half-windows
  (missing partner, or two candidate partners with different LocalSolution
  identity). Message: chart Center, the half-window found, expected partner
  window.
- **E7 DuplicateSectorTags.** Input LocalSolution has two sectors with
  identical exact (a,b,p). Message: chart Center, the tag, both sectors'
  leading coefficients.
- **E8 ErrorEstimateAbort.** Accumulated relative error estimate > 1 at any
  in-window eps order. Message: chart, eps order, estimate, contributions
  breakdown (incoming vs t-tail).
- **E9a BoundOutsideChart** (`|t| ≥ Radius`): chart Center, Radius, the
  offending bound (both coordinates). **E9b BoundaryEvaluationBlowup**
  (non-numeric or `|F| > 10^(WP/2)`): chart Center, t, the monomial (m,b,p),
  the value's magnitude — and NO Pade or re-summation rescue (F5).
- **E10 MalformedTag.** Non-exact a/b, negative or non-integer p, inexact T
  classification, float m reaching SectorMonomialIntegral. Message: the
  offending value and its Head.
- **E11 EmptyOrMalformedInput.** `segments === {}`, missing keys, Coeffs
  dimensions inconsistent with windows. A LocalSolution representing exact
  zero must say so explicitly (zero-dimension Coeffs with valid windows
  integrates to exact zero WITH those windows); a missing/empty decomposition
  is NOT silently treated as zero (old behavior at
  DiffExp/RegularizedIntegration.m:2207-2210, 2833-2836 — forbidden F7).
- **E12 InternalWindowMismatch.** Invariant I4/I5 violation — always a bug;
  message includes both window accounts.

### 5.1 Forbidden fallbacks — enumerated and banned

Every site where the old code degraded silently (or warned and continued)
maps to a hard rule here. Implementers: if you find yourself writing any of
the following, stop — it is a spec violation.

- **F1** Dropping resonant/divergent endpoint coefficients ("zero epsilon
  regulator; dropping contribution",
  DiffExp/RegularizedIntegration.m:1470-1477) → E2 after the 2.2.3 check.
- **F2** Unit-regulator substitution `b := 1` for divergent b = 0 content
  (DiffExp/RegularizedIntegration.m:654-660) → forbidden; b = 0 content is
  integrated by elementary antiderivatives only (no eps mixing —
  RegularizedIntegration.m:562-570's own lesson).
- **F3** Averaged-exponent explicit-log integration when sector structure is
  "unresolved" (DiffExp/RegularizedIntegration.m:1599-1631) → cannot arise
  (tags are exact data); any code path that would need it is E12.
- **F4** Salvage of unresolved offsets / "dropped N endpoint coefficient(s)"
  warnings-and-continue (DiffExp/RegularizedIntegration.m:1565-1596,
  1659-1689) → dies with the fitter; no equivalent exists.
- **F5** Pade rescue at boundary evaluation
  (DiffExp/RegularizedIntegration.m:1845-1870) → E9b.
- **F6** Warn-and-drop on non-integer eps powers / IntegerPower snapping
  (DiffExp/RegularizedIntegration.m:705-741, 273-284) → E10.
- **F7** Empty decomposition treated as zero integrand
  (DiffExp/RegularizedIntegration.m:2207-2210, 2833-2836) → E11.
- **F8** "Keep at least one output order" floor
  (FeynmanTrick/DiffExpIntegration.m:658-661) → windows may be empty; report
  them empty.
- **F9** `Quiet[Check[..., 0]]` around boundary/coefficient evaluation
  (FeynmanTrick/DiffExpIntegration.m:1037-1059, 1080-1093) → E9b/E10.
- **F10** Zero-padding windows to a common shape
  (FeynmanTrick/DiffExpIntegration.m:1258-1280 `ShiftRawBoundariesToFinite`)
  → honest min-windows only.
- **F11** Real-log split silently applied to b ≠ 0 interior crossings (the
  old warning-and-wrong path, DiffExp/RegularizedIntegration.m:507-530) →
  the exact phase-paired form of 2.2.4, or E3 when the prescription is
  missing.
- **F12** Absolute 10^-40 leading-coefficient trim
  (FeynmanTrick/DiffExpIntegration.m:1177-1190) → relative laurentLeadTol.
- **F13** Jacobian sampled at an interval midpoint
  (DiffExp/RegularizedIntegration.m:2186-2194) → exact affine |s|.
- **F14** Silent unsnapped-coordinate passthrough
  (DiffExp/RegularizedIntegration.m:2388-2398) → snap-or-E5.
- **F15** Iteration caps with warning-and-continue
  (RegularizeIntegrand maxIter, DiffExp/RegularizedIntegration.m:450-464) →
  no iteration exists in the new design; any loop bound that can be hit is
  an error.
- **F16** Per-eps-order processing that drops cross-order content (the
  collapsed-tower disease D2 generally) → all arithmetic is symbolic-eps
  EpsLaurent; there is no per-order code path to fall back to.

---

## 6. ABSORBED OLD CODE

This module + SectorSeries.m absorb the consumption side of
`DiffExp/RegularizedIntegration.m` (3015 lines, 15 exports; export list at
RegularizedIntegration.m:45-76). Disposition of every export:

| old export (file:line of ::usage) | disposition |
|---|---|
| `ApplyRegularizationStep` (RegularizedIntegration.m:45; impl 363-444) | DIES. The IBP regularization step is subsumed by the dimreg boundary rule applied to the closed form (cells 11/12). |
| `RegularizeIntegrand` (:47; impl 446-466) | DIES (same reason; its maxIter warning is F15). |
| `IntegrateSingularTerm` (:49; impl 482-487) | ABSORBED → `SectorMonomialIntegral` + `IntegrateLocalSolution` (closed forms, no regularization loop). |
| `IntegrateSingularTermLaurent` (:51; impl 489-703) | ABSORBED → same. Its three internal paths map: interior-split (507-530) → 2.2.4; meromorphic Hadamard split (571-652) → 2.2.4 b = 0 closed forms; subtraction dispatch (668-675) → dimreg rule. |
| `IntegrateDecomposition` (:53; impl 1880-1888) | DIES (the DecomposeSingularity input format is replaced by exact Sector data from SectorSeries). |
| `IntegrateDecompositionLaurent` (:55; impl 1890-1919) | ABSORBED → `IntegrateLocalSolution` (per-sector sum with EpsLaurent assembly). |
| `EvaluateLimitAtSingularity` (:57; impl 1921-1952) | ABSORBED → `EndpointSectorLimit` (2.4): the b ≠ 0 drop rule and the divergence gate live HERE (REVIEW-minimalism 12); API.m's EndpointLimit = SectorSeries combination + this function (API.md §7). |
| `EvaluateEndpointLimitSectors` (:59; impl 1965-2141) | DIES WITH THE FITTER. Its contract (drop rule applied per sector to the ABSOLUTE exponent, not to a collapsed exponent — RegularizedIntegration.m:1954-1964) is automatic in the sector-native representation. FT call site FeynmanTrick/DiffExpIntegration.m:1015 retargets to the API limit entry (shim audit R4). |
| `FitResidualEndpointSectors` (:61; impl 909-1331, plus clusterRootSpecs 860-907) | DIES WITH THE FITTER (Prony/Hankel recovery, candidate dominance gates, salvage). ~470 lines deleted by design; the N-root fitter (f48cd94) is superseded (RewritePlan section 10). |
| `IntegrateSegmentData` (:63; impl 2150-2211) | ABSORBED → per-segment step of `IntegrateSegmentedLine` (with F13 fixed: exact affine jacobian). |
| `IntegratePiecewiseSaved` (:65; impl 2218-2276) | ABSORBED → `IntegrateSegmentedLine` (overlap selection becomes the strict tiling contract E5). |
| `DefiniteIntegral` (:67; impl 2279-2282) | ABSORBED → thin API.m wrapper (IntegrateOverLine) over `IntegrateSegmentedLine`. |
| `IndefiniteIntegral` (:69; impl 2285-2353) | DROPPED v1 (DEC-24: API exposes IntegrateOverLine only; the 6 IndefiniteIntegral test sites are retargeted/retired at M6, when the disposition table is updated). No FT or example dependency; API.m's ToPiecewise-equivalents cover the use case. Record in Docs/ExportDisposition.md. |
| `DefiniteIntegralWithPrefactor` (:71; impl 2487-2547) | DIES (non-Laurent variant; superseded). |
| `DefiniteIntegralWithPrefactorLaurent` (:76; impl 2549-2625 + segment workers 2632-3010) | SPLIT: the power-law/rational prefactor multiplication and the pole-absorption-into-a (impl 2680-2838) are SectorSeries.m's rational-multiply (a-shift at the endpoint charts is exact tag arithmetic); the integration core is `IntegrateSegmentedLine`. The FT call site FeynmanTrick/DiffExpIntegration.m:875 gets the named replacement: API.m `IntegrateOverLine[segments, {0,1}, "Prefactor" -> <|"PowerAtLower" -> v1−1, "PowerAtUpper" -> v2−1|>]` = SectorSeries multiply then this module (shim audit R4). v1−1, v2−1 are INTEGER (FT propagator exponents): per DEC-22, v1 prefactors are rational in x times integer powers of eps only — eps-DEPENDENT closed-form prefactors (the Beta-function pin class) are handled at the FT layer, which already carries EpsPrefactors as integer shifts; an eps-dependent power is a loud API-level error (API test 17 replaced accordingly). The Gamma prefactor stays in the FT layer (FeynmanTrick/DiffExpIntegration.m:570-571). |

Supporting internals that DIE with their callers: the Laurent helper family
(LaurentAdd/Scale/Trim/Zero, RegularizedIntegration.m:176-223 — EpsSeries
owns this), `linearSeriesCombine` (119-174), `KnownSeriesDerivative`
(313-341), `EffectiveZeroExprQ` theta-branch probing (286-308 — Tolerances
owns zero-classification), `thetaRulesAtPoint`/`localSidePhase`/
`branchDirection` (1772-1806 — replaced by the Prescriptions-driven branch
rule shared with SectorSeries), `EvaluateIntegralAtPoint` (1809-1873 —
replaced by 2.2.5), `FiniteEndpointConstant` (1694-1699 — replaced by the
dimreg rule), segment/coordinate helpers (2355-2424 — replaced by the exact
affine contract).

### 6.1 Numerical lessons that MUST be preserved (read from the old code)

1. **Combine before integrating.** Individual IBP-weighted terms have
   endpoint poles that cancel only in the sum
   (FeynmanTrick/DiffExpIntegration.m:16-21, 532-545). New form: input
   precondition + the 2.2.3 object-level check.
2. **Never regulate regular powers of meromorphic terms.** The unit
   regulator leaks `c1(B Log B − A Log|A| − (B−A))` from eps^-1 into eps^0
   (old comment DiffExp/RegularizedIntegration.m:562-570; fixed by the
   Hadamard path in 9aeb300). New form: b = 0 content uses elementary
   antiderivatives exclusively; no formal eps ever introduced (F2).
3. **Interior-split real-log convention.** The negative arm of a
   meromorphic crossing must use Log|t|; complex Log leaks iπ (and −π² at
   the next order) (DiffExp/RegularizedIntegration.m:1833-1843;
   Tests/test_interior_singular_integration.m; commit 9aeb300,
   Docs/FeynmanTrickBoxFamilyStatus.md "Fixed this campaign" item 1). New
   form: 2.2.4 b = 0 closed forms + reality assert I7.
4. **The drop rule is per-sector on the ABSOLUTE exponent.** Collapsed
   exponents discard or keep whole towers wholesale — the D2 disease
   (DiffExp/RegularizedIntegration.m:1954-1964 comment; pentagon record,
   RewritePlan section 10). New form: tags exact, rule trivially per-sector.
5. **Eps-lookahead is consumed by endpoint poles.** Old
   `IntegrationPoleAllowance` default 4
   (FeynmanTrick/DiffExpIntegration.m:1421-1436;
   Docs/FeynmanTrickBananaStatus.md "Fixes Implemented" item 5) and the
   carry-full-depth runner lesson (Scripts/run_ft_stepwise.m:154-161;
   box_bubble needed 9 at L2, 11 at L1). New form: the exact window shift of
   2.1.3 + EpsWindow propagation; the static budget (RewritePlan 3.4) must
   dominate, and E4 is the backstop.
6. **Precision discipline at boundary evaluation.** Evaluate Log T / T^α
   from inputs raised to 2×WP under `$MinPrecision = WP`
   (numericAtActivePrecision pattern, DiffExp/RegularizedIntegration.m:82-112;
   ledger seed "input precision raise").
7. **Snap segment boundaries, exactly.**
   (DiffExp/RegularizedIntegration.m:2388-2398;
   FeynmanTrick/DiffExpIntegration.m:151-178) — new form: snap-or-E5, snapTol
   from Tolerances.m.
8. **Relative, not absolute, zero tests on Laurent leading coefficients**
   (FeynmanTrick/DiffExpIntegration.m:1177-1190 → laurentLeadTol; D3 cure).
9. **Dump-replay debuggability.** DIFFEXP_DUMP_LAURENT_DIR +
   Scripts/eval_dump_generic.m gave seconds-per-iteration debugging on the
   box campaign (memory brief). Keep the equivalent hook (2.3).
10. **Trim incomplete top orders, never emit them.** Partially combined top
    eps orders poisoned the fitter (FeynmanTrick/DiffExpIntegration.m:636-667;
    Docs/FeynmanTrickBananaStatus.md "Root cause 1"). New form: honest
    CompleteMax arithmetic makes "incomplete but present" unrepresentable.

---

## 7. DEPENDENCIES

May call (acyclic order Tolerances < Config < EpsSeries < SectorSeries <
Indicial < Solve < Transport/Integrate < API):

- **Tolerances.m** — laurentLeadTol (cancellation gate 2.2.3 + result trim
  2.3 step 6 — REVIEW-math D3/D16, REVIEW-minimalism 1/17), snapTol
  (tiling/endpoint snapping, E5; computed-value snapping only), matchTol
  (additivity spot-check I6), and `DE2Error` (section 5). NO chopFloor: it
  is the absolute coefficient noise floor, never a classification tolerance.
  No literal tolerance constants in this module.
- **Config.m** — WorkingPrecision, Verbosity via the validated accessor
  only; ``Config`PrintInfo`` for the gated info prints (2.2.3 step 3 —
  REVIEW-minimalism 5; no other print helper is defined here).
- **EpsSeries.m** — all eps-Laurent arithmetic incl. Laurent division
  window-shift semantics (cells 9/10) and window-carrying add/scale.
- **SectorSeries.m** — LocalSolution structural accessors/validators and the
  SHARED branch rule used for negative-arm evaluation (2.2.5 must be the
  same function or a re-export; divergence between Integrate's and
  SectorSeries' branch resolution is a spec violation).

MUST NOT call: Indicial.m, Solve.m, Transport.m (same level — segment data
arrives as input), API.m. MUST NOT read FT-layer state.

---

## 8. UNIT TESTS

File: Tests/diffexp2/test_integrate.m (kernel-run by the orchestrator).
`L := Log[T]` below; all numeric asserts at relative 10^-25 unless noted,
WP = 80. "window" asserts check Min and CompleteMax exactly.

1. **monomial_cell1_plain_power** — `SectorMonomialIntegral[2, 0, 0, 1/3, 4]`
   = `(1/3)^3/3 = 1/81` at eps^0 only; window [0,4].
2. **monomial_cell2_b0_log** — m=1, b=0, p=2, T=1/2:
   `T^2·eps^2·(Log^2 T/(2·2) − ... )` per cell 2 formula; closed form
   `eps^2 · T^2 (2 Log^2 T − 2 Log T + 1)/8`; assert eps^2 coefficient
   equals `(1/4)(2 Log^2 2 + 2 Log 2 + 1)/8` (note Log(1/2) = −Log 2);
   orders ≠ 2 zero; window Min = 2.
3. **monomial_cell7_regular_beps** — m=0, b=2, p=0, T=1/2:
   `T^(1+2eps)/(1+2eps)`; assert eps^0 = 1/2, eps^1 = Log(1/2) − 1,
   eps^2 = Log^2(1/2) − 2 Log(1/2) + 2 (general: T·Σ pattern); window
   Min = 0, no shift.
4. **monomial_cell9_single_pole** — m=−1, b=b0, p=0, T:
   `T^(b0·eps)/(b0·eps)`; assert eps^-1 = 1/b0, eps^0 = L,
   eps^1 = b0·L^2/2, eps^2 = b0^2·L^3/6; window Min = −1 and
   CompleteMax = kMaxOut (primitive level), and at sector level
   CompleteMax = κmax − 1 (the shift).
5. **monomial_cell10_p1_pole** — m=−1, b=b0, p=1:
   assert {eps^-1: −1/b0^2, eps^0: 0, eps^1: L^2/2, eps^2: b0·L^3/3}
   (derived: `eps·∂_β[T^β/β]`, β=b0·eps). Window Min = −1.
6. **monomial_cell10_p2_pole_depth3** — m=−1, b=b0, p=2 (the depth-(p+1)
   pin and the M3 `∫ x^(−1) log^k x` resonant-source class):
   assert {eps^-1: 1/b0^3, eps^0: 0, eps^1: 0, eps^2: L^3/6,
   eps^3: b0·L^4/8, eps^4: b0^2·L^5/20}
   (derived: `(eps^2/2)·∂_β^2[T^β/β]`). The zero coefficients at eps^0 and
   eps^1 are exact cancellations — assert EXACT zero, not just small.
7. **monomial_cell11_dimreg_negative_power** — m=−5/2, b=1, p=0, T=1/4:
   NO error; `F(0)=0`; assert eps^0 = T^(−3/2)/(−3/2) = −(2/3)·8 = −16/3;
   window Min = 0, no shift (the half-integer m also exercises exact
   non-integer tag handling — no snapping).
8. **error_cell3_b0_log_divergence** — m=−1, b=0, p=0, nonzero coefficient,
   endpoint at center: assert error E2 fires and its payload contains chart
   center, (a,b,p), n, k.
9. **error_cell4_b0_logchain_divergence** — m=−1, b=0, p=2: E2 fires (the
   "any p" requirement, finding 6(ii)).
10. **cancellation_gate_exact_and_chop** — three sub-cases on a merged
    (a=−1, b=0, p=0) sector at the lower endpoint: coefficient exactly 0
    (integrates cleanly, `CancelledDivergent` = 0); coefficient 10^-60 with
    scale O(1) at WP 80 (dropped, `CancelledDivergent` = 1, no error);
    coefficient 10^-3 (E2 fires).
11. **interior_pv_residue** — replicate
    Tests/test_interior_singular_integration.m:41-51 in the new
    representation: sector (a=−1, b=0, p=0), Coeffs c[0,0]=1, c[0,1]=1,
    bounds {−1/5, 1/4}: assert eps^0 = Log(B/|A|) + (B−A) with A=−1/5,
    B=1/4; Im exactly below 10^-25; other orders zero; window unshifted.
12. **interior_fp_double_pole** — sector (a=−2, b=0, p=0), c[0,0]=1,
    c[0,2]=1, bounds {−1/5, 1/4}: assert eps^0 = (−1/B + 1/A) + (B−A),
    real (Tests/test_interior_singular_integration.m:53-64);
    `FinitePartDrops` ≥ 1 recorded.
13. **interior_pv_logchain** — sector (a=−1, b=0, p=1), c[0,0]=1, bounds
    {−A0, B0}: assert value = eps·(Log^2 B0 − Log^2 A0)/2 (2.2.4 closed
    form with p=1: `eps^p(Log^(p+1)t2 − Log^(p+1)|t1|)/((p+1)p!)`), real.
14. **interior_idelta_pairing** — sector (a=−1, b=b0, p=0), c[0,0]=1,
    prescription σ=+1, bounds {−A0, B0}: assert eps^-1 coefficient is
    EXACTLY absent (window Min = 0 — the pairing cancels the pole and the
    window is NOT shifted, finding 7(v)); eps^0 = Log(B0/A0) − iπ;
    eps^1 = b0(Log^2 B0 − (Log A0 + iπ)^2)/2. Repeat with σ=−1 (conjugate
    phase). Missing prescription ⇒ E3 fires.
15. **pairing_assert_unpaired_half** — segment chain where a chart center
    image lies strictly inside (xLo,xHi) but only the upper half-window is
    supplied: E6 fires.
16. **mobius_rejected** — segment with ChartMap Type "Mobius": E1 fires
    naming the chart.
17. **window_honesty_pole_shift** — LocalSolution with two sectors,
    (0,0,0) and (−1, b0, 0), both with Coeffs window [0, K]: output window
    CompleteMax = K−1 and Min = −1; asserting order K available ⇒ E4.
18. **tiling_gap_error** — two segments with a 10^-3 gap: E5;
    with a 10^-50 mismatch (< snapTol at WP 80): snapped, no error.
19. **assembly_linearity_two_charts** — integrand `x` represented on two
    charts tiling [0,1] (centers 0 and 1, affine scales ±1): assembled
    eps^0 = 1/2 exactly; additivity invariant I6 must pass with the random
    split (run the call 3 times).
20. **beta_function_end_to_end** — `∫_0^1 x^(eps−1)(1−x) dx = 1/(eps(1+eps))`
    via two charts: chart at 0 with sector (a=−1, b=1, p=0),
    c[0,0]=1, c[0,1]=−1 over [0, 1/2]; chart at 1 (t = 1−x, scale −1) with
    sector (0,0,0), Coeffs c[k,n] = the eps-expansion coefficients of
    `t·(1−t)^(eps−1)` (eps-DEPENDENT b=0 coefficients), over [1/2, 1].
    Assert assembled orders {eps^-1: 1, eps^0: −1, eps^1: 1, eps^2: −1}
    to 10^-25; window Min = −1; exercises orientation (negative scale),
    pole cell, regular cells, assembly. (This is the FT `integrate`
    boundary case in miniature: prefactor (v1,v2) = (eps-shifted, 2).)
21. **reality_invariant_trips** — feed real data through a path that
    wrongly applies a complex log on the negative arm (test-only hook or
    deliberately corrupted prescription): I7/E-class must fire — guards the
    9aeb300 regression.
22. **ftc_crosscheck_vs_sectorseries** — for a random multi-sector
    LocalSolution (no divergent b=0 content), compare
    `IntegrateLocalSolution[ls, {T', T'+δ}]` against
    `δ·SectorSeriesEvaluate[ls, T'+δ/2]` with Richardson δ, δ/2 to
    O(δ^2)-extrapolated 10^-12 — pins Integrate's branch/normalization
    conventions to SectorSeries' evaluate (NOT always-on; unit only).

Integration-level (not unit; M4/M5 per RewritePlan): interior-pole crossing
on the box L1 line (pole at 1/4 — campaign pin), near-endpoint pole
box_triangle L3 geometry 0.9617, budget-formula reproduction of box_bubble's
9/11 (RewritePlan 3.4 validation unit), banana IntegrationPoleAllowance=4
record dominated.

---

## 9. LINE BUDGET

RewritePlan 3.2 allots **~500 lines**. Estimate:

| block | lines |
|---|---|
| case-table primitive (closed forms, dimreg rule, EpsLaurent assembly) | ~120 |
| chart-level integral (classification, per-sector assembly, cancellation gate, t-tail estimates) | ~140 |
| interior pairing (b=0 PV/FP + b≠0 phase-paired forms) | ~80 |
| segment assembly (tiling, jacobian, merge-halves, min-windows, trim) | ~80 |
| error helper + assert helpers (E1-E12, I1-I9) | ~60 |
| dump hook + verbosity | ~25 |
| total | ~505 |

If over budget, cut in this order (never cut error richness or asserts):
1. Move the FTC crosscheck helper (test 22 support) entirely into the test
   file (−15).
2. Demote the b ≠ 0 interior phase-paired path (2.2.4 last bullet) to loud
   error "b ≠ 0 interior crossing unsupported v1" with a LessonsLedger
   waiver — NO in-scope example needs it (box interior pole is b = 0
   IBP-induced; DE singularities of all 8 examples sit at the integration
   endpoints) (−50). Keep test 14 disabled with a pointer to the waiver.
3. Inline the monomial primitive into the sector loop and expose it only via
   a thin test shim (−20). Last resort: it is the unit-test surface.

What must NOT be traded for budget: the cancellation gate, the pairing
asserts, honest window arithmetic, the twelve-cell completeness, the
structured error payloads.

---

## 10. OPEN QUESTIONS

- **Q1 EpsLaurent field names.** This spec assumes
  `<|"Min", "CompleteMax", "Coefficients"|>` mirroring EpsWindow; EpsSeries.md
  owns the truth. Reconcile at the M0 spec review; pure rename.
- **Q2 Segment input shape.** `IntegrationSegment` (3.3) and the ChartMap
  canonical form are assumed; Transport.md/API.md own what TransportTo
  actually emits (old shape: 5-element segment lists,
  DiffExp/RegularizedIntegration.m:2147-2149). Reconcile field names and who
  merges lower/upper transport halves (this spec places the merge in 2.3
  step 4).
- **Q3 PV vs PV ∓ iπ·residue for b = 0 interior poles.** v1 pins the old
  real-PV/finite-part convention (campaign-validated on the Euclidean box,
  pin to 11 digits). For Minkowski lines a surviving b = 0 interior residue
  arguably should pick up ∓iπ from the chart prescription. Decision
  deferred: v1 = real PV + the I7 reality assert; revisit if an M5 example
  with Minkowski kinematics and a non-cancelling interior residue appears.
  (Note pentagon kinematics are Euclidean-region per Scripts/FTExamples.m.)
- **Q4 Even-m finite-part drops.** 2.2.4 adopts the old Hadamard
  finite-part for m < −1 interior content surviving the chop (old behavior,
  parity-safe). Alternative: treat non-cancelling even-m divergences as E2
  (they likely indicate upstream cancellation failure). v1 keeps FP +
  metadata counter; the M4 multisector closed-form gate should decide
  whether any real example ever exercises it with a nonzero coefficient.
- **Q5 Cancellation chop scale.** 2.2.3 uses
  `scale_k = max_n' |c[k,n',comp]|` of the merged combination. Alternative
  scales (per-sector max; cross-k max) change behavior only for
  pathologically scaled data; Tolerances.md should own the named semantics
  (`chopFloor` relative form) — flag for the Tolerances spec review.
- **Q6 ErrorEstimate granularity.** This spec emits per-eps-order estimates
  and adds them to the incoming LocalSolution estimates. Whether the FT
  layer wants per-(integral, order) maps in the final result shape is an
  API.md question; the data is available either way.
- **Q7 Algebraic (non-rational) a, b performance.** Cells select exactly on
  algebraic tags (R2: correctness fine, performance flagged). If exact
  algebraic `m+1` arithmetic in the geometric expansions is slow on a real
  example, the fix belongs in EpsSeries coefficient representation, not in
  a tolerance — record here so nobody "fixes" it with a numeric snap.
