# Performance Gap Analysis: old DiffExp (~30 s banana) vs DiffExp2 (~25–40 min)

Date: 2026-06-12.  Scope: structural comparison of the two engines on the FT
ladder workload (`Scripts/run_ft_stepwise.m` driving `DiffExp/` vs
`Scripts/run_ft_stepwise2.m` driving `DiffExp2/`).  Both runners carry the
same eps depth (full boundary depth, K up to ~8–9 at the deepest level) and
the same DivisionOrder 4 geometry.  The old engine runs at WP 500 /
ExpansionOrder 50; DiffExp2 at WP 100 / ExpansionOrder 40 — the old engine
is ~50–80x faster while paying ~5x more digits per arithmetic operation.
The gap is therefore **operation count and per-op machinery, not precision**.

Verdict up front: the single biggest structural difference is what each
engine does at a *regular interior chart*.  The old engine propagates the
**value vector** (one small solve per segment, homogeneous data cached
across eps orders, block-by-block); DiffExp2 builds a **full d-column
fundamental system per chart** and then matches weights — d coupled
recursions, an eps-Laurent Gaussian elimination, and a full-array linear
combination, per chart, even where the local structure is trivial.

---

## 1. What the old engine actually does (architecture audit)

### 1.1 Transport = value chaining; chart centers ARE the previous eval points
- `TransportTo` keeps a value object `Currbcs` and per segment runs
  `IntegrateSystem[Currbcs, CurrLine]`, evaluates at `CurrEvalPoint`, and
  *that value becomes the next boundary*: `Currbcs = {line /. x -> CurrEvalPoint, CurrEval}`
  (`DiffExp/Transport.m:1002, 1099–1101, 1171–1174`).  No fundamental
  matrix is ever transported across a boundary.
- In the Predivision strategy the *next regular chart is centered exactly at
  the evaluation point* (`currentCenter = CurrEvalPoint`,
  `DiffExp/Transport.m:1154–1165`; only when a pole is ahead does
  `FindNextCenterPointL/R` move the center to straddle it).  Consequence:
  the boundary data for the next segment sits at the next chart's **t = 0**,
  and `IntegrateSystem`'s `FixAt === 0` branch fixes the integration
  constants by *reading series coefficients at the origin* — no evaluation,
  no linear system beyond a per-block b×b solve
  (`DiffExp/Transport.m:296–307`).
- Step size: evaluation at `±RadiusOfConvergence/DivisionOrder` of the
  rescaled chart (`Transport.m:1041`), i.e. at 1/k of the distance to the
  nearest singularity (charts are rescaled so that distance = 1,
  `Mobius.m:38–69`).  With k = 4, nmax = 50 the per-step truncation tail is
  ~(1/4)^50 ≈ 1e-30.

### 1.2 Expansion-order policy
- **Fixed** `ExpansionOrder` (50 in the FT driver) per segment by default.
- Optional adaptivity exists in two modes, both *tail-driven*:
  `AccuracyGoalValidate -> "Before"` sizes the order from the magnitude of
  the largest expanded-matrix term at the two segment endpoints
  (`Transport.m:774–841`), `"After"` re-runs a segment at +25 orders if the
  decimated-series error probe missed the AccuracyGoal
  (`Transport.m:1191–1213`).  The 30 s banana run uses the fixed order.
- Error estimate "Fast" = re-evaluate the segment with the series order
  decreased by `ICurrEvalErrorSeriesDecrease = ceil(0.7·coupling)+2` and
  diff (`Transport.m:907–924`) — the same idea DiffExp2 already has in
  `SegmentErrorProbe`.

### 1.3 Recursion data layout: eps FIRST, numeric, per block
- The matrices are stored/expanded **per eps order** before any solving:
  `d<var>_<ord>.m` files, or `Series[..., {eps, 0, K}]` of the closed form
  once per line (`MatrixLoading.m:91–111, 323–331`).
- `IntegrateSystem` runs an **eps-order cascade**: at order k the source is
  `Σ_{l≥1} A_l · f_{k−l}` from already-solved lower orders
  (`Transport.m:239–271`).  Every per-order system is **numeric in eps**
  (eps eliminated entirely; coefficients are WP bignum floats).
- The system is split into **diagonal blocks** (sectors) by the coupling
  graph (`MatrixLoading.m:358–385`); blocks solve in lower-triangular
  order, with off-diagonal couplings folded into the per-order source.
  Block sizes for banana/sunrise are 1–3; d never appears as a coupled
  solve dimension.
- Per block, the homogeneous fundamental matrix (b columns of a
  b-dimensional recursion) is computed **once at eps^0 and cached across
  all eps orders** (`Transport.m:276–285`, `Recurrence.m:292–315`,
  LessonsLedger L42).  Each eps order then costs one *particular* recursion
  + a b×b constant fix.
- The recursion itself is the denominator-cleared polynomial recurrence
  (L41): `f_{m+1} = d0⁻¹/(m+1)·(Σ_j A_j f_{m−j} − Σ_i d_i (m−i+1) f_{m−i+1})`
  with `A_j` numeric matrices of the cleared numerator's t-degree dA
  (`Recurrence.m:241–273`), all scalars `N[..., WP]`.
- Per-segment matrix prep is *one rational substitution + one numeric
  series expansion*: `PrepareMatricesFrom1` substitutes the line relation
  into the per-eps-order factored matrices; `PrepareMatricesExpanded` does
  `Series[N[...], {x, 0, 50}]` (`MatrixLoading.m:260–296, 337–355`).
  Centers are arbitrary-precision **floats** — no exact-rational growth.

### 1.4 Cost model, old engine (banana L1: d = 5, blocks ≈ {2,1,1,1}, K = 5, nmax = 50, dA ≈ 8)
- homogeneous (once/segment): Σ_blocks b³·nmax·dA ≈ 11·50·8 ≈ **4.4e3** mults
- per eps order: source assembly ≈ d²·nmax·dA ≈ 1.0e4, particular ≈
  Σ b²·nmax·dA ≈ 3e3, constant fix ≈ b³ each ⇒ ×(K+1) ≈ **8e4** mults
- total ≈ **1e5 bignum(500-digit) multiply-adds per segment** ⇒ ~0.2–0.3 s,
  matching 30 s / ~100 segment+tile units end-to-end.

## 2. Where DiffExp2's cycles go at the SAME regular chart

Per chart, `TransportLine` runs (`DiffExp2/Transport.m:311–396`):
1. `PrepareChart`: exact rational shift `A(x0 + βt)` with `Cancel@Together`
   on d² entries, `ChartIndicial` (trivial for regular, but still d²
   valuation passes), exact `V/VInv` (`Solve.m:62–121`); memoized per
   center, but every *new* center pays it.
2. `prepareCleared`: exact `PolynomialLCM` over d² denominators + exact
   eps-Laurent expansion (`ratEpsList`) of every (t-degree, entry) pair
   (`Solve.m:212–243`).
3. `SolveHomogeneous`: **d columns**, each one a **d-dimensional framed
   recursion** over eps-frames of width `W = K + 3·Pmax + 2·cdMax + 7`
   (`Solve.m:474–479`); at a regular chart Pmax = cdMax = 0 ⇒ **W = K + 7**
   (= 12 at K = 5).  Per column ≈ nmax·dN·nsp·d²·W multiply-adds plus
   O(nmax·d·W²) frame inversions/convolutions in `blockSolveTPFrame`.
4. `ODEResidualCheck` over **all d columns** (d + 1 full evaluations).
5. Matching: d basis evaluations at the match point + `MatchWeights`
   (d³ EpsSeries Gauss steps, each O(W²)) + `CombineLocalSolutions`
   (full (W × nmax × d) array combine per column)
   (`Transport.m:232–271`, `SectorSeries.m:419–458`).

Cost model (banana L1: d = 5, K = 5 ⇒ W = 12, nmax = 40, dN ≈ 8, nsp ≈ 2):
- one recursion column ≈ 40·8·2·25·12 ≈ 1.9e6 (+0.3e6 frame ops) ≈ **2.2e6**
- basis = d columns ≈ **1.1e7**; matching/combine/residual ≈ +0.5e6
- ⇒ ~1.2e7 mixed exact-rational/bignum ops per chart ⇒ the measured
  **2–4 s/chart**.  Op-count ratio vs the old segment: **~100x**; measured
  ratio ~10–15x (the old engine pays ~5–10x per op at WP 500 and benefits
  from SeriesData C kernels).

Chart counts: ~15–30 per line, 2 lines per level, 4 levels ⇒ ~120–240
chart-solves ≈ 4–16 min, plus the per-tile integration phase
(`LineIntegral`: MultiplyRational + IntegrateLocalSolution per nonzero
coefficient per tile, `API.m:73–124`) of comparable magnitude at the deep
levels.

## 3. Ranked levers

### Lever 1 — value-vector propagation through regular interior charts (expected ×4–7 on chart cost; ×2.5–4 end-to-end)  [PROTOTYPED]
**Old-code evidence.** The old engine never transports a basis: `Currbcs`
chaining + centers at eval points + `FixAt === 0` coefficient matching
(§1.1).  The b-column block bases it does build are (a) per *block*, not
per system, and (b) computed once per segment at eps^0.

**The DiffExp2 form.** At a regular chart the transported solution is the
unique Cauchy solution with `f(0) = v`, where `v` = incoming object
evaluated *at this chart's center* — a point that lies well inside the
previous chart's disk by the predivision geometry: consecutive centers are
`step = s/k` apart while the previous radius ≥ `s(1+k)/k`, so
`step/radius ≤ 1/(1+k)` (= 1/5 at k = 4; per-step tail (1/5)^40 ≈ 1e-28,
*better* than the basis path's basis-evaluation at radius/k).  So:
evaluate the previous object at the new center, feed that EpsSeries vector
as `init` into ONE `runRecursion` (the init path already exists — it is
how `SolveHomogeneous` seeds its unit columns), assemble, done.  No
`MatchWeights`, no `CombineLocalSolutions`, no per-column residual checks.

**Invariants preserved.**
- *Exact sector tags*: a regular chart yields the single exact (0,0,0)
  sector (DEC-6) — same as the basis path after combination.
- *Honest eps-windows*: at a regular chart every recursion step n ≥ 1 is
  CASE T with `deltaList = n + 0·eps` — no eps-division, `TopValid` never
  erodes; the delivered window is capped at the incoming `CompleteMax`
  (`capWindow`), exactly the eps-first cascade's window preservation.  The
  frame uses the same deterministic work-window shape as `SolveHomogeneous`
  with Pmax = cdMax = 0.
- *Loud ambiguity errors*: there is no matching solve to go singular; the
  geometric precondition (center inside the previous disk) is asserted by
  `EvaluateLocalSolution`'s radius check, with a conservative 9/10-margin
  pre-check that falls back to the basis path (a performance choice, not a
  silent fix).
- *Always-on ODE residual self-check*: runs on the propagated solution
  itself (`SolveValueRegular` calls `ODEResidualCheck` unconditionally).

**What keeps the basis path**: singular charts (sector decomposition is
*required* there for limits/integration/crossing) and the first chart of a
line (plain-value boundary data is only valid at its anchor, t-window 0).
That is ~1–3 of ~15–30 charts per line.

**Arithmetic.** Per regular chart: basis ≈ d·R + matching ≈ 1.15·d·R;
value ≈ R·(K+5)/(K+7) + one residual check ≈ 0.9·R.  Factor ≈ **1.3·d**
(≈ 6.4x at d = 5, ≈ 9x at d = 7).  Over a line with 90% regular charts:
×4–6 on the transport phase.

**Prototype.** `DiffExp2/Solve.m` `SolveValueRegular` +
`DiffExp2/Transport.m` `TransportLine` branch, gated by
`DE2_VALUE_TRANSPORT=1` (default off; flag-off behavior bit-identical).
Validation plan in §5.

**Risks.** (i) The handoff point moves from matchPt to the center — both
are valid anchorings of the same analytic germ; end values should agree to
MatchTol; pins decide.  (ii) After a singular chart in a tight singularity
cluster the next center can fall outside the singular chart's disk — the
margin pre-check falls back to basis mode for that chart.  (iii) The
`capWindow` trim must keep `ValidateLocalSolution`'s dims invariant
(it slices full eps-rows only).  (iv) Promote to default only after the
oracle pins (bubble/sunrise/banana vs pySecDec references) pass under the
flag.

### Lever 2 — per-line chart-prep amortization: the matrix is the SAME on the whole line (expected ×1.5–2.5 on the post-L1 chart cost)
Only the center changes from chart to chart, yet DiffExp2 redoes per
center: the exact rational shift (d² gcds), the d² valuation scans of
`MatrixPoleData`, the `PolynomialLCM`, and the exact `ratEpsList`
eps-expansion of every cleared coefficient.  The old engine's analogue is
one substitution + one *numeric* series expansion (§1.3, last bullet) —
trivially cheap because its centers are floats.

Implementation sketch respecting invariants:
1. Per line (once): keep the cleared `den(x, eps)`, `num(x, eps)` of the
   *unshifted* system and the exact singular-factor roots
   (`FindSingularities` already computes these).
2. Per regular chart: certify regularity **exactly** — the center is not a
   root of any singular factor (exact comparison against the root list; no
   numerics in the *decision*).  This replaces `MatrixPoleData` on shifted
   entries entirely; the indicial data is the canonical trivial record.
3. Coefficient *values* (not structure): build `dL/NhatL` frame lists by
   Taylor-shifting the per-line coefficient polynomials to the center
   numerically at WP+20 (Horner in x0), never running `Cancel/Together` on
   float polynomials.  Numeric values in coefficients are already the
   steady state (the ByteCount gate numericizes recursion content); the
   I-6 rule "structure decided on exact data" is satisfied because the only
   structural decision — regular vs singular — is made exactly in step 2.
4. Singular charts keep the fully exact `PrepareChart`/`ChartIndicial`
   path: tags must stay exact there.

Risk: medium.  Float coefficient evaluation near a *barely-missed* factor
root could lose digits; the radius geometry bounds `|den(x0)|` away from 0
(nearest root is at distance ≥ radius), and the existing `E3` assert on
`d0` stays.  The decision logic stays exact, so no honesty loss.

### Lever 3 — eps-first cascade: QUANTIFIED, mostly NOT worth adopting
The question "solve the coupled eps-window W vs W separate per-order
systems with back-substitution":
- Coupled (current): every scalar op carries a width-W frame; matrix ops
  cost `nsp` sparse eps-slots each ⇒ per delivered eps order the overhead
  factor is `W·nsp / ((K+1)·L)` where L = eps-degree of the cleared system
  (≈ nsp).  With W = K+7: **(K+7)/(K+1) ≈ 1.2–2.0** — small.
- The old engine's *real* eps win was not the cascade itself but the
  **basis cache across eps orders** (L42): homogeneous work paid once, not
  ×(K+1).  Lever 1 eliminates the homogeneous basis on regular charts
  altogether, which subsumes this advantage.  A second per-order solver
  would violate the ONE-solver design (Config `$droppedReasons`,
  RewritePlan I2) for ≤2x — rejected.
- Residual micro-lever worth taking: in `blockSolveTPFrame`, when
  `deltaList` has a single eps^0 entry (every CASE-T step at a regular
  chart), `frInv` + two `frConv`s reduce to one scalar division — replace
  the O(W²) frame ops with an O(W) scale.  Expected ×1.1–1.3 on the
  recursion core; zero honesty impact (no eps-division either way).

### Lever 4 — data-plane slimming (expected ×1.2–1.5 global; ×1.5–2 on the tile phase)
- `ValidateLocalSolution` runs full `Position[...]` scans over the whole
  (W × nmax × d) coefficient array on **every** `EvaluateLocalSolution` /
  `MultiplyRational` / `CombineLocalSolutions` call (`SectorSeries.m:44–68`)
  — dozens of times per chart/tile.  Validate at construction through one
  `mkLocalSolution` choke point (the EpsSeries pattern, `EpsSeries.m:49`),
  keep an O(1) shape check on read.  Structure validation on every read is
  *not* one of the four core invariants; this is amortization, not
  relaxation.
- Tile integration (`IntegrateLocalSolution`): the per-(n, component)
  `cESc` EpsSeries construction + `esTimes` pair per monomial cell
  (`Integrate.m:162–177`) should run as one slab convolution per sector
  (the same transformation MultiplyRational already received in M5d).
- `TransportLine` keeps every chart's full LocalSolution; fine for memory,
  but `LineIntegral` re-sorts and re-validates them per master — hoist.

### Lever 5 — explicitly NOT copied from the old engine (invariant guard)
- Numeric eigen-snapping & rank tests (L44 sqrt-noise failure mode) — the
  exact indicial algebra stays.
- Silent strategy fallbacks and the singular-endpoint sector *fitting* —
  sector tags stay exact data.
- `Floor`-ed fractional eps powers, Mobius charts (DEC-18) — stay dropped.
- The old block-triangular `IntegrationSequence` is **subsumed** by lever 1
  for regular charts (a value solve is already matrix·vector, d²; block
  decomposition would only have shaved the d-column basis).  It could still
  shrink the d-column basis work at *singular* charts, but those are ~10%
  of charts — revisit only if singular charts dominate post-L1 profiles.

## 4. Combined expectation

| phase | now | after L1 | after L1+L2+L3micro+L4 |
|---|---|---|---|
| regular chart | 2–4 s | 0.4–0.8 s | **0.15–0.4 s** |
| line (≈20 charts) | 50–100 s | 10–20 s | 4–10 s |
| banana transport (8 lines) | ~8–14 min | ~1.5–3 min | ~0.5–1.5 min |
| tile/limit phase | ~10–25 min | unchanged | ~3–8 min (L4 tile item) |
| **end-to-end** | 25–40 min | — | **~4–9 min** |

The remaining ~10x to the old 30 s is per-op machinery: WL-interpreter
Association/frame-list arithmetic vs `SeriesData` C kernels, and the old
engine's pure-bignum steady state.  Once L1/L2 land, the natural follow-up
is packing the (by then fully numeric) frame lists of a chart into packed
arrays for the recursion inner loop — orthogonal to all invariants.

## 5. Prototype status and validation plan (not run here: no-kernel constraint)

Patch (this commit): `DiffExp2/Solve.m` (+`SolveValueRegular`, `capWindow`),
`DiffExp2/Transport.m` (`TransportLine` value-mode branch).  Inert unless
`DE2_VALUE_TRANSPORT=1`.

1. Baseline (flag off): `Tests/test_solve.m`, `Tests/test_transport.m`,
   `Tests/test_pins.m` — must be unchanged (the branch is dead code).
2. Flag on: same tests + `FT_EXAMPLES=bubble,sunrise`
   `Scripts/run_ft_stepwise2.m` — STEPWISE rows must match the flag-off run
   to MatchTol; then banana for the timing claim (expect ×4–6 transport).
3. `DEBUG_CHART=1` before/after per-chart timings; confirm the value-mode
   fallback never triggers on the banana lines except adjacent to the
   singular-target charts.
4. Promotion to default = a follow-up commit that flips the gate to a
   Config key after (1)–(3) hold.
