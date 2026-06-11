# Module spec: DiffExp2/EpsSeries.m

Status: M0 deliverable for Docs/RewritePlan.md section 3.2 (EpsSeries.m, ~300
lines).  This document is a contract: an implementation agent reading ONLY
this spec, Docs/RewritePlan.md, and the old code cited here must be able to
write the module and its unit tests.  Package context: `DiffExp2`EpsSeries``,
file `DiffExp2/EpsSeries.m`.

All old-code citations are against the frozen oracle tree (DiffExp/,
FeynmanTrick/) at the commit pinned by RewritePlan section 1 (f48cd94 plus
the post-baseline commits present in this working tree; line numbers were
verified against the current working tree on 2026-06-11).

---

## 1. PURPOSE

EpsSeries.m is the arithmetic foundation of DiffExp2: a single, honest
implementation of truncated eps-Laurent arrays (a coefficient list plus an
explicit completeness window `<|"Min" -> kmin, "CompleteMax" -> kmax|>`),
replacing the two divergent, silently-padding Laurent association layers that
the old code grew independently in FeynmanTrick/DiffExpIntegration.m and
DiffExp/RegularizedIntegration.m.  Every operation (add, scale, eps-power
shift, multiply, invert, divide, truncate, trim, tolerance-equality, and
expression conversion at the API boundary only) computes the exact window of
eps-orders on which its result is COMPLETE, and it is impossible to read,
pad, or fabricate a coefficient outside that window.  The core mathematical
subtlety this module owns is LAURENT DIVISION with explicit window-shift
semantics: dividing by a series whose eps^0 part vanishes (e.g. `b*eps`, or
the pseudo-resonance denominator `(b_i - b_j)*eps` of RewritePlan I2) shifts
both `Min` AND `CompleteMax` down — this is how the `1/(b eps)` integration
enhancement enters the pipeline, and it is a well-defined SHIFT, not an
inversion failure (math-lens review finding 6(iv), Docs/reviews/
rewrite_plan_review_3lens.json).  Coefficients are exact rationals/algebraics
or arbitrary-precision numbers; the module never forces `N`, `SetPrecision`,
or `Chop` on them.

---

## 2. PUBLIC SYMBOLS

All exported from `DiffExp2`EpsSeries``.  "EpsSeries" below means the
association defined in section 3.  "coefficient" means any Mathematica
expression FREE of the eps symbol (exact number, arbitrary-precision number,
or symbolic expression such as a rational function of a chart variable, a
`Logx`/theta carrier, or a symbolic indeterminate).

Constructors and validation:

- `ESNew[kmin_Integer, coeffs_List] -> EpsSeries`
  Builds the series Σ_{k=kmin}^{kmin+Length[coeffs]-1} coeffs[[k-kmin+1]] eps^k
  with `CompleteMax = kmin + Length[coeffs] - 1`.  `coeffs` must be a
  non-empty flat list of coefficients.  Validates structure (ERR-BAD-CONSTRUCT)
  including eps-freeness of every coefficient against the pinned eps symbol(s)
  (section 4, I-3).  By constructing through `ESNew` the caller ASSERTS the
  below-Min zero guarantee (section 3): all true coefficients below `kmin`
  are zero.

- `ESZero[kmax_Integer] -> EpsSeries`
  The zero series known complete through `kmax`:
  `<|"EpsWindow" -> <|"Min" -> kmax, "CompleteMax" -> kmax|>, "Coeffs" -> {0}|>`.
  This is the canonical all-zero representation (see ESTrim).  Replaces
  `LaurentZero` (FeynmanTrick/DiffExpIntegration.m:1154-1157,
  DiffExp/RegularizedIntegration.m:187-190) WITHOUT the old layer's
  freedom to make windows of zeros that read as data.

- `ESQ[expr_] -> True | False`
  Structural validator: association with exactly the keys of section 3,
  integer window fields, list length equal to window width >= 1.  Cheap; used
  by every public function's argument check.

Accessors:

- `ESWindow[s_?ESQ] -> <|"Min" -> _Integer, "CompleteMax" -> _Integer|>`
  Returns the EpsWindow association verbatim in the RewritePlan 3.1 shape, so
  it can be stored directly into a `LocalSolution["EpsWindow"]`.

- `ESMinPower[s_?ESQ] -> Integer` — `ESWindow[s]["Min"]`.

- `ESCompleteMax[s_?ESQ] -> Integer` — `ESWindow[s]["CompleteMax"]`.

- `ESCoefficient[s_?ESQ, k_Integer] -> coefficient`
  For `Min <= k <= CompleteMax`: the stored coefficient.  For
  `k < Min`: exact `0` — this is the below-Min zero GUARANTEE, not a
  fallback.  For `k > CompleteMax`: LOUD ERROR (ERR-WINDOW-READ).  There is
  deliberately NO default-value variant of this accessor (see F13).

- `ESCoefficientList[s_?ESQ, k1_Integer, k2_Integer] -> List`
  Coefficients of eps^k1 .. eps^k2 in order (length k2-k1+1), with below-Min
  entries filled by exact 0.  LOUD ERROR if `k2 > CompleteMax` (ERR-RANGE),
  if `k1 > k2` (ERR-RANGE), or if `k1 > Min` and any coefficient in
  `[Min, k1-1]` is non-negligible per `ESCoeffZeroQ` at the series' own scale
  (ERR-DROP-BELOW): slicing away a known nonzero pole must be a deliberate
  two-step (inspect `ESLeading`, then `ESShift`/handle the pole), never an
  implicit projection.  Replaces `LaurentToRange`/`LaurentToNonNegativeList`
  (FeynmanTrick/DiffExpIntegration.m:1192-1196,
  DiffExp/RegularizedIntegration.m:222-223), which silently zero-padded above
  the window and silently dropped negative powers.

- `ESLeading[s_?ESQ] -> {k_Integer, coefficient} | None`
  The lowest order whose coefficient is non-negligible (exact-nonzero, or
  numeric and >= laurentLeadTol relative to the series scale; symbolic
  not-provably-zero counts as nonzero).  `None` when every coefficient is
  negligible.  Returning `None` is honest data for the caller's decision, not
  a fallback; callers that require a nonzero lead must error themselves.

Arithmetic (all window rules normative; see section 3 table):

- `ESAdd[a_?ESQ, b_?ESQ, ...] -> EpsSeries`  (n-ary, folds binary)
  Sum.  `Min = Min[aMin, bMin]`, `CompleteMax = Min[aCM, bCM]`.  The result
  window is never the union: extending the shorter operand with zeros (old
  `LaurentAdd`, FeynmanTrick/DiffExpIntegration.m:1159-1170 takes
  `maxPower = Max[...]`; DiffExp/RegularizedIntegration.m:197-208 same) is
  FORBIDDEN (F2).  Stored coefficients above the result CompleteMax are
  dropped.  Coefficient normalization: numeric coefficients are combined
  with plain `Plus`; non-numeric ones with `Together@*Expand` (the old
  layer's anti-swell normalization, FeynmanTrick/DiffExpIntegration.m:1166).
  No `Chop`, no `N`.

- `ESScale[c_, s_?ESQ] -> EpsSeries`
  Multiply every coefficient by the eps-free scalar `c` (validated eps-free,
  ERR-BAD-CONSTRUCT).  Window unchanged — including when `c` is exactly 0
  (the result is a known-zero series on the same window).  Replaces
  `LaurentScale` (FeynmanTrick/DiffExpIntegration.m:1172-1175;
  DiffExp/RegularizedIntegration.m:192-195 — note the latter ran `PChop` over
  every coefficient; the new module does NOT chop).

- `ESShift[s_?ESQ, j_Integer] -> EpsSeries`
  Multiply by eps^j exactly: `Min += j`, `CompleteMax += j`, coefficients
  unchanged.  This is the series-level primitive behind the FT eps-prefactor
  bookkeeping (the `shift` argument of `MultiplyLaurentShifted`,
  FeynmanTrick/DiffExpIntegration.m:1241-1256, and the `J = eps^{k_j} I`
  prefactor convention documented at DiffExpIntegration.m:669-679; matrix
  side: FeynmanTrick/EpsPrefactors.m:194-209).  ESShift is the ONLY exported
  operation that raises CompleteMax, and it raises Min rigidly with it; it
  cannot pad.

- `ESTimes[a_?ESQ, b_?ESQ] -> EpsSeries`
  Cauchy product (convolution).  `Min = aMin + bMin`;
  `CompleteMax = Min[aCM + bMin, bCM + aMin]`.
  Derivation (normative): coefficient `c_k = Σ_i a_i b_{k-i}`; `c_k` is
  complete iff every term with a possibly-nonzero factor is known, i.e.
  `k - bMin <= aCM` (a-side) and `k - aMin <= bCM` (b-side).  Note the rule
  rewards trimming: raising an operand's `Min` past exact-zero leads widens
  the product window legitimately, because the below-Min zero guarantee is
  what makes those convolution terms drop.  Cost is O(width^2) on arrays of
  ~5-15 entries; not a performance concern (R6 concerns t-series, not eps).
  Replaces the convolution core of `MultiplyLaurentShifted`
  (FeynmanTrick/DiffExpIntegration.m:1241-1256) WITHOUT its caller-chosen
  `minOut`/`maxOut` output window (F3): the result window is computed, never
  requested.

- `ESInvert[d_?ESQ] -> EpsSeries`
  Multiplicative inverse `1/d`.  Let `L` = `ESLeading[d][[1]]` (LOUD ERROR
  ERR-DIV-ZERO when `ESLeading[d] === None`).  Then with
  `d = d_L eps^L (1 + u)`, `u` known on relative orders `1 .. dCM - L`:
  `Min = -L`, `CompleteMax = dCM - 2L`.
  Coefficients by the standard recursion `e_0 = 1/d_L`,
  `e_m = -(1/d_L) Σ_{j=1}^{m} d_{L+j} e_{m-j}` for relative orders
  `m = 1 .. (CompleteMax - Min) = dCM - L` (every `d_{L+j}` needed lies
  inside d's window, since `L + j <= L + dCM - L = dCM`).  Division never
  numericizes `1/d_L` (exact stays exact; `(b_i - b_j)` as a symbol stays
  symbolic).
  WINDOW-SHIFT SEMANTICS (the load-bearing case, math review finding 6(iv)):
  if `d` has vanishing eps^0 part — e.g. `d = b*eps` known complete through
  `dCM = K` — then `L = 1` and the result has `Min = -1`,
  `CompleteMax = K - 2`: BOTH edges shift down.  This is exact: the unknown
  tail of `d` at order `K+1` feeds `1/d` first at order `K-1` relative to the
  `eps^-1` lead.  It is NOT an error and NOT a special case in the code — the
  general window rule produces it.

- `ESDivide[a_?ESQ, b_?ESQ] -> EpsSeries`
  `ESTimes[a, ESInvert[b]]`.  Resulting window (derived, normative):
  `Min = aMin - L`, `CompleteMax = Min[aCM - L, aMin + bCM - 2L]` with `L`
  the leading index of `b`.  Two named consequences:
  (i) the `1/(b eps)` INTEGRATION ENHANCEMENT: `a` complete on `[0, K]`
  divided by `b1*eps` complete on `[1, K+1]` gives window `[-1, K-1]` — one
  full order of Laurent depth gained below, one order of completeness lost on
  top.  The lost top order must be paid for by the caller's order budget
  (RewritePlan 3.4 `enhancement(Lk)` term); EpsSeries makes the loss visible
  instead of laundering it (F11).
  (ii) the PSEUDO-RESONANCE shift of RewritePlan I2 / math review finding 3:
  the Frobenius denominator `(b_i - b_j)*eps` at an integer-spaced-a
  collision is exactly this division; it must run without complaint and
  produce the shifted window.  Solve.m's joint-solve construction then
  recombines to keep coefficients finite at eps=0; EpsSeries itself does no
  recombination.

- `ESTruncate[s_?ESQ, newCM_Integer] -> EpsSeries`
  Lower `CompleteMax` to `newCM`, dropping stored coefficients above.  LOUD
  ERROR if `newCM > CompleteMax` (ERR-TRUNCATE-EXTEND — there is no padding
  operation in this module) or if `newCM < Min` (use the canonical zero/
  rethink the call; also ERR-TRUNCATE-EXTEND, message distinguishes the two).

- `ESTrim[s_?ESQ] -> EpsSeries`
  Advance `Min` past leading coefficients that are negligible per
  `ESCoeffZeroQ` at the series' own scale (exact zeros always; numeric zeros
  at relative laurentLeadTol; symbolic not-provably-zero coefficients STOP
  the trim).  `CompleteMax` NEVER changes.  If everything is negligible the
  result is `ESZero[CompleteMax[s]]` — completeness is preserved, unlike old
  `LaurentTrim` which reset the window of an all-zero series to `[0, 0]`
  (FeynmanTrick/DiffExpIntegration.m:1186-1187;
  DiffExp/RegularizedIntegration.m:213) (F8).  After a numeric-zero advance
  the below-Min guarantee is "zero at laurentLeadTol", documented as such.
  ESTrim is never called implicitly by the arithmetic ops; window honesty
  must not depend on hidden tolerance decisions (the one exception is
  ESInvert/ESDivide, whose leading-index search uses the identical predicate
  and is specified above).

Predicates and comparison:

- `ESCoeffZeroQ[c_, scale_] -> True | False`
  THE zero test of the new core, exported for reuse (Solve.m/Transport.m
  matching solves use the same semantics — Lessons Ledger seed
  "numerical-zero leading-coefficient skipping (generalizes to matching
  solves)", RewritePlan section 5).  True iff `PossibleZeroQ[c]`, OR `c` is
  numeric (NumericQ) and `Abs[N[c, wp]] < laurentLeadTol * scale` where `wp`
  and `laurentLeadTol` come from Tolerances.m and `scale > 0`.  Symbolic
  non-numeric `c` that is not provably zero -> False, ALWAYS: this predicate
  never substitutes values for `Logx`/theta/indeterminates (the theta-aware,
  Logx-probed variant lives where those symbols have semantics — the old
  model is DiffExp/IntegrationStrategies/Recurrence.m:595-634; SectorSeries/
  Solve pre-probe and pass a numeric value here).  The test is RELATIVE by
  construction; there is no absolute-threshold code path (F9).  Scale
  convention used by ESTrim/ESInvert/ESLeading/ERR-DROP-BELOW:
  `scale = Max[Abs[N[c_k]]]` over the series' numeric coefficients (symbolic
  coefficients excluded from the scale); if that scale is 0, only exact
  zeros qualify.
  Calibration precedent the tolerance must respect: the box campaign's
  apparent-singularity fix used the relative threshold
  `Max[10^(-ChopPrecision/2), 10^-24]` and found the `10^-24` floor
  load-bearing (DiffExp/IntegrationStrategies/Recurrence.m:613-618;
  Docs/FeynmanTrickBoxFamilyStatus.md:116-119).  laurentLeadTol's derivation
  from WorkingPrecision is Tolerances.m's contract; this module only consumes
  the named value.

- `ESSameQ[a_?ESQ, b_?ESQ] -> True | False`
  Equality to tolerance on the SHARED complete range: for every
  `k in [Min[aMin, bMin], Min[aCM, bCM]]`, the difference
  `ESCoefficient[a,k] - ESCoefficient[b,k]` must satisfy `ESCoeffZeroQ` at
  matchTol (Tolerances.m) with scale = the max numeric coefficient magnitude
  of both series on that range (both all-zero -> True; scale 0 -> exact
  comparison).  Symbolic coefficient differences must be `PossibleZeroQ`.
  Window AGREEMENT is deliberately not required — two series may be equal as
  functions while known to different depths; consumers that need window
  equality check `ESWindow` themselves.  Never errors; returns False
  honestly.  If `Min[aCM,bCM] < Min[aMin,bMin]` the shared range collapses
  to the windows' guaranteed-zero region and the result is True only if both
  are all-zero there (i.e. trivially True); flagged in section 10.

Conversion (API boundary ONLY — see usage rule below):

- `ESFromExpression[expr_, epsSym_Symbol, kmax_Integer] -> EpsSeries`
  Exact Laurent expansion of a closed-form expression around `epsSym = 0`,
  complete through `kmax`.  Implementation: `Series[expr, {epsSym, 0, kmax}]`
  once, then exact `SeriesData` part extraction; coefficients passed through
  `Together` only, NEVER `N`/`SetPrecision` (Gamma-prefactor coefficients
  come out as exact `EulerGamma`/`Zeta` combinations — the FT gamma
  prefactor, FeynmanTrick/DiffExpIntegration.m:571,880, is the canonical
  client).  `Min` = the exact leading exponent reported by `SeriesData`.
  Eps-FREE input (where `Series` returns the expression itself, not a
  `SeriesData`) is legal and total: result window `[0, kmax]`, coefficients
  `{expr, 0, ..., 0}` (the constant-term analogue of old `SeriesAlways`,
  DiffExp/SeriesOps.m:168-171); in particular
  `ESFromExpression[0, epsSym, kmax] := ESZero[kmax]`.
  LOUD ERROR (ERR-EXPAND-FAIL) when: `Series` fails or returns a non-
  `SeriesData` result that still CONTAINS `epsSym` (the old fall-back-to-
  unexpanded branch); the leading power is not an integer
  (the old code FLOORED fractional powers,
  FeynmanTrick/DiffExpIntegration.m:1224 — forbidden, F7); any extracted
  coefficient still contains `epsSym`; or expr has an essential singularity.
  There is NO fall-back-to-unexpanded-expression branch (old
  `ExpandIBPCoeffLaurent` `Quiet[Check[..., $Failed]]` then "use expr as-is",
  FeynmanTrick/DiffExpIntegration.m:1207-1217 — forbidden, F6).

- `ESToExpression[s_?ESQ, epsSym_Symbol] -> expression`
  `Σ_k ESCoefficient[s,k] epsSym^k` over the stored window, returned as a
  plain expression.  Carries NO window metadata — therefore: USAGE RULE
  (normative for the whole library): ESToExpression/ESFromExpression may be
  called only at the API boundary (API.m input parsing and result output, and
  test harnesses).  No core module converts an EpsSeries to an expression
  mid-pipeline and re-expands it; that round trip is the D2 disease
  (RewritePlan section 2) in eps-direction miniature, and it launders
  completeness (the x-direction analogue — `Normal` + re-`Series` laundering
  one phantom top order — is catalogued at
  Docs/FeynmanTrickBoxFamilyStatus.md:149-150).

Mapping:

- `ESMap[f_, s_?ESQ] -> EpsSeries`
  Apply `f` to every stored coefficient; window unchanged; result
  re-validated by `ESNew` (so `f` introducing eps-dependence errors loudly).
  Intended uses: precision raising at the API boundary (the 2xWP input
  lesson, DiffExp/Transport.m:541-562 — the RAISING is API.m's job, ESMap is
  just the vehicle), component extraction when coefficients are vectors,
  explicit caller-owned chopping.  ESMap must not be used inside EpsSeries
  itself to apply numeric coercions (I-6).

Error context:

- `$ESErrorContext` (settable, default `"(no context)"`)
  A dynamically scoped string callers `Block` around EpsSeries calls, e.g.
  `Block[{DiffExp2`EpsSeries`$ESErrorContext = "chart x0=1, sector
  (a=-1,b=2,p=0), matching order k=3"}, ...]`.  Every EpsSeries error message
  embeds it (section 5), which is how this module — which knows nothing of
  charts or sectors — satisfies the library-wide rule that errors name
  (chart, sector, order).  SectorSeries/Solve/Transport/Integrate MUST set it
  around every EpsSeries call cluster; their specs carry the matching
  obligation.

Total: 22 exported symbols + 1 exported context variable.

---

## 3. DATA CONTRACTS

### 3.1 EpsSeries (produced/consumed by every symbol above)

```
EpsSeries = <|
  "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,   (* integers *)
  "Coeffs"    -> {c_kmin, c_{kmin+1}, ..., c_kmax}           (* length kmax-kmin+1 >= 1 *)
|>
```

The `"EpsWindow"` value is VERBATIM the EpsWindow association of RewritePlan
3.1 (`<|"Min" -> kmin, "CompleteMax" -> kmax|>`), so `ESWindow[s]` can be
stored into `LocalSolution["EpsWindow"]` without translation.

Semantics: the object represents
`f(eps) = Σ_{k=kmin}^{kmax} c_k eps^k + O(eps^{kmax+1})`
together with TWO guarantees:

- COMPLETENESS: every order `k <= CompleteMax` of the true function is
  exactly (up to coefficient arithmetic precision) the stored/derived value.
- BELOW-MIN ZERO GUARANTEE: the true function has no nonzero coefficient
  below `Min`.  Hence every EpsSeries is complete on ALL of
  `(-inf, CompleteMax]`; `Min` is simultaneously the storage offset and an
  exact zero certificate.  (After a numeric `ESTrim` advance the certificate
  is "zero at laurentLeadTol relative to the series scale"; exact inputs give
  exact certificates.)  This guarantee is what makes the convolution and
  addition window rules sound and makes trimming a legitimate window
  improvement.

Coefficient domain: any expression free of the eps symbol.  Exact integers,
rationals, algebraics and symbolic constants pass through ALL operations
structurally unharmed — the module contains no call to `N`, `SetPrecision`,
or `Chop` on stored coefficients (numericization happens only inside
tolerance PREDICATES, on local copies).  Arbitrary-precision numbers keep
their precision under Mathematica's own arithmetic rules; the module neither
raises nor floors precision (that is the API.m / Transport.m 2xWP +
$MinPrecision ledger lesson, out of scope here).

Canonical zero: `ESZero[kmax]`.  All-zero series produced by arithmetic keep
their honest window; `ESTrim` canonicalizes them to `ESZero[CompleteMax]`.

Negative `Min` is ordinary (Laurent!): in particular the homogeneous
true-resonance sectors of RewritePlan 3.1 have `kmin = -p` BY CONSTRUCTION
(invariant exemption; math review finding 14).  No EpsSeries code path may
assume `Min >= 0` (the old `LaurentToNonNegativeList` projection is exactly
the assumption to kill — F12).

### 3.2 Relation to Sector["Coeffs"]

RewritePlan 3.1 defines (verbatim):

```
Sector = <|
  "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
  "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT *)
|>
```

The eps direction (`k`) of that 3D array carries one EpsWindow shared by the
sector; SectorSeries.m slices per `(n, comp)` into EpsSeries objects (or
operates on the shared window using the rules in the table below — the RULES
are normative, the slicing is an implementation choice for SectorSeries).
EpsSeries.m is the single owner of the window algebra; no other module may
re-derive its own variant (the disease being replaced is precisely the two
divergent copies at FeynmanTrick/DiffExpIntegration.m:1143-1256 and
DiffExp/RegularizedIntegration.m:176-223).

### 3.3 Window arithmetic table (normative)

For inputs `a` with window `[aMin, aCM]`, `b`/`d` with `[bMin, bCM]`, leading
index `L = ESLeading[d][[1]]`:

| op                | result Min          | result CompleteMax                  |
|-------------------|---------------------|-------------------------------------|
| ESAdd[a,b]        | Min[aMin, bMin]     | Min[aCM, bCM]                       |
| ESScale[c,a]      | aMin                | aCM                                 |
| ESShift[a,j]      | aMin + j            | aCM + j                             |
| ESTimes[a,b]      | aMin + bMin         | Min[aCM + bMin, bCM + aMin]         |
| ESInvert[d]       | -L                  | dCM - 2L                            |
| ESDivide[a,d]     | aMin - L            | Min[aCM - L, aMin + dCM - 2L]       |
| ESTruncate[a,m]   | aMin                | m  (m <= aCM enforced)              |
| ESTrim[a]         | >= aMin (advanced)  | aCM (never changes)                 |

Closure: every rule yields `Min <= CompleteMax` (proofs one-line from
`aMin <= aCM`, `bMin <= bCM`, `L <= dCM`); asserted anyway (I-2).

Worked shift example (must appear as a code comment in the implementation):
`ESDivide` of `1 + O(eps^6)` (window `[0,5]`, coeffs `{1,0,0,0,0,0}`) by
`3 eps + O(eps^7)` (window `[1,6]`) gives `L = 1`, window
`[-1, Min[5-1, 0+6-2]] = [-1, 4]`, coefficients `{1/3, 0, 0, 0, 0, 0}` —
exact `1/3`, both edges shifted down by one.

---

## 4. INVARIANTS (always on, cheap)

- I-1 STRUCTURE: every EpsSeries returned by a public function satisfies
  `ESQ`: exactly the two keys, `EpsWindow` with exactly `Min`/`CompleteMax`
  integer entries, `Coeffs` a list of length `CompleteMax - Min + 1 >= 1`.
  Checked at construction (all public ops construct through one internal
  `mkSeries[kmin, kmax, coeffs]` that performs the check — single choke
  point).

- I-2 NON-EMPTY WINDOW: `Min <= CompleteMax` after every operation
  (mathematically guaranteed by the section 3.3 rules; asserted at the choke
  point; violation = internal bug, ERR-WINDOW-EMPTY).

- I-3 EPS-FREE COEFFICIENTS: at `ESNew`/`ESMap`/`ESScale`/`ESFromExpression`,
  every (new) coefficient is `FreeQ` of the pinned eps symbol(s).  The
  symbol set comes from Config.m's validated accessor (both `eps` and
  `\[Epsilon]` are accepted library-wide per RewritePlan 3.2 API.m); if
  Config is not yet initialized the check uses the literal global symbols
  `Global`eps` and `\[Epsilon]`.  Catches the D1 type-poisoning class where
  an un-substituted eps rides inside a "coefficient".

- I-4 NO PADDING: no exported function increases `CompleteMax` relative to
  the honest rule; `ESShift` is the only CompleteMax-raising operation and it
  moves `Min` rigidly with it.  Enforced by construction (there is no pad/
  extend/resize entry point) and by the unit suite (section 8, tests 4, 16).

- I-5 BELOW-MIN GUARANTEE PRESERVED: each op's rule keeps the zero
  certificate valid (add: min of Mins; times: sum of Mins; invert: lead is
  `1/d_L != 0` at `-L`; trim: only advances past negligible coefficients).
  Spot assert: after `ESTrim`, the first stored coefficient is
  non-negligible or the series is canonical zero.

- I-6 NO NUMERIC COERCION: the module never applies `N`, `SetPrecision`, or
  `Chop` to stored coefficients (predicates numericize LOCAL copies only).
  Enforced by code review + unit tests 13/14 (exact-rational passthrough,
  precision preservation); there is no runtime assert because exactness of
  every coefficient is the test's job, not a per-op scan.

These are all O(window width) or O(1); they stay on in production
(RewritePlan A2).

---

## 5. ERROR CONTRACT

Mechanism: one internal helper `esError[id_String, details__]` that prints
`"DiffExp2 EpsSeries error <id>: <details>; window=<Min,CompleteMax if
applicable>; context: " <> $ESErrorContext` and then `Abort[]`s — same
loudness class as old `DiffExp`Utilities`ReportError`
(DiffExp/Utilities.m:66).  If the orchestrator later standardizes a
library-wide error helper (open question O-1), only `esError`'s body changes.
Every message embeds `$ESErrorContext`, which callers must populate with
chart/sector/order names (section 2); EpsSeries supplies the local facts
(requested order, window, magnitudes, tolerance values) and the caller
context supplies the physics location.  NO EpsSeries error is maskable by
options; there is no `Quiet` mode.

Enumerated loud errors:

- ERR-BAD-CONSTRUCT — fires in `ESNew`/`ESScale`/`ESMap`/`ESZero` argument
  validation: non-integer window field, empty or non-list `Coeffs`, length
  mismatch, eps-dependent coefficient.  Message carries: which check failed,
  the offending index and coefficient (`InputForm`, truncated to ~200 chars),
  the eps symbol found (for the eps-free check), and context.

- ERR-WINDOW-READ — fires in `ESCoefficient[s, k]` with `k > CompleteMax`.
  Message carries: requested `k`, the window `[Min, CompleteMax]`, and
  context.  This is the single most important error in the module: it is the
  EpsWindow-propagation backstop of RewritePlan 3.4 ("EpsWindow propagation
  turns any miscount into a named error").

- ERR-RANGE — fires in `ESCoefficientList[s, k1, k2]` when `k1 > k2` or
  `k2 > CompleteMax`.  Carries: `k1`, `k2`, window, context.

- ERR-DROP-BELOW — fires in `ESCoefficientList[s, k1, k2]` when `k1 > Min`
  and some coefficient in `[Min, k1-1]` fails `ESCoeffZeroQ` at the series
  scale.  Carries: the orders being dropped, their numeric magnitudes (or
  `Symbolic`), the scale and laurentLeadTol used, window, context.  (Makes
  the old silent pole-projection of `LaurentToNonNegativeList` impossible —
  F12.)

- ERR-DIV-ZERO — fires in `ESInvert`/`ESDivide` when `ESLeading[denominator]
  === None`.  Carries: denominator window, max numeric coefficient
  magnitude, scale, laurentLeadTol, and context.  NOTE the asymmetry this
  error must respect: a denominator with vanishing eps^0 part but a
  non-negligible HIGHER coefficient (e.g. `b*eps`) is NOT this error — it is
  the window shift of section 2/3.3.  Only an entirely-negligible
  denominator fires.

- ERR-TRUNCATE-EXTEND — fires in `ESTruncate[s, m]` when `m > CompleteMax`
  (a padding attempt) or `m < Min` (an empty-window request).  Carries: `m`,
  window, which of the two violations, context.

- ERR-EXPAND-FAIL — fires in `ESFromExpression` per section 2 (Series
  failure / non-integer leading power / residual eps in coefficients /
  essential singularity).  Carries: the input expression (`InputForm`,
  truncated), `epsSym`, requested `kmax`, what `Series` returned (head, and
  leading power when fractional), context.

- ERR-WINDOW-EMPTY — internal choke-point assert violation (I-2); carries
  the op name and both input windows.  Unreachable by design; still loud.

### NO SILENT FALLBACKS — enumerated temptations, each forbidden

Each item names the old-code site that did it, and the rule that replaces it:

- F1 Out-of-window read returns 0.
  Old: `LaurentCoeff` (FeynmanTrick/DiffExpIntegration.m:1146-1152;
  DiffExp/RegularizedIntegration.m:179-185) returned 0 BOTH below MinPower
  and above MaxPower.  New: below `Min` is the zero guarantee (legitimate);
  above `CompleteMax` is ERR-WINDOW-READ.  No accessor with a default value
  exists.

- F2 Union-window addition (zero-padding the shorter operand).
  Old: `LaurentAdd` takes `maxPower = Max[...]`
  (FeynmanTrick/DiffExpIntegration.m:1159-1170;
  DiffExp/RegularizedIntegration.m:197-208).  New: `CompleteMax = Min[...]`,
  always.

- F3 Caller-chosen output windows in multiplication.
  Old: `MultiplyLaurentShifted[..., minOut, maxOut]`
  (FeynmanTrick/DiffExpIntegration.m:1241-1256) computed whatever range the
  caller asked for, complete or not.  New: `ESTimes` has no window
  parameters; the window is computed by the section 3.3 rule.

- F4 Zero-padding boundary vectors to a level-wide max after a uniform
  shift.  Old: `ShiftRawBoundariesToFinite`
  (FeynmanTrick/DiffExpIntegration.m:1258-1280, zero-reads at 1265-1271),
  flagged by RewritePlan A1 and the campaign defect list
  (Docs/FeynmanTrickBoxFamilyStatus.md:147-148).  New: per-master EpsSeries
  with individual honest windows; the uniform per-level shift is `ESShift`;
  the FT-side replacement carries `CompleteMax` metadata instead of pads
  (FT hardening contract, RewritePlan A1 "ShiftRawBoundaries zero-padding ->
  CompleteMaxPower metadata").

- F5 "Keep at least one output order" floor masking an empty window.
  Old: FeynmanTrick/DiffExpIntegration.m:658-661.  New: an empty resulting
  window in a CONSUMER's order computation must surface as that consumer's
  loud error; EpsSeries cannot represent an empty window (I-2) and provides
  no floor.

- F6 Conversion falling back to the unexpanded expression when `Series`
  fails.  Old: `ExpandIBPCoeffLaurent`
  (FeynmanTrick/DiffExpIntegration.m:1207-1217).  New: ERR-EXPAND-FAIL.

- F7 `Floor[]` of a non-integer leading power.
  Old: FeynmanTrick/DiffExpIntegration.m:1224 (`If[!IntegerQ[minPower],
  minPower = Floor[minPower]]`).  New: ERR-EXPAND-FAIL naming the fractional
  power.  (Fractional eps-powers do not exist in this library's domain; if
  one appears the input is wrong and must be said so.)

- F8 Trimming an all-zero series to window `[0, 0]`.
  Old: FeynmanTrick/DiffExpIntegration.m:1186-1187;
  DiffExp/RegularizedIntegration.m:213.  New: `ESTrim` returns
  `ESZero[CompleteMax]` — completeness survives.

- F9 Absolute zero thresholds.
  Old: `zeroCoeffQ`'s `10^-40` (FeynmanTrick/DiffExpIntegration.m:1177-1178),
  `RationalizationTolerance` fallbacks `10^-40`
  (FeynmanTrick/DiffExpIntegration.m:154-155, 182-183), `EffectiveZeroExprQ`
  absolute test (DiffExp/RegularizedIntegration.m:286-308).  New: ALL zero
  decisions go through `ESCoeffZeroQ` (relative, scale explicit,
  laurentLeadTol from Tolerances.m).  This is RewritePlan 3.2's "the FT
  LaurentTrim 10^-40 ABSOLUTE test becomes a relative one" realized.

- F10 Forcing exact data numeric.
  Old: `ApplyRegularizationStep` `SetPrecision`s `a`, `b`, `c` and the
  combination coefficients (DiffExp/RegularizedIntegration.m:373-380).  New:
  no numeric coercion anywhere in EpsSeries (I-6); exact rationals/algebraic
  `b`'s flow through division exactly.

- F11 Phantom top orders when dividing by an eps^0-vanishing denominator.
  Old: `ApplyRegularizationStep` extends `transformedMaxOrder` by one and
  reads the missing source order as 0
  (DiffExp/RegularizedIntegration.m:381, zero-reads at 388, 391), so after
  the `1/(b eps)` shift (411-416) the top output order silently contains
  only partial information.  New: the division window rule SHIFTS
  `CompleteMax` down; missing orders are the caller's budget problem
  (RewritePlan 3.4 enhancement term), surfaced by ERR-WINDOW-READ when
  under-budgeted, never fabricated.

- F12 Silent projection to non-negative eps powers at output.
  Old: `LaurentToNonNegativeList`
  (FeynmanTrick/DiffExpIntegration.m:1195-1196;
  DiffExp/RegularizedIntegration.m:222-223) dropped negative powers — i.e.
  dropped POLES — without a check.  New: ERR-DROP-BELOW; deliberate pole
  dropping requires the caller to inspect `ESLeading` and truncate/shift
  explicitly (the FT b!=0 endpoint drop RULE remains, but it lives in
  Integrate.m/FT as an explicit documented step, RewritePlan 3.3).

- F13 Default-to-zero evaluation chains.
  Old consumers wrap coefficient evaluation in `Quiet[Check[..., 0]]`
  (FeynmanTrick/DiffExpIntegration.m:1037-1059, and the
  Limit -> substitution -> 0 chain at 1080-1093, plus the
  `If[NumericQ[...], ... += ...]` silent drop at 1061-1063).  Those sites are
  FT-hardening work (RewritePlan A1), but EpsSeries must not make the
  pattern expressible: no accessor takes a default, no constructor accepts
  `$Failed`-like placeholders (ERR-BAD-CONSTRUCT catches them as
  eps-free-but-invalid only if non-expression heads appear — the validator
  explicitly rejects `$Failed`, `Indeterminate`, `ComplexInfinity`,
  `DirectedInfinity[_]`, `Overflow[]` coefficients, same pattern class as
  DiffExp/Utilities.m:129-130).

---

## 6. ABSORBED OLD CODE

EpsSeries.m replaces, with file:line:

1. The FT-layer Laurent association toolkit,
   FeynmanTrick/DiffExpIntegration.m:1143-1256:
   `LaurentMaxPower` (1143-1144) -> `ESCompleteMax`;
   `LaurentCoeff` (1146-1152) -> `ESCoefficient` (loud above window);
   `LaurentZero` (1154-1157) -> `ESZero`;
   `LaurentAdd` (1159-1170) -> `ESAdd` (min-window);
   `LaurentScale` (1172-1175) -> `ESScale` (no chop);
   `zeroCoeffQ` (1177-1178) -> `ESCoeffZeroQ` (relative);
   `LaurentTrim` (1180-1190) -> `ESTrim` (window-preserving);
   `LaurentToRange`/`LaurentToNonNegativeList` (1192-1196) ->
   `ESCoefficientList` (loud);
   `ExpandIBPCoeffLaurent` (1198-1239) -> `ESFromExpression` (loud; the
   dimension-variable substitution at 1200, 1209, 1216 stays in the FT layer
   — EpsSeries receives the already-substituted expression);
   `MultiplyLaurentShifted` (1241-1256) -> `ESShift` + `ESTimes`.

2. The boundary-shift bookkeeping `ShiftRawBoundariesToFinite`,
   FeynmanTrick/DiffExpIntegration.m:1258-1280 -> `ESTrim` + `ESShift` +
   per-master windows (consumer-side restructuring is FT hardening, but the
   primitives and the no-padding rule live here).

3. The DUPLICATE core-layer Laurent toolkit,
   DiffExp/RegularizedIntegration.m:176-223 (`LaurentMaxPower` 176-177,
   `LaurentCoeff` 179-185, `LaurentZero` 187-190, `LaurentScale` 192-195,
   `LaurentAdd` 197-208, `LaurentTrim` 210-220, `LaurentToNonNegativeList`
   222-223) and its zero predicate `EffectiveZeroExprQ` (286-308; the
   theta-branch probing half of that predicate, 291-305, moves to
   SectorSeries/Solve where theta has semantics — EpsSeries gets only the
   plain relative test).

4. The `{epsMinPower, gList}` offset-list convention of regularized
   integration, DiffExp/RegularizedIntegration.m:347-352 (format comment),
   and the eps-window halves of `ApplyRegularizationStep` (363-436): the
   `1/(b eps)` special case (411-416) and the geometric
   `1/(1 + a + b eps)` expansion (418-432) both become `ESDivide`; the
   per-order window bookkeeping of `IntegrateSingularTermLaurent`
   (489-495 epsMin/epsMax setup; 681-702 result-window assembly) becomes
   EpsSeries windows.  The INTEGRATION content of those functions (case
   analysis, boundary terms) moves to Integrate.m per its own spec; only the
   eps-Laurent arithmetic is absorbed here.

5. The honest-window embryo in the FT combine step,
   FeynmanTrick/DiffExpIntegration.m:635-667 (`completeMaxPower` trimming of
   incomplete top orders, comment 635-641) — the idea generalizes into the
   `CompleteMax` field; its floor (658-661) is deleted (F5); its
   trust warning (662-666) becomes ERR-WINDOW-READ at the point of use.

NOT absorbed (boundary statement): DiffExp/SeriesOps.m (218 lines) operates
on x-direction `SeriesData` and `Logx` towers — that is SectorSeries.m
territory; EpsSeries has no x/t concept.  FeynmanTrick/EpsPrefactors.m
(matrix-level prefactor inference, 33-185) stays in the FT layer; only the
series-level shift convention maps to `ESShift`.  The old core's per-eps-order
tower representation (eps orders as separate transported systems) is not
"absorbed" anywhere — it is the D2 disease the sector-native representation
deletes (RewritePlan section 2).

### Numerical lessons that MUST be preserved (read from the old code)

- N-1 RELATIVE leading-zero classification, with the campaign calibration.
  Misclassifying a leading Laurent coefficient shifts the entire series and,
  in the FT layer, the uniform per-level prefactor shift (legacy review
  finding 13a).  The absolute `10^-40` test
  (FeynmanTrick/DiffExpIntegration.m:1177-1178) failed exactly this way; the
  working fix in the box campaign was a RELATIVE threshold
  `Max[10^(-ChopPrecision/2), 10^-24]` judged against the series' own
  magnitude, with Logx/theta probing where those symbols occur
  (DiffExp/IntegrationStrategies/Recurrence.m:595-634;
  Docs/FeynmanTrickBoxFamilyStatus.md:116-124: "the 1e-24 floor is
  load-bearing").  `ESCoeffZeroQ` + Tolerances`laurentLeadTol must reproduce
  this behavior class on numeric coefficients; the probing variant is
  SectorSeries/Solve's obligation using this same predicate underneath.

- N-2 Cancellation residues sit far ABOVE chop scale.  At WP 300 the box L2
  apparent-singularity sources carried ~1e-29 RELATIVE residues
  (precision-tracked at ~271 digits, far above PChop and
  `10^(-ChopPrecision/2)`) — Docs/FeynmanTrickBoxFamilyStatus.md:104-108.
  Consequence for this module: the zero test must use the RELATIVE scale of
  the series at hand, and `laurentLeadTol` must be loose enough to absorb
  ~30 digits of cancellation at WP 300 (Tolerances.m owns the number; this
  spec pins the requirement and test 11 pins the behavior).

- N-3 Division by an eps^0-vanishing denominator costs one COMPLETE order
  per leading-index unit and gains one Laurent order below — never model it
  as "compute one more order from zeros" (the
  `ApplyRegularizationStep` extension, DiffExp/RegularizedIntegration.m:381,
  388, 391).  The cost must be paid upstream by the order budget
  (RewritePlan 3.4 `enhancement(Lk)`; empirical record: box_bubble needed
  9 orders at L2 / 11 at L1, Docs/FeynmanTrickBoxFamilyStatus.md:21-26, and
  the stepwise runner now carries full incoming depth because single-level
  cuts starved deep chains, Scripts/run_ft_stepwise.m:154-161).

- N-4 Together/Expand normalization on symbolic coefficients.  The old
  layers ran `Together[Expand[...]]` on every combined coefficient
  (FeynmanTrick/DiffExpIntegration.m:1166, 1252) — without it, convolution
  chains of x-rational IBP coefficients blow up structurally.  Keep for
  non-numeric coefficients; identity on plain numbers.

- N-5 Exactness is information.  Everywhere the old pipeline forced exact
  data numeric (`SetPrecision` in ApplyRegularizationStep,
  DiffExp/RegularizedIntegration.m:373-380; `numericAtActivePrecision`
  wrappers, FeynmanTrick/DiffExpIntegration.m:88-95) it later needed
  rationalization/snap tolerances to recover structure.  The sector-native
  design carries tags exactly (RewritePlan I1/3.1); EpsSeries is where that
  policy bites operationally: exact in -> exact out, no exceptions (I-6,
  test 13).

- N-6 Window truth beats output convenience.  The `LaurentToNonNegativeList`
  projection and the zero-padded boundary lists let callers consume
  fixed-shape lists, and every campaign defect in the "limit/direct paths
  lack the integrate path's incomplete-top-order trimming;
  ShiftRawBoundariesToFinite zero-pads unknown orders" cluster
  (Docs/FeynmanTrickBoxFamilyStatus.md:146-149) traces to shape convenience
  defeating completeness.  EpsSeries deliberately has NO fixed-shape output
  helper; consumers adapt to windows, not vice versa.

---

## 7. DEPENDENCIES

May call (acyclic order per task contract:
Tolerances < Config < EpsSeries < SectorSeries < Indicial < Solve <
Transport/Integrate < API):

- `DiffExp2`Tolerances``: `laurentLeadTol` (zero/lead classification),
  `matchTol` (ESSameQ), and the working-precision value used to numericize
  LOCAL copies inside predicates.  Exact symbol casing per
  Docs/specs/Tolerances.md; this spec uses the RewritePlan 3.2 names.
- `DiffExp2`Config``: ONLY the pinned eps symbol set for the I-3 eps-free
  validation (both `eps` and `\[Epsilon]`, RewritePlan 3.2 API.m note), with
  the literal-global fallback when Config is uninitialized.  `epsSym` is an
  EXPLICIT argument of the conversion functions, so Config is not needed for
  any arithmetic semantics.

Must NOT call: SectorSeries, Indicial, Solve, Transport, Integrate, API, any
`DiffExp`` (old) or `FeynmanTrick`` context.  EpsSeries knows nothing about
x/t, charts, sectors, Logx, theta symbols, or matrices.

Consumed by: SectorSeries.m (sector coefficient windows), Indicial.m
(eps-rational residue manipulation), Solve.m (Frobenius recursion at
symbolic eps — every recursion-denominator division is `ESDivide`),
Transport.m (eps-graded weight matching), Integrate.m (the
`(a+n+1+b eps)^-(j+1)` denominators), API.m (I/O conversion), and the
hardened FT layer through the API.

---

## 8. UNIT TESTS

File: `Tests/test_eps_series.m` (runs under the standard battery harness;
zero kernel-time dependence on the old library).  Where exact values are
stated they are closed-form expectations, compared with `===`/exact equality
unless noted.

1. `construct_validate_ok` — `ESNew[-2, {1, 0, 3/7}]` returns an ESQ-true
   object with `ESWindow -> <|"Min" -> -2, "CompleteMax" -> 0|>`,
   `ESCoefficient[s, -2] === 1`, `ESCoefficient[s, 0] === 3/7`.
2. `construct_validate_errors` — each of: `ESNew[1/2, {1}]`, `ESNew[0, {}]`,
   `ESNew[0, {eps}]` (eps symbol), `ESNew[0, {$Failed}]`,
   `ESNew[0, {Indeterminate}]` aborts with ERR-BAD-CONSTRUCT (harness runs
   each in a `CheckAbort`).
3. `coefficient_access_semantics` — for `s = ESNew[-1, {2, 5}]` (window
   `[-1,0]`): `ESCoefficient[s, -3] === 0` (guarantee, no error);
   `ESCoefficient[s, 0] === 5`; `ESCoefficient[s, 1]` aborts
   (ERR-WINDOW-READ) and the message contains "1", "-1", "0", and the
   blocked `$ESErrorContext` string.
4. `add_min_window_no_union` — `a = ESNew[0, {1,1,1,1,1,1}]` (CM 5),
   `b = ESNew[-1, {7, 2, 3, 4, 9}]` (CM 3):
   `ESAdd[a,b]` has window `[-1, 3]` (NOT 5) and coefficients
   `{7, 3, 4, 5, 10}` exactly.  Regression pin against old `LaurentAdd`
   max-union (FeynmanTrick/DiffExpIntegration.m:1162).
5. `mul_window_convolution` — `a = ESNew[0, {1,1,1}]` (=1+eps+eps^2, CM 2),
   `b = ESNew[1, {2,0,1}]` (=2eps+eps^3, CM 3): product window
   `[1, Min[2+1, 0+3]] = [1,3]`, coefficients `{2, 2, 3}` exactly.
6. `div_by_beps_shift` (THE math-review 6(iv) pin) —
   `one = ESFromExpression[1, eps, 5]` (window `[0,5]`),
   `d = ESFromExpression[3 eps, eps, 6]` (window `[1,6]`):
   `q = ESDivide[one, d]` has window `[-1, 4]` and
   `ESCoefficient[q,-1] === 1/3` (EXACT Rational, `Head === Rational`),
   all other window coefficients exactly 0.  Both `Min` and `CompleteMax`
   shifted down by 1 relative to the numerator.
7. `div_geometric_exact` — `ESDivide[ESFromExpression[1, eps, 4],
   ESFromExpression[1 - eps, eps, 4]]` = window `[0, 4]`, coefficients
   `{1,1,1,1,1}` exact integers.
8. `invert_beps_trailing_zeros` — `d = ESNew[1, {b0, 0, 0, 0}]` with
   `b0 = 5/3` (window `[1,4]`): `ESInvert[d]` window `[-1, 2]`
   (`dCM - 2L = 4 - 2`), coefficients `{3/5, 0, 0, 0}` exact.
9. `pseudo_resonance_symbolic_shift` — with symbols `b1, b2`:
   `den = ESNew[1, {b1 - b2, 0, 0}]` (window `[1,3]`),
   `num = ESFromExpression[1 + eps, eps, 2]`:
   `q = ESDivide[num, den]` runs WITHOUT error (division by an
   eps^0-vanishing denominator is a shift, not a failure), window
   `[-1, Min[2-1, 0+3-2]] = [-1, 1]`,
   `ESCoefficient[q, -1] === 1/(b1 - b2)` (`Together`-normalized symbolic).
10. `invert_zero_denominator_error` — `ESInvert[ESZero[5]]` aborts with
    ERR-DIV-ZERO; message contains the window `[5,5]` and tolerance info.
11. `relative_not_absolute_lead` (the LaurentTrim 10^-40 regression pin,
    FeynmanTrick/DiffExpIntegration.m:1177-1178) — at WP-300-class numbers:
    (a) `d = ESNew[0, {SetPrecision[10^-60, 300], 2, 1}]`: the lead of
    `ESTrim[d]`/`ESInvert[d]` is index 1 (the `2`), because `10^-60` is
    negligible RELATIVE to scale 2;
    (b) `d2 = ESNew[0, {SetPrecision[10^-60, 300],
    SetPrecision[10^-58, 300]}]`: the lead is index 0 — absolutely tiny but
    relatively O(10^-2), NOT skipped.  (An absolute 10^-40 test gets (b)
    wrong; this test kills it.)
12. `trim_preserves_zero_window` — `ESTrim[ESNew[-2, {0, 0, 0, 0}]]`
    (window `[-2,1]`) === `ESZero[1]`: `CompleteMax` 1 survives.  Regression
    pin against old `LaurentTrim` -> `[0,0]`
    (FeynmanTrick/DiffExpIntegration.m:1186-1187).
13. `exact_rational_passthrough` — chain
    `ESDivide[ESTimes[ESAdd[a, b], c], d]` over series built from
    `{1/3, 22/7, -5/11, ...}`: every coefficient of the result has head
    Integer or Rational (no Real anywhere), and one hand-computed
    coefficient matches exactly.  (No forced N — I-6.)
14. `precision_300_digits` — coefficients `N[1/3, 300]` etc.; after an
    add/mul/div chain every inexact result coefficient has
    `Precision >= 290`; `ESSameQ` against the exact-rational version of the
    same chain is True at matchTol.  (M1 gate: "Laurent division
    window-shift semantics vs exact rationals at 300 digits", RewritePlan M1.)
15. `shift_moves_both_edges` — `ESShift[ESNew[0, {1,2}], -3]` has window
    `[-3, -2]` and unchanged coefficients; `ESShift[s, 3]` then `ESShift[_, -3]`
    round-trips to `s` exactly.
16. `truncate_no_pad` — `ESTruncate[ESNew[0,{1,2,3}], 1]` -> window `[0,1]`,
    coeffs `{1,2}`; `ESTruncate[ESNew[0,{1}], 5]` aborts
    (ERR-TRUNCATE-EXTEND); `ESTruncate[ESNew[0,{1,2}], -1]` aborts.
17. `from_expression_exact_laurent` —
    `ESFromExpression[(1+eps)/(eps^2 (1-eps)), eps, 3]` -> window `[-2, 3]`,
    coefficients `{1, 2, 2, 2, 2, 2}` exact integers.
18. `from_expression_gamma_prefactor` —
    `ESFromExpression[Gamma[1+eps], eps, 2]` -> window `[0,2]`, coefficients
    `{1, -EulerGamma, EulerGamma^2/2 + Pi^2/12}` exact symbolic (no `N`;
    compare with `Simplify[... - expected] === 0`).  Pins the FT gamma
    prefactor path (FeynmanTrick/DiffExpIntegration.m:571, 880).
19. `from_expression_loud_failures` — `ESFromExpression[Exp[1/eps], eps, 3]`
    aborts (ERR-EXPAND-FAIL; old code would have silently fallen back,
    FeynmanTrick/DiffExpIntegration.m:1207-1217);
    `ESFromExpression[Sqrt[eps], eps, 3]` aborts naming the fractional
    leading power 1/2 (old code Floor-ed it, :1224).
20. `coefficient_list_contract` — `s = ESNew[-1, {4, 5, 6}]` (window
    `[-1,1]`): `ESCoefficientList[s, -2, 1] === {0, 4, 5, 6}`;
    `ESCoefficientList[s, 0, 2]` aborts (ERR-RANGE);
    `ESCoefficientList[s, 0, 1]` aborts (ERR-DROP-BELOW: dropping the
    nonzero `4` at k=-1); after `t = ESTrim[ESNew[-1, {0, 5, 6}]]`,
    `ESCoefficientList[t, 0, 1] === {5, 6}` (no error — the dropped range is
    genuinely zero).
21. `same_q_tolerance_and_windows` — series equal on shared range but with
    different CompleteMax compare True; perturbing one coefficient by
    `10 * matchTol * scale` flips to False; symbolic equal coefficients
    (`b1 - b2` vs `-(b2 - b1)`) compare True via PossibleZeroQ.
22. `error_context_propagation` — `Block[{$ESErrorContext = "chart x0=1,
    sector (a=-1,b=2,p=0), order k=7"}, ESCoefficient[ESNew[0,{1}], 3]]`
    aborts with a message containing that exact string.
23. `window_object_shape` — `ESWindow[ESNew[2, {1,2}]]` is EXACTLY
    `<|"Min" -> 2, "CompleteMax" -> 3|>` (key names and order verbatim per
    RewritePlan 3.1) so it can be stored into `LocalSolution["EpsWindow"]`
    unmodified.

Coverage notes: tests 4, 5, 6, 16 jointly enforce I-4 (no padding); 6, 8, 9
pin the window-shift semantics from three directions (numeric, trailing-zero
exact, symbolic pseudo-resonance); 11 and 12 are direct regression pins
against the two old-code Laurent layers.

---

## 9. LINE BUDGET

RewritePlan 3.2 allocates ~300 lines to EpsSeries.m (core total ~3.3k,
ceiling 3.5k).  Indicative breakdown: package boilerplate + usage messages
~45; validator/choke-point/`esError` ~35; accessors ~30; add/scale/shift
~35; times ~25; invert/divide ~30; truncate/trim ~25; predicates
(`ESCoeffZeroQ`, `ESSameQ`) ~35; conversions ~35; `ESMap` + context var ~10.
Total ~305.

If over budget, cut IN THIS ORDER (each cut must keep the error contract
intact):

1. Reduce `ESAdd` to binary only (callers fold) — saves ~5.
2. Drop `ESZero` (callers write `ESNew[kmax, {0}]`; keep the canonical-zero
   RULE in ESTrim) — saves ~8.
3. Drop `ESMap` (API.m reconstructs via `ESNew[ESMinPower[s],
   f /@ ESCoefficientList[s, ESMinPower[s], ESCompleteMax[s]]]`) — saves ~10.
4. Drop `ESCoefficientList`, moving the ERR-DROP-BELOW check to the single
   API.m output site (the only caller allowed to slice) — saves ~20.  This
   is the last acceptable cut; the arithmetic core, the window rules, the
   relative zero test, ESTrim, ESSameQ, the conversions, and ALL errors are
   not negotiable.

What may NOT be cut to save lines: any invariant, any error, the n-ary
window rules table, the exactness guarantee, or test coverage.

---

## 10. OPEN QUESTIONS

- O-1 Error mechanism unification: this spec uses print-then-`Abort[]`
  (matching old `ReportError`, DiffExp/Utilities.m:66) behind one internal
  helper.  If the orchestrator standardizes a library-wide helper (e.g. in
  Tolerances.m or Config.m), only `esError`'s body changes.  Decision owner:
  orchestrator at M1 merge.
- O-2 Exact casing/naming of Tolerances exports (`laurentLeadTol`,
  `matchTol`, working-precision accessor) — defer to Docs/specs/Tolerances.md;
  this spec binds to the RewritePlan 3.2 names semantically.
- O-3 `ESSameQ` degenerate shared range (one series' CompleteMax below the
  other's Min): spec'd as comparing over the guaranteed-zero region (usually
  trivially True).  If M3 matching tests want a loud error for "windows do
  not overlap meaningfully", add it then — do not add speculatively.
- O-4 Machine-precision coefficients: spec'd as ACCEPTED (the module is
  precision-agnostic; API.m owns the 2xWP raise per the ledger seed,
  DiffExp/Transport.m:541-562).  If M4 parity shows machine numbers leaking
  into the core, tighten ERR-BAD-CONSTRUCT to reject `MachinePrecision`
  coefficients.
- O-5 Coefficient normalization knob: `Together@*Expand` on non-numeric
  coefficients is spec'd unconditionally (N-4).  If the R6/M4 benchmark
  shows `Together` dominating on large x-rational coefficients, a
  SectorSeries-level batching strategy is the fix — do NOT add a "skip
  normalization" option here (it would make results
  representation-dependent).
- O-6 Whether SectorSeries stores per-(n,comp) EpsSeries objects or one
  shared window + raw 3D array using these rules (section 3.2 leaves it
  open).  EpsSeries' contract is identical either way; SectorSeries' spec
  decides the representation.
