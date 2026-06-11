# Module spec: DiffExp2/Indicial.m

Status: M0 deliverable (RewritePlan.md section 6, agent tasks 4-11).
Amended at the M0 spec-amendment pass per Docs/specs/DECISIONS-M0.md
(binding), REVIEW-math.md (D5, D6c, D8a, D17) and REVIEW-minimalism.md
(defects 4, 5, 14, 15, 33, 35).
Authority: Docs/RewritePlan.md (v2) is the execution contract; this spec
binds the implementation of `Indicial.m`.  An implementation agent must be
able to build the module from THIS file, RewritePlan.md, and the old code
cited here — nothing else.  Where this spec quotes RewritePlan section 3.1
data structures, the quote is verbatim and normative.

Context name: `DiffExp2`Indicial``.
Line budget: ~300 (RewritePlan 3.2; see section 9).
Position in the acyclic module order:
`Tolerances < Config < EpsSeries < SectorSeries < Indicial < Solve <
Transport/Integrate < API`.

---

## 1. PURPOSE

Indicial.m turns the exact eps-rational differential-equation matrix at a
chart center into the exact algebraic data that fixes the local solution
structure: it extracts the t-Laurent pole data of the chart-local matrix,
applies a ported Moser/shearing rank reduction (lattice saturation in the
style of old `FuchsianizeLocal`/`TrimFuchsianLattice`,
DiffExp/LocalSeries.m:15-16, 183-253) when the pole order exceeds one,
enforces the I1 contract — the characteristic polynomial of the residue
A_{-1}(eps) must factor exactly as a product of (lambda − a_i − b_i·eps)
over rationals/algebraics, LOUD ERROR otherwise — detects Jordan/confluent
structure by exact linear algebra over Q(eps) (never numerically), and
partitions the resulting {a, b, p} sector specs into true-resonance
log-chain families (same b, integer-spaced a) and pseudo-resonance
joint-solve families (integer-spaced a, different b, per I2).  Its output
is the sole source of sector tags in the library: tags are exact data
derived here once per chart and never reconstructed numerically (disease
D2's cure).  Indicial.m performs NO numerics, consumes NO tolerances, and
contains NO fallback of any kind: every chart is either classified exactly
or rejected with a loud, fully-attributed error.

---

## 2. PUBLIC SYMBOLS

All matrices/vectors below have entries that are exact rational functions
of their stated variables with exact numeric coefficients (Integer,
Rational, or exact algebraics: `Root`, `AlgebraicNumber`, Gaussian
rationals).  `t` is the chart coordinate with the (possible) singularity
at t = 0; the caller (Solve.m's `PrepareChart`, DEC-7) has already
applied the chart map.  `eps` is the dimensional regulator symbol.
Everything in this module is generic-eps algebra: valuations, ranks,
spectra, and nilpotency are computed over the field Q(alpha)(eps)
(alpha = algebraic numbers present in the input), never at a specialized
eps value.

PUBLIC SURFACE (REVIEW-minimalism defect 4: cut C-2 taken up front):
the exported symbols are `ChartIndicial`, `MatrixPoleData`,
`FuchsianReduce`, and `EpsDegenerateFamilies` (DEC-7).
`AffineSpectrum`, `JordanChains`, and `PartitionResonanceFamilies`
(2.4-2.6) are PACKAGE-PRIVATE — specified below for implementability
and tested through ChartIndicial's output fields.

### 2.1 ChartIndicial — the orchestrator

```
ChartIndicial[A_?MatrixQ, t_Symbol, eps_Symbol, chartRef_Association]
  -> IndicialData (section 3.6)
```

- `A`: d x d matrix, f'(t) = A(t, eps) . f(t) (d/dt convention), entries
  exact rational in (t, eps).  d >= 1.
- `chartRef`: error-attribution record (section 3.1).  Opaque to the
  algebra; every error payload embeds it.

Behavior:
1. Validate input exactness (error E1 otherwise; see section 5).
2. `pole = MatrixPoleData[A, t]`.
3. If `pole["PoleOrder"] <= 0`: return the trivial regular-chart
   IndicialData (one family whose single EigRecord is a = 0, b = 0,
   Multiplicity = d, BlockSizes = ConstantArray[1, d], Chains = the
   standard basis vectors; d sector specs, each
   `<|"a"->0,"b"->0,"p"->0|>`; identity gauge).  This is the SAME code
   path as singular charts — the degenerate case.  Per DEC-6
   (superseding the RewritePlan 3.1 one-sector phrasing): a regular
   chart is ONE (a=0,b=0,p=0) sector FAMILY with d-dimensional
   coefficient space; Solve.m's basis completeness and Transport.md's
   T-1 assert match this.
4. If `pole["PoleOrder"] == 1`: residue R = coefficient of t^-1 of A;
   skip reduction (identity gauge).
5. If `pole["PoleOrder"] >= 2`: `red = FuchsianReduce[A, t, eps,
   chartRef]` (errors E3/E4/E5 may fire); R = red["Residue"].
6. `spec = AffineSpectrum[R, eps, chartRef]` (error E2 may fire).
7. `spec = JordanChains[R, spec, eps, chartRef]` (error E6 may fire).
8. `fams = PartitionResonanceFamilies[spec, chartRef]`.
9. Assemble and return IndicialData; run the section-4 invariant asserts.

An apparent singularity (function analytic, matrix singular: the box L2
beta-root chart at t* = 7/11, Docs/FeynmanTrickBoxFamilyStatus.md:80-88)
goes through EXACTLY this path.  There is no "the solution is analytic
here" shortcut anywhere in this module.

### 2.2 MatrixPoleData — exact t-Laurent head of the matrix

```
MatrixPoleData[A_?MatrixQ, t_Symbol]
  -> <|"PoleOrder" -> r_Integer,
       "Coefficients" -> <|(-r) -> A_{-r}, ..., (-1) -> A_{-1}|>|>
```

- `r` = max over entries of (pole order in t at t = 0) computed exactly:
  per entry, `Cancel[Together[entry]]`, then r_entry =
  (t-min-degree of denominator) − (t-min-degree of numerator), i.e. the
  negative of the old `LocalZValuation` (DiffExp/LocalSeries.m:90-99);
  r = Max[0, max over entries].  Entries identically zero contribute
  nothing (old `ZeroQGeneric` exact test, LocalSeries.m:32-33).
- `Coefficients`: exact matrix Laurent coefficients of A in t for the
  pole part only, each a d x d matrix over Q(alpha)(eps).  Computed by
  the rational-Laurent recursion of old `LaurentCoefficientsRational`
  (LocalSeries.m:255-299) — NEVER by `Series`/`SeriesCoefficient` on
  potentially inexact objects, and never with asymptotic objects.
- If r computed from naive degrees but the candidate leading coefficient
  matrix A_{-r} is identically zero after cancellation, r is decremented
  until the leading coefficient is nonzero (cancellation is exact, so
  this terminates; a fully holomorphic matrix returns r = 0 and empty
  `Coefficients`).  The old behavior of bailing out with `$Failed` when
  "residue is identically zero" (DiffExp/IntegrationStrategies/
  Recurrence.m:129-132, ResonantRecurrence.m:41-44) is REPLACED: a zero
  residue with no deeper pole means the chart is REGULAR and is
  classified as such, loudly-correctly, not refused.

### 2.3 FuchsianReduce — ported Moser/shearing rank reduction

```
FuchsianReduce[A_?MatrixQ, t_Symbol, eps_Symbol, chartRef_Association]
  -> ReductionData (section 3.5)
```

Reduces a chart-local matrix with t-pole order r >= 2 to an equivalent
Fuchsian (simple-pole) system by a meromorphic gauge transform, exactly.
This is the port of `FuchsianizeLocal` + `TrimFuchsianLattice`
(DiffExp/LocalSeries.m:183-216, 218-253) with the call-site hardening of
`BuildFuchsianizedRecurrenceData`
(DiffExp/IntegrationStrategies/ResonantRecurrence.m:131-217) folded in.

ALGORITHM (concrete; implement exactly this):

Step 0 — theta form.  M(t, eps) := t * A(t, eps), entries
`Cancel[Together[...]]`-normalized.  Work in theta form throughout:
theta f = M f where theta = t d/dt.  The system is Fuchsian iff M is
holomorphic at t = 0 (matrix min t-valuation >= 0); then the d/dt residue
is M(0).  (Old convention: LocalSeries.m:51-53, 130-131.)

Step 1 — Moser pre-check.  Let v = min t-valuation of M (v <= -1 here,
since r = 1 - v >= 2).  Let L = (coefficient of t^v in M) = A_{-(1-v)}.
Require L nilpotent: `L^d == 0` exactly (MatrixPower + exact zero test).
If L is NOT nilpotent the singular point is genuinely irregular
(exponential local behavior; Moser's necessary condition fails): raise
E3.  No saturation is attempted.  This pre-check is NEW relative to the
old code (which just ran to MaxSteps); it converts the common irregular
case into an immediate, precisely-worded error.

Step 2 — lattice saturation loop (port of FuchsianizeLocal,
LocalSeries.m:188-216).  T := IdentityMatrix[d].  For step = 0, 1, ...,
MaxSteps-1 (MaxSteps = 200, the old default, LocalSeries.m:184):
  a. B := LocalConnectionMatrix(M, T) = Cancel[Together[
     Inverse[T] . (M . T - thetaT)]] with thetaT = t * D[T, t]
     (LocalSeries.m:130-131).
  b. If min t-valuation of B >= 0: saturation done, go to Step 3.
  c. Pick the column j with the most negative column min-valuation
     (LocalSeries.m:206-207; ties: smallest j).
  d. Adjoin u = B[[All, j]] to the lattice (port of
     AdjoinVectorToLattice, LocalSeries.m:133-181):
       - vals_i = t-valuation of u_i; m = Min[vals] (m < 0 here);
       - c_i = leading t-coefficient of u_i if vals_i == m else 0
         (LocalLeadingCoefficient, LocalSeries.m:101-116);
       - pivot = first i with c_i not identically zero.  c is by
         construction not all zero, so a missing pivot is an internal
         inconsistency: raise E6 (the old code Throw'd a "NoPivot"
         Failure, LocalSeries.m:151-158 — keep the loudness, name the
         chart);
       - Emat = IdentityMatrix[d]; Emat[[pivot, pivot]] = 1/c_pivot;
         Emat[[i, pivot]] = -c_i/c_pivot for i != pivot
         (entries are rational functions of eps and algebraics — exact);
       - uNew = Cancel /@ (Emat . u); U = IdentityMatrix[d] with column
         `pivot` replaced by uNew;
       - T := Cancel[Together[(T . Inverse[Emat]) . U]].
  If the loop exhausts MaxSteps without reaching min-valuation >= 0:
  raise E4 (non-termination).  NOTE: false success is impossible — the
  loop only returns when B is genuinely holomorphic, an exact test.
  CORRECTNESS PRECONDITION (the banana lesson): every quantity in this
  loop is EXACT.  The old code had to `Rationalize` numeric theta
  matrices before reduction because "exact pole-cancellation tests never
  fire on floats, so the trim pass ratcheted lattice columns by x once
  per pass up to x^50-scale entries while the precision collapsed to ~15
  digits" (ResonantRecurrence.m:143-154;
  Docs/FeynmanTrickBananaStatus.md:229-235).  DiffExp2 inputs are exact
  by construction, so the port ASSERTS exactness (E1) instead of
  rationalizing.  `Rationalize` must not appear in this module.

Step 3 — trim pass (port of TrimFuchsianLattice,
LocalSeries.m:223-253, WITH ONE MANDATORY AMENDMENT).  Purpose: remove
unnecessary column poles from the lattice so the gauge's t-degrees — and
hence the integer a-shifts and coefficient-window bloat seen by
downstream consumers — stay minimal.
For pass = 1..MaxPasses (= 50, old default, LocalSeries.m:219):
  for each column p WITH min t-valuation of T[[All, p]] < 0
  (PRECONDITION — see below): trial T' = T . D_p with D_p = identity
  except (p,p) -> t; recompute B' = LocalConnectionMatrix(M, T'); if B'
  still has min t-valuation >= 0, accept (T := T').  Stop when a full
  pass changes nothing.  All tests exact.
AMENDMENT RATIONALE (do not port verbatim): the old code trials EVERY
column (LocalSeries.m:234-247).  On systems with a decoupled direction,
multiplying a pole-FREE column by t preserves Fuchsianity while merely
shifting that solution's exponent down by 1, so the old loop accepts the
"trim" again every pass and ratchets until MaxPasses — in EXACT
arithmetic (hand-verifiable on test T-12's matrix: after saturation,
column 2 of T is pole-free and the old loop would drive the reduced
eigenvalue 0 -> −50 before the degree guard kills the run).  The banana
campaign saw the float-amplified version of this ratchet
(Docs/FeynmanTrickBananaStatus.md:229-235); exactness alone does NOT
fully remove it.  The precondition (only trial columns that actually
carry a pole — the documented purpose, LocalSeries.m:16 "removes
unnecessary column poles") removes it structurally.  (Trim is cut
candidate #1 under budget pressure; see section 9.)

Step 4 — degree guard (port of the degenerate-lattice guard,
ResonantRecurrence.m:179-205).  maxDeg = max over T entries of
Max[Exponent[num, t], Exponent[den, t]] (after Together).  Require
maxDeg <= 4 d + 8.  In exact arithmetic the numeric-ratcheting failure
mode this guard was written for cannot occur, so a violation here means
an internal bug or a pathological input: raise E5 (the old code's
PrintWarning + silent `$Failed`, ResonantRecurrence.m:197-204, is
FORBIDDEN).

Step 5 — finalize.  TInv := Cancel[Together[Inverse[T]]]; assert
T . TInv == Identity exactly.  B := LocalConnectionMatrix(M, T) (final,
exact).  Residue := B /. t -> 0 (B is holomorphic at 0; substitute and
Cancel — exact).  Return ReductionData (section 3.5) carrying T, TInv, B,
Residue, Steps, the original pole order, and `Trimmed`.

The reduced system the downstream solver actually solves is
theta g = B(t, eps) g, equivalently g' = (B/t) g — pole order exactly 1
(asserted).  The ORIGINAL system's solutions are f = T . g; the
composition contract is section 3.7.

### 2.4 AffineSpectrum — the I1 contract

```
AffineSpectrum[R_?MatrixQ, eps_Symbol, chartRef_Association]
  -> {EigRecord..} (section 3.3)
```

`R` = the (reduced) residue, d x d over Q(alpha)(eps).

ALGORITHM:
1. chi(lambda) = CharacteristicPolynomial[R, lambda], coefficients
   Cancel/Together-normalized in Q(alpha)(eps).
2. Clear eps-denominators: chiC = chi * (polynomial LCM of coefficient
   denominators in eps); chiC in Q(alpha)[lambda, eps].  (The cleared
   denominator is eps-only and cannot vanish identically; it carries no
   roots in lambda.)
3. `FactorList[chiC]` (with `Extension -> Automatic` when algebraic
   coefficients are present).  For each non-constant irreducible factor
   fct with exponent e:
   - if Exponent[fct, lambda] == 0: it is part of the cleared
     denominator content; skip;
   - if Exponent[fct, lambda] == 1: root = Together[-c0/c1] where
     c1, c0 are the lambda-coefficients;
   - if Exponent[fct, lambda] >= 2: attempt exact roots via
     `Solve[fct == 0, lambda]`; each solution is a candidate root
     (algebraic a, b arise here — e.g. eigenvalues +-Sqrt[2] (1 + eps)).
     If Solve cannot produce exact roots, raise E2 with the factor.
4. AFFINE VERIFICATION per root r(eps): compute b = D[r, eps]; require
   `D[b, eps] === 0` after Together/RootReduce (i.e. r is affine in eps);
   then a = (r /. eps -> 0), and assert
   `RootReduce[Together[r - (a + b*eps)]] === 0`.  Any failure: E2.
   The check is on the EIGENVALUES, never on the residue entries —
   "linear-in-eps residue does NOT imply affine eigenvalues"
   (RewritePlan section 2, I1; reviewer counterexample [[0,1],[eps,0]]
   -> +-Sqrt[eps]; Docs/reviews/rewrite_plan_review_3lens.json, math
   finding 2 and legacy finding 1).  Conversely, residues with
   eps-RATIONAL or eps-polynomial entries of any degree are accepted as
   long as every eigenvalue verifies affine.
5. Multiplicity of (a, b) = factor exponent x (number of times the root
   occurs in the factor) summed over factors; merge records with
   identical exact (a, b) (equality via RootReduce).
6. CLOSURE ASSERT (decisive, cheap): expand
   lc(eps) * Product_i (lambda - a_i - b_i eps)^{m_i} - chiC
   and require identically zero (Together/RootReduce), where lc is
   chiC's leading lambda-coefficient.  Sum_i m_i == d asserted.

### 2.5 JordanChains — exact confluence detection

```
JordanChains[R_?MatrixQ, spectrum:{__Association}, eps_Symbol,
             chartRef_Association]
  -> {EigRecord..}  (input records augmented with Jordan data)
```

For each distinct eigenvalue lam = a + b*eps with multiplicity m:
1. Nested nullities: for j = 1, 2, ... compute n_j = d - rank over
   Q(alpha)(eps) of (R - lam I)^j, until n_j == m (j is then the largest
   block size).  Rank by exact Gaussian elimination / `RowReduce` on the
   symbolic matrix with `Cancel`-normalized pivots and exact zero tests
   (`Together`/`RootReduce`-based `PossibleZeroQ` confirmed by exact
   cancellation).  NEVER pass a `Tolerance` option; NEVER numericize.
2. Block sizes: #blocks of size >= j equals n_j - n_{j-1} (n_0 = 0).
   If the nullity sequence is inconsistent (non-monotone increments,
   or n_j never reaches m for j <= m): raise E6 — with exact arithmetic
   this is an internal bug, not a data condition.
3. Chains: top vectors of length-k chains are exact null vectors of
   (R - lam I)^k independent of Null((R - lam I)^(k-1)) and of vectors
   already claimed at that level; chain members v_{q-1} = (R - lam I).v_q.
   This mirrors the level bookkeeping of `NumericJordanData`
   (ResonantRecurrence.m:266-339) executed in exact arithmetic.
   Normalization convention (for test reproducibility): each chain-top
   vector is scaled so its first nonzero component is 1, entries
   Cancel-normalized.
4. Verify each chain exactly: (R - lam I).v_q == v_{q-1}, and
   (R - lam I).v_1 == 0.  Failure: E6.

The eps -> 0 collision of EIGENVECTORS across different-b eigenvalues
(the log x = (x^(2 eps) - 1)/(2 eps) class) is NOT Jordan structure at
generic eps; it IS detected here, exactly, by the per-family
EpsZeroDegeneracy rank computation of section 2.6 step 7 and exported
through `EpsDegenerateFamilies` (2.7; DEC-7).  Division of labor
(REVIEW-math D8): Indicial DETECTS (exact rank at eps = 0), Transport
REMOVES (unimodular recombination, Transport.md 2.8).

### 2.6 PartitionResonanceFamilies — resonance / pseudo-resonance

```
PartitionResonanceFamilies[spectrum:{__Association},
                           chartRef_Association]
  -> {Family..} (section 3.4)
```

1. FAMILY RELATION: eigenvalue records i, j belong to the same family
   iff (a_i - a_j) is an exact integer (test on exact values; for
   algebraic a: `IntegerQ[RootReduce[a_i - a_j]]`).  b plays NO role in
   the family relation — pseudo-resonance (different b, integer-spaced
   a) deliberately lands in one family (RewritePlan section 2, I2).
   This is the transitive closure; with affine eigenvalues it is an
   equivalence relation per b-independent a-lattice.
2. Within a family, the SAME-b subgroups (exact b equality via
   RootReduce) are the true-resonance log-chain classes — the
   generalization of old `resonanceClasses`
   (ResonantRecurrence.m:440-460) from numbers to (a, b) pairs.
3. Per-solution log bound p (the generalized `maxLogPowers` formula,
   ResonantRecurrence.m:489-493): for the solution at eigenvalue
   (a, b), Jordan chain position q (0-based, q = 0 the eigenvector):
     p = q + Sum of multiplicities m(a', b) over same-b family members
             with a' > a, a' - a a positive integer.
   Cross-b family members contribute NOTHING to p: pseudo-resonant
   collisions are eps-Laurent shifts removed by the joint solve, not
   logs (I2: the eps-regular combination
   (S_i - S_j)/((b_i - b_j) eps) keeps tags exact and log-free).
4. Sector specs: one per solution, `<|"a" -> a, "b" -> b, "p" -> p|>`,
   p as in 3.  These p are ANSATZ BOUNDS (allocation for Solve.m's
   log-chain ansatz); realized Sector["p"] tags in the final
   LocalSolution may be smaller (Solve decides; see section 3.4 note).
5. Collision list (consumed by Solve.m's recursion bookkeeping and the
   3.4 budget asserts): for every ordered pair (i, j) in the family
   with n := a_j - a_i a non-negative integer and (a_i, b_i) != (a_j, b_j):
     <|"n" -> n, "LowerIdx" -> i, "UpperIdx" -> j,
       "DeltaB" -> b_i - b_j,
       "Type" -> If[b_i === b_j, "Log", "LaurentShift"]|>
   For n == 0 (equal a, different b) both ordered pairs satisfy the
   condition; record ONLY the canonically ordered one (i before j in
   the family's sector sort) — one collision, one record.  (REVIEW-math
   D6c, resolved by its second option per DEC-7: n = 0 collisions
   ARRIVE canonically ordered, and Solve.m's `PrepareChart` symmetrizes
   to both ordered pairs where its recursion bookkeeping needs them;
   this single canonical record and T-11 are normative.)
   "Log" entries (same b, n > 0) are where the recursion operator
   L_n = (lam_i + n) I - R is singular identically in eps -> log bump.
   "LaurentShift" entries are where det L_n carries the exact factor
   ((a_i + n - a_j) + (b_i - b_j) eps) = (b_i - b_j) eps -> the
   pseudo-resonant 1/eps that the joint solve must absorb
   (pseudoResonanceShift = 0 is ASSERTED downstream, RewritePlan 3.4).
6. Family record fields: section 3.4.  `JointSolve -> True` iff the
   family contains >= 2 distinct exact b values.  `LogMax` = max p over
   the family's sector specs.
7. `"EpsZeroDegeneracy" -> r0_Integer` per family (REVIEW-math D8a):
   assemble the family's chain-top vectors v_i(eps) (columns),
   normalize each by its eps-valuation (exact), substitute eps -> 0
   (exact), and set r0 = (number of columns) − rank over Q(alpha) of
   the resulting matrix.  r0 > 0 iff the eps -> 0 eigenvectors collide;
   r0 is the recursion depth Transport.m's RecombineBasis needs.
   Computed only for families with >= 2 distinct b (elsewhere
   r0 := 0).  All-exact, no tolerance (invariant I-8 applies).

Families are returned sorted by (min a at eps = 0, then b) for
deterministic output; sectors within a family sorted by (a, b, p).

### 2.7 EpsDegenerateFamilies — eps = 0 degeneracy selector (DEC-7)

```
EpsDegenerateFamilies[data_Association (* IndicialData, section 3.6 *)]
  -> {<|"FamilyIndex" -> _Integer, "EpsZeroDegeneracy" -> _Integer|>..}
```

Thin exact selector, no new mathematics: one record per Family whose
`EpsZeroDegeneracy` r0 is > 0 (same a mod Z, distinct b, colliding
eps = 0 eigenvectors — the 2.6 step 7 field), `FamilyIndex` indexing
into `data["Families"]`; empty list when no family is degenerate.
Transport.m's RecombineBasis is tag-driven off this output and performs
NO numeric degeneracy detection (DEC-7; REVIEW-math D8; REVIEW-
minimalism defect 14).

---

## 3. DATA CONTRACTS

### 3.1 ChartRef (consumed; constructed by Solve.m/Transport.m)

```
ChartRef = <|
  "Name"     -> _String,   (* e.g. "banana.L1.xx1=0" or "box.L2.seg12" —
                              unique, human-readable; appears verbatim in
                              every error *)
  "Center"   -> exact x0,  (* in the ORIGINAL variable; exact number or
                              algebraic; Infinity allowed for Mobius-mapped
                              charts *)
  "Variable" -> _Symbol    (* the original line/integration variable *)
|>
```
Extra keys are permitted and ignored.  Indicial never reads anything but
these three, and only for error payloads and logs.

### 3.2 SectorSpec (produced)

The tag part of RewritePlan 3.1's Sector, quoted verbatim minus the
coefficient array:

```
  Sector = <|
    "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
    "Coeffs" -> c[[k, n, comp]]
  |>
```

Indicial emits `SectorSpec = <|"a" -> _, "b" -> _, "p" -> _Integer|>` —
exactly the first three fields, same names, same types ("a" may also be
exact algebraic when chart centers are algebraic; RewritePlan R2 covers
this: performance flag, not correctness).  Solve.m attaches "Coeffs".
EpsWindow note (RewritePlan 3.1 INVARIANTS, verbatim): "homogeneous
true-resonance sectors have kmin = −p BY CONSTRUCTION (log-chain weights
are eps-independent under the (eps Logx)^p/p! normalization); window
arithmetic must account, not error."  The p emitted here is the number
that rule keys off.

### 3.3 EigRecord (produced by AffineSpectrum, augmented by JordanChains)

```
EigRecord = <|
  "a" -> exact, "b" -> exact,        (* eigenvalue = a + b*eps *)
  "Multiplicity" -> _Integer,        (* algebraic multiplicity *)
  (* added by JordanChains: *)
  "BlockSizes" -> {_Integer..},      (* Jordan block sizes, descending;
                                        Total == Multiplicity *)
  "Chains" -> {{v1, v2, ..., vk}..}  (* one list per block, ordered
       eigenvector FIRST (v1), each v an exact d-vector over
       Q(alpha)(eps); (R - lam I).v_{q+1} == v_q, (R - lam I).v1 == 0 *)
|>
```

### 3.4 Family (produced)

```
Family = <|
  "Sectors"    -> {SectorSpec..},    (* one per basis solution; Length ==
                                        sum of member multiplicities *)
  "Members"    -> {EigRecord..},     (* the family's eigenvalue records *)
  "Class"      -> "Single" | "Confluent" | "TrueResonant" | "Pseudo",
       (* Single: one eigenvalue, diagonalizable, no integer partner.
          Confluent: one eigenvalue with a Jordan block (p > 0 from
            chains alone).
          TrueResonant: >= 2 same-b members at integer-spaced a.
          Pseudo: >= 2 distinct b values present (may ALSO contain
            true-resonant same-b subgroups; Pseudo wins the label and
            JointSolve goes True). *)
  "JointSolve" -> True | False,      (* True iff >= 2 distinct b values:
                                        Solve.m must solve the family
                                        jointly, choosing the eps-regular
                                        combination at each collision
                                        (I2 spec decision) *)
  "LogMax"     -> _Integer,          (* max p over Sectors: Solve.m's
                                        log-chain ansatz depth *)
  "EpsZeroDegeneracy" -> _Integer,   (* >= 0; section 2.6 step 7
                                        (REVIEW-math D8a): number of
                                        chain-top columns minus their
                                        exact rank at eps = 0.  r0 > 0
                                        iff the eps -> 0 eigenvectors
                                        collide; r0 = RecombineBasis'
                                        recursion depth.  Always 0 when
                                        fewer than 2 distinct b. *)
  "Collisions" -> {<|"n" -> _Integer, "LowerIdx" -> _Integer,
                     "UpperIdx" -> _Integer, "DeltaB" -> exact,
                     "Type" -> "Log" | "LaurentShift"|>..}
|>
```

NORMATIVITY (REVIEW-minimalism defect 15; DEC-7): the Family and
EigRecord shapes above — including the collision keys
n/LowerIdx/UpperIdx/DeltaB/Type — are consumed VERBATIM by Solve.m's
`PrepareChart` (Solve's former From/To/Offset names are renamed to
these).  Solve.m itself builds the per-family spectral frame V/VInv,
CollisionDepth, and the ChartSystem from this output; Indicial emits
exactly what this section specifies, nothing more.

NOTE ON p SEMANTICS: `Sectors[..,"p"]` is the allocation bound for the
ansatz (old `MaxLogPowers` role, ResonantRecurrence.m:489-493).  Solve.m
constructs the log-chain ansatz to depth LogMax and may emit final
Sector tags with smaller realized p when top log weights are EXACTLY
zero; it must NOT trim on numerically-small coefficients (that is the
matching-side "numerical-zero leading-coefficient skipping" lesson,
which belongs to Transport.m, not to tag derivation).

### 3.5 ReductionData (produced by FuchsianReduce)

```
ReductionData = <|
  "PoleOrder"    -> _Integer,        (* ORIGINAL pole order r of A, >= 2;
                                        1 or 0 when reduction was skipped *)
  "Gauge"        -> T,               (* d x d exact rational in (t, eps);
                                        Identity when no reduction *)
  "GaugeInverse" -> TInv,            (* exact; T.TInv == Id asserted *)
  "ThetaMatrix"  -> B,               (* reduced system: theta g = B g,
                                        B holomorphic at t = 0, exact *)
  "Residue"      -> B /. t -> 0,     (* exact d x d over Q(alpha)(eps) *)
  "Steps"        -> _Integer,        (* saturation steps used *)
  "Trimmed"      -> True | False,
  "GaugeValuation" -> _Integer       (* min t-valuation over T entries:
                                        a lower bound on the integer
                                        a-shift composition introduces *)
|>
```

### 3.6 IndicialData (produced by ChartIndicial; consumed by Solve.m)

```
IndicialData = <|
  "Chart"      -> ChartRef,
  "Dimension"  -> d,
  "PoleData"   -> MatrixPoleData result (section 2.2),
  "Reduction"  -> ReductionData,     (* identity-gauge record when pole
                                        order <= 1 *)
  "Residue"    -> R,                 (* == Reduction["Residue"]; the
                                        matrix all spectrum data refers
                                        to: the REDUCED system's residue *)
  "Spectrum"   -> {EigRecord..},
  "Families"   -> {Family..},
  "Regular"    -> True | False       (* True iff PoleOrder <= 0: one
                                        family, d (0,0,0) sector specs
                                        (DEC-6) *)
|>
```

### 3.7 Gauge composition contract

(Binding on Solve.m; stated here because the gauge originates here.)

All sector tags, families, and collision data in IndicialData refer to
the REDUCED system theta g = B g.  Solve.m solves the reduced system and
composes f = T . g via SectorSeries.m's multiply-by-rational operation
(T entries are rational in (t, eps); t-poles of T sit at t = 0 by
construction and shift `a` by integers; T's eps-denominators — the
1/c_pivot factors of Step 2d — go through EpsSeries.m's Laurent-division
window-shift semantics, honestly).  Consequences, all binding:
- The final LocalSolution's sectors are the formal image of the reduced
  sectors under T: a-tags shift by integers >= GaugeValuation; b and the
  family structure (b-set, JointSolve, LogMax) are GAUGE-INVARIANT and
  must be identical before and after composition (cheap assert).
- The always-on ODE residual spot-check (RewritePlan 3.1 INVARIANTS)
  for the composed LocalSolution runs against the ORIGINAL matrix A,
  not against B.  This is what "downstream consumers see the ORIGINAL
  system's solutions" means operationally; the old equivalent was
  `ApplyGaugeToSolution` + `ValidateFiniteWidthSolution`
  (LocalSeries.m:886-944, 952-1022).
- Boundary/matching data are always in the ORIGINAL frame; Solve.m maps
  incoming values with TInv at most once per chart, exactly.
- Window bookkeeping: composition may LOWER kmin (eps-denominators in T)
  and lower the leading a; both are recorded through the standard
  EpsWindow/TWindow arithmetic — never padded, never clamped.

---

## 4. INVARIANTS (always on, cheap)

I-1  Input exactness: every matrix entry is free of Real/Complex inexact
     numbers and free of SeriesData heads (in t AND in eps).  Checked at
     every public entry point (E1 on failure).
I-2  Sum of multiplicities over Spectrum == d; sum of Length[Sectors]
     over Families == d; Total[BlockSizes] == Multiplicity per record.
I-3  Char-poly closure: lc * Prod (lambda - a_i - b_i eps)^{m_i} == chiC
     identically (section 2.4 step 6).
I-4  Chain identities: (R - lam I).v_{q+1} == v_q and
     (R - lam I).v_1 == 0, exactly, for every emitted chain.
I-5  Reduction validity: min t-valuation of ThetaMatrix >= 0;
     T . TInv == Identity exactly; recomputing
     LocalConnectionMatrix(M, T) reproduces ThetaMatrix exactly.
I-6  Regular chart <=> exactly d sector specs, all (0, 0, 0), one
     family (DEC-6, amending the RewritePlan 3.1 one-sector phrasing;
     scoped to ODE solutions — the products-with-IBP-coefficients
     exemption lives in SectorSeries.md, not here).
I-7  Every "a" with b == 0 in a Fuchsian-from-the-start FT chart is
     rational; algebraic a/b only ever enter via algebraic chart centers
     or algebraic char-poly factors — both exact (R2).
I-8  No tolerance, no `N[...]`, no `Rationalize`, no `Chop`, no
     `Tolerance ->` option anywhere in the module (grep-clean; a unit
     test greps the source — see T-21).
I-9  Determinism: identical input produces identical output (sorted
     families/sectors, canonicalized chain normalization, fixed
     tie-breaks).

All asserts raise through ``Tolerances`DE2Error`` (section 5.0; DEC-1)
with the chart name; none is downgradeable.

---

## 5. ERROR CONTRACT

### 5.0 Mechanism

All errors are raised via ``Tolerances`DE2Error[id, payload]`` (DEC-1;
REVIEW-math D17, REVIEW-minimalism defect 5): DE2Error prints a
one-line summary and performs
`Throw[Failure["DiffExp2", payload], "DiffExp2Error"]`; the only Catch
sites are API.m's entry points.  This module defines NO other error
mechanism.  The E-names below are the `id` argument; payload ALWAYS
contains:

```
<|"ID" -> <name below>, "Module" -> "Indicial",
  "Chart" -> chartRef["Name"], "Center" -> chartRef["Center"],
  "Variable" -> chartRef["Variable"], ...class-specific fields...|>
```

No error in this module is catchable-and-continuable by design: there is
no degraded result to continue WITH.  `Quiet`, `Check[..., fallback]`,
and silent `$Failed` returns are forbidden in this module (see 5.9).

### 5.1 E1 "InexactInput"

Fires when: any entry of an input matrix (A or R) contains an inexact
number, or a SeriesData head, or a non-rational non-algebraic
x/eps-dependence (e.g. Sqrt[t] — RewritePlan non-goal, loud at
LoadSystem but re-checked here).
Carries: the entry position {i, j}, the offending sub-expression, and —
for the SeriesData case — the hint that eps-TRUNCATED slice exports
(FeynmanTrick/MatrixExport.m:42-77 `ExportDiffExpMatrix`, the
d<var>_<k>.m format of DiffExp/MatrixLoading.m:29) cannot certify the I1
contract and that the exact full export (`ExportGeneralMatrix`,
MatrixExport.m:91-114, d<var>_full.m) is required (RewritePlan I1
PREREQUISITE).
Forbidden fallback: `Rationalize` of float entries (the old workaround,
ResonantRecurrence.m:151-154); per-eps-order re-expansion of SeriesData.

### 5.2 E2 "NonAffineEigenvalue" (the I1 contract violation)

Fires when: any root of the characteristic polynomial of the residue
fails the affine verification of section 2.4 step 4, or an irreducible
factor of lambda-degree >= 2 admits no exact root extraction.
Carries: the characteristic polynomial chiC (cleared form), the
offending irreducible factor, the non-affine root when one was computed,
and the text of risk R1 ("out-of-scope perturbative fallback documented"
— documented, NOT implemented).
Forbidden fallbacks (every tempting one, explicitly):
- numeric `Eigenvalues` + `Rationalize` snapping (old
  Recurrence.m:134-136, ResonantRecurrence.m:380-386);
- per-eps-order eigenvalue expansion ("solve at eps^0, correct
  perturbatively") — this is disease D2 reborn;
- treating the residue as if b = 0 ("keep towers explicit") — the old
  DecomposeSingularity salvage (Docs/FeynmanTrickBananaStatus.md:256-259)
  has no analogue here;
- accepting the affine PART of a mixed spectrum and dropping the rest.

### 5.3 E3 "IrregularSingularPoint"

Fires when: pole order r >= 2 and the leading Laurent coefficient
A_{-r} (equivalently the t^v coefficient of the theta matrix) is NOT
nilpotent (Moser necessary condition fails; local solutions contain
exp(c/t^(r-1)) behavior — genuinely out of scope, RewritePlan R11).
Carries: r, the leading coefficient matrix, and the nonzero entry
pattern of A_{-r}^d as the non-nilpotency witness (the optional
exact-eigenvalue witness is cut C-3, taken up front per
REVIEW-minimalism defect 4).
Forbidden fallback: attempting saturation anyway; truncating the pole.

### 5.4 E4 "RankReductionNonTermination"

Fires when: the saturation loop reaches MaxSteps (= 200) without
achieving min t-valuation >= 0.
Carries: MaxSteps, the final min t-valuation, the final gauge T's max
t-degree, and r.  (Old behavior: `Failure["FuchsianizationDidNotTerminate"]`
returned as a VALUE, LocalSeries.m:212-215, which
BuildFuchsianizedRecurrenceData converted to silent `$Failed`,
ResonantRecurrence.m:162-164 — both halves forbidden; the new module
throws.)
Forbidden fallback: returning the partially reduced system; retrying
with a larger bound silently.

### 5.5 E5 "DegenerateGaugeLattice"

Fires when: the Step-4 degree guard maxDeg > 4 d + 8 trips
(ResonantRecurrence.m:179-205 lesson).
Carries: maxDeg, the bound 4 d + 8, d, Steps.
Forbidden fallback: the old PrintWarning + `$Failed`
(ResonantRecurrence.m:197-204); accepting the oversized gauge.

### 5.6 E6 "InternalAlgebraInconsistency"

Fires when: pivot extraction finds no nonzero leading coefficient
(LocalSeries.m:151-158 class); Jordan nullity sequence inconsistent or
chains fail their exact identities (section 2.5); T.TInv != Id; the
char-poly closure assert I-3 fails.  These are bugs, not data: the
payload says so.
Carries: a stage tag ("Pivot" | "Nullity" | "Chain" | "GaugeInverse" |
"CharPolyClosure"), the eigenvalue (a, b) where applicable, and the
offending object.
Forbidden fallback: `LeastSquares`/`PseudoInverse` repair (the old
SolveRecurrenceStep chain, ResonantRecurrence.m:546-571 and 614-640,
including "numerically singular but no null space found - treat as
non-singular", :617-621 — none of that machinery may be ported).

### 5.7 E7 "BadShape"

Fires when: A is not square, d == 0, t === eps, chartRef missing a
required key.  Carries: the actual shapes/keys.

### 5.8 Loudness of WARNING-class events: there are none

This module defines NO warnings.  Everything is either fine or an error.
(The old code's verbosity-gated PrintWarning culture around this
functionality — e.g. Frobenius.m:42-45's "root of the indicial equation
is of degree greater than two ... not thoroughly tested" warning gated
by `IgnoreIndicialCheck` — is replaced by exact handling (algebraic
roots verified affine) or E2.)

### 5.9 Anti-fallback sweep (complete enumeration of tempting sites)

1. Inexact input -> Rationalize: FORBIDDEN (E1).  Old site:
   ResonantRecurrence.m:151-154.
2. Eigenvalues numerically + snap: FORBIDDEN (E2 path is exact).  Old
   sites: Recurrence.m:134-136; ResonantRecurrence.m:380-386, 423-425.
3. JordanDecomposition exact-fails -> numeric fallback chain:
   FORBIDDEN.  Old site: ResonantRecurrence.m:396-419 (three-stage
   fallback ending in ReportError).
4. Toleranced SVD null spaces / NumericJordanData: FORBIDDEN — exists
   only because old inputs were numeric.  Old site:
   ResonantRecurrence.m:219-351.
5. Applicability predicates returning False/$Failed under
   `Quiet[Check[...]]` so a dispatcher can try something else:
   FORBIDDEN — there is no dispatcher.  Old sites: Recurrence.m:103,
   170-172, 181-196, 198-201; ResonantRecurrence.m:16-49, 122-129;
   Dispatch.m:5-132.
6. IntegerQ on inexact eigenvalue differences (silently never resonant):
   moot under exactness, and the exactness assert is the guard.  Old
   sites: ResonantRecurrence.m:449 (correct only because inputs were
   snapped); the box-campaign defect record
   Docs/FeynmanTrickBoxFamilyStatus.md:109-111.
7. Zero residue -> "not applicable": FORBIDDEN; classify as regular
   (section 2.2).  Old sites: Recurrence.m:129-132,
   ResonantRecurrence.m:41-44.
8. Non-termination -> $Failed -> some other strategy: FORBIDDEN (E4).
   Old sites: LocalSeries.m:212-215 + ResonantRecurrence.m:162-164.
9. Degenerate gauge -> warn + $Failed: FORBIDDEN (E5).  Old site:
   ResonantRecurrence.m:197-204.
10. "Largest indicial root only" (solve for rMax and ignore the other
    roots): FORBIDDEN — all d roots are enumerated.  Old site:
    Frobenius.m:39-40.
11. Truncating the pole part / dropping t-Laurent heads below some
    size: FORBIDDEN; pole data is exact (section 2.2).
12. Accepting slice-truncated eps input "because the slices look
    exact": FORBIDDEN (E1 + the I1 prerequisite, RewritePlan I1;
    legacy review finding 1-2).

---

## 6. ABSORBED OLD CODE

Indicial.m replaces the following old functionality (the old library is
FROZEN as parity oracle; nothing below is edited, only superseded):

| Old code | Lines | Disposition |
|---|---|---|
| `PrepareSingularRecurrence` (residue extraction, numeric eigenvalues, diagonalizability rank test, positive-integer non-resonance test) | DiffExp/IntegrationStrategies/Recurrence.m:100-172 | absorbed: MatrixPoleData + AffineSpectrum + JordanChains, exact |
| `SingularRecurrenceApplicableQ` wrapper | Recurrence.m:198-201 | dropped — no applicability predicates |
| `GeneralSingularRecurrenceApplicableQ`, `ExpandedMatrixMinOrder` | DiffExp/IntegrationStrategies/ResonantRecurrence.m:16-49, 51-60 | absorbed: MatrixPoleData |
| `FuchsianizedSingularRecurrenceApplicableQ` | ResonantRecurrence.m:122-129 | absorbed: pole-order branch of ChartIndicial |
| `BuildFuchsianizedRecurrenceData` (rationalize-then-reduce, trim call, degree guard, gauge/inverse assembly) | ResonantRecurrence.m:131-217 | absorbed: FuchsianReduce (exactness assert replaces Rationalize; guard becomes E5) |
| `NumericJordanData` (toleranced generalized-eigenvector chains) | ResonantRecurrence.m:219-351 | replaced by exact JordanChains; survives only as the chain-bookkeeping blueprint |
| `ComputeResonanceStructure` (snapped spectrum, Jordan fallback chain, block structure, resonance classes, `maxLogPowers`, initial vectors) | ResonantRecurrence.m:353-527 | absorbed: AffineSpectrum + JordanChains + PartitionResonanceFamilies; `maxLogPowers` formula :489-493 generalized to (a,b) pairs (section 2.6 step 3); initial vectors :496-501 become exact Chains |
| `LocalZValuation`, `LocalLeadingCoefficient`, min-valuation helpers | DiffExp/LocalSeries.m:90-128 | ported (exact valuation arithmetic) |
| `LocalConnectionMatrix` | LocalSeries.m:130-131 | ported verbatim (exact) |
| `AdjoinVectorToLattice` | LocalSeries.m:133-181 | ported (Step 2d) |
| `FuchsianizeLocal` | LocalSeries.m:183-216 | ported (Step 2) + Moser pre-check (Step 1, new) |
| `TrimFuchsianLattice` | LocalSeries.m:218-253 | ported (Step 3); first cut candidate |
| `LaurentCoefficientsRational` | LocalSeries.m:255-299 | ported into MatrixPoleData (pole heads only; the full-window variant belongs to EpsSeries.m/Solve.m) |
| `ApplyGaugeToSolution`, `ValidateFiniteWidthSolution` | LocalSeries.m:886-944, 952-1022 | NOT in Indicial: the composition/validation CONTRACT (section 3.7) binds Solve.m/SectorSeries.m; cited here because the gauge originates here |
| Scalar indicial equation, largest-root choice, `IgnoreIndicialCheck` warning | DiffExp/Frobenius.m:32-45 | dropped (scalar Frobenius path deleted with the strategy stack); all-roots enumeration replaces largest-root |
| `DispatchStrategy` spectrum-class routing | DiffExp/IntegrationStrategies/Dispatch.m:5-132 | dropped — ONE solver consumes ONE classification (this module's output) |

The rest of LocalSeries.m (PrepareFiniteWidthData/RecursiveFiniteWidthSolve,
:424-726) is Solve.m's inventory, not Indicial's (RewritePlan I2).

### Numerical lessons that MUST survive the port

(Each verified in the old code and campaign docs; each becomes spec text
above and a test below.)

N-1 EXACT ARITHMETIC OR NOTHING in the lattice reduction.  Float inputs
    made the trim pass ratchet lattice columns to x^50-scale while
    precision collapsed (banana L1 fix 1,
    Docs/FeynmanTrickBananaStatus.md:229-235;
    ResonantRecurrence.m:143-154).  New form: E1 exactness assert; no
    Rationalize.  (Tests T-19, T-21.)
N-1b The trim ratchet is only PARTLY a float problem: trialing
    pole-free columns can ratchet in exact arithmetic too (spec
    discovery, verified by hand on T-12's matrix against
    LocalSeries.m:231-247's verbatim logic).  New form: the Step-3
    column-valuation precondition (section 2.3).  (Test T-12 fails
    against a verbatim port — it doubles as the regression test.)
N-2 Float noise splits degenerate eigenvalues by ~sqrt(noise); exact
    JordanDecomposition on rationalized noise returns ill-scaled
    eigenvectors (banana fix 2, FeynmanTrickBananaStatus.md:236-242;
    ResonantRecurrence.m:373-394).  The size-3 Jordan block at the
    banana L1 endpoints was MISSED entirely by the old exact-on-noise
    path.  New form: exactness makes the failure mode unconstructible;
    JordanChains is the only Jordan path.  (Test T-15.)
N-3 Rank tests cannot see near-parallel eigenvectors; degenerate snapped
    spectra must not be "diagonalized" (banana fix 3,
    FeynmanTrickBananaStatus.md:243-247; Recurrence.m:138-146).  New
    form: block structure from exact nullity jumps, never from
    eigenvector-matrix rank with a tolerance.  (Test T-7.)
N-4 IntegerQ on approximate numbers is always False — the old resonance
    guard silently let resonant cases through until the box campaign
    fixed it (Docs/FeynmanTrickBoxFamilyStatus.md:109-111).  New form:
    integer tests only ever on exact differences (section 2.6 step 1).
N-5 Gauge degree guard 4d + 8 (ResonantRecurrence.m:179-205): kept as a
    hard internal-consistency error E5.
N-6 Apparent singularities are real charts.  The box L2 beta-root chart
    at t* = 7/11 (eigenvalue −1 recorded by the campaign;
    FeynmanTrickBoxFamilyStatus.md:80-88, memory brief) broke FIVE
    compounding layers of the old particular-solution path precisely
    because "the function is analytic here" was implicitly assumed in
    bookkeeping.  Indicial treats it as any singular chart.  (Test T-17.)
N-7 The worst in-scope structure is known and bounded: banana L1
    endpoints xx1 in {0, 1} have a NILPOTENT DOUBLE POLE with 5-fold
    resonant residues (FeynmanTrickBananaStatus.md:222-227) and a
    genuine size-3 Jordan block (:236-242); the xx1 = 1/2 chart has
    residue eigenvalues {0,0,0,0,1} (:408-417).  In scope from day one
    (RewritePlan 3.2 Indicial; R11 evidence).  (Tests T-15, T-16.)
N-8 Pseudo-resonance is in scope from day one: the banana level-2 upper
    endpoint mixes x^0, x^(−1+eps), x^(2 eps) sectors
    (FeynmanTrickBananaStatus.md:343-346; segment-1 b values 2 and 1
    over a = −1, :81-84).  The family partition must put these in ONE
    joint-solve family.  (Test T-10.)
N-9 The exact full-eps matrix is a PREREQUISITE: slice exports are
    eps-truncated by SeriesCoefficient at export order
    (FeynmanTrick/MatrixExport.m:42-77; loaded per-order via
    DiffExp/MatrixLoading.m:29) and cannot certify I1.  Indicial
    refuses SeriesData input (E1) rather than guessing.

---

## 7. DEPENDENCIES

May call (per the acyclic order `Tolerances < Config < EpsSeries <
SectorSeries < Indicial`; REVIEW-minimalism defect 33, names per
DEC-1):
- `Tolerances.m`: ``Tolerances`DE2Error`` ONLY (sole Tolerances use —
  the I-8 "no numeric tolerances" invariant is about thresholds, not
  the error helper).  NO tolerance accessor is consumed — deliberate
  and load-bearing: this module has no numerics, hence no thresholds
  (invariant I-8).  Stating it here forecloses "just one little
  snapTol" drift.
- `Config.m`: ``Config`PrintInfo`` for verbosity-gated info logging.
  No configuration key changes Indicial's mathematics.
- `EpsSeries.m`, `SectorSeries.m`: NOT USED.  Allowed by the order but
  unnecessary — Indicial works on exact symbols, not truncated arrays.
  Gauge composition (which does need SectorSeries/EpsSeries) executes
  in Solve.m under the section 3.7 contract.

Consumed by: `Solve.m` (primary: `PrepareChart` makes one ChartIndicial
call per chart, caches the result on the chart, and assembles the
ChartSystem itself — DEC-7), `Transport.m` (chart classification for
segmentation/matching decisions; RecombineBasis is tag-driven off
`EpsDegenerateFamilies`, section 2.7), `API.m` (SolveAtPoint surface).

Wolfram built-ins relied on (implementation note): Cancel, Together,
Factor/FactorList (+ Extension), Solve (exact, on univariate factors),
RootReduce, CharacteristicPolynomial, RowReduce/NullSpace WITHOUT
Tolerance options, MatrixPower, Inverse, PolynomialLCM.

---

## 8. UNIT TESTS

File: `Tests/test_indicial.m` (battery-registered at M2).  "Assert"
means exact equality (SameQ after RootReduce/Together canonicalization)
unless stated.  Error tests assert via `Catch[..., "DiffExp2Error"]`
returning the `Failure["DiffExp2", ...]` with the named Module/ID
payload fields (DEC-1; REVIEW-math D17).  Tests T-15..T-17 need
vendored campaign matrices and the single kernel; all others are pure
closed-form and fast.

T-1  test_regular_chart: A = {{0, 1}, {1/(1−t), 0}} (holomorphic at 0).
     Assert: Regular -> True; one family;
     Sectors == ConstantArray[<|"a"->0,"b"->0,"p"->0|>, 2] (DEC-6 /
     REVIEW-math D5: d sector specs, one per basis solution);
     Gauge == IdentityMatrix[2]; PoleOrder == 0.
T-2  test_zero_residue_after_cancellation: A with entries that LOOK
     singular but cancel, e.g. A = {{(2 t)/(2 t^2) − 1/t, 1}, {0, 0}}
     == {{0,1},{0,0}}.  Assert: classified regular (NOT an error, NOT a
     refusal) — kills the old "residue identically zero -> $Failed"
     behavior (Recurrence.m:129-132).
T-3  test_simple_pole_two_families: A = {{(1+eps)/t, 1}, {0,
     (3−2 eps)/(2 t)}}.  Assert: PoleOrder == 1; Spectrum ==
     {(a,b,m)} = {(1,1,1), (3/2,−1,1)}; TWO families (a-difference 1/2
     not an integer), both Class -> "Single", JointSolve -> False, all
     p == 0.
T-4  test_i1_violation_sqrt_eps: A = {{0, 1}, {eps, 0}}/t (RewritePlan
     section 2 counterexample).  Assert: E2 thrown; payload contains
     ID -> "NonAffineEigenvalue", Chart name, the char poly
     lambda^2 − eps, and the irreducible factor.
T-5  test_i1_violation_rational_root: residue {{eps/(1+eps), 0},{0, 0}}
     (A = that /t).  Assert: E2 (root eps/(1+eps) is eps-rational but
     not affine) — proves the check is on eigenvalues, not on entry
     polynomial degree.
T-6  test_affine_pass_nonlinear_entries: residue {{2 eps, 1},
     {−eps^2, 0}} (char poly (lambda − eps)^2).  Assert: NO error;
     Spectrum == {(0, 1, m=2)}; BlockSizes == {2}; sectors
     {(0,1,0), (0,1,1)}; Class -> "Confluent"; LogMax == 1.  (Entries
     quadratic in eps are fine when eigenvalues are affine.)
T-7  test_jordan_vs_diagonalizable: residues diag(eps, eps) vs
     {{eps, 1}, {0, eps}}.  Assert: former gives BlockSizes {1,1},
     sectors {(0,1,0),(0,1,0)}, LogMax 0; latter gives BlockSizes {2},
     sectors {(0,1,0),(0,1,1)}, LogMax 1.  (Lesson N-3: structure from
     exact nullities.)
T-8  test_true_resonance_2f1_resonant: residue of
     Tests/Hypergeometric2F1_Resonant_Matrices/dz_0.m at z = 0
     (a=1/4, b=1/3, c=2; computed by hand from the committed file:
     z*A|_{z->0} == {{0, 0}, {1/12, −2}}).  Assert: eigenvalues
     {0, −2}, both b == 0; ONE family, Class -> "TrueResonant",
     JointSolve -> False; sector specs {(−2,0,1), (0,0,0)} (lower
     eigenvalue gets p = 0 + mult(0) = 1 per the section 2.6 formula);
     LogMax == 1; Collisions == {<|n->2, Type->"Log", DeltaB->0,...|>}.
T-9  test_nonresonant_2f1: residue of
     Tests/Hypergeometric2F1_Matrices/dz_0.m at z = 0 (c = 3/2;
     z*A|_{z->0} == {{0, 0}, {1/12, −3/2}}).  Assert: eigenvalues
     {0, −3/2}; two families; no collisions; all p == 0.
T-10 test_pseudo_resonance_banana_mixture: residue
     diag(0, −1 + eps, 2 eps) (the banana level-2 endpoint sector
     mixture x^0, x^(−1+eps), x^(2 eps);
     FeynmanTrickBananaStatus.md:343-346).  Assert: ONE family
     containing all three; Class -> "Pseudo"; JointSolve -> True;
     all p == 0; LogMax == 0; Collisions contains
     {n->1, DeltaB->1, Type->"LaurentShift"} ((−1,1)->(0,0)),
     {n->1, DeltaB->−1, Type->"LaurentShift"} ((−1,1)->(0,2)), and
     {n->0, DeltaB->−2, Type->"LaurentShift"} ((0,0)->(0,2)).
T-11 test_pseudo_resonance_same_a: residue diag(eps, 2 eps).  Assert:
     one family, JointSolve -> True, collision at n == 0 with
     DeltaB == −eps coefficient ... stated exactly: Collisions ==
     {<|n->0, DeltaB->−1, Type->"LaurentShift",...|>} (b_i − b_j with
     i the (0,1) record, j the (0,2) record); p == 0 everywhere.
T-12 test_rank_reduction_closed_form: A = {{0, t^−2}, {0, 0}}.
     Closed-form solution basis of f' = A f: {(1, 0), (1/t, −1)}.
     Assert: PoleOrder == 2; leading {{0,1},{0,0}} nilpotent (no E3);
     FuchsianReduce terminates with Steps >= 1; reduced ThetaMatrix
     holomorphic; Residue eigenvalues {0, 1} with b == 0; one
     TrueResonant family; AND the composition check: for each reduced
     basis solution g of theta g = B g (closed-form computable since B
     is exactly constant for this input), f = T.g satisfies f' == A.f
     identically in t — verifying section 3.7 with zero numerics.
T-13 test_irregular_loud: A = {{0, 1}, {1, 0}}/t^2.  Leading
     {{0,1},{1,0}} non-nilpotent.  Assert: E3 with PoleOrder == 2 and
     the chart name; no saturation iterations executed.
T-14 test_nontermination_loud: A = {{0, t^−2}, {t^−1, 0}} (nilpotent
     leading matrix but ramified-irregular, Katz rank 1/2: equivalent
     scalar equation t^3 y'' + 2 t^2 y' − y = 0).  Assert: E4 fires
     with Steps == 200 and the chart name.  False success impossible
     (section 2.3 Step 2 note) — the test pins the loop bound.
T-15 test_banana_l1_endpoint_chart [kernel, M2 gate]: load the vendored
     banana level-1 dxx1 full matrix (Tests/refs/, M0 task 14; OQ-1),
     chart xx1 = 0.  Assert: PoleOrder == 2; leading matrix nilpotent;
     reduction terminates; the reduced residue's spectrum is a single
     5-fold integer-spaced b == 0 ... (record: "5-fold resonant
     residues", FeynmanTrickBananaStatus.md:222-227); Jordan data
     contains a size-3 block (:236-242).  The exact (a, b, BlockSizes)
     tuple is FROZEN as a pin on the first kernel run (written to
     Tests/refs/banana_l1_endpoint_indicial.m) and asserted thereafter.
T-16 test_banana_l1_half_chart [kernel, M2 gate]: same matrix, chart
     xx1 = 1/2.  Assert: PoleOrder == 1; eigenvalues {0,0,0,0,1} with
     b == 0 (FeynmanTrickBananaStatus.md:408-417); one TrueResonant
     family containing {0, 1}; Jordan block data frozen as a pin on
     first run (the doc does not record the x = 1/2 block structure;
     OQ-3).
T-17 test_box_l2_apparent_chart [kernel, M2 gate]: vendored box L2
     matrix, chart t* = 7/11 (the beta-root apparent singularity,
     FeynmanTrickBoxFamilyStatus.md:80-88).  Assert: treated as a
     genuine singular chart; the spectrum contains eigenvalue −1
     (campaign record, memory brief / RewritePlan M2); the
     integer-spaced family containing {−1, 0} is TrueResonant; full
     spectrum frozen as a pin on first run.
T-18 test_constructed_spectrum_roundtrip: R = S.diag(1 + 2 eps, −1/2,
     3 eps, 1 + 2 eps).S^−1 with S a fixed unimodular integer matrix
     (R is diagonalizable by construction).  Assert: Spectrum ==
     {(1,2,m=2), (−1/2,0,1), (0,3,1)} exactly; I-3 closure assert
     passes; the (1,2) record has BlockSizes {1,1} (diagonalizable
     repeated eigenvalue — NOT a Jordan block, so it contributes two
     p == 0 sectors); family partition: a-values {1, 0} are
     integer-spaced, so {(1,2), (0,3)} form ONE family with
     Class -> "Pseudo", JointSolve -> True, all p == 0, and exactly
     one collision <|n -> 1, DeltaB -> 1, Type -> "LaurentShift"|>
     (lower (0,3), upper (1,2)); {(−1/2,0)} is its own family with
     Class -> "Single" ((−1/2) minus 0 or 1 is not an integer).  Also
     a determinism check: two calls give identical output (I-9).
T-19 test_inexact_input_loud: A = {{0.5/t, 0}, {0, 0}}.  Assert: E1
     naming entry {1,1}; the word "Rationalize" appears nowhere in the
     handling path (lesson N-1).
T-20 test_series_input_loud: A with a SeriesData entry (an eps-slice
     artifact).  Assert: E1 with the full-export hint (section 5.1,
     lesson N-9).
T-21 test_no_numerics_in_source: read DiffExp2/Indicial.m as text;
     assert no occurrence of `N[`, `Rationalize[`, `Chop[`,
     `Tolerance ->`, `Quiet[`, `Check[` anywhere (invariant I-8 made
     mechanical; the module defines no local error helper).  Calls to
     ``Tolerances`DE2Error`` are ALLOWED — I-8 is about numeric
     thresholds, not the error primitive (REVIEW-minimalism defect 33,
     name per DEC-1).
T-22 test_error_payload_completeness: for each of E1-E7 triggered
     synthetically, assert payload contains ID, Module, Chart,
     Center, Variable (DEC-1), and the class-specific fields of
     section 5.
T-23 test_eps_zero_degeneracy: residues diag(eps, 2 eps) vs
     {{0, 1}, {0, 2 eps}} (the log x class: the residue of Transport.md
     T12's matrix {{0, 1/x}, {0, 2 eps/x}}).  Both are ONE Pseudo
     family (JointSolve -> True, n == 0 collision).  Assert: the former
     has EpsZeroDegeneracy == 0 (eps -> 0 eigenvectors (1,0), (0,1)
     stay independent) and `EpsDegenerateFamilies` returns {}; the
     latter has EpsZeroDegeneracy == 1 (chain tops (1,0) and (1, 2 eps)
     collide at eps = 0) and `EpsDegenerateFamilies` returns
     {<|"FamilyIndex" -> 1, "EpsZeroDegeneracy" -> 1|>}.
     (REVIEW-minimalism defect 14's T-23, adapted to the integer r0
     field of REVIEW-math D8a per DEC-7.)

---

## 9. LINE BUDGET

RewritePlan 3.2 allots Indicial.m ~300 lines (module-map entry:
"Indicial.m (~300) exact residue from the FULL eps-rational matrix;
char-poly factorization contract (I1); Jordan/confluence -> {a,b,p}
specs; resonance/pseudo-resonance partitioning; HIGHER-ORDER POLES:
ported rank-reduction ... with loud non-termination error").  Estimated
spend: pole data ~35; reduction (saturation + adjoin + trim + guards)
~95; spectrum/I1 ~55; Jordan chains ~65; partitioning ~40 (+~10 for the
2.6 step 7 EpsZeroDegeneracy rank, REVIEW-math D8a); the 2.7
`EpsDegenerateFamilies` selector ~4 (DEC-7); orchestration, payload
assembly, asserts ~40 (the error PRIMITIVE is Tolerances.m's DE2Error,
DEC-1 — not budgeted here).  Pre-cut total ~344.

PRE-AUTHORIZED at M0 review (REVIEW-minimalism defect 4): cuts C-2 and
C-3 are TAKEN UP FRONT (−18; already reflected in sections 2 and 5.3),
compensating the D8a/DEC-7 additions (+~14).  RESTATED TOTAL: ~326.
C-1 (trim pass) stays IN until M2 evidence, because T-12/N-1b
regression-pins it.  Per defect 4's library-wide rule: if the
implementation lands >10% over the restated figure, STOP and report to
the orchestrator BEFORE writing more code.

C-1 Drop the trim pass (Step 3, ~25 lines).  Correctness-neutral; cost
    is larger integer a-shifts and wider t/eps windows downstream
    (honest-window arithmetic absorbs it).  Keep the degree guard.
    M2-EVIDENCE-ONLY (defect 4): taking it requires re-pointing T-12's
    N-1b regression role first.
C-2 TAKEN (defect 4): the standalone exports `AffineSpectrum`/
    `JordanChains`/`PartitionResonanceFamilies` are package-private
    (section 2), tested through ChartIndicial's output fields: −10.
C-3 TAKEN (defect 4): E3's witness is restricted to the entry pattern
    of A_{-r}^d, no eigenvalue computation in the error path
    (section 5.3): −8.
C-4 is withdrawn (REVIEW-minimalism defect 35): Solve.m consumes Chains
    verbatim (defect-15 resolution; section 3.4 NORMATIVITY), so chain
    construction cannot leave Indicial.  Remaining headroom comes from
    C-2/C-3 (pre-authorized, defect 4) and, on M2 evidence only, C-1.

What may NEVER be cut: the Moser pre-check, the I1 affine verification
+ closure assert, the exactness asserts, any error path, the collision
list (Solve.m's joint-solve and the 3.4 budget asserts consume it).

---

## 10. OPEN QUESTIONS

(Numbering gaps are deliberate: OQ-2 was CLOSED at the M0 amendment
pass by DEC-1 — the error primitive is ``Tolerances`DE2Error``, section
5.0; REVIEW-math D17, REVIEW-minimalism defect 5.  OQ-5 was DELETED as
subsumed by REVIEW-minimalism defect 15 / DEC-7 — Solve.m's
`PrepareChart` consumes the section 3.4 collision records VERBATIM,
see the NORMATIVITY note there.)

OQ-1 Vendored matrix paths for T-15/T-16/T-17: M0 task 14 vendors pin
     generators and reference data into Tests/refs/, but the banana L1
     dxx1 FULL matrix and the box L2 chart matrix are kernel-generated
     campaign artifacts.  The orchestrator must produce them during the
     M0 idle-kernel oracle pass (RewritePlan M0 KERNEL item) and fix
     the file names; this spec assumes
     Tests/refs/banana_l1_dxx1_full.m and Tests/refs/box_l2_full.m.
OQ-3 The banana xx1 = 1/2 chart's Jordan block structure is not
     recorded in the campaign docs (only the spectrum {0,0,0,0,1} is,
     FeynmanTrickBananaStatus.md:416-417); T-16 freezes it as a pin on
     first kernel run.  If the frozen structure contradicts the size-3
     block expectation IMPLIED nowhere for this chart, no spec change
     is needed — only T-15's endpoint chart carries the size-3 claim.
OQ-4 MaxSteps (200) and MaxPasses (50) are spec constants (old
     defaults, LocalSeries.m:184, 219).  Should they be Config keys?
     Current answer: NO — the kept-config table (Config.md per
     RewritePlan 3.2) does not include them, and a chart needing more
     than 200 saturation steps deserves a human, not a knob.  Revisit
     only with evidence.
OQ-6 Performance of exact Q(eps) linear algebra at d = 5 with deep
     rational entries (banana L1) is untested; R2 flags performance,
     not correctness.  Mitigation if slow: per-chart memoization (one
     ChartIndicial call per chart already), Cancel after every
     elimination step, and evaluation-interpolation for ranks ONLY as
     a verified-exact technique (reconstruct, then verify exactly) —
     any such optimization must keep invariant I-8's spirit: no
     unverified numerics.
OQ-7 Algebraic chart centers (complex singularities, RewritePlan 3.1
     Radius note) put algebraic numbers in EVERY entry; FactorList
     with Extension -> Automatic must cope.  If Mathematica's
     factorization gives trouble over towers (Q(i, Root[...])), the
     fallback is RootReduce-normalized Solve verification (section 2.4
     step 3 path for ALL factors) — slower, still exact.  Flag at M2
     if hit.
