# Adversarial spec review — MATHEMATICAL CORRECTNESS lens (M0 tasks 12-13)

Scope: Docs/specs/{Tolerances,Config,EpsSeries,SectorSeries,Indicial,Solve,
Transport,Integrate,API,FTShimContract}.md against Docs/RewritePlan.md and the
frozen oracle tree.  Every defect below was verified by recomputation or by
reading the cited old code; where a unit test's expected value is claimed
wrong, the correct value was derived independently.  Defects only.

Severity: BLOCKER = spec as written produces wrong mathematics or an
unimplementable/aborting contract on in-scope (campaign) cases.  MAJOR =
contract inconsistency between specs, or a formula/test wrong on reachable
inputs.  MINOR = local formula slip, undefined edge case, or wording that an
implementation agent could misread.

---

## BLOCKERS

### D1 (BLOCKER) — Tolerances.md §2 "LaurentLeadTol", §6, §8 T4; EpsSeries.md §6 N-1/N-2, §8 test 11
**Problem.** `LaurentLeadTol[chopDigits] = 10^(-Floor[chopDigits/2])`
contradicts the campaign calibration both specs claim to preserve.  The
verified in-tree fix is `tolRel = Max[10^(-ChopPrecision/2), 10^(-24)]`
(DiffExp/IntegrationStrategies/Recurrence.m:618, with the code comment:
cancellation residues "~1e-29 relative at WP 300 ... genuine leading content
is either O(1) relative or exactly zero"; Docs/FeynmanTrickBoxFamilyStatus.md
"the 1e-24 floor is load-bearing").  At default chopDigits = 250 the spec's
value is 10^-125: a 1e-29-relative cancellation residue classifies NONZERO —
the opposite of the load-bearing behavior.  Concretely, EpsSeries test 11(a)
is unsatisfiable as written: coefficient `10^-60` against scale `2` is
relative 5·10^-61, which is 54 decades ABOVE the upper ambiguity-band edge
(10^-115) of tol = 10^-125, so `NumericallyZeroQ`/`ESCoeffZeroQ` returns
False and the lead is index 0, not index 1 as the test asserts.  The
Tolerances rationale sentence "real signal at relative >= ~10^-40 empirically
(the campaign's constant)" conflates the FT ABSOLUTE 10^-40 test with a
relative claim and is unsupported; the verified Recurrence.m comment says
genuine leading content is O(1) relative.
**Amendment.** In Tolerances.md §2 replace the LaurentLeadTol entry with:
"`LaurentLeadTol[chopDigits_Integer] -> Max[10^(-Floor[chopDigits/2]), 10^-24]`
RELATIVE leading-coefficient zero test for eps-Laurent window trimming.  The
10^-24 floor is the campaign calibration
(DiffExp/IntegrationStrategies/Recurrence.m:613-618;
Docs/FeynmanTrickBoxFamilyStatus.md:116-124): incomplete cancellations leave
residues up to ~1e-29 relative at WP 300 which MUST classify zero; genuine
leading content is O(1) relative or exactly zero."  In §8 replace T4's
`LaurentLeadTol[250] === 10^-125` with `LaurentLeadTol[250] === 10^-24` and
add `LaurentLeadTol[40] === 10^-20` (the Max picks the derived value below
chopDigits 48).  In EpsSeries.md test 11(b) replace "under
`LaurentLeadTol[250] = 10^-125`" with "under `LaurentLeadTol[250] = 10^-24`".
Delete the Tolerances §2 sentence "real signal at relative >= ~10^-40
empirically (the campaign's constant); `chopDigits/2` (default 125 decades)
separates both by >= 80 decades" and replace with the Recurrence.m-comment
calibration quoted above.

### D2 (BLOCKER) — Tolerances.md §2 "NumericallyZeroQ", §3 consumer table vs EpsSeries.md §2 "ESCoeffZeroQ", §5
**Problem.** Two incompatible zero predicates are both declared THE only one.
Tolerances.md §3: EpsWindow "movement is decided in EpsSeries.m, but ONLY
through `NumericallyZeroQ` with `LaurentLeadTol` — no other predicate may
move a window boundary" (three-valued, aborts in a ±10-decade band).
EpsSeries.md defines and exports `ESCoeffZeroQ` (two-valued, no band, never
errors) and uses it for ESTrim/ESLeading/ESInvert/ERR-DROP-BELOW — every
window-moving decision.  Beyond the naming clash this is mathematically
load-bearing: with the corrected tol of D1 (10^-24) the ±10-decade band is
[10^-34, 10^-14], and the documented 1e-29-relative campaign residues land
INSIDE it — `NumericallyZeroQ` would abort on exactly the inputs the
calibration exists to pass.  The band is appropriate for matching/rank
decisions (where the two populations are separated by >100 decades), not for
cancellation-residue trimming (where the "zero" population extends up to
1e-29 relative).
**Amendment.** In Tolerances.md §2 NumericallyZeroQ and §3, replace the
"ONLY through NumericallyZeroQ" sentence with: "EpsWindow movement is decided
in EpsSeries.m through `EpsSeries`ESCoeffZeroQ` (two-valued, relative,
tolerance `Tol["LaurentLeadTol"]`); `NumericallyZeroQ` with its ambiguity
band is reserved for structure-creating decisions whose zero/nonzero
populations are provably separated by more than 2·$AmbiguityBandDecades:
matching pivots (MatchTol/RankTol) and rank/nullspace cuts.  ESCoeffZeroQ has
no ambiguity band BY DESIGN: cancellation residues (up to ~1e-29 relative,
D1) and genuine O(1)-relative leads cannot be separated by a 20-decade band
around any single threshold."  Update the §3 consumer table row for
EpsSeries accordingly.  Tolerances E5/T9 stand for the matching/rank uses.

### D3 (BLOCKER) — Integrate.md §2.2.3 step 3 (cancellation gate tolerance)
**Problem.** The b = 0 divergence gate drops a coefficient iff
`|c| <= chopFloor · scale_k` while citing "threshold class
`Max[10^(−ChopPrecision/2), 10^−24]`".  `ChopFloor` is the ABSOLUTE noise
floor `10^-chopDigits` (Tolerances.md §2): at the campaign settings (WP 300,
ChopPrecision 250, the FT pin `precision - 50` kept by
FTShimContract.md §8 step 3) the gate threshold is 10^-250, so the documented
~1e-29-relative cancellation residues do NOT drop and E2 fires on the
box/pentagon-class cancellations the gate exists to admit.  Wrong tolerance
name with materially wrong value.
**Amendment.** In §2.2.3 step 3 replace "`|c| ≤ chopFloor · scale_k` ...
`chopFloor` from Tolerances.m" with "`ESCoeffZeroQ[c, scale_k]` at
`Tol["LaurentLeadTol"]` (= Max[10^(-Floor[chopDigits/2]), 10^-24] per
Tolerances.md as amended by D1; the campaign class
Recurrence.m:613-618)".  In §7 change the Tolerances dependency from
"chopFloor (2.2.3)" to "laurentLeadTol (2.2.3, 2.3 step 6)".  Test 10's
values remain valid (at WP 80, tol = Max[10^-20, 10^-24] = 10^-20: 10^-60
drops, 10^-3 errors).

### D4 (BLOCKER) — Config.md §3.2 (AccuracyGoal/AccuracyGoalValidate defaults), §4, §5 E8, §8 C1/C14; Transport.md §2.10/§6.1/Q4
**Problem.** The schema defaults are mutually aborting: "AccuracyGoalValidate"
defaults to "Before" and "AccuracyGoal" defaults to "?", while invariant 4 /
E8 make non-None Validate with non-numeric AccuracyGoal an error.  Therefore
`LoadConfiguration[{}]` (test C1) aborts via the same rule C14 asserts.
Additionally the 3.2 row still describes the OLD adaptive order search
("averaging window ..., +/-10 steps, ...") as the DiffExp2 meaning, while
Transport.md §6.1 explicitly does NOT port it ("AccuracyGoalValidate
'Before' adaptive order search / 'After' redo | NOT ported; replaced by
E11").
**Amendment.** In Config.md §3.2 set `"AccuracyGoalValidate"` Default to
`None`.  Replace its meaning column with: "None: no per-segment accuracy
check; 'Before'/'After': accepted for compatibility, both map to the
DiffExp2 post-hoc per-segment accuracy check (Transport.md E11: segment
probe error must satisfy err <= 10^-DigitsNeeded, error names the suggested
ExpansionOrder); the old adaptive search (old Transport.m:776-841,
1191-1213) is NOT ported — ledger waiver."  C1 then passes; C14 unchanged
(explicitly setting "Before" with AccuracyGoal "?" still aborts).

### D5 (BLOCKER) — Indicial.md §2.1 step 3, §4 I-6, §8 T-1 (regular-chart sector count) vs §4 I-2 and Solve.md §4 I-2
**Problem.** For a regular chart the spec returns "one family, one sector
spec (0,0,0)" (step 3), asserts "exactly one family with exactly one sector
spec" (I-6), and T-1 asserts `Sectors == {<|a->0,b->0,p->0|>}` for a 2x2
system.  This contradicts Indicial's own I-2 ("sum of Length[Sectors] over
Families == d") and Solve.md I-2 (one column per spec, SystemSize columns).
A d-dimensional regular chart has residue 0 with eigenvalue 0 of
multiplicity d (BlockSizes {1,...,1}) and needs d sector specs.
**Amendment.** In §2.1 step 3 replace "(one family, one sector spec
`<|"a"->0,"b"->0,"p"->0|>`, identity gauge)" with "(one family whose single
EigRecord is a = 0, b = 0, Multiplicity = d, BlockSizes = ConstantArray[1,d],
Chains = the standard basis vectors; d sector specs, each
`<|"a"->0,"b"->0,"p"->0|>`; identity gauge)".  In I-6 replace "exactly one
sector spec (0, 0, 0)" with "exactly d sector specs, all (0, 0, 0), one
family".  In T-1 replace `Sectors == {<|a->0,b->0,p->0|>}` with
`Sectors == ConstantArray[<|"a"->0,"b"->0,"p"->0|>, 2]`.

### D6 (BLOCKER) — Solve.md §2/§3.2 vs Indicial.md §2.6/§3.6 vs Transport.md §2.6 step 1: the ChartSystem contract has no producer and a nonexistent entry point
**Problem.** (i) Transport.md §2.6 step 1 and §7 call
``DiffExp2`Solve`SolveChart[sys, chart, ...]``; Solve.md §2 exports "Exactly
three" symbols and SolveChart is not one of them.  Transport's `chart` (its
§3.2: geometry record, no matrix) is also a different object from Solve's
`ChartSystem` (§3.2: ThetaMatrix, Gauge, Families, Residue).  (ii) Solve's
ChartSystem requires per family `"V"`, `"VInv"`, `"J"`, `"CollisionDepth"`
and the certified det-V collision-factor property; Indicial's IndicialData
provides none of these (it emits per-eigenvalue `Chains`, no assembled V, no
inverse, no CollisionDepth, no det certification) — acknowledged in Solve
OQ1/Indicial OQ-5 but unresolved, and the M0 gate requires these specs to be
implementable as written.  (iii) Collision-record cardinality conflicts: for
n = 0 collisions Indicial §2.6 step 5 records "ONLY the canonically ordered
one" (T-11 asserts exactly one record for diag(eps, 2eps)) while Solve §3.2
demands "all ordered pairs with integer offset >= 0" (two records at n = 0).
(iv) Nobody owns applying the chart map to the loaded system, forming
theta-form B = t·A in chart coordinates, and calling ChartIndicial.
**Amendment.** (a) Add to Solve.md §2 a fourth export:
"`PrepareChart[sys_Association, chart_Association] -> ChartSystem` — applies
chart["Map"] to the loaded exact matrix, forms the theta matrix, calls
``Indicial`ChartIndicial``, and assembles the ChartSystem of 3.2: per family,
V = the d x (family size) columns extended to the full d x d matrix by
concatenating ALL families' chain vectors in family order (column order =
(family, root, chain position), eigenvector first); VInv = exact
`Inverse[V]`; J = the corresponding block-diagonal Jordan data;
CollisionDepth = the longest directed path in the family's collision DAG
restricted to Type == 'LaurentShift' edges; the det-V collision-factor
property is certified here by exact factorization of Det[V(eps)] (E2 on
failure).  +~40 lines against the §9 budget (OQ1 reserved +20; take the
remaining 20 from cut item 1)."  (b) Transport.md §2.6 step 1 becomes:
"`cs = DiffExp2`Solve`PrepareChart[sys, chart]`; `basis =
SolveHomogeneous[cs, req]`; `particular = SolveParticular[cs, source, req]`
when a source is present; metadata incl. CouplingDepth from
basis["Diagnostics"]."  (c) Indicial.md §2.6 step 5: delete the
n = 0 single-record rule and record BOTH ordered pairs (amend T-11 to assert
two records, `DeltaB -> -1` and `DeltaB -> +1`), OR Solve.md §3.2 must state
that n = 0 collisions arrive canonically ordered and are symmetrized by
PrepareChart — pick one and apply to both files.

### D7 (BLOCKER) — Transport.md §2.7 (MatchWeights column normalization eps^(-kmin_i))
**Problem.** "columns pre-normalized by eps^(-kmin_i) so every column's
window starts at eps^0 — note homogeneous true-resonance sectors have
kmin = -p BY CONSTRUCTION" conflates the COEFFICIENT-storage window
(kmin = −p) with the VALUE window of the evaluated column.  Per
SectorSeries.md §2.4.1 the evaluated value's MinPower is
min over sectors of (first present row + p): a true-resonance column with
storage kmin = −p has value MinPower = 0 ALREADY (the (eps Logx)^p/p!
normalization is exactly what makes log-chain values start at eps^0 —
RewritePlan 3.1 exemption).  Multiplying such a column's values by
eps^(-kmin_i) = eps^(+p) shifts genuinely O(eps^0) data to order p and makes
the graded solve wrong on EVERY chart with a log (2F1 resonant, banana, box
L2) — a wrong-mathematics path, not a style issue.  For compensated
pseudo-resonant columns the value window can start below 0 only with
numerically-cancelling content (Solve.md §3.7 step 3 / I-5), which must be
TRIMMED, not shifted.
**Amendment.** Replace the parenthetical in §2.7 with: "columns are the
EpsSeries VALUES of the basis solutions at the match point; before the
graded solve, each column is trimmed (`ESTrim`) and it is ASSERTED that the
trimmed window starts at eps^0 — negative-order content must be numerically
zero (it is the cancelling polar content of compensated joint columns,
Solve.md I-5); a negative-order coefficient failing `ESCoeffZeroQ` at
matchTol is error E5, gray zone is E7.  No eps-power renormalization of
columns is performed: true-resonance columns already have value MinPower = 0
under the (eps Logx)^p/p! normalization."  Align §3.5 ("Min -> 0") wording:
the assertion applies to ALL columns post-trim, not only "regular incoming
data".

### D8 (BLOCKER) — Transport.md §2.8/Q2 vs Indicial.md §2.5 vs Solve.md OQ4: "EpsDegenerateFamilies" is produced by no module
**Problem.** RecombineBasis is driven by `chart["EpsDegenerateFamilies"]`
"from Indicial.m metadata"; Indicial.md §2.5 explicitly declines to compute
it ("The eps -> 0 collision of EIGENVECTORS ... is NOT detected here"),
pointing at the family flags; Solve.md OQ4 claims recombination is
unnecessary ("the compensated basis already has ord_eps det = 0").  The
Solve claim is FALSE for the canonical in-scope class: for
A = {{0, 1/x}, {0, 2 eps/x}} (Transport T12) the two columns are exact
single-sector solutions — no Laurent division ever fires in the recursion,
no compensation is registered, and the eigenvectors (1,0) and (1,2eps)
collide at eps = 0, so det F(0) = 0 at every match point.  As specced, the
metadata is absent, RecombineBasis never fires, and E5 aborts on a valid
day-one system (banana endpoint data is the same class).
**Amendment.** (a) Indicial.md: add to §2.6 a step 7 and to §3.4 a field:
"`"EpsZeroDegeneracy" -> r0_Integer` per family: assemble the family's
chain-top vectors v_i(eps) (columns), normalize each by its eps-valuation
(exact), substitute eps -> 0 (exact), and set r0 = (number of columns) −
rank over Q(alpha) of the resulting matrix.  r0 > 0 iff the eps -> 0
eigenvectors collide; r0 is the recursion depth RecombineBasis needs.
Computed only for families with >= 2 distinct b (elsewhere r0 := 0).  Delete
the sentence claiming this is 'NOT detected here'; the division of labor is:
Indicial DETECTS (exact rank at eps = 0), Transport REMOVES (unimodular
recombination)."  (b) Transport.md §2.8: "EpsDegenerateFamilies" :=
families with `EpsZeroDegeneracy > 0`; recursion depth := r0.  (c) Solve.md
OQ4: replace the current answer with "RAW columns are not needed, but
recombination IS: joint compensation only fires when a Laurent division
occurs in the recursion; eps=0 eigenvector collisions without any division
(the log x class) leave det F(0) = 0 and are Transport's RecombineBasis job,
keyed on Indicial's EpsZeroDegeneracy."

---

## MAJORS

### D9 (MAJOR) — Transport.md §8 T12: the incoming data is not a solution of the test system
**Problem.** For A = {{0, 1/x}, {0, 2 eps/x}}, f2 must satisfy
f2' = 2 eps f2/x, so f2 = C·x^(2eps); the test's incoming vector
`((x^(2eps)-1)/(2 eps), x^(2eps)-1)` has second component x^(2eps) − 1,
which fails the ODE (its derivative is 2 eps x^(2eps-1), not
2 eps (x^(2eps)-1)/x).  With weights (0,1) the basis reproduces
B2 = ((x^(2eps)-1)/(2eps), x^(2eps)); as written the residual assert E6
fires and the test cannot pass.
**Amendment.** In T12 replace "incoming data = the regular solution
((x^(2eps)-1)/(2 eps), x^(2eps)-1)" with "incoming data = the regular
solution ((x^(2eps)-1)/(2 eps), x^(2eps))".

### D10 (MAJOR) — Solve.md §8 SU-03 and API.md §8 test 12: wrong 2F1 sector exponents for the committed matrices
**Problem.** Tests/Hypergeometric2F1_Matrices/dz_0.m (verified) is
A = {{0, 1}, {1/(12 z (1-z)), (18-19z)/(12(z-1)z)}}, so
z·A|_{z->0} = {{0,0},{1/12, -3/2}} with eigenvalues {0, −c} = {0, −3/2}
(Indicial.md T-9 states this correctly).  The scalar indicial roots
{0, 1−c} = {0, −1/2} are roots of the SCALAR equation; the (y, y')
companion VECTOR solutions carry tags {0, −3/2} (the y' component
dominates: y ~ z^(1-c) gives the vector ~ z^(-c)).  Solve SU-03 asserts
"roots {0, 1−c(eps)}" and API test 12 asserts "sector exponents EXACTLY
{0, -1/2} (indicial roots of 2F1 at z=0: 0 and 1 - c = -1/2)" — both wrong
for these matrices and inconsistent with Indicial T-9.
**Amendment.** Solve.md SU-03: replace "roots {0, 1−c(eps)}" with "residue
eigenvalues {0, −c} (companion-system tags; the scalar exponent 1−c appears
as component-1 content at column n = 1 of the a = −c sector)".  API.md test
12: replace "sector exponents EXACTLY {0, -1/2} (indicial roots ... 1 - c =
-1/2)" with "sector exponents EXACTLY {0, -3/2} (companion residue
eigenvalues {0, −c}; the scalar root 1−c = −1/2 is the a = −3/2 sector's
first component at n = 1), b == 0 for both".

### D11 (MAJOR) — SectorSeries.md §2.6 (re-expansion tail estimate): formula valid only in a unit-radius coordinate
**Problem.** `TailEstimate[k] = Max|c'[k,N,comp]| · q^(N+1)/(1−q)` with
q = 1/DivisionOrder.  The output coefficients c' are in the t' = t − D
coordinate whose radius is rho, and the evaluation point is r =
rho/DivisionOrder.  From the geometric model |c'_n| ~ A·rho^(−n) the true
tail at r is Σ_{n>N} |c'_n| r^n = |c'_N| · r^N · q/(1−q).  The spec's
expression equals this only when rho = 1; for rho > 1 it underestimates by
rho^N and for rho < 1 (small charts — common) it overestimates by rho^(−N),
making ErrorEstimate aborts/passes wrong by many orders of magnitude.
**Amendment.** Replace the formula with:
"`TailEstimate[k] = Max over comp of Abs[c'[k, N, comp]] * r^N * q/(1-q)`,
with `r = rho/DivisionOrder` (the design evaluation radius) and
`q = 1/DivisionOrder` (so q/(1−q) = 1/(DivisionOrder−1)).  Equivalently, in
the radius-normalized coordinate w = t'/rho the familiar form
|ĉ_N| q^(N+1)/(1−q) holds with ĉ_N = c'_N rho^N; the implementation must use
whichever form matches its stored coordinate."

### D12 (MAJOR) — Integrate.md §2.2.3/§4 I2 vs SectorSeries.md §2.2/§9 cut 4: the cancellation gate needs the integer-spaced merge, which I2 does not require and a sanctioned cut would remove
**Problem.** The gate sums "the offending coefficients across the master
combination" by relying on same-tag merging, and I2 only requires pairwise
distinct (a,b,p).  Distinct-tag sectors (−1,0,0) with c[0,0] = 1 and
(−2,0,0) with c[0,1] = −1 both carry the monomial t^(−1); they satisfy I2,
the gate inspects them separately, and E2 fires although the merged t^(−1)
coefficient is exactly 0.  Soundness requires the canonical form INCLUDING
the integer-spaced same-(b,p) merge (SectorSeries §2.2 bullet 2) — yet
SectorSeries §9 cut item 4 offers to cut exactly that merge, which would
silently break this gate.
**Amendment.** In Integrate.md §2.2.3 precondition add: "and integer-spaced
same-(b,p) sectors merged to minimal a (SectorSeries.md 2.2; assert at
entry: no two sectors share b, p with integer a-difference — violation is
E7-class).  The gate's per-monomial coefficients are well-defined only on
this canonical form."  In SectorSeries.md §9 cut item 4 append: "FORBIDDEN
for objects consumed by Integrate.m — Integrate 2.2.3's cancellation gate
and EndpointLimit's t^0 readout (API.md 2.5) require the merge; this cut may
only be taken together with moving the merge into Integrate.m's entry
validation."

### D13 (MAJOR) — API.md §2.5 (EndpointLimit): "constant of the (0,0,0)-sector" is ill-defined on canonical (merged) objects
**Problem.** After the integer-spaced merge a LocalSolution may contain NO
(0,0,0) sector even when the limit is finite and nonzero: a merged
(−1,0,0) sector holds the former (0,0,0) content at column n = 1.  As
written EndpointLimit returns nothing/zero for such inputs (the same
wholesale-drop disease class the spec cites).  RewritePlan 3.3's shorthand
is the source; the API spec must be the precise one.
**Amendment.** In §2.5 replace "the limit is the eps-Laurent constant of the
(a,b,p) = (0,0,0) sector" with: "the limit is the coefficient of the
monomial t^0: Σ over sectors with b == 0, p == 0 and integer a <= 0 of
c[k, n = −a, comp] (exact tag selection; on canonical merged input at most
one such sector exists per component).  Sectors with b != 0 are dropped
exactly (dimreg).  Divergence gate first: any b = 0 sector contributing a
monomial t^m with m < 0 (any p), or t^0 with p > 0, whose merged coefficient
fails the cancellation test (Integrate.md 2.2.3 semantics, laurentLeadTol)
is the loud error E26 naming (component, tag, eps order)."

### D14 (MAJOR) — Transport.md §5 E8 vs API.md §5 E19 vs Config.md §3.2: three different scopes for AbortOnAnalyticContinuationFail -> False
**Problem.** Transport E8 permits the degraded (honest-partial) result ONLY
"iff ... AND the chart is the FINAL one in the plan"; API E19 degrades to
warning + "MultivaluedFail" -> True with no positional restriction; Config's
key row says "False: downgrade to warning + flagged, FT pipeline mode" with
no restriction.  An implementation agent cannot satisfy all three.
**Amendment.** Adopt Transport's rule as normative (it is the one that makes
the pentagon configuration gap loud mid-path).  In API.md E19 append: "The
downgrade applies ONLY when the offending chart is the final chart of the
plan (Transport.md E8); a mid-path unprescribed/conflicting branch point
aborts regardless of this key."  In Config.md 3.2, AbortOnAnalyticContinuationFail
meaning column, append: "(False softens the FINAL chart only; mid-path
failures always abort — Transport.md E8)".

### D15 (MAJOR) — Tolerances.md §2/§8 T4 vs Transport.md §2.7/§5 E7: SnapTol === RankTol makes the pivot gray zone empty
**Problem.** Transport's three-valued pivot logic ("|pivot| < snapTol*scale
is snapped to exact 0; snapTol*scale <= |pivot| <= rankTol*scale is the gray
zone = LOUD ERROR") requires snapTol < rankTol, but Tolerances derives both
as 10^(-Floor[wp/2]) and T4 asserts `SnapTol[500] === RankTol[500]` — the
gray zone degenerates to a point and E7/L10's three-valued honesty is
vacuous.
**Amendment.** In Tolerances.md §2 redefine
`RankTol[wp_Integer] -> 10^(-Floor[wp/4])` with the note "rank/pivot cuts
must sit a band above the snap floor; the gray zone
[SnapTol*scale, RankTol*scale] (width wp/4 decades) is the loud-error
region of Transport.md E7".  Update T4 to `SnapTol[500] === 10^-250,
RankTol[500] === 10^-125` and the §6 absorption note for the fitter's
`10^(-prec/2)` cut ("tightened to wp/4 to open the gray zone; the fitter
itself is deleted").  Re-check Transport T15's geometric-mean construction
still lands strictly inside the band (it does, by construction).

### D16 (MAJOR) — Tolerances.md §3 (binding consumer table) contradicts four module specs and omits a needed export
**Problem.** The table is declared binding ("anything else is a spec
violation found in review"), and: (i) Transport.md §7 consumes `snapTol`
(pivot snap 2.7, E14, I2 geometry) and `chopReserve` (E3) — neither granted;
`chopReserve` is not even an exported Tolerances symbol.  (ii) EpsSeries.md
uses `matchTol` (ESSameQ; its §7 says so) — not granted.  (iii) Integrate.md
§7 uses `matchTol` (I6 additivity check) — not granted.  (iv) The table's
EpsSeries row grants "ChopFloor (coefficient chop after add/mul/div)", which
directly contradicts EpsSeries I-6 ("the module never applies N,
SetPrecision, or Chop to stored coefficients").
**Amendment.** In Tolerances.md §2 add
"`ChopReserve[wp_Integer, chopDigits_Integer] -> wp − chopDigits` (digit
reserve between working precision and the chop floor; consumed by
Transport.m DigitBudget E3; old reserve WP − ChopPrecision = 250,
State.m:109,135)."  Rewrite the §3 consumer table rows: "EpsSeries.m:
LaurentLeadTol (ESCoeffZeroQ window decisions), MatchTol (ESSameQ) — NO
ChopFloor: EpsSeries never chops (its I-6)"; "Transport.m: MatchTol,
RankTol, SnapTol, ResidTol, ChopReserve, $SafetyDigits,
$InputPrecisionFactor, expansion-order constants, EvalErrorSeriesDecrease";
"Integrate.m: LaurentLeadTol (cancellation gate + result trim), MatchTol
(I6), SnapTol (tiling/endpoint snapping), ChopFloor (pre-Pade-free numeric
hygiene only if used; else drop)".

### D17 (MAJOR) — error-mechanism layering: Tolerances.md §5 vs Indicial.md §5.0 vs EpsSeries §5 vs SectorSeries §5.1 vs Solve §5 vs Integrate §5 vs API §5
**Problem.** Tolerances.md says errors go through "the shared DiffExp2
loud-error primitive (... the API.m spec owns its exact form)" — but API.m
is the TOP of the acyclic order; the seven lower modules cannot call an
API-owned helper.  Meanwhile each spec pins a DIFFERENT mechanism:
print+Abort (Tolerances, EpsSeries), Message+Throw (SectorSeries),
Throw[Failure["DiffExp2", payload]] (Indicial), "single error head
DiffExp2::solve" (Solve), "structured Failure" (Integrate).  Unit tests are
written against the mechanism (CheckAbort vs Catch vs message assertions),
so this is not cosmetic: the M1-M3 test suites as specced cannot all pass
against one library.
**Amendment.** Add to Tolerances.md §2 (the bottom module): "`DE2Error[
module_String, tag_String, payload_Association]` — prints one line
'DiffExp2 <module> error <tag>: <payload summary>' and
`Throw[Failure["DiffExp2", Join[<|"Module"->module,"Error"->tag|>, payload]],
"DiffExp2"]`.  THE library-wide loud-error primitive; every module raises
through it; no module defines its own Abort/Message idiom."  In each other
spec's §5, replace the mechanism sentence with a reference to
`Tolerances`DE2Error` (payload fields unchanged), and in each §8 state that
error tests assert via `Catch[..., "DiffExp2"]` returning the Failure with
the named Module/Error fields.  Delete Tolerances' "API.m spec owns its
exact form" clause; resolve EpsSeries O-1, Indicial OQ-2 accordingly.

### D18 (MAJOR) — eps-Laurent value-object shape: four incompatible field sets across EpsSeries.md §3.1, SectorSeries.md §3.4, Integrate.md §3.4, API.md §3.8
**Problem.** The owner (EpsSeries) defines
`<|"EpsWindow" -> <|"Min","CompleteMax"|>, "Coeffs"|>`; SectorSeries'
EpsLaurentValue uses `"MinPower"/"CompleteMax"/"Coefficients"`; Integrate's
EpsLaurent uses `"Min"/"CompleteMax"/"Coefficients"`; API's LaurentValue
uses `"MinPower"/"Coefficients"` (+"CompleteMax").  All three consumers
defer via open questions, but the M0 gate is spec freeze; an implementation
agent reading any single pair produces incompatible code.
**Amendment.** Pin: the EpsSeries object (`"EpsWindow"`+`"Coeffs"`) is the
ONLY in-core shape; SectorSeries §3.4 and Integrate §3.4 are rewritten to
return EpsSeries objects with metadata keys ("PadeFallbacks",
"TailEstimates") carried in a sibling association, not inside the series;
API §3.8 LaurentValue is the FT-boundary shape only (`"MinPower"`,
`"Coefficients"`, `"CompleteMax"`), produced exclusively by API.m via
`ESMinPower`/`ESCoefficientList`/`ESCompleteMax`.  Mark Q1 in
SectorSeries/Integrate and API OQ as resolved by this rule.

### D19 (MAJOR) — Tolerances.md §2 SnapTol + Transport.md §5 E14/W5 + API.md §2.3: wp-derived SnapTol makes user-input target snapping dead and creates a degenerate-geometry failure mode
**Problem.** SnapTol[wp] = 10^(-Floor[wp/2]) (10^-250 at defaults).  A
user/FT numeric target near an exact singularity (machine input, error
~1e-16; e.g. `to = 0.3333333333333333` for 1/3) is neither equal nor within
snapTol, so E14 does not fire and W5 never snaps; the point is treated as
regular at distance ~1e-16 from a singularity — chart radius 1e-16,
DigitBudget/E3 blow-up or a nonsensical plan.  The old absolute 10^-10
(State.m:123) caught this class.  The wp/2 derivation is right for
INTERNALLY computed values (full precision) but wrong for external inputs
whose error is set by input resolution.
**Amendment.** In Transport.md E14/W5 and API.md W5 replace "within
snapTol" with "within `InputSnapTol[v]`, where `InputSnapTol[v_] :=
10^(-Floor[3*Min[Precision[v], Accuracy[v]]/4])` for inexact external input
v (a new Tolerances.m export; exact input uses exact membership only)".
Add the `InputSnapTol` entry to Tolerances.md §2 with the note: "external
numeric inputs carry input-resolution error; snapping them at the internal
SnapTol (wp/2 decades) is a dead test — the lesson behind the old absolute
10^-10, made input-aware instead of absolute."

### D20 (MAJOR) — API.md §2.3/§2.6 vs FTShimContract.md §7 G1 / Transport.md §2.1: "ExtraSingularFactors" is consumed but never accepted
**Problem.** Transport.md 2.1 reads `sys["ExtraSingularFactors"]` "(config;
the FT layer's IBP-coefficient denominators ...)", and the shim contract
(G1) requires a TransportTo/IntegrateOverLine option, flagged as an API gap
— but API.md §2 defines no such option and LoadSystem's SystemInfo has no
such field.  Without it, the segmentation alphabet excludes IBP-coefficient
poles, making SectorSeries `::interiorpole` reachable on FT lines (the
"Combination" rational multiply is not closed) — an algebra-closure hole,
not just plumbing.
**Amendment.** In API.md §2.3 and §2.6 add the option:
"`"ExtraSingularFactors" -> {poly..}` — exact polynomials in the pinned
variables, unioned EXACTLY into the segmentation alphabet for this call
(charts are placed at every real root in range; complex roots cap radii).
Loud error if any "Combination"/"Prefactor" denominator has a root not in
the alphabet (closure check; the implicit-widening alternative is rejected
per FTShimContract G1)."  Mirror one sentence in §3.4 SystemInfo ("the
per-call extra alphabet is NOT system state — it is an option").

---

## MINORS

### D21 (MINOR) — Solve.md §3.6 CASE R: missing d_0 factor
**Problem.** "the block-i rows of level ℓ read −N_i·û[n,ℓ]_i +
d_0·eps·û[n,ℓ+1]_i = R[n,ℓ]_i" drops d_0 on the first term; from
L_n = d_0(δ I − J) with δ ≡ 0 the correct row set is
−d_0·N_i·û[n,ℓ]_i + d_0·eps·û[n,ℓ+1]_i = R[n,ℓ]_i.  (The eigendirection
sentence "d_0·eps·û[n,ℓ+1]_eig = R_eig" is consistent only with the
corrected form.)
**Amendment.** Insert `d_0·` before `N_i·û[n,ℓ]_i`.

### D22 (MINOR) — Integrate.md §2.2.2: `kMaxOut = κmax + shift_s(k)` with `shift_s(k)` never defined
**Problem.** The per-monomial truncation argument is undefined, and a
literal reading (per-k CompleteMax κmax + shift before the eps^k multiply)
gives k-dependent windows inconsistent with the static account of §2.1.3.
**Amendment.** Replace `κmax + shift_s(k)` with: "`kMaxOut = K_s − k`, where
`K_s = κmax + p` for sectors that cannot hit the pole cell and
`K_s = κmax − 1` for pole-hit sectors (b ≠ 0 with a + n + 1 = 0 reachable),
so that after the eps^k shift every contribution shares the sector's
uniform CompleteMax K_s of 2.1.3."

### D23 (MINOR) — EpsSeries.md §2 ESFromExpression: leading exponent above kmax undefined
**Problem.** For input with exact leading eps-power > kmax (e.g.
`ESFromExpression[eps^7, eps, 3]`), `Series` returns a termless SeriesData
and "Min = the exact leading exponent reported by SeriesData" is undefined;
the spec covers only eps-free and zero inputs.
**Amendment.** Add to the ESFromExpression contract: "If `Series` returns a
SeriesData with no terms through kmax (leading exponent above the window),
the result is `ESZero[kmax]` — completeness through kmax and the below-Min
certificate both hold."

### D24 (MINOR) — EpsSeries.md §2: ESCoeffZeroQ signature has no tolerance argument but ESSameQ calls it "at matchTol"
**Problem.** `ESCoeffZeroQ[c_, scale_]` is hardwired to laurentLeadTol, yet
ESSameQ's contract says coefficient differences "must satisfy ESCoeffZeroQ
at matchTol" — uncallable as signed.
**Amendment.** Change the export to
`ESCoeffZeroQ[c_, scale_, tol_:Automatic]` with "Automatic ->
Tol[\"LaurentLeadTol\"]; ESSameQ passes Tol[\"MatchTol\"] explicitly".

### D25 (MINOR) — Transport.md §2.5: DigitBudget with the default AccuracyGoal "?" is undefined
**Problem.** `DigitBudget[accuracyGoal, segmentCount]` and the plan's
"DigitsNeeded" assume numeric AccuracyGoal; the Config default is "?".  E11
is conditioned on numeric AG but E3 and the plan field are not; the old code
silently skipped budgeting (Transport.m:776-777 class).
**Amendment.** Add to §2.5: "With AccuracyGoal == \"?\": DigitsNeeded :=
ChopDigits (the chop reserve is then the only budget), E3 still applies, and
E11 is skipped — recorded in the plan as `"BudgetMode" -> "Unvalidated"`,
never silently."

### D26 (MINOR) — SectorSeries.md §2.5 step 3 / §7: the `::geomambiguous` guard has no named tolerance
**Problem.** The pole-modulus-vs-Radius comparison errors "when the
difference is below the Tolerances guard", but no Tolerances export is named
and Tolerances §3 expects SectorSeries to use none beyond
ChopFloor/ResidTol.
**Amendment.** In both files name it: SectorSeries uses
`Tol["RankTol"]·Max[|t_i|, Radius]` as the ambiguity width for the
|t_i| vs Radius (and |D| vs Radius) comparisons; add "RankTol (geometry
ambiguity guard, ::geomambiguous)" to the SectorSeries row of Tolerances §3.

### D27 (MINOR) — Solve.md §5 E9: the "would integrate t^{-1}" clause describes a nonexistent code path
**Problem.** Solve never integrates in t; a source tag hitting δ = 0 is
CASE R (log bump), handled exactly, and b = 0 endpoint divergences belong
wholly to Integrate.m's table.  The clause "a source sector with a ≤ −1
integer and b = 0 colliding with nothing is flagged when its particular
would integrate t^{-1}" can only mislead an implementer into adding a
spurious divergence error inside Solve.
**Amendment.** Delete the clause; E9 ends at "...source in d/dt form
detected (caller contract breach is not detectable in general; the theta
normalization of 3.3 is the caller's obligation)".

### D28 (MINOR) — Transport.md §2.9(ii)/§7 vs Integrate.md §7: contradictory ownership of the crossing/branch rule
**Problem.** Transport declares "CrossingOperator is the single exported
authority that Integrate.m reuses" and its §7 allows Integrate to call it;
Integrate's §7 forbids calling Transport.m and instead binds its negative-arm
evaluation to SectorSeries' shared branch rule ("must be the same function
or a re-export").  Both cannot hold.
**Amendment.** Single owner: SectorSeries owns the scalar branch rule
(sigma substitution, its 2.4.1); Transport.md 2.9(ii) becomes "the
sector-level mixing operator is Transport-internal and is DEFINED as the
SectorSeries sigma rule applied tagwise; Integrate.m's interior pairing and
negative-arm boundary terms use the SectorSeries rule directly (its 2.2.5),
guaranteeing the conventions coincide"; delete "(Integrate may call
Transport's exported CrossingOperator)" from Transport §7.

### D29 (MINOR) — Tolerances.md §8 T9: "exclusive" band wording contradicts the stated comparisons
**Problem.** With `True iff r < tol/10^band` and `False iff r > tol*10^band`,
the endpoints 10^-135/10^-115 themselves fall into the ELSE branch (loud
error), i.e. the error band is CLOSED; T9's parenthetical calls it
"exclusive".  An implementer pinning boundary behavior to the comment gets
the opposite of the comparison spec.
**Amendment.** Replace "(closed-form: band is [tol/10^10, tol*10^10] =
[10^-135, 10^-115], exclusive)" with "(closed-form: the ERROR band is the
closed interval [tol/10^10, tol*10^10] = [10^-135, 10^-115]; the True/False
regions are the open exteriors)".
