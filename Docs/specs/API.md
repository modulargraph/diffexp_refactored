# Module spec: DiffExp2/API.m

> **PENDING AMENDMENTS (apply before M5 consumption; tracked, binding
> resolutions in DECISIONS-M0.md):** REVIEW-math D10 (test 12 exponents,
> DEC-17), D13 (EndpointLimit on merged objects), D14 (abort scope E19,
> DEC-15), D18 (value shape keys, DEC-13); REVIEW-minimalism 3 (PrepareChart/
> EpsDegenerateFamilies naming, DEC-7), 5 (DE2Error, DEC-1), 7 (Prefactor
> narrowing + test 17, DEC-22), 12 (EndpointLimit owner line); DEC-24
> (IndefiniteIntegral note). Already applied: D19/D20, min-8/9/11/13,
> DEC-10/11/14/23.


Status: M0 deliverable (RewritePlan section 3.2, API.m entry; section 6 M0
task 4-11).  This document is the FULL user-facing compatibility contract for
DiffExp2.  It is written to be implementable by an agent that reads ONLY this
spec, Docs/RewritePlan.md, and the old code cited here by file:line.  All old
code paths are relative to the repo root /Users/mhidding/Code/diffexp_refactored.
"Old Transport.m" means DiffExp/Transport.m (the frozen oracle), never the new
DiffExp2/Transport.m.

---

## 1. PURPOSE

API.m is the thin user-facing surface of DiffExp2: it ingests external inputs
(matrix directories, boundary conditions, lines/points, prefactor specs),
normalizes and validates them ONCE at the boundary (precision raising, eps
symbol aliasing, wildcard/closed-form/Laurent boundary-condition parsing), and
exposes seven entry points — LoadSystem, PrepareBoundaryConditions,
TransportTo, SolveAtPoint, EndpointLimit, IntegrateOverLine, ToPiecewise —
that delegate all mathematics to the lower modules (Transport.m, Integrate.m,
Solve.m, SectorSeries.m via the dependency order in section 7).  It owns the
result-object format (the association every test and the FT layer consumes:
"SeriesValues", "ErrorEstimates", "KinematicPoint", "SegmentData",
"EndpointIsSingularity", plus the new honest-window keys), and it is the
compatibility layer that mirrors the old DiffExp` entry points
(PrepareBoundaryConditions / IntegrateSystem / TransportTo / ToPiecewise,
old Transport.m:23-26) closely enough that the classic examples and tests
migrate with mechanical edits only, while making the FT-pipeline operations
(solve at singular points, endpoint limits, regularized integration over
[0,1]) first-class instead of reach-ins.

---

## 2. PUBLIC SYMBOLS

All symbols live in context DiffExp2` and are the ONLY exported symbols of
the library besides the Config.m configuration functions (LoadConfiguration,
UpdateConfiguration, CurrentConfiguration — specified in Docs/specs/Config.md;
API.m re-exposes them in DiffExp2` but their implementation and line budget
belong to Config.m).

Throughout: `x` denotes the configured LineParameter symbol; `eps` denotes
the canonical DiffExp2 regulator symbol.  EVERY input expression (matrix
entries, boundary conditions, prefactors, combinations) may use either
Global`eps or \[Epsilon]; both are normalized to the canonical symbol at
ingestion (old precedent: MatrixLoading.m:118 accepts both in closed-form
matrices).  Output never contains the regulator symbol — only Laurent
coefficient arrays with explicit windows.

### 2.1 LoadSystem

    LoadSystem[dir_String, opts___Rule]  ->  SystemInfo (Association, sec. 3.4)

Loads the DE system from a matrix directory, replacing any previously loaded
system wholesale (exactly-one-system invariant, section 4).  Accepted
directory contents, probed in this order:

  (a) EXACT FULL FORMAT (preferred; REQUIRED for any FT / Indicial /
      singular-chart use): one file `d<var>_full.m` per kinematic variable
      <var>, each a square matrix of rational functions of (<vars>, eps)
      — the ExportGeneralMatrix format (FeynmanTrick/MatrixExport.m:91-114).
      Mode "Exact".
  (b) CLOSED FORM: one file `d<var>_d.m` per variable, same mathematical
      content class as (a) (exact symbolic eps), the old closed-form path
      (MatrixLoading.m:42-46, 113-125).  Mode "Exact".
  (c) LEGACY SLICES: files `d<var>_<k>.m`, k = 0..kmax, per-eps-order
      Taylor slices (MatrixLoading.m:29, 91-111; written by
      ExportDiffExpMatrix, FeynmanTrick/MatrixExport.m:42-83).  Mode
      "LegacySlices".  Accepted ONLY for classic regular-point parity
      transport, with a loud one-time warning at load:
      "System loaded from eps-truncated slice files (orders 0..kmax).
      Exact sector spectra cannot be certified; SolveAtPoint at singular
      charts, EndpointLimit, IntegrateOverLine, and transport THROUGH or
      TO singular points are disabled for this system.  Re-export with
      ExportGeneralMatrix for full functionality."
      Rationale: RewritePlan I1 prerequisite; legacy review findings 1-2;
      math review finding 13.
  (d) CANONICAL d_1.m FORMAT (MatrixLoading.m:60-85): REJECTED with a loud
      error in v1; the format is DECLARED LEGACY-PARITY-ONLY (DEC-23).  Its
      only in-repo user is the five-point non-planar classic example whose
      alphabet carries square roots
      (Reference/Examples/FivePointNonPlanar_example.m:42-66, 117), which is
      out of scope v1 anyway (NON-GOALS).  Error text must name the file
      and say "canonical dlog format is not supported by DiffExp2 v1; run
      this example against Legacy/ (frozen old library)".

Validation performed at load (all from old LoadMatrices, kept):
  - at least one matrix file found, else loud error (MatrixLoading.m:31-33);
  - variable names parsed from filenames (MatrixLoading.m:35-50); collision
    with LineParameter is a loud error (MatrixLoading.m:52-54);
  - all matrices square and of equal dimension, else loud error naming every
    file and its dimensions (MatrixLoading.m:127-134, 132-134);
  - entries contain only heads {Association, List, Complex, Integer, Plus,
    Power, Rational, Symbol, Times}, else loud error naming the offending
    heads and file (MatrixLoading.m:148-157);
  - variables appearing in entries but lacking a d<var> file: loud error
    naming them (MatrixLoading.m:97-100, 118-120);
  - ANY half-integer or higher fractional power of an x-dependent quantity
    (Power[a_, b_] with Denominator[b] >= 2): loud error naming the factor
    and file — "matrices with irrational x-dependence are out of scope for
    DiffExp2 v1" (RewritePlan NON-GOALS; replaces the old sqrt machinery
    MatrixLoading.m:159-212 and State.m:227-230, which are ledger items);
  - irreducible singular factors computed and stored exactly
    (MatrixLoading.m:217-231), eps NOT set to zero when computing them in
    Exact mode (the eps-dependent denominators of FT systems are real
    apparent singularities; the old `/. \[Epsilon] -> 0` at
    MatrixLoading.m:228 loses them — DiffExp2 keeps both the eps->0 factor
    set for segmentation and the full set for Indicial).

Options: "Variables" -> {syms} (override filename parsing, old
FEC[Variables] behavior MatrixLoading.m:27-28).  Everything else (Verbosity,
WorkingPrecision, ...) comes from Config.m; LoadSystem takes NO numerical
options of its own.

Input form (DEC-23, binding): FILES ONLY in v1 — LoadSystem accepts a
directory path, never an in-memory matrix Association; unit tests write
temp files.  Together with (d) above this fixes the v1 ingestion surface:
exact full format and closed form are first-class, legacy slices are
parity-fenced, canonical dlog d_1.m is legacy-parity-only.

Old entry point mirrored: LoadConfiguration[... MatrixDirectory -> dir ...]
triggering LoadMatrices (DiffExp/DiffExp.m:171-179, MatrixLoading.m:21-232).
Intentional differences:
  - loading is an explicit verb, not a config side effect.  LoadConfiguration
    with MatrixDirectory still works (the API.m re-export wrapper forwards;
    see 2.8) so the classic call pattern survives unchanged;
  - the exact full format is new (consumed, not just written);
  - sqrt matrices and canonical d_1.m loudly rejected (were supported);
  - returns SystemInfo instead of nothing.

### 2.2 PrepareBoundaryConditions

    PrepareBoundaryConditions[bcs_List, anchor_Association | anchor_List]
        ->  BoundaryConditions (Association, sec. 3.5)

Mirrors old PrepareBoundaryConditions (old Transport.m:31-121).  `anchor` is
a point (<|t -> -1|>) or a line (<|t -> -1/x|>, {t -> endpoint x}); List
input is converted to Association and KeySorted (old Transport.m:36).  Each
entry of `bcs` (one per master integral, length must equal NumIntegrals —
loud error on mismatch, hardened from the old warning at Transport.m:43-45)
is one of:

  (i)  the string "?": per-integral wildcard — this integral's boundary
       data is ignored during constant fixing (old Transport.m:63-65,
       consumed at old Transport.m:288-294 where "?" positions are dropped
       from the matching rows).  Used by Banana_example.m:32-33 and
       test_topiecewise.m:37-38.
  (ii) a List of eps coefficients {c0, c1, ...}, first entry = coefficient
       of eps^0 (old Transport.m:56-62).  Too few coefficients for the
       requested window: loud error naming the integral index and the
       needed order (old Transport.m:59-61).
  (iii) a LaurentValue <|"MinPower" -> m, "Coefficients" -> {...}|> (NEW):
       exact Laurent boundary data.  This replaces the FT side-channel of
       separate EpsPrefactors arrays (Scripts/run_ft_stepwise.m:142-148).
  (iv) a closed-form expression in eps (and the anchor variables), e.g. the
       equal-mass banana Gamma-function expressions
       (Reference/Examples/Banana_example.m:34-40).  Auto-expanded in eps to
       the needed window (old Transport.m:67-74 via SeriesAlways).  Expansion
       starting below the window minimum is a loud error in the (ii)-list
       sense; the old "boundary conditions should start at finite order"
       check (Transport.m:72) generalizes to "start at or above the
       requested window minimum".
  (v)  an x-SeriesData asymptotic value, e.g. 0 + O[x]^(1/2)
       (Reference/Examples/MultiplePolylogarithms_example.m:99-102).
  (vi) expressions containing free symbolic indeterminates (any symbol other
       than x, Logx, eps, "?"): carried linearly through transport (old
       Transport.m:565); see TransportTo for the Pade consequence.

Semantics on a LINE anchor (asymptotic boundary conditions): after eps
expansion, each coefficient is evaluated on the line and its leading
asymptotics extracted (old Transport.m:91-106, LeadingCoefficientSeries with
2 orders); powers x^(a+b*eps) arising from the substitution (e.g.
(-1/t)^(1+3eps) /. t -> -1/x) are parsed into EXACT sector tags {a, b, p}
of the t=anchor chart wherever possible — the parser is SectorSeries.m's
ParseTaggedPower (SectorSeries.md 2.11; sole caller is this function;
REVIEW-minimalism 13), whose inert $FailedParse return selects the
documented minimum-contract branch below (data for a documented branch,
never a silent fallback) — this is the sector-native upgrade
of the old eps-expansion into Logx polynomials; the minimum contract is the
old behavior: expand, rewrite Log[a*x] -> Log[a] + Logx and Log[x] -> Logx
(old Transport.m:120), and zero the sub-leading-in-eps orders below the
expression's leading eps order with an Info print per zeroed order (old
Transport.m:76-84 — this print is MANDATORY; the zeroing is a documented
semantic, not a fallback).  On a POINT anchor, values must be x-free, else
loud error (old Transport.m:102-104).

Precision: any input below WorkingPrecision is raised to 2x WorkingPrecision
with a warning (old Transport.m:38-41; ledger lesson "2x WP raise").

Intentional differences: returns a named BoundaryConditions association
instead of the old positional {line, bcs1} pair; TransportTo accepts BOTH
(section 2.3) so old scripts run unchanged.

### 2.3 TransportTo

    TransportTo[bcs, target_Association | target_List]                 (to = 1)
    TransportTo[bcs, target, to_?NumericQ]
    TransportTo[bcs, target, to_?NumericQ, saveExpansions_?BooleanQ]
    TransportTo[bcs, target, opts___Rule]    (* "SaveExpansions" -> True etc. *)
        ->  TransportResult (Association, sec. 3.6)

Mirrors old TransportTo (old Transport.m:514-1247) including its exact
positional signature TransportTo[bcs, line, to:1, save:False] — the
positional boolean is used by Banana_example.m:52,144, test_banana_refactored.m:90,
test_topiecewise.m:60, and FT TransportLevel
(FeynmanTrick/DiffExpIntegration.m:440-446, 474-480) and MUST keep working.
The old 5th positional argument SampleAtList (old Transport.m:514,
870-897) is DROPPED — no in-repo user outside Transport.m itself; migration
note: use "SaveExpansions" -> True and evaluate the returned SegmentData /
ToPiecewise closures.

Accepted `bcs` forms (all old forms, old Transport.m:530-537 and 146-192):
  (1) BoundaryConditions association from PrepareBoundaryConditions;
  (2) the old positional pair {anchorAssociation, valuesMatrix} (FT
      TransportLevel passes {startPoint, bcs},
      FeynmanTrick/DiffExpIntegration.m:440-446);
  (3) a previous TransportResult association — CHAINING (old
      Transport.m:531-533 reads "KinematicPoint" and "SeriesValues").
      NEW: chaining also ingests "ErrorEstimates" so accumulated error is
      carried across calls; the old assoc-chaining path silently dropped it
      and re-derived errors from Accuracy (old Transport.m:582-587) — that
      Quiet block is deleted; if the chained result has no ErrorEstimates
      (EstimateError was False) the new accumulation starts at 0 and the
      result says so via "ErrorEstimates" -> Missing["NotComputed"];
  (4) the old 3-element form {anchor, values, errors} (old
      Transport.m:579-581).
  Chaining from a result with "MultivaluedFail" -> True is a loud error
  (old code silently allowed re-use of the partial result).

Target semantics (all old, must be preserved):
  - target a LINE (depends on x): transport along it from FixAt to `to`,
    where FixAt is derived from the bcs anchor by relating lines
    (LineSegmentation.m:148-..., CheckBoundaryConditionsAndReparametrize);
  - target a POINT, bcs anchored at a point: the straight line
    p1*(1-x) + x*p2 is constructed automatically and to = 1 (old
    Transport.m:547-550) — pattern TransportTo[Results1, <|t -> 1/2|>]
    (Banana_example.m:90, test_unequal_mass_full.m:58);
  - target a POINT, bcs anchored on a line: the target must lie on that
    line, found via RelateLines, else loud error "Endpoint does not lie on
    same line as the boundary conditions" (old Transport.m:551-554);
  - BIDIRECTIONAL: `to` may be on either side of FixAt; to < FixAt
    transports backwards (old Transport.m:691-695, 1112-1116; pattern
    TransportTo[Results2, <|t -> x|>, -1], test_banana_refactored.m:138).
  - multi-variable lines: <|psq -> 1/2, mm1 -> 1+x, mm2 -> 1+x/2, ...|>
    (Banana_example.m:128-132); the line must fix ALL system variables,
    else loud error naming the missing ones (old Transport.m:522-524);
  - the line must not lie ON a singularity: loud error naming the vanishing
    irreducible factors (old Transport.m:571-574);
  - non-linear lines (x appearing at powers other than +-1) rejected loudly
    (LineSegmentation.m:153-155).

Option (DEC-10; closes FTShimContract G1):
  `"ExtraSingularFactors" -> {poly..}` — exact polynomials in the pinned
  variables, unioned EXACTLY into the segmentation alphabet for this call
  (charts are placed at every real root in range; complex roots cap radii).
  Loud error if any "Combination"/"Prefactor" denominator has a root not in
  the alphabet (closure check; the implicit-widening alternative is rejected
  per FTShimContract G1).  Threaded to Transport.m segmentation; same option
  on IntegrateOverLine (2.6).

SINGULAR-ENDPOINT MODE (old Transport.m:605-613 detection, 1053-1096 final
segment; test_singular_endpoint.m): if `to` is a singular point of the
system on this line, the final chart is solved but NOT evaluated; the result
carries "EndpointIsSingularity" -> True and "LocalSolution" -> the
LocalSolution object at the target chart (RewritePlan 3.1) INSTEAD of
numbers in "SeriesValues" — intentional upgrade: the old code stuffed
analytic-continued SeriesData tables into "SeriesValues" (old
Transport.m:1084-1085); DiffExp2 returns the typed object that EndpointLimit
and SolveAtPoint consume.  For compatibility, "SeriesValues" at a singular
endpoint holds the per-(integral, eps-order) x-series evaluation forms of
that LocalSolution (Normal-izable, Logx-bearing), matching what
test_singular_endpoint.m:119-140 asserts.  Endpoint singularity detection is
EXACT: `to` is compared against the exact singularity roots; an INEXACT
numeric `to` within `InputSnapTol[to]` of an exact root is snapped WITH A
WARNING (W5) — the input's own accuracy, never wp, sets the snap scale
(Tolerances.md InputSnapTol + $NearSingularityGuardDecades entries; DEC-14;
REVIEW-math D19, REVIEW-minimalism 11), and an inexact `to` in the
suspicious band between InputSnapTol[to] and
10^(-$NearSingularityGuardDecades) is the LOUD Transport.md E14, never a
guess; exact `to` uses exact membership only.  The old numeric-chop
membership test (old Transport.m:605-609) is replaced (a chop miss there
silently attempted evaluation AT the singularity).

Eps-window honesty: "SeriesValues" is a NumIntegrals x (CompleteMax-Min+1)
matrix; column k is the coefficient of eps^(Window["Min"]+k-1).  For classic
systems with Min = 0 this is bit-identical in shape to the old
NumIntegrals x (EpsilonOrder+1) matrix, so existing consumers
(test_topiecewise.m:97, Banana_example.m:98-112, FivePointNonPlanar_example.m:198)
migrate without index changes.  "EpsilonOrder" is kept as an alias for
Window["CompleteMax"] (old key, old Transport.m:1244).

Symbolic indeterminates: if the bcs carry free symbols, transport is linear
in them; UsePade is force-disabled with a LOUD warning and recorded as
"PadeDisabled" -> True in the result (old Transport.m:565-569 prints Info
and mutates global config — the global mutation is kept for compat but the
result metadata is new); error estimates are computed per indeterminate (old
Transport.m:947-962, ComputeErrorsPerIndeterminate) and the reported
per-(integral, order) error is the max over indeterminate coefficients.

Underdetermined boundary data (including via "?" wildcards): named free
parameters Subscript[c, epsord, intind, i] are introduced with a loud
warning, exactly the old behavior (old Transport.m:393-410, 418-434),
never silently zeroed; Pade disabled as above.

Error estimation: "EstimateError" -> "Fast" (default) runs the two-point
full-vs-reduced-order probe per segment at both the matching point and the
evaluation point, avoiding x = 0 where logs vanish spuriously (old
Transport.m:905-925); errors accumulate additively across segments per
(integral, eps order) (old Transport.m:980-984); a segment error > 1 aborts
loudly (old Transport.m:991-993).  Output key "ErrorEstimates" has the same
shape as "SeriesValues" (consumed by Tests/test_multiple_polylogarithms.m:106-107
and Tests/test_five_point_nonplanar.m:249).  "EstimateError" -> False yields
"ErrorEstimates" -> Missing["NotComputed"] (old Transport.m:1238-1239),
NEVER a zero array.

Precision/budgeting lessons (implemented in Transport.m, surfaced here
because the API owns the user knobs): predivision two-pass segment counting
and DigitsNeeded = AccuracyGoal + Ceiling[Log10[#segments]] + safety (old
Transport.m:666-760, 758-760); AccuracyGoal numeric => always-on
per-segment accuracy validation (Transport.md E11; DEC-11 — the old
AccuracyGoalValidate "Before" adaptive expansion-order search, old
Transport.m:776-841, and "After" redo-segment, old Transport.m:1191-1213,
are NOT ported; ledger waiver "replaced by exact recursion + ErrorEstimate
gate"; REVIEW-minimalism 8); input precision raise to 2x WP (old
Transport.m:541-544, 559-562).

Old entry point mirrored: TransportTo (old Transport.m:514).  IntegrateSystem
as a public name is DROPPED (it was both engine and API, old
Transport.m:146-511); its two public uses migrate:
  - IntegrateSystem[bcs, singularLine] (Tests/test_decomposition.m:58)
    -> SolveAtPoint (section 2.4);
  - IntegrateSystem[line] general-solution mode (old Transport.m:24, 146
    default bcs = "?") -> SolveAtPoint["?", ...].

### 2.4 SolveAtPoint

    SolveAtPoint[bcs, point_Association, opts___Rule]  ->  LocalSolution
    SolveAtPoint["?", point_Association, opts___Rule]  ->  LocalSolution
        (general solution with symbolic constants)

NEW entry point (FT-pipeline first-class).  Constructs the LocalSolution
(RewritePlan 3.1) at `point` — which MAY be a singular point of the system —
by chart construction (Indicial.m -> Solve.m).  With bcs: integration
constants are fixed by matching the bcs at their anchor inside the chart
(the old FixAt = 0 asymptotic match, old Transport.m:296-307, and the
numeric-FixAt match, old Transport.m:317-326, become one object-level
matching in Transport.m/Solve.m).  With "?": all constants stay symbolic
Subscript[c, i] (old IntegrateSystem default, old Transport.m:146,
439-441).

Options: "Direction" -> +1|-1 (which side of the singular point the chart
is oriented toward; determines the iδ branch via the Prescriptions record).
Loud errors: LegacySlices mode (section 2.1c); indicial contract violation,
rank-reduction non-termination, missing/conflicting prescriptions — all
propagated from Indicial.m/Solve.m with the chart named (their specs own the
text; API adds the user-facing point coordinates).

Old functionality replaced: the FT pattern "solve on a line through/at a
singular point and inspect the local series" currently done by
IntegrateSystem + DecomposeSingularity
(Tests/test_decomposition.m:58; FeynmanTrick/DiffExpIntegration.m:822) —
DecomposeSingularity becomes unnecessary because the LocalSolution carries
exact sector tags (RewritePlan I2; SectorSeries.m absorbs
SingularityDecomposition.m's role per RewritePlan 3.2).

### 2.5 EndpointLimit

    EndpointLimit[src, opts___Rule]  ->  {LaurentValue ..}
        src: TransportResult with "EndpointIsSingularity" -> True,
             or LocalSolution

Computes lim_{t -> 0} (chart coordinate) of each component, per RewritePlan
3.3: the limit is the eps-Laurent constant of the (a,b,p) = (0,0,0) sector;
sectors with b != 0 are dropped EXACTLY (dimreg convention, even when a < 0
— stated, it is the drop rule); a b = 0 sector with a < 0, or with a = 0 and
p > 0, that does not cancel in the assembled combination is a LOUD
divergence error naming (component index, sector tag {a,b,p}, eps order).

Options:
  "Direction" -> +1|-1 : side of approach (resolves the branch via the
      chart's Prescriptions; mirrors the old direction argument of
      EvaluateEndpointLimitSectors, DiffExp/RegularizedIntegration.m:59,
      1965, called at FeynmanTrick/DiffExpIntegration.m:1014-1018).
  "Combination" -> coefficient matrix {{c_ij(x, eps)}} or vector, rational
      in (x, eps): the combination Σ_j c_ij f_j is formed AT THE OBJECT
      LEVEL (SectorSeries rational-multiply + add) BEFORE the limit, so
      exact pole cancellations happen symbolically.  This is the
      combine-before-limit hardening of A1: it replaces FT's
      EvaluateLimitFromTransport with its Quiet[Check[..., 0]] chains
      (FeynmanTrick/DiffExpIntegration.m:908-1110, silent zeros at
      1037-1059 and 1080-1093) — those fallbacks are FORBIDDEN here; any
      coefficient that cannot be evaluated/limited is a loud error naming
      the coefficient index and the chart.

Old functionality replaced: EvaluateEndpointLimitSectors
(DiffExp/RegularizedIntegration.m:1965-2140) including its salvage/drop
endings — the "endpoint tower with unresolved sectors" and "dropped
unresolved sector content" warnings (RegularizedIntegration.m:2125-2132)
become loud errors (they cannot occur when tags are exact, but the assert
stays); EvaluateLimitAtSingularity (RegularizedIntegration.m:57).

### 2.6 IntegrateOverLine

    IntegrateOverLine[result_Association, {x0_, x1_}, opts___Rule]
        ->  {LaurentValue ..}    (one per output component)

Integrates the transported solutions (result must carry SegmentData, i.e.
TransportTo was called with SaveExpansions; loud error otherwise, mirroring
old DefiniteIntegralWithPrefactorLaurent's "No segment data available",
DiffExp/RegularizedIntegration.m:2562) over [x0, x1] in main-line
coordinates, across interior singular points, with exact Laurent
bookkeeping.  Delegates to Integrate.m (per-sector closed forms, full
{b=0/b!=0} x {a+n+1 <0/=0/>0} x {p=0/p>0} case table, interior-pole PV/iδ
half-segment pairing) and SectorSeries.m (rational multiply, partial
fractions across charts).

Options:
  "Prefactor" -> <|"PowerAtLower" -> v1(eps), "PowerAtUpper" -> v2(eps),
      "RationalFactor" -> r(x, eps)|> : integrand prefactor
      (x - x0)^v1 (x1 - x)^v2 r(x, eps), matching the old prefactorSpec
      consumed by DefiniteIntegralWithPrefactor[Laurent]
      (DiffExp/RegularizedIntegration.m:71-76;
      FeynmanTrick/DiffExpIntegration.m:868-877).  Default: no prefactor.
  "Combination" -> as in 2.5, applied at object level before integration
      (exact cancellation; RewritePlan 3.3 "integrate" case).

Hard requirements:
  - charts in SegmentData must be affine (Mobius charts REJECTED loudly,
    RewritePlan 3.2 Integrate.m; old FT contract "Required for
    integration!", FeynmanTrick/DiffExpIntegration.m:352);
  - LegacySlices mode rejected (section 2.1c);
  - the SegmentData chain must cover [x0, x1] without gaps (assert; the old
    code silently integrated whatever segments existed);
  - window arithmetic per EpsSeries.m: the returned LaurentValue windows
    reflect every endpoint-pole enhancement (depth p+1 per RewritePlan 3.4)
    and carry honest CompleteMax.

Old functionality replaced: DefiniteIntegralWithPrefactorLaurent
(DiffExp/RegularizedIntegration.m:2549-...) and the supporting
DefiniteIntegral / IntegratePiecewiseSaved / IntegrateSegmentData /
IntegrateDecomposition[Laurent] / IntegrateSingularTerm[Laurent] /
RegularizeIntegrand stack (RegularizedIntegration.m:45-76 usage block) —
per-export delete-vs-reimplement dispositions live in
Docs/ExportDisposition.md; the API-level contract is: ONE entry point, no
Laurent side-channels (epsMinPower threading is inside LaurentValue).

### 2.7 ToPiecewise

    ToPiecewise[result_Association]                       ->  matrix of closures
    ToPiecewise[result_Association, ord_Integer]          ->  truncated variant

Mirrors old ToPiecewise (old Transport.m:1255-1313; usage
Banana_example.m:57, 152; Tests/test_topiecewise.m:67).  Consumes a
TransportResult with SegmentData and returns a NumIntegrals x WindowLength
matrix of single-argument closures: ToPiecewise[res][[i, j]][tval] evaluates
integral i, eps-column j at main-line parameter tval, choosing the covering
segment and evaluating its LocalSolution through SectorSeries (which owns
Pade acceleration; the old positional pade flag, old Transport.m:1255 and
its GetPade use at 1293-1298, is DROPPED — evaluation acceleration is a
Config concern (UsePade), not a per-call flag.  Migration note: delete the
second argument; behavior with UsePade -> True is equivalent).
Out-of-range tval: loud error naming the covered interval (the old
Piecewise default silently returned 0 outside all segment conditions — a
forbidden fallback).  Compressed segment data (SaveExpansionsCompress
workflow, old Transport.m:1278-1287) is supported by uncompressing
transparently; the file-hash export variant is kept only if budget allows
(section 9).

### 2.8 Configuration functions (implemented in Config.m, surfaced here)

    LoadConfiguration[optsList]      (resets to defaults, then applies)
    UpdateConfiguration[optsList | rules__]
    CurrentConfiguration[]

Contract and key table: Docs/specs/Config.md.  API-relevant guarantees:
LoadConfiguration with MatrixDirectory triggers LoadSystem (compat with
every classic example, e.g. Banana_example.m:24, 84, 119;
FeynmanTrick/DiffExpIntegration.m:367); DeltaPrescriptions parsing of
poly + I*\[Delta] forms (old DiffExp/DiffExp.m:140-148) and irreducibility
validation (old DiffExp/DiffExp.m:164-169) are Config.m's; unknown keys are
a loud error (the validated accessor kills the FEC-class silent-miss bug).

### 2.9 Compatibility matrix (every observed old call pattern)

Verified by grep over Reference/Examples/*.m, Tests/*.m, Scripts/*.m,
FeynmanTrick/*.m.  "KEEP" = must work unchanged (modulo context rename
DiffExp` -> DiffExp2`).  "BREAK" = intentionally broken, migration note
given.

| # | Pattern (with a real citation) | Status |
|---|---|---|
| 1 | PrepareBoundaryConditions[{...closed forms, "?"...}, <\|t -> -1/x\|>] (Banana_example.m:31-41) | KEEP |
| 2 | PrepareBoundaryConditions[series-valued bcs, {t -> endpoint x}] with 0 + O[x]^(1/2) entries (MultiplePolylogarithms_example.m:94-102) | KEEP |
| 3 | PrepareBoundaryConditions[numericValues, pointList] (FivePointNonPlanar_example.m:181) | KEEP |
| 4 | TransportTo[preparedBcs, <\|t -> -1\|>] point target (Banana_example.m:47) | KEEP |
| 5 | TransportTo[resultAssoc, <\|t -> x\|>, 10, True] chaining + to + save (Banana_example.m:52) | KEEP |
| 6 | TransportTo[resultAssoc, <\|t -> 1/2\|>] auto-line between points (Banana_example.m:90) | KEEP |
| 7 | TransportTo[bcs, multiVarLine, 1] (Banana_example.m:128-132) | KEEP |
| 8 | TransportTo[result, <\|t -> x\|>, -1] backward transport (test_banana_refactored.m:138) | KEEP |
| 9 | TransportTo[{startPoint, bcsMatrix}, <\|xx -> 0\|>, 1, True] singular endpoint, raw pair (DiffExpIntegration.m:440-446) | KEEP |
| 10 | TransportTo[bcs, line] default to = 1 (test_recurrence_speedup.m:42) | KEEP |
| 11 | result["SeriesValues"][[i, j]], ["ErrorEstimates"][[i, 1]], ["KinematicPoint"], ["NumIntegrals"], ["EndpointIsSingularity"], ["EpsilonOrder"], ["ComputationTime"] (test_multiple_polylogarithms.m:106-107, test_five_point_nonplanar.m:249, test_singular_endpoint.m:105-140, test_topiecewise.m:77) | KEEP (keys preserved; "EpsilonOrder" aliases Window CompleteMax) |
| 12 | ToPiecewise[result]; ToPiecewise[result][[i, j]][5] (Banana_example.m:57-59, 152-156) | KEEP |
| 13 | LoadConfiguration / UpdateConfiguration / CurrentConfiguration (Banana_example.m:24-27, 84-86, 120) | KEEP (Config.m) |
| 14 | Both eps and \[Epsilon] in matrices and bcs (MatrixLoading.m:118; Banana_example.m:34-40 uses \[Epsilon], FT uses eps) | KEEP |
| 15 | IntegrateSystem[result, singularLine] (test_decomposition.m:58) | BREAK -> SolveAtPoint[result, point]; the singular line's center becomes the point argument |
| 16 | IntegrateSystem[line] general solutions (old Transport.m:24) | BREAK -> SolveAtPoint["?", point] |
| 17 | TransportTo 5th argument SampleAtList (old Transport.m:514) | BREAK -> SaveExpansions + ToPiecewise (no in-repo user) |
| 18 | Positional SegmentData 5-tuples seg[[1..5]] (test_singular_endpoint.m:144-176; DiffExpIntegration.m:916-995) | BREAK -> named SegmentRecord (sec. 3.7); FT call sites are rewritten at M5 to EndpointLimit/IntegrateOverLine; test_singular_endpoint is retargeted in M6 |
| 19 | ToPiecewise[savedData, True] Pade flag (old Transport.m:1255) | BREAK -> Config UsePade (no in-repo user passes it) |
| 20 | DiffExp`Transport`TransportTo / DiffExp`RegularizedIntegration`* / DiffExp`SingularityDecomposition`DecomposeSingularity fully qualified reach-ins (test_transport_level.m:103; DiffExpIntegration.m:822, 875, 1015) | BREAK -> DiffExp2` public symbols per the FT shim contract (M0 task 3) |
| 21 | Five-point non-planar classic example (canonical d_1.m + sqrt alphabet, FivePointNonPlanar_example.m:42-66, 117) | BREAK -> runs against Legacy/ only (v1 NON-GOAL); M6 battery disposition |
| 22 | FT EpsPrefactors side-channel for Laurent boundary data (run_ft_stepwise.m:142-153) | BREAK -> LaurentValue bcs entries (sec. 2.2 iii) at M5 cutover |

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution (verbatim from RewritePlan 3.1)

    LocalSolution = <|
      "Center" -> exact x0, "ChartMap" -> affine (Mobius optional, see 3.2),
      "Radius" -> distance to nearest singularity IN THE COMPLEX PLANE
                  (complex singularities are real: pentagon/unequal-mass
                  lines have them; old code projects ghosts Re, Re±Im —
                  ledger item; new code uses true complex distance),
      "Sectors" -> { Sector.. },
      "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
      "TWindow"   -> <|"CompleteMax" -> nmax|>,   (* t-order honesty: top
                  (couplingDepth−1) orders of chained particular solutions
                  are degraded (old MaxCouplingOrder/ISafetyExpansionSubtract
                  lesson); tracked, not guessed *)
      "ErrorEstimate" -> per (eps-order) accumulated error (two-point
                  full-vs-reduced-order probe, additive across segments,
                  abort > 1 — ported old machinery, user-facing),
      "Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
                  "LeadingCoeffSign"|> with the derived chart Im-sign;
                  consistency-checked at construction (even multiplicity = no
                  constraint; conflict or missing prescription at a chart with
                  b!=0/p>0 sectors = LOUD ERROR; sqrt factors auto-prescribed
                  as in old State.m DEqnSquareRoots)
    |>

### 3.2 Sector (verbatim from RewritePlan 3.1)

    Sector = <|
      "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
      "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT:
                  solutions are vectors; matching and ODE residuals need all
                  components *)
    |>
    Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                    Σ_k eps^k Σ_n c[k,n,comp] t^n

### 3.3 EpsWindow (verbatim from RewritePlan 3.1)

    EpsWindow = <|"Min" -> kmin, "CompleteMax" -> kmax|>

Consumers check need <= CompleteMax and fail loudly naming (chart, sector,
order).  Exemption: homogeneous true-resonance sectors have kmin = -p by
construction; window arithmetic must account, not error (RewritePlan 3.1
invariants).

### 3.4 SystemInfo (returned by LoadSystem)

    <|
      "Variables" -> {sym..},          (* parsed or overridden *)
      "NumIntegrals" -> n_Integer,
      "Format" -> "Exact" | "LegacySlices",
      "SliceMax" -> k_Integer | None,  (* LegacySlices only: highest order *)
      "MatrixDirectory" -> dir_String,
      "IrreducibleFactors" -> {poly..},        (* eps -> 0 projection *)
      "IrreducibleFactorsExact" -> {poly..}    (* with eps, Exact mode *)
    |>

### 3.5 BoundaryConditions

    <|
      "Anchor" -> Association (line or point, KeySorted, precision >= 2x WP),
      "Values" -> {entry_1, ..., entry_NumIntegrals},
          entry = "?" | LaurentValue | <|"Sectors" -> {TaggedBC..}|>,
          TaggedBC = <|"a", "b", "p", "Coeffs" -> c[[k, n]]|>  (* parsed
              asymptotic data on a line anchor; on a point anchor every
              entry is a LaurentValue with numeric/symbolic coefficients *)
      "Indeterminates" -> {sym..},     (* free symbols, possibly {} *)
      "ErrorEstimates" -> matrix | 0,  (* carried accumulated errors *)
      "EpsWindow" -> EpsWindow
    |>

### 3.6 TransportResult

    <|
      "KinematicPoint" -> Association,           (* old key, old Transport.m:1235 *)
      "EndpointIsSingularity" -> True | False,   (* old key, :1236 *)
      "SeriesValues" -> matrix [i, k],           (* old key, :1237; column k =
          coefficient of eps^(Window Min + k - 1); numbers at a regular
          endpoint, x-series forms at a singular endpoint *)
      "LocalSolution" -> LocalSolution | Missing["RegularEndpoint"],  (* NEW *)
      "ErrorEstimates" -> matrix (same shape) | Missing["NotComputed"],  (* old key, :1238 *)
      "SegmentData" -> {SegmentRecord..} | Missing["NotRequested"],    (* old key, :1240 *)
      "EpsWindow" -> EpsWindow,                  (* NEW; honest *)
      "EpsilonOrder" -> Window CompleteMax,      (* old alias key, :1244 *)
      "NumIntegrals" -> n,                       (* old key, :1243 *)
      "PadeDisabled" -> True | False,            (* NEW *)
      "MultivaluedFail" -> True | False,         (* NEW; explicit *)
      "ComputationTime" -> seconds               (* old key, :1242 *)
    |>

### 3.7 SegmentRecord (replaces the old positional 5-tuple, old Transport.m:847-867)

    <|
      "Chart" -> LocalSolution,        (* with constants fixed *)
      "LineMap" -> x -> expr,          (* local t -> main-line x relation;
                                          old element 2 *)
      "MainInterval" -> {xlo, xhi},    (* old element 3 *)
      "LocalInterval" -> {tlo, thi}    (* old element 4 *)
    |>

The old element 5 (compressed per-(integral, eps-order) SeriesData) is
subsumed by "Chart"; compression (SaveExpansionsCompress) applies to the
whole record.

### 3.8 LaurentValue (FT-compatible)

    <|"MinPower" -> m_Integer, "Coefficients" -> {c_m, c_{m+1}, ...}|>

Identical to the FT layer's existing Laurent convention
(Scripts/run_ft_stepwise.m:47-53), with the DiffExp2 addition that any
LaurentValue produced by EndpointLimit/IntegrateOverLine also carries
"CompleteMax" -> kmax (honest top).  Consumers indexing beyond CompleteMax
is a loud error, never zero-fill.

---

## 4. INVARIANTS (always on, cheap)

1. EXACTLY-ONE-SYSTEM: every entry point except LoadSystem asserts a system
   is loaded ("No system loaded — call LoadSystem[dir]"); LoadSystem
   replaces ALL system state atomically (no partial mutation on error —
   validate fully before installing).
2. SHAPE: Dimensions[SeriesValues] === {NumIntegrals, CompleteMax - Min + 1};
   ErrorEstimates has the same dimensions or is Missing; asserted on every
   result construction and on every chained input.
3. WINDOW SANITY: Min <= CompleteMax on every EpsWindow/LaurentValue; window
   never empty; chaining never widens a window (the chained request must
   satisfy need <= incoming CompleteMax or error loudly naming the integral
   and order).
4. PRECISION FLOOR: after ingestion, every numeric input has precision
   >= 2x WorkingPrecision (raised with warning, old Transport.m:38-41,
   541-562); asserted before delegation.
5. ANCHOR CONSISTENCY: point anchors have x-free values (old
   Transport.m:102-104); line anchors and targets fix all system variables
   (old Transport.m:166-168, 522-524); BC entry count equals NumIntegrals.
6. REGULAR-CHART DEGENERACY (RewritePlan 3.1): a TransportResult at a
   regular endpoint has values free of Logx and fractional powers; asserted
   cheaply by head inspection.
7. SEGMENT TILING: SegmentData MainIntervals tile the transported range
   contiguously (no gaps/overlaps beyond shared endpoints); asserted at
   result construction and again by ToPiecewise/IntegrateOverLine.
8. LEGACY-MODE FENCE: SystemInfo["Format"] === "LegacySlices" implies the
   singular-operation gate (sections 2.1c, 5.E1) — checked at the TOP of
   SolveAtPoint, EndpointLimit, IntegrateOverLine, and inside TransportTo
   when the segmentation reports a singular chart on the path or a singular
   target.
9. ERROR MONOTONICITY: chained accumulated errors are non-negative and
   non-decreasing across a transport call.

---

## 5. ERROR CONTRACT

NO silent fallbacks anywhere.  Every error below is "loud" = raised through
the single error-reporting function with the named payload; warnings are
printed AND recorded in result metadata where indicated.  Message text must
carry, at minimum, the items listed in [brackets].

LoadSystem:
  E1  Singular operation on LegacySlices system [operation name; directory;
      slice max; remedy: ExportGeneralMatrix].  Fired by SolveAtPoint /
      EndpointLimit / IntegrateOverLine always, and by TransportTo when the
      path contains or targets a singular chart.
  E2  No matrix files in directory [directory] (MatrixLoading.m:31-33).
  E3  Mixed/ambiguous formats (both _full.m and slices present for the same
      variable) [variable; both filenames].  Old code had no such check.
  E4  Dimension mismatch [every filename with its dimensions]
      (MatrixLoading.m:127-134).
  E5  Unsupported heads [heads; filename] (MatrixLoading.m:148-157).
  E6  Fractional powers of x-dependent factors [factor; filename] — v1
      out-of-scope sqrt error (NON-GOALS; replaces MatrixLoading.m:159-212).
  E7  Canonical d_1.m present [filename; Legacy remedy] (sec. 2.1d).
  E8  Variable/LineParameter collision [symbol] (MatrixLoading.m:52-54).
  E9  Entry variables without derivative files [variables; filename]
      (MatrixLoading.m:97-100, 118-120).
  W1  WARNING + SystemInfo flag: LegacySlices mode loaded (sec. 2.1c).

PrepareBoundaryConditions / bcs ingestion:
  E10 BC count != NumIntegrals [given count; expected; system dir] —
      hardened from old warning (old Transport.m:43-45).
  E11 Too few eps coefficients [integral index; supplied window; needed
      window] (old Transport.m:59-61, generalized to Laurent windows).
  E12 Closed form starts below the requested window minimum [integral
      index; leading eps power; window min] (old Transport.m:72).
  E13 Point-anchored values depend on x [integral index]
      (old Transport.m:102-104).
  W2  WARNING: precision raised to 2x WP (old Transport.m:38-41).
  I1  INFO (mandatory print): sub-leading eps orders of a closed form
      assumed exactly zero on a line anchor [integral index; eps order]
      (old Transport.m:76-84).

TransportTo:
  E14 Target/line does not fix all variables [missing variables]
      (old Transport.m:522-524).
  E15 Line lies on a singularity [vanishing irreducible factors]
      (old Transport.m:571-574).
  E16 Endpoint not on the bcs line [target point; anchor line]
      (old Transport.m:551-554).
  E17 Non-linear line segment [offending powers]
      (LineSegmentation.m:153-155).
  E18 Boundary matching failure [integral indices; eps order; chart center;
      solve residual] (old Transport.m:376-388) — the old least-squares
      rescue (old Transport.m:348-374) is FORBIDDEN: a failed exact solve is
      an error, full context dumped to the diagnostic state, never a
      best-fit substitution.
  E19 Unprescribed branch point [chart center on the main line; the sector
      tags with b != 0 or p > 0 that need continuation; the
      DeltaPrescriptions remedy] — the old GiveMultivaluedError text (old
      Transport.m:470-478) and SingularityCheck (old Transport.m:480-492)
      unified and ON BY DEFAULT.  With config
      "AbortOnAnalyticContinuationFail" -> False (FT sets this,
      FeynmanTrick/DiffExpIntegration.m:411-414): degrade to WARNING +
      "MultivaluedFail" -> True in the result; chaining FROM such a result
      is error E25.
  E20 Prescription sign conflict / indeterminate at a chart [chart center;
      per-factor (factor, sign, multiplicity, leading-coeff sign) table]
      (RewritePlan 3.1 Prescriptions; old AnalyticContinuation.m:18-68
      semantics) — at segment 1 in the negative direction this replaces the
      old "integrate in the positive line direction" abort (old
      Transport.m:1026-1032).
  E21 Segment error estimate > 1 [segment index; chart center; worst
      (integral, eps order); estimate] (old Transport.m:991-993).
  E22 Window shortfall during marching [chart center; sector tag; needed
      eps order; available CompleteMax] (RewritePlan 3.1 invariant).
  E23 Chained error array window mismatch [shapes] — replaces the old
      silent PadRight zero-fill (old Transport.m:589-591), FORBIDDEN.
  E24 Singular target in LegacySlices mode = E1.
  E25 Chaining from a result with "MultivaluedFail" -> True [the failed
      chart] — old code silently reused partial results.
  W3  WARNING + "PadeDisabled" -> True: symbolic indeterminates or
      wildcard-induced free parameters disable Pade (old Transport.m:565-569,
      178-181, 406-410).
  W4  WARNING: free parameters Subscript[c, epsord, intind, i] introduced
      [integral indices; eps order; count] (old Transport.m:406-407,
      427-428).
  W5  WARNING: numeric target snapped to exact singularity [given value;
      exact root; snapTol] (replaces silent chop matching, old
      Transport.m:605-609).

SolveAtPoint: E1, E14, E19, E20 as above, plus errors propagated verbatim
from Indicial.m (char-poly factorization contract; rank-reduction
non-termination) and Solve.m, with the user-facing point coordinates
appended.

EndpointLimit:
  E26 Divergent b = 0 sector not cancelled in the assembled combination
      [component index; sector tag {a,b,p}; eps order] (RewritePlan 3.3 /
      Integrate.m b=0 rule).
  E27 Unresolved sector content [chart; sector tag] — replaces the old
      warn-and-drop (DiffExp/RegularizedIntegration.m:2125-2132), FORBIDDEN
      to drop.
  E28 Combination coefficient not rational in (x, eps) or unevaluable at
      the chart [coefficient index] — replaces FT's Quiet[Check[..., 0]]
      (FeynmanTrick/DiffExpIntegration.m:1037-1059, 1080-1093), FORBIDDEN.

IntegrateOverLine:
  E29 No SegmentData [remedy: SaveExpansions -> True]
      (old RegularizedIntegration.m:2562 precedent).
  E30 Mobius charts present [chart centers] (RewritePlan 3.2 Integrate.m;
      FeynmanTrick/DiffExpIntegration.m:352).
  E31 Coverage gap [requested interval; covered union].
  E32 Endpoint divergence per Integrate.m case table [endpoint; sector tag;
      eps order].
  E33 Interior-pole half-segment pairing violation [pole position; the two
      half-segments] (RewritePlan 3.2 Integrate.m: assert, not assume).

ToPiecewise:
  E34 No SegmentData (old Transport.m:1259-1261).
  E35 Evaluation outside covered range [tval; covered interval] — replaces
      the silent Piecewise-default 0, FORBIDDEN.

Enumerated tempting-fallback sites, all FORBIDDEN (A1):
  F1  treating slice dirs as exact (-> E1/W1);
  F2  least-squares boundary matching (-> E18);
  F3  zeroing unmatched/free constants (-> W4, symbols kept);
  F4  PadRight on chained error arrays (-> E23);
  F5  Quiet[Check[..., 0]] on endpoint evaluation/limits (-> E26-E28);
  F6  warn-and-drop unresolved endpoint content (-> E27);
  F7  Piecewise default 0 outside range (-> E35);
  F8  zero-filled ErrorEstimates when not computed (-> Missing);
  F9  numeric-chop singularity membership (-> W5 exact snap);
  F10 "keep at least one output order" floors and any zero-padded window
      (legacy review finding 13; windows are honest or error);
  F11 silently continuing after MultivaluedFail (-> E25);
  F12 evaluating AT a singular endpoint because detection missed (-> exact
      detection + singular-endpoint mode).

---

## 6. ABSORBED OLD CODE

| Old code (file:line) | Disposition |
|---|---|
| DiffExp/Transport.m:23-26 (exports) | names kept: PrepareBoundaryConditions, TransportTo, ToPiecewise; IntegrateSystem dropped (-> SolveAtPoint) |
| DiffExp/Transport.m:31-121 PrepareBoundaryConditions | API.m sec. 2.2 (parsing); the Logx rewrite (:120), leading-asymptotics extraction (:91-106), eps-expansion (:67-84), wildcard (:63-65), 2xWP raise (:38-41) all preserved; sector-tag parsing is the upgrade |
| DiffExp/Transport.m:146-511 IntegrateSystem | engine split: dispatch/recursion -> Solve.m; boundary-equation assembly (:288-345) -> Transport.m matching; "?" masking (:288-294), free-parameter introduction (:393-441), Pade disabling (:178-181) -> API semantics, sec. 2.3; least-squares rescue (:348-374) DELETED (E18); GeneralSolution crosscheck (:443-452) -> always-on ODE residual invariant (RewritePlan 3.1) |
| DiffExp/Transport.m:470-492 GiveMultivaluedError + SingularityCheck | E19, default ON (old default OFF via CrosscheckFlags "SingularityCheck" -> 0, DiffExp/State.m:205) |
| DiffExp/Transport.m:514-1247 TransportTo | marching/segmentation -> Transport.m; API keeps: signature (:514), assoc chaining (:530-537), error carry (:577-591 minus the Quiet block), auto-line between points (:547-554), bidirectionality (:691-695, 1112-1116), singular-endpoint mode (:605-613, 1053-1096), digit budgeting (:758-760), AccuracyGoalValidate (:776-841, 1191-1213), two-point error probe + per-indeterminate errors + abort>1 (:905-997), result keys (:1234-1246) |
| DiffExp/Transport.m:1255-1313 ToPiecewise | sec. 2.7; Pade flag (:1255, 1293-1298) dropped; compress handling (:1278-1287) kept |
| DiffExp/MatrixLoading.m:21-232 LoadMatrices | LoadSystem sec. 2.1; eps/\[Epsilon] (:118), validation (:127-164), irreducible factors (:217-231) kept; sqrt machinery (:166-212) -> E6 + ledger; canonical d_1 (:60-85) -> E7; per-chart expansion (PrepareMatrices* :252-355) -> Transport.m/SectorSeries.m, NOT API; InitializeIntegrationSequence (:357-385) -> Solve.m (block ordering + MaxCouplingOrder for TWindow) |
| DiffExp/DiffExp.m:86-197 config functions | Config.m (sec. 2.8); ChopPrecision sync (:121-130), DeltaPrescriptions parsing (:140-148, 164-169) per Config spec |
| DiffExp/RegularizedIntegration.m:1965-2140 EvaluateEndpointLimitSectors; :57 EvaluateLimitAtSingularity | EndpointLimit sec. 2.5; warn-and-drop (:2125-2132) -> E27 |
| DiffExp/RegularizedIntegration.m:2549-... DefiniteIntegralWithPrefactorLaurent and the :45-76 export stack | IntegrateOverLine sec. 2.6 (closed forms live in Integrate.m); per-export dispositions in Docs/ExportDisposition.md |
| FeynmanTrick/DiffExpIntegration.m:267-510 TransportLevel | NOT absorbed (stays FT) but its DiffExp surface contracts to: LoadSystem + UpdateConfiguration + two TransportTo singular-endpoint calls (:440-446, 474-480); the double LoadConfiguration + prescription re-add dance (:460-469) becomes unnecessary (config is not clobbered by transport) |
| FeynmanTrick/DiffExpIntegration.m:908-1110 EvaluateLimitFromTransport | replaced by EndpointLimit with "Combination" (combine-before-limit); silent zeros (:1037-1059, 1080-1093) -> E28 |
| FeynmanTrick/MatrixExport.m:42-83 / :91-114 | export side unchanged at M0; M5 switches the FT pipeline to ExportGeneralMatrix (RewritePlan I1/M5) |

Numerical lessons that MUST survive (ledger cross-references): 2x WP input
raise + $MinPrecision floors; digit budget += log10(segments); adaptive
expansion-order search semantics; two-point error probe avoiding x = 0;
per-indeterminate error tracking; abort at error > 1; coupling-depth t-order
degradation feeding TWindow (MatrixLoading.m:383, State.m:222); match point
at radius/DivisionOrder of BOTH adjacent charts (FT pins DivisionOrder = 4,
FeynmanTrick/DiffExpIntegration.m:272; classic default 3); banana classic
line needs UseMobius + UsePade + RadiusOfConvergence -> 10
(Banana_example.m:70-80) — so TransportTo must keep full Mobius support even
though IntegrateOverLine rejects it.

---

## 7. DEPENDENCIES

API.m may call: Tolerances.m, Config.m, EpsSeries.m, SectorSeries.m,
Indicial.m, Solve.m, Transport.m, Integrate.m — i.e. everything below it in
the acyclic order Tolerances < Config < EpsSeries < SectorSeries < Indicial
< Solve < Transport/Integrate < API.  Nothing calls API.m.  Typical
delegation: LoadSystem -> Config (validated keys) + Indicial-adjacent
factor extraction; PrepareBoundaryConditions -> EpsSeries (Laurent arrays) +
SectorSeries (tag parsing); TransportTo -> Transport; SolveAtPoint ->
Indicial + Solve (+ Transport for the chart frame); EndpointLimit ->
SectorSeries (combination) + Integrate (drop rule/divergence checks);
IntegrateOverLine -> SectorSeries + Integrate; ToPiecewise -> SectorSeries
(evaluation incl. Pade).  The FT layer calls ONLY DiffExp2 public symbols
(this spec) per the M0 shim contract.

---

## 8. UNIT TESTS

All tests are kernel tests (orchestrator-run).  "2F1 system" = the 2x2
companion system of Hypergeometric2F1[1/4, 1/3, 3/2, z]
(Reference/Examples/Hypergeometric2F1_example.m:26-30), exported in FULL
format dz_full.m for these tests (entries are eps-free — valid Exact input).
Closed-form references are computed in-test with Mathematica built-ins at
200 digits.

 1. test_api_loadsystem_full_format: write dz_full.m for the 2F1 matrix;
    LoadSystem returns Format -> "Exact", NumIntegrals -> 2, Variables ->
    {z}; loading the same matrix written with \[Epsilon] instead of eps
    gives the identical normalized system (compare entries exactly).
 2. test_api_loadsystem_legacy_slices: load Tests/Hypergeometric2F1_Matrices
    (slice files); assert W1 warning fired, Format -> "LegacySlices"; assert
    SolveAtPoint[..., <|z -> 0|>] raises E1 naming "SolveAtPoint" and the
    directory.
 3. test_api_loadsystem_rejects_sqrt: matrix containing Sqrt[1 - z] -> E6
    naming the factor 1 - z.
 4. test_api_loadsystem_rejects_canonical: directory with d_1.m -> E7.
 5. test_api_loadsystem_dimension_mismatch: 2x2 dz_full.m + 3x3 dw_full.m
    -> E4 naming both files and dims.
 6. test_api_bc_closed_form_banana: equal-mass banana bcs
    (Banana_example.m:31-41) on <|t -> -1/x|>; assert entry 4 parses to the
    (0,0,0) sector with eps^0 coefficient exactly 1 (eps^3 Gamma[eps]^3
    E^(3 EulerGamma eps) -> 1 + O(eps)); assert entry 3 produces sector tags
    {a,b} = {1,0}, {1,1}, {1,2}, {1,3} (from (-1/t)^(1+k eps) /. t -> -1/x);
    assert entries 1-2 are wildcard markers; assert I1 info fired for the
    zeroed sub-leading orders.
 7. test_api_bc_count_mismatch: 3 entries for the 2F1 system -> E10.
 8. test_api_bc_window_shortfall: 2 eps coefficients supplied,
    EpsilonOrder -> 4 -> E11 naming integral and needed order.
 9. test_api_transport_2f1_value: bcs at z = 1/2 from
    {Hypergeometric2F1[1/4,1/3,3/2,1/2], D[2F1]} (200 digits); TransportTo
    to <|z -> 3/4|>; assert |SeriesValues[[1,1]] −
    Hypergeometric2F1[1/4,1/3,3/2,3/4]| < 10^-20 (mirrors
    test_singular_endpoint.m:271-291); assert EndpointIsSingularity ===
    False and ErrorEstimates dims == SeriesValues dims.
10. test_api_transport_bidirectional_roundtrip: 1/2 -> 3/4 (forward), chain
    the result back to 1/2 (backward, to < FixAt); assert round-trip
    matches the starting values to 10^-18 and accumulated error
    monotonicity (invariant 9).
11. test_api_transport_chaining_error_carry: two-leg transport with
    EstimateError -> "Fast"; assert leg-2 ErrorEstimates >= leg-1
    ErrorEstimates entrywise (the old assoc-chaining dropped this; E23/
    carry semantics).
12. test_api_singular_endpoint_2f1: transport 1/2 -> 0; assert
    EndpointIsSingularity === True, LocalSolution present with sector
    exponents EXACTLY {0, -1/2} (indicial roots of 2F1 at z=0: 0 and
    1 - c = -1/2), b == 0 for both; evaluate the LocalSolution at
    z = 1/20 and assert agreement with Hypergeometric2F1 to 10^-12
    (mirrors test_singular_endpoint.m:144-198).
13. test_api_wildcard_constants_fixed: banana with "?" rows transported to
    t = -1; assert no free parameters introduced (W4 NOT fired) and values
    match the vendored oracle checkpoint (Tests/refs/, M0 task 15) to
    10^-25.
14. test_api_symbolic_indeterminates: replace one bc coefficient by symbol
    bc1; assert result is linear in bc1 (D[result, bc1] x-free and
    bc1-free), W3 fired, "PadeDisabled" -> True; substituting the true
    value reproduces test 9's numbers to 10^-18.
15. test_api_endpoint_limit_drop_rule: build (via SolveAtPoint on a
    diagonal exact system dz_full.m = DiagonalMatrix[{0, (-1+eps)/z,
    2 eps/z}]/...) a LocalSolution with sectors x^0 (constant 7/3),
    x^(-1+eps), x^(2 eps); EndpointLimit returns exactly 7/3 at eps^0 with
    b != 0 sectors dropped; flipping the x^(-1+eps) sector to b = 0 (matrix
    entry -1/z) must raise E26 naming {a,b,p} = {-1,0,0}.
16. test_api_endpoint_limit_combination: same system; "Combination" ->
    {{1, -1, 0}} where the two components have equal divergent b=0 parts
    that cancel exactly at object level; assert finite limit returned and
    that componentwise EndpointLimit on either alone raises E26
    (combine-before-limit).
17. test_api_integrate_over_line_beta: trivial 1x1 system f' = 0, f = 1 on
    [0,1]; IntegrateOverLine with Prefactor PowerAtLower -> -1+eps,
    PowerAtUpper -> -1+eps; assert the LaurentValue equals
    Series[Beta[eps, eps], {eps, 0, CompleteMax}] coefficientwise to
    10^-25: leading term 2/eps (Beta[eps,eps] = 2/eps - ... ; reference
    computed in-test by Series[Gamma[eps]^2/Gamma[2 eps], ...]).
18. test_api_integrate_requires_segments: TransportTo without
    SaveExpansions then IntegrateOverLine -> E29; with UseMobius -> True
    charts -> E30.
19. test_api_topiecewise_matches_transport: equal-mass banana, transport
    t = -1 -> 5 with SaveExpansions; ToPiecewise values at t = 3 and t = 2
    agree with direct TransportTo results to rel. 10^-10 (mirrors
    test_topiecewise.m:82-195); evaluation at t = 7 (outside) -> E35.
20. test_api_unprescribed_branch_point: equal-mass banana WITHOUT the
    t - 16 + I \[Delta] prescription, transport crossing t = 16 -> E19
    naming the chart and the offending sectors; with
    "AbortOnAnalyticContinuationFail" -> False assert warning +
    "MultivaluedFail" -> True, and chaining from it -> E25.
21. test_api_line_on_singularity: banana line <|t -> 16|> as a LINE through
    the singular surface (e.g. <|t -> 16 + 0 x|>) -> E15 naming t - 16.
22. test_api_eps_epsilon_alias_end_to_end: same transport with bcs written
    in eps vs \[Epsilon]; identical SeriesValues bitwise.
23. test_api_result_shape_invariants: for tests 9, 12, 19: assert invariant
    2 (shapes), 3 (window sanity), 7 (segment tiling) directly on the
    returned associations.

---

## 9. LINE BUDGET

Target ~250 lines (RewritePlan 3.2: "API.m (~250)"), hard ceiling within
the core total 3.5k.  This is achievable ONLY because API.m contains no
mathematics: parsing/validation/dispatch plus result assembly.  Counting
guide: LoadSystem ~70 (file probing + validations), PrepareBoundaryConditions
~60 (entry-type dispatch; tag parsing helpers live in SectorSeries),
TransportTo wrapper ~45 (input normalization, chaining, result assembly),
SolveAtPoint ~15, EndpointLimit ~20, IntegrateOverLine ~20, ToPiecewise ~20.
Config functions are counted in Config.m (~150).

If over budget, cut in this order (never cut error contracts):
  1. SaveExpansionsCompressDirectory file-hash export (old
     Transport.m:854-863): keep in-memory Compress only; waive in ledger.
  2. The old 3-element {anchor, values, errors} bcs form (sec. 2.3 form 4):
     chaining via result associations covers every live use.
  3. ToPiecewise truncation argument (ord): callers can truncate via
     Config "SaveExpansionsOrder".
  4. Move the closed-form-to-sector-tag parser entirely into SectorSeries.m
     (API keeps only the dispatch switch) and re-count.
NEVER cut: positional TransportTo compatibility, wildcards, singular-endpoint
mode, the result keys, or any E/W item of section 5.

---

## 10. OPEN QUESTIONS

1. Key naming: keep "SeriesValues"/"KinematicPoint" verbatim (this spec) vs
   introduce "Values"/"Point" with aliases.  Spec chose verbatim to make
   M4/M6 parity mechanical; the minimalism reviewer may prefer aliases ->
   decide before M1.
2. Per-master vs uniform EpsWindow in TransportResult: this spec mandates a
   uniform window (matrix shape) with per-master LaurentValues only in
   EndpointLimit/IntegrateOverLine outputs.  FT carries per-master MinPower
   today (run_ft_stepwise.m:142-148); if M5 finds genuinely ragged windows
   across masters of one level, the result needs a per-master window key —
   resolve at M5 cutover design.
3. Singular-endpoint "SeriesValues" compatibility forms: this spec keeps
   Logx-bearing x-series alongside "LocalSolution" for
   test_singular_endpoint.m-class consumers.  If the M6 battery disposition
   retargets that test to LocalSolution directly, the compat forms can be
   dropped (saves ~10 lines).
4. Canonical d_1.m rejection (E7) assumes the five-point example is the
   only user — confirm against Docs/ExportDisposition.md (M0 task 2) before
   freezing; if another consumer surfaces, E7 must point at a conversion
   path instead of Legacy.
5. Should LoadSystem also accept an in-memory matrix Association (no
   directory) for unit tests?  Tests above write temp files; an in-memory
   form would simplify them but adds surface.  Default: files only.
6. "AbortOnAnalyticContinuationFail" default: this spec follows the old
   effective default (abort/error; FT opts out).  Confirm Config.md agrees
   on the key name and default.
7. ErrorEstimate field of LocalSolution vs "ErrorEstimates" of the result:
   the result key is the per-(integral, eps-order) matrix; the LocalSolution
   field is per eps-order (RewritePlan 3.1).  The reduction (max over
   components vs per-component) at the seam needs one sentence in
   Transport.md — flag for the Transport spec agent.
