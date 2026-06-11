# DiffExp2: minimal sector-native rewrite plan (A + B heavy)

Status: v2, 2026-06-11 — adversarially reviewed (3-lens panel: mathematical
correctness, legacy numerical wisdom, execution realism; ~40 findings folded
in).  Owner directive: full hardening (A) plus sector-structured series
end-to-end (B heavy); a rewrite is approved if it keeps classic point-to-point
transport AND makes the Feynman-trick (FT) pipeline first-class with better
logic flow and less code.  This document is the execution contract: after
context compaction, implementation proceeds from this file (plus
memory/box-seg12-continuation.md and the survey/review reports referenced in
section 11).

---

## 1. Goals and non-goals

GOALS
- One minimal library (working name `DiffExp2`, context `DiffExp2``):
  (a) classic series transport of eps-expanded DE systems between regular
      points (original DiffExp use case, API-compatible where cheap), and
  (b) FT-pipeline needs natively: solve AT singular points, endpoint limits,
      integration over [0,1] across interior singular points, exact Laurent
      bookkeeping in eps.
- Sector-native representation: x^(a+b*eps) (eps Logx)^p tags carried as
  exact data, never reconstructed numerically.
- No silent fallbacks; honest completeness metadata in BOTH eps-order and
  t-order; one tolerance module.
- Core target ~3.1-3.3k lines (hard ceiling 3.5k); current DiffExp core is
  ~11k.  FT layer (~4.7k) retained with hardening (~ -1k).

NON-GOALS (v1)
- No new physics beyond the 8 examples + classic transport.
- Matrices with irrational (sqrt) x-dependence: OUT OF SCOPE v1 — loud error
  at LoadSystem.  (Old code's DEqnSquareRoots machinery is ledger-documented
  for a future v1.1; algebraic *exponents* a,b ARE supported.)
- No bug-for-bug compatibility; RESULTS parity on the pin suite only.
- Old library not deleted: frozen at f48cd94 as the parity oracle until the
  full ladder (incl. double_box_planar triage) is done.

---

## 2. Diagnosis

Empirical record: bubble, sunrise, banana, box (11 digits), box_bubble (all
orders) validated end-to-end; box_triangle NO VERDICT on the N-root run
(power loss mid-run); pentagon FAILS its pin at every order INCLUDING the
leading pole even with the working N-root fitter — its content loss happens
in the zero-regulator drop rule on collapsed towers, upstream of any fitting,
PLUS it emits "current point is not recognized as a branch point... add
DeltaPrescriptions" (a CONFIGURATION gap that no rewrite fixes by itself).
Every new structural configuration broke something; every break was silent.

Disease classes:
D1 SILENT DEGRADATION (shape fallbacks, zero-fills, IntegerQ-on-floats,
   silent drops, type poisoning, zero-padded windows).
D2 EXACT INFORMATION DISCARDED THEN RECONSTRUCTED (sector structure is
   indicial data; per-eps-order towers collapse it; Prony fitter + salvage +
   drop rules exist only to recover it — pentagon proves even the fixed
   fitter cannot, because the drops precede it).
D3 SCATTERED AD-HOC TOLERANCES.

Architectural insights (as amended by review):
I1 THE SECTOR SPECTRUM IS EXACT ALGEBRA — WITH A VERIFIABLE CONTRACT.
   At a singular chart the exponents are eigenvalues of the eps-dependent
   residue A_{-1}(eps).  Caution (reviewer counterexample [[0,1],[eps,0]] →
   ±sqrt(eps)): linear-in-eps residue does NOT imply affine eigenvalues.
   CONTRACT: Indicial.m requires the characteristic polynomial of
   A_{-1}(eps) to factor exactly as Π_i (λ − a_i − b_i·eps) over
   rationals/algebraics; LOUD ERROR otherwise (risk R1).  All 8 example
   families satisfy this (campaign-verified for box L1/L2; M2 verifies the
   rest).  PREREQUISITE: this needs the EXACT eps-rational matrix — the
   current per-order slice exports (d<var>_k.m, eps-TRUNCATED at export
   order) cannot certify it.  DiffExp2 LoadSystem consumes the exact full
   export (d<var>_full.m / ExportGeneralMatrix); slice dirs allowed only for
   legacy parity transport.  FT cutover (M5) switches MatrixExport to the
   full format.
I2 SYMBOLIC-eps FROBENIUS REMOVES THE RESONANCE INDUSTRY — WITH TWO NAMED
   RESIDUES OF IT.  Per-sector recursion denominators are
   (n + a_i − a_j + (b_i − b_j) eps).  (i) TRUE resonance (same b,
   integer-spaced a): explicit log-chains (confluent p-families), exact.
   (ii) PSEUDO-resonance (integer-spaced a, DIFFERENT b): denominator =
   (b_i−b_j)·eps exactly — a Laurent SHIFT.  Spec decision: solve
   integer-spaced-a families JOINTLY, choosing the eps-regular combination
   at each collision (log-chain-analogue construction keeping coefficients
   finite at eps=0); never silently lose window.  Campaign reality: banana
   segment-1 endpoint mixes x^0, x^(−1+eps), x^(2eps) — pseudo-resonance is
   in-scope from day one.  What this deletes: the strategy stack
   DiffExp/IntegrationStrategies/{Recurrence 870, ResonantRecurrence 1527,
   VOP 304, Default 197, Dispatch 132, Helpers 124, IntegrationStrategies 43}
   + Wronskian.m 112 + Frobenius.m 141 = ~3.45k lines, plus LocalSeries.m
   (1026: the finite-width machinery is subsumed; FuchsianizeLocal is ported
   — see 3.2 Indicial).  Replaced by ONE solver (~700-800 incl. the ported
   rank reduction).

---

## 3. Core design

### 3.1 The fundamental object

LocalSolution at expansion point x0 (chart coordinate t; chart map recorded):

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

  Sector = <|
    "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
    "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT:
                solutions are vectors; matching and ODE residuals need all
                components *)
  |>
  Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                  Σ_k eps^k Σ_n c[k,n,comp] t^n

INVARIANTS (always on, cheap):
- Regular chart ⇔ exactly one sector (a=0,b=0,p=0): classic transport is the
  degenerate case, same code path.
- Tags exact; coefficients are the only numerics.
- EpsWindow/TWindow are the truth; no padding; consumers check need <=
  CompleteMax and fail loudly naming (chart, sector, order).  EXEMPTION:
  homogeneous true-resonance sectors have kmin = −p BY CONSTRUCTION
  (log-chain weights are eps-independent under the (eps Logx)^p/p!
  normalization); window arithmetic must account, not error.
- ODE residual spot-check at one random interior point per chart.

### 3.2 Module map

  DiffExp2/
    Tolerances.m   (~100)  ALL thresholds, derived from WorkingPrecision,
                           named semantics (chopFloor, matchTol, snapTol,
                           rankTol, laurentLeadTol — the FT LaurentTrim
                           10^-40 ABSOLUTE test becomes a relative one here).
    Config.m       (~150)  option schema + validated accessor (the FEC-class
                           silent-miss bug is impossible by construction);
                           variable-context pinning; the KEPT config surface
                           (table in Docs/specs/Config.md per M0; live keys
                           enumerated by the legacy review: AccuracyGoal[+
                           Validate], ChopPrecision sync, DeltaPrescriptions,
                           DivisionOrder, EpsilonOrder, EstimateError,
                           ExpansionOrder, LineParameter, MatrixDirectory,
                           RadiusOfConvergence, RationalizationTolerance,
                           SegmentationStrategy(Predivision only v1),
                           UseMobius, UsePade, Variables, Verbosity,
                           WorkingPrecision, Crosscheck*, SaveExpansions*,
                           AbortOnAnalyticContinuationFail; each kept key has
                           a DiffExp2 meaning or a waiver).
    EpsSeries.m    (~300)  truncated eps-LAURENT arrays: add/mul and
                           LAURENT DIVISION with explicit window-shift
                           semantics (denominator with vanishing eps^0 part
                           shifts kmin AND CompleteMax — the 1/(b eps)
                           enhancement is a shift, not an inversion); honest
                           window arithmetic for all ops.
    SectorSeries.m (~400)  LocalSolution/Sector algebra: evaluate (with
                           branch rule), multiply by rational c(x,eps)
                           (closed; partial fractions split poles ACROSS
                           charts — absorbs old SingularityDecomposition.m's
                           role; the FT call site gets a named replacement
                           API), re-expand around a new center (truncation
                           contract: target order + tail bound (Δ/R)^N feeding
                           ErrorEstimate), differentiate.  Pade evaluation
                           DECISION: ported here as the evaluation
                           accelerator (~+100 lines, loud fallback — old
                           GetPade warns silently); the R6 benchmark compares
                           old-with-Pade vs new-with-Pade.
    Indicial.m     (~300)  exact residue from the FULL eps-rational matrix;
                           char-poly factorization contract (I1); Jordan/
                           confluence -> {a,b,p} specs; resonance/pseudo-
                           resonance partitioning; HIGHER-ORDER POLES:
                           ported rank-reduction (Moser/shearing a la
                           FuchsianizeLocal — banana L1 endpoints have a
                           NILPOTENT DOUBLE POLE, in-scope day one) with
                           loud non-termination error.
    Solve.m        (~700)  ONE solver: Frobenius per sector family at
                           symbolic eps (eps-Laurent coefficient arithmetic);
                           x-denominators of A cleared up front so recursion
                           coefficients are polynomial (the old rational-
                           recurrence fast path's lesson — this is the
                           performance backbone); true-resonance log-chains
                           explicit; pseudo-resonant families solved jointly
                           (I2 spec); inhomogeneous sources in the same
                           representation with the RESONANT-SOURCE LOG-BUMP
                           rule (source tag hits homogeneous tag exactly ->
                           p -> p+1 ansatz; the old "general solver returns
                           empty particular" hole gets a closed-form pin).
    Transport.m    (~600)  segmentation (exact singularity solve incl.
                           complex roots; radius = complex distance; chart
                           sizing with the GetCPL/GetCPR lesson: match point
                           at radius/DivisionOrder of BOTH adjacent charts);
                           predivision two-pass with digit budgeting
                           (DigitsNeeded += log10(#segments) + safety);
                           input precision raised to 2x WP, $MinPrecision
                           floors (ledger); marching loop with eps-GRADED
                           Laurent weight-matching: assert ord_eps det F = 0
                           or recombine basis unimodularly at eps=0
                           ((S_i−S_j)/((b_i−b_j)eps) columns) — regular
                           incoming data CAN have 1/eps weights in a naive
                           sector basis (log x = (x^(2eps)−1)/(2eps) class);
                           crossing operator = phase e^(±iπ(a+b eps)) TIMES
                           the unipotent log-chain mixing matrix
                           M_{p→p−j} = (±iπ eps)^j / j! (p>0 sectors populate
                           lower-p members), convention pinned to the
                           old principal-branch-far-side form (Logx →
                           Logx − 2πi θm class; 9aeb300 interior-split
                           lesson); two-point error probe per segment.
    Integrate.m    (~500)  per-sector exact antiderivatives with the FULL
                           case table {b=0 / b≠0} × {a+n+1 <0 / =0 / >0} ×
                           {p=0 / p>0}: denominators (a+n+1+b eps)^(j+1),
                           j=0..p (pole depth p+1, NOT a single 1/(b eps));
                           dimreg convention t^(a+n+1+b eps)|_{t=0} := 0
                           for b≠0 (stated, it IS the drop rule's integral
                           counterpart); b=0 divergence error for a+n+1 <= 0
                           (any p) UNLESS cancelled in the assembled
                           combination (checked at the object level);
                           interior-pole PV/iδ with the pairing of the two
                           half-segments enforced (assert, not assumed);
                           Mobius charts REJECTED loudly v1 (old FT requires
                           UseMobius->False for integration; same contract).
    API.m          (~250)  LoadSystem (full-format matrix dir, closed form,
                           legacy slice dir for parity transport only),
                           TransportTo (assoc chaining, BIDIRECTIONAL,
                           singular-endpoint mode returning the series,
                           "?" per-integral wildcard BCs, closed-form eps
                           expressions with auto-expansion, symbolic
                           indeterminate coefficients), SolveAtPoint,
                           EndpointLimit, IntegrateOverLine, SaveExpansions/
                           ToPiecewise equivalents, ErrorEstimates output,
                           both `eps` and `\[Epsilon]` accepted.  The full
                           compatibility contract = Docs/specs/API.md (M0),
                           driven by Reference/Examples + Tests usage.
  Core total target: ~3.3k (ceiling 3.5k).

RETAINED: FeynmanTrick/ (hardened, consumes DiffExp2 API; shim audit per
R4), Scripts/, Tests/, comparator + pySecDec infra.  FROZEN: DiffExp/ (old)
as oracle; open old-core defects are NOT debugged further — converted to
pins/ledger entries (R8).

### 3.3 FT boundary cases (unchanged from v1 of this plan)

limitUpper/limitLower: constant of the (0,0,0)-sector; b≠0 dropped exactly;
divergence loud.  integrate: rational-multiply + exact cancellation at
object level + Integrate.m closed forms.  direct: Laurent convolution with
honest windows.

### 3.4 Static order budget (amended per review)

  need(L0) = requested
  need(Lk) = need(L(k-1)) + prefactorShift(Lk)              [cumulative
             convention as in old code]
           + Max[0, -minEpsPower(IBP coeffs at Lk)]
           + enhancement(Lk)        [= max over endpoint pole sectors
             (a+n+1=0 reachable, b≠0) of (p+1); NOT a count — several p=0
             sectors give depth 1; one p=2 sector gives 3]
           + matchingShift(Lk)      [ord_eps det F at singular charts if
             not removed by basis recombination; target 0, asserted]
           + pseudoResonanceShift(Lk) [0 under the joint-solve spec;
             asserted]
           + 1 safety
EpsWindow propagation turns any miscount into a named error.  Diagnostic
env overrides retained.  Validation: the formula must reproduce the
campaign's measured needs (box_bubble: 9 at L2, 11 at L1) — an M5 unit.

---

## 4. Stage A hardening (retained FT layer; new core by construction)

A1 No silent fallbacks (survey-cataloged sites; the full list lives in the
   three survey reports + legacy-review finding 13: FEC-key bug, prefactor
   zero-fills -> validated accessor, reduction-residual assert, negative-
   index requests loud, EvaluateLimitFromTransport combine-before-limit and
   its Quiet[Check[...,0]] chains, expOrd=30, phantom-order laundering,
   ShiftRawBoundaries zero-padding -> CompleteMaxPower metadata, LaurentTrim
   absolute->relative tolerance, "keep at least one output order" floor
   deleted).
A2 Always-on invariants: reduction residual; Euclidean-reality flag; ODE
   spot-checks; structural type asserts.
A3 ReductionCache ON by default (dbox iteration cost).
A4 Battery + ladder as the only merge gate; add the four survey-sketched
   tests + per-module unit suites (M1-M4).

---

## 5. Lessons Ledger (M0 deliverable; seeds — reviewer-amended)

Docs/LessonsLedger.md: every item implemented-or-waived with reason.  Seeds:
input precision raise to 2x WP + $MinPrecision/$MaxExtraPrecision floors;
ChopPrecision/LinearSolveChopPrecision sync; expansion-order adaptive search
(AccuracyGoalValidate Before/After semantics); DivisionOrder = match point
at radius/k of BOTH adjacent charts (GetCPL/GetCPR), FT pins k=4; predivision
segment-count digit budgeting; coupling-depth t-order degradation
(MaxCouplingOrder, ISafety, ICurrEvalErrorSeriesDecrease formulas);
two-point error probe avoiding x=0 (log-vanishing trap), per-indeterminate
errors, abort>1; RoC chart rescaling keeps high-order coefficients O(1)
(banana REQUIRES RoC=10); Pade evaluation rationale + silent-fallback fix;
Mobius charts (which examples need them; banana classic line does);
complex-singularity ghost projection (Re, Re±Im, suppression rule);
delta-prescription derivation (multiplicity parity, leading-coeff sign,
conflicts, sqrt auto-prescriptions, SingularityCheck default-on);
interior-split real-log convention (9aeb300); principal-branch crossing
convention (sign +1 = NO replacement; sign −1 = Logx − 2πi θm SHIFT that
mixes log-chains binomially; e^(−2πi b) only for denominator > 2);
denominator-clearing recursion fast mode; segment block caching;
eps-window trimming semantics incl. uniform per-level shift; FIRE
variable-context pinning; numerical-zero leading-coefficient skipping
(generalizes to matching solves); the resonant-source empty-particular hole.

---

## 6. Milestones

M0 FOUNDATIONS (parallel agents + IDLE KERNEL pre-generates oracles)
   Agent tasks (one per agent, prompts enumerated by the execution review):
   (1) LessonsLedger.md (exhaustive, file:line; seeds above);
   (2) Docs/ExportDisposition.md: EVERY ::usage export of DiffExp/*.m (~60)
       + State.m (~90) classified kept-in-API / absorbed-into-<module> /
       dropped-with-reason; M5 may not start with FT-referenced symbols
       unclassified (known orphans: GetPade/SEval*, ToPiecewise,
       IntegrateSystem, PrepareBoundaryConditions, DecomposeSingularity*,
       FindMatrixSingularities, sqrt machinery, CrosscheckFlags, the 15
       RegularizedIntegration exports);
   (3) FT reach-in shim contract (State x27, Symbols x10, Utilities x5,
       Transport x2, RegularizedIntegration x2, SingularityDecomposition x1);
   (4-11) Docs/specs/<Module>.md for the 8 modules (signatures, invariants,
       error contracts per section 3);
   (12-13) adversarial spec reviews (math lens; minimalism/error-contract
       lens);
   (14) Tests/PINS.md + VENDOR the pin generators and frozen reference JSONs
       from /tmp into Tests/refs/ (power-loss lesson: /tmp is not storage);
   (15) Scripts/dump_transport_checkpoints.m (old TransportTo -> JSON values
       at every segment boundary; the missing parity harness);
   (16) pentagon triage: root-cause the unrecognized-branch-point warning
       against FTExamples.m DeltaPrescriptions; re-verify the pentagon pin
       sign via the box1l sign-protocol (R9).
   KERNEL (parallel): pre-generate oracle artifacts — checkpoint dumps for
   bubble/sunrise/2F1/banana lines, dispatch dumps, per-example STEPWISE
   logs — so M3/M4 gates need zero old-code kernel time.
   GATE: specs reviewed; ledger reviewed; disposition table complete;
   oracle artifacts on disk (in-repo).

M1 Tolerances + Config + EpsSeries (+units: Laurent division window-shift
   semantics vs exact rationals at 300 digits).  GATE: units green.

M2 Indicial + SectorSeries (+units).  Indicial vs known systems: box L1
   (b={-1,+1} class), box L2 apparent chart — CONCRETE pin: compute the
   indicials from the stored campaign matrices and assert the recorded
   beta-root chart values (t*=7/11, eigenvalue −1; memory brief) —
   banana L1 endpoint (nilpotent double pole through rank reduction;
   5-fold resonant residue), 2F1.  GATE: units green.

M3 Solve + single-chart validation: regular-point parity vs old (bubble
   segment, 1e-30, via dispatch-dump replay); 2F1 closed forms; box L2
   apparent chart (campaign J1/J2 pins to 20 digits); banana L1 endpoint
   (double pole + log chains); pseudo-resonance unit (x^0 vs x^(-1+eps) vs
   x^(2eps) family, banana-derived); resonant-source log-bump closed form
   (∫ x^(-1) log^k x class).  GATE: closed-form suite green.

M4 Transport + Integrate: classic parity on bubble/sunrise/2F1/banana lines
   vs M0 checkpoint dumps (1e-25; old config pinned per example, e.g.
   banana = UnequalMassConfiguration with UsePade, UseMobius, RoC=10 —
   parity examples run WITHOUT Mobius where Integrate is involved);
   multisector closed forms (now exact, no fit); interior-pole crossing
   (box L1 pole at 1/4, campaign pin); near-endpoint pole (box_triangle L3
   geometry 0.9617); benchmark gate: banana-line wall-clock within 2x of
   old-with-Pade.  GATE: parity + closed forms + benchmark.

M5 FT CUTOVER: MatrixExport -> full format; FT consumes DiffExp2 API (shim
   audit per disposition table); static budget ladder; EpsWindow end-to-end.
   Ladder with incremental commits per passing example, long runs in
   background with dump-on-failure: bubble, sunrise, banana, box,
   box_bubble, box_triangle, pentagon, double_box_planar — each vs pins
   (budget-formula unit: reproduce box_bubble's 9/11).
   GATE: >= 6/8 with MANDATORY triage of any failure into D2-class
   (rewrite should have fixed — blocker), prescription/config-class
   (FTExamples fix), or pin-class (R9 protocol).  dbox is stretch (R5).

M6 CLEANUP: Stage-A hardening edits to retained FT files (survey batches);
   battery DISPOSITION: every test retargeted-to-DiffExp2 / kept-on-Legacy
   (on-demand oracle only) / deleted (test_multisector_fit's fitter gates
   die with the old core; its closed forms were ported in M4); old DiffExp
   -> Legacy/ with a sweep of external dependents (README, AGENTS.md,
   Reference/Examples, Papers paths) (R10); docs rewrite; dead-code
   deletions; final tag.
   GATE: battery (retargeted) + full ladder green against the new core.

Effort estimate (amended): M0 one parallel-agent day + idle-kernel oracle
generation; M1-M4 2-4 sessions; M5 2-4 sessions (kernel-bound: full ladder
up to ~2.5h/pass, several passes expected); M6 one session.

---

## 7. Agent execution model

ORCHESTRATOR owns the single kernel (all validation serial), merges,
commits.  REVIEW agents: M0 ledger/specs/dispositions; per-milestone
adversarial diff review BEFORE kernel time.  IMPLEMENTATION agents: one
module per agent, code + unit tests, never run wolframscript; orchestrator
runs tests and hands failures back (SendMessage) or fixes small breaks
directly (campaign-calibrated: the 553-line N-root delivery needed exactly
one orchestrator-fixed gate bug).  Kernel discipline: pre-generate oracles
in M0; batch gates; long runs in background with dumps.

## 8. Risk register

R1 Indicial contract violation (eigenvalues not affine in eps): loud error;
   out-of-scope perturbative fallback documented.
R2 Algebraic a/b: exact-symbolic supported; performance flagged not
   correctness.
R3 Ill-conditioned weight matching near chart edges: GetCPL/GetCPR geometry
   + residual assert.
R4 FT hidden reach-ins: disposition table + shim audit gate M5.
R5 dbox genuinely new structure: triage as capability work, not failure.
R6 Performance: eps-Laurent ARRAY arithmetic + cleared-denominator
   recursion + Pade port; M4 benchmark gate vs old-with-Pade.
R7 Single-kernel throttling: accepted; mitigated by M0 oracle pre-generation
   and background ladder runs.
R8 In-flight old-core defect work freezes at f48cd94; open defects become
   pins/ledger entries only (pentagon drop-rule, box_triangle near-endpoint
   leak, seg12 thread); M5 failures must be attributable via the M0
   pentagon-triage and pin-protocol work.
R9 Pin correctness (loop_package sign conventions): the box1l sign protocol
   re-run for pentagon/box_triangle pins BEFORE M5 judgments.
R10 External dependents on old API/paths: M6 sweep item.
R11 Non-Fuchsian charts beyond the ported rank reduction: loud error naming
   the chart; evidence so far: banana L1 double pole is the worst in-scope
   case and is covered.

## 9. Pin suite

As v1 of this plan, PLUS: all pins/generators vendored in-repo (Tests/refs/,
Tests/PINS.md) at M0; box L2 apparent-chart indicial + pointwise values from
the campaign memory; banana L1 endpoint closed forms; pseudo-resonance and
resonant-source closed-form units; pentagon pin held PROVISIONAL until the
R9 sign protocol passes; box_triangle loop_package pin likewise.

## 10. Baseline at plan time

master = f48cd94.  Validated: bubble, sunrise, banana, box (11 digits),
box_bubble (all orders).  box_triangle: NO VERDICT (N-root run lost to
power-cut; pre-N-root runs were complex-corrupt from the L3→L2 step,
interior pole at 0.9617 near endpoint).  pentagon: FAILS pin at all orders
incl. leading pole with the working N-root fitter; failure sits in the
zero-regulator drop rule on collapsed towers (D2) AND an unresolved
DeltaPrescriptions configuration warning (M0 task 16 triages).  N-root
fitter itself: 12/12 units, battery 18/18, no regressions — committed as
the best-possible patch on the old representation and superseded by this
plan at M5.  Env-gated debug instrumentation intentionally retained.

## 11. Referenced artifacts

memory/box-seg12-continuation.md (campaign mechanisms, pins, debug envs);
the three survey reports (FT-layer refactor batches; DiffExp dead code;
Scripts/Tests/docs hygiene incl. new-test sketches) and the three plan
reviews (math/legacy/execution) — task outputs under the session transcript
dir; full review text preserved at
/private/tmp/claude-501/.../tasks/wl6pozq2m.output (vendor key findings into
Docs/specs/ during M0 — /tmp is not storage).
