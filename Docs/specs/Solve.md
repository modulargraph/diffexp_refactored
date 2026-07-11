# Module spec: DiffExp2/Solve.m (~700 lines; restated ~725 at M0 — section 9)

Status: M0 deliverable (RewritePlan task 4-11).  This spec is a contract: an
implementation agent reading ONLY this file, Docs/RewritePlan.md, and the
cited old code must be able to write Solve.m and its unit tests.  Context:
`DiffExp2`Solve``.  All old-code citations are against the frozen oracle at
f48cd94 (working tree paths under DiffExp/).

---

## 1. PURPOSE

Solve.m is THE one local solver of DiffExp2: given a chart system assembled
by its own PrepareChart (section 2.5, DEC-7) from Indicial.m's exact output
(theta-form matrix, exact in eps, at worst a simple pole at the chart origin
after rank reduction, plus the exact spectral data of the residue), it
constructs (a) the homogeneous fundamental system — exactly one
LocalSolution column per indicial sector spec — and (b) particular solutions
for inhomogeneous sources given in the same sector-native representation, by
a single Frobenius recursion at SYMBOLIC eps whose coefficients are truncated
eps-Laurent arrays (EpsSeries.m).  The x-denominators of the matrix are
cleared up front so the recursion has finite width in the chart variable (the
old rational-recurrence fast path's lesson — the performance backbone, R6);
the only divisions anywhere are EpsSeries divisions by indicial offsets
(n + a_i − a_j + (b_i − b_j)·eps), so true resonance (identically-zero
offset) becomes an explicit log-chain, pseudo-resonance (offset = Δb·eps)
becomes an exact Laurent shift handled by the JOINT family construction, and
the resonant-source log-bump (p → p+1) closes the old "empty particular"
hole.  Solve.m replaces the entire old strategy stack
(Recurrence/ResonantRecurrence/VOP/Default/Dispatch/Helpers + Frobenius +
Wronskian + the finite-width core of LocalSeries) with one code path; the
regular chart is the degenerate case (single sector (0,0,0), all divisions
Taylor), not a separate strategy.

---

## 2. PUBLIC SYMBOLS

Exactly five exported symbols (M0 review: SolveChart added per
REVIEW-minimalism defect 3, ratified DEC-9; PrepareChart added per
REVIEW-math D6, ratified DEC-7).  Everything else is `Private`` (unit tests
may reach into the private context; nothing else may).

### 2.1 SolveHomogeneous

```
SolveHomogeneous[chart_Association, req_Association] -> fundamental_Association
```

- `chart`: a ChartSystem (section 3.2).
- `req`: <|"EpsWindow" -> <|"Min" -> kmin_Integer, "CompleteMax" -> kmax_Integer|>,
          "TOrder" -> nmax_Integer (>= 0)|>
  The eps-window and t-order the CALLER needs.  Solve widens its internal
  work window (section 3.6) so that the DELIVERED windows meet the request;
  this is possible for homogeneous solutions because the matrix data is
  exact in eps.
- Returns a FundamentalSystem (section 3.4): one LocalSolution column per
  indicial sector spec, ordered (family, descending a, chain position),
  normalized per section 3.5, with Diagnostics.
- Side effect: runs the ODE-residual spot check (section 2.3) on every
  column before returning; a failed check is a loud error, never a warning.
- Idempotent.  The exact prepared data (cleared polynomials, spectral
  frame V/VInv/J) lives in the ChartSystem assembled ONCE per chart by
  PrepareChart (2.5); the eps-expanded work-window arrays are computed per
  solve call.  The old memo-table plumbing is CUT (§9 cut 3, taken at M0
  review to fund SolveChart — REVIEW-minimalism defects 3-4): on the
  Transport path SolveChart performs the homogeneous and particular solves
  in one call, so nothing is recomputed there; results must be
  bit-identical regardless of call pattern.

### 2.2 SolveParticular

```
SolveParticular[chart_Association, source_Association, req_Association]
  -> LocalSolution
```

- `chart`: same ChartSystem.  `req`: same shape as above.
- `source`: a SourceData (section 3.3) — sector-native, THETA-form
  (section 3.3 normalization), with honest EpsWindow/TWindow.
- Returns ONE LocalSolution: the particular solution with the canonical
  kernel choice (all homogeneous freedom set to zero, section 3.5) and with
  pseudo-resonant compensation terms included (section 3.7), windows honest:
  EpsWindow.CompleteMax = source CompleteMax minus the recorded
  enhancement shifts (resonant 1/eps hits and pseudo-resonant collisions on
  the solved path), TWindow.CompleteMax = Min[req TOrder, source
  TWindow.CompleteMax].  If the delivered window is smaller than req, that
  is NOT an error here — the honest window is returned and the CONSUMER
  (Transport/Integrate) fails loudly on need > CompleteMax (RewritePlan 3.1
  invariant); but an EMPTY delivered window (CompleteMax < Min) IS a loud
  error at this level (E6).
- A structurally zero source (empty "Sectors" list) returns the zero
  particular with FULL requested windows (old behavior at
  Recurrence.m:318-322, 846-850, made structural instead of PChop-numeric).
- Side effect: ODE-residual spot check with the source term included.

### 2.3 ODEResidualCheck

```
ODEResidualCheck[chart_Association, sol_Association, source_:None,
                 probe_:Automatic] -> resid_Real
```

- Evaluates r(t0) = theta f(t0) − B(t0,eps)·f(t0) − s(t0) at one interior
  probe point t0 (Automatic: a random real t0 in (Radius/4, 3·Radius/4),
  drawn once per chart with a fixed seed derived from the chart center so
  runs are reproducible), for every delivered eps-order in the window, on
  the principal sheet (t0 > 0 real; DiffExp2 carries NO theta/branch
  symbols — coefficients are complex numbers, so a single evaluation
  suffices, unlike the old both-branches workaround at
  ResonantRecurrence.m:1003-1015).
- Returns the max relative residual (relative to the evaluated solution
  scale at the same eps-order; scale floor 1).  ERRORS loudly (E7) when the
  residual exceeds `residualTol` from Tolerances.m.  Exported so Transport
  can re-check after re-expansion; called internally by 2.1/2.2 (always-on
  invariant, RewritePlan 3.1).
- `sol` may be a single LocalSolution or a FundamentalSystem (checks every
  column, homogeneous: source = None).

### 2.4 SolveChart

```
SolveChart[cs_Association, req_Association, source_:None]
  -> <|"Basis" -> FundamentalSystem, "Particular" -> LocalSolution | None,
       "CouplingDepth" -> _Integer,
       "IntegrationSequence" -> _Association (* multi-SCC only *)|>
```

The Transport-facing wrapper (`cs` = a ChartSystem from PrepareChart, 2.5)
implements the old `InitializeIntegrationSequence` role as an exact,
certified SCC decomposition.  The graph is built from `ThetaOriginal`, the
original physical theta-form matrix `t beta A(x0 + beta t)`, BEFORE any
full-system rank-reduction gauge or spectral frame can mix components.  An
equation-row/component-column entry gives the dependency edge `c -> r` iff
`Cancel[Together[ThetaOriginal[[r,c]]]]` is exactly nonzero.  The input must
be exact.  In particular, eps-only and analytic-regulator-only entries remain
edges: eps and all analytic regulators are left symbolic and there is no
zero-point, generic-point, or other regulator specialization.

The resulting certificate has the implemented schema

```
<|"Schema" -> "DiffExp2.IntegrationSequence/v1",
  "Exact" -> True, "RegulatorSpecialization" -> None,
  "MatrixHash" -> ..., "PatternHash" -> ...,
  "NonzeroPattern" -> {{0|1, ...}, ...},
  "DependencyEdges" -> {{sourceComponent, targetEquation}, ...},
  "Components" -> {{masterIndex, ...}, ...},
  "VertexToComponent" -> {...},
  "CondensationEdges" -> {{sourceBlock, targetBlock}, ...},
  "TopologicalOrder" -> Range[numberOfBlocks],
  "BlockDepths" -> {...}, "CouplingDepth" -> _Integer,
  "Certificate" -> <|"Partition" -> True,
    "StrongConnectivity" -> True,
    "AcyclicTopologicalOrder" -> True,
    "EdgeCoverage" -> True|>|>
```

Every SCC's master indices are sorted.  SCCs are renumbered in a
deterministic source-first topological order, with the smallest master index
as the ready-queue tie break.  Construction rechecks the partition, forward
and reverse reachability inside every SCC, topological order of every
condensation edge, and coverage of every original dependency edge.
`CouplingDepth` counts blocks on the longest condensation-DAG path.  The
one-SCC direct path reports `CouplingDepth -> 0`; a multi-SCC result returns
the certificate as `IntegrationSequence`.

Before executing a multi-SCC plan, Solve prepares the diagonal block systems
and audits every block which is the TARGET of a condensation edge.  At the
guarded work order

```
workTOrder = req["TOrder"] + 2 + 2 CouplingDepth
```

it evaluates the exact `recurrencePoleDepth` of each target's cleared
finite-width recurrence.  If all recorded depths are zero, execution remains
blockwise.  If any target has positive depth, propagating an upstream source
would repeatedly consume future eps orders, so the whole chart is coarsened:
the lazy envelope fully prepares the registered original physical system at
the same affine chart with `UseSCCSkeleton -> False`, then runs one monolithic
strict recurrence.  This is the same configured recurrence solver acting on
the same exact system, not an alternate-solver fallback.  The exact SCC
certificate remains the planner metadata even though its blocks are not used
as the execution units.

The coarsened basis and particular diagnostics contain

```
<|"SCCSolved" -> False,
  "SCCExecutionMode" -> "MonolithicStrict",
  "SCCExecutionReason" -> "DownstreamRecurrencePoleDepth",
  "SCCExecutionPlan" -> <|
    "Mode" -> "MonolithicStrict",
    "Reason" -> "DownstreamRecurrencePoleDepth",
    "WorkTOrder" -> workTOrder,
    "OffendingTargets" ->
      {<|"TargetBlock" -> _Integer, "IncomingBlocks" -> {...},
          "Indices" -> {...},
          "RecurrencePoleDepth" -> _Integer|>, ...}|>,
  "SCCCertificate" -> IntegrationSequence,
  "CouplingDepth" -> IntegrationSequence["CouplingDepth"]|>
```

The monolithic basis additionally carries an explicit
`"SCCRecombineGroups"` list derived from the ACTUAL lazily prepared full
ChartSystem and its returned basis specs.  Transport never tries to infer
global epsilon-degeneracy data from the deliberately incomplete SCC envelope.

Full-system preparation is a performance candidate, not a new capability
requirement.  Solve catches only a Failure raised by that preparation attempt.
It then runs the strict block path and records

```
<|"SCCExecutionMode" -> "BlockSequentialStrict",
  "SCCExecutionReason" -> "CoarseningDeclined",
  "SCCCoarseningCandidate" -> theMonolithicPlan,
  "CoarseningDeclined" -> thePreparationFailure, ...|>
```

A subsequent block preparation or solve Failure is not caught or converted.
Thus the heuristic cannot make a block-solvable chart into a capability
regression.

The monolithic core retains the caller's public window caps and the mandatory
full original-system residual check.

For the depth-zero/blockwise execution mode, the following is the implemented
contract:

1. Each diagonal block is prepared from the corresponding principal
   submatrix of the registered ORIGINAL physical system, using its original
   variable and the same affine chart.  It is not reconstructed as a
   synthetic `ThetaOriginal/t` system.  Consequently each block receives its
   own complete rank reduction and spectral preparation while retaining a
   stable physical-submatrix clearing identity across chart centers.
2. Each diagonal homogeneous system is solved strictly.  For every basis
   column, downstream responses are then solved in source-first block order.
   The source from block `u` to block `v` is the exact original-frame product
   `ThetaOriginal[[component[v], component[u]]] . solution[u]`, formed by
   `SectorSeries` rational multiplication.  Fractional powers, log sectors,
   eps windows, prescriptions, and center-pole tag shifts are therefore
   preserved.  No coefficient is sampled to decide that a coupling vanishes.
3. An external theta-form source is projected by physical component, combined
   with the incoming off-diagonal sources, and solved in the same order.  SCC
   orchestration never pre-applies a block gauge.  `SolveParticular` applies
   the target block's `T^-1` exactly once at its public boundary and gauges
   the answer back once, as required by the rank-reduced SourceData contract
   in 3.3.
4. Intermediate solves use deterministic work halos.  Their Taylor order is
   the guarded `workTOrder` above.  The initial upper-eps halo is the
   longest exact DAG-path cost, including eps-pole valuations of the
   off-diagonal coupling, spectral transforms, rank-reduction gauges,
   finite-width recurrence denominators, target dimension, and collision
   depth.  An external particular also funds its first block solve.  Honest
   delivered-window deficits trigger bounded monotone widening on the same
   strict recurrence path while unused certified source halo remains; no
   missing frame is treated as zero.  A finite external source that cannot
   fund the requested top returns its nonempty honest lower window with a
   `"RequestedWindowNotReached"` diagnostic (requested, delivered, source,
   and work tops plus attempt count).  Only an empty delivered window is E6.
5. The persistent homogeneous cache contains only the assembled parent SCC
   basis.  Diagonal-block cache entries live in a dynamically scoped scratch
   cache and are discarded after assembly, so one multi-block chart consumes
   one persistent slot.  In blockwise mode with `RecurrenceBackend -> Cpp`,
   diagonal homogeneous requests are submitted as one native batch only when
   the grouped spectral transform is enabled, no analytic-regulator variables
   are configured, and the block count fits `HomogeneousCacheCapacity[]`.
   The batch additionally enforces one parent solve tag and exact one-time
   response consumption.  If a guard is false, the same configured backend
   runs the block solves sequentially; it does not select a fallback
   recurrence solver.  `MonolithicStrict` mode invokes that configured backend on the
   single lazily full-prepared chart.
6. The assembled solutions are restored to original component order, then
   capped to the caller's eps `CompleteMax` and public `TOrder`.  Only after
   that cap does `ODEResidualCheck` verify the FULL original theta-form system
   (and the full external source, when present).  Thus diagonal-block checks
   are not accepted as a substitute for the final recombination proof.

Successful blockwise basis and particular results explicitly record
`"SCCExecutionMode" -> "BlockSequentialStrict"`.  In that mode, `Basis`
contains all downstream particular responses needed
to make each seeded block column solve the full homogeneous system.  In
`MonolithicStrict` mode it is the ordinary full-system strict basis.  `Particular` is
`None` exactly when no external `source` was supplied; otherwise it is the
canonical particular produced by the selected execution mode, including the
structurally zero case.

### 2.5 PrepareChart

```
PrepareChart[sys_Association, chart_Association] -> ChartSystem
```

(DEC-7; REVIEW-math D6.)  Applies the exact affine chart to the loaded
matrix and forms `ThetaOriginal = t beta A(x0 + beta t)`.  It always builds
and certifies the `IntegrationSequence` of 2.4 at this point, before any
global indicial reduction.  With `"UseSCCSkeleton" -> True` and more than
one SCC, PrepareChart returns the lazy SCC envelope of 3.2 without calling a
full-system ``Indicial`ChartIndicial``; each diagonal physical submatrix is
fully prepared only when the SCC solver needs it.  Transport requests this
mode.  The execution audit of 2.4 may later full-prepare that original system
lazily when a downstream target has positive recurrence pole depth.  With one
SCC, or without the skeleton flag, PrepareChart constructs the full
ChartSystem immediately as before.

For a full ChartSystem it calls ``Indicial`ChartIndicial`` and assembles the
data of 3.2: each Jordan block is first multiplied by
the unique nonnegative epsilon monomial that clears the lowest valuation of
all entries in that block (one COMMON factor for every member, preserving
the unit Jordan superdiagonal); V = the full d x d matrix formed by
concatenating these epsilon-primitive chain vectors in family order (column
order = (family, root in descending a, chain position), eigenvector first —
matching the column order of 3.4/3.5); VInv = exact `Inverse[V]`
(Q(alpha)(eps) arithmetic); J = the corresponding block-diagonal Jordan
data; per family, CollisionDepth is the executed homogeneous Taylor-layer
budget defined below (parallel CASE-P targets at one `n` count only by the
largest target Jordan block, while distinct positive `n` layers compose); the
exact nonsingularity of V is certified here; the epsilon monomial inserted
into its determinant by primitive normalization is explicitly accounted by
`SpectralBlockEpsShifts` (E2 on failure).  n = 0 collisions arriving from Indicial
canonically ordered (one record per unordered pair, Indicial.md 2.6
step 5) are SYMMETRIZED here: the reversed pair (DeltaB negated) is added
when absent, so the Family record of 3.2 always carries all ordered pairs
with integer offset >= 0.  +~40 lines against the §9 budget (the former
OQ1 reserved +20 for the det-V certification; the remaining 20 funded by
taking part of §9 cut 1).

NOT exported (deliberately): no SolveGeneral convenience wrapper (Transport
assembles general solutions from columns + particular + matching weights);
no strategy/dispatch entry point; no per-eps-order entry point (the old
epsord-loop re-entry pattern of DispatchStrategy at
IntegrationStrategies/Dispatch.m:5 dies — symbolic eps means ONE solve per
chart).

---

## 3. DATA CONTRACTS

### 3.1 LocalSolution / Sector / EpsWindow (RewritePlan 3.1, verbatim)

Solve.m produces LocalSolution objects exactly as defined in RewritePlan 3.1:

```
LocalSolution = <|
  "Center" -> exact x0, "ChartMap" -> affine (Mobius optional, see 3.2),
  "Radius" -> distance to nearest singularity IN THE COMPLEX PLANE ...,
  "Sectors" -> { Sector.. },
  "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "TWindow"   -> <|"CompleteMax" -> nmax|>,
  "ErrorEstimate" -> per (eps-order) accumulated error ...,
  "Prescriptions" -> LIST of <|"Factor", "Sign", "Multiplicity",
              "LeadingCoeffSign"|> ...
|>

Sector = <|
  "a" -> rational, "b" -> exact rational/algebraic, "p" -> integer >= 0,
  "Coeffs" -> c[[k, n, comp]]    (* eps-order x t-power x COMPONENT *)
|>
Value(x)_comp = Σ_sectors t^a t^(b eps) (eps Log t)^p / p! ·
                Σ_k eps^k Σ_n c[k,n,comp] t^n
```

Solve.m fills: "Sectors", "EpsWindow", "TWindow"; copies "Center",
"ChartMap", "Radius", "Prescriptions" through from the ChartSystem
unchanged; initializes "ErrorEstimate" to the zero map (Transport
accumulates).  The (eps Log t)^p/p! normalization is load-bearing: under it,
true-resonance log-chain weights are eps-independent and the stored
coefficient of the p-th member carries an exact eps^(−p) monomial factor
(kmin = −p BY CONSTRUCTION — RewritePlan 3.1 invariant exemption).
"Coeffs" arrays are EpsSeries-backed: per (n, comp) a truncated eps-Laurent
array at WorkingPrecision with an explicit honest window (EpsSeries.m
contract).  Tags a, b, p are EXACT (rational/algebraic); coefficients are
the only numerics.

### 3.2 ChartSystem (assembled by PrepareChart, 2.5)

RESOLUTION (M0 review; REVIEW-minimalism defect 15 + REVIEW-math D6,
adapted per DEC-7): (i) Solve.m consumes Indicial's Family/EigRecord
shapes VERBATIM (Docs/specs/Indicial.md 3.3/3.4 are normative; collision
keys n/LowerIdx/UpperIdx/DeltaB/Type — this spec's former From/To/Offset
names are renamed to those).  (ii) Solve.m BUILDS the spectral frame
itself, once per chart inside PrepareChart: V = the full d x d matrix of
all families' epsilon-primitive chain vectors in spec column order, VInv =
exact Inverse[V]
(Q(alpha)(eps) arithmetic), certifies exact nonsingularity, and records the
known eps -> 0 determinant valuation inserted by the blockwise primitive
shifts (the E2 third check; failure stays E2).
(iii) CollisionDepth is computed from executed homogeneous CASE-P layers:
for each source root, discard `n=0` LaurentShift records, group its remaining
records by positive Taylor order `n`, sum the maximum target BlockSize in
each group, and take the maximum over source roots.  Thus parallel divisions
at one `n` are not added, while divisions at distinct Taylor layers compose.
(iv) ChartSystem is ASSEMBLED BY PrepareChart (DEC-7 — this supersedes
the review's Transport-assembly variant): chart geometry fields from
Transport's chart record (Transport.md 3.2) + IndicialData fields, with
"ThetaMatrix"/"Gauge"/"GaugeInverse"/"Residue" lifted from
IndicialData["Reduction"].  (v) `ThetaOriginal` and its exact
`IntegrationSequence` are recorded before that reduction.  The full form is:

```
ChartSystem = <|
  "ChartVar" -> t (symbol),
  "Center" -> exact x0, "ChartMap" -> chart map record,
  "Radius" -> exact/numeric complex-plane distance,
  "Prescriptions" -> (RewritePlan 3.1 list; passed through),
  "SystemSize" -> d_Integer,
  "ThetaOriginal" -> t beta A(x0 + beta t), exact in the original
                   physical master basis,
  "IntegrationSequence" -> the exact certificate of 2.4,
  "ThetaMatrix" -> B(t,eps): d x d, exact rational in (t, eps),
                   HOLOMORPHIC at t = 0   (* theta f = B f; B = t·A after
                   Indicial's rank reduction; Solve NEVER sees a pole order
                   > 1 in A — that is Indicial's contract, E1 enforces *),
  "Gauge" -> T(t,eps), "GaugeInverse" -> T^(-1)(t,eps): exact rational
                   matrices; Identity when no rank reduction was needed.
                   Original-frame solutions f = T·g where g solves
                   theta g = B g + T^(-1)·s   (old pattern:
                   ResonantRecurrence.m:1489 bVecG = TInv.bVec and
                   :1497 fGeneral = T.gGeneral),
  "Residue" -> M0(eps) = B(0,eps) exact,
  "SpectralBlockEpsShifts" -> one nonnegative integer per Jordan block,
         recording the common eps monomial applied to every chain member,
  "V" -> V(eps), "VInv" -> V^(-1)(eps): full d x d, exact in eps,
         M0·V = V·J, BUILT BY PrepareChart from epsilon-primitive versions
         of the families' Chains
         (column order = (family, root in descending a, chain position),
         eigenvector first),
  "J" -> block-diagonal Jordan data in the same column order: per root i
         the block λ_i·I_{q_i} + N_i (N_i = upper shift; chain vectors
         are UNIT vectors in the J-frame),
  "Families" -> { Family.. }
|>

Family = <|
  "Roots" -> {<|"a" -> a_i, "b" -> b_i, "BlockSize" -> q_i|> ..}
             (* integer-spaced a within the family: a_i - a_j ∈ Z;
                Σ q_i over ALL families = SystemSize *),
  "ColumnRange" -> the family's contiguous column indices in V/J,
  "Collisions" -> {<|"n" -> n_ij = a_j - a_i (>= 0), "LowerIdx" -> i,
                     "UpperIdx" -> j, "DeltaB" -> b_i - b_j,
                     "Type" -> (Indicial.md 3.4 vocabulary)|> ..}
                  (* all ordered pairs with integer offset >= 0 AFTER
                     PrepareChart's n = 0 symmetrization (2.5); entries
                     with DeltaB == 0 and n > 0 are TRUE resonances;
                     DeltaB != 0 (Type == "LaurentShift") are PSEUDO
                     resonances (any n >= 0) *),
  "CollisionDepth" -> Max over source roots i of
                      Sum over executed n >= 1 layers of
                      Max[BlockSize of each LaurentShift target at that n]
                      (integer >= 0; computed by PrepareChart; n=0 records
                      are diagnostic degeneracies, not homogeneous
                      recursion layers)
|>
```

When `UseSCCSkeleton -> True` discovers more than one SCC, the lazy form is
instead

```
SCCEnvelope = <|
  "ChartVar", "Center", "ChartMap", "Radius", "Prescriptions",
  "SystemHash", "SystemClearKey", "SystemSize",
  "ThetaOriginal" -> exact original theta-form matrix,
  "IntegrationSequence" -> exact certificate of 2.4,
  "IndicialData" -> <|"Regular" -> True|False,
                       "Dimension" -> d, "SCCSkeleton" -> True|>,
  "SCCSkeleton" -> True, "ChartSystemKind" -> "SCCEnvelope"
|>
```

The envelope deliberately has no global `ThetaMatrix`, gauge, residue,
spectral frame, or families.  It retains `SystemClearKey` so diagonal blocks
can be rebuilt from the registered original physical matrix, and
`SystemHash` supplies their shared parent solve-cache tag.  A regular value
transport may materialize a transient identity gauge/frame, but the cached
envelope remains skeletal.

Requirements certified at PrepareChart assembly (violations are E2):
char-poly factorization already certified by Indicial (I1 contract);
Det[V(eps)] is not identically zero and the monomial valuation deliberately
inserted by primitive normalization is `Sum[q_i shift_i]`;
Σ BlockSizes = SystemSize; ThetaMatrix entries have no
t-pole (exact check on the rational form).  These spectral requirements
apply to a full ChartSystem and independently to every prepared diagonal SCC;
the lazy envelope instead certifies the exact graph invariants in 2.4 and its
ordinary/singular classification.  SolveHomogeneous/SolveParticular re-assert
the full-chart invariants cheaply on receipt (Σ BlockSizes, t-pole).

### 3.3 SourceData (consumed)

```
SourceData = <|
  "Sectors" -> { Sector.. }       (* same Sector shape as 3.1 *),
  "EpsWindow" -> <|"Min", "CompleteMax"|>,
  "TWindow" -> <|"CompleteMax"|>
|>
```

A LocalSolution is accepted wherever a SourceData is (only these three keys
are read).  NORMALIZATION CONTRACT: the source is in THETA form — Solve
accepts the source in the original physical frame.  On a rank-reduced chart
`f = T.g`, SolveParticular forms `T^-1.s` exactly once before solving
`theta g = B_reduced.g + T^-1.s`, and gauges the result back with `T`.  A
caller must neither pre-apply `T^-1` nor add another factor of `t`.  A caller
holding a d/dt-form source b (f' = A f + b) must convert s = t·b first (an
exact a -> a+1 tag shift via SectorSeries.m).  This is stated to kill the
off-by-one the old code encoded as `s = bLeadPow + 1`
(Recurrence.m:641).  Source sectors may have any exact tags, including
p > 0 and eps-Laurent coefficient windows with kmin < 0 (IBP poles).

### 3.4 FundamentalSystem (produced)

```
FundamentalSystem = <|
  "Columns" -> { LocalSolution.. }   (* length = SystemSize *),
  "Specs"   -> { <|"a", "b", "p", "Family" -> fIdx, "Root" -> i,
                   "ChainPos" -> q (0-based)|> .. }  (* parallel to Columns;
                   p = structural max log power actually populated *),
  "Diagnostics" -> <|
     "LogCeilings" -> per column the precomputed cutoff P (3.6),
     "PseudoCollisionsHit" -> per column the list of (n, i->j) collisions
        actually exercised,
     "WindowExtensions" -> per column the by-construction kmin extension
        (true-resonance -p and/or collision depth),
     "ResidualProbe" -> <|"t" -> t0, "MaxResidual" -> r|>,
     "CouplingDepth" -> longest chain in the block dependency DAG (set by
        SolveChart, 2.4; 0 for a direct SolveHomogeneous call — Transport
        reads it from here, REVIEW-math D6(b))
  |>
|>
```

### 3.5 Normalization and canonical choices (pinned here, once)

- Homogeneous column for (family F, root i, chain position q): in the
  J-frame the initial data is the exact monomial ladder
  û[n=0, ℓ] = ê^{(i)}_{q−ℓ+1} · eps^{−ℓ} for ℓ = 0..q, zero for ℓ > q,
  where ê^{(i)}_m is the m-th unit chain vector of block i (ê_1 = the
  eigenvector).  Equivalently: the column equals the CLASSICAL Frobenius
  solution t^{λ_i}(v log^q t/q! + ...) rewritten in the (eps Log t)^ℓ/ℓ!
  basis; its top log member has the eps-independent weight ê_1 stored at
  eps-order k = −q.  In the original frame the leading t^0 coefficient is
  the corresponding exact column of V(eps) (Indicial normalizes V: first
  nonzero component of each chain's eigenvector = 1, exact).  This is the
  old convention "initial condition at the highest log power is the
  eigenvector" (ResonantRecurrence.m:496-501) made exact.
- ALL kernel freedom met later in a recursion (true-resonant steps n > 0,
  and the n = 0 chain descent) is fixed to ZERO along the colliding
  eigendirections.  This is the canonical-basis choice; any other valid
  choice differs by multiples of other columns and would break parity
  pinning.  (The old code instead resolved n = 0 freedom from next-step
  solvability with pseudo-inverses — ResonantRecurrence.m:684-733; in exact
  arithmetic the constraint system is solved exactly and the REMAINING
  freedom, which the old code left to minimum-norm accidents, is pinned to
  zero.)
- Particular solutions: kernel components zero at every resonant step (the
  source fixes everything else).  "The" particular is thereby unique.
- Solutions are constructed on the PRINCIPAL SHEET (t real > 0, Log t
  real).  Solve performs no analytic continuation; crossing operators are
  Transport.m's job.

### 3.6 The recursion (normative)

Setup, once per chart (the exact spectral part lives in PrepareChart's
ChartSystem; the eps-expansion below happens once per SolveChart call —
see 2.1, memo table cut at M0):

1. CLEAR DENOMINATORS (the performance backbone; old fast path:
   Recurrence.m:13-77 RationalizeAMatrixCore, :206-230, :427-453; lesson
   restated at LocalSeries.m:396-422 ClearZDenominators).  Compute
   D(t,eps) = lcm of the t-denominators of B's entries, normalized by its
   content so that d_0(eps) := D(0,eps) is an eps-rational function with
   d_0(0) ≠ 0 (E3 if impossible — e.g. a factor (t − c·eps) whose root
   degenerates onto the chart origin at eps = 0; old exact check:
   LocalSeries.m:403-412).  Set N(t,eps) = D·B (polynomial in t, exact in
   eps).  Write D = Σ_{j=0}^{dD} d_j(eps) t^j, N = Σ_{j=0}^{dN} N_j(eps) t^j.
   Note N_0 = d_0·M0.  NO benefit heuristic, NO timeout, NO "direct series
   mode" alternative (forbidden fallback F2): the input is exact rational
   by the LoadSystem contract, so clearing always succeeds.

2. SPECTRAL FRAME: N̂_j := VInv·N_j·V, once per chart in the GLOBAL frame
   of 3.2 (V/VInv built by PrepareChart; J block-diagonal in (family,
   root, chain) column order).  N̂_0 = d_0·J(eps).  Expand every d_j(eps),
   N̂_j(eps), and V/VInv entry ONCE into EpsSeries arrays at
   WorkingPrecision over a deterministic WORK WINDOW.  A single worst entry
   is not a sufficient recurrence-pole budget: a negative-valued multiplier
   at Taylor offset `j` can be reused every `j` steps.  Define, for
   `1 <= j <= J = Min[nmax, Max[dD,dN]]`,

   ```
   stepCost[j] = Max[0, -Min[finite eps valuations of d_j and N̂_j entries]]
   depth[0] = 0
   depth[n] = Max[0, Max_{1 <= j <= Min[n,J]}
                         (depth[n-j] + stepCost[j])]
   recurrencePoleDepth = depth[nmax].
   ```

   Missing/structurally zero coefficients contribute cost zero.  This is a
   longest-path dynamic program on the executed Taylor dependency DAG; it is
   tighter than `nmax` times the worst pole when a costly multiplier first
   occurs at `j > 1`, and unlike the former one-hit bound it covers repeated
   composition.

   For a homogeneous solve let `Pmax` be the largest log ceiling below,
   `Cmax = Max[Family["CollisionDepth"]]`, `RPD` the depth above,
   `SPD = Max[0,-Min valuation(V)]`, and `GPD` the pole depth of exactly the
   gauge `t`-Laurent coefficients used by `applyGauge` (orders
   `t^gv .. t^(gv+nmax)`, where `gv` is the minimum `t` valuation of T).
   The implemented frame is exactly

   ```
   base0    = Min[reqMin,0] - Pmax - Cmax - 2
   top0     = reqMax + Pmax + Cmax + 2
              - Min[0, reqMin - Pmax - 2]
   workBase = Min[base0, reqMin - Pmax - Cmax - RPD] - SPD
   workTop  = Max[top0, reqMax + Pmax + Cmax + RPD] + SPD + GPD.
   ```

   `SPD + GPD` pays the honest-window losses of the final, separately
   applied back-transforms û -> V.û -> T.V.û; it is private scratch
   headroom and does not extend the delivered request.  `SPD` also lowers
   the shared physical frame base, whereas T is multiplied in its own local
   epsilon frame.

   For a particular source sector, the analogous executed-layer budget is
   computed directly from that tag, including `n=0`: group all different-b
   target blocks by `n = aTarget-aSource >= 0`, take the maximum target block
   size per `n`, and sum the layers.  With this `PD`, that sector's log
   ceiling `P`, source window `[sMin,sMax]`, and
   `SFD = poleDepth(VInv)+poleDepth(V)`, SolveParticular uses

   ```
   sourceTop0 = sMax + P + PD + 2 - Min[0, sMin - P - 2]
   workBase   = Min[Min[sMin,0] - P - PD - 4,
                    sMin - P - PD - RPD] - SFD
   workTop    = Max[sourceTop0, sMax + P + PD + RPD] + SFD.
   ```

   (`Gauge == Identity` is currently required for a nonzero particular
   source.)

   This is the
   old eps work-buffer rule (LocalSeries.m:444-455, :529-532;
   Docs/RecursiveFuchsianSolver.md "Epsilon Poles") made DETERMINISTIC: the
   buffer is computed exactly from the family data up front — there is no
   "repeat with a larger buffer until stable" loop (forbidden F6-analogue).

3. LOG CEILINGS (exact, up front):
   - column (root i, chain pos q):  P = q + Σ_{j: b_j = b_i, a_j > a_i} q_j
     (the old MaxLogPowers formula, ResonantRecurrence.m:489-493).
   - particular from source sector (a_σ, b_σ, p_σ):
     P = p_σ + Σ_{j: b_j = b_σ, a_j − a_σ ∈ Z≥0} q_j
     (each resonant passage bumps by at most the block size; this REPLACES
     the old kMaxInitial guess + grow-and-retry loop,
     ResonantRecurrence.m:1191-1195, :1242-1254 — forbidden F6).
   Trailing log members whose final coefficient arrays all classify zero
   under laurentLeadTol (relative, Tolerances.md as amended — the
   Max[10^(-Floor[chopDigits/2]), 10^-24] campaign calibration, DEC-2;
   tested through Tolerances`NumericallyZeroQ per DEC-3, scale = the
   column's own coefficient scale) are dropped at assembly, the drop
   recorded in Diagnostics; this is the ONE numeric structural decision in
   Solve.m and it is guarded by the residual check.

RECURSION (per column / per source-sector tag).  Ansatz in the J-frame with
tag exponent λ = a + b·eps (a, b from the column's root or the source
sector's tag):

```
g = Σ_{ℓ=0}^{P} t^λ (eps Log t)^ℓ/ℓ! · Σ_{n>=0} t^n û[n,ℓ],
û[n,ℓ] ∈ (EpsSeries)^d
```

Using theta(t^{λ+n} w_ℓ) = (λ+n) t^{λ+n} w_ℓ + eps·t^{λ+n} w_{ℓ−1} with
w_ℓ = (eps Log t)^ℓ/ℓ!, the cleared equation D·theta g = N·g + D·ŝ
(ŝ = V^(-1)·s in the J-frame, per-sector EpsSeries arrays β̂[n,ℓ]) gives,
at each (n, ℓ):

```
L_n(eps)·û[n,ℓ]  =  R[n,ℓ] − d_0(eps)·eps·û[n,ℓ+1]

L_n(eps) = d_0(eps)·[ (λ+n)·I − J(eps) ]

R[n,ℓ] =   Σ_{j=1}^{min(n,dN)}  N̂_j(eps)·û[n−j,ℓ]
         − Σ_{j=1}^{min(n,dD)}  d_j(eps)·[ (λ+n−j)·û[n−j,ℓ]
                                           + eps·û[n−j,ℓ+1] ]
         + Σ_{j=0}^{min(n,dD)}  d_j(eps)·β̂[n−j,ℓ]      (particular only)
```

with û[·, P+1] := 0 and û[m, ·] := 0 for m < 0.  Per t-order the cost is
O((dD + dN)·d²) EpsSeries multiply-adds — finite width, NOT an O(n)
convolution (this is exactly the old denominator-cleared recursion,
Recurrence.m:254-262 and :494-511, in theta form with the log tower of
:356-395 folded in).

WHAT IS SOLVED AT EACH n: one finite linear block over the unknowns
{û[n,ℓ] : ℓ = 0..P}.  L_n is block-diagonal over the family's Jordan
blocks; for block i (size q_i) the diagonal operator is
d_0·(δ_i(n,eps)·I_{q_i} − N_i) with the INDICIAL OFFSET

```
δ_i(n, eps) = (a + n − a_i) + (b − b_i)·eps .
```

Block-by-block, top-down in ℓ:

- CASE T (Taylor): a + n − a_i ≠ 0.  Invert exactly:
  (δ·I − N_i)^(-1) = Σ_{m=0}^{q_i−1} N_i^m / δ^{m+1}.  All divisions are
  EpsSeries divisions by δ_i (unit at eps = 0): no window effect.  This is
  the generic step and the ONLY step a regular chart ever executes
  (δ = n, V = I, J = 0 — the classic-transport degenerate path).

- CASE P (pseudo-resonant): a + n = a_i, b ≠ b_i, so δ_i = (b − b_i)·eps
  EXACTLY.  The division is an honest EpsSeries LAURENT division: the
  quotient's window shifts down by one per power of δ_i (EpsSeries.m
  window-shift contract).  Record the shift; then apply the JOINT
  compensation rule of section 3.7.

- CASE R (true-resonant): a + n = a_i AND b = b_i, so δ_i ≡ 0 identically
  in eps.  No division happens.  The block-i rows of level ℓ read
  −d_0·N_i·û[n,ℓ]_i + d_0·eps·û[n,ℓ+1]_i = R[n,ℓ]_i: they DETERMINE the
  next-higher log member û[n,ℓ+1]_i (the LOG BUMP, p → p+1: along the
  eigendirection the equation is d_0·eps·û[n,ℓ+1]_eig = R_eig, an EXACT
  monomial eps-division — the normalization's eps^{-1}, accounted as
  by-construction kmin extension, not as data loss for homogeneous columns;
  for particulars it is a real honest window shift of one order per hit,
  the integral counterpart of Integrate.m's (a+n+1+b·eps)^{p+1}
  denominators).  Implementation: solve the block-i ℓ-LADDER at this n as
  ONE small exact linear system over {û[n,ℓ]_i : ℓ = 0..P} (the
  SolveSingularNBlock structure, ResonantRecurrence.m:891-935, but exact,
  square once the kernel directions are pinned to zero per 3.5 — no
  PseudoInverse, no projection).  If the ladder is INCONSISTENT at the
  precomputed ceiling P, that is a loud error E5 (the ceiling formula is
  then provably wrong — a bug, not a runtime condition to retry).

RESONANT-SOURCE LOG-BUMP RULE (the I2/finding-9 decision, normative): when a
source sector tag exactly hits a homogeneous tag — b_σ = b_i and
a_i − a_σ ∈ Z≥0 — the particular ansatz MUST include log members up to the
ceiling of step 3 (p_σ + Σ q_j), and the bump materializes through CASE R
above at n = a_i − a_σ.  The old general solver silently returned an EMPTY
particular on exactly these inputs (Docs/FeynmanTrickBoxFamilyStatus.md:
139-142; mechanism: ResonantRecurrence.m:1178-1183 returns a zero particular
when SourceBasePowers comes back {}, and PowerCoefficient silently maps
extraction failures to 0 at :1116-1139).  There is no old-code parity oracle
for this path; the M3 closed-form pin is test SU-08 (section 8).

PARTICULAR ASSEMBLY: one recursion per source-sector tag class (source
sectors sharing (b_σ, a_σ mod 1) align with one family's frame; tags
aligning with NO family run pure CASE T).  Sum of the per-tag particulars;
transform each source into the spectral frame with VInv, solve there,
back-transform û -> V·û, then apply the gauge: f = T·g via SectorSeries
rational-multiply.  The scratch window pays the exact epsilon-pole depths
of both VInv and V, in their actual order, before returning the honest
minimum window over contributing pieces.  CASE-P source hits use the same
registered homogeneous-target compensation and two-probe certificate as
homogeneous columns; their finite source window still records every real
collision-order loss described below.

### 3.7 Joint pseudo-resonant construction (the I2 spec decision, normative)

At a CASE P collision (column from root σ, colliding block i, offset n), the
quotient û[n,·]_i = R/((b−b_i)·eps)·(chain-corrected) acquires a finite
POLAR PART γ(eps) = Σ_{k<0} γ_k eps^k along the block-i chain directions
(a size-`q` target can consume `q` orders at that layer; distinct executed
Taylor layers can compose).  The joint rule:

0. The spectral seeds are the epsilon-primitive Jordan blocks constructed by
   PrepareChart (2.5).  Indicial's deterministic first-entry normalization
   may be meromorphic when affine eigenvalues coalesce at eps=0; leaving that
   arbitrary pole in a source seed makes a value-finiteness certificate
   indistinguishable from a failed collision.  Multiplying a whole block by
   eps^s changes only the basis normalization, moves the inverse factor into
   VInv, and preserves every Jordan identity.  Consequently any negative
   order inspected below was created by the recurrence/collision rather than
   inherited from an arbitrary eigenvector normalization.
1. Keep the full Laurent quotient in û (the recursion continues with it;
   EpsSeries windows track honestly).
2. REGISTER the compensation term  −Σ γ-polar(eps) × (family column of root
   i, matching chain member)  onto the object under construction: extra
   Sectors with root-i tags whose coefficients are (−γ-polar) × the
   already-built column-i coefficient arrays (the EpsSeries-weighted
   combination `SectorSeries`CombineLocalSolutions`, SectorSeries.md 2.10 —
   REVIEW-minimalism defect 6).  Columns are built
   in DESCENDING-a order within each family precisely so column i exists
   when the collision fires.  CASE-P divisions at one Taylor order act on
   disjoint target blocks in parallel; the work budget charges their maximum
   BlockSize once.  Only distinct increasing `n` layers compose, exactly as
   recorded by `CollisionDepth` (or the per-source `PD` for particulars).
3. The delivered object (column or particular) therefore has VALUE finite
   at eps = 0 through its honest complete window: the negative eps-orders
   introduced by each registered collision cancel.  This is ASSERTED at
   two deterministic positive interior probes in the gauge-reduced physical
   frame `g = GaugeInverse.f`: each registered spectral collision order is
   mapped through the exact epsilon valuations of its target column in `V`,
   and the corresponding reduced-component polar coefficients are checked.
   Applying `VInv` in the certificate is forbidden: when the spectral frame
   has an epsilon-singular inverse (banana L2 at x = 0 is the minimal
   example), a finite physical value has legitimate polar spectral
   coordinates.  Probe radii
   are chosen deterministically inside the disk so the unmatched top Taylor
   tail of a finite-order compensated target lies below `laurentLeadTol`.
   Every checked coefficient includes its own significance uncertainty.
   This value check, not the ODE residual (which any added homogeneous column
   also satisfies), is the runtime discharge of invariant I-5.  Failure is
   E5-class, loud.
4. Window accounting: per-sector arrays carry kmin lowered by the collision
   depth BY CONSTRUCTION (exemption class, like kmin = −p); for HOMOGENEOUS
   columns the internal work buffer (3.6 step 2) restores the delivered
   CompleteMax to the request — pseudoResonanceShift = 0 as RewritePlan 3.4
   asserts.  For PARTICULARS the source window is finite and the shift is
   REAL: CompleteMax drops by one per collision on the solved path,
   recorded in Diagnostics and in the delivered EpsWindow.  Never padded,
   never hidden.

Banana reality check (why this is day-one scope): the banana L2/L1 endpoint
data mixes x^0, x^(−1+eps), x^(2eps) — one family, offsets {0, 1}, three
distinct b values (Docs/FeynmanTrickBananaStatus.md:340-348; RewritePlan I2).

---

## 4. INVARIANTS (always on, cheap)

- I-1 TAGS EXACT: every Sector has exact rational/algebraic a, b and
  integer p ≥ 0; coefficients are the only numerics.  Checked at assembly.
- I-2 COMPLETENESS OF THE BASIS: number of delivered columns =
  SystemSize = Σ BlockSizes; specs enumerate every (root, chain position)
  exactly once.
- I-3 REGULAR CHART DEGENERACY: a regular chart (M0 = 0; per DEC-6 one
  (0,0,0) FAMILY whose single root has Multiplicity = d, BlockSizes =
  ConstantArray[1, d], V = I — d columns, d sector specs) produces exactly
  one sector (a=0, b=0, p=0) per column and the SAME code path executes
  (no SolveSimple-style special
  case — the old special case was the NormalizeLogPower hole's host,
  Default.m:6-18 + Docs/FeynmanTrickBoxFamilyStatus.md:57-73).
- I-4 WINDOW EXEMPTIONS ARE STRUCTURAL: kmin extensions appear ONLY as
  (a) −p on true-resonance members, (b) −collisionDepth on compensated
  objects; both recorded in Diagnostics["WindowExtensions"], and window
  arithmetic accounts for them rather than erroring (RewritePlan 3.1
  EXEMPTION clause).
- I-5 EPS-REGULARITY OF COMPENSATED OBJECTS: at the probe point, the
  assembled value's eps-orders below the object's honest content (the
  negative orders introduced by collisions/log-chains that must cancel
  across sectors) vanish to matchTol relative.  (Homogeneous joint columns
  and compensated particulars only; single-sector objects skip.)
- I-6 ODE RESIDUAL: every delivered column and particular passes
  ODEResidualCheck at one random interior probe per chart (RewritePlan 3.1
  invariant).  Never sampled at t = 0 (logs vanish there spuriously — the
  old two-point-probe lesson, RewritePlan §5).
- I-7 WINDOW HONESTY: EpsWindow/TWindow of every output are the computed
  honest windows; no zero-padding anywhere; every Laurent-division shift on
  finite-window data is visible in the delivered window.
- I-8 STRUCTURAL DECISIONS ON EXACT DATA ONLY: collision detection,
  resonance classification (IntegerQ on exact rational differences, exact
  b-equality), log ceilings, and block structure are computed from exact
  tags/spectral data BEFORE any numeric coefficient exists.  The
  IntegerQ-on-floats / snapped-spectrum disease class
  (Recurrence.m:649-663; ResonantRecurrence.m:379-394;
  Docs/FeynmanTrickBoxFamilyStatus.md items 3-4) is impossible by
  construction, and any code path that would re-derive structure from
  numeric coefficients is a review-rejectable defect.

---

## 5. ERROR CONTRACT

All errors are raised via ``Tolerances`DE2Error[id, payload]`` (DEC-1: the
library-wide primitive — prints a one-line summary and
`Throw[Failure["DiffExp2", payload], "DiffExp2Error"]`, caught only at
API.m entry points); this module defines no other error mechanism
(REVIEW-math D17 / REVIEW-minimalism defect 5, as superseded by DEC-1).
`id` is the E-tag below; the payload always carries "ID", "Module" ->
"Solve", and AT MINIMUM: the chart identification (Center, ChartVar,
ChartMap summary), the family roots involved, the sector tags (a, b, p),
the t-order n and log level ℓ where applicable, and the
requested-vs-available windows where applicable.  Enumerated conditions:

- E1 NON-FUCHSIAN INPUT: ThetaMatrix has a t-pole (exact rational check at
  receipt).  Message: chart + the offending entry's pole order.  (Indicial
  owns rank reduction; Solve refuses rather than reducing again.  Old
  analogue: ResonantRecurrence.m:1465-1472.)
- E2 SPECTRAL DATA MISMATCH: Σ BlockSizes ≠ SystemSize; or M0·V − V·J ≠ 0
  (exact); or V is identically singular; or a recorded primitive block
  shift does not clear all negative valuations with minimum zero.  Message:
  chart + family roots + which check failed.
- E3 DEGENERATE CLEARED DENOMINATOR: d_0(eps)|_{eps=0} = 0 after content
  normalization (a singularity location degenerating onto the chart origin
  at eps = 0).  Message: chart + the offending denominator factor.  (Old
  exact analogue: LocalSeries.m:403-412.)
- E4 WORK-WINDOW OVERFLOW: an EpsSeries operation requests content outside
  the precomputed work window (would indicate the deterministic buffer
  formula of 3.6 is wrong).  Message: chart + operation + window arithmetic.
- E5 RESONANT LADDER INCONSISTENT / EPS-REGULARITY FAILURE: the CASE R block
  system has no solution at the precomputed ceiling P, or invariant I-5
  fails.  Message: chart + family + n + ℓ + tags + (for I-5) the residual
  eps-orders and magnitudes.  NO retry, NO ceiling growth.
- E6 EMPTY DELIVERED WINDOW: a particular's honest EpsWindow comes out with
  CompleteMax < Min after the recorded shifts (the source cannot support
  even one output order).  Message: chart + source tags + source window +
  the shift ledger.  (The old "keep at least one output order" floor
  masked exactly this — RewritePlan A1; deleted.)
- E7 ODE RESIDUAL FAILURE: ODEResidualCheck above residualTol.  Message:
  chart + column spec / particular tags + probe t0 + per-eps-order residual
  table.
- E8 NON-FINITE COEFFICIENT: any Indeterminate/ComplexInfinity/overflow in
  a coefficient array, at the moment of creation.  Message: chart + tag +
  (k, n, comp).  (Replaces TrimNonFiniteSeriesTails' trim-and-warn,
  ResonantRecurrence.m:62-120 — forbidden F7.)
- E9 MALFORMED REQUEST/SOURCE: req windows inverted (Min > CompleteMax);
  source sector tags non-exact (inexact a/b heads, negative or non-integer
  p); source TWindow/EpsWindow inconsistent with its Coeffs dimensions.
  The theta-form normalization contract (3.3) is NOT mechanically
  detectable and is not guessed at: a d/dt-form source produces wrong
  (shifted-tag) results that the always-on ODE residual check (E7, which
  includes the source term) catches — that is the designed backstop, and a
  code comment at the SourceData validator must say so.  (REVIEW-math D27
  / REVIEW-minimalism defect 32: Solve never integrates in t; b = 0
  endpoint divergences are wholly Integrate.m's table.)

NO SILENT FALLBACKS.  Every place a fallback might be tempting, explicitly
forbidden (each cites the old site that must NOT be reproduced):

- F1 NO strategy dispatch ladder, no "if this solver declines, try that
  one" (Dispatch.m:87-131; Recurrence.m:855-861 resonance handoff;
  Default.m:45-47 pivot fallback to VOPAlt; Wronskian.m:85-108 pivot
  retry).  One solver, one path; any "not applicable" condition is a typed
  loud error (E1-E3).
- F2 NO rationalization escape hatches: no TimeConstrained with a {False}
  fallback, no "degree too high, use series mode" benefit heuristic, no
  direct-series O(n²) alternate mode (Recurrence.m:25-74 incl. the
  5-second timeout at :73, :66-68 degree heuristic, :294-307 mode switch).
- F3 NO "divisor == 0 -> coefficient := 0" (Recurrence.m:505-507, :521-527
  — commented "Should not happen ... but be safe").  Exact offsets make
  the zero case a STRUCTURAL branch (CASE R), never a numeric surprise.
- F4 NO numeric rescue chain LinearSolve -> LeastSquares -> PseudoInverse,
  no projection of an inconsistent rhs onto the column space, no
  minimum-norm kernel choices (SafeNumericLinearSolve
  ResonantRecurrence.m:529-587; SolveRecurrenceStep :599-640;
  SolveSingularNBlock's PseudoInverse :926; LocalSeries.m:572-580).
  Inconsistency = E5.
- F5 NO silent empty particular: a particular may be zero ONLY for a
  structurally empty source (3.3); never because power extraction failed
  (ResonantRecurrence.m:1178-1183 empty SourceBasePowers -> zero
  particular; PowerCoefficient -> 0 on $Failed at :1116-1139, the
  box-campaign eps^0 deficit's root mechanism per
  Docs/FeynmanTrickBoxFamilyStatus.md:90-114).  In the sector-native
  representation "leading power" is a tag, not an inference — there is
  nothing to extract.
- F6 NO grow-and-retry loops for log ceilings or eps buffers
  (ResonantRecurrence.m:1242-1254 kMax growth; the
  Docs/RecursiveFuchsianSolver.md "repeat with a larger epsilon buffer"
  loop): both ceilings are computed exactly up front (3.6); a miss is a
  bug surfaced by E4/E5.
- F7 NO trimming of non-finite coefficient tails with a warning
  (TrimNonFiniteSeriesTails, ResonantRecurrence.m:62-120).  E8 instead.
- F8 NO type-normalization warn-and-continue: the old code warned when a
  particular failed to normalize to SeriesData and proceeded
  (Recurrence.m:574-579, :777-787).  In DiffExp2 coefficients are EpsSeries
  arrays end-to-end; a type mismatch is a hard error at the boundary it
  occurs.
- F9 NO "accept unless residual PROVEN above tolerance" significance-zero
  logic (ResonantRecurrence.m:1232-1239): residual comparisons use fixed
  WorkingPrecision values and relative tolerances from Tolerances.m; an
  undetermined comparison is itself E8-class.
- F10 NO zero-filling, no phantom orders, no "keep at least one order"
  floors anywhere in window assembly (RewritePlan A1 catalogue).
- F11 NO fallback to any Wronskian/Frobenius-style scalar-reduction path —
  those modules are deleted, and the old refusal messages
  (Dispatch.m:110-118, ResonantRecurrence.m:1359-1368) become the
  unconditional behavior.
- F12 NO numeric leading-power detection / noise-skipping inside Solve
  (Recurrence.m:587-639 with the load-bearing 1e-24 floor at :618).  The
  representation makes it meaningless here; the LESSON (relative, not
  absolute, judgments on cancellation residue) lives on in Tolerances.m
  and in the matching solves (RewritePlan §5 seed).

---

## 6. ABSORBED OLD CODE

Solve.m replaces the following (RewritePlan I2: the ~3.45k-line strategy
stack + Wronskian + Frobenius, plus LocalSeries' finite-width core):

| Old code | Lines | Disposition |
|---|---|---|
| DiffExp/IntegrationStrategies/Recurrence.m | 870 | replaced entirely |
| DiffExp/IntegrationStrategies/ResonantRecurrence.m | 1527 | replaced entirely |
| DiffExp/IntegrationStrategies/VOP.m | 304 | deleted (scalar-reduction VOP concept dropped) |
| DiffExp/IntegrationStrategies/Default.m | 197 | deleted (SolveSimple/SolveDefault) |
| DiffExp/IntegrationStrategies/Dispatch.m | 132 | deleted (no dispatch) |
| DiffExp/IntegrationStrategies/Helpers.m | 124 | deleted |
| DiffExp/IntegrationStrategies.m | 43 | deleted (loader) |
| DiffExp/Frobenius.m | 141 | deleted (scalar Frobenius + order reduction) |
| DiffExp/Wronskian.m | 112 | deleted (n-th order combination + Wronskian inversion) |
| DiffExp/LocalSeries.m (finite-width core + recursion doc) | part of 1026 | subsumed: ClearZDenominators :396-422, PrepareFiniteWidthData :431-472, RecursiveFiniteWidthSolve :500-720, SolveLocalFuchsianSeries :728-770; FuchsianizeLocal/TrimFuchsianLattice (:188+, :15-16) are ported by INDICIAL, not here |
| Docs/RecursiveFuchsianSolver.md | — | the recursion's design note; 3.6 is its production form |

Specific functionality mapping with file:line:

- Denominator-cleared finite-width recursion (homogeneous):
  Recurrence.m:241-273 (core at :254-262) and singular variant :475-545
  (rational-mode recursion :491-511) -> section 3.6 recursion.  The
  eigenbasis transform-once pattern (:449-452 rHat = PInv.R.P, back at
  :533-540) -> the J-frame.
- Log-tower particular at regular points with top-down ℓ sweep:
  Recurrence.m:324-407 (ordering rationale in the comment :345-348) ->
  the ℓ = P..0 sweep.
- Singular particular with source leading-power handling:
  Recurrence.m:549-790 -> SolveParticular; its five campaign bug fixes
  (Docs/FeynmanTrickBoxFamilyStatus.md:90-124) are obsoleted by
  construction (tags replace power inference) — see lessons below.
- Resonance structure / Jordan chains / MaxLogPowers:
  ResonantRecurrence.m:362-527 -> Indicial supplies exact J/V; the
  MaxLogPowers formula :489-493 -> log ceilings (3.6 step 3); initial
  vectors :496-501 -> normalization (3.5).
- Resonant fundamental matrix with kernel-freedom resolution:
  ResonantRecurrence.m:646-854 (n = 0 descent :684-733, n >= 1 :738-833)
  -> CASE R ladder + canonical zero kernel (3.5).
- Coupled (n, all-ℓ) block solve at singular steps:
  SolveSingularNBlock ResonantRecurrence.m:891-935 -> CASE R ladder
  system (exact, square).
- Unified particular over source power classes:
  ResonantRecurrence.m:1141-1276 -> per-source-tag recursion (3.6).
- Fuchsianized-frame wrapping: ResonantRecurrence.m:1459-1527 (TInv·source
  at :1489, T·solution at :1497) -> Gauge/GaugeInverse contract (3.2/3.3).
- Residual self-check: CheckParticularResidual ResonantRecurrence.m:987-1023
  and the gated DebugCheckBlockSolution :1391-1457 -> ODEResidualCheck
  (always-on, not env-gated).
- Config key `UseRationalRecurrence` (State.m:26, :127; Dispatch.m:19-57;
  Transport.m:141): MEANING ABSORBED — the cleared-denominator recursion is
  unconditional; the key is dropped (Config.md owns the waiver entry).

NUMERICAL LESSONS THAT MUST SURVIVE (read from the old code, normative):

- L1 Finite width beats convolution: per-order cost must be O(dD + dN), not
  O(n) (Recurrence.m:254-262, :494-511; R6 benchmark mechanism per the
  execution review finding 7).
- L2 Transform to the spectral frame once, back-substitute per order;
  never re-decompose per step (Recurrence.m:449-452, :533-540).
- L3 The ℓ-sweep is top-down because theta of w_{ℓ+1} feeds level ℓ at the
  same t-power (Recurrence.m:345-348, :368-378; ComputeRecurrenceRHS
  ResonantRecurrence.m:865-882).
- L4 n = 0 Jordan-chain freedom must be consistent with the next chain
  equation — the old minimum-norm choice silently broke it (banana bug 4,
  Docs/FeynmanTrickBananaStatus.md:248-253;
  ResonantRecurrence.m:684-733).  In DiffExp2 the exact ladder solve plus
  the pinned-zero canonical freedom makes this structural; the unit test
  SU-05 locks it.
- L5 Structural decisions need exact arithmetic: the old code learned to
  rationalize before lattice reduction (ResonantRecurrence.m:144-154,
  banana bug 1) and to distrust snapped numeric spectra
  (:227-351, :379-394, banana bugs 2-3).  DiffExp2 generalizes: structure
  from exact data only (invariant I-8).
- L6 Precision floors: coefficient arrays live at fixed WorkingPrecision;
  significance-collapsed inputs are re-fixed to WP, never allowed to veto
  a solve or poison comparisons (SafeNumericLinearSolve's lesson at
  ResonantRecurrence.m:536-545, kept; its rescue ladder, dropped per F4).
- L7 Cancellation residue sits far above absolute chop scales; all
  tolerance judgments relative (Recurrence.m:594-618 incl. the 1e-24 floor
  rationale at :612-618) — folded into Tolerances.m's
  laurentLeadTol/matchTol semantics (relative, with the load-bearing
  10^-24 floor; REVIEW-minimalism defect 1, DEC-2); inside Solve only the
  assembly-time log-member trim (3.6 step 3) and the invariant probes
  consume tolerances.
- L8 Zero-source short-circuit must be structural, not numeric-chop
  (old PChop test Recurrence.m:318, :846 -> empty Sectors list, 2.2).
- L9 Type poisoning is silent death: every old warn-and-continue
  normalization site (Recurrence.m:556-585, :770-789) becomes a typed
  error (F8); coefficient containers are one type end-to-end.
- L10 The resonant-source hole: the empty-particular bypass
  (Docs/FeynmanTrickBoxFamilyStatus.md:139-142) is closed by the log-bump
  rule with closed-form pin SU-08 — the M3 gate item RewritePlan names.

---

## 7. DEPENDENCIES

Acyclic order (RewritePlan):
Tolerances < Config < EpsSeries < SectorSeries < Indicial < Solve <
Transport/Integrate < API.

Solve.m MAY call:
- Tolerances.m: residualTol, matchTol, rankTol, laurentLeadTol
  (assembly-time log-member trim, 3.6 step 3 — the binding consumer table,
  Tolerances.md §3 as amended by REVIEW-minimalism defect 17), the
  NumericallyZeroQ predicate (DEC-3), and DE2Error (the library-wide
  loud-error primitive, DEC-1).  Named semantics only; no literal
  thresholds in Solve.m.
- Config.m: WorkingPrecision, Verbosity (via the validated accessor);
  Solve takes windows/orders as EXPLICIT arguments — it never reads
  ExpansionOrder/EpsilonOrder-style globals (the old
  ctx["ExpansionOrder"]/FEC pattern dies with the dispatch layer).
- EpsSeries.m: all coefficient arithmetic — add, mul, scalar mul, Taylor
  and LAURENT division with window-shift semantics, polar-part extraction
  via the blessed ESCoefficientList idiom (EpsSeries.md §2,
  REVIEW-minimalism defect 36 — not a separate export), window queries.
  (If matrix-valued EpsSeries helpers are needed, they belong in
  EpsSeries.m, not here.)
- SectorSeries.m: LocalSolution assembly, rational multiply (gauge T,
  theta-form conversions, block-coupling sources in 2.4),
  CombineLocalSolutions (SectorSeries.md 2.10: EpsSeries-weighted
  combinations — compensation terms), evaluation at a point (residual
  probe), differentiate.
- Indicial.m: ChartIndicial, called ONLY from PrepareChart (2.5, DEC-7);
  Solve never re-derives spectral data — it assembles V/VInv/J from the
  delivered Chains (3.2 resolution).

Solve.m MUST NOT call Transport.m, Integrate.m, API.m, or any FeynmanTrick/
symbol, and must not read or write any global mutable state (no
DiffExp`State` analogue).

---

## 8. UNIT TESTS (Tests/test_solve_*.m; names normative)

Closed-form expectations are exact unless noted.  WP = 50 digits unless a
test says otherwise; comparisons relative 1e-40 where exact values exist.
Error tests assert via `Catch[..., "DiffExp2Error"]` returning the
`Failure["DiffExp2", ...]` with the named ID/Module payload fields (DEC-1).

> COVERAGE RULE (M0 review): every error ID and warning ID of section 5
> has at least one unit test that triggers it and asserts its required
> payload fields, OR a one-line waiver in the test file naming the ID and
> why it is untestable in isolation (e.g. "unreachable by construction,
> guarded by assert X").  The IDs currently missing tests are enumerated
> in Docs/specs/REVIEW-minimalism.md defect 25 (for Solve: E2, E4, E5, E6,
> E8); the implementation agent closes the list before the module's
> milestone gate.

- SU-01 `regular_exponential`: chart B = t·[[1]] (i.e. f' = f), regular.
  Assert: one column, single sector (0,0,0), c[0,n,1] = 1/n! exactly
  (rationals to WP) for n ≤ 12; same code path marker as singular charts
  (no special-case branch taken — assert via Diagnostics).
- SU-02 `regular_log_source`: regular scalar chart B = 0, theta-form source
  s = t·(eps Log t) (tag a=1, b=0, p=1).  Closed form:
  f = t·eps·(Log t − 1) = t[(eps Log t)/1! − eps].  Assert both sector
  members exactly; locks the top-down ℓ sweep (L3).
- SU-03 `frobenius_2f1`: the 2F1 companion system at x = 0 (chart matrices
  from Tests/Hypergeometric2F1_Matrices, exact load); residue eigenvalues
  {0, −c} (companion-system tags, DEC-17; for the committed dz_0.m this is
  {0, −3/2} — the scalar exponent 1−c appears as component-1 content at
  column n = 1 of the a = −c sector).  Assert the n ≤ 8 coefficients
  against the exact hypergeometric ratios (Pochhammer products) for both
  columns (M3 gate item).
- SU-04 `singular_nonresonant_2x2`: B with M0 = diag(−1+eps, 2eps) and
  N(t) = M0 + t·[[0,0],[1,0]].  Column 1 closed form (worked in 3.7's
  class): component 2 coefficient at n = 1 is −1/eps exactly (window
  [−1, ...]); joint compensation adds the (0, 2, 0)-tag sector with
  coefficient +1/eps × column 2.  Assert: tags exactly
  {(−1,1,0), (0,2,0)}; per-sector kmin = −1; I-5 probe cancellation; eps^0
  value of component 2 at t = 1/3 equals Log[1/3] to 40 digits (the
  (t^{2eps} − t^{eps})/eps -> Log t limit).  THE pseudo-resonance pin.
- SU-05 `true_resonance_jordan_chain`: B = M0 = [[0,1],[0,0]] (nilpotent,
  one root (0,0), q = 2).  Assert column 1 = (1,0) (single p=0 sector);
  column 2 sectors: p=1 member with coefficient (1,0) stored at k = −1
  (kmin = −p), p=0 member (0,1); n = 0 descent exact (L4); residual probe
  exact zero.
- SU-06 `integer_spaced_resonance_log_bump_homogeneous`: 2x2 with
  M0 = diag(0,2), N(t) = M0 + (t + t²)·[[0,0],[1,0]].  Worked closed form
  (derived in this spec's design record): column for root 0 =
  (1, −t + t²·Log t + O(t³)), i.e. p=1 sector appears at n = 2 with
  û[2,1] = (0, 1/eps).  Assert coefficients exactly through n = 3 and the
  structural p = 1.
- SU-07 `pseudo_resonance_three_sector_banana_class`: 3x3 diagonal-plus-
  coupling system with roots {(−1,1), (0,0), (0,2)} (the banana
  segment-1 endpoint mixture x^(−1+eps), x^0, x^(2eps) —
  Docs/FeynmanTrickBananaStatus.md:340-348).  Assert: one family,
  collision table has the (−1,1)->(0,0) and (−1,1)->(0,2) offset-1 pairs
  and the (0,0)<->(0,2) offset-0 pair; all three columns delivered with
  requested CompleteMax (homogeneous shift = 0, RewritePlan 3.4 assert);
  I-5 holds.  The exact seven-root banana-L1 companion pin has 24 collision
  records but only one executed positive-`n` layer: assert
  `CollisionDepth == 1` (ten parallel CASE-P visits must not be summed).
- SU-08 `resonant_source_log_bump_particular` (M3 closed-form pin, L10):
  scalar B = 0 with source tag (0,0,k): particular =
  (1/eps)·(eps Log t)^{k+1}/(k+1)! exactly, for k = 0,1,2.  And the offset
  variant B = 2·(scalar), source t²: particular = t²·Log t =
  (1/eps)·t²·(eps Log t).  Assert exact coefficients, p = k+1 structural,
  and the recorded one-order window shift of the particular.
- SU-09 `pseudo_resonant_source`: scalar B = 2eps, source tag (0,0,0)
  constant 1 with source window [0,K].  Particular: c = −1/(2eps), window
  [−1, K−1] — assert the shifted honest window (NOT K), the compensation
  +1/(2eps)·t^{2eps} sector, and eps^0 value Log t at the probe (exact
  (t^{2eps}−1)/(2eps) class).
- SU-10 `empty_source_zero_particular`: empty Sectors source -> zero
  particular, FULL requested window, no error (L8).
- SU-11 `window_honesty_no_padding`: source with CompleteMax = K and a
  CASE P collision: delivered CompleteMax = K−1; then a consumer-style
  read at order K must fail loudly naming (chart, sector, order) — i.e.
  the object carries no phantom order K (RewritePlan 3.1 consumer rule;
  the assertion here is on the produced metadata).
- SU-12 `residual_hook_catches_corruption`: corrupt one coefficient of a
  delivered column by 1e-10 relative and call ODEResidualCheck: assert E7
  fires and the message names the chart, the sector tag, and the probe.
- SU-13 `e3_degenerate_denominator`: B with denominator (t − eps): assert
  E3 names the factor.
- SU-14 `e1_nonfuchsian_refused`: ThetaMatrix with a 1/t entry: assert E1.
- SU-15 `box_l2_apparent_chart` (M3 gate; the box L2 J1/J2 pointwise
  20-digit pins were flagged NOT FOUND at M0 and are regenerated at M3 by
  the orchestrator's kernel queue per DEC-25 — Tests/PINS.md): the box L2
  apparent-singularity chart at
  t* = 7/11 (campaign: memory brief + Docs/FeynmanTrickBoxFamilyStatus.md:
  75-128).  Assert the J1/J2 pointwise pins to 20 digits and that the
  apparent chart's particular reproduces analyticity (the assembled
  general solution's a<0 content cancels — this exercises CASE T/P/R
  together on real campaign data).
- SU-16 `banana_l1_endpoint_chains` (M3 gate; fixture per Tests/PINS.md):
  the banana L1 endpoint system AFTER Indicial's rank reduction (nilpotent
  double pole upstream; 5-fold resonant residue, genuine size-3 Jordan
  block — Docs/FeynmanTrickBananaStatus.md:220-264).  Assert: chain
  structure (log² members present), gauge-wrapped output residual = 0 at
  probe, endpoint tower values against the vendored closed forms.
- SU-17 `dispatch_dump_replay_parity` (M3 gate): replay a vendored
  bubble-segment dispatch dump (old DEBUG_DUMP_DISPATCH_DIR capture,
  Dispatch.m:10-16; oracle values pre-generated in M0) through
  SolveHomogeneous/SolveParticular; agreement 1e-30 per RewritePlan M3.
- SU-18 `finite_width_scaling_guard`: structural performance assert — for
  a degree-3 rational 4x4 chart at TOrder 200, count EpsSeries
  multiply-adds per t-order (instrumented private counter): must be
  bounded by a constant in n (L1/R6 mechanism), not grow linearly.
- SU-25 `epsilon_primitive_casep_frame`: the exact epsilon-regular 3x3
  residue has raw source eigenvector `(1,1/eps,0)`, a coalescing partner,
  and an offset-one pseudo target `(0,1,1)`.  Assert PrepareChart records
  block shifts `{0,0,1}`, the source seed becomes `(eps,1,0)`, the CASE-P
  quotient is regular, and the uncompromised two-probe certificate passes.

---

## 9. LINE BUDGET

RewritePlan 3.2: Solve.m ~700 lines (the I2 figure ~700-800 INCLUDES the
ported rank reduction, which lives in Indicial.m; Solve proper is 700).
Indicative breakdown: clearing + spectral-frame prep + work-window math
~110; the recursion core (shared by homogeneous and particular — ONE
routine with β̂ = 0 for homogeneous) ~190; CASE R ladder + joint CASE P
compensation ~130; assembly/normalization/window bookkeeping ~110;
ODEResidualCheck + invariant probes ~70; errors/validation ~60; headers
~30.

M0 RESTATEMENT: +40 PrepareChart (DEC-7/REVIEW-math D6: 20 = the former
OQ1 det-V-certification reserve realized, 20 funded by taking cut 1 in
part) and +35 SolveChart (REVIEW-minimalism defect 3/DEC-9: funded by
taking cut 3, the memoization plumbing, ~−30).  Restated total: ~725; the
net +25 over the plan figure is DEC-7/DEC-9-imposed scope, accounted
against the defect-4 library headroom (restated totals ~3455 of 3500).
Landing >10% over 725 stops and reports to the orchestrator BEFORE
writing more code.

If over budget, cut in this order:
1. (TAKEN IN PART at M0: −20 funds PrepareChart) Move any matrix-EpsSeries
   convenience helpers into EpsSeries.m (they are generic; ~30-50 lines).
2. Trim Diagnostics to the fields the 3.4 budget validation and M5 ladder
   actually consume (LogCeilings, WindowExtensions, ResidualProbe,
   CouplingDepth).
3. (TAKEN at M0: funds SolveChart) The memoization plumbing is dropped —
   correctness-neutral; PrepareChart's explicit ChartSystem handoff (2.5)
   plus SolveChart's one-call solve replace its role on the Transport
   path.
4. Fold SU-18's instrumentation behind a single env-gated counter.

> PRE-AUTHORIZED at M0 review (REVIEW-minimalism.md defect 4), to be taken
> at implementation time WITHOUT further sign-off, freeing ~200 lines
> library-wide to fund the missing operations:
> - SectorSeries.md §9 cuts 1, 3, 5 are taken up front (Main-coordinates
>   option, multi-point evaluation form, public PadeEvaluate): -55; the
>   SectorSeries target is restated as 400 + 75 Pade = 475 TOTAL.
> - Indicial.md §9 cuts C-2 and C-3 are taken up front: -18 (C-1 trim-pass
>   stays IN until M2 evidence, because T-12/N-1b regression-pins it).
> - Transport.md §9 cut 3 (SavedCharts via caller callback): -30.
> - Tolerances.md: the dead adaptive-search constants are DELETED, not
>   folded (defect 8): -10.
> - EpsSeries.md §9 cuts 1-2 (binary ESAdd, drop ESZero): -13.
> - Integrate.md §9 cut 2 (b != 0 interior phase-paired path demoted to a
>   loud not-implemented error + ledger waiver) is pre-authorized but NOT
>   taken by default; it is the designated reserve (-50) if the post-fund
>   total exceeds 3450.
> Any module landing >10% over its restated figure stops and reports to
> the orchestrator BEFORE writing more code.

NEVER cut (load-bearing for correctness or the M5 gates): the joint
pseudo-resonant construction (3.7), the CASE R exact ladder + log bump, the
always-on residual check, the loud-error contract, the cleared-denominator
finite width (R6 benchmark fails without it).

---

## 10. OPEN QUESTIONS

- OQ1 — DELETED at M0 review (subsumed by the §3.2 RESOLUTION and
  PrepareChart, 2.5: Solve consumes Indicial's shapes verbatim, builds
  epsilon-primitive V/VInv itself, and certifies exact nonsingularity here;
  REVIEW-minimalism defect 15, REVIEW-math D6, DEC-7).
- OQ2 Algebraic (non-rational) a, b (R2): collision tests remain exact, but
  EpsSeries arrays with algebraic-number tag arithmetic in δ-divisions may
  be slow; if an example hits this, the division denominators
  (n + Δa + Δb·eps) should be numericized at WP AFTER the exact zero/
  collision classification.  Decide at first occurrence; flagged
  performance-only per R2.
- OQ3 The compensation construction's eps-regularity for collision chains
  of depth ≥ 2 is discharged at runtime by invariant I-5 rather than by a
  written proof.  If a depth-2 case fails I-5 in M3, the construction
  order (descending-a, polar-part-only γ) needs the full-Laurent-γ variant
  — both are valid bases; the spec pins polar-part-only as canonical and
  M3 revisits only on evidence.
- OQ4 — RESOLVED at M0 review (REVIEW-math D8; DEC-7): RAW columns are
  not needed, but recombination IS: joint compensation only fires when a
  Laurent division occurs in the recursion; eps=0 eigenvector collisions
  without any division (the log x class) leave det F(0) = 0 and are
  Transport's RecombineBasis job, keyed on Indicial's exported
  EpsDegenerateFamilies (DEC-7 canonical name).
- OQ5 The probe-point policy for charts whose Radius is not yet final at
  solve time (Transport may shrink radii after product operations per the
  SectorSeries multiply contract): current spec draws the probe from the
  chart's stated Radius; if Transport re-checks after shrinking, it calls
  ODEResidualCheck with an explicit probe argument.
- OQ6 Budget interplay: the deterministic work-window formula (3.6 step 2)
  double-counts nothing with RewritePlan 3.4's static budget, but the
  3.4 validation unit (box_bubble 9/11) should also record the Solve-level
  Diagnostics so miscounts are attributable; needs a one-line hook in the
  M5 ladder runner — owner: Transport/API specs.
