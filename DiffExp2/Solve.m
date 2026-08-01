(* DiffExp2/Solve.m — THE one local solver: symbolic-eps Frobenius with
   EpsSeries coefficient arithmetic.  Spec: Docs/specs/Solve.md (binding);
   decisions: DECISIONS-M0.md.
   True resonance -> explicit log-chain ladder (CASE R); pseudo-resonance ->
   exact Laurent shift (CASE P); regular charts run the SAME recursion
   (CASE T only, V = I).
   CASE-P quotients keep their honest Laurent windows internally, then the
   polar coefficient is cancelled by the already-built target-family
   column (the joint construction of Solve.md 3.7).  The delivered
   multi-sector column is therefore regular as a VALUE while its sector
   representation honestly retains the compensating Laurent pieces. *)

BeginPackage["DiffExp2`Solve`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`", "DiffExp2`CppBackend`"}];

PrepareChart::usage = "PrepareChart[sys, chart] applies the chart map, runs ChartIndicial, and assembles the ChartSystem (theta matrix, gauge, V/VInv spectral frame, families).";
SolveHomogeneous::usage = "SolveHomogeneous[chartSystem, req] gives the FundamentalSystem: one LocalSolution column per indicial sector spec.";
PrewarmHomogeneousBatch::usage = "PrewarmHomogeneousBatch[{chartSystem1, ...}, req] batches the C++ recurrence work for several boundary-independent chart bases into one native task pool, verifies and assembles each result through the ordinary SolveHomogeneous path, and populates its memo cache. All chart systems must belong to one system/configuration; there is no Wolfram fallback.";
HomogeneousCacheCapacity::usage = "HomogeneousCacheCapacity[] gives the bounded number of verified chart bases retained by SolveHomogeneous. It is intended for transport schedulers that must preflight a complete prewarm before submitting any batch.";
SolveParticular::usage = "SolveParticular[chartSystem, source, req] gives THE particular solution (canonical kernel choice) for a sector-native theta-form source supplied in the original physical frame.";
SolveChart::usage = "SolveChart[chartSystem, req, source] gives <|\"Basis\", \"Particular\", \"CouplingDepth\"|>.";
SolveValueRegular::usage = "SolveValueRegular[chartSystem, req, vals] propagates an incoming value vector (one EpsSeries per component at chart center t = 0) through a regular chart with one d-dimensional recursion, without constructing a basis or matching matrix. The delivered epsilon window is capped by the incoming window. Non-regular charts fail loudly.";
SolveNativeValueRegular::usage = "SolveNativeValueRegular[chartSystem, req, vals] propagates a regular center-value vector as one retained C++ local solution. It accepts an SCC skeleton by materializing the exact identity physical frame for this single value recursion, caps delivery to the honest incoming epsilon window, and never returns a coefficient tensor.";
SolveNativeLocalFamily::usage = "SolveNativeLocalFamily[chartSystem, req, <|\"a\"->a,\"b\"->b,\"p\"->p|>, init] runs one uncompensated homogeneous family through the persistent C++ solver and returns an opaque native handle record, never a Wolfram coefficient tensor. init is the same (p+1)-by-d EpsSeries ladder accepted by the framed recurrence. This narrow migration seam requires an identity gauge, grouped native assembly, no unresolved analytic regulators, and no pseudo-resonant family collisions; general transport continues to use SolveHomogeneous/SolveParticular.";
PrepareSCCCouplingMatrix::usage = "PrepareSCCCouplingMatrix[sccChartSystem, sourceBlock, targetBlock, sourceShape, serialization] prepares one exact cross-SCC ThetaOriginal block as a deterministic JSON-ready sparse rational-multiplier matrix. serialization is Automatic (the active C++ serialization Block) or the exact field <|\"domain\"->...,\"symbols\"->{...}|>. Signed epsilon shifts are preserved; execution later proves the requested/work halo contract.";
PrepareNativeRationalRow::usage = "PrepareNativeRationalRow[chartSystem, sourceShape, cvec, physicalVariable, serialization] prepares cvec(center+scale t,eps) as the strict JSON-ready rational row consumed by CppBackend`ApplyPersistentRationalRow. sourceShape supplies the retained local's exact EpsWindow, TWindow, and Dimension; serialization is Automatic or <|\"domain\"->\"acb\"|\"rational\",\"symbols\"->{}|>. The affine dx Jacobian is deliberately not included.";
PrepareNativeSCCComposite::usage = "PrepareNativeSCCComposite[sccChartSystem, req] captures (without executing) the ordinary grouped native homogeneous requests for every supported diagonal SCC block, prepares their strict typed persistent composite manifest with one full original-master physical q/C owner, and returns the opaque C++ SCC handle record. This first slice is an explicit preparation API only; SolveHomogeneous does not dispatch through it.";
SolveNativeSCCBasisColumn::usage = "SolveNativeSCCBasisColumn[sccChartSystem, req, seedBlock, seedLocalComponent:1] executes one strict regular exact-Rational or Acb block-DAG SCC basis column, or a certified exact-Rational/Acb regular-singular Jordan column, through an already captured persistent composite and returns an opaque native local handle record without coefficient tensors. Singular Acb execution is restricted to exact schedules without CASE-P collisions. seedBlock and seedLocalComponent are one-based; the three-argument scalar-v1 call is unchanged. This explicit migration seam is not yet used by SolveHomogeneous or transport.";
SolveNativeSCCBasis::usage = "SolveNativeSCCBasis[sccChartSystem, req, threads:Automatic] executes the complete physical SCC basis as one ordered native column batch, retaining every column atomically and returning opaque handles sorted by physical basis index. No coefficient tensor crosses the bridge.";
SolveNativeRegularBasis::usage = "SolveNativeRegularBasis[chartSystem, req, threads:Automatic, forceMonolithic:False, equationOwner:Automatic] returns a complete retained basis for any regular chart. Multi-block SCC envelopes normally use the ordered native SCC batch; forceMonolithic=True instead solves the full physical frame. When equationOwner is supplied, the framed fallback is published under that exact frame-independent physical owner. No coefficient tensor crosses the bridge.";
PrepareNativeRegularBasisOwner::usage = "PrepareNativeRegularBasisOwner[chartSystem,req] preserves the legacy primitive framed-chart owner path. PrepareNativeRegularBasisOwner[chartSystem,req,anchor] instead prepares a lightweight frame-independent ordinary q/C owner in anchor's existing native session and returns its compact causal value-solver prototype without constructing a recurrence frame.";
WithNativeSCCCompositeCacheReservation::usage = "WithNativeSCCCompositeCacheReservation[count,expr] evaluates expr with capacity reserved for exactly count additional live native SCC composite owners. The reservation is dynamically scoped, never evicts a live public handle, and leaves ordinary direct preparation subject to the default bounded capacity.";
ClearSolveCaches::usage = "ClearSolveCaches[] empties the PrepareChart, exact-SCC-structure, exact-clearing, physical-cleared-ODE, rational-multiplier, SolveHomogeneous, and native SCC composite memo caches, then closes persistent native sessions. Called by API`LoadSystem; the SolveHomogeneous cache additionally self-flushes whenever the chart's SystemHash changes and is entry-capped.";
DropWolframPreparationCaches::usage = "DropWolframPreparationCaches[] drops only Wolfram-side chart/operator/multiplier preparation memo state while preserving every retained native session, chart, SCC, local, match, and tile-plan handle.";
ODEResidualCheck::usage = "ODEResidualCheck[chartSystem, sol, source, probe] checks the theta-form ODE residual at an interior probe point; loud error above ResidTol.";

Begin["`Private`"];

err[id_, cs_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Solve",
    "Chart" -> ToString[Lookup[cs, "Center", "?"], InputForm]|>, payload]];
cfg = DiffExp2`Config`CFG;
numMag = DiffExp2`Tolerances`NumericMagnitude;
(* zeroQ: exact-first.  Together-canonical rational functions over the
   Gaussian rationals in non-numeric symbols are a zero-DECISION domain:
   === 0 decides.  Only forms outside it (inexact numbers, radicals/Root,
   numeric symbol constants like Pi) fall through to PossibleZeroQ, whose
   numeric ztest meprec-storms and can mis-answer on giant exact rationals.
   zeroCanQ takes input the caller already Together/Cancel-canonicalized. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && TrueQ[PossibleZeroQ[c]]);
zeroQ[e_] := zeroCanQ[Together[e]];

esNew = DiffExp2`EpsSeries`ESNew; esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd; esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes; esDiv = DiffExp2`EpsSeries`ESDivide;
esShift = DiffExp2`EpsSeries`ESShift; esCoeff = DiffExp2`EpsSeries`ESCoefficient;
esMin = DiffExp2`EpsSeries`ESMinPower; esCM = DiffExp2`EpsSeries`ESCompleteMax;
esFrom = DiffExp2`EpsSeries`ESFromExpression; esTrim = DiffExp2`EpsSeries`ESTrim;

(* Exact algebraic coefficients are grounded at 2x WP before the framed
   recurrence, but their later cancellations can be stored as 0``A rather
   than exact 0.  A structural first-nonzero scan must not mistake those
   resolved zeros for physical epsilon support: doing so retains dozens of
   adaptive scratch rows and can collapse a later honest Cauchy window.

   Use Chop only as a cheap candidate filter, then require the FULL numeric
   uncertainty enclosure to lie below the configured absolute ChopFloor.
   Thus 0``17 at WP500 is preserved, while 0``1000 is safely exactified.
   Symbolic analytic regulators and all resolved material values pass
   through unchanged. *)
certifiedFrameZeroQ[z_?InexactNumberQ] := Module[
  {floor = DiffExp2`Tolerances`Tol["ChopFloor"]},
  Chop[z, floor] === 0 &&
    TrueQ[Last[DiffExp2`Tolerances`NumericMagnitudeBounds[z, 20]] <= floor]];
certifiedFrameZeroQ[z_] := z === 0;
certifiedFrameChop[z_] := If[certifiedFrameZeroQ[z], 0, z];

(* exact t-Laurent coefficient of a rational function at t = 0 (local copy
   of the Indicial recursion; entries exact) *)
(* contract: p is 0 or a polynomial part of a CANCELED fraction, so the
   syntactic test in zeroCanQ is complete *)
polyMinDeg[p_, t_] := If[zeroCanQ[p], Infinity, Exponent[p, t, Min]];
tVal[e_, t_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    polyMinDeg[Numerator[c], t] - polyMinDeg[Denominator[c], t]]];
tLaurent[e_, t_, k_Integer] := Module[{c, v, num, den, nc, dc, ord, csr},
  c = Cancel[Together[e]];
  If[zeroCanQ[c] || k < (v = tVal[c, t]), Return[0]];
  num = Numerator[c]; den = Denominator[c];
  num = Cancel[num/t^polyMinDeg[num, t]]; den = Cancel[den/t^polyMinDeg[den, t]];
  ord = k - v;
  nc = Table[Coefficient[num, t, j], {j, 0, ord}];
  dc = Table[Coefficient[den, t, j], {j, 0, ord}];
  csr = Table[0, {ord + 1}];
  csr[[1]] = Cancel[Together[nc[[1]]/dc[[1]]]];
  Do[csr[[m + 1]] = Cancel[Together[
      (nc[[m + 1]] - Sum[dc[[j + 1]]*csr[[m - j + 1]], {j, 1, m}])/dc[[1]]]],
    {m, 1, ord}];
  csr[[ord + 1]]];

(* ---- PrepareChart ---- *)

(* A system-level exact clearing pair is expensive but affine-chart
   independent.  Register the exact matrix once and compute the pair lazily
   only if a regular chart asks for it.  The SHA-256 key is accompanied by an
   exact equality check, so a hash collision is loud rather than contaminating
   another system. *)
$pcCache = <||>; $systemClearRegistry = <||>; $globalClearedCache = <||>;
$chartClearedCache = <||>; $exactSCCStructureCache = <||>;
$suppressIntermediateODEResidualChecks = False;
systemClearKey[sys_Association] := Module[{vtag = Replace[
    sys["Variable"], s_Symbol :> {Context[s], SymbolName[s]}]},
  {Hash[{sys["Matrix"], sys["Variable"]}, "SHA256"],
    Dimensions[sys["Matrix"]], Sequence @@ vtag}];
registerSystemClearInput[sys_Association] := Module[
  {key = systemClearKey[sys], old},
  old = If[KeyExistsQ[$systemClearRegistry, key],
    $systemClearRegistry[key], None];
  If[old =!= None &&
      !(old["Variable"] === sys["Variable"] && old["Matrix"] === sys["Matrix"]),
    err["E6", <|"Center" -> "system-clear-registry"|>,
      <|"Key" -> key, "Detail" -> "system clearing cache key collision"|>]];
  If[old === None, AssociateTo[$systemClearRegistry, key ->
    <|"Variable" -> sys["Variable"], "Matrix" -> sys["Matrix"]|>]];
  key];

(* Exact replacement for the old InitializeIntegrationSequence sparsity
   analysis.  Equation row r depends on component c when B[[r,c]] is a
   nonzero element of the full exact coefficient field, so the dependency
   edge is c -> r.  In particular, epsilon- or analytic-regulator-only
   couplings are edges: no specialization is permitted while constructing
   the graph.  Components are returned in a deterministic source-first
   topological order and accompanied by a certificate which is rechecked
   here before any block solve can consume it. *)
deterministicTopologicalOrder[components_List, edges_List] := Module[
  {n = Length[components], indegree, outgoing, ready, order = {}, u},
  indegree = Table[Count[edges, {_, i}], {i, n}];
  outgoing = Table[Sort[Cases[edges, {i, j_} :> j]], {i, n}];
  ready = SortBy[Select[Range[n], indegree[[#]] === 0 &],
    First[components[[#]]] &];
  While[ready =!= {},
    u = First[ready]; ready = Rest[ready]; AppendTo[order, u];
    Do[
      indegree[[v]]--;
      If[indegree[[v]] === 0,
        ready = SortBy[Append[ready, v], First[components[[#]]] &]],
      {v, outgoing[[u]]}]];
  If[Length[order] =!= n,
    err["E6", <|"Center" -> "integration-sequence"|>, <|
      "CondensationEdges" -> edges,
      "Detail" -> "SCC condensation graph is cyclic"|>]];
  order];

(* Private instrumentation seam used by the focused cache contract.  It
   counts actual graph/certificate constructions, not cache lookups. *)
$exactSCCStructureCoreCalls = 0;
exactSCCStructure[m_?MatrixQ] := Module[
  {d = Length[m], normalized, pattern, edgePairs, graph, rawComponents,
   rawVertexToComponent, rawCondensation, topo, oldToNew, components,
   vertexToComponent, condensation, outgoing, depths, certificateQ},
  $exactSCCStructureCoreCalls++;
  If[d === 0 || Dimensions[m] =!= {d, d},
    err["E6", <|"Center" -> "integration-sequence"|>, <|
      "Dimensions" -> Dimensions[m],
      "Detail" -> "SCC analysis requires a nonempty square matrix"|>]];
  If[!FreeQ[m, _?InexactNumberQ],
    err["E1", <|"Center" -> "integration-sequence"|>, <|
      "Detail" -> "SCC decisions require the full exact matrix"|>]];
  normalized = Map[Cancel[Together[#]] &, m, {2}];
  pattern = Map[If[zeroCanQ[#], 0, 1] &, normalized, {2}];
  edgePairs = Flatten[Table[
    If[pattern[[r, c]] === 1, {{c, r}}, {}], {r, d}, {c, d}], 2];
  graph = Graph[Range[d], DirectedEdge @@@ edgePairs,
    DirectedEdges -> True];
  (* On the supported Wolfram kernel ConnectedComponents of a directed graph
     is the strong-component operation (WeaklyConnectedComponents is the
     distinct weak operation).  Its component numbering is not an API
     contract.  Number the condensation
     deterministically after constructing all cross-component edges. *)
  rawComponents = Sort /@ ConnectedComponents[graph];
  rawVertexToComponent = ConstantArray[0, d];
  Do[Scan[(rawVertexToComponent[[#]] = i) &, rawComponents[[i]]],
    {i, Length[rawComponents]}];
  rawCondensation = Sort[DeleteDuplicates[Select[
    ({rawVertexToComponent[[#[[1]]]], rawVertexToComponent[[#[[2]]]]} & /@
      edgePairs), #[[1]] =!= #[[2]] &]]];
  topo = deterministicTopologicalOrder[rawComponents, rawCondensation];
  oldToNew = AssociationThread[topo, Range[Length[topo]]];
  components = rawComponents[[topo]];
  vertexToComponent = oldToNew /@ rawVertexToComponent;
  condensation = Sort[DeleteDuplicates[
    ({oldToNew[#[[1]]], oldToNew[#[[2]]]} & /@ rawCondensation)]];
  outgoing = Table[Sort[Cases[condensation, {i, j_} :> j]],
    {i, Length[components]}];
  (* CouplingDepth counts blocks on the longest dependency chain, matching
     classic MaxCouplingOrder and Transport's decrement heuristic.  The
     one-SCC SolveChart fast path reports 0 separately. *)
  depths = ConstantArray[1, Length[components]];
  Do[Do[depths[[v]] = Max[depths[[v]], depths[[u]] + 1],
      {v, outgoing[[u]]}], {u, Length[components]}];
  certificateQ =
    Sort[Flatten[components]] === Range[d] &&
    Total[Length /@ components] === d &&
    AllTrue[Range[Length[components]], Function[i, Module[
      {sub = Subgraph[graph, components[[i]]], root = First[components[[i]]]},
      Sort[VertexOutComponent[sub, root]] === components[[i]] &&
        Sort[VertexOutComponent[ReverseGraph[sub], root]] ===
          components[[i]]]]] &&
    AllTrue[condensation, #[[1]] < #[[2]] &] &&
    AllTrue[edgePairs, Function[e,
      With[{a = vertexToComponent[[e[[1]]]],
          b = vertexToComponent[[e[[2]]]]},
        a === b || MemberQ[condensation, {a, b}]]]];
  If[!TrueQ[certificateQ],
    err["E6", <|"Center" -> "integration-sequence"|>, <|
      "Components" -> components, "CondensationEdges" -> condensation,
      "Detail" -> "exact SCC/topological certificate failed"|>]];
  <|"Schema" -> "DiffExp2.IntegrationSequence/v1",
    "Exact" -> True, "RegulatorSpecialization" -> None,
    "MatrixHash" -> Hash[normalized, "SHA256"],
    "PatternHash" -> Hash[pattern, "SHA256"],
    "NonzeroPattern" -> pattern, "DependencyEdges" -> edgePairs,
    "Components" -> components,
    "VertexToComponent" -> vertexToComponent,
    "CondensationEdges" -> condensation,
    "TopologicalOrder" -> Range[Length[components]],
    "BlockDepths" -> depths,
    "CouplingDepth" -> If[depths === {}, 0, Max[depths]],
    "Certificate" -> <|"Partition" -> True,
      "StrongConnectivity" -> True, "AcyclicTopologicalOrder" -> True,
      "EdgeCoverage" -> True|>|>];

(* An invertible affine substitution x -> x0 + beta t and multiplication by
   the nonzero monomial beta t are injective on the exact rational-function
   field.  They therefore preserve the entrywise structural-zero pattern.
   Build the SCC graph/certificate once from the original exact system and
   reuse it at every chart of that system/level.  The SHA key is only an
   index: registerSystemClearInput has just rechecked the complete exact input
   against $systemClearRegistry, so a hash collision is loud without doing a
   second full-matrix comparison here.  MatrixHash is chart-local metadata
   and is replaced by prepareChartCore after ThetaOriginal has been formed. *)
exactSCCStructureForSystem[sys_Association, key_] := Module[
  {structure},
  If[!KeyExistsQ[$systemClearRegistry, key],
    err["E6", <|"Center" -> "integration-sequence-cache"|>, <|
      "Key" -> key,
      "Detail" -> "exact SCC structure key was not collision-certified"|>]];
  If[KeyExistsQ[$exactSCCStructureCache, key],
    Return[$exactSCCStructureCache[key], Module]];
  structure = exactSCCStructure[sys["Matrix"]];
  AssociateTo[$exactSCCStructureCache, key -> structure];
  structure];

familyCollisionDepth[roots_List, collisions_List] := Module[{perRoot},
  perRoot = Table[Module[{events, layers},
    (* Homogeneous recursion starts at n=1; same-a n=0 symmetrization is
       diagnostic degeneracy, not an executed division layer. *)
    events = Select[collisions,
      #["Type"] === "LaurentShift" && #["LowerIdx"] === i && #["n"] >= 1 &];
    layers = GatherBy[events, #["n"] &];
    Total[Map[Max[roots[[#[[All, "UpperIdx"]]]][[All, "BlockSize"]]] &, layers]]],
    {i, Length[roots]}];
  If[perRoot === {}, 0, Max[perRoot]]];

(* Jordan chains are defined only up to one nonzero scalar in eps per
   block.  Indicial's deterministic first-entry normalization can make that
   scalar meromorphic when several affine eigenvalues coalesce at eps=0.
   A CASE-P certificate is a statement about poles CREATED by the collision;
   it must not be asked to cancel an arbitrary pole already present in the
   chosen seed vector.  Put every block in its epsilon-primitive lattice by
   clearing the common lowest valuation.  The same monomial multiplies every
   member of a Jordan chain, so (R-lambda I)v_j=v_{j-1} and the unit Jordan
   superdiagonal are unchanged.  Different blocks may be scaled
   independently because the spectral Jordan matrix is block diagonal. *)
epsilonPrimitiveJordanChain[chain_List, eps_Symbol] := Module[
  {vals, shift, normalized, normalizedVals},
  vals = Select[tVal[#, eps] & /@ Flatten[chain], IntegerQ];
  If[vals === {},
    err["E2", <|"Center" -> "spectral-normalization"|>,
      <|"Detail" -> "Jordan chain has no nonzero entry"|>]];
  shift = Max[0, -Min[vals]];
  normalized = Map[Cancel[Together[eps^shift*#]] &, chain, {2}];
  normalizedVals = Select[tVal[#, eps] & /@ Flatten[normalized], IntegerQ];
  If[AnyTrue[normalizedVals, # < 0 &] || Min[normalizedVals] =!= 0,
    err["E2", <|"Center" -> "spectral-normalization"|>,
      <|"Shift" -> shift, "Valuations" -> normalizedVals,
        "Detail" -> "epsilon-primitive Jordan block invariant failed"|>]];
  {normalized, shift}];

PrepareChart[sys_Association, chart_Association] := Module[
  {sysClearKey = registerSystemClearInput[sys], pcKey},
  pcKey = {sysClearKey, chart["Center"], Lookup[chart, "Scale", 1],
    Lookup[chart, "Radius", None], Lookup[chart, "LocalRadius", None],
    Lookup[chart, "Prescriptions", {}],
    TrueQ[Lookup[chart, "UseSCCSkeleton", False]]};
  If[KeyExistsQ[$pcCache, pcKey], Return[$pcCache[pcKey]]];
  $pcCache[pcKey] = prepareChartCore[sys, chart, sysClearKey]];

prepareChartCore[sys_Association, chart_Association, sysClearKey_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t, x0, beta, A, Achart, idata, d,
   cols, V, VInv, fams, detV, colIdx, blockEpsShifts,
   integrationSequence, thetaOriginal, regularOriginal, apparent,
   baseReduction, composedGauge, composedGaugeInverse},
  t = chart["ChartVar"]; x0 = chart["Center"]; beta = Lookup[chart, "Scale", 1];
  A = sys["Matrix"];
  Achart = Map[Cancel[Together[#]] &,
    beta*(A /. sys["Variable"] -> x0 + beta*t), {2}];
  thetaOriginal = Map[Cancel[Together[t*#]] &, Achart, {2}];
  If[zeroCanQ[Cancel[Together[beta]]],
    err["E1", chart, <|
      "Detail" -> "affine chart Scale must be nonzero for SCC analysis"|>]];
  (* Preserve the previous exact-input contract.  The reusable graph is
     proved from the exact original system, while ThetaOriginal remains the
     chart-local matrix used by recurrences and residual checks. *)
  If[!FreeQ[thetaOriginal, _?InexactNumberQ],
    err["E1", chart, <|
      "Detail" -> "SCC decisions require the full exact chart matrix"|>]];
  (* The integration sequence belongs to the original master basis.  A
     full-system rank-reduction gauge or spectral frame may mix components
     and must not be allowed to densify/change the physical dependency DAG.
     Every diagonal SCC is prepared independently later. *)
  integrationSequence = Join[
    exactSCCStructureForSystem[sys, sysClearKey],
    <|"MatrixHash" -> Hash[thetaOriginal, "SHA256"]|>];
  d = Length[Achart];
  If[TrueQ[Lookup[chart, "UseSCCSkeleton", False]] &&
      Length[integrationSequence["Components"]] > 1,
    (* A full spectral preparation before preparing every diagonal SCC defeats
       the integration sequence (the double-box full ChartIndicial alone is
       minute-scale).  Transport and the SCC solver need only the original
       exact theta matrix, chart metadata, dimension, and the ordinary/singular
       classification.  Each diagonal block below receives its own complete
       PrepareChart, including any rank reduction and spectral frame. *)
    regularOriginal = AllTrue[Flatten[Achart],
      With[{v = tVal[#, t]}, v === Infinity || TrueQ[v >= 0]] &] &&
      !TrueQ[DiffExp2`Indicial`EpsilonCoalescingDenominatorQ[
        Achart, t, eps]];
    Return[<|"ChartVar" -> t, "Center" -> x0,
      "ChartMap" -> <|"Center" -> x0, "Scale" -> beta|>,
      "Radius" -> Lookup[chart, "LocalRadius", chart["Radius"]],
      "SystemHash" -> Hash[sys["Matrix"]],
      "SystemClearKey" -> sysClearKey,
      "Prescriptions" -> Lookup[chart, "Prescriptions", {}],
      "SystemSize" -> d, "ThetaOriginal" -> thetaOriginal,
      "IntegrationSequence" -> integrationSequence,
      "IndicialData" -> <|"Regular" -> regularOriginal,
        "Dimension" -> d, "SCCSkeleton" -> True|>,
      "SCCSkeleton" -> True,
      "ChartSystemKind" -> "SCCEnvelope"|>, Module]];
  apparent = DiffExp2`Indicial`EpsilonCoalescingApparentReduce[
    Achart, t, eps,
    <|"Name" -> Lookup[chart, "Name", "chart@" <> ToString[x0, InputForm]],
      "Center" -> x0, "Variable" -> sys["Variable"]|>];
  idata = DiffExp2`Indicial`ChartIndicial[apparent["Matrix"], t, eps,
    <|"Name" -> Lookup[chart, "Name", "chart@" <> ToString[x0, InputForm]],
      "Center" -> x0, "Variable" -> sys["Variable"]|>];
  If[TrueQ[apparent["Applied"]],
    baseReduction = idata["Reduction"];
    composedGauge = Map[Cancel[Together[#]] &,
      apparent["Gauge"] . baseReduction["Gauge"], {2}];
    composedGaugeInverse = Map[Cancel[Together[#]] &,
      baseReduction["GaugeInverse"] . apparent["GaugeInverse"], {2}];
    If[!AllTrue[Flatten[Map[Cancel[Together[#]] &,
          composedGauge . composedGaugeInverse -
            IdentityMatrix[Length[Achart]], {2}]], # === 0 &] ||
        !AllTrue[Flatten[Map[Cancel[Together[#]] &,
          composedGaugeInverse . composedGauge -
            IdentityMatrix[Length[Achart]], {2}]], # === 0 &],
      err["E6", chart, <|
        "Stage" -> "EpsilonCoalescingApparentGaugeComposition",
        "Detail" -> "composed chart gauge is not an exact two-sided inverse"|>]];
    baseReduction = Join[baseReduction, <|
      "Gauge" -> composedGauge,
      "GaugeInverse" -> composedGaugeInverse,
      "GaugeInverseCertified" -> True,
      "GaugeInverseCertificateSchema" ->
        "diffexp2-indicial-exact-gauge-inverse-v1",
      "ApparentReduction" -> KeyDrop[apparent, {
        "Matrix", "Gauge", "GaugeInverse"}],
      "ReductionMethod" -> StringRiffle[
        DeleteCases[{
          "EpsilonCoalescingProjectorGauge",
          Lookup[baseReduction, "ReductionMethod", None]}, None], "+"]|>];
    idata = Join[idata, <|"Reduction" -> baseReduction|>]];
  d = idata["Dimension"];
  cols = {}; fams = {}; colIdx = 0; blockEpsShifts = {};
  Do[Module[{members, roots = {}, colRange0 = colIdx + 1, collisions = {}},
    members = Reverse[SortBy[fam["Members"], {#["a"], #["b"]} &]];
    Do[
      Do[Module[{chain = mem["Chains"][[ci]], normalized},
        normalized = epsilonPrimitiveJordanChain[chain, eps];
        chain = normalized[[1]];
        AppendTo[blockEpsShifts, normalized[[2]]];
        AppendTo[roots, <|"a" -> mem["a"], "b" -> mem["b"],
          "BlockSize" -> Length[chain]|>];
        Do[AppendTo[cols, chain[[q]]]; colIdx++, {q, Length[chain]}]],
        {ci, Length[mem["Chains"]]}],
      {mem, members}];
    Do[Module[{ri = roots[[i]], rj = roots[[j]], nn},
      If[i =!= j,
        nn = RootReduce[Together[rj["a"] - ri["a"]]];
        If[IntegerQ[nn] && nn >= 0 &&
            !(zeroQ[ri["a"] - rj["a"]] && zeroQ[ri["b"] - rj["b"]]),
          AppendTo[collisions, <|"n" -> nn, "LowerIdx" -> i, "UpperIdx" -> j,
            "DeltaB" -> RootReduce[Together[ri["b"] - rj["b"]]],
            "Type" -> If[zeroQ[ri["b"] - rj["b"]], "Log", "LaurentShift"]|>]]]],
      {i, Length[roots]}, {j, Length[roots]}];
    AppendTo[fams, <|"Roots" -> roots,
      "ColumnRange" -> Range[colRange0, colIdx],
      "Collisions" -> collisions,
      "CollisionDepth" -> familyCollisionDepth[roots, collisions]|>]],
    {fam, idata["Families"]}];
  V = Transpose[cols];
  detV = Together[Det[V]];
  If[zeroCanQ[detV], err["E2", chart, <|"Detail" -> "spectral frame V is singular"|>]];
  VInv = Map[Cancel[Together[#]] &, Inverse[V], {2}];
  <|"ChartVar" -> t, "Center" -> x0,
    "ChartMap" -> <|"Center" -> x0, "Scale" -> beta|>,
    (* LocalSolution evaluation is in t, not the physical line variable. *)
    "Radius" -> Lookup[chart, "LocalRadius", chart["Radius"]],
    "SystemHash" -> Hash[sys["Matrix"]],   (* solve-cache flush tag *)
    "SystemClearKey" -> sysClearKey,
    "Prescriptions" -> Lookup[chart, "Prescriptions", {}],
    "SystemSize" -> d,
    "ThetaMatrix" -> idata["Reduction"]["ThetaMatrix"],
    (* original theta matrix t*A(chart): the residual check runs on the
       GAUGED-BACK solution f = T.g, which solves the ORIGINAL system -
       checking it against the reduced B would fail every gauge chart *)
    "ThetaOriginal" -> thetaOriginal,
    "Gauge" -> idata["Reduction"]["Gauge"],
    "GaugeInverse" -> idata["Reduction"]["GaugeInverse"],
    "IntegrationSequence" -> integrationSequence,
    "Residue" -> idata["Residue"],
    "SpectralBlockEpsShifts" -> blockEpsShifts,
    "V" -> V, "VInv" -> VInv, "Families" -> fams,
    "IndicialData" -> idata|>];

(* block list in column order: one entry per Jordan block *)
blockList[cs_] := Flatten[Map[Function[fam,
  Module[{off = First[fam["ColumnRange"]] - 1, bl = {}},
    Do[AppendTo[bl, <|"a" -> r["a"], "b" -> r["b"], "q" -> r["BlockSize"],
        "Cols" -> Range[off + 1, off + r["BlockSize"]]|>];
      off += r["BlockSize"],
      {r, fam["Roots"]}];
    bl]], cs["Families"]], 1];

(* log ceiling for tag (a, b): own p0 + same-b block sizes at higher
   integer-spaced a (the generalized MaxLogPowers formula) *)
logCeiling[cs_, aT_, bT_, p0_, includeEqual_:False] := p0 + Total[Map[
  Module[{dd = RootReduce[Together[#["a"] - aT]]},
    If[zeroQ[#["b"] - bT] && IntegerQ[dd] &&
        (TrueQ[dd > 0] || (includeEqual && dd === 0)), #["q"], 0]] &,
  blockList[cs]]];

(* CASE-P divisions at one Taylor order act on disjoint spectral blocks in
   parallel; only distinct increasing n layers compose.  A size-q Jordan
   target costs q epsilon orders at that layer. *)
pseudoDepthForTag[cs_, aT_, bT_, nStart_Integer] := Module[{events, groups},
  events = Select[blockList[cs], Function[blk, Module[{nn =
      RootReduce[Together[blk["a"] - aT]]},
    IntegerQ[nn] && nn >= nStart && !zeroQ[bT - blk["b"]]]]];
  groups = GatherBy[events,
    RootReduce[Together[#["a"] - aT]] &];
  Total[Map[Max[#[[All, "q"]]] &, groups]]];

(* ---- the recursion ---- *)

(* ===================== PACKED-FRAME RECURSION CORE =====================
   Coefficient vectors are PLAIN LISTS over a fixed eps-frame
   [fb, fb+W-1] (index i <-> eps^(fb+i-1)); products via ListConvolve.
   This bypasses the per-op EpsSeries object layer (Association churn +
   per-coefficient Together) that cost ~18s/chart on 3-master systems.
   Honest windows: the frame is sized by the deterministic work-window
   formula; a per-(n,log,component) integer shadow follows completeness
   through the actual dependency DAG.  EpsSeries objects appear only at
   the boundaries (sources in, LocalSolutions out). *)

(* fast exact eps-expansion of a RATIONAL expr into a frame list;
   ByteCount-gated numericization at 2x WP *)
ratEpsList[expr_, eps_, fb_, W_] := Module[
  {c, num, den, vn, vd, v, rel, out, top},
  out = ConstantArray[0, W];
  c = Cancel[Together[expr]];
  If[zeroCanQ[c], Return[out]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min];
  v = vn - vd;
  If[v < fb, err["E4", <|"Center" -> "(frame)"|>,
    <|"Valuation" -> v, "FrameBase" -> fb,
      "Detail" -> "eps-valuation below the work frame (buffer formula violated)"|>]];
  num = num/eps^vn // Cancel; den = den/eps^vd // Cancel;
  top = fb + W - 1;
  If[v > top, err["E4", <|"Center" -> "(frame)"|>,
    <|"Valuation" -> v, "FrameTop" -> top,
      "Detail" -> "eps-valuation above the work frame (buffer formula violated)"|>]];
  rel = top - v;
  Module[{ncl = Table[Coefficient[num, eps, j], {j, 0, rel}],
    dcl = Table[Coefficient[den, eps, j], {j, 0, rel}], csr},
    csr = ConstantArray[0, rel + 1];
    csr[[1]] = Together[ncl[[1]]/dcl[[1]]];
    Do[csr[[m + 1]] = Together[
        (ncl[[m + 1]] - Sum[dcl[[j + 1]]*csr[[m - j + 1]], {j, 1, m}])/dcl[[1]]],
      {m, 1, rel}];
    csr = Map[preparedEpsCoefficient, csr];
    Do[out[[v - fb + 1 + m]] = csr[[m + 1]], {m, 0, rel}]];
  out];

(* Adaptive singular homogeneous solves use a deliberately narrow lower
   frame first.  A genuine lower underflow throws this private retry signal;
   the same recurrence is then rerun on a wider rectangle.  Outside that
   narrowly scoped Block, the established public E4 contract is unchanged. *)
$adaptiveLowerFrameProbe = False;
$cancellationAuditBase = None;
$disableAdaptiveLowerFrames = False;   (* private parity/benchmark seam *)
$disableRationalDenominatorFusion = False;  (* private parity/benchmark seam *)
$disableGroupedSpectralTransform = False;   (* private parity/benchmark seam *)
lowerFrameUnderflow[cs_, payload_Association] :=
  If[TrueQ[$adaptiveLowerFrameProbe],
    Throw[Failure["AdaptiveLowerFrame", Join[
      <|"ID" -> "AdaptiveLowerFrame"|>, payload]],
      "DiffExp2AdaptiveLowerFrame"],
    err["E4", cs, payload]];

(* One (d x W) epsilon-frame shift.  This is shared by the ordinary sparse
   polynomial path and the grouped rational path below, so both retain the
   same strict lower-frame contract. *)
shiftFrameBlock[M_, s_Integer, fb_Integer, W_Integer, cs_] := Which[
  s === 0, M,
  s > 0 && s < W,
    Join[ConstantArray[0, {Length[M], s}], M[[All, 1 ;; W - s]], 2],
  s >= W, ConstantArray[0, Dimensions[M]],
  -s < W,
    If[AnyTrue[Flatten[M[[All, 1 ;; -s]]], # =!= 0 &],
      lowerFrameUnderflow[cs, <|"FrameBase" -> fb, "Shift" -> s,
        "Detail" -> "epsilon shift would discard nonzero lower-frame content"|>]];
    Join[M[[All, -s + 1 ;; W]], ConstantArray[0, {Length[M], -s}], 2],
  True,
    If[AnyTrue[Flatten[M], # =!= 0 &],
      lowerFrameUnderflow[cs, <|"FrameBase" -> fb, "Shift" -> s,
        "Detail" -> "epsilon shift lies wholly below the work frame"|>]];
    ConstantArray[0, Dimensions[M]]];

matrixShiftProduct[A_, M_, s_Integer, fb_Integer, W_Integer, cs_] :=
  If[TrueQ[$adaptiveLowerFrameProbe] && s < 0,
    shiftFrameBlock[A . M, s, fb, W, cs],
    A . shiftFrameBlock[M, s, fb, W, cs]];

exactNonRationalNumberQ[e_] := NumericQ[e] && Precision[e] === Infinity &&
  !IntegerQ[e] && Head[e] =!= Rational;

$numericizeAllPreparedNumbers = False;
$disablePreparedDirectNumericization = False;  (* private parity seam *)
$preparedDirectGuardDigits = 32;
preparedNumericizationRequiredQ[e_] :=
  exactNonRationalNumberQ[e] ||
  (TrueQ[$numericizeAllPreparedNumbers] && e =!= 0 && NumericQ[e] &&
    Precision[e] === Infinity);

groupedEpsExactSafeQ[e_] := e === 0 ||
  (* A Rational native session is an exact coefficient-field computation,
     not an Acb computation whose decimal inputs happen to have zero radius.
     In particular, the CASE-P Rational shadow is prepared under
     $cppExactDomain after the production Acb owner has already been built.
     Preserve every exact Q coefficient there regardless of expression
     size; reconstructing a rounded 2x-WP Real later would be neither exact
     nor compatible with the shadow's provenance contract. *)
  (TrueQ[$cppExactDomain] &&
    (IntegerQ[e] || Head[e] === Rational)) ||
  (ByteCount[e] <= 500 && !preparedNumericizationRequiredQ[e]);

finitePreparedInexactScalarQ[value_, requiredPrecision_] :=
  InexactNumberQ[value] && NumberQ[value] &&
    FreeQ[value, Indeterminate | ComplexInfinity | DirectedInfinity[_]] &&
    NumberQ[Precision[value]] && Precision[value] >= requiredPrecision;

(* Exact planning has already classified the epsilon valuation, P/Q grouping,
   spectral data and frame bounds before coefficients reach this boundary.
   For the production Acb owner, try numerical evaluation of the original
   exact scalar first: explicitly canonicalizing thousands of distinct
   algebraic coefficients with RootReduce can dominate singular-chart
   preparation.  A direct result is only authoritative when it retains the
   full requested precision and is nonzero.  Zero or cancellation-degraded
   results fall back to exact canonicalization below. *)
preparedDirectEpsCoefficient[e_, wp2_Integer] := Module[{candidate},
  candidate = Quiet[Check[
    Block[{$MaxExtraPrecision = Max[$MaxExtraPrecision,
        DiffExp2`Tolerances`$MaxExtraPrecisionValue, 2 wp2]},
      N[e, wp2 + $preparedDirectGuardDigits]], $Failed]];
  If[candidate === $Failed ||
      !finitePreparedInexactScalarQ[candidate, wp2] ||
      TrueQ[candidate == 0],
    $Failed,
    N[candidate, wp2]]];

preparedEpsCoefficient[e_] := Module[{wp2, direct, canonical, grounded},
  If[groupedEpsExactSafeQ[e], Return[e]];
  If[!NumericQ[e], Return[e]];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*
    DiffExp2`Config`CFG["WorkingPrecision"];
  If[!TrueQ[$cppExactDomain] &&
      !TrueQ[$disablePreparedDirectNumericization],
    direct = preparedDirectEpsCoefficient[e, wp2];
    If[direct =!= $Failed, Return[direct]]];
  canonical = RootReduce[e];
  If[canonical === 0, Return[0]];
  Block[{$MaxExtraPrecision = Max[$MaxExtraPrecision,
      DiffExp2`Tolerances`$MaxExtraPrecisionValue, 2 wp2]},
    grounded = N[canonical, wp2 + $preparedDirectGuardDigits]];
  (* Never manufacture precision with SetPrecision, and never Chop an exact
     coefficient which exact canonicalization proved nonzero.  EncodeScalar
     will carry the value's actual Accuracy into the Arb input radius. *)
  If[finitePreparedInexactScalarQ[grounded, wp2] &&
      !TrueQ[grounded == 0],
    N[grounded, wp2],
    canonical]];

(* Canonical epsilon-rational entry.  eps^Valuation is kept outside P/Q;
   Q is normalized to Q(0)=1 so denominator equality is an exact grouping
   key and coefficient division is causal on the finite frame. *)
epsRationalData[expr_, eps_, fb_Integer, top_Integer, cs_] := Module[
  {c, num, den, vn, vd, v, num0, den0, q0, p, q, qd},
  c = Cancel[Together[expr]];
  If[zeroCanQ[c], Return[<|"Zero" -> True, "Valuation" -> Infinity|>]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min];
  v = vn - vd;
  If[v < fb,
    err["E4", cs, <|"Valuation" -> v, "FrameBase" -> fb,
      "Detail" -> "eps-valuation below the work frame (buffer formula violated)"|>]];
  If[v > top,
    err["E4", cs, <|"Valuation" -> v, "FrameTop" -> top,
      "Detail" -> "eps-valuation above the work frame (buffer formula violated)"|>]];
  num0 = Cancel[num/eps^vn]; den0 = Cancel[den/eps^vd];
  q0 = Coefficient[den0, eps, 0];
  If[zeroCanQ[q0],
    err["E5", cs, <|"Denominator" -> den0,
      "Detail" -> "epsilon-unit denominator has zero constant coefficient"|>]];
  p = Expand[Cancel[Together[num0/q0]]];
  q = Expand[Cancel[Together[den0/q0]]];
  If[!PolynomialQ[p, eps] || !PolynomialQ[q, eps],
    err["E5", cs, <|"Expression" -> c,
      "Detail" -> "prepared epsilon-rational entry is not polynomial-over-polynomial"|>]];
  qd = Exponent[q, eps];
  <|"Zero" -> False, "Valuation" -> v, "P" -> p, "Q" -> q,
    "DenominatorDegree" -> qd|>];

(* Frame-independent epsilon-rational data for the physical cleared ODE.
   Unlike recurrence preparation this payload is an exact equation owner: it
   must not acquire upper-frame boundary corrections or reject a valuation
   merely because one particular local solve requested a narrower slab. *)
physicalEpsRationalData[expr_, eps_, cs_] := Module[
  {c, num, den, vn, vd, v, num0, den0, q0, p, q},
  c = Cancel[Together[expr]];
  If[zeroCanQ[c], Return[<|"Zero" -> True|>]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min];
  v = vn - vd;
  If[!IntegerQ[v],
    err["E5", cs, <|"Expression" -> c,
      "Detail" -> "physical ODE coefficient has a nonintegral epsilon valuation"|>]];
  num0 = Cancel[num/eps^vn]; den0 = Cancel[den/eps^vd];
  q0 = Coefficient[den0, eps, 0];
  If[zeroCanQ[q0],
    err["E5", cs, <|"Denominator" -> den0,
      "Detail" -> "physical ODE epsilon denominator is not a causal unit"|>]];
  p = Expand[Cancel[Together[num0/q0]]];
  q = Expand[Cancel[Together[den0/q0]]];
  If[!PolynomialQ[p, eps] || !PolynomialQ[q, eps] ||
      zeroCanQ[Coefficient[p, eps, 0]] ||
      !TrueQ[Coefficient[q, eps, 0] === 1],
    err["E5", cs, <|"Expression" -> c,
      "Detail" -> "physical ODE coefficient did not normalize to eps^v P/Q with P(0)!=0 and Q(0)=1"|>]];
  <|"Zero" -> False, "Valuation" -> v, "P" -> p, "Q" -> q|>];

$physicalClearedODECache = <||>;
$physicalClearedODECacheMax = 256;

(* Capture the exact equation in the delivered master basis, before any
   spectral V/VInv recurrence frame is applied:

                    q(t,eps) theta f = C(t,eps) f,
                    C = q ThetaOriginal.

   The full exact signature is retained beside its hash cache index, so a
   collision can never substitute another physical equation. *)
physicalClearedODEData[cs_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], theta,
   useGlobalAffineClear, signature, key, cached, symbolic, den, denCoeffs,
   denContent, cMatrix, qDegree, cDegree, qCoeffs, cCoeffs, qData, cData,
   identity, data},
  theta = Lookup[cs, "ThetaOriginal", None];
  If[!MatrixQ[theta] || Dimensions[theta] =!=
      {cs["SystemSize"], cs["SystemSize"]},
    err["E5", cs, <|
      "Detail" -> "physical residual ownership requires the exact original-master theta matrix"|>]];
  (* Regular recurrence preparation already proves and caches the exact
     global clearing

                         D(x,eps) f'(x) = N(x,eps) f(x)

     once per input system.  The physical owner used to ignore that proof
     and recompute an LCM/GCD from every affine copy

                   t beta A(x0 + beta t),

     which is mathematically the same clearing but can become catastrophically
     expensive as exact rational chart centers acquire large denominators.
     Affine substitution is an automorphism of the coefficient field, so use
     the existing global-clear path for an identity-frame regular chart and
     convert its already-certified dExpr/NhatExpr coefficients directly to
     the physical q/C record.  Singular/nonidentity and standalone fixture
     records retain the legacy local clearing below. *)
  useGlobalAffineClear =
    !TrueQ[$disableGlobalClearedHoist] &&
    KeyExistsQ[cs, "SystemClearKey"] &&
    KeyExistsQ[cs, "Center"] &&
    AssociationQ[Lookup[cs, "ChartMap", None]] &&
    KeyExistsQ[cs["ChartMap"], "Scale"];
  signature = If[useGlobalAffineClear,
    {"physical-global-affine-clear-v2", cs["SystemClearKey"], t,
      cs["Center"], cs["ChartMap", "Scale"]},
    {t, theta}];
  key = Hash[signature, "SHA256"];
  cached = Lookup[$physicalClearedODECache, key, None];
  If[AssociationQ[cached] && SameQ[cached["Signature"], signature],
    Return[cached["Data"], Module]];
  If[cached =!= None,
    err["E5", cs, <|
      "Detail" -> "physical cleared-ODE cache-key collision with unequal exact input"|>]];
  If[useGlobalAffineClear,
    (* q(t) theta = C(t) is captured in the original-master frame.  Its
       affine-global construction is valid at ordinary and singular centers;
       only recurrence preparation needs q(0,eps=0) != 0. *)
    symbolic = affinePhysicalClearedFromGlobal[cs];
    qCoeffs = symbolic["dExpr"];
    cCoeffs = symbolic["NhatExpr"],
    den = Together[PolynomialLCM @@
      (Denominator[Together[#]] & /@ Flatten[theta])];
    denCoeffs = Select[CoefficientList[den, t], !zeroCanQ[#] &];
    If[denCoeffs === {},
      err["E5", cs, <|
        "Detail" -> "physical cleared-ODE denominator is identically zero"|>]];
    denContent = If[Length[denCoeffs] === 1, First[denCoeffs],
      Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]];
    den = Cancel[Together[den/denContent]];
    cMatrix = Map[Cancel[Together[den*#]] &, theta, {2}];
    If[!PolynomialQ[den, t] ||
        !AllTrue[Flatten[cMatrix], PolynomialQ[#, t] &],
      err["E5", cs, <|
        "Detail" -> "physical cleared ODE is not polynomial in the local chart coordinate"|>]];
    qDegree = Max[0, Exponent[den, t]];
    cDegree = Max[0, Max[Exponent[#, t] & /@ Flatten[cMatrix]]];
    qCoeffs = Table[Cancel[Together[Coefficient[den, t, j]]],
      {j, 0, qDegree}];
    cCoeffs = Table[Map[Cancel[Together[#]] &,
        Map[Coefficient[#, t, j] &, cMatrix, {2}], {2}],
      {j, 0, cDegree}]];
  qData = physicalEpsRationalData[#, eps, cs] & /@ qCoeffs;
  cData = Map[physicalEpsRationalData[#, eps, cs] &, cCoeffs, {3}];
  identity = "de2-physical-ode-" <>
    IntegerString[Hash[{"physical-original-master", qData, cData},
      "SHA256"], 16, 64];
  data = <|"Identity" -> identity, "Q" -> qData, "C" -> cData|>;
  If[Length[$physicalClearedODECache] >= $physicalClearedODECacheMax,
    KeyDropFrom[$physicalClearedODECache,
      First[Keys[$physicalClearedODECache]]]];
  AssociateTo[$physicalClearedODECache, key -> <|
    "Signature" -> signature, "Data" -> data|>];
  data];

physicalEpsRationalExpression[data_Association, eps_, cs_] := Module[
  {zero = Lookup[data, "Zero", None], valuation, p, q},
  If[TrueQ[zero], Return[0, Module]];
  valuation = Lookup[data, "Valuation", None];
  p = Lookup[data, "P", None];
  q = Lookup[data, "Q", None];
  If[zero =!= False || !IntegerQ[valuation] ||
      !PolynomialQ[p, eps] || !PolynomialQ[q, eps] ||
      zeroCanQ[Coefficient[p, eps, 0]] ||
      !TrueQ[Coefficient[q, eps, 0] === 1],
    err["E6", cs, <|"PhysicalCoefficient" -> data,
      "Detail" ->
        "retained regular equation contains a malformed epsilon-rational coefficient"|>]];
  Cancel[Together[eps^valuation p/q]]];

(* A frame-independent regular owner has already paid the exact polynomial
   LCM/GCD clear needed to publish q(t,eps) theta f = C(t,eps) f.  Rebuild the
   ordinary identity-frame recurrence directly from that immutable q/C
   record instead of clearing the same chart a second time when value
   transport falls back to a basis solve. *)
regularClearedSymbolicFromPhysicalData[cs_Association,
    data_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], d = cs["SystemSize"],
   qData, cData, expectedIdentity, dExpr, nExpr},
  qData = Lookup[data, "Q", None];
  cData = Lookup[data, "C", None];
  expectedIdentity = If[ListQ[qData] && ListQ[cData],
    "de2-physical-ode-" <> IntegerString[
      Hash[{"physical-original-master", qData, cData}, "SHA256"],
      16, 64], None];
  If[!ListQ[qData] || qData === {} || !ListQ[cData] || cData === {} ||
      !AllTrue[qData, AssociationQ] ||
      !AllTrue[cData, MatrixQ[#] && Dimensions[#] === {d, d} &&
          AllTrue[Flatten[#], AssociationQ] &] ||
      Lookup[data, "Identity", None] =!= expectedIdentity,
    err["E6", cs, <|"PhysicalData" -> KeyTake[data, {"Identity"}],
      "Detail" ->
        "retained regular equation q/C data failed its exact content identity or dimensions"|>]];
  dExpr = physicalEpsRationalExpression[#, eps, cs] & /@ qData;
  nExpr = Map[physicalEpsRationalExpression[#, eps, cs] &,
    cData, {3}];
  <|"dExpr" -> dExpr, "NhatExpr" -> nExpr,
    "dD" -> Length[dExpr] - 1, "dN" -> Length[nExpr] - 1|>];

(* Q=1 entries need no second Together/series-division pass after
   epsRationalData has already canonicalized them. *)
polynomialEpsListFromData[data_Association, eps_, fb_Integer, W_Integer] :=
 Module[{out = ConstantArray[0, W], v = data["Valuation"], pCoeffs, kmax},
  pCoeffs = CoefficientList[data["P"], eps];
  kmax = Min[Length[pCoeffs] - 1, fb + W - 1 - v];
  Do[out[[v - fb + 1 + k]] = preparedEpsCoefficient[pCoeffs[[k + 1]]],
    {k, 0, kmax}];
  out];

(* Let Mbar be the Laurent series of one rational entry truncated exactly to
   [fb, top], as ratEpsList stores it.  Instead of storing all O(W) entries of
   Mbar, store B = Q Mbar.  B consists of the finite numerator plus at most
   deg(Q) boundary coefficients above top.  Causal division by Q after the
   matrix application therefore reproduces Mbar.U, including the old upper
   cutoff when U has negative epsilon powers. *)
rationalBoundaryData[rec_List, qCoeffs_List, eps_, top_Integer] := Module[
  {v = rec[[4]], p = rec[[5]], qd = Length[qCoeffs] - 1, rel,
   pCoeffs, csr, coeffAt, low, boundary, n0, reason},
  rel = top - v;
  pCoeffs = CoefficientList[p, eps];
  csr = ConstantArray[0, rel + 1];
  Do[
    n0 = If[k + 1 <= Length[pCoeffs], pCoeffs[[k + 1]], 0];
    csr[[k + 1]] = Together[(n0 - Sum[
      qCoeffs[[i + 1]]*csr[[k - i + 1]], {i, 1, Min[k, qd]}])/
        qCoeffs[[1]]],
    {k, 0, rel}];
  coeffAt[k_Integer] := If[0 <= k <= rel, csr[[k + 1]], 0];
  low = Table[{v + k, pCoeffs[[k + 1]]},
    {k, 0, Min[Length[pCoeffs] - 1, rel]}];
  boundary = Table[{top + h, Together[Sum[
      qCoeffs[[i + 1]]*coeffAt[rel + h - i], {i, h, qd}]]},
    {h, 1, qd}];
  reason = Which[
    !AllTrue[pCoeffs, groupedEpsExactSafeQ], "NumeratorByteCount",
    !AllTrue[csr, groupedEpsExactSafeQ], "ImpulseByteCount",
    !AllTrue[boundary[[All, 2]], groupedEpsExactSafeQ], "BoundaryByteCount",
    True, None];
  <|"Pairs" -> Select[Join[low, boundary], #[[2]] =!= 0 &],
    "Safe" -> (reason === None), "RouteReason" -> reason|>];

buildRationalMatrixGroup[records_List, eps_, top_Integer, d_Integer] := Module[
  {q = records[[1, 1]], qCoeffs, boundaryData, bad, decorated,
   shifts, mats, rules},
  qCoeffs = CoefficientList[q, eps];
  If[!groupedEpsExactSafeQ[q] || !AllTrue[qCoeffs, groupedEpsExactSafeQ],
    Return[<|"GroupedSafe" -> False, "Records" -> records,
      "RouteReason" -> "DenominatorByteCount"|>]];
  boundaryData = rationalBoundaryData[#, qCoeffs, eps, top] & /@ records;
  bad = SelectFirst[boundaryData, !TrueQ[#["Safe"]] &, None];
  If[bad =!= None,
    Return[<|"GroupedSafe" -> False, "Records" -> records,
      "RouteReason" -> bad["RouteReason"]|>]];
  decorated = MapThread[{#1[[2]], #1[[3]], #2["Pairs"]} &,
    {records, boundaryData}];
  shifts = Union @@ Map[#[[3, All, 1]] &, decorated];
  mats = Table[
    rules = Map[Function[rr, With[
        {z = SelectFirst[rr[[3]], #[[1]] === s &, None]},
        If[z === None, Nothing, {rr[[1]], rr[[2]]} -> z[[2]]]]],
      decorated];
    {s, SparseArray[rules, {d, d}]},
    {s, shifts}];
  <|"GroupedSafe" -> True, "Denominator" -> q,
    "DenominatorCoefficients" -> qCoeffs,
    "NumeratorSp" -> mats, "EntryCount" -> Length[records]|>];

(* Hybrid preparation: epsilon-Laurent polynomials stay on the established
   sparse-shift path.  Only entries with a genuinely nonconstant epsilon
   denominator are grouped, avoiding a regression on charts such as banana
   x=0 where Nhat is already finite. *)
prepareNhatHybrid[exprs_List, eps_, fb_Integer, W_Integer, cs_] := Module[
  {top = fb + W - 1, nj = Length[exprs], d = Length[First[exprs]],
   polySp = {}, ratGroups = {}, valuations, stats = {}, analyzed,
   polyL, idxs, records, builtGroups, groups, unsafeGroups,
   legacyRecords, legacyIndices},
  valuations = ConstantArray[Infinity, {nj, d, d}];
  Do[
    analyzed = Table[
      epsRationalData[exprs[[j, r, c]], eps, fb, top, cs],
      {r, d}, {c, d}];
    valuations[[j]] = Map[#["Valuation"] &, analyzed, {2}];
    records = Flatten[Table[
      If[TrueQ[analyzed[[r, c, "Zero"]]] ||
          analyzed[[r, c, "DenominatorDegree"]] === 0,
        Nothing,
        {analyzed[[r, c, "Q"]], r, c,
          analyzed[[r, c, "Valuation"]], analyzed[[r, c, "P"]]}],
      {r, d}, {c, d}], 1];
    builtGroups = If[records === {}, {},
      Map[buildRationalMatrixGroup[#, eps, top, d] &,
        GatherBy[records, First]]];
    groups = Select[builtGroups, TrueQ[#["GroupedSafe"]] &];
    unsafeGroups = Select[builtGroups, !TrueQ[#["GroupedSafe"]] &];
    legacyRecords = If[unsafeGroups === {}, {},
      Flatten[unsafeGroups[[All, "Records"]], 1]];
    legacyIndices = Association[Table[
      ((rec[[2]] - 1)*d + rec[[3]]) -> True, {rec, legacyRecords}]];
    polyL = Table[Which[
      TrueQ[analyzed[[r, c, "Zero"]]], ConstantArray[0, W],
      analyzed[[r, c, "DenominatorDegree"]] === 0,
        polynomialEpsListFromData[analyzed[[r, c]], eps, fb, W],
      KeyExistsQ[legacyIndices, (r - 1)*d + c],
        ratEpsList[exprs[[j, r, c]], eps, fb, W],
      True, ConstantArray[0, W]],
      {r, d}, {c, d}];
    idxs = Select[Range[W], Module[{i = #},
      AnyTrue[Flatten[polyL, 1], #[[i]] =!= 0 &]] &];
    AppendTo[polySp,
      Table[{i + fb - 1, Map[#[[i]] &, polyL, {2}]}, {i, idxs}]];
    AppendTo[ratGroups, groups];
    AppendTo[stats, <|
      "PolynomialShifts" -> Length[idxs],
      "LegacySparseShifts" -> Length[idxs],
      "RationalEntries" -> Length[records],
      "GroupedRationalEntries" -> Total[Lookup[groups, "EntryCount", 0]],
      "RationalGroups" -> Length[groups],
      "RationalNumeratorShifts" -> Total[
        Length[#["NumeratorSp"]] & /@ groups],
      "LegacyRationalEntries" -> Length[legacyRecords],
      "LegacyRationalGroups" -> Length[unsafeGroups],
      "LegacyRouteReasons" -> If[unsafeGroups === {}, <||>,
        Counts[unsafeGroups[[All, "RouteReason"]]]],
      "LegacyRationalShiftUpperBound" -> If[records === {}, 0,
        Length[Union @@ (Range[#[[4]], top] & /@ records)]],
      "FrameWidth" -> W|>],
    {j, nj}];
  <|"PolynomialSp" -> polySp, "RationalGroups" -> ratGroups,
    "Valuations" -> valuations, "Stats" -> stats|>];

(* B.U for one grouped rational matrix Mbar, where B = Q Mbar includes the
   finite-top boundary correction.  Keep this operation per (lag, group): a
   lower-frame witness must be raised before any cancellation with another
   lag is allowed. *)
rationalMatrixGroupNumerator[group_Association, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {rhs = ConstantArray[0, Dimensions[U]]},
  (* During an adaptive probe, apply the coefficient matrix before the
     lower-frame shift so banana's rank-one/square-zero annihilation is
     visible before the underflow decision.  Terminal/ordinary runs retain
     the established operation order bit-for-bit. *)
  Do[rhs += matrixShiftProduct[sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, group["NumeratorSp"]}];
  rhs];

(* Causal finite-frame division by the normalized epsilon unit Q.  Q(0) is
   nonzero, so this map is linear and lower triangular. *)
divideRationalMatrixRHS[rhs_List, q_List, W_Integer] := Module[
  {qd = Length[q] - 1, y = ConstantArray[0, Dimensions[rhs]]},
  y[[All, 1]] = rhs[[All, 1]]/q[[1]];
  Do[y[[All, i]] = (rhs[[All, i]] - Sum[
      q[[h + 1]]*y[[All, i - h]], {h, 1, Min[i - 1, qd]}])/q[[1]],
    {i, 2, W}];
  y];

(* Unfused single-lag form retained as the equivalence/benchmark seam. *)
applyRationalMatrixGroups[groups_List, U_List, fb_Integer, W_Integer, cs_] :=
 Module[{acc = ConstantArray[0, Dimensions[U]], rhs},
  Do[
    rhs = rationalMatrixGroupNumerator[group, U, fb, W, cs];
    acc += divideRationalMatrixRHS[
      rhs, group["DenominatorCoefficients"], W],
    {group, groups}];
  acc];

(* White-box benchmark/equivalence seam used by the focused tests. *)
applyPreparedNhat[polySp_List, groups_List, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {acc = ConstantArray[0, Dimensions[U]]},
  Do[acc += matrixShiftProduct[sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, polySp}];
  If[groups =!= {},
    acc += applyRationalMatrixGroups[groups, U, fb, W, cs]];
  acc];

(* One dense epsilon-rational matrix uses the same exact denominator groups
   and finite-top boundary corrections as Nhat.  The wrapper records the
   actual frame because adaptive homogeneous columns can finish on different
   lower rectangles. *)
preparedFramedMatrixValuations[polySp_List, groups_List, d_Integer] := Module[
  {vals = ConstantArray[Infinity, {d, d}], record},
  record[sp_List] := Do[
    If[sp[[2, r, c]] =!= 0,
      vals[[r, c]] = Min[vals[[r, c]], sp[[1]]]],
    {r, d}, {c, d}];
  Scan[record, polySp];
  Scan[Function[group, Scan[record, group["NumeratorSp"]]], groups];
  vals];

prepareFramedMatrix[M_?MatrixQ, eps_Symbol, fb_Integer, W_Integer, cs_] :=
 Module[{d = Length[M], p, polySp, groups},
  If[TrueQ[Normal[M] === IdentityMatrix[d]],
    Return[<|"PolynomialSp" -> {{0, IdentityMatrix[d]}},
      "RationalGroups" -> {}, "Valuations" ->
        Table[If[r === c, 0, Infinity], {r, d}, {c, d}],
      "FrameBase" -> fb, "FrameWidth" -> W, "Identity" -> True,
      "Stats" -> <|"PolynomialShifts" -> 1,
        "LegacySparseShifts" -> 1, "RationalEntries" -> 0,
        "GroupedRationalEntries" -> 0, "RationalGroups" -> 0,
        "RationalNumeratorShifts" -> 0, "LegacyRationalEntries" -> 0,
        "LegacyRationalGroups" -> 0, "LegacyRouteReasons" -> <||>,
        "LegacyRationalShiftUpperBound" -> 0, "FrameWidth" -> W|>|>]];
  p = prepareNhatHybrid[{M}, eps, fb, W, cs];
  polySp = First[p["PolynomialSp"]];
  groups = First[p["RationalGroups"]];
  <|"PolynomialSp" -> polySp, "RationalGroups" -> groups,
    (* Match the legacy spectral contract: valuation belongs to the actual
       prepared/Chopped finite frame, not the pre-numericization expression. *)
    "Valuations" -> preparedFramedMatrixValuations[polySp, groups, d],
    "FrameBase" -> fb, "FrameWidth" -> W, "Identity" -> False,
    "Stats" -> First[p["Stats"]]|>];

(* frConv checked every nonzero V_(r,c) product before summing a row.  Keep
   that strict witness even though the grouped application below performs a
   matrix product: exact cancellation between columns must not conceal a
   product whose leading epsilon power lies below the work frame. *)
framedMatrixLowerWitness[valuations_?MatrixQ, U_List,
    fb_Integer, cs_] := Module[{uVal, d = Length[valuations], v, lead},
  uVal = frameValuation[#, fb] & /@ U;
  Do[
    v = valuations[[r, c]];
    If[IntegerQ[v] && IntegerQ[uVal[[c]]],
      lead = v + uVal[[c]];
      If[lead < fb,
        lowerFrameUnderflow[cs, <|"FrameBase" -> fb,
          "Row" -> r, "Column" -> c, "MatrixValuation" -> v,
          "InputValuation" -> uVal[[c]], "ProductLead" -> lead,
          "Detail" ->
            "framed matrix product would discard nonzero lower-epsilon content"|>]]],
    {r, d}, {c, d}];
  uVal];

(* The entrywise witness has proved every active contribution in-frame.
   Multiply before shifting so content in a column unused by this sparse
   coefficient matrix cannot trigger a false lower-frame refusal.  The
   ordinary shift guard remains as a second exact assertion. *)
matrixShiftProductAfterWitness[A_, M_, s_Integer,
    fb_Integer, W_Integer, cs_] :=
  shiftFrameBlock[A . M, s, fb, W, cs];

rationalMatrixGroupNumeratorAfterWitness[group_Association, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {rhs = ConstantArray[0, Dimensions[U]]},
  Do[rhs += matrixShiftProductAfterWitness[
      sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, group["NumeratorSp"]}];
  rhs];

applyPreparedFramedMatrix[prep_Association, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {vals = prep["Valuations"], d = Length[U],
   acc = ConstantArray[0, Dimensions[U]], rhs},
  If[prep["FrameBase"] =!= fb || prep["FrameWidth"] =!= W ||
      Dimensions[U] =!= {d, W} || Dimensions[vals] =!= {d, d},
    err["E5", cs, <|"FrameBase" -> fb, "FrameWidth" -> W,
      "InputDimensions" -> Dimensions[U],
      "ValuationDimensions" -> Dimensions[vals],
      "Detail" -> "prepared framed matrix metadata/dimensions mismatch"|>]];
  If[TrueQ[Lookup[prep, "Identity", False]], Return[U]];
  framedMatrixLowerWitness[vals, U, fb, cs];
  Do[acc += matrixShiftProductAfterWitness[
      sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, prep["PolynomialSp"]}];
  Do[
    rhs = rationalMatrixGroupNumeratorAfterWitness[
      group, U, fb, W, cs];
    acc += divideRationalMatrixRHS[
      rhs, group["DenominatorCoefficients"], W],
    {group, prep["RationalGroups"]}];
  acc];

(* framed product: a, b based at fb -> product re-framed at fb.  A nonzero
   product below the physical frame is a budgeting failure, never a
   sanctioned truncation.  Exact zero and one-term frames occur pervasively
   in the spectral V/VInv and gauge transforms.  Their product is exactly a
   zero or an epsilon shift-and-scale, so do not pay for a dense W-by-W
   ListConvolve at every Taylor row/basis column.  Support classification is
   structural (=!= 0) only; genuinely multi-term frames retain the original
   convolution path unchanged. *)
frConv[a_, b_, fb_, W_] := Module[{sa, sb, ia, ib, lead, shiftMono},
  sa = Select[Range[W], a[[#]] =!= 0 &];
  sb = Select[Range[W], b[[#]] =!= 0 &];
  If[sa === {} || sb === {}, Return[ConstantArray[0, W], Module]];
  ia = First[sa]; ib = First[sb];
  lead = 2 fb + ia + ib - 2;
  If[lead < fb,
    lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
      <|"FrameBase" -> fb, "ProductLead" -> lead,
        "Detail" -> "framed convolution would discard nonzero lower-epsilon content"|>]];
  shiftMono[mono_, pos_, other_] := Module[
    {out = ConstantArray[0, W], shift = fb + pos - 1, dstLo, dstHi},
    dstLo = Max[1, 1 + shift];
    dstHi = Min[W, W + shift];
    If[dstLo <= dstHi,
      out[[dstLo ;; dstHi]] =
        mono*other[[dstLo - shift ;; dstHi - shift]]];
    out];
  If[Length[sa] === 1, Return[shiftMono[a[[ia]], ia, b], Module]];
  If[Length[sb] === 1, Return[shiftMono[b[[ib]], ib, a], Module]];
  Take[ListConvolve[a, b, {1, -1}, 0], {1 - fb, W - fb}]];

(* framed inverse: leading nonzero at frame position pos0 (eps-power
   fb+pos0-1); the inverse leads at the negated power (must lie in-frame) *)
frInv[dlist_, fb_, W_] := Module[{pos0, drel, m, e, out, lead, startIdx},
  pos0 = SelectFirst[Range[W], dlist[[#]] =!= 0 &, None];
  If[pos0 === None,
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "framed inverse of an identically zero series"|>]];
  lead = fb + pos0 - 1;
  startIdx = -lead - fb + 1;   (* frame index of eps^(-lead) *)
  If[startIdx < 1,
    err["E4", <|"Center" -> "(frame)"|>,
      <|"Lead" -> lead, "FrameBase" -> fb,
        "Detail" -> "inverse leads below the work frame"|>]];
  If[startIdx > W,
    err["E4", <|"Center" -> "(frame)"|>,
      <|"Lead" -> lead, "FrameTop" -> fb + W - 1,
        "Detail" -> "inverse leads above the work frame"|>]];
  drel = dlist[[pos0 ;;]];
  m = W - startIdx;            (* relative orders we can store *)
  e = ConstantArray[0, m + 1];
  e[[1]] = If[NumericQ[drel[[1]]], 1/drel[[1]], Together[1/drel[[1]]]];
  Do[e[[k + 1]] = -e[[1]]*Sum[
      If[j + 1 <= Length[drel], drel[[j + 1]], 0]*e[[k - j + 1]], {j, 1, k}],
    {k, 1, m}];
  out = ConstantArray[0, W];
  out[[startIdx ;; startIdx + m]] = e;
  out];

(* x eps and / eps as frame shifts *)
frShiftUp[a_] := Prepend[Most[a], 0];
frDivEps[a_] := Module[{},
  If[First[a] =!= 0,
    lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
      <|"Detail" -> "division by epsilon would discard the lowest framed coefficient"|>]];
  Append[Rest[a], 0]];

epsValuation[e_, eps_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    Exponent[Numerator[c], eps, Min] - Exponent[Denominator[c], eps, Min]]];

(* Symbolic finite-width data.  D is defined only up to multiplication by
   a nonzero function of eps.  Removing its polynomial content in t
   before the d0(eps=0) check: a global 1/eps coefficient is a legitimate
   Laurent multiplier, while a moving factor such as (t-eps) still fails
   E3 because its primitive d0 vanishes at eps=0. *)
$clearedSymbolicLegacyCache = <||>;
$clearedSymbolicLegacyCacheMax = 32;

clearedSymbolicLegacyCore[cs_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], den,
   denCoeffs, denContent, num, d0, dD, dN, dExpr, NhatExpr,
   phaseQ, phaseTime, phase, identityVQ, rightPolynomial,
   transformedPolynomial, coefficientLists},
  phaseQ = Environment["DEBUG_SOLVE_PHASES"] === "1";
  phaseTime = SessionTime[];
  phase[label_String] := If[phaseQ, Module[{now = SessionTime[]},
    Print["SOLVEPHASE center=", cs["Center"], " phase=cleared-", label,
      " dt=", N[now - phaseTime, 6]];
    phaseTime = now]];
  den = Together[PolynomialLCM @@ (Denominator[Together[#]] & /@ Flatten[cs["ThetaMatrix"]])];
  phase["denominator-lcm"];
  denCoeffs = Select[CoefficientList[den, t], !zeroCanQ[#] &];
  If[denCoeffs === {},
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator is identically zero"|>]];
  denContent = If[Length[denCoeffs] === 1, First[denCoeffs],
    Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]];
  den = Cancel[Together[den/denContent]];
  phase["denominator-content"];
  num = Map[Cancel[Together[#*den]] &, cs["ThetaMatrix"], {2}];
  phase["numerators"];
  d0 = Together[den /. t -> 0];
  If[zeroQ[d0 /. eps -> 0],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  dExpr = Table[Cancel[Together[Coefficient[den, t, j]]], {j, 0, dD}];
  phase["denominator-coefficients"];
  identityVQ = !TrueQ[$disableIdentityNhatShortcut] &&
    TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    TrueQ[Normal[cs["V"]] === IdentityMatrix[cs["SystemSize"]]] &&
    TrueQ[Normal[cs["VInv"]] === IdentityMatrix[cs["SystemSize"]]];
  NhatExpr = Which[
    identityVQ,
      Table[Map[Cancel[Together[#]] &,
        Map[Coefficient[#, t, j] &, num, {2}], {2}], {j, 0, dN}],
    TrueQ[$disablePolynomialNhatTransform],
      Table[Module[{Nj},
        Nj = Map[Coefficient[#, t, j] &, num, {2}];
        Map[Cancel[Together[#]] &, cs["VInv"] . Nj . cs["V"], {2}]],
        {j, 0, dN}],
    True,
      (* V and VInv are residue-frame matrices and therefore t-independent.
         Linearity of exact coefficient extraction then permits ONE
         polynomial similarity transform instead of dN+1 scalar transforms:
           [t^j](VInv.num.V) = VInv.[t^j]num.V.
         Keep the multiplication explicitly right-first; the final
         coefficient canonicalization remains identical to the legacy
         contract used by epsilon valuation and preparation. *)
      If[!FreeQ[{cs["V"], cs["VInv"]}, t],
        err["E5", cs, <|
          "Detail" -> "spectral frame depends on the chart variable during polynomial Nhat transform"|>]];
      If[!AllTrue[Flatten[num], PolynomialQ[#, t] &],
        err["E5", cs, <|
          "Detail" -> "cleared numerator is not polynomial in the chart variable"|>]];
      rightPolynomial = num . cs["V"];
      phase["nhat-polynomial-right"];
      transformedPolynomial = cs["VInv"] . rightPolynomial;
      phase["nhat-polynomial-left"];
      transformedPolynomial = Map[Cancel[Together[#]] &,
        transformedPolynomial, {2}];
      phase["nhat-polynomial-canonicalize"];
      If[!AllTrue[Flatten[transformedPolynomial], PolynomialQ[#, t] &],
        err["E5", cs, <|
          "Detail" -> "spectrally transformed numerator is not polynomial in the chart variable"|>]];
      coefficientLists = Map[CoefficientList[#, t] &,
        transformedPolynomial, {2}];
      phase["nhat-coefficient-lists"];
      Table[Map[Cancel[Together[#]] &,
        Map[If[j + 1 <= Length[#], #[[j + 1]], 0] &,
          coefficientLists, {2}], {2}], {j, 0, dN}]
    ];
  phase[If[identityVQ, "nhat-identity-final",
    If[TrueQ[$disablePolynomialNhatTransform], "nhat-legacy-final",
      "nhat-final-canonicalize"]]];
  <|"dExpr" -> dExpr, "NhatExpr" -> NhatExpr,
    "dD" -> dD, "dN" -> dN|>];

(* Singular SCC planning asks for the same exact spectral and physical
   operators from several independent proof layers: seed-halo costing,
   homogeneous frame planning, and request capture.  The operator is a pure
   function of the exact reduced theta matrix and spectral frame, so retain a
   small content-checked cache instead of repeating algebraic coefficient
   extraction and RootReduce work at every layer. *)
clearedSymbolicLegacy[cs_Association] := Module[
  {signature, key, cached, data},
  signature = {"cleared-symbolic-legacy-v1",
    cs["ChartVar"], cs["SystemSize"], cs["ThetaMatrix"],
    cs["V"], cs["VInv"],
    TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    TrueQ[$disableIdentityNhatShortcut],
    TrueQ[$disablePolynomialNhatTransform]};
  key = Hash[signature, "SHA256"];
  cached = Lookup[$clearedSymbolicLegacyCache, key, None];
  If[AssociationQ[cached] &&
      SameQ[Lookup[cached, "Signature", None], signature],
    Return[cached["Data"], Module]];
  If[cached =!= None,
    err["E6", cs, <|"CacheKey" -> key,
      "Detail" ->
        "legacy cleared-symbolic cache key collided with unequal exact input"|>]];
  data = clearedSymbolicLegacyCore[cs];
  If[Length[$clearedSymbolicLegacyCache] >=
      $clearedSymbolicLegacyCacheMax,
    KeyDropFrom[$clearedSymbolicLegacyCache,
      First[Keys[$clearedSymbolicLegacyCache]]]];
  AssociateTo[$clearedSymbolicLegacyCache, key ->
    <|"Signature" -> signature, "Data" -> data|>];
  data];

$disableGlobalClearedHoist = False;
$disableIdentityNhatShortcut = False;
$disablePolynomialNhatTransform = False;

regularIdentityFrameQ[cs_Association] :=
  TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
  TrueQ[cs["Gauge"] === IdentityMatrix[cs["SystemSize"]]] &&
  TrueQ[Normal[cs["V"]] === IdentityMatrix[cs["SystemSize"]]] &&
  TrueQ[Normal[cs["VInv"]] === IdentityMatrix[cs["SystemSize"]]];

(* Exact clearing in the original system variable.  Affine substitution is a
   field automorphism, so the LCM of canceled entry denominators may be shifted
   to every regular chart instead of recomputed after the algebraic shift. *)
globalClearedSystem[systemKey_] := Module[
  {cached, input, x, A, den, denCoeffs, denContent, num,
   numCoefficientLists, phaseQ, phaseTime},
  cached = If[KeyExistsQ[$globalClearedCache, systemKey],
    $globalClearedCache[systemKey], None];
  If[cached =!= None, Return[cached, Module]];
  input = If[KeyExistsQ[$systemClearRegistry, systemKey],
    $systemClearRegistry[systemKey], None];
  If[input === None,
    err["E6", <|"Center" -> "global-system-clear"|>,
      <|"Key" -> systemKey,
        "Detail" -> "chart references an unregistered exact system"|>]];
  x = input["Variable"];
  phaseQ = Environment["DEBUG_SOLVE_PHASES"] === "1";
  phaseTime = SessionTime[];
  A = Map[Cancel[Together[#]] &, input["Matrix"], {2}];
  den = Together[PolynomialLCM @@
    (Denominator[Together[#]] & /@ Flatten[A])];
  If[phaseQ, Print["SOLVEPHASE phase=global-denominator-lcm dt=",
    N[SessionTime[] - phaseTime, 6]]];
  denCoeffs = Select[CoefficientList[den, x], !zeroCanQ[#] &];
  If[denCoeffs === {},
    err["E3", <|"Center" -> "global-system-clear"|>,
      <|"Denominator" -> den,
        "Detail" -> "global cleared denominator is identically zero"|>]];
  denContent = If[Length[denCoeffs] === 1, First[denCoeffs],
    Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]];
  den = Cancel[Together[den/denContent]];
  num = Map[Cancel[Together[#*den]] &, A, {2}];
  If[!PolynomialQ[den, x] ||
      !AllTrue[Flatten[num], PolynomialQ[#, x] &],
    err["E5", <|"Center" -> "global-system-clear"|>, <|
      "Detail" ->
        "globally cleared denominator or numerator is not polynomial in the system variable"|>]];
  denCoeffs = CoefficientList[den, x];
  numCoefficientLists = Map[CoefficientList[#, x] &, num, {2}];
  cached = <|"Variable" -> x, "Denominator" -> den,
    "Numerator" -> num,
    "DenominatorCoefficientList" -> denCoeffs,
    "NumeratorCoefficientLists" -> numCoefficientLists,
    "Dimension" -> Length[A]|>;
  AssociateTo[$globalClearedCache, systemKey -> cached];
  cached];

trimExactPolynomialCoefficientList[list_List] := Module[{last},
  last = SelectFirst[Reverse[Range[Length[list]]],
    !zeroCanQ[list[[#]]] &, None];
  If[last === None, {0}, Take[list, last]]];

affinePolynomialCoefficientTransform[degree_Integer?NonNegative,
    center_, beta_] := Table[
  If[m < j, 0,
    If[j === 0, 1, beta^j] Binomial[m, j]
      If[m === j, 1, center^(m - j)]],
  {j, 0, degree}, {m, 0, degree}];

affineTranslatePolynomialCoefficientList[list_List, transform_List] :=
  trimExactPolynomialCoefficientList[
    Take[transform, {1, Length[list]}, {1, Length[list]}] . list];

polynomialCoefficientListValuation[list_List] := Module[{first},
  first = SelectFirst[Range[Length[list]], !zeroCanQ[list[[#]]] &,
    Infinity];
  first - 1];

affinePhysicalClearedFromGlobal[cs_Association] := Module[
  {systemKey = cs["SystemClearKey"], key, cached, global, x,
   t = cs["ChartVar"], center = cs["Center"], beta = cs["ChartMap", "Scale"],
   affineContentInvariantQ, den, denCoeffs, denContent, num, dD, dN,
   dExpr, NhatExpr, coefficientLists, activePairEntries, centerPower,
   transformDegree, coefficientTransform, numCoefficientLists,
   activeCoefficientLists,
   phaseQ, phaseTime, phase},
  key = {systemKey, center, beta, t};
  cached = If[KeyExistsQ[$chartClearedCache, key],
    $chartClearedCache[key], None];
  If[cached =!= None, Return[cached, Module]];
  global = globalClearedSystem[systemKey];
  x = global["Variable"];
  phaseQ = Environment["DEBUG_SOLVE_PHASES"] === "1";
  phaseTime = SessionTime[];
  phase[label_String] := If[phaseQ, Module[{now = SessionTime[]},
    Print["SOLVEPHASE center=", center, " phase=hoist-", label,
      " dt=", N[now - phaseTime, 6]];
    phaseTime = now]];
  (* globalClearedSystem has already removed the content of D(x,eps) in x.
     For coefficient-independent exact numeric c and nonzero beta,
     x -> c + beta t is an invertible coefficient-field automorphism.  It
     therefore preserves primitivity up to a unit: a nonunit epsilon/parameter
     factor common to every t coefficient after the substitution would pull
     back to the same common factor of every x coefficient before it.  Do not
     re-prove that invariant with a PolynomialGCD at every ordinary chart.
     Apart from being redundant, the proof can take minutes once an atlas has
     accumulated large rational centers.  A chart-constant unit introduced by
     the affine map need not be removed because multiplying both sides of
     D theta = N by that unit leaves the recurrence and physical differential
     equation unchanged.

     Symbolic geometry may itself contain eps or an analytic regulator and
     can introduce genuine nonunit content (for example beta=eps).  Keep the
     old exact GCD fallback for that out-of-contract case. *)
  If[zeroCanQ[beta],
    err["E3", cs, <|
      "Detail" -> "affine-clearing scale must be nonzero"|>]];
  affineContentInvariantQ = AllTrue[{center, beta},
    NumericQ[#] && FreeQ[#, _?InexactNumberQ] &];
  (* For ordinary exact geometry, translate the already-cleared global
     coefficient lists directly.  If p(x)=Sum[p_m x^m,m], then

       [t^j] p(center+beta t) =
         beta^j Sum[Binomial[m,j] center^(m-j) p_m,{m,j,degree}].

     This is the same field automorphism used by the expression path below,
     but it performs the polynomial expansion once globally and only a
     triangular linear transform at each chart.  In particular it avoids
     rebuilding and re-expanding all 23^2 large rational expressions at each
     Henn level-3 center. *)
  If[affineContentInvariantQ,
    transformDegree = Max[
      Length[global["DenominatorCoefficientList"]] - 1,
      Max[Length /@ Flatten[
        global["NumeratorCoefficientLists"], 1]] - 1];
    coefficientTransform = affinePolynomialCoefficientTransform[
      transformDegree, center, beta];
    dExpr = affineTranslatePolynomialCoefficientList[
      global["DenominatorCoefficientList"], coefficientTransform];
    numCoefficientLists = Map[
      affineTranslatePolynomialCoefficientList[#, coefficientTransform] &,
      global["NumeratorCoefficientLists"], {2}];
    numCoefficientLists = Map[Prepend[beta*#, 0] &,
      numCoefficientLists, {2}];
    activeCoefficientLists = Prepend[
      Select[Flatten[numCoefficientLists, 1],
        AnyTrue[#, !zeroCanQ[#] &] &], dExpr];
    centerPower = Min[
      polynomialCoefficientListValuation /@ activeCoefficientLists];
    If[!IntegerQ[centerPower] || centerPower < 0,
      err["E3", cs, <|"CenterPower" -> centerPower,
        "Detail" ->
          "affine physical q/C pair has an invalid common center valuation"|>]];
    If[centerPower > 0,
      dExpr = Drop[dExpr, centerPower];
      numCoefficientLists = Map[
        If[AllTrue[#, zeroCanQ], {0}, Drop[#, centerPower]] &,
        numCoefficientLists, {2}]];
    dExpr = trimExactPolynomialCoefficientList[dExpr];
    numCoefficientLists = Map[trimExactPolynomialCoefficientList,
      numCoefficientLists, {2}];
    dD = Length[dExpr] - 1;
    dN = Max[0,
      Max[Length /@ Flatten[numCoefficientLists, 1]] - 1];
    NhatExpr = Table[
      Map[If[j + 1 <= Length[#], #[[j + 1]], 0] &,
        numCoefficientLists, {2}],
      {j, 0, dN}];
    phase["coefficient-translation"];
    cached = <|"dExpr" -> dExpr, "NhatExpr" -> NhatExpr,
      "dD" -> dD, "dN" -> dN|>;
    AssociateTo[$chartClearedCache, key -> cached];
    Return[cached, Module]];
  den = Cancel[Together[global["Denominator"] /. x -> center + beta*t]];
  If[zeroCanQ[den],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "affine-shifted global denominator is identically zero"|>]];
  denContent = If[affineContentInvariantQ, 1,
    denCoeffs = Select[CoefficientList[den, t], !zeroCanQ[#] &];
    If[denCoeffs === {},
      err["E3", cs, <|"Denominator" -> den,
        "Detail" -> "affine-shifted global denominator is identically zero"|>]];
    If[Length[denCoeffs] === 1, First[denCoeffs],
      Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]]];
  den = Cancel[Together[den/denContent]];
  num = Map[Cancel[Together[beta*t*(# /. x -> center + beta*t)/denContent]] &,
    global["Numerator"], {2}];
  (* D f' = N f becomes D theta f = beta t N f.  At a singular chart
     center, the affine denominator D can contain t, which is then a common
     factor of both sides solely because of the explicit theta multiplier.
     The legacy local construction starts from theta=beta t A and cancels
     this factor automatically.  The hoisted path must do the same: retaining
     it makes q(0)=C(0)=0 and fabricates a singular/resonant Taylor layer.

     Global primitivity and affine invertibility prove that t is the only new
     possible common polynomial factor, so compare exact center valuations;
     do not reintroduce the expensive full polynomial GCD that this path was
     designed to avoid. *)
  activePairEntries = Prepend[
    Select[Flatten[num], !zeroCanQ[#] &], den];
  centerPower = Min[Exponent[#, t, Min] & /@ activePairEntries];
  If[!IntegerQ[centerPower] || centerPower < 0,
    err["E3", cs, <|"CenterPower" -> centerPower,
      "Detail" -> "affine physical q/C pair has an invalid common center valuation"|>]];
  If[centerPower > 0,
    den = Cancel[Together[den/t^centerPower]];
    num = Map[Cancel[Together[#/t^centerPower]] &, num, {2}]];
  phase["affine-pair"];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  dExpr = Table[Cancel[Together[Coefficient[den, t, j]]], {j, 0, dD}];
  (* num is already an entrywise-canonical polynomial matrix from the exact
     global clearing and affine substitution above.  Asking Coefficient for
     every {Taylor lag,row,column} independently makes Wolfram re-expand the
     same polynomial and run thousands of unrelated rational GCDs.  This was
     the dominant Henn level-3 anchor cost (hours for a 23-by-23 system).
     Extract every entry once and index its exact coefficient list.  Unlike
     the spectral legacy path there is no subsequent matrix product here, so
     the coefficients need no second Together/Cancel canonicalization. *)
  If[!AllTrue[Flatten[num], PolynomialQ[#, t] &],
    err["E5", cs, <|
      "Detail" ->
        "affine globally cleared numerator is not polynomial in the chart variable"|>]];
  coefficientLists = Map[CoefficientList[#, t] &, num, {2}];
  phase["coefficient-lists"];
  NhatExpr = Table[
    Map[If[j + 1 <= Length[#], #[[j + 1]], 0] &,
      coefficientLists, {2}],
    {j, 0, dN}];
  phase["coefficients"];
  cached = <|"dExpr" -> dExpr, "NhatExpr" -> NhatExpr,
    "dD" -> dD, "dN" -> dN|>;
  AssociateTo[$chartClearedCache, key -> cached];
  cached];

regularClearedFromGlobal[cs_Association] := Module[
  {data = affinePhysicalClearedFromGlobal[cs], d0},
  d0 = First[data["dExpr"]];
  If[zeroQ[d0 /. DiffExp2`Config`CanonicalEps[] -> 0],
    err["E3", cs, <|"Denominator" -> d0,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  data];

clearedSymbolic[cs_Association] := If[
  !TrueQ[$disableGlobalClearedHoist] && regularIdentityFrameQ[cs] &&
    KeyExistsQ[cs, "SystemClearKey"],
  regularClearedFromGlobal[cs],
  clearedSymbolicLegacy[cs]];

(* The ordinary Frobenius recurrence diagonalizes the residue first and
   therefore prepares V^-1 N(t) V.  When V degenerates at eps=0, that
   representation can turn an eps-integral differential equation into a
   recurrence with a repeatable artificial eps pole.  The confluent path
   keeps the same exact t-Fuchsian gauge but prepares the reduced physical
   basis instead.  V is still used once for homogeneous seeds and V^-1 once
   for a sourced n=0 resonance; neither is allowed to contaminate every
   later Taylor lag. *)
clearedPhysicalSymbolic[cs_Association] := Module[
  {d = cs["SystemSize"], identityGaugeQ},
  identityGaugeQ =
    TrueQ[Normal[cs["Gauge"]] === IdentityMatrix[d]] &&
    TrueQ[cs["ThetaMatrix"] === cs["ThetaOriginal"]];
  If[!TrueQ[$disableGlobalClearedHoist] && identityGaugeQ &&
      KeyExistsQ[cs, "SystemClearKey"],
    (* A singular SCC block is prepared as its own registered exact system.
       When rank reduction is the identity, its reduced physical recurrence
       is therefore precisely the affine image of that registered system.
       Reuse the global q/C proof instead of rebuilding the same 12-by-12
       PolynomialLCM/GCD locally merely because the chart is singular. *)
    affinePhysicalClearedFromGlobal[cs],
    clearedSymbolicLegacy[Join[cs, <|
      "V" -> IdentityMatrix[d], "VInv" -> IdentityMatrix[d]|>]]]];

epsilonRegularPrincipalCertificate[cs_Association,
    nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], blocks, physical, dVals, nVals,
   d0Val, integralOperators, integralSeeds, positiveTaylor,
   originalPoleDepth, eligible, reason},
  blocks = blockList[cs];
  physical = clearedPhysicalSymbolic[cs];
  dVals = Select[epsValuation[#, eps] & /@ physical["dExpr"],
    IntegerQ];
  nVals = Select[Map[epsValuation[#, eps] &,
      Flatten[physical["NhatExpr"]]], IntegerQ];
  d0Val = epsValuation[First[physical["dExpr"]], eps];
  integralOperators = IntegerQ[d0Val] && d0Val === 0 &&
    AllTrue[Join[dVals, nVals], # >= 0 &];
  integralSeeds = matrixEpsPoleDepth[cs["V"]] === 0;
  originalPoleDepth = recurrencePoleDepth[clearedSymbolic[cs], nmax];
  (* Positive Taylor layers are normally solved in the epsilon-integral
     reduced physical frame.  A T/P/R layer is exceptional only at that one
     n: RecurrenceSolver routes the complete affine-Jordan schedule through
     its retained exact spectral principal/source transaction, then returns
     to the physical recurrence.  That transaction supports arbitrary
     certified Jordan chains and propagates CASE-P hits to the composite
     compensator, so neither a nontrivial chain nor a positive resonance
     justifies applying V^-1 N_j V (and its artificial epsilon pole) at every
     later Taylor layer. *)
  positiveTaylor = AllTrue[Flatten[Table[
      !zeroCanQ[Cancel[Together[target["a"] + n - other["a"]]]],
      {target, blocks}, {n, 1, nmax}, {other, blocks}]], TrueQ];
  eligible = !TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    nmax >= 1 && originalPoleDepth > 0 && integralOperators &&
    integralSeeds;
  reason = Which[
    TrueQ[Lookup[cs["IndicialData"], "Regular", False]], "regular-chart",
    nmax < 1, "no-positive-taylor-layer",
    originalPoleDepth === 0, "spectral-recurrence-already-epsilon-integral",
    !integralOperators, "nonintegral-reduced-operator",
    !integralSeeds, "nonintegral-homogeneous-seed-frame",
    True, "certified-affine-jordan-epsilon-regular-principal-recurrence"];
  <|"Eligible" -> eligible, "Reason" -> reason,
    "Symbolic" -> physical, "D0Valuation" -> d0Val,
    "MinimumOperatorValuation" -> If[Join[dVals, nVals] === {},
      Infinity, Min[Join[dVals, nVals]]],
    "OriginalRecurrencePoleDepth" -> originalPoleDepth,
    "SpectralSeedPoleDepth" -> matrixEpsPoleDepth[cs["V"]],
    "InitialSourcePoleDepth" -> inverseSpectralTransformPoleDepth[cs],
    "PositiveTaylorNonresonant" -> positiveTaylor,
    "ExceptionalLayerPolicy" ->
      "exact-affine-jordan-spectral-transaction"|>];

(* A negative-valued j-th Taylor multiplier can be used again after every
   j Taylor steps.  The old one-hit epsPoleDepth misses this composition.
   This longest-path bound is deterministic and tighter than nmax times the
   worst pole when the pole occurs only at j>1. *)
recurrencePoleDepth[data_Association, nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], dD = data["dD"],
   dN = data["dN"], maxJ, stepCost, depth, phaseQ, phaseTime},
  phaseQ = Environment["DEBUG_SOLVE_PHASES"] === "1";
  phaseTime = SessionTime[];
  If[nmax <= 0, Return[0]];
  maxJ = Min[nmax, Max[dD, dN]];
  If[maxJ <= 0, Return[0]];
  stepCost = Table[Module[{vals = {}},
      If[j <= dD, AppendTo[vals, epsValuation[data["dExpr"][[j + 1]], eps]]];
      If[j <= dN, vals = Join[vals,
        Flatten[Map[epsValuation[#, eps] &, data["NhatExpr"][[j + 1]], {2}]]]];
      vals = Select[vals, IntegerQ];
      If[vals === {}, 0, Max[0, -Min[vals]]]],
    {j, 1, maxJ}];
  If[phaseQ, Print["SOLVEPHASE phase=pole-depth-valuations dt=",
    N[SessionTime[] - phaseTime, 6], " maxLag=", maxJ]];
  phaseTime = SessionTime[];
  depth = ConstantArray[0, nmax + 1];
  Do[depth[[n + 1]] = Max[Join[{0}, Table[
      depth[[n - j + 1]] + stepCost[[j]], {j, 1, Min[n, maxJ]}]]],
    {n, 1, nmax}];
  If[phaseQ, Print["SOLVEPHASE phase=pole-depth-dp dt=",
    N[SessionTime[] - phaseTime, 6], " result=", depth[[-1]]]];
  depth[[-1]]];

(* Deepest pole of ONE prepared Taylor multiplier, without composing it
   along the recurrence.  prepareCleared materializes every dExpr/NhatExpr
   lag before runRecursion starts, so the adaptive first rectangle must
   contain even a deep high-j operator (for example t^5/eps^10 when j=1 is
   pole-free).  The longer recurrencePoleDepth remains the terminal bound. *)
recurrenceSingleUsePoleDepth[data_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], vals},
  vals = Join[
    Map[epsValuation[#, eps] &, data["dExpr"]],
    Flatten[Map[epsValuation[#, eps] &, data["NhatExpr"], {3}]]];
  vals = Select[vals, IntegerQ];
  If[vals === {}, 0, Max[0, -Min[vals]]]];

(* Exact epsilon valuation of the t-Laurent coefficients that applyGauge
   will use.  Evaluating T at a generic numeric t is not sufficient: a
   mixed denominator such as eps+t can have finite generic valuation while
   its expansion at the chart origin contains increasingly deep eps poles.
   Keeping the coefficient expressions also lets applyGauge and the work
   budget use precisely the same structural valuation. *)
gaugeCoefficientData[cs_, nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"],
   T = cs["Gauge"], gv, coeffs, vals},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[<|"TMin" -> 0, "EpsValuation" -> 0, "Coefficients" -> None|>]];
  gv = Min[tVal[#, t] & /@ Flatten[T]];
  (* applyGauge convolves output orders np=0..nmax with
     m=gv..gv+nmax.  Budget and expand exactly that same coefficient range,
     including the allowed (though uncommon) gv>0 case. *)
  coeffs = Table[Map[tLaurent[#, t, m] &, T, {2}],
    {m, gv, gv + nmax}];
  vals = Select[Map[epsValuation[#, eps] &, Flatten[coeffs]], IntegerQ];
  <|"TMin" -> gv,
    "EpsValuation" -> If[vals === {}, 0, Min[0, Min[vals]]],
    "Coefficients" -> coeffs|>];

(* The native planner needs only a safe lower bound on the epsilon
   valuation of the finite gauge-coefficient range; it does not consume the
   coefficients themselves.  When the denominator has only a monomial
   epsilon content, division in t is by an epsilon-independent power series.
   Every requested quotient coefficient is then an epsilon-independent
   linear combination of numerator coefficients of no higher t order.  The
   minimum numerator valuation is therefore an exact, conservative bound,
   obtained without materializing nmax+1 Together[] expressions.  Mixed
   denominators such as eps+t retain the full coefficient calculation: their
   t expansion can create progressively deeper epsilon poles. *)
gaugeCoefficientEpsValuation[cs_, nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"],
   T = cs["Gauge"], entryBound, bounds},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[0, Module]];
  entryBound[entry_] := Module[
    {c = Cancel[Together[entry]], num, den, denEpsMin, reducedDen,
     numTMin, reducedNum, degree, coeffs, vals},
    If[zeroCanQ[c], Return[Infinity, Module]];
    num = Numerator[c]; den = Denominator[c];
    If[!PolynomialQ[num, {t, eps}] || !PolynomialQ[den, {t, eps}],
      Return[$Failed, Module]];
    denEpsMin = Exponent[den, eps, Min];
    reducedDen = Cancel[den/eps^denEpsMin];
    If[!FreeQ[reducedDen, eps], Return[$Failed, Module]];
    numTMin = polyMinDeg[num, t];
    reducedNum = Cancel[num/t^numTMin];
    degree = Min[nmax, Exponent[reducedNum, t]];
    coeffs = Table[Coefficient[reducedNum, t, j], {j, 0, degree}];
    vals = Select[Map[epsValuation[#, eps] &, coeffs], IntegerQ];
    If[vals === {}, Infinity, Min[vals] - denEpsMin]];
  bounds = entryBound /@ Flatten[T];
  If[MemberQ[bounds, $Failed],
    gaugeCoefficientData[cs, nmax]["EpsValuation"],
    bounds = Select[bounds, IntegerQ];
    If[bounds === {}, 0, Min[0, Min[bounds]]]]];

matrixEpsPoleDepth[m_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], vVals, vMin},
  vVals = Select[Map[epsValuation[#, eps] &, Flatten[m]], IntegerQ];
  vMin = If[vVals === {}, 0, Min[0, Min[vVals]]];
  -vMin];

spectralTransformPoleDepth[cs_] := matrixEpsPoleDepth[cs["V"]];
inverseSpectralTransformPoleDepth[cs_] := matrixEpsPoleDepth[cs["VInv"]];

(* assembleSolution applies V and then T as two finite Cauchy products.
   Their honest-window losses therefore add, even if entries of the formal
   product T.V would cancel.  This is private work headroom only: the
   delivered window is still capped to the caller's request. *)
finalTransformPoleDepth[cs_, nmax_Integer] :=
  spectralTransformPoleDepth[cs] -
    gaugeCoefficientEpsValuation[cs, nmax];

prepareCleared[cs_, fb_, W_, symbolic_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], data, dD, dN, dL, nhatPrep,
   d0Expr, d0InvScalar, ratGroups, ratDenominators, ratDenominatorKeys,
   numericGeometryQ},
  data = If[symbolic === Automatic, clearedSymbolic[cs], symbolic];
  dD = data["dD"]; dN = data["dN"];
  (* Exact algebraic centers/scales are valuable to the planner, but keeping
     their transformed scalar operators exact through dozens of recurrence
     steps causes severe rational/algebraic expression swell.  On a regular
     chart all structural decisions are already complete, so ground every
     exact numeric operator coefficient once at 2x WP.  Symbolic analytic
     regulators are not NumberQ and remain exact. *)
  numericGeometryQ = TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    AnyTrue[{cs["Center"], cs["ChartMap", "Scale"]},
      exactNonRationalNumberQ];
  Block[{$numericizeAllPreparedNumbers = numericGeometryQ},
    dL = Map[ratEpsList[#, eps, fb, W] &, data["dExpr"]];
    nhatPrep = prepareNhatHybrid[data["NhatExpr"], eps, fb, W, cs]];
  (* Cross-lag fusion keys.  epsRationalData already normalized every Q to
     Q(0)=1, so SameQ of the exact denominator expression is the canonical
     equivalence relation.  Store the small integer index once rather than
     rediscovering it in every Taylor recurrence. *)
  ratGroups = nhatPrep["RationalGroups"];
  ratDenominators = DeleteDuplicatesBy[
    Flatten[ratGroups, 1], #["Denominator"] &];
  ratDenominators = KeyTake[#,
      {"Denominator", "DenominatorCoefficients"}] & /@ ratDenominators;
  ratDenominatorKeys = Lookup[ratDenominators, "Denominator", {}];
  ratGroups = Map[Function[groups, Map[Function[group,
      Append[group, "DenominatorIndex" -> SelectFirst[
        Range[Length[ratDenominatorKeys]],
        ratDenominatorKeys[[#]] === group["Denominator"] &]]], groups]],
    ratGroups];
  (* Sparse-offset forms.  Epsilon-Laurent-polynomial Nhat entries use the
     established shift/matrix representation.  Nonconstant denominators use
     finite grouped numerators plus a causal coefficient recurrence. *)
  (* d0 is an epsilon unit by the cleared-system contract.  In the common
     (and banana) case it is actually epsilon-independent.  Cache the
     reciprocal of the already-expanded coefficient so the recursion keeps
     ratEpsList's exact/2x-WP numericization choice and never constructs a
     dense framed inverse for a scalar.  A genuinely epsilon-dependent d0
     retains the general frInv/frConv path. *)
  d0Expr = Cancel[Together[data["dExpr"][[1]]]];
  d0InvScalar = If[FreeQ[d0Expr, eps],
    Module[{i0 = 1 - fb, c0},
      If[i0 < 1 || i0 > W,
        err["E4", cs, <|"FrameBase" -> fb, "FrameTop" -> fb + W - 1,
          "Detail" -> "epsilon-independent d0 lies outside the work frame"|>]];
      c0 = dL[[1, i0]];
      If[zeroCanQ[c0],
        err["E5", cs, <|"Detail" -> "cleared d0 is zero"|>]];
      If[NumericQ[c0], 1/c0, Together[1/c0]]],
    None];
  Module[{dSp},
    dSp = Table[Module[{lst = dL[[j + 1]]},
      Table[{i + fb - 1, lst[[i]]}, {i, Select[Range[W], lst[[#]] =!= 0 &]}]],
      {j, 0, dD}];
    <|"dL" -> dL, "dD" -> dD, "dN" -> dN,
      "d0InvScalar" -> d0InvScalar,
      "dSp" -> dSp,
      "NhatSp" -> nhatPrep["PolynomialSp"],
      "NhatRationalGroups" -> ratGroups,
      "NhatRationalDenominators" -> ratDenominators,
      "NhatValuations" -> nhatPrep["Valuations"],
      "NhatStats" -> nhatPrep["Stats"]|>]];

prepareEpsilonRegularCleared[cs_Association, fb_Integer, W_Integer,
    certificate_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], physical, blocks, jordan,
   d0, spectralExpr, spectral, source},
  If[!TrueQ[Lookup[certificate, "Eligible", False]],
    err["E6", cs, <|"Certificate" -> KeyDrop[certificate, "Symbolic"],
      "Detail" -> "epsilon-regular recurrence preparation requires an exact eligible certificate"|>]];
  physical = prepareCleared[cs, fb, W, certificate["Symbolic"]];
  blocks = blockList[cs];
  jordan = ConstantArray[0, {cs["SystemSize"], cs["SystemSize"]}];
  Do[
    Do[jordan[[column, column]] =
        Cancel[Together[block["a"] + block["b"] eps]],
      {column, block["Cols"]}];
    Do[jordan[[block["Cols", position],
          block["Cols", position + 1]]] = 1,
      {position, 1, Length[block["Cols"]] - 1}],
    {block, blocks}];
  d0 = First[certificate["Symbolic", "dExpr"]];
  spectralExpr = Map[Cancel[Together[d0 #]] &, jordan, {2}];
  spectral = prepareNhatHybrid[{spectralExpr}, eps, fb, W, cs];
  source = prepareFramedMatrix[cs["VInv"], eps, fb, W, cs];
  Join[physical, <|
    "EpsilonRegularPrincipal" -> True,
    "EpsilonRegularCertificate" -> KeyDrop[certificate, "Symbolic"],
    "SpectralPrincipal" -> <|
      "PolynomialSp" -> First[spectral["PolynomialSp"]],
      "RationalGroups" -> First[spectral["RationalGroups"]],
      "Valuations" -> First[spectral["Valuations"]]|>,
    "SpectralSource" -> source|>]];

(* EpsSeries -> frame list *)
esToFrame[s_, fb_, W_] := Module[{out = ConstantArray[0, W], m1, m2},
  m1 = esMin[s]; m2 = esCM[s];
  Do[If[1 <= k - fb + 1 <= W, out[[k - fb + 1]] = esCoeff[s, k]], {k, m1, m2}];
  out];

(* Completeness shadow helpers.  Infinity means a structurally exact zero,
   which never constrains a sum.  Finite integers are the highest eps power
   known complete for one framed scalar. *)
validShift[Infinity, _Integer, _Integer] := Infinity;
validShift[k_Integer, s_Integer, top_Integer] := Min[top, k + s];
validMin[x_List] := If[x === {}, Infinity, Min @@ x];
frameValuation[x_List, fb_Integer] := Module[{p},
  p = SelectFirst[Range[Length[x]], x[[#]] =!= 0 &, None];
  If[p === None, Infinity, fb + p - 1]];

(* The old validity loop visited every nonzero epsilon coefficient of every
   Nhat entry.  validShift is monotone in the shift, so its minimum is exactly
   the entry's epsilon valuation.  This O(d^2) update is coefficient/window
   identical and avoids restoring an O(W) loop beside the grouped product. *)
updateNhatValidity[accValid_List, valuations_?MatrixQ,
    inputValid_List, frameTop_Integer] := Table[
  Min[accValid[[r]], validMin[Table[
    If[IntegerQ[valuations[[r, c]]],
      validShift[inputValid[[c]], valuations[[r, c]], frameTop],
      Infinity],
    {c, Length[inputValid]}]]],
  {r, Length[accValid]}];

(* Exact O(W) division by the affine spectral offset dA + dB eps.
   CASE T (dA != 0) is a coefficient recurrence; CASE P is an exact
   Laurent shift.  No tolerance-based structural decisions occur here. *)
affineDivFrame[rhs_List, dA_, dB_, fb_Integer, W_Integer] := Module[
  {top = fb + W - 1, invA, invB, y, i},
  If[Length[rhs] =!= W,
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "affine frame divisor received the wrong width"|>]];
  Which[
    !zeroCanQ[dA],
      If[!(fb <= 0 <= top),
        err["E4", <|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "FrameTop" -> top,
            "Detail" -> "CASE-T affine inverse exponent zero is outside the work frame"|>]];
      invA = If[NumericQ[dA], 1/dA, Together[1/dA]];
      If[zeroCanQ[dB],
        invA*rhs,
        y = ConstantArray[0, W];
        y[[1]] = invA*rhs[[1]];
        Do[y[[i]] = invA*(rhs[[i]] - dB*y[[i - 1]]), {i, 2, W}];
        y],
    !zeroCanQ[dB],
      If[!(fb <= -1 && top >= 1),
        err["E4", <|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "FrameTop" -> top,
            "Detail" -> "CASE-P epsilon inverse is outside the work frame"|>]];
      invB = If[NumericQ[dB], 1/dB, Together[1/dB]];
      invB*frDivEps[rhs],
    True,
      err["E5", <|"Center" -> "(frame)"|>,
        <|"Detail" -> "zero affine divisor belongs to CASE R"|>]]];

applyInvD0Frame[v_List, invD0_, d0InvScalar_, fb_Integer, W_Integer] :=
  If[d0InvScalar === None, frConv[v, invD0, fb, W], d0InvScalar*v];

(* Solve d0(eps) ((dA+dB eps) I - N_q) y = rhs.  N_q has a unit
   superdiagonal, so the affine/Jordan solve is bottom-up and d0^-1 is
   applied ONCE per row (not once per Jordan step).  This replaces the old
   dense powers of a truncated framed inverse with q linear recurrences. *)
blockSolveTPFrame[rhs_List, dA_, dB_, invD0_, d0InvScalar_, q_Integer,
    fb_Integer, W_Integer] := Module[{top = fb + W - 1, z, r, j, v},
  If[q < 1 || Length[rhs] =!= q || !AllTrue[rhs, ListQ[#] && Length[#] === W &],
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "malformed affine/Jordan block frame"|>]];
  If[zeroCanQ[dA] && zeroCanQ[dB],
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "zero affine divisor belongs to CASE R"|>]];
  If[zeroCanQ[dA],
    (* Preserve the strict legacy underflow contract before bottom-up sums
       can conceal an individually out-of-frame rhs_j/eps^j term.  The old
       implementation also built inv^q even for an exact-zero RHS. *)
    If[fb > -q || top < 1,
      lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
        <|"FrameBase" -> fb, "FrameTop" -> top, "JordanSize" -> q,
          "Detail" -> "CASE-P Jordan inverse exceeds the lower work frame"|>]];
    Do[
      v = frameValuation[rhs[[j]], fb];
      If[IntegerQ[v] && v - j < fb,
        lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "RHSComponent" -> j,
            "RHSValuation" -> v, "RequiredShift" -> -j,
            "Detail" -> "CASE-P Jordan solve would discard nonzero lower-frame content"|>]],
      {j, q}]];
  z = Table[ConstantArray[0, W], {q}];
  z[[q]] = affineDivFrame[rhs[[q]], dA, dB, fb, W];
  Do[z[[r]] = affineDivFrame[rhs[[r]] + z[[r + 1]], dA, dB, fb, W],
    {r, q - 1, 1, -1}];
  applyInvD0Frame[#, invD0, d0InvScalar, fb, W] & /@ z];

(* the recursion on frames.  UValid mirrors U without an eps axis:
   UValid[[n+1,l+1,r]] is the honest complete top of that scalar, or
   Infinity for a structural exact zero.  This shadow follows the actual
   dependency DAG; pseudo/log divisions at parallel blocks or log levels
   therefore take a maximum rather than being globally counted. *)
runRecursionFramedSrc[args___] := runRecursion[args];
runRecursionWolfram[cs_, prep_, aT_, bT_, P_, nmax_, srcHat_, fb_, W_, init_] := Module[
  {d = cs["SystemSize"], dL = prep["dL"],
   dD = prep["dD"], dN = prep["dN"], invD0, d0InvScalar, blocks, U,
   hits = {}, n0,
   UValid, topValid, frameTop = fb + W - 1, lamL, zeroV, srcData},
  d0InvScalar = Lookup[prep, "d0InvScalar", None];
  invD0 = If[d0InvScalar === None, frInv[dL[[1]], fb, W], None];
  blocks = blockList[cs];
  zeroV = ConstantArray[0, W];
  lamL[m_] := Module[{l = zeroV}, l[[-fb + 1]] = Together[aT + m];
    If[-fb + 2 <= W, l[[-fb + 2]] = Together[bT]]; l];
  srcData[n_, l_] := If[srcHat === None, None,
    Module[{sv = srcHat[n, l]},
      Which[
        sv === None, None,
        AssociationQ[sv] && KeyExistsQ[sv, "Frames"] && KeyExistsQ[sv, "Validity"],
          {sv["Frames"], sv["Validity"]},
        (* already frame lists (vectors of plain lists of length W) *)
        VectorQ[First[sv], AtomQ] || (ListQ[First[sv]] && Length[First[sv]] === W &&
          !AssociationQ[First[sv]]), {sv, ConstantArray[frameTop, d]},
        True, {Map[esToFrame[#, fb, W] &, sv], esCM /@ sv}]]];
  U = Table[Table[Table[zeroV, {d}], {l, 0, P + 1}], {n, 0, nmax}];
  UValid = ConstantArray[Infinity, {nmax + 1, P + 2, d}];
  n0 = If[init === None, 0, 1];
  If[init =!= None,
    Do[
      U[[1, l + 1]] = Map[esToFrame[#, fb, W] &, init[[l + 1]]];
      (* Seed from the input contracts, including finite zero series. *)
      UValid[[1, l + 1]] = esCM /@ init[[l + 1]],
      {l, 0, Min[P, Length[init] - 1]}]];
  (* shift a (d x W) block by s eps-slots (s > 0: content moves UP) *)
  Module[{shB, dSp = prep["dSp"], NhatSp = prep["NhatSp"],
      NhatRationalGroups = prep["NhatRationalGroups"],
      NhatRationalDenominators = prep["NhatRationalDenominators"],
      NhatValuations = prep["NhatValuations"], fuseRational},
  fuseRational = NhatRationalDenominators =!= {} &&
    !TrueQ[$disableRationalDenominatorFusion];
  shB[M_, s_] := shiftFrameBlock[M, s, fb, W, cs];
  Do[Module[{RData, R, RValid, rhsFull, rhsValid},
    RData = Table[Module[{acc = ConstantArray[0, {d, W}],
        accValid = ConstantArray[Infinity, d], rationalRHS, rationalUsed,
        denominatorIndex},
      rationalRHS = If[fuseRational,
        Table[ConstantArray[0, {d, W}], {Length[NhatRationalDenominators]}], {}];
      rationalUsed = If[fuseRational,
        ConstantArray[False, Length[NhatRationalDenominators]], {}];
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]]},
          Do[
            acc += matrixShiftProduct[
              sp[[2]], Ublk, sp[[1]], fb, W, cs],
            {sp, NhatSp[[j + 1]]}];
          If[NhatRationalGroups[[j + 1]] =!= {},
            If[fuseRational,
              Do[
                denominatorIndex = group["DenominatorIndex"];
                (* Accumulate only after the per-contribution numerator
                   application has performed its own lower-frame check. *)
                rationalRHS[[denominatorIndex]] +=
                  rationalMatrixGroupNumerator[group, Ublk, fb, W, cs];
                rationalUsed[[denominatorIndex]] = True,
                {group, NhatRationalGroups[[j + 1]]}],
              acc += applyRationalMatrixGroups[
                NhatRationalGroups[[j + 1]], Ublk, fb, W, cs]]];
          accValid = updateNhatValidity[accValid,
            NhatValuations[[j + 1]],
            UValid[[n - j + 1, l + 1]], frameTop]]],
        {j, 1, Min[n, dN]}];
      If[fuseRational,
        Do[If[TrueQ[rationalUsed[[denominatorIndex]]],
          acc += divideRationalMatrixRHS[
            rationalRHS[[denominatorIndex]],
            NhatRationalDenominators[[denominatorIndex,
              "DenominatorCoefficients"]], W]],
          {denominatorIndex, Length[NhatRationalDenominators]}]];
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]], term, a0, termValid},
          (* lam*U + eps*U_above: lam = (aT + n - j) + bT*eps *)
          a0 = Together[aT + n - j];
          term = a0*Ublk + shB[Together[bT]*Ublk, 1] +
            shB[U[[n - j + 1, l + 2]], 1];
          termValid = Table[validMin[{
              If[zeroCanQ[a0], Infinity, UValid[[n - j + 1, l + 1, r]]],
              If[zeroCanQ[bT], Infinity,
                validShift[UValid[[n - j + 1, l + 1, r]], 1, frameTop]],
              validShift[UValid[[n - j + 1, l + 2, r]], 1, frameTop]}],
            {r, d}];
          Do[
            acc -= sp[[2]]*shB[term, sp[[1]]];
            Do[accValid[[r]] = Min[accValid[[r]],
                validShift[termValid[[r]], sp[[1]], frameTop]], {r, d}],
            {sp, dSp[[j + 1]]}]]],
        {j, 1, Min[n, dD]}];
      If[srcHat =!= None,
        Do[If[n - j >= 0,
          Module[{sd = srcData[n - j, l], sv, svValid},
            If[sd =!= None,
              sv = sd[[1]]; svValid = sd[[2]];
              Do[
                acc += sp[[2]]*shB[sv, sp[[1]]];
                Do[accValid[[r]] = Min[accValid[[r]],
                    validShift[svValid[[r]], sp[[1]], frameTop]], {r, d}],
                {sp, dSp[[j + 1]]}]]]],
          {j, 0, Min[n, dD]}]];
      {acc, accValid}], {l, 0, P}];
    R = RData[[All, 1]];
    RValid = RData[[All, 2]];
    (* solve top-down in l *)
    Do[Module[{},
      rhsFull = R[[l + 1]];
      rhsValid = RValid[[l + 1]];
      Do[
        rhsFull -= sp[[2]]*shB[shB[U[[n + 1, l + 2]], 1], sp[[1]]];
        Do[rhsValid[[r]] = Min[rhsValid[[r]],
            validShift[UValid[[n + 1, l + 2, r]], 1 + sp[[1]], frameTop]],
          {r, d}],
        {sp, dSp[[1]]}];
      rhsFull = Table[rhsFull[[r]], {r, d}];
      Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"],
          dA, dB, solvedValid},
        dA = Together[aT + n - aI]; dB = Together[bT - bI];
        Which[
          !zeroCanQ[dA],
          U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], dA, dB,
              invD0, d0InvScalar, q, fb, W]]];
          solvedValid = Table[validMin[rhsValid[[colsB[[r ;; q]]]]], {r, q}];
          UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
            Thread[colsB -> solvedValid]],
          !zeroCanQ[dB],
          (U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], dA, dB,
              invD0, d0InvScalar, q, fb, W]]];
          solvedValid = Table[validMin[Table[
              validShift[rhsValid[[colsB[[r + m]]]], -(m + 1), frameTop],
              {m, 0, q - r}]], {r, q}];
          UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
            Thread[colsB -> solvedValid]];
          If[l == 0, AppendTo[hits, <|"n" -> n, "Cols" -> colsB,
            "DeltaB" -> dB, "FrameBase" -> fb,
            "GammaFrames" -> U[[n + 1, l + 1, colsB]],
            "GammaValidity" -> UValid[[n + 1, l + 1, colsB]]|>]]),
          True, Null]],
        {blk, blocks}]],
      {l, P, 0, -1}];
    (* CASE R ladder *)
    Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"]},
      If[zeroQ[aT + n - aI] && zeroQ[bT - bI],
        Module[{Rt, RtValid, assigned = Association[], assignedValid = Association[]},
          Rt = Table[Table[
            applyInvD0Frame[R[[l + 1, colsB[[r]]]], invD0, d0InvScalar,
              fb, W], {r, q}], {l, 0, P}];
          RtValid = RValid[[All, colsB]];
          Do[If[l + 1 <= P,
            assigned[{l + 1, q}] = frDivEps[Rt[[l + 1, q]]];
            assignedValid[{l + 1, q}] =
              validShift[RtValid[[l + 1, q]], -1, frameTop]],
            {l, 0, P - 1}];
          Do[
            Do[Module[{upr = Lookup[assigned, Key[{l + 1, r}], zeroV], uprValid},
              uprValid = Lookup[assignedValid, Key[{l + 1, r}], Infinity];
              If[!KeyExistsQ[assigned, {l, r + 1}],
                assigned[{l, r + 1}] = frShiftUp[upr] - Rt[[l + 1, r]];
                assignedValid[{l, r + 1}] = Min[
                  validShift[uprValid, 1, frameTop], RtValid[[l + 1, r]]]]],
              {r, 1, q - 1}],
            {l, P, 0, -1}];
          Do[If[KeyExistsQ[assigned, {l, r}],
            U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
              colsB[[r]] -> assigned[{l, r}]];
            UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
              colsB[[r]] -> assignedValid[{l, r}]]],
            {l, 0, P}, {r, 1, q}]]]],
      {blk, blocks}]],
    {n, n0, nmax}]];
  topValid = Module[{f = Select[Flatten[UValid], IntegerQ]},
    If[f === {}, frameTop, Min @@ f]];
  <|"U" -> U, "Validity" -> UValid, "Hits" -> hits, "P" -> P,
    "FrameBase" -> fb, "TopValid" -> topValid|>];

(* ---- compiled framed-recurrence seam ---------------------------------
   Wolfram retains every exact structural decision.  In particular the
   T/P/R schedule is serialized explicitly; the Acb backend never decides
   resonance from a numerical enclosure.  A forced Cpp backend is strict:
   unsupported coefficient fields or bridge failures are loud and never
   fall back to runRecursionWolfram. *)
$cppExactDomain = False;  (* focused exact-parity seam; Acb is production *)
$cppBuildRequestOnly = False;
$cppSerializationDomain = "acb";
$cppSerializationSymbols = {};
$cppStaticRecordOverride = None;
$cppUsePersistentSessions = True;

(* Input serialization and native arithmetic have different contracts.
   $InputPrecisionFactor protects decimal/algebraic input encoding; it must
   not silently double the precision of every Arb multiplication.  Native
   output is requested at WorkingPrecision+20 decimal digits, so retain that
   complete output width plus a fixed 32-bit arithmetic guard.  Arb parses
   the more precise serialized input into an enclosing ball at this work
   precision, preserving rigor while avoiding an accidental 2x-WP runtime. *)
$cppNativeOutputGuardDigits = 20;
$cppNativeArithmeticGuardBits = 32;

cppNativeOutputDigits[wp_Integer?Positive] :=
  wp + $cppNativeOutputGuardDigits;

cppNativePrecisionBits[wp_Integer?Positive] := Max[64,
  Ceiling[cppNativeOutputDigits[wp]*Log[2, 10]] +
    $cppNativeArithmeticGuardBits];

(* Schema-2 stores the exact SCC certificate with every retained native
   operator.  Wolfram's certificate is one-based and counts vertices on the
   longest path; the C++ protocol is zero-based and counts coupling edges. *)
cppPersistentSCC[cs_Association] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None]},
  If[!AssociationQ[seq],
    err["E6", cs, <|"Detail" ->
      "persistent C++ chart is missing its exact SCC certificate"|>]];
  <|"components" -> ((# - 1) & /@ seq["Components"]),
    "structural_edges" -> ((# - 1) & /@ seq["DependencyEdges"]),
    "condensation_edges" -> ((# - 1) & /@ seq["CondensationEdges"]),
    "topological_order" -> (seq["TopologicalOrder"] - 1),
    "coupling_depth" -> Max[0, seq["CouplingDepth"] - 1]|>];

cppPersistentPrescription[record_Association] := <|
  "factor_exact" -> ToString[
    Lookup[record, "ExactFactor", record["Factor"]], InputForm],
  "sign" -> record["Sign"],
  "multiplicity" -> record["Multiplicity"],
  "leading_coefficient_sign" -> record["LeadingCoeffSign"]|>;

cppExactRealAlgebraicQ[value_] := Module[{canonical},
  If[!FreeQ[value, _?InexactNumberQ], Return[False, Module]];
  canonical = Quiet[Check[RootReduce[value], $Failed]];
  canonical =!= $Failed && NumericQ[canonical] &&
    TrueQ[Quiet[Check[Element[canonical, Algebraics], False]]] &&
    TrueQ[Quiet[Check[FullSimplify[Im[canonical] == 0], False]]]];

cppExactAlgebraicTruthQ[condition_] := TrueQ[Quiet[Check[
  FullSimplify[RootReduce[condition]], False]]];

cppPersistentGeometry[cs_Association] := Module[
  {digits = Max[80, 2 DiffExp2`Config`CFG["WorkingPrecision"]],
   center = RootReduce[cs["Center"]],
   scale = RootReduce[cs["ChartMap", "Scale"]], radius, encoded},
  radius = If[cs["Radius"] === Infinity, Infinity,
    RootReduce[cs["Radius"]]];
  encoded[value_, label_] := Module[{result =
      DiffExp2`CppBackend`EncodeScalar[value, digits]},
    If[FailureQ[result], err["E5", cs, <|"Field" -> label,
      "Value" -> value, "Detail" ->
        "persistent chart geometry has no rigorous Acb specialization"|>]];
    result];
  Join[<|
    "center_exact" -> ToString[center, InputForm],
    "scale_exact" -> ToString[scale, InputForm],
    "center_numeric" -> encoded[center, "center"],
    "scale_numeric" -> encoded[scale, "scale"],
    "infinite_radius" -> TrueQ[cs["Radius"] === Infinity],
    "prescriptions" ->
      (cppPersistentPrescription /@ Lookup[cs, "Prescriptions", {}])|>,
  If[TrueQ[cs["Radius"] === Infinity], <||>,
    <|"radius_exact" -> ToString[radius, InputForm],
      "radius_numeric" -> encoded[radius, "radius"]|>]]];

(* Exact producer certificates for the deliberately narrow first composite
   slice.  These are structural predicates, never numerical enclosure tests.
   C++ binds the exact V identity below to the retained assembly operator and
   binds the separately prepared VInv entries to the target-source transform;
   the inverse and determinant facts remain collision-bound exact producer
   certificates. *)
sccExactIdentityMatrixQ[matrix_, dimension_Integer] := Module[{delta},
  If[!MatrixQ[matrix] || Dimensions[matrix] =!= {dimension, dimension} ||
      !FreeQ[matrix, _?InexactNumberQ], Return[False, Module]];
  delta = Map[Cancel[Together[#]] &,
    Normal[matrix] - IdentityMatrix[dimension], {2}];
  AllTrue[Flatten[delta], # === 0 &]];

sccSpectralFrameCertificate[cs_Association] := Module[
  {dimension = Lookup[cs, "SystemSize", None],
   t = Lookup[cs, "ChartVar", None],
   eps = DiffExp2`Config`CanonicalEps[], v, vInv, det,
   detValuation, left, right, identityV, fail},
  fail[detail_, extra_:<||>] := Join[<|"admissible" -> False,
      "identity_v" -> False, "detail" -> detail|>, extra];
  If[!IntegerQ[dimension] || dimension < 1 || !MatchQ[t, _Symbol],
    Return[fail[
      "spectral-frame certification requires a positive dimension and exact chart variable"],
      Module]];
  v = Quiet[Check[Normal[Lookup[cs, "V", None]], $Failed]];
  vInv = Quiet[Check[Normal[Lookup[cs, "VInv", None]], $Failed]];
  If[v === $Failed || vInv === $Failed || !MatrixQ[v] ||
      !MatrixQ[vInv] || Dimensions[v] =!= {dimension, dimension} ||
      Dimensions[vInv] =!= {dimension, dimension},
    Return[fail["V/VInv are not exact square matrices of the block dimension"],
      Module]];
  If[!FreeQ[{v, vInv}, _?InexactNumberQ],
    Return[fail["V/VInv contain inexact coefficients"], Module]];
  If[!FreeQ[{v, vInv}, t],
    Return[fail[
      "persistent CompositeSCC spectral frames must be independent of the local chart variable"],
      Module]];
  left = Quiet[Check[Map[Cancel[Together[#]] &,
      v . vInv - IdentityMatrix[dimension], {2}], $Failed]];
  right = Quiet[Check[Map[Cancel[Together[#]] &,
      vInv . v - IdentityMatrix[dimension], {2}], $Failed]];
  If[left === $Failed || right === $Failed ||
      !AllTrue[Flatten[left], # === 0 &] ||
      !AllTrue[Flatten[right], # === 0 &],
    Return[fail["V and VInv are not exact two-sided inverses"], Module]];
  det = Quiet[Check[Cancel[Together[Det[v]]], $Failed]];
  detValuation = If[det === $Failed, $Failed,
    Quiet[Check[epsValuation[det, eps], $Failed]]];
  If[det === $Failed || !IntegerQ[detValuation],
    Return[fail[
      "persistent CompositeSCC requires a finite Laurent determinant valuation",
      <|"det_epsilon_valuation" -> detValuation|>], Module]];
  identityV = sccExactIdentityMatrixQ[v, dimension] &&
    sccExactIdentityMatrixQ[vInv, dimension];
  <|"admissible" -> True, "identity_v" -> identityV,
    "dimension" -> dimension,
    "det_epsilon_valuation" -> detValuation,
    "v_exact_identity" ->
      DiffExp2`SectorSeries`ExactExpressionIdentity[v, t],
    "vinv_exact_identity" ->
      DiffExp2`SectorSeries`ExactExpressionIdentity[vInv, t],
    "det_exact_identity" ->
      DiffExp2`SectorSeries`ExactExpressionIdentity[det, t]|>];

sccGaugeFrameCertificate[cs_Association] := Module[
  {dimension = Lookup[cs, "SystemSize", None],
   t = Lookup[cs, "ChartVar", None], gauge, inverse, reduction,
   identityGauge, gaugeIdentity, inverseIdentity, pairIdentity,
   matrixIdentity, fail},
  fail[detail_, extra_:<||>] := Join[<|"admissible" -> False,
      "identity_gauge" -> False, "detail" -> detail|>, extra];
  (* Preserve SparseArray structure here.  The SCC gauge is commonly an
     identity plus one rational-function column; densifying it before the
     exact two-sided product caused multi-gigabyte intermediate expansion. *)
  gauge = Quiet[Check[Lookup[cs, "Gauge", None], $Failed]];
  inverse = Quiet[Check[Lookup[cs, "GaugeInverse", None], $Failed]];
  If[!IntegerQ[dimension] || dimension < 1 || !MatchQ[t, _Symbol] ||
      gauge === $Failed || inverse === $Failed ||
      !MatrixQ[gauge] || !MatrixQ[inverse] ||
      Dimensions[gauge] =!= {dimension, dimension} ||
      Dimensions[inverse] =!= {dimension, dimension},
    Return[fail["Gauge/GaugeInverse are not exact square matrices"],
      Module]];
  If[!FreeQ[{gauge, inverse}, _?InexactNumberQ],
    Return[fail["Gauge/GaugeInverse contain inexact coefficients"], Module]];
  reduction = Lookup[Lookup[cs, "IndicialData", <||>],
    "Reduction", <||>];
  If[!TrueQ[Lookup[reduction, "GaugeInverseCertified", False]] ||
      Lookup[reduction, "GaugeInverseCertificateSchema", None] =!=
        "diffexp2-indicial-exact-gauge-inverse-v1",
    Return[fail[
      "Gauge/GaugeInverse lack the exact indicial producer certificate"],
      Module]];
  (* FuchsianReduce proved T.TInv==I over the exact rational-function
     field before retaining this chart.  For a square matrix this proves
     invertibility and the reverse identity.  Repeating dense products here
     caused multi-gigabyte expression swell on the L1 banana gauge. *)
  identityGauge = sccExactIdentityMatrixQ[gauge, dimension] &&
    sccExactIdentityMatrixQ[inverse, dimension];
  matrixIdentity[matrix_, role_] := ExportString[<|
      "schema" -> "diffexp2-exact-gauge-matrix-v1", "role" -> role,
      "dimension" -> dimension,
      "entries" -> Map[
        DiffExp2`SectorSeries`ExactExpressionIdentity[#, t] &,
        Normal[matrix], {2}]|>, "RawJSON", "Compact" -> True];
  gaugeIdentity = matrixIdentity[gauge, "Gauge"];
  inverseIdentity = matrixIdentity[inverse, "GaugeInverse"];
  pairIdentity = "de2-exact-two-sided-gauge-" <>
    IntegerString[Hash[{gaugeIdentity, inverseIdentity}, "SHA256"], 16, 64];
  <|"admissible" -> True, "identity_gauge" -> identityGauge,
    "dimension" -> dimension,
    "gauge_exact_identity" -> gaugeIdentity,
    "gauge_inverse_exact_identity" -> inverseIdentity,
    "gauge_det_exact_identity" -> pairIdentity|>];

sccNoFamilyCollisionQ[cs_Association] := Module[
  {families = Lookup[cs, "Families", None], collisions, depths},
  If[!ListQ[families], Return[False, Module]];
  collisions = Flatten[Lookup[families, "Collisions", {}], 1];
  depths = Lookup[families, "CollisionDepth", {}];
  collisions === {} && AllTrue[depths, # === 0 &]];

sccNativeBlockCapabilities[cs_Association,
    spectralFrame_:Automatic, gaugeFrameInput_:Automatic] := Module[
  {dimension = Lookup[cs, "SystemSize", 0], frame = spectralFrame,
   gaugeFrame = gaugeFrameInput},
  If[frame === Automatic, frame = sccSpectralFrameCertificate[cs]];
  If[gaugeFrame === Automatic, gaugeFrame = sccGaugeFrameCertificate[cs]];
  <|"regular" -> TrueQ[Lookup[
      Lookup[cs, "IndicialData", <||>], "Regular", False]],
    "identity_gauge" ->
      (sccExactIdentityMatrixQ[Lookup[cs, "Gauge", None], dimension] &&
       sccExactIdentityMatrixQ[
         Lookup[cs, "GaugeInverse", None], dimension]),
    "exact_gauge" -> TrueQ[Lookup[gaugeFrame, "admissible", False]],
    "identity_v" -> TrueQ[Lookup[frame, "identity_v", False]],
    "epsilon_unimodular_v" ->
      TrueQ[Lookup[frame, "admissible", False]],
    "no_pseudo" -> sccNoFamilyCollisionQ[cs]|>];


sccBlockPrincipalMatrixRecord[cs_Association] := Module[
  {clearKey = Lookup[cs, "SystemClearKey", None], parentInput, matrix,
   variable, dimension = Lookup[cs, "SystemSize", None]},
  If[!IntegerQ[Lookup[cs, "SCCBlock", None]] || clearKey === None ||
      !KeyExistsQ[$systemClearRegistry, clearKey],
    err["E6", cs, <|"SystemClearKey" -> clearKey,
      "Detail" -> "native SCC principal metadata requires a registered diagonal block chart"|>]];
  parentInput = $systemClearRegistry[clearKey];
  matrix = Lookup[parentInput, "Matrix", None];
  variable = Lookup[parentInput, "Variable", None];
  If[!MatchQ[variable, _Symbol] ||
      Dimensions[matrix] =!= {dimension, dimension},
    err["E6", cs, <|"Dimensions" -> Dimensions[matrix],
      "Expected" -> {dimension, dimension},
      "Detail" -> "native SCC principal matrix registry entry is malformed"|>]];
  sccExactMatrixRecord[matrix, variable, cs]];

$cppRegularValueAggregationGuardDigits = 0;
$cppRegularValueMatchingCertificationDigits = Automatic;

(* Regular value recentering is an accuracy-producing operation and must
   follow the explicit per-level matching target.  ResidTol is only a
   WP-derived ODE spot-check threshold; using Max[ResidTol, MatchTol] here
   silently made every producer retry above Floor[WP/10] ineffective. *)
cppRegularValueRelativeAccuracyMaxExact[] := ToString[
  If[IntegerQ[$cppRegularValueMatchingCertificationDigits],
    10^-$cppRegularValueMatchingCertificationDigits,
    DiffExp2`Tolerances`Tol["MatchTol"]]/
    10^(DiffExp2`Tolerances`$SafetyDigits +
      $cppRegularValueAggregationGuardDigits), InputForm];

cppPersistentMetadata[cs_Association, fb_Integer, W_Integer] := Module[
  {systemIdentity, chartIdentity, chartAnalytic},
  (* SolveCacheTag joins diagonal SCC blocks back under the parent level's
     one native session.  Ordinary charts use their collision-certified
     exact system key.  The complete prepared operator is independently
     stored in the chart signature, so this grouping cannot alias unequal
     recurrence tensors. *)
  systemIdentity = Lookup[cs, "SolveCacheTag",
    Lookup[cs, "SystemClearKey", Lookup[cs, "SystemHash", None]]];
  chartIdentity = {Lookup[cs, "SystemClearKey", None],
    Lookup[cs, "Center", None], Lookup[cs, "ChartMap", None],
    Lookup[cs, "SCCBlock", None], fb, W};
  chartAnalytic = <|
    "PrescriptionIdentity" ->
      ToString[Lookup[cs, "Prescriptions", {}], InputForm],
    "Prescriptions" ->
      (cppPersistentPrescription /@ Lookup[cs, "Prescriptions", {}]),
    "regular_value_relative_accuracy_max_exact" ->
      cppRegularValueRelativeAccuracyMaxExact[],
    "geometry" -> cppPersistentGeometry[cs]|>;
  If[IntegerQ[Lookup[cs, "SCCBlock", None]],
    chartAnalytic = Join[chartAnalytic, <|
      "principal_matrix" -> sccBlockPrincipalMatrixRecord[cs],
      "native_scc_capabilities" -> sccNativeBlockCapabilities[cs]|>]];
  <|"SystemIdentity" -> systemIdentity,
    "ChartIdentity" -> chartIdentity,
    "SessionAnalytic" -> <|
      "RegulatorSymbols" -> (SymbolName /@ $cppSerializationSymbols),
      "Policy" -> "exact-structure-prescription-specialized"|>,
    "ChartAnalytic" -> chartAnalytic,
    "SCC" -> cppPersistentSCC[cs]|>];

(* Exact metadata for session-owned native LocalSolutions.  Unlike the Acb
   recurrence payload, these fields are structural facts: canonical forms,
   integer/sign predicates, chart geometry, and branch prescriptions are
   never inferred from a numerical enclosure. *)
cppExactTruth[value_] := Which[TrueQ[value], "yes",
  TrueQ[!value], "no", True, "unknown"];
cppExactSign[value_] := Which[TrueQ[value === -1], "negative",
  TrueQ[value === 0], "zero", TrueQ[value === 1], "positive",
  True, "unknown"];

cppNativeExactDescriptor[value_, inputDigits_Integer, cs_] := Module[
  {canonical, rationalQ, zeroFact, integerFact, signFact, specialization},
  If[!FreeQ[value, _?InexactNumberQ] ||
      !FreeQ[value, DiffExp2`Config`CanonicalEps[]],
    err["E5", cs, <|"Tag" -> value, "Detail" ->
      "native local sector tags must be exact and epsilon-independent"|>]];
  canonical = Quiet[Check[RootReduce[Together[value]], $Failed]];
  If[canonical === $Failed || !NumericQ[canonical],
    err["E5", cs, <|"Tag" -> value, "Detail" ->
      "native local sector tag contains an unresolved analytic regulator"|>]];
  If[!zeroCanQ[RootReduce[Im[canonical]]],
    err["E5", cs, <|"Tag" -> value, "Detail" ->
      "native local sector tags must be exact real scalars"|>]];
  canonical = RootReduce[Re[canonical]];
  rationalQ = IntegerQ[canonical] || Head[canonical] === Rational;
  If[!rationalQ && !TrueQ[Quiet[FullSimplify[
      Element[canonical, Algebraics]]]],
    err["E5", cs, <|"Tag" -> canonical, "Detail" ->
      "native local exact tag lies outside the rational/algebraic descriptor domains"|>]];
  zeroFact = zeroCanQ[canonical];
  integerFact = Quiet[Check[FullSimplify[Element[canonical, Integers]],
    Indeterminate]];
  signFact = Quiet[Check[Sign[canonical], Indeterminate]];
  specialization = DiffExp2`CppBackend`EncodeScalar[canonical, inputDigits];
  If[FailureQ[specialization],
    err["E5", cs, <|"Tag" -> canonical,
      "BackendFailure" -> specialization,
      "Detail" -> "native local exact tag has no Acb specialization"|>]];
  If[rationalQ,
    <|"domain" -> "rational", "canonical" ->
        ToString[canonical, InputForm],
      "is_zero" -> cppExactTruth[zeroFact],
      "is_integer" -> cppExactTruth[integerFact],
      "sign" -> cppExactSign[signFact],
      "specialization" -> specialization|>,
    <|"domain" -> "algebraic", "canonical" ->
        ToString[canonical, InputForm],
      "is_zero" -> cppExactTruth[zeroFact],
      "is_integer" -> cppExactTruth[integerFact],
      "sign" -> cppExactSign[signFact],
      "specialization" -> specialization|>]];

cppNativeChartMetadata[cs_Association, inputDigits_Integer] := Module[
  {center = cs["Center"], scale = cs["ChartMap", "Scale"],
   radius = cs["Radius"], encodedRadius, realRadius},
  If[!FreeQ[{center, scale}, _?InexactNumberQ] ||
      !AllTrue[{center, scale}, NumericQ],
    err["E5", cs, <|"Geometry" -> {center, scale}, "Detail" ->
      "native local chart center and scale must be exact numeric scalars"|>]];
  If[!zeroCanQ[RootReduce[Im[center]]] ||
      !zeroCanQ[RootReduce[Im[scale]]] || zeroCanQ[scale],
    err["E5", cs, <|"Geometry" -> {center, scale}, "Detail" ->
      "native local chart requires real center and nonzero real scale"|>]];
  If[radius === Infinity,
    <|"center_exact" -> ToString[RootReduce[Re[center]], InputForm],
      "scale_exact" -> ToString[RootReduce[Re[scale]], InputForm],
      "infinite_radius" -> True|>,
    If[!NumericQ[radius] || !TrueQ[N[radius, 30] > 0] ||
        !zeroCanQ[RootReduce[Im[radius]]],
      err["E5", cs, <|"Radius" -> radius, "Detail" ->
        "native local chart radius must be a positive real scalar"|>]];
    (* Retain an exact rational radius as the native scalar string itself.
       Besides avoiding an unnecessary Acb pair, this makes local metadata
       byte-for-byte comparable with a composite's exact retained geometry.
       Algebraic/non-rational callers keep the ordinary enclosing encoding. *)
    realRadius = RootReduce[Re[radius]];
    encodedRadius = If[IntegerQ[realRadius] || Head[realRadius] === Rational,
      ToString[realRadius, InputForm],
      DiffExp2`CppBackend`EncodeScalar[radius, inputDigits]];
    If[FailureQ[encodedRadius],
      err["E5", cs, <|"Radius" -> radius,
        "BackendFailure" -> encodedRadius,
        "Detail" -> "native local chart radius is not Acb-encodable"|>]];
    <|"center_exact" -> ToString[RootReduce[Re[center]], InputForm],
      "scale_exact" -> ToString[RootReduce[Re[scale]], InputForm],
      (* The numerical radius constructs ChartGeometry, while this exact
         identity binds algebraic radii to a retained SCC parent.  Rational
         radii keep both fields too, so one metadata schema covers both. *)
      "radius_exact" -> ToString[realRadius, InputForm],
      "radius" -> encodedRadius, "infinite_radius" -> False|>]];

cppNativePrescriptions[cs_Association] := Map[Function[record, Module[
    {factor = Lookup[record, "ExactFactor", Lookup[record, "Factor", None]],
     sign = Lookup[record, "Sign", None],
     multiplicity = Lookup[record, "Multiplicity", None],
     leading = Lookup[record, "LeadingCoeffSign", None]},
    If[factor === None || !FreeQ[factor, _?InexactNumberQ] ||
        !MemberQ[{-1, 1}, sign] || !IntegerQ[multiplicity] ||
        multiplicity < 1 || !MemberQ[{-1, 1}, leading],
      err["E5", cs, <|"Prescription" -> record, "Detail" ->
        "native local analytic-continuation prescription is malformed"|>]];
    <|"factor_exact" -> ToString[factor, InputForm], "sign" -> sign,
      "multiplicity" -> multiplicity,
      "leading_coefficient_sign" -> leading|>]],
  Lookup[cs, "Prescriptions", {}]];

cppNativeLocalMetadata[cs_Association, a_, b_, p_Integer,
    inputDigits_Integer, checkpointIdentity_String] := <|
  "chart" -> cppNativeChartMetadata[cs, inputDigits],
  "tag" -> <|"a" -> cppNativeExactDescriptor[a, inputDigits, cs],
    "b" -> cppNativeExactDescriptor[b, inputDigits, cs],
    (* The native parser binds a,b to the recurrence targets.  p is already
       a strict nonnegative integer in run[\"p\"]; retain its canonical fact
       beside them for checkpoint/audit consumers. *)
    "p" -> <|"domain" -> "integer", "canonical" -> ToString[p]|>|>,
  "prescriptions" -> cppNativePrescriptions[cs],
  "checkpoint_identity" -> checkpointIdentity|>;

cppScalar[e_, digits_Integer, cs_] := Module[{encoded},
  encoded = If[$cppSerializationDomain === "symbolic",
    DiffExp2`CppBackend`EncodeSymbolicScalar[e, $cppSerializationSymbols],
    DiffExp2`CppBackend`EncodeScalar[e, digits]];
  If[FailureQ[encoded],
    err["E5", cs, <|"Coefficient" -> e,
      "BackendFailure" -> encoded,
      "Detail" -> "C++ recurrence cannot represent a prepared coefficient"|>]];
  If[$cppSerializationDomain === "rational",
    If[!((IntegerQ[e] || Head[e] === Rational) && TrueQ[Im[e] === 0]),
      err["E5", cs, <|"Coefficient" -> e,
        "Detail" -> "exact C++ parity mode requires rational coefficients"|>]];
    First[encoded],
    encoded]];

cppValidity[k_] := If[k === Infinity, Null, k];

cppMatrixShift[sp_List, digits_Integer, cs_] := Module[
  {shift = sp[[1]], matrix = Normal[sp[[2]]], d, entries},
  d = Length[matrix];
  entries = Flatten[Table[
    If[matrix[[r, c]] =!= 0,
      {{r - 1, c - 1, cppScalar[matrix[[r, c]], digits, cs]}}, {}],
    {r, d}, {c, d}], 2];
  <|"s" -> shift, "e" -> entries|>];

cppPhysicalRational[data_Association, digits_Integer, cs_] := If[
  TrueQ[data["Zero"]],
  <|"zero" -> True|>,
  <|"zero" -> False, "valuation" -> data["Valuation"],
    "numerator" -> (cppScalar[#, digits, cs] & /@
      CoefficientList[data["P"], DiffExp2`Config`CanonicalEps[]]),
    "denominator" -> (cppScalar[#, digits, cs] & /@
      CoefficientList[data["Q"], DiffExp2`Config`CanonicalEps[]])|>];

cppPhysicalODEPayload[data_Association, ownerIdentity_String,
    digits_Integer, cs_] := Module[{d = cs["SystemSize"], c},
  c = Map[Function[lag, Flatten[Table[
      If[TrueQ[lag[[r, col, "Zero"]]], Nothing,
        <|"r" -> r - 1, "c" -> col - 1,
          "v" -> cppPhysicalRational[lag[[r, col]], digits, cs]|>],
      {r, d}, {col, d}]]], data["C"]];
  <|"schema" -> "diffexp2-physical-cleared-ode-v1",
    "basis" -> "physical-original-master",
    "theta_coordinate" -> "local-t",
    "owner_signature_identity" -> ownerIdentity,
    "payload_identity" -> data["Identity"],
    "q" -> (cppPhysicalRational[#, digits, cs] & /@ data["Q"]),
    "c" -> c|>];

$cppStaticOperatorCache = <||>;
$cppStaticOperatorCacheMax = 1024;

(* Decimal encoding of every prepared lag dominated the Wolfram side of
   repeated persistent solves if rebuilt per homogeneous column/SCC source.
   Cache the complete schema-1-compatible static payload once per exact
   prepared frame.  The hash is only an index and every hit compares the full
   structural signature; the opaque token then gives CppBackend a lightweight
   identity for the already collision-certified payload. *)
cppStaticOperatorPayload[cs_, prep_, blocks_List, fb_Integer, W_Integer,
    vPrep_, inputDigits_Integer, precisionBits_Integer,
    physicalDataOverride_:Automatic] := Module[
  {d = cs["SystemSize"], signature, key, cached, dLags, denominators,
   nLags, assembly = Null, assemblyGroups, assemblyBase, payload, record,
   physicalData, physicalPayload, token, assemblyExactIdentity,
   epsilonRegularQ, spectralPrincipal = Null, spectralSource = Null,
   encodeAdditionalMatrix, spectralSourceExactIdentity},
  physicalData = If[AssociationQ[physicalDataOverride],
    physicalDataOverride, physicalClearedODEData[cs]];
  (* Validate an owner-supplied record before it can enter either the static
     operator identity or the native physical-equation payload. *)
  If[AssociationQ[physicalDataOverride],
    regularClearedSymbolicFromPhysicalData[cs, physicalData]];
  assemblyExactIdentity = If[AssociationQ[vPrep],
    DiffExp2`SectorSeries`ExactExpressionIdentity[
      Normal[cs["V"]], cs["ChartVar"]], None];
  signature = {$cppSerializationDomain, $cppSerializationSymbols,
    inputDigits, precisionBits, d, fb, W, prep, blocks,
    If[AssociationQ[vPrep], vPrep, Automatic], assemblyExactIdentity,
    physicalData,
    cfg["ChopPrecision"]};
  key = Hash[signature, "SHA256"];
  token = "de2-operator-" <> IntegerString[key, 16, 64];
  cached = Lookup[$cppStaticOperatorCache, key, None];
  If[AssociationQ[cached] && SameQ[cached["Signature"], signature],
    Return[cached, Module]];
  If[cached =!= None,
    err["E5", cs, <|"Detail" ->
      "C++ static operator cache-key collision with unequal full identity"|>]];
  dLags = Map[Function[lag, Map[Function[sp,
      <|"s" -> sp[[1]], "v" -> cppScalar[sp[[2]], inputDigits, cs]|>], lag]],
    prep["dSp"]];
  denominators = Map[Function[den,
      cppScalar[#, inputDigits, cs] & /@ den["DenominatorCoefficients"]],
    prep["NhatRationalDenominators"]];
  encodeAdditionalMatrix[prepared_Association, identity_:False,
      exactIdentity_:None] := Module[{groups, base, encoded},
    groups = Lookup[prepared, "RationalGroups", {}];
    base = Length[denominators];
    denominators = Join[denominators,
      Map[(cppScalar[#, inputDigits, cs] & /@
          #["DenominatorCoefficients"]) &, groups]];
    encoded = <|
      "poly" -> (cppMatrixShift[#, inputDigits, cs] & /@
        prepared["PolynomialSp"]),
      "rat" -> MapIndexed[Function[{group, idx}, <|
        "q" -> base + First[idx] - 1,
        "num" -> (cppMatrixShift[#, inputDigits, cs] & /@
          group["NumeratorSp"])|>], groups],
      "val" -> (cppValidity /@ Flatten[prepared["Valuations"]])|>;
    If[BooleanQ[identity], AssociateTo[encoded, "identity" -> identity]];
    If[StringQ[exactIdentity],
      AssociateTo[encoded, "exact_identity" -> exactIdentity]];
    encoded];
  epsilonRegularQ = TrueQ[Lookup[prep,
    "EpsilonRegularPrincipal", False]];
  If[epsilonRegularQ,
    spectralPrincipal = encodeAdditionalMatrix[
      prep["SpectralPrincipal"]];
    spectralSourceExactIdentity =
      DiffExp2`SectorSeries`ExactExpressionIdentity[
        Normal[cs["VInv"]], cs["ChartVar"]];
    spectralSource = encodeAdditionalMatrix[prep["SpectralSource"],
      TrueQ[prep["SpectralSource", "Identity"]],
      spectralSourceExactIdentity]];
  If[AssociationQ[vPrep],
    assemblyGroups = vPrep["RationalGroups"];
    assemblyBase = Length[denominators];
    denominators = Join[denominators,
      Map[(cppScalar[#, inputDigits, cs] & /@
          #["DenominatorCoefficients"]) &, assemblyGroups]];
    assembly = <|
      "identity" -> TrueQ[Lookup[vPrep, "Identity", False]],
      "exact_identity" -> assemblyExactIdentity,
      "poly" -> (cppMatrixShift[#, inputDigits, cs] & /@
        vPrep["PolynomialSp"]),
      "rat" -> MapIndexed[Function[{group, idx}, <|
        "q" -> assemblyBase + First[idx] - 1,
        "num" -> (cppMatrixShift[#, inputDigits, cs] & /@
          group["NumeratorSp"])|>], assemblyGroups],
      "val" -> (cppValidity /@ Flatten[vPrep["Valuations"]])|>];
  nLags = MapThread[Function[{poly, groups, vals}, <|
      "poly" -> (cppMatrixShift[#, inputDigits, cs] & /@ poly),
      "rat" -> Map[Function[group, <|
        "q" -> group["DenominatorIndex"] - 1,
        "num" -> (cppMatrixShift[#, inputDigits, cs] & /@
          group["NumeratorSp"])|>], groups],
      "val" -> (cppValidity /@ Flatten[vals])|>],
    {prep["NhatSp"], prep["NhatRationalGroups"],
      prep["NhatValuations"]}];
  physicalPayload = cppPhysicalODEPayload[
    physicalData, token, inputDigits, cs];
  payload = <|"domain" -> $cppSerializationDomain,
    "symbols" -> (SymbolName /@ $cppSerializationSymbols),
    "precision_bits" -> precisionBits,
    "d" -> d, "fb" -> fb, "w" -> W,
    "d_lags" -> dLags, "denominators" -> denominators,
    "nhat_lags" -> nLags,
    "epsilon_regular_principal" -> epsilonRegularQ,
    "spectral_principal" -> spectralPrincipal,
    "spectral_source" -> spectralSource,
    "d0_inverse" -> If[prep["d0InvScalar"] === None, Null,
      cppScalar[prep["d0InvScalar"], inputDigits, cs]],
    "blocks" -> Map[(#["Cols"] - 1) &, blocks],
    "assembly" -> assembly, "physical_ode" -> physicalPayload,
    "chop_digits" -> cfg["ChopPrecision"]|>;
  (* Stable content identity prevents an evicted Wolfram serialization cache
     entry from duplicating an already-retained native chart.  The digest is
     only an index: both caches retain and compare the complete certificate. *)
  record = <|"Signature" -> signature, "Payload" -> payload,
    "Token" -> token|>;
  If[Length[$cppStaticOperatorCache] >= $cppStaticOperatorCacheMax,
    KeyDropFrom[$cppStaticOperatorCache,
      First[Keys[$cppStaticOperatorCache]]]];
  AssociateTo[$cppStaticOperatorCache, key -> record];
  record];

cppRunRecursionCore[cs_, prep_, aT_, bT_, P_Integer, nmax_Integer,
    srcHat_, fb_Integer, W_Integer, init_, vPrep_:Automatic,
    responseOverride_:Automatic] := Module[
  {d = cs["SystemSize"], wp = cfg["WorkingPrecision"], inputDigits,
   outputDigits, precisionBits, blocks, schedule, staticRecord,
   initFrames, initValidity, request, response, decodedU,
   decodedValidity, hits, decode, topValid,
   timingQ = Environment["DEBUG_CPP_RECURRENCE"] === "1", t0, tRequest,
   tCall, tDone,
   assembledData = None, sourceRecords, sourcePayload},
  t0 = SessionTime[];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  outputDigits = cppNativeOutputDigits[wp];
  precisionBits = cppNativePrecisionBits[wp];
  blocks = blockList[cs];
  schedule = Table[Map[Function[blk, Module[{dA, dB, kind},
      dA = Together[aT + n - blk["a"]];
      dB = Together[bT - blk["b"]];
      kind = Which[!zeroCanQ[dA], "T", !zeroCanQ[dB], "P", True, "R"];
      <|"case" -> kind, "da" -> cppScalar[dA, inputDigits, cs],
        "db" -> cppScalar[dB, inputDigits, cs]|>]], blocks],
    {n, 0, nmax}];
  staticRecord = If[AssociationQ[$cppStaticRecordOverride],
    $cppStaticRecordOverride,
    cppStaticOperatorPayload[cs, prep, blocks, fb, W,
      vPrep, inputDigits, precisionBits]];
  initFrames = Table[If[init =!= None && l + 1 <= Length[init],
      Map[esToFrame[#, fb, W] &, init[[l + 1]]],
      ConstantArray[0, {d, W}]], {l, 0, P}];
  initValidity = Table[If[init =!= None && l + 1 <= Length[init],
      esCM /@ init[[l + 1]], ConstantArray[Infinity, d]], {l, 0, P}];
  sourcePayload = If[srcHat === None, Null,
    sourceRecords = Table[Module[{sv = srcHat[n, l], frames, validity},
      Which[
        sv === None,
          <|"Present" -> False, "Frames" -> ConstantArray[0, {d, W}],
            "Validity" -> ConstantArray[Infinity, d]|>,
        AssociationQ[sv] && KeyExistsQ[sv, "Frames"] &&
            KeyExistsQ[sv, "Validity"],
          <|"Present" -> True, "Frames" -> sv["Frames"],
            "Validity" -> sv["Validity"]|>,
        ListQ[sv] && Length[sv] === d &&
            AllTrue[sv, ListQ[#] && Length[#] === W &],
          <|"Present" -> True, "Frames" -> sv,
            "Validity" -> ConstantArray[fb + W - 1, d]|>,
        ListQ[sv] && Length[sv] === d &&
            AllTrue[sv, DiffExp2`EpsSeries`ESQ],
          <|"Present" -> True,
            "Frames" -> (esToFrame[#, fb, W] & /@ sv),
            "Validity" -> (esCM /@ sv)|>,
        True,
          err["E5", cs, <|"SourceIndex" -> {n, l},
            "Detail" -> "C++ recurrence received a malformed materialized source frame"|>]]],
      {n, 0, nmax}, {l, 0, P}];
    <|"frames" -> (cppScalar[#, inputDigits, cs] & /@
        Flatten[Map[# ["Frames"] &, sourceRecords, {2}]]),
      "validity" -> (cppValidity /@
        Flatten[Map[# ["Validity"] &, sourceRecords, {2}]]),
      "present" -> Flatten[Map[TrueQ[# ["Present"]] &, sourceRecords, {2}]]|>];
  request = Join[staticRecord["Payload"], <|
    "schema" -> 1, "output_digits" -> outputDigits,
    "nmax" -> nmax, "p" -> P,
    "has_initial" -> (init =!= None),
    "adaptive_probe" -> TrueQ[$adaptiveLowerFrameProbe],
    "cancellation_audit_base" -> If[
      IntegerQ[$cancellationAuditBase], $cancellationAuditBase, Null],
    "a_target" -> cppScalar[aT, inputDigits, cs],
    "b_target" -> cppScalar[bT, inputDigits, cs],
    "a_shift_min" -> 0,
    "a_shifts" -> Table[cppScalar[Together[aT + m], inputDigits, cs],
      {m, 0, nmax}],
    "schedule" -> schedule,
    "initial" -> (cppScalar[#, inputDigits, cs] & /@ Flatten[initFrames]),
    "initial_validity" -> (cppValidity /@ Flatten[initValidity]),
    "source" -> sourcePayload,
    "return_u" -> !AssociationQ[vPrep]|>];
  If[TrueQ[$cppBuildRequestOnly], Return[request, Module]];
  tRequest = SessionTime[];
  response = If[AssociationQ[responseOverride], responseOverride,
    If[TrueQ[$cppUsePersistentSessions] &&
        Environment["DE2_CPP_PERSISTENT"] =!= "0",
      DiffExp2`CppBackend`RunPersistentRequest[request,
        Append[cppPersistentMetadata[cs, fb, W],
          "PreparedToken" -> staticRecord["Token"]]],
      DiffExp2`CppBackend`RunRequest[request]]];
  tCall = SessionTime[];
  If[FailureQ[response],
    err["E5", cs, <|"BackendFailure" -> response,
      "Detail" -> "compiled recurrence backend call failed"|>]];
  If[Lookup[response, "status", "error"] =!= "ok",
    If[TrueQ[$adaptiveLowerFrameProbe] &&
        Lookup[response, "id", ""] === "E4" &&
        AnyTrue[{"below the work frame", "lower-frame", "lowest framed",
          "epsilon shift", "framed convolution", "pseudo Jordan"},
          StringContainsQ[Lookup[response, "detail", ""], #] &],
      Throw[Failure["AdaptiveLowerFrame", <|
        "ID" -> "AdaptiveLowerFrame", "FrameBase" -> fb,
        "Shift" -> Lookup[response, "shift", 0],
        "Detail" -> Lookup[response, "detail",
          "compiled recurrence lower-frame underflow"]|>],
        "DiffExp2AdaptiveLowerFrame"]];
    err[If[MemberQ[{"E4", "E5", "E6"}, Lookup[response, "id", ""]],
        response["id"], "E5"], cs,
      <|"BackendID" -> Lookup[response, "id", "CPP"],
        "FrameBase" -> Lookup[response, "frame_base", fb],
        "Shift" -> Lookup[response, "shift", 0],
        "Detail" -> Lookup[response, "detail",
          "compiled recurrence returned an error"]|>]];
  decode = DiffExp2`CppBackend`DecodeScalar[#, outputDigits] &;
  If[KeyExistsQ[response, "assembled"],
    Module[{amin = response["assembled", "min"],
        amax = response["assembled", "max"], coeffs},
      coeffs = DiffExp2`CppBackend`DecodeScalars[
        response["assembled", "coefficients"], outputDigits];
      If[AnyTrue[coeffs, FailureQ],
        err["E5", cs, <|"Detail" ->
          "compiled recurrence returned an undecodable assembled coefficient"|>]];
      assembledData = <|"Min" -> amin, "CompleteMax" -> amax,
        "Coefficients" -> ArrayReshape[coeffs,
          {P + 1, amax - amin + 1, nmax + 1, d}]|>],
    decodedU = DiffExp2`CppBackend`DecodeScalars[
      response["u"], outputDigits];
    If[AnyTrue[decodedU, FailureQ],
      err["E5", cs, <|"Detail" ->
        "compiled recurrence returned an undecodable coefficient"|>]];
    decodedU = ArrayReshape[decodedU, {nmax + 1, P + 2, d, W}];
    decodedValidity = ArrayReshape[
      Replace[response["validity"], Null -> Infinity, {1}],
      {nmax + 1, P + 2, d}]];
  hits = Map[Function[hit, Module[{hf, hv, cols},
      cols = hit["cols"] + 1;
      hf = ArrayReshape[DiffExp2`CppBackend`DecodeScalars[
        hit["frames"], outputDigits], {Length[cols], W}];
      hv = Replace[hit["validity"], Null -> Infinity, {1}];
      <|"n" -> hit["n"], "Cols" -> cols,
        "DeltaB" -> decode[hit["delta_b"]], "FrameBase" -> fb,
        "GammaFrames" -> hf, "GammaValidity" -> hv|>]],
    response["hits"]];
  topValid = Replace[response["top_valid"], Null -> Infinity];
  tDone = SessionTime[];
  If[timingQ, Print["CPPREC ", <|
    "Chart" -> ToString[Lookup[cs, "Center", "?"], InputForm],
    "d" -> d, "nmax" -> nmax, "P" -> P, "W" -> W,
    "PrepareSeconds" -> N[tRequest - t0, 6],
    "BridgeSeconds" -> N[tCall - tRequest, 6],
    "DecodeSeconds" -> N[tDone - tCall, 6],
    "KernelMilliseconds" -> response["elapsed_ms"]|>]];
  Join[<|"Hits" -> hits,
    "P" -> P, "FrameBase" -> fb, "TopValid" -> topValid,
    "BackendDiagnostics" -> Join[<|"Backend" -> "Cpp",
      "KernelMilliseconds" -> response["elapsed_ms"]|>,
      If[AssociationQ[Lookup[response, "persistent", None]],
        <|"Persistent" -> response["persistent"]|>, <||>]]|>,
    If[assembledData === None,
      <|"U" -> decodedU, "Validity" -> decodedValidity|>,
      <|"CppAssembled" -> assembledData|>]]];

cppRegulatorSymbols[cs_, prep_, aT_, bT_, P_Integer, nmax_Integer,
    srcHat_, init_, vPrep_] := Module[{sourceProbe, data, symbols},
  sourceProbe = If[srcHat === None, {},
    Table[Quiet[srcHat[n, l]], {n, 0, nmax}, {l, 0, P}]];
  data = {aT, bT, prep["dL"], prep["NhatSp"],
    prep["NhatRationalDenominators"], init, sourceProbe,
    If[AssociationQ[vPrep], {vPrep["PolynomialSp"],
      vPrep["RationalGroups"]}, {}]};
  symbols = DeleteDuplicates[Cases[data,
    s_Symbol /; Context[s] === "Global`" &&
      s =!= DiffExp2`Config`CanonicalEps[] &&
      s =!= Lookup[cs, "ChartVar", None] && !NumericQ[s], Infinity]];
  SortBy[symbols, SymbolName]];

cppRunRecursion[cs_, prep_, aT_, bT_, P_Integer, nmax_Integer,
    srcHat_, fb_Integer, W_Integer, init_, vPrep_:Automatic,
    responseOverride_:Automatic] := Module[{symbols, domain},
  symbols = cppRegulatorSymbols[cs, prep, aT, bT, P, nmax,
    srcHat, init, vPrep];
  domain = Which[symbols =!= {}, "symbolic",
    TrueQ[$cppExactDomain], "rational", True, "acb"];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    cppRunRecursionCore[cs, prep, aT, bT, P, nmax, srcHat, fb, W,
      init, vPrep, responseOverride]]];

(* A Wolfram kernel cannot advance two TransportLine evaluations while a
   LibraryLink call is blocking.  Homogeneous chart bases are nevertheless
   independent of the incoming boundary, so the two endpoint arms can
   collect their already-prepared recurrence requests and execute those in
   one native task pool.  Capture/injection is deliberately below every
   exact structural decision: the normal solve is run twice (request pass,
   then verified assembly pass), and only the native response is supplied
   on the second pass.  Thus resonance, frame, analytic-regulator, residual,
   and strict no-fallback semantics remain owned by the ordinary path. *)
$cppHomogeneousBatchCapture = False;
$cppHomogeneousBatchInjection = None;
$cppHomogeneousBatchInjectionUses = 0;
$cppHomogeneousBatchTag = "DiffExp2CppHomogeneousBatch";

cppConfiguredThreads[count_Integer] := Module[{threads},
  threads = Quiet[Check[ToExpression[Environment["DE2_CPP_THREADS"]], 4]];
  If[!IntegerQ[threads] || threads < 1, threads = 4];
  Min[threads, Max[1, count]]];

cppBatchRecurrences[cs_, prep_, tasks_List, nmax_Integer, fb_Integer,
    W_Integer, vPrep_] := Module[
  {requests, threads, response, results, started = SessionTime[], wp,
   inputDigits, precisionBits, staticRecord, metadata, symbols, domain},
  threads = cppConfiguredThreads[Length[tasks]];
  (* One native batch has one coefficient field.  Derive that field from the
     union of every run before entering the serialization Block and keep the
     Block live through request construction, persistent preparation, replay,
     and decoded assembly.  In particular, never recompute a symbolic static
     operator after cppRunRecursion's per-call Block has unwound to Acb. *)
  symbols = SortBy[DeleteDuplicates[Flatten[
    cppRegulatorSymbols[cs, prep, # ["a"], # ["b"], # ["P"], nmax,
      None, # ["Init"], vPrep] & /@ tasks]], SymbolName];
  domain = Which[symbols =!= {}, "symbolic",
    TrueQ[$cppExactDomain], "rational", True, "acb"];
  wp = cfg["WorkingPrecision"];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  precisionBits = cppNativePrecisionBits[wp];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    (* Build/certify the immutable operator once.  Every task then shares the
       same association and token instead of re-hashing prep/vPrep. *)
    staticRecord = cppStaticOperatorPayload[cs, prep, blockList[cs], fb, W,
      vPrep, inputDigits, precisionBits];
    Block[{$cppStaticRecordOverride = staticRecord},
      requests = Map[Function[task,
        Block[{$cppBuildRequestOnly = True},
          cppRunRecursionCore[cs, prep, task["a"], task["b"], task["P"],
            nmax, None, fb, W, task["Init"], vPrep]]], tasks];
      metadata = Append[cppPersistentMetadata[cs, fb, W],
        "PreparedToken" -> staticRecord["Token"]];
      If[TrueQ[$cppHomogeneousBatchCapture],
        Throw[<|"Requests" -> requests,
          (* Exact producer-side task identities survive static-payload
             compaction.  Later persistent SCC planning must not recover
             roots or log ceilings from Acb strings/midpoints. *)
          "TaskMetadata" ->
            (KeyTake[#, {"a", "b", "P"}] & /@ tasks),
          "Metadata" -> metadata|>,
          $cppHomogeneousBatchTag]];
      response = If[AssociationQ[$cppHomogeneousBatchInjection],
        If[Lookup[$cppHomogeneousBatchInjection, "Requests", None] =!= requests ||
            !ListQ[Lookup[$cppHomogeneousBatchInjection, "Results", None]],
          err["E5", cs, <|
            "Detail" -> "prewarmed C++ recurrence request did not reproduce exactly during verified assembly"|>]];
        $cppHomogeneousBatchInjectionUses++;
        <|"status" -> "ok",
          "results" -> $cppHomogeneousBatchInjection["Results"]|>,
        (* Capture/replay above never touches the native session. *)
        If[TrueQ[$cppUsePersistentSessions] &&
            Environment["DE2_CPP_PERSISTENT"] =!= "0",
          DiffExp2`CppBackend`RunPersistentRequests[
            requests, metadata, threads],
          DiffExp2`CppBackend`RunRequest[
            <|"batch" -> requests, "threads" -> threads|>]]];
      If[FailureQ[response] || Lookup[response, "status", "error"] =!= "ok" ||
          !ListQ[Lookup[response, "results", None]] ||
          Length[response["results"]] =!= Length[tasks],
        err["E5", cs, <|"BackendResponse" -> response,
          "Detail" -> "compiled recurrence batch call failed"|>]];
      results = MapThread[Function[{task, raw},
        cppRunRecursionCore[cs, prep, task["a"], task["b"], task["P"],
          nmax, None, fb, W, task["Init"], vPrep, raw]],
        {tasks, response["results"]}];
      If[Environment["DEBUG_CPP_RECURRENCE"] === "1",
        Print["CPPBATCH ", <|"Tasks" -> Length[tasks], "Threads" -> threads,
          "FrameBase" -> fb, "FrameWidth" -> W,
          "Seconds" -> N[SessionTime[] - started, 6]|>]];
      results]]];

runRecursion[cs_, prep_, aT_, bT_, P_, nmax_, srcHat_, fb_, W_, init_,
    vPrep_:Automatic] :=
  If[cfg["RecurrenceBackend"] === "Cpp",
    cppRunRecursion[cs, prep, aT, bT, P, nmax, srcHat, fb, W, init, vPrep],
    runRecursionWolfram[cs, prep, aT, bT, P, nmax, srcHat, fb, W, init]];

(* ---- assembly ---- *)

(* frame U -> original-frame LocalSolution with sectors (a, b, l) *)
assembleSolution[cs_, aT_, bT_, rec_, nmax_, vPrep_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], d = cs["SystemSize"], U = rec["U"],
   UValid = rec["Validity"], P = rec["P"], fb = rec["FrameBase"], W,
   frameTop, kminO, kmaxO, VexpL, VVal, WL, WValid, secs, ls, firstNZ, vp,
   assembled},
  If[KeyExistsQ[rec, "CppAssembled"],
    assembled = rec["CppAssembled"];
    kminO = assembled["Min"];
    kmaxO = assembled["CompleteMax"];
    secs = Table[<|"a" -> aT, "b" -> bT, "p" -> l,
      "Coeffs" -> assembled["Coefficients"][[l + 1]]|>,
      {l, 0, rec["P"]}];
    ls = <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"], "Sectors" -> secs,
      "EpsWindow" -> <|"Min" -> kminO, "CompleteMax" -> kmaxO|>,
      "TWindow" -> <|"CompleteMax" -> nmax|>,
      "ErrorEstimate" -> ConstantArray[0, kmaxO - kminO + 1],
      "Prescriptions" -> cs["Prescriptions"]|>;
    ls = applyGauge[cs, ls, nmax];
    Return[DiffExp2`SectorSeries`CanonicalizeLocalSolution[ls], Module]];
  W = Length[U[[1, 1, 1]]];
  frameTop = fb + W - 1;
  If[TrueQ[$disableGroupedSpectralTransform],
    VexpL = Map[ratEpsList[Together[#], eps, fb, W] &, cs["V"], {2}];
    VVal = Map[frameValuation[#, fb] &, VexpL, {2}];
    WL = Table[
      Table[Sum[frConv[VexpL[[r, c]], U[[n + 1, l + 1, c]], fb, W],
        {c, d}], {r, d}], {n, 0, nmax}, {l, 0, P}],
    vp = If[AssociationQ[vPrep], vPrep,
      prepareFramedMatrix[cs["V"], eps, fb, W, cs]];
    VVal = vp["Valuations"];
    WL = Table[applyPreparedFramedMatrix[
      vp, U[[n + 1, l + 1]], fb, W, cs],
      {n, 0, nmax}, {l, 0, P}]];
  WValid = Table[
    Table[validMin[Table[
      If[VVal[[r, c]] === Infinity, Infinity,
        validShift[UValid[[n + 1, l + 1, c]], VVal[[r, c]], frameTop]],
      {c, d}]], {r, d}],
    {n, 0, nmax}, {l, 0, P}];
  kmaxO = Module[{f = Select[Flatten[WValid], IntegerQ]},
    If[f === {}, frameTop, Min @@ f]];
  (* content min: an all-zero result has no first coefficient and is
     represented canonically at its honest complete top. *)
  firstNZ = SelectFirst[Range[W],
    Function[i, AnyTrue[Flatten[WL, 2],
      Function[fr, !certifiedFrameZeroQ[fr[[i]]]]]], None];
  kminO = If[firstNZ === None, kmaxO, fb + firstNZ - 1];
  If[firstNZ =!= None && kmaxO < kminO,
    err["E6", cs, <|"FrameBase" -> fb, "FrameTop" -> frameTop,
      "ContentMin" -> kminO, "CompleteMax" -> kmaxO,
      "Tag" -> {aT, bT},
      "Detail" -> "recursion completeness is exhausted below the first stored coefficient"|>]];
  secs = If[firstNZ === None,
    Table[<|"a" -> aT, "b" -> bT, "p" -> l,
      "Coeffs" -> ConstantArray[0, {1, nmax + 1, d}]|>, {l, 0, P}],
    Table[<|"a" -> aT, "b" -> bT, "p" -> l,
      "Coeffs" -> Table[Table[Table[certifiedFrameChop[
        WL[[n + 1, l + 1, r, k - fb + 1]]], {r, d}],
        {n, 0, nmax}], {k, kminO, kmaxO}]|>, {l, 0, P}]];
  ls = <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"], "Sectors" -> secs,
    "EpsWindow" -> <|"Min" -> kminO, "CompleteMax" -> kmaxO|>,
    "TWindow" -> <|"CompleteMax" -> nmax|>,
    "ErrorEstimate" -> ConstantArray[0, kmaxO - kminO + 1],
    "Prescriptions" -> cs["Prescriptions"]|>;
  ls = applyGauge[cs, ls, nmax];
  DiffExp2`SectorSeries`CanonicalizeLocalSolution[ls]];

(* gauge composition f = T.g: T exact rational, poles only at t = 0.
   Frame-list implementation: the k-axis (eps) stays a plain index, the
   t-convolution runs per (np, r) with ListConvolve on eps-frame lists.
   A negative eps-lead in T pays the window top (same rule as the V
   multiply in assembleSolution). *)
applyGauge[cs_, ls_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], T = cs["Gauge"], d,
   gv, kmin, kmax, newSecs, vT, fbT, WG, TexpL, gaugeData, Tcoeffs},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[ls]];
  d = Length[T];
  gaugeData = gaugeCoefficientData[cs, nmax];
  gv = gaugeData["TMin"];
  Tcoeffs = gaugeData["Coefficients"];
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  vT = gaugeData["EpsValuation"];
  fbT = vT; WG = (kmax - kmin) + (-vT) + 1;
  (* TexpL[m - gv + 1][r][c] = eps-frame list (base fbT, width WG) of the
     t^m Laurent coefficient of T_rc *)
  TexpL = Map[
    If[zeroCanQ[#], ConstantArray[0, WG], ratEpsList[#, eps, fbT, WG]] &,
    Tcoeffs, {3}];
  newSecs = Map[Module[{arr = #["Coeffs"], aS = #["a"], outF, kminN, kmaxN,
      ncolsS = Dimensions[#["Coeffs"]][[2]], topValidG, firstNZ},
    (* arr columns as eps-frame lists on [kmin, kmin + WG - 1] (pad top) *)
    outF = Table[ConstantArray[0, WG], {nmax + 1}, {d}];
    Do[Module[{n = np - (m - gv)},
      If[0 <= n <= nmax && n < ncolsS,
        Do[Module[{arrCol = PadRight[arr[[All, n + 1, c]], WG], tl},
          Do[
            tl = TexpL[[m - gv + 1, r, c]];
            If[!AllTrue[tl, # === 0 &],
              outF[[np + 1, r]] += Take[
                ListConvolve[tl, arrCol, {1, -1}, 0], WG]],
            {r, d}]],
          {c, d}]]],
      {np, 0, nmax}, {m, gv, gv + nmax}];
    topValidG = kmax + vT;  (* pay |vT| at the top (vT <= 0) *)
    firstNZ = SelectFirst[Range[WG],
      AnyTrue[Flatten[outF[[All, All, #]]],
        Function[z, !certifiedFrameZeroQ[z]]] &, None];
    kminN = If[firstNZ === None, topValidG, kmin + fbT + firstNZ - 1];
    If[firstNZ =!= None && topValidG < kminN,
      err["E6", cs, <|"InputWindow" -> {kmin, kmax},
        "GaugeValuation" -> vT, "ContentMin" -> kminN,
        "CompleteMax" -> topValidG, "Tag" -> {aS, #["b"], #["p"]},
        "Detail" -> "gauge multiplication exhausted completeness below the first stored coefficient"|>]];
    kmaxN = topValidG;
    <|"a" -> Together[aS + gv], "b" -> #["b"], "p" -> #["p"],
      "Coeffs" -> Table[Table[Table[
        Module[{idx = k - (kmin + fbT) + 1},
          certifiedFrameChop[
            If[1 <= idx <= WG, outF[[n + 1, r, idx]], 0]]], {r, d}],
        {n, 0, nmax}], {k, kminN, kmaxN}],
      "KMin" -> kminN, "KMax" -> kmaxN|>] &, ls["Sectors"]];
  Module[{kminA = Min[#["KMin"] & /@ newSecs], kmaxA = Min[#["KMax"] & /@ newSecs]},
    newSecs = Map[Module[{sh = #["KMin"] - kminA},
      <|"a" -> #["a"], "b" -> #["b"], "p" -> #["p"],
        "Coeffs" -> Table[
          If[kminA + i - 1 < #["KMin"],
            Map[0 &, #["Coeffs"][[1]], {2}],
            #["Coeffs"][[i - sh]]],
          {i, 1, kmaxA - kminA + 1}]|>] &, newSecs];
    Join[ls, <|"Sectors" -> newSecs,
      "EpsWindow" -> <|"Min" -> kminA, "CompleteMax" -> kmaxA|>,
      "ErrorEstimate" -> ConstantArray[0, kmaxA - kminA + 1]|>]]];

(* Joint CASE-P compensation (Solve.md 3.7).  At a collision the raw
   quotient coefficient gamma along target block r multiplies the target
   root's leading chain vector.  Its negative-epsilon polar part is not a
   new solution: it must be cancelled by -gammaPolar times the already
   built target-family column.  Zeros above eps^-1 in that polynomial are
   structural and may be certified through the product's honest top.

   Target columns are kept in their wider work windows until every source
   column has been compensated.  A depth-one collision consumes one upper
   order in the target product; parallel collisions consume it only once.
   The deterministic CollisionDepth work budget guarantees the requested
   top survives. *)
compensatePseudoColumn[cs_, raw_Association, hits_List, workColumns_List,
    reqMax_Integer] := Module[{out = raw, records = {}},
  Do[Module[{fb = hit["FrameBase"], gammas = hit["GammaFrames"], cols = hit["Cols"]},
    Do[Module[{gamma = gammas[[r]], negIdx, targetIndex, target, wMin,
        targetMin, targetMax, desiredMax, wMax, weight},
      negIdx = Select[Range[Length[gamma]],
        fb + # - 1 < 0 && gamma[[#]] =!= 0 &];
      If[negIdx =!= {},
        targetIndex = cols[[r]];
        If[targetIndex > Length[workColumns],
          err["E5", cs, <|"TargetColumn" -> targetIndex,
            "BuiltColumns" -> Length[workColumns], "Hit" -> KeyDrop[hit,
              {"GammaFrames", "GammaValidity"}],
            "Detail" -> "pseudo-resonant compensation target is not yet built (family ordering violated)"|>]];
        target = workColumns[[targetIndex]];
        wMin = fb + First[negIdx] - 1;
        targetMin = target["EpsWindow", "Min"];
        targetMax = target["EpsWindow", "CompleteMax"];
        desiredMax = Min[out["EpsWindow", "CompleteMax"], targetMax + wMin];
        If[desiredMax < reqMax,
          err["E4", cs, <|"RequestedCompleteMax" -> reqMax,
            "AvailableCompleteMax" -> desiredMax, "WeightMin" -> wMin,
            "TargetWindow" -> target["EpsWindow"],
            "Detail" -> "pseudo-resonant compensation exhausted the deterministic work window"|>]];
        wMax = desiredMax - targetMin;
        weight = esNew[wMin, Table[
          If[k < 0 && 1 <= k - fb + 1 <= Length[gamma],
            -gamma[[k - fb + 1]], 0], {k, wMin, wMax}]];
        out = DiffExp2`SectorSeries`CombineLocalSolutions[
          {1, weight}, {out, target}];
        AppendTo[records, <|"n" -> hit["n"], "TargetColumn" -> targetIndex,
          "DeltaB" -> hit["DeltaB"], "WeightWindow" ->
            <|"Min" -> wMin, "CompleteMax" -> wMax|>|>]]],
      {r, Length[cols]}]],
    {hit, hits}];
  {out, records}];

pseudoHitPolarOrders[hit_Association] := If[KeyExistsQ[hit, "PolarOrders"],
  hit["PolarOrders"],
  Module[{fb = hit["FrameBase"], gammas = hit["GammaFrames"],
      cols = hit["Cols"]},
    DeleteDuplicates[Flatten[Table[
      Map[{cols[[r]], fb + # - 1} &,
        Select[Range[Length[gammas[[r]]]],
          fb + # - 1 < 0 && gammas[[r, #]] =!= 0 &]],
      {r, Length[cols]}], 1]]]];

publicPseudoHit[hit_Association] := Append[
  KeyDrop[hit, {"GammaFrames", "GammaValidity", "FrameBase"}],
  "PolarOrders" -> pseudoHitPolarOrders[hit]];

(* Spectral-coordinate probe retained as a diagnostic/test seam.  It is NOT
   the CASE-P certificate: when VInv is epsilon-singular, a finite physical
   solution can have legitimate polar spectral coordinates. *)
spectralProbeValue[cs_, ls_, t0_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], value, fcomps,
   H, hTop, out},
  value = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
    "UsePade" -> False, "ComputeTailEstimates" -> False]["Value"];
  fcomps = Table[esNew[esMin[value],
      Table[esCoeff[value, k][[c]], {k, esMin[value], esCM[value]}]],
    {c, cs["SystemSize"]}];
  H = Map[Cancel[Together[#]] &,
    cs["VInv"] . (cs["GaugeInverse"] /. t -> t0), {2}];
  hTop = Max[0, -1 - Min[esMin /@ fcomps]];
  out = Table[Module[{acc = None},
      Do[If[!zeroCanQ[H[[r, c]]],
        Module[{term = esTimes[esFrom[H[[r, c]], eps, hTop], fcomps[[c]]]},
          acc = If[acc === None, term, esAdd[acc, term]]]],
        {c, cs["SystemSize"]}];
      If[acc === None, esZero[esCM[value]], acc]],
    {r, cs["SystemSize"]}];
  out];

(* Gauge-reduced physical value g = GaugeInverse.f.  Unlike the individual
   coordinates of u = VInv.g, this frame does not acquire artificial Laurent
   poles from an epsilon-singular spectral inverse. *)
reducedProbeValue[cs_, ls_, t0_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], value, fcomps,
   H, hTop, out},
  value = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
    "UsePade" -> False, "ComputeTailEstimates" -> False]["Value"];
  fcomps = Table[esNew[esMin[value],
      Table[esCoeff[value, k][[c]], {k, esMin[value], esCM[value]}]],
    {c, cs["SystemSize"]}];
  H = Map[Cancel[Together[#]] &,
    cs["GaugeInverse"] /. t -> t0, {2}];
  hTop = Max[0, -1 - Min[esMin /@ fcomps]];
  out = Table[Module[{acc = None},
      Do[If[!zeroCanQ[H[[r, c]]],
        Module[{term = esTimes[esFrom[H[[r, c]], eps, hTop], fcomps[[c]]]},
          acc = If[acc === None, term, esAdd[acc, term]]]],
        {c, cs["SystemSize"]}];
      If[acc === None, esZero[esCM[value]], acc]],
    {r, cs["SystemSize"]}];
  out];

(* Map every registered spectral collision order through the exact V column.
   The certificate checks the corresponding polar orders in the reduced
   physical components.  Requiring the VInv coordinates themselves to be
   pole-free is wrong when VInv has epsilon poles: the banana L2 x=0 frame is
   the minimal counterexample. *)
reducedPolarOrders[cs_, hits_List] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], V = cs["V"], pairs},
  pairs = DeleteDuplicates[Flatten[pseudoHitPolarOrders /@ hits, 1]];
  DeleteDuplicates[Flatten[Table[Module[
      {target = pair[[1]], k = pair[[2]], v},
      Flatten[Table[
        v = epsValuation[V[[c, target]], eps];
        If[IntegerQ[v] && k + v < 0,
          Table[{c, ko, target}, {ko, k + v, -1}], {}],
        {c, cs["SystemSize"]}], 1]],
    {pair, pairs}], 1]]];

certifyPseudoCompensation[cs_, ls_, hits_List, label_] := Module[
  {polarPairs, probes, nmax, minA, tailOrder, fracRaw, frac, wp2,
   ltol = DiffExp2`Tolerances`Tol["LaurentLeadTol"]},
  If[hits === {}, Return[True]];
  polarPairs = reducedPolarOrders[cs, hits];
  If[polarPairs === {}, Return[True]];
  (* A compensated target whose integer a is shifted relative to the raw
     sector has one unmatched top Taylor tail at finite TOrder.  Probe deeply
     enough inside the disk that this known truncation tail lies below the
     structural floor; the two probes remain deterministic and nonzero. *)
  nmax = ls["TWindow", "CompleteMax"];
  minA = Min[Floor[#["a"]] & /@ ls["Sectors"]];
  tailOrder = Max[1, nmax + minA];
  fracRaw = N[(ltol/100)^(1/tailOrder), 40];
  frac = Min[1/5, Rationalize[fracRaw, Abs[fracRaw]/50]];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"];
  probes = N[{cs["Radius"]*frac/2, 2*cs["Radius"]*frac/3}, wp2];
  Do[Module[{gv = reducedProbeValue[cs, ls, t0]},
    Do[Module[{r = pair[[1]], k = pair[[2]], target = pair[[3]], allVals,
        scale, z, mag, uncertainty, acc},
      If[k > esCM[gv[[r]]],
        (* The source may not certify this higher polar order.  Its honest
           CompleteMax is retained, so skip rather than infer a value. *)
        Continue[]];
      allVals = Table[If[k < esMin[gv[[c]]], 0, esCoeff[gv[[c]], k]],
        {c, Length[gv]}];
      scale = Max[1, Sequence @@
        (numMag[#, 20] & /@ Select[Flatten[allVals], NumericQ])];
      z = If[k < esMin[gv[[r]]], 0, esCoeff[gv[[r]], k]];
      uncertainty = If[InexactNumberQ[z],
        acc = Accuracy[z];
        If[NumericQ[acc] && acc =!= Infinity, 10^-Floor[acc], 0], 0];
      Which[
        NumericQ[z],
          mag = numMag[z, 20];
          If[TrueQ[mag + uncertainty > ltol*scale],
            err["E5", cs, <|"Certificate" -> label, "Probe" -> t0,
              "EpsOrder" -> k, "TargetColumn" -> target,
              "ReducedComponent" -> r,
              "Residual" -> z, "ResidualUncertainty" -> uncertainty,
              "Scale" -> scale, "Tolerance" -> ltol,
              "Detail" -> "joint pseudo-resonant polar value did not cancel"|>]],
        !zeroQ[z],
          err["E5", cs, <|"Certificate" -> label, "Probe" -> t0,
            "EpsOrder" -> k, "TargetColumn" -> target,
            "ReducedComponent" -> r, "Residual" -> z,
            "Detail" -> "joint pseudo-resonant symbolic polar value did not cancel"|>],
        True, Null]],
      {pair, polarPairs}]],
    {t0, probes}];
  True];

(* ---- exact component/source algebra for SCC orchestration ------------
   Couplings are multiplied in the ORIGINAL theta-form frame.  These helpers
   deliberately reuse SectorSeries`MultiplyRational, so fractional powers,
   logarithms, epsilon Laurent windows, analytic prescriptions, and center
   pole shifts follow the same implementation as every other solver source.
   No coefficient is sampled or specialized to decide whether it is zero. *)
sccStructuralZeroQ[e_] :=
  FreeQ[e, _?InexactNumberQ] && zeroQ[e];

sccSymbolIdentity[s_Symbol] := <|"context" -> Context[s],
  "name" -> SymbolName[s]|>;

(* Coupling preparation must serialize in exactly the field retained by its
   parent native session.  Inference from Config`Variables is unsound after
   chart pruning: it can reorder or enlarge the field relative to the live
   session.  The present native codec exports bare names, so non-Global
   contexts, duplicate bare names, and collisions with the chart/epsilon
   binders are rejected until an explicit context-to-native-name map exists. *)
sccSerializationField[spec_, cs_Association] := Module[
  {field, domain, symbols, names, binders,
   eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"]},
  field = If[spec === Automatic,
    <|"domain" -> $cppSerializationDomain,
      "symbols" -> $cppSerializationSymbols|>, spec];
  If[!AssociationQ[field] ||
      Sort[Keys[field]] =!= Sort[{"domain", "symbols"}],
    err["E5", cs, <|"Serialization" -> field,
      "Detail" -> "SCC coupling serialization must contain exactly domain and symbols"|>]];
  domain = field["domain"];
  symbols = field["symbols"];
  If[!MemberQ[{"acb", "rational", "symbolic"}, domain] ||
      !ListQ[symbols] || !AllTrue[symbols,
        MatchQ[#, _Symbol] && !NumericQ[#] &],
    err["E5", cs, <|"Serialization" -> field,
      "Detail" -> "SCC coupling serialization field is malformed"|>]];
  If[Length[DeleteDuplicates[symbols, SameQ]] =!= Length[symbols],
    err["E5", cs, <|"Symbols" -> (sccSymbolIdentity /@ symbols),
      "Detail" -> "SCC coupling serialization contains duplicate exact symbols"|>]];
  names = SymbolName /@ symbols;
  If[Length[DeleteDuplicates[names]] =!= Length[names],
    err["E5", cs, <|"Symbols" -> (sccSymbolIdentity /@ symbols),
      "Detail" -> "native bare regulator names are duplicated across symbol contexts"|>]];
  If[!AllTrue[symbols, Context[#] === "Global`" &],
    err["E5", cs, <|"Symbols" -> (sccSymbolIdentity /@ symbols),
      "Detail" -> "native bare-name symbolic serialization cannot preserve non-Global regulator contexts"|>]];
  binders = SymbolName /@ {eps, t};
  If[Intersection[names, binders] =!= {},
    err["E5", cs, <|"Symbols" -> (sccSymbolIdentity /@ symbols),
      "Binders" -> (sccSymbolIdentity /@ {eps, t}),
      "Detail" -> "native regulator name collides with the epsilon or chart-variable binder"|>]];
  If[(domain === "symbolic") =!= (symbols =!= {}),
    err["E5", cs, <|"Serialization" -> field,
      "Detail" -> "symbolic serialization requires a nonempty symbol list and numeric/rational serialization requires none"|>]];
  <|"domain" -> domain, "symbols" -> symbols,
    "symbol_identities" -> (sccSymbolIdentity /@ symbols)|>];

cppPreparedRationalMultiplierJSON[prepared_Association,
    inputDigits_Integer, cs_Association, includeAnalyticTail_:False] := Module[
  {encoded = <|
      (* epsilon_shift is deliberately signed.  scc.prepare retains it and the
         later execution work contract proves that work_min supplies its halo. *)
      "epsilon_shift" -> prepared["EpsilonShift"],
      "center_pole_order" -> prepared["CenterPoleOrder"],
      "exact_identity" -> prepared["ExactIdentity"],
      "proven_zero" -> TrueQ[prepared["ProvenZero"]]|>, rationals,
   kernels = Lookup[prepared, "TaylorKernels", None], kernelCount},
  If[ListQ[kernels],
    encoded = Append[encoded, "kernels" ->
      Map[cppScalar[#, inputDigits, cs] &, kernels, {2}]]];
  If[TrueQ[includeAnalyticTail],
    rationals = Lookup[prepared, "AnalyticRationals", None];
    kernelCount = Lookup[prepared, "EpsilonKernelCount",
      If[ListQ[kernels], Length[kernels], None]];
    If[!IntegerQ[kernelCount] || kernelCount < 1 ||
        !ListQ[rationals] || Length[rationals] =!= kernelCount ||
        !AllTrue[rationals, AssociationQ[#] &&
          Sort[Keys[#]] === Sort[{"NumeratorCoefficients",
            "DenominatorCoefficients"}] &&
          ListQ[# ["NumeratorCoefficients"]] &&
          # ["NumeratorCoefficients"] =!= {} &&
          ListQ[# ["DenominatorCoefficients"]] &&
          # ["DenominatorCoefficients"] =!= {} &],
      err["E5", cs, <|"PreparedMultiplier" -> prepared,
        "Detail" -> "native rational-row tail certification requires one nonempty analytic numerator/denominator pair per epsilon kernel"|>]];
    encoded = Append[encoded, "analytic_coefficients" -> Map[
      <|"numerator" -> (cppScalar[#, inputDigits, cs] & /@
            # ["NumeratorCoefficients"]),
        "denominator" -> (cppScalar[#, inputDigits, cs] & /@
          # ["DenominatorCoefficients"])|> &, rationals]]];
  If[!KeyExistsQ[encoded, "kernels"] &&
      !KeyExistsQ[encoded, "analytic_coefficients"],
    err["E5", cs, <|"PreparedMultiplier" -> prepared,
      "Detail" -> "native rational multiplier needs Taylor kernels or a compact analytic rational source"|>]];
  encoded];

(* Prepare the user-facing coefficient row in precisely the same finite
   multiplier representation used by native SCC propagation.  The retained
   local, rather than a requested outer window, owns the rectangle: after a
   match or epsilon shift its honest CompleteMax may be smaller. *)
PrepareNativeRationalRow[cs_Association, sourceShape_Association,
    cvec_List, physicalVar_Symbol, serialization_:Automatic] := Module[
  {d = Lookup[cs, "SystemSize", None], epsWindow, tWindow, dimension,
   shape, field, domain, symbols, inputDigits, t, localExpressions,
   prepared, active, encodedEntries, identityPayload, identity},
  epsWindow = Lookup[sourceShape, "EpsWindow", None];
  tWindow = Lookup[sourceShape, "TWindow", None];
  dimension = Lookup[sourceShape, "Dimension", None];
  If[!IntegerQ[d] || d < 1 || Length[cvec] =!= d ||
      !AssociationQ[Lookup[cs, "ChartMap", None]] ||
      !MatchQ[Lookup[cs, "ChartVar", None], _Symbol] ||
      !AssociationQ[epsWindow] ||
      Sort[Keys[epsWindow]] =!= Sort[{"Min", "CompleteMax"}] ||
      !IntegerQ[epsWindow["Min"]] ||
      !IntegerQ[epsWindow["CompleteMax"]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"] ||
      !AssociationQ[tWindow] ||
      Sort[Keys[tWindow]] =!= {"CompleteMax"} ||
      !IntegerQ[tWindow["CompleteMax"]] ||
      tWindow["CompleteMax"] < 0 || dimension =!= d,
    err["E8", cs, <|"SourceShape" -> sourceShape,
      "CoefficientCount" -> Length[cvec], "Dimension" -> d,
      "Detail" -> "native rational-row preparation requires the retained local's exact epsilon/Taylor windows and one coefficient per physical component"|>]];
  field = sccSerializationField[serialization, cs];
  domain = field["domain"];
  symbols = field["symbols"];
  If[!MemberQ[{"acb", "rational"}, domain] || symbols =!= {},
    err["E5", cs, <|"Serialization" -> field,
      "Detail" -> "retained rational-row application supports only specialized Acb or exact Rational sessions; unresolved analytic regulators must be specialized before transport"|>]];
  t = cs["ChartVar"];
  shape = <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"],
    "Prescriptions" -> Lookup[cs, "Prescriptions", {}],
    "EpsWindow" -> epsWindow, "TWindow" -> tWindow,
    "Dimension" -> d|>;
  localExpressions = Map[Cancel[Together[# /. physicalVar ->
      cs["Center"] + cs["ChartMap", "Scale"]*t]] &, cvec];
  prepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
      shape, #, t, "SerializationDomain" -> domain,
      "PrepareTaylorKernels" -> False] & /@
    localExpressions;
  active = Select[Range[d],
    !TrueQ[prepared[[#]]["ProvenZero"]] &];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*
    cfg["WorkingPrecision"];
  encodedEntries = Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    Map[<|"column" -> # - 1,
        "multiplier" -> cppPreparedRationalMultiplierJSON[
          prepared[[#]], inputDigits, cs, True]|> &, active]];
  identityPayload = <|
    "schema" -> "diffexp2-prepared-rational-local-row-identity-v1",
    "chart" -> <|
      "center" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
        cs["Center"], t],
      "scale" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
        cs["ChartMap", "Scale"], t]|>,
    "source_shape" -> <|"epsilon_min" -> epsWindow["Min"],
      "epsilon_complete_max" -> epsWindow["CompleteMax"],
      "taylor_complete_max" -> tWindow["CompleteMax"],
      "dimension" -> d|>,
    "physical_variable" -> sccSymbolIdentity[physicalVar],
    "serialization" -> <|"domain" -> domain, "symbols" -> {}|>,
    "entries" -> MapIndexed[<|
        "column" -> First[#2] - 1,
        "physical_exact_identity" ->
          DiffExp2`SectorSeries`ExactExpressionIdentity[#1, physicalVar],
        "local_exact_identity" -> prepared[[First[#2]]]["ExactIdentity"],
        "proven_zero" -> TrueQ[prepared[[First[#2]]]["ProvenZero"]],
        "epsilon_shift" -> prepared[[First[#2]]]["EpsilonShift"],
        "center_pole_order" ->
          prepared[[First[#2]]]["CenterPoleOrder"]|> &, cvec]|>;
  identity = ExportString[identityPayload, "RawJSON", "Compact" -> True];
  <|"schema" -> "diffexp2-prepared-rational-local-row-v1",
    "columns" -> d, "exact_identity" -> identity,
    "entries" -> encodedEntries|>];

(* THE canonical exact cell record shared by the future parent system/theta
   manifests and every sparse cross-edge binding.  ProvenZero is decided in
   the exact Wolfram field; an inexact or specialized numerical zero remains
   active. *)
sccExactCellRecord[entry_, var_Symbol] := Module[
  {canonical = Together[entry]},
  <|"exact" ->
      DiffExp2`SectorSeries`ExactExpressionIdentity[canonical, var],
    "proven_zero" -> sccStructuralZeroQ[canonical]|>];

sccExactMatrixRecord[matrix_, var_Symbol, cs_Association] := Module[
  {dense = Normal[matrix], dims, cache = <||>, record},
  dims = Dimensions[dense];
  If[Length[dims] =!= 2 || dims[[1]] =!= dims[[2]],
    err["E6", cs, <|"Dimensions" -> dims,
      "Detail" -> "parent exact matrix record requires one square D by D matrix"|>]];
  record[entry_] := Module[
    {canonical = Together[entry], signature, key, cached, cell},
    signature = {canonical, sccSymbolIdentity[var]};
    key = Hash[signature, "SHA256"];
    cached = Lookup[cache, key, None];
    If[AssociationQ[cached] &&
        SameQ[cached["Signature"], signature],
      Return[cached["Cell"], Module]];
    If[cached =!= None,
      err["E6", cs, <|"CacheKey" -> key,
        "Detail" -> "parent exact-cell cache key collided with an unequal full signature"|>]];
    cell = sccExactCellRecord[canonical, var];
    AssociateTo[cache, key -> <|"Signature" -> signature,
      "Cell" -> cell|>];
    cell];
  Map[record, dense, {2}]];

(* Built exactly once for the top-level scc.prepare manifest.  Coupling
   groups carry only their indexed cells; duplicating both D by D records in
   every group would weaken ownership and bloat the handoff. *)
sccParentExactRecords[cs_Association] := Module[
  {clearKey = Lookup[cs, "SystemClearKey", None], parentInput,
   originalMatrix, originalVariable, thetaMatrix, thetaVariable,
   d = Lookup[cs, "SystemSize", None]},
  If[clearKey === None || !KeyExistsQ[$systemClearRegistry, clearKey],
    err["E6", cs, <|"SystemClearKey" -> clearKey,
      "Detail" -> "parent exact records require the registered original physical matrix"|>]];
  parentInput = $systemClearRegistry[clearKey];
  originalMatrix = Lookup[parentInput, "Matrix", None];
  originalVariable = Lookup[parentInput, "Variable", None];
  thetaMatrix = Lookup[cs, "ThetaOriginal", None];
  thetaVariable = Lookup[cs, "ChartVar", None];
  If[!MatchQ[originalVariable, _Symbol] ||
      !MatchQ[thetaVariable, _Symbol] ||
      Dimensions[originalMatrix] =!= {d, d} ||
      Dimensions[thetaMatrix] =!= {d, d},
    err["E6", cs, <|"SystemDimensions" -> Dimensions[originalMatrix],
      "ThetaDimensions" -> Dimensions[thetaMatrix], "Expected" -> {d, d},
      "Detail" -> "parent exact records do not bind two full D by D matrices and explicit variables"|>]];
  <|"exact_system_record" ->
      sccExactMatrixRecord[originalMatrix, originalVariable, cs],
    "exact_theta_record" ->
      sccExactMatrixRecord[thetaMatrix, thetaVariable, cs]|>];

sccPrescriptionIdentity[records_List, var_Symbol] := Map[
  Function[record, <|
    "factor" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
      Lookup[record, "ExactFactor", Lookup[record, "Factor", None]], var],
    "sign" -> Lookup[record, "Sign", None],
    "multiplicity" -> Lookup[record, "Multiplicity", None],
    "leading_coefficient_sign" ->
      Lookup[record, "LeadingCoeffSign", None]|>], records];

PrepareSCCCouplingMatrix[cs_Association, sourceBlock_Integer,
    targetBlock_Integer, sourceShape_Association,
    serialization_:Automatic] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], components, nb,
   sourceVertices, targetVertices, sourceDimension, checkedShape,
   preparationShape, matrix, expectedEdges, rawEntries, actualEdges,
   serializationField, symbols, domain, inputDigits, encodedEntries,
   matrixIdentity, matrixIdentityPayload,
   clearKey = Lookup[cs, "SystemClearKey", None], parentInput,
   originalMatrix, originalVariable, shapeCenter, shapeScale, shapeMap,
   requiredShapeKeys, missingShapeKeys, epsWindow, tWindow, sourceNCols,
   parentDimension = Lookup[cs, "SystemSize", None],
   t = Lookup[cs, "ChartVar", None]},
  If[!AssociationQ[seq] || !MatchQ[t, _Symbol] ||
      !MatrixQ[Lookup[cs, "ThetaOriginal", None]] ||
      !AssociationQ[Lookup[cs, "ChartMap", None]] ||
      !AllTrue[{"Center", "Radius", "Prescriptions"},
        KeyExistsQ[cs, #] &] || !KeyExistsQ[cs["ChartMap"], "Scale"],
    err["E6", cs, <|"Detail" ->
      "SCC coupling preparation requires an exact SCC chart envelope"|>]];
  serializationField = sccSerializationField[serialization, cs];
  domain = serializationField["domain"];
  symbols = serializationField["symbols"];
  components = seq["Components"];
  nb = Length[components];
  If[!Between[sourceBlock, {1, nb}] || !Between[targetBlock, {1, nb}] ||
      sourceBlock === targetBlock ||
      !MemberQ[seq["CondensationEdges"], {sourceBlock, targetBlock}],
    err["E8", cs, <|"SourceBlock" -> sourceBlock,
      "TargetBlock" -> targetBlock,
      "CondensationEdges" -> seq["CondensationEdges"],
      "Detail" -> "requested SCC coupling is not a certified cross-block condensation edge"|>]];
  sourceVertices = components[[sourceBlock]];
  targetVertices = components[[targetBlock]];
  If[clearKey === None || !KeyExistsQ[$systemClearRegistry, clearKey],
    err["E6", cs, <|"SystemClearKey" -> clearKey,
      "Detail" -> "SCC coupling preparation cannot bind the registered original physical matrix"|>]];
  parentInput = $systemClearRegistry[clearKey];
  originalMatrix = Lookup[parentInput, "Matrix", None];
  originalVariable = Lookup[parentInput, "Variable", None];
  If[!MatchQ[originalVariable, _Symbol] ||
      Dimensions[originalMatrix] =!= {parentDimension, parentDimension},
    err["E6", cs, <|"OriginalVariable" -> originalVariable,
      "OriginalDimensions" -> Dimensions[originalMatrix],
      "Expected" -> {parentDimension, parentDimension},
      "Detail" -> "registered original SCC parent matrix/variable is malformed"|>]];
  requiredShapeKeys = {"Center", "ChartMap", "Radius",
    "Prescriptions", "EpsWindow", "TWindow"};
  missingShapeKeys = Select[requiredShapeKeys,
    !KeyExistsQ[sourceShape, #] &];
  If[missingShapeKeys =!= {},
    err["E8", cs, <|"MissingKeys" -> missingShapeKeys,
      "Detail" -> "SCC coupling source shape must carry full chart geometry, prescriptions, and windows"|>]];
  shapeMap = sourceShape["ChartMap"];
  epsWindow = sourceShape["EpsWindow"];
  tWindow = sourceShape["TWindow"];
  If[!AssociationQ[shapeMap] || !KeyExistsQ[shapeMap, "Scale"] ||
      !AssociationQ[epsWindow] ||
      !AllTrue[{"Min", "CompleteMax"}, KeyExistsQ[epsWindow, #] &] ||
      !IntegerQ[epsWindow["Min"]] ||
      !IntegerQ[epsWindow["CompleteMax"]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"] ||
      !AssociationQ[tWindow] ||
      !KeyExistsQ[tWindow, "CompleteMax"] ||
      !IntegerQ[tWindow["CompleteMax"]] ||
      tWindow["CompleteMax"] < 0 ||
      !ListQ[sourceShape["Prescriptions"]],
    err["E8", cs, <|"SourceGeometry" ->
        KeyTake[sourceShape, requiredShapeKeys],
      "Detail" -> "SCC coupling source geometry/window contract is malformed"|>]];
  shapeCenter = sourceShape["Center"];
  shapeScale = shapeMap["Scale"];
  If[!SameQ[sourceShape["Radius"], cs["Radius"]],
    err["E8", cs, <|"SourceRadius" -> sourceShape["Radius"],
      "ChartRadius" -> cs["Radius"],
      "Detail" -> "SCC coupling source shape belongs to a different chart radius"|>]];
  If[!SameQ[shapeCenter, cs["Center"]],
    err["E8", cs, <|"SourceCenter" -> shapeCenter,
      "ChartCenter" -> cs["Center"],
      "Detail" -> "SCC coupling source shape belongs to a different exact chart center"|>]];
  If[KeyExistsQ[shapeMap, "Center"] &&
      !SameQ[shapeMap["Center"], shapeCenter],
    err["E8", cs, <|"SourceCenter" -> shapeCenter,
      "ChartMapCenter" -> shapeMap["Center"],
      "Detail" -> "SCC coupling source top-level and ChartMap centers disagree"|>]];
  If[!SameQ[shapeScale, cs["ChartMap", "Scale"]],
    err["E8", cs, <|"SourceScale" -> shapeScale,
      "ChartScale" -> cs["ChartMap", "Scale"],
      "Detail" -> "SCC coupling source shape belongs to a different exact chart scale"|>]];
  If[!SameQ[sourceShape["Prescriptions"], cs["Prescriptions"]],
    err["E8", cs, <|"SourcePrescriptions" -> sourceShape["Prescriptions"],
      "ChartPrescriptions" -> cs["Prescriptions"],
      "Detail" -> "SCC coupling source shape has different analytic-continuation prescriptions"|>]];
  If[KeyExistsQ[sourceShape, "Sectors"],
    checkedShape = DiffExp2`SectorSeries`ValidateLocalSolution[sourceShape];
    sourceDimension = Dimensions[
      First[checkedShape["Sectors"]]["Coeffs"]][[3]];
    sourceNCols = Dimensions[
      First[checkedShape["Sectors"]]["Coeffs"]][[2]];
    If[tWindow["CompleteMax"] =!= sourceNCols - 1 ||
        (KeyExistsQ[sourceShape, "Dimension"] &&
          sourceShape["Dimension"] =!= sourceDimension),
      err["E8", cs, <|"SourceDimension" -> sourceDimension,
        "DeclaredDimension" -> Lookup[sourceShape, "Dimension", None],
        "TaylorColumns" -> sourceNCols, "TWindow" -> tWindow,
        "Detail" -> "SCC coupling LocalSolution has a contradictory dimension or Taylor window"|>]];
    preparationShape = <|"Center" -> checkedShape["Center"],
      "ChartMap" -> checkedShape["ChartMap"],
      "Radius" -> checkedShape["Radius"],
      "Prescriptions" -> checkedShape["Prescriptions"],
      "EpsWindow" -> checkedShape["EpsWindow"],
      "TWindow" -> checkedShape["TWindow"],
      "Dimension" -> sourceDimension|>,
    sourceDimension = Lookup[sourceShape, "Dimension", None];
    If[!IntegerQ[sourceDimension] || sourceDimension < 1,
      err["E8", cs, <|"SourceDimension" -> sourceDimension,
        "Detail" -> "SCC coupling metadata shape requires an explicit positive integer Dimension"|>]];
    preparationShape = sourceShape];
  If[sourceDimension =!= Length[sourceVertices],
    err["E8", cs, <|"SourceDimension" -> sourceDimension,
      "ExpectedDimension" -> Length[sourceVertices],
      "SourceVertices" -> sourceVertices,
      "Detail" -> "SCC coupling source shape has the wrong component dimension"|>]];
  matrix = Map[Cancel[Together[#]] &,
    cs["ThetaOriginal"][[targetVertices, sourceVertices]], {2}];
  expectedEdges = Sort[Select[seq["DependencyEdges"],
    MemberQ[sourceVertices, #[[1]]] &&
      MemberQ[targetVertices, #[[2]]] &]];
  rawEntries = Flatten[Table[
    Module[{coefficient = matrix[[row, column]], originalCoefficient,
        prepared, edge, originalCell, thetaCell},
      If[sccStructuralZeroQ[coefficient], {},
        edge = {sourceVertices[[column]], targetVertices[[row]]};
        originalCoefficient = originalMatrix[[
          targetVertices[[row]], sourceVertices[[column]]]];
        originalCell = sccExactCellRecord[
          originalCoefficient, originalVariable];
        thetaCell = sccExactCellRecord[coefficient, t];
        If[TrueQ[originalCell["proven_zero"]] ||
            TrueQ[thetaCell["proven_zero"]],
          err["E6", cs, <|"GlobalEdge" -> edge,
            "OriginalCell" -> originalCell, "ThetaCell" -> thetaCell,
            "Detail" -> "active SCC coupling edge is marked proven zero in its exact parent cell record"|>]];
        prepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
          preparationShape, coefficient, t,
          "SerializationDomain" -> domain,
          "PrepareTaylorKernels" -> False];
        If[TrueQ[prepared["ProvenZero"]] ||
            prepared["ExactIdentity"] =!= thetaCell["exact"],
          err["E6", cs, <|"GlobalEdge" -> edge,
            "PreparedIdentity" -> prepared["ExactIdentity"],
            "PreparedProvenZero" -> prepared["ProvenZero"],
            "ThetaCell" -> thetaCell,
            "Detail" -> "active prepared multiplier does not equal its nonzero exact ThetaOriginal cell"|>]];
        {<|"GlobalEdge" -> edge, "Row" -> row - 1,
          "Column" -> column - 1,
          "SourceVertex" -> sourceVertices[[column]] - 1,
          "TargetVertex" -> targetVertices[[row]] - 1,
          "ExactOriginalEntry" -> originalCell["exact"],
          "ExactThetaEntry" -> thetaCell["exact"],
          "Prepared" -> prepared|>}]],
    {row, Length[targetVertices]}, {column, Length[sourceVertices]}], 2];
  actualEdges = Sort[Lookup[rawEntries, "GlobalEdge", {}]];
  If[actualEdges =!= expectedEdges,
    err["E6", cs, <|"SourceBlock" -> sourceBlock,
      "TargetBlock" -> targetBlock, "ExpectedEdges" -> expectedEdges,
      "PreparedEdges" -> actualEdges,
      "Detail" -> "prepared coupling nonzeros do not bind one-to-one to the exact SCC dependency graph"|>]];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*
    cfg["WorkingPrecision"];
  encodedEntries = Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    Map[<|"row" -> #["Row"], "column" -> #["Column"],
        "source_vertex" -> #["SourceVertex"],
        "target_vertex" -> #["TargetVertex"],
        "exact_original_entry" -> #["ExactOriginalEntry"],
        "exact_theta_entry" -> #["ExactThetaEntry"],
        "multiplier" -> cppPreparedRationalMultiplierJSON[
          #["Prepared"], inputDigits, cs, True]|> &, rawEntries]];
  (* A compact JSON string keeps the native field scalar while making its
     collision certificate fully structural.  No matrix InputForm string is
     trusted: every retained entry is the same context-explicit AST identity
     used in the full parent D by D records. *)
  matrixIdentityPayload = <|
    "schema" -> "diffexp2-scc-coupling-v1",
    "source_block" -> sourceBlock - 1,
    "target_block" -> targetBlock - 1,
    "source_vertices" -> (sourceVertices - 1),
    "target_vertices" -> (targetVertices - 1),
    "rows" -> Length[targetVertices],
    "columns" -> Length[sourceVertices],
    "source_shape" -> <|
      "center" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
        shapeCenter, t],
      "scale" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
        shapeScale, t],
      "radius" -> DiffExp2`SectorSeries`ExactExpressionIdentity[
        sourceShape["Radius"], t],
      "prescriptions" -> sccPrescriptionIdentity[
        sourceShape["Prescriptions"], t],
      "eps_min" -> epsWindow["Min"],
      "eps_complete_max" -> epsWindow["CompleteMax"],
      "t_complete_max" -> tWindow["CompleteMax"],
      "dimension" -> sourceDimension|>,
    "serialization" -> <|"domain" -> domain,
      "symbols" -> serializationField["symbol_identities"]|>,
    "entries" -> Map[<|"row" -> #["row"],
        "column" -> #["column"],
        "source_vertex" -> #["source_vertex"],
        "target_vertex" -> #["target_vertex"],
        "exact_original_entry" -> #["exact_original_entry"],
        "exact_theta_entry" -> #["exact_theta_entry"],
        "multiplier_exact_identity" -> #["multiplier", "exact_identity"],
        "epsilon_shift" -> #["multiplier", "epsilon_shift"],
        "center_pole_order" -> #["multiplier", "center_pole_order"],
        "proven_zero" -> #["multiplier", "proven_zero"]|> &,
      encodedEntries]|>;
  matrixIdentity = ExportString[matrixIdentityPayload,
    "RawJSON", "Compact" -> True];
  <|"source_block" -> sourceBlock - 1,
    "target_block" -> targetBlock - 1,
    "source_vertices" -> (sourceVertices - 1),
    "target_vertices" -> (targetVertices - 1),
    "rows" -> Length[targetVertices], "columns" -> Length[sourceVertices],
    "entries" -> encodedEntries, "exact_identity" -> matrixIdentity,
    "domain" -> domain, "symbols" -> (SymbolName /@ symbols)|>];

sccSourceZeroQ[None] := True;
sccSourceZeroQ[source_Association] :=
  Lookup[source, "Sectors", {}] === {} ||
    AllTrue[Flatten[Lookup[source, "Sectors", {}][[All, "Coeffs"]]],
      sccStructuralZeroQ];

sccReframeLocalSolution[ls_Association, cs_Association] := Join[ls, <|
  "Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
  "Radius" -> cs["Radius"], "Prescriptions" -> cs["Prescriptions"]|>];

sccSelectComponents[ls_Association, indices_List] := Module[{secs},
  secs = Map[Append[#, "Coeffs" -> # ["Coeffs"][[All, All, indices]]] &,
    ls["Sectors"]];
  Join[ls, <|"Sectors" -> secs|>]];

sccEmbedComponents[ls_Association, indices_List, dimension_Integer] := Module[
  {secs},
  secs = Map[Function[sec, Module[{arr = sec["Coeffs"], out},
      out = ConstantArray[0,
        {Length[arr], Length[arr[[1]]], dimension}];
      Do[out[[All, All, indices[[j]]]] = arr[[All, All, j]],
        {j, Length[indices]}];
      Append[sec, "Coeffs" -> out]]], ls["Sectors"]];
  Join[ls, <|"Sectors" -> secs|>]];

sccSourceToLocalSolution[cs_Association, source_Association,
    dimension_Integer] := Module[{secs = Lookup[source, "Sectors", {}],
   win = source["EpsWindow"], tw = source["TWindow"], dims},
  If[secs === {},
    Return[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"],
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> ConstantArray[0,
          {1, tw["CompleteMax"] + 1, dimension}]|>},
      "EpsWindow" -> <|"Min" -> win["CompleteMax"],
        "CompleteMax" -> win["CompleteMax"]|>,
      "TWindow" -> tw, "ErrorEstimate" -> {0},
      "Prescriptions" -> cs["Prescriptions"]|>, Module]];
  dims = Dimensions[First[secs]["Coeffs"]];
  If[Length[dims] =!= 3 || Last[dims] =!= dimension ||
      !AllTrue[secs, Dimensions[# ["Coeffs"]][[3]] === dimension &],
    err["E8", cs, <|"SourceDimensions" -> (Dimensions[# ["Coeffs"]] & /@ secs),
      "ExpectedComponents" -> dimension,
      "Detail" -> "SCC source component dimension mismatch"|>]];
  <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"], "Sectors" -> secs,
    "EpsWindow" -> win, "TWindow" -> tw,
    "ErrorEstimate" -> ConstantArray[0,
      win["CompleteMax"] - win["Min"] + 1],
    "Prescriptions" -> cs["Prescriptions"]|>];

sccLocalSolutionSource[ls_Association] := KeyTake[ls,
  {"Sectors", "EpsWindow", "TWindow"}];

sccProjectSource[source_Association, indices_List] := Module[{secs},
  If[Lookup[source, "Sectors", {}] === {}, Return[source, Module]];
  secs = Map[Append[#, "Coeffs" -> # ["Coeffs"][[All, All, indices]]] &,
    source["Sectors"]];
  Join[source, <|"Sectors" -> If[
    AllTrue[Flatten[secs[[All, "Coeffs"]]], sccStructuralZeroQ], {}, secs]|>]];

sccMatrixTimesLocalSolution[cs_Association, matrix_?MatrixQ,
    ls_Association] := Module[
  {nr = Length[matrix], nc, inputComponents, terms = {}, t = cs["ChartVar"],
   coeff, one, product},
  nc = If[nr === 0, 0, Length[First[matrix]]];
  inputComponents = Dimensions[First[ls["Sectors"]]["Coeffs"]][[3]];
  If[nc =!= inputComponents,
    err["E8", cs, <|"MatrixDimensions" -> Dimensions[matrix],
      "InputComponents" -> inputComponents,
      "Detail" -> "SCC coupling matrix/component mismatch"|>]];
  Do[
    coeff = Cancel[Together[matrix[[r, c]]]];
    If[!sccStructuralZeroQ[coeff],
      one = sccSelectComponents[ls, {c}];
      product = DiffExp2`SectorSeries`MultiplyRational[one, coeff, t];
      AppendTo[terms, sccEmbedComponents[product, {r}, nr]]],
    {r, nr}, {c, nc}];
  If[terms === {}, None,
    If[Length[terms] === 1, First[terms],
      DiffExp2`SectorSeries`CombineLocalSolutions[
        ConstantArray[1, Length[terms]], terms]]]];

sccCombineSources[cs_Association, sources_List, dimension_Integer] := Module[
  {active, lss, combined},
  active = Select[sources, AssociationQ[#] && !sccSourceZeroQ[#] &];
  If[active === {}, Return[None, Module]];
  lss = sccSourceToLocalSolution[cs, #, dimension] & /@ active;
  combined = If[Length[lss] === 1, First[lss],
    DiffExp2`SectorSeries`CombineLocalSolutions[
      ConstantArray[1, Length[lss]], lss]];
  sccLocalSolutionSource[combined]];

(* ---- public functions ---- *)

(* SolveHomogeneous memo: the fundamental basis depends ONLY on the
   ChartSystem and the request window — solveHomogeneousCore reads cs
   (blocks/families, the prepared cleared system, V/gauge) and req (TOrder,
   EpsWindow); boundary values never enter (they meet the basis later in
   MatchWeights).  The lo/hi endpoint transports of one level share the
   anchor chart (and any coinciding centers), so the second direction's
   solve is a pure replay: memoize on (Hash[cs], request window, WP).
   PrepareChart already dedups cs itself ($pcCache), so shared charts hash
   identically.  WorkingPrecision is in the key because ratEpsList
   numericizes large coefficients at 2x WP.  Memory policy: the cache holds
   ONE system's charts (flushed when SystemHash changes — each FT level is
   a fresh system), is entry-capped as a runaway guard, and is cleared by
   ClearSolveCaches[] from API`LoadSystem.  A memo hit skips the
   ODEResidualCheck rerun: the identical result already passed it when
   computed. *)
$shCache = <||>; $shSysTag = None; $shCacheMax = 64;
HomogeneousCacheCapacity[] := $shCacheMax;
homogeneousCacheKey[cs_Association, req_Association] :=
  {Hash[cs], req["TOrder"], req["EpsWindow", "Min"],
    req["EpsWindow", "CompleteMax"], cfg["WorkingPrecision"],
    cfg["ChopPrecision"],
    cfg["RecurrenceBackend"],
    TrueQ[$cppExactDomain],
    TrueQ[$disableAdaptiveLowerFrames],
    TrueQ[$disableRationalDenominatorFusion],
    TrueQ[$disableGroupedSpectralTransform],
    TrueQ[$disablePolynomialNhatTransform]};
SolveHomogeneous[cs_Association, req_Association] := Module[
  {tag = Lookup[cs, "SolveCacheTag", Lookup[cs, "SystemHash", None]], key},
  key = homogeneousCacheKey[cs, req];
  If[tag =!= $shSysTag, $shCache = <||>; $shSysTag = tag];
  If[KeyExistsQ[$shCache, key], Return[$shCache[key]]];
  If[Length[$shCache] >= $shCacheMax, $shCache = <||>];
  $shCache[key] = If[TrueQ[Lookup[cs, "SCCSkeleton", False]] &&
      AssociationQ[Lookup[cs, "IntegrationSequence", None]] &&
      Length[cs["IntegrationSequence", "Components"]] > 1,
    (* Diagonal block bases are implementation scratch for this assembled
       parent basis.  Isolate their memo entries dynamically so a six-block
       chart still occupies exactly one persistent cache slot and cannot
       evict an already-prewarmed arm. *)
    Block[{$shCache = <||>, $shSysTag = tag},
      Module[{blocks, execution, fullcs, fs},
        blocks = sccBlockChartSystem[cs, #] & /@
          Range[Length[cs["IntegrationSequence", "Components"]]];
        execution = sccExecutionPlan[cs, req, blocks];
        If[execution["Mode"] === "MonolithicStrict",
          fullcs = sccTryMonolithicChartSystem[cs];
          If[FailureQ[fullcs],
            fs = solveSCCBasis[cs, req, blocks];
            sccAnnotateDeclinedBasis[fs, cs, execution, fullcs],
            fs = solveHomogeneousCore[fullcs, req];
            sccAnnotateCoarsenedBasis[fs, cs, execution, fullcs]],
          solveSCCBasis[cs, req, blocks]]]],
    solveHomogeneousCore[cs, req]]];

PrewarmHomogeneousBatch[chartSystems_List, req_Association] := Module[
  {systems, tag, keys, uncached, captures, requests, lengths, response,
   results, offsets, assembled, threads, started = SessionTime[], allSystems,
   sccSystems, ordinaryInput},
  If[chartSystems === {}, Return[{}, Module]];
  If[!AllTrue[chartSystems, AssociationQ],
    err["E6", <|"Center" -> "homogeneous-batch"|>, <|
      "Detail" -> "homogeneous batch entries must be prepared ChartSystem associations"|>]];
  If[cfg["RecurrenceBackend"] =!= "Cpp",
    err["E6", First[chartSystems], <|
      "Detail" -> "PrewarmHomogeneousBatch requires RecurrenceBackend -> Cpp; no alternate solver is selected"|>]];
  If[TrueQ[$disableGroupedSpectralTransform],
    err["E6", First[chartSystems], <|
      "Detail" -> "homogeneous native batching requires the production grouped spectral-transform path"|>]];
  allSystems = DeleteDuplicatesBy[chartSystems,
    homogeneousCacheKey[#, req] &];
  sccSystems = Select[allSystems,
    TrueQ[Lookup[#, "SCCSkeleton", False]] &];
  (* A skeleton expands into several independently prepared native requests;
     the existing capture/injection protocol is one request-pool per original
     ChartSystem.  Prewarm those systems through their ordinary strict SCC
     solve (still cached), then batch only genuine one-frame systems. *)
  If[sccSystems =!= {},
    Block[{$cppHomogeneousBatchCapture = False},
      SolveHomogeneous[#, req] & /@ sccSystems]];
  ordinaryInput = Select[allSystems,
    !TrueQ[Lookup[#, "SCCSkeleton", False]] &];
  If[ordinaryInput === {},
    Return[SolveHomogeneous[#, req] & /@ allSystems, Module]];
  tag = Lookup[First[ordinaryInput], "SolveCacheTag",
    Lookup[First[ordinaryInput], "SystemHash", None]];
  If[!AllTrue[ordinaryInput,
      Lookup[#, "SolveCacheTag", Lookup[#, "SystemHash", None]] === tag &],
    err["E6", First[ordinaryInput], <|
      "Detail" -> "one homogeneous native batch cannot mix different differential systems"|>]];
  If[tag =!= $shSysTag, $shCache = <||>; $shSysTag = tag];
  (* SolveHomogeneous's key is the authoritative identity.  This also drops
     the shared anchor chart appearing in both endpoint plans. *)
  systems = ordinaryInput;
  keys = homogeneousCacheKey[#, req] & /@ systems;
  uncached = Pick[systems, Not /@ (KeyExistsQ[$shCache, #] & /@ keys)];
  If[uncached === {}, Return[SolveHomogeneous[#, req] & /@ systems, Module]];
  If[Length[$shCache] + Length[uncached] > $shCacheMax,
    err["E6", First[uncached], <|
      "Charts" -> Length[uncached], "CacheLimit" -> $shCacheMax,
      "Detail" -> "homogeneous batch would exceed the bounded chart cache; submit smaller waves"|>]];
  captures = Map[Function[cs,
    Module[{captured = Catch[
        Block[{$cppHomogeneousBatchCapture = True},
          SolveHomogeneous[cs, req]], $cppHomogeneousBatchTag]},
      If[!AssociationQ[captured] ||
          !ListQ[Lookup[captured, "Requests", None]] ||
          !AssociationQ[Lookup[captured, "Metadata", None]],
        err["E5", cs, <|
          "Detail" -> "homogeneous solve did not yield a capturable C++ recurrence batch"|>]];
      captured]], uncached];
  lengths = Length[# ["Requests"]] & /@ captures;
  If[AnyTrue[lengths, # < 1 &],
    err["E5", First[uncached], <|
      "Detail" -> "captured homogeneous recurrence batch was empty"|>]];
  requests = Flatten[# ["Requests"] & /@ captures, 1];
  threads = cppConfiguredThreads[Length[requests]];
  response = If[TrueQ[$cppUsePersistentSessions] &&
      Environment["DE2_CPP_PERSISTENT"] =!= "0",
    DiffExp2`CppBackend`RunPersistentRequestGroups[captures, threads],
    DiffExp2`CppBackend`RunRequest[
      <|"batch" -> requests, "threads" -> threads|>]];
  If[FailureQ[response] || Lookup[response, "status", "error"] =!= "ok" ||
      !ListQ[Lookup[response, "results", None]] ||
      Length[response["results"]] =!= Length[requests],
    err["E5", First[uncached], <|"BackendResponse" -> response,
      "Detail" -> "combined homogeneous recurrence batch failed"|>]];
  results = response["results"];
  offsets = FoldList[Plus, 0, lengths];
  assembled = MapIndexed[Function[{cs, idx}, Module[
      {i = First[idx], injection, uses, solved},
      injection = <|"Requests" -> captures[[i, "Requests"]],
        "Results" -> Take[results, {offsets[[i]] + 1, offsets[[i + 1]]}]|>;
      uses = 0;
      solved = Block[{$cppHomogeneousBatchInjection = injection,
          $cppHomogeneousBatchInjectionUses = 0},
        Module[{value = SolveHomogeneous[cs, req]},
          uses = $cppHomogeneousBatchInjectionUses; value]];
      If[uses =!= 1,
        err["E5", cs, <|"InjectionUses" -> uses,
          "Detail" -> "prewarmed recurrence response was not consumed exactly once"|>]];
      solved]], uncached];
  If[Environment["DEBUG_CPP_RECURRENCE"] === "1",
    Print["CPPARM BATCH ", <|"Charts" -> Length[uncached],
      "Tasks" -> Length[requests], "Threads" -> threads,
      "Seconds" -> N[SessionTime[] - started, 6]|>]];
  (* Return in the caller's deduplicated order; cached entries and the newly
     verified entries are indistinguishable at this point. *)
  SolveHomogeneous[#, req] & /@ allSystems];

ClearSolveCaches[] := ($pcCache = <||>; $shCache = <||>; $shSysTag = None;
  $systemClearRegistry = <||>; $globalClearedCache = <||>;
  $chartClearedCache = <||>; $clearedSymbolicLegacyCache = <||>;
  $exactSCCStructureCache = <||>;
  $physicalClearedODECache = <||>; $cppStaticOperatorCache = <||>;
  $nativeSCCCompositeCache = <||>;
  $nativeSCCCompositeReservedCapacity = 0;
  $homogeneousFramePlanOverride = None;
  $cppHomogeneousFrameOverride = None;
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
  DiffExp2`CppBackend`ClearPersistentSessions[];);

DropWolframPreparationCaches[] := Module[{},
  $pcCache = <||>; $shCache = <||>; $shSysTag = None;
  $systemClearRegistry = <||>; $globalClearedCache = <||>;
  $chartClearedCache = <||>; $clearedSymbolicLegacyCache = <||>;
  $exactSCCStructureCache = <||>;
  $physicalClearedODECache = <||>; $cppStaticOperatorCache = <||>;
  $nativeSCCCompositeCache = <||>;
  $nativeSCCCompositeReservedCapacity = 0;
  $homogeneousFramePlanOverride = None;
  $cppHomogeneousFrameOverride = None;
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
  Null];

(* The exact homogeneous work rectangle is needed both by the ordinary
   column builder and by persistent SCC preparation.  Keeping the arithmetic
   in one helper lets a composite take the union of heterogeneous Jordan
   block rectangles before any native operator is serialized.  A temporary
   override passes the already-computed plan into the capture call, avoiding
   a second denominator clear/pole-depth scan. *)
$homogeneousFramePlanOverride = None;
$cppHomogeneousFrameOverride = None;

homogeneousFramePlan[cs_Association, req_Association] := Module[
  {blocks = blockList[cs], fams = cs["Families"], nmax, reqMin, reqMax,
   pMax, cdMax, symbolic, poleDepth, singleUseDepth, spectralDepth,
   transformDepth, wideTop, terminalFb, startFb, adaptiveQ,
   epsilonRegular, epsilonRegularQ, initialSourceDepth},
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"];
  reqMax = req["EpsWindow", "CompleteMax"];
  pMax = Max[0, Max[Table[
      logCeiling[cs, b["a"], b["b"], b["q"] - 1], {b, blocks}]]];
  cdMax = Max[0, Max[#["CollisionDepth"] & /@ fams]];
  (* The hybrid principal solve is currently a native recurrence capability.
     Keep the Wolfram recurrence on its original spectral operator so it
     remains an independent oracle for parity tests. *)
  epsilonRegular = If[cfg["RecurrenceBackend"] === "Cpp",
    epsilonRegularPrincipalCertificate[cs, nmax],
    <|"Eligible" -> False,
      "Reason" -> "epsilon-regular-principal-is-native-only"|>];
  epsilonRegularQ = TrueQ[epsilonRegular["Eligible"]];
  symbolic = If[epsilonRegularQ, epsilonRegular["Symbolic"],
    clearedSymbolic[cs]];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  singleUseDepth = recurrenceSingleUsePoleDepth[symbolic];
  initialSourceDepth = If[epsilonRegularQ,
    epsilonRegular["InitialSourcePoleDepth"], 0];
  spectralDepth = spectralTransformPoleDepth[cs];
  transformDepth = finalTransformPoleDepth[cs, nmax];
  wideTop = reqMax + pMax + cdMax + 2 -
    Min[0, reqMin - pMax - 2];
  wideTop = Max[wideTop,
    reqMax + pMax + cdMax + poleDepth] + transformDepth;
  terminalFb = Min[Min[reqMin, 0] - pMax - cdMax - 2,
      reqMin - pMax - cdMax - poleDepth - initialSourceDepth] -
    spectralDepth;
  startFb = Min[Min[reqMin, 0] - pMax - cdMax - 2,
      reqMin - pMax - cdMax - singleUseDepth - initialSourceDepth] -
    spectralDepth;
  adaptiveQ = !TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    startFb > terminalFb && !TrueQ[$disableAdaptiveLowerFrames];
  If[!adaptiveQ, startFb = terminalFb];
  <|"Identity" -> {cs, req}, "PMax" -> pMax,
    "CollisionDepthMax" -> cdMax, "Symbolic" -> symbolic,
    "PoleDepth" -> poleDepth, "SingleUseDepth" -> singleUseDepth,
    "EpsilonRegularPrincipal" -> epsilonRegularQ,
    "EpsilonRegularCertificate" -> KeyDrop[epsilonRegular, "Symbolic"],
    "InitialSourceDepth" -> initialSourceDepth,
    "SpectralDepth" -> spectralDepth,
    "TransformDepth" -> transformDepth, "FrameTop" -> wideTop,
    "TerminalFrameBase" -> terminalFb, "StartFrameBase" -> startFb,
    "Adaptive" -> adaptiveQ|>];

homogeneousFramePlanFor[cs_Association, req_Association] :=
  If[AssociationQ[$homogeneousFramePlanOverride] &&
      SameQ[Lookup[$homogeneousFramePlanOverride, "Identity", None],
        {cs, req}],
    $homogeneousFramePlanOverride,
    homogeneousFramePlan[cs, req]];

makeRegularOwnerPhysicalPreparation[cs_Association, req_Association,
    equationIdentity_String, physicalData_Association] := Module[
  {core},
  (* Validate before retention; the identity below then binds the complete
     immutable hint, not merely its compact physical-payload digest. *)
  regularClearedSymbolicFromPhysicalData[cs, physicalData];
  core = <|
    "Schema" -> "DiffExp2RegularOwnerPhysicalPreparation/v1",
    "EquationIdentity" -> equationIdentity,
    "PhysicalPayloadIdentity" -> physicalData["Identity"],
    "Geometry" -> cppPersistentGeometry[cs],
    "Request" -> req,
    "PhysicalData" -> physicalData|>;
  Append[core, "Identity" -> ("de2-regular-owner-preparation-" <>
    IntegerString[Hash[core, "SHA256"], 16, 64])]];

regularOwnerPhysicalPreparationData[cs_Association, req_Association,
    preparation_Association] := Module[{core, data, expectedIdentity},
  core = KeyDrop[preparation, "Identity"];
  data = Lookup[preparation, "PhysicalData", None];
  expectedIdentity = "de2-regular-owner-preparation-" <>
    IntegerString[Hash[core, "SHA256"], 16, 64];
  If[Sort[Keys[preparation]] =!= Sort[{
        "Schema", "EquationIdentity", "PhysicalPayloadIdentity",
        "Geometry", "Request", "PhysicalData", "Identity"}] ||
      Lookup[preparation, "Schema", None] =!=
        "DiffExp2RegularOwnerPhysicalPreparation/v1" ||
      !StringQ[Lookup[preparation, "EquationIdentity", None]] ||
      Lookup[preparation, "Geometry", None] =!= cppPersistentGeometry[cs] ||
      Lookup[preparation, "Request", None] =!= req ||
      !AssociationQ[data] ||
      Lookup[preparation, "PhysicalPayloadIdentity", None] =!=
        Lookup[data, "Identity", None] ||
      Lookup[preparation, "Identity", None] =!= expectedIdentity,
    err["E6", cs, <|
      "OwnerPreparation" -> KeyDrop[preparation, "PhysicalData"],
      "Detail" ->
        "regular owner physical preparation failed its exact request, geometry, or content identity"|>]];
  regularClearedSymbolicFromPhysicalData[cs, data];
  data];

(* Prepare the complete immutable operator once.  Chart-only owner
   construction and every dynamic local run share this exact capture, so a
   later value or basis solve cannot silently drift to a different frame,
   physical q/C payload, or analytic chart identity. *)
prepareNativeLocalFamilyShared[cs_Association, req_Association,
    tag_Association, initPrototype_List,
    ownerPhysicalPreparation_:Automatic] := Module[
  {d = Lookup[cs, "SystemSize", 0], a, b, p, nmax, reqMin, reqMax,
   matchingFamilies, matchingBlocks, allowedP, blocks, fams, pMax,
   pBudget, cdMax, symbolic,
   poleDepth, spectralDepth, transformDepth, wideTop, fb, W, prep, vPrep,
   symbols, domain, wp, inputDigits, precisionBits, staticRecord,
   ownerPhysicalData,
   persistentMetadata, shared, ownerTimingQ, ownerPhaseStarted, ownerPhase},
  ownerTimingQ = Environment["DE2_NATIVE_OWNER_TIMING"] === "1";
  ownerPhaseStarted = SessionTime[];
  ownerPhase[label_String] := If[ownerTimingQ,
    Module[{now = SessionTime[]},
      Print["DE2 NATIVE OWNER PHASE center=", InputForm[cs["Center"]],
        " phase=", label, " elapsedMs=",
        N[1000 (now - ownerPhaseStarted), 8], " frame={", fb, ",",
        W, "} dimension=", d];
      ownerPhaseStarted = now]];
  If[cfg["RecurrenceBackend"] =!= "Cpp" ||
      !TrueQ[$cppUsePersistentSessions] ||
      Environment["DE2_CPP_PERSISTENT"] === "0",
    err["E5", cs, <|"Detail" ->
      "SolveNativeLocalFamily requires the persistent C++ recurrence backend; no fallback is permitted"|>]];
  If[TrueQ[Lookup[cs, "SCCSkeleton", False]],
    err["E6", cs, <|"Detail" ->
      "native local family solve does not yet orchestrate an SCC skeleton"|>]];
  If[d < 1 || !MatrixQ[Lookup[cs, "Gauge", None]] ||
      cs["Gauge"] =!= IdentityMatrix[d],
    err["E5", cs, <|"Detail" ->
      "native local family solve requires Gauge === IdentityMatrix[d]; a native sequential V-then-Gauge assembly chain is not implemented yet"|>]];
  If[TrueQ[$disableGroupedSpectralTransform],
    err["E5", cs, <|"Detail" ->
      "native local family solve requires the grouped native assembly path"|>]];
  If[!AllTrue[{"a", "b", "p"}, KeyExistsQ[tag, #] &],
    err["E8", cs, <|"Tag" -> tag, "Detail" ->
      "native local family tag requires exact a, b, and p fields"|>]];
  {a, b, p} = Lookup[tag, {"a", "b", "p"}];
  If[!IntegerQ[p] || p < 0,
    err["E8", cs, <|"Tag" -> tag, "Detail" ->
      "native local family p must be a nonnegative exact integer"|>]];
  If[!FreeQ[{a, b}, _?InexactNumberQ] ||
      !FreeQ[{a, b}, DiffExp2`Config`CanonicalEps[]] ||
      !AllTrue[{a, b}, NumericQ],
    err["E5", cs, <|"Tag" -> tag, "Detail" ->
      "native local family a and b must be exact, epsilon-independent, and fully specialized"|>]];
  If[!AllTrue[{"TOrder", "EpsWindow"}, KeyExistsQ[req, #] &] ||
      !AssociationQ[Lookup[req, "EpsWindow", None]] ||
      !AllTrue[{"Min", "CompleteMax"},
        KeyExistsQ[req["EpsWindow"], #] &],
    err["E8", cs, <|"Request" -> req, "Detail" ->
      "native local family request has no complete Taylor/epsilon window"|>]];
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"];
  reqMax = req["EpsWindow", "CompleteMax"];
  If[!IntegerQ[nmax] || nmax < 0 || !IntegerQ[reqMin] ||
      !IntegerQ[reqMax] || reqMin > reqMax,
    err["E8", cs, <|"Request" -> req, "Detail" ->
      "native local family request windows must be finite ordered integers"|>]];
  If[Length[initPrototype] =!= p + 1 || !AllTrue[initPrototype,
      ListQ[#] && Length[#] === d &&
        AllTrue[#, DiffExp2`EpsSeries`ESQ] &],
    err["E8", cs, <|"Tag" -> tag, "InitialDimensions" ->
      Quiet[Check[Dimensions[initPrototype], Missing["Ragged"]]],
      "Expected" -> {p + 1, d}, "Detail" ->
      "native local family initial ladder must contain p+1 rows of d EpsSeries values"|>]];

  fams = Lookup[cs, "Families", {}];
  matchingFamilies = Select[fams, AnyTrue[Lookup[#, "Roots", {}],
      zeroQ[a - # ["a"]] && zeroQ[b - # ["b"]] &] &];
  If[matchingFamilies === {},
    err["E8", cs, <|"Tag" -> tag, "Detail" ->
      "native local family tag is not an exact indicial root of this chart"|>]];
  (* local.solve currently assembles recurrence output but does not expose
     RecurrenceResult::hits for the Wolfram compensation transaction.  Be
     conservative across the whole prepared frame: a hit against any family
     must keep using the ordinary compensated path until compensation itself
     is session-owned. *)
  If[AnyTrue[fams, Lookup[#, "Collisions", {}] =!= {} &],
    err["E5", cs, <|"Tag" -> tag, "Detail" ->
      "native local handle cannot yet represent Wolfram pseudo-resonant family compensation anywhere in the prepared frame"|>]];

  blocks = blockList[cs];
  matchingBlocks = Select[blocks,
    zeroQ[a - # ["a"]] && zeroQ[b - # ["b"]] &];
  allowedP = DeleteDuplicates[Flatten[Table[
      logCeiling[cs, block["a"], block["b"], qpos],
      {block, matchingBlocks}, {qpos, 0, block["q"] - 1}]]];
  If[!MemberQ[allowedP, p],
    err["E8", cs, <|"Tag" -> tag, "AllowedLogPowers" -> allowedP,
      "Detail" ->
        "native local family p is not an exact log ceiling of the selected indicial root"|>]];
  pMax = Max[Join[{0}, Table[
      logCeiling[cs, block["a"], block["b"], block["q"] - 1],
      {block, blocks}]]];
  pBudget = Max[p, pMax];
  cdMax = Max[Join[{0}, Lookup[fams, "CollisionDepth", {}]]];
  ownerPhysicalData = If[AssociationQ[ownerPhysicalPreparation],
    regularOwnerPhysicalPreparationData[
      cs, req, ownerPhysicalPreparation], Automatic];
  symbolic = If[AssociationQ[ownerPhysicalData],
    regularClearedSymbolicFromPhysicalData[cs, ownerPhysicalData],
    clearedSymbolic[cs]];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  spectralDepth = spectralTransformPoleDepth[cs];
  transformDepth = finalTransformPoleDepth[cs, nmax];
  wideTop = reqMax + pBudget + cdMax + 2 -
    Min[0, reqMin - pBudget - 2];
  wideTop = Max[wideTop,
    reqMax + pBudget + cdMax + poleDepth] + transformDepth;
  fb = Min[Min[reqMin, 0] - pBudget - cdMax - 2,
      reqMin - pBudget - cdMax - poleDepth] - spectralDepth;
  W = wideTop - fb + 1;
  ownerPhase["preflightAndFramePlan"];
  prep = prepareCleared[cs, fb, W, symbolic];
  ownerPhase["prepareCleared"];
  vPrep = prepareFramedMatrix[cs["V"],
    DiffExp2`Config`CanonicalEps[], fb, W, cs];
  ownerPhase["prepareFramedMatrix"];
  symbols = cppRegulatorSymbols[cs, prep, a, b, p, nmax,
    None, initPrototype, vPrep];
  If[symbols =!= {},
    err["E5", cs, <|"Tag" -> tag,
      "RegulatorSymbols" -> (SymbolName /@ symbols), "Detail" ->
      "native local handle rejects unresolved analytic regulators; specialize them before solving"|>]];

  domain = If[TrueQ[$cppExactDomain], "rational", "acb"];
  wp = cfg["WorkingPrecision"];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  precisionBits = cppNativePrecisionBits[wp];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    staticRecord = cppStaticOperatorPayload[cs, prep, blocks, fb, W,
      vPrep, inputDigits, precisionBits, ownerPhysicalData];
    ownerPhase["buildAndSerializeStaticOperator"]];
  persistentMetadata = Append[cppPersistentMetadata[cs, fb, W],
    "PreparedToken" -> staticRecord["Token"]];
  ownerPhase["metadataAndIdentity"];
  shared = <|"ChartSystem" -> cs, "RequestSpec" -> req,
    "Dimension" -> d,
    "Tag" -> <|"a" -> a, "b" -> b, "p" -> p|>,
    "TaylorOrder" -> nmax,
    "RequestedMin" -> reqMin, "RequestedMax" -> reqMax,
    "FrameBase" -> fb, "FrameWidth" -> W,
    "PreparedCleared" -> prep, "PreparedAssembly" -> vPrep,
    "SerializationDomain" -> domain, "InputDigits" -> inputDigits,
    "StaticRecord" -> staticRecord,
    "PersistentMetadata" -> persistentMetadata|>;
  shared];

prepareNativeLocalFamilyDynamic[shared_Association, init_List,
    retainValueSolverPrototype_:False] := Module[
  {cs = shared["ChartSystem"], req = shared["RequestSpec"],
   d = shared["Dimension"], normalizedTag = shared["Tag"],
   a, b, p, nmax = shared["TaylorOrder"],
   reqMin = shared["RequestedMin"], reqMax = shared["RequestedMax"],
   fb = shared["FrameBase"], W = shared["FrameWidth"],
   prep = shared["PreparedCleared"], vPrep = shared["PreparedAssembly"],
   domain = shared["SerializationDomain"],
   inputDigits = shared["InputDigits"],
   staticRecord = shared["StaticRecord"],
   persistentMetadata = shared["PersistentMetadata"], symbols, request,
   checkpointIdentity, localMetadata, prepared, valueRunKeys},
  {a, b, p} = Lookup[normalizedTag, {"a", "b", "p"}];
  If[Length[init] =!= p + 1 || !AllTrue[init,
      ListQ[#] && Length[#] === d &&
        AllTrue[#, DiffExp2`EpsSeries`ESQ] &],
    err["E8", cs, <|"Tag" -> normalizedTag,
      "InitialDimensions" ->
        Quiet[Check[Dimensions[init], Missing["Ragged"]]],
      "Expected" -> {p + 1, d}, "Detail" ->
        "native local family initial ladder must contain p+1 rows of d EpsSeries values"|>]];
  (* The expensive chart-wide regulator audit was completed once while the
     immutable preparation was captured.  A reused preparation can vary only
     its initial ladder, so audit that small dynamic payload directly. *)
  symbols = DeleteDuplicates[Cases[init,
    s_Symbol /; Context[s] === "Global`" &&
      s =!= DiffExp2`Config`CanonicalEps[] &&
      s =!= Lookup[cs, "ChartVar", None] && !NumericQ[s], Infinity]];
  If[symbols =!= {},
    err["E5", cs, <|"Tag" -> normalizedTag,
      "RegulatorSymbols" -> (SymbolName /@ SortBy[symbols, SymbolName]),
      "Detail" ->
        "native local handle rejects unresolved analytic regulators; specialize them before solving"|>]];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    request = Block[{$cppStaticRecordOverride = staticRecord,
        $cppBuildRequestOnly = True},
      cppRunRecursionCore[cs, prep, a, b, p, nmax, None, fb, W,
        init, vPrep]];
    checkpointIdentity = "de2-native-local-" <>
      IntegerString[Hash[{persistentMetadata["SystemIdentity"],
        persistentMetadata["ChartIdentity"], cs["ChartMap"],
        cs["Radius"], cs["Prescriptions"], {a, b, p}, req, init,
        staticRecord["Token"]}, "SHA256"],
        16, 64];
    localMetadata = cppNativeLocalMetadata[cs, a, b, p, inputDigits,
      checkpointIdentity]];
  prepared = <|"Dimension" -> d, "Tag" -> normalizedTag,
    "RequestedMin" -> reqMin, "RequestedMax" -> reqMax,
    "Request" -> request, "PersistentMetadata" -> persistentMetadata,
    "LocalMetadata" -> localMetadata,
    "CheckpointIdentity" -> checkpointIdentity|>;
  If[!TrueQ[retainValueSolverPrototype], Return[prepared, Module]];
  valueRunKeys = {"nmax", "p", "has_initial", "adaptive_probe",
    "a_target", "b_target", "a_shift_min", "a_shifts", "schedule",
    "initial", "initial_validity", "source", "return_u"};
  Append[prepared, "ValueSolver" -> <|
      "schema" -> "diffexp2-native-regular-value-solver-prototype-v1",
      "run" -> KeyTake[request, valueRunKeys],
      "metadata" -> localMetadata,
      "tail_proxy_max_exact" -> ToString[
        DiffExp2`Tolerances`Tol["LaurentLeadTol"]/100, InputForm],
      "relative_accuracy_max_exact" ->
        cppRegularValueRelativeAccuracyMaxExact[]|>]];

prepareNativeLocalFamilyRun[cs_Association, req_Association,
    tag_Association, init_List, retainValueSolverPrototype_:False] := Module[
  {shared = prepareNativeLocalFamilyShared[cs, req, tag, init]},
  prepareNativeLocalFamilyDynamic[
    shared, init, retainValueSolverPrototype]];

nativeLocalFamilyFinalize[cs_Association, req_Association,
    prepared_Association, response_, expectedSession_:Automatic,
    expectedChart_:Automatic, expectedOperator_:Automatic] := Module[
  {backendID, result, reqMax = prepared["RequestedMax"],
   checkpointIdentity = prepared["CheckpointIdentity"],
   normalizedTag = prepared["Tag"], forbiddenPayloadKeys},
  If[FailureQ[response],
    err["E5", cs, <|"BackendFailure" -> response, "Detail" ->
      "persistent native local solve failed"|>]];
  If[!AssociationQ[response] ||
      Lookup[response, "status", "error"] =!= "ok",
    backendID = Lookup[response, "id", "E5"];
    err[If[MemberQ[{"E4", "E5", "E6"}, backendID], backendID, "E5"],
      cs, <|"BackendID" -> backendID,
        "Detail" -> Lookup[response, "detail",
          "persistent native local solve returned an error"]|>]];
  forbiddenPayloadKeys = Intersection[Keys[response],
    {"assembled", "coefficients", "u", "validity"}];
  If[!StringQ[Lookup[response, "session", None]] ||
      !StringQ[Lookup[response, "local", None]] ||
      !StringQ[Lookup[response, "chart", None]] ||
      (expectedSession =!= Automatic &&
        response["session"] =!= expectedSession) ||
      (expectedChart =!= Automatic && response["chart"] =!= expectedChart) ||
      Lookup[response, "source_operator_identity", None] =!=
        If[expectedOperator === Automatic,
          prepared["PersistentMetadata", "PreparedToken"],
          expectedOperator] ||
      Lookup[response, "dimension", None] =!= prepared["Dimension"] ||
      !IntegerQ[Lookup[response, "epsilon_min", None]] ||
      !IntegerQ[Lookup[response, "epsilon_max", None]] ||
      response["epsilon_min"] > response["epsilon_max"] ||
      response["epsilon_max"] < reqMax ||
      Lookup[response, "taylor_complete_max", None] =!= req["TOrder"] ||
      Lookup[response, "checkpoint_identity", None] =!=
        checkpointIdentity ||
      !AssociationQ[Lookup[response, "metadata", None]] ||
      !TrueQ[Lookup[response, "native_retained", False]] ||
      Lookup[response, "json_coefficients", None] =!= 0 ||
      Lookup[response, "pseudo_hit_count", None] =!= 0 ||
      forbiddenPayloadKeys =!= {},
    If[StringQ[Lookup[response, "local", None]],
      Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[response]]];
    err["E6", cs, <|"BackendResponse" -> response,
      "ExpectedCheckpointIdentity" -> checkpointIdentity,
      "ForbiddenPayloadKeys" -> forbiddenPayloadKeys,
      "Detail" ->
        "native local family summary violated its exact identity or opaque retained-local contract"|>]];
  result = <|"Type" -> "DiffExp2NativeLocalFamily",
    "Session" -> response["session"], "Local" -> response["local"],
    "NativeChart" -> response["chart"],
    "Tag" -> normalizedTag,
    "Chart" -> <|"Center" -> cs["Center"],
      "ChartMap" -> cs["ChartMap"], "Radius" -> cs["Radius"],
      "Prescriptions" -> cs["Prescriptions"]|>,
    "EpsWindow" -> <|"Min" -> response["epsilon_min"],
      "CompleteMax" -> response["epsilon_max"]|>,
    "TWindow" -> <|"CompleteMax" ->
      response["taylor_complete_max"]|>,
    "CheckpointIdentity" -> checkpointIdentity,
    "NativeSummary" -> KeyDrop[response,
      {"status", "session", "local", "chart", "metadata"}]|>;
  result];

(* First production seam for a session-owned native LocalSolution.  It is
   intentionally narrower than SolveHomogeneous: no SCC orchestration,
   pseudo-resonant compensation, or rank-reduction gauge is hidden behind
   the opaque handle.  Those operations still require Wolfram tensors until
   their native counterparts preserve the same sequential completeness
   contract. *)
SolveNativeLocalFamily[cs_Association, req_Association,
    tag_Association, init_List, retainValueSolverPrototype_:False,
    equationOwner_:Automatic] := Module[
  {prepared, response, result, expectedSession, expectedChart,
   expectedOperator},
  prepared = prepareNativeLocalFamilyRun[cs, req, tag, init,
    retainValueSolverPrototype];
  response = DiffExp2`CppBackend`RunPersistentLocalSolve[
    prepared["Request"], prepared["PersistentMetadata"],
    prepared["LocalMetadata"], equationOwner];
  expectedSession = If[equationOwner === Automatic, Automatic,
    Lookup[equationOwner, "Session", Lookup[equationOwner, "session", None]]];
  expectedChart = If[equationOwner === Automatic, Automatic,
    Lookup[equationOwner, "NativeEquationOwner",
      Lookup[equationOwner, "EquationOwner", None]]];
  expectedOperator = If[equationOwner === Automatic, Automatic,
    Lookup[equationOwner, "EquationIdentity",
      Lookup[equationOwner, "ChartIdentity", None]]];
  result = nativeLocalFamilyFinalize[cs, req, prepared, response,
    expectedSession, expectedChart, expectedOperator];
  If[!TrueQ[retainValueSolverPrototype], Return[result, Module]];
  Join[result, <|"NativeValueSolver" -> prepared["ValueSolver"]|>]];

solveHomogeneousCore[cs_Association, req_Association] := Module[
  {d = cs["SystemSize"], blocks = blockList[cs], nmax, reqMax,
   columns = {}, workColumns = {}, specs = {}, hitsAll = {}, compAll = {},
   fams = cs["Families"], colCursor = 0, wideTop,
   symbolic, singleUseDepth,
   startFb, terminalFb, adaptiveQ, prepCache = <||>, prepFor,
   vPrepCache = <||>, vPrepFor,
   adaptiveDiags = {}, certs = {}, cppTimingQ, solveStarted, residualStarted,
   knownAdaptiveFb, cppBatchQ, cppTasks = {}, cppBatchResults = {},
   cppBatchCursor = 0, taskCursor, framePlan, forcedFrame},
  cppTimingQ = Environment["DEBUG_CPP_RECURRENCE"] === "1";
  solveStarted = SessionTime[];
  nmax = req["TOrder"];
  reqMax = req["EpsWindow", "CompleteMax"];
  framePlan = homogeneousFramePlanFor[cs, req];
  symbolic = framePlan["Symbolic"];
  singleUseDepth = framePlan["SingleUseDepth"];
  wideTop = framePlan["FrameTop"];
  terminalFb = framePlan["TerminalFrameBase"];
  startFb = framePlan["StartFrameBase"];
  adaptiveQ = TrueQ[framePlan["Adaptive"]];
  forcedFrame = $cppHomogeneousFrameOverride;
  If[AssociationQ[forcedFrame],
    If[!IntegerQ[Lookup[forcedFrame, "FrameBase", None]] ||
        !IntegerQ[Lookup[forcedFrame, "FrameTop", None]] ||
        (forcedFrame["FrameBase"] > terminalFb &&
          forcedFrame["FrameBase"] > startFb) ||
        forcedFrame["FrameTop"] < wideTop,
      err["E6", cs, <|"RequiredFrame" -> {terminalFb, wideTop},
        "ForcedFrame" -> forcedFrame,
        "Detail" -> "forced homogeneous capture frame does not contain the exact chart work rectangle"|>]];
    terminalFb = forcedFrame["FrameBase"];
    startFb = terminalFb;
    wideTop = forcedFrame["FrameTop"];
    adaptiveQ = False];
  (* Once any column proves that a wider lower rectangle is required, every
     later column starts from that already-certified width.  A wider frame is
     algebraically identical and avoids repeating the same failed narrow
     probes (especially costly across a LibraryLink boundary). *)
  knownAdaptiveFb = startFb;
  prepFor[ff_Integer] := If[KeyExistsQ[prepCache, ff], prepCache[ff],
    AssociateTo[prepCache, ff -> If[
      TrueQ[Lookup[framePlan, "EpsilonRegularPrincipal", False]],
      prepareEpsilonRegularCleared[cs, ff, wideTop - ff + 1,
        Join[framePlan["EpsilonRegularCertificate"],
          <|"Eligible" -> True, "Symbolic" -> symbolic|>]],
      prepareCleared[cs, ff, wideTop - ff + 1, symbolic]]];
    prepCache[ff]];
  vPrepFor[ff_Integer, ww_Integer] := Module[{key = {ff, ww}},
    If[KeyExistsQ[vPrepCache, key], vPrepCache[key],
      AssociateTo[vPrepCache, key -> prepareFramedMatrix[
        cs["V"], DiffExp2`Config`CanonicalEps[], ff, ww, cs]];
      vPrepCache[key]]];
  cppBatchQ = cfg["RecurrenceBackend"] === "Cpp" &&
    !TrueQ[$disableGroupedSpectralTransform];
  If[cppBatchQ,
    taskCursor = 0;
    Do[Module[{fam = fams[[fi]], root, q, taskP, taskInit},
      Do[
        root = fam["Roots"][[ri]];
        q = root["BlockSize"];
        Do[
          taskP = logCeiling[cs, root["a"], root["b"], qpos];
          taskInit = Table[
            Module[{vv = Table[esZero[wideTop], {d}]},
              If[l <= qpos,
                vv[[taskCursor + (qpos + 1 - l)]] =
                  esShift[esNew[0, PadRight[{1}, wideTop + 1 + l]], -l]];
              vv], {l, 0, taskP}];
          AppendTo[cppTasks, <|"a" -> root["a"], "b" -> root["b"],
            "P" -> taskP, "Init" -> taskInit|>],
          {qpos, 0, q - 1}];
        taskCursor += q,
        {ri, Length[fam["Roots"]]}]],
      {fi, Length[fams]}];
    cppBatchResults = cppBatchRecurrences[cs, prepFor[terminalFb],
      cppTasks, nmax, terminalFb, wideTop - terminalFb + 1,
      vPrepFor[terminalFb, wideTop - terminalFb + 1]];
    knownAdaptiveFb = terminalFb];
  Do[Module[{fIdx = fi, fam = fams[[fi]]},
    Do[Module[{blk, root, q},
      (* identify the block for this root via the running cursor *)
      root = fam["Roots"][[ri]]; q = root["BlockSize"];
      Do[Module[{P, init, rec, ls, comp, fbRun = startFb, WRun,
          attempts = 0, underflow, used, nextUsed, tColumn, tRec, tAssembly,
          tCompensation, tCertificate, certificate},
        tColumn = SessionTime[];
        fbRun = knownAdaptiveFb;
        P = logCeiling[cs, root["a"], root["b"], qpos];
        (* init ladder: U[0, l] = e_{qpos+1-l} eps^{-l}, l = 0..qpos *)
        init = Table[
          Module[{vv = Table[esZero[wideTop], {d}]},
            If[l <= qpos,
              vv[[colCursor + (qpos + 1 - l)]] =
                esShift[esNew[0, PadRight[{1}, wideTop + 1 + l]], -l]];
            vv],
          {l, 0, P}];
        If[cppBatchQ,
          cppBatchCursor++;
          rec = cppBatchResults[[cppBatchCursor]];
          fbRun = terminalFb;
          WRun = wideTop - terminalFb + 1;
          attempts = 1,
          While[True,
            WRun = wideTop - fbRun + 1;
            attempts++;
            underflow = Catch[
              Block[{$adaptiveLowerFrameProbe = adaptiveQ && fbRun > terminalFb},
                runRecursion[cs, prepFor[fbRun], root["a"], root["b"],
                  P, nmax, None, fbRun, WRun, init, Automatic]],
              "DiffExp2AdaptiveLowerFrame"];
            If[!FailureQ[underflow], rec = underflow; Break[]];
            (* Monotone geometric widening.  Since terminalFb is exactly the
               previous scalar bound, termination and the old error behavior
               are preserved even for scalar/idempotent repeated poles. *)
            used = startFb - fbRun;
            nextUsed = If[used === 0, Max[1, singleUseDepth], 2 used];
            fbRun = Max[terminalFb, startFb - nextUsed]]];
        knownAdaptiveFb = Min[knownAdaptiveFb, fbRun];
        tRec = SessionTime[];
        AppendTo[adaptiveDiags, <|
          "Tag" -> {root["a"], root["b"], qpos},
          "Adaptive" -> adaptiveQ, "Attempts" -> attempts,
          "FrameBase" -> fbRun, "FrameTop" -> wideTop,
          "FrameWidth" -> WRun,
          "TerminalFrameBase" -> terminalFb,
          "TerminalFrameWidth" -> wideTop - terminalFb + 1|>];
        ls = assembleSolution[cs, root["a"], root["b"], rec, nmax,
          If[TrueQ[$disableGroupedSpectralTransform], Automatic,
            vPrepFor[fbRun, WRun]]];
        tAssembly = SessionTime[];
        comp = compensatePseudoColumn[cs, ls, rec["Hits"], workColumns, reqMax];
        ls = comp[[1]];
        tCompensation = SessionTime[];
        certificate = certifyPseudoCompensation[cs, ls, rec["Hits"],
          <|"Kind" -> "Homogeneous", "Tag" ->
            {root["a"], root["b"], qpos}|>];
        AppendTo[certs, certificate];
        tCertificate = SessionTime[];
        If[cppTimingQ, Print["CPPCOL ", <|
          "Tag" -> {root["a"], root["b"], qpos},
          "RecurrenceSeconds" -> N[tRec - tColumn, 6],
          "AssemblySeconds" -> N[tAssembly - tRec, 6],
          "CompensationSeconds" -> N[tCompensation - tAssembly, 6],
          "CertificateSeconds" -> N[tCertificate - tCompensation, 6]|>]];
        compAll = Join[compAll, Map[Append[#, "SourceColumn" ->
          (Length[workColumns] + 1)] &, comp[[2]]]];
        hitsAll = Join[hitsAll, publicPseudoHit /@ rec["Hits"]];
        (* Work-frame headroom is scratch space for Cauchy products and
           collision divisions, not a delivered accuracy claim. *)
        If[ls["EpsWindow", "CompleteMax"] < reqMax,
          err["E6", cs, <|"RequestedCompleteMax" -> reqMax,
            "DeliveredWindow" -> ls["EpsWindow"],
            "Tag" -> {root["a"], root["b"], qpos},
            "Detail" -> "homogeneous work budget did not reach the requested epsilon order"|>]];
        AppendTo[workColumns, ls];
        ls = capWindow[cs, ls, reqMax];
        AppendTo[columns, ls];
        AppendTo[specs, <|"a" -> root["a"], "b" -> root["b"], "p" -> qpos,
          "Family" -> fIdx, "Root" -> ri, "ChainPos" -> qpos|>]],
        {qpos, 0, q - 1}];
      colCursor += q],
      {ri, Length[fam["Roots"]]}]],
    {fi, Length[fams]}];
  Module[{fs = <|"Columns" -> columns, "Specs" -> specs,
      "Diagnostics" -> <|"PseudoCollisionsHit" -> hitsAll,
        "PseudoCompensations" -> compAll,
        "PseudoCollisionsCompensated" -> And @@ certs,
        "AdaptiveLowerFrames" -> adaptiveDiags,
        "CouplingDepth" -> 0|>|>},
    residualStarted = SessionTime[];
    ODEResidualCheck[cs, fs];
    If[cppTimingQ, Print["CPPSOLVE ", <|
      "BeforeResidualSeconds" -> N[residualStarted - solveStarted, 6],
      "ResidualSeconds" -> N[SessionTime[] - residualStarted, 6]|>]];
    fs]];

(* A caller supplies a theta-form source in the ORIGINAL physical frame,
   theta f = B_original.f + s_f.  On a rank-reduced chart f = T.g, hence

       theta g = B_reduced.g + T^-1.s_f.

   This transform is deliberately performed once, at the SolveParticular
   boundary.  In particular there is no additional factor of t (the source
   is already in theta form), and SCC/block callers must not pre-transform
   their physical coupling sources.

   Use SectorSeries' rational multiplier so center t-poles become exact
   shifts of the sector's a tag (analytic factors remain Taylor series),
   while b/log tags and honest epsilon windows follow the existing closed
   algebra.  Structurally zero source components
   are skipped: a pole in an unused T^-1 entry must not consume completeness.
   The result is still SourceData (only the three source keys are exposed). *)
sourceLocalSolution[cs_, source_] := Module[{ls, d = cs["SystemSize"], ncomp},
  ls = <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"], "Sectors" -> source["Sectors"],
    "EpsWindow" -> source["EpsWindow"], "TWindow" -> source["TWindow"],
    "ErrorEstimate" -> ConstantArray[0,
      source["EpsWindow", "CompleteMax"] - source["EpsWindow", "Min"] + 1],
    "Prescriptions" -> cs["Prescriptions"]|>;
  ls = DiffExp2`SectorSeries`ValidateLocalSolution[ls];
  ncomp = Dimensions[First[ls["Sectors"]]["Coeffs"]][[3]];
  If[ncomp =!= d,
    err["E9", cs, <|"SourceComponents" -> ncomp, "SystemSize" -> d,
      "Detail" -> "particular source component dimension does not match the chart system"|>]];
  ls];

reduceParticularSource[cs_, source_] := Module[
  {d = cs["SystemSize"], t = cs["ChartVar"], TInv = cs["GaugeInverse"],
   ls, active, terms, out},
  (* SameQ identity is load-bearing parity: ordinary charts retain their
     SourceData byte-for-byte and pay no wrapper/canonicalization churn. *)
  If[TInv === IdentityMatrix[d], Return[source]];
  ls = sourceLocalSolution[cs, source];
  active = Table[AnyTrue[ls["Sectors"], Function[sec,
      AnyTrue[Flatten[sec["Coeffs"][[All, All, c]]],
        !(FreeQ[#, _?InexactNumberQ] && zeroQ[#]) &]]], {c, d}];
  terms = Flatten[Table[
    If[zeroQ[TInv[[r, c]]] || !TrueQ[active[[c]]], {},
      {DiffExp2`SectorSeries`MultiplyRational[
        Append[ls, "Sectors" -> Map[Function[sec,
          Append[sec, "Coeffs" -> Map[
            Function[v, ReplacePart[ConstantArray[0, d], r -> v[[c]]]],
            sec["Coeffs"], {2}]]], ls["Sectors"]]],
        TInv[[r, c]], t]}],
    {r, d}, {c, d}], 2];
  (* A nonempty SourceData can still have structurally zero coefficient
     slabs.  Its value is frame-independent, so preserve the input window. *)
  If[terms === {}, Return[source]];
  out = If[Length[terms] === 1, First[terms],
    DiffExp2`SectorSeries`CombineLocalSolutions[
      ConstantArray[1, Length[terms]], terms]];
  KeyTake[out, {"Sectors", "EpsWindow", "TWindow"}]];

SolveParticular[cs_Association, source_Association, req_Association] := Module[
  {d = cs["SystemSize"], nmax, reqMin, reqMax, parts = {}, wideTop, prep,
   eps = DiffExp2`Config`CanonicalEps[], hitsAll = {}, compAll = {}, certs = {},
   homTargets = None, symbolic, poleDepth, spectralDepth,
   inverseSpectralDepth, sourceFrameDepth, singleUseDepth, reducedSource,
   vPrepCache = <||>, vPrepFor, blockSystems, execution, fullcs, result},
  If[TrueQ[Lookup[cs, "SCCSkeleton", False]] &&
      AssociationQ[Lookup[cs, "IntegrationSequence", None]] &&
      Length[cs["IntegrationSequence", "Components"]] > 1,
    blockSystems = sccBlockChartSystem[cs, #] & /@
      Range[Length[cs["IntegrationSequence", "Components"]]];
    execution = sccExecutionPlan[cs, req, blockSystems];
    If[execution["Mode"] === "MonolithicStrict",
      fullcs = sccTryMonolithicChartSystem[cs];
      If[FailureQ[fullcs],
        result = solveSCCParticular[cs, source, req, blockSystems];
        Return[sccAnnotateDeclinedParticular[
          result, cs, execution, fullcs], Module],
        result = SolveParticular[fullcs, source, req];
        Return[sccAnnotateCoarsenedParticular[
          result, cs, execution], Module]],
      Return[solveSCCParticular[cs, source, req, blockSystems], Module]]];
  nmax = Min[req["TOrder"], source["TWindow", "CompleteMax"]];
  reqMin = req["EpsWindow", "Min"]; reqMax = req["EpsWindow", "CompleteMax"];
  If[source["Sectors"] === {},
    Return[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"],
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> Table[0, {reqMax - reqMin + 1}, {nmax + 1}, {d}]|>},
      "EpsWindow" -> <|"Min" -> reqMin, "CompleteMax" -> reqMax|>,
      "TWindow" -> <|"CompleteMax" -> nmax|>,
      "ErrorEstimate" -> ConstantArray[0, reqMax - reqMin + 1],
      "Prescriptions" -> cs["Prescriptions"]|>]];
  reducedSource = reduceParticularSource[cs, source];
  (* Rational source multiplication currently preserves the finite Taylor
     width, but derive the recursion bound from the transformed contract:
     this stays correct if that algebra later returns a stricter window. *)
  nmax = Min[req["TOrder"], reducedSource["TWindow", "CompleteMax"]];
  symbolic = clearedSymbolic[cs];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  singleUseDepth = recurrenceSingleUsePoleDepth[symbolic];
  spectralDepth = spectralTransformPoleDepth[cs];
  inverseSpectralDepth = inverseSpectralTransformPoleDepth[cs];
  sourceFrameDepth = spectralDepth + inverseSpectralDepth;
  vPrepFor[ff_Integer, ww_Integer] := Module[{key = {ff, ww}},
    If[KeyExistsQ[vPrepCache, key], vPrepCache[key],
      AssociateTo[vPrepCache, key -> prepareFramedMatrix[
        cs["V"], eps, ff, ww, cs]];
      vPrepCache[key]]];
  Do[Module[{aS = sec["a"], bS = sec["b"], pS = sec["p"], arr = sec["Coeffs"],
      P, srcMin, srcMax, VInvExp, VInvVal, bHat, bHatValid, srcFn, rec, ls,
      wideTop2, prep2, pseudoDepth, comp, desiredMax, assemblyPrep},
    P = logCeiling[cs, aS, bS, pS, True];  (* sources: Z>=0 incl. the same-a hit *)
    pseudoDepth = pseudoDepthForTag[cs, aS, bS, 0];
    srcMin = reducedSource["EpsWindow", "Min"];
    srcMax = reducedSource["EpsWindow", "CompleteMax"];
    wideTop2 = srcMax + P + pseudoDepth + 2 - Min[0, srcMin - P - 2];
    wideTop2 = Max[wideTop2,
      srcMax + P + pseudoDepth + poleDepth + 2 singleUseDepth];
    wideTop2 += sourceFrameDepth;
    Module[{fb2 = Min[Min[srcMin, 0] - P - pseudoDepth - 2 - 2,
          srcMin - P - pseudoDepth - poleDepth - 2 singleUseDepth], Wd2},
      (* The source first crosses VInv, then the solved coefficient crosses
         V on assembly.  Both matrices live in the shared frame, so their
         lower pole depths and upper scratch losses add. *)
      fb2 -= sourceFrameDepth;
      Wd2 = wideTop2 - fb2 + 1;
      prep2 = prepareCleared[cs, fb2, Wd2, symbolic];
      VInvExp = Map[ratEpsList[Together[#], eps, fb2, Wd2] &, cs["VInv"], {2}];
      VInvVal = Map[frameValuation[#, fb2] &, VInvExp, {2}];
      (* J-frame source as FRAME LISTS: bHat[n+1][r] *)
      bHat = Table[
        Module[{vc = Table[
            Module[{fl = ConstantArray[0, Wd2]},
              Do[
                If[1 <= k - fb2 + 1 <= Wd2,
                  fl[[k - fb2 + 1]] =
                    arr[[k - srcMin + 1, n + 1, c]]],
                {k, srcMin, srcMax}];
              fl],
            {c, d}]},
          Table[
            certifiedFrameChop /@
              Sum[frConv[VInvExp[[r, c]], vc[[c]], fb2, Wd2],
                {c, d}],
            {r, d}]],
        {n, 0, Min[nmax, Length[arr[[1]]] - 1]}];
      (* The source is finite data: zero frame slots above srcMax are
         unknown, not an extension of completeness.  Transform its actual
         top through the exact VInv entry valuations. *)
      bHatValid = Table[Table[validMin[Table[
          If[VInvVal[[r, c]] === Infinity, Infinity,
            validShift[srcMax, VInvVal[[r, c]], fb2 + Wd2 - 1]],
          {c, d}]], {r, d}], {Length[bHat]}];
      (* hand frame lists straight through: wrap as a function returning
         pre-framed vectors via a marker *)
      (* fresh formal names: bare Private n/l carry values from package
         evaluation and corrupt Function parameter lists *)
      srcFn = With[{bb = bHat, bv = bHatValid, pp = pS},
        Function[{srcN, srcL},
          If[srcL === pp && srcN + 1 <= Length[bb],
            <|"Frames" -> bb[[srcN + 1]],
              "Validity" -> bv[[srcN + 1]]|>, None]]];
      (* The native recurrence already supports applying V before returning
         its result.  Particular solves previously omitted this payload,
         decoded the much wider J-frame U tensor, and then performed the
         same V convolution in Wolfram.  Supplying the existing prepared V
         keeps the exact frame/window decisions unchanged while returning
         only the assembled physical coefficients. *)
      assemblyPrep = If[TrueQ[$disableGroupedSpectralTransform], Automatic,
        vPrepFor[fb2, Wd2]];
      rec = runRecursionFramedSrc[
        cs, prep2, aS, bS, P, nmax, srcFn, fb2, Wd2, None,
        assemblyPrep]];
    hitsAll = Join[hitsAll, publicPseudoHit /@ rec["Hits"]];
    ls = assembleSolution[cs, aS, bS, rec, nmax, assemblyPrep];
    desiredMax = Min[ls["EpsWindow", "CompleteMax"], reqMax];
    If[rec["Hits"] =!= {},
      If[homTargets === None,
        Module[{hreq = Join[req, <|"EpsWindow" ->
            Join[req["EpsWindow"], <|"CompleteMax" ->
              Max[reqMax, reducedSource["EpsWindow", "CompleteMax"]]|>]|>]},
          homTargets = SolveHomogeneous[cs, hreq]["Columns"]]];
      comp = compensatePseudoColumn[cs, ls, rec["Hits"], homTargets, desiredMax];
      ls = comp[[1]];
      compAll = Join[compAll, Map[Append[#, "SourceSector" ->
        {aS, bS, pS}] &, comp[[2]]]];
      AppendTo[certs, certifyPseudoCompensation[cs, ls, rec["Hits"],
        <|"Kind" -> "Particular", "Tag" -> {aS, bS, pS}|>]],
      AppendTo[certs, True]];
    AppendTo[parts, ls]],
    {sec, reducedSource["Sectors"]}];
  Module[{ls = If[Length[parts] === 1, First[parts],
      DiffExp2`SectorSeries`CombineLocalSolutions[
        ConstantArray[1, Length[parts]], parts]]},
    ls = capWindow[cs, ls, reqMax];
    ls = Join[ls, <|"Diagnostics" -> <|
      "PseudoCollisionsHit" -> hitsAll,
      "PseudoCompensations" -> compAll,
      "PseudoCollisionsCompensated" -> And @@ certs|>|>];
    ODEResidualCheck[cs, ls, source];
    ls]];

(* ---- exact SCC/block integration sequence ---------------------------- *)
sccBlockChartSystem[cs_Association, block_Integer] := Module[
  {seq = cs["IntegrationSequence"], indices, t = cs["ChartVar"],
   parentInput, sys, chart, sub, clearKey = cs["SystemClearKey"]},
  indices = seq["Components"][[block]];
  If[!KeyExistsQ[$systemClearRegistry, clearKey],
    err["E6", cs, <|"SystemClearKey" -> clearKey,
      "Detail" -> "SCC envelope references an unregistered original system"|>]];
  parentInput = $systemClearRegistry[clearKey];
  (* Prepare each diagonal block from the original physical system and the
     original affine chart.  Besides avoiding a synthetic theta/t system,
     this gives every center of the same block one shared global-clearing
     cache key. *)
  sys = <|"Variable" -> parentInput["Variable"],
    "Matrix" -> parentInput["Matrix"][[indices, indices]]|>;
  chart = <|"Center" -> cs["ChartMap", "Center"],
    "Scale" -> cs["ChartMap", "Scale"],
    "Radius" -> cs["Radius"], "LocalRadius" -> cs["Radius"],
    "ChartVar" -> t,
    "Name" -> "scc" <> ToString[block] <> "@" <>
      ToString[cs["Center"], InputForm],
    "Prescriptions" -> cs["Prescriptions"]|>;
  sub = PrepareChart[sys, chart];
  (* Use the parent's collision-certified exact system key, matching ordinary
     charts of the same physical system.  This lets a regular anchor and a
     singular composite receiving chart coexist in one retained session;
     the old machine-sized SystemHash tag unnecessarily split that atlas. *)
  Join[sub, <|"SolveCacheTag" -> Lookup[cs, "SystemClearKey",
      Lookup[cs, "SystemHash", None]],
    "SCCParentHash" -> seq["MatrixHash"], "SCCBlock" -> block,
    "SCCIndices" -> indices|>]];

(* SCC propagation is advantageous only while an off-diagonal source can be
   carried through a target block at finite epsilon width.  A negative
   epsilon recurrence multiplier at such a target makes every upstream
   column fund future epsilon orders repeatedly.  Coarsen that chart to one
   original-system recurrence instead: this is the same strict solver and
   exact matrix, not a fallback to a different integration algorithm.  The
   exact SCC certificate remains attached as planning/diagnostic metadata. *)
sccMonolithicChartSystem[cs_Association] := Module[
  {clearKey = cs["SystemClearKey"], parentInput, chart, full},
  If[!KeyExistsQ[$systemClearRegistry, clearKey],
    err["E6", cs, <|"SystemClearKey" -> clearKey,
      "Detail" -> "SCC envelope references an unregistered original system"|>]];
  parentInput = $systemClearRegistry[clearKey];
  chart = <|"Center" -> cs["ChartMap", "Center"],
    "Scale" -> cs["ChartMap", "Scale"],
    "Radius" -> cs["Radius"], "LocalRadius" -> cs["Radius"],
    "ChartVar" -> cs["ChartVar"],
    "Name" -> "scc-monolithic@" <> ToString[cs["Center"], InputForm],
    "Prescriptions" -> cs["Prescriptions"],
    "UseSCCSkeleton" -> False|>;
  full = PrepareChart[parentInput, chart];
  Join[full, <|"SolveCacheTag" -> Lookup[cs, "SystemClearKey",
    Lookup[cs, "SystemHash", None]]|>]];

(* A monolithic full-frame preparation is only a performance candidate.
   Catch precisely that preparation attempt so an unsupported global
   indicial frame can decline coarsening without swallowing any later strict
   block solve failure.  Its printed DE2 error is replaced by the retained
   Failure object in the execution diagnostics. *)
sccTryMonolithicChartSystem[cs_Association] :=
  Block[{Print = Function[Null]},
    Quiet[Catch[sccMonolithicChartSystem[cs], "DiffExp2Error"]]];

sccExecutionPlan[cs_Association, req_Association,
    blockSystems_List] := Module[
  {seq = cs["IntegrationSequence"], targets, nmax, records, offending,
   incoming},
  targets = Sort[DeleteDuplicates[
    If[seq["CondensationEdges"] === {}, {},
      seq["CondensationEdges"][[All, 2]]]]];
  nmax = sccWorkTOrder[cs, req];
  incoming[block_] := Sort[Cases[
    seq["CondensationEdges"], {source_, block} :> source]];
  records = Table[<|"TargetBlock" -> block,
      "IncomingBlocks" -> incoming[block],
      "Indices" -> seq["Components"][[block]],
      "RecurrencePoleDepth" -> recurrencePoleDepth[
        clearedSymbolic[blockSystems[[block]]], nmax]|>,
    {block, targets}];
  offending = Select[records, # ["RecurrencePoleDepth"] > 0 &];
  If[offending === {},
    <|"Mode" -> "BlockSequentialStrict", "WorkTOrder" -> nmax,
      "DownstreamTargets" -> records|>,
    <|"Mode" -> "MonolithicStrict",
      "Reason" -> "DownstreamRecurrencePoleDepth",
      "WorkTOrder" -> nmax, "OffendingTargets" -> offending|>]];

sccCoarseningDiagnostics[cs_Association, execution_Association] := <|
  "SCCSolved" -> False,
  "SCCExecutionMode" -> "MonolithicStrict",
  "SCCExecutionReason" -> execution["Reason"],
  "SCCExecutionPlan" -> execution,
  "SCCCertificate" -> cs["IntegrationSequence"],
  "CouplingDepth" -> cs["IntegrationSequence", "CouplingDepth"]|>;

sccAnnotateCoarsenedBasis[fs_Association, cs_Association,
    execution_Association, fullcs_Association] := Module[{groups},
  (* Recombination must use the actual monolithic indicial frame.  The lazy
     SCC envelope intentionally has no global Families and therefore cannot
     serve as Transport's fallback source of degeneracy groups. *)
  groups = sccBlockRecombineRecords[
    Join[fullcs, <|"SCCBlock" -> "MonolithicStrict"|>], fs, 0];
  Join[fs, <|"Diagnostics" -> Join[
    Lookup[fs, "Diagnostics", <||>],
    sccCoarseningDiagnostics[cs, execution],
    <|"SCCRecombineGroups" -> groups|>]|>]];

sccAnnotateCoarsenedParticular[ls_Association, cs_Association,
    execution_Association] := Join[ls, <|"Diagnostics" -> Join[
      Lookup[ls, "Diagnostics", <||>],
      sccCoarseningDiagnostics[cs, execution]]|>];

sccDeclinedDiagnostics[cs_Association, execution_Association,
    failure_?FailureQ] := <|
  "SCCExecutionMode" -> "BlockSequentialStrict",
  "SCCExecutionReason" -> "CoarseningDeclined",
  "SCCCoarseningCandidate" -> execution,
  "CoarseningDeclined" -> failure,
  "SCCCertificate" -> cs["IntegrationSequence"],
  "CouplingDepth" -> cs["IntegrationSequence", "CouplingDepth"]|>;

sccAnnotateDeclinedBasis[fs_Association, cs_Association,
    execution_Association, failure_?FailureQ] := Join[fs, <|
  "Diagnostics" -> Join[Lookup[fs, "Diagnostics", <||>],
    sccDeclinedDiagnostics[cs, execution, failure]]|>];

sccAnnotateDeclinedParticular[ls_Association, cs_Association,
    execution_Association, failure_?FailureQ] := Join[ls, <|
  "Diagnostics" -> Join[Lookup[ls, "Diagnostics", <||>],
    sccDeclinedDiagnostics[cs, execution, failure]]|>];

sccBlockCouplingSource[cs_Association, target_Integer,
    state_Association] := Module[
  {components = cs["IntegrationSequence", "Components"], targetIndices,
   pieces = {}, sourceBlock, mat, product},
  targetIndices = components[[target]];
  Do[
    sourceBlock = key;
    mat = cs["ThetaOriginal"][[targetIndices, components[[sourceBlock]]]];
    If[AnyTrue[Flatten[mat], !sccStructuralZeroQ[Cancel[Together[#]]] &],
      product = sccMatrixTimesLocalSolution[cs, mat, state[key]];
      If[AssociationQ[product],
        AppendTo[pieces, sccLocalSolutionSource[product]]]],
    {key, Sort[Keys[state]]}];
  sccCombineSources[cs, pieces, Length[targetIndices]]];

sccCombineBlockState[cs_Association, state_Association] := Module[
  {components = cs["IntegrationSequence", "Components"], terms},
  terms = Map[Function[key,
    sccEmbedComponents[state[key], components[[key]], cs["SystemSize"]]],
    Sort[Keys[state]]];
  If[terms === {},
    err["E6", cs, <|"Detail" -> "SCC block state is empty"|>]];
  If[Length[terms] === 1, First[terms],
    DiffExp2`SectorSeries`CombineLocalSolutions[
      ConstantArray[1, Length[terms]], terms]]];

sccBlockRecombineRecords[subcs_Association, fs_Association,
    columnOffset_Integer] := Module[{records = {}, specs = fs["Specs"], cols},
  Do[
    cols = Select[Range[Length[specs]],
      specs[[#]]["Family"] === rec["FamilyIndex"] &&
        specs[[#]]["ChainPos"] === 0 &];
    Do[If[Length[group] >= 2,
      AppendTo[records, <|"Columns" -> (columnOffset + group),
        "EpsZeroDegeneracy" -> rec["EpsZeroDegeneracy"],
        "SCCBlock" -> subcs["SCCBlock"]|>]],
      {group, GatherBy[cols, Together[specs[[#]]["a"]] &]}],
    {rec, DiffExp2`Indicial`EpsDegenerateFamilies[subcs["IndicialData"]]}];
  records];

sccParticularBlockCost[target_Association, nmax_Integer:0] := Module[
  {certificate = epsilonRegularPrincipalCertificate[target, nmax],
   symbolic},
  symbolic = If[TrueQ[certificate["Eligible"]],
    certificate["Symbolic"], clearedSymbolic[target]];
  spectralTransformPoleDepth[target] +
    inverseSpectralTransformPoleDepth[target] +
    (* A sourced resonant layer crosses the exact V^-1/d0 boundary before
       the spectral solve.  Keep one complete coefficient beyond the
       valuation shift so the round trip can certify, rather than merely
       store, the requested physical top order. *)
    If[TrueQ[certificate["Eligible"]], 1, 0] +
    matrixEpsPoleDepth[target["Gauge"]] +
    matrixEpsPoleDepth[target["GaugeInverse"]] +
    recurrencePoleDepth[symbolic, nmax] +
    2 recurrenceSingleUsePoleDepth[symbolic] +
    (* The canonical particular/log ladder can consume one honest source
       order per target Jordan position even for a regular zero diagonal. *)
    target["SystemSize"] +
    Max[0, Sequence @@ (# ["CollisionDepth"] & /@ target["Families"])] ];

sccSeedWorkHalos[cs_Association, blockSystems_List, nmax_Integer:0] := Module[
  {seq = cs["IntegrationSequence"], nb, edges, targets, targetCosts,
   edgeCost, outgoing, path, seedGaugeCosts, topological},
  nb = Length[seq["Components"]];
  edges = seq["CondensationEdges"];
  targets = DeleteDuplicates[If[edges === {}, {}, edges[[All, 2]]]];
  (* The target recurrence/gauge cost is independent of the incoming source
     edge.  Dense SCC condensations used to rebuild the same exact algebraic
     clear once for every predecessor, which turned one singular chart into
     minutes of identical PolynomialGCD/coefficient work. *)
  targetCosts = AssociationMap[
    sccParticularBlockCost[blockSystems[[#]], nmax] &, targets];
  edgeCost[{u_, v_}] := Module[{mat},
    mat = cs["ThetaOriginal"][[seq["Components"][[v]],
      seq["Components"][[u]]]];
    matrixEpsPoleDepth[mat] + targetCosts[v]];
  outgoing = Table[Cases[edges,
      {source_, target_} /; source === block :> target],
    {block, nb}];
  topological = seq["TopologicalOrder"];
  If[Sort[topological] =!= Range[nb],
    err["E6", cs, <|"TopologicalOrder" -> topological,
      "Detail" -> "SCC seed-work planning requires a complete topological order"|>]];
  path = ConstantArray[0, nb];
  Do[If[outgoing[[u]] =!= {},
    path[[u]] = Max[(edgeCost[{u, #}] + path[[#]]) & /@ outgoing[[u]]]],
    {u, Reverse[topological]}];
  seedGaugeCosts = matrixEpsPoleDepth[# ["Gauge"]] & /@ blockSystems;
  MapThread[Plus, {seedGaugeCosts, path}]];

sccInitialWorkHalo[cs_Association, blockSystems_List, nmax_Integer:0] :=
  Max[0, Sequence @@ sccSeedWorkHalos[cs, blockSystems, nmax]];

(* An external source, unlike a homogeneous seed column, must first pass
   through the particular kernel of the block where it enters.  Fund that
   first solve as well as the downstream coupling path.  The adaptive audit
   below remains authoritative: this is only the deterministic first width. *)
sccParticularInitialWorkHalo[cs_Association, blockSystems_List,
    nmax_Integer:0] :=
  sccInitialWorkHalo[cs, blockSystems, nmax] +
    Max[0, Sequence @@ (sccParticularBlockCost[#, nmax] & /@ blockSystems)];

sccWorkTOrder[cs_Association, req_Association] :=
  req["TOrder"] + 2 + 2 cs["IntegrationSequence", "CouplingDepth"];

(* ---- strict native composite SCC preparation ------------------------
   This is a preparation/capture seam only.  It deliberately does not alter
   SolveHomogeneous's production dispatch and it never selects the existing
   monolithic coarsening candidate. *)
$nativeSCCCompositeCache = <||>;
$nativeSCCCompositeCacheMax = 32;
$nativeSCCCompositeReservedCapacity = 0;

(* A singular composite starts from a compact, cancellation-audited epsilon
   rectangle.  These dynamically scoped values are retry coordinates, not
   user configuration: the public default is deliberately small, while an
   explicit frame-exhaustion witness may widen one edge without ever falling
   back to the often enormous scalar terminal bound. *)
$nativeSCCFrameLowerExtra = Automatic;
$nativeSCCFrameTopHalo = Automatic;
$nativeSCCFrameRetryLimit = 5;
$nativeSCCFrameLowerExtraMax = 16;
$nativeSCCFrameTopHaloMax = 32;

sccNativeDefaultFrameLowerExtra[cs_Association] := Min[8, Max[4,
  2 Lookup[Lookup[cs, "IntegrationSequence", <||>],
    "CouplingDepth", 0]]];

sccNativeEffectiveFrameLowerExtra[cs_Association] := If[
  IntegerQ[$nativeSCCFrameLowerExtra] &&
    0 <= $nativeSCCFrameLowerExtra <= $nativeSCCFrameLowerExtraMax,
  $nativeSCCFrameLowerExtra, sccNativeDefaultFrameLowerExtra[cs]];

sccNativeEffectiveFrameTopHalo[] := If[
  IntegerQ[$nativeSCCFrameTopHalo] &&
    1 <= $nativeSCCFrameTopHalo <= $nativeSCCFrameTopHaloMax,
  $nativeSCCFrameTopHalo, 8];

nativeSCCCompositeEffectiveCapacity[] := Max[
  $nativeSCCCompositeCacheMax, $nativeSCCCompositeReservedCapacity];

nativeSCCCompositeCacheAdmissionQ[] :=
  Length[$nativeSCCCompositeCache] < nativeSCCCompositeEffectiveCapacity[];

nativeSCCCompositeRequireAdmission[cs_Association] :=
  If[!nativeSCCCompositeCacheAdmissionQ[],
    err["E6", cs, <|"Capacity" ->
        nativeSCCCompositeEffectiveCapacity[],
      "DefaultCapacity" -> $nativeSCCCompositeCacheMax,
      "ReservedCapacity" -> $nativeSCCCompositeReservedCapacity,
      "Detail" -> "native SCC composite cache capacity is exhausted; clear solver caches before preparing another public handle"|>]];

(* Native atlas construction must retain every receiving equation owner until
   the exact tile plan has taken strong ownership.  That legitimate scoped
   demand can exceed the conservative direct-call bound, but it must not turn
   the cache into an unbounded store or evict handles which have already been
   published.  Reserve count slots from the occupancy at scope entry; nested
   schedulers retain the larger enclosing ceiling.  Block restores the prior
   policy on success, Throw, or other nonlocal exit. *)
SetAttributes[WithNativeSCCCompositeCacheReservation, HoldRest];
WithNativeSCCCompositeCacheReservation[count_Integer?NonNegative,
    expression_] := Block[{$nativeSCCCompositeReservedCapacity = Max[
      $nativeSCCCompositeReservedCapacity,
      Length[$nativeSCCCompositeCache] + count]}, expression];
WithNativeSCCCompositeCacheReservation[count_, expression_] :=
  err["E6", <||>, <|"Reservation" -> count,
    "Detail" -> "native SCC composite cache reservation must be a nonnegative integer"|>];

$nativeSCCColumnRunKeys = {"nmax", "p", "has_initial",
  "adaptive_probe", "cancellation_audit_base", "a_target", "b_target",
  "a_shift_min", "a_shifts", "schedule", "initial",
  "initial_validity", "source", "return_u"};

(* The compact execution descriptor lives inside the cache record whose
   complete parent signature is compared with SameQ.  Binding its exact
   public handle as well prevents a later explicit column solve from
   consuming run records belonging to another native composite without
   duplicating the (large) full signature or static operator tensors. *)
sccNativeCompositeExecutionDescriptorQ[entry_Association,
    signature_] := Module[{descriptor = Lookup[entry, "Execution", None]},
  AssociationQ[descriptor] &&
    SameQ[Lookup[entry, "Signature", None], signature] &&
    SameQ[Lookup[descriptor, "PublicHandle", None],
      Lookup[entry, "Result", None]] &&
    ListQ[Lookup[descriptor, "BlockDimensions", None]] &&
    ListQ[Lookup[descriptor, "Runs", None]] &&
    ListQ[Lookup[descriptor, "TaskMetadata", None]] &&
    ListQ[Lookup[descriptor, "ColumnPlans", None]] &&
    AssociationQ[Lookup[descriptor, "RationalShadowDecision", None]] &&
    AssociationQ[Lookup[descriptor, "Contract", None]] &&
    StringQ[Lookup[descriptor, "ParentIdentity", None]] &&
    ListQ[Lookup[descriptor, "Components", None]] &&
    ListQ[Lookup[descriptor, "CondensationEdges", None]] &&
    ListQ[Lookup[descriptor, "TopologicalOrder", None]] &&
    IntegerQ[Lookup[descriptor, "InputDigits", None]] &&
    AssociationQ[Lookup[descriptor, "NativeStatistics", None]]];

sccNativeCompositeStatisticsQ[stats_, handle_Association] :=
  AssociationQ[stats] && Lookup[stats, "status", "error"] === "ok" &&
    Lookup[stats, "session", None] === Lookup[handle, "Session", None] &&
    Lookup[stats, "scc", None] === Lookup[handle, "SCC", None] &&
    Lookup[stats, "key", None] === Lookup[handle, "Key", None];

sccAllSameQ[values_List] := values =!= {} &&
  AllTrue[Rest[values], SameQ[#, First[values]] &];

sccNativeCompositeCacheSignature[cs_Association, req_Association] := {
  cs, req, cfg["WorkingPrecision"], cfg["ChopPrecision"],
  cfg["Variables"], cfg["RecurrenceBackend"],
  DiffExp2`Tolerances`$InputPrecisionFactor,
  cppNativePrecisionBits[cfg["WorkingPrecision"]],
  TrueQ[$cppExactDomain], TrueQ[$cppUsePersistentSessions],
  TrueQ[$numericizeAllPreparedNumbers],
  TrueQ[$disablePreparedDirectNumericization],
  TrueQ[$disableGlobalClearedHoist],
  TrueQ[$disableIdentityNhatShortcut],
  TrueQ[$disableAdaptiveLowerFrames],
  TrueQ[$adaptiveLowerFrameProbe],
  TrueQ[$disableRationalDenominatorFusion],
  TrueQ[$disableGroupedSpectralTransform],
  TrueQ[$disablePolynomialNhatTransform],
  sccNativeEffectiveFrameLowerExtra[cs],
  sccNativeEffectiveFrameTopHalo[],
  Environment["DE2_SCC_DIAGNOSTIC_START_FRAME"],
  Environment["DE2_SCC_DIAGNOSTIC_STRICT_PROBE"],
  Environment["DE2_SCC_DIAGNOSTIC_LOWER_EXTRA"],
  Environment["DE2_SCC_DIAGNOSTIC_TOP_HALO"],
  Environment["DE2_CPP_PERSISTENT"]};

sccCaptureHomogeneousGroup[blockcs_Association,
    workReq_Association, framePlan_:Automatic,
    forcedFrame_:None] := Module[{captured},
  captured = Catch[
    Block[{$homogeneousFramePlanOverride =
          If[AssociationQ[framePlan], framePlan, None],
        $cppHomogeneousFrameOverride = forcedFrame,
        $cppHomogeneousBatchCapture = True,
        $cppHomogeneousBatchInjection = None,
        $cppHomogeneousBatchInjectionUses = 0,
        $cppBuildRequestOnly = True},
      solveHomogeneousCore[blockcs, workReq]],
    $cppHomogeneousBatchTag];
  If[!AssociationQ[captured] ||
      !ListQ[Lookup[captured, "Requests", None]] ||
      Lookup[captured, "Requests", {}] === {} ||
      !AllTrue[captured["Requests"], AssociationQ] ||
      !ListQ[Lookup[captured, "TaskMetadata", None]] ||
      Length[captured["TaskMetadata"]] =!=
        Length[captured["Requests"]] ||
      !AllTrue[captured["TaskMetadata"], AssociationQ[#] &&
        Sort[Keys[#]] === Sort[{"a", "b", "P"}] &] ||
      !AssociationQ[Lookup[captured, "Metadata", None]],
    err["E5", blockcs, <|
      "Detail" -> "diagonal SCC homogeneous solve did not yield one nonempty captured native request group"|>]];
  captured];

sccCapturedCompositeContract[captures_List, cs_Association,
    expectedNMax_Integer] := Module[
  {requests, frameRecords, auditRecords, fieldRecords, sessionRecords,
   domain, names, symbols, serialization},
  If[captures === {} || !AllTrue[captures, AssociationQ],
    err["E5", cs, <|"Detail" ->
      "native SCC capture produced no diagonal request groups"|>]];
  requests = Flatten[Lookup[captures, "Requests", {}], 1];
  If[requests === {} || !AllTrue[requests, AssociationQ],
    err["E5", cs, <|"Detail" ->
      "native SCC capture produced no recurrence requests"|>]];
  frameRecords = Map[{
      Lookup[#, "fb", None], Lookup[#, "w", None],
      Lookup[#, "nmax", None]} &, requests];
  If[!sccAllSameQ[frameRecords] ||
      !MatchQ[First[frameRecords], {_Integer, _Integer, _Integer}] ||
      First[frameRecords][[2]] < 1 ||
      First[frameRecords][[3]] =!= expectedNMax,
    err["E6", cs, <|"CapturedFrames" -> DeleteDuplicates[frameRecords],
      "ExpectedNMax" -> expectedNMax,
      "Detail" -> "native SCC first slice requires identical fb/w/nmax across every captured block request"|>]];
  auditRecords = Lookup[requests, "cancellation_audit_base", None];
  If[!sccAllSameQ[auditRecords] ||
      !MatchQ[First[auditRecords], _Integer | Null],
    err["E6", cs, <|"CapturedAuditBases" ->
        DeleteDuplicates[auditRecords],
      "Detail" -> "native SCC capture requires one identical optional cancellation-audit base across every block request"|>]];
  fieldRecords = Map[{
      Lookup[#, "domain", None], Lookup[#, "symbols", None],
      Lookup[#, "precision_bits", None],
      Lookup[#, "output_digits", None]} &, requests];
  If[!sccAllSameQ[fieldRecords] ||
      !MemberQ[{"acb", "rational", "symbolic"},
        First[fieldRecords][[1]]] ||
      !ListQ[First[fieldRecords][[2]]] ||
      !AllTrue[First[fieldRecords][[2]], StringQ] ||
      !IntegerQ[First[fieldRecords][[3]]] ||
      First[fieldRecords][[3]] < 1 ||
      !IntegerQ[First[fieldRecords][[4]]] ||
      First[fieldRecords][[4]] < 1,
    err["E5", cs, <|"CapturedFields" -> DeleteDuplicates[fieldRecords],
      "Detail" -> "native SCC first slice requires one identical captured scalar domain, symbol list, and precision"|>]];
  domain = First[fieldRecords][[1]];
  names = First[fieldRecords][[2]];
  sessionRecords = Map[{
      Lookup[#["Metadata"], "SystemIdentity", None],
      Lookup[#["Metadata"], "SessionAnalytic", None]} &, captures];
  If[!sccAllSameQ[sessionRecords] ||
      !AllTrue[captures,
        Lookup[Lookup[#["Metadata"], "SessionAnalytic", <||>],
          "RegulatorSymbols", None] === names &],
    err["E5", cs, <|"Detail" ->
      "captured diagonal blocks do not belong to one persistent session analytic identity"|>]];
  symbols = (Symbol["Global`" <> #] & /@ names);
  serialization = sccSerializationField[
    <|"domain" -> domain, "symbols" -> symbols|>, cs];
  <|"FrameBase" -> First[frameRecords][[1]],
    "FrameWidth" -> First[frameRecords][[2]],
    "NMax" -> First[frameRecords][[3]],
    "CancellationAuditBase" -> First[auditRecords],
    "Domain" -> domain, "SymbolNames" -> names,
    "Serialization" -> KeyTake[serialization, {"domain", "symbols"}]|>];

sccNativeSourceShape[cs_Association, dimension_Integer,
    fb_Integer, completeMax_Integer, tOrder_Integer] := <|
  "Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
  "Radius" -> cs["Radius"], "Prescriptions" -> cs["Prescriptions"],
  "EpsWindow" -> <|"Min" -> fb, "CompleteMax" -> completeMax|>,
  "TWindow" -> <|"CompleteMax" -> tOrder|>,
  "Dimension" -> dimension|>;

(* Native block state is g_B = V_B u_B because the retained diagonal chart
   assembles V_B before materializing its LocalSolution.  Cross-SCC matrices
   are intentionally kept in the physical g basis.  Prepare V_T^-1 once per
   target block so C_TS g_S is summed in that physical basis and only then
   converted to the target recurrence basis.  No source-side V_S belongs in
   this payload. *)
PrepareSCCSpectralSourceTransform[blockcs_Association,
    sourceShape_Association, serialization_:Automatic] := Module[
  {certificate = sccSpectralFrameCertificate[blockcs], dimension,
   t = Lookup[blockcs, "ChartVar", None], vInv, field, domain, symbols,
   inputDigits, rawEntries, encodedEntries, activeRows, activeColumns,
   identityPayload, identity, epsWindow, tWindow},
  If[!TrueQ[Lookup[certificate, "admissible", False]],
    err["E6", blockcs, <|"SpectralFrame" -> certificate,
      "Detail" -> "target spectral source transform is not exact Laurent-unimodular"|>]];
  dimension = certificate["dimension"];
  epsWindow = Lookup[sourceShape, "EpsWindow", None];
  tWindow = Lookup[sourceShape, "TWindow", None];
  If[Lookup[sourceShape, "Dimension", None] =!= dimension ||
      !AssociationQ[epsWindow] ||
      !AllTrue[{"Min", "CompleteMax"}, KeyExistsQ[epsWindow, #] &] ||
      !IntegerQ[epsWindow["Min"]] ||
      !IntegerQ[epsWindow["CompleteMax"]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"] ||
      !AssociationQ[tWindow] ||
      !IntegerQ[Lookup[tWindow, "CompleteMax", None]] ||
      tWindow["CompleteMax"] < 0,
    err["E8", blockcs, <|"SourceShape" -> sourceShape,
      "ExpectedDimension" -> dimension,
      "Detail" -> "target spectral source transform requires the complete native work rectangle"|>]];
  field = sccSerializationField[serialization, blockcs];
  domain = field["domain"];
  symbols = field["symbols"];
  vInv = Normal[blockcs["VInv"]];
  rawEntries = Flatten[Table[Module[
      {entry = Cancel[Together[vInv[[row, column]]]], prepared,
       exactEntry},
      exactEntry = DiffExp2`SectorSeries`ExactExpressionIdentity[entry, t];
      prepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
        sourceShape, entry, t, "SerializationDomain" -> domain,
        "PrepareTaylorKernels" -> False];
      If[TrueQ[prepared["ProvenZero"]], {},
        If[prepared["ExactIdentity"] =!= exactEntry ||
            prepared["CenterPoleOrder"] =!= 0,
          err["E6", blockcs, <|"Row" -> row, "Column" -> column,
            "Prepared" -> prepared, "ExactEntry" -> exactEntry,
            "Detail" -> "t-independent VInv entry did not prepare as a center-regular exact multiplier"|>]];
        {<|"Row" -> row - 1, "Column" -> column - 1,
          "ExactEntry" -> exactEntry, "Prepared" -> prepared|>}]],
    {row, dimension}, {column, dimension}], 2];
  activeRows = DeleteDuplicates[Lookup[rawEntries, "Row", {}]];
  activeColumns = DeleteDuplicates[Lookup[rawEntries, "Column", {}]];
  If[Sort[activeRows] =!= Range[0, dimension - 1] ||
      Sort[activeColumns] =!= Range[0, dimension - 1],
    err["E6", blockcs, <|"ActiveRows" -> activeRows,
      "ActiveColumns" -> activeColumns,
      "Detail" -> "certified VInv has an empty structural row or column"|>]];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*
    cfg["WorkingPrecision"];
  encodedEntries = Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    Map[<|"row" -> #["Row"], "column" -> #["Column"],
        "exact_entry" -> #["ExactEntry"],
        "multiplier" -> cppPreparedRationalMultiplierJSON[
          #["Prepared"], inputDigits, blockcs, True]|> &, rawEntries]];
  identityPayload = <|
    "schema" -> "diffexp2-scc-spectral-source-transform-identity-v1",
    "state_basis" -> "reduced-g-after-spectral-assembly",
    "target_recurrence_basis" -> "spectral-u",
    "dimension" -> dimension,
    "identity" -> certificate["identity_v"],
    "epsilon_unimodular" -> True,
    "det_epsilon_valuation" ->
      certificate["det_epsilon_valuation"],
    "v_exact_identity" -> certificate["v_exact_identity"],
    "vinv_exact_identity" -> certificate["vinv_exact_identity"],
    "det_exact_identity" -> certificate["det_exact_identity"],
    "source_window" -> <|"epsilon_min" -> epsWindow["Min"],
      "epsilon_complete_max" -> epsWindow["CompleteMax"],
      "taylor_complete_max" -> tWindow["CompleteMax"]|>,
    "serialization" -> <|"domain" -> domain,
      "symbols" -> field["symbol_identities"]|>,
    "entries" -> Map[<|"row" -> #["row"],
        "column" -> #["column"],
        "exact_entry" -> #["exact_entry"],
        "epsilon_shift" -> #["multiplier", "epsilon_shift"],
        "center_pole_order" ->
          #["multiplier", "center_pole_order"]|> &, encodedEntries]|>;
  identity = ExportString[identityPayload, "RawJSON", "Compact" -> True];
  <|"schema" -> "diffexp2-scc-spectral-source-transform-v1",
    "rows" -> dimension, "columns" -> dimension,
    "identity" -> certificate["identity_v"],
    "epsilon_unimodular" -> True,
    "det_epsilon_valuation" -> certificate["det_epsilon_valuation"],
    "v_exact_identity" -> certificate["v_exact_identity"],
    "vinv_exact_identity" -> certificate["vinv_exact_identity"],
    "det_exact_identity" -> certificate["det_exact_identity"],
    "exact_identity" -> identity, "domain" -> domain,
    "symbols" -> (SymbolName /@ symbols),
    "entries" -> encodedEntries|>];

PrepareSCCGaugeMultiplier[sourceShape_Association, entry_, t_Symbol,
    exactIdentity_String] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], epsWindow, tWindow, frameWidth,
   tmax, valuation, centerPoleOrder, coefficients, valuations, jmin,
   jtop, epsilonFrames, kernels},
  epsWindow = sourceShape["EpsWindow"];
  tWindow = sourceShape["TWindow"];
  frameWidth = epsWindow["CompleteMax"] - epsWindow["Min"] + 1;
  tmax = tWindow["CompleteMax"];
  valuation = tVal[entry, t];
  If[valuation === Infinity,
    Return[<|"EpsilonShift" -> 0, "CenterPoleOrder" -> 0,
      "TaylorKernels" -> ConstantArray[0, {frameWidth, tmax + 1}],
      "ExactIdentity" -> exactIdentity, "ProvenZero" -> True|>, Module]];
  centerPoleOrder = Max[0, -valuation];
  (* The native gauge schema consumes only this finite kernel rectangle.
     Expand in t first so rational functions of t never survive into a long
     epsilon quotient recurrence. *)
  coefficients = Table[tLaurent[entry, t, n - centerPoleOrder],
    {n, 0, tmax}];
  valuations = Select[epsValuation[#, eps] & /@ coefficients, IntegerQ];
  If[valuations === {},
    Return[<|"EpsilonShift" -> 0,
      "CenterPoleOrder" -> centerPoleOrder,
      "TaylorKernels" -> ConstantArray[0, {frameWidth, tmax + 1}],
      (* The exact entry may first appear above the retained Taylor window;
         keep its structural row/column and honest zero finite kernel. *)
      "ExactIdentity" -> exactIdentity, "ProvenZero" -> False|>, Module]];
  jmin = Min[valuations];
  jtop = jmin + frameWidth - 1;
  epsilonFrames = Map[Function[coefficient, Module[{v},
      If[zeroCanQ[coefficient], Return[ConstantArray[0, frameWidth], Module]];
      v = epsValuation[coefficient, eps];
      If[!IntegerQ[v],
        err["E4", sourceShape, <|"Entry" -> exactIdentity,
          "Detail" -> "gauge Taylor coefficient has no finite epsilon valuation"|>]];
      If[v > jtop, ConstantArray[0, frameWidth],
        ratEpsList[coefficient, eps, jmin, frameWidth]]]], coefficients];
  kernels = Transpose[epsilonFrames];
  <|"EpsilonShift" -> jmin, "CenterPoleOrder" -> centerPoleOrder,
    "TaylorKernels" -> kernels, "ExactIdentity" -> exactIdentity,
    "ProvenZero" -> False|>];

PrepareSCCGaugeTransform[blockcs_Association, role_String,
    sourceShape_Association, serialization_:Automatic,
    certificateInput_:Automatic] := Module[
  {certificate = certificateInput, dimension,
   t = Lookup[blockcs, "ChartVar", None], matrix, field, domain, symbols,
   inputDigits, rawEntries, encodedEntries, activeRows, activeColumns,
   identityPayload, identity, epsWindow, tWindow, nonzeroRules},
  If[certificate === Automatic,
    certificate = sccGaugeFrameCertificate[blockcs]];
  If[!MemberQ[{"to_physical", "to_reduced"}, role] ||
      !TrueQ[Lookup[certificate, "admissible", False]],
    err["E6", blockcs, <|"Role" -> role, "GaugeFrame" -> certificate,
      "Detail" -> "SCC gauge transform requires an exact two-sided Gauge/GaugeInverse certificate"|>]];
  dimension = certificate["dimension"];
  epsWindow = Lookup[sourceShape, "EpsWindow", None];
  tWindow = Lookup[sourceShape, "TWindow", None];
  If[Lookup[sourceShape, "Dimension", None] =!= dimension ||
      !AssociationQ[epsWindow] || !AssociationQ[tWindow] ||
      !IntegerQ[Lookup[epsWindow, "Min", None]] ||
      !IntegerQ[Lookup[epsWindow, "CompleteMax", None]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"] ||
      !IntegerQ[Lookup[tWindow, "CompleteMax", None]] ||
      tWindow["CompleteMax"] < 0,
    err["E8", blockcs, <|"SourceShape" -> sourceShape,
      "Detail" -> "SCC gauge transform requires the complete native work rectangle"|>]];
  field = sccSerializationField[serialization, blockcs];
  domain = field["domain"]; symbols = field["symbols"];
  matrix = If[role === "to_physical",
    blockcs["Gauge"], blockcs["GaugeInverse"]];
  nonzeroRules = Select[Most[ArrayRules[SparseArray[matrix]]],
    MatchQ[First[#], {_Integer, _Integer}] &&
      !TrueQ[PossibleZeroQ[Last[#]]] &];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC GAUGE PREP start block=", Lookup[blockcs, "SCCBlock", None],
      " role=", role, " entries=", Length[nonzeroRules],
      " eps=", {epsWindow["Min"], epsWindow["CompleteMax"]},
      " tmax=", tWindow["CompleteMax"], " memory=", MemoryInUse[]]];
  rawEntries = Map[Function[rule, Module[
      {position = First[rule], entry = Cancel[Together[Last[rule]]],
       row, column, prepared, exactEntry},
      row = First[position]; column = Last[position];
      exactEntry = DiffExp2`SectorSeries`ExactExpressionIdentity[entry, t];
      (* The native consumer can reconstruct this finite Taylor rectangle
         directly from the complete low-degree rational functions.  Keeping
         the custom expanded preparation here made one nontrivial endpoint
         spend minutes constructing and serializing hundreds of redundant
         exact Taylor coefficients. *)
      prepared = DiffExp2`SectorSeries`PrepareRationalMultiplier[
        sourceShape, entry, t,
        "SerializationDomain" -> domain,
        "PrepareTaylorKernels" -> False];
      If[TrueQ[prepared["ProvenZero"]], {},
        If[prepared["ExactIdentity"] =!= exactEntry,
          err["E6", blockcs, <|"Role" -> role, "Row" -> row,
            "Column" -> column,
            "Detail" -> "prepared gauge multiplier changed its exact entry identity"|>]];
        <|"Row" -> row - 1, "Column" -> column - 1,
          "ExactEntry" -> exactEntry, "Prepared" -> prepared|>]]],
    nonzeroRules];
  rawEntries = DeleteCases[rawEntries, {}];
  activeRows = DeleteDuplicates[Lookup[rawEntries, "Row", {}]];
  activeColumns = DeleteDuplicates[Lookup[rawEntries, "Column", {}]];
  If[Sort[activeRows] =!= Range[0, dimension - 1] ||
      Sort[activeColumns] =!= Range[0, dimension - 1],
    err["E6", blockcs, <|"Role" -> role,
      "Detail" -> "certified gauge transform has an empty structural row or column"|>]];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"];
  encodedEntries = Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = symbols},
    Map[<|"row" -> #["Row"], "column" -> #["Column"],
        "exact_entry" -> #["ExactEntry"],
        "multiplier" -> cppPreparedRationalMultiplierJSON[
          #["Prepared"], inputDigits, blockcs, True]|> &, rawEntries]];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC GAUGE PREP done block=", Lookup[blockcs, "SCCBlock", None],
      " role=", role, " memory=", MemoryInUse[]]];
  identityPayload = <|
    "schema" -> "diffexp2-scc-gauge-transform-identity-v1",
    "role" -> role, "dimension" -> dimension,
    "identity" -> certificate["identity_gauge"],
    "gauge_exact_identity" -> certificate["gauge_exact_identity"],
    "gauge_inverse_exact_identity" -> certificate["gauge_inverse_exact_identity"],
    "gauge_det_exact_identity" -> certificate["gauge_det_exact_identity"],
    "source_window" -> <|"epsilon_min" -> epsWindow["Min"],
      "epsilon_complete_max" -> epsWindow["CompleteMax"],
      "taylor_complete_max" -> tWindow["CompleteMax"]|>,
    "entries" -> Map[<|"row" -> #["row"], "column" -> #["column"],
        "exact_entry" -> #["exact_entry"],
        "epsilon_shift" -> #["multiplier", "epsilon_shift"],
        "center_pole_order" -> #["multiplier", "center_pole_order"]|> &,
      encodedEntries]|>;
  identity = ExportString[identityPayload, "RawJSON", "Compact" -> True];
  <|"schema" -> "diffexp2-scc-gauge-transform-v1", "role" -> role,
    "rows" -> dimension, "columns" -> dimension,
    "identity" -> certificate["identity_gauge"],
    "gauge_exact_identity" -> certificate["gauge_exact_identity"],
    "gauge_inverse_exact_identity" -> certificate["gauge_inverse_exact_identity"],
    "gauge_det_exact_identity" -> certificate["gauge_det_exact_identity"],
    "exact_identity" -> identity, "domain" -> domain,
    "symbols" -> (SymbolName /@ symbols), "entries" -> encodedEntries|>];

(* Retain the exact affine-Jordan proof used to classify every singular Acb
   recurrence.  This record contains only exact Rational task data and the
   prepared Jordan-chain partition; no serialized Acb coefficient or
   midpoint participates.  Rational C++ composites independently reconstruct
   the same certificate from their prepared operator and cross-check it,
   while Acb composites reuse this retained proof for T/P/R validation. *)
sccExactAffineJordanIndicialRecord[blockcs_Association] := Module[
  {dimension = Lookup[blockcs, "SystemSize", None], blocks,
   malformed},
  blocks = blockList[blockcs];
  malformed = Select[blocks, Function[block,
      !ListQ[Lookup[block, "Cols", None]] ||
      Lookup[block, "Cols", {}] === {} ||
      !AllTrue[Lookup[block, "Cols", {}], IntegerQ]]];
  If[!IntegerQ[dimension] || dimension < 1 || blocks === {} ||
      malformed =!= {} ||
      Sort[Flatten[Lookup[blocks, "Cols", {}]]] =!= Range[dimension],
    err["E6", blockcs, <|"UnsupportedJordanBlocks" -> malformed,
      "Detail" -> "native SCC affine-Jordan root/partition metadata is malformed"|>]];
  (* Symbolic and non-Rational algebraic composites remain preparable, but
     cannot advertise this deliberately Rational certificate/execution
     scope. *)
  If[!AllTrue[Flatten[Lookup[blocks, {"a", "b"}, {}]],
      IntegerQ[#] || Head[#] === Rational &],
    Return[None, Module]];
  <|"schema" -> "diffexp2-exact-affine-jordan-indicial-v1",
    "dimension" -> dimension,
    "blocks" -> MapIndexed[Function[{block, index}, <|
        "block" -> First[index] - 1,
        "columns" -> (block["Cols"] - 1),
        "a" -> ToString[Cancel[Together[block["a"]]], InputForm],
        "b" -> ToString[Cancel[Together[block["b"]]], InputForm]|>],
      blocks]|>];

sccCapturedBlockRecord[parentSystemRecord_List, parentGeometry_Association,
    seq_Association, blockcs_Association, captured_Association,
    capabilities_Association, block_Integer,
    sourceTransform_, gaugeTransforms_] := Module[
  {vertices = seq["Components"][[block]], expectedPrincipal,
   analytic, principal, capturedCapabilities, capturedGeometry,
   capturedIdentity, exactIndicial},
  expectedPrincipal = parentSystemRecord[[vertices, vertices]];
  analytic = Lookup[captured["Metadata"], "ChartAnalytic", None];
  If[!AssociationQ[analytic],
    err["E5", blockcs, <|"Detail" ->
      "captured diagonal block is missing persistent chart analytic metadata"|>]];
  principal = Lookup[analytic, "principal_matrix", None];
  capturedCapabilities = Lookup[
    analytic, "native_scc_capabilities", None];
  capturedGeometry = Lookup[analytic, "geometry", None];
  If[!SameQ[principal, expectedPrincipal] ||
      !SameQ[capturedCapabilities, capabilities] ||
      !SameQ[capturedGeometry, parentGeometry],
    err["E6", blockcs, <|
      "Block" -> block, "ExpectedPrincipal" -> expectedPrincipal,
      "CapturedPrincipal" -> principal,
      "ExpectedCapabilities" -> capabilities,
      "CapturedCapabilities" -> capturedCapabilities,
      "Detail" -> "captured diagonal chart metadata does not bind its indexed parent principal block and geometry"|>]];
  If[!KeyExistsQ[captured["Metadata"], "ChartIdentity"],
    err["E6", blockcs, <|"Block" -> block,
      "Detail" -> "captured diagonal chart is missing its exact principal identity"|>]];
  capturedIdentity = ToString[
    captured["Metadata", "ChartIdentity"], InputForm];
  If[!StringQ[capturedIdentity] || StringLength[capturedIdentity] === 0,
    err["E6", blockcs, <|"Block" -> block,
      "Detail" -> "captured diagonal chart has no exact principal identity"|>]];
  exactIndicial = sccExactAffineJordanIndicialRecord[blockcs];
  Join[<|"block" -> block - 1, "vertices" -> (vertices - 1),
    "regular" -> capabilities["regular"],
    "identity_gauge" -> capabilities["identity_gauge"],
    "exact_gauge" -> capabilities["exact_gauge"],
    "identity_v" -> capabilities["identity_v"],
    "epsilon_unimodular_v" -> capabilities["epsilon_unimodular_v"],
    "no_pseudo" -> capabilities["no_pseudo"]|>,
    If[AssociationQ[sourceTransform],
      <|"source_transform" -> sourceTransform|>, <||>],
    If[AssociationQ[gaugeTransforms], gaugeTransforms, <||>],
    If[AssociationQ[exactIndicial],
      <|"exact_affine_jordan_indicial" -> exactIndicial|>, <||>]]];

sccNativeCompositeIdentity[parent_Association, blocks_List,
    couplings_List, domain_String, symbolNames_List] := Module[
  {couplingIdentities, identity},
  couplingIdentities = Map[KeyTake[#,
      {"source_block", "target_block", "source_vertices",
       "target_vertices", "rows", "columns", "exact_identity",
       "domain", "symbols"}] &, couplings];
  identity = Quiet[Check[ExportString[<|
      "schema" -> "diffexp2-native-scc-composite-v1",
      "parent" -> parent,
      "blocks" -> blocks,
      "couplings" -> couplingIdentities,
      "serialization" -> <|"domain" -> domain,
        "symbols" -> symbolNames|>|>,
    "RawJSON", "Compact" -> True], $Failed]];
  identity];

(* Domain-independent identity for an exact Rational SCC shadow and its Acb
   specialization.  Numeric coefficient encodings and session handles are
   deliberately excluded; the complete original/theta matrices, geometry,
   work rectangle, Jordan partitions, spectral frames, and coupling
   expressions remain proof obligations. *)
sccRationalShadowIdentity[parent_Association, blocks_List,
    couplings_List] := Module[{blockProofs, couplingProofs, payload, encoded},
  blockProofs = Map[<|
      "block" -> # ["block"], "vertices" -> # ["vertices"],
      "regular" -> # ["regular"],
      "identity_gauge" -> # ["identity_gauge"],
      "exact_gauge" -> # ["exact_gauge"],
      "identity_v" -> # ["identity_v"],
      "epsilon_unimodular_v" -> # ["epsilon_unimodular_v"],
      "no_pseudo" -> # ["no_pseudo"],
      "exact_affine_jordan_indicial" ->
        Lookup[#, "exact_affine_jordan_indicial", Null],
      "source_transform" -> If[AssociationQ[
          Lookup[#, "source_transform", None]],
        Join[KeyTake[# ["source_transform"], {"identity",
          "det_epsilon_valuation", "v_exact_identity",
          "vinv_exact_identity", "det_exact_identity"}],
          <|"entries" -> Map[<|"row" -> # ["row"],
              "column" -> # ["column"],
              "exact_entry" -> # ["exact_entry"],
              "epsilon_shift" -> # ["multiplier", "epsilon_shift"],
              "center_pole_order" ->
                # ["multiplier", "center_pole_order"]|> &,
            # ["source_transform", "entries"]]|>], Null],
      "to_physical" -> If[AssociationQ[Lookup[#, "to_physical", None]],
        Join[KeyTake[# ["to_physical"], {"role", "identity",
          "gauge_exact_identity", "gauge_inverse_exact_identity",
          "gauge_det_exact_identity"}], <|"entries" -> Map[
            <|"row" -> # ["row"], "column" -> # ["column"],
              "exact_entry" -> # ["exact_entry"],
              "epsilon_shift" -> # ["multiplier", "epsilon_shift"],
              "center_pole_order" -> # ["multiplier", "center_pole_order"]|> &,
            # ["to_physical", "entries"]]|>], Null],
      "to_reduced" -> If[AssociationQ[Lookup[#, "to_reduced", None]],
        Join[KeyTake[# ["to_reduced"], {"role", "identity",
          "gauge_exact_identity", "gauge_inverse_exact_identity",
          "gauge_det_exact_identity"}], <|"entries" -> Map[
            <|"row" -> # ["row"], "column" -> # ["column"],
              "exact_entry" -> # ["exact_entry"],
              "epsilon_shift" -> # ["multiplier", "epsilon_shift"],
              "center_pole_order" -> # ["multiplier", "center_pole_order"]|> &,
            # ["to_reduced", "entries"]]|>], Null]|> &, blocks];
  couplingProofs = Map[<|
      "source_block" -> # ["source_block"],
      "target_block" -> # ["target_block"],
      "source_vertices" -> # ["source_vertices"],
      "target_vertices" -> # ["target_vertices"],
      "rows" -> # ["rows"], "columns" -> # ["columns"],
      "entries" -> Map[<|"row" -> # ["row"],
          "column" -> # ["column"],
          "source_vertex" -> # ["source_vertex"],
          "target_vertex" -> # ["target_vertex"],
          "exact_original_entry" -> # ["exact_original_entry"],
          "exact_theta_entry" -> # ["exact_theta_entry"],
          "epsilon_shift" -> # ["multiplier", "epsilon_shift"],
          "center_pole_order" ->
            # ["multiplier", "center_pole_order"],
          "multiplier_exact_identity" ->
            # ["multiplier", "exact_identity"],
          "proven_zero" -> # ["multiplier", "proven_zero"]|> &,
        # ["entries"]]|> &, couplings];
  payload = <|"schema" -> "diffexp2-scc-rational-shadow-v1",
    "parent" -> KeyTake[parent, {"dimension", "exact_system_record",
      "exact_theta_record", "chart", "scc", "execution",
      "work_contract"}], "blocks" -> blockProofs,
    "couplings" -> couplingProofs|>;
  encoded = ExportString[payload, "RawJSON", "Compact" -> True];
  If[!StringQ[encoded], Return[$Failed, Module]];
  "de2-rational-shadow-" <>
    IntegerString[Hash[encoded, "SHA256"], 16, 64]];

(* The composite owner must retain one equation in the full original-master
   basis.  Diagonal chart payloads are intentionally not reused: their q/C
   records have block dimension and cannot certify a sourced propagation
   after the completed parent has been assembled.  The already collision-
   checked composite exact identity is also the physical owner signature. *)
sccParentPhysicalODEPayload[cs_Association, ownerIdentity_String,
    serialization_Association, inputDigits_Integer] := Module[
  {field, data, payload},
  If[StringLength[ownerIdentity] == 0 || inputDigits < 1,
    err["E6", cs, <|"Detail" ->
      "parent physical ODE capture requires a nonempty composite identity and positive input precision"|>]];
  field = sccSerializationField[serialization, cs];
  data = physicalClearedODEData[cs];
  payload = Block[{
      $cppSerializationDomain = field["domain"],
      $cppSerializationSymbols = field["symbols"]},
    cppPhysicalODEPayload[data, ownerIdentity, inputDigits, cs]];
  If[!AssociationQ[payload] ||
      Lookup[payload, "owner_signature_identity", None] =!=
        ownerIdentity ||
      !StringQ[Lookup[payload, "payload_identity", None]] ||
      Length[Lookup[payload, "c", {}]] < 1,
    err["E6", cs, <|"Detail" ->
      "parent physical q/C payload did not retain its full composite owner identity"|>]];
  payload];

(* Acb is the execution field, not the provenance field.  When the exact
   full-parent q/C equation is rational, retain a second canonical Rational
   encoding beside the ball payload.  Terminal adjoints and other structural
   consumers must not have to reconstruct exact zero/resonance facts from
   decimal Arb enclosures or depend on a separate Rational-shadow solve having
   happened earlier in the process. *)
sccParentRationalShadowPhysicalODEPayload[cs_Association,
    ownerIdentity_String, serialization_Association,
    inputDigits_Integer] := Module[{field, data, rationalCoefficientQ},
  field = sccSerializationField[serialization, cs];
  If[field["domain"] =!= "acb" || field["symbols"] =!= {},
    Return[None, Module]];
  data = physicalClearedODEData[cs];
  rationalCoefficientQ[entry_Association] :=
    TrueQ[Lookup[entry, "Zero", False]] ||
      (PolynomialQ[Lookup[entry, "P", $Failed],
          DiffExp2`Config`CanonicalEps[]] &&
       PolynomialQ[Lookup[entry, "Q", $Failed],
          DiffExp2`Config`CanonicalEps[]] &&
       AllTrue[Join[
          CoefficientList[entry["P"],
            DiffExp2`Config`CanonicalEps[]],
          CoefficientList[entry["Q"],
            DiffExp2`Config`CanonicalEps[]]],
         IntegerQ[#] || Head[#] === Rational &]);
  If[!AllTrue[Join[data["Q"], Flatten[data["C"]]],
      AssociationQ[#] && rationalCoefficientQ[#] &],
    Return[None, Module]];
  Block[{$cppSerializationDomain = "rational",
      $cppSerializationSymbols = {}},
    cppPhysicalODEPayload[data, ownerIdentity, inputDigits, cs]]];

(* The terminal singular tail theorem must be applied after the same exact
   block gauge which Fuchsianizes the diagonal SCCs.  Cross-block entries can
   still have negative epsilon valuation, although their condensation graph
   is acyclic.  A diagonal epsilon shearing g_i=eps^s_i h_i makes those
   couplings causal without multiplying the scalar q by a common epsilon
   monomial (which would destroy the formal-unit q(0,eps) premise).

   The difference constraints

                   s_source >= s_target - valuation(C_target,source)

   are solved as longest paths from a zero super-source.  A change on the
   d-th pass proves a positive cycle and is rejected: such a cycle would be a
   genuine two-sided epsilon recurrence, not the finite DAG reservoir used by
   this bridge. *)
sccRationalShadowSingularTailPayload[cs_Association, blockSystems_List,
    ownerIdentity_String, serialization_Association,
    inputDigits_Integer] := Module[
  {field, eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"],
   d = cs["SystemSize"], components, gauge, gaugeInverse, reducedTheta,
   reducedSystem, reducedData, positions, constraints, shifts, changed,
   candidate,
   shearedTheta, shearedSystem, shearedData, rationalCoefficientQ,
   q0, valuations, payload, fail},
  fail[detail_String] := Module[{},
    If[Environment["DE2_DIAGNOSTIC_SINGULAR_TAIL_FRAME"] === "1",
      Print["DE2 SINGULAR TAIL FRAME unavailable: ", detail]];
    None];
  field = sccSerializationField[serialization, cs];
  If[!MemberQ[{"rational", "acb"}, field["domain"]] ||
      field["symbols"] =!= {},
    Return[fail["execution field is not exact Rational-specializable"],
      Module]];
  components = Lookup[
    Lookup[cs, "IntegrationSequence", <||>], "Components", {}];
  If[!ListQ[components] || Length[components] =!= Length[blockSystems] ||
      Sort[Flatten[components]] =!= Range[d],
    err["E6", cs, <|"Detail" ->
      "singular-tail block gauge does not cover the full parent exactly"|>]];
  gauge = IdentityMatrix[d];
  gaugeInverse = IdentityMatrix[d];
  MapThread[Function[{vertices, block},
      gauge[[vertices, vertices]] = block["Gauge"];
      gaugeInverse[[vertices, vertices]] = block["GaugeInverse"]],
    {components, blockSystems}];
  reducedTheta = Map[Cancel[Together[#]] &,
    gaugeInverse . cs["ThetaOriginal"] . gauge -
      gaugeInverse . (t D[gauge, t]), {2}];
  reducedSystem = KeyDrop[
    Join[cs, <|"ThetaOriginal" -> reducedTheta|>], {"SystemClearKey"}];
  reducedData = physicalClearedODEData[reducedSystem];
  q0 = First[reducedData["Q"]];
  If[TrueQ[Lookup[q0, "Zero", False]] ||
      Lookup[q0, "Valuation", None] =!= 0,
    Return[fail["Fuchsian block gauge did not produce a formal-unit q0"],
      Module]];
  positions = Position[reducedData["C"],
    entry_Association /; !TrueQ[Lookup[entry, "Zero", False]],
    {3}];
  constraints = ({#[[2]], #[[3]],
        Lookup[Extract[reducedData["C"], #], "Valuation", None]} &) /@
    positions;
  If[!AllTrue[constraints,
      MatchQ[#, {_Integer, _Integer, _Integer}] &],
    Return[fail["reduced q/C valuations are malformed"], Module]];
  shifts = ConstantArray[0, d];
  Do[
    changed = False;
    Scan[Function[constraint,
      candidate = shifts[[constraint[[1]]]] - constraint[[3]];
      If[shifts[[constraint[[2]]]] < candidate,
        shifts[[constraint[[2]]]] = candidate;
        changed = True]], constraints];
    If[!changed, Break[]];
    If[iteration === d,
      Return[fail["negative epsilon valuations contain a positive cycle"],
        Module]],
    {iteration, d}];
  shifts -= Min[shifts];
  shearedTheta = MapIndexed[
    Cancel[Together[
      eps^(shifts[[#2[[2]]]] - shifts[[#2[[1]]]]) #1]] &,
    reducedTheta, {2}];
  shearedSystem = KeyDrop[
    Join[cs, <|"ThetaOriginal" -> shearedTheta|>], {"SystemClearKey"}];
  shearedData = physicalClearedODEData[shearedSystem];
  q0 = First[shearedData["Q"]];
  valuations = Join[
    Cases[shearedData["Q"],
      entry_Association /; !TrueQ[Lookup[entry, "Zero", False]] :>
        Lookup[entry, "Valuation", None]],
    Cases[Flatten[shearedData["C"]],
      entry_Association /; !TrueQ[Lookup[entry, "Zero", False]] :>
        Lookup[entry, "Valuation", None]]];
  If[TrueQ[Lookup[q0, "Zero", False]] ||
      Lookup[q0, "Valuation", None] =!= 0 ||
      !AllTrue[valuations, IntegerQ[#] && # >= 0 &],
    Return[fail[
      "epsilon shearing did not produce a causal formal-unit q/C equation"],
      Module]];
  rationalCoefficientQ[entry_Association] :=
    TrueQ[Lookup[entry, "Zero", False]] ||
      (PolynomialQ[Lookup[entry, "P", $Failed], eps] &&
       PolynomialQ[Lookup[entry, "Q", $Failed], eps] &&
       AllTrue[Join[CoefficientList[entry["P"], eps],
           CoefficientList[entry["Q"], eps]],
         IntegerQ[#] || Head[#] === Rational &]);
  If[!AllTrue[Join[shearedData["Q"], Flatten[shearedData["C"]]],
      AssociationQ[#] && rationalCoefficientQ[#] &],
    Return[fail["sheared q/C coefficients are not exact Rational"],
      Module]];
  payload = Block[{$cppSerializationDomain = "rational",
      $cppSerializationSymbols = {}},
    cppPhysicalODEPayload[
      shearedData, ownerIdentity, inputDigits, shearedSystem]];
  <|"schema" -> "diffexp2-scc-singular-tail-frame-v1",
    "equation" -> payload, "epsilon_shifts" -> shifts|>];

PrepareNativeSCCComposite[cs_Association, req_Association] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], epsWindow,
   requestedMin, requestedMax, publicTOrder, workTOrder, plannedTop,
   signature, cacheKey, cached, cacheStats, blockSystems,
   capabilities, badBlocks, captures, capturedContract, fb, width, workTop,
   parentRecords, parentGeometry, parent, blockRecords, serialization,
   couplings, identity, manifest, prepared, result, scale, radius,
   center, missingReq, components, condensation, executionDescriptor,
   inputDigits, runRecords, taskRecords, columnPlans, blockDimensions,
   framePlans, forcedFrame, physicalPayload, spectralFrames,
   sourceTransforms, gaugeTransforms, gaugeFrames, gaugePrepStart,
   gaugePrepMemory, gaugeProbeRecord, rationalShadowDecision,
   seedWorkHalos, blockRequiredTops, workReqs, reservoirMax,
   rationalShadowPhysicalPayload, rationalShadowTailPayload,
   diagnosticStartFrame, singularCompositeQ, compactFrameQ,
   coreFb, diagnosticLowerExtra, diagnosticStrictProbe,
   diagnosticTopHalo, singularTailTOrderHalo},
  missingReq = Select[{"TOrder", "EpsWindow"},
    !KeyExistsQ[req, #] &];
  epsWindow = Lookup[req, "EpsWindow", None];
  If[missingReq =!= {} || !AssociationQ[epsWindow] ||
      !AllTrue[{"Min", "CompleteMax"},
        KeyExistsQ[epsWindow, #] &] ||
      !IntegerQ[Lookup[req, "TOrder", None]] || req["TOrder"] < 0 ||
      !IntegerQ[Lookup[epsWindow, "Min", None]] ||
      !IntegerQ[Lookup[epsWindow, "CompleteMax", None]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"],
    err["E6", cs, <|"Request" -> req,
      "Detail" -> "native SCC preparation requires an ordered integer epsilon window and nonnegative Taylor order"|>]];
  components = If[AssociationQ[seq],
    Lookup[seq, "Components", {}], {}];
  If[!ListQ[components] || Length[components] < 1 ||
      (Length[components] === 1 &&
        TrueQ[Lookup[Lookup[cs, "IndicialData", <||>],
          "Regular", False]]),
    err["E6", cs, <|"Detail" ->
      "native SCC preparation requires a multi-block exact SCC certificate or a singular one-block certificate"|>]];
  If[cfg["RecurrenceBackend"] =!= "Cpp" ||
      TrueQ[$disableGroupedSpectralTransform] ||
      !TrueQ[$cppUsePersistentSessions] ||
      Environment["DE2_CPP_PERSISTENT"] === "0",
    err["E6", cs, <|"Detail" ->
      "native SCC preparation requires the persistent grouped C++ recurrence backend"|>]];
  If[DownValues[DiffExp2`CppBackend`PreparePersistentSCC] === {} ||
      DownValues[DiffExp2`CppBackend`PersistentSCCStatistics] === {},
    err["E5", cs, <|"Detail" ->
      "CppBackend persistent SCC prepare/statistics bridge is not available"|>]];
  center = Lookup[cs, "Center", None];
  scale = Lookup[Lookup[cs, "ChartMap", <||>], "Scale", None];
  radius = Lookup[cs, "Radius", None];
  If[!FreeQ[{center, scale, radius}, _?InexactNumberQ] ||
      !cppExactRealAlgebraicQ[center] ||
      !cppExactRealAlgebraicQ[scale] ||
      !cppExactAlgebraicTruthQ[scale != 0] ||
      radius === Infinity ||
      !cppExactRealAlgebraicQ[radius] ||
      !cppExactAlgebraicTruthQ[radius > 0],
    err["E6", cs, <|"Scale" -> scale, "Radius" -> radius,
      "Detail" -> "native SCC first slice requires exact real algebraic geometry, a nonzero scale, and a positive finite radius"|>]];

  signature = sccNativeCompositeCacheSignature[cs, req];
  cacheKey = Hash[signature, "SHA256"];
  cached = Lookup[$nativeSCCCompositeCache, cacheKey, None];
  If[AssociationQ[cached] && SameQ[cached["Signature"], signature],
    If[!sccNativeCompositeExecutionDescriptorQ[cached, signature],
      err["E6", cs, <|"CacheKey" -> cacheKey,
        "Detail" -> "native SCC cached execution descriptor is not bound to its complete composite signature and public handle"|>]];
    cacheStats = DiffExp2`CppBackend`PersistentSCCStatistics[
      cached["Result"]];
    If[sccNativeCompositeStatisticsQ[cacheStats, cached["Result"]],
      cached["Execution"] = Join[cached["Execution"],
        <|"NativeStatistics" -> cacheStats|>];
      AssociateTo[$nativeSCCCompositeCache, cacheKey -> cached];
      Return[cached["Result"], Module]];
    (* Direct native release/session clearing can invalidate an otherwise
       exact Solve cache entry.  Drop only that proved-stale entry and rebuild
       it; never evict a live public handle as a cache policy. *)
    KeyDropFrom[$nativeSCCCompositeCache, cacheKey];
    cached = None];
  If[cached =!= None,
    err["E6", cs, <|"CacheKey" -> cacheKey,
      "Detail" -> "native SCC composite cache key collided with an unequal full signature"|>]];
  nativeSCCCompositeRequireAdmission[cs];

  condensation = seq["CondensationEdges"];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC BLOCK SYSTEMS start center=", center,
      " memory=", MemoryInUse[]]];
  blockSystems = Map[Function[block,
      If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
        Print["DE2 SCC BLOCK SYSTEM start center=", center,
          " block=", block, " indices=", components[[block]],
          " memory=", MemoryInUse[]]];
      Module[{system = sccBlockChartSystem[cs, block]},
        If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
          Print["DE2 SCC BLOCK SYSTEM done center=", center,
            " block=", block,
            " regular=", Lookup[
              Lookup[system, "IndicialData", <||>], "Regular", None],
            " memory=", MemoryInUse[]]];
        system]],
    Range[Length[components]]];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC BLOCK SYSTEMS done center=", center,
      " memory=", MemoryInUse[]]];
  spectralFrames = sccSpectralFrameCertificate /@ blockSystems;
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC SPECTRAL FRAMES done center=", center,
      " memory=", MemoryInUse[]]];
  gaugeFrames = sccGaugeFrameCertificate /@ blockSystems;
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC GAUGE FRAMES done center=", center,
      " identity=", Lookup[gaugeFrames, "identity_gauge", None],
      " memory=", MemoryInUse[]]];
  capabilities = MapThread[sccNativeBlockCapabilities,
    {blockSystems, spectralFrames, gaugeFrames}];
  badBlocks = Select[Range[Length[blockSystems]],
    Function[block, Module[{record = capabilities[[block]]},
      (* "regular" is a classification, not an admission predicate.  The
         retained C++ chart proves the complete exact affine-Jordan
         indicial operator for a singular block and revalidates every T/P/R
         schedule at execution.  Wolfram must therefore admit both Boolean
         classes here while continuing to require the exact producer facts
         which the current composite representation actually depends on. *)
      !MemberQ[{True, False}, Lookup[record, "regular", None]] ||
        !TrueQ[Lookup[record, "exact_gauge", False]] ||
        !TrueQ[Lookup[record, "epsilon_unimodular_v", False]] ||
        !MemberQ[{True, False},
          Lookup[record, "no_pseudo", None]]]]];
  If[badBlocks =!= {},
    err["E6", cs, <|"UnsupportedBlocks" -> Map[
        <|"Block" -> #, "Capabilities" -> capabilities[[#]],
          "SpectralFrame" -> spectralFrames[[#]]|> &,
        badBlocks],
      "Detail" -> "native SCC preparation requires regular or exact affine-Jordan diagonal blocks with exact invertible Gauge/GaugeInverse and an exact t-independent epsilon-unimodular V/VInv pair; no_pseudo is retained provenance, not an admission decision"|>]];
  requestedMin = epsWindow["Min"];
  reservoirMax = epsWindow["CompleteMax"];
  requestedMax = Lookup[req, "RequiredCompleteMax", reservoirMax];
  If[!IntegerQ[requestedMax] ||
      !(requestedMin <= requestedMax <= reservoirMax),
    err["E6", cs, <|"RequestedWindow" -> epsWindow,
      "RequiredCompleteMax" -> requestedMax,
      "Detail" -> "native SCC required epsilon edge must lie inside its private reservoir"|>]];
  singularCompositeQ = !TrueQ[$disableAdaptiveLowerFrames] &&
    AnyTrue[blockSystems,
      !TrueQ[Lookup[Lookup[#, "IndicialData", <||>],
        "Regular", False]] &];
  publicTOrder = req["TOrder"];
  (* Preserve an explicitly requested private singular-tail reservoir, but
     do not impose it on every singular chart.  Terminal matching decides
     whether its Frobenius bridge actually contracts before a bounded
     Taylor-order retry asks for more rows.  Structurally noncontracting
     charts must keep the established exact physical fallback at public
     order. *)
  singularTailTOrderHalo =
    Lookup[req, "SingularTailTOrderHalo", 0];
  If[!IntegerQ[singularTailTOrderHalo] ||
      singularTailTOrderHalo < 0,
    err["E6", cs, <|
      "SingularTailTOrderHalo" -> singularTailTOrderHalo,
      "Detail" ->
        "native SCC singular-tail Taylor halo must be a nonnegative integer"|>]];
  If[!TrueQ[singularCompositeQ],
    singularTailTOrderHalo = 0];
  workTOrder = sccWorkTOrder[cs, req] + singularTailTOrderHalo;
  seedWorkHalos = sccSeedWorkHalos[cs, blockSystems, workTOrder];
  blockRequiredTops = reservoirMax + seedWorkHalos;
  diagnosticStartFrame =
    Environment["DE2_SCC_DIAGNOSTIC_START_FRAME"] === "1";
  compactFrameQ = diagnosticStartFrame || singularCompositeQ;
  diagnosticStrictProbe =
    Environment["DE2_SCC_DIAGNOSTIC_STRICT_PROBE"] === "1";
  diagnosticLowerExtra = Switch[
    Environment["DE2_SCC_DIAGNOSTIC_LOWER_EXTRA"],
    "1", 1, "2", 2, "4", 4, "8", 8, "16", 16,
    _, If[diagnosticStartFrame, 0,
      sccNativeEffectiveFrameLowerExtra[cs]]];
  diagnosticTopHalo = Switch[
    Environment["DE2_SCC_DIAGNOSTIC_TOP_HALO"],
    "4", 4, "8", 8, "16", 16, "32", 32, "64", 64,
    _, If[diagnosticStartFrame, None,
      If[singularCompositeQ, sccNativeEffectiveFrameTopHalo[], None]]];
  If[IntegerQ[diagnosticTopHalo],
    seedWorkHalos = Min[#, diagnosticTopHalo] & /@ seedWorkHalos;
    blockRequiredTops = reservoirMax + seedWorkHalos];
  plannedTop = Max[blockRequiredTops];
  workReqs = Map[Function[requiredTop, Join[req, <|
      "EpsWindow" -> Join[epsWindow, <|"CompleteMax" -> requiredTop|>],
      "TOrder" -> workTOrder|>]], blockRequiredTops];
  framePlans = MapThread[homogeneousFramePlan, {blockSystems, workReqs}];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC FRAME PLAN center=", center,
      " requested=", {requestedMin, requestedMax},
      " reservoirMax=", reservoirMax,
      " plannedTop=", plannedTop,
      " seedWorkHalos=", seedWorkHalos,
      " blockRequiredTops=", blockRequiredTops,
      " workTOrder=", workTOrder,
      " blocks=", KeyTake[#, {"PMax", "CollisionDepthMax",
          "PoleDepth", "SingleUseDepth", "SpectralDepth",
          "TransformDepth", "EpsilonRegularPrincipal",
          "InitialSourceDepth", "FrameTop", "TerminalFrameBase"}] & /@
        framePlans]];
  coreFb = Min[Lookup[framePlans, "StartFrameBase"]];
  fb = If[compactFrameQ, coreFb,
    Min[Lookup[framePlans, "TerminalFrameBase"]]];
  If[compactFrameQ,
    fb = Max[Min[Lookup[framePlans, "TerminalFrameBase"]],
      fb - diagnosticLowerExtra]];
  workTop = Max[Lookup[framePlans, "FrameTop"]];
  If[IntegerQ[diagnosticTopHalo],
    workTop = Max[plannedTop,
      Max[Map[# ["FrameTop"] - # ["PoleDepth"] &,
        framePlans]]];
    framePlans = Map[Join[#, <|"FrameTop" -> workTop|>] &,
      framePlans]];
  forcedFrame = <|"FrameBase" -> fb, "FrameTop" -> workTop|>;
  captures = Block[{$shCache = <||>, $shSysTag = None,
      $cppStaticOperatorCache = <||>,
      $cancellationAuditBase = If[
        compactFrameQ && fb < coreFb, coreFb, None],
      $adaptiveLowerFrameProbe = compactFrameQ &&
        (singularCompositeQ || diagnosticStrictProbe)},
    MapThread[sccCaptureHomogeneousGroup[
        #1, #2, #3, forcedFrame] &,
      {blockSystems, workReqs, framePlans}]];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC CAPTURES done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  capturedContract = sccCapturedCompositeContract[
    captures, cs, workTOrder];
  If[capturedContract["FrameBase"] =!= fb ||
      capturedContract["FrameWidth"] =!= workTop - fb + 1,
    err["E6", cs, <|"PlannedFrame" -> forcedFrame,
      "CapturedFrame" -> {capturedContract["FrameBase"],
        capturedContract["FrameWidth"]},
      "Detail" -> "native SCC capture did not preserve the exact union work rectangle"|>]];
  width = capturedContract["FrameWidth"];
  If[!(fb <= requestedMin <= requestedMax <= workTop),
    err["E6", cs, <|"FrameBase" -> fb, "FrameWidth" -> width,
      "RequestedWindow" -> epsWindow,
      "Detail" -> "captured native SCC work frame does not contain the requested epsilon window"|>]];

  parentRecords = sccParentExactRecords[cs];
  parentGeometry = cppPersistentGeometry[cs];
  parent = Join[<|"dimension" -> cs["SystemSize"]|>,
    parentRecords, <|
      "chart" -> parentGeometry,
      "scc" -> cppPersistentSCC[cs],
      "execution" -> <|"mode" -> "BlockSequentialStrict",
        "work_t_order" -> workTOrder|>,
      "work_contract" -> <|"work_min" -> fb,
        "cancellation_audit_base" ->
          capturedContract["CancellationAuditBase"],
        "requested_min" -> requestedMin,
        "requested_max" -> requestedMax,
        "work_complete_max" -> workTop,
        "public_t_order" -> publicTOrder,
        "singular_tail_t_order_halo" ->
          singularTailTOrderHalo,
        "wolfram_coupling_depth" -> seq["CouplingDepth"]|>|>];
  serialization = capturedContract["Serialization"];
  sourceTransforms = MapThread[
    If[TrueQ[Lookup[#3, "identity_v", False]], None,
      PrepareSCCSpectralSourceTransform[#1,
        sccNativeSourceShape[cs, #2, fb, workTop, workTOrder],
        serialization]] &,
    {blockSystems, Length /@ components, spectralFrames, framePlans}];
  gaugePrepStart = AbsoluteTime[];
  gaugePrepMemory = MemoryInUse[];
  gaugeTransforms = MapThread[
    If[TrueQ[Lookup[#3, "identity_gauge", False]], None,
      <|"to_physical" -> PrepareSCCGaugeTransform[#1, "to_physical",
          sccNativeSourceShape[cs, #2, fb, workTop, workTOrder],
          serialization, #4],
        "to_reduced" -> PrepareSCCGaugeTransform[#1, "to_reduced",
          sccNativeSourceShape[cs, #2, fb, workTop, workTOrder],
          serialization, #4]|>] &,
    {blockSystems, Length /@ components, capabilities, gaugeFrames}];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC TRANSFORMS done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  If[Environment["DE2_SCC_GAUGE_PREP_PROBE"] === "1" &&
      AnyTrue[gaugeTransforms, AssociationQ],
    gaugeProbeRecord = <|
      "schema" -> "diffexp2-scc-gauge-prep-probe-v1",
      "elapsed_seconds" -> N[AbsoluteTime[] - gaugePrepStart],
      "memory_before" -> gaugePrepMemory,
      "memory_after" -> MemoryInUse[],
      "max_memory_used" -> MaxMemoryUsed[],
      "certified" -> AllTrue[gaugeFrames,
        TrueQ[Lookup[#, "admissible", False]] &],
      "transforms" -> Map[If[AssociationQ[#],
          AssociationMap[<|"entries" -> Length[Lookup[#, "entries", {}]],
              "kernel_shapes" -> DeleteDuplicates[
                Dimensions[Lookup[Lookup[#, "multiplier", <||>],
                    "kernels", {}]] & /@ Lookup[#, "entries", {}]]|> &,
            #], <||>] &, gaugeTransforms]|>;
    Print["DE2 SCC GAUGE PREP PROBE ", InputForm[gaugeProbeRecord]];
    err["E6", cs, <|"GaugePrepProbe" -> gaugeProbeRecord,
      "Detail" -> "requested SCC gauge preparation probe completed before block capture"|>]];
  blockRecords = MapThread[
    sccCapturedBlockRecord[parentRecords["exact_system_record"],
      parentGeometry, seq, #1, #2, #3, #4, #5, #6] &,
    {blockSystems, captures, capabilities,
      Range[Length[blockSystems]], sourceTransforms, gaugeTransforms}];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC BLOCK RECORDS done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  couplings = Map[
    PrepareSCCCouplingMatrix[cs, #[[1]], #[[2]],
      sccNativeSourceShape[cs, Length[components[[#[[1]]]]],
      fb, workTop, workTOrder], serialization] &,
    condensation];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC COUPLINGS done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  identity = sccNativeCompositeIdentity[parent, blockRecords, couplings,
    capturedContract["Domain"], capturedContract["SymbolNames"]];
  If[!StringQ[identity],
    err["E6", cs, <|"Detail" ->
      "native SCC exact parent identity could not be serialized"|>]];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*
    cfg["WorkingPrecision"];
  physicalPayload = sccParentPhysicalODEPayload[
    cs, identity, serialization, inputDigits];
  rationalShadowPhysicalPayload =
    sccParentRationalShadowPhysicalODEPayload[
      cs, identity, serialization, inputDigits];
  rationalShadowTailPayload = If[
    TrueQ[singularCompositeQ] &&
      (serialization["domain"] === "rational" ||
       AssociationQ[rationalShadowPhysicalPayload]),
    sccRationalShadowSingularTailPayload[
      cs, blockSystems, identity, serialization, inputDigits],
    None];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC PHYSICAL PAYLOAD done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  manifest = <|"identity" -> identity,
    "rational_shadow_identity" ->
      sccRationalShadowIdentity[parent, blockRecords, couplings],
    "parent" -> parent,
    "blocks" -> blockRecords, "couplings" -> couplings,
    "physical_ode" -> physicalPayload|>;
  If[AssociationQ[rationalShadowPhysicalPayload],
    AssociateTo[manifest,
      "rational_shadow_physical_ode" ->
        rationalShadowPhysicalPayload]];
  If[AssociationQ[rationalShadowTailPayload],
    AssociateTo[manifest,
      "rational_shadow_singular_tail" ->
        rationalShadowTailPayload]];
  runRecords = Map[Function[capture,
      KeyTake[#, $nativeSCCColumnRunKeys] & /@ capture["Requests"]],
    captures];
  If[!AllTrue[Flatten[runRecords, 1], AssociationQ[#] &&
        Sort[Keys[#]] === Sort[$nativeSCCColumnRunKeys] &],
    err["E6", cs, <|"Detail" ->
      "captured native SCC run records are incomplete after static-payload compaction"|>]];
  taskRecords = Lookup[captures, "TaskMetadata", None];
  If[!ListQ[taskRecords] || Length[taskRecords] =!= Length[captures] ||
      !(And @@ MapThread[Length[#1] === Length[#2] &,
        {taskRecords, runRecords}]) ||
      !AllTrue[Flatten[taskRecords, 1], AssociationQ[#] &&
        Sort[Keys[#]] === Sort[{"a", "b", "P"}] &],
    err["E6", cs, <|"Detail" ->
      "captured native SCC exact task metadata is incomplete after static-payload compaction"|>]];
  blockDimensions = Lookup[blockSystems, "SystemSize", None];
  columnPlans = If[MemberQ[{"rational", "acb"},
        capturedContract["Domain"]] &&
      capturedContract["SymbolNames"] === {},
    sccNativeCapturedColumnPlans[cs, blockSystems,
      captures, capturedContract, inputDigits],
    (* Symbolic composites remain valid prepared objects.  Rational and Acb
       plans above are derived only from exact task/indicial metadata; Acb
       encodings are compared solely with independently encoded exact
       zero/one values and never inspected through a midpoint. *)
    ConstantArray[{}, Length[blockSystems]]];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC COLUMN PLANS done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  rationalShadowDecision = sccNativeRationalShadowDecision[
    capturedContract["Domain"], blockRecords, couplings, columnPlans];
  prepared = DiffExp2`CppBackend`PreparePersistentSCC[captures, manifest];
  If[Environment["DE2_SCC_GAUGE_TIMING"] === "1",
    Print["DE2 SCC NATIVE PREPARE done center=", center,
      " t=", AbsoluteTime[], " memory=", MemoryInUse[]]];
  If[FailureQ[prepared] || !AssociationQ[prepared] ||
      !StringQ[Lookup[prepared, "Session", None]] ||
      !StringQ[Lookup[prepared, "SCC", None]] ||
      !StringQ[Lookup[prepared, "Key", None]],
    err["E5", cs, <|"BackendFailure" -> prepared,
      "Detail" -> "persistent native SCC preparation failed"|>]];
  result = prepared;
  cacheStats = DiffExp2`CppBackend`PersistentSCCStatistics[result];
  If[FailureQ[cacheStats] ||
      !sccNativeCompositeStatisticsQ[cacheStats, result],
    Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[result]];
    err["E5", cs, <|"BackendFailure" -> cacheStats,
      "Detail" -> "new persistent native SCC could not be inspected before caching"|>]];
  executionDescriptor = <|
    "PublicHandle" -> result, "BlockDimensions" -> blockDimensions,
    "Runs" -> runRecords, "TaskMetadata" -> taskRecords,
    "ColumnPlans" -> columnPlans,
    "RationalShadowDecision" -> rationalShadowDecision,
    "Contract" -> capturedContract,
    "ParentIdentity" -> identity, "Components" -> components,
    "CondensationEdges" -> condensation,
    "TopologicalOrder" -> seq["TopologicalOrder"],
    "InputDigits" -> inputDigits, "NativeStatistics" -> cacheStats|>;
  AssociateTo[$nativeSCCCompositeCache, cacheKey -> <|
    "Signature" -> signature, "Result" -> result,
    "Execution" -> executionDescriptor|>];
  result];

sccNativeCanonicalScalarScheduleQ[schedule_, nmax_Integer] :=
  ListQ[schedule] && Length[schedule] === nmax + 1 &&
    And @@ MapIndexed[Function[{row, index}, Module[
      {n = First[index] - 1, step},
      If[!ListQ[row] || Length[row] =!= 1 ||
          !AssociationQ[First[row]], Return[False, Module]];
      step = First[row];
      Sort[Keys[step]] === Sort[{"case", "da", "db"}] &&
        Lookup[step, "case", None] === If[n === 0, "R", "T"] &&
        Lookup[step, "da", None] === ToString[n, InputForm] &&
        Lookup[step, "db", None] === "0"]], schedule];

sccNativeCanonicalRegularBlockScheduleQ[schedule_, nmax_Integer,
    dimension_Integer] :=
  dimension > 0 && ListQ[schedule] && Length[schedule] === nmax + 1 &&
    And @@ MapIndexed[Function[{row, index}, Module[
      {n = First[index] - 1},
      ListQ[row] && Length[row] === dimension &&
        AllTrue[row, Function[step,
          AssociationQ[step] &&
            Sort[Keys[step]] === Sort[{"case", "da", "db"}] &&
            Lookup[step, "case", None] === If[n === 0, "R", "T"] &&
            Lookup[step, "da", None] === ToString[n, InputForm] &&
            Lookup[step, "db", None] === "0"]]]], schedule];

(* Certify the exact request that the ordinary homogeneous builder captured;
   execution never fabricates a seed schedule independently of that builder.
   This v1 domain is intentionally only the exact eps^0 scalar family. *)
sccNativeCanonicalScalarHomogeneousRequestQ[request_Association,
    contract_Association] := Module[
  {fb = Lookup[contract, "FrameBase", None],
   width = Lookup[contract, "FrameWidth", None],
   nmax = Lookup[contract, "NMax", None], frameTop, unitIndex, initial},
  If[!MatchQ[{fb, width, nmax}, {_Integer, _Integer, _Integer}] ||
      width < 1 || nmax < 0 ||
      Sort[Keys[request]] =!= Sort[$nativeSCCColumnRunKeys],
    Return[False, Module]];
  frameTop = fb + width - 1;
  unitIndex = 1 - fb;
  initial = Lookup[request, "initial", None];
  Lookup[request, "nmax", None] === nmax &&
    Lookup[request, "p", None] === 0 &&
    TrueQ[Lookup[request, "has_initial", False]] &&
    TrueQ[!Lookup[request, "adaptive_probe", True]] &&
    Lookup[request, "a_target", None] === "0" &&
    Lookup[request, "b_target", None] === "0" &&
    Lookup[request, "a_shift_min", None] === 0 &&
    Lookup[request, "a_shifts", None] ===
      (ToString[#, InputForm] & /@ Range[0, nmax]) &&
    sccNativeCanonicalScalarScheduleQ[
      Lookup[request, "schedule", None], nmax] &&
    ListQ[initial] && Length[initial] === width &&
    1 <= unitIndex <= width &&
    And @@ MapIndexed[
      If[First[#2] === unitIndex, #1 === "1", #1 === "0"] &,
      initial] &&
    Lookup[request, "initial_validity", None] === {frameTop} &&
    Lookup[request, "source", Missing["Source"]] === Null &&
    TrueQ[!Lookup[request, "return_u", True]]];

(* A multidimensional regular block is captured as one ordinary homogeneous
   request per local basis component.  Derive the component from the honest
   flattened initial tensor instead of trusting request order.  Returning a
   one-based component (or None) keeps the selection proof separate from the
   later C++ provenance check. *)
sccNativeCanonicalRegularBlockComponent[request_Association,
    contract_Association, dimension_Integer] := Module[
  {fb = Lookup[contract, "FrameBase", None],
   width = Lookup[contract, "FrameWidth", None],
   nmax = Lookup[contract, "NMax", None], frameTop, unitIndex, initial,
   validity, zeroFrame, unitFrame, frames, selected},
  If[dimension < 1 ||
      !MatchQ[{fb, width, nmax}, {_Integer, _Integer, _Integer}] ||
      width < 1 || nmax < 0 ||
      Sort[Keys[request]] =!= Sort[$nativeSCCColumnRunKeys],
    Return[None, Module]];
  frameTop = fb + width - 1;
  unitIndex = 1 - fb;
  initial = Lookup[request, "initial", None];
  validity = Lookup[request, "initial_validity", None];
  If[Lookup[request, "nmax", None] =!= nmax ||
      Lookup[request, "p", None] =!= 0 ||
      !TrueQ[Lookup[request, "has_initial", False]] ||
      !TrueQ[!Lookup[request, "adaptive_probe", True]] ||
      Lookup[request, "a_target", None] =!= "0" ||
      Lookup[request, "b_target", None] =!= "0" ||
      Lookup[request, "a_shift_min", None] =!= 0 ||
      Lookup[request, "a_shifts", None] =!=
        (ToString[#, InputForm] & /@ Range[0, nmax]) ||
      !sccNativeCanonicalRegularBlockScheduleQ[
        Lookup[request, "schedule", None], nmax, dimension] ||
      !ListQ[initial] || Length[initial] =!= dimension width ||
      !ListQ[validity] || validity =!= ConstantArray[frameTop, dimension] ||
      !(1 <= unitIndex <= width) ||
      Lookup[request, "source", Missing["Source"]] =!= Null ||
      !TrueQ[!Lookup[request, "return_u", True]],
    Return[None, Module]];
  zeroFrame = ConstantArray["0", width];
  unitFrame = ReplacePart[zeroFrame, unitIndex -> "1"];
  frames = Partition[initial, width];
  If[!AllTrue[frames, # === zeroFrame || # === unitFrame &],
    Return[None, Module]];
  selected = Select[Range[dimension], frames[[#]] === unitFrame &];
  If[Length[selected] === 1, First[selected], None]];

(* Acb recurrence payloads encode even exact integers as {real,imag} pairs.
   Certify the same canonical regular unit-column shape without interpreting
   an interval midpoint as an exact value: the encodings below are generated
   independently from exact integers in the selected coefficient domain, and
   the native executor rechecks them as exact Arb objects. *)
sccNativeCanonicalEncodedRegularBlockComponent[
    request_Association, task_Association, contract_Association,
    dimension_Integer, zero_, one_, encodedATarget_, encodedBTarget_,
    encodedAShifts_List] := Module[
  {fb = Lookup[contract, "FrameBase", None],
   width = Lookup[contract, "FrameWidth", None],
   nmax = Lookup[contract, "NMax", None], frameTop, unitIndex,
   initial, validity, schedule, zeroFrame, unitFrame, frames, selected},
  If[dimension < 1 ||
      !MatchQ[{fb, width, nmax}, {_Integer, _Integer, _Integer}] ||
      width < 1 || nmax < 0 || Length[encodedAShifts] =!= nmax + 1 ||
      Sort[Keys[request]] =!= Sort[$nativeSCCColumnRunKeys] ||
      Sort[Keys[task]] =!= Sort[{"a", "b", "P"}] ||
      !FreeQ[{task["a"], task["b"]}, _?InexactNumberQ] ||
      task["P"] =!= 0,
    Return[None, Module]];
  frameTop = fb + width - 1;
  unitIndex = 1 - fb;
  initial = Lookup[request, "initial", None];
  validity = Lookup[request, "initial_validity", None];
  schedule = Lookup[request, "schedule", None];
  If[Lookup[request, "nmax", None] =!= nmax ||
      Lookup[request, "p", None] =!= 0 ||
      !TrueQ[Lookup[request, "has_initial", False]] ||
      TrueQ[Lookup[request, "adaptive_probe", True]] ||
      Lookup[request, "a_target", None] =!= encodedATarget ||
      Lookup[request, "b_target", None] =!= encodedBTarget ||
      Lookup[request, "a_shift_min", None] =!= 0 ||
      Lookup[request, "a_shifts", None] =!= encodedAShifts ||
      !ListQ[schedule] || Length[schedule] =!= nmax + 1 ||
      !AllTrue[schedule, ListQ[#] && Length[#] === dimension &&
          AllTrue[#, AssociationQ[#] &&
            Sort[Keys[#]] === Sort[{"case", "da", "db"}] &&
            MemberQ[{"R", "T", "P"}, Lookup[#, "case", None]] &&
            MatchQ[Lookup[#, "da", None], {_String, _String}] &&
            MatchQ[Lookup[#, "db", None], {_String, _String}] &] &] ||
      !ListQ[initial] || Length[initial] =!= dimension width ||
      !ListQ[validity] ||
      validity =!= ConstantArray[frameTop, dimension] ||
      !(1 <= unitIndex <= width) ||
      Lookup[request, "source", Missing["Source"]] =!= Null ||
      TrueQ[Lookup[request, "return_u", True]],
    Return[None, Module]];
  zeroFrame = ConstantArray[zero, width];
  unitFrame = ReplacePart[zeroFrame, unitIndex -> one];
  frames = Partition[initial, width];
  If[!AllTrue[frames, # === zeroFrame || # === unitFrame &],
    Return[None, Module]];
  selected = Select[Range[dimension], frames[[#]] === unitFrame &];
  If[Length[selected] === 1, First[selected], None]];

(* A singular Jordan seed has the same unique eps^0 unit in its log-zero
   row, followed by the canonical eps-shifted upper-chain ladder.  The C++
   certificate verifies that complete ladder; Wolfram uses only the exact
   unit position to map a captured request to its local basis component. *)
sccNativeCapturedJordanComponent[request_Association,
    task_Association, contract_Association, dimension_Integer,
    zero_, one_] := Module[
  {fb = Lookup[contract, "FrameBase", None],
   width = Lookup[contract, "FrameWidth", None],
   nmax = Lookup[contract, "NMax", None], frameTop, unitIndex, p,
   initial, validity, frames, logZeroFrames, zeroFrame, unitFrame,
   selected},
  p = Lookup[request, "p", None];
  If[dimension < 1 ||
      !MatchQ[{fb, width, nmax, p},
        {_Integer, _Integer, _Integer, _Integer}] ||
      width < 1 || nmax < 0 || p < 0 ||
      Sort[Keys[request]] =!= Sort[$nativeSCCColumnRunKeys] ||
      Sort[Keys[task]] =!= Sort[{"a", "b", "P"}] ||
      Lookup[task, "P", None] =!= p,
    Return[None, Module]];
  frameTop = fb + width - 1;
  unitIndex = 1 - fb;
  initial = Lookup[request, "initial", None];
  validity = Lookup[request, "initial_validity", None];
  If[Lookup[request, "nmax", None] =!= nmax ||
      !TrueQ[Lookup[request, "has_initial", False]] ||
      Lookup[request, "a_shift_min", None] =!= 0 ||
      !ListQ[initial] ||
      Length[initial] =!= (p + 1) dimension width ||
      !ListQ[validity] ||
      validity =!= ConstantArray[frameTop, (p + 1) dimension] ||
      !(1 <= unitIndex <= width) ||
      Lookup[request, "source", Missing["Source"]] =!= Null ||
      TrueQ[Lookup[request, "return_u", True]],
    Return[None, Module]];
  frames = Partition[initial, width];
  logZeroFrames = Take[frames, dimension];
  zeroFrame = ConstantArray[zero, width];
  unitFrame = ReplacePart[zeroFrame, unitIndex -> one];
  If[!AllTrue[logZeroFrames,
      # === zeroFrame || # === unitFrame &],
    Return[None, Module]];
  selected = Select[Range[dimension],
    logZeroFrames[[#]] === unitFrame &];
  If[Length[selected] === 1, First[selected], None]];

sccNativeParticularRunForTag[blockcs_Association,
    contract_Association, a_, b_, sourceP_Integer,
    inputDigits_Integer, domain_String] := Module[
  {dimension = blockcs["SystemSize"], blocks = blockList[blockcs],
   nmax = contract["NMax"], width = contract["FrameWidth"], p,
   encode, zero, schedule},
  If[sourceP < 0,
    err["E6", blockcs, <|"SourceLogPower" -> sourceP,
      "Detail" -> "native SCC source log power must be nonnegative"|>]];
  p = logCeiling[blockcs, a, b, sourceP, True];
  If[!MemberQ[{"rational", "acb"}, domain],
    err["E6", blockcs, <|"Domain" -> domain,
      "Detail" -> "native singular SCC column plans require a Rational or Acb coefficient field"|>]];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    encode[value_] := cppScalar[Cancel[Together[value]],
      inputDigits, blockcs];
    zero = encode[0];
    schedule = Table[Map[Function[blk, Module[{dA, dB, kind},
          dA = Cancel[Together[a + n - blk["a"]]];
          dB = Cancel[Together[b - blk["b"]]];
          kind = Which[!zeroCanQ[dA], "T", !zeroCanQ[dB], "P",
            True, "R"];
          <|"case" -> kind, "da" -> encode[dA],
            "db" -> encode[dB]|>]], blocks],
      {n, 0, nmax}];
    <|"nmax" -> nmax, "p" -> p, "has_initial" -> False,
      "adaptive_probe" -> True,
      "cancellation_audit_base" -> Lookup[
        contract, "CancellationAuditBase", Null],
      "a_target" -> encode[a], "b_target" -> encode[b],
      "a_shift_min" -> 0,
      "a_shifts" -> Table[encode[a + n], {n, 0, nmax}],
      "schedule" -> schedule,
      "initial" -> ConstantArray[zero, (p + 1) dimension width],
      "initial_validity" -> ConstantArray[Null, (p + 1) dimension],
      "source" -> Null, "return_u" -> False|>]];

(* Build the complete cheap run plan once per captured basis column.  Along
   a DAG the affine sector tag is preserved by the admitted couplings, while
   the required log ceiling can grow at each resonant target.  Propagate the
   exact ceiling in topological order and retain the resulting requests. *)
sccNativeCapturedColumnPlans[cs_Association, blockSystems_List,
    captures_List, contract_Association, inputDigits_Integer] := Module[
  {seq, dimensions, requests, tasks, zero, one, plans, topological,
   edges, domain},
  seq = cs["IntegrationSequence"];
  topological = seq["TopologicalOrder"];
  edges = seq["CondensationEdges"];
  dimensions = Lookup[blockSystems, "SystemSize", None];
  requests = Lookup[captures, "Requests", None];
  tasks = Lookup[captures, "TaskMetadata", None];
  domain = Lookup[contract, "Domain", None];
  If[!MemberQ[{"rational", "acb"}, domain] ||
      Lookup[contract, "SymbolNames", None] =!= {},
    err["E6", cs, <|"Contract" -> contract,
      "Detail" -> "native singular SCC column planning requires exact tasks in a Rational or Acb field without regulator symbols"|>]];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    zero = cppScalar[0, inputDigits, First[blockSystems]];
    one = cppScalar[1, inputDigits, First[blockSystems]]];
  plans = Table[Module[{blockRequests = requests[[seedBlock]],
      blockTasks = tasks[[seedBlock]], blockPlans},
    blockPlans = MapThread[Function[{rawRequest, task}, Module[
        {request = KeyTake[rawRequest, $nativeSCCColumnRunKeys],
         component, reachable = <||>, targetRecords = {},
         predecessors, sourceP, run, a, b},
        component = sccNativeCapturedJordanComponent[request, task,
          contract, dimensions[[seedBlock]], zero, one];
        If[component === None,
          err["E6", blockSystems[[seedBlock]], <|
            "SeedBlock" -> seedBlock,
            "Task" -> task,
            "RunShape" -> <|
              "nmax" -> Lookup[request, "nmax", None],
              "p" -> Lookup[request, "p", None],
              "initial_length" -> Length[Lookup[request, "initial", {}]],
              "validity_length" ->
                Length[Lookup[request, "initial_validity", {}]]|>,
            "Detail" -> "captured singular SCC request has no unique canonical log-zero unit component"|>]];
        a = task["a"]; b = task["b"];
        AssociateTo[reachable, seedBlock -> task["P"]];
        Do[If[target =!= seedBlock,
          predecessors = Select[
            First /@ Select[edges, Last[#] === target &],
            KeyExistsQ[reachable, #] &];
          If[predecessors =!= {},
            sourceP = Max[Lookup[reachable, predecessors]];
            run = sccNativeParticularRunForTag[
              blockSystems[[target]], contract, a, b, sourceP,
              inputDigits, domain];
            AssociateTo[reachable, target -> run["p"]];
            AppendTo[targetRecords,
              <|"Block" -> target, "Run" -> run|>]]],
          {target, topological}];
        <|"SeedLocalComponent" -> component,
          "SeedRun" -> KeyTake[request, $nativeSCCColumnRunKeys],
          "Tag" -> <|"a" -> a, "b" -> b,
            "p" -> task["P"]|>,
          "Targets" -> targetRecords|>]],
      {blockRequests, blockTasks}];
    blockPlans = SortBy[blockPlans, # ["SeedLocalComponent"] &];
    If[Lookup[blockPlans, "SeedLocalComponent", {}] =!=
        Range[dimensions[[seedBlock]]],
      err["E6", blockSystems[[seedBlock]], <|
        "DerivedLocalComponents" ->
          Lookup[blockPlans, "SeedLocalComponent", {}],
        "Detail" -> "captured singular SCC requests are not one complete canonical local basis"|>]];
    blockPlans], {seedBlock, Length[blockSystems]}];
  plans];

(* Select the exact Rational shadow before an Acb basis batch only from exact
   producer metadata.  This is a deterministic route-selection certificate,
   not a sampled numerical heuristic:

   - a CASE-P step in a captured seed schedule is checked before its Acb
     recurrence starts, and the complete basis submits every captured seed;
   - with identity gauge/spectral frames, a genuine polar cross-SCC
     multiplier shifts source support by -p.  The exact Rational executor is
     closed under those tag shifts and its later specialization proves the
     same physical equation, geometry, and public work window.  Different
     DAG paths could cancel a shifted sector, so this condition deliberately
     certifies that the Rational route is safe, not that Acb must fail.

   Nonidentity frames and target-only CASE-P schedules remain on the strict
   runtime gate: proving their composed tag support requires exact matrix
   composition, so this producer never guesses from individual entries. *)
sccNativeSeedRunContainsPseudoQ[columnPlans_List] := AnyTrue[
  Flatten[columnPlans, 1], Function[plan,
    AssociationQ[plan] && AnyTrue[
      Flatten[Lookup[Lookup[plan, "SeedRun", <||>], "schedule", {}], 1],
      AssociationQ[#] && Lookup[#, "case", None] === "P" &]]];

sccNativeActivePolarCouplingQ[coupling_Association] := AnyTrue[
  Lookup[coupling, "entries", {}], Function[entry,
    AssociationQ[entry] &&
      !TrueQ[Lookup[Lookup[entry, "multiplier", <||>],
        "proven_zero", True]] &&
      IntegerQ[Lookup[Lookup[entry, "multiplier", <||>],
        "center_pole_order", None]] &&
      Lookup[entry["multiplier"], "center_pole_order", 0] > 0]];

sccNativeRationalShadowDecision[domain_, blockRecords_List,
    couplings_List, columnPlans_List] := Module[
  {seedCaseP, identityFrames, polarCoupling, reason},
  If[domain =!= "acb", Return[<|
    "RequiresRationalShadow" -> False,
    "Certificate" -> "non-acb-domain"|>, Module]];
  seedCaseP = sccNativeSeedRunContainsPseudoQ[columnPlans];
  identityFrames = blockRecords =!= {} && AllTrue[blockRecords,
    AssociationQ[#] && TrueQ[Lookup[#, "identity_gauge", False]] &&
      TrueQ[Lookup[#, "identity_v", False]] &];
  polarCoupling = identityFrames &&
    AnyTrue[couplings, sccNativeActivePolarCouplingQ];
  reason = Which[
    seedCaseP, "exact-schedule-acb-case-p-compensation",
    polarCoupling,
      "exact-tagged-acb-polar-cross-coupling",
    True, "runtime-exact-schedule-and-tag-gate-required"];
  (* CASE-P and identity-frame polar couplings are more than local
     coefficient hazards.  Even when an Acb column solve succeeds, a later
     singular match must compose its right normal frame with the matching
     weights without losing their shared algebraic correlation.  The exact
     Rational shadow is the retained authority for that composition; a
     finite Acb fallback can manufacture a weight-radius residual which no
     Taylor-order or working-precision retry can reduce. *)
  <|"RequiresRationalShadow" -> TrueQ[seedCaseP || polarCoupling],
    "Certificate" -> reason,
    "RationalShadowFallback" -> TrueQ[seedCaseP || polarCoupling],
    "SeedCaseP" -> TrueQ[seedCaseP],
    "IdentityFrames" -> TrueQ[identityFrames],
    "PolarCrossCoupling" -> TrueQ[polarCoupling]|>];

sccNativeCachedRationalShadowDecision[cs_Association, req_Association,
    prepared_Association] := Module[
  {signature, cacheKey, cached, execution, decision},
  signature = sccNativeCompositeCacheSignature[cs, req];
  cacheKey = Hash[signature, "SHA256"];
  cached = Lookup[$nativeSCCCompositeCache, cacheKey, None];
  If[!AssociationQ[cached] ||
      !sccNativeCompositeExecutionDescriptorQ[cached, signature] ||
      !SameQ[Lookup[cached, "Result", None], prepared],
    err["E6", cs, <|"CacheKey" -> cacheKey,
      "Detail" -> "native SCC Rational-shadow decision is not collision-bound to the prepared public handle"|>]];
  execution = cached["Execution"];
  decision = Lookup[execution, "RationalShadowDecision", None];
  If[!AssociationQ[decision] ||
      !BooleanQ[Lookup[decision, "RequiresRationalShadow", None]] ||
      !StringQ[Lookup[decision, "Certificate", None]],
    err["E6", cs, <|"Decision" -> decision,
      "Detail" -> "cached native SCC Rational-shadow decision is malformed"|>]];
  decision];

sccNativeColumnRun[request_Association] :=
  KeyTake[request, $nativeSCCColumnRunKeys];

sccNativeParticularRunTemplate[request_Association] := Module[
  {run = sccNativeColumnRun[request], initial},
  initial = Lookup[run, "initial", {}];
  Join[run, <|"p" -> 0, "has_initial" -> False,
    "adaptive_probe" -> False, "a_target" -> "0",
    "b_target" -> "0", "a_shift_min" -> 0,
    "initial" -> ConstantArray["0", Length[initial]],
    "initial_validity" -> {Null}, "source" -> Null,
    "return_u" -> False|>]];

sccNativeRegularBlockParticularRunTemplate[request_Association,
    dimension_Integer] := Module[
  {run = sccNativeColumnRun[request], initial},
  initial = Lookup[run, "initial", {}];
  Join[run, <|"p" -> 0, "has_initial" -> False,
    "adaptive_probe" -> False, "a_target" -> "0",
    "b_target" -> "0", "a_shift_min" -> 0,
    "initial" -> ConstantArray["0", Length[initial]],
    "initial_validity" -> ConstantArray[Null, dimension],
    "source" -> Null, "return_u" -> False|>]];

sccNativeEncodedRegularBlockParticularRunTemplate[
    request_Association, dimension_Integer, zero_] := Module[
  {run = sccNativeColumnRun[request], initial},
  initial = Lookup[run, "initial", {}];
  Join[run, <|"p" -> 0, "has_initial" -> False,
    "adaptive_probe" -> False, "a_target" -> zero,
    "b_target" -> zero, "a_shift_min" -> 0,
    "initial" -> ConstantArray[zero, Length[initial]],
    "initial_validity" -> ConstantArray[Null, dimension],
    "source" -> Null, "return_u" -> False|>]];

sccNativeRegularBlockStatisticsQ[stats_, handle_Association,
    blockDimensions_List, domain_String:"rational"] := Module[
  {blockCharts = Lookup[stats, "block_charts", None],
   evidence = Lookup[stats, "capability_evidence", None], expected},
  expected = Switch[domain,
    "rational", "exact-rational-regular-block-dag-column-v2",
    "acb", "acb-regular-block-dag-column-v2",
    _, Return[False, Module]];
  sccNativeCompositeStatisticsQ[stats, handle] &&
    TrueQ[Lookup[stats, "execution_implemented", False]] &&
    Lookup[stats, "execution_scope", None] === expected &&
    TrueQ[Lookup[stats, "regular_block_dag_column_execution", False]] &&
    TrueQ[!Lookup[stats, "scalar_block_dag_column_execution", True]] &&
    TrueQ[!Lookup[stats, "general_scc_execution", True]] &&
    Lookup[stats, "execution_mode", None] === "BlockSequentialStrict" &&
    Lookup[stats, "dimension", None] === Total[blockDimensions] &&
    Lookup[stats, "blocks", None] === Length[blockDimensions] &&
    ListQ[blockCharts] && Length[blockCharts] === Length[blockDimensions] &&
    And @@ MapThread[Function[{record, block},
      AssociationQ[record] && Lookup[record, "block", None] === block - 1 &&
        Lookup[record, "dimension", None] === blockDimensions[[block]] &&
        StringQ[Lookup[record, "chart", None]] &&
        StringQ[Lookup[record, "principal_identity", None]]],
      {blockCharts, Range[Length[blockDimensions]]}] &&
    AssociationQ[evidence] &&
    Lookup[evidence, "identity_v", None] ===
      "native-retained-spectral-assembly-and-target-inverse" &&
    Lookup[evidence, "regular", None] ===
      "collision-bound-producer-certificate" &&
    (Lookup[evidence, "identity_gauge", None] ===
        "collision-bound-producer-certificate" ||
      (Lookup[evidence, "identity_gauge", None] ===
          "optional-fast-path-exact-directional-gauge-frame" &&
       Lookup[evidence, "exact_gauge", None] ===
          "retained-gauge-and-gauge-inverse-multiplier-certificates")) &&
    Lookup[evidence, "no_pseudo", None] ===
      "collision-bound-producer-certificate"];

sccNativeSingularBlockStatisticsQ[stats_, handle_Association,
    blockDimensions_List, scalarShape_, domain_String:"rational"] := Module[
  {blockCharts = Lookup[stats, "block_charts", None],
   evidence = Lookup[stats, "capability_evidence", None],
   expected, expectedNoPseudo, expectedPseudoExecution},
  expected = Switch[domain,
    "rational", If[TrueQ[scalarShape],
      "exact-rational-regular-singular-scalar-block-dag-column-v1",
      "exact-rational-regular-singular-jordan-block-dag-column-v2"],
    "acb", If[TrueQ[scalarShape],
      "acb-regular-singular-scalar-block-dag-column-v1",
      "acb-regular-singular-jordan-block-dag-column-v1"],
    _, Return[False, Module]];
  expectedNoPseudo = If[domain === "acb",
    "runtime-exact-schedule-case-p-gate",
    "producer-provenance-only-execution-revalidated-by-exact-schedule-certificate"];
  expectedPseudoExecution = If[domain === "acb",
    "exact-rational-certificate-case-p-rejected-for-acb",
    "exact-rational-joint-compensation-and-formal-overlap-certificate"];
  sccNativeCompositeStatisticsQ[stats, handle] &&
    TrueQ[Lookup[stats, "execution_implemented", False]] &&
    Lookup[stats, "execution_scope", None] === expected &&
    TrueQ[Lookup[stats,
      "regular_singular_jordan_block_dag_column_execution", False]] &&
    SameQ[TrueQ[Lookup[stats,
      "regular_singular_scalar_block_dag_column_execution", False]],
      TrueQ[scalarShape]] &&
    TrueQ[!Lookup[stats, "general_scc_execution", True]] &&
    Lookup[stats, "execution_mode", None] === "BlockSequentialStrict" &&
    Lookup[stats, "dimension", None] === Total[blockDimensions] &&
    Lookup[stats, "blocks", None] === Length[blockDimensions] &&
    ListQ[blockCharts] && Length[blockCharts] === Length[blockDimensions] &&
    And @@ MapThread[Function[{record, block},
      AssociationQ[record] && Lookup[record, "block", None] === block - 1 &&
        Lookup[record, "dimension", None] === blockDimensions[[block]] &&
        StringQ[Lookup[record, "chart", None]] &&
        StringQ[Lookup[record, "principal_identity", None]] &&
        AssociationQ[Lookup[record,
          "exact_affine_jordan_indicial", None]] &&
        Lookup[Lookup[record, "exact_affine_jordan_indicial", <||>],
          "dimension", None] === blockDimensions[[block]]],
      {blockCharts, Range[Length[blockDimensions]]}] &&
    AssociationQ[evidence] &&
    Lookup[evidence, "identity_v", None] ===
      "native-retained-spectral-assembly-and-target-inverse" &&
    Lookup[evidence, "regular_or_regular_singular", None] ===
      "collision-bound-producer-certificate" &&
    (Lookup[evidence, "identity_gauge", None] ===
        "collision-bound-producer-certificate" ||
      (Lookup[evidence, "identity_gauge", None] ===
          "optional-fast-path-exact-directional-gauge-frame" &&
       Lookup[evidence, "exact_gauge", None] ===
          "retained-gauge-and-gauge-inverse-multiplier-certificates")) &&
    Lookup[evidence, "jordan_indicial", None] ===
      "retained-exact-rational-full-matrix-certificate" &&
    Lookup[evidence, "no_pseudo", None] === expectedNoPseudo &&
    Lookup[evidence, "pseudo_schedule_execution", None] ===
      expectedPseudoExecution &&
    Lookup[evidence, "resonance_schedule", None] ===
      "retained-affine-jordan-verified-exact-captured-run"];

sccNativePseudoDiagnosticsQ[response_Association] := Module[
  {records = Lookup[response, "block_diagnostics", None]},
  ListQ[records] && records =!= {} && AllTrue[records, Function[record,
    Module[{hits, compensations, depth, uncompensated},
      If[!AssociationQ[record], Return[False, Module]];
      hits = Lookup[record, "pseudo_hit_count", None];
      compensations = Lookup[record, "pseudo_compensation_count", None];
      depth = Lookup[record, "max_pseudo_depth", None];
      uncompensated = Lookup[record,
        "uncompensated_pseudo_hit_count", None];
      And[IntegerQ[hits], hits >= 0,
        IntegerQ[compensations], compensations >= 0,
        IntegerQ[depth], depth >= 0,
        uncompensated === 0,
        If[hits === 0, True,
          (* A certified CASE-P hit need not materialize a compensation
             term: the exact Rational compensator deliberately omits a
             polar weight when every retained negative-epsilon gamma
             coefficient is exactly zero.  The positive depth proves that
             a genuine hit was inspected, while the exact value certificate
             and zero unresolved-hit count prove that omitting the term was
             not an approximate or fallback decision. *)
          depth > 0 &&
            TrueQ[Lookup[record, "pseudo_value_certified", False]]]]]]]];

sccNativeCanonicalJSONValue[value_Association] := Association@Map[
  Function[key, key -> sccNativeCanonicalJSONValue[value[key]]],
  Sort[Keys[value]]];
sccNativeCanonicalJSONValue[value_List] :=
  sccNativeCanonicalJSONValue /@ value;
sccNativeCanonicalJSONValue[value_] := value;

sccNativeBlockProvenanceIdentityQ[identity_, schema_String,
    parentIdentity_String, basisIndex_Integer, localComponent_,
    submittedSeed_Association, submittedTargets_List] := Module[
  {record, componentBindingQ},
  If[!StringQ[identity] || StringLength[identity] === 0,
    Return[False, Module]];
  record = Quiet[Check[ImportString[identity, "RawJSON"], $Failed]];
  componentBindingQ = AssociationQ[record] &&
    If[IntegerQ[localComponent],
      Lookup[record, "seed_local_component", None] ===
        localComponent - 1,
      !KeyExistsQ[record, "seed_local_component"]];
  AssociationQ[record] &&
    Lookup[record, "schema", None] === schema &&
    Lookup[record, "scc_exact_identity", None] === parentIdentity &&
    Lookup[record, "basis_index", None] === basisIndex &&
    componentBindingQ &&
    SameQ[sccNativeCanonicalJSONValue[Lookup[record, "seed", None]],
      sccNativeCanonicalJSONValue[submittedSeed]] &&
    SameQ[sccNativeCanonicalJSONValue[Lookup[record, "targets", None]],
      sccNativeCanonicalJSONValue[submittedTargets]]];

(* Public SCC-column summaries are bounded references to exact records retained
   by the native SCC/local owner and its checkpoint.  The FNV values are useful
   diagnostics only; exact authority comes from the already-validated live
   handles, checkpoint identity, and native column binding. *)
sccNativeColumnReferenceQ[reference_, scc_String, seedBlock_Integer,
    basisIndex_Integer] := Module[{diagnostics, sccFingerprint,
    columnFingerprint, sccBytes, columnBytes},
  If[!AssociationQ[reference], Return[False, Module]];
  diagnostics = Lookup[reference, "identity_diagnostics", None];
  If[!AssociationQ[diagnostics], Return[False, Module]];
  sccFingerprint = Lookup[diagnostics,
    "scc_exact_identity_fingerprint", None];
  columnFingerprint = Lookup[diagnostics,
    "exact_column_identity_fingerprint", None];
  sccBytes = Lookup[diagnostics, "scc_exact_identity_bytes", None];
  columnBytes = Lookup[diagnostics, "exact_column_identity_bytes", None];
  Sort[Keys[reference]] === Sort[{"schema", "authority", "scc",
      "seed_block", "basis_index", "identity_diagnostics"}] &&
    Sort[Keys[diagnostics]] === Sort[{"algorithm",
      "scc_exact_identity_fingerprint", "scc_exact_identity_bytes",
      "exact_column_identity_fingerprint",
      "exact_column_identity_bytes"}] &&
    Lookup[reference, "schema", None] ===
      "diffexp2-retained-scc-column-reference-v1" &&
    Lookup[reference, "authority", None] ===
      "retained-native-exact-column-owner" &&
    Lookup[reference, "scc", None] === scc &&
    Lookup[reference, "seed_block", None] === seedBlock &&
    Lookup[reference, "basis_index", None] === basisIndex &&
    Lookup[diagnostics, "algorithm", None] === "fnv1a64-v1" &&
    StringQ[sccFingerprint] &&
    StringMatchQ[sccFingerprint,
      RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
    StringQ[columnFingerprint] &&
    StringMatchQ[columnFingerprint,
      RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
    IntegerQ[sccBytes] && sccBytes > 0 &&
    IntegerQ[columnBytes] && columnBytes > 0];

sccNativeReachableTargetBlocks[components_List, edges_List,
    topological_List, seedBlock_Integer, cs_Association] := Module[
  {count = Length[components], positions, reachable, outgoing},
  If[Sort[topological] =!= Range[count] ||
      !AllTrue[edges, MatchQ[#, {_Integer, _Integer}] &&
        1 <= #[[1]] <= count && 1 <= #[[2]] <= count &],
    err["E6", cs, <|"TopologicalOrder" -> topological,
      "CondensationEdges" -> edges,
      "Detail" -> "cached native SCC execution graph is malformed"|>]];
  positions = AssociationThread[topological, Range[count]];
  If[!AllTrue[edges,
      positions[#[[1]]] < positions[#[[2]]] &],
    err["E6", cs, <|"TopologicalOrder" -> topological,
      "CondensationEdges" -> edges,
      "Detail" -> "cached native SCC condensation edges violate their exact topological order"|>]];
  reachable = ConstantArray[False, count];
  reachable[[seedBlock]] = True;
  outgoing[block_Integer] := Cases[edges, {block, target_} :> target];
  Do[If[TrueQ[reachable[[block]]],
    Scan[(reachable[[#]] = True) &, outgoing[block]]],
    {block, topological}];
  Select[topological, # =!= seedBlock && TrueQ[reachable[[#]]] &]];

sccNativeBuildColumnRequest[cs_Association, req_Association,
    seedBlock_Integer, seedLocalComponent_Integer:1] := Module[
  {prepared, signature, cacheKey, cached, execution, blockDimensions,
   runRecords, taskRecords, contract, domain, canonicalRequests, badBlocks,
   components, edges,
   topological, targetBlocks, seedRun, targetRuns, inputDigits,
   checkpointIdentity, checkpointDigest, seed, targets, stats,
   scalarExecution,
   componentMaps, seedRequestIndex, targetCanonicalRequests,
   selectedSeedLocalComponent, expectedBasisIndex, expectedCapability,
   expectedProvenanceSchema, columnPlans, selectedPlan,
   singularExecution, seedTag, targetTags, plannedTargetBlocks,
   provenanceLocalComponent, encodedZero, encodedOne, blockCharts,
   blockRegularity, capturedSingularExecution,
   advertisedSingularExecution},
  If[DownValues[DiffExp2`CppBackend`RunPersistentSCCColumn] === {},
    err["E5", cs, <|"Detail" ->
      "CppBackend persistent SCC column bridge is not available"|>]];
  prepared = PrepareNativeSCCComposite[cs, req];
  signature = sccNativeCompositeCacheSignature[cs, req];
  cacheKey = Hash[signature, "SHA256"];
  cached = Lookup[$nativeSCCCompositeCache, cacheKey, None];
  If[!AssociationQ[cached] ||
      !sccNativeCompositeExecutionDescriptorQ[cached, signature] ||
      !SameQ[Lookup[cached, "Result", None], prepared],
    err["E6", cs, <|"CacheKey" -> cacheKey,
      "Detail" -> "native SCC execution data is absent or not collision-bound to the prepared public handle"|>]];
  execution = cached["Execution"];
  blockDimensions = execution["BlockDimensions"];
  runRecords = execution["Runs"];
  taskRecords = execution["TaskMetadata"];
  contract = execution["Contract"];
  domain = Lookup[contract, "Domain", None];
  components = execution["Components"];
  edges = execution["CondensationEdges"];
  topological = execution["TopologicalOrder"];
  columnPlans = execution["ColumnPlans"];
  inputDigits = execution["InputDigits"];
  stats = Lookup[execution, "NativeStatistics", None];
  (* execution_scope is an executor admission result, not the mathematical
     chart classification.  In particular an exact regular-singular capture
     can temporarily advertise "unsupported" when one later native
     prerequisite (such as a polar cross-SCC coupling) is unavailable.  Do
     not then reinterpret honest p>0 Jordan seeds as regular unit columns.
     The retained block regularity is exact producer metadata, and the
     nonempty captured column plans independently certify the seed shapes.
     The strict capability predicate below remains authoritative before any
     request is submitted to C++. *)
  blockCharts = If[AssociationQ[stats],
    Lookup[stats, "block_charts", None], None];
  blockRegularity = If[ListQ[blockCharts],
    Lookup[blockCharts, "regular", None], None];
  capturedSingularExecution = ListQ[blockRegularity] &&
    Length[blockRegularity] === Length[blockDimensions] &&
    AllTrue[blockRegularity, MemberQ[{True, False}, #] &] &&
    MemberQ[blockRegularity, False] && ListQ[columnPlans] &&
    Length[columnPlans] === Length[blockDimensions] &&
    AllTrue[columnPlans, ListQ[#] && # =!= {} &];
  advertisedSingularExecution = AssociationQ[stats] &&
    MemberQ[{
      "exact-rational-regular-singular-scalar-block-dag-column-v1",
      "exact-rational-regular-singular-jordan-block-dag-column-v2",
      "acb-regular-singular-scalar-block-dag-column-v1",
      "acb-regular-singular-jordan-block-dag-column-v1"},
      Lookup[stats, "execution_scope", None]];
  singularExecution = capturedSingularExecution ||
    advertisedSingularExecution;
  scalarExecution = AllTrue[blockDimensions, # === 1 &];
  If[Length[blockDimensions] < 1 ||
      Length[runRecords] =!= Length[blockDimensions] ||
      Length[taskRecords] =!= Length[blockDimensions] ||
      Length[columnPlans] =!= Length[blockDimensions] ||
      Length[components] =!= Length[blockDimensions] ||
      !AllTrue[blockDimensions, IntegerQ[#] && # > 0 &] ||
      !(And @@ MapThread[
        ListQ[#1] && Length[#1] === #2 &, {components, blockDimensions}]) ||
      Total[blockDimensions] =!= cs["SystemSize"] ||
      !MemberQ[{"rational", "acb"}, domain] ||
      Lookup[contract, "SymbolNames", None] =!= {},
    err["E6", cs, <|"Contract" -> contract,
      "BlockDimensions" -> blockDimensions,
      "Detail" -> "native SCC column requires one or more exact-rational or Acb blocks with dimensions matching the parent SCC partition and no regulator field"|>]];
  If[advertisedSingularExecution && !MemberQ[
      Switch[domain,
        "rational", {
          "exact-rational-regular-singular-scalar-block-dag-column-v1",
          "exact-rational-regular-singular-jordan-block-dag-column-v2"},
        "acb", {
          "acb-regular-singular-scalar-block-dag-column-v1",
          "acb-regular-singular-jordan-block-dag-column-v1"},
        _, {}], Lookup[stats, "execution_scope", None]],
    err["E6", cs, <|"Contract" -> contract,
      "NativeStatistics" -> stats,
      "Detail" -> "native singular SCC capability does not match the captured scalar domain"|>]];
  If[!IntegerQ[inputDigits] || inputDigits < 1,
    err["E6", cs, <|"InputDigits" -> inputDigits,
      "Detail" -> "cached native SCC local-metadata precision is malformed"|>]];
  If[seedBlock < 1 || seedBlock > Length[blockDimensions],
    err["E6", cs, <|"SeedBlock" -> seedBlock,
      "Blocks" -> Length[blockDimensions],
      "Detail" -> "native SCC seed block is outside the one-based block range"|>]];
  If[seedLocalComponent < 1 ||
      seedLocalComponent > blockDimensions[[seedBlock]],
    err["E6", cs, <|"SeedBlock" -> seedBlock,
      "SeedLocalComponent" -> seedLocalComponent,
      "BlockDimension" -> blockDimensions[[seedBlock]],
      "Detail" -> "native SCC seed component is outside the one-based local block range"|>]];
  If[!(And @@ MapThread[Function[{runs, tasks, dimension},
      ListQ[runs] && Length[runs] === dimension &&
        AllTrue[runs, AssociationQ] && ListQ[tasks] &&
        Length[tasks] === dimension && AllTrue[tasks, AssociationQ]],
      {runRecords, taskRecords, blockDimensions}]),
    err["E6", cs, <|"Detail" ->
      "native SCC column requires one captured homogeneous request per local block component"|>]];
  targetBlocks = sccNativeReachableTargetBlocks[components, edges,
    topological, seedBlock, cs];
  If[singularExecution,
    If[!AllTrue[columnPlans, ListQ] ||
        Length[columnPlans[[seedBlock]]] =!=
          blockDimensions[[seedBlock]],
      err["E6", cs, <|"SeedBlock" -> seedBlock,
        "Detail" -> "cached singular SCC column plans do not cover the selected block dimension"|>]];
    selectedPlan = Select[columnPlans[[seedBlock]],
      Lookup[#, "SeedLocalComponent", None] === seedLocalComponent &];
    If[Length[selectedPlan] =!= 1,
      err["E6", cs, <|"SeedBlock" -> seedBlock,
        "SeedLocalComponent" -> seedLocalComponent,
        "Detail" -> "cached singular SCC column plan does not select exactly one canonical seed"|>]];
    selectedPlan = First[selectedPlan];
    If[!AssociationQ[Lookup[selectedPlan, "Tag", None]] ||
        !ListQ[Lookup[selectedPlan, "Targets", None]] ||
        !AllTrue[selectedPlan["Targets"], AssociationQ[#] &&
          IntegerQ[Lookup[#, "Block", None]] &&
          AssociationQ[Lookup[#, "Run", None]] &],
      err["E6", cs, <|"SeedBlock" -> seedBlock,
        "Detail" -> "cached singular SCC column plan is malformed"|>]];
    seedRun = selectedPlan["SeedRun"];
    plannedTargetBlocks = Lookup[selectedPlan["Targets"], "Block", {}];
    If[plannedTargetBlocks =!= targetBlocks,
      err["E6", cs, <|"SeedBlock" -> seedBlock,
        "ExpectedTargets" -> targetBlocks,
        "PlannedTargets" -> plannedTargetBlocks,
        "Detail" -> "cached singular SCC column plan differs from the exact condensation reachability order"|>]];
    targetRuns = Lookup[selectedPlan["Targets"], "Run", {}];
    selectedSeedLocalComponent = seedLocalComponent;
    expectedBasisIndex =
      components[[seedBlock, selectedSeedLocalComponent]] - 1;
    seedTag = selectedPlan["Tag"];
    targetTags = Map[<|"a" -> seedTag["a"],
        "b" -> seedTag["b"], "p" -> # ["p"]|> &, targetRuns];
    expectedCapability = Switch[domain,
      "rational", If[scalarExecution,
        "exact-rational-regular-singular-scalar-block-dag-column-v1",
        "exact-rational-regular-singular-jordan-block-dag-column-v2"],
      "acb", If[scalarExecution,
        "acb-regular-singular-scalar-block-dag-column-v1",
        "acb-regular-singular-jordan-block-dag-column-v1"]];
    expectedProvenanceSchema = Switch[domain,
      "rational", If[scalarExecution,
        "diffexp2-native-scc-regular-singular-scalar-column-v1",
        "diffexp2-native-scc-regular-singular-jordan-column-v2"],
      "acb", If[scalarExecution,
        "diffexp2-native-scc-acb-regular-singular-scalar-column-v1",
        "diffexp2-native-scc-acb-regular-singular-jordan-column-v1"]],

    If[domain === "acb",
      Block[{$cppSerializationDomain = "acb",
          $cppSerializationSymbols = {}},
        encodedZero = cppScalar[0, inputDigits, cs];
        encodedOne = cppScalar[1, inputDigits, cs];
      componentMaps = MapThread[Function[{runs, tasks, dimension},
          MapThread[Function[{run, task},
            sccNativeCanonicalEncodedRegularBlockComponent[
              run, task, contract, dimension, encodedZero, encodedOne,
              cppScalar[task["a"], inputDigits, cs],
              cppScalar[task["b"], inputDigits, cs],
              Table[cppScalar[task["a"] + n, inputDigits, cs],
                {n, 0, contract["NMax"]}]]], {runs, tasks}]],
        {runRecords, taskRecords, blockDimensions}]];
      badBlocks = Select[Range[Length[componentMaps]], Function[block,
        MemberQ[componentMaps[[block]], None] ||
          Sort[componentMaps[[block]]] =!=
            Range[blockDimensions[[block]]]]];
      If[badBlocks =!= {},
        err["E6", cs, <|"UnsupportedBlocks" -> badBlocks,
          "DerivedLocalComponents" -> componentMaps,
          "Detail" -> "captured regular Acb block requests are not one complete exact-encoded eps^0 local basis"|>]];
      seedRequestIndex = First@FirstPosition[
        componentMaps[[seedBlock]], seedLocalComponent];
      selectedSeedLocalComponent =
        componentMaps[[seedBlock, seedRequestIndex]];
      seedRun = sccNativeColumnRun[
        runRecords[[seedBlock, seedRequestIndex]]];
      targetCanonicalRequests = MapThread[
        Function[{runs, localComponents},
          runs[[First@FirstPosition[localComponents, 1]]]],
        {runRecords, componentMaps}];
      targetRuns = MapThread[
        sccNativeEncodedRegularBlockParticularRunTemplate[
          #1, #2, encodedZero] &,
        {targetCanonicalRequests[[targetBlocks]],
         blockDimensions[[targetBlocks]]}];
      expectedBasisIndex =
        components[[seedBlock, selectedSeedLocalComponent]] - 1;
      expectedCapability = If[scalarExecution,
        "acb-regular-scalar-block-dag-column-v1",
        "acb-regular-block-dag-column-v2"];
      expectedProvenanceSchema = If[scalarExecution,
        "diffexp2-native-scc-acb-regular-scalar-column-v1",
        "diffexp2-native-scc-acb-regular-column-v2"],

      (* Keep exact regular scalar v1 on its original request/template and
         checkpoint path. In particular, the optional component is not added
         to its wire identity or returned record. *)
      If[scalarExecution,
      canonicalRequests = First /@ runRecords;
      badBlocks = Select[Range[Length[canonicalRequests]],
        !sccNativeCanonicalScalarHomogeneousRequestQ[
          canonicalRequests[[#]], contract] &];
      If[badBlocks =!= {},
        err["E6", cs, <|"UnsupportedBlocks" -> badBlocks,
          "Detail" -> "captured scalar block request is not the canonical exact eps^0 regular homogeneous run"|>]];
      seedRun = sccNativeColumnRun[canonicalRequests[[seedBlock]]];
      targetRuns = sccNativeParticularRunTemplate /@
        canonicalRequests[[targetBlocks]];
      selectedSeedLocalComponent = 1;
      expectedBasisIndex = First[components[[seedBlock]]] - 1;
      expectedCapability =
        "exact-rational-regular-scalar-block-dag-column-v1";
      expectedProvenanceSchema =
        "diffexp2-native-scc-column-v1",
      componentMaps = MapThread[Function[{runs, dimension},
          sccNativeCanonicalRegularBlockComponent[#, contract, dimension] & /@
            runs], {runRecords, blockDimensions}];
      badBlocks = Select[Range[Length[componentMaps]], Function[block,
        MemberQ[componentMaps[[block]], None] ||
          Sort[componentMaps[[block]]] =!= Range[blockDimensions[[block]]]]];
      If[badBlocks =!= {},
        err["E6", cs, <|"UnsupportedBlocks" -> badBlocks,
          "DerivedLocalComponents" -> componentMaps,
          "Detail" -> "captured regular block requests are not a complete permutation of canonical exact eps^0 local unit columns"|>]];
      seedRequestIndex = First@FirstPosition[
        componentMaps[[seedBlock]], seedLocalComponent];
      selectedSeedLocalComponent =
        componentMaps[[seedBlock, seedRequestIndex]];
      seedRun = sccNativeColumnRun[
        runRecords[[seedBlock, seedRequestIndex]]];
      targetCanonicalRequests = MapThread[
        Function[{runs, localComponents},
          runs[[First@FirstPosition[localComponents, 1]]]],
        {runRecords, componentMaps}];
      targetRuns = MapThread[
        sccNativeRegularBlockParticularRunTemplate,
        {targetCanonicalRequests[[targetBlocks]],
         blockDimensions[[targetBlocks]]}];
      expectedBasisIndex =
        components[[seedBlock, selectedSeedLocalComponent]] - 1;
      expectedCapability =
        "exact-rational-regular-block-dag-column-v2";
      expectedProvenanceSchema =
        "diffexp2-native-scc-column-v2"]];
    seedTag = <|"a" -> 0, "b" -> 0, "p" -> 0|>;
    targetTags = Map[<|"a" -> 0, "b" -> 0,
        "p" -> # ["p"]|> &, targetRuns]];
  provenanceLocalComponent = If[scalarExecution, None,
    selectedSeedLocalComponent];
  checkpointDigest = Hash[{cacheKey,
      execution["ParentIdentity"], prepared, seedBlock, seedRun,
      MapThread[Rule, {targetBlocks, targetRuns}]}, "SHA256"];
  checkpointIdentity = "de2-native-scc-column-" <>
    IntegerString[checkpointDigest, 16, 64];
  seed = <|"block" -> seedBlock - 1, "run" -> seedRun,
    "metadata" -> cppNativeLocalMetadata[cs,
      seedTag["a"], seedTag["b"], seedTag["p"], inputDigits,
      checkpointIdentity <> ":seed:" <> ToString[seedBlock - 1]]|>;
  targets = MapThread[Function[{block, run, tag}, <|
      "block" -> block - 1, "run" -> run,
      "metadata" -> cppNativeLocalMetadata[cs,
        tag["a"], tag["b"], tag["p"], inputDigits,
        checkpointIdentity <> ":target:" <> ToString[block - 1]]|>],
    {targetBlocks, targetRuns, targetTags}];
  If[FailureQ[stats] ||
      !sccNativeCompositeStatisticsQ[stats, prepared] ||
      !TrueQ[Lookup[stats, "execution_implemented", False]] ||
      Lookup[stats, "execution_scope", None] =!=
        expectedCapability ||
      If[singularExecution,
        !sccNativeSingularBlockStatisticsQ[
          stats, prepared, blockDimensions, scalarExecution, domain],
        !scalarExecution &&
          !sccNativeRegularBlockStatisticsQ[
            stats, prepared, blockDimensions, domain]],
    err["E6", cs, <|"NativeStatistics" -> stats,
      "ExpectedCapability" -> expectedCapability,
      "Detail" -> "retained native SCC does not advertise the required strict regular or exact affine-Jordan block-DAG column capability"|>]];
  <|"Prepared" -> prepared, "Seed" -> seed, "Targets" -> targets,
    "CheckpointIdentity" -> checkpointIdentity,
    "ExpectedCapability" -> expectedCapability,
    "ExpectedProvenanceSchema" -> expectedProvenanceSchema,
    "ExpectedBasisIndex" -> expectedBasisIndex,
    "ParentIdentity" -> execution["ParentIdentity"],
    "ProvenanceLocalComponent" -> provenanceLocalComponent,
    "ScalarExecution" -> scalarExecution,
    "SingularExecution" -> singularExecution,
    "SelectedSeedLocalComponent" -> selectedSeedLocalComponent,
    "SeedBlock" -> seedBlock|>];

sccNativeFinalizeColumn[cs_Association, req_Association,
    spec_Association, response_] := Module[
  {prepared = spec["Prepared"], seed = spec["Seed"],
   targets = spec["Targets"],
   checkpointIdentity = spec["CheckpointIdentity"],
   expectedCapability = spec["ExpectedCapability"],
   expectedProvenanceSchema = spec["ExpectedProvenanceSchema"],
   expectedBasisIndex = spec["ExpectedBasisIndex"],
   parentIdentity = spec["ParentIdentity"],
   provenanceLocalComponent = spec["ProvenanceLocalComponent"],
   scalarExecution = spec["ScalarExecution"],
   singularExecution = spec["SingularExecution"],
   selectedSeedLocalComponent = spec["SelectedSeedLocalComponent"],
   seedBlock = spec["SeedBlock"], backendID, provenance, requiredMax,
   forbiddenPayloadKeys, result},
  requiredMax = Lookup[req, "RequiredCompleteMax",
    req["EpsWindow", "CompleteMax"]];
  If[FailureQ[response],
    err["E5", cs, <|"BackendFailure" -> response,
      "Detail" -> "persistent native SCC basis-column solve failed"|>]];
  If[!AssociationQ[response] ||
      Lookup[response, "status", "error"] =!= "ok",
    backendID = Lookup[response, "id", "E5"];
    err[If[MemberQ[{"E4", "E5", "E6"}, backendID], backendID, "E5"],
      cs, <|"BackendID" -> backendID,
        "Detail" -> Lookup[response, "detail",
          "persistent native SCC basis-column solve returned an error"]|>]];
  provenance = Lookup[response, "column_provenance", None];
  forbiddenPayloadKeys = Intersection[Keys[response],
    {"assembled", "coefficients", "u", "validity"}];
  If[!StringQ[Lookup[response, "session", None]] ||
      !StringQ[Lookup[response, "local", None]] ||
      Lookup[response, "session", None] =!= prepared["Session"] ||
      Lookup[response, "scc", None] =!= prepared["SCC"] ||
      Lookup[response, "chart", None] =!= prepared["SCC"] ||
      Lookup[response, "dimension", None] =!= cs["SystemSize"] ||
      !IntegerQ[Lookup[response, "epsilon_min", None]] ||
      Lookup[response, "epsilon_min", 1] >
        requiredMax ||
      !IntegerQ[Lookup[response, "epsilon_max", None]] ||
      !TrueQ[requiredMax <= response["epsilon_max"]] ||
      Lookup[response, "taylor_complete_max", None] =!= req["TOrder"] ||
      Lookup[response, "checkpoint_identity", None] =!= checkpointIdentity ||
      !TrueQ[Lookup[response, "native_retained", False]] ||
      Lookup[response, "json_coefficients", None] =!= 0 ||
      Lookup[response, "execution_capability", None] =!=
        expectedCapability ||
      Lookup[response, "pseudo_hit_count", None] =!= 0 ||
      (singularExecution && !sccNativePseudoDiagnosticsQ[response]) ||
      forbiddenPayloadKeys =!= {} || !AssociationQ[provenance] ||
      !sccNativeColumnReferenceQ[provenance, prepared["SCC"],
        seedBlock - 1, expectedBasisIndex],
    Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[response]];
    err["E6", cs, <|"BackendResponse" -> response,
      "ForbiddenPayloadKeys" -> forbiddenPayloadKeys,
      "Detail" -> "native SCC basis-column summary violated its opaque retained-local or exact-provenance contract"|>]];
  result = <|"Type" -> "DiffExp2NativeSCCBasisColumn",
    "Session" -> response["session"], "Local" -> response["local"],
    "NativeSCC" -> response["scc"], "NativeChart" -> response["chart"],
    "SeedBlock" -> seedBlock,
    "BasisIndex" -> Lookup[provenance, "basis_index"] + 1,
    "Chart" -> <|"Center" -> cs["Center"],
      "ChartMap" -> cs["ChartMap"], "Radius" -> cs["Radius"],
      "Prescriptions" -> cs["Prescriptions"]|>,
    "EpsWindow" -> <|"Min" -> response["epsilon_min"],
      "CompleteMax" -> response["epsilon_max"]|>,
    "TWindow" -> <|"CompleteMax" ->
      response["taylor_complete_max"]|>,
    "CheckpointIdentity" -> checkpointIdentity,
    "ColumnProvenance" -> provenance,
    "NativeSummary" -> KeyDrop[response,
      {"status", "session", "local", "scc", "chart", "metadata",
       "column_provenance"}]|>;
  If[scalarExecution, result,
    Join[result,
      <|"SeedLocalComponent" -> selectedSeedLocalComponent|>]]];

sccNativeDiagnosticFrameOverrideQ[] := AnyTrue[{
    "DE2_SCC_DIAGNOSTIC_START_FRAME",
    "DE2_SCC_DIAGNOSTIC_STRICT_PROBE",
    "DE2_SCC_DIAGNOSTIC_LOWER_EXTRA",
    "DE2_SCC_DIAGNOSTIC_TOP_HALO"},
  StringQ[Environment[#]] &];

sccNativeFrameFailureKind[failure_] := Module[
  {text = ToLowerCase[ToString[failure, InputForm]]},
  Which[
    AnyTrue[{
        "below the work frame",
        "below the target frame",
        "lower-frame",
        "retained lower epsilon frame",
        "insufficient lower halo",
        "lowest framed",
        "cancellation guard",
        "epsilon shift would discard nonzero",
        "epsilon shift lies wholly below",
        "pseudo jordan solve would discard",
        "canonical jordan seed exceeds"},
      StringContainsQ[text, #] &], "Lower",
    AnyTrue[{
        "above the work frame",
        "above the retained epsilon reservoir",
        "above retained epsilon reservoir",
        "needs coefficients above the retained epsilon reservoir",
        "insufficient upper halo",
        "retained upper epsilon",
        "polar frame incomplete"},
      StringContainsQ[text, #] &], "Upper",
    True, None]];

sccNativeRethrowFrameFailure[cs_Association, failure_,
    attempts_List] := Module[{payload},
  payload = If[FailureQ[failure] &&
      AssociationQ[Quiet[Check[failure[[2]], None]]],
    failure[[2]], <|"ID" -> "E5",
      "Detail" -> ToString[failure, InputForm]|>];
  err[Lookup[payload, "ID", "E5"], cs,
    Join[KeyDrop[payload, {"ID", "Module", "Chart"}],
      If[Length[attempts] > 1,
        <|"SCCFrameAttempts" -> attempts,
          "Detail" -> Lookup[payload, "Detail",
            "bounded native SCC epsilon-frame retries were exhausted"]|>,
        <||>]]]];

SetAttributes[sccNativeWithFrameRetries, HoldRest];
sccNativeWithFrameRetries[cs_Association, req_Association,
    expression_] := Module[
  {automaticQ = !sccNativeDiagnosticFrameOverrideQ[],
   lowerExtra = sccNativeEffectiveFrameLowerExtra[cs],
   topHalo = sccNativeEffectiveFrameTopHalo[],
   attempt = 0, result = None, kind, attempts = {}, next},
  If[!automaticQ, Return[expression, Module]];
  WithNativeSCCCompositeCacheReservation[$nativeSCCFrameRetryLimit,
    While[attempt < $nativeSCCFrameRetryLimit,
      attempt++;
      (* DE2Error prints before throwing.  Intermediate bounded probes are
         intentionally silent; a terminal failure is reissued once below
         with the complete compact-attempt history. *)
      result = Block[{Print = (Null &),
          $nativeSCCFrameLowerExtra = lowerExtra,
          $nativeSCCFrameTopHalo = topHalo},
        Catch[expression, "DiffExp2Error"]];
      If[!FailureQ[result],
        If[Environment["DE2_SCC_FRAME_RETRY_TIMING"] === "1" &&
            attempts =!= {},
          Print["DE2 SCC FRAME RETRY recovered attempts=",
            Append[attempts, <|"Attempt" -> attempt,
              "LowerExtra" -> lowerExtra, "TopHalo" -> topHalo,
              "Result" -> "ok"|>]]];
        Return[result, Module]];
      kind = sccNativeFrameFailureKind[result];
      AppendTo[attempts, <|"Attempt" -> attempt,
        "LowerExtra" -> lowerExtra, "TopHalo" -> topHalo,
        "FailureEdge" -> kind|>];
      next = Switch[kind,
        "Lower", If[lowerExtra < $nativeSCCFrameLowerExtraMax,
          Min[$nativeSCCFrameLowerExtraMax, 2 lowerExtra], None],
        "Upper", If[topHalo < $nativeSCCFrameTopHaloMax,
          Min[$nativeSCCFrameTopHaloMax, 2 topHalo], None],
        _, None];
      If[!IntegerQ[next],
        Return[sccNativeRethrowFrameFailure[
          cs, result, attempts], Module]];
      If[kind === "Lower", lowerExtra = next, topHalo = next]]];
  sccNativeRethrowFrameFailure[cs, result, attempts]];

solveNativeSCCBasisColumnAttempt[cs_Association, req_Association,
    seedBlock_Integer, seedLocalComponent_Integer:1] := Module[
  {spec, response},
  spec = sccNativeBuildColumnRequest[
    cs, req, seedBlock, seedLocalComponent];
  response = DiffExp2`CppBackend`RunPersistentSCCColumn[
    spec["Prepared"], spec["Seed"], spec["Targets"],
    spec["CheckpointIdentity"]];
  sccNativeFinalizeColumn[cs, req, spec, response]];

SolveNativeSCCBasisColumn[cs_Association, req_Association,
    seedBlock_Integer, seedLocalComponent_Integer:1] :=
  sccNativeWithFrameRetries[cs, req,
    solveNativeSCCBasisColumnAttempt[
      cs, req, seedBlock, seedLocalComponent]];

solveNativeSCCBasisAttempt[cs_Association, req_Association,
    threads_:Automatic] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], components,
   seeds, specs,
   prepared, workerCount, requestColumns, batch, raw, columns, cleanup,
   forbiddenPayloadKeys},
  components = If[AssociationQ[seq],
    Lookup[seq, "Components", {}], {}];
  If[!ListQ[components] || Length[components] < 1 ||
      (Length[components] === 1 &&
        TrueQ[Lookup[Lookup[cs, "IndicialData", <||>],
          "Regular", False]]),
    err["E6", cs, <|"Detail" ->
      "native SCC basis batch requires a multi-block exact SCC certificate or a singular one-block certificate"|>]];
  If[DownValues[DiffExp2`CppBackend`RunPersistentSCCColumns] === {},
    err["E5", cs, <|"Detail" ->
      "CppBackend persistent SCC column-batch bridge is not available"|>]];
  seeds = Flatten[Table[{block, component},
      {block, Length[components]},
      {component, Length[components[[block]]]}], 1];
  specs = sccNativeBuildColumnRequest[cs, req, #[[1]], #[[2]]] & /@
    seeds;
  If[Length[specs] =!= cs["SystemSize"] ||
      !sccAllSameQ[Lookup[specs, "Prepared", None]],
    err["E6", cs, <|"Detail" ->
      "native SCC basis columns did not bind one complete shared composite"|>]];
  prepared = First[Lookup[specs, "Prepared"]];
  workerCount = Which[
    threads === Automatic, cppConfiguredThreads[Length[specs]],
    IntegerQ[threads] && threads > 0, Min[threads, Length[specs]],
    True, err["E6", cs, <|"Threads" -> threads,
      "Detail" -> "native SCC basis batch thread count must be a positive integer or Automatic"|>]];
  requestColumns = Map[<|"Seed" -> # ["Seed"],
      "Targets" -> # ["Targets"],
      "CheckpointIdentity" -> # ["CheckpointIdentity"]|> &, specs];
  batch = DiffExp2`CppBackend`RunPersistentSCCColumns[
    prepared, requestColumns, workerCount];
  If[FailureQ[batch] || !AssociationQ[batch] ||
      Lookup[batch, "status", "error"] =!= "ok",
    err["E5", cs, <|"BackendFailure" -> batch,
      "Detail" -> "persistent native SCC basis batch failed"|>]];
  raw = Lookup[batch, "results", None];
  forbiddenPayloadKeys = Intersection[Keys[batch],
    {"assembled", "coefficients", "u", "validity"}];
  If[!ListQ[raw] || Length[raw] =!= Length[specs] ||
      Lookup[batch, "session", None] =!= prepared["Session"] ||
      Lookup[batch, "scc", None] =!= prepared["SCC"] ||
      Lookup[batch, "columns", None] =!= Length[specs] ||
      !TrueQ[Lookup[batch, "atomic_retention", False]] ||
      Lookup[batch, "json_coefficients", None] =!= 0 ||
      forbiddenPayloadKeys =!= {},
    If[ListQ[raw], Scan[
      If[AssociationQ[#] && StringQ[Lookup[#, "local", None]],
        Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]]] &, raw]];
    err["E6", cs, <|"BackendResponse" -> batch,
      "ForbiddenPayloadKeys" -> forbiddenPayloadKeys,
      "Detail" -> "native SCC basis batch violated its ordered opaque-retention contract"|>]];
  cleanup[] := Scan[
    If[AssociationQ[#] && StringQ[Lookup[#, "local", None]],
      Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]]] &, raw];
  columns = Catch[
    MapThread[sccNativeFinalizeColumn[cs, req, #1, #2] &,
      {specs, raw}],
    "DiffExp2Error", Function[{failure, tag},
      cleanup[]; Throw[failure, tag]]];
  columns = SortBy[columns, # ["BasisIndex"] &];
  If[Lookup[columns, "BasisIndex", {}] =!= Range[cs["SystemSize"]],
    cleanup[];
    err["E6", cs, <|"BasisIndices" ->
      Lookup[columns, "BasisIndex", {}],
      "Detail" -> "native SCC basis batch did not cover every physical basis index exactly once"|>]];
  <|"Type" -> "DiffExp2NativeSCCBasis",
    "Session" -> prepared["Session"], "NativeSCC" -> prepared["SCC"],
    "Columns" -> columns, "Dimension" -> cs["SystemSize"],
    "Chart" -> <|"Center" -> cs["Center"],
      "ChartMap" -> cs["ChartMap"], "Radius" -> cs["Radius"],
      "Prescriptions" -> cs["Prescriptions"]|>,
    "EpsWindow" -> <|"Min" -> Min[
        Lookup[Lookup[columns, "EpsWindow"], "Min"]],
      "CompleteMax" -> Min[
        Lookup[Lookup[columns, "EpsWindow"], "CompleteMax"]]|>,
    "TWindow" -> <|"CompleteMax" -> req["TOrder"]|>,
    "NativeSummary" -> KeyDrop[batch, {"results"}]|>];

SolveNativeSCCBasis[cs_Association, req_Association,
    threads_:Automatic] :=
  sccNativeWithFrameRetries[cs, req,
    solveNativeSCCBasisAttempt[cs, req, threads]];

SolveNativeRegularBasis[cs_Association, req_Association,
    threads_:Automatic, forceMonolithic_:False,
    equationOwner_:Automatic] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], d = cs["SystemSize"],
   epsWindow = Lookup[req, "EpsWindow", None], epsMin, epsMax,
   physical, tag = <|"a" -> 0, "b" -> 0, "p" -> 0|>, unitInitials,
   shared, prepared, workerCount, batch, raw, columns, cleanup, sessions,
   charts, forbiddenPayloadKeys, expectedLocalChart,
   expectedLocalOperator, ownerPhysicalPreparation},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "SolveNativeRegularBasis requires a regular chart"|>]];
  If[threads =!= Automatic && !(IntegerQ[threads] && threads > 0),
    err["E6", cs, <|"Threads" -> threads,
      "Detail" -> "native regular basis thread count must be a positive integer or Automatic"|>]];
  If[equationOwner === Automatic && !TrueQ[forceMonolithic] &&
      AssociationQ[seq] &&
      Length[Lookup[seq, "Components", {}]] > 1,
    Return[SolveNativeSCCBasis[cs, req, threads], Module]];
  If[!AssociationQ[epsWindow] ||
      !IntegerQ[Lookup[epsWindow, "Min", None]] ||
      !IntegerQ[Lookup[epsWindow, "CompleteMax", None]],
    err["E8", cs, <|"Request" -> req,
      "Detail" -> "native regular basis requires a finite epsilon window"|>]];
  epsMin = epsWindow["Min"];
  epsMax = epsWindow["CompleteMax"];
  If[epsMin > 0 || epsMax < 0,
    err["E6", cs, <|"EpsWindow" -> epsWindow,
      "Detail" -> "native regular basis normalization requires eps^0 inside the requested window"|>]];
  physical = regularPhysicalChartSystem[cs];
  unitInitials = Table[{Table[DiffExp2`EpsSeries`ESNew[epsMin,
      Table[If[component === basisIndex && power === 0, 1, 0],
        {power, epsMin, epsMax}]], {component, d}]},
    {basisIndex, d}];
  ownerPhysicalPreparation = If[
    AssociationQ[equationOwner] &&
      Lookup[equationOwner, "OwnerKind", None] ===
        "FrameIndependentRegularPhysicalEquation" &&
      AssociationQ[Lookup[equationOwner, "PhysicalPreparation", None]],
    equationOwner["PhysicalPreparation"], Automatic];
  If[AssociationQ[ownerPhysicalPreparation] &&
      (Lookup[ownerPhysicalPreparation, "EquationIdentity", None] =!=
          Lookup[equationOwner, "EquationIdentity", None] ||
       Lookup[ownerPhysicalPreparation, "PhysicalPayloadIdentity", None] =!=
          Lookup[equationOwner, "PhysicalPayloadIdentity", None]),
    err["E6", cs, <|
      "Detail" ->
        "regular basis owner preparation is not bound to the supplied physical equation owner"|>]];
  (* The immutable chart preparation is the dominant Wolfram-side cost for
     large regular receivers.  Capture it once, then serialize only the d
     distinct unit-column runs and checkpoint identities. *)
  shared = prepareNativeLocalFamilyShared[
    physical, req, tag, First[unitInitials], ownerPhysicalPreparation];
  prepared = prepareNativeLocalFamilyDynamic[shared, #, False] & /@
    unitInitials;
  If[Length[prepared] =!= d ||
      !Apply[SameQ, Lookup[prepared, "PersistentMetadata"]],
    err["E6", cs, <|"Detail" ->
      "regular basis dynamic runs did not retain one exact immutable preparation"|>]];
  workerCount = Which[
    threads === Automatic, cppConfiguredThreads[d],
    IntegerQ[threads] && threads > 0, Min[threads, d],
    True, err["E6", cs, <|"Threads" -> threads,
      "Detail" ->
        "native regular basis batch thread count must be a positive integer or Automatic"|>]];
  If[DownValues[DiffExp2`CppBackend`RunPersistentLocalSolves] === {},
    err["E5", cs, <|"Detail" ->
      "CppBackend persistent retained-local batch bridge is not available"|>]];
  batch = DiffExp2`CppBackend`RunPersistentLocalSolves[
    Lookup[prepared, "Request"],
    First[Lookup[prepared, "PersistentMetadata"]],
    Lookup[prepared, "LocalMetadata"], workerCount, equationOwner];
  If[FailureQ[batch] || !AssociationQ[batch] ||
      Lookup[batch, "status", "error"] =!= "ok",
    err["E5", cs, <|"BackendFailure" -> batch,
      "Detail" -> "persistent native regular basis batch failed"|>]];
  raw = Lookup[batch, "results", None];
  forbiddenPayloadKeys = Intersection[Keys[batch],
    {"assembled", "coefficients", "u", "validity"}];
  If[!ListQ[raw] || Length[raw] =!= d ||
      !StringQ[Lookup[batch, "session", None]] ||
      !StringQ[Lookup[batch, "chart", None]] ||
      Lookup[batch, "attempted", None] =!= d ||
      Lookup[batch, "succeeded", None] =!= d ||
      Lookup[batch, "failed", None] =!= 0 ||
      Lookup[batch, "requested_threads", None] =!= workerCount ||
      !IntegerQ[Lookup[batch, "worker_threads", None]] ||
      !TrueQ[1 <= batch["worker_threads"] <= workerCount] ||
      !TrueQ[Lookup[batch, "atomic_retention", False]] ||
      Lookup[batch, "json_coefficients", None] =!= 0 ||
      forbiddenPayloadKeys =!= {},
    If[ListQ[raw], Scan[
      If[AssociationQ[#] && StringQ[Lookup[#, "local", None]],
        Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]]] &, raw]];
    err["E6", cs, <|"BackendResponse" -> batch,
      "ForbiddenPayloadKeys" -> forbiddenPayloadKeys,
      "Detail" ->
        "native regular basis batch violated its ordered opaque-retention contract"|>]];
  cleanup[] := Scan[
    If[AssociationQ[#] && StringQ[Lookup[#, "local", None]],
      Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]]] &, raw];
  expectedLocalChart = If[equationOwner === Automatic, batch["chart"],
    Lookup[equationOwner, "NativeEquationOwner",
      Lookup[equationOwner, "EquationOwner", None]]];
  expectedLocalOperator = If[equationOwner === Automatic,
    First[prepared]["PersistentMetadata", "PreparedToken"],
    Lookup[equationOwner, "EquationIdentity",
      Lookup[equationOwner, "ChartIdentity", None]]];
  columns = Catch[
    MapThread[Function[{spec, response, basisIndex},
      Join[nativeLocalFamilyFinalize[physical, req, spec, response,
          batch["session"], expectedLocalChart, expectedLocalOperator],
        <|"Type" -> "DiffExp2NativeRegularBasisColumn",
          "BasisIndex" -> basisIndex|>]],
      {prepared, raw, Range[d]}],
    "DiffExp2Error", Function[{failure, errorTag},
      cleanup[]; Throw[failure, errorTag]]];
  sessions = DeleteDuplicates[Lookup[columns, "Session", None]];
  charts = DeleteDuplicates[Lookup[columns, "NativeChart", None]];
  If[Length[columns] =!= d || sessions === {None} ||
      Length[sessions] =!= 1 || charts === {None} || Length[charts] =!= 1 ||
      Lookup[columns, "BasisIndex", {}] =!= Range[d],
    cleanup[];
    err["E6", cs, <|"Sessions" -> sessions, "NativeCharts" -> charts,
      "BasisIndices" -> Lookup[columns, "BasisIndex", {}],
      "Detail" -> "retained monolithic regular basis did not bind one complete native chart"|>]];
  <|"Type" -> "DiffExp2NativeRegularBasis",
    "Session" -> First[sessions], "NativeChart" -> First[charts],
    "Columns" -> columns, "Dimension" -> d,
    "Chart" -> <|"Center" -> cs["Center"],
      "ChartMap" -> cs["ChartMap"], "Radius" -> cs["Radius"],
      "Prescriptions" -> cs["Prescriptions"]|>,
    "EpsWindow" -> epsWindow,
    "TWindow" -> <|"CompleteMax" -> req["TOrder"]|>,
    "NativeSummary" -> Join[KeyDrop[batch, {"results"}], <|
      "execution_capability" ->
        "retained-regular-monolithic-unit-basis-v2",
      "columns" -> d|>]|>];

(* Exact tile geometry needs every receiving equation owner before marching,
   but it does not need a local-solution slab.  Capture one validated ordinary
   (0,0,0) run to define the immutable value-solver prototype, then prepare
   only its physical monolithic chart.  The content token binds subsequent
   value or fallback-basis runs to that same native owner.  Singular receivers
   continue to use SCC composites. *)
PrepareNativeRegularBasisOwner[cs_Association, req_Association] := Module[
  {d = cs["SystemSize"],
   epsWindow = Lookup[req, "EpsWindow", None], epsMin, epsMax, physical,
   unitValues, prepared, owner, expectedIdentity, ownerTimingQ,
   registrationStarted},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "PrepareNativeRegularBasisOwner requires a regular chart"|>]];
  If[!AssociationQ[epsWindow] ||
      !IntegerQ[Lookup[epsWindow, "Min", None]] ||
      !IntegerQ[Lookup[epsWindow, "CompleteMax", None]],
    err["E8", cs, <|"Request" -> req,
      "Detail" -> "native regular owner requires a finite epsilon window"|>]];
  epsMin = epsWindow["Min"];
  epsMax = epsWindow["CompleteMax"];
  If[epsMin > 0 || epsMax < 0,
    err["E6", cs, <|"EpsWindow" -> epsWindow,
      "Detail" -> "native regular owner normalization requires eps^0 inside the requested window"|>]];
  physical = regularPhysicalChartSystem[cs];
  unitValues = Table[DiffExp2`EpsSeries`ESNew[epsMin,
      Table[If[component === 1 && power === 0, 1, 0],
        {power, epsMin, epsMax}]], {component, d}];
  prepared = prepareNativeLocalFamilyRun[physical, req,
    <|"a" -> 0, "b" -> 0, "p" -> 0|>, {unitValues}, True];
  ownerTimingQ = Environment["DE2_NATIVE_OWNER_TIMING"] === "1";
  registrationStarted = SessionTime[];
  owner = DiffExp2`CppBackend`PreparePersistentChart[
    prepared["Request"], prepared["PersistentMetadata"]];
  If[ownerTimingQ, Print[
    "DE2 NATIVE OWNER PHASE center=", InputForm[physical["Center"]],
    " phase=cppOwnerRegistration elapsedMs=",
    N[1000 (SessionTime[] - registrationStarted), 8],
    " dimension=", d]];
  If[FailureQ[owner] || !AssociationQ[owner],
    err["E5", cs, <|"BackendFailure" -> owner, "Detail" ->
      "persistent native regular chart-owner preparation failed"|>]];
  expectedIdentity = prepared["PersistentMetadata", "PreparedToken"];
  If[!StringQ[Lookup[owner, "Session", None]] ||
      !StringQ[Lookup[owner, "Chart", None]] ||
      !StringQ[Lookup[owner, "ChartIdentity", None]] ||
      owner["ChartIdentity"] =!= expectedIdentity,
    err["E6", cs, <|"PreparedOwner" -> owner,
      "ExpectedChartIdentity" -> expectedIdentity,
      "Detail" -> "prepared native regular owner lost its exact session/chart identity"|>]];
  <|"Type" -> "DiffExp2NativeRegularBasisOwner",
    "Session" -> owner["Session"],
    "NativeChart" -> owner["Chart"],
    "ChartIdentity" -> owner["ChartIdentity"],
    "ValueSolver" -> prepared["ValueSolver"]|>];

(* Production deferred transport can bind ordinary receiving equations before
   any recurrence frame exists.  The anchor has already selected and opened
   the exact persistent scalar session, so this form serializes only physical
   q/C, exact geometry, and one compact local-metadata prototype.  Keep the
   two-argument definition above unchanged for compatibility and explicit
   framed-owner tests. *)
PrepareNativeRegularBasisOwner[cs_Association, req_Association,
    anchor_Association] := Module[
  {d = cs["SystemSize"], epsWindow = Lookup[req, "EpsWindow", None],
   physical, physicalData, geometry, sessionStats, domain, wp,
   inputDigits, identityRecord, identity, checkpointIdentity,
   physicalPayload, localMetadata, owner, relativeAccuracy,
   physicalPreparation},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "PrepareNativeRegularBasisOwner requires a regular chart"|>]];
  If[!AssociationQ[epsWindow] ||
      !IntegerQ[Lookup[epsWindow, "Min", None]] ||
      !IntegerQ[Lookup[epsWindow, "CompleteMax", None]] ||
      epsWindow["Min"] > epsWindow["CompleteMax"] ||
      !IntegerQ[Lookup[req, "TOrder", None]] || req["TOrder"] < 0,
    err["E8", cs, <|"Request" -> req,
      "Detail" -> "lightweight regular owner requires finite ordered epsilon and Taylor windows"|>]];
  If[epsWindow["Min"] > 0 || epsWindow["CompleteMax"] < 0,
    err["E6", cs, <|"EpsWindow" -> epsWindow,
      "Detail" -> "native regular owner normalization requires eps^0 inside the requested window"|>]];
  sessionStats = DiffExp2`CppBackend`PersistentSessionCounters[anchor];
  domain = If[AssociationQ[sessionStats],
    Lookup[sessionStats, "domain", None], None];
  If[!MemberQ[{"acb", "rational"}, domain],
    err["E5", cs, <|"BackendResponse" -> sessionStats,
      "Detail" -> "lightweight regular owner requires an existing Acb or Rational anchor session"|>]];
  physical = regularPhysicalChartSystem[cs];
  physicalData = physicalClearedODEData[physical];
  geometry = cppPersistentGeometry[physical];
  wp = cfg["WorkingPrecision"];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  identityRecord = <|
    "schema" -> "diffexp2-frame-independent-regular-equation-identity-v1",
    "domain" -> domain, "dimension" -> d,
    "physical_identity" -> physicalData["Identity"],
    "center" -> RootReduce[physical["Center"]],
    "chart_map" -> RootReduce[physical["ChartMap", "Scale"]],
    "radius" -> If[physical["Radius"] === Infinity, Infinity,
      RootReduce[physical["Radius"]]],
    "prescriptions" -> Lookup[physical, "Prescriptions", {}]|>;
  identity = "de2-equation-" <>
    IntegerString[Hash[identityRecord, "SHA256"], 16, 64];
  checkpointIdentity = "de2-native-physical-value-prototype-" <>
    IntegerString[Hash[{identityRecord, req}, "SHA256"], 16, 64];
  relativeAccuracy = cppRegularValueRelativeAccuracyMaxExact[];
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    physicalPayload = cppPhysicalODEPayload[
      physicalData, identity, inputDigits, physical];
    localMetadata = cppNativeLocalMetadata[physical, 0, 0, 0,
      inputDigits, checkpointIdentity]];
  owner = DiffExp2`CppBackend`PreparePersistentRegularEquationOwner[
    anchor, <|
      "capability" ->
        "frame-independent-regular-physical-equation-owner-v1",
      "key" -> ("regular-equation:" <>
        IntegerString[Hash[identityRecord, "SHA256"], 16, 64]),
      "identity" -> identity, "dimension" -> d,
      "geometry" -> geometry, "physical_ode" -> physicalPayload,
      "relative_accuracy_max_exact" -> relativeAccuracy|>];
  If[FailureQ[owner] || !AssociationQ[owner] ||
      Lookup[owner, "Session", None] =!= Lookup[anchor, "Session",
        Lookup[anchor, "session", None]] ||
      !StringQ[Lookup[owner, "EquationOwner", None]] ||
      !StringStartsQ[owner["EquationOwner"], "eq:"] ||
      Lookup[owner, "EquationIdentity", None] =!= identity ||
      Lookup[owner, "PhysicalPayloadIdentity", None] =!=
        physicalData["Identity"],
    err["E5", cs, <|"BackendFailure" -> owner,
      "Detail" -> "frame-independent regular equation-owner preparation lost its exact session or q/C identity"|>]];
  physicalPreparation = makeRegularOwnerPhysicalPreparation[
    physical, req, identity, physicalData];
  <|"Type" -> "DiffExp2NativeRegularBasisOwner",
    "OwnerKind" -> "FrameIndependentRegularPhysicalEquation",
    "Session" -> owner["Session"],
    "NativeEquationOwner" -> owner["EquationOwner"],
    "ChartIdentity" -> owner["EquationIdentity"],
    "EquationIdentity" -> owner["EquationIdentity"],
    "PhysicalPayloadIdentity" -> owner["PhysicalPayloadIdentity"],
    "PhysicalPreparation" -> physicalPreparation,
    "ValueSolver" -> <|
      "schema" -> "diffexp2-native-ordinary-physical-value-solver-v1",
      "taylor_complete_max" -> req["TOrder"],
      "metadata" -> localMetadata,
      "relative_accuracy_max_exact" -> relativeAccuracy|>|>];

sccAggregateDiagnostics[records_List, seq_Association,
    recombine_List] := Module[{diags, getLists, compensated},
  diags = Select[records, AssociationQ];
  getLists[key_] := Flatten[
    Map[If[ListQ[#], #, {#}] &, Lookup[diags, key, {}]], 1];
  compensated = And @@ (TrueQ[Lookup[#,
      "PseudoCollisionsCompensated", True]] & /@ diags);
  <|"PseudoCollisionsHit" -> getLists["PseudoCollisionsHit"],
    "PseudoCompensations" -> getLists["PseudoCompensations"],
    "PseudoCollisionsCompensated" -> compensated,
    "AdaptiveLowerFrames" -> getLists["AdaptiveLowerFrames"],
    "SCCSolved" -> True,
    "SCCExecutionMode" -> "BlockSequentialStrict",
    "SCCCertificate" -> seq,
    "SCCRecombineGroups" -> recombine,
    "CouplingDepth" -> seq["CouplingDepth"]|>];

solveSCCBasisAtTop[cs_Association, req_Association, blockSystems_List,
    workTop_Integer] := Module[
  {seq = cs["IntegrationSequence"], nb, workReq, blockBases, columns = {},
   specs = {}, diagnostics = {}, recombine = {}, familyOffset = 0,
   columnOffset = 0, state, source, part, full, minDelivered, workTOrder},
  nb = Length[seq["Components"]];
  (* Block particulars can be individually ill-scaled even though their
     assembled physical column is well-conditioned.  Extra Taylor guard rows
     keep the mandatory intermediate residual probe clear of their finite
     tail.  Triangular recurrences make the retained rows independent of this
     halo; the final full-system proof is rerun after capping to the public N. *)
  workTOrder = sccWorkTOrder[cs, req];
  workReq = Join[req, <|"EpsWindow" ->
    Join[req["EpsWindow"], <|"CompleteMax" -> workTop|>],
    "TOrder" -> workTOrder|>];
  blockBases = Block[{$suppressIntermediateODEResidualChecks = True},
    If[cfg["RecurrenceBackend"] === "Cpp" &&
        !TrueQ[$disableGroupedSpectralTransform] &&
        cfg["Variables"] === {} &&
        Length[blockSystems] <= HomogeneousCacheCapacity[],
      (* Diagonal SCC homogeneous problems are independent and share the
         parent's SolveCacheTag.  Submit their native recurrence payloads as
         one batch; the surrounding transient cache then supplies the
         ordinary SolveHomogeneous reads below without persistent block
         entries.  Their per-block ODE proofs are subsumed by the mandatory
         original-system proof after all columns are assembled. *)
      PrewarmHomogeneousBatch[blockSystems, workReq],
      SolveHomogeneous[#, workReq] & /@ blockSystems]];
  Do[
    recombine = Join[recombine,
      sccBlockRecombineRecords[blockSystems[[block]], blockBases[[block]],
        columnOffset]];
    diagnostics = Append[diagnostics,
      blockBases[[block]]["Diagnostics"]];
    Do[
      state = <|block -> sccReframeLocalSolution[
        blockBases[[block]]["Columns"][[localColumn]], cs]|>;
      Do[
        source = sccBlockCouplingSource[cs, target, state];
        If[AssociationQ[source] && !sccSourceZeroQ[source],
          part = Block[{$suppressIntermediateODEResidualChecks = True},
            SolveParticular[blockSystems[[target]], source, workReq]];
          diagnostics = Append[diagnostics,
            Lookup[part, "Diagnostics", <||>]];
          AssociateTo[state, target -> sccReframeLocalSolution[part, cs]]],
        {target, block + 1, nb}];
      full = sccCombineBlockState[cs, state];
      AppendTo[columns, full];
      AppendTo[specs, Join[blockBases[[block]]["Specs"][[localColumn]], <|
        "Family" -> (familyOffset +
          blockBases[[block]]["Specs"][[localColumn]]["Family"]),
        "SCCBlock" -> block|>]],
      {localColumn, Length[blockBases[[block]]["Columns"]]}];
    familyOffset += Length[blockSystems[[block]]["Families"]];
    columnOffset += Length[blockBases[[block]]["Columns"]],
    {block, nb}];
  minDelivered = Min[# ["EpsWindow", "CompleteMax"] & /@ columns];
  <|"Columns" -> columns, "Specs" -> specs,
    "DiagnosticsRecords" -> diagnostics,
    "RecombineGroups" -> recombine,
    "MinDelivered" -> minDelivered|>];

solveSCCBasis[cs_Association, req_Association,
    blockSystems_List] := Module[
  {seq = cs["IntegrationSequence"], requested = req["EpsWindow", "CompleteMax"],
   workTop, attempt = 0, maxAttempts, built, deficit, columns, fs},
  workTop = requested + sccInitialWorkHalo[
    cs, blockSystems, sccWorkTOrder[cs, req]];
  maxAttempts = Length[blockSystems] + 3;
  While[True,
    attempt++;
    built = solveSCCBasisAtTop[cs, req, blockSystems, workTop];
    deficit = requested - built["MinDelivered"];
    If[deficit <= 0, Break[]];
    If[attempt >= maxAttempts,
      err["E6", cs, <|"RequestedCompleteMax" -> requested,
        "DeliveredCompleteMax" -> built["MinDelivered"],
        "WorkCompleteMax" -> workTop, "Attempts" -> attempt,
        "Detail" -> "SCC source propagation exhausted the deterministic epsilon halo"|>]];
    (* Every retry stays on the same strict recurrence path.  The measured
       honest deficit is monotone and contains no assumed-zero coefficients. *)
    workTop += Max[1, deficit]];
  columns = capTWindow[
      capWindow[cs, #, requested], req["TOrder"]] & /@ built["Columns"];
  fs = <|"Columns" -> columns, "Specs" -> built["Specs"],
    "Diagnostics" -> sccAggregateDiagnostics[
      built["DiagnosticsRecords"], seq, built["RecombineGroups"]]|>;
  ODEResidualCheck[cs, fs];
  fs];

sccZeroParticular[cs_Association, source_Association,
    req_Association] := Module[{nmax, min, max, d = cs["SystemSize"]},
  nmax = Min[req["TOrder"], source["TWindow", "CompleteMax"]];
  min = req["EpsWindow", "Min"]; max = req["EpsWindow", "CompleteMax"];
  <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"],
    "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
      "Coeffs" -> ConstantArray[0, {max - min + 1, nmax + 1, d}]|>},
    "EpsWindow" -> <|"Min" -> min, "CompleteMax" -> max|>,
    "TWindow" -> <|"CompleteMax" -> nmax|>,
    "ErrorEstimate" -> ConstantArray[0, max - min + 1],
    "Prescriptions" -> cs["Prescriptions"]|>];

solveSCCParticularAtTop[cs_Association, source_Association,
    req_Association, blockSystems_List, workTop_Integer] := Module[
  {seq = cs["IntegrationSequence"], nb, state = <||>, external, coupling,
   total, part, full, workReq},
  nb = Length[seq["Components"]];
  workReq = Join[req, <|"EpsWindow" ->
    Join[req["EpsWindow"], <|"CompleteMax" -> workTop|>],
    "TOrder" -> sccWorkTOrder[cs, req]|>];
  Do[
    external = sccProjectSource[source, seq["Components"][[block]]];
    coupling = sccBlockCouplingSource[cs, block, state];
    total = sccCombineSources[cs, {external, coupling},
      Length[seq["Components"][[block]]]];
    If[AssociationQ[total] && !sccSourceZeroQ[total],
      part = Block[{$suppressIntermediateODEResidualChecks = True},
        SolveParticular[blockSystems[[block]], total, workReq]];
      AssociateTo[state, block -> sccReframeLocalSolution[part, cs]]],
    {block, nb}];
  If[state === <||>, Return[<|"Solution" ->
      sccZeroParticular[cs, source, workReq],
    "MinDelivered" -> workTop|>, Module]];
  full = sccCombineBlockState[cs, state];
  <|"Solution" -> full,
    "MinDelivered" -> full["EpsWindow", "CompleteMax"]|>];

solveSCCParticular[cs_Association, source_Association, req_Association,
    blockSystems_List] := Module[
  {requested = req["EpsWindow", "CompleteMax"], workTop, attempt = 0,
   maxAttempts, built, deficit, full, sourceTop,
   requestedWindowNotReached = None},
  (* Validate the physical source dimension before slicing it. *)
  sccSourceToLocalSolution[cs, source, cs["SystemSize"]];
  If[sccSourceZeroQ[source], Return[sccZeroParticular[cs, source, req], Module]];
  workTop = requested + sccParticularInitialWorkHalo[
    cs, blockSystems, sccWorkTOrder[cs, req]];
  sourceTop = source["EpsWindow", "CompleteMax"];
  maxAttempts = Length[blockSystems] + 3;
  While[True,
    attempt++;
    built = solveSCCParticularAtTop[
      cs, source, req, blockSystems, workTop];
    deficit = requested - built["MinDelivered"];
    If[deficit <= 0, Break[]];
    (* An inhomogeneous finite source is allowed to deliver a narrower honest
       epsilon window than requested.  Retry only while unused source halo can
       actually fund a wider solve; asking beyond its certified top cannot
       manufacture future coefficients. *)
    If[attempt >= maxAttempts || workTop >= sourceTop,
      requestedWindowNotReached = <|
        "RequestedCompleteMax" -> requested,
        "DeliveredCompleteMax" -> built["MinDelivered"],
        "SourceCompleteMax" -> sourceTop,
        "WorkCompleteMax" -> workTop, "Attempts" -> attempt|>;
      Break[]];
    workTop += Min[Max[1, deficit], sourceTop - workTop]];
  full = capTWindow[
    capWindow[cs, built["Solution"], requested], req["TOrder"]];
  If[full["EpsWindow", "CompleteMax"] < full["EpsWindow", "Min"],
    err["E6", cs, <|"EpsWindow" -> full["EpsWindow"],
      "Detail" -> "SCC particular propagation produced an empty epsilon window"|>]];
  full = Join[full, <|"Diagnostics" -> Join[
    Lookup[full, "Diagnostics", <||>], <|
      "SCCSolved" -> True,
      "SCCExecutionMode" -> "BlockSequentialStrict",
      "SCCCertificate" -> cs["IntegrationSequence"],
      "CouplingDepth" -> cs["IntegrationSequence", "CouplingDepth"]|>]|>];
  If[AssociationQ[requestedWindowNotReached],
    full = Join[full, <|"Diagnostics" -> Join[
      Lookup[full, "Diagnostics", <||>],
      <|"RequestedWindowNotReached" -> requestedWindowNotReached|>]|>]];
  ODEResidualCheck[cs, full, source];
  full];

SolveChart[cs_Association, req_Association, source_:None] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], basis, part, multiQ},
  multiQ = AssociationQ[seq] && Length[seq["Components"]] > 1;
  If[!multiQ,
    basis = SolveHomogeneous[cs, req];
    part = If[source === None, None, SolveParticular[cs, source, req]];
    Return[<|"Basis" -> basis, "Particular" -> part,
      "CouplingDepth" -> 0|>, Module]];
  basis = SolveHomogeneous[cs, req];
  part = If[source === None, None, SolveParticular[cs, source, req]];
  <|"Basis" -> basis, "Particular" -> part,
    "CouplingDepth" -> seq["CouplingDepth"],
    "IntegrationSequence" -> seq|>];

(* ---- value-vector propagation (regular charts; prototype) ----
   The incoming VALUE at the chart center is the t^0 Cauchy datum of the
   transported solution: ONE runRecursion with init = vals replaces the
   d-column basis + MatchWeights + CombineLocalSolutions of the basis
   path (the old engine's Currbcs chaining, DiffExp/Transport.m).
   Regular charts only: the single (0,0,0) family makes every n >= 1
   step CASE T (dA = n != 0).  Negative-valued Taylor multipliers can
   nevertheless lower the honest CompleteMax and are included in the work
   frame.  Invariants kept:
   exact (0,0,0) sector tag, honest window (capped at the incoming
   CompleteMax), no matching solve to go ambiguous, and the always-on
   ODE residual check runs on the propagated solution itself. *)
regularPhysicalChartSystem[cs_Association] := Module[
  {d = cs["SystemSize"]},
  If[TrueQ[Lookup[cs, "SCCSkeleton", False]] &&
      TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    Join[cs, <|"SCCSkeleton" -> False,
      (* Recombining the exact SCC skeleton must not invent a second session
         identity.  Use the same collision-certified parent system key as
         ordinary and composite charts in a mixed retained atlas. *)
      "SolveCacheTag" -> Lookup[cs, "SystemClearKey",
        Lookup[cs, "SystemHash", None]],
      "ThetaMatrix" -> cs["ThetaOriginal"],
      "Gauge" -> IdentityMatrix[d],
      "GaugeInverse" -> IdentityMatrix[d],
      "Residue" -> ConstantArray[0, {d, d}],
      "V" -> IdentityMatrix[d], "VInv" -> IdentityMatrix[d],
      "Families" -> {<|
        "Roots" -> Table[<|"a" -> 0, "b" -> 0,
          "BlockSize" -> 1|>, d],
        "ColumnRange" -> Range[d], "Collisions" -> {},
        "CollisionDepth" -> 0|>}|>],
    cs]];

SolveNativeValueRegular[cs_Association, req_Association,
    vals_List] := Module[
  {physical, d = cs["SystemSize"], vMin, vCM, requestedMin,
   requestedMax, deliveredMax, nativeReq, result},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "SolveNativeValueRegular requires a regular chart"|>]];
  If[Length[vals] =!= d || !AllTrue[vals, DiffExp2`EpsSeries`ESQ],
    err["E8", cs, <|"Components" -> Length[vals], "Dimension" -> d,
      "Detail" -> "native value vector must be d EpsSeries components"|>]];
  If[!AssociationQ[Lookup[req, "EpsWindow", None]] ||
      !IntegerQ[Lookup[req, "TOrder", None]] || req["TOrder"] < 0,
    err["E8", cs, <|"Request" -> req,
      "Detail" -> "native value transport requires finite Taylor and epsilon windows"|>]];
  vMin = Min[esMin /@ vals];
  vCM = Min[esCM /@ vals];
  requestedMin = Lookup[req["EpsWindow"], "Min", None];
  requestedMax = Lookup[req["EpsWindow"], "CompleteMax", None];
  If[!IntegerQ[requestedMin] || !IntegerQ[requestedMax] ||
      requestedMin > requestedMax,
    err["E8", cs, <|"Request" -> req,
      "Detail" -> "native value transport epsilon window is malformed"|>]];
  deliveredMax = Min[requestedMax, vCM];
  If[deliveredMax < requestedMin,
    err["E6", cs, <|"IncomingCompleteMax" -> vCM,
      "RequestedMin" -> requestedMin,
      "Detail" -> "incoming value has no honest overlap with the requested epsilon window"|>]];
  nativeReq = Join[req, <|"EpsWindow" -> <|
      "Min" -> Min[requestedMin, vMin],
      "CompleteMax" -> deliveredMax|>|>];
  physical = regularPhysicalChartSystem[cs];
  result = SolveNativeLocalFamily[physical, nativeReq,
    <|"a" -> 0, "b" -> 0, "p" -> 0|>, {vals}];
  Join[result, <|"Type" -> "DiffExp2NativeValueRegular",
    "IncomingEpsWindow" -> <|"Min" -> vMin,
      "CompleteMax" -> vCM|>|>]];

SolveValueRegular[cs_Association, req_Association, vals_List] := Module[
  {d = cs["SystemSize"], nmax, vMin, vCM, fb, wideTop, Wd, prep, rec, ls,
   symbolic, poleDepth, numericInputQ, phaseQ, phaseTime, phase,
   preparedNums, vPrep},
  If[TrueQ[Lookup[cs, "SCCSkeleton", False]] &&
      TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    (* At an ordinary point the exact residue is zero, so the physical frame
       itself is already a valid identity spectral/gauge frame.  Materialize
       that frame only for this value recursion: the cached SCC envelope stays
       skeletal and never pays global ChartIndicial/spectral preparation. *)
    Return[SolveValueRegular[regularPhysicalChartSystem[cs], req, vals],
      Module]];
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "SolveValueRegular requires a regular chart (pole order 0); singular charts keep the basis+matching path"|>]];
  If[Length[vals] =!= d || !AllTrue[vals, DiffExp2`EpsSeries`ESQ],
    err["E8", cs, <|"Components" -> Length[vals], "Dimension" -> d,
      "Detail" -> "value vector must be d EpsSeries components"|>]];
  nmax = req["TOrder"];
  vMin = Min[esMin /@ vals]; vCM = Min[esCM /@ vals];
  numericInputQ = AnyTrue[vals, Function[v,
    !FreeQ[Table[esCoeff[v, k], {k, esMin[v], esCM[v]}],
      _?InexactNumberQ]]];
  phaseQ = Environment["DEBUG_SOLVE_PHASES"] === "1";
  phaseTime = SessionTime[];
  phase[label_String, meta_:None] := If[phaseQ, Module[{now = SessionTime[]},
    Print["SOLVEPHASE center=", cs["Center"], " phase=", label,
      " dt=", N[now - phaseTime, 6],
      If[meta === None, "", " meta=" <> ToString[meta, InputForm]]];
    phaseTime = now]];
  If[phaseQ,
    Print["SOLVEPHASE center=", cs["Center"], " phase=start",
      " nmax=", nmax, " inputNumeric=", numericInputQ]];
  If[AllTrue[vals, Function[v,
      AllTrue[Table[esCoeff[v, k], {k, esMin[v], esCM[v]}], # === 0 &]]],
    (* exactly-zero incoming value: the propagated solution is zero *)
    Return[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"],
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> Table[0, {vCM - vMin + 1}, {nmax + 1}, {d}]|>},
      "EpsWindow" -> <|"Min" -> vMin, "CompleteMax" -> vCM|>,
      "TWindow" -> <|"CompleteMax" -> nmax|>,
      "ErrorEstimate" -> ConstantArray[0, vCM - vMin + 1],
      "Prescriptions" -> cs["Prescriptions"]|>]];
  (* frame: the SolveHomogeneous work-window shape with Pmax = cdMax = 0
     (regular chart), based on the INCOMING window instead of req *)
  symbolic = clearedSymbolic[cs];
  phase["cleared-symbolic"];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  phase["pole-depth"];
  fb = Min[vMin, 0] - 2;
  wideTop = vCM + 2 - Min[0, vMin - 2];
  fb = Min[fb, vMin - poleDepth];
  wideTop = Max[wideTop, vCM + poleDepth];
  Wd = wideTop - fb + 1;
  prep = prepareCleared[cs, fb, Wd, symbolic];
  If[phaseQ,
    preparedNums = Cases[{prep["dL"], prep["NhatSp"],
        prep["NhatRationalDenominators"]}, _?NumericQ, Infinity]];
  phase["prepare-cleared", <|"fb" -> fb, "width" -> Wd,
    "exactNumericCount" -> If[phaseQ,
      Count[preparedNums, z_ /; Precision[z] === Infinity], Missing["NotScanned"]],
    "inexactCount" -> If[phaseQ,
      Count[preparedNums, _?InexactNumberQ], Missing["NotScanned"]]|>];
  vPrep = If[TrueQ[$disableGroupedSpectralTransform], Automatic,
    prepareFramedMatrix[cs["V"], DiffExp2`Config`CanonicalEps[],
      fb, Wd, cs]];
  rec = runRecursion[cs, prep, 0, 0, 0, nmax, None, fb, Wd, {vals},
    If[cfg["RecurrenceBackend"] === "Cpp", vPrep, Automatic]];
  phase["run-recursion"];
  ls = assembleSolution[cs, 0, 0, rec, nmax, vPrep];
  ls = capWindow[cs, ls, vCM];
  phase["assemble-cap"];
  ODEResidualCheck[cs, ls];
  phase["ode-residual"];
  ls];

(* honesty cap: nothing above the incoming CompleteMax is complete (the
   frame's buffer slots above it absorb eps-shifts but are not delivered) *)
capWindow[cs_, ls_, capCM_] := Module[
  {kmin = ls["EpsWindow", "Min"], kmax = ls["EpsWindow", "CompleteMax"],
   keep, ncols, ncomp},
  If[kmax <= capCM, Return[ls]];
  If[capCM < kmin,
    {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
    Return[Join[ls, <|
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> ConstantArray[0, {1, ncols, ncomp}]|>},
      "EpsWindow" -> <|"Min" -> capCM, "CompleteMax" -> capCM|>,
      "ErrorEstimate" -> {0}|>], Module]];
  keep = capCM - kmin + 1;
  Join[ls, <|
    "Sectors" -> Map[Append[#, "Coeffs" -> Take[#["Coeffs"], keep]] &,
      ls["Sectors"]],
    "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> capCM|>,
    "ErrorEstimate" -> If[ListQ[ls["ErrorEstimate"]],
      Take[PadRight[ls["ErrorEstimate"], kmax - kmin + 1], keep],
      ls["ErrorEstimate"]]|>]];

capTWindow[ls_Association, cap_Integer] := Module[
  {current = ls["TWindow", "CompleteMax"]},
  If[current <= cap, Return[ls, Module]];
  Join[ls, <|
    "Sectors" -> Map[Append[#, "Coeffs" ->
      #["Coeffs"][[All, 1 ;; cap + 1, All]]] &, ls["Sectors"]],
    "TWindow" -> <|"CompleteMax" -> cap|>|>]];

(* ---- residual check ---- *)

(* Residual-check numeric handoff (the Transport`numHandoff policy): the
   evaluated f / theta-f VALUES feeding the d^2 esTimes grid numericize at
   2x WP.  Exact evaluation outputs otherwise compound into exact-rational
   giants across the grid (measured: the check was ~3.1 s of a 13 s chart
   at d = 7, Docs/SpeedIdeas.md §2).  Values only — the residual compare
   below is numeric by construction and uses componentwise uncertainty
   bounds; windows untouched. *)
numV[s_] := esNew[esMin[s],
  N[Table[esCoeff[s, k], {k, esMin[s], esCM[s]}],
    DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]];

ODEResidualCheck[cs_Association, sol_Association, source_:None, probe_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], t0, sols, maxRel = 0,
   Bsub, bMinVal, bt0Cache = <||>, bt0For, win, rtol, formalData,
   formalSymbols, candidateRules, regularSpecializationQ, sampleCandidates,
   sampleRules, sampleResiduals, requiredSamples = 5,
   regulatorProbeValues},
  (* SCC diagonal bases and coupling particulars are implementation
     intermediates.  Their assembled columns receive one mandatory residual
     proof against the original full system; repeating the same proof on
     every block/source is mathematically redundant and dominated L3 wall
     time.  The dynamic flag is scoped only around those internal calls. *)
  If[TrueQ[$suppressIntermediateODEResidualChecks], Return[0, Module]];
  rtol = DiffExp2`Tolerances`Tol["ResidTol"];
  sols = If[KeyExistsQ[sol, "Columns"], sol["Columns"], {sol}];
  (* Exact analytic-regulator coefficients remain formal throughout the
     recurrence.  The ordinary numeric residual handoff below intentionally
     N[]s coefficient arrays; on a formal field that turns harmless Taylor
     truncation tails into nonnumeric expressions such as 10^-30 rho^n.
     Spot-check such a solution at several deterministic small rational field
     specializations instead.  Candidate rules that land on a coefficient
     pole are skipped.  This specializes only the independent residual check
     -- never the solve, epsilon lattice, indicial decisions, or delivered
     symbolic coefficients.  Like the ordinary one-t probe below this is a
     diagnostic spot-check, not an identity proof in the regulator field. *)
  formalData = {Lookup[cs, "ThetaOriginal", cs["ThetaMatrix"]],
    Map[{#["a"], #["b"], #["Coeffs"]} &,
      Flatten[Lookup[sols, "Sectors", {}]]],
    If[source === None, {}, Lookup[source, "Sectors", {}]]};
  formalSymbols = SortBy[DeleteDuplicates[Cases[formalData,
    s_Symbol /; Context[s] === "Global`" && s =!= eps && s =!= t &&
      !NumericQ[s], Infinity]], SymbolName];
  If[formalSymbols =!= {},
    (* Include separated O(1) points rather than only reciprocals tending to
       zero: a valid rational coefficient can have a pole near zero and make
       a finite Taylor residual spuriously ill-conditioned at every such
       sample. The first legacy point is retained as a pole-skip regression. *)
    regulatorProbeValues = {1/43, -1, 1, -2, 2, -1/2, 1/2,
      -3/2, 3/2, -1/3, 1/3, -3, 3, -2/3, 2/3};
    candidateRules[sample_Integer] := MapIndexed[
      (#1 -> If[sample <= Length[regulatorProbeValues],
        regulatorProbeValues[[1 + Mod[sample - 1 + 2 (First[#2] - 1),
          Length[regulatorProbeValues]]]],
        (-1)^(sample + First[#2])*
          Prime[10 + sample + 2 First[#2]]/
          Prime[11 + sample + 3 First[#2]]]) &, formalSymbols];
    regularSpecializationQ[rules_List] := Module[{specialized},
      specialized = Quiet[Check[formalData /. rules, $Failed]];
      specialized =!= $Failed &&
        FreeQ[specialized, Indeterminate | ComplexInfinity |
          DirectedInfinity[___] | _Failure]];
    sampleCandidates = Select[Range[1, 128],
      regularSpecializationQ[candidateRules[#]] &, requiredSamples];
    If[Length[sampleCandidates] < requiredSamples,
      err["E7", cs, <|"Regulators" -> formalSymbols,
        "RegularSamplesFound" -> Length[sampleCandidates],
        "SamplesRequired" -> requiredSamples,
        "Detail" -> "could not find enough regular rational specializations for the symbolic ODE residual spot-check"|>]];
    sampleRules = candidateRules /@ Take[sampleCandidates, requiredSamples];
    sampleResiduals = Table[With[{rules = sampleRules[[sample]]},
      ODEResidualCheck[cs /. rules, sol /. rules,
        If[source === None, None, source /. rules], probe]],
      {sample, requiredSamples}];
    Return[Max[sampleResiduals], Module]];
  (* truncation-aware probe: the common sector prefactor t^a cancels in the
     relative residual, so the degree-N tail scales as (t0/R)^(N+1). *)
  t0 = If[probe === Automatic,
    Module[{nmaxS, frac, fracRaw, raw},
      nmaxS = Min[#["TWindow", "CompleteMax"] & /@ sols];
      fracRaw = N[(rtol/100)^(1/(nmaxS + 1)), 30];
      (* A fixed absolute Rationalize tolerance rounded the small-order
         probe to exactly zero (for example nmax=0..2 at ResidTol=10^-10),
         sending singular sectors through the forbidden origin limit. *)
      frac = Min[1/4, Rationalize[fracRaw, Abs[fracRaw]/50]];
      raw = cs["Radius"]*frac*(2 + Mod[Hash[cs["Center"]], 5])/10;
      (* a SIMPLE RATIONAL probe: algebraic radii otherwise drag every
         downstream evaluation into exact algebraic arithmetic *)
      Rationalize[N[raw, 20], N[raw, 20]/50]],
    probe];
  (* B(t0) is column-independent: hoist the d^2 exact substitutions out of
     the per-column loop (they dominated the check at d = 7: ~7x recompute)
     and memoize the eps-frame expansion per (k1, k2) window — columns
     share windows in the common case.  Values are identical per column;
     pure cost hoist, no honesty change. *)
  Bsub = Map[Together[# /. t -> t0] &,
    Lookup[cs, "ThetaOriginal", cs["ThetaMatrix"]], {2}];
  bMinVal = Module[{v = Select[
      epsValuation[#, eps] & /@ Flatten[Bsub], IntegerQ]},
    If[v === {}, 0, Min[0, Min[v]]]];
  bt0For[k1_, k2_] := Module[
    {opTop = Max[0, k2 - k1], key},
    key = {bMinVal, opTop};
    If[!KeyExistsQ[bt0Cache, key],
      bt0Cache[key] = Map[Module[
          {v = epsValuation[#, eps], fl},
        If[v === Infinity || v > opTop, esZero[opTop],
          fl = ratEpsList[#, eps, bMinVal, opTop - bMinVal + 1];
          (* Leading frame zeros are padding, not part of the operator's
             mathematical Laurent support. *)
          esTrim[esNew[bMinVal, fl]]]] &, Bsub, {2}]];
    bt0Cache[key]];
  Do[Module[{f, df, lhs, rhs, srcv, k1, k2, scale, Bt0},
    If[Environment["DEBUG_RESID"] === "1",
      Print["RESID col window=", ls["EpsWindow"], " tags=",
        {#["a"], #["b"], #["p"]} & /@ ls["Sectors"]]];
    f = numV[DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
      "UsePade" -> False, "ComputeTailEstimates" -> False]["Value"]];
    df = numV[DiffExp2`SectorSeries`EvaluateLocalSolution[
      DiffExp2`SectorSeries`DifferentiateLocalSolution[ls], t0,
      "UsePade" -> False, "ComputeTailEstimates" -> False]["Value"]];
    (* Coefficients below an EpsSeries Min are certified zero.  Residuals
       must therefore start at the UNION of supports, not their intersection;
       otherwise a low-order B.f term can be hidden by a zero derivative. *)
    k1 = Min[esMin[f], esMin[df]];
    k2 = Min[esCM[f], esCM[df]];
    Bt0 = bt0For[esMin[f], k2];
    (* theta f = t f' *)
    Module[{thetaF, Bf, resid},
      thetaF = esScale[t0, df];
      Bf = Module[{fv = f, comps},
        comps = Length[esCoeff[fv, esMin[fv]]];
        Table[Module[{s = None},
          Do[Module[{cESr},
            If[!AllTrue[Bt0[[r, c]]["Coeffs"], # === 0 &],
              cESr = esNew[esMin[fv],
                Table[esCoeff[fv, k][[c]], {k, esMin[fv], esCM[fv]}]];
              Module[{term = esTimes[Bt0[[r, c]], cESr]},
                s = If[s === None, term, esAdd[s, term]]]]],
            {c, comps}];
          If[s === None, esZero[k2], s]], {r, comps}]];
      srcv = If[source =!= None && AssociationQ[source] && KeyExistsQ[source, "Sectors"],
        Module[{sLS = Join[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
            "Radius" -> cs["Radius"], "ErrorEstimate" -> 0,
            "Prescriptions" -> cs["Prescriptions"],
            "EpsWindow" -> source["EpsWindow"], "TWindow" -> source["TWindow"]|>,
            <|"Sectors" -> source["Sectors"]|>]},
          DiffExp2`SectorSeries`EvaluateLocalSolution[sLS, t0,
            "UsePade" -> False, "ComputeTailEstimates" -> False]["Value"]],
        None];
      Do[Module[{tf = esCoeff[thetaF, k], bfk, sv, residualEntries,
          badEntries, rv, sc},
        bfk = Table[If[esMin[Bf[[r]]] <= k <= esCM[Bf[[r]]], esCoeff[Bf[[r]], k], 0],
          {r, Length[Bf]}];
        sv = If[srcv =!= None && esMin[srcv] <= k <= esCM[srcv],
          esCoeff[srcv, k], ConstantArray[0, Length[Bf]]];
        residualEntries = Flatten[{tf - bfk - sv}];
        badEntries = Select[residualEntries,
          !NumericQ[#] &&
            !(FreeQ[#, _?InexactNumberQ] && TrueQ[PossibleZeroQ[#]]) &];
        If[badEntries =!= {},
          err["E7", cs, <|"Probe" -> t0, "EpsOrder" -> k,
            "Residual" -> First[badEntries],
            "Detail" -> "ODE residual contains uncertified nonnumeric content"|>]];
        rv = Select[residualEntries, NumericQ];
        sc = Max[1, Sequence @@ (numMag[#, 20] & /@
          Select[Flatten[{tf, bfk, sv}], NumericQ])];
        Module[{rk = Max[0, Sequence @@
              (Last[DiffExp2`Tolerances`NumericMagnitudeBounds[#, 20]] & /@ rv)]/sc},
          If[Environment["DEBUG_RESID"] === "1" && rk > 10^-12,
            Print["RESID k=", k, " rel=", N[rk, 4]]];
          maxRel = Max[maxRel, rk]]],
        {k, Min[k1, esMin[thetaF], Min[esMin /@ Bf],
            If[srcv === None, k1, esMin[srcv]]],
          Min[k2, esCM[thetaF], Min[esCM /@ Bf],
            If[srcv === None, k2, esCM[srcv]]]}]]],
    {ls, sols}];
  If[maxRel > rtol,
    err["E7", cs, <|"Probe" -> t0, "MaxRelativeResidual" -> N[maxRel, 6],
      "ResidTol" -> rtol, "Detail" -> "ODE residual check failed"|>]];
  maxRel];

End[];
EndPackage[];
