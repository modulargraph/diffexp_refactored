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
SolveValueRegular::usage = "SolveValueRegular[chartSystem, req, vals] propagates an incoming VALUE vector (one EpsSeries per component: the solution value AT THE CHART CENTER t = 0) through a REGULAR chart with ONE d-dimensional recursion (init = vals); no basis, no matching. The delivered eps-window is capped by the incoming window. Loud error on non-regular charts. (Value-transport prototype; see Docs/PerfGapAnalysis.md lever 1.)";
SolveNativeValueRegular::usage = "SolveNativeValueRegular[chartSystem, req, vals] propagates a regular center-value vector as one retained C++ local solution. It accepts an SCC skeleton by materializing the exact identity physical frame for this single value recursion, caps delivery to the honest incoming epsilon window, and never returns a coefficient tensor.";
SolveNativeLocalFamily::usage = "SolveNativeLocalFamily[chartSystem, req, <|\"a\"->a,\"b\"->b,\"p\"->p|>, init] runs one uncompensated homogeneous family through the persistent C++ solver and returns an opaque native handle record, never a Wolfram coefficient tensor. init is the same (p+1)-by-d EpsSeries ladder accepted by the framed recurrence. This narrow migration seam requires an identity gauge, grouped native assembly, no unresolved analytic regulators, and no pseudo-resonant family collisions; general transport continues to use SolveHomogeneous/SolveParticular.";
PrepareSCCCouplingMatrix::usage = "PrepareSCCCouplingMatrix[sccChartSystem, sourceBlock, targetBlock, sourceShape, serialization] prepares one exact cross-SCC ThetaOriginal block as a deterministic JSON-ready sparse rational-multiplier matrix. serialization is Automatic (the active C++ serialization Block) or the exact field <|\"domain\"->...,\"symbols\"->{...}|>. Signed epsilon shifts are preserved; execution later proves the requested/work halo contract.";
PrepareNativeSCCComposite::usage = "PrepareNativeSCCComposite[sccChartSystem, req] captures (without executing) the ordinary grouped native homogeneous requests for every supported diagonal SCC block, prepares their strict typed persistent composite manifest, and returns the opaque C++ SCC handle record. This first slice is an explicit preparation API only; SolveHomogeneous does not dispatch through it.";
SolveNativeSCCBasisColumn::usage = "SolveNativeSCCBasisColumn[sccChartSystem, req, seedBlock, seedLocalComponent:1] executes one strict regular exact-Rational or Acb block-DAG SCC basis column, or an exact-Rational regular-singular Jordan column, through an already captured persistent composite and returns an opaque native local handle record without coefficient tensors. seedBlock and seedLocalComponent are one-based; the three-argument scalar-v1 call is unchanged. This explicit migration seam is not yet used by SolveHomogeneous or transport.";
SolveNativeSCCBasis::usage = "SolveNativeSCCBasis[sccChartSystem, req, threads:Automatic] executes the complete physical SCC basis as one ordered native column batch, retaining every column atomically and returning opaque handles sorted by physical basis index. No coefficient tensor crosses the bridge.";
SolveNativeRegularBasis::usage = "SolveNativeRegularBasis[chartSystem, req, threads:Automatic] returns a complete retained basis for any regular chart. Multi-block SCC envelopes use the ordered native SCC batch; a single strongly connected block uses the same retained full-system recurrence with exact eps^0 unit seeds. No coefficient tensor crosses the bridge.";
ClearSolveCaches::usage = "ClearSolveCaches[] empties the PrepareChart, exact-SCC-structure, exact-clearing, rational-multiplier, SolveHomogeneous, and native SCC composite memo caches, then closes persistent native sessions. Called by API`LoadSystem; the SolveHomogeneous cache additionally self-flushes whenever the chart's SystemHash changes and is entry-capped.";
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
   integrationSequence, thetaOriginal, regularOriginal},
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
      With[{v = tVal[#, t]}, v === Infinity || TrueQ[v >= 0]] &];
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
  idata = DiffExp2`Indicial`ChartIndicial[Achart, t, eps,
    <|"Name" -> Lookup[chart, "Name", "chart@" <> ToString[x0, InputForm]],
      "Center" -> x0, "Variable" -> sys["Variable"]|>];
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
preparedNumericizationRequiredQ[e_] :=
  exactNonRationalNumberQ[e] ||
  (TrueQ[$numericizeAllPreparedNumbers] && e =!= 0 && NumericQ[e] &&
    Precision[e] === Infinity);

groupedEpsExactSafeQ[e_] := e === 0 ||
  (ByteCount[e] <= 500 && !preparedNumericizationRequiredQ[e]);

preparedEpsCoefficient[e_] := Module[{wp2},
  If[groupedEpsExactSafeQ[e], Return[e]];
  If[!NumericQ[e], Return[e]];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*
    DiffExp2`Config`CFG["WorkingPrecision"];
  Block[{$MaxExtraPrecision = Max[$MaxExtraPrecision,
      DiffExp2`Tolerances`$MaxExtraPrecisionValue, 2 wp2]},
    Chop[SetPrecision[N[RootReduce[e], wp2], wp2],
      DiffExp2`Tolerances`Tol["ChopFloor"]]]];

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
clearedSymbolicLegacy[cs_] := Module[
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
  {cached, input, x, A, den, denCoeffs, denContent, num, phaseQ, phaseTime},
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
  cached = <|"Variable" -> x, "Denominator" -> den,
    "Numerator" -> num, "Dimension" -> Length[A]|>;
  AssociateTo[$globalClearedCache, systemKey -> cached];
  cached];

regularClearedFromGlobal[cs_Association] := Module[
  {systemKey = cs["SystemClearKey"], key, cached, global, x,
   t = cs["ChartVar"], center = cs["Center"], beta = cs["ChartMap", "Scale"],
   den, denCoeffs, denContent, num, d0, dD, dN, dExpr, NhatExpr,
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
  den = Cancel[Together[global["Denominator"] /. x -> center + beta*t]];
  denCoeffs = Select[CoefficientList[den, t], !zeroCanQ[#] &];
  If[denCoeffs === {},
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "affine-shifted global denominator is identically zero"|>]];
  denContent = If[Length[denCoeffs] === 1, First[denCoeffs],
    Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]];
  den = Cancel[Together[den/denContent]];
  num = Map[Cancel[Together[beta*t*(# /. x -> center + beta*t)/denContent]] &,
    global["Numerator"], {2}];
  phase["affine-pair"];
  d0 = Together[den /. t -> 0];
  If[zeroQ[d0 /. DiffExp2`Config`CanonicalEps[] -> 0],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  dExpr = Table[Cancel[Together[Coefficient[den, t, j]]], {j, 0, dD}];
  NhatExpr = Table[
    Map[Cancel[Together[#]] &, Map[Coefficient[#, t, j] &, num, {2}], {2}],
    {j, 0, dN}];
  phase["coefficients"];
  cached = <|"dExpr" -> dExpr, "NhatExpr" -> NhatExpr,
    "dD" -> dD, "dN" -> dN|>;
  AssociateTo[$chartClearedCache, key -> cached];
  cached];

clearedSymbolic[cs_Association] := If[
  !TrueQ[$disableGlobalClearedHoist] && regularIdentityFrameQ[cs] &&
    KeyExistsQ[cs, "SystemClearKey"],
  regularClearedFromGlobal[cs],
  clearedSymbolicLegacy[cs]];

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
    gaugeCoefficientData[cs, nmax]["EpsValuation"];

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

cppPersistentGeometry[cs_Association] := Join[<|
    "center_exact" -> ToString[cs["Center"], InputForm],
    "scale_exact" -> ToString[cs["ChartMap", "Scale"], InputForm],
    "infinite_radius" -> TrueQ[cs["Radius"] === Infinity],
    "prescriptions" ->
      (cppPersistentPrescription /@ Lookup[cs, "Prescriptions", {}])|>,
  If[TrueQ[cs["Radius"] === Infinity], <||>,
    <|"radius_exact" -> ToString[cs["Radius"], InputForm]|>]];

(* Exact producer certificates for the deliberately narrow first composite
   slice.  These are structural predicates, never numerical enclosure tests.
   C++ independently proves identity_v from the retained assembly operator;
   the other facts remain collision-bound until native execution rechecks
   them. *)
sccExactIdentityMatrixQ[matrix_, dimension_Integer] := Module[{delta},
  If[!MatrixQ[matrix] || Dimensions[matrix] =!= {dimension, dimension} ||
      !FreeQ[matrix, _?InexactNumberQ], Return[False, Module]];
  delta = Map[Cancel[Together[#]] &,
    Normal[matrix] - IdentityMatrix[dimension], {2}];
  AllTrue[Flatten[delta], # === 0 &]];

sccNoFamilyCollisionQ[cs_Association] := Module[
  {families = Lookup[cs, "Families", None], collisions, depths},
  If[!ListQ[families], Return[False, Module]];
  collisions = Flatten[Lookup[families, "Collisions", {}], 1];
  depths = Lookup[families, "CollisionDepth", {}];
  collisions === {} && AllTrue[depths, # === 0 &]];

sccNativeBlockCapabilities[cs_Association] := Module[
  {dimension = Lookup[cs, "SystemSize", 0]},
  <|"regular" -> TrueQ[Lookup[
      Lookup[cs, "IndicialData", <||>], "Regular", False]],
    "identity_gauge" ->
      (sccExactIdentityMatrixQ[Lookup[cs, "Gauge", None], dimension] &&
       sccExactIdentityMatrixQ[
         Lookup[cs, "GaugeInverse", None], dimension]),
    "identity_v" ->
      (sccExactIdentityMatrixQ[Lookup[cs, "V", None], dimension] &&
       sccExactIdentityMatrixQ[Lookup[cs, "VInv", None], dimension]),
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

$cppStaticOperatorCache = <||>;
$cppStaticOperatorCacheMax = 1024;

(* Decimal encoding of every prepared lag dominated the Wolfram side of
   repeated persistent solves if rebuilt per homogeneous column/SCC source.
   Cache the complete schema-1-compatible static payload once per exact
   prepared frame.  The hash is only an index and every hit compares the full
   structural signature; the opaque token then gives CppBackend a lightweight
   identity for the already collision-certified payload. *)
cppStaticOperatorPayload[cs_, prep_, blocks_List, fb_Integer, W_Integer,
    vPrep_, inputDigits_Integer, precisionBits_Integer] := Module[
  {d = cs["SystemSize"], signature, key, cached, dLags, denominators,
   nLags, assembly = Null, assemblyGroups, assemblyBase, payload, record},
  signature = {$cppSerializationDomain, $cppSerializationSymbols,
    inputDigits, precisionBits, d, fb, W, prep, blocks,
    If[AssociationQ[vPrep], vPrep, Automatic], cfg["ChopPrecision"]};
  key = Hash[signature, "SHA256"];
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
  If[AssociationQ[vPrep],
    assemblyGroups = vPrep["RationalGroups"];
    assemblyBase = Length[denominators];
    denominators = Join[denominators,
      Map[(cppScalar[#, inputDigits, cs] & /@
          #["DenominatorCoefficients"]) &, assemblyGroups]];
    assembly = <|
      "identity" -> TrueQ[Lookup[vPrep, "Identity", False]],
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
  payload = <|"domain" -> $cppSerializationDomain,
    "symbols" -> (SymbolName /@ $cppSerializationSymbols),
    "precision_bits" -> precisionBits,
    "d" -> d, "fb" -> fb, "w" -> W,
    "d_lags" -> dLags, "denominators" -> denominators,
    "nhat_lags" -> nLags,
    "d0_inverse" -> If[prep["d0InvScalar"] === None, Null,
      cppScalar[prep["d0InvScalar"], inputDigits, cs]],
    "blocks" -> Map[(#["Cols"] - 1) &, blocks],
    "assembly" -> assembly, "chop_digits" -> cfg["ChopPrecision"]|>;
  (* Stable content identity prevents an evicted Wolfram serialization cache
     entry from duplicating an already-retained native chart.  The digest is
     only an index: both caches retain and compare the complete certificate. *)
  record = <|"Signature" -> signature, "Payload" -> payload,
    "Token" -> ("de2-operator-" <> IntegerString[key, 16, 64])|>;
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
  outputDigits = wp + 20;
  precisionBits = Ceiling[inputDigits*Log[2, 10]] + 32;
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
  precisionBits = Ceiling[inputDigits*Log[2, 10]] + 32;
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
    inputDigits_Integer, cs_Association] := <|
  (* epsilon_shift is deliberately signed.  scc.prepare retains it and the
     later execution work contract proves that work_min supplies its halo. *)
  "epsilon_shift" -> prepared["EpsilonShift"],
  "center_pole_order" -> prepared["CenterPoleOrder"],
  "kernels" -> Map[cppScalar[#, inputDigits, cs] &,
    prepared["TaylorKernels"], {2}],
  "exact_identity" -> prepared["ExactIdentity"],
  "proven_zero" -> TrueQ[prepared["ProvenZero"]]|>;

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
          preparationShape, coefficient, t];
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
          #["Prepared"], inputDigits, cs]|> &, rawEntries]];
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
  $chartClearedCache = <||>; $exactSCCStructureCache = <||>;
  $cppStaticOperatorCache = <||>; $nativeSCCCompositeCache = <||>;
  $homogeneousFramePlanOverride = None;
  $cppHomogeneousFrameOverride = None;
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
  DiffExp2`CppBackend`ClearPersistentSessions[];);

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
   transformDepth, wideTop, terminalFb, startFb, adaptiveQ},
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"];
  reqMax = req["EpsWindow", "CompleteMax"];
  pMax = Max[0, Max[Table[
      logCeiling[cs, b["a"], b["b"], b["q"] - 1], {b, blocks}]]];
  cdMax = Max[0, Max[#["CollisionDepth"] & /@ fams]];
  symbolic = clearedSymbolic[cs];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  singleUseDepth = recurrenceSingleUsePoleDepth[symbolic];
  spectralDepth = spectralTransformPoleDepth[cs];
  transformDepth = finalTransformPoleDepth[cs, nmax];
  wideTop = reqMax + pMax + cdMax + 2 -
    Min[0, reqMin - pMax - 2];
  wideTop = Max[wideTop,
    reqMax + pMax + cdMax + poleDepth] + transformDepth;
  terminalFb = Min[Min[reqMin, 0] - pMax - cdMax - 2,
      reqMin - pMax - cdMax - poleDepth] - spectralDepth;
  startFb = Min[Min[reqMin, 0] - pMax - cdMax - 2,
      reqMin - pMax - cdMax - singleUseDepth] - spectralDepth;
  adaptiveQ = !TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    startFb > terminalFb && !TrueQ[$disableAdaptiveLowerFrames];
  If[!adaptiveQ, startFb = terminalFb];
  <|"Identity" -> {cs, req}, "PMax" -> pMax,
    "CollisionDepthMax" -> cdMax, "Symbolic" -> symbolic,
    "PoleDepth" -> poleDepth, "SingleUseDepth" -> singleUseDepth,
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

(* First production seam for a session-owned native LocalSolution.  It is
   intentionally narrower than SolveHomogeneous: no SCC orchestration,
   pseudo-resonant compensation, or rank-reduction gauge is hidden behind
   the opaque handle.  Those operations still require Wolfram tensors until
   their native counterparts preserve the same sequential completeness
   contract. *)
SolveNativeLocalFamily[cs_Association, req_Association,
    tag_Association, init_List] := Module[
  {d = Lookup[cs, "SystemSize", 0], a, b, p, nmax, reqMin, reqMax,
   matchingFamilies, matchingBlocks, allowedP, blocks, fams, pMax,
   pBudget, cdMax, symbolic,
   poleDepth, spectralDepth, transformDepth, wideTop, fb, W, prep, vPrep,
   symbols, domain, wp, inputDigits, precisionBits, staticRecord, request,
   persistentMetadata, checkpointIdentity, localMetadata, response,
   backendID},
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
  If[Length[init] =!= p + 1 || !AllTrue[init,
      ListQ[#] && Length[#] === d &&
        AllTrue[#, DiffExp2`EpsSeries`ESQ] &],
    err["E8", cs, <|"Tag" -> tag, "InitialDimensions" ->
      Quiet[Check[Dimensions[init], Missing["Ragged"]]],
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
  symbolic = clearedSymbolic[cs];
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
  prep = prepareCleared[cs, fb, W, symbolic];
  vPrep = prepareFramedMatrix[cs["V"],
    DiffExp2`Config`CanonicalEps[], fb, W, cs];
  symbols = cppRegulatorSymbols[cs, prep, a, b, p, nmax,
    None, init, vPrep];
  If[symbols =!= {},
    err["E5", cs, <|"Tag" -> tag,
      "RegulatorSymbols" -> (SymbolName /@ symbols), "Detail" ->
      "native local handle rejects unresolved analytic regulators; specialize them before solving"|>]];

  domain = If[TrueQ[$cppExactDomain], "rational", "acb"];
  wp = cfg["WorkingPrecision"];
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  precisionBits = Ceiling[inputDigits*Log[2, 10]] + 32;
  Block[{$cppSerializationDomain = domain,
      $cppSerializationSymbols = {}},
    staticRecord = cppStaticOperatorPayload[cs, prep, blocks, fb, W,
      vPrep, inputDigits, precisionBits];
    request = Block[{$cppStaticRecordOverride = staticRecord,
        $cppBuildRequestOnly = True},
      cppRunRecursionCore[cs, prep, a, b, p, nmax, None, fb, W,
        init, vPrep]];
    persistentMetadata = Append[cppPersistentMetadata[cs, fb, W],
      "PreparedToken" -> staticRecord["Token"]];
    checkpointIdentity = "de2-native-local-" <>
      IntegerString[Hash[{persistentMetadata["SystemIdentity"],
        persistentMetadata["ChartIdentity"], cs["ChartMap"],
        cs["Radius"], cs["Prescriptions"], {a, b, p}, req, init,
        staticRecord["Token"]}, "SHA256"],
        16, 64];
    localMetadata = cppNativeLocalMetadata[cs, a, b, p, inputDigits,
      checkpointIdentity];
    response = DiffExp2`CppBackend`RunPersistentLocalSolve[
      request, persistentMetadata, localMetadata]];
  If[FailureQ[response],
    err["E5", cs, <|"BackendFailure" -> response, "Detail" ->
      "persistent native local solve failed"|>]];
  If[Lookup[response, "status", "error"] =!= "ok",
    backendID = Lookup[response, "id", "E5"];
    err[If[MemberQ[{"E4", "E5", "E6"}, backendID], backendID, "E5"],
      cs, <|"BackendID" -> backendID,
        "Detail" -> Lookup[response, "detail",
          "persistent native local solve returned an error"]|>]];
  If[Lookup[response, "epsilon_max", reqMin - 1] < reqMax,
    Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[response]];
    err["E6", cs, <|"Tag" -> tag, "RequestedCompleteMax" -> reqMax,
      "AvailableCompleteMax" -> Lookup[response, "epsilon_max",
        Missing["NotAvailable"]],
      "Detail" ->
        "native local family work budget did not reach the requested epsilon order"|>]];
  <|"Type" -> "DiffExp2NativeLocalFamily",
    "Session" -> response["session"], "Local" -> response["local"],
    "NativeChart" -> response["chart"],
    "Tag" -> <|"a" -> a, "b" -> b, "p" -> p|>,
    "Chart" -> <|"Center" -> cs["Center"],
      "ChartMap" -> cs["ChartMap"], "Radius" -> cs["Radius"],
      "Prescriptions" -> cs["Prescriptions"]|>,
    "EpsWindow" -> <|"Min" -> response["epsilon_min"],
      "CompleteMax" -> response["epsilon_max"]|>,
    "TWindow" -> <|"CompleteMax" -> response["taylor_complete_max"]|>,
    "CheckpointIdentity" -> checkpointIdentity,
    "NativeSummary" -> KeyDrop[response,
      {"status", "session", "local", "chart", "metadata"}]|>];

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
        forcedFrame["FrameBase"] > terminalFb ||
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
    AssociateTo[prepCache, ff -> prepareCleared[
      cs, ff, wideTop - ff + 1, symbolic]]; prepCache[ff]];
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
  Join[sub, <|"SolveCacheTag" -> Lookup[cs, "SystemHash", None],
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
  Join[full, <|"SolveCacheTag" -> Lookup[cs, "SystemHash", None]|>]];

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

sccParticularBlockCost[target_Association, nmax_Integer:0] :=
  spectralTransformPoleDepth[target] +
    inverseSpectralTransformPoleDepth[target] +
    matrixEpsPoleDepth[target["Gauge"]] +
    matrixEpsPoleDepth[target["GaugeInverse"]] +
    recurrencePoleDepth[clearedSymbolic[target], nmax] +
    2 recurrenceSingleUsePoleDepth[clearedSymbolic[target]] +
    (* The canonical particular/log ladder can consume one honest source
       order per target Jordan position even for a regular zero diagonal. *)
    target["SystemSize"] +
    Max[0, Sequence @@ (# ["CollisionDepth"] & /@ target["Families"])];

sccInitialWorkHalo[cs_Association, blockSystems_List, nmax_Integer:0] := Module[
  {seq = cs["IntegrationSequence"], nb, edgeCost, incoming, path},
  nb = Length[seq["Components"]];
  edgeCost[{u_, v_}] := Module[{mat, target = blockSystems[[v]]},
    mat = cs["ThetaOriginal"][[seq["Components"][[v]],
      seq["Components"][[u]]]];
    matrixEpsPoleDepth[mat] + sccParticularBlockCost[target, nmax]];
  incoming = Table[Cases[seq["CondensationEdges"], {u_, v} :> u],
    {v, nb}];
  path = ConstantArray[0, nb];
  Do[If[incoming[[v]] =!= {},
    path[[v]] = Max[(path[[#]] + edgeCost[{#, v}]) & /@ incoming[[v]]]],
    {v, nb}];
  If[path === {}, 0, Max[path]]];

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
$nativeSCCColumnRunKeys = {"nmax", "p", "has_initial",
  "adaptive_probe", "a_target", "b_target", "a_shift_min",
  "a_shifts", "schedule", "initial", "initial_validity", "source",
  "return_u"};

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
  TrueQ[$cppExactDomain], TrueQ[$cppUsePersistentSessions],
  TrueQ[$numericizeAllPreparedNumbers],
  TrueQ[$disableGlobalClearedHoist],
  TrueQ[$disableIdentityNhatShortcut],
  TrueQ[$disableAdaptiveLowerFrames],
  TrueQ[$adaptiveLowerFrameProbe],
  TrueQ[$disableRationalDenominatorFusion],
  TrueQ[$disableGroupedSpectralTransform],
  TrueQ[$disablePolynomialNhatTransform],
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
  {requests, frameRecords, fieldRecords, sessionRecords, domain, names,
   symbols, serialization},
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
    "Domain" -> domain, "SymbolNames" -> names,
    "Serialization" -> KeyTake[serialization, {"domain", "symbols"}]|>];

sccNativeSourceShape[cs_Association, dimension_Integer,
    fb_Integer, completeMax_Integer, tOrder_Integer] := <|
  "Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
  "Radius" -> cs["Radius"], "Prescriptions" -> cs["Prescriptions"],
  "EpsWindow" -> <|"Min" -> fb, "CompleteMax" -> completeMax|>,
  "TWindow" -> <|"CompleteMax" -> tOrder|>,
  "Dimension" -> dimension|>;

sccCapturedBlockRecord[parentSystemRecord_List, parentGeometry_Association,
    seq_Association, blockcs_Association, captured_Association,
    capabilities_Association, block_Integer] := Module[
  {vertices = seq["Components"][[block]], expectedPrincipal,
   analytic, principal, capturedCapabilities, capturedGeometry,
   capturedIdentity},
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
  <|"block" -> block - 1, "vertices" -> (vertices - 1),
    "regular" -> capabilities["regular"],
    "identity_gauge" -> capabilities["identity_gauge"],
    "identity_v" -> capabilities["identity_v"],
    "no_pseudo" -> capabilities["no_pseudo"]|>];

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

PrepareNativeSCCComposite[cs_Association, req_Association] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], epsWindow,
   requestedMin, requestedMax, publicTOrder, workTOrder, plannedTop,
   workReq, signature, cacheKey, cached, cacheStats, blockSystems,
   capabilities, badBlocks, captures, capturedContract, fb, width, workTop,
   parentRecords, parentGeometry, parent, blockRecords, serialization,
   couplings, identity, manifest, prepared, result, scale, radius,
   center, missingReq, components, condensation, executionDescriptor,
   inputDigits, runRecords, taskRecords, columnPlans, blockDimensions,
   framePlans, forcedFrame},
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
  If[!TrueQ[Lookup[cs, "SCCSkeleton", False]] ||
      !AssociationQ[seq] ||
      Length[Lookup[seq, "Components", {}]] <= 1,
    err["E6", cs, <|"Detail" ->
      "native SCC preparation requires a multi-block SCC skeleton envelope"|>]];
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
      !(IntegerQ[center] || Head[center] === Rational) ||
      !((IntegerQ[scale] || Head[scale] === Rational) && scale =!= 0) ||
      radius === Infinity ||
      !((IntegerQ[radius] || Head[radius] === Rational) && radius > 0),
    err["E6", cs, <|"Scale" -> scale, "Radius" -> radius,
      "Detail" -> "native SCC first slice requires a rational center, nonzero rational scale, and positive finite rational radius"|>]];

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
  If[Length[$nativeSCCCompositeCache] >= $nativeSCCCompositeCacheMax,
    err["E6", cs, <|"Capacity" -> $nativeSCCCompositeCacheMax,
      "Detail" -> "native SCC composite cache capacity is exhausted; clear solver caches before preparing another public handle"|>]];

  components = seq["Components"];
  condensation = seq["CondensationEdges"];
  blockSystems = sccBlockChartSystem[cs, #] & /@
    Range[Length[components]];
  capabilities = sccNativeBlockCapabilities /@ blockSystems;
  badBlocks = Select[Range[Length[blockSystems]],
    Function[block, Module[{record = capabilities[[block]]},
      (* "regular" is a classification, not an admission predicate.  The
         retained C++ chart proves the complete exact affine-Jordan
         indicial operator for a singular block and revalidates every T/P/R
         schedule at execution.  Wolfram must therefore admit both Boolean
         classes here while continuing to require the exact producer facts
         which the current composite representation actually depends on. *)
      !MemberQ[{True, False}, Lookup[record, "regular", None]] ||
        !TrueQ[Lookup[record, "identity_gauge", False]] ||
        !TrueQ[Lookup[record, "identity_v", False]] ||
        !MemberQ[{True, False},
          Lookup[record, "no_pseudo", None]]]]];
  If[badBlocks =!= {},
    err["E6", cs, <|"UnsupportedBlocks" -> Map[
        <|"Block" -> #, "Capabilities" -> capabilities[[#]]|> &,
        badBlocks],
      "Detail" -> "native SCC preparation requires regular or exact affine-Jordan diagonal blocks with exact identity Gauge/GaugeInverse and V/VInv; no_pseudo is retained provenance, not an admission decision"|>]];
  requestedMin = epsWindow["Min"];
  requestedMax = epsWindow["CompleteMax"];
  publicTOrder = req["TOrder"];
  workTOrder = sccWorkTOrder[cs, req];
  plannedTop = requestedMax + sccInitialWorkHalo[
    cs, blockSystems, workTOrder];
  workReq = Join[req, <|
    "EpsWindow" -> Join[epsWindow, <|"CompleteMax" -> plannedTop|>],
    "TOrder" -> workTOrder|>];
  framePlans = homogeneousFramePlan[#, workReq] & /@ blockSystems;
  fb = Min[Lookup[framePlans, "TerminalFrameBase"]];
  workTop = Max[Lookup[framePlans, "FrameTop"]];
  forcedFrame = <|"FrameBase" -> fb, "FrameTop" -> workTop|>;
  captures = Block[{$shCache = <||>, $shSysTag = None,
      $cppStaticOperatorCache = <||>},
    MapThread[sccCaptureHomogeneousGroup[
        #1, workReq, #2, forcedFrame] &,
      {blockSystems, framePlans}]];
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
        "requested_min" -> requestedMin,
        "requested_max" -> requestedMax,
        "work_complete_max" -> workTop,
        "public_t_order" -> publicTOrder,
        "wolfram_coupling_depth" -> seq["CouplingDepth"]|>|>];
  blockRecords = MapThread[
    sccCapturedBlockRecord[parentRecords["exact_system_record"],
      parentGeometry, seq, #1, #2, #3, #4] &,
    {blockSystems, captures, capabilities,
      Range[Length[blockSystems]]}];
  serialization = capturedContract["Serialization"];
  couplings = Map[
    PrepareSCCCouplingMatrix[cs, #[[1]], #[[2]],
      sccNativeSourceShape[cs, Length[components[[#[[1]]]]],
        fb, workTop, workTOrder], serialization] &,
    condensation];
  identity = sccNativeCompositeIdentity[parent, blockRecords, couplings,
    capturedContract["Domain"], capturedContract["SymbolNames"]];
  If[!StringQ[identity],
    err["E6", cs, <|"Detail" ->
      "native SCC exact parent identity could not be serialized"|>]];
  manifest = <|"identity" -> identity, "parent" -> parent,
    "blocks" -> blockRecords, "couplings" -> couplings|>;
  inputDigits = DiffExp2`Tolerances`$InputPrecisionFactor*
    cfg["WorkingPrecision"];
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
  columnPlans = If[capturedContract["Domain"] === "rational" &&
      capturedContract["SymbolNames"] === {},
    sccNativeCapturedColumnPlans[cs, blockSystems,
      captures, capturedContract, inputDigits],
    (* Acb/symbolic composites remain valid prepared objects.  Exact SCC
       column execution is deliberately Rational-only, so never infer its
       unit lattice or T/P/R plan from numerical/string specializations. *)
    ConstantArray[{}, Length[blockSystems]]];
  prepared = DiffExp2`CppBackend`PreparePersistentSCC[captures, manifest];
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
    "ColumnPlans" -> columnPlans, "Contract" -> capturedContract,
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
    dimension_Integer, zero_, one_, integers_List] := Module[
  {fb = Lookup[contract, "FrameBase", None],
   width = Lookup[contract, "FrameWidth", None],
   nmax = Lookup[contract, "NMax", None], frameTop, unitIndex,
   initial, validity, schedule, zeroFrame, unitFrame, frames, selected},
  If[dimension < 1 ||
      !MatchQ[{fb, width, nmax}, {_Integer, _Integer, _Integer}] ||
      width < 1 || nmax < 0 || Length[integers] =!= nmax + 1 ||
      Sort[Keys[request]] =!= Sort[$nativeSCCColumnRunKeys] ||
      Sort[Keys[task]] =!= Sort[{"a", "b", "P"}] ||
      !TrueQ[zeroCanQ[task["a"]]] ||
      !TrueQ[zeroCanQ[task["b"]]] || task["P"] =!= 0,
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
      Lookup[request, "a_target", None] =!= zero ||
      Lookup[request, "b_target", None] =!= zero ||
      Lookup[request, "a_shift_min", None] =!= 0 ||
      Lookup[request, "a_shifts", None] =!= integers ||
      !ListQ[schedule] || Length[schedule] =!= nmax + 1 ||
      !(And @@ MapIndexed[Function[{row, index}, Module[
          {n = First[index] - 1},
          ListQ[row] && Length[row] === dimension &&
            AllTrue[row, AssociationQ[#] &&
              Sort[Keys[#]] === Sort[{"case", "da", "db"}] &&
              Lookup[#, "case", None] === If[n === 0, "R", "T"] &&
              Lookup[#, "da", None] === integers[[n + 1]] &&
              Lookup[#, "db", None] === zero &]]], schedule]) ||
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
    zero_String, one_String] := Module[
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
      TrueQ[Lookup[request, "adaptive_probe", True]] ||
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
    inputDigits_Integer] := Module[
  {dimension = blockcs["SystemSize"], blocks = blockList[blockcs],
   nmax = contract["NMax"], width = contract["FrameWidth"], p,
   encode, zero, schedule},
  If[sourceP < 0,
    err["E6", blockcs, <|"SourceLogPower" -> sourceP,
      "Detail" -> "native SCC source log power must be nonnegative"|>]];
  p = logCeiling[blockcs, a, b, sourceP, True];
  Block[{$cppSerializationDomain = "rational",
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
      "adaptive_probe" -> False,
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
   edges},
  seq = cs["IntegrationSequence"];
  topological = seq["TopologicalOrder"];
  edges = seq["CondensationEdges"];
  dimensions = Lookup[blockSystems, "SystemSize", None];
  requests = Lookup[captures, "Requests", None];
  tasks = Lookup[captures, "TaskMetadata", None];
  Block[{$cppSerializationDomain = "rational",
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
              inputDigits];
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
    Lookup[evidence, "identity_v", None] === "native-retained-assembly" &&
    Lookup[evidence, "regular", None] ===
      "collision-bound-producer-certificate" &&
    Lookup[evidence, "identity_gauge", None] ===
      "collision-bound-producer-certificate" &&
    Lookup[evidence, "no_pseudo", None] ===
      "collision-bound-producer-certificate"];

sccNativeSingularBlockStatisticsQ[stats_, handle_Association,
    blockDimensions_List, scalarShape_] := Module[
  {blockCharts = Lookup[stats, "block_charts", None],
   evidence = Lookup[stats, "capability_evidence", None],
   expected = If[TrueQ[scalarShape],
     "exact-rational-regular-singular-scalar-block-dag-column-v1",
     "exact-rational-regular-singular-jordan-block-dag-column-v2"]},
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
    Lookup[evidence, "identity_v", None] === "native-retained-assembly" &&
    Lookup[evidence, "regular_or_regular_singular", None] ===
      "collision-bound-producer-certificate" &&
    Lookup[evidence, "identity_gauge", None] ===
      "collision-bound-producer-certificate" &&
    Lookup[evidence, "jordan_indicial", None] ===
      "retained-exact-rational-full-matrix-certificate" &&
    Lookup[evidence, "no_pseudo", None] ===
      "producer-provenance-only-execution-revalidated-by-exact-schedule-certificate" &&
    Lookup[evidence, "pseudo_schedule_execution", None] ===
      "exact-rational-joint-compensation-and-formal-overlap-certificate" &&
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
          compensations > 0 && depth > 0 &&
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
   provenanceLocalComponent, encodedZero, encodedOne, encodedIntegers},
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
  singularExecution = AssociationQ[stats] &&
    MemberQ[{
      "exact-rational-regular-singular-scalar-block-dag-column-v1",
      "exact-rational-regular-singular-jordan-block-dag-column-v2"},
      Lookup[stats, "execution_scope", None]];
  scalarExecution = AllTrue[blockDimensions, # === 1 &];
  If[Length[blockDimensions] < 2 ||
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
      "Detail" -> "native SCC column requires two or more exact-rational or regular-Acb blocks with dimensions matching the parent SCC partition and no regulator field"|>]];
  If[singularExecution && domain =!= "rational",
    err["E6", cs, <|"Contract" -> contract,
      "Detail" -> "native singular SCC execution remains exact-rational because pseudo/resonance decisions require exact affine data"|>]];
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
    expectedCapability = If[scalarExecution,
      "exact-rational-regular-singular-scalar-block-dag-column-v1",
      "exact-rational-regular-singular-jordan-block-dag-column-v2"];
    expectedProvenanceSchema = If[scalarExecution,
      "diffexp2-native-scc-regular-singular-scalar-column-v1",
      "diffexp2-native-scc-regular-singular-jordan-column-v2"],

    If[domain === "acb",
      Block[{$cppSerializationDomain = "acb",
          $cppSerializationSymbols = {}},
        encodedZero = cppScalar[0, inputDigits, cs];
        encodedOne = cppScalar[1, inputDigits, cs];
        encodedIntegers = Table[cppScalar[n, inputDigits, cs],
          {n, 0, contract["NMax"]}]];
      componentMaps = MapThread[Function[{runs, tasks, dimension},
          MapThread[sccNativeCanonicalEncodedRegularBlockComponent[
              #1, #2, contract, dimension, encodedZero, encodedOne,
              encodedIntegers] &, {runs, tasks}]],
        {runRecords, taskRecords, blockDimensions}];
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
          stats, prepared, blockDimensions, scalarExecution],
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
   seedBlock = spec["SeedBlock"], backendID, provenance,
   forbiddenPayloadKeys, result},
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
        req["EpsWindow", "CompleteMax"] ||
      Lookup[response, "epsilon_max", None] =!=
        req["EpsWindow", "CompleteMax"] ||
      Lookup[response, "taylor_complete_max", None] =!= req["TOrder"] ||
      Lookup[response, "checkpoint_identity", None] =!= checkpointIdentity ||
      !TrueQ[Lookup[response, "native_retained", False]] ||
      Lookup[response, "json_coefficients", None] =!= 0 ||
      Lookup[response, "execution_capability", None] =!=
        expectedCapability ||
      Lookup[response, "pseudo_hit_count", None] =!= 0 ||
      (singularExecution && !sccNativePseudoDiagnosticsQ[response]) ||
      forbiddenPayloadKeys =!= {} || !AssociationQ[provenance] ||
      Lookup[provenance, "scc", None] =!= prepared["SCC"] ||
      Lookup[provenance, "scc_exact_identity", None] =!= parentIdentity ||
      Lookup[provenance, "seed_block", None] =!= seedBlock - 1 ||
      Lookup[provenance, "basis_index", None] =!= expectedBasisIndex ||
      !StringQ[Lookup[provenance, "exact_column_identity", None]] ||
      StringLength[Lookup[provenance, "exact_column_identity", ""]] === 0 ||
      !sccNativeBlockProvenanceIdentityQ[
        Lookup[provenance, "exact_column_identity", None],
        expectedProvenanceSchema, parentIdentity, expectedBasisIndex,
        provenanceLocalComponent, seed, targets],
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

SolveNativeSCCBasisColumn[cs_Association, req_Association,
    seedBlock_Integer, seedLocalComponent_Integer:1] := Module[
  {spec, response},
  spec = sccNativeBuildColumnRequest[
    cs, req, seedBlock, seedLocalComponent];
  response = DiffExp2`CppBackend`RunPersistentSCCColumn[
    spec["Prepared"], spec["Seed"], spec["Targets"],
    spec["CheckpointIdentity"]];
  sccNativeFinalizeColumn[cs, req, spec, response]];

SolveNativeSCCBasis[cs_Association, req_Association,
    threads_:Automatic] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], seeds, specs,
   prepared, workerCount, requestColumns, batch, raw, columns, cleanup,
   forbiddenPayloadKeys},
  If[!AssociationQ[seq] || Length[Lookup[seq, "Components", {}]] < 2,
    err["E6", cs, <|"Detail" ->
      "native SCC basis batch requires a multi-block SCC skeleton"|>]];
  If[DownValues[DiffExp2`CppBackend`RunPersistentSCCColumns] === {},
    err["E5", cs, <|"Detail" ->
      "CppBackend persistent SCC column-batch bridge is not available"|>]];
  seeds = Flatten[Table[{block, component},
      {block, Length[seq["Components"]]},
      {component, Length[seq["Components"][[block]]]}], 1];
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
    "EpsWindow" -> req["EpsWindow"],
    "TWindow" -> <|"CompleteMax" -> req["TOrder"]|>,
    "NativeSummary" -> KeyDrop[batch, {"results"}]|>];

SolveNativeRegularBasis[cs_Association, req_Association,
    threads_:Automatic] := Module[
  {seq = Lookup[cs, "IntegrationSequence", None], d = cs["SystemSize"],
   epsWindow = Lookup[req, "EpsWindow", None], epsMin, epsMax,
   physical, columns = {}, unitValues, built, cleanup, sessions, charts},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "SolveNativeRegularBasis requires a regular chart"|>]];
  If[threads =!= Automatic && !(IntegerQ[threads] && threads > 0),
    err["E6", cs, <|"Threads" -> threads,
      "Detail" -> "native regular basis thread count must be a positive integer or Automatic"|>]];
  If[AssociationQ[seq] && Length[Lookup[seq, "Components", {}]] > 1,
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
  cleanup[] := Scan[
    Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &, columns];
  Catch[
    Do[
      unitValues = Table[DiffExp2`EpsSeries`ESNew[epsMin,
          Table[If[component === basisIndex && power === 0, 1, 0],
            {power, epsMin, epsMax}]], {component, d}];
      built = SolveNativeLocalFamily[physical, req,
        <|"a" -> 0, "b" -> 0, "p" -> 0|>, {unitValues}];
      AppendTo[columns, Join[built,
        <|"Type" -> "DiffExp2NativeRegularBasisColumn",
          "BasisIndex" -> basisIndex|>]],
      {basisIndex, d}],
    "DiffExp2Error", Function[{failure, tag},
      cleanup[]; Throw[failure, tag]]];
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
    "NativeSummary" -> <|
      "execution_capability" ->
        "retained-regular-monolithic-unit-basis-v1",
      "columns" -> d, "worker_threads" -> 1,
      "json_coefficients" -> 0|>|>];

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
      "SolveCacheTag" -> Lookup[cs, "SystemHash", None],
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
