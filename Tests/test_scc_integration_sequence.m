(* Focused permanent contracts for exact SCC/block-sequential local solving.
   This file intentionally exercises only public solver entry points, except
   for narrow diagnostic/cache seams whose behavior is part of the tests. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps;
rho = Global`rho; sigma = Global`sigma;
req[min_, max_, n_] := <|"EpsWindow" -> <|"Min" -> min,
  "CompleteMax" -> max|>, "TOrder" -> n|>;
chart[name_, radius_:1] := <|"ChartVar" -> t, "Center" -> 0,
  "Scale" -> 1, "Radius" -> radius, "LocalRadius" -> radius,
  "Name" -> name, "Prescriptions" -> {}|>;
skeletonChart[name_, radius_:1] :=
  Append[chart[name, radius], "UseSCCSkeleton" -> True];
monoChartSystem[sys_, ch_] := Module[{cs},
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys,
    KeyDrop[ch, "UseSCCSkeleton"]]];
  If[FailureQ[cs], cs, KeyDrop[cs, "IntegrationSequence"]]];
evalValue[ls_, point_] := DiffExp2`SectorSeries`EvaluateLocalSolution[
  ls, point, "UsePade" -> False,
  "ComputeTailEstimates" -> False]["Value"];
epsCoefficient[ls_, point_, k_] :=
  DiffExp2`EpsSeries`ESCoefficient[evalValue[ls, point], k];
maxNumericMagnitude[data_] := Max[0, Sequence @@ Abs[Flatten[N[data, 60]]]];
tagList[ls_] := {# ["a"], # ["b"], # ["p"]} & /@ ls["Sectors"];

catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 3, "RecurrenceBackend" -> "Wolfram",
  "Variables" -> {}}]];

(* ------------------------------------------------------------------ *)
(* Exact three-block chain: f1'=0, f2'=f1, f3'=f2.  The source-first
   fundamental matrix is exp(t L), and CouplingDepth counts vertices. *)
chainSystem = <|"Matrix" -> {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}},
  "Variable" -> x|>;
chainCS = catchDE2[DiffExp2`Solve`PrepareChart[
  chainSystem, skeletonChart["scc-three-chain", 2]]];
chainReq = req[0, 2, 4];
chainSolved = catchDE2[DiffExp2`Solve`SolveChart[chainCS, chainReq]];
assert["scc_exact_three_chain_graph",
  !FailureQ[chainCS] && TrueQ[chainCS["SCCSkeleton"]] &&
  chainCS["IntegrationSequence", "Components"] === {{1}, {2}, {3}} &&
  chainCS["IntegrationSequence", "CondensationEdges"] === {{1, 2}, {2, 3}} &&
  chainCS["IntegrationSequence", "CouplingDepth"] === 3];
assert["scc_exact_three_chain_basis",
  !FailureQ[chainSolved] && chainSolved["CouplingDepth"] === 3 &&
  chainSolved["Basis", "Diagnostics", "SCCSolved"] === True &&
  chainSolved["Basis", "Diagnostics", "SCCExecutionMode"] ===
    "BlockSequentialStrict" &&
  chainSolved["Basis", "Columns"][[1, "Sectors", 1, "Coeffs", 1,
      1 ;; 3]] === {{1, 0, 0}, {0, 1, 0}, {0, 0, 1/2}} &&
  chainSolved["Basis", "Columns"][[2, "Sectors", 1, "Coeffs", 1,
      1 ;; 2]] === {{0, 1, 0}, {0, 0, 1}} &&
  chainSolved["Basis", "Columns"][[3, "Sectors", 1, "Coeffs", 1,
      1]] === {0, 0, 1}];

(* A regular SCC envelope must support the value-vector fast path without
   materializing a global indicial/spectral frame. *)
chainValues = {
  DiffExp2`EpsSeries`ESNew[0, {1, 0, 0}],
  DiffExp2`EpsSeries`ESNew[0, {0, 0, 0}],
  DiffExp2`EpsSeries`ESNew[0, {0, 0, 0}]};
chainValueSolution = catchDE2[DiffExp2`Solve`SolveValueRegular[
  chainCS, chainReq, chainValues]];
assert["scc_regular_envelope_value_recursion",
  !FailureQ[chainValueSolution] && !KeyExistsQ[chainCS, "V"] &&
  maxNumericMagnitude[epsCoefficient[chainValueSolution, 1/3, 0] -
      {1, 1/3, 1/18}] < 10^-70];

(* ------------------------------------------------------------------ *)
(* An edge proportional only to epsilon is still an exact graph edge. *)
epsEdgeSystem = <|"Matrix" -> {{0, 0}, {eps, 0}}, "Variable" -> x|>;
epsEdgeCS = catchDE2[DiffExp2`Solve`PrepareChart[
  epsEdgeSystem, skeletonChart["scc-eps-edge"]]];
epsEdgeSolved = catchDE2[DiffExp2`Solve`SolveChart[
  epsEdgeCS, req[0, 2, 3]]];
assert["scc_epsilon_only_edge_is_structural",
  !FailureQ[epsEdgeSolved] &&
  epsEdgeCS["IntegrationSequence", "Components"] === {{1}, {2}} &&
  MemberQ[epsEdgeCS["IntegrationSequence", "DependencyEdges"], {1, 2}] &&
  maxNumericMagnitude[
    epsCoefficient[epsEdgeSolved["Basis", "Columns"][[1]], 1/4, 1] -
      {0, 1/4}] < 10^-70];

(* Block preparation at a nontrivial affine chart must use the original
   physical principal submatrix.  Its theta matrix is exactly the parent
   slice; no synthetic theta/t re-charting is allowed. *)
affineSystem = <|"Matrix" -> {
    {(1 + eps)/(1 - x/2), 0, 0},
    {(2 - eps)/(1 + x), (-1 + 2 eps)/(2 - x), 0},
    {(3 + eps)/(3 - x), (1 - eps)/(1 + x/2),
      (2 + eps)/(1 - x/3)}}, "Variable" -> x|>;
affineChart = <|"ChartVar" -> t, "Center" -> 1/3, "Scale" -> 2/5,
  "Radius" -> 3/2, "LocalRadius" -> 3/2,
  "Name" -> "scc-affine-block-principal", "Prescriptions" -> {},
  "UseSCCSkeleton" -> True|>;
affineCS = catchDE2[DiffExp2`Solve`PrepareChart[
  affineSystem, affineChart]];
affineBlocks = If[FailureQ[affineCS], {}, Table[
  catchDE2[DiffExp2`Solve`Private`sccBlockChartSystem[affineCS, block]],
  {block, Length[affineCS["IntegrationSequence", "Components"]]}]];
affineBlockParity = If[FailureQ[affineCS] ||
    AnyTrue[affineBlocks, FailureQ], False,
  And @@ Table[Module[{indices =
      affineCS["IntegrationSequence", "Components"][[block]], diff},
    diff = Map[Cancel[Together[#]] &,
      affineBlocks[[block, "ThetaOriginal"]] -
        affineCS["ThetaOriginal"][[indices, indices]], {2}];
    diff === ConstantArray[0, Dimensions[diff]] &&
      affineBlocks[[block, "ChartMap"]] === affineCS["ChartMap"]],
    {block, Length[affineBlocks]}]];
assert["scc_affine_block_uses_original_principal_system",
  TrueQ[affineBlockParity]];

(* ------------------------------------------------------------------ *)
(* A 1/eps off-diagonal coupling needs one future source order.  The SCC
   halo must reproduce the strict monolithic result at the requested top. *)
haloSystem = <|"Matrix" -> {{0, 0}, {1/eps, 0}}, "Variable" -> x|>;
haloChart = skeletonChart["scc-particular-halo"];
haloCS = catchDE2[DiffExp2`Solve`PrepareChart[haloSystem, haloChart]];
haloFullCS = monoChartSystem[haloSystem, haloChart];
haloReq = req[0, 1, 4];
haloSource = <|"Sectors" -> {<|"a" -> 1, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[
      {If[n === 0, {1, 2, 3}[[k + 1]], 0], 0},
      {k, 0, 2}, {n, 0, 4}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 4|>|>;
haloSCC = catchDE2[DiffExp2`Solve`SolveParticular[
  haloCS, haloSource, haloReq]];
haloFull = catchDE2[DiffExp2`Solve`SolveParticular[
  haloFullCS, haloSource, haloReq]];
assert["scc_inverse_epsilon_particular_halo_parity",
  !FailureQ[haloSCC] && !FailureQ[haloFull] &&
  haloSCC["EpsWindow"] === haloFull["EpsWindow"] ===
    <|"Min" -> -1, "CompleteMax" -> 1|> &&
  tagList[haloSCC] === tagList[haloFull] &&
  maxNumericMagnitude[
    epsCoefficient[haloSCC, 1/3, -1] -
      epsCoefficient[haloFull, 1/3, -1]] < 10^-70 &&
  maxNumericMagnitude[
    epsCoefficient[haloSCC, 1/3, 1] -
      epsCoefficient[haloFull, 1/3, 1]] < 10^-70];

(* A finite source is allowed to miss the requested upper window.  Its
   honest nonempty [-1,-1] result is diagnostic, not an E6. *)
haloLowSource = Join[haloSource, <|"Sectors" -> {<|"a" -> 1,
    "b" -> 0, "p" -> 0, "Coeffs" -> Table[
      {If[k === 0 && n === 0, 1, 0], 0},
      {k, 0, 0}, {n, 0, 4}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 0|>|>];
haloLowReq = req[-2, 1, 4];
haloLowSCC = catchDE2[DiffExp2`Solve`SolveParticular[
  haloCS, haloLowSource, haloLowReq]];
haloLowFull = catchDE2[DiffExp2`Solve`SolveParticular[
  haloFullCS, haloLowSource, haloLowReq]];
assert["scc_finite_source_honest_lower_window",
  !FailureQ[haloLowSCC] && !FailureQ[haloLowFull] &&
  haloLowSCC["EpsWindow"] === <|"Min" -> -1, "CompleteMax" -> -1|> &&
  haloLowFull["EpsWindow"] === haloLowSCC["EpsWindow"] &&
  AssociationQ[haloLowSCC["Diagnostics", "RequestedWindowNotReached"]]];

(* ------------------------------------------------------------------ *)
(* The target SCC is itself rank-reduced.  The source is physical theta
   data; SolveParticular must apply the target GaugeInverse exactly once. *)
rankSystem = <|"Matrix" -> {{0, 0, 0}, {1/x, 0, 1/x^2}, {0, 1, 0}},
  "Variable" -> x|>;
rankCS = catchDE2[DiffExp2`Solve`PrepareChart[
  rankSystem, skeletonChart["scc-rank-target"]]];
rankBlocks = If[FailureQ[rankCS], {}, Table[
  catchDE2[DiffExp2`Solve`Private`sccBlockChartSystem[rankCS, block]],
  {block, Length[rankCS["IntegrationSequence", "Components"]]}]];
rankSource = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0,
      {0, 1, 0}, {0, 0, 0}], {k, 0, 3}, {n, 0, 5}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 3|>,
  "TWindow" -> <|"CompleteMax" -> 5|>|>;
rankParticular = catchDE2[DiffExp2`Solve`SolveParticular[
  rankCS, rankSource, req[0, 2, 5]]];
assert["scc_rank_reduced_target_source_frame",
  !FailureQ[rankParticular] &&
  rankCS["IntegrationSequence", "Components"] === {{1}, {2, 3}} &&
  Length[rankBlocks] === 2 &&
  rankBlocks[[2, "Gauge"]] =!= IdentityMatrix[2] &&
  maxNumericMagnitude[epsCoefficient[rankParticular, 1/3, 0] -
      {0, -1, -1/3}] < 10^-60];

(* ------------------------------------------------------------------ *)
(* The residual probe must remain sensitive for negative leading powers. *)
negativeACS = catchDE2[DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{-3/x}}, "Variable" -> x|>,
  chart["scc-negative-a-residual"]]];
negativeABasis = catchDE2[DiffExp2`Solve`SolveHomogeneous[
  negativeACS, req[0, 0, 5]]];
negativeAColumn = If[FailureQ[negativeABasis], <||>,
  negativeABasis["Columns"][[1]]];
negativeAGood = If[AssociationQ[negativeAColumn],
  catchDE2[DiffExp2`Solve`ODEResidualCheck[negativeACS, negativeAColumn]],
  negativeABasis];
negativeABadColumn = If[AssociationQ[negativeAColumn],
  ReplacePart[negativeAColumn,
    {"Sectors", 1, "Coeffs", 1, 4, 1} -> 1], <||>];
negativeABad = If[AssociationQ[negativeABadColumn],
  catchDE2[DiffExp2`Solve`ODEResidualCheck[negativeACS,
    negativeABadColumn]], negativeABasis];
assert["scc_negative_a_residual_corruption_is_loud",
  !FailureQ[negativeAGood] && FailureQ[negativeABad] &&
  negativeABad["ID"] === "E7"];

(* ------------------------------------------------------------------ *)
(* Block scratch cache entries must not consume the persistent parent cap. *)
oldCacheCap = DiffExp2`Solve`Private`$shCacheMax;
DiffExp2`Solve`Private`$shCacheMax = 1;
DiffExp2`Solve`Private`$shCache = <||>;
DiffExp2`Solve`Private`$shSysTag = None;
cacheBasis1 = catchDE2[DiffExp2`Solve`SolveHomogeneous[chainCS, chainReq]];
parentCacheKey = DiffExp2`Solve`Private`homogeneousCacheKey[
  chainCS, chainReq];
cacheBasis2 = catchDE2[DiffExp2`Solve`SolveHomogeneous[chainCS, chainReq]];
assert["scc_parent_cache_survives_subblock_capacity",
  !FailureQ[cacheBasis1] && SameQ[cacheBasis1, cacheBasis2] &&
  Length[DiffExp2`Solve`Private`$shCache] === 1 &&
  KeyExistsQ[DiffExp2`Solve`Private`$shCache, parentCacheKey]];
DiffExp2`Solve`Private`$shCacheMax = oldCacheCap;
DiffExp2`Solve`Private`$shCache = <||>;
DiffExp2`Solve`Private`$shSysTag = None;

(* ------------------------------------------------------------------ *)
(* A downstream diagonal recurrence with negative epsilon valuation must
   select one strict monolithic recurrence.  Recombination records come from
   that actual full frame and must work even when passed the SCC envelope. *)
monoSystem = <|"Matrix" -> {{eps/x, 0},
    {1/x, 2 eps/x + 1/(eps (1 - x))}}, "Variable" -> x|>;
monoCS = catchDE2[DiffExp2`Solve`PrepareChart[
  monoSystem, skeletonChart["scc-monolithic-rule", 1/2]]];
monoReq = req[0, 2, 4];
monoSolved = catchDE2[DiffExp2`Solve`SolveChart[monoCS, monoReq]];
monoDiagnostics = If[FailureQ[monoSolved], <||>,
  monoSolved["Basis", "Diagnostics"]];
monoGroups = Lookup[monoDiagnostics, "SCCRecombineGroups", Missing[]];
monoRecombined = If[FailureQ[monoSolved], monoSolved,
  catchDE2[DiffExp2`Transport`Private`recombineDegenerate[
    monoCS, monoSolved["Basis", "Columns"], monoSolved["Basis", "Specs"],
    monoDiagnostics]]];
assert["scc_monolithic_execution_rule_and_plan",
  !FailureQ[monoSolved] && monoDiagnostics["SCCSolved"] === False &&
  monoDiagnostics["SCCExecutionMode"] === "MonolithicStrict" &&
  monoDiagnostics["SCCExecutionReason"] ===
    "DownstreamRecurrencePoleDepth" &&
  monoDiagnostics["SCCExecutionPlan", "Mode"] === "MonolithicStrict" &&
  AllTrue[monoDiagnostics["SCCExecutionPlan", "OffendingTargets"],
    # ["RecurrencePoleDepth"] > 0 &]];
assert["scc_monolithic_explicit_recombination_groups",
  ListQ[monoGroups] && monoGroups =!= {} &&
  AllTrue[monoGroups, AssociationQ[#] && Length[# ["Columns"]] >= 2 &&
      IntegerQ[# ["EpsZeroDegeneracy"]] &&
      # ["SCCBlock"] === "MonolithicStrict" &] &&
  ListQ[monoRecombined] &&
  monoRecombined =!= monoSolved["Basis", "Columns"]];

(* If only monolithic preparation is unsupported, coarsening is declined and
   the already-supported strict block path remains authoritative.  Use a
   candidate whose block path is independently supported; the recombination
   fixture above has a genuine block-local CASE-P obstruction and must stay
   loud if its monolithic preparation is unavailable. *)
declineSystem = <|"Matrix" -> {{0, 0}, {1, 1/(eps (1 - x))}},
  "Variable" -> x|>;
declineCS = catchDE2[DiffExp2`Solve`PrepareChart[
  declineSystem, skeletonChart["scc-monolithic-decline"]]];
declineReq = req[0, 1, 2];
DiffExp2`Solve`Private`$shCache = <||>;
DiffExp2`Solve`Private`$shSysTag = None;
declinedSolved = Block[{
    DiffExp2`Solve`Private`sccTryMonolithicChartSystem =
      Function[Failure["DiffExp2", <|"ID" -> "E3",
        "Detail" -> "injected monolithic preparation refusal"|>]]},
  catchDE2[DiffExp2`Solve`SolveChart[declineCS, declineReq]]];
assert["scc_monolithic_preparation_decline_uses_strict_blocks",
  !FailureQ[declinedSolved] &&
  declinedSolved["Basis", "Diagnostics", "SCCSolved"] === True &&
  declinedSolved["Basis", "Diagnostics", "SCCExecutionMode"] ===
    "BlockSequentialStrict" &&
  declinedSolved["Basis", "Diagnostics", "SCCExecutionReason"] ===
    "CoarseningDeclined" &&
  FailureQ[declinedSolved["Basis", "Diagnostics",
    "CoarseningDeclined"]] &&
  AssociationQ[declinedSolved["Basis", "Diagnostics",
    "SCCCoarseningCandidate"]]];

(* ------------------------------------------------------------------ *)
(* Mixed analytic-regulator blocks are intentionally solved sequentially;
   the optional C++ backend must retain their exact formal field. *)
If[TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  catchDE2[DiffExp2`Config`UpdateConfiguration[{
    "RecurrenceBackend" -> "Cpp", "Variables" -> {rho, sigma}}]];
  regulatorSystem = <|"Matrix" -> {{rho, 0}, {0, sigma}},
    "Variable" -> x|>;
  regulatorCS = catchDE2[DiffExp2`Solve`PrepareChart[
    regulatorSystem, skeletonChart["scc-cpp-mixed-regulators"]]];
  regulatorSolved = catchDE2[DiffExp2`Solve`SolveHomogeneous[
    regulatorCS, req[0, 0, 2]]];
  assert["scc_cpp_mixed_regulators_sequential_exact",
    !FailureQ[regulatorSolved] &&
    regulatorSolved["Diagnostics", "SCCExecutionMode"] ===
      "BlockSequentialStrict" &&
    regulatorSolved["Columns"][[1, "Sectors", 1, "Coeffs", 1,
        All, 1]] === {1, rho, rho^2/2} &&
    regulatorSolved["Columns"][[2, "Sectors", 1, "Coeffs", 1,
        All, 2]] === {1, sigma, sigma^2/2}];

  (* Disabling the grouped spectral transform also disables native block
     batching, but each strict C++ block solve must remain available. *)
  catchDE2[DiffExp2`Config`UpdateConfiguration[{
    "RecurrenceBackend" -> "Cpp", "Variables" -> {}}]];
  groupedOffSystem = <|"Matrix" -> {{1 + eps, 0}, {1 + x, 2 - eps}},
    "Variable" -> x|>;
  groupedOffCS = catchDE2[DiffExp2`Solve`PrepareChart[
    groupedOffSystem, skeletonChart["scc-cpp-grouped-off"]]];
  DiffExp2`Solve`Private`$shCache = <||>;
  DiffExp2`Solve`Private`$shSysTag = None;
  groupedOffSolved = Block[{
      DiffExp2`Solve`Private`$disableGroupedSpectralTransform = True},
    catchDE2[DiffExp2`Solve`SolveHomogeneous[
      groupedOffCS, req[0, 2, 5]]]];
  assert["scc_cpp_grouped_transform_off_sequential_fallback",
    !FailureQ[groupedOffSolved] &&
    Length[groupedOffSolved["Columns"]] === 2 &&
    groupedOffSolved["Diagnostics", "SCCExecutionMode"] ===
      "BlockSequentialStrict"];

  (* More blocks than the native prewarm capacity take the sequential C++
     route; this is a supported execution rule, not a capacity error. *)
  manyBlockDimension = 65;
  manyBlockSystem = <|"Matrix" -> ConstantArray[0,
      {manyBlockDimension, manyBlockDimension}], "Variable" -> x|>;
  manyBlockCS = catchDE2[DiffExp2`Solve`PrepareChart[
    manyBlockSystem, skeletonChart["scc-cpp-over-capacity"]]];
  DiffExp2`Solve`Private`$shCache = <||>;
  DiffExp2`Solve`Private`$shSysTag = None;
  manyBlockSolved = catchDE2[DiffExp2`Solve`SolveHomogeneous[
    manyBlockCS, req[0, 0, 1]]];
  assert["scc_cpp_more_than_capacity_uses_sequential_blocks",
    !FailureQ[manyBlockSolved] &&
    Length[manyBlockCS["IntegrationSequence", "Components"]] === 65 &&
    Length[manyBlockSolved["Columns"]] === 65 &&
    manyBlockSolved["Diagnostics", "SCCExecutionMode"] ===
      "BlockSequentialStrict"];

  (* The persistent cache owns assembled parent bases, not block scratch.
     Two prewarm waves must retain all four parent entries under cap eight. *)
  prewarmSystem = <|"Matrix" -> {{1 + x + eps, 0},
      {1 + x^2, -1 + 2 x + eps}}, "Variable" -> x|>;
  prewarmChart[center_, name_] := <|"ChartVar" -> t,
    "Center" -> center, "Scale" -> 1, "Radius" -> 1,
    "LocalRadius" -> 1, "Name" -> name, "Prescriptions" -> {},
    "UseSCCSkeleton" -> True|>;
  prewarmSystems = MapIndexed[catchDE2[DiffExp2`Solve`PrepareChart[
      prewarmSystem, prewarmChart[#1,
        "scc-prewarm-parent-" <> ToString[First[#2]]]]] &,
    {0, 1/5, 2/5, 3/5}];
  prewarmReq = req[0, 1, 4];
  prewarmOldCap = DiffExp2`Solve`Private`$shCacheMax;
  DiffExp2`Solve`Private`$shCacheMax = 8;
  DiffExp2`Solve`Private`$shCache = <||>;
  DiffExp2`Solve`Private`$shSysTag = None;
  prewarmFirst = catchDE2[DiffExp2`Solve`PrewarmHomogeneousBatch[
    prewarmSystems[[1 ;; 2]], prewarmReq]];
  prewarmSecond = catchDE2[DiffExp2`Solve`PrewarmHomogeneousBatch[
    prewarmSystems[[3 ;; 4]], prewarmReq]];
  prewarmParentKeys = DiffExp2`Solve`Private`homogeneousCacheKey[
      #, prewarmReq] & /@ prewarmSystems;
  assert["scc_cpp_prewarm_retains_parent_entries_across_waves",
    !FailureQ[prewarmFirst] && !FailureQ[prewarmSecond] &&
    Length[prewarmFirst] === 2 && Length[prewarmSecond] === 2 &&
    Length[DiffExp2`Solve`Private`$shCache] === 4 &&
    AllTrue[prewarmParentKeys,
      KeyExistsQ[DiffExp2`Solve`Private`$shCache, #] &]];
  DiffExp2`Solve`Private`$shCacheMax = prewarmOldCap,
  Print["  SKIP: scc_cpp_mixed_regulators_sequential_exact ",
    "(compiled backend unavailable)"];
  Scan[assert[#, True] &, {
    "scc_cpp_mixed_regulators_sequential_exact",
    "scc_cpp_grouped_transform_off_sequential_fallback",
    "scc_cpp_more_than_capacity_uses_sequential_blocks",
    "scc_cpp_prewarm_retains_parent_entries_across_waves"}]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
