(* DiffExp2/NativeTransport.m -- persistent native chart-atlas and arm
   orchestration.  Wolfram owns exact input/planning metadata; no retained
   local coefficient tensor is materialized here. *)

BeginPackage["DiffExp2`NativeTransport`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Solve`", "DiffExp2`Transport`",
   "DiffExp2`CppBackend`"}];

PrepareNativeRegularIndependentArms::usage =
  "PrepareNativeRegularIndependentArms[sys,boundary,lowerPlan,upperPlan] prepares one shared retained anchor, every regular receiving basis, and one exact lower/upper native tile plan. It returns only opaque native locals/bases and exact atlas metadata.";
RunNativeRegularIndependentArms::usage =
  "RunNativeRegularIndependentArms[atlas,cvec,var] marches both prepared regular arms using retained plan-derived matches, materializes each receiving local in C++, applies cvec(center+scale t,eps) natively, and integrates every planned tile. The current bridge executes arms sequentially; the retained atlas is compatible with the native concurrent-arm runner.";
ReleaseNativeRegularIndependentArms::usage =
  "ReleaseNativeRegularIndependentArms[runOrAtlas] releases public line, projected-local, materialized-local, match, basis, anchor, and tile-plan handles created by the explicit native regular-arm seam. Prepared chart/SCC caches remain session-owned and reusable.";
NativeRegularLineIntegral::usage =
  "NativeRegularLineIntegral[sys,boundary,from,{lo,hi},cvec] runs the explicit persistent-native regular-arm seam and returns an association whose Value is the compatibility EpsSeries integral. The anchor must lie strictly inside {lo,hi}; unsupported singular/nonrational geometry fails loudly without fallback.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "NativeTransport"|>, payload]];
cfg = DiffExp2`Config`CFG;
esQ = DiffExp2`EpsSeries`ESQ;
esNew = DiffExp2`EpsSeries`ESNew;
esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax;
esAdd = DiffExp2`EpsSeries`ESAdd;
esScale = DiffExp2`EpsSeries`ESScale;

exactRationalQ[value_] := IntegerQ[value] || Head[value] === Rational;

exactRationalString[value_, label_String] := Module[{canonical},
  canonical = Quiet[Check[RootReduce[value], value]];
  If[!exactRationalQ[canonical],
    err["E6", <|"Field" -> label, "Value" -> value,
      "Detail" -> "the current persistent exact-path protocol requires rational real geometry"|>]];
  ToString[canonical, InputForm]];

sameExactQ[left_, right_] := SameQ[left, right] ||
  TrueQ[Quiet[Check[PossibleZeroQ[RootReduce[left - right]], False]]];

inClosedArmQ[value_, from_, to_] := Module[{lo = Min[from, to], hi = Max[from, to]},
  TrueQ[lo <= value <= hi]];

nativeBoundaryValues[boundary_, dimension_Integer] := Module[{rows},
  Which[
    ListQ[boundary] && Length[boundary] === dimension &&
        AllTrue[boundary, esQ], boundary,
    ListQ[boundary] && Length[boundary] === dimension &&
        AllTrue[boundary, ListQ] &&
        SameQ @@ (Length /@ boundary) && Length[First[boundary]] > 0,
      rows = boundary;
      esNew[0, #] & /@ rows,
    True,
      err["E8", <|"Dimension" -> dimension,
        "BoundaryShape" -> Quiet[Check[Dimensions[boundary], "ragged"]],
        "Detail" -> "native regular transport requires d EpsSeries components or a rectangular d-by-epsilon coefficient table"|>]]];

normalizeSharedAnchor[lower_Association, upper_Association] := Module[
  {lowerAnchor, upperAnchor, commonRadius, commonMatchRadius, scale,
   common, lower2, upper2, roc = cfg["RadiusOfConvergence"]},
  If[lower["Direction"] =!= -1 || upper["Direction"] =!= 1 ||
      !sameExactQ[lower["From"], upper["From"]] ||
      !TrueQ[lower["To"] < lower["From"] < upper["To"]],
    err["E8", <|"Lower" -> KeyTake[lower, {"From", "To", "Direction"}],
      "Upper" -> KeyTake[upper, {"From", "To", "Direction"}],
      "Detail" -> "native independent arms require lower < shared anchor < upper"|>]];
  lowerAnchor = First[lower["Charts"]];
  upperAnchor = First[upper["Charts"]];
  If[!sameExactQ[lowerAnchor["Center"], upperAnchor["Center"]] ||
      !SameQ[Lookup[lowerAnchor, "Prescriptions", {}],
        Lookup[upperAnchor, "Prescriptions", {}]],
    err["E8", <|"LowerAnchor" -> lowerAnchor,
      "UpperAnchor" -> upperAnchor,
      "Detail" -> "independent arm anchors disagree in center or branch prescription"|>]];
  If[!AllTrue[{
        lowerAnchor["Radius"], upperAnchor["Radius"],
        lowerAnchor["MatchRadius"], upperAnchor["MatchRadius"], roc},
      exactRationalQ],
    err["E6", <|"LowerAnchor" -> lowerAnchor,
      "UpperAnchor" -> upperAnchor,
      "Detail" -> "shared-anchor normalization currently requires rational radius and scale data"|>]];
  (* SegmentLine caps anchor geometry by each arm length.  The larger of the
     two conservative caps is still no larger than the true singularity
     radius and gives one chart capable of serving both arms. *)
  commonRadius = Max[lowerAnchor["Radius"], upperAnchor["Radius"]];
  commonMatchRadius = Max[lowerAnchor["MatchRadius"],
    upperAnchor["MatchRadius"]];
  scale = Together[commonMatchRadius/roc];
  common = Join[lowerAnchor, <|"Radius" -> commonRadius,
    "MatchRadius" -> commonMatchRadius, "Scale" -> scale,
    "LocalRadius" -> Together[commonRadius/scale]|>];
  lower2 = Join[lower, <|"Charts" -> ReplacePart[lower["Charts"], 1 -> common]|>];
  upper2 = Join[upper, <|"Charts" -> ReplacePart[upper["Charts"], 1 -> common]|>];
  DiffExp2`Transport`ValidatePlan[lower2];
  DiffExp2`Transport`ValidatePlan[upper2];
  {lower2, upper2}];

nativeBranchSheets[plan_Association] := Module[
  {records, grouped, factorString},
  records = Flatten[Lookup[plan["Charts"], "Prescriptions", {}], 1];
  factorString[record_Association] := ToString[
    Lookup[record, "ExactFactor", Lookup[record, "Factor", None]],
    InputForm];
  grouped = GatherBy[records, factorString];
  Map[Module[{signs = DeleteDuplicates[Lookup[#, "Sign", None]], factor},
      factor = factorString[First[#]];
      If[Length[signs] =!= 1 || !MemberQ[{-1, 1}, First[signs]],
        err["E7", <|"Factor" -> factor, "Prescriptions" -> #,
          "Detail" -> "native path topology found conflicting branch-sheet signs"|>]];
      <|"factor_exact" -> factor, "sign" -> First[signs]|>] &,
    grouped]];

nativeComplexProjections[plan_Association] := Module[
  {singularities = plan["Singularities"], all, real, projected, from,
   to, nonreal, data, pairs, relevant},
  all = Lookup[singularities, "All", {}];
  real = Lookup[singularities, "Real", {}];
  projected = Lookup[singularities, "Projected", {}];
  from = plan["From"]; to = plan["To"];
  (* Exact pair membership avoids classifying a nearly-real algebraic root
     from a floating midpoint. *)
  nonreal = Select[all, Function[root,
    !AnyTrue[real, Function[r, sameExactQ[root, r]]]]];
  data = Map[Function[root, Module[{re, im, h},
      re = Quiet[RootReduce[(root + Conjugate[root])/2]];
      im = Quiet[RootReduce[(root - Conjugate[root])/(2 I)]];
      h = If[TrueQ[Re[N[im, 60]] < 0], -im, im];
      {re, Quiet[RootReduce[h]]}]], nonreal];
  pairs = DeleteDuplicatesBy[data, ToString[#, InputForm] &];
  relevant = Select[pairs, Function[pair,
    AnyTrue[{pair[[1]] - pair[[2]], pair[[1]], pair[[1]] + pair[[2]]},
      inClosedArmQ[#, from, to] &]]];
  Map[Function[pair, Module[{re = pair[[1]], h = pair[[2]], flags},
    If[!exactRationalQ[re] || !exactRationalQ[h] || !TrueQ[h > 0],
      err["E6", <|"Projection" -> pair,
        "Detail" -> "an on-arm complex projection is not representable by the current rational native path protocol"|>]];
    flags = {AnyTrue[projected, sameExactQ[#, re - h] &],
      AnyTrue[projected, sameExactQ[#, re] &],
      AnyTrue[projected, sameExactQ[#, re + h] &]};
    <|"source_identity" -> ("complex-projection:" <>
        IntegerString[Hash[{re, h}, "SHA256"], 16, 64]),
      "real_part_exact" -> exactRationalString[re, "projection real part"],
      "imaginary_magnitude_exact" ->
        exactRationalString[h, "projection imaginary magnitude"],
      "retain_minus_imaginary" -> First[flags],
      "retain_real_part" -> flags[[2]],
      "retain_plus_imaginary" -> Last[flags]|>]], relevant]];

nativeArmTopology[plan_Association] := Module[
  {singularities = plan["Singularities"], from = plan["From"],
   to = plan["To"], real, boundary},
  real = Select[Lookup[singularities, "Real", {}],
    inClosedArmQ[#, from, to] &];
  boundary = DeleteDuplicates@Select[
    Lookup[singularities, "ProjectionWaypoints", {}],
    inClosedArmQ[#, from, to] &];
  <|"singular_points" ->
      (exactRationalString[#, "real singular point"] & /@ real),
    "boundary_points" ->
      (exactRationalString[#, "projection waypoint"] & /@ boundary),
    "complex_projections" -> nativeComplexProjections[plan],
    "branch_sheets" -> nativeBranchSheets[plan]|>];

nativeBasisOwner[basis_Association] := Module[{owner},
  owner = Lookup[basis, "NativeSCC",
    Lookup[basis, "NativeChart", None]];
  If[!StringQ[owner] || StringLength[owner] == 0,
    err["E6", <|"Basis" -> KeyTake[basis,
        {"Type", "Session", "NativeSCC", "NativeChart"}],
      "Detail" -> "retained regular basis exposes no native chart/SCC owner"|>]];
  owner];

nativeArmRequest[plan_Association, owners_List] := Module[{},
  If[Length[owners] =!= Length[plan["Charts"]] ||
      !AllTrue[owners, StringQ],
    err["E6", <|"Owners" -> owners,
      "ChartCount" -> Length[plan["Charts"]],
      "Detail" -> "native arm owner list does not cover its exact chart chain"|>]];
  <|"from_exact" -> exactRationalString[plan["From"], "arm start"],
    "to_exact" -> exactRationalString[plan["To"], "arm endpoint"],
    "charts" -> owners, "topology" -> nativeArmTopology[plan]|>];

nativeCheckpointIdentity[prefix_String, payload_] := prefix <>
  IntegerString[Hash[payload, "SHA256"], 16, 64];

nativeLocalStatistics[local_Association] := Module[{stats},
  stats = DiffExp2`CppBackend`PersistentLocalStatistics[local];
  If[FailureQ[stats] || !AssociationQ[stats],
    err["E5", <|"BackendFailure" -> stats,
      "Detail" -> "could not inspect a retained native local"|>]];
  stats];

nativeLocalShape[local_Association] := Module[{stats = nativeLocalStatistics[local]},
  <|"EpsWindow" -> <|"Min" -> stats["epsilon_min"],
      "CompleteMax" -> stats["epsilon_max"]|>,
    "TWindow" -> <|"CompleteMax" -> stats["taylor_complete_max"]|>,
    "Dimension" -> stats["dimension"]|>];

nativeIdentityLattice[dimension_Integer, min_Integer, max_Integer,
    identity_String] := Module[{frame},
  If[!(min <= 0 <= max),
    err["E6", <|"EpsilonWindow" -> {min, max},
      "Detail" -> "ordinary-basis exact lattice requires epsilon^0 in the matching window"|>]];
  frame[row_, column_] := <|"min" -> min, "max" -> max,
    "coefficients" -> Table[
      ToString[Boole[row === column && power === 0], InputForm],
      {power, min, max}]|>;
  <|"schema" -> "diffexp2-exact-evaluated-epsilon-lattice-v1",
    "identity" -> identity,
    "evaluated_basis" -> Table[frame[row, column],
      {row, dimension}, {column, dimension}]|>];

nativeMatchWindow[current_Association, basis_Association,
    requested_Association] := Module[{stats, bstats, min, max},
  stats = nativeLocalStatistics[current];
  bstats = nativeLocalStatistics /@ basis["Columns"];
  min = Max[requested["Min"], stats["epsilon_min"],
    Max[Lookup[bstats, "epsilon_min"]]];
  max = Min[requested["CompleteMax"], stats["epsilon_max"],
    Min[Lookup[bstats, "epsilon_max"]]];
  If[min > max,
    err["E6", <|"Requested" -> requested,
      "Incoming" -> KeyTake[stats, {"epsilon_min", "epsilon_max"}],
      "Basis" -> (KeyTake[#, {"epsilon_min", "epsilon_max"}] & /@ bstats),
      "Detail" -> "retained match has no honest common epsilon window"|>]];
  <|"min" -> min, "max" -> max,
    "required_complete_max" -> max|>];

nativeMatchPolicy[atlas_Association, arm_String, index_Integer,
    basis_Association, current_Association, maxSteps_Integer] := Module[
  {window, checkpoint, policy, latticeIdentity, dimension = basis["Dimension"]},
  window = nativeMatchWindow[current, basis, atlas["Request", "EpsWindow"]];
  checkpoint = nativeCheckpointIdentity["de2-native-planned-match-", {
    atlas["Plan", "TilePlan"], arm, index,
    Lookup[basis["Columns"], "CheckpointIdentity"],
    Lookup[current, "CheckpointIdentity",
      Lookup[current, "checkpoint_identity", None]], window}];
  policy = <|"epsilon" -> window,
    "checkpoint_identity" -> checkpoint|>;
  If[atlas["Domain"] === "acb",
    latticeIdentity = nativeCheckpointIdentity[
      "de2-native-ordinary-lattice-", {checkpoint, dimension, window}];
    policy = Join[policy, <|
      "exact_lattice" -> nativeIdentityLattice[dimension,
        window["min"], window["max"], latticeIdentity],
      "refinement" -> <|"relative_tolerance" ->
          ("1e-" <> ToString[DiffExp2`Tolerances`Tol["MatchDigits"]]),
        "max_steps" -> maxSteps|>|>]];
  policy];

Options[PrepareNativeRegularIndependentArms] = {"Threads" -> Automatic};

PrepareNativeRegularIndependentArms[sys_Association, boundary_,
    lowerPlan_Association, upperPlan_Association, OptionsPattern[]] := Module[
  {plans, lower, upper, dimension = Length[sys["Matrix"]], values,
   epsMin, epsMax, req, anchorSystem, anchor, prepareArm, lowerData,
   upperData, sessions, anchorOwner, lowerOwners, upperOwners, planIdentity,
   nativePlan, anchorStats, sessionInfo, sessionStats, domain,
   threads = OptionValue["Threads"]},
  plans = normalizeSharedAnchor[lowerPlan, upperPlan];
  {lower, upper} = plans;
  values = nativeBoundaryValues[boundary, dimension];
  epsMin = Min[0, Min[esMin /@ values]];
  epsMax = Min[cfg["EpsilonOrder"], Min[esCM /@ values]];
  If[epsMax < 0,
    err["E6", <|"BoundaryWindow" -> {epsMin, epsMax},
      "Detail" -> "native regular basis normalization has no epsilon^0 coefficient"|>]];
  req = <|"EpsWindow" -> <|"Min" -> epsMin,
      "CompleteMax" -> epsMax|>,
    "TOrder" -> cfg["ExpansionOrder"]|>;
  anchorSystem = DiffExp2`Solve`PrepareChart[sys, First[lower["Charts"]]];
  If[!TrueQ[Lookup[anchorSystem["IndicialData"], "Regular", False]],
    err["E8", <|"Center" -> anchorSystem["Center"],
      "Detail" -> "native independent-arm anchor is not regular"|>]];
  anchor = DiffExp2`Solve`SolveNativeValueRegular[
    anchorSystem, req, values];
  anchorStats = nativeLocalStatistics[anchor];
  sessionInfo = DiffExp2`CppBackend`PersistentSessionInformation[];
  sessionStats = Lookup[sessionInfo, anchor["Session"], None];
  domain = If[AssociationQ[sessionStats],
    Lookup[sessionStats, "domain", None], None];
  If[!MemberQ[{"acb", "rational"}, domain],
    err["E5", <|"Domain" -> domain,
      "Detail" -> "native regular arm requires Acb or Rational retained locals"|>]];
  prepareArm[plan_Association] := Module[{systems, bases},
    systems = Prepend[
      DiffExp2`Solve`PrepareChart[sys, #] & /@ Rest[plan["Charts"]],
      anchorSystem];
    If[!AllTrue[systems,
        TrueQ[Lookup[# ["IndicialData"], "Regular", False]] &],
      err["E8", <|"Centers" -> Lookup[systems, "Center"],
        "Detail" -> "explicit native regular-arm seam encountered a singular chart"|>]];
    bases = Prepend[
      DiffExp2`Solve`SolveNativeRegularBasis[#, req, threads] & /@
        Rest[systems], None];
    <|"Plan" -> plan, "ChartSystems" -> systems,
      "Bases" -> bases|>];
  lowerData = prepareArm[lower];
  upperData = prepareArm[upper];
  sessions = DeleteDuplicates@Join[{anchor["Session"]},
    Cases[Join[Rest[lowerData["Bases"]], Rest[upperData["Bases"]]],
      b_Association :> b["Session"]]];
  If[Length[sessions] =!= 1,
    err["E6", <|"Sessions" -> sessions,
      "Detail" -> "prepared native atlas was split across solver sessions"|>]];
  anchorOwner = anchor["NativeChart"];
  lowerOwners = Prepend[nativeBasisOwner /@ Rest[lowerData["Bases"]],
    anchorOwner];
  upperOwners = Prepend[nativeBasisOwner /@ Rest[upperData["Bases"]],
    anchorOwner];
  planIdentity = nativeCheckpointIdentity["de2-native-independent-arms-", {
    lower, upper, lowerOwners, upperOwners, req,
    cfg["DivisionOrder"]}];
  nativePlan = DiffExp2`CppBackend`CreatePersistentTilePlan[anchor,
    nativeArmRequest[lower, lowerOwners],
    nativeArmRequest[upper, upperOwners], planIdentity,
    cfg["DivisionOrder"]];
  If[FailureQ[nativePlan] || !AssociationQ[nativePlan],
    err["E5", <|"BackendFailure" -> nativePlan,
      "Detail" -> "persistent native independent-arm planning failed"|>]];
  <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
    "Session" -> First[sessions], "Domain" -> domain,
    "Request" -> req, "Anchor" -> anchor, "Plan" -> nativePlan,
    "Lower" -> lowerData, "Upper" -> upperData,
    "PlanCheckpointIdentity" -> planIdentity|>];

nativeApplyRow[atlas_Association, chartSystem_Association,
    current_Association, cvec_List, var_Symbol, arm_String,
    tile_Integer] := Module[{shape, row, checkpoint, projected},
  shape = nativeLocalShape[current];
  row = DiffExp2`Solve`PrepareNativeRationalRow[chartSystem, shape,
    cvec, var, <|"domain" -> atlas["Domain"], "symbols" -> {}|>];
  checkpoint = nativeCheckpointIdentity["de2-native-integrand-local-", {
    atlas["PlanCheckpointIdentity"], arm, tile,
    Lookup[current, "CheckpointIdentity",
      Lookup[current, "checkpoint_identity", None]],
    row["exact_identity"]}];
  projected = DiffExp2`CppBackend`ApplyPersistentRationalRow[
    current, row, checkpoint];
  If[FailureQ[projected] || !AssociationQ[projected],
    err["E5", <|"BackendFailure" -> projected,
      "Arm" -> arm, "Tile" -> tile,
      "Detail" -> "native rational-row projection failed"|>]];
  projected];

nativeIntegrateTile[atlas_Association, arm_String, tile_Integer,
    local_Association, certifyTail_] := Module[{shape, epsilon, checkpoint, line},
  shape = nativeLocalShape[local];
  epsilon = <|"min" -> shape["EpsWindow", "Min"],
    "max" -> shape["EpsWindow", "CompleteMax"]|>;
  checkpoint = nativeCheckpointIdentity["de2-native-tile-integral-", {
    atlas["PlanCheckpointIdentity"], arm, tile,
    Lookup[local, "CheckpointIdentity",
      Lookup[local, "checkpoint_identity", None]], epsilon}];
  line = DiffExp2`CppBackend`RunPersistentTileIntegral[
    atlas["Plan"], arm, tile, local, epsilon, checkpoint, certifyTail];
  If[FailureQ[line] || !AssociationQ[line],
    err["E5", <|"BackendFailure" -> line, "Arm" -> arm,
      "Tile" -> tile, "Detail" -> "persistent native tile integration failed"|>]];
  line];

Options[RunNativeRegularIndependentArms] = {
  "CertifyTail" -> False, "MaxRefinementSteps" -> 2};

RunNativeRegularIndependentArms[atlas_Association, cvec_List,
    var_Symbol, OptionsPattern[]] := Module[
  {certify = OptionValue["CertifyTail"],
   maxSteps = OptionValue["MaxRefinementSteps"], runArm},
  If[Lookup[atlas, "Type", None] =!=
      "DiffExp2NativeRegularIndependentArmAtlas" ||
      !BooleanQ[certify] || !IntegerQ[maxSteps] || maxSteps < 0,
    err["E8", <|"AtlasType" -> Lookup[atlas, "Type", None],
      "CertifyTail" -> certify, "MaxRefinementSteps" -> maxSteps,
      "Detail" -> "native arm execution options or atlas are malformed"|>]];
  runArm[arm_String, data_Association] := Module[
    {current = atlas["Anchor"], lines = {}, matches = {}, locals = {},
     projected = {}, basis, policy, match, next, scalar, chartCount},
    chartCount = Length[data["ChartSystems"]];
    Do[
      scalar = nativeApplyRow[atlas, data["ChartSystems"][[tile]],
        current, cvec, var, arm, tile];
      AppendTo[projected, scalar];
      AppendTo[lines, nativeIntegrateTile[
        atlas, arm, tile, scalar, certify]];
      If[tile < chartCount,
        basis = data["Bases"][[tile + 1]];
        policy = nativeMatchPolicy[atlas, arm, tile, basis,
          current, maxSteps];
        match = DiffExp2`CppBackend`RunPersistentPlannedMatch[
          atlas["Plan"], arm, tile, basis["Columns"], current, policy];
        If[FailureQ[match] || !AssociationQ[match],
          err["E5", <|"BackendFailure" -> match, "Arm" -> arm,
            "Match" -> tile,
            "Detail" -> "persistent native planned match failed"|>]];
        next = DiffExp2`CppBackend`MaterializePersistentLocalMatch[match,
          nativeCheckpointIdentity["de2-native-receiving-local-", {
            policy["checkpoint_identity"], arm, tile}]];
        If[FailureQ[next] || !AssociationQ[next],
          err["E5", <|"BackendFailure" -> next, "Arm" -> arm,
            "Match" -> tile,
            "Detail" -> "persistent native match materialization failed"|>]];
        AppendTo[matches, match]; AppendTo[locals, next]; current = next],
      {tile, chartCount}];
    <|"Arm" -> arm, "FinalLocal" -> current, "Lines" -> lines,
      "Matches" -> matches, "MaterializedLocals" -> locals,
      "ProjectedLocals" -> projected|>];
  <|"Type" -> "DiffExp2NativeRegularIndependentArmRun",
    "Atlas" -> atlas, "Lower" -> runArm["lower", atlas["Lower"]],
    "Upper" -> runArm["upper", atlas["Upper"]]|>];

ReleaseNativeRegularIndependentArms[obj_Association] := Module[
  {atlas, run, lines = {}, projected = {}, materialized = {}, matches = {},
   bases, locals, responses = {}, failures, releaseAll, releaseOKQ},
  {atlas, run} = Which[
    Lookup[obj, "Type", None] ===
        "DiffExp2NativeRegularIndependentArmRun", {obj["Atlas"], obj},
    Lookup[obj, "Type", None] ===
        "DiffExp2NativeRegularIndependentArmAtlas", {obj, None},
    True, err["E8", <|"Type" -> Lookup[obj, "Type", None],
      "Detail" -> "native arm release requires an atlas or completed arm run"|>]];
  If[AssociationQ[run],
    lines = Join[run["Lower", "Lines"], run["Upper", "Lines"]];
    projected = Join[run["Lower", "ProjectedLocals"],
      run["Upper", "ProjectedLocals"]];
    materialized = Join[run["Lower", "MaterializedLocals"],
      run["Upper", "MaterializedLocals"]];
    matches = Join[run["Lower", "Matches"], run["Upper", "Matches"]]];
  bases = Flatten[Map[Lookup[#, "Columns", {}] &,
    Join[Rest[atlas["Lower", "Bases"]],
      Rest[atlas["Upper", "Bases"]]]], 1];
  releaseAll[fn_, items_List, key_] := Scan[Function[item,
    AppendTo[responses, Quiet[fn[item]]]],
    DeleteDuplicatesBy[Select[items, AssociationQ],
      Lookup[#, key, Lookup[#, ToUpperCase[StringTake[key, 1]] <>
          StringDrop[key, 1], None]] &]];
  releaseAll[DiffExp2`CppBackend`ReleasePersistentLineIntegral,
    lines, "line"];
  locals = Join[projected, materialized, bases, {atlas["Anchor"]}];
  releaseAll[DiffExp2`CppBackend`ReleasePersistentLocal, locals, "local"];
  releaseAll[DiffExp2`CppBackend`ReleasePersistentLocalMatch,
    matches, "match"];
  AppendTo[responses,
    Quiet[DiffExp2`CppBackend`ReleasePersistentTilePlan[atlas["Plan"]]]];
  releaseOKQ[response_] := AssociationQ[response] &&
    Lookup[response, "status", "error"] === "ok";
  failures = Select[responses, !TrueQ[releaseOKQ[#]] &];
  <|"Released" -> Length[responses] - Length[failures],
    "Failures" -> failures|>];

decodeLineValue[line_Association, outputDigits_Integer] := Module[
  {exported, value, coefficients, decoded},
  exported = DiffExp2`CppBackend`ExportPersistentLineIntegral[line,
    Lookup[line, "checkpoint_identity",
      Lookup[line, "CheckpointIdentity", ""]], outputDigits];
  If[FailureQ[exported] || !AssociationQ[exported],
    err["E5", <|"BackendFailure" -> exported,
      "Detail" -> "could not export final native line-integral compatibility value"|>]];
  value = exported["value"];
  coefficients = value["coefficients"];
  decoded = DiffExp2`CppBackend`DecodeScalars[coefficients, outputDigits];
  If[FailureQ[decoded],
    err["E5", <|"BackendFailure" -> decoded,
      "Detail" -> "could not decode final native line-integral coefficients"|>]];
  esNew[value["min"], decoded]];

Options[NativeRegularLineIntegral] = {"Threads" -> Automatic,
  "CertifyTail" -> False, "MaxRefinementSteps" -> 2,
  "RetainNativeState" -> True};

NativeRegularLineIntegral[sys_Association, boundary_, from_, {lo_, hi_},
    cvec_List, OptionsPattern[]] := Module[
  {lower, upper, atlas, run, digits, lowerValues, upperValues,
   lowerTotal, upperTotal, value, result, retain},
  retain = OptionValue["RetainNativeState"];
  If[!BooleanQ[retain],
    err["E8", <|"RetainNativeState" -> retain,
      "Detail" -> "RetainNativeState must be True or False"|>]];
  If[!TrueQ[lo < from < hi],
    err["E8", <|"Range" -> {lo, hi}, "Anchor" -> from,
      "Detail" -> "explicit native independent-arm line integral requires an interior anchor"|>]];
  lower = DiffExp2`Transport`SegmentLine[sys, {from, lo}];
  upper = DiffExp2`Transport`SegmentLine[sys, {from, hi}];
  atlas = PrepareNativeRegularIndependentArms[sys, boundary, lower, upper,
    "Threads" -> OptionValue["Threads"]];
  run = RunNativeRegularIndependentArms[atlas, cvec, sys["Variable"],
    "CertifyTail" -> OptionValue["CertifyTail"],
    "MaxRefinementSteps" -> OptionValue["MaxRefinementSteps"]];
  digits = cfg["WorkingPrecision"];
  lowerValues = decodeLineValue[#, digits] & /@ run["Lower", "Lines"];
  upperValues = decodeLineValue[#, digits] & /@ run["Upper", "Lines"];
  lowerTotal = Fold[esAdd, First[lowerValues], Rest[lowerValues]];
  upperTotal = Fold[esAdd, First[upperValues], Rest[upperValues]];
  value = esAdd[esScale[-1, lowerTotal], upperTotal];
  result = <|"Type" -> "DiffExp2NativeRegularLineIntegral",
    "Value" -> value, "Atlas" -> atlas, "Run" -> run,
    "CompatibilityExports" -> Length[lowerValues] + Length[upperValues]|>;
  If[retain, result,
    ReleaseNativeRegularIndependentArms[run];
    KeyDrop[Append[result, "ReleasedNativeState" -> True],
      {"Atlas", "Run"}]]];

End[];
EndPackage[];
