(* DiffExp2/NativeTransport.m -- persistent native chart-atlas and arm
   orchestration.  Wolfram owns exact input/planning metadata; no retained
   local coefficient tensor is materialized here. *)

BeginPackage["DiffExp2`NativeTransport`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Solve`", "DiffExp2`Transport`",
   "DiffExp2`CppBackend`"}];

PrepareNativeRegularIndependentArms::usage =
  "PrepareNativeRegularIndependentArms[sys,boundary,lowerPlan,upperPlan] prepares one shared retained anchor, every regular receiving basis, and one exact lower/upper native tile plan. Option Integrand->{cvec,var} derives and solves the honest epsilon halo required by polar coefficient rows. It returns only opaque native locals/bases and exact atlas metadata.";
RunNativeRegularIndependentArms::usage =
  "RunNativeRegularIndependentArms[atlas,cvec,var] precomputes one exact rational integrand row per tile, then marches the lower and upper arms concurrently in one persistent C++ request. Matching remains vector-valued, row projection is hidden, every tile and both arm sums remain native, and only the two final locals plus lower/upper/combined line handles are published.";
NativeRegularIndependentArmPlansSupportedQ::usage =
  "NativeRegularIndependentArmPlansSupportedQ[lower,upper] is the side-effect-free public-dispatch eligibility predicate for the current exact-rational native path protocol. It accepts only a shared interior anchor, entirely regular charts, and topology/geometry representable by the retained native tile planner.";
ReleaseNativeRegularIndependentArms::usage =
  "ReleaseNativeRegularIndependentArms[runOrAtlas] releases the public final-local, line-aggregate, basis, anchor, and tile-plan handles created by the explicit native regular-arm seam. Hidden matches, projected locals, and per-tile lines are reclaimed through their retained owner chains. Prepared chart/SCC caches remain session-owned and reusable.";
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
esZero = DiffExp2`EpsSeries`ESZero;
esTruncate = DiffExp2`EpsSeries`ESTruncate;

exactRationalQ[value_] := IntegerQ[value] || Head[value] === Rational;

nativeRationalEpsilonShift[expression_, physicalVar_Symbol] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], canonical, numerator,
   denominator, numeratorValuation, denominatorValuation},
  canonical = Together[expression];
  If[FreeQ[canonical, _?InexactNumberQ] &&
      TrueQ[PossibleZeroQ[canonical]], Return[0, Module]];
  numerator = Numerator[canonical];
  denominator = Denominator[canonical];
  If[!PolynomialQ[numerator, {physicalVar, eps}] ||
      !PolynomialQ[denominator, {physicalVar, eps}] ||
      !FreeQ[denominator, _?InexactNumberQ],
    err["E5", <|"Expression" -> expression,
      "Detail" -> "native regular-arm integrands must be rational in the physical variable and epsilon with an exact denominator"|>]];
  numeratorValuation = Exponent[numerator, eps, Min];
  denominatorValuation = Exponent[denominator, eps, Min];
  If[!IntegerQ[numeratorValuation] || !IntegerQ[denominatorValuation],
    err["E5", <|"Expression" -> expression,
      "Detail" -> "could not determine the exact epsilon valuation of a native regular-arm integrand"|>]];
  numeratorValuation - denominatorValuation];

nativeIntegrandMinimumShift[cvec_List, physicalVar_Symbol,
    dimension_Integer] := Module[{shifts},
  If[Length[cvec] =!= dimension,
    err["E8", <|"CoefficientCount" -> Length[cvec],
      "Dimension" -> dimension,
      "Detail" -> "native regular-arm integration requires one coefficient per physical component"|>]];
  shifts = nativeRationalEpsilonShift[#, physicalVar] & /@ cvec;
  If[shifts === {}, 0, Min[shifts]]];

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

(* This predicate deliberately performs no native preparation and catches no
   recurrence result.  It is the public API's pre-selection boundary: cases
   outside the current rational tile protocol may continue through the
   established Wolfram orchestration, while every failure after this
   predicate selected the persistent path remains loud.  Keep the checks in
   lockstep with normalizeSharedAnchor/nativeArmTopology. *)
nativeTopologyProtocolQ[plan_Association] := Quiet[Check[Module[
  {from = Lookup[plan, "From", None], to = Lookup[plan, "To", None],
   singularities, all, real, projected, nonreal, data, pairs, relevant,
   realOnArm, boundaryOnArm, records, grouped},
  If[!exactRationalQ[from] || !exactRationalQ[to], Return[False, Module]];
  singularities = Lookup[plan, "Singularities", None];
  If[!AssociationQ[singularities], Return[False, Module]];
  all = Lookup[singularities, "All", {}];
  real = Lookup[singularities, "Real", {}];
  projected = Lookup[singularities, "Projected", {}];
  If[!ListQ[all] || !ListQ[real] || !ListQ[projected],
    Return[False, Module]];
  realOnArm = Select[real, inClosedArmQ[#, from, to] &];
  boundaryOnArm = Select[
    Lookup[singularities, "ProjectionWaypoints", {}],
    inClosedArmQ[#, from, to] &];
  If[!AllTrue[Join[realOnArm, boundaryOnArm], exactRationalQ],
    Return[False, Module]];
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
  If[!AllTrue[relevant,
      exactRationalQ[#[[1]]] && exactRationalQ[#[[2]]] &&
        TrueQ[#[[2]] > 0] &], Return[False, Module]];
  records = Flatten[Lookup[Lookup[plan, "Charts", {}],
    "Prescriptions", {}], 1];
  If[!AllTrue[records, AssociationQ[#] &&
      MemberQ[{-1, 1}, Lookup[#, "Sign", None]] &],
    Return[False, Module]];
  grouped = GatherBy[records, ToString[
      Lookup[#, "ExactFactor", Lookup[#, "Factor", None]], InputForm] &];
  AllTrue[grouped,
    Length[DeleteDuplicates[Lookup[#, "Sign", None]]] === 1 &]
  ], False]];

NativeRegularIndependentArmPlansSupportedQ[lower_Association,
    upper_Association] := Quiet[Check[Module[
  {lowerAnchor, upperAnchor, charts, geometryKeys =
      {"Center", "Radius", "MatchRadius", "Scale", "LocalRadius"},
   geometry},
  If[Lookup[lower, "Direction", None] =!= -1 ||
      Lookup[upper, "Direction", None] =!= 1 ||
      !sameExactQ[Lookup[lower, "From", None],
        Lookup[upper, "From", None]] ||
      !TrueQ[Lookup[lower, "To", None] < Lookup[lower, "From", None] <
        Lookup[upper, "To", None]], Return[False, Module]];
  If[!ListQ[Lookup[lower, "Charts", None]] ||
      !ListQ[Lookup[upper, "Charts", None]] ||
      lower["Charts"] === {} || upper["Charts"] === {},
    Return[False, Module]];
  lowerAnchor = First[lower["Charts"]];
  upperAnchor = First[upper["Charts"]];
  If[!AssociationQ[lowerAnchor] || !AssociationQ[upperAnchor] ||
      !sameExactQ[Lookup[lowerAnchor, "Center", None],
        Lookup[upperAnchor, "Center", None]] ||
      !SameQ[Lookup[lowerAnchor, "Prescriptions", {}],
        Lookup[upperAnchor, "Prescriptions", {}]] ||
      !exactRationalQ[cfg["RadiusOfConvergence"]],
    Return[False, Module]];
  charts = Join[lower["Charts"], upper["Charts"]];
  If[!AllTrue[charts, AssociationQ[#] &&
      Lookup[#, "Singular", Missing["Absent"]] === False &],
    Return[False, Module]];
  geometry = Flatten[Map[Lookup[#, geometryKeys, None] &, charts]];
  If[!AllTrue[geometry, exactRationalQ], Return[False, Module]];
  If[!AllTrue[Cases[charts,
      chart_Association /; KeyExistsQ[chart, "IncomingMatchPoint"] :>
        chart["IncomingMatchPoint"]], exactRationalQ],
    Return[False, Module]];
  nativeTopologyProtocolQ[lower] && nativeTopologyProtocolQ[upper]
  ], False]];
NativeRegularIndependentArmPlansSupportedQ[___] := False;

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

nativeOpaqueCheckpoint[handle_Association, label_String] := Module[{value},
  value = Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]];
  If[!StringQ[value] || StringLength[value] == 0,
    err["E6", <|"Field" -> label,
      "Detail" -> "retained native source has no checkpoint identity"|>]];
  value];

nativeOpaqueLocalHandleQ[handle_, session_String] :=
  AssociationQ[handle] &&
  Lookup[handle, "session", Lookup[handle, "Session", None]] === session &&
  StringQ[Lookup[handle, "local", Lookup[handle, "Local", None]]] &&
  StringLength[Lookup[handle, "local", Lookup[handle, "Local", ""]]] > 0 &&
  StringQ[Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]]] &&
  StringLength[Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", ""]]] > 0;

nativeOpaqueLineHandleQ[handle_, session_String] :=
  AssociationQ[handle] &&
  Lookup[handle, "session", Lookup[handle, "Session", None]] === session &&
  StringQ[Lookup[handle, "line", Lookup[handle, "Line", None]]] &&
  StringLength[Lookup[handle, "line", Lookup[handle, "Line", ""]]] > 0 &&
  StringQ[Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]]] &&
  StringLength[Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", ""]]] > 0;

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

(* Solve's opaque Wolfram records already retain the honest solve rectangle.
   Prefer it over a native stats round-trip; top_valid can only narrow a later
   live local, so preparing against this declared rectangle remains safe. *)
nativePreparationShape[local_Association, dimension_Integer] := Module[
  {epsWindow = Lookup[local, "EpsWindow", None],
   tWindow = Lookup[local, "TWindow", None]},
  If[AssociationQ[epsWindow] && AssociationQ[tWindow] &&
      IntegerQ[Lookup[epsWindow, "Min", None]] &&
      IntegerQ[Lookup[epsWindow, "CompleteMax", None]] &&
      epsWindow["Min"] <= epsWindow["CompleteMax"] &&
      IntegerQ[Lookup[tWindow, "CompleteMax", None]] &&
      tWindow["CompleteMax"] >= 0,
    <|"EpsWindow" -> KeyTake[epsWindow, {"Min", "CompleteMax"}],
      "TWindow" -> <|"CompleteMax" -> tWindow["CompleteMax"]|>,
      "Dimension" -> dimension|>,
    nativeLocalShape[local]]];

(* A materialized receiving local owns the intersection of its incoming and
   basis-column rectangles.  Preparing a multiplier against the union
   envelope of the receiving basis is therefore conservative: the native row
   parser accepts excess epsilon/Taylor kernels, while it rejects a row that
   is even one coefficient too short for the live local. *)
nativeBasisEnvelopeShape[basis_Association] := Module[
  {columns, shapes, dimensions, epsMin, epsMax, tMax,
   dimension = Lookup[basis, "Dimension", None]},
  columns = Lookup[basis, "Columns", {}];
  If[columns === {} || !AllTrue[columns, AssociationQ] ||
      !IntegerQ[dimension] || dimension < 1,
    err["E6", <|"Basis" -> KeyTake[basis, {"Type", "Dimension"}],
      "Detail" -> "native integrand-row preparation requires a nonempty retained receiving basis"|>]];
  shapes = nativePreparationShape[#, dimension] & /@ columns;
  dimensions = DeleteDuplicates[Lookup[shapes, "Dimension", None]];
  If[Length[dimensions] =!= 1 || First[dimensions] =!= dimension,
    err["E6", <|"Dimensions" -> dimensions,
      "ExpectedDimension" -> dimension,
      "Detail" -> "retained receiving-basis columns disagree in dimension"|>]];
  epsMin = Min[Lookup[Lookup[shapes, "EpsWindow"], "Min"]];
  epsMax = Max[Lookup[Lookup[shapes, "EpsWindow"], "CompleteMax"]];
  tMax = Max[Lookup[Lookup[shapes, "TWindow"], "CompleteMax"]];
  <|"EpsWindow" -> <|"Min" -> epsMin, "CompleteMax" -> epsMax|>,
    "TWindow" -> <|"CompleteMax" -> tMax|>,
    "Dimension" -> First[dimensions]|>];

nativePreparedArmRows[atlas_Association, data_Association, cvec_List,
    var_Symbol] := Module[{systems, bases, shapes, rows},
  systems = data["ChartSystems"];
  bases = Rest[data["Bases"]];
  shapes = Prepend[nativeBasisEnvelopeShape /@ bases,
    nativePreparationShape[atlas["Anchor"], atlas["Dimension"]]];
  If[Length[systems] =!= Length[shapes],
    err["E6", <|"ChartCount" -> Length[systems],
      "ShapeCount" -> Length[shapes],
      "Detail" -> "native arm row envelopes do not reproduce its chart chain"|>]];
  rows = MapThread[DiffExp2`Solve`PrepareNativeRationalRow[
      #1, #2, cvec, var,
      <|"domain" -> atlas["Domain"], "symbols" -> {}|>] &,
    {systems, shapes}];
  If[!AllTrue[rows, AssociationQ],
    err["E5", <|"Detail" ->
      "native arm integrand-row preparation returned a malformed row"|>]];
  rows];

nativeArmExecution[atlas_Association, data_Association, cvec_List,
    var_Symbol] := Module[{bases, rows},
  bases = Rest[data["Bases"]];
  rows = nativePreparedArmRows[atlas, data, cvec, var];
  <|"receiving_basis" -> Lookup[bases, "Columns", {}],
    "integrand_rows" -> rows|>];

nativePreparedRowsMinimumShift[rows_List] := Module[
  {entries, multipliers, shifts},
  entries = Flatten[Lookup[rows, "entries", {}], 1];
  If[entries === {}, Return[0, Module]];
  If[!AllTrue[entries, AssociationQ],
    err["E5", <|"Detail" ->
      "prepared native integrand rows contain malformed entries"|>]];
  multipliers = Lookup[entries, "multiplier", None];
  shifts = If[AllTrue[multipliers, AssociationQ],
    Lookup[multipliers, "epsilon_shift", None], {None}];
  If[!AllTrue[shifts, IntegerQ],
    err["E5", <|"Detail" ->
      "prepared native integrand rows do not expose signed epsilon shifts"|>]];
  Min[shifts]];

nativeReleaseResponseHandles[atlas_Association, response_] := Module[
  {session = Lookup[atlas, "Session", None], arms, summaries, lines,
   locals, sourceRecords, sourceTokens, localToken, lineToken, unique},
  If[!AssociationQ[response], Return[Null, Module]];
  If[!StringQ[session] || StringLength[session] == 0,
    Return[Null, Module]];
  arms = Lookup[response, "arms", <||>];
  summaries = If[AssociationQ[arms], Select[Values[arms], AssociationQ], {}];
  lines = Join[
    Cases[{Lookup[response, "combined_line_result", None]}, _Association],
    Cases[Lookup[summaries, "line_result", None], _Association]];
  locals = Cases[Lookup[summaries, "final_local", None], _Association];
  sourceRecords = Join[{Lookup[atlas, "Anchor", None]},
    Flatten[Map[Lookup[#, "Columns", {}] &,
      Join[Rest[Lookup[Lookup[atlas, "Lower", <||>], "Bases", {}]],
        Rest[Lookup[Lookup[atlas, "Upper", <||>], "Bases", {}]]]], 1]];
  localToken[record_Association] := Lookup[record, "local",
    Lookup[record, "Local", None]];
  lineToken[record_Association] := Lookup[record, "line",
    Lookup[record, "Line", None]];
  sourceTokens = DeleteDuplicates@Cases[
    localToken /@ Select[sourceRecords, AssociationQ], _String];
  unique[items_List, key_String] := DeleteDuplicatesBy[items,
    Lookup[#, key, ToString[#, InputForm]] &];
  Scan[Function[line, Module[{token = lineToken[line]},
      If[StringQ[token] && StringLength[token] > 0,
        Quiet[DiffExp2`CppBackend`RunRequest[<|"schema" -> 2,
          "op" -> "integration.release", "session" -> session,
          "line" -> token|>]]]]], unique[lines, "line"]];
  Scan[Function[local, Module[{token = localToken[local]},
      If[StringQ[token] && StringLength[token] > 0 &&
          !MemberQ[sourceTokens, token],
        Quiet[DiffExp2`CppBackend`RunRequest[<|"schema" -> 2,
          "op" -> "local.release", "session" -> session,
          "local" -> token|>]]]]], unique[locals, "local"]];
  Null];

Options[PrepareNativeRegularIndependentArms] = {
  "Threads" -> Automatic, "Integrand" -> Automatic};

PrepareNativeRegularIndependentArms[sys_Association, boundary_,
    lowerPlan_Association, upperPlan_Association, OptionsPattern[]] := Module[
  {plans, lower, upper, dimension = Length[sys["Matrix"]], values,
   epsMin, epsMax, req, anchorSystem, anchor = None, prepareArm, lowerData,
   upperData, sessions, anchorOwner, lowerOwners, upperOwners, planIdentity,
   nativePlan = None, sessionInfo, sessionStats, domain, integrand,
   preparedShift, halo, targetMax = cfg["EpsilonOrder"], availableMax,
   cleanup, output, preparedBases = {},
   threads = OptionValue["Threads"]},
  plans = normalizeSharedAnchor[lowerPlan, upperPlan];
  {lower, upper} = plans;
  values = nativeBoundaryValues[boundary, dimension];
  integrand = OptionValue["Integrand"];
  preparedShift = Which[
    integrand === Automatic, 0,
    MatchQ[integrand, {_List, _Symbol}],
      nativeIntegrandMinimumShift[integrand[[1]], integrand[[2]],
        dimension],
    True, err["E8", <|"Integrand" -> integrand,
      "Detail" -> "Integrand must be Automatic or {coefficientVector, physicalVariable}"|>]];
  halo = Max[0, -preparedShift];
  epsMin = Min[0, Min[esMin /@ values]];
  availableMax = Min[esCM /@ values];
  epsMax = targetMax + halo;
  If[availableMax < epsMax,
    err["E6", <|"BoundaryCompleteMax" -> availableMax,
      "RequestedCompleteMax" -> targetMax,
      "IntegrandEpsilonShift" -> preparedShift,
      "RequiredSolveCompleteMax" -> epsMax,
      "Detail" -> "native regular-arm boundary data does not contain the epsilon halo required by its integrand pole"|>]];
  req = <|"EpsWindow" -> <|"Min" -> epsMin,
      "CompleteMax" -> epsMax|>,
    "TOrder" -> cfg["ExpansionOrder"]|>;
  cleanup[] := Module[{basisLocals, locals},
    If[AssociationQ[nativePlan],
      Quiet[DiffExp2`CppBackend`ReleasePersistentTilePlan[nativePlan]]];
    basisLocals = If[preparedBases === {}, {},
      Flatten[Lookup[preparedBases, "Columns", {}], 1]];
    locals = DeleteDuplicatesBy[
      Join[Select[basisLocals, AssociationQ],
        Cases[{anchor}, _Association]],
      Lookup[#, "local", Lookup[#, "Local", ToString[#, InputForm]]] &];
    Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &, locals];
    Null];
  output = Catch[
  anchorSystem = DiffExp2`Solve`PrepareChart[sys, First[lower["Charts"]]];
  If[!TrueQ[Lookup[anchorSystem["IndicialData"], "Regular", False]],
    err["E8", <|"Center" -> anchorSystem["Center"],
      "Detail" -> "native independent-arm anchor is not regular"|>]];
  anchor = DiffExp2`Solve`SolveNativeValueRegular[
    anchorSystem, req, values];
  sessionInfo = DiffExp2`CppBackend`PersistentSessionInformation[];
  sessionStats = Lookup[sessionInfo, anchor["Session"], None];
  domain = If[AssociationQ[sessionStats],
    Lookup[sessionStats, "domain", None], None];
  If[!MemberQ[{"acb", "rational"}, domain],
    err["E5", <|"Domain" -> domain,
      "Detail" -> "native regular arm requires Acb or Rational retained locals"|>]];
  prepareArm[plan_Association] := Module[{systems, bases, built},
    systems = Prepend[
      DiffExp2`Solve`PrepareChart[sys, #] & /@ Rest[plan["Charts"]],
      anchorSystem];
    If[!AllTrue[systems,
        TrueQ[Lookup[# ["IndicialData"], "Regular", False]] &],
      err["E8", <|"Centers" -> Lookup[systems, "Center"],
        "Detail" -> "explicit native regular-arm seam encountered a singular chart"|>]];
    bases = Prepend[Map[
      Function[system,
        built = DiffExp2`Solve`SolveNativeRegularBasis[
          system, req, threads];
        AppendTo[preparedBases, built];
        built], Rest[systems]], None];
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
    "Dimension" -> dimension,
    "Request" -> req, "Anchor" -> anchor, "Plan" -> nativePlan,
    "Lower" -> lowerData, "Upper" -> upperData,
    "TargetCompleteMax" -> targetMax,
    "PreparedIntegrandEpsilonShift" -> preparedShift,
    "PlanCheckpointIdentity" -> planIdentity|>,
  "DiffExp2Error", Function[{failure, tag},
    cleanup[]; Throw[failure, tag]]];
  output];

Options[RunNativeRegularIndependentArms] = {
  "CertifyTail" -> False, "MaxRefinementSteps" -> 2};

RunNativeRegularIndependentArms[atlas_Association, cvec_List,
    var_Symbol, OptionsPattern[]] := Module[
  {certify = OptionValue["CertifyTail"],
   maxSteps = OptionValue["MaxRefinementSteps"], lower, upper, epsilon,
   refinement, checkpointRoot, response, nativeArms, armRecord,
   backendID, sourceIdentities, rowShift, targetMax, sourceMin,
   requiredSolveMax, availableSolveMax, runMax, result},
  If[Lookup[atlas, "Type", None] =!=
      "DiffExp2NativeRegularIndependentArmAtlas" ||
      !BooleanQ[certify] || !IntegerQ[maxSteps] ||
      !TrueQ[0 <= maxSteps <= 32],
    err["E8", <|"AtlasType" -> Lookup[atlas, "Type", None],
      "CertifyTail" -> certify, "MaxRefinementSteps" -> maxSteps,
      "Detail" -> "native arm execution options or atlas are malformed"|>]];
  lower = nativeArmExecution[atlas, atlas["Lower"], cvec, var];
  upper = nativeArmExecution[atlas, atlas["Upper"], cvec, var];
  rowShift = Min[nativePreparedRowsMinimumShift[
      lower["integrand_rows"]],
    nativePreparedRowsMinimumShift[upper["integrand_rows"]]];
  targetMax = Lookup[atlas, "TargetCompleteMax", None];
  sourceMin = atlas["Request", "EpsWindow", "Min"];
  availableSolveMax = atlas["Request", "EpsWindow", "CompleteMax"];
  If[!IntegerQ[targetMax] || !IntegerQ[sourceMin] ||
      !IntegerQ[availableSolveMax],
    err["E6", <|"Detail" ->
      "native arm atlas has no finite source/output epsilon contract"|>]];
  requiredSolveMax = targetMax + Max[0, -rowShift];
  If[requiredSolveMax > availableSolveMax,
    err["E6", <|"IntegrandEpsilonShift" -> rowShift,
      "AvailableSolveCompleteMax" -> availableSolveMax,
      "RequiredSolveCompleteMax" -> requiredSolveMax,
      "Detail" -> "native arm atlas was prepared without enough epsilon halo for this integrand; prepare it with the Integrand option"|>]];
  (* Matches consume the source work frame; projected tile rows may start
     below it (epsilon poles) or entirely above the requested public order.
     The latter still needs one nonempty retained interval so the compatibility
     boundary can return the rigorously zero truncation through targetMax. *)
  (* A positive multiplier shift can move the actual first nonzero source
     row above targetMax, especially after exact leading-zero trimming.  Use
     the certified source complete edge, not its declared lower bound, so the
     projected row still has a nonempty retained interval whose lower-edge
     zero guarantee proves the public truncation is zero. *)
  runMax = availableSolveMax + Max[0, rowShift];
  epsilon = <|"min" -> sourceMin + Min[0, rowShift],
    "max" -> runMax,
    "match_required_complete_max" -> requiredSolveMax,
    "required_complete_max" -> targetMax|>;
  refinement = <|"relative_tolerance" ->
      ("1e-" <> ToString[DiffExp2`Tolerances`Tol["MatchDigits"]]),
    "max_steps" -> maxSteps|>;
  sourceIdentities = <|
    "anchor" -> nativeOpaqueCheckpoint[atlas["Anchor"], "anchor"],
    "lower_basis" -> MapIndexed[
      nativeOpaqueCheckpoint[#1, "lower basis " <>
        ToString[First[#2]]] &,
      Flatten[lower["receiving_basis"], 1]],
    "upper_basis" -> MapIndexed[
      nativeOpaqueCheckpoint[#1, "upper basis " <>
        ToString[First[#2]]] &,
      Flatten[upper["receiving_basis"], 1]]|>;
  checkpointRoot = nativeCheckpointIdentity["de2-native-arm-run-", {
    atlas["PlanCheckpointIdentity"], atlas["Domain"], sourceIdentities,
    epsilon, refinement, certify,
    Lookup[lower["integrand_rows"], "exact_identity"],
    Lookup[upper["integrand_rows"], "exact_identity"]}];
  response = DiffExp2`CppBackend`RunPersistentNativeArms[
    atlas["Plan"], atlas["Anchor"],
    <|"lower" -> lower, "upper" -> upper|>, epsilon, checkpointRoot,
    refinement, certify];
  result = Catch[
  If[FailureQ[response] || !AssociationQ[response] ||
      Lookup[response, "status", "error"] =!= "ok",
    backendID = If[AssociationQ[response],
      Lookup[response, "id", "E5"], "E5"];
    err[If[MemberQ[{"E3", "E4", "E5", "E6", "E7", "E8"},
        backendID], backendID, "E5"],
      <|"BackendFailure" -> response,
        "Detail" -> "persistent concurrent native arm execution failed"|>]];
  nativeArms = Lookup[response, "arms", None];
  If[!AssociationQ[nativeArms] ||
      Sort[Keys[nativeArms]] =!= Sort[{"lower", "upper"}] ||
      Lookup[response, "session", None] =!= atlas["Session"] ||
      !TrueQ[Lookup[response, "native_retained", False]] ||
      Lookup[response, "json_coefficients", None] =!= 0 ||
      !TrueQ[Lookup[response, "atomic_publication", False]] ||
      !nativeOpaqueLineHandleQ[
        Lookup[response, "combined_line_result", None], atlas["Session"]],
    err["E5", <|"BackendResponse" -> response,
      "Detail" -> "concurrent native arm response violated its opaque atomic-publication contract"|>]];
  armRecord[name_String] := Module[{raw = nativeArms[name]},
    If[!AssociationQ[raw] ||
        !nativeOpaqueLocalHandleQ[
          Lookup[raw, "final_local", None], atlas["Session"]] ||
        !nativeOpaqueLineHandleQ[
          Lookup[raw, "line_result", None], atlas["Session"]] ||
        !IntegerQ[Lookup[raw, "matches", None]] ||
        !IntegerQ[Lookup[raw, "tiles", None]],
      err["E5", <|"Arm" -> name, "BackendResponse" -> raw,
        "Detail" -> "concurrent native arm summary is malformed"|>]];
    <|"Arm" -> name, "FinalLocal" -> raw["final_local"],
      "Line" -> raw["line_result"], "Matches" -> raw["matches"],
      "Tiles" -> raw["tiles"],
      "ElapsedMilliseconds" -> Lookup[raw, "elapsed_ms", None]|>];
  <|"Type" -> "DiffExp2NativeRegularIndependentArmRun",
    "Atlas" -> atlas, "Lower" -> armRecord["lower"],
    "Upper" -> armRecord["upper"],
    "CombinedLine" -> response["combined_line_result"],
    "NativeSummary" -> KeyDrop[response,
      {"status", "arms", "combined_line_result"}]|>,
  "DiffExp2Error", Function[{failure, tag},
    nativeReleaseResponseHandles[atlas, response]; Throw[failure, tag]]];
  result];

ReleaseNativeRegularIndependentArms[obj_Association] := Module[
  {atlas, run, lines = {}, runLocals = {}, matches = {}, bases, locals,
   responses = {}, failures, releaseAll, releaseOKQ},
  {atlas, run} = Which[
    Lookup[obj, "Type", None] ===
        "DiffExp2NativeRegularIndependentArmRun", {obj["Atlas"], obj},
    Lookup[obj, "Type", None] ===
        "DiffExp2NativeRegularIndependentArmAtlas", {obj, None},
    True, err["E8", <|"Type" -> Lookup[obj, "Type", None],
      "Detail" -> "native arm release requires an atlas or completed arm run"|>]];
  If[AssociationQ[run],
    If[KeyExistsQ[run, "CombinedLine"],
      (* Current atomic whole-arm result: intermediate matches, projected
         locals, and tile lines were never entered in a public registry. *)
      lines = {run["CombinedLine"], run["Lower", "Line"],
        run["Upper", "Line"]};
      runLocals = {run["Lower", "FinalLocal"],
        run["Upper", "FinalLocal"]},
      (* Accept a pre-cutover run record long enough to release it safely. *)
      lines = Join[run["Lower", "Lines"], run["Upper", "Lines"]];
      runLocals = Join[run["Lower", "ProjectedLocals"],
        run["Upper", "ProjectedLocals"],
        run["Lower", "MaterializedLocals"],
        run["Upper", "MaterializedLocals"]];
      matches = Join[run["Lower", "Matches"],
        run["Upper", "Matches"]]]];
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
  locals = Join[runLocals, bases, {atlas["Anchor"]}];
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

decodeLineValue[line_Association, outputDigits_Integer,
    targetCompleteMax_Integer] := Module[
  {exported, value, coefficients, decoded, series},
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
  series = esNew[value["min"], decoded];
  Which[
    esMin[series] > targetCompleteMax, esZero[targetCompleteMax],
    esCM[series] > targetCompleteMax,
      esTruncate[series, targetCompleteMax],
    esCM[series] < targetCompleteMax,
      err["E6", <|"NativeWindow" ->
          <|"Min" -> esMin[series], "CompleteMax" -> esCM[series]|>,
        "TargetCompleteMax" -> targetCompleteMax,
        "Detail" -> "final native line integral did not cover the requested epsilon order"|>],
    True, series]];

Options[NativeRegularLineIntegral] = {"Threads" -> Automatic,
  "CertifyTail" -> False, "MaxRefinementSteps" -> 2,
  "RetainNativeState" -> True, "ArmPlans" -> Automatic};

NativeRegularLineIntegral[sys_Association, boundary_, from_, {lo_, hi_},
    cvec_List, OptionsPattern[]] := Module[
  {lower, upper, atlas, run = None, digits, value, result, retain, output,
   armPlans = OptionValue["ArmPlans"]},
  retain = OptionValue["RetainNativeState"];
  If[!BooleanQ[retain],
    err["E8", <|"RetainNativeState" -> retain,
      "Detail" -> "RetainNativeState must be True or False"|>]];
  If[!TrueQ[lo < from < hi],
    err["E8", <|"Range" -> {lo, hi}, "Anchor" -> from,
      "Detail" -> "explicit native independent-arm line integral requires an interior anchor"|>]];
  If[armPlans === Automatic,
    lower = DiffExp2`Transport`SegmentLine[sys, {from, lo}];
    upper = DiffExp2`Transport`SegmentLine[sys, {from, hi}],
    If[!MatchQ[armPlans, {_Association, _Association}],
      err["E8", <|"ArmPlans" -> armPlans,
        "Detail" -> "ArmPlans must be Automatic or an exact {lowerPlan,upperPlan} pair"|>]];
    {lower, upper} = armPlans;
    If[!sameExactQ[Lookup[lower, "From", None], from] ||
        !sameExactQ[Lookup[lower, "To", None], lo] ||
        !sameExactQ[Lookup[upper, "From", None], from] ||
        !sameExactQ[Lookup[upper, "To", None], hi],
      err["E8", <|"Range" -> {lo, hi}, "Anchor" -> from,
        "LowerPlan" -> KeyTake[lower, {"From", "To"}],
        "UpperPlan" -> KeyTake[upper, {"From", "To"}],
        "Detail" -> "precomputed native arm plans do not bind the requested anchor and range"|>]]];
  atlas = PrepareNativeRegularIndependentArms[sys, boundary, lower, upper,
    "Threads" -> OptionValue["Threads"],
    "Integrand" -> {cvec, sys["Variable"]}];
  output = Catch[
    run = RunNativeRegularIndependentArms[atlas, cvec, sys["Variable"],
      "CertifyTail" -> OptionValue["CertifyTail"],
      "MaxRefinementSteps" -> OptionValue["MaxRefinementSteps"]];
    digits = cfg["WorkingPrecision"];
    (* The signed -lower+upper aggregate is already assembled under the native
       honest Laurent-window rules.  Export it once, at the public compatibility
       boundary, instead of serializing every physical tile. *)
    value = decodeLineValue[run["CombinedLine"], digits,
      run["Atlas", "TargetCompleteMax"]];
    result = <|"Type" -> "DiffExp2NativeRegularLineIntegral",
      "Value" -> value, "Atlas" -> atlas, "Run" -> run,
      "CompatibilityExports" -> 1|>;
    If[retain, result,
      ReleaseNativeRegularIndependentArms[run];
      KeyDrop[Append[result, "ReleasedNativeState" -> True],
        {"Atlas", "Run"}]],
    "DiffExp2Error", Function[{failure, tag},
      Quiet[ReleaseNativeRegularIndependentArms[
        If[AssociationQ[run], run, atlas]]];
      Throw[failure, tag]]];
  output];

End[];
EndPackage[];
