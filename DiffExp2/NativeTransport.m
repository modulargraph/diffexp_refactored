(* DiffExp2/NativeTransport.m -- persistent native chart-atlas and arm
   orchestration.  Wolfram owns exact input/planning metadata; no retained
   local coefficient tensor is materialized here. *)

BeginPackage["DiffExp2`NativeTransport`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Solve`", "DiffExp2`Transport`",
   "DiffExp2`CppBackend`"}];

PrepareNativeRegularIndependentArms::usage =
  "PrepareNativeRegularIndependentArms[sys,boundary,lowerPlan,upperPlan] prepares one shared regular retained anchor, dispatches every non-anchor chart strictly to a regular basis or supported exact affine-Jordan SCC basis, and creates one exact lower/upper native tile plan. Unsupported singular charts fail loudly without fallback. Option Integrand->{cvec,var}, or mutually exclusive Integrands->{cvecs,var}, derives the honest global epsilon halo required by polar coefficient rows. TargetCompleteMax->Automatic uses the configured EpsilonOrder; a nonnegative integer selects the public atlas target before the coefficient halo is added to the source solve. It returns only opaque native locals/bases and exact atlas metadata.";
RunNativeRegularIndependentArms::usage =
  "RunNativeRegularIndependentArms[atlas,cvec,var] precomputes one exact rational integrand row per tile, then marches regular or supported exact affine-Jordan singular receiving charts on the lower and upper arms concurrently in one persistent C++ request. Matching remains vector-valued, row projection is hidden, every tile and both arm sums remain native, and only the two final locals plus lower/upper/combined line handles are published.";
RunNativeTransportObservableBatch::usage =
  "RunNativeTransportObservableBatch[atlas,observables,var] marches the retained lower/upper atlas exactly once, then contracts every ordered integrate, limitLower, and limitUpper observable without rematching. Each observable contains Operation, Identity, CheckpointIdentity, CoefficientVector, and Epsilon; integrate observables may additionally contain TailPolicy. Results are opaque retained line/endpoint handles and preserve request order.";
SaveNativeTransportObservableBatchCheckpoint::usage =
  "SaveNativeTransportObservableBatchCheckpoint[batch,path,identity] atomically saves one completed retained observable batch through the schema-2 native checkpoint protocol. The returned exact manifest binds the stable line, endpoint, and transport-state handles together with the pre-save transport-arm march counter.";
RestoreNativeTransportObservableBatchCheckpoint::usage =
  "RestoreNativeTransportObservableBatchCheckpoint[manifest] restores one completed retained observable batch without replaying chart preparation, matching, or transport. It validates the stable public handle manifest and transport-arm march counter before returning an exportable opaque batch.";
ReleaseNativeTransportObservableBatch::usage =
  "ReleaseNativeTransportObservableBatch[batch] releases every public line, endpoint, and retained transport-state token produced by RunNativeTransportObservableBatch, then releases its atlas anchor, bases, and tile plan. Strongly owned hidden matches and locals are reclaimed without exposing coefficient tensors.";
ExportNativeTransportObservableBatch::usage =
  "ExportNativeTransportObservableBatch[batch,outputDigits] performs the sole compatibility export for each retained line or endpoint result in request order and returns EpsSeries values. Integrate results additionally retain the validated native Scope, ErrorGuarantee, and compact ErrorEnvelope diagnostics; endpoint result shape is unchanged. No chart, match, local-sector, or tile coefficient tensor is serialized.";
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
nativeNonemptyStringQ[value_] := StringQ[value] &&
  StringLength[StringTrim[value]] > 0;

(* The retained exact path planner deliberately accepts geometry over Q
   only.  SegmentLine can nevertheless produce an ordinary chart whose
   center and local convergence radius are rational while its positive
   affine scale is a real algebraic number.  A lower dyadic floor is a safe
   bridge for precisely that case: it keeps the center, orientation, and
   branch data exact while making the represented physical disk strictly
   smaller.  No approximate comparison or RootApproximant inference enters
   this certificate. *)
$nativeScaleFloorBits = 64;

exactAlgebraicTruthQ[condition_] := TrueQ[Quiet[Check[
  FullSimplify[RootReduce[condition]], False]]];

exactRealAlgebraicQ[value_] := Module[{canonical},
  If[!FreeQ[value, _?InexactNumberQ], Return[False, Module]];
  canonical = Quiet[Check[RootReduce[value], $Failed]];
  canonical =!= $Failed && NumericQ[canonical] &&
    TrueQ[Quiet[Check[Element[canonical, Algebraics], False]]] &&
    exactAlgebraicTruthQ[Im[canonical] == 0]];

nativeInwardScaleFloor[scale_] := Module[{canonical, floor},
  canonical = Quiet[Check[RootReduce[scale], $Failed]];
  If[canonical === $Failed || !exactRealAlgebraicQ[canonical] ||
      exactRationalQ[canonical] ||
      !exactAlgebraicTruthQ[canonical > 0], Return[$Failed, Module]];
  floor = Floor[2^$nativeScaleFloorBits*canonical]/
    2^$nativeScaleFloorBits;
  If[exactRationalQ[floor] &&
      exactAlgebraicTruthQ[0 < floor < canonical], floor, $Failed]];

nativeScaleBridgePrerequisiteQ[chart_Association] := Module[
  {scale = Lookup[chart, "Scale", None],
   center = Lookup[chart, "Center", None],
   localRadius = Lookup[chart, "LocalRadius", None],
   radius = Lookup[chart, "Radius", None],
   matchRadius = Lookup[chart, "MatchRadius", None],
   roc = cfg["RadiusOfConvergence"], floor},
  If[exactRationalQ[scale],
    Return[exactRationalQ[center] && exactRationalQ[localRadius] &&
      exactRationalQ[radius] && exactRationalQ[matchRadius] &&
      scale =!= 0 && localRadius > 0 && radius > 0 && matchRadius > 0,
      Module]];
  floor = nativeInwardScaleFloor[scale];
  Lookup[chart, "Singular", Missing["Absent"]] === False &&
    exactRationalQ[center] && exactRationalQ[localRadius] &&
    localRadius > 0 && exactRationalQ[roc] && roc > 0 &&
    floor =!= $Failed && exactRealAlgebraicQ[radius] &&
    exactRealAlgebraicQ[matchRadius] &&
    exactAlgebraicTruthQ[radius == scale*localRadius] &&
    exactAlgebraicTruthQ[matchRadius == scale*roc]];

bridgeNativeRegularChartScale[chart_Association, index_Integer] := Module[
  {scale = Lookup[chart, "Scale", None], center, localRadius, radius,
   matchRadius, roc = cfg["RadiusOfConvergence"], nativeScale,
   nativeRadius, nativeMatchRadius, certificate},
  If[exactRationalQ[scale], Return[chart, Module]];
  center = Lookup[chart, "Center", None];
  localRadius = Lookup[chart, "LocalRadius", None];
  radius = Lookup[chart, "Radius", None];
  matchRadius = Lookup[chart, "MatchRadius", None];
  If[Lookup[chart, "Singular", Missing["Absent"]] =!= False,
    err["E6", <|"ChartIndex" -> index, "Center" -> center,
      "Scale" -> scale,
      "Detail" -> "nonrational native scale bridging is supported only for explicitly regular charts"|>]];
  If[!exactRationalQ[center] || !exactRationalQ[localRadius] ||
      !TrueQ[localRadius > 0] || !exactRationalQ[roc] ||
      !TrueQ[roc > 0],
    err["E6", <|"ChartIndex" -> index, "Center" -> center,
      "Scale" -> scale, "LocalRadius" -> localRadius,
      "Detail" -> "an algebraic native scale bridge requires an unchanged rational center and positive rational local radius"|>]];
  nativeScale = nativeInwardScaleFloor[scale];
  If[nativeScale === $Failed || !exactRealAlgebraicQ[radius] ||
      !exactRealAlgebraicQ[matchRadius] ||
      !exactAlgebraicTruthQ[radius == scale*localRadius] ||
      !exactAlgebraicTruthQ[matchRadius == scale*roc],
    err["E6", <|"ChartIndex" -> index, "Center" -> center,
      "Scale" -> scale, "Radius" -> radius,
      "MatchRadius" -> matchRadius, "LocalRadius" -> localRadius,
      "Detail" -> "regular algebraic scale geometry lacks the exact affine Radius=Scale LocalRadius and MatchRadius=Scale RadiusOfConvergence certificate"|>]];
  nativeRadius = Together[nativeScale*localRadius];
  nativeMatchRadius = Together[nativeScale*roc];
  If[!AllTrue[{nativeScale, nativeRadius, nativeMatchRadius},
        exactRationalQ] ||
      !exactAlgebraicTruthQ[0 < nativeRadius < radius] ||
      !exactAlgebraicTruthQ[0 < nativeMatchRadius < matchRadius],
    err["E6", <|"ChartIndex" -> index, "Center" -> center,
      "Scale" -> scale, "NativeScale" -> nativeScale,
      "Detail" -> "the exact lower rational scale did not produce a strict inward physical chart"|>]];
  certificate = <|
    "Schema" -> "diffexp2-native-inward-rational-scale-v1",
    "FloorBits" -> $nativeScaleFloorBits,
    "OriginalScale" -> scale, "NativeScale" -> nativeScale,
    "OriginalPhysicalRadius" -> radius,
    "NativePhysicalRadius" -> nativeRadius,
    "LocalRadius" -> localRadius,
    "CenterPreserved" -> center|>;
  Join[chart, <|"Scale" -> nativeScale, "Radius" -> nativeRadius,
    "MatchRadius" -> nativeMatchRadius,
    "NativeRationalScaleBridge" -> certificate|>]];

bridgeNativeRegularPlanScales[plan_Association] := Module[
  {charts = Lookup[plan, "Charts", None], bridged, records},
  If[!ListQ[charts] || charts === {},
    err["E8", <|"Detail" ->
      "native algebraic-scale normalization requires a nonempty chart chain"|>]];
  bridged = MapIndexed[bridgeNativeRegularChartScale[#1, First[#2]] &,
    charts];
  records = Cases[bridged,
    chart_Association /; KeyExistsQ[chart, "NativeRationalScaleBridge"] :>
      chart["NativeRationalScaleBridge"]];
  If[records === {}, plan,
    (* IncomingMatchPoint was built for the original affine scales.  The
       retained C++ planner owns fresh exact handoffs for the shrunken
       charts; never let stale Wolfram points masquerade as that proof. *)
    bridged = KeyDrop[#, {"IncomingMatchPoint", "SymmetricMatch"}] & /@
      bridged;
    Join[plan, <|"Charts" -> bridged,
      "NativeRationalScaleBridges" -> records,
      "NativeHandoffAuthority" -> "CppExactTilePlanner"|>]]];

nativePlannerScale[chart_Association] := Module[
  {scale = Lookup[chart, "Scale", None]},
  If[exactRationalQ[scale], scale, nativeInwardScaleFloor[scale]]];

nativeBridgedPlanPreflightQ[plan_Association] := Quiet[Check[Module[
  {charts = Lookup[plan, "Charts", None],
   from = Lookup[plan, "From", None], to = Lookup[plan, "To", None],
   dir = Lookup[plan, "Direction", None], k = cfg["DivisionOrder"],
   scales, reaches, centers, finalRadius},
  If[!ListQ[charts] || charts === {} || !MemberQ[{-1, 1}, dir] ||
      !IntegerQ[k] || k < 2, Return[False, Module]];
  scales = nativePlannerScale /@ charts;
  If[!AllTrue[scales, exactRationalQ[#] && # =!= 0 &] ||
      Length[DeleteDuplicates[Sign /@ scales]] =!= 1 ||
      !AllTrue[Lookup[charts, "LocalRadius", None],
        exactRationalQ[#] && # > 0 &], Return[False, Module]];
  centers = Lookup[charts, "Center", None];
  If[!AllTrue[centers, exactRationalQ] ||
      !exactAlgebraicTruthQ[First[centers] == from] ||
      !AllTrue[centers,
        exactAlgebraicTruthQ[dir*(# - from) >= 0] &&
          exactAlgebraicTruthQ[dir*(to - #) >= 0] &],
    Return[False, Module]];
  reaches = MapThread[
    Abs[#1]*Min[1/k, #2/2] &,
    {scales, Lookup[charts, "LocalRadius"]}];
  If[!And @@ Table[exactAlgebraicTruthQ[
        0 < dir*(centers[[i + 1]] - centers[[i]]) <=
          reaches[[i]] + reaches[[i + 1]]],
      {i, Length[centers] - 1}], Return[False, Module]];
  finalRadius = Abs[Last[scales]]*Last[Lookup[charts, "LocalRadius"]];
  exactAlgebraicTruthQ[Abs[to - Last[centers]] < finalRadius]
  ], False]];

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
  {lowerAnchor, upperAnchor, charts, commonRadius, commonMatchRadius,
   commonScale, commonAnchor, lowerOverlay, upperOverlay,
   roc = cfg["RadiusOfConvergence"]},
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
      !exactRationalQ[roc] ||
      !AllTrue[{Lookup[lowerAnchor, "Radius", None],
          Lookup[upperAnchor, "Radius", None],
          Lookup[lowerAnchor, "MatchRadius", None],
          Lookup[upperAnchor, "MatchRadius", None]}, exactRationalQ],
    Return[False, Module]];
  charts = Join[lower["Charts"], upper["Charts"]];
  If[!AllTrue[charts, AssociationQ[#] &&
      Lookup[#, "Singular", Missing["Absent"]] === False &],
    Return[False, Module]];
  If[!AllTrue[charts, nativeScaleBridgePrerequisiteQ],
    Return[False, Module]];
  commonRadius = Max[lowerAnchor["Radius"], upperAnchor["Radius"]];
  commonMatchRadius = Max[lowerAnchor["MatchRadius"],
    upperAnchor["MatchRadius"]];
  commonScale = Together[commonMatchRadius/roc];
  commonAnchor = Join[lowerAnchor, <|"Radius" -> commonRadius,
    "MatchRadius" -> commonMatchRadius, "Scale" -> commonScale,
    "LocalRadius" -> Together[commonRadius/commonScale]|>];
  lowerOverlay = Join[lower, <|"Charts" ->
    ReplacePart[lower["Charts"], 1 -> commonAnchor]|>];
  upperOverlay = Join[upper, <|"Charts" ->
    ReplacePart[upper["Charts"], 1 -> commonAnchor]|>];
  nativeTopologyProtocolQ[lower] && nativeTopologyProtocolQ[upper] &&
    nativeBridgedPlanPreflightQ[lowerOverlay] &&
    nativeBridgedPlanPreflightQ[upperOverlay]
  ], False]];
NativeRegularIndependentArmPlansSupportedQ[___] := False;

nativeBasisOwner[basis_Association] := Module[{owner},
  owner = Lookup[basis, "NativeSCC",
    Lookup[basis, "NativeChart", None]];
  If[!StringQ[owner] || StringLength[owner] == 0,
    err["E6", <|"Basis" -> KeyTake[basis,
        {"Type", "Session", "NativeSCC", "NativeChart"}],
      "Detail" -> "retained receiving basis exposes no native chart/SCC owner"|>]];
  owner];

nativeReceivingBasis[system_Association, req_Association, threads_] := Module[
  {regular = TrueQ[Lookup[
      Lookup[system, "IndicialData", <||>], "Regular", False]],
   sequence, components, built, expectedTypes},
  If[regular,
    built = DiffExp2`Solve`SolveNativeRegularBasis[
      system, req, threads];
    expectedTypes = {
      "DiffExp2NativeRegularBasis", "DiffExp2NativeSCCBasis"},
    sequence = Lookup[system, "IntegrationSequence", None];
    components = If[AssociationQ[sequence],
      Lookup[sequence, "Components", None], None];
    If[!ListQ[components] || components === {} ||
        !AllTrue[components, ListQ] ||
        Sort[Flatten[components]] =!= Range[system["SystemSize"]],
      err["E8", <|"Center" -> Lookup[system, "Center", None],
        "IntegrationSequence" -> sequence,
        "Detail" -> "singular native receiving chart has no complete exact SCC certificate; no alternate basis solver is selected"|>]];
    (* SolveNativeSCCBasis is the strict admission gate for singular data.
       It accepts only a certified exact affine-Jordan block DAG (including
       one SCC) and propagates every unsupported/CASE-P failure unchanged. *)
    built = DiffExp2`Solve`SolveNativeSCCBasis[
      system, req, threads];
    expectedTypes = {"DiffExp2NativeSCCBasis"}];
  If[!AssociationQ[built] ||
      !MemberQ[expectedTypes, Lookup[built, "Type", None]] ||
      Lookup[built, "Dimension", None] =!= system["SystemSize"] ||
      !ListQ[Lookup[built, "Columns", None]] ||
      Length[built["Columns"]] =!= system["SystemSize"] ||
      Lookup[built["Columns"], "BasisIndex", {}] =!=
        Range[system["SystemSize"]],
    If[AssociationQ[built] && ListQ[Lookup[built, "Columns", None]],
      Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
        built["Columns"]]];
    err["E6", <|"Center" -> Lookup[system, "Center", None],
      "Basis" -> If[AssociationQ[built],
        KeyTake[built, {"Type", "Dimension", "Session",
          "NativeChart", "NativeSCC"}], built],
      "Detail" -> "native receiving-chart dispatcher did not return one complete ordered opaque basis"|>]];
  built];

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
  "Threads" -> Automatic, "Integrand" -> Automatic,
  "Integrands" -> Automatic, "TargetCompleteMax" -> Automatic};

PrepareNativeRegularIndependentArms[sys_Association, boundary_,
    lowerPlan_Association, upperPlan_Association, OptionsPattern[]] := Module[
  {plans, lower, upper, dimension = Length[sys["Matrix"]], values,
   epsMin, epsMax, req, anchorSystem, anchor = None, prepareArm, lowerData,
   upperData, sessions, anchorOwner, lowerOwners, upperOwners, planIdentity,
   nativePlan = None, geometryAudit = None, sessionInfo, sessionStats,
   domain, integrand,
   preparedShift, halo, targetMax, targetOption =
     OptionValue["TargetCompleteMax"], availableMax,
   cleanup, output, preparedBases = {},
   containsSingularReceivingCharts = False,
   integrands = OptionValue["Integrands"],
   threads = OptionValue["Threads"]},
  If[targetOption =!= Automatic &&
      (!IntegerQ[targetOption] || targetOption < 0),
    err["E8", <|"TargetCompleteMax" -> targetOption,
      "Detail" -> "TargetCompleteMax must be Automatic or a nonnegative integer"|>]];
  targetMax = If[targetOption === Automatic,
    cfg["EpsilonOrder"], targetOption];
  plans = normalizeSharedAnchor[lowerPlan, upperPlan];
  plans = bridgeNativeRegularPlanScales /@ plans;
  {lower, upper} = plans;
  values = nativeBoundaryValues[boundary, dimension];
  integrand = OptionValue["Integrand"];
  If[integrand =!= Automatic && integrands =!= Automatic,
    err["E8", <|"Integrand" -> integrand, "Integrands" -> integrands,
      "Detail" -> "Integrand and Integrands are mutually exclusive native halo declarations"|>]];
  preparedShift = Which[
    integrand === Automatic && integrands === Automatic, 0,
    MatchQ[integrand, {_List, _Symbol}],
      nativeIntegrandMinimumShift[integrand[[1]], integrand[[2]],
        dimension],
    MatchQ[integrands, {_List, _Symbol}] &&
        AllTrue[integrands[[1]], ListQ],
      If[integrands[[1]] === {}, 0,
        Min[nativeIntegrandMinimumShift[#, integrands[[2]],
            dimension] & /@ DeleteDuplicates[integrands[[1]], SameQ]]],
    True, err["E8", <|"Integrand" -> integrand,
      "Integrands" -> integrands,
      "Detail" -> "Integrand must be Automatic or {coefficientVector,physicalVariable}; Integrands must be Automatic or {coefficientVectors,physicalVariable}"|>]];
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
      "Detail" -> "native independent-arm atlas requires Acb or Rational retained locals"|>]];
  prepareArm[plan_Association] := Module[{systems, bases, built, kinds},
    systems = Prepend[
      DiffExp2`Solve`PrepareChart[sys, #] & /@ Rest[plan["Charts"]],
      anchorSystem];
    kinds = Prepend[Map[
      If[TrueQ[Lookup[Lookup[#, "IndicialData", <||>],
          "Regular", False]], "Regular", "SingularSCC"] &,
      Rest[systems]], "Anchor"];
    If[MemberQ[kinds, "SingularSCC"],
      containsSingularReceivingCharts = True];
    bases = Prepend[Map[
      Function[system,
        built = nativeReceivingBasis[system, req, threads];
        AppendTo[preparedBases, built];
        built], Rest[systems]], None];
    <|"Plan" -> plan, "ChartSystems" -> systems,
      "Bases" -> bases, "BasisKinds" -> kinds|>];
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
  If[FailureQ[nativePlan] || !AssociationQ[nativePlan] ||
      Lookup[nativePlan, "status", "error"] =!= "ok" ||
      !nativeNonemptyStringQ[Lookup[nativePlan, "session", None]] ||
      !nativeNonemptyStringQ[Lookup[nativePlan, "tile_plan", None]] ||
      !AssociationQ[Lookup[nativePlan, "lower", None]] ||
      !AssociationQ[Lookup[nativePlan, "upper", None]],
    err["E5", <|"BackendFailure" -> nativePlan,
      "Detail" -> "persistent native independent-arm planning failed"|>]];
  (* A successful tile.plan response is published only after C++
     validate_exact_arm_plan has proved every exact match and tile endpoint
     strictly inside its retained rational chart.  Each bridge certificate
     above proves 0 < native Scale < original Scale with center and positive
     LocalRadius unchanged, so those native disks are exact subsets of the
     original algebraic disks.  This composes the two proofs without a
     duplicate Wolfram geometry walk. *)
  geometryAudit = <|
    "Schema" -> "diffexp2-native-exact-geometry-proof-v1",
    "BackendExactPlanValidated" -> True,
    "NativeDisksContainedInOriginal" -> True,
    "HandoffAuthority" -> "CppExactTilePlanner",
    "LowerBridgeCount" -> Length[Lookup[lower,
      "NativeRationalScaleBridges", {}]],
    "UpperBridgeCount" -> Length[Lookup[upper,
      "NativeRationalScaleBridges", {}]]|>;
  <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
    "Session" -> First[sessions], "Domain" -> domain,
    "Dimension" -> dimension,
    "Request" -> req, "Anchor" -> anchor, "Plan" -> nativePlan,
    "Lower" -> lowerData, "Upper" -> upperData,
    "NativeGeometryAudit" -> geometryAudit,
    "ContainsSingularReceivingCharts" ->
      containsSingularReceivingCharts,
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

nativeTransportStateHandleQ[state_, session_String, arm_String] :=
  AssociationQ[state] &&
  Lookup[state, "session", Lookup[state, "Session", None]] === session &&
  Lookup[state, "arm", Lookup[state, "Arm", None]] === arm &&
  nativeNonemptyStringQ[Lookup[state, "transport_state",
    Lookup[state, "TransportState", None]]] &&
  nativeNonemptyStringQ[Lookup[state, "checkpoint_identity",
    Lookup[state, "CheckpointIdentity", None]]] &&
  nativeNonemptyStringQ[Lookup[state, "provenance_identity",
    Lookup[state, "ProvenanceIdentity", None]]];

nativeOpaqueEndpointHandleQ[handle_, session_String] :=
  AssociationQ[handle] &&
  Lookup[handle, "session", Lookup[handle, "Session", None]] === session &&
  nativeNonemptyStringQ[Lookup[handle, "endpoint",
    Lookup[handle, "Endpoint", None]]] &&
  nativeNonemptyStringQ[Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]]];

nativeOrderedObservableOutputs[raw_, expected_List, session_String,
    kind_String] := Module[{records, handlesOK},
  records = If[AssociationQ[raw], Lookup[raw,
    If[kind === "line", "lines", "endpoints"], None], None];
  If[!ListQ[records],
    err["E5", <|"Kind" -> kind, "BackendResponse" -> raw,
      "Detail" -> "native observable response has no ordered result list"|>]];
  handlesOK = Switch[kind,
    "line", AllTrue[records, nativeOpaqueLineHandleQ[#, session] &],
    "endpoint", AllTrue[records,
      nativeOpaqueEndpointHandleQ[#, session] &],
    _, False];
  If[Length[records] =!= Length[expected] ||
      !handlesOK ||
      Lookup[records, "request_index"] =!= Range[0, Length[expected] - 1] ||
      Lookup[records, "observable_identity"] =!=
        Lookup[expected, "Identity"] ||
      Lookup[records, "checkpoint_identity"] =!=
        Lookup[expected, "CheckpointIdentity"],
    err["E5", <|"Kind" -> kind, "BackendResponse" -> raw,
      "ExpectedIdentities" -> Lookup[expected, "Identity"],
      "Detail" -> "native observable outputs changed request order, identity, checkpoint, session, or opaque handle shape"|>]];
  records];

nativeBatchObservable[raw_, dimension_Integer] := Module[
  {baseKeys = {"Operation", "Identity", "CheckpointIdentity",
      "CoefficientVector", "Epsilon"}, keys, operation, epsilon,
   tailPolicy},
  If[!AssociationQ[raw],
    err["E8", <|"Observable" -> raw,
      "Detail" -> "native transport batch observables must be associations"|>]];
  operation = Lookup[raw, "Operation", None];
  keys = Sort[Keys[raw]];
  If[!MemberQ[{"integrate", "limitLower", "limitUpper"}, operation] ||
      (operation === "integrate" && keys =!= Sort[baseKeys] &&
        keys =!= Sort[Append[baseKeys, "TailPolicy"]]) ||
      (operation =!= "integrate" && keys =!= Sort[baseKeys]),
    err["E8", <|"Operation" -> operation, "Keys" -> Keys[raw],
      "Detail" -> "native transport observables require exactly Operation, Identity, CheckpointIdentity, CoefficientVector, and Epsilon; only integrate accepts optional TailPolicy"|>]];
  epsilon = raw["Epsilon"];
  tailPolicy = Lookup[raw, "TailPolicy", "stored"];
  If[!nativeNonemptyStringQ[raw["Identity"]] ||
      !nativeNonemptyStringQ[raw["CheckpointIdentity"]] ||
      !ListQ[raw["CoefficientVector"]] ||
      Length[raw["CoefficientVector"]] =!= dimension ||
      !AssociationQ[epsilon] ||
      Sort[Keys[epsilon]] =!= Sort[{"Min", "Max",
        "RequiredCompleteMax"}] ||
      !AllTrue[Lookup[epsilon, {"Min", "Max",
          "RequiredCompleteMax"}], IntegerQ] ||
      !TrueQ[epsilon["Min"] <= epsilon["RequiredCompleteMax"] <=
        epsilon["Max"]] ||
      (operation === "integrate" &&
        !MemberQ[{"stored", "attempt", "require"}, tailPolicy]),
    err["E8", <|"Identity" -> Lookup[raw, "Identity", None],
      "Operation" -> operation,
      "Detail" -> "native transport observable identities, coefficient dimension, epsilon window, or tail policy are malformed"|>]];
  <|"Operation" -> operation, "Identity" -> raw["Identity"],
    "CheckpointIdentity" -> raw["CheckpointIdentity"],
    "CoefficientVector" -> raw["CoefficientVector"],
    "Epsilon" -> epsilon, "TailPolicy" -> tailPolicy|>];

nativePrepareBatchObservable[observable_Association,
    prepareRows_] := Module[{lowerRows = {}, upperRows = {}, relevant,
   shift, coefficients = observable["CoefficientVector"]},
  Switch[observable["Operation"],
    "integrate",
      lowerRows = prepareRows["lower", coefficients];
      upperRows = prepareRows["upper", coefficients],
    "limitLower",
      lowerRows = prepareRows["lower", coefficients],
    "limitUpper",
      upperRows = prepareRows["upper", coefficients]];
  relevant = Join[lowerRows, upperRows];
  If[relevant === {},
    err["E6", <|"Identity" -> observable["Identity"],
      "Detail" -> "native transport observable prepared no retained chart row"|>]];
  shift = nativePreparedRowsMinimumShift[relevant];
  Join[observable, <|"LowerRows" -> lowerRows,
    "UpperRows" -> upperRows, "MinimumEpsilonShift" -> shift|>]];

nativeBatchReleasePublished[states_, pair_, lowerEndpoints_,
    upperEndpoints_] := Module[{lines, endpoints},
  lines = If[AssociationQ[pair], Lookup[pair, "lines", {}], {}];
  endpoints = Join[
    If[AssociationQ[lowerEndpoints],
      Lookup[lowerEndpoints, "endpoints", {}], {}],
    If[AssociationQ[upperEndpoints],
      Lookup[upperEndpoints, "endpoints", {}], {}]];
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLineIntegral[#]] &,
    Select[lines, AssociationQ]];
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentEndpoint[#]] &,
    Select[endpoints, AssociationQ]];
  If[AssociationQ[states],
    Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentTransportArm[#]] &,
      Select[Values[states], AssociationQ]]];
  Null];

Options[RunNativeTransportObservableBatch] = {
  "MaxRefinementSteps" -> 2};

RunNativeTransportObservableBatch[atlas_Association, observables_List,
    var_Symbol, OptionsPattern[]] := Module[
  {maxSteps = OptionValue["MaxRefinementSteps"], normalized, identities,
   checkpoints, prepared, sourceMin, availableMax, projectedRequired,
   statePublicRequired,
   matchRequired, preparedShift, declaredShift, epsilon, refinement,
   lowerBasis, upperBasis, checkpointRoot, march, states = <||>,
   lowerState, upperState, integrates, lowerLimits, upperLimits,
   pair = None, lowerEndpoints = None, upperEndpoints = None,
   pairRecords = {}, lowerRecords = {}, upperRecords = {},
   pairByIdentity = <||>, lowerByIdentity = <||>, upperByIdentity = <||>,
   resultRecords, prepareRows, output},
  If[Lookup[atlas, "Type", None] =!=
        "DiffExp2NativeRegularIndependentArmAtlas" ||
      !IntegerQ[maxSteps] || !TrueQ[0 <= maxSteps <= 32],
    err["E8", <|"AtlasType" -> Lookup[atlas, "Type", None],
      "MaxRefinementSteps" -> maxSteps,
      "Detail" -> "native transport observable batch received a malformed atlas or refinement bound"|>]];
  normalized = nativeBatchObservable[#, atlas["Dimension"]] & /@
    observables;
  identities = Lookup[normalized, "Identity"];
  checkpoints = Lookup[normalized, "CheckpointIdentity"];
  If[DuplicateFreeQ[identities] =!= True ||
      DuplicateFreeQ[checkpoints] =!= True,
    err["E8", <|"Identities" -> identities,
      "CheckpointIdentities" -> checkpoints,
      "Detail" -> "native transport observable identities and checkpoint identities must each be pairwise unique"|>]];
  If[normalized === {},
    Return[<|"Type" -> "DiffExp2NativeTransportObservableBatch",
      "Atlas" -> atlas, "States" -> <||>, "Results" -> {},
      "NativeMarches" -> 0, "CompatibilityExports" -> 0|>, Module]];
  prepareRows[arm_String, coefficients_List] :=
    prepareRows[arm, coefficients] = Switch[arm,
      "lower", nativePreparedArmRows[atlas, atlas["Lower"],
        coefficients, var],
      "upper", nativePreparedArmRows[atlas, atlas["Upper"],
        coefficients, var]];
  prepared = nativePrepareBatchObservable[#, prepareRows] & /@ normalized;
  sourceMin = atlas["Request", "EpsWindow", "Min"];
  availableMax = atlas["Request", "EpsWindow", "CompleteMax"];
  projectedRequired = Max[Lookup[Lookup[prepared, "Epsilon"],
    "RequiredCompleteMax"]];
  (* The retained state gate belongs to the unprojected source.  A polar row
     may legitimately deliver only eps^-1 from a source that starts at eps^0;
     do not make that projected top the source state's lower completeness
     edge.  Positive shifts are deliberately conservative: they do not buy a
     lower source solve order, so every requested nonnegative public top must
     still fit the prepared source atlas. *)
  statePublicRequired = Max[sourceMin, projectedRequired];
  matchRequired = Max[statePublicRequired,
    Max[(#["Epsilon", "RequiredCompleteMax"] -
          #["MinimumEpsilonShift"]) & /@ prepared]];
  preparedShift = Min[Lookup[prepared, "MinimumEpsilonShift"]];
  declaredShift = Lookup[atlas, "PreparedIntegrandEpsilonShift", 0];
  If[!IntegerQ[sourceMin] || !IntegerQ[availableMax] ||
      !IntegerQ[declaredShift] || preparedShift < declaredShift ||
      statePublicRequired > availableMax || matchRequired > availableMax,
    err["E6", <|"PreparedMinimumEpsilonShift" -> preparedShift,
      "AtlasDeclaredMinimumEpsilonShift" -> declaredShift,
      "ProjectedRequiredCompleteMax" -> projectedRequired,
      "StatePublicRequiredCompleteMax" -> statePublicRequired,
      "MatchRequiredCompleteMax" -> matchRequired,
      "AvailableSolveCompleteMax" -> availableMax,
      "Detail" -> "native atlas does not contain the global epsilon halo required by this observable batch; prepare it with Integrands"|>]];
  epsilon = <|"min" -> sourceMin, "max" -> availableMax,
    "required_complete_max" -> statePublicRequired,
    "match_required_complete_max" -> matchRequired|>;
  refinement = <|"relative_tolerance" ->
      ("1e-" <> ToString[DiffExp2`Tolerances`Tol["MatchDigits"]]),
    "max_steps" -> maxSteps|>;
  lowerBasis = Lookup[Rest[atlas["Lower", "Bases"]], "Columns", {}];
  upperBasis = Lookup[Rest[atlas["Upper", "Bases"]], "Columns", {}];
  checkpointRoot = nativeCheckpointIdentity["de2-native-observable-batch-", {
    atlas["PlanCheckpointIdentity"], atlas["Domain"],
    nativeOpaqueCheckpoint[atlas["Anchor"], "anchor"], epsilon,
    refinement, Map[KeyTake[#, {"Operation", "Identity",
      "CheckpointIdentity", "Epsilon", "TailPolicy",
      "MinimumEpsilonShift"}] &,
      prepared],
    Map[Lookup[#, "exact_identity"] &,
      Flatten[Join[Lookup[prepared, "LowerRows"],
        Lookup[prepared, "UpperRows"]], 1]]}];
  output = Catch[
    march = DiffExp2`CppBackend`RunPersistentTransportArms[
      atlas["Plan"], atlas["Anchor"], <|
        "lower" -> <|"receiving_basis" -> lowerBasis|>,
        "upper" -> <|"receiving_basis" -> upperBasis|>|>,
      epsilon, checkpointRoot <> ":march", refinement];
    states = If[AssociationQ[march] &&
        AssociationQ[Lookup[march, "states", None]],
      march["states"], <||>];
    If[FailureQ[march] || !AssociationQ[march] ||
        Lookup[march, "status", "error"] =!= "ok" ||
        Lookup[march, "session", None] =!= atlas["Session"] ||
        !TrueQ[Lookup[march, "native_retained", False]] ||
        Lookup[march, "json_coefficients", None] =!= 0,
      err["E5", <|"BackendFailure" -> march,
        "Detail" -> "persistent two-arm observable-batch march failed"|>]];
    lowerState = Lookup[states, "lower", <||>];
    upperState = Lookup[states, "upper", <||>];
    If[Sort[Keys[states]] =!= {"lower", "upper"} ||
        !nativeTransportStateHandleQ[lowerState, atlas["Session"],
          "lower"] ||
        !nativeTransportStateHandleQ[upperState, atlas["Session"],
          "upper"],
      err["E5", <|"BackendResponse" -> march,
        "Detail" -> "two-arm observable-batch march did not return exact lower/upper retained states"|>]];
    integrates = Select[prepared, #["Operation"] === "integrate" &];
    lowerLimits = Select[prepared, #["Operation"] === "limitLower" &];
    upperLimits = Select[prepared, #["Operation"] === "limitUpper" &];
    If[integrates =!= {},
      pair = DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
        lowerState, upperState,
        Map[<|"Identity" -> #["Identity"],
            "CheckpointIdentity" -> #["CheckpointIdentity"],
            "LowerIntegrandRows" -> #["LowerRows"],
            "UpperIntegrandRows" -> #["UpperRows"],
            "Epsilon" -> #["Epsilon"],
            "TailPolicy" -> #["TailPolicy"]|> &, integrates],
        checkpointRoot <> ":integrals"];
      If[FailureQ[pair] || !AssociationQ[pair] ||
          Lookup[pair, "status", "error"] =!= "ok",
        err["E5", <|"BackendFailure" -> pair,
          "Detail" -> "native paired observable contraction failed"|>]];
      pairRecords = nativeOrderedObservableOutputs[
        pair, integrates, atlas["Session"], "line"];
      pairByIdentity = AssociationThread[Lookup[integrates, "Identity"],
        pairRecords]];
    If[lowerLimits =!= {},
      lowerEndpoints =
        DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
          lowerState,
          Map[<|"Identity" -> #["Identity"],
              "CheckpointIdentity" -> #["CheckpointIdentity"],
              "IntegrandRow" -> Last[#["LowerRows"]],
              "Epsilon" -> #["Epsilon"]|> &, lowerLimits],
          checkpointRoot <> ":lower-endpoints"];
      If[FailureQ[lowerEndpoints] || !AssociationQ[lowerEndpoints] ||
          Lookup[lowerEndpoints, "status", "error"] =!= "ok",
        err["E5", <|"BackendFailure" -> lowerEndpoints,
          "Detail" -> "native lower endpoint observable batch failed"|>]];
      lowerRecords = nativeOrderedObservableOutputs[
        lowerEndpoints, lowerLimits, atlas["Session"], "endpoint"];
      lowerByIdentity = AssociationThread[Lookup[lowerLimits, "Identity"],
        lowerRecords]];
    If[upperLimits =!= {},
      upperEndpoints =
        DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
          upperState,
          Map[<|"Identity" -> #["Identity"],
              "CheckpointIdentity" -> #["CheckpointIdentity"],
              "IntegrandRow" -> Last[#["UpperRows"]],
              "Epsilon" -> #["Epsilon"]|> &, upperLimits],
          checkpointRoot <> ":upper-endpoints"];
      If[FailureQ[upperEndpoints] || !AssociationQ[upperEndpoints] ||
          Lookup[upperEndpoints, "status", "error"] =!= "ok",
        err["E5", <|"BackendFailure" -> upperEndpoints,
          "Detail" -> "native upper endpoint observable batch failed"|>]];
      upperRecords = nativeOrderedObservableOutputs[
        upperEndpoints, upperLimits, atlas["Session"], "endpoint"];
      upperByIdentity = AssociationThread[Lookup[upperLimits, "Identity"],
        upperRecords]];
    resultRecords = MapIndexed[Function[{observable, position},
      Join[<|"RequestIndex" -> First[position] - 1,
        "Operation" -> observable["Operation"],
        "Identity" -> observable["Identity"],
        "CheckpointIdentity" -> observable["CheckpointIdentity"],
        "Epsilon" -> observable["Epsilon"]|>,
       <|If[observable["Operation"] === "integrate", "Line",
          "Endpoint"] -> Switch[observable["Operation"],
            "integrate", pairByIdentity[observable["Identity"]],
            "limitLower", lowerByIdentity[observable["Identity"]],
            "limitUpper", upperByIdentity[observable["Identity"]]]|>]],
      prepared];
    <|"Type" -> "DiffExp2NativeTransportObservableBatch",
      "Atlas" -> atlas, "States" -> states,
      "Results" -> resultRecords, "NativeMarches" -> 2,
      "CompatibilityExports" -> 0,
      "NativeSummary" -> <|
        "March" -> KeyDrop[march, {"status", "states"}],
        "Pair" -> If[AssociationQ[pair],
          KeyDrop[pair, {"status", "lines"}], None],
        "LowerEndpoints" -> If[AssociationQ[lowerEndpoints],
          KeyDrop[lowerEndpoints, {"status", "endpoints"}], None],
        "UpperEndpoints" -> If[AssociationQ[upperEndpoints],
          KeyDrop[upperEndpoints, {"status", "endpoints"}], None]|>|>,
    "DiffExp2Error", Function[{failure, tag},
      nativeBatchReleasePublished[states, pair, lowerEndpoints,
        upperEndpoints]; Throw[failure, tag]]];
  output];

nativeDecodeBatchExport[exported_Association, outputDigits_Integer,
    requiredCompleteMax_Integer] := Module[
  {value, coefficients, decoded, series},
  If[Lookup[exported, "status", "error"] =!= "ok" ||
      !AssociationQ[Lookup[exported, "value", None]],
    err["E5", <|"BackendFailure" -> exported,
      "Detail" -> "native observable compatibility export failed"|>]];
  value = exported["value"];
  coefficients = Lookup[value, "coefficients", None];
  If[!IntegerQ[Lookup[value, "min", None]] ||
      !IntegerQ[Lookup[value, "max", None]] ||
      !ListQ[coefficients] ||
      Length[coefficients] =!= value["max"] - value["min"] + 1,
    err["E5", <|"BackendResponse" -> exported,
      "Detail" -> "native observable export returned a malformed epsilon vector"|>]];
  decoded = DiffExp2`CppBackend`DecodeScalars[coefficients, outputDigits];
  If[FailureQ[decoded] || !ListQ[decoded],
    err["E5", <|"BackendFailure" -> decoded,
      "Detail" -> "native observable export coefficients could not be decoded"|>]];
  series = esNew[value["min"], decoded];
  If[esCM[series] < requiredCompleteMax,
    err["E6", <|"NativeWindow" -> <|"Min" -> esMin[series],
        "CompleteMax" -> esCM[series]|>,
      "RequiredCompleteMax" -> requiredCompleteMax,
      "Detail" -> "native observable export does not cover its required epsilon maximum"|>]];
  series];

nativeLineExportCertification[exported_Association] := Module[
  {scope = Lookup[exported, "scope", None],
   guarantee = Lookup[exported, "error_guarantee", None], value,
   envelope, allowedGuarantees = {"none", "advisory", "certified"},
   envelopeKeys = {"min", "max", "guarantee",
     "absolute_upper_approx", "bound_encoding", "provenance"}, bounds},
  value = Lookup[exported, "value", None];
  If[!MemberQ[{"stored_truncation",
        "full_local_with_certified_tail"}, scope] ||
      !MemberQ[allowedGuarantees, guarantee] || !AssociationQ[value],
    err["E5", <|"Scope" -> scope, "ErrorGuarantee" -> guarantee,
      "Detail" -> "native line export omitted or malformed its integration certification summary"|>]];
  envelope = Lookup[value, "error", None];
  If[envelope === None,
    If[guarantee =!= "none" ||
        scope === "full_local_with_certified_tail",
      err["E5", <|"Scope" -> scope, "ErrorGuarantee" -> guarantee,
        "Detail" -> "native line export certification requires a returned error envelope"|>]],
    If[!AssociationQ[envelope] ||
        Sort[Keys[envelope]] =!= Sort[envelopeKeys] ||
        !IntegerQ[Lookup[envelope, "min", None]] ||
        !IntegerQ[Lookup[envelope, "max", None]] ||
        envelope["min"] > envelope["max"] ||
        Lookup[envelope, "guarantee", None] =!= guarantee ||
        !StringQ[Lookup[envelope, "bound_encoding", None]] ||
        StringLength[envelope["bound_encoding"]] === 0 ||
        !StringQ[Lookup[envelope, "provenance", None]],
      err["E5", <|"Scope" -> scope, "ErrorGuarantee" -> guarantee,
        "ErrorEnvelope" -> envelope,
        "Detail" -> "native line export returned a malformed error envelope"|>]];
    bounds = envelope["absolute_upper_approx"];
    If[!ListQ[bounds] ||
        Length[bounds] =!= envelope["max"] - envelope["min"] + 1 ||
        !AllTrue[bounds, NumberQ[#] && TrueQ[# >= 0] &] ||
        !IntegerQ[Lookup[value, "min", None]] ||
        !IntegerQ[Lookup[value, "max", None]] ||
        envelope["min"] < value["min"] ||
        envelope["max"] > value["max"],
      err["E5", <|"ValueWindow" ->
          Lookup[value, {"min", "max"}, None],
        "ErrorEnvelope" -> envelope,
        "Detail" -> "native line export error bounds do not match a finite subwindow of the returned epsilon vector"|>]]];
  If[scope === "full_local_with_certified_tail" &&
      guarantee =!= "certified",
    err["E5", <|"Scope" -> scope, "ErrorGuarantee" -> guarantee,
      "Detail" -> "full-local line scope lacks a certified native error guarantee"|>]];
  <|"Scope" -> scope, "ErrorGuarantee" -> guarantee,
    "ErrorEnvelope" -> envelope|>];

ExportNativeTransportObservableBatch[batch_Association,
    outputDigits_Integer] := Module[{results, exported},
  If[Lookup[batch, "Type", None] =!=
        "DiffExp2NativeTransportObservableBatch" || outputDigits < 1,
    err["E8", <|"Type" -> Lookup[batch, "Type", None],
      "OutputDigits" -> outputDigits,
      "Detail" -> "native observable export requires a completed batch and positive output digits"|>]];
  results = Lookup[batch, "Results", {}];
  exported = Map[Function[result, Module[{raw, decoded},
    raw = Switch[result["Operation"],
      "integrate",
        DiffExp2`CppBackend`ExportPersistentLineIntegral[
          result["Line"], result["CheckpointIdentity"], outputDigits],
      "limitLower" | "limitUpper",
        DiffExp2`CppBackend`ExportPersistentEndpoint[
          result["Endpoint"], result["CheckpointIdentity"],
          outputDigits]];
    decoded = Append[result, "Value" -> nativeDecodeBatchExport[raw,
      outputDigits, result["Epsilon", "RequiredCompleteMax"]]];
    If[result["Operation"] === "integrate",
      Join[decoded, nativeLineExportCertification[raw]], decoded]]],
    results];
  Join[batch, <|"ExportedResults" -> exported,
    "CompatibilityExports" -> Length[exported]|>]];

nativeObservableCheckpointResult[result_Association,
    session_String] := Module[
  {operation, kind, handle, token, checkpoint, provenance},
  operation = Lookup[result, "Operation", None];
  {kind, handle} = Switch[operation,
    "integrate", {"line", Lookup[result, "Line", None]},
    "limitLower" | "limitUpper",
      {"endpoint", Lookup[result, "Endpoint", None]},
    _, err["E8", <|"Operation" -> operation,
      "Detail" -> "completed native checkpoint contains an unsupported observable operation"|>]];
  checkpoint = Lookup[result, "CheckpointIdentity", None];
  token = Switch[kind,
    "line",
      If[!nativeOpaqueLineHandleQ[handle, session],
        err["E8", <|"Result" -> result,
          "Detail" -> "completed native checkpoint contains a malformed line handle"|>]];
      Lookup[handle, "line", Lookup[handle, "Line", None]],
    "endpoint",
      If[!nativeOpaqueEndpointHandleQ[handle, session],
        err["E8", <|"Result" -> result,
          "Detail" -> "completed native checkpoint contains a malformed endpoint handle"|>]];
      Lookup[handle, "endpoint", Lookup[handle, "Endpoint", None]]];
  If[Lookup[handle, "checkpoint_identity",
        Lookup[handle, "CheckpointIdentity", None]] =!= checkpoint,
    err["E8", <|"Result" -> result,
      "Detail" -> "observable and retained-result checkpoint identities differ"|>]];
  provenance = Lookup[handle, "provenance_identity",
    Lookup[handle, "ProvenanceIdentity", None]];
  If[!nativeNonemptyStringQ[provenance],
    err["E8", <|"Result" -> result,
      "Detail" -> "completed native checkpoint result has no exact provenance identity"|>]];
  <|"RequestIndex" -> result["RequestIndex"],
    "Operation" -> operation, "Identity" -> result["Identity"],
    "CheckpointIdentity" -> checkpoint, "Epsilon" -> result["Epsilon"],
    "Kind" -> kind, "Handle" -> token,
    "ProvenanceIdentity" -> provenance|>];

nativeObservableCheckpointManifestQ[manifest_] := Module[
  {keys = {"Schema", "Path", "CheckpointIdentity", "ManifestIdentity",
      "TransportArmMarches", "StateHandles", "Results"}, core, results,
   resultKeys = {"RequestIndex", "Operation", "Identity",
      "CheckpointIdentity", "Epsilon", "Kind", "Handle",
      "ProvenanceIdentity"}},
  If[!AssociationQ[manifest] || Sort[Keys[manifest]] =!= Sort[keys],
    Return[False, Module]];
  core = KeyDrop[manifest, "ManifestIdentity"];
  results = manifest["Results"];
  TrueQ[manifest["Schema"] ===
      "DiffExp2.NativeTransportObservableCheckpoint/v1"] &&
    StringQ[manifest["Path"]] && StringLength[manifest["Path"]] > 0 &&
    nativeNonemptyStringQ[manifest["CheckpointIdentity"]] &&
    manifest["ManifestIdentity"] === nativeCheckpointIdentity[
      "de2-native-observable-checkpoint-manifest-", core] &&
    IntegerQ[manifest["TransportArmMarches"]] &&
    manifest["TransportArmMarches"] >= 2 &&
    AssociationQ[manifest["StateHandles"]] &&
    Sort[Keys[manifest["StateHandles"]]] === {"lower", "upper"} &&
    AllTrue[Values[manifest["StateHandles"]], AssociationQ[#] &&
      Sort[Keys[#]] === Sort[{"Handle", "CheckpointIdentity",
        "ProvenanceIdentity"}] &&
      nativeNonemptyStringQ[#["Handle"]] &&
      StringStartsQ[#["Handle"], "transport:"] &&
      nativeNonemptyStringQ[#["CheckpointIdentity"]] &&
      nativeNonemptyStringQ[#["ProvenanceIdentity"]] &] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]], "Handle"]] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]],
      "CheckpointIdentity"]] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]],
      "ProvenanceIdentity"]] &&
    ListQ[results] && results =!= {} &&
    AllTrue[results, AssociationQ[#] &&
      Sort[Keys[#]] === Sort[resultKeys] &] &&
    Lookup[results, "RequestIndex"] === Range[0, Length[results] - 1] &&
    DuplicateFreeQ[Lookup[results, "Identity"]] &&
    DuplicateFreeQ[Lookup[results, "CheckpointIdentity"]] &&
    DuplicateFreeQ[Lookup[results, "Handle"]] &&
    DuplicateFreeQ[Lookup[results, "ProvenanceIdentity"]] &&
    AllTrue[results, Function[result,
      nativeNonemptyStringQ[result["Identity"]] &&
      nativeNonemptyStringQ[result["CheckpointIdentity"]] &&
      nativeNonemptyStringQ[result["Handle"]] &&
      nativeNonemptyStringQ[result["ProvenanceIdentity"]] &&
      AssociationQ[result["Epsilon"]] &&
      Sort[Keys[result["Epsilon"]]] ===
        Sort[{"Min", "Max", "RequiredCompleteMax"}] &&
      AllTrue[Lookup[result["Epsilon"],
          {"Min", "Max", "RequiredCompleteMax"}], IntegerQ] &&
      result["Epsilon", "Min"] <=
        result["Epsilon", "RequiredCompleteMax"] <=
        result["Epsilon", "Max"] &&
      Switch[result["Operation"],
        "integrate", result["Kind"] === "line" &&
          StringStartsQ[result["Handle"], "line:"],
        "limitLower" | "limitUpper", result["Kind"] === "endpoint" &&
          StringStartsQ[result["Handle"], "e:"],
        _, False]]]];

SaveNativeTransportObservableBatchCheckpoint[batch_Association,
    path_String, identity_String] := Module[
  {atlas, session, results, states, stateHandles, stats, saved, expanded,
   core},
  atlas = Lookup[batch, "Atlas", None];
  results = Lookup[batch, "Results", None];
  states = Lookup[batch, "States", None];
  session = If[AssociationQ[atlas], Lookup[atlas, "Session", None], None];
  If[Lookup[batch, "Type", None] =!=
        "DiffExp2NativeTransportObservableBatch" ||
      !nativeNonemptyStringQ[session] || !ListQ[results] || results === {} ||
      !AssociationQ[states] || Sort[Keys[states]] =!= {"lower", "upper"} ||
      !nativeNonemptyStringQ[identity],
    err["E8", <|"Detail" ->
      "native observable checkpoint save requires one completed nonempty two-arm batch and a nonempty identity"|>]];
  If[!nativeTransportStateHandleQ[states["lower"], session, "lower"] ||
      !nativeTransportStateHandleQ[states["upper"], session, "upper"],
    err["E8", <|"Detail" ->
      "native observable checkpoint save received malformed lower/upper transport states"|>]];
  stateHandles = AssociationMap[<|
      "Handle" -> Lookup[states[#], "transport_state",
        Lookup[states[#], "TransportState", None]],
      "CheckpointIdentity" -> Lookup[states[#], "checkpoint_identity",
        Lookup[states[#], "CheckpointIdentity", None]],
      "ProvenanceIdentity" -> Lookup[states[#], "provenance_identity",
        Lookup[states[#], "ProvenanceIdentity", None]]|> &,
    {"lower", "upper"}];
  stats = DiffExp2`CppBackend`RunRequest[<|"schema" -> 2,
    "op" -> "session.stats", "session" -> session|>];
  If[FailureQ[stats] || !AssociationQ[stats] ||
      Lookup[stats, "status", "error"] =!= "ok" ||
      !IntegerQ[Lookup[stats, "transport_arm_marches", None]],
    err["E5", <|"BackendFailure" -> stats,
      "Detail" -> "could not inspect native transport counters before checkpoint save"|>]];
  expanded = ExpandFileName[path];
  If[!DirectoryQ[DirectoryName[expanded]],
    Quiet[Check[CreateDirectory[DirectoryName[expanded],
      CreateIntermediateDirectories -> True], Null]]];
  If[!DirectoryQ[DirectoryName[expanded]],
    err["E5", <|"Path" -> expanded,
      "Detail" -> "native observable checkpoint directory could not be created"|>]];
  saved = DiffExp2`CppBackend`SavePersistentCheckpoint[
    session, expanded, identity];
  If[FailureQ[saved] || !AssociationQ[saved] ||
      Lookup[saved, "status", "error"] =!= "ok" ||
      Lookup[saved, "path", None] =!= expanded ||
      Lookup[saved, "checkpoint_identity", None] =!= identity ||
      !TrueQ[Lookup[saved, "atomic", False]],
    err["E5", <|"BackendFailure" -> saved,
      "Detail" -> "schema-2 native observable checkpoint save failed"|>]];
  core = <|"Schema" ->
      "DiffExp2.NativeTransportObservableCheckpoint/v1",
    "Path" -> expanded, "CheckpointIdentity" -> identity,
    "TransportArmMarches" -> stats["transport_arm_marches"],
    "StateHandles" -> stateHandles,
    "Results" -> (nativeObservableCheckpointResult[#, session] & /@
      results)|>;
  Append[core, "ManifestIdentity" -> nativeCheckpointIdentity[
    "de2-native-observable-checkpoint-manifest-", core]]];

RestoreNativeTransportObservableBatchCheckpoint[manifest_Association] :=
 Module[{restored, session = None, close, expectedLines, expectedEndpoints,
   expectedStates, stats, results, restoredHandles, restoredRecordMap,
   lineMap, endpointMap, stateMap, restoredResultIdentitiesQ,
   restoredStateIdentitiesQ},
  If[!nativeObservableCheckpointManifestQ[manifest],
    err["E8", <|"Detail" ->
      "native observable checkpoint manifest is malformed or self-inconsistent"|>]];
  If[!FileExistsQ[manifest["Path"]],
    err["E8", <|"Path" -> manifest["Path"],
      "Detail" -> "native observable checkpoint file does not exist"|>]];
  restored = DiffExp2`CppBackend`RestorePersistentCheckpoint[
    manifest["Path"], manifest["CheckpointIdentity"]];
  If[FailureQ[restored] || !AssociationQ[restored] ||
      Lookup[restored, "status", "error"] =!= "ok" ||
      !nativeNonemptyStringQ[Lookup[restored, "session", None]],
    err["E5", <|"BackendFailure" -> restored,
      "Detail" -> "schema-2 native observable checkpoint restore failed"|>]];
  session = restored["session"];
  close[] := Quiet[DiffExp2`CppBackend`ClosePersistentSession[session]];
  (* Lookup[{}, key] is Missing[KeyAbsent, key], not {}.  A completed FT
     level may legitimately contain only line integrals (or only endpoint
     limits), so preserve the empty visibility class explicitly. *)
  expectedLines = Lookup[
    Select[manifest["Results"], # ["Kind"] === "line" &], "Handle", {}];
  expectedEndpoints = Lookup[
    Select[manifest["Results"], # ["Kind"] === "endpoint" &], "Handle", {}];
  expectedStates = Lookup[Values[manifest["StateHandles"]], "Handle"];
  restoredHandles[collection_, key_String] := If[ListQ[collection],
    Map[If[AssociationQ[#], Lookup[#, key, None], #] &, collection], {}];
  restoredRecordMap[collection_, key_String] := If[ListQ[collection] &&
      AllTrue[collection, AssociationQ],
    Association@Map[Lookup[#, key, None] -> # &, collection], <||>];
  lineMap = restoredRecordMap[Lookup[restored, "line_results", {}],
    "line"];
  endpointMap = restoredRecordMap[Lookup[restored, "endpoints", {}],
    "endpoint"];
  stateMap = restoredRecordMap[Lookup[restored, "transport_states", {}],
    "transport_state"];
  restoredResultIdentitiesQ = AllTrue[manifest["Results"],
    Function[result, Module[{record = Lookup[
        If[result["Kind"] === "line", lineMap, endpointMap],
        result["Handle"], None]},
      AssociationQ[record] &&
        Lookup[record, "checkpoint_identity", None] ===
          result["CheckpointIdentity"] &&
        Lookup[record, "provenance_identity", None] ===
          result["ProvenanceIdentity"]]]];
  restoredStateIdentitiesQ = AllTrue[Values[manifest["StateHandles"]],
    Function[state, Module[{record = Lookup[stateMap,
        state["Handle"], None]},
      AssociationQ[record] &&
        Lookup[record, "checkpoint_identity", None] ===
          state["CheckpointIdentity"] &&
        Lookup[record, "provenance_identity", None] ===
          state["ProvenanceIdentity"]]]];
  If[Lookup[restored, "checkpoint_identity", None] =!=
        manifest["CheckpointIdentity"] ||
      !TrueQ[Lookup[restored, "replayed_wolfram_preprocessing", True] ===
        False] ||
      Sort[restoredHandles[Lookup[restored, "line_results", {}],
          "line"]] =!= Sort[expectedLines] ||
      Sort[restoredHandles[Lookup[restored, "endpoints", {}],
          "endpoint"]] =!= Sort[expectedEndpoints] ||
      Sort[restoredHandles[Lookup[restored, "transport_states", {}],
          "transport_state"]] =!=
        Sort[expectedStates] || !TrueQ[restoredResultIdentitiesQ] ||
      !TrueQ[restoredStateIdentitiesQ],
    close[];
    err["E5", <|
      "RestoredLineHandles" -> restoredHandles[
        Lookup[restored, "line_results", {}], "line"],
      "RestoredEndpointHandles" -> restoredHandles[
        Lookup[restored, "endpoints", {}], "endpoint"],
      "RestoredTransportStateHandles" -> restoredHandles[
        Lookup[restored, "transport_states", {}], "transport_state"],
      "Detail" -> "restored native checkpoint visibility differs from the completed observable manifest"|>]];
  stats = DiffExp2`CppBackend`RunRequest[<|"schema" -> 2,
    "op" -> "session.stats", "session" -> session|>];
  If[FailureQ[stats] || !AssociationQ[stats] ||
      Lookup[stats, "status", "error"] =!= "ok" ||
      Lookup[stats, "transport_arm_marches", None] =!=
        manifest["TransportArmMarches"],
    close[];
    err["E5", <|"BackendFailure" -> stats,
      "ExpectedTransportArmMarches" -> manifest["TransportArmMarches"],
      "Detail" -> "native checkpoint restore changed the transport-arm march counter"|>]];
  results = Map[Function[result,
    Join[KeyTake[result, {"RequestIndex", "Operation", "Identity",
        "CheckpointIdentity", "Epsilon"}],
      <|If[result["Kind"] === "line", "Line", "Endpoint"] ->
        <|"session" -> session, result["Kind"] -> result["Handle"],
          "request_index" -> result["RequestIndex"],
          "observable_identity" -> result["Identity"],
          "checkpoint_identity" -> result["CheckpointIdentity"],
          "provenance_identity" -> result["ProvenanceIdentity"]|>|>]],
    manifest["Results"]];
  <|"Type" -> "DiffExp2NativeTransportObservableBatch",
    "Atlas" -> None, "States" -> <||>, "Results" -> results,
    "NativeMarches" -> 0,
    "RestoredNativeMarches" -> manifest["TransportArmMarches"],
    "CompatibilityExports" -> 0, "RestoredCheckpoint" -> True,
    "RestoredSession" -> session, "CheckpointManifest" -> manifest|>];

ReleaseNativeTransportObservableBatch[batch_Association] := Module[
  {results, lines, endpoints, states, responses = {}, atlasResponse,
   failures, releaseOKQ, restoredSession, closed},
  If[Lookup[batch, "Type", None] =!=
      "DiffExp2NativeTransportObservableBatch",
    err["E8", <|"Type" -> Lookup[batch, "Type", None],
      "Detail" -> "native transport batch release requires one completed observable batch"|>]];
  If[TrueQ[Lookup[batch, "RestoredCheckpoint", False]],
    restoredSession = Lookup[batch, "RestoredSession", None];
    If[!nativeNonemptyStringQ[restoredSession],
      err["E8", <|"Detail" ->
        "restored native transport batch has no exact session token"|>]];
    closed = DiffExp2`CppBackend`ClosePersistentSession[restoredSession];
    Return[<|"Released" -> If[AssociationQ[closed] &&
          Lookup[closed, "status", "error"] === "ok", 1, 0],
      "Failures" -> If[AssociationQ[closed] &&
          Lookup[closed, "status", "error"] === "ok", {}, {closed}]|>,
      Module]];
  results = Lookup[batch, "Results", {}];
  lines = Cases[Lookup[results, "Line", None], _Association];
  endpoints = Cases[Lookup[results, "Endpoint", None], _Association];
  states = If[AssociationQ[Lookup[batch, "States", None]],
    Values[batch["States"]], {}];
  Scan[AppendTo[responses,
      Quiet[DiffExp2`CppBackend`ReleasePersistentLineIntegral[#]]] &,
    DeleteDuplicatesBy[lines, Lookup[#, "line", None] &]];
  Scan[AppendTo[responses,
      Quiet[DiffExp2`CppBackend`ReleasePersistentEndpoint[#]]] &,
    DeleteDuplicatesBy[endpoints, Lookup[#, "endpoint", None] &]];
  Scan[AppendTo[responses,
      Quiet[DiffExp2`CppBackend`ReleasePersistentTransportArm[#]]] &,
    DeleteDuplicatesBy[Select[states, AssociationQ],
      Lookup[#, "transport_state", None] &]];
  atlasResponse = ReleaseNativeRegularIndependentArms[batch["Atlas"]];
  releaseOKQ[response_] := AssociationQ[response] &&
    Lookup[response, "status", "error"] === "ok";
  failures = Select[responses, !TrueQ[releaseOKQ[#]] &];
  If[!AssociationQ[atlasResponse] ||
      Lookup[atlasResponse, "Failures", {"malformed"}] =!= {},
    AppendTo[failures, atlasResponse]];
  <|"Released" -> Count[responses, _?releaseOKQ] +
      Lookup[atlasResponse, "Released", 0],
    "Failures" -> failures|>];

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
