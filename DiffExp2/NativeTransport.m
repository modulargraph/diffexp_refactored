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
  "RunNativeTransportObservableBatch[atlas,observables,var] marches the retained lower/upper atlas exactly once, then contracts every ordered integrate, limitLower, and limitUpper observable without rematching. ObservableContractionChunkSize defaults to 1 so prepared rows are sent to the native backend in bounded-memory deterministic chunks. Each observable contains Operation, Identity, CheckpointIdentity, CoefficientVector, and Epsilon; integrate observables may additionally contain TailPolicy. Results are opaque retained line/endpoint handles and preserve request order.";
RunNativeTransportObservableBatchOwned::usage =
  "RunNativeTransportObservableBatchOwned[atlasSymbol,observables,var] is the ownership-taking batch entry point. It compacts atlasSymbol after preparing the observable rows, before allocating native transport states, so large Wolfram ChartSystems are no longer retained at the march peak.";
SaveNativeTransportObservableBatchCheckpoint::usage =
  "SaveNativeTransportObservableBatchCheckpoint[batch,path,identity] atomically saves one completed retained observable batch through the schema-2 native checkpoint protocol. The returned compact manifest binds the stable line, endpoint, and transport-state handles by SHA-256 provenance digests together with the pre-save transport-arm march counter.";
RestoreNativeTransportObservableBatchCheckpoint::usage =
  "RestoreNativeTransportObservableBatchCheckpoint[manifest] restores one completed retained observable batch without replaying chart preparation, matching, or transport. It validates either a compact v2 digest manifest or a legacy v1 full-provenance manifest and the transport-arm march counter before returning an exportable opaque batch.";
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

nativeStageTiming[fields___] := If[
  Environment["DE2_NATIVE_STAGE_TIMING"] === "1",
  Print["DE2 NATIVE STAGE ", fields, " t=", SessionTime[],
    " memory=", MemoryInUse[]]];

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

nativeAlgebraicRationalApprox[value_] := Module[
  {canonical, digits, approximation},
  canonical = Quiet[Check[RootReduce[value], $Failed]];
  If[canonical === $Failed || !exactRealAlgebraicQ[canonical],
    Return[$Failed, Module]];
  If[exactRationalQ[canonical], Return[canonical, Module]];
  digits = Min[300, Max[80, 2 cfg["WorkingPrecision"]]];
  approximation = Quiet[Check[
    Rationalize[N[canonical, digits], 0], $Failed]];
  If[exactRationalQ[approximation], approximation, $Failed]];

nativeCertifiedInexactRadiusFloor[value_] := Module[
  {precision = Precision[value], midpoint, exponent, precisionBits,
   uncertainty, lower, floor},
  If[Head[value] =!= Real || !TrueQ[0 < value < Infinity] ||
      !NumberQ[precision] || precision === Infinity || precision < 32,
    Return[$Failed, Module]];
  midpoint = Quiet[Check[Rationalize[value, 0], $Failed]];
  exponent = Quiet[Check[Last[MantissaExponent[Abs[value], 2]],
      $Failed]];
  If[!exactRationalQ[midpoint] || !IntegerQ[exponent],
    Return[$Failed, Module]];
  (* A Mathematica p-digit Real carries an uncertainty of roughly one last
     binary place.  Sixteen ulps is deliberately conservative and gives an
     exact interval certificate independent of later SetPrecision calls. *)
  precisionBits = Floor[precision*Log[2, 10]];
  uncertainty = 2^(exponent - precisionBits + 4);
  lower = midpoint - uncertainty;
  floor = Floor[2^$nativeScaleFloorBits*lower]/
    2^$nativeScaleFloorBits;
  If[!exactRationalQ[uncertainty] || !exactRationalQ[floor] ||
      !TrueQ[0 < floor <= lower < midpoint],
    Return[$Failed, Module]];
  <|"Floor" -> floor, "Midpoint" -> midpoint,
    "Uncertainty" -> uncertainty, "PrecisionDigits" -> precision,
    "PrecisionBits" -> precisionBits|>];

nativeScaleBridgePrerequisiteQ[chart_Association] := Module[
  {scale = Lookup[chart, "Scale", None],
   center = Lookup[chart, "Center", None],
   localRadius = Lookup[chart, "LocalRadius", None],
   radius = Lookup[chart, "Radius", None],
   matchRadius = Lookup[chart, "MatchRadius", None],
   roc = cfg["RadiusOfConvergence"], floor, radiusFloor},
  If[exactRationalQ[scale],
    If[exactRationalQ[radius],
      Return[exactRationalQ[center] && exactRationalQ[localRadius] &&
        exactRationalQ[matchRadius] && scale =!= 0 &&
        localRadius > 0 && radius > 0 && matchRadius > 0, Module]];
    radiusFloor = nativeInwardScaleFloor[radius];
    Return[exactRationalQ[center] && scale =!= 0 &&
      exactRealAlgebraicQ[localRadius] &&
      exactRealAlgebraicQ[radius] && exactRationalQ[matchRadius] &&
      radiusFloor =!= $Failed &&
      exactAlgebraicTruthQ[radius == Abs[scale]*localRadius] &&
      TrueQ[0 < matchRadius < radiusFloor], Module]];
  floor = nativeInwardScaleFloor[scale];
  BooleanQ[Lookup[chart, "Singular", Missing["Absent"]]] &&
    exactRationalQ[center] && exactRationalQ[localRadius] &&
    localRadius > 0 && exactRationalQ[roc] && roc > 0 &&
    floor =!= $Failed && exactRealAlgebraicQ[radius] &&
    exactRealAlgebraicQ[matchRadius] &&
    exactAlgebraicTruthQ[radius == scale*localRadius] &&
    exactAlgebraicTruthQ[matchRadius == scale*roc]];

bridgeNativeRegularChartScale[chart_Association, index_Integer] := Module[
  {scale = Lookup[chart, "Scale", None], center, localRadius, radius,
   matchRadius, roc = cfg["RadiusOfConvergence"], nativeScale,
   nativeRadius, nativeMatchRadius, certificate, radiusFloor,
   canonicalRadius, originalLocalRadius},
  If[exactRationalQ[scale],
    radius = Lookup[chart, "Radius", None];
    If[exactRationalQ[radius], Return[chart, Module]];
    center = Lookup[chart, "Center", None];
    localRadius = Lookup[chart, "LocalRadius", None];
    originalLocalRadius = localRadius;
    matchRadius = Lookup[chart, "MatchRadius", None];
    canonicalRadius = Quiet[Check[RootReduce[radius], $Failed]];
    If[canonicalRadius =!= $Failed &&
        exactRealAlgebraicQ[canonicalRadius],
      nativeRadius = nativeInwardScaleFloor[canonicalRadius];
      If[!exactRationalQ[center] || scale === 0 ||
          !exactRealAlgebraicQ[localRadius] ||
          !exactAlgebraicTruthQ[
            canonicalRadius == Abs[scale]*localRadius] ||
          nativeRadius === $Failed || !exactRationalQ[matchRadius] ||
          !TrueQ[0 < matchRadius < nativeRadius],
        err["E6", <|"ChartIndex" -> index, "Center" -> center,
          "Scale" -> scale, "Radius" -> radius,
          "MatchRadius" -> matchRadius,
          "Detail" -> "exact algebraic native radius has no strict positive inward rational disk containing its match disk"|>]];
      localRadius = Together[nativeRadius/Abs[scale]];
      certificate = <|
        "Schema" ->
          "diffexp2-native-inward-rational-algebraic-radius-v1",
        "BridgeKind" -> "ExactAlgebraicPhysicalRadius",
        "FloorBits" -> $nativeScaleFloorBits,
        "OriginalPhysicalRadius" -> canonicalRadius,
        "NativePhysicalRadius" -> nativeRadius,
        "OriginalLocalRadius" -> originalLocalRadius,
        "NativeLocalRadius" -> localRadius,
        "ScalePreserved" -> scale, "CenterPreserved" -> center,
        "MatchRadiusPreserved" -> matchRadius,
        "PrescriptionsPreserved" -> True,
        "ExactContainmentProved" -> True|>;
      Return[Join[chart, <|"Radius" -> nativeRadius,
        "LocalRadius" -> localRadius,
        "NativeRationalScaleBridge" -> certificate|>], Module]];
    radiusFloor = nativeCertifiedInexactRadiusFloor[radius];
    If[!exactRationalQ[center] || scale === 0 ||
        radiusFloor === $Failed || !exactRationalQ[matchRadius] ||
        !TrueQ[0 < matchRadius < radiusFloor["Floor"]],
      err["E6", <|"ChartIndex" -> index, "Center" -> center,
        "Scale" -> scale, "Radius" -> radius,
        "MatchRadius" -> matchRadius,
        "Detail" -> "inexact native radius has no certified positive inward rational disk containing its exact match disk"|>]];
    nativeRadius = radiusFloor["Floor"];
    localRadius = Together[nativeRadius/Abs[scale]];
    certificate = <|
      "Schema" -> "diffexp2-native-inward-rational-radius-v1",
      "FloorBits" -> $nativeScaleFloorBits,
      "OriginalRadiusMidpoint" -> radiusFloor["Midpoint"],
      "OriginalRadiusUncertainty" -> radiusFloor["Uncertainty"],
      "OriginalRadiusPrecisionDigits" -> radiusFloor["PrecisionDigits"],
      "OriginalRadiusPrecisionBits" -> radiusFloor["PrecisionBits"],
      "CertifiedOriginalRadiusLower" ->
        radiusFloor["Midpoint"] - radiusFloor["Uncertainty"],
      "NativePhysicalRadius" -> nativeRadius,
      "NativeLocalRadius" -> localRadius,
      "ScalePreserved" -> scale, "CenterPreserved" -> center,
      "PrescriptionsPreserved" -> True|>;
    Return[Join[chart, <|"Radius" -> nativeRadius,
      "LocalRadius" -> localRadius,
      "NativeRationalScaleBridge" -> certificate|>], Module]];
  center = Lookup[chart, "Center", None];
  localRadius = Lookup[chart, "LocalRadius", None];
  radius = Lookup[chart, "Radius", None];
  matchRadius = Lookup[chart, "MatchRadius", None];
  If[!BooleanQ[Lookup[chart, "Singular", Missing["Absent"]]],
    err["E6", <|"ChartIndex" -> index, "Center" -> center,
      "Scale" -> scale,
      "Detail" -> "nonrational native scale bridging requires an explicit regular/singular chart classification"|>]];
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
    "Singular" -> chart["Singular"],
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
    dimension_Integer] := Module[{active, shifts},
  If[Length[cvec] =!= dimension,
    err["E8", <|"CoefficientCount" -> Length[cvec],
      "Dimension" -> dimension,
      "Detail" -> "native regular-arm integration requires one coefficient per physical component"|>]];
  (* A structurally exact zero has no valuation and must not force the
     minimum back to eps^0.  Inexact zeros remain active so an uncertain
     coefficient cannot silently relax the honest source window. *)
  active = Select[cvec, !TrueQ[FreeQ[#, _?InexactNumberQ] &&
        PossibleZeroQ[Together[#]]] &];
  shifts = nativeRationalEpsilonShift[#, physicalVar] & /@ active;
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
   normalization, common, lower2, upper2,
   roc = cfg["RadiusOfConvergence"]},
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
  If[!AllTrue[{lowerAnchor["Radius"], upperAnchor["Radius"]},
        exactRealAlgebraicQ] ||
      !AllTrue[{lowerAnchor["Radius"], upperAnchor["Radius"]},
        exactAlgebraicTruthQ[# > 0] &] ||
      !AllTrue[{lowerAnchor["MatchRadius"],
          upperAnchor["MatchRadius"], roc}, exactRationalQ],
    err["E6", <|"LowerAnchor" -> lowerAnchor,
      "UpperAnchor" -> upperAnchor,
      "Detail" -> "shared-anchor normalization requires positive exact real algebraic radii and rational match/scale data"|>]];
  (* SegmentLine caps anchor geometry by each arm length.  The larger of the
     two conservative caps is still no larger than the true singularity
     radius and gives one chart capable of serving both arms. *)
  commonRadius = Max[lowerAnchor["Radius"], upperAnchor["Radius"]];
  commonMatchRadius = Max[lowerAnchor["MatchRadius"],
    upperAnchor["MatchRadius"]];
  scale = Together[commonMatchRadius/roc];
  normalization = <|
    "Schema" -> "diffexp2-native-shared-anchor-normalization-v1",
    "IncomingLowerPhysicalRadius" -> lowerAnchor["Radius"],
    "IncomingUpperPhysicalRadius" -> upperAnchor["Radius"],
    "IncomingLowerMatchRadius" -> lowerAnchor["MatchRadius"],
    "IncomingUpperMatchRadius" -> upperAnchor["MatchRadius"],
    "SelectedPhysicalRadius" -> commonRadius,
    "SelectedMatchRadius" -> commonMatchRadius|>;
  common = Join[lowerAnchor, <|"Radius" -> commonRadius,
    "MatchRadius" -> commonMatchRadius, "Scale" -> scale,
    "LocalRadius" -> Together[commonRadius/scale],
    "SharedAnchorNormalization" -> normalization|>];
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
  Map[Function[pair, Module[{re = pair[[1]], h = pair[[2]], flags,
      nativeRe, nativeH},
    nativeRe = nativeAlgebraicRationalApprox[re];
    nativeH = nativeAlgebraicRationalApprox[h];
    If[nativeRe === $Failed || nativeH === $Failed ||
        !TrueQ[nativeH > 0],
      err["E6", <|"Projection" -> pair,
        "Detail" -> "an on-arm complex projection is not representable by the current rational native path protocol"|>]];
    flags = {AnyTrue[projected, sameExactQ[#, re - h] &],
      AnyTrue[projected, sameExactQ[#, re] &],
      AnyTrue[projected, sameExactQ[#, re + h] &]};
    <|"source_identity" -> ("complex-projection:" <>
        IntegerString[Hash[{re, h}, "SHA256"], 16, 64]),
      "real_part_exact" ->
        exactRationalString[nativeRe, "projection real part"],
      "imaginary_magnitude_exact" ->
        exactRationalString[nativeH, "projection imaginary magnitude"],
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
  If[!AllTrue[relevant, Module[{re, h},
      re = nativeAlgebraicRationalApprox[#[[1]]];
      h = nativeAlgebraicRationalApprox[#[[2]]];
      re =!= $Failed && h =!= $Failed && TrueQ[h > 0]] &],
    Return[False, Module]];
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
          Lookup[upperAnchor, "Radius", None]},
        exactRealAlgebraicQ] ||
      !AllTrue[{Lookup[lowerAnchor, "Radius", None],
          Lookup[upperAnchor, "Radius", None]},
        exactAlgebraicTruthQ[# > 0] &] ||
      !AllTrue[{Lookup[lowerAnchor, "MatchRadius", None],
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
  lowerOverlay = bridgeNativeRegularPlanScales[lowerOverlay];
  upperOverlay = bridgeNativeRegularPlanScales[upperOverlay];
  nativeTopologyProtocolQ[lower] && nativeTopologyProtocolQ[upper] &&
    nativeBridgedPlanPreflightQ[lowerOverlay] &&
    nativeBridgedPlanPreflightQ[upperOverlay]
  ], False]];
NativeRegularIndependentArmPlansSupportedQ[___] := False;

nativeBasisOwner[basis_Association] := Module[{owner},
  owner = Lookup[basis, "NativeSCC",
    Lookup[basis, "NativeChart", Lookup[basis, "SCC", None]]];
  If[!StringQ[owner] || StringLength[owner] == 0,
    err["E6", <|"Basis" -> KeyTake[basis,
        {"Type", "Session", "NativeSCC", "NativeChart"}],
      "Detail" -> "retained receiving basis exposes no native chart/SCC owner"|>]];
  owner];

SetAttributes[nativeCatchDE2, HoldFirst];
nativeCatchDE2[expression_] := Catch[expression, "DiffExp2Error"];

(* An Acb regular-singular basis is first attempted as a capability probe.
   CASE-P and exact multi-tag sources deliberately fail that probe before the
   Rational shadow is selected.  DE2Error prints before throwing, so a plain
   Catch makes a successfully recovered probe look like a terminal error.
   Buffer Print only across this narrow probe.  The caller discards precisely
   the expected DE2Error record and replays every other progress/debug print;
   unrecognized failures replay the complete buffer before remaining loud. *)
SetAttributes[nativeCatchDE2Buffered, HoldFirst];
nativeCatchDE2Buffered[expression_] := Module[{records = {}, result},
  result = Block[{Print = Function[Null,
      AppendTo[records, HoldComplete[##]], HoldAllComplete]},
    nativeCatchDE2[expression]];
  {result, records}];

nativeReplayPrintRecord[HoldComplete[args___]] := Print[args];
nativeReplayPrintRecords[records_List] :=
  Scan[nativeReplayPrintRecord, records];

$nativeAcbRationalShadowFailureNeedles = {
  "native Acb regular-singular SCC execution rejects exact CASE-P collisions",
  "native Acb regular-singular SCC execution requires the exact Rational shadow for a multi-tag gauge source",
  "requires the exact Rational shadow"};

nativeAcbRationalShadowTriggerTextQ[text_String] := AnyTrue[
  $nativeAcbRationalShadowFailureNeedles, StringContainsQ[text, #] &];

nativeAcbShadowTriggerPrintRecordQ[record_] :=
  MatchQ[record,
    HoldComplete["DiffExp2 error ", _String, ": ", _String]] &&
  !FreeQ[record,
    text_String /; nativeAcbRationalShadowTriggerTextQ[text], Infinity];

nativeAcbCasePFailureQ[failure_] := Module[{},
  If[!FailureQ[failure], Return[False, Module]];
  !FreeQ[Unevaluated[failure],
    text_String /; nativeAcbRationalShadowTriggerTextQ[text], Infinity]];

nativeShadowColumn[source_Association, targetSCC_Association,
    shadowIdentity_String, system_Association] := Module[
  {checkpointIdentity, response, provenance, result},
  checkpointIdentity = nativeCheckpointIdentity[
    "de2-native-acb-rational-shadow-", {
      source["CheckpointIdentity"], shadowIdentity,
      targetSCC["Session"], targetSCC["SCC"]}];
  response = DiffExp2`CppBackend`SpecializePersistentRationalSCCColumn[
    source, targetSCC, shadowIdentity, checkpointIdentity];
  If[FailureQ[response] || !AssociationQ[response] ||
      Lookup[response, "status", "error"] =!= "ok" ||
      Lookup[response, "session", None] =!= targetSCC["Session"] ||
      Lookup[response, "scc", None] =!= targetSCC["SCC"] ||
      Lookup[response, "json_coefficients", None] =!= 0 ||
      !StringQ[Lookup[response, "local", None]],
    If[AssociationQ[response] &&
        Lookup[response, "status", "error"] === "ok" &&
        StringQ[Lookup[response, "local", None]],
      Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[response]]];
    err["E5", <|"BackendFailure" -> response,
      "Detail" -> "exact Rational CASE-P shadow could not be specialized into its paired Acb SCC"|>]];
  provenance = Lookup[response, "column_provenance", None];
  If[!AssociationQ[provenance] ||
      Lookup[provenance, "basis_index", None] =!=
        source["BasisIndex"] - 1 ||
      !StringQ[Lookup[provenance, "exact_column_identity", None]],
    Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[response]];
    err["E6", <|"BackendResponse" -> response,
      "Detail" -> "specialized CASE-P shadow lost its exact physical-column provenance"|>]];
  result = <|"Type" -> "DiffExp2NativeSCCBasisColumn",
    "Session" -> response["session"], "Local" -> response["local"],
    "NativeSCC" -> response["scc"],
    "NativeChart" -> response["chart"],
    "SeedBlock" -> source["SeedBlock"],
    "BasisIndex" -> source["BasisIndex"],
    "Chart" -> <|"Center" -> system["Center"],
      "ChartMap" -> system["ChartMap"], "Radius" -> system["Radius"],
      "Prescriptions" -> system["Prescriptions"]|>,
    "EpsWindow" -> <|"Min" -> response["epsilon_min"],
      "CompleteMax" -> response["epsilon_max"]|>,
    "TWindow" -> <|"CompleteMax" ->
      response["taylor_complete_max"]|>,
    "CheckpointIdentity" -> checkpointIdentity,
    "ColumnProvenance" -> provenance,
    "NativeSummary" -> KeyDrop[response,
      {"status", "session", "local", "scc", "chart", "metadata",
       "column_provenance"}]|>;
  If[KeyExistsQ[source, "SeedLocalComponent"],
    Append[result, "SeedLocalComponent" -> source["SeedLocalComponent"]],
    result]];

nativeRationalShadowBasis[system_Association, req_Association, threads_,
    targetSCC_Association] := Module[
  {targetStats, shadowIdentity, rationalBasis, imported = {}, cleanup,
   result},
  targetStats = DiffExp2`CppBackend`PersistentSCCStatistics[targetSCC];
  shadowIdentity = If[AssociationQ[targetStats],
    Lookup[targetStats, "rational_shadow_identity", None], None];
  If[!StringQ[shadowIdentity] || StringLength[shadowIdentity] == 0,
    err["E6", <|"BackendStatistics" -> targetStats,
      "Detail" -> "Acb CASE-P SCC exposes no domain-independent Rational-shadow identity"|>]];
  (* The exact-domain shadow is a second live composite owner, distinct from
     targetSCC.  Reserve it at the point where the producer certificate says
     it is actually required; this also covers deferred streamed bases after
     the per-arm owner-preparation scope has ended. *)
  rationalBasis =
    DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[1,
      Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
        nativeCatchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
          system, req, threads]]]];
  If[FailureQ[rationalBasis] || !AssociationQ[rationalBasis] ||
      Lookup[rationalBasis, "Type", None] =!= "DiffExp2NativeSCCBasis" ||
      Lookup[rationalBasis, "Dimension", None] =!= system["SystemSize"],
    If[FailureQ[rationalBasis],
      Throw[rationalBasis, "DiffExp2Error"],
      err["E6", <|"RationalBasis" -> rationalBasis,
        "Detail" -> "exact CASE-P shadow did not produce one complete Rational SCC basis"|>]]];
  cleanup[] := Module[{},
    Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
      imported];
    Quiet[DiffExp2`CppBackend`ClosePersistentSession[rationalBasis]];
    Null];
  result = Catch[
    Do[AppendTo[imported,
      nativeShadowColumn[column, targetSCC, shadowIdentity, system]],
      {column, rationalBasis["Columns"]}];
    If[Lookup[imported, "BasisIndex", {}] =!=
        Range[system["SystemSize"]],
      err["E6", <|"BasisIndices" -> Lookup[imported, "BasisIndex", {}],
        "Detail" -> "specialized CASE-P shadow basis is not in complete physical order"|>]];
    <|"Type" -> "DiffExp2NativeSCCBasis",
      "Session" -> targetSCC["Session"],
      "NativeSCC" -> targetSCC["SCC"],
      "Columns" -> imported, "Dimension" -> system["SystemSize"],
      "Chart" -> <|"Center" -> system["Center"],
        "ChartMap" -> system["ChartMap"], "Radius" -> system["Radius"],
        "Prescriptions" -> system["Prescriptions"]|>,
      "EpsWindow" -> req["EpsWindow"],
      "TWindow" -> <|"CompleteMax" -> req["TOrder"]|>,
      "NativeSummary" -> <|
        "specialization_capability" ->
          "exact-rational-shadow-to-acb-local-v1",
        "rational_shadow_identity" -> shadowIdentity|>|>,
    "DiffExp2Error", Function[{failure, tag}, cleanup[];
      Throw[failure, tag]]];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[rationalBasis]];
  result];

nativeReceivingBasis[system_Association, req_Association, threads_] := Module[
  {regular = TrueQ[Lookup[
      Lookup[system, "IndicialData", <||>], "Regular", False]],
   sequence, components, built, expectedTypes, targetSCC, attempt,
   probePrints, shadowDecision},
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
    targetSCC = DiffExp2`Solve`PrepareNativeSCCComposite[system, req];
    shadowDecision =
      DiffExp2`Solve`Private`sccNativeCachedRationalShadowDecision[
        system, req, targetSCC];
    If[TrueQ[shadowDecision["RequiresRationalShadow"]],
      (* The collision-bound producer certificate selects the exact route for
         a CASE-P schedule or identity-frame polar source.  Build the shadow
         directly: no speculative Acb columns are executed or discarded. *)
      built = nativeRationalShadowBasis[
        system, req, threads, targetSCC];
      built = Join[built, <|"NativeSummary" -> Join[
        Lookup[built, "NativeSummary", <||>], <|
          "selection_capability" ->
            "producer-certified-proactive-rational-shadow-v1",
          "selection_certificate" ->
            shadowDecision["Certificate"]|>]|>],
      {attempt, probePrints} = nativeCatchDE2Buffered[
        DiffExp2`Solve`SolveNativeSCCBasis[system, req, threads]];
      built = Which[
        !FailureQ[attempt],
          nativeReplayPrintRecords[probePrints]; attempt,
        nativeAcbCasePFailureQ[attempt],
          nativeReplayPrintRecords[
            Select[probePrints, !nativeAcbShadowTriggerPrintRecordQ[#] &]];
          nativeRationalShadowBasis[system, req, threads, targetSCC],
        True,
          nativeReplayPrintRecords[probePrints];
          Throw[attempt, "DiffExp2Error"]]];
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

(* Keep this predicate identical to the owner/basis dispatch below.  A
   singular receiving chart always owns an SCC composite; a regular chart
   does so exactly when its certified condensation has more than one block.
   Monolithic regular charts instead retain one prepare-only physical chart
   owner and its immutable value-run prototype. *)
nativeReceivingSystemUsesSCCCompositeQ[system_Association] := Module[
  {regular = TrueQ[Lookup[
      Lookup[system, "IndicialData", <||>], "Regular", False]],
   sequence, components},
  If[!regular, Return[True, Module]];
  sequence = Lookup[system, "IntegrationSequence", None];
  components = If[AssociationQ[sequence],
    Lookup[sequence, "Components", {}], {}];
  ListQ[components] && Length[components] > 1];

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

$nativeObservableCheckpointSchemaV1 =
  "DiffExp2.NativeTransportObservableCheckpoint/v1";
$nativeObservableCheckpointSchemaV2 =
  "DiffExp2.NativeTransportObservableCheckpoint/v2";

(* Native provenance identities are deliberately complete, recursively
   auditable JSON records.  They belong in the atomic native checkpoint, but
   copying them into the Wolfram MX sidecar can duplicate tens of megabytes
   per result.  The compact manifest binds each restored record by SHA-256;
   the full string is recovered from, and validated against, the native file
   before an opaque handle is reconstructed. *)
nativeProvenanceSHA256[identity_String] :=
  IntegerString[Hash[identity, "SHA256"], 16, 64];

nativeProvenanceSHA256Q[value_] := StringQ[value] &&
  StringMatchQ[value, RegularExpression["[0-9a-f]{64}"]];

nativeCheckpointProvenanceKey[schema_] := Switch[schema,
  $nativeObservableCheckpointSchemaV1, "ProvenanceIdentity",
  $nativeObservableCheckpointSchemaV2, "ProvenanceSHA256",
  _, None];

nativeCheckpointProvenanceReference[identity_String, schema_] :=
  Switch[schema,
    $nativeObservableCheckpointSchemaV1, identity,
    $nativeObservableCheckpointSchemaV2, nativeProvenanceSHA256[identity],
    _, None];

nativeCheckpointProvenanceReferenceQ[value_, schema_] := Switch[schema,
  $nativeObservableCheckpointSchemaV1, nativeNonemptyStringQ[value],
  $nativeObservableCheckpointSchemaV2, nativeProvenanceSHA256Q[value],
  _, False];

nativeCheckpointProvenanceMatchesQ[identity_, reference_, schema_] :=
  nativeNonemptyStringQ[identity] &&
    nativeCheckpointProvenanceReference[identity, schema] === reference;

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

nativeArmRowRecipes[atlas_Association, data_Association] := Module[
  {systems, bases, shapes, charts},
  systems = data["ChartSystems"];
  bases = Rest[data["Bases"]];
  shapes = Prepend[
    If[AllTrue[bases, AssociationQ],
      nativeBasisEnvelopeShape /@ bases,
      ConstantArray[<|"EpsWindow" -> atlas["Request", "EpsWindow"],
        "TWindow" -> <|"CompleteMax" -> atlas["Request", "TOrder"]|>,
        "Dimension" -> atlas["Dimension"]|>, Length[bases]]],
    nativePreparationShape[atlas["Anchor"], atlas["Dimension"]]];
  If[Length[systems] =!= Length[shapes],
    err["E6", <|"ChartCount" -> Length[systems],
      "ShapeCount" -> Length[shapes],
      "Detail" -> "native arm row envelopes do not reproduce its chart chain"|>]];
  (* PrepareNativeRationalRow consumes only this small geometric chart
     projection.  Capturing it before atlas compaction lets observable rows
     be constructed one tile at a time after the native march without
     retaining the full chart systems or solved bases. *)
  charts = KeyTake[#, {"SystemSize", "ChartMap", "ChartVar", "Center",
        "Radius", "Prescriptions"}] & /@ systems;
  MapThread[<|"Chart" -> #1, "Shape" -> #2|> &, {charts, shapes}]];

nativePrepareArmRecipeRow[recipe_Association, cvec_List,
    var_Symbol, domain_String] := Module[{row},
  If[Sort[Keys[recipe]] =!= {"Chart", "Shape"} ||
      !AssociationQ[recipe["Chart"]] ||
      !AssociationQ[recipe["Shape"]],
    err["E6", <|"Detail" ->
      "lazy native integrand-row recipe is malformed"|>]];
  row = DiffExp2`Solve`PrepareNativeRationalRow[
    recipe["Chart"], recipe["Shape"], cvec, var,
    <|"domain" -> domain, "symbols" -> {}|>];
  If[!AssociationQ[row],
    err["E5", <|"Detail" ->
      "lazy native integrand-row preparation returned a malformed row"|>]];
  row];

nativeDropRationalMultiplierPreparationCache[] := Module[{},
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
  Null];

nativePreparedArmRows[atlas_Association, data_Association, cvec_List,
    var_Symbol] := Module[{recipes, rows},
  recipes = nativeArmRowRecipes[atlas, data];
  rows = nativePrepareArmRecipeRow[#, cvec, var, atlas["Domain"]] & /@
    recipes;
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
  "Integrands" -> Automatic, "TargetCompleteMax" -> Automatic,
  "RequiredTargetCompleteMax" -> Automatic,
  "DeferReceivingBases" -> False};

PrepareNativeRegularIndependentArms[sys_Association, boundary_,
    lowerPlan_Association, upperPlan_Association, OptionsPattern[]] := Module[
  {plans, lower, upper, dimension = Length[sys["Matrix"]], values,
   epsMin, epsMax, req, anchorSystem, anchor = None, prepareArm, lowerData,
   upperData, sessions, anchorOwner, lowerOwners, upperOwners, planIdentity,
   nativePlan = None, geometryAudit = None, sessionInfo, sessionStats,
   domain, integrand,
   preparedShift, halo, targetMax, targetOption =
     OptionValue["TargetCompleteMax"], requiredTargetMax,
   requiredTargetOption = OptionValue["RequiredTargetCompleteMax"],
   availableMax,
   minimumSolveMax, matchEpsilonPadding,
   cleanup, output, preparedBases = {}, preparedOwners = {},
   containsSingularReceivingCharts = False,
   integrands = OptionValue["Integrands"],
   threads = OptionValue["Threads"],
   deferReceivingBases = TrueQ[OptionValue["DeferReceivingBases"]]},
  If[targetOption =!= Automatic &&
      (!IntegerQ[targetOption] || targetOption < 0),
    err["E8", <|"TargetCompleteMax" -> targetOption,
      "Detail" -> "TargetCompleteMax must be Automatic or a nonnegative integer"|>]];
  targetMax = If[targetOption === Automatic,
    cfg["EpsilonOrder"], targetOption];
  requiredTargetMax = If[requiredTargetOption === Automatic,
    targetMax, requiredTargetOption];
  If[!IntegerQ[requiredTargetMax] || requiredTargetMax < 0 ||
      requiredTargetMax > targetMax,
    err["E8", <|"TargetCompleteMax" -> targetMax,
      "RequiredTargetCompleteMax" -> requiredTargetOption,
      "Detail" -> "RequiredTargetCompleteMax must be an integer between zero and TargetCompleteMax"|>]];
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
  minimumSolveMax = requiredTargetMax + halo;
  If[availableMax < minimumSolveMax,
    err["E6", <|"BoundaryCompleteMax" -> availableMax,
      "RequestedCompleteMax" -> targetMax,
      "IntegrandEpsilonShift" -> preparedShift,
      "RequiredSolveCompleteMax" -> minimumSolveMax,
      "Detail" -> "native regular-arm boundary data does not contain the epsilon halo required by its integrand pole"|>]];
  (* Retain the already available boundary reservoir as internal match work.
     Singularly perturbed receiving bases can have a positive determinant
     valuation even after exact negative-power saturation.  That valuation
     consumes solve coefficients, but must not inflate the public or residual
     target: those remain bound to TargetCompleteMax below. *)
  epsMax = availableMax;
  matchEpsilonPadding = epsMax - minimumSolveMax;
  req = <|"EpsWindow" -> <|"Min" -> epsMin,
      "CompleteMax" -> epsMax|>,
    (* CompleteMax is the private reservoir available to a match.  The
       receiving equation only has to deliver the state edge below; keeping
       those two contracts separate prevents a wider boundary reservoir from
       becoming a recursive public basis demand. *)
    "RequiredCompleteMax" -> minimumSolveMax,
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
    Scan[Function[record,
      If[StringQ[Lookup[record, "SCC", None]],
        Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[record]]]],
      preparedOwners];
    Null];
  output = Catch[
  nativeStageTiming["anchor-prepare-start"];
  anchorSystem = DiffExp2`Solve`PrepareChart[sys, First[lower["Charts"]]];
  nativeStageTiming["anchor-prepare-done"];
  If[!TrueQ[Lookup[anchorSystem["IndicialData"], "Regular", False]],
    err["E8", <|"Center" -> anchorSystem["Center"],
      "Detail" -> "native independent-arm anchor is not regular"|>]];
  anchor = DiffExp2`Solve`SolveNativeValueRegular[
    anchorSystem, req, values];
  nativeStageTiming["anchor-solve-done boundaryMax=", availableMax,
    " request=", req["EpsWindow"], " anchor=",
    Lookup[anchor, "EpsWindow", None]];
  sessionInfo = DiffExp2`CppBackend`PersistentSessionInformation[];
  sessionStats = Lookup[sessionInfo, anchor["Session"], None];
  domain = If[AssociationQ[sessionStats],
    Lookup[sessionStats, "domain", None], None];
  If[!MemberQ[{"acb", "rational"}, domain],
    err["E5", <|"Domain" -> domain,
      "Detail" -> "native independent-arm atlas requires Acb or Rational retained locals"|>]];
  anchorOwner = anchor["NativeChart"];
  prepareArm[plan_Association] := Module[
    {systems, bases, built, prepared, kinds, ownerRecords, owners,
     sccOwnerCount},
    systems = Prepend[
      MapIndexed[Function[{chart, index},
        nativeStageTiming["chart-prepare-start index=", First[index]];
        prepared = DiffExp2`Solve`PrepareChart[sys, chart];
        nativeStageTiming["chart-prepare-done index=", First[index]];
        prepared], Rest[plan["Charts"]]],
      anchorSystem];
    kinds = Prepend[Map[
      If[TrueQ[Lookup[Lookup[#, "IndicialData", <||>],
          "Regular", False]], "Regular", "SingularSCC"] &,
      Rest[systems]], "Anchor"];
    If[MemberQ[kinds, "SingularSCC"],
      containsSingularReceivingCharts = True];
    sccOwnerCount = Count[Rest[systems],
      system_ /; nativeReceivingSystemUsesSCCCompositeQ[system]];
    If[deferReceivingBases,
      ownerRecords =
        DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[
          sccOwnerCount,
          MapIndexed[Function[{system, index},
            nativeStageTiming["owner-prepare-start index=", First[index]];
            built = If[TrueQ[Lookup[
                  Lookup[system, "IndicialData", <||>],
                  "Regular", False]],
              DiffExp2`Solve`PrepareNativeRegularBasisOwner[system, req],
              DiffExp2`Solve`PrepareNativeSCCComposite[system, req]];
            AppendTo[preparedOwners, built];
            nativeStageTiming["owner-prepare-done index=", First[index],
              " center=", InputForm[Lookup[system, "Center", None]],
              " regular=", TrueQ[Lookup[
                Lookup[system, "IndicialData", <||>],
                "Regular", False]],
              " ownerType=", Lookup[built, "Type", None]];
            built], Rest[systems]]];
      owners = Prepend[nativeBasisOwner /@ ownerRecords, anchorOwner];
      bases = ConstantArray[None, Length[systems]],
      ownerRecords = {};
      bases = Prepend[
        DiffExp2`Solve`WithNativeSCCCompositeCacheReservation[
          (* Every eager target may select one exact Rational shadow.  The
             native shadow session is closed immediately, but its Solve-side
             descriptor remains as a deliberately stale cache record until
             the next cache drop, so fund both records for the whole map. *)
          2 sccOwnerCount,
          MapIndexed[Function[{system, index},
            nativeStageTiming["basis-solve-start index=", First[index]];
            built = nativeReceivingBasis[system, req, threads];
            AppendTo[preparedBases, built];
            nativeStageTiming["basis-solve-done index=", First[index]];
            built], Rest[systems]]], None];
      owners = Prepend[nativeBasisOwner /@ Rest[bases], anchorOwner]];
    <|"Plan" -> plan, "ChartSystems" -> systems,
      "Bases" -> bases, "BasisKinds" -> kinds,
      "Owners" -> owners, "OwnerRecords" -> ownerRecords|>];
  lowerData = prepareArm[lower];
  upperData = prepareArm[upper];
  sessions = DeleteDuplicates@Join[{anchor["Session"]},
    Cases[Join[Rest[lowerData["Bases"]], Rest[upperData["Bases"]],
      lowerData["OwnerRecords"], upperData["OwnerRecords"]],
      b_Association :> b["Session"]]];
  If[Length[sessions] =!= 1,
    err["E6", <|"Sessions" -> sessions,
      "Detail" -> "prepared native atlas was split across solver sessions"|>]];
  lowerOwners = lowerData["Owners"];
  upperOwners = upperData["Owners"];
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
     above preserves the center and proves that the resulting physical disk
     is a strict exact subset of the original algebraic disk: either Scale is
     floored with LocalRadius fixed, or Radius/LocalRadius are floored with
     Scale fixed.  This composes the two proofs without a duplicate Wolfram
     geometry walk. *)
  geometryAudit = <|
    "Schema" -> "diffexp2-native-exact-geometry-proof-v1",
    "BackendExactPlanValidated" -> True,
    "NativeDisksContainedInOriginal" -> True,
    "HandoffAuthority" -> "CppExactTilePlanner",
    "LowerBridgeCount" -> Length[Lookup[lower,
      "NativeRationalScaleBridges", {}]],
    "UpperBridgeCount" -> Length[Lookup[upper,
      "NativeRationalScaleBridges", {}]]|>;
  (* The tile plan now strongly owns every equation owner.  Regular owners
     were prepare-only, so there is no disposable local slab to release. *)
  lowerData = KeyDrop[lowerData, "OwnerRecords"];
  upperData = KeyDrop[upperData, "OwnerRecords"];
  <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
    "Session" -> First[sessions], "Domain" -> domain,
    "Dimension" -> dimension,
    "Request" -> req, "Anchor" -> anchor, "Plan" -> nativePlan,
    "Lower" -> lowerData, "Upper" -> upperData,
    "NativeGeometryAudit" -> geometryAudit,
    "ContainsSingularReceivingCharts" ->
      containsSingularReceivingCharts,
    "DeferredReceivingBases" -> deferReceivingBases,
    "TargetCompleteMax" -> targetMax,
    "RequiredTargetCompleteMax" -> requiredTargetMax,
    "MatchEpsilonPadding" -> matchEpsilonPadding,
    "Threads" -> threads,
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
   tailPolicy, divergentCancellation, integrateKeys},
  If[!AssociationQ[raw],
    err["E8", <|"Observable" -> raw,
      "Detail" -> "native transport batch observables must be associations"|>]];
  operation = Lookup[raw, "Operation", None];
  keys = Sort[Keys[raw]];
  integrateKeys = Join[baseKeys,
    Select[{"TailPolicy", "DivergentCancellation"},
      KeyExistsQ[raw, #] &]];
  If[!MemberQ[{"integrate", "limitLower", "limitUpper"}, operation] ||
      (operation === "integrate" && keys =!= Sort[integrateKeys]) ||
      (operation =!= "integrate" && keys =!= Sort[baseKeys]),
    err["E8", <|"Operation" -> operation, "Keys" -> Keys[raw],
      "Detail" -> "native transport observables require exactly Operation, Identity, CheckpointIdentity, CoefficientVector, and Epsilon; only integrate accepts optional TailPolicy and DivergentCancellation"|>]];
  epsilon = raw["Epsilon"];
  tailPolicy = Lookup[raw, "TailPolicy", "stored"];
  divergentCancellation = Lookup[raw, "DivergentCancellation", None];
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
        (!MemberQ[{"stored", "attempt", "require"}, tailPolicy] ||
         (divergentCancellation =!= None &&
          (!AssociationQ[divergentCancellation] ||
           Sort[Keys[divergentCancellation]] =!=
             Sort[{"Mode", "RelativeTolerance", "Provenance"}] ||
           Lookup[divergentCancellation, "Mode", None] =!=
             "bounded-relative-acb" ||
           !nativeNonemptyStringQ[Lookup[divergentCancellation,
             "RelativeTolerance", None]] ||
           !nativeNonemptyStringQ[Lookup[divergentCancellation,
             "Provenance", None]])))),
    err["E8", <|"Identity" -> Lookup[raw, "Identity", None],
      "Operation" -> operation,
      "Detail" -> "native transport observable identities, coefficient dimension, epsilon window, or tail policy are malformed"|>]];
  Join[<|"Operation" -> operation, "Identity" -> raw["Identity"],
    "CheckpointIdentity" -> raw["CheckpointIdentity"],
    "CoefficientVector" -> raw["CoefficientVector"],
    "Epsilon" -> epsilon, "TailPolicy" -> tailPolicy|>,
    If[divergentCancellation === None, <||>,
      <|"DivergentCancellation" -> divergentCancellation|>]]];

nativePrepareBatchObservable[observable_Association, var_Symbol,
    dimension_Integer] := Module[
  {shift, coefficients = observable["CoefficientVector"]},
  (* This is deliberately a metadata-only prepass.  Eagerly preparing every
     chart row here used to retain both the rows and SectorSeries multiplier
     cache entries throughout the march, defeating native tile streaming. *)
  shift = nativeIntegrandMinimumShift[coefficients, var, dimension];
  Append[observable, "MinimumEpsilonShift" -> shift]];

nativePairStreamObservableMetadata[observable_Association] :=
  Join[<|"Identity" -> observable["Identity"],
    "CheckpointIdentity" -> observable["CheckpointIdentity"],
    "Epsilon" -> observable["Epsilon"],
    "TailPolicy" -> "stored"|>,
    If[KeyExistsQ[observable, "DivergentCancellation"],
      <|"DivergentCancellation" ->
        observable["DivergentCancellation"]|>, <||>]];

nativeContractStoredPairObservableStreamed[lowerState_Association,
    upperState_Association, observable_Association,
    lowerRecipes_List, upperRecipes_List, var_Symbol, domain_String,
    checkpointRoot_String] := Module[
  {stream = None, response, row = None, completed = False, output,
   tag = Unique["nativePairTileStream"], cleanup},
  cleanup[] := Module[{},
    Clear[row];
    nativeDropRationalMultiplierPreparationCache[];
    If[AssociationQ[stream] && !completed,
      Quiet[DiffExp2`CppBackend`AbortPersistentTransportPairObservableStream[
        stream]]];
    Null];
  output = Catch[
    CheckAbort[Catch[
      stream =
        DiffExp2`CppBackend`BeginPersistentTransportPairObservableStream[
          lowerState, upperState,
          nativePairStreamObservableMetadata[observable],
          checkpointRoot];
      If[FailureQ[stream] || !AssociationQ[stream] ||
          Lookup[stream, "status", "error"] =!= "ok",
        Throw[stream, tag]];
      Do[
        nativeStageTiming["paired-row-prepare-start side=lower tile=", tile];
        row = nativePrepareArmRecipeRow[lowerRecipes[[tile]],
          observable["CoefficientVector"], var, domain];
        nativeStageTiming["paired-row-prepare-done side=lower tile=", tile];
        response =
          DiffExp2`CppBackend`AddPersistentTransportPairObservableStreamTile[
            stream, "lower", tile - 1, row];
        nativeStageTiming["paired-tile-add-done side=lower tile=", tile];
        Clear[row];
        nativeDropRationalMultiplierPreparationCache[];
        If[FailureQ[response] || !AssociationQ[response] ||
            Lookup[response, "status", "error"] =!= "ok",
          Throw[response, tag]],
        {tile, Length[lowerRecipes]}];
      Do[
        nativeStageTiming["paired-row-prepare-start side=upper tile=", tile];
        row = nativePrepareArmRecipeRow[upperRecipes[[tile]],
          observable["CoefficientVector"], var, domain];
        nativeStageTiming["paired-row-prepare-done side=upper tile=", tile];
        response =
          DiffExp2`CppBackend`AddPersistentTransportPairObservableStreamTile[
            stream, "upper", tile - 1, row];
        nativeStageTiming["paired-tile-add-done side=upper tile=", tile];
        Clear[row];
        nativeDropRationalMultiplierPreparationCache[];
        If[FailureQ[response] || !AssociationQ[response] ||
            Lookup[response, "status", "error"] =!= "ok",
          Throw[response, tag]],
        {tile, Length[upperRecipes]}];
      nativeStageTiming["paired-stream-finish-start"];
      response =
        DiffExp2`CppBackend`FinishPersistentTransportPairObservableStream[
          stream];
      nativeStageTiming["paired-stream-finish-done"];
      If[AssociationQ[response] &&
          Lookup[response, "status", "error"] === "ok",
        completed = True];
      response,
      tag], cleanup[]; Abort[]],
    "DiffExp2Error", Function[{failure, errorTag},
      cleanup[]; Throw[failure, errorTag]]];
  cleanup[];
  output];

nativeBatchReleasePublished[states_, pairs_, lowerEndpointBatches_,
    upperEndpointBatches_] := Module[
  {pairList, lowerList, upperList, lines, endpoints},
  pairList = Select[If[ListQ[pairs], pairs, {pairs}], AssociationQ];
  lowerList = Select[If[ListQ[lowerEndpointBatches],
    lowerEndpointBatches, {lowerEndpointBatches}], AssociationQ];
  upperList = Select[If[ListQ[upperEndpointBatches],
    upperEndpointBatches, {upperEndpointBatches}], AssociationQ];
  lines = Flatten[Lookup[pairList, "lines", {}], 1];
  endpoints = Flatten[Join[Lookup[lowerList, "endpoints", {}],
    Lookup[upperList, "endpoints", {}]], 1];
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLineIntegral[#]] &,
    Select[lines, AssociationQ]];
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentEndpoint[#]] &,
    Select[endpoints, AssociationQ]];
  If[AssociationQ[states],
    Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentTransportArm[#]] &,
      Select[Values[states], AssociationQ]]];
  Null];

nativeChunkedObservableSummary[responses_List, chunkSize_Integer,
    resultKey_String] := If[responses === {}, None, <|
  "ChunkSize" -> chunkSize,
  "ChunkCount" -> Length[responses],
  "Chunks" -> (KeyDrop[#, {"status", resultKey}] & /@ responses)|>];

nativeStreamTransportArm[atlas_Association, data_Association,
    arm_String, epsilon_Association, checkpointRoot_String,
    refinement_Association] := Module[
  {systems = Rest[data["ChartSystems"]], current = atlas["Anchor"],
   tiles = {atlas["Anchor"]}, hopEpsilon, basis, response, next, output},
  hopEpsilon = <|"min" -> epsilon["min"], "max" -> epsilon["max"],
    "required_complete_max" -> epsilon["match_required_complete_max"]|>;
  output = Catch[
    Do[
      nativeStageTiming["stream-basis-start arm=", arm,
        " index=", index];
      basis = nativeReceivingBasis[systems[[index]], atlas["Request"],
        Lookup[atlas, "Threads", Automatic]];
      nativeStageTiming["stream-basis-done arm=", arm,
        " index=", index,
        " center=", InputForm[Lookup[systems[[index]], "Center", None]],
        " regular=", TrueQ[Lookup[
          Lookup[systems[[index]], "IndicialData", <||>],
          "Regular", False]],
        " elapsedMs=", Lookup[
          Lookup[basis, "NativeSummary", <||>], "elapsed_ms", None],
        " workers=", Lookup[
          Lookup[basis, "NativeSummary", <||>], "worker_threads", None],
        " capability=", Lookup[
          Lookup[basis, "NativeSummary", <||>],
          "execution_capability", Lookup[
            Lookup[basis, "NativeSummary", <||>],
            "selection_capability", None]]];
      response = DiffExp2`CppBackend`ConsumePersistentTransportHop[
        atlas["Plan"], arm, index, basis["Columns"], current,
        hopEpsilon, checkpointRoot, refinement];
      If[FailureQ[response] || !AssociationQ[response] ||
          Lookup[response, "status", "error"] =!= "ok",
        If[AssociationQ[basis],
          Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
            Lookup[basis, "Columns", {}]]];
        err["E5", <|"Arm" -> arm, "Match" -> index,
          "BackendFailure" -> response,
          "Detail" -> "streamed native transport hop failed"|>]];
      next = Lookup[response, "next_local", None];
      If[AssociationQ[next] && !KeyExistsQ[next, "session"],
        next = Append[next, "session" -> atlas["Session"]]];
      If[!nativeOpaqueLocalHandleQ[next, atlas["Session"]] ||
          !ListQ[Lookup[response, "consumed_basis_handles", None]],
        err["E5", <|"Arm" -> arm, "Match" -> index,
          "BackendResponse" -> response,
          "Detail" -> "streamed native transport hop returned a malformed consumed-local result"|>]];
      AppendTo[tiles, next];
      current = next;
      basis = None;
      ClearSystemCache[];
      nativeStageTiming["stream-hop-done arm=", arm,
        " index=", index],
      {index, Length[systems]}];
    <|"tile_sources" -> tiles|>,
    "DiffExp2Error", Function[{failure, tag},
      Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
        Rest[tiles]];
      Throw[failure, tag]]];
  output];

Options[RunNativeTransportObservableBatch] = {
  "MaxRefinementSteps" -> 2,
  "ObservableContractionChunkSize" -> 1,
  "ConsumeReceivingBases" -> False,
  "AtlasCompactor" -> Identity};

RunNativeTransportObservableBatch[atlasIn_Association, observables_List,
    var_Symbol, OptionsPattern[]] := Module[
  {atlas = atlasIn,
   maxSteps = OptionValue["MaxRefinementSteps"], normalized, identities,
   contractionChunkSize =
     OptionValue["ObservableContractionChunkSize"],
   checkpoints, prepared, sourceMin, availableMax, projectedRequired,
   contractionRequired, statePublicRequired,
   matchRequired, preparedShift, declaredShift, epsilon, refinement,
   lowerBasis, upperBasis, checkpointRoot, march, states = <||>,
   lowerState, upperState, integrates, lowerLimits, upperLimits,
   pairResponses = {}, lowerEndpointResponses = {},
   upperEndpointResponses = {}, chunkResponse, chunkRecords,
   observableChunks, chunk, chunkIndices, integrateIndices,
   lowerLimitIndices, upperLimitIndices, legacyObservables,
   endpointObservables, lowerRowRecipes, upperRowRecipes, row,
   pairRecords = {}, lowerRecords = {}, upperRecords = {},
   pairByIdentity = <||>, lowerByIdentity = <||>, upperByIdentity = <||>,
   resultRecords, output, cacheMemoryBefore,
   cacheMemoryAfter, workingAtlas = atlas,
   atlasCompactor = OptionValue["AtlasCompactor"],
   consumeReceivingBases = TrueQ[OptionValue["ConsumeReceivingBases"]],
   deferredReceivingBases, streamedArms = <||>},
  If[Lookup[atlas, "Type", None] =!=
        "DiffExp2NativeRegularIndependentArmAtlas" ||
      !IntegerQ[maxSteps] || !TrueQ[0 <= maxSteps <= 32],
    err["E8", <|"AtlasType" -> Lookup[atlas, "Type", None],
      "MaxRefinementSteps" -> maxSteps,
      "Detail" -> "native transport observable batch received a malformed atlas or refinement bound"|>]];
  If[!IntegerQ[contractionChunkSize] ||
      !TrueQ[1 <= contractionChunkSize <= 32],
    err["E8", <|"ObservableContractionChunkSize" ->
      contractionChunkSize,
      "Detail" -> "native observable contraction chunk size must be an integer from 1 through 32"|>]];
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
  (* Capture only the chart geometry and honest local rectangle needed to
     prepare one projected row.  The complete chart systems/bases may then
     be compacted before the march, while no prepared multiplier or row is
     retained in Wolfram. *)
  lowerRowRecipes = nativeArmRowRecipes[atlas, atlas["Lower"]];
  upperRowRecipes = nativeArmRowRecipes[atlas, atlas["Upper"]];
  prepared = nativePrepareBatchObservable[#, var, atlas["Dimension"]] & /@
    normalized;
  deferredReceivingBases = TrueQ[
    Lookup[atlas, "DeferredReceivingBases", False]];
  workingAtlas = If[deferredReceivingBases, atlas, Join[atlas, <|
    "Lower" -> KeyDrop[atlas["Lower"], "ChartSystems"],
    "Upper" -> KeyDrop[atlas["Upper"], "ChartSystems"]|>]];
  If[Head[atlasCompactor] =!= Function && atlasCompactor =!= Identity,
    err["E8", <|"Detail" ->
      "AtlasCompactor must be Identity or a one-argument Function"|>]];
  atlasCompactor[workingAtlas];
  (* Keep the large input association in exactly one short-lived local.
     Using the pattern argument throughout this Module substitutes the full
     atlas into the evaluated body, so neither the caller compactor nor cache
     clearing can reclaim ChartSystems before the native march allocation. *)
  atlas = workingAtlas;
  If[Head[atlasCompactor] === Function,
    nativeStageTiming["atlas-owner-compacted"]];
  (* The compact row recipes are self-contained and every chart/SCC/local
     needed by the march is strongly retained in C++.  Drop duplicate
     Wolfram preparation memo state before the two-arm allocation peak;
     rows themselves will be prepared only after the march, one tile at a
     time.  Never close a session or release an opaque native owner here. *)
  cacheMemoryBefore = MemoryInUse[];
  If[!deferredReceivingBases,
    DiffExp2`Solve`DropWolframPreparationCaches[]];
  ClearSystemCache[];
  Share[];
  cacheMemoryAfter = MemoryInUse[];
  nativeStageTiming["post-atlas-cache-drop before=", cacheMemoryBefore,
    " after=", cacheMemoryAfter,
    " released=", cacheMemoryBefore - cacheMemoryAfter];
  sourceMin = workingAtlas["Request", "EpsWindow", "Min"];
  availableMax = workingAtlas["Request", "EpsWindow", "CompleteMax"];
  projectedRequired = Max[Lookup[Lookup[prepared, "Epsilon"],
    "RequiredCompleteMax"]];
  (* Retain the unprojected source through every downstream contraction.
     Multiplication by a row beginning at eps^s needs source order D-s to
     deliver D.  An integrated monomial may additionally begin at eps^-1,
     so reserve one more source order for integrate observables.  The exact
     ladder preplanner already charges this same 1-s loss; keeping it here
     prevents the arm state from discarding that reserved coefficient before
     paired contraction. *)
  contractionRequired = Max[Map[
    #1["Epsilon", "RequiredCompleteMax"] -
      #1["MinimumEpsilonShift"] +
      If[#1["Operation"] === "integrate", 1, 0] &,
    prepared]];
  statePublicRequired = Max[sourceMin, contractionRequired];
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
      "ContractionRequiredCompleteMax" -> contractionRequired,
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
  lowerBasis = If[deferredReceivingBases, {},
    Lookup[Rest[workingAtlas["Lower", "Bases"]], "Columns", {}]];
  upperBasis = If[deferredReceivingBases, {},
    Lookup[Rest[workingAtlas["Upper", "Bases"]], "Columns", {}]];
  (* Bind the exact physical observable before any row exists.  Each native
     add_tile additionally binds the actual prepared-row identity, so the
     final line provenance proves both this physical request and every lazy
     chart specialization without requiring an eager all-row identity list. *)
  checkpointRoot = nativeCheckpointIdentity["de2-native-observable-batch-", <|
    "schema" ->
      "diffexp2-native-observable-lazy-row-checkpoint-identity-v1",
    "plan_checkpoint_identity" ->
      workingAtlas["PlanCheckpointIdentity"],
    "domain" -> workingAtlas["Domain"],
    "anchor_checkpoint_identity" ->
      nativeOpaqueCheckpoint[workingAtlas["Anchor"], "anchor"],
    "epsilon" -> epsilon, "refinement" -> refinement,
    "observables" -> Map[KeyTake[#, {"Operation", "Identity",
      "CheckpointIdentity", "CoefficientVector", "Epsilon",
      "TailPolicy", "MinimumEpsilonShift",
      "DivergentCancellation"}] &, prepared]|>];
  output = Catch[
    march = If[deferredReceivingBases,
      streamedArms = <|
        "lower" -> nativeStreamTransportArm[workingAtlas,
          workingAtlas["Lower"], "lower", epsilon,
          checkpointRoot <> ":march", refinement],
        "upper" -> nativeStreamTransportArm[workingAtlas,
          workingAtlas["Upper"], "upper", epsilon,
          checkpointRoot <> ":march", refinement]|>;
      nativeStageTiming["transport-state-publish-start"];
      DiffExp2`CppBackend`PublishPersistentConsumedTransportStates[
        workingAtlas["Plan"], workingAtlas["Anchor"], streamedArms,
        epsilon, checkpointRoot <> ":march", refinement],
      If[consumeReceivingBases,
        DiffExp2`CppBackend`RunPersistentConsumingTransportArms,
        DiffExp2`CppBackend`RunPersistentTransportArms][
        workingAtlas["Plan"], workingAtlas["Anchor"], <|
          "lower" -> <|"receiving_basis" -> lowerBasis|>,
          "upper" -> <|"receiving_basis" -> upperBasis|>|>,
        epsilon, checkpointRoot <> ":march", refinement]];
    nativeStageTiming["transport-state-publish-done"];
    states = If[AssociationQ[march] &&
        AssociationQ[Lookup[march, "states", None]],
      march["states"], <||>];
    If[FailureQ[march] || !AssociationQ[march] ||
        Lookup[march, "status", "error"] =!= "ok" ||
        Lookup[march, "session", None] =!= workingAtlas["Session"] ||
        !TrueQ[Lookup[march, "native_retained", False]] ||
        Lookup[march, "json_coefficients", None] =!= 0,
      err["E5", <|"BackendFailure" -> march,
        "Detail" -> "persistent two-arm observable-batch march failed"|>]];
    lowerState = Lookup[states, "lower", <||>];
    upperState = Lookup[states, "upper", <||>];
    If[Sort[Keys[states]] =!= {"lower", "upper"} ||
        !nativeTransportStateHandleQ[lowerState, workingAtlas["Session"],
          "lower"] ||
        !nativeTransportStateHandleQ[upperState, workingAtlas["Session"],
          "upper"],
      err["E5", <|"BackendResponse" -> march,
        "Detail" -> "two-arm observable-batch march did not return exact lower/upper retained states"|>]];
    If[consumeReceivingBases || deferredReceivingBases,
      workingAtlas = Join[workingAtlas, <|
        "Lower" -> Join[KeyDrop[workingAtlas["Lower"], "ChartSystems"],
          <|"Bases" -> {None}|>],
        "Upper" -> Join[KeyDrop[workingAtlas["Upper"], "ChartSystems"],
          <|"Bases" -> {None}|>]|>];
      atlas = workingAtlas;
      atlasCompactor[workingAtlas];
      If[deferredReceivingBases,
        DiffExp2`Solve`DropWolframPreparationCaches[]]];
    integrateIndices = Select[Range[Length[prepared]],
      prepared[[#, "Operation"]] === "integrate" &];
    lowerLimitIndices = Select[Range[Length[prepared]],
      prepared[[#, "Operation"]] === "limitLower" &];
    upperLimitIndices = Select[Range[Length[prepared]],
      prepared[[#, "Operation"]] === "limitUpper" &];
    integrates = prepared[[integrateIndices]];
    lowerLimits = prepared[[lowerLimitIndices]];
    upperLimits = prepared[[upperLimitIndices]];
    If[integrates =!= {},
      nativeStageTiming["paired-contraction-start observables=",
        Length[integrates]];
      observableChunks = Partition[integrateIndices,
        UpTo[contractionChunkSize]];
      Do[
        chunkIndices = observableChunks[[chunkIndex]];
        chunk = prepared[[chunkIndices]];
        If[Length[chunk] === 1 &&
            First[chunk]["TailPolicy"] === "stored",
          (* The default singleton path never owns more than one prepared
             row: prepare, synchronously add, clear the row and multiplier
             memo, then advance to the next tile. *)
          chunkResponse = nativeContractStoredPairObservableStreamed[
            lowerState, upperState, First[chunk], lowerRowRecipes,
            upperRowRecipes, var, workingAtlas["Domain"],
            checkpointRoot <> ":integrals:chunk:" <>
              ToString[chunkIndex]],
          (* Tail certification still uses the legacy atomic batch API.
             Likewise, an explicit chunk size above one preserves its
             requested batched behavior.  Only this bounded chunk is eager. *)
          legacyObservables = Map[Function[observable,
            Module[{lowerRows, upperRows},
              lowerRows = Map[Function[recipe,
                row = nativePrepareArmRecipeRow[recipe,
                  observable["CoefficientVector"], var,
                  workingAtlas["Domain"]];
                nativeDropRationalMultiplierPreparationCache[];
                row], lowerRowRecipes];
              upperRows = Map[Function[recipe,
                row = nativePrepareArmRecipeRow[recipe,
                  observable["CoefficientVector"], var,
                  workingAtlas["Domain"]];
                nativeDropRationalMultiplierPreparationCache[];
                row], upperRowRecipes];
              Join[<|"Identity" -> observable["Identity"],
                "CheckpointIdentity" ->
                  observable["CheckpointIdentity"],
                "LowerIntegrandRows" -> lowerRows,
                "UpperIntegrandRows" -> upperRows,
                "Epsilon" -> observable["Epsilon"],
                "TailPolicy" -> observable["TailPolicy"]|>,
                If[KeyExistsQ[observable, "DivergentCancellation"],
                  <|"DivergentCancellation" ->
                    observable["DivergentCancellation"]|>, <||>]]]],
            chunk];
          chunkResponse =
            DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
              lowerState, upperState, legacyObservables,
              checkpointRoot <> ":integrals:chunk:" <>
                ToString[chunkIndex]];
          Clear[legacyObservables, row];
          nativeDropRationalMultiplierPreparationCache[]];
        (* Record even a malformed current response before validation: a
           backend failure may still have published handles that the batch
           failure path must release. *)
        If[AssociationQ[chunkResponse],
          AppendTo[pairResponses, chunkResponse]];
        If[FailureQ[chunkResponse] || !AssociationQ[chunkResponse] ||
            Lookup[chunkResponse, "status", "error"] =!= "ok",
          err["E5", <|"BackendFailure" -> chunkResponse,
            "ChunkIndex" -> chunkIndex,
            "Detail" -> "native paired observable contraction chunk failed"|>]];
        chunkRecords = nativeOrderedObservableOutputs[
          chunkResponse, chunk, workingAtlas["Session"], "line"];
        pairRecords = Join[pairRecords, chunkRecords],
        {chunkIndex, Length[observableChunks]}];
      nativeStageTiming["paired-contraction-done"];
      pairByIdentity = AssociationThread[Lookup[integrates, "Identity"],
        pairRecords]];
    If[lowerLimits =!= {},
      observableChunks = Partition[lowerLimitIndices,
        UpTo[contractionChunkSize]];
      Do[
        chunkIndices = observableChunks[[chunkIndex]];
        chunk = prepared[[chunkIndices]];
        endpointObservables = Map[Function[observable,
          row = nativePrepareArmRecipeRow[Last[lowerRowRecipes],
            observable["CoefficientVector"], var,
            workingAtlas["Domain"]];
          nativeDropRationalMultiplierPreparationCache[];
          <|"Identity" -> observable["Identity"],
            "CheckpointIdentity" -> observable["CheckpointIdentity"],
            "IntegrandRow" -> row,
            "Epsilon" -> observable["Epsilon"]|>], chunk];
        chunkResponse =
          DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
            lowerState, endpointObservables,
            checkpointRoot <> ":lower-endpoints:chunk:" <>
              ToString[chunkIndex]];
        Clear[endpointObservables, row];
        nativeDropRationalMultiplierPreparationCache[];
        If[AssociationQ[chunkResponse],
          AppendTo[lowerEndpointResponses, chunkResponse]];
        If[FailureQ[chunkResponse] || !AssociationQ[chunkResponse] ||
            Lookup[chunkResponse, "status", "error"] =!= "ok",
          err["E5", <|"BackendFailure" -> chunkResponse,
            "ChunkIndex" -> chunkIndex,
            "Detail" -> "native lower endpoint observable chunk failed"|>]];
        chunkRecords = nativeOrderedObservableOutputs[
          chunkResponse, chunk, workingAtlas["Session"], "endpoint"];
        lowerRecords = Join[lowerRecords, chunkRecords],
        {chunkIndex, Length[observableChunks]}];
      lowerByIdentity = AssociationThread[Lookup[lowerLimits, "Identity"],
        lowerRecords]];
    If[upperLimits =!= {},
      observableChunks = Partition[upperLimitIndices,
        UpTo[contractionChunkSize]];
      Do[
        chunkIndices = observableChunks[[chunkIndex]];
        chunk = prepared[[chunkIndices]];
        endpointObservables = Map[Function[observable,
          row = nativePrepareArmRecipeRow[Last[upperRowRecipes],
            observable["CoefficientVector"], var,
            workingAtlas["Domain"]];
          nativeDropRationalMultiplierPreparationCache[];
          <|"Identity" -> observable["Identity"],
            "CheckpointIdentity" -> observable["CheckpointIdentity"],
            "IntegrandRow" -> row,
            "Epsilon" -> observable["Epsilon"]|>], chunk];
        chunkResponse =
          DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
            upperState, endpointObservables,
            checkpointRoot <> ":upper-endpoints:chunk:" <>
              ToString[chunkIndex]];
        Clear[endpointObservables, row];
        nativeDropRationalMultiplierPreparationCache[];
        If[AssociationQ[chunkResponse],
          AppendTo[upperEndpointResponses, chunkResponse]];
        If[FailureQ[chunkResponse] || !AssociationQ[chunkResponse] ||
            Lookup[chunkResponse, "status", "error"] =!= "ok",
          err["E5", <|"BackendFailure" -> chunkResponse,
            "ChunkIndex" -> chunkIndex,
            "Detail" -> "native upper endpoint observable chunk failed"|>]];
        chunkRecords = nativeOrderedObservableOutputs[
          chunkResponse, chunk, workingAtlas["Session"], "endpoint"];
        upperRecords = Join[upperRecords, chunkRecords],
        {chunkIndex, Length[observableChunks]}];
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
      "Atlas" -> workingAtlas, "States" -> states,
      "Results" -> resultRecords, "NativeMarches" -> 2,
      "CompatibilityExports" -> 0,
      "NativeSummary" -> <|
        "March" -> KeyDrop[march, {"status", "states"}],
        "PairCalls" -> Length[pairResponses],
        "Pair" -> nativeChunkedObservableSummary[
          pairResponses, contractionChunkSize, "lines"],
        "LowerEndpointCalls" -> Length[lowerEndpointResponses],
        "LowerEndpoints" -> nativeChunkedObservableSummary[
          lowerEndpointResponses, contractionChunkSize, "endpoints"],
        "UpperEndpointCalls" -> Length[upperEndpointResponses],
        "UpperEndpoints" -> nativeChunkedObservableSummary[
          upperEndpointResponses, contractionChunkSize, "endpoints"]|>|>,
    "DiffExp2Error", Function[{failure, tag},
      Clear[row, legacyObservables, endpointObservables];
      nativeDropRationalMultiplierPreparationCache[];
      nativeBatchReleasePublished[states, pairResponses,
        lowerEndpointResponses, upperEndpointResponses];
      If[AssociationQ[streamedArms],
        Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
          Flatten[Rest[Lookup[#, "tile_sources", {}]] & /@
            Select[Values[streamedArms], AssociationQ], 1]]];
      Throw[failure, tag]]];
  output];

SetAttributes[RunNativeTransportObservableBatchOwned, HoldFirst];
Options[RunNativeTransportObservableBatchOwned] = {
  "MaxRefinementSteps" -> 2,
  "ObservableContractionChunkSize" -> 1};
RunNativeTransportObservableBatchOwned[owner_Symbol, observables_List,
    var_Symbol, OptionsPattern[]] :=
  RunNativeTransportObservableBatch[owner, observables, var,
    "MaxRefinementSteps" -> OptionValue["MaxRefinementSteps"],
    "ObservableContractionChunkSize" ->
      OptionValue["ObservableContractionChunkSize"],
    "ConsumeReceivingBases" -> True,
    "AtlasCompactor" -> Function[compactAtlas, owner = compactAtlas]];

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

nativeObservableCheckpointResult[result_Association, session_String,
    schema_:$nativeObservableCheckpointSchemaV2] := Module[
  {operation, kind, handle, token, checkpoint, provenance, provenanceKey},
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
  provenanceKey = nativeCheckpointProvenanceKey[schema];
  If[provenanceKey === None,
    err["E8", <|"Schema" -> schema,
      "Detail" -> "native checkpoint result requested an unsupported manifest schema"|>]];
  Join[<|"RequestIndex" -> result["RequestIndex"],
    "Operation" -> operation, "Identity" -> result["Identity"],
    "CheckpointIdentity" -> checkpoint, "Epsilon" -> result["Epsilon"],
    "Kind" -> kind, "Handle" -> token|>,
    <|provenanceKey ->
      nativeCheckpointProvenanceReference[provenance, schema]|>]];

nativeObservableCheckpointManifestQ[manifest_] := Module[
  {keys = {"Schema", "Path", "CheckpointIdentity", "ManifestIdentity",
      "TransportArmMarches", "StateHandles", "Results"}, core, results,
   schema, provenanceKey,
   resultKeys = {"RequestIndex", "Operation", "Identity",
      "CheckpointIdentity", "Epsilon", "Kind", "Handle"}},
  If[!AssociationQ[manifest] || Sort[Keys[manifest]] =!= Sort[keys],
    Return[False, Module]];
  schema = manifest["Schema"];
  provenanceKey = nativeCheckpointProvenanceKey[schema];
  If[provenanceKey === None, Return[False, Module]];
  AppendTo[resultKeys, provenanceKey];
  core = KeyDrop[manifest, "ManifestIdentity"];
  results = manifest["Results"];
  TrueQ[MemberQ[{$nativeObservableCheckpointSchemaV1,
        $nativeObservableCheckpointSchemaV2}, schema]] &&
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
        provenanceKey}] &&
      nativeNonemptyStringQ[#["Handle"]] &&
      StringStartsQ[#["Handle"], "transport:"] &&
      nativeNonemptyStringQ[#["CheckpointIdentity"]] &&
      nativeCheckpointProvenanceReferenceQ[#[provenanceKey], schema] &] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]], "Handle"]] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]],
      "CheckpointIdentity"]] &&
    DuplicateFreeQ[Lookup[Values[manifest["StateHandles"]],
      provenanceKey]] &&
    ListQ[results] && results =!= {} &&
    AllTrue[results, AssociationQ[#] &&
      Sort[Keys[#]] === Sort[resultKeys] &] &&
    Lookup[results, "RequestIndex"] === Range[0, Length[results] - 1] &&
    DuplicateFreeQ[Lookup[results, "Identity"]] &&
    DuplicateFreeQ[Lookup[results, "CheckpointIdentity"]] &&
    DuplicateFreeQ[Lookup[results, "Handle"]] &&
    DuplicateFreeQ[Lookup[results, provenanceKey]] &&
    AllTrue[results, Function[result,
      nativeNonemptyStringQ[result["Identity"]] &&
      nativeNonemptyStringQ[result["CheckpointIdentity"]] &&
      nativeNonemptyStringQ[result["Handle"]] &&
      nativeCheckpointProvenanceReferenceQ[
        result[provenanceKey], schema] &&
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
   core, schema = $nativeObservableCheckpointSchemaV2,
   provenanceKey},
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
  provenanceKey = nativeCheckpointProvenanceKey[schema];
  stateHandles = AssociationMap[Function[side, Module[{provenance},
      provenance = Lookup[states[side], "provenance_identity",
        Lookup[states[side], "ProvenanceIdentity", None]];
      Join[<|
        "Handle" -> Lookup[states[side], "transport_state",
          Lookup[states[side], "TransportState", None]],
        "CheckpointIdentity" -> Lookup[states[side], "checkpoint_identity",
          Lookup[states[side], "CheckpointIdentity", None]]|>,
        <|provenanceKey -> nativeCheckpointProvenanceReference[
          provenance, schema]|>]]], {"lower", "upper"}];
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
  core = <|"Schema" -> schema,
    "Path" -> expanded, "CheckpointIdentity" -> identity,
    "TransportArmMarches" -> stats["transport_arm_marches"],
    "StateHandles" -> stateHandles,
    "Results" -> (nativeObservableCheckpointResult[#, session, schema] & /@
      results)|>;
  Append[core, "ManifestIdentity" -> nativeCheckpointIdentity[
    "de2-native-observable-checkpoint-manifest-", core]]];

RestoreNativeTransportObservableBatchCheckpoint[manifest_Association] :=
 Module[{restored, session = None, close, expectedLines, expectedEndpoints,
   expectedStates, stats, results, restoredHandles, restoredRecordMap,
   lineMap, endpointMap, stateMap, restoredResultIdentitiesQ,
   restoredStateIdentitiesQ, schema, provenanceKey},
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
  schema = manifest["Schema"];
  provenanceKey = nativeCheckpointProvenanceKey[schema];
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
        nativeCheckpointProvenanceMatchesQ[
          Lookup[record, "provenance_identity", None],
          result[provenanceKey], schema]]]];
  restoredStateIdentitiesQ = AllTrue[Values[manifest["StateHandles"]],
    Function[state, Module[{record = Lookup[stateMap,
        state["Handle"], None]},
      AssociationQ[record] &&
        Lookup[record, "checkpoint_identity", None] ===
          state["CheckpointIdentity"] &&
        nativeCheckpointProvenanceMatchesQ[
          Lookup[record, "provenance_identity", None],
          state[provenanceKey], schema]]]];
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
  results = Map[Function[result, Module[{record},
    record = Lookup[If[result["Kind"] === "line", lineMap, endpointMap],
      result["Handle"]];
    Join[KeyTake[result, {"RequestIndex", "Operation", "Identity",
        "CheckpointIdentity", "Epsilon"}],
      <|If[result["Kind"] === "line", "Line", "Endpoint"] ->
        <|"session" -> session, result["Kind"] -> result["Handle"],
          "request_index" -> result["RequestIndex"],
          "observable_identity" -> result["Identity"],
          "checkpoint_identity" -> result["CheckpointIdentity"],
          "provenance_identity" -> record["provenance_identity"]|>|>]]],
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
    Select[Join[Rest[atlas["Lower", "Bases"]],
      Rest[atlas["Upper", "Bases"]]], AssociationQ]], 1];
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
