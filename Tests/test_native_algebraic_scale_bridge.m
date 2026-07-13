(* Exact no-FIRE regression for the regular algebraic-scale native bridge. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 80,
  "ChopPrecision" -> 40,
  "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0,
  "DivisionOrder" -> 3,
  "RadiusOfConvergence" -> 1,
  "Variables" -> {},
  "Verbosity" -> 0}];

x = Global`x;
goldenScale = RootReduce[(Sqrt[5] - 1)/2];
nativeScale = Floor[2^64*goldenScale]/2^64;
exactFloorQ = IntegerQ[Numerator[nativeScale]] &&
  IntegerQ[Denominator[nativeScale]] &&
  TrueQ[FullSimplify[RootReduce[0 < nativeScale < goldenScale]]];

system = DiffExp2`LoadSystem[<|
  "Matrix" -> {{0}}, "Variable" -> x|>];
singularities = <|"All" -> {}, "Real" -> {}, "Projected" -> {},
  "ProjectionWaypoints" -> {}|>;
minusDeltaPrescription = <|"Factor" -> 1 + x,
  "ExactFactor" -> 1 + x, "Sign" -> -1, "Multiplicity" -> 1,
  "LeadingCoeffSign" -> 1|>;

anchor = <|"Center" -> 0, "Singular" -> False,
  "Radius" -> 24/23, "MatchRadius" -> 24/23,
  "Scale" -> 24/23, "LocalRadius" -> 1,
  "ChartVar" -> Global`t, "UseSCCSkeleton" -> True,
  "Name" -> "bubble-anchor", "Prescriptions" -> {}|>;
endpoint = <|"Center" -> 1/2, "Singular" -> False,
  "Radius" -> goldenScale, "MatchRadius" -> goldenScale,
  "Scale" -> goldenScale, "LocalRadius" -> 1,
  "IncomingMatchPoint" -> 1/3, "SymmetricMatch" -> False,
  "ChartVar" -> Global`t, "UseSCCSkeleton" -> True,
  "Name" -> "bubble-algebraic-endpoint",
  "Prescriptions" -> {minusDeltaPrescription}|>;

lowerPlan = <|"From" -> 0, "To" -> -1/2, "Direction" -> -1,
  "Charts" -> {anchor}, "SegmentCount" -> 1,
  "EndpointIsSingular" -> False, "DigitsNeeded" -> 20,
  "Singularities" -> singularities|>;
upperPlan = <|"From" -> 0, "To" -> 1/2, "Direction" -> 1,
  "Charts" -> {anchor, endpoint}, "SegmentCount" -> 2,
  "EndpointIsSingular" -> False, "DigitsNeeded" -> 20,
  "Singularities" -> singularities|>;

eligibleQ = TrueQ[
  DiffExp2`NativeTransport`NativeRegularIndependentArmPlansSupportedQ[
    lowerPlan, upperPlan]];
singularBridge = catchDE2[
  DiffExp2`NativeTransport`Private`bridgeNativeRegularChartScale[
    Join[endpoint, <|"Singular" -> True|>], 2]];
singularScopeQ = AssociationQ[singularBridge] &&
  TrueQ[singularBridge["Singular"]] &&
  singularBridge["Center"] === 1/2 &&
  singularBridge["Scale"] === nativeScale &&
  singularBridge["Prescriptions"] === endpoint["Prescriptions"] &&
  TrueQ[FullSimplify[RootReduce[
    0 < singularBridge["Radius"] < endpoint["Radius"]]]];

inexactRadiusSource = Sqrt[41/4];
inexactRadiusChart = <|"Center" -> 0, "Singular" -> False,
  "Radius" -> N[inexactRadiusSource, 80], "MatchRadius" -> 1/11,
  "Scale" -> 1/11, "LocalRadius" -> N[11 inexactRadiusSource, 80],
  "UseSCCSkeleton" -> True, "Name" -> "certified-inexact-radius",
  "Prescriptions" -> {minusDeltaPrescription}|>;
inexactRadiusBridge = catchDE2[
  DiffExp2`NativeTransport`Private`bridgeNativeRegularChartScale[
    inexactRadiusChart, 3]];
inexactRadiusCertificate = If[AssociationQ[inexactRadiusBridge],
  Lookup[inexactRadiusBridge, "NativeRationalScaleBridge", <||>], <||>];
inexactRadiusBridgeQ = AssociationQ[inexactRadiusBridge] &&
  Lookup[inexactRadiusCertificate, "Schema", None] ===
    "diffexp2-native-inward-rational-radius-v1" &&
  inexactRadiusBridge["Center"] === inexactRadiusChart["Center"] &&
  inexactRadiusBridge["Scale"] === inexactRadiusChart["Scale"] &&
  inexactRadiusBridge["MatchRadius"] ===
    inexactRadiusChart["MatchRadius"] &&
  inexactRadiusBridge["Prescriptions"] ===
    inexactRadiusChart["Prescriptions"] &&
  IntegerQ[Numerator[inexactRadiusBridge["Radius"]]] &&
  IntegerQ[Denominator[inexactRadiusBridge["Radius"]]] &&
  IntegerQ[Numerator[inexactRadiusBridge["LocalRadius"]]] &&
  IntegerQ[Denominator[inexactRadiusBridge["LocalRadius"]]] &&
  TrueQ[0 < inexactRadiusBridge["Radius"] < inexactRadiusSource] &&
  TrueQ[inexactRadiusBridge["Radius"] +
      inexactRadiusCertificate["OriginalRadiusUncertainty"] <=
    inexactRadiusCertificate["OriginalRadiusMidpoint"]] &&
  inexactRadiusBridge["LocalRadius"] ===
    Together[inexactRadiusBridge["Radius"]/
      Abs[inexactRadiusChart["Scale"]]];

(* An exact algebraic physical radius is a stronger input than the interval
   certificate above.  Keep the match disk deliberately close to the true
   radius: the bridge must choose a strict inward rational floor without
   shrinking past the exact rational handoff disk. *)
exactRadiusSource = RootReduce[Sqrt[41/4]];
exactRadiusMatch = 16/5;
exactRadiusChart = <|"Center" -> 0, "Singular" -> False,
  "Radius" -> exactRadiusSource, "MatchRadius" -> exactRadiusMatch,
  "Scale" -> exactRadiusMatch,
  "LocalRadius" -> Together[exactRadiusSource/exactRadiusMatch],
  "UseSCCSkeleton" -> True, "Name" -> "exact-algebraic-radius",
  "Prescriptions" -> {minusDeltaPrescription}|>;
exactRadiusBridge = catchDE2[
  DiffExp2`NativeTransport`Private`bridgeNativeRegularChartScale[
    exactRadiusChart, 4]];
exactRadiusCertificate = If[AssociationQ[exactRadiusBridge],
  Lookup[exactRadiusBridge, "NativeRationalScaleBridge", <||>], <||>];
exactRadiusBridgeQ = AssociationQ[exactRadiusBridge] &&
  Lookup[exactRadiusCertificate, "Schema", None] ===
    "diffexp2-native-inward-rational-algebraic-radius-v1" &&
  Lookup[exactRadiusCertificate, "BridgeKind", None] ===
    "ExactAlgebraicPhysicalRadius" &&
  Lookup[exactRadiusCertificate, "FloorBits", None] === 64 &&
  Lookup[exactRadiusCertificate, "OriginalPhysicalRadius", None] ===
    exactRadiusSource &&
  Lookup[exactRadiusCertificate, "NativePhysicalRadius", None] ===
    exactRadiusBridge["Radius"] &&
  Lookup[exactRadiusCertificate, "OriginalLocalRadius", None] ===
    exactRadiusChart["LocalRadius"] &&
  Lookup[exactRadiusCertificate, "NativeLocalRadius", None] ===
    exactRadiusBridge["LocalRadius"] &&
  TrueQ[Lookup[exactRadiusCertificate, "ExactContainmentProved", False]] &&
  exactRadiusBridge["Center"] === exactRadiusChart["Center"] &&
  exactRadiusBridge["Scale"] === exactRadiusChart["Scale"] &&
  exactRadiusBridge["MatchRadius"] === exactRadiusMatch &&
  exactRadiusBridge["Prescriptions"] ===
    exactRadiusChart["Prescriptions"] &&
  IntegerQ[Numerator[exactRadiusBridge["Radius"]]] &&
  IntegerQ[Denominator[exactRadiusBridge["Radius"]]] &&
  IntegerQ[Numerator[exactRadiusBridge["LocalRadius"]]] &&
  IntegerQ[Denominator[exactRadiusBridge["LocalRadius"]]] &&
  TrueQ[FullSimplify[RootReduce[
    0 < exactRadiusMatch < exactRadiusBridge["Radius"] <
      exactRadiusSource]]] &&
  exactRadiusBridge["LocalRadius"] ===
    Together[exactRadiusBridge["Radius"]/exactRadiusMatch];

(* SegmentLine may cap the same retained anchor differently on its two
   arms.  Normalization must select the larger exact disk, after which both
   copies receive one identical conservative rational bridge. *)
sharedLowerRadius = RootReduce[Sqrt[10]];
sharedUpperRadius = exactRadiusSource;
sharedLowerMatch = 31/10;
sharedUpperMatch = 16/5;
sharedSelectedMatch = Max[sharedLowerMatch, sharedUpperMatch];
sharedLowerAnchor = <|"Center" -> 0, "Singular" -> False,
  "Radius" -> sharedLowerRadius, "MatchRadius" -> sharedLowerMatch,
  "Scale" -> sharedLowerMatch,
  "LocalRadius" -> Together[sharedLowerRadius/sharedLowerMatch],
  "ChartVar" -> Global`t, "UseSCCSkeleton" -> True,
  "Name" -> "exact-shared-lower", "Prescriptions" -> {}|>;
sharedUpperAnchor = Join[sharedLowerAnchor, <|
  "Radius" -> sharedUpperRadius,
  "MatchRadius" -> sharedUpperMatch, "Scale" -> sharedUpperMatch,
  "LocalRadius" -> Together[sharedUpperRadius/sharedUpperMatch],
  "Name" -> "exact-shared-upper"|>];
sharedLowerPlan = <|"From" -> 0, "To" -> -1, "Direction" -> -1,
  "Charts" -> {sharedLowerAnchor}, "SegmentCount" -> 1,
  "EndpointIsSingular" -> False, "DigitsNeeded" -> 20,
  "Singularities" -> singularities|>;
sharedUpperPlan = <|"From" -> 0, "To" -> 1, "Direction" -> 1,
  "Charts" -> {sharedUpperAnchor}, "SegmentCount" -> 1,
  "EndpointIsSingular" -> False, "DigitsNeeded" -> 20,
  "Singularities" -> singularities|>;
sharedNormalized = catchDE2[
  DiffExp2`NativeTransport`Private`normalizeSharedAnchor[
    sharedLowerPlan, sharedUpperPlan]];
sharedNormalizationQ = ListQ[sharedNormalized] &&
  Length[sharedNormalized] === 2 &&
  Module[{lower = First[sharedNormalized]["Charts"][[1]],
      upper = Last[sharedNormalized]["Charts"][[1]], provenance},
    provenance = Lookup[lower, "SharedAnchorNormalization", <||>];
    lower === upper && lower["Center"] === 0 &&
      lower["Radius"] === sharedUpperRadius &&
      lower["MatchRadius"] === sharedSelectedMatch &&
      lower["Scale"] === sharedSelectedMatch &&
      lower["LocalRadius"] ===
        Together[sharedUpperRadius/sharedSelectedMatch] &&
      AssociationQ[provenance] &&
      Lookup[provenance, "Schema", None] ===
        "diffexp2-native-shared-anchor-normalization-v1" &&
      Lookup[provenance, "IncomingLowerPhysicalRadius", None] ===
        sharedLowerRadius &&
      Lookup[provenance, "IncomingUpperPhysicalRadius", None] ===
        sharedUpperRadius &&
      Lookup[provenance, "IncomingLowerMatchRadius", None] ===
        sharedLowerMatch &&
      Lookup[provenance, "IncomingUpperMatchRadius", None] ===
        sharedUpperMatch &&
      Lookup[provenance, "SelectedPhysicalRadius", None] ===
        sharedUpperRadius &&
      Lookup[provenance, "SelectedMatchRadius", None] ===
        sharedSelectedMatch];
sharedEligibleQ = TrueQ[
  DiffExp2`NativeTransport`NativeRegularIndependentArmPlansSupportedQ[
    sharedLowerPlan, sharedUpperPlan]];

(* Exact-equality overlap is not stable under two strict inward floors.  The
   public eligibility gate must reject this atlas before native preparation. *)
tightChart1 = Join[endpoint, <|"Center" -> 1/3,
  "Scale" -> Sqrt[2], "Radius" -> Sqrt[2],
  "MatchRadius" -> Sqrt[2], "Prescriptions" -> {}|>];
tightChart2 = Join[endpoint, <|"Center" -> 1,
  "Scale" -> 2 - Sqrt[2], "Radius" -> 2 - Sqrt[2],
  "MatchRadius" -> 2 - Sqrt[2], "Prescriptions" -> {}|>];
tightLowerPlan = Join[lowerPlan, <|"Charts" -> {anchor},
  "SegmentCount" -> 1|>];
tightUpperPlan = Join[upperPlan, <|"To" -> 1,
  "Charts" -> {anchor, tightChart1, tightChart2},
  "SegmentCount" -> 3|>];
tightOriginalEqualityQ = TrueQ[FullSimplify[RootReduce[
  1 - 1/3 == Sqrt[2]/3 + (2 - Sqrt[2])/3]]];
tightRejectedQ = tightOriginalEqualityQ &&
  !TrueQ[DiffExp2`NativeTransport`NativeRegularIndependentArmPlansSupportedQ[
      tightLowerPlan, tightUpperPlan]] &&
  !TrueQ[DiffExp2`NativeTransport`Private`nativeBridgedPlanPreflightQ[
    tightUpperPlan]];

projectionRoot = Sqrt[2] + I*Sqrt[3];
projectionPlan = <|"From" -> -1, "To" -> 4,
  "Singularities" -> <|
    "All" -> {projectionRoot, Conjugate[projectionRoot]},
    "Real" -> {},
    "Projected" -> {Sqrt[2] - Sqrt[3], Sqrt[2], Sqrt[2] + Sqrt[3]},
    "ProjectionWaypoints" -> {}|>|>;
projectionRecords = catchDE2[
  DiffExp2`NativeTransport`Private`nativeComplexProjections[
    projectionPlan]];
projectionBridgeQ = ListQ[projectionRecords] &&
  Length[projectionRecords] === 1 &&
  Module[{record = First[projectionRecords], re, h},
    re = ToExpression[record["real_part_exact"], InputForm];
    h = ToExpression[record["imaginary_magnitude_exact"], InputForm];
    (IntegerQ[re] || Head[re] === Rational) &&
      (IntegerQ[h] || Head[h] === Rational) && h > 0 &&
      TrueQ[record["retain_minus_imaginary"]] &&
      TrueQ[record["retain_real_part"]] &&
      TrueQ[record["retain_plus_imaginary"]] &&
      Abs[N[re, 70] - N[Sqrt[2], 70]] < 10^-60 &&
      Abs[N[h, 70] - N[Sqrt[3], 70]] < 10^-60];

exactNativeValue[string_String] := ToExpression[string, InputForm];
originalPhysicalRadius[chart_Association] := Module[
  {certificate = Lookup[chart, "NativeRationalScaleBridge", None]},
  If[AssociationQ[certificate], certificate["OriginalPhysicalRadius"],
    Abs[chart["Scale"]]*chart["LocalRadius"]]];
insideBothDisksQ[point_, chart_Association] :=
  TrueQ[FullSimplify[RootReduce[
    Abs[point - chart["Center"]] <
      Abs[chart["Scale"]]*chart["LocalRadius"]]]] &&
  TrueQ[FullSimplify[RootReduce[
    Abs[point - chart["Center"]] < originalPhysicalRadius[chart]]]];
exactArmGeometryQ[armStats_Association, plan_Association] := Module[
  {charts = plan["Charts"], nativeCharts = armStats["charts"],
   matches = armStats["matches"], tiles = armStats["tiles"],
   chartBindingsQ, matchesQ, tilesQ},
  If[Length[nativeCharts] =!= Length[charts] ||
      Length[matches] + 1 =!= Length[charts] ||
      Length[tiles] =!= Length[charts], Return[False, Module]];
  chartBindingsQ = And @@ Table[
    exactNativeValue[nativeCharts[[i, "center_exact"]]] ===
        charts[[i, "Center"]] &&
      exactNativeValue[nativeCharts[[i, "scale_exact"]]] ===
        charts[[i, "Scale"]] &&
      exactNativeValue[nativeCharts[[i, "radius_exact"]]] ===
        charts[[i, "LocalRadius"]], {i, Length[charts]}];
  matchesQ = And @@ Table[Module[{point = exactNativeValue[
        matches[[i, "physical_exact"]]]},
      matches[[i, "producing_chart_index"]] === i - 1 &&
      matches[[i, "receiving_chart_index"]] === i &&
      insideBothDisksQ[point, charts[[i]]] &&
      insideBothDisksQ[point, charts[[i + 1]]]],
    {i, Length[matches]}];
  tilesQ = And @@ Table[Module[{begin = exactNativeValue[
        tiles[[i, "physical_begin_exact"]]], end = exactNativeValue[
        tiles[[i, "physical_end_exact"]]]},
      tiles[[i, "chart_index"]] === i - 1 &&
      insideBothDisksQ[begin, charts[[i]]] &&
      insideBothDisksQ[end, charts[[i]]]],
    {i, Length[tiles]}];
  chartBindingsQ && matchesQ && tilesQ &&
    exactNativeValue[First[tiles]["physical_begin_exact"]] ===
      plan["From"] &&
    exactNativeValue[Last[tiles]["physical_end_exact"]] === plan["To"]];

runDomain[exactDomain_] := Module[
  {atlas, stats, normalized, certificate, match, matchPoint, audit,
   geometryQ, released, after, session, domain, checks},
  atlas = Block[{DiffExp2`Solve`Private`$cppExactDomain = exactDomain},
    catchDE2[
      DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
        system, {DiffExp2`EpsSeries`ESNew[0, {1}]},
        lowerPlan, upperPlan, "Threads" -> 2]]];
  If[FailureQ[atlas] || !AssociationQ[atlas],
    Return[<|"OK" -> False, "Atlas" -> atlas|>, Module]];
  session = atlas["Session"];
  domain = If[TrueQ[exactDomain], "rational", "acb"];
  stats = DiffExp2`CppBackend`PersistentTilePlanStatistics[atlas["Plan"]];
  normalized = Last[atlas["Upper", "Plan", "Charts"]];
  certificate = Lookup[normalized, "NativeRationalScaleBridge", <||>];
  audit = Lookup[atlas, "NativeGeometryAudit", <||>];
  geometryQ = AssociationQ[stats] &&
    exactArmGeometryQ[stats["lower"], atlas["Lower", "Plan"]] &&
    exactArmGeometryQ[stats["upper"], atlas["Upper", "Plan"]];
  match = If[AssociationQ[stats] &&
      ListQ[Lookup[Lookup[stats, "upper", <||>], "matches", None]],
    First[stats["upper", "matches"]], <||>];
  matchPoint = If[StringQ[Lookup[match, "physical_exact", None]],
    ToExpression[match["physical_exact"], InputForm], Indeterminate];
  checks = {
    atlas["Domain"] === domain,
    normalized["Center"] === endpoint["Center"],
    normalized["Prescriptions"] === endpoint["Prescriptions"],
    normalized["Scale"] === nativeScale,
    normalized["LocalRadius"] === 1,
    normalized["Radius"] === nativeScale,
    normalized["MatchRadius"] === nativeScale,
    !KeyExistsQ[normalized, "IncomingMatchPoint"],
    !KeyExistsQ[normalized, "SymmetricMatch"],
    KeyExistsQ[endpoint, "IncomingMatchPoint"],
    Lookup[certificate, "Schema", None] ===
      "diffexp2-native-inward-rational-scale-v1",
    Lookup[certificate, "FloorBits", None] === 64,
    Lookup[certificate, "OriginalScale", None] === goldenScale,
    Lookup[certificate, "NativeScale", None] === nativeScale,
    Lookup[certificate, "OriginalPhysicalRadius", None] === goldenScale,
    Lookup[certificate, "NativePhysicalRadius", None] === nativeScale,
    TrueQ[FullSimplify[RootReduce[
      0 < normalized["Radius"] < endpoint["Radius"]]]],
    AssociationQ[stats] && Lookup[stats, "status", "error"] === "ok",
    Lookup[stats["upper", "charts"][[2]], "center_exact", None] === "1/2",
    Lookup[stats["upper", "charts"][[2]], "scale_exact", None] ===
      ToString[nativeScale, InputForm],
    Lookup[stats["upper", "charts"][[2]], "radius_exact", None] === "1",
    IntegerQ[matchPoint] || Head[matchPoint] === Rational,
    TrueQ[0 < matchPoint < 1/2],
    TrueQ[Abs[matchPoint] < 24/23],
    TrueQ[Abs[matchPoint - 1/2] < nativeScale],
    TrueQ[FullSimplify[RootReduce[
      Abs[matchPoint - 1/2] < goldenScale]]],
    Lookup[audit, "Schema", None] ===
      "diffexp2-native-exact-geometry-proof-v1",
    TrueQ[Lookup[audit, "BackendExactPlanValidated", False]],
    TrueQ[Lookup[audit, "NativeDisksContainedInOriginal", False]],
    Lookup[audit, "HandoffAuthority", None] === "CppExactTilePlanner",
    Lookup[audit, "LowerBridgeCount", None] === 0,
    Lookup[audit, "UpperBridgeCount", None] === 1,
    geometryQ};
  released = DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    atlas];
  after = Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    session, <||>];
  <|"OK" -> And @@ checks && AssociationQ[released] &&
      Lookup[released, "Failures", {"missing"}] === {} &&
      Lookup[after, "tile_plans", -1] === 0 &&
      Lookup[after, "locals", -1] === 0,
    "Checks" -> checks, "Atlas" -> KeyTake[atlas,
      {"Type", "Domain", "NativeGeometryAudit"}],
    "Stats" -> If[AssociationQ[stats],
      KeyTake[stats, {"status", "lower_matches", "upper_matches",
        "lower_tiles", "upper_tiles"}], stats],
    "Released" -> released, "After" -> after|>];

runExactSharedAnchor[] := Module[
  {atlas, stats, lowerChart, upperChart, certificate, audit, session,
   normalization, released, after, geometryQ, checks},
  atlas = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
    catchDE2[
      DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
        system, {DiffExp2`EpsSeries`ESNew[0, {1}]},
        sharedLowerPlan, sharedUpperPlan, "Threads" -> 2]]];
  If[FailureQ[atlas] || !AssociationQ[atlas],
    Return[<|"OK" -> False, "Atlas" -> atlas|>, Module]];
  session = atlas["Session"];
  stats = DiffExp2`CppBackend`PersistentTilePlanStatistics[atlas["Plan"]];
  lowerChart = First[atlas["Lower", "Plan", "Charts"]];
  upperChart = First[atlas["Upper", "Plan", "Charts"]];
  certificate = Lookup[lowerChart, "NativeRationalScaleBridge", <||>];
  normalization = Lookup[lowerChart, "SharedAnchorNormalization", <||>];
  audit = Lookup[atlas, "NativeGeometryAudit", <||>];
  geometryQ = AssociationQ[stats] &&
    exactArmGeometryQ[stats["lower"], atlas["Lower", "Plan"]] &&
    exactArmGeometryQ[stats["upper"], atlas["Upper", "Plan"]];
  checks = {
    atlas["Domain"] === "rational",
    lowerChart === upperChart,
    lowerChart["Center"] === 0,
    lowerChart["Scale"] === sharedSelectedMatch,
    lowerChart["MatchRadius"] === sharedSelectedMatch,
    IntegerQ[Numerator[lowerChart["Radius"]]],
    IntegerQ[Denominator[lowerChart["Radius"]]],
    IntegerQ[Numerator[lowerChart["LocalRadius"]]],
    IntegerQ[Denominator[lowerChart["LocalRadius"]]],
    lowerChart["LocalRadius"] ===
      Together[lowerChart["Radius"]/sharedSelectedMatch],
    TrueQ[FullSimplify[RootReduce[
      sharedSelectedMatch < lowerChart["Radius"] < sharedUpperRadius]]],
    Lookup[normalization, "Schema", None] ===
      "diffexp2-native-shared-anchor-normalization-v1",
    Lookup[normalization, "IncomingLowerPhysicalRadius", None] ===
      sharedLowerRadius,
    Lookup[normalization, "IncomingUpperPhysicalRadius", None] ===
      sharedUpperRadius,
    Lookup[normalization, "IncomingLowerMatchRadius", None] ===
      sharedLowerMatch,
    Lookup[normalization, "IncomingUpperMatchRadius", None] ===
      sharedUpperMatch,
    Lookup[normalization, "SelectedPhysicalRadius", None] ===
      sharedUpperRadius,
    Lookup[normalization, "SelectedMatchRadius", None] ===
      sharedSelectedMatch,
    Lookup[certificate, "Schema", None] ===
      "diffexp2-native-inward-rational-algebraic-radius-v1",
    Lookup[certificate, "BridgeKind", None] ===
      "ExactAlgebraicPhysicalRadius",
    Lookup[certificate, "OriginalPhysicalRadius", None] ===
      sharedUpperRadius,
    Lookup[certificate, "NativePhysicalRadius", None] ===
      lowerChart["Radius"],
    Lookup[certificate, "OriginalLocalRadius", None] ===
      Together[sharedUpperRadius/sharedSelectedMatch],
    Lookup[certificate, "NativeLocalRadius", None] ===
      lowerChart["LocalRadius"],
    TrueQ[Lookup[certificate, "ExactContainmentProved", False]],
    AssociationQ[stats] && Lookup[stats, "status", "error"] === "ok",
    Lookup[audit, "Schema", None] ===
      "diffexp2-native-exact-geometry-proof-v1",
    TrueQ[Lookup[audit, "BackendExactPlanValidated", False]],
    TrueQ[Lookup[audit, "NativeDisksContainedInOriginal", False]],
    Lookup[audit, "LowerBridgeCount", None] === 1,
    Lookup[audit, "UpperBridgeCount", None] === 1,
    geometryQ};
  released = DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    atlas];
  after = Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    session, <||>];
  <|"OK" -> And @@ checks && AssociationQ[released] &&
      Lookup[released, "Failures", {"missing"}] === {} &&
      Lookup[after, "tile_plans", -1] === 0 &&
      Lookup[after, "locals", -1] === 0,
    "Checks" -> checks,
    "Atlas" -> KeyTake[atlas,
      {"Type", "Domain", "NativeGeometryAudit"}],
    "Stats" -> If[AssociationQ[stats],
      KeyTake[stats, {"status", "lower_matches", "upper_matches",
        "lower_tiles", "upper_tiles"}], stats],
    "Released" -> released, "After" -> after|>];

rationalResult = runDomain[True];
DiffExp2`Solve`ClearSolveCaches[];
acbResult = runDomain[False];
DiffExp2`Solve`ClearSolveCaches[];
sharedAnchorResult = runExactSharedAnchor[];
DiffExp2`Solve`ClearSolveCaches[];

(* A backend error association from tile.plan must become Failure at this
   public seam instead of flowing onward as if it were an opaque plan. *)
emptyTopology = <|"singular_points" -> {}, "boundary_points" -> {},
  "complex_projections" -> {}, "branch_sheets" -> {}|>;
dummyArm = <|"from_exact" -> "0", "to_exact" -> "1",
  "charts" -> {"c:missing"}, "topology" -> emptyTopology|>;
loudFailure = DiffExp2`CppBackend`CreatePersistentTilePlan[
  "missing-session", dummyArm, dummyArm, "status-gate-fixture", 3];
loudData = If[FailureQ[loudFailure], loudFailure[[2]], <||>];
loudStatusQ = FailureQ[loudFailure] &&
  Lookup[loudData, "Operation", None] === "tile.plan" &&
  AssociationQ[Lookup[loudData, "BackendResponse", None]] &&
  Lookup[loudData["BackendResponse"], "status", "ok"] === "error";

ok = exactFloorQ && eligibleQ && singularScopeQ &&
  inexactRadiusBridgeQ && exactRadiusBridgeQ &&
  sharedNormalizationQ && sharedEligibleQ && tightRejectedQ &&
  projectionBridgeQ &&
  TrueQ[rationalResult["OK"]] &&
  TrueQ[acbResult["OK"]] && TrueQ[sharedAnchorResult["OK"]] &&
  loudStatusQ;

If[TrueQ[ok],
  Print["PASS: exact regular algebraic-scale native bridge"],
  Print["FAIL: ", InputForm[<|"ExactFloor" -> exactFloorQ,
    "Eligible" -> eligibleQ, "SingularScope" -> singularScopeQ,
    "InexactRadiusBridge" -> inexactRadiusBridgeQ,
    "ExactRadiusBridge" -> exactRadiusBridgeQ,
    "SharedNormalization" -> sharedNormalizationQ,
    "SharedEligible" -> sharedEligibleQ,
    "TightRejected" -> tightRejectedQ,
    "Rational" -> rationalResult,
    "Acb" -> acbResult, "SharedAnchor" -> sharedAnchorResult,
    "LoudFailure" -> loudFailure,
    "LoudStatus" -> loudStatusQ|>]];
  Exit[1]];
