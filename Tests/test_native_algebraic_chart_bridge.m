(* A rational affine chart is only a topology surrogate.  Genuine algebraic
   prepared geometry and every physical/local handoff identity must remain
   exact and receive a rigorous Acb specialization. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
  "Variables" -> {}}]];

roc = DiffExp2`Config`CFG["RadiusOfConvergence"];
center = (3 - Sqrt[3])/6;
scale = Sqrt[3]/2;
localRadius = 1/2;
physicalRadius = RootReduce[Abs[scale] localRadius];
matchPoint = RootReduce[center/2];

anchor = <|"Center" -> 0, "Scale" -> 1, "Radius" -> 1,
  "LocalRadius" -> 1, "MatchRadius" -> roc,
  "Singular" -> False, "Prescriptions" -> {}|>;
algebraic = <|"Center" -> center, "Scale" -> scale,
  "Radius" -> physicalRadius, "LocalRadius" -> localRadius,
  "MatchRadius" -> Min[physicalRadius, Abs[scale] roc],
  "Singular" -> False, "Prescriptions" -> {},
  "IncomingMatchPoint" -> matchPoint|>;
plan = <|"From" -> 0, "To" -> 1/3,
  "Charts" -> {anchor, algebraic}|>;

bridged = catchDE2[
  DiffExp2`NativeTransport`Private`bridgeNativeRegularPlanScales[plan]];
planning = If[AssociationQ[bridged],
  bridged["Charts"][[2, "NativeRationalPlanningGeometry"]], $Failed];
certificate = If[AssociationQ[bridged],
  bridged["Charts"][[2, "NativeRationalScaleBridge"]], $Failed];
certified = If[AssociationQ[bridged], catchDE2[
  DiffExp2`NativeTransport`Private`nativeCertifiedArmGeometry[bridged]],
  bridged];
request = If[AssociationQ[bridged], catchDE2[
  DiffExp2`NativeTransport`Private`nativeArmRequest[
    bridged, {"c:anchor", "c:algebraic"}]], bridged];
baseRequest = <|"from_exact" -> "0", "to_exact" -> "1/3",
  "charts" -> {"c:anchor"}, "topology" -> <||>|>;
capturedSingle = Block[{
    DiffExp2`CppBackend`RunRequest = Function[raw, raw]},
  DiffExp2`CppBackend`CreatePersistentArmTilePlan[
    "s:test", request, "algebraic-single-test", 3]];
capturedPair = Block[{
    DiffExp2`CppBackend`RunRequest = Function[raw, raw]},
  DiffExp2`CppBackend`CreatePersistentTilePlan[
    "s:test", baseRequest, request, "algebraic-pair-test", 3]];
capturedPairRaw = If[FailureQ[capturedPair],
  Lookup[capturedPair[[2]], "BackendResponse", capturedPair], capturedPair];
halfProtocol = KeyDrop[request, "certified_geometry"];
rejectedHalfProtocol = Block[{
    DiffExp2`CppBackend`RunRequest = Function[raw, raw]},
  DiffExp2`CppBackend`CreatePersistentArmTilePlan[
    "s:test", halfProtocol, "algebraic-half-test", 3]];

fakePrepared = <|"Center" -> center,
  "ChartMap" -> <|"Scale" -> scale|>, "Radius" -> localRadius,
  "Prescriptions" -> {}|>;
persistentGeometry = catchDE2[
  DiffExp2`Solve`Private`cppPersistentGeometry[fakePrepared]];

match = If[AssociationQ[certified], First[certified["matches"]], $Failed];
tiles = If[AssociationQ[certified], certified["tiles"], $Failed];
encodedScalarQ[record_] := AssociationQ[record] &&
  StringQ[Lookup[record, "exact", None]] &&
  MatchQ[Lookup[record, "value", None], {_String, _String}] &&
  MemberQ[{-1, 0, 1}, Lookup[record, "sign", None]];

ok = AssociationQ[bridged] && AssociationQ[planning] &&
  AssociationQ[certificate] && AssociationQ[certified] &&
  AssociationQ[request] && AssociationQ[capturedSingle] &&
  AssociationQ[capturedPairRaw] && FailureQ[rejectedHalfProtocol] &&
  AssociationQ[persistentGeometry] &&
  RootReduce[bridged["Charts"][[2, "Center"]] - center] === 0 &&
  RootReduce[bridged["Charts"][[2, "Scale"]] - scale] === 0 &&
  RootReduce[bridged["Charts"][[2, "Radius"]] - physicalRadius] === 0 &&
  AllTrue[Lookup[planning,
    {"Center", "Scale", "Radius", "LocalRadius", "MatchRadius"}],
    IntegerQ[#] || Head[#] === Rational &] &&
  TrueQ[RootReduce[
    certificate["CenterDisplacement"] + planning["Radius"] <
      physicalRadius]] &&
  TrueQ[certificate["ExactContainmentProved"]] &&
  Lookup[bridged, "NativeHandoffAuthority", None] ===
    "CppExactTilePlanner" &&
  !KeyExistsQ[bridged["Charts"][[2]], "IncomingMatchPoint"] &&
  RootReduce[bridged["Charts"][[2,
      "NativeCertifiedIncomingMatchPoint"]] - matchPoint] === 0 &&
  Lookup[certified, "schema", None] ===
    "diffexp2-wolfram-certified-algebraic-arm-v1" &&
  Length[certified["charts"]] === 2 && Length[certified["matches"]] === 1 &&
  Length[tiles] === 2 &&
  AllTrue[Lookup[match,
    {"physical", "producing_local", "receiving_local"}], encodedScalarQ] &&
  AllTrue[Flatten[Lookup[tiles,
      {"physical_begin", "physical_end", "local_begin", "local_end"}]],
    encodedScalarQ] &&
  request["planning_charts"][[2, "center_exact"]] ===
    ToString[planning["Center"], InputForm] &&
  request["certified_geometry", "exact_identity"] ===
    certified["exact_identity"] &&
  capturedSingle["op"] === "tile.plan_arm" &&
  capturedSingle["arm"] === request &&
  capturedPairRaw["op"] === "tile.plan" &&
  capturedPairRaw["lower"] === baseRequest &&
  capturedPairRaw["upper"] === request &&
  persistentGeometry["center_exact"] === ToString[RootReduce[center], InputForm] &&
  persistentGeometry["scale_exact"] === ToString[RootReduce[scale], InputForm] &&
  MatchQ[persistentGeometry["center_numeric"], {_String, _String}] &&
  MatchQ[persistentGeometry["scale_numeric"], {_String, _String}] &&
  MatchQ[persistentGeometry["radius_numeric"], {_String, _String}];

Print[If[ok, "PASS", "FAIL"],
  ": algebraic charts retain exact geometry behind rational topology surrogates"];
If[!ok, Print[InputForm[{
  "Bridged" -> bridged, "Planning" -> planning,
  "Certificate" -> certificate, "Certified" -> certified,
  "Request" -> request, "CapturedSingle" -> capturedSingle,
  "CapturedPair" -> capturedPair, "CapturedPairRaw" -> capturedPairRaw,
  "RejectedHalfProtocol" -> rejectedHalfProtocol,
  "PersistentGeometry" -> persistentGeometry}]]];
Exit[If[ok, 0, 1]];
