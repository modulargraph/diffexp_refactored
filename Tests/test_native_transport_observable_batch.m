(* Focused level-facing retained transport-observable batch smoke. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];
nativeOKQ[value_] := AssociationQ[value] &&
  Lookup[value, "status", "error"] === "ok";

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 1, "DivisionOrder" -> 2,
  "Variables" -> {}, "Verbosity" -> 0}];

x = Global`x;
system = DiffExp2`LoadSystem[<|"Matrix" -> {{0}}, "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {0, -1/4}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {0, 1/4}];
one = DiffExp2`EpsSeries`ESNew[0, {1, 0}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrands" -> {{{1}, {1/Global`eps}}, x},
    "TargetCompleteMax" -> 0]];

epsilon = <|"Min" -> -1, "Max" -> 0,
  "RequiredCompleteMax" -> -1|>;
observable[operation_, identity_, coefficients_:{1}] := <|
  "Operation" -> operation, "Identity" -> identity,
  "CheckpointIdentity" -> identity <> ":checkpoint",
  "CoefficientVector" -> coefficients, "Epsilon" -> epsilon|>;
observables = {
  Append[observable["integrate", "integral"], "TailPolicy" -> "stored"],
  observable["limitLower", "lower-limit"],
  observable["limitUpper", "upper-limit"],
  Append[observable["integrate", "polar-integral", {1/Global`eps}],
    "TailPolicy" -> "stored"]};
sessionStats[] := If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];

beforeInvalid = sessionStats[];
duplicate = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, {observable["limitLower", "duplicate"],
      observable["limitUpper", "duplicate"]}, x]]];
afterInvalid = sessionStats[];
assert["native_observable_batch_rejects_duplicate_identity_before_march",
  FailureQ[duplicate] &&
    Lookup[beforeInvalid, "transport_states", -1] ===
      Lookup[afterInvalid, "transport_states", -2] &&
    Lookup[beforeInvalid, "transport_arm_marches", -1] ===
      Lookup[afterInvalid, "transport_arm_marches", -2]];

run = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, observables, x, "MaxRefinementSteps" -> 1]]];
results = If[AssociationQ[run], Lookup[run, "Results", {}], {}];
assert["native_observable_batch_marches_once_and_preserves_request_order",
  AssociationQ[run] &&
    Lookup[run, "Type", None] ===
      "DiffExp2NativeTransportObservableBatch" &&
    Lookup[run, "NativeMarches", 0] === 2 &&
    Lookup[results, "RequestIndex"] === {0, 1, 2, 3} &&
    Lookup[results, "Identity"] ===
      {"integral", "lower-limit", "upper-limit", "polar-integral"} &&
    atlas["TargetCompleteMax"] === 0 &&
    atlas["Request", "EpsWindow", "CompleteMax"] === 1 &&
    run["States", "lower", "epsilon", "required_complete_max"] === 0 &&
    Sort[Keys[Lookup[run, "States", <||>]]] === {"lower", "upper"}];

exportedRun = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[run, 50]],
  run];
exportedResults = If[AssociationQ[exportedRun],
  Lookup[exportedRun, "ExportedResults", {}], {}];
lineValue = If[Length[exportedResults] === 4,
  exportedResults[[1, "Value"]], None];
lowerValue = If[Length[exportedResults] === 4,
  exportedResults[[2, "Value"]], None];
upperValue = If[Length[exportedResults] === 4,
  exportedResults[[3, "Value"]], None];
polarValue = If[Length[exportedResults] === 4,
  exportedResults[[4, "Value"]], None];
assert["native_observable_batch_contracts_integrals_polar_order_and_endpoints",
  AssociationQ[exportedRun] &&
    Lookup[exportedRun, "CompatibilityExports", 0] === 4 &&
    AllTrue[{lineValue, lowerValue, upperValue, polarValue},
      DiffExp2`EpsSeries`ESQ] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[lineValue, 0] - 1/2,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[lowerValue, 0] - 1,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[upperValue, 0] - 1,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[polarValue, -1] - 1/2,
      30]] < 10^-20]];

released = If[AssociationQ[exportedRun],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[exportedRun],
  exportedRun];
afterRelease = sessionStats[];
assert["native_observable_batch_releases_complete_owner_closure",
  AssociationQ[released] && Lookup[released, "Failures", {"missing"}] === {} &&
    Lookup[afterRelease, "line_results", -1] === 0 &&
    Lookup[afterRelease, "endpoints", -1] === 0 &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "matches", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
