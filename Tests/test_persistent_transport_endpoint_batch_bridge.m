(* Focused Wolfram bridge smoke for retained-arm endpoint batches. *)

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
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30,
  "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0,
  "DivisionOrder" -> 2,
  "Variables" -> {},
  "Verbosity" -> 0}];

x = Global`x;
system = DiffExp2`LoadSystem[<|"Matrix" -> {{0}}, "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {0, -1/4}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {0, 1/4}];
one = DiffExp2`EpsSeries`ESNew[-2, {0, 0, 1}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrand" -> {{1}, x}]];

rows = If[FailureQ[atlas], {}, catchDE2[
  DiffExp2`NativeTransport`Private`nativePreparedArmRows[
    atlas, atlas["Lower"], {1}, x]]];
endpointRow = If[ListQ[rows] && rows =!= {}, Last[rows], <||>];
receivingBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Lower", "Bases"]], "Columns", {}]];
transportEpsilon = If[FailureQ[atlas], <||>, <|
  "min" -> atlas["Request", "EpsWindow", "Min"],
  "max" -> atlas["Request", "EpsWindow", "CompleteMax"],
  "required_complete_max" -> atlas["TargetCompleteMax"],
  "match_required_complete_max" -> atlas["TargetCompleteMax"]|>];
state = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`RunPersistentTransportArm[
    atlas["Plan"], "lower", atlas["Anchor"], receivingBasis,
    transportEpsilon, "transport-endpoint-bridge-state",
    <|"relative_tolerance" -> "1e-20", "max_steps" -> 1|>]];

stateRequired = If[nativeOKQ[state],
  state["epsilon", "required_complete_max"], 0];
epsilon = <|"Min" -> 0, "Max" -> 0,
  "RequiredCompleteMax" -> 0|>;
observable[identity_, checkpoint_] := <|
  "Identity" -> identity,
  "CheckpointIdentity" -> checkpoint,
  "IntegrandRow" -> endpointRow,
  "Epsilon" -> epsilon,
  "PublicationRelativeTolerance" -> "1e-8"|>;
sessionStats[] := If[!nativeOKQ[state], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    state["session"], <||>]];
mutationCounters[stats_Association] := KeyTake[stats, {
  "transport_endpoint_batches", "transport_endpoint_rows", "endpoints",
  "endpoint_limits", "local_matches", "pending_endpoint_limits"}];

before = sessionStats[];
missingProvenance = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
    KeyDrop[state, "provenance_identity"], {}, "missing-provenance"],
  state];
duplicateIdentity = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    observable["duplicate", "duplicate-checkpoint-1"],
    observable["duplicate", "duplicate-checkpoint-2"]},
    "duplicate-identity"], state];
badRow = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    Join[observable["bad-row", "bad-row-checkpoint"], <|
      "IntegrandRow" -> Join[endpointRow, <|"columns" -> 0|>]|>]},
    "bad-row"], state];
forbiddenPoint = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    Append[observable["point", "point-checkpoint"], "Point" -> 0]},
    "forbidden-point"], state];
missingPublicationTolerance = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    KeyDrop[observable["missing-publication",
      "missing-publication-checkpoint"],
      "PublicationRelativeTolerance"]},
    "missing-publication"], state];
overTarget = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    Join[observable["over-target", "over-target-checkpoint"], <|
      "Epsilon" -> <|"Min" -> 0, "Max" -> stateRequired + 1,
        "RequiredCompleteMax" -> stateRequired + 1|>|>]},
    "over-target"], state];
afterInvalid = sessionStats[];

assert["transport_endpoint_bridge_prepares_one_retained_final_row",
  nativeOKQ[state] && AssociationQ[endpointRow]];
assert["transport_endpoint_bridge_rejects_invalid_records_before_native_call",
    FailureQ[missingProvenance] && FailureQ[duplicateIdentity] &&
    FailureQ[badRow] && FailureQ[forbiddenPoint] &&
    FailureQ[missingPublicationTolerance] && FailureQ[overTarget] &&
    SameQ[mutationCounters[before], mutationCounters[afterInvalid]]];

zero = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[
    state, {}, "transport-endpoint-zero"], state];
afterZero = sessionStats[];
assert["transport_endpoint_bridge_zero_is_allocation_free",
  nativeOKQ[zero] && Lookup[zero, "observables", -1] === 0 &&
    Lookup[zero, "endpoints", Missing["Absent"]] === {} &&
    Lookup[afterZero, "transport_endpoint_batches", -1] ===
      Lookup[before, "transport_endpoint_batches", 0] + 1 &&
    Lookup[afterZero, "transport_endpoint_rows", -1] ===
      Lookup[before, "transport_endpoint_rows", 0] &&
    Lookup[afterZero, "endpoints", -1] ===
      Lookup[before, "endpoints", 0]];

oneResult = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    observable["one-observable", "one-observable-checkpoint"]},
    "transport-endpoint-one"], state];
oneEndpoints = If[nativeOKQ[oneResult],
  Lookup[oneResult, "endpoints", {}], {}];
oneStats = If[Length[oneEndpoints] === 1,
  DiffExp2`CppBackend`PersistentEndpointStatistics[First[oneEndpoints]],
  <||>];
assert["transport_endpoint_bridge_one_returns_usable_endpoint",
  nativeOKQ[oneResult] && Length[oneEndpoints] === 1 &&
    Lookup[First[oneEndpoints], "session", None] === state["session"] &&
    nativeOKQ[oneStats] &&
    Lookup[oneResult, "no_projected_local_publication", False] === True &&
    Lookup[oneResult, "json_coefficients", -1] === 0];

manyIdentities = {"observable-z", "observable-a"};
manyResult = If[nativeOKQ[state],
  DiffExp2`CppBackend`RunPersistentTransportEndpointBatch[state, {
    observable[manyIdentities[[1]], "many-checkpoint-1"],
    observable[manyIdentities[[2]], "many-checkpoint-2"]},
    "transport-endpoint-many", 2], state];
manyEndpoints = If[nativeOKQ[manyResult],
  Lookup[manyResult, "endpoints", {}], {}];
afterMany = sessionStats[];
assert["transport_endpoint_bridge_many_preserves_order_and_counters",
  nativeOKQ[manyResult] && Length[manyEndpoints] === 2 &&
    Lookup[manyEndpoints, "request_index"] === {0, 1} &&
    Lookup[manyEndpoints, "observable_identity"] === manyIdentities &&
    DuplicateFreeQ[Lookup[manyEndpoints, "endpoint"]] &&
    Lookup[manyResult, "requested_observable_threads", 0] === 2 &&
    Lookup[manyResult, "observable_worker_threads", 0] === 2 &&
    Lookup[afterMany, "transport_endpoint_batches", -1] ===
      Lookup[before, "transport_endpoint_batches", 0] + 3 &&
    Lookup[afterMany, "transport_endpoint_rows", -1] ===
      Lookup[before, "transport_endpoint_rows", 0] + 3 &&
    Lookup[afterMany, "endpoints", -1] ===
      Lookup[before, "endpoints", 0] + 3 &&
    Lookup[afterMany, "local_matches", -1] ===
      Lookup[before, "local_matches", 0]];

endpointHandles = Join[oneEndpoints, manyEndpoints];
endpointReleases = DiffExp2`CppBackend`ReleasePersistentEndpoint /@
  endpointHandles;
stateRelease = If[nativeOKQ[state],
  DiffExp2`CppBackend`ReleasePersistentTransportArm[state], state];
atlasRelease = If[FailureQ[atlas], atlas,
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas]];
afterRelease = If[nativeOKQ[state],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    state["session"], <||>], <||>];

assert["transport_endpoint_bridge_releases_results_and_state",
  AllTrue[endpointReleases, nativeOKQ] && nativeOKQ[stateRelease] &&
    AssociationQ[atlasRelease] && atlasRelease["Failures"] === {} &&
    Lookup[afterRelease, "endpoints", -1] === 0 &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
