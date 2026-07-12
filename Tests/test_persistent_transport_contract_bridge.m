(* Focused Wolfram bridge coverage for retained-arm 0/1/N contraction. *)

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
    transportEpsilon, "transport-contract-bridge-state",
    <|"relative_tolerance" -> "1e-20", "max_steps" -> 1|>]];

tileCount = If[nativeOKQ[state], Lookup[state, "tiles", 0], 0];
epsilon = <|"Min" -> 0, "Max" -> 0,
  "RequiredCompleteMax" -> 0|>;
observable[identity_, checkpoint_] := <|
  "Identity" -> identity,
  "CheckpointIdentity" -> checkpoint,
  "IntegrandRows" -> rows,
  "Epsilon" -> epsilon,
  "TailPolicy" -> "stored"|>;
sessionStats[] := If[!nativeOKQ[state], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    state["session"], <||>]];
mutationCounters[stats_Association] := KeyTake[stats, {
  "transport_contractions", "transport_observables", "line_results",
  "line_integrations", "local_matches", "pending_local_solves",
  "pending_matches", "pending_line_integrations"}];

before = sessionStats[];
missingProvenance = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[
    KeyDrop[state, "provenance_identity"], {}, "missing-provenance"],
  state];
duplicateIdentity = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[state, {
    observable["duplicate", "duplicate-checkpoint-1"],
    observable["duplicate", "duplicate-checkpoint-2"]},
    "duplicate-identity"], state];
badRows = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[state, {
    Join[observable["bad-rows", "bad-rows-checkpoint"],
      <|"IntegrandRows" -> {}|>]}, "bad-rows"], state];
badTail = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[state, {
    Join[observable["bad-tail", "bad-tail-checkpoint"],
      <|"TailPolicy" -> "certified"|>]}, "bad-tail"], state];
afterInvalid = sessionStats[];

assert["transport_contract_smallest_retained_arm_prepared",
  nativeOKQ[state] && tileCount === 1 && Length[rows] === 1];
assert["transport_contract_bridge_rejects_bad_records_before_native_call",
  FailureQ[missingProvenance] && FailureQ[duplicateIdentity] &&
    FailureQ[badRows] && FailureQ[badTail] &&
    SameQ[mutationCounters[before], mutationCounters[afterInvalid]]];

zero = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[
    state, {}, "transport-contract-zero"], state];
afterZero = sessionStats[];
assert["transport_contract_zero_observables_is_allocation_free",
  nativeOKQ[zero] && Lookup[zero, "observables", -1] === 0 &&
    Lookup[zero, "lines", Missing["Absent"]] === {} &&
    Lookup[afterZero, "transport_contractions", -1] ===
      Lookup[before, "transport_contractions", 0] + 1 &&
    Lookup[afterZero, "transport_observables", -1] ===
      Lookup[before, "transport_observables", 0] &&
    Lookup[afterZero, "line_results", -1] ===
      Lookup[before, "line_results", 0] &&
    Lookup[afterZero, "local_matches", -1] ===
      Lookup[before, "local_matches", 0]];

oneResult = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[state, {
    observable["one-observable", "one-observable-checkpoint"]},
    "transport-contract-one"], state];
oneLines = If[nativeOKQ[oneResult], Lookup[oneResult, "lines", {}], {}];
oneStats = If[Length[oneLines] === 1,
  DiffExp2`CppBackend`PersistentLineIntegralStatistics[First[oneLines]],
  <||>];
assert["transport_contract_one_returns_directly_usable_opaque_line",
  nativeOKQ[oneResult] && Length[oneLines] === 1 &&
    Lookup[First[oneLines], "session", None] === state["session"] &&
    nativeOKQ[oneStats] &&
    Lookup[oneResult, "no_rematching", False] === True &&
    Lookup[oneResult, "json_coefficients", -1] === 0];

manyIdentities = {"observable-z", "observable-a"};
manyResult = If[nativeOKQ[state],
  DiffExp2`CppBackend`ContractPersistentTransportObservables[state, {
    observable[manyIdentities[[1]], "many-checkpoint-1"],
    observable[manyIdentities[[2]], "many-checkpoint-2"]},
    "transport-contract-many"], state];
manyLines = If[nativeOKQ[manyResult], Lookup[manyResult, "lines", {}], {}];
afterMany = sessionStats[];
assert["transport_contract_many_preserves_request_order_without_rematching",
  nativeOKQ[manyResult] && Length[manyLines] === 2 &&
    Lookup[manyLines, "request_index"] === {0, 1} &&
    Lookup[manyLines, "observable_identity"] === manyIdentities &&
    DuplicateFreeQ[Lookup[manyLines, "line"]] &&
    Lookup[afterMany, "transport_contractions", -1] ===
      Lookup[before, "transport_contractions", 0] + 3 &&
    Lookup[afterMany, "transport_observables", -1] ===
      Lookup[before, "transport_observables", 0] + 3 &&
    Lookup[afterMany, "line_integrations", -1] ===
      Lookup[before, "line_integrations", 0] + 3 tileCount &&
    Lookup[afterMany, "local_matches", -1] ===
      Lookup[before, "local_matches", 0]];

lineHandles = Join[oneLines, manyLines];
lineReleases = DiffExp2`CppBackend`ReleasePersistentLineIntegral /@
  lineHandles;
finalLocal = If[nativeOKQ[state], Lookup[state, "final_local", <||>], <||>];
anchorLocalHandle = If[FailureQ[atlas], None,
  Lookup[atlas["Anchor"], "local",
    Lookup[atlas["Anchor"], "Local", None]]];
finalLocalRelease = If[AssociationQ[finalLocal] &&
    Lookup[finalLocal, "local",
      Lookup[finalLocal, "Local", None]] =!= anchorLocalHandle,
  DiffExp2`CppBackend`ReleasePersistentLocal[finalLocal], None];
stateRelease = If[nativeOKQ[state],
  DiffExp2`CppBackend`ReleasePersistentTransportArm[state], state];
atlasRelease = If[FailureQ[atlas], atlas,
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas]];
afterRelease = If[nativeOKQ[state],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    state["session"], <||>], <||>];

assert["transport_contract_bridge_releases_public_results_and_state",
  AllTrue[lineReleases, nativeOKQ] &&
    (finalLocalRelease === None || nativeOKQ[finalLocalRelease]) &&
    nativeOKQ[stateRelease] &&
    AssociationQ[atlasRelease] && atlasRelease["Failures"] === {} &&
    Lookup[afterRelease, "line_results", -1] === 0 &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
