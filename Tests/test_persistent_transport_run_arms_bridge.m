(* Focused Wolfram bridge smoke for atomic retained two-arm transport. *)

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
  "DivisionOrder" -> 3,
  "Variables" -> {},
  "Verbosity" -> 0}];

x = Global`x; eps = Global`eps;
lambda = 1/2 + eps/3;
system = DiffExp2`LoadSystem[<|
  "Matrix" -> {{lambda/x, 1}, {0, 0}}, "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 0}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 1}];
zero = DiffExp2`EpsSeries`ESNew[-3, {0, 0, 0, 0}];
one = DiffExp2`EpsSeries`ESNew[-3, {0, 0, 0, 1}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one, zero}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrand" -> {{1, 0}, x}]];

lowerBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Lower", "Bases"]], "Columns", {}]];
upperBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Upper", "Bases"]], "Columns", {}]];
arms = <|
  "lower" -> <|"receiving_basis" -> lowerBasis|>,
  "upper" -> <|"receiving_basis" -> upperBasis|>|>;
epsilon = If[FailureQ[atlas], <||>, <|
  "min" -> atlas["Request", "EpsWindow", "Min"],
  "max" -> atlas["Request", "EpsWindow", "CompleteMax"],
  "required_complete_max" -> atlas["TargetCompleteMax"],
  "match_required_complete_max" -> atlas["TargetCompleteMax"]|>];
refinement = <|"relative_tolerance" -> "1e-20", "max_steps" -> 1|>;
sessionStats[] := If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];
publicationStats[stats_Association] := KeyTake[stats, {
  "locals", "matches", "line_results", "transport_states",
  "local_matches", "transport_arm_marches", "pending_local_solves",
  "pending_matches", "pending_transport_states"}];

before = sessionStats[];
badArms = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`RunPersistentTransportArms[
    atlas["Plan"], atlas["Anchor"], KeyDrop[arms, "upper"], epsilon,
    "transport-run-arms-bad", refinement]];
afterBad = sessionStats[];
assert["transport_run_arms_bridge_rejects_nonexact_arm_schema",
  FailureQ[badArms] &&
    SameQ[publicationStats[before], publicationStats[afterBad]]];

result = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`RunPersistentTransportArms[
    atlas["Plan"], atlas["Anchor"], arms, epsilon,
    "transport-run-arms-bridge", refinement]];
afterRun = sessionStats[];
states = If[nativeOKQ[result], Lookup[result, "states", <||>], <||>];
lowerState = Lookup[states, "lower", <||>];
upperState = Lookup[states, "upper", <||>];
lowerStats = If[AssociationQ[lowerState],
  DiffExp2`CppBackend`PersistentTransportArmStatistics[lowerState],
  lowerState];
upperStats = If[AssociationQ[upperState],
  DiffExp2`CppBackend`PersistentTransportArmStatistics[upperState],
  upperState];

assert["transport_run_arms_bridge_returns_only_two_overlapping_states",
  nativeOKQ[result] && Sort[Keys[states]] === {"lower", "upper"} &&
    Lookup[result, "workers", 0] === 2 &&
    Lookup[result, "max_parallel_arms", 0] === 2 &&
    Lookup[result, "worker_overlap", False] === True &&
    Lookup[result, "public_result_tokens", None] ===
      "transport_states_only" &&
    Lookup[result, "dependency_only_final_locals", False] === True &&
    Lookup[lowerState, "matches", -1] === 1 &&
    Lookup[upperState, "matches", -1] === 1 &&
    Lookup[Lookup[lowerState, "final_local", <||>], "public_token", True]
      === False &&
    Lookup[Lookup[upperState, "final_local", <||>], "public_token", True]
      === False && nativeOKQ[lowerStats] && nativeOKQ[upperStats]];

assert["transport_run_arms_bridge_has_two_state_only_publication",
  Lookup[afterRun, "locals", -1] === Lookup[before, "locals", 0] &&
    Lookup[afterRun, "matches", -1] === Lookup[before, "matches", 0] &&
    Lookup[afterRun, "line_results", -1] ===
      Lookup[before, "line_results", 0] &&
    Lookup[afterRun, "transport_states", -1] ===
      Lookup[before, "transport_states", 0] + 2 &&
    Lookup[afterRun, "local_matches", -1] ===
      Lookup[before, "local_matches", 0] + 2 &&
    Lookup[afterRun, "transport_arm_marches", -1] ===
      Lookup[before, "transport_arm_marches", 0] + 2 &&
    Lookup[afterRun, "pending_local_solves", -1] === 0 &&
    Lookup[afterRun, "pending_matches", -1] === 0 &&
    Lookup[afterRun, "pending_transport_states", -1] === 0];

stateReleases = If[nativeOKQ[result],
  DiffExp2`CppBackend`ReleasePersistentTransportArm /@
    {lowerState, upperState}, {result}];
atlasRelease = If[FailureQ[atlas], atlas,
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas]];
afterRelease = sessionStats[];
assert["transport_run_arms_bridge_releases_states_without_final_local_tokens",
  AllTrue[stateReleases, nativeOKQ] && AssociationQ[atlasRelease] &&
    Lookup[atlasRelease, "Failures", {"missing"}] === {} &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "matches", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0 &&
    Lookup[afterRelease, "line_results", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
