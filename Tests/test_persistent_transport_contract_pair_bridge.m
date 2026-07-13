(* Focused Wolfram bridge coverage for retained paired-arm contraction. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label,
    If[detail === None, "", ": " <> ToString[detail, InputForm]]]];
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
(* Pair integration consumes one primitive epsilon halo even for this
   epsilon-independent fixture.  Retain that certified zero coefficient so
   the requested epsilon^0 result remains globally complete. *)
one = DiffExp2`EpsSeries`ESNew[-2, {0, 0, 1, 0}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrand" -> {{1}, x}, "TargetCompleteMax" -> 1]];

lowerRows = If[FailureQ[atlas], {}, catchDE2[
  DiffExp2`NativeTransport`Private`nativePreparedArmRows[
    atlas, atlas["Lower"], {1}, x]]];
upperRows = If[FailureQ[atlas], {}, catchDE2[
  DiffExp2`NativeTransport`Private`nativePreparedArmRows[
    atlas, atlas["Upper"], {1}, x]]];
lowerBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Lower", "Bases"]], "Columns", {}]];
upperBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Upper", "Bases"]], "Columns", {}]];
transportEpsilon = If[FailureQ[atlas], <||>, <|
  "min" -> atlas["Request", "EpsWindow", "Min"],
  "max" -> atlas["Request", "EpsWindow", "CompleteMax"],
  "required_complete_max" -> atlas["TargetCompleteMax"],
  "match_required_complete_max" -> atlas["TargetCompleteMax"]|>];
transport = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`RunPersistentTransportArms[
    atlas["Plan"], atlas["Anchor"], <|
      "lower" -> <|"receiving_basis" -> lowerBasis|>,
      "upper" -> <|"receiving_basis" -> upperBasis|>|>,
    transportEpsilon, "transport-pair-bridge-states",
    <|"relative_tolerance" -> "1e-20", "max_steps" -> 1|>]];
states = If[nativeOKQ[transport], Lookup[transport, "states", <||>], <||>];
lowerState = Lookup[states, "lower", <||>];
upperState = Lookup[states, "upper", <||>];

epsilon = <|"Min" -> 0, "Max" -> 0,
  "RequiredCompleteMax" -> 0|>;
boundedCancellation = <|
  "Mode" -> "bounded-relative-acb",
  "RelativeTolerance" ->
    "1.0000000000000000000000000000000000000000000000000e-24",
  "Provenance" -> ExportString[<|
      "schema" -> "feynman-trick-divergent-cancellation-v1",
      "producer" -> "DiffExp2`Tolerances`Tol",
      "key" -> "LaurentLeadTol",
      "exact_value" -> "1/1000000000000000000000000"|>,
    "RawJSON", "Compact" -> True]|>;
observable[identity_, checkpoint_, tail_:Automatic] := Module[{result},
  result = <|"Identity" -> identity,
    "CheckpointIdentity" -> checkpoint,
    "LowerIntegrandRows" -> lowerRows,
    "UpperIntegrandRows" -> upperRows,
    "Epsilon" -> epsilon|>;
  If[tail === Automatic, result,
    Append[result, "TailPolicy" -> tail]]];
sessionStats[] := If[!nativeOKQ[transport], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    transport["session"], <||>]];
mutationCounters[stats_Association] := KeyTake[stats, {
  "locals", "matches", "line_results", "line_integrations",
  "transport_contractions", "transport_observables",
  "transport_pair_contractions", "transport_pair_observables",
  "pending_line_integrations"}];

beforeInvalid = sessionStats[];
extraSign = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {Append[observable["extra-sign", "extra-sign-checkpoint"],
      "Signs" -> {-1, 1}]}, "extra-sign-root"], transport];
missingUpper = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {KeyDrop[observable["missing-upper", "missing-upper-checkpoint"],
      "UpperIntegrandRows"]}, "missing-upper-root"], transport];
duplicateIdentity = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {observable["duplicate", "duplicate-checkpoint-1"],
     observable["duplicate", "duplicate-checkpoint-2"]},
    "duplicate-identity-root"], transport];
duplicateCheckpoint = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {observable["duplicate-a", "duplicate-checkpoint"],
     observable["duplicate-b", "duplicate-checkpoint"]},
    "duplicate-checkpoint-root"], transport];
overRequired = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {Join[observable["over-required", "over-required-checkpoint"], <|
      "Epsilon" -> <|"Min" -> 0, "Max" -> 2,
        "RequiredCompleteMax" -> 2|>|>]}, "over-required-root"],
  transport];
swapped = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    upperState, lowerState, {}, "swapped-root"], transport];
tampered = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    Join[lowerState, <|"provenance_identity" -> "tampered"|>],
    upperState, {}, "tampered-root"], transport];
afterInvalid = sessionStats[];

assert["transport_pair_bridge_rejects_malformed_signs_and_duplicates",
  FailureQ[extraSign] && FailureQ[missingUpper] &&
    FailureQ[duplicateIdentity] && FailureQ[duplicateCheckpoint] &&
    (FailureQ[overRequired] ||
      Lookup[overRequired, "status", None] === "error") &&
    FailureQ[swapped]];
assert["transport_pair_bridge_rejects_native_provenance_tamper_atomically",
  AssociationQ[tampered] && Lookup[tampered, "status", "ok"] === "error" &&
    SameQ[mutationCounters[beforeInvalid],
      mutationCounters[afterInvalid]]];

beforeZero = sessionStats[];
zero = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState, {}, "transport-pair-zero"], transport];
afterZero = sessionStats[];
assert["transport_pair_bridge_zero_is_line_allocation_free",
  nativeOKQ[zero] && Lookup[zero, "observables", -1] === 0 &&
    Lookup[zero, "lines", Missing["Absent"]] === {} &&
    Lookup[zero, "combination", None] === "negative-lower-plus-upper" &&
    Lookup[zero, "max_parallel_arms", -1] === 0 &&
    Lookup[afterZero, "line_results", -1] ===
      Lookup[beforeZero, "line_results", 0] &&
    Lookup[afterZero, "line_integrations", -1] ===
      Lookup[beforeZero, "line_integrations", 0] &&
    Lookup[afterZero, "transport_pair_contractions", -1] ===
      Lookup[beforeZero, "transport_pair_contractions", 0] + 1 &&
    Lookup[afterZero, "transport_contractions", -1] ===
      Lookup[beforeZero, "transport_contractions", 0] + 2];

oneResult = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {Append[observable["one-observable", "one-observable-checkpoint"],
      "DivergentCancellation" -> boundedCancellation]},
    "transport-pair-one"], transport];
oneLines = If[nativeOKQ[oneResult], Lookup[oneResult, "lines", {}], {}];
oneLineStats = If[Length[oneLines] === 1,
  DiffExp2`CppBackend`PersistentLineIntegralStatistics[First[oneLines]],
  <||>];
oneExport = If[Length[oneLines] === 1,
  DiffExp2`CppBackend`ExportPersistentLineIntegral[
    First[oneLines], First[oneLines]["checkpoint_identity"], 50],
  <||>];
oneDecoded = If[nativeOKQ[oneExport],
  DiffExp2`CppBackend`DecodeScalars[
    oneExport["value", "coefficients"], 50], {}];
oneEpsilonMin = If[nativeOKQ[oneExport],
  oneExport["value", "min"], Missing["NoMinimum"]];
oneEpsilonZero = If[ListQ[oneDecoded] && IntegerQ[oneEpsilonMin] &&
    1 <= 1 - oneEpsilonMin <= Length[oneDecoded],
  oneDecoded[[1 - oneEpsilonMin]], Missing["NoEpsilonZero"]];
assert["transport_pair_bridge_one_defaults_tail_and_fixes_minus_plus",
  nativeOKQ[oneResult] && Length[oneLines] === 1 &&
    Lookup[oneResult, "combination", None] ===
      "negative-lower-plus-upper" &&
    Lookup[oneResult, "max_parallel_arms", 0] === 1 &&
    Lookup[oneResult, "no_remarching", False] === True &&
    Lookup[oneResult, "no_rematching", False] === True &&
    nativeOKQ[oneLineStats] && nativeOKQ[oneExport] &&
    Lookup[Lookup[oneLineStats, "diagnostics", <||>],
      "divergent_cancellation_mode", None] === "bounded-relative-acb" &&
    Lookup[Lookup[oneLineStats, "diagnostics", <||>],
      "divergent_cancellation_provenance", None] ===
        boundedCancellation["Provenance"] &&
    NumberQ[oneEpsilonZero] &&
    TrueQ[Abs[N[oneEpsilonZero - 1/2, 30]] < 10^-20],
  {KeyTake[oneResult, {"status", "observables", "combination",
      "max_parallel_arms", "no_remarching", "no_rematching"}],
   KeyTake[oneLineStats, {"status"}], KeyTake[oneExport, {"status"}],
   oneEpsilonZero}];

manyIdentities = {"many-z", "many-a"};
manyResult = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ContractPersistentTransportPairObservables[
    lowerState, upperState,
    {observable[manyIdentities[[1]], "many-checkpoint-1", "stored"],
     observable[manyIdentities[[2]], "many-checkpoint-2"]},
    "transport-pair-many"], transport];
manyLines = If[nativeOKQ[manyResult], Lookup[manyResult, "lines", {}], {}];
afterMany = sessionStats[];
assert["transport_pair_bridge_many_preserves_order_and_opaque_outputs",
  nativeOKQ[manyResult] && Length[manyLines] === 2 &&
    Lookup[manyLines, "request_index"] === {0, 1} &&
    Lookup[manyLines, "observable_identity"] === manyIdentities &&
    DuplicateFreeQ[Lookup[manyLines, "line"]] &&
    Lookup[manyResult, "max_parallel_arms", 0] === 1 &&
    Lookup[afterMany, "transport_pair_observables", -1] === 3 &&
    Lookup[afterMany, "pending_line_integrations", -1] === 0,
  {KeyTake[manyResult, {"status", "observables", "combination",
      "max_parallel_arms"}],
   Lookup[manyLines, {"request_index", "observable_identity", "line"}],
   KeyTake[afterMany, {"transport_pair_observables",
      "pending_line_integrations"}]}];

lineHandles = Join[oneLines, manyLines];
lineReleases = DiffExp2`CppBackend`ReleasePersistentLineIntegral /@
  lineHandles;
stateReleases = If[nativeOKQ[transport],
  DiffExp2`CppBackend`ReleasePersistentTransportArm /@
    {lowerState, upperState}, {transport}];
atlasRelease = If[FailureQ[atlas], atlas,
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas]];
afterRelease = sessionStats[];
assert["transport_pair_bridge_releases_lines_states_and_atlas",
  AllTrue[lineReleases, nativeOKQ] &&
    AllTrue[stateReleases, nativeOKQ] && AssociationQ[atlasRelease] &&
    Lookup[atlasRelease, "Failures", {"missing"}] === {} &&
    Lookup[afterRelease, "line_results", -1] === 0 &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "matches", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
