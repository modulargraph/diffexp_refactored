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
    "Integrands" -> {{{1 + Global`eps},
      {1/Global`eps + 1}}, x},
    "TargetCompleteMax" -> 0]];

epsilon = <|"Min" -> -1, "Max" -> 0,
  "RequiredCompleteMax" -> -1|>;
observable[operation_, identity_, coefficients_:{1 + Global`eps}] := <|
  "Operation" -> operation, "Identity" -> identity,
  "CheckpointIdentity" -> identity <> ":checkpoint",
  "CoefficientVector" -> coefficients, "Epsilon" -> epsilon|>;
observables = {
  Append[observable["integrate", "integral"], "TailPolicy" -> "require"],
  observable["limitLower", "lower-limit"],
  observable["limitUpper", "upper-limit"],
  Append[observable["integrate", "polar-integral",
      {1/Global`eps + 1}],
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
    (* The polar integral reserves one source order for its possible
       eps^-1 primitive before the retained arm state is contracted. *)
    run["States", "lower", "epsilon", "required_complete_max"] === 1 &&
    Sort[Keys[Lookup[run, "States", <||>]]] === {"lower", "upper"}];

checkpointPath = FileNameJoin[{$TemporaryDirectory,
  "de2-native-observable-batch-" <> ToString[$ProcessID] <> ".checkpoint"}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
checkpointManifest = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`SaveNativeTransportObservableBatchCheckpoint[
    run, checkpointPath, "native-observable-batch-roundtrip"]], run];
assert["native_observable_batch_checkpoint_is_atomic_and_manifest_bound",
  AssociationQ[checkpointManifest] && FileExistsQ[checkpointPath] &&
    Lookup[checkpointManifest, "Schema", None] ===
      "DiffExp2.NativeTransportObservableCheckpoint/v2" &&
    Lookup[checkpointManifest, "TransportArmMarches", -1] === 2 &&
    Lookup[Lookup[checkpointManifest, "Results", {}], "Identity"] ===
      Lookup[results, "Identity"] &&
    AllTrue[Values[Lookup[checkpointManifest, "StateHandles", <||>]],
      AssociationQ[#] && Sort[Keys[#]] ===
        Sort[{"Handle", "CheckpointIdentity", "ProvenanceSHA256"}] &&
        StringMatchQ[# ["ProvenanceSHA256"],
          RegularExpression["[0-9a-f]{64}"]] &] &&
    DuplicateFreeQ[Lookup[checkpointManifest["Results"],
      "ProvenanceSHA256"]] &&
    AllTrue[Lookup[checkpointManifest["Results"], "ProvenanceSHA256"],
      StringMatchQ[#, RegularExpression["[0-9a-f]{64}"]] &]];

(* Preserve one exact v1 fixture before releasing the live session.  V2 is
   what new sidecars persist; the legacy manifest remains loadable against
   the same schema-2 native checkpoint. *)
legacySchema = "DiffExp2.NativeTransportObservableCheckpoint/v1";
legacyStateHandles = AssociationMap[Function[side, <|
    "Handle" -> run["States", side, "transport_state"],
    "CheckpointIdentity" -> run["States", side, "checkpoint_identity"],
    "ProvenanceIdentity" -> run["States", side, "provenance_identity"]|>],
  {"lower", "upper"}];
legacyCore = <|"Schema" -> legacySchema,
  "Path" -> checkpointManifest["Path"],
  "CheckpointIdentity" -> checkpointManifest["CheckpointIdentity"],
  "TransportArmMarches" -> checkpointManifest["TransportArmMarches"],
  "StateHandles" -> legacyStateHandles,
  "Results" ->
    (DiffExp2`NativeTransport`Private`nativeObservableCheckpointResult[
        #, run["Atlas", "Session"], legacySchema] & /@ results)|>;
legacyManifest = Append[legacyCore, "ManifestIdentity" ->
  DiffExp2`NativeTransport`Private`nativeCheckpointIdentity[
    "de2-native-observable-checkpoint-manifest-", legacyCore]];
assert["compact checkpoint manifest replaces recursive provenance with digests",
  ByteCount[checkpointManifest] < ByteCount[legacyManifest]/10 &&
    AllTrue[Join[Keys /@ Values[checkpointManifest["StateHandles"]],
        Keys /@ checkpointManifest["Results"]],
      FreeQ[#, "ProvenanceIdentity"] &]];

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
integralCertification = If[Length[exportedResults] === 4,
  KeyTake[exportedResults[[1]],
    {"Scope", "ErrorGuarantee", "ErrorEnvelope"}], None];
polarCertification = If[Length[exportedResults] === 4,
  KeyTake[exportedResults[[4]],
    {"Scope", "ErrorGuarantee", "ErrorEnvelope"}], None];
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
assert["native_observable_batch_exports_compact_line_certification_only_for_integrals",
  AssociationQ[integralCertification] &&
    integralCertification["Scope"] ===
      "full_local_with_certified_tail" &&
    integralCertification["ErrorGuarantee"] === "certified" &&
    AssociationQ[integralCertification["ErrorEnvelope"]] &&
    integralCertification["ErrorEnvelope", "guarantee"] ===
      "certified" &&
    ListQ[integralCertification["ErrorEnvelope",
      "absolute_upper_approx"]] &&
    integralCertification["ErrorEnvelope",
      "absolute_upper_approx"] =!= {} &&
    polarCertification === <|"Scope" -> "stored_truncation",
      "ErrorGuarantee" -> "none", "ErrorEnvelope" -> None|> &&
    AllTrue[exportedResults[[{2, 3}]],
      Intersection[Keys[#],
        {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] === {} &]];
malformedCertification = catchDE2[
  DiffExp2`NativeTransport`Private`nativeLineExportCertification[<|
    "scope" -> "full_local_with_certified_tail",
    "error_guarantee" -> "certified",
    "value" -> <|"min" -> 0, "max" -> 0,
      "error" -> <|"min" -> 0, "max" -> 0,
        "guarantee" -> "advisory", "absolute_upper_approx" -> {1.},
        "bound_encoding" -> "approximate-double",
        "provenance" -> "tampered"|>|>|>]];
assert["native_observable_batch_rejects_malformed_certification_metadata",
  FailureQ[malformedCertification]];

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

tamperedResults = checkpointManifest["Results"];
tamperedResults[[{1, 4}, "Handle"]] =
  Reverse[tamperedResults[[{1, 4}, "Handle"]]];
tamperedCore = Join[KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"Results" -> tamperedResults|>];
tamperedManifest = Append[tamperedCore, "ManifestIdentity" ->
  ("de2-native-observable-checkpoint-manifest-" <>
    IntegerString[Hash[tamperedCore, "SHA256"], 16, 64])];
tamperedRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedManifest]];
assert["checkpoint_restore_rejects_swapped_per_handle_observable_identity",
  FailureQ[tamperedRestore]];

tamperedStates = checkpointManifest["StateHandles"];
tamperedStates["lower", "CheckpointIdentity"] =
  "tampered-lower-state-checkpoint";
tamperedStateCore = Join[KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"StateHandles" -> tamperedStates|>];
tamperedStateManifest = Append[tamperedStateCore, "ManifestIdentity" ->
  ("de2-native-observable-checkpoint-manifest-" <>
    IntegerString[Hash[tamperedStateCore, "SHA256"], 16, 64])];
tamperedStateRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedStateManifest]];
assert["checkpoint_restore_binds_each_transport_state_checkpoint_and_provenance",
  FailureQ[tamperedStateRestore]];

tamperedProvenanceResults = checkpointManifest["Results"];
tamperedProvenanceResults[[1, "ProvenanceSHA256"]] =
  StringRepeat["0", 64];
tamperedProvenanceCore = Join[
  KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"Results" -> tamperedProvenanceResults|>];
tamperedProvenanceManifest = Append[tamperedProvenanceCore,
  "ManifestIdentity" ->
    DiffExp2`NativeTransport`Private`nativeCheckpointIdentity[
      "de2-native-observable-checkpoint-manifest-",
      tamperedProvenanceCore]];
tamperedProvenanceRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedProvenanceManifest]];
assert["checkpoint_restore_rejects_recomputed_manifest_with_wrong_provenance_digest",
  FailureQ[tamperedProvenanceRestore]];

restoredRun = If[AssociationQ[checkpointManifest], catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    checkpointManifest]], checkpointManifest];
restoredSession = If[AssociationQ[restoredRun],
  Lookup[restoredRun, "RestoredSession", None], None];
restoredStatsBefore = If[StringQ[restoredSession],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    restoredSession, <||>], <||>];
restoredExport = If[AssociationQ[restoredRun], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    restoredRun, 50]], restoredRun];
restoredStatsAfter = If[StringQ[restoredSession],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    restoredSession, <||>], <||>];
assert["checkpoint_restore_export_does_not_repeat_transport_arm_marches",
  AssociationQ[restoredExport] &&
    Lookup[restoredStatsBefore, "transport_arm_marches", -1] === 2 &&
    Lookup[restoredStatsAfter, "transport_arm_marches", -2] ===
      Lookup[restoredStatsBefore, "transport_arm_marches", -1] &&
    Lookup[Lookup[restoredExport, "ExportedResults", {}], "Value"] ===
      Lookup[exportedResults, "Value"] &&
    (KeyTake[#, {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] & /@
      Select[restoredExport["ExportedResults"],
        #["Operation"] === "integrate" &]) ===
      (KeyTake[#, {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] & /@
        Select[exportedResults, #["Operation"] === "integrate" &])];
restoredReleased = If[AssociationQ[restoredExport],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    restoredExport], restoredExport];
assert["restored_native_observable_batch_closes_only_its_restored_session",
  AssociationQ[restoredReleased] &&
    Lookup[restoredReleased, "Failures", {"missing"}] === {} &&
    !KeyExistsQ[DiffExp2`CppBackend`PersistentSessionInformation[],
      restoredSession]];

legacyRestored = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    legacyManifest]];
legacyExport = If[AssociationQ[legacyRestored], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    legacyRestored, 50]], legacyRestored];
assert["checkpoint_restore_remains_backward_compatible_with_v1_full_provenance_manifest",
  AssociationQ[legacyExport] &&
    Lookup[Lookup[legacyExport, "ExportedResults", {}], "Value"] ===
      Lookup[exportedResults, "Value"]];
legacyReleased = If[AssociationQ[legacyExport],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    legacyExport], legacyExport];
assert["legacy_manifest_restore_closes_its_restored_session",
  AssociationQ[legacyReleased] &&
    Lookup[legacyReleased, "Failures", {"missing"}] === {}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
