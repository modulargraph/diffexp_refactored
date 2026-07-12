(* Focused production-loader contract for completed schema-2 native FT state.
   This is definitions-only: the persistent batch round-trip itself is covered
   by test_native_transport_observable_batch.m. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
tmpDir = CreateDirectory[FileNameJoin[{$TemporaryDirectory,
  "DiffExp2_ft_native_resume_" <> ToString[$ProcessID]}]];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> "Cpp"];
SetEnvironment["FT_LADDER_CHECKPOINT_DIR" -> tmpDir];
SetEnvironment["FT_RESUME_LADDER_CHECKPOINT" -> ""];
SetEnvironment["FT_ALLOW_STALE_LADDER_CHECKPOINT" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label,
    If[detail === None, "", ": " <> ToString[InputForm[detail]]]]];

name = "native-resume-fixture";
prepKey = 271828;
z = Global`zNativeResume;
mastersHere = {{2, 0}};
mastersBelow = {{1, 1}};
levelData = <|"Masters" -> mastersHere, "FeynmanParameter" -> z,
  "CombinedPositions" -> {1, 2}, "DiffMatrix" -> {{0}}|>;
data = <|"NumLevels" -> 1, "Levels" -> <|
  0 -> <|"Masters" -> mastersBelow|>, 1 -> levelData|>|>;
requests = FeynmanTrick`LevelReduction`BoundaryRequestRecords[
  mastersBelow, levelData["CombinedPositions"]];
batch = <|"Schema" -> "FeynmanTrick.LevelIBPBatch/v1",
  "Key" -> StringRepeat["a", 64],
  "PayloadKey" -> StringRepeat["b", 64],
  "KeyRecord" -> {"native-resume-fixture"}, "UpperLevel" -> 1,
  "MastersAbove" -> mastersHere, "BoundaryRequests" -> requests,
  "CoefficientVectors" ->
    Association[requests[[1, "NeededVec"]] -> {1}]|>;
nativePlan = ft2BuildNativeEpsilonPlan[data, epsOrder,
  levelEpsilonHalos, Identity, <|1 -> batch|>];
requiredRaw = nativePlan["Levels", 1, "RequiredRawTop"];
deepBoundaryFixture = <|
  "BoundaryValues" -> {Join[{1}, ConstantArray[0,
      nativePlan["DeepRequiredRawTop"]]]},
  "EpsPrefactors" -> {0},
  "RequestedEpsilonOrder" -> nativePlan["DeepRequiredRawTop"],
  "WorkingEpsilonOrder" -> nativePlan["DeepRequiredRawTop"]|>;
execution = ft2FinalizeNativeEpsilonPlan[
  nativePlan, deepBoundaryFixture, nativePlan["DeepRequiredRawTop"]];
boundaryValues = {Join[{1}, ConstantArray[0, requiredRaw + 1]]};
boundaryPrefactors = {0};
basis = ft2NormalizeEpsilonBasis[
  {{0}}, boundaryValues, boundaryPrefactors, Global`eps,
  nativePlan["Levels", 1, "Gauge"]];
system = <|"Variable" -> z, "Matrix" -> basis["Matrix"],
  "SingularFactors" -> {}, "SingularFactorsExact" -> {}|>;
entries = ft2PrepareBoundaryEntries[1, batch, boundaryPrefactors,
  z, Global`eps, Identity];
ledger = ft2NativeEpsilonLedger[entries, boundaryValues,
  nativePlan["Levels", 1, "RequiredOutputRawTop"]];
delta = {{z, 1}, {1 - z, 1}};
configuration = <|"WorkingPrecision" -> wp,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> ledger["TargetCompleteMax"],
  "DivisionOrder" -> divisionOrder,
  "RadiusOfConvergence" -> radiusOfConvergence,
  "StepDivisionOrder" -> stepDivisionOrder,
  "RecurrenceBackend" -> "Cpp",
  "SingularMatchPrecondition" -> singularMatchPrecondition,
  "ValueTransportMode" -> Environment["DE2_VALUE_TRANSPORT"],
  "CppThreads" -> cppArmThreadBudget|>;
contract = ft2NativeTransportContract[name, 1, prepKey, system,
  boundaryValues, boundaryPrefactors, entries, ledger, configuration,
  delta, {}, execution["Identity"]];
atlasIdentity = "native-resume-atlas";
payloadIdentity = "native-resume-payload";
audit = <|"Schema" -> "FeynmanTrick.NativeObservableBatch/v1",
  "BatchKey" -> batch["Key"], "BatchPayloadKey" -> batch["PayloadKey"],
  "RequestIdentities" -> Lookup[entries, "RequestIdentity"],
  "CoefficientIdentities" -> Lookup[entries, "CoefficientIdentity"],
  "ObservableIdentities" -> Lookup[entries, "Identity"],
  "ObservableCheckpointIdentities" ->
    Lookup[entries, "CheckpointIdentity"],
  "DeltaPrescriptions" -> delta,
  "DeltaPrescriptionIdentity" -> ft2CanonicalIdentity[
    "ft2-delta-prescriptions-", delta],
  "ExtraSingularFactorsIdentity" -> ft2CanonicalIdentity[
    "ft2-extra-singular-factors-", {}],
  "AtlasPlanIdentity" -> atlasIdentity,
  "NativeBatchPayloadIdentity" -> payloadIdentity,
  "NativeEpsilonPlanIdentity" -> execution["Identity"],
  "SourceCompleteMax" -> ledger["SourceCompleteMax"],
  "TargetCompleteMax" -> ledger["TargetCompleteMax"],
  "DeliverableCompleteMax" -> ledger["DeliverableCompleteMax"],
  "RequiredRawTop" -> ledger["DownstreamRawTop"],
  "CoefficientHalo" -> ledger["CoefficientHalo"],
  "IntegrationHalo" -> ledger["IntegrationHalo"],
  "MatchEpsilonPadding" ->
    ledger["SourceCompleteMax"] - ledger["CoefficientHalo"] -
      ledger["TargetCompleteMax"]|>;
checkpointIdentity = "native-resume-state";
nativeStateFile = FileNameJoin[{tmpDir,
  name <> "_level1_native_transport.de2cp"}];
resumeCore = <|"Schema" -> "FeynmanTrick.NativeTransportResume/v1",
  "ContractIdentity" -> contract["Identity"],
  "AtlasPlanIdentity" -> atlasIdentity,
  "NativeBatchPayloadIdentity" -> payloadIdentity,
  "CheckpointIdentity" -> checkpointIdentity,
  "State" -> <|"CheckpointIdentity" -> checkpointIdentity,
    "Path" -> nativeStateFile|>|>;
resumeRecord = Append[resumeCore, "Identity" -> ft2CanonicalIdentity[
  "ft2-native-transport-resume-", resumeCore]];

Export[nativeStateFile, "opaque fixture", "Text"];
nativeSidecar = FileNameJoin[{tmpDir,
  name <> "_level1_native_transport.mx"}];
nativePayload = <|"Kind" -> "NativeTransport", "Example" -> name,
  "Level" -> 1, "PrepKey" -> prepKey, "System" -> system,
  "Variable" -> z, "BoundaryValues" -> boundaryValues,
  "BoundaryPrefactors" -> boundaryPrefactors,
  "EpsilonBasis" -> basis["CheckpointRecord"],
  "MastersHere" -> mastersHere, "MastersBelow" -> mastersBelow,
  "Requests" -> requests,
  "Reductions" -> Association[requests[[1, "NeededVec"]] -> 1],
  "ExtraSingularFactors" -> {}, "Anchor" -> anchor,
  "WorkingPrecision" -> wp, "DivisionOrder" -> divisionOrder,
  "RadiusOfConvergence" -> radiusOfConvergence,
  "ValueTransportMode" -> Environment["DE2_VALUE_TRANSPORT"],
  "RecurrenceBackend" -> "Cpp",
  "SingularMatchPrecondition" -> singularMatchPrecondition,
  "DeltaPrescriptionSign" -> deltaPrescriptionSign,
  "EpsilonOrder" -> epsOrder,
  "BoundaryExtraOrder" -> boundaryExtraOrder,
  "LevelEpsilonHalos" -> levelEpsilonHalos,
  "ExpansionOrder" -> expansionOrder,
  "RequestedEpsilonOrder" -> requiredRaw,
  "NativeLedger" -> ledger, "NativeObservableBatch" -> audit,
  "NativeTransportContract" -> contract,
  "NativeTransportCheckpoint" -> resumeRecord,
  "NativeEpsilonPlan" -> execution["Record"],
  "NativeEpsilonPlanIdentity" -> execution["Identity"],
  "Tainted" -> False|>;
writeNative = saveLadderCheckpoint[nativeSidecar, nativePayload];
loadedNative = loadLadderCheckpoint[
  nativeSidecar, name, data, prepKey, nativePlan];
assert["schema-2 native sidecar passes the normal production loader",
  writeNative === nativeSidecar && AssociationQ[loadedNative] &&
    loadedNative["Kind"] === "NativeTransport" &&
    loadedNative["NativeTransportContract"] === contract,
  loadedNative];

boundaryAudit = Join[audit, <|
  "SourceCompleteMax" -> requiredRaw,
  "TargetCompleteMax" -> requiredRaw,
  "DeliverableCompleteMax" -> requiredRaw,
  "RequiredRawTop" -> requiredRaw,
  "CoefficientHalo" -> 0, "IntegrationHalo" -> 0,
  "MatchEpsilonPadding" -> 0|>];
boundaryFile = FileNameJoin[{tmpDir, name <> "_level1_boundary.mx"}];
saveLadderCheckpoint[boundaryFile, <|"Kind" -> "Boundary",
  "Example" -> name, "Level" -> 1, "PrepKey" -> prepKey,
  "BoundaryValues" -> {Take[First[boundaryValues], requiredRaw + 1]},
  "BoundaryPrefactors" -> {0}, "MastersHere" -> mastersHere,
  "Anchor" -> anchor, "WorkingPrecision" -> wp,
  "RecurrenceBackend" -> "Cpp", "DeltaPrescriptionSign" ->
    deltaPrescriptionSign, "EpsilonOrder" -> epsOrder,
  "BoundaryExtraOrder" -> boundaryExtraOrder,
  "LevelEpsilonHalos" -> levelEpsilonHalos,
  "SourceExpansionOrder" -> expansionOrder,
  "RequestedEpsilonOrder" -> requiredRaw,
  "RequiredRawTop" -> requiredRaw,
  "PreservedRawCompleteMax" -> requiredRaw,
  "BoundaryShift" -> 0,
  "PreservedSourceCompleteMax" -> requiredRaw,
  "NativeObservableBatch" -> boundaryAudit,
  "NativeEpsilonPlan" -> execution["Record"],
  "NativeEpsilonPlanIdentity" -> execution["Identity"],
  "Tainted" -> False|>];
discovered = ft2DiscoveredLadderCheckpoint[
  name, data, prepKey, nativePlan];
assert["automatic startup prefers completed native state over same-level boundary data",
  AssociationQ[discovered] && discovered["Kind"] === "NativeTransport",
  discovered];

tamperedFile = FileNameJoin[{tmpDir,
  name <> "_level1_native_transport_tampered.mx"}];
saveLadderCheckpoint[tamperedFile, Join[nativePayload, <|
  "NativeTransportContract" -> Join[contract, <|
    "Identity" -> "tampered"|>]|>]];
assert["native loader rejects a tampered exact transport contract",
  loadLadderCheckpoint[tamperedFile, name, data, prepKey,
    nativePlan] === $Failed];

Clear[Global`$FT2LadderCheckpoint]; Get[nativeSidecar];
Global`$FT2LadderCheckpoint = Join[Global`$FT2LadderCheckpoint,
  <|"SourceFingerprint" -> "stale-native-source"|>];
DumpSave[nativeSidecar, Global`$FT2LadderCheckpoint];
Clear[Global`$FT2LadderCheckpoint];
assert["native state rejects stale provenance even when legacy stale opt-in is enabled",
  loadLadderCheckpoint[nativeSidecar, name, data, prepKey,
    nativePlan] === $Failed];

Quiet[DeleteDirectory[tmpDir, DeleteContents -> True]];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];
SetEnvironment["FT_LADDER_CHECKPOINT_DIR" -> None];
SetEnvironment["FT_RESUME_LADDER_CHECKPOINT" -> None];
SetEnvironment["FT_ALLOW_STALE_LADDER_CHECKPOINT" -> None];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
