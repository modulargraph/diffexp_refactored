(* Static/package-level contract for the DiffExp2 Feynman-trick facade.
   PipelinePlan must not launch wolframscript or FIRE. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

plan = FeynmanTrick`PipelinePlan["bubble",
  "WorkingPrecision" -> 250,
  "MatchingDigits" -> 18,
  "ExpansionOrder" -> 40,
  "CppThreads" -> 6,
  "FIREPath" -> FileNameJoin[{$TemporaryDirectory, "fire-facade-test"}],
  "LevelEpsilonHalos" -> {0, 4}];
minusPlan = FeynmanTrick`PipelinePlan["bubble",
  "DeltaPrescriptionSign" -> -1];
levelSignPlan = FeynmanTrick`PipelinePlan["bubble",
  "DeltaPrescriptionSign" -> -1,
  "LevelDeltaPrescriptionSigns" -> {1}];
invalidSignPlan = FeynmanTrick`PipelinePlan["bubble",
  "DeltaPrescriptionSign" -> 0];
invalidLevelSignPlan = FeynmanTrick`PipelinePlan["bubble",
  "LevelDeltaPrescriptionSigns" -> {1, 0}];
migrationPlan = FeynmanTrick`PipelinePlan["bubble",
  "MigrateLegacyPreparation" -> True];
nativeCheckpointPlan = FeynmanTrick`PipelinePlan["bubble",
  "SaveNativeTransportCheckpoint" -> True];
basisOnlyValuePlan = FeynmanTrick`PipelinePlan["bubble",
  "ValueTransport" -> True, "NativeValueHopExecution" -> False];
anchorVectorPlan = FeynmanTrick`PipelinePlan["banana_unequal",
  "FixedParameterValues" -> {1/5, 3/10, 2/5}];
invalidAnchorEndpointPlan = FeynmanTrick`PipelinePlan["bubble",
  "FixedParameterValues" -> {0}];
invalidInexactAnchorPlan = FeynmanTrick`PipelinePlan["bubble",
  "FixedParameterValues" -> {0.2}];

SetEnvironment["DE2_VALUE_TRANSPORT" -> None];
SetEnvironment["DE2_NATIVE_VALUE_HOP_EXECUTION" -> None];
defaultValueTransportRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["DE2_VALUE_TRANSPORT" -> "0"];
legacyValueTransportRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["DE2_VALUE_TRANSPORT" -> "1"];
SetEnvironment["DE2_NATIVE_VALUE_HOP_EXECUTION" -> "0"];
basisOnlyValuePlanRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["DE2_VALUE_TRANSPORT" -> "2"];
invalidValueTransportRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["DE2_VALUE_TRANSPORT" -> None];
SetEnvironment["DE2_NATIVE_VALUE_HOP_EXECUTION" -> "2"];
invalidNativeValueHopRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["DE2_NATIVE_VALUE_HOP_EXECUTION" -> None];

SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> "+1"];
plusRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> "-1"];
minusRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> "0"];
invalidRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> None];
SetEnvironment["FT_LEVEL_DELTA_PRESCRIPTION_SIGNS" -> "+1,-1,+1"];
levelSignRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_LEVEL_DELTA_PRESCRIPTION_SIGNS" -> "1,0"];
invalidLevelSignRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_LEVEL_DELTA_PRESCRIPTION_SIGNS" -> None];
SetEnvironment["FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT" -> "2"];
invalidNativeCheckpointRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT" -> None];
SetEnvironment["FT_FIXED_PARAMETER_VALUES" -> "1/5, 3/10"];
anchorVectorRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_FIXED_PARAMETER_VALUES" -> "1/5, 1"];
invalidAnchorVectorRunnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_FIXED_PARAMETER_VALUES" -> None];

rootFacadeDefinitions = DownValues[FeynmanTrick`PipelinePlan];
Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]];

assert["pipeline plan is an association", AssociationQ[plan]];
assert["Feynman-trick release version is exported",
  StringQ[FeynmanTrick`$FeynmanTrickVersion] &&
  StringQ[FeynmanTrick`$FeynmanTrickVersion::usage]];
assert["root Feynman-trick loader is idempotent",
  DownValues[FeynmanTrick`PipelinePlan] === rootFacadeDefinitions &&
  First[$ContextPath] === "FeynmanTrick`"];
assert["supported example registry is centralized",
  ContainsAll[FeynmanTrick`SupportedExamples[],
    {"bubble", "banana_unequal", "banana4_unequal", "kite"}]];
assert["unknown examples fail before launching a child",
  FailureQ[FeynmanTrick`PipelinePlan["definitely_not_an_example"]]];
assert["invalid FIRE paths fail in the parent plan",
  FailureQ[FeynmanTrick`PipelinePlan["bubble", "FIREPath" -> ""]] &&
  FailureQ[FeynmanTrick`PipelinePlan["bubble", "FIREPath" -> 17]]];
assert["malformed typed plans fail before process construction",
  FailureQ[FeynmanTrick`RunIntegrationPipeline[
    <|"Schema" -> "FeynmanTrick.PipelinePlan/v1"|>]]];
assert["root facade and native level-reduction seam load independently",
  Context[FeynmanTrick`RunIntegrationPipeline] === "FeynmanTrick`" &&
  Context[FeynmanTrick`LevelReduction`PrepareLevelIBPBatch] ===
    "FeynmanTrick`LevelReduction`" &&
  Length[DownValues[FeynmanTrick`RunIntegrationPipeline]] > 0 &&
  Length[DownValues[
    FeynmanTrick`LevelReduction`PrepareLevelIBPBatch]] > 0];
assert["pipeline plan schema",
  plan["Schema"] === "FeynmanTrick.PipelinePlan/v1"];
assert["C++ recurrence is the facade default",
  plan["Settings", "RecurrenceBackend"] === "Cpp" &&
  plan["Environment", "DE2_RECURRENCE_BACKEND"] === "Cpp"];
assert["delta prescription sign defaults to deterministic +1",
  plan["Settings", "DeltaPrescriptionSign"] === 1 &&
    plan["Environment", "FT_DELTA_PRESCRIPTION_SIGN"] === "1" &&
    AssociationQ[plusRunnerSettings] &&
    plusRunnerSettings["DeltaPrescriptionSign"] === 1];
assert["delta prescription sign carries strict -1 through plan and runner",
  AssociationQ[minusPlan] &&
    minusPlan["Settings", "DeltaPrescriptionSign"] === -1 &&
    minusPlan["Environment", "FT_DELTA_PRESCRIPTION_SIGN"] === "-1" &&
    AssociationQ[minusRunnerSettings] &&
    minusRunnerSettings["DeltaPrescriptionSign"] === -1 &&
    minusPlan["Settings"] =!= plan["Settings"] &&
    minusPlan["Environment"] =!= plan["Environment"]];
assert["delta prescription sign rejects zero in both public seams",
  FailureQ[invalidSignPlan] && FailureQ[invalidRunnerSettings]];
assert["per-level delta signs round-trip independently of the deepest rim",
  AssociationQ[levelSignPlan] &&
    levelSignPlan["Settings", "DeltaPrescriptionSign"] === -1 &&
    levelSignPlan["Settings", "LevelDeltaPrescriptionSigns"] === {1} &&
    levelSignPlan["Environment", "FT_DELTA_PRESCRIPTION_SIGN"] === "-1" &&
    levelSignPlan["Environment",
      "FT_LEVEL_DELTA_PRESCRIPTION_SIGNS"] === "1" &&
    AssociationQ[levelSignRunnerSettings] &&
    levelSignRunnerSettings["LevelDeltaPrescriptionSigns"] ===
      {1, -1, 1}];
assert["per-level delta signs reject non-unit entries",
  FailureQ[invalidLevelSignPlan] &&
    FailureQ[invalidLevelSignRunnerSettings]];
assert["legacy preparation migration is explicit in plan and environment",
  AssociationQ[migrationPlan] &&
    TrueQ[migrationPlan["Settings", "MigrateLegacyPreparation"]] &&
    migrationPlan["Environment", "FT_MIGRATE_LEGACY_PREP"] === "1" &&
    plan["Environment", "FT_MIGRATE_LEGACY_PREP"] === "0"];
assert["heavy native transport snapshots are explicit opt-in acceleration",
  plan["Settings", "SaveNativeTransportCheckpoint"] === False &&
    plan["Environment", "FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT"] === "0" &&
    AssociationQ[nativeCheckpointPlan] &&
    TrueQ[nativeCheckpointPlan["Settings",
      "SaveNativeTransportCheckpoint"]] &&
    nativeCheckpointPlan["Environment",
      "FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT"] === "1" &&
    FailureQ[invalidNativeCheckpointRunnerSettings]];
assert["fast transport settings are explicit",
  plan["Environment", "DE2_VALUE_TRANSPORT"] === "1" &&
    plan["Environment", "DE2_NATIVE_VALUE_HOP_EXECUTION"] === "1" &&
    plan["Environment", "FT_CPP_BATCH_ENDPOINT_ARMS"] === "1" &&
    plan["Environment", "DE2_CPP_THREADS"] === "6" &&
    plan["Settings", "MatchingDigits"] === 18 &&
    plan["Environment", "FT_MATCH_DIGITS"] === "18"];
assert["direct runners share the release value-transport default",
  AssociationQ[defaultValueTransportRunnerSettings] &&
    TrueQ[defaultValueTransportRunnerSettings["ValueTransport"]] &&
    AssociationQ[legacyValueTransportRunnerSettings] &&
    legacyValueTransportRunnerSettings["ValueTransport"] === False &&
    AssociationQ[basisOnlyValuePlanRunnerSettings] &&
    TrueQ[basisOnlyValuePlanRunnerSettings["ValueTransport"]] &&
    basisOnlyValuePlanRunnerSettings["NativeValueHopExecution"] === False &&
    FailureQ[invalidValueTransportRunnerSettings] &&
    FailureQ[invalidNativeValueHopRunnerSettings]];
assert["value-aware planning is independent of optional native value-hop execution",
  AssociationQ[basisOnlyValuePlan] &&
    TrueQ[basisOnlyValuePlan["Settings", "ValueTransport"]] &&
    basisOnlyValuePlan["Settings", "NativeValueHopExecution"] === False &&
    basisOnlyValuePlan["Environment", "DE2_VALUE_TRANSPORT"] === "1" &&
    basisOnlyValuePlan["Environment",
      "DE2_NATIVE_VALUE_HOP_EXECUTION"] === "0"];
assert["FIRE installation path is explicit and reproducible",
  plan["Environment", "FT_FIRE_PATH"] ===
    ExpandFileName[FileNameJoin[{$TemporaryDirectory, "fire-facade-test"}]] &&
  plan["Settings", "FIREPath"] === plan["Environment", "FT_FIRE_PATH"]];
assert["epsilon halos serialize deterministically",
  plan["Environment", "FT_LEVEL_EPS_HALOS"] === "0,4"];
assert["exact per-level anchors round-trip through plan and runner",
  AssociationQ[anchorVectorPlan] &&
    anchorVectorPlan["Settings", "FixedParameterValues"] ===
      {1/5, 3/10, 2/5} &&
    anchorVectorPlan["Environment", "FT_FIXED_PARAMETER_VALUES"] ===
      "1/5,3/10,2/5" &&
    AssociationQ[anchorVectorRunnerSettings] &&
    anchorVectorRunnerSettings["FixedParameterValues"] === {1/5, 3/10}];
assert["per-level anchors reject endpoints and inexact values",
  FailureQ[invalidAnchorEndpointPlan] &&
    FailureQ[invalidInexactAnchorPlan] &&
    FailureQ[invalidAnchorVectorRunnerSettings]];
assert["argv does not use a shell",
  ListQ[plan["Command"]] && plan["Command"][[1]] === "wolframscript" &&
    plan["Command"][[2]] === "-file"];

checkpoint = CreateTemporary[];
resume = FeynmanTrick`PipelinePlan["bubble",
  "ResumeFrom" -> checkpoint];
assert["resume checkpoint is wired into the typed plan",
  resume["ResumeFrom"] === ExpandFileName[checkpoint] &&
  resume["Environment", "FT_RESUME_LADDER_CHECKPOINT"] ===
    ExpandFileName[checkpoint] &&
  resume["CheckpointDirectory"] === DirectoryName[ExpandFileName[checkpoint]]];
DeleteFile[checkpoint];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
