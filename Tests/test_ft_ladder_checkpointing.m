(* Focused tests for the stepwise runner's resumable endpoint-arm snapshots.
   The runner is loaded definition-only: no FIRE preparation or transport is
   performed here. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
tmpDir = CreateDirectory[FileNameJoin[{$TemporaryDirectory,
  "DiffExp2_ft_checkpoint_test_" <> ToString[$ProcessID]}]];
SetEnvironment["FT_LADDER_CHECKPOINT_DIR" -> tmpDir];
SetEnvironment["FT_RESUME_LADDER_CHECKPOINT" -> ""];
SetEnvironment["FT_ALLOW_STALE_LADDER_CHECKPOINT" -> "0"];

Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passes = 0; failures = 0;
pass[name_] := (Print["  [PASS] ", name]; passes++);
fail[name_, detail_] := (Print["  [FAIL] ", name, ": ", detail]; failures++);
test[name_, condition_, detail_: None] :=
  If[TrueQ[condition], pass[name], fail[name, detail]];

runnerSource = Import[
  FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}], "Text"];
lowerPos = First@First@StringPosition[runnerSource,
  "Print[\"FTLADDER TRANSPORT ARM level=\", level, \" endpoint=lower\"]"];
lowerSavePos = First@First@StringPosition[runnerSource,
  "(* This write must finish before the expensive upper solve starts. *)"];
upperPos = First@First@StringPosition[runnerSource,
  "Print[\"FTLADDER TRANSPORT ARM level=\", level, \" endpoint=upper\"]"];
test["lower arm is saved before upper transport starts",
  lowerPos < lowerSavePos < upperPos, {lowerPos, lowerSavePos, upperPos}];
test["resume computes only missing endpoint arms",
  StringContainsQ[runnerSource, "needLo && !AssociationQ[trLoCache]"] &&
    StringContainsQ[runnerSource, "needHi && !AssociationQ[trHiCache]"]];
test["one-kernel C++ arm batching requires both missing arms",
  StringContainsQ[runnerSource,
    "needLo && needHi && !AssociationQ[trLoCache] &&"] &&
    StringContainsQ[runnerSource,
      "!AssociationQ[trHiCache]"] &&
    StringContainsQ[runnerSource,
      "DiffExp2`Solve`PrewarmHomogeneousBatch[roundSystems, armReq]"]];
test["arm batching only fills genuinely idle native workers",
  StringContainsQ[runnerSource,
    "Length[A] < cppArmThreadBudget"] &&
    StringContainsQ[runnerSource,
      "Length[roundSystems] === 2"] &&
    StringContainsQ[runnerSource,
      "roundSystems[[1]] =!= roundSystems[[2]]"]];
test["arm batching preflights the complete bounded cache",
  StringContainsQ[runnerSource,
    "DiffExp2`Solve`HomogeneousCacheCapacity[]"] &&
    StringContainsQ[runnerSource,
      "Length[armUniqueCharts] > armCacheCapacity"] &&
    StringContainsQ[runnerSource,
      "FTLADDER CPP ARM BATCH SKIP"]];
test["arm batching opens no Wolfram subkernels",
  And @@ (!StringContainsQ[runnerSource, #] & /@
    {"ParallelSubmit", "ParallelMap", "LaunchKernels", "ParallelNeeds"})];
test["native arm prewarm does not replace synchronous arm checkpoints",
  StringContainsQ[runnerSource,
    "This prewarm is pure cache state: it never marks an arm complete."] &&
    StringContainsQ[runnerSource,
      "saveTransportProgress[]];\n    If[needHi"]];
test["checkpoint replacement requests atomic overwrite",
  StringContainsQ[runnerSource,
    "RenameFile[tmp, file, OverwriteTarget -> True]"]];

name = "checkpoint-fixture";
prepKey = 112358;
z = Global`z;
mastersHere = {{2, 0}};
mastersBelow = {{1, 1}};
levelData = <|"Masters" -> mastersHere, "FeynmanParameter" -> z,
  "CombinedPositions" -> {1, 2}|>;
data = <|"NumLevels" -> 1, "Levels" -> <|
  0 -> <|"Masters" -> mastersBelow|>, 1 -> levelData|>|>;
requests = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
  mastersBelow, levelData["CombinedPositions"]];
reductions = Association[requests[[1, "NeededVec"]] -> 1];
low = <|"Charts" -> {<|"Arm" -> "lower"|>},
  "Final" -> <|"Endpoint" -> 0|>|>;
high = <|"Charts" -> {<|"Arm" -> "upper"|>},
  "Final" -> <|"Endpoint" -> 1|>|>;
file = FileNameJoin[{tmpDir, "fixture_level1_transport.mx"}];

basePayload = <|
  "Kind" -> "Transport", "Example" -> name, "Level" -> 1,
  "PrepKey" -> prepKey,
  "System" -> <|"Variable" -> z, "Matrix" -> {{0}}|>,
  "Variable" -> z, "BoundaryValues" -> {{1}},
  "BoundaryPrefactors" -> {0}, "MastersHere" -> mastersHere,
  "MastersBelow" -> mastersBelow, "Requests" -> requests,
  "Reductions" -> reductions, "ExtraSingularFactors" -> {},
  "Anchor" -> anchor, "WorkingPrecision" -> wp,
  "RecurrenceBackend" -> recurrenceBackend,
  "DivisionOrder" -> divisionOrder,
  "RadiusOfConvergence" -> radiusOfConvergence,
  "ValueTransportMode" -> Environment["DE2_VALUE_TRANSPORT"],
  "SingularMatchPrecondition" -> singularMatchPrecondition,
  "EpsilonOrder" -> epsOrder, "BoundaryExtraOrder" -> boundaryExtraOrder,
  "LevelEpsilonHalos" -> levelEpsilonHalos,
  "ExpansionOrder" -> expansionOrder,
  "RequestedEpsilonOrder" -> requestedEpsilonOrder[1]|>;

writeResult = saveLadderCheckpoint[file, Join[basePayload, <|
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["lower-arm checkpoint writes", writeResult === file && FileExistsQ[file],
  writeResult];
loaded = loadLadderCheckpoint[file, name, data, prepKey];
test["lower-only checkpoint is resumable",
  AssociationQ[loaded] && AssociationQ[loaded["TransportLow"]] &&
    loaded["TransportHigh"] === None, loaded];

(* Updating the same destination exercises atomic OverwriteTarget replacement. *)
writeResult = saveLadderCheckpoint[file, Join[basePayload, <|
  "ChartCache" -> high["Charts"], "TransportLow" -> None,
  "TransportHigh" -> high, "CompletedArms" -> {"Upper"}|>]];
loaded = loadLadderCheckpoint[file, name, data, prepKey];
test["upper-only checkpoint atomically replaces and resumes",
  writeResult === file && AssociationQ[loaded] &&
    loaded["TransportLow"] === None &&
    AssociationQ[loaded["TransportHigh"]], loaded];
test["atomic writer leaves no temporary snapshots",
  FileNames["*.tmp-*.mx", tmpDir] === {}, FileNames["*", tmpDir]];

saveLadderCheckpoint[file, Join[basePayload, <|
  "ChartCache" -> Join[low["Charts"], high["Charts"]],
  "TransportLow" -> low, "TransportHigh" -> high,
  "CompletedArms" -> {"Lower", "Upper"}|>]];
loaded = loadLadderCheckpoint[file, name, data, prepKey];
test["completed two-arm checkpoint remains resumable",
  AssociationQ[loaded] && AssociationQ[loaded["TransportLow"]] &&
    AssociationQ[loaded["TransportHigh"]], loaded];

badArmsFile = FileNameJoin[{tmpDir, "bad-arms.mx"}];
saveLadderCheckpoint[badArmsFile, Join[basePayload, <|
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Upper"}|>]];
test["inconsistent completed-arm metadata is rejected",
  loadLadderCheckpoint[badArmsFile, name, data, prepKey] === $Failed];

badChartsFile = FileNameJoin[{tmpDir, "bad-charts.mx"}];
saveLadderCheckpoint[badChartsFile, Join[basePayload, <|
  "ChartCache" -> {}, "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["chart cache inconsistent with cached arms is rejected",
  loadLadderCheckpoint[badChartsFile, name, data, prepKey] === $Failed];

configFile = FileNameJoin[{tmpDir, "bad-config.mx"}];
saveLadderCheckpoint[configFile, Join[basePayload, <|
  "DivisionOrder" -> divisionOrder + 1,
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["transport configuration fingerprint remains enforced",
  loadLadderCheckpoint[configFile, name, data, prepKey] === $Failed];

backendFile = FileNameJoin[{tmpDir, "bad-recurrence-backend.mx"}];
saveLadderCheckpoint[backendFile, Join[basePayload, <|
  "RecurrenceBackend" -> If[recurrenceBackend === "Cpp", "Wolfram", "Cpp"],
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["recurrence backend mismatch is rejected",
  loadLadderCheckpoint[backendFile, name, data, prepKey] === $Failed];

preconditionFile = FileNameJoin[{tmpDir, "bad-precondition-mode.mx"}];
saveLadderCheckpoint[preconditionFile, Join[basePayload, <|
  "SingularMatchPrecondition" -> !singularMatchPrecondition,
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["singular precondition mode mismatch is rejected",
  loadLadderCheckpoint[preconditionFile, name, data, prepKey] === $Failed];

highOrderTransportFile = FileNameJoin[{tmpDir, "high-order-transport.mx"}];
saveLadderCheckpoint[highOrderTransportFile, Join[basePayload, <|
  "ExpansionOrder" -> expansionOrder + 10,
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["higher-order transport checkpoint remains order-specific",
  loadLadderCheckpoint[highOrderTransportFile, name, data, prepKey] === $Failed];

lowOrderTransportFile = FileNameJoin[{tmpDir, "low-order-transport.mx"}];
saveLadderCheckpoint[lowOrderTransportFile, Join[basePayload, <|
  "ExpansionOrder" -> expansionOrder - 10,
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["lower-order transport checkpoint cannot downgrade the run",
  loadLadderCheckpoint[lowOrderTransportFile, name, data, prepKey] === $Failed];

highOrderBoundaryFile = FileNameJoin[{tmpDir, "high-order-boundary.mx"}];
saveLadderCheckpoint[highOrderBoundaryFile, <|
  "Kind" -> "Boundary", "Example" -> name, "Level" -> 1,
  "PrepKey" -> prepKey, "BoundaryValues" -> {{1}},
  "BoundaryPrefactors" -> {0}, "MastersHere" -> mastersHere,
  "Anchor" -> anchor, "WorkingPrecision" -> wp,
  "RecurrenceBackend" -> recurrenceBackend,
  "EpsilonOrder" -> epsOrder, "BoundaryExtraOrder" -> boundaryExtraOrder,
  "LevelEpsilonHalos" -> levelEpsilonHalos,
  "SourceExpansionOrder" -> expansionOrder + 10,
  "RequestedEpsilonOrder" -> requestedEpsilonOrder[1]|>];
loaded = loadLadderCheckpoint[highOrderBoundaryFile, name, data, prepKey];
test["higher-order boundary checkpoint may seed a lower-order run",
  AssociationQ[loaded] &&
    loaded["SourceExpansionOrder"] === expansionOrder + 10, loaded];

lowOrderBoundaryFile = FileNameJoin[{tmpDir, "low-order-boundary.mx"}];
saveLadderCheckpoint[lowOrderBoundaryFile, <|
  "Kind" -> "Boundary", "Example" -> name, "Level" -> 1,
  "PrepKey" -> prepKey, "BoundaryValues" -> {{1}},
  "BoundaryPrefactors" -> {0}, "MastersHere" -> mastersHere,
  "Anchor" -> anchor, "WorkingPrecision" -> wp,
  "RecurrenceBackend" -> recurrenceBackend,
  "EpsilonOrder" -> epsOrder, "BoundaryExtraOrder" -> boundaryExtraOrder,
  "LevelEpsilonHalos" -> levelEpsilonHalos,
  "SourceExpansionOrder" -> expansionOrder - 10,
  "RequestedEpsilonOrder" -> requestedEpsilonOrder[1]|>];
test["lower-order boundary checkpoint cannot downgrade the run",
  loadLadderCheckpoint[lowOrderBoundaryFile, name, data, prepKey] === $Failed];

(* Rewrite a valid snapshot's source provenance to exercise stale opt-in. *)
staleFile = FileNameJoin[{tmpDir, "stale.mx"}];
saveLadderCheckpoint[staleFile, Join[basePayload, <|
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
Clear[Global`$FT2LadderCheckpoint]; Get[staleFile];
Global`$FT2LadderCheckpoint = Join[Global`$FT2LadderCheckpoint,
  <|"SourceFingerprint" -> "deliberately-stale"|>];
DumpSave[staleFile, Global`$FT2LadderCheckpoint];
Clear[Global`$FT2LadderCheckpoint];
allowStaleLadderCheckpoint = False;
test["stale source checkpoint is rejected by default",
  loadLadderCheckpoint[staleFile, name, data, prepKey] === $Failed];
allowStaleLadderCheckpoint = True;
loaded = loadLadderCheckpoint[staleFile, name, data, prepKey];
test["explicit stale opt-in still works for a partial arm",
  AssociationQ[loaded] && TrueQ[loaded["Tainted"]], loaded];

Quiet[DeleteDirectory[tmpDir, DeleteContents -> True]];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["FT_LADDER_CHECKPOINT_DIR" -> None];
SetEnvironment["FT_RESUME_LADDER_CHECKPOINT" -> None];
SetEnvironment["FT_ALLOW_STALE_LADDER_CHECKPOINT" -> None];

Print["\n", passes, " passed, ", failures, " failed."];
If[failures > 0, Quit[1], Quit[0]];
