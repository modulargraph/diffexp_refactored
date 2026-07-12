(* Focused tests for the stepwise runner's Wolfram endpoint-arm snapshots
   and the retained-native rejection boundary.  The runner is loaded
   definition-only: no FIRE preparation or transport is performed here. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
(* Legacy partial-arm snapshots are an explicit Wolfram-backend feature.
   The retained C++ observable batch owns opaque native state and resumes
   only from completed numeric Boundary checkpoints. *)
SetEnvironment["DE2_RECURRENCE_BACKEND" -> "Wolfram"];
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
test["runner opens no Wolfram subkernels",
  And @@ (!StringContainsQ[runnerSource, #] & /@
    {"ParallelSubmit", "ParallelMap", "LaunchKernels", "ParallelNeeds"})];
test["retained native branch uses the observable batch dispatcher",
  StringContainsQ[runnerSource,
    "nativeDispatch = ft2RunNativeBoundaryDispatch["] &&
    StringContainsQ[runnerSource,
      "rawES = nativeDispatch[\"Values\"]"]];
test["legacy synchronous arm checkpoints remain in the Wolfram branch",
  StringContainsQ[runnerSource,
    "This write must finish before the expensive upper solve starts."] &&
    StringContainsQ[runnerSource,
      "saveTransportProgress[]];\n    If[needHi"]];
test["checkpoint replacement requests atomic overwrite",
  StringContainsQ[runnerSource,
    "RenameFile[tmp, file, OverwriteTarget -> True]"]];
test["prepared snapshot schema invalidates legacy reduction keys",
  $ftPrepCacheVersion === 2, $ftPrepCacheVersion];
test["epsilon-basis transport bumps the checkpoint schema",
  $ftLadderCheckpointVersion === 2, $ftLadderCheckpointVersion];

(* A two-link raw 1/eps chain is the minimal form of the double-box L3
   failure.  The exact diagonal basis must remove both poles, shift the
   incoming finite arrays without inventing coefficients, and remain
   idempotent when a transport checkpoint is resumed. *)
chainEps = Global`eps;
chainZ = Global`zBasis;
chainMatrix = {
  {0, 1/(chainEps*chainZ), 0},
  {0, 0, 1/(chainEps*chainZ)},
  {0, 0, 0}};
chainBoundary = Table[
  Table[Global`bc[i, k], {k, 0, 4}], {i, 1, 3}];
chainInputPrefactors = {2, 2, 2};
chainBasis = ft2NormalizeEpsilonBasis[
  chainMatrix, chainBoundary, chainInputPrefactors, chainEps];
test["raw 1/eps chain receives the exact relative normalization",
  AssociationQ[chainBasis] &&
    chainBasis["CheckpointRecord", "CanonicalPrefactors"] === {0, -1, -2} &&
    chainBasis["BoundaryPrefactors"] === {4, 3, 2} &&
    chainBasis["BoundaryShifts"] === {2, 1, 0}, chainBasis];
test["normalized chain matrix is exactly pole-free",
  AssociationQ[chainBasis] &&
    And @@ MapThread[TrueQ[PossibleZeroQ[Together[#1 - #2]]] &,
      {Flatten[chainBasis["Matrix"]],
       Flatten[{{0, 1/chainZ, 0}, {0, 0, 1/chainZ}, {0, 0, 0}}]}] &&
    !FeynmanTrick`EpsPrefactors`CheckEpsPoles[
      chainBasis["Matrix"], chainEps], chainBasis];
test["basis conversion retains, but never widens, the common boundary window",
  AssociationQ[chainBasis] &&
    chainBasis["InputCompleteMax"] === 4 &&
    chainBasis["CompleteMax"] === 4 &&
    Length /@ chainBasis["BoundaryValues"] === {5, 5, 5} &&
    chainBasis["BoundaryValues"] === {
      {0, 0, Global`bc[1, 0], Global`bc[1, 1], Global`bc[1, 2]},
      {0, Global`bc[2, 0], Global`bc[2, 1], Global`bc[2, 2], Global`bc[2, 3]},
      chainBoundary[[3]]}, chainBasis];
chainPhysicalBefore = MapThread[
  DiffExp2`EpsSeries`ESShift[DiffExp2`EpsSeries`ESNew[0, #1], -#2] &,
  {chainBoundary, chainInputPrefactors}];
chainPhysicalAfter = MapThread[
  DiffExp2`EpsSeries`ESShift[DiffExp2`EpsSeries`ESNew[0, #1], -#2] &,
  {chainBasis["BoundaryValues"], chainBasis["BoundaryPrefactors"]}];
test["physical boundary is unchanged on every certified handoff order",
  And @@ MapThread[DiffExp2`EpsSeries`ESSameQ,
    {chainPhysicalBefore, chainPhysicalAfter}],
  {DiffExp2`EpsSeries`ESWindow /@ chainPhysicalBefore,
   DiffExp2`EpsSeries`ESWindow /@ chainPhysicalAfter}];
test["relative basis shifts honestly reduce physical upper windows",
  DiffExp2`EpsSeries`ESCompleteMax /@ chainPhysicalBefore === {2, 2, 2} &&
    DiffExp2`EpsSeries`ESCompleteMax /@ chainPhysicalAfter === {0, 1, 2},
  {DiffExp2`EpsSeries`ESWindow /@ chainPhysicalBefore,
   DiffExp2`EpsSeries`ESWindow /@ chainPhysicalAfter}];
chainReductionCoefficients = {1 + chainEps, 2 - chainEps, 3 + chainEps^2};
chainReducedValue[rows_, prefactors_] := Module[{terms},
  terms = MapThread[
    DiffExp2`EpsSeries`ESTimes[
      DiffExp2`EpsSeries`ESFromExpression[
        Together[#1/chainEps^#3], chainEps, 4],
      DiffExp2`EpsSeries`ESNew[0, #2]] &,
    {chainReductionCoefficients, rows, prefactors}];
  Fold[DiffExp2`EpsSeries`ESAdd, First[terms], Rest[terms]]];
chainReducedBefore = chainReducedValue[chainBoundary, chainInputPrefactors];
chainReducedAfter = chainReducedValue[
  chainBasis["BoundaryValues"], chainBasis["BoundaryPrefactors"]];
test["IBP reduction is unchanged on its certified normalized-basis window",
  DiffExp2`EpsSeries`ESSameQ[chainReducedBefore, chainReducedAfter] &&
    DiffExp2`EpsSeries`ESCompleteMax[chainReducedAfter] <=
      DiffExp2`EpsSeries`ESCompleteMax[chainReducedBefore],
  {DiffExp2`EpsSeries`ESWindow[chainReducedBefore],
   DiffExp2`EpsSeries`ESWindow[chainReducedAfter]}];
chainBasisAgain = ft2NormalizeEpsilonBasis[chainMatrix,
  chainBasis["BoundaryValues"], chainBasis["BoundaryPrefactors"], chainEps];
test["checkpoint resume does not apply the diagonal shift twice",
  AssociationQ[chainBasisAgain] &&
    chainBasisAgain["BoundaryShifts"] === {0, 0, 0} &&
    chainBasisAgain["BoundaryValues"] === chainBasis["BoundaryValues"] &&
    chainBasisAgain["BoundaryPrefactors"] === chainBasis["BoundaryPrefactors"] &&
    chainBasisAgain["CheckpointRecord"] === chainBasis["CheckpointRecord"],
  chainBasisAgain];
test["nonuniform boundary windows fail instead of being padded",
  FailureQ[ft2NormalizeEpsilonBasis[{{0, 1/chainEps}, {0, 0}},
    {{1, 2}, {3}}, {0, 0}, chainEps]]];
noPoleBoundary = {{Global`u0, Global`u1}, {Global`v0, Global`v1}};
noPoleBasis = ft2NormalizeEpsilonBasis[
  {{chainZ, 0}, {0, -chainZ}}, noPoleBoundary, {3, 3}, chainEps];
test["pole-free levels with a common incoming prefactor remain exact no-ops",
  AssociationQ[noPoleBasis] &&
    noPoleBasis["Matrix"] === {{chainZ, 0}, {0, -chainZ}} &&
    noPoleBasis["BoundaryValues"] === noPoleBoundary &&
    noPoleBasis["BoundaryPrefactors"] === {3, 3} &&
    noPoleBasis["BoundaryShifts"] === {0, 0}, noPoleBasis];
test["IBP reductions are expressed in the active normalized basis",
  StringContainsQ[runnerSource,
    "eps^currentPrefactors[[j]]"] &&
    StringContainsQ[runnerSource,
      "currentPrefactors = epsilonBasis[\"BoundaryPrefactors\"]"]];
projectedRunnerPrescriptions = levelDeltaPrescriptions[chainZ,
  <|"SingularFactors" -> {}|>, {chainZ - 1/3 + chainEps}];
test["runner projects extra factors before constructing delta prescriptions",
  AnyTrue[projectedRunnerPrescriptions,
    FreeQ[First[#], chainEps] &&
      TrueQ[PossibleZeroQ[First[#] /. chainZ -> 1/3]] && Last[#] === 1 &] &&
    AllTrue[projectedRunnerPrescriptions, FreeQ[First[#], chainEps] &],
  projectedRunnerPrescriptions];

cacheTopology = FeynmanTrick`FIREInterface`DefineTopology[
  "snapshot-key", {Global`l1}, {},
  {1 - Global`l1^2, 2 - Global`l1^2}, {}];
cacheTopology["WorkDirectory"] = tmpDir;
cacheTopology["ProblemNumber"] = 41;
cacheTopology["StartFileReady"] = True;
cacheTopology["SetupFingerprintRecord"] = <|
  "Schema" -> "FeynmanTrick.FIRESetup/v1",
  "StartFileSHA256" -> "snapshot-start",
  "Restrictions" -> {{-1, -1}}|>;
oldAutoDetectRestrictions =
  FeynmanTrick`Private`$FTConfig["AutoDetectRestrictions"];
prepIdentityKey = ftPrepKey["snapshot-key", cacheTopology, {{1, 2}}];
FeynmanTrick`SetFTOption[
  "AutoDetectRestrictions", !TrueQ[oldAutoDetectRestrictions]];
prepChangedOptionKey = ftPrepKey[
  "snapshot-key", cacheTopology, {{1, 2}}];
FeynmanTrick`SetFTOption[
  "AutoDetectRestrictions", oldAutoDetectRestrictions];
test["prepared snapshot key covers setup-affecting FIRE options",
  prepIdentityKey =!= prepChangedOptionKey,
  {prepIdentityKey, prepChangedOptionKey}];
prepRuntimeA = Block[{
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord},
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[] :=
    <|"Runtime" -> "A"|>;
  ftPrepKey["snapshot-key", cacheTopology, {{1, 2}}]];
prepRuntimeB = Block[{
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord},
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[] :=
    <|"Runtime" -> "B"|>;
  ftPrepKey["snapshot-key", cacheTopology, {{1, 2}}]];
test["prepared snapshot key covers FIRE runtime identity",
  prepRuntimeA =!= prepRuntimeB, {prepRuntimeA, prepRuntimeB}];
cacheData = <|"NumLevels" -> 1, "Levels" -> <|
  0 -> <|"Masters" -> {{1, 1}}|>,
  1 -> <|"Masters" -> {{2, 0}}, "CombinedPositions" -> {1, 2},
    "Topology" -> cacheTopology, "Computed" -> True,
    "DiffMatrix" -> {{0}}|>|>|>;
cacheNeeded = {2, 0};
cacheExpectedKey =
  FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    cacheTopology, cacheNeeded];
test["prepared snapshots use FIREInterface's hardened reduction key",
  requiredReductionKeys[cacheData] === {cacheExpectedKey},
  requiredReductionKeys[cacheData]];
oldReductionCache = FeynmanTrick`FIREInterface`Private`$ReductionCache;
FeynmanTrick`FIREInterface`Private`$ReductionCache = Association[
  cacheExpectedKey -> <|"Reduction" -> 1, "Masters" -> {{2, 0}}|>];
prepSnapshotFile = FileNameJoin[{tmpDir, "hardened-prep.mx"}];
prepSnapshotWrite = savePreparedFT[
  prepSnapshotFile, 24680, cacheData];
FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
prepSnapshotLoad = loadPreparedFT[prepSnapshotFile, 24680];
test["hardened reduction keys survive prepared snapshot round-trip",
  prepSnapshotWrite === prepSnapshotFile &&
    AssociationQ[prepSnapshotLoad] &&
    preparedReductionCacheQ[prepSnapshotLoad,
      FeynmanTrick`FIREInterface`Private`$ReductionCache],
  {prepSnapshotWrite, prepSnapshotLoad}];
FeynmanTrick`FIREInterface`Private`$ReductionCache = oldReductionCache;

name = "checkpoint-fixture";
prepKey = 112358;
z = Global`z;
mastersHere = {{2, 0}};
mastersBelow = {{1, 1}};
levelData = <|"Masters" -> mastersHere, "FeynmanParameter" -> z,
  "CombinedPositions" -> {1, 2}|>;
data = <|"NumLevels" -> 1, "Levels" -> <|
  0 -> <|"Masters" -> mastersBelow|>, 1 -> levelData|>|>;
requests = FeynmanTrick`LevelReduction`BoundaryRequestRecords[
  mastersBelow, levelData["CombinedPositions"]];
reductions = Association[requests[[1, "NeededVec"]] -> 1];
low = <|"Charts" -> {<|"Arm" -> "lower"|>},
  "Final" -> <|"Endpoint" -> 0|>|>;
high = <|"Charts" -> {<|"Arm" -> "upper"|>},
  "Final" -> <|"Endpoint" -> 1|>|>;
file = FileNameJoin[{tmpDir, "fixture_level1_transport.mx"}];
fixtureBasis = ft2NormalizeEpsilonBasis[{{0}}, {{1}}, {0}, Global`eps];

basePayload = <|
  "Kind" -> "Transport", "Example" -> name, "Level" -> 1,
  "PrepKey" -> prepKey,
  "System" -> <|"Variable" -> z, "Matrix" -> fixtureBasis["Matrix"],
    "SingularFactors" -> {z}, "SingularFactorsExact" -> {z + Global`eps}|>,
  "Variable" -> z, "BoundaryValues" -> {{1}},
  "BoundaryPrefactors" -> {0},
  "EpsilonBasis" -> fixtureBasis["CheckpointRecord"],
  "MastersHere" -> mastersHere,
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
    loaded["TransportHigh"] === None &&
    loaded["System", "SingularFactors"] === {z} &&
    loaded["System", "SingularFactorsExact"] === {z + Global`eps}, loaded];

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
test["native backend rejects legacy partial-arm transport snapshots",
  Block[{recurrenceBackend = "Cpp"},
    loadLadderCheckpoint[file, name, data, prepKey] === $Failed]];

badArmsFile = FileNameJoin[{tmpDir, "bad-arms.mx"}];
saveLadderCheckpoint[badArmsFile, Join[basePayload, <|
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Upper"}|>]];
test["inconsistent completed-arm metadata is rejected",
  loadLadderCheckpoint[badArmsFile, name, data, prepKey] === $Failed];

badBasisFile = FileNameJoin[{tmpDir, "bad-epsilon-basis.mx"}];
saveLadderCheckpoint[badBasisFile, Join[basePayload, <|
  "EpsilonBasis" -> Join[fixtureBasis["CheckpointRecord"],
    <|"Prefactors" -> {1}|>],
  "ChartCache" -> low["Charts"], "TransportLow" -> low,
  "TransportHigh" -> None, "CompletedArms" -> {"Lower"}|>]];
test["inconsistent epsilon-basis checkpoint metadata is rejected",
  loadLadderCheckpoint[badBasisFile, name, data, prepKey] === $Failed];

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
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];

Print["\n", passes, " passed, ", failures, " failed."];
If[failures > 0, Quit[1], Quit[0]];
