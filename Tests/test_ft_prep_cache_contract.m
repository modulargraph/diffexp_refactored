(* Definitions-only tests for the stable FIRE/Feynman-trick preparation
   cache contract.  This test reads source files and snapshot payloads but
   never calls SetupFIRE, FIRE7, or a transport solver. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];

tmpDir = CreateDirectory[FileNameJoin[{$TemporaryDirectory,
  "DiffExp2_ft_prep_contract_test_" <> ToString[$ProcessID]}]];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
SetEnvironment["FT_PREP_CACHE_DIR" -> tmpDir];
SetEnvironment["FT_MIGRATE_LEGACY_PREP" -> "0"];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passes = 0; failures = 0;
pass[name_] := (Print["  [PASS] ", name]; passes++);
fail[name_, detail_] := (Print["  [FAIL] ", name, ": ", detail]; failures++);
test[name_, condition_, detail_:None] :=
  If[TrueQ[condition], pass[name], fail[name, detail]];

test["preparation-only runner mode is strict and explicit",
  TrueQ[ft2ParsePreparationOnly["0"] === False] &&
    TrueQ[ft2ParsePreparationOnly["1"] === True] &&
    FailureQ[ft2ParsePreparationOnly["yes"]]];

inputTopology = FeynmanTrick`FIREInterface`DefineTopology[
  "prep-contract", {Global`l1}, {},
  {1 - Global`l1^2, 2 - Global`l1^2}, {}];
sequence = {{1, 2}};
baseContract = ftPrepContractRecord[
  "prep-contract", inputTopology, sequence];
baseKey = ftPrepContractKey[baseContract];

test["prepared snapshot schema is the exact v3 contract",
  $ftPrepCacheVersion === 3 &&
    baseContract["Schema"] === "FeynmanTrick.PreparedFTContract/v3",
  {$ftPrepCacheVersion, baseContract}];
test["preparation manifest is the narrow reviewed three-file set",
  $ftPrepPreparationSourcePaths === {
    "FeynmanTrick/PropagatorAlgebra.m",
    "FeynmanTrick/FIREInterface.m",
    "FeynmanTrick/FeynmanTrickIteration.m"},
  $ftPrepPreparationSourcePaths];
test["preparation source identities contain no absolute paths",
  AllTrue[$ftPrepPreparationSourceIdentities, Function[identity,
    StringQ[identity["RelativePath"]] &&
      !StringStartsQ[identity["RelativePath"], "/"] &&
      !StringContainsQ[identity["RelativePath"], repoRoot] &&
      StringLength[identity["SHA256"]] === 64]],
  $ftPrepPreparationSourceIdentities];

excludedPaths = {
  "FeynmanTrick/DiffExp2Pipeline.m",
  "Scripts/run_ft_stepwise2.m",
  "FeynmanTrick/BoundaryConditions.m",
  "FeynmanTrick/LevelReduction.m"
};
allVirtualSources = Association@Map[
  # -> <|"RelativePath" -> #, "SHA256" -> IntegerString[
    Hash[{"fixture", #}, "SHA256"], 16, 64]|> &,
  Join[$ftPrepPreparationSourcePaths, excludedPaths]];
virtualKey[sourceMap_Association] := Block[{
    $ftPrepPreparationSourceIdentities =
      ftPrepSelectPreparationSourceIdentities[sourceMap]},
  ftPrepKey["prep-contract", inputTopology, sequence]];
virtualBaseKey = virtualKey[allVirtualSources];
excludedKeys = Table[
  virtualKey[Join[allVirtualSources, Association[
    path -> <|"RelativePath" -> path,
      "SHA256" -> IntegerString[Hash[{"edited", path}, "SHA256"],
        16, 64]|>]]],
  {path, excludedPaths}];
test["facade runner boundary and LevelReduction edits preserve prep identity",
  AllTrue[excludedKeys, # === virtualBaseKey &],
  Thread[excludedPaths -> excludedKeys]];

includedKeys = Table[
  virtualKey[Join[allVirtualSources, Association[
    path -> <|"RelativePath" -> path,
      "SHA256" -> IntegerString[Hash[{"edited", path}, "SHA256"],
        16, 64]|>]]],
  {path, $ftPrepPreparationSourcePaths}];
test["each preparation-module edit invalidates prep identity",
  AllTrue[includedKeys, # =!= virtualBaseKey &] &&
    DuplicateFreeQ[includedKeys],
  Thread[$ftPrepPreparationSourcePaths -> includedKeys]];

orchestrationOnlyKey = Block[{
    $ftLadderSourceFingerprint = Hash["orchestration-only-edit", "SHA256"],
    deltaPrescriptionSign = -1,
    wp = 777, expansionOrder = 123},
  ftPrepKey["prep-contract", inputTopology, sequence]];
test["ladder provenance and numerical transport settings do not enter prep key",
  orchestrationOnlyKey === baseKey,
  {baseKey, orchestrationOnlyKey}];

changedTopology = inputTopology;
changedTopology["Propagators"] = {1 - Global`l1^2, 3 - Global`l1^2};
test["exact topology and sequence changes invalidate prep identity",
  ftPrepKey["prep-contract", changedTopology, sequence] =!= baseKey &&
    ftPrepKey["prep-contract", inputTopology, {{2, 1}}] =!= baseKey];

oldFixedParameter = FeynmanTrick`Private`$FTConfig["FixedParameterValue"];
oldFixedParameterValues =
  FeynmanTrick`Private`$FTConfig["FixedParameterValues"];
FeynmanTrick`SetFTOption["FixedParameterValue", 7/13];
changedConfigKey = ftPrepKey["prep-contract", inputTopology, sequence];
FeynmanTrick`SetFTOption["FixedParameterValue", oldFixedParameter];
test["evaluated preparation configuration invalidates prep identity",
  changedConfigKey =!= baseKey, {baseKey, changedConfigKey}];

FeynmanTrick`SetFTOption["FixedParameterValues", {1/5}];
changedAnchorVectorKey = ftPrepKey["prep-contract", inputTopology, sequence];
FeynmanTrick`SetFTOption[
  "FixedParameterValues", oldFixedParameterValues];
test["per-level anchor vector invalidates prep identity",
  changedAnchorVectorKey =!= baseKey,
  {baseKey, changedAnchorVectorKey}];

runtimeA = Block[{
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord},
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[] :=
    <|"Runtime" -> "A"|>;
  ftPrepKey["prep-contract", inputTopology, sequence]];
runtimeB = Block[{
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord},
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[] :=
    <|"Runtime" -> "B"|>;
  ftPrepKey["prep-contract", inputTopology, sequence]];
test["exact FIRE runtime invalidates prep identity",
  runtimeA =!= runtimeB, {runtimeA, runtimeB}];

(* Synthetic v2 migration uses the same full structural, retained FIRE-setup,
   and reduction-key checks as a real snapshot.  The v1 fixture below exercises
   the separate explicit tuple-key re-map without inventing setup provenance. *)
runtime = baseContract["FIRERuntime"];
preparedTopology = inputTopology;
preparedTopology["WorkDirectory"] = tmpDir;
preparedTopology["ProblemNumber"] = 73;
preparedTopology["StartFileReady"] = True;
preparedTopology["SetupFingerprintRecord"] = Join[<|
  "Schema" -> "FeynmanTrick.FIRESetup/v1",
  "AutoDetectRestrictions" -> baseContract[
    "PreparationConfiguration", "AutoDetectRestrictions"],
  "StartFileSHA256" -> "fixture-start"|>, runtime];
preparedData = <|
  "TopTopology" -> inputTopology,
  "CombinationSequence" -> sequence,
  "NumericalPoint" -> {},
  "NumLevels" -> 1,
  "FixedParamValue" -> baseContract[
    "PreparationConfiguration", "FixedParameterValue"],
  "Levels" -> <|
    0 -> <|"Masters" -> {{1, 1}}|>,
    1 -> <|"Topology" -> preparedTopology, "Masters" -> {{2, 0}},
      "CombinedPositions" -> {1, 2}, "Computed" -> True,
      "DiffMatrix" -> {{0}}|>
  |>
|>;
neededKey = First[requiredReductionKeys[preparedData]];
legacyCache = Association[neededKey -> <|
  "Reduction" -> Global`G[1, {2, 0}], "Masters" -> {{2, 0}}|>];
FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
hydrationCalls = {};
hydratedBatches = Block[{
    FeynmanTrick`LevelReduction`PrepareLevelIBPBatch},
  FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[data_, level_] := Module[
    {key = First[requiredReductionKeys[data]]},
    AppendTo[hydrationCalls, level];
    AssociateTo[FeynmanTrick`FIREInterface`Private`$ReductionCache,
      key -> <|"Reduction" -> Global`G[1, {2, 0}],
        "Masters" -> {{2, 0}}|>];
    <|"UpperLevel" -> level|>];
  hydratePreparedReductionCache[preparedData]];
test["fresh preparation hydrates exact boundary reductions before snapshot",
  AssociationQ[hydratedBatches] && hydrationCalls === {1} &&
    preparedReductionCacheQ[preparedData,
      FeynmanTrick`FIREInterface`Private`$ReductionCache],
  {hydratedBatches, hydrationCalls,
    FeynmanTrick`FIREInterface`Private`$ReductionCache}];
legacyPayload = <|"Version" -> 2, "Key" -> 424242,
  "FTData" -> preparedData, "ReductionCache" -> legacyCache|>;
legacyFile = FileNameJoin[{tmpDir, "prep-contract_legacy.mx"}];
Global`$FT2PreparedSnapshot = legacyPayload;
DumpSave[legacyFile, Global`$FT2PreparedSnapshot];
Clear[Global`$FT2PreparedSnapshot];
targetFile = ftPrepFile["prep-contract", baseKey];
migrated = migrateLegacyPreparedFT[
  "prep-contract", targetFile, baseContract];
Clear[Global`$FT2PreparedSnapshot]; Get[targetFile];
migratedPayload = Global`$FT2PreparedSnapshot;
Clear[Global`$FT2PreparedSnapshot];
test["explicit v2 migration rewrites one fully validated snapshot as v3",
  AssociationQ[migrated] && FileExistsQ[targetFile] &&
    migratedPayload["Version"] === 3 &&
    migratedPayload["Contract"] === baseContract &&
    migratedPayload["Provenance", "Kind"] === "LegacyV2Validated" &&
    migratedPayload["Key"] === baseKey,
  migratedPayload];
badLegacy = Join[legacyPayload, <|"FTData" -> Join[preparedData,
  <|"CombinationSequence" -> {{2, 1}}|>]|>];
test["legacy migration rejects structurally stale prepared data",
  !legacyPreparedFTPayloadCompatibleQ[badLegacy, baseContract]];

v1Name = "prep-v1";
v1Contract = ftPrepContractRecord[v1Name, inputTopology, sequence];
v1Topology = KeyDrop[preparedTopology, {"SetupFingerprintRecord",
  "SetupFingerprint"}];
v1Data = Join[preparedData, <|"Levels" -> <|
  0 -> preparedData["Levels", 0],
  1 -> Join[preparedData["Levels", 1],
    <|"Topology" -> v1Topology|>]|>|>];
v1OldKey = Join[legacyV1TopologyIdentity[v1Topology], {{2, 0}}];
v1Payload = <|"Version" -> 1, "Key" -> 131313,
  "FTData" -> v1Data, "ReductionCache" -> Association[
    v1OldKey -> <|"Reduction" -> Global`G[1, {2, 0}],
      "Masters" -> {{2, 0}}|>]|>;
v1File = FileNameJoin[{tmpDir, v1Name <> "_legacy.mx"}];
Global`$FT2PreparedSnapshot = v1Payload;
DumpSave[v1File, Global`$FT2PreparedSnapshot];
Clear[Global`$FT2PreparedSnapshot];
v1TargetFile = ftPrepFile[v1Name, ftPrepContractKey[v1Contract]];
v1Migrated = migrateLegacyPreparedFT[
  v1Name, v1TargetFile, v1Contract];
Clear[Global`$FT2PreparedSnapshot]; Get[v1TargetFile];
v1MigratedPayload = Global`$FT2PreparedSnapshot;
Clear[Global`$FT2PreparedSnapshot];
test["explicit v1 migration exactly re-keys legacy tuple reductions",
  AssociationQ[v1Migrated] &&
    v1MigratedPayload["Provenance", "Kind"] === "LegacyV1Rekeyed" &&
    TrueQ[v1MigratedPayload["Provenance",
      "VerifiedFIRESetupProvenance"] === False] &&
    preparedReductionCacheQ[v1Data,
      v1MigratedPayload["ReductionCache"]],
  v1MigratedPayload];

FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
v1Reloaded = loadPreparedFT[v1TargetFile, v1Contract];
fireLaunchObserved = False;
v1ExactHit = Block[{
    FeynmanTrick`FIREInterface`Private`preparedTopologyCompatibleQ,
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6},
  FeynmanTrick`FIREInterface`Private`preparedTopologyCompatibleQ[___] :=
    (fireLaunchObserved = True; False);
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[___] :=
    (fireLaunchObserved = True; $Failed);
  FeynmanTrick`FIREInterface`Private`runFIRE6[___] :=
    (fireLaunchObserved = True; 99);
  FeynmanTrick`FIREInterface`ReduceIntegrals[v1Topology, {{2, 0}}]
];
test["migrated v1 exact cache hit cannot enter any FIRE launch seam",
  AssociationQ[v1Reloaded] && !TrueQ[fireLaunchObserved] &&
    AssociationQ[v1ExactHit] &&
    v1ExactHit[{2, 0}] === Global`G[1, {2, 0}],
  {fireLaunchObserved, v1ExactHit}];
v1CollisionCache = Association[
  Join[legacyV1TopologyIdentity[v1Topology], {{2}}] -> <|
    "Reduction" -> 1, "Masters" -> {{2, 0}}|>,
  Join[legacyV1TopologyIdentity[v1Topology], {{2, 0}}] -> <|
    "Reduction" -> 2, "Masters" -> {{2, 0}}|>];
test["v1 re-key rejects nonidentical normalized-key collisions",
  FailureQ[legacyV1RekeyReductionCache[v1Data, v1CollisionCache]]];
v1IncompleteCache = Association[
  Join[legacyV1TopologyIdentity[v1Topology], {{1, 0}}] -> <|
    "Reduction" -> 1, "Masters" -> {{2, 0}}|>];
test["v1 re-key rejects caches missing a current exact boundary key",
  FailureQ[legacyV1RekeyReductionCache[v1Data, v1IncompleteCache]]];

SetEnvironment["FT_MIGRATE_LEGACY_PREP" -> "1"];
migrationSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
SetEnvironment["FT_MIGRATE_LEGACY_PREP" -> "0"];
defaultMigrationSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
test["legacy migration is an explicit opt-in runner setting",
  TrueQ[migrationSettings["MigrateLegacyPreparation"]] &&
    TrueQ[defaultMigrationSettings["MigrateLegacyPreparation"] === False],
  {migrationSettings, defaultMigrationSettings}];

Quiet[DeleteDirectory[tmpDir, DeleteContents -> True]];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["FT_PREP_CACHE_DIR" -> None];
SetEnvironment["FT_MIGRATE_LEGACY_PREP" -> None];

Print["\n", passes, " passed, ", failures, " failed."];
If[failures > 0, Quit[1], Quit[0]];
