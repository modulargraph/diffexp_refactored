(* Process-free call-count/parity tests for the invocation-local level IBP
   bundle.  ReduceIntegrals is replaced by an exact fake, so no FIRE process,
   database, or license seat is used. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"],
  {General::shdw, Symbol::shdw}];
FeynmanTrick`SetFTOption["Verbosity", 0];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

Module[{topology, ftData, changedFtData, transport, calls = {}, baseline,
        baselineCalls, batch, optimized, optimizedCalls, keyAgain, changedKey,
        keyAfterGlobalOptionChange, badBatch, tamperedBatch,
        oldReductionCacheOption, oldAutoDetectRestrictions, dumpSource},
  topology = FeynmanTrick`FIREInterface`DefineTopology[
    "level_batch_stub", {Global`l1}, {},
    {1 - Global`l1^2, 2 - Global`l1^2}, {}];
  topology["ProblemNumber"] = 31;
  topology["WorkDirectory"] = FileNameJoin[{$TemporaryDirectory,
    "unused_level_batch_stub"}];
  topology["StartFileReady"] = True;
  topology["SetupFingerprintRecord"] = <|
    "Schema" -> "FeynmanTrick.FIRESetup/v1",
    "StartFileSHA256" -> "level-batch-start",
    "Propagators" -> topology["Propagators"],
    "Restrictions" -> {{-1, -1}},
    "AutoDetectRestrictions" -> False|>;

  (* The lower master is direct for CombinedPositions {1,2}. *)
  ftData = <|"Levels" -> <|
    0 -> <|"Masters" -> {{0, 0}}|>,
    1 -> <|
      "Masters" -> {{0, 0}},
      "CombinedPositions" -> {1, 2},
      "FeynmanParameter" -> Global`x,
      "Topology" -> topology,
      "DiffMatrix" -> {{0}}
    |>
  |>|>;
  transport = <|
    "SegmentData" -> {},
    "EpsPrefactorsAbove" -> {0},
    "EpsilonOrder" -> 1,
    "BoundaryValuesAbove" -> {{1, 0}}
  |>;
  oldReductionCacheOption = Lookup[
    FeynmanTrick`Private`$FTConfig, "ReductionCache", True];
  FeynmanTrick`SetFTOption["ReductionCache", False];

  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] := (
      AppendTo[calls, integrals];
      AssociationMap[2 Global`G[1, {0, 0}] &, integrals]);

    baseline = {
      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
        ftData, 1],
      FeynmanTrick`DiffExpIntegration`Private`RequiredTransportEpsilonOrder[
        ftData, 1, 1, {0}],
      FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
        ftData, 0, transport, 1]
    };
    baselineCalls = calls;

    calls = {};
    batch = FeynmanTrick`DiffExpIntegration`Private`PrepareLevelIBPBatch[
      ftData, 1];
    optimized = {
      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
        ftData, 1, batch],
      FeynmanTrick`DiffExpIntegration`Private`RequiredTransportEpsilonOrder[
        ftData, 1, 1, {0}, batch],
      FeynmanTrick`DiffExpIntegration`ComputeLevelBoundary[
        ftData, 0, transport, 1, batch]
    };
    optimizedCalls = calls;
  ];

  assert["independent consumers previously request the same FIRE batch thrice",
    baselineCalls === {{{0, 0}}, {{0, 0}}, {{0, 0}}}];
  assert["explicit level bundle performs one FIRE reduction request",
    optimizedCalls === {{{0, 0}}}];
  assert["bundled singular-factor budget and boundary outputs are exact parity",
    optimized === baseline && AssociationQ[optimized[[3]]] &&
      optimized[[3, "BoundaryValues"]] === {{2, 0}}];

  keyAgain =
    FeynmanTrick`DiffExpIntegration`Private`levelIBPBatchSpec[
      ftData, 1]["Key"];
  changedFtData = ftData;
  changedFtData["Levels"][1]["Topology"]["SetupFingerprintRecord"]
    ["StartFileSHA256"] = "changed-level-batch-start";
  changedKey =
    FeynmanTrick`DiffExpIntegration`Private`levelIBPBatchSpec[
      changedFtData, 1]["Key"];
  assert["level bundle key is deterministic for identical exact inputs",
    batch["Key"] === keyAgain];
  assert["level bundle key changes with exact topology content",
    batch["Key"] =!= changedKey];
  oldAutoDetectRestrictions = Lookup[
    FeynmanTrick`Private`$FTConfig, "AutoDetectRestrictions", False];
  FeynmanTrick`SetFTOption[
    "AutoDetectRestrictions", !TrueQ[oldAutoDetectRestrictions]];
  keyAfterGlobalOptionChange =
    FeynmanTrick`DiffExpIntegration`Private`levelIBPBatchSpec[
      ftData, 1]["Key"];
  FeynmanTrick`SetFTOption[
    "AutoDetectRestrictions", oldAutoDetectRestrictions];
  assert["bundle key uses setup-time restrictions, not current global state",
    batch["Key"] === keyAfterGlobalOptionChange];

  badBatch = Join[batch, <|"Key" -> "stale"|>];
  assert["mismatched level bundle is rejected instead of reused",
    FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
      ftData, 1, badBatch] === $Failed];
  tamperedBatch = batch;
  tamperedBatch["CoefficientVectors"][{0, 0}] = {3};
  assert["same-shape coefficient payload mutation is rejected",
    FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
      ftData, 1, tamperedBatch] === $Failed];
  dumpSource = Import[
    FileNameJoin[{repoRoot, "Scripts", "dump_transport_checkpoints.m"}],
    "Text"];
  assert["cache-disabled checkpoint dumper threads one level batch",
    StringContainsQ[dumpSource,
      "Private`PrepareLevelIBPBatch[\n        ftData, level]"] &&
    StringContainsQ[dumpSource,
      "currentPrefactors, levelIBPBatch"] &&
    StringContainsQ[dumpSource,
      "ftData, level, levelIBPBatch"] &&
    StringContainsQ[dumpSource,
      "dtcFTEpsOrder, levelIBPBatch"]];
  FeynmanTrick`SetFTOption["ReductionCache", oldReductionCacheOption];
];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
