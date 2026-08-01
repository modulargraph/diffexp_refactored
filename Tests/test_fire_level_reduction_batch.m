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

Module[{topology, ftData, changedFtData, calls = {}, baseline, baselineBatch,
        baselineCalls, batch, optimized, optimizedCalls, keyAgain, changedKey,
        keyAfterGlobalOptionChange, badBatch, tamperedBatch,
        poleBatch, poleBudget, regulatorBatch, fractionalBatch,
        fractionalBudget, singularBatch, singularFactors,
        singularFactorsCached, singularFactorCacheSize,
        oldReductionCacheOption, oldAutoDetectRestrictions,
        levelReductionSource, runnerSource},
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
  oldReductionCacheOption = Lookup[
    FeynmanTrick`Private`$FTConfig, "ReductionCache", True];
  FeynmanTrick`SetFTOption["ReductionCache", False];

  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] := (
      AppendTo[calls, integrals];
      AssociationMap[2 Global`G[1, {0, 0}] &, integrals]);

    baselineBatch = FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[
      ftData, 1];
    baseline = {
      FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
        ftData, 1],
      FeynmanTrick`LevelReduction`RequiredTransportEpsilonOrder[
        ftData, 1, 1, {0}],
      baselineBatch["Reductions"],
      baselineBatch["CoefficientVectors"]
    };
    baselineCalls = calls;

    calls = {};
    batch = FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[
      ftData, 1];
    optimized = {
      FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
        ftData, 1, batch],
      FeynmanTrick`LevelReduction`RequiredTransportEpsilonOrder[
        ftData, 1, 1, {0}, batch],
      batch["Reductions"],
      batch["CoefficientVectors"]
    };
    optimizedCalls = calls;
  ];

  assert["independent consumers previously request the same FIRE batch thrice",
    baselineCalls === {{{0, 0}}, {{0, 0}}, {{0, 0}}}];
  assert["explicit level bundle performs one FIRE reduction request",
    optimizedCalls === {{{0, 0}}}];
  assert["bundled factors, epsilon budget, and boundary reduction are exact parity",
    optimized === baseline && optimized[[1]] === {} &&
      optimized[[2]] === 1 &&
      optimized[[4]][{0, 0}] === {2}];

  (* FIRE already publishes one rational numerator times one denominator
     power.  Extract its factors directly, but still remove an exact
     numerator/denominator cancellation and memoize the deterministic batch
     result so a private epsilon retry cannot repeat the factorization. *)
  FeynmanTrick`LevelReduction`Private`$levelIBPSingularFactorCache = <||>;
  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] :=
      AssociationMap[
        ((Global`x^2 - 1)/
          ((Global`x - 1) (Global`x - 2))) Global`G[1, {0, 0}] &,
        integrals];
    singularBatch =
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, 1];
  ];
  singularFactors =
    FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
      ftData, 1, singularBatch];
  singularFactorCacheSize = Length[
    FeynmanTrick`LevelReduction`Private`$levelIBPSingularFactorCache];
  singularFactorsCached = Block[{
      FeynmanTrick`LevelReduction`Private`coefficientSingularFactors},
    FeynmanTrick`LevelReduction`Private`coefficientSingularFactors[___] :=
      Failure["UnexpectedSingularFactorRecomputation", <||>];
    FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
      ftData, 1, singularBatch]];
  assert["single-rational factor extraction removes exact cancellations",
    AssociationQ[singularBatch] &&
      Length[singularFactors] === 1 &&
      TrueQ[PossibleZeroQ[Expand[First[singularFactors] -
        (Global`x - 2)]]]];
  assert["singular factors are memoized by exact batch payload",
    singularFactorsCached === singularFactors &&
      singularFactorCacheSize === 1 &&
      Length[
        FeynmanTrick`LevelReduction`Private`$levelIBPSingularFactorCache]
        === 1];

  (* Epsilon budgeting is exact bookkeeping, not numerical cleanup.  A tiny
     but nonzero exact pole must deepen the requested frame, analytic
     regulators must remain symbolic, and unsupported fractional epsilon
     frames must fail rather than being floored. *)
  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] :=
      AssociationMap[
        Global`G[1, {0, 0}]/(10^100 (Global`d - 4)) &,
        integrals];
    poleBatch = FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, 1];
    poleBudget = FeynmanTrick`LevelReduction`RequiredTransportEpsilonOrder[
      ftData, 1, 1, {0}, poleBatch];
  ];
  assert["tiny exact epsilon pole is never chopped from the budget",
    AssociationQ[poleBatch] && poleBudget === 2];

  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] :=
      AssociationMap[
        (1 + Global`analyticRegulator) Global`G[1, {0, 0}]/
          (Global`d - 4) &,
        integrals];
    regulatorBatch =
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, 1];
  ];
  assert["analytic regulator remains exact in the shared coefficient vector",
    AssociationQ[regulatorBatch] &&
      regulatorBatch["CoefficientVectors"][{0, 0}] ===
        {(1 + Global`analyticRegulator)/(Global`d - 4)}];

  Block[{FeynmanTrick`FIREInterface`ReduceIntegrals},
    FeynmanTrick`FIREInterface`ReduceIntegrals[_, integrals_List] :=
      AssociationMap[
        Sqrt[FeynmanTrick`FTeps] Global`G[1, {0, 0}] &,
        integrals];
    fractionalBatch =
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, 1];
    fractionalBudget =
      FeynmanTrick`LevelReduction`RequiredTransportEpsilonOrder[
        ftData, 1, 1, {0}, fractionalBatch];
  ];
  assert["fractional epsilon frame is rejected instead of silently floored",
    AssociationQ[fractionalBatch] && fractionalBudget === $Failed];

  keyAgain =
    FeynmanTrick`LevelReduction`LevelIBPBatchSpec[
      ftData, 1]["Key"];
  changedFtData = ftData;
  changedFtData["Levels"][1]["Topology"]["SetupFingerprintRecord"]
    ["StartFileSHA256"] = "changed-level-batch-start";
  changedKey = FeynmanTrick`LevelReduction`LevelIBPBatchSpec[
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
    FeynmanTrick`LevelReduction`LevelIBPBatchSpec[ftData, 1]["Key"];
  FeynmanTrick`SetFTOption[
    "AutoDetectRestrictions", oldAutoDetectRestrictions];
  assert["bundle key uses setup-time restrictions, not current global state",
    batch["Key"] === keyAfterGlobalOptionChange];

  badBatch = Join[batch, <|"Key" -> "stale"|>];
  assert["mismatched level bundle is rejected instead of reused",
    FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
      ftData, 1, badBatch] === $Failed];
  tamperedBatch = batch;
  tamperedBatch["CoefficientVectors"][{0, 0}] = {3};
  assert["same-shape coefficient payload mutation is rejected",
    FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
      ftData, 1, tamperedBatch] === $Failed];
  levelReductionSource = Import[FileNameJoin[{
    repoRoot, "FeynmanTrick", "LevelReduction.m"}], "Text"];
  runnerSource = Import[FileNameJoin[{
    repoRoot, "Scripts", "run_ft_stepwise2.m"}], "Text"];
  assert["release level-reduction seam has no legacy package dependency",
    !StringContainsQ[levelReductionSource, "DiffExp`"] &&
      !StringContainsQ[runnerSource,
        "FeynmanTrick`DiffExpIntegration`"] &&
      StringContainsQ[runnerSource,
        "FeynmanTrick`LevelReduction`PrepareLevelIBPBatch"]];
  FeynmanTrick`SetFTOption["ReductionCache", oldReductionCacheOption];
];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
