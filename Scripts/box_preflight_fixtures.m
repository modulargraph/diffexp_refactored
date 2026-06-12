(* Box preflight, phase A (kernel + FIRE): build box ftData exactly the way
   Scripts/run_ft_stepwise2.m does, then vendor the per-level static
   artifacts to Tests/refs/bench/box_L<n>.m so every later audit/iteration
   runs WITHOUT FIRE (the dump-replay lesson).  Fixture contents per level:
   the seam-normalized DiffMatrix, CollectLevelIBPSingularFactors, masters,
   boundary requests, and the normalized reductions.  No transports, no
   integration: cheap static artifacts only. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

example = "box";
anchor = 11/23;

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];
(* private FIRE work dir: the default $TemporaryDirectory/FeynmanTrick is
   shared across agent kernels *)
Module[{wd = FileNameJoin[{$TemporaryDirectory,
    "FT_boxpre_" <> ToString[$ProcessID]}]},
  If[!DirectoryQ[wd], CreateDirectory[wd, CreateIntermediateDirectories -> True]];
  FeynmanTrick`SetFTOption["WorkDirectory", wd]];
FeynmanTrick`SetFTOption["DimensionExpression", FTExampleDimension[example]];

topology = FTExampleTopology[example, "step"];
If[topology === $Failed, Print["TOPOLOGY FAIL"]; Quit[1]];
ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  topology, FTExampleSequence[example], {}];
outputDir = FileNameJoin[{$TemporaryDirectory,
  "FT2pre_" <> example <> "_" <> ToString[$ProcessID]}];
If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
Print["RunFullIteration start t=", SessionTime[]];
ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
If[ftData === $Failed, Print["ITERATION FAIL"]; Quit[1]];
Print["RunFullIteration done t=", SessionTime[]];

nLevels = ftData["NumLevels"];
Print["NumLevels = ", nLevels];

(* the runner's seam normalization, verbatim (run_ft_stepwise2.m:113-117) *)
ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
dimExpr = FeynmanTrick`Private`DimensionExpression[];
normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
  ftEps -> Global`eps);

benchDir = FileNameJoin[{repoRoot, "Tests", "refs", "bench"}];
If[!DirectoryQ[benchDir], CreateDirectory[benchDir,
  CreateIntermediateDirectories -> True]];

Do[Module[
  {levelData = ftData["Levels"][level], levelBelow = ftData["Levels"][level - 1],
   var, A, mastersHere, mastersBelow, requests, neededVecs, reductions,
   extraFacs, fixture, path},
  var = levelData["FeynmanParameter"];
  A = normalizeFT[levelData["DiffMatrix"]];
  mastersHere = levelData["Masters"];
  mastersBelow = levelBelow["Masters"];
  requests = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
    mastersBelow, levelData["CombinedPositions"]];
  neededVecs = DeleteDuplicates[#["NeededVec"] & /@ requests];
  reductions = FeynmanTrick`FIREInterface`ReduceIntegrals[
    levelData["Topology"], neededVecs];
  If[reductions === $Failed, Print["FIRE FAIL level ", level]; Quit[1]];
  extraFacs = normalizeFT[
    FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
      ftData, level]];
  fixture = <|
    "Example" -> example,
    "Level" -> level,
    "Variable" -> var,
    "Matrix" -> A,
    "ExtraSingularFactors" -> extraFacs,
    "Masters" -> mastersHere,
    "MastersBelow" -> mastersBelow,
    "CombinedPositions" -> levelData["CombinedPositions"],
    "Requests" -> requests,
    "Reductions" -> Association[KeyValueMap[#1 -> normalizeFT[#2] &, reductions]],
    "Anchor" -> anchor,
    "FTSettings" -> <|"WorkingPrecision" -> 120, "ExpansionOrder" -> 40,
      "DivisionOrder" -> 4, "StepDivisionOrder" -> 4,
      "EpsilonOrder" -> Max[level + 4, 1]|>,
    "GeneratedBy" -> "Scripts/box_preflight_fixtures.m"|>;
  path = FileNameJoin[{benchDir, example <> "_L" <> ToString[level] <> ".m"}];
  Put[fixture, path];
  Print["LEVEL ", level, " var=", var, " d=", Length[A],
    " masters=", InputForm[mastersHere],
    " below=", InputForm[mastersBelow]];
  Print["  cases: ", InputForm[{#["MasterVec"], #["Case"], #["Vi"], #["Vj"],
    #["NeededVec"]} & /@ requests]];
  Print["  extraFacs: ", InputForm[extraFacs]];
  Print["  matrix denominator factors: ", InputForm[DeleteDuplicates[
    Select[Flatten[Map[FactorList[Denominator[Together[#]]][[All, 1]] &,
      Flatten[A]]], !FreeQ[#, var] &]]]];
  Print["  wrote ", path]],
  {level, nLevels, 1, -1}];

Print["FIXTURES DONE t=", SessionTime[]];
Quit[0];
