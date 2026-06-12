(* Scripts/gen_bench_fixtures.m — ONE-TIME generator for the per-chart
   benchmark fixtures consumed by Scripts/bench_chart.m.

   Runs the FT layer (FeynmanTrick + FIRE) exactly like
   Scripts/run_ft_stepwise2.m up to RunFullIteration, then freezes the
   LEVEL-1 system of each example to a plain .m file so the benchmark
   never needs the FT layer again:

     Tests/refs/bench/sunrise_L1.m   (d = 3, variable xx1)
     Tests/refs/bench/banana_L1.m    (d = 7, variable xx1)

   Each fixture is an Association:
     <|"Example", "Level", "Matrix" (exact eps-rational, FT symbols
       normalized: dimVar -> DimensionExpression, FTeps -> Global`eps),
       "Variable", "ExtraSingularFactors"
       (CollectLevelIBPSingularFactors, normalized), "Masters",
       "NumLevels"|>

   The generator aborts loudly if any non-System/Global symbol leaks
   into the fixture (the unexported-context silent no-op trap).

   Run (under the shared kernel lock):
     env WolframKernel='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel' \
       wolframscript -file Scripts/gen_bench_fixtures.m
   Expected: ~5-10 min (FIRE reductions dominate). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", 120];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];

outDir = FileNameJoin[{repoRoot, "Tests", "refs", "bench"}];
If[!DirectoryQ[outDir],
  CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

genOne[name_String] := Module[
  {topology, ftData, outputDir, levelData, var, A, extraFacs,
   ftEps, dimVar, dimExpr, normalizeFT, fixture, leaks, path},
  Print["FIXTURE ", name, " start t=", SessionTime[]];
  FeynmanTrick`SetFTOption["DimensionExpression", FTExampleDimension[name]];
  topology = FTExampleTopology[name, "step"];
  If[topology === $Failed, Print["TOPOLOGY FAIL ", name]; Quit[1]];
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, FTExampleSequence[name], {}];
  outputDir = FileNameJoin[{$TemporaryDirectory,
    "FT2bench_" <> name <> "_" <> ToString[$ProcessID]}];
  If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
  ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
  If[ftData === $Failed, Print["ITERATION FAIL ", name]; Quit[1]];
  Print["  iteration done t=", SessionTime[]];
  (* symbol normalization at the seam — the run_ft_stepwise2.m normalizeFT *)
  ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  dimExpr = FeynmanTrick`Private`DimensionExpression[];
  normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
    ftEps -> Global`eps);
  levelData = ftData["Levels"][1];
  var = levelData["FeynmanParameter"];
  A = normalizeFT[levelData["DiffMatrix"]];
  extraFacs = normalizeFT[
    FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[ftData, 1]];
  Print["  d=", Length[A], " var=", var,
    " extraFacs=", Length[extraFacs], " t=", SessionTime[]];
  fixture = <|"Example" -> name, "Level" -> 1,
    "Matrix" -> A, "Variable" -> var,
    "ExtraSingularFactors" -> extraFacs,
    "Masters" -> levelData["Masters"],
    "NumLevels" -> ftData["NumLevels"]|>;
  (* loud context-leak gate: a FeynmanTrick`... symbol in the fixture would
     re-parse as a fresh silent symbol in the FT-free benchmark kernel *)
  leaks = DeleteDuplicates[Cases[fixture,
    s_Symbol /; !MemberQ[{"System`", "Global`"}, Context[s]] :>
      Context[s] <> SymbolName[s], {0, Infinity}, Heads -> True]];
  If[leaks =!= {}, Print["CONTEXT LEAK in ", name, ": ", leaks]; Quit[1]];
  path = FileNameJoin[{outDir, name <> "_L1.m"}];
  Put[fixture, path];
  Print["  written ", path, " (", FileByteCount[path], " bytes) t=",
    SessionTime[]]];

genOne["sunrise"];
genOne["banana"];
Print["DONE t=", SessionTime[]];
Quit[0];
