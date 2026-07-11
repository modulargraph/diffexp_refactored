(* Compare the chart-basis portion of established sequential lower/upper
   endpoint transport with the one-kernel C++ homogeneous-basis prewarm.
   Planning and boundary matching are deliberately outside the timing: the
   optimization changes only PrepareChart/SolveHomogeneous scheduling.  This
   uses the committed banana L1 matrix fixture; no FIRE process or FT
   preparation is involved.

   Environment:
     ARM_BENCH_WP=100
     ARM_BENCH_EXPANSION_ORDER=20
     ARM_BENCH_EPSILON_ORDER=5
     ARM_BENCH_DIVISION_ORDER=4
     ARM_BENCH_SCALAR=0
     DE2_CPP_THREADS=4
*)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envInt[name_, default_Integer] := Module[{raw = Environment[name], value},
  value = If[StringQ[raw], Quiet[Check[ToExpression[raw], $Failed]], $Failed];
  If[IntegerQ[value], value, default]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

wp = envInt["ARM_BENCH_WP", 100];
expansionOrder = envInt["ARM_BENCH_EXPANSION_ORDER", 20];
epsilonOrder = envInt["ARM_BENCH_EPSILON_ORDER", 5];
divisionOrder = envInt["ARM_BENCH_DIVISION_ORDER", 4];
scalarFixture = Environment["ARM_BENCH_SCALAR"] === "1";
anchor = 11/23;

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["compiled C++ backend is unavailable"]; Exit[2]];
If[FailureQ[catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> wp, "ExpansionOrder" -> expansionOrder,
    "EpsilonOrder" -> epsilonOrder, "DivisionOrder" -> divisionOrder,
    "StepDivisionOrder" -> divisionOrder,
    "RecurrenceBackend" -> "Cpp", "Variables" -> {}}]]], Exit[2]];

fixture = If[scalarFixture,
  <|"Matrix" -> {{(1 - 2 Global`eps)/Global`x +
      (1 + Global`eps)/(Global`x - 1) + 1/(Global`x + 2)}},
    "Variable" -> Global`x|>,
  Get[FileNameJoin[
    {repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]]];
makeSystem[] := catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fixture["Matrix"], "Variable" -> fixture["Variable"]|>]];
makePlans[sys_] := catchDE2[{
  DiffExp2`Transport`SegmentLine[sys, {anchor, 0}],
  DiffExp2`Transport`SegmentLine[sys, {anchor, 1}]}];
armReq = <|"EpsWindow" -> <|"Min" -> 0,
    "CompleteMax" -> epsilonOrder|>, "TOrder" -> expansionOrder|>;

sysSequential = makeSystem[];
plansSequential = makePlans[sysSequential];
If[FailureQ[sysSequential] || FailureQ[plansSequential], Exit[2]];
{sequentialSeconds, sequentialResult} = AbsoluteTiming[catchDE2[Module[
  {systems},
  systems = DiffExp2`Solve`PrepareChart[sysSequential, #] & /@
    Join[plansSequential[[1, "Charts"]], plansSequential[[2, "Charts"]]];
  DiffExp2`Solve`SolveHomogeneous[#, armReq] & /@ systems]]];
If[FailureQ[sequentialResult],
  Print["sequential chart solve failed: ", sequentialResult]; Exit[2]];

sysBatched = makeSystem[];
plansBatched = makePlans[sysBatched];
If[FailureQ[sysBatched] || FailureQ[plansBatched], Exit[2]];
{batchedSeconds, batchedResult} = AbsoluteTiming[catchDE2[Module[
  {lo = plansBatched[[1, "Charts"]], hi = plansBatched[[2, "Charts"]],
   loSystems, hiSystems, rounds, systems},
  loSystems = ConstantArray[None, Length[lo]];
  hiSystems = ConstantArray[None, Length[hi]];
  rounds = Max[Length[lo], Length[hi]];
  Do[
    systems = {};
    If[i <= Length[lo],
      loSystems[[i]] = DiffExp2`Solve`PrepareChart[sysBatched, lo[[i]]];
      AppendTo[systems, loSystems[[i]]]];
    If[i <= Length[hi],
      hiSystems[[i]] = DiffExp2`Solve`PrepareChart[sysBatched, hi[[i]]];
      AppendTo[systems, hiSystems[[i]]]];
    If[Length[systems] === 2 && systems[[1]] =!= systems[[2]],
      DiffExp2`Solve`PrewarmHomogeneousBatch[systems, armReq]],
    {i, rounds}];
  DiffExp2`Solve`SolveHomogeneous[#, armReq] & /@
    Join[loSystems, hiSystems]]]];
If[FailureQ[batchedResult],
  Print["batched chart solve failed: ", batchedResult]; Exit[2]];

sectors[result_] := Flatten[#["Sectors"] & /@ result["Columns"]];
tags[result_] := Map[{#["a"], #["b"], #["p"]} &, sectors[result]];
coefficients[result_] := Flatten[
  Lookup[sectors[result], "Coeffs"]];
structureParity = And @@ MapThread[tags[#1] === tags[#2] &,
  {sequentialResult, batchedResult}];
maxDifference = If[structureParity,
  Max[Abs[N[Flatten[MapThread[
    coefficients[#1] - coefficients[#2] &,
    {sequentialResult, batchedResult}]], 50]]], Infinity];

Print["ARM_BENCH ", ExportString[<|
  "Fixture" -> If[scalarFixture, "scalar-two-endpoint", "banana-L1"],
  "WorkingPrecision" -> wp,
  "ExpansionOrder" -> expansionOrder,
  "EpsilonOrder" -> epsilonOrder,
  "DivisionOrder" -> divisionOrder,
  "Threads" -> Environment["DE2_CPP_THREADS"],
  "LowerCharts" -> Length[plansBatched[[1, "Charts"]]],
  "UpperCharts" -> Length[plansBatched[[2, "Charts"]]],
  "SequentialSeconds" -> sequentialSeconds,
  "BatchedSeconds" -> batchedSeconds,
  "Speedup" -> sequentialSeconds/batchedSeconds,
  "StructureParity" -> structureParity,
  "MaxCoefficientDifference" -> maxDifference|>,
  "RawJSON", "Compact" -> True]];
If[!TrueQ[structureParity] || !TrueQ[maxDifference < 10^-60], Exit[1]];
