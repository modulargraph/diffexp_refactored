(* Focused regression for the singular banana chart whose residue
   eigenvectors degenerate at eps=0.  The native recurrence may use the
   spectral frame at n=0, but every positive Taylor layer must be solved in
   the epsilon-regular reduced physical frame. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["FAIL: compiled DiffExp2 backend is not available"];
  Exit[1]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
req[min_, max_, n_] := <|"EpsWindow" -> <|"Min" -> min,
  "CompleteMax" -> max|>, "TOrder" -> n|>;
tags[ls_] := {#["a"], #["b"], #["p"]} & /@ ls["Sectors"];

fix = Get[FileNameJoin[
  {repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
bananaSystem[] := catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
bananaChart[sys_] := Module[{plan},
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, 0}]];
  SelectFirst[Reverse[plan["Charts"]], TrueQ[#["Singular"]] &]];

solve[backend_] := Module[{sys, cs, result, seconds},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 5, "DivisionOrder" -> 4,
    "StepDivisionOrder" -> 4, "RecurrenceBackend" -> backend,
    "Variables" -> {}}]];
  sys = bananaSystem[];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, bananaChart[sys]]];
  seconds = First@AbsoluteTiming[
    result = catchDE2[DiffExp2`Solve`SolveChart[cs, req[0, 5, 10]]]];
  {result, seconds}];

{wolfram, wolframSeconds} = solve["Wolfram"];
{native, nativeSeconds} = solve["Cpp"];
structure = !FailureQ[wolfram] && !FailureQ[native] &&
  Length[wolfram["Basis", "Columns"]] ===
    Length[native["Basis", "Columns"]] &&
  And @@ MapThread[(#1["EpsWindow"] === #2["EpsWindow"] &&
      tags[#1] === tags[#2]) &,
    {wolfram["Basis", "Columns"], native["Basis", "Columns"]}];
difference = If[structure,
  Max[Abs[Flatten[MapThread[Flatten[N[
      #1["Sectors"][[All, "Coeffs"]] -
      #2["Sectors"][[All, "Coeffs"]], 70]] &,
    {wolfram["Basis", "Columns"], native["Basis", "Columns"]}]]]],
  Infinity];

Print["banana epsilon-regular N10 seconds: Wolfram=",
  N[wolframSeconds, 5], " Cpp=", N[nativeSeconds, 5],
  " maxdiff=", difference];
If[structure && difference < 10^-80,
  Print["PASS: epsilon-regular banana singular recurrence parity"];
  Exit[0],
  Print["FAIL: epsilon-regular banana singular recurrence parity"];
  If[FailureQ[wolfram], Print["Wolfram failure: ", wolfram]];
  If[FailureQ[native], Print["Cpp failure: ", native]];
  Exit[1]];
