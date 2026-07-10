(* Production strict-recurrence lower-frame benchmark.

   Compares the adaptive singular homogeneous rectangle with the former
   scalar terminal rectangle on the actual banana L1 endpoint charts.

   Environment: ALF_ORDERS=10,20 (default), ALF_WP=100, ALF_EPS_ORDER=2. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envOr[name_, def_] := Module[{v = Environment[name]}, If[StringQ[v], v, def]];
orders = ToExpression /@ StringSplit[envOr["ALF_ORDERS", "10,20"], ","];
wp = ToExpression[envOr["ALF_WP", "100"]];
epsOrder = ToExpression[envOr["ALF_EPS_ORDER", "2"]];
If[!AllTrue[orders, IntegerQ[#] && # >= 10 &], Print["invalid ALF_ORDERS"]; Quit[2]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[e_] := Quiet[Catch[e, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> wp, "ExpansionOrder" -> Max[orders],
  "EpsilonOrder" -> epsOrder, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];

fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
If[FailureQ[sys], Print[sys]; Quit[1]];

req[n_] := <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> epsOrder|>,
  "TOrder" -> n|>;
emit[a_] := Print["ALF ", ExportString[a, "RawJSON", "Compact" -> True]];

failed = False;
Do[
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, endpoint}]];
  chart = SelectFirst[Reverse[plan["Charts"]], TrueQ[# ["Singular"]] &];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart]];
  Do[
    DiffExp2`Solve`ClearSolveCaches[];
    ta = First@AbsoluteTiming[
      adaptive = Block[{DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = False},
        catchDE2[DiffExp2`Solve`SolveChart[cs, req[n]]]]];
    DiffExp2`Solve`ClearSolveCaches[];
    tl = First@AbsoluteTiming[
      legacy = Block[{DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = True},
        catchDE2[DiffExp2`Solve`SolveChart[cs, req[n]]]]];
    parity = !FailureQ[adaptive] && !FailureQ[legacy] &&
      adaptive["Basis", "Columns"] === legacy["Basis", "Columns"] &&
      adaptive["Basis", "Specs"] === legacy["Basis", "Specs"];
    If[!TrueQ[parity], failed = True];
    rows = adaptive["Basis", "Diagnostics", "AdaptiveLowerFrames"];
    emit[<|"Endpoint" -> endpoint, "TaylorOrder" -> n,
      "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
      "ExactParity" -> parity,
      "AdaptiveSeconds" -> N[ta, 8], "TerminalSeconds" -> N[tl, 8],
      "Speedup" -> N[tl/ta, 8],
      "AdaptiveFrameWidths" -> rows[[All, "FrameWidth"]],
      "TerminalFrameWidths" -> rows[[All, "TerminalFrameWidth"]],
      "Attempts" -> rows[[All, "Attempts"]]|>],
    {n, orders}],
  {endpoint, {0, 1}}];

If[failed, Quit[1], Quit[0]];
