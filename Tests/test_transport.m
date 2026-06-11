(* Smoke tests for DiffExp2/Transport.m (M4 parity suite grows at M5). *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];
passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100,
  "ExpansionOrder" -> 30, "EpsilonOrder" -> 2}]];
x = Global`x; eps = Global`eps;

(* TT1: regular transport f' = f from 0 to 1/2: exact e^(1/2) *)
sys = <|"Matrix" -> {{1}}, "Variable" -> x, "SingularFactors" -> {}|>;
plan = DiffExp2`Transport`SegmentLine[sys, {0, 1/2}];
res = catchDE2[DiffExp2`Transport`TransportLine[sys, {{1, 0, 0}}, plan]];
assert["tt1_regular_exponential",
  !FailureQ[res] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res["Value"], 0][[1]] - Exp[1/2], 30]] < 10^-25];

(* TT2: singular endpoint x^eps: transport (1/2)^eps from 1/2 to 0 *)
sys2 = <|"Matrix" -> {{eps/x}}, "Variable" -> x, "SingularFactors" -> {x}|>;
plan2 = DiffExp2`Transport`SegmentLine[sys2, {1/2, 0}];
bvals = Transpose[Table[{SeriesCoefficient[(1/2)^eps, {eps, 0, k}]}, {k, 0, 2}]];
res2 = catchDE2[DiffExp2`Transport`TransportLine[sys2, bvals, plan2]];
assert["tt2_singular_endpoint_weight",
  !FailureQ[res2] && res2["EndpointIsSingular"] &&
  Module[{f = res2["Final"], sec},
    sec = SelectFirst[f["Sectors"], PossibleZeroQ[#["b"] - 1] &];
    sec =!= Missing["NotFound"] &&
    Abs[N[sec["Coeffs"][[-f["EpsWindow", "Min"] + 1, 1, 1]] - 1, 20]] < 10^-20]];

(* TT3: regular line with a far singularity: f' = f/(x-2): f = c(x-2);
   f(0) = 1/2 -> c = -1/4 -> f(1) = 1/4 *)
sys3 = <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x, "SingularFactors" -> {x - 2}|>;
plan3 = DiffExp2`Transport`SegmentLine[sys3, {0, 1}];
res3 = catchDE2[DiffExp2`Transport`TransportLine[sys3, {{1/2, 0, 0}}, plan3]];
assert["tt3_far_singularity",
  !FailureQ[res3] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res3["Value"], 0][[1]] - 1/4, 30]] < 10^-25];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
