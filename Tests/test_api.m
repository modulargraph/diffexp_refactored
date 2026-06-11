(* End-to-end API tests: the mini-FT closed-form pin. *)
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
  "ExpansionOrder" -> 40, "EpsilonOrder" -> 3}]];
x = Global`x; eps = Global`eps;
esC = DiffExp2`EpsSeries`ESCoefficient;
eqN[a_, b_, tol_] := TrueQ[Abs[N[a - b, 30]] < tol];

(* THE mini-FT pin: f' = (eps/x) f, f(11/23) = (11/23)^eps;
   Int_0^1 f dx = 1/(1+eps) = 1 - eps + eps^2 - eps^3 *)
sys = DiffExp2`API`LoadSystem[<|"Matrix" -> {{eps/x}}, "Variable" -> x|>];
assert["load_system", sys["SingularFactors"] === {x}];
bvals = Transpose[Table[{SeriesCoefficient[(11/23)^eps, {eps, 0, k}]}, {k, 0, 3}]];
li = catchDE2[DiffExp2`API`LineIntegral[sys, bvals, 11/23, {0, 1}, {1}]];
assert["line_integral_x_to_eps",
  !FailureQ[li] &&
  eqN[esC[li, 0], 1, 10^-20] &&
  eqN[esC[li, 1], -1, 10^-18] &&
  eqN[esC[li, 2], 1, 10^-16]];

(* singular-endpoint limit: 2-component f = {x^eps, const}:
   A = diag(eps/x, 0): limit at 0 of {1,1}.f = the constant *)
sys2 = DiffExp2`API`LoadSystem[<|"Matrix" -> {{eps/x, 0}, {0, 0}}, "Variable" -> x|>];
bv2 = Transpose[Table[{SeriesCoefficient[(11/23)^eps, {eps, 0, k}],
  If[k === 0, 7, 0]}, {k, 0, 3}]];
tr2 = catchDE2[DiffExp2`API`TransportEndpoint[sys2, bv2, 11/23, 0]];
lim2 = catchDE2[DiffExp2`API`EndpointLimitValues[tr2, {1, 1}]];
assert["endpoint_limit_drop_rule",
  !FailureQ[lim2] && eqN[esC[lim2, 0], 7, 10^-20]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
