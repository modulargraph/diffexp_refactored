(* DiffExpIntegrate must extend its integration rules when the input
   carries deeper Logx towers than IMaxLogOrder (the auto-extension relies
   on the exported SeriesOps`NormalizeLogPower; an unexported symbol made
   the detection silently return no integer powers, so x^-1 Logx^k sources
   integrated as if constant - the box-family endpoint-segment bug). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get["DiffExp.m"], {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
test[name_String, ok_] := If[TrueQ[ok],
  passed++; Print["  PASS: ", name],
  failed++; Print["  FAIL: ", name]
];

x = DiffExp`Symbols`x;
LL = DiffExp`Symbols`Logx;
wp = 60;
DiffExp`State`DiffExpConfiguration["WorkingPrecision"] = wp;
DiffExp`Integration`UpdateIntReps[1];

(* NormalizeLogPower is exported and rationalizes numeric powers *)
test["NormalizeLogPower exported, integer passthrough",
  DiffExp`SeriesOps`NormalizeLogPower[3] === 3];
test["NormalizeLogPower rationalizes numeric power",
  DiffExp`SeriesOps`NormalizeLogPower[N[2, 30]] === 2];

(* integrate -Logx^2/(2x): must give -Logx^3/6 and raise IMaxLogOrder *)
src = SeriesData[x, 0, {SetPrecision[-1/2, wp]*LL^2}, -1, 6, 1];
out = DiffExp`Integration`DiffExpIntegrate[src, x];
c0 = DiffExp`SeriesOps`LogxCoeffNS[
  DiffExp`SingularityDecomposition`Private`GetCoefficientAtPower[out, 0], 3];
test["x^-1 Logx^2 integrates to Logx^3/6",
  NumericQ[N[c0]] && Abs[N[c0 + 1/6, 20]] < 10^-25];
test["IMaxLogOrder extended to 2", DiffExp`State`IMaxLogOrder >= 2];

(* the endpoint ladder: Y_n = Int[-Y_{n-1}/x], Y_0 = 1: Y_n = (-Logx)^n/n! *)
coupling = SeriesData[x, 0, {SetPrecision[-1, wp]}, -1, 10, 1];
y = SeriesData[x, 0, {SetPrecision[1, wp]}, 0, 10, 1];
Do[
  y = DiffExp`Integration`DiffExpIntegrate[
    DiffExp`SeriesOps`SExpand[coupling*y], x],
  {n, 1, 5}
];
truth = (-Log[SetPrecision[1/10, wp]])^5/120;
got = Normal[y] /. LL -> Log[SetPrecision[1/10, wp]] /. x -> SetPrecision[1/10, wp];
test["five-fold 1/x ladder reproduces (-Logx)^5/5!",
  Abs[N[got - truth, 20]] < 10^-18];

Print["\nResults: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
