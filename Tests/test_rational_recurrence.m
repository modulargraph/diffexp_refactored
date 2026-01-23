(* Test script for the rational recurrence method *)
(* Compares UseRationalRecurrence -> True against the default method *)
(* using the Hypergeometric 2F1 test case *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Rational Recurrence Method Test"];
Print["===========================================\n"];

Get["DiffExp.m"];

(* Hypergeometric parameters *)
a = 1/4;
b = 1/3;
c = 3/2;

(* Boundary conditions at z0 *)
z0 = 1/100;
y2F1AtStart = N[Hypergeometric2F1[a, b, c, z0], 200];
yPrimeAtStart = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z0, 200];

BCs = {
  <|z -> z0|>,
  {
    {y2F1AtStart},
    {yPrimeAtStart}
  }
};

(* Expected values at z = 1/2 *)
y2F1Exact = N[Hypergeometric2F1[a, b, c, 1/2], 200];
yPrimeExact = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 1/2, 200];

(* === Test 1: Default strategy === *)
Print["Test 1: Running with default strategy..."];
LoadConfiguration[{
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 200,
  ExpansionOrder -> 60,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0,
  IntegrationStrategy -> "Default",
  "UseRationalRecurrence" -> False
}];

result1 = TransportTo[BCs, <|z -> 1/2|>];
val1Y = result1["SeriesValues"][[1, 1]];
val1Yp = result1["SeriesValues"][[2, 1]];

diff1Y = Abs[val1Y - y2F1Exact];
diff1Yp = Abs[val1Yp - yPrimeExact];
Print["  2F1 accuracy (default): ", Floor[-Log10[diff1Y]], " digits"];
Print["  y'  accuracy (default): ", Floor[-Log10[diff1Yp]], " digits"];

(* === Test 2: Rational recurrence === *)
Print["\nTest 2: Running with rational recurrence..."];
LoadConfiguration[{
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 200,
  ExpansionOrder -> 60,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0,
  IntegrationStrategy -> "Default",
  "UseRationalRecurrence" -> True
}];

result2 = TransportTo[BCs, <|z -> 1/2|>];
val2Y = result2["SeriesValues"][[1, 1]];
val2Yp = result2["SeriesValues"][[2, 1]];

diff2Y = Abs[val2Y - y2F1Exact];
diff2Yp = Abs[val2Yp - yPrimeExact];
Print["  2F1 accuracy (recurrence): ", Floor[-Log10[diff2Y]], " digits"];
Print["  y'  accuracy (recurrence): ", Floor[-Log10[diff2Yp]], " digits"];

(* === Test 3: Compare methods === *)
Print["\nTest 3: Agreement between methods..."];
diffMethodsY = Abs[val1Y - val2Y];
diffMethodsYp = Abs[val1Yp - val2Yp];
Print["  2F1 diff between methods: ", ScientificForm[diffMethodsY, 3]];
Print["  y'  diff between methods: ", ScientificForm[diffMethodsYp, 3]];

(* === Summary === *)
Print["\n==========================================="];
tolerance = 10^-40;
test1Pass = diff2Y < tolerance && diff2Yp < tolerance;
test2Pass = diffMethodsY < tolerance && diffMethodsYp < tolerance;

numPass = Count[{test1Pass, test2Pass}, True];
Print["Results:"];
If[test1Pass,
  Print["  [PASS] Rational recurrence gives correct 2F1 results (>40 digits)"],
  Print["  [FAIL] Rational recurrence accuracy insufficient: ", ScientificForm[Max[diff2Y, diff2Yp], 3]]
];
If[test2Pass,
  Print["  [PASS] Methods agree to within tolerance"],
  Print["  [FAIL] Methods disagree: ", ScientificForm[Max[diffMethodsY, diffMethodsYp], 3]]
];
Print[numPass, " / 2 tests passed"];
If[numPass === 2, Print["All tests PASSED!"], Print["Some tests FAILED."]];
Print["==========================================="];
