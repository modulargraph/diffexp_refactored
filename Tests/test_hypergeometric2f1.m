(* Test: Hypergeometric 2F1 via DiffExp *)
(*
   Tests DiffExp transport for the hypergeometric ODE:
   z(1-z)y'' + [c - (a+b+1)z]y' - ab*y = 0

   System: Y = {y, y'}, dY/dz = M(z) * Y
   Parameters: a = 1/4, b = 1/3, c = 3/2

   Boundary conditions from gamma function values:
   - At z=0: 2F1(a,b;c;0) = 1, y'(0) = ab/c
   - Start at z=0.01 (near regular singular point z=0)

   Target: z = 1/2 (within convergence radius |z| < 1)
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Hypergeometric 2F1 Test"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Hypergeometric parameters *)
a = 1/4;
b = 1/3;
c = 3/2;

Print["Hypergeometric parameters:"];
Print["a = ", a, ", b = ", b, ", c = ", c];
Print["c - a - b = ", c - a - b, " > 0 (convergent at z=1)\n"];

(* Show gamma function formula *)
Print["At z=1: 2F1(a,b;c;1) = Gamma[c]*Gamma[c-a-b] / (Gamma[c-a]*Gamma[c-b])"];
z1Value = Gamma[c]*Gamma[c-a-b] / (Gamma[c-a]*Gamma[c-b]);
Print["                     = ", N[z1Value, 15], "\n"];

(* Configuration for 2F1 *)
Config2F1 = {
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 200,
  ExpansionOrder -> 60,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0  (* Only eps^0 order for pure 2F1 *)
};

Print["Loading 2F1 configuration..."];
LoadConfiguration[Config2F1];
Print["System loaded: ", DiffExp`State`NumIntegrals, " integrals\n"];

(* Boundary conditions near z=0 *)
(* At z=0: 2F1(a,b;c;0) = 1, derivative = ab/c *)
z0 = 1/100;  (* Start near z=0 *)
Print["=== Boundary conditions near z=0 ==="];
Print["At z=0 exactly: 2F1 = 1, y' = ab/c = ", a*b/c];

(* Use Mathematica's high-precision values at z0 *)
y2F1AtStart = N[Hypergeometric2F1[a, b, c, z0], 200];
yPrimeAtStart = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z0, 200];

Print["At z = ", z0, ":"];
Print["  Y1 = ", NumberForm[y2F1AtStart, 15]];
Print["  Y2 = ", NumberForm[yPrimeAtStart, 15], "\n"];

(* Boundary conditions format with EpsilonOrder = 0 *)
BCs = {
  <|z -> z0|>,
  {
    {y2F1AtStart},
    {yPrimeAtStart}
  }
};

(* Target: transport to z = 1/2 *)
Print["=== Transporting to z = 1/2 ==="];

Result = TransportTo[BCs, <|z -> 1/2|>];

Print["\nTransport complete!"];
Print["Result point: ", Result["KinematicPoint"]];

(* Extract results *)
y2F1Result = Result["SeriesValues"][[1, 1]];
yPrimeResult = Result["SeriesValues"][[2, 1]];

Print["\n=== Comparison with Mathematica ==="];

(* Expected values at z = 1/2 *)
y2F1Exact = N[Hypergeometric2F1[a, b, c, 1/2], 200];
yPrimeExact = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 1/2, 200];

Print["2F1(1/4, 1/3; 3/2; 1/2):"];
Print["  DiffExp  = ", NumberForm[y2F1Result, 20]];
Print["  Exact    = ", NumberForm[y2F1Exact, 20]];
diff1 = Abs[y2F1Result - y2F1Exact];
Print["  Accuracy = ", ScientificForm[diff1, 3]];

Print["\nDerivative y'(1/2):"];
Print["  DiffExp  = ", NumberForm[yPrimeResult, 20]];
Print["  Exact    = ", NumberForm[yPrimeExact, 20]];
diff2 = Abs[yPrimeResult - yPrimeExact];
Print["  Accuracy = ", ScientificForm[diff2, 3]];

(* Test pass/fail - require at least 40 digits of accuracy *)
tolerance = 10^-40;
test1Pass = diff1 < tolerance;
test2Pass = diff2 < tolerance;

Print["\n==========================================="];
If[test1Pass && test2Pass,
  Print["TEST PASSED"];
  Print["  2F1 correct to ", Floor[-Log10[diff1]], " digits"];
  Print["  y'  correct to ", Floor[-Log10[diff2]], " digits"];,
  Print["TEST FAILED"];
  If[!test1Pass, Print["  - 2F1 mismatch: ", ScientificForm[diff1, 3]]];
  If[!test2Pass, Print["  - y' mismatch: ", ScientificForm[diff2, 3]]];
];
Print["==========================================="];
