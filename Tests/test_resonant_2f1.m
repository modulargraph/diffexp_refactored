(* Test: Resonant Hypergeometric 2F1 via DiffExp *)
(*
   Tests the general singular recurrence strategy for the hypergeometric ODE
   with RESONANT eigenvalues (c = integer).

   z(1-z)y'' + [c - (a+b+1)z]y' - ab*y = 0

   System: Y = {y, y'}, dY/dz = M(z) * Y
   Parameters: a = 1/4, b = 1/3, c = 2

   Residue matrix M_0 = lim_{z->0} z*A(z) has eigenvalues 0 and -(c-1) = -1.
   Difference: 0 - (-1) = 1 (positive integer!) -> RESONANT

   The second Frobenius solution involves log(z) terms.
   The general singular recurrence must handle this resonance.

   Boundary conditions from the regular solution:
   - At z=0: 2F1(a,b;c;0) = 1 (regular solution, no logs)
   - y'(0) = ab/c = 1/24

   Strategy: expand around z=0 (the singular point), transport to z=1/2.
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Resonant Hypergeometric 2F1 Test"];
Print["(c=2, eigenvalue difference = 1)"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Hypergeometric parameters *)
a = 1/4;
b = 1/3;
c = 2;

Print["Hypergeometric parameters:"];
Print["a = ", a, ", b = ", b, ", c = ", c];
Print["Eigenvalues of residue at z=0: {0, ", 1-c, "} -> difference = ", c-1, " (resonant!)"];
Print["c - a - b = ", c - a - b, " > 0 (convergent at z=1)\n"];

(* Configuration for resonant 2F1 *)
ConfigResonant2F1 = {
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Resonant_Matrices"}] <> "/",
  Verbosity -> 3,
  UseMobius -> True,
  UsePade -> True,
  UseRationalRecurrence -> True,
  WorkingPrecision -> 200,
  ExpansionOrder -> 60,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0
};

Print["Loading resonant 2F1 configuration..."];
LoadConfiguration[ConfigResonant2F1];
Print["System loaded: ", DiffExp`State`NumIntegrals, " integrals\n"];

(* Boundary conditions: provide the regular 2F1 solution as BC *)
(* The solution is expanded around z=0 via the line z = x *)
Print["=== Setting up boundary conditions at z=0 ==="];
Print["Using the regular 2F1 solution (no logs):"];
Print["  y(z) = 2F1(1/4, 1/3; 2; z)"];
Print["  y'(z) = d/dz 2F1(1/4, 1/3; 2; z)\n"];

(* BCs as closed-form expressions in z *)
bcExpressions = {
  Hypergeometric2F1[a, b, c, z],
  D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z
};

(* Prepare BCs with line z = x (expansion around z=0) *)
Print["Preparing boundary conditions (series expansion around z=0)..."];
PreparedBCs = PrepareBoundaryConditions[bcExpressions, <|z -> x|>];
Print["Boundary conditions prepared.\n"];

(* Transport to z = 1/2 *)
Print["=== Transporting from z=0 to z=1/2 ==="];
Print["(This expansion goes through the resonant singular point z=0)"];

Result = TransportTo[PreparedBCs, <|z -> 1/2|>];

Print["\nTransport complete!"];
Print["Result point: ", Result["KinematicPoint"]];

(* Extract results *)
y2F1Result = Result["SeriesValues"][[1, 1]];
yPrimeResult = Result["SeriesValues"][[2, 1]];

Print["\n=== Comparison with Mathematica ==="];

(* Expected values at z = 1/2 *)
y2F1Exact = N[Hypergeometric2F1[a, b, c, 1/2], 200];
yPrimeExact = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 1/2, 200];

Print["2F1(1/4, 1/3; 2; 1/2):"];
Print["  DiffExp  = ", NumberForm[y2F1Result, 20]];
Print["  Exact    = ", NumberForm[y2F1Exact, 20]];
diff1 = Abs[y2F1Result - y2F1Exact];
Print["  Accuracy = ", ScientificForm[diff1, 3]];

Print["\nDerivative y'(1/2):"];
Print["  DiffExp  = ", NumberForm[yPrimeResult, 20]];
Print["  Exact    = ", NumberForm[yPrimeExact, 20]];
diff2 = Abs[yPrimeResult - yPrimeExact];
Print["  Accuracy = ", ScientificForm[diff2, 3]];

(* Test pass/fail - require at least 30 digits of accuracy *)
tolerance = 10^-30;
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
