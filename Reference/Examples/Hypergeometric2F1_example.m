(* Hypergeometric 2F1 Example *)
(*
   The Gauss hypergeometric function 2F1(a,b;c;z) satisfies:
   z(1-z)y'' + [c - (a+b+1)z]y' - ab*y = 0

   Converting to first-order system with Y = {y, y'}:
   dY/dz = M(z) * Y
   M(z) = {{0, 1}, {ab/(z(1-z)), ((a+b+1)z-c)/(z(1-z))}}

   For DiffExp (system: eps * dY/dz = M * Y), we put M at eps^0.
   This means the constraint 0 = M_0 * Y_0 must be satisfied,
   which requires the matrix structure to allow non-trivial Y_0.
*)

(* Parameters for the hypergeometric function *)
a = 1/4;
b = 1/3;
c = 3/2;

Print["Hypergeometric parameters:"];
Print["a = ", a];
Print["b = ", b];
Print["c = ", c];
Print["c - a - b = ", c - a - b, " (must be > 0 for convergence at z=1)\n"];

(* The differential equation matrix *)
M = {
  {0, 1},
  {a*b/(z*(1-z)), ((a+b+1)*z - c)/(z*(1-z))}
};

Print["Differential equation matrix M(z):"];
Print[M // MatrixForm // Simplify];

(* Output directory - clean up and create fresh files *)
outputDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "Tests", "Hypergeometric2F1_Matrices"}];
Print["\nOutput directory: ", outputDir];

(* Delete any existing files *)
Quiet[DeleteFile[FileNameJoin[{outputDir, #}]] & /@
  {"dz_d.m", "dz_0.m", "dz_1.m", "dz_2.m", "dz_3.m", "dz_4.m"}];

(* Create matrices - put full M at eps^0, zeros at higher orders *)
M0 = M // Simplify;
MZero = {{0, 0}, {0, 0}};

Export[FileNameJoin[{outputDir, "dz_0.m"}], M0];
Export[FileNameJoin[{outputDir, "dz_1.m"}], MZero];
Export[FileNameJoin[{outputDir, "dz_2.m"}], MZero];
Export[FileNameJoin[{outputDir, "dz_3.m"}], MZero];
Export[FileNameJoin[{outputDir, "dz_4.m"}], MZero];

Print["Matrices exported successfully."];

(* Verify the exported matrix *)
Print["\n=== Exported dz_0.m ==="];
Print[Import[FileNameJoin[{outputDir, "dz_0.m"}]] // MatrixForm];

(* Boundary conditions at z = z0 *)
z0 = 1/10;  (* Start at z = 0.1 to avoid singularity at z = 0 *)
Print["\n=== Boundary conditions at z = ", z0, " ==="];

y2F1val = N[Hypergeometric2F1[a, b, c, z0], 50];
yPrimeval = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z0, 50];

Print["Y1 = 2F1(", a, ",", b, ";", c, ";", z0, ") = ", NumberForm[y2F1val, 20]];
Print["Y2 = d/dz 2F1 at z=", z0, " = ", NumberForm[yPrimeval, 20]];

(* Test point z = 1/2 *)
Print["\n=== Expected values at z = 1/2 ==="];
valueAtHalf = N[Hypergeometric2F1[a, b, c, 1/2], 50];
derivAtHalf = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 1/2, 50];
Print["2F1(", a, ",", b, ";", c, ";1/2) = ", NumberForm[valueAtHalf, 20]];
Print["d/dz 2F1 at z=1/2 = ", NumberForm[derivAtHalf, 20]];
