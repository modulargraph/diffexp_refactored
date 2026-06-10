(* Print the exact (eps-dependent) residue-matrix eigenvalues of an exported
   level matrix at given singular points: rho_i(eps) = a_i + b_i*eps.
   These are the true local sector exponents the transport must realize. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

matFile = envOrDefault["MAT_FILE",
  "/var/folders/r9/t81xsghn0r198pktdqz_0v700000gn/T/FT_stepwise_banana_82095/Level_1_Matrices/dxx1_full.m"];
pointsRaw = envOrDefault["POINTS", "0,1"];
points = ToExpression /@ StringSplit[pointsRaw, ","];

mat = Get[matFile];
Print["matrix dims=", Dimensions[mat]];
vars = Variables[Level[mat, {-1}]];
Print["vars=", vars];

(* identify the line variable (xx1-like) and eps symbol *)
xVar = SelectFirst[vars, StringMatchQ[SymbolName[#], "xx" ~~ DigitCharacter ...] &, First[vars]];
epsCandidates = Select[vars, StringMatchQ[SymbolName[#], "eps" | "FTeps" | "ep"] &];
dCandidates = Select[vars, SymbolName[#] === "d" &];
Print["xVar=", xVar, " epsCandidates=", epsCandidates, " dCandidates=", dCandidates];

eps = Symbol["epsLocal"];
matEps = If[Length[dCandidates] > 0,
  mat /. dCandidates[[1]] -> 2 - 2*eps,
  If[Length[epsCandidates] > 0, mat /. epsCandidates[[1]] -> eps, mat]
];

Do[
  Module[{doublePole, residue, evs},
    doublePole = Map[
      Function[entry, SeriesCoefficient[entry, {xVar, pt, -2}]],
      matEps, {2}
    ] // Simplify;
    residue = Map[
      Function[entry, SeriesCoefficient[entry, {xVar, pt, -1}]],
      matEps, {2}
    ] // Simplify;
    If[!TrueQ[doublePole == ConstantArray[0, Dimensions[matEps]]],
      Print["POINT ", pt, " has DOUBLE POLE content: ",
        InputForm[doublePole /. eps -> Symbol["e"]]];
    ];
    evs = Sort[Eigenvalues[residue] // Simplify // Factor];
    Print["POINT ", pt, " exact exponents rho_i(eps) = ", InputForm[evs]];
  ],
  {pt, points}
];

Quit[0];
