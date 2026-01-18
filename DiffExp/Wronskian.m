(* ::Package:: *)

(* DiffExp Wronskian Subpackage *)
(* This package provides Wronskian inverse and related functions *)

BeginPackage["DiffExp`Wronskian`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Wronskian functions *)
MatrixLogxInverse::usage = "MatrixLogxInverse[Mat_] computes inverse of matrix with Logx terms.";
NullSpaceTryAgainOnFail::usage = "NullSpaceTryAgainOnFail[ex_,r___] computes nullspace with retry on fail.";
CombineDifferentialEquationsHomogeneous::usage = "CombineDifferentialEquationsHomogeneous[Amat_,topind_] combines differential equations.";

Begin["`Private`"];

(* Inverse of matrix containing Logx terms *)
MatrixLogxInverse[Mat_] := Module[
  {MaxLogxPower = (Mat // Dimensions // First) - 1, BB, AAA, AA},

  AAA = Mat // DiffExp`SeriesOps`SExpand;
  AA[0] = AAA /. DiffExp`Symbols`Logx -> 0;
  BB[0] = Inverse[AA[0]];

  Do[
    AA[mm] = Coefficient[AAA, DiffExp`Symbols`Logx^mm];
    BB[mm] = -DiffExp`SeriesOps`MatrixMultiplySExpand[BB[0],
      (Sum[AA[jj] . BB[mm - jj], {jj, 1, mm}])];,
    {mm, 1, MaxLogxPower}
  ];

  Sum[BB[mm] DiffExp`Symbols`Logx^mm, {mm, 0, MaxLogxPower}] // DiffExp`SeriesOps`SExpand
];

(* NullSpace with retry on failure *)
NullSpaceTryAgainOnFail[ex_, r___] := Module[{Res, ValidateTerm},
  Res = NullSpace[ex, r];

  If[Quiet[Length[Cases[Res, SeriesData[DiffExp`Symbols`x, _, List[], k_, k_, _] /; k < 0, Infinity]] > 0],
    DiffExp`Utilities`PrintWarning["Encountered a problem while determining the nullspace of a matrix. Try setting the option \"HomogeneousSolve\" -> \"DontExpand\" "];
    Global`DebugData = {ex, r};
    Abort[];
  ];

  Res
];

(* Derives an n-th order differential equation for a single integral *)
CombineDifferentialEquationsHomogeneous[Amat_, topind_: 1] := Module[
  {Amats, n = Amat // Length, Solns, MtildeMat},

  DiffExp`Utilities`PrintDebug["Getting higher order derivatives.."][2];

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    Amats = DiffExp`Utilities`PChop@NestList[
      If[DiffExp`State`FEC["VerbosityDebug"] >= 4, EchoTiming, Identity][
        DiffExp`SeriesOps`SExpand[(DiffExp`SeriesOps`SD[#, DiffExp`Symbols`x] +
          DiffExp`SeriesOps`MatrixMultiplySExpand[#, Amat])]
      ] &,
      IdentityMatrix[n], n
    ];
    Solns = Amats[[All, topind]] // Transpose //
      NullSpaceTryAgainOnFail[#, Method -> "DivisionFreeRowReduction",
        Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;,

    Amats = DiffExp`Utilities`PChop@NestList[Together[(D[#, DiffExp`Symbols`x] + # . Amat)] &, IdentityMatrix[n], n];
    Solns = DiffExp`SeriesOps`DiffExpSeries[Amats[[All, topind]] // Transpose // NullSpace[#] & // Together];
  ];

  If[Length[Solns] > 1,
    Global`DebugData = {Amats[[All, topind]]};
    DiffExp`Utilities`ReportError["Found multiple solutions while combining the differential equations. Use the option IntegrationStrategy -> \"VariationOfParameters\", or choose a different line."]
  ];

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    MtildeMat = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[
      Amats[[All, topind]][[Range[n], Range[n]]]];,
    MtildeMat = Amats[[All, topind]][[Range[n], Range[n]]];
  ];

  {Solns[[1]]/(Solns[[1]] // Last), MtildeMat}
];

End[];

EndPackage[];
