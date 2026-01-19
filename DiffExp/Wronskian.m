(* ::Package:: *)

(* DiffExp Wronskian Subpackage *)
(* This package provides Wronskian inverse and related functions *)

BeginPackage["DiffExp`Wronskian`", {"DiffExp`Symbols`", "DiffExp`State`", "DiffExp`Utilities`", "DiffExp`SeriesOps`"}];

(* Wronskian functions *)
MatrixLogxInverse::usage = "MatrixLogxInverse[Mat_] computes inverse of matrix with Logx terms.";
NullSpaceTryAgainOnFail::usage = "NullSpaceTryAgainOnFail[ex_,r___] computes nullspace with retry on fail.";
CombineDifferentialEquationsHomogeneous::usage = "CombineDifferentialEquationsHomogeneous[Amat_,topind_] combines differential equations.";
CombineDifferentialEquationsWithPivotSelection::usage = "CombineDifferentialEquationsWithPivotSelection[Amat_] tries multiple pivot indices to find an invertible Mtilde.";

Begin["`Private`"];

(* Inverse of matrix containing Logx terms *)
MatrixLogxInverse[Mat_] := Module[
  {maxLogPower = (Mat // Dimensions // First) - 1, inverseCoeffs, matrixExpanded, logCoeffs},

  matrixExpanded = Mat // DiffExp`SeriesOps`SExpand;
  logCoeffs[0] = matrixExpanded /. DiffExp`Symbols`Logx -> 0;
  inverseCoeffs[0] = Inverse[logCoeffs[0]];

  Do[
    logCoeffs[mm] = Coefficient[matrixExpanded, DiffExp`Symbols`Logx^mm];
    inverseCoeffs[mm] = -DiffExp`SeriesOps`MatrixMultiplySExpand[inverseCoeffs[0],
      (Sum[logCoeffs[jj] . inverseCoeffs[mm - jj], {jj, 1, mm}])];,
    {mm, 1, maxLogPower}
  ];

  Sum[inverseCoeffs[mm] DiffExp`Symbols`Logx^mm, {mm, 0, maxLogPower}] // DiffExp`SeriesOps`SExpand
];

(* NullSpace with retry on failure *)
NullSpaceTryAgainOnFail[ex_, r___] := Module[{nullSpaceResult},
  nullSpaceResult = NullSpace[ex, r];

  If[Quiet[Length[Cases[nullSpaceResult, SeriesData[DiffExp`Symbols`x, _, List[], k_, k_, _] /; k < 0, Infinity]] > 0],
    DiffExp`Utilities`PrintWarning["Encountered a problem while determining the nullspace of a matrix. Try setting the option \"HomogeneousSolve\" -> \"DontExpand\" "];
    DiffExp`State`LastErrorContext = {ex, r};
    Abort[];
  ];

  nullSpaceResult
];

(* Derives an n-th order differential equation for a single integral *)
CombineDifferentialEquationsHomogeneous[Amat_, topind_: 1] := Module[
  {Amats, matrixSize = Amat // Length, Solns, MtildeMat},

  DiffExp`Utilities`PrintDebug["Getting higher order derivatives.."][2];

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    Amats = DiffExp`Utilities`PChop@NestList[
      If[DiffExp`State`FEC["VerbosityDebug"] >= 4, EchoTiming, Identity][
        DiffExp`SeriesOps`SExpand[(DiffExp`SeriesOps`SD[#, DiffExp`Symbols`x] +
          DiffExp`SeriesOps`MatrixMultiplySExpand[#, Amat])]
      ] &,
      IdentityMatrix[matrixSize], matrixSize
    ];
    Solns = Amats[[All, topind]] // Transpose //
      NullSpaceTryAgainOnFail[#, Method -> "DivisionFreeRowReduction",
        Tolerance -> 10^-DiffExp`State`LinearSolveChopPrecisionVal] &;,

    Amats = DiffExp`Utilities`PChop@NestList[Together[(D[#, DiffExp`Symbols`x] + # . Amat)] &, IdentityMatrix[matrixSize], matrixSize];
    Solns = DiffExp`SeriesOps`DiffExpSeries[Amats[[All, topind]] // Transpose // NullSpace[#] & // Together];
  ];

  If[Length[Solns] > 1,
    DiffExp`Utilities`PrintDebug["Pivot ", topind, " yields singular Mtilde (", Length[Solns], " null space solutions). Trying next pivot..."][2];
    Return[$Failed]
  ];

  If[DiffExp`State`FEC["HomogeneousSolve"] === "Expand",
    MtildeMat = DiffExp`Utilities`PChop@DiffExp`SeriesOps`DiffExpSeries[
      Amats[[All, topind]][[Range[matrixSize], Range[matrixSize]]]];,
    MtildeMat = Amats[[All, topind]][[Range[matrixSize], Range[matrixSize]]];
  ];

  {Solns[[1]]/(Solns[[1]] // Last), MtildeMat}
];

(* Wrapper that tries multiple pivot indices to find an invertible Mtilde *)
(* Returns {HomogeneousEquation, MtildeMat, pivotIndex} on success, or $NeedsFallback if all pivots fail *)
CombineDifferentialEquationsWithPivotSelection[Amat_] := Module[
  {matrixSize = Amat // Length, result, pivotIndex},

  (* Use Catch/Throw for early return from Do loop *)
  Catch[
    (* Try each pivot index from 1 to matrixSize *)
    Do[
      DiffExp`Utilities`PrintDebug["Trying pivot index ", pivotIndex, " of ", matrixSize][3];
      result = CombineDifferentialEquationsHomogeneous[Amat, pivotIndex];

      If[result =!= $Failed,
        DiffExp`Utilities`PrintDebug["Successfully found invertible Mtilde with pivot ", pivotIndex][2];
        Throw[Append[result, pivotIndex], "PivotFound"]
      ];
      , {pivotIndex, matrixSize}
    ];

    (* All pivots failed - return indicator for fallback to VOPAlt *)
    DiffExp`Utilities`PrintInfo["All pivot indices yield singular Mtilde. Falling back to VOPAlt strategy."][2];
    $NeedsFallback

    , "PivotFound"  (* Catch tag *)
  ]
];

End[];

EndPackage[];
