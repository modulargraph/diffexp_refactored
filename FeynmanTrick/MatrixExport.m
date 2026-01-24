(* ::Package:: *)
(* MatrixExport - Export differential matrices in DiffExp format *)

BeginPackage["FeynmanTrick`MatrixExport`", {"FeynmanTrick`"}];

ExportDiffExpMatrix::usage =
  "ExportDiffExpMatrix[matrix, variable, directory, epsOrder] exports a matrix \
with eps dependence to DiffExp order-by-order format: d{variable}_{0..epsOrder}.m. \
The eps symbol is taken from FTConfiguration.";

ExportGeneralMatrix::usage =
  "ExportGeneralMatrix[matrix, variable, directory] exports the full matrix \
(with symbolic eps) to a single file d{variable}_full.m.";

LoadGeneralMatrix::usage =
  "LoadGeneralMatrix[variable, directory] loads a general matrix from \
d{variable}_full.m.";

Begin["`Private`"];

(* ============================================================ *)
(* ExportDiffExpMatrix                                           *)
(* Exports in DiffExp order-by-order format                      *)
(* Files: d{var}_0.m, d{var}_1.m, ..., d{var}_n.m               *)
(* Each file contains a nested list (matrix) with rational       *)
(* function entries in the line parameter variable               *)
(* ============================================================ *)

ExportDiffExpMatrix[matrix_, variable_, directory_String, epsOrder_Integer:4] :=
Module[{eps, varName, fileName, epsMat, k},
  eps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  varName = ToString[variable];

  If[!DirectoryQ[directory],
    CreateDirectory[directory, CreateIntermediateDirectories -> True]
  ];

  Do[
    (* Extract coefficient of eps^k from each matrix entry *)
    epsMat = Map[
      Function[entry,
        SeriesCoefficient[entry, {eps, 0, k}] // Normal // Together
      ],
      matrix,
      {2}
    ];

    fileName = FileNameJoin[{directory,
      "d" <> varName <> "_" <> ToString[k] <> ".m"
    }];

    Put[epsMat, fileName];

    If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 2,
      Print["Exported: ", fileName];
    ];
  , {k, 0, epsOrder}];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Exported DiffExp matrices to: ", directory,
          " (orders 0-", epsOrder, ")"];
  ];
];


(* ============================================================ *)
(* ExportGeneralMatrix                                           *)
(* Exports full matrix with symbolic eps to a single file        *)
(* ============================================================ *)

ExportGeneralMatrix[matrix_, variable_, directory_String] :=
Module[{varName, fileName},
  varName = ToString[variable];

  If[!DirectoryQ[directory],
    CreateDirectory[directory, CreateIntermediateDirectories -> True]
  ];

  fileName = FileNameJoin[{directory,
    "d" <> varName <> "_full.m"
  }];

  Put[matrix, fileName];

  If[FeynmanTrick`Private`$FTConfig["Verbosity"] >= 1,
    Print["Exported general matrix to: ", fileName];
  ];
];


(* ============================================================ *)
(* LoadGeneralMatrix                                             *)
(* Loads a general matrix from file                              *)
(* ============================================================ *)

LoadGeneralMatrix[variable_, directory_String] :=
Module[{varName, fileName},
  varName = ToString[variable];
  fileName = FileNameJoin[{directory,
    "d" <> varName <> "_full.m"
  }];

  If[!FileExistsQ[fileName],
    Print["Error: File not found: ", fileName];
    Return[$Failed];
  ];

  Get[fileName]
];


End[];
EndPackage[];
