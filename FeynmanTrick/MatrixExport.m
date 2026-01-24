(* ::Package:: *)
(* MatrixExport - Export differential matrices in DiffExp format *)

BeginPackage["FeynmanTrick`MatrixExport`", {"FeynmanTrick`"}];

ExportDiffExpMatrix::usage =
  "ExportDiffExpMatrix[matrix, variable, directory, epsOrder] exports a matrix \
with eps dependence to DiffExp order-by-order format: d{variable}_{0..epsOrder}.m. \
The dimension variable d is substituted as d = 4 - 2*eps before expansion.";

ExportGeneralMatrix::usage =
  "ExportGeneralMatrix[matrix, variable, directory] exports the full matrix \
(with symbolic eps) to a single file d{variable}_full.m.";

LoadGeneralMatrix::usage =
  "LoadGeneralMatrix[variable, directory] loads a general matrix from \
d{variable}_full.m.";

Begin["`Private`"];

(* Helper: strip contexts from symbols for clean file output *)
(* Substitutes all non-Global, non-System symbols with Global equivalents *)
stripContexts[expr_] :=
Module[{allSyms, nonGlobal, rules},
  allSyms = Cases[expr, _Symbol, Infinity] // DeleteDuplicates;
  nonGlobal = Select[allSyms,
    (Context[#] =!= "Global`" && Context[#] =!= "System`") &
  ];
  rules = Rule[#, Symbol[SymbolName[#]]] & /@ nonGlobal;
  expr /. rules
];


(* ============================================================ *)
(* ExportDiffExpMatrix                                           *)
(* Exports in DiffExp order-by-order format                      *)
(* Files: d{var}_0.m, d{var}_1.m, ..., d{var}_n.m               *)
(* Each file contains a nested list (matrix) with rational       *)
(* function entries in the line parameter variable               *)
(* ============================================================ *)

ExportDiffExpMatrix[matrix_, variable_, directory_String, epsOrder_Integer:4] :=
Module[{eps, dimVar, varName, fileName, epsMat, epsMatrix, k},
  eps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  varName = If[Head[variable] === Symbol, SymbolName[variable], ToString[variable]];

  If[!DirectoryQ[directory],
    CreateDirectory[directory, CreateIntermediateDirectories -> True]
  ];

  (* Substitute d -> 4 - 2*eps to express matrix in terms of eps *)
  epsMatrix = matrix /. dimVar -> (4 - 2*eps);

  Do[
    (* Extract coefficient of eps^k from each matrix entry *)
    epsMat = Map[
      Function[entry,
        SeriesCoefficient[entry, {eps, 0, k}] // Normal // Together
      ],
      epsMatrix,
      {2}
    ];

    (* Strip contexts for clean output *)
    epsMat = stripContexts[epsMat];

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
Module[{varName, fileName, exportMatrix, eps, dimVar},
  varName = If[Head[variable] === Symbol, SymbolName[variable], ToString[variable]];
  eps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];

  If[!DirectoryQ[directory],
    CreateDirectory[directory, CreateIntermediateDirectories -> True]
  ];

  (* Substitute d -> 4 - 2*eps and strip contexts *)
  exportMatrix = matrix /. dimVar -> (4 - 2*eps);
  exportMatrix = stripContexts[exportMatrix];

  fileName = FileNameJoin[{directory,
    "d" <> varName <> "_full.m"
  }];

  Put[exportMatrix, fileName];

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
  varName = If[Head[variable] === Symbol, SymbolName[variable], ToString[variable]];
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
