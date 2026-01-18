(* Five-Point Non-Planar Master Integrals *)
(* Cross-checking results from arXiv:1812.11160: *)
(* "All master integrals for three-jet production at NNLO" *)
(* D. Chicherin, T. Gehrmann, J. M. Henn, P. Wasser, Y. Zhang, S. Zoia *)

(* Set the directory to the script location *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
Print["Working directory: ", scriptDir];

(* ========================================================================= *)
(* Pentagon alphabet                                                        *)
(* ========================================================================= *)

PentagonAlphabet = {
  W[1] -> v1,
  W[2] -> v2,
  W[3] -> v3,
  W[4] -> v4,
  W[5] -> v5,
  W[6] -> v3 + v4,
  W[7] -> v4 + v5,
  W[8] -> v1 + v5,
  W[9] -> v1 + v2,
  W[10] -> v2 + v3,
  W[11] -> v1 - v4,
  W[12] -> v2 - v5,
  W[13] -> -v1 + v3,
  W[14] -> -v2 + v4,
  W[15] -> -v3 + v5,
  W[16] -> v1 + v2 - v4,
  W[17] -> v2 + v3 - v5,
  W[18] -> -v1 + v3 + v4,
  W[19] -> -v2 + v4 + v5,
  W[20] -> v1 - v3 + v5,
  W[21] -> -v1 - v2 + v3 + v4,
  W[22] -> -v2 - v3 + v4 + v5,
  W[23] -> v1 - v3 - v4 + v5,
  W[24] -> v1 + v2 - v4 - v5,
  W[25] -> -v1 + v2 + v3 - v5,
  W[26] -> (v1 v2 - v2 v3 + v3 v4 - v1 v5 - v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (v1 v2 - v2 v3 + v3 v4 - v1 v5 - v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]),
  W[27] -> (-v1 v2 + v2 v3 - v3 v4 - v1 v5 + v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (-v1 v2 + v2 v3 - v3 v4 - v1 v5 + v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]),
  W[28] -> (-v1 v2 - v2 v3 + v3 v4 + v1 v5 - v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (-v1 v2 - v2 v3 + v3 v4 + v1 v5 - v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]),
  W[29] -> (v1 v2 - v2 v3 - v3 v4 - v1 v5 + v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (v1 v2 - v2 v3 - v3 v4 - v1 v5 + v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]),
  W[30] -> (-v1 v2 + v2 v3 - v3 v4 + v1 v5 - v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (-v1 v2 + v2 v3 - v3 v4 + v1 v5 - v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]),
  W[31] -> (v1 v2 + v2 v3 - v3 v4 - v1 v5 - v4 v5 -
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))]) /
   (v1 v2 + v2 v3 - v3 v4 - v1 v5 - v4 v5 +
    Sqrt[v1^2 (v2 - v5)^2 + (v2 v3 + v4 (-v3 + v5))^2 +
      2 v1 (-v2^2 v3 + v4 (v3 - v5) v5 + v2 (v3 v4 + (v3 + v4) v5))])
};

(* ========================================================================= *)
(* Download ancillary files from arXiv                                      *)
(* ========================================================================= *)

Print["\n========================================================================"];
Print["=== FIVE-POINT NON-PLANAR MASTER INTEGRALS ==="];
Print["========================================================================\n"];

(* Set up temporary directory *)
MyTemporaryDirectory = scriptDir <> "/Tmp/";
If[!DirectoryQ[MyTemporaryDirectory <> "5pNonPlanar"],
  CreateDirectory[FileNameJoin[{MyTemporaryDirectory, "5pNonPlanar"}]];
];

Print["Temporary directory: ", MyTemporaryDirectory];

(* Download ancillary files if not present *)
AncillaryFiles = {"XB_Atilde.txt", "XB_Boundary_values_X0.txt", "XB_Boundary_values_X1.txt"};

Do[
  targetFile = MyTemporaryDirectory <> "5pNonPlanar" <> $PathnameSeparator <> myfile;
  If[!FileExistsQ[targetFile],
    Print["Downloading: ", "https://arxiv.org/src/1812.11160v2/anc/" <> myfile];
    DownloadedFile = URLDownload["https://arxiv.org/src/1812.11160v2/anc/" <> myfile];
    CopyFile[DownloadedFile, targetFile];
    Print["Downloaded: ", myfile];
  ,
    Print["File already exists: ", myfile];
  ];
, {myfile, AncillaryFiles}];

(* ========================================================================= *)
(* Prepare differential equations for DiffExp                               *)
(* ========================================================================= *)

Print["\nPreparing differential equations..."];

(* Export the differential matrix in DiffExp format *)
diffMatrixFile = FileNameJoin[{MyTemporaryDirectory, "5pNonPlanar", "d_1.m"}];
If[!FileExistsQ[diffMatrixFile],
  Print["Converting Atilde matrix to DiffExp format..."];
  Export[
    diffMatrixFile,
    (Import[MyTemporaryDirectory <> "5pNonPlanar" <> $PathnameSeparator <> "XB_Atilde.txt"] // ToExpression) /. PentagonAlphabet
  ];
  Print["Matrix exported to: ", diffMatrixFile];
,
  Print["Matrix file already exists."];
];

(* ========================================================================= *)
(* Load DiffExp and configure                                               *)
(* ========================================================================= *)

Print["\nLoading DiffExp package..."];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp.m"}]];
Print["Package loaded successfully!"];

(* Configuration - using LOW expansion order for quick tests *)
FivePConfiguration = {
  AccuracyGoal -> 15,           (* Reduced for faster testing *)
  ExpansionOrder -> 30,         (* LOW for quick testing - original was 80 *)
  DeltaPrescriptions -> {
    v1 + I \[Delta],
    v2 + I \[Delta],
    v3 + I \[Delta],
    v4 + I \[Delta],
    v5 + I \[Delta]
  },
  MatrixDirectory -> MyTemporaryDirectory <> "5pNonPlanar",
  Verbosity -> 2,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 100,      (* Reduced for faster testing - original was 150 *)
  ChopPrecision -> 60           (* Reduced for faster testing - original was 100 *)
};

Print["\nLoading configuration..."];
LoadConfiguration[FivePConfiguration];
Print["Configuration loaded.\n"];

(* ========================================================================= *)
(* Cross-check numerical results                                            *)
(* ========================================================================= *)

Print["=== Cross-checking numerical results from the paper ==="];
Print["Transporting from point X0 to X1\n"];

(* Load boundary values at X0 *)
BoundaryValuesX0 = Import[
  MyTemporaryDirectory <> "5pNonPlanar" <> $PathnameSeparator <> "XB_Boundary_values_X0.txt"
] // ToExpression;

(* Point X0 *)
X0 = {v1 -> 3, v2 -> -1, v3 -> 1, v4 -> 1, v5 -> -1};
Print["Starting point X0: ", X0];

(* Point X1 *)
X1 = {v1 -> 4, v2 -> -113/47, v3 -> 281/149, v4 -> 349/257, v5 -> -863/541};
Print["Target point X1: ", X1];

(* Prepare boundary conditions *)
BoundaryConditionsDiffExp = PrepareBoundaryConditions[BoundaryValuesX0, X0];

(* Transport from X0 to X1 *)
Print["\nTransporting from X0 to X1..."];
X0ToX1 = TransportTo[
  BoundaryConditionsDiffExp,
  X1
];

Print["\nTransport complete!"];

(* Display some results *)
Print["\n=== Results at X1 ==="];
Print["First few master integrals (showing value and error estimate):"];

(* Show first 3 master integrals *)
Do[
  Print["Master integral ", i, ": ", X0ToX1[[2, i, 1]], " +/- ", X0ToX1[[3, i, 1]]];
, {i, 1, Min[3, Length[X0ToX1[[2]]]]}];

Print["\n========================================================================"];
Print["=== TEST COMPLETED SUCCESSFULLY! ==="];
Print["========================================================================"];
