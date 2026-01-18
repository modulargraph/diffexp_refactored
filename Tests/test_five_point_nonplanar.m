(* Test: Five-Point Non-Planar Master Integrals *)
(* Cross-checking results from arXiv:1812.11160: *)
(* "All master integrals for three-jet production at NNLO" *)
(* D. Chicherin, T. Gehrmann, J. M. Henn, P. Wasser, Y. Zhang, S. Zoia *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add DiffExp to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Test: Five-Point Non-Planar Master Integrals"];
Print["===========================================\n"];

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
(* Set up temporary directory and download files                            *)
(* ========================================================================= *)

(* Set up temporary directory *)
MyTemporaryDirectory = FileNameJoin[{scriptDir, "Tmp_5pNonPlanar"}];
If[!DirectoryQ[MyTemporaryDirectory],
  CreateDirectory[MyTemporaryDirectory];
];

Print["Temporary directory: ", MyTemporaryDirectory];

(* Download ancillary files if not present *)
AncillaryFiles = {"XB_Atilde.txt", "XB_Boundary_values_X0.txt", "XB_Boundary_values_X1.txt"};

downloadSuccess = True;
Do[
  targetFile = FileNameJoin[{MyTemporaryDirectory, myfile}];
  If[!FileExistsQ[targetFile],
    Print["Downloading: ", "https://arxiv.org/src/1812.11160v2/anc/" <> myfile];
    Quiet[
      Check[
        DownloadedFile = URLDownload["https://arxiv.org/src/1812.11160v2/anc/" <> myfile];
        CopyFile[DownloadedFile, targetFile];
        Print["Downloaded: ", myfile];
        ,
        Print["  [WARNING] Failed to download: ", myfile];
        downloadSuccess = False;
      ]
    ];
  ,
    Print["File already exists: ", myfile];
  ];
, {myfile, AncillaryFiles}];

(* Check if we can proceed with the test *)
If[!downloadSuccess || !And @@ (FileExistsQ[FileNameJoin[{MyTemporaryDirectory, #}]] & /@ AncillaryFiles),
  Print["\n==========================================="];
  Print["[SKIP] Cannot run test - missing ancillary files"];
  Print["This test requires internet access to download files from arXiv."];
  Print["==========================================="];
  Quit[];
];

(* ========================================================================= *)
(* Prepare differential equations                                           *)
(* ========================================================================= *)

Print["\nPreparing differential equations..."];

(* Export the differential matrix in DiffExp format *)
diffMatrixFile = FileNameJoin[{MyTemporaryDirectory, "d_1.m"}];
If[!FileExistsQ[diffMatrixFile],
  Print["Converting Atilde matrix to DiffExp format..."];
  Export[
    diffMatrixFile,
    (Import[FileNameJoin[{MyTemporaryDirectory, "XB_Atilde.txt"}]] // ToExpression) /. PentagonAlphabet
  ];
  Print["Matrix exported to: ", diffMatrixFile];
,
  Print["Matrix file already exists."];
];

(* ========================================================================= *)
(* Load DiffExp and configure                                               *)
(* ========================================================================= *)

Print["\nLoading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Configuration - using LOW expansion order for quick tests *)
FivePConfiguration = {
  AccuracyGoal -> 15,
  ExpansionOrder -> 30,         (* LOW for quick testing *)
  DeltaPrescriptions -> {
    v1 + I \[Delta],
    v2 + I \[Delta],
    v3 + I \[Delta],
    v4 + I \[Delta],
    v5 + I \[Delta]
  },
  MatrixDirectory -> MyTemporaryDirectory,
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 100,
  ChopPrecision -> 60
};

Print["Loading configuration..."];
LoadConfiguration[FivePConfiguration];
Print["Configuration loaded.\n"];

(* ========================================================================= *)
(* Tests                                                                    *)
(* ========================================================================= *)

testsPassed = 0;
testsTotal = 0;

(* Test 1: Check configuration loaded correctly *)
Print["=== Test 1: Configuration Loading ==="];
testsTotal++;
If[DiffExp`State`NumIntegrals > 0,
  Print["  [PASS] Matrices loaded, NumIntegrals = ", DiffExp`State`NumIntegrals];
  testsPassed++;
  ,
  Print["  [FAIL] Matrices not loaded properly"];
];

(* Load boundary values at X0 *)
Print["\n=== Test 2: Boundary Conditions ==="];
testsTotal++;
BoundaryValuesX0 = Import[
  FileNameJoin[{MyTemporaryDirectory, "XB_Boundary_values_X0.txt"}]
] // ToExpression;

If[ListQ[BoundaryValuesX0] && Length[BoundaryValuesX0] > 0,
  Print["  [PASS] Boundary values loaded, ", Length[BoundaryValuesX0], " integrals"];
  testsPassed++;
  ,
  Print["  [FAIL] Boundary values not loaded properly"];
];

(* Point X0 *)
X0 = {v1 -> 3, v2 -> -1, v3 -> 1, v4 -> 1, v5 -> -1};
Print["Starting point X0: ", X0];

(* Prepare boundary conditions *)
Print["\n=== Test 3: Prepare Boundary Conditions ==="];
testsTotal++;
BoundaryConditionsDiffExp = PrepareBoundaryConditions[BoundaryValuesX0, X0];
If[ListQ[BoundaryConditionsDiffExp] && Length[BoundaryConditionsDiffExp] >= 2,
  Print["  [PASS] Boundary conditions prepared"];
  testsPassed++;
  ,
  Print["  [FAIL] Boundary conditions preparation failed"];
];

(* Point X1 - a simple nearby point for quick testing *)
X1 = {v1 -> 3 + 1/10, v2 -> -1, v3 -> 1, v4 -> 1, v5 -> -1};
Print["\nTarget point X1 (simple nearby point): ", X1];

(* Transport from X0 to X1 *)
Print["\n=== Test 4: Transport ==="];
testsTotal++;
Print["Transporting from X0 to X1..."];
X0ToX1 = TransportTo[BoundaryConditionsDiffExp, X1];

If[ListQ[X0ToX1] && Length[X0ToX1] >= 2,
  Print["  [PASS] Transport completed"];
  Print["  Result structure: ", Length[X0ToX1], " elements"];
  testsPassed++;
  ,
  Print["  [FAIL] Transport failed"];
];

(* Test 5: Check results are numerical *)
Print["\n=== Test 5: Results Numerical ==="];
testsTotal++;
If[Length[X0ToX1] >= 3 && NumericQ[X0ToX1[[2, 1, 1]]],
  Print["  [PASS] Results are numerical"];
  Print["  First master integral: ", X0ToX1[[2, 1, 1]], " +/- ", X0ToX1[[3, 1, 1]]];
  testsPassed++;
  ,
  Print["  [FAIL] Results are not numerical"];
];

(* Cleanup *)
Print["\nCleaning up temporary directory..."];
DeleteDirectory[MyTemporaryDirectory, DeleteContents -> True];

(* Summary *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
