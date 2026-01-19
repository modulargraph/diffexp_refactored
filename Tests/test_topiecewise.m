(* Test ToPiecewise: compare piecewise evaluation vs manual transport *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Testing ToPiecewise: Piecewise vs Manual Transport"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading refactored DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Equal mass banana configuration *)
Print["=== Equal Mass Banana Configuration ==="];
EqualMassConfiguration = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True
};

Print["Loading configuration..."];
LoadConfiguration[EqualMassConfiguration];

(* Boundary conditions *)
Print["\n=== Setting up boundary conditions ==="];
EqualMassBoundaryConditions = {
  "?",
  "?",
  \[Epsilon] (1 + 3 \[Epsilon]) (1 + 4 \[Epsilon]) * (
    -4 E^(3 EulerGamma \[Epsilon]) Gamma[\[Epsilon]]^3/t +
    6 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^2 Gamma[\[Epsilon]]^3/Gamma[-2 \[Epsilon]] +
    8 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 2 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^3 Gamma[\[Epsilon]] Gamma[2 \[Epsilon]]/Gamma[-3 \[Epsilon]] +
    3 E^(3 EulerGamma \[Epsilon]) (-1/t)^(1 + 3 \[Epsilon]) \[Epsilon] Gamma[-\[Epsilon]]^4 Gamma[3 \[Epsilon]]/Gamma[-4 \[Epsilon]]
  ),
  E^(3 EulerGamma \[Epsilon]) \[Epsilon]^3 Gamma[\[Epsilon]]^3
};

Print["Preparing boundary conditions..."];
PreparedBCs = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
Print["Boundary conditions prepared\n"];

(* Transport to t = -1 first *)
Print["=== Transporting to t = -1 ==="];
ResultsAtMinus1 = TransportTo[PreparedBCs, <|t -> -1|>];
Print["Transport to t = -1 complete\n"];

(* ============================================================ *)
(* Transport from t = -1 to t = 5 with SaveExpansions = True *)
(* ============================================================ *)
Print["=== Transporting from t = -1 to t = 5 (with SaveExpansions=True) ==="];
ResultsWithSave = TransportTo[ResultsAtMinus1, <|t -> x|>, 5, True];
Print["Transport complete!\n"];

(* ============================================================ *)
(* Use ToPiecewise to create piecewise functions *)
(* ============================================================ *)
Print["=== Creating piecewise functions with ToPiecewise ==="];
PiecewiseFunctions = ToPiecewise[ResultsWithSave];
Print["Piecewise functions created."];
Print["Dimensions: ", Dimensions[PiecewiseFunctions], "\n"];

(* ============================================================ *)
(* Manually transport to t = 3 *)
(* ============================================================ *)
Print["=== Manually transporting to t = 3 ==="];
ManualResultsAt3 = TransportTo[ResultsAtMinus1, <|t -> x|>, 3];
Print["Manual transport to t = 3 complete."];
Print["Result point: ", ManualResultsAt3[[1]], "\n"];

(* ============================================================ *)
(* Compare: Evaluate piecewise at t = 3 vs manual transport *)
(* ============================================================ *)
Print["=== Comparing piecewise evaluation vs manual transport at t = 3 ===\n"];

evalPoint = 3;
testsPassed = 0;
testsTotal = 0;

(* Evaluate piecewise functions at t = 3 *)
Print["Evaluating piecewise functions at t = ", evalPoint, "..."];
PiecewiseValuesAt3 = Table[
  PiecewiseFunctions[[i, j]][evalPoint],
  {i, Dimensions[PiecewiseFunctions][[1]]},
  {j, Dimensions[PiecewiseFunctions][[2]]}
];

(* Get manual transport values *)
ManualValuesAt3 = ManualResultsAt3[[2]];

Print["\n--- Comparison of values ---"];
Print["(Integral index, Epsilon order): Piecewise vs Manual"];

maxRelDiff = 0;
maxAbsDiff = 0;
Do[
  piecewiseVal = PiecewiseValuesAt3[[i, j]] // N[#, 20] &;
  manualVal = ManualValuesAt3[[i, j]] // N[#, 20] &;

  absDiff = Abs[piecewiseVal - manualVal];
  If[absDiff > maxAbsDiff, maxAbsDiff = absDiff];

  If[Abs[manualVal] > 10^-50,
    relDiff = absDiff / Abs[manualVal];
    If[relDiff > maxRelDiff, maxRelDiff = relDiff];
  ];

  , {i, Length[PiecewiseValuesAt3]}, {j, Length[PiecewiseValuesAt3[[1]]]}
];

Print["\nMaximum absolute difference: ", N[maxAbsDiff, 5]];
Print["Maximum relative difference: ", N[maxRelDiff, 5]];

(* Test 1: Check that piecewise functions return numerical values *)
testsTotal++;
numericCheck = And @@ Flatten[Table[
  NumericQ[N[PiecewiseValuesAt3[[i, j]]]],
  {i, Length[PiecewiseValuesAt3]}, {j, Length[PiecewiseValuesAt3[[1]]]}
]];
If[numericCheck,
  Print["\n  [PASS] Piecewise functions return numerical values"];
  testsPassed++;
  ,
  Print["\n  [FAIL] Piecewise functions contain non-numerical values"];
];

(* Test 2: Check that the differences are small (< 10^-10) *)
testsTotal++;
If[maxRelDiff < 10^-10,
  Print["  [PASS] Piecewise values match manual transport (rel. diff < 10^-10)"];
  testsPassed++;
  ,
  Print["  [FAIL] Piecewise values differ from manual transport"];
  Print["\n  Piecewise values at t=3:"];
  Print[TableForm[N[PiecewiseValuesAt3, 10]]];
  Print["\n  Manual transport values at t=3:"];
  Print[TableForm[N[ManualValuesAt3, 10]]];
];

(* Test 3: Check dimensions match *)
testsTotal++;
If[Dimensions[PiecewiseValuesAt3] === Dimensions[ManualValuesAt3],
  Print["  [PASS] Dimensions match: ", Dimensions[PiecewiseValuesAt3]];
  testsPassed++;
  ,
  Print["  [FAIL] Dimensions mismatch"];
  Print["    Piecewise: ", Dimensions[PiecewiseValuesAt3]];
  Print["    Manual: ", Dimensions[ManualValuesAt3]];
];

(* Test 4: Test at another point within the range (t = 2) *)
(* Note: t=0 and t=4 are singularities, so we use t=2 which is safe *)
Print["\n=== Testing at t = 2 (non-singular point) ==="];
evalPoint2 = 2;

Print["Manually transporting to t = ", evalPoint2, "..."];
ManualResultsAt2 = TransportTo[ResultsAtMinus1, <|t -> x|>, evalPoint2];
ManualValuesAt2 = ManualResultsAt2[[2]];

PiecewiseValuesAt2 = Table[
  PiecewiseFunctions[[i, j]][evalPoint2],
  {i, Dimensions[PiecewiseFunctions][[1]]},
  {j, Dimensions[PiecewiseFunctions][[2]]}
];

maxRelDiff2 = 0;
Do[
  piecewiseVal = PiecewiseValuesAt2[[i, j]] // N[#, 20] &;
  manualVal = ManualValuesAt2[[i, j]] // N[#, 20] &;

  If[Abs[manualVal] > 10^-50,
    relDiff = Abs[(piecewiseVal - manualVal) / manualVal];
    If[NumericQ[relDiff] && relDiff > maxRelDiff2, maxRelDiff2 = relDiff];
  ];

  , {i, Length[PiecewiseValuesAt2]}, {j, Length[PiecewiseValuesAt2[[1]]]}
];

Print["Maximum relative difference at t=2: ", N[maxRelDiff2, 5]];

testsTotal++;
If[maxRelDiff2 < 10^-10,
  Print["  [PASS] Piecewise values match manual transport at t=2 (rel. diff < 10^-10)"];
  testsPassed++;
  ,
  Print["  [FAIL] Piecewise values differ from manual transport at t=2"];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
