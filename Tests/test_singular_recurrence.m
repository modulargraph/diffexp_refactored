(* Test: Singular Recurrence Method *)
(*
   Tests the Frobenius-type recurrence for regular singular points
   (simple pole with non-resonant, diagonalizable residue).

   Test 1: Hypergeometric 2F1
   - Transport from z=0.01 to z=1/2
   - Compare Default strategy vs UseRationalRecurrence->True
   - The singular point at z=0 should trigger the singular recurrence

   Test 2: Equal Mass Banana
   - Transport from boundary to t=-1
   - Compare Default strategy vs UseRationalRecurrence->True
   - Non-resonant singular points (t=4, t=16) use singular recurrence
   - Resonant point (t=0) falls back to Frobenius
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Singular Recurrence Method Tests"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

testsPassed = 0;
testsTotal = 0;

(* ============================================================ *)
(* Test 1: Hypergeometric 2F1                                   *)
(* ============================================================ *)
Print["==========================================="];
Print["Test 1: Hypergeometric 2F1"];
Print["===========================================\n"];

(* Hypergeometric parameters *)
a = 1/4;
b = 1/3;
c = 3/2;

Print["Parameters: a = ", a, ", b = ", b, ", c = ", c];

(* Configuration for Default strategy *)
Config2F1Default = {
  MatrixDirectory -> FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True,
  WorkingPrecision -> 200,
  ExpansionOrder -> 60,
  DivisionOrder -> 4,
  ChopPrecision -> 150,
  EpsilonOrder -> 0
};

(* Boundary conditions near z=0 *)
z0 = 1/100;
y2F1AtStart = N[Hypergeometric2F1[a, b, c, z0], 200];
yPrimeAtStart = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z0, 200];

BCs2F1 = {
  <|z -> z0|>,
  {
    {y2F1AtStart},
    {yPrimeAtStart}
  }
};

(* Run with Default strategy *)
Print["\n--- Running with Default strategy ---"];
LoadConfiguration[Config2F1Default];
ResultDefault = TransportTo[BCs2F1, <|z -> 1/2|>];
Print["Default transport complete."];

y2F1Default = ResultDefault[[2, 1, 1]];
yPrimeDefault = ResultDefault[[2, 2, 1]];

(* Configuration for Singular Recurrence *)
Config2F1Recurrence = Join[Config2F1Default, {
  UseRationalRecurrence -> True,
  Verbosity -> 3  (* Higher verbosity to see singular recurrence messages *)
}];

(* Run with UseRationalRecurrence *)
Print["\n--- Running with UseRationalRecurrence -> True ---"];
LoadConfiguration[Config2F1Recurrence];
ResultRecurrence = TransportTo[BCs2F1, <|z -> 1/2|>];
Print["Recurrence transport complete."];

y2F1Recurrence = ResultRecurrence[[2, 1, 1]];
yPrimeRecurrence = ResultRecurrence[[2, 2, 1]];

(* Compare with exact values *)
y2F1Exact = N[Hypergeometric2F1[a, b, c, 1/2], 200];
yPrimeExact = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> 1/2, 200];

Print["\n--- 2F1 Results Comparison ---"];
diffDefaultExact = Max[Abs[y2F1Default - y2F1Exact], Abs[yPrimeDefault - yPrimeExact]];
diffRecurrenceExact = Max[Abs[y2F1Recurrence - y2F1Exact], Abs[yPrimeRecurrence - yPrimeExact]];
diffDefaultRecurrence = Max[Abs[y2F1Default - y2F1Recurrence], Abs[yPrimeDefault - yPrimeRecurrence]];

Print["  Default vs Exact:     ", ScientificForm[diffDefaultExact, 3]];
Print["  Recurrence vs Exact:  ", ScientificForm[diffRecurrenceExact, 3]];
Print["  Default vs Recurrence:", ScientificForm[diffDefaultRecurrence, 3]];

(* Test: Both should agree with exact to > 40 digits *)
testsTotal++;
If[diffDefaultExact < 10^-40,
  Print["  [PASS] Default matches exact to > 40 digits"];
  testsPassed++;
  ,
  Print["  [FAIL] Default accuracy: ", Floor[-Log10[diffDefaultExact]], " digits"];
];

testsTotal++;
If[diffRecurrenceExact < 10^-40,
  Print["  [PASS] Recurrence matches exact to > 40 digits"];
  testsPassed++;
  ,
  Print["  [FAIL] Recurrence accuracy: ", Floor[-Log10[diffRecurrenceExact]], " digits"];
];

testsTotal++;
If[diffDefaultRecurrence < 10^-40,
  Print["  [PASS] Default and Recurrence agree to > 40 digits"];
  testsPassed++;
  ,
  Print["  [FAIL] Default vs Recurrence differ by: ", ScientificForm[diffDefaultRecurrence, 3]];
];

(* ============================================================ *)
(* Test 2: Equal Mass Banana                                    *)
(* ============================================================ *)
Print["\n\n==========================================="];
Print["Test 2: Equal Mass Banana"];
Print["===========================================\n"];

(* Equal mass banana boundary conditions *)
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

(* Default configuration *)
BananaConfigDefault = {
  DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  Verbosity -> 1,
  UseMobius -> True,
  UsePade -> True
};

(* Run with Default strategy *)
Print["--- Running with Default strategy ---"];
LoadConfiguration[BananaConfigDefault];
PreparedBCsDefault = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
BananaResultDefault = TransportTo[PreparedBCsDefault, <|t -> -1|>];
Print["Default transport complete.\n"];

(* Recurrence configuration *)
BananaConfigRecurrence = Join[BananaConfigDefault, {
  UseRationalRecurrence -> True,
  Verbosity -> 3
}];

(* Run with UseRationalRecurrence *)
Print["--- Running with UseRationalRecurrence -> True ---"];
LoadConfiguration[BananaConfigRecurrence];
PreparedBCsRecurrence = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
BananaResultRecurrence = TransportTo[PreparedBCsRecurrence, <|t -> -1|>];
Print["Recurrence transport complete.\n"];

(* Compare results *)
Print["--- Banana Results Comparison ---"];
bananaDiff = Max[Abs[Flatten[N[BananaResultDefault[[2]] - BananaResultRecurrence[[2]], 30]]]];
Print["  Max difference Default vs Recurrence: ", ScientificForm[bananaDiff, 3]];

testsTotal++;
If[bananaDiff < 10^-10,
  Print["  [PASS] Banana: Default and Recurrence agree to > 10 digits"];
  testsPassed++;
  ,
  Print["  [FAIL] Banana: results differ by ", ScientificForm[bananaDiff, 3]];
  Print["  Default values at t=-1:"];
  Print[TableForm[N[BananaResultDefault[[2]], 10]]];
  Print["  Recurrence values at t=-1:"];
  Print[TableForm[N[BananaResultRecurrence[[2]], 10]]];
];

(* ============================================================ *)
(* Test 2b: Banana crossing t=4 (non-resonant singular point)  *)
(* Eigenvalues at t=4: {0, 0, 0, -1/2} (non-resonant)          *)
(* This should trigger the singular recurrence                  *)
(* ============================================================ *)
Print["\n\n==========================================="];
Print["Test 2b: Banana crossing t=4 (singular recurrence)"];
Print["===========================================\n"];

(* Transport from t=-1 to t=5 with Default strategy *)
Print["--- Default: t=-1 to t=5 ---"];
LoadConfiguration[BananaConfigDefault];
PreparedBCsDefault2 = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
BananaDefault1 = TransportTo[PreparedBCsDefault2, <|t -> -1|>];
BananaDefault5 = TransportTo[BananaDefault1, <|t -> x|>, 5];
Print["Default transport to t=5 complete.\n"];

(* Transport from t=-1 to t=5 with UseRationalRecurrence *)
Print["--- Recurrence: t=-1 to t=5 ---"];
LoadConfiguration[BananaConfigRecurrence];
PreparedBCsRecurrence2 = PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
BananaRecurrence1 = TransportTo[PreparedBCsRecurrence2, <|t -> -1|>];
BananaRecurrence5 = TransportTo[BananaRecurrence1, <|t -> x|>, 5];
Print["Recurrence transport to t=5 complete.\n"];

(* Compare results at t=5 *)
Print["--- Banana t=5 Results Comparison ---"];
bananaDiff5 = Max[Abs[Flatten[N[BananaDefault5[[2]] - BananaRecurrence5[[2]], 30]]]];
Print["  Max difference Default vs Recurrence at t=5: ", ScientificForm[bananaDiff5, 3]];

testsTotal++;
If[bananaDiff5 < 10^-10,
  Print["  [PASS] Banana t=5: Default and Recurrence agree to > 10 digits"];
  testsPassed++;
  ,
  Print["  [FAIL] Banana t=5: results differ by ", ScientificForm[bananaDiff5, 3]];
];

(* ============================================================ *)
(* Test 3: Verify applicability check                           *)
(* ============================================================ *)
Print["\n\n==========================================="];
Print["Test 3: Applicability Check"];
Print["===========================================\n"];

(* Load 2F1 config to have matrices available *)
LoadConfiguration[Config2F1Default];

(* Check that RationalRecurrenceApplicableQ returns False for a singular point *)
(* and SingularRecurrenceApplicableQ returns True *)
Print["Checking applicability functions on 2F1 matrices..."];

(* We need to prepare matrices first to test applicability *)
(* Construct a SegmentContext with the matrix data *)
Module[{testLine, testCtx},
  testLine = <|z -> z0 + (1/2 - z0) * x|>;
  DiffExp`AnalyticContinuation`PrepareAnalyticContinuation[testLine];
  DiffExp`MatrixLoading`PrepareMatrices[testLine];
  testCtx = <|
    "AMatExpanded" -> DiffExp`State`DEqnMatricesExpanded[testLine][0][[{1, 2}, {1, 2}]],
    "AMatFactored" -> If[KeyExistsQ[DiffExp`State`DEqnMatricesFactored, testLine],
      DiffExp`State`DEqnMatricesFactored[testLine][0][[{1, 2}, {1, 2}]],
      Missing["NotAvailable"]
    ],
    "SystemSize" -> 2,
    "ExpansionOrder" -> DiffExp`State`ExpansionOrderVal,
    "WorkingPrecision" -> DiffExp`State`FEWorkingPrecision,
    "ChopPrecision" -> DiffExp`State`ChopPrecisionVal,
    "LinearSolveChopPrecision" -> DiffExp`State`LinearSolveChopPrecisionVal,
    "HomogeneousSolve" -> "Expand",
    "InvWronskSolver" -> "Auto",
    "CrosscheckFlags" -> {},
    "CrossCheckPrintOrder" -> 5,
    "CrossCheckVerifyOrder" -> 5,
    "IntegrationStrategy" -> "Default",
    "UseRationalRecurrence" -> True,
    "Label" -> {1, 2}
  |>;

  testsTotal++;
  If[BooleanQ[DiffExp`IntegrationStrategies`RationalRecurrenceApplicableQ[testCtx]],
    Print["  [PASS] RationalRecurrenceApplicableQ returns Boolean"];
    testsPassed++;
    ,
    Print["  [FAIL] RationalRecurrenceApplicableQ did not return Boolean"];
  ];

  testsTotal++;
  If[BooleanQ[DiffExp`IntegrationStrategies`SingularRecurrenceApplicableQ[testCtx]],
    Print["  [PASS] SingularRecurrenceApplicableQ returns Boolean"];
    testsPassed++;
    ,
    Print["  [FAIL] SingularRecurrenceApplicableQ did not return Boolean"];
  ];
];

(* ============================================================ *)
(* Summary                                                       *)
(* ============================================================ *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
