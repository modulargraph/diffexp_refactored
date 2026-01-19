(* Test the refactored DiffExp package with the equal mass banana example *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Testing Refactored DiffExp with Equal Mass Banana"];
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

Print["\nCurrent configuration:"];
Print[CurrentConfiguration[]];

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

(* Transport to t = -1 *)
Print["=== Transporting to t = -1 ==="];
Results1 = TransportTo[PreparedBCs, <|t -> -1|>];
Print["Transport to t = -1 complete"];
Print["Result point: ", Results1[[1]]];

(* Check that result has the expected structure *)
testsPassed = 0;
testsTotal = 0;

testsTotal++;
If[Head[Results1] === List && Length[Results1] >= 2,
  Print["  [PASS] TransportTo returned expected structure"];
  testsPassed++;
  ,
  Print["  [FAIL] TransportTo returned unexpected structure"];
];

testsTotal++;
If[AssociationQ[Results1[[1]]] && KeyExistsQ[Results1[[1]], t],
  Print["  [PASS] Result contains point with key 't'"];
  testsPassed++;
  ,
  Print["  [FAIL] Result missing point with key 't'"];
];

testsTotal++;
If[Results1[[1]][t] === -1,
  Print["  [PASS] Result evaluated at t = -1"];
  testsPassed++;
  ,
  Print["  [FAIL] Result not at expected point t = -1, got: ", Results1[[1]][t]];
];

(* Transport from t = -1 to t = 5, saving expansions *)
Print["\n=== Transporting from t = -1 to t = 5 (with save) ==="];
Results2 = TransportTo[Results1, <|t -> x|>, 5, True];
Print["Transport complete!"];

testsTotal++;
If[Head[Results2] === List && Length[Results2] >= 2,
  Print["  [PASS] Transport with save returned list"];
  testsPassed++;
  ,
  Print["  [FAIL] Transport with save returned unexpected structure"];
];

(* Test ToPiecewise *)
Print["\n=== Testing ToPiecewise ==="];
ResultsForEval = ToPiecewise[Results2];

testsTotal++;
If[Head[ResultsForEval] === List && MatrixQ[ResultsForEval],
  Print["  [PASS] ToPiecewise returned matrix"];
  testsPassed++;
  ,
  Print["  [FAIL] ToPiecewise returned unexpected structure"];
];

(* Evaluate at t = 3 *)
Print["\n=== Evaluating at t = 3 ==="];
evalPoint = 3;
Print["Master integral values at t = ", evalPoint, ":"];
evaluatedResults = Table[
  ResultsForEval[[i, j]][evalPoint] // N[#, 10] &,
  {i, 1, 4}, {j, 1, Min[5, Dimensions[ResultsForEval][[2]]]}
];
Print[TableForm[evaluatedResults]];

testsTotal++;
(* Check that results are numerical - use N to convert arbitrary precision to machine numbers *)
numericCheck = And @@ Flatten[Table[
  NumericQ[evaluatedResults[[i, j]]] || NumericQ[N[evaluatedResults[[i, j]]]],
  {i, Length[evaluatedResults]}, {j, Length[evaluatedResults[[1]]]}
]];
If[numericCheck,
  Print["  [PASS] Results are numerical"];
  testsPassed++;
  ,
  Print["  [FAIL] Results contain non-numerical values"];
];

(* Round-trip test: transport back to t = -1 and compare with original *)
Print["\n=== Round-trip Test: Transporting back to t = -1 ==="];
Results3 = TransportTo[Results2[[1]], <|t -> x|>, -1];
Print["Transport back complete!"];
Print["Result point: ", Results3[[1]]];

testsTotal++;
If[Results3[[1]][t] === -1,
  Print["  [PASS] Round-trip returned to t = -1"];
  testsPassed++;
  ,
  Print["  [FAIL] Round-trip did not return to t = -1, got: ", Results3[[1]][t]];
];

(* Compare values from original transport (Results1) with round-trip (Results3) *)
Print["\n=== Comparing original vs round-trip values ==="];
originalValues = Results1[[2]];
roundTripValues = Results3[[2]];

(* Calculate maximum relative difference *)
maxRelDiff = 0;
Do[
  orig = originalValues[[i, j]];
  roundTrip = roundTripValues[[i, j]];
  If[NumericQ[orig] && NumericQ[roundTrip] && Abs[orig] > 10^-50,
    relDiff = Abs[(orig - roundTrip)/orig];
    If[NumericQ[relDiff] && relDiff > maxRelDiff,
      maxRelDiff = relDiff;
    ];
  ];
  , {i, Length[originalValues]}, {j, Length[originalValues[[1]]]}
];

Print["Maximum relative difference: ", N[maxRelDiff, 5]];

testsTotal++;
(* Allow for some numerical error - should be very small *)
If[maxRelDiff < 10^-10,
  Print["  [PASS] Round-trip values match original (rel. diff < 10^-10)"];
  testsPassed++;
  ,
  Print["  [FAIL] Round-trip values differ significantly from original"];
  Print["  Original values at t=-1:"];
  Print[TableForm[N[originalValues, 10]]];
  Print["  Round-trip values at t=-1:"];
  Print[TableForm[N[roundTripValues, 10]]];
];

Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
