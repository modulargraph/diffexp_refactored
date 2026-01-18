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
  DiffExp`DeltaPrescriptions -> {t - 16 + I * \[Delta]},
  DiffExp`MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_EqualMass_Matrices"}] <> "/",
  DiffExp`Verbosity -> 1,
  DiffExp`UseMobius -> True,
  DiffExp`UsePade -> True
};

Print["Loading configuration..."];
DiffExp`LoadConfiguration[EqualMassConfiguration];

Print["\nCurrent configuration:"];
Print[DiffExp`CurrentConfiguration[]];

(* Boundary conditions *)
Print["\n=== Setting up boundary conditions ==="];
EqualMassBoundaryConditions = {
  "?",
  "?",
  DiffExp`\[Epsilon] (1 + 3 DiffExp`\[Epsilon]) (1 + 4 DiffExp`\[Epsilon]) * (
    -4 E^(3 EulerGamma DiffExp`\[Epsilon]) Gamma[DiffExp`\[Epsilon]]^3/t +
    6 E^(3 EulerGamma DiffExp`\[Epsilon]) (-1/t)^(1 + DiffExp`\[Epsilon]) DiffExp`\[Epsilon] Gamma[-DiffExp`\[Epsilon]]^2 Gamma[DiffExp`\[Epsilon]]^3/Gamma[-2 DiffExp`\[Epsilon]] +
    8 E^(3 EulerGamma DiffExp`\[Epsilon]) (-1/t)^(1 + 2 DiffExp`\[Epsilon]) DiffExp`\[Epsilon] Gamma[-DiffExp`\[Epsilon]]^3 Gamma[DiffExp`\[Epsilon]] Gamma[2 DiffExp`\[Epsilon]]/Gamma[-3 DiffExp`\[Epsilon]] +
    3 E^(3 EulerGamma DiffExp`\[Epsilon]) (-1/t)^(1 + 3 DiffExp`\[Epsilon]) DiffExp`\[Epsilon] Gamma[-DiffExp`\[Epsilon]]^4 Gamma[3 DiffExp`\[Epsilon]]/Gamma[-4 DiffExp`\[Epsilon]]
  ),
  E^(3 EulerGamma DiffExp`\[Epsilon]) DiffExp`\[Epsilon]^3 Gamma[DiffExp`\[Epsilon]]^3
};

Print["Preparing boundary conditions..."];
PreparedBCs = DiffExp`PrepareBoundaryConditions[EqualMassBoundaryConditions, <|t -> -1/x|>];
Print["Boundary conditions prepared\n"];

(* Transport to t = -1 *)
Print["=== Transporting to t = -1 ==="];
Results1 = DiffExp`TransportTo[PreparedBCs, <|t -> -1|>];
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
Results2 = DiffExp`TransportTo[Results1, <|t -> x|>, 5, True];
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
ResultsForEval = DiffExp`ToPiecewise[Results2];

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
(* Check that results are numerical *)
If[And @@ Flatten[NumericQ /@ evaluatedResults],
  Print["  [PASS] Results are numerical"];
  testsPassed++;
  ,
  Print["  [FAIL] Results contain non-numerical values"];
];

Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
