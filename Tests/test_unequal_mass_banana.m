(* Test the refactored DiffExp package with the unequal mass banana example *)
(* This is a simpler test that verifies configuration loading *)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Testing Refactored DiffExp with Unequal Mass Banana"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading refactored DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

(* Test configuration for unequal mass banana *)
(* Use DiffExp`State` symbols for configuration keys to match the package's internal symbols *)
Print["=== Unequal Mass Banana Configuration ==="];
UnequalMassConfiguration = {
  DiffExp`State`ChopPrecision -> 500,
  DiffExp`State`DivisionOrder -> 4,
  DiffExp`State`ExpansionOrder -> 70,
  DiffExp`State`MatrixDirectory -> FileNameJoin[{scriptDir, "Banana_Matrices"}] <> "/",
  DiffExp`State`RadiusOfConvergence -> 10,
  DiffExp`State`UseMobius -> True,
  DiffExp`State`UsePade -> True,
  WorkingPrecision -> 1000,
  DiffExp`State`Verbosity -> 1
};

Print["Loading configuration..."];
LoadConfiguration[UnequalMassConfiguration];

Print["\nCurrent configuration:"];
Print[CurrentConfiguration[]];

(* Test Results *)
testsPassed = 0;
testsTotal = 0;

(* Check that matrices were loaded *)
testsTotal++;
If[DiffExp`State`NumIntegrals > 0,
  Print["  [PASS] Matrices loaded, NumIntegrals = ", DiffExp`State`NumIntegrals];
  testsPassed++;
  ,
  Print["  [FAIL] Matrices not loaded properly"];
];

(* Check configuration values *)
testsTotal++;
If[DiffExp`State`ChopPrecisionVal === 500,
  Print["  [PASS] ChopPrecision set correctly"];
  testsPassed++;
  ,
  Print["  [FAIL] ChopPrecision not set correctly, got: ", DiffExp`State`ChopPrecisionVal];
];

testsTotal++;
If[DiffExp`State`ExpansionOrderVal === 70,
  Print["  [PASS] ExpansionOrder set correctly"];
  testsPassed++;
  ,
  Print["  [FAIL] ExpansionOrder not set correctly, got: ", DiffExp`State`ExpansionOrderVal];
];

testsTotal++;
If[DiffExp`State`UseMobiusVal === True,
  Print["  [PASS] UseMobius set correctly"];
  testsPassed++;
  ,
  Print["  [FAIL] UseMobius not set correctly"];
];

(* Check that external scales include the kinematic variables *)
testsTotal++;
expectedVars = {psq, mm1, mm2, mm3, mm4};
If[ContainsAll[DiffExp`State`ExternalScalesVal, expectedVars],
  Print["  [PASS] ExternalScales contains expected kinematic variables"];
  testsPassed++;
  ,
  Print["  [FAIL] ExternalScales missing variables, got: ", DiffExp`State`ExternalScalesVal];
];

Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
