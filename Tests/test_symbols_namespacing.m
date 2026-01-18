(* Test script for DiffExp Symbols namespacing *)
(* This test verifies that symbols are correctly shared across subpackages *)

(* Set the directory to the script location *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
Print["Working directory: ", scriptDir];

(* Load the subpackages in order *)
Print["\n=== Loading Symbols subpackage ==="];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp", "Symbols.m"}]];
Print["Symbols loaded."];

Print["\n=== Loading State subpackage ==="];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp", "State.m"}]];
Print["State loaded."];

Print["\n=== Loading Utilities subpackage ==="];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp", "Utilities.m"}]];
Print["Utilities loaded."];

Print["\n=== Loading SeriesOps subpackage ==="];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp", "SeriesOps.m"}]];
Print["SeriesOps loaded."];

Print["\n=== Loading Integration subpackage ==="];
Get[FileNameJoin[{ParentDirectory[scriptDir], "DiffExp", "Integration.m"}]];
Print["Integration loaded."];

(* Test 1: Check that eps evaluates to \[Epsilon] *)
Print["\n=== Test 1: eps alias ==="];
Print["DiffExp`Symbols`eps = ", DiffExp`Symbols`eps];
testResult1 = (DiffExp`Symbols`eps === DiffExp`Symbols`\[Epsilon]);
Print["Test 1 (eps === \[Epsilon]): ", If[testResult1, "PASSED", "FAILED"]];

(* Test 2: Check that x defaults to Global`x *)
Print["\n=== Test 2: Default line parameter ==="];
Print["DiffExp`Symbols`x = ", DiffExp`Symbols`x];
testResult2 = (DiffExp`Symbols`x === Global`x);
Print["Test 2 (x === Global`x): ", If[testResult2, "PASSED", "FAILED"]];

(* Test 3: Check that DependsQ works with the x symbol *)
Print["\n=== Test 3: DependsQ function ==="];
testExpr = 3 Global`x^2 + 5 Global`x;
testResult3 = DiffExp`Utilities`DependsQ[testExpr, Global`x];
Print["DependsQ[", testExpr, ", x] = ", testResult3];
Print["Test 3 (DependsQ detects x): ", If[testResult3, "PASSED", "FAILED"]];

(* Test 4: Check SExpand works *)
Print["\n=== Test 4: SExpand function ==="];
testSeries = Series[Exp[Global`x], {Global`x, 0, 5}];
expandedSeries = DiffExp`SeriesOps`SExpand[testSeries];
testResult4 = Head[expandedSeries] === SeriesData;
Print["SExpand[Series[Exp[x],{x,0,5}]] = ", expandedSeries];
Print["Test 4 (SExpand returns SeriesData): ", If[testResult4, "PASSED", "FAILED"]];

(* Test 5: Check UpdateIntReps works *)
Print["\n=== Test 5: UpdateIntReps function ==="];
DiffExp`Integration`UpdateIntReps[3];
testResult5 = Length[DiffExp`Integration`IntReps] > 0;
Print["IntReps has ", Length[DiffExp`Integration`IntReps], " rules"];
Print["Test 5 (IntReps populated): ", If[testResult5, "PASSED", "FAILED"]];

(* Test 6: Check DiffExpIntegrate works with simple series *)
Print["\n=== Test 6: DiffExpIntegrate function ==="];
testSeries6 = Series[Global`x^2, {Global`x, 0, 10}];
integratedSeries = DiffExp`Integration`DiffExpIntegrate[testSeries6];
(* The integral of x^2 should be x^3/3 *)
testResult6a = Head[integratedSeries] === SeriesData;
(* Check the coefficient of x^3 is 1/3 *)
testResult6b = (SeriesCoefficient[integratedSeries, 3] === 1/3);
Print["DiffExpIntegrate[x^2 series] = ", integratedSeries];
Print["Test 6a (returns SeriesData): ", If[testResult6a, "PASSED", "FAILED"]];
Print["Test 6b (x^3 coeff is 1/3): ", If[testResult6b, "PASSED", "FAILED"]];

(* Test 7: Check integration with Logx *)
Print["\n=== Test 7: Integration with Logx ==="];
(* Integrate x^2 * Log[x] should give x^3(Log[x]/3 - 1/9) *)
testSeries7 = Series[Global`x^2 Log[Global`x], {Global`x, 0, 10}];
integratedSeries7 = DiffExp`Integration`DiffExpIntegrate[testSeries7];
testResult7 = Head[integratedSeries7] === SeriesData;
Print["DiffExpIntegrate[x^2 Log[x] series] = ", integratedSeries7];
Print["Test 7 (integration with Log): ", If[testResult7, "PASSED", "FAILED"]];

(* Summary *)
Print["\n=== Summary ==="];
allPassed = testResult1 && testResult2 && testResult3 && testResult4 && testResult5 && testResult6a && testResult6b && testResult7;
Print["All tests: ", If[allPassed, "PASSED", "SOME FAILED"]];
