(* Test package loading for refactored DiffExp *)

(* Set path to subpackages *)
$Path = Prepend[$Path, FileNameJoin[{DirectoryName[$InputFileName], "..", "DiffExp"}]];

Print["Testing refactored DiffExp package loading..."];
Print[""];

(* Test 1: Load subpackages individually *)
Print["Test 1: Loading individual subpackages..."];

testsPassed = 0;
testsTotal = 0;

(* Load Symbols *)
testsTotal++;
Get["Symbols.m"];
If[ValueQ[DiffExp`Symbols`\[Epsilon]::usage],
  Print["  [PASS] Symbols.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Symbols.m failed to load"];
];

(* Load State *)
testsTotal++;
Get["State.m"];
If[ValueQ[DiffExp`State`DiffExpConfiguration],
  Print["  [PASS] State.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] State.m failed to load"];
];

(* Load Utilities *)
testsTotal++;
Get["Utilities.m"];
If[ValueQ[DiffExp`Utilities`PrintInfo::usage],
  Print["  [PASS] Utilities.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Utilities.m failed to load"];
];

(* Load SeriesOps *)
testsTotal++;
Get["SeriesOps.m"];
If[ValueQ[DiffExp`SeriesOps`SExpand::usage],
  Print["  [PASS] SeriesOps.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] SeriesOps.m failed to load"];
];

(* Load Integration *)
testsTotal++;
Get["Integration.m"];
If[ValueQ[DiffExp`Integration`DiffExpIntegrate::usage],
  Print["  [PASS] Integration.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Integration.m failed to load"];
];

(* Load Pade *)
testsTotal++;
Get["Pade.m"];
If[ValueQ[DiffExp`Pade`GetPade::usage],
  Print["  [PASS] Pade.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Pade.m failed to load"];
];

(* Load Mobius *)
testsTotal++;
Get["Mobius.m"];
If[ValueQ[DiffExp`Mobius`GetMobius::usage],
  Print["  [PASS] Mobius.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Mobius.m failed to load"];
];

(* Load AnalyticContinuation *)
testsTotal++;
Get["AnalyticContinuation.m"];
If[ValueQ[DiffExp`AnalyticContinuation`PrepareAnalyticContinuation::usage],
  Print["  [PASS] AnalyticContinuation.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] AnalyticContinuation.m failed to load"];
];

(* Load LineSegmentation *)
testsTotal++;
Get["LineSegmentation.m"];
If[ValueQ[DiffExp`LineSegmentation`RelateLines::usage],
  Print["  [PASS] LineSegmentation.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] LineSegmentation.m failed to load"];
];

(* Load Frobenius *)
testsTotal++;
Get["Frobenius.m"];
If[ValueQ[DiffExp`Frobenius`FrobeniusSolutions::usage],
  Print["  [PASS] Frobenius.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Frobenius.m failed to load"];
];

(* Load LocalSeries *)
testsTotal++;
Get["LocalSeries.m"];
If[ValueQ[DiffExp`LocalSeries`RecursiveFiniteWidthSolve::usage],
  Print["  [PASS] LocalSeries.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] LocalSeries.m failed to load"];
];

(* Load Wronskian *)
testsTotal++;
Get["Wronskian.m"];
If[ValueQ[DiffExp`Wronskian`MatrixLogxInverse::usage],
  Print["  [PASS] Wronskian.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Wronskian.m failed to load"];
];

(* Load MatrixLoading *)
testsTotal++;
Get["MatrixLoading.m"];
If[ValueQ[DiffExp`MatrixLoading`LoadMatrices::usage],
  Print["  [PASS] MatrixLoading.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] MatrixLoading.m failed to load"];
];

(* Load IntegrationStrategies *)
testsTotal++;
Get["IntegrationStrategies.m"];
If[ValueQ[DiffExp`IntegrationStrategies`DispatchStrategy::usage],
  Print["  [PASS] IntegrationStrategies.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] IntegrationStrategies.m failed to load"];
];

(* Load Transport *)
testsTotal++;
Get["Transport.m"];
If[ValueQ[DiffExp`Transport`TransportTo::usage],
  Print["  [PASS] Transport.m loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] Transport.m failed to load"];
];

Print[""];
Print["Test 2: Verifying key symbols..."];

(* Test that x defaults to Global`x *)
testsTotal++;
If[DiffExp`Symbols`x === Global`x,
  Print["  [PASS] x defaults to Global`x"];
  testsPassed++;
  ,
  Print["  [FAIL] x does not default to Global`x, got: ", DiffExp`Symbols`x];
];

(* Test that eps is alias for epsilon *)
testsTotal++;
If[DiffExp`Symbols`eps === DiffExp`Symbols`\[Epsilon],
  Print["  [PASS] eps is alias for \[Epsilon]"];
  testsPassed++;
  ,
  Print["  [FAIL] eps is not alias for \[Epsilon]"];
];

(* Test default configuration *)
testsTotal++;
If[KeyExistsQ[DiffExp`State`DiffExpConfiguration, DiffExp`State`ExpansionOrder],
  Print["  [PASS] Default configuration has ExpansionOrder"];
  testsPassed++;
  ,
  Print["  [FAIL] Default configuration missing ExpansionOrder"];
];

Print[""];
Print["Test 3: Testing basic function availability..."];

(* Test SExpand *)
testsTotal++;
testSeries = Series[1/(1-Global`x), {Global`x, 0, 3}];
result = DiffExp`SeriesOps`SExpand[testSeries];
If[Head[result] === SeriesData,
  Print["  [PASS] SExpand works on series"];
  testsPassed++;
  ,
  Print["  [FAIL] SExpand failed"];
];

(* Test DiffExpIntegrate *)
testsTotal++;
DiffExp`Integration`UpdateIntReps[1];
intResult = DiffExp`Integration`DiffExpIntegrate[testSeries];
If[Head[intResult] === SeriesData,
  Print["  [PASS] DiffExpIntegrate works"];
  testsPassed++;
  ,
  Print["  [FAIL] DiffExpIntegrate failed"];
];

Print[""];
Print["==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED!"];
];
Print["==========================================="];
