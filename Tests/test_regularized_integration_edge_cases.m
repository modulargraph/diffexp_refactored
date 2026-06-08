(* Focused tests for singular regularized integration edge cases. *)

SetDirectory[DirectoryName[$InputFileName]];
$Path = Prepend[$Path, FileNameJoin[{ParentDirectory[Directory[]], "DiffExp"}]];

Get["DiffExp.m"];

DiffExp`State`DiffExpConfiguration[RationalizationTolerance] = 10^-30;
DiffExp`State`DiffExpConfiguration[ChopPrecision] = 80;
DiffExp`State`DiffExpConfiguration[WorkingPrecision] = 100;

testsTotal = 0;
testsPassed = 0;

pass[name_] := (testsPassed++; Print["  [PASS] ", name]);
fail[name_, got_] := Print["  [FAIL] ", name, ": ", got];
closeQ[a_, b_, tol_:10^-40] := TrueQ[Abs[N[a - b, 80]] < tol];

x = DiffExp`Symbols`x;
logx = DiffExp`Symbols`Logx;

Print["==========================================="];
Print["Regularized Integration Edge Case Tests"];
Print["==========================================="];

testsTotal++;
pureLogSeries = SeriesData[x, 0, {logx}, 0, 1, 1];
pureLogDecomp = DiffExp`SingularityDecomposition`DecomposeSingularity[{pureLogSeries}];
If[Length[pureLogDecomp] == 1,
  pass["pure Logx series is not dropped by decomposition"],
  fail["pure Logx series is not dropped by decomposition", pureLogDecomp]
];

testsTotal++;
directNonzeroInterval =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 1, 0, {1, 0}, {1, 2}
  ];
If[directNonzeroInterval["MinPower"] == 0 &&
    Length[directNonzeroInterval["Coefficients"]] >= 1 &&
    closeQ[directNonzeroInterval["Coefficients"][[1]], Log[2]],
  pass["no analytic-regularization pole away from x=0"],
  fail["no analytic-regularization pole away from x=0", directNonzeroInterval]
];

testsTotal++;
topOrderPole =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    -1, 1, 0, {1}, {0, -1}
  ];
If[topOrderPole["MinPower"] == -1 &&
    Length[topOrderPole["Coefficients"]] >= 2 &&
    closeQ[topOrderPole["Coefficients"][[1]], 1] &&
    closeQ[topOrderPole["Coefficients"][[2]], I Pi],
  pass["a=-1 regularization preserves finite part with one input order"],
  fail["a=-1 regularization preserves finite part with one input order", topOrderPole]
];

testsTotal++;
explicitLogIntegral =
  DiffExp`RegularizedIntegration`IntegrateSingularTermLaurent[
    0, 0, 0, {pureLogSeries}, {0, 1}
  ];
If[explicitLogIntegral["MinPower"] == 0 &&
    Length[explicitLogIntegral["Coefficients"]] >= 1 &&
    closeQ[explicitLogIntegral["Coefficients"][[1]], -1],
  pass["explicit Logx in Taylor coefficient integrates by parts correctly"],
  fail["explicit Logx in Taylor coefficient integrates by parts correctly", explicitLogIntegral]
];

testsTotal++;
negativeSegmentSeries = Table[
  SeriesData[x, 0, {logx^k/k!}, -1, 0, 1],
  {k, 0, 2}
];
negativeSegmentData = {
  {<|x -> x|>, x -> x, {-1, 0}, {-1, 0}, {negativeSegmentSeries}}
};
prefactorPhaseResult =
  DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
    negativeSegmentData,
    {-1, 0},
    <|"PowerAtLower" -> 0, "PowerAtUpper" -> 1,
      "RationalFactor" -> 1, "Variable" -> x|>,
    0
  ][[1]];
If[prefactorPhaseResult["MinPower"] == 0 &&
    Length[prefactorPhaseResult["Coefficients"]] >= 1 &&
    closeQ[prefactorPhaseResult["Coefficients"][[1]], -1],
  pass["upper endpoint prefactor keeps negative-side phase"],
  fail["upper endpoint prefactor keeps negative-side phase", prefactorPhaseResult]
];

Print["==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"],
  Print["Some tests FAILED"];
  Exit[1]
];
Print["==========================================="];
