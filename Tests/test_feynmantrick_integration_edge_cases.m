(* Focused tests for FeynmanTrick integration edge cases. *)

SetDirectory[DirectoryName[$InputFileName]];
repoDir = ParentDirectory[Directory[]];
SetDirectory[repoDir];

Get["FeynmanTrick/FeynmanTrick.m"];
$Path = Prepend[$Path, FileNameJoin[{repoDir, "DiffExp"}]];
Get["DiffExp.m"];

FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["DimensionExpression", 2 - 2*FeynmanTrick`FTeps];

DiffExp`State`DiffExpConfiguration[RationalizationTolerance] = 10^-30;
DiffExp`State`DiffExpConfiguration[ChopPrecision] = 80;
DiffExp`State`DiffExpConfiguration[WorkingPrecision] = 100;

testsTotal = 0;
testsPassed = 0;

pass[name_] := (testsPassed++; Print["  [PASS] ", name]);
fail[name_, got_] := Print["  [FAIL] ", name, ": ", got];
closeQ[a_, b_, tol_:10^-35] := TrueQ[Abs[N[a - b, 80]] < tol];

Print["==========================================="];
Print["FeynmanTrick Integration Edge Case Tests"];
Print["==========================================="];

Block[{DiffExp`Symbols`x},
  xLocal = DiffExp`Symbols`x;
  logx = DiffExp`Symbols`Logx;
  thetaP = DiffExp`Symbols`\[Theta]p;
  thetaM = DiffExp`Symbols`\[Theta]m;

  testsTotal++;
  param = Global`xxThetaTest;
  straddlingTransport = <|
    "SegmentData" -> {
      {
        <|param -> 1/2 + xLocal|>,
        xLocal -> 1/2 + xLocal,
        {0, 1},
        {-1/2, 1/2},
        {{
          SeriesData[xLocal, 0, {thetaP - thetaM}, 0, 1, 1]
        }}
      }
    },
    "NumIntegrals" -> 1,
    "EpsilonOrder" -> 0,
    "SnapValues" -> {0, 1}
  |>;
  straddlingResult = FeynmanTrick`DiffExpIntegration`IntegrateCombinedMasters[
    straddlingTransport,
    {1},
    1,
    1,
    0,
    {0},
    param
  ];
  If[ListQ[straddlingResult] && Length[straddlingResult] == 1 &&
      closeQ[straddlingResult[[1]], 0],
    pass["theta branches are resolved pointwise on sign-straddling segments"],
    fail["theta branches are resolved pointwise on sign-straddling segments",
      straddlingResult]
  ];

  testsTotal++;
  negativeACData = {
    {
      <|param -> 1 + xLocal|>,
      xLocal -> 1 + xLocal,
      {0, 1/2},
      {-1, -1/2},
      {{
        SeriesData[xLocal, 0, {logx - 2 Pi I thetaM}, 0, 1, 1]
      }}
    }
  };
  negativeACResult = DiffExp`RegularizedIntegration`DefiniteIntegralWithPrefactorLaurent[
    negativeACData,
    {0, 1/2},
    <|"PowerAtLower" -> 0, "PowerAtUpper" -> 0,
      "RationalFactor" -> 1, "Variable" -> param|>,
    0
  ][[1]];
  negativeACExpected = Log[2]/2 - 1/2 - I Pi/2;
  If[negativeACResult["MinPower"] == 0 &&
      Length[negativeACResult["Coefficients"]] >= 1 &&
      closeQ[negativeACResult["Coefficients"][[1]], negativeACExpected],
    pass["negative local Logx uses DiffExp analytic-continuation constants"],
    fail["negative local Logx uses DiffExp analytic-continuation constants",
      negativeACResult]
  ];
];

Print["==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"],
  Print["Some tests FAILED"];
  Exit[1]
];
Print["==========================================="];
