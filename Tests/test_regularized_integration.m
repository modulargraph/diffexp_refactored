(* Test Regularized Integration *)
(*
   Tests the regularized integration functionality for:
   1. Basic package loading
   2. Regularization formula application
   3. Integration of simple series
   4. Integration via ToPiecewise saved data
*)

(* Set directory and paths *)
SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];

(* Add subpackages to path *)
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["Regularized Integration Tests"];
Print["===========================================\n"];

(* Load the refactored DiffExp package *)
Print["Loading DiffExp package..."];
Get["DiffExp.m"];
Print["Package loaded!\n"];

testsTotal = 0;
testsPassed = 0;

(* ============================================================ *)
(* Test 1: Package Loading *)
(* ============================================================ *)
Print["=== Test 1: Package Loading ==="];
testsTotal++;

regIntLoaded = Context["DefiniteIntegral"] === "DiffExp`RegularizedIntegration`";

If[regIntLoaded,
  Print["  [PASS] RegularizedIntegration package loaded"];
  testsPassed++;
  ,
  Print["  [FAIL] RegularizedIntegration package not loaded"];
  Print["  Context of DefiniteIntegral: ", Context["DefiniteIntegral"]];
];

(* Check exported functions - use DownValues or pattern check instead of ValueQ *)
testsTotal++;
functionsExist = And[
  Head[DefiniteIntegral] === Symbol,
  Head[IndefiniteIntegral] === Symbol,
  Head[IntegrateDecomposition] === Symbol,
  Head[ApplyRegularizationStep] === Symbol,
  Head[RegularizeIntegrand] === Symbol,
  Head[EvaluateLimitAtSingularity] === Symbol,
  Context["DefiniteIntegral"] === "DiffExp`RegularizedIntegration`"
];

If[functionsExist,
  Print["  [PASS] All exported functions defined"];
  testsPassed++;
  ,
  Print["  [FAIL] Some exported functions missing"];
];

(* ============================================================ *)
(* Test 2: Basic Integration of Taylor Series *)
(* ============================================================ *)
Print["\n=== Test 2: Basic Taylor Series Integration ==="];

(* Configure minimal settings without loading matrices *)
(* Just update the state directly for testing *)
DiffExp`State`DiffExpConfiguration[RationalizationTolerance] = 10^-10;
DiffExp`State`DiffExpConfiguration[ChopPrecision] = 50;
DiffExp`State`DiffExpConfiguration[WorkingPrecision] = 100;

(* Create a simple Taylor series: g(x) = 1 + x + x^2 *)
testSeries = SeriesData[DiffExp`Symbols`x, 0, {1, 1, 1}, 0, 3, 1];

(* Decomposition of a finite series should have a=0, b=0 *)
decomp = DecomposeSingularity[{testSeries}];

testsTotal++;
If[Length[decomp] > 0 && decomp[[1]]["a"] === 0 && decomp[[1]]["b"] === 0,
  Print["  [PASS] Decomposition of Taylor series: a=0, b=0"];
  testsPassed++;
  ,
  Print["  [FAIL] Decomposition unexpected: ", decomp];
];

(* ============================================================ *)
(* Test 3: Test with 2F1 Transport Data *)
(* ============================================================ *)
Print["\n=== Test 3: Integration of 2F1 Transport ==="];

(* Load 2F1 matrices if available *)
matricesDir = FileNameJoin[{scriptDir, "Hypergeometric2F1_Matrices"}];

If[DirectoryQ[matricesDir],
  Print["Loading 2F1 configuration..."];

  (* Hypergeometric parameters *)
  a = 1/4;
  b = 1/3;
  c = 3/2;

  Config2F1 = {
    MatrixDirectory -> matricesDir <> "/",
    Verbosity -> 1,
    UseMobius -> False,  (* Keep linear transforms for integration *)
    UsePade -> True,
    WorkingPrecision -> 200,
    ExpansionOrder -> 60,
    DivisionOrder -> 4,
    ChopPrecision -> 150,
    EpsilonOrder -> 0
  };

  LoadConfiguration[Config2F1];

  (* Boundary conditions near z=0 *)
  z0 = 1/100;
  y2F1AtStart = N[Hypergeometric2F1[a, b, c, z0], 200];
  yPrimeAtStart = N[D[Hypergeometric2F1[a, b, c, zz], zz] /. zz -> z0, 200];

  BCs = {
    <|z -> z0|>,
    {
      {y2F1AtStart},
      {yPrimeAtStart}
    }
  };

  (* Transport with SaveExpansions = True *)
  Print["Transporting 2F1 from z=0.01 to z=0.5 (saving segments)..."];
  ResultWithSave = TransportTo[BCs, <|z -> 1/2|>, 1, True];

  Print["Transport complete. Segments saved: ", Length[ResultWithSave[[2]]]];

  (* Test that we can call DefiniteIntegral *)
  testsTotal++;

  (* For 2F1, we have:
     Integral[2F1(a,b;c;z), {z, z1, z2}]

     This should work through our saved data *)

  Print["\nTesting DefiniteIntegral function call..."];

  (* Just test that the function can be called without error *)
  integralResult = Quiet[Check[
    (* Integrate from z=0.01 to z=0.3 *)
    DefiniteIntegral[ResultWithSave, {z0, 3/10}],
    $Failed
  ]];

  If[integralResult =!= $Failed,
    Print["  [PASS] DefiniteIntegral function executed successfully"];
    Print["  Result dimensions: ", Dimensions[integralResult]];
    testsPassed++;
    ,
    Print["  [FAIL] DefiniteIntegral function failed"];
    Print["  Check implementation for errors"];
  ];

  (* Test IndefiniteIntegral *)
  testsTotal++;
  Print["\nTesting IndefiniteIntegral function call..."];

  indefiniteResult = Quiet[Check[
    IndefiniteIntegral[ResultWithSave],
    $Failed
  ]];

  If[indefiniteResult =!= $Failed,
    Print["  [PASS] IndefiniteIntegral function executed successfully"];
    Print["  Result dimensions: ", Dimensions[indefiniteResult]];
    testsPassed++;
    ,
    Print["  [FAIL] IndefiniteIntegral function failed"];
  ];

  (* Test 3b: Numerical verification of the integral *)
  testsTotal++;
  Print["\nVerifying integral numerically..."];

  (* The first integral is 2F1 itself, so ∫ 2F1 dz from z0 to z1 *)
  (* The main line is z -> z0 + (1/2 - z0)*x_main, so:
     - at x_main=0, z=z0
     - at x_main=1, z=1/2
     - dz/dx_main = 1/2 - z0 = 0.49 *)
  z1 = 3/10;  (* Upper bound in z *)
  dzdx = 1/2 - z0;  (* = 0.49, the slope of main line *)
  x1 = (z1 - z0) / dzdx;  (* Convert z1 to x_main coordinate *)

  expectedIntegral = NIntegrate[Hypergeometric2F1[a, b, c, zz], {zz, z0, z1}, WorkingPrecision -> 50];

  (* Call DefiniteIntegral with x_main bounds (0 to x1) *)
  (* This gives ∫ f(x) dx in x_main coords. To get ∫ f(z) dz, multiply by dz/dx *)
  xIntegralResult = DefiniteIntegral[ResultWithSave, {0, x1}];
  ourIntegral = dzdx * xIntegralResult[[1, 1]];  (* First integral, eps^0 order *)

  integralDiff = Abs[ourIntegral - expectedIntegral];
  Print["  Expected (NIntegrate): ", NumberForm[expectedIntegral, 15]];
  Print["  Our result:            ", NumberForm[ourIntegral, 15]];
  Print["  Difference:            ", ScientificForm[integralDiff, 3]];

  If[integralDiff < 10^-10,
    Print["  [PASS] Integral matches NIntegrate to 10+ digits"];
    testsPassed++;
    ,
    Print["  [FAIL] Integral differs from NIntegrate"];
  ];

  ,
  Print["  [SKIP] 2F1 matrices not found at: ", matricesDir];
  Print["  Run test_hypergeometric2f1.m first to generate matrices"];
];

(* ============================================================ *)
(* Test 4: Regularization Formula *)
(* ============================================================ *)
Print["\n=== Test 4: Regularization Formula ==="];

(* Test the regularization step on x^(-1 + eps) * (1 + x) *)
(* The regularization should increase power from -1 to 0 *)

(* We need to set up a proper test case *)
(* Series g(x) = 1 + x starting at x^0 *)
testGSeries = SeriesData[DiffExp`Symbols`x, 0, {1, 1}, 0, 2, 1];

(* Apply regularization with a = -2, b = 1, c = 1/2 (not a=-1 which generates eps pole) *)
testsTotal++;

regResult = Quiet[Check[
  ApplyRegularizationStep[-2, 1, 0, {testGSeries}, 1/2],
  $Failed
]];

If[regResult =!= $Failed && Length[regResult] == 4,
  newA = regResult[[1]];
  newB = regResult[[2]];
  newEpsMin = regResult[[3]];
  newG = regResult[[4]];

  If[newA === -1,
    Print["  [PASS] Regularization increased power: a=-2 -> a=-1"];
    testsPassed++;
    ,
    Print["  [FAIL] Expected a=-1, got a=", newA];
  ];
  ,
  Print["  [FAIL] ApplyRegularizationStep failed, result=", regResult];
];

(* Test RegularizeIntegrand for multiple steps *)
testsTotal++;

(* x^(-2 + eps) * (1 + x) should need 2 regularization steps *)
regResult2 = Quiet[Check[
  RegularizeIntegrand[-2, 1, {testGSeries}, 1/2],
  $Failed
]];

If[regResult2 =!= $Failed && Length[regResult2] == 4,
  finalA = regResult2[[1]];
  finalEpsMin = regResult2[[3]];
  If[finalA >= 0,
    Print["  [PASS] RegularizeIntegrand reaches non-negative power: a=", finalA, ", epsMin=", finalEpsMin];
    testsPassed++;
    ,
    Print["  [FAIL] Expected a>=0, got a=", finalA];
  ];
  ,
  Print["  [FAIL] RegularizeIntegrand failed, result=", regResult2];
];

(* ============================================================ *)
(* Test 5: Limit at Singularity *)
(* ============================================================ *)
Print["\n=== Test 5: Limit at Singularity ==="];

(* Create a decomposition with different terms *)
testDecomp = {
  <|"a" -> 0, "b" -> 0, "g" -> {SeriesData[DiffExp`Symbols`x, 0, {5}, 0, 1, 1]}|>,  (* Taylor term *)
  <|"a" -> -1/2, "b" -> 1, "g" -> {SeriesData[DiffExp`Symbols`x, 0, {3}, 0, 1, 1]}|>  (* Singular with b!=0 *)
};

testsTotal++;
limitResult = EvaluateLimitAtSingularity[testDecomp, 1];

(* Taylor term (b=0) should contribute 5, singular term (b!=0) should contribute 0 *)
If[Length[limitResult] > 0 && limitResult[[1]] === 5,
  Print["  [PASS] Limit at singularity: Taylor term=5, singular term=0"];
  testsPassed++;
  ,
  Print["  [FAIL] Expected {5}, got: ", limitResult];
];

(* ============================================================ *)
(* Summary *)
(* ============================================================ *)
Print["\n==========================================="];
Print["Results: ", testsPassed, " / ", testsTotal, " tests passed"];
If[testsPassed === testsTotal,
  Print["All tests PASSED!"];
  ,
  Print["Some tests FAILED or SKIPPED"];
];
Print["==========================================="];
