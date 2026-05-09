(* Test script for DiffExp`LocalSeries` *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["LocalSeries Recursive Solver Tests"];
Print["===========================================\n"];

Get["DiffExp.m"];

testsPassed = 0;
testsTotal = 0;

pass[name_, condition_] := Module[{},
  testsTotal++;
  If[TrueQ[condition],
    testsPassed++;
    Print["  [PASS] ", name],
    Print["  [FAIL] ", name]
  ];
];

minMatrixValuation[m_, var_] := Module[{vals},
  vals = DeleteCases[DiffExp`LocalSeries`LocalZValuation[#, var] & /@ Flatten[m], Infinity];
  If[vals === {}, Infinity, Min[vals]]
];

Clear[zz, ee, aa, cc];

Print["\n--- Rational Laurent coefficients ---"];

lcZ = DiffExp`LocalSeries`LaurentCoefficientsRational[
  zz^-2/(1 - zz),
  zz,
  {-3, 3}
];
pass[
  "z-Laurent recurrence extracts negative and non-negative powers",
  (Lookup[lcZ, #] & /@ Range[-3, 3]) === {0, 1, 1, 1, 1, 1, 1}
];

lcEps = DiffExp`LocalSeries`LaurentCoefficientsRational[
  1/(ee*(1 - ee)),
  ee,
  {-2, 3}
];
pass[
  "epsilon Laurent recurrence handles a pole",
  (Lookup[lcEps, #] & /@ Range[-2, 3]) === {0, 1, 1, 1, 1, 1}
];

Print["\n--- Fuchsian lattice saturation ---"];

Mbad = {{0, -1/zz}, {0, 0}};
fb = DiffExp`LocalSeries`FuchsianizeLocal[Mbad, zz, "MaxSteps" -> 20];

pass[
  "apparent higher theta-pole is Fuchsianized",
  ListQ[fb] && Length[fb] === 2 && minMatrixValuation[fb[[2]], zz] >= 0
];

pass[
  "returned gauge reproduces the transformed connection",
  ListQ[fb] && DiffExp`LocalSeries`LocalConnectionMatrix[Mbad, fb[[1]], zz] === fb[[2]]
];

Print["\n--- Denominator clearing ---"];

aa = 23/10;
cc = 17/10;
prep = DiffExp`LocalSeries`PrepareFiniteWidthData[
  {{aa + cc*zz/(1 - zz)}},
  zz,
  ee,
  {0, 0}
];

pass[
  "regular denominator is normalized to q(0)=1",
  prep["q"] === 1 - zz
];

pass[
  "cleared numerator has finite z-width",
  prep["CAssoc"][{0, 0}] === {{aa}} && prep["CAssoc"][{1, 0}] === {{cc - aa}}
];

Print["\n--- Scalar finite-width recurrence ---"];

qAssoc = <|
  {0, 0} -> 1,
  {1, 0} -> -1
|>;

CAssoc = <|
  {0, 0} -> {{aa}},
  {1, 0} -> {{cc - aa}}
|>;

sol = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qAssoc,
  CAssoc,
  1,
  aa,
  {6, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

gotScalar = Table[
  DiffExp`LocalSeries`FiniteWidthCoefficient[sol, 0, n, 0][[1]],
  {n, 0, 6}
];
expectedScalar = Table[Pochhammer[cc, n]/n!, {n, 0, 6}];

pass[
  "scalar recurrence matches z^a (1-z)^(-c)",
  Simplify[gotScalar == expectedScalar]
];

resScalar = DiffExp`LocalSeries`ValidateFiniteWidthSolution[sol, qAssoc, CAssoc];
pass[
  "scalar recurrence residual vanishes",
  PossibleZeroQ[resScalar["MaxResidual"]] && AllTrue[Flatten[Values[resScalar["Residuals"]]], PossibleZeroQ]
];

localRun = DiffExp`LocalSeries`SolveLocalFuchsianSeries[
  {{aa + cc*zz/(1 - zz)}},
  zz,
  ee,
  aa,
  {4, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>,
  "ReturnOriginalSolution" -> False
];

pass[
  "end-to-end local wrapper reproduces the finite-width scalar solution",
  AssociationQ[localRun] &&
    Table[
      DiffExp`LocalSeries`FiniteWidthCoefficient[localRun["RegularSolution"], 0, n, 0][[1]],
      {n, 0, 4}
    ] === expectedScalar[[1 ;; 5]]
];

Print["\n--- epsilon-dependent exponents and logs ---"];

qLog = <|{0, 0} -> 1|>;
LogC = <|{0, 1} -> {{2}}|>;

solLog = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qLog,
  LogC,
  1,
  0,
  {0, 0, 3, 3},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

gotLog = Table[
  DiffExp`LocalSeries`FiniteWidthCoefficient[solLog, k, 0, ell][[1]],
  {k, 0, 3}, {ell, 0, 3}
];
expectedLog = Table[
  If[ell == k, 2^k/k!, 0],
  {k, 0, 3}, {ell, 0, 3}
];

pass[
  "theta y = 2 eps y generates the expected log tower",
  gotLog === expectedLog
];

resLog = DiffExp`LocalSeries`ValidateFiniteWidthSolution[solLog, qLog, LogC];
pass[
  "log tower residual vanishes",
  PossibleZeroQ[resLog["MaxResidual"]] && AllTrue[Flatten[Values[resLog["Residuals"]]], PossibleZeroQ]
];

Print["\n--- negative epsilon powers ---"];

bad = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  <|{0, 0} -> 1|>,
  <|{0, -1} -> {{1}}|>,
  1,
  0,
  {0, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

pass[
  "non-Laurent z^(1/eps) behavior is rejected",
  MatchQ[bad, _Failure]
];

nil = {{0, 1}, {0, 0}};
solNil = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  <|{0, 0} -> 1|>,
  <|{0, -1} -> nil|>,
  2,
  0,
  {0, -1, 0, 1},
  "Prescribed" -> <|{0, 0, 0} -> {0, 1}|>
];

pass[
  "nilpotent epsilon pole gives a finite Laurent-log tail",
  DiffExp`LocalSeries`FiniteWidthCoefficient[solNil, -1, 0, 1] === {1, 0} &&
    DiffExp`LocalSeries`FiniteWidthCoefficient[solNil, 0, 0, 0] === {0, 1}
];

Print["\n--- Gauge multiplication ---"];

TAssoc = DiffExp`LocalSeries`RationalMatrixZEpsLaurentAssoc[
  {{1/zz}},
  zz,
  {-1, 0},
  ee,
  {0, 0}
];

fSol = DiffExp`LocalSeries`ApplyGaugeToSolution[
  sol,
  TAssoc,
  {-1, 5},
  {0, 0}
];

pass[
  "gauge multiplication shifts z powers correctly",
  DiffExp`LocalSeries`FiniteWidthCoefficient[fSol, 0, -1, 0][[1]] === 1 &&
    DiffExp`LocalSeries`FiniteWidthCoefficient[fSol, 0, 0, 0][[1]] === cc
];

Print["\n--- Implementation hygiene ---"];

source = Import[FileNameJoin[{parentDir, "DiffExp", "LocalSeries.m"}], "Text"];
pass[
  "LocalSeries does not call the built-in asymptotic constructor",
  !StringContainsQ[source, RegularExpression["\\bSeries\\s*\\["]]
];
pass[
  "LocalSeries does not construct built-in asymptotic objects",
  !StringContainsQ[source, RegularExpression["\\bSeriesData\\s*\\["]]
];

Print["\n==========================================="];
Print["LocalSeries tests: ", testsPassed, " / ", testsTotal, " passed"];
Print["==========================================="];

If[testsPassed === testsTotal, Quit[0], Quit[1]];
