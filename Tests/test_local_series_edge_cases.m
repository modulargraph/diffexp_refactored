(* Edge-case test suite for DiffExp`LocalSeries` *)
(* Covers regular singular, infinity, logs, epsilon poles, nonlinear epsilon dependence, and irregular rejection. *)

SetDirectory[DirectoryName[$InputFileName]];
scriptDir = Directory[];
parentDir = ParentDirectory[scriptDir];
$Path = Prepend[$Path, FileNameJoin[{parentDir, "DiffExp"}]];

Print["==========================================="];
Print["LocalSeries Edge-Case Test Suite"];
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

zeroResidualQ[sol_, q_, c_] := Module[{res},
  res = DiffExp`LocalSeries`ValidateFiniteWidthSolution[sol, q, c];
  PossibleZeroQ[res["MaxResidual"]] &&
    AllTrue[Flatten[Values[res["Residuals"]]], PossibleZeroQ]
];

Clear[zz, ee, aa, cc];

Print["\n--- Regular singular and fractional powers ---"];

aa = 7/3;
cc = 5/4;
qFinite = <|{0, 0} -> 1, {1, 0} -> -1|>;
cFinite = <|{0, 0} -> {{aa}}, {1, 0} -> {{cc - aa}}|>;

solFractional = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qFinite,
  cFinite,
  1,
  aa,
  {5, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

pass[
  "fractional exponent sector gives expected coefficients",
  Table[DiffExp`LocalSeries`FiniteWidthCoefficient[solFractional, 0, n, 0][[1]], {n, 0, 5}] ===
    Table[Pochhammer[cc, n]/n!, {n, 0, 5}]
];
pass["fractional exponent residual vanishes", zeroResidualQ[solFractional, qFinite, cFinite]];

Print["\n--- Infinity coordinate ---"];

aa = 11/5;
cc = 9/7;
qInfinity = <|{0, 0} -> 1, {1, 0} -> -1|>;
cInfinity = <|{0, 0} -> {{-aa}}, {1, 0} -> {{aa + cc}}|>;

solInfinity = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qInfinity,
  cInfinity,
  1,
  -aa,
  {5, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

pass[
  "z=1/x infinity sector has expected regular tail",
  Table[DiffExp`LocalSeries`FiniteWidthCoefficient[solInfinity, 0, n, 0][[1]], {n, 0, 5}] ===
    Table[Pochhammer[cc, n]/n!, {n, 0, 5}]
];
pass["infinity residual vanishes", zeroResidualQ[solInfinity, qInfinity, cInfinity]];

Print["\n--- Jordan and resonant logarithms ---"];

nj = {{0, 1}, {0, 0}};
qJordan = <|{0, 0} -> 1|>;
cJordan = <|{0, 0} -> nj|>;

solJordan = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qJordan,
  cJordan,
  2,
  0,
  {0, 0, 0, 1},
  "Prescribed" -> <|{0, 0, 0} -> {0, 1}|>
];

pass[
  "nilpotent residue generates one logarithm",
  DiffExp`LocalSeries`FiniteWidthCoefficient[solJordan, 0, 0, 0] === {0, 1} &&
    DiffExp`LocalSeries`FiniteWidthCoefficient[solJordan, 0, 0, 1] === {1, 0}
];
pass["Jordan residual vanishes", zeroResidualQ[solJordan, qJordan, cJordan]];

qResonant = <|{0, 0} -> 1|>;
cResonant = <|
  {0, 0} -> {{0, 0}, {0, 1}},
  {1, 0} -> {{0, 0}, {1, 0}}
|>;

solResonant = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qResonant,
  cResonant,
  2,
  0,
  {2, 0, 0, 1},
  "Prescribed" -> <|{0, 0, 0} -> {1, 0}|>
];

pass[
  "integer resonance between distinct exponents forces a logarithm",
  DiffExp`LocalSeries`FiniteWidthCoefficient[solResonant, 0, 1, 1] === {0, 1}
];
pass["integer-resonant residual vanishes", zeroResidualQ[solResonant, qResonant, cResonant]];

Print["\n--- epsilon-dependent exponents beyond linear order ---"];

qEps = <|{0, 0} -> 1|>;
cEps = <|{0, 1} -> {{2}}, {0, 2} -> {{3}}|>;
solEps = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qEps,
  cEps,
  1,
  0,
  {0, 0, 3, 3},
  "Prescribed" -> <|
    {0, 0, 0} -> {1},
    {1, 0, 0} -> {0},
    {2, 0, 0} -> {0},
    {3, 0, 0} -> {0}
  |>
];

expectedEpsLog = {
  {1, 0, 0, 0},
  {0, 2, 0, 0},
  {0, 3, 2, 0},
  {0, 0, 6, 4/3}
};

pass[
  "nonlinear epsilon exponent expands into mixed log powers",
  Table[DiffExp`LocalSeries`FiniteWidthCoefficient[solEps, k, 0, ell][[1]], {k, 0, 3}, {ell, 0, 3}] === expectedEpsLog
];
pass["nonlinear epsilon residual vanishes", zeroResidualQ[solEps, qEps, cEps]];

Print["\n--- epsilon poles ---"];

badEps = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  <|{0, 0} -> 1|>,
  <|{0, -1} -> {{1}}|>,
  1,
  0,
  {0, 0, 0, 0},
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

pass["semisimple epsilon pole is rejected as non-Laurent", MatchQ[badEps, _Failure]];

nEps = {{0, 1}, {0, 0}};
solNilEps = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  <|{0, 0} -> 1|>,
  <|{0, -1} -> nEps|>,
  2,
  0,
  {0, -1, 0, 1},
  "Prescribed" -> <|{0, 0, 0} -> {0, 1}|>
];

pass[
  "nilpotent epsilon pole remains finite Laurent-log",
  DiffExp`LocalSeries`FiniteWidthCoefficient[solNilEps, -1, 0, 1] === {1, 0} &&
    DiffExp`LocalSeries`FiniteWidthCoefficient[solNilEps, 0, 0, 0] === {0, 1}
];

Print["\n--- epsilon-dependent finite-width denominator ---"];

qMixed = <|{0, 0} -> 1, {1, -1} -> 1|>;
cMixed = <|{0, 0} -> {{3/2}}, {1, -1} -> {{3/2}}|>;
solMixed = DiffExp`LocalSeries`RecursiveFiniteWidthSolve[
  qMixed,
  cMixed,
  1,
  3/2,
  {4, 0, 1, 0},
  "EpsilonWorkMax" -> 2,
  "EquationEpsilonMin" -> -1,
  "Prescribed" -> <|{0, 0, 0} -> {1}|>
];

pass[
  "finite-width recurrence accepts epsilon-pole denominator with work buffer",
  AssociationQ[solMixed] && zeroResidualQ[solMixed, qMixed, cMixed]
];

Print["\n--- Fuchsianization and irregular rejection ---"];

badBasis = {{0, -1/zz}, {0, 0}};
fb = DiffExp`LocalSeries`FuchsianizeLocal[badBasis, zz, "MaxSteps" -> 20];
pass[
  "apparent higher pole can be removed by lattice saturation",
  ListQ[fb] && DiffExp`LocalSeries`LocalZValuation[fb[[2, 1, 2]], zz] >= 0
];

irregular = DiffExp`LocalSeries`FuchsianizeLocal[{{1/zz}}, zz, "MaxSteps" -> 8];
pass["true irregular scalar pole does not falsely Fuchsianize", MatchQ[irregular, _Failure]];

Print["\n==========================================="];
Print["LocalSeries edge-case tests: ", testsPassed, " / ", testsTotal, " passed"];
If[testsPassed === testsTotal, Print["All tests PASSED!"], Print["Some tests FAILED!"]];
Print["==========================================="];

If[testsPassed === testsTotal, Quit[0], Quit[1]];
