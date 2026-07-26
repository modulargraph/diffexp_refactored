(* Focused tests for canonical dlog transport, including algebraic letters.
   This path is intentionally separate from the rational LoadSystem solver. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repoRoot, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[
  Catch[expression, "DiffExp2Error"]];

x = Global`x;

regularSystem = DiffExp2`LoadCanonicalSystem[<|
  "ConstantMatrices" -> {{{1}}},
  "Letters" -> {x},
  "Variables" -> {x}
|>];
assert["canonical_system_schema",
  regularSystem["Schema"] === "DiffExp2.CanonicalSystem/v1" &&
  regularSystem["Dimension"] === 1 &&
  regularSystem["LetterCount"] === 1];

regularResult = catchDE2[DiffExp2`TransportCanonicalLine[
  regularSystem,
  {{1, 0, 0}},
  {x -> 1},
  {x -> 2},
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 2,
  "WorkingPrecision" -> 50,
  "UsePade" -> True,
  "ChartCenters" -> {1/2},
  "ChartBoundaries" -> {0, 1}
]];
regularExpected = {1, Log[2], Log[2]^2/2};
assert["canonical_regular_transport",
  !FailureQ[regularResult] &&
  Max[Abs[N[
    First[regularResult["Value"]] - regularExpected, 40]]] < 10^-18 &&
  regularResult["PadeFallbacks"] === 0];

implicitGeometry = catchDE2[DiffExp2`TransportCanonicalLine[
  regularSystem,
  {{1, 0}},
  {x -> 1},
  {x -> 2},
  "ExpansionOrder" -> 20,
  "EpsilonOrder" -> 1,
  "WorkingPrecision" -> 50,
  "UsePade" -> True
]];
assert["canonical_chart_geometry_is_explicit",
  FailureQ[implicitGeometry] && implicitGeometry["ID"] === "E6"];

algebraicSystem = DiffExp2`LoadCanonicalSystem[<|
  "ConstantMatrices" -> {{{1}}},
  "Letters" -> {1 + Sqrt[x]},
  "Variables" -> {x}
|>];
algebraicResult = catchDE2[DiffExp2`TransportCanonicalLine[
  algebraicSystem,
  {{1, 0, 0}},
  {x -> 1},
  {x -> 4},
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 2,
  "WorkingPrecision" -> 50,
  "UsePade" -> True,
  "ChartCenters" -> {1/2},
  "ChartBoundaries" -> {0, 1}
]];
algebraicLog = Log[3/2];
assert["canonical_algebraic_letter_transport",
  !FailureQ[algebraicResult] &&
  Max[Abs[N[
    First[algebraicResult["Value"]] -
      {1, algebraicLog, algebraicLog^2/2}, 40]]] < 10^-18];

singularSystem = DiffExp2`LoadCanonicalSystem[<|
  "ConstantMatrices" -> {
    {{0, 1}, {0, 0}},
    {{1, 0}, {0, 0}}
  },
  "Letters" -> {x, 1 + x},
  "Variables" -> {x}
|>];
singularResult = catchDE2[DiffExp2`TransportCanonicalLine[
  singularSystem,
  {
    {1, 0, 0},
    {0, 0, 0}
  },
  {x -> 0},
  {x -> 1},
  "ExpansionOrder" -> 30,
  "EpsilonOrder" -> 2,
  "WorkingPrecision" -> 50,
  "UsePade" -> True,
  "ChartCenters" -> {0},
  "ChartBoundaries" -> {0, 1}
]];
assert["canonical_regular_singular_start",
  !FailureQ[singularResult] &&
  TrueQ[singularResult["Charts"][[1, "SingularStart"]]] &&
  Max[Abs[N[
    singularResult["Value"][[1]] - regularExpected, 40]]] < 10^-14 &&
  Max[Abs[singularResult["Value"][[2]]]] === 0];

incompatibleBoundary = catchDE2[DiffExp2`TransportCanonicalLine[
  singularSystem,
  {
    {0, 0},
    {1, 0}
  },
  {x -> 0},
  {x -> 1},
  "ExpansionOrder" -> 20,
  "EpsilonOrder" -> 1,
  "WorkingPrecision" -> 50,
  "UsePade" -> True,
  "ChartCenters" -> {0},
  "ChartBoundaries" -> {0, 1}
]];
assert["canonical_incompatible_singular_boundary_loud",
  FailureQ[incompatibleBoundary] &&
  incompatibleBoundary["ID"] === "E7"];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
