(* Focused property/parity coverage for the Euclidean numerical tadpole
   Laurent recurrence.  This test does not invoke FIRE or native transport. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = Prepend[$Path, repoRoot];
Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label,
    If[detail === None, "", ": "], If[detail === None, "", detail]]];

eps = FeynmanTrick`FTeps;
dimension = 4 - 2*eps;
precision = 100;
order = 8;
FeynmanTrick`SetFTOption["WorkingPrecision", precision];
FeynmanTrick`SetFTOption["DimensionExpression", dimension];
FeynmanTrick`SetFTOption["Verbosity", 0];

fast[u_, f_, v_, loops_, n_:order, d_:dimension, wp_:precision] :=
  FeynmanTrick`BoundaryConditions`Private`fastNumericalTadpoleBoundary[
    u, f, v, loops, n, eps, d, wp];
legacy[u_, f_, v_, loops_, n_:order, d_:dimension, wp_:precision] :=
  FeynmanTrick`BoundaryConditions`Private`seriesTadpoleBoundary[
    u, f, v, loops, n, eps, d, wp];

compareCase[label_String, u_, f_, v_Integer, loops_Integer,
    expectedMin_Integer, expectedLeading_] := Module[
  {fastResult, legacyResult, relativeError, scale},
  fastResult = fast[u, f, v, loops];
  legacyResult = legacy[u, f, v, loops];
  scale = Max[1, Max[Abs[N[Last[legacyResult], 80]]]];
  relativeError = Max[Abs[N[Last[fastResult] - Last[legacyResult], 80]]]/scale;
  assert[label <> " Laurent window",
    First[fastResult] === expectedMin &&
      Length[Last[fastResult]] === order - expectedMin + 1,
    {First[fastResult], Length[Last[fastResult]]}];
  assert[label <> " leading coefficient",
    Abs[N[First[Last[fastResult]] - expectedLeading, 80]] < 10^-75,
    First[Last[fastResult]]];
  assert[label <> " Series parity",
    relativeError < 10^-75, relativeError];
];

u = N[13/10, precision];
f = N[7/5, precision];

(* a=v-L d0/2=-3, b=2: the box-bubble deepest-boundary shape. *)
compareCase["negative integer shift", u, f, 1, 2, -1,
  -f^3/(12*u^5)];

(* a=0, b=1. *)
compareCase["zero integer shift", u, f, 2, 1, -1, u^-2];

(* a=2, b=1. *)
compareCase["positive integer shift", u, f, 4, 1, 0, 1/(6*f^2)];

positiveResult = fast[u, f, 4, 1];
coefficientPrecisions = Precision /@
  Select[Last[positiveResult], InexactNumberQ];
assert["fast path restores configured output precision",
  Min[coefficientPrecisions] >= precision, Min[coefficientPrecisions]];

assert["negative F is branch-sensitive and rejected by the fast path",
  fast[u, -f, 1, 2] === $Failed];
assert["complex U is rejected by the fast path",
  fast[u + I/10, f, 1, 2] === $Failed];
assert["non-affine dimension is rejected by the fast path",
  fast[u, f, 1, 2, order, dimension + eps^2] === $Failed];
assert["inexact dimension is rejected by the fast path",
  fast[u, f, 1, 2, order, N[4, precision] - 2*eps] === $Failed];
assert["noninteger Gamma base shift is rejected by the fast path",
  fast[u, f, 1, 1, order, 3 - 2*eps] === $Failed];
assert["nonpositive v is rejected by the fast path",
  fast[u, f, 0, 2] === $Failed];

fallbackExpected = legacy[u, -f, 1, 2];
fallbackActual =
  FeynmanTrick`BoundaryConditions`EvaluateTadpoleBoundary[u, -f, 1, 2, order];
assert["public evaluator preserves the branch-sensitive Series fallback",
  fallbackActual === fallbackExpected];

boxU = N[0.10041327902827258006, 300];
boxF = N[0.00361394742846773332, 300];
fastTiming = First@AbsoluteTiming[
  boxResult = fast[boxU, boxF, 1, 2, 20, dimension, 300]];
assert["order-20 box-bubble-shaped recurrence completes",
  First[boxResult] === -1 && Length[Last[boxResult]] === 22];
Print["  INFO: order-20 WP300 fast seconds = ", fastTiming];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
