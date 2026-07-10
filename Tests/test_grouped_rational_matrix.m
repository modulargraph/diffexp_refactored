(* Focused tests for DiffExp2/Solve.m grouped epsilon-rational Nhat products. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

eps = Global`eps; rho = Global`rho; x = Global`x; t = Global`t;
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];
csTest = <|"Center" -> "grouped-test"|>;

legacySp[m_, fb_Integer, W_Integer] := Module[{frames, idxs},
  frames = Map[DiffExp2`Solve`Private`ratEpsList[#, eps, fb, W] &, m, {2}];
  idxs = Select[Range[W], Module[{i = #},
    AnyTrue[Flatten[frames, 1], #[[i]] =!= 0 &]] &];
  Table[{i + fb - 1, Map[#[[i]] &, frames, {2}]}, {i, idxs}]];

legacyApply[sp_, u_, fb_Integer, W_Integer] := Module[
  {acc = ConstantArray[0, Dimensions[u]]},
  Do[acc += term[[2]] . DiffExp2`Solve`Private`shiftFrameBlock[
      u, term[[1]], fb, W, csTest],
    {term, sp}];
  acc];

legacyValidity[sp_, inputValid_List, frameTop_Integer] := Module[
  {acc = ConstantArray[Infinity, Length[inputValid]], d = Length[inputValid]},
  Do[Do[If[term[[2, r, c]] =!= 0,
      acc[[r]] = Min[acc[[r]],
        DiffExp2`Solve`Private`validShift[inputValid[[c]], term[[1]], frameTop]]],
    {r, d}, {c, d}], {term, sp}];
  acc];

(* Mixed Laurent-polynomial, affine-denominator, quadratic-denominator, and
   symbolic analytic-regulator coefficients.  Scalar/sign variants of the
   affine denominator must canonicalize into one group. *)
fb = -3; W = 13; top = fb + W - 1;
mixed = {
  {eps^-2 (1 + eps + 2 eps^2)/(1 + 5 eps),
    (2 - eps)/(1 - eps + eps^2), 1 + eps},
  {eps (1 + rho eps)/(2 + 10 eps), 0, eps^-1 (3 + eps)},
  {0, eps^2 (1 + 3 eps)/(1 - eps + eps^2),
    -3 (1 + rho eps)/(-3 - 15 eps)}
};
hybrid = catchDE2[DiffExp2`Solve`Private`prepareNhatHybrid[
  {mixed}, eps, fb, W, csTest]];
SeedRandom[20260710];
u = Join[ConstantArray[0, {3, 3}], RandomInteger[{-3, 3}, {3, W - 3}], 2];
oldSp = legacySp[mixed, fb, W];
oldValue = legacyApply[oldSp, u, fb, W];
newValue = DiffExp2`Solve`Private`applyPreparedNhat[
  hybrid["PolynomialSp"][[1]], hybrid["RationalGroups"][[1]],
  u, fb, W, csTest];
assert["grouped_mixed_exact_frame_equivalence",
  AllTrue[Flatten[oldValue - newValue], Together[#] === 0 &]];

qGroups = hybrid["RationalGroups"][[1]];
assert["grouped_normalized_denominators_and_finite_boundary",
  Length[qGroups] === 2 &&
  Sort[#["EntryCount"] & /@ qGroups] === {2, 3} &&
  Max[Flatten[#["NumeratorSp"][[All, 1]] & /@ qGroups]] > top &&
  Max[Length[#["NumeratorSp"]] & /@ qGroups] <
    hybrid["Stats"][[1, "LegacyRationalShiftUpperBound"]]];

inputValid = {4, 7, Infinity};
oldValid = legacyValidity[oldSp, inputValid, top];
newValid = DiffExp2`Solve`Private`updateNhatValidity[
  ConstantArray[Infinity, 3], hybrid["Valuations"][[1]], inputValid, top];
assert["grouped_complete_max_shadow_equivalence",
  oldValid === newValid];

(* The original numerator may extend beyond the frame top.  The grouped
   object must represent Q times the STORED rational series, not silently
   restore those discarded numerator coefficients. *)
fbTiny = -1; WTiny = 2; topTiny = 0;
tiny = {{(1 + 2 eps + 3 eps^2 + 4 eps^3)/(1 + eps)}};
tinyHybrid = DiffExp2`Solve`Private`prepareNhatHybrid[
  {tiny}, eps, fbTiny, WTiny, csTest];
tinyU = {{5, 7}};
tinyOld = legacyApply[legacySp[tiny, fbTiny, WTiny], tinyU, fbTiny, WTiny];
tinyNew = DiffExp2`Solve`Private`applyPreparedNhat[
  tinyHybrid["PolynomialSp"][[1]], tinyHybrid["RationalGroups"][[1]],
  tinyU, fbTiny, WTiny, csTest];
assert["grouped_numerator_beyond_frame_top",
  tinyOld === tinyNew &&
  Max[tinyHybrid["RationalGroups"][[1, 1, "NumeratorSp", All, 1]]] > topTiny];

(* A mixed exact/inexact grouped recurrence can spoil the exact finite-top
   cancellation.  Large denominator/impulse data must therefore route
   through ratEpsList's legacy dense representation as one whole group. *)
fbHuge = -1; WHuge = 5;
huge = {{1/(1 + 10^1200 eps)}};
hugeHybrid = DiffExp2`Solve`Private`prepareNhatHybrid[
  {huge}, eps, fbHuge, WHuge, csTest];
hugeU = {{1, 0, 0, 0, 0}};  (* eps^-1 *)
hugeOldSp = legacySp[huge, fbHuge, WHuge];
hugeOld = legacyApply[hugeOldSp, hugeU, fbHuge, WHuge];
hugeNew = DiffExp2`Solve`Private`applyPreparedNhat[
  hugeHybrid["PolynomialSp"][[1]], hugeHybrid["RationalGroups"][[1]],
  hugeU, fbHuge, WHuge, csTest];
hugeStats = hugeHybrid["Stats"][[1]];
hugeStoredCoeffs = Flatten[#[[2]] & /@ hugeHybrid["PolynomialSp"][[1]]];
assert["grouped_large_exact_denominator_routes_legacy",
  hugeNew === hugeOld && Last[First[hugeNew]] === 0 &&
  hugeHybrid["RationalGroups"][[1]] === {} &&
  hugeStats["GroupedRationalEntries"] === 0 &&
  hugeStats["LegacyRationalEntries"] === 1 &&
  hugeStats["LegacyRationalGroups"] === 1 &&
  Lookup[hugeStats["LegacyRouteReasons"], "DenominatorByteCount", 0] === 1 &&
  Max[ByteCount /@ hugeStoredCoeffs] < 1000];

(* Pure Q=1 input stays entirely on the old sparse path. *)
poly = {{eps^-2 (1 + eps + eps^3), 2 - eps}, {0, eps^2}};
polyHybrid = DiffExp2`Solve`Private`prepareNhatHybrid[
  {poly}, eps, fb, W, csTest];
polySp = legacySp[poly, fb, W];
assert["grouped_q1_uses_legacy_sparse_path",
  polyHybrid["RationalGroups"] === {{}} &&
  polyHybrid["PolynomialSp"][[1]] === polySp &&
  polyHybrid["Stats"][[1, "RationalEntries"]] === 0];

(* Property sweep: exact equality for common and distinct denominators,
   valuations on both sides of zero, and quadratic coefficient recurrences. *)
randomParity = And @@ Table[
  SeedRandom[9000 + seed];
  Module[{fbR = -4, WR = 16, denoms, mr, hr, ur, oldr, newr},
    denoms = {1, 1 + 2 eps, 1 - eps + eps^2};
    mr = Table[Module[{v = RandomInteger[{-2, 2}], q, p},
      q = RandomChoice[denoms];
      p = RandomChoice[{-3, -2, -1, 1, 2, 3}] +
        Sum[RandomInteger[{-3, 3}] eps^k, {k, 1, RandomInteger[{0, 3}]}];
      Cancel[eps^v p/q]], {3}, {3}];
    hr = DiffExp2`Solve`Private`prepareNhatHybrid[
      {mr}, eps, fbR, WR, csTest];
    ur = Join[ConstantArray[0, {3, 3}],
      RandomInteger[{-2, 2}, {3, WR - 3}], 2];
    oldr = legacyApply[legacySp[mr, fbR, WR], ur, fbR, WR];
    newr = DiffExp2`Solve`Private`applyPreparedNhat[
      hr["PolynomialSp"][[1]], hr["RationalGroups"][[1]],
      ur, fbR, WR, csTest];
    oldr === newr],
  {seed, 1, 30}];
assert["grouped_random_exact_property_sweep_30",
  randomParity];

(* Both implementations retain the same loud lower-frame refusal. *)
fbBad = -2; WBad = 8;
badM = {{eps^-1/(1 + 3 eps)}};
badH = DiffExp2`Solve`Private`prepareNhatHybrid[
  {badM}, eps, fbBad, WBad, csTest];
badU = {ConstantArray[1, WBad]};
badOld = catchDE2[legacyApply[legacySp[badM, fbBad, WBad], badU, fbBad, WBad]];
badNew = catchDE2[DiffExp2`Solve`Private`applyPreparedNhat[
  badH["PolynomialSp"][[1]], badH["RationalGroups"][[1]],
  badU, fbBad, WBad, csTest]];
assert["grouped_lower_frame_guard_parity",
  FailureQ[badOld] && badOld["ID"] === "E4" &&
  FailureQ[badNew] && badNew["ID"] === "E4"];

(* Cross-lag fusion: two B_j operators with the same normalized Q carry
   distinct finite-top boundary corrections.  Summing their B_j.U_j right
   sides before ONE causal division must be exactly identical to the old
   sum of two separately divided results. *)
fbFuse = -2; WFuse = 5; topFuse = fbFuse + WFuse - 1;
lagMatrices = {
  {{(1 + 2 eps + 3 eps^2 + 4 eps^3 + 5 eps^4)/(1 + eps)}},
  {{eps^-1 (2 - rho eps + 3 eps^2 + 7 eps^4)/(1 + eps)}}
};
lagHybrid = DiffExp2`Solve`Private`prepareNhatHybrid[
  lagMatrices, eps, fbFuse, WFuse, csTest];
lagGroups = lagHybrid["RationalGroups"];
lagU = {{{0, 0, 2, -1, 3}}, {{0, 0, -4, 5, 1}}};
lagOld = Total[MapThread[
  DiffExp2`Solve`Private`applyRationalMatrixGroups[
    #1, #2, fbFuse, WFuse, csTest] &, {lagGroups, lagU}]];
lagRHS = Total[MapThread[
  DiffExp2`Solve`Private`rationalMatrixGroupNumerator[
    First[#1], #2, fbFuse, WFuse, csTest] &, {lagGroups, lagU}]];
lagFused = DiffExp2`Solve`Private`divideRationalMatrixRHS[
  lagRHS, lagGroups[[1, 1, "DenominatorCoefficients"]], WFuse];
assert["grouped_cross_lag_same_q_finite_top_exact",
  AllTrue[Flatten[lagOld - lagFused], Together[#] === 0 &] &&
  AllTrue[lagGroups,
    Max[#[[1, "NumeratorSp", All, 1]]] > topFuse &]];

(* Per-contribution lower guards are intentionally stronger than checking
   the final fused sum.  These two contributions cancel algebraically, but
   each would discard a nonzero eps^-1 input and must remain loud. *)
cancelMatrices = {{{eps^-1/(1 + eps)}}, {{-eps^-1/(1 + eps)}}};
cancelHybrid = DiffExp2`Solve`Private`prepareNhatHybrid[
  cancelMatrices, eps, -1, 5, csTest];
cancelGroups = cancelHybrid["RationalGroups"];
cancelU = {{1, 2, 3, 4, 5}};
cancelFused = catchDE2[Module[{rr = ConstantArray[0, {1, 5}]},
  Do[rr += DiffExp2`Solve`Private`rationalMatrixGroupNumerator[
      First[gg], cancelU, -1, 5, csTest], {gg, cancelGroups}];
  DiffExp2`Solve`Private`divideRationalMatrixRHS[
    rr, cancelGroups[[1, 1, "DenominatorCoefficients"]], 5]]];
assert["grouped_cross_lag_cancellation_keeps_each_underflow_witness",
  FailureQ[cancelFused] && cancelFused["ID"] === "E4"];

(* End-to-end production path.  The first scalar system has one Q shared by
   Taylor lags 1 and 2; the second has distinct Qs.  Toggle the private seam
   WITHOUT clearing the homogeneous cache to also pin its memo-key mode. *)
fusionChart[name_] := <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 1, "LocalRadius" -> 1, "Singular" -> True, "Name" -> name|>;
fusionReq = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TOrder" -> 6|>;
fusionSystemSolve[A_, name_] := Module[{cs, fused, unfused, prep, sym},
  DiffExp2`Solve`ClearSolveCaches[];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[
    <|"Matrix" -> A, "Variable" -> x|>, fusionChart[name]]];
  If[FailureQ[cs], Return[{cs, cs, cs}]];
  sym = DiffExp2`Solve`Private`clearedSymbolic[cs];
  prep = catchDE2[DiffExp2`Solve`Private`prepareCleared[cs, -6, 15, sym]];
  fused = Block[{
      DiffExp2`Solve`Private`$disableRationalDenominatorFusion = False},
    catchDE2[DiffExp2`Solve`SolveChart[cs, fusionReq]]];
  (* no cache clear: fusion mode must be part of the memo key *)
  unfused = Block[{
      DiffExp2`Solve`Private`$disableRationalDenominatorFusion = True},
    catchDE2[DiffExp2`Solve`SolveChart[cs, fusionReq]]];
  {fused, unfused, prep}];

sameQRun = fusionSystemSolve[
  {{eps/x + (1 + x)/(1 + 2 eps)}}, "grouped_fusion_same_q"];
assert["grouped_cross_lag_same_q_end_to_end_and_cache_key",
  !AnyTrue[sameQRun, FailureQ] &&
  sameQRun[[1, "Basis", "Columns"]] ===
    sameQRun[[2, "Basis", "Columns"]] &&
  Length[sameQRun[[3, "NhatRationalDenominators"]]] === 1 &&
  Count[Flatten[sameQRun[[3, "NhatRationalGroups"]], 1],
    g_Association /; KeyExistsQ[g, "DenominatorIndex"]] >= 2];

distinctQRun = fusionSystemSolve[
  {{eps/x + 1/(1 + 2 eps) + x/(1 - eps + eps^2)}},
  "grouped_fusion_distinct_q"];
assert["grouped_cross_lag_distinct_q_stays_separate",
  !AnyTrue[distinctQRun, FailureQ] &&
  distinctQRun[[1, "Basis", "Columns"]] ===
    distinctQRun[[2, "Basis", "Columns"]] &&
  Length[distinctQRun[[3, "NhatRationalDenominators"]]] === 2];

(* Inhomogeneous recurrence shares the same grouped operator path.  A
   nonresonant fractional-power source exercises source validity, the
   particular solve, and homogeneous compensation-target cache behavior. *)
sourceFusion = <|"Sectors" -> {<|"a" -> 1/2, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[{If[k === 0 && n === 0, 1, 0]},
      {k, 0, 2}, {n, 0, 6}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 6|>|>;
sourceCS = catchDE2[DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{eps/x + (1 + x)/(1 + 2 eps)}}, "Variable" -> x|>,
  fusionChart["grouped_fusion_source"]]];
DiffExp2`Solve`ClearSolveCaches[];
sourceFused = Block[{
    DiffExp2`Solve`Private`$disableRationalDenominatorFusion = False},
  catchDE2[DiffExp2`Solve`SolveChart[sourceCS, fusionReq, sourceFusion]]];
DiffExp2`Solve`ClearSolveCaches[];
sourceUnfused = Block[{
    DiffExp2`Solve`Private`$disableRationalDenominatorFusion = True},
  catchDE2[DiffExp2`Solve`SolveChart[sourceCS, fusionReq, sourceFusion]]];
assert["grouped_cross_lag_inhomogeneous_exact_parity",
  !FailureQ[sourceFused] && !FailureQ[sourceUnfused] &&
  sourceFused["Basis", "Columns"] === sourceUnfused["Basis", "Columns"] &&
  sourceFused["Particular"] === sourceUnfused["Particular"]];

(* Real banana L1 endpoint: the x=1 frame activates one nonconstant group
   (six entries with Q=1+5 eps), while x=0 stays polynomial.  SolveChart's
   always-on residual certificate and delivered windows exercise the path
   through pseudo/fractional analytic-regulator sectors, not just a helper. *)
fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
planU = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, 1}]];
chartU = SelectFirst[Reverse[planU["Charts"]], TrueQ[#["Singular"]] &];
csU = catchDE2[DiffExp2`Solve`PrepareChart[sys, chartU]];
symbolicU = DiffExp2`Solve`Private`clearedSymbolic[csU];
prepU = catchDE2[DiffExp2`Solve`Private`prepareCleared[
  csU, -20, 45, symbolicU]];
statsU = prepU["NhatStats"];
assert["grouped_banana_upper_metadata",
  Max[statsU[[All, "RationalEntries"]]] === 6 &&
  Max[statsU[[All, "LegacyRationalEntries"]]] === 0 &&
  Max[statsU[[All, "RationalGroups"]]] === 1 &&
  Max[statsU[[All, "RationalNumeratorShifts"]]] === 6 &&
  Max[statsU[[All, "LegacyRationalShiftUpperBound"]]] >= 25];

reqU = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TOrder" -> 10|>;
DiffExp2`Solve`ClearSolveCaches[];
solU = Block[{
    DiffExp2`Solve`Private`$disableRationalDenominatorFusion = False},
  catchDE2[DiffExp2`Solve`SolveChart[csU, reqU]]];
(* no clear: the mode-key regression is also exercised on the real chart *)
solUUnfused = Block[{
    DiffExp2`Solve`Private`$disableRationalDenominatorFusion = True},
  catchDE2[DiffExp2`Solve`SolveChart[csU, reqU]]];
assert["grouped_banana_upper_local_solution_windows_and_residual",
  !FailureQ[solU] && !FailureQ[solUUnfused] &&
  solU["Basis", "Columns"] === solUUnfused["Basis", "Columns"] &&
  solU["Basis", "Specs"] === solUUnfused["Basis", "Specs"] &&
  Length[solU["Basis", "Columns"]] === 7 &&
  AllTrue[solU["Basis", "Columns"],
    #["EpsWindow", "CompleteMax"] === 2 &] &&
  TrueQ[solU["Basis", "Diagnostics", "PseudoCollisionsCompensated"]] &&
  Sort[DeleteDuplicates[solU["Basis", "Specs"][[All, "b"]]]] === {0, 1, 2, 3}];

planL = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, 0}]];
chartL = SelectFirst[Reverse[planL["Charts"]], TrueQ[#["Singular"]] &];
csL = catchDE2[DiffExp2`Solve`PrepareChart[sys, chartL]];
prepL = catchDE2[DiffExp2`Solve`Private`prepareCleared[
  csL, -20, 45, DiffExp2`Solve`Private`clearedSymbolic[csL]]];
assert["grouped_banana_lower_stays_polynomial",
  Max[prepL["NhatStats"][[All, "RationalEntries"]]] === 0 &&
  Max[prepL["NhatStats"][[All, "PolynomialShifts"]]] <= 5];

(* A stable benchmark seam: consumers can report these counts without
   materializing the legacy rational tail. *)
Print["GROUPED_BENCH ", <|
  "FrameWidth" -> 45,
  "UpperLegacyShiftUpperBound" ->
    Max[statsU[[All, "LegacyRationalShiftUpperBound"]]],
  "UpperGroupedNumeratorShifts" ->
    Max[statsU[[All, "RationalNumeratorShifts"]]],
  "LowerRationalEntries" ->
    Max[prepL["NhatStats"][[All, "RationalEntries"]]]|>];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["FAILED"]; Quit[1], Print["All tests PASSED!"]; Quit[0]];
