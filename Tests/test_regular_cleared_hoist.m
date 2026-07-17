(* Exact system-global clearing hoist: regular affine chart equivalence. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100,
  "ExpansionOrder" -> 12, "EpsilonOrder" -> 3}]];
x = Global`x; eps = Global`eps; t = Global`t;

chart[c_, beta_, name_] := <|"ChartVar" -> t, "Center" -> c,
  "Scale" -> beta, "Radius" -> 1, "LocalRadius" -> 1,
  "Singular" -> False, "Name" -> name, "Prescriptions" -> {}|>;
dpoly[data_] := Sum[data["dExpr"][[j + 1]]*t^j,
  {j, 0, data["dD"]}];
npoly[data_] := Sum[data["NhatExpr"][[j + 1]]*t^j,
  {j, 0, data["dN"]}];
sameClearedEquationQ[a_, b_] := AllTrue[
  Flatten[Map[Cancel[Together[#]] &,
    dpoly[a]*npoly[b] - dpoly[b]*npoly[a], {2}]], # === 0 &];

SeedRandom[271828];
randomParity = And @@ Table[Module[
    {a1, a2, p, den, A, sys, beta, cs, legacy, hoist, legacyDense,
     depthLegacy, depthHoist},
    a1 = RandomInteger[{2, 7}]; a2 = RandomInteger[{8, 13}];
    p[r_, c_] := RandomInteger[{-4, 4}] +
      RandomInteger[{-3, 3}] x + RandomInteger[{-2, 2}] x^2 +
      eps*(RandomInteger[{-3, 3}] + RandomInteger[{-2, 2}] x);
    den[r_, c_] := (x + If[OddQ[r + c], a1, a2])*
      (1 + (r + c) eps);
    A = Table[Cancel[Together[p[r, c]/den[r, c]]], {r, 2}, {c, 2}];
    sys = <|"Matrix" -> A, "Variable" -> x|>;
    beta = If[OddQ[seed], 2/5, (1 + Sqrt[2])/5];
    cs = catchDE2[DiffExp2`Solve`PrepareChart[sys,
      chart[1/3, beta, "hoist-random-" <> ToString[seed]]]];
    If[FailureQ[cs], Return[False, Module]];
    legacy = Block[{DiffExp2`Solve`Private`$disableGlobalClearedHoist = True},
      catchDE2[DiffExp2`Solve`Private`clearedSymbolic[cs]]];
    hoist = catchDE2[DiffExp2`Solve`Private`clearedSymbolic[cs]];
    legacyDense = Block[{
        DiffExp2`Solve`Private`$disableGlobalClearedHoist = True,
        DiffExp2`Solve`Private`$disableIdentityNhatShortcut = True},
      catchDE2[DiffExp2`Solve`Private`clearedSymbolic[cs]]];
    If[AnyTrue[{legacy, hoist, legacyDense}, FailureQ], Return[False, Module]];
    depthLegacy = DiffExp2`Solve`Private`recurrencePoleDepth[legacy, 12];
    depthHoist = DiffExp2`Solve`Private`recurrencePoleDepth[hoist, 12];
    sameClearedEquationQ[legacy, hoist] &&
      sameClearedEquationQ[legacy, legacyDense] &&
      depthLegacy === depthHoist && depthLegacy ===
        DiffExp2`Solve`Private`recurrencePoleDepth[legacyDense, 12]],
  {seed, 1, 8}];
assert["regular_global_clear_random_affine_exact_equivalence_8",
  randomParity];

(* An invertible affine map can introduce a coefficient-field unit even
   though it cannot introduce epsilon/regulator content.  The hoist
   deliberately retains that unit instead of spending a polynomial GCD to
   reproduce the legacy normalization. *)
unitSys = <|"Matrix" -> {{1/(x + 1)}}, "Variable" -> x|>;
unitChart = catchDE2[DiffExp2`Solve`PrepareChart[unitSys,
  chart[1, 2, "affine-unit-content"]]];
unitLegacy = Block[{DiffExp2`Solve`Private`$disableGlobalClearedHoist = True},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolic[unitChart]]];
unitHoist = catchDE2[
  DiffExp2`Solve`Private`clearedSymbolic[unitChart]];
assert["regular_global_clear_retains_only_harmless_affine_unit",
  !FailureQ[unitLegacy] && !FailureQ[unitHoist] &&
    sameClearedEquationQ[unitLegacy, unitHoist] &&
    Cancel[Together[dpoly[unitHoist]/dpoly[unitLegacy]]] === 2 &&
    DiffExp2`Solve`Private`recurrencePoleDepth[unitLegacy, 12] ===
      DiffExp2`Solve`Private`recurrencePoleDepth[unitHoist, 12]];

(* Analytic regulators stay in the exact coefficient field: the hoist must
   neither numericize them nor treat them as chart geometry. *)
rho = Global`rho;
regSys = <|"Matrix" -> {
    {(1 + rho*x)/(eps*(2 + x)), (rho + eps*x)/(3 - x)},
    {eps/(2 + x), rho*x/(3 - x)}}, "Variable" -> x|>;
regChart = catchDE2[DiffExp2`Solve`PrepareChart[regSys,
  chart[2/7, (1 + Sqrt[3])/9, "analytic-regulator"]]];
regLegacy = Block[{DiffExp2`Solve`Private`$disableGlobalClearedHoist = True},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolic[regChart]]];
regHoist = catchDE2[DiffExp2`Solve`Private`clearedSymbolic[regChart]];
assert["regular_global_clear_preserves_exact_analytic_regulator",
  !FailureQ[regLegacy] && !FailureQ[regHoist] &&
    sameClearedEquationQ[regLegacy, regHoist] &&
    !FreeQ[regHoist, rho] && FreeQ[regHoist, _?InexactNumberQ] &&
    DiffExp2`Solve`Private`recurrencePoleDepth[regLegacy, 12] ===
      DiffExp2`Solve`Private`recurrencePoleDepth[regHoist, 12]];

(* Singular/gauged charts are deliberately outside the optimized dispatch;
   their legacy cleared construction must remain byte-for-byte identical. *)
gaugeSys = <|"Matrix" -> {{0, x^-2}, {0, 0}}, "Variable" -> x|>;
gaugeChart = catchDE2[DiffExp2`Solve`PrepareChart[gaugeSys,
  chart[0, 1, "singular-gauged-legacy"]]];
gaugeWrapped = catchDE2[DiffExp2`Solve`Private`clearedSymbolic[gaugeChart]];
gaugeLegacy = catchDE2[
  DiffExp2`Solve`Private`clearedSymbolicLegacy[gaugeChart]];
assert["singular_gauged_clearing_path_unchanged",
  !FailureQ[gaugeWrapped] && !FailureQ[gaugeLegacy] &&
    !TrueQ[DiffExp2`Solve`Private`regularIdentityFrameQ[gaugeChart]] &&
    gaugeWrapped === gaugeLegacy];

(* Pole depth cannot be cached by system alone.  The same A=x/eps has a
   t^2/eps first lag at center 0 but a t/eps first lag at center 1. *)
counterSys = <|"Matrix" -> {{x/eps}}, "Variable" -> x|>;
counter0 = catchDE2[DiffExp2`Solve`PrepareChart[counterSys,
  chart[0, 1, "pole-depth-center-0"]]];
counter1 = catchDE2[DiffExp2`Solve`PrepareChart[counterSys,
  chart[1, 1, "pole-depth-center-1"]]];
counterData0 = catchDE2[DiffExp2`Solve`Private`clearedSymbolic[counter0]];
counterData1 = catchDE2[DiffExp2`Solve`Private`clearedSymbolic[counter1]];
counterDepth0 = DiffExp2`Solve`Private`recurrencePoleDepth[counterData0, 6];
counterDepth1 = DiffExp2`Solve`Private`recurrencePoleDepth[counterData1, 6];
assert["regular_pole_depth_is_chart_specific_counterexample",
  counterDepth0 === 3 && counterDepth1 === 6];

(* Cache keys include the exact system and affine chart, and ClearSolveCaches
   owns every new cache's lifetime. *)
cacheCountsBefore = Length /@ {
  DiffExp2`Solve`Private`$systemClearRegistry,
  DiffExp2`Solve`Private`$globalClearedCache,
  DiffExp2`Solve`Private`$chartClearedCache};
DiffExp2`Solve`ClearSolveCaches[];
cacheCountsAfter = Length /@ {
  DiffExp2`Solve`Private`$systemClearRegistry,
  DiffExp2`Solve`Private`$globalClearedCache,
  DiffExp2`Solve`Private`$chartClearedCache};
assert["regular_global_clear_cache_keys_and_lifetime",
  Min[cacheCountsBefore] >= 2 && cacheCountsAfter === {0, 0, 0}];

(* PrepareChart's own cache must not alias systems that happen to have the
   same expression matrix but use distinct independent-variable symbols. *)
y = Global`y;
sameMatrixX = <|"Matrix" -> {{1/(1 + x)}}, "Variable" -> x|>;
sameMatrixY = <|"Matrix" -> {{1/(1 + x)}}, "Variable" -> y|>;
sameChartX = catchDE2[DiffExp2`Solve`PrepareChart[sameMatrixX,
  chart[2, 1, "cache-variable-x"]]];
sameChartY = catchDE2[DiffExp2`Solve`PrepareChart[sameMatrixY,
  chart[2, 1, "cache-variable-y"]]];
assert["prepare_chart_cache_separates_system_variables",
  !FailureQ[sameChartX] && !FailureQ[sameChartY] &&
    sameChartX["SystemClearKey"] =!= sameChartY["SystemClearKey"] &&
    sameChartX["ThetaMatrix"] =!= sameChartY["ThetaMatrix"]];

Print["Regular cleared hoist tests: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1]];
