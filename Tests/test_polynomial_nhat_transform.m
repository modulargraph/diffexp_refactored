(* Exact polynomial-first spectral transform for singular/gauged charts. *)
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
t = Global`t; x = Global`x; eps = Global`eps; rho = Global`rho;

npoly[data_] := Sum[data["NhatExpr"][[j + 1]]*t^j,
  {j, 0, data["dN"]}];
entrywiseExactQ[a_, b_] := Dimensions[a] === Dimensions[b] &&
  AllTrue[Flatten[Map[Cancel[Together[#]] &, a - b, {2}]], # === 0 &];
sameNhatQ[a_, b_] := a["dN"] === b["dN"] &&
  entrywiseExactQ[npoly[a], npoly[b]];
canonicalNhatQ[data_] := AllTrue[Flatten[data["NhatExpr"]],
  Cancel[Together[#]] === # &];
syntheticCS[num_, V_, name_] := <|
  "ChartVar" -> t, "Center" -> name, "SystemSize" -> Length[V],
  "ThetaMatrix" -> num, "V" -> V,
  "VInv" -> Map[Cancel[Together[#]] &, Inverse[V], {2}],
  "IndicialData" -> <|"Regular" -> False|>|>;

SeedRandom[314159];
randomParity = And @@ Table[Module[
    {rad, left, right, base, V, coeff, num, cs, old, new},
    rad = If[OddQ[seed], 1, Sqrt[2]];
    left = {{1, RandomInteger[{-2, 2}]}, {0, 1}};
    right = {{1, 0}, {RandomInteger[{-2, 2}], 1}};
    (* det(base) = rad eps, hence the exact inverse has an epsilon pole. *)
    base = {{rad, rad}, {rho, rho + eps}};
    V = Map[Cancel[Together[#]] &, left . base . right, {2}];
    coeff[r_, c_, j_] := RandomInteger[{-3, 3}] +
      eps RandomInteger[{-2, 2}] + rho RandomInteger[{-2, 2}] +
      If[j === 1 && r === c, RandomInteger[{1, 3}]/eps, 0] +
      If[EvenQ[seed] && j === 2 && r === 1 && c === 2, Sqrt[3]/5, 0];
    num = Table[Sum[coeff[r, c, j]*t^j, {j, 0, 3}], {r, 2}, {c, 2}];
    cs = syntheticCS[num, V, "random-" <> ToString[seed]];
    old = Block[{
        DiffExp2`Solve`Private`$disablePolynomialNhatTransform = True},
      catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[cs]]];
    new = Block[{
        DiffExp2`Solve`Private`$disablePolynomialNhatTransform = False},
      catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[cs]]];
    !FailureQ[old] && !FailureQ[new] &&
      DiffExp2`Solve`Private`matrixEpsPoleDepth[cs["VInv"]] >= 1 &&
      old["dExpr"] === new["dExpr"] && sameNhatQ[old, new] &&
      canonicalNhatQ[new] && FreeQ[new, _?InexactNumberQ] &&
      DiffExp2`Solve`Private`recurrencePoleDepth[old, 12] ===
        DiffExp2`Solve`Private`recurrencePoleDepth[new, 12]],
  {seed, 1, 12}];
assert["polynomial_nhat_random_exact_eps_pole_inverse_parity_12",
  randomParity];

(* Pin radical-field parity separately: operation order may change the
   syntax, but every coefficient difference must cancel exactly. *)
radV = {{Sqrt[2], Sqrt[2]}, {rho, rho + eps}};
radNum = {{(1 + rho/eps) + (2 + Sqrt[3]) t + eps t^3,
    (rho - eps) t + Sqrt[5] t^2},
  {1/eps + rho t^2, 2 - rho + (1 + eps) t^3}};
radCS = syntheticCS[radNum, radV, "radical-parity"];
radOld = Block[{
    DiffExp2`Solve`Private`$disablePolynomialNhatTransform = True},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[radCS]]];
radNew = Block[{
    DiffExp2`Solve`Private`$disablePolynomialNhatTransform = False},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[radCS]]];
assert["polynomial_nhat_radical_entrywise_cancel_together_parity",
  !FailureQ[radOld] && !FailureQ[radNew] && sameNhatQ[radOld, radNew] &&
    (!FreeQ[radNew, Sqrt[2]] || !FreeQ[radNew, Sqrt[3]] ||
      !FreeQ[radNew, Sqrt[5]]) && !FreeQ[radNew, rho] &&
    canonicalNhatQ[radNew]];

(* The optimization is valid only for the residue frame built by
   PrepareChart.  A chart-variable-dependent synthetic frame is rejected
   loudly rather than silently taking a different algebraic route. *)
badCS = Join[radCS, <|"Center" -> "t-dependent-frame",
  "V" -> (radCS["V"] /. rho -> rho + t)|>];
bad = Block[{
    DiffExp2`Solve`Private`$disablePolynomialNhatTransform = False},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[badCS]]];
assert["polynomial_nhat_t_dependent_frame_precondition_loud",
  FailureQ[bad] && bad["ID"] === "E5"];

(* A real rank-reduced/gauged ChartSystem exercises the same optimized
   dispatch without moving V through the t-dependent gauge. *)
gaugeChart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 1, "LocalRadius" -> 1, "Singular" -> True,
  "Name" -> "polynomial-nhat-gauge", "Prescriptions" -> {}|>;
gaugeCS = catchDE2[DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0, x^-2}, {0, 0}}, "Variable" -> x|>, gaugeChart]];
gaugeOld = Block[{
    DiffExp2`Solve`Private`$disablePolynomialNhatTransform = True},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[gaugeCS]]];
gaugeNew = Block[{
    DiffExp2`Solve`Private`$disablePolynomialNhatTransform = False},
  catchDE2[DiffExp2`Solve`Private`clearedSymbolicLegacy[gaugeCS]]];
assert["polynomial_nhat_rank_reduced_gauge_exact_parity",
  !FailureQ[gaugeCS] && !FailureQ[gaugeOld] && !FailureQ[gaugeNew] &&
    gaugeCS["Gauge"] =!= IdentityMatrix[2] && sameNhatQ[gaugeOld, gaugeNew] &&
    DiffExp2`Solve`Private`recurrencePoleDepth[gaugeOld, 12] ===
      DiffExp2`Solve`Private`recurrencePoleDepth[gaugeNew, 12]];

Print["Polynomial Nhat transform tests: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1]];
