(* End-to-end API tests: the mini-FT closed-form pin. *)
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
  "ExpansionOrder" -> 40, "EpsilonOrder" -> 3}]];
x = Global`x; eps = Global`eps;
esC = DiffExp2`EpsSeries`ESCoefficient;
eqN[a_, b_, tol_] := TrueQ[Abs[N[a - b, 30]] < tol];

algCoord = (-74641 + 23*Sqrt[57983929])/133584;
numCoord = DiffExp2`API`Private`numericLocalCoordinate[algCoord];
symCoord = DiffExp2`API`Private`numericLocalCoordinate[Global`rho + Sqrt[2]];
assert["integration_algebraic_coordinate_grounded_at_2x_wp",
  InexactNumberQ[numCoord] && Precision[numCoord] >= 190 &&
  Abs[N[numCoord - algCoord, 180]] < 10^-180 &&
  symCoord === Global`rho + Sqrt[2]];

(* THE mini-FT pin: f' = (eps/x) f, f(11/23) = (11/23)^eps;
   Int_0^1 f dx = 1/(1+eps) = 1 - eps + eps^2 - eps^3 *)
sys = DiffExp2`API`LoadSystem[<|"Matrix" -> {{eps/x}}, "Variable" -> x|>];
assert["load_system", sys["SingularFactors"] === {x}];

(* LoadSystem retains the full exact matrix/factor data for solving, while
   exposing a separate eps-zero alphabet to the line planner.  In
   particular, the eps-only denominator valuation must not erase x=0. *)
epsAlphabetMatrix = {{1/(eps*(x + eps))}};
epsAlphabetSys = DiffExp2`API`LoadSystem[<|
  "Matrix" -> epsAlphabetMatrix, "Variable" -> x|>];
assert["load_system_epsilon_zero_alphabet_keeps_exact_matrix",
  epsAlphabetSys["Matrix"] === epsAlphabetMatrix &&
  epsAlphabetSys["SingularFactorsExact"] === {x + eps} &&
  epsAlphabetSys["SingularFactors"] === {x}];
epsAlphabetSings = DiffExp2`Transport`FindSingularities[epsAlphabetSys];
assert["load_system_epsilon_zero_alphabet_drives_planner",
  Keys[epsAlphabetSings["Factors"]] === {x} &&
  epsAlphabetSings["All"] === {0}];
epsMovingOnlySys = DiffExp2`API`LoadSystem[<|
  "Matrix" -> {{1/(x + eps)}}, "Variable" -> x|>];
epsMovingEndpoint = catchDE2[DiffExp2`API`TransportEndpoint[
  epsMovingOnlySys, {{1, 0, 0, 0}}, 1/2, 0]];
epsMovingLimit = If[FailureQ[epsMovingEndpoint], epsMovingEndpoint,
  catchDE2[DiffExp2`API`EndpointLimitValues[epsMovingEndpoint, {1}]]];
assert["apparent_moving_matrix_pole_is_exactly_desingularized",
  !FailureQ[epsMovingLimit] &&
  eqN[esC[epsMovingLimit, 0], 0, 10^-30] &&
  eqN[esC[epsMovingLimit, 1], 2, 10^-30] &&
  eqN[esC[epsMovingLimit, 2], -4, 10^-29] &&
  eqN[esC[epsMovingLimit, 3], 8, 10^-28]];
epsMovingNonApparentSys = DiffExp2`API`LoadSystem[<|
  "Matrix" -> {{eps/(x + eps)}}, "Variable" -> x|>];
epsMovingNonApparent = catchDE2[DiffExp2`API`TransportEndpoint[
  epsMovingNonApparentSys, {{1, 0, 0, 0}}, 1/2, 0]];
assert["nonapparent_moving_matrix_pole_at_projected_center_is_loud",
  FailureQ[epsMovingNonApparent] && epsMovingNonApparent["ID"] === "E3"];
bvals = Transpose[Table[{SeriesCoefficient[(11/23)^eps, {eps, 0, k}]}, {k, 0, 3}]];
li = catchDE2[DiffExp2`API`LineIntegral[sys, bvals, 11/23, {0, 1}, {1}]];
assert["line_integral_x_to_eps",
  !FailureQ[li] &&
  (* ~12-digit pipeline precision in v1 (recorded niggle); the pin is the
     VALUE {1,-1,1,-1}, not yet the precision *)
  eqN[esC[li, 0], 1, 10^-10] &&
  eqN[esC[li, 1], -1, 10^-10] &&
  eqN[esC[li, 2], 1, 10^-9]];

(* Tiling ownership must stay inside the same R/2 envelope certified by
   SegmentErrorProbe.  The old 0.9 R clamp assigned the small left chart far
   beyond that bound for this unequal-radius pair; the half-disks meet only
   at 1/4, which is therefore the unique safe breakpoint. *)
mockTiles = catchDE2[DiffExp2`API`Private`halfRadiusTiles[{
  <|"Chart" -> <|"Center" -> -1/9, "Radius" -> 13/18|>|>,
  <|"Chart" -> <|"Center" -> 4/3, "Radius" -> 13/6|>|>}, -1/9, 4/3]];
assert["line_integral_tiles_use_certified_half_radius_overlap",
  !FailureQ[mockTiles] && Length[mockTiles] === 2 &&
  mockTiles[[1, 3]] === 1/4 && mockTiles[[2, 2]] === 1/4];
mockTileGap = catchDE2[DiffExp2`API`Private`halfRadiusTiles[{
  <|"Chart" -> <|"Center" -> 0, "Radius" -> 1|>|>,
  <|"Chart" -> <|"Center" -> 2, "Radius" -> 1|>|>}, 0, 2]];
assert["line_integral_half_radius_tiling_hole_loud",
  FailureQ[mockTileGap] && mockTileGap["ID"] === "E9"];
mockTwoSidedTiles = catchDE2[DiffExp2`API`Private`halfRadiusTiles[
  Map[<|"Chart" -> <|"Center" -> #[[1]], "Radius" -> #[[2]]|>|> &,
    {{4/3, 11/3}, {7/3, 1}, {7/3, 1}, {17/6, 3/2}}],
  4/3, 17/6]];
assert["line_integral_two_sided_duplicate_anchor_tiles_monotone",
  !FailureQ[mockTwoSidedTiles] &&
  And @@ Thread[Differences[Join[{mockTwoSidedTiles[[1, 2]]},
      mockTwoSidedTiles[[All, 3]]]] >= 0]];

(* Cancellation must happen before the endpoint divergence gate.  Each
   component is 1/x and has a divergent integral on [0,1], whereas the
   requested scalar combination is identically zero. *)
sysCancel = DiffExp2`API`LoadSystem[
  <|"Matrix" -> {{-1/x, 0}, {0, -1/x}}, "Variable" -> x|>];
bvCancel = {{2, 0, 0, 0}, {2, 0, 0, 0}};
liBadDim = catchDE2[DiffExp2`API`LineIntegral[
  sysCancel, bvCancel, 1/2, {0, 1}, {1}]];
assert["line_integral_coefficient_dimension_loud",
  FailureQ[liBadDim] && liBadDim["ID"] === "E9"];
liOne = catchDE2[DiffExp2`API`LineIntegral[
  sysCancel, bvCancel, 1/2, {0, 1}, {1, 0}]];
assert["line_integral_individual_endpoint_divergence_loud",
  FailureQ[liOne] && liOne["ID"] === "E2"];
liCancel = catchDE2[DiffExp2`API`LineIntegral[
  sysCancel, bvCancel, 1/2, {0, 1}, {1, -1}]];
assert["line_integral_cross_component_endpoint_cancellation",
  !FailureQ[liCancel] &&
  And @@ Table[eqN[esC[liCancel, k], 0, 10^-30], {k, 0, 3}]];
liCancelLaurent = catchDE2[DiffExp2`API`LineIntegral[
  sysCancel, bvCancel, 1/2, {0, 1}, {1/eps, -1/eps}]];
assert["line_integral_laurent_endpoint_cancellation",
  !FailureQ[liCancelLaurent] &&
  And @@ Table[eqN[esC[liCancelLaurent, k], 0, 10^-30], {k, -1, 2}]];
trCancel = catchDE2[DiffExp2`API`TransportEndpoint[
  sysCancel, bvCancel, 1/2, 0]];
limCancel = If[FailureQ[trCancel], trCancel,
  catchDE2[DiffExp2`API`EndpointLimitValues[trCancel, {1, -1}]]];
assert["endpoint_limit_cross_component_cancellation",
  !FailureQ[limCancel] &&
  And @@ Table[eqN[esC[limCancel, k], 0, 10^-30], {k, 0, 3}]];

(* singular-endpoint limit: 2-component f = {x^eps, const}:
   A = diag(eps/x, 0): limit at 0 of {1,1}.f = the constant *)
sys2 = DiffExp2`API`LoadSystem[<|"Matrix" -> {{eps/x, 0}, {0, 0}}, "Variable" -> x|>];
bv2 = Transpose[Table[{SeriesCoefficient[(11/23)^eps, {eps, 0, k}],
  If[k === 0, 7, 0]}, {k, 0, 3}]];
tr2 = catchDE2[DiffExp2`API`TransportEndpoint[sys2, bv2, 11/23, 0]];
lim2 = catchDE2[DiffExp2`API`EndpointLimitValues[tr2, {1, 1}]];
assert["endpoint_limit_drop_rule",
  !FailureQ[lim2] && eqN[esC[lim2, 0], 7, 10^-20]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
