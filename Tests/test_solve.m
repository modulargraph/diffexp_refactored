(* Focused closed-form unit tests for DiffExp2/Solve.m. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100}]];

pc = DiffExp2`Solve`PrepareChart;
sc = DiffExp2`Solve`SolveChart;
sp = DiffExp2`Solve`SolveParticular;
t = Global`t; eps = Global`eps; x = Global`x;
req[kmin_, kmax_, nord_] := <|"EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "TOrder" -> nord|>;
mkSys[A_] := <|"Matrix" -> A, "Variable" -> x|>;
mkChart[c_, r_, name_] := <|"ChartVar" -> t, "Center" -> c, "Radius" -> r, "Name" -> name|>;
rowOf[ls_, sec_, k_] := sec["Coeffs"][[k - ls["EpsWindow", "Min"] + 1]];
tagsOf[ls_] := {#["a"], #["b"], #["p"]} & /@ ls["Sectors"];

resolvedFrameZero = SetPrecision[0., 250];
unresolvedFrameZero = 0``17;
assert["su00_frame_chop_exactifies_only_certified_centered_zero",
  DiffExp2`Solve`Private`certifiedFrameChop[resolvedFrameZero] === 0 &&
  SameQ[DiffExp2`Solve`Private`certifiedFrameChop[unresolvedFrameZero],
    unresolvedFrameZero]];
assert["su00_frame_chop_uses_configured_absolute_floor",
  DiffExp2`Solve`Private`certifiedFrameChop[
    SetPrecision[10^-80, 200]] === 0 &&
  DiffExp2`Solve`Private`certifiedFrameChop[
    SetPrecision[10^-40, 200]] =!= 0];
assert["su00_frame_chop_preserves_symbolic_regulator",
  DiffExp2`Solve`Private`certifiedFrameChop[Global`rho] === Global`rho];

(* SU-01: f' = f *)
cs1 = pc[mkSys[{{1}}], mkChart[0, 10, "su01"]];
r1 = sc[cs1, req[0, 2, 8]];
col1 = r1["Basis"]["Columns"][[1]];
assert["su01_regular_exponential",
  tagsOf[col1] === {{0, 0, 0}} &&
  rowOf[col1, col1["Sectors"][[1]], 0][[;; 6, 1]] === Table[1/n!, {n, 0, 5}]];

(* Residual comparison starts at the union of supports: the derivative of
   this wrong constant is zero, but B.f has a nonzero eps^0 term. *)
badConst1 = <|"Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 10,
  "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {1}, {0}],
      {k, 0, 2}, {n, 0, 4}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 4|>,
  "ErrorEstimate" -> {0, 0, 0}, "Prescriptions" -> {}|>;
badResid1 = catchDE2[DiffExp2`Solve`ODEResidualCheck[cs1, badConst1]];
assert["su01_residual_checks_certified_zero_lower_orders",
  FailureQ[badResid1] && badResid1["ID"] === "E7"];

(* Raw Abs can collapse this nonzero coefficient because its centered-zero
   imaginary component has poor accuracy.  The ODE proof must stay loud. *)
complexBad1 =
  4.4267459561002104836`0.2683567342206439 +
    0``-0.022790862816694443*I;
badComplex1 = Join[badConst1, <|"Sectors" -> {
  <|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {complexBad1}, {0}],
      {k, 0, 2}, {n, 0, 4}]|>}|>];
badComplexResid1 = catchDE2[
  DiffExp2`Solve`ODEResidualCheck[cs1, badComplex1]];
assert["su01_complex_residual_magnitude_not_poisoned",
  FailureQ[badComplexResid1] && badComplexResid1["ID"] === "E7"];

mkBadEps1[value_] := Join[badConst1, <|"Sectors" -> {
  <|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 1 && n === 0, {value}, {0}],
      {k, 0, 2}, {n, 0, 4}]|>}|>];
badDiagonalResid1 = catchDE2[DiffExp2`Solve`ODEResidualCheck[cs1,
  mkBadEps1[(5 + 5 I)*10^-9]]];
badUncertainZeroResid1 = catchDE2[DiffExp2`Solve`ODEResidualCheck[cs1,
  mkBadEps1[0``5]]];
assert["su01_residual_uses_true_modulus_and_uncertainty_upper_bound",
  FailureQ[badDiagonalResid1] && badDiagonalResid1["ID"] === "E7" &&
  FailureQ[badUncertainZeroResid1] &&
  badUncertainZeroResid1["ID"] === "E7"];

(* A negative-valued Taylor multiplier composes at every recurrence step.
   f' = f/(eps(1-eps)) has f = exp[t/(eps(1-eps))], hence its t^4
   coefficient starts at eps^-4 and remains complete through eps^2 when
   the exact normalized datum is solved to the requested order. *)
cs1p = pc[mkSys[{{1/(eps*(1 - eps))}}], mkChart[0, 2, "su01_repeated_pole"]];
r1p = catchDE2[sc[cs1p, req[0, 2, 4]]];
assert["su01_repeated_negative_multiplier_budget",
  !FailureQ[r1p] && Module[{col = First[r1p["Basis"]["Columns"]], sec},
    sec = First[col["Sectors"]];
    col["EpsWindow"] === <|"Min" -> -4, "CompleteMax" -> 2|> &&
    rowOf[col, sec, -4][[5, 1]] === 1/24 &&
    rowOf[col, sec, 2][[5, 1]] === 7/2]];

(* Both sides of the physical frame fail through the named E4 path, rather
   than leaking negative-length Table/Part diagnostics. *)
e1hi = catchDE2[DiffExp2`Solve`Private`ratEpsList[eps^5, eps, -1, 3]];
e1inv = catchDE2[DiffExp2`Solve`Private`frInv[{1, 0}, -3, 2]];
assert["su01_work_frame_range_guards",
  FailureQ[e1hi] && e1hi["ID"] === "E4" &&
  FailureQ[e1inv] && e1inv["ID"] === "E4"];

(* Exact zero/monomial frame products must be coefficient-identical to the
   dense convolution they replace, including negative shifts and upper
   truncation.  The strict lower-frame guard remains active. *)
fb1conv = -2; W1conv = 7;
legacyConv1[a_, b_] :=
  Take[ListConvolve[a, b, {1, -1}, 0], {1 - fb1conv, W1conv - fb1conv}];
dense1conv = {0, 0, 3, -2, 5, 0, 7};
negMono1conv = {0, 11, 0, 0, 0, 0, 0};
posMono1conv = {0, 0, 0, 0, 13, 0, 0};
zero1conv = ConstantArray[0, W1conv];
fastNeg1conv = DiffExp2`Solve`Private`frConv[
  negMono1conv, dense1conv, fb1conv, W1conv];
fastPos1conv = DiffExp2`Solve`Private`frConv[
  dense1conv, posMono1conv, fb1conv, W1conv];
fastZero1conv = DiffExp2`Solve`Private`frConv[
  zero1conv, dense1conv, fb1conv, W1conv];
badConv1 = catchDE2[DiffExp2`Solve`Private`frConv[
  ReplacePart[zero1conv, 1 -> 1], ReplacePart[zero1conv, 1 -> 1],
  fb1conv, W1conv]];
assert["su01_zero_monomial_frame_convolution_exact",
  fastNeg1conv === legacyConv1[negMono1conv, dense1conv] &&
  fastPos1conv === legacyConv1[dense1conv, posMono1conv] &&
  fastZero1conv === zero1conv &&
  FailureQ[badConv1] && badConv1["ID"] === "E4"];

(* SU-24: exact affine/Jordan frame solves.  testFrame24 and the local
   legacy copy exist only to pin the private recurrence seam: the new
   implementation must agree with the former inverse-power construction
   wherever the latter's finite upper frame contains all needed terms. *)
testFrame24[fb_Integer, W_Integer, terms_List] := ReplacePart[
  ConstantArray[0, W], Map[(#[[1]] - fb + 1) -> #[[2]] &, terms]];
legacyBlockSolve24[rhs_, dA_, dB_, invD0_, q_, fb_, W_] := Module[
  {delta, inv, z},
  delta = testFrame24[fb, W, Select[{{0, dA}, {1, dB}}, #[[2]] =!= 0 &]];
  inv = DiffExp2`Solve`Private`frInv[delta, fb, W];
  z = Table[Module[{acc = ConstantArray[0, W], pw = inv},
      Do[
        If[r + m <= q,
          acc += DiffExp2`Solve`Private`frConv[pw, rhs[[r + m]], fb, W]];
        If[m < q - 1,
          pw = DiffExp2`Solve`Private`frConv[pw, inv, fb, W]],
        {m, 0, q - 1}];
      DiffExp2`Solve`Private`frConv[acc, invD0, fb, W]],
    {r, q}];
  z];

(* d0^-1 is applied once after the Jordan solve, not at every traversal
   of the superdiagonal. *)
bt24scalar = DiffExp2`Solve`Private`blockSolveTPFrame[
  {{0}, {1}}, 3, 0, None, 1/2, 2, 0, 1];
assert["su24_case_t_d0_scalar_applied_once",
  bt24scalar === {{1/18}, {1/6}}];

fb24t = -3; W24t = 14; q24t = 4;
rhs24t = Table[testFrame24[fb24t, W24t,
    Table[{k, (7 + r + 2 k)/(11 + r + k)}, {k, -3, 2}]], {r, q24t}];
delta24t = testFrame24[fb24t, W24t, {{0, 3}, {1, -2}}];
invD024t = testFrame24[fb24t, W24t, {{0, 3/5}}];
fast24t = DiffExp2`Solve`Private`blockSolveTPFrame[
  rhs24t, 3, -2, None, 3/5, q24t, fb24t, W24t];
legacy24t = legacyBlockSolve24[
  rhs24t, 3, -2, invD024t, q24t, fb24t, W24t];
op24t = Table[(5/3) (
    DiffExp2`Solve`Private`frConv[delta24t, fast24t[[r]], fb24t, W24t] -
      If[r < q24t, fast24t[[r + 1]], ConstantArray[0, W24t]]),
  {r, q24t}];
assert["su24_case_t_affine_q4_operator_identity",
  op24t === rhs24t];
(* rhs24t starts at eps^-3, while the legacy inverse frame stops at the
   common frame top.  Its first W-3 slots are therefore the certified
   comparison region; the final three require inverse coefficients beyond
   that old frame and deliberately are not used as equivalence evidence. *)
assert["su24_case_t_affine_q4_legacy_safe_interior",
  fast24t[[All, 1 ;; W24t - 3]] === legacy24t[[All, 1 ;; W24t - 3]]];

fb24p = -5; W24p = 10; q24p = 5;
zero24p = ConstantArray[0, W24p];
rhs24p = ReplacePart[ConstantArray[zero24p, q24p],
  {q24p -> testFrame24[fb24p, W24p, {{0, 1}}]}];
fast24p = DiffExp2`Solve`Private`blockSolveTPFrame[
  rhs24p, 0, -2, None, 1/3, q24p, fb24p, W24p];
expected24p = Table[testFrame24[fb24p, W24p,
    {{-(q24p - r + 1), (1/3) (-2)^(-(q24p - r + 1))}}],
  {r, q24p}];
delta24p = testFrame24[fb24p, W24p, {{1, -2}}];
op24p = Table[3 (
    DiffExp2`Solve`Private`frConv[delta24p, fast24p[[r]], fb24p, W24p] -
      If[r < q24p, fast24p[[r + 1]], zero24p]), {r, q24p}];
assert["su24_case_p_q5_signs_and_full_frame_identity",
  fast24p === expected24p && op24p === rhs24p];

fb24bad = -4; W24bad = 9;
zero24bad = ConstantArray[0, W24bad];
rhs24bad = ReplacePart[ConstantArray[zero24bad, q24p],
  {q24p -> testFrame24[fb24bad, W24bad, {{0, 1}}]}];
bad24p = catchDE2[DiffExp2`Solve`Private`blockSolveTPFrame[
  rhs24bad, 0, -2, None, 1/3, q24p, fb24bad, W24bad]];
assert["su24_case_p_q5_insufficient_lower_halo_loud",
  FailureQ[bad24p] && bad24p["ID"] === "E4"];

(* An epsilon-dependent d0 cannot use the scalar shortcut.  This ODE has
   the normalized exact solution (1+eps+eps^2+x)/(1+eps+eps^2), so only
   its x^1 coefficient is nonzero and its epsilon sequence is periodic. *)
cs24d0 = pc[mkSys[{{1/(1 + eps + eps^2 + x)}}],
  mkChart[0, 1/2, "su24_eps_dependent_d0"]];
r24d0 = catchDE2[sc[cs24d0, req[0, 5, 4]]];
assert["su24_eps_dependent_d0_fallback_end_to_end",
  !FailureQ[r24d0] && Module[{col = First[r24d0["Basis", "Columns"]], sec},
    sec = First[col["Sectors"]];
    Table[rowOf[col, sec, k][[2, 1]], {k, 0, 5}] ===
      {1, -1, 0, 1, -1, 0} &&
    AllTrue[Flatten[Table[rowOf[col, sec, k][[3 ;; 5, 1]], {k, 0, 5}]],
      # === 0 &]]];

(* The automatic residual probe must remain strictly inside a singular
   chart even at very low Taylor order; the former fixed 10^-3
   rationalization tolerance rounded it to the origin. *)
cs1low = pc[mkSys[{{eps/x}}], mkChart[0, 1, "su01_low_order_probe"]];
r1low = catchDE2[sc[cs1low, req[0, 1, 0]]];
assert["su01_low_order_residual_probe_nonzero", !FailureQ[r1low]];

(* SU-04: singular nonresonant 2x2: theta f = (M0 + t N1) f,
   M0 = diag(-1+eps, 2eps), N1 = {{0,0},{1,0}}.
   A = (M0 + t N1)/t.  Column 1 (root a=-1,b=1): the second component obeys
   theta f2 = 2eps f2 + t f1 with f1 = t^(-1+eps): particular component
   f2 = t^eps/(1 - eps) (from (a+n) offset: n=1: delta = (−1+1) + (1−2)eps = −eps;
   recompute: lambda+n − lambda_2 = (−1+eps+1) − 2eps = −eps... 1/(−eps)?? ) *)
A4 = (DiagonalMatrix[{-1 + eps, 2 eps}] + t*{{0, 0}, {1, 0}})/t;
cs4 = pc[mkSys[A4 /. t -> x], mkChart[0, 1, "su04"]];
r4 = catchDE2[sc[cs4, req[-2, 3, 6]]];
assert["su04_pseudo_resonant_2x2_runs",
  !FailureQ[r4] && Length[r4["Basis"]["Columns"]] === 2 &&
  Length[r4["Basis"]["Diagnostics"]["PseudoCollisionsHit"]] >= 1];
(* the collision quotient: column with tag (-1,1): component 2 at n=1 gets
   1/((b−b2)eps) = 1/(−eps): Laurent content at eps^-1 *)
col4 = SelectFirst[r4["Basis"]["Columns"], MemberQ[tagsOf[#], {-1, 1, 0}] &];
assert["su04_laurent_window",
  col4["EpsWindow", "Min"] <= -1];
ev4 = DiffExp2`SectorSeries`EvaluateLocalSolution[col4, 1/3,
  "UsePade" -> False]["Value"];
assert["su04_joint_pseudo_compensation",
  MemberQ[tagsOf[col4], {0, 2, 0}] &&
  Length[r4["Basis"]["Diagnostics"]["PseudoCompensations"]] >= 1 &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[ev4, -1][[2]], 50]] < 10^-40 &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[ev4, 0][[2]] - Log[1/3], 50]] <
    10^-40];
badCol4 = Join[col4, <|"Sectors" -> Map[Function[sec,
    If[{sec["a"], sec["b"], sec["p"]} === {0, 2, 0},
      Join[sec, <|"Coeffs" -> ReplacePart[sec["Coeffs"],
        {-1 - col4["EpsWindow", "Min"] + 1, 1, 2} ->
          sec["Coeffs"][[
            -1 - col4["EpsWindow", "Min"] + 1, 1, 2]] + 1/10]|>],
      sec]], col4["Sectors"]]|>];
badCert4 = catchDE2[DiffExp2`Solve`Private`certifyPseudoCompensation[
  cs4, badCol4, r4["Basis", "Diagnostics", "PseudoCollisionsHit"],
  "su04_corrupted"]];
assert["su04_pseudo_compensation_certificate_loud",
  FailureQ[badCert4] && badCert4["ID"] === "E5"];

(* SU-05: true resonance Jordan chain: B = M0 = {{0,1},{0,0}} *)
cs5 = pc[mkSys[{{0, 1}, {0, 0}}/x], mkChart[0, 1, "su05"]];
r5 = catchDE2[sc[cs5, req[0, 2, 4]]];
assert["su05_jordan_chain",
  !FailureQ[r5] && Length[r5["Basis"]["Columns"]] === 2 &&
  Sort[#["ChainPos"] & /@ r5["Basis"]["Specs"]] === {0, 1}];
(* chain-pos-1 column = log solution: value f = (e1 Log t + e2)?? theta f = M0 f:
   f = {Log t, 1}: in our basis: sector p=1 row k=-1 carries the top weight *)
col5 = r5["Basis"]["Columns"][[SelectFirst[Range[2],
  r5["Basis"]["Specs"][[#]]["ChainPos"] === 1 &]]];
assert["su05_log_member_window",
  col5["EpsWindow", "Min"] === -1 &&
  MemberQ[tagsOf[col5], {0, 0, 1}]];

(* SU-06: integer-spaced resonance log bump: M0 = diag(0,2),
   N(t) = M0 + (t + t^2) {{0,0},{1,0}}: f1 = 1 (const, a=0 root);
   theta f2 = 2 f2 + (t + t^2) f1: n=1: (0+1-2) = -1: f2_1 = -(1);
   n=2: delta = 0 RESONANT: log bump: t^2 coefficient: eps-division *)
A6 = (DiagonalMatrix[{0, 2}] + (t + t^2)*{{0, 0}, {1, 0}})/t;
cs6 = pc[mkSys[A6 /. t -> x], mkChart[0, 1, "su06"]];
r6 = catchDE2[sc[cs6, req[0, 2, 5]]];
assert["su06_log_bump_runs", !FailureQ[r6]];
(* the a=0 column must carry a p=1 sector (the bump) with nonzero content *)
col6 = r6["Basis"]["Columns"][[SelectFirst[Range[2],
  r6["Basis"]["Specs"][[#]]["a"] === 0 &]]];
assert["su06_log_bump_sector",
  MemberQ[tagsOf[col6], {0, 0, 1}] || MemberQ[tagsOf[col6], {2, 0, 1}]];

(* SU-07: banana pseudo-resonance class: diag(0, -1+eps, 2eps), pure diagonal *)
cs7 = pc[mkSys[DiagonalMatrix[{0, -1 + eps, 2 eps}]/x], mkChart[0, 1, "su07"]];
r7 = catchDE2[sc[cs7, req[0, 2, 4]]];
assert["su07_banana_class",
  !FailureQ[r7] && Length[r7["Basis"]["Columns"]] === 3 &&
  First[cs7["Families"]]["CollisionDepth"] === 1 &&
  (* diagonal system: no actual coupling, columns are pure monomials *)
  AllTrue[r7["Basis"]["Columns"], Length[#["Sectors"]] === 1 &] &&
  (* CASE-P operators are visited, but their RHS is structurally zero. *)
  AllTrue[r7["Basis"]["Columns"], #["EpsWindow", "CompleteMax"] === 2 &]];

(* Exact L1 banana indicial pattern: 24 pair records collapse to one
   executed n=1 dependency layer; the ten CASE-P visits are parallel. *)
cs7b = pc[mkSys[DiagonalMatrix[{3 eps, 2 eps, 2 eps, 2 eps, 0,
    -1 + eps, -1 + eps}]/x], mkChart[0, 1, "su07_banana_l1"]];
r7b = catchDE2[sc[cs7b, req[0, 2, 3]]];
assert["su07_banana_l1_parallel_collision_budget",
  Length[First[cs7b["Families"]]["Collisions"]] === 24 &&
  First[cs7b["Families"]]["CollisionDepth"] === 1 &&
  !FailureQ[r7b] &&
  Length[r7b["Basis"]["Diagnostics"]["PseudoCollisionsHit"]] === 10 &&
  AllTrue[r7b["Basis"]["Columns"],
    #["EpsWindow", "CompleteMax"] === 2 &]];

(* SU-08: resonant-source log bump (the old "empty particular" hole):
   scalar theta f = 0·f + source t^0 (i.e. theta f = s, s = 1 at n=0):
   exact solution f = Log t = (eps Log t)/eps: sector (0,0,1) row k=-1. *)
cs8 = pc[mkSys[{{0}}/x (* theta f = 0 *) * 0 + {{0}}], mkChart[0, 1, "su08"]];
(* build theta-form source: sector (0,0,0), coeff 1 at n=0, window [0,2] *)
src8 = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[{If[k === 0 && n === 0, 1, 0]}, {k, 0, 2}, {n, 0, 3}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 3|>|>;
p8 = catchDE2[sp[cs8, src8, req[0, 2, 3]]];
assert["su08_resonant_source_log_bump",
  !FailureQ[p8] &&
  MemberQ[tagsOf[p8], {0, 0, 1}] &&
  Module[{sec = SelectFirst[p8["Sectors"], #["p"] === 1 &]},
    rowOf[p8, sec, -1][[1, 1]] === 1]];

(* SU-09: pseudo-resonant source: scalar theta f = 2eps f + 1.  The raw
   particular -1/(2eps) is jointly completed by +t^(2eps)/(2eps), giving
   the canonical eps-regular value (t^(2eps)-1)/(2eps). *)
cs9 = pc[mkSys[{{2 eps}}/x], mkChart[0, 1, "su09"]];
p9 = catchDE2[sp[cs9, src8, req[0, 2, 3]]];
ev9 = If[FailureQ[p9], None,
  DiffExp2`SectorSeries`EvaluateLocalSolution[p9, 1/3,
    "UsePade" -> False]["Value"]];
assert["su09_pseudo_resonant_source",
  !FailureQ[p9] && p9["EpsWindow", "Min"] === -1 &&
  p9["EpsWindow", "CompleteMax"] === 1 &&
  MemberQ[tagsOf[p9], {0, 2, 0}] &&
  p9["Diagnostics", "PseudoCollisionsCompensated"] === True &&
  Length[p9["Diagnostics", "PseudoCompensations"]] >= 1 &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[ev9, -1][[1]], 50]] < 10^-40 &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[ev9, 0][[1]] - Log[1/3], 50]] <
    10^-40];

(* Parallel pseudo hits at the same n are one dependency layer, not three
   sequential losses.  theta f = diag(eps,2eps,3eps).f + {1,1,1}. *)
cs9p = pc[mkSys[DiagonalMatrix[{eps, 2 eps, 3 eps}]/x],
  mkChart[0, 1, "su09_parallel"]];
src9p = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[
      If[k === 0 && n === 0, {1, 1, 1}, {0, 0, 0}],
      {k, 0, 3}, {n, 0, 2}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 3|>,
  "TWindow" -> <|"CompleteMax" -> 2|>|>;
p9p = catchDE2[sp[cs9p, src9p, req[0, 3, 2]]];
assert["su09_parallel_pseudo_one_order_loss",
  !FailureQ[p9p] && p9p["EpsWindow"] === <|"Min" -> -1, "CompleteMax" -> 2|> &&
  rowOf[p9p, First[p9p["Sectors"]], -1][[1]] === {-1, -1/2, -1/3}];

(* A size-2 pseudo-resonant Jordan block has pole depth exactly q=2. *)
cs9q = pc[mkSys[{{2 eps, 1}, {0, 2 eps}}/x], mkChart[0, 1, "su09_q2"]];
src9q = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[
      If[k === 0 && n === 0, {0, 1}, {0, 0}],
      {k, 0, 3}, {n, 0, 2}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 3|>,
  "TWindow" -> <|"CompleteMax" -> 2|>|>;
p9q = catchDE2[sp[cs9q, src9q, req[0, 3, 2]]];
assert["su09_pseudo_jordan_q2_window",
  !FailureQ[p9q] && p9q["EpsWindow"] === <|"Min" -> -2, "CompleteMax" -> 1|> &&
  rowOf[p9q, First[p9q["Sectors"]], -2][[1, 1]] === 1/4];

(* Particular work frames include the full pseudo-Jordan pole depth. *)
j5 = DiagonalMatrix[ConstantArray[2 eps, 5]] +
  DiagonalMatrix[ConstantArray[1, 4], 1];
cs9q5 = pc[mkSys[j5/x], mkChart[0, 1, "su09_q5_particular"]];
src9q5 = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{0, 0, 0, 0, 1}, {0, 0, 0, 0, 0}}}|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 0|>,
  "TWindow" -> <|"CompleteMax" -> 1|>|>;
p9q5 = catchDE2[sp[cs9q5, src9q5, req[0, 0, 1]]];
assert["su09_pseudo_jordan_q5_particular_budget",
  !FailureQ[p9q5] &&
  p9q5["EpsWindow"] === <|"Min" -> -5, "CompleteMax" -> -5|> &&
  rowOf[p9q5, First[p9q5["Sectors"]], -5][[1, 1]] === -1/32];

(* SU-10: empty source -> zero particular with full windows *)
p10 = sp[cs9, <|"Sectors" -> {}, "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 5|>,
  "TWindow" -> <|"CompleteMax" -> 8|>|>, req[0, 5, 8]];
assert["su10_empty_source",
  p10["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 5|> &&
  AllTrue[Flatten[p10["Sectors"][[1]]["Coeffs"]], # === 0 &]];

(* SU-13: degenerate cleared denominator (t - eps class) *)
cs13 = catchDE2[Module[{csX = pc[mkSys[{{1/(x - eps)}}], mkChart[0, 1, "su13"]]},
  sc[csX, req[0, 1, 3]]]];
assert["su13_degenerate_denominator",
  FailureQ[cs13] && (cs13["ID"] === "E3" || cs13["ID"] === "E2" || cs13["ID"] === "E1")];

(* SU-14 analogue: non-Fuchsian refusal happens in Indicial (E3/E4 there);
   PrepareChart surfaces it *)
e14 = catchDE2[pc[mkSys[{{0, 1}, {1, 0}}/x^2], mkChart[0, 1, "su14"]]];
assert["su14_irregular_surfaced", FailureQ[e14] && e14["ID"] === "E3"];

(* 2F1 resonant chart (the committed matrices' residue, closed form):
   z A = {{0,0},{1/12,-2}}: basis exponents {0,-2}, TrueResonant *)
cs2f1 = pc[mkSys[{{0, 0}, {1/(12 x), -2/x}}], mkChart[0, 1, "su2f1"]];
r2f1 = catchDE2[sc[cs2f1, req[0, 2, 6]]];
assert["su2f1_true_resonant_basis",
  !FailureQ[r2f1] && Length[r2f1["Basis"]["Columns"]] === 2];

(* SU-18: gauge chart end-to-end (T-12 rank reduction, closed form
   f = (c1 - c2/x, c2)): the residual check must verify the GAUGED-BACK
   solution against the ORIGINAL system *)
cs18 = pc[mkSys[{{0, x^-2}, {0, 0}}], mkChart[0, 1, "su18"]];
r18 = catchDE2[sc[cs18, req[0, 2, 6]]];
assert["su18_gauge_chart_closed_form",
  !FailureQ[r18] && Length[r18["Basis"]["Columns"]] === 2 &&
  Module[{vals},
    vals = Map[Module[{ev = DiffExp2`SectorSeries`EvaluateLocalSolution[#, 1/3,
        "UsePade" -> False]["Value"]},
      N[DiffExp2`EpsSeries`ESCoefficient[ev, 0], 10]] &,
      r18["Basis"]["Columns"]];
    (* span of (1,0) and (-3,1) at x=1/3: check both columns lie in it *)
    AllTrue[vals, Abs[#[[1]] + 3*#[[2]] - (#[[1]] + 3*#[[2]])] < 10^-8 &] &&
    Module[{m = vals},
      Abs[Det[m]] > 10^-8 &&
      AllTrue[m, NumericQ[#[[1]]] && NumericQ[#[[2]]] &]]]];

(* SU-19: the canonical spectral lattice clears arbitrary poles from each
   Jordan block before recurrence.  The same depth moves into VInv, where
   source/matching transforms must still budget it honestly; homogeneous
   columns themselves now start finite instead of inheriting a normalization
   pole.  Gauge poles remain an independent final-transform budget below. *)
cs19v = pc[mkSys[{{1, 0}, {eps^-8, 0}}/x],
  mkChart[0, 1, "su19_spectral_eps_pole"]];
r19v = catchDE2[sc[cs19v, req[0, 2, 6]]];
assert["su19_spectral_primitive_lattice_budget",
  !FailureQ[r19v] &&
  cs19v["SpectralBlockEpsShifts"] === {8, 0} &&
  cs19v["V"] === {{eps^8, 0}, {1, 1}} &&
  DiffExp2`Solve`Private`spectralTransformPoleDepth[cs19v] === 0 &&
  DiffExp2`Solve`Private`inverseSpectralTransformPoleDepth[cs19v] === 8 &&
  DiffExp2`Solve`Private`finalTransformPoleDepth[cs19v, 6] === 0 &&
  (#["EpsWindow", "CompleteMax"] & /@ r19v["Basis"]["Columns"]) === {2, 2} &&
  Min[#["EpsWindow", "Min"] & /@ r19v["Basis"]["Columns"]] === 0];

cs19t = pc[mkSys[{{0, 1/(eps^8*x^2)}, {0, 0}}],
  mkChart[0, 1, "su19_gauge_eps_pole"]];
r19t = catchDE2[sc[cs19t, req[0, 2, 6]]];
assert["su19_gauge_transform_budget",
  !FailureQ[r19t] &&
  DiffExp2`Solve`Private`spectralTransformPoleDepth[cs19t] === 0 &&
  DiffExp2`Solve`Private`gaugeCoefficientData[cs19t, 6]["EpsValuation"] === -8 &&
  DiffExp2`Solve`Private`finalTransformPoleDepth[cs19t, 6] === 8 &&
  AllTrue[r19t["Basis"]["Columns"],
    #["EpsWindow"] === <|"Min" -> -8, "CompleteMax" -> 2|> &]];

(* SU-22: two successive CASE-P layers.  The a=0 column first collides
   with (1,eps), then its eps^-1 descendant collides with (2,2eps), so the
   raw depth reaches eps^-2.  Joint target-column compensation must cancel
   both polar orders as values. *)
B22 = DiagonalMatrix[{0, 1 + eps, 2 + 2 eps}] +
  t*{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
cs22 = pc[mkSys[(B22 /. t -> x)/x], mkChart[0, 1, "su22_depth_two"]];
r22 = catchDE2[sc[cs22, req[-3, 3, 4]]];
col22 = If[FailureQ[r22], None,
  SelectFirst[r22["Basis", "Columns"], MemberQ[tagsOf[#], {0, 0, 0}] &]];
ev22 = If[col22 === None, None,
  DiffExp2`SectorSeries`EvaluateLocalSolution[col22, 1/3,
    "UsePade" -> False]["Value"]];
assert["su22_depth_two_pseudo_compensation",
  !FailureQ[r22] && col22 =!= None &&
  r22["Basis", "Diagnostics", "PseudoCollisionsCompensated"] === True &&
  col22["EpsWindow", "Min"] === -2 &&
  AllTrue[Flatten[Table[
    DiffExp2`EpsSeries`ESCoefficient[ev22, k], {k, -2, -1}]],
    PossibleZeroQ]];

(* SU-23: an epsilon-singular spectral inverse does not make its individual
   coordinates a valid CASE-P regularity certificate.  This is the exact
   three-root structure of the banana L2 x=0 endpoint: the compensated
   physical collision direction is regular, while VInv maps its finite
   eps^0 value back to a legitimate eps^-1 spectral coordinate. *)
r23c = -264/529 + (396/529)*eps;
M023 = {{-1 + eps, 0, 0}, {-2, -1 + 3 eps, 0},
  {0, 0, -2 + 2 eps}};
cs23 = pc[mkSys[M023/x + {{0, 1, 0}, {0, 0, r23c}, {0, 0, 0}}],
  mkChart[0, 1, "su23_singular_vinv_certificate"]];
r23 = catchDE2[sc[cs23, req[0, 3, 8]]];
col23 = If[FailureQ[r23], None, SelectFirst[r23["Basis", "Columns"],
  MemberQ[tagsOf[#], {-2, 2, 0}] &]];
oldFrame23 = If[col23 === None, None,
  DiffExp2`Solve`Private`spectralProbeValue[cs23, col23, 1/3]];
assert["su23_collision_certificate_uses_reduced_physical_frame",
  !FailureQ[r23] && col23 =!= None &&
  r23["Basis", "Diagnostics", "PseudoCollisionsCompensated"] === True &&
  cs23["VInv"] === {{-1/eps, 1, 0}, {1/eps, 0, 0}, {0, 0, 1}} &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[oldFrame23[[1]], -1], 30]] > 10^-6];

(* SU-25: a coalescing pair can make Indicial's deterministic first-entry
   normalization meromorphic even though the exact system is epsilon-
   regular.  Here the raw source eigenvector is (1,1/eps,0), while the
   n=1 CASE-P target is (0,1,1).  Without a blockwise primitive scaling the
   collision's eps^-1 certificate sees the source seed's unrelated t/eps in
   reduced component 2 and fails.  Multiplying the whole one-vector Jordan
   block by eps gives the equally valid seed (eps,1,0); the transformed
   coupling also gains eps, so the exact quotient is regular rather than
   hiding a pole behind a relaxed certificate. *)
M025 = {{1 - 2 eps, 0, 0}, {-1, 1 - eps, 1}, {0, 0, 2 - eps}};
N125 = {{0, 0, 0}, {1, 0, 0}, {1, 0, 0}};
chain25 = {{1, 1/eps}, {0, 1/eps}};
primitive25 = DiffExp2`Solve`Private`epsilonPrimitiveJordanChain[
  chain25, eps];
lambda25 = 1 - 2 eps;
Vchain25 = Transpose[chain25];
Rchain25 = Map[Cancel[Together[#]] &,
  Vchain25.{{lambda25, 1}, {0, lambda25}}.Inverse[Vchain25], {2}];
assert["su25_primitive_scaling_is_common_over_jordan_chain",
  primitive25 === {{{eps, 1}, {0, 1}}, 1} &&
  Map[Cancel[Together[#]] &,
    Rchain25.Transpose[primitive25[[1]]] -
      Transpose[primitive25[[1]]].{{lambda25, 1}, {0, lambda25}},
    {2}] === ConstantArray[0, {2, 2}]];
cs25 = pc[mkSys[(M025 + x*N125)/x],
  mkChart[0, 1, "su25_primitive_casep_frame"]];
r25 = catchDE2[sc[cs25, req[0, 2, 4]]];
col25 = If[FailureQ[r25], None, SelectFirst[r25["Basis", "Columns"],
  MemberQ[tagsOf[#], {1, -2, 0}] &]];
assert["su25_casep_uses_epsilon_primitive_jordan_blocks",
  !FailureQ[r25] && col25 =!= None &&
  cs25["SpectralBlockEpsShifts"] === {0, 0, 1} &&
  cs25["V"] === {{0, 0, eps}, {1, 1, 1}, {1, 0, 0}} &&
  r25["Basis", "Diagnostics", "PseudoCollisionsCompensated"] === True &&
  r25["Basis", "Diagnostics", "PseudoCollisionsHit"] === {
    <|"n" -> 1, "Cols" -> {1}, "DeltaB" -> -1,
      "PolarOrders" -> {}|>} &&
  col25["EpsWindow", "Min"] === 0 &&
  Module[{sec = First[col25["Sectors"]]},
    rowOf[col25, sec, 0][[1 ;; 2]] ===
      {{0, 1, 0}, {0, -1, -1}} &&
    rowOf[col25, sec, 1][[1]] === {1, 0, 0}]];

(* SU-20: a finite V can have an eps-singular inverse (det V ~ eps^8).
   Particular-source framing crosses VInv before the recurrence and must
   budget that depth independently of the final V transform. *)
cs20 = pc[mkSys[{{0, eps^-8}, {0, 1}}/x],
  mkChart[0, 1, "su20_particular_vinv_eps_pole"]];
src20 = <|"Sectors" -> {<|"a" -> 1/2, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {1, 0}, {0, 0}],
      {k, 0, 10}, {n, 0, 6}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 10|>,
  "TWindow" -> <|"CompleteMax" -> 6|>|>;
p20 = catchDE2[sp[cs20, src20, req[0, 2, 6]]];
assert["su20_particular_inverse_spectral_budget",
  DiffExp2`Solve`Private`spectralTransformPoleDepth[cs20] === 0 &&
  DiffExp2`Solve`Private`inverseSpectralTransformPoleDepth[cs20] === 8 &&
  !FailureQ[p20] && p20["EpsWindow", "CompleteMax"] === 2];

(* SU-21: gauge coefficient extraction and convolution must use the same
   m=gv..gv+nmax range.  The old gv..nmax range was empty for gv>nmax and
   applyGauge then indexed past it. *)
gd21 = DiffExp2`Solve`Private`gaugeCoefficientData[
  <|"ChartVar" -> t, "Gauge" -> t^2 IdentityMatrix[1],
    "SystemSize" -> 1|>, 1];
assert["su21_positive_gauge_valuation_range",
  gd21["TMin"] === 2 && gd21["EpsValuation"] === 0 &&
  Length[gd21["Coefficients"]] === 2 &&
  gd21["Coefficients"][[1]] === {{1}} &&
  gd21["Coefficients"][[2]] === {{0}}];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
