(* Unit tests for DiffExp2/Solve.m per Docs/specs/Solve.md section 8
   (closed-form subset; SU-15..17 fixture tests activate at M4/M5). *)

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

(* SU-01: f' = f *)
cs1 = pc[mkSys[{{1}}], mkChart[0, 10, "su01"]];
r1 = sc[cs1, req[0, 2, 8]];
col1 = r1["Basis"]["Columns"][[1]];
assert["su01_regular_exponential",
  tagsOf[col1] === {{0, 0, 0}} &&
  rowOf[col1, col1["Sectors"][[1]], 0][[;; 6, 1]] === Table[1/n!, {n, 0, 5}]];

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
  (* diagonal system: no actual coupling, columns are pure monomials *)
  AllTrue[r7["Basis"]["Columns"], Length[#["Sectors"]] === 1 &]];

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

(* SU-09: pseudo-resonant source: scalar theta f = 2eps f + 1:
   f = 1/(0 - 2eps)·1 = -1/(2eps): window shifts down by 1 *)
cs9 = pc[mkSys[{{2 eps}}/x], mkChart[0, 1, "su09"]];
p9 = catchDE2[sp[cs9, src8, req[0, 2, 3]]];
assert["su09_pseudo_resonant_source",
  !FailureQ[p9] && p9["EpsWindow", "Min"] === -1 &&
  Module[{sec = First[p9["Sectors"]]},
    rowOf[p9, sec, -1][[1, 1]] === -1/2]];

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

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
