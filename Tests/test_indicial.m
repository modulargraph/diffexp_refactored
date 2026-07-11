(* Focused unit tests for DiffExp2/Indicial.m. Optional campaign-matrix pins
   are skipped loudly when their development fixtures are absent. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
failsWith[f_, id_String] := FailureQ[f] && f["ID"] === id && f["Module"] === "Indicial";

catchDE2[DiffExp2`Config`LoadConfiguration[{}]];

ci = DiffExp2`Indicial`ChartIndicial;
ref = <|"Name" -> "test.chart", "Center" -> 0, "Variable" -> Global`x|>;
t = Global`t; eps = Global`eps;
specOf[fam_] := {#["a"], #["b"], #["p"]} & /@ fam["Sectors"];

(* T-1 *)
r1 = ci[{{0, 1}, {1/(1 - t), 0}}, t, eps, ref];
assert["test_regular_chart",
  r1["Regular"] === True && Length[r1["Families"]] === 1 &&
  r1["Families"][[1]]["Sectors"] === ConstantArray[<|"a" -> 0, "b" -> 0, "p" -> 0|>, 2] &&
  r1["Reduction"]["Gauge"] === IdentityMatrix[2] &&
  r1["PoleData"]["PoleOrder"] === 0];

(* T-2 *)
r2 = catchDE2[ci[{{(2 t)/(2 t^2) - 1/t, 1}, {0, 0}}, t, eps, ref]];
assert["test_zero_residue_after_cancellation",
  !FailureQ[r2] && r2["Regular"] === True];

(* T-3 *)
r3 = ci[{{(1 + eps)/t, 1}, {0, (3 - 2 eps)/(2 t)}}, t, eps, ref];
assert["test_simple_pole_two_families",
  r3["PoleData"]["PoleOrder"] === 1 &&
  Sort[{#["a"], #["b"], #["Multiplicity"]} & /@ r3["Spectrum"]] ===
    Sort[{{1, 1, 1}, {3/2, -1, 1}}] &&
  Length[r3["Families"]] === 2 &&
  AllTrue[r3["Families"], #["Class"] === "Single" && #["JointSolve"] === False &],
  AllTrue[Flatten[specOf /@ r3["Families"], 1], #[[3] ] === 0 &]];

(* T-4 *)
e4 = catchDE2[ci[{{0, 1/t}, {eps/t, 0}}, t, eps, ref]];
assert["test_i1_violation_sqrt_eps",
  failsWith[e4, "E2"] && e4["Chart"] === "test.chart"];

(* T-5 *)
e5 = catchDE2[ci[{{eps/((1 + eps) t), 0}, {0, 0}}, t, eps, ref]];
assert["test_i1_violation_rational_root", failsWith[e5, "E2"]];

(* T-6 *)
r6 = ci[{{2 eps/t, 1/t}, {-eps^2/t, 0}}, t, eps, ref];
assert["test_affine_pass_nonlinear_entries",
  {#["a"], #["b"], #["Multiplicity"]} & /@ r6["Spectrum"] === {{0, 1, 2}} &&
  r6["Spectrum"][[1]]["BlockSizes"] === {2} &&
  specOf[r6["Families"][[1]]] === {{0, 1, 0}, {0, 1, 1}} &&
  r6["Families"][[1]]["Class"] === "Confluent" &&
  r6["Families"][[1]]["LogMax"] === 1];

(* T-7 *)
r7a = ci[{{eps/t, 0}, {0, eps/t}}, t, eps, ref];
r7b = ci[{{eps/t, 1/t}, {0, eps/t}}, t, eps, ref];
assert["test_jordan_vs_diagonalizable",
  r7a["Spectrum"][[1]]["BlockSizes"] === {1, 1} &&
  specOf[r7a["Families"][[1]]] === {{0, 1, 0}, {0, 1, 0}} &&
  r7a["Families"][[1]]["LogMax"] === 0 &&
  r7b["Spectrum"][[1]]["BlockSizes"] === {2} &&
  specOf[r7b["Families"][[1]]] === {{0, 1, 0}, {0, 1, 1}} &&
  r7b["Families"][[1]]["LogMax"] === 1];

(* T-8: 2F1 resonant residue {{0,0},{1/12,-2}} *)
r8 = ci[{{0, 0}, {1/(12 t), -2/t}}, t, eps, ref];
f8 = r8["Families"][[1]];
assert["test_true_resonance_2f1_resonant",
  Length[r8["Families"]] === 1 && f8["Class"] === "TrueResonant" &&
  f8["JointSolve"] === False &&
  specOf[f8] === {{-2, 0, 1}, {0, 0, 0}} &&
  f8["LogMax"] === 1 &&
  Length[f8["Collisions"]] === 1 &&
  f8["Collisions"][[1]]["n"] === 2 &&
  f8["Collisions"][[1]]["Type"] === "Log" &&
  f8["Collisions"][[1]]["DeltaB"] === 0];

(* T-9: 2F1 nonresonant residue {{0,0},{1/12,-3/2}} *)
r9 = ci[{{0, 0}, {1/(12 t), -3/(2 t)}}, t, eps, ref];
assert["test_nonresonant_2f1",
  Length[r9["Families"]] === 2 &&
  AllTrue[r9["Families"], #["Collisions"] === {} &] &&
  AllTrue[Flatten[specOf /@ r9["Families"], 1], #[[3]] === 0 &]];

(* T-10: banana pseudo-resonance mixture diag(0, -1+eps, 2eps) *)
r10 = ci[DiagonalMatrix[{0, -1 + eps, 2 eps}]/t, t, eps, ref];
f10 = r10["Families"][[1]];
colKeys[c_] := {c["n"], c["DeltaB"], c["Type"]};
assert["test_pseudo_resonance_banana_mixture",
  Length[r10["Families"]] === 1 && f10["Class"] === "Pseudo" &&
  f10["JointSolve"] === True && f10["LogMax"] === 0 &&
  AllTrue[f10["Sectors"], #["p"] === 0 &] &&
  Sort[colKeys /@ f10["Collisions"]] ===
    Sort[{{1, 1, "LaurentShift"}, {1, -1, "LaurentShift"}, {0, -2, "LaurentShift"}}]];

(* T-11: same a, different b *)
r11 = ci[DiagonalMatrix[{eps, 2 eps}]/t, t, eps, ref];
f11 = r11["Families"][[1]];
assert["test_pseudo_resonance_same_a",
  f11["JointSolve"] === True &&
  Length[f11["Collisions"]] === 1 &&
  colKeys[f11["Collisions"][[1]]] === {0, -1, "LaurentShift"} &&
  AllTrue[f11["Sectors"], #["p"] === 0 &]];

(* T-12: rank reduction closed form A = {{0, t^-2}, {0, 0}} *)
r12 = ci[{{0, t^-2}, {0, 0}}, t, eps, ref];
red12 = r12["Reduction"];
assert["test_rank_reduction_closed_form",
  r12["PoleData"]["PoleOrder"] === 2 &&
  red12["Steps"] >= 1 &&
  Sort[{#["a"], #["b"]} & /@ r12["Spectrum"]] === {{0, 0}, {1, 0}} &&
  Length[r12["Families"]] === 1 &&
  r12["Families"][[1]]["Class"] === "TrueResonant"];
(* composition check: B is constant for this input; closed-form reduced basis
   g with theta g = B g via matrix power t^B: verify f = T.g solves f' = A f *)
Module[{T12 = red12["Gauge"], B12 = red12["ThetaMatrix"], g1, g2, ok},
  ok = matZero = AllTrue[Flatten[#], PossibleZeroQ[Together[#]] &] &;
  (* B constant: residue == ThetaMatrix *)
  If[AllTrue[Flatten[B12 - red12["Residue"]], PossibleZeroQ[Together[#]] &],
    Module[{R12 = red12["Residue"], sols},
      (* solutions g = t^R . c for constant theta-matrix *)
      sols = Transpose[MatrixExp[R12*Log[t]]];
      assert["test_rank_reduction_composition",
        AllTrue[sols, Module[{f = T12 . #},
          AllTrue[Flatten[{D[f, t] - {{0, t^-2}, {0, 0}} . f}],
            PossibleZeroQ[Together[FullSimplify[#]]] &]] &]]],
    assert["test_rank_reduction_composition", False]]];

(* T-13: irregular loud *)
e13 = catchDE2[ci[{{0, 1/t^2}, {1/t^2, 0}}, t, eps, ref]];
assert["test_irregular_loud",
  failsWith[e13, "E3"] && e13["PoleOrder"] === 2];

(* T-14: ramified-irregular non-termination *)
e14 = catchDE2[ci[{{0, t^-2}, {t^-1, 0}}, t, eps, ref]];
assert["test_nontermination_loud",
  failsWith[e14, "E4"] && e14["MaxSteps"] === 200];

(* T-18: constructed diagonalizable spectrum roundtrip + determinism *)
S18 = {{1, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}};
R18 = S18 . DiagonalMatrix[{1 + 2 eps, -1/2, 3 eps, 1 + 2 eps}] . Inverse[S18];
r18 = ci[R18/t, t, eps, ref];
r18b = ci[R18/t, t, eps, ref];
assert["test_constructed_spectrum_roundtrip",
  Sort[{#["a"], #["b"], #["Multiplicity"]} & /@ r18["Spectrum"]] ===
    Sort[{{1, 2, 2}, {-1/2, 0, 1}, {0, 3, 1}}] &&
  Module[{rec = SelectFirst[r18["Spectrum"], #["a"] === 1 &]},
    rec["BlockSizes"] === {1, 1}] &&
  Length[r18["Families"]] === 2 &&
  Module[{pf = SelectFirst[r18["Families"], #["Class"] === "Pseudo" &]},
    pf =!= Missing["NotFound"] && pf["JointSolve"] === True &&
    AllTrue[pf["Sectors"], #["p"] === 0 &] &&
    Length[pf["Collisions"]] === 1 &&
    colKeys[pf["Collisions"][[1]]] === {1, 1, "LaurentShift"}] &&
  Module[{sf = SelectFirst[r18["Families"], #["Class"] === "Single" &]},
    sf =!= Missing["NotFound"] &&
    sf["Members"][[1]]["a"] === -1/2] &&
  r18 === r18b];

(* T-19: inexact input loud *)
e19 = catchDE2[ci[{{0.5/t, 0}, {0, 0}}, t, eps, ref]];
assert["test_inexact_input_loud", failsWith[e19, "E1"]];

(* T-20: SeriesData input loud with full-export hint *)
e20 = catchDE2[ci[{{Series[1/(1 - eps), {eps, 0, 2}], 0}, {0, 0}}/t, t, eps, ref]];
assert["test_series_input_loud",
  failsWith[e20, "E1"] && StringContainsQ[e20["Detail"], "full export"]];

(* T-21: no numerics in source *)
src21 = ReadString[FileNameJoin[{repoRoot, "DiffExp2", "Indicial.m"}]];
assert["test_no_numerics_in_source",
  !StringContainsQ[src21, "N["] &&
  !StringContainsQ[src21, "Rationalize["] &&
  !StringContainsQ[src21, "Chop["] &&
  !StringContainsQ[src21, "Tolerance ->"] &&
  !StringContainsQ[src21, "Quiet["] &&
  !StringContainsQ[src21, "Check["]];

(* T-22: payload completeness for triggered errors *)
assert["test_error_payload_completeness",
  AllTrue[{e4, e13, e14, e19, e20},
    #["Module"] === "Indicial" && #["Chart"] === "test.chart" &&
    !MatchQ[#["Center"], _Missing] && !MatchQ[#["Variable"], _Missing] &]];

(* T-23: eps-zero degeneracy *)
r23a = ci[DiagonalMatrix[{eps, 2 eps}]/t, t, eps, ref];
r23b = ci[{{0, 1/t}, {0, 2 eps/t}}, t, eps, ref];
assert["test_eps_zero_degeneracy",
  r23a["Families"][[1]]["EpsZeroDegeneracy"] === 0 &&
  DiffExp2`Indicial`EpsDegenerateFamilies[r23a] === {} &&
  r23b["Families"][[1]]["EpsZeroDegeneracy"] === 1 &&
  DiffExp2`Indicial`EpsDegenerateFamilies[r23b] ===
    {<|"FamilyIndex" -> 1, "EpsZeroDegeneracy" -> 1|>}];

(* T-15..T-17 placeholders: vendored campaign matrices *)
If[FileExistsQ[FileNameJoin[{repoRoot, "Tests", "refs", "banana_l1_dxx1_full.m"}]],
  Print["  (T-15/T-16 banana matrices found - implement pins)"],
  Print["  SKIP: T-15..T-17 (optional development campaign matrices are not shipped)"]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
