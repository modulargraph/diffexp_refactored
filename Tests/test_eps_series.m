(* Focused unit tests for DiffExp2/EpsSeries.m. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
failsWith[f_, id_String] := FailureQ[f] && f["ID"] === id;

(* tolerances live: wp 300 per tests 11/14 *)
catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 300}]];

es = DiffExp2`EpsSeries`ESNew;
zero = DiffExp2`EpsSeries`ESZero;
coeff = DiffExp2`EpsSeries`ESCoefficient;
clist = DiffExp2`EpsSeries`ESCoefficientList;
lead = DiffExp2`EpsSeries`ESLeading;
add = DiffExp2`EpsSeries`ESAdd;
times = DiffExp2`EpsSeries`ESTimes;
inv = DiffExp2`EpsSeries`ESInvert;
div = DiffExp2`EpsSeries`ESDivide;
shift = DiffExp2`EpsSeries`ESShift;
trunc = DiffExp2`EpsSeries`ESTruncate;
trim = DiffExp2`EpsSeries`ESTrim;
trimThrough = DiffExp2`EpsSeries`ESTrimThrough;
sameq = DiffExp2`EpsSeries`ESSameQ;
fromx = DiffExp2`EpsSeries`ESFromExpression;
win = DiffExp2`EpsSeries`ESWindow;
esq = DiffExp2`EpsSeries`ESQ;

(* 1 *)
s1 = es[-2, {1, 0, 3/7}];
assert["construct_validate_ok",
  esq[s1] && win[s1] === <|"Min" -> -2, "CompleteMax" -> 0|> &&
  coeff[s1, -2] === 1 && coeff[s1, 0] === 3/7];

(* 2 *)
assert["construct_validate_errors",
  failsWith[catchDE2[es[1/2, {1}]], "ERR-BAD-CONSTRUCT"] &&
  failsWith[catchDE2[es[0, {}]], "ERR-BAD-CONSTRUCT"] &&
  failsWith[catchDE2[es[0, {Global`eps}]], "ERR-BAD-CONSTRUCT"] &&
  failsWith[catchDE2[es[0, {$Failed}]], "ERR-BAD-CONSTRUCT"] &&
  failsWith[catchDE2[es[0, {Indeterminate}]], "ERR-BAD-CONSTRUCT"]];

(* 3 *)
s3 = es[-1, {2, 5}];
e3 = Block[{DiffExp2`EpsSeries`$ESErrorContext = "ctx-three"},
  catchDE2[coeff[s3, 1]]];
assert["coefficient_access_semantics",
  coeff[s3, -3] === 0 && coeff[s3, 0] === 5 &&
  failsWith[e3, "ERR-WINDOW-READ"] && e3["Requested"] === 1 &&
  e3["Min"] === -1 && e3["CompleteMax"] === 0 && e3["Context"] === "ctx-three"];

(* 4 *)
s4 = add[es[0, {1, 1, 1, 1, 1, 1}], es[-1, {7, 2, 3, 4, 9}]];
assert["add_min_window_no_union",
  win[s4] === <|"Min" -> -1, "CompleteMax" -> 3|> &&
  s4["Coeffs"] === {7, 3, 4, 5, 10}];

(* 5 *)
s5 = times[es[0, {1, 1, 1}], es[1, {2, 0, 1}]];
assert["mul_window_convolution",
  win[s5] === <|"Min" -> 1, "CompleteMax" -> 3|> &&
  s5["Coeffs"] === {2, 2, 3}];

(* 6 — THE math-review 6(iv) pin *)
one6 = fromx[1, Global`eps, 5];
d6 = fromx[3 Global`eps, Global`eps, 6];
q6 = div[one6, d6];
assert["div_by_beps_shift",
  win[one6] === <|"Min" -> 0, "CompleteMax" -> 5|> &&
  win[d6] === <|"Min" -> 1, "CompleteMax" -> 6|> &&
  win[q6] === <|"Min" -> -1, "CompleteMax" -> 4|> &&
  coeff[q6, -1] === 1/3 && Head[coeff[q6, -1]] === Rational &&
  AllTrue[Range[0, 4], coeff[q6, #] === 0 &]];

(* 7 *)
q7 = div[fromx[1, Global`eps, 4], fromx[1 - Global`eps, Global`eps, 4]];
assert["div_geometric_exact",
  win[q7] === <|"Min" -> 0, "CompleteMax" -> 4|> &&
  q7["Coeffs"] === {1, 1, 1, 1, 1}];

(* 8 *)
q8 = inv[es[1, {5/3, 0, 0, 0}]];
assert["invert_beps_trailing_zeros",
  win[q8] === <|"Min" -> -1, "CompleteMax" -> 2|> &&
  q8["Coeffs"] === {3/5, 0, 0, 0}];

(* 9 *)
q9 = div[fromx[1 + Global`eps, Global`eps, 2], es[1, {Global`b1 - Global`b2, 0, 0}]];
assert["pseudo_resonance_symbolic_shift",
  win[q9] === <|"Min" -> -1, "CompleteMax" -> 1|> &&
  PossibleZeroQ[coeff[q9, -1] - 1/(Global`b1 - Global`b2)]];

(* 10 *)
e10 = catchDE2[inv[zero[5]]];
assert["invert_zero_denominator_error",
  failsWith[e10, "ERR-DIV-ZERO"] && e10["Min"] === 5 && e10["CompleteMax"] === 5];

(* 11 — relative-not-absolute lead (the 10^-40 regression pin) *)
d11a = es[0, {SetPrecision[10^-60, 300], 2, 1}];
d11b = es[0, {SetPrecision[10^-60, 300], SetPrecision[10^-58, 300]}];
assert["relative_not_absolute_lead",
  lead[d11a][[1]] === 1 && lead[d11a][[2]] === 2 &&
  lead[d11b][[1]] === 0];

(* A low-accuracy centered-zero imaginary component must not collapse the
   trailing coefficient that defines the relative series scale. *)
complexScale11 =
  4.4267459561002104836`0.2683567342206439 +
    0``-0.022790862816694443*I;
d11c = es[0, {10^-30, complexScale11}];
assert["complex_scale_does_not_poison_leading_order",
  lead[d11c][[1]] === 1 && lead[d11c][[2]] === complexScale11];
diagonal11 = es[0, {(8 + 8 I)*10^-25, 1}];
assert["true_complex_modulus_preserves_laurent_lead",
  DiffExp2`EpsSeries`ESMinPower[trim[diagonal11]] === 0];

(* 12 *)
assert["trim_preserves_zero_window",
  trim[es[-2, {0, 0, 0, 0}]] === zero[1]];
resolvedCenteredZero = SetPrecision[0., 250];
assert["trim_resolved_all_centered_zero_uses_unit_decision_scale",
  trim[es[-54, ConstantArray[resolvedCenteredZero, 13]]] === zero[-42]];
unresolvedCenteredZero = 0``17;
assert["trim_unresolved_centered_zero_is_not_discarded",
  DiffExp2`EpsSeries`ESMinPower[
    trim[es[-1, {unresolvedCenteredZero}]]] === -1];
privateReservoir12 = es[-1, {N[10^-8, 100], N[2, 100],
    N[3, 100], N[10^80, 100]}];
assert["private reservoir cannot change the public Laurent lead",
  DiffExp2`EpsSeries`ESMinPower[trim[privateReservoir12]] === 2 &&
    DiffExp2`EpsSeries`ESMinPower[
      trimThrough[privateReservoir12, 1]] === -1 &&
    coeff[trimThrough[privateReservoir12, 1], -1] ===
      coeff[privateReservoir12, -1] &&
    DiffExp2`EpsSeries`ESCompleteMax[
      trimThrough[privateReservoir12, 1]] === 2];
assert["all-zero public prefix advances without discarding private suffix",
  win[trimThrough[es[-2, {0, 0, 7, 9}], -1]] ===
      <|"Min" -> 0, "CompleteMax" -> 1|> &&
    trimThrough[es[-2, {0, 0, 7, 9}], -1]["Coeffs"] === {7, 9}];

(* 13 — exact passthrough chain; hand-computed pin *)
q13 = div[times[add[es[0, {1/3, 22/7}], es[0, {-5/11, 2}]], es[0, {1, 1}]], es[0, {2, 1}]];
assert["exact_rational_passthrough",
  AllTrue[q13["Coeffs"], MemberQ[{Integer, Rational}, Head[#]] &] &&
  coeff[q13, 0] === -2/33 && coeff[q13, 1] === 587/231];

(* 14 — 300-digit chain vs exact *)
n14 = {N[1/3, 300], N[22/7, 300], N[-5/11, 300], N[2, 300]};
qn = div[times[add[es[0, n14[[{1, 2}]]], es[0, n14[[{3, 4}]]]], es[0, {N[1, 300], N[1, 300]}]],
  es[0, {N[2, 300], N[1, 300]}]];
assert["precision_300_digits",
  AllTrue[Select[qn["Coeffs"], InexactNumberQ], Precision[#] >= 290 &] &&
  sameq[qn, q13]];

(* 15 *)
s15 = es[0, {1, 2}];
assert["shift_moves_both_edges",
  win[shift[s15, -3]] === <|"Min" -> -3, "CompleteMax" -> -2|> &&
  shift[s15, -3]["Coeffs"] === {1, 2} &&
  shift[shift[s15, 3], -3] === s15];

(* 16 *)
assert["truncate_no_pad",
  trunc[es[0, {1, 2, 3}], 1] === es[0, {1, 2}] &&
  failsWith[catchDE2[trunc[es[0, {1}], 5]], "ERR-TRUNCATE-EXTEND"] &&
  failsWith[catchDE2[trunc[es[0, {1, 2}], -1]], "ERR-TRUNCATE-EXTEND"]];

(* 17 *)
s17 = fromx[(1 + Global`eps)/(Global`eps^2 (1 - Global`eps)), Global`eps, 3];
assert["from_expression_exact_laurent",
  win[s17] === <|"Min" -> -2, "CompleteMax" -> 3|> &&
  s17["Coeffs"] === {1, 2, 2, 2, 2, 2}];
s17uncertain = fromx[0``17, Global`eps, 3];
assert["from_expression_inexact_zero_is_not_structural_zero",
  win[s17uncertain] === <|"Min" -> 0, "CompleteMax" -> 3|> &&
  InexactNumberQ[coeff[s17uncertain, 0]]];

(* 18 *)
s18 = fromx[Gamma[1 + Global`eps], Global`eps, 2];
assert["from_expression_gamma_prefactor",
  win[s18] === <|"Min" -> 0, "CompleteMax" -> 2|> &&
  coeff[s18, 0] === 1 &&
  Simplify[coeff[s18, 1] + EulerGamma] === 0 &&
  Simplify[coeff[s18, 2] - (EulerGamma^2/2 + Pi^2/12)] === 0];

(* 19 *)
assert["from_expression_loud_failures",
  failsWith[catchDE2[fromx[Exp[1/Global`eps], Global`eps, 3]], "ERR-EXPAND-FAIL"] &&
  Module[{f = catchDE2[fromx[Sqrt[Global`eps], Global`eps, 3]]},
    failsWith[f, "ERR-EXPAND-FAIL"] && f["LeadingPower"] === 1/2]];

(* 20 *)
s20 = es[-1, {4, 5, 6}];
t20 = trim[es[-1, {0, 5, 6}]];
assert["coefficient_list_contract",
  clist[s20, -2, 1] === {0, 4, 5, 6} &&
  failsWith[catchDE2[clist[s20, 0, 2]], "ERR-RANGE"] &&
  failsWith[catchDE2[clist[s20, 0, 1]], "ERR-DROP-BELOW"] &&
  clist[t20, 0, 1] === {5, 6}];
assert["polar_slice_idiom",
  clist[es[0, {7, 8}], 0, -1] === {} &&
  clist[s20, -1, -1] === {4}];

(* 21 *)
a21 = es[0, {N[1, 300], N[2, 300], N[3, 300]}];
b21 = es[0, {N[1, 300], N[2, 300]}];
b21p = es[0, {N[1, 300] + 10*10^-150*2, N[2, 300]}];
assert["same_q_tolerance_and_windows",
  sameq[a21, b21] &&
  !sameq[a21, b21p] &&
  sameq[es[1, {Global`b1 - Global`b2}], es[1, {-(Global`b2 - Global`b1)}]]];

(* 22 *)
e22 = Block[{DiffExp2`EpsSeries`$ESErrorContext =
    "chart x0=1, sector (a=-1,b=2,p=0), order k=7"},
  catchDE2[coeff[es[0, {1}], 3]]];
assert["error_context_propagation",
  e22["Context"] === "chart x0=1, sector (a=-1,b=2,p=0), order k=7"];

(* 23 *)
assert["window_object_shape",
  win[es[2, {1, 2}]] === <|"Min" -> 2, "CompleteMax" -> 3|>];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
