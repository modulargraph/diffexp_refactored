(* Unit tests for DiffExp2/Integrate.m: the case-table closed forms. *)
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
smi = DiffExp2`Integrate`SectorMonomialIntegral;
esC = DiffExp2`EpsSeries`ESCoefficient;
esMn = DiffExp2`EpsSeries`ESMinPower;
eqQ[x_, y_] := TrueQ[PossibleZeroQ[Simplify[x - y]]];

(* cell 1: Int_0^T t^m = T^(m+1)/(m+1) *)
r1 = smi[2, 0, 0, 1/2, 3];
assert["cell1", eqQ[esC[r1, 0], (1/2)^3/3] && esC[r1, 1] === 0];

(* cell 2: Int_0^T t (eps Log t) = eps T^2 (Log T/2 - 1/4) *)
r2 = smi[1, 0, 1, 1/3, 3];
assert["cell2", esMn[r2] === 1 &&
  eqQ[esC[r2, 1], (1/3)^2*(Log[1/3]/2 - 1/4)]];

(* cell 7: Int_0^T t^(b eps) = T^(1 + b eps)/(1 + b eps): orders 0..2 *)
r7 = smi[0, 2, 0, 1/2, 2];
ex7 = Normal[Series[(1/2)^(1 + 2 Global`e)/(1 + 2 Global`e), {Global`e, 0, 2}]];
assert["cell7", AllTrue[Range[0, 2],
  eqQ[esC[r7, #], SeriesCoefficient[ex7, {Global`e, 0, #}]] &]];

(* cell 9: Int_0^T t^(-1 + b eps) = T^(b eps)/(b eps): Laurent shift by 1 *)
r9 = smi[-1, 3, 0, 1/4, 2];
ex9 = Normal[Series[(1/4)^(3 Global`e)/(3 Global`e), {Global`e, 0, 2}]];
assert["cell9", esMn[r9] === -1 &&
  AllTrue[Range[-1, 2],
    eqQ[esC[r9, #], SeriesCoefficient[ex9, {Global`e, 0, #}]] &]];

(* cell 10: p=1 pole depth 2 *)
r10 = smi[-1, 2, 1, 1/2, 1];
ex10 = Normal[Series[Integrate[Global`u^(-1 + 2 Global`e)*Global`e*Log[Global`u],
  {Global`u, 0, 1/2}, Assumptions -> Global`e > 0], {Global`e, 0, 1}]];
assert["cell10", esMn[r10] === -1 &&
  AllTrue[Range[-1, 1],
    eqQ[esC[r10, #], SeriesCoefficient[ex10, {Global`e, 0, #}]] &]];

(* cell 11: analytic regularization: Int_0^T t^(-2 + b eps): NO error *)
r11 = smi[-2, 1, 0, 1/2, 1];
ex11 = Normal[Series[(1/2)^(-1 + Global`e)/(-1 + Global`e), {Global`e, 0, 1}]];
assert["cell11", AllTrue[Range[0, 1],
  eqQ[esC[r11, #], SeriesCoefficient[ex11, {Global`e, 0, #}]] &]];

(* divergent cell 3: loud *)
e3 = catchDE2[smi[-1, 0, 0, 1/2, 1]];
assert["cell3_loud", FailureQ[e3] && e3["ID"] === "E2"];

(* IntegrateLocalSolution: f = 1 + t over [0, 1/2]: 1/2 + 1/8 *)
mk1[secs_, kmin_, kmax_, ncols_] := <|"Center" -> 0,
  "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>, "Radius" -> 1,
  "Sectors" -> secs, "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
  "TWindow" -> <|"CompleteMax" -> ncols - 1|>,
  "ErrorEstimate" -> ConstantArray[0, kmax - kmin + 1], "Prescriptions" -> {}|>;
ls1 = mk1[{<|"a" -> 0, "b" -> 0, "p" -> 0,
  "Coeffs" -> {{{1}, {1}, {0}}}|>}, 0, 0, 3];
ri1 = DiffExp2`Integrate`IntegrateLocalSolution[ls1, {0, 1/2}];
assert["ils_regular", eqQ[esC[ri1["Values"][[1]], 0], 1/2 + 1/8]];

(* interior PV: f = 1/t over [-1/4, 1/2]: PV = Log[(1/2)/(1/4)] = Log 2 *)
ls2 = mk1[{<|"a" -> -1, "b" -> 0, "p" -> 0,
  "Coeffs" -> {{{1}, {0}, {0}}}|>}, 0, 0, 3];
ri2 = DiffExp2`Integrate`IntegrateLocalSolution[ls2, {-1/4, 1/2}];
assert["ils_interior_pv", eqQ[esC[ri2["Values"][[1]], 0], Log[2]]];

(* EndpointSectorLimit: drop rule *)
ls3 = mk1[{<|"a" -> 0, "b" -> 0, "p" -> 0, "Coeffs" -> {{{7}, {1}, {0}}}|>,
  <|"a" -> -1, "b" -> 2, "p" -> 0, "Coeffs" -> {{{5}, {0}, {0}}}|>}, 0, 0, 3];
lim3 = DiffExp2`Integrate`EndpointSectorLimit[ls3];
assert["endpoint_limit_drop", esC[lim3[[1]], 0] === 7];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
