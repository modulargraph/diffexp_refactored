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

(* Primitive work depth is relative to the coefficient-window width.  A
   one-row negative Laurent window times the m=-1 dimensional pole is a
   valid eps^-6 result, not a request to truncate the primitive at -3. *)
ls1neg = mk1[{<|"a" -> -1, "b" -> 1, "p" -> 0,
  "Coeffs" -> {{{1}}}|>}, -5, -5, 1];
ri1neg = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  ls1neg, {0, 1/2}]];
assert["ils_negative_absolute_window_primitive_budget",
  !FailureQ[ri1neg] &&
  ri1neg["Values"][[1]]["EpsWindow"] ===
    <|"Min" -> -6, "CompleteMax" -> -6|> &&
  esC[ri1neg["Values"][[1]], -6] === 1];

(* interior PV: f = 1/t over [-1/4, 1/2]: PV = Log[(1/2)/(1/4)] = Log 2 *)
ls2 = mk1[{<|"a" -> -1, "b" -> 0, "p" -> 0,
  "Coeffs" -> {{{1}, {0}, {0}}}|>}, 0, 0, 3];
ri2 = DiffExp2`Integrate`IntegrateLocalSolution[ls2, {-1/4, 1/2}];
assert["ils_interior_pv", eqQ[esC[ri2["Values"][[1]], 0], Log[2]]];

(* An exact-zero fractional-a tower is structurally inactive and therefore
   does not require either the b=0 real-PV parity rule or branch data. *)
ls2zeroFrac = mk1[{
  <|"a" -> 1/2, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{0}, {0}, {0}}}|>,
  <|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{1}, {0}, {0}}}|>
  }, 0, 0, 3];
ri2zeroFrac = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  ls2zeroFrac, {-1/4, 1/2}]];
assert["ils_interior_zero_fractional_sector_inactive",
  !FailureQ[ri2zeroFrac] &&
    eqQ[esC[ri2zeroFrac["Values"][[1]], 0], 3/4]];

(* Interior i-delta pairing: the m=-1,b!=0 cell is a SINGLE regular
   epsilon series.  Exercise the generic p formula at several log depths;
   in every case the separate-arm eps^-1 pole is absent. *)
unitRows[kmax_Integer, ncols_Integer] := Table[
  Table[{Boole[k === 0 && n === 0]}, {n, 0, ncols - 1}], {k, 0, kmax}];
phasePrescription[sg_] := {<|"Factor" -> Global`t, "Sign" -> sg,
  "Multiplicity" -> 1, "LeadingCoeffSign" -> 1|>};
nearQ[x_, y_] := TrueQ[Abs[N[x - y, 60]] < 10^-50];
A0 = 1/5; B0 = 1/4; b0 = 2; kPhase = 3;
Do[
  lsPhase = Join[mk1[{<|"a" -> -1, "b" -> b0, "p" -> pp,
      "Coeffs" -> unitRows[kPhase, 3]|>}, 0, kPhase, 3],
    <|"Prescriptions" -> phasePrescription[1]|>];
  riPhase = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
    lsPhase, {-A0, B0}]];
  leadPhase = (Log[B0]^(pp + 1) - (Log[A0] + I*Pi)^(pp + 1))/
    (pp!*(pp + 1));
  nextPhase = b0*(Log[B0]^(pp + 2) - (Log[A0] + I*Pi)^(pp + 2))/
    (pp!*(pp + 2));
  assert["ils_interior_idelta_p" <> ToString[pp],
    !FailureQ[riPhase] && esMn[riPhase["Values"][[1]]] === pp &&
      nearQ[esC[riPhase["Values"][[1]], pp], leadPhase] &&
      nearQ[esC[riPhase["Values"][[1]], pp + 1], nextPhase]],
  {pp, 0, 3}];

(* The sign choice is physical branch data, not a silent +i0 default. *)
lsPhaseMinus = Join[mk1[{<|"a" -> -1, "b" -> b0, "p" -> 0,
    "Coeffs" -> unitRows[kPhase, 3]|>}, 0, kPhase, 3],
  <|"Prescriptions" -> phasePrescription[-1]|>];
riPhaseMinus = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  lsPhaseMinus, {-A0, B0}]];
assert["ils_interior_idelta_sigma_minus",
  !FailureQ[riPhaseMinus] && esMn[riPhaseMinus["Values"][[1]]] === 0 &&
    nearQ[esC[riPhaseMinus["Values"][[1]], 0], Log[B0/A0] + I*Pi]];

lsPhaseMissing = mk1[{<|"a" -> -1, "b" -> b0, "p" -> 0,
  "Coeffs" -> unitRows[kPhase, 3]|>}, 0, kPhase, 3];
riPhaseMissing = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  lsPhaseMissing, {-A0, B0}]];
assert["ils_interior_idelta_missing_prescription_loud",
  FailureQ[riPhaseMissing] && riPhaseMissing["ID"] === "E3"];

(* Away from m=-1, the same branch-resolved antiderivative difference
   supplies the full phase (including b eps) without a pole special case. *)
lsPhaseRegular = Join[mk1[{<|"a" -> 0, "b" -> b0, "p" -> 0,
    "Coeffs" -> unitRows[kPhase, 3]|>}, 0, kPhase, 3],
  <|"Prescriptions" -> phasePrescription[1]|>];
riPhaseRegular = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  lsPhaseRegular, {-A0, B0}]];
regular1 = b0*(B0*(Log[B0] - 1) +
  A0*(Log[A0] + I*Pi - 1));
assert["ils_interior_idelta_regular_cell_phase",
  !FailureQ[riPhaseRegular] &&
    nearQ[esC[riPhaseRegular["Values"][[1]], 0], A0 + B0] &&
    nearQ[esC[riPhaseRegular["Values"][[1]], 1], regular1]];

(* EndpointSectorLimit: drop rule *)
ls3 = mk1[{<|"a" -> 0, "b" -> 0, "p" -> 0, "Coeffs" -> {{{7}, {1}, {0}}}|>,
  <|"a" -> -1, "b" -> 2, "p" -> 0, "Coeffs" -> {{{5}, {0}, {0}}}|>}, 0, 0, 3];
lim3 = DiffExp2`Integrate`EndpointSectorLimit[ls3];
assert["endpoint_limit_drop", esC[lim3[[1]], 0] === 7];

(* Canonical integer-a merging may bury the finite t^0 coefficient at
   n=-a.  The negative-power cell is zero and the positive cell vanishes. *)
ls4 = mk1[{
  <|"a" -> -1, "b" -> 0, "p" -> 0, "Coeffs" -> {{{0}, {2}, {9}}}|>,
  <|"a" -> 0, "b" -> 0, "p" -> 0, "Coeffs" -> {{{5}, {0}, {0}}}|>
  }, 0, 0, 3];
lim4 = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls4]];
assert["endpoint_limit_buried_finite_after_integer_merge",
  !FailureQ[lim4] && esC[lim4[[1]], 0] === 7];

(* A finite cell in a higher-a sector must survive even when shifting it
   into a far lower-a sector would exceed the fixed Taylor width. *)
ls4wide = mk1[{
  <|"a" -> -5, "b" -> 0, "p" -> 0, "Coeffs" -> {{{0}, {0}, {0}}}|>,
  <|"a" -> 0, "b" -> 0, "p" -> 0, "Coeffs" -> {{{23}, {0}, {0}}}|>
  }, 0, 0, 3];
lim4wide = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[
  DiffExp2`SectorSeries`CanonicalizeLocalSolution[ls4wide]]];
assert["endpoint_limit_no_fixed_width_canonicalization_loss",
  FailureQ[lim4wide] && lim4wide["ID"] === "E10"];

(* With enough Taylor columns the same towers merge without loss and the
   buried t^0 coefficient is certified. *)
ls4covered = mk1[{
  <|"a" -> -5, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{0}, {0}, {0}, {0}, {0}, {0}}}|>,
  <|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{23}, {0}, {0}, {0}, {0}, {0}}}|>
  }, 0, 0, 6];
lim4covered = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[
  DiffExp2`SectorSeries`CanonicalizeLocalSolution[ls4covered]]];
assert["endpoint_limit_covered_far_shift_finite",
  !FailureQ[lim4covered] && esC[lim4covered[[1]], 0] === 23];

(* If a nonzero negative-a tower has not been expanded through its possible
   finite t^0 cell, cancellation of the visible poles is insufficient to
   certify a limit. *)
ls4short = mk1[{
  <|"a" -> -5, "b" -> 0, "p" -> 0, "Coeffs" -> {{{1}, {0}, {0}}}|>,
  <|"a" -> -5, "b" -> 0, "p" -> 0, "Coeffs" -> {{{-1}, {0}, {0}}}|>
  }, 0, 0, 3];
lim4short = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls4short]];
assert["endpoint_limit_incomplete_negative_tower_loud",
  FailureQ[lim4short] && lim4short["ID"] === "E10"];
int4short = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[ls4short, {0, 1/2}]];
assert["endpoint_integral_incomplete_negative_tower_loud",
  FailureQ[int4short] && int4short["ID"] === "E10"];

(* The t^-1 cells cancel after integer-spaced sectors are merged; their
   finite t^0 cells add. *)
ls5 = mk1[{
  <|"a" -> -2, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{0}, {1}, {2}}}|>,
  <|"a" -> -1, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{-1}, {5}, {0}}}|>
  }, 0, 0, 3];
lim5 = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls5]];
assert["endpoint_limit_merged_negative_power_cancellation",
  !FailureQ[lim5] && esC[lim5[[1]], 0] === 7];

(* Positive powers times logs vanish; they are unrelated to the divergent
   t^0 Log[t] cell and must not cause a false positive. *)
ls6 = mk1[{<|"a" -> 1, "b" -> 0, "p" -> 2,
  "Coeffs" -> {{{11}, {13}, {0}}}|>}, 0, 0, 3];
lim6 = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls6]];
assert["endpoint_limit_positive_power_log_vanishes",
  !FailureQ[lim6] && esC[lim6[[1]], 0] === 0];

(* Non-cancelling b=0 endpoint divergences remain loud. *)
ls7 = mk1[{<|"a" -> 0, "b" -> 0, "p" -> 1,
  "Coeffs" -> {{{1}, {0}, {0}}}|>}, 0, 0, 3];
lim7 = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls7]];
assert["endpoint_limit_log_divergence_loud",
  FailureQ[lim7] && lim7["ID"] === "E2"];

ls8 = mk1[{<|"a" -> -1, "b" -> 0, "p" -> 0,
  "Coeffs" -> {{{1}, {17}, {19}}}|>}, 0, 0, 3];
lim8 = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[ls8]];
assert["endpoint_limit_negative_power_divergence_loud",
  FailureQ[lim8] && lim8["ID"] === "E2"];

mkNearCancelledDivergence[delta_] := mk1[{
  <|"a" -> -2, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{0}, {1}, {0}}}|>,
  <|"a" -> -1, "b" -> 0, "p" -> 0,
    "Coeffs" -> {{{-1 + delta}, {0}, {0}}}|>
  }, 0, 0, 3];
lsDiagonalDiv = mkNearCancelledDivergence[(8 + 8 I)*10^-25];
lsUncertainDiv = mkNearCancelledDivergence[0``17];
limDiagonalDiv = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[lsDiagonalDiv]];
intDiagonalDiv = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  lsDiagonalDiv, {0, 1/2}]];
limUncertainDiv = catchDE2[DiffExp2`Integrate`EndpointSectorLimit[lsUncertainDiv]];
intUncertainDiv = catchDE2[DiffExp2`Integrate`IntegrateLocalSolution[
  lsUncertainDiv, {0, 1/2}]];
assert["endpoint_cancellation_true_modulus_and_uncertainty_loud",
  And @@ (FailureQ[#] && #["ID"] === "E2" & /@
    {limDiagonalDiv, intDiagonalDiv, limUncertainDiv, intUncertainDiv})];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
