(* Smoke tests for DiffExp2/Transport.m (M4 parity suite grows at M5). *)
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
  "ExpansionOrder" -> 30, "EpsilonOrder" -> 2}]];
x = Global`x; eps = Global`eps;

(* TT1: regular transport f' = f from 0 to 1/2: exact e^(1/2) *)
sys = <|"Matrix" -> {{1}}, "Variable" -> x, "SingularFactors" -> {}|>;
plan = DiffExp2`Transport`SegmentLine[sys, {0, 1/2}];
res = catchDE2[DiffExp2`Transport`TransportLine[sys, {{1, 0, 0}}, plan]];
assert["tt1_regular_exponential",
  !FailureQ[res] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res["Value"], 0][[1]] - Exp[1/2], 30]] < 10^-25];

(* TT2: singular endpoint x^eps: transport (1/2)^eps from 1/2 to 0 *)
sys2 = <|"Matrix" -> {{eps/x}}, "Variable" -> x, "SingularFactors" -> {x}|>;
plan2 = DiffExp2`Transport`SegmentLine[sys2, {1/2, 0}];
bvals = Transpose[Table[{SeriesCoefficient[(1/2)^eps, {eps, 0, k}]}, {k, 0, 2}]];
res2 = catchDE2[DiffExp2`Transport`TransportLine[sys2, bvals, plan2]];
assert["tt2_singular_endpoint_weight",
  !FailureQ[res2] && res2["EndpointIsSingular"] &&
  Module[{f = res2["Final"], sec},
    sec = SelectFirst[f["Sectors"], PossibleZeroQ[#["b"] - 1] &];
    sec =!= Missing["NotFound"] &&
    Abs[N[sec["Coeffs"][[-f["EpsWindow", "Min"] + 1, 1, 1]] - 1, 20]] < 10^-20]];

(* The same singular endpoint approached on its negative local arm needs an
   explicit sheet even though no far-side crossing follows.  The leading
   x^eps weight carries both the affine chart scale and the prescribed
   negative-axis phase. *)
negativeEndpointResult[sg_] := Module[{p, r},
  catchDE2[DiffExp2`Config`UpdateConfiguration[{
    "DeltaPrescriptions" -> {{x, sg}}}]];
  p = DiffExp2`Transport`SegmentLine[sys2, {-1, 0}];
  r = catchDE2[DiffExp2`Transport`TransportLine[
    sys2, {{1, 0, 0}}, p]];
  r];
res2negPlus = negativeEndpointResult[1];
res2negMinus = negativeEndpointResult[-1];
negativeEndpointLeading[r_] := Module[{f = r["Final"], sec},
  sec = SelectFirst[f["Sectors"], PossibleZeroQ[# ["b"] - 1] &];
  Table[sec["Coeffs"][[k - f["EpsWindow", "Min"] + 1, 1, 1]],
    {k, 0, 2}]];
assert["tt2_negative_arm_endpoint_plus_minus_i0",
  !FailureQ[res2negPlus] && !FailureQ[res2negMinus] &&
    Module[{cp = negativeEndpointLeading[res2negPlus],
        cm = negativeEndpointLeading[res2negMinus], qp, qm},
      qp = Log[2] - I Pi; qm = Log[2] + I Pi;
      And @@ Table[
        Abs[N[cp[[k + 1]] - qp^k/k!, 50]] < 10^-15 &&
        Abs[N[cm[[k + 1]] - qm^k/k!, 50]] < 10^-15,
        {k, 0, 2}]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {}}]];
res2negMissing = catchDE2[DiffExp2`Transport`TransportLine[
  sys2, {{1, 0, 0}}, DiffExp2`Transport`SegmentLine[sys2, {-1, 0}]]];
assert["tt2_negative_arm_endpoint_missing_prescription_loud",
  FailureQ[res2negMissing] && res2negMissing["ID"] === "E8"];

(* TT3: regular line with a far singularity: f' = f/(x-2): f = c(x-2);
   f(0) = 1/2 -> c = -1/4 -> f(1) = 1/4 *)
sys3 = <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x, "SingularFactors" -> {x - 2}|>;
plan3 = DiffExp2`Transport`SegmentLine[sys3, {0, 1}];
res3 = catchDE2[DiffExp2`Transport`TransportLine[sys3, {{1/2, 0, 0}}, plan3]];
assert["tt3_far_singularity",
  !FailureQ[res3] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res3["Value"], 0][[1]] - 1/4, 30]] < 10^-25];

(* TT3a: RadiusOfConvergence is an affine coordinate normalization, not a
   no-op.  The physical line and answer stay fixed while x=center+Scale t,
   and LocalSolution's radius is expressed in t. *)
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "RadiusOfConvergence" -> 10}]];
plan3a = DiffExp2`Transport`SegmentLine[sys3, {0, 1}];
res3a = catchDE2[DiffExp2`Transport`TransportLine[
  sys3, {{1/2, 0, 0}}, plan3a]];
assert["tt3a_affine_roc_scale_and_local_radius",
  plan3a["Charts"][[1, "Scale"]] === 1/5 &&
  plan3a["Charts"][[1, "Radius"]] === 2 &&
  res3a["Charts"][[1, "LocalSolution", "Radius"]] === 10];
assert["tt3a_affine_roc_transport_invariant",
  !FailureQ[res3a] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res3a["Value"], 0][[1]] - 1/4,
    30]] < 10^-25];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "RadiusOfConvergence" -> 1}]];

(* TT3-eps: planning uses the first nonzero eps coefficient of each exact
   singular factor, after stripping any overall eps valuation.  This is a
   planner projection only: it prevents regulator-dependent moving roots
   from entering line geometry while retaining their eps=0 limit exactly. *)
epsValuationFactors3 = DiffExp2`Transport`EpsilonZeroSingularFactors[
  {eps*(x + eps)}, x];
epsValuationSings3 = DiffExp2`Transport`FindSingularities[<|
  "Matrix" -> {{0}}, "Variable" -> x, "SingularFactors" -> {},
  "ExtraSingularFactors" -> {eps*(x + eps)}|>];
assert["tt3_eps_overall_valuation_keeps_first_nonzero_coefficient",
  epsValuationFactors3 === {x} &&
  Keys[epsValuationSings3["Factors"]] === {x} &&
  epsValuationSings3["All"] === {0}];

epsShiftSings3 = DiffExp2`Transport`FindSingularities[<|
  "Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x + eps}|>];
assert["tt3_eps_shifted_root_projects_exactly_to_zero",
  Keys[epsShiftSings3["Factors"]] === {x} &&
  epsShiftSings3["All"] === {0} && epsShiftSings3["Real"] === {0}];

assert["tt3_epsilon_alias_projects_to_the_same_alphabet",
  DiffExp2`Transport`EpsilonZeroSingularFactors[
    {x + Global`\[Epsilon]}, x] === {x}];
assert["tt3_non_epsilon_regulator_is_not_specialized",
  DiffExp2`Transport`EpsilonZeroSingularFactors[
    {x + Global`etaReg}, x] === {x + Global`etaReg}];

epsLegacyPrecedence3 = DiffExp2`Transport`FindSingularities[<|
  "Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x - 1},
  "SingularFactorsExact" -> {x + eps}|>];
assert["tt3_exact_factor_metadata_wins_over_legacy_planner_field",
  Keys[epsLegacyPrecedence3["Factors"]] === {x} &&
    epsLegacyPrecedence3["All"] === {0}];

projectedPrescription3[sg_] := Module[{p, final},
  catchDE2[DiffExp2`Config`UpdateConfiguration[{
    "DeltaPrescriptions" -> {{x + eps, sg}}}]];
  p = catchDE2[DiffExp2`Transport`SegmentLine[<|
    "Matrix" -> {{0}}, "Variable" -> x, "SingularFactors" -> {},
    "ExtraSingularFactors" -> {x + eps}|>, {1/2, 0}]];
  final = If[FailureQ[p], None,
    SelectFirst[p["Charts"], # ["Center"] === 0 &, None]];
  final];
epsPrescriptionPlus3 = projectedPrescription3[1];
epsPrescriptionMinus3 = projectedPrescription3[-1];
assert["tt3_projected_extra_factor_keeps_plus_minus_idelta_sheet",
  And @@ MapThread[Function[{chart, sg},
    AssociationQ[chart] && Length[chart["Prescriptions"]] === 1 &&
      chart["Prescriptions"][[1, "Factor"]] === x &&
      chart["Prescriptions"][[1, "ExactFactor"]] === x + eps &&
      chart["Prescriptions"][[1, "Sign"]] === sg &&
      chart["Prescriptions"][[1, "Multiplicity"]] === 1 &&
      chart["Prescriptions"][[1, "LeadingCoeffSign"]] === 1],
    {{epsPrescriptionPlus3, epsPrescriptionMinus3}, {1, -1}}]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {}}]];

epsOutsideSys3 = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x - 2 + eps}|>;
epsOutsideSings3 = DiffExp2`Transport`FindSingularities[epsOutsideSys3];
epsOutsidePlan3 = catchDE2[
  DiffExp2`Transport`SegmentLine[epsOutsideSys3, {0, 1}]];
assert["tt3_eps_moving_root_outside_line_uses_limit_root",
  Keys[epsOutsideSings3["Factors"]] === {x - 2} &&
  epsOutsideSings3["All"] === {2} &&
  !FailureQ[epsOutsidePlan3] &&
  AllTrue[epsOutsidePlan3["Charts"], !TrueQ[# ["Singular"]] &]];

epsInfinitySings3 = DiffExp2`Transport`FindSingularities[<|
  "Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {1 + eps*x}|>];
assert["tt3_eps_degree_drop_sends_root_to_infinity",
  DiffExp2`Transport`EpsilonZeroSingularFactors[{1 + eps*x}, x] === {} &&
  epsInfinitySings3["All"] === {} &&
  epsInfinitySings3["Factors"] === <||>];

noEpsSings3 = DiffExp2`Transport`FindSingularities[<|
  "Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x^2 - 1}|>];
assert["tt3_epsilon_independent_factor_parity",
  DiffExp2`Transport`EpsilonZeroSingularFactors[{x^2 - 1}, x] ===
    {x^2 - 1} &&
  Keys[noEpsSings3["Factors"]] === {x^2 - 1} &&
  Sort[noEpsSings3["All"]] === {-1, 1}];

(* TT3b: faithfully retain the old DiffExp Re(z)+/-Im(z) real projection
   waypoints while true complex roots alone determine convergence radii.
   The ghosts are regular centers, never fake singular charts. *)
sys3b = <|"Matrix" -> {{1/(x^2 + 1)}}, "Variable" -> x,
  "SingularFactors" -> {x^2 + 1}|>;
sings3b = DiffExp2`Transport`FindSingularities[sys3b];
plan3b = DiffExp2`Transport`SegmentLine[sys3b, {-3/2, 3/2}];
centers3b = # ["Center"] & /@ plan3b["Charts"];
assert["tt3b_complex_projection_waypoints_restored",
  sings3b["Projected"] === {-1, 0, 1} &&
  And @@ (MemberQ[centers3b, #] & /@ {-1, 0, 1})];
assert["tt3b_ghosts_regular_true_radius_preserved",
  DiffExp2`Transport`ChartRadius[0, sings3b["All"]] === 1 &&
  AllTrue[Select[plan3b["Charts"], MemberQ[{-1, 0, 1}, #["Center"]] &],
    !TrueQ[# ["Singular"]] &]];
ghostOne3b = SelectFirst[plan3b["Charts"], # ["Center"] === 1 &];
assert["tt3b_true_radius_and_projected_match_scale_are_distinct",
  ghostOne3b["Radius"] === Sqrt[2] &&
  ghostOne3b["MatchRadius"] === 1 && ghostOne3b["Scale"] === 1 &&
  ghostOne3b["LocalRadius"] === Sqrt[2]];
sys3bSuppressed = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x^2 + 1, x - 1/2}|>;
sings3bSuppressed = DiffExp2`Transport`FindSingularities[sys3bSuppressed];
assert["tt3b_occupied_complex_projection_strip_suppresses_ghost",
  MemberQ[sings3bSuppressed["Projected"], -1] &&
  MemberQ[sings3bSuppressed["Projected"], 0] &&
  MemberQ[sings3bSuppressed["Projected"], 1/2] &&
  !MemberQ[sings3bSuppressed["Projected"], 1]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{"DivisionOrder" -> 2}]];
sys3bK2 = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {(x - Sqrt[2])^2 + 3}|>;
plan3bK2Forced = DiffExp2`Transport`SegmentLine[sys3bK2, {-1/2, 2}];
plan3bK2Endpoint = DiffExp2`Transport`SegmentLine[
  sys3bK2, {7/5, 22661/10000}];
assert["tt3b_k2_algebraic_projection_respects_true_half_radius",
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[plan3bK2Forced]]] &&
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[plan3bK2Endpoint]]] &&
  Module[{last = Last[plan3bK2Endpoint["Charts"]]},
    TrueQ[Abs[N[plan3bK2Endpoint["To"] - last["Center"], 40]] <=
      N[last["Radius"], 40]/2]]];
sys3bProjectedEndpoint = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {(x + 3)^2 + 9}|>;
plan3bProjectedEndpointR = catchDE2[DiffExp2`Transport`SegmentLine[
  sys3bProjectedEndpoint, {-5, 0}]];
plan3bProjectedEndpointL = catchDE2[DiffExp2`Transport`SegmentLine[
  sys3bProjectedEndpoint, {1, 0}]];
assert["tt3b_regular_projected_endpoint_terminates_both_directions",
  !FailureQ[plan3bProjectedEndpointR] &&
  !FailureQ[plan3bProjectedEndpointL] &&
  Last[plan3bProjectedEndpointR["Charts"]]["Center"] === 0 &&
  Last[plan3bProjectedEndpointL["Charts"]]["Center"] === 0 &&
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[
    plan3bProjectedEndpointR]]] &&
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[
    plan3bProjectedEndpointL]]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{"DivisionOrder" -> 3}]];
symPairs3b = Select[Range[2, Length[plan3b["Charts"]]],
  TrueQ[Lookup[plan3b["Charts"][[#]], "SymmetricMatch", False]] &];
assert["tt3b_classic_both_sides_plus_minus_one_over_k",
  symPairs3b =!= {} && AllTrue[symPairs3b, Function[i,
    Module[{a = plan3b["Charts"][[i - 1]], b = plan3b["Charts"][[i]],
        mp = plan3b["Charts"][[i, "IncomingMatchPoint"]]},
      Together[(mp - a["Center"])/a["MatchRadius"]] === 1/3 &&
      Together[(mp - b["Center"])/b["MatchRadius"]] === -1/3]]]];
plan3bReverse = DiffExp2`Transport`SegmentLine[sys3b, {3/2, -3/2}];
symPairs3bReverse = Select[Range[2, Length[plan3bReverse["Charts"]]],
  TrueQ[Lookup[plan3bReverse["Charts"][[#]], "SymmetricMatch", False]] &];
assert["tt3b_classic_symmetric_match_reverse_direction",
  symPairs3bReverse =!= {} && AllTrue[symPairs3bReverse, Function[i,
    Module[{a = plan3bReverse["Charts"][[i - 1]],
        b = plan3bReverse["Charts"][[i]],
        mp = plan3bReverse["Charts"][[i, "IncomingMatchPoint"]]},
      Together[(mp - a["Center"])/a["MatchRadius"]] === -1/3 &&
      Together[(mp - b["Center"])/b["MatchRadius"]] === 1/3]]]];

(* TT3c: Config's sign-aware DeltaPrescriptions must reach the singular
   chart records.  This was parsed but silently discarded by SegmentLine. *)
catchDE2[DiffExp2`Config`UpdateConfiguration[{"DeltaPrescriptions" ->
  {{x, 1}, {x - 1/2, -1}, {1 - x, 1}}}]];
sys3c = <|"Matrix" -> {{eps/(x (x - 1/2) (1 - x))}}, "Variable" -> x,
  "SingularFactors" -> {x, x - 1/2, 1 - x}|>;
plan3c = DiffExp2`Transport`SegmentLine[sys3c, {11/23, 1}];
chartHalf3c = SelectFirst[plan3c["Charts"], # ["Center"] === 1/2 &];
chartOne3c = SelectFirst[plan3c["Charts"], # ["Center"] === 1 &];
assert["tt3c_delta_prescriptions_reach_charts",
  AssociationQ[chartHalf3c] && AssociationQ[chartOne3c] &&
  Length[chartHalf3c["Prescriptions"]] === 1 &&
  Length[chartOne3c["Prescriptions"]] === 1];
assert["tt3c_plus_minus_im_sign_derived",
  DiffExp2`SectorSeries`ChartImSign[chartHalf3c] === -1 &&
  DiffExp2`SectorSeries`ChartImSign[chartOne3c] === -1];
catchDE2[DiffExp2`Config`UpdateConfiguration[{"DeltaPrescriptions" -> {}}]];

(* TT4 (regression, banana level-1 geometry): an interior SINGULAR chart
   whose radius dwarfs the producing chart's.  Singular factors x and
   (x - 1/2) put the anchor 11/23 at radius 1/46 while the chart at 1/2
   has radius 1/2: the naive incoming match point 1/2 - (1/2)/k lies far
   BEHIND the anchor, outside its disk.  Contract: success with finite
   values OR the loud E8 — never the raw "point outside the chart radius"
   evaluation error.  Exact value: y' = (eps/x) y, y(11/23) = (11/23)^eps
   -> y = x^eps -> y(1) = 1 at every eps order. *)
sys4 = <|"Matrix" -> {{eps/x}}, "Variable" -> x, "SingularFactors" -> {x}|>;
bv4 = Transpose[Table[{SeriesCoefficient[(11/23)^eps, {eps, 0, kk}]}, {kk, 0, 2}]];
res4 = catchDE2[DiffExp2`API`TransportEndpoint[sys4, bv4, 11/23, 1,
  "ExtraSingularFactors" -> {x - 1/2}]];
assert["tt4_singular_handoff_contract",
  If[FailureQ[res4], res4["ID"] === "E8", True]];
assert["tt4_singular_handoff_value",
  !FailureQ[res4] && Module[{v = res4["Value"]},
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, 0][[1]] - 1, 30]] < 10^-10 &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, 1][[1]], 30]] < 10^-10 &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, 2][[1]], 30]] < 10^-10]];
plan4 = DiffExp2`Transport`SegmentLine[
  Join[sys4, <|"ExtraSingularFactors" -> {x - 1/2}|>], {11/23, 1}];
assert["tt4_plan_marks_half_singular",
  AnyTrue[plan4["Charts"], TrueQ[#["Singular"]] && #["Center"] === 1/2 &]];
assert["tt4_validateplan_accepts_segmentline_plan",
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[plan4]]] &&
  plan4["SegmentCount"] === Length[plan4["Charts"]]];

(* TT5: the mirror direction — singular ENDPOINT at 0 from the same
   squeezed anchor (the trLo hazard).  Same contract; the eps^0 weight of
   the b = 1 sector at the 0-chart is exactly 1 (the eps^0 component of
   x^eps is constant, so no truncation enters that order). *)
res5 = catchDE2[DiffExp2`API`TransportEndpoint[sys4, bv4, 11/23, 0,
  "ExtraSingularFactors" -> {x - 1/2}]];
assert["tt5_singular_endpoint_contract",
  If[FailureQ[res5], res5["ID"] === "E8", True]];
assert["tt5_singular_endpoint_weight",
  !FailureQ[res5] && res5["EndpointIsSingular"] &&
  Module[{f = res5["Final"], sec},
    sec = SelectFirst[f["Sectors"], PossibleZeroQ[#["b"] - 1] &];
    sec =!= Missing["NotFound"] &&
    Abs[N[sec["Coeffs"][[-f["EpsWindow", "Min"] + 1, 1, 1]] - 1, 20]] < 10^-20]];

(* TT5b: branch phase is applied exactly once across an interior singularity.
   With +(i0), the negative arm has Log[-t]+i Pi.  Rightward transport from
   1/4 to 3/4 starts on that arm and ends at equal positive distance, so the
   normalized exact factor is Exp[-i Pi eps].
   The reverse path is its inverse.  This catches a former rightward-only
   double application of ApplyCrossing. *)
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {{x - 1/2, 1}}}]];
sys5b = <|"Matrix" -> {{eps/(x - 1/2)}}, "Variable" -> x,
  "SingularFactors" -> {x - 1/2}|>;
res5bRight = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5b, {{1, 0, 0}}, 1/4, 3/4]];
res5bLeft = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5b, {{1, 0, 0}}, 3/4, 1/4]];
assert["tt5b_rightward_crossing_phase_once",
  !FailureQ[res5bRight] && Module[{v = res5bRight["Value"],
      q = -I Pi},
    And @@ Table[Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, k][[1]] -
      q^k/k!, 50]] < 10^-8, {k, 0, 2}]]];
assert["tt5b_leftward_crossing_phase_once",
  !FailureQ[res5bLeft] && Module[{v = res5bLeft["Value"],
      q = I Pi},
    And @@ Table[Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, k][[1]] -
      q^k/k!, 50]] < 10^-8, {k, 0, 2}]]];

(* End-to-end noninteger-a crossing: this is the direct fractional-power
   analogue of the b eps phase above, and pins both physical sheets. *)
sys5bFrac = <|"Matrix" -> {{1/(2 (x - 1/2))}}, "Variable" -> x,
  "SingularFactors" -> {x - 1/2}|>;
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {{x - 1/2, 1}}}]];
res5bFracPlus = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5bFrac, {{1, 0, 0}}, 1/4, 3/4]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {{x - 1/2, -1}}}]];
res5bFracMinus = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5bFrac, {{1, 0, 0}}, 1/4, 3/4]];
assert["tt5b_fractional_crossing_plus_minus_i0",
  !FailureQ[res5bFracPlus] && !FailureQ[res5bFracMinus] &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[
      res5bFracPlus["Value"], 0][[1]] + I, 40]] < 10^-15 &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[
      res5bFracMinus["Value"], 0][[1]] - I, 40]] < 10^-15];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {}}]];
res5bMissing = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5b, {{1, 0, 0}}, 1/4, 3/4]];
assert["tt5b_rightward_crossing_missing_prescription_loud",
  FailureQ[res5bMissing] && res5bMissing["ID"] === "E8"];

(* TT5c: when the mandatory regular chart after a singular crossing is
   clipped to a nearby endpoint, its incoming point must be recomputed from
   the ACTUAL receiver.  Retaining the natural pre-clip +R/k point can put it
   outside the receiver disk (endpoint 1/10), or near 0.97 R and silently
   outrun the half-radius error model (endpoint 51/100). *)
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {{x, 1}}}]];
sys5c = <|"Matrix" -> {{eps/x}}, "Variable" -> x,
  "SingularFactors" -> {x}|>;
plan5cNear = DiffExp2`Transport`SegmentLine[sys5c, {-1, 1/10}];
last5cNear = Last[plan5cNear["Charts"]];
assert["tt5c_clipped_receiver_recomputes_shared_match_point",
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[plan5cNear]]] &&
  TrueQ[PossibleZeroQ[Together[
    last5cNear["IncomingMatchPoint"] - last5cNear["Center"] +
      last5cNear["MatchRadius"]/3]]]];
res5c = catchDE2[DiffExp2`API`TransportEndpoint[
  sys5c, {{1, 0, 0}}, -1, 51/100]];
assert["tt5c_clipped_receiver_transport_matches_xeps_ratio",
  !FailureQ[res5c] && Module[{v = res5c["Value"],
      q = Log[51/100] - I Pi},
    And @@ Table[Abs[N[DiffExp2`EpsSeries`ESCoefficient[v, k][[1]] -
      q^k/k!, 50]] < 10^-8, {k, 0, 2}]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {}}]];

(* TT6: ValidatePlan rejects a hand-built plan whose singular handoff is
   geometrically impossible (producing margin-disk and approach interval
   disjoint: gap 1/46 >= (9/10)/46 + 1/1000): loud E8 with the chain
   geometry, not a deep evaluation error. *)
plan6 = <|"From" -> 11/23, "To" -> 1, "Direction" -> 1,
  "Charts" -> {
    <|"Center" -> 11/23, "Singular" -> False, "Radius" -> 1/46,
      "Name" -> "tt6anchor"|>,
    <|"Center" -> 1/2, "Singular" -> True, "Radius" -> 1/1000,
      "Name" -> "tt6sing"|>}|>;
res6 = catchDE2[DiffExp2`Transport`ValidatePlan[plan6]];
assert["tt6_validateplan_e8_on_impossible_handoff",
  FailureQ[res6] && res6["ID"] === "E8"];
plan6b = <|"From" -> 0, "To" -> 2, "Direction" -> 1,
  "Charts" -> {
    <|"Center" -> 0, "Singular" -> False, "Radius" -> 1,
      "Name" -> "tt6bproducer"|>,
    <|"Center" -> 2, "Singular" -> False, "Radius" -> 1/10,
      "IncomingMatchPoint" -> 1/2, "Name" -> "tt6breceiver"|>}|>;
res6b = catchDE2[DiffExp2`Transport`ValidatePlan[plan6b]];
assert["tt6_validateplan_checks_regular_receiving_disk",
  FailureQ[res6b] && res6b["ID"] === "E8"];
plan6c = <|"From" -> 0, "To" -> 1, "Direction" -> 1,
  "Charts" -> {
    <|"Center" -> 0, "Singular" -> False, "Radius" -> 2,
      "Name" -> "tt6cproducer"|>,
    <|"Center" -> 1, "Singular" -> False, "Radius" -> 2,
      "IncomingMatchPoint" -> 7/4, "Name" -> "tt6creceiver"|>}|>;
res6c = catchDE2[DiffExp2`Transport`ValidatePlan[plan6c]];
assert["tt6_validateplan_checks_half_radius_probe_envelope",
  FailureQ[res6c] && res6c["ID"] === "E8"];
tinyPlan6c = <|"From" -> 1, "To" -> 1 + 10^-101,
  "Direction" -> 1, "EndpointIsSingular" -> False,
  "Charts" -> {<|"Center" -> 1, "Singular" -> False,
    "Radius" -> 10^-100, "Name" -> "tt6ctiny"|>}|>;
assert["tt6_exact_distance_preserves_tiny_scale_and_anchor_zero",
  DiffExp2`Transport`Private`numericDistance[
    1 + 10^-100, 1, 40] > 0 &&
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[tinyPlan6c]]]];
plan6d = plan4;
singularPos6d = FirstPosition[plan6d["Charts"],
  c_ /; TrueQ[c["Singular"]], Missing["NotFound"]];
If[singularPos6d =!= Missing["NotFound"],
  charts6d = plan6d["Charts"];
  charts6d[[First[singularPos6d]]] = Join[
    charts6d[[First[singularPos6d]]], <|"IncomingMatchPoint" ->
      charts6d[[First[singularPos6d], "IncomingMatchPoint"]] + 1/1000|>];
  plan6d = Join[plan6d, <|"Charts" -> charts6d|>]];
res6d = catchDE2[DiffExp2`Transport`ValidatePlan[plan6d]];
assert["tt6_validateplan_recomputes_canonical_singular_handoff",
  singularPos6d =!= Missing["NotFound"] && FailureQ[res6d] &&
  res6d["ID"] === "E8"];

(* TT7: bounded recombination safety.  The same-a/different-b mock family
   is eligible for the legacy divided-difference heuristic when hit-free,
   but an uncompensated pseudo-collision diagnostic must preserve the
   original Laurent basis exactly. *)
mkMockColumn[b_, vec_] := <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 1,
  "Sectors" -> {<|"a" -> 0, "b" -> b, "p" -> 0,
    "Coeffs" -> Table[{If[k === 0, #, 0] & /@ vec}, {k, 0, 2}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 0|>,
  "ErrorEstimate" -> {0, 0, 0}, "Prescriptions" -> {}|>;
mockBasis = {mkMockColumn[0, {1, 0}], mkMockColumn[2, {1, 2}]};
mockSpecs = {
  <|"a" -> 0, "b" -> 0, "Family" -> 1, "ChainPos" -> 0|>,
  <|"a" -> 0, "b" -> 2, "Family" -> 1, "ChainPos" -> 0|>};
mockCS = <|"Center" -> 0, "IndicialData" -> <|"Families" -> {
  <|"EpsZeroDegeneracy" -> 1|>}|>|>;
mockHitFree = DiffExp2`Transport`Private`recombineDegenerate[
  mockCS, mockBasis, mockSpecs, <|"PseudoCollisionsHit" -> {}|>];
mockGuarded = DiffExp2`Transport`Private`recombineDegenerate[
  mockCS, mockBasis, mockSpecs,
  <|"PseudoCollisionsHit" -> {<|"n" -> 1, "Cols" -> {1}|>}|>];
mockCompensated = DiffExp2`Transport`Private`recombineDegenerate[
  mockCS, mockBasis, mockSpecs,
  <|"PseudoCollisionsHit" -> {<|"n" -> 1, "Cols" -> {1}|>},
    "PseudoCollisionsCompensated" -> True|>];
mockMissing = DiffExp2`Transport`Private`recombineDegenerate[
  mockCS, mockBasis, mockSpecs, <||>];
assert["tt7_hit_free_recombination_still_active",
  mockHitFree =!= mockBasis && mockHitFree[[2]]["EpsWindow", "Min"] === -1];
assert["tt7_pseudo_hit_basis_not_recombined",
  mockGuarded === mockBasis && mockMissing === mockBasis];
assert["tt7_compensated_hit_still_recombines_same_a_pair",
  mockCompensated =!= mockBasis &&
  mockCompensated[[2]]["EpsWindow", "Min"] === -1];

(* TT8: Laurent-field matching with a unit determinant after cancellation.
     F = {{1/eps, 1}, {1/eps+1, 1}}, w = {eps, 1}
   gives analytic v = {2, 2+eps}.  This exercises negative entry windows,
   the expected weights, and the always-on residual assertion. *)
esT[k_, cs_] := DiffExp2`EpsSeries`ESNew[k, cs];
invE = esT[-1, {1, 0, 0, 0, 0}];
one = esT[0, {1, 0, 0, 0}];
invEPlusOne = esT[-1, {1, 1, 0, 0, 0}];
fLaurent = {{invE, one}, {invEPlusOne, one}};
vLaurent = {esT[0, {2, 0, 0, 0}], esT[0, {2, 1, 0, 0}]};
wLaurent = catchDE2[DiffExp2`Transport`MatchWeights[
  fLaurent, vLaurent, "tt8_laurent"]];
assert["tt8_laurent_matching_weights",
  !FailureQ[wLaurent] &&
  DiffExp2`EpsSeries`ESCoefficient[wLaurent[[1]], 1] === 1 &&
  DiffExp2`EpsSeries`ESCoefficient[wLaurent[[1]], 0] === 0 &&
  DiffExp2`EpsSeries`ESCoefficient[wLaurent[[2]], 0] === 1];
laurentOut = If[FailureQ[wLaurent], wLaurent,
  Table[DiffExp2`EpsSeries`ESAdd[
    DiffExp2`EpsSeries`ESTimes[fLaurent[[c, 1]], wLaurent[[1]]],
    DiffExp2`EpsSeries`ESTimes[fLaurent[[c, 2]], wLaurent[[2]]]], {c, 2}]];
assert["tt8_laurent_matching_output",
  !FailureQ[laurentOut] && And @@ MapThread[
    DiffExp2`EpsSeries`ESSameQ, {laurentOut, vLaurent}]];

(* Full Laurent pivoting may exchange basis columns.  Swapping TT8's
   columns forces the lowest-order first pivot into original column 2; the
   returned weights must still be mapped to original basis order. *)
fLaurentSwapped = fLaurent[[All, {2, 1}]];
wLaurentSwapped = catchDE2[DiffExp2`Transport`MatchWeights[
  fLaurentSwapped, vLaurent, "tt8_laurent_column_pivot"]];
assert["tt8_laurent_full_column_pivot_maps_weights_back",
  !FailureQ[wLaurentSwapped] &&
  DiffExp2`EpsSeries`ESCoefficient[wLaurentSwapped[[1]], 0] === 1 &&
  DiffExp2`EpsSeries`ESCoefficient[wLaurentSwapped[[2]], 1] === 1];

(* TT9: exercise the residual checker directly with a deliberately wrong
   weight vector.  The error must name the failed order/component. *)
badWeights = {esT[1, {2, 0, 0}], esT[0, {1, 0, 0}]};
badMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    fLaurent, vLaurent, badWeights, "tt9_bad_match"]];
assert["tt9_matching_residual_loud",
  FailureQ[badMatch] && badMatch["ID"] === "E6" &&
  IntegerQ[badMatch["EpsOrder"]] && IntegerQ[badMatch["Component"]]];

(* TT10: the matching seam must preserve the repository's pinned 2x-WP
   input headroom.  Banana's honest Laurent basis is ill-conditioned and
   can consume far more than the old ad-hoc 20 guard digits. *)
handoffProbe = DiffExp2`Transport`Private`numHandoff[1/3];
assert["tt10_matching_handoff_uses_2x_wp",
  Precision[handoffProbe] >=
    DiffExp2`Tolerances`$InputPrecisionFactor*
      DiffExp2`Config`CFG["WorkingPrecision"] - 1];

(* TT11: row elimination trims structural zeros at LaurentLeadTol, so the
   residual proof must use that same floor when it is looser than MatchTol.
   A resolved residual below the floor passes; one above remains loud. *)
unitF = {{esT[0, {1, 0, 0}]}};
unitV = {esT[0, {1, 0, 0}]};
smallFloorWeights = {esT[0, {N[1 + 10^-26, 50], 0, 0}]};
largeFloorWeights = {esT[0, {N[1 + 10^-22, 50], 0, 0}]};
smallFloorMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, smallFloorWeights, "tt11_small_floor"]];
largeFloorMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, largeFloorWeights, "tt11_large_floor"]];
diagonalFloorWeights = {esT[0,
  {1 + (8 + 8 I)*10^-25, 0, 0}]};
diagonalFloorMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, diagonalFloorWeights, "tt11_true_complex_modulus"]];
diagonalOrderWeights = {esT[0,
  {1, (8 + 8 I)*10^-25, 0}]};
diagonalOrderMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, diagonalOrderWeights, "tt11_true_modulus_order1"]];
insideDiagonalWeights = {esT[0,
  {1, (7 + 7 I)*10^-25, 0}]};
insideDiagonalMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, insideDiagonalWeights, "tt11_inside_true_modulus"]];
uncertainZeroWeights = {esT[0, {1, 0``17, 0}]};
uncertainZeroMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    unitF, unitV, uncertainZeroWeights, "tt11_uncertain_zero"]];
assert["tt11_residual_below_laurent_floor_accepted", smallFloorMatch === Null];
assert["tt11_residual_above_laurent_floor_loud",
  FailureQ[largeFloorMatch] && largeFloorMatch["ID"] === "E6" &&
  largeFloorMatch["EffectiveTolerance"] ===
    DiffExp2`Tolerances`Tol["LaurentLeadTol"] &&
  AssociationQ[largeFloorMatch[[2]]] &&
  KeyExistsQ[largeFloorMatch[[2]], "ResidualUncertainty"]];
assert["tt11_residual_uses_true_complex_modulus",
  FailureQ[diagonalFloorMatch] && diagonalFloorMatch["ID"] === "E6" &&
  FailureQ[diagonalOrderMatch] && diagonalOrderMatch["ID"] === "E6" &&
  diagonalOrderMatch["EpsOrder"] === 1 &&
  insideDiagonalMatch === Null &&
  FailureQ[uncertainZeroMatch] && uncertainZeroMatch["ID"] === "E6"];

(* TT12: matrix dimension alone must not relax the structural proof.  Even
   in a 7-column identity system only one term is active here; accepting the
   deliberately wrong 4e-24 weight would hide a real matching error. *)
id7F = Table[esT[0, {KroneckerDelta[r, c], 0, 0}], {r, 7}, {c, 7}];
id7V = Table[esT[0, {KroneckerDelta[r, 1], 0, 0}], {r, 7}];
wrong7W = Join[{esT[0, {N[1 + 4*10^-24, 80], 0, 0}]},
  Table[esT[0, {0, 0, 0}], {6}]];
wrong7Match = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    id7F, id7V, wrong7W, "tt12_dimension_not_a_tolerance"]];
assert["tt12_dimension_does_not_relax_residual",
  FailureQ[wrong7Match] && wrong7Match["ID"] === "E6" &&
  wrong7Match["EffectiveTolerance"] ===
    DiffExp2`Tolerances`Tol["LaurentLeadTol"]];

(* TT13: a symbolic residual with an inexact coefficient is not certified
   merely because NumericQ is false; without a parameter domain it must be
   rejected rather than silently accepted. *)
symF = {{esT[0, {Global`yy, 0, 0}]}};
symV = {esT[0, {0, 0, 0}]};
symW = {esT[0, {0.01`30, 0, 0}]};
symMatch = catchDE2[
  DiffExp2`Transport`Private`matchingResidualAssert[
    symF, symV, symW, "tt13_inexact_symbolic"]];
assert["tt13_inexact_symbolic_residual_loud",
  FailureQ[symMatch] && symMatch["ID"] === "E6"];

(* TT14: matching input elimination owns a RankTol decision domain.  At high
   WP, resolved coefficients far below the fixed LaurentLeadTol structural
   floor are still rank information.  After input classification, row
   cancellation is suffix-invariant and only a certified centered zero may
   advance the formal epsilon valuation. *)
catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 500}]];
rankDelta = 10^-40;
rankOne = esT[0, {1, 0, 0, 0}];
rankZero = esT[0, {0, 0, 0, 0}];
rankOff = esT[0, {rankDelta, 1, 0, 0}];
rankF = {{rankOne, rankOff}, {rankZero, rankOne}};
rankV = {rankOff, rankOne};
rankTrimKept = DiffExp2`Transport`Private`mwTrim[
  {0, {rankDelta, 1, 0}}, "tt14_rank_keep"];
rankTrimDropped = DiffExp2`Transport`Private`mwTrim[
  {0, {10^-140, 1, 0}}, "tt14_rank_drop"];
rankComplexScale =
  4.4267459561002104836`0.2683567342206439 +
    0``-0.022790862816694443*I;
rankTrimComplex = DiffExp2`Transport`Private`mwTrim[
  {0, {10^-130, rankComplexScale}}, "tt14_complex_scale"];
polarInputTrim = DiffExp2`Transport`Private`mwInputTrim[
  {-1, {10^-40, 1, 0}}, "tt14_polar_input"];
finiteRankInput = DiffExp2`Transport`Private`mwInputTrim[
  {0, {10^-40, 1, 0}}, "tt14_finite_rank_input"];
growingSchurTrim = DiffExp2`Transport`Private`mwCancellationTrim[
  {0, {10^-40, 10^300, 0}}, "tt14_growing_schur"];
suffixLead = N[10^-40, 500];
suffixTail = N[10^300, 500];
suffixSolo = DiffExp2`Transport`Private`mwCancellationTrim[
  {2, {suffixLead}}, "tt14_suffix_solo"];
suffixGrown = DiffExp2`Transport`Private`mwCancellationTrim[
  {2, {suffixLead, suffixTail}}, "tt14_suffix_grown"];
exactCancellation = DiffExp2`Transport`Private`mwCancellationTrim[
  {2, {0, suffixLead}}, "tt14_exact_cancellation"];
uncertainCancellation = catchDE2[
  DiffExp2`Transport`Private`mwCancellationTrim[
    {2, {0``10, suffixLead}}, "tt14_uncertain_cancellation"]];
rankW = catchDE2[DiffExp2`Transport`MatchWeights[
  rankF, rankV, "tt14_rank_matching"]];
assert["tt14_matching_trim_uses_ranktol_not_laurent_floor",
  rankDelta < DiffExp2`Tolerances`Tol["LaurentLeadTol"] &&
  rankTrimKept[[1]] === 0 && rankTrimDropped[[1]] === 1 &&
  rankTrimComplex[[1]] === 1 &&
  polarInputTrim[[1]] === 0 && finiteRankInput[[1]] === 0 &&
  growingSchurTrim[[1]] === 0 &&
  suffixSolo[[1]] === 2 && suffixGrown[[1]] === 2 &&
  First[suffixGrown[[2]]] === suffixLead &&
  exactCancellation[[1]] === 3 &&
  FailureQ[uncertainCancellation] && uncertainCancellation["ID"] === "E5" &&
  !FailureQ[rankW] &&
  DiffExp2`EpsSeries`ESCoefficient[rankW[[1]], 0] === 0 &&
  DiffExp2`EpsSeries`ESCoefficient[rankW[[2]], 0] === 1];
catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100}]];

(* The value-vector transport prototype is intentionally environment-gated.
   Keep each A/B run local even when this test file is launched from a shell
   that already defines DE2_VALUE_TRANSPORT.  SetEnvironment accepts None to
   remove a variable, which faithfully restores an initially-unset state. *)
SetAttributes[withValueTransport, HoldRest];
withValueTransport[enabled_, expr_] := Module[
  {old = Quiet[Environment["DE2_VALUE_TRANSPORT"]]},
  Internal`WithLocalSettings[
    SetEnvironment["DE2_VALUE_TRANSPORT" -> If[TrueQ[enabled], "1", "0"]],
    expr,
    SetEnvironment["DE2_VALUE_TRANSPORT" -> If[StringQ[old], old, None]]]];

numericMagnitude[z_] := Max[0, Sequence @@
  (DiffExp2`Tolerances`NumericMagnitude[#, 40] & /@
    Select[Flatten[{z}], NumericQ])];
vectorESCloseQ[a_, b_, tol_] := Module[{k1, k2},
  If[!And[DiffExp2`EpsSeries`ESQ[a], DiffExp2`EpsSeries`ESQ[b]],
    Return[False, Module]];
  k1 = Min[DiffExp2`EpsSeries`ESMinPower[a],
    DiffExp2`EpsSeries`ESMinPower[b]];
  k2 = Min[DiffExp2`EpsSeries`ESCompleteMax[a],
    DiffExp2`EpsSeries`ESCompleteMax[b]];
  k2 >= k1 && AllTrue[Range[k1, k2], Function[k,
    Module[{av = DiffExp2`EpsSeries`ESCoefficient[a, k],
        bv = DiffExp2`EpsSeries`ESCoefficient[b, k], flat},
      flat = Flatten[{av, bv}];
      AllTrue[flat, NumericQ] &&
        numericMagnitude[av - bv] <=
          tol*Max[1, numericMagnitude[av], numericMagnitude[bv]]]]]];

cachedBasisCenters[] := SortBy[
  (First[#["Columns"]]["Center"] & /@
    Values[DiffExp2`Solve`Private`$shCache]), N[#, 20] &];

(* TT15: regular multi-chart A/B.  The system is genuinely epsilon
   dependent (including a rational 1/(2-eps) entry), and the boundary starts
   at eps^1.  Thus both paths must preserve the nontrivial honest window
   [1,2], not merely agree at eps^0.  A structural complex pair supplies
   regular projected waypoints; at EO100 value-aware planning should need no
   fundamental basis at all, including at the center-evaluable anchor. *)
DiffExp2`Config`UpdateConfiguration[{"ExpansionOrder" -> 100}];
sys15 = <|"Matrix" -> {{1 + eps, eps/3}, {-eps/5, 1/(2 - eps)}},
  "Variable" -> x, "SingularFactors" -> {x^2 + 1}|>;
plan15 = withValueTransport[False,
  DiffExp2`Transport`SegmentLine[sys15, {0, 3/2}]];
plan15Value = withValueTransport[True,
  DiffExp2`Transport`SegmentLine[sys15, {0, 3/2}]];
bv15 = {{0, 1, -2}, {0, 3, 1/2}};
DiffExp2`Solve`ClearSolveCaches[];
basis15 = withValueTransport[False,
  catchDE2[DiffExp2`Transport`TransportLine[sys15, bv15, plan15]]];
basisCenters15 = cachedBasisCenters[];
DiffExp2`Solve`ClearSolveCaches[];
value15 = withValueTransport[True,
  catchDE2[DiffExp2`Transport`TransportLine[sys15, bv15, plan15Value]]];
valueCenters15 = cachedBasisCenters[];
assert["tt15_value_mode_regular_multichart_exercised",
  Length[plan15["Charts"]] >= 3 &&
  Length[plan15Value["Charts"]] >= Length[plan15["Charts"]] &&
  AllTrue[plan15Value["Charts"], !TrueQ[#["Singular"]] &] &&
  Length[basisCenters15] === Length[plan15["Charts"]] &&
  valueCenters15 === {}];
assert["tt15_value_mode_preserves_epsilon_window",
  !FailureQ[basis15] && !FailureQ[value15] &&
  basis15["Value", "EpsWindow"] ===
    <|"Min" -> 1, "CompleteMax" -> 2|> &&
  value15["Value", "EpsWindow"] === basis15["Value", "EpsWindow"] &&
  value15["Final", "EpsWindow"] === basis15["Final", "EpsWindow"]];
assert["tt15_value_mode_matches_basis_mode",
  !FailureQ[basis15] && !FailureQ[value15] &&
  vectorESCloseQ[value15["Value"], basis15["Value"], 10^-12]];

(* TT16: singular-adjacent branch coverage using the squeezed-anchor chain
   from TT4.  The extra x=1/2 singularity is structural geometry rather than
   a pole of the scalar ODE, so crossing is single-valued and has an exact
   reference.  Value mode must still use a basis at the marked singular
   chart.  Value-aware planning refines any receding regular hop that would
   exceed the strict center margin, while the singular chart itself remains a
   basis solve with the same crossing prescription. *)
(* EO50 gives the production structural center margin (~0.309). *)
DiffExp2`Config`UpdateConfiguration[{"ExpansionOrder" -> 50}];
sys16 = Join[sys4, <|"ExtraSingularFactors" -> {x - 1/2}|>];
plan16Value = withValueTransport[True,
  DiffExp2`Transport`SegmentLine[sys16, {11/23, 1}]];
DiffExp2`Solve`ClearSolveCaches[];
basis16 = withValueTransport[False,
  catchDE2[DiffExp2`Transport`TransportLine[sys16, bv4, plan4]]];
basisCenters16 = cachedBasisCenters[];
DiffExp2`Solve`ClearSolveCaches[];
value16 = withValueTransport[True,
  catchDE2[DiffExp2`Transport`TransportLine[sys16, bv4, plan16Value]]];
valueCenters16 = cachedBasisCenters[];
singPos16 = FirstPosition[plan16Value["Charts"],
  c_ /; TrueQ[c["Singular"]] && c["Center"] === 1/2, Missing["NotFound"]];
assert["tt16_value_mode_singular_crossing_chain_exercised",
  singPos16 =!= Missing["NotFound"] &&
  First[singPos16] < Length[plan16Value["Charts"]] &&
  !TrueQ[plan16Value["Charts"][[First[singPos16] + 1, "Singular"]]] &&
  plan16Value["Charts"][[First[singPos16] + 1, "Center"]] > 1/2];
assert["tt16_value_mode_singular_and_margin_fallbacks",
  !FailureQ[value16] && MemberQ[valueCenters16, 1/2] &&
  Length[valueCenters16] <= Length[basisCenters16] &&
  Length[basisCenters16] === Length[plan4["Charts"]]];
assert["tt16_value_mode_matches_basis_across_crossing",
  !FailureQ[basis16] && !FailureQ[value16] &&
  value16["Value", "EpsWindow"] === basis16["Value", "EpsWindow"] &&
  vectorESCloseQ[value16["Value"], basis16["Value"], 10^-10]];
DiffExp2`Config`UpdateConfiguration[{"ExpansionOrder" -> 30}];

(* TT17: a finite-accuracy boundary can have enough digits for the requested
   result but too little headroom for another center-to-center value solve.
   N[..., 2 WP] cannot manufacture those digits.  The quality gate is
   relative for large values and absolute below unit scale, matching the ODE
   residual's Max[1, scale] convention; an insufficient handoff must choose
   the ordinary basis path before the recurrence, not fail afterward. *)
qualityHigh = {esT[0, {N[1, 40], N[10^40, 30], N[10^-40, 30]}]};
qualityLow = {esT[0, {N[1, 10], 0, 0}]};
assert["tt17_value_handoff_significance_gate_scale_aware",
  DiffExp2`Transport`Private`valueHandoffAccurateQ[qualityHigh] &&
  !DiffExp2`Transport`Private`valueHandoffAccurateQ[qualityLow]];
margin17 = DiffExp2`Transport`Private`valueCenterMargin[50];
marginTarget17 = N[(DiffExp2`Tolerances`Tol["LaurentLeadTol"]/100)^(1/51), 30];
assert["tt17_value_center_margin_uses_structural_tail_floor",
  Abs[margin17 - marginTarget17] < 10^-25 && 3/10 < margin17 < 8/25 &&
  margin17^51 <= N[DiffExp2`Tolerances`Tol["LaurentLeadTol"]/100, 20]*
    (1 + 10^-20)];

(* End-to-end witness: f1' = K f2, f2' = 0 and boundary
   f1(0)=1-K h, f2(0)=1 make f1(h)=1 by a 70-decade cancellation at the
   second chart center h.  The boundary still has 80 digits of relative
   precision and matches safely at its natural scale, while the center value
   has only about ten absolute digits and must trigger the basis choice. *)
sys17 = <|"Matrix" -> {{0, 10^70}, {0, 0}}, "Variable" -> x,
  "SingularFactors" -> {}|>;
plan17 = DiffExp2`Transport`SegmentLine[sys17, {0, 1/2}];
h17 = plan17["Charts"][[2, "Center"]];
bv17 = {{N[1 - 10^70*h17, 80], 0, 0}, {N[1, 80], 0, 0}};
DiffExp2`Solve`ClearSolveCaches[];
value17 = withValueTransport[True,
  catchDE2[DiffExp2`Transport`TransportLine[sys17, bv17, plan17]]];
valueCenters17 = cachedBasisCenters[];
assert["tt17_low_significance_value_handoff_uses_basis",
  !FailureQ[value17] && Length[plan17["Charts"]] >= 2 &&
  valueCenters17 === {h17}];

(* TT18: a genuine epsilon-adic rank defect is repaired by a column-lattice
   saturation, while an ordinary resolved but ill-conditioned constant
   matrix is not.  The model F={{1,1},{0,eps}} requires exactly the
   replacement (col2-col1)/eps and then admits regular matching weights. *)
satF18 = {{esT[0, {1, 0, 0, 0}], esT[0, {1, 0, 0, 0}]},
  {esT[0, {0, 0, 0, 0}], esT[0, {0, 1, 0, 0}]}};
satV18 = {esT[0, {1, 0, 0}], esT[0, {1, 0, 0}]};
satPlan18 = catchDE2[
  DiffExp2`Transport`Private`mwSaturationPlan[satF18, "tt18_eps_lattice"]];
satW18 = If[FailureQ[satPlan18], satPlan18,
  catchDE2[DiffExp2`Transport`MatchWeights[
    satPlan18["Matrix"], satV18, "tt18_eps_lattice"]]];
illF18 = {{esT[0, {N[1, 80], 0, 0}],
    esT[0, {N[10^70, 80], 0, 0}]},
  {esT[0, {0, 0, 0}], esT[0, {N[1, 80], 0, 0}]}};
illPlan18 = catchDE2[
  DiffExp2`Transport`Private`mwSaturationPlan[illF18, "tt18_constant_scale"]];
(* A zero component row can carry arbitrary-precision centered zeros after
   cancellation.  Its central-value scale is exactly zero; that is an
   absolute-zero case, not evidence for a spurious column valuation.  The
   old scale==0 refusal assigned col2 valuation -1, then generated the no-op
   pair col2 -> eps col2 -> (eps col2)/eps forever. *)
centeredZero18 = 0``200;
zeroRowF18 = {
  {esT[0, {1, 0, 0, 0}],
    esT[-1, ConstantArray[centeredZero18, 4]]},
  {esT[0, {0, 0, 0, 0}], esT[-1, {0, 1, 0, 0}]}};
zeroRowPlan18 = catchDE2[DiffExp2`Transport`Private`mwSaturationPlan[
  zeroRowF18, "tt18_centered_zero_row"]];
tinyResolved18 = N[10^-70, 200];
tinyTrim18 = catchDE2[DiffExp2`Transport`Private`mwInputTrim[
  {0, {tinyResolved18, 0, 0}}, "tt18_tiny_resolved"]];
(* Determinant valuation is part of the saturation certificate.  Pin the
   polynomial-time Newton implementation against the former exhaustive
   Leibniz sum at n=4, then exercise the n=16 kite-sized case that must never
   attempt 16! permutations. *)
legacyDet18[M_] := Module[{n = Length[M], terms, fac, out},
  terms = DeleteCases[Table[
    fac = Table[M[[r, perm[[r]]]], {r, n}];
    If[AnyTrue[fac, DiffExp2`Transport`Private`mwZeroQ], Nothing,
      DiffExp2`Transport`Private`mwScaleBy[Signature[perm],
        Fold[DiffExp2`Transport`Private`mwMul, First[fac], Rest[fac]]]],
    {perm, Permutations[Range[n]]}], Nothing];
  out = Fold[DiffExp2`Transport`Private`mwAdd, First[terms], Rest[terms]];
  DiffExp2`Transport`Private`mwCancellationTrim[
    out, "tt18_legacy_det", terms]];
SeedRandom[1804];
detSmallFrames18 = Table[{0, RandomInteger[{-3, 3}, 5]}, {4}, {4}];
detSmallOld18 = legacyDet18[detSmallFrames18];
detSmallNew18 = catchDE2[DiffExp2`Transport`Private`mwDetFrame[
  detSmallFrames18, "tt18_newton_det_small"]];
detLargeFrames18 = Table[If[r === c,
    If[r === 1, {1, {2, 0, 0, 0}}, {0, {1, 0, 0, 0, 0}}],
    {4, {0}}], {r, 16}, {c, 16}];
detLargeNew18 = catchDE2[DiffExp2`Transport`Private`mwDetFrame[
  detLargeFrames18, "tt18_newton_det_kite16"]];
assert["tt18_epsilon_lattice_saturates_to_regular_weights",
  !FailureQ[satPlan18] && satPlan18["InitialShifts"] === {0, 0} &&
  satPlan18["Steps"] === 1 && !FailureQ[satW18] &&
  AllTrue[satW18, DiffExp2`EpsSeries`ESMinPower[#] >= 0 &]];
assert["tt18_constant_ill_conditioning_is_not_epsilon_degeneracy",
  !FailureQ[illPlan18] && illPlan18["InitialShifts"] === {0, 0} &&
  illPlan18["Steps"] === 0];
assert["tt18_centered_zero_row_does_not_create_noop_saturation_cycle",
  !FailureQ[zeroRowPlan18] &&
  zeroRowPlan18["InitialShifts"] === {0, 0} &&
  zeroRowPlan18["Steps"] === 0];
assert["tt18_zero_scale_fix_preserves_resolved_tiny_rank_data",
  !FailureQ[tinyTrim18] && tinyTrim18[[1]] === 0 &&
  tinyTrim18[[2, 1]] === tinyResolved18];
assert["tt18_polynomial_determinant_matches_exhaustive_small_case",
  !FailureQ[detSmallNew18] && detSmallNew18 === detSmallOld18];
assert["tt18_polynomial_determinant_scales_to_kite_dimension_16",
  !FailureQ[detLargeNew18] && detLargeNew18[[1]] === 1 &&
  First[detLargeNew18[[2]]] === 2];

(* TT19: after epsilon-lattice saturation, normalize an ill-conditioned
   constant match frame by the epsilon-independent right action
   P=F_eps0^-1.  The witness combines both a 10^70 column scale and a
   10^70 cancellation.  Its essential 10^-70 entry must not be thresholded;
   the re-evaluated LocalSolutions must have identity leading values without
   changing the honest epsilon window. *)
kappa19 = N[10^70, 200];
basis19 = {mkMockColumn[0, {kappa19, 0}],
  mkMockColumn[0, {kappa19, 1}]};
basisValues19[bb_List] := Module[{ev},
  ev = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#, 0,
      "UsePade" -> False, "ImSign" -> 1]["Value"] &, bb];
  Table[DiffExp2`EpsSeries`ESNew[
      DiffExp2`EpsSeries`ESMinPower[ev[[i]]],
      Table[DiffExp2`EpsSeries`ESCoefficient[ev[[i]], k][[c]],
        {k, DiffExp2`EpsSeries`ESMinPower[ev[[i]]],
          DiffExp2`EpsSeries`ESCompleteMax[ev[[i]]]}]],
    {c, 2}, {i, 2}]];
pre19 = catchDE2[DiffExp2`Transport`Private`mwConstantMatchPrecondition[
  basis19, basisValues19[basis19], basisValues19, 2,
  "tt19_constant_frame"]];
assert["tt19_constant_match_frame_normalized_without_window_loss",
  !FailureQ[pre19] && TrueQ[pre19["Applied"]] &&
  Min[# ["EpsWindow", "CompleteMax"] & /@ pre19["Basis"]] === 2 &&
  DiffExp2`Transport`Private`mwIdentityMatrixQ[
    DiffExp2`Transport`Private`mwLeadingValueMatrix[
      pre19["Matrix"], "tt19_check"], "tt19_check"]];
lowKappa19 = N[10^70, 80];
lowBasis19 = {mkMockColumn[0, {lowKappa19, 0}],
  mkMockColumn[0, {lowKappa19, 1}]};
lowPre19 = catchDE2[
  DiffExp2`Transport`Private`mwConstantMatchPrecondition[
    lowBasis19, basisValues19[lowBasis19], basisValues19, 2,
    "tt19_underresolved_frame"]];
assert["tt19_constant_match_frame_underresolved_identity_loud",
  FailureQ[lowPre19] && lowPre19["ID"] === "E5"];

(* The crossing operator is a whole-tagged-solution pullback, not only the
   b=1 scalar phase exercised by TT5b.  Pin the opposite prescription sign,
   a fractional a-sector, and the p=2 log-chain mixing by comparing the
   transformed positive-arm value to the original negative-arm value. *)
mkCross19[a_, p_, sigma_] := <|
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 1,
  "Sectors" -> {<|"a" -> a, "b" -> 0, "p" -> p,
    "Coeffs" -> Table[{{If[k === 0, 1, 0]}}, {k, 0, 2}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 0|>, "ErrorEstimate" -> {0, 0, 0},
  "Prescriptions" -> {<|"Factor" -> x, "Sign" -> sigma,
    "Multiplicity" -> 1, "LeadingCoeffSign" -> 1|>}|>;
crossPullback19[ls_, sigma_] := {
  DiffExp2`SectorSeries`EvaluateLocalSolution[ls, -1/4,
    "UsePade" -> False, "ImSign" -> sigma]["Value"],
  DiffExp2`SectorSeries`EvaluateLocalSolution[
    DiffExp2`Transport`ApplyCrossing[ls, sigma], 1/4,
    "UsePade" -> False, "ImSign" -> sigma]["Value"]};
fracCross19 = crossPullback19[mkCross19[1/2, 0, -1], -1];
logCross19 = crossPullback19[mkCross19[0, 2, 1], 1];
crossExact19[pair_] := Module[{kmin, kmax},
  kmin = Min[DiffExp2`EpsSeries`ESMinPower /@ pair];
  kmax = Min[DiffExp2`EpsSeries`ESCompleteMax /@ pair];
  And @@ Flatten[Table[
    Together[Expand[#]] === 0 & /@
      (DiffExp2`EpsSeries`ESCoefficient[pair[[1]], k] -
        DiffExp2`EpsSeries`ESCoefficient[pair[[2]], k]),
    {k, kmin, kmax}]]];
assert["tt19_crossing_fractional_sigma_minus_and_log_chain",
  crossExact19[fracCross19] && crossExact19[logCross19]];

(* Cancellation trimming may discard only a certified centered zero.  A
   resolved tiny nonzero remains formal pivot data even below 10^-24, while
   a centered value whose uncertainty exceeds the cancellation contract is
   loud rather than becoming a pivot. *)
resolvedTiny19 = N[10^-40, 80];
underresolved19 = 0``10;
ambiguous19 = catchDE2[DiffExp2`Transport`Private`mwCenteredZeroQ[
  underresolved19, "tt19_uncertainty_overlap"]];
assert["tt19_centered_zero_gate_keeps_tiny_and_rejects_overlap",
  !DiffExp2`Transport`Private`mwCenteredZeroQ[
    resolvedTiny19, "tt19_resolved_tiny"] &&
  FailureQ[ambiguous19] && ambiguous19["ID"] === "E5"];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
