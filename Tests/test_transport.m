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

(* TT3: regular line with a far singularity: f' = f/(x-2): f = c(x-2);
   f(0) = 1/2 -> c = -1/4 -> f(1) = 1/4 *)
sys3 = <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x, "SingularFactors" -> {x - 2}|>;
plan3 = DiffExp2`Transport`SegmentLine[sys3, {0, 1}];
res3 = catchDE2[DiffExp2`Transport`TransportLine[sys3, {{1/2, 0, 0}}, plan3]];
assert["tt3_far_singularity",
  !FailureQ[res3] &&
  Abs[N[DiffExp2`EpsSeries`ESCoefficient[res3["Value"], 0][[1]] - 1/4, 30]] < 10^-25];

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
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[plan4]]]];

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

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
