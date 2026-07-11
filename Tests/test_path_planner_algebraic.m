(* Exact algebraic path-planner regression and performance tests. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

x = Global`x;
DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 250, "ChopPrecision" -> 25,
  "ExpansionOrder" -> 50, "EpsilonOrder" -> 4,
  "DivisionOrder" -> 3, "RadiusOfConvergence" -> 10}];

(* This is the degree-eight Landau factor on the Euclidean unequal-mass
   banana deformation {1,1,1,1}->{2,3/2,4/3,1}, p^2=-1.  Before the
   numeric-key/exact-confirmation planner, ProjectComplexRoots plus repeated
   exact Abs[Root-center] reductions took 268.8 seconds for this system. *)
landau8 = 17845920000 + 33872256000*x + 31460140800*x^2 +
  18746691840*x^3 + 7155462240*x^4 + 1701573696*x^5 +
  234203760*x^6 + 6424176*x^7 + 279841*x^8;
bananaFactors = {3 + x, 30 + 11*x, landau8, 2 + x, 1 + x};
bananaSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> bananaFactors|>;

(* Independent numerical transcription of the old DiffExp projection rule.
   It is only an oracle for the projected values; production roots and
   membership remain exact. *)
numericProjectionReference[roots_List] := Module[
  {rows, raw, ordered},
  rows = ({Re[#], Abs[Im[#]]} &) /@ N[roots, 100];
  raw = Flatten[Map[Function[row, Module[
      {re = row[[1]], h = row[[2]], leftOccupied, rightOccupied},
      If[TrueQ[h == 0], {re},
        leftOccupied = AnyTrue[rows,
          TrueQ[re - h < #[[1]] < re] &];
        rightOccupied = AnyTrue[rows,
          TrueQ[re < #[[1]] < re + h] &];
        Join[If[leftOccupied, {}, {re - h}], {re},
          If[rightOccupied, {}, {re + h}]]]]], rows]];
  ordered = Sort[raw];
  Fold[If[#1 === {} || Abs[#2 - Last[#1]] > 10^-70,
      Append[#1, #2], #1] &, {}, ordered]];

findTime = First@AbsoluteTiming[
  bananaSings = TimeConstrained[
    DiffExp2`Transport`FindSingularities[bananaSys], 10, $Failed]];
If[bananaSings === $Failed,
  assert["planner_degree8_projection_fast", False];
  Print["Path planner algebraic tests: timed out in FindSingularities"];
  Exit[1]];
refProjected = numericProjectionReference[bananaSings["All"]];
gotProjected = Re[N[bananaSings["Projected"], 90]];
assert["planner_exact_roots_and_real_subset",
  Length[bananaSings["All"]] === 12 &&
  Length[bananaSings["Real"]] === 4 &&
  FreeQ[bananaSings["All"], _?InexactNumberQ]];
assert["planner_complex_projection_plus_minus_im_parity",
  Length[gotProjected] === Length[refProjected] === 12 &&
  Max[Abs[gotProjected - N[refProjected, 90]]] < 10^-70 &&
  FreeQ[bananaSings["Projected"], _?InexactNumberQ]];
assert["planner_degree8_projection_fast",
  findTime < 5];

SetEnvironment["DE2_VALUE_TRANSPORT" -> "1"];
planTime = First@AbsoluteTiming[
  bananaPlan = TimeConstrained[
    DiffExp2`Transport`SegmentLine[bananaSys, {0, 1}], 10, $Failed]];
If[bananaPlan === $Failed,
  assert["planner_unequal_banana_total_fast", False];
  Print["Path planner algebraic tests: timed out in SegmentLine"];
  Exit[1]];
assert["planner_unequal_banana_geometry_regression",
  bananaPlan["SegmentCount"] === 3 &&
  Count[bananaPlan["Charts"], c_ /; TrueQ[c["Singular"]]] === 0 &&
  Length[bananaPlan["Singularities", "ProjectionWaypoints"]] === 12 &&
  !FailureQ[Catch[
    DiffExp2`Transport`ValidatePlan[bananaPlan], "DiffExp2Error"]]];
assert["planner_unequal_banana_total_fast",
  planTime < 5];

SetEnvironment["DE2_VALUE_TRANSPORT" -> "0"];
classicPlan = DiffExp2`Transport`SegmentLine[
  <|"Matrix" -> {{0}}, "Variable" -> x,
    "SingularFactors" -> {x^2 + 1}|>, {-3/2, 3/2}];
SetEnvironment["DE2_VALUE_TRANSPORT" -> "1"];
symPairs = Select[Range[2, Length[classicPlan["Charts"]]],
  TrueQ[Lookup[classicPlan["Charts"][[#]], "SymmetricMatch", False]] &];
assert["planner_classic_coupled_plus_minus_one_third",
  symPairs =!= {} && AllTrue[symPairs, Function[i, Module[
    {a = classicPlan["Charts"][[i - 1]], b = classicPlan["Charts"][[i]],
      mp = classicPlan["Charts"][[i, "IncomingMatchPoint"]]},
    Together[(mp - a["Center"])/a["MatchRadius"]] === 1/3 &&
    Together[(mp - b["Center"])/b["MatchRadius"]] === -1/3]]]];

(* Distinct roots far below the numeric prefilter scale must remain distinct:
   a numeric key only selects candidates, exact reduction owns merging. *)
nearSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x, 10^80*x - 1}|>;
nearSings = DiffExp2`Transport`FindSingularities[nearSys];
assert["planner_nearby_exact_roots_never_numeric_merged",
  Length[nearSings["All"]] === 2 &&
  Sort[nearSings["All"]] === {0, 10^-80}];
nearRouteSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x - 1, 10^80*(x - 1) - 1}|>;
SetEnvironment["DE2_VALUE_TRANSPORT" -> "0"];
nearRoute = TimeConstrained[
  DiffExp2`Transport`SegmentLine[nearRouteSys,
    {1 - 10^-80, 1 + 2*10^-80}], 10, $Failed];
SetEnvironment["DE2_VALUE_TRANSPORT" -> "1"];
nearCenters = If[AssociationQ[nearRoute],
  Select[nearRoute["Charts"], TrueQ[# ["Singular"]] &][[All, "Center"]],
  {}];
assert["planner_segmentline_orders_ultraclose_exact_poles",
  AssociationQ[nearRoute] && nearCenters === {1, 1 + 10^-80} &&
    !FailureQ[Catch[
      DiffExp2`Transport`ValidatePlan[nearRoute], "DiffExp2Error"]]];
closeAlgebraicDistance = DiffExp2`Transport`Private`numericDistance[
  Sqrt[2 + 10^-200], Sqrt[2], 40];
closeAlgebraicReference = N[
  10^-200/(Sqrt[2 + 10^-200] + Sqrt[2]), 60];
assert["planner_close_algebraic_distance_exact_fallback",
  NumericQ[closeAlgebraicDistance] && closeAlgebraicDistance > 0 &&
  Abs[N[closeAlgebraicDistance/closeAlgebraicReference - 1, 30]] < 10^-35];

(* Conversely, two exact representations of the same algebraic pole merge
   only after their numeric keys nominate an exact RootReduce comparison. *)
dupSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x^2 - 2, x - Sqrt[2]}|>;
dupSings = DiffExp2`Transport`FindSingularities[dupSys];
assert["planner_algebraic_duplicate_exactly_confirmed",
  Length[dupSings["All"]] === 2 &&
  Count[dupSings["All"], r_ /;
    TrueQ[PossibleZeroQ[RootReduce[r - Sqrt[2]]]]] === 1];

(* Expensive Root moduli may be numerical, but are evaluated from exact
   endpoints at WorkingPrecision.  Cheap quadratic geometry stays exact. *)
landauOnly = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {landau8}|>;
landauSings = DiffExp2`Transport`FindSingularities[landauOnly];
landauRadius = DiffExp2`Transport`ChartRadius[0, landauSings["All"]];
landauRadiusRef = Min[Abs[N[landauSings["All"], 270]]];
assert["planner_high_degree_radius_wp_numeric",
  InexactNumberQ[landauRadius] && Precision[landauRadius] >= 240 &&
  Abs[N[landauRadius - landauRadiusRef, 220]] < 10^-210];
quadSings = DiffExp2`Transport`FindSingularities[
  <|"Matrix" -> {{0}}, "Variable" -> x,
    "SingularFactors" -> {x^2 + 1}|>];
assert["planner_cheap_quadratic_radius_remains_exact",
  DiffExp2`Transport`ChartRadius[0, quadSings["All"]] === 1];

(* A high-degree real algebraic endpoint must still be recognized EXACTLY,
   excluded from its own radius alphabet, and installed as a singular chart.
   The numeric-radius optimization is not an endpoint-membership shortcut. *)
quintic = x^5 - x - 1;
quinticSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {quintic}|>;
quinticSings = DiffExp2`Transport`FindSingularities[quinticSys];
quinticEndpoint = First[quinticSings["Real"]];
quinticRadius = DiffExp2`Transport`ChartRadius[
  quinticEndpoint, quinticSings["All"]];
quinticPlan = DiffExp2`Transport`SegmentLine[
  quinticSys, {0, quinticEndpoint}];
assert["planner_high_degree_singular_endpoint_and_radius",
  quinticPlan["EndpointIsSingular"] &&
  TrueQ[Last[quinticPlan["Charts"]]["Singular"]] &&
  Last[quinticPlan["Charts"]]["Center"] === quinticEndpoint &&
  NumericQ[quinticRadius] && TrueQ[quinticRadius > 0] &&
  !FailureQ[Catch[
    DiffExp2`Transport`ValidatePlan[quinticPlan], "DiffExp2Error"]]];

(* Analytic-continuation metadata is orthogonal to the geometry speed path:
   exact real singular charts still receive the configured +/-Im side. *)
DiffExp2`Config`UpdateConfiguration[{
  "DeltaPrescriptions" -> {{x - 1/2, -1}}}];
prescribedSys = <|"Matrix" -> {{0}}, "Variable" -> x,
  "SingularFactors" -> {x - 1/2}|>;
prescribedPlan = DiffExp2`Transport`SegmentLine[prescribedSys, {0, 1}];
halfChart = SelectFirst[prescribedPlan["Charts"],
  TrueQ[# ["Singular"]] && # ["Center"] === 1/2 &];
assert["planner_delta_prescription_sign_preserved",
  AssociationQ[halfChart] && Length[halfChart["Prescriptions"]] === 1 &&
  DiffExp2`SectorSeries`ChartImSign[halfChart] === -1];

Print["Path planner algebraic tests: ", passed, " passed, ", failed,
  " failed. findTime=", N[findTime, 6], " planTime=", N[planTime, 6]];
If[failed > 0, Exit[1], Exit[0]];
