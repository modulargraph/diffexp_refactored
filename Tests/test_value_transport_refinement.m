(* Focused value-aware segmentation and anchor propagation tests. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
SetAttributes[withValueTransport, HoldRest];
withValueTransport[enabled_, expr_] := Module[
  {old = Quiet[Environment["DE2_VALUE_TRANSPORT"]]},
  Internal`WithLocalSettings[
    SetEnvironment["DE2_VALUE_TRANSPORT" -> If[TrueQ[enabled], "1", "0"]],
    expr,
    SetEnvironment["DE2_VALUE_TRANSPORT" -> If[StringQ[old], old, None]]]];

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 80,
  "ExpansionOrder" -> 60, "EpsilonOrder" -> 0, "DivisionOrder" -> 3}]];
x = Global`x; eps = Global`eps;

(* A classic two-center step has center ratio 1/2 and is unsafe under the
   EO60 value margin.  Value-aware planning bisects it, while flag-off
   planning remains the exact classic plan. *)
sys = <|"Matrix" -> {{1/(x - 2)}}, "Variable" -> x,
  "SingularFactors" -> {x - 2}|>;
planOff = withValueTransport[False,
  catchDE2[DiffExp2`Transport`SegmentLine[sys, {0, 1}]]];
planOn = withValueTransport[True,
  catchDE2[DiffExp2`Transport`SegmentLine[sys, {0, 1}]]];
centersOff = # ["Center"] & /@ planOff["Charts"];
centersOn = # ["Center"] & /@ planOn["Charts"];
margin = DiffExp2`Transport`Private`valueCenterMargin[60];

assert["vt_refine_flag_off_classic_plan_unchanged",
  !FailureQ[planOff] && centersOff === {0, 1} &&
    planOff["SegmentCount"] === 2];
assert["vt_refine_exact_centers_and_updated_count",
  !FailureQ[planOn] && centersOn === {0, 1/2, 1} &&
    FreeQ[centersOn, _?InexactNumberQ] &&
    planOn["SegmentCount"] === Length[planOn["Charts"]] === 3];
assert["vt_refine_strict_center_margin",
  AllTrue[Range[2, Length[planOn["Charts"]]], Function[i,
    Module[{a = planOn["Charts"][[i - 1]], b = planOn["Charts"][[i]]},
      TrueQ[Abs[N[b["Center"] - a["Center"], 40]] <
        margin*N[a["Radius"], 40]]]]]];
assert["vt_refine_validateplan_and_fallback_matches",
  !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[planOn]]] &&
    AllTrue[Range[2, Length[planOn["Charts"]]], Function[i,
      Module[{a = planOn["Charts"][[i - 1]], b = planOn["Charts"][[i]],
          mp = planOn["Charts"][[i, "IncomingMatchPoint"]]},
        TrueQ[Abs[N[mp - a["Center"], 40]] <= N[a["Radius"], 40]/3] &&
          TrueQ[Abs[N[mp - b["Center"], 40]] <= N[b["Radius"], 40]/3]]]]];

(* With exact center data, every regular chart including the first anchor is
   a one-vector solve.  The homogeneous basis cache therefore stays empty.
   Running the same refined plan with the flag disabled exercises its stored
   certified basis-fallback match points. *)
DiffExp2`Solve`ClearSolveCaches[];
valueRes = withValueTransport[True,
  catchDE2[DiffExp2`Transport`TransportLine[sys, {{1}}, planOn]]];
valueBasisCacheCount = Length[DiffExp2`Solve`Private`$shCache];
DiffExp2`Solve`ClearSolveCaches[];
basisRes = withValueTransport[False,
  catchDE2[DiffExp2`Transport`TransportLine[sys, {{1}}, planOn]]];
tileIntegral = catchDE2[DiffExp2`API`LineIntegral[sys, {{1}}, 0, {0, 1}, {1},
  "PrecomputedCharts" -> valueRes["Charts"]]];
assert["vt_refine_first_anchor_and_all_regular_charts_use_values",
  !FailureQ[valueRes] && valueBasisCacheCount === 0];
assert["vt_refine_basis_fallback_matches_same_solution",
  !FailureQ[basisRes] &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[valueRes["Value"], 0][[1]] -
      DiffExp2`EpsSeries`ESCoefficient[basisRes["Value"], 0][[1]], 30]] < 10^-20];
assert["vt_refine_half_radius_tiles_cover_integral",
  !FailureQ[tileIntegral] &&
    Abs[N[DiffExp2`EpsSeries`ESCoefficient[tileIntegral, 0] - 3/4, 30]] < 10^-20];

(* SegmentErrorProbe's authoritative full-vs-reduced definition is unchanged:
   only the advisory tail scans inside its two evaluations are suppressed. *)
probeLS = valueRes["Charts"][[1, "LocalSolution"]];
probePoint = probeLS["Radius"]/5;
probeDec = DiffExp2`Tolerances`EvalErrorSeriesDecrease[1];
probeFull = DiffExp2`SectorSeries`EvaluateLocalSolution[probeLS, probePoint,
  "UsePade" -> False];
probeReduced = DiffExp2`SectorSeries`EvaluateLocalSolution[probeLS, probePoint,
  "UsePade" -> False, "TOrderReduction" -> probeDec];
probeLegacy = Table[Module[
    {kf = DiffExp2`EpsSeries`ESCoefficient[probeFull["Value"], k], kr},
    kr = If[DiffExp2`EpsSeries`ESMinPower[probeReduced["Value"]] <= k <=
        DiffExp2`EpsSeries`ESCompleteMax[probeReduced["Value"]],
      DiffExp2`EpsSeries`ESCoefficient[probeReduced["Value"], k], 0*kf];
    Max[0, Sequence @@ (Last[
      DiffExp2`Transport`Private`numMagBounds[#, 20]] & /@
        Select[Flatten[{kf - kr}], NumericQ])]],
  {k, DiffExp2`EpsSeries`ESMinPower[probeFull["Value"]],
    DiffExp2`EpsSeries`ESCompleteMax[probeFull["Value"]]}];
probeFast = DiffExp2`Transport`SegmentErrorProbe[probeLS, probePoint, 0];
assert["vt_tail_fastpath_segment_probe_semantics_unchanged",
  probeFast === probeLegacy];

(* Refinement may add regular charts but must not add, remove, or reorder
   singular charts.  The existing singular FixWithin records still validate. *)
sysSing = <|"Matrix" -> {{eps/x}}, "Variable" -> x,
  "SingularFactors" -> {x}|>;
sysSingExtra = Join[sysSing, <|"ExtraSingularFactors" -> {x - 1/2}|>];
planSingOff = withValueTransport[False,
  catchDE2[DiffExp2`Transport`SegmentLine[sysSingExtra, {11/23, 1}]]];
planSingOn = withValueTransport[True,
  catchDE2[DiffExp2`Transport`SegmentLine[sysSingExtra, {11/23, 1}]]];
singularCenters[plan_] := # ["Center"] & /@
  Select[plan["Charts"], TrueQ[# ["Singular"]] &];
assert["vt_refine_singular_chain_unchanged_and_valid",
  !FailureQ[planSingOff] && !FailureQ[planSingOn] &&
    singularCenters[planSingOn] === singularCenters[planSingOff] &&
    !FailureQ[catchDE2[DiffExp2`Transport`ValidatePlan[planSingOn]]]];

Print["Value transport refinement tests: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1]];
