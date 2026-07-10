(* A/B benchmark for the dense spectral V transform on the real banana L1
   singular endpoints.  One Wolfram kernel, cold Solve caches per mode. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

env[name_, default_] := Module[{v = Environment[name]},
  If[StringQ[v] && StringLength[StringTrim[v]] > 0, v, default]];
orders = ToExpression["{" <> env["GST_ORDERS", "10,20"] <> "}"];
wp = ToExpression[env["GST_WORKING_PRECISION", "100"]];
epsOrder = ToExpression[env["GST_EPS_ORDER", "2"]];
endpointNames = StringTrim /@ StringSplit[
  env["GST_ENDPOINTS", "upper"], ","];
mode = env["GST_MODE", "both"];
reps = ToExpression[env["GST_REPS", "2"]];
warmup = ToExpression[env["GST_WARMUP", "1"]];
If[!AllTrue[orders, IntegerQ[#] && # >= 1 &] ||
    !AllTrue[endpointNames, MemberQ[{"lower", "upper"}, #] &] ||
    !MemberQ[{"both", "legacy", "grouped"}, mode] ||
    !IntegerQ[reps] || reps < 1 || !MemberQ[{0, 1}, warmup],
  Print["invalid GST_ORDERS, GST_ENDPOINTS, GST_MODE, GST_REPS, or GST_WARMUP"];
  Exit[2]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[e_] := Quiet[Catch[e, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> wp, "ExpansionOrder" -> First[orders],
  "EpsilonOrder" -> epsOrder, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];
fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
If[FailureQ[sys], Print["load failed: ", sys]; Exit[1]];

endpointData[name_] := Module[{to, plan, chart, cs},
  to = If[name === "lower", 0, 1];
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, to}]];
  chart = SelectFirst[Reverse[plan["Charts"]], TrueQ[#["Singular"]] &];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart]];
  <|"Name" -> name, "Chart" -> chart, "ChartSystem" -> cs|>];
endpoints = endpointData /@ endpointNames;

sameSolutionQ[a_, b_] := AssociationQ[a] && AssociationQ[b] &&
  !FailureQ[a] && !FailureQ[b] &&
  a["Basis", "Columns"] === b["Basis", "Columns"] &&
  a["Basis", "Specs"] === b["Basis", "Specs"] &&
  a["Basis", "Diagnostics"] === b["Basis", "Diagnostics"];

structureSignature[s_] := If[!AssociationQ[s] || FailureQ[s], $Failed,
  <|"Specs" -> s["Basis", "Specs"],
    "Columns" -> Map[Function[col, <|
      "EpsWindow" -> col["EpsWindow"], "TWindow" -> col["TWindow"],
      "Tags" -> Lookup[col["Sectors"], {"a", "b", "p"}],
      "Dimensions" -> (Dimensions[#] & /@ Lookup[col["Sectors"], "Coeffs"])|>],
      s["Basis", "Columns"]]|>];

coefficientVector[s_] := Flatten[
  Lookup[#, "Coeffs"] & /@ Flatten[Lookup[
    s["Basis", "Columns"], "Sectors"], 1]];

numericDifference[a_, b_] := Module[{aa, bb, diffs, scale, digits},
  If[sameSolutionQ[a, b], Return[0]];
  aa = coefficientVector[a]; bb = coefficientVector[b];
  If[Length[aa] =!= Length[bb], Return[Infinity]];
  digits = Min[50, wp];
  diffs = Quiet[Check[N[aa - bb, digits], $Failed]];
  If[diffs === $Failed || !AllTrue[diffs, NumericQ], Return[Infinity]];
  scale = Max[1, Sequence @@ (Abs /@ Quiet[N[Join[aa, bb], digits]])];
  Max[0, Sequence @@ (Abs /@ diffs)]/scale];

runMode[which_String, cs_, req_] := Module[{result, seconds},
  DiffExp2`Solve`ClearSolveCaches[];
  seconds = First[AbsoluteTiming[result = Block[{
      DiffExp2`Solve`Private`$disableGroupedSpectralTransform =
        (which === "legacy")},
    catchDE2[DiffExp2`Solve`SolveChart[cs, req]]]]];
  {seconds, result}];

Do[
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> wp, "ExpansionOrder" -> n,
    "EpsilonOrder" -> epsOrder, "DivisionOrder" -> 3,
    "StepDivisionOrder" -> 3, "Variables" -> {}}]];
  req = <|"EpsWindow" -> <|"Min" -> 0,
      "CompleteMax" -> epsOrder|>, "TOrder" -> n|>;
  Do[
    cs = ep["ChartSystem"];
    legacy = Missing["NotRun"]; grouped = Missing["NotRun"];
    legacyTimes = {}; groupedTimes = {};
    If[warmup === 1,
      Scan[(runMode[#, cs, req];) &,
        If[mode === "both", {"legacy", "grouped"}, {mode}]]];
    runOrder = Which[
      mode === "both", Flatten[Table[
        If[OddQ[r], {"legacy", "grouped"}, {"grouped", "legacy"}],
        {r, reps}]],
      True, ConstantArray[mode, reps]];
    Do[
      measured = runMode[which, cs, req];
      If[which === "legacy",
        AppendTo[legacyTimes, measured[[1]]]; legacy = measured[[2]],
        AppendTo[groupedTimes, measured[[1]]]; grouped = measured[[2]]],
      {which, runOrder}];
    legacyTime = If[legacyTimes === {}, Missing["NotRun"], Median[legacyTimes]];
    groupedTime = If[groupedTimes === {}, Missing["NotRun"], Median[groupedTimes]];
    reference = If[AssociationQ[grouped], grouped, legacy];
    diags = If[!AssociationQ[reference], {},
      reference["Basis", "Diagnostics", "AdaptiveLowerFrames"]];
    frames = If[diags === {}, {},
      DeleteDuplicates[Lookup[#, {"FrameBase", "FrameWidth"}] & /@ diags]];
    stats = If[frames === {}, {}, Map[Function[fw,
      DiffExp2`Solve`Private`prepareFramedMatrix[
        cs["V"], DiffExp2`Config`CanonicalEps[], fw[[1]], fw[[2]], cs][
          "Stats"]], frames]];
    Print["SPECTRAL_BENCH ", <|
      "Endpoint" -> ep["Name"], "TOrder" -> n,
      "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
      "Repetitions" -> reps, "Warmup" -> warmup,
      "LegacySamples" -> legacyTimes, "GroupedSamples" -> groupedTimes,
      "LegacySeconds" -> legacyTime, "GroupedSeconds" -> groupedTime,
      "Speedup" -> If[NumericQ[legacyTime] && NumericQ[groupedTime],
        legacyTime/groupedTime, Missing["NotAvailable"]],
      "PairedSpeedups" -> If[mode === "both",
        MapThread[#1/#2 &, {legacyTimes, groupedTimes}],
        Missing["NotAvailable"]],
      "ExactParity" -> If[mode === "both",
        sameSolutionQ[legacy, grouped], Missing["NotAvailable"]],
      "StructuralParity" -> If[mode === "both",
        structureSignature[legacy] === structureSignature[grouped],
        Missing["NotAvailable"]],
      "MaxRelativeCoefficientDifference" -> If[mode === "both",
        numericDifference[legacy, grouped], Missing["NotAvailable"]],
      "Frames" -> frames, "PreparedStats" -> stats|>];
    If[mode === "both" &&
        (!TrueQ[structureSignature[legacy] === structureSignature[grouped]] ||
          !TrueQ[numericDifference[legacy, grouped] <=
            DiffExp2`Tolerances`Tol["ResidTol"]]), Exit[1]],
    {ep, endpoints}],
  {n, orders}];

Exit[0];
