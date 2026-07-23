(* Focused persistent-native lower/upper arm orchestration test. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 100,
  "ChopPrecision" -> 50,
  "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2,
  "DivisionOrder" -> 3,
  "Verbosity" -> 0}];

sys = DiffExp2`API`LoadSystem[<|
  "Matrix" -> {{0}}, "Variable" -> Global`x|>];

result = Catch[
  DiffExp2`NativeTransport`NativeRegularLineIntegral[
    sys, {{1, 0, 0, 0}}, 0, {-1, 1},
    {1/Global`eps + Global`x^2}, "Threads" -> 10],
  "DiffExp2Error"];

released = If[AssociationQ[result],
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    result["Run"]], <||>];
sessionStats = If[AssociationQ[result],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    result["Atlas", "Session"], <||>], <||>];

positiveResult = Catch[
  DiffExp2`NativeTransport`NativeRegularLineIntegral[
    sys, {{0, 1, 0}}, 0, {-1, 1},
    {Global`eps^5}, "Threads" -> 10],
  "DiffExp2Error"];
positiveReleased = If[AssociationQ[positiveResult],
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    positiveResult["Run"]], <||>];
positiveSessionStats = If[AssociationQ[positiveResult],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    positiveResult["Atlas", "Session"], <||>], <||>];

ok = !FailureQ[result] && AssociationQ[result] &&
  Lookup[result, "Type", None] ===
    "DiffExp2NativeRegularLineIntegral" &&
  TrueQ[Abs[N[
      DiffExp2`EpsSeries`ESCoefficient[result["Value"], -1] - 2,
      40]] < 10^-30] &&
  TrueQ[Abs[N[
      DiffExp2`EpsSeries`ESCoefficient[result["Value"], 0] - 2/3,
      40]] < 10^-30] &&
  result["Value", "EpsWindow"] ===
    <|"Min" -> -1, "CompleteMax" -> 2|> &&
  result["Atlas", "TargetCompleteMax"] === 2 &&
  result["Atlas", "Request", "EpsWindow", "CompleteMax"] === 3 &&
  result["Atlas", "PreparedIntegrandEpsilonShift"] === -1 &&
  result["Atlas", "RegularValueAggregationGuardDigits"] ===
    Ceiling[Log10[Max[
      Length[result["Atlas", "Lower", "Plan", "Tiles"]] +
      Length[result["Atlas", "Upper", "Plan", "Tiles"]], 1]]] &&
  result["Run", "Lower", "Matches"] > 0 &&
  result["Run", "Upper", "Matches"] > 0 &&
  result["Run", "Lower", "Tiles"] ===
    Length[result["Atlas", "Lower", "ChartSystems"]] &&
  result["Run", "Upper", "Tiles"] ===
    Length[result["Atlas", "Upper", "ChartSystems"]] &&
  result["CompatibilityExports"] === 1 &&
  Lookup[result["Run", "NativeSummary"], "json_coefficients", None] === 0 &&
  TrueQ[Lookup[result["Run", "NativeSummary"],
    "atomic_publication", False]] &&
  TrueQ[Lookup[result["Run", "NativeSummary"], "worker_overlap", False]] &&
  Lookup[released, "Failures", {"missing"}] === {} &&
  Lookup[sessionStats, "locals", -1] === 0 &&
  Lookup[sessionStats, "matches", -1] === 0 &&
  Lookup[sessionStats, "tile_plans", -1] === 0 &&
  Lookup[sessionStats, "line_results", -1] === 0 &&
  AssociationQ[positiveResult] &&
  positiveResult["Value"] ===
    DiffExp2`EpsSeries`ESZero[2] &&
  positiveResult["Atlas", "PreparedIntegrandEpsilonShift"] === 5 &&
  positiveResult["Atlas", "Request", "EpsWindow", "CompleteMax"] === 2 &&
  Lookup[positiveReleased, "Failures", {"missing"}] === {} &&
  Lookup[positiveSessionStats, "locals", -1] === 0 &&
  Lookup[positiveSessionStats, "matches", -1] === 0 &&
  Lookup[positiveSessionStats, "tile_plans", -1] === 0 &&
  Lookup[positiveSessionStats, "line_results", -1] === 0;

DiffExp2`CppBackend`ClearPersistentSessions[];

If[TrueQ[ok],
  Print["PASS: persistent native regular independent arms"],
  Print["FAIL: negative=", InputForm[If[FailureQ[result], result,
      KeyTake[result, {"Type", "Value", "CompatibilityExports"}]]],
    "; positive=", InputForm[If[FailureQ[positiveResult], positiveResult,
      KeyTake[positiveResult,
        {"Type", "Value", "CompatibilityExports"}]]]];
  Exit[1]];
