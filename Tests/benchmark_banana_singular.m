(* Opt-in order-50 timing gate for the two singular banana endpoint bases.
   The full production gate invokes this before starting its wall timer so a
   future failure identifies endpoint/Fuchsian recurrence regressions without
   changing the five-minute end-to-end contract. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["FAIL: compiled DiffExp2 backend is not available"];
  Exit[1]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
req[min_, max_, n_] := <|"EpsWindow" -> <|"Min" -> min,
  "CompleteMax" -> max|>, "TOrder" -> n|>;

limitText = Environment["BANANA_SINGULAR_MAX_SECONDS"];
limit = If[StringQ[limitText] &&
    StringMatchQ[StringTrim[limitText], DigitCharacter ..],
  ToExpression[StringTrim[limitText]], 60];
If[!IntegerQ[limit] || limit < 1,
  Print["FAIL: BANANA_SINGULAR_MAX_SECONDS must be a positive integer"];
  Exit[2]];

config = catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 500,
  "LinearSolveChopPrecision" -> 20,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 5,
  "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3,
  "RecurrenceBackend" -> "Cpp",
  "Variables" -> {}}]];
If[FailureQ[config], Print["FAIL: configuration: ", config]; Exit[1]];

fixture = Get[FileNameJoin[
  {repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
system = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fixture["Matrix"],
    "Variable" -> fixture["Variable"]|>]];
If[FailureQ[system], Print["FAIL: fixture system: ", system]; Exit[1]];

solveEndpoint[endpoint_] := Module[{plan, chart, prepared, result, seconds},
  plan = catchDE2[DiffExp2`Transport`SegmentLine[
    system, {11/23, endpoint}]];
  If[FailureQ[plan], Return[{plan, Infinity}]];
  chart = SelectFirst[Reverse[plan["Charts"]],
    TrueQ[#1["Singular"]] &, Missing["NoSingularEndpoint"]];
  If[MissingQ[chart], Return[{chart, Infinity}]];
  prepared = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  If[FailureQ[prepared], Return[{prepared, Infinity}]];
  seconds = First@AbsoluteTiming[
    result = catchDE2[DiffExp2`Solve`SolveChart[
      prepared, req[0, 5, 50]]]];
  {result, seconds}
];

{lower, lowerSeconds} = solveEndpoint[0];
DiffExp2`Solve`ClearSolveCaches[];
{upper, upperSeconds} = solveEndpoint[1];
totalSeconds = lowerSeconds + upperSeconds;

valid[result_] := AssociationQ[result] &&
  AssociationQ[Lookup[result, "Basis", None]] &&
  Length[Lookup[result["Basis"], "Columns", {}]] === 7;

Print["banana singular order-50 seconds: lower=",
  N[lowerSeconds, 6], " upper=", N[upperSeconds, 6],
  " total=", N[totalSeconds, 6], " ceiling=", limit];
If[valid[lower] && valid[upper] && totalSeconds <= limit,
  Print["PASS: banana singular endpoint timing regression"];
  Exit[0],
  If[!valid[lower], Print["lower failure: ", lower]];
  If[!valid[upper], Print["upper failure: ", upper]];
  If[NumericQ[totalSeconds] && totalSeconds > limit,
    Print["timing failure: ", N[totalSeconds, 6],
      " seconds exceeds ", limit]];
  Exit[1]];
