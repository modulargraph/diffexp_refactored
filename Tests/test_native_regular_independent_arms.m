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
    sys, {{1, 0, 0}}, 0, {-1, 1},
    {1 + Global`eps*Global`x^2}, "Threads" -> 10],
  "DiffExp2Error"];

ok = !FailureQ[result] && AssociationQ[result] &&
  Lookup[result, "Type", None] ===
    "DiffExp2NativeRegularLineIntegral" &&
  TrueQ[Abs[N[
      DiffExp2`EpsSeries`ESCoefficient[result["Value"], 0] - 2,
      40]] < 10^-30] &&
  TrueQ[Abs[N[
      DiffExp2`EpsSeries`ESCoefficient[result["Value"], 1] - 2/3,
      40]] < 10^-30] &&
  Length[result["Run", "Lower", "Matches"]] > 0 &&
  Length[result["Run", "Upper", "Matches"]] > 0 &&
  result["CompatibilityExports"] ===
    Length[result["Run", "Lower", "Lines"]] +
      Length[result["Run", "Upper", "Lines"]] &&
  AllTrue[Join[result["Run", "Lower", "ProjectedLocals"],
      result["Run", "Upper", "ProjectedLocals"]],
    Lookup[#, "json_coefficients", None] === 0 &];

DiffExp2`CppBackend`ClearPersistentSessions[];

If[TrueQ[ok],
  Print["PASS: persistent native regular independent arms"],
  Print["FAIL: ", InputForm[result]];
  Exit[1]];
