(* Focused Wolfram producer -> persistent C++ CASE-P compensation smoke.
   A scalar (a,b)=(1,0) seed enters a size-2 (1,1) Jordan block, creating
   eps^-2 compensation, and both exact tags then cross a third SCC block. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps;
request = <|"EpsWindow" -> <|"Min" -> -2, "CompleteMax" -> 0|>,
  "TOrder" -> 0|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-scc-pseudo-compensation",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;
system = <|"Matrix" -> {
    {1/x, 0, 0, 0},
    {0, (1 + eps)/x, 1/x, 0},
    {1/x, 1, (1 + eps)/x, 0},
    {0, 1/x, 0, (3/2 + 3 eps)/x}}, "Variable" -> x|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  stats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  solved = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 1]]];
  evaluated = If[FailureQ[solved], solved,
    DiffExp2`CppBackend`EvaluatePersistentLocal[solved,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  {cs, prepared, stats, solved, evaluated}];

{cs, prepared, stats, solved, evaluated} = result;
diagnostics = If[AssociationQ[solved],
  Lookup[Lookup[solved, "NativeSummary", <||>],
    "block_diagnostics", {}], {}];
collisionRecords = Select[diagnostics,
  AssociationQ[#] && Lookup[#, "block", None] === 1 &&
    Lookup[#, "source_b", None] === "0" &];
tailTags = Sort@DeleteDuplicates[Lookup[
  Select[diagnostics, AssociationQ[#] &&
    Lookup[#, "block", None] === 2 &], "source_b", {}]];
decoded = If[AssociationQ[evaluated] &&
    Lookup[evaluated, "status", "error"] === "ok",
  DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60], {}];

collisionOk = Length[collisionRecords] === 1 && With[
  {record = First[collisionRecords]},
  Lookup[record, "pseudo_hit_count", None] === 1 &&
    Lookup[record, "pseudo_compensation_count", None] === 2 &&
    Lookup[record, "max_pseudo_depth", None] === 2 &&
    TrueQ[Lookup[record, "pseudo_value_certified", False]] &&
    Lookup[record, "uncompensated_pseudo_hit_count", None] === 0];

ok = !AnyTrue[result, FailureQ] &&
  Lookup[cs["IntegrationSequence"], "Components", None] ===
    {{1}, {2, 3}, {4}} &&
  Lookup[stats, "execution_scope", None] ===
    "exact-rational-regular-singular-jordan-block-dag-column-v2" &&
  Lookup[stats["capability_evidence"],
    "pseudo_schedule_execution", None] ===
      "exact-rational-joint-compensation-and-formal-overlap-certificate" &&
  Lookup[solved["NativeSummary"], "pseudo_hit_count", None] === 0 &&
  collisionOk && tailTags === {"0", "1"} &&
  Lookup[evaluated, "status", "error"] === "ok" &&
  Lookup[evaluated["value"], "min", None] === -2 &&
  Length[decoded] >= 8 &&
  AllTrue[Take[decoded, 8],
    !FailureQ[#] && NumericQ[#] && Abs[N[#, 50]] < 10^-30 &];

If[AssociationQ[solved] && StringQ[Lookup[solved, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[solved]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": Wolfram-wired exact SCC CASE-P compensation"];
If[!ok, Print[InputForm[result]]];
Exit[If[ok, 0, 1]];
