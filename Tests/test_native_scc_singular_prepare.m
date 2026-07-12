(* Focused Wolfram-to-native smoke for singular Jordan SCC preparation.
   The regular flag is a chart classification: a non-regular diagonal block
   is admissible when the retained C++ operator proves its exact affine
   Jordan indicial form. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps;
lambda = 1/2 + eps/3;
request = <|"EpsWindow" -> <|"Min" -> -3, "CompleteMax" -> 0|>,
  "TOrder" -> 0|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-singular-jordan-scc-prepare",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  system = <|"Matrix" -> {
      {lambda/x, 1/x, 0},
      {1, lambda/x, 0},
      {1/x, 0, lambda/x}}, "Variable" -> x|>;
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  stats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  solved = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 1, 2]]];
  evaluated = If[FailureQ[solved], solved,
    DiffExp2`CppBackend`EvaluatePersistentLocal[solved,
      <|"exact" -> "1/4"|>, <|"tail_estimate" -> False|>, 60]];
  {cs, prepared, stats, solved, evaluated}];

{cs, prepared, stats, solved, evaluated} = result;
blocks = If[AssociationQ[stats], Lookup[stats, "block_charts", {}], {}];
firstCertificate = If[ListQ[blocks] && Length[blocks] >= 1 &&
    AssociationQ[blocks[[1]]],
  Lookup[blocks[[1]], "exact_affine_jordan_indicial", <||>], <||>];

ok = !AnyTrue[{cs, prepared, stats, solved, evaluated}, FailureQ] &&
  AssociationQ[cs] && AssociationQ[prepared] && AssociationQ[stats] &&
  Lookup[cs["IntegrationSequence"], "Components", None] ===
    {{1, 2}, {3}} &&
  Lookup[stats, "status", "error"] === "ok" &&
  TrueQ[Lookup[stats, "execution_implemented", False]] &&
  Lookup[stats, "execution_scope", None] ===
    "exact-rational-regular-singular-jordan-block-dag-column-v2" &&
  Lookup[firstCertificate, "dimension", None] === 2 &&
  Lookup[firstCertificate, "max_jordan_size", None] === 2 &&
  Lookup[solved, "Type", None] === "DiffExp2NativeSCCBasisColumn" &&
  Lookup[solved, "SeedLocalComponent", None] === 2 &&
  Lookup[solved, "BasisIndex", None] === 2 &&
  Lookup[solved["NativeSummary"], "execution_capability", None] ===
    "exact-rational-regular-singular-jordan-block-dag-column-v2" &&
  Lookup[evaluated, "status", "error"] === "ok" &&
  Lookup[evaluated["value"], "dimension", None] === 3;

If[AssociationQ[solved] && StringQ[Lookup[solved, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[solved]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

(* Preparation itself remains field-generic. Exact column-plan construction
   must not inspect Acb midpoint strings as a Rational unit lattice. *)
acbResult = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  acbSystem = <|"Matrix" -> {
      {lambda/x, 1/x, 0},
      {1, lambda/x, 0},
      {1/x, 0, lambda/x}}, "Variable" -> x|>;
  acbCs = catchDE2[DiffExp2`Solve`PrepareChart[acbSystem, chart]];
  acbPrepared = If[FailureQ[acbCs], acbCs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[
      acbCs, request]]];
  acbStats = If[FailureQ[acbPrepared], acbPrepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[acbPrepared]];
  {acbCs, acbPrepared, acbStats}];
{acbCs, acbPrepared, acbStats} = acbResult;
acbOk = !AnyTrue[acbResult, FailureQ] &&
  Lookup[acbStats, "status", "error"] === "ok" &&
  Lookup[acbStats, "execution_implemented", True] === False;
If[AssociationQ[acbPrepared] &&
    StringQ[Lookup[acbPrepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[acbPrepared]]];
DiffExp2`Solve`ClearSolveCaches[];
ok = TrueQ[ok] && TrueQ[acbOk];

Print[If[ok, "PASS", "FAIL"],
  ": native singular Jordan SCC preparation and column execution"];
If[!ok, Print[InputForm[{
  cs, prepared, stats, solved, evaluated, acbResult}]]];
Exit[If[ok, 0, 1]];
