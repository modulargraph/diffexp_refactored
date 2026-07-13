(* A native executor may reject a regular-singular SCC for a later
   capability reason.  Its "unsupported" execution_scope must not cause the
   Wolfram request builder to reinterpret an exact p>0 Jordan seed as a
   regular eps^0 unit column. *)

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
  "TOrder" -> 2|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-singular-classifier-polar-edge",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  (* The diagonal {1,2} block has an exact p>0 Jordan seed. *)
  system = <|"Matrix" -> {
      {lambda/x, 1/x, 0},
      {1, lambda/x, 0},
      {1/x, 0, lambda/x}}, "Variable" -> x|>;
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  stats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  (* Classifier regression, independent of which coupling capabilities the
     current native build implements: retain the real collision-bound block
     metadata but make the executor admission result unsupported for this
     one call. *)
  unsupportedStats = If[AssociationQ[stats], Join[stats, <|
    "execution_implemented" -> False, "execution_scope" -> "unsupported",
    "regular_singular_scalar_block_dag_column_execution" -> False,
    "regular_singular_jordan_block_dag_column_execution" -> False|>],
    stats];
  savedStatsDownValues =
    DownValues[DiffExp2`CppBackend`PersistentSCCStatistics];
  statsWasProtected = MemberQ[
    Attributes[DiffExp2`CppBackend`PersistentSCCStatistics], Protected];
  If[statsWasProtected,
    Unprotect[DiffExp2`CppBackend`PersistentSCCStatistics]];
  DownValues[DiffExp2`CppBackend`PersistentSCCStatistics] = {};
  DiffExp2`CppBackend`PersistentSCCStatistics[_] := unsupportedStats;
  solved = If[FailureQ[prepared], prepared,
    Block[{Print = Function[Null]},
      catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
        cs, request, 1, 2]]]];
  DownValues[DiffExp2`CppBackend`PersistentSCCStatistics] =
    savedStatsDownValues;
  If[statsWasProtected,
    Protect[DiffExp2`CppBackend`PersistentSCCStatistics]];
  {cs, prepared, stats, solved}];

{cs, prepared, stats, solved} = result;
payload = If[FailureQ[solved], solved[[2]], <||>];
ok = !AnyTrue[Take[result, 3], FailureQ] && FailureQ[solved] &&
  Lookup[stats, "execution_scope", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  Lookup[payload, "ID", None] === "E6" &&
  Lookup[payload, "ExpectedCapability", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  StringContainsQ[Lookup[payload, "Detail", ""],
    "does not advertise the required strict regular or exact affine-Jordan"] &&
  !StringContainsQ[Lookup[payload, "Detail", ""],
    "not one complete exact-encoded eps^0 local basis"];

If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": unsupported native scope preserves exact singular seed classification"];
If[!ok, Print[InputForm[{
  "Stats" -> If[AssociationQ[stats],
    KeyTake[stats, {"execution_implemented", "execution_scope",
      "block_charts"}], stats],
  "Solved" -> solved}]]];
Exit[If[ok, 0, 1]];
