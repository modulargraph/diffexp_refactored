(* One focused end-to-end check for the production numeric coefficient field:
   a regular Acb SCC column must execute and remain retained, rather than
   stopping after composite preparation or returning coefficient slabs. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

decodeValue[evaluated_] := Module[
  {decoded, min, max, dimension},
  If[!AssociationQ[evaluated] ||
      Lookup[evaluated, "status", "error"] =!= "ok",
    Return[{}, Module]];
  min = Lookup[evaluated["value"], "min", 1];
  max = Lookup[evaluated["value"], "max", 0];
  dimension = Lookup[evaluated["value"], "dimension", 0];
  decoded = DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60];
  If[!IntegerQ[min] || !IntegerQ[max] || min > max ||
      dimension =!= 3 || !ListQ[decoded] ||
      Length[decoded] =!= (max - min + 1) dimension ||
      !AllTrue[decoded, !FailureQ[#] && NumericQ[#] &],
    Return[{}, Module]];
  <|"Min" -> min,
    "Table" -> ArrayReshape[decoded, {max - min + 1, dimension}]|>];

coefficient[value_Association, power_Integer] := Module[
  {min = value["Min"], table = value["Table"]},
  If[1 <= power - min + 1 <= Length[table],
    table[[power - min + 1]], {}]];
coefficient[_, _] := {};

x = Global`x; t = Global`t; eps = Global`eps;
request = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 1|>,
  "TOrder" -> 1|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-acb-scc-one-plus-two-column",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 1, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  system = <|"Matrix" -> {{0, 0, 0}, {eps, 0, 1}, {0, 1, 0}},
    "Variable" -> x|>;
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  stats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  propagated = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 1]]];
  evaluated = If[FailureQ[propagated], propagated,
    DiffExp2`CppBackend`EvaluatePersistentLocal[propagated,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  {cs, prepared, stats, propagated, evaluated}];

{cs, prepared, stats, propagated, evaluated} = result;
value = decodeValue[evaluated];
eps0 = coefficient[value, 0];
eps1 = coefficient[value, 1];

ok = !AnyTrue[result, FailureQ] && AssociationQ[stats] &&
  Lookup[cs["IntegrationSequence"], "Components", None] ===
    {{1}, {2, 3}} &&
  TrueQ[Lookup[stats, "execution_implemented", False]] &&
  Lookup[stats, "execution_scope", None] ===
    "acb-regular-block-dag-column-v2" &&
  Length[eps0] === 3 && Length[eps1] === 3 &&
  Max[Abs[N[eps0 - {1, 0, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[eps1 - {0, 1/2, 0}, 50]]] < 10^-40 &&
  Lookup[propagated, "SeedLocalComponent", None] === 1 &&
  Lookup[propagated, "BasisIndex", None] === 1 &&
  Lookup[propagated["NativeSummary"], "execution_capability", None] ===
    "acb-regular-block-dag-column-v2" &&
  Lookup[propagated["NativeSummary"], "json_coefficients", None] === 0;

If[AssociationQ[propagated] &&
    StringQ[Lookup[propagated, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[propagated]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": retained regular Acb SCC basis column executes natively"];
If[!ok, Print[InputForm[result]]];
Exit[If[ok, 0, 1]];
