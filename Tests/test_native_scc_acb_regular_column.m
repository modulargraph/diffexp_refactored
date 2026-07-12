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

decodeValue[evaluated_, expectedDimension_:3] := Module[
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
      dimension =!= expectedDimension || !ListQ[decoded] ||
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
  basis = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
      cs, request, 2]]];
  valueLocal = If[FailureQ[basis], basis,
    catchDE2[DiffExp2`Solve`SolveNativeValueRegular[cs, request,
      {DiffExp2`EpsSeries`ESNew[0, {1, 0}],
       DiffExp2`EpsSeries`ESNew[0, {0, 0}],
       DiffExp2`EpsSeries`ESNew[0, {0, 0}]}]]];
  valueEvaluation = If[FailureQ[valueLocal], valueLocal,
    DiffExp2`CppBackend`EvaluatePersistentLocal[valueLocal,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  propagated = If[FailureQ[basis], basis, First[basis["Columns"]]];
  evaluated = If[FailureQ[propagated], propagated,
    DiffExp2`CppBackend`EvaluatePersistentLocal[propagated,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  monoSystem = <|"Matrix" -> {{0, 1}, {1, 0}}, "Variable" -> x|>;
  monoChart = Join[chart, <|
    "Name" -> "native-acb-monolithic-regular-basis"|>];
  monoCs = catchDE2[DiffExp2`Solve`PrepareChart[monoSystem, monoChart]];
  monoBasis = If[FailureQ[monoCs], monoCs,
    catchDE2[DiffExp2`Solve`SolveNativeRegularBasis[
      monoCs, request, 2]]];
  monoEvaluations = If[FailureQ[monoBasis], {monoBasis},
    DiffExp2`CppBackend`EvaluatePersistentLocal[#,
      <|"exact" -> "0"|>, <|"tail_estimate" -> False|>, 60] & /@
      monoBasis["Columns"]];
  {cs, prepared, stats, basis, valueLocal, valueEvaluation,
    propagated, evaluated, monoCs, monoBasis, monoEvaluations}];

{cs, prepared, stats, basis, valueLocal, valueEvaluation,
  propagated, evaluated, monoCs, monoBasis, monoEvaluations} = result;
value = decodeValue[evaluated];
transportedValue = decodeValue[valueEvaluation];
eps0 = coefficient[value, 0];
eps1 = coefficient[value, 1];
transported0 = coefficient[transportedValue, 0];
transported1 = coefficient[transportedValue, 1];
monoValues = decodeValue[#, 2] & /@ monoEvaluations;
mono0 = coefficient[#, 0] & /@ monoValues;

ok = !AnyTrue[result, FailureQ] && AssociationQ[stats] &&
  Lookup[cs["IntegrationSequence"], "Components", None] ===
    {{1}, {2, 3}} &&
  TrueQ[Lookup[stats, "execution_implemented", False]] &&
  Lookup[stats, "execution_scope", None] ===
    "acb-regular-block-dag-column-v2" &&
  Lookup[basis, "Type", None] === "DiffExp2NativeSCCBasis" &&
  Lookup[basis, "Dimension", None] === 3 &&
  Lookup[Lookup[basis, "Columns", {}], "BasisIndex", {}] ===
    Range[3] &&
  Lookup[basis["NativeSummary"], "columns", None] === 3 &&
  Lookup[basis["NativeSummary"], "worker_threads", None] === 2 &&
  TrueQ[Lookup[basis["NativeSummary"], "atomic_retention", False]] &&
  Lookup[valueLocal, "Type", None] === "DiffExp2NativeValueRegular" &&
  Lookup[valueLocal, "Session", None] === Lookup[basis, "Session", None] &&
  Length[transported0] === 3 && Length[transported1] === 3 &&
  Max[Abs[N[transported0 - {1, 0, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[transported1 - {0, 1/2, 0}, 50]]] < 10^-40 &&
  Lookup[monoCs["IntegrationSequence"], "Components", None] ===
    {{1, 2}} &&
  Lookup[monoBasis, "Type", None] === "DiffExp2NativeRegularBasis" &&
  Lookup[monoBasis, "Dimension", None] === 2 &&
  Lookup[monoBasis["NativeSummary"], "execution_capability", None] ===
    "retained-regular-monolithic-unit-basis-v1" &&
  Length[mono0] === 2 &&
  Max[Abs[N[mono0[[1]] - {1, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[mono0[[2]] - {0, 1}, 50]]] < 10^-40 &&
  Length[eps0] === 3 && Length[eps1] === 3 &&
  Max[Abs[N[eps0 - {1, 0, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[eps1 - {0, 1/2, 0}, 50]]] < 10^-40 &&
  Lookup[propagated, "SeedLocalComponent", None] === 1 &&
  Lookup[propagated, "BasisIndex", None] === 1 &&
  Lookup[propagated["NativeSummary"], "execution_capability", None] ===
    "acb-regular-block-dag-column-v2" &&
  Lookup[propagated["NativeSummary"], "json_coefficients", None] === 0;

If[AssociationQ[basis] && ListQ[Lookup[basis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    basis["Columns"]]];
If[AssociationQ[valueLocal] && StringQ[Lookup[valueLocal, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[valueLocal]]];
If[AssociationQ[monoBasis] &&
    ListQ[Lookup[monoBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    monoBasis["Columns"]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": retained regular Acb SCC basis column executes natively"];
If[!ok, Print[InputForm[{
  "FailurePositions" -> Select[Range[Length[result]],
    FailureQ[result[[#]]] &],
  "BasisType" -> If[AssociationQ[basis], Lookup[basis, "Type", None],
    Head[basis]],
  "ValueType" -> If[AssociationQ[valueLocal],
    Lookup[valueLocal, "Type", None], Head[valueLocal]],
  "ValueSession" -> If[AssociationQ[valueLocal],
    Lookup[valueLocal, "Session", None], None],
  "BasisSession" -> If[AssociationQ[basis],
    Lookup[basis, "Session", None], None],
  "Transported" -> {transported0, transported1},
  "MonolithicType" -> If[AssociationQ[monoBasis],
    Lookup[monoBasis, "Type", None], Head[monoBasis]],
  "MonolithicValues" -> mono0,
  "BasisColumn" -> {eps0, eps1}}]]];
Exit[If[ok, 0, 1]];
