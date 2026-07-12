(* One deliberately tiny end-to-end smoke for the explicit multidimensional
   retained SCC column seam.  It covers a scalar seed feeding a two-component
   target and selection of the second local column of that target block. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

decodeValue[evaluated_] := Module[
  {decoded, min, max, dimension, table},
  If[!AssociationQ[evaluated] ||
      Lookup[evaluated, "status", "error"] =!= "ok",
    Return[<||>, Module]];
  min = Lookup[evaluated["value"], "min", 1];
  max = Lookup[evaluated["value"], "max", 0];
  dimension = Lookup[evaluated["value"], "dimension", 0];
  decoded = DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60];
  If[!IntegerQ[min] || !IntegerQ[max] || min > max ||
      !IntegerQ[dimension] || dimension < 1 || !ListQ[decoded] ||
      Length[decoded] =!= (max - min + 1) dimension ||
      !AllTrue[decoded, !FailureQ[#] && NumericQ[#] &],
    Return[<||>, Module]];
  table = ArrayReshape[decoded, {max - min + 1, dimension}];
  <|"Min" -> min, "Max" -> max, "Dimension" -> dimension,
    "Table" -> table|>];

coefficient[value_Association, power_Integer] := Module[
  {min = Lookup[value, "Min", 1], max = Lookup[value, "Max", 0]},
  If[min <= power <= max, value["Table"][[power - min + 1]], {}]];

x = Global`x; t = Global`t; eps = Global`eps;
request = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 1|>,
  "TOrder" -> 1|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-scc-one-plus-two-block-column",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
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
  propagated = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 1]]];
  propagatedEvaluation = If[FailureQ[propagated], propagated,
    DiffExp2`CppBackend`EvaluatePersistentLocal[propagated,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  vectorSeed = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 2, 2]]];
  vectorEvaluation = If[FailureQ[vectorSeed], vectorSeed,
    DiffExp2`CppBackend`EvaluatePersistentLocal[vectorSeed,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  {cs, prepared, propagated, propagatedEvaluation,
    vectorSeed, vectorEvaluation}];

{cs, prepared, propagated, propagatedEvaluation,
  vectorSeed, vectorEvaluation} = result;
propagatedValue = decodeValue[propagatedEvaluation];
vectorValue = decodeValue[vectorEvaluation];
propagated0 = coefficient[propagatedValue, 0];
propagated1 = coefficient[propagatedValue, 1];
vector0 = coefficient[vectorValue, 0];
vector1 = coefficient[vectorValue, 1];

ok = !AnyTrue[{cs, prepared, propagated, propagatedEvaluation,
      vectorSeed, vectorEvaluation}, FailureQ] &&
  AssociationQ[cs] &&
  Lookup[cs["IntegrationSequence"], "Components", None] ===
    {{1}, {2, 3}} &&
  Length[propagated0] === 3 && Length[propagated1] === 3 &&
  Length[vector0] === 3 && Length[vector1] === 3 &&
  Max[Abs[N[propagated0 - {1, 0, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[propagated1 - {0, 1/2, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[vector0 - {0, 1/2, 1}, 50]]] < 10^-40 &&
  Max[Abs[N[vector1, 50]]] < 10^-40 &&
  Lookup[propagated, "Type", None] ===
    "DiffExp2NativeSCCBasisColumn" &&
  Lookup[propagated, "SeedBlock", None] === 1 &&
  Lookup[propagated, "SeedLocalComponent", None] === 1 &&
  Lookup[propagated, "BasisIndex", None] === 1 &&
  Lookup[vectorSeed, "SeedBlock", None] === 2 &&
  Lookup[vectorSeed, "SeedLocalComponent", None] === 2 &&
  Lookup[vectorSeed, "BasisIndex", None] === 3 &&
  Lookup[propagated["NativeSummary"], "execution_capability", None] ===
    "exact-rational-regular-block-dag-column-v2" &&
  Lookup[vectorSeed["NativeSummary"], "execution_capability", None] ===
    "exact-rational-regular-block-dag-column-v2" &&
  Lookup[propagated["NativeSummary"], "json_coefficients", None] === 0 &&
  Lookup[vectorSeed["NativeSummary"], "json_coefficients", None] === 0;

If[AssociationQ[propagated] &&
    StringQ[Lookup[propagated, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[propagated]]];
If[AssociationQ[vectorSeed] && StringQ[Lookup[vectorSeed, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[vectorSeed]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": native retained one-plus-two-dimensional SCC basis columns"];
Exit[If[ok, 0, 1]];
