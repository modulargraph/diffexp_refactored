(* One deliberately tiny end-to-end smoke for the explicit retained SCC
   column seam.  Keep this separate from the broad recurrence suite: it is
   intended to run once after rebuilding the native library. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps;
request = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 1|>,
  "TOrder" -> 1|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-scc-two-scalar-column", "Prescriptions" -> {},
  "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 1, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  system = <|"Matrix" -> {{0, 0}, {eps, 0}}, "Variable" -> x|>;
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  solved = If[FailureQ[prepared], prepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      cs, request, 1]]];
  evaluated = If[FailureQ[solved], solved,
    DiffExp2`CppBackend`EvaluatePersistentLocal[solved,
      <|"exact" -> "1/2"|>, <|"tail_estimate" -> False|>, 60]];
  {cs, prepared, solved, evaluated}];

{cs, prepared, solved, evaluated} = result;
decoded = If[AssociationQ[evaluated] &&
    Lookup[evaluated, "status", "error"] === "ok",
  DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60], {}];
valueMin = If[AssociationQ[evaluated],
  Lookup[evaluated["value"], "min", 1], 1];
valueMax = If[AssociationQ[evaluated],
  Lookup[evaluated["value"], "max", 0], 0];
valueDimension = If[AssociationQ[evaluated],
  Lookup[evaluated["value"], "dimension", 0], 0];
valueTable = If[valueDimension === 2 && IntegerQ[valueMin] &&
    IntegerQ[valueMax] && valueMin <= valueMax &&
    ListQ[decoded] &&
    AllTrue[decoded, !FailureQ[#] && NumericQ[#] &] &&
    Length[decoded] === (valueMax - valueMin + 1) valueDimension,
  ArrayReshape[decoded, {valueMax - valueMin + 1, valueDimension}], {}];
eps0 = If[valueTable === {} || !(valueMin <= 0 <= valueMax), {},
  valueTable[[1 - valueMin]]];
eps1 = If[valueTable === {} || !(valueMin <= 1 <= valueMax), {},
  valueTable[[2 - valueMin]]];
provenance = If[AssociationQ[solved],
  Lookup[solved, "ColumnProvenance", <||>], <||>];
nativeSummary = If[AssociationQ[solved],
  Lookup[solved, "NativeSummary", <||>], <||>];

ok = !AnyTrue[{cs, prepared, solved, evaluated}, FailureQ] &&
  AssociationQ[prepared] && AssociationQ[solved] &&
  AssociationQ[evaluated] &&
  Lookup[evaluated, "status", "error"] === "ok" &&
  Length[eps0] === 2 && Length[eps1] === 2 &&
  Max[Abs[N[eps0 - {1, 0}, 50]]] < 10^-40 &&
  Max[Abs[N[eps1 - {0, 1/2}, 50]]] < 10^-40 &&
  Lookup[solved, "Type", None] === "DiffExp2NativeSCCBasisColumn" &&
  Lookup[solved, "SeedBlock", None] === 1 &&
  Lookup[solved, "BasisIndex", None] === 1 &&
  AssociationQ[provenance] &&
  Lookup[provenance, "schema", None] ===
    "diffexp2-retained-scc-column-reference-v1" &&
  Lookup[provenance, "authority", None] ===
    "retained-native-exact-column-owner" &&
  Sort[Keys[provenance]] === Sort[{"schema", "authority", "scc",
    "seed_block", "basis_index", "identity_diagnostics"}] &&
  Lookup[provenance, "seed_block", None] === 0 &&
  Lookup[provenance, "basis_index", None] === 0 &&
  AssociationQ[Lookup[provenance, "identity_diagnostics", None]] &&
  Sort[Keys[provenance["identity_diagnostics"]]] === Sort[{
    "algorithm", "scc_exact_identity_fingerprint",
    "scc_exact_identity_bytes", "exact_column_identity_fingerprint",
    "exact_column_identity_bytes"}] &&
  Lookup[provenance["identity_diagnostics"], "algorithm", None] ===
    "fnv1a64-v1" &&
  StringMatchQ[Lookup[provenance["identity_diagnostics"],
      "scc_exact_identity_fingerprint", ""],
    RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
  StringMatchQ[Lookup[provenance["identity_diagnostics"],
      "exact_column_identity_fingerprint", ""],
    RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
  IntegerQ[Lookup[provenance["identity_diagnostics"],
    "exact_column_identity_bytes", None]] &&
  Lookup[provenance["identity_diagnostics"],
    "exact_column_identity_bytes", 0] > 0 &&
  FreeQ[Keys[provenance],
    "scc_exact_identity" | "exact_column_identity"] &&
  StringLength[ExportString[
    provenance, "RawJSON", "Compact" -> True]] < 1024 &&
  AssociationQ[nativeSummary] &&
  Lookup[nativeSummary, "json_coefficients", None] === 0 &&
  Intersection[Keys[nativeSummary],
    {"assembled", "coefficients", "u", "validity"}] === {};

If[AssociationQ[solved] && StringQ[Lookup[solved, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[solved]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": native retained two-scalar SCC basis column"];
Exit[If[ok, 0, 1]];
