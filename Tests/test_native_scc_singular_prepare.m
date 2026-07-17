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

decodeEvaluation[evaluated_] := Module[{decoded},
  If[!AssociationQ[evaluated] ||
      Lookup[evaluated, "status", "error"] =!= "ok",
    Return[{}, Module]];
  decoded = DiffExp2`CppBackend`DecodeScalars[
    evaluated["value", "coefficients"], 60];
  If[ListQ[decoded] &&
      AllTrue[decoded, !FailureQ[#] && NumericQ[#] &], decoded, {}]];

decodeEvaluationFrame[evaluated_] := Module[
  {decoded, min, max, dimension},
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
  <|"Min" -> min, "Max" -> max, "Dimension" -> dimension,
    "Table" -> ArrayReshape[decoded,
      {max - min + 1, dimension}]|>];

exactTagEntryQ[entry_] := Module[{run, metadata, tag, a, b, p},
  If[!AssociationQ[entry], Return[False, Module]];
  run = Lookup[entry, "run", <||>];
  metadata = Lookup[entry, "metadata", <||>];
  tag = Lookup[metadata, "tag", <||>];
  {a, b, p} = Lookup[tag, {"a", "b", "p"}, <||>];
  AssociationQ[run] && AssociationQ[a] && AssociationQ[b] &&
    AssociationQ[p] &&
    Lookup[a, "domain", None] === "rational" &&
    Lookup[a, "canonical", None] === "1/2" &&
    Lookup[b, "domain", None] === "rational" &&
    Lookup[b, "canonical", None] === "1/3" &&
    Lookup[p, "domain", None] === "integer" &&
    Lookup[p, "canonical", None] ===
      ToString[Lookup[run, "p", None]]];

x = Global`x; t = Global`t; eps = Global`eps;
lambda = 1/2 + eps/3;
request = <|"EpsWindow" -> <|"Min" -> -3, "CompleteMax" -> 0|>,
  "TOrder" -> 2|>;
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
rationalValue = decodeEvaluation[evaluated];
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
  Lookup[solved["TWindow"], "CompleteMax", None] === 2 &&
  Lookup[solved["NativeSummary"], "execution_capability", None] ===
    "exact-rational-regular-singular-jordan-block-dag-column-v2" &&
  Lookup[evaluated, "status", "error"] === "ok" &&
  Lookup[evaluated["value"], "dimension", None] === 3;

If[AssociationQ[solved] && StringQ[Lookup[solved, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[solved]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

(* Acb column plans are built from exact task/indicial metadata and exact
   encoded zero/one.  Exercise the ordered batch so the Jordan seed and its
   downstream source edge both execute in ComplexBall arithmetic. *)
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
  acbBasis = If[FailureQ[acbPrepared], acbPrepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
      acbCs, request, 2]]];
  acbColumn = If[FailureQ[acbBasis], acbBasis,
    acbBasis["Columns"][[2]]];
  If[!FailureQ[acbColumn],
    DiffExp2`Solve`DropWolframPreparationCaches[];
    ClearSystemCache[];
    Share[]];
  acbEvaluated = If[FailureQ[acbColumn], acbColumn,
    DiffExp2`CppBackend`EvaluatePersistentLocal[acbColumn,
      <|"exact" -> "1/4"|>, <|"tail_estimate" -> False|>, 60]];
  {acbCs, acbPrepared, acbStats, acbBasis, acbColumn, acbEvaluated}];
{acbCs, acbPrepared, acbStats, acbBasis, acbColumn, acbEvaluated} =
  acbResult;
acbValue = decodeEvaluation[acbEvaluated];
rationalFrame = decodeEvaluationFrame[evaluated];
acbFrame = decodeEvaluationFrame[acbEvaluated];
acbCommonMin = If[AssociationQ[rationalFrame] && AssociationQ[acbFrame],
  Max[Lookup[rationalFrame, "Min", 1], Lookup[acbFrame, "Min", 1]], 1];
acbCommonMax = If[AssociationQ[rationalFrame] && AssociationQ[acbFrame],
  Min[Lookup[rationalFrame, "Max", 0], Lookup[acbFrame, "Max", 0]], 0];
rationalCommonValue = If[acbCommonMin <= acbCommonMax &&
    Lookup[rationalFrame, "Dimension", 0] ===
      Lookup[acbFrame, "Dimension", -1],
  Flatten[rationalFrame["Table"][[
    acbCommonMin - rationalFrame["Min"] + 1 ;;
      acbCommonMax - rationalFrame["Min"] + 1]]], {}];
acbCommonValue = If[rationalCommonValue === {}, {},
  Flatten[acbFrame["Table"][[
    acbCommonMin - acbFrame["Min"] + 1 ;;
      acbCommonMax - acbFrame["Min"] + 1]]]];
acbMaximumDifference = If[Length[rationalCommonValue] > 0 &&
    Length[rationalCommonValue] === Length[acbCommonValue],
  Max[Abs[N[acbCommonValue - rationalCommonValue, 50]]], Infinity];
acbProvenance = If[AssociationQ[acbColumn],
  Lookup[acbColumn, "ColumnProvenance", <||>], <||>];
acbIdentityDiagnostics = If[AssociationQ[acbProvenance],
  Lookup[acbProvenance, "identity_diagnostics", <||>], <||>];
acbDiagnostics = If[AssociationQ[acbColumn],
  Lookup[acbColumn["NativeSummary"], "block_diagnostics", {}], {}];
acbSourceRecords = Select[acbDiagnostics,
  AssociationQ[#] && Lookup[#, "block", None] === 1 &&
    Lookup[#, "role", None] === "particular-tag" &&
    Lookup[#, "predecessors", None] === {0} &&
    Lookup[#, "source_sectors", 0] > 0 &];
acbOk = !AnyTrue[acbResult, FailureQ] &&
  Lookup[acbStats, "status", "error"] === "ok" &&
  TrueQ[Lookup[acbStats, "execution_implemented", False]] &&
  Lookup[acbStats, "execution_scope", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  Lookup[acbBasis, "Type", None] === "DiffExp2NativeSCCBasis" &&
  Lookup[acbBasis, "Dimension", None] === 3 &&
  Lookup[Lookup[acbBasis, "Columns", {}], "BasisIndex", {}] ===
    Range[3] &&
  Lookup[acbBasis["NativeSummary"], "worker_threads", None] === 2 &&
  TrueQ[Lookup[acbBasis["NativeSummary"], "atomic_retention", False]] &&
  Lookup[acbColumn, "SeedLocalComponent", None] === 2 &&
  Lookup[acbColumn["TWindow"], "CompleteMax", None] === 2 &&
  Lookup[acbColumn["NativeSummary"], "execution_capability", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  Length[acbSourceRecords] === 1 &&
  Lookup[First[acbSourceRecords], "result_taylor_max", 0] >= 2 &&
  AssociationQ[acbProvenance] &&
  Lookup[acbProvenance, "schema", None] ===
    "diffexp2-retained-scc-column-reference-v1" &&
  Lookup[acbProvenance, "authority", None] ===
    "retained-native-exact-column-owner" &&
  Lookup[acbProvenance, "scc", None] === acbPrepared["SCC"] &&
  Lookup[acbProvenance, "seed_block", None] ===
    acbColumn["SeedBlock"] - 1 &&
  Lookup[acbProvenance, "basis_index", None] === 1 &&
  AssociationQ[acbIdentityDiagnostics] &&
  Lookup[acbIdentityDiagnostics, "algorithm", None] ===
    "fnv1a64-v1" &&
  StringMatchQ[Lookup[acbIdentityDiagnostics,
      "scc_exact_identity_fingerprint", ""],
    RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
  StringMatchQ[Lookup[acbIdentityDiagnostics,
      "exact_column_identity_fingerprint", ""],
    RegularExpression["^fnv1a64:[0-9a-f]{16}$"]] &&
  IntegerQ[Lookup[acbIdentityDiagnostics,
    "exact_column_identity_bytes", None]] &&
  Lookup[acbIdentityDiagnostics, "exact_column_identity_bytes", 0] > 0 &&
  FreeQ[Keys[acbProvenance],
    "scc_exact_identity" | "exact_column_identity"] &&
  StringLength[ExportString[
    acbProvenance, "RawJSON", "Compact" -> True]] < 1024 &&
  Length[rationalCommonValue] > 0 &&
  Length[acbCommonValue] === Length[rationalCommonValue] &&
  acbMaximumDifference < 10^-40;
If[AssociationQ[acbBasis] && ListQ[Lookup[acbBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    acbBasis["Columns"]]];
If[AssociationQ[acbPrepared] &&
    StringQ[Lookup[acbPrepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[acbPrepared]]];
DiffExp2`Solve`ClearSolveCaches[];

(* A genuine cross-block affine collision is certified from exact tags and
   schedules, while its coefficients stay in ComplexBall arithmetic.  Each
   diagonal block can truthfully be no-pseudo in isolation, so this exercises
   the composite runtime CASE-P gate and the ball-certified compensation. *)
pseudoRequest = <|
  "EpsWindow" -> <|"Min" -> -2, "CompleteMax" -> 0|>,
  "TOrder" -> 0|>;
pseudoResult = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  pseudoSystem = <|"Matrix" -> {
      {1/x, 0, 0, 0},
      {0, (1 + eps)/x, 1/x, 0},
      {1/x, 1, (1 + eps)/x, 0},
      {0, 1/x, 0, (3/2 + 3 eps)/x}}, "Variable" -> x|>;
  pseudoChart = Join[chart, <|
    "Name" -> "native-acb-scc-case-p-certified"|>];
  pseudoCs = catchDE2[DiffExp2`Solve`PrepareChart[
    pseudoSystem, pseudoChart]];
  pseudoPrepared = If[FailureQ[pseudoCs], pseudoCs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[
      pseudoCs, pseudoRequest]]];
  pseudoStats = If[FailureQ[pseudoPrepared], pseudoPrepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[pseudoPrepared]];
  pseudoSolved = If[FailureQ[pseudoPrepared], pseudoPrepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasisColumn[
      pseudoCs, pseudoRequest, 1]]];
  {pseudoCs, pseudoPrepared, pseudoStats, pseudoSolved}];
{pseudoCs, pseudoPrepared, pseudoStats, pseudoSolved} = pseudoResult;
pseudoBlocks = If[AssociationQ[pseudoStats],
  Lookup[pseudoStats, "block_charts", {}], {}];
pseudoDiagnostics = If[AssociationQ[pseudoSolved],
  Lookup[pseudoSolved["NativeSummary"], "block_diagnostics", {}], {}];
pseudoCollisionRecords = Select[pseudoDiagnostics,
  AssociationQ[#] && Lookup[#, "block", None] === 1 &&
    Lookup[#, "source_b", None] === "0" &];
pseudoCollisionOk = Length[pseudoCollisionRecords] === 1 && With[
  {record = First[pseudoCollisionRecords]},
  Lookup[record, "pseudo_hit_count", None] === 1 &&
    Lookup[record, "pseudo_compensation_count", None] === 2 &&
    Lookup[record, "max_pseudo_depth", None] === 2 &&
    TrueQ[Lookup[record, "pseudo_value_certified", False]] &&
    Lookup[record, "uncompensated_pseudo_hit_count", None] === 0];
pseudoOk = !AnyTrue[Take[pseudoResult, 3], FailureQ] &&
  AssociationQ[pseudoSolved] &&
  TrueQ[Lookup[pseudoStats, "execution_implemented", False]] &&
  Lookup[pseudoStats, "execution_scope", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  AllTrue[pseudoBlocks,
    AssociationQ[#] && TrueQ[Lookup[#, "no_pseudo", False]] &] &&
  Lookup[pseudoSolved["NativeSummary"], "execution_capability", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  TrueQ[pseudoCollisionOk];

automaticAcbBasis = catchDE2[
  DiffExp2`NativeTransport`Private`nativeReceivingBasis[
    pseudoCs, pseudoRequest, 2]];
automaticAcbOk = AssociationQ[automaticAcbBasis] &&
  Lookup[automaticAcbBasis, "Type", None] ===
    "DiffExp2NativeSCCBasis" &&
  Lookup[automaticAcbBasis, "Dimension", None] === 4 &&
  Lookup[Lookup[automaticAcbBasis, "Columns", {}],
    "BasisIndex", {}] === Range[4] &&
  AllTrue[Lookup[automaticAcbBasis, "Columns", {}],
    Lookup[Lookup[#, "NativeSummary", <||>],
      "execution_capability", None] ===
        "acb-regular-singular-jordan-block-dag-column-v1" &];
If[AssociationQ[automaticAcbBasis] &&
    ListQ[Lookup[automaticAcbBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    automaticAcbBasis["Columns"]]];

(* The exact Rational CASE-P solver is a proof-producing shadow for this one
   SCC only.  Import its completed columns into the already prepared Acb SCC,
   close the Rational session, and verify that the Acb locals remain live. *)
shadowResult = Block[{DiffExp2`Solve`Private`$cppExactDomain = True},
  shadowCs = catchDE2[DiffExp2`Solve`PrepareChart[
    pseudoSystem, pseudoChart]];
  shadowPrepared = If[FailureQ[shadowCs], shadowCs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[
      shadowCs, pseudoRequest]]];
  shadowStats = If[FailureQ[shadowPrepared], shadowPrepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[shadowPrepared]];
  shadowBasis = If[FailureQ[shadowPrepared], shadowPrepared,
    catchDE2[DiffExp2`Solve`SolveNativeSCCBasis[
      shadowCs, pseudoRequest, 2]]];
  {shadowCs, shadowPrepared, shadowStats, shadowBasis}];
shadowIdentity = If[AssociationQ[pseudoStats],
  Lookup[pseudoStats, "rational_shadow_identity", None], None];
shadowImports = If[AssociationQ[shadowBasis] &&
    ListQ[Lookup[shadowBasis, "Columns", None]] &&
    StringQ[shadowIdentity],
  MapIndexed[DiffExp2`CppBackend`SpecializePersistentRationalSCCColumn[
      #1, pseudoPrepared, shadowIdentity,
      "native-scc-shadow-import:" <> ToString[First[#2]]] &,
    shadowBasis["Columns"]], {}];
badShadowImport = If[AssociationQ[shadowBasis] &&
    ListQ[Lookup[shadowBasis, "Columns", None]],
  DiffExp2`CppBackend`SpecializePersistentRationalSCCColumn[
    First[shadowBasis["Columns"]], pseudoPrepared, "mismatched-shadow",
    "native-scc-shadow-import:bad"], None];
shadowClosed = If[AssociationQ[shadowPrepared],
  DiffExp2`CppBackend`ClosePersistentSession[shadowPrepared], None];
shadowEvaluated = Map[
  DiffExp2`CppBackend`EvaluatePersistentLocal[#,
      <|"exact" -> "1/4"|>, <|"tail_estimate" -> False|>, 60] &,
  shadowImports];
pseudoEvaluated = If[AssociationQ[pseudoSolved],
  DiffExp2`CppBackend`EvaluatePersistentLocal[pseudoSolved,
    <|"exact" -> "1/4"|>, <|"tail_estimate" -> False|>, 60],
  pseudoSolved];
pseudoFrame = decodeEvaluationFrame[pseudoEvaluated];
shadowFrames = decodeEvaluationFrame /@ shadowEvaluated;
shadowFrame = If[shadowFrames === {}, <||>, First[shadowFrames]];
shadowSourceSummaries = If[AssociationQ[shadowBasis],
  Lookup[Lookup[shadowBasis, "Columns", {}], "NativeSummary", {}], {}];
shadowSourceMaxima =
  Lookup[shadowSourceSummaries, "epsilon_max", Missing["epsilon_max"]];
shadowSourceTopValid =
  Lookup[shadowSourceSummaries, "top_valid", Missing["top_valid"]];
shadowFrameMaxima =
  Lookup[shadowFrames, "Max", Missing["frame_max"]];
shadowReservoirOk = Length[shadowFrames] === 4 &&
  shadowFrameMaxima === shadowSourceMaxima &&
  And @@ Thread[shadowSourceTopValid >= shadowSourceMaxima] &&
  AllTrue[shadowFrameMaxima,
    # > Lookup[pseudoRequest["EpsWindow"], "CompleteMax", Infinity] &];
commonMin = If[AssociationQ[pseudoFrame] && AssociationQ[shadowFrame],
  Max[Lookup[pseudoFrame, "Min", 1], Lookup[shadowFrame, "Min", 1]], 1];
commonMax = If[AssociationQ[pseudoFrame] && AssociationQ[shadowFrame],
  Min[Lookup[pseudoFrame, "Max", 0], Lookup[shadowFrame, "Max", 0]], 0];
pseudoCommonValue = If[commonMin <= commonMax &&
    Lookup[pseudoFrame, "Dimension", 0] ===
      Lookup[shadowFrame, "Dimension", -1],
  Flatten[pseudoFrame["Table"][[
    commonMin - pseudoFrame["Min"] + 1 ;;
      commonMax - pseudoFrame["Min"] + 1]]], {}];
shadowCommonValue = If[pseudoCommonValue === {}, {},
  Flatten[shadowFrame["Table"][[
    commonMin - shadowFrame["Min"] + 1 ;;
      commonMax - shadowFrame["Min"] + 1]]]];
pseudoShadowMaximumDifference = If[Length[pseudoCommonValue] > 0 &&
    Length[pseudoCommonValue] === Length[shadowCommonValue],
  Max[Abs[N[pseudoCommonValue - shadowCommonValue, 50]]], Infinity];
pseudoShadowParityOk = Length[pseudoCommonValue] > 0 &&
  Length[pseudoCommonValue] === Length[shadowCommonValue] &&
  pseudoShadowMaximumDifference < 10^-40;
shadowImportOk = Length[shadowImports] === 4 &&
  AllTrue[shadowImports, AssociationQ[#] &&
      Lookup[#, "status", "error"] === "ok" &&
      Lookup[#, "specialization_capability", None] ===
        "exact-rational-shadow-to-acb-local-v1" &&
      Lookup[#, "rational_shadow_identity", None] === shadowIdentity &] &&
  AssociationQ[badShadowImport] &&
  Lookup[badShadowImport, "status", "ok"] === "error" &&
  AssociationQ[shadowClosed] &&
  AllTrue[shadowEvaluated, AssociationQ[#] &&
      Lookup[#, "status", "error"] === "ok" &&
      Lookup[Lookup[#, "value", <||>], "dimension", None] === 4 &] &&
  TrueQ[shadowReservoirOk] &&
  TrueQ[pseudoShadowParityOk];
If[AssociationQ[pseudoSolved] &&
    StringQ[Lookup[pseudoSolved, "Local", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[pseudoSolved]]];
Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
  shadowImports];
If[AssociationQ[pseudoPrepared] &&
    StringQ[Lookup[pseudoPrepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[pseudoPrepared]]];
DiffExp2`Solve`ClearSolveCaches[];
ok = TrueQ[ok] && TrueQ[acbOk] && TrueQ[pseudoOk] &&
  TrueQ[automaticAcbOk] && TrueQ[shadowImportOk];

Print[If[ok, "PASS", "FAIL"],
  ": Rational/Acb singular Jordan SCC parity, source propagation, and ball-certified CASE-P"];
If[!ok, Print[InputForm[{
  "RationalFailure" -> Select[result, FailureQ],
  "AcbFailure" -> Select[acbResult, FailureQ],
  "RationalValueLength" -> Length[rationalValue],
  "AcbValueLength" -> Length[acbValue],
  "RationalFrame" -> KeyDrop[rationalFrame, "Table"],
  "AcbFrame" -> KeyDrop[acbFrame, "Table"],
  "AcbCommonFrame" -> {acbCommonMin, acbCommonMax},
  "AcbCommonValueLengths" ->
    {Length[rationalCommonValue], Length[acbCommonValue]},
  "AcbMaximumDifference" -> acbMaximumDifference,
  "AcbStats" -> If[AssociationQ[acbStats],
    KeyTake[acbStats, {"execution_implemented", "execution_scope",
      "capability_evidence"}], acbStats],
  "AcbBasis" -> If[AssociationQ[acbBasis],
    KeyTake[acbBasis, {"Type", "Dimension", "NativeSummary"}],
    acbBasis],
  "PseudoStats" -> If[AssociationQ[pseudoStats],
    KeyTake[pseudoStats, {"execution_implemented", "execution_scope",
      "block_charts"}], pseudoStats],
  "PseudoSolved" -> pseudoSolved,
  "PseudoCollisionRecords" -> pseudoCollisionRecords,
  "AutomaticAcbBasis" -> If[AssociationQ[automaticAcbBasis],
    KeyTake[automaticAcbBasis,
      {"Type", "Dimension", "NativeSummary"}], automaticAcbBasis],
  "ShadowFailures" -> Select[shadowResult, FailureQ],
  "ShadowIdentity" -> shadowIdentity,
  "ShadowImports" -> (If[AssociationQ[#],
      KeyTake[#, {"status", "session", "local", "scc",
        "specialization_capability", "rational_shadow_identity"}], #] & /@
    shadowImports),
  "BadShadowImport" -> badShadowImport,
  "ShadowClosed" -> shadowClosed,
  "ShadowEvaluatedStatus" -> Lookup[shadowEvaluated, "status", None],
  "PseudoEvaluatedStatus" -> Lookup[pseudoEvaluated, "status", None],
  "PseudoFrame" -> KeyDrop[pseudoFrame, "Table"],
  "ShadowFrame" -> KeyDrop[shadowFrame, "Table"],
  "ShadowReservoirOk" -> shadowReservoirOk,
  "CommonFrame" -> {commonMin, commonMax},
  "PseudoShadowCommonValueLengths" ->
    {Length[pseudoCommonValue], Length[shadowCommonValue]},
  "PseudoShadowMaximumDifference" ->
    pseudoShadowMaximumDifference }]]];
Exit[If[ok, 0, 1]];
