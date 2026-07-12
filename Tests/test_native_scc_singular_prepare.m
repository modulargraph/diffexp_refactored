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
acbIdentity = If[AssociationQ[acbColumn], Quiet[Check[
    ImportString[Lookup[acbColumn["ColumnProvenance"],
      "exact_column_identity", ""], "RawJSON"], $Failed]], $Failed];
acbDiagnostics = If[AssociationQ[acbColumn],
  Lookup[acbColumn["NativeSummary"], "block_diagnostics", {}], {}];
acbSourceRecords = Select[acbDiagnostics,
  AssociationQ[#] && Lookup[#, "block", None] === 1 &&
    Lookup[#, "role", None] === "particular" &&
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
  AssociationQ[acbIdentity] &&
  exactTagEntryQ[Lookup[acbIdentity, "seed", None]] &&
  AllTrue[Lookup[acbIdentity, "targets", {}], exactTagEntryQ] &&
  Length[rationalValue] > 0 && Length[acbValue] === Length[rationalValue] &&
  Max[Abs[N[acbValue - rationalValue, 50]]] < 10^-40;
If[AssociationQ[acbBasis] && ListQ[Lookup[acbBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    acbBasis["Columns"]]];
If[AssociationQ[acbPrepared] &&
    StringQ[Lookup[acbPrepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[acbPrepared]]];
DiffExp2`Solve`ClearSolveCaches[];

(* A genuine cross-block affine collision is outside this Acb slice.  Each
   diagonal block can truthfully be no-pseudo in isolation, so C++ must still
   revalidate the exact source/target schedule and reject CASE-P rather than
   running approximate compensation or dropping to another integration
   path. *)
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
    "Name" -> "native-acb-scc-case-p-rejected"|>];
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
pseudoOk = !AnyTrue[Take[pseudoResult, 3], FailureQ] &&
  FailureQ[pseudoSolved] &&
  TrueQ[Lookup[pseudoStats, "execution_implemented", False]] &&
  Lookup[pseudoStats, "execution_scope", None] ===
    "acb-regular-singular-jordan-block-dag-column-v1" &&
  AllTrue[pseudoBlocks,
    AssociationQ[#] && TrueQ[Lookup[#, "no_pseudo", False]] &] &&
  Lookup[pseudoSolved[[2]], "ID", None] === "E5" &&
  StringContainsQ[
    Lookup[pseudoSolved[[2]], "Detail", ""], "CASE-P"];

automaticShadowBasis = catchDE2[
  DiffExp2`NativeTransport`Private`nativeReceivingBasis[
    pseudoCs, pseudoRequest, 2]];
automaticShadowOk = AssociationQ[automaticShadowBasis] &&
  Lookup[automaticShadowBasis, "Type", None] ===
    "DiffExp2NativeSCCBasis" &&
  Lookup[automaticShadowBasis, "Dimension", None] === 4 &&
  Lookup[Lookup[automaticShadowBasis, "NativeSummary", <||>],
    "specialization_capability", None] ===
      "exact-rational-shadow-to-acb-local-v1" &&
  Lookup[Lookup[automaticShadowBasis, "Columns", {}],
    "BasisIndex", {}] === Range[4];
If[AssociationQ[automaticShadowBasis] &&
    ListQ[Lookup[automaticShadowBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    automaticShadowBasis["Columns"]]];

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
      Lookup[Lookup[#, "value", <||>], "dimension", None] === 4 &];
Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
  shadowImports];
If[AssociationQ[pseudoPrepared] &&
    StringQ[Lookup[pseudoPrepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[pseudoPrepared]]];
DiffExp2`Solve`ClearSolveCaches[];
ok = TrueQ[ok] && TrueQ[acbOk] && TrueQ[pseudoOk] &&
  TrueQ[automaticShadowOk] && TrueQ[shadowImportOk];

Print[If[ok, "PASS", "FAIL"],
  ": Rational/Acb singular Jordan SCC parity, source propagation, and CASE-P guard"];
If[!ok, Print[InputForm[{
  "RationalFailure" -> Select[result, FailureQ],
  "AcbFailure" -> Select[acbResult, FailureQ],
  "RationalValueLength" -> Length[rationalValue],
  "AcbValueLength" -> Length[acbValue],
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
  "AutomaticShadowBasis" -> If[AssociationQ[automaticShadowBasis],
    KeyTake[automaticShadowBasis,
      {"Type", "Dimension", "NativeSummary"}], automaticShadowBasis],
  "ShadowFailures" -> Select[shadowResult, FailureQ],
  "ShadowIdentity" -> shadowIdentity,
  "ShadowImports" -> (If[AssociationQ[#],
      KeyTake[#, {"status", "session", "local", "scc",
        "specialization_capability", "rational_shadow_identity"}], #] & /@
    shadowImports),
  "BadShadowImport" -> badShadowImport,
  "ShadowClosed" -> shadowClosed,
  "ShadowEvaluatedStatus" -> Lookup[shadowEvaluated, "status", None]}]]];
Exit[If[ok, 0, 1]];
