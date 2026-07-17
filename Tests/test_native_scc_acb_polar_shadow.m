(* All diagonal SCCs may be ordinary while a polar cross-edge makes the
   composite regular-singular.  Exact tags and schedules now let Acb execute
   the normal path and certify its coefficient balls directly.  If a ball is
   ambiguous, nativeReceivingBasis must still catch that narrowly identified
   failure, solve the Rational shadow, and specialize the completed columns
   back into the paired Acb owner without printing a false terminal error. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t;
request = <|"EpsWindow" -> <|"Min" -> -3, "CompleteMax" -> 0|>,
  "TOrder" -> 2|>;
chart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> "native-acb-all-regular-polar-shadow",
  "Prescriptions" -> {}, "UseSCCSkeleton" -> True|>;

result = Block[{DiffExp2`Solve`Private`$cppExactDomain = False},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
    "Variables" -> {}}]];
  DiffExp2`Solve`ClearSolveCaches[];
  (* theta[y2] == y1/x shifts the exact source tag from a=0 to a=-1. *)
  system = <|"Matrix" -> {{0, 0}, {1/x^2, 0}}, "Variable" -> x|>;
  cs = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
  prepared = If[FailureQ[cs], cs,
    catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs, request]]];
  stats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  decision = If[FailureQ[prepared], prepared,
    DiffExp2`Solve`Private`sccNativeCachedRationalShadowDecision[
      cs, request, prepared]];
  capturedPrints = {};
  basis = If[FailureQ[prepared], prepared,
    Block[{Print = Function[Null,
        AppendTo[capturedPrints, HoldComplete[##]], HoldAllComplete]},
      catchDE2[DiffExp2`NativeTransport`Private`nativeReceivingBasis[
        cs, request, 2]]]];
  afterStats = If[FailureQ[prepared], prepared,
    DiffExp2`CppBackend`PersistentSCCStatistics[prepared]];
  {cs, prepared, stats, decision, basis, capturedPrints, afterStats}];

{cs, prepared, stats, decision, basis, capturedPrints, afterStats} = result;
forcedFallbackProbeCount = 0;
forcedFallbackBasis = If[FailureQ[prepared], prepared,
  Block[{
      DiffExp2`NativeTransport`Private`nativeCatchDE2Buffered =
        Function[Null,
          forcedFallbackProbeCount++;
          {Failure["DiffExp2", <|"Detail" ->
              "native Acb CASE-P coefficient is numerically ambiguous; requires the exact Rational shadow"|>],
            {}}, HoldAllComplete]},
    catchDE2[DiffExp2`NativeTransport`Private`nativeReceivingBasis[
      cs, request, 2]]]];
falseTerminalPrints = Select[capturedPrints,
  MatchQ[#, HoldComplete["DiffExp2 error ", _String, ": ", _String]] &];
triggerPrintRecognized =
  DiffExp2`NativeTransport`Private`nativeAcbShadowTriggerPrintRecordQ[
    HoldComplete["DiffExp2 error ", "E5", ": ",
      "Detail=\"native Acb regular-singular SCC execution requires the exact Rational shadow for a multi-tag gauge source\""]];
unrelatedPrintPreserved =
  !DiffExp2`NativeTransport`Private`nativeAcbShadowTriggerPrintRecordQ[
    HoldComplete["DiffExp2 error ", "E6", ": ",
      "Detail=\"unrelated terminal failure\""]];
syntheticPolarCoupling = {<|"entries" -> {<|"multiplier" -> <|
  "proven_zero" -> False, "center_pole_order" -> 1|>|>}|>};
nonIdentityDecision =
  DiffExp2`Solve`Private`sccNativeRationalShadowDecision["acb",
    {<|"identity_gauge" -> True, "identity_v" -> False|>},
    syntheticPolarCoupling, {}];
seedCasePDecision =
  DiffExp2`Solve`Private`sccNativeRationalShadowDecision["acb",
    {<|"identity_gauge" -> False, "identity_v" -> False|>}, {},
    {{<|"SeedRun" -> <|"schedule" ->
      {{<|"case" -> "P"|>}}|>|>}}];
ok = !AnyTrue[Take[result, 3], FailureQ] &&
  Lookup[cs["IntegrationSequence"], "Components", None] === {{1}, {2}} &&
  Lookup[stats, "execution_scope", None] ===
    "acb-regular-singular-scalar-block-dag-column-v1" &&
  AllTrue[Lookup[stats, "block_charts", {}],
    TrueQ[Lookup[#, "regular", False]] &] &&
  AssociationQ[decision] &&
  !TrueQ[Lookup[decision, "RequiresRationalShadow", True]] &&
  Lookup[decision, "Certificate", None] ===
    "exact-tagged-acb-polar-cross-coupling" &&
  TrueQ[Lookup[decision, "RationalShadowFallback", False]] &&
  AssociationQ[basis] &&
  Lookup[basis, "Type", None] === "DiffExp2NativeSCCBasis" &&
  Lookup[basis, "Dimension", None] === 2 &&
  Lookup[Lookup[basis, "Columns", {}], "BasisIndex", {}] === {1, 2} &&
  AllTrue[Lookup[basis, "Columns", {}],
    Lookup[Lookup[#, "NativeSummary", <||>],
      "execution_capability", None] ===
        "acb-regular-singular-scalar-block-dag-column-v1" &] &&
  AssociationQ[afterStats] &&
  Lookup[afterStats, "scc_column_solves", None] === 2 &&
  forcedFallbackProbeCount === 1 &&
  AssociationQ[forcedFallbackBasis] &&
  Lookup[forcedFallbackBasis, "Type", None] ===
    "DiffExp2NativeSCCBasis" &&
  Lookup[forcedFallbackBasis, "Dimension", None] === 2 &&
  Lookup[Lookup[forcedFallbackBasis, "Columns", {}],
    "BasisIndex", {}] === {1, 2} &&
  Lookup[Lookup[forcedFallbackBasis, "NativeSummary", <||>],
    "specialization_capability", None] ===
      "exact-rational-shadow-to-acb-local-v1" &&
  !TrueQ[Lookup[nonIdentityDecision,
    "RequiresRationalShadow", True]] &&
  Lookup[nonIdentityDecision, "Certificate", None] ===
    "runtime-exact-schedule-and-tag-gate-required" &&
  !TrueQ[Lookup[nonIdentityDecision,
    "RationalShadowFallback", False]] &&
  !TrueQ[Lookup[seedCasePDecision,
    "RequiresRationalShadow", True]] &&
  Lookup[seedCasePDecision, "Certificate", None] ===
    "exact-schedule-acb-case-p-compensation" &&
  TrueQ[Lookup[seedCasePDecision,
    "RationalShadowFallback", False]] &&
  falseTerminalPrints === {} && TrueQ[triggerPrintRecognized] &&
  TrueQ[unrelatedPrintPreserved];

If[AssociationQ[basis] && ListQ[Lookup[basis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    basis["Columns"]]];
If[AssociationQ[forcedFallbackBasis] &&
    ListQ[Lookup[forcedFallbackBasis, "Columns", None]],
  Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
    forcedFallbackBasis["Columns"]]];
If[AssociationQ[prepared] && StringQ[Lookup[prepared, "SCC", None]],
  Quiet[DiffExp2`CppBackend`ReleasePersistentSCC[prepared]]];
DiffExp2`Solve`ClearSolveCaches[];

Print[If[ok, "PASS", "FAIL"],
  ": all-regular polar Acb SCC uses direct balls with exact fallback"];
If[!ok, Print[InputForm[{
  "Failures" -> Select[result, FailureQ],
  "Stats" -> If[AssociationQ[stats],
    KeyTake[stats, {"execution_scope", "block_charts",
      "rational_shadow_identity"}], stats],
  "Decision" -> decision,
  "AfterStats" -> If[AssociationQ[afterStats],
    KeyTake[afterStats, {"scc_column_solves", "scc_column_solve_ms"}],
    afterStats],
  "NonIdentityDecision" -> nonIdentityDecision,
  "SeedCasePDecision" -> seedCasePDecision,
  "Basis" -> If[AssociationQ[basis],
    KeyTake[basis, {"Type", "Dimension", "NativeSummary"}], basis],
  "ForcedFallbackProbeCount" -> forcedFallbackProbeCount,
  "ForcedFallbackBasis" -> If[AssociationQ[forcedFallbackBasis],
    KeyTake[forcedFallbackBasis,
      {"Type", "Dimension", "NativeSummary"}], forcedFallbackBasis],
  "CapturedPrints" -> capturedPrints,
  "TriggerPrintRecognized" -> triggerPrintRecognized,
  "UnrelatedPrintPreserved" -> unrelatedPrintPreserved}]]];
Exit[If[ok, 0, 1]];
