(* Focused level-facing retained transport-observable batch smoke. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];
nativeOKQ[value_] := AssociationQ[value] &&
  Lookup[value, "status", "error"] === "ok";

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 1, "DivisionOrder" -> 2,
  "Variables" -> {}, "Verbosity" -> 0}];

x = Global`x;

(* A singular basis can store all columns one epsilon row above the physical
   source obtained after inverse Laurent matching weights.  Row preparation
   must use the atlas work floor, while retaining the basis's public upper
   edge rather than any private matching reservoir. *)
rowFixtureChart = <|"SystemSize" -> 1, "ChartVar" -> Global`t,
  "Center" -> 0, "ChartMap" -> <|"Center" -> 0, "Scale" -> 1|>,
  "Radius" -> 2, "Prescriptions" -> {}|>;
rowFixtureAtlas = <|"Dimension" -> 1,
  "Request" -> <|"EpsWindow" -> <|"Min" -> 0,
      "CompleteMax" -> 2|>, "TOrder" -> 4|>,
  "Anchor" -> <|"EpsWindow" -> <|"Min" -> 0,
      "CompleteMax" -> 2|>,
    "TWindow" -> <|"CompleteMax" -> 4|>|>|>;
rowFixtureBasis = <|"Dimension" -> 1, "Columns" -> {
    <|"EpsWindow" -> <|"Min" -> 1, "CompleteMax" -> 2|>,
      "TWindow" -> <|"CompleteMax" -> 4|>|>}|>;
rowFixtureRecipes =
  DiffExp2`NativeTransport`Private`nativeArmRowRecipes[
    rowFixtureAtlas,
    <|"ChartSystems" -> {rowFixtureChart, rowFixtureChart},
      "Bases" -> {None, rowFixtureBasis}|>];
rowFixturePrepared =
  DiffExp2`NativeTransport`Private`nativePrepareArmRecipeRow[
    Last[rowFixtureRecipes], {1/(1 - x/2)}, x, "rational"];
assert["singular row recipe restores the physical atlas epsilon floor",
  rowFixtureRecipes[[2, "Shape", "EpsWindow"]] ===
      <|"Min" -> 0, "CompleteMax" -> 2|> &&
    Length[rowFixturePrepared["entries"][[1, "multiplier",
      "analytic_coefficients"]]] === 3];

(* The banana4 final tile is a concrete counterexample to using the basis
   envelope alone: its live source is eps[-5,4], the integrand begins at
   eps^-4, and the regulated primitive consumes through eps^1.  The lazy
   producer must therefore emit ten multiplier coefficients even though its
   pre-march recipe had only nine rows. *)
bananaConsumerRecipe =
  DiffExp2`NativeTransport`Private`nativeConsumerRowRecipe[
    <|"Chart" -> rowFixtureChart,
      "Shape" -> <|"EpsWindow" -> <|"Min" -> -1,
          "CompleteMax" -> 7|>,
        "TWindow" -> <|"CompleteMax" -> 4|>,
        "Dimension" -> 1|>|>,
    <|"Min" -> -5, "CompleteMax" -> 4|>,
    <|"Identity" -> "banana4-consumer-width",
      "MinimumEpsilonShift" -> -4,
      "Epsilon" -> <|"Min" -> -1, "Max" -> 0,
        "RequiredCompleteMax" -> 0|>|>, 1];
bananaConsumerPrepared =
  DiffExp2`NativeTransport`Private`nativePrepareArmRecipeRow[
    bananaConsumerRecipe,
    {DiffExp2`Config`CanonicalEps[]^-4/(1 - x/2)}, x, "rational"];
assert["regulated consumer widens the exact banana4 multiplier prefix by one",
  bananaConsumerRecipe["Shape", "EpsWindow"] ===
      <|"Min" -> -1, "CompleteMax" -> 8|> &&
    Length[bananaConsumerPrepared["entries"][[1, "multiplier",
      "analytic_coefficients"]]] === 10];

system = DiffExp2`LoadSystem[<|"Matrix" -> {{0}}, "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {0, -1/4}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {0, 1/4}];
one = DiffExp2`EpsSeries`ESNew[0, {1, 0}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrands" -> {{{1 + Global`eps},
      {1/Global`eps + 1}}, x},
    "TargetCompleteMax" -> 0]];

epsilon = <|"Min" -> -1, "Max" -> 0,
  "RequiredCompleteMax" -> -1|>;
observable[operation_, identity_, coefficients_:{1 + Global`eps}] := <|
  "Operation" -> operation, "Identity" -> identity,
  "CheckpointIdentity" -> identity <> ":checkpoint",
  "CoefficientVector" -> coefficients, "Epsilon" -> epsilon|>;
observables = {
  Append[observable["integrate", "integral"], "TailPolicy" -> "require"],
  observable["limitLower", "lower-limit"],
  observable["limitUpper", "upper-limit"],
  Append[observable["integrate", "polar-integral",
      {1/Global`eps + 1}],
    "TailPolicy" -> "stored"],
  Append[observable["integrate", "second-integral"],
    "TailPolicy" -> "stored"]};
sessionStats[] := If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];
compactCounters = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`PersistentSessionCounters[atlas]];
assert["native_session_counters_are_fixed_size_and_summary_free",
  AssociationQ[compactCounters] &&
    Lookup[compactCounters, "status", "error"] === "ok" &&
    Lookup[compactCounters, "scope", None] ===
      "fixed-session-counters" &&
    Lookup[compactCounters, "session", None] === atlas["Session"] &&
    ByteCount[compactCounters] < 10000 &&
    FreeQ[Keys[compactCounters],
      "chart_stats" | "local_stats" | "match_stats" |
      "transport_state_stats" | "retained_derivation"]];

(* A streamed hop certifies every retained row which a later observable can
   consume.  Only the gap above that required edge is private match work. *)
capturedHopEpsilon = None;
streamFixtureAnchor = <|"session" -> "stream-fixture-session",
  "local" -> "l:anchor", "checkpoint_identity" -> "anchor-checkpoint"|>;
streamFixture = catchDE2[Block[{
    DiffExp2`NativeTransport`Private`nativeReceivingBasis =
      Function[{receiver, request, threads, consume, equationOwner},
        <|"Columns" -> {<|"session" -> "stream-fixture-session",
            "local" -> "l:basis",
            "checkpoint_identity" -> "basis-checkpoint"|>}|>],
    DiffExp2`CppBackend`ConsumePersistentTransportHop =
      Function[{plan, arm, index, basis, incoming, epsilon, root,
          refinement},
        capturedHopEpsilon = epsilon;
        <|"status" -> "ok", "next_local" ->
            <|"local" -> "l:next",
              "checkpoint_identity" -> "next-checkpoint"|>,
          "consumed_basis_handles" -> {"l:basis"}|>]},
  DiffExp2`NativeTransport`Private`nativeStreamTransportArm[
    <|"Session" -> "stream-fixture-session",
      "Anchor" -> streamFixtureAnchor, "Plan" -> <||>, "Threads" -> 1,
      "Request" -> <|"EpsWindow" -> <|"Min" -> 0,
          "CompleteMax" -> 8|>, "RequiredCompleteMax" -> 2|>|>,
    <|"ChartSystems" -> {<||>, <|"Center" -> 1,
          "IndicialData" -> <|"Regular" -> True|>|>},
      "ValueSolvers" -> {None}, "OwnerRecords" -> {}|>,
    "lower", <|"min" -> 0, "max" -> 8,
      "required_complete_max" -> 7,
      "match_required_complete_max" -> 7|>,
    "stream-fixture-root", <|"relative_tolerance" -> "1e-8",
      "max_steps" -> 2|>]]];
assert["streamed hop certifies the observable-consumed reservoir below the private work top",
  AssociationQ[streamFixture] &&
    Length[streamFixture["tile_sources"]] === 2 &&
    capturedHopEpsilon === <|"min" -> 0, "max" -> 8,
      "required_complete_max" -> 7|>];

beforeInvalid = sessionStats[];
duplicate = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, {observable["limitLower", "duplicate"],
      observable["limitUpper", "duplicate"]}, x]]];
afterInvalid = sessionStats[];
assert["native_observable_batch_rejects_duplicate_identity_before_march",
  FailureQ[duplicate] &&
    Lookup[beforeInvalid, "transport_states", -1] ===
      Lookup[afterInvalid, "transport_states", -2] &&
    Lookup[beforeInvalid, "transport_arm_marches", -1] ===
      Lookup[afterInvalid, "transport_arm_marches", -2]];

(* A definitions-only failure fixture proves that completed chunks and even
   a malformed current response are both reclaimed before the two retained
   states.  No real native owner is allocated by this block. *)
mockAtlas = <|
  "Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
  "Dimension" -> 1, "PlanCheckpointIdentity" -> "mock-plan-checkpoint",
  "Domain" -> "acb", "Session" -> "mock-session",
  "Plan" -> <|"plan" -> "mock-plan"|>,
  "Anchor" -> <|"checkpoint_identity" -> "mock-anchor-checkpoint"|>,
  "PreparedIntegrandEpsilonShift" -> 0,
  "Request" -> <|"EpsWindow" -> <|"Min" -> 0,
      "CompleteMax" -> 1|>|>,
  "Lower" -> <|"Bases" -> {None}, "ChartSystems" -> {None}|>,
  "Upper" -> <|"Bases" -> {None}, "ChartSystems" -> {None}|>|>;
mockIntegrates = Table[
  <|"Operation" -> "integrate", "Identity" -> "mock-" <> ToString[i],
    "CheckpointIdentity" -> "mock-checkpoint-" <> ToString[i],
    "CoefficientVector" -> {1},
    "Epsilon" -> <|"Min" -> 0, "Max" -> 0,
      "RequiredCompleteMax" -> 0|>, "TailPolicy" -> "stored"|>,
  {i, 3}];
mockPairCalls = 0; mockPairRoots = {}; mockPairIdentities = {};
mockReleasedLines = {}; mockReleasedStates = {}; mockMarchCalls = 0;
mockFailure = catchDE2[Block[{
    DiffExp2`NativeTransport`Private`nativeArmRowRecipes =
      Function[{ignoredAtlas, ignoredArm}, {<|"MockRecipe" -> True|>}],
    DiffExp2`NativeTransport`Private`nativePrepareBatchObservable =
      Function[{raw, ignoredVar, ignoredDimension},
        Append[raw, "MinimumEpsilonShift" -> 0]],
    DiffExp2`Solve`DropWolframPreparationCaches = Function[Null, Null],
    DiffExp2`CppBackend`RunPersistentTransportArms =
      Function[{plan, anchor, arms, epsWindow, root, refinement},
        mockMarchCalls++;
        <|"status" -> "ok", "session" -> "mock-session",
          "native_retained" -> True, "json_coefficients" -> 0,
          "states" -> <|
            "lower" -> <|"session" -> "mock-session", "arm" -> "lower",
              "transport_state" -> "transport:mock-lower",
              "checkpoint_identity" -> "mock-state-lower-checkpoint",
              "provenance_identity" -> "mock-state-lower-provenance",
              "tiles" -> 1,
              "tile_source_epsilon" -> {
                <|"min" -> 0, "max" -> 1|>}|>,
            "upper" -> <|"session" -> "mock-session", "arm" -> "upper",
              "transport_state" -> "transport:mock-upper",
              "checkpoint_identity" -> "mock-state-upper-checkpoint",
              "provenance_identity" -> "mock-state-upper-provenance",
              "tiles" -> 1,
              "tile_source_epsilon" -> {
                <|"min" -> 0, "max" -> 1|>}|>|>|>],
    DiffExp2`NativeTransport`Private`nativeContractStoredPairObservableStreamed =
      Function[{lowerState, upperState, request, lowerRecipes,
          upperRecipes, ignoredVar, ignoredDomain, root,
          lowerSourceEpsilon, upperSourceEpsilon}, Module[
        {call, line},
        mockPairCalls++; call = mockPairCalls;
        AppendTo[mockPairRoots, root];
        AppendTo[mockPairIdentities, {request["Identity"]}];
        line = <|"session" -> "mock-session",
          "line" -> "line:mock-" <> ToString[call],
          "checkpoint_identity" -> request["CheckpointIdentity"],
          "provenance_identity" -> "mock-line-provenance-" <>
            ToString[call], "request_index" -> 0,
          "observable_identity" -> request["Identity"]|>;
        <|"status" -> If[call === 2, "error", "ok"],
          "lines" -> {line}|>]],
    DiffExp2`CppBackend`ReleasePersistentLineIntegral =
      Function[handle, AppendTo[mockReleasedLines, handle["line"]];
        <|"status" -> "ok"|>],
    DiffExp2`CppBackend`ReleasePersistentTransportArm =
      Function[handle,
        AppendTo[mockReleasedStates, handle["transport_state"]];
        <|"status" -> "ok"|>]},
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    mockAtlas, mockIntegrates, x]]];
assert["native_observable_chunk_failure_releases_prior_and_current_lines",
  FailureQ[mockFailure] && mockMarchCalls === 1 && mockPairCalls === 2 &&
    mockPairIdentities === {{"mock-1"}, {"mock-2"}} &&
    DuplicateFreeQ[mockPairRoots] &&
    StringEndsQ[mockPairRoots[[1]], ":integrals:chunk:1"] &&
    StringEndsQ[mockPairRoots[[2]], ":integrals:chunk:2"] &&
    Sort[mockReleasedLines] === {"line:mock-1", "line:mock-2"} &&
    Sort[mockReleasedStates] ===
      {"transport:mock-lower", "transport:mock-upper"}];

beforeRun = sessionStats[];
run = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, observables, x, "MaxRefinementSteps" -> 1]]];
afterRun = sessionStats[];
results = If[AssociationQ[run], Lookup[run, "Results", {}], {}];
assert["native_observable_batch_marches_once_and_preserves_request_order",
  AssociationQ[run] &&
    Lookup[run, "Type", None] ===
      "DiffExp2NativeTransportObservableBatch" &&
    Lookup[run, "NativeMarches", 0] === 2 &&
    Lookup[results, "RequestIndex"] === {0, 1, 2, 3, 4} &&
    Lookup[results, "Identity"] ===
      {"integral", "lower-limit", "upper-limit", "polar-integral",
       "second-integral"} &&
    run["NativeSummary", "PairCalls"] === 3 &&
    run["NativeSummary", "Pair", "ChunkCount"] === 3 &&
    Lookup[afterRun, "transport_pair_contractions", -1] ===
      Lookup[beforeRun, "transport_pair_contractions", 0] + 3 &&
    Lookup[afterRun, "transport_arm_marches", -1] ===
      Lookup[beforeRun, "transport_arm_marches", 0] + 2 &&
    atlas["TargetCompleteMax"] === 0 &&
    atlas["Request", "EpsWindow", "CompleteMax"] === 1 &&
    (* The polar integral reserves one source order for its possible
       eps^-1 primitive before the retained arm state is contracted. *)
    run["States", "lower", "epsilon", "required_complete_max"] === 1 &&
    Sort[Keys[Lookup[run, "States", <||>]]] === {"lower", "upper"}];

checkpointPath = FileNameJoin[{$TemporaryDirectory,
  "de2-native-observable-batch-" <> ToString[$ProcessID] <> ".checkpoint"}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
checkpointManifest = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`SaveNativeTransportObservableBatchCheckpoint[
    run, checkpointPath, "native-observable-batch-roundtrip"]], run];
assert["native_observable_batch_checkpoint_is_atomic_and_manifest_bound",
  AssociationQ[checkpointManifest] && FileExistsQ[checkpointPath] &&
    Lookup[checkpointManifest, "Schema", None] ===
      "DiffExp2.NativeTransportObservableCheckpoint/v2" &&
    Lookup[checkpointManifest, "TransportArmMarches", -1] === 2 &&
    Lookup[Lookup[checkpointManifest, "Results", {}], "Identity"] ===
      Lookup[results, "Identity"] &&
    AllTrue[Values[Lookup[checkpointManifest, "StateHandles", <||>]],
      AssociationQ[#] && Sort[Keys[#]] ===
        Sort[{"Handle", "CheckpointIdentity", "ProvenanceSHA256"}] &&
        StringMatchQ[# ["ProvenanceSHA256"],
          RegularExpression["[0-9a-f]{64}"]] &] &&
    DuplicateFreeQ[Lookup[checkpointManifest["Results"],
      "ProvenanceSHA256"]] &&
    AllTrue[Lookup[checkpointManifest["Results"], "ProvenanceSHA256"],
      StringMatchQ[#, RegularExpression["[0-9a-f]{64}"]] &]];

(* Preserve one exact v1 fixture before releasing the live session.  V2 is
   what new sidecars persist; the legacy manifest remains loadable against
   the same schema-2 native checkpoint. *)
legacySchema = "DiffExp2.NativeTransportObservableCheckpoint/v1";
legacyStateHandles = AssociationMap[Function[side, <|
    "Handle" -> run["States", side, "transport_state"],
    "CheckpointIdentity" -> run["States", side, "checkpoint_identity"],
    "ProvenanceIdentity" -> run["States", side, "provenance_identity"]|>],
  {"lower", "upper"}];
legacyCore = <|"Schema" -> legacySchema,
  "Path" -> checkpointManifest["Path"],
  "CheckpointIdentity" -> checkpointManifest["CheckpointIdentity"],
  "TransportArmMarches" -> checkpointManifest["TransportArmMarches"],
  "StateHandles" -> legacyStateHandles,
  "Results" ->
    (DiffExp2`NativeTransport`Private`nativeObservableCheckpointResult[
        #, run["Atlas", "Session"], legacySchema] & /@ results)|>;
legacyManifest = Append[legacyCore, "ManifestIdentity" ->
  DiffExp2`NativeTransport`Private`nativeCheckpointIdentity[
    "de2-native-observable-checkpoint-manifest-", legacyCore]];
assert["compact checkpoint manifest replaces recursive provenance with digests",
  ByteCount[checkpointManifest] < ByteCount[legacyManifest]/2 &&
    AllTrue[Join[Keys /@ Values[checkpointManifest["StateHandles"]],
        Keys /@ checkpointManifest["Results"]],
      FreeQ[#, "ProvenanceIdentity"] &]];

exportedRun = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[run, 50]],
  run];
exportedResults = If[AssociationQ[exportedRun],
  Lookup[exportedRun, "ExportedResults", {}], {}];
lineValue = If[Length[exportedResults] === 5,
  exportedResults[[1, "Value"]], None];
lowerValue = If[Length[exportedResults] === 5,
  exportedResults[[2, "Value"]], None];
upperValue = If[Length[exportedResults] === 5,
  exportedResults[[3, "Value"]], None];
polarValue = If[Length[exportedResults] === 5,
  exportedResults[[4, "Value"]], None];
secondLineValue = If[Length[exportedResults] === 5,
  exportedResults[[5, "Value"]], None];
integralCertification = If[Length[exportedResults] === 5,
  KeyTake[exportedResults[[1]],
    {"Scope", "ErrorGuarantee", "ErrorEnvelope"}], None];
polarCertification = If[Length[exportedResults] === 5,
  KeyTake[exportedResults[[4]],
    {"Scope", "ErrorGuarantee", "ErrorEnvelope"}], None];
assert["native_observable_batch_contracts_integrals_polar_order_and_endpoints",
  AssociationQ[exportedRun] &&
    Lookup[exportedRun, "CompatibilityExports", 0] === 5 &&
    AllTrue[{lineValue, lowerValue, upperValue, polarValue,
        secondLineValue},
      DiffExp2`EpsSeries`ESQ] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[lineValue, 0] - 1/2,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[lowerValue, 0] - 1,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[upperValue, 0] - 1,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[polarValue, -1] - 1/2,
      30]] < 10^-20] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[secondLineValue, 0] - 1/2,
      30]] < 10^-20]];
assert["native_observable_batch_exports_compact_line_certification_only_for_integrals",
  AssociationQ[integralCertification] &&
    integralCertification["Scope"] ===
      "full_local_with_certified_tail" &&
    integralCertification["ErrorGuarantee"] === "certified" &&
    AssociationQ[integralCertification["ErrorEnvelope"]] &&
    integralCertification["ErrorEnvelope", "guarantee"] ===
      "certified" &&
    ListQ[integralCertification["ErrorEnvelope",
      "absolute_upper_approx"]] &&
    integralCertification["ErrorEnvelope",
      "absolute_upper_approx"] =!= {} &&
    polarCertification === <|"Scope" -> "stored_truncation",
      "ErrorGuarantee" -> "none", "ErrorEnvelope" -> None|> &&
    AllTrue[exportedResults[[{2, 3}]],
      Intersection[Keys[#],
        {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] === {} &]];
malformedCertification = catchDE2[
  DiffExp2`NativeTransport`Private`nativeLineExportCertification[<|
    "scope" -> "full_local_with_certified_tail",
    "error_guarantee" -> "certified",
    "value" -> <|"min" -> 0, "max" -> 0,
      "error" -> <|"min" -> 0, "max" -> 0,
        "guarantee" -> "advisory", "absolute_upper_approx" -> {1.},
        "bound_encoding" -> "approximate-double",
        "provenance" -> "tampered"|>|>|>]];
assert["native_observable_batch_rejects_malformed_certification_metadata",
  FailureQ[malformedCertification]];

released = If[AssociationQ[exportedRun],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[exportedRun],
  exportedRun];
afterRelease = sessionStats[];
assert["native_observable_batch_releases_complete_owner_closure",
  AssociationQ[released] && Lookup[released, "Failures", {"missing"}] === {} &&
    Lookup[afterRelease, "line_results", -1] === 0 &&
    Lookup[afterRelease, "endpoints", -1] === 0 &&
    Lookup[afterRelease, "transport_states", -1] === 0 &&
    Lookup[afterRelease, "locals", -1] === 0 &&
    Lookup[afterRelease, "matches", -1] === 0 &&
    Lookup[afterRelease, "tile_plans", -1] === 0];

tamperedResults = checkpointManifest["Results"];
tamperedResults[[{1, 4}, "Handle"]] =
  Reverse[tamperedResults[[{1, 4}, "Handle"]]];
tamperedCore = Join[KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"Results" -> tamperedResults|>];
tamperedManifest = Append[tamperedCore, "ManifestIdentity" ->
  ("de2-native-observable-checkpoint-manifest-" <>
    IntegerString[Hash[tamperedCore, "SHA256"], 16, 64])];
tamperedRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedManifest]];
assert["checkpoint_restore_rejects_swapped_per_handle_observable_identity",
  FailureQ[tamperedRestore]];

tamperedStates = checkpointManifest["StateHandles"];
tamperedStates["lower", "CheckpointIdentity"] =
  "tampered-lower-state-checkpoint";
tamperedStateCore = Join[KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"StateHandles" -> tamperedStates|>];
tamperedStateManifest = Append[tamperedStateCore, "ManifestIdentity" ->
  ("de2-native-observable-checkpoint-manifest-" <>
    IntegerString[Hash[tamperedStateCore, "SHA256"], 16, 64])];
tamperedStateRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedStateManifest]];
assert["checkpoint_restore_binds_each_transport_state_checkpoint_and_provenance",
  FailureQ[tamperedStateRestore]];

tamperedProvenanceResults = checkpointManifest["Results"];
tamperedProvenanceResults[[1, "ProvenanceSHA256"]] =
  StringRepeat["0", 64];
tamperedProvenanceCore = Join[
  KeyDrop[checkpointManifest, "ManifestIdentity"],
  <|"Results" -> tamperedProvenanceResults|>];
tamperedProvenanceManifest = Append[tamperedProvenanceCore,
  "ManifestIdentity" ->
    DiffExp2`NativeTransport`Private`nativeCheckpointIdentity[
      "de2-native-observable-checkpoint-manifest-",
      tamperedProvenanceCore]];
tamperedProvenanceRestore = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    tamperedProvenanceManifest]];
assert["checkpoint_restore_rejects_recomputed_manifest_with_wrong_provenance_digest",
  FailureQ[tamperedProvenanceRestore]];

restoredRun = If[AssociationQ[checkpointManifest], catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    checkpointManifest]], checkpointManifest];
restoredSession = If[AssociationQ[restoredRun],
  Lookup[restoredRun, "RestoredSession", None], None];
restoredStatsBefore = If[StringQ[restoredSession],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    restoredSession, <||>], <||>];
restoredExport = If[AssociationQ[restoredRun], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    restoredRun, 50]], restoredRun];
restoredStatsAfter = If[StringQ[restoredSession],
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    restoredSession, <||>], <||>];
assert["checkpoint_restore_export_does_not_repeat_transport_arm_marches",
  AssociationQ[restoredExport] &&
    Lookup[restoredStatsBefore, "transport_arm_marches", -1] === 2 &&
    Lookup[restoredStatsAfter, "transport_arm_marches", -2] ===
      Lookup[restoredStatsBefore, "transport_arm_marches", -1] &&
    Lookup[Lookup[restoredExport, "ExportedResults", {}], "Value"] ===
      Lookup[exportedResults, "Value"] &&
    (KeyTake[#, {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] & /@
      Select[restoredExport["ExportedResults"],
        #["Operation"] === "integrate" &]) ===
      (KeyTake[#, {"Scope", "ErrorGuarantee", "ErrorEnvelope"}] & /@
        Select[exportedResults, #["Operation"] === "integrate" &])];
restoredReleased = If[AssociationQ[restoredExport],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    restoredExport], restoredExport];
assert["restored_native_observable_batch_closes_only_its_restored_session",
  AssociationQ[restoredReleased] &&
    Lookup[restoredReleased, "Failures", {"missing"}] === {} &&
    !KeyExistsQ[DiffExp2`CppBackend`PersistentSessionInformation[],
      restoredSession]];

legacyRestored = catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    legacyManifest]];
legacyExport = If[AssociationQ[legacyRestored], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    legacyRestored, 50]], legacyRestored];
assert["checkpoint_restore_remains_backward_compatible_with_v1_full_provenance_manifest",
  AssociationQ[legacyExport] &&
    Lookup[Lookup[legacyExport, "ExportedResults", {}], "Value"] ===
      Lookup[exportedResults, "Value"]];
legacyReleased = If[AssociationQ[legacyExport],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    legacyExport], legacyExport];
assert["legacy_manifest_restore_closes_its_restored_session",
  AssociationQ[legacyReleased] &&
    Lookup[legacyReleased, "Failures", {"missing"}] === {}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];

streamEndpoint = <|"Center" -> 1/4, "Singular" -> False,
  "Radius" -> 1/2, "MatchRadius" -> 1/2, "Scale" -> 1/2,
  "LocalRadius" -> 1, "IncomingMatchPoint" -> 1/8,
  "SymmetricMatch" -> False, "ChartVar" -> Global`t,
  "UseSCCSkeleton" -> True, "Name" -> "stream-endpoint",
  "Prescriptions" -> {}|>;
streamUpperPlan = Join[upperPlan, <|
  "Charts" -> {First[upperPlan["Charts"]], streamEndpoint},
  "SegmentCount" -> 2|>];
ownerAtlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, streamUpperPlan, "Threads" -> 1,
    "Integrand" -> {{1 + Global`eps}, x},
    "TargetCompleteMax" -> 0, "DeferReceivingBases" -> True]];
ownerRun = If[FailureQ[ownerAtlas], ownerAtlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatchOwned[
    ownerAtlas, {observable["limitLower", "owned-lower"]}, x,
    "MaxRefinementSteps" -> 1]]];
ownerCounters = If[AssociationQ[ownerRun],
  DiffExp2`CppBackend`PersistentSessionCounters[ownerRun["Atlas"]],
  ownerRun];
assert["owned deferred observable batch streams and compacts caller atlas",
  AssociationQ[ownerRun] && AssociationQ[ownerAtlas] &&
    !KeyExistsQ[ownerAtlas["Lower"], "ChartSystems"] &&
    !KeyExistsQ[ownerAtlas["Upper"], "ChartSystems"] &&
    !KeyExistsQ[ownerRun["Atlas", "Lower"], "ChartSystems"] &&
    !KeyExistsQ[ownerRun["Atlas", "Upper"], "ChartSystems"]];
assert["owned deferred stream measures physical_selection_and_releases_owner",
  AssociationQ[ownerCounters] &&
    Lookup[ownerCounters,
      "transport_physical_value_hop_attempts", -1] === 1 &&
    Lookup[ownerCounters, "transport_physical_value_hop_successes", -1] +
      Lookup[ownerCounters,
        "transport_physical_value_hop_ineligible", -1] === 1 &&
    Lookup[ownerCounters, "transport_framed_basis_hops", -1] ===
      Lookup[ownerCounters,
        "transport_physical_value_hop_ineligible", -2] &&
    Lookup[ownerCounters, "regular_equation_owners", -1] === 0 &&
    !KeyExistsQ[ownerRun["Atlas", "Lower"], "OwnerRecords"] &&
    !KeyExistsQ[ownerRun["Atlas", "Upper"], "OwnerRecords"]];
If[AssociationQ[ownerRun],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[ownerRun]];

releasedOwnerAtlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, streamUpperPlan, "Threads" -> 1,
    "Integrand" -> {{1 + Global`eps}, x},
    "TargetCompleteMax" -> 0, "DeferReceivingBases" -> True]];
releasedOwnerSession = If[AssociationQ[releasedOwnerAtlas],
  releasedOwnerAtlas["Session"], None];
releasedOwnerResponse = If[AssociationQ[releasedOwnerAtlas],
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    releasedOwnerAtlas], releasedOwnerAtlas];
releasedOwnerCounters = If[StringQ[releasedOwnerSession],
  DiffExp2`CppBackend`PersistentSessionCounters[releasedOwnerSession],
  releasedOwnerResponse];
assert["unexecuted deferred atlas release drops_public_equation_owner",
  AssociationQ[releasedOwnerResponse] &&
    Lookup[releasedOwnerResponse, "Failures", {"missing"}] === {} &&
    AssociationQ[releasedOwnerCounters] &&
    Lookup[releasedOwnerCounters, "regular_equation_owners", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
