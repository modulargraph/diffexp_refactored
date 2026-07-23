(* Focused production-runner test for the retained native observable cutover.
   Loads definitions only and replaces the six native seams, so no FIRE job,
   chart solve, or persistent C++ session is started. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> "Cpp"];

Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label, If[detail === None, "", ": "],
    If[detail === None, "", detail]]];

assert["downstream publication digits remain independent of guarded producer digits",
  ft2DownstreamPublicationDigits[
      5, <|4 -> 8, 5 -> 10, 6 -> 12|>] === 8 &&
    ft2DownstreamPublicationDigits[
      1, <|1 -> 10, 2 -> 12|>] === matchDigits];
assert["internal matching certification consumes at most one downstream safety guard",
  ft2LevelMatchingCertificationDigits[12, 8] === 10 &&
    ft2LevelMatchingCertificationDigits[14, 12] === 14 &&
    ft2LevelMatchingCertificationDigits[8, 8] === 8];

reservoirBackendFailure = <|
  "reason" -> "acb_match_residual_inconclusive",
  "retryable_epsilon_reservoir" -> True,
  "retryable_matching_clearance" -> False,
  "required_additional_epsilon_orders" -> 3|>;
clearanceBackendFailure = <|
  "reason" -> "acb_match_residual_inconclusive",
  "retryable_epsilon_reservoir" -> False,
  "retryable_matching_clearance" -> True,
  "required_additional_epsilon_orders" -> 0,
  "residual" -> <|"complete_through_required" -> True,
    "scope" -> "stored-taylor-truncation",
    "coefficient_verdicts" -> <|
      "pass" -> 113, "fail" -> 0, "inconclusive" -> 139|>,
    "required_coefficient_verdicts" -> <|
      "pass" -> 100, "fail" -> 0, "inconclusive" -> 5|>|>|>;
reservoirFailure = Failure["DiffExp2", <|
  "BackendFailure" -> reservoirBackendFailure|>];
clearanceFailure = Failure["DiffExp2", <|
  "BackendFailure" -> clearanceBackendFailure|>];
propagatedBackendFailure = Join[clearanceBackendFailure, <|
  "retryable_matching_clearance" -> False,
  "retryable_propagated_enclosure" -> True,
  "normal_frame_attempt" -> <|
    "physical_clearance_source" -> "propagated-enclosure",
    "physical_clearance_source_probes" -> <|
      "basis" -> <|"verdict" -> "pass"|>,
      "weights" -> <|"verdict" -> "inconclusive"|>,
      "incoming" -> <|"verdict" -> "inconclusive"|>|>|>|>];
propagatedFailure = Failure["DiffExp2", <|
  "BackendFailure" -> propagatedBackendFailure|>];
continuityBackendFailure = <|
  "reason" -> "acb_match_residual_inconclusive",
  "retryable_epsilon_reservoir" -> False,
  "retryable_matching_clearance" -> True,
  "required_additional_epsilon_orders" -> 0,
  "residual" -> <|
    "status" -> "materialized-continuity-inconclusive",
    "complete_through_required" -> True,
    "scope" -> "materialized-continuity-clearance",
    "coefficient_verdicts" -> <|
      "pass" -> 0, "fail" -> 0, "inconclusive" -> 1|>,
    "required_coefficient_verdicts" -> <|
      "pass" -> 0, "fail" -> 0, "inconclusive" -> 1|>|>|>;
continuityFailure = Failure["DiffExp2", <|
  "BackendFailure" -> continuityBackendFailure|>];
reservoirRetry = ft2NativeMatchingReservoirRetry[reservoirFailure, 3];
clearanceRetry = ft2NativeMatchingClearanceRetry[
  clearanceFailure, 3, 50];
continuityRetry = ft2NativeMatchingClearanceRetry[
  continuityFailure, 1, 50];
assert["complete matching clearance retries Taylor order, not epsilon width",
  ft2NativeMatchingReservoirRetry[clearanceFailure, 3] === None &&
    ft2NativeMatchingClearanceRetry[
      reservoirFailure, 3, 50] === None &&
    ft2NativeMatchingReservoirRetryQ[reservoirRetry] &&
    reservoirRetry[[2, "AdditionalOrders"]] === 3 &&
    ft2NativeMatchingClearanceRetryQ[clearanceRetry] &&
    ft2NativeMatchingRetryQ[clearanceRetry] &&
    clearanceRetry[[2, "Level"]] === 3 &&
    clearanceRetry[[2, "CurrentExpansionOrder"]] === 50 &&
    clearanceRetry[[2, "AdditionalOrders"]] === 50 &&
    clearanceRetry[[2, "ResidualVerdicts", "inconclusive"]] === 5,
  {reservoirRetry, clearanceRetry}];
assert["materialized handoff clearance permits one failing-level Taylor probe",
  ft2NativeMatchingReservoirRetry[continuityFailure, 1] === None &&
    ft2NativeMatchingClearanceRetryQ[continuityRetry] &&
    continuityRetry[[2, "Level"]] === 1 &&
    continuityRetry[[2, "CurrentExpansionOrder"]] === 50 &&
    continuityRetry[[2, "AdditionalOrders"]] === 50 &&
    continuityRetry[[2, "ResidualVerdicts", "inconclusive"]] === 1,
  continuityRetry];
assert["producer enclosure is not misclassified as a Taylor-order retry",
  ft2NativeMatchingReservoirRetry[propagatedFailure, 3] === None &&
    ft2NativeMatchingClearanceRetry[
      propagatedFailure, 3, 50] === None,
  {ft2NativeMatchingReservoirRetry[propagatedFailure, 3],
   ft2NativeMatchingClearanceRetry[propagatedFailure, 3, 50]}];
producerRetry = ft2NativeMatchingProducerRetry[
  propagatedFailure, 3, 4, 8];
assert["propagated enclosure retries only the immediately preceding producer",
  ft2NativeMatchingProducerRetryQ[producerRetry] &&
    ft2NativeMatchingRetryQ[producerRetry] &&
    producerRetry[[2, "Level"]] === 3 &&
    producerRetry[[2, "ProducerLevel"]] === 4 &&
    producerRetry[[2, "NumLevels"]] === 4 &&
    producerRetry[[2, "CurrentMatchingDigits"]] === 8 &&
    producerRetry[[2, "AdditionalOrders"]] ===
      DiffExp2`Tolerances`$SafetyDigits &&
    ft2NativeMatchingProducerRetry[
      propagatedFailure, 4, 4, 8] === None &&
    ft2NativeMatchingProducerRetry[
      clearanceFailure, 3, 4, 8] === None,
  producerRetry];
weightOnlyPropagatedFailure = Failure["DiffExp2", <|
  "BackendFailure" -> ReplacePart[propagatedBackendFailure,
    {"normal_frame_attempt", "physical_clearance_source_probes",
      "incoming", "verdict"} -> "pass"]|>];
assert["weight-only propagated enclosure is not sent back to the producer",
  ft2NativeMatchingProducerRetry[
    weightOnlyPropagatedFailure, 3, 4, 8] === None];
algebraicWeightPropagatedFailure = Failure["DiffExp2", <|
  "BackendFailure" -> Join[
    weightOnlyPropagatedFailure[[2, "BackendFailure"]], <|
      "weight_shadow_retry" ->
        "upstream-accuracy-required-nonrational-chart"|>]|>];
assert["nonrational weight enclosure uses guarded upstream accuracy",
  ft2NativeMatchingProducerRetryQ[
    ft2NativeMatchingProducerRetry[
      algebraicWeightPropagatedFailure, 3, 4, 8]]];
assert["producer precision is bound to producer checkpoints, not consumers",
  ft2CheckpointExpectedMatchingDigits[
      "Boundary", 3, <|4 -> 10|>] === 10 &&
    ft2CheckpointExpectedMatchingDigits[
      "NativeTransport", 3, <|4 -> 10|>] === matchDigits &&
    ft2CheckpointExpectedMatchingDigits[
      "NativeTransport", 4, <|4 -> 10|>] === 10];
assert["checkpoint publication identity follows the downstream consumer",
  ft2CheckpointExpectedPublicationDigits[
      "Boundary", 3, <|3 -> 10, 4 -> 12|>] === 10 &&
    ft2CheckpointExpectedPublicationDigits[
      "NativeTransport", 4, <|3 -> 10, 4 -> 12|>] === 10 &&
    ft2CheckpointExpectedPublicationDigits[
      "NativeTransport", 3, <|3 -> 10, 4 -> 12|>] === matchDigits];
assert["producer digit budget retains its safety margin up the ladder",
  ft2RaiseMatchingProducerDigits[
      <|4 -> matchDigits + 2|>, 3, matchDigits + 2, 4] ===
    <|4 -> matchDigits + 4, 3 -> matchDigits + 2|> &&
    ft2RaiseMatchingProducerDigits[
      <|3 -> matchDigits + 2, 4 -> matchDigits + 4|>,
      3, matchDigits + 4, 4] ===
        <|3 -> matchDigits + 4, 4 -> matchDigits + 6|>];
assert["banana4 final producer seed carries one measured digit plus guards",
  ft2RaiseMatchingProducerDigits[
      <||>, 2, matchDigits + 1, 4] ===
    <|2 -> matchDigits + 1, 3 -> matchDigits + 3,
      4 -> matchDigits + 5|>];

terminalOutputFailure = Failure["DiffExp2", <|
  "BackendFailure" -> <|
    "reason" -> "terminal_output_ball_inconclusive",
    "retryable_level_accuracy" -> True,
    "failure_functional" -> 0, "failure_epsilon" -> 2|>|>];
terminalOutputRetry = ft2NativeTerminalOutputProducerRetry[
  terminalOutputFailure, 5, 6, 10];
measuredTerminalOutputFailure = Failure["DiffExp2", <|
  "BackendFailure" -> Join[
    terminalOutputFailure[[2, "BackendFailure"]],
    <|"required_additional_digits" -> 1|>]|>];
measuredTerminalOutputRetry = ft2NativeTerminalOutputProducerRetry[
  measuredTerminalOutputFailure, 5, 6, 10];
assert[
  "wide terminal boundary retries its own level as the next consumer's producer",
  ft2NativeMatchingProducerRetryQ[terminalOutputRetry] &&
    terminalOutputRetry[[2, "Level"]] === 4 &&
    terminalOutputRetry[[2, "ProducerLevel"]] === 5 &&
    terminalOutputRetry[[2, "NumLevels"]] === 6 &&
    terminalOutputRetry[[2, "CurrentMatchingDigits"]] === 10 &&
    terminalOutputRetry[[2, "AdditionalOrders"]] ===
      DiffExp2`Tolerances`$SafetyDigits &&
    measuredTerminalOutputRetry[[2, "AdditionalOrders"]] === 1 &&
    ft2NativeTerminalOutputProducerRetry[
      terminalOutputFailure, 1, 6, 10] === None,
  <|"fallback" -> terminalOutputRetry,
    "measured" -> measuredTerminalOutputRetry|>];

assert["matching Taylor progress requires fewer inconclusive coefficients",
  ft2NativeMatchingClearanceProgressQ[
      <|"pass" -> 113, "fail" -> 0, "inconclusive" -> 139|>,
      <|"pass" -> 133, "fail" -> 0, "inconclusive" -> 119|>] &&
    !ft2NativeMatchingClearanceProgressQ[
      <|"pass" -> 113, "fail" -> 0, "inconclusive" -> 139|>,
      <|"pass" -> 113, "fail" -> 0, "inconclusive" -> 139|>]];

matchingTaylorRetryCalls = {};
matchingTaylorRetryResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[matchingTaylorRetryCalls, taylorOrders];
      If[Length[matchingTaylorRetryCalls] === 1,
        Failure["FeynmanTrickNativeMatchingTaylor", <|
          "Level" -> 3, "AdditionalOrders" -> expansionOrder,
          "CurrentExpansionOrder" -> expansionOrder,
          "ResidualVerdicts" -> <|
            "pass" -> 10, "fail" -> 0, "inconclusive" -> 4|>,
          "MatchingTaylorOrders" -> taylorOrders|>], True]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries["matching-taylor-retry-fixture"]];
assert["matching retry driver doubles only the failing level Taylor order",
  matchingTaylorRetryResult === True &&
    matchingTaylorRetryCalls ===
      {<||>, <|3 -> 2 expansionOrder|>},
  matchingTaylorRetryCalls];

matchingTaylorStallCalls = {};
matchingTaylorStallResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[matchingTaylorStallCalls, taylorOrders];
      Failure["FeynmanTrickNativeMatchingTaylor", <|
        "Level" -> 2, "AdditionalOrders" -> expansionOrder,
        "CurrentExpansionOrder" -> Lookup[taylorOrders, 2, expansionOrder],
        "ResidualVerdicts" -> <|
          "pass" -> 113, "fail" -> 0, "inconclusive" -> 139|>,
        "MatchingTaylorOrders" -> taylorOrders|>]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries["matching-taylor-stall-fixture"]];
assert["matching retry stops after a Taylor order increase makes no progress",
  FailureQ[matchingTaylorStallResult] &&
    matchingTaylorStallResult[[1]] ===
      "FeynmanTrickNativeMatchingTaylorStalled" &&
    matchingTaylorStallCalls === {<||>, <|2 -> 2 expansionOrder|>},
  {matchingTaylorStallResult, matchingTaylorStallCalls}];

continuityTaylorStallCalls = {};
continuityTaylorStallResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[continuityTaylorStallCalls, taylorOrders];
      Failure["FeynmanTrickNativeMatchingTaylor", <|
        "Level" -> 1, "AdditionalOrders" -> expansionOrder,
        "CurrentExpansionOrder" -> Lookup[taylorOrders, 1, expansionOrder],
        "ResidualVerdicts" -> <|
          "pass" -> 0, "fail" -> 0, "inconclusive" -> 1|>,
        "MatchingTaylorOrders" -> taylorOrders,
        "BackendFailure" -> continuityBackendFailure|>]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries[
    "materialized-continuity-stall-fixture"]];
assert[
  "unchanged materialized continuity clearance stops after exactly one probe",
  FailureQ[continuityTaylorStallResult] &&
    continuityTaylorStallResult[[1]] ===
      "FeynmanTrickNativeMatchingTaylorStalled" &&
    continuityTaylorStallCalls === {<||>, <|1 -> 2 expansionOrder|>} &&
    continuityTaylorStallResult[[2, "BackendFailure", "residual",
      "scope"]] === "materialized-continuity-clearance",
  {continuityTaylorStallResult, continuityTaylorStallCalls}];

matchingProducerRetryCalls = {};
matchingProducerRetryResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[matchingProducerRetryCalls, levelDigits];
      If[Length[matchingProducerRetryCalls] === 1,
        Failure["FeynmanTrickNativeMatchingProducer", <|
          "Level" -> 3, "ProducerLevel" -> 4,
          "NumLevels" -> 4,
          "CurrentMatchingDigits" -> matchDigits,
          "AdditionalOrders" -> DiffExp2`Tolerances`$SafetyDigits,
          "MatchingDigitsByLevel" -> levelDigits|>], True]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries["matching-producer-retry-fixture"]];
assert["matching retry driver tightens only the preceding producer level",
  matchingProducerRetryResult === True &&
    matchingProducerRetryCalls ===
      {<||>, <|4 -> matchDigits +
        DiffExp2`Tolerances`$SafetyDigits|>},
  matchingProducerRetryCalls];

seededRetryCalls = {};
seededRetryResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[seededRetryCalls,
        {epsilonHalos, taylorOrders, levelDigits}];
      If[Length[seededRetryCalls] === 1,
        Failure["FeynmanTrickNativeMatchingProducer", <|
          "Level" -> 3, "ProducerLevel" -> 4,
          "NumLevels" -> 4,
          "CurrentMatchingDigits" -> matchDigits,
          "AdditionalOrders" -> DiffExp2`Tolerances`$SafetyDigits,
          "MatchingDigitsByLevel" -> levelDigits|>], True]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries[
    "seeded-matching-retry-fixture", None, <|2 -> 5|>,
    <|1 -> 75|>, <||>]];
assert[
  "seeded timing-gate halos survive producer-accuracy retries",
  seededRetryResult === True &&
    seededRetryCalls === {
      {<|2 -> 5|>, <|1 -> 75|>, <||>},
      {<|2 -> 5|>, <|1 -> 75|>,
        <|4 -> matchDigits +
          DiffExp2`Tolerances`$SafetyDigits|>}},
  seededRetryCalls];

preapprovedHighDigitCalls = {};
preapprovedHighDigitResult = Block[{
    runExample = Function[
      {runName, familyRequest, epsilonHalos, taylorOrders, levelDigits},
      AppendTo[preapprovedHighDigitCalls, levelDigits];
      If[Length[preapprovedHighDigitCalls] === 1,
        Failure["FeynmanTrickNativeMatchingProducer", <|
          "Level" -> 2, "ProducerLevel" -> 3,
          "NumLevels" -> 6,
          "CurrentMatchingDigits" -> matchDigits + 2,
          "AdditionalOrders" -> 1,
          "MatchingDigitsByLevel" -> levelDigits|>], True]],
    DiffExp2`Solve`ClearSolveCaches = Function[{}, Null]},
  ft2RunExampleWithMatchingRetries[
    "preapproved-high-digit-retry-fixture", None, <||>, <||>,
    <|3 -> matchDigits + 2, 4 -> matchDigits + 8,
      5 -> matchDigits + 11, 6 -> matchDigits + 13|>]];
assert[
  "retry ceiling checks only newly raised digits, not preapproved seeds",
  preapprovedHighDigitResult === True &&
    preapprovedHighDigitCalls === {
      <|3 -> matchDigits + 2, 4 -> matchDigits + 8,
        5 -> matchDigits + 11, 6 -> matchDigits + 13|>,
      <|3 -> matchDigits + 3, 4 -> matchDigits + 8,
        5 -> matchDigits + 11, 6 -> matchDigits + 13|>},
  {preapprovedHighDigitResult, preapprovedHighDigitCalls}];

x = Global`xNativeFT;
epsilon = Global`eps;
request[index_, case_, vi_, vj_, needed_] := <|
  "MasterIndex" -> index, "MasterVec" -> {vi, vj},
  "Case" -> case, "Vi" -> vi, "Vj" -> vj,
  "NeededVec" -> needed|>;
requests = {
  request[1, "integrate", 2, 3, {1}],
  request[2, "limitLower", 0, 1, {2}],
  request[3, "limitUpper", 1, 0, {3}],
  request[4, "direct", 0, 0, {4}],
  request[5, "integrate", 1, 1, {5}],
  request[6, "integrate", 1, 1, {6}]};
vectors = Association[{
  {1} -> {2, 3/epsilon},
  {2} -> {epsilon^3, epsilon},
  {3} -> {1, epsilon^2},
  {4} -> {epsilon, 2},
  {5} -> {0, 0},
  {6} -> {epsilon^21, epsilon^20}}];
batch = <|"Schema" -> "FeynmanTrick.LevelIBPBatch/v1",
  "Key" -> StringRepeat["a", 64],
  "PayloadKey" -> StringRepeat["b", 64],
  "KeyRecord" -> {"synthetic"}, "UpperLevel" -> 3,
  "MastersAbove" -> {{1, 0}, {0, 1}},
  "BoundaryRequests" -> requests,
  "CoefficientVectors" -> vectors|>;
normalizeIdentity = Function[value, value];
entries = ft2PrepareBoundaryEntries[
  3, batch, {1, 0}, x, epsilon, normalizeIdentity];
badCaseRequests = ReplacePart[requests, 2 ->
  Join[requests[[2]], <|"Case" -> "limitUpper"|>]];
badIndexRequests = ReplacePart[requests, 3 ->
  Join[requests[[3]], <|"MasterIndex" -> 99|>]];
assert["batch schema, key, master index, and case contracts fail loudly",
  FailureQ[ft2PrepareBoundaryEntries[3,
      Join[batch, <|"Schema" -> "wrong"|>],
      {1, 0}, x, epsilon, normalizeIdentity]] &&
    FailureQ[ft2PrepareBoundaryEntries[3,
      Join[batch, <|"BoundaryRequests" -> badCaseRequests|>],
      {1, 0}, x, epsilon, normalizeIdentity]] &&
    FailureQ[ft2PrepareBoundaryEntries[3,
      Join[batch, <|"BoundaryRequests" -> badIndexRequests|>],
      {1, 0}, x, epsilon, normalizeIdentity]]];

expectedWeight = x*(1 - x)^2;
assert["coefficient vectors include basis, beta, and x weights exactly",
  ListQ[entries] && Length[entries] === 6 &&
    And @@ MapThread[TrueQ[PossibleZeroQ[Together[#1 - #2]]] &,
      {entries[[1, "CoefficientVector"]],
       {24 expectedWeight/epsilon, 36 expectedWeight/epsilon}}] &&
    entries[[2, "CoefficientVector"]] === {epsilon^2, epsilon} &&
    entries[[4, "CoefficientVector"]] === {1, 2}];
assert["all-zero coefficient vectors are proved before halo accounting",
  TrueQ[entries[[5, "ProvenZero"]]] &&
    !TrueQ[entries[[1, "ProvenZero"]]] &&
    !TrueQ[entries[[6, "ProvenZero"]]]];

sourceRows = {
  Join[{1}, ConstantArray[0, 10]],
  Join[{2}, ConstantArray[0, 10]]};
ledger = ft2NativeEpsilonLedger[entries, sourceRows, 7, 3];
assert["epsilon ledger separates coefficient and integration halos",
  AssociationQ[ledger] &&
    KeyTake[ledger, {"AvailableSourceCompleteMax", "SourceCompleteMax",
      "CoefficientHalo",
      "IntegrationHalo", "PublicTargetCompleteMax",
      "TargetCompleteMax",
      "DeliverableCompleteMax",
      "DownstreamRawTop"}] ===
      <|"AvailableSourceCompleteMax" -> 10,
        "SourceCompleteMax" -> 10, "CoefficientHalo" -> 1,
        "IntegrationHalo" -> 1, "PublicTargetCompleteMax" -> 4,
        "TargetCompleteMax" -> 8, "DeliverableCompleteMax" -> 7,
        "DownstreamRawTop" -> 3|>,
  ledger];
assert["observable minima do not discount independently required raw depth",
  ledger["OutputMinimums"] ===
    <|1 -> -2, 2 -> 1, 3 -> -1, 4 -> 0, 6 -> 19|> &&
    !KeyExistsQ[ledger["OutputMinimums"], 5],
  ledger["OutputMinimums"]];
terminalLedger = ft2NativeEpsilonLedger[entries, sourceRows, 0, 0];
assert["terminal raw completeness is exactly epsOrder, independent of poles",
  AssociationQ[terminalLedger] &&
    terminalLedger["SourceCompleteMax"] === 10 &&
    terminalLedger["TargetCompleteMax"] === 1 &&
    terminalLedger["DeliverableCompleteMax"] === 0 &&
    terminalLedger["DownstreamRawTop"] === 0 &&
    terminalLedger["DeliverableCompleteMax"] >=
      terminalLedger["DownstreamRawTop"],
  terminalLedger];
runnerSource = Import[
  FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}], "Text"];
assert["native floor excludes the consumable boundary-extra reservoir",
  nativeRequiredRawTop[0] === epsOrder &&
    nativeRequiredRawTop[2] ===
      Max[epsOrder + 2 + levelEpsilonHalo[2], 1] &&
    requestedEpsilonOrder[2] >= nativeRequiredRawTop[2]];
assert["Cpp handoff retains the common certified raw edge and records its width",
  StringContainsQ[runnerSource,
      "needTop = Min[Min[kmaxAvail], nextReq + boundaryExtraOrder]"] &&
    StringContainsQ[runnerSource,
      "\"PreservedRawCompleteMax\" -> If[recurrenceBackend === \"Cpp\""] &&
    StringContainsQ[runnerSource,
      "\"PreservedSourceCompleteMax\" -> If["]];

newCounts[] := <|"segment" -> 0, "prepare" -> 0, "run" -> 0,
  "export" -> 0, "releaseBatch" -> 0, "releaseAtlas" -> 0|>;
fixtureCertifiedEnvelope = <|"guarantee" -> "certified",
  "absolute_upper_approx" -> {1.*^-30},
  "bound_encoding" -> "approximate-double",
  "provenance" -> "definitions-only-fixture"|>;
fixtureExportResults[results_List] := MapIndexed[Function[{observable, pos},
  Module[{exported = Append[observable,
      "Value" -> DiffExp2`EpsSeries`ESZero[
        observable["Epsilon", "Max"]]]},
    Switch[observable["Operation"],
      "integrate",
        If[First[pos] === 1,
          Join[exported, <|
            "Scope" -> "full_local_with_certified_tail",
            "ErrorGuarantee" -> "certified",
            "ErrorEnvelope" -> fixtureCertifiedEnvelope|>],
          Join[exported, <|"Scope" -> "stored_truncation",
            "ErrorGuarantee" -> "none", "ErrorEnvelope" -> None|>]],
      _, exported]]], results];
counts = newCounts[];
capturedPrepare = None;
capturedObservables = None;
capturedMatchingCertificationDigits = None;
capturedPublicationDigits = None;
capturedMaxExtraPrecision = None;

mixed = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; <|"Path" -> path|>],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, requiredTargetMax, threads},
      counts["prepare"]++;
      capturedMaxExtraPrecision = $MaxExtraPrecision;
      capturedPrepare = <|"Boundary" -> boundary,
        "CoefficientVectors" -> coefficientVectors,
        "TargetCompleteMax" -> targetMax,
        "RequiredTargetCompleteMax" -> requiredTargetMax,
        "Threads" -> threads|>;
      <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
        "PlanCheckpointIdentity" -> "synthetic-atlas-plan"|>],
    ft2NativeRun = Function[
      {atlas, observables, variable, matchingCertificationDigits,
       publicationDigits},
      counts["run"]++; capturedObservables = observables;
      capturedMatchingCertificationDigits =
        matchingCertificationDigits;
      capturedPublicationDigits = publicationDigits;
      <|"Type" -> "DiffExp2NativeTransportObservableBatch",
        "Atlas" -> atlas, "NativeMarches" -> 2,
        "Results" -> observables|>],
    ft2NativeExport = Function[{nativeBatch, digits},
      counts["export"]++;
      Join[nativeBatch, <|"CompatibilityExports" ->
          Length[nativeBatch["Results"]],
        "ExportedResults" ->
          fixtureExportResults[nativeBatch["Results"]]|>]],
    ft2NativeReleaseBatch = Function[nativeBatch,
      counts["releaseBatch"]++;
      <|"Released" -> 1, "Failures" -> {}|>],
    ft2NativeReleaseAtlas = Function[atlas,
      counts["releaseAtlas"]++;
      <|"Released" -> 1, "Failures" -> {}|>]},
  ft2RunNativeBoundaryDispatch[
    <|"Matrix" -> IdentityMatrix[2], "Variable" -> x|>,
    sourceRows, entries, ledger, x, 11/23, {x, 1 - x},
    {{x, 1}, {1 - x, 1}}, 6, 50, 8, 8]];

assert["mixed dispatch plans both arms but invokes one native batch/export",
  AssociationQ[mixed] &&
    counts === <|"segment" -> 2, "prepare" -> 1, "run" -> 1,
      "export" -> 1, "releaseBatch" -> 1, "releaseAtlas" -> 0|> &&
    mixed["NativeBatchCalls"] === 1 &&
    mixed["NativeMarches"] === 2 &&
    mixed["CompatibilityExports"] === 4 &&
    capturedMatchingCertificationDigits === 8 &&
    capturedPublicationDigits === 8,
  {counts, mixed}];
integrateObservables = Select[capturedObservables,
  # ["Operation"] === "integrate" &];
limitObservables = Select[capturedObservables,
  MemberQ[{"limitLower", "limitUpper"}, # ["Operation"]] &];
divergentPolicies = Lookup[integrateObservables,
  "DivergentCancellation", Missing["NotFound"]];
divergentProvenance = Quiet[Check[
    ImportString[#, "RawJSON"] & /@
      Lookup[divergentPolicies, "Provenance"], $Failed]];
defaultMatchTolExact = ToString[
  DiffExp2`Tolerances`Tol["MatchTol"], InputForm];
lowTargetCancellationPolicy = Block[{
    DiffExp2`Tolerances`Tol = Function[key, Switch[key,
      "LaurentLeadTol", 10^-24, "MatchTol", 10^-8,
      _, $Failed]]}, ft2DivergentCancellationPolicy[]];
lowTargetCancellationProvenance = Quiet[Check[ImportString[
    lowTargetCancellationPolicy["Provenance"], "RawJSON"], $Failed]];
assert["atlas preparation receives every active non-direct vector and exact target",
  AssociationQ[capturedPrepare] &&
    capturedMaxExtraPrecision >=
      DiffExp2`Tolerances`$MaxExtraPrecisionValue &&
    capturedPrepare["TargetCompleteMax"] === 8 &&
    capturedPrepare["RequiredTargetCompleteMax"] === 4 &&
    Length[capturedPrepare["CoefficientVectors"]] === 4 &&
    AllTrue[capturedPrepare["Boundary"],
      DiffExp2`EpsSeries`ESMinPower[#] === -1 &&
        DiffExp2`EpsSeries`ESCompleteMax[#] === 10 &],
  {capturedMaxExtraPrecision, capturedPrepare}];
assert["observable order and operation-specific epsilon windows are stable",
  Lookup[capturedObservables, "Operation"] ===
    {"integrate", "limitLower", "limitUpper", "integrate"} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Min"] ===
      {-2, 1, -1, 7} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Max"] ===
      {7, 7, 7, 7} &&
    Lookup[Lookup[capturedObservables, "Epsilon"],
      "RequiredCompleteMax"] === {7, 7, 7, 7} &&
    Lookup[Select[capturedObservables,
        #["Operation"] === "integrate" &], "TailPolicy"] ===
      {"stored", "stored"} &&
    Keys[capturedObservables[[1]]] ===
      {"Operation", "Identity", "CheckpointIdentity",
       "CoefficientVector", "Epsilon", "TailPolicy",
       "DivergentCancellation"}];
assert["FT integrate observables alone carry one explicit bounded cancellation policy",
  Length[divergentPolicies] === 2 &&
    Lookup[divergentPolicies, "Mode"] ===
      {"bounded-relative-acb", "bounded-relative-acb"} &&
    SameQ @@ Lookup[divergentPolicies, "RelativeTolerance"] &&
    StringQ[divergentPolicies[[1, "RelativeTolerance"]]] &&
    StringMatchQ[divergentPolicies[[1, "RelativeTolerance"]],
      NumberString ~~ "e-24"] &&
    SameQ @@ Lookup[divergentPolicies, "Provenance"] &&
    FreeQ[limitObservables, "DivergentCancellation"] &&
    ListQ[divergentProvenance] &&
    divergentProvenance === ConstantArray[<|
        "schema" -> "feynman-trick-divergent-cancellation-v2",
        "producer" -> "DiffExp2`Tolerances`Tol",
        "formula" ->
          "Max[LaurentLeadTol,MatchTol/10^SafetyDigits]",
        "laurent_lead_tol_exact" ->
          "1/1000000000000000000000000",
        "match_tol_exact" -> defaultMatchTolExact,
        "safety_digits" -> 2,
        "effective_exact_value" ->
          "1/1000000000000000000000000"|>, 2],
  {divergentPolicies, divergentProvenance, limitObservables}];
assert["low matching targets retain two guarded digits beyond the requested result",
  AssociationQ[lowTargetCancellationPolicy] &&
    StringMatchQ[lowTargetCancellationPolicy["RelativeTolerance"],
      NumberString ~~ "e-10"] &&
    AssociationQ[lowTargetCancellationProvenance] &&
    lowTargetCancellationProvenance["effective_exact_value"] ===
      "1/10000000000" &&
    lowTargetCancellationProvenance["match_tol_exact"] ===
      "1/100000000" &&
    lowTargetCancellationProvenance["safety_digits"] === 2,
  {lowTargetCancellationPolicy, lowTargetCancellationProvenance}];
assert["high-shift integration is marched rather than unsafely pruned",
  capturedObservables[[-1, "Identity"]] === entries[[6, "Identity"]] &&
    capturedObservables[[-1, "Epsilon", "Min"]] === 7];
assert["direct and proven-zero results merge in original master order",
  Length[mixed["Values"]] === 6 &&
    DiffExp2`EpsSeries`ESCoefficient[mixed["Values"][[4]], 0] === 5 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[4]]] === 7 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[5]]] === 7 &&
    DiffExp2`EpsSeries`ESMinPower[mixed["Values"][[1]]] === 7,
  DiffExp2`EpsSeries`ESWindow /@ mixed["Values"]];
assert["dispatch preserves per-master certification and explicit non-applicability",
  Length[mixed["Certifications"]] === 6 &&
    mixed["Certifications"][[1]] === <|
      "Applicability" -> "applicable", "Operation" -> "integrate",
      "Scope" -> "full_local_with_certified_tail",
      "ErrorGuarantee" -> "certified",
      "ErrorEnvelope" -> fixtureCertifiedEnvelope|> &&
    mixed["Certifications"][[2]] ===
      ft2NotApplicableCertification["limitLower", "endpoint-limit"] &&
    mixed["Certifications"][[3]] ===
      ft2NotApplicableCertification["limitUpper", "endpoint-limit"] &&
    mixed["Certifications"][[4]] ===
      ft2NotApplicableCertification["direct", "direct"] &&
    mixed["Certifications"][[5]] ===
      ft2NotApplicableCertification["integrate", "proven-zero"] &&
    mixed["Certifications"][[6]] === <|
      "Applicability" -> "applicable", "Operation" -> "integrate",
      "Scope" -> "stored_truncation", "ErrorGuarantee" -> "none",
      "ErrorEnvelope" -> Null|>,
  mixed["Certifications"]];
assert["checkpoint audit record binds requests, coefficients, prescriptions, atlas, and payload",
  ft2NativeCheckpointRecordQ[mixed["CheckpointRecord"]] &&
    mixed["CheckpointRecord", "AtlasPlanIdentity"] ===
      "synthetic-atlas-plan" &&
    Length[mixed["CheckpointRecord", "RequestIdentities"]] === 6 &&
    Length[mixed["CheckpointRecord", "CoefficientIdentities"]] === 6,
  mixed["CheckpointRecord"]];

checkpointEvents = {};
publishedResume = None;
publishedAudit = None;
checkpointContractIdentity = "synthetic-native-contract";
checkpointSpec = <|"Mode" -> "Save",
  "Path" -> FileNameJoin[{$TemporaryDirectory,
    "synthetic-native-state.checkpoint"}],
  "ContractIdentity" -> checkpointContractIdentity,
  "Publish" -> Function[{resumeRecord, auditRecord},
    AppendTo[checkpointEvents, "publish"];
    publishedResume = resumeRecord; publishedAudit = auditRecord;
    "published-before-export"]|>;
checkpointed = Block[{
    ft2NativeSegmentLine = Function[{system, path}, <|"Path" -> path|>],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, requiredTargetMax, threads},
      <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
        "PlanCheckpointIdentity" -> "checkpoint-atlas-plan"|>],
    ft2NativeRun = Function[
      {atlas, observables, variable, matchingCertificationDigits,
       publicationDigits},
      AppendTo[checkpointEvents, "run"];
      <|"Type" -> "DiffExp2NativeTransportObservableBatch",
        "Atlas" -> atlas, "NativeMarches" -> 2,
        "Results" -> observables|>],
    ft2NativeSaveCheckpoint = Function[{nativeBatch, path, identity},
      AppendTo[checkpointEvents, "save"];
      <|"CheckpointIdentity" -> identity,
        "TransportArmMarches" -> 2|>],
    ft2NativeExport = Function[{nativeBatch, digits},
      AppendTo[checkpointEvents, "export"];
      Join[nativeBatch, <|"CompatibilityExports" ->
          Length[nativeBatch["Results"]],
        "ExportedResults" ->
          fixtureExportResults[nativeBatch["Results"]]|>]],
    ft2NativeReleaseBatch = Function[nativeBatch,
      AppendTo[checkpointEvents, "release"];
      <|"Released" -> 1, "Failures" -> {}|>],
    ft2NativeReleaseAtlas = Function[atlas,
      <|"Released" -> 1, "Failures" -> {}|>]},
  ft2RunNativeBoundaryDispatch[
    <|"Matrix" -> IdentityMatrix[2], "Variable" -> x|>,
    sourceRows, entries, ledger, x, 11/23, {x, 1 - x},
    {{x, 1}, {1 - x, 1}}, 6, 50, 8, 8,
    "synthetic-epsilon-plan",
    checkpointSpec]];
assert["native sidecar publisher runs synchronously after schema-2 save and before export",
  AssociationQ[checkpointed] &&
    checkpointEvents === {"run", "save", "publish", "export", "release"} &&
    ft2NativeTransportResumeRecordQ[publishedResume] &&
    ft2NativeCheckpointRecordQ[publishedAudit] &&
    checkpointed["NativeTransportCheckpoint"] === publishedResume,
  {checkpointEvents, publishedResume}];

restoreCounts = <|"segment" -> 0, "prepare" -> 0, "run" -> 0,
  "restore" -> 0, "export" -> 0|>;
restoredDispatch = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      restoreCounts["segment"]++; $Failed],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, requiredTargetMax, threads},
      restoreCounts["prepare"]++; $Failed],
    ft2NativeRun = Function[
      {atlas, observables, variable, matchingCertificationDigits,
       publicationDigits},
      restoreCounts["run"]++; $Failed],
    ft2NativeRestoreCheckpoint = Function[manifest,
      restoreCounts["restore"]++;
      <|"Type" -> "DiffExp2NativeTransportObservableBatch",
        "Atlas" -> None, "NativeMarches" -> 0,
        "RestoredNativeMarches" -> 2,
        "Results" -> capturedObservables|>],
    ft2NativeExport = Function[{nativeBatch, digits},
      restoreCounts["export"]++;
      Join[nativeBatch, <|"CompatibilityExports" ->
          Length[nativeBatch["Results"]],
        "ExportedResults" ->
          fixtureExportResults[nativeBatch["Results"]]|>]],
    ft2NativeReleaseBatch = Function[nativeBatch,
      <|"Released" -> 1, "Failures" -> {}|>],
    ft2NativeReleaseAtlas = Function[atlas,
      <|"Released" -> 1, "Failures" -> {}|>]},
  ft2RunNativeBoundaryDispatch[
    <|"Matrix" -> IdentityMatrix[2], "Variable" -> x|>,
    sourceRows, entries, ledger, x, 11/23, {x, 1 - x},
    {{x, 1}, {1 - x, 1}}, 6, 50, 8, 8,
    "synthetic-epsilon-plan",
    <|"Mode" -> "Restore", "Record" -> publishedResume,
      "ContractIdentity" -> checkpointContractIdentity|>]];
assert["resume dispatch restores and exports without replanning or remarching",
  AssociationQ[restoredDispatch] &&
    restoreCounts === <|"segment" -> 0, "prepare" -> 0, "run" -> 0,
      "restore" -> 1, "export" -> 1|> &&
    restoredDispatch["NativeBatchCalls"] === 0 &&
    restoredDispatch["NativeMarches"] === 0 &&
    TrueQ[restoredDispatch["RestoredNativeTransport"]] &&
    And @@ MapThread[DiffExp2`EpsSeries`ESSameQ,
      {restoredDispatch["Values"], checkpointed["Values"]}],
  {restoreCounts, restoredDispatch}];
assert["fresh and restored dispatch preserve identical certification records",
  restoredDispatch["Certifications"] === checkpointed["Certifications"] &&
    restoredDispatch["Certifications"] === mixed["Certifications"],
  {restoredDispatch["Certifications"], checkpointed["Certifications"]}];

directRequests = {
  request[1, "direct", 0, 0, {14}],
  request[2, "integrate", 1, 1, {15}]};
directBatch = <|"Schema" -> "FeynmanTrick.LevelIBPBatch/v1",
  "Key" -> StringRepeat["c", 64],
  "PayloadKey" -> StringRepeat["d", 64],
  "KeyRecord" -> {"direct-only"}, "UpperLevel" -> 3,
  "MastersAbove" -> {{1, 0}, {0, 1}},
  "BoundaryRequests" -> directRequests,
  "CoefficientVectors" -> Association[{
    {14} -> {epsilon, 2}, {15} -> {0, 0}}]|>;
directEntries = ft2PrepareBoundaryEntries[
  3, directBatch, {1, 0}, x, epsilon, normalizeIdentity];
directLedger = ft2NativeEpsilonLedger[directEntries, sourceRows, 4, 4];
counts = newCounts[];
directOnly = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; $Failed],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, requiredTargetMax, threads},
      counts["prepare"]++; $Failed],
    ft2NativeRun = Function[
      {atlas, observables, variable, matchingCertificationDigits,
       publicationDigits},
      counts["run"]++; $Failed],
    ft2NativeExport = Function[{nativeBatch, digits},
      counts["export"]++; $Failed],
    ft2NativeReleaseBatch = Function[nativeBatch,
      counts["releaseBatch"]++; $Failed],
    ft2NativeReleaseAtlas = Function[atlas,
      counts["releaseAtlas"]++; $Failed]},
  ft2RunNativeBoundaryDispatch[
    <|"Matrix" -> IdentityMatrix[2], "Variable" -> x|>,
    sourceRows, directEntries, directLedger, x, 11/23, {},
    {{x, 1}}, 6, 50, 8, 8]];
assert["direct-only and proven-zero work skips atlas, march, export, and release",
  AssociationQ[directOnly] && Total[Values[counts]] === 0 &&
    directOnly["NativeBatchCalls"] === 0 &&
    directOnly["NativeMarches"] === 0 &&
    directOnly["CompatibilityExports"] === 0 &&
    DiffExp2`EpsSeries`ESCoefficient[
      directOnly["Values"][[1]], 0] === 5 &&
    directLedger["IntegrationHalo"] === 0 &&
    directOnly["Certifications"] === {
      ft2NotApplicableCertification["direct", "direct"],
      ft2NotApplicableCertification["integrate", "proven-zero"]} &&
    ft2NativeCheckpointRecordQ[directOnly["CheckpointRecord"]],
  {counts, directOnly}];

counts = newCounts[];
malformedBatch = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; <|"Path" -> path|>],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, requiredTargetMax, threads},
      counts["prepare"]++;
      <|"Type" -> "DiffExp2NativeRegularIndependentArmAtlas",
        "PlanCheckpointIdentity" -> "malformed-run-atlas"|>],
    ft2NativeRun = Function[
      {atlas, observables, variable, matchingCertificationDigits,
       publicationDigits},
      counts["run"]++;
      <|"Type" -> "MalformedPublishedBatch", "Atlas" -> atlas|>],
    ft2NativeExport = Function[{nativeBatch, digits},
      counts["export"]++; $Failed],
    ft2NativeReleaseBatch = Function[nativeBatch,
      counts["releaseBatch"]++; $Failed],
    ft2NativeReleaseAtlas = Function[atlas,
      counts["releaseAtlas"]++;
      <|"Released" -> 1, "Failures" -> {}|>]},
  ft2RunNativeBoundaryDispatch[
    <|"Matrix" -> IdentityMatrix[2], "Variable" -> x|>,
    sourceRows, entries, ledger, x, 11/23, {x, 1 - x},
    {{x, 1}, {1 - x, 1}}, 6, 50, 8, 8]];
assert["malformed published batch falls back to releasing its atlas owner",
  FailureQ[malformedBatch] &&
    counts === <|"segment" -> 2, "prepare" -> 1, "run" -> 1,
      "export" -> 0, "releaseBatch" -> 0, "releaseAtlas" -> 1|>,
  {counts, malformedBatch}];

printedRows = printRows["certification-fixture", 0, {{1, 1}},
  {mixed["Values"][[1]]}, {0}, {mixed["Certifications"][[1]]}];
certifiedStepRow = ft2StepwiseRow["certification-fixture", 0, {1, 1},
  mixed["Values"][[1]], 0, mixed["Certifications"][[1]]];
storedStepRow = ft2StepwiseRow["certification-fixture", 0, {2, 0},
  mixed["Values"][[6]], 0, mixed["Certifications"][[6]]];
directStepRow = ft2StepwiseRow["certification-fixture", 0, {0, 0},
  mixed["Values"][[4]], 0, mixed["Certifications"][[4]]];
certifiedFinalRow = ft2FinalRow["certification-fixture",
  mixed["Values"][[1]], mixed["Certifications"][[1]]];
syntheticOutput = StringRiffle[{
  ft2OutputLine["STEPWISE ", certifiedStepRow],
  ft2OutputLine["STEPWISE ", storedStepRow],
  ft2OutputLine["STEPWISE ", directStepRow],
  ft2OutputLine["FINAL ", certifiedFinalRow]}, "\n"];
parsedPipeline = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  <|"Schema" -> "definitions-only-plan"|>,
  <|"StandardOutput" -> syntheticOutput, "StandardError" -> "",
    "ExitCode" -> 0|>];
positiveRaw = DiffExp2`EpsSeries`ESNew[-1,
  {11, 22, 3 + 4 I, 5 - 6 I}];
negligiblePoleRaw = DiffExp2`EpsSeries`ESNew[-2,
  {10^-900, 0, 39.65552683429765`20}];
catastrophicPoleRaw = DiffExp2`EpsSeries`ESNew[-2,
  {10^17, 10^41, 10^67}];
doubleBoxRemnantRaw = DiffExp2`EpsSeries`ESNew[-5,
  {0``13, 12, 2.626, -77.8, -201.7, -247.0}];
doubleBoxBadPoleRaw = DiffExp2`EpsSeries`ESNew[-5,
  {1, 12, 2.626, -77.8, -201.7, -247.0}];
positiveCertification =
  ft2NotApplicableCertification["direct", "direct"];
legacyStepRow = Block[{epsOrder = 0},
  ft2StepwiseRow["positive-output-fixture", 0, {1}, positiveRaw, 0,
    positiveCertification]];
legacyFinalRow = Block[{epsOrder = 0},
  ft2FinalRow["positive-output-fixture", positiveRaw,
    positiveCertification]];
positiveStepRow = Block[{epsOrder = 2},
  ft2StepwiseRow["positive-output-fixture", 0, {1}, positiveRaw, 0,
    positiveCertification]];
positiveFinalRow = Block[{epsOrder = 2},
  ft2FinalRow["positive-output-fixture", positiveRaw,
    positiveCertification]];
incompleteRaw = DiffExp2`EpsSeries`ESNew[0, {1, 2}];
incompleteStepRow = Block[{epsOrder = 2},
  ft2StepwiseRow["incomplete-output-fixture", 0, {1}, incompleteRaw, 0,
    positiveCertification]];
incompleteFinalRow = Block[{epsOrder = 2},
  ft2FinalRow["incomplete-output-fixture", incompleteRaw,
    positiveCertification]];
negligiblePoleAudit = ft2PretrimFinalPoleAudit[
  "banana4", {{1}}, {negligiblePoleRaw}];
catastrophicPoleAudit = ft2PretrimFinalPoleAudit[
  "banana4", {{1}}, {catastrophicPoleRaw}];
doubleBoxRemnantAudit = ft2PretrimFinalPoleAudit[
  "double_box_planar", {ConstantArray[1, 7]},
  {doubleBoxRemnantRaw}];
doubleBoxBadPoleAudit = ft2PretrimFinalPoleAudit[
  "double_box_planar", {ConstantArray[1, 7]},
  {doubleBoxBadPoleRaw}];
doubleBoxFloored = ft2ApplyExpectedFinalPoleFloor[
  "double_box_planar", doubleBoxRemnantRaw];
positiveSyntheticOutput = StringRiffle[{
  ft2OutputLine["STEPWISE ", positiveStepRow],
  ft2OutputLine["FINAL ", positiveFinalRow]}, "\n"];
positiveParsedPipeline =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    <|"Schema" -> "positive-definitions-only-plan"|>,
    <|"StandardOutput" -> positiveSyntheticOutput,
      "StandardError" -> "", "ExitCode" -> 0|>];
complexCoefficientQ[value_, re_, im_] := AssociationQ[value] &&
  TrueQ[Abs[N[value["Re"] - re, 30]] < 10^-20] &&
  TrueQ[Abs[N[value["Im"] - im, 30]] < 10^-20];
assert["STEPWISE printer retains one compact certification per master",
  ListQ[printedRows] && Length[printedRows] === 1 &&
    printedRows[[1, "Certification"]] === mixed["Certifications"][[1]] &&
    StringStartsQ[ft2OutputLine["STEPWISE ", printedRows[[1]]],
      "STEPWISE {"] &&
    FailureQ[printRows["bad-length", 0, {{1}},
      {mixed["Values"][[1]]}, {0}, {}]],
  printedRows];
assert["facade parser preserves certified stored and not-applicable records",
  AssociationQ[parsedPipeline] && parsedPipeline["Status"] === "Succeeded" &&
    Length[parsedPipeline["Stepwise"]] === 3 &&
    parsedPipeline["Stepwise"][[1, "Certification", "Scope"]] ===
      "full_local_with_certified_tail" &&
    parsedPipeline["Stepwise"][[1, "Certification", "ErrorEnvelope",
      "guarantee"]] === "certified" &&
    parsedPipeline["Stepwise"][[2, "Certification"]] === <|
      "Applicability" -> "applicable", "Operation" -> "integrate",
      "Scope" -> "stored_truncation", "ErrorGuarantee" -> "none",
      "ErrorEnvelope" -> Null|> &&
    parsedPipeline["Stepwise"][[3, "Certification"]] ===
      ft2NotApplicableCertification["direct", "direct"] &&
    parsedPipeline["Final", "Certification"] ===
      parsedPipeline["Stepwise"][[1, "Certification"]],
  parsedPipeline];
assert["order-zero output retains the exact singleton-compatible shape",
  AssociationQ[legacyStepRow] && AssociationQ[legacyFinalRow] &&
    legacyStepRow["Coefficients"][[All, 1]] === {-1, 0} &&
    Keys[legacyFinalRow] ===
      {"Example", "Finite", "RawMinPower", "Certification"} &&
    !KeyExistsQ[legacyFinalRow, "Coefficients"] &&
    TrueQ[Abs[N[legacyFinalRow["Finite"] - 22, 30]] < 10^-20],
  {legacyStepRow, legacyFinalRow}];
assert["positive epsilon requests expose every STEPWISE and FINAL coefficient",
  AssociationQ[positiveStepRow] && AssociationQ[positiveFinalRow] &&
    positiveStepRow["Coefficients"][[All, 1]] === {-1, 0, 1, 2} &&
    positiveFinalRow["Coefficients"][[All, 1]] === {-1, 0, 1, 2} &&
    TrueQ[Abs[N[positiveFinalRow["Finite"] - 22, 30]] < 10^-20] &&
    complexCoefficientQ[
      positiveStepRow["Coefficients"][[3, 2]], 3, 4] &&
    complexCoefficientQ[
      positiveFinalRow["Coefficients"][[4, 2]], 5, -6],
  {positiveStepRow, positiveFinalRow}];
assert["output completeness fails explicitly below the requested epsilon top",
  FailureQ[incompleteStepRow] && FailureQ[incompleteFinalRow] &&
    incompleteStepRow["RequestedCompleteMax"] === 2 &&
    incompleteStepRow["AvailableCompleteMax"] === 1 &&
    incompleteFinalRow["RequestedCompleteMax"] === 2 &&
    incompleteFinalRow["AvailableCompleteMax"] === 1,
  {incompleteStepRow, incompleteFinalRow}];
assert["pre-trim pole audit accepts harmless absolute remnants but rejects catastrophic poles",
  TrueQ[negligiblePoleAudit] &&
    FailureQ[catastrophicPoleAudit] &&
    catastrophicPoleAudit[[1]] ===
      "FeynmanTrickUnexpectedFinalPole",
  {negligiblePoleAudit, catastrophicPoleAudit}];
assert["proven double-box pole floor removes only audited sub-floor remnants",
  TrueQ[doubleBoxRemnantAudit] &&
    FailureQ[doubleBoxBadPoleAudit] &&
    doubleBoxBadPoleAudit[[1]] ===
      "FeynmanTrickUnexpectedFinalPole" &&
    DiffExp2`EpsSeries`ESMinPower[doubleBoxFloored] === -4 &&
    DiffExp2`EpsSeries`ESCoefficient[doubleBoxFloored, -4] === 12,
  {doubleBoxRemnantAudit, doubleBoxBadPoleAudit, doubleBoxFloored}];
assert["facade parser preserves positive-order Laurent rows and complex encoding",
  AssociationQ[positiveParsedPipeline] &&
    positiveParsedPipeline["Status"] === "Succeeded" &&
    positiveParsedPipeline["Stepwise"][[1, "Coefficients"]][[All, 1]] ===
      {-1, 0, 1, 2} &&
    positiveParsedPipeline["Final", "Coefficients"][[All, 1]] ===
      {-1, 0, 1, 2} &&
    complexCoefficientQ[
      positiveParsedPipeline["Final", "Coefficients"][[3, 2]], 3, 4],
  positiveParsedPipeline];

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
