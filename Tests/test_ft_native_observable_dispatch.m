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
    ft2NativeRun = Function[{atlas, observables, variable},
      counts["run"]++; capturedObservables = observables;
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
    {{x, 1}, {1 - x, 1}}, 6, 50]];

assert["mixed dispatch plans both arms but invokes one native batch/export",
  AssociationQ[mixed] &&
    counts === <|"segment" -> 2, "prepare" -> 1, "run" -> 1,
      "export" -> 1, "releaseBatch" -> 1, "releaseAtlas" -> 0|> &&
    mixed["NativeBatchCalls"] === 1 &&
    mixed["NativeMarches"] === 2 &&
    mixed["CompatibilityExports"] === 4,
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
    ft2NativeRun = Function[{atlas, observables, variable},
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
    {{x, 1}, {1 - x, 1}}, 6, 50, "synthetic-epsilon-plan",
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
    ft2NativeRun = Function[{atlas, observables, variable},
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
    {{x, 1}, {1 - x, 1}}, 6, 50, "synthetic-epsilon-plan",
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
    ft2NativeRun = Function[{atlas, observables, variable},
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
    {{x, 1}}, 6, 50]];
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
    ft2NativeRun = Function[{atlas, observables, variable},
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
    {{x, 1}, {1 - x, 1}}, 6, 50]];
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
