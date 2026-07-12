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
ledger = ft2NativeEpsilonLedger[entries, sourceRows, 7];
assert["epsilon ledger separates coefficient and integration halos",
  AssociationQ[ledger] &&
    KeyTake[ledger, {"AvailableSourceCompleteMax", "SourceCompleteMax",
      "CoefficientHalo",
      "IntegrationHalo", "TargetCompleteMax",
      "DeliverableCompleteMax",
      "DownstreamRawTop"}] ===
      <|"AvailableSourceCompleteMax" -> 10,
        "SourceCompleteMax" -> 10, "CoefficientHalo" -> 1,
        "IntegrationHalo" -> 1, "TargetCompleteMax" -> 9,
        "DeliverableCompleteMax" -> 8, "DownstreamRawTop" -> 7|>,
  ledger];
assert["observable minima do not discount independently required raw depth",
  ledger["OutputMinimums"] ===
    <|1 -> -2, 2 -> 1, 3 -> -1, 4 -> 0, 6 -> 19|> &&
    !KeyExistsQ[ledger["OutputMinimums"], 5],
  ledger["OutputMinimums"]];
terminalLedger = ft2NativeEpsilonLedger[entries, sourceRows, 0];
assert["terminal raw completeness is exactly epsOrder, independent of poles",
  AssociationQ[terminalLedger] &&
    terminalLedger["SourceCompleteMax"] === 10 &&
    terminalLedger["TargetCompleteMax"] === 9 &&
    terminalLedger["DeliverableCompleteMax"] === 8 &&
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
  StringContainsQ[runnerSource, "needTop = Min[kmaxAvail]"] &&
    StringContainsQ[runnerSource,
      "\"PreservedRawCompleteMax\" -> If[recurrenceBackend === \"Cpp\""] &&
    StringContainsQ[runnerSource,
      "\"PreservedSourceCompleteMax\" -> If["]];

newCounts[] := <|"segment" -> 0, "prepare" -> 0, "run" -> 0,
  "export" -> 0, "releaseBatch" -> 0, "releaseAtlas" -> 0|>;
counts = newCounts[];
capturedPrepare = None;
capturedObservables = None;

mixed = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; <|"Path" -> path|>],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, threads},
      counts["prepare"]++;
      capturedPrepare = <|"Boundary" -> boundary,
        "CoefficientVectors" -> coefficientVectors,
        "TargetCompleteMax" -> targetMax, "Threads" -> threads|>;
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
        "ExportedResults" -> Map[
          Append[#, "Value" -> DiffExp2`EpsSeries`ESZero[
            #["Epsilon", "RequiredCompleteMax"]]] &,
          nativeBatch["Results"]]|>]],
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
assert["atlas preparation receives every active non-direct vector and exact target",
  AssociationQ[capturedPrepare] &&
    capturedPrepare["TargetCompleteMax"] === 9 &&
    Length[capturedPrepare["CoefficientVectors"]] === 4 &&
    AllTrue[capturedPrepare["Boundary"],
      DiffExp2`EpsSeries`ESMinPower[#] === -1 &&
        DiffExp2`EpsSeries`ESCompleteMax[#] === 10 &],
  capturedPrepare];
assert["observable order and operation-specific epsilon windows are stable",
  Lookup[capturedObservables, "Operation"] ===
    {"integrate", "limitLower", "limitUpper", "integrate"} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Min"] ===
      {-2, 1, -1, 8} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Max"] ===
      {8, 8, 8, 8} &&
    Lookup[Lookup[capturedObservables, "Epsilon"],
      "RequiredCompleteMax"] === {8, 8, 8, 8} &&
    Keys[capturedObservables[[1]]] ===
      {"Operation", "Identity", "CheckpointIdentity",
       "CoefficientVector", "Epsilon", "TailPolicy"}];
assert["high-shift integration is marched rather than unsafely pruned",
  capturedObservables[[-1, "Identity"]] === entries[[6, "Identity"]] &&
    capturedObservables[[-1, "Epsilon", "Min"]] === 8];
assert["direct and proven-zero results merge in original master order",
  Length[mixed["Values"]] === 6 &&
    DiffExp2`EpsSeries`ESCoefficient[mixed["Values"][[4]], 0] === 5 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[4]]] === 8 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[5]]] === 8 &&
    DiffExp2`EpsSeries`ESMinPower[mixed["Values"][[1]]] === 8,
  DiffExp2`EpsSeries`ESWindow /@ mixed["Values"]];
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
       targetMax, threads},
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
        "ExportedResults" -> Map[
          Append[#, "Value" -> DiffExp2`EpsSeries`ESZero[
            #["Epsilon", "RequiredCompleteMax"]]] &,
          nativeBatch["Results"]]|>]],
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
       targetMax, threads}, restoreCounts["prepare"]++; $Failed],
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
        "ExportedResults" -> Map[
          Append[#, "Value" -> DiffExp2`EpsSeries`ESZero[
            #["Epsilon", "RequiredCompleteMax"]]] &,
          nativeBatch["Results"]]|>]],
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
directLedger = ft2NativeEpsilonLedger[directEntries, sourceRows, 4];
counts = newCounts[];
directOnly = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; $Failed],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, threads}, counts["prepare"]++; $Failed],
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
    ft2NativeCheckpointRecordQ[directOnly["CheckpointRecord"]],
  {counts, directOnly}];

counts = newCounts[];
malformedBatch = Block[{
    ft2NativeSegmentLine = Function[{system, path},
      counts["segment"]++; <|"Path" -> path|>],
    ft2NativePrepare = Function[
      {system, boundary, lower, upper, coefficientVectors, variable,
       targetMax, threads},
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

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
