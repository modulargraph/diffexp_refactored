(* Focused production-runner test for the retained native observable cutover.
   Loads definitions only and replaces the six native seams, so no FIRE job,
   chart solve, or persistent C++ session is started. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> "Cpp"];
SetEnvironment["FT_INTEGRATION_POLE_ALLOWANCE" -> "2"];

Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label, If[detail === None, "", ": "],
    If[detail === None, "", detail]]];

x = Global`xNativeFT;
epsilon = Global`eps;
requests = {
  <|"Case" -> "integrate", "Vi" -> 2, "Vj" -> 3,
    "NeededVec" -> {1}|>,
  <|"Case" -> "limitLower", "Vi" -> 1, "Vj" -> 0,
    "NeededVec" -> {2}|>,
  <|"Case" -> "limitUpper", "Vi" -> 0, "Vj" -> 1,
    "NeededVec" -> {3}|>,
  <|"Case" -> "direct", "Vi" -> 0, "Vj" -> 0,
    "NeededVec" -> {4}|>,
  <|"Case" -> "integrate", "Vi" -> 1, "Vj" -> 1,
    "NeededVec" -> {5}|>};
vectors = Association[{
  {1} -> {2, 3/epsilon},
  {2} -> {epsilon^3, epsilon},
  {3} -> {1, epsilon^2},
  {4} -> {epsilon, 2},
  {5} -> {0, 0}}];
batch = <|"Key" -> "synthetic-batch",
  "PayloadKey" -> "synthetic-payload",
  "BoundaryRequests" -> requests,
  "CoefficientVectors" -> vectors|>;
normalizeIdentity = Function[value, value];
entries = ft2PrepareBoundaryEntries[
  3, batch, {1, 0}, x, epsilon, normalizeIdentity];

expectedWeight = x*(1 - x)^2;
assert["coefficient vectors include basis, beta, and x weights exactly",
  ListQ[entries] && Length[entries] === 5 &&
    And @@ MapThread[TrueQ[PossibleZeroQ[Together[#1 - #2]]] &,
      {entries[[1, "CoefficientVector"]],
       {24 expectedWeight/epsilon, 36 expectedWeight/epsilon}}] &&
    entries[[2, "CoefficientVector"]] === {epsilon^2, epsilon} &&
    entries[[4, "CoefficientVector"]] === {1, 2}];
assert["all-zero coefficient vectors are proved before halo accounting",
  TrueQ[entries[[5, "ProvenZero"]]] &&
    !TrueQ[entries[[1, "ProvenZero"]]]];

sourceRows = {
  Join[{1}, ConstantArray[0, 8]],
  Join[{2}, ConstantArray[0, 8]]};
ledger = ft2NativeEpsilonLedger[entries, sourceRows, 7];
assert["epsilon ledger separates coefficient and integration halos",
  AssociationQ[ledger] &&
    KeyTake[ledger, {"SourceCompleteMax", "CoefficientHalo",
      "IntegrationHalo", "TargetCompleteMax",
      "DeliverableCompleteMax", "PlannedBoundaryShift",
      "DownstreamRawTop"}] ===
      <|"SourceCompleteMax" -> 8, "CoefficientHalo" -> 1,
        "IntegrationHalo" -> 2, "TargetCompleteMax" -> 7,
        "DeliverableCompleteMax" -> 5,
        "PlannedBoundaryShift" -> 3, "DownstreamRawTop" -> 4|>,
  ledger];
assert["observable minima use shift-HI only for integration",
  ledger["OutputMinimums"] ===
    <|1 -> -3, 2 -> 1, 3 -> -1, 4 -> 0|> &&
    !KeyExistsQ[ledger["OutputMinimums"], 5],
  ledger["OutputMinimums"]];

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
      <|"Type" -> "SyntheticAtlas"|>],
    ft2NativeRun = Function[{atlas, observables, variable},
      counts["run"]++; capturedObservables = observables;
      <|"Type" -> "SyntheticBatch", "NativeMarches" -> 2,
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
    sourceRows, entries, ledger, x, 11/23, {x, 1 - x}, 6, 50]];

assert["mixed dispatch plans both arms but invokes one native batch/export",
  AssociationQ[mixed] &&
    counts === <|"segment" -> 2, "prepare" -> 1, "run" -> 1,
      "export" -> 1, "releaseBatch" -> 1, "releaseAtlas" -> 0|> &&
    mixed["NativeBatchCalls"] === 1 &&
    mixed["NativeMarches"] === 2 &&
    mixed["CompatibilityExports"] === 3,
  {counts, mixed}];
assert["atlas preparation receives every active non-direct vector and exact target",
  AssociationQ[capturedPrepare] &&
    capturedPrepare["TargetCompleteMax"] === 7 &&
    Length[capturedPrepare["CoefficientVectors"]] === 3 &&
    AllTrue[capturedPrepare["Boundary"],
      DiffExp2`EpsSeries`ESMinPower[#] === -2 &&
        DiffExp2`EpsSeries`ESCompleteMax[#] === 8 &],
  capturedPrepare];
assert["observable order and operation-specific epsilon windows are stable",
  Lookup[capturedObservables, "Operation"] ===
    {"integrate", "limitLower", "limitUpper"} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Min"] ===
      {-3, 1, -1} &&
    Lookup[Lookup[capturedObservables, "Epsilon"], "Max"] ===
      {5, 5, 5} &&
    Lookup[Lookup[capturedObservables, "Epsilon"],
      "RequiredCompleteMax"] === {5, 5, 5} &&
    Keys[capturedObservables[[1]]] ===
      {"Operation", "Identity", "CheckpointIdentity",
       "CoefficientVector", "Epsilon", "TailPolicy"}];
assert["direct and proven-zero results merge in original master order",
  Length[mixed["Values"]] === 5 &&
    DiffExp2`EpsSeries`ESCoefficient[mixed["Values"][[4]], 0] === 5 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[4]]] === 4 &&
    DiffExp2`EpsSeries`ESCompleteMax[mixed["Values"][[5]]] === 5,
  DiffExp2`EpsSeries`ESWindow /@ mixed["Values"]];

directEntries = entries[[{4, 5}]];
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
    sourceRows, directEntries, directLedger, x, 11/23, {}, 6, 50]];
assert["direct-only and proven-zero work skips atlas, march, export, and release",
  AssociationQ[directOnly] && Total[Values[counts]] === 0 &&
    directOnly["NativeBatchCalls"] === 0 &&
    directOnly["NativeMarches"] === 0 &&
    directOnly["CompatibilityExports"] === 0 &&
    DiffExp2`EpsSeries`ESCoefficient[
      directOnly["Values"][[1]], 0] === 5,
  {counts, directOnly}];

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];
SetEnvironment["FT_INTEGRATION_POLE_ALLOWANCE" -> None];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
