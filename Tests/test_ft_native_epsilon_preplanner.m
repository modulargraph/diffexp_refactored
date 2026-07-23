(* Definitions-only proof fixture for the exact native full-ladder epsilon
   planner.  No FIRE process, chart solve, or persistent C++ session runs. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> "Cpp"];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> "1"];

Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label,
    If[detail === None, "", ": "], If[detail === None, "", detail]]];

eps = Global`eps;
x = Global`xNativePlan;
request[index_, case_, vi_, vj_, needed_] := <|
  "MasterIndex" -> index, "MasterVec" -> {vi, vj},
  "Case" -> case, "Vi" -> vi, "Vj" -> vj,
  "NeededVec" -> needed|>;

requests = {
  request[1, "integrate", 1, 1, {1}],
  request[2, "direct", 0, 0, {2}],
  request[3, "integrate", 1, 1, {3}],
  request[4, "limitLower", 0, 1, {4}]};
batch = <|
  "Schema" -> "FeynmanTrick.LevelIBPBatch/v1",
  "Key" -> StringRepeat["a", 64],
  "PayloadKey" -> StringRepeat["b", 64],
  "KeyRecord" -> {"native-epsilon-plan-fixture"},
  "UpperLevel" -> 1,
  "MastersAbove" -> {{1, 0}, {0, 1}},
  "BoundaryRequests" -> requests,
  "CoefficientVectors" -> Association[{
    {1} -> {eps^-5, 0},
    {2} -> {0, eps^-2},
    {3} -> {0, 0},
    {4} -> {eps^2, 0}}]|>;
ftData = <|"NumLevels" -> 1, "Levels" -> <|
  0 -> <|"Masters" -> {{1, 1}, {0, 0}, {1, 1}, {0, 1}}|>,
  1 -> <|"Masters" -> {{1, 0}, {0, 1}},
    "FeynmanParameter" -> x,
    "CombinedPositions" -> {1, 2},
    "DiffMatrix" -> {{0, 1/(eps*x)}, {0, 0}}|>|>|>;

plannerMatrixNormalizeCalls = 0;
plan = Block[{FeynmanTrick`LevelReduction`PrepareLevelIBPBatch =
      Function[Null, $Failed]},
  ft2BuildNativeEpsilonPlan[ftData, 0, {0},
    Function[value,
      If[value === ftData["Levels", 1, "DiffMatrix"],
        plannerMatrixNormalizeCalls++];
      value],
    <|1 -> batch|>]];
levelPlan = If[AssociationQ[plan], plan["Levels"][1], <||>];

assert["planner normalizes the exact relative diagonal gauge",
  AssociationQ[plan] &&
    plannerMatrixNormalizeCalls === 1 &&
    levelPlan["Gauge", "RelativePrefactors"] === {1, 0} &&
    levelPlan["Gauge", "Record", "PoleFree"] === True,
  {plannerMatrixNormalizeCalls, plan}];

runtimeNormalizeCalls = 0;
runtimeMatrix = ft2RuntimeLevelMatrix[
  ftData["Levels", 1], levelPlan,
  Function[value, runtimeNormalizeCalls++; value]];
assert["native runtime reuses the planner-owned normalized level matrix",
  plannerMatrixNormalizeCalls === 1 &&
    runtimeNormalizeCalls === 0 &&
    runtimeMatrix === levelPlan["Gauge", "Matrix"],
  {plannerMatrixNormalizeCalls, runtimeNormalizeCalls, runtimeMatrix}];

fallbackNormalizeCalls = 0;
fallbackMatrix = ft2RuntimeLevelMatrix[
  ftData["Levels", 1], None,
  Function[value, fallbackNormalizeCalls++; value]];
assert["non-native runtime retains the direct normalization path",
  fallbackNormalizeCalls === 1 &&
    fallbackMatrix === ftData["Levels", 1, "DiffMatrix"],
  {fallbackNormalizeCalls, fallbackMatrix}];

tamperedRuntimeLevel = Association[levelPlan];
tamperedRuntimeGauge = Association[levelPlan["Gauge"]];
tamperedRuntimeGauge["Matrix"] = IdentityMatrix[2];
tamperedRuntimeLevel["Gauge"] = tamperedRuntimeGauge;
assert["runtime rejects a planner matrix whose exact gauge hash changed",
  FailureQ[ft2RuntimeLevelMatrix[
    ftData["Levels", 1], tamperedRuntimeLevel,
    Function[value, value]]],
  tamperedRuntimeLevel];

tamperedRuntimeData = Association[ftData["Levels", 1]];
tamperedRuntimeData["DiffMatrix"] = IdentityMatrix[2];
assert["runtime rejects same-size raw input changed after planning",
  FailureQ[ft2RuntimeLevelMatrix[
    tamperedRuntimeData, levelPlan, Function[value, value]]],
  tamperedRuntimeData];

runnerSource = Import[
  FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}], "Text"];
assert["production level loop consumes the checked planned matrix path",
  StringCount[runnerSource,
      "A = ft2RuntimeLevelMatrix[levelData, plannedLevel, normalizeFT]"] ===
      1 &&
    !StringContainsQ[runnerSource,
      "A = normalizeFT[levelData[\"DiffMatrix\"]]"],
  "the production call site must not bypass ft2RuntimeLevelMatrix"];
assert["mixed direct/integrate losses include the one native primitive row",
  levelPlan["Record", "RelativeMinimumEpsilonShifts"] ===
      {-6, -2, 0, 1} &&
    levelPlan["Record", "EntryLosses"] === <|1 -> 7, 2 -> 2, 4 -> 0|> &&
    levelPlan["Record", "IntrinsicLoss"] === 7,
  levelPlan];
assert["proven-zero rows do not enter the level loss",
  levelPlan["Record", "ProvenZero"] === {False, False, True, False} &&
    !KeyExistsQ[levelPlan["Record", "EntryLosses"], 3],
  levelPlan];
assert["automatic raw requirement covers a loss greater than four",
  plan["Record", "DeepRequiredRawTop"] === 7 &&
    levelPlan["RequiredOutputRawTop"] === 0 &&
    levelPlan["RequiredRawTop"] === 7,
  plan["Record"]];

runtimeEntries = ft2PrepareBoundaryEntries[
  1, batch, {4, 3}, x, eps, Function[value, value]];
runtimeParity = ft2ValidateNativePlanRuntimeLevel[
  levelPlan, {4, 3}, runtimeEntries, eps];
sourceRows = Table[Join[{i}, ConstantArray[0, 10]], {i, 1, 2}];
runtimeLedger = ft2NativeEpsilonLedger[runtimeEntries, sourceRows, 0, 0];
assert["common prefactor shifts cancel instead of becoming recursive halos",
  AssociationQ[runtimeParity] && runtimeParity["CommonOffset"] === 3 &&
    Lookup[runtimeEntries, "MinimumEpsilonShift"] === {-9, -5, 0, -2} &&
    AssociationQ[runtimeLedger] &&
    runtimeLedger["DeliverableCompleteMax"] === 0,
  {runtimeParity, runtimeLedger}];

deepBoundary[requested_Integer, prefactor_Integer] := Module[
  {working = requested + prefactor},
  <|"BoundaryValues" -> ConstantArray[
      ConstantArray[0, working + 1], 2],
    "EpsPrefactors" -> ConstantArray[prefactor, 2],
    "WorkingEpsilonOrder" -> working,
    "RequestedEpsilonOrder" -> requested|>];
tooShallow = ft2FinalizeNativeEpsilonPlan[
  plan, deepBoundary[6, 3], 6];
execution = ft2FinalizeNativeEpsilonPlan[
  plan, deepBoundary[7, 3], 7];
assert["deepest API padding covers the gauge offset without charging it twice",
  FailureQ[tooShallow] && AssociationQ[execution] &&
    execution["DeepGaugeOffset"] === 3 &&
    execution["DeepRequiredSourceCompleteMax"] === 10 &&
    execution["DeepRequestedBoundaryOrder"] === 7 &&
    execution["DeepBoundaryCompleteMax"] === 10 &&
    execution["Record", "DeepBoundaryWorkingEpsilonOrder"] === 10 &&
    execution["Record", "DeepRequestedBoundarySurplus"] === 0 &&
    execution["Record", "DeepSourceSurplus"] === 0 &&
    ft2NativeEpsilonExecutionRecordQ[
      execution["Record"], execution["Identity"], plan],
  {tooShallow, execution}];
tamperedExecution = execution["Record"];
tamperedExecution["DeepGaugeOffset"] = 4;
assert["execution-plan checkpoint tampering is rejected",
  !ft2NativeEpsilonExecutionRecordQ[
    tamperedExecution, execution["Identity"], plan]];

tamperedPlan = Association[plan];
tamperedLevels = Association[plan["Levels"]];
tamperedLevel = Association[plan["Levels"][1]];
tamperedLevel["RequiredRawTop"] = 8;
AssociateTo[tamperedLevels, 1 -> tamperedLevel];
AssociateTo[tamperedPlan, "Levels" -> tamperedLevels];
assert["runtime plan data must exactly match its checkpoint record",
  !ft2NativeEpsilonPlanQ[tamperedPlan]];

plusPrescriptions = Block[{deltaPrescriptionSign = 1},
  levelDeltaPrescriptions[x, <|"SingularFactors" -> {}|>, {}]];
minusPrescriptions = Block[{deltaPrescriptionSign = -1},
  levelDeltaPrescriptions[x, <|"SingularFactors" -> {}|>, {}]];
assert["+1 and -1 delta rims propagate into distinct checkpoint identities",
  AllTrue[plusPrescriptions, Last[#] === 1 &] &&
  AllTrue[minusPrescriptions, Last[#] === -1 &] &&
  ft2CanonicalIdentity["ft2-delta-prescriptions-", plusPrescriptions] =!=
      ft2CanonicalIdentity["ft2-delta-prescriptions-", minusPrescriptions],
  {plusPrescriptions, minusPrescriptions}];

profileTmp = CreateDirectory[FileNameJoin[{$TemporaryDirectory,
  "DiffExp2_matching_halo_profile_" <> ToString[$ProcessID]}]];
profileContract = ft2MatchingHaloProfileContract[
  "matching-profile-fixture", "prepared-fixture", plan];
profileFile = Block[{matchingHaloProfileRoot = profileTmp},
  ft2MatchingHaloProfileFile[profileContract]];
profileSaved = ft2SaveMatchingHaloProfile[
  profileFile, profileContract, <|1 -> 2|>];
profileLoaded = ft2LoadMatchingHaloProfile[
  profileFile, profileContract];
profileSavedLower = ft2SaveMatchingHaloProfile[
  profileFile, profileContract, <|1 -> 1|>];
profileLoadedAfterLower = ft2LoadMatchingHaloProfile[
  profileFile, profileContract];
runLevelMergedBounds = ft2MergeMatchingHaloBounds[1,
  {profileLoadedAfterLower, <|1 -> 1|>}];
assert["learned matching-halo profiles round-trip under an exact contract",
  ft2MatchingHaloProfileContractQ[profileContract] &&
    AssociationQ[profileSaved] && profileLoaded === {2} &&
    AssociationQ[profileSavedLower] && profileLoadedAfterLower === {2},
  {profileContract, profileSaved, profileLoaded,
    profileSavedLower, profileLoadedAfterLower}];
assert["run-level merge keeps a loaded higher halo over an explicit lower bound",
  runLevelMergedBounds === <|1 -> 2|>, runLevelMergedBounds];
assert["profile writes release the same-filesystem serialization lock",
  !DirectoryQ[profileFile <> ".lock"], profileFile <> ".lock"];

mergedProfileBounds = ft2MergeMatchingHaloBounds[3,
  {<|2 -> 5|>, <|1 -> 1, 2 -> 2|>}];
assert["explicit direct-call halos and learned profiles merge only upward",
  mergedProfileBounds === <|1 -> 1, 2 -> 5, 3 -> 0|> &&
    FailureQ[ft2MergeMatchingHaloBounds[3,
      {<|1 -> -1|>, <||>}]] &&
    FailureQ[ft2MergeMatchingHaloBounds[3,
      {<|1 -> $ft2MatchingHaloProfileMax + 1|>, <||>}]],
  mergedProfileBounds];

changedProfileContract = Block[{matchDigits = matchDigits + 1},
  ft2MatchingHaloProfileContract[
    "matching-profile-fixture", "prepared-fixture", plan]];
changedProfileFile = Block[{matchingHaloProfileRoot = profileTmp},
  ft2MatchingHaloProfileFile[changedProfileContract]];
assert["matching-affecting configuration changes select a distinct profile",
  ft2MatchingHaloProfileContractQ[changedProfileContract] &&
    changedProfileContract["Identity"] =!= profileContract["Identity"] &&
    changedProfileFile =!= profileFile &&
    ft2LoadMatchingHaloProfile[profileFile, changedProfileContract] === {0},
  {profileContract, changedProfileContract}];

sourceIndependentProfileContract =
  Block[{$ftLadderSourceFingerprint = "unrelated-source-edit"},
    ft2MatchingHaloProfileContract[
      "matching-profile-fixture", "prepared-fixture", plan]];
assert[
  "unrelated source edits do not discard a certified matching-halo lower bound",
  sourceIndependentProfileContract === profileContract,
  {profileContract, sourceIndependentProfileContract}];

noncanonicalBasePlan = Join[profileContract["Record", "BasePlanRecord"],
  <|"MatchingPrivateHalos" -> {1}|>];
noncanonicalRecord = Join[profileContract["Record"], <|
  "BasePlanRecord" -> noncanonicalBasePlan,
  "BasePlanIdentity" -> ft2CanonicalIdentity[
    "ft2-native-epsilon-plan-", noncanonicalBasePlan]|>];
noncanonicalContract = <|"Record" -> noncanonicalRecord,
  "Identity" -> ft2CanonicalIdentity[
    "ft2-matching-halo-profile-contract-", noncanonicalRecord],
  "NumLevels" -> 1|>;
assert["profile contracts reject a nonzero or noncanonical base plan",
  !ft2MatchingHaloProfileContractQ[noncanonicalContract],
  noncanonicalContract];

Clear[Global`$FT2MatchingHaloProfile]; Get[profileFile];
Global`$FT2MatchingHaloProfile = Join[Global`$FT2MatchingHaloProfile,
  <|"MatchingPrivateHalos" -> {3}|>];
DumpSave[profileFile, Global`$FT2MatchingHaloProfile];
Clear[Global`$FT2MatchingHaloProfile];
tamperedProfileLoad = ft2LoadMatchingHaloProfile[
  profileFile, profileContract];
assert["tampered profile bounds are rejected instead of becoming trusted",
  tamperedProfileLoad === {0}, tamperedProfileLoad];
assert["profile publication is atomic and leaves no temporary snapshots",
  FileNames["*.tmp-*.mx", profileTmp] === {},
  FileNames["*", profileTmp]];
assert["the learned-profile cache has an explicit environment opt-out",
  StringContainsQ[runnerSource,
    "FT_DISABLE_MATCHING_HALO_PROFILE"]];
Quiet[DeleteDirectory[profileTmp, DeleteContents -> True]];

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> None];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
