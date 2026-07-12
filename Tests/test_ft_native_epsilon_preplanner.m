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

plan = Block[{FeynmanTrick`LevelReduction`PrepareLevelIBPBatch =
      Function[Null, $Failed]},
  ft2BuildNativeEpsilonPlan[ftData, 0, {0}, Function[value, value],
    <|1 -> batch|>]];
levelPlan = If[AssociationQ[plan], plan["Levels"][1], <||>];

assert["planner normalizes the exact relative diagonal gauge",
  AssociationQ[plan] &&
    levelPlan["Gauge", "RelativePrefactors"] === {1, 0} &&
    levelPlan["Gauge", "Record", "PoleFree"] === True,
  plan];
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
runtimeLedger = ft2NativeEpsilonLedger[runtimeEntries, sourceRows, 0];
assert["common prefactor shifts cancel instead of becoming recursive halos",
  AssociationQ[runtimeParity] && runtimeParity["CommonOffset"] === 3 &&
    Lookup[runtimeEntries, "MinimumEpsilonShift"] === {-9, -5, 0, -2} &&
    AssociationQ[runtimeLedger] &&
    runtimeLedger["DeliverableCompleteMax"] === 0,
  {runtimeParity, runtimeLedger}];

tooShallow = ft2FinalizeNativeEpsilonPlan[plan, {4, 3}, 9];
execution = ft2FinalizeNativeEpsilonPlan[plan, {4, 3}, 10];
assert["deepest prefactors are charged exactly once through their gauge offset",
  FailureQ[tooShallow] && AssociationQ[execution] &&
    execution["DeepGaugeOffset"] === 3 &&
    execution["DeepRequiredBoundaryOrder"] === 10 &&
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

SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> None];
SetEnvironment["DE2_RECURRENCE_BACKEND" -> None];
SetEnvironment["FT_DELTA_PRESCRIPTION_SIGN" -> None];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
