(* Process-free production-runner contract for selected and All custom families. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
environmentNames = {
  "FT_RUNNER_DEFINITIONS_ONLY", "FT_FAMILY_REQUEST_FILE",
  "FT_FAMILY_REQUEST_ID", "FT_EXAMPLES"
};
savedEnvironment = AssociationMap[Environment, environmentNames];
SetEnvironment[{
  "FT_RUNNER_DEFINITIONS_ONLY" -> "1",
  "FT_FAMILY_REQUEST_FILE" -> None,
  "FT_FAMILY_REQUEST_ID" -> None,
  "FT_EXAMPLES" -> None
}];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

requestDirectory = FileNameJoin[{$TemporaryDirectory,
  "ft-runner-family-contract-" <> ToString[$ProcessID]}];
If[DirectoryQ[requestDirectory],
  DeleteDirectory[requestDirectory, DeleteContents -> True]];

rawFamily = <|
  "Name" -> "runner_family_contract",
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {Global`p},
  "Propagators" -> {Global`l^2 + 1, (Global`l + Global`p)^2 + 2},
  "Replacements" -> {Global`p^2 -> Global`s},
  "NumericalPoint" -> {Global`s -> -1},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps,
  "CombinationSequence" -> {{1, 2}}
|>;
targets = {{1, 1}, {2, 0}};
family = FeynmanTrick`CreateFamily[rawFamily, targets];
request = FeynmanTrick`PipelineRequest`CreatePipelineRequest[family];
requestFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
  request, requestDirectory];
written = FeynmanTrick`PipelineRequest`WritePipelineRequest[
  request, requestFile];
loaded = ft2ResolveFamilyRequest[requestFile, request["RequestID"]];

assert["paired exact WXF request loads before execution",
  StringQ[written] && AssociationQ[loaded] && SameQ[loaded, request]];
assert["request file and identity are a strict pair",
  FailureQ[ft2ResolveFamilyRequest[requestFile, ""]] &&
    FailureQ[ft2ResolveFamilyRequest["", request["RequestID"]]] &&
    ft2ResolveFamilyRequest["", ""] === None];
assert["identity mismatch is rejected",
  FailureQ[ft2ResolveFamilyRequest[requestFile, "ft-request-wrong"]]];

allFamily = FeynmanTrick`CreateFamily[rawFamily, All];
allRequest = FeynmanTrick`PipelineRequest`CreatePipelineRequest[allFamily];
allFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
  allRequest, requestDirectory];
FeynmanTrick`PipelineRequest`WritePipelineRequest[allRequest, allFile];
loadedAll = ft2ResolveFamilyRequest[allFile, allRequest["RequestID"]];
assert["execution-ready All request enters production validation exactly",
  AssociationQ[loadedAll] && SameQ[loadedAll, allRequest] &&
    loadedAll["ExecutionPolicy", "MasterDiscovery"] === "AtExecution"];

symbolicMassFamily = FeynmanTrick`CreateFamily[
  Join[rawFamily, <|
    "Propagators" -> {
      Global`m1sq + Global`l^2,
      Global`m2sq + (Global`l + Global`p)^2},
    "Replacements" -> {Global`u -> Global`s},
    "NumericalPoint" -> {
      Global`u -> 99, Global`s -> Global`q, Global`q -> -1,
      Global`m1sq -> 3/2, Global`m2sq -> 5/4}
  |>], All];
symbolicMassRequest =
  FeynmanTrick`PipelineRequest`CreatePipelineRequest[symbolicMassFamily];
symbolicMassDiscoveryTopology =
  ft2AllDiscoveryTopology[symbolicMassRequest];
assert["All discovery freezes chained coefficients without rewriting replacement left-hand sides",
  FreeQ[symbolicMassDiscoveryTopology["Propagators"],
    Global`m1sq | Global`m2sq] &&
    symbolicMassDiscoveryTopology["Replacements"] ===
      {Global`u -> -1}];

allResolution = ft2ValidateAllDiscoveryMasters[
  allRequest, {{2, 0}, {1, 1}, {2, 0}}];
assert["runner canonicalizes discovered masters and binds their identities",
  AssociationQ[allResolution] &&
    allResolution["Masters"] === {{1, 1}, {2, 0}} &&
    Lookup[allResolution["OutputRequests"], "RequestOrdinal"] === {1, 2} &&
    AllTrue[allResolution["OutputRequests"],
      #["SelectionRequestID"] === allResolution["SelectionRequestID"] &&
      #["ResolutionID"] === allResolution["ResolutionID"] &]];
allResolutionLine = ft2OutputLine["OUTPUT_RESOLUTION ", allResolution];
assert["runner emits one machine-readable All discovery manifest",
  StringQ[allResolutionLine] &&
    StringStartsQ[allResolutionLine, "OUTPUT_RESOLUTION {"]];
setupWithoutAddedNumerators = Join[ft2AllDiscoveryTopology[allRequest], <|
  "OriginalPropagators" ->
    ft2AllDiscoveryTopology[allRequest]["Propagators"],
  "OriginalNumPropagators" -> 2,
  "NumPropagators" -> 2,
  "NumeratorPositions" -> {},
  "StartFileReady" -> True
|>];
setupWithAddedNumerator = Join[setupWithoutAddedNumerators, <|
  "NumPropagators" -> 3,
  "NumeratorPositions" -> {3}
|>];
assert["FIRE-added numerator slots fail before basis discovery",
  AssociationQ[ft2ValidateAllDiscoverySetup[
    allRequest, setupWithoutAddedNumerators]] &&
  FailureQ[ft2ValidateAllDiscoverySetup[
    allRequest, setupWithAddedNumerator]] &&
  FailureQ[ft2ValidateAllDiscoverySetup[allRequest,
    Join[setupWithoutAddedNumerators, <|"Replacements" -> {}|>]]]];
assert["merged-position numerator masters fail closed",
  FailureQ[ft2ValidateAllDiscoveryMasters[allRequest, {{-1, 1}}]] &&
    FailureQ[ft2ValidateAllDiscoveryMasters[allRequest, {{1, 1, 0}}]]];
setupCalls = 0; basisCalls = 0;
mockDiscovered = Block[{
    ft2AllSetupFIRE = Function[topology,
      setupCalls++;
      Join[topology, <|
        "OriginalNumPropagators" -> topology["NumPropagators"],
        "OriginalPropagators" -> topology["Propagators"],
        "NumeratorPositions" -> {}, "StartFileReady" -> True|>]],
    ft2AllFindBasis = Function[topology,
      basisCalls++;
      Join[topology, <|"Masters" -> {{2, 0}, {1, 1}, {2, 0}}|>]]},
  ft2DiscoverAllOutputSelection[allRequest]];
assert["execution-time discovery orchestrates one SetupFIRE and FindBasis seam",
  setupCalls === 1 && basisCalls === 1 && SameQ[mockDiscovered, allResolution]];

prescribedFamily = FeynmanTrick`CreateFamily[
  Join[rawFamily, <|"Prescriptions" -> {{Global`x, 1}}|>], targets];
prescribedRequest =
  FeynmanTrick`PipelineRequest`CreatePipelineRequest[prescribedFamily];
prescribedFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
  prescribedRequest, requestDirectory];
FeynmanTrick`PipelineRequest`WritePipelineRequest[
  prescribedRequest, prescribedFile];
assert["unwired analytic prescriptions fail instead of being ignored",
  FailureQ[ft2ResolveFamilyRequest[
    prescribedFile, prescribedRequest["RequestID"]]]];

assumedFamily = FeynmanTrick`CreateFamily[
  Join[rawFamily, <|
    "KinematicAssumptions" -> Global`s < 0|>], targets];
assumedRequest =
  FeynmanTrick`PipelineRequest`CreatePipelineRequest[assumedFamily];
assumedFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
  assumedRequest, requestDirectory];
FeynmanTrick`PipelineRequest`WritePipelineRequest[
  assumedRequest, assumedFile];
assert["unwired kinematic assumptions fail instead of being ignored",
  FailureQ[ft2ResolveFamilyRequest[
    assumedFile, assumedRequest["RequestID"]]]];

topologyOnlyPrescription = Join[family, <|
  "Topology" -> Join[family["Topology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>],
  "TopTopology" -> Join[family["TopTopology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>]
|>];
forgedPrescriptionCore =
  FeynmanTrick`PipelineRequest`Private`requestCore[
    topologyOnlyPrescription];
forgedPrescriptionRequest = Append[forgedPrescriptionCore,
  "RequestID" ->
    FeynmanTrick`PipelineRequest`Private`requestIdentity[
      forgedPrescriptionCore]];
forgedPrescriptionFile = FileNameJoin[
  {requestDirectory, "forged-topology-prescription.wxf"}];
Export[forgedPrescriptionFile, forgedPrescriptionRequest, "WXF"];
assert["production request load rejects topology-only prescriptions",
  FailureQ[ft2ResolveFamilyRequest[
    forgedPrescriptionFile, forgedPrescriptionRequest["RequestID"]]]];

threeFamily = FeynmanTrick`CreateFamily[<|
  "Name" -> "runner_three_sequence",
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {Global`p},
  "Propagators" -> {
    Global`l^2 + 1, (Global`l + Global`p)^2 + 2,
    (Global`l - Global`p)^2 + 3},
  "Replacements" -> {Global`p^2 -> -1},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps,
  "CombinationSequence" -> {{1, 2}, {1, 3}}
|>, All];
invalidSequenceFamily = Join[threeFamily, <|
  "CombinationSequence" -> {{1, 2}, {1, 2}}|>];
forgedSequenceCore =
  FeynmanTrick`PipelineRequest`Private`requestCore[
    invalidSequenceFamily];
forgedSequenceRequest = Append[forgedSequenceCore,
  "RequestID" ->
    FeynmanTrick`PipelineRequest`Private`requestIdentity[
      forgedSequenceCore]];
forgedSequenceFile = FileNameJoin[
  {requestDirectory, "forged-invalid-sequence.wxf"}];
Export[forgedSequenceFile, forgedSequenceRequest, "WXF"];
assert["production request load rejects repeated eliminated merge positions",
  FailureQ[ft2ResolveFamilyRequest[
    forgedSequenceFile, forgedSequenceRequest["RequestID"]]]];

runName = "family_" <> StringTake[request["RequestID"], -16];
FeynmanTrick`SetFTOption["DimensionExpression", family["Dimension"]];
contract = ftPrepCustomContractRecord[runName, request];
reverseFamily = FeynmanTrick`CreateFamily[rawFamily, Reverse[targets]];
reverseRequest =
  FeynmanTrick`PipelineRequest`CreatePipelineRequest[reverseFamily];
reverseContract = ftPrepCustomContractRecord[
  "family_" <> StringTake[reverseRequest["RequestID"], -16],
  reverseRequest];
allRunName = "family_" <> StringTake[allRequest["RequestID"], -16];
allContract = ftPrepCustomContractRecord[
  allRunName, allRequest, allResolution];

assert["custom preparation contract contains every mathematical input",
  contract["FamilyID"] === family["FamilyID"] &&
    contract["PipelineRequestID"] === request["RequestID"] &&
    contract["NumericalPoint"] === family["NumericalPoint"] &&
    contract["DimensionExpression"] === family["Dimension"] &&
    contract["OutputSelection"] === targets &&
    !AnyTrue[contract["CustomFamilySources"], FailureQ]];
assert["ordered L0 targets participate in preparation identity",
  ftPrepContractKey[contract] =!= ftPrepContractKey[reverseContract]];
assert["resolved All basis and original selection are both preparation-bound",
  AssociationQ[allContract] &&
    allContract["OutputSelectionMode"] === "AllResolved" &&
    allContract["OutputSelection"] === All &&
    allContract["ResolvedOutputSelection"] === allResolution["Masters"] &&
    allContract["ResolvedOutputRequests"] === allResolution["OutputRequests"] &&
    allContract["AllSelectionRequestID"] ===
      allResolution["SelectionRequestID"] &&
    allContract["OutputResolutionID"] === allResolution["ResolutionID"] &&
    ftPrepContractKey[allContract] =!= ftPrepContractKey[contract]];
assert["All preparation rejects a missing or stale resolution",
  FailureQ[ftPrepCustomContractRecord[allRunName, allRequest]] &&
    FailureQ[ftPrepCustomContractRecord[allRunName, allRequest,
      Join[allResolution, <|"ResolutionID" -> "stale"|>]]]];

ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  family["Topology"], family["CombinationSequence"],
  family["NumericalPoint"], "OutputIntegrals" -> targets];
ftData = Join[ftData, <|
  "FamilyID" -> family["FamilyID"],
  "PipelineRequestID" -> request["RequestID"]|>];
runtime =
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[];
setup = Join[runtime, <|
  "AutoDetectRestrictions" ->
    FeynmanTrick`Private`$FTConfig["AutoDetectRestrictions"]|>];
ftData["Levels"][1] = <|
  "Topology" -> <|"SetupFingerprintRecord" -> setup|>,
  "Masters" -> {{1}}, "DiffMatrix" -> {{0}}, "Computed" -> True|>;
assert["prepared custom data proves request, point, dimension, and L0 order",
  preparedFTDataMatchesContractQ[ftData, contract]];
assert["prepared-data request mismatch is rejected",
  !preparedFTDataMatchesContractQ[
    Join[ftData, <|"PipelineRequestID" -> "wrong"|>], contract]];

allFTData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  allFamily["Topology"], allFamily["CombinationSequence"],
  allFamily["NumericalPoint"],
  "OutputIntegrals" -> allResolution["Masters"]];
allFTData = Join[allFTData, <|
  "FamilyID" -> allFamily["FamilyID"],
  "PipelineRequestID" -> allRequest["RequestID"],
  "AllSelectionRequestID" -> allResolution["SelectionRequestID"],
  "OutputResolutionID" -> allResolution["ResolutionID"],
  "ResolvedOutputRequests" -> allResolution["OutputRequests"]
|>];
allFTData["Levels"][1] = <|
  "Topology" -> <|"SetupFingerprintRecord" -> setup|>,
  "Masters" -> {{1}}, "DiffMatrix" -> {{0}}, "Computed" -> True|>;
assert["prepared All data proves the exact resolved basis and identities",
  preparedFTDataMatchesContractQ[allFTData, allContract] &&
    !preparedFTDataMatchesContractQ[
      Join[allFTData, <|"OutputResolutionID" -> "stale"|>], allContract]];

allDiscoveryContract = ft2AllDiscoveryContractRecord[allRunName, allRequest];
allDiscoveryFile = FileNameJoin[{requestDirectory, "all-resolution.wxf"}];
allDiscoveryWrite = ft2SaveAllDiscoveryCache[
  allDiscoveryFile, allRequest, allDiscoveryContract, allResolution];
assert["All discovery cache round-trips an exact source/runtime contract",
  StringQ[allDiscoveryWrite] &&
    allDiscoveryContract["WolframRuntime"] ===
      <|"Version" -> $Version, "SystemID" -> $SystemID|> &&
    MemberQ[Lookup[allDiscoveryContract["DiscoverySources"],
      "RelativePath", None], "Scripts/run_ft_stepwise2.m"] &&
    !AnyTrue[allDiscoveryContract["DiscoverySources"], FailureQ] &&
    SameQ[ft2LoadAllDiscoveryCache[
      allDiscoveryFile, allRequest, allDiscoveryContract], allResolution]];
assert["All discovery cache rejects a changed discovery contract",
  ft2LoadAllDiscoveryCache[allDiscoveryFile, allRequest,
    Join[allDiscoveryContract, <|"DimensionExpression" -> 6 - 2 FeynmanTrick`FTeps|>]] === $Failed &&
  ft2AllDiscoveryContractIdentity[allDiscoveryContract] =!=
    ft2AllDiscoveryContractIdentity[Join[allDiscoveryContract, <|
      "WolframRuntime" -> <|"Version" -> "changed", "SystemID" -> $SystemID|>
    |>]]];
discoveryCalls = 0;
cachedResolutionPair = Block[{
    prepCacheRoot = requestDirectory,
    forcePrepRebuild = False,
    ft2DiscoverAllOutputSelection = Function[theRequest,
      discoveryCalls++; allResolution]},
  {ft2ResolveAllOutputSelection[allRunName, allRequest],
   ft2ResolveAllOutputSelection[allRunName, allRequest]}];
assert["validated All discovery is reused without rerunning FIRE",
  discoveryCalls === 1 && cachedResolutionPair ===
    {allResolution, allResolution}];

$ft2ActivePipelineRequestID = request["RequestID"];
$ft2ActiveFamilyID = family["FamilyID"];
$ft2ActiveAllSelectionRequestID = None;
$ft2ActiveOutputResolutionID = None;
nativeContract = ft2NativeTransportContract[
  runName, 1, ftPrepContractKey[contract], <|"Matrix" -> {{0}}|>,
  {{0}}, {0}, {}, <||>, <||>, {}, {}, "native-plan-test"];
assert["native transport contract is request-bound",
  ft2NativeTransportContractQ[nativeContract] &&
    nativeContract["Record", "PipelineRequestID"] === request["RequestID"] &&
    nativeContract["Record", "FamilyID"] === family["FamilyID"]];

$ft2ActivePipelineRequestID = allRequest["RequestID"];
$ft2ActiveFamilyID = allFamily["FamilyID"];
$ft2ActiveAllSelectionRequestID = allResolution["SelectionRequestID"];
$ft2ActiveOutputResolutionID = allResolution["ResolutionID"];
allNativeContract = ft2NativeTransportContract[
  allRunName, 1, ftPrepContractKey[allContract], <|"Matrix" -> {{0}}|>,
  {{0}}, {0}, {}, <||>, <||>, {}, {}, "native-plan-all-test"];
assert["native/checkpoint contracts bind the resolved All basis identity",
  ft2NativeTransportContractQ[allNativeContract] &&
    ft2CheckpointRequestMetadata[] === <|
      "PipelineRequestID" -> allRequest["RequestID"],
      "FamilyID" -> allFamily["FamilyID"],
      "AllSelectionRequestID" -> allResolution["SelectionRequestID"],
      "OutputResolutionID" -> allResolution["ResolutionID"]|> &&
    allNativeContract["Record", "PipelineRequestID"] ===
      allRequest["RequestID"] &&
    allNativeContract["Record", "AllSelectionRequestID"] ===
      allResolution["SelectionRequestID"] &&
    allNativeContract["Record", "OutputResolutionID"] ===
      allResolution["ResolutionID"]];

certification = ft2NotApplicableCertification["direct", "test"];
rawValues = {
  DiffExp2`EpsSeries`ESNew[0, {7}],
  DiffExp2`EpsSeries`ESNew[0, {11}]
};
finalRows = ft2CustomFinalRows[
  family["Name"], family["FamilyID"], request["RequestID"],
  targets, rawValues, {certification, certification},
  request["OutputRequests"]];
assert["custom FINAL rows retain exact target order and identities",
  AssociationQ /@ finalRows === {True, True} &&
    Lookup[finalRows, "Master"] === targets &&
    Lookup[finalRows, "RequestOrdinal"] === {1, 2} &&
    Lookup[finalRows, "RequestID"] ===
      Lookup[request["OutputRequests"], "RequestID"] &&
    Max[Abs[Lookup[finalRows, "Finite"] - {7, 11}]] == 0];

allFinalRows = ft2CustomFinalRows[
  allFamily["Name"], allFamily["FamilyID"], allRequest["RequestID"],
  allResolution["Masters"], rawValues,
  {certification, certification}, allResolution["OutputRequests"]];
assert["resolved All FINAL rows retain selection and resolution identities",
  AssociationQ /@ allFinalRows === {True, True} &&
    Lookup[allFinalRows, "Master"] === allResolution["Masters"] &&
    Lookup[allFinalRows, "SelectionRequestID"] ===
      ConstantArray[allResolution["SelectionRequestID"], 2] &&
    Lookup[allFinalRows, "ResolutionID"] ===
      ConstantArray[allResolution["ResolutionID"], 2]];

If[DirectoryQ[requestDirectory],
  DeleteDirectory[requestDirectory, DeleteContents -> True]];
KeyValueMap[
  If[StringQ[#2], SetEnvironment[#1 -> #2], SetEnvironment[#1 -> None]] &,
  savedEnvironment];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
