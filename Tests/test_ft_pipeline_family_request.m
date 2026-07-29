(* Process-free contract tests for exact custom-family pipeline requests. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

requestDirectory = FileNameJoin[{$TemporaryDirectory,
  "ft-pipeline-request-test-" <> ToString[$ProcessID]}];
If[DirectoryQ[requestDirectory],
  DeleteDirectory[requestDirectory, DeleteContents -> True]];

rawFamily = <|
  "Name" -> "two_point_api_test",
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {Global`p},
  "Propagators" -> {
    Global`l^2 + 1,
    (Global`l + Global`p)^2 + 2
  },
  "Replacements" -> {Global`p^2 -> -1},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps,
  "CombinationSequence" -> {{1, 2}}
|>;

family = FeynmanTrick`CreateFamily[rawFamily];
orderedTargets = {{1, 1}, {2, 0}};
reverseTargets = Reverse[orderedTargets];
plan = FeynmanTrick`PipelinePlan[family, orderedTargets,
  "RequestDirectory" -> requestDirectory];
epsilonPlan = FeynmanTrick`PipelinePlan[family, orderedTargets,
  "EpsilonOrder" -> 2,
  "RequestDirectory" -> requestDirectory];
singlePlan = FeynmanTrick`PipelinePlan[family, First[orderedTargets],
  "RequestDirectory" -> requestDirectory];
samePlan = FeynmanTrick`PipelinePlan[family, orderedTargets,
  "RequestDirectory" -> requestDirectory];
reversePlan = FeynmanTrick`PipelinePlan[family, reverseTargets,
  "RequestDirectory" -> requestDirectory];
allPlan = FeynmanTrick`PipelinePlan[family, All,
  "RequestDirectory" -> requestDirectory];
registryPlan = FeynmanTrick`PipelinePlan["bubble"];

numeratorFamily = FeynmanTrick`CreateFamily[<|
  "Name" -> "declared_numerator_api_test",
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {Global`p},
  "Propagators" -> {
    Global`l^2 + 1,
    (Global`l + Global`p)^2 + 2,
    (Global`l - Global`p)^2 + 3
  },
  "Replacements" -> {Global`p^2 -> -1},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps,
  "NumeratorPositions" -> {3},
  "CombinationSequence" -> {{1, 2}}
|>, {{1, 1, -1}}];
numeratorPlan = FeynmanTrick`PipelinePlan[
  numeratorFamily, "RequestDirectory" -> requestDirectory];
wrongLengthAnchorPlan = FeynmanTrick`PipelinePlan[
  family, orderedTargets,
  "FixedParameterValues" -> {1/5, 3/10},
  "RequestDirectory" -> requestDirectory];
wrongLengthPrescriptionPlan = FeynmanTrick`PipelinePlan[
  family, orderedTargets,
  "LevelDeltaPrescriptionSigns" -> {1, -1},
  "RequestDirectory" -> requestDirectory];

assert["custom family plan is canonical and content addressed",
  AssociationQ[plan] && plan["InputKind"] === "Family" &&
    StringStartsQ[plan["RequestID"], "ft-request-"] &&
    plan["Environment", "FT_FAMILY_REQUEST_ID"] === plan["RequestID"] &&
    plan["Environment", "FT_FAMILY_REQUEST_FILE"] === plan["RequestFile"] &&
    FileNameTake[plan["RequestFile"]] === plan["RequestID"] <> ".wxf"];
assert["PipelinePlan remains process- and write-free",
  !FileExistsQ[plan["RequestFile"]]];
assert["canonical declared-numerator family is production-plannable",
  AssociationQ[numeratorPlan] &&
    numeratorPlan["Request", "Family", "NumeratorPositions"] === {3} &&
    numeratorPlan["Request", "Family", "CombinationSequence"] ===
      {{1, 2}} &&
    Lookup[numeratorPlan["OutputRequests"], "IndexVector"] ===
      {{1, 1, -1}}];
assert["custom-family anchor length is rejected before execution",
  FailureQ[wrongLengthAnchorPlan]];
assert["custom-family per-level prescription length is rejected before execution",
  FailureQ[wrongLengthPrescriptionPlan]];
assert["base plans clear inherited optional and custom-family mode variables",
  AssociationQ[registryPlan] &&
    Lookup[registryPlan["Environment"], {
      "FT_FAMILY_REQUEST_FILE", "FT_FAMILY_REQUEST_ID",
      "FT_RESUME_LADDER_CHECKPOINT", "FT_STOP_AFTER_BOUNDARY_LEVEL"},
      Missing["Absent"]] === {"", "", "", ""} &&
    Lookup[plan["Environment"], {
      "FT_FAMILY_REQUEST_FILE", "FT_FAMILY_REQUEST_ID"}] ===
      {plan["RequestFile"], plan["RequestID"]} &&
    Lookup[plan["Environment"], {
      "FT_RESUME_LADDER_CHECKPOINT", "FT_STOP_AFTER_BOUNDARY_LEVEL"},
      Missing["Absent"]] === {"", ""}];
assert["identical exact requests have identical identities",
  samePlan["RequestID"] === plan["RequestID"] &&
    samePlan["Request"] === plan["Request"]];
assert["ordered output selection participates in request identity",
  Lookup[plan["OutputRequests"], "IndexVector"] === orderedTargets &&
    Lookup[reversePlan["OutputRequests"], "IndexVector"] === reverseTargets &&
    reversePlan["RequestID"] =!= plan["RequestID"]];
assert["All is explicit, execution-time, and potentially expensive",
  allPlan["Request", "OutputMode"] === "AllPendingDiscovery" &&
    allPlan["ExecutionPolicy", "MasterDiscovery"] === "AtExecution" &&
    allPlan["ExecutionPolicy", "CostClass"] === "PotentiallyExpensive" &&
    allPlan["ExecutionPolicy", "ExecutionReady"] === True &&
    !FileExistsQ[allPlan["RequestFile"]]];

allResolution =
  FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
    allPlan["Request"], {{2, 0}, {1, 1}, {2, 0}}];
assert["All resolution sorts and deduplicates FIRE masters deterministically",
  AssociationQ[allResolution] &&
    allResolution["Masters"] === {{1, 1}, {2, 0}} &&
    SameQ[allResolution,
      FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
        allPlan["Request"], {{1, 1}, {2, 0}}]] &&
    Lookup[allResolution["OutputRequests"], "RequestOrdinal"] === {1, 2} &&
    DuplicateFreeQ[Lookup[allResolution["OutputRequests"], "RequestID"]] &&
    AllTrue[allResolution["OutputRequests"],
      #["SelectionRequestID"] ===
        allPlan["OutputRequests"][[1, "RequestID"]] &&
      #["ResolutionID"] === allResolution["ResolutionID"] &]];
assert["All resolution preserves the original content-addressed selection",
  allResolution["PipelineRequestID"] === allPlan["RequestID"] &&
    allPlan["Request", "OutputRequests"] === allPlan["OutputRequests"] &&
    allPlan["OutputRequests"][[1, "IndexVector"]] === All &&
    TrueQ[FeynmanTrick`PipelineRequest`ResolvedAllOutputSelectionQ[
      allPlan["Request"], allResolution]]];
assert["unsupported discovered arity and merged numerators fail closed",
  FailureQ[FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
    allPlan["Request"], {{1, 1, 0}}]] &&
  FailureQ[FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
    allPlan["Request"], {{-1, 1}}]]];

singlePropagatorFamily = FeynmanTrick`CreateFamily[<|
  "Name" -> "single_propagator_api_test",
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {},
  "Propagators" -> {Global`l^2 + 1},
  "Replacements" -> {},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps
|>, {{1}}];
assert["zero-level single-propagator ladders fail at planning time",
  FailureQ[FeynmanTrick`PipelinePlan[
    singlePropagatorFamily, {{1}},
    "RequestDirectory" -> requestDirectory]]];

assert["the built-in runner is request-aware for explicit targets",
  plan["ExecutionPolicy", "RunnerSupport"] ===
      "BuiltInRequestAwareRunner" &&
    plan["ExecutionPolicy", "ExecutionReady"] === True];
assert["planning All never launches a child or writes its request",
  !FileExistsQ[plan["RequestFile"]] &&
    !FileExistsQ[allPlan["RequestFile"]]];

fakeRunner = CreateTemporary[];
readyPlan = FeynmanTrick`PipelinePlan[family, orderedTargets,
  "Runner" -> fakeRunner,
  "RequestAwareRunner" -> True,
  "RequestDirectory" -> requestDirectory];
preparedPlan =
  FeynmanTrick`DiffExp2Pipeline`Private`preparePlanRequest[readyPlan];
roundTrip = If[AssociationQ[preparedPlan],
  FeynmanTrick`PipelineRequest`ReadPipelineRequest[
    preparedPlan["RequestFile"], preparedPlan["RequestID"]], $Failed];
assert["request-aware child handoff materializes exact WXF atomically",
  AssociationQ[preparedPlan] && AssociationQ[roundTrip] &&
    SameQ[roundTrip, preparedPlan["Request"]]];
preparedAllPlan =
  FeynmanTrick`DiffExp2Pipeline`Private`preparePlanRequest[allPlan];
assert["execution preparation materializes an execution-ready All request",
  AssociationQ[preparedAllPlan] &&
    SameQ[
      FeynmanTrick`PipelineRequest`ReadPipelineRequest[
        preparedAllPlan["RequestFile"], preparedAllPlan["RequestID"]],
      preparedAllPlan["Request"]]];
tampered = Join[preparedPlan["Request"], <|"RequestID" -> "tampered"|>];
assert["request identity validation rejects tampering",
  !TrueQ[FeynmanTrick`PipelineRequest`PipelineRequestQ[tampered]]];
staleFamily = Join[family, <|
  "FamilyID" -> "ft-family-stale",
  "FamilyIdentity" -> "ft-family-stale"|>];
assert["request creation rejects stale internal family identities",
  FailureQ[FeynmanTrick`PipelineRequest`CreatePipelineRequest[staleFamily]]];
topologyOnlyPrescription = Join[family, <|
  "Topology" -> Join[family["Topology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>],
  "TopTopology" -> Join[family["TopTopology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>]
|>];
assert["request validation rejects topology-only mathematical fields",
  FailureQ[FeynmanTrick`PipelineRequest`CreatePipelineRequest[
    topologyOnlyPrescription]]];
assumedFamily = FeynmanTrick`CreateFamily[
  Join[rawFamily, <|"KinematicAssumptions" -> (Global`s < 0)|>],
  orderedTargets];
assert["unwired kinematic assumptions fail at planning time",
  AssociationQ[assumedFamily] &&
    FailureQ[FeynmanTrick`PipelinePlan[
      assumedFamily, orderedTargets,
      "RequestDirectory" -> requestDirectory]]];
threeFamily = FeynmanTrick`CreateFamily[<|
  "Name" -> "three_sequence_validation",
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
assert["request validation replays active-position sequence invariants",
  FailureQ[FeynmanTrick`PipelineRequest`CreatePipelineRequest[
    invalidSequenceFamily]]];
assert["plan duplicate metadata is checked before child launch",
  FailureQ[FeynmanTrick`DiffExp2Pipeline`Private`preparePlanRequest[
    Join[readyPlan, <|"FamilyID" -> "ft-family-stale"|>]]]];

validCertification = <|
  "Applicability" -> "not-applicable",
  "Operation" -> "direct",
  "Reason" -> "process-free-test"
|>;
validFinalRow[thePlan_, ordinal_Integer, finite_Integer] := Module[
  {outputRequest = thePlan["OutputRequests"][[ordinal]], row,
   epsilonOrder = thePlan["Settings", "EpsilonOrder"]},
  row = <|
    "Example" -> thePlan["FamilyName"],
    "RequestOrdinal" -> ordinal,
    "RequestID" -> outputRequest["RequestID"],
    "PhysicalIntegralID" -> outputRequest["PhysicalIntegralID"],
    "Master" -> outputRequest["IndexVector"],
    "PipelineRequestID" -> thePlan["RequestID"],
    "FamilyID" -> thePlan["FamilyID"],
    "Finite" -> finite,
    "RawMinPower" -> 0,
    "Certification" -> validCertification
  |>;
  If[epsilonOrder > 0,
    Append[row, "Coefficients" ->
      Table[{power, If[power === 0, finite, 0]},
        {power, 0, epsilonOrder}]],
    row]
];
finalOutput[rows_List] := StringRiffle[
  ("FINAL " <> ExportString[#, "RawJSON", "Compact" -> True] &) /@ rows,
  "\n"];
resolutionOutput[resolution_Association] :=
  "OUTPUT_RESOLUTION " <>
    ExportString[resolution, "RawJSON", "Compact" -> True];

oneRows = {validFinalRow[singlePlan, 1, 7]};
manyRows = {
  validFinalRow[plan, 1, 7],
  validFinalRow[plan, 2, 11]
};
oneProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[oneRows]|>;
manyProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[manyRows]|>;
reorderedProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[Reverse[manyRows]]|>;
malformedProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[{KeyDrop[First[manyRows], "FamilyID"],
    Last[manyRows]}]|>;
missingPayloadProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[{
    KeyDrop[First[manyRows], {"RawMinPower", "Certification"}],
    Last[manyRows]}]|>;
nonNumericPayloadProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[{
    Join[First[manyRows], <|"Finite" -> "7"|>],
    Last[manyRows]}]|>;
invalidJSONProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> "FINAL {not-json}"|>;
validPlusInvalidJSONProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[manyRows] <> "\nFINAL {not-json}"|>;
invalidResolutionJSONProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> "OUTPUT_RESOLUTION {not-json}\n" <>
    finalOutput[manyRows]|>;
stoppedReorderedProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[Reverse[manyRows]] <>
    "\nSTOPPED_AFTER_BOUNDARY_LEVEL 0"|>;
epsilonRows = {
  validFinalRow[epsilonPlan, 1, 19],
  validFinalRow[epsilonPlan, 2, 23]
};
epsilonProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[epsilonRows]|>;
epsilonMissingCoefficientsProcess = <|
  "ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[{
    KeyDrop[First[epsilonRows], "Coefficients"], Last[epsilonRows]}]|>;
epsilonTruncatedCoefficientsProcess = <|
  "ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[{
    Join[First[epsilonRows],
      <|"Coefficients" -> {{0, 19}, {1, 0}}|>],
    Last[epsilonRows]}]|>;
positiveMinEpsilonRows = ReplacePart[epsilonRows, 1 ->
  Join[First[epsilonRows], <|
    "Finite" -> 0, "RawMinPower" -> 1,
    "Coefficients" -> {{1, 29}, {2, 31}}|>]];
positiveMinEpsilonProcess = <|
  "ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[positiveMinEpsilonRows]|>;
positiveMinNonzeroFiniteProcess = <|
  "ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[
    ReplacePart[positiveMinEpsilonRows, 1 ->
      Join[First[positiveMinEpsilonRows], <|"Finite" -> 1|>]]]|>;
validAllFinalRow[ordinal_Integer, finite_Integer] := Module[
  {outputRequest = allResolution["OutputRequests"][[ordinal]]},
  <|
    "Example" -> allPlan["FamilyName"],
    "RequestOrdinal" -> ordinal,
    "RequestID" -> outputRequest["RequestID"],
    "PhysicalIntegralID" -> outputRequest["PhysicalIntegralID"],
    "Master" -> outputRequest["IndexVector"],
    "PipelineRequestID" -> allPlan["RequestID"],
    "FamilyID" -> allPlan["FamilyID"],
    "SelectionRequestID" -> allResolution["SelectionRequestID"],
    "ResolutionID" -> allResolution["ResolutionID"],
    "Finite" -> finite,
    "RawMinPower" -> 0,
    "Certification" -> validCertification
  |>
];
allRows = {validAllFinalRow[1, 13], validAllFinalRow[2, 17]};
allProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> resolutionOutput[allResolution] <> "\n" <>
    finalOutput[allRows]|>;
allReorderedProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> resolutionOutput[allResolution] <> "\n" <>
    finalOutput[Reverse[allRows]]|>;
allTamperedResolutionProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> resolutionOutput[allResolution] <> "\n" <>
    finalOutput[
    ReplacePart[allRows, 2 -> Join[allRows[[2]],
      <|"ResolutionID" -> "stale"|>]]]|>;
allMissingManifestProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> finalOutput[allRows]|>;
staleManifest = Join[allResolution, <|"ResolutionID" -> "stale"|>];
allStaleManifestProcess = <|"ExitCode" -> 0, "StandardError" -> "",
  "StandardOutput" -> resolutionOutput[staleManifest] <> "\n" <>
    finalOutput[allRows]|>;
oneResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  singlePlan, oneProcess];
manyResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  plan, manyProcess];
reorderedResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  plan, reorderedProcess];
malformedResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  plan, malformedProcess];
missingPayloadResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    plan, missingPayloadProcess];
nonNumericPayloadResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    plan, nonNumericPayloadProcess];
invalidJSONResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  plan, invalidJSONProcess];
validPlusInvalidJSONResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    plan, validPlusInvalidJSONProcess];
invalidResolutionJSONResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    plan, invalidResolutionJSONProcess];
stoppedReorderedResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    plan, stoppedReorderedProcess];
epsilonResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  epsilonPlan, epsilonProcess];
epsilonMissingCoefficientsResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    epsilonPlan, epsilonMissingCoefficientsProcess];
epsilonTruncatedCoefficientsResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    epsilonPlan, epsilonTruncatedCoefficientsProcess];
positiveMinEpsilonResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    epsilonPlan, positiveMinEpsilonProcess];
positiveMinNonzeroFiniteResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    epsilonPlan, positiveMinNonzeroFiniteProcess];
allResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  allPlan, allProcess];
allReorderedResult = FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
  allPlan, allReorderedProcess];
allTamperedResolutionResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    allPlan, allTamperedResolutionProcess];
allMissingManifestResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    allPlan, allMissingManifestProcess];
allStaleManifestResult =
  FeynmanTrick`DiffExp2Pipeline`Private`pipelineResult[
    allPlan, allStaleManifestProcess];
assert["one-output result preserves Final compatibility",
  Length[oneResult["Outputs"]] === 1 &&
    oneResult["Final"] === First[oneResult["Outputs"]] &&
    oneResult["Final", "Finite"] === 7];
assert["multiple outputs remain ordered without inventing one Final",
  Lookup[manyResult["Outputs"], "RequestOrdinal"] === {1, 2} &&
    MissingQ[manyResult["Final"]]];
assert["reordered child FINAL rows fail closed",
  reorderedResult["Status"] === "Failed" &&
    FailureQ[reorderedResult["OutputValidation"]] &&
    reorderedResult["Outputs"] === {} &&
    MissingQ[reorderedResult["Final"]] &&
    Length[reorderedResult["UnvalidatedOutputs"]] === 2];
assert["malformed child FINAL rows fail closed",
  malformedResult["Status"] === "Failed" &&
    FailureQ[malformedResult["OutputValidation"]] &&
    malformedResult["Outputs"] === {} &&
    MissingQ[malformedResult["Final"]]];
assert["custom FINAL numerical payload and certification are mandatory",
  missingPayloadResult["Status"] === "Failed" &&
    FailureQ[missingPayloadResult["OutputValidation"]] &&
    nonNumericPayloadResult["Status"] === "Failed" &&
    FailureQ[nonNumericPayloadResult["OutputValidation"]]];
assert["missing or invalid family FINAL output fails by count",
  invalidJSONResult["Status"] === "Failed" &&
    FailureQ[invalidJSONResult["OutputValidation"]] &&
    invalidJSONResult["Outputs"] === {}];
assert["malformed tagged JSON cannot hide beside an otherwise valid contract",
  validPlusInvalidJSONResult["Status"] === "Failed" &&
    FailureQ[validPlusInvalidJSONResult["OutputValidation"]] &&
    invalidResolutionJSONResult["Status"] === "Failed" &&
    FailureQ[invalidResolutionJSONResult["OutputResolutionValidation"]]];
assert["stopped marker cannot bypass FINAL identity validation",
  stoppedReorderedResult["Status"] === "Failed" &&
    FailureQ[stoppedReorderedResult["OutputValidation"]] &&
    stoppedReorderedResult["Outputs"] === {}];
assert["positive EpsilonOrder requires a complete ordered Laurent payload",
  epsilonResult["Status"] === "Succeeded" &&
    epsilonResult["Outputs"] === epsilonRows &&
    epsilonMissingCoefficientsResult["Status"] === "Failed" &&
    FailureQ[epsilonMissingCoefficientsResult["OutputValidation"]] &&
    epsilonTruncatedCoefficientsResult["Status"] === "Failed" &&
    FailureQ[epsilonTruncatedCoefficientsResult["OutputValidation"]]];
assert["positive raw minimum powers retain an honest zero finite coefficient",
  positiveMinEpsilonResult["Status"] === "Succeeded" &&
    positiveMinEpsilonResult["Outputs"] === positiveMinEpsilonRows &&
    positiveMinNonzeroFiniteResult["Status"] === "Failed" &&
    FailureQ[positiveMinNonzeroFiniteResult["OutputValidation"]]];
assert["All FINAL rows expose a validated dynamic output contract",
  allResult["Status"] === "Succeeded" &&
    allResult["Outputs"] === allRows &&
    allResult["OutputResolutionValidation"] === allResolution &&
    allResult["OutputResolution"] === allResolution &&
    allResult["ResolvedOutputRequests"] === allResolution["OutputRequests"] &&
    MissingQ[allResult["Final"]]];
assert["All FINAL order and resolution identities fail closed",
  allReorderedResult["Status"] === "Failed" &&
    FailureQ[allReorderedResult["OutputValidation"]] &&
    allTamperedResolutionResult["Status"] === "Failed" &&
    FailureQ[allTamperedResolutionResult["OutputValidation"]] &&
    allMissingManifestResult["Status"] === "Failed" &&
    FailureQ[allMissingManifestResult["OutputResolutionValidation"]] &&
    allStaleManifestResult["Status"] === "Failed" &&
    FailureQ[allStaleManifestResult["OutputResolutionValidation"]]];

If[FileExistsQ[fakeRunner], DeleteFile[fakeRunner]];
If[DirectoryQ[requestDirectory],
  DeleteDirectory[requestDirectory, DeleteContents -> True]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
