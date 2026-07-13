(* ::Package:: *)
(* Exact, content-addressed handoff records for custom FT family runs. *)

BeginPackage["FeynmanTrick`PipelineRequest`", {
  "FeynmanTrick`FamilySpec`"
}];

CreatePipelineRequest::usage =
  "CreatePipelineRequest[familySpec] builds an exact, content-addressed " <>
  "request envelope for a canonical FeynmanTrick family specification.";
PipelineRequestQ::usage =
  "PipelineRequestQ[request] returns True only for an internally consistent " <>
  "FeynmanTrick pipeline request envelope.";
PipelineRequestPath::usage =
  "PipelineRequestPath[request, directory] returns the deterministic WXF " <>
  "path for a valid content-addressed request.";
WritePipelineRequest::usage =
  "WritePipelineRequest[request, file] atomically writes a valid request as " <>
  "WXF, or reuses an identical existing content-addressed file.";
ReadPipelineRequest::usage =
  "ReadPipelineRequest[file, expectedID] reads and validates an exact WXF " <>
  "request. expectedID may be Automatic.";
CreateResolvedAllOutputSelection::usage =
  "CreateResolvedAllOutputSelection[request, masters] validates FIRE's L0 " <>
  "master list for an All request, sorts and deduplicates it deterministically, " <>
  "and returns an exact resolution record with ordered output-request identities.";
ResolvedAllOutputSelectionQ::usage =
  "ResolvedAllOutputSelectionQ[request, resolution] verifies that an All-master " <>
  "resolution record is exact, canonical, and bound to the original request.";

Begin["`Private`"];

$requestSchema = "FeynmanTrick.PipelineRequest/v1";
$requestPrefix = "ft-request-";

requestFailure[detail_String, data_:<||>] :=
  Failure["FeynmanTrickPipelineRequest", Join[<|"Detail" -> detail|>, data]];

stableDigest[expr_] := IntegerString[Hash[expr, "SHA256"], 16, 64];

exactExpressionQ[expr_] := FreeQ[Unevaluated[expr], _Real];

$familyDefinitionKeys = {
  "LoopMomenta", "ExternalMomenta", "Propagators", "Replacements",
  "NumericalPoint", "Dimension", "EliminatedPositions",
  "AnalyticPrescription", "Prescriptions", "KinematicAssumptions"
};

canonicalFamilyIdentityQ[family_Association] := Module[
  {definition, topology, expectedID, requiredDefinitionKeys,
   expectedFamilyKeys},
  definition = Lookup[family, "Definition", None];
  topology = Lookup[family, "Topology", None];
  requiredDefinitionKeys = Take[$familyDefinitionKeys, 7];
  expectedFamilyKeys = {
    "Schema", "FamilyID", "FamilyIdentity", "InputKind", "Name",
    "Definition", "Topology", "TopTopology", "NumPropagators",
    "EliminatedPositions", "NumericalPoint", "Dimension",
    "CombinationSequence", "OutputIntegralSource", "OutputIntegralMode",
    "OutputIntegrals", "L0OutputRequests"
  };
  If[!AssociationQ[definition] || !AssociationQ[topology] ||
      Sort[Keys[family]] =!= Sort[expectedFamilyKeys] ||
      !AllTrue[requiredDefinitionKeys, KeyExistsQ[definition, #] &] ||
      !AllTrue[Keys[definition], MemberQ[$familyDefinitionKeys, #] &],
    Return[False, Module]];
  expectedID = "ft-family-" <>
    stableDigest[{"FeynmanTrick.FamilyDefinition/v1", definition}];
  MemberQ[{"RawFamily", "Topology"}, Lookup[family, "InputKind", None]] &&
    StringQ[Lookup[family, "Name", None]] &&
    StringMatchQ[family["Name"],
      RegularExpression["[A-Za-z][A-Za-z0-9_]*"]] &&
    Lookup[topology, "StartFileReady", False] === False &&
    Lookup[topology, "ProblemNumber", 0] === 0 &&
    Lookup[topology, "WorkDirectory", ""] === "" &&
    Lookup[topology, "Masters", {}] === {} &&
    Lookup[topology, "MasterRules", {}] === {} &&
    !AnyTrue[{
      "SetupFingerprint", "SetupFingerprintRecord", "NumeratorPositions",
      "OriginalPropagators", "OriginalNumPropagators"},
      KeyExistsQ[topology, #] &] &&
    Lookup[family, "FamilyID", None] === expectedID &&
    Lookup[family, "FamilyIdentity", None] === expectedID &&
    Lookup[family, "TopTopology", None] === topology &&
    Lookup[family, "Name", None] === Lookup[topology, "Name", None] &&
    Lookup[family, "NumPropagators", None] ===
      Length[definition["Propagators"]] &&
    Lookup[topology, "NumPropagators", None] ===
      Lookup[family, "NumPropagators", None] &&
    Lookup[family, "NumericalPoint", None] ===
      definition["NumericalPoint"] &&
    Lookup[family, "Dimension", None] === definition["Dimension"] &&
    Lookup[family, "EliminatedPositions", None] ===
      definition["EliminatedPositions"] &&
    (* Optional mathematical fields cannot be smuggled into only the topology:
       the exact Definition is the sole family-identity source. *)
    AllTrue[$familyDefinitionKeys,
      Lookup[topology, #, Missing["Absent"]] ===
        Lookup[definition, #, Missing["Absent"]] &]
];

canonicalCombinationSequenceQ[sequence_, n_Integer,
    eliminated_List] := Module[{active, pair},
  If[!ListQ[sequence] || Length[sequence] > Max[0, n - Length[eliminated] - 1] ||
      !AllTrue[sequence, MatchQ[#, {_Integer, _Integer}] &],
    Return[False, Module]];
  active = Complement[Range[n], eliminated];
  Do[
    pair = sequence[[step]];
    If[pair[[1]] === pair[[2]] ||
        !MemberQ[active, pair[[1]]] || !MemberQ[active, pair[[2]]],
      Return[False, Module]];
    active = DeleteCases[active, pair[[2]]],
    {step, Length[sequence]}];
  True
];

explicitOutputRequestQ[request_, familyID_String, output_List,
    ordinal_Integer] := Module[{expectedID},
  expectedID = "ft-l0-" <> stableDigest[{
    "FeynmanTrick.L0OutputRequest/v1", familyID, output}];
  AssociationQ[request] &&
    Lookup[request, "Schema", None] ===
      "FeynmanTrick.L0OutputRequest/v1" &&
    Lookup[request, "RequestID", None] === expectedID &&
    Lookup[request, "RequestIdentity", None] === expectedID &&
    Lookup[request, "PhysicalIntegralID", None] === expectedID &&
    Lookup[request, "RequestOrdinal", None] === ordinal &&
    Lookup[request, "IndexVector", None] === output &&
    TrueQ[Lookup[request, "Resolved", False]]
];

allOutputRequestQ[request_, familyID_String] := Module[{expectedID},
  expectedID = "ft-l0-all-" <> stableDigest[{
    "FeynmanTrick.L0OutputRequest/v1", familyID, All}];
  AssociationQ[request] &&
    Lookup[request, "Schema", None] ===
      "FeynmanTrick.L0OutputRequest/v1" &&
    Lookup[request, "RequestID", None] === expectedID &&
    Lookup[request, "RequestIdentity", None] === expectedID &&
    Lookup[request, "PhysicalIntegralID", None] === expectedID &&
    Lookup[request, "RequestOrdinal", None] === 1 &&
    Lookup[request, "IndexVector", None] === All &&
    Lookup[request, "Selection", None] === All &&
    TrueQ[Lookup[request, "Resolved", True] === False]
];

familyRequestShapeQ[family_] := Module[
  {mode, outputs, requests, n, eliminated, sequence, validatedOutputs,
   replayed},
  If[!AssociationQ[family] ||
      Lookup[family, "Schema", None] =!= "FeynmanTrick.FamilySpec/v1" ||
      !StringQ[Lookup[family, "FamilyID", None]] ||
      !IntegerQ[Lookup[family, "NumPropagators", None]] ||
      !AssociationQ[Lookup[family, "Topology", None]] ||
      !ListQ[Lookup[family, "CombinationSequence", None]] ||
      !exactExpressionQ[family] || !canonicalFamilyIdentityQ[family],
    Return[False, Module]];
  n = family["NumPropagators"];
  eliminated = family["EliminatedPositions"];
  sequence = family["CombinationSequence"];
  If[!ListQ[eliminated] || !DuplicateFreeQ[eliminated] ||
      !AllTrue[eliminated, IntegerQ[#] && 1 <= # <= n &] ||
      !canonicalCombinationSequenceQ[sequence, n, eliminated],
    Return[False, Module]];
  mode = Lookup[family, "OutputIntegralMode", None];
  outputs = Lookup[family, "OutputIntegrals", None];
  requests = Lookup[family, "L0OutputRequests", None];
  If[!MemberQ[{"Explicit", "AllPendingDiscovery"}, mode] ||
      !ListQ[requests] || requests === {}, Return[False, Module]];
  (* Replay the authoritative process-free constructor.  InputKind and the
     descriptive target-source tag can legitimately differ because this
     replay starts from the already canonical unprepared topology; every
     mathematical field, topology field, sequence, target and identity must
     reproduce byte-for-byte. *)
  replayed = Quiet[Check[
    FeynmanTrick`FamilySpec`CreateFamily[
      family["Topology"], If[mode === "AllPendingDiscovery", All, outputs],
      "Name" -> family["Name"],
      "CombinationSequence" -> sequence], $Failed]];
  If[!AssociationQ[replayed] ||
      KeyDrop[replayed, {"InputKind", "OutputIntegralSource"}] =!=
        KeyDrop[family, {"InputKind", "OutputIntegralSource"}],
    Return[False, Module]];
  Switch[mode,
    "Explicit",
      validatedOutputs = Quiet[
        FeynmanTrick`FamilySpec`NormalizeOutputIntegrals[
          outputs, n, eliminated]];
      If[validatedOutputs =!= $Failed,
        validatedOutputs = Quiet[
          FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence[
            validatedOutputs, sequence, eliminated]]];
      validatedOutputs === outputs && ListQ[outputs] && outputs =!= {} &&
        Length[outputs] === Length[requests] &&
        AllTrue[outputs, ListQ[#] && Length[#] === n &&
          AllTrue[#, IntegerQ] &] &&
        And @@ MapThread[
          Function[{request, output, ordinal},
            explicitOutputRequestQ[
              request, family["FamilyID"], output, ordinal]],
          {requests, outputs, Range[Length[outputs]]}],
    "AllPendingDiscovery",
      outputs === All && Length[requests] === 1 &&
        allOutputRequestQ[First[requests], family["FamilyID"]],
    _, False
  ]
];

executionPolicy["Explicit"] := <|
  "MasterDiscovery" -> "NotRequired",
  "ExecutionReady" -> True,
  "CostClass" -> "RequestedIntegrals"
|>;
executionPolicy["AllPendingDiscovery"] := <|
  "MasterDiscovery" -> "AtExecution",
  "ExecutionReady" -> True,
  "CostClass" -> "PotentiallyExpensive"
|>;

requestCore[family_Association] := <|
  "Schema" -> $requestSchema,
  "Family" -> family,
  "OutputMode" -> family["OutputIntegralMode"],
  "OutputRequests" -> family["L0OutputRequests"],
  "ExecutionPolicy" -> executionPolicy[family["OutputIntegralMode"]]
|>;

requestIdentity[core_Association] := $requestPrefix <>
  stableDigest[{"FeynmanTrick.PipelineRequestIdentity/v1", core}];

CreatePipelineRequest[family_Association] := Module[{core},
  If[!familyRequestShapeQ[family],
    Return[requestFailure[
      "the family specification is malformed or contains inexact data"],
      Module]];
  core = requestCore[family];
  Append[core, "RequestID" -> requestIdentity[core]]
];

CreatePipelineRequest[input_] := requestFailure[
  "CreatePipelineRequest requires a canonical FamilySpec/v1 association",
  <|"Input" -> HoldForm[input]|>];

resolvedAllPhysicalIntegralID[familyID_String, vector_List] :=
  "ft-l0-" <> stableDigest[{
    "FeynmanTrick.L0OutputRequest/v1", familyID, vector}];

resolvedAllRequestRecord[familyID_String, selectionRequestID_String,
    resolutionID_String, vector_List, ordinal_Integer] := Module[
  {requestID, physicalID},
  physicalID = resolvedAllPhysicalIntegralID[familyID, vector];
  requestID = "ft-l0-resolved-" <> stableDigest[{
    "FeynmanTrick.ResolvedL0OutputRequest/v1", resolutionID,
    ordinal, vector}];
  <|
    "Schema" -> "FeynmanTrick.ResolvedL0OutputRequest/v1",
    "RequestID" -> requestID,
    "RequestIdentity" -> requestID,
    "PhysicalIntegralID" -> physicalID,
    "RequestOrdinal" -> ordinal,
    "IndexVector" -> vector,
    "Resolved" -> True,
    "SelectionRequestID" -> selectionRequestID,
    "ResolutionID" -> resolutionID
  |>
];

createResolvedAllOutputSelection[request_Association, masters_] := Module[
  {family, n, normalized, canonical, validated, selectionRequestID,
   resolutionID, outputRequests},
  If[Lookup[request, "OutputMode", None] =!= "AllPendingDiscovery" ||
      !ListQ[masters] || masters === {},
    Return[requestFailure[
      "All-master discovery must return a nonempty list of index vectors",
      <|"Masters" -> masters|>], Module]];
  family = request["Family"];
  n = family["NumPropagators"];
  If[!AllTrue[masters,
      ListQ[#] && Length[#] === n && AllTrue[#, IntegerQ] &],
    Return[requestFailure[
      "discovered masters contain an unsupported arity or noninteger index",
      <|"ExpectedArity" -> n, "Masters" -> masters|>], Module]];
  normalized = Quiet[
    FeynmanTrick`FamilySpec`NormalizeOutputIntegrals[
      masters, n, family["EliminatedPositions"]]];
  If[normalized === $Failed,
    Return[requestFailure[
      "discovered masters violate the family's eliminated-position contract",
      <|"Masters" -> masters,
        "EliminatedPositions" -> family["EliminatedPositions"]|>], Module]];
  validated = Quiet[
    FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence[
      normalized, family["CombinationSequence"],
      family["EliminatedPositions"]]];
  If[validated === $Failed,
    Return[requestFailure[
      "discovered masters contain numerator powers at a Feynman-trick merge position",
      <|"Masters" -> masters,
        "CombinationSequence" -> family["CombinationSequence"]|>], Module]];
  canonical = Sort[DeleteDuplicates[validated]];
  selectionRequestID = request["OutputRequests"][[1, "RequestID"]];
  resolutionID = "ft-all-resolution-" <> stableDigest[{
    "FeynmanTrick.ResolvedAllOutputSelection/v1", request["RequestID"],
    selectionRequestID, canonical}];
  outputRequests = MapIndexed[
    resolvedAllRequestRecord[family["FamilyID"], selectionRequestID,
      resolutionID, #1, First[#2]] &,
    canonical];
  <|
    "Schema" -> "FeynmanTrick.ResolvedAllOutputSelection/v1",
    "PipelineRequestID" -> request["RequestID"],
    "SelectionRequestID" -> selectionRequestID,
    "FamilyID" -> family["FamilyID"],
    "OutputMode" -> "AllResolved",
    "Masters" -> canonical,
    "OutputRequests" -> outputRequests,
    "ResolutionID" -> resolutionID
  |>
];

CreateResolvedAllOutputSelection[request_, masters_] := Module[{},
  If[!TrueQ[PipelineRequestQ[request]],
    Return[requestFailure[
      "All-master resolution requires a valid pipeline request"], Module]];
  createResolvedAllOutputSelection[request, masters]
];

ResolvedAllOutputSelectionQ[request_, resolution_] := Module[{expected},
  If[!TrueQ[PipelineRequestQ[request]] || !AssociationQ[resolution] ||
      !exactExpressionQ[resolution] ||
      Sort[Keys[resolution]] =!= Sort[{
        "Schema", "PipelineRequestID", "SelectionRequestID", "FamilyID",
        "OutputMode", "Masters", "OutputRequests", "ResolutionID"}],
    Return[False, Module]];
  expected = createResolvedAllOutputSelection[
    request, Lookup[resolution, "Masters", None]];
  AssociationQ[expected] && SameQ[expected, resolution]
];

PipelineRequestQ[request_] := Module[{core, expectedKeys},
  If[!AssociationQ[request] || !exactExpressionQ[request], Return[False, Module]];
  expectedKeys = {
    "Schema", "Family", "OutputMode", "OutputRequests",
    "ExecutionPolicy", "RequestID"
  };
  If[Sort[Keys[request]] =!= Sort[expectedKeys] ||
      Lookup[request, "Schema", None] =!= $requestSchema ||
      !familyRequestShapeQ[Lookup[request, "Family", None]],
    Return[False, Module]];
  core = KeyTake[request, Most[expectedKeys]];
  request["OutputMode"] === request["Family", "OutputIntegralMode"] &&
    request["OutputRequests"] === request["Family", "L0OutputRequests"] &&
    request["ExecutionPolicy"] === executionPolicy[request["OutputMode"]] &&
    request["RequestID"] === requestIdentity[core]
];

PipelineRequestPath[request_, directory_String] := Module[{dir},
  If[!TrueQ[PipelineRequestQ[request]],
    Return[requestFailure["cannot form a path for an invalid request"], Module]];
  If[StringLength[StringTrim[directory]] === 0,
    Return[requestFailure["request directory must be a nonempty string"], Module]];
  dir = ExpandFileName[StringTrim[directory]];
  FileNameJoin[{dir, request["RequestID"] <> ".wxf"}]
];

PipelineRequestPath[request_, directory_] := requestFailure[
  "request directory must be a nonempty string",
  <|"Directory" -> directory|>];

ReadPipelineRequest[file_String, expectedID_:Automatic] := Module[{path, value},
  path = ExpandFileName[file];
  If[!FileExistsQ[path],
    Return[requestFailure["pipeline request file does not exist",
      <|"File" -> path|>], Module]];
  value = Quiet[Check[Import[path, "WXF"], $Failed]];
  If[!TrueQ[PipelineRequestQ[value]],
    Return[requestFailure["pipeline request file is malformed or has a stale digest",
      <|"File" -> path|>], Module]];
  If[expectedID =!= Automatic && value["RequestID"] =!= expectedID,
    Return[requestFailure["pipeline request identity does not match the expected identity",
      <|"File" -> path, "ExpectedRequestID" -> expectedID,
        "ActualRequestID" -> value["RequestID"]|>], Module]];
  value
];

ReadPipelineRequest[file_, expectedID_:Automatic] := requestFailure[
  "pipeline request file must be a path string",
  <|"File" -> file, "ExpectedRequestID" -> expectedID|>];

WritePipelineRequest[request_, file_String] := Module[
  {path, directory, existing, tmp, written, loaded},
  If[!TrueQ[PipelineRequestQ[request]],
    Return[requestFailure["refusing to write an invalid pipeline request"], Module]];
  path = ExpandFileName[file];
  directory = DirectoryName[path];
  If[!DirectoryQ[directory],
    Quiet[Check[
      CreateDirectory[directory, CreateIntermediateDirectories -> True],
      Return[requestFailure["could not create the pipeline request directory",
        <|"Directory" -> directory|>], Module]]]];
  If[FileExistsQ[path],
    existing = ReadPipelineRequest[path, request["RequestID"]];
    If[AssociationQ[existing] && SameQ[existing, request], Return[path, Module]];
    Return[requestFailure[
      "an incompatible file already occupies the content-addressed request path",
      <|"File" -> path|>], Module]
  ];
  tmp = path <> ".tmp-" <> ToString[$ProcessID] <> "-" <>
    IntegerString[Hash[{AbsoluteTime[], $ProcessID}], 16];
  written = Quiet[Check[Export[tmp, request, "WXF"], $Failed]];
  If[written === $Failed || !FileExistsQ[tmp],
    Return[requestFailure["could not write the pipeline request",
      <|"File" -> path|>], Module]];
  loaded = ReadPipelineRequest[tmp, request["RequestID"]];
  If[!AssociationQ[loaded] || !SameQ[loaded, request],
    Quiet[DeleteFile[tmp]];
    Return[requestFailure["the written pipeline request did not round-trip exactly",
      <|"File" -> path|>], Module]
  ];
  If[FileExistsQ[path],
    existing = ReadPipelineRequest[path, request["RequestID"]];
    Quiet[DeleteFile[tmp]];
    If[AssociationQ[existing] && SameQ[existing, request], Return[path, Module]];
    Return[requestFailure[
      "a concurrent incompatible request occupied the content-addressed path",
      <|"File" -> path|>], Module]
  ];
  If[Quiet[Check[RenameFile[tmp, path]; True, False]], path,
    Quiet[If[FileExistsQ[tmp], DeleteFile[tmp]]];
    requestFailure["could not publish the pipeline request atomically",
      <|"File" -> path|>]]
];

WritePipelineRequest[request_, file_] := requestFailure[
  "pipeline request file must be a path string", <|"File" -> file|>];

End[];
EndPackage[];
