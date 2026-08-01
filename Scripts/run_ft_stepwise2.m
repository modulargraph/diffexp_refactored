(* DiffExp2 FT cutover runner: the M5 ladder.
   Mirrors Scripts/run_ft_stepwise.m but the transport/integration chain is
   DiffExp2 (sector-native, no fitting): per level the in-memory exact
   DiffMatrix is loaded directly (no slice export round-trip), boundary
   C++ mode prepares both retained arms once and evaluates every level
   observable in one native batch; explicit Wolfram mode retains the legacy
   LineIntegral / EndpointLimitValues / direct chain.  Output:
   STEPWISE/FINAL rows compatible with the old comparator. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
(* FeynmanTrick.m loads the root DiffExp2 umbrella (and therefore every
   implementation module) exactly once.  Reloading the implementation here
   would reset configuration and risks same-name context capture. *)
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

ft2ParsePreparationOnly[value_String] := Switch[value,
  "0", False,
  "1", True,
  _, Failure["FeynmanTrickPreparationOnly", <|
    "Detail" -> "FT_PREPARATION_ONLY must be 0 or 1",
    "Value" -> value|>]];
ft2PreparationOnly =
  ft2ParsePreparationOnly[envOrDefault["FT_PREPARATION_ONLY", "0"]];
If[FailureQ[ft2PreparationOnly],
  Print["Invalid Feynman-trick runner environment: ", ft2PreparationOnly];
  Exit[2]];

(* The package facade and this script share one strict parser/default set.
   In particular, Cpp is now the release default; the Wolfram recurrence is
   still available only by explicit selection. *)
runnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
If[FailureQ[runnerSettings],
  Print["Invalid Feynman-trick runner environment: ", runnerSettings];
  Exit[2]];

(* A custom family is handed to the child as exact WXF plus its independent
   content identity.  Read and validate the pair before any FIRE preparation
   can start.  Registry execution is the exact historical no-variable path. *)
ft2LoadFamilyRequest[file_String, expectedID_String] := Module[
  {request, family, activeDenominatorCount, unwiredFields},
  request = FeynmanTrick`PipelineRequest`ReadPipelineRequest[file, expectedID];
  If[FailureQ[request], Return[request, Module]];
  family = request["Family"];
  If[!MemberQ[{"Explicit", "AllPendingDiscovery"},
        Lookup[request, "OutputMode", None]] ||
      !TrueQ[Lookup[Lookup[request, "ExecutionPolicy", <||>],
        "ExecutionReady", False]],
    Return[Failure["FeynmanTrickFamilyRequest", <|
      "Detail" -> "the production runner requires explicit targets or an execution-ready All-master request",
      "RequestID" -> Lookup[request, "RequestID", None],
      "OutputMode" -> Lookup[request, "OutputMode", None]|>], Module]];
  activeDenominatorCount = family["NumPropagators"] -
    Length[Union[family["EliminatedPositions"],
      family["NumeratorPositions"]]];
  If[activeDenominatorCount <= 1 ||
      Length[family["CombinationSequence"]] =!=
        activeDenominatorCount - 1,
    Return[Failure["FeynmanTrickFamilyRequest", <|
      "Detail" -> "the production ladder requires a complete nonempty combination sequence ending in one active denominator",
      "ActiveDenominators" -> activeDenominatorCount,
      "NumeratorPositions" -> family["NumeratorPositions"],
      "CombinationSequence" -> family["CombinationSequence"]|>], Module]];
  unwiredFields = Select[
    {"AnalyticPrescription", "Prescriptions", "KinematicAssumptions"},
    KeyExistsQ[family["Definition"], #] &];
  If[unwiredFields =!= {},
    Return[Failure["FeynmanTrickFamilyRequest", <|
      "Detail" -> "custom branch or kinematic-assumption fields are not yet wired into the production runner and cannot be ignored safely",
      "UnsupportedFields" -> unwiredFields,
      "RequestID" -> request["RequestID"]|>], Module]];
  request
];

ft2ResolveFamilyRequest[file_, expectedID_] := Which[
  file === "" && expectedID === "", None,
  !StringQ[file] || !StringQ[expectedID] || file === "" || expectedID === "",
    Failure["FeynmanTrickFamilyRequest", <|
      "Detail" -> "FT_FAMILY_REQUEST_FILE and FT_FAMILY_REQUEST_ID must be supplied together"|>],
  True, ft2LoadFamilyRequest[ExpandFileName[file], expectedID]
];

ft2FamilyRequestFile = envOrDefault["FT_FAMILY_REQUEST_FILE", ""];
ft2FamilyRequestID = envOrDefault["FT_FAMILY_REQUEST_ID", ""];
ft2FamilyRequest = ft2ResolveFamilyRequest[
  ft2FamilyRequestFile, ft2FamilyRequestID];
If[FailureQ[ft2FamilyRequest],
  Print["Invalid custom Feynman-trick request: ", ft2FamilyRequest];
  Exit[2]];
ft2FamilyRunName = If[AssociationQ[ft2FamilyRequest],
  "family_" <> StringTake[ft2FamilyRequest["RequestID"], -16], None];
If[AssociationQ[ft2FamilyRequest] &&
    envOrDefault["FT_EXAMPLES", "bubble"] =!= ft2FamilyRunName,
  Print["Invalid custom Feynman-trick request: FT_EXAMPLES must equal the content-addressed internal run name ",
    ft2FamilyRunName];
  Exit[2]];

(* Restore the exact FIRE7 backend contract serialized by PipelinePlan before
   any family preparation or reduction-cache identity is computed. *)
Scan[(FeynmanTrick`SetFTOption[#, runnerSettings[#]]) &, {
  "FIREPath", "FIREBackend", "FIRECalc", "FIREModularWorkers",
  "FIREUseMultiprime", "FIREPrimeLimit", "FIREKeepModularTables",
  "FIREDimensionSeparated", "FIREMultiprimeWidth", "FIREMPIExecutable",
  "FIREBasisProbeCount", "FIREModularCacheDirectory",
  "FIRETimeoutSeconds"
}];
If[runnerSettings["FixedParameterValues"] =!= Automatic,
  FeynmanTrick`SetFTOption[
    "FixedParameterValues", runnerSettings["FixedParameterValues"]]];

singularMatchPrecondition = runnerSettings["SingularMatchPrecondition"];
valueTransport = runnerSettings["ValueTransport"];
valueTransportMode = If[TrueQ[valueTransport], "1", "0"];
nativeValueHopExecution = runnerSettings["NativeValueHopExecution"];
nativeValueHopExecutionMode =
  If[TrueQ[nativeValueHopExecution], "1", "0"];
recurrenceBackend = runnerSettings["RecurrenceBackend"];
cppBatchEndpointArms = runnerSettings["BatchEndpointArms"];
cppArmThreadBudget = runnerSettings["CppThreads"];
observableContractionThreadsText = envOrDefault[
  "FT_OBSERVABLE_CONTRACTION_THREADS",
  ToString[Min[2, cppArmThreadBudget]]];
If[!StringMatchQ[observableContractionThreadsText,
    DigitCharacter..] ||
    !TrueQ[1 <= FromDigits[observableContractionThreadsText] <= 32],
  Print["Invalid FT_OBSERVABLE_CONTRACTION_THREADS: expected an integer from 1 through 32"];
  Exit[2]];
observableContractionThreads =
  FromDigits[observableContractionThreadsText];
observableContractionChunkSizeText = envOrDefault[
  "FT_OBSERVABLE_CONTRACTION_CHUNK_SIZE",
  ToString[observableContractionThreads]];
If[!StringMatchQ[observableContractionChunkSizeText,
    DigitCharacter..] ||
    !TrueQ[1 <= FromDigits[observableContractionChunkSizeText] <= 32],
  Print["Invalid FT_OBSERVABLE_CONTRACTION_CHUNK_SIZE: expected an integer from 1 through 32"];
  Exit[2]];
observableContractionChunkSize =
  FromDigits[observableContractionChunkSizeText];
deltaPrescriptionSign = runnerSettings["DeltaPrescriptionSign"];
configuredLevelDeltaPrescriptionSigns =
  runnerSettings["LevelDeltaPrescriptionSigns"];
levelDeltaPrescriptionSigns = configuredLevelDeltaPrescriptionSigns;
DiffExp2`Transport`Private`$enableSingularMatchPrecondition =
  singularMatchPrecondition;
If[singularMatchPrecondition,
  Print["DE2 singular match precondition enabled"]];
wp = runnerSettings["WorkingPrecision"];
matchDigits = runnerSettings["MatchingDigits"];
matchingCertificationSafetyDigitsText = envOrDefault[
  "FT_MATCHING_CERTIFICATION_SAFETY_DIGITS",
  ToString[DiffExp2`Tolerances`$SafetyDigits]];
If[!StringMatchQ[matchingCertificationSafetyDigitsText,
    DigitCharacter..],
  Print["Invalid FT_MATCHING_CERTIFICATION_SAFETY_DIGITS: expected a nonnegative integer"];
  Exit[2]];
matchingCertificationSafetyDigits =
  FromDigits[matchingCertificationSafetyDigitsText];
matchingCertificationMaximumDigitsText = envOrDefault[
  "FT_MATCHING_CERTIFICATION_MAX_DIGITS", ToString[wp]];
If[!StringMatchQ[matchingCertificationMaximumDigitsText,
    DigitCharacter..] ||
    !TrueQ[1 <= FromDigits[matchingCertificationMaximumDigitsText] <= wp],
  Print["Invalid FT_MATCHING_CERTIFICATION_MAX_DIGITS: expected an integer from 1 through FT_WORKING_PRECISION"];
  Exit[2]];
matchingCertificationMaximumDigits =
  FromDigits[matchingCertificationMaximumDigitsText];
ft2ParseLevelMatchingCertificationDigits[value_String,
    maximum_Integer] := Module[{tokens, parts, rules},
  If[StringLength[StringTrim[value]] === 0, Return[<||>, Module]];
  tokens = StringTrim /@ StringSplit[value, ","];
  parts = (StringTrim /@ StringSplit[#, ":"]) & /@ tokens;
  If[AnyTrue[parts,
      Length[#] =!= 2 ||
        !And @@ (StringMatchQ[#, DigitCharacter..] & /@ #) &],
    Return[Failure["FeynmanTrickMatchingCertificationDigits", <|
      "Detail" ->
        "FT_MATCHING_CERTIFICATION_DIGITS_BY_LEVEL must be a comma-separated level:digits map",
      "Value" -> value|>], Module]];
  rules = (FromDigits[#[[1]]] -> FromDigits[#[[2]]]) & /@ parts;
  If[Length[DeleteDuplicates[First /@ rules]] =!= Length[rules] ||
      AnyTrue[rules,
        First[#] < 1 || Last[#] < 1 || Last[#] > maximum &],
    Return[Failure["FeynmanTrickMatchingCertificationDigits", <|
      "Detail" ->
        "per-level certification entries require unique positive levels and digits within working precision",
      "Value" -> value, "WorkingPrecision" -> maximum|>], Module]];
  Association[rules]
];
matchingCertificationDigitsByLevel =
  ft2ParseLevelMatchingCertificationDigits[
    envOrDefault[
      "FT_MATCHING_CERTIFICATION_DIGITS_BY_LEVEL", ""], wp];
If[FailureQ[matchingCertificationDigitsByLevel],
  Print["Invalid FT_MATCHING_CERTIFICATION_DIGITS_BY_LEVEL: ",
    matchingCertificationDigitsByLevel];
  Exit[2]];
ft2MatchingCertificationComputationDigits[
    certificationByLevel_Association] :=
  Select[Map[
    Min[wp, Max[matchDigits,
      # + matchingCertificationSafetyDigits]] &,
    certificationByLevel], # > matchDigits &];
matchingCertificationComputationDigitsByLevel =
  ft2MatchingCertificationComputationDigits[
    matchingCertificationDigitsByLevel];
epsOrder = runnerSettings["EpsilonOrder"];
expansionOrder = runnerSettings["ExpansionOrder"];
boundaryExtraOrder = runnerSettings["BoundaryExtraOrder"];
divisionOrder = runnerSettings["DivisionOrder"];
stopAfterBoundaryLevel = runnerSettings["StopAfterBoundaryLevel"];
radiusOfConvergence = runnerSettings["RadiusOfConvergence"];
stepDivisionOrder = runnerSettings["StepDivisionOrder"];
If[runnerSettings["RequestedStepDivisionOrder"] =!= divisionOrder,
  Print["FT_STEP_DIVISION_ORDER=",
    runnerSettings["RequestedStepDivisionOrder"],
    " overridden by classic coupled segmentation; using FT_DIVISION_ORDER=",
    divisionOrder, " for both placement and +/-1/k matching"]];
levelEpsilonHalos = runnerSettings["LevelEpsilonHalos"];
levelEpsilonHalo[level_Integer] := If[1 <= level <= Length[levelEpsilonHalos],
  levelEpsilonHalos[[level]], 0];
ft2UserRawFloor[epsilonOrder_Integer, halos_List,
    level_Integer] := If[level === 0, epsilonOrder,
  Max[epsilonOrder + level +
    If[1 <= level <= Length[halos], halos[[level]], 0], 1]];
requestedEpsilonOrder[level_Integer] := Max[
  epsOrder + level + boundaryExtraOrder + levelEpsilonHalo[level], 1];
nativeRequiredRawTop[lowerLevel_Integer] :=
  ft2UserRawFloor[epsOrder, levelEpsilonHalos, lowerLevel];

(* Old FeynmanTrick/DiffExp prescribed both endpoints and every matrix/IBP
   segmentation factor with a consistent +i delta side.  DiffExp2's config
   is reset at every level, so rebuild that effective list explicitly here;
   otherwise Transport receives an empty Prescriptions record and loses the
   branch sheet after crossing an interior singularity. *)
levelDeltaPrescriptionSign[level_Integer] := If[
  ListQ[levelDeltaPrescriptionSigns] &&
    1 <= level <= Length[levelDeltaPrescriptionSigns],
  levelDeltaPrescriptionSigns[[level]],
  deltaPrescriptionSign
];

levelDeltaPrescriptionSignsForCount[nLevels_Integer] := Which[
  configuredLevelDeltaPrescriptionSigns === Automatic,
    ConstantArray[deltaPrescriptionSign, nLevels],
  ListQ[configuredLevelDeltaPrescriptionSigns] &&
      Length[configuredLevelDeltaPrescriptionSigns] === nLevels &&
      AllTrue[configuredLevelDeltaPrescriptionSigns,
        MemberQ[{-1, 1}, #] &],
    configuredLevelDeltaPrescriptionSigns,
  True,
    $Failed
];

levelDeltaPrescriptions[level_Integer, var_Symbol, sys_Association,
    extra_List] := Module[
  {raw, factors, projectedExtra, sign, endpointEquivalentQ},
  sign = levelDeltaPrescriptionSign[level];
  projectedExtra = DiffExp2`Transport`EpsilonZeroSingularFactors[extra, var];
  raw = Join[Lookup[sys, "SingularFactors", {}], projectedExtra];
  factors = Flatten[Map[Module[{fl = FactorList[Factor[Numerator[Together[#]]]]},
      First /@ Select[fl, !FreeQ[First[#], var] &]] &, raw]];
  endpointEquivalentQ[factor_] :=
    AnyTrue[{var, 1 - var},
      TrueQ[PossibleZeroQ[Expand[factor - #]]] ||
        TrueQ[PossibleZeroQ[Expand[factor + #]]] &];
  factors = Select[factors, !endpointEquivalentQ[#] &];
  factors = DeleteDuplicates[factors,
    TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
      TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &];
  (* Keep the endpoint polynomials in their physical orientation.  In
     particular, factoring 1-var produces the canonical factor var-1 and a
     numeric unit -1; discarding that unit would put both endpoints on the
     same local rim.  Config's side-aware canonicalization must see 1-var so
     it can turn {1-var,sign} into {var-1,-sign}. *)
  Join[{{var, sign}, {1 - var, sign}}, {#, sign} & /@ factors]];

inputPrecision = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
prepCacheRoot = runnerSettings["PrepCacheRoot"];
forcePrepRebuild = runnerSettings["ForcePrepRebuild"];
migrateLegacyPrep = runnerSettings["MigrateLegacyPreparation"];
resumeLadderFile = runnerSettings["ResumeCheckpoint"];
ladderCheckpointDir = runnerSettings["CheckpointDirectory"];
allowStaleLadderCheckpoint = runnerSettings["AllowStaleCheckpoint"];
saveNativeTransportCheckpoint =
  runnerSettings["SaveNativeTransportCheckpoint"];
matchingHaloProfileEnabled =
  envOrDefault["FT_DISABLE_MATCHING_HALO_PROFILE", "0"] =!= "1";
matchingHaloProfileRoot = FileNameJoin[
  {prepCacheRoot, "MatchingHaloProfiles"}];

(* RunFullIteration is FIRE-dominated but independent of DiffExp2's
   transport settings.  Persist both the populated ftData and FIRE's
   in-memory reduction cache so a fresh Wolfram process can resume at the
   transport ladder.  Preparation has a deliberately narrow, exact contract:
   topology, combination sequence, dimension, preparation-affecting options,
   FIRE/Wolfram runtime, and only the three modules that construct levels,
   FIRE problems, bases, reductions, and differential matrices.  Paths in the
   source manifest are repository-relative; runner, facade, boundary,
   transport, and export-only edits therefore cannot invalidate expensive
   FIRE data.  The complete contract is stored beside the digest and compared
   exactly on load, so the digest is never the sole stale-data guard.

   FT_REBUILD_PREP=1 forces a rebuild.  Legacy snapshots have no exact
   preparation-source manifest and are rejected by default.  The explicit
   FT_MIGRATE_LEGACY_PREP=1 escape hatch accepts exactly one candidate only
   after exact input/level/reduction-key validation and rewrites it immediately
   as v3; ambiguous candidates remain rejected.  V2 also proves its retained
   FIRE setup record.  V1 cannot: its old tuple keys are uniquely bound to the
   retained level topologies and re-keyed with the current exact fallback
   setup record, while provenance remains explicitly marked unverified. *)
$ftPrepCacheVersion = 3;
$ftPrepContractSchema = "FeynmanTrick.PreparedFTContract/v3";
$ftPrepPreparationSourcePaths = {
  "FeynmanTrick/PropagatorAlgebra.m",
  "FeynmanTrick/FIREInterface.m",
  "FeynmanTrick/FeynmanTrickIteration.m"
};

(* Why this manifest is intentionally only three files:
   - FeynmanTrick.m supplies configuration plumbing, but every value that can
     affect preparation is evaluated into PreparationConfiguration or
     DimensionExpression below.
   - MatrixExport.m only serializes diagnostic artifacts after a level has
     already been computed; it does not mutate the retained FTData.
   - LevelReduction.m classifies transport-time boundary requests.  It does
     not construct FTData, and every request must still find its exact,
     setup-fingerprinted reduction key before a snapshot can load.
   - BoundaryConditions.m, DiffExp2Pipeline.m, and this runner act only after
     preparation.  Their broader provenance belongs to ladder checkpoints,
     not the FIRE preparation cache. *)

ftPrepRelativeSourceIdentity[relativePath_String] := Module[{path},
  path = FileNameJoin[Prepend[FileNameSplit[relativePath], repoRoot]];
  If[!FileExistsQ[path], Return[Failure["FeynmanTrickPreparationSource", <|
    "Detail" -> "preparation source is missing",
    "RelativePath" -> relativePath|>], Module]];
  <|"RelativePath" -> relativePath,
    "SHA256" -> IntegerString[FileHash[path, "SHA256"], 16, 64]|>];

$ftPrepPreparationSourceIdentities =
  ftPrepRelativeSourceIdentity /@ $ftPrepPreparationSourcePaths;

$ftPrepCustomFamilySourcePaths = {
  "FeynmanTrick/FamilySpec.m",
  "FeynmanTrick/PipelineRequest.m"
};
$ftPrepCustomFamilySourceIdentities =
  ftPrepRelativeSourceIdentity /@ $ftPrepCustomFamilySourcePaths;
$ft2AllDiscoveryRunnerSourceIdentity =
  ftPrepRelativeSourceIdentity["Scripts/run_ft_stepwise2.m"];

ftPrepSelectPreparationSourceIdentities[sourceMap_Association] :=
  Lookup[sourceMap, $ftPrepPreparationSourcePaths,
    Missing["PreparationSourceAbsent"]];

ftPrepConfigurationRecord[] := Module[{cfg = FeynmanTrick`Private`$FTConfig},
  <|
    "DimensionVariable" -> Lookup[cfg, "DimensionVariable", Missing["Unset"]],
    "EpsilonSymbol" -> Lookup[cfg, "EpsilonSymbol", Missing["Unset"]],
    "FixedParameterValue" ->
      Lookup[cfg, "FixedParameterValue", Missing["Unset"]],
    "FixedParameterValues" ->
      Lookup[cfg, "FixedParameterValues", Automatic],
    "AutoDetectRestrictions" ->
      Lookup[cfg, "AutoDetectRestrictions", Missing["Unset"]]
  |>];

ftPrepContractRecord[name_String, topology_Association, sequence_List] := <|
  "Schema" -> $ftPrepContractSchema,
  "CacheVersion" -> $ftPrepCacheVersion,
  "Example" -> name,
  "Topology" -> topology,
  "CombinationSequence" -> sequence,
  "NumericalPoint" -> {},
  "DimensionExpression" -> FeynmanTrick`Private`DimensionExpression[],
  "PreparationConfiguration" -> ftPrepConfigurationRecord[],
  "WolframRuntime" -> <|"Version" -> $Version, "SystemID" -> $SystemID|>,
  "FIRERuntime" ->
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[],
  "PreparationSources" -> $ftPrepPreparationSourceIdentities
|>;

ftPrepCustomContractRecord[name_String, request_Association,
    resolution_:None] := Module[
  {family = request["Family"], contract, mode, selectionFields},
  contract = ftPrepContractRecord[
    name, family["Topology"], family["CombinationSequence"]];
  mode = request["OutputMode"];
  selectionFields = Switch[mode,
    "Explicit",
      If[resolution =!= None,
        Return[Failure["FeynmanTrickPreparationContract", <|
          "Detail" -> "explicit output selection cannot carry an All-master resolution"|>],
          Module]];
      <|
        "OutputSelectionMode" -> "Explicit",
        "OutputSelection" -> family["OutputIntegrals"]
      |>,
    "AllPendingDiscovery",
      If[!TrueQ[
          FeynmanTrick`PipelineRequest`ResolvedAllOutputSelectionQ[
            request, resolution]],
        Return[Failure["FeynmanTrickPreparationContract", <|
          "Detail" -> "All output selection requires a valid resolved master contract"|>],
          Module]];
      <|
        (* Keep the user's original All selection explicit while binding the
           discovered ordered basis into every derived preparation identity. *)
        "OutputSelectionMode" -> "AllResolved",
        "OutputSelection" -> All,
        "AllSelectionRequestID" -> resolution["SelectionRequestID"],
        "OutputResolutionID" -> resolution["ResolutionID"],
        "ResolvedOutputSelection" -> resolution["Masters"],
        "ResolvedOutputRequests" -> resolution["OutputRequests"]
      |>,
    _,
      Return[Failure["FeynmanTrickPreparationContract", <|
        "Detail" -> "unsupported custom-family output mode",
        "OutputMode" -> mode|>], Module]
  ];
  Join[contract, <|
    "CustomFamilyPreparationSchema" ->
      "FeynmanTrick.CustomFamilyPreparation/v1",
    "FamilyID" -> family["FamilyID"],
    "PipelineRequestID" -> request["RequestID"],
    "NumericalPoint" -> family["NumericalPoint"],
    "DimensionExpression" -> family["Dimension"],
    "CustomFamilySources" -> $ftPrepCustomFamilySourceIdentities
  |>, selectionFields]
];

$ft2AllDiscoveryCacheSchema = "FeynmanTrick.AllMasterDiscoveryCache/v2";
$ft2AllDiscoveryContractSchema = "FeynmanTrick.AllMasterDiscoveryContract/v2";

ft2AllDiscoveryTopology[request_Association] := Module[
  {family = request["Family"], topology, propagators, replacements},
  topology = family["Topology"];
  propagators =
    FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
      topology["Propagators"], family["NumericalPoint"]];
  replacements =
    FeynmanTrick`FIREInterface`Private`applyFIRENumericalPointToReplacements[
      topology["Replacements"], family["NumericalPoint"]];
  If[propagators === $Failed || replacements === $Failed,
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "NumericalPoint substitutions are cyclic or do not reach a fixed point"|>],
      Module]];
  topology = Join[topology, <|
    "Name" -> family["Name"] <> "_L0_all_" <>
      StringTake[request["RequestID"], -12],
    "Propagators" -> propagators,
    "Replacements" -> replacements
  |>];
  topology
];

ft2AllDiscoveryContractRecord[name_String, request_Association] := Module[
  {family = request["Family"], fireSource, discoveryTopology},
  discoveryTopology = ft2AllDiscoveryTopology[request];
  If[FailureQ[discoveryTopology], Return[discoveryTopology, Module]];
  fireSource = SelectFirst[$ftPrepPreparationSourceIdentities,
    Lookup[#, "RelativePath", None] === "FeynmanTrick/FIREInterface.m" &,
    Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "FIREInterface source identity is absent"|>]];
  <|
    "Schema" -> $ft2AllDiscoveryContractSchema,
    "Example" -> name,
    "PipelineRequestID" -> request["RequestID"],
    "SelectionRequestID" -> request["OutputRequests"][[1, "RequestID"]],
    "FamilyID" -> family["FamilyID"],
    "DiscoveryTopology" -> discoveryTopology,
    "CombinationSequence" -> family["CombinationSequence"],
    "NumericalPoint" -> family["NumericalPoint"],
    "DimensionExpression" -> family["Dimension"],
    "PreparationConfiguration" -> ftPrepConfigurationRecord[],
    "WolframRuntime" -> <|"Version" -> $Version, "SystemID" -> $SystemID|>,
    "FIRERuntime" ->
      FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[],
    "DiscoverySources" -> Join[
      {fireSource, $ft2AllDiscoveryRunnerSourceIdentity},
      $ftPrepCustomFamilySourceIdentities]
  |>
];

ft2AllDiscoveryContractIdentity[contract_Association] :=
  "ft-all-discovery-contract-" <>
    IntegerString[Hash[contract, "SHA256"], 16, 64];

ft2AllDiscoveryCacheFile[name_String, contract_Association] :=
  FileNameJoin[{prepCacheRoot, name <> "_all_resolution_" <>
    StringTake[ft2AllDiscoveryContractIdentity[contract], -32] <> ".wxf"}];

ft2AllDiscoveryCachePayloadQ[payload_, request_Association,
    contract_Association] := Module[{resolution},
  If[!AssociationQ[payload] || Sort[Keys[payload]] =!= Sort[{
      "Schema", "Contract", "ContractIdentity", "Resolution"}] ||
      Lookup[payload, "Schema", None] =!= $ft2AllDiscoveryCacheSchema ||
      Lookup[payload, "Contract", None] =!= contract ||
      Lookup[payload, "ContractIdentity", None] =!=
        ft2AllDiscoveryContractIdentity[contract],
    Return[False, Module]];
  resolution = Lookup[payload, "Resolution", None];
  TrueQ[FeynmanTrick`PipelineRequest`ResolvedAllOutputSelectionQ[
    request, resolution]]
];

ft2LoadAllDiscoveryCache[file_String, request_Association,
    contract_Association] := Module[{payload},
  If[!FileExistsQ[file], Return[$Failed, Module]];
  payload = Quiet[Check[Import[file, "WXF"], $Failed]];
  If[ft2AllDiscoveryCachePayloadQ[payload, request, contract],
    payload["Resolution"], $Failed]
];

ft2SaveAllDiscoveryCache[file_String, request_Association,
    contract_Association, resolution_Association] := Module[
  {payload, directory, tmp, wrote, loaded},
  payload = <|
    "Schema" -> $ft2AllDiscoveryCacheSchema,
    "Contract" -> contract,
    "ContractIdentity" -> ft2AllDiscoveryContractIdentity[contract],
    "Resolution" -> resolution
  |>;
  If[!ft2AllDiscoveryCachePayloadQ[payload, request, contract],
    Return[$Failed, Module]];
  directory = DirectoryName[file];
  If[!DirectoryQ[directory],
    Quiet[Check[
      CreateDirectory[directory, CreateIntermediateDirectories -> True],
      Return[$Failed, Module]]]];
  tmp = file <> ".tmp-" <> ToString[$ProcessID];
  If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
  wrote = Quiet[Check[Export[tmp, payload, "WXF"], $Failed]];
  If[wrote === $Failed || !FileExistsQ[tmp], Return[$Failed, Module]];
  loaded = Quiet[Check[Import[tmp, "WXF"], $Failed]];
  If[!ft2AllDiscoveryCachePayloadQ[loaded, request, contract],
    Quiet[DeleteFile[tmp]];
    Return[$Failed, Module]];
  If[!Quiet[Check[
      RenameFile[tmp, file, OverwriteTarget -> True]; True, False]],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    Return[$Failed, Module]];
  file
];

(* Overridable seams let focused tests exercise every discovery guard without
   starting FIRE.  Production calls the package functions directly. *)
ft2AllSetupFIRE[topology_Association] :=
  FeynmanTrick`FIREInterface`SetupFIRE[topology];
ft2AllFindBasis[topology_Association] :=
  FeynmanTrick`FIREInterface`FindBasis[topology];

ft2ValidateAllDiscoverySetup[request_Association, setupTopology_] := Module[
  {family = request["Family"], expectedTopology, expectedN, actualN,
   numerators, expectedNumerators},
  If[!AssociationQ[setupTopology],
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "SetupFIRE failed during L0 All-master discovery"|>], Module]];
  expectedN = family["NumPropagators"];
  expectedTopology = ft2AllDiscoveryTopology[request];
  If[FailureQ[expectedTopology], Return[expectedTopology, Module]];
  actualN = Lookup[setupTopology, "NumPropagators", None];
  numerators = Lookup[setupTopology, "NumeratorPositions", {}];
  expectedNumerators = family["NumeratorPositions"];
  If[actualN =!= expectedN || numerators =!= expectedNumerators ||
      Lookup[setupTopology, "OriginalNumPropagators", expectedN] =!= expectedN,
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "FIRE changed the declared irreducible-numerator slots",
      "ExpectedArity" -> expectedN,
      "ActualArity" -> actualN,
      "ExpectedNumeratorPositions" -> expectedNumerators,
      "NumeratorPositions" -> numerators|>], Module]];
  If[!TrueQ[Lookup[setupTopology, "StartFileReady", False]] ||
      Lookup[setupTopology, "Name", None] =!= expectedTopology["Name"] ||
      Lookup[setupTopology, "LoopMomenta", None] =!=
        expectedTopology["LoopMomenta"] ||
      Lookup[setupTopology, "ExternalMomenta", None] =!=
        expectedTopology["ExternalMomenta"] ||
      Lookup[setupTopology, "OriginalPropagators", None] =!=
        expectedTopology["Propagators"] ||
      Lookup[setupTopology, "Replacements", None] =!=
        expectedTopology["Replacements"] ||
      Lookup[setupTopology, "EliminatedPositions", {}] =!=
        expectedTopology["EliminatedPositions"],
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "SetupFIRE returned stale or mathematically different L0 topology metadata"|>],
      Module]];
  setupTopology
];

ft2ValidateAllDiscoveryMasters[request_Association, masters_] :=
  FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
    request, masters];

ft2DiscoverAllOutputSelection[request_Association] := Module[
  {topology, setupTopology, basisTopology, validatedSetup, resolution},
  topology = ft2AllDiscoveryTopology[request];
  If[FailureQ[topology], Return[topology, Module]];
  setupTopology = Quiet[Check[ft2AllSetupFIRE[topology], $Failed]];
  validatedSetup = ft2ValidateAllDiscoverySetup[request, setupTopology];
  If[FailureQ[validatedSetup], Return[validatedSetup, Module]];
  basisTopology = Quiet[Check[ft2AllFindBasis[validatedSetup], $Failed]];
  If[!AssociationQ[basisTopology] ||
      KeyDrop[basisTopology, "Masters"] =!=
        KeyDrop[validatedSetup, "Masters"],
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "FindBasis failed or returned a different L0 setup during All-master discovery"|>],
      Module]];
  resolution = ft2ValidateAllDiscoveryMasters[
    request, Lookup[basisTopology, "Masters", None]];
  resolution
];

ft2ResolveAllOutputSelection[name_String, request_Association] := Module[
  {contract, file, cached, resolution, saved},
  contract = ft2AllDiscoveryContractRecord[name, request];
  If[FailureQ[contract], Return[contract, Module]];
  If[AnyTrue[Lookup[contract, "DiscoverySources", {}], FailureQ],
    Return[Failure["FeynmanTrickAllMasterDiscovery", <|
      "Detail" -> "All-master discovery source identity is incomplete",
      "Contract" -> contract|>], Module]];
  file = ft2AllDiscoveryCacheFile[name, contract];
  cached = If[forcePrepRebuild, $Failed,
    ft2LoadAllDiscoveryCache[file, request, contract]];
  If[AssociationQ[cached],
    Print["FTPREP ALL DISCOVERY CACHE HIT ", file];
    Return[cached, Module]];
  Print["FTPREP ALL DISCOVERY CACHE MISS ", file];
  resolution = ft2DiscoverAllOutputSelection[request];
  If[FailureQ[resolution], Return[resolution, Module]];
  saved = ft2SaveAllDiscoveryCache[file, request, contract, resolution];
  If[saved === $Failed,
    Print["FTPREP ALL DISCOVERY CACHE WRITE FAILED ", file],
    Print["FTPREP ALL DISCOVERY CACHE WRITE ", file]];
  resolution
];

ftPrepContractKey[contract_Association] := Hash[contract, "SHA256"];

$ftPrepProvenanceSchema = "FeynmanTrick.PreparedFTProvenance/v1";

ftPrepFreshProvenance[] := <|
  "Schema" -> $ftPrepProvenanceSchema,
  "Kind" -> "PreparedV3",
  "SourceSnapshotVersion" -> $ftPrepCacheVersion,
  "VerifiedFIRESetupProvenance" -> True,
  "VerifiedPreparationSourceProvenance" -> True
|>;

ftPrepProvenanceQ[provenance_] := AssociationQ[provenance] &&
  Lookup[provenance, "Schema", None] === $ftPrepProvenanceSchema &&
  Switch[Lookup[provenance, "Kind", None],
    "PreparedV3",
      Lookup[provenance, "SourceSnapshotVersion", None] === 3 &&
        TrueQ[Lookup[provenance, "VerifiedFIRESetupProvenance", False]] &&
        TrueQ[Lookup[provenance,
          "VerifiedPreparationSourceProvenance", False]],
    "LegacyV2Validated",
      Lookup[provenance, "SourceSnapshotVersion", None] === 2 &&
        TrueQ[Lookup[provenance, "VerifiedFIRESetupProvenance", False]] &&
        TrueQ[Lookup[provenance,
          "VerifiedPreparationSourceProvenance", True] === False],
    "LegacyV1Rekeyed",
      Lookup[provenance, "SourceSnapshotVersion", None] === 1 &&
        TrueQ[Lookup[provenance,
          "VerifiedFIRESetupProvenance", True] === False] &&
        TrueQ[Lookup[provenance,
          "VerifiedPreparationSourceProvenance", True] === False] &&
        Lookup[provenance, "ReductionKeyMigration", None] ===
          "LegacyTupleToFallbackSetup/v1" &&
        IntegerQ[Lookup[provenance, "RekeyedEntries", None]] &&
        Lookup[provenance, "RekeyedEntries", -1] >= 0 &&
        StringQ[Lookup[provenance, "SourceSnapshotSHA256", None]],
    _, False];

$ftLadderCheckpointVersion = 2;
$ft2ActivePipelineRequestID = None;
$ft2ActiveFamilyID = None;
$ft2ActiveAllSelectionRequestID = None;
$ft2ActiveOutputResolutionID = None;
$ft2CheckpointRequestMetadataKeys = {
  "PipelineRequestID", "FamilyID", "AllSelectionRequestID",
  "OutputResolutionID"
};

ft2CheckpointRequestMetadata[] :=
  If[StringQ[$ft2ActivePipelineRequestID] && StringQ[$ft2ActiveFamilyID],
    Join[<|"PipelineRequestID" -> $ft2ActivePipelineRequestID,
      "FamilyID" -> $ft2ActiveFamilyID|>,
      If[StringQ[$ft2ActiveAllSelectionRequestID] &&
          StringQ[$ft2ActiveOutputResolutionID],
        <|
          "AllSelectionRequestID" -> $ft2ActiveAllSelectionRequestID,
          "OutputResolutionID" -> $ft2ActiveOutputResolutionID
        |>, <||>]],
    <||>];

$ftLadderSourcePaths = Sort[DeleteDuplicates[Join[
    FileNames["*.m", FileNameJoin[{repoRoot, "DiffExp2"}], Infinity],
    FileNames["*.m", FileNameJoin[{repoRoot, "FeynmanTrick"}], Infinity],
    Select[FileNames["*", FileNameJoin[{repoRoot, "cpp"}], Infinity],
      FileType[#] === File &],
    {FileNameJoin[{repoRoot, "CMakeLists.txt"}]},
    {FileNameJoin[{repoRoot, "FeynmanTrick.m"}]},
    {ExpandFileName[$InputFileName]}]]];
$ftLadderSourceFingerprint = Hash[
  ({#, FileHash[#, "SHA256"]} & /@ $ftLadderSourcePaths), "SHA256"];

preparedFTDataQ[data_] := AssociationQ[data] &&
  IntegerQ[Lookup[data, "NumLevels", None]] &&
  AssociationQ[Lookup[data, "Levels", None]] &&
  AllTrue[Range[data["NumLevels"]], Function[level,
    Module[{ld = Lookup[data["Levels"], level, <||>], masters, mat},
      masters = Lookup[ld, "Masters", {}]; mat = Lookup[ld, "DiffMatrix", {}];
      TrueQ[Lookup[ld, "Computed", False]] && masters =!= {} &&
        MatrixQ[mat] && Dimensions[mat] === {Length[masters], Length[masters]}]]];

ftPrepRuntimeRecordCompatibleQ[stored_, expected_Association] :=
  AssociationQ[stored] &&
    KeyTake[stored, Keys[expected]] === expected;

preparedFTDataMatchesContractQ[data_, contract_Association,
    provenance_:Automatic] := Module[
  {nLevels, levels, config, runtime, provenanceRecord, legacyV1Q,
   customQ, levelZero, outputMode, customOutputMatchQ},
  If[!preparedFTDataQ[data] ||
      Lookup[contract, "Schema", None] =!= $ftPrepContractSchema,
    Return[False, Module]];
  provenanceRecord = Replace[provenance, Automatic :> ftPrepFreshProvenance[]];
  If[!ftPrepProvenanceQ[provenanceRecord], Return[False, Module]];
  legacyV1Q = Lookup[provenanceRecord, "Kind", None] ===
    "LegacyV1Rekeyed";
  nLevels = Lookup[data, "NumLevels", None];
  levels = Lookup[data, "Levels", <||>];
  config = Lookup[contract, "PreparationConfiguration", <||>];
  runtime = Lookup[contract, "FIRERuntime", <||>];
  customQ = KeyExistsQ[contract, "PipelineRequestID"];
  levelZero = Lookup[levels, 0, <||>];
  outputMode = Lookup[contract, "OutputSelectionMode", None];
  customOutputMatchQ = If[!customQ, True, Switch[outputMode,
    "Explicit",
      Lookup[levelZero, "Masters", Missing["Absent"]] ===
        Lookup[contract, "OutputSelection", Missing["ContractAbsent"]],
    "AllResolved",
      Lookup[contract, "OutputSelection", Missing["ContractAbsent"]] === All &&
      Lookup[levelZero, "Masters", Missing["Absent"]] ===
        Lookup[contract, "ResolvedOutputSelection",
          Missing["ContractAbsent"]] &&
      Lookup[data, "AllSelectionRequestID", Missing["Absent"]] ===
        Lookup[contract, "AllSelectionRequestID",
          Missing["ContractAbsent"]] &&
      Lookup[data, "OutputResolutionID", Missing["Absent"]] ===
        Lookup[contract, "OutputResolutionID",
          Missing["ContractAbsent"]] &&
      Lookup[data, "ResolvedOutputRequests", Missing["Absent"]] ===
        Lookup[contract, "ResolvedOutputRequests",
          Missing["ContractAbsent"]],
    _, False]];
  TrueQ[
    Lookup[data, "TopTopology", Missing["Absent"]] ===
      Lookup[contract, "Topology", Missing["Absent"]] &&
    Lookup[data, "CombinationSequence", Missing["Absent"]] ===
      Lookup[contract, "CombinationSequence", Missing["Absent"]] &&
    Lookup[data, "NumericalPoint", Missing["Absent"]] ===
      Lookup[contract, "NumericalPoint", Missing["Absent"]] &&
    Lookup[data, "FixedParamValue", Missing["Absent"]] ===
      Lookup[config, "FixedParameterValue", Missing["Unset"]] &&
    Lookup[data, "FixedParamValues",
      ConstantArray[
        Lookup[data, "FixedParamValue", Missing["Absent"]],
        nLevels]] ===
      Replace[
        Lookup[config, "FixedParameterValues", Automatic],
        Automatic :> ConstantArray[
          Lookup[config, "FixedParameterValue", Missing["Unset"]],
          nLevels]] &&
    (!customQ || (
      Lookup[contract, "CustomFamilyPreparationSchema", None] ===
        "FeynmanTrick.CustomFamilyPreparation/v1" &&
      Lookup[data, "FamilyID", Missing["Absent"]] ===
        Lookup[contract, "FamilyID", Missing["ContractAbsent"]] &&
      Lookup[data, "PipelineRequestID", Missing["Absent"]] ===
        Lookup[contract, "PipelineRequestID", Missing["ContractAbsent"]] &&
      TrueQ[customOutputMatchQ] &&
      Lookup[Lookup[data, "TopTopology", <||>],
        "Dimension", Missing["Absent"]] ===
        Lookup[contract, "DimensionExpression", Missing["ContractAbsent"]])) &&
    nLevels === Length[Lookup[contract, "CombinationSequence", {}]] &&
    AssociationQ[runtime] && runtime =!= <||> &&
    AllTrue[Range[nLevels], Function[level,
      Module[{topology, setup},
        topology = Lookup[levels[level], "Topology", <||>];
        setup = Lookup[topology, "SetupFingerprintRecord", None];
        If[legacyV1Q,
          (* Historical v1 topologies contain no SetupFingerprintRecord.
             Never fabricate one: exact reductions are re-keyed below with
             FIREInterface's full fallback topology/config record, and any
             later cache miss remains unable to launch FIRE. *)
          !AssociationQ[setup] &&
            Lookup[
              FeynmanTrick`FIREInterface`Private`fallbackSetupFingerprintRecord[
                topology],
              {"AutoDetectRestrictions", "DimensionVariable"},
              Missing["Unset"]] ===
            Lookup[config,
              {"AutoDetectRestrictions", "DimensionVariable"},
              Missing["Unset"]],
          ftPrepRuntimeRecordCompatibleQ[setup, runtime] &&
            Lookup[setup, "AutoDetectRestrictions", Missing["Unset"]] ===
              Lookup[config, "AutoDetectRestrictions", Missing["Unset"]]
        ]
      ]]]
  ]
];

requiredReductionKeys[data_] := Flatten[Table[
  Module[{ld = data["Levels"][level], below = data["Levels"][level - 1],
      topo, reqs},
    topo = ld["Topology"];
    reqs = FeynmanTrick`LevelReduction`BoundaryRequestRecords[
      below["Masters"], ld["CombinedPositions"]];
    Map[Function[req,
      FeynmanTrick`FIREInterface`Private`reductionCacheKey[topo,
        FeynmanTrick`FIREInterface`Private`normalizeIntegralIndex[
          topo, req["NeededVec"]]]], reqs]],
  {level, 1, data["NumLevels"]}], 1];

preparedReductionCacheQ[data_, rc_] := AssociationQ[rc] &&
  AllTrue[requiredReductionKeys[data], Function[key,
    Module[{entry = Lookup[rc, Key[key], Missing["NotCached"]]},
      AssociationQ[entry] && KeyExistsQ[entry, "Reduction"] &&
        KeyExistsQ[entry, "Masters"]]]];

hydratePreparedReductionCache[data_Association] := Module[
  {nLevels = Lookup[data, "NumLevels", 0], batches},
  If[!IntegerQ[nLevels] || nLevels < 0,
    Return[Failure["PreparedFTReductionHydration", <|
      "Detail" -> "prepared data has an invalid level count",
      "NumLevels" -> nLevels|>], Module]];
  batches = AssociationMap[
    FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[data, #] &,
    Range[nLevels]];
  If[AnyTrue[Values[batches], !AssociationQ[#] &],
    Return[Failure["PreparedFTReductionHydration", <|
      "Detail" -> "an exact boundary-reduction batch failed",
      "FailedLevels" -> Keys@Select[batches, !AssociationQ[#] &]|>],
      Module]];
  If[!preparedReductionCacheQ[data,
      FeynmanTrick`FIREInterface`Private`$ReductionCache],
    Return[Failure["PreparedFTReductionHydration", <|
      "Detail" ->
        "exact boundary batches did not populate every snapshot reduction key"
      |>], Module]];
  batches
];

ftPrepKey[name_String, topology_Association, sequence_List] :=
  ftPrepContractKey[ftPrepContractRecord[name, topology, sequence]];

ftPrepFile[name_, key_] := FileNameJoin[{prepCacheRoot,
  name <> "_" <> IntegerString[Abs[key], 16] <> ".mx"}];

savePreparedFT[file_String, contract_Association, data_,
    provenance_:Automatic] := Module[
  {payload, tmp, ok, key = ftPrepContractKey[contract], provenanceRecord},
  provenanceRecord = Replace[provenance, Automatic :> ftPrepFreshProvenance[]];
  If[!ftPrepProvenanceQ[provenanceRecord] ||
      !preparedFTDataMatchesContractQ[data, contract, provenanceRecord],
    Print["FTPREP CACHE CONTRACT MISMATCH; not saving ", file];
    Return[$Failed, Module]];
  If[!preparedReductionCacheQ[data,
      FeynmanTrick`FIREInterface`Private`$ReductionCache],
    Print["FTPREP CACHE INCOMPLETE; not saving ", file];
    Return[$Failed, Module]];
  If[!DirectoryQ[DirectoryName[file]],
    CreateDirectory[DirectoryName[file], CreateIntermediateDirectories -> True]];
  payload = <|
    "Version" -> $ftPrepCacheVersion, "Key" -> key,
    "Contract" -> contract, "Provenance" -> provenanceRecord,
    "FTData" -> data,
    "ReductionCache" -> KeyTake[
      FeynmanTrick`FIREInterface`Private`$ReductionCache,
      requiredReductionKeys[data]]|>;
  Global`$FT2PreparedSnapshot = payload;
  tmp = file <> ".tmp-" <> ToString[$ProcessID] <> ".mx";
  If[FileExistsQ[tmp], DeleteFile[tmp]];
  ok = Quiet[Check[DumpSave[tmp, Global`$FT2PreparedSnapshot]; True, False]];
  Clear[Global`$FT2PreparedSnapshot];
  If[!TrueQ[ok], Print["FTPREP CACHE WRITE FAILED ", file];
    Return[$Failed, Module]];
  If[FileExistsQ[file], DeleteFile[file]];
  Quiet[Check[RenameFile[tmp, file],
    Print["FTPREP CACHE RENAME FAILED ", file]; Return[$Failed, Module]]];
  Print["FTPREP CACHE WRITE ", file];
  file];

readPreparedFTPayload[file_String] := Module[{payload, ok},
  If[!FileExistsQ[file], Return[$Failed, Module]];
  Clear[Global`$FT2PreparedSnapshot];
  ok = Quiet[Check[Get[file]; True, False]];
  If[!TrueQ[ok], Clear[Global`$FT2PreparedSnapshot];
    Return[$Failed, Module]];
  payload = Global`$FT2PreparedSnapshot;
  Clear[Global`$FT2PreparedSnapshot];
  If[AssociationQ[payload], payload, $Failed]];

loadPreparedFT[file_String, contract_Association] := Module[
  {payload, data, reductionCache, provenance,
   key = ftPrepContractKey[contract], reason},
  If[!FileExistsQ[file], Return[$Failed, Module]];
  payload = readPreparedFTPayload[file];
  If[!AssociationQ[payload],
    Print["FTPREP CACHE REJECT ", file, ": unreadable payload"];
    Return[$Failed, Module]];
  data = Lookup[payload, "FTData", None];
  reductionCache = Lookup[payload, "ReductionCache", <||>];
  provenance = Lookup[payload, "Provenance", None];
  reason = Which[
    Lookup[payload, "Version", None] =!= $ftPrepCacheVersion,
      "unsupported snapshot version",
    Lookup[payload, "Contract", None] =!= contract,
      "preparation contract differs",
    Lookup[payload, "Key", None] =!= key,
      "contract digest differs",
    !ftPrepProvenanceQ[provenance],
      "preparation provenance is absent or malformed",
    !preparedFTDataMatchesContractQ[data, contract, provenance],
      "prepared level data does not match the exact contract",
    !preparedReductionCacheQ[data, reductionCache],
      "required exact reduction keys are absent or malformed",
    True, None];
  If[StringQ[reason],
    Print["FTPREP CACHE REJECT ", file, ": ", reason];
    Return[$Failed, Module]];
  FeynmanTrick`FIREInterface`Private`$ReductionCache =
    Join[FeynmanTrick`FIREInterface`Private`$ReductionCache,
      reductionCache];
  Print["FTPREP CACHE HIT ", file];
  data];

legacyV1TopologyIdentity[topology_Association] := Lookup[topology,
  {"WorkDirectory", "Name", "ProblemNumber", "NumPropagators"},
  Missing["Absent"]];

legacyV1RekeyReductionCache[data_Association, oldCache_Association] := Module[
  {levels, topologies, topologyIdentities, rekeyed = <||>, result,
   oldKey, entry, matches, topology, integral, newKey, existing},
  levels = Lookup[data, "Levels", <||>];
  topologies = Lookup[levels[#], "Topology", <||>] & /@
    Range[Lookup[data, "NumLevels", 0]];
  If[!AllTrue[topologies, AssociationQ],
    Return[Failure["LegacyPreparationMigration", <|
      "Detail" -> "v1 prepared levels do not all contain topologies"|>],
      Module]];
  topologyIdentities = legacyV1TopologyIdentity /@ topologies;
  If[!DuplicateFreeQ[topologyIdentities],
    Return[Failure["LegacyPreparationMigration", <|
      "Detail" -> "v1 level topology identities are ambiguous",
      "Identities" -> topologyIdentities|>], Module]];
  result = Catch[
    Scan[Function[rule,
      oldKey = First[rule]; entry = Last[rule];
      If[!MatchQ[oldKey, {_, _String, _Integer, _Integer, _List}] ||
          !AssociationQ[entry] ||
          !KeyExistsQ[entry, "Reduction"] ||
          !ListQ[Lookup[entry, "Masters", None]],
        Throw[Failure["LegacyPreparationMigration", <|
          "Detail" -> "v1 reduction entry or tuple key is malformed",
          "Key" -> oldKey|>], "LegacyV1Rekey"]];
      matches = Pick[topologies,
        SameQ[#, Take[oldKey, 4]] & /@ topologyIdentities, True];
      If[Length[matches] =!= 1,
        Throw[Failure["LegacyPreparationMigration", <|
          "Detail" -> "v1 reduction key does not identify exactly one level",
          "Key" -> oldKey, "MatchingLevels" -> Length[matches]|>],
          "LegacyV1Rekey"]];
      topology = First[matches];
      integral =
        FeynmanTrick`FIREInterface`Private`normalizeIntegralIndex[
          topology, Last[oldKey]];
      If[integral === $Failed,
        Throw[Failure["LegacyPreparationMigration", <|
          "Detail" -> "v1 reduction integral cannot be normalized",
          "Key" -> oldKey|>], "LegacyV1Rekey"]];
      newKey = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
        topology, integral];
      existing = Lookup[rekeyed, Key[newKey], Missing["Absent"]];
      If[MissingQ[existing],
        AssociateTo[rekeyed, newKey -> entry],
        If[existing =!= entry,
          Throw[Failure["LegacyPreparationMigration", <|
            "Detail" -> "nonidentical v1 entries collide after exact re-key",
            "OldKey" -> oldKey, "NewKey" -> newKey|>],
            "LegacyV1Rekey"]]]
    ], Normal[oldCache]];
    rekeyed,
    "LegacyV1Rekey"];
  If[FailureQ[result], Return[result, Module]];
  If[!preparedReductionCacheQ[data, result],
    Return[Failure["LegacyPreparationMigration", <|
      "Detail" ->
        "re-keyed v1 cache does not cover every current exact boundary key"|>],
      Module]];
  result
];

legacyPreparedFTCandidate[payload_, contract_Association,
    sourceFile_String] := Module[
  {version, data, reductionCache, provenance},
  If[!AssociationQ[payload] ||
      !IntegerQ[Lookup[payload, "Key", None]], Return[$Failed, Module]];
  version = Lookup[payload, "Version", None];
  data = Lookup[payload, "FTData", None];
  reductionCache = Lookup[payload, "ReductionCache", <||>];
  Switch[version,
    2,
      provenance = <|
        "Schema" -> $ftPrepProvenanceSchema,
        "Kind" -> "LegacyV2Validated",
        "SourceSnapshotVersion" -> 2,
        "VerifiedFIRESetupProvenance" -> True,
        "VerifiedPreparationSourceProvenance" -> False|>;
      If[preparedFTDataMatchesContractQ[data, contract, provenance] &&
          preparedReductionCacheQ[data, reductionCache],
        <|"Data" -> data, "ReductionCache" -> reductionCache,
          "Provenance" -> provenance|>, $Failed],
    1,
      If[!preparedFTDataQ[data] || !AssociationQ[reductionCache],
        Return[$Failed, Module]];
      reductionCache = legacyV1RekeyReductionCache[data, reductionCache];
      If[FailureQ[reductionCache], Return[$Failed, Module]];
      provenance = <|
        "Schema" -> $ftPrepProvenanceSchema,
        "Kind" -> "LegacyV1Rekeyed",
        "SourceSnapshotVersion" -> 1,
        "VerifiedFIRESetupProvenance" -> False,
        "VerifiedPreparationSourceProvenance" -> False,
        "ReductionKeyMigration" -> "LegacyTupleToFallbackSetup/v1",
        "RekeyedEntries" -> Length[reductionCache],
        "SourceSnapshotSHA256" -> IntegerString[
          FileHash[sourceFile, "SHA256"], 16, 64]|>;
      If[preparedFTDataMatchesContractQ[data, contract, provenance] &&
          preparedReductionCacheQ[data, reductionCache],
        <|"Data" -> data, "ReductionCache" -> reductionCache,
          "Provenance" -> provenance|>, $Failed],
    _, $Failed]
];

legacyPreparedFTPayloadCompatibleQ[payload_, contract_Association] := Module[
  {version = Lookup[payload, "Version", None], provenance, data,
   reductionCache},
  (* Definitions-only compatibility seam retained for focused v2 tests.  V1
     validation additionally needs the source file hash and exact key re-map,
     so it is exercised through legacyPreparedFTCandidate. *)
  If[version =!= 2, Return[False, Module]];
  provenance = <|"Schema" -> $ftPrepProvenanceSchema,
    "Kind" -> "LegacyV2Validated", "SourceSnapshotVersion" -> 2,
    "VerifiedFIRESetupProvenance" -> True,
    "VerifiedPreparationSourceProvenance" -> False|>;
  data = Lookup[payload, "FTData", None];
  reductionCache = Lookup[payload, "ReductionCache", <||>];
  TrueQ[IntegerQ[Lookup[payload, "Key", None]] &&
    preparedFTDataMatchesContractQ[data, contract, provenance] &&
    preparedReductionCacheQ[data, reductionCache]]
];

migrateLegacyPreparedFT[name_String, targetFile_String,
    contract_Association] := Module[
  {candidates, valid = {}, payload, candidate, sourceFile,
   oldReductionCache, wrote},
  candidates = DeleteCases[
    Sort[FileNames[name <> "_*.mx", prepCacheRoot]], targetFile];
  Scan[Function[file,
    payload = readPreparedFTPayload[file];
    If[AssociationQ[payload] &&
        MemberQ[{1, 2}, Lookup[payload, "Version", None]],
      candidate = legacyPreparedFTCandidate[payload, contract, file];
      If[AssociationQ[candidate],
        AppendTo[valid, Join[<|"SourceFile" -> file|>, candidate]],
        Print["FTPREP LEGACY REJECT ", file,
          ": exact input/level/reduction validation failed"]]]],
    candidates];
  If[valid === {},
    Print["FTPREP LEGACY MIGRATION NONE for ", name];
    Return[$Failed, Module]];
  If[Length[valid] =!= 1,
    Print["FTPREP LEGACY MIGRATION REJECT ", name,
      ": ambiguous validated snapshots ", Lookup[valid, "SourceFile"]];
    Return[$Failed, Module]];
  candidate = First[valid];
  sourceFile = candidate["SourceFile"];
  oldReductionCache = FeynmanTrick`FIREInterface`Private`$ReductionCache;
  FeynmanTrick`FIREInterface`Private`$ReductionCache = Join[
    oldReductionCache, candidate["ReductionCache"]];
  wrote = savePreparedFT[targetFile, contract, candidate["Data"],
    candidate["Provenance"]];
  If[wrote === $Failed,
    FeynmanTrick`FIREInterface`Private`$ReductionCache = oldReductionCache;
    Print["FTPREP LEGACY MIGRATION WRITE FAILED ", sourceFile];
    Return[$Failed, Module]];
  Print["FTPREP LEGACY MIGRATED ", sourceFile, " -> ", targetFile];
  candidate["Data"]
];

saveLadderCheckpoint[file_, payload_] := Module[
  {tmp, saved, wrote, renamed},
  If[ladderCheckpointDir === "", Return[Null, Module]];
  If[!DirectoryQ[DirectoryName[file]],
    CreateDirectory[DirectoryName[file], CreateIntermediateDirectories -> True]];
  saved = Join[payload, ft2CheckpointRequestMetadata[], <|
    "CheckpointVersion" -> $ftLadderCheckpointVersion,
    "SourceFingerprint" -> $ftLadderSourceFingerprint|>];
  Global`$FT2LadderCheckpoint = saved;
  tmp = file <> ".tmp-" <> ToString[$ProcessID] <> ".mx";
  If[FileExistsQ[tmp], DeleteFile[tmp]];
  wrote = Quiet[Check[
    DumpSave[tmp, Global`$FT2LadderCheckpoint]; FileExistsQ[tmp], False]];
  Clear[Global`$FT2LadderCheckpoint];
  If[!TrueQ[wrote],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    Print["FTLADDER CHECKPOINT WRITE FAILED ", file];
    Return[$Failed, Module]];
  (* tmp lives beside the destination, so the overwrite is one filesystem
     rename: a killed upper-arm transport can never expose a missing or
     half-written lower-arm checkpoint. *)
  renamed = Quiet[Check[
    RenameFile[tmp, file, OverwriteTarget -> True]; True, False]];
  If[!TrueQ[renamed],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    Print["FTLADDER CHECKPOINT RENAME FAILED ", file];
    Return[$Failed, Module]];
  Print["FTLADDER CHECKPOINT ", file];
  file];

ladderCheckpointReject[file_, detail_] :=
  (Print["FTLADDER RESUME REJECT ", file, ": ", detail]; $Failed);

(* The matrix fixes only a relative epsilon gauge.  Capture that exact gauge
   once so the full-ladder planner and the runtime normalization use the same
   pole-free basis without repeating MatrixPoleOrders/FindEpsPrefactors. *)
ft2DiagonalRelativeEpsilonGauge[matrix_, epsSymbol_Symbol] := Module[
  {d, rawPoleOrders, hasRawPoles, canonical, relative, normalized,
   normalizedPoleOrders, record},
  d = Length[matrix];
  If[!MatrixQ[matrix] || Dimensions[matrix] =!= {d, d} || d === 0,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "level matrix must be a nonempty square matrix",
      "Dimensions" -> Dimensions[matrix]|>], Module]];
  rawPoleOrders =
    FeynmanTrick`EpsPrefactors`MatrixPoleOrders[matrix, epsSymbol];
  hasRawPoles = Max[Flatten[rawPoleOrders]] > 0;
  canonical = If[hasRawPoles,
    FeynmanTrick`EpsPrefactors`FindEpsPrefactors[matrix, epsSymbol],
    ConstantArray[0, d]];
  If[!ListQ[canonical] || Length[canonical] =!= d ||
      !AllTrue[canonical, IntegerQ],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "could not determine an exact integer epsilon basis",
      "Candidate" -> canonical|>], Module]];
  relative = canonical - Min[canonical];
  normalized = FeynmanTrick`EpsPrefactors`ApplyEpsPrefactors[
    matrix, relative, epsSymbol];
  normalizedPoleOrders = If[hasRawPoles,
    FeynmanTrick`EpsPrefactors`MatrixPoleOrders[normalized, epsSymbol],
    rawPoleOrders];
  If[Max[Flatten[normalizedPoleOrders]] > 0,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "diagonal epsilon normalization did not remove every matrix pole",
      "CanonicalPrefactors" -> canonical,
      "RelativePrefactors" -> relative,
      "RemainingPoleOrders" -> normalizedPoleOrders|>], Module]];
  record = <|
    "Schema" -> "FeynmanTrick.RelativeEpsilonGauge/v1",
    "InputMatrixHash" -> Hash[matrix, "SHA256"],
    "CanonicalPrefactors" -> canonical,
    "RelativePrefactors" -> relative,
    "RawPoleOrders" -> rawPoleOrders,
    "NormalizedPoleOrders" -> normalizedPoleOrders,
    "NormalizedMatrixHash" -> Hash[normalized, "SHA256"],
    "PoleFree" -> True|>;
  <|"Matrix" -> normalized, "CanonicalPrefactors" -> canonical,
    "RelativePrefactors" -> relative, "RawPoleOrders" -> rawPoleOrders,
    "NormalizedPoleOrders" -> normalizedPoleOrders,
    "Record" -> record,
    "Identity" -> ft2CanonicalIdentity["ft2-relative-epsilon-gauge-",
      record]|>];

(* A dimension-dependent FIRE basis can contain a simple apparent pole whose
   position tends to the FT endpoint as eps -> 0.  Expanding that pole as
   1/(x-r(eps)) would spend one negative epsilon order per Taylor order and,
   unless the pole is apparent, can miss eps Log[eps] terms.  Ask Indicial
   for the narrow exact projector-gauge certificate first.  When it applies,
   normalize the resulting matrix in epsilon and retain the complete
   physical-from-transport basis map:

       I = G K,       L = D K,       I = G D^-1 L.

   Observable rows and numerical boundaries use this same exact map.  The
   ordinary diagonal-only record remains byte-for-byte unchanged for every
   level without a certified apparent pole. *)
ft2RelativeEpsilonGauge[matrix_, epsSymbol_Symbol,
    physicalVar_:None] := Module[
  {apparent, diagonal, d, relative, physicalFromRelative,
   relativeFromPhysical, normalized, identity, connectionResidual,
   record, chartRef},
  If[!MatchQ[physicalVar, _Symbol],
    Return[ft2DiagonalRelativeEpsilonGauge[matrix, epsSymbol], Module]];
  chartRef = <|"Name" -> "ft-apparent@" <>
      ToString[physicalVar, InputForm],
    "Center" -> 0, "Variable" -> physicalVar|>;
  apparent = DiffExp2`Indicial`EpsilonCoalescingApparentReduce[
    matrix, physicalVar, epsSymbol, chartRef];
  If[!TrueQ[Lookup[apparent, "Applied", False]],
    Return[ft2DiagonalRelativeEpsilonGauge[matrix, epsSymbol], Module]];
  diagonal = ft2DiagonalRelativeEpsilonGauge[
    apparent["Matrix"], epsSymbol];
  If[FailureQ[diagonal], Return[diagonal, Module]];
  d = Length[matrix];
  relative = diagonal["RelativePrefactors"];
  physicalFromRelative = Map[Cancel[Together[#]] &,
    apparent["Gauge"] .
      DiagonalMatrix[epsSymbol^-relative], {2}];
  relativeFromPhysical = Map[Cancel[Together[#]] &,
    DiagonalMatrix[epsSymbol^relative] .
      apparent["GaugeInverse"], {2}];
  If[!And @@ Flatten[Map[TrueQ[PossibleZeroQ[#]] &,
      physicalFromRelative . relativeFromPhysical -
        IdentityMatrix[d], {2}]] ||
      !And @@ Flatten[Map[TrueQ[PossibleZeroQ[#]] &,
      relativeFromPhysical . physicalFromRelative -
        IdentityMatrix[d], {2}]],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" ->
        "composite apparent/epsilon basis map is not an exact two-sided inverse"|>],
      Module]];
  normalized = diagonal["Matrix"];
  connectionResidual = Map[Cancel[Together[#]] &,
    relativeFromPhysical .
      (matrix . physicalFromRelative -
        D[physicalFromRelative, physicalVar]) - normalized, {2}];
  If[!And @@ Flatten[Map[TrueQ[PossibleZeroQ[#]] &,
      connectionResidual, {2}]],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" ->
        "composite apparent/epsilon basis does not reproduce the normalized connection"|>],
      Module]];
  record = <|
    "Schema" -> "FeynmanTrick.RelativeEpsilonGauge/v2",
    "InputMatrixHash" -> Hash[matrix, "SHA256"],
    "CanonicalPrefactors" -> ConstantArray[0, d],
    "RelativePrefactors" -> ConstantArray[0, d],
    "RawPoleOrders" -> diagonal["RawPoleOrders"],
    "NormalizedPoleOrders" -> diagonal["NormalizedPoleOrders"],
    "NormalizedMatrixHash" -> Hash[normalized, "SHA256"],
    "PoleFree" -> True,
    "CompositeApparent" -> True,
    "PhysicalVariable" -> ToString[physicalVar, InputForm],
    "ApparentReduction" -> KeyDrop[apparent, {
      "Matrix", "Gauge", "GaugeInverse"}],
    "DiagonalGaugeIdentity" -> diagonal["Identity"],
    "PhysicalFromRelativeHash" ->
      Hash[physicalFromRelative, "SHA256"],
    "RelativeFromPhysicalHash" ->
      Hash[relativeFromPhysical, "SHA256"]|>;
  identity = ft2CanonicalIdentity[
    "ft2-relative-epsilon-gauge-", record];
  <|"Matrix" -> normalized,
    "CanonicalPrefactors" -> ConstantArray[0, d],
    "RelativePrefactors" -> ConstantArray[0, d],
    "RawPoleOrders" -> diagonal["RawPoleOrders"],
    "NormalizedPoleOrders" -> diagonal["NormalizedPoleOrders"],
    "PhysicalFromRelative" -> physicalFromRelative,
    "RelativeFromPhysical" -> relativeFromPhysical,
    "CompositeApparent" -> True,
    "Record" -> record, "Identity" -> identity|>];

ft2CompositeBoundaryTransform[gauge_Association,
    boundaryValues_List, boundaryPrefactors_List,
    epsSymbol_Symbol, physicalVar_Symbol, anchor_] := Module[
  {d = Length[boundaryValues], inputTop, inputSeries, transform,
   outputs, coefficient, coefficientValuation, coefficientSeries,
   term, out, commonOffset, shifted, commonTop, rows, record},
  inputTop = First[Length /@ boundaryValues] - 1;
  inputSeries = MapThread[
    DiffExp2`EpsSeries`ESShift[
      DiffExp2`EpsSeries`ESNew[0, #1], -#2] &,
    {boundaryValues, boundaryPrefactors}];
  transform = Map[Cancel[Together[# /. physicalVar -> anchor]] &,
    gauge["RelativeFromPhysical"], {2}];
  If[!FreeQ[transform, physicalVar],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" ->
        "composite boundary transform retained the Feynman parameter at the anchor"|>],
      Module]];
  outputs = Table[
    out = None;
    Do[
      coefficient = transform[[i, j]];
      If[!TrueQ[PossibleZeroQ[coefficient]],
        coefficientValuation = ft2ExactEpsilonValuation[
          coefficient, physicalVar, epsSymbol];
        If[FailureQ[coefficientValuation],
          Return[coefficientValuation, Module]];
        coefficientSeries = catch2[
          DiffExp2`EpsSeries`ESFromExpression[
            coefficient, epsSymbol,
            inputTop + coefficientValuation]];
        If[FailureQ[coefficientSeries],
          Return[coefficientSeries, Module]];
        term = DiffExp2`EpsSeries`ESTimes[
          coefficientSeries, inputSeries[[j]]];
        out = If[out === None, term,
          DiffExp2`EpsSeries`ESAdd[out, term]]],
      {j, d}];
    If[out === None,
      Return[Failure["FeynmanTrickEpsilonBasis", <|
        "Detail" ->
          "composite relative-from-physical transform has an identically zero row",
        "Row" -> i|>], Module]];
    out,
    {i, d}];
  commonOffset = Max[0,
    -Min[DiffExp2`EpsSeries`ESMinPower /@ outputs]];
  shifted = DiffExp2`EpsSeries`ESShift[#, commonOffset] & /@ outputs;
  commonTop = Min[DiffExp2`EpsSeries`ESCompleteMax /@ shifted];
  If[!IntegerQ[commonTop] || commonTop < 0,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" ->
        "composite boundary transform has no finite nonnegative common epsilon window",
      "CommonOffset" -> commonOffset,
      "Windows" -> ({
          DiffExp2`EpsSeries`ESMinPower[#],
          DiffExp2`EpsSeries`ESCompleteMax[#]} & /@ shifted)|>],
      Module]];
  rows = Table[
    Table[DiffExp2`EpsSeries`ESCoefficient[shifted[[i]], k],
      {k, 0, commonTop}],
    {i, d}];
  record = <|
    "Schema" -> "FeynmanTrick.CompositeApparentBoundary/v1",
    "GaugeIdentity" -> gauge["Identity"],
    "Anchor" -> anchor,
    "InputPrefactors" -> boundaryPrefactors,
    "InputCompleteMax" -> inputTop,
    "CommonOffset" -> commonOffset,
    "CompleteMax" -> commonTop|>;
  <|"BoundaryValues" -> rows,
    "BoundaryPrefactors" -> ConstantArray[commonOffset, d],
    "BoundaryShifts" -> ConstantArray[commonOffset, d] -
      boundaryPrefactors,
    "InputPrefactors" -> boundaryPrefactors,
    "InputCompleteMax" -> inputTop,
    "CompleteMax" -> commonTop,
    "BoundaryTransformRecord" -> record|>];

(* FIRE returns differential equations in its physical master basis I.  Some
   otherwise regular systems contain epsilon poles in that basis.  Transport
   the exactly equivalent basis J_i = eps^k_i I_i instead:

       A_J = D A_I D^-1,  D = DiagonalMatrix[eps^k_i].

   FindEpsPrefactors fixes only the relative k_i.  Its common offset is chosen
   here so that every conversion from the incoming finite boundary basis is a
   nonnegative coefficient shift.  The plain DiffExp2 boundary seam has one
   common complete upper order, so we deliberately retain exactly the input
   width.  Relative shifts can consume physical upper orders; they must be
   supplied by explicit lookahead/halos and are never filled with assumed
   zeros. *)
ft2NormalizeEpsilonBasis[matrix_, boundaryValues_List,
    boundaryPrefactors_List, epsSymbol_Symbol,
    suppliedGauge_:Automatic, physicalVar_:None,
    anchor_:None] := Module[
  {d, widths, inputTop, gauge, rawPoleOrders, canonical, relative,
   matrixHash, gaugeInputHash, gaugeNormalizedHash,
   commonOffset, effective, boundaryShifts, normalized,
   normalizedPoleOrders, shiftedBoundary, record},
  d = Length[matrix];
  ft2NativeStageTiming["epsilon-basis dimension=", d, " begin"];
  If[!MatrixQ[matrix] || Dimensions[matrix] =!= {d, d} || d === 0,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "level matrix must be a nonempty square matrix",
      "Dimensions" -> Dimensions[matrix]|>], Module]];
  If[Length[boundaryValues] =!= d ||
      Length[boundaryPrefactors] =!= d ||
      !AllTrue[boundaryValues, ListQ] ||
      !AllTrue[boundaryPrefactors, IntegerQ],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "boundary values/prefactors do not match the level matrix",
      "Dimension" -> d, "BoundaryRows" -> Length[boundaryValues],
      "Prefactors" -> boundaryPrefactors|>], Module]];
  ft2NativeStageTiming["epsilon-basis shape-ready"];
  widths = Length /@ boundaryValues;
  If[MemberQ[widths, 0] || Length[DeleteDuplicates[widths]] =!= 1,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "epsilon boundary rows must have one nonempty common window",
      "Widths" -> widths|>], Module]];
  inputTop = First[widths] - 1;
  ft2NativeStageTiming["epsilon-basis boundary-ready top=", inputTop];
  gauge = If[suppliedGauge === Automatic,
    ft2RelativeEpsilonGauge[matrix, epsSymbol, physicalVar], suppliedGauge];
  If[FailureQ[gauge], Return[gauge, Module]];
  ft2NativeStageTiming["epsilon-basis gauge-ready supplied=",
    suppliedGauge =!= Automatic];
  matrixHash = Hash[matrix, "SHA256"];
  gaugeInputHash =
    Lookup[Lookup[gauge, "Record", <||>], "InputMatrixHash", None];
  gaugeNormalizedHash =
    Lookup[Lookup[gauge, "Record", <||>], "NormalizedMatrixHash", None];
  (* The native planner owns and reuses the already gauge-normalized matrix
     to avoid repeating expensive d -> 4-2 eps rational canonicalization.
     The same exact gauge therefore legitimately arrives bound either to its
     raw input hash or to its normalized output hash.  Every other hash is
     rejected, and the normalized payload is checked again below. *)
  If[!AssociationQ[gauge] ||
      !MemberQ[{gaugeInputHash, gaugeNormalizedHash}, matrixHash] ||
      !TrueQ[Lookup[Lookup[gauge, "Record", <||>], "PoleFree", False]] ||
      Lookup[gauge, "Identity", None] =!=
        ft2CanonicalIdentity["ft2-relative-epsilon-gauge-",
          Lookup[gauge, "Record", None]],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "supplied relative epsilon gauge does not match the level matrix"|>],
    Module]];
  ft2NativeStageTiming["epsilon-basis gauge-identity-ready"];
  If[TrueQ[Lookup[gauge, "CompositeApparent", False]],
    If[!MatchQ[physicalVar, _Symbol] ||
        !TrueQ[NumericQ[anchor]] ||
        Lookup[gauge["Record"], "PhysicalFromRelativeHash", None] =!=
          Hash[Lookup[gauge, "PhysicalFromRelative", None], "SHA256"] ||
        Lookup[gauge["Record"], "RelativeFromPhysicalHash", None] =!=
          Hash[Lookup[gauge, "RelativeFromPhysical", None], "SHA256"],
      Return[Failure["FeynmanTrickEpsilonBasis", <|
        "Detail" ->
          "composite apparent gauge lacks a bound physical variable, anchor, or exact transform hashes"|>],
        Module]];
    With[{transformed = ft2CompositeBoundaryTransform[
        gauge, boundaryValues, boundaryPrefactors,
        epsSymbol, physicalVar, anchor]},
      If[FailureQ[transformed], Return[transformed, Module]];
      record = Join[<|
        "Schema" -> "FeynmanTrick.EpsilonBasis/v2",
        "CanonicalPrefactors" -> gauge["CanonicalPrefactors"],
        "RelativePrefactors" -> gauge["RelativePrefactors"],
        "Prefactors" -> transformed["BoundaryPrefactors"],
        "RawPoleOrders" -> gauge["RawPoleOrders"],
        "NormalizedPoleOrders" -> gauge["NormalizedPoleOrders"],
        "CompleteMax" -> transformed["CompleteMax"],
        "NormalizedMatrixHash" ->
          Hash[gauge["Matrix"], "SHA256"],
        "RelativeGaugeIdentity" -> gauge["Identity"],
        "CompositeApparent" -> True,
        "BoundaryTransformRecord" ->
          transformed["BoundaryTransformRecord"],
        "PoleFree" -> True|>];
      Return[Join[<|"Matrix" -> gauge["Matrix"]|>,
        transformed, <|"CheckpointRecord" -> record|>], Module]]];
  canonical = gauge["CanonicalPrefactors"];
  relative = gauge["RelativePrefactors"];
  rawPoleOrders = gauge["RawPoleOrders"];
  normalizedPoleOrders = gauge["NormalizedPoleOrders"];
  normalized = gauge["Matrix"];
  If[Length[canonical] =!= d || Length[relative] =!= d ||
      !AllTrue[Join[canonical, relative], IntegerQ] ||
      Min[relative] =!= 0 ||
      Hash[normalized, "SHA256"] =!=
        gauge["Record", "NormalizedMatrixHash"],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "relative epsilon gauge has inconsistent dimensions or hashes"|>],
    Module]];
  ft2NativeStageTiming["epsilon-basis normalized-hash-ready"];
  commonOffset = Max[boundaryPrefactors - canonical];
  effective = canonical + commonOffset;
  boundaryShifts = effective - boundaryPrefactors;
  If[!AllTrue[boundaryShifts, IntegerQ[#] && # >= 0 &],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "epsilon basis would require an unknown negative boundary shift",
      "InputPrefactors" -> boundaryPrefactors,
      "EffectivePrefactors" -> effective|>], Module]];
  shiftedBoundary = MapThread[
    Function[{shift, row},
      Table[If[n < shift, 0, row[[n - shift + 1]]],
        {n, 0, inputTop}]],
    {boundaryShifts, boundaryValues}];
  ft2NativeStageTiming["epsilon-basis boundary-shift-ready"];
  record = <|
    "Schema" -> "FeynmanTrick.EpsilonBasis/v1",
    "CanonicalPrefactors" -> canonical,
    "RelativePrefactors" -> relative,
    "Prefactors" -> effective,
    "RawPoleOrders" -> rawPoleOrders,
    "NormalizedPoleOrders" -> normalizedPoleOrders,
    "CompleteMax" -> inputTop,
    "NormalizedMatrixHash" -> Hash[normalized, "SHA256"],
    "RelativeGaugeIdentity" -> gauge["Identity"],
    "PoleFree" -> True|>;
  <|
    "Matrix" -> normalized,
    "BoundaryValues" -> shiftedBoundary,
    "BoundaryPrefactors" -> effective,
    "BoundaryShifts" -> boundaryShifts,
    "InputPrefactors" -> boundaryPrefactors,
    "InputCompleteMax" -> inputTop,
    "CompleteMax" -> inputTop,
    "CheckpointRecord" -> record|>];

ft2EpsilonBasisCheckpointQ[payload_Association] := Module[
  {record = Lookup[payload, "EpsilonBasis", None], system, values,
   prefactors, matrix, widths, schema},
  system = Lookup[payload, "System", None];
  values = Lookup[payload, "BoundaryValues", None];
  prefactors = Lookup[payload, "BoundaryPrefactors", None];
  If[!AssociationQ[record] || !AssociationQ[system] || !ListQ[values] ||
      !ListQ[prefactors] || !AllTrue[values, ListQ], Return[False, Module]];
  widths = Length /@ values;
  If[MemberQ[widths, 0] || Length[DeleteDuplicates[widths]] =!= 1,
    Return[False, Module]];
  matrix = Lookup[system, "Matrix", None];
  schema = Lookup[record, "Schema", None];
  MemberQ[{"FeynmanTrick.EpsilonBasis/v1",
      "FeynmanTrick.EpsilonBasis/v2"}, schema] &&
    TrueQ[If[schema === "FeynmanTrick.EpsilonBasis/v2",
      Lookup[record, "CompositeApparent", False] &&
        AssociationQ[Lookup[record, "BoundaryTransformRecord", None]],
      !TrueQ[Lookup[record, "CompositeApparent", False]]]] &&
    TrueQ[Lookup[record, "PoleFree", False]] &&
    Lookup[record, "Prefactors", None] === prefactors &&
    Lookup[record, "CompleteMax", None] === First[widths] - 1 &&
    Lookup[record, "NormalizedMatrixHash", None] ===
      Hash[matrix, "SHA256"]];

(* The explicit Wolfram backend updates a legacy Transport checkpoint after
   each successful endpoint arm.  Cpp instead saves a completed schema-2
   NativeTransport session plus an atomic MX sidecar before compatibility
   export; restoring that state reuses its retained line/endpoint results and
   never interprets partial Wolfram arms as native ownership.  A Boundary
   checkpoint starts directly at the named lower level.  Legacy stale data
   requires explicit opt-in, while native state always requires exact current
   source provenance. *)
ft2CheckpointExpectedMatchingDigits[kind_String, level_Integer,
    matchingDigitsByLevel_Association] := If[kind === "Boundary",
  Lookup[matchingDigitsByLevel, level + 1, matchDigits],
  Lookup[matchingDigitsByLevel, level, matchDigits]];

(* MatchingDigits names the level whose solve created a checkpoint, whereas
   PublicationDigits names the lower consumer whose requested accuracy was
   used when the numerical boundary/transport result was certified.  Both
   are part of the checkpoint's numerical identity. *)
ft2CheckpointExpectedPublicationDigits[kind_String, level_Integer,
    matchingDigitsByLevel_Association] := If[kind === "Boundary",
  Lookup[matchingDigitsByLevel, level, matchDigits],
  Lookup[matchingDigitsByLevel, level - 1, matchDigits]];

(* A level publishes the boundary consumed by level-1.  Its internal match
   solve may carry guarded extra digits, but endpoint acceptance belongs to
   the downstream consumer and must not inherit that private target. *)
ft2DownstreamPublicationDigits[level_Integer,
    matchingDigitsByLevel_Association] :=
  Lookup[matchingDigitsByLevel, level - 1, matchDigits];

ft2LevelMatchingCertificationDigits[level_Integer,
    levelMatchingDigits_Integer,
    downstreamPublicationDigits_Integer] :=
  Min[Lookup[matchingCertificationDigitsByLevel, level,
      matchingCertificationMaximumDigits],
    levelMatchingDigits,
    downstreamPublicationDigits +
      matchingCertificationSafetyDigits];

ft2FixedParameterValues[data_Association] := Module[
  {nLevels = Lookup[data, "NumLevels", None], values},
  values = Lookup[data, "FixedParamValues",
    If[IntegerQ[nLevels],
      ConstantArray[
        Lookup[data, "FixedParamValue", Missing["Absent"]], nLevels],
      $Failed]];
  If[!IntegerQ[nLevels] || nLevels < 1 ||
      !ListQ[values] || Length[values] =!= nLevels ||
      !AllTrue[values,
        (IntegerQ[#] || Head[#] === Rational) && 0 < # < 1 &],
    $Failed, values]
];

ft2LevelAnchor[data_Association, level_Integer] := Module[
  {values = ft2FixedParameterValues[data]},
  If[values === $Failed || !TrueQ[1 <= level <= Length[values]],
    $Failed, values[[level]]]
];

loadLadderCheckpoint[file_, name_, data_, prepKey_, nativePlan_:None,
    matchingDigitsByLevel_:<||>,
    publicationDigitsByLevel_:Automatic] := Module[
  {payload, ok, kind, level, levelData, belowData, mastersHere, mastersBelow,
   currentRequests, savedEO, savedFingerprint, stale, needInt, needLo, needHi,
   cachedArms, recordedArms, expectedCharts, boundaryWidths, boundaryShift,
   requiredRaw, savedRequiredRaw, preservedRaw, preservedSource, nativeRecord,
   nativePlanRecord, nativePlanIdentity, transportKindQ,
   nativeTransportRecord, nativeContract, storedRequestMetadata,
   expectedRequestMetadata, expectedMatchingDigits,
   expectedPublicationDigits, expectedAnchor,
   effectivePublicationDigitsByLevel},
  effectivePublicationDigitsByLevel = If[
    publicationDigitsByLevel === Automatic,
    matchingDigitsByLevel, publicationDigitsByLevel];
  If[!AssociationQ[matchingDigitsByLevel] ||
      !AssociationQ[effectivePublicationDigitsByLevel],
    Return[ladderCheckpointReject[file,
      "matching/publication digit profiles must be Associations"], Module]];
  If[!FileExistsQ[file],
    Return[ladderCheckpointReject[file, "file does not exist"], Module]];
  Clear[Global`$FT2LadderCheckpoint];
  ok = Quiet[Check[Get[file]; True, False]];
  If[!TrueQ[ok], Clear[Global`$FT2LadderCheckpoint];
    Return[ladderCheckpointReject[file, "could not load MX payload"], Module]];
  payload = Global`$FT2LadderCheckpoint;
  Clear[Global`$FT2LadderCheckpoint];
  If[!AssociationQ[payload],
    Return[ladderCheckpointReject[file, "payload is not an Association"], Module]];
  storedRequestMetadata = KeyTake[payload,
    $ft2CheckpointRequestMetadataKeys];
  expectedRequestMetadata = ft2CheckpointRequestMetadata[];
  If[storedRequestMetadata =!= expectedRequestMetadata,
    Return[ladderCheckpointReject[file,
      "custom family/request/output-resolution identity does not match"],
      Module]];
  If[KeyExistsQ[payload, "CheckpointVersion"] &&
      payload["CheckpointVersion"] =!= $ftLadderCheckpointVersion,
    Return[ladderCheckpointReject[file, "unsupported checkpoint version"], Module]];
  kind = Lookup[payload, "Kind", "Transport"];
  If[!MemberQ[{"Transport", "NativeTransport", "Boundary"}, kind],
    Return[ladderCheckpointReject[file, "unknown checkpoint kind"], Module]];
  If[kind === "Transport" && recurrenceBackend === "Cpp",
    Return[ladderCheckpointReject[file,
      "legacy partial-arm Transport snapshots cannot resume the retained native observable batch; resume from a numeric Boundary checkpoint"],
      Module]];
  If[kind === "NativeTransport" && recurrenceBackend =!= "Cpp",
    Return[ladderCheckpointReject[file,
      "native retained Transport snapshots require the Cpp recurrence backend"],
      Module]];
  transportKindQ = MemberQ[{"Transport", "NativeTransport"}, kind];
  level = Lookup[payload, "Level", None];
  If[!IntegerQ[level] || level < 1 || level > data["NumLevels"],
    Return[ladderCheckpointReject[file, "invalid resume level"], Module]];
  expectedAnchor = ft2LevelAnchor[data, level];
  If[expectedAnchor === $Failed,
    Return[ladderCheckpointReject[file,
      "prepared FT data has an invalid per-level anchor vector"], Module]];
  expectedMatchingDigits = ft2CheckpointExpectedMatchingDigits[
    kind, level, matchingDigitsByLevel];
  expectedPublicationDigits = ft2CheckpointExpectedPublicationDigits[
    kind, level, effectivePublicationDigitsByLevel];
  If[Lookup[payload, "Example", None] =!= name,
    Return[ladderCheckpointReject[file, "example does not match"], Module]];
  If[Lookup[payload, "WorkingPrecision", None] =!= wp,
    Return[ladderCheckpointReject[file, "WorkingPrecision does not match"], Module]];
  If[Lookup[payload, "MatchingDigits", matchDigits] =!=
      expectedMatchingDigits,
    Return[ladderCheckpointReject[file, "MatchingDigits does not match"], Module]];
  If[Lookup[payload, "PublicationDigits", matchDigits] =!=
      expectedPublicationDigits,
    Return[ladderCheckpointReject[file,
      "PublicationDigits does not match"], Module]];
  If[KeyExistsQ[payload, "RecurrenceBackend"] &&
      payload["RecurrenceBackend"] =!= recurrenceBackend,
    Return[ladderCheckpointReject[file,
      "recurrence backend mode does not match"], Module]];
  If[Lookup[payload, "Anchor", None] =!= expectedAnchor,
    Return[ladderCheckpointReject[file, "anchor does not match"], Module]];
  If[KeyExistsQ[payload, "PrepKey"] && payload["PrepKey"] =!= prepKey,
    Return[ladderCheckpointReject[file, "prepared FT data does not match"], Module]];
  If[KeyExistsQ[payload, "EpsilonOrder"] &&
      payload["EpsilonOrder"] =!= epsOrder,
    Return[ladderCheckpointReject[file, "EpsilonOrder does not match"], Module]];
  If[KeyExistsQ[payload, "BoundaryExtraOrder"] &&
      payload["BoundaryExtraOrder"] =!= boundaryExtraOrder,
    Return[ladderCheckpointReject[file, "BoundaryExtraOrder does not match"], Module]];
  If[KeyExistsQ[payload, "LevelEpsilonHalos"] &&
      payload["LevelEpsilonHalos"] =!= levelEpsilonHalos,
    Return[ladderCheckpointReject[file, "level epsilon halos do not match"], Module]];
  If[Lookup[payload, "DeltaPrescriptionSign", 1] =!=
      deltaPrescriptionSign,
    Return[ladderCheckpointReject[file,
      "delta prescription sign does not match"], Module]];
  If[Lookup[payload, "LevelDeltaPrescriptionSigns",
        ConstantArray[
          Lookup[payload, "DeltaPrescriptionSign", 1],
          data["NumLevels"]]] =!=
      levelDeltaPrescriptionSignsForCount[data["NumLevels"]],
    Return[ladderCheckpointReject[file,
      "per-level delta prescription signs do not match"], Module]];
  levelData = data["Levels"][level];
  mastersHere = levelData["Masters"];
  If[Lookup[payload, "MastersHere", mastersHere] =!= mastersHere,
    Return[ladderCheckpointReject[file, "level masters do not match"], Module]];
  If[!ListQ[Lookup[payload, "BoundaryValues", None]] ||
      Length[payload["BoundaryValues"]] =!= Length[mastersHere] ||
      !ListQ[Lookup[payload, "BoundaryPrefactors", None]] ||
      Length[payload["BoundaryPrefactors"]] =!= Length[mastersHere],
    Return[ladderCheckpointReject[file, "boundary vector has the wrong dimension"],
      Module]];
  savedEO = If[transportKindQ, Lookup[payload, "ExpansionOrder", None],
    Lookup[payload, "SourceExpansionOrder", None]];
  If[!IntegerQ[savedEO] || savedEO < 1,
    Return[ladderCheckpointReject[file, "missing source ExpansionOrder"], Module]];
  (* Transport checkpoints contain order-specific chart series and therefore
     require exact order parity.  Boundary checkpoints contain only endpoint
     Laurent data: a higher-order source may safely seed a lower-order run,
     but lower-order data must never silently downgrade the requested order. *)
  If[transportKindQ && expansionOrder =!= savedEO,
    Return[ladderCheckpointReject[file,
      "transport ExpansionOrder does not match"], Module]];
  If[kind === "Boundary" && savedEO < expansionOrder,
    Return[ladderCheckpointReject[file,
      "boundary source ExpansionOrder is lower than requested"], Module]];
  If[transportKindQ,
    If[KeyExistsQ[payload, "DivisionOrder"] &&
        payload["DivisionOrder"] =!= divisionOrder,
      Return[ladderCheckpointReject[file,
        "DivisionOrder does not match"], Module]];
    If[KeyExistsQ[payload, "RadiusOfConvergence"] &&
        payload["RadiusOfConvergence"] =!= radiusOfConvergence,
      Return[ladderCheckpointReject[file,
        "RadiusOfConvergence does not match"], Module]];
    If[KeyExistsQ[payload, "ValueTransportMode"] &&
        payload["ValueTransportMode"] =!= valueTransportMode,
      Return[ladderCheckpointReject[file,
        "DE2_VALUE_TRANSPORT mode does not match"], Module]];
    If[KeyExistsQ[payload, "NativeValueHopExecutionMode"] &&
        payload["NativeValueHopExecutionMode"] =!=
          nativeValueHopExecutionMode,
      Return[ladderCheckpointReject[file,
        "DE2_NATIVE_VALUE_HOP_EXECUTION mode does not match"], Module]];
    If[KeyExistsQ[payload, "SingularMatchPrecondition"] &&
        payload["SingularMatchPrecondition"] =!= singularMatchPrecondition,
      Return[ladderCheckpointReject[file,
        "singular match precondition mode does not match"], Module]];
    belowData = data["Levels"][level - 1];
    mastersBelow = belowData["Masters"];
    currentRequests =
      FeynmanTrick`LevelReduction`BoundaryRequestRecords[
        mastersBelow, levelData["CombinedPositions"]];
    If[Lookup[payload, "Variable", None] =!= levelData["FeynmanParameter"] ||
        Lookup[payload, "MastersBelow", None] =!= mastersBelow ||
        Lookup[payload, "Requests", None] =!= currentRequests,
      Return[ladderCheckpointReject[file,
        "transport level metadata does not match prepared FT data"], Module]];
    requiredRaw = If[kind === "NativeTransport" &&
        ft2NativeEpsilonPlanQ[nativePlan],
      nativePlan["Levels"][level]["RequiredRawTop"],
      requestedEpsilonOrder[level]];
    If[Lookup[payload, "RequestedEpsilonOrder", None] =!= requiredRaw,
      Return[ladderCheckpointReject[file,
        "requested epsilon window does not match"], Module]];
    If[!AssociationQ[Lookup[payload, "System", None]] ||
        !ft2EpsilonBasisCheckpointQ[payload] ||
        !AssociationQ[Lookup[payload, "Reductions", None]] ||
        !AllTrue[currentRequests,
          KeyExistsQ[payload["Reductions"], #["NeededVec"]] &] ||
        !ListQ[Lookup[payload, "ExtraSingularFactors", None]],
      Return[ladderCheckpointReject[file,
        "transport payload is incomplete or has inconsistent epsilon-basis metadata"],
        Module]];
    If[kind === "Transport",
      If[!ListQ[Lookup[payload, "ChartCache", None]] ||
          !AllTrue[{Lookup[payload, "TransportLow", None],
              Lookup[payload, "TransportHigh", None]},
            (# === None || AssociationQ[#]) &],
        Return[ladderCheckpointReject[file,
          "legacy transport payload has malformed arm snapshots"], Module]];
      needInt = AnyTrue[currentRequests, #["Case"] === "integrate" &];
      needLo = needInt || AnyTrue[currentRequests, #["Case"] === "limitLower" &];
      needHi = needInt || AnyTrue[currentRequests, #["Case"] === "limitUpper" &];
      cachedArms = Pick[{"Lower", "Upper"},
        AssociationQ /@ {payload["TransportLow"], payload["TransportHigh"]}];
      recordedArms = Lookup[payload, "CompletedArms", cachedArms];
      If[recordedArms =!= cachedArms,
        Return[ladderCheckpointReject[file,
          "completed-arm metadata does not match transport payload"], Module]];
      expectedCharts = Join[
        If[AssociationQ[payload["TransportLow"]],
          Lookup[payload["TransportLow"], "Charts", {}], {}],
        If[AssociationQ[payload["TransportHigh"]],
          Lookup[payload["TransportHigh"], "Charts", {}], {}]];
      If[payload["ChartCache"] =!= expectedCharts,
        Return[ladderCheckpointReject[file,
          "chart cache does not match completed transport arms"], Module]];
      If[(needLo || needHi) && cachedArms === {},
        Print["FTLADDER RESUME transport has no completed endpoint arms; ",
          "both required arms will be computed"]],
      nativeTransportRecord = Lookup[payload,
        "NativeTransportCheckpoint", None];
      nativeContract = Lookup[payload, "NativeTransportContract", None];
      nativeRecord = Lookup[payload, "NativeObservableBatch", None];
      nativePlanRecord = Lookup[payload, "NativeEpsilonPlan", None];
      nativePlanIdentity = Lookup[payload,
        "NativeEpsilonPlanIdentity", None];
      If[!ft2NativeTransportResumeRecordQ[nativeTransportRecord] ||
          !ft2NativeTransportContractQ[nativeContract] ||
          nativeTransportRecord["ContractIdentity"] =!=
            nativeContract["Identity"] ||
          nativeContract["Record", "SourceFingerprint"] =!=
            $ftLadderSourceFingerprint ||
          nativeContract["Record", "Example"] =!= name ||
          nativeContract["Record", "Level"] =!= level ||
          nativeContract["Record", "PrepKey"] =!= prepKey ||
          !ft2NativeCheckpointRecordQ[nativeRecord] ||
          !AssociationQ[Lookup[payload, "NativeLedger", None]] ||
          !ft2NativeEpsilonExecutionRecordQ[nativePlanRecord,
            nativePlanIdentity, nativePlan] ||
          nativeContract["Record", "NativeEpsilonPlanIdentity"] =!=
            nativePlanIdentity ||
          nativeRecord["NativeEpsilonPlanIdentity"] =!= nativePlanIdentity ||
          nativeRecord["AtlasPlanIdentity"] =!=
            nativeTransportRecord["AtlasPlanIdentity"] ||
          !StringQ[Lookup[nativeTransportRecord["State"], "Path", None]] ||
          !FileExistsQ[nativeTransportRecord["State", "Path"]],
        Return[ladderCheckpointReject[file,
          "native transport checkpoint contract, batch, or epsilon-plan identity is inconsistent"],
          Module]]]];
    If[kind === "Boundary" && recurrenceBackend === "Cpp",
      boundaryWidths = If[AllTrue[payload["BoundaryValues"], ListQ],
        Length /@ payload["BoundaryValues"], {}];
      boundaryShift = If[payload["BoundaryPrefactors"] =!= {} &&
          SameQ @@ payload["BoundaryPrefactors"] &&
          IntegerQ[First[payload["BoundaryPrefactors"]]] &&
          First[payload["BoundaryPrefactors"]] >= 0,
        First[payload["BoundaryPrefactors"]], None];
      requiredRaw = If[ft2NativeEpsilonPlanQ[nativePlan],
        nativePlan["Levels"][level]["RequiredRawTop"],
        nativeRequiredRawTop[level]];
      savedRequiredRaw = Lookup[payload, "RequiredRawTop", None];
      preservedRaw = Lookup[payload, "PreservedRawCompleteMax", None];
      preservedSource =
        Lookup[payload, "PreservedSourceCompleteMax", None];
      nativeRecord = Lookup[payload, "NativeObservableBatch", None];
      nativePlanRecord = Lookup[payload, "NativeEpsilonPlan", None];
      nativePlanIdentity = Lookup[payload,
        "NativeEpsilonPlanIdentity", None];
      If[boundaryWidths === {} || MemberQ[boundaryWidths, 0] ||
          Length[DeleteDuplicates[boundaryWidths]] =!= 1 ||
          !IntegerQ[boundaryShift] ||
          Lookup[payload, "BoundaryShift", None] =!= boundaryShift ||
          !IntegerQ[savedRequiredRaw] ||
          savedRequiredRaw < requiredRaw ||
          Lookup[payload, "RequestedEpsilonOrder", None] =!=
            savedRequiredRaw ||
          !IntegerQ[preservedRaw] || preservedRaw < savedRequiredRaw ||
          !IntegerQ[preservedSource] ||
          preservedSource =!= preservedRaw + boundaryShift ||
          preservedSource =!= First[boundaryWidths] - 1 ||
          !ft2NativeCheckpointRecordQ[nativeRecord] ||
          nativeRecord["RequiredRawTop"] =!= savedRequiredRaw ||
          nativeRecord["DeliverableCompleteMax"] =!= preservedRaw ||
          (ft2NativeEpsilonPlanQ[nativePlan] &&
            (!AssociationQ[nativePlanRecord] ||
              !StringQ[nativePlanIdentity] ||
              nativePlanIdentity =!= ft2CanonicalIdentity[
                "ft2-native-epsilon-execution-plan-",
                nativePlanRecord] ||
              Lookup[nativeRecord, "NativeEpsilonPlanIdentity", None] =!=
                nativePlanIdentity)),
        Return[ladderCheckpointReject[file,
          "native Boundary checkpoint has an insufficient or internally inconsistent required floor, preserved source width, shift, or observable-batch identity"],
          Module]],
      If[kind === "Boundary" &&
          KeyExistsQ[payload, "RequestedEpsilonOrder"] &&
          payload["RequestedEpsilonOrder"] =!=
            requestedEpsilonOrder[level],
        Return[ladderCheckpointReject[file,
          "requested epsilon window does not match"], Module]]];
  savedFingerprint = Lookup[payload, "SourceFingerprint", Missing["Unversioned"]];
  stale = savedFingerprint =!= $ftLadderSourceFingerprint ||
    TrueQ[Lookup[payload, "Tainted", False]];
  If[kind === "NativeTransport" && stale,
    Return[ladderCheckpointReject[file,
      "native retained transport state requires exact current source provenance"],
      Module]];
  If[stale && !allowStaleLadderCheckpoint,
    Return[ladderCheckpointReject[file,
      "source provenance is stale/unversioned; set FT_ALLOW_STALE_LADDER_CHECKPOINT=1 to opt in"],
      Module]];
  If[stale,
    Print["FTLADDER RESUME WARNING stale/unversioned checkpoint explicitly accepted"]];
  Join[payload, <|"Kind" -> kind, "Tainted" -> stale,
    "SourceExpansionOrder" -> savedEO|>]];

ft2DiscoveredLadderCheckpoint[name_String, data_, prepKey_, nativePlan_,
    matchingDigitsByLevel_:<||>,
    publicationDigitsByLevel_:Automatic] :=
 Module[{records, record, loaded, levelFromFile, add,
    effectivePublicationDigitsByLevel},
  effectivePublicationDigitsByLevel = If[
    publicationDigitsByLevel === Automatic,
    matchingDigitsByLevel, publicationDigitsByLevel];
  If[!AssociationQ[matchingDigitsByLevel] ||
      !AssociationQ[effectivePublicationDigitsByLevel],
    Return[$Failed, Module]];
  If[ladderCheckpointDir === "" || !DirectoryQ[ladderCheckpointDir],
    Return[None, Module]];
  levelFromFile[file_String, suffix_String] := Module[{matches},
    matches = StringCases[FileNameTake[file],
      RegularExpression["^.*_level([0-9]+)_" <> suffix <> "\\.mx$"] ->
        "$1"];
    If[Length[matches] === 1, FromDigits[First[matches]], Infinity]];
  add[pattern_String, kind_String, suffix_String, phase_Integer] :=
    Map[<|"File" -> #, "Kind" -> kind,
        "Level" -> levelFromFile[#, suffix], "Phase" -> phase,
        "Date" -> Quiet[Check[AbsoluteTime[FileDate[#]], 0]]|> &,
      FileNames[pattern, ladderCheckpointDir]];
  records = Join[
    add[name <> "_level*_native_transport.mx", "NativeTransport",
      "native_transport", 2],
    If[recurrenceBackend === "Cpp", {},
      add[name <> "_level*_transport.mx", "Transport", "transport", 1]],
    add[name <> "_level*_boundary.mx", "Boundary", "boundary", 0]];
  records = SortBy[Select[records, IntegerQ[#["Level"]] &],
    {#["Level"], -#["Phase"], -#["Date"]} &];
  Do[
    loaded = loadLadderCheckpoint[record["File"], name, data, prepKey,
      nativePlan, matchingDigitsByLevel,
      effectivePublicationDigitsByLevel];
    If[AssociationQ[loaded],
      Print["FTLADDER AUTO RESUME kind=", loaded["Kind"],
        " level=", loaded["Level"], " file=", record["File"]];
      Return[loaded, Module]],
    {record, records}];
  None];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", wp];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds",
  runnerSettings["FIRETimeoutSeconds"]];

eps = Global`eps;
esC = DiffExp2`EpsSeries`ESCoefficient;
esMn = DiffExp2`EpsSeries`ESMinPower;
esCMx = DiffExp2`EpsSeries`ESCompleteMax;
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

ft2EpsilonTrimAmbiguityFailureQ[failure_] :=
  FailureQ[failure] &&
  failure[[1]] === "DiffExp2" &&
  AssociationQ[failure[[2]]] &&
  Lookup[failure[[2]], "ID", None] === "E5" &&
  Lookup[failure[[2]], "Module", None] === "Tolerances" &&
  Lookup[failure[[2]], "Detail", None] ===
    "cannot classify coefficient against scale within the ambiguity band";

cleanNumber[value_] := Module[{n = N[value, 50], re, im},
  If[!NumericQ[n], Return[ToString[InputForm[value]]]];  (* JSON-safe *)
  re = Re[n]; im = Im[n];
  If[Abs[im] < 10^-30, re, <|"Re" -> re, "Im" -> im|>]];

ft2NotApplicableCertification[operation_String, reason_String] := <|
  "Applicability" -> "not-applicable", "Operation" -> operation,
  "Reason" -> reason|>;

ft2OutputCertificationQ[record_] := AssociationQ[record] && Switch[
  Lookup[record, "Applicability", None],
  "applicable",
    Sort[Keys[record]] === Sort[{"Applicability", "Operation", "Scope",
        "ErrorGuarantee", "ErrorEnvelope"}] &&
      Lookup[record, "Operation", None] === "integrate" &&
      (Which[
        Lookup[record, "Scope", None] ===
            "full_local_with_certified_tail",
          Lookup[record, "ErrorGuarantee", None] === "certified" &&
            AssociationQ[Lookup[record, "ErrorEnvelope", None]] &&
            Lookup[record["ErrorEnvelope"], "guarantee", None] ===
              "certified",
        Lookup[record, "Scope", None] === "stored_truncation",
          Lookup[record, "ErrorGuarantee", None] === "none" &&
            Lookup[record, "ErrorEnvelope", None] === Null,
        True, False]),
  "not-applicable",
    Sort[Keys[record]] ===
      Sort[{"Applicability", "Operation", "Reason"}] &&
      StringQ[Lookup[record, "Operation", None]] &&
      StringQ[Lookup[record, "Reason", None]],
  _, False];

ft2NativeIntegrationCertification[exportedResult_] := Module[
  {scope, guarantee, envelope, record},
  If[!AssociationQ[exportedResult] ||
      !And @@ (KeyExistsQ[exportedResult, #] & /@
        {"Scope", "ErrorGuarantee", "ErrorEnvelope"}),
    Return[ft2NativeFailure[
      "native integral export omitted its certification record"], Module]];
  scope = exportedResult["Scope"];
  guarantee = exportedResult["ErrorGuarantee"];
  envelope = exportedResult["ErrorEnvelope"];
  (* None is the native in-kernel sentinel.  JSON has no representation for
     an arbitrary Wolfram symbol, so expose the same absence as JSON null. *)
  If[envelope === None, envelope = Null];
  record = <|"Applicability" -> "applicable",
    "Operation" -> "integrate", "Scope" -> scope,
    "ErrorGuarantee" -> guarantee, "ErrorEnvelope" -> envelope|>;
  If[ft2OutputCertificationQ[record], record,
    ft2NativeFailure[
      "native integral export returned malformed certification metadata",
      <|"Certification" -> record|>]]];

ft2CertificationForEntry[entry_Association, exportedResult_] := Which[
  TrueQ[Lookup[entry, "ProvenZero", False]],
    ft2NotApplicableCertification[Lookup[entry, "Case", "unknown"],
      "proven-zero"],
  Lookup[entry, "Case", None] === "integrate",
    ft2NativeIntegrationCertification[exportedResult],
  MemberQ[{"limitLower", "limitUpper"}, Lookup[entry, "Case", None]],
    ft2NotApplicableCertification[entry["Case"], "endpoint-limit"],
  Lookup[entry, "Case", None] === "direct",
    ft2NotApplicableCertification["direct", "direct"],
  True,
    ft2NativeFailure["cannot classify native output certification",
      <|"Entry" -> entry|>]];

ft2RequestedCoefficientRows[raw_, requestedTop_Integer,
    outputKind_String] := Module[{rowMin, availableTop, outputTop},
  If[requestedTop < 0,
    Return[Failure["FeynmanTrickOutputCompleteness", <|
      "Detail" -> outputKind <>
        " output requires a nonnegative integer requested epsilon order",
      "RequestedCompleteMax" -> requestedTop|>], Module]];
  rowMin = esMn[raw];
  availableTop = esCMx[raw];
  If[availableTop < requestedTop,
    Return[Failure["FeynmanTrickOutputCompleteness", <|
      "Detail" -> outputKind <>
        " output does not cover the requested epsilon order",
      "RequestedCompleteMax" -> requestedTop,
      "AvailableCompleteMax" -> availableTop,
      "RawMinPower" -> rowMin|>], Module]];
  outputTop = Min[requestedTop, availableTop];
  Table[{power, cleanNumber[esC[raw, power]]},
    {power, rowMin, outputTop}]];

ft2RequestedCoefficientRows[raw_, requestedTop_, outputKind_String] :=
  Failure["FeynmanTrickOutputCompleteness", <|
    "Detail" -> outputKind <>
      " output requires a nonnegative integer requested epsilon order",
    "RequestedCompleteMax" -> requestedTop|>];

ft2ExpectedPoleFreeExampleQ[example_String] := MemberQ[
  {"bubble", "sunrise", "banana", "banana_unequal", "banana4",
    "banana4_unequal", "kite", "pentagon_massive"}, example];

ft2ExpectedFinalMinimumPower[example_String] := Which[
  ft2ExpectedPoleFreeExampleQ[example], 0,
  example === "double_box_planar", -4,
  True, Missing["UnknownFinalPoleFloor", example]];

ft2FinalPoleNegligibleQ[coefficient_] := Module[{bounds},
  If[FreeQ[coefficient, _?InexactNumberQ] &&
      TrueQ[PossibleZeroQ[coefficient]], Return[True, Module]];
  If[!NumericQ[coefficient], Return[False, Module]];
  bounds = catch2[
    DiffExp2`Tolerances`NumericMagnitudeBounds[coefficient, 30]];
  ListQ[bounds] && Length[bounds] === 2 &&
    TrueQ[Last[bounds] <= 10^-8]];

ft2UnexpectedFinalPoleRows[example_String, raw_] := Module[
  {minimum = esMn[raw],
   expectedMinimum = ft2ExpectedFinalMinimumPower[example], inspected},
  If[MissingQ[expectedMinimum] || minimum >= expectedMinimum, Return[{}]];
  inspected = Table[{power, esC[raw, power]},
    {power, minimum, expectedMinimum - 1}];
  Select[inspected, !ft2FinalPoleNegligibleQ[Last[#]] &]];

ft2PretrimFinalPoleAudit[example_String, masters_List,
    rawValues_List] := Module[{bad, expectedMinimum},
  expectedMinimum = ft2ExpectedFinalMinimumPower[example];
  If[MissingQ[expectedMinimum], Return[True, Module]];
  If[Length[masters] =!= Length[rawValues],
    Return[Failure["FeynmanTrickUnexpectedFinalPole", <|
      "Detail" ->
        "pre-trim final pole audit received inconsistent master/value counts",
      "Example" -> example, "Masters" -> masters,
      "ValueCount" -> Length[rawValues]|>], Module]];
  bad = MapThread[Function[{master, raw}, Module[{rows},
      rows = ft2UnexpectedFinalPoleRows[example, raw];
      If[rows === {}, Nothing, <|
        "Master" -> master,
        "NegativeCoefficients" ->
          ({First[#], cleanNumber[N[Last[#], Min[50, wp]]]} & /@ rows)|>]]],
    {masters, rawValues}];
  If[bad === {}, True,
    Failure["FeynmanTrickUnexpectedFinalPole", <|
      "Detail" ->
        "a physical example has unresolved coefficients below its proven final pole floor before ESTrimThrough; refusing to let trimming hide them",
      "Example" -> example, "ExpectedMinimumPower" -> expectedMinimum,
      "AbsoluteTolerance" -> 10^-8,
      "Failures" -> bad|>]]];

ft2ApplyExpectedFinalPoleFloor[example_String, raw_] := Module[
  {expectedMinimum = ft2ExpectedFinalMinimumPower[example],
   completeMax = esCMx[raw]},
  If[MissingQ[expectedMinimum] || esMn[raw] >= expectedMinimum,
    Return[raw, Module]];
  If[completeMax < expectedMinimum,
    Return[Failure["FeynmanTrickUnexpectedFinalPole", <|
      "Detail" ->
        "the retained final epsilon window ends below its proven pole floor",
      "Example" -> example, "ExpectedMinimumPower" -> expectedMinimum,
      "AvailableCompleteMax" -> completeMax|>], Module]];
  (* The audit above proves that every discarded numerical enclosure is
     below the public absolute tolerance.  The example-specific pole theorem
     supplies the stronger exact statement: these coefficients are
     structural zeros, so advance Min instead of carrying numerical
     cancellation remnants into the published Laurent window. *)
  DiffExp2`EpsSeries`ESNew[expectedMinimum,
    Table[esC[raw, power], {power, expectedMinimum, completeMax}]]];

ft2StepwiseRow[example_, level_, master_, raw_, prefactor_,
    certification_] := Module[{rowMin, coefficients},
  If[!ft2OutputCertificationQ[certification],
    Return[Failure["FeynmanTrickOutputCertification", <|
      "Detail" -> "STEPWISE certification is malformed",
      "Certification" -> certification|>], Module]];
  rowMin = esMn[raw];
  coefficients = ft2RequestedCoefficientRows[raw, epsOrder, "STEPWISE"];
  If[FailureQ[coefficients], Return[coefficients, Module]];
  <|"Example" -> example, "Level" -> level, "Master" -> master,
    "EpsPrefactor" -> prefactor, "RawMinPower" -> rowMin,
    "Coefficients" -> coefficients,
    "Certification" -> certification|>];

ft2FinalRow[example_, raw_, certification_] := Module[
  {coefficients, result, unexpectedPoles},
  If[!ft2OutputCertificationQ[certification],
    Return[Failure["FeynmanTrickOutputCertification", <|
      "Detail" -> "FINAL certification is malformed",
      "Certification" -> certification|>], Module]];
  unexpectedPoles = ft2UnexpectedFinalPoleRows[example, raw];
  If[unexpectedPoles =!= {},
    Return[Failure["FeynmanTrickUnexpectedFinalPole", <|
      "Detail" -> "a pole-free physical example has unresolved negative epsilon coefficients; refusing to publish an unchecked FINAL value",
      "Example" -> example, "Tolerance" -> 10^-8,
      "NegativeCoefficients" ->
        ({First[#], cleanNumber[Last[#]]} & /@ unexpectedPoles),
      "IntegrationCertification" -> certification|>], Module]];
  coefficients = ft2RequestedCoefficientRows[raw, epsOrder, "FINAL"];
  If[FailureQ[coefficients], Return[coefficients, Module]];
  result = <|"Example" -> example,
    "Finite" -> cleanNumber[esC[raw, 0]],
    "RawMinPower" -> esMn[raw],
    "Certification" -> certification|>;
  (* EpsilonOrder==0 retains the exact historical FINAL association.  A
     positive request adds the complete Laurent row without replacing the
     long-standing Finite compatibility field. *)
  If[epsOrder > 0,
    Append[result, "Coefficients" -> coefficients], result]];

ft2CustomFinalRows[example_, familyID_String, pipelineRequestID_String,
    masters_List, rawValues_List, certifications_List,
    outputRequests_List] := Module[{n, rows},
  n = Length[masters];
  If[n === 0 || Length[rawValues] =!= n ||
      Length[certifications] =!= n || Length[outputRequests] =!= n ||
      !And @@ MapThread[
        Function[{request, master, ordinal},
          AssociationQ[request] &&
            Lookup[request, "IndexVector", None] === master &&
            Lookup[request, "RequestOrdinal", None] === ordinal &&
            StringQ[Lookup[request, "RequestID", None]] &&
            StringQ[Lookup[request, "PhysicalIntegralID", None]]],
        {outputRequests, masters, Range[n]}],
    Return[Failure["FeynmanTrickOutputCertification", <|
      "Detail" -> "custom FINAL values do not match the ordered output-request contract",
      "Masters" -> masters,
      "OutputRequests" -> outputRequests|>], Module]];
  rows = MapThread[ft2FinalRow,
    {ConstantArray[example, n], rawValues, certifications}];
  If[AnyTrue[rows, FailureQ],
    Return[First[Select[rows, FailureQ]], Module]];
  MapThread[
    Function[{row, master, request}, Join[row, <|
        "Master" -> master,
        "PipelineRequestID" -> pipelineRequestID,
        "FamilyID" -> familyID,
        "RequestID" -> request["RequestID"],
        "RequestOrdinal" -> request["RequestOrdinal"],
        "PhysicalIntegralID" -> request["PhysicalIntegralID"]|>,
      KeyTake[request, {"SelectionRequestID", "ResolutionID"}]]],
    {rows, masters, outputRequests}]
];

ft2OutputLine[prefix_String, row_Association] := Module[{json},
  json = Quiet[Check[ExportString[
    row /. x_Rational :> N[x, 50], "RawJSON", "Compact" -> True],
    $Failed]];
  If[StringQ[json], prefix <> json,
    Failure["FeynmanTrickOutputCertification", <|
      "Detail" -> "output row is not JSON serializable", "Row" -> row|>]]];

printRows[example_, level_, masters_, rawES_List, prefactors_,
    certifications_List] := Module[{rows, lines},
  If[Length[masters] =!= Length[rawES] ||
      Length[masters] =!= Length[prefactors] ||
      Length[masters] =!= Length[certifications],
    Return[Failure["FeynmanTrickOutputCertification", <|
      "Detail" -> "STEPWISE row inputs have inconsistent lengths"|>],
      Module]];
  rows = MapThread[ft2StepwiseRow,
    {ConstantArray[example, Length[masters]],
      ConstantArray[level, Length[masters]], masters, rawES, prefactors,
      certifications}];
  If[AnyTrue[rows, FailureQ], Return[First[Select[rows, FailureQ]], Module]];
  lines = ft2OutputLine["STEPWISE ", #] & /@ rows;
  If[AnyTrue[lines, FailureQ], Return[First[Select[lines, FailureQ]], Module]];
  Scan[Print, lines];
  rows];

(* combined endpoint limit: lim Sum_j c_j(x) f_j(x) at the chart center.
   Assemble the scalar combination first: testing each master separately
   rejects legitimate cancellations between endpoint-singular terms. *)
limitCombined[tres_, cvec_, var_] := Module[
  {ls = tres["Final"], pieces = {}, active = {}, scalar, scale},
  scale = ls["ChartMap", "Scale"];
  Do[Module[{cc = cvec[[j]], lsP, lsM},
    If[!PossibleZeroQ[Together[cc]],
      If[Environment["DEBUG_LI"] === "1",
        Print["      limit j=", j, " mul start t=", SessionTime[]]];
      lsP = Join[ls, <|"Sectors" -> Map[
        Join[#, <|"Coeffs" -> #["Coeffs"][[All, All, {j}]]|>] &,
        ls["Sectors"]]|>];
      lsM = DiffExp2`SectorSeries`MultiplyRational[lsP,
        Together[cc /. var -> ls["Center"] + scale*Global`t], Global`t];
      AppendTo[pieces, lsM]; AppendTo[active, j]]],
    {j, Length[cvec]}];
  If[pieces === {}, Return[DiffExp2`EpsSeries`ESZero[
    ls["EpsWindow", "CompleteMax"]], Module]];
  scalar = If[Length[pieces] === 1, First[pieces],
    DiffExp2`SectorSeries`CombineLocalSolutions[
      ConstantArray[1, Length[pieces]], pieces]];
  If[Environment["DEBUG_LI"] === "1",
    Print["      combined limit components=", active,
      " lim start t=", SessionTime[]]];
  DiffExp2`Integrate`EndpointSectorLimit[scalar][[1]]];

(* Exact, runner-local construction for the retained native observable seam.
   Keeping this data preparation separate from runExample makes it possible
   to test the production contract without starting FIRE. *)
ft2NativeFailure[detail_String, data_:<||>] :=
  Failure["FeynmanTrickNativeBoundary", Join[<|"Detail" -> detail|>, data]];

ft2NativeMatchingReservoirRetry[failure_?FailureQ,
    level_Integer] := Module[{data, backend, additional},
  data = Quiet[Check[failure[[2]], <||>]];
  backend = If[AssociationQ[data],
    Lookup[data, "BackendFailure", None], None];
  additional = If[AssociationQ[backend] &&
      Lookup[backend, "reason", None] ===
        "acb_match_residual_inconclusive" &&
      TrueQ[Lookup[backend, "retryable_epsilon_reservoir", False]],
    Lookup[backend, "required_additional_epsilon_orders", None], None];
  If[IntegerQ[additional] && additional > 0,
    Failure["FeynmanTrickNativeMatchingReservoir", <|
      "Detail" -> "native matching needs a wider private epsilon reservoir",
      "Level" -> level, "AdditionalOrders" -> additional,
      "BackendFailure" -> backend|>], None]];

ft2NativeMatchingReservoirRetry[_, _Integer] := None;

ft2NativeMatchingReservoirRetryQ[failure_] := FailureQ[failure] &&
  Quiet[Check[failure[[1]] ===
    "FeynmanTrickNativeMatchingReservoir", False]];

(* A long ordinary arm can lose one private high coefficient at each
   finite-Laurent materialization.  If a wider retry certifies the previous
   obstruction and the failure moves to a later match, adding one order at a
   time would replay the whole arm once per chart.  Exponentially back off the
   private halo in that demonstrated-progress case.  This is conservative:
   the halo is not a publication request, and every retry must still pass the
   same coefficientwise residual certificates. *)
ft2NativeMatchingReservoirBackoff[current_Integer, additional_Integer,
    previousBackend_, currentBackend_] := Module[
  {previousArm, currentArm, previousMatch, currentMatch},
  If[current < 0 || additional < 1,
    Return[$Failed, Module]];
  previousArm = If[AssociationQ[previousBackend],
    Lookup[previousBackend, "arm", None], None];
  currentArm = If[AssociationQ[currentBackend],
    Lookup[currentBackend, "arm", None], None];
  previousMatch = If[AssociationQ[previousBackend],
    Lookup[previousBackend, "match", None], None];
  currentMatch = If[AssociationQ[currentBackend],
    Lookup[currentBackend, "match", None], None];
  If[StringQ[previousArm] && currentArm === previousArm &&
      IntegerQ[previousMatch] && IntegerQ[currentMatch] &&
      currentMatch > previousMatch,
    Max[current + additional, 2 current],
    current + additional]];

(* A terminal factorized observable can consume more epsilon orders than the
   ordinary one-row primitive.  That is a property of the producing level's
   matched functional, not of the lower consumer's matching solve.  Learn the
   observed producer loss separately; increasing the consumer halo here would
   raise the requested handoff edge and create a moving goalpost. *)
ft2NativeHandoffProducerRetry[lowerLevel_Integer, producerLevel_Integer,
    requiredTop_Integer, availableTop_Integer,
    producerSourceTop_Integer, baseIntrinsicLoss_Integer,
    currentProducerPrivateLoss_Integer] := Module[
  {observedLoss, requiredPrivateLoss, additional},
  If[lowerLevel < 1 || producerLevel =!= lowerLevel + 1 ||
      availableTop >= requiredTop, Return[None, Module]];
  observedLoss = producerSourceTop - availableTop;
  requiredPrivateLoss = Max[0, observedLoss - baseIntrinsicLoss];
  additional = requiredPrivateLoss - currentProducerPrivateLoss;
  If[observedLoss < 0 || additional < 1,
    Return[Failure["FeynmanTrickNativeProducerReservoirStalled", <|
      "Detail" ->
        "short native handoff did not imply a larger producer-loss bound",
      "Level" -> producerLevel, "ConsumerLevel" -> lowerLevel,
      "RequiredCompleteMax" -> requiredTop,
      "AvailableCompleteMax" -> availableTop,
      "ProducerSourceCompleteMax" -> producerSourceTop,
      "BaseIntrinsicLoss" -> baseIntrinsicLoss,
      "CurrentProducerPrivateLoss" -> currentProducerPrivateLoss|>],
      Module]];
  Failure["FeynmanTrickNativeProducerReservoir", <|
    "Detail" ->
      "native level handoff needs a wider producer epsilon reservoir",
    "Level" -> producerLevel, "ConsumerLevel" -> lowerLevel,
    "AdditionalOrders" -> additional,
    "RequiredCompleteMax" -> requiredTop,
    "AvailableCompleteMax" -> availableTop,
    "ProducerSourceCompleteMax" -> producerSourceTop,
    "ObservedProducerLoss" -> observedLoss,
    "BaseIntrinsicLoss" -> baseIntrinsicLoss,
    "CurrentProducerPrivateLoss" -> currentProducerPrivateLoss,
    "RequiredProducerPrivateLoss" -> requiredPrivateLoss|>]];

ft2NativeProducerReservoirRetryQ[failure_] := FailureQ[failure] &&
  Quiet[Check[failure[[1]] ===
    "FeynmanTrickNativeProducerReservoir", False]];

(* A complete epsilon residual with an inconclusive accuracy verdict is not a
   reservoir deficit.  A stored-Taylor residual directly motivates a larger
   Taylor work order.  An independently materialized continuity failure does
   not identify the missing resource by itself, but permits one bounded Taylor
   probe.  The progress guard below accepts only a later handoff failure,
   fewer inconclusive coefficients, or a material reduction in the same
   handoff's physical midpoint defect. *)
ft2NativeMatchingClearanceRetry[failure_?FailureQ,
    level_Integer, currentExpansionOrder_Integer] := Module[
  {data, backend, residual, verdicts, progress},
  data = Quiet[Check[failure[[2]], <||>]];
  backend = If[AssociationQ[data],
    Lookup[data, "BackendFailure", None], None];
  residual = If[AssociationQ[backend],
    Lookup[backend, "residual", None], None];
  verdicts = If[AssociationQ[residual],
    Lookup[residual, "required_coefficient_verdicts",
      Lookup[residual, "coefficient_verdicts", None]], None];
  progress = ft2NativeMatchingClearanceProgressRecord[backend, verdicts];
  If[currentExpansionOrder > 0 && AssociationQ[backend] &&
      AssociationQ[residual] &&
      Lookup[backend, "reason", None] ===
        "acb_match_residual_inconclusive" &&
      !TrueQ[Lookup[backend, "retryable_epsilon_reservoir", False]] &&
      TrueQ[Lookup[backend, "retryable_matching_clearance", False]] &&
      TrueQ[Lookup[residual, "complete_through_required", False]] &&
      MemberQ[{"stored-taylor-truncation",
          "materialized-continuity-clearance"},
        Lookup[residual, "scope", None]],
    Failure["FeynmanTrickNativeMatchingTaylor", <|
      "Detail" ->
        "native matching clearance permits one larger finite-Taylor overlap probe",
      "Level" -> level,
      "CurrentExpansionOrder" -> currentExpansionOrder,
      "AdditionalOrders" -> currentExpansionOrder,
      "ResidualVerdicts" -> verdicts,
      "ProgressRecord" -> progress,
      "BackendFailure" -> backend|>], None]];

ft2NativeMatchingClearanceRetry[_, _Integer, _] := None;

ft2NativeMatchingClearanceRetryQ[failure_] := FailureQ[failure] &&
  Quiet[Check[failure[[1]] ===
    "FeynmanTrickNativeMatchingTaylor", False]];

(* If the zero-radius physical midpoint passes but the full Acb ball does not,
   neither a wider epsilon frame nor more Taylor terms repairs the cause.  The
   uncertainty belongs to the boundary produced by the next-higher FT level.
   Preserve the consumer's requested tolerance and ask that producer for the
   standard safety margin instead. *)
ft2NativeMatchingProducerRetry[failure_?FailureQ, level_Integer,
    nLevels_Integer, currentProducerDigits_Integer] := Module[
  {data, backend, normalFrame, source, sourceProbes, incomingProbe,
   basisProbe, weightsProbe, upstreamWeightQ, additional},
  data = Quiet[Check[failure[[2]], <||>]];
  backend = If[AssociationQ[data],
    Lookup[data, "BackendFailure", None], None];
  normalFrame = If[AssociationQ[backend],
    Lookup[backend, "normal_frame_attempt", None], None];
  source = If[AssociationQ[normalFrame],
    Lookup[normalFrame, "physical_clearance_source", None], None];
  sourceProbes = If[AssociationQ[normalFrame],
    Lookup[normalFrame, "physical_clearance_source_probes", None], None];
  incomingProbe = If[AssociationQ[sourceProbes],
    Lookup[sourceProbes, "incoming", None], None];
  basisProbe = If[AssociationQ[sourceProbes],
    Lookup[sourceProbes, "basis", None], None];
  weightsProbe = If[AssociationQ[sourceProbes],
    Lookup[sourceProbes, "weights", None], None];
  upstreamWeightQ =
    AssociationQ[backend] &&
      Lookup[backend, "weight_shadow_retry", None] ===
        "upstream-accuracy-required-nonrational-chart" &&
      AssociationQ[basisProbe] &&
      Lookup[basisProbe, "verdict", None] === "pass" &&
      AssociationQ[weightsProbe] &&
      Lookup[weightsProbe, "verdict", None] === "inconclusive";
  additional = DiffExp2`Tolerances`$SafetyDigits;
  If[level < nLevels && currentProducerDigits > 0 && additional > 0 &&
      AssociationQ[backend] &&
      Lookup[backend, "reason", None] ===
        "acb_match_residual_inconclusive" &&
      !TrueQ[Lookup[backend, "retryable_epsilon_reservoir", False]] &&
      !TrueQ[Lookup[backend, "retryable_matching_clearance", False]] &&
      TrueQ[Lookup[backend, "retryable_propagated_enclosure", False]] &&
      source === "propagated-enclosure" &&
      ((AssociationQ[incomingProbe] &&
          Lookup[incomingProbe, "verdict", None] === "inconclusive") ||
        TrueQ[upstreamWeightQ]),
    Failure["FeynmanTrickNativeMatchingProducer", <|
      "Detail" ->
        "native matching needs a tighter immediately preceding boundary producer",
      "Level" -> level, "ProducerLevel" -> level + 1,
      "NumLevels" -> nLevels,
      "CurrentMatchingDigits" -> currentProducerDigits,
      "AdditionalOrders" -> additional,
      "BackendFailure" -> backend|>], None]];

ft2NativeMatchingProducerRetry[_, _Integer, _Integer, _Integer] := None;

(* A terminal boundary contraction is itself the producer for the next lower
   FT level.  If its rigorous output ball misses the requested publication
   width, raise this level and every upstream producer by the same guarded
   accuracy ladder used for propagated matching enclosures. *)
ft2NativeTerminalOutputProducerRetry[failure_?FailureQ,
    level_Integer, nLevels_Integer,
    currentLevelDigits_Integer] := Module[
  {data, backend, reportedAdditional, additional},
  data = Quiet[Check[failure[[2]], <||>]];
  backend = If[AssociationQ[data],
    Lookup[data, "BackendFailure", None], None];
  reportedAdditional = If[AssociationQ[backend],
    Lookup[backend, "required_additional_digits", None], None];
  additional = If[
    IntegerQ[reportedAdditional] && reportedAdditional > 0,
    reportedAdditional,
    DiffExp2`Tolerances`$SafetyDigits];
  If[level >= 1 && level <= nLevels && currentLevelDigits > 0 &&
      additional > 0 && AssociationQ[backend] &&
      Lookup[backend, "reason", None] ===
        "terminal_output_ball_inconclusive" &&
      TrueQ[Lookup[backend, "retryable_level_accuracy", False]],
    Failure["FeynmanTrickNativeMatchingProducer", <|
      "Detail" ->
        "native terminal boundary needs a tighter producer accuracy ladder",
      "Level" -> level - 1, "ProducerLevel" -> level,
      "NumLevels" -> nLevels,
      "CurrentMatchingDigits" -> currentLevelDigits,
      "AdditionalOrders" -> additional,
      "ProducerRetryKind" -> "terminal-output",
      "ProgressRecord" ->
        ft2NativeTerminalOutputProgressRecord[backend],
      "BackendFailure" -> backend|>], None]];

ft2NativeTerminalOutputProducerRetry[
    _, _Integer, _Integer, _Integer] := None;

ft2NativeMatchingProducerRetryQ[failure_] := FailureQ[failure] &&
  Quiet[Check[failure[[1]] ===
    "FeynmanTrickNativeMatchingProducer", False]];

ft2NativeTerminalOutputRadiusExponent[value_] := Module[{parts},
  parts = If[ListQ[value], value, {value}];
  parts = Replace[parts, {
      text_String /; StringMatchQ[text,
          ("-" | "+" | "") ~~ DigitCharacter..] :>
        Switch[StringTake[text, 1],
          "-", -FromDigits[StringDrop[text, 1]],
          "+", FromDigits[StringDrop[text, 1]],
          _, FromDigits[text]],
      "zero" -> -Infinity
    }, {1}];
  If[parts === {} || !AllTrue[parts, NumericQ], None, Max[parts]]
];

(* Producer precision is a useful retry resource only if it actually tightens
   the first unpublished terminal coefficient.  Track the evaluator's ordered
   failure position and its binary radius.  A later coefficient, or a smaller
   radius at the same coefficient, is progress; an unchanged or wider ball is
   a structural floor and must not consume the whole producer-digit budget. *)
ft2NativeTerminalOutputProgressRecord[backend_] := Module[
  {conditioning, entries, functional, epsilonPower, entry, radius},
  If[!AssociationQ[backend], Return[None, Module]];
  functional = Lookup[backend, "failure_functional", None];
  epsilonPower = Lookup[backend, "failure_epsilon", None];
  conditioning = Lookup[backend, "conditioning", <||>];
  entries = If[AssociationQ[conditioning],
    Lookup[conditioning, "entries", {}], {}];
  entry = If[ListQ[entries],
    SelectFirst[entries,
      AssociationQ[#] &&
        Lookup[#, "power", None] === epsilonPower &&
        Lookup[#, "component", None] === functional &, None],
    None];
  radius = If[AssociationQ[entry],
    ft2NativeTerminalOutputRadiusExponent[
      Lookup[entry, "combined_radius2exp", None]], None];
  <|
    "Scope" -> Lookup[backend, "scope", None],
    "Functional" -> functional,
    "EpsilonPower" -> epsilonPower,
    "CombinedRadius2Exponent" -> radius
  |>
];

ft2NativeTerminalOutputProgressQ[previous_, current_] := Module[
  {previousPosition, currentPosition, previousRadius, currentRadius},
  If[!AssociationQ[previous] || !AssociationQ[current] ||
      Lookup[previous, "Scope", None] =!=
        Lookup[current, "Scope", None],
    Return[False, Module]];
  previousPosition =
    Lookup[previous, {"Functional", "EpsilonPower"}, None];
  currentPosition =
    Lookup[current, {"Functional", "EpsilonPower"}, None];
  If[MatchQ[previousPosition, {_Integer, _Integer}] &&
      MatchQ[currentPosition, {_Integer, _Integer}] &&
      (currentPosition[[1]] > previousPosition[[1]] ||
        (currentPosition[[1]] === previousPosition[[1]] &&
          currentPosition[[2]] > previousPosition[[2]])),
    Return[True, Module]];
  previousRadius =
    Lookup[previous, "CombinedRadius2Exponent", None];
  currentRadius =
    Lookup[current, "CombinedRadius2Exponent", None];
  currentPosition === previousPosition &&
    NumericQ[previousRadius] && NumericQ[currentRadius] &&
    TrueQ[currentRadius < previousRadius]
];

(* A matching residual can require a guarded accuracy ladder: an L-level
   producer requested at D digits must itself receive a boundary with the
   safety margin still intact.  Terminal-output publication is different:
   advance one producer seam at a time.  If that producer cannot publish the
   tighter boundary, its own typed terminal failure identifies the next seam.
   Blindly raising the entire ladder can activate wider non-monotone
   enclosures before the first requested producer has been tested. *)
ft2RaiseMatchingProducerDigits[current_Association,
    producerLevel_Integer, requestedDigits_Integer,
    nLevels_Integer] := Module[{guard, raised},
  guard = DiffExp2`Tolerances`$SafetyDigits;
  If[producerLevel < 1 || producerLevel > nLevels ||
      requestedDigits < matchDigits || guard < 1,
    Return[$Failed, Module]];
  raised = Association[Table[
    upperLevel -> requestedDigits +
      (upperLevel - producerLevel) guard,
    {upperLevel, producerLevel, nLevels}]];
  Merge[{current, raised}, Max]];

(* The producer budget must scale with ladder depth.  A fixed allowance of
   four safety guards is not enough for a seven-level ladder: even one
   downstream retry has to preserve a guard at every upstream handoff.  The
   retry loop remains bounded separately by maxAttempts, while this limit
   permits one full guarded correction beyond the deepest possible chain. *)
ft2MatchingProducerDigitLimit[nLevels_Integer] := Module[{guard},
  guard = DiffExp2`Tolerances`$SafetyDigits;
  If[nLevels < 1 || guard < 1 || !IntegerQ[wp] || !IntegerQ[matchDigits],
    Return[$Failed, Module]];
  Min[wp, matchDigits + (nLevels + 1) guard]];

ft2NativeMatchingRetryQ[failure_] :=
  ft2NativeMatchingReservoirRetryQ[failure] ||
  ft2NativeProducerReservoirRetryQ[failure] ||
  ft2NativeMatchingClearanceRetryQ[failure] ||
  ft2NativeMatchingProducerRetryQ[failure];

ft2NativeMatchingClearanceProgressRecord[backend_, verdicts_] := Module[
  {normalFrame, midpointResidual, ratio, geometry},
  normalFrame = If[AssociationQ[backend],
    Lookup[backend, "normal_frame_attempt", backend], <||>];
  midpointResidual = If[AssociationQ[normalFrame],
    Lookup[normalFrame, "physical_midpoint_residual", <||>], <||>];
  ratio = If[AssociationQ[midpointResidual],
    Lookup[midpointResidual,
      "maximum_required_residual_to_scale_upper_approx", None], None];
  If[!NumericQ[ratio] || !TrueQ[ratio >= 0], ratio = None];
  geometry = If[AssociationQ[backend],
    Lookup[backend, "geometry", <||>], <||>];
  <|
    "ResidualVerdicts" -> verdicts,
    "Arm" -> If[AssociationQ[backend],
      Lookup[backend, "arm", None], None],
    "Match" -> If[AssociationQ[backend],
      Lookup[backend, "match", None], None],
    "PhysicalPoint" -> If[AssociationQ[geometry],
      Lookup[geometry, "physical_exact", None], None],
    "MidpointResidualRatio" -> ratio
  |>
];

ft2NativeMatchingClearancePosition[record_] := Module[
  {arm, match, rank},
  If[!AssociationQ[record], Return[None, Module]];
  arm = Lookup[record, "Arm", None];
  match = Lookup[record, "Match", None];
  rank = Switch[arm, "lower", 0, "upper", 1, _, None];
  If[IntegerQ[rank] && IntegerQ[match], {rank, match}, None]
];

ft2NativeMatchingClearanceProgressQ[previous_, current_] := Module[
  {previousVerdicts, currentVerdicts, previousInconclusive,
   currentInconclusive, previousPosition, currentPosition,
   previousRatio, currentRatio},
  If[!AssociationQ[previous] || !AssociationQ[current],
    Return[False, Module]];
  previousVerdicts = Lookup[
    previous, "ResidualVerdicts", previous];
  currentVerdicts = Lookup[
    current, "ResidualVerdicts", current];
  previousInconclusive = If[AssociationQ[previousVerdicts],
    Lookup[previousVerdicts, "inconclusive", None], None];
  currentInconclusive = If[AssociationQ[currentVerdicts],
    Lookup[currentVerdicts, "inconclusive", None], None];
  If[IntegerQ[previousInconclusive] &&
      IntegerQ[currentInconclusive] &&
      currentInconclusive < previousInconclusive,
    Return[True, Module]];
  previousPosition =
    ft2NativeMatchingClearancePosition[previous];
  currentPosition =
    ft2NativeMatchingClearancePosition[current];
  If[ListQ[previousPosition] && ListQ[currentPosition] &&
      (currentPosition[[1]] > previousPosition[[1]] ||
        (currentPosition[[1]] === previousPosition[[1]] &&
          currentPosition[[2]] > previousPosition[[2]])),
    Return[True, Module]];
  previousRatio = Lookup[
    previous, "MidpointResidualRatio", None];
  currentRatio = Lookup[
    current, "MidpointResidualRatio", None];
  ListQ[previousPosition] && currentPosition === previousPosition &&
    NumericQ[previousRatio] && NumericQ[currentRatio] &&
    TrueQ[previousRatio > 0] &&
    TrueQ[currentRatio < 99 previousRatio/100]
];

ft2CanonicalIdentity[prefix_String, value_] := prefix <>
  IntegerString[Hash[value, "SHA256"], 16, 64];

(* Matching-reservoir deficits are structural lower bounds for one exact
   prepared ladder/configuration.  Persist them separately from numerical
   Boundary and NativeTransport checkpoints: a learned profile may seed a
   new plan, but it is never evidence that numerical state can be restored.
   The contract includes the zero-private-halo base plan, so a profile cannot
   cross a changed matrix, FIRE batch, public epsilon request, or level loss.
   It deliberately excludes the whole-tree source fingerprint: a learned
   lower bound remains safe across unrelated implementation edits because
   every new run still performs the current matching certificates and may
   monotonically discover more clearance. *)
$ft2MatchingHaloProfileMax = 64;
ft2NormalizeMatchingHaloBounds[bounds_, nLevels_Integer] := Module[{values},
  If[nLevels < 1, Return[Failure["FeynmanTrickMatchingHaloProfile", <|
    "Detail" -> "matching-halo profile requires a positive level count"|>],
    Module]];
  values = Which[
    AssociationQ[bounds] && AllTrue[Keys[bounds],
        IntegerQ[#] && 1 <= # <= nLevels &],
      Lookup[bounds, Range[nLevels], 0],
    ListQ[bounds] && Length[bounds] === nLevels, bounds,
    bounds === Automatic || bounds === None, ConstantArray[0, nLevels],
    True, Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "matching-halo bounds must be a level association or full list",
      "Bounds" -> bounds, "NumLevels" -> nLevels|>], Module]];
  If[!AllTrue[values, IntegerQ[#] &&
      0 <= # <= $ft2MatchingHaloProfileMax &],
    Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "matching-halo bounds must be bounded nonnegative integers",
      "Bounds" -> values|>], Module]];
  values];

ft2MergeMatchingHaloBounds[nLevels_Integer, candidates_List] := Module[
  {normalized},
  normalized = ft2NormalizeMatchingHaloBounds[#, nLevels] & /@ candidates;
  If[AnyTrue[normalized, FailureQ],
    Return[First[Select[normalized, FailureQ]], Module]];
  Association@MapIndexed[First[#2] -> #1 &,
    Apply[Max, Transpose[normalized], {1}]]];

ft2MatchingHaloProfileContract[name_String, prepKey_,
    basePlan_Association, anchors_List:{11/23}] := Module[
  {record, identity},
  If[!ft2NativeEpsilonPlanQ[basePlan] ||
      basePlan["Record", "MatchingPrivateHalos"] =!=
        ConstantArray[0, basePlan["NumLevels"]] ||
      basePlan["Record", "ProducerPrivateLosses"] =!=
        ConstantArray[0, basePlan["NumLevels"]],
    Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "matching-halo profile contract requires the exact zero-private-halo base plan"|>],
    Module]];
  record = Join[<|
    "Schema" -> "FeynmanTrick.MatchingHaloProfileContract/v3",
    "Example" -> name,
    "PrepKey" -> prepKey,
    "MatchingProofPolicy" ->
      "exact-match-reservoir-and-staged-terminal-producer-v3",
    "BasePlanIdentity" -> basePlan["Identity"],
    "BasePlanRecord" -> basePlan["Record"],
    "Configuration" -> <|
      "WorkingPrecision" -> wp,
      "MatchingDigits" -> matchDigits,
      "MatchingCertificationSafetyDigits" ->
        matchingCertificationSafetyDigits,
      "MatchingCertificationMaximumDigits" ->
        matchingCertificationMaximumDigits,
      "MatchingCertificationDigitsByLevel" ->
        matchingCertificationDigitsByLevel,
      "MatchingCertificationComputationDigitsByLevel" ->
        matchingCertificationComputationDigitsByLevel,
      "ExpansionOrder" -> expansionOrder,
      "EpsilonOrder" -> epsOrder,
      "BoundaryExtraOrder" -> boundaryExtraOrder,
      "DivisionOrder" -> divisionOrder,
      "StepDivisionOrder" -> stepDivisionOrder,
      "RadiusOfConvergence" -> radiusOfConvergence,
      "RecurrenceBackend" -> recurrenceBackend,
      "SingularMatchPrecondition" -> singularMatchPrecondition,
      "DeltaPrescriptionSign" -> deltaPrescriptionSign,
      "LevelDeltaPrescriptionSigns" -> levelDeltaPrescriptionSigns,
      "ValueTransportMode" -> valueTransportMode,
      "NativeValueHopExecutionMode" -> nativeValueHopExecutionMode,
      "CppThreads" -> cppArmThreadBudget,
      "MatchingTaylorProgressPolicy" ->
        "later-position-or-fewer-inconclusive-or-1pct-midpoint-v1",
      "Anchor" -> anchors|>|>, ft2CheckpointRequestMetadata[]];
  identity = ft2CanonicalIdentity[
    "ft2-matching-halo-profile-contract-", record];
  <|"Record" -> record, "Identity" -> identity,
    "NumLevels" -> basePlan["NumLevels"]|>];

ft2MatchingHaloProfileContractQ[contract_] := AssociationQ[contract] &&
  Sort[Keys[contract]] === Sort[{"Record", "Identity", "NumLevels"}] &&
  AssociationQ[contract["Record"]] &&
  Lookup[contract["Record"], "Schema", None] ===
    "FeynmanTrick.MatchingHaloProfileContract/v3" &&
  Lookup[contract["Record"], "MatchingProofPolicy", None] ===
    "exact-match-reservoir-and-staged-terminal-producer-v3" &&
  IntegerQ[contract["NumLevels"]] && contract["NumLevels"] >= 1 &&
  With[{base = Lookup[contract["Record"], "BasePlanRecord", <||>]},
    AssociationQ[base] &&
      Lookup[base, "Schema", None] ===
        "FeynmanTrick.NativeEpsilonPlan/v1" &&
      Lookup[base, "NumLevels", None] === contract["NumLevels"] &&
      ListQ[Lookup[base, "Levels", None]] &&
      Length[base["Levels"]] === contract["NumLevels"] &&
      Lookup[base, "MatchingPrivateHalos", None] ===
        ConstantArray[0, contract["NumLevels"]] &&
      Lookup[base, "ProducerPrivateLosses", None] ===
        ConstantArray[0, contract["NumLevels"]]] &&
  Lookup[contract["Record"], "BasePlanIdentity", None] ===
    ft2CanonicalIdentity["ft2-native-epsilon-plan-",
      Lookup[contract["Record"], "BasePlanRecord", None]] &&
  contract["Identity"] === ft2CanonicalIdentity[
    "ft2-matching-halo-profile-contract-", contract["Record"]];

ft2MatchingHaloProfileFile[contract_Association] :=
  FileNameJoin[{matchingHaloProfileRoot,
    "matching_halos_" <> StringTake[contract["Identity"], -64] <> ".mx"}];

ft2MatchingHaloProfileQ[profile_, contract_Association] := Module[
  {core, schema, expectedKeys},
  schema = If[AssociationQ[profile],
    Lookup[profile, "Schema", None], None];
  expectedKeys = Switch[schema,
    "FeynmanTrick.MatchingHaloProfile/v2",
      {"Schema", "Contract", "ContractIdentity", "NumLevels",
        "MatchingPrivateHalos", "ProducerPrivateLosses", "Identity"},
    "FeynmanTrick.MatchingHaloProfile/v3",
      {"Schema", "Contract", "ContractIdentity", "NumLevels",
        "MatchingPrivateHalos", "ProducerPrivateLosses",
        "ProducerMatchingDigitExtras", "Identity"},
    _, {}];
  If[!ft2MatchingHaloProfileContractQ[contract] ||
      !AssociationQ[profile] ||
      expectedKeys === {} ||
      Sort[Keys[profile]] =!= Sort[expectedKeys],
    Return[False, Module]];
  core = KeyDrop[profile, "Identity"];
  TrueQ[MemberQ[{"FeynmanTrick.MatchingHaloProfile/v2",
        "FeynmanTrick.MatchingHaloProfile/v3"}, schema] &&
    profile["Contract"] === contract["Record"] &&
    profile["ContractIdentity"] === contract["Identity"] &&
    profile["NumLevels"] === contract["NumLevels"] &&
    ListQ[profile["MatchingPrivateHalos"]] &&
    Length[profile["MatchingPrivateHalos"]] === contract["NumLevels"] &&
    AllTrue[profile["MatchingPrivateHalos"],
      IntegerQ[#] && 0 <= # <= $ft2MatchingHaloProfileMax &] &&
    ListQ[profile["ProducerPrivateLosses"]] &&
    Length[profile["ProducerPrivateLosses"]] === contract["NumLevels"] &&
    AllTrue[profile["ProducerPrivateLosses"],
      IntegerQ[#] && 0 <= # <= $ft2MatchingHaloProfileMax &] &&
    (schema === "FeynmanTrick.MatchingHaloProfile/v2" ||
      (ListQ[profile["ProducerMatchingDigitExtras"]] &&
       Length[profile["ProducerMatchingDigitExtras"]] ===
         contract["NumLevels"] &&
       AllTrue[profile["ProducerMatchingDigitExtras"],
         IntegerQ[#] && 0 <= # <= $ft2MatchingHaloProfileMax &])) &&
    profile["Identity"] === ft2CanonicalIdentity[
      "ft2-matching-halo-profile-", core]]];

ft2LoadMatchingHaloProfile[file_String, contract_Association] := Module[
  {ok, profile},
  If[!FileExistsQ[file], Return[ConstantArray[0,
    contract["NumLevels"]], Module]];
  Clear[Global`$FT2MatchingHaloProfile];
  ok = Quiet[Check[Get[file]; True, False]];
  profile = If[TrueQ[ok], Global`$FT2MatchingHaloProfile, None];
  Clear[Global`$FT2MatchingHaloProfile];
  If[!TrueQ[ok] || !ft2MatchingHaloProfileQ[profile, contract],
    Print["FTLADDER MATCH PROFILE REJECT ", file,
      ": malformed, tampered, or contract-mismatched"];
    Return[ConstantArray[0, contract["NumLevels"]], Module]];
  Print["FTLADDER MATCH PROFILE HIT ", file,
    " privateHalos=", profile["MatchingPrivateHalos"]];
  profile["MatchingPrivateHalos"]];

ft2LoadProducerLossProfile[file_String, contract_Association] := Module[
  {ok, profile},
  If[!FileExistsQ[file], Return[ConstantArray[0,
    contract["NumLevels"]], Module]];
  Clear[Global`$FT2MatchingHaloProfile];
  ok = Quiet[Check[Get[file]; True, False]];
  profile = If[TrueQ[ok], Global`$FT2MatchingHaloProfile, None];
  Clear[Global`$FT2MatchingHaloProfile];
  If[!TrueQ[ok] || !ft2MatchingHaloProfileQ[profile, contract],
    Return[ConstantArray[0, contract["NumLevels"]], Module]];
  profile["ProducerPrivateLosses"]];

ft2LoadProducerDigitProfile[file_String, contract_Association] := Module[
  {ok, profile},
  If[!FileExistsQ[file], Return[ConstantArray[0,
    contract["NumLevels"]], Module]];
  Clear[Global`$FT2MatchingHaloProfile];
  ok = Quiet[Check[Get[file]; True, False]];
  profile = If[TrueQ[ok], Global`$FT2MatchingHaloProfile, None];
  Clear[Global`$FT2MatchingHaloProfile];
  If[!TrueQ[ok] || !ft2MatchingHaloProfileQ[profile, contract],
    Return[ConstantArray[0, contract["NumLevels"]], Module]];
  Lookup[profile, "ProducerMatchingDigitExtras",
    ConstantArray[0, contract["NumLevels"]]]];

ft2AcquireMatchingHaloProfileLock[file_String] := Module[
  {lock = file <> ".lock", acquired, age, attempts = 200},
  If[!DirectoryQ[DirectoryName[file]],
    Quiet[Check[CreateDirectory[DirectoryName[file],
      CreateIntermediateDirectories -> True], Null]]];
  Do[
    acquired = Quiet[Check[CreateDirectory[lock]; True, False]];
    If[TrueQ[acquired], Return[lock, Module]];
    If[DirectoryQ[lock],
      age = Quiet[Check[AbsoluteTime[] - AbsoluteTime[FileDate[lock]], 0]];
      (* A profile write is a tiny local DumpSave plus rename.  Ten minutes
         cannot be a live writer; it is a lock left by a killed process. *)
      If[NumericQ[age] && age > 600,
        Quiet[Check[DeleteDirectory[lock, DeleteContents -> True], Null]]]];
    Pause[0.05],
    {attempts}];
  Failure["FeynmanTrickMatchingHaloProfile", <|
    "Detail" -> "timed out acquiring the matching-halo profile lock",
    "File" -> file, "Lock" -> lock|>]];

ft2SaveMatchingHaloProfileUnlocked[file_String, contract_Association,
    bounds_, producerLossBounds_:Automatic,
    producerDigitBounds_:Automatic] := Module[
  {existing, existingProducer, existingProducerDigits, merged,
   mergedProducer, mergedProducerDigits, core, profile, tmp, wrote, renamed},
  existing = ft2LoadMatchingHaloProfile[file, contract];
  existingProducer = ft2LoadProducerLossProfile[file, contract];
  existingProducerDigits = ft2LoadProducerDigitProfile[file, contract];
  merged = ft2MergeMatchingHaloBounds[contract["NumLevels"],
    {existing, bounds}];
  If[FailureQ[merged], Return[merged, Module]];
  mergedProducer = ft2MergeMatchingHaloBounds[contract["NumLevels"],
    {existingProducer, producerLossBounds}];
  If[FailureQ[mergedProducer], Return[mergedProducer, Module]];
  mergedProducerDigits = ft2MergeMatchingHaloBounds[
    contract["NumLevels"],
    {existingProducerDigits, producerDigitBounds}];
  If[FailureQ[mergedProducerDigits],
    Return[mergedProducerDigits, Module]];
  core = <|
    "Schema" -> "FeynmanTrick.MatchingHaloProfile/v3",
    "Contract" -> contract["Record"],
    "ContractIdentity" -> contract["Identity"],
    "NumLevels" -> contract["NumLevels"],
    "MatchingPrivateHalos" -> Lookup[merged,
      Range[contract["NumLevels"]], 0],
    "ProducerPrivateLosses" -> Lookup[mergedProducer,
      Range[contract["NumLevels"]], 0],
    "ProducerMatchingDigitExtras" -> Lookup[mergedProducerDigits,
      Range[contract["NumLevels"]], 0]|>;
  profile = Append[core, "Identity" -> ft2CanonicalIdentity[
    "ft2-matching-halo-profile-", core]];
  If[!DirectoryQ[DirectoryName[file]],
    Quiet[Check[CreateDirectory[DirectoryName[file],
      CreateIntermediateDirectories -> True], Null]]];
  tmp = file <> ".tmp-" <> ToString[$ProcessID] <> "-" <>
    IntegerString[Hash[{AbsoluteTime[], RandomInteger[]}], 16] <> ".mx";
  If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
  Global`$FT2MatchingHaloProfile = profile;
  wrote = Quiet[Check[
    DumpSave[tmp, Global`$FT2MatchingHaloProfile]; FileExistsQ[tmp], False]];
  Clear[Global`$FT2MatchingHaloProfile];
  If[!TrueQ[wrote],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "could not write the temporary matching-halo profile",
      "File" -> file|>], Module]];
  renamed = Quiet[Check[
    RenameFile[tmp, file, OverwriteTarget -> True]; True, False]];
  If[!TrueQ[renamed],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "could not atomically publish the matching-halo profile",
      "File" -> file|>], Module]];
  Print["FTLADDER MATCH PROFILE ", file,
    " privateHalos=", profile["MatchingPrivateHalos"],
    " producerLosses=", profile["ProducerPrivateLosses"],
    " producerDigitExtras=", profile["ProducerMatchingDigitExtras"]];
  profile];

ft2SaveMatchingHaloProfile[file_String, contract_Association,
    bounds_, producerLossBounds_:Automatic,
    producerDigitBounds_:Automatic] := Module[{lock},
  If[!ft2MatchingHaloProfileContractQ[contract],
    Return[Failure["FeynmanTrickMatchingHaloProfile", <|
      "Detail" -> "cannot save a matching-halo profile under a malformed contract"|>],
    Module]];
  lock = ft2AcquireMatchingHaloProfileLock[file];
  If[FailureQ[lock], Return[lock, Module]];
  Internal`WithLocalSettings[
    Null,
    ft2SaveMatchingHaloProfileUnlocked[
      file, contract, bounds, producerLossBounds, producerDigitBounds],
    If[DirectoryQ[lock], Quiet[Check[
      DeleteDirectory[lock, DeleteContents -> True], Null]]]]];

ft2CanonicalBatchKeyQ[value_] := StringQ[value] &&
  StringMatchQ[value, RegularExpression["[0-9a-f]{64}"]];

ft2NativeCheckpointRecordQ[record_] := AssociationQ[record] &&
  Lookup[record, "Schema", None] ===
    "FeynmanTrick.NativeObservableBatch/v1" &&
  ft2CanonicalBatchKeyQ[Lookup[record, "BatchKey", None]] &&
  ft2CanonicalBatchKeyQ[Lookup[record, "BatchPayloadKey", None]] &&
  ListQ[Lookup[record, "RequestIdentities", None]] &&
  ListQ[Lookup[record, "CoefficientIdentities", None]] &&
  Length[record["RequestIdentities"]] ===
    Length[record["CoefficientIdentities"]] &&
  AllTrue[Join[record["RequestIdentities"],
      record["CoefficientIdentities"],
      Lookup[record, "ObservableIdentities", {}],
      Lookup[record, "ObservableCheckpointIdentities", {}],
      {Lookup[record, "DeltaPrescriptionIdentity", None],
       Lookup[record, "ExtraSingularFactorsIdentity", None],
       Lookup[record, "NativeBatchPayloadIdentity", None]}],
    StringQ[#] && StringLength[#] > 0 &] &&
  ListQ[Lookup[record, "DeltaPrescriptions", None]] &&
  Lookup[record, "DeltaPrescriptionIdentity", None] ===
    ft2CanonicalIdentity["ft2-delta-prescriptions-",
      record["DeltaPrescriptions"]] &&
  (Lookup[record, "AtlasPlanIdentity", None] === None ||
    (StringQ[record["AtlasPlanIdentity"]] &&
      StringLength[record["AtlasPlanIdentity"]] > 0)) &&
  (Lookup[record, "NativeEpsilonPlanIdentity", None] === None ||
    (StringQ[record["NativeEpsilonPlanIdentity"]] &&
      StringLength[record["NativeEpsilonPlanIdentity"]] > 0)) &&
  AllTrue[Lookup[record, {"SourceCompleteMax", "TargetCompleteMax",
      "RequiredTargetCompleteMax",
      "DeliverableCompleteMax", "RequiredRawTop",
      "CoefficientHalo", "IntegrationHalo",
      "MatchEpsilonPadding"},
    None], IntegerQ] &&
  With[{source = record["SourceCompleteMax"],
      target = record["TargetCompleteMax"],
      requiredTarget = record["RequiredTargetCompleteMax"],
      deliverable = record["DeliverableCompleteMax"],
      required = record["RequiredRawTop"],
      coefficientHalo = record["CoefficientHalo"],
      integrationHalo = record["IntegrationHalo"],
      matchPadding = record["MatchEpsilonPadding"],
      hasRequiredSolve =
        KeyExistsQ[record, "RequiredSolveCompleteMax"],
      hasMaximumSourceLoss =
        KeyExistsQ[record, "MaximumSourceLoss"]},
    coefficientHalo >= 0 && matchPadding >= 0 &&
      MemberQ[{0, 1}, integrationHalo] &&
      source >= target >= deliverable >= required &&
      target >= requiredTarget &&
      If[hasRequiredSolve && hasMaximumSourceLoss,
        With[{requiredSolve = record["RequiredSolveCompleteMax"],
            maximumSourceLoss = record["MaximumSourceLoss"]},
          IntegerQ[requiredSolve] && IntegerQ[maximumSourceLoss] &&
            maximumSourceLoss >= 0 &&
            requiredSolve === requiredTarget + maximumSourceLoss &&
            requiredSolve + matchPadding === source &&
            source >= requiredSolve],
        !hasRequiredSolve && !hasMaximumSourceLoss &&
          requiredTarget + matchPadding === source - coefficientHalo]];

ft2NativeTransportContract[name_String, level_Integer, prepKey_, sys_,
    boundaryValues_, boundaryPrefactors_, entries_List, ledger_Association,
    configuration_Association, deltaPrescriptions_List,
    extraSingularFactors_List, nativePlanIdentity_String] := Module[
  {record},
  record = Join[<|
    "Schema" -> "FeynmanTrick.NativeTransportContract/v1",
    "SourceFingerprint" -> $ftLadderSourceFingerprint,
    "Example" -> name, "Level" -> level, "PrepKey" -> prepKey,
    "SystemIdentity" -> ft2CanonicalIdentity[
      "ft2-native-system-", sys],
    "BoundaryIdentity" -> ft2CanonicalIdentity[
      "ft2-native-boundary-", {boundaryValues, boundaryPrefactors}],
    "BatchIdentity" -> ft2CanonicalIdentity[
      "ft2-native-batch-", Map[KeyTake[#,
        {"MasterIndex", "Case", "RequestIdentity", "CoefficientIdentity",
          "Identity", "CheckpointIdentity"}] &, entries]],
    "LedgerIdentity" -> ft2CanonicalIdentity[
      "ft2-native-ledger-", ledger],
    "ConfigurationIdentity" -> ft2CanonicalIdentity[
      "ft2-native-configuration-", configuration],
    "BranchIdentity" -> ft2CanonicalIdentity[
      "ft2-native-branch-",
      {deltaPrescriptions, extraSingularFactors}],
    "NativeEpsilonPlanIdentity" -> nativePlanIdentity|>,
    ft2CheckpointRequestMetadata[]];
  <|"Record" -> record, "Identity" -> ft2CanonicalIdentity[
    "ft2-native-transport-contract-", record]|>];

ft2NativeTransportContractQ[contract_] := AssociationQ[contract] &&
  Sort[Keys[contract]] === {"Identity", "Record"} &&
  AssociationQ[contract["Record"]] &&
  Lookup[contract["Record"], "Schema", None] ===
    "FeynmanTrick.NativeTransportContract/v1" &&
  StringQ[Lookup[contract["Record"], "NativeEpsilonPlanIdentity", None]] &&
  contract["Identity"] === ft2CanonicalIdentity[
    "ft2-native-transport-contract-", contract["Record"]];

ft2NativeTransportResumeRecordQ[record_] := Module[{core},
  If[!AssociationQ[record] || Sort[Keys[record]] =!= Sort[{
      "Schema", "ContractIdentity", "AtlasPlanIdentity",
      "NativeBatchPayloadIdentity", "CheckpointIdentity",
      "State", "Identity"}], Return[False, Module]];
  core = KeyDrop[record, "Identity"];
  TrueQ[record["Schema"] ===
      "FeynmanTrick.NativeTransportResume/v1"] &&
    AllTrue[Lookup[record, {"ContractIdentity", "AtlasPlanIdentity",
        "NativeBatchPayloadIdentity", "CheckpointIdentity"}],
      StringQ[#] && StringLength[#] > 0 &] &&
    AssociationQ[record["State"]] &&
    Lookup[record["State"], "CheckpointIdentity", None] ===
      record["CheckpointIdentity"] &&
    record["Identity"] === ft2CanonicalIdentity[
      "ft2-native-transport-resume-", core]];

ft2ExactEpsilonValuation[expression_, physicalVar_Symbol,
    epsSymbol_Symbol] := Module[
  {canonical = Together[expression], numerator, denominator,
   numeratorValuation, denominatorValuation},
  If[!FreeQ[canonical, _?InexactNumberQ],
    Return[ft2NativeFailure["native boundary coefficients must be exact",
      <|"Expression" -> canonical|>], Module]];
  If[TrueQ[PossibleZeroQ[canonical]], Return[0, Module]];
  numerator = Numerator[canonical];
  denominator = Denominator[canonical];
  If[!PolynomialQ[numerator, {physicalVar, epsSymbol}] ||
      !PolynomialQ[denominator, {physicalVar, epsSymbol}],
    Return[ft2NativeFailure[
      "native boundary coefficients must be rational in the Feynman parameter and epsilon",
      <|"Expression" -> canonical|>], Module]];
  numeratorValuation = Exponent[numerator, epsSymbol, Min];
  denominatorValuation = Exponent[denominator, epsSymbol, Min];
  If[!IntegerQ[numeratorValuation] || !IntegerQ[denominatorValuation],
    Return[ft2NativeFailure[
      "could not determine an exact integer epsilon valuation",
      <|"Expression" -> canonical|>], Module]];
  numeratorValuation - denominatorValuation];

(* Exact native primitive bound: a monomial chart primitive can acquire only
   one epsilon pole, at alpha0==0 on one center endpoint.  Normalized log
   chains still begin at -1, while paired m==-1 definite primitives cancel
   that pole.  Coefficient poles are accounted independently by HC. *)
ft2NativeIntegrationHalo[entries_List] := If[
  AnyTrue[entries, Lookup[#, "Case", None] === "integrate" &], 1, 0];

ft2PrepareBoundaryEntries[level_Integer, batch_Association,
    currentPrefactors_List, physicalVar_Symbol, epsSymbol_Symbol,
    normalize_, epsilonGauge_:Automatic] := Module[
  {requests = Lookup[batch, "BoundaryRequests", None],
   vectors = Lookup[batch, "CoefficientVectors", None],
   batchKey = Lookup[batch, "Key", Missing["NoBatchKey"]],
   payloadKey = Lookup[batch, "PayloadKey", Missing["NoPayloadKey"]],
   dimension = Length[currentPrefactors], entries,
   physicalFromRelative},
  physicalFromRelative = If[
    AssociationQ[epsilonGauge] &&
      TrueQ[Lookup[epsilonGauge, "CompositeApparent", False]],
    Lookup[epsilonGauge, "PhysicalFromRelative", None], None];
  If[physicalFromRelative =!= None &&
      (!MatrixQ[physicalFromRelative] ||
        Dimensions[physicalFromRelative] =!= {dimension, dimension}),
    Return[ft2NativeFailure[
      "composite apparent coefficient transform has the wrong dimension",
      <|"Expected" -> {dimension, dimension},
        "Actual" -> Quiet[Check[
          Dimensions[physicalFromRelative], Missing["Malformed"]]]|>],
      Module]];
  If[Lookup[batch, "Schema", None] =!=
        "FeynmanTrick.LevelIBPBatch/v1" ||
      Lookup[batch, "UpperLevel", None] =!= level ||
      !ft2CanonicalBatchKeyQ[batchKey] ||
      !ft2CanonicalBatchKeyQ[payloadKey] ||
      !ListQ[Lookup[batch, "KeyRecord", None]] ||
      !ListQ[Lookup[batch, "MastersAbove", None]] ||
      Length[batch["MastersAbove"]] =!= dimension ||
      !ListQ[requests] || requests === {} ||
      !AssociationQ[vectors] || dimension === 0 ||
      Lookup[requests, "MasterIndex", {}] =!= Range[Length[requests]],
    Return[ft2NativeFailure[
      "level IBP batch schema, keys, master indices, or coefficient vectors are invalid",
      <|"Schema" -> Lookup[batch, "Schema", None],
        "UpperLevel" -> Lookup[batch, "UpperLevel", None],
        "BatchKey" -> batchKey, "PayloadKey" -> payloadKey,
        "MasterIndices" -> If[ListQ[requests],
          Lookup[requests, "MasterIndex", {}], None]|>], Module]];
  entries = MapIndexed[Function[{request, position}, Module[
      {masterIndex = First[position], needed, raw, base, coefficients,
       nonzeroCoefficients, provenZero, shifts, case, expectedCase, vi,
       vj, requestIdentity, coefficientIdentity, identityHash},
      needed = Lookup[request, "NeededVec", Missing["NoNeededVector"]];
      case = Lookup[request, "Case", None];
      vi = Lookup[request, "Vi", None];
      vj = Lookup[request, "Vj", None];
      expectedCase = If[IntegerQ[vi] && IntegerQ[vj], Which[
        vi > 0 && vj > 0, "integrate",
        vi > 0 && vj === 0, "limitUpper",
        vi === 0 && vj > 0, "limitLower",
        True, "direct"], None];
      If[!MemberQ[{"integrate", "limitLower", "limitUpper", "direct"},
          case] || case =!= expectedCase ||
          Lookup[request, "MasterIndex", None] =!= masterIndex ||
          !ListQ[Lookup[request, "MasterVec", None]] ||
          !ListQ[needed] ||
          !KeyExistsQ[vectors, needed],
        Return[ft2NativeFailure["malformed batched boundary request",
          <|"MasterIndex" -> masterIndex, "Request" -> request|>], Module]];
      raw = vectors[needed];
      If[!ListQ[raw] || Length[raw] =!= dimension,
        Return[ft2NativeFailure[
          "batched coefficient vector has the wrong dimension",
          <|"MasterIndex" -> masterIndex, "Dimension" -> dimension|>],
          Module]];
      base = normalize /@ raw;
      If[physicalFromRelative =!= None,
        base = Map[Together, base . physicalFromRelative]];
      base = MapThread[Together[#1/epsSymbol^#2] &,
        {base, currentPrefactors}];
      coefficients = If[case === "integrate",
        (Together[
            Gamma[vi + vj]/(Gamma[vi]*Gamma[vj])*
            physicalVar^(vi - 1)*(1 - physicalVar)^(vj - 1)*#] & /@ base),
        base];
      (* Zero entries have infinite epsilon valuation and therefore do not
         participate in a nonzero row minimum.  Retain the conventional zero
         shift only for a row proved identically zero. *)
      nonzeroCoefficients = Select[coefficients,
        !TrueQ[PossibleZeroQ[Together[#]]] &];
      provenZero = nonzeroCoefficients === {};
      shifts = ft2ExactEpsilonValuation[#, physicalVar, epsSymbol] & /@
        nonzeroCoefficients;
      If[AnyTrue[shifts, FailureQ],
        Return[First[Select[shifts, FailureQ]], Module]];
      requestIdentity = ft2CanonicalIdentity["ft2-request-",
        KeyTake[request, {"MasterIndex", "MasterVec", "Vi", "Vj",
          "Case", "NeededVec"}]];
      coefficientIdentity = ft2CanonicalIdentity["ft2-coefficients-",
        Together /@ coefficients];
      identityHash = ft2CanonicalIdentity["",
        {level, masterIndex, batchKey, payloadKey, requestIdentity,
          coefficientIdentity}];
      <|"MasterIndex" -> masterIndex, "Case" -> case,
        "Vi" -> vi, "Vj" -> vj, "NeededVec" -> needed,
        "BatchKey" -> batchKey, "BatchPayloadKey" -> payloadKey,
        "RequestIdentity" -> requestIdentity,
        "CoefficientIdentity" -> coefficientIdentity,
        "CoefficientVector" -> coefficients,
        "ProvenZero" -> provenZero,
        "MinimumEpsilonShift" -> If[provenZero, 0, Min[shifts]],
        "Identity" -> ("ft2-level-" <> ToString[level] <> "-master-" <>
          ToString[masterIndex] <> "-" <> identityHash),
        "CheckpointIdentity" -> ("ft2-level-" <> ToString[level] <>
          "-observable-checkpoint-" <> identityHash)|>]], requests];
  If[AnyTrue[entries, FailureQ], First[Select[entries, FailureQ]], entries]];

(* Exact full-ladder epsilon planning in the relative matrix gauge.  If the
   incoming finite representation carries one common prefactor q, then its
   source edge is S=D+q while every prepared row shift is s=sbar-q.  Hence

                  S + s - delta = D + sbar - delta,

   so q cancels and must never be recursively charged as a new halo. *)
ft2BuildNativeEpsilonPlan[ftData_Association, epsilonOrder_Integer,
    halos_List, normalize_, suppliedBatches_:Automatic,
    matchingPrivateHalos_:Automatic,
    producerPrivateLosses_:Automatic] := Module[
  {levels = Lookup[ftData, "Levels", None], nLevels, previousRequired,
   previousPublicRequired,
   levelRecords = {}, runtimeLevels = <||>, levelData, matrix, gauge,
   batch, entries, active, entryLosses, baseIntrinsicLoss, intrinsicLoss,
   matchingSolveLoss, matchingHalos, producerLosses, producerPrivateLoss,
   userFloor, publicRequired,
   required, record, identity},
  nLevels = Lookup[ftData, "NumLevels", If[AssociationQ[levels],
    Length[Select[Keys[levels], IntegerQ[#] && # > 0 &]], None]];
  If[!AssociationQ[levels] || !IntegerQ[nLevels] || nLevels < 1 ||
      epsilonOrder < 0 || !AllTrue[halos, IntegerQ[#] && # >= 0 &],
    Return[ft2NativeFailure[
      "native epsilon preplanner received invalid levels, epsilon order, or level halos"],
      Module]];
  matchingHalos = Which[
    matchingPrivateHalos === Automatic, ConstantArray[0, nLevels],
    ListQ[matchingPrivateHalos] &&
        Length[matchingPrivateHalos] === nLevels &&
        AllTrue[matchingPrivateHalos, IntegerQ[#] && # >= 0 &],
      matchingPrivateHalos,
    AssociationQ[matchingPrivateHalos] &&
        AllTrue[Range[nLevels],
          IntegerQ[Lookup[matchingPrivateHalos, #, 0]] &&
            Lookup[matchingPrivateHalos, #, 0] >= 0 &],
      Lookup[matchingPrivateHalos, Range[nLevels], 0],
    True, Return[ft2NativeFailure[
      "native epsilon preplanner received invalid private matching halos",
      <|"MatchingPrivateHalos" -> matchingPrivateHalos,
        "NumLevels" -> nLevels|>], Module]];
  producerLosses = Which[
    producerPrivateLosses === Automatic, ConstantArray[0, nLevels],
    ListQ[producerPrivateLosses] &&
        Length[producerPrivateLosses] === nLevels &&
        AllTrue[producerPrivateLosses, IntegerQ[#] && # >= 0 &],
      producerPrivateLosses,
    AssociationQ[producerPrivateLosses] &&
        AllTrue[Range[nLevels],
          IntegerQ[Lookup[producerPrivateLosses, #, 0]] &&
            Lookup[producerPrivateLosses, #, 0] >= 0 &],
      Lookup[producerPrivateLosses, Range[nLevels], 0],
    True, Return[ft2NativeFailure[
      "native epsilon preplanner received invalid private producer losses",
      <|"ProducerPrivateLosses" -> producerPrivateLosses,
        "NumLevels" -> nLevels|>], Module]];
  previousRequired = epsilonOrder;
  previousPublicRequired = epsilonOrder;
  Do[
    If[!KeyExistsQ[levels, level],
      Return[ft2NativeFailure[
        "native epsilon preplanner is missing a positive FT level",
        <|"Level" -> level|>], Module]];
    levelData = levels[level];
    matrix = normalize[Lookup[levelData, "DiffMatrix", None]];
    gauge = ft2RelativeEpsilonGauge[matrix, Global`eps,
      Lookup[levelData, "FeynmanParameter", Missing["NoVariable"]]];
    If[FailureQ[gauge], Return[gauge, Module]];
    batch = If[suppliedBatches === Automatic,
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, level],
      If[AssociationQ[suppliedBatches] &&
          KeyExistsQ[suppliedBatches, level], suppliedBatches[level],
        $Failed]];
    If[batch === $Failed || !AssociationQ[batch],
      Return[ft2NativeFailure[
        "native epsilon preplanner could not obtain the exact level IBP batch",
        <|"Level" -> level|>], Module]];
    entries = ft2PrepareBoundaryEntries[level, batch,
      gauge["RelativePrefactors"],
      Lookup[levelData, "FeynmanParameter", Missing["NoVariable"]],
      Global`eps, normalize, gauge];
    If[FailureQ[entries], Return[entries, Module]];
    active = Select[entries,
      !TrueQ[Lookup[#, "ProvenZero", False]] &];
    entryLosses = Association@Map[Function[entry,
      entry["MasterIndex"] -> Max[0,
        If[entry["Case"] === "integrate", 1, 0] -
          entry["MinimumEpsilonShift"]]], active];
    baseIntrinsicLoss = If[entryLosses === <||>, 0,
      Max[Values[entryLosses]]];
    producerPrivateLoss = producerLosses[[level]];
    intrinsicLoss = baseIntrinsicLoss + producerPrivateLoss;
    (* Matching is an internal change of basis, not a physical epsilon-order
       consumer.  Its factorization/refinement halo is private to the match
       transaction and the residual certificate is authoritative only through
       RequiredOutputRawTop.  Charging elimination depth recursively here was
       the architectural error that expanded a public order-0 banana request
       to order 33. *)
    matchingSolveLoss = matchingHalos[[level]];
    userFloor = ft2UserRawFloor[epsilonOrder, halos, level];
    publicRequired = Max[userFloor,
      previousPublicRequired + intrinsicLoss];
    required = Max[userFloor,
      previousRequired + intrinsicLoss + matchingSolveLoss];
    record = <|
      "Schema" -> "FeynmanTrick.NativeEpsilonPlanLevel/v1",
      "Level" -> level,
      "GaugeIdentity" -> gauge["Identity"],
      "GaugeRecord" -> gauge["Record"],
      "RelativeGauge" -> gauge["RelativePrefactors"],
      "BatchKey" -> Lookup[batch, "Key", None],
      "BatchPayloadKey" -> Lookup[batch, "PayloadKey", None],
      "RequestIdentities" -> Lookup[entries, "RequestIdentity"],
      "RelativeCoefficientIdentities" ->
        Lookup[entries, "CoefficientIdentity"],
      "ProvenZero" -> Lookup[entries, "ProvenZero"],
      "RelativeMinimumEpsilonShifts" ->
        Lookup[entries, "MinimumEpsilonShift"],
      "EntryLosses" -> entryLosses,
      "BaseIntrinsicLoss" -> baseIntrinsicLoss,
      "ProducerPrivateLoss" -> producerPrivateLoss,
      "IntrinsicLoss" -> intrinsicLoss,
      "MatchingSolveLoss" -> matchingSolveLoss,
      "UserRawFloor" -> userFloor,
      "RequiredOutputPublicRawTop" -> previousPublicRequired,
      "RequiredOutputRawTop" -> previousRequired,
      "RequiredPublicRawTop" -> publicRequired,
      "RequiredRawTop" -> required|>;
    AppendTo[levelRecords, record];
    AssociateTo[runtimeLevels, level -> <|
      "Record" -> record, "Gauge" -> gauge, "Batch" -> batch,
      "RelativeEntries" -> entries,
      "RequiredOutputPublicRawTop" -> previousPublicRequired,
      "RequiredOutputRawTop" -> previousRequired,
      "RequiredPublicRawTop" -> publicRequired,
      "RequiredRawTop" -> required|>];
    previousPublicRequired = publicRequired;
    previousRequired = required,
    {level, 1, nLevels}];
  record = <|
    "Schema" -> "FeynmanTrick.NativeEpsilonPlan/v1",
    "EpsilonOrder" -> epsilonOrder,
    "LevelEpsilonHalos" -> halos,
    "MatchingPrivateHalos" -> matchingHalos,
    "ProducerPrivateLosses" -> producerLosses,
    "NumLevels" -> nLevels,
    "Levels" -> levelRecords,
    "DeepRequiredPublicRawTop" -> previousPublicRequired,
    "DeepRequiredRawTop" -> previousRequired|>;
  identity = ft2CanonicalIdentity["ft2-native-epsilon-plan-", record];
  <|"Schema" -> record["Schema"], "Identity" -> identity,
    "Record" -> record, "Levels" -> runtimeLevels,
    "NumLevels" -> nLevels,
    "DeepRequiredPublicRawTop" -> previousPublicRequired,
    "DeepRequiredRawTop" -> previousRequired|>];

(* The native epsilon planner has already normalized every exact level matrix
   and retained the result inside its gauge.  Re-evaluating the same
   d -> 4-2 eps substitution at runtime is not a harmless lookup: for the
   15-dimensional banana4 level it can repeat tens of minutes of symbolic
   rational-function canonicalization before the first stage marker.  Reuse
   the planner-owned matrix after checking its exact gauge binding and hash.
   The non-native backend keeps its historical normalization path. *)
ft2RuntimeLevelMatrix[levelData_Association, plannedLevel_,
    normalize_] := Module[
  {raw = Lookup[levelData, "DiffMatrix", None], gauge, matrix, record,
   levelRecord},
  If[!AssociationQ[plannedLevel], Return[normalize[raw], Module]];
  gauge = Lookup[plannedLevel, "Gauge", None];
  levelRecord = Lookup[plannedLevel, "Record", None];
  If[!AssociationQ[gauge] || !AssociationQ[levelRecord],
    Return[ft2NativeFailure[
      "native runtime level lost its planned epsilon gauge"], Module]];
  matrix = Lookup[gauge, "Matrix", None];
  record = Lookup[gauge, "Record", None];
  If[!MatrixQ[raw] || !MatrixQ[matrix] ||
      Dimensions[matrix] =!= Dimensions[raw] ||
      !AssociationQ[record] ||
      Lookup[levelRecord, "GaugeIdentity", None] =!=
        Lookup[gauge, "Identity", None] ||
      Lookup[gauge, "Identity", None] =!=
        ft2CanonicalIdentity["ft2-relative-epsilon-gauge-", record] ||
      Lookup[record, "NormalizedMatrixHash", None] =!=
        Hash[matrix, "SHA256"],
    Return[ft2NativeFailure[
      "native runtime level matrix does not match its planned epsilon gauge"],
      Module]];
  matrix];

ft2NativeEpsilonPlanQ[plan_] := Module[
  {record, levels, nLevels, levelRecords, deepRequired, deepPublic,
   matchingHalos, producerLosses},
  If[!AssociationQ[plan], Return[False, Module]];
  record = Lookup[plan, "Record", None];
  levels = Lookup[plan, "Levels", None];
  nLevels = Lookup[plan, "NumLevels", None];
  If[Lookup[plan, "Schema", None] =!=
        "FeynmanTrick.NativeEpsilonPlan/v1" ||
      !AssociationQ[record] || !AssociationQ[levels] ||
      !IntegerQ[nLevels] || nLevels < 1 ||
      Lookup[record, "Schema", None] =!=
        "FeynmanTrick.NativeEpsilonPlan/v1" ||
      Lookup[record, "NumLevels", None] =!= nLevels ||
      Lookup[plan, "Identity", None] =!=
        ft2CanonicalIdentity["ft2-native-epsilon-plan-", record],
    Return[False, Module]];
  levelRecords = Lookup[record, "Levels", None];
  deepRequired = Lookup[record, "DeepRequiredRawTop", None];
  deepPublic = Lookup[record, "DeepRequiredPublicRawTop", None];
  matchingHalos = Lookup[record, "MatchingPrivateHalos", None];
  producerLosses = Lookup[record, "ProducerPrivateLosses", None];
  If[!ListQ[levelRecords] || Length[levelRecords] =!= nLevels ||
      Sort[Keys[levels]] =!= Range[nLevels] ||
      !IntegerQ[deepRequired] || !IntegerQ[deepPublic] ||
      deepPublic > deepRequired ||
      !ListQ[matchingHalos] || Length[matchingHalos] =!= nLevels ||
      !AllTrue[matchingHalos, IntegerQ[#] && # >= 0 &] ||
      !ListQ[producerLosses] || Length[producerLosses] =!= nLevels ||
      !AllTrue[producerLosses, IntegerQ[#] && # >= 0 &] ||
      Lookup[plan, "DeepRequiredPublicRawTop", None] =!= deepPublic ||
      Lookup[plan, "DeepRequiredRawTop", None] =!= deepRequired ||
      Lookup[Last[levelRecords], "RequiredPublicRawTop", None] =!=
        deepPublic ||
      Lookup[Last[levelRecords], "RequiredRawTop", None] =!= deepRequired,
    Return[False, Module]];
  AllTrue[Range[nLevels], Function[level,
    With[{runtime = levels[level], saved = levelRecords[[level]]},
      AssociationQ[runtime] && AssociationQ[saved] &&
        Lookup[saved, "Schema", None] ===
          "FeynmanTrick.NativeEpsilonPlanLevel/v1" &&
        Lookup[saved, "Level", None] === level &&
        Lookup[runtime, "Record", None] === saved &&
        AssociationQ[Lookup[runtime, "Gauge", None]] &&
        AssociationQ[Lookup[runtime, "Batch", None]] &&
        ListQ[Lookup[runtime, "RelativeEntries", None]] &&
        Lookup[runtime, "RequiredOutputRawTop", None] ===
          Lookup[saved, "RequiredOutputRawTop", None] &&
        Lookup[runtime, "RequiredOutputPublicRawTop", None] ===
          Lookup[saved, "RequiredOutputPublicRawTop", None] &&
        Lookup[runtime, "RequiredPublicRawTop", None] ===
          Lookup[saved, "RequiredPublicRawTop", None] &&
        Lookup[runtime, "RequiredRawTop", None] ===
          Lookup[saved, "RequiredRawTop", None] &&
        IntegerQ[Lookup[saved, "RequiredOutputPublicRawTop", None]] &&
        IntegerQ[Lookup[saved, "RequiredOutputRawTop", None]] &&
        IntegerQ[Lookup[saved, "RequiredPublicRawTop", None]] &&
        IntegerQ[Lookup[saved, "RequiredRawTop", None]] &&
        Lookup[saved, "RequiredOutputPublicRawTop", 1] <=
          Lookup[saved, "RequiredOutputRawTop", 0] &&
        Lookup[saved, "RequiredPublicRawTop", 1] <=
          Lookup[saved, "RequiredRawTop", 0] &&
        Lookup[saved, "MatchingSolveLoss", None] ===
          matchingHalos[[level]] &&
        Lookup[saved, "ProducerPrivateLoss", None] ===
          producerLosses[[level]] &&
        IntegerQ[Lookup[saved, "BaseIntrinsicLoss", None]] &&
        Lookup[saved, "BaseIntrinsicLoss", -1] >= 0 &&
        Lookup[saved, "IntrinsicLoss", None] ===
          Lookup[saved, "BaseIntrinsicLoss", 0] +
            Lookup[saved, "ProducerPrivateLoss", 0] &&
        IntegerQ[Lookup[saved, "IntrinsicLoss", None]] &&
        Lookup[saved, "IntrinsicLoss", -1] >= 0 &&
        IntegerQ[Lookup[saved, "MatchingSolveLoss", None]] &&
        Lookup[saved, "MatchingSolveLoss", -1] >= 0]]]
  ];

(* DeepestLevelBoundary's order argument is the requested PHYSICAL Laurent
   top B.  Its returned finite rows are padded through

                         S = B + Max[p_i],

   and report that edge as WorkingEpsilonOrder.  Keep B and S distinct: the
   relative-gauge offset is charged against S, not against B a second time. *)
ft2DeepBoundaryWindow[deepBoundary_Association,
    requestedBoundaryOrder_Integer] := Module[
  {values = Lookup[deepBoundary, "BoundaryValues", None],
   prefactors = Lookup[deepBoundary, "EpsPrefactors", None], widths,
   completeMax, workingOrder =
     Lookup[deepBoundary, "WorkingEpsilonOrder", None],
   reportedRequested =
     Lookup[deepBoundary, "RequestedEpsilonOrder", None], expectedMax},
  If[requestedBoundaryOrder < 0 || !ListQ[values] || values === {} ||
      !AllTrue[values, ListQ] || !ListQ[prefactors] ||
      Length[prefactors] =!= Length[values] ||
      !AllTrue[prefactors, IntegerQ[#] && # >= 0 &],
    Return[ft2NativeFailure[
      "deepest boundary returned malformed finite rows or epsilon prefactors",
      <|"RequestedBoundaryOrder" -> requestedBoundaryOrder,
        "Prefactors" -> prefactors|>], Module]];
  widths = Length /@ values;
  If[MemberQ[widths, 0] || Length[DeleteDuplicates[widths]] =!= 1,
    Return[ft2NativeFailure[
      "deepest boundary rows do not have one common nonempty epsilon window",
      <|"Widths" -> widths|>], Module]];
  completeMax = First[widths] - 1;
  expectedMax = requestedBoundaryOrder + Max[prefactors];
  If[reportedRequested =!= requestedBoundaryOrder ||
      workingOrder =!= completeMax || completeMax =!= expectedMax,
    Return[ft2NativeFailure[
      "deepest boundary order metadata does not match its returned finite window",
      <|"RequestedBoundaryOrder" -> requestedBoundaryOrder,
        "ReportedRequestedEpsilonOrder" -> reportedRequested,
        "WorkingEpsilonOrder" -> workingOrder,
        "ReturnedCompleteMax" -> completeMax,
        "ExpectedCompleteMax" -> expectedMax|>], Module]];
  <|"BoundaryPrefactors" -> prefactors,
    "RequestedBoundaryOrder" -> requestedBoundaryOrder,
    "CompleteMax" -> completeMax,
    "WorkingEpsilonOrder" -> workingOrder|>];

ft2FinalizeNativeEpsilonPlan[plan_Association,
    deepBoundary_Association, requestedBoundaryOrder_Integer] := Module[
  {deepLevel, relative, window, deepPrefactors, gaugeOffset,
   requiredSourceCompleteMax, sourceCompleteMax, record, identity},
  If[!ft2NativeEpsilonPlanQ[plan],
    Return[ft2NativeFailure[
      "cannot finalize a malformed native epsilon plan"], Module]];
  window = ft2DeepBoundaryWindow[deepBoundary, requestedBoundaryOrder];
  If[FailureQ[window], Return[window, Module]];
  deepLevel = plan["Levels"][plan["NumLevels"]];
  relative = deepLevel["Gauge", "RelativePrefactors"];
  deepPrefactors = window["BoundaryPrefactors"];
  If[Length[deepPrefactors] =!= Length[relative] ||
      !AllTrue[deepPrefactors, IntegerQ],
    Return[ft2NativeFailure[
      "deepest boundary prefactors do not match the planned relative gauge",
      <|"DeepPrefactors" -> deepPrefactors,
        "RelativeGauge" -> relative|>], Module]];
  gaugeOffset = Max[deepPrefactors - relative];
  requiredSourceCompleteMax =
    plan["DeepRequiredRawTop"] + gaugeOffset;
  sourceCompleteMax = window["CompleteMax"];
  If[sourceCompleteMax < requiredSourceCompleteMax,
    Return[ft2NativeFailure[
      "deepest returned boundary window is below its exact planned gauge requirement",
      <|"RequestedBoundaryOrder" -> requestedBoundaryOrder,
        "SourceCompleteMax" -> sourceCompleteMax,
        "RequiredSourceCompleteMax" -> requiredSourceCompleteMax|>],
    Module]];
  record = <|
    "Schema" -> "FeynmanTrick.NativeEpsilonExecutionPlan/v2",
    "BasePlanIdentity" -> plan["Identity"],
    "BasePlanRecord" -> plan["Record"],
    "DeepBoundaryPrefactors" -> deepPrefactors,
    "DeepGaugeOffset" -> gaugeOffset,
    "DeepRequiredSourceCompleteMax" -> requiredSourceCompleteMax,
    "DeepRequestedBoundaryOrder" -> requestedBoundaryOrder,
    "DeepBoundaryCompleteMax" -> sourceCompleteMax,
    "DeepBoundaryWorkingEpsilonOrder" ->
      window["WorkingEpsilonOrder"],
    "DeepRequestedBoundarySurplus" ->
      requestedBoundaryOrder - plan["DeepRequiredRawTop"],
    "DeepSourceSurplus" ->
      sourceCompleteMax - requiredSourceCompleteMax|>;
  identity = ft2CanonicalIdentity[
    "ft2-native-epsilon-execution-plan-", record];
  <|"Record" -> record, "Identity" -> identity,
    "DeepGaugeOffset" -> gaugeOffset,
    "DeepRequiredSourceCompleteMax" -> requiredSourceCompleteMax,
    "DeepRequestedBoundaryOrder" -> requestedBoundaryOrder,
    "DeepBoundaryCompleteMax" -> sourceCompleteMax|>];

ft2NativeEpsilonExecutionRecordQ[record_, identity_, plan_] := Module[
  {relative, deepPrefactors, gaugeOffset, requiredSourceCompleteMax,
   requestedBoundaryOrder, sourceCompleteMax, workingOrder},
  If[!ft2NativeEpsilonPlanQ[plan] || !AssociationQ[record] ||
      Lookup[record, "Schema", None] =!=
        "FeynmanTrick.NativeEpsilonExecutionPlan/v2" ||
      identity =!= ft2CanonicalIdentity[
        "ft2-native-epsilon-execution-plan-", record] ||
      Lookup[record, "BasePlanIdentity", None] =!= plan["Identity"] ||
      Lookup[record, "BasePlanRecord", None] =!= plan["Record"],
    Return[False, Module]];
  relative = plan["Levels"][plan["NumLevels"]]["Gauge",
    "RelativePrefactors"];
  deepPrefactors = Lookup[record, "DeepBoundaryPrefactors", None];
  If[!ListQ[deepPrefactors] || Length[deepPrefactors] =!= Length[relative] ||
      !AllTrue[deepPrefactors, IntegerQ], Return[False, Module]];
  gaugeOffset = Max[deepPrefactors - relative];
  requiredSourceCompleteMax =
    plan["DeepRequiredRawTop"] + gaugeOffset;
  requestedBoundaryOrder =
    Lookup[record, "DeepRequestedBoundaryOrder", None];
  sourceCompleteMax = Lookup[record, "DeepBoundaryCompleteMax", None];
  workingOrder =
    Lookup[record, "DeepBoundaryWorkingEpsilonOrder", None];
  TrueQ[Lookup[record, "DeepGaugeOffset", None] === gaugeOffset &&
    Lookup[record, "DeepRequiredSourceCompleteMax", None] ===
      requiredSourceCompleteMax &&
    IntegerQ[requestedBoundaryOrder] && requestedBoundaryOrder >= 0 &&
    IntegerQ[sourceCompleteMax] &&
    sourceCompleteMax === workingOrder &&
    sourceCompleteMax === requestedBoundaryOrder + Max[deepPrefactors] &&
    sourceCompleteMax >= requiredSourceCompleteMax &&
    Lookup[record, "DeepRequestedBoundarySurplus", None] ===
      requestedBoundaryOrder - plan["DeepRequiredRawTop"] &&
    Lookup[record, "DeepSourceSurplus", None] ===
      sourceCompleteMax - requiredSourceCompleteMax]
  ];

ft2ValidateNativePlanRuntimeLevel[planned_Association,
    currentPrefactors_List, entries_List, epsSymbol_Symbol] := Module[
  {relative = planned["Gauge", "RelativePrefactors"], offsets,
   commonOffset, relativeEntries = planned["RelativeEntries"],
   coefficientParity, expectedShifts},
  If[Length[currentPrefactors] =!= Length[relative] ||
      !AllTrue[currentPrefactors, IntegerQ] ||
      Length[entries] =!= Length[relativeEntries],
    Return[ft2NativeFailure[
      "runtime epsilon basis does not match its planned level"], Module]];
  offsets = currentPrefactors - relative;
  If[!SameQ @@ offsets,
    Return[ft2NativeFailure[
      "runtime epsilon basis differs from the relative gauge by a noncommon shift",
      <|"RuntimePrefactors" -> currentPrefactors,
        "RelativeGauge" -> relative|>], Module]];
  commonOffset = First[offsets];
  coefficientParity = And @@ MapThread[Function[{runtime, plannedEntry},
    And @@ MapThread[TrueQ[PossibleZeroQ[Together[
        #1*epsSymbol^commonOffset - #2]]] &,
      {runtime["CoefficientVector"],
       plannedEntry["CoefficientVector"]}]],
    {entries, relativeEntries}];
  expectedShifts = MapThread[If[TrueQ[#2], 0, #1 - commonOffset] &,
    {Lookup[relativeEntries, "MinimumEpsilonShift"],
     Lookup[relativeEntries, "ProvenZero"]}];
  If[Lookup[entries, "BatchKey"] =!= Lookup[relativeEntries, "BatchKey"] ||
      Lookup[entries, "BatchPayloadKey"] =!=
        Lookup[relativeEntries, "BatchPayloadKey"] ||
      Lookup[entries, "RequestIdentity"] =!=
        Lookup[relativeEntries, "RequestIdentity"] ||
      Lookup[entries, "ProvenZero"] =!=
        Lookup[relativeEntries, "ProvenZero"] ||
      Lookup[entries, "MinimumEpsilonShift"] =!=
        expectedShifts ||
      !TrueQ[coefficientParity],
    Return[ft2NativeFailure[
      "runtime FIRE rows do not reproduce the planned common-shift invariant",
      <|"CommonOffset" -> commonOffset,
        "BatchKeyParity" ->
          (Lookup[entries, "BatchKey"] ===
            Lookup[relativeEntries, "BatchKey"]),
        "BatchPayloadParity" ->
          (Lookup[entries, "BatchPayloadKey"] ===
            Lookup[relativeEntries, "BatchPayloadKey"]),
        "RequestParity" ->
          (Lookup[entries, "RequestIdentity"] ===
            Lookup[relativeEntries, "RequestIdentity"]),
        "ZeroParity" ->
          (Lookup[entries, "ProvenZero"] ===
            Lookup[relativeEntries, "ProvenZero"]),
        "ShiftParity" ->
          (Lookup[entries, "MinimumEpsilonShift"] ===
            expectedShifts),
        "CoefficientParity" -> coefficientParity|>], Module]];
  <|"CommonOffset" -> commonOffset,
    "PlannedIntrinsicLoss" -> planned["Record", "IntrinsicLoss"]|>];

ft2NativeEpsilonLedger[entries_List, currentBCs_List,
    downstreamFiniteTop_Integer,
    downstreamPublicFiniteTop_Integer] := Module[
  {widths, availableSourceMax, active, nonDirect,
   coefficientShift, coefficientHalo, integrationHalo, entrySourceLosses,
   maximumSourceLoss, targetMax, requiredSolveMax,
   publicTargetMax, deliverableMax, publicDeliverableMax,
   maximumDeliverableMax, outputMins, capacityByMaster,
   activeCapacities},
  widths = If[AllTrue[currentBCs, ListQ], Length /@ currentBCs, {}];
  If[widths === {} || MemberQ[widths, 0] ||
      Length[DeleteDuplicates[widths]] =!= 1,
    Return[ft2NativeFailure[
      "native epsilon ledger requires one nonempty common source window",
      <|"Widths" -> widths|>], Module]];
  availableSourceMax = First[widths] - 1;
  active = Select[entries, !TrueQ[Lookup[#, "ProvenZero", False]] &];
  nonDirect = Select[active, #["Case"] =!= "direct" &];
  coefficientShift = If[nonDirect === {}, 0,
    Min[Lookup[nonDirect, "MinimumEpsilonShift"]]];
  integrationHalo = ft2NativeIntegrationHalo[active];
  (* The padded native source already carries the one global primitive order
     whenever any integrate row is present.  Compute the remaining atlas halo
     from each row's combined loss before taking a maximum.  Taking the worst
     coefficient pole and the integration halo independently can combine two
     different rows (for example, a limit row at -6 and an integrate row at
     -5) and invent one unavailable epsilon order even though every row has
     sufficient exact capacity. *)
  entrySourceLosses = Association@Map[Function[entry,
    entry["MasterIndex"] -> Max[0,
      If[entry["Case"] === "integrate", 1, 0] -
        entry["MinimumEpsilonShift"]]], nonDirect];
  maximumSourceLoss = If[entrySourceLosses === <||>, 0,
    Max[Values[entrySourceLosses]]];
  coefficientHalo = Max[0, maximumSourceLoss - integrationHalo];
  (* downstreamFiniteTop is the independently planned public output edge.
     Keep it distinct from the source reservoir: integration needs one state
     coefficient beyond that edge, while every remaining source coefficient
     is private matching work.  Promoting the whole reservoir to targetMax
     makes MatchEpsilonPadding identically zero and moves the goalpost every
     time the planner supplies more data. *)
  deliverableMax = downstreamFiniteTop;
  publicDeliverableMax = downstreamPublicFiniteTop;
  targetMax = deliverableMax + integrationHalo;
  publicTargetMax = publicDeliverableMax + integrationHalo;
  requiredSolveMax = publicDeliverableMax + maximumSourceLoss;
  outputMins = Association@Map[Function[entry,
    entry["MasterIndex"] -> If[entry["Case"] === "integrate",
      entry["MinimumEpsilonShift"] - integrationHalo,
      entry["MinimumEpsilonShift"]]], active];
  capacityByMaster = Association@Map[Function[entry,
    entry["MasterIndex"] -> (availableSourceMax +
      entry["MinimumEpsilonShift"] -
      If[entry["Case"] === "integrate", 1, 0])], active];
  activeCapacities = Values[capacityByMaster];
  maximumDeliverableMax = If[activeCapacities === {}, availableSourceMax,
    Min[Prepend[activeCapacities, availableSourceMax]]];
  If[downstreamPublicFiniteTop < 0 ||
      downstreamPublicFiniteTop > downstreamFiniteTop ||
      targetMax < 0 || publicTargetMax < 0 ||
      maximumDeliverableMax < deliverableMax ||
      requiredSolveMax > availableSourceMax,
    Return[ft2NativeFailure[
      "source epsilon depth cannot cover the downstream raw boundary window",
      <|"AvailableSourceCompleteMax" -> availableSourceMax,
        "RequiredRawTop" -> downstreamFiniteTop,
        "CoefficientHalo" -> coefficientHalo,
        "IntegrationHalo" -> integrationHalo,
        "PublicTargetCompleteMax" -> publicTargetMax,
        "RequiredSolveCompleteMax" -> requiredSolveMax,
        "TargetCompleteMax" -> targetMax,
        "PublicDeliverableCompleteMax" -> publicDeliverableMax,
        "DeliverableCompleteMax" -> deliverableMax,
        "MaximumDeliverableCompleteMax" -> maximumDeliverableMax,
        "CapacityByMaster" -> capacityByMaster|>], Module]];
  <|"AvailableSourceCompleteMax" -> availableSourceMax,
    "SourceCompleteMax" -> availableSourceMax,
    "CoefficientMinimumShift" -> coefficientShift,
    "CoefficientHalo" -> coefficientHalo,
    "IntegrationHalo" -> integrationHalo,
    "EntrySourceLosses" -> entrySourceLosses,
    "MaximumSourceLoss" -> maximumSourceLoss,
    "PublicTargetCompleteMax" -> publicTargetMax,
    "RequiredSolveCompleteMax" -> requiredSolveMax,
    "TargetCompleteMax" -> targetMax,
    "PublicDeliverableCompleteMax" -> publicDeliverableMax,
    "DeliverableCompleteMax" -> deliverableMax,
    "MaximumDeliverableCompleteMax" -> maximumDeliverableMax,
    "OutputMinimums" -> outputMins,
    "CapacityByMaster" -> capacityByMaster,
    "DownstreamFiniteTop" -> downstreamFiniteTop,
    "DownstreamReservoirRawTop" -> downstreamFiniteTop,
    "DownstreamRawTop" -> downstreamPublicFiniteTop|>];

ft2DirectBoundaryValue[entry_Association, currentBCs_List,
    physicalVar_Symbol, anchor_, epsSymbol_Symbol,
  completeMax_Integer] := Module[
  {coefficients = entry["CoefficientVector"], out = None, coefficient,
   coefficientShift, coefficientSeries, boundarySeries, term},
  Do[
    coefficient = Together[coefficients[[j]] /. physicalVar -> anchor];
    If[!FreeQ[coefficient, physicalVar],
      Return[ft2NativeFailure[
        "direct coefficient retained the Feynman parameter at the anchor",
        <|"MasterIndex" -> entry["MasterIndex"]|>], Module]];
    If[!TrueQ[PossibleZeroQ[coefficient]],
      coefficientShift = ft2ExactEpsilonValuation[
        coefficient, physicalVar, epsSymbol];
      If[FailureQ[coefficientShift], Return[coefficientShift, Module]];
      coefficientSeries = If[coefficientShift > completeMax,
        DiffExp2`EpsSeries`ESZero[completeMax],
        catch2[DiffExp2`EpsSeries`ESFromExpression[
          coefficient, epsSymbol, completeMax]]];
      If[FailureQ[coefficientSeries], Return[coefficientSeries, Module]];
      boundarySeries = DiffExp2`EpsSeries`ESNew[0, currentBCs[[j]]];
      term = DiffExp2`EpsSeries`ESTimes[coefficientSeries, boundarySeries];
      out = If[out === None, term,
        DiffExp2`EpsSeries`ESAdd[out, term]]],
    {j, Length[coefficients]}];
  If[out === None, out = DiffExp2`EpsSeries`ESZero[completeMax]];
  If[esCMx[out] < completeMax,
    Return[ft2NativeFailure[
      "direct convolution lacks the requested complete epsilon top",
      <|"MasterIndex" -> entry["MasterIndex"],
        "AvailableCompleteMax" -> esCMx[out],
        "RequiredCompleteMax" -> completeMax|>], Module]];
  If[esCMx[out] > completeMax,
    out = DiffExp2`EpsSeries`ESTruncate[out, completeMax]];
  out];

(* Overridable seams for the focused definitions-only structural test. *)
ft2NativeSegmentLine[sys_, path_] :=
  DiffExp2`Transport`SegmentLine[
    sys, path, "ValueTailContract" -> "NativeCertified"];
ft2NativePrepare[sys_, boundary_, lower_, upper_, coefficientVectors_,
    integrandRequiredMaxima_, physicalVar_, targetMax_,
    requiredTargetMax_, threads_,
    matchingCertificationDigits_] :=
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    sys, boundary, lower, upper, "Threads" -> threads,
    "Integrands" -> {coefficientVectors, physicalVar},
    "TargetCompleteMax" -> targetMax,
    "RequiredTargetCompleteMax" -> requiredTargetMax,
    "IntegrandRequiredCompleteMaxima" -> integrandRequiredMaxima,
    "MatchingCertificationDigits" -> matchingCertificationDigits,
    "DeferReceivingBases" -> True];
SetAttributes[ft2NativeRun, HoldFirst];
ft2NativeRun[atlas_Symbol, observables_, physicalVar_,
    matchingCertificationDigits_Integer,
    publicationDigits_Integer] :=
  DiffExp2`NativeTransport`RunNativeTransportObservableBatchOwned[
    atlas, observables, physicalVar,
    "ObservableContractionChunkSize" ->
      observableContractionChunkSize,
    "ObservableContractionThreads" ->
      observableContractionThreads,
    "MatchingCertificationDigits" -> matchingCertificationDigits,
    "PublicationDigits" -> publicationDigits];
ft2NativeExport[batch_, digits_] :=
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    batch, digits];
ft2NativeSaveCheckpoint[batch_, path_, identity_] :=
  DiffExp2`NativeTransport`SaveNativeTransportObservableBatchCheckpoint[
    batch, path, identity];
ft2NativeRestoreCheckpoint[manifest_] :=
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    manifest];
ft2NativeReleaseBatch[batch_] :=
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[batch];
ft2NativeReleaseAtlas[atlas_] :=
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas];

ft2NativeStageTiming[fields___] := If[
  Environment["DE2_NATIVE_STAGE_TIMING"] === "1",
  Print["FTLADDER NATIVE STAGE ", fields, " t=", SessionTime[],
    " memory=", MemoryInUse[]]];

(* Native transport is exact-singleton strict everywhere else; only FT
   integrate observables carry this explicit, identity-bound bounded-relative
   policy.  Keep two guarded digits beyond a deliberately low matching target;
   otherwise retain the fixed structural Laurent floor. *)
ft2DivergentCancellationPolicy[] := Module[
  {laurentTol, matchTol, guardDigits = DiffExp2`Tolerances`$SafetyDigits,
   tol, decimal, provenance},
  laurentTol = Together[DiffExp2`Tolerances`Tol["LaurentLeadTol"]];
  matchTol = Together[DiffExp2`Tolerances`Tol["MatchTol"]];
  tol = Max[laurentTol, matchTol/10^guardDigits];
  If[!FreeQ[tol, _?InexactNumberQ] || !TrueQ[0 < tol < 1],
    Return[ft2NativeFailure[
      "effective divergent-cancellation tolerance must be one exact number strictly between zero and one",
      <|"LaurentLeadTol" -> laurentTol, "MatchTol" -> matchTol,
        "GuardDigits" -> guardDigits, "EffectiveTolerance" -> tol|>],
      Module]];
  decimal = ToString[ScientificForm[N[tol, 50], 50,
      NumberFormat -> (Row[{#1, "e", #3}] &)], OutputForm];
  provenance = ExportString[<|
      "schema" -> "feynman-trick-divergent-cancellation-v2",
      "producer" -> "DiffExp2`Tolerances`Tol",
      "formula" -> "Max[LaurentLeadTol,MatchTol/10^SafetyDigits]",
      "laurent_lead_tol_exact" -> ToString[laurentTol, InputForm],
      "match_tol_exact" -> ToString[matchTol, InputForm],
      "safety_digits" -> guardDigits,
      "effective_exact_value" -> ToString[tol, InputForm]|>,
    "RawJSON", "Compact" -> True];
  <|"Mode" -> "bounded-relative-acb",
    "RelativeTolerance" -> decimal,
    "Provenance" -> provenance|>];

ft2RunNativeBoundaryDispatch[sys_Association, currentBCs_List,
    entries_List, ledger_Association, physicalVar_Symbol, anchor_,
    extraSingularFactors_List, deltaPrescriptions_List, threads_Integer,
    outputDigits_Integer, matchingCertificationDigits_Integer,
    publicationDigits_Integer,
    nativePlanIdentity_:None,
    checkpointSpec_:None] :=
 Block[{$MaxExtraPrecision = Max[$MaxExtraPrecision,
     DiffExp2`Tolerances`$MaxExtraPrecisionValue, 2 outputDigits]},
 Module[
  {deliverableMax = ledger["DeliverableCompleteMax"],
   publicDeliverableMax = ledger["PublicDeliverableCompleteMax"],
   publicRequiredTop = ledger["DownstreamRawTop"],
   integrationHalo = ledger["IntegrationHalo"], directRequiredTop =
     ledger["DeliverableCompleteMax"], provenZeroEntries, activeEntries,
   directEntries, nonDirectEntries, nativeEntries, directValues = <||>,
   nativeValues = <||>, zeroValues = <||>, transportSystem, lowerPlan,
   upperPlan, paddedBoundary, observables, atlas = None, batch = None,
   exported = None, exportedResults = {}, exportedValues,
   cleanupResult = None, result,
   dispatchTag = Unique["ft2NativeDispatch"], values, releaseOKQ,
   directRules, masterIndices, batchKeys, batchPayloadKeys,
   identities, checkpointIdentities, requestIdentities,
   coefficientIdentities,
   prescriptionIdentity, extraFactorsIdentity, atlasPlanIdentity = None,
   nativeBatchPayloadIdentity, publishedBatchQ, checkpointMode = None,
   checkpointContractIdentity = None, checkpointIdentity = None,
   nativeCheckpointState = None, nativeResumeRecord = None,
   restoredNativeQ = False, resumeCore, checkpointAuditRecord = None,
   publishResult = None, makeCheckpointAuditRecord, nativeBatchMatchesQ,
   nativeResultByMaster = <||>, certifications,
   divergentCancellation, integrandRequiredMaxima},
  masterIndices = Lookup[entries, "MasterIndex", {}];
  batchKeys = DeleteDuplicates[Lookup[entries, "BatchKey", {}]];
  batchPayloadKeys =
    DeleteDuplicates[Lookup[entries, "BatchPayloadKey", {}]];
  identities = Lookup[entries, "Identity", {}];
  checkpointIdentities = Lookup[entries, "CheckpointIdentity", {}];
  requestIdentities = Lookup[entries, "RequestIdentity", {}];
  coefficientIdentities = Lookup[entries, "CoefficientIdentity", {}];
  If[masterIndices =!= Range[Length[entries]] ||
      Length[batchKeys] =!= 1 || Length[batchPayloadKeys] =!= 1 ||
      !ft2CanonicalBatchKeyQ[First[batchKeys]] ||
      !ft2CanonicalBatchKeyQ[First[batchPayloadKeys]] ||
      !DuplicateFreeQ[identities] ||
      !DuplicateFreeQ[checkpointIdentities] ||
      !AllTrue[Join[identities, checkpointIdentities, requestIdentities,
          coefficientIdentities],
        StringQ[#] && StringLength[#] > 0 &] ||
      !AllTrue[Lookup[entries, "Case", {}],
        MemberQ[{"integrate", "limitLower", "limitUpper", "direct"}, #] &] ||
      !(nativePlanIdentity === None ||
        (StringQ[nativePlanIdentity] &&
          StringLength[nativePlanIdentity] > 0)) ||
      !(checkpointSpec === None || AssociationQ[checkpointSpec]),
    Return[ft2NativeFailure[
      "native boundary dispatch received inconsistent master, batch, case, or identity metadata",
      <|"MasterIndices" -> masterIndices, "BatchKeys" -> batchKeys,
        "BatchPayloadKeys" -> batchPayloadKeys|>], Module]];
  If[AssociationQ[checkpointSpec],
    checkpointMode = Lookup[checkpointSpec, "Mode", None];
    checkpointContractIdentity = Lookup[checkpointSpec,
      "ContractIdentity", None];
    If[!MemberQ[{"Save", "Restore"}, checkpointMode] ||
        !StringQ[checkpointContractIdentity] ||
        StringLength[checkpointContractIdentity] === 0 ||
        (checkpointMode === "Save" &&
          (Sort[Keys[checkpointSpec]] =!=
              Sort[{"Mode", "Path", "ContractIdentity", "Publish"}] ||
            !StringQ[checkpointSpec["Path"]] ||
            StringLength[checkpointSpec["Path"]] === 0 ||
            Head[checkpointSpec["Publish"]] =!= Function)) ||
        (checkpointMode === "Restore" &&
          (Sort[Keys[checkpointSpec]] =!=
              Sort[{"Mode", "Record", "ContractIdentity"}] ||
            !ft2NativeTransportResumeRecordQ[
              checkpointSpec["Record"]] ||
            checkpointSpec["Record", "ContractIdentity"] =!=
              checkpointContractIdentity)),
      Return[ft2NativeFailure[
        "native checkpoint specification is malformed or contract-mismatched"],
        Module]]];
  prescriptionIdentity = ft2CanonicalIdentity[
    "ft2-delta-prescriptions-", deltaPrescriptions];
  extraFactorsIdentity = ft2CanonicalIdentity[
    "ft2-extra-singular-factors-", extraSingularFactors];
  makeCheckpointAuditRecord[] := <|
    "Schema" -> "FeynmanTrick.NativeObservableBatch/v1",
    "BatchKey" -> First[batchKeys],
    "BatchPayloadKey" -> First[batchPayloadKeys],
    "RequestIdentities" -> Lookup[entries, "RequestIdentity"],
    "CoefficientIdentities" -> Lookup[entries, "CoefficientIdentity"],
    "ObservableIdentities" -> If[nativeEntries === {}, {},
      Lookup[nativeEntries, "Identity"]],
    "ObservableCheckpointIdentities" -> If[nativeEntries === {}, {},
      Lookup[nativeEntries, "CheckpointIdentity"]],
    "DeltaPrescriptions" -> deltaPrescriptions,
    "DeltaPrescriptionIdentity" -> prescriptionIdentity,
    "ExtraSingularFactorsIdentity" -> extraFactorsIdentity,
    "AtlasPlanIdentity" -> atlasPlanIdentity,
    "NativeBatchPayloadIdentity" -> nativeBatchPayloadIdentity,
    "NativeEpsilonPlanIdentity" -> nativePlanIdentity,
    "SourceCompleteMax" -> ledger["SourceCompleteMax"],
    "TargetCompleteMax" -> ledger["TargetCompleteMax"],
    "RequiredTargetCompleteMax" ->
      ledger["PublicDeliverableCompleteMax"],
    "RequiredSolveCompleteMax" ->
      ledger["RequiredSolveCompleteMax"],
    "DeliverableCompleteMax" -> deliverableMax,
    "RequiredRawTop" -> ledger["DownstreamRawTop"],
    "CoefficientHalo" -> ledger["CoefficientHalo"],
    "IntegrationHalo" -> integrationHalo,
    "MaximumSourceLoss" -> ledger["MaximumSourceLoss"],
    "MatchEpsilonPadding" -> If[AssociationQ[atlas],
      Lookup[atlas, "MatchEpsilonPadding",
        ledger["SourceCompleteMax"] -
          ledger["RequiredSolveCompleteMax"]],
      ledger["SourceCompleteMax"] -
        ledger["RequiredSolveCompleteMax"]]|>;
  provenZeroEntries = Select[entries,
    TrueQ[Lookup[#, "ProvenZero", False]] &];
  activeEntries = Select[entries,
    !TrueQ[Lookup[#, "ProvenZero", False]] &];
  directEntries = Select[activeEntries, #["Case"] === "direct" &];
  nonDirectEntries = Map[
    Join[#, <|"DeclaredOutputMin" ->
        ledger["OutputMinimums"][#["MasterIndex"]],
      "OutputMin" -> Min[
        ledger["OutputMinimums"][#["MasterIndex"]], deliverableMax]|>] &,
    Select[activeEntries, #["Case"] =!= "direct" &]];
  (* Do not prune a nonzero integral merely because its coefficient shift is
     above the requested window: endpoint primitive poles are not certified
     by that shift.  Clamping Min to D asks native transport to prove the
     zero row if the whole observable really starts later. *)
  nativeEntries = nonDirectEntries;
  If[checkpointMode === "Restore" && nativeEntries === {},
    Return[ft2NativeFailure[
      "native transport restore was requested for a level with no retained transport observables"],
      Module]];
  directRules = Map[Function[entry, Module[{value},
      value = ft2DirectBoundaryValue[entry, currentBCs, physicalVar,
        anchor, Global`eps, directRequiredTop];
      If[FailureQ[value], value, entry["MasterIndex"] -> value]]],
    directEntries];
  If[AnyTrue[directRules, FailureQ],
    Return[First[Select[directRules, FailureQ]], Module]];
  directValues = Association[directRules];
  Do[AssociateTo[zeroValues, entry["MasterIndex"] ->
      DiffExp2`EpsSeries`ESZero[deliverableMax]],
    {entry, provenZeroEntries}];
  nativeBatchPayloadIdentity = ft2CanonicalIdentity[
    "ft2-native-observable-payload-",
    {First[batchKeys], First[batchPayloadKeys],
      matchingCertificationDigits, publicationDigits,
      Map[KeyTake[#, {"MasterIndex", "Case", "RequestIdentity",
          "CoefficientIdentity", "Identity", "CheckpointIdentity"}] &,
        entries],
      KeyTake[ledger, {"SourceCompleteMax", "CoefficientHalo",
        "IntegrationHalo", "TargetCompleteMax",
        "DeliverableCompleteMax", "DownstreamRawTop"}],
      prescriptionIdentity, extraFactorsIdentity, nativePlanIdentity}];
  If[nativeEntries =!= {},
    transportSystem = Join[sys, <|"ExtraSingularFactors" ->
      Select[extraSingularFactors, !FreeQ[#, physicalVar] &]|>];
    If[checkpointMode =!= "Restore",
      ft2NativeStageTiming["lower-plan-start"];
      lowerPlan = catch2[
        ft2NativeSegmentLine[transportSystem, {anchor, 0}]];
      ft2NativeStageTiming["lower-plan-done"];
      ft2NativeStageTiming["upper-plan-start"];
      upperPlan = catch2[
        ft2NativeSegmentLine[transportSystem, {anchor, 1}]];
      ft2NativeStageTiming["upper-plan-done"];
      If[FailureQ[lowerPlan] || FailureQ[upperPlan],
        Return[First[Select[{lowerPlan, upperPlan}, FailureQ]], Module]];
      ft2NativeStageTiming["plans-ready charts=",
        Length /@ {lowerPlan["Charts"], upperPlan["Charts"]}]];
      If[Environment["FT_DUMP_NATIVE_PLANS"] === "1",
        Print["FTLADDER NATIVE PLAN DUMP lower=", InputForm[
          KeyTake[#, {"Center", "Radius", "MatchRadius", "Scale",
              "LocalRadius", "Singular", "IncomingMatchPoint"}] & /@
            lowerPlan["Charts"]]];
        Print["FTLADDER NATIVE PLAN DUMP upper=", InputForm[
          KeyTake[#, {"Center", "Radius", "MatchRadius", "Scale",
              "LocalRadius", "Singular", "IncomingMatchPoint"}] & /@
            upperPlan["Charts"]]];
        Return[ft2NativeFailure[
          "native plan dump requested; stopped before atlas preparation"],
          Module]];
    paddedBoundary = If[integrationHalo === 0,
      DiffExp2`EpsSeries`ESNew[0, #] & /@ currentBCs,
      DiffExp2`EpsSeries`ESNew[-integrationHalo,
        Join[ConstantArray[0, integrationHalo], #]] & /@ currentBCs];
    divergentCancellation = ft2DivergentCancellationPolicy[];
    If[FailureQ[divergentCancellation], Return[divergentCancellation, Module]];
    observables = Map[Function[entry, Module[{observable},
      observable = <|"Operation" -> entry["Case"],
        "Identity" -> entry["Identity"],
        "CheckpointIdentity" -> entry["CheckpointIdentity"],
        "CoefficientVector" -> entry["CoefficientVector"],
        "Epsilon" -> <|"Min" -> entry["OutputMin"],
          "Max" -> deliverableMax,
          (* Max retains the private reservoir needed by later arithmetic.
             Publication and completeness certification stop at the
             independently planned downstream public edge.  Otherwise a
             guard coefficient can reject an already certified boundary and
             turn every added producer digit into a moving goalpost. *)
          "RequiredCompleteMax" -> publicDeliverableMax|>|>;
      If[entry["Case"] === "integrate",
        (* Production transport returns the honest stored Taylor truncation.
           Full-local tail models are an optional certification product and
           can multiply memory across every projected tile; they are not
           needed for the finite-order FT ladder or analytic regularization. *)
        Join[observable, <|"TailPolicy" -> "stored",
          "DivergentCancellation" -> divergentCancellation|>],
        observable]]],
      nativeEntries];
    integrandRequiredMaxima = Map[
      #["Epsilon", "RequiredCompleteMax"] +
        If[#["Operation"] === "integrate", 1, 0] &,
      observables];
    nativeBatchMatchesQ[candidate_, expectedAtlas_] :=
      AssociationQ[candidate] &&
      Lookup[candidate, "Type", None] ===
        "DiffExp2NativeTransportObservableBatch" &&
      Lookup[candidate, "Atlas", None] === expectedAtlas &&
      ListQ[Lookup[candidate, "Results", None]] &&
      Length[candidate["Results"]] === Length[observables] &&
      Lookup[candidate["Results"], "Operation"] ===
        Lookup[observables, "Operation"] &&
      Lookup[candidate["Results"], "Identity"] ===
        Lookup[observables, "Identity"] &&
      Lookup[candidate["Results"], "CheckpointIdentity"] ===
        Lookup[observables, "CheckpointIdentity"] &&
      Lookup[candidate["Results"], "Epsilon"] ===
        Lookup[observables, "Epsilon"];
    result = Catch[Internal`WithLocalSettings[
      Null,
      If[checkpointMode === "Restore",
        resumeCore = checkpointSpec["Record"];
        nativeResumeRecord = resumeCore;
        atlasPlanIdentity = resumeCore["AtlasPlanIdentity"];
        nativeBatchPayloadIdentity = ft2CanonicalIdentity[
          "ft2-native-observable-payload-",
          {nativeBatchPayloadIdentity, atlasPlanIdentity,
            matchingCertificationDigits, publicationDigits,
            Map[KeyTake[#, {"Operation", "Identity",
                "CheckpointIdentity", "Epsilon", "TailPolicy",
                "DivergentCancellation"}] &,
              observables]}];
        checkpointIdentity = ft2CanonicalIdentity[
          "ft2-native-transport-state-",
          {checkpointContractIdentity, atlasPlanIdentity,
            nativeBatchPayloadIdentity,
            Map[KeyTake[#, {"Operation", "Identity",
                "CheckpointIdentity", "Epsilon", "TailPolicy",
                "DivergentCancellation"}] &,
              observables]}];
        If[resumeCore["NativeBatchPayloadIdentity"] =!=
              nativeBatchPayloadIdentity ||
            resumeCore["CheckpointIdentity"] =!= checkpointIdentity,
          Throw[ft2NativeFailure[
            "native checkpoint identities do not match the reconstructed observable batch"],
            dispatchTag]];
        nativeCheckpointState = resumeCore["State"];
        batch = catch2[ft2NativeRestoreCheckpoint[nativeCheckpointState]];
        restoredNativeQ = True,
        ft2NativeStageTiming["atlas-prepare-start"];
        atlas = catch2[ft2NativePrepare[transportSystem, paddedBoundary,
          lowerPlan, upperPlan, Lookup[nativeEntries, "CoefficientVector"],
          integrandRequiredMaxima, physicalVar,
          ledger["TargetCompleteMax"],
          ledger["PublicDeliverableCompleteMax"], threads,
          matchingCertificationDigits]];
        If[FailureQ[atlas] || !AssociationQ[atlas] ||
            Lookup[atlas, "Type", None] =!=
              "DiffExp2NativeRegularIndependentArmAtlas" ||
            !StringQ[Lookup[atlas, "PlanCheckpointIdentity", None]] ||
            StringLength[atlas["PlanCheckpointIdentity"]] === 0,
          Throw[If[FailureQ[atlas], atlas,
            ft2NativeFailure["native atlas preparation returned a malformed result",
              <|"Result" -> atlas|>]], dispatchTag]];
        ft2NativeStageTiming["atlas-prepare-done"];
        atlasPlanIdentity = atlas["PlanCheckpointIdentity"];
        nativeBatchPayloadIdentity = ft2CanonicalIdentity[
          "ft2-native-observable-payload-",
          {nativeBatchPayloadIdentity, atlasPlanIdentity,
            matchingCertificationDigits, publicationDigits,
            Map[KeyTake[#, {"Operation", "Identity",
                "CheckpointIdentity", "Epsilon", "TailPolicy",
                "DivergentCancellation"}] &,
              observables]}];
        ft2NativeStageTiming["observable-run-start"];
        batch = catch2[ft2NativeRun[
          atlas, observables, physicalVar,
          matchingCertificationDigits, publicationDigits]];
        ft2NativeStageTiming["observable-run-done"];
        If[checkpointMode === "Save",
          If[FailureQ[batch] || !nativeBatchMatchesQ[batch, atlas],
            Throw[If[FailureQ[batch], batch,
              ft2NativeFailure[
                "native observable batch was malformed before checkpoint save",
                <|"Result" -> batch|>]], dispatchTag]];
          checkpointIdentity = ft2CanonicalIdentity[
            "ft2-native-transport-state-",
            {checkpointContractIdentity, atlasPlanIdentity,
              nativeBatchPayloadIdentity,
              Map[KeyTake[#, {"Operation", "Identity",
                  "CheckpointIdentity", "Epsilon", "TailPolicy",
                  "DivergentCancellation"}] &,
                observables]}];
          nativeCheckpointState = catch2[ft2NativeSaveCheckpoint[
            batch, checkpointSpec["Path"], checkpointIdentity]];
          If[FailureQ[nativeCheckpointState],
            Throw[nativeCheckpointState, dispatchTag]];
          resumeCore = <|
            "Schema" -> "FeynmanTrick.NativeTransportResume/v1",
            "ContractIdentity" -> checkpointContractIdentity,
            "AtlasPlanIdentity" -> atlasPlanIdentity,
            "NativeBatchPayloadIdentity" -> nativeBatchPayloadIdentity,
            "CheckpointIdentity" -> checkpointIdentity,
            "State" -> nativeCheckpointState|>;
          nativeResumeRecord = Append[resumeCore,
            "Identity" -> ft2CanonicalIdentity[
              "ft2-native-transport-resume-", resumeCore]];
          checkpointAuditRecord = makeCheckpointAuditRecord[];
          If[!ft2NativeTransportResumeRecordQ[nativeResumeRecord] ||
              !ft2NativeCheckpointRecordQ[checkpointAuditRecord],
            Throw[ft2NativeFailure[
              "native checkpoint save produced an inconsistent resume or audit manifest"],
              dispatchTag]];
          publishResult = checkpointSpec["Publish"][
            nativeResumeRecord, checkpointAuditRecord];
          If[publishResult === $Failed || FailureQ[publishResult],
            Throw[ft2NativeFailure[
              "native checkpoint state was saved but its atomic ladder sidecar could not be published",
              <|"PublishResult" -> publishResult|>], dispatchTag]]]];
      If[FailureQ[batch] || !nativeBatchMatchesQ[batch, atlas],
        Throw[If[FailureQ[batch], batch,
          ft2NativeFailure["native observable batch returned a malformed result",
            <|"Result" -> batch|>]], dispatchTag]];
      exported = catch2[ft2NativeExport[batch, outputDigits]];
      If[FailureQ[exported] || !AssociationQ[exported] ||
          Lookup[exported, "Type", None] =!=
            "DiffExp2NativeTransportObservableBatch",
        Throw[If[FailureQ[exported], exported,
          ft2NativeFailure["native observable export returned a malformed result",
            <|"Result" -> exported|>]], dispatchTag]];
      exported,
      publishedBatchQ[candidate_] := AssociationQ[candidate] &&
        Lookup[candidate, "Type", None] ===
          "DiffExp2NativeTransportObservableBatch" &&
        Lookup[candidate, "Atlas", None] === atlas;
      cleanupResult = Which[
        publishedBatchQ[batch], catch2[ft2NativeReleaseBatch[batch]],
        AssociationQ[atlas], catch2[ft2NativeReleaseAtlas[atlas]],
        True, <|"Released" -> 0, "Failures" -> {}|>]], dispatchTag];
    releaseOKQ[release_] := AssociationQ[release] &&
      Lookup[release, "Failures", {"malformed"}] === {};
    If[FailureQ[result], Return[result, Module]];
    If[!releaseOKQ[cleanupResult],
      Return[ft2NativeFailure["native observable owner cleanup failed",
        <|"ReleaseResult" -> cleanupResult|>], Module]];
    exportedResults = Lookup[exported, "ExportedResults", None];
    If[!ListQ[exportedResults] ||
        Length[exportedResults] =!= Length[nativeEntries] ||
        Lookup[exportedResults, "Identity"] =!=
          Lookup[nativeEntries, "Identity"] ||
        !AllTrue[Lookup[exportedResults, "Value", {}],
          DiffExp2`EpsSeries`ESQ],
      Return[ft2NativeFailure[
        "native observable export changed request order or value shape"],
        Module]];
    exportedValues = Lookup[exportedResults, "Value"];
    If[AnyTrue[exportedValues, esCMx[#] < publicRequiredTop &],
      Return[ft2NativeFailure[
        "native observable export is incomplete at the required downstream raw edge"],
      Module]];
    With[{commonTop = Min[deliverableMax, Min[esCMx /@ exportedValues]]},
      exportedValues = Map[Which[
          esMn[#] > commonTop, DiffExp2`EpsSeries`ESZero[commonTop],
          esCMx[#] > commonTop,
            DiffExp2`EpsSeries`ESTruncate[#, commonTop],
          True, #] &,
        exportedValues]];
    nativeValues = AssociationThread[
      Lookup[nativeEntries, "MasterIndex"],
      exportedValues]];
    nativeResultByMaster = AssociationThread[
      Lookup[nativeEntries, "MasterIndex"], exportedResults];
  values = Map[Function[masterIndex,
      Lookup[Join[directValues, zeroValues, nativeValues], masterIndex,
        Missing["MissingBoundaryValue", masterIndex]]],
    Lookup[entries, "MasterIndex"]];
  If[AnyTrue[values, MissingQ],
    Return[ft2NativeFailure[
      "native/direct result merge did not cover every lower master",
      <|"Missing" -> Cases[values, _Missing]|>], Module]];
  With[{commonTop = Min[esCMx /@ values]},
    If[commonTop < publicRequiredTop,
      Return[ft2NativeFailure[
        "native/direct result merge lost the required downstream epsilon edge"],
        Module]];
    values = Map[Which[
        esMn[#] > commonTop, DiffExp2`EpsSeries`ESZero[commonTop],
        esCMx[#] > commonTop,
          DiffExp2`EpsSeries`ESTruncate[#, commonTop],
        True, #] &, values]];
  certifications = Map[Function[entry,
      ft2CertificationForEntry[entry,
        Lookup[nativeResultByMaster, entry["MasterIndex"], None]]],
    entries];
  If[AnyTrue[certifications, FailureQ],
    Return[First[Select[certifications, FailureQ]], Module]];
  If[checkpointAuditRecord === None,
    checkpointAuditRecord = makeCheckpointAuditRecord[]];
  <|"Values" -> values, "Certifications" -> certifications,
    "NativeBatchCalls" -> If[nativeEntries === {} || restoredNativeQ, 0, 1],
    "NativeMarches" -> If[AssociationQ[exported],
      Lookup[exported, "NativeMarches", Missing["NotReported"]], 0],
    "RestoredNativeTransport" -> restoredNativeQ,
    "NativeTransportCheckpoint" -> nativeResumeRecord,
    "CompatibilityExports" -> If[AssociationQ[exported],
      Lookup[exported, "CompatibilityExports", Length[nativeEntries]], 0],
    "CheckpointRecord" -> checkpointAuditRecord|>]];

runExample[name_String, familyRequest_:None,
    matchingPrivateHalos_:Automatic,
    matchingTaylorOrders_:Automatic,
    matchingDigitsByLevel_:Automatic,
    producerPrivateLosses_:Automatic] := Module[
  {topology, sequence, prepContract, prepKey, prepFile, ftData, outputDir, nLevels,
   boundaryOrder, deepBoundary, currentBCs, currentPrefactors,
   resumeCheckpoint = None, startLevel, finalRaw = None,
   finalCertifications = None, ftEps, dimVar,
   dimExpr, normalizeFT, nativeEpsilonPlan = None,
   baseNativeEpsilonPlan = None, matchingHaloProfileContract = None,
   matchingHaloProfileFile = "", loadedMatchingPrivateHalos = <||>,
   loadedProducerPrivateLosses = <||>,
   loadedProducerMatchingDigitExtras = {}, loadedProducerMatchingDigits = <||>,
   effectiveMatchingPrivateHalos = <||>,
   effectiveProducerPrivateLosses = <||>,
   effectiveMatchingTaylorOrders = <||>,
   effectiveMatchingDigitsByLevel = <||>,
   effectiveMatchingComputationDigitsByLevel = <||>,
   baseNativeBatches = <||>,
   nativeEpsilonExecution = None, initialDeepPrefactors,
   deepRelativeGauge, deepGaugeOffset, deepBoundaryWindow,
   deepRequiredSourceCompleteMax, deepBoundaryDeficit,
   discoveredCheckpoint, customQ, family, outputName, numericalPoint,
   outputIntegrals, outputRequests, dimensionExpression, familyID,
   pipelineRequestID, outputMode = None, outputResolution = None,
   outputResolutionLine, pretrimFinalPoleAudit, fixedParameterValues,
   prepBatches, prepSaveResult},
  customQ = familyRequest =!= None;
  If[customQ &&
      !TrueQ[FeynmanTrick`PipelineRequest`PipelineRequestQ[familyRequest]],
    Print["CUSTOM FAMILY REQUEST INVALID ", familyRequest];
    Return[$Failed]];
  If[customQ,
    family = familyRequest["Family"];
    If[name =!= "family_" <> StringTake[familyRequest["RequestID"], -16],
      Print["CUSTOM FAMILY REQUEST REJECT: internal run name is not content-addressed"];
      Return[$Failed]];
    outputMode = familyRequest["OutputMode"];
    If[!MemberQ[{"Explicit", "AllPendingDiscovery"}, outputMode],
      Print["CUSTOM FAMILY REQUEST REJECT: unsupported output mode"];
      Return[$Failed]];
    If[AnyTrue[{"AnalyticPrescription", "Prescriptions"},
        KeyExistsQ[family["Definition"], #] &],
      Print["CUSTOM FAMILY REQUEST REJECT: analytic prescriptions are not yet wired"];
      Return[$Failed]];
    outputName = family["Name"];
    numericalPoint = family["NumericalPoint"];
    outputIntegrals = If[outputMode === "Explicit",
      family["OutputIntegrals"], All];
    outputRequests = If[outputMode === "Explicit",
      familyRequest["OutputRequests"], {}];
    dimensionExpression = family["Dimension"];
    familyID = family["FamilyID"];
    pipelineRequestID = familyRequest["RequestID"];
    topology = family["Topology"];
    sequence = family["CombinationSequence"],
    outputName = name;
    numericalPoint = {};
    outputIntegrals = Automatic;
    outputRequests = Automatic;
    dimensionExpression = FTExampleDimension[name];
    familyID = None;
    pipelineRequestID = None;
    topology = FTExampleTopology[name, "step"];
    If[topology === $Failed, Return[$Failed]];
    sequence = FTExampleSequence[name]
  ];
  Print["EXAMPLE ", outputName];
  FeynmanTrick`SetFTOption["DimensionExpression", dimensionExpression];
  If[customQ && outputMode === "AllPendingDiscovery",
    outputResolution = ft2ResolveAllOutputSelection[name, familyRequest];
    If[FailureQ[outputResolution],
      Print["CUSTOM FAMILY ALL DISCOVERY FAIL ", outputResolution];
      Return[$Failed]];
    outputIntegrals = outputResolution["Masters"];
    outputRequests = outputResolution["OutputRequests"];
    Print["CUSTOM FAMILY ALL RESOLVED count=", Length[outputIntegrals],
      " identity=", outputResolution["ResolutionID"]];
    outputResolutionLine = ft2OutputLine[
      "OUTPUT_RESOLUTION ", outputResolution];
    If[FailureQ[outputResolutionLine],
      Print["CUSTOM FAMILY ALL RESOLUTION OUTPUT FAIL ",
        outputResolutionLine];
      Return[$Failed]];
    Print[outputResolutionLine]];
  $ft2ActivePipelineRequestID = pipelineRequestID;
  $ft2ActiveFamilyID = familyID;
  $ft2ActiveAllSelectionRequestID = If[AssociationQ[outputResolution],
    outputResolution["SelectionRequestID"], None];
  $ft2ActiveOutputResolutionID = If[AssociationQ[outputResolution],
    outputResolution["ResolutionID"], None];
  prepContract = If[customQ,
    ftPrepCustomContractRecord[name, familyRequest, outputResolution],
    ftPrepContractRecord[name, topology, sequence]];
  If[!AssociationQ[prepContract] ||
      AnyTrue[Join[
        Lookup[prepContract, "PreparationSources", {}],
        Lookup[prepContract, "CustomFamilySources", {}]], FailureQ],
    Print["FTPREP CONTRACT BUILD FAILED ", prepContract];
    Return[$Failed]];
  prepKey = ftPrepContractKey[prepContract];
  prepFile = ftPrepFile[name, prepKey];
  ftData = If[forcePrepRebuild, $Failed,
    loadPreparedFT[prepFile, prepContract]];
  If[ftData === $Failed && !forcePrepRebuild && migrateLegacyPrep && !customQ,
    ftData = migrateLegacyPreparedFT[name, prepFile, prepContract]];
  If[ftData === $Failed,
    Print["FTPREP CACHE MISS ", prepFile];
    ftData = If[customQ,
      FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
        topology, sequence, numericalPoint,
        "OutputIntegrals" -> outputIntegrals],
      FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
        topology, sequence, {}]];
    If[customQ && AssociationQ[ftData],
      ftData = Join[ftData, <|
          "FamilyID" -> familyID,
          "PipelineRequestID" -> pipelineRequestID|>,
        If[AssociationQ[outputResolution], <|
          "AllSelectionRequestID" -> outputResolution["SelectionRequestID"],
          "OutputResolutionID" -> outputResolution["ResolutionID"],
          "ResolvedOutputRequests" -> outputResolution["OutputRequests"]
        |>, <||>]]];
    If[ftData === $Failed, Return[$Failed]];
    outputDir = FileNameJoin[{$TemporaryDirectory,
      "FT2_" <> name <> "_" <> ToString[$ProcessID]}];
    If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
    CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
    ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
    If[preparedFTDataQ[ftData],
      (* RunFullIteration prepares the differential systems, while the exact
         boundary reductions are first requested by the native epsilon
         preplanner.  Hydrate them before writing the snapshot so a learned
         one-order retry can actually reload the expensive preparation. *)
      prepBatches = hydratePreparedReductionCache[ftData];
      If[FailureQ[prepBatches],
        Print["FTPREP REDUCTION HYDRATION FAIL ", prepBatches];
        Return[$Failed, Module]];
      prepSaveResult = savePreparedFT[
        prepFile, prepContract, ftData];
      If[prepSaveResult === $Failed && TrueQ[ft2PreparationOnly],
        Return[$Failed, Module]]]];
  If[ftData === $Failed, Return[$Failed]];
  If[TrueQ[ft2PreparationOnly],
    Print["FTPREP ONLY COMPLETE ", prepFile];
    Return[True, Module]];
  nLevels = ftData["NumLevels"];
  levelDeltaPrescriptionSigns =
    levelDeltaPrescriptionSignsForCount[nLevels];
  If[levelDeltaPrescriptionSigns === $Failed,
      Print[
        "FTLADDER DELTA PRESCRIPTION CONTRACT FAIL: expected ",
        nLevels, " per-level signs, received ",
        InputForm[configuredLevelDeltaPrescriptionSigns]];
      Return[$Failed, Module]];
  fixedParameterValues = ft2FixedParameterValues[ftData];
  If[fixedParameterValues === $Failed,
    Print["FTLADDER ANCHOR CONTRACT FAIL data=",
      KeyTake[ftData, {"NumLevels", "FixedParamValue",
        "FixedParamValues"}]];
    Return[$Failed, Module]];
  effectiveMatchingTaylorOrders = Which[
    matchingTaylorOrders === Automatic || matchingTaylorOrders === None,
      <||>,
    AssociationQ[matchingTaylorOrders] &&
        AllTrue[Keys[matchingTaylorOrders],
          IntegerQ[#] && 1 <= # <= nLevels &] &&
        AllTrue[Values[matchingTaylorOrders],
          IntegerQ[#] && # >= expansionOrder &],
      matchingTaylorOrders,
    True,
      Print["FTLADDER MATCH TAYLOR RETRY CONTRACT FAIL orders=",
        matchingTaylorOrders];
      Return[$Failed, Module]];
  effectiveMatchingDigitsByLevel = Which[
    matchingDigitsByLevel === Automatic || matchingDigitsByLevel === None,
      <||>,
    AssociationQ[matchingDigitsByLevel] &&
        AllTrue[Keys[matchingDigitsByLevel],
          IntegerQ[#] && 1 <= # <= nLevels &] &&
        AllTrue[Values[matchingDigitsByLevel],
          IntegerQ[#] && # >= matchDigits && # <= wp &],
      matchingDigitsByLevel,
    True,
      Print["FTLADDER MATCH PRODUCER RETRY CONTRACT FAIL digits=",
        matchingDigitsByLevel];
      Return[$Failed, Module]];
  ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  dimExpr = FeynmanTrick`Private`DimensionExpression[];
  normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
    ftEps -> Global`eps);
  If[recurrenceBackend === "Cpp",
    baseNativeEpsilonPlan = ft2BuildNativeEpsilonPlan[
      ftData, epsOrder, levelEpsilonHalos, normalizeFT, Automatic, <||>];
    If[FailureQ[baseNativeEpsilonPlan],
      Print["FTLADDER NATIVE BASE EPSILON PLAN FAIL ",
        baseNativeEpsilonPlan];
      Return[$Failed]];
    matchingHaloProfileContract = ft2MatchingHaloProfileContract[
      name, prepKey, baseNativeEpsilonPlan, fixedParameterValues];
    If[FailureQ[matchingHaloProfileContract],
      Print["FTLADDER MATCH PROFILE CONTRACT FAIL ",
        matchingHaloProfileContract];
      Return[$Failed]];
    matchingHaloProfileFile = ft2MatchingHaloProfileFile[
      matchingHaloProfileContract];
    loadedMatchingPrivateHalos = If[matchingHaloProfileEnabled,
      ft2LoadMatchingHaloProfile[matchingHaloProfileFile,
        matchingHaloProfileContract],
      ConstantArray[0, baseNativeEpsilonPlan["NumLevels"]]];
    loadedProducerPrivateLosses = If[matchingHaloProfileEnabled,
      ft2LoadProducerLossProfile[matchingHaloProfileFile,
        matchingHaloProfileContract],
      ConstantArray[0, baseNativeEpsilonPlan["NumLevels"]]];
    loadedProducerMatchingDigitExtras = If[matchingHaloProfileEnabled,
      ft2LoadProducerDigitProfile[matchingHaloProfileFile,
        matchingHaloProfileContract],
      ConstantArray[0, baseNativeEpsilonPlan["NumLevels"]]];
    loadedProducerMatchingDigits = Select[
      AssociationThread[
        Range[baseNativeEpsilonPlan["NumLevels"]],
        matchDigits + loadedProducerMatchingDigitExtras],
      # > matchDigits &];
    effectiveMatchingDigitsByLevel = Merge[
      {effectiveMatchingDigitsByLevel,
        loadedProducerMatchingDigits}, Max];
    (* A per-level certification override may need additional arithmetic
       guard digits without asking the preceding level to publish a tighter
       boundary.  Keep that internal computation budget separate from the
       learned producer/publication ladder.  The propagated Acb boundary
       still carries its honest enclosure at the public target. *)
    effectiveMatchingComputationDigitsByLevel = Merge[
      {effectiveMatchingDigitsByLevel,
        matchingCertificationComputationDigitsByLevel}, Max];
    effectiveMatchingPrivateHalos = ft2MergeMatchingHaloBounds[
      baseNativeEpsilonPlan["NumLevels"],
      {loadedMatchingPrivateHalos, matchingPrivateHalos}];
    If[FailureQ[effectiveMatchingPrivateHalos],
      Print["FTLADDER MATCH PROFILE MERGE FAIL ",
        effectiveMatchingPrivateHalos];
      Return[$Failed]];
    effectiveProducerPrivateLosses = ft2MergeMatchingHaloBounds[
      baseNativeEpsilonPlan["NumLevels"],
      {loadedProducerPrivateLosses, producerPrivateLosses}];
    If[FailureQ[effectiveProducerPrivateLosses],
      Print["FTLADDER PRODUCER LOSS PROFILE MERGE FAIL ",
        effectiveProducerPrivateLosses];
      Return[$Failed]];
    If[matchingHaloProfileEnabled &&
        AssociationQ[matchingDigitsByLevel] &&
        matchingDigitsByLevel =!= <||>,
      With[{savedProfile = ft2SaveMatchingHaloProfile[
          matchingHaloProfileFile, matchingHaloProfileContract,
          effectiveMatchingPrivateHalos,
          effectiveProducerPrivateLosses,
          AssociationMap[
            Max[0, Lookup[effectiveMatchingDigitsByLevel, #,
                matchDigits] - matchDigits] &,
            Range[baseNativeEpsilonPlan["NumLevels"]]]]},
        If[FailureQ[savedProfile],
          Print["FTLADDER PRODUCER DIGIT PROFILE SAVE WARNING ",
            savedProfile]]]];
    baseNativeBatches = AssociationMap[
      baseNativeEpsilonPlan["Levels"][#]["Batch"] &,
      Range[baseNativeEpsilonPlan["NumLevels"]]];
    nativeEpsilonPlan = If[
      Values[effectiveMatchingPrivateHalos] ===
          ConstantArray[0, baseNativeEpsilonPlan["NumLevels"]] &&
        Values[effectiveProducerPrivateLosses] ===
          ConstantArray[0, baseNativeEpsilonPlan["NumLevels"]],
      baseNativeEpsilonPlan,
      ft2BuildNativeEpsilonPlan[
        ftData, epsOrder, levelEpsilonHalos, normalizeFT,
        baseNativeBatches, effectiveMatchingPrivateHalos,
        effectiveProducerPrivateLosses]];
    If[FailureQ[nativeEpsilonPlan],
      Print["FTLADDER NATIVE EPSILON PLAN FAIL ", nativeEpsilonPlan];
      Return[$Failed]];
    Print["FTLADDER NATIVE EPSILON PLAN identity=",
      nativeEpsilonPlan["Identity"],
      " losses=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "IntrinsicLoss"],
      " matchLosses=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "MatchingSolveLoss"],
      " producerLosses=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "ProducerPrivateLoss"],
      " required=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "RequiredRawTop"]]];
  If[resumeLadderFile =!= "",
    resumeCheckpoint = loadLadderCheckpoint[ExpandFileName[resumeLadderFile],
      name, ftData, prepKey, nativeEpsilonPlan,
      effectiveMatchingComputationDigitsByLevel,
      effectiveMatchingDigitsByLevel];
    If[resumeCheckpoint === $Failed, Return[$Failed]],
    discoveredCheckpoint = ft2DiscoveredLadderCheckpoint[
      name, ftData, prepKey, nativeEpsilonPlan,
      effectiveMatchingComputationDigitsByLevel,
      effectiveMatchingDigitsByLevel];
    If[AssociationQ[discoveredCheckpoint],
      resumeCheckpoint = discoveredCheckpoint]];
  If[AssociationQ[resumeCheckpoint],
    startLevel = resumeCheckpoint["Level"];
    currentBCs = resumeCheckpoint["BoundaryValues"];
    currentPrefactors = resumeCheckpoint["BoundaryPrefactors"];
    If[recurrenceBackend === "Cpp",
      nativeEpsilonExecution = <|
        "Record" -> resumeCheckpoint["NativeEpsilonPlan"],
        "Identity" -> resumeCheckpoint["NativeEpsilonPlanIdentity"]|>];
    Print["FTLADDER RESUME kind=", resumeCheckpoint["Kind"],
      " level=", startLevel,
      " savedEO=", resumeCheckpoint["SourceExpansionOrder"],
      " requestedLowerEO=", expansionOrder],
    startLevel = nLevels;
    boundaryOrder = If[recurrenceBackend === "Cpp",
      nativeEpsilonPlan["DeepRequiredRawTop"] + boundaryExtraOrder,
      requestedEpsilonOrder[nLevels]];
    deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
      ftData, boundaryOrder,
      "DeltaPrescriptionSign" -> deltaPrescriptionSign];
    If[!AssociationQ[deepBoundary], Return[$Failed]];
    If[recurrenceBackend === "Cpp",
      deepBoundaryWindow =
        ft2DeepBoundaryWindow[deepBoundary, boundaryOrder];
      If[FailureQ[deepBoundaryWindow],
        Print["FTLADDER NATIVE DEEP WINDOW FAIL ", deepBoundaryWindow];
        Return[$Failed]];
      initialDeepPrefactors = deepBoundaryWindow["BoundaryPrefactors"];
      deepRelativeGauge = nativeEpsilonPlan["Levels"][nLevels]
        ["Gauge"]["RelativePrefactors"];
      If[Length[initialDeepPrefactors] =!= Length[deepRelativeGauge] ||
          !AllTrue[initialDeepPrefactors, IntegerQ],
        Print["FTLADDER NATIVE DEEP PREFAC FAIL prefactors=",
          initialDeepPrefactors, " gauge=", deepRelativeGauge];
        Return[$Failed]];
      deepGaugeOffset = Max[initialDeepPrefactors - deepRelativeGauge];
      deepRequiredSourceCompleteMax =
        nativeEpsilonPlan["DeepRequiredRawTop"] + deepGaugeOffset;
      deepBoundaryDeficit = Max[0, deepRequiredSourceCompleteMax -
        deepBoundaryWindow["CompleteMax"]];
      If[deepBoundaryDeficit > 0,
        Print["FTLADDER NATIVE DEEP RETRY requested=", boundaryOrder,
          " sourceAvailable=", deepBoundaryWindow["CompleteMax"],
          " sourceRequired=", deepRequiredSourceCompleteMax,
          " deficit=", deepBoundaryDeficit,
          " gaugeOffset=", deepGaugeOffset];
        boundaryOrder += deepBoundaryDeficit;
        deepBoundary =
          FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
            ftData, boundaryOrder,
            "DeltaPrescriptionSign" -> deltaPrescriptionSign];
        If[!AssociationQ[deepBoundary],
          Print["FTLADDER NATIVE DEEP RETRY FAIL"];
          Return[$Failed]];
        deepBoundaryWindow =
          ft2DeepBoundaryWindow[deepBoundary, boundaryOrder];
        If[FailureQ[deepBoundaryWindow] ||
            deepBoundaryWindow["BoundaryPrefactors"] =!=
              initialDeepPrefactors,
          Print["FTLADDER NATIVE DEEP RETRY FAIL window=",
            deepBoundaryWindow];
          Return[$Failed]]];
      nativeEpsilonExecution = ft2FinalizeNativeEpsilonPlan[
        nativeEpsilonPlan, deepBoundary, boundaryOrder];
      If[FailureQ[nativeEpsilonExecution],
        Print["FTLADDER NATIVE EPSILON FINALIZE FAIL ",
          nativeEpsilonExecution];
        Return[$Failed]]];
    (* coefficients are the only numerics (tags stay exact): numericize the
       deep boundary at 2x WP so the chain runs at arbitrary precision instead
       of exact-symbolic (Log/Gamma giants grind the Laurent-field algebra) *)
    currentBCs = N[deepBoundary["BoundaryValues"], inputPrecision];
    currentPrefactors = deepBoundary["EpsPrefactors"];
    If[FailureQ[printRows[outputName, nLevels,
      ftData["Levels"][nLevels]["Masters"],
      Table[DiffExp2`EpsSeries`ESShift[
        DiffExp2`EpsSeries`ESNew[0, currentBCs[[i]]], -currentPrefactors[[i]]],
        {i, Length[currentBCs]}],
      currentPrefactors,
      ConstantArray[ft2NotApplicableCertification[
        "deepestBoundary", "initial-boundary"], Length[currentBCs]]]],
      Return[$Failed]]];

  Module[{abortRes},
  abortRes = Catch[
  Do[Module[
    {levelData = ftData["Levels"][level], levelBelow = ftData["Levels"][level - 1],
     var, A, sys, mastersBelow, mastersHere, requests, reductions,
     extraFacs, rawES, trimmedRawES, trimFailureIndex,
     rawMin, shift, kmaxAvail, nextReq, needTop,
     trLoCache, trHiCache, chartCache,
     resumeTransport, resumeNativeTransport, levelExpansionOrder,
     requestedLevelExpansionOrder,
     levelEpsilonGauge,
     needInt, needLo, needHi,
     transportCheckpointFile, saveTransportProgress, completedArms,
     transportSys = None, planLo = None, planHi = None, armReq,
     loPlanCharts, hiPlanCharts, armRounds, armBatchResult,
     armUniqueCharts, armCacheCapacity, levelIBPBatch, rawExtraFacs,
     epsilonBasis, epsilonBasisRecord, nativeEntries = None,
     nativeLedger = None, nativeDispatch = None, downstreamFiniteTop,
     esCMxLevel, configResult, deltaPrescriptions, plannedLevel = None,
     runtimePlanCheck, nativeTransportContract = None,
     nativeCheckpointSpec = None, nativeStateFile, nativeSidecarFile,
     nativeConfigurationRecord, rowCertifications,
     downstreamPublicFiniteTop, levelMatchingDigits,
     downstreamPublicationDigits, levelMatchingCertificationDigits,
     anchor},
    If[recurrenceBackend === "Cpp",
      plannedLevel = nativeEpsilonPlan["Levels"][level]];
    var = levelData["FeynmanParameter"];
    anchor = fixedParameterValues[[level]];
    ft2NativeStageTiming["level=", level,
      " runtime-matrix-acquire-begin"];
    A = ft2RuntimeLevelMatrix[levelData, plannedLevel, normalizeFT];
    If[FailureQ[A],
      Print["FTLADDER RUNTIME MATRIX FAIL level=", level, " ", A];
      Throw[$Failed, "FT2Abort"]];
    ft2NativeStageTiming["level=", level,
      " runtime-matrix-acquire-ready planned=",
      AssociationQ[plannedLevel]];
    Print["LEVEL ", level, " var=", var, " d=", Length[A]];
    mastersHere = levelData["Masters"];
    mastersBelow = levelBelow["Masters"];
    resumeTransport = AssociationQ[resumeCheckpoint] &&
      resumeCheckpoint["Kind"] === "Transport" &&
      level === resumeCheckpoint["Level"];
    resumeNativeTransport = AssociationQ[resumeCheckpoint] &&
      resumeCheckpoint["Kind"] === "NativeTransport" &&
      level === resumeCheckpoint["Level"];
    ft2NativeStageTiming["level=", level,
      " epsilon-basis-gauge-extract-begin"];
    levelEpsilonGauge = If[AssociationQ[plannedLevel],
      plannedLevel["Gauge"], Automatic];
    ft2NativeStageTiming["level=", level,
      " epsilon-basis-gauge-extract-ready"];
    epsilonBasis = ft2NormalizeEpsilonBasis[
      A, currentBCs, currentPrefactors, Global`eps,
      levelEpsilonGauge, var, anchor];
    If[FailureQ[epsilonBasis],
      Print["FTLADDER EPS BASIS FAIL level=", level, " ", epsilonBasis];
      Throw[$Failed, "FT2Abort"]];
    A = epsilonBasis["Matrix"];
    currentBCs = epsilonBasis["BoundaryValues"];
    currentPrefactors = epsilonBasis["BoundaryPrefactors"];
    epsilonBasisRecord = epsilonBasis["CheckpointRecord"];
    ft2NativeStageTiming["level=", level, " epsilon-basis-ready"];
    If[resumeTransport || resumeNativeTransport,
      If[Lookup[resumeCheckpoint, "EpsilonBasis", None] =!=
          epsilonBasisRecord,
        Print["FTLADDER RESUME REJECT: epsilon basis does not match level matrix"];
        Throw[$Failed, "FT2Abort"]],
      If[Max[Flatten[epsilonBasisRecord["RawPoleOrders"]]] > 0 ||
          Max[epsilonBasis["BoundaryShifts"]]> 0,
        Print["FTLADDER EPS BASIS level=", level,
          " canonical=", epsilonBasisRecord["CanonicalPrefactors"],
          " prefactors=", currentPrefactors,
          " boundaryShifts=", epsilonBasis["BoundaryShifts"],
          " completeMax=", epsilonBasis["CompleteMax"],
          " poleFree=", epsilonBasisRecord["PoleFree"]]]];
    requestedLevelExpansionOrder = Lookup[
      effectiveMatchingTaylorOrders, level, expansionOrder];
    levelMatchingDigits = Lookup[
      effectiveMatchingComputationDigitsByLevel, level, matchDigits];
    downstreamPublicationDigits = ft2DownstreamPublicationDigits[
      level, effectiveMatchingDigitsByLevel];
    levelMatchingCertificationDigits =
      ft2LevelMatchingCertificationDigits[
        level, levelMatchingDigits, downstreamPublicationDigits];
    levelExpansionOrder = If[resumeTransport || resumeNativeTransport,
      Max[resumeCheckpoint["SourceExpansionOrder"],
        requestedLevelExpansionOrder],
      requestedLevelExpansionOrder];
    (* CoefficientVectors are part of the production transport contract.
       Rebuild/revalidate the one batched FIRE payload for every resume.  A
       schema-2 native state is restored only after these exact request and
       coefficient identities reproduce its saved contract; legacy partial
       arms remain confined to the explicit Wolfram branch. *)
    levelIBPBatch = If[AssociationQ[plannedLevel],
      plannedLevel["Batch"],
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, level]];
    If[levelIBPBatch === $Failed,
      Print["FIRE FAIL"]; Return[$Failed, Module]];
    requests = levelIBPBatch["BoundaryRequests"];
    reductions = levelIBPBatch["Reductions"];
    ft2NativeStageTiming["level=", level, " ibp-batch-ready"];
    rowCertifications = Map[Function[request,
      Switch[Lookup[request, "Case", None],
        "integrate", ft2NotApplicableCertification["integrate",
          "native-certification-unavailable"],
        "limitLower", ft2NotApplicableCertification["limitLower",
          "endpoint-limit"],
        "limitUpper", ft2NotApplicableCertification["limitUpper",
          "endpoint-limit"],
        "direct", ft2NotApplicableCertification["direct", "direct"],
        _, ft2NotApplicableCertification["unknown", "unclassified"]]],
      requests];
    rawExtraFacs =
      FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
        ftData, level, levelIBPBatch];
    If[rawExtraFacs === $Failed,
      Print["FIRE BATCH FAIL"]; Return[$Failed, Module]];
    extraFacs = normalizeFT[rawExtraFacs];
    ft2NativeStageTiming["level=", level,
      " singular-factors-ready"];
    If[resumeTransport || resumeNativeTransport,
      sys = resumeCheckpoint["System"];
      If[Lookup[sys, "Variable", None] =!= var ||
          Lookup[sys, "Matrix", None] =!= A,
        Print["FTLADDER RESUME REJECT: saved system does not match level matrix"];
        Throw[$Failed, "FT2Abort"]];
      If[Lookup[resumeCheckpoint, "Requests", None] =!= requests ||
          Lookup[resumeCheckpoint, "Reductions", None] =!=
            AssociationMap[normalizeFT, reductions] ||
          Lookup[resumeCheckpoint, "ExtraSingularFactors", None] =!=
            extraFacs,
        Print["FTLADDER RESUME REJECT: saved transport does not match the revalidated level IBP batch"];
        Throw[$Failed, "FT2Abort"]];
      If[resumeTransport, DiffExp2`Solve`ClearSolveCaches[]],
      sys = catch2[DiffExp2`API`LoadSystem[
        <|"Matrix" -> A, "Variable" -> var|>]];
      If[FailureQ[sys], Print["LOAD FAIL ", sys]; Return[$Failed, Module]]];
    If[recurrenceBackend === "Cpp",
      nativeEntries = ft2PrepareBoundaryEntries[level, levelIBPBatch,
        currentPrefactors, var, Global`eps, normalizeFT,
        levelEpsilonGauge];
      If[FailureQ[nativeEntries],
        Print["FTLADDER NATIVE COEFFICIENT FAIL level=", level, " ",
          nativeEntries];
        Throw[$Failed, "FT2Abort"]];
      runtimePlanCheck = ft2ValidateNativePlanRuntimeLevel[
        plannedLevel, currentPrefactors, nativeEntries, Global`eps];
      If[FailureQ[runtimePlanCheck],
        Print["FTLADDER NATIVE PLAN PARITY FAIL level=", level, " ",
          runtimePlanCheck];
        Throw[$Failed, "FT2Abort"]];
      downstreamFiniteTop = plannedLevel["RequiredOutputRawTop"];
      downstreamPublicFiniteTop =
        plannedLevel["RequiredOutputPublicRawTop"];
      nativeLedger = ft2NativeEpsilonLedger[
        nativeEntries, currentBCs, downstreamFiniteTop,
        downstreamPublicFiniteTop];
      If[FailureQ[nativeLedger],
        Print["FTLADDER NATIVE EPSILON FAIL level=", level, " ",
          nativeLedger];
        Throw[$Failed, "FT2Abort"]];
      ft2NativeStageTiming["level=", level,
        " epsilon-ledger=", InputForm[nativeLedger],
        " boundary-width=", Length /@ currentBCs,
        " prefactors=", currentPrefactors];
      esCMxLevel = nativeLedger["TargetCompleteMax"],
      esCMxLevel = requestedEpsilonOrder[level]];
    (* configure DiffExp2 for this level *)
    deltaPrescriptions =
      levelDeltaPrescriptions[level, var, sys, extraFacs];
    configResult = catch2[DiffExp2`Config`LoadConfiguration[{
      "WorkingPrecision" -> wp, "ExpansionOrder" -> levelExpansionOrder,
      "LinearSolveChopPrecision" -> levelMatchingDigits,
      "EpsilonOrder" -> esCMxLevel,
      "DivisionOrder" -> divisionOrder,
      "RadiusOfConvergence" -> radiusOfConvergence,
      (* The restored classic predivision planner couples placement and
         matching: adjacent regular segments meet at +1/k and -1/k. *)
      "StepDivisionOrder" -> stepDivisionOrder,
      "RecurrenceBackend" -> recurrenceBackend,
      "DeltaPrescriptions" -> deltaPrescriptions,
      "Variables" -> {}}]];
    If[FailureQ[configResult],
      Print["CONFIG FAIL level=", level, " ", configResult];
      Throw[$Failed, "FT2Abort"]];
    deltaPrescriptions = configResult["DeltaPrescriptions"];
    If[recurrenceBackend === "Cpp",
      nativeConfigurationRecord = <|
        "WorkingPrecision" -> wp,
        "MatchingDigits" -> levelMatchingDigits,
        "MatchingCertificationDigits" ->
          levelMatchingCertificationDigits,
        "PublicationDigits" -> downstreamPublicationDigits,
        "ExpansionOrder" -> levelExpansionOrder,
        "EpsilonOrder" -> esCMxLevel,
        "DivisionOrder" -> divisionOrder,
        "RadiusOfConvergence" -> radiusOfConvergence,
        "StepDivisionOrder" -> stepDivisionOrder,
        "RecurrenceBackend" -> recurrenceBackend,
        "SingularMatchPrecondition" -> singularMatchPrecondition,
        "ValueTransportMode" -> valueTransportMode,
        "NativeValueHopExecutionMode" -> nativeValueHopExecutionMode,
        "CppThreads" -> cppArmThreadBudget|>;
      nativeTransportContract = ft2NativeTransportContract[
        name, level, prepKey, sys, currentBCs, currentPrefactors,
        nativeEntries, nativeLedger, nativeConfigurationRecord,
        deltaPrescriptions, extraFacs, nativeEpsilonExecution["Identity"]];
      If[!ft2NativeTransportContractQ[nativeTransportContract],
        Print["FTLADDER NATIVE CHECKPOINT CONTRACT FAIL level=", level];
        Throw[$Failed, "FT2Abort"]];
      nativeStateFile = If[ladderCheckpointDir === "", "",
        FileNameJoin[{ladderCheckpointDir,
          name <> "_level" <> ToString[level] <>
            "_native_transport.de2cp"}]];
      nativeSidecarFile = If[ladderCheckpointDir === "", "",
        FileNameJoin[{ladderCheckpointDir,
          name <> "_level" <> ToString[level] <>
            "_native_transport.mx"}]];
      nativeCheckpointSpec = Which[
        resumeNativeTransport,
          If[nativeTransportContract =!=
              resumeCheckpoint["NativeTransportContract"],
            Print["FTLADDER NATIVE CHECKPOINT CONTRACT MISMATCH level=",
              level];
            Throw[$Failed, "FT2Abort"]];
          <|"Mode" -> "Restore",
            "Record" -> resumeCheckpoint["NativeTransportCheckpoint"],
            "ContractIdentity" -> nativeTransportContract["Identity"]|>,
        saveNativeTransportCheckpoint && nativeStateFile =!= "",
          (* A completed native-state snapshot is an optional acceleration
             point, not the default durability mechanism.  It retains every
             chart/local needed to reconstruct the live C++ session and can
             be much larger than the numerical result.  The ordinary
             boundary checkpoint written immediately after export remains
             enabled independently and is the memory-bounded default. *)
          <|"Mode" -> "Save", "Path" -> nativeStateFile,
            "ContractIdentity" -> nativeTransportContract["Identity"],
            "Publish" -> Function[{resumeRecord, auditRecord},
              saveLadderCheckpoint[nativeSidecarFile, <|
                "Kind" -> "NativeTransport", "Example" -> name,
                "Level" -> level, "PrepKey" -> prepKey,
                "System" -> sys, "Variable" -> var,
                "BoundaryValues" -> currentBCs,
                "BoundaryPrefactors" -> currentPrefactors,
                "EpsilonBasis" -> epsilonBasisRecord,
                "MastersHere" -> mastersHere,
                "MastersBelow" -> mastersBelow,
                "Requests" -> requests,
                "Reductions" -> AssociationMap[normalizeFT, reductions],
                "ExtraSingularFactors" -> extraFacs,
                "Anchor" -> anchor, "WorkingPrecision" -> wp,
                "MatchingDigits" -> levelMatchingDigits,
                "PublicationDigits" -> downstreamPublicationDigits,
                "DivisionOrder" -> divisionOrder,
                "RadiusOfConvergence" -> radiusOfConvergence,
                "ValueTransportMode" -> valueTransportMode,
                "NativeValueHopExecutionMode" ->
                  nativeValueHopExecutionMode,
                "RecurrenceBackend" -> recurrenceBackend,
                "SingularMatchPrecondition" ->
                  singularMatchPrecondition,
                "DeltaPrescriptionSign" -> deltaPrescriptionSign,
                "LevelDeltaPrescriptionSigns" ->
                  levelDeltaPrescriptionSigns,
                "EpsilonOrder" -> epsOrder,
                "BoundaryExtraOrder" -> boundaryExtraOrder,
                "LevelEpsilonHalos" -> levelEpsilonHalos,
                "ExpansionOrder" -> levelExpansionOrder,
                "RequestedEpsilonOrder" ->
                  plannedLevel["RequiredRawTop"],
                "NativeLedger" -> nativeLedger,
                "NativeObservableBatch" -> auditRecord,
                "NativeTransportContract" -> nativeTransportContract,
                "NativeTransportCheckpoint" -> resumeRecord,
                "NativeEpsilonPlan" -> nativeEpsilonExecution["Record"],
                "NativeEpsilonPlanIdentity" ->
                  nativeEpsilonExecution["Identity"],
                "Tainted" -> False|>]]|>,
        True, None];
      ft2NativeStageTiming["level=", level,
        " native-dispatch-ready"];
      nativeDispatch = ft2RunNativeBoundaryDispatch[
        sys, currentBCs, nativeEntries, nativeLedger, var, anchor,
        extraFacs, deltaPrescriptions, cppArmThreadBudget, inputPrecision,
        levelMatchingCertificationDigits,
        downstreamPublicationDigits,
        nativeEpsilonExecution["Identity"], nativeCheckpointSpec];
      If[FailureQ[nativeDispatch],
        Print["FTLADDER NATIVE BATCH FAIL level=", level, " ",
          nativeDispatch];
        With[{reservoirRetry = ft2NativeMatchingReservoirRetry[
            nativeDispatch, level]},
          With[{producerRetry = Replace[
              ft2NativeMatchingProducerRetry[
                nativeDispatch, level, nLevels,
                Lookup[effectiveMatchingDigitsByLevel, level + 1,
                  matchDigits]],
              None :> ft2NativeTerminalOutputProducerRetry[
                nativeDispatch, level, nLevels,
                levelMatchingDigits]]},
          With[{retry = Which[
              ft2NativeMatchingReservoirRetryQ[reservoirRetry],
                reservoirRetry,
              ft2NativeMatchingProducerRetryQ[producerRetry],
                producerRetry,
              True,
                ft2NativeMatchingClearanceRetry[
                  nativeDispatch, level, levelExpansionOrder]]},
          Throw[If[ft2NativeMatchingRetryQ[retry],
            Failure[retry[[1]], Join[retry[[2]], <|
              "MatchingPrivateHalos" -> effectiveMatchingPrivateHalos,
              "ProducerPrivateLosses" ->
                effectiveProducerPrivateLosses,
              "MatchingTaylorOrders" -> effectiveMatchingTaylorOrders,
              "MatchingDigitsByLevel" ->
                effectiveMatchingDigitsByLevel,
              "MatchingHaloProfileEnabled" -> matchingHaloProfileEnabled,
              "MatchingHaloProfileContract" ->
                matchingHaloProfileContract,
              "MatchingHaloProfileFile" -> matchingHaloProfileFile|>]],
            $Failed], "FT2Abort"]]]]];
      If[!resumeNativeTransport && nativeStateFile =!= "" &&
          AssociationQ[nativeDispatch["NativeTransportCheckpoint"]],
        If[!ft2NativeTransportResumeRecordQ[
            nativeDispatch["NativeTransportCheckpoint"]],
          Print["FTLADDER NATIVE CHECKPOINT MANIFEST FAIL level=", level];
          Throw[$Failed, "FT2Abort"]];
        Print["FTLADDER NATIVE CHECKPOINT READY level=", level,
          " file=", nativeSidecarFile]];
      rawES = nativeDispatch["Values"];
      rowCertifications = nativeDispatch["Certifications"];
      Print["FTLADDER NATIVE BATCH level=", level,
        " requests=", Length[nativeEntries],
        " batchCalls=", nativeDispatch["NativeBatchCalls"],
        " armMarches=", nativeDispatch["NativeMarches"],
        " restored=", nativeDispatch["RestoredNativeTransport"],
        " exports=", nativeDispatch["CompatibilityExports"],
        " HC=", nativeLedger["CoefficientHalo"],
        " HI=", nativeLedger["IntegrationHalo"],
        " sourceAvailable=", nativeLedger["AvailableSourceCompleteMax"],
        " sourceRequired=", nativeLedger["SourceCompleteMax"],
        " atlasTop=", nativeLedger["TargetCompleteMax"],
        " rawTop=", nativeLedger["DeliverableCompleteMax"],
        " t=", SessionTime[]],
    (* One pair of endpoint transports per level serves every master.  Each
       arm is checkpointed synchronously before the next arm begins. *)
    needInt = AnyTrue[requests, #["Case"] === "integrate" &];
    needLo = needInt || AnyTrue[requests, #["Case"] === "limitLower" &];
    needHi = needInt || AnyTrue[requests, #["Case"] === "limitUpper" &];
    transportCheckpointFile = FileNameJoin[{ladderCheckpointDir,
      name <> "_level" <> ToString[level] <> "_transport.mx"}];
    If[resumeTransport,
      trLoCache = resumeCheckpoint["TransportLow"];
      trHiCache = resumeCheckpoint["TransportHigh"];
      chartCache = resumeCheckpoint["ChartCache"];
      Print["FTLADDER REUSE TRANSPORT level=", level,
        " expansionOrder=", levelExpansionOrder,
        " cachedArms=", Pick[{"lower", "upper"},
          AssociationQ /@ {trLoCache, trHiCache}],
        " charts=", Length[chartCache]],
      trLoCache = None; trHiCache = None; chartCache = {}];
    saveTransportProgress[] := Module[{saved},
      If[ladderCheckpointDir === "", Return[Null, Module]];
      completedArms = Pick[{"Lower", "Upper"},
        AssociationQ /@ {trLoCache, trHiCache}];
      chartCache = Join[
        If[AssociationQ[trLoCache], trLoCache["Charts"], {}],
        If[AssociationQ[trHiCache], trHiCache["Charts"], {}]];
      saved = saveLadderCheckpoint[transportCheckpointFile, <|
        "Kind" -> "Transport", "Example" -> name, "Level" -> level,
        "PrepKey" -> prepKey, "System" -> sys,
        "Variable" -> var, "BoundaryValues" -> currentBCs,
        "BoundaryPrefactors" -> currentPrefactors,
        "EpsilonBasis" -> epsilonBasisRecord,
        "MastersHere" -> mastersHere, "MastersBelow" -> mastersBelow,
        "Requests" -> requests,
        "Reductions" -> AssociationMap[normalizeFT, reductions],
        "ExtraSingularFactors" -> extraFacs, "ChartCache" -> chartCache,
        "TransportLow" -> trLoCache, "TransportHigh" -> trHiCache,
        "CompletedArms" -> completedArms,
        "Anchor" -> anchor, "WorkingPrecision" -> wp,
        "MatchingDigits" -> levelMatchingDigits,
        "PublicationDigits" -> downstreamPublicationDigits,
        "DivisionOrder" -> divisionOrder,
        "RadiusOfConvergence" -> radiusOfConvergence,
        "ValueTransportMode" -> valueTransportMode,
        "NativeValueHopExecutionMode" -> nativeValueHopExecutionMode,
        "RecurrenceBackend" -> recurrenceBackend,
        "SingularMatchPrecondition" -> singularMatchPrecondition,
        "DeltaPrescriptionSign" -> deltaPrescriptionSign,
        "LevelDeltaPrescriptionSigns" ->
          levelDeltaPrescriptionSigns,
        "EpsilonOrder" -> epsOrder,
        "BoundaryExtraOrder" -> boundaryExtraOrder,
        "LevelEpsilonHalos" -> levelEpsilonHalos,
        "ExpansionOrder" -> levelExpansionOrder,
        "Tainted" -> If[AssociationQ[resumeCheckpoint],
          TrueQ[Lookup[resumeCheckpoint, "Tainted", False]], False],
        "RequestedEpsilonOrder" -> esCMxLevel|>];
      If[saved === $Failed, Throw[$Failed, "FT2Abort"]];
      saved];
    (* One Wolfram kernel cannot evaluate two marching loops concurrently,
       but their homogeneous chart bases do not depend on the incoming
       boundary.  For the C++ backend, collect one lower/upper chart pair at
       a time into a single native request pool.  The ordinary solve is then
       replayed with those responses, including all residual certificates,
       and its memo cache makes the subsequent marches cheap.  Pair-sized
       waves bound bridge memory and DE2_CPP_THREADS remains the sole native
       worker budget (no second Wolfram kernel/license and no oversubscription).

       This prewarm is pure cache state: it never marks an arm complete.
       The lower transport result is still synchronously checkpointed before
       the upper march starts, and a resume still computes only a missing
       arm.  Value-vector transport has boundary-dependent recurrences, so it
       intentionally keeps the established sequential path. *)
    If[cppBatchEndpointArms && recurrenceBackend === "Cpp" &&
        !valueTransport &&
        Length[A] < cppArmThreadBudget &&
        needLo && needHi && !AssociationQ[trLoCache] &&
        !AssociationQ[trHiCache],
      transportSys = Join[sys, <|"ExtraSingularFactors" ->
        Select[extraFacs, !FreeQ[#, var] &]|>];
      planLo = catch2[DiffExp2`Transport`SegmentLine[
        transportSys, {anchor, 0}]];
      planHi = catch2[DiffExp2`Transport`SegmentLine[
        transportSys, {anchor, 1}]];
      If[FailureQ[planLo] || FailureQ[planHi],
        Print["TRANSPORT PLAN FAIL ", {planLo, planHi}];
        Throw[$Failed, "FT2Abort"]];
      loPlanCharts = planLo["Charts"];
      hiPlanCharts = planHi["Charts"];
      armReq = <|"EpsWindow" -> <|"Min" -> 0,
          "CompleteMax" -> esCMxLevel|>,
        "TOrder" -> levelExpansionOrder|>;
      armRounds = Max[Length[loPlanCharts], Length[hiPlanCharts]];
      (* Preflight the complete prewarm before submitting its first wave.
         Otherwise a long pair of arms can fill the bounded homogeneous
         cache, abort a later wave, or clear all prewarmed entries when the
         ordinary march asks for its first uncached tail chart.  The chart
         count is a conservative upper bound (shared anchors are removed). *)
      armUniqueCharts = DeleteDuplicates[Join[loPlanCharts, hiPlanCharts]];
      armCacheCapacity = DiffExp2`Solve`HomogeneousCacheCapacity[];
      If[Length[armUniqueCharts] > armCacheCapacity,
        Print["FTLADDER CPP ARM BATCH SKIP level=", level,
          " uniqueCharts=", Length[armUniqueCharts],
          " cacheCapacity=", armCacheCapacity],
        Print["FTLADDER CPP ARM BATCH level=", level,
          " lowerCharts=", Length[loPlanCharts],
          " upperCharts=", Length[hiPlanCharts],
          " rounds=", armRounds];
        Do[Module[{roundCharts, roundSystems},
          roundCharts = Join[
            If[ri <= Length[loPlanCharts], {loPlanCharts[[ri]]}, {}],
            If[ri <= Length[hiPlanCharts], {hiPlanCharts[[ri]]}, {}]];
          roundSystems = catch2[
            DiffExp2`Solve`PrepareChart[transportSys, #] & /@ roundCharts];
          If[FailureQ[roundSystems],
            Print["CPP ARM PREP FAIL round=", ri, " ", roundSystems];
            Throw[$Failed, "FT2Abort"]];
          (* A single tail chart, or the identical shared anchor, has no idle
             sibling work to fill.  Let the ordinary lower-first march own it
             instead of paying the two-pass collection overhead. *)
          If[Length[roundSystems] === 2 &&
              roundSystems[[1]] =!= roundSystems[[2]],
            armBatchResult = catch2[
              DiffExp2`Solve`PrewarmHomogeneousBatch[roundSystems, armReq]];
            If[FailureQ[armBatchResult],
              Print["CPP ARM BATCH FAIL round=", ri, " ", armBatchResult];
              Throw[$Failed, "FT2Abort"]]]],
          {ri, armRounds}]]];
    If[needLo && !AssociationQ[trLoCache],
      Print["FTLADDER TRANSPORT ARM level=", level, " endpoint=lower"];
      trLoCache = catch2[If[AssociationQ[planLo],
        DiffExp2`Transport`TransportLine[
          transportSys, currentBCs, planLo],
        DiffExp2`API`TransportEndpoint[
          sys, currentBCs, anchor, 0,
          "ExtraSingularFactors" -> extraFacs]]];
      If[FailureQ[trLoCache],
        Print["TRANSPORT FAIL lower ", trLoCache];
        Throw[$Failed, "FT2Abort"]];
      (* This write must finish before the expensive upper solve starts. *)
      saveTransportProgress[]];
    If[needHi && !AssociationQ[trHiCache],
      Print["FTLADDER TRANSPORT ARM level=", level, " endpoint=upper"];
      trHiCache = catch2[If[AssociationQ[planHi],
        DiffExp2`Transport`TransportLine[
          transportSys, currentBCs, planHi],
        DiffExp2`API`TransportEndpoint[
          sys, currentBCs, anchor, 1,
          "ExtraSingularFactors" -> extraFacs]]];
      If[FailureQ[trHiCache],
        Print["TRANSPORT FAIL upper ", trHiCache];
        Throw[$Failed, "FT2Abort"]];
      saveTransportProgress[]];
    chartCache = Join[
      If[AssociationQ[trLoCache], trLoCache["Charts"], {}],
      If[AssociationQ[trHiCache], trHiCache["Charts"], {}]];
    (* Direct-only levels have no endpoint arm to trigger the write. *)
    If[!needLo && !needHi && !resumeTransport, saveTransportProgress[]];
    If[Environment["DEBUG_WINPROG"] === "1",
      Print["WINPROG level=", level,
        " reqEpsOrder=", requestedEpsilonOrder[level],
        " basisPrefactors=", Min @@ currentPrefactors, "..",
          Max @@ currentPrefactors,
        " basisCompleteMax=", epsilonBasis["CompleteMax"],
        " trLoFinalWin=", If[trLoCache === None, None,
          trLoCache["Final"]["EpsWindow"]],
        " trHiFinalWin=", If[trHiCache === None, None,
          trHiCache["Final"]["EpsWindow"]]]];
    (* per lower master: dispatch the boundary case *)
    rawES = Table[Module[
      {req = requests[[mi]], expr2, cvecBase, cvec, case, res2, vi, vj, gammaFac},
      case = req["Case"]; vi = req["Vi"]; vj = req["Vj"];
      Print["  master ", mi, " case=", case, " vi=", vi, " vj=", vj,
        " t=", SessionTime[]];
      If[!KeyExistsQ[reductions, req["NeededVec"]],
        Print["MISSING REDUCTION ", req["NeededVec"]]; Return[$Failed, Module]];
      expr2 = normalizeFT[reductions[req["NeededVec"]]];
      cvecBase = Table[
        Together[Coefficient[expr2, Global`G[1, mastersHere[[j]]]]]/
          eps^currentPrefactors[[j]],
        {j, Length[mastersHere]}];
      Switch[case,
        "integrate",
        Module[{w},
          gammaFac = Gamma[vi + vj]/(Gamma[vi]*Gamma[vj]);
          cvec = Table[Together[
            var^(vi - 1)*(1 - var)^(vj - 1)*cvecBase[[j]]], {j, Length[mastersHere]}];
          w = catch2[DiffExp2`API`LineIntegral[sys, currentBCs, anchor, {0, 1},
            cvec, "ExtraSingularFactors" -> extraFacs,
            "PrecomputedCharts" -> chartCache]];
          If[FailureQ[w], Print["INTEGRATE FAIL master ", mi, ": ", w];
            Return[$Failed, Module]];
          DiffExp2`EpsSeries`ESScale[gammaFac, w]],
        "limitUpper",
        Module[{tr = trHiCache},
          If[TrueQ[tr["EndpointIsSingular"]],
            limitCombined[tr, cvecBase, var],
            Module[{vv = tr["Value"], out = None},
              Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
                Module[{cES = DiffExp2`EpsSeries`ESFromExpression[
                    Together[cvecBase[[j]] /. var -> 1], eps, esCMx[vv]],
                    comp},
                  comp = DiffExp2`EpsSeries`ESNew[esMn[vv],
                    Table[esC[vv, k][[j]], {k, esMn[vv], esCMx[vv]}]];
                  Module[{term = DiffExp2`EpsSeries`ESTimes[cES, comp]},
                    out = If[out === None, term,
                      DiffExp2`EpsSeries`ESAdd[out, term]]]]],
                {j, Length[mastersHere]}];
              If[out === None,
                DiffExp2`EpsSeries`ESZero[esCMx[vv]], out]]]],
        "limitLower",
        Module[{tr = trLoCache},
          If[TrueQ[tr["EndpointIsSingular"]],
            limitCombined[tr, cvecBase, var],
            Module[{vv = tr["Value"], out = None},
              Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
                Module[{cES = DiffExp2`EpsSeries`ESFromExpression[
                    Together[cvecBase[[j]] /. var -> 0], eps, esCMx[vv]],
                    comp},
                  comp = DiffExp2`EpsSeries`ESNew[esMn[vv],
                    Table[esC[vv, k][[j]], {k, esMn[vv], esCMx[vv]}]];
                  Module[{term = DiffExp2`EpsSeries`ESTimes[cES, comp]},
                    out = If[out === None, term,
                      DiffExp2`EpsSeries`ESAdd[out, term]]]]],
                {j, Length[mastersHere]}];
              If[out === None,
                DiffExp2`EpsSeries`ESZero[esCMx[vv]], out]]]],
        "direct",
        Module[{out = None, kmax = requestedEpsilonOrder[level]},
          Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
            Module[{cES, bES, term},
              cES = DiffExp2`EpsSeries`ESFromExpression[
                Together[cvecBase[[j]]], eps, kmax];
              bES = DiffExp2`EpsSeries`ESNew[0, currentBCs[[j]]];
              term = DiffExp2`EpsSeries`ESTimes[cES, bES];
              out = If[out === None, term, DiffExp2`EpsSeries`ESAdd[out, term]]]],
            {j, Length[mastersHere]}];
          If[out === None, DiffExp2`EpsSeries`ESZero[kmax], out]]]],
      {mi, Length[mastersBelow]}];
    ];
    If[MemberQ[rawES, $Failed], Throw[$Failed, "FT2Abort"]];
    If[level === 1,
      pretrimFinalPoleAudit = ft2PretrimFinalPoleAudit[
        name, mastersBelow, rawES];
      If[FailureQ[pretrimFinalPoleAudit],
        Print["FTLADDER PRETRIM FINAL POLE FAIL ",
          pretrimFinalPoleAudit];
        Throw[$Failed, "FT2Abort"]];
      rawES = ft2ApplyExpectedFinalPoleFloor[name, #] & /@ rawES;
      If[AnyTrue[rawES, FailureQ],
        Print["FTLADDER FINAL POLE FLOOR FAIL ",
          First[Select[rawES, FailureQ]]];
        Throw[$Failed, "FT2Abort"]]];
    (* Matching halos extend a private high-order reservoir.  Its remote
       coefficients can be much larger than the physical pole/finite prefix
       and therefore must not participate in ESTrim's relative leading-order
       decision.  Otherwise merely increasing a halo changes already computed
       low epsilon coefficients into zeros. *)
    (* Trimming is an optimization, not permission to decide a numerically
       ambiguous Laurent coefficient.  The strict EpsSeries classifier must
       remain loud for algebraic callers; at this handoff, conservatively
       retain the untrimmed row as nonzero and pay the extra epsilon order.
       Suppress the classifier's preliminary Print because the structured
       warning below owns the row/master context. *)
    trimmedRawES = Table[Module[{candidate},
      candidate = Block[{Print = (Null &)}, catch2[
        If[recurrenceBackend === "Cpp",
          DiffExp2`EpsSeries`ESTrimThrough[
            rawES[[i]], downstreamPublicFiniteTop],
          DiffExp2`EpsSeries`ESTrim[rawES[[i]]]]]];
      If[ft2EpsilonTrimAmbiguityFailureQ[candidate],
        Print["FTLADDER EPSILON TRIM AMBIGUOUS RETAIN level=", level,
          " row=", i, " master=", mastersBelow[[i]],
          " rawWindow=", {esMn[rawES[[i]]], esCMx[rawES[[i]]]},
          " publicDecisionTop=", downstreamPublicFiniteTop,
          " coefficient=", Lookup[candidate[[2]], "Coefficient", None],
          " scale=", Lookup[candidate[[2]], "Scale", None]];
        rawES[[i]],
        candidate]],
      {i, Length[rawES]}];
    If[AnyTrue[trimmedRawES, FailureQ],
      trimFailureIndex = First[
        FirstPosition[trimmedRawES, _?FailureQ]];
      Print["FTLADDER EPSILON TRIM FAIL level=", level,
        " row=", trimFailureIndex,
        " master=", mastersBelow[[trimFailureIndex]],
        " rawWindow=", {esMn[rawES[[trimFailureIndex]]],
          esCMx[rawES[[trimFailureIndex]]]},
        " publicDecisionTop=", downstreamPublicFiniteTop,
        " failure=", trimmedRawES[[trimFailureIndex]]];
      Throw[$Failed, "FT2Abort"]];
    rawES = trimmedRawES;
    If[FailureQ[printRows[outputName, level - 1, mastersBelow, rawES,
        ConstantArray[0, Length[mastersBelow]], rowCertifications]],
      Print["FTLADDER OUTPUT CERTIFICATION FAIL level=", level - 1];
      Throw[$Failed, "FT2Abort"]];
    (* shift to finite for the next level's transport *)
    rawMin = Min[esMn /@ rawES];
    shift = Max[0, -rawMin];
    kmaxAvail = esCMx /@ rawES;
    If[level > 1,
      If[recurrenceBackend === "Cpp",
        nextReq = nativeEpsilonPlan["Levels"][level - 1]
          ["RequiredRawTop"];
        (* Carry a bounded certified reservoir to the next level.  The
           independently planned floor is mandatory; BoundaryExtraOrder is
           the explicit budget for later matching/CASE-P lattice losses.
           Propagating every deep-level halo can turn a modest public request
           into dozens of unnecessary epsilon orders and an OOM contraction. *)
        needTop = Min[Min[kmaxAvail], nextReq + boundaryExtraOrder],
        nextReq = requestedEpsilonOrder[level - 1];
        needTop = nextReq - shift];
      If[(recurrenceBackend === "Cpp" && needTop < nextReq) ||
          AnyTrue[kmaxAvail, # < needTop &],
        Print["FTLADDER INCOMPLETE level=", level - 1,
          " requiredTop=", If[recurrenceBackend === "Cpp",
            nextReq, needTop], " availableTops=", kmaxAvail,
          " shift=", shift];
        If[recurrenceBackend === "Cpp" && needTop < nextReq,
          With[{retry = ft2NativeHandoffProducerRetry[
              level - 1, level, nextReq, Min[kmaxAvail],
              nativeLedger["AvailableSourceCompleteMax"] -
                runtimePlanCheck["CommonOffset"],
              plannedLevel["Record", "BaseIntrinsicLoss"],
              plannedLevel["Record", "ProducerPrivateLoss"]]},
            Throw[If[ft2NativeProducerReservoirRetryQ[retry],
              Failure[retry[[1]], Join[retry[[2]], <|
                "MatchingPrivateHalos" ->
                  effectiveMatchingPrivateHalos,
                "ProducerPrivateLosses" ->
                  effectiveProducerPrivateLosses,
                "MatchingTaylorOrders" ->
                  effectiveMatchingTaylorOrders,
                "MatchingDigitsByLevel" ->
                  effectiveMatchingDigitsByLevel,
                "MatchingHaloProfileEnabled" ->
                  matchingHaloProfileEnabled,
                "MatchingHaloProfileContract" ->
                  matchingHaloProfileContract,
                "MatchingHaloProfileFile" ->
                  matchingHaloProfileFile|>]],
              $Failed], "FT2Abort"]],
          Throw[$Failed, "FT2Abort"]]];
      (* Numericize only at a genuine level handoff: exact Log-trees from
         tile antiderivatives otherwise compound into symbolic giants that
         grind the next level's recursion (meprec storms). *)
      currentBCs = Table[Module[{r = rawES[[i]]},
        Table[N[esC[r, k], inputPrecision], {k, -shift, needTop}]],
        {i, Length[rawES]}];
      If[!AllTrue[currentBCs,
          Length[#] === If[recurrenceBackend === "Cpp",
            needTop + shift + 1, nextReq + 1] &],
        Print["FTLADDER INTERNAL ERROR: nonuniform boundary width at level ",
          level - 1, " lengths=", Length /@ currentBCs,
          " expected=", If[recurrenceBackend === "Cpp",
            needTop + shift + 1, nextReq + 1]];
        Throw[$Failed, "FT2Abort"]];
      currentPrefactors = ConstantArray[shift, Length[rawES]];
      If[ladderCheckpointDir =!= "",
        saveLadderCheckpoint[FileNameJoin[{ladderCheckpointDir,
            name <> "_level" <> ToString[level - 1] <> "_boundary.mx"}], <|
          "Kind" -> "Boundary", "Example" -> name, "Level" -> level - 1,
          "PrepKey" -> prepKey, "BoundaryValues" -> currentBCs,
          "BoundaryPrefactors" -> currentPrefactors,
          "MastersHere" -> mastersBelow,
          "Anchor" -> If[level > 1,
            fixedParameterValues[[level - 1]], None],
          "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
          "MatchingDigits" -> levelMatchingDigits,
          "PublicationDigits" -> downstreamPublicationDigits,
          "RecurrenceBackend" -> recurrenceBackend,
          "DeltaPrescriptionSign" -> deltaPrescriptionSign,
          "LevelDeltaPrescriptionSigns" ->
            levelDeltaPrescriptionSigns,
          "BoundaryExtraOrder" -> boundaryExtraOrder,
          "LevelEpsilonHalos" -> levelEpsilonHalos,
          "SourceExpansionOrder" -> levelExpansionOrder,
          "RequestedEpsilonOrder" -> nextReq,
          "RequiredRawTop" -> If[recurrenceBackend === "Cpp",
            nextReq, None],
          "PreservedRawCompleteMax" -> If[recurrenceBackend === "Cpp",
            needTop, None],
          "BoundaryShift" -> If[recurrenceBackend === "Cpp", shift, None],
          "PreservedSourceCompleteMax" -> If[
            recurrenceBackend === "Cpp", needTop + shift, None],
          "NativeObservableBatch" -> If[recurrenceBackend === "Cpp",
            Join[nativeDispatch["CheckpointRecord"],
              <|"RequiredRawTop" -> nextReq,
                "DeliverableCompleteMax" -> needTop|>], None],
          "NativeEpsilonPlan" -> If[recurrenceBackend === "Cpp",
            nativeEpsilonExecution["Record"], None],
          "NativeEpsilonPlanIdentity" -> If[
            recurrenceBackend === "Cpp",
            nativeEpsilonExecution["Identity"], None],
          "Tainted" -> If[AssociationQ[resumeCheckpoint],
            TrueQ[Lookup[resumeCheckpoint, "Tainted", False]], False]|>]];
      If[IntegerQ[stopAfterBoundaryLevel] &&
          level - 1 === stopAfterBoundaryLevel,
        Print["STOPPED_AFTER_BOUNDARY_LEVEL ", stopAfterBoundaryLevel];
        Throw["Stopped", "FT2Abort"]],
      (* At the terminal level there is no downstream halo to populate. *)
      If[AnyTrue[kmaxAvail, # < epsOrder &],
        Print["FTLADDER INCOMPLETE FINAL requiredTop=", epsOrder,
          " availableTops=", kmaxAvail];
        Throw[$Failed, "FT2Abort"]]];
    finalRaw = rawES;
    finalCertifications = rowCertifications],
    {level, startLevel, 1, -1}], "FT2Abort"];
  If[FailureQ[abortRes], Return[abortRes]];
  If[abortRes === $Failed, Return[$Failed]];
  If[abortRes === "Stopped", Return[True]]];
  If[finalRaw === None || MemberQ[finalRaw, $Failed], Return[$Failed]];

  If[customQ,
    Module[{finalRows, finalLines},
      finalRows = ft2CustomFinalRows[
        outputName, familyID, pipelineRequestID,
        ftData["Levels"][0]["Masters"], finalRaw, finalCertifications,
        outputRequests];
      If[FailureQ[finalRows], Print["FINAL OUTPUT FAIL ", finalRows];
        Return[$Failed, Module]];
      finalLines = ft2OutputLine["FINAL ", #] & /@ finalRows;
      If[AnyTrue[finalLines, FailureQ],
        Print["FINAL OUTPUT FAIL ", First[Select[finalLines, FailureQ]]];
        Return[$Failed, Module]];
      Scan[Print, finalLines]],
    (* Preserve the registry FINAL association and serialization exactly. *)
    Module[{finalRow, finalLine},
      finalRow = ft2FinalRow[name, finalRaw[[1]],
        First[finalCertifications]];
      If[FailureQ[finalRow], Print["FINAL OUTPUT FAIL ", finalRow];
        Return[$Failed, Module]];
      finalLine = ft2OutputLine["FINAL ", finalRow];
      If[FailureQ[finalLine], Print["FINAL OUTPUT FAIL ", finalLine];
        Return[$Failed, Module]];
      Print[finalLine]]];
  True];

ft2RunExampleWithMatchingRetries[name_String,
    familyRequest_:None, initialMatchingPrivateHalos_:<||>,
    initialMatchingTaylorOrders_:<||>,
    initialMatchingDigitsByLevel_:<||>,
    initialProducerPrivateLosses_:<||>] := Module[
  {matchingPrivateHalos = initialMatchingPrivateHalos,
   matchingTaylorOrders = initialMatchingTaylorOrders,
   matchingDigitsByLevel = initialMatchingDigitsByLevel,
   producerPrivateLosses = initialProducerPrivateLosses,
   result, data, level, additional, current, merged, updated,
   currentExpansionOrder, currentTaylorOrders, profileEnabled, profileContract,
   profileFile, profileSave, nLevels, attempt = 0, maxAttempts = 12,
   residualVerdicts, currentTaylorProgress, previousTaylorProgress,
   matchingTaylorProgress = <||>,
   matchingReservoirProgress = <||>,
   producerLevel, currentMatchingDigits, currentMatchingDigitsByLevel,
   currentProducerLosses, mergedProducerLosses,
   raisedMatchingDigits, preapprovedMatchingDigits, overLimitLevels,
   backendFailure, previousBackendFailure,
   producerRetryKind, currentProducerProgress, previousProducerProgress,
   terminalOutputProducerProgress = <||>,
   maxHalo = $ft2MatchingHaloProfileMax,
   maxTaylorOrder = 4 expansionOrder,
   maxProducerDigits},
  If[!And @@ (AssociationQ /@
        {matchingPrivateHalos, matchingTaylorOrders,
          matchingDigitsByLevel, producerPrivateLosses}),
    Return[$Failed, Module]];
  While[attempt < maxAttempts,
    ++attempt;
    result = runExample[name, familyRequest, matchingPrivateHalos,
      matchingTaylorOrders, matchingDigitsByLevel,
      producerPrivateLosses];
    If[!ft2NativeMatchingRetryQ[result], Return[result, Module]];
    data = result[[2]];
    level = Lookup[data, "Level", None];
    additional = Lookup[data, "AdditionalOrders", None];
    If[!IntegerQ[level] ||
        !IntegerQ[additional] || additional < 1,
      Return[$Failed, Module]];
    If[ft2NativeProducerReservoirRetryQ[result],
      currentProducerLosses = Lookup[
        data, "ProducerPrivateLosses", producerPrivateLosses];
      profileContract = Lookup[data, "MatchingHaloProfileContract", None];
      nLevels = Lookup[profileContract, "NumLevels", None];
      If[!ft2MatchingHaloProfileContractQ[profileContract] ||
          !IntegerQ[nLevels] || level > nLevels,
        Print["FTLADDER PRODUCER LOSS RETRY PROFILE CONTRACT FAIL"];
        Return[$Failed, Module]];
      mergedProducerLosses = ft2MergeMatchingHaloBounds[
        nLevels, {producerPrivateLosses, currentProducerLosses}];
      If[FailureQ[mergedProducerLosses],
        Print["FTLADDER PRODUCER LOSS RETRY PROFILE FAIL ",
          mergedProducerLosses];
        Return[$Failed, Module]];
      updated = Lookup[mergedProducerLosses, level, 0] + additional;
      If[updated > maxHalo,
        Print["FTLADDER PRODUCER LOSS RETRY EXHAUSTED level=",
          level, " requestedPrivateLoss=", updated,
          " limit=", maxHalo];
        Return[$Failed, Module]];
      AssociateTo[mergedProducerLosses, level -> updated];
      producerPrivateLosses = mergedProducerLosses;
      current = Lookup[data, "MatchingPrivateHalos",
        matchingPrivateHalos];
      merged = ft2MergeMatchingHaloBounds[
        nLevels, {matchingPrivateHalos, current}];
      If[FailureQ[merged], Return[$Failed, Module]];
      matchingPrivateHalos = merged;
      profileEnabled = TrueQ[Lookup[data,
        "MatchingHaloProfileEnabled", False]];
      profileFile = Lookup[data, "MatchingHaloProfileFile", ""];
      If[profileEnabled && StringQ[profileFile] &&
          StringLength[profileFile] > 0,
        profileSave = ft2SaveMatchingHaloProfile[
          profileFile, profileContract, matchingPrivateHalos,
          producerPrivateLosses];
        If[FailureQ[profileSave],
          Print["FTLADDER PRODUCER LOSS PROFILE SAVE WARNING ",
            profileSave]]];
      Print["FTLADDER NATIVE PRODUCER RESERVOIR RETRY level=", level,
        " additional=", additional,
        " producerLosses=", producerPrivateLosses];
      DiffExp2`Solve`ClearSolveCaches[];
      Continue[]];
    If[ft2NativeMatchingProducerRetryQ[result],
      producerLevel = Lookup[data, "ProducerLevel", None];
      producerRetryKind = Lookup[data, "ProducerRetryKind", None];
      nLevels = Lookup[data, "NumLevels", None];
      currentMatchingDigits = Lookup[
        data, "CurrentMatchingDigits", None];
      currentMatchingDigitsByLevel = Lookup[
        data, "MatchingDigitsByLevel", matchingDigitsByLevel];
      If[!IntegerQ[producerLevel] || producerLevel <= level ||
          !IntegerQ[nLevels] || producerLevel > nLevels ||
          !IntegerQ[currentMatchingDigits] ||
          currentMatchingDigits < matchDigits ||
          !AssociationQ[currentMatchingDigitsByLevel],
        Return[$Failed, Module]];
      maxProducerDigits = ft2MatchingProducerDigitLimit[nLevels];
      If[!IntegerQ[maxProducerDigits],
        Return[$Failed, Module]];
      If[producerRetryKind === "terminal-output",
        currentProducerProgress =
          Lookup[data, "ProgressRecord", None];
        previousProducerProgress = Lookup[
          terminalOutputProducerProgress, producerLevel, None];
        If[AssociationQ[previousProducerProgress] &&
            !ft2NativeTerminalOutputProgressQ[
              previousProducerProgress, currentProducerProgress],
          Print[
            "FTLADDER NATIVE MATCH PRODUCER RETRY STALLED level=",
            level, " producerLevel=", producerLevel,
            " previousProgress=", previousProducerProgress,
            " currentProgress=", currentProducerProgress];
          Return[
            Failure[
              "FeynmanTrickNativeMatchingProducerStalled", <|
                "Detail" ->
                  "raising producer digits did not move the terminal failure later or tighten its ball",
                "Level" -> level,
                "ProducerLevel" -> producerLevel,
                "PreviousProgressRecord" ->
                  previousProducerProgress,
                "ProgressRecord" -> currentProducerProgress,
                "BackendFailure" ->
                  Lookup[data, "BackendFailure", None]|>],
            Module]];
        If[AssociationQ[currentProducerProgress],
          AssociateTo[terminalOutputProducerProgress,
            producerLevel -> currentProducerProgress]]];
      updated = currentMatchingDigits + additional;
      preapprovedMatchingDigits = Merge[
        {matchingDigitsByLevel,
          currentMatchingDigitsByLevel}, Max];
      raisedMatchingDigits =
        If[producerRetryKind === "terminal-output",
          Merge[{preapprovedMatchingDigits,
            <|producerLevel -> updated|>}, Max],
          ft2RaiseMatchingProducerDigits[
            preapprovedMatchingDigits,
            producerLevel, updated, nLevels]];
      overLimitLevels = If[AssociationQ[raisedMatchingDigits],
        Select[Keys[raisedMatchingDigits],
          raisedMatchingDigits[#] >
              Lookup[preapprovedMatchingDigits, #, matchDigits] &&
            raisedMatchingDigits[#] > maxProducerDigits &],
        {None}];
      If[!AssociationQ[raisedMatchingDigits] ||
          overLimitLevels =!= {},
        Print["FTLADDER NATIVE MATCH PRODUCER RETRY EXHAUSTED level=",
          level, " producerLevel=", producerLevel,
          " requestedMatchingDigits=", updated,
          " raisedLevelDigits=", raisedMatchingDigits,
          " newlyOverLimitLevels=", overLimitLevels,
          " limit=", maxProducerDigits];
        Return[$Failed, Module]];
      matchingDigitsByLevel = raisedMatchingDigits;
      Print["FTLADDER NATIVE MATCH PRODUCER RETRY level=", level,
        " producerLevel=", producerLevel,
        " previousMatchingDigits=", currentMatchingDigits,
        " nextMatchingDigits=", updated,
        " levelDigits=", matchingDigitsByLevel];
      DiffExp2`Solve`ClearSolveCaches[];
      Continue[]];
    If[level < 1, Return[$Failed, Module]];
    If[ft2NativeMatchingClearanceRetryQ[result],
      currentExpansionOrder = Lookup[
        data, "CurrentExpansionOrder", None];
      currentTaylorOrders = Lookup[
        data, "MatchingTaylorOrders", matchingTaylorOrders];
      If[!IntegerQ[currentExpansionOrder] || currentExpansionOrder < 1 ||
          !AssociationQ[currentTaylorOrders],
        Return[$Failed, Module]];
      residualVerdicts = Lookup[data, "ResidualVerdicts", None];
      currentTaylorProgress = Lookup[data, "ProgressRecord",
        ft2NativeMatchingClearanceProgressRecord[
          Lookup[data, "BackendFailure", None], residualVerdicts]];
      previousTaylorProgress = Lookup[
        matchingTaylorProgress, level, None];
      If[AssociationQ[previousTaylorProgress] &&
          !ft2NativeMatchingClearanceProgressQ[
            previousTaylorProgress, currentTaylorProgress],
        Print["FTLADDER NATIVE MATCH TAYLOR RETRY STALLED level=",
          level, " previousProgress=", previousTaylorProgress,
          " currentProgress=", currentTaylorProgress,
          " currentExpansionOrder=", currentExpansionOrder];
        Return[Failure["FeynmanTrickNativeMatchingTaylorStalled", <|
          "Detail" -> "raising the finite-Taylor order did not move the failure later, reduce inconclusive coefficients, or materially reduce the same-handoff midpoint defect",
          "Level" -> level,
          "CurrentExpansionOrder" -> currentExpansionOrder,
          "PreviousProgressRecord" -> previousTaylorProgress,
          "ProgressRecord" -> currentTaylorProgress,
          "PreviousResidualVerdicts" ->
            Lookup[previousTaylorProgress,
              "ResidualVerdicts", previousTaylorProgress],
          "ResidualVerdicts" -> residualVerdicts,
          "BackendFailure" -> Lookup[data, "BackendFailure", None]|>],
          Module]];
      If[AssociationQ[currentTaylorProgress],
        AssociateTo[
          matchingTaylorProgress, level -> currentTaylorProgress]];
      updated = currentExpansionOrder + additional;
      If[updated > maxTaylorOrder,
        Print["FTLADDER NATIVE MATCH TAYLOR RETRY EXHAUSTED level=",
          level, " requestedExpansionOrder=", updated,
          " limit=", maxTaylorOrder];
        Return[$Failed, Module]];
      matchingTaylorOrders = Merge[
        {matchingTaylorOrders, currentTaylorOrders}, Max];
      AssociateTo[matchingTaylorOrders, level -> updated];
      Print["FTLADDER NATIVE MATCH TAYLOR RETRY level=", level,
        " previousExpansionOrder=", currentExpansionOrder,
        " nextExpansionOrder=", updated,
        " levelOrders=", matchingTaylorOrders];
      DiffExp2`Solve`ClearSolveCaches[];
      Continue[]];
    current = Lookup[data, "MatchingPrivateHalos", matchingPrivateHalos];
    profileContract = Lookup[data, "MatchingHaloProfileContract", None];
    nLevels = Lookup[profileContract, "NumLevels", None];
    If[!ft2MatchingHaloProfileContractQ[profileContract] ||
        !IntegerQ[nLevels] || level > nLevels,
      Print["FTLADDER NATIVE MATCH RETRY PROFILE CONTRACT FAIL"];
      Return[$Failed, Module]];
    merged = ft2MergeMatchingHaloBounds[
      nLevels,
      {matchingPrivateHalos, current}];
    If[FailureQ[merged],
      Print["FTLADDER NATIVE MATCH RETRY PROFILE FAIL ", merged];
      Return[$Failed, Module]];
    backendFailure = Lookup[data, "BackendFailure", None];
    previousBackendFailure = Lookup[
      matchingReservoirProgress, level, None];
    updated = ft2NativeMatchingReservoirBackoff[
      Lookup[merged, level, 0], additional,
      previousBackendFailure, backendFailure];
    If[updated === $Failed, Return[$Failed, Module]];
    If[AssociationQ[backendFailure],
      AssociateTo[matchingReservoirProgress,
        level -> backendFailure]];
    If[updated > maxHalo,
      Print["FTLADDER NATIVE MATCH RETRY EXHAUSTED level=", level,
        " requestedPrivateHalo=", updated, " limit=", maxHalo];
      Return[$Failed, Module]];
    AssociateTo[merged, level -> updated];
    matchingPrivateHalos = merged;
    profileEnabled = TrueQ[Lookup[data,
      "MatchingHaloProfileEnabled", False]];
    profileFile = Lookup[data, "MatchingHaloProfileFile", ""];
    If[profileEnabled && ft2MatchingHaloProfileContractQ[profileContract] &&
        StringQ[profileFile] && StringLength[profileFile] > 0,
      profileSave = ft2SaveMatchingHaloProfile[
        profileFile, profileContract, matchingPrivateHalos,
        producerPrivateLosses];
      If[FailureQ[profileSave],
        Print["FTLADDER MATCH PROFILE SAVE WARNING ", profileSave]]];
    Print["FTLADDER NATIVE MATCH RETRY level=", level,
      " additional=", additional,
      " privateHalos=", matchingPrivateHalos];
    DiffExp2`Solve`ClearSolveCaches[]];
  Print["FTLADDER NATIVE MATCH RETRY EXHAUSTED attempts=", maxAttempts,
    " privateHalos=", matchingPrivateHalos,
    " producerLosses=", producerPrivateLosses,
    " taylorOrders=", matchingTaylorOrders,
    " levelDigits=", matchingDigitsByLevel];
  $Failed];

(* Let the focused checkpoint tests load these definitions without starting
   FIRE or terminating their Wolfram kernel. *)
If[envOrDefault["FT_RUNNER_DEFINITIONS_ONLY", "0"] =!= "1",
  If[AssociationQ[ft2FamilyRequest],
    If[ft2RunExampleWithMatchingRetries[
        ft2FamilyRunName, ft2FamilyRequest] === $Failed,
      Print["FAILED ", ft2FamilyRunName]; Quit[1]],
    requested = StringTrim /@
      StringSplit[envOrDefault["FT_EXAMPLES", "bubble"], ","];
    Do[
      If[ft2RunExampleWithMatchingRetries[name] === $Failed,
        Print["FAILED ", name]; Quit[1]],
      {name, requested}]];
  Quit[0]];
