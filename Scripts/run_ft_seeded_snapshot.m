(* Diagnostic runner for an explicitly selected prepared Feynman-trick
   snapshot.  This bypasses only FIRE preparation-cache source provenance;
   the snapshot must still match the requested built-in topology, combination
   sequence, dimension, and preparation configuration.  Transport and native
   checkpoints retain their ordinary validation.

   Required environment:
     FT_SEEDED_PREP_SNAPSHOT=/absolute/path/example_digest.mx
     FT_EXAMPLES=one-built-in-example

   A validated FT_FAMILY_REQUEST_FILE/FT_FAMILY_REQUEST_ID pair may select a
   custom family instead of FT_EXAMPLES.  This is useful when a source-only
   preparation-contract change invalidates the cache filename but the exact
   topology, output selection, FIRE runtime, and prepared reductions remain
   unchanged.

   This runner deliberately performs one attempt.  A structured epsilon retry
   is returned to the caller instead of silently changing the experiment. *)

seededSnapshotPath = Environment["FT_SEEDED_PREP_SNAPSHOT"];
If[!StringQ[seededSnapshotPath] ||
    StringLength[StringTrim[seededSnapshotPath]] == 0,
  Print["FT_SEEDED_PREP_SNAPSHOT is required"];
  Exit[2]];
seededSnapshotPath = ExpandFileName[seededSnapshotPath];
If[!FileExistsQ[seededSnapshotPath],
  Print["seeded preparation snapshot does not exist: ",
    seededSnapshotPath];
  Exit[2]];

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

seededCustomQ = AssociationQ[ft2FamilyRequest];
seededExample = If[seededCustomQ, ft2FamilyRunName,
  StringTrim[Environment["FT_EXAMPLES"]]];
If[!seededCustomQ &&
    (seededExample == "" || StringContainsQ[seededExample, ","]),
  Print["FT_EXAMPLES must select exactly one built-in example when no custom family request is supplied"];
  Exit[2]];
seededPrivateHaloText =
  StringTrim[envOrDefault["FT_SEEDED_MATCHING_PRIVATE_HALOS", ""]];
seededPrivateHaloParts =
  StringTrim /@ StringSplit[seededPrivateHaloText, ","];
seededPrivateHalos = Which[
  seededPrivateHaloText === "", Automatic,
  seededPrivateHaloParts =!= {} &&
      AllTrue[seededPrivateHaloParts,
        StringMatchQ[#, RegularExpression["[0-9]+"]] &],
    FromDigits /@ seededPrivateHaloParts,
  True, $Failed];
If[seededPrivateHalos === $Failed,
  Print["invalid FT_SEEDED_MATCHING_PRIVATE_HALOS: ",
    seededPrivateHaloText];
  Exit[2]];

seededPayload = readPreparedFTPayload[seededSnapshotPath];
If[!AssociationQ[seededPayload] ||
    Lookup[seededPayload, "Version", None] =!= $ftPrepCacheVersion ||
    !AssociationQ[Lookup[seededPayload, "Contract", None]] ||
    !AssociationQ[Lookup[seededPayload, "FTData", None]] ||
    !AssociationQ[Lookup[seededPayload, "ReductionCache", None]] ||
    !preparedReductionCacheQ[
      seededPayload["FTData"], seededPayload["ReductionCache"]],
  Print["seeded preparation snapshot is malformed or incomplete"];
  Exit[2]];

seededContractKeys = Join[{
  "Example", "Topology", "CombinationSequence", "DimensionExpression",
  "NumericalPoint", "PreparationConfiguration", "FIRERuntime"},
  If[seededCustomQ, {
    "CustomFamilyPreparationSchema", "FamilyID", "PipelineRequestID",
    "OutputSelectionMode", "OutputSelection", "AllSelectionRequestID",
    "OutputResolutionID", "ResolvedOutputSelection",
    "ResolvedOutputRequests"}, {}]];
seededPrepKey = Lookup[seededPayload, "Key", None];
If[!IntegerQ[seededPrepKey],
  Print["seeded preparation snapshot has no integer preparation key"];
  Exit[2]];

(* Boundary checkpoints produced from this explicit snapshot bind its
   preparation key. Reuse that exact key as part of the same stale-source
   opt-in; topology/configuration parity is still checked below. *)
ftPrepContractKey[contract_Association] := seededPrepKey;

seededLoadPreparedFT[file_String, contract_Association] := Module[
  {savedContract = seededPayload["Contract"]},
  If[Lookup[savedContract, seededContractKeys,
        Missing["Absent"]] =!=
      Lookup[contract, seededContractKeys, Missing["Absent"]],
    Print["FTPREP SEEDED SNAPSHOT CONTRACT MISMATCH ",
      seededSnapshotPath];
    Return[$Failed, Module]];
  FeynmanTrick`FIREInterface`Private`$ReductionCache =
    Join[FeynmanTrick`FIREInterface`Private`$ReductionCache,
      seededPayload["ReductionCache"]];
  Print["FTPREP SEEDED SNAPSHOT ", seededSnapshotPath,
    " requestedCache=", file];
  seededPayload["FTData"]];

loadPreparedFT[file_String, contract_Association] :=
  seededLoadPreparedFT[file, contract];

seededResult = If[seededCustomQ,
  runExample[seededExample, ft2FamilyRequest, seededPrivateHalos],
  runExample[seededExample, None, seededPrivateHalos]];
Print["SEEDED RESULT ", InputForm[seededResult]];
Exit[If[TrueQ[seededResult], 0, 1]];
