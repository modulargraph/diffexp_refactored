(* DiffExp2 FT cutover runner: the M5 ladder.
   Mirrors Scripts/run_ft_stepwise.m but the transport/integration chain is
   DiffExp2 (sector-native, no fitting): per level the in-memory exact
   DiffMatrix is loaded directly (no slice export round-trip), boundary
   values are transported from the anchor 11/23 to both endpoints, and the
   boundary cases run through LineIntegral / EndpointLimitValues / direct
   convolution.  Output: STEPWISE/FINAL rows compatible with the old
   comparator. *)

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

(* The package facade and this script share one strict parser/default set.
   In particular, Cpp is now the release default; the Wolfram recurrence is
   still available only by explicit selection. *)
runnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
If[FailureQ[runnerSettings],
  Print["Invalid Feynman-trick runner environment: ", runnerSettings];
  Exit[2]];

singularMatchPrecondition = runnerSettings["SingularMatchPrecondition"];
recurrenceBackend = runnerSettings["RecurrenceBackend"];
cppBatchEndpointArms = runnerSettings["BatchEndpointArms"];
cppArmThreadBudget = runnerSettings["CppThreads"];
DiffExp2`Transport`Private`$enableSingularMatchPrecondition =
  singularMatchPrecondition;
If[singularMatchPrecondition,
  Print["DE2 singular match precondition enabled"]];
wp = runnerSettings["WorkingPrecision"];
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
requestedEpsilonOrder[level_Integer] := Max[
  epsOrder + level + boundaryExtraOrder + levelEpsilonHalo[level], 1];

(* Old FeynmanTrick/DiffExp prescribed both endpoints and every matrix/IBP
   segmentation factor with a consistent +i delta side.  DiffExp2's config
   is reset at every level, so rebuild that effective list explicitly here;
   otherwise Transport receives an empty Prescriptions record and loses the
   branch sheet after crossing an interior singularity. *)
levelDeltaPrescriptions[var_Symbol, sys_Association, extra_List] := Module[
  {raw, factors},
  raw = Join[{var, 1 - var}, Lookup[sys, "SingularFactors", {}], extra];
  factors = Flatten[Map[Module[{fl = FactorList[Factor[Numerator[Together[#]]]]},
      First /@ Select[fl, !FreeQ[First[#], var] &]] &, raw]];
  factors = DeleteDuplicates[factors,
    TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
      TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &];
  {#, 1} & /@ factors];

anchor = 11/23;
inputPrecision = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
prepCacheRoot = runnerSettings["PrepCacheRoot"];
forcePrepRebuild = runnerSettings["ForcePrepRebuild"];
resumeLadderFile = runnerSettings["ResumeCheckpoint"];
ladderCheckpointDir = runnerSettings["CheckpointDirectory"];
allowStaleLadderCheckpoint = runnerSettings["AllowStaleCheckpoint"];

(* RunFullIteration is FIRE-dominated but independent of DiffExp2's
   transport settings.  Persist both the populated ftData and FIRE's
   in-memory reduction cache so a fresh Wolfram process can resume at the
   transport ladder.  The key covers topology/sequence/dimension and every
   FeynmanTrick source file; edits to this runner or DiffExp2 intentionally
   do not invalidate preparation.  FT_REBUILD_PREP=1 forces a rebuild. *)
$ftPrepCacheVersion = 1;
$ftPrepSourceFingerprint = Hash[
  ({#, FileHash[#, "SHA256"]} & /@ Sort[FileNames["*.m",
    FileNameJoin[{repoRoot, "FeynmanTrick"}], Infinity]]), "SHA256"];
$ftLadderCheckpointVersion = 1;
$ftLadderSourceFingerprint = Hash[
  ({#, FileHash[#, "SHA256"]} & /@ Sort[Join[
    FileNames["*.m", FileNameJoin[{repoRoot, "DiffExp2"}], Infinity],
    Select[FileNames["*", FileNameJoin[{repoRoot, "cpp"}], Infinity],
      FileType[#] === File &],
    {FileNameJoin[{repoRoot, "CMakeLists.txt"}]},
    {ExpandFileName[$InputFileName]}]]), "SHA256"];

preparedFTDataQ[data_] := AssociationQ[data] &&
  IntegerQ[Lookup[data, "NumLevels", None]] &&
  AssociationQ[Lookup[data, "Levels", None]] &&
  AllTrue[Range[data["NumLevels"]], Function[level,
    Module[{ld = Lookup[data["Levels"], level, <||>], masters, mat},
      masters = Lookup[ld, "Masters", {}]; mat = Lookup[ld, "DiffMatrix", {}];
      TrueQ[Lookup[ld, "Computed", False]] && masters =!= {} &&
        MatrixQ[mat] && Dimensions[mat] === {Length[masters], Length[masters]}]]];

requiredReductionKeys[data_] := Flatten[Table[
  Module[{ld = data["Levels"][level], below = data["Levels"][level - 1],
      topo, reqs},
    topo = ld["Topology"];
    reqs = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
      below["Masters"], ld["CombinedPositions"]];
    Map[Function[req, {topo["WorkDirectory"], topo["Name"],
      topo["ProblemNumber"], topo["NumPropagators"],
      PadRight[req["NeededVec"], topo["NumPropagators"], 0]}], reqs]],
  {level, 1, data["NumLevels"]}], 1];

preparedReductionCacheQ[data_, rc_] := AssociationQ[rc] &&
  AllTrue[requiredReductionKeys[data], Function[key,
    Module[{entry = Lookup[rc, Key[key], Missing["NotCached"]]},
      AssociationQ[entry] && KeyExistsQ[entry, "Reduction"] &&
        KeyExistsQ[entry, "Masters"]]]];

ftPrepKey[name_, topology_, sequence_] := Hash[{
  $ftPrepCacheVersion, $Version, name, topology, sequence,
  FeynmanTrick`Private`DimensionExpression[],
  Lookup[FeynmanTrick`Private`$FTConfig,
    {"DimensionVariable", "EpsilonSymbol", "FixedParameterValue"}],
  $ftPrepSourceFingerprint}, "SHA256"];

ftPrepFile[name_, key_] := FileNameJoin[{prepCacheRoot,
  name <> "_" <> IntegerString[Abs[key], 16] <> ".mx"}];

savePreparedFT[file_, key_, data_] := Module[{payload, tmp, ok},
  If[!preparedFTDataQ[data], Return[$Failed, Module]];
  If[!preparedReductionCacheQ[data,
      FeynmanTrick`FIREInterface`Private`$ReductionCache],
    Print["FTPREP CACHE INCOMPLETE; not saving ", file];
    Return[$Failed, Module]];
  If[!DirectoryQ[DirectoryName[file]],
    CreateDirectory[DirectoryName[file], CreateIntermediateDirectories -> True]];
  payload = <|
    "Version" -> $ftPrepCacheVersion, "Key" -> key, "FTData" -> data,
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

loadPreparedFT[file_, key_] := Module[{payload, ok},
  If[!FileExistsQ[file], Return[$Failed, Module]];
  Clear[Global`$FT2PreparedSnapshot];
  ok = Quiet[Check[Get[file]; True, False]];
  If[!TrueQ[ok], Clear[Global`$FT2PreparedSnapshot];
    Return[$Failed, Module]];
  payload = Global`$FT2PreparedSnapshot;
  Clear[Global`$FT2PreparedSnapshot];
  If[!AssociationQ[payload] || Lookup[payload, "Version", None] =!=
      $ftPrepCacheVersion || Lookup[payload, "Key", None] =!= key ||
      !preparedFTDataQ[Lookup[payload, "FTData", None]] ||
      !preparedReductionCacheQ[Lookup[payload, "FTData", None],
        Lookup[payload, "ReductionCache", <||>]],
    Return[$Failed, Module]];
  FeynmanTrick`FIREInterface`Private`$ReductionCache =
    Join[FeynmanTrick`FIREInterface`Private`$ReductionCache,
      Lookup[payload, "ReductionCache", <||>]];
  Print["FTPREP CACHE HIT ", file];
  payload["FTData"]];

saveLadderCheckpoint[file_, payload_] := Module[
  {tmp, saved, wrote, renamed},
  If[ladderCheckpointDir === "", Return[Null, Module]];
  If[!DirectoryQ[DirectoryName[file]],
    CreateDirectory[DirectoryName[file], CreateIntermediateDirectories -> True]];
  saved = Join[payload, <|
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

(* A transport checkpoint is updated after each successful endpoint arm.
   In particular, the lower arm reaches disk before the upper arm starts, so
   an interrupted upper solve can resume without repeating the lower one.
   The same schema accepts either arm (or both) and computes only what is
   missing before replaying boundary assembly.  A boundary checkpoint starts
   directly at the named lower level.  Unversioned/stale files require an
   explicit opt-in: they may contain mathematically valid data, but source
   provenance cannot be proved. *)
loadLadderCheckpoint[file_, name_, data_, prepKey_] := Module[
  {payload, ok, kind, level, levelData, belowData, mastersHere, mastersBelow,
   currentRequests, savedEO, savedFingerprint, stale, needInt, needLo, needHi,
   cachedArms, recordedArms, expectedCharts},
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
  If[KeyExistsQ[payload, "CheckpointVersion"] &&
      payload["CheckpointVersion"] =!= $ftLadderCheckpointVersion,
    Return[ladderCheckpointReject[file, "unsupported checkpoint version"], Module]];
  kind = Lookup[payload, "Kind", "Transport"];
  If[!MemberQ[{"Transport", "Boundary"}, kind],
    Return[ladderCheckpointReject[file, "unknown checkpoint kind"], Module]];
  level = Lookup[payload, "Level", None];
  If[!IntegerQ[level] || level < 1 || level > data["NumLevels"],
    Return[ladderCheckpointReject[file, "invalid resume level"], Module]];
  If[Lookup[payload, "Example", None] =!= name,
    Return[ladderCheckpointReject[file, "example does not match"], Module]];
  If[Lookup[payload, "WorkingPrecision", None] =!= wp,
    Return[ladderCheckpointReject[file, "WorkingPrecision does not match"], Module]];
  If[KeyExistsQ[payload, "RecurrenceBackend"] &&
      payload["RecurrenceBackend"] =!= recurrenceBackend,
    Return[ladderCheckpointReject[file,
      "recurrence backend mode does not match"], Module]];
  If[Lookup[payload, "Anchor", None] =!= anchor,
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
  savedEO = If[kind === "Transport", Lookup[payload, "ExpansionOrder", None],
    Lookup[payload, "SourceExpansionOrder", None]];
  If[!IntegerQ[savedEO] || savedEO < 1,
    Return[ladderCheckpointReject[file, "missing source ExpansionOrder"], Module]];
  (* Transport checkpoints contain order-specific chart series and therefore
     require exact order parity.  Boundary checkpoints contain only endpoint
     Laurent data: a higher-order source may safely seed a lower-order run,
     but lower-order data must never silently downgrade the requested order. *)
  If[kind === "Transport" && expansionOrder =!= savedEO,
    Return[ladderCheckpointReject[file,
      "transport ExpansionOrder does not match"], Module]];
  If[kind === "Boundary" && savedEO < expansionOrder,
    Return[ladderCheckpointReject[file,
      "boundary source ExpansionOrder is lower than requested"], Module]];
  If[kind === "Transport",
    If[KeyExistsQ[payload, "DivisionOrder"] &&
        payload["DivisionOrder"] =!= divisionOrder,
      Return[ladderCheckpointReject[file,
        "DivisionOrder does not match"], Module]];
    If[KeyExistsQ[payload, "RadiusOfConvergence"] &&
        payload["RadiusOfConvergence"] =!= radiusOfConvergence,
      Return[ladderCheckpointReject[file,
        "RadiusOfConvergence does not match"], Module]];
    If[KeyExistsQ[payload, "ValueTransportMode"] &&
        payload["ValueTransportMode"] =!= Environment["DE2_VALUE_TRANSPORT"],
      Return[ladderCheckpointReject[file,
        "DE2_VALUE_TRANSPORT mode does not match"], Module]];
    If[KeyExistsQ[payload, "SingularMatchPrecondition"] &&
        payload["SingularMatchPrecondition"] =!= singularMatchPrecondition,
      Return[ladderCheckpointReject[file,
        "singular match precondition mode does not match"], Module]];
    belowData = data["Levels"][level - 1];
    mastersBelow = belowData["Masters"];
    currentRequests =
      FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
        mastersBelow, levelData["CombinedPositions"]];
    If[Lookup[payload, "Variable", None] =!= levelData["FeynmanParameter"] ||
        Lookup[payload, "MastersBelow", None] =!= mastersBelow ||
        Lookup[payload, "Requests", None] =!= currentRequests,
      Return[ladderCheckpointReject[file,
        "transport level metadata does not match prepared FT data"], Module]];
    If[Lookup[payload, "RequestedEpsilonOrder", None] =!=
        requestedEpsilonOrder[level],
      Return[ladderCheckpointReject[file,
        "requested epsilon window does not match"], Module]];
    If[!AssociationQ[Lookup[payload, "System", None]] ||
        !AssociationQ[Lookup[payload, "Reductions", None]] ||
        !AllTrue[currentRequests,
          KeyExistsQ[payload["Reductions"], #["NeededVec"]] &] ||
        !ListQ[Lookup[payload, "ExtraSingularFactors", None]] ||
        !ListQ[Lookup[payload, "ChartCache", None]] ||
        !AllTrue[{Lookup[payload, "TransportLow", None],
            Lookup[payload, "TransportHigh", None]},
          (# === None || AssociationQ[#]) &],
      Return[ladderCheckpointReject[file,
        "transport payload is incomplete"], Module]];
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
    If[KeyExistsQ[payload, "RequestedEpsilonOrder"] &&
        payload["RequestedEpsilonOrder"] =!=
          requestedEpsilonOrder[level],
      Return[ladderCheckpointReject[file,
        "requested epsilon window does not match"], Module]]];
  savedFingerprint = Lookup[payload, "SourceFingerprint", Missing["Unversioned"]];
  stale = savedFingerprint =!= $ftLadderSourceFingerprint ||
    TrueQ[Lookup[payload, "Tainted", False]];
  If[stale && !allowStaleLadderCheckpoint,
    Return[ladderCheckpointReject[file,
      "source provenance is stale/unversioned; set FT_ALLOW_STALE_LADDER_CHECKPOINT=1 to opt in"],
      Module]];
  If[stale,
    Print["FTLADDER RESUME WARNING stale/unversioned checkpoint explicitly accepted"]];
  Join[payload, <|"Kind" -> kind, "Tainted" -> stale,
    "SourceExpansionOrder" -> savedEO|>]];

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", wp];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds",
  runnerSettings["FIRETimeoutSeconds"]];
persistentFIREStorage = runnerSettings["PersistentFIREStorage"];
fireStorageDirectory = runnerSettings["FIREStorageDirectory"];
FeynmanTrick`SetFTOption["PersistentFIREStorage", persistentFIREStorage];
If[fireStorageDirectory =!= "",
  FeynmanTrick`SetFTOption["FIREStorageDirectory",
    ExpandFileName[fireStorageDirectory]]
];
If[persistentFIREStorage,
  Print["Persistent FIRE storage enabled",
    If[fireStorageDirectory === "", "", ": " <> ExpandFileName[fireStorageDirectory]]]
];

eps = Global`eps;
esC = DiffExp2`EpsSeries`ESCoefficient;
esMn = DiffExp2`EpsSeries`ESMinPower;
esCMx = DiffExp2`EpsSeries`ESCompleteMax;
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

cleanNumber[value_] := Module[{n = N[value, 50], re, im},
  If[!NumericQ[n], Return[ToString[InputForm[value]]]];  (* JSON-safe *)
  re = Re[n]; im = Im[n];
  If[Abs[im] < 10^-30, re, <|"Re" -> re, "Im" -> im|>]];

printRows[example_, level_, masters_, rawES_List, prefactors_] := Module[{},
  Do[Module[{r = rawES[[i]], rowMin},
    rowMin = esMn[r];
    Print["STEPWISE ", ExportString[<|
      "Example" -> example, "Level" -> level, "Master" -> masters[[i]],
      "EpsPrefactor" -> prefactors[[i]], "RawMinPower" -> rowMin,
      "Coefficients" -> Table[{p, cleanNumber[esC[r, p]]},
        {p, rowMin, Min[0, esCMx[r]]}]|> /. x_Rational :> N[x, 50],
      "RawJSON", "Compact" -> True]]],
    {i, Length[masters]}]];

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

runExample[name_String] := Module[
  {topology, sequence, prepKey, prepFile, ftData, outputDir, nLevels,
   boundaryOrder, deepBoundary, currentBCs, currentPrefactors,
   resumeCheckpoint = None, startLevel, finalRaw = None},
  Print["EXAMPLE ", name];
  FeynmanTrick`SetFTOption["DimensionExpression", FTExampleDimension[name]];
  topology = FTExampleTopology[name, "step"];
  If[topology === $Failed, Return[$Failed]];
  sequence = FTExampleSequence[name];
  prepKey = ftPrepKey[name, topology, sequence];
  prepFile = ftPrepFile[name, prepKey];
  ftData = If[forcePrepRebuild, $Failed, loadPreparedFT[prepFile, prepKey]];
  If[ftData === $Failed,
    Print["FTPREP CACHE MISS ", prepFile];
    ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
      topology, sequence, {}];
    outputDir = FileNameJoin[{$TemporaryDirectory,
      "FT2_" <> name <> "_" <> ToString[$ProcessID]}];
    If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
    CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
    ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
    If[preparedFTDataQ[ftData], savePreparedFT[prepFile, prepKey, ftData]]];
  If[ftData === $Failed, Return[$Failed]];
  nLevels = ftData["NumLevels"];
  If[resumeLadderFile =!= "",
    resumeCheckpoint = loadLadderCheckpoint[ExpandFileName[resumeLadderFile],
      name, ftData, prepKey];
    If[resumeCheckpoint === $Failed, Return[$Failed]]];
  If[AssociationQ[resumeCheckpoint],
    startLevel = resumeCheckpoint["Level"];
    currentBCs = resumeCheckpoint["BoundaryValues"];
    currentPrefactors = resumeCheckpoint["BoundaryPrefactors"];
    Print["FTLADDER RESUME kind=", resumeCheckpoint["Kind"],
      " level=", startLevel,
      " savedEO=", resumeCheckpoint["SourceExpansionOrder"],
      " requestedLowerEO=", expansionOrder],
    startLevel = nLevels;
    boundaryOrder = requestedEpsilonOrder[nLevels];
    deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
      ftData, boundaryOrder];
    If[!AssociationQ[deepBoundary], Return[$Failed]];
    (* coefficients are the only numerics (tags stay exact): numericize the
       deep boundary at 2x WP so the chain runs at arbitrary precision instead
       of exact-symbolic (Log/Gamma giants grind the Laurent-field algebra) *)
    currentBCs = N[deepBoundary["BoundaryValues"], inputPrecision];
    currentPrefactors = deepBoundary["EpsPrefactors"];
    printRows[name, nLevels, ftData["Levels"][nLevels]["Masters"],
      Table[DiffExp2`EpsSeries`ESShift[
        DiffExp2`EpsSeries`ESNew[0, currentBCs[[i]]], -currentPrefactors[[i]]],
        {i, Length[currentBCs]}],
      currentPrefactors]];

  Module[{abortRes},
  abortRes = Catch[
  Do[Module[
    {levelData = ftData["Levels"][level], levelBelow = ftData["Levels"][level - 1],
     var, A, sys, mastersBelow, mastersHere, requests, neededVecs, reductions,
     extraFacs, rawES, rawMin, shift, kmaxAvail, nextReq, needTop, ftEps, dimVar,
     dimExpr, normalizeFT, trLoCache, trHiCache, chartCache,
     resumeTransport, levelExpansionOrder, needInt, needLo, needHi,
     transportCheckpointFile, saveTransportProgress, completedArms,
     transportSys = None, planLo = None, planHi = None, armReq,
     loPlanCharts, hiPlanCharts, armRounds, armBatchResult,
     armUniqueCharts, armCacheCapacity},
    var = levelData["FeynmanParameter"];
    (* normalize FT-layer symbols at the seam: dimension d -> 2-2eps form,
       FT epsilon symbol -> the DiffExp2 canonical Global`eps *)
    ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
    dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
    dimExpr = FeynmanTrick`Private`DimensionExpression[];
    normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
      ftEps -> Global`eps);
    A = normalizeFT[levelData["DiffMatrix"]];
    Print["LEVEL ", level, " var=", var, " d=", Length[A]];
    mastersHere = levelData["Masters"];
    mastersBelow = levelBelow["Masters"];
    resumeTransport = AssociationQ[resumeCheckpoint] &&
      resumeCheckpoint["Kind"] === "Transport" &&
      level === resumeCheckpoint["Level"];
    levelExpansionOrder = If[resumeTransport,
      resumeCheckpoint["SourceExpansionOrder"], expansionOrder];
    If[resumeTransport,
      sys = resumeCheckpoint["System"];
      If[Lookup[sys, "Variable", None] =!= var ||
          Lookup[sys, "Matrix", None] =!= A,
        Print["FTLADDER RESUME REJECT: saved system does not match level matrix"];
        Throw[$Failed, "FT2Abort"]];
      requests = resumeCheckpoint["Requests"];
      reductions = resumeCheckpoint["Reductions"];
      extraFacs = resumeCheckpoint["ExtraSingularFactors"];
      DiffExp2`Solve`ClearSolveCaches[],
      sys = catch2[DiffExp2`API`LoadSystem[
        <|"Matrix" -> A, "Variable" -> var|>]];
      If[FailureQ[sys], Print["LOAD FAIL ", sys]; Return[$Failed, Module]];
      requests = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
        mastersBelow, levelData["CombinedPositions"]];
      neededVecs = DeleteDuplicates[#["NeededVec"] & /@ requests];
      reductions = FeynmanTrick`FIREInterface`ReduceIntegrals[
        levelData["Topology"], neededVecs];
      If[reductions === $Failed, Print["FIRE FAIL"]; Return[$Failed, Module]];
      extraFacs = normalizeFT[
        FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
          ftData, level]]];
    (* configure DiffExp2 for this level *)
    catch2[DiffExp2`Config`LoadConfiguration[{
      "WorkingPrecision" -> wp, "ExpansionOrder" -> levelExpansionOrder,
      "EpsilonOrder" -> (esCMxLevel = requestedEpsilonOrder[level]),
      "DivisionOrder" -> divisionOrder,
      "RadiusOfConvergence" -> radiusOfConvergence,
      (* The restored classic predivision planner couples placement and
         matching: adjacent regular segments meet at +1/k and -1/k. *)
      "StepDivisionOrder" -> stepDivisionOrder,
      "RecurrenceBackend" -> recurrenceBackend,
      "DeltaPrescriptions" -> levelDeltaPrescriptions[var, sys, extraFacs],
      "Variables" -> {}}]];
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
        "MastersHere" -> mastersHere, "MastersBelow" -> mastersBelow,
        "Requests" -> requests,
        "Reductions" -> AssociationMap[normalizeFT, reductions],
        "ExtraSingularFactors" -> extraFacs, "ChartCache" -> chartCache,
        "TransportLow" -> trLoCache, "TransportHigh" -> trHiCache,
        "CompletedArms" -> completedArms,
        "Anchor" -> anchor, "WorkingPrecision" -> wp,
        "DivisionOrder" -> divisionOrder,
        "RadiusOfConvergence" -> radiusOfConvergence,
        "ValueTransportMode" -> Environment["DE2_VALUE_TRANSPORT"],
        "RecurrenceBackend" -> recurrenceBackend,
        "SingularMatchPrecondition" -> singularMatchPrecondition,
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
        Environment["DE2_VALUE_TRANSPORT"] =!= "1" &&
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
        " bcMinPow=", Min @@ currentPrefactors, "..", Max @@ currentPrefactors,
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
    If[MemberQ[rawES, $Failed], Throw[$Failed, "FT2Abort"]];
    rawES = DiffExp2`EpsSeries`ESTrim /@ rawES;
    printRows[name, level - 1, mastersBelow, rawES,
      ConstantArray[0, Length[mastersBelow]]];
    (* shift to finite for the next level's transport *)
    rawMin = Min[esMn /@ rawES];
    shift = Max[0, -rawMin];
    kmaxAvail = esCMx /@ rawES;
    If[level > 1,
      nextReq = requestedEpsilonOrder[level - 1];
      needTop = nextReq - shift;
      If[AnyTrue[kmaxAvail, # < needTop &],
        Print["FTLADDER INCOMPLETE level=", level - 1,
          " requiredTop=", needTop, " availableTops=", kmaxAvail,
          " shift=", shift];
        Throw[$Failed, "FT2Abort"]];
      (* Numericize only at a genuine level handoff: exact Log-trees from
         tile antiderivatives otherwise compound into symbolic giants that
         grind the next level's recursion (meprec storms). *)
      currentBCs = Table[Module[{r = rawES[[i]]},
        Table[N[esC[r, k], inputPrecision], {k, -shift, needTop}]],
        {i, Length[rawES]}];
      If[!AllTrue[currentBCs, Length[#] === nextReq + 1 &],
        Print["FTLADDER INTERNAL ERROR: nonuniform boundary width at level ",
          level - 1, " lengths=", Length /@ currentBCs,
          " expected=", nextReq + 1];
        Throw[$Failed, "FT2Abort"]];
      currentPrefactors = ConstantArray[shift, Length[rawES]];
      If[ladderCheckpointDir =!= "",
        saveLadderCheckpoint[FileNameJoin[{ladderCheckpointDir,
            name <> "_level" <> ToString[level - 1] <> "_boundary.mx"}], <|
          "Kind" -> "Boundary", "Example" -> name, "Level" -> level - 1,
          "PrepKey" -> prepKey, "BoundaryValues" -> currentBCs,
          "BoundaryPrefactors" -> currentPrefactors,
          "MastersHere" -> mastersBelow, "Anchor" -> anchor,
          "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
          "RecurrenceBackend" -> recurrenceBackend,
          "BoundaryExtraOrder" -> boundaryExtraOrder,
          "LevelEpsilonHalos" -> levelEpsilonHalos,
          "SourceExpansionOrder" -> levelExpansionOrder,
          "RequestedEpsilonOrder" -> nextReq,
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
    finalRaw = rawES],
    {level, startLevel, 1, -1}], "FT2Abort"];
  If[abortRes === $Failed, Return[$Failed]];
  If[abortRes === "Stopped", Return[True]]];
  If[finalRaw === None || MemberQ[finalRaw, $Failed], Return[$Failed]];

  Print["FINAL ", ExportString[<|
    "Example" -> name,
    "Finite" -> cleanNumber[esC[finalRaw[[1]], 0]],
    "RawMinPower" -> esMn[finalRaw[[1]]]|> /. x_Rational :> N[x, 50],
    "RawJSON", "Compact" -> True]];
  True];

(* Let the focused checkpoint tests load these definitions without starting
   FIRE or terminating their Wolfram kernel. *)
If[envOrDefault["FT_RUNNER_DEFINITIONS_ONLY", "0"] =!= "1",
  requested = StringTrim /@
    StringSplit[envOrDefault["FT_EXAMPLES", "bubble"], ","];
  Do[
    If[runExample[name] === $Failed, Print["FAILED ", name]; Quit[1]],
    {name, requested}];
  Quit[0]];
