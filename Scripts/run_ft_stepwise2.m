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

(* The package facade and this script share one strict parser/default set.
   In particular, Cpp is now the release default; the Wolfram recurrence is
   still available only by explicit selection. *)
runnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];
If[FailureQ[runnerSettings],
  Print["Invalid Feynman-trick runner environment: ", runnerSettings];
  Exit[2]];

FeynmanTrick`SetFTOption["FIREPath", runnerSettings["FIREPath"]];

singularMatchPrecondition = runnerSettings["SingularMatchPrecondition"];
recurrenceBackend = runnerSettings["RecurrenceBackend"];
cppBatchEndpointArms = runnerSettings["BatchEndpointArms"];
cppArmThreadBudget = runnerSettings["CppThreads"];
deltaPrescriptionSign = runnerSettings["DeltaPrescriptionSign"];
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
levelDeltaPrescriptions[var_Symbol, sys_Association, extra_List] := Module[
  {raw, factors, projectedExtra},
  projectedExtra = DiffExp2`Transport`EpsilonZeroSingularFactors[extra, var];
  raw = Join[{var, 1 - var}, Lookup[sys, "SingularFactors", {}],
    projectedExtra];
  factors = Flatten[Map[Module[{fl = FactorList[Factor[Numerator[Together[#]]]]},
      First /@ Select[fl, !FreeQ[First[#], var] &]] &, raw]];
  factors = DeleteDuplicates[factors,
    TrueQ[PossibleZeroQ[Expand[#1 - #2]]] ||
      TrueQ[PossibleZeroQ[Expand[#1 + #2]]] &];
  {#, deltaPrescriptionSign} & /@ factors];

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
   FeynmanTrick source file, while each exact reduction key embeds the
   setup-time FIRE fingerprint.  Edits to this runner or DiffExp2 intentionally
   do not invalidate preparation.  FT_REBUILD_PREP=1 forces a rebuild. *)
$ftPrepCacheVersion = 2;
$ftPrepSourceFingerprint = Hash[
  ({#, FileHash[#, "SHA256"]} & /@ Sort[FileNames["*.m",
    FileNameJoin[{repoRoot, "FeynmanTrick"}], Infinity]]), "SHA256"];
$ftLadderCheckpointVersion = 2;
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

ftPrepKey[name_, topology_, sequence_] := Hash[{
  $ftPrepCacheVersion, $Version, name, topology, sequence,
  FeynmanTrick`Private`DimensionExpression[],
  Lookup[FeynmanTrick`Private`$FTConfig,
    {"DimensionVariable", "EpsilonSymbol", "FixedParameterValue",
      "AutoDetectRestrictions"}],
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[],
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

(* The matrix fixes only a relative epsilon gauge.  Capture that exact gauge
   once so the full-ladder planner and the runtime normalization use the same
   pole-free basis without repeating MatrixPoleOrders/FindEpsPrefactors. *)
ft2RelativeEpsilonGauge[matrix_, epsSymbol_Symbol] := Module[
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
    suppliedGauge_:Automatic] := Module[
  {d, widths, inputTop, gauge, rawPoleOrders, canonical, relative,
   commonOffset, effective, boundaryShifts, normalized,
   normalizedPoleOrders, shiftedBoundary, record},
  d = Length[matrix];
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
  widths = Length /@ boundaryValues;
  If[MemberQ[widths, 0] || Length[DeleteDuplicates[widths]] =!= 1,
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "epsilon boundary rows must have one nonempty common window",
      "Widths" -> widths|>], Module]];
  inputTop = First[widths] - 1;
  gauge = If[suppliedGauge === Automatic,
    ft2RelativeEpsilonGauge[matrix, epsSymbol], suppliedGauge];
  If[FailureQ[gauge], Return[gauge, Module]];
  If[!AssociationQ[gauge] ||
      Lookup[Lookup[gauge, "Record", <||>], "InputMatrixHash", None] =!=
        Hash[matrix, "SHA256"] ||
      !TrueQ[Lookup[Lookup[gauge, "Record", <||>], "PoleFree", False]] ||
      Lookup[gauge, "Identity", None] =!=
        ft2CanonicalIdentity["ft2-relative-epsilon-gauge-",
          Lookup[gauge, "Record", None]],
    Return[Failure["FeynmanTrickEpsilonBasis", <|
      "Detail" -> "supplied relative epsilon gauge does not match the level matrix"|>],
    Module]];
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
   prefactors, matrix, widths},
  system = Lookup[payload, "System", None];
  values = Lookup[payload, "BoundaryValues", None];
  prefactors = Lookup[payload, "BoundaryPrefactors", None];
  If[!AssociationQ[record] || !AssociationQ[system] || !ListQ[values] ||
      !ListQ[prefactors] || !AllTrue[values, ListQ], Return[False, Module]];
  widths = Length /@ values;
  If[MemberQ[widths, 0] || Length[DeleteDuplicates[widths]] =!= 1,
    Return[False, Module]];
  matrix = Lookup[system, "Matrix", None];
  TrueQ[Lookup[record, "Schema", None] ===
      "FeynmanTrick.EpsilonBasis/v1"] &&
    TrueQ[Lookup[record, "PoleFree", False]] &&
    Lookup[record, "Prefactors", None] === prefactors &&
    Lookup[record, "CompleteMax", None] === First[widths] - 1 &&
    Lookup[record, "NormalizedMatrixHash", None] ===
      Hash[matrix, "SHA256"]];

(* A transport checkpoint is updated after each successful endpoint arm.
   In particular, the lower arm reaches disk before the upper arm starts, so
   an interrupted upper solve can resume without repeating the lower one.
   The same schema accepts either arm (or both) and computes only what is
   missing before replaying boundary assembly.  A boundary checkpoint starts
   directly at the named lower level.  Unversioned/stale files require an
   explicit opt-in: they may contain mathematically valid data, but source
   provenance cannot be proved. *)
loadLadderCheckpoint[file_, name_, data_, prepKey_, nativePlan_:None] := Module[
  {payload, ok, kind, level, levelData, belowData, mastersHere, mastersBelow,
   currentRequests, savedEO, savedFingerprint, stale, needInt, needLo, needHi,
   cachedArms, recordedArms, expectedCharts, boundaryWidths, boundaryShift,
   requiredRaw, preservedRaw, preservedSource, nativeRecord,
   nativePlanRecord, nativePlanIdentity},
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
  If[kind === "Transport" && recurrenceBackend === "Cpp",
    Return[ladderCheckpointReject[file,
      "legacy partial-arm Transport snapshots cannot resume the retained native observable batch; resume from a numeric Boundary checkpoint"],
      Module]];
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
  If[Lookup[payload, "DeltaPrescriptionSign", 1] =!=
      deltaPrescriptionSign,
    Return[ladderCheckpointReject[file,
      "delta prescription sign does not match"], Module]];
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
      FeynmanTrick`LevelReduction`BoundaryRequestRecords[
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
        !ft2EpsilonBasisCheckpointQ[payload] ||
        !AssociationQ[Lookup[payload, "Reductions", None]] ||
        !AllTrue[currentRequests,
          KeyExistsQ[payload["Reductions"], #["NeededVec"]] &] ||
        !ListQ[Lookup[payload, "ExtraSingularFactors", None]] ||
        !ListQ[Lookup[payload, "ChartCache", None]] ||
        !AllTrue[{Lookup[payload, "TransportLow", None],
            Lookup[payload, "TransportHigh", None]},
          (# === None || AssociationQ[#]) &],
      Return[ladderCheckpointReject[file,
        "transport payload is incomplete or has inconsistent epsilon-basis metadata"],
        Module]];
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
    If[recurrenceBackend === "Cpp",
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
          Lookup[payload, "RequiredRawTop", None] =!= requiredRaw ||
          Lookup[payload, "RequestedEpsilonOrder", None] =!= requiredRaw ||
          !IntegerQ[preservedRaw] || preservedRaw < requiredRaw ||
          !IntegerQ[preservedSource] ||
          preservedSource =!= preservedRaw + boundaryShift ||
          preservedSource =!= First[boundaryWidths] - 1 ||
          !ft2NativeCheckpointRecordQ[nativeRecord] ||
          nativeRecord["RequiredRawTop"] =!= requiredRaw ||
          nativeRecord["DeliverableCompleteMax"] =!= preservedRaw ||
          (ft2NativeEpsilonPlanQ[nativePlan] &&
            (!ft2NativeEpsilonExecutionRecordQ[nativePlanRecord,
                nativePlanIdentity, nativePlan] ||
              Lookup[nativeRecord, "NativeEpsilonPlanIdentity", None] =!=
                nativePlanIdentity)),
        Return[ladderCheckpointReject[file,
          "native Boundary checkpoint has an inconsistent required floor, preserved source width, shift, or observable-batch identity"],
          Module]],
      If[KeyExistsQ[payload, "RequestedEpsilonOrder"] &&
          payload["RequestedEpsilonOrder"] =!=
            requestedEpsilonOrder[level],
        Return[ladderCheckpointReject[file,
          "requested epsilon window does not match"], Module]]]];
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

(* Exact, runner-local construction for the retained native observable seam.
   Keeping this data preparation separate from runExample makes it possible
   to test the production contract without starting FIRE. *)
ft2NativeFailure[detail_String, data_:<||>] :=
  Failure["FeynmanTrickNativeBoundary", Join[<|"Detail" -> detail|>, data]];

ft2CanonicalIdentity[prefix_String, value_] := prefix <>
  IntegerString[Hash[value, "SHA256"], 16, 64];

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
      "DeliverableCompleteMax", "RequiredRawTop",
      "CoefficientHalo", "IntegrationHalo"},
    None], IntegerQ] &&
  With[{source = record["SourceCompleteMax"],
      target = record["TargetCompleteMax"],
      deliverable = record["DeliverableCompleteMax"],
      required = record["RequiredRawTop"],
      coefficientHalo = record["CoefficientHalo"],
      integrationHalo = record["IntegrationHalo"]},
    coefficientHalo >= 0 && MemberQ[{0, 1}, integrationHalo] &&
      target === source - coefficientHalo &&
      source >= target >= deliverable >= required];

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
    normalize_] := Module[
  {requests = Lookup[batch, "BoundaryRequests", None],
   vectors = Lookup[batch, "CoefficientVectors", None],
   batchKey = Lookup[batch, "Key", Missing["NoBatchKey"]],
   payloadKey = Lookup[batch, "PayloadKey", Missing["NoPayloadKey"]],
   dimension = Length[currentPrefactors], entries},
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
      base = MapThread[Together[normalize[#1]/epsSymbol^#2] &,
        {raw, currentPrefactors}];
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
    halos_List, normalize_, suppliedBatches_:Automatic] := Module[
  {levels = Lookup[ftData, "Levels", None], nLevels, previousRequired,
   levelRecords = {}, runtimeLevels = <||>, levelData, matrix, gauge,
   batch, entries, active, entryLosses, intrinsicLoss, userFloor,
   required, record, identity},
  nLevels = Lookup[ftData, "NumLevels", If[AssociationQ[levels],
    Length[Select[Keys[levels], IntegerQ[#] && # > 0 &]], None]];
  If[!AssociationQ[levels] || !IntegerQ[nLevels] || nLevels < 1 ||
      epsilonOrder < 0 || !AllTrue[halos, IntegerQ[#] && # >= 0 &],
    Return[ft2NativeFailure[
      "native epsilon preplanner received invalid levels, epsilon order, or level halos"],
      Module]];
  previousRequired = epsilonOrder;
  Do[
    If[!KeyExistsQ[levels, level],
      Return[ft2NativeFailure[
        "native epsilon preplanner is missing a positive FT level",
        <|"Level" -> level|>], Module]];
    levelData = levels[level];
    matrix = normalize[Lookup[levelData, "DiffMatrix", None]];
    gauge = ft2RelativeEpsilonGauge[matrix, Global`eps];
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
      Global`eps, normalize];
    If[FailureQ[entries], Return[entries, Module]];
    active = Select[entries,
      !TrueQ[Lookup[#, "ProvenZero", False]] &];
    entryLosses = Association@Map[Function[entry,
      entry["MasterIndex"] -> Max[0,
        If[entry["Case"] === "integrate", 1, 0] -
          entry["MinimumEpsilonShift"]]], active];
    intrinsicLoss = If[entryLosses === <||>, 0,
      Max[Values[entryLosses]]];
    userFloor = ft2UserRawFloor[epsilonOrder, halos, level];
    required = Max[userFloor, previousRequired + intrinsicLoss];
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
      "IntrinsicLoss" -> intrinsicLoss,
      "UserRawFloor" -> userFloor,
      "RequiredOutputRawTop" -> previousRequired,
      "RequiredRawTop" -> required|>;
    AppendTo[levelRecords, record];
    AssociateTo[runtimeLevels, level -> <|
      "Record" -> record, "Gauge" -> gauge, "Batch" -> batch,
      "RelativeEntries" -> entries,
      "RequiredOutputRawTop" -> previousRequired,
      "RequiredRawTop" -> required|>];
    previousRequired = required,
    {level, 1, nLevels}];
  record = <|
    "Schema" -> "FeynmanTrick.NativeEpsilonPlan/v1",
    "EpsilonOrder" -> epsilonOrder,
    "LevelEpsilonHalos" -> halos,
    "NumLevels" -> nLevels,
    "Levels" -> levelRecords,
    "DeepRequiredRawTop" -> previousRequired|>;
  identity = ft2CanonicalIdentity["ft2-native-epsilon-plan-", record];
  <|"Schema" -> record["Schema"], "Identity" -> identity,
    "Record" -> record, "Levels" -> runtimeLevels,
    "NumLevels" -> nLevels,
    "DeepRequiredRawTop" -> previousRequired|>];

ft2NativeEpsilonPlanQ[plan_] := Module[
  {record, levels, nLevels, levelRecords, deepRequired},
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
  If[!ListQ[levelRecords] || Length[levelRecords] =!= nLevels ||
      Sort[Keys[levels]] =!= Range[nLevels] ||
      !IntegerQ[deepRequired] ||
      Lookup[plan, "DeepRequiredRawTop", None] =!= deepRequired ||
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
        Lookup[runtime, "RequiredRawTop", None] ===
          Lookup[saved, "RequiredRawTop", None] &&
        IntegerQ[Lookup[saved, "RequiredOutputRawTop", None]] &&
        IntegerQ[Lookup[saved, "RequiredRawTop", None]] &&
        IntegerQ[Lookup[saved, "IntrinsicLoss", None]] &&
        Lookup[saved, "IntrinsicLoss", -1] >= 0]]]
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
    downstreamFiniteTop_Integer] := Module[
  {widths, availableSourceMax, active, nonDirect,
   coefficientShift, coefficientHalo, integrationHalo, targetMax,
   deliverableMax, outputMins, capacityByMaster, activeCapacities},
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
  coefficientHalo = Max[0, -coefficientShift];
  integrationHalo = ft2NativeIntegrationHalo[active];
  (* downstreamFiniteTop is an independently required raw floor F.  A pole
     that happens to appear (or cancel) never buys completeness.  Preserve
     the honest source reservoir instead: with available work edge S,
     Prepare's maximal public target is T=S-HC.  Per-entry exact capacities
     are S+s for direct/limits and S+s-1 for native integration. *)
  targetMax = availableSourceMax - coefficientHalo;
  outputMins = Association@Map[Function[entry,
    entry["MasterIndex"] -> If[entry["Case"] === "integrate",
      entry["MinimumEpsilonShift"] - integrationHalo,
      entry["MinimumEpsilonShift"]]], active];
  capacityByMaster = Association@Map[Function[entry,
    entry["MasterIndex"] -> (availableSourceMax +
      entry["MinimumEpsilonShift"] -
      If[entry["Case"] === "integrate", 1, 0])], active];
  activeCapacities = Values[capacityByMaster];
  deliverableMax = If[activeCapacities === {}, availableSourceMax,
    Min[Prepend[activeCapacities, availableSourceMax]]];
  If[downstreamFiniteTop < 0 || targetMax < 0 ||
      deliverableMax < downstreamFiniteTop,
    Return[ft2NativeFailure[
      "source epsilon depth cannot cover the downstream raw boundary window",
      <|"AvailableSourceCompleteMax" -> availableSourceMax,
        "RequiredRawTop" -> downstreamFiniteTop,
        "CoefficientHalo" -> coefficientHalo,
        "IntegrationHalo" -> integrationHalo,
        "TargetCompleteMax" -> targetMax,
        "DeliverableCompleteMax" -> deliverableMax,
        "CapacityByMaster" -> capacityByMaster|>], Module]];
  <|"AvailableSourceCompleteMax" -> availableSourceMax,
    "SourceCompleteMax" -> availableSourceMax,
    "CoefficientMinimumShift" -> coefficientShift,
    "CoefficientHalo" -> coefficientHalo,
    "IntegrationHalo" -> integrationHalo,
    "TargetCompleteMax" -> targetMax,
    "DeliverableCompleteMax" -> deliverableMax,
    "OutputMinimums" -> outputMins,
    "CapacityByMaster" -> capacityByMaster,
    "DownstreamFiniteTop" -> downstreamFiniteTop,
    "DownstreamRawTop" -> downstreamFiniteTop|>];

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
  DiffExp2`Transport`SegmentLine[sys, path];
ft2NativePrepare[sys_, boundary_, lower_, upper_, coefficientVectors_,
    physicalVar_, targetMax_, threads_] :=
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    sys, boundary, lower, upper, "Threads" -> threads,
    "Integrands" -> {coefficientVectors, physicalVar},
    "TargetCompleteMax" -> targetMax];
ft2NativeRun[atlas_, observables_, physicalVar_] :=
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, observables, physicalVar];
ft2NativeExport[batch_, digits_] :=
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    batch, digits];
ft2NativeReleaseBatch[batch_] :=
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[batch];
ft2NativeReleaseAtlas[atlas_] :=
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[atlas];

ft2RunNativeBoundaryDispatch[sys_Association, currentBCs_List,
    entries_List, ledger_Association, physicalVar_Symbol, anchor_,
    extraSingularFactors_List, deltaPrescriptions_List, threads_Integer,
    outputDigits_Integer, nativePlanIdentity_:None] :=
 Module[
  {deliverableMax = ledger["DeliverableCompleteMax"],
   integrationHalo = ledger["IntegrationHalo"], directRequiredTop =
     ledger["DeliverableCompleteMax"], provenZeroEntries, activeEntries,
   directEntries, nonDirectEntries, nativeEntries, directValues = <||>,
   nativeValues = <||>, zeroValues = <||>, transportSystem, lowerPlan,
   upperPlan, paddedBoundary, observables, atlas = None, batch = None,
   exported = None, exportedResults, exportedValues,
   cleanupResult = None, result,
   dispatchTag = Unique["ft2NativeDispatch"], values, releaseOKQ,
   directRules, masterIndices, batchKeys, batchPayloadKeys,
   identities, checkpointIdentities, requestIdentities,
   coefficientIdentities,
   prescriptionIdentity, extraFactorsIdentity, atlasPlanIdentity = None,
   nativeBatchPayloadIdentity, publishedBatchQ},
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
          StringLength[nativePlanIdentity] > 0)),
    Return[ft2NativeFailure[
      "native boundary dispatch received inconsistent master, batch, case, or identity metadata",
      <|"MasterIndices" -> masterIndices, "BatchKeys" -> batchKeys,
        "BatchPayloadKeys" -> batchPayloadKeys|>], Module]];
  prescriptionIdentity = ft2CanonicalIdentity[
    "ft2-delta-prescriptions-", deltaPrescriptions];
  extraFactorsIdentity = ft2CanonicalIdentity[
    "ft2-extra-singular-factors-", extraSingularFactors];
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
    lowerPlan = catch2[
      ft2NativeSegmentLine[transportSystem, {anchor, 0}]];
    upperPlan = catch2[
      ft2NativeSegmentLine[transportSystem, {anchor, 1}]];
    If[FailureQ[lowerPlan] || FailureQ[upperPlan],
      Return[First[Select[{lowerPlan, upperPlan}, FailureQ]], Module]];
    paddedBoundary = If[integrationHalo === 0,
      DiffExp2`EpsSeries`ESNew[0, #] & /@ currentBCs,
      DiffExp2`EpsSeries`ESNew[-integrationHalo,
        Join[ConstantArray[0, integrationHalo], #]] & /@ currentBCs];
    observables = Map[Function[entry, Module[{observable},
      observable = <|"Operation" -> entry["Case"],
        "Identity" -> entry["Identity"],
        "CheckpointIdentity" -> entry["CheckpointIdentity"],
        "CoefficientVector" -> entry["CoefficientVector"],
        "Epsilon" -> <|"Min" -> entry["OutputMin"],
          "Max" -> deliverableMax,
          "RequiredCompleteMax" -> deliverableMax|>|>;
      If[entry["Case"] === "integrate",
        Append[observable, "TailPolicy" -> "stored"], observable]]],
      nativeEntries];
    result = Catch[Internal`WithLocalSettings[
      Null,
      atlas = catch2[ft2NativePrepare[transportSystem, paddedBoundary,
        lowerPlan, upperPlan, Lookup[nativeEntries, "CoefficientVector"],
        physicalVar, ledger["TargetCompleteMax"], threads]];
      If[FailureQ[atlas] || !AssociationQ[atlas] ||
          Lookup[atlas, "Type", None] =!=
            "DiffExp2NativeRegularIndependentArmAtlas" ||
          !StringQ[Lookup[atlas, "PlanCheckpointIdentity", None]] ||
          StringLength[atlas["PlanCheckpointIdentity"]] === 0,
        Throw[If[FailureQ[atlas], atlas,
          ft2NativeFailure["native atlas preparation returned a malformed result",
            <|"Result" -> atlas|>]], dispatchTag]];
      atlasPlanIdentity = atlas["PlanCheckpointIdentity"];
      nativeBatchPayloadIdentity = ft2CanonicalIdentity[
        "ft2-native-observable-payload-",
        {nativeBatchPayloadIdentity, atlasPlanIdentity,
          Map[KeyTake[#, {"Operation", "Identity",
              "CheckpointIdentity", "Epsilon", "TailPolicy"}] &,
            observables]}];
      batch = catch2[ft2NativeRun[atlas, observables, physicalVar]];
      If[FailureQ[batch] || !AssociationQ[batch] ||
          Lookup[batch, "Type", None] =!=
            "DiffExp2NativeTransportObservableBatch" ||
          Lookup[batch, "Atlas", None] =!= atlas,
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
    exportedValues = Map[If[esCMx[#] > deliverableMax,
        DiffExp2`EpsSeries`ESTruncate[#, deliverableMax], #] &,
      Lookup[exportedResults, "Value"]];
    If[AnyTrue[exportedValues, esCMx[#] < deliverableMax &],
      Return[ft2NativeFailure[
        "native observable export is incomplete at the common retained raw edge"],
        Module]];
    nativeValues = AssociationThread[
      Lookup[nativeEntries, "MasterIndex"],
      exportedValues]];
  values = Map[Function[masterIndex,
      Lookup[Join[directValues, zeroValues, nativeValues], masterIndex,
        Missing["MissingBoundaryValue", masterIndex]]],
    Lookup[entries, "MasterIndex"]];
  If[AnyTrue[values, MissingQ],
    Return[ft2NativeFailure[
      "native/direct result merge did not cover every lower master",
      <|"Missing" -> Cases[values, _Missing]|>], Module]];
  <|"Values" -> values,
    "NativeBatchCalls" -> If[nativeEntries === {}, 0, 1],
    "NativeMarches" -> If[AssociationQ[exported],
      Lookup[exported, "NativeMarches", Missing["NotReported"]], 0],
    "CompatibilityExports" -> If[AssociationQ[exported],
      Lookup[exported, "CompatibilityExports", Length[nativeEntries]], 0],
    "CheckpointRecord" -> <|
      "Schema" -> "FeynmanTrick.NativeObservableBatch/v1",
      "BatchKey" -> First[batchKeys],
      "BatchPayloadKey" -> First[batchPayloadKeys],
      "RequestIdentities" -> Lookup[entries, "RequestIdentity"],
      "CoefficientIdentities" -> Lookup[entries, "CoefficientIdentity"],
      "ObservableIdentities" -> If[nativeEntries === {}, {},
        Lookup[nativeEntries, "Identity"]],
      "ObservableCheckpointIdentities" ->
        If[nativeEntries === {}, {},
          Lookup[nativeEntries, "CheckpointIdentity"]],
      "DeltaPrescriptions" -> deltaPrescriptions,
      "DeltaPrescriptionIdentity" -> prescriptionIdentity,
      "ExtraSingularFactorsIdentity" -> extraFactorsIdentity,
      "AtlasPlanIdentity" -> atlasPlanIdentity,
      "NativeBatchPayloadIdentity" -> nativeBatchPayloadIdentity,
      "NativeEpsilonPlanIdentity" -> nativePlanIdentity,
      "SourceCompleteMax" -> ledger["SourceCompleteMax"],
      "TargetCompleteMax" -> ledger["TargetCompleteMax"],
      "DeliverableCompleteMax" -> deliverableMax,
      "RequiredRawTop" -> ledger["DownstreamRawTop"],
      "CoefficientHalo" -> ledger["CoefficientHalo"],
      "IntegrationHalo" -> integrationHalo|>|>];

runExample[name_String] := Module[
  {topology, sequence, prepKey, prepFile, ftData, outputDir, nLevels,
   boundaryOrder, deepBoundary, currentBCs, currentPrefactors,
   resumeCheckpoint = None, startLevel, finalRaw = None, ftEps, dimVar,
   dimExpr, normalizeFT, nativeEpsilonPlan = None,
   nativeEpsilonExecution = None, initialDeepPrefactors,
   deepRelativeGauge, deepGaugeOffset, deepBoundaryWindow,
   deepRequiredSourceCompleteMax, deepBoundaryDeficit},
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
  ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
  dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  dimExpr = FeynmanTrick`Private`DimensionExpression[];
  normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
    ftEps -> Global`eps);
  If[recurrenceBackend === "Cpp",
    nativeEpsilonPlan = ft2BuildNativeEpsilonPlan[
      ftData, epsOrder, levelEpsilonHalos, normalizeFT];
    If[FailureQ[nativeEpsilonPlan],
      Print["FTLADDER NATIVE EPSILON PLAN FAIL ", nativeEpsilonPlan];
      Return[$Failed]];
    Print["FTLADDER NATIVE EPSILON PLAN identity=",
      nativeEpsilonPlan["Identity"],
      " losses=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "IntrinsicLoss"],
      " required=", Lookup[nativeEpsilonPlan["Record", "Levels"],
        "RequiredRawTop"]]];
  If[resumeLadderFile =!= "",
    resumeCheckpoint = loadLadderCheckpoint[ExpandFileName[resumeLadderFile],
      name, ftData, prepKey, nativeEpsilonPlan];
    If[resumeCheckpoint === $Failed, Return[$Failed]]];
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
      ftData, boundaryOrder];
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
            ftData, boundaryOrder];
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
    printRows[name, nLevels, ftData["Levels"][nLevels]["Masters"],
      Table[DiffExp2`EpsSeries`ESShift[
        DiffExp2`EpsSeries`ESNew[0, currentBCs[[i]]], -currentPrefactors[[i]]],
        {i, Length[currentBCs]}],
      currentPrefactors]];

  Module[{abortRes},
  abortRes = Catch[
  Do[Module[
    {levelData = ftData["Levels"][level], levelBelow = ftData["Levels"][level - 1],
     var, A, sys, mastersBelow, mastersHere, requests, reductions,
     extraFacs, rawES, rawMin, shift, kmaxAvail, nextReq, needTop,
     trLoCache, trHiCache, chartCache,
     resumeTransport, levelExpansionOrder, needInt, needLo, needHi,
     transportCheckpointFile, saveTransportProgress, completedArms,
     transportSys = None, planLo = None, planHi = None, armReq,
     loPlanCharts, hiPlanCharts, armRounds, armBatchResult,
     armUniqueCharts, armCacheCapacity, levelIBPBatch, rawExtraFacs,
     epsilonBasis, epsilonBasisRecord, nativeEntries = None,
     nativeLedger = None, nativeDispatch = None, downstreamFiniteTop,
     esCMxLevel, configResult, deltaPrescriptions, plannedLevel = None,
     runtimePlanCheck},
    var = levelData["FeynmanParameter"];
    If[recurrenceBackend === "Cpp",
      plannedLevel = nativeEpsilonPlan["Levels"][level]];
    A = normalizeFT[levelData["DiffMatrix"]];
    Print["LEVEL ", level, " var=", var, " d=", Length[A]];
    mastersHere = levelData["Masters"];
    mastersBelow = levelBelow["Masters"];
    resumeTransport = AssociationQ[resumeCheckpoint] &&
      resumeCheckpoint["Kind"] === "Transport" &&
      level === resumeCheckpoint["Level"];
    epsilonBasis = ft2NormalizeEpsilonBasis[
      A, currentBCs, currentPrefactors, Global`eps,
      If[AssociationQ[plannedLevel], plannedLevel["Gauge"], Automatic]];
    If[FailureQ[epsilonBasis],
      Print["FTLADDER EPS BASIS FAIL level=", level, " ", epsilonBasis];
      Throw[$Failed, "FT2Abort"]];
    A = epsilonBasis["Matrix"];
    currentBCs = epsilonBasis["BoundaryValues"];
    currentPrefactors = epsilonBasis["BoundaryPrefactors"];
    epsilonBasisRecord = epsilonBasis["CheckpointRecord"];
    If[resumeTransport,
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
    levelExpansionOrder = If[resumeTransport,
      resumeCheckpoint["SourceExpansionOrder"], expansionOrder];
    (* CoefficientVectors are part of the production transport contract.
       Rebuild/revalidate the one batched FIRE payload even when the Wolfram
       backend resumes a legacy arm snapshot; native Transport snapshots were
       rejected at the loader and can never be reinterpreted as retained
       state. *)
    levelIBPBatch = If[AssociationQ[plannedLevel],
      plannedLevel["Batch"],
      FeynmanTrick`LevelReduction`PrepareLevelIBPBatch[ftData, level]];
    If[levelIBPBatch === $Failed,
      Print["FIRE FAIL"]; Return[$Failed, Module]];
    requests = levelIBPBatch["BoundaryRequests"];
    reductions = levelIBPBatch["Reductions"];
    rawExtraFacs =
      FeynmanTrick`LevelReduction`CollectLevelIBPSingularFactors[
        ftData, level, levelIBPBatch];
    If[rawExtraFacs === $Failed,
      Print["FIRE BATCH FAIL"]; Return[$Failed, Module]];
    extraFacs = normalizeFT[rawExtraFacs];
    If[resumeTransport,
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
      DiffExp2`Solve`ClearSolveCaches[],
      sys = catch2[DiffExp2`API`LoadSystem[
        <|"Matrix" -> A, "Variable" -> var|>]];
      If[FailureQ[sys], Print["LOAD FAIL ", sys]; Return[$Failed, Module]]];
    If[recurrenceBackend === "Cpp",
      nativeEntries = ft2PrepareBoundaryEntries[level, levelIBPBatch,
        currentPrefactors, var, Global`eps, normalizeFT];
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
      nativeLedger = ft2NativeEpsilonLedger[
        nativeEntries, currentBCs, downstreamFiniteTop];
      If[FailureQ[nativeLedger],
        Print["FTLADDER NATIVE EPSILON FAIL level=", level, " ",
          nativeLedger];
        Throw[$Failed, "FT2Abort"]];
      esCMxLevel = nativeLedger["TargetCompleteMax"],
      esCMxLevel = requestedEpsilonOrder[level]];
    (* configure DiffExp2 for this level *)
    deltaPrescriptions = levelDeltaPrescriptions[var, sys, extraFacs];
    configResult = catch2[DiffExp2`Config`LoadConfiguration[{
      "WorkingPrecision" -> wp, "ExpansionOrder" -> levelExpansionOrder,
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
      nativeDispatch = ft2RunNativeBoundaryDispatch[
        sys, currentBCs, nativeEntries, nativeLedger, var, anchor,
        extraFacs, deltaPrescriptions, cppArmThreadBudget, wp,
        nativeEpsilonExecution["Identity"]];
      If[FailureQ[nativeDispatch],
        Print["FTLADDER NATIVE BATCH FAIL level=", level, " ",
          nativeDispatch];
        Throw[$Failed, "FT2Abort"]];
      rawES = nativeDispatch["Values"];
      Print["FTLADDER NATIVE BATCH level=", level,
        " requests=", Length[nativeEntries],
        " batchCalls=", nativeDispatch["NativeBatchCalls"],
        " armMarches=", nativeDispatch["NativeMarches"],
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
        "DivisionOrder" -> divisionOrder,
        "RadiusOfConvergence" -> radiusOfConvergence,
        "ValueTransportMode" -> Environment["DE2_VALUE_TRANSPORT"],
        "RecurrenceBackend" -> recurrenceBackend,
        "SingularMatchPrecondition" -> singularMatchPrecondition,
        "DeltaPrescriptionSign" -> deltaPrescriptionSign,
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
    rawES = DiffExp2`EpsSeries`ESTrim /@ rawES;
    printRows[name, level - 1, mastersBelow, rawES,
      ConstantArray[0, Length[mastersBelow]]];
    (* shift to finite for the next level's transport *)
    rawMin = Min[esMn /@ rawES];
    shift = Max[0, -rawMin];
    kmaxAvail = esCMx /@ rawES;
    If[level > 1,
      If[recurrenceBackend === "Cpp",
        nextReq = nativeEpsilonPlan["Levels"][level - 1]
          ["RequiredRawTop"];
        needTop = Min[kmaxAvail],
        nextReq = requestedEpsilonOrder[level - 1];
        needTop = nextReq - shift];
      If[(recurrenceBackend === "Cpp" && needTop < nextReq) ||
          AnyTrue[kmaxAvail, # < needTop &],
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
          "MastersHere" -> mastersBelow, "Anchor" -> anchor,
          "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
          "RecurrenceBackend" -> recurrenceBackend,
          "DeltaPrescriptionSign" -> deltaPrescriptionSign,
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
            nativeDispatch["CheckpointRecord"], None],
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
