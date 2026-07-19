(* A consumed terminal CASE-P Rational-shadow basis is private proof state:
   the live exact witness must remain forbidden before matching, while the
   completed terminal factorization must checkpoint and restore for later
   observable contraction. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 80,
  "ChopPrecision" -> 40, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0, "DivisionOrder" -> 3,
  "Variables" -> {}, "Verbosity" -> 0}];

x = Global`x; eps = Global`eps;
(* A regular two-block DAG deterministically retains the terminal exact-right
   factorization.  Forcing its basis through the Rational-shadow importer
   exercises the same live-only witness ownership as CASE-P, without making
   the lifecycle test depend on a particular singular-frame choice.  The
   singular-SCC tests cover the CASE-P schedule itself. *)
system = DiffExp2`LoadSystem[<|
  "Matrix" -> {{(1/2 + eps/3)/x, 1}, {0, 0}},
  "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 1/4}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 3/4}];
(* The physical target ends at eps^0, but this CASE-P frame has weights down
   to eps^-2.  Retain the structured retry's three private source orders so
   the terminal residual has a common complete window through the target. *)
zero = DiffExp2`EpsSeries`ESNew[-2, {0, 0, 0, 0, 0, 0}];
one = DiffExp2`EpsSeries`ESNew[-2, {0, 0, 1, 0, 0, 0}];
observable = <|"Operation" -> "limitUpper",
  "Identity" -> "casep-terminal-limit",
  "CheckpointIdentity" -> "casep-terminal-limit:checkpoint",
  "CoefficientVector" -> {1, 0},
  "Epsilon" -> <|"Min" -> -2, "Max" -> 0,
    "RequiredCompleteMax" -> 0|>|>;

forcedShadowCalls = 0;
Clear[forcedReceivingBasis];
forcedReceivingBasis[receiver_, request_, threads_, rest___] :=
  Module[{target},
    forcedShadowCalls++;
    target = DiffExp2`Solve`PrepareNativeSCCComposite[
      receiver, request];
    DiffExp2`NativeTransport`Private`nativeRationalShadowBasis[
      receiver, request, threads, target]];

atlas = Block[{
    DiffExp2`NativeTransport`Private`nativeReceivingBasis =
      forcedReceivingBasis},
  catchDE2[
    DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
      system, {one, zero}, lowerPlan, upperPlan,
      "Threads" -> 1, "Integrands" -> {{{1, 0}}, x},
      "TargetCompleteMax" -> 0, "DeferReceivingBases" -> False]]];
(* The eager plan is now bound to the SCC owners selected by the forced
   shadow bases.  Drop those first basis handles and stream fresh ones during
   the consuming run, exactly as production does, so no unrelated public
   witness can account for the checkpoint result. *)
eagerColumns = If[AssociationQ[atlas], Flatten[
  Lookup[Join[Rest[atlas["Lower", "Bases"]],
      Rest[atlas["Upper", "Bases"]]], "Columns", {}], 1], {}];
Scan[Quiet[DiffExp2`CppBackend`ReleasePersistentLocal[#]] &,
  eagerColumns];
If[AssociationQ[atlas],
  atlas = Join[atlas, <|
    "DeferredReceivingBases" -> True,
    "Lower" -> Join[atlas["Lower"], <|
      "Bases" -> ConstantArray[None,
        Length[atlas["Lower", "ChartSystems"]]]|>],
    "Upper" -> Join[atlas["Upper"], <|
      "Bases" -> ConstantArray[None,
        Length[atlas["Upper", "ChartSystems"]]]|>]|>]];
forcedShadowCalls = 0;

run = If[FailureQ[atlas], atlas,
  Block[{DiffExp2`NativeTransport`Private`nativeReceivingBasis =
      forcedReceivingBasis},
    catchDE2[
      DiffExp2`NativeTransport`RunNativeTransportObservableBatchOwned[
        atlas, {observable}, x, "MaxRefinementSteps" -> 1]]]];
stateStats = If[AssociationQ[run],
  Map[
    DiffExp2`CppBackend`PersistentTransportArmStatistics,
    run["States"]], run];
assert["forced CASE-P shadow reaches a consumed terminal match",
  AssociationQ[run] && forcedShadowCalls === 2 &&
    AllTrue[Values[stateStats],
      AssociationQ[#] &&
        Lookup[#, "status", "error"] === "ok" &] &&
    AllTrue[Values[stateStats],
      TrueQ[Lookup[#, "terminal_factorized_match", False]] &]];

checkpointPath = FileNameJoin[{$TemporaryDirectory,
  "diffexp2-consumed-casep-terminal-" <>
    ToString[$ProcessID] <> ".de2cp"}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
manifest = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`SaveNativeTransportObservableBatchCheckpoint[
    run, checkpointPath, "consumed-casep-terminal-roundtrip"]], run];
assert["consumed CASE-P terminal factorization checkpoints atomically",
  AssociationQ[manifest] &&
    Lookup[manifest, "Schema", None] ===
      "DiffExp2.NativeTransportObservableCheckpoint/v2" &&
    FileExistsQ[checkpointPath]];

liveExport = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    run, 50]], run];
restored = If[AssociationQ[manifest], catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    manifest]], manifest];
restoredExport = If[AssociationQ[restored], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    restored, 50]], restored];
assert["restored CASE-P terminal result is byte-identical and does not remarch",
  AssociationQ[liveExport] && AssociationQ[restoredExport] &&
    Lookup[Lookup[liveExport, "ExportedResults", {}], "Value"] ===
      Lookup[Lookup[restoredExport, "ExportedResults", {}], "Value"] &&
    Lookup[restoredExport, "NativeMarches", -1] === 0];

restoredRelease = If[AssociationQ[restored],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    restored], restored];
liveRelease = If[AssociationQ[run],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[
    run], run];
assert["live and restored CASE-P terminal owners release cleanly",
  AssociationQ[restoredRelease] && AssociationQ[liveRelease] &&
    Lookup[restoredRelease, "Failures", {None}] === {} &&
    Lookup[liveRelease, "Failures", {None}] === {}];

If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0,
  Print["Diagnostics: ", InputForm[<|
    "ForcedShadowCalls" -> forcedShadowCalls,
    "RunType" -> If[AssociationQ[run], Lookup[run, "Type", None], run],
    "StateStats" -> If[AssociationQ[stateStats],
      Map[
        KeyTake[#, {"status", "arm", "matches",
          "terminal_factorized_match"}] &, stateStats], stateStats],
    "Manifest" -> If[AssociationQ[manifest],
      KeyTake[manifest, {"Schema", "Path", "CheckpointIdentity"}],
      manifest],
    "LiveValues" -> If[AssociationQ[liveExport],
      Lookup[Lookup[liveExport, "ExportedResults", {}], "Value"],
      liveExport],
    "RestoredType" -> If[AssociationQ[restored],
      Lookup[restored, "Type", None], restored],
    "RestoredNativeMarches" -> If[AssociationQ[restoredExport],
      Lookup[restoredExport, "NativeMarches", None], restoredExport],
    "RestoredValues" -> If[AssociationQ[restoredExport],
      Lookup[Lookup[restoredExport, "ExportedResults", {}], "Value"],
      restoredExport]|>]];
  Exit[1],
  Exit[0]];
