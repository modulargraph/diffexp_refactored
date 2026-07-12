(* Regression for checkpoint manifests with one empty public result class.
   Production FT integration-only levels have lines and no endpoints. *)

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
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0, "DivisionOrder" -> 2,
  "Variables" -> {}, "Verbosity" -> 0}];

x = Global`x;
system = DiffExp2`LoadSystem[<|"Matrix" -> {{0}}, "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {0, -1/4}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {0, 1/4}];
one = DiffExp2`EpsSeries`ESNew[0, {1}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one}, lowerPlan, upperPlan, "Threads" -> 1,
    "Integrands" -> {{{1}}, x}, "TargetCompleteMax" -> 0]];
observable = <|"Operation" -> "integrate", "Identity" -> "line-only",
  "CheckpointIdentity" -> "line-only:checkpoint",
  "CoefficientVector" -> {1},
  "Epsilon" -> <|"Min" -> 0, "Max" -> 0,
    "RequiredCompleteMax" -> 0|>, "TailPolicy" -> "stored"|>;
run = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeTransportObservableBatch[
    atlas, {observable}, x, "MaxRefinementSteps" -> 1]]];

checkpointPath = FileNameJoin[{$TemporaryDirectory,
  "de2-native-single-visibility-" <> ToString[$ProcessID] <> ".checkpoint"}];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
manifest = If[AssociationQ[run], catchDE2[
  DiffExp2`NativeTransport`SaveNativeTransportObservableBatchCheckpoint[
    run, checkpointPath, "single-visibility-roundtrip"]], run];
assert["integration-only manifest has one line and no endpoint",
  AssociationQ[manifest] &&
    Lookup[Lookup[manifest, "Results", {}], "Kind", {}] === {"line"}];

restored = If[AssociationQ[manifest], catchDE2[
  DiffExp2`NativeTransport`RestoreNativeTransportObservableBatchCheckpoint[
    manifest]], manifest];
exported = If[AssociationQ[restored], catchDE2[
  DiffExp2`NativeTransport`ExportNativeTransportObservableBatch[
    restored, 40]], restored];
value = If[AssociationQ[exported] &&
    Length[Lookup[exported, "ExportedResults", {}]] === 1,
  First[exported["ExportedResults"]]["Value"], None];
assert["integration-only visibility restores and exports without replay",
  AssociationQ[exported] &&
    TrueQ[Lookup[exported, "RestoredCheckpoint", False]] &&
    Lookup[exported, "RestoredNativeMarches", -1] === 2 &&
    Lookup[exported, "NativeMarches", -1] === 0 &&
    DiffExp2`EpsSeries`ESQ[value] &&
    TrueQ[Abs[N[DiffExp2`EpsSeries`ESCoefficient[value, 0] - 1/2,
      25]] < 10^-20]];

If[AssociationQ[exported],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[exported]];
If[AssociationQ[run],
  DiffExp2`NativeTransport`ReleaseNativeTransportObservableBatch[run]];
If[FileExistsQ[checkpointPath], DeleteFile[checkpointPath]];
DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
