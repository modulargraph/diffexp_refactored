(* Process-free tests for the opt-in persistent FIRE database seam.  FIRE is
   stubbed throughout: this test must never start FIRE6 or a Fermat worker. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0;
failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]
];

beginSession = FeynmanTrick`FIREInterface`Private`beginFIREStorageSession;
resetAttempt = FeynmanTrick`FIREInterface`Private`resetFIREStorageAttempt;
publishSession = FeynmanTrick`FIREInterface`Private`publishFIREStorageSession;
abortSession = FeynmanTrick`FIREInterface`Private`abortFIREStorageSession;
safeRead = FeynmanTrick`FIREInterface`Private`safeReadString;
cachedReduction = FeynmanTrick`FIREInterface`Private`cachedReduction;
oldConfig = FeynmanTrick`FTConfiguration[];
root = FileNameJoin[{$TemporaryDirectory,
  "ft-fire-storage-test-" <> ToString[$ProcessID] <> "-" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
work = FileNameJoin[{root, "work"}];
cacheRoot = FileNameJoin[{root, "cache"}];
fireRoot = FileNameJoin[{root, "fake-fire"}];
fireBin = FileNameJoin[{fireRoot, "bin", "FIRE6"}];
CreateDirectory[work, CreateIntermediateDirectories -> True];
CreateDirectory[DirectoryName[fireBin], CreateIntermediateDirectories -> True];
Export[fireBin, "fake FIRE build A\n", "Text"];

topology = <|
  "Name" -> "persistent-test",
  "ProblemNumber" -> 7,
  "NumPropagators" -> 2,
  "WorkDirectory" -> work,
  "Variables" -> {Global`x},
  "StartFileReady" -> True
|>;
startFile = FileNameJoin[{work, topology["Name"] <> ".start"}];
startContent = "fake start content A\n";
Export[startFile, startContent, "Text"];

FeynmanTrick`SetFTOption["FIREPath", fireRoot];
FeynmanTrick`SetFTOption["WorkDirectory", work];
FeynmanTrick`SetFTOption["FIREStorageDirectory", cacheRoot];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;

Print["=== Persistent FIRE storage lifecycle ==="];
FeynmanTrick`SetFTOption["PersistentFIREStorage", False];
assert["persistent storage is opt-in",
  beginSession[topology, fireBin] === None && !DirectoryQ[cacheRoot]];
FeynmanTrick`SetFTOption["PersistentFIREStorage", True];
session1 = beginSession[topology, fireBin];
assert["session starts with a content key", AssociationQ[session1] &&
  StringLength[session1["Key"]] === 64];
assert["first attempt starts from an empty database",
  resetAttempt[session1] &&
  DirectoryQ[session1["AttemptStorage"]]];
marker = FileNameJoin[{session1["AttemptStorage"], "marker.txt"}];
Export[marker, "generation one", "Text"];
assert["successful attempt publishes", publishSession[session1]];
keyDir = session1["KeyDirectory"];
pointer1 = StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]];
assert["CURRENT names an immutable generation",
  StringMatchQ[pointer1, RegularExpression["gen-[0-9a-f]{32}"]] &&
  FileExistsQ[FileNameJoin[{keyDir, "generations", pointer1,
    "storage", "marker.txt"}]]];

session2 = beginSession[topology, fireBin];
assert["next attempt clones the committed generation",
  resetAttempt[session2] &&
  FileExistsQ[FileNameJoin[{session2["AttemptStorage"], "marker.txt"}]]];
abortSession[session2];
assert["abort removes attempt and lock without moving CURRENT",
  !DirectoryQ[session2["AttemptDirectory"]] &&
  !DirectoryQ[session2["LockDirectory"]] &&
  StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]] === pointer1];

session3 = beginSession[topology, fireBin];
contended = beginSession[topology, fireBin];
assert["same-key writers are serialized", contended === $Failed];
abortSession[session3];

Export[startFile, startContent <> "changed\n", "Text"];
changedSession = beginSession[topology, fireBin];
assert["start-file content changes the storage key",
  AssociationQ[changedSession] && changedSession["Key"] =!= session1["Key"]];
abortSession[changedSession];
Export[startFile, startContent, "Text"];
Export[fireBin, "fake FIRE build B\n", "Text"];
changedBuildSession = beginSession[topology, fireBin];
assert["FIRE build content changes the storage key",
  AssociationQ[changedBuildSession] &&
  changedBuildSession["Key"] =!= session1["Key"]];
abortSession[changedBuildSession];
Export[fireBin, "fake FIRE build A\n", "Text"];

Print["=== Stubbed reduction publication and cache ==="];
requested = {2, 1};
master = {1, 1};
pointerBeforeSuccess = StringTrim[
  safeRead[FileNameJoin[{keyDir, "CURRENT"}]]];
configObserved = "";
result = Block[{
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6,
    FIRE`Tables2Rules,
    FIRE`Tables2Masters},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6[
      _String, dir_String, config_String, session_] := Module[{table},
    configObserved = safeRead[FileNameJoin[{dir, config <> ".config"}]];
    If[AssociationQ[session],
      If[!resetAttempt[session], Return[125, Module]];
      Export[FileNameJoin[{session["AttemptStorage"], "second.txt"}],
        "generation two", "Text"]
    ];
    table = FileNameJoin[{dir, topology["Name"] <> "_reduce.tables"}];
    Export[table, "stub table", "Text"];
    0
  ];
  FIRE`Tables2Rules[_String] := {
    Global`G[7, requested] -> 3 Global`G[7, master]};
  FIRE`Tables2Masters[_String] := {{7, master}};
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
    topology, {requested}]
];
pointerAfterSuccess = StringTrim[
  safeRead[FileNameJoin[{keyDir, "CURRENT"}]]];
masterCache = cachedReduction[topology, master];
assert["config opts into FIRE storage and keepall",
  StringContainsQ[configObserved, "#storage           !"] &&
  StringContainsQ[configObserved, "#keepall"]];
assert["parsed reduction is returned",
  AssociationQ[result] &&
  result["Reductions"][requested] === 3 Global`G[1, master]];
assert["successful parse advances CURRENT",
  pointerAfterSuccess =!= pointerBeforeSuccess &&
  FileExistsQ[FileNameJoin[{keyDir, "generations", pointerAfterSuccess,
    "storage", "second.txt"}]]];
assert["returned master identity is cached",
  AssociationQ[masterCache] &&
  masterCache["Reduction"] === Global`G[1, master]];
assert["publication removes lock and attempt",
  !DirectoryQ[FileNameJoin[{keyDir, "LOCK"}]] &&
  FileNames["attempt-*", FileNameJoin[{keyDir, "attempts"}]] === {}];

retrySession = beginSession[topology, fireBin];
retryAttempt = 0;
retrySawCleanBase = False;
retryBadFile = FileNameJoin[{retrySession["AttemptStorage"], "partial.txt"}];
retryExit = Block[{
    FeynmanTrick`FIREInterface`Private`cleanupFIREProcesses,
    FeynmanTrick`FIREInterface`Private`cleanupFIREWorkDir,
    FeynmanTrick`FIREInterface`Private`runFIRE6Once},
  FeynmanTrick`FIREInterface`Private`cleanupFIREProcesses[___] := Null;
  FeynmanTrick`FIREInterface`Private`cleanupFIREWorkDir[___] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6Once[___] := Module[{},
    retryAttempt++;
    If[retryAttempt === 1,
      Export[retryBadFile, "partial failed mutation", "Text"];
      1,
      retrySawCleanBase = !FileExistsQ[retryBadFile] &&
        FileExistsQ[FileNameJoin[{retrySession["AttemptStorage"],
          "second.txt"}]];
      0
    ]
  ];
  FeynmanTrick`FIREInterface`Private`runFIRE6[
    fireBin, work, "retry-stub", retrySession]
];
assert["retry discards partial mutation and reclones CURRENT",
  retryExit === 0 && retryAttempt === 2 && retrySawCleanBase];
abortSession[retrySession];

Print["=== Failure rollback ==="];
pointerBeforeFailure = pointerAfterSuccess;
parseFailure = Block[{
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6,
    FIRE`Tables2Rules,
    FIRE`Tables2Masters},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6[
      _String, dir_String, _String, session_] := Module[{},
    If[!resetAttempt[session], Return[125, Module]];
    Export[FileNameJoin[{session["AttemptStorage"], "uncommitted.txt"}],
      "must be discarded", "Text"];
    Export[FileNameJoin[{dir, topology["Name"] <> "_reduce.tables"}],
      "invalid stub table", "Text"];
    0
  ];
  FIRE`Tables2Rules[_String] := $Failed;
  FIRE`Tables2Masters[_String] := {{7, master}};
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[topology, {{3, 1}}]
];
assert["parse failure is loud", parseFailure === $Failed];
assert["parse failure preserves CURRENT and cleans attempt",
  StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]] ===
    pointerBeforeFailure &&
  !DirectoryQ[FileNameJoin[{keyDir, "LOCK"}]] &&
  FileNames["attempt-*", FileNameJoin[{keyDir, "attempts"}]] === {}];

compatibilityFailure = Block[{
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6,
    FIRE`Tables2Rules,
    FIRE`Tables2Masters},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6[
      _String, dir_String, _String, session_] := Module[{},
    If[!resetAttempt[session], Return[125, Module]];
    Export[FileNameJoin[{session["AttemptStorage"], "wrong-key.txt"}],
      "must be discarded", "Text"];
    Export[startFile, startContent <> "changed during run\n", "Text"];
    Export[FileNameJoin[{dir, topology["Name"] <> "_reduce.tables"}],
      "valid stub table", "Text"];
    0
  ];
  FIRE`Tables2Rules[_String] := {
    Global`G[7, {5, 1}] -> Global`G[7, master]};
  FIRE`Tables2Masters[_String] := {{7, master}};
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[topology, {{5, 1}}]
];
Export[startFile, startContent, "Text"];
assert["compatibility mutation prevents publication",
  compatibilityFailure === $Failed &&
  StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]] ===
    pointerBeforeFailure &&
  !DirectoryQ[FileNameJoin[{keyDir, "LOCK"}]] &&
  FileNames["attempt-*", FileNameJoin[{keyDir, "attempts"}]] === {}];

runFailure = Block[{
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6[
      _String, _String, _String, session_] := Module[{},
    resetAttempt[session];
    Export[FileNameJoin[{session["AttemptStorage"], "failed-run.txt"}],
      "must be discarded", "Text"];
    1
  ];
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[topology, {{4, 1}}]
];
assert["FIRE exit failure is loud", runFailure === $Failed];
assert["FIRE exit failure preserves CURRENT and cleans attempt",
  StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]] ===
    pointerBeforeFailure &&
  !DirectoryQ[FileNameJoin[{keyDir, "LOCK"}]] &&
  FileNames["attempt-*", FileNameJoin[{keyDir, "attempts"}]] === {}];

Export[FileNameJoin[{work, topology["Name"] <> "_reduce.tables"}],
  "stale table from an earlier success", "Text"];
noOutputFailure = Block[{
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
    FeynmanTrick`FIREInterface`Private`runFIRE6},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
  FeynmanTrick`FIREInterface`Private`runFIRE6[
      _String, _String, _String, session_] := Module[{},
    resetAttempt[session];
    Export[FileNameJoin[{session["AttemptStorage"], "no-output.txt"}],
      "must be discarded", "Text"];
    0
  ];
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[topology, {{6, 1}}]
];
assert["zero exit without a fresh table cannot publish stale output",
  noOutputFailure === $Failed &&
  StringTrim[safeRead[FileNameJoin[{keyDir, "CURRENT"}]]] ===
    pointerBeforeFailure &&
  !DirectoryQ[FileNameJoin[{keyDir, "LOCK"}]] &&
  FileNames["attempt-*", FileNameJoin[{keyDir, "attempts"}]] === {}];

FeynmanTrick`Private`$FTConfig = oldConfig;
FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
If[DirectoryQ[root], DeleteDirectory[root, DeleteContents -> True]];

Print["Results: ", passed, " passed, ", failed, " failed."];
If[failed > 0, Exit[1]];
