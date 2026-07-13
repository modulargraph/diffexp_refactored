(* Focused FIRE 7 lock and completed-cache publication contracts.

   This suite does not invoke FIRE.  It uses short-lived worker kernels to
   exercise real interprocess lock contention, and a fake MPI launcher to
   distinguish a trusted completed cache from an unverified final artifact. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];

Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick", "FIRE7Runner.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

makeExecutable[path_String] :=
  Lookup[RunProcess[{"/bin/chmod", "+x", path}], "ExitCode", -1] === 0;

waitForFile[path_String, seconds_:10] := Module[{deadline = AbsoluteTime[] + seconds},
  While[AbsoluteTime[] < deadline,
    If[FileExistsQ[path], Return[True, Module]];
    Pause[.05]];
  False
];

waitForProcess[process_ProcessObject, seconds_:15] := Module[
  {deadline = AbsoluteTime[] + seconds},
  While[AbsoluteTime[] < deadline,
    If[ProcessStatus[process] =!= "Running", Return[True, Module]];
    Pause[.05]];
  False
];

stopProcess[process_ProcessObject] := If[ProcessStatus[process] === "Running",
  Quiet[KillProcess[process]];
  waitForProcess[process, 5]];

readWXF[path_String] := If[FileExistsQ[path],
  Quiet[Check[Import[path, "WXF"], $Failed]], $Failed];

sha256[path_String] :=
  IntegerString[FileHash[path, "SHA256"], 16, 64];

tempRoot = FileNameJoin[{$TemporaryDirectory,
  "ft_fire7_lock_cache_" <> ToString[$ProcessID] <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
If[DirectoryQ[tempRoot], DeleteDirectory[tempRoot, DeleteContents -> True]];
CreateDirectory[tempRoot, CreateIntermediateDirectories -> True];

runnerPath = FileNameJoin[{repoRoot, "FeynmanTrick", "FIRE7Runner.m"}];

(* ---------------------------------------------------------------------- *)
(* A creator may have made lock/ but not yet atomically published owner.wxf.
   Missing or unreadable metadata in such a fresh lock is indeterminate, not
   proof that the lock is stale.                                           *)
(* ---------------------------------------------------------------------- *)

freshOwnerMetadataCase[tag_String, ownerKind_String] := Module[
  {jobDir, lockDir, ownerPath, readyPath, outcomePath, workerPath,
   workerText, process, ready, preserved, staleLocks},
  jobDir = FileNameJoin[{tempRoot, "fresh-owner-" <> tag}];
  lockDir = FileNameJoin[{jobDir, "lock"}];
  ownerPath = FileNameJoin[{lockDir, "owner.wxf"}];
  readyPath = FileNameJoin[{jobDir, "worker-ready"}];
  outcomePath = FileNameJoin[{jobDir, "worker-outcome.wxf"}];
  workerPath = FileNameJoin[{jobDir, "worker.m"}];
  CreateDirectory[lockDir, CreateIntermediateDirectories -> True];
  If[ownerKind === "Unreadable",
    Export[ownerPath, "this is deliberately not WXF\n", "Text"]];
  workerText = StringJoin[
    "Quiet[Get[", ToString[runnerPath, InputForm], "]];\n",
    "Export[", ToString[readyPath, InputForm], ", \"ready\\n\", \"Text\"];\n",
    "Pause[.10];\n",
    "result = FeynmanTrick`FIRE7Runner`Private`acquireLock[",
      ToString[jobDir, InputForm], "];\n",
    "Export[", ToString[outcomePath, InputForm], ", result, \"WXF\"];\n",
    "If[StringQ[result], Pause[.25]; ",
      "FeynmanTrick`FIRE7Runner`Private`releaseLock[result]];\n"];
  Export[workerPath, workerText, "Text"];
  process = StartProcess[{"wolframscript", "-file", workerPath}];
  ready = waitForFile[readyPath, 10];
  If[ready, Pause[.45]];
  staleLocks = FileNames["lock.stale-*", jobDir];
  preserved = ready && DirectoryQ[lockDir] && staleLocks === {};
  stopProcess[process];
  preserved
];

assert["a fresh lock with missing owner metadata is not stolen immediately",
  freshOwnerMetadataCase["missing", "Missing"]];
assert["a fresh lock with unreadable owner metadata is not stolen immediately",
  freshOwnerMetadataCase["unreadable", "Unreadable"]];

(* ---------------------------------------------------------------------- *)
(* Same-job contenders serialize: the second kernel waits while the live
   owner holds the lock, then enters the critical section after release.    *)
(* ---------------------------------------------------------------------- *)

contentionDir = FileNameJoin[{tempRoot, "contention"}];
CreateDirectory[contentionDir, CreateIntermediateDirectories -> True];
eventsPath = FileNameJoin[{contentionDir, "events.log"}];
readyPath = FileNameJoin[{contentionDir, "worker-ready"}];
outcomePath = FileNameJoin[{contentionDir, "worker-outcome.wxf"}];
workerPath = FileNameJoin[{contentionDir, "worker.m"}];

appendEvent[path_String, event_String] := Module[{stream = OpenAppend[path]},
  WriteString[stream, event, "\n"];
  Close[stream]];

ownerLock = FeynmanTrick`FIRE7Runner`Private`acquireLock[contentionDir];
appendEvent[eventsPath, "A-enter"];
workerText = StringJoin[
  "Quiet[Get[", ToString[runnerPath, InputForm], "]];\n",
  "appendEvent[path_String, event_String] := Module[{s = OpenAppend[path]}, ",
    "WriteString[s, event, \"\\n\"]; Close[s]];\n",
  "Export[", ToString[readyPath, InputForm], ", \"ready\\n\", \"Text\"];\n",
  "result = FeynmanTrick`FIRE7Runner`Private`withJobLock[",
    ToString[contentionDir, InputForm], ",\n",
  "  appendEvent[", ToString[eventsPath, InputForm], ", \"B-enter\"];\n",
  "  Pause[.10];\n",
  "  appendEvent[", ToString[eventsPath, InputForm], ", \"B-exit\"];\n",
  "  \"B-complete\"];\n",
  "Export[", ToString[outcomePath, InputForm], ", result, \"WXF\"];\n"];
Export[workerPath, workerText, "Text"];
contender = StartProcess[{"wolframscript", "-file", workerPath}];
contenderReady = waitForFile[readyPath, 10];
If[contenderReady, Pause[.40]];
eventsWhileOwned = If[FileExistsQ[eventsPath],
  Select[StringSplit[ReadString[eventsPath], "\n"], # =!= "" &], {}];
appendEvent[eventsPath, "A-exit"];
If[StringQ[ownerLock],
  FeynmanTrick`FIRE7Runner`Private`releaseLock[ownerLock]];
contenderFinished = waitForProcess[contender, 15];
If[!contenderFinished, stopProcess[contender]];
contentionResult = readWXF[outcomePath];
contentionEvents = If[FileExistsQ[eventsPath],
  Select[StringSplit[ReadString[eventsPath], "\n"], # =!= "" &], {}];

assert["same-job contender does not enter while the first owner is live",
  StringQ[ownerLock] && contenderReady && eventsWhileOwned === {"A-enter"}];
assert["same-job contention serializes both critical sections",
  contenderFinished && contentionResult === "B-complete" &&
    contentionEvents === {"A-enter", "A-exit", "B-enter", "B-exit"}];

(* ---------------------------------------------------------------------- *)
(* Cache reuse is a completed-publication contract, not a nonempty-file
   contract.  Only a valid marker for this job and this final SHA may hit.  *)
(* ---------------------------------------------------------------------- *)

configDir = FileNameJoin[{tempRoot, "config"}];
CreateDirectory[configDir, CreateIntermediateDirectories -> True];
Export[FileNameJoin[{configDir, "cache.config"}], StringJoin[
  "#variables x,d\n",
  "#start\n",
  "#problem 1 cache.start\n",
  "#integrals cache.m\n",
  "#output cache.tables\n"], "Text"];
Export[FileNameJoin[{configDir, "cache.m"}], "{{1,{2,1}}}\n", "Text"];
Export[FileNameJoin[{configDir, "cache.start"}], "fake start\n", "Text"];

fireRoot = FileNameJoin[{tempRoot, "fire7"}];
fireBin = FileNameJoin[{fireRoot, "bin"}];
CreateDirectory[fireBin, CreateIntermediateDirectories -> True];
Export[FileNameJoin[{fireRoot, "FIRE7.m"}], "(* fake FIRE7 *)\n", "Text"];
Export[FileNameJoin[{fireRoot, "paths.inc"}],
  "RESULTING_MPRIME_COUNT=16\n", "Text"];
Do[Export[FileNameJoin[{fireBin, name}], "fake binary\n", "Text"],
  {name, {"FIRE7p", "FIRE7mp", "FIRE7_MPI", "reconstruct"}}];

invocationLog = FileNameJoin[{tempRoot, "fake-mpi-invocations.log"}];
fakeMPI = FileNameJoin[{tempRoot, "fake-mpirun"}];
Export[fakeMPI, StringJoin[
  "#!/bin/sh\n",
  "printf 'run\\n' >> ", ToString[invocationLog, InputForm], "\n",
  "mkdir -p artifacts\n",
  "printf 'fresh reconstructed table\\n' > ",
    "artifacts/cache_x_d_0.tables\n",
  "exit 0\n"], "Text"];
assert["fake MPI launcher is executable", makeExecutable[fakeMPI]];

baseSettings = <|
  "MPIExecutable" -> fakeMPI,
  "Workers" -> 1,
  "Calc" -> "flint",
  "UseMultiprime" -> True,
  "PrimeLimit" -> 16,
  "DimensionSeparated" -> False,
  "MultiprimeWidth" -> 16,
  "TimeoutSeconds" -> 10
|>;

identityManifest = FeynmanTrick`FIRE7Runner`ParseConfig[configDir, "cache"];
launcherJobA = FeynmanTrick`FIRE7Runner`Private`modularJob[
  "Reduction", fireRoot, configDir, identityManifest, baseSettings];
fakeMPI2 = FileNameJoin[{tempRoot, "fake-mpirun-copy"}];
Export[fakeMPI2, ReadString[fakeMPI], "Text"];
makeExecutable[fakeMPI2];
launcherJobB = FeynmanTrick`FIRE7Runner`Private`modularJob[
  "Reduction", fireRoot, configDir, identityManifest,
  Join[baseSettings, <|"MPIExecutable" -> fakeMPI2|>]];
Export[fakeMPI2, ReadString[fakeMPI] <> "# changed launcher\n", "Text"];
makeExecutable[fakeMPI2];
launcherJobC = FeynmanTrick`FIRE7Runner`Private`modularJob[
  "Reduction", fireRoot, configDir, identityManifest,
  Join[baseSettings, <|"MPIExecutable" -> fakeMPI2|>]];
assert["modular cache identity binds the MPI launcher path and content",
  AssociationQ[launcherJobA] && AssociationQ[launcherJobB] &&
    AssociationQ[launcherJobC] &&
    DuplicateFreeQ[Lookup[{launcherJobA, launcherJobB, launcherJobC}, "ID"]]];

invocationCount[] := If[FileExistsQ[invocationLog],
  Count[StringSplit[ReadString[invocationLog], "\n"], "run"], 0];

completeRecord[job_Association, final_String] := <|
  "Schema" -> FeynmanTrick`FIRE7Runner`Private`$jobSchema <> "/complete",
  "JobID" -> job["ID"],
  "FinalTable" -> FileNameTake[final],
  "FinalSHA256" -> sha256[final],
  "Completed" -> DateString["ISODateTime"]
|>;

runCacheCase[tag_String, markerKind_String] := Module[
  {settings, manifest, runManifest, job, jobDir, final, markerPath, record,
   canonical, before, result, after, canonicalText},
  settings = Join[baseSettings, <|
    "CacheDirectory" -> FileNameJoin[{tempRoot, "cache-" <> tag}]|>];
  manifest = FeynmanTrick`FIRE7Runner`ParseConfig[configDir, "cache"];
  runManifest = FeynmanTrick`FIRE7Runner`Private`jobManifest[manifest];
  job = FeynmanTrick`FIRE7Runner`Private`modularJob[
    "Reduction", fireRoot, configDir, manifest, settings];
  jobDir = job["Directory"];
  final = FeynmanTrick`FIRE7Runner`ExpectedReconstructedFile[
    jobDir, runManifest];
  CreateDirectory[DirectoryName[final], CreateIntermediateDirectories -> True];
  Export[final, "stale-" <> tag <> "\n", "Text"];
  markerPath = FileNameJoin[{jobDir, "complete.wxf"}];
  record = completeRecord[job, final];
  Switch[markerKind,
    "Valid", Export[markerPath, record, "WXF"],
    "WrongJob", Export[markerPath,
      ReplacePart[record, "JobID" -> StringRepeat["0", 64]], "WXF"],
    "WrongSHA", Export[markerPath,
      ReplacePart[record, "FinalSHA256" -> StringRepeat["0", 64]], "WXF"],
    "Unreadable", Export[markerPath, "not WXF\n", "Text"],
    "Missing", Null];
  canonical = FileNameJoin[{configDir, "cache.tables"}];
  If[FileExistsQ[canonical], DeleteFile[canonical]];
  before = invocationCount[];
  result = FeynmanTrick`FIRE7Runner`RunModular[
    fireRoot, configDir, "cache", settings];
  after = invocationCount[];
  canonicalText = If[FileExistsQ[canonical], ReadString[canonical], ""];
  <|"Result" -> result, "Before" -> before, "After" -> after,
    "CanonicalText" -> canonicalText, "Final" -> final,
    "Marker" -> markerPath, "Job" -> job|>
];

validCase = runCacheCase["valid", "Valid"];
assert["matching JobID and final SHA permit completed-cache reuse",
  AssociationQ[validCase["Result"]] &&
    TrueQ[Lookup[validCase["Result"], "Success", False]] &&
    TrueQ[Lookup[validCase["Result"], "Reused", False]] &&
    validCase["After"] === validCase["Before"] &&
    validCase["CanonicalText"] === "stale-valid\n"];

missingCase = runCacheCase["missing", "Missing"];
assert["a nonempty final without a complete marker is not reused",
  AssociationQ[missingCase["Result"]] &&
    TrueQ[Lookup[missingCase["Result"], "Success", False]] &&
    !TrueQ[Lookup[missingCase["Result"], "Reused", True]] &&
    missingCase["After"] === missingCase["Before"] + 1 &&
    missingCase["CanonicalText"] === "fresh reconstructed table\n"];

wrongJobCase = runCacheCase["wrong-job", "WrongJob"];
assert["a complete marker for another JobID is not reused",
  AssociationQ[wrongJobCase["Result"]] &&
    TrueQ[Lookup[wrongJobCase["Result"], "Success", False]] &&
    !TrueQ[Lookup[wrongJobCase["Result"], "Reused", True]] &&
    wrongJobCase["After"] === wrongJobCase["Before"] + 1 &&
    wrongJobCase["CanonicalText"] === "fresh reconstructed table\n"];

wrongSHACase = runCacheCase["wrong-sha", "WrongSHA"];
assert["a complete marker whose SHA does not match the final is not reused",
  AssociationQ[wrongSHACase["Result"]] &&
    TrueQ[Lookup[wrongSHACase["Result"], "Success", False]] &&
    !TrueQ[Lookup[wrongSHACase["Result"], "Reused", True]] &&
    wrongSHACase["After"] === wrongSHACase["Before"] + 1 &&
    wrongSHACase["CanonicalText"] === "fresh reconstructed table\n"];

unreadableCase = runCacheCase["unreadable", "Unreadable"];
assert["an unreadable complete marker does not make a stale final reusable",
  AssociationQ[unreadableCase["Result"]] &&
    TrueQ[Lookup[unreadableCase["Result"], "Success", False]] &&
    !TrueQ[Lookup[unreadableCase["Result"], "Reused", True]] &&
    unreadableCase["After"] === unreadableCase["Before"] + 1 &&
    unreadableCase["CanonicalText"] === "fresh reconstructed table\n"];

If[DirectoryQ[tempRoot], DeleteDirectory[tempRoot, DeleteContents -> True]];

Print["\nFIRE7 lock/cache contract tests: ", passed, " passed, ",
  failed, " failed."];
If[failed > 0, Exit[1], Exit[0]];
