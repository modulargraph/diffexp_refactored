(* Regression tests for safe, invocation-local FIRE reuse and the FIRE
   process wrapper.  Fake runners are used throughout, so no FIRE executable,
   database, or license seat is used. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

Module[{dir, topology, unverifiedTopology, otherTopology, verifiedTopology,
        relocatedTopology,
        changedSetupTopology, oldCache, oldReductionCacheOption,
        oldDimensionVariable, cachedHitResult, callsBeforeCachedHit,
        blockedResult, callsBeforeBlockedMiss,
        runCalls = 0, fireRequests = {}, results, first, second, key1, key2,
        unverifiedKey, otherKey, verifiedKey, relocatedKey, changedSetupKey},
  dir = FileNameJoin[{$TemporaryDirectory,
    "ft_fire_inmemory_cache_" <> ToString[$ProcessID]}];
  If[DirectoryQ[dir], DeleteDirectory[dir, DeleteContents -> True]];
  CreateDirectory[dir, CreateIntermediateDirectories -> True];

  topology = FeynmanTrick`FIREInterface`DefineTopology[
    "cache_stub", {Global`l1}, {},
    {1 - Global`l1^2, 2 - Global`l1^2}, {}];
  topology["WorkDirectory"] = dir;
  topology["ProblemNumber"] = 17;
  topology["StartFileReady"] = True;
  Export[FileNameJoin[{dir, "cache_stub.start"}], "fake start", "Text"];
  topology["SetupFingerprintRecord"] = Join[
    <|
      "Schema" -> "FeynmanTrick.FIRESetup/v1",
      "VerifiedStartFile" -> True,
      "StartFileSHA256" ->
        FeynmanTrick`FIREInterface`Private`fileSHA256[
          FileNameJoin[{dir, "cache_stub.start"}]],
      "AutoDetectRestrictions" -> False
    |>,
    FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[]
  ];

  oldCache = FeynmanTrick`FIREInterface`Private`$ReductionCache;
  oldReductionCacheOption = Lookup[
    FeynmanTrick`Private`$FTConfig, "ReductionCache", True];
  FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
  FeynmanTrick`SetFTOption["ReductionCache", True];

  results = Block[{
      FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
      FeynmanTrick`FIREInterface`Private`runFIRE6,
      FIRE`Tables2Rules, FIRE`Tables2Masters},
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
    FeynmanTrick`FIREInterface`Private`runFIRE6[
        _String, runDir_String, _String] := Module[{requestFile, tableFile},
      runCalls++;
      requestFile = FileNameJoin[{runDir, "cache_stub_reduce.m"}];
      AppendTo[fireRequests, Get[requestFile]];
      tableFile = FileNameJoin[{runDir, "cache_stub_reduce.tables"}];
      Export[tableFile, "fake FIRE table", "Text"];
      0];
    FIRE`Tables2Rules[_String] := {
      Global`G[17, {2, 0}] -> 3 Global`G[17, {1, 0}]
    };
    FIRE`Tables2Masters[_String] := {{17, {1, 0}}};

    {
      FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
        topology, {{2}, {2, 0}, {2}}],
      FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
        topology, {{1}, {1, 0}}]
    }
  ];

  {first, second} = results;
  assert["normalized duplicate requests are sent to FIRE once",
    fireRequests === {{{17, {2, 0}}}}];
  assert["reported master identity avoids a second FIRE invocation",
    runCalls === 1];
  assert["deduplicated input reduction is expanded to original keys",
    first["Reductions"][{2}] === 3 Global`G[1, {1, 0}] &&
    first["Reductions"][{2, 0}] === 3 Global`G[1, {1, 0}]];
  assert["cached master identity preserves both public index forms",
    second["Reductions"][{1}] === Global`G[1, {1, 0}] &&
    second["Reductions"][{1, 0}] === Global`G[1, {1, 0}] &&
    MatchQ[second["TablesFile"], Missing["Cached"]]];

  key1 = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    topology, {1, 0}];
  key2 = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    topology, {1, 0}];
  unverifiedTopology = KeyDrop[topology, "SetupFingerprintRecord"];
  unverifiedKey = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    unverifiedTopology, {1, 0}];
  otherTopology = unverifiedTopology;
  otherTopology["ProblemNumber"] = 18;
  otherKey = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    otherTopology, {1, 0}];
  assert["reduction cache key is deterministic for one prepared topology",
    key1 === key2];
  assert["unverified reduction cache key is scoped to the FIRE problem",
    unverifiedKey =!= otherKey];

  verifiedTopology = topology;
  verifiedKey = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    verifiedTopology, {1, 0}];
  relocatedTopology = verifiedTopology;
  relocatedTopology["WorkDirectory"] = dir <> "_relocated";
  relocatedTopology["ProblemNumber"] = 99;
  relocatedKey = FeynmanTrick`FIREInterface`Private`reductionCacheKey[
    relocatedTopology, {1, 0}];
  changedSetupTopology = relocatedTopology;
  changedSetupTopology["SetupFingerprintRecord"]["Restrictions"] =
    {{-1, -1}, {1, 0}};
  changedSetupKey =
    FeynmanTrick`FIREInterface`Private`reductionCacheKey[
      changedSetupTopology, {1, 0}];
  assert["verified setup content survives transient relocation",
    verifiedKey === relocatedKey];
  assert["verified setup restriction change invalidates reduction cache",
    relocatedKey =!= changedSetupKey];

  FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
  FeynmanTrick`FIREInterface`Private`cacheReduction[
    topology, {9, 0}, 1/(Global`d - 4), {{1, 0}}];
  oldDimensionVariable =
    FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  FeynmanTrick`SetFTOption["DimensionVariable", Global`dd];
  callsBeforeCachedHit = runCalls;
  cachedHitResult = Block[{
      FeynmanTrick`FIREInterface`Private`runFIRE6},
    FeynmanTrick`FIREInterface`Private`runFIRE6[___] := (runCalls++; 0);
    FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
      topology, {{9, 0}}]
  ];
  FeynmanTrick`SetFTOption["DimensionVariable", oldDimensionVariable];
  assert["cached 1/(d-4) rejects changed dimension before interpretation",
    cachedHitResult === $Failed && runCalls === callsBeforeCachedHit];

  FeynmanTrick`FIREInterface`Private`$ReductionCache = <||>;
  oldDimensionVariable =
    FeynmanTrick`Private`$FTConfig["DimensionVariable"];
  FeynmanTrick`SetFTOption["DimensionVariable", Global`differentD];
  callsBeforeBlockedMiss = runCalls;
  blockedResult = Block[{
      FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
      FeynmanTrick`FIREInterface`Private`runFIRE6},
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := Null;
    FeynmanTrick`FIREInterface`Private`runFIRE6[___] := (runCalls++; 0);
    FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
      topology, {{2, 0}}]
  ];
  FeynmanTrick`SetFTOption["DimensionVariable", oldDimensionVariable];
  assert["cache miss refuses changed FIRE runtime identity before launch",
    blockedResult === $Failed && runCalls === callsBeforeBlockedMiss];

  FeynmanTrick`FIREInterface`Private`$ReductionCache = oldCache;
  FeynmanTrick`SetFTOption["ReductionCache", oldReductionCacheOption];
  Quiet[DeleteDirectory[dir, DeleteContents -> True]];
];

Module[{dir, quickBin, slowBin, pidFile, oldTimeout, chmodResults,
        quickExit, quickSeconds, quickLog, slowExit, slowSeconds, slowLog,
        pids = {}, processGoneQ, orphanFree},
  dir = FileNameJoin[{$TemporaryDirectory,
    "ft_fire_runner_polling_" <> ToString[$ProcessID]}];
  If[DirectoryQ[dir], DeleteDirectory[dir, DeleteContents -> True]];
  CreateDirectory[dir, CreateIntermediateDirectories -> True];
  quickBin = FileNameJoin[{dir, "fake_FIRE6_quick"}];
  slowBin = FileNameJoin[{dir, "fake_FIRE6_slow"}];
  pidFile = FileNameJoin[{dir, "slow_pids.txt"}];

  Export[quickBin,
    "#!/bin/sh\nprintf 'quick fake FIRE completed\\n'\nexit 0\n", "Text"];
  Export[slowBin, StringJoin[
    "#!/bin/sh\n",
    "sleep 10 &\n",
    "child=$!\n",
    "printf '%s %s\\n' \"$$\" \"$child\" > \"", pidFile, "\"\n",
    "wait \"$child\"\n"
  ], "Text"];
  chmodResults = Lookup[
    RunProcess[{"/bin/chmod", "+x", #}], "ExitCode", -1] & /@
      {quickBin, slowBin};
  assert["fake FIRE runners are executable", chmodResults === {0, 0}];

  oldTimeout = Lookup[
    FeynmanTrick`Private`$FTConfig, "FIRETimeoutSeconds", 600];
  FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 5];
  quickSeconds = First@AbsoluteTiming[
    quickExit =
      FeynmanTrick`FIREInterface`Private`runFIRE6Once[
        quickBin, dir, "quick"]];
  quickLog = FeynmanTrick`FIREInterface`Private`safeReadString[
    FileNameJoin[{dir, "fire_stdout.log"}]];
  assert["short FIRE process does not pay the former two-second poll",
    quickExit === 0 && quickSeconds < 1.5];
  assert["short FIRE process output is preserved",
    StringContainsQ[quickLog, "quick fake FIRE completed"]];

  FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1];
  slowSeconds = First@AbsoluteTiming[
    slowExit =
      FeynmanTrick`FIREInterface`Private`runFIRE6Once[
        slowBin, dir, "slow"]];
  slowLog = FeynmanTrick`FIREInterface`Private`safeReadString[
    FileNameJoin[{dir, "fire_stdout.log"}]];
  If[FileExistsQ[pidFile],
    pids = Quiet[ToExpression /@ StringSplit[
      StringTrim[Import[pidFile, "Text"]]]]];
  processGoneQ[pid_] := Lookup[
    RunProcess[{"/bin/kill", "-0", ToString[pid]}],
    "ExitCode", 0] =!= 0;
  (* A terminated grandchild can remain visible as an init-owned zombie for
     a scheduler tick after the monitored parent has been explicitly reaped.
     Allow bounded OS reaping time, but a genuinely leaked ten-second worker
     still fails this assertion. *)
  orphanFree = False;
  Do[
    orphanFree = Length[pids] === 2 && AllTrue[pids, processGoneQ];
    If[TrueQ[orphanFree], Break[]];
    Pause[0.05],
    {40}
  ];
  assert["FIRE timeout keeps exit code and exact diagnostic",
    slowExit === 124 &&
      StringContainsQ[slowLog, "FIRE6 timeout after 1s"] &&
      slowSeconds < 4.5];
  assert["FIRE timeout cleans parent and child processes", orphanFree];
  assert["FIRE runner removes isolated attempt directories",
    FileNames["run_*", dir] === {}];

  FeynmanTrick`SetFTOption["FIRETimeoutSeconds", oldTimeout];
  Quiet[DeleteDirectory[dir, DeleteContents -> True]];
];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
