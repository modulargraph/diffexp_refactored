(* Process-free regression tests for safe, invocation-local FIRE reuse.
   A fake FIRE runner records the exact request file; no FIRE executable,
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

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
