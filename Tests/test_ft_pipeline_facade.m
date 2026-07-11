(* Static/package-level contract for the DiffExp2 Feynman-trick facade.
   PipelinePlan must not launch wolframscript or FIRE. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

plan = FeynmanTrick`PipelinePlan["bubble",
  "WorkingPrecision" -> 250,
  "ExpansionOrder" -> 40,
  "CppThreads" -> 6,
  "LevelEpsilonHalos" -> {0, 4}];

assert["pipeline plan is an association", AssociationQ[plan]];
assert["root facade and legacy implementation keep separate contexts",
  Context[FeynmanTrick`RunIntegrationPipeline] === "FeynmanTrick`" &&
  Context[FeynmanTrick`DiffExpIntegration`RunIntegrationPipeline] ===
    "FeynmanTrick`DiffExpIntegration`" &&
  Length[DownValues[FeynmanTrick`RunIntegrationPipeline]] > 0 &&
  Length[DownValues[
    FeynmanTrick`DiffExpIntegration`RunIntegrationPipeline]] > 0];
assert["pipeline plan schema",
  plan["Schema"] === "FeynmanTrick.PipelinePlan/v1"];
assert["C++ recurrence is the facade default",
  plan["Settings", "RecurrenceBackend"] === "Cpp" &&
  plan["Environment", "DE2_RECURRENCE_BACKEND"] === "Cpp"];
assert["fast transport settings are explicit",
  plan["Environment", "DE2_VALUE_TRANSPORT"] === "1" &&
  plan["Environment", "FT_CPP_BATCH_ENDPOINT_ARMS"] === "1" &&
  plan["Environment", "DE2_CPP_THREADS"] === "6"];
assert["epsilon halos serialize deterministically",
  plan["Environment", "FT_LEVEL_EPS_HALOS"] === "0,4"];
assert["argv does not use a shell",
  ListQ[plan["Command"]] && plan["Command"][[1]] === "wolframscript" &&
    plan["Command"][[2]] === "-file"];

checkpoint = CreateTemporary[];
resume = FeynmanTrick`PipelinePlan["bubble",
  "ResumeFrom" -> checkpoint];
assert["resume checkpoint is wired into the typed plan",
  resume["ResumeFrom"] === ExpandFileName[checkpoint] &&
  resume["Environment", "FT_RESUME_LADDER_CHECKPOINT"] ===
    ExpandFileName[checkpoint] &&
  resume["CheckpointDirectory"] === DirectoryName[ExpandFileName[checkpoint]]];
DeleteFile[checkpoint];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
