(* Real subprocess smoke test for the DiffExp2 Feynman-trick facade.  The
   child is a tiny fake runner: no FIRE setup or recurrence solve occurs. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

runner = FileNameJoin[{$TemporaryDirectory,
  "ft-pipeline-fake-runner-" <> ToString[$ProcessID] <> ".m"}];
Export[runner, StringRiffle[{
  "Print[\"STEPWISE \", ExportString[<|\"Example\" -> Environment[\"FT_EXAMPLES\"], \"Level\" -> 1|>, \"RawJSON\", \"Compact\" -> True]];",
  "Print[\"FINAL \", ExportString[<|\"Example\" -> Environment[\"FT_EXAMPLES\"], \"Value\" -> 7, \"FIREPath\" -> Environment[\"FT_FIRE_PATH\"]|>, \"RawJSON\", \"Compact\" -> True]];"
}, "\n"], "Text"];

result = FeynmanTrick`RunIntegrationPipeline["bubble",
  "Runner" -> runner,
  "WorkingDirectory" -> repoRoot,
  "FIREPath" -> FileNameJoin[{$TemporaryDirectory, "ft-pipeline-fire"}],
  "PreparedCacheDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-prep"}],
  "CheckpointDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-checkpoints"}]];

planned = FeynmanTrick`PipelinePlan["bubble",
  "Runner" -> runner,
  "WorkingDirectory" -> repoRoot,
  "PreparedCacheDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-plan-process-prep"}],
  "CheckpointDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-plan-process-checkpoints"}]];
plannedResult = FeynmanTrick`RunIntegrationPipeline[planned];

If[!AssociationQ[result] || Lookup[result, "Status", "Failed"] =!= "Succeeded",
  Print["  DIAGNOSTIC: ", InputForm[result]]];

assert["pipeline child returns the typed result schema",
  AssociationQ[result] &&
    result["Schema"] === "FeynmanTrick.PipelineResult/v1" &&
    result["Status"] === "Succeeded" && result["ExitCode"] === 0];
assert["pipeline child receives the planned environment",
  result["Final", "Example"] === "bubble" &&
    result["Final", "Value"] === 7 &&
    result["Final", "FIREPath"] ===
      ExpandFileName[FileNameJoin[{$TemporaryDirectory, "ft-pipeline-fire"}]]];
assert["pipeline result retains parsed stepwise rows",
  Length[result["Stepwise"]] === 1 &&
    result["Stepwise"][[1, "Level"]] === 1];
assert["an existing typed plan can be executed",
  plannedResult["Schema"] === "FeynmanTrick.PipelineResult/v1" &&
    plannedResult["Status"] === "Succeeded" &&
    plannedResult["Final", "Example"] === "bubble"];

If[FileExistsQ[runner], DeleteFile[runner]];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
