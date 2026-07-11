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
  "Print[\"FINAL \", ExportString[<|\"Example\" -> Environment[\"FT_EXAMPLES\"], \"Value\" -> 7|>, \"RawJSON\", \"Compact\" -> True]];"
}, "\n"], "Text"];

result = FeynmanTrick`RunIntegrationPipeline["bubble",
  "Runner" -> runner,
  "WorkingDirectory" -> repoRoot,
  "PreparedCacheDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-prep"}],
  "CheckpointDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-checkpoints"}]];

If[!AssociationQ[result] || Lookup[result, "Status", "Failed"] =!= "Succeeded",
  Print["  DIAGNOSTIC: ", InputForm[result]]];

assert["pipeline child returns the typed result schema",
  AssociationQ[result] &&
    result["Schema"] === "FeynmanTrick.PipelineResult/v1" &&
    result["Status"] === "Succeeded" && result["ExitCode"] === 0];
assert["pipeline child receives the planned environment",
  result["Final", "Example"] === "bubble" &&
    result["Final", "Value"] === 7];
assert["pipeline result retains parsed stepwise rows",
  Length[result["Stepwise"]] === 1 &&
    result["Stepwise"][[1, "Level"]] === 1];

If[FileExistsQ[runner], DeleteFile[runner]];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
