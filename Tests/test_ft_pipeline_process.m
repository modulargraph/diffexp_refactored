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
  "epsOrder = ToExpression[Environment[\"FT_EPS_ORDER\"]];",
  "coefficients = Table[{power, 7 + power}, {power, 0, epsOrder}];",
  "step = <|\"Example\" -> Environment[\"FT_EXAMPLES\"], \"Level\" -> 1, \"RawMinPower\" -> 0, \"Coefficients\" -> coefficients|>;",
  "final = <|\"Example\" -> Environment[\"FT_EXAMPLES\"], \"Finite\" -> 7, \"RawMinPower\" -> 0, \"FIREPath\" -> Environment[\"FT_FIRE_PATH\"]|>;",
  "If[epsOrder > 0, final = Append[final, \"Coefficients\" -> coefficients]];",
  "Print[\"STEPWISE \", ExportString[step, \"RawJSON\", \"Compact\" -> True]];",
  "Print[\"FINAL \", ExportString[final, \"RawJSON\", \"Compact\" -> True]];"
}, "\n"], "Text"];

result = FeynmanTrick`RunIntegrationPipeline["bubble",
  "Runner" -> runner,
  "WorkingDirectory" -> repoRoot,
  "FIREPath" -> FileNameJoin[{$TemporaryDirectory, "ft-pipeline-fire"}],
  "PreparedCacheDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-prep"}],
  "CheckpointDirectory" -> FileNameJoin[{$TemporaryDirectory,
    "ft-pipeline-process-checkpoints"}],
  "EpsilonOrder" -> 2];

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
    result["Final", "Finite"] === 7 &&
    result["Final", "FIREPath"] ===
      ExpandFileName[FileNameJoin[{$TemporaryDirectory, "ft-pipeline-fire"}]]];
assert["pipeline result retains every requested positive epsilon coefficient",
  Length[result["Stepwise"]] === 1 &&
    result["Stepwise"][[1, "Level"]] === 1 &&
    result["Plan", "Settings", "EpsilonOrder"] === 2 &&
    result["Stepwise"][[1, "Coefficients"]][[All, 1]] === {0, 1, 2} &&
    result["Final", "Coefficients"][[All, 1]] === {0, 1, 2}];
assert["an existing typed plan can be executed",
  plannedResult["Schema"] === "FeynmanTrick.PipelineResult/v1" &&
    plannedResult["Status"] === "Succeeded" &&
    plannedResult["Final", "Example"] === "bubble" &&
    !KeyExistsQ[plannedResult["Final"], "Coefficients"]];

If[FileExistsQ[runner], DeleteFile[runner]];
Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
