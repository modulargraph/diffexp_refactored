(* FeynmanTrick/DiffExp2Pipeline.m -- release facade for the tested
   Scripts/run_ft_stepwise2.m ladder.

   The current minimal slice deliberately executes the mature ladder in a
   clean Wolfram subprocess.  This gives package users typed plans, stable
   result records, checkpoint/resume, and safe argv/environment construction
   now, while the numerical level loop is moved out of the script in later
   slices.  No shell command is constructed. *)

BeginPackage["FeynmanTrick`DiffExp2Pipeline`"];

PipelinePlan::usage = "PipelinePlan[example, opts] returns the reproducible command/environment/checkpoint plan for the DiffExp2 Feynman-trick ladder without running it.";
RunIntegrationPipeline::usage = "RunIntegrationPipeline[example, opts] runs the DiffExp2 Feynman-trick ladder with the C++ recurrence backend by default and returns a PipelineResult (or PipelineProcess when Asynchronous -> True).";
ResumeIntegrationPipeline::usage = "ResumeIntegrationPipeline[example, checkpoint, opts] resumes the DiffExp2 ladder from an atomic transport or boundary checkpoint.";
RunnerSettingsFromEnvironment::usage = "RunnerSettingsFromEnvironment[] parses and validates the runner environment. It is shared by the package facade and Scripts/run_ft_stepwise2.m.";

Begin["`Private`"];

$pipelinePlanSchema = "FeynmanTrick.PipelinePlan/v1";
$pipelineResultSchema = "FeynmanTrick.PipelineResult/v1";
$pipelineProcessSchema = "FeynmanTrick.PipelineProcess/v1";
$packageRoot = ParentDirectory[DirectoryName[$InputFileName]];

failure[detail_String, data_:<||>] := Failure["FeynmanTrickPipeline",
  Join[<|"Detail" -> detail|>, data]];

envOrDefault[name_String, default_String] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0,
    StringTrim[value], default]];

parseInteger[name_String, text_String, min_Integer] := Module[{trim = StringTrim[text], value},
  If[!StringMatchQ[trim, RegularExpression["[+-]?[0-9]+"]],
    Return[failure["environment value must be an integer",
      <|"Name" -> name, "Value" -> text|>], Module]];
  value = ToExpression[trim];
  If[!IntegerQ[value] || value < min,
    failure["environment integer is below its allowed minimum",
      <|"Name" -> name, "Value" -> value, "Minimum" -> min|>], value]];

parsePositiveRational[name_String, text_String] := Module[{trim = StringTrim[text], value},
  If[!StringMatchQ[trim, RegularExpression["[0-9]+(/[0-9]+)?"]],
    Return[failure["environment value must be a positive integer or rational",
      <|"Name" -> name, "Value" -> text|>], Module]];
  value = ToExpression[trim];
  If[!(Head[value] === Rational || IntegerQ[value]) || value <= 0,
    failure["environment value must be positive",
      <|"Name" -> name, "Value" -> value|>], value]];

parseFlag[name_String, text_String] := Switch[StringTrim[text],
  "0", False,
  "1", True,
  _, failure["environment flag must be 0 or 1",
    <|"Name" -> name, "Value" -> text|>]];

parseHalos[text_String] := Module[{parts, values},
  parts = StringTrim /@ StringSplit[text, ","];
  If[parts === {""}, Return[{}, Module]];
  values = Map[parseInteger["FT_LEVEL_EPS_HALOS", #, 0] &, parts];
  If[AnyTrue[values, FailureQ], First[Select[values, FailureQ]], values]];

RunnerSettingsFromEnvironment[] := Module[
  {backend, threads, wp, epsOrder, expansionOrder, boundaryExtraOrder,
   divisionOrder, requestedStepDivisionOrder, radius, halos, stop, singular,
   batch, rebuild, allowStale, fireTimeout, resume, checkpointDir, prepRoot,
   values},
  backend = envOrDefault["DE2_RECURRENCE_BACKEND", "Cpp"];
  If[!MemberQ[{"Cpp", "Wolfram"}, backend],
    Return[failure["DE2_RECURRENCE_BACKEND must be Cpp or Wolfram",
      <|"Value" -> backend|>], Module]];
  values = {
    threads = parseInteger["DE2_CPP_THREADS",
      envOrDefault["DE2_CPP_THREADS", "4"], 1],
    wp = parseInteger["FT_WORKING_PRECISION",
      envOrDefault["FT_WORKING_PRECISION", "500"], 20],
    epsOrder = parseInteger["FT_EPS_ORDER",
      envOrDefault["FT_EPS_ORDER", "0"], 0],
    expansionOrder = parseInteger["FT_EXPANSION_ORDER",
      envOrDefault["FT_EXPANSION_ORDER", "50"], 10],
    boundaryExtraOrder = parseInteger["FT_BOUNDARY_EXTRA_ORDER",
      envOrDefault["FT_BOUNDARY_EXTRA_ORDER", "4"], 0],
    fireTimeout = parseInteger["FT_FIRE_TIMEOUT_SECONDS",
      envOrDefault["FT_FIRE_TIMEOUT_SECONDS", "1800"], 1],
    divisionOrder = parseInteger["FT_DIVISION_ORDER",
      envOrDefault["FT_DIVISION_ORDER", "3"], 2],
    radius = parsePositiveRational["FT_RADIUS_OF_CONVERGENCE",
      envOrDefault["FT_RADIUS_OF_CONVERGENCE", "1"]],
    singular = parseFlag["DE2_SINGULAR_MATCH_PRECONDITION",
      envOrDefault["DE2_SINGULAR_MATCH_PRECONDITION", "0"]],
    batch = parseFlag["FT_CPP_BATCH_ENDPOINT_ARMS",
      envOrDefault["FT_CPP_BATCH_ENDPOINT_ARMS", "1"]],
    rebuild = parseFlag["FT_REBUILD_PREP",
      envOrDefault["FT_REBUILD_PREP", "0"]],
    allowStale = parseFlag["FT_ALLOW_STALE_LADDER_CHECKPOINT",
      envOrDefault["FT_ALLOW_STALE_LADDER_CHECKPOINT", "0"]],
    halos = parseHalos[envOrDefault["FT_LEVEL_EPS_HALOS", "0"]]
  };
  If[AnyTrue[values, FailureQ], Return[First[Select[values, FailureQ]], Module]];
  requestedStepDivisionOrder = parseInteger["FT_STEP_DIVISION_ORDER",
    envOrDefault["FT_STEP_DIVISION_ORDER", ToString[divisionOrder, InputForm]], 1];
  If[FailureQ[requestedStepDivisionOrder], Return[requestedStepDivisionOrder, Module]];
  stop = envOrDefault["FT_STOP_AFTER_BOUNDARY_LEVEL", ""];
  stop = If[stop === "", Missing["NotSet"],
    parseInteger["FT_STOP_AFTER_BOUNDARY_LEVEL", stop, 0]];
  If[FailureQ[stop], Return[stop, Module]];
  resume = envOrDefault["FT_RESUME_LADDER_CHECKPOINT", ""];
  checkpointDir = envOrDefault["FT_LADDER_CHECKPOINT_DIR",
    If[resume === "", "", DirectoryName[ExpandFileName[resume]]]];
  prepRoot = envOrDefault["FT_PREP_CACHE_DIR",
    FileNameJoin[{$TemporaryDirectory, "DiffExp2_FT_Prepared"}]];
  <|
    "RecurrenceBackend" -> backend, "CppThreads" -> threads,
    "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
    "ExpansionOrder" -> expansionOrder,
    "BoundaryExtraOrder" -> boundaryExtraOrder,
    "FIRETimeoutSeconds" -> fireTimeout,
    "DivisionOrder" -> divisionOrder,
    "RequestedStepDivisionOrder" -> requestedStepDivisionOrder,
    (* The classic planner couples placement and +/-1/k matching. *)
    "StepDivisionOrder" -> divisionOrder,
    "RadiusOfConvergence" -> radius, "LevelEpsilonHalos" -> halos,
    "StopAfterBoundaryLevel" -> stop,
    "SingularMatchPrecondition" -> singular,
    "BatchEndpointArms" -> (batch && backend === "Cpp"),
    "PrepCacheRoot" -> prepRoot, "ForcePrepRebuild" -> rebuild,
    "ResumeCheckpoint" -> resume, "CheckpointDirectory" -> checkpointDir,
    "AllowStaleCheckpoint" -> allowStale
  |>];

Options[PipelinePlan] = {
  "WorkingPrecision" -> 500,
  "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 0,
  "BoundaryExtraOrder" -> 4,
  "LevelEpsilonHalos" -> {0},
  "DivisionOrder" -> 3,
  "RadiusOfConvergence" -> 1,
  "RecurrenceBackend" -> "Cpp",
  "CppThreads" -> Automatic,
  "ValueTransport" -> True,
  "BatchEndpointArms" -> True,
  "SingularMatchPrecondition" -> False,
  "PreparedCacheDirectory" -> Automatic,
  "CheckpointDirectory" -> Automatic,
  "ResumeFrom" -> None,
  "RebuildPreparation" -> False,
  "AllowStaleCheckpoint" -> False,
  "StopAfterBoundaryLevel" -> None,
  "FIRETimeoutSeconds" -> 1800,
  "Runner" -> Automatic,
  "WolframScript" -> Automatic,
  "WorkingDirectory" -> Automatic,
  "ExtraEnvironment" -> <||>,
  "Asynchronous" -> False
};

boolQ[value_] := MemberQ[{True, False}, value];
boolString[value_] := If[TrueQ[value], "1", "0"];
inputString[value_] := ToString[value, InputForm];

resolveExecutable[value_] := Module[{candidate},
  candidate = Replace[value, Automatic -> "wolframscript"];
  If[!StringQ[candidate] || StringLength[StringTrim[candidate]] === 0,
    Return[failure["\"WolframScript\" must name an executable",
      <|"Value" -> value|>], Module]];
  candidate = StringTrim[candidate];
  (* RunProcess resolves a bare executable name through the operating-system
     process search path.  An explicit path can and should be checked now. *)
  If[MemberQ[{"", "."}, DirectoryName[candidate]], candidate,
    candidate = ExpandFileName[candidate];
    If[FileExistsQ[candidate], candidate,
      failure["the explicit wolframscript executable does not exist",
        <|"Executable" -> candidate|>]]]];

validatePlanOptions[settings_Association] := Module[{checks},
  checks = {
    IntegerQ[settings["WorkingPrecision"]] && settings["WorkingPrecision"] >= 20,
    IntegerQ[settings["ExpansionOrder"]] && settings["ExpansionOrder"] >= 10,
    IntegerQ[settings["EpsilonOrder"]] && settings["EpsilonOrder"] >= 0,
    IntegerQ[settings["BoundaryExtraOrder"]] && settings["BoundaryExtraOrder"] >= 0,
    ListQ[settings["LevelEpsilonHalos"]] &&
      AllTrue[settings["LevelEpsilonHalos"], IntegerQ[#] && # >= 0 &],
    IntegerQ[settings["DivisionOrder"]] && settings["DivisionOrder"] >= 2,
    (IntegerQ[settings["RadiusOfConvergence"]] ||
      Head[settings["RadiusOfConvergence"]] === Rational) &&
      settings["RadiusOfConvergence"] > 0,
    MemberQ[{"Cpp", "Wolfram"}, settings["RecurrenceBackend"]],
    IntegerQ[settings["CppThreads"]] && settings["CppThreads"] >= 1,
    And @@ (boolQ[settings[#]] & /@ {
      "ValueTransport", "BatchEndpointArms", "SingularMatchPrecondition",
      "RebuildPreparation", "AllowStaleCheckpoint", "Asynchronous"}),
    IntegerQ[settings["FIRETimeoutSeconds"]] &&
      settings["FIRETimeoutSeconds"] >= 1,
    AssociationQ[settings["ExtraEnvironment"]] &&
      AllTrue[Keys[settings["ExtraEnvironment"]], StringQ] &&
      AllTrue[Values[settings["ExtraEnvironment"]], StringQ]
  };
  If[And @@ checks, True,
    failure["one or more pipeline options are invalid", <|"Settings" -> settings|>]]];

PipelinePlan[example_String, OptionsPattern[]] := Module[
  {threads, runner, executable, workdir, prep, checkpoint, resume, stop,
   settings, valid, env},
  If[StringLength[StringTrim[example]] === 0 || StringContainsQ[example, ","],
    Return[failure["example must be one nonempty registry name",
      <|"Example" -> example|>], Module]];
  threads = Replace[OptionValue["CppThreads"], Automatic :>
    Max[1, Min[10, $ProcessorCount]]];
  runner = ExpandFileName[Replace[OptionValue["Runner"], Automatic :>
    FileNameJoin[{$packageRoot, "Scripts", "run_ft_stepwise2.m"}]]];
  workdir = ExpandFileName[Replace[OptionValue["WorkingDirectory"],
    Automatic -> $packageRoot]];
  prep = ExpandFileName[Replace[OptionValue["PreparedCacheDirectory"],
    Automatic :> FileNameJoin[{$TemporaryDirectory, "DiffExp2_FT_Prepared"}]]];
  checkpoint = ExpandFileName[Replace[OptionValue["CheckpointDirectory"],
    Automatic :> FileNameJoin[{$TemporaryDirectory,
      "DiffExp2_FT_Checkpoints", example}]]];
  resume = OptionValue["ResumeFrom"];
  If[resume =!= None,
    If[!StringQ[resume] || !FileExistsQ[resume],
      Return[failure["resume checkpoint does not exist",
        <|"Checkpoint" -> resume|>], Module]];
    resume = ExpandFileName[resume];
    If[OptionValue["CheckpointDirectory"] === Automatic,
      checkpoint = DirectoryName[resume]]];
  stop = OptionValue["StopAfterBoundaryLevel"];
  settings = <|
    "WorkingPrecision" -> OptionValue["WorkingPrecision"],
    "ExpansionOrder" -> OptionValue["ExpansionOrder"],
    "EpsilonOrder" -> OptionValue["EpsilonOrder"],
    "BoundaryExtraOrder" -> OptionValue["BoundaryExtraOrder"],
    "LevelEpsilonHalos" -> OptionValue["LevelEpsilonHalos"],
    "DivisionOrder" -> OptionValue["DivisionOrder"],
    "RadiusOfConvergence" -> OptionValue["RadiusOfConvergence"],
    "RecurrenceBackend" -> OptionValue["RecurrenceBackend"],
    "CppThreads" -> threads,
    "ValueTransport" -> OptionValue["ValueTransport"],
    "BatchEndpointArms" -> OptionValue["BatchEndpointArms"],
    "SingularMatchPrecondition" -> OptionValue["SingularMatchPrecondition"],
    "RebuildPreparation" -> OptionValue["RebuildPreparation"],
    "AllowStaleCheckpoint" -> OptionValue["AllowStaleCheckpoint"],
    "FIRETimeoutSeconds" -> OptionValue["FIRETimeoutSeconds"],
    "Asynchronous" -> OptionValue["Asynchronous"],
    "ExtraEnvironment" -> OptionValue["ExtraEnvironment"]
  |>;
  valid = validatePlanOptions[settings];
  If[FailureQ[valid], Return[valid, Module]];
  If[stop =!= None && !(IntegerQ[stop] && stop >= 0),
    Return[failure["\"StopAfterBoundaryLevel\" must be None or a nonnegative integer",
      <|"Value" -> stop|>], Module]];
  If[!FileExistsQ[runner],
    Return[failure["pipeline runner does not exist", <|"Runner" -> runner|>], Module]];
  If[!DirectoryQ[workdir],
    Return[failure["pipeline working directory does not exist",
      <|"WorkingDirectory" -> workdir|>], Module]];
  executable = resolveExecutable[OptionValue["WolframScript"]];
  If[FailureQ[executable], Return[executable, Module]];
  env = <|
    "FT_EXAMPLES" -> example,
    "DE2_RECURRENCE_BACKEND" -> settings["RecurrenceBackend"],
    "DE2_CPP_THREADS" -> inputString[threads],
    "DE2_VALUE_TRANSPORT" -> boolString[settings["ValueTransport"]],
    "DE2_SINGULAR_MATCH_PRECONDITION" ->
      boolString[settings["SingularMatchPrecondition"]],
    "FT_CPP_BATCH_ENDPOINT_ARMS" ->
      boolString[settings["BatchEndpointArms"]],
    "FT_WORKING_PRECISION" -> inputString[settings["WorkingPrecision"]],
    "FT_EXPANSION_ORDER" -> inputString[settings["ExpansionOrder"]],
    "FT_EPS_ORDER" -> inputString[settings["EpsilonOrder"]],
    "FT_BOUNDARY_EXTRA_ORDER" -> inputString[settings["BoundaryExtraOrder"]],
    "FT_LEVEL_EPS_HALOS" -> StringRiffle[
      inputString /@ settings["LevelEpsilonHalos"], ","],
    "FT_DIVISION_ORDER" -> inputString[settings["DivisionOrder"]],
    "FT_STEP_DIVISION_ORDER" -> inputString[settings["DivisionOrder"]],
    "FT_RADIUS_OF_CONVERGENCE" ->
      inputString[settings["RadiusOfConvergence"]],
    "FT_PREP_CACHE_DIR" -> prep,
    "FT_LADDER_CHECKPOINT_DIR" -> checkpoint,
    "FT_REBUILD_PREP" -> boolString[settings["RebuildPreparation"]],
    "FT_ALLOW_STALE_LADDER_CHECKPOINT" ->
      boolString[settings["AllowStaleCheckpoint"]],
    "FT_FIRE_TIMEOUT_SECONDS" -> inputString[settings["FIRETimeoutSeconds"]]
  |>;
  If[resume =!= None, env["FT_RESUME_LADDER_CHECKPOINT"] = resume];
  If[stop =!= None,
    env["FT_STOP_AFTER_BOUNDARY_LEVEL"] = inputString[stop]];
  (* ExtraEnvironment is for additive diagnostics (DEBUG_* etc.).  Canonical
     plan settings win on duplicate keys so the visible Settings record can
     never disagree with the process environment actually executed. *)
  env = Join[settings["ExtraEnvironment"], env];
  <|
    "Schema" -> $pipelinePlanSchema, "Example" -> example,
    "Command" -> {executable, "-file", runner},
    "WorkingDirectory" -> workdir, "Environment" -> env,
    "Settings" -> KeyDrop[settings, {"ExtraEnvironment"}],
    "PreparedCacheDirectory" -> prep,
    "CheckpointDirectory" -> checkpoint, "ResumeFrom" -> resume
  |>];

jsonRows[output_String, prefix_String] := DeleteCases[
  Map[Function[line, If[StringStartsQ[line, prefix],
      Quiet[Check[ImportString[StringDrop[line, StringLength[prefix]],
        "RawJSON"], Missing["InvalidJSON", line]]], Nothing]],
    StringSplit[output, {"\r\n", "\n", "\r"}]],
  _Missing];

pipelineResult[plan_Association, process_Association] := Module[
  {stdout, stderr, exit, finals, steps, stopped, status},
  stdout = Lookup[process, "StandardOutput", ""];
  stderr = Lookup[process, "StandardError", ""];
  exit = Lookup[process, "ExitCode", Missing["NotAvailable"]];
  steps = jsonRows[stdout, "STEPWISE "];
  finals = jsonRows[stdout, "FINAL "];
  stopped = StringContainsQ[stdout, "STOPPED_AFTER_BOUNDARY_LEVEL"];
  status = Which[
    exit =!= 0, "Failed",
    stopped, "Stopped",
    finals =!= {}, "Succeeded",
    True, "Incomplete"];
  <|
    "Schema" -> $pipelineResultSchema, "Status" -> status,
    "ExitCode" -> exit, "Final" -> If[finals === {},
      Missing["NotProduced"], Last[finals]],
    "Stepwise" -> steps, "Plan" -> plan,
    "StandardOutput" -> stdout, "StandardError" -> stderr
  |>];

processCommand[plan_Association] := Module[{envProgram, assignments},
  envProgram = If[$OperatingSystem === "Windows", None, "/usr/bin/env"];
  If[!StringQ[envProgram] || !FileExistsQ[envProgram],
    Return[failure[
      "the pipeline facade currently requires /usr/bin/env for inherited-environment execution",
      <|"OperatingSystem" -> $OperatingSystem|>], Module]];
  assignments = KeyValueMap[#1 <> "=" <> #2 &,
    KeySort[plan["Environment"]]];
  Join[{envProgram}, assignments, plan["Command"]]];

runPlan[plan_, asynchronous_] := Module[{process, command},
  If[FailureQ[plan], Return[plan, Module]];
  command = processCommand[plan];
  If[FailureQ[command], Return[command, Module]];
  If[TrueQ[asynchronous],
    process = Quiet[Check[StartProcess[command,
      ProcessDirectory -> plan["WorkingDirectory"]], $Failed]];
    If[Head[process] =!= ProcessObject,
      Return[failure["could not start the pipeline subprocess",
        <|"Plan" -> plan|>], Module]];
    <|"Schema" -> $pipelineProcessSchema, "Status" -> "Running",
      "Process" -> process, "Plan" -> plan|>,
    process = Quiet[Check[RunProcess[command, All,
      ProcessDirectory -> plan["WorkingDirectory"]], $Failed]];
    If[!AssociationQ[process],
      Return[failure["pipeline subprocess could not be executed",
        <|"Plan" -> plan|>], Module]];
    pipelineResult[plan, process]]];

Options[RunIntegrationPipeline] = Options[PipelinePlan];
RunIntegrationPipeline[example_String, opts:OptionsPattern[]] := Module[{plan},
  plan = PipelinePlan[example, opts];
  runPlan[plan, OptionValue["Asynchronous"]]];

Options[ResumeIntegrationPipeline] = Options[PipelinePlan];
ResumeIntegrationPipeline[example_String, checkpoint_String,
    opts:OptionsPattern[]] := Module[{rules, plan},
  rules = DeleteCases[{opts},
    HoldPattern[Rule["ResumeFrom", _] | RuleDelayed["ResumeFrom", _]]];
  plan = PipelinePlan[example,
    Sequence @@ Join[rules, {"ResumeFrom" -> checkpoint}]];
  runPlan[plan, OptionValue["Asynchronous"]]];

End[];
EndPackage[];
