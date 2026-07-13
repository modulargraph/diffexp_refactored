(* FeynmanTrick/DiffExp2Pipeline.m -- release facade for the tested
   Scripts/run_ft_stepwise2.m ladder.

   The current minimal slice deliberately executes the mature ladder in a
   clean Wolfram subprocess.  This gives package users typed plans, stable
   result records, checkpoint/resume, and safe argv/environment construction
   now, while the numerical level loop is moved out of the script in later
   slices.  No shell command is constructed. *)

BeginPackage["FeynmanTrick`DiffExp2Pipeline`"];

PipelinePlan::usage = "PipelinePlan[example, opts] returns the reproducible command/environment/checkpoint plan for the DiffExp2 Feynman-trick ladder without running it. PipelinePlan[familySpec, integrals, opts] additionally creates an exact content-addressed custom-family request; integrals may be one index vector, an ordered list, or All for deterministic L0 master discovery at execution.";
RunIntegrationPipeline::usage = "RunIntegrationPipeline[example, opts] runs the DiffExp2 Feynman-trick ladder with the C++ recurrence backend by default. RunIntegrationPipeline[familySpec, integrals, opts] executes selected targets or execution-time All-master discovery with the built-in runner, or with an external runner explicitly declared request-aware. RunIntegrationPipeline[plan] executes an existing PipelinePlan. It returns a PipelineResult (or PipelineProcess when Asynchronous -> True).";
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

parseDeltaPrescriptionSign[name_String, text_String] := Module[
  {trim = StringTrim[text], value},
  If[!StringMatchQ[trim, RegularExpression["[+-]?1"]],
    Return[failure["environment delta-prescription sign must be +1 or -1",
      <|"Name" -> name, "Value" -> text|>], Module]];
  value = ToExpression[trim];
  If[MemberQ[{-1, 1}, value], value,
    failure["environment delta-prescription sign must be +1 or -1",
      <|"Name" -> name, "Value" -> text|>]]];

parseHalos[text_String] := Module[{parts, values},
  parts = StringTrim /@ StringSplit[text, ","];
  If[parts === {""}, Return[{}, Module]];
  values = Map[parseInteger["FT_LEVEL_EPS_HALOS", #, 0] &, parts];
  If[AnyTrue[values, FailureQ], First[Select[values, FailureQ]], values]];

RunnerSettingsFromEnvironment[] := Module[
  {backend, threads, wp, epsOrder, expansionOrder, boundaryExtraOrder,
   divisionOrder, requestedStepDivisionOrder, radius, halos, stop, singular,
   deltaPrescriptionSign, batch, rebuild, migrateLegacyPrep, allowStale,
   saveNativeTransportCheckpoint,
   fireTimeout, firePath,
   resume, checkpointDir, prepRoot, values},
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
    deltaPrescriptionSign = parseDeltaPrescriptionSign[
      "FT_DELTA_PRESCRIPTION_SIGN",
      envOrDefault["FT_DELTA_PRESCRIPTION_SIGN", "1"]],
    batch = parseFlag["FT_CPP_BATCH_ENDPOINT_ARMS",
      envOrDefault["FT_CPP_BATCH_ENDPOINT_ARMS", "1"]],
    rebuild = parseFlag["FT_REBUILD_PREP",
      envOrDefault["FT_REBUILD_PREP", "0"]],
    migrateLegacyPrep = parseFlag["FT_MIGRATE_LEGACY_PREP",
      envOrDefault["FT_MIGRATE_LEGACY_PREP", "0"]],
    allowStale = parseFlag["FT_ALLOW_STALE_LADDER_CHECKPOINT",
      envOrDefault["FT_ALLOW_STALE_LADDER_CHECKPOINT", "0"]],
    saveNativeTransportCheckpoint = parseFlag[
      "FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT",
      envOrDefault["FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT", "0"]],
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
  firePath = ExpandFileName[envOrDefault["FT_FIRE_PATH",
    FileNameJoin[{$packageRoot, "Dependencies", "fire", "FIRE6"}]]];
  <|
    "RecurrenceBackend" -> backend, "CppThreads" -> threads,
    "WorkingPrecision" -> wp, "EpsilonOrder" -> epsOrder,
    "ExpansionOrder" -> expansionOrder,
    "BoundaryExtraOrder" -> boundaryExtraOrder,
    "FIRETimeoutSeconds" -> fireTimeout,
    "FIREPath" -> firePath,
    "DivisionOrder" -> divisionOrder,
    "RequestedStepDivisionOrder" -> requestedStepDivisionOrder,
    (* The classic planner couples placement and +/-1/k matching. *)
    "StepDivisionOrder" -> divisionOrder,
    "RadiusOfConvergence" -> radius, "LevelEpsilonHalos" -> halos,
    "StopAfterBoundaryLevel" -> stop,
    "SingularMatchPrecondition" -> singular,
    "DeltaPrescriptionSign" -> deltaPrescriptionSign,
    "BatchEndpointArms" -> (batch && backend === "Cpp"),
    "PrepCacheRoot" -> prepRoot, "ForcePrepRebuild" -> rebuild,
    "MigrateLegacyPreparation" -> migrateLegacyPrep,
    "ResumeCheckpoint" -> resume, "CheckpointDirectory" -> checkpointDir,
    "AllowStaleCheckpoint" -> allowStale,
    "SaveNativeTransportCheckpoint" -> saveNativeTransportCheckpoint
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
  "DeltaPrescriptionSign" -> 1,
  "PreparedCacheDirectory" -> Automatic,
  "FIREPath" -> Automatic,
  "CheckpointDirectory" -> Automatic,
  "ResumeFrom" -> None,
  "RebuildPreparation" -> False,
  "MigrateLegacyPreparation" -> False,
  "AllowStaleCheckpoint" -> False,
  "SaveNativeTransportCheckpoint" -> False,
  "StopAfterBoundaryLevel" -> None,
  "FIRETimeoutSeconds" -> 1800,
  "Runner" -> Automatic,
  "WolframScript" -> Automatic,
  "WorkingDirectory" -> Automatic,
  "RequestAwareRunner" -> False,
  "RequestDirectory" -> Automatic,
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
    MemberQ[{-1, 1}, settings["DeltaPrescriptionSign"]],
    IntegerQ[settings["CppThreads"]] && settings["CppThreads"] >= 1,
    And @@ (boolQ[settings[#]] & /@ {
      "ValueTransport", "BatchEndpointArms", "SingularMatchPrecondition",
      "RebuildPreparation", "MigrateLegacyPreparation",
      "AllowStaleCheckpoint", "SaveNativeTransportCheckpoint",
      "Asynchronous"}),
    IntegerQ[settings["FIRETimeoutSeconds"]] &&
      settings["FIRETimeoutSeconds"] >= 1,
    AssociationQ[settings["ExtraEnvironment"]] &&
      AllTrue[Keys[settings["ExtraEnvironment"]], StringQ] &&
      AllTrue[Values[settings["ExtraEnvironment"]], StringQ],
    StringQ[settings["FIREPath"]] &&
      StringLength[StringTrim[settings["FIREPath"]]] > 0
  };
  If[And @@ checks, True,
    failure["one or more pipeline options are invalid", <|"Settings" -> settings|>]]];

buildPipelinePlan[example_String, registryQ_, OptionsPattern[PipelinePlan]] := Module[
  {threads, runner, executable, workdir, prep, firePath, firePathOption,
   checkpoint, resume, stop,
   settings, valid, env},
  If[StringLength[StringTrim[example]] === 0 || StringContainsQ[example, ","],
    Return[failure["example must be one nonempty registry name",
      <|"Example" -> example|>], Module]];
  If[TrueQ[registryQ] &&
      !MemberQ[FeynmanTrick`SupportedExamples[], example],
    Return[failure["example is not present in the Feynman-trick registry",
      <|"Example" -> example,
        "SupportedExamples" -> FeynmanTrick`SupportedExamples[]|>], Module]];
  threads = Replace[OptionValue["CppThreads"], Automatic :>
    Max[1, Min[10, $ProcessorCount]]];
  runner = ExpandFileName[Replace[OptionValue["Runner"], Automatic :>
    FileNameJoin[{$packageRoot, "Scripts", "run_ft_stepwise2.m"}]]];
  workdir = ExpandFileName[Replace[OptionValue["WorkingDirectory"],
    Automatic -> $packageRoot]];
  prep = ExpandFileName[Replace[OptionValue["PreparedCacheDirectory"],
    Automatic :> FileNameJoin[{$TemporaryDirectory, "DiffExp2_FT_Prepared"}]]];
  firePathOption = Replace[OptionValue["FIREPath"],
    Automatic :> FileNameJoin[{$packageRoot, "Dependencies", "fire", "FIRE6"}]];
  If[!StringQ[firePathOption] || StringLength[StringTrim[firePathOption]] === 0,
    Return[failure["\"FIREPath\" must be a nonempty path string",
      <|"Value" -> OptionValue["FIREPath"]|>], Module]];
  firePath = ExpandFileName[StringTrim[firePathOption]];
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
    "DeltaPrescriptionSign" -> OptionValue["DeltaPrescriptionSign"],
    "RebuildPreparation" -> OptionValue["RebuildPreparation"],
    "MigrateLegacyPreparation" ->
      OptionValue["MigrateLegacyPreparation"],
    "AllowStaleCheckpoint" -> OptionValue["AllowStaleCheckpoint"],
    "SaveNativeTransportCheckpoint" ->
      OptionValue["SaveNativeTransportCheckpoint"],
    "FIREPath" -> firePath,
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
    (* /usr/bin/env inherits the parent process.  Empty canonical assignments
       are therefore part of a reproducible registry plan: a stale custom
       request or checkpoint control variable must not change its mode. *)
    "FT_FAMILY_REQUEST_FILE" -> "",
    "FT_FAMILY_REQUEST_ID" -> "",
    "FT_RESUME_LADDER_CHECKPOINT" ->
      If[resume === None, "", resume],
    "FT_STOP_AFTER_BOUNDARY_LEVEL" ->
      If[stop === None, "", inputString[stop]],
    "DE2_RECURRENCE_BACKEND" -> settings["RecurrenceBackend"],
    "DE2_CPP_THREADS" -> inputString[threads],
    "DE2_VALUE_TRANSPORT" -> boolString[settings["ValueTransport"]],
    "DE2_SINGULAR_MATCH_PRECONDITION" ->
      boolString[settings["SingularMatchPrecondition"]],
    "FT_DELTA_PRESCRIPTION_SIGN" ->
      inputString[settings["DeltaPrescriptionSign"]],
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
    "FT_FIRE_PATH" -> firePath,
    "FT_LADDER_CHECKPOINT_DIR" -> checkpoint,
    "FT_REBUILD_PREP" -> boolString[settings["RebuildPreparation"]],
    "FT_MIGRATE_LEGACY_PREP" ->
      boolString[settings["MigrateLegacyPreparation"]],
    "FT_ALLOW_STALE_LADDER_CHECKPOINT" ->
      boolString[settings["AllowStaleCheckpoint"]],
    "FT_SAVE_NATIVE_TRANSPORT_CHECKPOINT" ->
      boolString[settings["SaveNativeTransportCheckpoint"]],
    "FT_FIRE_TIMEOUT_SECONDS" -> inputString[settings["FIRETimeoutSeconds"]]
  |>;
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

PipelinePlan[example_String, opts:OptionsPattern[]] :=
  buildPipelinePlan[example, True, opts];

familyPipelinePlan[family_Association, targets_,
    opts:OptionsPattern[PipelinePlan]] := Module[
  {canonical, request, requestDirectoryOption, requestDirectory,
   requestFile, executionName, plan, runnerSupport, policy, activeCount},
  canonical = If[targets === Automatic,
    FeynmanTrick`FamilySpec`CreateFamily[family],
    FeynmanTrick`FamilySpec`CreateFamily[family, targets]];
  If[canonical === $Failed || !AssociationQ[canonical],
    Return[failure["custom family normalization failed"], Module]];
  If[AnyTrue[{"AnalyticPrescription", "Prescriptions",
        "KinematicAssumptions"},
      KeyExistsQ[canonical["Definition"], #] &],
    Return[failure[
      "custom analytic-prescription and kinematic-assumption fields are not yet wired into the production runner and cannot be ignored safely",
      <|"FamilyID" -> canonical["FamilyID"]|>], Module]];
  activeCount = canonical["NumPropagators"] -
    Length[canonical["EliminatedPositions"]];
  If[activeCount <= 1 ||
      Length[canonical["CombinationSequence"]] =!= activeCount - 1,
    Return[failure[
      "the production ladder requires a complete nonempty combination sequence ending in one active propagator",
      <|"FamilyID" -> canonical["FamilyID"],
        "ActivePropagators" -> activeCount,
        "CombinationSequence" -> canonical["CombinationSequence"]|>],
      Module]];
  request = FeynmanTrick`PipelineRequest`CreatePipelineRequest[canonical];
  If[FailureQ[request], Return[request, Module]];
  If[!boolQ[OptionValue["RequestAwareRunner"]],
    Return[failure["\"RequestAwareRunner\" must be True or False",
      <|"Value" -> OptionValue["RequestAwareRunner"]|>], Module]];
  requestDirectoryOption = OptionValue["RequestDirectory"];
  If[requestDirectoryOption === Automatic,
    requestDirectoryOption = FileNameJoin[{
      $TemporaryDirectory, "DiffExp2_FT_Requests"}]];
  If[!StringQ[requestDirectoryOption] ||
      StringLength[StringTrim[requestDirectoryOption]] === 0,
    Return[failure["\"RequestDirectory\" must be a nonempty path string",
      <|"Value" -> OptionValue["RequestDirectory"]|>], Module]];
  requestDirectory = ExpandFileName[StringTrim[requestDirectoryOption]];
  requestFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
    request, requestDirectory];
  If[FailureQ[requestFile], Return[requestFile, Module]];
  executionName = "family_" <> StringTake[request["RequestID"], -16];
  plan = buildPipelinePlan[executionName, False, opts];
  If[FailureQ[plan], Return[plan, Module]];
  runnerSupport = Which[
    OptionValue["Runner"] === Automatic,
      "BuiltInRequestAwareRunner",
    TrueQ[OptionValue["RequestAwareRunner"]],
      "ExplicitRequestAwareRunner",
    True,
      "UnverifiedExternalRunner"
  ];
  policy = Join[request["ExecutionPolicy"], <|
    "RunnerSupport" -> runnerSupport,
    "ExecutionReady" -> (
      TrueQ[request["ExecutionPolicy", "ExecutionReady"]] &&
      MemberQ[{"BuiltInRequestAwareRunner", "ExplicitRequestAwareRunner"},
        runnerSupport])|>];
  plan["Environment"] = Join[plan["Environment"], <|
    "FT_FAMILY_REQUEST_FILE" -> requestFile,
    "FT_FAMILY_REQUEST_ID" -> request["RequestID"]|>];
  Join[plan, <|
    "InputKind" -> "Family",
    "FamilyName" -> canonical["Name"],
    "FamilyID" -> canonical["FamilyID"],
    "Request" -> request,
    "RequestID" -> request["RequestID"],
    "RequestDirectory" -> requestDirectory,
    "RequestFile" -> requestFile,
    "OutputRequests" -> request["OutputRequests"],
    "ExecutionPolicy" -> policy
  |>]
];

PipelinePlan[family_Association, targets:(All | _List),
    opts:OptionsPattern[]] := familyPipelinePlan[family, targets, opts];

PipelinePlan[family_Association, opts:OptionsPattern[]] :=
  familyPipelinePlan[family, Automatic, opts];

jsonRowParse[output_String, prefix_String] := Module[
  {taggedLines, parsed, invalidMask},
  taggedLines = Select[StringSplit[output, {"\r\n", "\n", "\r"}],
    StringStartsQ[#, prefix] &];
  parsed = Map[
    Quiet[Check[
      ImportString[StringDrop[#, StringLength[prefix]], "RawJSON"],
      $Failed]] &,
    taggedLines];
  invalidMask = SameQ[#, $Failed] & /@ parsed;
  <|
    "Rows" -> Pick[parsed, Not /@ invalidMask],
    "Errors" -> Pick[taggedLines, invalidMask]
  |>
];

jsonRows[output_String, prefix_String] :=
  jsonRowParse[output, prefix]["Rows"];

jsonRealNumberQ[value_] := IntegerQ[value] || Head[value] === Real;

jsonNumericValueQ[value_] := jsonRealNumberQ[value] ||
  (AssociationQ[value] &&
    Sort[Keys[value]] === Sort[{"Re", "Im"}] &&
    jsonRealNumberQ[value["Re"]] && jsonRealNumberQ[value["Im"]]);

jsonNumericZeroQ[value_] := jsonNumericValueQ[value] && Which[
  jsonRealNumberQ[value], TrueQ[value == 0],
  AssociationQ[value],
    TrueQ[value["Re"] == 0] && TrueQ[value["Im"] == 0],
  True, False];

outputCertificationQ[record_] := AssociationQ[record] && Switch[
  Lookup[record, "Applicability", None],
  "applicable",
    Sort[Keys[record]] === Sort[{"Applicability", "Operation", "Scope",
        "ErrorGuarantee", "ErrorEnvelope"}] &&
      Lookup[record, "Operation", None] === "integrate" && Which[
        Lookup[record, "Scope", None] ===
            "full_local_with_certified_tail",
          Lookup[record, "ErrorGuarantee", None] === "certified" &&
            AssociationQ[Lookup[record, "ErrorEnvelope", None]] &&
            Lookup[record["ErrorEnvelope"], "guarantee", None] ===
              "certified",
        Lookup[record, "Scope", None] === "stored_truncation",
          Lookup[record, "ErrorGuarantee", None] === "none" &&
            Lookup[record, "ErrorEnvelope", None] === Null,
        True, False],
  "not-applicable",
    Sort[Keys[record]] ===
        Sort[{"Applicability", "Operation", "Reason"}] &&
      StringQ[Lookup[record, "Operation", None]] &&
      StringQ[Lookup[record, "Reason", None]],
  _, False];

validateCustomFinalPayload[plan_Association, row_] := Module[
  {epsilonOrder, mode, expectedKeys, rawMin, coefficients, finite,
   coefficientPayloadQ, finiteConsistencyQ},
  If[!AssociationQ[row],
    Return[failure["custom FINAL row must be a JSON object",
      <|"Row" -> row|>], Module]];
  epsilonOrder = Lookup[Lookup[plan, "Settings", <||>],
    "EpsilonOrder", Missing["Absent"]];
  If[!IntegerQ[epsilonOrder] || epsilonOrder < 0,
    Return[failure[
      "family plan has no valid EpsilonOrder for FINAL validation",
      <|"EpsilonOrder" -> epsilonOrder|>], Module]];
  mode = Lookup[Lookup[plan, "Request", <||>], "OutputMode", None];
  expectedKeys = {
    "Example", "Finite", "RawMinPower", "Certification", "Master",
    "PipelineRequestID", "FamilyID", "RequestID", "RequestOrdinal",
    "PhysicalIntegralID"};
  If[mode === "AllPendingDiscovery",
    expectedKeys = Join[expectedKeys,
      {"SelectionRequestID", "ResolutionID"}]];
  If[epsilonOrder > 0, AppendTo[expectedKeys, "Coefficients"]];
  If[Sort[Keys[row]] =!= Sort[expectedKeys],
    Return[failure["custom FINAL row has the wrong payload fields",
      <|"ExpectedKeys" -> Sort[expectedKeys],
        "ActualKeys" -> Sort[Keys[row]]|>], Module]];
  rawMin = row["RawMinPower"];
  finite = row["Finite"];
  If[Lookup[row, "Example", None] =!= Lookup[plan, "FamilyName", None] ||
      !IntegerQ[rawMin] ||
      !jsonNumericValueQ[finite] ||
      !TrueQ[outputCertificationQ[row["Certification"]]],
    Return[failure[
      "custom FINAL row has malformed numerical or certification payload",
      <|"ExpectedExample" -> Lookup[plan, "FamilyName", None],
        "Row" -> row|>], Module]];
  If[epsilonOrder === 0,
    If[rawMin > 0 && !jsonNumericZeroQ[finite],
      Return[failure[
        "custom FINAL finite coefficient conflicts with its positive minimum power",
        <|"RawMinPower" -> rawMin, "Finite" -> finite|>], Module]];
    Return[True, Module]];
  coefficients = row["Coefficients"];
  finiteConsistencyQ = If[rawMin <= 0,
    ListQ[coefficients] && Length[coefficients] >= 1 - rawMin &&
      SameQ[coefficients[[1 - rawMin, 2]], finite],
    jsonNumericZeroQ[finite]];
  coefficientPayloadQ = ListQ[coefficients] &&
    Length[coefficients] === Max[0, epsilonOrder - rawMin + 1] &&
    And @@ MapIndexed[
      Function[{coefficient, index},
        ListQ[coefficient] && Length[coefficient] === 2 &&
          coefficient[[1]] === rawMin + First[index] - 1 &&
          jsonNumericValueQ[coefficient[[2]]]],
      coefficients] &&
    TrueQ[finiteConsistencyQ];
  If[!TrueQ[coefficientPayloadQ],
    Return[failure[
      "custom FINAL coefficients do not cover the requested epsilon window",
      <|"RawMinPower" -> rawMin, "EpsilonOrder" -> epsilonOrder,
        "Coefficients" -> coefficients|>], Module]];
  True
];

validateCustomFinalPayloads[plan_Association, finals_List] := Module[
  {validations},
  validations = validateCustomFinalPayload[plan, #] & /@ finals;
  If[AnyTrue[validations, FailureQ],
    failure["child FINAL rows have malformed result payloads",
      <|"PayloadFailures" -> Select[validations, FailureQ]|>],
    True]
];

validateExplicitFamilyFinalRows[plan_Association, finals_List] := Module[
  {requests, pipelineRequestID, familyID, n, payloadValidation},
  If[Lookup[plan, "InputKind", "Registry"] =!= "Family", Return[True, Module]];
  requests = Lookup[plan, "OutputRequests", None];
  pipelineRequestID = Lookup[plan, "RequestID", None];
  familyID = Lookup[plan, "FamilyID", None];
  If[!ListQ[requests],
    Return[failure["family plan has no ordered output-request contract"],
      Module]];
  n = Length[requests];
  If[Length[finals] =!= n,
    Return[failure["child FINAL row count does not match the family request",
      <|"ExpectedCount" -> n, "ActualCount" -> Length[finals]|>], Module]];
  payloadValidation = validateCustomFinalPayloads[plan, finals];
  If[FailureQ[payloadValidation], Return[payloadValidation, Module]];
  If[!And @@ MapThread[
      Function[{row, request, ordinal},
        AssociationQ[row] && AssociationQ[request] &&
          Lookup[row, "RequestOrdinal", None] === ordinal &&
          Lookup[request, "RequestOrdinal", None] === ordinal &&
          Lookup[row, "RequestID", None] ===
            Lookup[request, "RequestID", Missing["Absent"]] &&
          Lookup[row, "PhysicalIntegralID", None] ===
            Lookup[request, "PhysicalIntegralID", Missing["Absent"]] &&
          Lookup[row, "Master", None] ===
            Lookup[request, "IndexVector", Missing["Absent"]] &&
          Lookup[row, "PipelineRequestID", None] === pipelineRequestID &&
          Lookup[row, "FamilyID", None] === familyID],
      {finals, requests, Range[n]}],
    Return[failure[
      "child FINAL rows are malformed, reordered, or belong to another request",
      <|"ExpectedOutputRequests" -> requests,
        "ActualFinalRows" -> finals|>], Module]];
  True
];

validateFamilyResolutionRows[plan_Association, rows_List] := Module[
  {request, mode, resolution},
  If[Lookup[plan, "InputKind", "Registry"] =!= "Family",
    Return[If[rows === {}, True,
      failure["registry child emitted an unexpected custom output resolution",
        <|"ResolutionRows" -> rows|>]], Module]];
  request = Lookup[plan, "Request", None];
  mode = Lookup[request, "OutputMode", None];
  Switch[mode,
    "Explicit",
      If[rows === {}, True,
        failure["explicit-target child emitted an unexpected All-master resolution",
          <|"ResolutionRows" -> rows|>]],
    "AllPendingDiscovery",
      If[Length[rows] =!= 1 || !AssociationQ[First[rows]],
        Return[failure[
          "All-master child must emit exactly one output-resolution manifest",
          <|"ResolutionRows" -> rows|>], Module]];
      resolution = First[rows];
      If[TrueQ[
          FeynmanTrick`PipelineRequest`ResolvedAllOutputSelectionQ[
            request, resolution]], resolution,
        failure[
          "All-master child output-resolution manifest is stale or malformed",
          <|"Resolution" -> resolution|>]],
    _, failure["family plan has an unsupported output-selection mode",
      <|"OutputMode" -> mode|>]
  ]
];

validateAllFamilyFinalRows[plan_Association, finals_List,
    declaredResolution_Association] := Module[
  {request, resolution, requests, masters, pipelineRequestID, familyID, n,
   payloadValidation},
  request = Lookup[plan, "Request", None];
  pipelineRequestID = Lookup[plan, "RequestID", None];
  familyID = Lookup[plan, "FamilyID", None];
  If[finals === {} || !AllTrue[finals, AssociationQ],
    Return[failure[
      "All-master child produced no complete ordered FINAL contract"], Module]];
  payloadValidation = validateCustomFinalPayloads[plan, finals];
  If[FailureQ[payloadValidation], Return[payloadValidation, Module]];
  masters = Lookup[finals, "Master", Missing["Absent"]];
  resolution =
    FeynmanTrick`PipelineRequest`CreateResolvedAllOutputSelection[
      request, masters];
  If[FailureQ[resolution],
    Return[failure[
      "All-master FINAL rows do not form a supported canonical resolution",
      <|"ResolutionFailure" -> resolution,
        "ActualFinalRows" -> finals|>], Module]];
  If[resolution =!= declaredResolution,
    Return[failure[
      "All-master FINAL rows do not match the child discovery manifest",
      <|"DeclaredResolution" -> declaredResolution,
        "DerivedResolution" -> resolution|>], Module]];
  (* Creation accepts raw FIRE order and duplicates so the runner can
     canonicalize discovery.  FINAL is stricter: it must already be exactly
     the sorted, duplicate-free order recorded by that resolution. *)
  If[masters =!= resolution["Masters"],
    Return[failure[
      "All-master FINAL rows are reordered or duplicated",
      <|"ExpectedMasters" -> resolution["Masters"],
        "ActualMasters" -> masters|>], Module]];
  requests = resolution["OutputRequests"];
  n = Length[requests];
  If[Length[finals] =!= n || !And @@ MapThread[
      Function[{row, outputRequest, ordinal},
        Lookup[row, "RequestOrdinal", None] === ordinal &&
          Lookup[outputRequest, "RequestOrdinal", None] === ordinal &&
          Lookup[row, "RequestID", None] === outputRequest["RequestID"] &&
          Lookup[row, "PhysicalIntegralID", None] ===
            outputRequest["PhysicalIntegralID"] &&
          Lookup[row, "Master", None] === outputRequest["IndexVector"] &&
          Lookup[row, "PipelineRequestID", None] === pipelineRequestID &&
          Lookup[row, "FamilyID", None] === familyID &&
          Lookup[row, "SelectionRequestID", None] ===
            resolution["SelectionRequestID"] &&
          Lookup[row, "ResolutionID", None] === resolution["ResolutionID"]],
      {finals, requests, Range[n]}],
    Return[failure[
      "All-master FINAL rows are malformed or resolution-identity mismatched",
      <|"ExpectedResolution" -> resolution,
        "ActualFinalRows" -> finals|>], Module]];
  resolution
];

validateFamilyFinalRows[plan_Association, finals_List,
    resolutionValidation_:True] := Module[{mode},
  If[Lookup[plan, "InputKind", "Registry"] =!= "Family", Return[True, Module]];
  mode = Lookup[Lookup[plan, "Request", <||>], "OutputMode", None];
  Switch[mode,
    "Explicit", validateExplicitFamilyFinalRows[plan, finals],
    "AllPendingDiscovery",
      If[AssociationQ[resolutionValidation],
        validateAllFamilyFinalRows[plan, finals, resolutionValidation],
        failure["All-master FINAL validation has no valid discovery manifest"]],
    _, failure["family plan has an unsupported output-selection mode",
      <|"OutputMode" -> mode|>]
  ]
];

pipelineResult[plan_Association, process_Association] := Module[
  {stdout, stderr, exit, finals, steps, stopped, status, outputValidation,
   validatedFinals, outputResolution, resolvedOutputRequests,
   resolutionRows, resolutionValidation, finalParse, resolutionParse,
   familyQ, finalParseValidation},
  stdout = Lookup[process, "StandardOutput", ""];
  stderr = Lookup[process, "StandardError", ""];
  exit = Lookup[process, "ExitCode", Missing["NotAvailable"]];
  steps = jsonRows[stdout, "STEPWISE "];
  finalParse = jsonRowParse[stdout, "FINAL "];
  resolutionParse = jsonRowParse[stdout, "OUTPUT_RESOLUTION "];
  finals = finalParse["Rows"];
  resolutionRows = resolutionParse["Rows"];
  stopped = StringContainsQ[stdout, "STOPPED_AFTER_BOUNDARY_LEVEL"];
  familyQ = Lookup[plan, "InputKind", "Registry"] === "Family";
  (* A successful custom child owes exactly the ordered FINAL contract.
     Every tagged row must parse: discarding a malformed extra line would let
     a partial or ambiguous child contract appear successful.  Registry runs
     retain their historical permissive parsing. *)
  finalParseValidation = If[familyQ && finalParse["Errors"] =!= {},
    failure["child emitted malformed FINAL JSON",
      <|"MalformedLines" -> finalParse["Errors"]|>], True];
  resolutionValidation = If[
    familyQ && resolutionParse["Errors"] =!= {},
    failure["child emitted malformed OUTPUT_RESOLUTION JSON",
      <|"MalformedLines" -> resolutionParse["Errors"]|>],
    validateFamilyResolutionRows[plan, resolutionRows]];
  outputValidation = Which[
    FailureQ[finalParseValidation], finalParseValidation,
    FailureQ[resolutionValidation], resolutionValidation,
    stopped && finals === {}, True,
    True, validateFamilyFinalRows[plan, finals, resolutionValidation]
  ];
  validatedFinals = If[FailureQ[outputValidation], {}, finals];
  outputResolution = If[AssociationQ[resolutionValidation] &&
      Lookup[resolutionValidation, "Schema", None] ===
        "FeynmanTrick.ResolvedAllOutputSelection/v1",
    resolutionValidation, Missing["NotApplicable"]];
  resolvedOutputRequests = Which[
    AssociationQ[outputResolution], outputResolution["OutputRequests"],
    Lookup[plan, "InputKind", "Registry"] === "Family" &&
      Lookup[Lookup[plan, "Request", <||>], "OutputMode", None] === "Explicit",
        Lookup[plan, "OutputRequests", {}],
    True, {}
  ];
  status = Which[
    exit =!= 0, "Failed",
    FailureQ[outputValidation], "Failed",
    stopped, "Stopped",
    validatedFinals =!= {}, "Succeeded",
    True, "Incomplete"];
  <|
    "Schema" -> $pipelineResultSchema, "Status" -> status,
    "ExitCode" -> exit,
    "Outputs" -> validatedFinals,
    "Final" -> Which[
      FailureQ[outputValidation], Missing["InvalidOutputs"],
      validatedFinals === {}, Missing["NotProduced"],
      Length[validatedFinals] === 1, First[validatedFinals],
      True, Missing["MultipleOutputs", Length[validatedFinals]]],
    "OutputValidation" -> outputValidation,
    "OutputResolutionValidation" -> resolutionValidation,
    "OutputResolution" -> outputResolution,
    "ResolvedOutputRequests" -> resolvedOutputRequests,
    "UnvalidatedOutputs" -> If[FailureQ[outputValidation], finals, {}],
    "UnvalidatedOutputResolutions" ->
      If[FailureQ[resolutionValidation], resolutionRows, {}],
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

preparePlanRequest[plan_Association] := Module[
  {request, requestFile, expectedFile, written, policy, family,
   executionName},
  If[Lookup[plan, "InputKind", "Registry"] =!= "Family", Return[plan, Module]];
  request = Lookup[plan, "Request", None];
  requestFile = Lookup[plan, "RequestFile", None];
  policy = Lookup[plan, "ExecutionPolicy", <||>];
  If[!TrueQ[Lookup[policy, "ExecutionReady", False]],
    Return[failure[
      If[Lookup[request, "OutputMode", None] === "AllPendingDiscovery",
        "All-integral discovery requires the built-in runner or an external Runner declared request-aware",
        "custom-family execution through an external Runner requires RequestAwareRunner -> True"],
      <|"ExecutionPolicy" -> policy,
        "RequestID" -> Lookup[request, "RequestID", None]|>], Module]];
  If[!TrueQ[FeynmanTrick`PipelineRequest`PipelineRequestQ[request]] ||
      !StringQ[requestFile],
    Return[failure["custom-family PipelinePlan request metadata is malformed"],
      Module]];
  family = request["Family"];
  executionName = "family_" <> StringTake[request["RequestID"], -16];
  If[Lookup[plan, "RequestID", None] =!= request["RequestID"] ||
      Lookup[plan, "FamilyID", None] =!= family["FamilyID"] ||
      Lookup[plan, "FamilyName", None] =!= family["Name"] ||
      Lookup[plan, "OutputRequests", None] =!= request["OutputRequests"] ||
      Lookup[plan, "Example", None] =!= executionName ||
      Lookup[plan["Environment"], "FT_EXAMPLES", None] =!= executionName,
    Return[failure[
      "custom-family PipelinePlan duplicates do not match its validated request",
      <|"ExpectedRequestID" -> request["RequestID"],
        "ExpectedFamilyID" -> family["FamilyID"],
        "ExpectedExecutionName" -> executionName|>], Module]];
  expectedFile = FeynmanTrick`PipelineRequest`PipelineRequestPath[
    request, Lookup[plan, "RequestDirectory", ""]];
  If[FailureQ[expectedFile] || ExpandFileName[requestFile] =!= expectedFile ||
      Lookup[plan["Environment"], "FT_FAMILY_REQUEST_ID", None] =!=
        request["RequestID"] ||
      Lookup[plan["Environment"], "FT_FAMILY_REQUEST_FILE", None] =!=
        expectedFile,
    Return[failure[
      "custom-family PipelinePlan request path or identity was modified",
      <|"ExpectedRequestFile" -> expectedFile,
        "RequestFile" -> requestFile|>], Module]];
  written = FeynmanTrick`PipelineRequest`WritePipelineRequest[
    request, expectedFile];
  If[FailureQ[written], written, plan]
];

runPlan[plan_, asynchronous_] := Module[{preparedPlan, process, command},
  If[FailureQ[plan], Return[plan, Module]];
  preparedPlan = preparePlanRequest[plan];
  If[FailureQ[preparedPlan], Return[preparedPlan, Module]];
  command = processCommand[preparedPlan];
  If[FailureQ[command], Return[command, Module]];
  If[TrueQ[asynchronous],
    process = Quiet[Check[StartProcess[command,
      ProcessDirectory -> preparedPlan["WorkingDirectory"]], $Failed]];
    If[Head[process] =!= ProcessObject,
      Return[failure["could not start the pipeline subprocess",
        <|"Plan" -> preparedPlan|>], Module]];
    <|"Schema" -> $pipelineProcessSchema, "Status" -> "Running",
      "Process" -> process, "Plan" -> preparedPlan|>,
    process = Quiet[Check[RunProcess[command, All,
      ProcessDirectory -> preparedPlan["WorkingDirectory"]], $Failed]];
    If[!AssociationQ[process],
      Return[failure["pipeline subprocess could not be executed",
        <|"Plan" -> preparedPlan|>], Module]];
    pipelineResult[preparedPlan, process]]];

Options[RunIntegrationPipeline] = Options[PipelinePlan];
RunIntegrationPipeline[example_String, opts:OptionsPattern[]] := Module[{plan},
  plan = PipelinePlan[example, opts];
  runPlan[plan, OptionValue["Asynchronous"]]];

RunIntegrationPipeline[family_Association, targets:(All | _List),
    opts:OptionsPattern[]] := Module[{plan},
  plan = PipelinePlan[family, targets, opts];
  runPlan[plan, OptionValue["Asynchronous"]]];

RunIntegrationPipeline[plan_Association, opts:OptionsPattern[]] := Module[
  {schema, asynchronous, required, missing, optionRules = {opts}, familyPlan},
  schema = Lookup[plan, "Schema", Missing["NotAvailable"]];
  If[schema === "FeynmanTrick.FamilySpec/v1",
    familyPlan = PipelinePlan[plan, opts];
    Return[runPlan[familyPlan, OptionValue["Asynchronous"]], Module]];
  If[optionRules =!= {},
    Return[failure[
      "options cannot be applied while executing an existing PipelinePlan; rebuild the plan instead",
      <|"Options" -> optionRules|>], Module]];
  If[schema =!= $pipelinePlanSchema,
    Return[failure["RunIntegrationPipeline requires a PipelinePlan/v1 record",
      <|"Schema" -> schema|>], Module]];
  required = {"Command", "WorkingDirectory", "Environment", "Settings"};
  missing = Select[required, !KeyExistsQ[plan, #] &];
  If[missing =!= {} || !MatchQ[plan["Command"], {_String ..}] ||
      !StringQ[plan["WorkingDirectory"]] ||
      !AssociationQ[plan["Environment"]] || !AssociationQ[plan["Settings"]],
    Return[failure["PipelinePlan record is malformed",
      <|"MissingKeys" -> missing|>], Module]];
  asynchronous = TrueQ[Lookup[Lookup[plan, "Settings", <||>],
    "Asynchronous", False]];
  runPlan[plan, asynchronous]];

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
