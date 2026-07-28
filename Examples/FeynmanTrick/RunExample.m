(* Public command-line launcher for the documented Feynman-trick examples.
   Numerical execution goes through FeynmanTrick`PipelinePlan and
   FeynmanTrick`RunIntegrationPipeline; the internal ladder driver is not an
   example-facing entry point. *)

repoRoot = ParentDirectory[
  ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

profiles = <|
  "bubble" -> <|
    "CheckpointSlug" -> "bubble",
    "Options" -> {
      "WorkingPrecision" -> 300,
      "ExpansionOrder" -> 40,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 10,
      "LevelEpsilonHalos" -> {0},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>,
  "sunrise" -> <|
    "CheckpointSlug" -> "sunrise",
    "Options" -> {
      "WorkingPrecision" -> 300,
      "ExpansionOrder" -> 50,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 10,
      "LevelEpsilonHalos" -> {0, 4},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>,
  "banana_unequal" -> <|
    "CheckpointSlug" -> "banana-unequal",
    "Options" -> {
      "WorkingPrecision" -> 500,
      "ExpansionOrder" -> 70,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 10,
      "LevelEpsilonHalos" -> {0, 4, 7},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>,
  "box_bubble" -> <|
    "CheckpointSlug" -> "box-bubble",
    "Options" -> {
      "WorkingPrecision" -> 300,
      "ExpansionOrder" -> 40,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 16,
      "LevelEpsilonHalos" -> {0},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>,
  "kite" -> <|
    "CheckpointSlug" -> "kite",
    "Options" -> {
      "WorkingPrecision" -> 300,
      "ExpansionOrder" -> 40,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 16,
      "LevelEpsilonHalos" -> {0, 4, 7, 7},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>,
  "banana4_unequal" -> <|
    "CheckpointSlug" -> "banana4-unequal",
    "Options" -> {
      "WorkingPrecision" -> 500,
      "ExpansionOrder" -> 70,
      "EpsilonOrder" -> 0,
      "BoundaryExtraOrder" -> 10,
      "LevelEpsilonHalos" -> {0, 4, 7, 7},
      "DivisionOrder" -> 3,
      "RadiusOfConvergence" -> 1}|>
|>;

environmentOrDefault[name_String, default_String] := Module[
  {value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0,
    StringTrim[value], default]];

positiveIntegerEnvironment[name_String, default_Integer] := Module[
  {text = environmentOrDefault[name, ToString[default]], value},
  If[!StringMatchQ[text, RegularExpression["[0-9]+"]],
    Return[Failure["ExampleConfiguration",
      <|"Detail" -> "environment value must be a positive integer",
        "Name" -> name, "Value" -> text|>], Module]];
  value = ToExpression[text];
  If[IntegerQ[value] && value >= 1, value,
    Failure["ExampleConfiguration",
      <|"Detail" -> "environment value must be a positive integer",
        "Name" -> name, "Value" -> text|>]]];

flagEnvironment[name_String, default_] := Switch[
  environmentOrDefault[name, If[TrueQ[default], "1", "0"]],
  "0", False,
  "1", True,
  value_, Failure["ExampleConfiguration",
    <|"Detail" -> "environment flag must be 0 or 1",
      "Name" -> name, "Value" -> value|>]];

buildExamplePlan[name_String] := Module[
  {profile, cacheRoot, threads, rebuild, firePath, commonOptions,
   extraEnvironment = <||>, plan},
  If[!KeyExistsQ[profiles, name],
    Return[Failure["ExampleConfiguration",
      <|"Detail" -> "unknown documented Feynman-trick example",
        "Example" -> name, "SupportedExamples" -> Keys[profiles]|>], Module]];
  profile = profiles[name];
  cacheRoot = ExpandFileName[environmentOrDefault[
    "DIFFEXP2_CACHE_DIR", FileNameJoin[{$HomeDirectory, ".cache", "diffexp2"}]]];
  threads = positiveIntegerEnvironment["DE2_CPP_THREADS", 4];
  rebuild = flagEnvironment["FT_REBUILD_PREP", False];
  If[FailureQ[threads], Return[threads, Module]];
  If[FailureQ[rebuild], Return[rebuild, Module]];
  firePath = environmentOrDefault["FT_FIRE_PATH", ""];
  commonOptions = {
    "RecurrenceBackend" -> "Cpp",
    "CppThreads" -> threads,
    "ValueTransport" -> True,
    "NativeValueHopExecution" -> True,
    "BatchEndpointArms" -> True,
    "PreparedCacheDirectory" -> FileNameJoin[{cacheRoot, "fire"}],
    "CheckpointDirectory" -> FileNameJoin[{
      cacheRoot, "ladder", profile["CheckpointSlug"]}],
    "RebuildPreparation" -> rebuild,
    "WorkingDirectory" -> repoRoot,
    "EchoOutput" -> True};
  Scan[Function[variable,
    With[{value = Environment[variable]},
      If[StringQ[value] && StringLength[StringTrim[value]] > 0,
        AssociateTo[extraEnvironment,
          variable -> StringTrim[value]]]]],
    {"DE2_NATIVE_STAGE_TIMING"}];
  If[extraEnvironment =!= <||>,
    commonOptions = Append[commonOptions,
      "ExtraEnvironment" -> extraEnvironment]];
  If[firePath =!= "",
    commonOptions = Append[commonOptions, "FIREPath" -> firePath]];
  plan = FeynmanTrick`PipelinePlan[name,
    Sequence @@ Join[profile["Options"], commonOptions]];
  plan
];

(* wolframscript -script keeps this filename as the first element. *)
arguments = Rest[$ScriptCommandLine];
If[arguments === {"--check"},
  plans = AssociationMap[buildExamplePlan, Keys[profiles]];
  failedPlans = Select[plans, FailureQ];
  If[Length[failedPlans] > 0,
    Print["Feynman-trick example plan check failed: ",
      InputForm[failedPlans]];
    Exit[1]];
  Print["Feynman-trick example plans valid: ",
    StringRiffle[Keys[plans], ", "]];
  Exit[0]];

planOnly = Length[arguments] === 2 && First[arguments] === "--plan";
If[planOnly, arguments = Rest[arguments]];
If[Length[arguments] =!= 1,
  Print["Usage: wolframscript -script ", $InputFileName,
    " [--plan] <example>"];
  Print["       wolframscript -script ", $InputFileName, " --check"];
  Exit[2]];

example = First[arguments];
plan = buildExamplePlan[example];
If[FailureQ[plan],
  Print["Could not build Feynman-trick example plan: ", InputForm[plan]];
  Exit[2]];

If[planOnly,
  Print[InputForm[plan]];
  Exit[0]];

If[!DirectoryQ[plan["PreparedCacheDirectory"]],
  CreateDirectory[plan["PreparedCacheDirectory"],
    CreateIntermediateDirectories -> True]];
If[!DirectoryQ[plan["CheckpointDirectory"]],
  CreateDirectory[plan["CheckpointDirectory"],
    CreateIntermediateDirectories -> True]];
result = FeynmanTrick`RunIntegrationPipeline[plan];
If[FailureQ[result],
  Print["Feynman-trick pipeline failed before producing a result: ",
    InputForm[result]];
  Exit[1]];

Print["Feynman-trick pipeline status: ", result["Status"],
  " (exit ", result["ExitCode"], ")"];
If[MemberQ[{"Succeeded", "Stopped"}, result["Status"]],
  Exit[0]];
exitCode = result["ExitCode"];
Exit[If[IntegerQ[exitCode] && exitCode =!= 0, exitCode, 1]];
