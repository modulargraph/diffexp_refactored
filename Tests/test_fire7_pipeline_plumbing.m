(* Process-free contract for FIRE7 option, environment, and runner plumbing. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

fireEnvironmentNames = {
  "FT_FIRE_PATH", "FT_FIRE_BACKEND", "FT_FIRE_CALC",
  "FT_FIRE_MODULAR_WORKERS", "FT_FIRE_USE_MULTIPRIME",
  "FT_FIRE_PRIME_LIMIT", "FT_FIRE_KEEP_MODULAR_TABLES",
  "FT_FIRE_DIMENSION_SEPARATED", "FT_FIRE_MULTIPRIME_WIDTH",
  "FT_FIRE_MPI_EXECUTABLE", "FT_FIRE_BASIS_PROBE_COUNT",
  "FT_FIRE_MODULAR_CACHE_DIR", "FT_RUNNER_DEFINITIONS_ONLY"
};
clearFireEnvironment[] := Scan[SetEnvironment[# -> None] &,
  fireEnvironmentNames];
clearFireEnvironment[];

defaultPlan = FeynmanTrick`PipelinePlan["bubble"];
defaultFirePath = ExpandFileName[FileNameJoin[{
  repoRoot, "Dependencies", "fire", "FIRE7", "FIRE7"}]];
defaultModularCache = ExpandFileName[FileNameJoin[{
  $TemporaryDirectory, "DiffExp2_FIRE7_Modular"}]];
defaultWorkers = Max[1, Min[10, $ProcessorCount]];

assert["FIRE7 modular defaults are explicit in the process-free plan",
  AssociationQ[defaultPlan] &&
    defaultPlan["Settings", "FIREPath"] === defaultFirePath &&
    defaultPlan["Settings", "FIREBackend"] === "Modular" &&
    defaultPlan["Settings", "FIRECalc"] === "flint" &&
    defaultPlan["Settings", "FIREModularWorkers"] === defaultWorkers &&
    TrueQ[defaultPlan["Settings", "FIREUseMultiprime"]] &&
    defaultPlan["Settings", "FIREPrimeLimit"] === 127 &&
    TrueQ[defaultPlan["Settings", "FIREKeepModularTables"]] &&
    defaultPlan["Settings", "FIREDimensionSeparated"] === False &&
    defaultPlan["Settings", "FIREMultiprimeWidth"] === 16 &&
    defaultPlan["Settings", "FIREMPIExecutable"] === Automatic &&
    defaultPlan["Settings", "FIREBasisProbeCount"] === 2 &&
    defaultPlan["Settings", "FIREModularCacheDirectory"] ===
      defaultModularCache];

assert["package configuration agrees with the facade FIRE7 defaults",
  With[{configuration = FeynmanTrick`FTConfiguration[]},
    configuration["FIREPath"] === defaultFirePath &&
      configuration["FIREBackend"] === "Modular" &&
      configuration["FIRECalc"] === "flint" &&
      configuration["FIREModularWorkers"] === defaultWorkers &&
      TrueQ[configuration["FIREUseMultiprime"]] &&
      TrueQ[configuration["FIREKeepModularTables"]] &&
      configuration["FIREModularCacheDirectory"] === defaultModularCache]];

customCache = ExpandFileName[FileNameJoin[{
  $TemporaryDirectory, "fire7-pipeline-plumbing-custom"}]];
customPlan = FeynmanTrick`PipelinePlan["bubble",
  "FIREPath" -> FileNameJoin[{$TemporaryDirectory, "fire7-custom"}],
  "FIREBackend" -> "Classical",
  "FIRECalc" -> "fermat",
  "FIREModularWorkers" -> 3,
  "FIREUseMultiprime" -> False,
  "FIREPrimeLimit" -> 19,
  "FIREKeepModularTables" -> False,
  "FIREDimensionSeparated" -> True,
  "FIREMultiprimeWidth" -> 7,
  "FIREMPIExecutable" -> "custom-mpirun",
  "FIREBasisProbeCount" -> 4,
  "FIREModularCacheDirectory" -> customCache,
  "ExtraEnvironment" -> <|
    "FT_FIRE_BACKEND" -> "stale-backend",
    "FT_FIRE_CALC" -> "stale-calc",
    "FT_FIRE_MODULAR_WORKERS" -> "999",
    "DEBUG_FIRE7_PLUMBING" -> "1"|>];

assert["all FIRE7 options serialize to canonical runner environment fields",
  AssociationQ[customPlan] &&
    customPlan["Environment", "FT_FIRE_BACKEND"] === "Classical" &&
    customPlan["Environment", "FT_FIRE_CALC"] === "fermat" &&
    customPlan["Environment", "FT_FIRE_MODULAR_WORKERS"] === "3" &&
    customPlan["Environment", "FT_FIRE_USE_MULTIPRIME"] === "0" &&
    customPlan["Environment", "FT_FIRE_PRIME_LIMIT"] === "19" &&
    customPlan["Environment", "FT_FIRE_KEEP_MODULAR_TABLES"] === "0" &&
    customPlan["Environment", "FT_FIRE_DIMENSION_SEPARATED"] === "1" &&
    customPlan["Environment", "FT_FIRE_MULTIPRIME_WIDTH"] === "7" &&
    customPlan["Environment", "FT_FIRE_MPI_EXECUTABLE"] ===
      "custom-mpirun" &&
    customPlan["Environment", "FT_FIRE_BASIS_PROBE_COUNT"] === "4" &&
    customPlan["Environment", "FT_FIRE_MODULAR_CACHE_DIR"] === customCache];

assert["canonical FIRE7 fields override ExtraEnvironment while additions survive",
  customPlan["Environment", "FT_FIRE_BACKEND"] === "Classical" &&
    customPlan["Environment", "FT_FIRE_CALC"] === "fermat" &&
    customPlan["Environment", "FT_FIRE_MODULAR_WORKERS"] === "3" &&
    customPlan["Environment", "DEBUG_FIRE7_PLUMBING"] === "1"];

invalidPlans = {
  FeynmanTrick`PipelinePlan["bubble", "FIREBackend" -> "Unknown"],
  FeynmanTrick`PipelinePlan["bubble", "FIRECalc" -> "bad calc"],
  FeynmanTrick`PipelinePlan["bubble", "FIREModularWorkers" -> 0],
  FeynmanTrick`PipelinePlan["bubble", "FIREModularWorkers" -> 11],
  FeynmanTrick`PipelinePlan["bubble", "FIREUseMultiprime" -> 1],
  FeynmanTrick`PipelinePlan["bubble", "FIREPrimeLimit" -> 0],
  FeynmanTrick`PipelinePlan["bubble", "FIREPrimeLimit" -> 128],
  FeynmanTrick`PipelinePlan["bubble", "FIREKeepModularTables" -> 1],
  FeynmanTrick`PipelinePlan["bubble", "FIREDimensionSeparated" -> 1],
  FeynmanTrick`PipelinePlan["bubble", "FIREMultiprimeWidth" -> 0],
  FeynmanTrick`PipelinePlan["bubble", "FIREMPIExecutable" -> ""],
  FeynmanTrick`PipelinePlan["bubble", "FIREBasisProbeCount" -> 1],
  FeynmanTrick`PipelinePlan["bubble", "FIREBasisProbeCount" -> 128],
  FeynmanTrick`PipelinePlan["bubble", "FIREModularCacheDirectory" -> ""]
};
assert["invalid FIRE7 plan options fail before process construction",
  AllTrue[invalidPlans, FailureQ]];

runnerEnvironment = <|
  "FT_FIRE_PATH" -> FileNameJoin[{$TemporaryDirectory, "runner-fire7"}],
  "FT_FIRE_BACKEND" -> "Classical",
  "FT_FIRE_CALC" -> "fermat",
  "FT_FIRE_MODULAR_WORKERS" -> "5",
  "FT_FIRE_USE_MULTIPRIME" -> "0",
  "FT_FIRE_PRIME_LIMIT" -> "23",
  "FT_FIRE_KEEP_MODULAR_TABLES" -> "0",
  "FT_FIRE_DIMENSION_SEPARATED" -> "1",
  "FT_FIRE_MULTIPRIME_WIDTH" -> "9",
  "FT_FIRE_MPI_EXECUTABLE" -> "runner-mpirun",
  "FT_FIRE_BASIS_PROBE_COUNT" -> "3",
  "FT_FIRE_MODULAR_CACHE_DIR" -> customCache,
  "FT_RUNNER_DEFINITIONS_ONLY" -> "1"
|>;
Scan[SetEnvironment, Normal[runnerEnvironment]];
runnerSettings =
  FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[];

assert["runner environment parser restores the complete typed FIRE7 contract",
  AssociationQ[runnerSettings] &&
    runnerSettings["FIREBackend"] === "Classical" &&
    runnerSettings["FIRECalc"] === "fermat" &&
    runnerSettings["FIREModularWorkers"] === 5 &&
    runnerSettings["FIREUseMultiprime"] === False &&
    runnerSettings["FIREPrimeLimit"] === 23 &&
    runnerSettings["FIREKeepModularTables"] === False &&
    runnerSettings["FIREDimensionSeparated"] === True &&
    runnerSettings["FIREMultiprimeWidth"] === 9 &&
    runnerSettings["FIREMPIExecutable"] === "runner-mpirun" &&
    runnerSettings["FIREBasisProbeCount"] === 3 &&
    runnerSettings["FIREModularCacheDirectory"] === customCache];

Quiet[Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]],
  {General::shdw, Symbol::shdw}];
runnerConfiguration = FeynmanTrick`FTConfiguration[];
assert["definitions-only runner installs FIRE7 settings before execution",
  And @@ (#1 === runnerConfiguration[#2] & @@@ {
    {runnerSettings["FIREPath"], "FIREPath"},
    {"Classical", "FIREBackend"}, {"fermat", "FIRECalc"},
    {5, "FIREModularWorkers"}, {False, "FIREUseMultiprime"},
    {23, "FIREPrimeLimit"}, {False, "FIREKeepModularTables"},
    {True, "FIREDimensionSeparated"}, {9, "FIREMultiprimeWidth"},
    {"runner-mpirun", "FIREMPIExecutable"}, {3, "FIREBasisProbeCount"},
    {customCache, "FIREModularCacheDirectory"}})];

unsafeEnvironmentCases = {
  {"FT_FIRE_MODULAR_WORKERS", "11"},
  {"FT_FIRE_PRIME_LIMIT", "128"},
  {"FT_FIRE_BASIS_PROBE_COUNT", "1"},
  {"FT_FIRE_BASIS_PROBE_COUNT", "128"}
};
unsafeEnvironmentRejected = True;
Do[
  clearFireEnvironment[];
  SetEnvironment[case[[1]] -> case[[2]]];
  unsafeEnvironmentRejected = unsafeEnvironmentRejected &&
    FailureQ[FeynmanTrick`DiffExp2Pipeline`RunnerSettingsFromEnvironment[]],
  {case, unsafeEnvironmentCases}];
assert["runner parser rejects unsafe modular settings and upper bounds",
  unsafeEnvironmentRejected];
clearFireEnvironment[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
