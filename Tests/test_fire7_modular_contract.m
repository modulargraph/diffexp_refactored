(* Focused FIRE 7 modular-runner and FIRE-interface contracts.

   This suite never invokes FIRE.  The only child processes are short fake
   shell launchers used to verify exact-artifact and watchdog semantics. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

makeExecutable[path_String] :=
  Lookup[RunProcess[{"/bin/chmod", "+x", path}], "ExitCode", -1] === 0;

processGoneOrZombieQ[pid_Integer] := Module[{probe, state},
  probe = Quiet[Check[
    RunProcess[{"/bin/ps", "-o", "state=", "-p", ToString[pid]}],
    <||>]];
  state = StringTrim[Lookup[probe, "StandardOutput", ""]];
  Lookup[probe, "ExitCode", 1] =!= 0 || state === "" ||
    StringStartsQ[state, "Z"]
];
processGoneOrZombieQ[_] := False;

tempRoot = FileNameJoin[{$TemporaryDirectory,
  "ft_fire7_contract_" <> ToString[$ProcessID] <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
If[DirectoryQ[tempRoot], DeleteDirectory[tempRoot, DeleteContents -> True]];
CreateDirectory[tempRoot, CreateIntermediateDirectories -> True];

(* ---------------------------------------------------------------------- *)
(* Exact config manifest and FIRE 7 output-name contracts.                *)
(* ---------------------------------------------------------------------- *)

configDir = FileNameJoin[{tempRoot, "config"}];
CreateDirectory[configDir, CreateIntermediateDirectories -> True];
Export[FileNameJoin[{configDir, "contract.config"}], StringJoin[
  "#threads 1\n",
  "#fthreads 1\n",
  "#variables x,m,d\n",
  "#start\n",
  "#problem 11 alpha.start\n",
  "#problem 12 beta.start\n",
  "#integrals contract.m\n",
  "#output contract.tables\n"], "Text"];
Export[FileNameJoin[{configDir, "contract.m"}],
  "{{11,{1,0}},{12,{0,1}}}\n", "Text"];
Export[FileNameJoin[{configDir, "alpha.start"}], "alpha start\n", "Text"];
Export[FileNameJoin[{configDir, "beta.start"}], "beta start\n", "Text"];

manifest = FeynmanTrick`FIRE7Runner`ParseConfig[configDir, "contract"];
expectedManifest = <|
  "ConfigPath" -> FileNameJoin[{configDir, "contract.config"}],
  "ConfigFile" -> "contract.config",
  "ConfigStem" -> "contract",
  "Variables" -> {"x", "m", "d"},
  "OutputFile" -> "contract.tables",
  "IntegralsFile" -> "contract.m",
  "MasterFiles" -> {},
  "Problems" -> {
    <|"ProblemNumber" -> 11, "StartFile" -> "alpha.start"|>,
    <|"ProblemNumber" -> 12, "StartFile" -> "beta.start"|>
  }
|>;
assert["ParseConfig preserves the exact d-last manifest",
  manifest === expectedManifest];
Export[FileNameJoin[{configDir, "contract-masters.config"}], StringJoin[
  "#threads 1\n",
  "#fthreads 1\n",
  "#variables x,m,d\n",
  "#start\n",
  "#problem 11 alpha.start\n",
  "#preferred contract.preferred\n",
  "#integrals contract.m\n",
  "#output contract-masters.tables\n"], "Text"];
Export[FileNameJoin[{configDir, "contract.preferred"}],
  "{{11,{1,0}},{11,{0,1}}}\n", "Text"];
masterManifest = FeynmanTrick`FIRE7Runner`ParseConfig[
  configDir, "contract-masters"];
assert["ParseConfig carries explicit master files as exact job inputs",
  AssociationQ[masterManifest] &&
    masterManifest["MasterFiles"] === {"contract.preferred"} &&
    MemberQ[
      FeynmanTrick`FIRE7Runner`Private`inputPaths[
        configDir, masterManifest],
      FileNameJoin[{configDir, "contract.preferred"}]]];
writerTopology = <|
  "Name" -> "alpha",
  "ProblemNumber" -> 11,
  "NumPropagators" -> 2,
  "Variables" -> {Global`x, Global`m},
  "Masters" -> {{1, 0}, {0, 1}}|>;
writerResult =
  FeynmanTrick`FIREInterface`Private`writeSingleProblemConfig[
    writerTopology, configDir, "writer.config",
    "contract.m", "writer.tables"];
writerManifest = FeynmanTrick`FIRE7Runner`ParseConfig[
  configDir, "writer"];
assert["single-problem reductions publish the selected basis as preferred",
  AssociationQ[writerResult] &&
    AssociationQ[writerManifest] &&
    writerManifest["MasterFiles"] === {"writer.preferred"} &&
    Get[FileNameJoin[{configDir, "writer.preferred"}]] ===
      {{11, {1, 0}}, {11, {0, 1}}}];
Export[FileNameJoin[{configDir, "unsafe.config"}], StringJoin[
  "#variables x,a/../../victim,d\n#start\n",
  "#problem 11 alpha.start\n#integrals contract.m\n",
  "#output contract.tables\n"], "Text"];
unsafeManifest = FeynmanTrick`FIRE7Runner`ParseConfig[
  configDir, "unsafe"];
Export[FileNameJoin[{configDir, "duplicate-vars.config"}], StringJoin[
  "#variables x,x,d\n#start\n",
  "#problem 11 alpha.start\n#integrals contract.m\n",
  "#output contract.tables\n"], "Text"];
duplicateVariableManifest = FeynmanTrick`FIRE7Runner`ParseConfig[
  configDir, "duplicate-vars"];
assert["ParseConfig rejects path-bearing and duplicate variable tokens",
  FailureQ[unsafeManifest] && FailureQ[duplicateVariableManifest]];
Export[FileNameJoin[{configDir, "colliding-output.config"}], StringJoin[
  "#variables x,d\n#start\n",
  "#problem 11 collision.tables\n#integrals contract.m\n",
  "#output collision.tables\n"], "Text"];
collidingOutputManifest = FeynmanTrick`FIRE7Runner`ParseConfig[
  configDir, "colliding-output"];
assert["ParseConfig rejects output/input collisions and wrong input suffixes",
  FailureQ[collidingOutputManifest]];
assert["reconstructed table name includes variables in exact d-last order",
  FeynmanTrick`FIRE7Runner`ExpectedReconstructedFile[
    configDir, manifest] ===
      FileNameJoin[{configDir, "contract_x_m_d_0.tables"}]];
assert["prime table name includes every supplied sample coordinate",
  FeynmanTrick`FIRE7Runner`ExpectedPrimeFile[
    configDir, manifest, {41, 53, 65, 1}] ===
      FileNameJoin[{configDir, "contract_41_53_65_1.tables"}]];

(* ---------------------------------------------------------------------- *)
(* Exact modular argv contract.                                           *)
(* ---------------------------------------------------------------------- *)

fireRoot = FileNameJoin[{tempRoot, "fire7"}];
fireBinDir = FileNameJoin[{fireRoot, "bin"}];
CreateDirectory[fireBinDir, CreateIntermediateDirectories -> True];
Export[FileNameJoin[{fireRoot, "FIRE7.m"}], "(* fake FIRE7 source *)\n", "Text"];
Export[FileNameJoin[{fireRoot, "paths.inc"}],
  "RESULTING_MPRIME_COUNT=16\n", "Text"];
Export[FileNameJoin[{fireBinDir, "FIRE7p"}], "fake prime worker\n", "Text"];
Export[FileNameJoin[{fireBinDir, "FIRE7mp"}],
  "fake multipoint worker\n", "Text"];
Export[FileNameJoin[{fireBinDir, "FIRE7_MPI"}], "fake mpi worker\n", "Text"];
Export[FileNameJoin[{fireBinDir, "reconstruct"}], "fake reconstruct\n", "Text"];

fakeMPI = FileNameJoin[{tempRoot, "fake_mpirun"}];
Export[fakeMPI, "#!/bin/sh\nexit 0\n", "Text"];
assert["fake MPI launcher is executable", makeExecutable[fakeMPI]];

defaultSettings = <|
  "MPIExecutable" -> fakeMPI,
  "Workers" -> 3,
  "Calc" -> "flint",
  "UseMultiprime" -> True,
  "PrimeLimit" -> 127,
  "DimensionSeparated" -> False,
  "MultiprimeWidth" -> 16,
  "TimeoutSeconds" -> 5
|>;
assert["same-job lock wait exceeds the configured FIRE job timeout",
  FeynmanTrick`FIRE7Runner`Private`lockWaitSeconds[
    defaultSettings, "Reduction"] === 70 &&
    FeynmanTrick`FIRE7Runner`Private`lockWaitSeconds[
      Join[defaultSettings, <|"BasisProbeCount" -> 3|>], "Basis"] === 75 &&
    FailureQ[FeynmanTrick`FIRE7Runner`Private`lockWaitSeconds[
      Join[defaultSettings, <|"LockTimeoutSeconds" -> 69|>],
      "Reduction"]]];

oversizedStem = "oversized";
Export[FileNameJoin[{configDir, oversizedStem <> ".config"}], StringJoin[
  "#threads 1\n#fthreads 1\n#variables ",
  StringRiffle[Join[
    ("v" <> ToString[#] & /@ Range[16]), {"d"}], ","], "\n",
  "#start\n#problem 11 alpha.start\n",
  "#integrals contract.m\n#output oversized.tables\n"], "Text"];
oversizedRun = FeynmanTrick`FIRE7Runner`RunModular[
  fireRoot, configDir, oversizedStem, defaultSettings];
assert["public modular runner rejects more than 16 variables before launch",
  FailureQ[oversizedRun] &&
    oversizedRun["Detail"] ===
      "FIRE7 modular reconstruction supports at most 16 variables including the dimension" &&
    oversizedRun["VariableCount"] === 17];

command = FeynmanTrick`FIRE7Runner`ModularCommand[
  fireRoot, manifest, defaultSettings];
expectedCommand = {
  ExpandFileName[fakeMPI], "--oversubscribe", "-np", "4",
  ExpandFileName[FileNameJoin[{fireBinDir, "FIRE7_MPI"}]],
  "--calc", "flint", "--zippel", "--reconstruct",
  "--geometric", "--early_abortion",
  "--rational_reconstruction_limit", "127",
  "--multitables", "--big_primes", "--config", "contract"
};
assert["ModularCommand uses one coordinator plus the configured workers",
  command === expectedCommand];
assert["default modular command enables reconstruction and multiprime Zippel",
  AssociationQ[manifest] &&
    And @@ (MemberQ[command, #] & /@
      {"--zippel", "--geometric", "--reconstruct", "--multitables",
        "--big_primes", "--early_abortion"})];
assert["default modular command avoids unsafe or unsupported policy flags",
  FreeQ[command, Alternatives[
    "--delete_tables", "--reserve", "--no_integrity", "--last_separated"]]];

adaptiveIdentity = FeynmanTrick`FIRE7Runner`Private`jobIdentity[
  "Reduction", fireRoot, configDir, manifest, defaultSettings];
legacyIdentityRecord = ReplacePart[adaptiveIdentity["Record"],
  "SamplingPolicy" -> "geometric-zippel/v1"];
legacyIdentityID = IntegerString[
  Hash[legacyIdentityRecord, "SHA256"], 16, 64];
assert["adaptive-limit policy is versioned into modular cache identity",
  AssociationQ[adaptiveIdentity] &&
    adaptiveIdentity["Record", "SamplingPolicy"] ===
      "geometric-zippel-adaptive-limits/v2" &&
    adaptiveIdentity["ID"] =!= legacyIdentityID];

singlePrimeCommand = FeynmanTrick`FIRE7Runner`ModularCommand[
  fireRoot, manifest, Join[defaultSettings, <|"UseMultiprime" -> False|>]];
assert["single-prime modular mode keeps independent geometric bases",
  MemberQ[singlePrimeCommand, "--big_primes"] &&
    !MemberQ[singlePrimeCommand, "--multitables"]];

(* ---------------------------------------------------------------------- *)
(* Exit zero is insufficient: only a nonempty exact reconstructed table   *)
(* may be published.                                                      *)
(* ---------------------------------------------------------------------- *)

missingMPI = FileNameJoin[{tempRoot, "fake_mpirun_missing"}];
Export[missingMPI, StringJoin[
  "#!/bin/sh\n",
  "mkdir -p artifacts\n",
  "printf 'wrong artifact\\n' > artifacts/contract.tables\n",
  "exit 0\n"], "Text"];
assert["missing-artifact fake launcher is executable",
  makeExecutable[missingMPI]];
missingSettings = Join[defaultSettings, <|
  "MPIExecutable" -> missingMPI,
  "CacheDirectory" -> FileNameJoin[{tempRoot, "cache-missing"}]
|>];
missingResult = FeynmanTrick`FIRE7Runner`RunModular[
  fireRoot, configDir, "contract", missingSettings];
assert["exit zero with only a canonical-looking decoy is rejected",
  FailureQ[missingResult] &&
    missingResult["Detail"] ===
      "FIRE7 modular reconstruction did not produce the exact final table" &&
    missingResult["ExitCode"] === 0 &&
    StringEndsQ[missingResult["ExpectedTable"],
      "contract_x_m_d_0.tables"] &&
    !FileExistsQ[FileNameJoin[{configDir, "contract.tables"}]]];

emptyMPI = FileNameJoin[{tempRoot, "fake_mpirun_empty"}];
Export[emptyMPI, StringJoin[
  "#!/bin/sh\n",
  "mkdir -p artifacts\n",
  ": > artifacts/contract_x_m_d_0.tables\n",
  "exit 0\n"], "Text"];
assert["empty-artifact fake launcher is executable", makeExecutable[emptyMPI]];
emptySettings = Join[defaultSettings, <|
  "MPIExecutable" -> emptyMPI,
  "CacheDirectory" -> FileNameJoin[{tempRoot, "cache-empty"}]
|>];
emptyResult = FeynmanTrick`FIRE7Runner`RunModular[
  fireRoot, configDir, "contract", emptySettings];
assert["exit zero with an empty exact-name artifact is rejected",
  FailureQ[emptyResult] &&
    emptyResult["Detail"] ===
      "FIRE7 modular reconstruction did not produce the exact final table" &&
    emptyResult["ExitCode"] === 0 &&
    !FileExistsQ[emptyResult["ExpectedTable"]] &&
    Length[FileNames[FileNameTake[emptyResult["ExpectedTable"]] <>
      ".empty-*", DirectoryName[emptyResult["ExpectedTable"]]]] === 1 &&
    !FileExistsQ[FileNameJoin[{configDir, "contract.tables"}]]];

(* ---------------------------------------------------------------------- *)
(* The watchdog owns only the launched process tree.                      *)
(* ---------------------------------------------------------------------- *)

watchDir = FileNameJoin[{tempRoot, "watchdog"}];
CreateDirectory[watchDir, CreateIntermediateDirectories -> True];
treePIDFile = FileNameJoin[{watchDir, "tree-pids.txt"}];
slowTree = FileNameJoin[{watchDir, "slow-tree"}];
Export[slowTree, StringJoin[
  "#!/bin/sh\n",
  "/bin/sleep 20 &\n",
  "child=$!\n",
  "printf '%s %s\\n' \"$$\" \"$child\" > \"", treePIDFile, "\"\n",
  "wait \"$child\"\n"], "Text"];
assert["watchdog fake process tree is executable", makeExecutable[slowTree]];

decoy = StartProcess[{"/bin/sleep", "20"}];
decoyPID = FeynmanTrick`FIRE7Runner`Private`processPID[decoy];
watchResult = FeynmanTrick`FIRE7Runner`RunSupervised[
  {slowTree}, watchDir, 1, FileNameJoin[{watchDir, "supervised"}]];
treePIDs = If[FileExistsQ[treePIDFile],
  Quiet[Check[ToExpression /@ StringSplit[
    StringTrim[ReadString[treePIDFile]]], {}]], {}];
treeGone = False;
Do[
  treeGone = Length[treePIDs] === 2 &&
    AllTrue[treePIDs, processGoneOrZombieQ];
  If[treeGone, Break[]];
  Pause[.05],
  {40}
];
decoySurvived = IntegerQ[decoyPID] &&
  ProcessStatus[decoy] === "Running" && !processGoneOrZombieQ[decoyPID];
assert["watchdog reports a bounded timeout for its launched tree",
  AssociationQ[watchResult] && watchResult["ExitCode"] === 124 &&
    TrueQ[watchResult["TimedOut"]] && treeGone];
assert["watchdog leaves an unrelated decoy process alive", decoySurvived];
If[ProcessStatus[decoy] === "Running", Quiet[KillProcess[decoy]]];

(* ---------------------------------------------------------------------- *)
(* Private FIRE-interface symbol and basis-sector contracts.              *)
(* ---------------------------------------------------------------------- *)

chainVariables = FeynmanTrick`FIREInterface`Private`extractVariables[
  {Global`l^2 - Global`a},
  {Global`a -> Global`b, Global`b -> Global`x + 2 Global`m,
    Global`p^2 -> Global`s},
  {Global`l}, {Global`p}];
assert["replacement closure removes intermediates and keeps final coefficients",
  Sort[SymbolName /@ chainVariables] === {"m", "s", "x"} &&
    FreeQ[chainVariables, Alternatives[Global`a, Global`b,
      Global`l, Global`p]]];

contextLoop = Symbol["FIRE7ContractLoop`ell"];
contextExternal = Symbol["FIRE7ContractExternal`mom"];
contextMass = Symbol["FIRE7ContractKinematics`mass"];
contextInvariant = Symbol["FIRE7ContractKinematics`inv"];
contextTopology = <|
  "Propagators" -> {contextLoop^2 - contextMass},
  "Replacements" -> {contextExternal^2 -> contextInvariant},
  "LoopMomenta" -> {contextLoop},
  "ExternalMomenta" -> {contextExternal}
|>;
contextRules =
  FeynmanTrick`FIREInterface`Private`buildFIRESubstitution[contextTopology];
assert["FIRE substitution canonicalizes all non-Global input symbols",
  Sort[contextRules /. Rule[left_, right_] :>
      {Context[left] <> SymbolName[left], Context[right] <> SymbolName[right]}] ===
    Sort[{
      {"FIRE7ContractLoop`ell", "Global`ell"},
      {"FIRE7ContractExternal`mom", "Global`mom"},
      {"FIRE7ContractKinematics`mass", "Global`mass"},
      {"FIRE7ContractKinematics`inv", "Global`inv"}
    }]];

dimensionCollision = Symbol["FIRE7ContractKinematics`d"];
dimensionCollisionResult = Block[{Print = Function[Null, Null, HoldAll]},
  FeynmanTrick`FIREInterface`Private`buildFIRESubstitution[<|
    "Propagators" -> {contextLoop^2 - dimensionCollision},
    "Replacements" -> {}, "LoopMomenta" -> {contextLoop},
    "ExternalMomenta" -> {}|>]];
assert["configured dimension basename collisions are rejected",
  dimensionCollisionResult === $Failed];
assert["finite-field setup accepts exact numbers and rejects machine reals",
  FeynmanTrick`FIREInterface`Private`exactFIREInputQ[{1/3, 2}] &&
    !FeynmanTrick`FIREInterface`Private`exactFIREInputQ[{1/3, 2.0}] &&
    !FeynmanTrick`FIREInterface`Private`exactFIREInputQ[Indeterminate] &&
    !FeynmanTrick`FIREInterface`Private`exactFIREInputQ[ComplexInfinity]];
pointA = Symbol["FIRE7ContractPoint`a"];
pointB = Symbol["FIRE7ContractPoint`b"];
pointSP = Symbol["FIRE7ContractPoint`sp"];
pointRules = {pointA -> pointB, pointB -> -3/2};
assert["NumericalPoint resolves exact chains without rewriting rule left-hand sides",
  FeynmanTrick`FIREInterface`Private`validFIRENumericalPointQ[pointRules] &&
    FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
      {pointA + 1}, pointRules] === {-1/2} &&
    FeynmanTrick`FIREInterface`Private`applyFIRENumericalPointToReplacements[
      {pointSP -> pointA + 1}, pointRules] === {pointSP -> -1/2}];
assert["NumericalPoint rejects inexact and cyclic assignments",
  !FeynmanTrick`FIREInterface`Private`validFIRENumericalPointQ[
    {pointA -> 1.0}] &&
    FeynmanTrick`FIREInterface`Private`applyFIRENumericalPoint[
      pointA, {pointA -> pointB, pointB -> pointA}] === $Failed];

loadRoot = FileNameJoin[{tempRoot, "load-provenance"}];
CreateDirectory[loadRoot, CreateIntermediateDirectories -> True];
loadSource = FileNameJoin[{loadRoot, "FIRE7.m"}];
Export[loadSource, "FIRE`SaveStart[___] := Null;\n", "Text"];
oldFirePath = FeynmanTrick`FTConfiguration[]["FIREPath"];
FeynmanTrick`SetFTOption["FIREPath", loadRoot];
Unprotect[FIRE`SaveStart]; Clear[FIRE`SaveStart];
FeynmanTrick`FIREInterface`Private`$LoadedFIRE7Source = None;
FeynmanTrick`FIREInterface`Private`$LoadedFIRE7SourceSHA256 = None;
firstLoad = FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[];
Export[loadSource,
  "FIRE`SaveStart[___] := Null;\n(* changed in place *)\n", "Text"];
secondLoad = Block[{Print = Function[Null, Null, HoldAll]},
  FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[]];
changedFingerprint =
  FeynmanTrick`FIREInterface`Private`currentFIRERuntimeFingerprintRecord[];
assert["loaded FIRE7 source is bound to its content hash, not only its path",
  firstLoad === True && secondLoad === $Failed &&
    MissingQ[changedFingerprint["FIRESourceSHA256"]]];
Unprotect[FIRE`SaveStart]; Clear[FIRE`SaveStart];
FeynmanTrick`FIREInterface`Private`$LoadedFIRE7Source = None;
FeynmanTrick`FIREInterface`Private`$LoadedFIRE7SourceSHA256 = None;
FeynmanTrick`SetFTOption["FIREPath", oldFirePath];

collisionA = Symbol["FIRE7ContractA`same"];
collisionB = Symbol["FIRE7ContractB`same"];
collisionResult = Block[{Print = Function[Null, Null, HoldAll]},
  FeynmanTrick`FIREInterface`Private`extractVariables[
    {collisionA + collisionB}, {}, {}, {}]];
assert["context basename collisions are rejected before FIRE serialization",
  collisionResult === $Failed];

oldDimension = FeynmanTrick`FTConfiguration[]["DimensionVariable"];
oldBackend = FeynmanTrick`FTConfiguration[]["FIREBackend"];
FeynmanTrick`SetFTOption["DimensionVariable", Global`d];
FeynmanTrick`SetFTOption["FIREBackend", "Modular"];
orderedVariables = FeynmanTrick`FIREInterface`Private`orderedFIREVariables[{
  <|"Variables" -> {Global`m, Global`d, Global`x, Global`m}|>
}];
invalidVariable = Symbol["Global`BadName"];
invalidOrdered = Block[{Print = Function[Null, Null, HoldAll]},
  FeynmanTrick`FIREInterface`Private`orderedFIREVariables[{
    <|"Variables" -> {Global`x, invalidVariable}|>
  }]];
tooManyVariables = Block[{Print = Function[Null, Null, HoldAll]},
  FeynmanTrick`FIREInterface`Private`orderedFIREVariables[{
    <|"Variables" ->
      (Symbol["Global`v" <> ToString[#]] & /@ Range[16])|>
  }]];
FeynmanTrick`SetFTOption["FIREBackend", "Classical"];
classicalManyVariables =
  FeynmanTrick`FIREInterface`Private`orderedFIREVariables[{
    <|"Variables" ->
      (Symbol["Global`v" <> ToString[#]] & /@ Range[16])|>
  }];
FeynmanTrick`SetFTOption["DimensionVariable", oldDimension];
FeynmanTrick`SetFTOption["FIREBackend", oldBackend];
assert["FIRE variables are stable-deduplicated with d forced last",
  orderedVariables === {Global`m, Global`x, Global`d}];
assert["unsafe FIRE variable spellings are rejected",
  invalidOrdered === $Failed &&
    FeynmanTrick`FIREInterface`Private`safeFIREVariableNameQ["x1"] &&
    !FeynmanTrick`FIREInterface`Private`safeFIREVariableNameQ["X"] &&
    !FeynmanTrick`FIREInterface`Private`safeFIREVariableNameQ["x_y"] &&
    !FeynmanTrick`FIREInterface`Private`safeFIREVariableNameQ["1x"]];
assert["modular reconstruction enforces FIRE7's 16-variable array bound",
  tooManyVariables === $Failed && Length[classicalManyVariables] === 17];

coverageTable = FileNameJoin[{tempRoot, "raw-coverage.tables"}];
Put[{
  {{101, {{201, 3}}}, {102, {}}},
  {{101, {7, {2, 1}}}, {102, {7, {1, 1}}},
   {201, {7, {0, 1}}}}
}, coverageTable];
assert["raw prime-table coverage maps relation IDs without coefficients",
  Sort[FeynmanTrick`FIREInterface`Private`rawFIRETableCoverage[
    coverageTable]] === Sort[{{7, {2, 1}}, {7, {1, 1}}}]];

parserSentinel = FileNameJoin[{tempRoot, "parser-sentinel.tables"}];
Export[parserSentinel, "nonempty\n", "Text"];
malformedRHSResult = Block[{FIRE`Tables2Rules, FIRE`Tables2Masters},
  FIRE`Tables2Rules[___] := {
    Global`G[7, {2, 1}] -> Global`G[7, 1]};
  FIRE`Tables2Masters[___] := {};
  FeynmanTrick`FIREInterface`Private`parseAndValidateFIRETable[
    parserSentinel, {{7, {2, 1}}}, {7}]];
assert["exact table parser rejects malformed master calls on rule right-hand sides",
  malformedRHSResult === $Failed];
duplicateIdentityResult = Block[{FIRE`Tables2Rules, FIRE`Tables2Masters},
  FIRE`Tables2Rules[___] := {
    Global`G[7, {2, 1}] -> Global`G[7, {2, 1}],
    Global`G[7, {2, 1}] -> Global`G[7, {2, 1}]};
  FIRE`Tables2Masters[___] := {{7, {2, 1}}};
  FeynmanTrick`FIREInterface`Private`parseAndValidateFIRETable[
    parserSentinel, {{7, {2, 1}}}, {7}]];
conflictingDuplicateResult = Block[{FIRE`Tables2Rules, FIRE`Tables2Masters},
  FIRE`Tables2Rules[___] := {
    Global`G[7, {2, 1}] -> Global`G[7, {2, 1}],
    Global`G[7, {2, 1}] -> 2 Global`G[7, {2, 1}]};
  FIRE`Tables2Masters[___] := {{7, {2, 1}}};
  FeynmanTrick`FIREInterface`Private`parseAndValidateFIRETable[
    parserSentinel, {{7, {2, 1}}}, {7}]];
assert["exact duplicate FIRE7 identity rules are canonicalized while conflicts fail",
  AssociationQ[duplicateIdentityResult] &&
    duplicateIdentityResult["Rules"] === {
      Global`G[7, {2, 1}] -> Global`G[7, {2, 1}]} &&
    conflictingDuplicateResult === $Failed];
assert["basis master validation binds each problem's requested index arity",
  FeynmanTrick`FIREInterface`Private`validFIREPairsForRequestsQ[
    {{7, {1, 1}}}, {{7, {2, 1}}}] &&
    !FeynmanTrick`FIREInterface`Private`validFIREPairsForRequestsQ[
      {{7, {1, 1, 0}}}, {{7, {2, 1}}}]];

sectorTopology = <|
  "NumPropagators" -> 5,
  "EliminatedPositions" -> {2, 2, 0, 99},
  "NumeratorPositions" -> {4, 4, -1},
  "BasisProbeIntegrals" -> {
    {3, 0, 1, -2, 1},
    {3, 0, 1, -2, 1},
    {1, 1, 1, 0, 1},
    {1, 0, 1, 1, 1}}
|>;
sectors = FeynmanTrick`FIREInterface`Private`basisSectors[sectorTopology];
sectorSeeds =
  FeynmanTrick`FIREInterface`Private`basisSeedIntegrals[sectorTopology];
expectedSectors = {
  {0, 0, 0, 0, 1},
  {0, 0, 1, 0, 0},
  {0, 0, 1, 0, 1},
  {1, 0, 0, 0, 0},
  {1, 0, 0, 0, 1},
  {1, 0, 1, 0, 0},
  {1, 0, 1, 0, 1}
};
assert["basis sectors enumerate only active denominator corners",
  sectors === expectedSectors &&
    AllTrue[sectors, #[[2]] === 0 && #[[4]] === 0 &] &&
    FreeQ[sectors, ConstantArray[0, 5]]];
assert["basis discovery adds one-dot and one-ISP first-shell probes",
  Length[sectorSeeds] === 27 && SubsetQ[sectorSeeds, sectors] &&
    MemberQ[sectorSeeds, {2, 0, 1, 0, 1}] &&
    MemberQ[sectorSeeds, {1, 0, 2, 0, 1}] &&
    MemberQ[sectorSeeds, {1, 0, 1, 0, 2}] &&
    MemberQ[sectorSeeds, {1, 0, 1, -1, 1}] &&
    MemberQ[sectorSeeds, {3, 0, 1, -2, 1}] &&
    AllTrue[sectorSeeds, #[[2]] === 0 && #[[4]] <= 0 &]];
assert["basis-sector pruning returns no request when every slot is inactive",
  FeynmanTrick`FIREInterface`Private`basisSectors[<|
    "NumPropagators" -> 3,
    "EliminatedPositions" -> {1},
    "NumeratorPositions" -> {2, 3}
  |>] === {}];

If[DirectoryQ[tempRoot],
  Quiet[DeleteDirectory[tempRoot, DeleteContents -> True]]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
