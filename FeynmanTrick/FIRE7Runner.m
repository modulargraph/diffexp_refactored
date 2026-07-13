(* ::Package:: *)
(* FIRE7Runner -- isolated process, modular-cache, and artifact mechanics.

   This module deliberately knows nothing about Feynman integral semantics.
   FIREInterface remains responsible for parsing and validating tables. *)

BeginPackage["FeynmanTrick`FIRE7Runner`"];

ParseConfig::usage = "ParseConfig[dir, stem] parses the generated FIRE 7 config into an exact run manifest.";
ExpectedReconstructedFile::usage = "ExpectedReconstructedFile[dir, manifest] returns FIRE7_MPI's exact reconstructed table path.";
ExpectedPrimeFile::usage = "ExpectedPrimeFile[dir, manifest, values] returns FIRE7p's exact sampled table path.";
ModularCommand::usage = "ModularCommand[fireRoot, manifest, settings] constructs the FIRE7_MPI argv vector without invoking a shell.";
RunSupervised::usage = "RunSupervised[argv, dir, timeout, logBase] runs one process tree with a bounded targeted watchdog.";
RunModular::usage = "RunModular[fireRoot, dir, stem, settings] resumes or performs exact FIRE7_MPI reconstruction and publishes the canonical table.";
RunBasisProbes::usage = "RunBasisProbes[fireRoot, dir, stem, settings] produces independent FIRE7p finite-field tables for basis discovery.";

Begin["`Private`"];

$jobSchema = "FeynmanTrick.FIRE7ModularJob/v2";
$jobCacheVersion = "v2";
$reductionSamplingPolicy = "geometric-zippel-adaptive-limits/v2";
$lockAcquireTimeoutSeconds = 60;
$lockOwnerGraceSeconds = 30;
If[!AssociationQ[$OwnedLockNonces], $OwnedLockNonces = <||>];

lockWaitSeconds[settings_Association, operation_String:"Reduction"] := Module[
  {jobTimeout, stages, minimum, configured},
  jobTimeout = Lookup[settings, "TimeoutSeconds", 1800];
  stages = If[operation === "Basis",
    Lookup[settings, "BasisProbeCount", 2], 2];
  minimum = If[IntegerQ[jobTimeout] && IntegerQ[stages] && stages >= 1,
    stages jobTimeout + 60, $Failed];
  configured = Replace[Lookup[settings, "LockTimeoutSeconds", Automatic],
    Automatic :> If[IntegerQ[minimum],
      Max[$lockAcquireTimeoutSeconds, minimum], $Failed]];
  If[IntegerQ[jobTimeout] && jobTimeout >= 1 &&
      IntegerQ[minimum] && IntegerQ[configured] && configured >= minimum,
    configured,
    failure["FIRE lock timeout is too short for all stages of one job",
      <|"JobTimeoutSeconds" -> jobTimeout,
        "StageCount" -> stages, "RequiredLockSeconds" -> minimum,
        "LockTimeoutSeconds" -> configured|>]]
];

failure[detail_String, data_:<||>] := Failure["FIRE7Runner",
  Join[<|"Detail" -> detail|>, data]];

$modularVariableLimit = 16;

validateModularVariableCount[manifest_Association] := Module[{count},
  count = Length[Lookup[manifest, "Variables", {}]];
  If[1 <= count <= $modularVariableLimit, True,
    failure[
      "FIRE7 modular reconstruction supports at most 16 variables including the dimension",
      <|"VariableCount" -> count,
        "MaximumVariableCount" -> $modularVariableLimit|>]]
];

safeReadString[path_String] := Module[{text},
  If[!FileExistsQ[path], Return["", Module]];
  text = Quiet[Check[ReadString[path], ""]];
  If[StringQ[text], text, ""]
];

nonemptyFileQ[path_String] := FileExistsQ[path] &&
  Quiet[Check[FileByteCount[path] > 0, False]];

fileSHA256[path_String] := If[FileExistsQ[path],
  IntegerString[FileHash[path, "SHA256"], 16, 64], Missing["NotFound"]];

safeLeafNameQ[name_String] := name =!= "" &&
  FileNameTake[name] === name && !MemberQ[{".", ".."}, name] &&
  !StringContainsQ[name, {"\n", "\r", FromCharacterCode[0]}];

safeVariableNameQ[name_String] :=
  StringMatchQ[name, RegularExpression["[a-z][a-z0-9]*"]];

configPath[dir_String, stem_String] := Module[{leaf},
  leaf = If[StringEndsQ[stem, ".config"], stem, stem <> ".config"];
  If[!safeLeafNameQ[leaf], Return[$Failed, Module]];
  FileNameJoin[{dir, leaf}]
];

lineValue[lines_List, key_String] := Module[{prefix, line},
  prefix = "#" <> key;
  line = SelectFirst[lines,
    Function[item, Module[{trim = StringTrim[item], words},
      words = StringSplit[trim, Whitespace];
      words =!= {} && First[words] === prefix]], Missing["NotFound"]];
  If[MissingQ[line], Return[Missing["NotFound"], Module]];
  StringTrim[StringDrop[StringTrim[line], StringLength[prefix]]]
];

ParseConfig[dir_String, stem_String] := Module[
  {path, text, lines, variablesText, variables, output, integrals,
   problemLines, problems, words, configLeaf, inputLeaves},
  path = configPath[dir, stem];
  If[path === $Failed || !nonemptyFileQ[path],
    Return[failure["FIRE config is missing or empty",
      <|"Directory" -> dir, "Config" -> stem|>], Module]];
  text = safeReadString[path];
  lines = StringSplit[text, {"\r\n", "\n", "\r"}];
  variablesText = lineValue[lines, "variables"];
  output = lineValue[lines, "output"];
  integrals = lineValue[lines, "integrals"];
  If[AnyTrue[{variablesText, output, integrals}, MissingQ],
    Return[failure["FIRE config lacks #variables, #integrals, or #output",
      <|"ConfigPath" -> path|>], Module]];
  variables = StringTrim /@ StringSplit[variablesText, ","];
  If[variables === {} || MemberQ[variables, ""] ||
      !DuplicateFreeQ[variables] || !AllTrue[variables, safeVariableNameQ],
    Return[failure[
      "FIRE config variables must be unique lowercase alphanumeric names",
      <|"ConfigPath" -> path, "Variables" -> variables|>], Module]];
  problemLines = Select[lines, Function[item,
    With[{w = StringSplit[StringTrim[item], Whitespace]},
      Length[w] >= 1 && First[w] === "#problem"]]];
  problems = Map[Function[line,
    words = StringSplit[StringTrim[line], Whitespace];
    If[Length[words] =!= 3 ||
       !StringMatchQ[words[[2]], RegularExpression["[0-9]+"]],
      Return[failure["malformed #problem line",
        <|"ConfigPath" -> path, "Line" -> line|>], Module]];
    <|"ProblemNumber" -> ToExpression[words[[2]]],
      "StartFile" -> words[[3]]|>], problemLines];
  If[FailureQ[problems], Return[problems, Module]];
  If[problems === {},
    Return[failure["FIRE config has no #problem entries",
      <|"ConfigPath" -> path|>], Module]];
  If[!DuplicateFreeQ[Lookup[problems, "ProblemNumber"]],
    Return[failure["FIRE config has duplicate problem numbers",
      <|"ConfigPath" -> path|>], Module]];
  configLeaf = FileNameTake[path];
  inputLeaves = Join[{configLeaf, integrals},
    Lookup[problems, "StartFile"]];
  If[!AllTrue[Join[{output}, inputLeaves], safeLeafNameQ],
    Return[failure["generated FIRE inputs must use safe leaf filenames",
      <|"ConfigPath" -> path|>], Module]];
  If[!StringEndsQ[output, ".tables"] ||
      !StringEndsQ[integrals, ".m"] ||
      !AllTrue[Lookup[problems, "StartFile"],
        StringEndsQ[#, ".start"] &],
    Return[failure["FIRE config requires .tables, .m, and .start filenames",
      <|"ConfigPath" -> path, "Output" -> output,
        "Integrals" -> integrals|>], Module]];
  If[MemberQ[inputLeaves, output] || !DuplicateFreeQ[inputLeaves],
    Return[failure["FIRE output and input filenames must be distinct",
      <|"ConfigPath" -> path, "Output" -> output,
        "Inputs" -> inputLeaves|>], Module]];
  <|
    "ConfigPath" -> path,
    "ConfigFile" -> configLeaf,
    "ConfigStem" -> FileBaseName[configLeaf],
    "Variables" -> variables,
    "OutputFile" -> output,
    "IntegralsFile" -> integrals,
    "Problems" -> problems
  |>
];

tableBase[manifest_Association] := FileBaseName[manifest["OutputFile"]];

tableDirectory[dir_String, manifest_Association] := Module[{relative},
  relative = DirectoryName[manifest["OutputFile"]];
  If[MemberQ[{"", "."}, relative], dir,
    FileNameJoin[{dir, relative}]]
];

ExpectedReconstructedFile[dir_String, manifest_Association] :=
  FileNameJoin[{tableDirectory[dir, manifest], tableBase[manifest] <> "_" <>
    StringRiffle[manifest["Variables"], "_"] <> "_0.tables"}];

ExpectedPrimeFile[dir_String, manifest_Association, values_List] :=
  FileNameJoin[{tableDirectory[dir, manifest], tableBase[manifest] <> "_" <>
    StringRiffle[ToString[#, InputForm] & /@ values, "_"] <> ".tables"}];

installedMultiprimeWidth[fireRoot_String] := Module[{path, text, hit},
  path = FileNameJoin[{fireRoot, "paths.inc"}];
  text = safeReadString[path];
  hit = StringCases[text,
    RegularExpression["(?m)^RESULTING_MPRIME_COUNT=([0-9]+)\\s*$"] :> "$1"];
  If[Length[hit] === 1, ToExpression[First[hit]], Missing["Unavailable"]]
];

resolveExecutable[command_String] := Module[{explicit, candidates, found},
  explicit = !MemberQ[{"", "."}, DirectoryName[command]];
  If[explicit,
    Return[If[FileExistsQ[ExpandFileName[command]],
      ExpandFileName[command], Missing["NotFound"]], Module]];
  candidates = FileNameJoin[{#, command}] & /@
    StringSplit[Environment["PATH"],
      If[$OperatingSystem === "Windows", ";", ":"]];
  found = SelectFirst[candidates, FileExistsQ, Missing["NotFound"]];
  If[StringQ[found], ExpandFileName[found], found]
];

resolvedMPIExecutable[settings_Association] := Module[{configured},
  configured = Replace[Lookup[settings, "MPIExecutable", "mpirun"],
    Automatic -> "mpirun"];
  If[StringQ[configured], resolveExecutable[configured], Missing["NotFound"]]
];

inputPaths[dir_String, manifest_Association] := Join[
  {manifest["ConfigPath"], FileNameJoin[{dir, manifest["IntegralsFile"]}]},
  FileNameJoin[{dir, #}] & /@ Lookup[manifest["Problems"], "StartFile"]
];

jobIdentity[operation_String, fireRoot_String, dir_String,
    manifest_Association, settings_Association] := Module[
  {paths, useMultiprime, mpiLauncher, record, expectedInputs,
   expectedRuntime},
  paths = inputPaths[dir, manifest];
  If[!AllTrue[paths, nonemptyFileQ],
    Return[failure["one or more FIRE job inputs are missing or empty",
      <|"Inputs" -> paths|>], Module]];
  useMultiprime = TrueQ[Lookup[settings, "UseMultiprime", True]];
  mpiLauncher = If[operation === "Reduction",
    resolvedMPIExecutable[settings], Missing["Inactive"]];
  If[operation === "Reduction" && !StringQ[mpiLauncher],
    Return[failure["MPI launcher is missing",
      <|"MPIExecutable" -> Lookup[settings, "MPIExecutable", "mpirun"]|>],
      Module]];
  record = <|
    "Schema" -> $jobSchema,
    "Operation" -> operation,
    "InputFiles" -> AssociationThread[FileNameTake /@ paths,
      fileSHA256 /@ paths],
    "Variables" -> manifest["Variables"],
    "Problems" -> manifest["Problems"],
    "OutputFile" -> manifest["OutputFile"],
    "Calc" -> Lookup[settings, "Calc", "flint"],
    "UseMultiprime" -> useMultiprime,
    "PrimeLimit" -> Lookup[settings, "PrimeLimit", 127],
    "DimensionSeparated" ->
      TrueQ[Lookup[settings, "DimensionSeparated", False]],
    "SamplingPolicy" -> If[operation === "Reduction",
      $reductionSamplingPolicy, Missing["Inactive"]],
    "VariableBasePolicy" -> If[operation === "Reduction",
      "big-primes/v1", Missing["Inactive"]],
    "FIRE7SourceSHA256" -> fileSHA256[
      FileNameJoin[{fireRoot, "FIRE7.m"}]],
    "FIRE7PrimeSHA256" -> If[MemberQ[{"Basis", "Reduction"}, operation],
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7p"}]],
      Missing["Inactive"]],
    "FIRE7MultiprimeSHA256" -> If[operation === "Reduction" && useMultiprime,
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7mp"}]],
      Missing["Inactive"]],
    "FIRE7MPISHA256" -> If[operation === "Reduction",
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7_MPI"}]],
      Missing["Inactive"]],
    "ReconstructSHA256" -> If[operation === "Reduction",
      fileSHA256[FileNameJoin[{fireRoot, "bin", "reconstruct"}]],
      Missing["Inactive"]],
    "MPILauncherPath" -> mpiLauncher,
    "MPILauncherSHA256" -> If[operation === "Reduction",
      fileSHA256[mpiLauncher], Missing["Inactive"]],
    "PathsIncSHA256" -> If[useMultiprime,
      fileSHA256[FileNameJoin[{fireRoot, "paths.inc"}]],
      Missing["Inactive"]],
    "InstalledMultiprimeWidth" -> If[useMultiprime,
      installedMultiprimeWidth[fireRoot], Missing["Inactive"]],
    "BasisProbeCount" -> If[operation === "Basis",
      Lookup[settings, "BasisProbeCount", 2], Missing["Inactive"]],
    "BasisProbePolicy" -> If[operation === "Basis",
      "independent-sha256/v1", Missing["Inactive"]]
  |>;
  expectedInputs = Lookup[settings, "ExpectedInputHashes", <||>];
  If[!AssociationQ[expectedInputs] ||
      !And @@ KeyValueMap[
        Lookup[record["InputFiles"], #1, Missing["Absent"]] === #2 &,
        expectedInputs],
    Return[failure["FIRE inputs do not match the prepared topology snapshot",
      <|"ExpectedInputHashes" -> expectedInputs,
        "ActualInputHashes" -> record["InputFiles"]|>], Module]];
  expectedRuntime = Lookup[settings, "ExpectedRuntimeHashes", <||>];
  If[!AssociationQ[expectedRuntime] ||
      !And @@ KeyValueMap[
        Lookup[record, #1, Missing["Absent"]] === #2 &,
        expectedRuntime],
    Return[failure["FIRE runtime does not match the prepared topology snapshot",
      <|"ExpectedRuntimeHashes" -> expectedRuntime|>], Module]];
  <|"Record" -> record,
    "ID" -> IntegerString[Hash[record, "SHA256"], 16, 64]|>
];

runtimeIdentityMatchesQ[job_Association, fireRoot_String] := Module[
  {identity, operation, useMultiprime, mpiLauncher, current},
  identity = job["Identity"];
  operation = identity["Operation"];
  useMultiprime = TrueQ[identity["UseMultiprime"]];
  mpiLauncher = Lookup[identity, "MPILauncherPath", Missing["Inactive"]];
  current = <|
    "FIRE7SourceSHA256" -> fileSHA256[FileNameJoin[{fireRoot, "FIRE7.m"}]],
    "FIRE7PrimeSHA256" -> If[
      MemberQ[{"Basis", "Reduction"}, operation],
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7p"}]],
      Missing["Inactive"]],
    "FIRE7MultiprimeSHA256" -> If[
      operation === "Reduction" && useMultiprime,
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7mp"}]],
      Missing["Inactive"]],
    "FIRE7MPISHA256" -> If[operation === "Reduction",
      fileSHA256[FileNameJoin[{fireRoot, "bin", "FIRE7_MPI"}]],
      Missing["Inactive"]],
    "ReconstructSHA256" -> If[operation === "Reduction",
      fileSHA256[FileNameJoin[{fireRoot, "bin", "reconstruct"}]],
      Missing["Inactive"]],
    "MPILauncherPath" -> mpiLauncher,
    "MPILauncherSHA256" -> If[operation === "Reduction" &&
        StringQ[mpiLauncher], fileSHA256[mpiLauncher], Missing["Inactive"]],
    "PathsIncSHA256" -> If[useMultiprime,
      fileSHA256[FileNameJoin[{fireRoot, "paths.inc"}]],
      Missing["Inactive"]],
    "InstalledMultiprimeWidth" -> If[useMultiprime,
      installedMultiprimeWidth[fireRoot], Missing["Inactive"]]
  |>;
  And @@ KeyValueMap[Lookup[identity, #1, Missing["Absent"]] === #2 &,
    current]
];

atomicPublish[source_String, destination_String] := Module[{temp, result},
  If[!nonemptyFileQ[source], Return[$Failed, Module]];
  If[!DirectoryQ[DirectoryName[destination]],
    CreateDirectory[DirectoryName[destination],
      CreateIntermediateDirectories -> True]];
  temp = destination <> ".tmp-" <> StringReplace[CreateUUID[], "-" -> ""];
  result = Quiet[Check[
    CopyFile[source, temp, OverwriteTarget -> True];
    RenameFile[temp, destination, OverwriteTarget -> True];
    destination,
    $Failed]];
  If[FileExistsQ[temp], Quiet[DeleteFile[temp]]];
  result
];

atomicExportWXF[path_String, expression_] := Module[{temp, result},
  If[!DirectoryQ[DirectoryName[path]],
    CreateDirectory[DirectoryName[path],
      CreateIntermediateDirectories -> True]];
  temp = path <> ".tmp-" <> StringReplace[CreateUUID[], "-" -> ""];
  result = Quiet[Check[
    If[Export[temp, expression, "WXF"] === $Failed,
      Return[$Failed, Module]];
    RenameFile[temp, path, OverwriteTarget -> True];
    path,
    $Failed]];
  If[FileExistsQ[temp], Quiet[DeleteFile[temp]]];
  result
];

safeImportWXF[path_String] := If[FileExistsQ[path],
  Quiet[Check[Import[path, "WXF"], $Failed]], $Failed];

jobManifest[manifest_Association] := ReplacePart[manifest,
  "OutputFile" -> FileNameJoin[{"artifacts", manifest["OutputFile"]}]];

copyJobInputs[dir_String, jobDir_String, manifest_Association,
    expectedHashes_Association] := Module[
  {paths = inputPaths[dir, manifest], target, expected, sourceBefore,
   sourceAfter, configTarget, text, lines},
  If[!DirectoryQ[jobDir],
    CreateDirectory[jobDir, CreateIntermediateDirectories -> True]];
  If[!DirectoryQ[FileNameJoin[{jobDir, "artifacts"}]],
    CreateDirectory[FileNameJoin[{jobDir, "artifacts"}],
      CreateIntermediateDirectories -> True]];
  Do[
    target = FileNameJoin[{jobDir, FileNameTake[path]}];
    expected = Lookup[expectedHashes, FileNameTake[path], Missing["Absent"]];
    sourceBefore = fileSHA256[path];
    If[MissingQ[expected] || sourceBefore =!= expected,
      Return[$Failed, Module]];
    If[!FileExistsQ[target] || fileSHA256[target] =!= expected,
      Quiet[Check[CopyFile[path, target, OverwriteTarget -> True],
        Return[$Failed, Module]]]];
    sourceAfter = fileSHA256[path];
    If[sourceAfter =!= expected || fileSHA256[target] =!= expected,
      Return[$Failed, Module]],
    {path, paths}];
  (* FIRE 7 calls create_directories on the output parent without guarding
     the empty-path case.  Give every cached job an explicit artifact parent
     while retaining the original config in the content identity. *)
  configTarget = FileNameJoin[{jobDir, manifest["ConfigFile"]}];
  text = safeReadString[configTarget];
  lines = StringSplit[text, {"\r\n", "\n", "\r"}];
  lines = Replace[lines, line_ /;
      StringStartsQ[StringTrim[line], "#output"] :>
        "#output            " <> jobManifest[manifest]["OutputFile"], {1}];
  If[Export[configTarget, StringRiffle[lines, "\n"] <> "\n", "Text"] === $Failed,
    Return[$Failed, Module]];
  True
];

processPID[process_ProcessObject] :=
  Replace[process, ProcessObject[data_] :> Lookup[data, "PID", Missing["PID"]]];

childPIDs[pid_Integer] := Module[{result, text},
  result = Quiet[Check[RunProcess[{"/usr/bin/pgrep", "-P", ToString[pid]}],
    <||>]];
  text = Lookup[result, "StandardOutput", ""];
  If[!StringQ[text] || StringTrim[text] === "", {},
    Select[Quiet[Check[ToExpression /@ StringSplit[StringTrim[text]], {}]],
      IntegerQ]]
];

descendantPIDs[pid_Integer] := Module[{children},
  children = childPIDs[pid];
  DeleteDuplicates[Join[children, Flatten[descendantPIDs /@ children]]]
];

signalPIDs[pids_List, signal_String] := If[pids =!= {},
  Quiet[Check[RunProcess[Join[{"/bin/kill", "-" <> signal},
    ToString /@ Select[pids, IntegerQ]]], Null]]];

pidGoneOrZombieQ[pid_Integer] := Module[{probe, state},
  probe = Quiet[Check[RunProcess[
    {"/bin/ps", "-o", "state=", "-p", ToString[pid]}], <||>]];
  state = StringTrim[Lookup[probe, "StandardOutput", ""]];
  Lookup[probe, "ExitCode", 1] =!= 0 || state === "" ||
    StringStartsQ[state, "Z"]
];

waitForPIDsToStop[pids_List, seconds_:2] := Module[{remaining, deadline},
  remaining = DeleteDuplicates[Select[pids, IntegerQ]];
  deadline = AbsoluteTime[] + seconds;
  While[remaining =!= {} && AbsoluteTime[] < deadline,
    remaining = Select[remaining, !pidGoneOrZombieQ[#] &];
    If[remaining =!= {}, Pause[.05]]];
  remaining
];

drainStream[input_, output_] := Module[{chunk},
  chunk = Quiet[Check[ReadString[input, EndOfBuffer], ""]];
  If[StringQ[chunk] && chunk =!= "", WriteString[output, chunk]];
];

RunSupervised[argv_List, dir_String, timeout_Integer, logBase_String] := Module[
  {process, pid, descendants = {}, outPath, errPath, outStream, errStream,
   start, timedOut = False, status, exitCode, elapsed, poll = .05,
   timeoutStream, survivors = {}},
  If[!AllTrue[argv, StringQ] || argv === {},
    Return[failure["process argv must be a nonempty string vector"], Module]];
  If[!DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  outPath = logBase <> ".stdout.log";
  errPath = logBase <> ".stderr.log";
  outStream = OpenWrite[outPath, PageWidth -> Infinity];
  errStream = OpenWrite[errPath, PageWidth -> Infinity];
  process = Quiet[Check[StartProcess[argv, ProcessDirectory -> dir], $Failed]];
  If[process === $Failed,
    Close[outStream]; Close[errStream];
    Return[failure["failed to start process", <|"Argv" -> argv|>], Module]];
  pid = processPID[process];
  start = AbsoluteTime[];
  CheckAbort[
    While[ProcessStatus[process] === "Running",
      drainStream[process["StandardOutput"], outStream];
      drainStream[process["StandardError"], errStream];
      elapsed = AbsoluteTime[] - start;
      If[elapsed >= timeout,
        timedOut = True;
        If[IntegerQ[pid], descendants = descendantPIDs[pid]];
        (* Signal the launcher/coordinator first so it can tear down ranks. *)
        Quiet[KillProcess[process]];
        Pause[.25];
        signalPIDs[descendants, "TERM"];
        Pause[1.0];
        signalPIDs[Join[If[IntegerQ[pid], {pid}, {}], descendants], "KILL"];
        survivors = waitForPIDsToStop[
          Join[If[IntegerQ[pid], {pid}, {}], descendants], 2];
        If[survivors =!= {},
          signalPIDs[survivors, "KILL"];
          survivors = waitForPIDsToStop[survivors, 2]];
        Break[]];
      Pause[poll];
      If[elapsed > 2, poll = .25];
    ],
    If[IntegerQ[pid], descendants = descendantPIDs[pid]];
    Quiet[KillProcess[process]];
    Pause[.1];
    signalPIDs[descendants, "TERM"];
    Pause[.5];
    signalPIDs[Join[If[IntegerQ[pid], {pid}, {}], descendants], "KILL"];
    survivors = waitForPIDsToStop[
      Join[If[IntegerQ[pid], {pid}, {}], descendants], 2];
    If[survivors =!= {},
      signalPIDs[survivors, "KILL"];
      survivors = waitForPIDsToStop[survivors, 2]];
    Close[outStream]; Close[errStream]; Abort[]
  ];
  drainStream[process["StandardOutput"], outStream];
  drainStream[process["StandardError"], errStream];
  Close[outStream]; Close[errStream];
  If[timedOut,
    timeoutStream = OpenAppend[errPath];
    WriteString[timeoutStream,
      "FIRE7 timeout after " <> ToString[timeout] <> "s\n"];
    Close[timeoutStream];
    exitCode = 124,
    status = ProcessStatus[process];
    exitCode = Quiet[Check[ProcessInformation[process, "ExitCode"], $Failed]];
    If[!IntegerQ[exitCode] && status === "Finished",
      exitCode = Lookup[ProcessInformation[process], "ExitCode", $Failed]]
  ];
  <|"ExitCode" -> exitCode, "TimedOut" -> timedOut,
    "PID" -> pid, "Argv" -> argv, "Directory" -> dir,
    "StandardOutputFile" -> outPath, "StandardErrorFile" -> errPath,
    "ElapsedSeconds" -> N[AbsoluteTime[] - start],
    "SurvivingPIDs" -> survivors|>
];

processAliveQ[pid_Integer] := Lookup[
  Quiet[Check[RunProcess[{"/bin/kill", "-0", ToString[pid]}], <||>]],
  "ExitCode", 1] === 0;

lockAgeSeconds[lock_String] := Quiet[Check[
  Max[0, AbsoluteTime[] - AbsoluteTime[FileDate[lock]]], 0]];

validLockOwnerQ[owner_] := AssociationQ[owner] &&
  IntegerQ[Lookup[owner, "PID", Missing["PID"]]] &&
  StringQ[Lookup[owner, "Host", Missing["Host"]]] &&
  StringQ[Lookup[owner, "Nonce", Missing["Nonce"]]];

tryCreateLock[jobDir_String] := Module[
  {lock, candidate, owner, result = $Failed},
  lock = FileNameJoin[{jobDir, "lock"}];
  candidate = FileNameJoin[{jobDir,
    ".lock-candidate-" <> StringReplace[CreateUUID[], "-" -> ""]}];
  If[!Quiet[Check[CreateDirectory[candidate]; True, False]],
    Return[$Failed, Module]];
  owner = <|"PID" -> $ProcessID, "Host" -> $MachineName,
    "Created" -> DateString["ISODateTime"], "Nonce" -> CreateUUID[]|>;
  If[atomicExportWXF[FileNameJoin[{candidate, "owner.wxf"}], owner] =!= $Failed,
    result = Quiet[Check[RenameDirectory[candidate, lock]; lock, $Failed]]];
  If[StringQ[result], $OwnedLockNonces[lock] = owner["Nonce"]];
  If[DirectoryQ[candidate],
    Quiet[DeleteDirectory[candidate, DeleteContents -> True]]];
  result
];

acquireLock[jobDir_String] := Module[
  {lock, acquired, owner, ownerPID, ownerHost, age, stale, stalePath,
   deadline},
  If[!DirectoryQ[jobDir],
    CreateDirectory[jobDir, CreateIntermediateDirectories -> True]];
  lock = FileNameJoin[{jobDir, "lock"}];
  deadline = AbsoluteTime[] + $lockAcquireTimeoutSeconds;
  While[AbsoluteTime[] < deadline,
    acquired = tryCreateLock[jobDir];
    If[StringQ[acquired], Return[acquired, Module]];
    If[!DirectoryQ[lock], Pause[.05]; Continue[]];
    owner = safeImportWXF[FileNameJoin[{lock, "owner.wxf"}]];
    age = lockAgeSeconds[lock];
    stale = False;
    If[validLockOwnerQ[owner],
      ownerPID = owner["PID"];
      ownerHost = owner["Host"];
      stale = Which[
        ownerHost === $MachineName, !processAliveQ[ownerPID],
        (* Without a distributed heartbeat there is no sound way to prove a
           foreign owner dead.  Never steal its valid lease automatically. *)
        True, False],
      (* A just-created directory may precede visibility of its owner record
         on a network filesystem.  Missing/corrupt metadata is indeterminate
         during a grace interval, never immediate evidence of staleness. *)
      stale = age >= $lockOwnerGraceSeconds];
    If[stale,
      stalePath = lock <> ".stale-" <>
        StringReplace[CreateUUID[], "-" -> ""];
      If[Quiet[Check[RenameDirectory[lock, stalePath]; True, False]],
        Continue[]]];
    Pause[.05]
  ];
  failure["timed out acquiring modular FIRE job lock",
    <|"JobDirectory" -> jobDir,
      "TimeoutSeconds" -> $lockAcquireTimeoutSeconds|>]
];

releaseLock[lock_String] := Module[{expected, owner},
  expected = Lookup[$OwnedLockNonces, lock, Missing["NotOwned"]];
  owner = safeImportWXF[FileNameJoin[{lock, "owner.wxf"}]];
  If[StringQ[expected] && AssociationQ[owner] &&
      Lookup[owner, "Nonce", Missing["Nonce"]] === expected &&
      DirectoryQ[lock],
    Quiet[DeleteDirectory[lock, DeleteContents -> True]]];
  KeyDropFrom[$OwnedLockNonces, lock]
];

SetAttributes[withJobLock, HoldRest];
withJobLock[jobDir_String, body_] := Module[{lock, result},
  lock = acquireLock[jobDir];
  If[FailureQ[lock], Return[lock, Module]];
  Internal`WithLocalSettings[
    Null,
    result = Quiet[Check[body, $Failed]],
    releaseLock[lock]
  ];
  result
];

modularJob[operation_String, fireRoot_String, dir_String,
    manifest_Association, settings_Association] := Module[
  {identity, cacheRoot, jobDir},
  identity = jobIdentity[operation, fireRoot, dir, manifest, settings];
  If[FailureQ[identity], Return[identity, Module]];
  cacheRoot = ExpandFileName[Lookup[settings, "CacheDirectory",
    FileNameJoin[{$TemporaryDirectory, "FeynmanTrick_FIRE7_Modular"}]]];
  jobDir = FileNameJoin[{cacheRoot, $jobCacheVersion,
    StringTake[identity["ID"], 2],
    identity["ID"]}];
  If[!DirectoryQ[jobDir],
    CreateDirectory[jobDir, CreateIntermediateDirectories -> True]];
  <|"ID" -> identity["ID"], "Identity" -> identity["Record"],
    "Directory" -> jobDir|>
];

ModularCommand[fireRoot_String, manifest_Association,
    settings_Association] := Module[
  {mpi, fireMPI, workers, calc, useMultiprime, limit, command},
  mpi = resolvedMPIExecutable[settings];
  fireMPI = ExpandFileName[FileNameJoin[{fireRoot, "bin", "FIRE7_MPI"}]];
  workers = Lookup[settings, "Workers", 1];
  calc = Lookup[settings, "Calc", "flint"];
  useMultiprime = TrueQ[Lookup[settings, "UseMultiprime", True]];
  limit = Lookup[settings, "PrimeLimit", 127];
  If[MissingQ[mpi] || !FileExistsQ[fireMPI] ||
      !IntegerQ[workers] || workers < 1 || workers > 10 ||
      !StringQ[calc] || !IntegerQ[limit] || limit < 1 || limit > 127,
    Return[failure["invalid FIRE7_MPI executable or modular settings",
      <|"MPIExecutable" -> mpi, "FIRE7MPI" -> fireMPI,
        "Workers" -> workers, "Calc" -> calc, "PrimeLimit" -> limit|>],
      Module]];
  command = {mpi, "--oversubscribe", "-np", ToString[workers + 1],
    fireMPI, "--calc", calc, "--zippel", "--reconstruct",
    "--geometric", "--early_abortion",
    "--rational_reconstruction_limit", ToString[limit]};
  If[useMultiprime, command = Append[command, "--multitables"]];
  command = Append[command, "--big_primes"];
  If[TrueQ[Lookup[settings, "DimensionSeparated", False]],
    command = Append[command, "--last_separated"]];
  Join[command, {"--config", manifest["ConfigStem"]}]
];

validateMultiprimeBuild[fireRoot_String, settings_Association] := Module[
  {enabled, expected, installed},
  enabled = TrueQ[Lookup[settings, "UseMultiprime", True]];
  If[!enabled, Return[True, Module]];
  expected = Lookup[settings, "MultiprimeWidth", 16];
  installed = installedMultiprimeWidth[fireRoot];
  If[!IntegerQ[expected] || expected < 1 || !IntegerQ[installed] ||
      expected =!= installed,
    failure["FIRE7 multiprime build width does not match configuration",
      <|"ConfiguredWidth" -> expected, "InstalledWidth" -> installed|>],
    True]
];

completeMarkerPath[jobDir_String] := FileNameJoin[{jobDir, "complete.wxf"}];

completedJobQ[job_Association, final_String] := Module[{record},
  If[!nonemptyFileQ[final], Return[False, Module]];
  record = safeImportWXF[completeMarkerPath[job["Directory"]]];
  AssociationQ[record] &&
    Lookup[record, "Schema", None] === $jobSchema <> "/complete" &&
    Lookup[record, "JobID", None] === job["ID"] &&
    Lookup[record, "FinalTable", None] === FileNameTake[final] &&
    Lookup[record, "FinalSHA256", None] === fileSHA256[final]
];

quarantineFile[path_String, tag_String] := Module[{target},
  If[!FileExistsQ[path], Return[True, Module]];
  target = path <> "." <> tag <> "-" <>
    StringReplace[CreateUUID[], "-" -> ""];
  Quiet[Check[RenameFile[path, target]; True, False]]
];

quarantineIncompleteResult[job_Association, final_String] := And[
  quarantineFile[completeMarkerPath[job["Directory"]], "invalid"],
  quarantineFile[final, "unverified"]
];

writeCompleteMarker[job_Association, final_String] := Module[{record},
  record = <|"Schema" -> $jobSchema <> "/complete",
    "JobID" -> job["ID"], "FinalTable" -> FileNameTake[final],
    "FinalSHA256" -> fileSHA256[final],
    "Completed" -> DateString["ISODateTime"]|>;
  atomicExportWXF[completeMarkerPath[job["Directory"]], record]
];

settingValidatorQ[settings_Association, key_String, value_] := Module[
  {validator = Lookup[settings, key, Automatic]},
  If[validator === Automatic, Return[True, Module]];
  TrueQ[Quiet[Check[validator[value], False]]]
];

rationalReconstructionCommand[fireRoot_String, final_String,
    settings_Association] := Module[{executable, calc, limit},
  executable = ExpandFileName[
    FileNameJoin[{fireRoot, "bin", "reconstruct"}]];
  calc = Lookup[settings, "Calc", "flint"];
  limit = Lookup[settings, "PrimeLimit", 127];
  If[!FileExistsQ[executable] || !StringQ[calc] ||
      !IntegerQ[limit] || limit < 1 || limit > 127,
    Return[failure["invalid FIRE7 rational-reconstruction settings",
      <|"Executable" -> executable, "Calc" -> calc,
        "PrimeLimit" -> limit|>], Module]];
  {executable, "--method", "rational", "--calc", calc,
    ExpandFileName[final], ToString[limit]}
];

cleanupModularIntermediates[jobDir_String, final_String] := Module[{files},
  files = Join[
    FileNames["*.tables", FileNameJoin[{jobDir, "artifacts"}]],
    FileNames["*.limits", FileNameJoin[{jobDir, "artifacts"}]]];
  Scan[Function[path, If[ExpandFileName[path] =!= ExpandFileName[final],
    Quiet[DeleteFile[path]]]], files]
];

RunModular[fireRoot_String, dir_String, stem_String,
    settings_Association] := Module[
  {manifest, variableCountOK, runManifest, widthOK, job, jobDir, final,
   canonical, command, timeout, lockTimeout,
   attemptDir, run, reconstructionCommand, reconstructionRun = Missing["Unused"],
   published, currentIdentity},
  manifest = ParseConfig[dir, stem];
  If[FailureQ[manifest], Return[manifest, Module]];
  variableCountOK = validateModularVariableCount[manifest];
  If[FailureQ[variableCountOK], Return[variableCountOK, Module]];
  runManifest = jobManifest[manifest];
  widthOK = validateMultiprimeBuild[fireRoot, settings];
  If[FailureQ[widthOK], Return[widthOK, Module]];
  job = modularJob["Reduction", fireRoot, dir, manifest, settings];
  If[FailureQ[job], Return[job, Module]];
  jobDir = job["Directory"];
  final = ExpectedReconstructedFile[jobDir, runManifest];
  canonical = FileNameJoin[{dir, manifest["OutputFile"]}];
  timeout = Lookup[settings, "TimeoutSeconds", 1800];
  lockTimeout = lockWaitSeconds[settings, "Reduction"];
  If[FailureQ[lockTimeout], Return[lockTimeout, Module]];
  Block[{$lockAcquireTimeoutSeconds = lockTimeout},
  withJobLock[jobDir,
    currentIdentity = jobIdentity["Reduction", fireRoot, dir, manifest,
      settings];
    If[FailureQ[currentIdentity] || currentIdentity["ID"] =!= job["ID"] ||
        !runtimeIdentityMatchesQ[job, fireRoot],
      Return[failure["FIRE job identity changed while waiting for its lock",
        <|"JobDirectory" -> jobDir|>], Module]];
    If[copyJobInputs[dir, jobDir, manifest,
        job["Identity"]["InputFiles"]] =!= True,
      Return[failure["FIRE job inputs changed after identity construction",
        <|"JobDirectory" -> jobDir|>], Module]];
    If[completedJobQ[job, final],
      If[!settingValidatorQ[settings, "FinalValidator", final],
        quarantineIncompleteResult[job, final],
        published = atomicPublish[final, canonical];
        Return[If[published === $Failed,
          failure["failed to publish reconstructed FIRE table"],
          <|"Success" -> True, "Reused" -> True, "Job" -> job,
            "FinalTable" -> final, "CanonicalTable" -> canonical|>], Module]]];
    (* An exact-name artifact without a matching atomic completion record is
       not trusted.  Preserve it for diagnosis, but remove it from the path
       that the fresh MPI/reconstruction attempt is required to create. *)
    If[FileExistsQ[final] || FileExistsQ[completeMarkerPath[jobDir]],
      If[!quarantineIncompleteResult[job, final],
        Return[failure["could not quarantine an unverified FIRE artifact",
          <|"ExpectedTable" -> final|>], Module]]];
    command = ModularCommand[fireRoot, manifest, settings];
    If[FailureQ[command], Return[command, Module]];
    attemptDir = FileNameJoin[{jobDir, "attempts",
      StringReplace[CreateUUID[], "-" -> ""]}];
    CreateDirectory[attemptDir, CreateIntermediateDirectories -> True];
    Export[FileNameJoin[{attemptDir, "argv.wxf"}], command, "WXF"];
    run = RunSupervised[command, jobDir, timeout,
      FileNameJoin[{attemptDir, "fire7_mpi"}]];
    Export[FileNameJoin[{attemptDir, "outcome.wxf"}], run, "WXF"];
    If[FailureQ[run], Return[run, Module]];
    If[Lookup[run, "ExitCode", $Failed] === 0 && !nonemptyFileQ[final],
      (* FIRE7_MPI 7.1 can finish finite-field/Thiele sampling yet omit the
         final one-variable rational table.  Complete that documented stage
         with FIRE7's own standalone reconstructor; this is not a Classical
         fallback and it reuses the preserved modular samples. *)
      If[FileExistsQ[final], quarantineFile[final, "empty"]];
      reconstructionCommand = rationalReconstructionCommand[
        fireRoot, final, settings];
      If[!FailureQ[reconstructionCommand],
        Export[FileNameJoin[{attemptDir, "reconstruct-argv.wxf"}],
          reconstructionCommand, "WXF"];
        reconstructionRun = RunSupervised[reconstructionCommand, jobDir,
          timeout, FileNameJoin[{attemptDir, "fire7_rational"}]];
        Export[FileNameJoin[{attemptDir, "reconstruct-outcome.wxf"}],
          reconstructionRun, "WXF"],
        reconstructionRun = reconstructionCommand]];
    (* Exit zero alone is insufficient for either stage. *)
    If[Lookup[run, "ExitCode", $Failed] =!= 0 ||
        FailureQ[reconstructionRun] ||
        (AssociationQ[reconstructionRun] &&
          Lookup[reconstructionRun, "ExitCode", $Failed] =!= 0) ||
        !nonemptyFileQ[final],
      Return[failure["FIRE7 modular reconstruction did not produce the exact final table",
        <|"ExitCode" -> Lookup[run, "ExitCode", $Failed],
          "ExpectedTable" -> final, "JobDirectory" -> jobDir,
          "Attempt" -> attemptDir,
          "ReconstructionOutcome" -> reconstructionRun|>], Module]];
    If[!runtimeIdentityMatchesQ[job, fireRoot],
      quarantineFile[final, "runtime-changed"];
      Return[failure["FIRE runtime changed during modular execution",
        <|"JobDirectory" -> jobDir|>], Module]];
    If[!settingValidatorQ[settings, "FinalValidator", final],
      quarantineFile[final, "semantically-invalid"];
      Return[failure["FIRE7 reconstructed table failed semantic validation",
        <|"ExpectedTable" -> final, "Attempt" -> attemptDir|>], Module]];
    If[writeCompleteMarker[job, final] === $Failed,
      Return[failure["failed to atomically mark FIRE reconstruction complete",
        <|"ExpectedTable" -> final|>], Module]];
    published = atomicPublish[final, canonical];
    If[published === $Failed,
      Return[failure["failed to publish reconstructed FIRE table",
        <|"Source" -> final, "Destination" -> canonical|>], Module]];
    If[!TrueQ[Lookup[settings, "KeepTables", True]],
      cleanupModularIntermediates[jobDir, final]];
    <|"Success" -> True, "Reused" -> False, "Job" -> job,
      "FinalTable" -> final, "CanonicalTable" -> canonical,
      "Process" -> run, "ReconstructionProcess" -> reconstructionRun|>
  ]]
];

probeValues[count_Integer, probe_Integer, jobID_String] := Join[
  Table[17 + Mod[Hash[{
      "FeynmanTrick.FIRE7BasisProbe/v2", jobID, probe, index}, "SHA256"],
    2147483000], {index, count}],
  {probe}
];

probeMarkerPath[table_String] := table <> ".complete.wxf";

completedProbeQ[job_Association, table_String, probe_Integer,
    values_List] := Module[{record},
  If[!nonemptyFileQ[table], Return[False, Module]];
  record = safeImportWXF[probeMarkerPath[table]];
  AssociationQ[record] &&
    Lookup[record, "Schema", None] === $jobSchema <> "/probe" &&
    Lookup[record, "JobID", None] === job["ID"] &&
    Lookup[record, "Probe", None] === probe &&
    Lookup[record, "Values", None] === values &&
    Lookup[record, "Table", None] === FileNameTake[table] &&
    Lookup[record, "TableSHA256", None] === fileSHA256[table]
];

writeProbeMarker[job_Association, table_String, probe_Integer,
    values_List] := atomicExportWXF[probeMarkerPath[table], <|
  "Schema" -> $jobSchema <> "/probe", "JobID" -> job["ID"],
  "Probe" -> probe, "Values" -> values,
  "Table" -> FileNameTake[table], "TableSHA256" -> fileSHA256[table],
  "Completed" -> DateString["ISODateTime"]|>];

quarantineProbe[table_String] := And[
  quarantineFile[probeMarkerPath[table], "invalid"],
  quarantineFile[table, "unverified"]
];

RunBasisProbes[fireRoot_String, dir_String, stem_String,
    settings_Association] := Module[
  {manifest, variableCountOK, runManifest, job, jobDir, firePrime, count,
   calc, timeout, lockTimeout, tables = {},
   probeRecords = {}, values, expected, command, attemptDir, run, probes,
   groupValid, currentIdentity},
  manifest = ParseConfig[dir, stem];
  If[FailureQ[manifest], Return[manifest, Module]];
  variableCountOK = validateModularVariableCount[manifest];
  If[FailureQ[variableCountOK], Return[variableCountOK, Module]];
  runManifest = jobManifest[manifest];
  job = modularJob["Basis", fireRoot, dir, manifest, settings];
  If[FailureQ[job], Return[job, Module]];
  jobDir = job["Directory"];
  firePrime = ExpandFileName[FileNameJoin[{fireRoot, "bin", "FIRE7p"}]];
  count = Length[manifest["Variables"]];
  calc = Lookup[settings, "Calc", "flint"];
  timeout = Lookup[settings, "TimeoutSeconds", 1800];
  lockTimeout = lockWaitSeconds[settings, "Basis"];
  probes = Lookup[settings, "BasisProbeCount", 2];
  If[!FileExistsQ[firePrime] || !IntegerQ[timeout] || timeout < 1 ||
      !IntegerQ[probes] || probes < 2 || probes > 127,
    Return[failure["invalid FIRE7p executable or basis-probe settings",
      <|"FIRE7p" -> firePrime, "TimeoutSeconds" -> timeout,
        "ProbeCount" -> probes|>], Module]];
  If[FailureQ[lockTimeout], Return[lockTimeout, Module]];
  Block[{$lockAcquireTimeoutSeconds = lockTimeout},
  withJobLock[jobDir,
    currentIdentity = jobIdentity["Basis", fireRoot, dir, manifest, settings];
    If[FailureQ[currentIdentity] || currentIdentity["ID"] =!= job["ID"] ||
        !runtimeIdentityMatchesQ[job, fireRoot],
      Return[failure["FIRE basis identity changed while waiting for its lock",
        <|"JobDirectory" -> jobDir|>], Module]];
    If[copyJobInputs[dir, jobDir, manifest,
        job["Identity"]["InputFiles"]] =!= True,
      Return[failure["FIRE basis inputs changed after identity construction"],
        Module]];
    Do[
      values = probeValues[count, probe, job["ID"]];
      expected = ExpectedPrimeFile[jobDir, runManifest, values];
      If[!completedProbeQ[job, expected, probe, values],
        If[FileExistsQ[expected] || FileExistsQ[probeMarkerPath[expected]],
          If[!quarantineProbe[expected],
            Return[failure["could not quarantine an unverified basis probe",
              <|"Probe" -> probe, "ExpectedTable" -> expected|>], Module]]];
        command = {firePrime, "--calc", calc,
          "--variables", StringRiffle[ToString /@ values, "_"],
          "--config", manifest["ConfigStem"]};
        attemptDir = FileNameJoin[{jobDir, "attempts",
          "probe" <> ToString[probe] <> "-" <>
            StringReplace[CreateUUID[], "-" -> ""]}];
        CreateDirectory[attemptDir, CreateIntermediateDirectories -> True];
        Export[FileNameJoin[{attemptDir, "argv.wxf"}], command, "WXF"];
        run = RunSupervised[command, jobDir, timeout,
          FileNameJoin[{attemptDir, "fire7p"}]];
        Export[FileNameJoin[{attemptDir, "outcome.wxf"}], run, "WXF"];
        If[FailureQ[run] || Lookup[run, "ExitCode", $Failed] =!= 0 ||
            !nonemptyFileQ[expected],
          Return[failure["FIRE7p basis probe failed",
            <|"Probe" -> probe, "ExpectedTable" -> expected,
              "Outcome" -> run|>], Module]]];
      AppendTo[tables, expected];
      AppendTo[probeRecords, {probe, values, expected}],
      {probe, probes}];
    If[!runtimeIdentityMatchesQ[job, fireRoot],
      Scan[quarantineProbe[#[[3]]] &, probeRecords];
      Return[failure["FIRE runtime changed during basis probes",
        <|"JobDirectory" -> jobDir|>], Module]];
    groupValid = settingValidatorQ[settings, "BasisSetValidator", tables];
    If[!groupValid,
      Scan[quarantineProbe[#[[3]]] &, probeRecords];
      Return[failure["FIRE7 basis probes failed semantic validation",
        <|"Tables" -> tables|>], Module]];
    If[!AllTrue[probeRecords,
        writeProbeMarker[job, #[[3]], #[[1]], #[[2]]] =!= $Failed &],
      Return[failure["failed to atomically mark FIRE7 basis probes complete",
        <|"Tables" -> tables|>], Module]];
    <|"Success" -> True, "Job" -> job, "Manifest" -> manifest,
      "Tables" -> tables|>
  ]]
];

End[];
EndPackage[];
