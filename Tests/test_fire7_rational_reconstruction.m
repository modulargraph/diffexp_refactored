(* Focused contract tests for FIRE7_MPI -> standalone rational reconstruction.

   FIRE is never invoked.  Short fake executables model the exact artifact
   boundary so stale cache files, exit codes, and reconstructed output names
   can be checked deterministically. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick", "FIRE7Runner.m"}]],
  {General::shdw, Symbol::shdw}];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

makeExecutable[path_String] :=
  Lookup[RunProcess[{"/bin/chmod", "+x", path}], "ExitCode", -1] === 0;

nonemptyFileQ[path_String] := FileExistsQ[path] &&
  Quiet[Check[FileByteCount[path] > 0, False]];

fileSHA256[path_String] := IntegerString[FileHash[path, "SHA256"], 16, 64];

readString[path_String] := If[FileExistsQ[path],
  Quiet[Check[ReadString[path], ""]], ""];

readArgv[path_String] := Module[{text = readString[path]},
  If[StringTrim[text] === "", {},
    StringSplit[StringTrim[text], {"\r\n", "\n", "\r"}]]
];

invocationCount[path_String] := Length[
  Select[StringSplit[readString[path], {"\r\n", "\n", "\r"}],
    # === "invoked" &]];

resolvedFromJob[path_String, jobDir_String] := ExpandFileName[
  If[StringStartsQ[path, "/"], path, FileNameJoin[{jobDir, path}]]];

tempRoot = FileNameJoin[{$TemporaryDirectory,
  "ft_fire7_rational_" <> ToString[$ProcessID] <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
If[DirectoryQ[tempRoot], DeleteDirectory[tempRoot, DeleteContents -> True]];
CreateDirectory[tempRoot, CreateIntermediateDirectories -> True];

mpiScript = StringJoin[
  "#!/bin/sh\n",
  "printf 'invoked\\n' >> mpi.invocations\n",
  "mkdir -p artifacts\n",
  "printf 'prime table 1\\n' > artifacts/contract_x_m_d_1.tables\n",
  "printf 'prime table 2\\n' > artifacts/contract_x_m_d_2.tables\n",
  "exit 0\n"];

successfulReconstruct = StringJoin[
  "#!/bin/sh\n",
  "printf 'invoked\\n' >> reconstruct.invocations\n",
  "printf '%s\\n' \"$@\" > reconstruct.argv\n",
  "mkdir -p \"$(dirname \"$5\")\"\n",
  "printf 'fresh reconstructed table\\n' > \"$5\"\n",
  "exit 0\n"];

failedReconstruct = StringJoin[
  "#!/bin/sh\n",
  "printf 'invoked\\n' >> reconstruct.invocations\n",
  "printf '%s\\n' \"$@\" > reconstruct.argv\n",
  "exit 9\n"];

missingReconstruct = StringJoin[
  "#!/bin/sh\n",
  "printf 'invoked\\n' >> reconstruct.invocations\n",
  "printf '%s\\n' \"$@\" > reconstruct.argv\n",
  "exit 0\n"];

wrongReconstruct = StringJoin[
  "#!/bin/sh\n",
  "printf 'invoked\\n' >> reconstruct.invocations\n",
  "printf '%s\\n' \"$@\" > reconstruct.argv\n",
  "mkdir -p artifacts\n",
  "printf 'wrong reconstructed table\\n' > artifacts/contract_x_m_d_wrong.tables\n",
  "exit 0\n"];

makeCase[name_String, reconstructScript_String] := Module[
  {caseRoot, inputDir, fireRoot, binDir, cacheDir, mpi, reconstruct,
   executablePaths, manifest, settings},
  caseRoot = FileNameJoin[{tempRoot, name}];
  inputDir = FileNameJoin[{caseRoot, "input"}];
  fireRoot = FileNameJoin[{caseRoot, "fire7"}];
  binDir = FileNameJoin[{fireRoot, "bin"}];
  cacheDir = FileNameJoin[{caseRoot, "cache"}];
  CreateDirectory[inputDir, CreateIntermediateDirectories -> True];
  CreateDirectory[binDir, CreateIntermediateDirectories -> True];
  Export[FileNameJoin[{inputDir, "contract.config"}], StringJoin[
    "#threads 1\n",
    "#fthreads 1\n",
    "#variables x,m,d\n",
    "#start\n",
    "#problem 1 contract.start\n",
    "#integrals contract.m\n",
    "#output contract.tables\n"], "Text"];
  Export[FileNameJoin[{inputDir, "contract.m"}],
    "{{1,{2,1}}}\n", "Text"];
  Export[FileNameJoin[{inputDir, "contract.start"}],
    "fake start\n", "Text"];
  Export[FileNameJoin[{fireRoot, "FIRE7.m"}],
    "(* fake FIRE7 source *)\n", "Text"];
  Export[FileNameJoin[{fireRoot, "paths.inc"}],
    "RESULTING_MPRIME_COUNT=16\n", "Text"];
  Do[Export[FileNameJoin[{binDir, executable}],
      "#!/bin/sh\nexit 0\n", "Text"],
    {executable, {"FIRE7_MPI", "FIRE7p", "FIRE7mp"}}];
  mpi = FileNameJoin[{caseRoot, "fake_mpirun"}];
  reconstruct = FileNameJoin[{binDir, "reconstruct"}];
  Export[mpi, mpiScript, "Text"];
  Export[reconstruct, reconstructScript, "Text"];
  executablePaths = Join[{mpi, reconstruct},
    FileNameJoin[{binDir, #}] & /@ {"FIRE7_MPI", "FIRE7p", "FIRE7mp"}];
  If[!AllTrue[executablePaths, makeExecutable],
    Return[$Failed, Module]];
  manifest = FeynmanTrick`FIRE7Runner`ParseConfig[inputDir, "contract"];
  settings = <|
    "MPIExecutable" -> mpi,
    "Workers" -> 2,
    "Calc" -> "flint",
    "UseMultiprime" -> True,
    "PrimeLimit" -> 7,
    "DimensionSeparated" -> False,
    "MultiprimeWidth" -> 16,
    "TimeoutSeconds" -> 10,
    "CacheDirectory" -> cacheDir
  |>;
  <|"Root" -> caseRoot, "InputDirectory" -> inputDir,
    "FireRoot" -> fireRoot, "CacheDirectory" -> cacheDir,
    "Manifest" -> manifest, "Settings" -> settings|>
];

jobLayout[case_Association] := Module[{job, runManifest, final},
  job = FeynmanTrick`FIRE7Runner`Private`modularJob[
    "Reduction", case["FireRoot"], case["InputDirectory"],
    case["Manifest"], case["Settings"]];
  runManifest = FeynmanTrick`FIRE7Runner`Private`jobManifest[
    case["Manifest"]];
  final = FeynmanTrick`FIRE7Runner`ExpectedReconstructedFile[
    job["Directory"], runManifest];
  <|"Job" -> job, "JobDirectory" -> job["Directory"],
    "Final" -> final,
    "Canonical" -> FileNameJoin[{case["InputDirectory"],
      case["Manifest"]["OutputFile"]}],
    "Complete" -> FileNameJoin[{job["Directory"], "complete.wxf"}]|>
];

runCase[case_Association] := FeynmanTrick`FIRE7Runner`RunModular[
  case["FireRoot"], case["InputDirectory"], "contract",
  case["Settings"]];

(* ---------------------------------------------------------------------- *)
(* A zero-exit MPI run with only prime tables must hand off to the exact  *)
(* official standalone rational reconstruction command.                  *)
(* ---------------------------------------------------------------------- *)

successCase = makeCase["success", successfulReconstruct];
successResult = runCase[successCase];
successLayout = jobLayout[successCase];
successArgs = readArgv[FileNameJoin[{successLayout["JobDirectory"],
  "reconstruct.argv"}]];
successComplete = Quiet[Check[Import[successLayout["Complete"], "WXF"],
  $Failed]];

assert["MPI prime tables are finalized by standalone reconstruction",
  AssociationQ[successResult] && TrueQ[successResult["Success"]] &&
    !TrueQ[successResult["Reused"]] &&
    readString[successLayout["Final"]] ===
      "fresh reconstructed table\n" &&
    readString[successLayout["Canonical"]] ===
      "fresh reconstructed table\n"];

assert["standalone reconstruction uses FIRE7's official rational argv",
  Length[successArgs] === 6 &&
    Take[successArgs, 4] === {"--method", "rational", "--calc", "flint"} &&
    resolvedFromJob[successArgs[[5]], successLayout["JobDirectory"]] ===
      ExpandFileName[successLayout["Final"]] &&
    Last[successArgs] === "7"];

assert["completion marker binds reuse to the exact final table hash",
  AssociationQ[successComplete] &&
    successComplete["JobID"] === successLayout["Job"]["ID"] &&
    successComplete["FinalTable"] === FileNameTake[successLayout["Final"]] &&
    successComplete["FinalSHA256"] === fileSHA256[successLayout["Final"]]];

successReuse = runCase[successCase];
assert["a validated completed reconstruction is reused without new workers",
  AssociationQ[successReuse] && TrueQ[successReuse["Success"]] &&
    TrueQ[successReuse["Reused"]] &&
    invocationCount[FileNameJoin[{successLayout["JobDirectory"],
      "mpi.invocations"}]] === 1 &&
    invocationCount[FileNameJoin[{successLayout["JobDirectory"],
      "reconstruct.invocations"}]] === 1];

(* ---------------------------------------------------------------------- *)
(* Reconstruction failure and a zero exit without the exact output are   *)
(* both hard failures; neither publishes nor marks the job complete.      *)
(* ---------------------------------------------------------------------- *)

failureCase = makeCase["reconstruct-failure", failedReconstruct];
failureResult = runCase[failureCase];
failureLayout = jobLayout[failureCase];
assert["nonzero standalone reconstruction fails the modular job",
  FailureQ[failureResult] &&
    invocationCount[FileNameJoin[{failureLayout["JobDirectory"],
      "reconstruct.invocations"}]] === 1 &&
    !FileExistsQ[failureLayout["Canonical"]] &&
    !FileExistsQ[failureLayout["Complete"]]];

wrongCase = makeCase["wrong-output", wrongReconstruct];
Export[FileNameJoin[{wrongCase["InputDirectory"], "contract.tables"}],
  "stale canonical table\n", "Text"];
wrongResult = runCase[wrongCase];
wrongLayout = jobLayout[wrongCase];
assert["exit zero with only a wrong reconstructed filename is rejected",
  FailureQ[wrongResult] &&
    invocationCount[FileNameJoin[{wrongLayout["JobDirectory"],
      "reconstruct.invocations"}]] === 1 &&
    !nonemptyFileQ[wrongLayout["Final"]] &&
    readString[wrongLayout["Canonical"]] === "stale canonical table\n" &&
    !FileExistsQ[wrongLayout["Complete"]]];

(* ---------------------------------------------------------------------- *)
(* A cache final without a valid marker is not evidence of completion.   *)
(* In particular, an old exact-name file cannot make a no-output          *)
(* reconstruction appear successful.                                    *)
(* ---------------------------------------------------------------------- *)

staleCase = makeCase["stale-unmarked", missingReconstruct];
staleLayout = jobLayout[staleCase];
CreateDirectory[DirectoryName[staleLayout["Final"]],
  CreateIntermediateDirectories -> True];
Export[staleLayout["Final"], "stale cached table\n", "Text"];
staleResult = runCase[staleCase];
assert["an unmarked stale exact-name table is neither reused nor accepted",
  FailureQ[staleResult] &&
    invocationCount[FileNameJoin[{staleLayout["JobDirectory"],
      "mpi.invocations"}]] === 1 &&
    invocationCount[FileNameJoin[{staleLayout["JobDirectory"],
      "reconstruct.invocations"}]] === 1 &&
    !FileExistsQ[staleLayout["Canonical"]] &&
    !FileExistsQ[staleLayout["Complete"]]];

(* A marker with the right job ID and filename but the wrong content hash
   must also force a fresh MPI/reconstruction attempt. *)
badHashCase = makeCase["bad-marker-hash", successfulReconstruct];
badHashLayout = jobLayout[badHashCase];
CreateDirectory[DirectoryName[badHashLayout["Final"]],
  CreateIntermediateDirectories -> True];
Export[badHashLayout["Final"], "unverified cached table\n", "Text"];
Export[badHashLayout["Complete"], <|
  "Schema" -> FeynmanTrick`FIRE7Runner`Private`$jobSchema <> "/complete",
  "JobID" -> badHashLayout["Job"]["ID"],
  "FinalTable" -> FileNameTake[badHashLayout["Final"]],
  "FinalSHA256" -> StringRepeat["0", 64],
  "Completed" -> DateString["ISODateTime"]|>, "WXF"];
badHashResult = runCase[badHashCase];
assert["a completion marker with the wrong table hash is not reused",
  AssociationQ[badHashResult] && TrueQ[badHashResult["Success"]] &&
    !TrueQ[badHashResult["Reused"]] &&
    invocationCount[FileNameJoin[{badHashLayout["JobDirectory"],
      "mpi.invocations"}]] === 1 &&
    invocationCount[FileNameJoin[{badHashLayout["JobDirectory"],
      "reconstruct.invocations"}]] === 1 &&
    readString[badHashLayout["Canonical"]] ===
      "fresh reconstructed table\n"];

If[DirectoryQ[tempRoot],
  Quiet[DeleteDirectory[tempRoot, DeleteContents -> True]]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
