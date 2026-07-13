(* Focused, process-only contracts for FIRE 7 Classical output publication.

   Short fake executables stand in for FIRE.  The test deliberately plants a
   stale canonical table and has the isolated child emit unrelated .tables
   files, so process exit status alone cannot satisfy the contract. *)

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

readString[path_String] := If[FileExistsQ[path],
  Quiet[Check[ReadString[path], ""]], ""];

makeCase[root_String, stem_String, output_String,
    script_String] := Module[{dir, binary},
  dir = FileNameJoin[{root, stem}];
  CreateDirectory[dir, CreateIntermediateDirectories -> True];
  Export[FileNameJoin[{dir, stem <> ".config"}], StringJoin[
    "#threads 1\n",
    "#fthreads 1\n",
    "#variables x,d\n",
    "#start\n",
    "#problem 1 ", stem, ".start\n",
    "#integrals ", stem, ".m\n",
    "#output ", output, "\n"], "Text"];
  Export[FileNameJoin[{dir, stem <> ".start"}], "fake start\n", "Text"];
  Export[FileNameJoin[{dir, stem <> ".m"}], "{{1,{1}}}\n", "Text"];
  binary = FileNameJoin[{dir, "fake_fire7"}];
  Export[binary, script, "Text"];
  If[!makeExecutable[binary], Return[$Failed, Module]];
  <|"Directory" -> dir, "Stem" -> stem, "Output" -> output,
    "Canonical" -> FileNameJoin[{dir, output}], "Binary" -> binary|>
];

tempRoot = FileNameJoin[{$TemporaryDirectory,
  "ft_fire7_classical_exact_" <> ToString[$ProcessID] <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
If[DirectoryQ[tempRoot], DeleteDirectory[tempRoot, DeleteContents -> True]];
CreateDirectory[tempRoot, CreateIntermediateDirectories -> True];

oldBackend = Lookup[FeynmanTrick`Private`$FTConfig,
  "FIREBackend", "Modular"];
oldVerbosity = Lookup[FeynmanTrick`Private`$FTConfig, "Verbosity", 1];
oldTimeout = Lookup[FeynmanTrick`Private`$FTConfig,
  "FIRETimeoutSeconds", 600];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 10];

(* Exit zero without the exact #output artifact is a failure.  A stale file
   in the caller's directory must not turn that isolated run into success,
   and a differently named table emitted by the child must not be copied. *)
missingCase = makeCase[tempRoot, "missing", "declared_missing.tables",
  StringJoin[
    "#!/bin/sh\n",
    "printf 'unrelated child table\\n' > unrelated_missing.tables\n",
    "exit 0\n"]];
Export[missingCase["Canonical"], "stale canonical table\n", "Text"];
missingExit = FeynmanTrick`FIREInterface`Private`runFIRE6Once[
  missingCase["Binary"], missingCase["Directory"], missingCase["Stem"]];

assert["Classical exit zero without config-declared output is rejected",
  missingExit =!= 0];
assert["a stale canonical table is not published or accepted as fresh",
  missingExit =!= 0 &&
    (!FileExistsQ[missingCase["Canonical"]] ||
      readString[missingCase["Canonical"]] === "stale canonical table\n")];
assert["an unrelated failed-run table is not copied from isolation",
  !FileExistsQ[FileNameJoin[{missingCase["Directory"],
    "unrelated_missing.tables"}]]];

(* The Classical dispatcher must publish only the exact newly generated
   #output artifact.  Other .tables files in the isolated run are ignored. *)
successCase = makeCase[tempRoot, "success", "declared_success.tables",
  StringJoin[
    "#!/bin/sh\n",
    "printf 'fresh declared table\\n' > artifacts/declared_success.tables\n",
    "printf 'unrelated child table\\n' > unrelated_success.tables\n",
    "exit 0\n"]];
Export[successCase["Canonical"], "old canonical table\n", "Text"];
FeynmanTrick`SetFTOption["FIREBackend", "Classical"];
successExit = FeynmanTrick`FIREInterface`Private`runFIRE6[
  successCase["Binary"], successCase["Directory"], successCase["Stem"]];

assert["Classical dispatcher publishes the exact newly generated output",
  successExit === 0 &&
    readString[successCase["Canonical"]] === "fresh declared table\n"];
assert["successful Classical run does not copy unrelated tables",
  !FileExistsQ[FileNameJoin[{successCase["Directory"],
    "unrelated_success.tables"}]]];

FeynmanTrick`SetFTOption["FIREBackend", oldBackend];
FeynmanTrick`SetFTOption["Verbosity", oldVerbosity];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", oldTimeout];

If[DirectoryQ[tempRoot],
  Quiet[DeleteDirectory[tempRoot, DeleteContents -> True]]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
