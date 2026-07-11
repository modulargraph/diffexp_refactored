(* Real FIRE cold/warm parity microbenchmark for PersistentFIREStorage.

   STATUS: this is a frozen NO-GO reproducer.  FIRE 6.5.2 warm reuse fails
   exact round-2 parity; success would indicate that FIRE or the integration
   strategy changed and deserves a fresh audit.  The expected current exit is 1.

   The topology is a fixed Euclidean one-loop massive bubble.  Two disjoint
   request batches mimic successive basis-closure rounds.  ReductionCache is
   disabled so the second batch can only benefit from FIRE's persisted sector
   database, not from FeynmanTrick's exact in-memory integral cache.

   Optional environment:
     FT_FIRE_PATH                  FIRE6 installation root
     FT_FIRE_STORAGE_BENCH_ROOT    artifact directory
     FT_FIRE_STORAGE_BENCH_KEEP=1  preserve artifacts after success
     FT_FIRE_STORAGE_BENCH_EXPLICIT_ENVELOPE=1
                                   diagnostic: explicitly union both rounds
*)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

envOrDefault[name_String, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]
];

firePath = ExpandFileName[envOrDefault["FT_FIRE_PATH",
  FileNameJoin[{repoRoot, "Dependencies", "fire", "FIRE6"}]]];
benchRoot = ExpandFileName[envOrDefault["FT_FIRE_STORAGE_BENCH_ROOT",
  FileNameJoin[{$TemporaryDirectory,
    "ft-fire-storage-real-" <> ToString[$ProcessID] <> "-" <>
      StringReplace[CreateUUID[], "-" -> ""]}]]];
keepArtifacts = envOrDefault["FT_FIRE_STORAGE_BENCH_KEEP", "0"] === "1";
explicitWarmEnvelope =
  envOrDefault["FT_FIRE_STORAGE_BENCH_EXPLICIT_ENVELOPE", "0"] === "1";
offDir = FileNameJoin[{benchRoot, "off"}];
onDir = FileNameJoin[{benchRoot, "on"}];
storageDir = FileNameJoin[{benchRoot, "storage"}];
evidenceDir = FileNameJoin[{benchRoot, "evidence"}];

If[DirectoryQ[benchRoot], DeleteDirectory[benchRoot, DeleteContents -> True]];
CreateDirectory[offDir, CreateIntermediateDirectories -> True];
CreateDirectory[onDir, CreateIntermediateDirectories -> True];
CreateDirectory[evidenceDir, CreateIntermediateDirectories -> True];

assertions = <||>;
assert[label_String, condition_] := (assertions[label] = TrueQ[condition]);
readText[path_String] := If[FileExistsQ[path], Quiet[ReadString[path]], ""];
storageSnapshot[] := Module[{files, databaseFiles},
  files = Select[FileNames["*", storageDir, Infinity], FileType[#] === File &];
  databaseFiles = Select[files,
    MemberQ[{"tmp", "kch"}, ToLowerCase[FileExtension[#]]] &];
  <|"Files" -> Length[files],
    "DatabaseFiles" -> Length[databaseFiles],
    "Bytes" -> Total[Quiet[FileByteCount /@ files]]|>
];
currentPointer[] := Module[{pointers = FileNames["CURRENT", storageDir, Infinity]},
  If[Length[pointers] === 1, StringTrim[readText[First[pointers]]],
    "NotAvailable"]
];
generationPointerQ[value_] := StringQ[value] &&
  StringMatchQ[value, RegularExpression["gen-[0-9a-f]{32}"]];
canonicalMasters[result_] := If[AssociationQ[result],
  Sort[DeleteDuplicates[Lookup[result, "Masters", {}]]], "NotAvailable"];
sameExactReductions[left_Association, right_Association] := Module[{keys},
  keys = Sort[Keys[left]];
  keys === Sort[Keys[right]] && AllTrue[keys,
    Together[Expand[left[#] - right[#]]] === 0 &]
];
sameExactReductions[_, _] := False;
reductionRows[result_] := If[AssociationQ[result],
  ({#, ToString[result["Reductions"][#], InputForm]} &) /@
    Keys[result["Reductions"]], {}];
restrictDetailedResult[result_, keys_List] := If[AssociationQ[result],
  Join[result, <|"Reductions" -> Association[
    (# -> result["Reductions"][#]) & /@ keys]|>], result];
fireReportedSeconds[log_String] := Module[{matches},
  matches = StringCases[log,
    RegularExpression["Total time: ([0-9.]+)"] -> "$1"];
  If[matches === {}, "NotReported", ToExpression[Last[matches]]]
];
snapshotRunArtifacts[dir_String, label_String] := Scan[
  Function[path, If[FileExistsQ[path],
    CopyFile[path, FileNameJoin[{evidenceDir,
      label <> "-" <> FileNameTake[path]}], OverwriteTarget -> True]]],
  FileNameJoin[{dir, #}] & /@ {
    "fire_stdout.log", "persistent_bubble_reduce.config",
    "persistent_bubble_reduce.m", "persistent_bubble_reduce.tables"}
];

If[!FileExistsQ[FileNameJoin[{firePath, "bin", "FIRE6"}]],
  Print["FIREPERSIST ERROR missing FIRE6 binary under ", firePath];
  Quit[2]
];

FeynmanTrick`SetFTOption["FIREPath", firePath];
FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 5];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["ReductionCache", False];
FeynmanTrick`SetFTOption["PersistentFIREStorage", False];

topology = FeynmanTrick`FIREInterface`DefineTopology[
  "persistent_bubble",
  {Global`l},
  {Global`p},
  {1 - Global`l^2, 2 - (Global`l + Global`p)^2},
  {Global`p^2 -> -3}
];
offTopology = FeynmanTrick`FIREInterface`SetupFIRE[topology, offDir];
If[offTopology === $Failed,
  Print["FIREPERSIST ERROR SetupFIRE failed"];
  Quit[2]
];

(* Reuse the exact .start bytes and problem number in the persistence-on arm,
   so the only experimental variable is FIRE's storage configuration. *)
CopyFile[
  FileNameJoin[{offDir, offTopology["Name"] <> ".start"}],
  FileNameJoin[{onDir, offTopology["Name"] <> ".start"}],
  OverwriteTarget -> True
];
onTopology = offTopology;
onTopology["WorkDirectory"] = onDir;

round1 = {{1, 1}, {2, 1}, {1, 2}};
round2 = {{3, 1}, {2, 2}, {1, 3}};
warmInput = If[explicitWarmEnvelope,
  DeleteDuplicates[Join[round1, round2]], round2];

{offTime1, offResult1} = AbsoluteTiming[
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[offTopology, round1]];
offLog1 = readText[FileNameJoin[{offDir, "fire_stdout.log"}]];
snapshotRunArtifacts[offDir, "off-round1"];
{offTime2, offResult2} = AbsoluteTiming[
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[offTopology, round2]];
offLog2 = readText[FileNameJoin[{offDir, "fire_stdout.log"}]];
snapshotRunArtifacts[offDir, "off-round2"];

FeynmanTrick`SetFTOption["PersistentFIREStorage", True];
FeynmanTrick`SetFTOption["FIREStorageDirectory", storageDir];
{onTime1, onResult1} = AbsoluteTiming[
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[onTopology, round1]];
onLog1 = readText[FileNameJoin[{onDir, "fire_stdout.log"}]];
snapshotRunArtifacts[onDir, "on-cold-round1"];
coldStorage = storageSnapshot[];
coldPointer = currentPointer[];
{onTime2, onFullResult2} = AbsoluteTiming[
  FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[onTopology, warmInput]];
onResult2 = restrictDetailedResult[onFullResult2, round2];
onLog2 = readText[FileNameJoin[{onDir, "fire_stdout.log"}]];
snapshotRunArtifacts[onDir, "on-warm-round2"];
warmStorage = storageSnapshot[];
warmPointer = currentPointer[];

assert["all reductions completed",
  AllTrue[{offResult1, offResult2, onResult1, onResult2}, AssociationQ]];
assert["round 1 reductions agree exactly",
  AssociationQ[offResult1] && AssociationQ[onResult1] &&
  sameExactReductions[offResult1["Reductions"], onResult1["Reductions"]]];
assert["round 2 reductions agree exactly",
  AssociationQ[offResult2] && AssociationQ[onResult2] &&
  sameExactReductions[offResult2["Reductions"], onResult2["Reductions"]]];
assert["round 1 master sets agree exactly",
  AssociationQ[offResult1] && AssociationQ[onResult1] &&
  canonicalMasters[offResult1] === canonicalMasters[onResult1]];
assert["round 2 master sets agree exactly",
  AssociationQ[offResult2] && AssociationQ[onResult2] &&
  canonicalMasters[offResult2] === canonicalMasters[onResult2]];
assert["cold run published nonempty storage",
  coldStorage["DatabaseFiles"] > 0 && coldStorage["Bytes"] > 0 &&
  generationPointerQ[coldPointer]];
assert["warm run consumed existing storage",
  StringContainsQ[onLog2, "Storage directory located"] &&
  StringContainsQ[onLog2, "Copying file"]];
assert["second successful run advanced CURRENT",
  generationPointerQ[coldPointer] && generationPointerQ[warmPointer] &&
  coldPointer =!= warmPointer];
assert["persistence-off configs did not use storage",
  !StringContainsQ[offLog1 <> offLog2, "Storage path:"]];

summary = <|
  "Topology" -> "fixed Euclidean one-loop massive bubble",
  "Round1" -> round1,
  "Round2" -> round2,
  "WarmInput" -> warmInput,
  "Seconds" -> <|
    "OffRound1" -> offTime1,
    "OffRound2" -> offTime2,
    "OnColdRound1" -> onTime1,
    "OnWarmRound2" -> onTime2|>,
  "FIRESeconds" -> <|
    "OffRound1" -> fireReportedSeconds[offLog1],
    "OffRound2" -> fireReportedSeconds[offLog2],
    "OnColdRound1" -> fireReportedSeconds[onLog1],
    "OnWarmRound2" -> fireReportedSeconds[onLog2]|>,
  "MasterSets" -> <|
    "OffRound1" -> canonicalMasters[offResult1],
    "OffRound2" -> canonicalMasters[offResult2],
    "OnColdRound1" -> canonicalMasters[onResult1],
    "OnWarmRound2" -> canonicalMasters[onResult2]|>,
  "Reductions" -> <|
    "OffRound1" -> reductionRows[offResult1],
    "OffRound2" -> reductionRows[offResult2],
    "OnColdRound1" -> reductionRows[onResult1],
    "OnWarmRound2" -> reductionRows[onResult2]|>,
  "StorageAfterCold" -> coldStorage,
  "StorageAfterWarm" -> warmStorage,
  "ColdPointer" -> coldPointer,
  "WarmPointer" -> warmPointer,
  "WarmLogStorageLocated" -> StringContainsQ[onLog2,
    "Storage directory located"],
  "WarmLogCopiedFiles" -> Length[StringCases[onLog2, "Copying file"]],
  "Assertions" -> assertions,
  "Artifacts" -> benchRoot
|>;
Print["FIREPERSIST ", ExportString[summary, "RawJSON", "Compact" -> True]];
Scan[Print[If[TrueQ[Last[#]], "  PASS: ", "  FAIL: "], First[#]] &,
  Normal[assertions]];

success = AllTrue[Values[assertions], TrueQ];
If[success && !keepArtifacts,
  Quiet[DeleteDirectory[benchRoot, DeleteContents -> True]]
];
If[success, Quit[0], Quit[1]];
