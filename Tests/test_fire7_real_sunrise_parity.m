(* Real FIRE 7 sunrise smoke/parity test.  Unlike the focused contract suites
   this invokes the locally built FIRE7/FIRE7p/FIRE7mp/FIRE7_MPI/reconstruct
   binaries.  Keep this test opt-in: it launches both exact and modular FIRE. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Quiet[Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]],
  {General::shdw, Symbol::shdw}];

fireRoot = FileNameJoin[{repoRoot, "Dependencies", "fire", "FIRE7", "FIRE7"}];
required = Join[{FileNameJoin[{fireRoot, "FIRE7.m"}]},
  FileNameJoin[{fireRoot, "bin", #}] & /@
    {"FIRE7", "FIRE7p", "FIRE7mp", "FIRE7_MPI", "reconstruct"}];
If[!AllTrue[required, FileExistsQ],
  Print["SKIP: local FIRE 7 build is unavailable."];
  Quit[0]];

workRoot = FileNameJoin[{$TemporaryDirectory,
  "ft_fire7_real_sunrise_" <> ToString[$ProcessID] <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""]}];
CreateDirectory[workRoot, CreateIntermediateDirectories -> True];

FeynmanTrick`SetFTOption["FIREPath", fireRoot];
FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 1];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 600];
FeynmanTrick`SetFTOption["FIREModularWorkers", 2];
FeynmanTrick`SetFTOption["FIREUseMultiprime", True];
FeynmanTrick`SetFTOption["FIREPrimeLimit", 127];
FeynmanTrick`SetFTOption["FIREMultiprimeWidth", 16];
FeynmanTrick`SetFTOption["FIREModularCacheDirectory",
  FileNameJoin[{workRoot, "modular-cache"}]];

target = {3, 1, 1};

makeSunrise[] := FeynmanTrick`FIREInterface`DefineTopology[
  "sunrise", {Global`l1, Global`l2}, {Global`p},
  {
    1 - Global`l1^2,
    1 - Global`l2^2,
    1 - (-Global`l1 - Global`l2 + Global`p)^2
  },
  {Global`p^2 -> -1}];

runBackend[backend_String] := Module[
  {topology, setup, basis, detailed, seconds, reduction},
  FeynmanTrick`SetFTOption["FIREBackend", backend];
  If[FeynmanTrick`FIREInterface`ClearFIREState[] === $Failed,
    Return[$Failed, Module]];
  topology = makeSunrise[];
  {seconds, setup} = AbsoluteTiming[
    FeynmanTrick`FIREInterface`SetupFIRE[topology,
      FileNameJoin[{workRoot, ToLowerCase[backend]}]]];
  If[!AssociationQ[setup], Return[$Failed, Module]];
  basis = FeynmanTrick`FIREInterface`FindBasis[setup];
  If[!AssociationQ[basis], Return[$Failed, Module]];
  detailed = FeynmanTrick`FIREInterface`ReduceIntegralsDetailed[
    basis, {target}];
  If[!AssociationQ[detailed], Return[$Failed, Module]];
  reduction = detailed["Reductions"][target];
  <|
    "Backend" -> backend,
    "SetupSeconds" -> seconds,
    "Masters" -> Sort[basis["Masters"]],
    "ReductionMasters" -> Sort[detailed["Masters"]],
    "TargetIsNonMaster" -> !MemberQ[basis["Masters"], target],
    "ReductionIsExact" -> FreeQ[reduction, _Real],
    "Reduction" -> reduction
  |>
];

classical = runBackend["Classical"];
modular = runBackend["Modular"];
ok = AssociationQ[classical] && AssociationQ[modular] &&
  TrueQ[classical["TargetIsNonMaster"]] &&
  TrueQ[modular["TargetIsNonMaster"]] &&
  TrueQ[classical["ReductionIsExact"]] &&
  TrueQ[modular["ReductionIsExact"]] &&
  classical["Masters"] === modular["Masters"] &&
  classical["ReductionMasters"] === modular["ReductionMasters"] &&
  TrueQ[Together[classical["Reduction"] - modular["Reduction"]] === 0];

Print["TARGET    ", InputForm[target]];
Print["CLASSICAL ", InputForm[classical]];
Print["MODULAR   ", InputForm[modular]];
Print[If[ok, "PASS", "FAIL"],
  ": real FIRE7 Classical/Modular sunrise parity"];

If[DirectoryQ[workRoot],
  Quiet[DeleteDirectory[workRoot, DeleteContents -> True]]];
If[ok, Quit[0], Quit[1]];
