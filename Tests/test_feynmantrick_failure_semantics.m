(* Process-free regression tests for mandatory FIRE/iteration failures. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"],
  {General::shdw, Symbol::shdw}];
FeynmanTrick`SetFTOption["Verbosity", 0];

passed = 0; failed = 0;
assert[label_, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

(* A failed mandatory level must prevent both its export and every later
   level.  Localized stubs keep this test independent of FIRE. *)
Module[{topology, ftData, result, computeCalls = {}, exportCalls = {}},
  topology = FeynmanTrick`FIREInterface`DefineTopology[
    "failure_semantics", {Global`l1}, {},
    {1 - Global`l1^2, 2 - Global`l1^2, 3 - Global`l1^2}, {}];
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}, {1, 3}}, {}];

  result = Block[{
      FeynmanTrick`FeynmanTrickIteration`ComputeLevelData,
      FeynmanTrick`FeynmanTrickIteration`ExportLevel},
    FeynmanTrick`FeynmanTrickIteration`ComputeLevelData[
        data_Association, level_Integer] := Module[{out = data},
      AppendTo[computeCalls, level];
      out["Levels"][level]["Computed"] = False;
      out["Levels"][level]["Error"] = "injected failure";
      out];
    FeynmanTrick`FeynmanTrickIteration`ExportLevel[
        _Association, level_Integer, ___] := AppendTo[exportCalls, level];
    FeynmanTrick`FeynmanTrickIteration`RunFullIteration[
      ftData, FileNameJoin[{$TemporaryDirectory, "unused_ft_failure_export"}]]
  ];

  assert["full_iteration_returns_failed_on_noncomputed_level",
    result === $Failed];
  assert["full_iteration_stops_at_first_failed_level",
    computeCalls === {1}];
  assert["full_iteration_does_not_export_failed_level",
    exportCalls === {}];
];

(* A stale/partial .tables file must never turn a nonzero FIRE exit into a
   successful result.  Stub both the process runner and the table reader;
   the reader must remain untouched. *)
Module[{dir, topology, partialTable, result, runCalls = 0, readCalls = 0},
  dir = FileNameJoin[{$TemporaryDirectory,
    "ft_fire_nonzero_" <> ToString[$ProcessID]}];
  If[DirectoryQ[dir], DeleteDirectory[dir, DeleteContents -> True]];
  CreateDirectory[dir, CreateIntermediateDirectories -> True];
  partialTable = FileNameJoin[{dir, "nonzero_stub.tables"}];
  Export[partialTable, "partial stale output", "Text"];

  topology = FeynmanTrick`FIREInterface`DefineTopology[
    "nonzero_stub", {Global`l1}, {}, {1 - Global`l1^2}, {}];
  topology["WorkDirectory"] = dir;
  topology["ProblemNumber"] = 71;
  topology["StartFileReady"] = True;

  result = Block[{FeynmanTrick`FIREInterface`Private`runFIRE6,
      FeynmanTrick`FIREInterface`Private`preparedTopologyCompatibleQ,
      FeynmanTrick`FIREInterface`Private`preparedStartHashes,
      FeynmanTrick`FIREInterface`Private`preparedRunnerRuntimeHashes,
      FeynmanTrick`FIREInterface`Private`ensureFIRELoaded,
      FIRE`Tables2Masters},
    FeynmanTrick`FIREInterface`Private`preparedTopologyCompatibleQ[___] := True;
    FeynmanTrick`FIREInterface`Private`preparedStartHashes[___] := <||>;
    FeynmanTrick`FIREInterface`Private`preparedRunnerRuntimeHashes[___] := <||>;
    FeynmanTrick`FIREInterface`Private`ensureFIRELoaded[] := True;
    FeynmanTrick`FIREInterface`Private`runFIRE6[___] :=
      (runCalls++; 134);
    FIRE`Tables2Masters[___] := (readCalls++; {{71, {1}}});
    FeynmanTrick`FIREInterface`FindBasis[topology]
  ];

  assert["find_basis_rejects_nonzero_fire_exit", result === $Failed];
  assert["find_basis_attempts_fire_once", runCalls === 1];
  assert["find_basis_ignores_partial_tables_after_failure", readCalls === 0];
  Quiet[DeleteDirectory[dir, DeleteContents -> True]];
];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
