(* Once-per-exact-system SCC structure cache contracts. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2, "RecurrenceBackend" -> "Wolfram",
  "Variables" -> {Global`rho}}]];

x = Global`x; t = Global`t; eps = Global`eps; rho = Global`rho;
chart[center_, scale_, name_] := <|"ChartVar" -> t,
  "Center" -> center, "Scale" -> scale, "Radius" -> 1,
  "LocalRadius" -> 1, "Name" -> name, "Prescriptions" -> {},
  "UseSCCSkeleton" -> True|>;

(* The off-diagonal edges cannot be discovered by specializing eps or rho:
   both remain nonzero elements of the exact coefficient field. *)
sys = <|"Matrix" -> {
    {(1 + x)/(2 - x), 0, 0, 0},
    {eps/(3 + x), (2 + x)/(4 - x), 0, 0},
    {0, rho/(5 + x), (3 + x)/(6 - x), 1/(7 + x)},
    {0, 0, 1/(8 + x), (4 + x)/(9 - x)}},
  "Variable" -> x|>;

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`Solve`Private`$exactSCCStructureCoreCalls = 0;
charts = {
  chart[0, 1, "scc-cache-a"],
  chart[1/5, 2/3, "scc-cache-b"],
  chart[-2/7, -3/5, "scc-cache-c"]};
systems = catchDE2[DiffExp2`Solve`PrepareChart[sys, #]] & /@ charts;
cachedCoreCalls = DiffExp2`Solve`Private`$exactSCCStructureCoreCalls;

assert["scc_structure_core_once_across_affine_charts",
  AllTrue[systems, !FailureQ[#] &] && cachedCoreCalls === 1 &&
    Length[DiffExp2`Solve`Private`$exactSCCStructureCache] === 1];

assert["scc_structure_preserves_epsilon_and_regulator_edges",
  AllTrue[systems, Function[cs,
    cs["IntegrationSequence", "Components"] === {{1}, {2}, {3, 4}} &&
    cs["IntegrationSequence", "CondensationEdges"] === {{1, 2}, {2, 3}} &&
    MemberQ[cs["IntegrationSequence", "DependencyEdges"], {1, 2}] &&
    MemberQ[cs["IntegrationSequence", "DependencyEdges"], {2, 3}]]]];

(* Direct chart-local analysis is the uncached reference implementation.
   The complete certificate, including the chart-specific MatrixHash, must
   be byte-for-byte identical at every center and signed nonzero scale. *)
localStructures = catchDE2[
    DiffExp2`Solve`Private`exactSCCStructure[# ["ThetaOriginal"]]] & /@
  systems;
assert["scc_cached_structure_matches_chart_local_reference",
  And @@ MapThread[#1["IntegrationSequence"] === #2 &,
    {systems, localStructures}]];

(* The independent-variable identity is part of the exact system identity,
   even when the stored matrix expression happens to be syntactically the
   same. *)
y = Global`y;
callsBeforeDistinctVariable =
  DiffExp2`Solve`Private`$exactSCCStructureCoreCalls;
distinctVariable = catchDE2[DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> sys["Matrix"], "Variable" -> y|>,
  chart[0, 1, "scc-cache-distinct-variable"]]];
assert["scc_structure_cache_keys_full_system_identity",
  !FailureQ[distinctVariable] &&
    DiffExp2`Solve`Private`$exactSCCStructureCoreCalls ===
      callsBeforeDistinctVariable + 1 &&
    Length[DiffExp2`Solve`Private`$exactSCCStructureCache] === 2];

(* ClearSolveCaches owns this cache's lifetime. *)
DiffExp2`Solve`ClearSolveCaches[];
callsBeforeRebuild = DiffExp2`Solve`Private`$exactSCCStructureCoreCalls;
rebuilt = catchDE2[DiffExp2`Solve`PrepareChart[
  sys, chart[1/9, 4/7, "scc-cache-after-clear"]]];
assert["scc_structure_cache_cleared_and_rebuilt",
  !FailureQ[rebuilt] &&
    DiffExp2`Solve`Private`$exactSCCStructureCoreCalls ===
      callsBeforeRebuild + 1 &&
    Length[DiffExp2`Solve`Private`$exactSCCStructureCache] === 1];

Print["SCC structure cache tests: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1]];
