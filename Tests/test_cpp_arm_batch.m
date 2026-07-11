(* One-kernel endpoint-arm recurrence batching.  This does not run the FT
   ladder: two independent prepared charts stand in for the lower/upper arm
   fronts, exercising the exact request capture, one native task pool,
   verified replay, memo population, strictness, and analytic regulators. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps; rho = Global`rho;
req = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TOrder" -> 8|>;
chart[name_, center_] := <|"ChartVar" -> t, "Center" -> center,
  "Scale" -> 1, "Radius" -> 2, "LocalRadius" -> 2,
  "Name" -> name, "Prescriptions" -> {}|>;
coeffs[fs_] := Flatten[
  Flatten[#["Sectors"] & /@ fs["Columns"]][[All, "Coeffs"]]];
clearHomogeneous[] := (
  DiffExp2`Solve`Private`$shCache = <||>;
  DiffExp2`Solve`Private`$shSysTag = None;);

oldThreads = Environment["DE2_CPP_THREADS"];
SetEnvironment["DE2_CPP_THREADS" -> "2"];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2, "RecurrenceBackend" -> "Cpp",
  "Variables" -> {}}]];

sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> {{1}}, "Variable" -> x|>]];
csLo = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart["batch-lo", 0]]];
csHi = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart["batch-hi", 1/3]]];

clearHomogeneous[];
sequential = catchDE2[
  DiffExp2`Solve`SolveHomogeneous[#, req] & /@ {csLo, csHi}];
sequentialCoeffs = If[FailureQ[sequential], {}, coeffs /@ sequential];

clearHomogeneous[];
batched = catchDE2[
  DiffExp2`Solve`PrewarmHomogeneousBatch[{csLo, csHi}, req]];
assert["cpp_arm_batch_returns_both_verified_chart_bases",
  ListQ[batched] && Length[batched] === 2 &&
    AllTrue[batched, AssociationQ[#] && KeyExistsQ[#, "Columns"] &]];
assert["cpp_arm_batch_populates_the_authoritative_homogeneous_cache",
  Length[DiffExp2`Solve`Private`$shCache] === 2 &&
    SameQ[batched,
      DiffExp2`Solve`SolveHomogeneous[#, req] & /@ {csLo, csHi}]];
batchCoeffs = If[FailureQ[batched], {}, coeffs /@ batched];
assert["cpp_arm_batch_matches_sequential_cpp_coefficients",
  Length[batchCoeffs] === 2 && Length[sequentialCoeffs] === 2 &&
    Max[Abs[N[Flatten[batchCoeffs - sequentialCoeffs], 70]]] < 10^-80];

(* The anchor chart appears in both endpoint plans.  Deduplication must use
   the same key as SolveHomogeneous and must not execute/cache it twice. *)
clearHomogeneous[];
duplicate = catchDE2[
  DiffExp2`Solve`PrewarmHomogeneousBatch[{csLo, csLo}, req]];
assert["cpp_arm_batch_deduplicates_shared_anchor_chart",
  ListQ[duplicate] && Length[duplicate] === 1 &&
    Length[DiffExp2`Solve`Private`$shCache] === 1];

(* An extra analytic regulator is retained over FLINT's exact rational
   coefficient field.  The native dispatcher deliberately serializes such
   a batch to one worker until allocator-level parallelism is certified. *)
regRate = (1 + rho)/(2 - rho);
regSys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> {{regRate}}, "Variable" -> x|>]];
regLo = catchDE2[DiffExp2`Solve`PrepareChart[
  regSys, chart["batch-reg-lo", 0]]];
regHi = catchDE2[DiffExp2`Solve`PrepareChart[
  regSys, chart["batch-reg-hi", 1/5]]];
regBatch = catchDE2[
  DiffExp2`Solve`PrewarmHomogeneousBatch[{regLo, regHi}, req]];
regCoeffs = If[FailureQ[regBatch], {},
  Table[regBatch[[i, "Columns", 1, "Sectors", 1, "Coeffs",
    1, All, 1]], {i, 2}]];
assert["cpp_arm_batch_preserves_exact_analytic_regulator_field",
  Length[regCoeffs] === 2 && FreeQ[regCoeffs, _?InexactNumberQ] &&
    And @@ Flatten[Table[
      Cancel[Together[regCoeffs[[i, n + 1]] - regRate^n/n!]] === 0,
      {i, 2}, {n, 0, req["TOrder"]}]]];

(* The batching seam is C++-only and must never become a hidden alternate
   solver selection. *)
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Wolfram"}]];
wrongBackend = catchDE2[
  DiffExp2`Solve`PrewarmHomogeneousBatch[{regLo, regHi}, req]];
assert["cpp_arm_batch_refuses_wolfram_fallback",
  FailureQ[wrongBackend] && wrongBackend["ID"] === "E6"];

If[!StringQ[oldThreads],
  SetEnvironment["DE2_CPP_THREADS" -> None],
  SetEnvironment["DE2_CPP_THREADS" -> oldThreads]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
