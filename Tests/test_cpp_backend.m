(* Focused parity and strictness tests for the optional FLINT/Arb C++
   recurrence backend. Build first with:
     cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
     cmake --build build -j4
   Set DE2_REQUIRE_CPP=1 to make an unavailable library a hard failure. *)

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
req[min_, max_, n_] := <|"EpsWindow" -> <|"Min" -> min,
  "CompleteMax" -> max|>, "TOrder" -> n|>;
chart[name_] := <|"ChartVar" -> t, "Center" -> 0, "Radius" -> 1,
  "Name" -> name|>;
tags[ls_] := {#["a"], #["b"], #["p"]} & /@ ls["Sectors"];

info = DiffExp2`CppBackend`BackendInformation[];
assert["cpp_backend_librarylink_flint_available",
  AssociationQ[info] && info["schema"] === 1 &&
  StringQ[info["flint"]]];

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100}]];
assert["cpp_backend_default_preserves_wolfram_reference",
  DiffExp2`Config`CFG["RecurrenceBackend"] === "Wolfram"];

(* Ordinary f'=f: exercise exact structural tags/windows and numerical
   coefficient parity through compiled recurrence+V assembly+ODE proof. *)
DiffExp2`Solve`ClearSolveCaches[];
csExp = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{1}}, "Variable" -> x|>, chart["cpp-exp"]];
rExpW = catchDE2[DiffExp2`Solve`SolveChart[csExp, req[0, 2, 10]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Cpp"}]];
rExpC = catchDE2[DiffExp2`Solve`SolveChart[csExp, req[0, 2, 10]]];
expCoeffs = If[FailureQ[rExpC], {},
  rExpC["Basis", "Columns"][[1, "Sectors", 1, "Coeffs", 1, All, 1]]];
assert["cpp_regular_exponential_parity",
  !FailureQ[rExpW] && !FailureQ[rExpC] &&
  rExpW["Basis", "Columns"][[1, "EpsWindow"]] ===
    rExpC["Basis", "Columns"][[1, "EpsWindow"]] &&
  tags[rExpC["Basis", "Columns"][[1]]] === {{0, 0, 0}} &&
  Max[Abs[N[expCoeffs - Table[1/n!, {n, 0, 10}], 70]]] < 10^-70];
assert["cpp_backend_mode_is_part_of_homogeneous_cache_key",
  Length[DiffExp2`Solve`Private`$shCache] === 2];

(* Inexact Cauchy data may legitimately carry less than the configured
   2x-WP input budget. The bridge must serialize the reliable midpoint plus
   its uncertainty (never Indeterminate padding), and the decoded result must
   retain that honest Accuracy rather than being blindly stamped WP+20. *)
lowPi = N[Pi, 50];
rLowPi = catchDE2[DiffExp2`Solve`SolveValueRegular[csExp, req[0, 0, 4],
  {DiffExp2`EpsSeries`ESNew[0, {lowPi}]}]];
lowPiCoeffs = If[FailureQ[rLowPi], {},
  rLowPi[["Sectors", 1, "Coeffs", 1, All, 1]]];
assert["cpp_lower_precision_handoff_preserves_uncertainty_without_indeterminate",
  !FailureQ[rLowPi] &&
  Max[Abs[N[lowPiCoeffs - Table[lowPi/n!, {n, 0, 4}], 40]]] < 10^-40 &&
  40 < Accuracy[First[lowPiCoeffs]] < 60];

rMachineComplex = catchDE2[DiffExp2`Solve`SolveValueRegular[
  csExp, req[0, 0, 4],
  {DiffExp2`EpsSeries`ESNew[0, {1.25 + 0.5 I}]}]];
assert["cpp_machine_complex_handoff_uses_arb_interval_input",
  !FailureQ[rMachineComplex] &&
  NumericQ[rMachineComplex[["Sectors", 1, "Coeffs", 1, 1, 1]]]];

(* True resonant Jordan chain: C++ receives the exact R schedule from
   Wolfram and must retain the logarithmic member and Laurent window. *)
DiffExp2`Solve`ClearSolveCaches[];
csJordan = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0, 1}, {0, 0}}/x, "Variable" -> x|>,
  chart["cpp-jordan"]];
rJordan = catchDE2[DiffExp2`Solve`SolveChart[csJordan, req[0, 2, 4]]];
jordanLog = If[FailureQ[rJordan], None,
  SelectFirst[rJordan["Basis", "Columns"], MemberQ[tags[#], {0, 0, 1}] &, None]];
assert["cpp_true_resonant_jordan_log_chain",
  !FailureQ[rJordan] && Length[rJordan["Basis", "Columns"]] === 2 &&
  AssociationQ[jordanLog] && jordanLog["EpsWindow", "Min"] === -1];

(* Inhomogeneous true resonance: theta f=1 gives Log[t]. This enters C++
   as a materialized source tensor and proves that source validity is not
   inferred from numerical zeros. *)
csSource = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0}}, "Variable" -> x|>, chart["cpp-source"]];
source = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[{If[k === 0 && n === 0, 1, 0]},
      {k, 0, 2}, {n, 0, 3}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TWindow" -> <|"CompleteMax" -> 3|>|>;
rSource = catchDE2[DiffExp2`Solve`SolveParticular[
  csSource, source, req[0, 2, 3]]];
sourceLog = If[FailureQ[rSource], None,
  SelectFirst[rSource["Sectors"], #["p"] === 1 &, None]];
assert["cpp_inhomogeneous_resonant_source_tensor",
  !FailureQ[rSource] && AssociationQ[sourceLog] &&
  rSource["EpsWindow", "Min"] === -1 &&
  Abs[N[sourceLog["Coeffs"][[1, 1, 1]] - 1, 70]] < 10^-70];

(* Extra symbolic analytic regulators use FLINT's exact multivariate
   rational-function field. Epsilon remains the independent Laurent axis;
   rho is neither sampled nor numerically specialized by the solve. *)
regRate = (1 + rho)/(2 - rho);
csReg = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{regRate}}, "Variable" -> x|>, chart["cpp-regulator"]];
rReg = catchDE2[DiffExp2`Solve`SolveChart[csReg, req[0, 1, 4]]];
regCoeffs = If[FailureQ[rReg], {},
  rReg["Basis", "Columns"][[1, "Sectors", 1, "Coeffs", 1, All, 1]]];
assert["cpp_symbolic_analytic_regulator_exact_field_parity",
  !FailureQ[rReg] &&
  And @@ MapThread[Cancel[Together[#1 - #2]] === 0 &,
    {regCoeffs, Table[regRate^n/n!, {n, 0, 4}]}] &&
  FreeQ[regCoeffs, _?InexactNumberQ]];

(* The symbolic residual certificate specializes only its independent
   probes.  Corrupt an exact rho-dependent Taylor coefficient to prove
   those probes remain a loud correctness check rather than a bypass. *)
badReg = If[FailureQ[rReg], <||>, Module[{col, old},
  col = rReg["Basis", "Columns"][[1]];
  old = col["Sectors"][[1, "Coeffs", 1, 2, 1]];
  (* This corruption vanishes at the first two regular O(1) probes. The
     expanded multi-probe check must still reject it. *)
  ReplacePart[col, {"Sectors", 1, "Coeffs", 1, 2, 1} ->
    old + (rho + 1) (rho - 1)/10]]];
badRegResidual = If[FailureQ[rReg], rReg,
  catchDE2[DiffExp2`Solve`ODEResidualCheck[csReg, badReg]]];
assert["cpp_symbolic_analytic_regulator_residual_rejects_corruption",
  !FailureQ[rReg] && FailureQ[badRegResidual] &&
  badRegResidual["ID"] === "E7"];

(* A valid regulator pole may coincide with an early deterministic residual
   specialization. Such a candidate must be skipped rather than corrupting
   the exact solve or causing a spurious SectorSeries failure. *)
poleRate = 1/(rho - 1/43);
csRegPole = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{poleRate}}, "Variable" -> x|>, chart["cpp-regulator-pole"]];
rRegPole = catchDE2[DiffExp2`Solve`SolveChart[csRegPole, req[0, 1, 4]]];
poleCoeffs = If[FailureQ[rRegPole], {},
  rRegPole["Basis", "Columns"][[1, "Sectors", 1, "Coeffs", 1, All, 1]]];
assert["cpp_symbolic_residual_skips_regulator_poles",
  !FailureQ[rRegPole] &&
  And @@ MapThread[Cancel[Together[#1 - #2]] === 0 &,
    {poleCoeffs, Table[poleRate^n/n!, {n, 0, 4}]}]];

(* Real singular banana fixture: tags, honest windows, pseudo-hit metadata,
   and every coefficient agree with the Wolfram oracle. N=10 keeps this a
   focused test; the N=50 speed gate lives in the benchmark script. *)
fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
bananaSystem[] := catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
bananaChart[sys_, endpoint_] := Module[{plan},
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, endpoint}]];
  SelectFirst[Reverse[plan["Charts"]], TrueQ[#["Singular"]] &]];
solveBanana[backend_, endpoint_] := Module[{sys, cs, result, seconds},
  catchDE2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
    "EpsilonOrder" -> 5, "DivisionOrder" -> 4,
    "StepDivisionOrder" -> 4, "RecurrenceBackend" -> backend,
    "Variables" -> {}}]];
  sys = bananaSystem[];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, bananaChart[sys, endpoint]]];
  seconds = First@AbsoluteTiming[
    result = catchDE2[DiffExp2`Solve`SolveChart[cs, req[0, 5, 10]]]];
  {result, seconds}];
{bananaW, bananaWT} = solveBanana["Wolfram", 0];
{bananaC, bananaCT} = solveBanana["Cpp", 0];
bananaStructure = !FailureQ[bananaW] && !FailureQ[bananaC] &&
  Length[bananaW["Basis", "Columns"]] === Length[bananaC["Basis", "Columns"]] &&
  And @@ MapThread[(#1["EpsWindow"] === #2["EpsWindow"] &&
      tags[#1] === tags[#2]) &, {bananaW["Basis", "Columns"],
        bananaC["Basis", "Columns"]}];
bananaDifference = If[bananaStructure,
  Max[Abs[Flatten[MapThread[Flatten[N[
      #1["Sectors"][[All, "Coeffs"]] -
      #2["Sectors"][[All, "Coeffs"]], 70]] &,
    {bananaW["Basis", "Columns"], bananaC["Basis", "Columns"]}]]]], Infinity];
assert["cpp_banana_singular_structure_and_coefficient_parity",
  bananaStructure && bananaDifference < 10^-80];
Print["  banana N10 seconds: Wolfram=", N[bananaWT, 5],
  " Cpp=", N[bananaCT, 5], " maxdiff=", bananaDifference];

(* The upper endpoint is the production grouped-rational case Q=1+5 eps;
   the lower endpoint above stays polynomial. Exercise native grouped
   denominator lags, compact V assembly, and pseudo sectors together. *)
{bananaUpperW, bananaUpperWT} = solveBanana["Wolfram", 1];
{bananaUpperC, bananaUpperCT} = solveBanana["Cpp", 1];
bananaUpperStructure = !FailureQ[bananaUpperW] && !FailureQ[bananaUpperC] &&
  Length[bananaUpperW["Basis", "Columns"]] ===
    Length[bananaUpperC["Basis", "Columns"]] &&
  And @@ MapThread[(#1["EpsWindow"] === #2["EpsWindow"] &&
      tags[#1] === tags[#2]) &, {bananaUpperW["Basis", "Columns"],
        bananaUpperC["Basis", "Columns"]}];
bananaUpperDifference = If[bananaUpperStructure,
  Max[Abs[Flatten[MapThread[Flatten[N[
      #1["Sectors"][[All, "Coeffs"]] -
      #2["Sectors"][[All, "Coeffs"]], 70]] &,
    {bananaUpperW["Basis", "Columns"],
      bananaUpperC["Basis", "Columns"]}]]]], Infinity];
assert["cpp_banana_upper_grouped_rational_parity",
  bananaUpperStructure && bananaUpperDifference < 10^-80];
Print["  banana upper N10 seconds: Wolfram=", N[bananaUpperWT, 5],
  " Cpp=", N[bananaUpperCT, 5], " maxdiff=", bananaUpperDifference];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
