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

(* Large recurrence responses must use the bounded bulk parser without
   changing exact values, symbolic regulator identity, Acb accuracy, or the
   scalar decoder's diagnostics on malformed/mixed public input. *)
decodePrecision = 120;
exactEncoded = {"0", "-7/13", "123456789012345678901234567890/37",
  "(1 + Global`rho)/(2 - Global`rho)"};
exactScalarDecoded =
  DiffExp2`CppBackend`DecodeScalar[#, decodePrecision] & /@ exactEncoded;
exactBulkDecoded = DiffExp2`CppBackend`DecodeScalars[
  exactEncoded, decodePrecision];
assert["cpp_bulk_decode_exact_and_symbolic_parity",
  SameQ[exactBulkDecoded, exactScalarDecoded]];

acbEncoded = {
  {"1.2345678901234567890123456789e+3", "-2.5e-1", "-250", "zero"},
  {"0", "3.141592653589793238462643383279e+0", "zero", "-180"},
  {"-9.5e-20", "0", "-120", "zero"}};
acbScalarDecoded =
  DiffExp2`CppBackend`DecodeScalar[#, decodePrecision] & /@ acbEncoded;
acbBulkDecoded = DiffExp2`CppBackend`DecodeScalars[
  acbEncoded, decodePrecision];
assert["cpp_bulk_decode_acb_value_precision_accuracy_parity",
  SameQ[acbBulkDecoded, acbScalarDecoded]];

malformedAcbEncoded = Join[acbEncoded[[1 ;; 2]],
  {{"not-a-decimal", "0", "zero", "zero"}}];
malformedAcbDecoded = DiffExp2`CppBackend`DecodeScalars[
  malformedAcbEncoded, decodePrecision];
assert["cpp_bulk_decode_malformed_acb_falls_back_elementwise",
  SameQ[malformedAcbDecoded[[1 ;; 2]], acbScalarDecoded[[1 ;; 2]]] &&
  FailureQ[malformedAcbDecoded[[3]]]];

mixedEncoded = {"5/11", First[acbEncoded], {"malformed"}};
mixedBulkDecoded = DiffExp2`CppBackend`DecodeScalars[
  mixedEncoded, decodePrecision];
mixedScalarDecoded =
  DiffExp2`CppBackend`DecodeScalar[#, decodePrecision] & /@ mixedEncoded;
assert["cpp_bulk_decode_mixed_records_preserves_scalar_fallback",
  SameQ[mixedBulkDecoded, mixedScalarDecoded] &&
  FailureQ[Last[mixedBulkDecoded]]];

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

(* A production L2 handoff exposed an edge BigReal whose nominal final
   decimal cells are Indeterminate to RealDigits.  The bridge must shorten
   the midpoint and enlarge (never shrink) its Arb interval to cover the
   omitted decimal places. *)
uncertainTailWitness =
  0.99999999999999018886135733657392966823529940615593`207.40966701248306;
uncertainTailInitialDigits = Floor[Precision[uncertainTailWitness]];
uncertainTailRaw = Quiet[RealDigits[
  Abs[N[uncertainTailWitness, uncertainTailInitialDigits]], 10,
  uncertainTailInitialDigits]];
uncertainTailRecord = DiffExp2`CppBackend`Private`decimalRecord[
  uncertainTailWitness, 1000];
uncertainTailEncoded = DiffExp2`CppBackend`Private`arbInexactString[
  uncertainTailWitness, 1000];
assert["cpp_uncertain_decimal_tail_retries_with_covering_interval",
  MatchQ[uncertainTailRaw, {_List, _Integer}] &&
  !VectorQ[First[uncertainTailRaw], IntegerQ] &&
  AssociationQ[uncertainTailRecord] &&
  uncertainTailRecord["Digits"] < uncertainTailInitialDigits &&
  StringQ[uncertainTailEncoded] &&
  StringStartsQ[uncertainTailEncoded, "["] &&
  uncertainTailRecord["DecimalErrorExponent"] >=
    -Floor[Accuracy[uncertainTailWitness]]];

(* A high-WP banana boundary can carry an essentially zero imaginary
   cancellation remnant with less than one reliable relative digit.
   RealDigits has no honest midpoint digit to publish, so the bridge must
   send a zero-centred enclosing ball instead of rejecting the coefficient
   or silently declaring it exact zero. *)
subdigitRemnant =
  -4.16139585655806706564333246`0.3512002696020475*^-1997;
subdigitEncoded = DiffExp2`CppBackend`Private`arbInexactString[
  subdigitRemnant, 1000];
assert["cpp_subdigit_remnant_uses_zero_centered_covering_ball",
  StringQ[subdigitEncoded] &&
  subdigitEncoded === "[0 +/- 2e-1996]"];

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

(* Rank-reduced inhomogeneous parity, including the source-frame seam.
   Fuchsian reduction produces T^-1=diag(t/(1+t),1), so the physical
   theta source must be rationally multiplied before either backend sees
   its materialized recurrence tensor. *)
csRankSource = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0, (1 + x)/x^2}, {0, 0}}, "Variable" -> x|>,
  <|"ChartVar" -> t, "Center" -> 0, "Radius" -> 1/2,
    "Name" -> "cpp-rank-source"|>];
rankSource = <|"Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {1, 0}, {0, 0}],
      {k, 0, 3}, {n, 0, 5}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 3|>,
  "TWindow" -> <|"CompleteMax" -> 5|>|>;
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Wolfram"}]];
rankSourceW = catchDE2[DiffExp2`Solve`SolveParticular[
  csRankSource, rankSource, req[0, 2, 5]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Cpp"}]];
rankSourceC = catchDE2[DiffExp2`Solve`SolveParticular[
  csRankSource, rankSource, req[0, 2, 5]]];
rankSourceDiff = If[FailureQ[rankSourceW] || FailureQ[rankSourceC] ||
    tags[rankSourceW] =!= tags[rankSourceC], Infinity,
  Max[Abs[Flatten[N[rankSourceW["Sectors"][[All, "Coeffs"]] -
    rankSourceC["Sectors"][[All, "Coeffs"]], 70]]]]];
assert["cpp_rank_reduced_nonmonomial_source_wolfram_parity",
  !FailureQ[rankSourceW] && !FailureQ[rankSourceC] &&
  csRankSource["GaugeInverse"] === DiagonalMatrix[{t/(1 + t), 1}] &&
  rankSourceW["EpsWindow"] === rankSourceC["EpsWindow"] &&
  tags[rankSourceW] === tags[rankSourceC] && rankSourceDiff < 10^-80];

(* Particular solves must apply a nonidentity epsilon-rational spectral V in
   the native recurrence before returning the physical tensor.  This is the
   direct parity pin for the compact native-assembly path. *)
csSpectralSource = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0, eps^-8}, {0, 1}}/x, "Variable" -> x|>,
  chart["cpp-particular-spectral-v"]];
spectralSource = <|"Sectors" -> {<|"a" -> 1/2, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {1, 0}, {0, 0}],
      {k, 0, 10}, {n, 0, 6}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 10|>,
  "TWindow" -> <|"CompleteMax" -> 6|>|>;
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Wolfram"}]];
spectralSourceW = catchDE2[DiffExp2`Solve`SolveParticular[
  csSpectralSource, spectralSource, req[0, 2, 6]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Cpp"}]];
spectralSourceC = catchDE2[DiffExp2`Solve`SolveParticular[
  csSpectralSource, spectralSource, req[0, 2, 6]]];
spectralSourceDiff = If[FailureQ[spectralSourceW] ||
    FailureQ[spectralSourceC] || tags[spectralSourceW] =!= tags[spectralSourceC],
  Infinity,
  Max[Abs[Flatten[N[spectralSourceW["Sectors"][[All, "Coeffs"]] -
    spectralSourceC["Sectors"][[All, "Coeffs"]], 70]]]]];
assert["cpp_particular_native_nonidentity_spectral_v_parity",
  !FailureQ[spectralSourceW] && !FailureQ[spectralSourceC] &&
  csSpectralSource["V"] =!= IdentityMatrix[2] &&
  spectralSourceW["EpsWindow"] === spectralSourceC["EpsWindow"] &&
  tags[spectralSourceW] === tags[spectralSourceC] &&
  spectralSourceDiff < 10^-80];

catchDE2[DiffExp2`Config`UpdateConfiguration[{
  "Variables" -> {rho}, "RecurrenceBackend" -> "Wolfram"}]];
csSymbolicSpectralSource = DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> {{0, (1 + rho)/eps^2}, {0, 1}}/x,
    "Variable" -> x|>, chart["cpp-particular-symbolic-spectral-v"]];
symbolicSpectralSource = <|"Sectors" -> {<|
    "a" -> 1/2, "b" -> 0, "p" -> 0,
    "Coeffs" -> Table[If[k === 0 && n === 0, {0, 1}, {0, 0}],
      {k, 0, 6}, {n, 0, 4}]|>},
  "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 6|>,
  "TWindow" -> <|"CompleteMax" -> 4|>|>;
symbolicSpectralSourceW = catchDE2[DiffExp2`Solve`SolveParticular[
  csSymbolicSpectralSource, symbolicSpectralSource, req[0, 1, 4]]];
catchDE2[DiffExp2`Config`UpdateConfiguration[
  {"RecurrenceBackend" -> "Cpp"}]];
symbolicSpectralSourceC = catchDE2[DiffExp2`Solve`SolveParticular[
  csSymbolicSpectralSource, symbolicSpectralSource, req[0, 1, 4]]];
symbolicSpectralParity = !FailureQ[symbolicSpectralSourceW] &&
  !FailureQ[symbolicSpectralSourceC] &&
  symbolicSpectralSourceW["EpsWindow"] ===
    symbolicSpectralSourceC["EpsWindow"] &&
  tags[symbolicSpectralSourceW] === tags[symbolicSpectralSourceC] &&
  AllTrue[Flatten[symbolicSpectralSourceW["Sectors"][[All, "Coeffs"]] -
    symbolicSpectralSourceC["Sectors"][[All, "Coeffs"]]],
    Cancel[Together[#]] === 0 &];
assert["cpp_particular_native_symbolic_regulator_spectral_v_parity",
  symbolicSpectralParity &&
  csSymbolicSpectralSource["V"] =!= IdentityMatrix[2]];
catchDE2[DiffExp2`Config`UpdateConfiguration[{"Variables" -> {}}]];

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
