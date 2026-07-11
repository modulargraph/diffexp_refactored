(* Unit tests for DiffExp2/Config.m per Docs/specs/Config.md section 8.
   E12 waiver: untestable until API.m LoadSystem exists (M5); the
   adopt-or-error rule is enforced through a Config-side check then. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "Tolerances.m"}]];
Get[FileNameJoin[{repoRoot, "DiffExp2", "Config.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
failsWith[f_, id_String] := FailureQ[f] && f["ID"] === id && f["Module"] === "Config";

cfg = DiffExp2`Config`CFG;
upd = DiffExp2`Config`UpdateConfiguration;
load = DiffExp2`Config`LoadConfiguration;
cc = DiffExp2`Config`CurrentConfiguration;
tol = DiffExp2`Tolerances`Tol;

(* C17 (first half) — must run before any load *)
assert["test_cfg_before_load_loud",
  failsWith[catchDE2[cfg["WorkingPrecision"]], "E11"] &&
  DiffExp2`Config`ConfiguredQ[] === False];

(* C1 *)
r = catchDE2[load[{}]];
assert["test_defaults_complete",
  !FailureQ[r] &&
  Sort[Keys[cc[]]] === Sort[Keys[DiffExp2`Config`ConfigSchema[]]] &&
  cfg["WorkingPrecision"] === 500 && cfg["ChopPrecision"] === 250 &&
  cfg["LinearSolveChopPrecision"] === 250 && cfg["DivisionOrder"] === 3 &&
  cfg["EpsilonOrder"] === 4 && cfg["EstimateError"] === "Fast" &&
  cfg["ExpansionOrder"] === 50 && cfg["RadiusOfConvergence"] === 1 &&
  cfg["SegmentationStrategy"] === "Predivision" && cfg["UsePade"] === False &&
  cfg["Verbosity"] === 1 && cfg["AccuracyGoal"] === "?" &&
  cfg["AccuracyGoalValidate"] === False &&
  cfg["StrictMatchingUncertainty"] === False &&
  cfg["AbortOnAnalyticContinuationFail"] === True];
assert["test_configured_q_after_load", DiffExp2`Config`ConfiguredQ[] === True];

(* C2 *)
assert["test_read_unknown_key_loud",
  failsWith[catchDE2[cfg["RationalizationTollerance"]], "E1"]];

(* C3 *)
assert["test_symbol_key_canonicalization",
  cfg[System`WorkingPrecision] === 500 &&
  cfg[Foo`Bar`WorkingPrecision] === 500];

(* C4 *)
e1 = catchDE2[cfg[Foo`Bar`NoSuchKey]];
assert["test_symbol_key_unknown_context_diagnostic",
  failsWith[e1, "E1"] && StringContainsQ[e1["Detail"], "Foo`Bar`"]];

(* C5 *)
e2 = catchDE2[upd[{"WorkingPrecision" -> 300, "NoSuchKey" -> 1}]];
assert["test_write_unknown_key_loud_atomic",
  failsWith[e2, "E2"] && cfg["WorkingPrecision"] === 500];

(* C6 *)
catchDE2[upd[{"WorkingPrecision" -> 300}]];
assert["test_chop_sync_automatic",
  cfg["ChopPrecision"] === 150 && cfg["LinearSolveChopPrecision"] === 150 &&
  tol["MatchTol"] === 10^-150];
catchDE2[upd[{"ChopPrecision" -> 200, "LinearSolveChopPrecision" -> 180,
  "WorkingPrecision" -> 500}]];
assert["test_chop_sync_explicit",
  tol["ChopFloor"] === 10^-200 && tol["MatchTol"] === 10^-180];
catchDE2[upd[{"ChopPrecision" -> 100}]];
assert["test_chop_sync_resets_lscp", cfg["LinearSolveChopPrecision"] === 100];

(* C7 *)
assert["test_chop_ge_wp_loud",
  failsWith[catchDE2[upd[{"WorkingPrecision" -> 100, "ChopPrecision" -> 100}]], "E4"]];

(* C8 *)
assert["test_estimate_error_strict",
  failsWith[catchDE2[upd[{"EstimateError" -> "Slow"}]], "E3"] &&
  failsWith[catchDE2[upd[{"EstimateError" -> "False"}]], "E3"]];

(* C9 *)
catchDE2[load[{"Variables" -> {Global`s, Global`m},
  "DeltaPrescriptions" -> {Global`m^2 - Global`s - I*Global`\[Delta], {Global`s - 4 Global`m^2, 1}}}]];
dp = cfg["DeltaPrescriptions"];
(* equivalence up to canonical orientation: {p,sg} == {-p,-sg} *)
sameRx[pair_, ref_] := PossibleZeroQ[pair[[1]] - ref[[1]]] && pair[[2]] === ref[[2]] ||
  PossibleZeroQ[pair[[1]] + ref[[1]]] && pair[[2]] === -ref[[2]];
assert["test_delta_prescriptions_parse",
  Length[dp] === 2 &&
  AnyTrue[dp, sameRx[#, {Global`m^2 - Global`s, -1}] &] &&
  AnyTrue[dp, sameRx[#, {Global`s - 4 Global`m^2, 1}] &]];
assert["test_delta_prescriptions_reject_nonlinear_delta",
  failsWith[catchDE2[upd[{"DeltaPrescriptions" -> {
    Global`s + I*Global`\[Delta] + Global`\[Delta]^2}}]], "E7"]];
assert["test_delta_prescriptions_reject_zero_or_nonpolynomial_factor",
  failsWith[catchDE2[upd[{"DeltaPrescriptions" -> {
      I*Global`\[Delta]}}]], "E7"] &&
  failsWith[catchDE2[upd[{"DeltaPrescriptions" -> {
      {Sin[Global`s], 1}}}]], "E7"] &&
  failsWith[catchDE2[upd[{"DeltaPrescriptions" -> {
      {Global`s + Global`\[Delta], 1}}}]], "E7"]];

(* C10 *)
assert["test_delta_prescriptions_reducible_loud",
  failsWith[catchDE2[load[{"Variables" -> {Global`s},
    "DeltaPrescriptions" -> {(Global`s - 1)^2 (Global`s - 2) + I*Global`\[Delta]}}]], "E6"]];

(* C11 *)
assert["test_line_parameter_clash_loud",
  failsWith[catchDE2[load[{"Variables" -> {Global`s}, "LineParameter" -> Global`s}]], "E5"]];

(* C12 *)
catchDE2[load[{"Variables" -> {Foo`Bar`s, Global`m}}]];
assert["test_variable_pinning",
  AllTrue[cfg["Variables"], Context[#] === "Global`" &] &&
  (SymbolName /@ cfg["Variables"]) === {"s", "m"}];

(* C13 *)
assert["test_dropped_keys_dedicated_errors",
  AllTrue[{"IntegrationStrategy", "UseRationalRecurrence", "InvWronskSolver",
    "HomogeneousSolve", "KeepMatrixExpansions", "Parallel", "IgnoreIndicialCheck"},
    Module[{f = catchDE2[upd[{# -> True}]]},
      failsWith[f, "E9"] && StringContainsQ[f["Detail"], "removed"]] &] &&
  failsWith[catchDE2[upd[{"SegmentationStrategy" -> "Dynamic"}]], "E10"]];

(* C14 *)
catchDE2[load[{}]];
assert["test_accuracy_goal_cross_field",
  failsWith[catchDE2[upd[{"AccuracyGoalValidate" -> True}]], "E8"] &&
  failsWith[catchDE2[upd[{"AccuracyGoalValidate" -> "Before"}]], "E8"] &&
  !FailureQ[catchDE2[upd[{"AccuracyGoal" -> 30, "AccuracyGoalValidate" -> True}]]] &&
  cfg["AccuracyGoalValidate"] === True];

(* C15 *)
catchDE2[upd[{"Verbosity" -> 3}]];
catchDE2[load[{"WorkingPrecision" -> 300}]];
assert["test_load_resets_then_applies",
  cfg["Verbosity"] === 1 && cfg["WorkingPrecision"] === 300];

(* C16 *)
catchDE2[load[{}]];
assert["test_rationalization_tolerance_automatic",
  cfg["RationalizationTolerance"] === 10^-250 && tol["SnapTol"] === 10^-250];
catchDE2[upd[{"RationalizationTolerance" -> 10^-30}]];
assert["test_rationalization_tolerance_explicit", tol["SnapTol"] === 10^-30];

(* C18 *)
exports = {"LoadConfiguration", "UpdateConfiguration", "CurrentConfiguration",
  "CFG", "ConfigSchema", "CanonicalKey", "PinnedVariable", "ConfiguredQ",
  "PrintInfo", "PrintWarning", "EpsSymbols", "CanonicalEps"};
assert["test_exports_visible_cross_context",
  Block[{$ContextPath = {"System`"}},
    DiffExp2`Config`CanonicalKey[System`WorkingPrecision] === "WorkingPrecision"] &&
  AllTrue[exports, MemberQ[Names["DiffExp2`Config`*"], #] &] &&
  AllTrue[exports,
    StringQ[ToExpression["DiffExp2`Config`" <> # <> "::usage"]] &]];

(* C19 — static: no raw store reads outside Config.m *)
de2files = FileNames["*.m", FileNameJoin[{repoRoot, "DiffExp2"}]];
offenders = Select[de2files,
  FileBaseName[#] =!= "Config" &&
  StringContainsQ[ReadString[#], "FEC[" | "DiffExpConfiguration"] &];
assert["test_no_raw_store_reads_in_tree", offenders === {}];

(* C20 *)
assert["test_dead_keys_dropped",
  AllTrue[{"CrosscheckLevel", "CrosscheckFlags", "LogFile", "UseMobius"},
    failsWith[catchDE2[upd[{# -> 1}]], "E9"] &] &&
  failsWith[catchDE2[upd[{"UseMobius" -> False}]], "E9"]];

(* C21 — sign-aware dedup (DEC-16) *)
catchDE2[load[{"Variables" -> {Global`s, Global`m},
  "DeltaPrescriptions" -> {{Global`m^2 - Global`s, -1}, {Global`s - Global`m^2, +1}}}]];
dp2 = cfg["DeltaPrescriptions"];
assert["test_delta_prescriptions_sign_aware_collapse", Length[dp2] === 1];
assert["test_delta_prescriptions_sign_aware_conflict",
  failsWith[catchDE2[load[{"Variables" -> {Global`s, Global`m},
    "DeltaPrescriptions" -> {{Global`m^2 - Global`s, -1}, {Global`s - Global`m^2, -1}}}]], "E7"]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
