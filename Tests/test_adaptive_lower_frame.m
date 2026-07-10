(* Focused strict-recurrence lower-frame adaptation tests. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[e_] := Quiet[Catch[e, "DiffExp2Error"]];

x = Global`x; t = Global`t; eps = Global`eps; delta = Global`delta;
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];

chart[name_] := <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 1,
  "Radius" -> 1, "LocalRadius" -> 1, "Singular" -> True,
  "Name" -> name|>;
req[n_] := <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> 2|>,
  "TOrder" -> n|>;

solveOnce[A_, name_, n_, disable_] := Module[{cs, out},
  DiffExp2`Solve`ClearSolveCaches[];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[
    <|"Matrix" -> A, "Variable" -> x|>, chart[name]]];
  If[FailureQ[cs], Return[cs]];
  out = Block[{DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = disable},
    catchDE2[DiffExp2`Solve`SolveChart[cs, req[n]]]];
  out];

sameBasisQ[a_, b_] := !FailureQ[a] && !FailureQ[b] &&
  a["Basis", "Columns"] === b["Basis", "Columns"] &&
  a["Basis", "Specs"] === b["Basis", "Specs"];
adaptiveRows[s_] := s["Basis", "Diagnostics", "AdaptiveLowerFrames"];

(* A genuine scalar repeated pole cannot be trimmed.  It must widen to the
   old scalar terminal rectangle and remain coefficient-identical. *)
scalarA = {{1/x + 1/eps}};
scalarAdaptive = solveOnce[scalarA, "adaptive_scalar", 10, False];
scalarLegacy = solveOnce[scalarA, "adaptive_scalar", 10, True];
assert["adaptive_scalar_repeated_pole_terminal_retry",
  sameBasisQ[scalarAdaptive, scalarLegacy] &&
  First[adaptiveRows[scalarAdaptive]]["Attempts"] > 1 &&
  First[adaptiveRows[scalarAdaptive]]["FrameWidth"] ===
    First[adaptiveRows[scalarAdaptive]]["TerminalFrameWidth"] &&
  First[scalarAdaptive["Basis", "Columns"]]["EpsWindow", "Min"] === -10];

(* Rank one is not enough: an idempotent pole map repeats forever. *)
idempotent = {{1, 0}, {0, 0}};
idemA = IdentityMatrix[2]/x + idempotent/eps;
idemAdaptive = solveOnce[idemA, "adaptive_idempotent", 10, False];
idemLegacy = solveOnce[idemA, "adaptive_idempotent", 10, True];
assert["adaptive_rank_one_idempotent_does_not_saturate",
  sameBasisQ[idemAdaptive, idemLegacy] &&
  AnyTrue[adaptiveRows[idemAdaptive],
    #["FrameWidth"] === #["TerminalFrameWidth"] && #["Attempts"] > 1 &]];

(* Preparation expands every Taylor lag.  A pole-free j=1 must not hide a
   much deeper j=5 pole from the first adaptive frame.  At n=10 the j=5
   multiplier composes twice, so the first single-use frame widens once to
   the old terminal depth and remains coefficient-identical. *)
highLagA = {{1/x + x^4/eps^10}};
highLagAdaptive = solveOnce[highLagA, "adaptive_high_lag", 10, False];
highLagLegacy = solveOnce[highLagA, "adaptive_high_lag", 10, True];
assert["adaptive_high_lag_deep_pole_preparation_and_retry",
  sameBasisQ[highLagAdaptive, highLagLegacy] &&
  First[adaptiveRows[highLagAdaptive]]["Attempts"] === 2 &&
  First[adaptiveRows[highLagAdaptive]]["FrameWidth"] ===
    First[adaptiveRows[highLagAdaptive]]["TerminalFrameWidth"] &&
  First[highLagAdaptive["Basis", "Columns"]]["EpsWindow", "Min"] === -20];

(* The banana-like square-zero map creates one pole layer, then annihilates
   it.  Matrix-before-shift must keep the narrow first rectangle. *)
nilpotent = {{0, 1}, {0, 0}};
nilA = IdentityMatrix[2]/x + nilpotent/eps;
nilAdaptive = solveOnce[nilA, "adaptive_square_zero", 10, False];
nilLegacy = solveOnce[nilA, "adaptive_square_zero", 10, True];
assert["adaptive_square_zero_saturates_at_depth_one",
  sameBasisQ[nilAdaptive, nilLegacy] &&
  AllTrue[adaptiveRows[nilAdaptive], #["Attempts"] === 1 &] &&
  AllTrue[adaptiveRows[nilAdaptive],
    #["FrameWidth"] < #["TerminalFrameWidth"] &] &&
  Min[# ["EpsWindow", "Min"] & /@ nilAdaptive["Basis", "Columns"]] === -1];

(* The forced-terminal parity seam changes the work-frame diagnostics and
   therefore belongs in the homogeneous memo key.  Do not clear between the
   two calls: this is the cache regression itself. *)
cacheCS = catchDE2[DiffExp2`Solve`PrepareChart[
  <|"Matrix" -> nilA, "Variable" -> x|>, chart["adaptive_cache_mode"]]];
DiffExp2`Solve`ClearSolveCaches[];
cacheAdaptive = Block[{
    DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = False},
  catchDE2[DiffExp2`Solve`SolveChart[cacheCS, req[10]]]];
cacheTerminal = Block[{
    DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = True},
  catchDE2[DiffExp2`Solve`SolveChart[cacheCS, req[10]]]];
assert["adaptive_private_mode_is_part_of_memo_key",
  sameBasisQ[cacheAdaptive, cacheTerminal] &&
  AllTrue[adaptiveRows[cacheAdaptive], TrueQ[#["Adaptive"]] &] &&
  AllTrue[adaptiveRows[cacheTerminal], !TrueQ[#["Adaptive"]] &] &&
  AllTrue[adaptiveRows[cacheTerminal],
    #["FrameWidth"] === #["TerminalFrameWidth"] &]];

(* Analytic regulators are coefficient-field symbols, never structural zero.
   Exercise the exact lower-bound seam directly: a symbolic delta in the
   discarded slot must request widening.  A full SolveChart residual probe
   is intentionally numerical and does not accept a free symbolic parameter. *)
regSignal = Catch[
  Block[{DiffExp2`Solve`Private`$adaptiveLowerFrameProbe = True},
    DiffExp2`Solve`Private`shiftFrameBlock[
      {{delta, 0}}, -1, -1, 2, <|"Center" -> "delta-control"|>]],
  "DiffExp2AdaptiveLowerFrame"];
assert["adaptive_analytic_regulator_is_nonzero",
  FailureQ[regSignal] && regSignal["ID"] === "AdaptiveLowerFrame"];

(* A rational tail begins far above eps^0, then alternating negative shifts
   pull it back into the delivered window.  The production optimization
   deliberately retains the full proven UValid/upper halo, so parity is exact. *)
tailMap = {{0, eps^8/(1 - eps)}, {1/eps, 0}};
tailA = IdentityMatrix[2]/x + tailMap;
tailAdaptive = solveOnce[tailA, "adaptive_delayed_tail", 10, False];
tailLegacy = solveOnce[tailA, "adaptive_delayed_tail", 10, True];
assert["adaptive_delayed_rational_tail_upper_halo",
  sameBasisQ[tailAdaptive, tailLegacy]];

(* Real banana endpoint charts.  N=10 is the focused parity gate; N=20 is
   reported by Scripts/bench_adaptive_lower_frame.m / manual benchmark so
   this unit stays bounded. *)
fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
sysB = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
bananaPair[endpoint_] := Module[{plan, ch, cs, a, b},
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sysB, {11/23, endpoint}]];
  ch = SelectFirst[Reverse[plan["Charts"]], TrueQ[# ["Singular"]] &];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sysB, ch]];
  DiffExp2`Solve`ClearSolveCaches[];
  a = Block[{DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = False},
    catchDE2[DiffExp2`Solve`SolveChart[cs, req[10]]]];
  DiffExp2`Solve`ClearSolveCaches[];
  b = Block[{DiffExp2`Solve`Private`$disableAdaptiveLowerFrames = True},
    catchDE2[DiffExp2`Solve`SolveChart[cs, req[10]]]];
  {a, b}];

banana0 = bananaPair[0];
banana1 = bananaPair[1];
assert["adaptive_banana_x0_n10_exact_parity_and_narrower",
  sameBasisQ @@ banana0 && AllTrue[adaptiveRows[banana0[[1]]],
    #["FrameWidth"] < #["TerminalFrameWidth"] &]];
assert["adaptive_banana_x1_n10_exact_parity_and_narrower",
  sameBasisQ @@ banana1 && AllTrue[adaptiveRows[banana1[[1]]],
    #["FrameWidth"] < #["TerminalFrameWidth"] &]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["FAILED"]; Quit[1], Print["All tests PASSED!"]; Quit[0]];
