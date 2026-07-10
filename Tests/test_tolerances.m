(* Unit tests for DiffExp2/Tolerances.m per Docs/specs/Tolerances.md section 8. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "Tolerances.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

(* Helper: run expr, return the Failure thrown via DE2Error or the value. *)
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
failsWith[f_, id_String] := FailureQ[f] && f["ID"] === id && f["Module"] === "Tolerances";

With[{ctx = "DiffExp2`Tolerances`"},
  T = Symbol[ctx <> #] &;
];
cd = DiffExp2`Tolerances`ChopDigits;
tol = DiffExp2`Tolerances`Tol;
nzq = DiffExp2`Tolerances`NumericallyZeroQ;
mag = DiffExp2`Tolerances`NumericMagnitude;
magBounds = DiffExp2`Tolerances`NumericMagnitudeBounds;
install = DiffExp2`Tolerances`InstallToleranceState;

(* T1 *)
assert["test_chop_digits_default", cd[500] === 250 && cd[100] === 50];

(* T2 *)
assert["test_chop_floor_exact",
  DiffExp2`Tolerances`ChopFloor[250] === 10^-250 &&
  Head[DiffExp2`Tolerances`ChopFloor[250]] === Rational];

(* T8 — must run BEFORE any install *)
assert["test_tol_before_install_loud", failsWith[catchDE2[tol["ChopFloor"]], "E1"]];

(* T3 — install the record Config would produce at wp=500 defaults *)
rec500 = <|"WorkingPrecision" -> 500, "ChopDigits" -> 250, "MatchDigits" -> 250,
  "ChopFloor" -> 10^-250, "MatchTol" -> 10^-250, "SnapTol" -> 10^-250,
  "RankTol" -> 10^-125, "LaurentLeadTol" -> 10^-24, "ResidTol" -> 10^-50|>;
catchDE2[install[rec500]];
assert["test_match_tol_sync_default",
  tol["MatchTol"] === 10^-250 && tol["ChopFloor"] === 10^-250];

(* T4 *)
assert["test_snap_rank_derivation",
  DiffExp2`Tolerances`SnapTol[500] === 10^-250 &&
  DiffExp2`Tolerances`RankTol[500] === 10^-125 &&
  DiffExp2`Tolerances`GeomGuardTol[500] === 10^-250 &&
  DiffExp2`Tolerances`ChopReserve[500, 250] === 250 &&
  DiffExp2`Tolerances`LaurentLeadTol[250] === 10^-24 &&
  DiffExp2`Tolerances`LaurentLeadTol[40] === 10^-20 &&
  DiffExp2`Tolerances`ResidTol[500] === 10^-50];

(* T5 — failed installs must not clobber the installed state *)
assert["test_install_schema_loud_missing",
  failsWith[catchDE2[install[KeyDrop[rec500, "RankTol"]]], "E3"]];
assert["test_install_schema_loud_extra",
  failsWith[catchDE2[install[Join[rec500, <|"Foo" -> 1/2|>]]], "E3"]];
assert["test_install_schema_loud_realtype",
  failsWith[catchDE2[install[Join[rec500, <|"ChopDigits" -> 250.0|>]]], "E3"]];

(* T6 *)
assert["test_install_chop_ge_wp_loud",
  failsWith[catchDE2[install[Join[rec500, <|"ChopDigits" -> 500, "MatchDigits" -> 500|>]]], "E4"]];

(* T7 *)
assert["test_tol_unknown_key_loud", failsWith[catchDE2[tol["ChopPrecision"]], "E2"]];
assert["test_tol_state_survived_bad_installs", tol["RankTol"] === 10^-125];

(* T9 *)
assert["test_numerically_zero_band_true", nzq[10^-300, 1, 10^-125, "t"] === True];
assert["test_numerically_zero_band_false", nzq[10^-40, 1, 10^-125, "t"] === False];
e5 = catchDE2[nzq[10^-125, 1, 10^-125, "t"]];
assert["test_numerically_zero_band_aborts",
  failsWith[e5, "E5"] && e5["Context"] === "t"];
assert["test_numerically_zero_band_edges",
  nzq[10^-130, 1, 10^-125, "t"] === True && nzq[10^-120, 1, 10^-125, "t"] === False];
assert["test_numerically_zero_floor_binary",
  nzq[10^-29, 1, 10^-24, "t", 4, True] === True &&
  nzq[10^-24, 1, 10^-24, "t", 4, True] === True &&
  nzq[10^-22, 1, 10^-24, "t", 4, True] === False];
rankFloorE5 = catchDE2[nzq[10^-24, 1,
  DiffExp2`Tolerances`RankTol[96], "rank-at-floor", 4, False]];
assert["test_ranktol_equal_to_floor_remains_ternary",
  DiffExp2`Tolerances`RankTol[96] === 10^-24 &&
  failsWith[rankFloorE5, "E5"] && rankFloorE5["BinaryFloor"] === False];

(* A zero imaginary component with poor accuracy must not poison the
   magnitude of a clearly nonzero real pivot.  Raw Abs collapses this exact
   Wolfram shape to 0``negative-accuracy; the classifier remains nonzero. *)
complexPivot =
  4.4267459561002104836`0.2683567342206439 +
    0``-0.022790862816694443*I;
complexScale = 11485.2588999987283699804`1.5130093642078184;
assert["test_complex_zero_imaginary_does_not_poison_rank_magnitude",
  PossibleZeroQ[Abs[N[complexPivot, 250]]] &&
  TrueQ[First[magBounds[complexPivot, 250]] > 0] &&
  nzq[complexPivot, complexScale, 10^-125, "complex-pivot"] === False];

(* The cure must retain the TRUE modulus, not replace it by an infinity
   norm: both resolved components contribute across a binary threshold. *)
diagonalFloor = (8 + 8 I)*10^-25;
assert["test_complex_magnitude_is_stable_true_modulus",
  TrueQ[PossibleZeroQ[mag[3 + 4 I, 50] - 5]] &&
  TrueQ[mag[diagonalFloor, 100] > 10^-24] &&
  nzq[diagonalFloor, 1, 10^-24, "diagonal-floor", 4, True] === False];

underresolvedRank = catchDE2[nzq[0``125, 1, 10^-125,
  "underresolved-rank", 4, False]];
assert["test_componentwise_uncertainty_bounds",
  failsWith[underresolvedRank, "E5"] &&
  nzq[0``140, 1, 10^-125, "resolved-rank-zero", 4, False] === True &&
  nzq[0``17, 1, 10^-24, "underresolved-floor", 4, True] === False];

(* T10 *)
assert["test_numerically_zero_exact_and_symbolic",
  nzq[0, 0, 10^-125, "t"] === True &&
  nzq[2 - 2, 1, 10^-125, "t"] === True &&
  nzq[someUndefinedSymbol, 1, 10^-125, "t"] === False &&
  nzq[10^-300, 0, 10^-125, "t"] === False];

(* T11 — the FT LaurentTrim regression in miniature *)
assert["test_laurent_trim_lesson_relative",
  nzq[10^-45, 10^-44, DiffExp2`Tolerances`LaurentLeadTol[250], "lead"] === False];

(* T12 *)
assert["test_operational_constants_pinned",
  DiffExp2`Tolerances`$SafetyDigits === 2 &&
  DiffExp2`Tolerances`$MinExpansionOrder === 10 &&
  DiffExp2`Tolerances`$InputPrecisionFactor === 2 &&
  DiffExp2`Tolerances`$MaxExtraPrecisionValue === 1000 &&
  DiffExp2`Tolerances`$AmbiguityBandDecades === 4 &&
  DiffExp2`Tolerances`$NearSingularityGuardDecades === 6 &&
  DiffExp2`Tolerances`EvalErrorSeriesDecrease[1] === 3 &&
  DiffExp2`Tolerances`EvalErrorSeriesDecrease[5] === 6];
assert["test_adaptive_search_symbols_absent",
  !MemberQ[Names["DiffExp2`Tolerances`*"], "$SafetyExpansionSubtract"] &&
  !MemberQ[Names["DiffExp2`Tolerances`*"], "$ExpansionOrdersAveraging"]];

(* T13 — export visibility from a foreign context *)
exportNames = {"ChopDigits", "ChopFloor", "MatchTol", "SnapTol", "InputSnapTol",
  "RankTol", "GeomGuardTol", "ChopReserve", "LaurentLeadTol", "ResidTol",
  "NumericMagnitude", "NumericMagnitudeBounds", "NumericallyZeroQ", "DE2Error",
  "InstallToleranceState", "Tol",
  "ToleranceStateInstalledQ", "EvalErrorSeriesDecrease", "$SafetyDigits",
  "$InputPrecisionFactor", "$MaxExtraPrecisionValue", "$MinExpansionOrder",
  "$AmbiguityBandDecades", "$NearSingularityGuardDecades"};
assert["test_exports_visible_cross_context",
  Block[{$ContextPath = {"System`"}},
    DiffExp2`Tolerances`ChopDigits[500] === 250] &&
  AllTrue[exportNames, MemberQ[Names["DiffExp2`Tolerances`*"], #] &] &&
  AllTrue[exportNames,
    StringQ[ToExpression["DiffExp2`Tolerances`" <> # <> "::usage"]] &]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
