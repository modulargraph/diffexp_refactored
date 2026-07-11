(* Unit tests for DiffExp2/SectorSeries.m per Docs/specs/SectorSeries.md
   section 8.  Cut items per the taken M0 cuts: t17 (multi-point form cut 3),
   t23/t29 Coordinates->Main variants (cut 1: chart-coordinate checks). *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];
failsWith[f_, id_String] := FailureQ[f] && f["ID"] === id && f["Module"] === "SectorSeries";

catchDE2[DiffExp2`Config`LoadConfiguration[{"WorkingPrecision" -> 100}]];

val = DiffExp2`SectorSeries`ValidateLocalSolution;
canon = DiffExp2`SectorSeries`CanonicalizeLocalSolution;
imsign = DiffExp2`SectorSeries`ChartImSign;
eval = DiffExp2`SectorSeries`EvaluateLocalSolution;
mul = DiffExp2`SectorSeries`MultiplyRational;
reex = DiffExp2`SectorSeries`ReexpandLocalSolution;
diff = DiffExp2`SectorSeries`DifferentiateLocalSolution;
sdec = DiffExp2`SectorSeries`SectorDecomposition;
combine = DiffExp2`SectorSeries`CombineLocalSolutions;
parse = DiffExp2`SectorSeries`ParseTaggedPower;
esC = DiffExp2`EpsSeries`ESCoefficient;
esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax;

(* helper: build a LocalSolution. sectors: {a,b,p} -> rows assoc k -> {col vals} *)
mkls[sectors_List, kmin_, kmax_, ncols_, opts___] := Module[{base},
  base = <|"Center" -> 0, "ChartMap" -> {"Affine", 0, 1}, "Radius" -> 1,
    "Sectors" -> Map[<|"a" -> #[[1, 1]], "b" -> #[[1, 2]], "p" -> #[[1, 3]],
      "Coeffs" -> Table[
        Module[{rows = #[[2]]},
          Table[{Lookup[rows, k, ConstantArray[0, ncols]][[n + 1]]}, {n, 0, ncols - 1}]],
        {k, kmin, kmax}][[All, All]]|> &, sectors],
    "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
    "TWindow" -> <|"CompleteMax" -> ncols - 1|>,
    "ErrorEstimate" -> ConstantArray[0, kmax - kmin + 1],
    "Prescriptions" -> {}, opts|>;
  base];

(* t01 *)
ls1 = mkls[{{{0, 0, 0}, <|0 -> {1, 2, 3}|>}, {{-1, 1, 0}, <|0 -> {4, 5, 6}|>}}, 0, 0, 3];
assert["t01_validate_roundtrip", val[ls1] === ls1];

(* t02 *)
ls2 = mkls[{{{0.5, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1];
assert["t02_validate_inexact_tag", failsWith[catchDE2[val[ls2]], "badtag"]];

(* t03: mismatched dims *)
ls3 = ls1; ls3["Sectors"] = ReplacePart[ls1["Sectors"],
  {2, Key["Coeffs"]} -> {{{1, 1}, {2, 2}, {3, 3}}}];
assert["t03_validate_shape", failsWith[catchDE2[val[ls3]], "dims"]];

(* t04 *)
ls4 = ls1; ls4["Sectors"] = ReplacePart[ls1["Sectors"],
  {1, Key["Coeffs"], 1, 1, 1} -> SeriesData[Global`x, 0, {1}, 0, 2, 1]];
assert["t04_validate_badcoeff", failsWith[catchDE2[val[ls4]], "badcoeff"]];

(* t05: exp *)
ls5 = mkls[{{{0, 0, 0}, <|0 -> Table[1/n!, {n, 0, 20}]|>}}, 0, 0, 21];
r5 = eval[ls5, 1/10, "UsePade" -> False];
assert["t05_eval_regular_exp",
  Abs[N[esC[r5["Value"], 0][[1]] - Exp[1/10], 30]] < 10^-20];
r5NoTails = eval[ls5, 1/10, "UsePade" -> False,
  "ComputeTailEstimates" -> False];
assert["t05_eval_tail_fastpath_value_parity_and_default_unchanged",
  r5NoTails["Value"] === r5["Value"] &&
    r5NoTails["PadeFallbacks"] === r5["PadeFallbacks"] &&
    ListQ[r5["TailEstimates"]] && Max[r5["TailEstimates"]] > 0 &&
    r5NoTails["TailEstimates"] === Missing["NotComputed"]];
tailScanDefault = Catch[Block[{
    DiffExp2`SectorSeries`Private`numMagBounds =
      Function[{z, wp}, Throw["tail-scan", "tail-scan"]]},
    eval[ls5, 1/10, "UsePade" -> False]], "tail-scan"];
tailScanOff = Catch[Block[{
    DiffExp2`SectorSeries`Private`numMagBounds =
      Function[{z, wp}, Throw["tail-scan", "tail-scan"]]},
    eval[ls5, 1/10, "UsePade" -> False,
      "ComputeTailEstimates" -> False]], "tail-scan"];
assert["t05_eval_tail_fastpath_skips_magnitude_scan",
  tailScanDefault === "tail-scan" && AssociationQ[tailScanOff] &&
    tailScanOff["Value"] === r5["Value"]];

ls5uncertain = mkls[{{{0, 0, 0}, <|0 -> {0``17}, 1 -> {1}|>}},
  0, 1, 1];
r5uncertain = eval[ls5uncertain, 1/3, "UsePade" -> False];
assert["t05_inexact_zero_row_does_not_shift_honest_epsilon_min",
  r5uncertain["Value", "EpsWindow", "Min"] === 0];

(* t06: branch half-integer *)
ls6 = mkls[{{{1/2, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1];
r6p = eval[ls6, -1/4, "UsePade" -> False, "ImSign" -> 1];
r6m = eval[ls6, -1/4, "UsePade" -> False, "ImSign" -> -1];
assert["t06_eval_branch_halfinteger",
  Simplify[esC[r6p["Value"], 0][[1]] - I/2] === 0 &&
  Simplify[esC[r6m["Value"], 0][[1]] + I/2] === 0];

(* t07: log chain at negative point *)
ls7 = mkls[{{{0, 0, 1}, <|0 -> {1}|>}}, 0, 0, 1];
r7 = eval[ls7, -1/4, "UsePade" -> False, "ImSign" -> 1];
assert["t07_eval_branch_logchain",
  Simplify[esC[r7["Value"], 1][[1]] - (-Log[4] + I Pi)] === 0];

(* t08: b-eps tower *)
ls8 = mkls[{{{0, 2, 0}, <|0 -> {1}|>}}, 0, 5, 1];
r8 = eval[ls8, 1/3, "UsePade" -> False];
assert["t08_eval_beps_tower",
  esCM[r8["Value"]] === 5 &&
  AllTrue[Range[0, 5],
    Simplify[esC[r8["Value"], #][[1]] - (2 Log[1/3])^#/#!] === 0 &]];
r8p = eval[ls8, -1/3, "UsePade" -> False, "ImSign" -> 1];
r8m = eval[ls8, -1/3, "UsePade" -> False, "ImSign" -> -1];
assert["t08_eval_beps_tower_negative_plus_minus_i0",
  AllTrue[Range[0, 5], Function[q,
    Simplify[esC[r8p["Value"], q][[1]] -
      (2 (Log[1/3] + I Pi))^q/q!] === 0 &&
    Simplify[esC[r8m["Value"], q][[1]] -
      (2 (Log[1/3] - I Pi))^q/q!] === 0]]];

(* t09: resonant kmin = -p *)
ls9 = mkls[{{{0, 0, 2}, <|-2 -> {1}|>}}, -2, -2, 1];
r9 = eval[ls9, 1/2, "UsePade" -> False];
assert["t09_eval_resonant_kmin",
  esMin[r9["Value"]] === 0 &&
  Simplify[esC[r9["Value"], 0][[1]] - Log[1/2]^2/2] === 0];

(* t10 *)
assert["t10_eval_radius_error",
  failsWith[catchDE2[eval[ls5, 1, "UsePade" -> False]], "radius"]];

(* t11 *)
ls11 = mkls[{{{0, 1, 0}, <|0 -> {1}|>}}, 0, 0, 1];
assert["t11_eval_origin_error",
  failsWith[catchDE2[eval[ls11, 0, "UsePade" -> False]], "originlimit"]];

(* t12 *)
ls12a = mkls[{{{1/2, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1];
ls12b = mkls[{{{0, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1];
assert["t12_eval_branch_missing",
  failsWith[catchDE2[eval[ls12a, -1/10, "UsePade" -> False]], "branchmissing"] &&
  !FailureQ[catchDE2[eval[ls12b, -1/10, "UsePade" -> False]]]];

(* A single-valued negative integer power needs no prescription.  In
   particular, the internal None sentinel must never leak into its phase. *)
ls12c = mkls[{{{3, 0, 0}, <|0 -> {2}|>}}, 0, 0, 1];
r12c = catchDE2[eval[ls12c, -1/10, "UsePade" -> False]];
assert["t12_eval_negative_integer_without_prescription",
  !FailureQ[r12c] && Simplify[esC[r12c["Value"], 0][[1]] + 1/500] === 0 &&
    FreeQ[r12c, None]];

(* t13 *)
mkPresc[entries_] := Append[ls12b, "Prescriptions" ->
  Map[<|"Factor" -> #[[1]], "Sign" -> #[[2]], "Multiplicity" -> #[[3]],
    "LeadingCoeffSign" -> #[[4]]|> &, entries]];
assert["t13_chart_imsign",
  imsign[mkPresc[{{Global`x, -1, 1, -1}}]] === 1 &&
  imsign[mkPresc[{{Global`x, -1, 3, -1}}]] === 1 &&
  imsign[mkPresc[{{Global`x, -1, 1, -1}, {Global`y, 1, 2, 1}}]] === 1 &&
  failsWith[catchDE2[imsign[mkPresc[{{Global`x, -1, 1, -1}, {Global`y, -1, 1, 1}}]]],
    "branchconflict"]];

(* t14: Pade load-bearing pin *)
ls14 = mkls[{{{0, 0, 0}, <|0 -> Table[(-1)^n, {n, 0, 10}]|>}}, 0, 0, 11];
r14p = eval[ls14, 3/4, "UsePade" -> True];
r14d = eval[ls14, 3/4, "UsePade" -> False];
assert["t14_pade_exact_rational",
  Abs[N[esC[r14d["Value"], 0][[1]] - 4/7, 20]] >= 10^-3 &&
  Simplify[esC[r14p["Value"], 0][[1]] - 4/7] === 0];

(* t15: loud fallback via hook *)
Block[{DiffExp2`SectorSeries`$ForcePadeFail = True},
  r15 = eval[ls14, 3/4, "UsePade" -> True]];
assert["t15_pade_loud_fallback",
  Length[r15["PadeFallbacks"]] === 1 &&
  esC[r15["Value"], 0] === esC[r14d["Value"], 0]];

(* t16: indeterminate + Pade *)
ls16 = mkls[{{{0, 0, 0}, <|0 -> {Global`c1, 2 Global`c1}|>}}, 0, 0, 2];
r16 = catchDE2[eval[ls16, 1/2, "UsePade" -> True]];
r16b = eval[ls16, 1/2, "UsePade" -> False];
assert["t16_pade_indeterminate",
  failsWith[r16, "indetpade"] &&
  Simplify[esC[r16b["Value"], 0][[1]] - 2 Global`c1] === 0];

(* t18: TOrderReduction seam *)
ls18 = mkls[{{{0, 0, 0}, <|0 -> Table[1/n!, {n, 0, 17}]|>}}, 0, 0, 18];
r18a = eval[ls5, 1/10, "UsePade" -> False, "TOrderReduction" -> 3];
r18b = eval[ls18, 1/10, "UsePade" -> False];
assert["t18_torder_reduction",
  esC[r18a["Value"], 0] === esC[r18b["Value"], 0]];

(* t19: center pole + far pole fold *)
ls19 = mkls[{{{0, 0, 0}, <|0 -> ConstantArray[0, 6] + UnitVector[6, 1]|>}}, 0, 0, 6];
r19 = mul[ls19, 1/(Global`t (1 - Global`t)), Global`t];
assert["t19_mul_center_pole",
  Length[r19["Sectors"]] === 1 &&
  {r19["Sectors"][[1]]["a"], r19["Sectors"][[1]]["b"], r19["Sectors"][[1]]["p"]} === {-1, 0, 0} &&
  Flatten[r19["Sectors"][[1]]["Coeffs"][[1, All, 1]]] === ConstantArray[1, 6]];
r19alg = mul[ls19, 1/(1 - (1 + Sqrt[2]) Global`t/10), Global`t];
assert["t19_mul_algebraic_taylor_coefficients_grounded_after_classification",
  InexactNumberQ[r19alg["Sectors"][[1]]["Coeffs"][[1, 2, 1]]] &&
  Precision[r19alg["Sectors"][[1]]["Coeffs"][[1, 2, 1]]] >= 190 &&
  Abs[N[r19alg["Sectors"][[1]]["Coeffs"][[1, 2, 1]] -
      (1 + Sqrt[2])/10, 180]] < 10^-180];

(* t20: eps-Laurent shift + geometric *)
r20a = mul[ls19, 1/(2 Global`eps), Global`t];
r20b = mul[mkls[{{{0, 0, 0}, <|0 -> {1, 0, 0, 0}|>}}, 0, 3, 4],
  1/(1 + Global`eps), Global`t];
assert["t20_mul_eps_laurent_shift",
  r20a["EpsWindow"] === <|"Min" -> -1, "CompleteMax" -> -1|> &&
  Flatten[r20a["Sectors"][[1]]["Coeffs"][[1, All, 1]]] === {1/2, 0, 0, 0, 0, 0} &&
  r20b["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 3|> &&
  Table[r20b["Sectors"][[1]]["Coeffs"][[k + 1, 1, 1]], {k, 0, 3}] === {1, -1, 1, -1}];

(* Prepared rational kernels are coefficient-independent.  Two different
   LocalSolutions with the same structural chart/window shape must share one
   cache entry and remain byte-for-byte equal to an uncached application. *)
DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
cacheInputA = mkls[{{{0, 0, 0}, <|0 -> {1, 2, 3, 4},
    1 -> {2, 0, 1, 0}, 2 -> {0, 1, 0, 1}, 3 -> {3, 2, 1, 0}|>}}, 0, 3, 4];
cacheInputB = mkls[{{{0, 0, 0}, <|0 -> {4, 3, 2, 1},
    1 -> {1, 1, 0, 0}, 2 -> {2, 0, 2, 0}, 3 -> {0, 1, 1, 0}|>}}, 0, 3, 4];
cacheMultiplier = (1 + Global`eps Global`t)/(1 - Global`t + Global`eps);
cacheResultA = mul[cacheInputA, cacheMultiplier, Global`t];
cacheResultB = mul[cacheInputB, cacheMultiplier, Global`t];
cacheEntriesAfterTwo = Length[
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache];
DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
cacheResultBUncached = mul[cacheInputB, cacheMultiplier, Global`t];
assert["t20_mul_prepared_kernel_cache_exact_parity",
  cacheEntriesAfterTwo === 1 && cacheResultB === cacheResultBUncached &&
  AssociationQ[cacheResultA]];
DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache = <||>;
contextGlobalResult = mul[cacheInputA, 1/(1 - Global`t), Global`t];
contextOtherResult = mul[cacheInputA, 1/(1 - Global`t), Other`t];
contextCacheEntries = Length[
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache];
assert["t20_mul_prepared_cache_keys_full_symbol_context",
  contextCacheEntries === 2 &&
  contextGlobalResult =!= contextOtherResult];
DiffExp2`Solve`ClearSolveCaches[];
assert["t20_mul_prepared_cache_clears_with_system_caches",
  DiffExp2`SectorSeries`Private`$multiplyRationalPreparedCache === <||>];

(* t21 *)
assert["t21_mul_interior_pole",
  failsWith[catchDE2[mul[ls19, 1/(Global`t - 1/2), Global`t]], "interiorpole"]];

(* t22 *)
assert["t22_mul_nonrational",
  failsWith[catchDE2[mul[ls19, Sqrt[Global`t], Global`t]], "nonrational"]];

(* t24: reexpand 1/t *)
ls24 = mkls[{{{-1, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1];
r24 = reex[ls24, 1/4, 10];
assert["t24_reexpand_singular_source",
  Length[r24["Sectors"]] === 1 &&
  {r24["Sectors"][[1]]["a"], r24["Sectors"][[1]]["b"], r24["Sectors"][[1]]["p"]} === {0, 0, 0} &&
  r24["Radius"] === 1/4 && r24["Prescriptions"] === {} &&
  Table[Simplify[r24["Sectors"][[1]]["Coeffs"][[1, m + 1, 1]] - (-1)^m 4^(m + 1)],
    {m, 0, 10}] === ConstantArray[0, 11]];
r24negative = catchDE2[reex[ls24, -1/4, 4]];
assert["t24_reexpand_negative_integer_pole_without_prescription",
  !FailureQ[r24negative] && FreeQ[r24negative, None] &&
    Flatten[r24negative["Sectors"][[1]]["Coeffs"][[1, All, 1]]] ===
      {-4, -16, -64, -256, -1024}];

ls24scaled = Join[ls24, <|"Center" -> 2,
  "ChartMap" -> <|"Center" -> 2, "Scale" -> 1/5|>|>];
r24scaled = reex[ls24scaled, 1/4, 4];
assert["t24_reexpand_composes_affine_chart_scale",
  r24scaled["Center"] === 41/20 &&
  r24scaled["ChartMap"] === <|"Center" -> 41/20, "Scale" -> 1/5|> &&
  r24scaled["Radius"] === 1/4];

(* t25: b-collapse *)
ls25 = mkls[{{{0, 1, 0}, <|0 -> {1}|>}}, 0, 3, 1];
r25 = reex[ls25, 1/2, 4];
assert["t25_reexpand_beps_collapse",
  Length[r25["Sectors"]] === 1 &&
  r25["Sectors"][[1]]["b"] === 0 &&
  AllTrue[Range[0, 3],
    Simplify[r25["Sectors"][[1]]["Coeffs"][[# + 1, 1, 1]] - Log[1/2]^#/#!] === 0 &] &&
  Table[Simplify[r25["Sectors"][[1]]["Coeffs"][[2, m + 1, 1]]], {m, 1, 4}] ===
    {2, -2, 8/3, -4}];

(* t26 *)
assert["t26_reexpand_errors",
  failsWith[catchDE2[reex[ls24, 0, 5]], "resingular"] &&
  failsWith[catchDE2[reex[ls24, 2, 5]], "radius"]];

(* t27: tagged derivative *)
ls27 = mkls[{{{1/2, 3, 1}, <|0 -> {1}|>}}, 0, 1, 1];
r27 = diff[ls27];
sec27a = SelectFirst[r27["Sectors"], #["p"] === 1 &];
sec27b = SelectFirst[r27["Sectors"], #["p"] === 0 &];
assert["t27_diff_tagged",
  sec27a["a"] === -1/2 && sec27a["b"] === 3 &&
  sec27a["Coeffs"][[1, 1, 1]] === 1/2 && sec27a["Coeffs"][[2, 1, 1]] === 3 &&
  sec27b["a"] === -1/2 && sec27b["Coeffs"][[2, 1, 1]] === 1 &&
  sec27b["Coeffs"][[1, 1, 1]] === 0];

(* t28: resonant window derivative (input carries headroom row for the shift) *)
ls9b = mkls[{{{0, 0, 2}, <|-2 -> {1}|>}}, -2, -1, 1];
r28 = diff[ls9b];
sec28 = SelectFirst[r28["Sectors"], #["p"] === 1 &];
assert["t28_diff_resonant_window",
  sec28 =!= Missing["NotFound"] && sec28["a"] === -1 &&
  sec28["Coeffs"][[2, 1, 1]] === 1];

(* t30: SectorDecomposition banana mix pin *)
ls30 = mkls[{{{0, 0, 0}, <|0 -> {7}|>}, {{-1, 1, 0}, <|0 -> {11}|>},
  {{0, 2, 0}, <|0 -> {13}|>}}, 0, 0, 1];
r30 = sdec[ls30];
assert["t30_sector_decomposition_banana_mix",
  ({#["a"], #["b"], #["p"]} & /@ r30["Sectors"]) ===
    {{-1, 1, 0}, {0, 0, 0}, {0, 2, 0}} &&
  (r30["Sectors"][[All, "Coeffs"]] // Flatten) === {11, 7, 13} &&
  r30["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 0|>];

(* t31: canonicalize merge + F10 *)
ls31 = mkls[{{{-1, 1, 0}, <|0 -> {1, 2, 0}|>}, {{0, 1, 0}, <|0 -> {5, 7, 0}|>},
  {{0, 2, 0}, <|0 -> {1, 0, 0}|>}, {{0, 2, 0}, <|0 -> {2, 0, 0}|>},
  {{3, 5, 0}, <|0 -> {0, 0, 0}|>},
  {{2, 7, 0}, <|0 -> {10^-200, 0, 0}|>}}, 0, 0, 3];
r31 = canon[ls31];
tags31 = {#["a"], #["b"], #["p"]} & /@ r31["Sectors"];
m31 = SelectFirst[r31["Sectors"], #["b"] === 1 &];
assert["t31_canonicalize_merge",
  tags31 === {{-1, 1, 0}, {0, 2, 0}, {2, 7, 0}} &&
  Flatten[m31["Coeffs"][[1, All, 1]]] === {1, 7, 7} &&
  Flatten[SelectFirst[r31["Sectors"], #["b"] === 2 &]["Coeffs"][[1, All, 1]]] === {3, 0, 0}];

(* Do not erase a higher-a tower when its known shifted content would fall
   outside the fixed Taylor slab.  Here the zero lower tower is dropped and
   the finite coefficient remains at its original tag. *)
ls31wide = mkls[{{{-5, 1, 0}, <|0 -> {0, 0, 0}|>},
  {{0, 1, 0}, <|0 -> {23, 0, 0}|>}}, 0, 0, 3];
r31wide = canon[ls31wide];
assert["t31_canonicalize_preserves_overflowing_integer_shift",
  ({#["a"], #["b"], #["p"]} & /@ r31wide["Sectors"]) === {{0, 1, 0}} &&
  First[r31wide["Sectors"]]["Coeffs"][[1, 1, 1]] === 23];

(* t32: window request error from the evaluation result *)
r32 = eval[ls8, 1/3, "UsePade" -> False];
assert["t32_window_request_error",
  FailureQ[catchDE2[DiffExp2`EpsSeries`ESCoefficient[r32["Value"], 6]]]];

(* t34: CombineLocalSolutions *)
lsA = mkls[{{{0, 1, 0}, <|0 -> {1, 2}|>}}, 0, 2, 2];
lsB = mkls[{{{0, 1, 0}, <|0 -> {10, 0}|>}, {{1, 0, 0}, <|0 -> {3, 0}|>}}, 0, 2, 2];
wA = DiffExp2`EpsSeries`ESNew[0, {2, 0, 0}];
wB = DiffExp2`EpsSeries`ESNew[-1, {1, 0, 0}];  (* Laurent weight: shifts kmin *)
r34 = combine[{wA, wB}, {lsA, lsB}];
sec34 = SelectFirst[r34["Sectors"], #["b"] === 1 && #["a"] === 0 &];
assert["t34_combine",
  r34["EpsWindow", "Min"] === -1 &&
  sec34["Coeffs"][[1, 1, 1]] === 10 &&   (* row -1: w_B(-1)*lsB row 0 *)
  sec34["Coeffs"][[2, 1, 1]] === 2];     (* row 0: w_A(0)*lsA row 0 *)

(* Error propagation uses the stable true modulus: a centered-zero
   imaginary component cannot erase a real weight, and two resolved
   components must not be reduced to an infinity norm. *)
ls34err = Join[mkls[{{{0, 0, 0}, <|0 -> {1}|>}}, 0, 0, 1],
  <|"ErrorEstimate" -> {1}|>];
complexWeight34 =
  4.4267459561002104836`0.2683567342206439 +
    0``-0.022790862816694443*I;
r34complex = combine[{complexWeight34}, {ls34err}];
r34diagonal = combine[{(8 + 8 I)*10^-25}, {ls34err}];
assert["t34_combine_error_uses_stable_true_modulus",
  TrueQ[First[r34complex["ErrorEstimate"]] > 4] &&
  TrueQ[PossibleZeroQ[First[r34diagonal["ErrorEstimate"]] -
    Sqrt[128]*10^-25]]];

(* t36: plain scalar weights are exact constants, so they preserve the
   complete input window; an exact-zero weight is non-constraining. *)
ls36a = mkls[{{{0, 0, 0}, <|0 -> {2}|>}}, 0, 3, 1];
ls36b = mkls[{{{0, 0, 0}, <|0 -> {5}|>}}, 0, 2, 1];
r36a = combine[{1, -1}, {ls36a, ls36b}];
r36b = combine[{0, 1}, {ls36b, ls36a}];
r36tiny = combine[{DiffExp2`EpsSeries`ESNew[0, {10.^-80, 0, 0, 0}]},
  {ls36a}];
r36tinyPlain = combine[{10.^-80}, {ls36a}];
r36finiteZero = combine[{DiffExp2`EpsSeries`ESZero[1], 1},
  {ls36a, ls36b}];
ls36t2 = mkls[{{{0, 0, 0}, <|0 -> {2, 3, 4}|>}}, 0, 3, 3];
ls36t1 = mkls[{{{0, 0, 0}, <|0 -> {5, 6}|>}}, 0, 3, 2];
r36t = combine[{1, 1}, {ls36t2, ls36t1}];
assert["t36_exact_scalar_weights_preserve_window",
  r36a["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 2|> &&
  esC[eval[r36a, 1/3, "UsePade" -> False]["Value"], 0][[1]] === -3 &&
  r36b["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 3|> &&
  esC[eval[r36b, 1/3, "UsePade" -> False]["Value"], 0][[1]] === 2 &&
  esC[eval[r36tiny, 1/3, "UsePade" -> False]["Value"], 0][[1]] =!= 0 &&
  esC[eval[r36tinyPlain, 1/3, "UsePade" -> False]["Value"], 0][[1]] =!= 0 &&
  r36finiteZero["EpsWindow"] === <|"Min" -> 0, "CompleteMax" -> 1|> &&
  r36t["TWindow", "CompleteMax"] === 1 &&
  Dimensions[First[r36t["Sectors"]]["Coeffs"]][[2]] === 2];

(* t35: ParseTaggedPower *)
p35a = parse[3 Global`x^(1 - 2 Global`eps) Log[Global`x]^2, Global`x, Global`eps];
p35b = parse[Global`c0 Global`x^(Global`eps), Global`x, Global`eps];
p35c = parse[Log[2 Global`x], Global`x, Global`eps];
p35d = parse[Exp[Global`x], Global`x, Global`eps];
assert["t35_parse_tagged_power",
  p35a === <|"a" -> 1, "b" -> -2, "p" -> 2, "Coefficient" -> 3|> &&
  p35b === <|"a" -> 0, "b" -> 1, "p" -> 0, "Coefficient" -> Global`c0|> &&
  p35c === DiffExp2`SectorSeries`$FailedParse &&
  p35d === DiffExp2`SectorSeries`$FailedParse];

(* t33: export hygiene *)
exports33 = {"ValidateLocalSolution", "CanonicalizeLocalSolution", "ChartImSign",
  "EvaluateLocalSolution", "MultiplyRational", "ReexpandLocalSolution",
  "DifferentiateLocalSolution", "SectorDecomposition", "CombineLocalSolutions",
  "ParseTaggedPower"};
assert["t33_export_hygiene",
  AllTrue[exports33, MemberQ[Names["DiffExp2`SectorSeries`*"], #] &] &&
  AllTrue[exports33,
    StringQ[ToExpression["DiffExp2`SectorSeries`" <> # <> "::usage"]] &]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
