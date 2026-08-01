(* Focused Wolfram-side smoke for the frame-independent physical q/C owner. *)
repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

t = Global`t;
eps = DiffExp2`Config`CanonicalEps[];
theta = {{t (1 + eps)/((1 - eps) (1 - 2 t))}};
cs = <|"ChartVar" -> t, "Center" -> 0, "SystemSize" -> 1,
  "ThetaOriginal" -> theta|>;

DiffExp2`Solve`ClearSolveCaches[];
data1 = catchDE2[DiffExp2`Solve`Private`physicalClearedODEData[cs]];
data2 = catchDE2[DiffExp2`Solve`Private`physicalClearedODEData[cs]];

assert["physical_ode_exact_capture_succeeds",
  AssociationQ[data1] && data1 === data2];
assert["physical_ode_clearing_removes_global_epsilon_content",
  MemberQ[{1, -1}, data1["Q"][[1, "P"]]] &&
  data1["Q"][[1, "Zero"]] === False &&
  data1["Q"][[1, "Valuation"]] === 0 &&
  data1["Q"][[1, "Q"]] === 1 &&
  data1["Q"][[2, "Zero"]] === False &&
  data1["Q"][[2, "Valuation"]] === 0 &&
  data1["Q"][[2, "P"]] === -2 data1["Q"][[1, "P"]] &&
  data1["Q"][[2, "Q"]] === 1];
assert["physical_ode_retains_rational_epsilon_C",
  TrueQ[data1["C"][[1, 1, 1, "Zero"]]] &&
  data1["C"][[2, 1, 1]] === <|"Zero" -> False,
    "Valuation" -> 0,
    "P" -> data1["Q"][[1, "P"]] (1 + eps),
    "Q" -> 1 - eps|>];
assert["physical_ode_cache_full_identity_hit",
  Length[DiffExp2`Solve`Private`$physicalClearedODECache] === 1];

physicalEntry[e_Association] := If[TrueQ[e["Zero"]], 0,
  eps^e["Valuation"]*e["P"]/e["Q"]];
physicalQPolynomial[d_Association] := Sum[
  physicalEntry[d["Q"][[j]]]*t^(j - 1), {j, Length[d["Q"]]}];
physicalCPolynomial[d_Association] := Sum[
  Map[physicalEntry, d["C"][[j]], {2}]*t^(j - 1),
  {j, Length[d["C"]]}];
samePhysicalEquationQ[a_Association, b_Association] := AllTrue[
  Flatten[Map[Cancel[Together[#]] &,
    physicalQPolynomial[a]*physicalCPolynomial[b] -
      physicalQPolynomial[b]*physicalCPolynomial[a], {2}]],
  # === 0 &];

(* A certified polynomial q/C pair must be encodable without first dividing
   by q and asking the generic rational-matrix clearer to rediscover it. *)
pairQ = 1 + t + t^2;
pairC = {{eps t, 0}, {(1 + eps)/(1 - eps) + t^2, 2 t}};
pairCs = <|"ChartVar" -> t, "Center" -> 0, "SystemSize" -> 2,
  "ThetaOriginal" -> Map[Cancel[Together[#/pairQ]] &, pairC, {2}]|>;
pairDirect = catchDE2[
  DiffExp2`Solve`Private`physicalPolynomialPairODEData[
    pairQ, pairC, pairCs]];
pairLegacy = catchDE2[
  DiffExp2`Solve`Private`physicalClearedODEData[pairCs]];
assert["physical_ode_precleared_polynomial_pair_matches_legacy_equation",
  AssociationQ[pairDirect] && AssociationQ[pairLegacy] &&
    samePhysicalEquationQ[pairDirect, pairLegacy] &&
    Lookup[First[pairDirect["Q"]], "Valuation", None] === 0];

(* A rational SCC gauge may put a small t-denominator into C while q remains
   the already-cleared parent polynomial.  Clear that denominator on the
   pair itself; the result must be polynomial and represent exactly the same
   equation as the independent whole-matrix clearer. *)
rationalPairQ = 1 - t;
rationalPairC = {{(1 + eps)/(eps (1 + t)), 0},
  {t/(1 + t)^2, (1 - t)/(1 + t)}};
rationalPairCs = <|"ChartVar" -> t, "Center" -> 0, "SystemSize" -> 2,
  "ThetaOriginal" ->
    Map[Cancel[Together[#/rationalPairQ]] &, rationalPairC, {2}]|>;
rationalPairCleared = catchDE2[
  DiffExp2`Solve`Private`physicalClearRationalPair[
    rationalPairQ, rationalPairC, rationalPairCs]];
rationalPairDirect = If[AssociationQ[rationalPairCleared], catchDE2[
    DiffExp2`Solve`Private`physicalPolynomialPairODEData[
      rationalPairCleared["QExpr"], rationalPairCleared["CMatrix"],
      rationalPairCs]], rationalPairCleared];
rationalPairLegacy = catchDE2[
  DiffExp2`Solve`Private`physicalClearedODEData[rationalPairCs]];
assert["physical_ode_gauge_local_rational_pair_matches_legacy_equation",
  AssociationQ[rationalPairCleared] &&
    PolynomialQ[rationalPairCleared["QExpr"], t] &&
    AllTrue[Flatten[rationalPairCleared["CMatrix"]],
      PolynomialQ[#, t] &] &&
    AssociationQ[rationalPairDirect] && AssociationQ[rationalPairLegacy] &&
    samePhysicalEquationQ[rationalPairDirect, rationalPairLegacy] &&
    Lookup[First[rationalPairDirect["Q"]], "Valuation", None] === 0];

centerPolePair = catchDE2[
  DiffExp2`Solve`Private`physicalClearRationalPair[
    1, {{1/t}}, <|"ChartVar" -> t, "Center" -> 0,
      "SystemSize" -> 1|>]];
assert["physical_ode_gauge_local_clear_preserves_genuine_center_pole",
  AssociationQ[centerPolePair] &&
    TrueQ[centerPolePair["GenuineCenterPole"]] &&
    centerPolePair["CenterPower"] === 0 &&
    centerPolePair["QCenterPower"] === 1];

(* Exercise the actual singular-tail dispatcher with a rational gauge.  The
   parent equation is (1-t) theta f=f and f=(1+t)^-1 g.  Its reduced C has a
   denominator 1+t, so this must select the gauge-local pair clear while
   retaining a causal q(0)=1 tail payload. *)
rationalTailCs = <|"ChartVar" -> t, "Center" -> 0, "SystemSize" -> 1,
  "ThetaOriginal" -> {{1/(1 - t)}},
  "IntegrationSequence" -> <|"Components" -> {{1}}|>|>;
rationalTailBlocks = {<|"Gauge" -> {{1/(1 + t)}},
  "GaugeInverse" -> {{1 + t}}|>};
rationalTail = catchDE2[
  DiffExp2`Solve`Private`sccRationalShadowSingularTailPayload[
    rationalTailCs, rationalTailBlocks, "rational-tail-owner",
    <|"domain" -> "rational", "symbols" -> {}|>, 80]];
assert["physical_ode_singular_tail_uses_gauge_local_rational_clear",
  AssociationQ[rationalTail] &&
    rationalTail["schema"] ===
      "diffexp2-scc-singular-tail-frame-v1" &&
    rationalTail["epsilon_shifts"] === {0} &&
    rationalTail["equation", "q"][[1, "valuation"]] === 0];

genuineRationalPoleCs = <|"ChartVar" -> t, "Center" -> 0,
  "SystemSize" -> 2, "ThetaOriginal" -> {{0, 1}, {0, 0}},
  "IntegrationSequence" -> <|"Components" -> {{1}, {2}}|>|>;
genuineRationalPoleBlocks = {
  <|"Gauge" -> {{t}}, "GaugeInverse" -> {{1/t}}|>,
  <|"Gauge" -> {{1}}, "GaugeInverse" -> {{1}}|>};
genuineRationalPoleTail = catchDE2[
  DiffExp2`Solve`Private`sccRationalShadowSingularTailPayload[
    genuineRationalPoleCs, genuineRationalPoleBlocks,
    "rational-center-pole-owner",
    <|"domain" -> "rational", "symbols" -> {}|>, 80]];
assert["physical_ode_gauge_local_clear_rejects_proved_rational_center_pole",
  genuineRationalPoleTail === None];

(* A regular chart belongs to an exact registered input system.  Its local
   physical equation is just the affine image of the global cleared pair.
   The fast path may retain a harmless coefficient-field unit that the
   independent legacy LCM/GCD construction removes, so compare q1 C2=q2 C1
   rather than a particular scalar normalization. *)
x = Global`x;
system = <|"Matrix" ->
    {{(1 + eps)/((1 - eps) (1 - 2 x)), 1/(1 - x)},
     {0, eps/(2 - x)}}, "Variable" -> x|>;
chart = <|"ChartVar" -> t, "Center" -> 1/3, "Scale" -> 2/5,
  "Radius" -> 1, "LocalRadius" -> 1, "Prescriptions" -> {}|>;
chartSystem = catchDE2[DiffExp2`Solve`PrepareChart[system, chart]];
physicalChart = If[AssociationQ[chartSystem],
  DiffExp2`Solve`Private`regularPhysicalChartSystem[chartSystem],
  chartSystem];
hoistedData = If[AssociationQ[physicalChart],
  catchDE2[DiffExp2`Solve`Private`physicalClearedODEData[physicalChart]],
  physicalChart];
legacyData = If[AssociationQ[physicalChart],
  Block[{DiffExp2`Solve`Private`$disableGlobalClearedHoist = True},
    catchDE2[
      DiffExp2`Solve`Private`physicalClearedODEData[physicalChart]]],
  physicalChart];
assert["physical_ode_global_affine_clear_equals_legacy_local_equation",
  AssociationQ[hoistedData] && AssociationQ[legacyData] &&
  samePhysicalEquationQ[hoistedData, legacyData] &&
  Length[DiffExp2`Solve`Private`$globalClearedCache] === 1 &&
  Length[DiffExp2`Solve`Private`$chartClearedCache] === 1];

(* Physical ownership is in the original-master frame, so the same global
   affine q/C pair remains valid when the recurrence chart is singular and
   has a nonidentity gauge or spectral frame. *)
singularSystem = <|"Matrix" -> {{(1 + eps)/x}}, "Variable" -> x|>;
singularChart = <|"ChartVar" -> t, "Center" -> 0, "Scale" -> 2,
  "Radius" -> 1, "LocalRadius" -> 1, "Prescriptions" -> {}|>;
singularChartSystem = catchDE2[
  DiffExp2`Solve`PrepareChart[singularSystem, singularChart]];
singularHoistedData = If[AssociationQ[singularChartSystem],
  catchDE2[
    DiffExp2`Solve`Private`physicalClearedODEData[singularChartSystem]],
  singularChartSystem];
singularLegacyData = If[AssociationQ[singularChartSystem],
  Block[{DiffExp2`Solve`Private`$disableGlobalClearedHoist = True},
    catchDE2[
      DiffExp2`Solve`Private`physicalClearedODEData[
        singularChartSystem]]],
  singularChartSystem];
assert["physical_ode_global_affine_clear_covers_singular_chart",
  AssociationQ[singularHoistedData] &&
  AssociationQ[singularLegacyData] &&
  !TrueQ[
    DiffExp2`Solve`Private`regularIdentityFrameQ[singularChartSystem]] &&
  samePhysicalEquationQ[singularHoistedData, singularLegacyData]];
assert["physical_ode_singular_affine_pair_cancels_artificial_center_factor",
  AssociationQ[singularHoistedData] &&
  samePhysicalEquationQ[singularHoistedData, singularLegacyData] &&
  Length[singularHoistedData["Q"]] === 1 &&
  Length[singularHoistedData["C"]] === 1 &&
  !TrueQ[singularHoistedData["Q"][[1, "Zero"]]]];

owner = "de2-operator-wolfram-physical-smoke";
payload = Block[{
    DiffExp2`Solve`Private`$cppSerializationDomain = "rational",
    DiffExp2`Solve`Private`$cppSerializationSymbols = {}},
  catchDE2[DiffExp2`Solve`Private`cppPhysicalODEPayload[
    data1, owner, 80, cs]]];
assert["physical_ode_payload_protocol_and_owner_tokens",
  AssociationQ[payload] &&
  payload["schema"] === "diffexp2-physical-cleared-ode-v1" &&
  payload["basis"] === "physical-original-master" &&
  payload["theta_coordinate"] === "local-t" &&
  payload["owner_signature_identity"] === owner &&
  payload["payload_identity"] === data1["Identity"]];
assert["physical_ode_payload_normalized_Q0_and_sparse_C",
  payload["q"][[1, "denominator"]] === {"1"} &&
  payload["q"][[2, "denominator"]] === {"1"} &&
  payload["c"][[1]] === {} && Length[payload["c"][[2]]] === 1 &&
  payload["c"][[2, 1, "r"]] === 0 &&
  payload["c"][[2, 1, "c"]] === 0 &&
  payload["c"][[2, 1, "v", "denominator"]] === {"1", "-1"}];

DiffExp2`Solve`ClearSolveCaches[];
assert["physical_ode_cache_owned_by_ClearSolveCaches",
  Length[DiffExp2`Solve`Private`$physicalClearedODECache] === 0];

Print["\nPhysical ODE payload smoke: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1], Exit[0]];
