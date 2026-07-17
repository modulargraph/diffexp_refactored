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
