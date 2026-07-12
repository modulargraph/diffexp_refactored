(* Focused Wolfram bridge test for combine-before-endpoint semantics.  Two
   retained masters are individually t^-1 divergent, while the exact row
   {1,-1} vanishes before the native endpoint gate. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  If[Environment["DE2_REQUIRE_CPP"] === "1", Exit[1], Exit[0]]];

nativeOKQ[value_] := AssociationQ[value] &&
  Lookup[value, "status", "error"] === "ok";
request[value_Association] := DiffExp2`CppBackend`RunRequest[value];

prescription = <|"factor_exact" -> "x", "sign" -> -1,
  "multiplicity" -> 1, "leading_coefficient_sign" -> 1|>;

prepareChart[session_, identity_, center_, prescriptions_] := request[<|
  "schema" -> 2, "op" -> "chart.prepare", "session" -> session,
  "key" -> identity, "identity" -> identity,
  "analytic" -> <|
    "geometry" -> <|"center_exact" -> center,
      "scale_exact" -> "1", "radius_exact" -> "2",
      "infinite_radius" -> False,
      "prescriptions" -> prescriptions|>,
    "principal_matrix" -> {
      {<|"exact" -> "0", "proven_zero" -> True|>,
       <|"exact" -> "0", "proven_zero" -> True|>},
      {<|"exact" -> "0", "proven_zero" -> True|>,
       <|"exact" -> "0", "proven_zero" -> True|>}},
    "native_scc_capabilities" -> <|"regular" -> False,
      "identity_gauge" -> True, "identity_v" -> True,
      "no_pseudo" -> True|>|>,
  "scc" -> <|"components" -> {{0}, {1}},
    "structural_edges" -> {}, "condensation_edges" -> {},
    "topological_order" -> {0, 1}, "coupling_depth" -> 0|>,
  "problem" -> <|"domain" -> "rational", "d" -> 2,
    "fb" -> 0, "w" -> 3,
    "d_lags" -> {{<|"s" -> 0, "v" -> "1"|>}},
    "denominators" -> {},
    "nhat_lags" -> {<|"poly" -> {}, "rat" -> {},
      "val" -> {Null, Null, Null, Null}|>},
    "d0_inverse" -> "1", "blocks" -> {{0}, {1}},
    "assembly" -> <|"identity" -> True, "poly" -> {},
      "rat" -> {}, "val" -> {0, Null, Null, 0}|>,
    "chop_digits" -> 0|>|>];

solveDivergentLocal[session_, chart_] := request[<|
  "schema" -> 2, "op" -> "local.solve", "session" -> session,
  "chart" -> chart,
  "run" -> <|"nmax" -> 2, "p" -> 0, "has_initial" -> True,
    "adaptive_probe" -> False, "a_target" -> "-1",
    "b_target" -> "0", "a_shift_min" -> 0,
    "a_shifts" -> {"-1", "0", "1"},
    "schedule" -> {
      {<|"case" -> "R", "da" -> "0", "db" -> "0"|>,
       <|"case" -> "R", "da" -> "0", "db" -> "0"|>},
      {<|"case" -> "R", "da" -> "0", "db" -> "0"|>,
       <|"case" -> "R", "da" -> "0", "db" -> "0"|>},
      {<|"case" -> "T", "da" -> "1", "db" -> "0"|>,
       <|"case" -> "T", "da" -> "1", "db" -> "0"|>}},
    "initial" -> {"1", "0", "0", "1", "0", "0"},
    "initial_validity" -> {2, 2}, "source" -> Null,
    "return_u" -> False|>,
  "metadata" -> <|
    "chart" -> <|"center_exact" -> "0", "scale_exact" -> "1",
      "radius" -> "2", "infinite_radius" -> False|>,
    "tag" -> <|
      "a" -> <|"domain" -> "rational", "canonical" -> "-1"|>,
      "b" -> <|"domain" -> "rational", "canonical" -> "0"|>|>,
    "prescriptions" -> {prescription},
    "checkpoint_identity" -> "weighted-endpoint-source"|>|>];

constantMultiplier[value_, identity_] := <|
  "epsilon_shift" -> 0, "center_pole_order" -> 0,
  "kernels" -> {{value, "0", "0"}, {"0", "0", "0"},
    {"0", "0", "0"}},
  "exact_identity" -> identity, "proven_zero" -> False|>;

preparedRow[left_, right_, identity_] := <|
  "schema" -> "diffexp2-prepared-rational-local-row-v1",
  "columns" -> 2, "exact_identity" -> identity,
  "entries" -> DeleteCases[{
    If[left === "0", Nothing,
      <|"column" -> 0,
        "multiplier" -> constantMultiplier[left, identity <> ":0"]|>],
    If[right === "0", Nothing,
      <|"column" -> 1,
        "multiplier" -> constantMultiplier[right, identity <> ":1"]|>]},
    Nothing]|>;

checkpointPath = FileNameJoin[{$TemporaryDirectory,
  "diffexp2-weighted-planned-endpoint-" <> ToString[$ProcessID] <> ".de2cp"}];
Quiet[DeleteFile[checkpointPath]];

created = request[<|"schema" -> 2, "op" -> "session.create",
  "domain" -> "rational", "precision_bits" -> 256,
  "output_digits" -> 40, "chart_capacity" -> 4,
  "local_capacity" -> 8, "endpoint_capacity" -> 4,
  "tile_plan_capacity" -> 2|>];
session = Lookup[created, "session", None];
anchor = If[StringQ[session],
  prepareChart[session, "weighted-endpoint-anchor", "-1/2", {}], <||>];
finalChart = If[StringQ[session],
  prepareChart[session, "weighted-endpoint-final", "0", {prescription}],
  <||>];
source = If[nativeOKQ[finalChart],
  solveDivergentLocal[session, finalChart["chart"]], <||>];

arm = If[nativeOKQ[anchor] && nativeOKQ[finalChart], <|
  "from_exact" -> "-1/2", "to_exact" -> "0",
  "charts" -> {anchor["chart"], finalChart["chart"]},
  "topology" -> <|"singular_points" -> {"0"},
    "boundary_points" -> {}, "complex_projections" -> {},
    "branch_sheets" -> {<|"factor_exact" -> "x", "sign" -> -1|>}|>|>,
  <||>];
plan = If[nativeOKQ[source] && AssociationQ[arm] && arm =!= <||>,
  DiffExp2`CppBackend`CreatePersistentArmTilePlan[
    source, arm, "weighted-endpoint-plan", 3], <||>];

cancellation = <|"mode" -> "exact-coefficient-field"|>;
before = If[StringQ[session], request[<|"schema" -> 2,
  "op" -> "session.stats", "session" -> session|>], <||>];
firstOnly = If[nativeOKQ[plan] && nativeOKQ[source],
  DiffExp2`CppBackend`RunPersistentWeightedPlannedEndpointLimit[
    plan, "upper", source, preparedRow["1", "0", "first-only"],
    "weighted-first-only", cancellation], <||>];
afterFirst = If[StringQ[session], request[<|"schema" -> 2,
  "op" -> "session.stats", "session" -> session|>], <||>];
secondOnly = If[nativeOKQ[plan] && nativeOKQ[source],
  DiffExp2`CppBackend`RunPersistentWeightedPlannedEndpointLimit[
    plan, "upper", source, preparedRow["0", "1", "second-only"],
    "weighted-second-only", cancellation], <||>];
afterSecond = If[StringQ[session], request[<|"schema" -> 2,
  "op" -> "session.stats", "session" -> session|>], <||>];
weighted = If[nativeOKQ[plan] && nativeOKQ[source],
  DiffExp2`CppBackend`RunPersistentWeightedPlannedEndpointLimit[
    plan, "upper", source, preparedRow["1", "-1", "difference"],
    "weighted-cancel", cancellation], <||>];
afterWeighted = If[StringQ[session], request[<|"schema" -> 2,
  "op" -> "session.stats", "session" -> session|>], <||>];

releasedSource = If[nativeOKQ[source],
  DiffExp2`CppBackend`ReleasePersistentLocal[source], <||>];
releasedPlan = If[nativeOKQ[plan],
  DiffExp2`CppBackend`ReleasePersistentTilePlan[plan], <||>];
If[nativeOKQ[anchor], request[<|"schema" -> 2, "op" -> "chart.release",
  "session" -> session, "chart" -> anchor["chart"]|>]];
If[nativeOKQ[finalChart], request[<|"schema" -> 2,
  "op" -> "chart.release", "session" -> session,
  "chart" -> finalChart["chart"]|>]];

ownedStats = If[nativeOKQ[weighted],
  DiffExp2`CppBackend`PersistentEndpointStatistics[weighted], <||>];
exported = If[nativeOKQ[weighted],
  DiffExp2`CppBackend`ExportPersistentEndpoint[weighted,
    "weighted-cancel:weighted-endpoint", 40], <||>];
saved = If[nativeOKQ[weighted],
  DiffExp2`CppBackend`SavePersistentCheckpoint[
    weighted, checkpointPath, "weighted-endpoint-roundtrip"], <||>];

If[StringQ[session], request[<|"schema" -> 2, "op" -> "session.close",
  "session" -> session|>]];
restored = If[nativeOKQ[saved],
  DiffExp2`CppBackend`RestorePersistentCheckpoint[
    checkpointPath, "weighted-endpoint-roundtrip"], <||>];
restoredSession = Lookup[restored, "session", None];
restoredHandle = If[nativeOKQ[weighted], weighted["endpoint"], None];
restoredEndpoint = If[StringQ[restoredSession] && StringQ[restoredHandle],
  <|"session" -> restoredSession, "endpoint" -> restoredHandle|>, <||>];
restoredStats = If[AssociationQ[restoredEndpoint] && restoredEndpoint =!= <||>,
  DiffExp2`CppBackend`PersistentEndpointStatistics[restoredEndpoint], <||>];
restoredExport = If[AssociationQ[restoredEndpoint] && restoredEndpoint =!= <||>,
  DiffExp2`CppBackend`ExportPersistentEndpoint[restoredEndpoint,
    "weighted-cancel:weighted-endpoint", 40], <||>];

composition = Lookup[weighted, "weighted_composition", <||>];
zeroCoefficientQ[response_] := Module[{coefficients},
  coefficients = Lookup[Lookup[response, "value", <||>],
    "coefficients", {}];
  ListQ[coefficients] && coefficients =!= {} &&
    AllTrue[Flatten[coefficients], # === "0" || # === "zero" &]];

ok = nativeOKQ[created] && nativeOKQ[anchor] && nativeOKQ[finalChart] &&
  nativeOKQ[source] && nativeOKQ[plan] &&
  AssociationQ[firstOnly] && Lookup[firstOnly, "status", "ok"] === "error" &&
  Lookup[firstOnly, "id", None] === "E2" &&
  AssociationQ[secondOnly] && Lookup[secondOnly, "status", "ok"] === "error" &&
  Lookup[secondOnly, "id", None] === "E2" &&
  Lookup[before, "locals", None] === 1 &&
  Lookup[afterFirst, "locals", None] === 1 &&
  Lookup[afterSecond, "locals", None] === 1 &&
  Lookup[afterFirst, "endpoints", None] === 0 &&
  Lookup[afterSecond, "endpoints", None] === 0 &&
  nativeOKQ[weighted] && Lookup[weighted, "native_retained", False] &&
  Lookup[weighted, "json_coefficients", None] === 0 &&
  Lookup[weighted, "dimension", None] === 1 &&
  Lookup[weighted, "effective_rim", None] === -1 &&
  Lookup[weighted, "approach_direction", None] === -1 &&
  Lookup[afterWeighted, "locals", None] === 1 &&
  Lookup[afterWeighted, "endpoints", None] === 1 &&
  Lookup[composition, "capability", None] ===
    "retained-native-weighted-plan-bound-endpoint-v1" &&
  Lookup[composition, "order", None] ===
    "rational-row-before-endpoint-gate" &&
  Lookup[composition, "row_checkpoint_identity", None] ===
    "weighted-cancel:weighted-row" &&
  Lookup[composition, "endpoint_checkpoint_identity", None] ===
    "weighted-cancel:weighted-endpoint" &&
  TrueQ[Lookup[composition, "projected_local_public_token_released", False]] &&
  nativeOKQ[releasedSource] && nativeOKQ[releasedPlan] &&
  nativeOKQ[ownedStats] && Lookup[ownedStats, "effective_rim", None] === -1 &&
  nativeOKQ[exported] && zeroCoefficientQ[exported] &&
  nativeOKQ[saved] && nativeOKQ[restored] &&
  Lookup[restored, "locals", Missing["Absent"]] === {} &&
  Lookup[restored, "tile_plans", Missing["Absent"]] === {} &&
  nativeOKQ[restoredStats] && Lookup[restoredStats, "effective_rim", None] === -1 &&
  nativeOKQ[restoredExport] && zeroCoefficientQ[restoredExport];

If[!ok,
  Print["fixture: ", <|"created" -> created, "anchor" -> anchor,
    "final" -> finalChart, "source" -> source, "plan" -> plan|>];
  Print["individual: ", <|"first" -> firstOnly,
    "second" -> secondOnly|>];
  Print["weighted: ", weighted];
  Print["stats: ", <|"before" -> before, "after_first" -> afterFirst,
    "after_second" -> afterSecond, "after_weighted" -> afterWeighted|>];
  Print["owned/export: ", <|"stats" -> ownedStats,
    "export" -> exported, "save" -> saved|>];
  Print["restore: ", <|"response" -> restored,
    "stats" -> restoredStats, "export" -> restoredExport|>]];

If[AssociationQ[restoredEndpoint] && restoredEndpoint =!= <||>,
  Quiet[DiffExp2`CppBackend`ReleasePersistentEndpoint[restoredEndpoint]]];
If[StringQ[restoredSession], request[<|"schema" -> 2,
  "op" -> "session.close", "session" -> restoredSession|>]];
Quiet[DeleteFile[checkpointPath]];

Print[If[ok, "PASS", "FAIL"],
  ": weighted retained row cancels before planned endpoint gate"];
Exit[If[ok, 0, 1]];
