(* Focused Wolfram smoke for full-parent physical q/C SCC manifest capture. *)
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
theta = {{t/(1 - eps), 0},
  {eps t/(1 - eps), 2 t/(1 - eps)}};
cs = <|"ChartVar" -> t, "Center" -> 0, "SystemSize" -> 2,
  "ThetaOriginal" -> theta|>;
owner = ExportString[<|"schema" -> "diffexp2-native-scc-composite-v1",
    "fixture" -> "full-parent-physical-owner"|>,
  "RawJSON", "Compact" -> True];
serialization = <|"domain" -> "rational", "symbols" -> {}|>;

DiffExp2`Solve`ClearSolveCaches[];
payload = catchDE2[
  DiffExp2`Solve`Private`sccParentPhysicalODEPayload[
    cs, owner, serialization, 80]];
cLags = Lookup[payload, "c", {}];
linearEntries = If[Length[cLags] >= 2, cLags[[2]], {}];
crossEntry = SelectFirst[linearEntries,
  Lookup[#, "r", None] === 1 && Lookup[#, "c", None] === 0 &,
  <||>];

assert["scc_parent_physical_payload_captured",
  AssociationQ[payload] &&
  payload["schema"] === "diffexp2-physical-cleared-ode-v1"];
assert["scc_parent_owner_is_exact_composite_identity",
  payload["owner_signature_identity"] === owner];
assert["scc_parent_payload_is_full_dimension",
  Length[cLags] === 2 && Length[linearEntries] === 3 &&
  Sort[Lookup[linearEntries, {"r", "c"}]] ===
    {{0, 0}, {1, 0}, {1, 1}}];
assert["scc_parent_payload_retains_epsilon_rational_cross_edge",
  payload["q"] === {<|"zero" -> False, "valuation" -> 0,
      "numerator" -> {"1"}, "denominator" -> {"1"}|>} &&
  Lookup[Lookup[crossEntry, "v", <||>], "valuation", None] === 1 &&
  Lookup[Lookup[crossEntry, "v", <||>], "denominator", None] ===
    {"1", "-1"}];
assert["scc_parent_payload_is_not_a_diagonal_block_certificate",
  AnyTrue[linearEntries,
    Lookup[#, "r"] === 1 && Lookup[#, "c"] === 0 &]];

Print["\nSCC physical ODE manifest smoke: ", passed, " passed, ", failed,
  " failed."];
If[failed > 0, Exit[1], Exit[0]];
