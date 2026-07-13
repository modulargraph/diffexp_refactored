(* Process-free contract tests for exact family and L0-target ingestion. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repoRoot, "FeynmanTrick.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]
];

raw = <|
  "LoopMomenta" -> {Global`l},
  "ExternalMomenta" -> {Global`p},
  "Propagators" -> {
    1 - Global`l^2,
    3/2 - (Global`l - Global`p)^2,
    4/3 - (Global`l - 2 Global`p)^2
  },
  "Replacements" -> {Global`p^2 -> Global`s},
  "NumericalPoint" -> {Global`s -> -1},
  "Dimension" -> 4 - 2 FeynmanTrick`FTeps
|>;

defaultFamily = FeynmanTrick`CreateFamily[raw];
assert["raw family canonical schema",
  AssociationQ[defaultFamily] &&
  defaultFamily["Schema"] === "FeynmanTrick.FamilySpec/v1" &&
  defaultFamily["InputKind"] === "RawFamily" &&
  defaultFamily["NumPropagators"] === 3 &&
  defaultFamily["NumericalPoint"] === {Global`s -> -1} &&
  defaultFamily["Definition", "NumericalPoint"] === {Global`s -> -1} &&
  defaultFamily["Dimension"] === 4 - 2 FeynmanTrick`FTeps &&
  defaultFamily["Definition", "Dimension"] ===
    4 - 2 FeynmanTrick`FTeps &&
  defaultFamily["Topology", "Dimension"] ===
    4 - 2 FeynmanTrick`FTeps];
assert["automatic left-to-right merge sequence",
  defaultFamily["CombinationSequence"] === {{1, 2}, {1, 3}}];
assert["automatic scalar target",
  defaultFamily["OutputIntegrals"] === {{1, 1, 1}} &&
  defaultFamily["L0OutputRequests"][[1, "IndexVector"]] === {1, 1, 1}];

singleFamily = FeynmanTrick`CreateFamily[raw, {2, 0, 1}];
orderedFamily = FeynmanTrick`CreateFamily[
  raw, {{2, 0, 1}, {1, 1, 1}, {2, 0, 1}}];
reorderedFamily = FeynmanTrick`CreateFamily[
  raw, {{1, 1, 1}, {2, 0, 1}}];
assert["single vector normalized to ordered list",
  singleFamily["OutputIntegrals"] === {{2, 0, 1}}];
assert["multiple targets preserve order and multiplicity",
  orderedFamily["OutputIntegrals"] ===
    {{2, 0, 1}, {1, 1, 1}, {2, 0, 1}} &&
  orderedFamily["L0OutputRequests"][[All, "RequestOrdinal"]] ===
    {1, 2, 3}];
assert["request identity is stable and content based",
  orderedFamily["L0OutputRequests"][[1, "RequestID"]] ===
    orderedFamily["L0OutputRequests"][[3, "RequestID"]] &&
  orderedFamily["L0OutputRequests"][[1, "RequestID"]] ===
    reorderedFamily["L0OutputRequests"][[2, "RequestID"]] &&
  defaultFamily["FamilyID"] === reorderedFamily["FamilyID"]];

otherPointFamily = FeynmanTrick`CreateFamily[
  Join[raw, <|"NumericalPoint" -> {Global`s -> -2}|>]];
otherDimensionFamily = FeynmanTrick`CreateFamily[
  Join[raw, <|"Dimension" -> 2 - 2 FeynmanTrick`FTeps|>]];
resolvedDimensionFamily = FeynmanTrick`CreateFamily[
  KeyDrop[raw, "Dimension"]];
assert["semantic family identity includes numerical point and dimension",
  defaultFamily["FamilyID"] =!= otherPointFamily["FamilyID"] &&
  defaultFamily["FamilyID"] =!= otherDimensionFamily["FamilyID"]];
assert["missing Dimension resolves to explicit configured semantics",
  resolvedDimensionFamily["Dimension"] ===
    FeynmanTrick`Private`DimensionExpression[] &&
  resolvedDimensionFamily["Definition", "Dimension"] ===
    FeynmanTrick`Private`DimensionExpression[]];

symbolicMassRaw = Join[raw, <|
  "Propagators" -> {
    Global`m1sq - Global`l^2,
    Global`m2sq - (Global`l - Global`p)^2,
    Global`m3sq - (Global`l - 2 Global`p)^2
  },
  "NumericalPoint" -> {
    Global`s -> -1, Global`m1sq -> 2,
    Global`m2sq -> 3/2, Global`m3sq -> 4/3
  }
|>];
symbolicMassFamily = FeynmanTrick`CreateFamily[symbolicMassRaw];
symbolicMassLevels =
  FeynmanTrick`FeynmanTrickIteration`BuildAllLevels[
    FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
      symbolicMassFamily["Topology"],
      symbolicMassFamily["CombinationSequence"],
      symbolicMassFamily["NumericalPoint"]]];
assert["NumericalPoint freezes symbols occurring directly in propagators",
  AssociationQ[symbolicMassLevels] &&
    FreeQ[Lookup[symbolicMassLevels["Levels", 1], "Propagators", {}],
      Global`m1sq | Global`m2sq | Global`m3sq] &&
    FreeQ[Lookup[symbolicMassLevels["Levels", 2], "Propagators", {}],
      Global`m1sq | Global`m2sq | Global`m3sq]];

allFamily = FeynmanTrick`CreateFamily[raw, All];
assert["All remains an explicit deferred discovery request",
  allFamily["OutputIntegrals"] === All &&
  allFamily["OutputIntegralMode"] === "AllPendingDiscovery" &&
  allFamily["L0OutputRequests"][[1, "IndexVector"]] === All &&
  allFamily["L0OutputRequests"][[1, "Resolved"]] === False];

topology = FeynmanTrick`FIREInterface`DefineTopology[
  "existing", raw["LoopMomenta"], raw["ExternalMomenta"],
  raw["Propagators"], raw["Replacements"]];
topology["Dimension"] = raw["Dimension"];
topology["NumericalPoint"] = raw["NumericalPoint"];
fromTopology = FeynmanTrick`CreateFamily[
  topology, "CombinationSequence" -> {{2, 3}, {1, 2}},
  "OutputIntegrals" -> {{0, 1, 2}, {1, 0, 0}}];
assert["existing topology accepted with exact metadata",
  fromTopology["InputKind"] === "Topology" &&
  fromTopology["Topology", "Name"] === "existing" &&
  fromTopology["CombinationSequence"] === {{2, 3}, {1, 2}} &&
  fromTopology["OutputIntegrals"] === {{0, 1, 2}, {1, 0, 0}}];

badLength = Quiet[FeynmanTrick`CreateFamily[raw, {1, 1}],
  FeynmanTrick`FamilySpec`NormalizeOutputIntegrals::length];
badInteger = Quiet[FeynmanTrick`CreateFamily[raw, {1, 1/2, 1}],
  FeynmanTrick`FamilySpec`NormalizeOutputIntegrals::integer];
badExact = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"Replacements" -> {Global`p^2 -> -1.0}|>]],
  FeynmanTrick`FamilySpec`CreateFamily::inexact];
badDelayed = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"Replacements" -> {Global`p^2 :> -1}|>]],
  FeynmanTrick`FamilySpec`CreateFamily::rules];
badDuplicateRule = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"Replacements" -> {
    Global`p^2 -> -1, Global`p^2 -> -2}|>]],
  FeynmanTrick`FamilySpec`CreateFamily::rules];
badNumericalPoint = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"NumericalPoint" -> {Global`s -> -1.0}|>]],
  FeynmanTrick`FamilySpec`CreateFamily::inexact];
badDimension = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"Dimension" -> Indeterminate|>]],
  {FeynmanTrick`FamilySpec`CreateFamily::dimension, Power::infy}];
badName = Quiet[FeynmanTrick`CreateFamily[
  raw, "Name" -> "not-fire-safe"],
  FeynmanTrick`FamilySpec`CreateFamily::name];
badMomenta = Quiet[FeynmanTrick`CreateFamily[
  Join[raw, <|"ExternalMomenta" -> {Global`p, Global`l}|>]],
  FeynmanTrick`FamilySpec`CreateFamily::momenta];
badSequence = Quiet[FeynmanTrick`CreateFamily[
  raw, "CombinationSequence" -> {{1, 2}, {2, 3}}],
  FeynmanTrick`FamilySpec`CreateFamily::pair];
badRepeatedSequence = Quiet[FeynmanTrick`CreateFamily[
  raw, "CombinationSequence" -> {{1, 2}, {1, 2}}],
  FeynmanTrick`FamilySpec`CreateFamily::pair];
assert["malformed target arity fails loudly", badLength === $Failed];
assert["noninteger target fails loudly", badInteger === $Failed];
assert["inexact family data rejected", badExact === $Failed];
assert["RuleDelayed replacements rejected", badDelayed === $Failed];
assert["duplicate replacement left-hand sides rejected",
  badDuplicateRule === $Failed];
assert["inexact numerical point rejected", badNumericalPoint === $Failed];
assert["nonfinite dimension rejected", badDimension === $Failed];
assert["unsafe FIRE name rejected", badName === $Failed];
assert["momentum symbols must be unique and disjoint",
  badMomenta === $Failed];
assert["merge sequence cannot reuse eliminated position",
  badSequence === $Failed];
assert["same-length merge sequence cannot reuse a removed position",
  badRepeatedSequence === $Failed];

topologyOnlyPrescriptionFamily = Join[defaultFamily, <|
  "Topology" -> Join[defaultFamily["Topology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>],
  "TopTopology" -> Join[defaultFamily["TopTopology"], <|
    "Prescriptions" -> {{Global`x, 1}}|>]
|>];
badTopologyOnlyPrescription = Quiet[
  FeynmanTrick`CreateFamily[topologyOnlyPrescriptionFamily],
  FeynmanTrick`FamilySpec`CreateFamily::canonical];
prescribedFamily = FeynmanTrick`CreateFamily[
  Join[raw, <|"Prescriptions" -> {{Global`x, 1}}|>]];
replayedPrescribedFamily = FeynmanTrick`CreateFamily[prescribedFamily];
assert["canonical re-ingestion rejects topology-only prescriptions",
  badTopologyOnlyPrescription === $Failed];
assert["canonical re-ingestion preserves consistently duplicated prescriptions",
  AssociationQ[replayedPrescribedFamily] &&
    replayedPrescribedFamily["Definition", "Prescriptions"] ===
      {{Global`x, 1}} &&
    replayedPrescribedFamily["Topology", "Prescriptions"] ===
      {{Global`x, 1}}];

eliminatedTopology = topology;
eliminatedTopology["EliminatedPositions"] = {2};
eliminatedAutomatic = FeynmanTrick`CreateFamily[eliminatedTopology];
eliminatedExplicit = FeynmanTrick`CreateFamily[
  eliminatedTopology, {{2, 0, 1}, {0, 0, 3}}];
badEliminatedTarget = Quiet[FeynmanTrick`CreateFamily[
  eliminatedTopology, {1, 1, 1}],
  FeynmanTrick`FamilySpec`NormalizeOutputIntegrals::sector];
badEliminatedSequence = Quiet[FeynmanTrick`CreateFamily[
  eliminatedTopology, "CombinationSequence" -> {{1, 2}}],
  FeynmanTrick`FamilySpec`CreateFamily::pair];
assert["automatic family target zeros eliminated positions",
  eliminatedAutomatic["OutputIntegrals"] === {{1, 0, 1}}];
assert["automatic merge chain skips eliminated positions",
  eliminatedAutomatic["CombinationSequence"] === {{1, 3}} &&
  eliminatedAutomatic["FamilyID"] =!= defaultFamily["FamilyID"]];
assert["explicit zero-sector targets survive canonical forwarding",
  eliminatedExplicit["OutputIntegrals"] === {{2, 0, 1}, {0, 0, 3}}];
assert["nonzero eliminated target rejected",
  badEliminatedTarget === $Failed];
assert["explicit merge cannot target an eliminated position",
  badEliminatedSequence === $Failed];

badMergedNumerator = Quiet[FeynmanTrick`CreateFamily[
  raw, {2, 0, -1}],
  FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence::merge];
spectatorNumerator = FeynmanTrick`CreateFamily[
  raw, {2, 0, -1}, "CombinationSequence" -> {{1, 2}}];
assert["negative index at a merged position rejected",
  badMergedNumerator === $Failed];
assert["negative index at a never-merged spectator remains supported",
  spectatorNumerator["OutputIntegrals"] === {{2, 0, -1}}];

preparedTopology = topology;
preparedTopology["StartFileReady"] = True;
preparedTopology["WorkDirectory"] = "/tmp/stale";
badPrepared = Quiet[FeynmanTrick`CreateFamily[preparedTopology],
  FeynmanTrick`FamilySpec`CreateFamily::prepared];
assert["prepared FIRE topology state rejected", badPrepared === $Failed];

legacyIteration =
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}, {1, 3}}, {}];
legacyEliminatedIteration =
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    eliminatedTopology, {{1, 3}}, {}];
multiIteration =
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}, {1, 3}}, {},
    "OutputIntegrals" -> {{0, 1, 2}, {2, 0, 1}}];
eliminatedIteration =
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    eliminatedTopology, eliminatedAutomatic["CombinationSequence"], {},
    "OutputIntegrals" -> eliminatedExplicit["OutputIntegrals"]];
badIterationEliminated = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    eliminatedTopology, {{1, 3}}, {},
    "OutputIntegrals" -> {1, 1, 1}],
  FeynmanTrick`FamilySpec`NormalizeOutputIntegrals::sector];
badIterationEliminatedSequence = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    eliminatedTopology, {{1, 2}}, {},
    "OutputIntegrals" -> {1, 0, 1}],
  FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence::sequence];
badIterationNumerator = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}, {1, 3}}, {},
    "OutputIntegrals" -> {2, 0, -1}],
  FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence::merge];
badIterationLength = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}}, {}, "OutputIntegrals" -> {1, 1}],
  FeynmanTrick`FamilySpec`NormalizeOutputIntegrals::length];
badIterationAll = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}}, {}, "OutputIntegrals" -> All],
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration::alltargets];
badIterationRepeatedSequence = Quiet[
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, {{1, 2}, {1, 2}}, {},
    "OutputIntegrals" -> {1, 1, 1}],
  FeynmanTrick`FamilySpec`ValidateOutputIntegralsForSequence::sequence];
assert["DefineFTIteration Automatic preserves legacy scalar target",
  legacyIteration["Levels", 0, "Masters"] === {{1, 1, 1}} &&
  legacyEliminatedIteration["Levels", 0, "Masters"] === {{1, 0, 1}}];
assert["DefineFTIteration preserves explicit L0 target order",
  multiIteration["Levels", 0, "Masters"] ===
    {{0, 1, 2}, {2, 0, 1}}];
assert["canonical eliminated targets forward to DefineFTIteration",
  eliminatedIteration["Levels", 0, "Masters"] ===
    {{2, 0, 1}, {0, 0, 3}}];
assert["DefineFTIteration rejects nonzero eliminated target",
  badIterationEliminated === $Failed];
assert["DefineFTIteration rejects a merge through an eliminated position",
  badIterationEliminatedSequence === $Failed];
assert["DefineFTIteration rejects numerator at merged position",
  badIterationNumerator === $Failed];
assert["DefineFTIteration rejects malformed target arity",
  badIterationLength === $Failed];
assert["DefineFTIteration defers unsupported All discovery loudly",
  badIterationAll === $Failed];
assert["DefineFTIteration rejects a repeated removed merge position",
  badIterationRepeatedSequence === $Failed];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
