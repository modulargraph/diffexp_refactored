(* End-to-end mixed regular/singular persistent native atlas. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp",
  "WorkingPrecision" -> 80,
  "ChopPrecision" -> 40,
  "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0,
  "DivisionOrder" -> 3,
  "Variables" -> {},
  "Verbosity" -> 0}];

x = Global`x; eps = Global`eps;
lambda = 1/2 + eps/3;
system = DiffExp2`LoadSystem[<|
  (* Two exact SCCs at every chart: a singular scalar affine family and a
     regular scalar family with a pole-free source edge. *)
  "Matrix" -> {{lambda/x, 1}, {0, 0}},
  "Variable" -> x|>];
lowerPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 0}];
upperPlan = DiffExp2`Transport`SegmentLine[system, {1/2, 1}];

zero = DiffExp2`EpsSeries`ESNew[-3, {0, 0, 0, 0}];
one = DiffExp2`EpsSeries`ESNew[-3, {0, 0, 0, 1}];
atlas = catchDE2[
  DiffExp2`NativeTransport`PrepareNativeRegularIndependentArms[
    system, {one, zero}, lowerPlan, upperPlan, "Threads" -> 2,
    "Integrand" -> {{1, 0}, x}]];

receivingBases = If[FailureQ[atlas], {}, Join[
  Rest[atlas["Lower", "Bases"]], Rest[atlas["Upper", "Bases"]]]];
sccBases = Select[receivingBases,
  Lookup[#, "Type", None] === "DiffExp2NativeSCCBasis" &];
singularBasis = If[FailureQ[atlas], atlas,
  SelectFirst[Rest[atlas["Lower", "Bases"]],
    Lookup[#, "Type", None] === "DiffExp2NativeSCCBasis" &,
    Missing["NotFound"]]];
regularSCCBasis = If[FailureQ[atlas], atlas,
  SelectFirst[Rest[atlas["Upper", "Bases"]],
    Lookup[#, "Type", None] === "DiffExp2NativeSCCBasis" &,
    Missing["NotFound"]]];
planStats = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`PersistentTilePlanStatistics[atlas["Plan"]]];
before = If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];

columns = If[AssociationQ[singularBasis],
  Lookup[singularBasis, "Columns", {}], {}];
forbidden = {"assembled", "coefficients", "u", "validity", "Sectors"};
opaqueBasisQ = AssociationQ[singularBasis] &&
  Lookup[singularBasis, "Dimension", None] === 2 &&
  Lookup[columns, "BasisIndex", {}] === {1, 2} &&
  AllTrue[columns, AssociationQ[#] &&
      Lookup[#, "Session", None] === atlas["Session"] &&
      Lookup[#, "NativeSCC", None] === singularBasis["NativeSCC"] &&
      StringQ[Lookup[#, "Local", None]] &&
      Intersection[Keys[#], forbidden] === {} &];

lowerCharts = If[AssociationQ[planStats],
  Lookup[Lookup[planStats, "lower", <||>], "charts", {}], {}];
upperCharts = If[AssociationQ[planStats],
  Lookup[Lookup[planStats, "upper", <||>], "charts", {}], {}];
planOwnsSCCQ = Length[lowerCharts] >= 2 &&
  Lookup[Last[lowerCharts], "chart", None] ===
    Lookup[singularBasis, "NativeSCC", None] &&
  TrueQ[Lookup[Last[lowerCharts], "singular_center", False]] &&
  Length[upperCharts] >= 2 &&
  Lookup[Last[upperCharts], "chart", None] ===
    Lookup[regularSCCBasis, "NativeSCC", None] &&
  !TrueQ[Lookup[Last[upperCharts], "singular_center", True]];

(* Remove the public SCC registry token.  The immutable tile plan and the
   already-retained opaque columns must keep typed strong ownership. *)
sccReleases = If[sccBases === {}, {Missing["NoSCCBasis"]},
  DiffExp2`CppBackend`ReleasePersistentSCC[<|
      "Session" -> #["Session"], "SCC" -> #["NativeSCC"]|>] & /@
    DeleteDuplicatesBy[sccBases, Lookup[#, "NativeSCC", None] &]];
planAfterSCCRelease = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`PersistentTilePlanStatistics[atlas["Plan"]]];
columnAfterSCCRelease = If[columns === {}, Missing["NoColumn"],
  DiffExp2`CppBackend`PersistentLocalStatistics[First[columns]]];

strongOwnershipQ = Length[sccReleases] === 2 &&
  AllTrue[sccReleases, AssociationQ[#] &&
      Lookup[#, "status", "error"] === "ok" &] &&
  AssociationQ[planAfterSCCRelease] &&
  Lookup[planAfterSCCRelease, "status", "error"] === "ok" &&
  AssociationQ[columnAfterSCCRelease] &&
  Lookup[columnAfterSCCRelease, "status", "error"] === "ok";

lowerReceivingBasis = If[FailureQ[atlas], {},
  Lookup[Rest[atlas["Lower", "Bases"]], "Columns", {}]];
transportEpsilon = If[FailureQ[atlas], <||>, <|
  "min" -> atlas["Request", "EpsWindow", "Min"],
  "max" -> atlas["Request", "EpsWindow", "CompleteMax"],
  "required_complete_max" -> atlas["TargetCompleteMax"],
  "match_required_complete_max" -> atlas["TargetCompleteMax"]|>];
transport = If[FailureQ[atlas], atlas,
  DiffExp2`CppBackend`RunPersistentTransportArm[
    atlas["Plan"], "lower", atlas["Anchor"], lowerReceivingBasis,
    transportEpsilon, "mixed-atlas-retained-lower",
    <|"relative_tolerance" -> "1e-25", "max_steps" -> 2|>]];
transportFinalRelease = If[AssociationQ[transport] &&
    Lookup[transport, "status", "error"] === "ok",
  DiffExp2`CppBackend`ReleasePersistentLocal[transport["final_local"]],
  transport];
transportStats = If[AssociationQ[transport] &&
    Lookup[transport, "status", "error"] === "ok",
  DiffExp2`CppBackend`PersistentTransportArmStatistics[transport],
  transport];
transportRelease = If[AssociationQ[transport] &&
    Lookup[transport, "status", "error"] === "ok",
  DiffExp2`CppBackend`ReleasePersistentTransportArm[transport], transport];
transportAfterRelease = If[AssociationQ[transport] &&
    Lookup[transport, "status", "error"] === "ok",
  DiffExp2`CppBackend`PersistentTransportArmStatistics[transport], transport];
afterTransport = If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];
transportQ = AssociationQ[transport] &&
  Lookup[transport, "status", "error"] === "ok" &&
  Lookup[transport, "native_retained", False] === True &&
  Lookup[transport, "json_coefficients", None] === 0 &&
  Lookup[transport, "matches", None] === 1 &&
  Lookup[transport, "tiles", None] === 2 &&
  AssociationQ[transportFinalRelease] &&
  Lookup[transportFinalRelease, "status", "error"] === "ok" &&
  AssociationQ[transportStats] &&
  Lookup[transportStats, "status", "error"] === "ok" &&
  Lookup[transportStats, "final_local", <||>]["local"] ===
    transport["final_local", "local"] &&
  AssociationQ[transportRelease] &&
  Lookup[transportRelease, "status", "error"] === "ok" &&
  AssociationQ[transportAfterRelease] &&
  Lookup[transportAfterRelease, "status", "ok"] === "error" &&
  Lookup[afterTransport, "locals", -1] === 5 &&
  Lookup[afterTransport, "transport_states", -1] === 0;

run = If[FailureQ[atlas], atlas, catchDE2[
  DiffExp2`NativeTransport`RunNativeRegularIndependentArms[
    atlas, {1, 0}, x]]];
exported = If[!AssociationQ[run] ||
    Lookup[run, "Type", None] =!=
      "DiffExp2NativeRegularIndependentArmRun", run,
  DiffExp2`CppBackend`ExportPersistentLineIntegral[
    run["CombinedLine"],
    Lookup[run["CombinedLine"], "checkpoint_identity", ""], 60]];
decoded = If[!AssociationQ[exported] ||
    Lookup[exported, "status", "error"] =!= "ok", exported,
  DiffExp2`CppBackend`DecodeScalars[
    exported["value", "coefficients"], 60]];
epsilonMin = If[AssociationQ[exported] &&
    Lookup[exported, "status", "error"] === "ok",
  exported["value", "min"], Missing["NoMinimum"]];
epsilonZero = If[ListQ[decoded] && IntegerQ[epsilonMin] &&
    1 <= 1 - epsilonMin <= Length[decoded],
  decoded[[1 - epsilonMin]], Missing["NoEpsilonZero"]];
expected = 2 Sqrt[2]/3;
runQ = AssociationQ[run] &&
  Lookup[run, "Type", None] ===
    "DiffExp2NativeRegularIndependentArmRun" &&
  run["Lower", "Matches"] === 1 && run["Upper", "Matches"] === 1 &&
  AssociationQ[exported] && Lookup[exported, "status", "error"] === "ok" &&
  ListQ[decoded] && NumberQ[epsilonZero] &&
  TrueQ[Abs[N[epsilonZero - expected, 30]] < 10^-5];

afterRun = If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];
nativeExecutionQ = runQ &&
  Lookup[afterRun, "local_matches", -1] === 3 &&
  Lookup[afterRun, "transport_states", -1] === 0 &&
  Lookup[afterRun, "line_integrations", 0] > 0 &&
  Lookup[afterRun, "line_exports", -1] === 1;

released = If[FailureQ[atlas], atlas,
  DiffExp2`NativeTransport`ReleaseNativeRegularIndependentArms[
    If[AssociationQ[run] && Lookup[run, "Type", None] ===
        "DiffExp2NativeRegularIndependentArmRun", run, atlas]]];
after = If[FailureQ[atlas], <||>,
  Lookup[DiffExp2`CppBackend`PersistentSessionInformation[],
    atlas["Session"], <||>]];

preparedQ = AssociationQ[atlas] &&
  Lookup[atlas, "Type", None] ===
    "DiffExp2NativeRegularIndependentArmAtlas" &&
  TrueQ[Lookup[atlas, "ContainsSingularReceivingCharts", False]] &&
  atlas["Lower", "BasisKinds"] === {"Anchor", "SingularSCC"} &&
  atlas["Upper", "BasisKinds"] === {"Anchor", "Regular"} &&
  Length[sccBases] === 2 &&
  AllTrue[sccBases, Lookup[#, "Session", None] === atlas["Session"] &] &&
  Lookup[before, "tile_plans", 0] === 1 &&
  Lookup[before, "scc_charts", 0] === 2 &&
  Lookup[before, "locals", 0] === 5 &&
  Lookup[before, "local_matches", -1] === 0 &&
  Lookup[before, "line_integrations", -1] === 0;

cleanupQ = AssociationQ[released] &&
  Lookup[released, "Failures", {"missing"}] === {} &&
  Lookup[after, "locals", -1] === 0 &&
  Lookup[after, "matches", -1] === 0 &&
  Lookup[after, "tile_plans", -1] === 0 &&
  Lookup[after, "line_results", -1] === 0 &&
  Lookup[after, "scc_charts", -1] === 0 &&
  Lookup[after, "local_matches", -1] === 3 &&
  Lookup[after, "line_integrations", 0] > 0;

DiffExp2`Solve`ClearSolveCaches[];
closedQ = DiffExp2`CppBackend`PersistentSessionInformation[] === <||>;

ok = preparedQ && opaqueBasisQ && planOwnsSCCQ && strongOwnershipQ &&
  transportQ && nativeExecutionQ && cleanupQ && closedQ;

If[TrueQ[ok],
  Print["PASS: mixed regular/singular native atlas execution"],
  Print["FAIL: ", InputForm[<|
    "Atlas" -> If[AssociationQ[atlas],
      KeyTake[atlas, {"Type", "Session", "ContainsSingularReceivingCharts"}],
      atlas],
    "Basis" -> If[AssociationQ[singularBasis],
      KeyTake[singularBasis,
        {"Type", "Session", "NativeSCC", "Dimension"}], singularBasis],
    "PlanStats" -> If[AssociationQ[planStats],
      KeyTake[planStats, {"status", "session", "tile_plan",
        "lower_matches", "upper_matches"}], planStats],
    "Before" -> before, "SCCReleases" -> sccReleases,
    "PlanAfterSCCRelease" -> If[AssociationQ[planAfterSCCRelease],
      KeyTake[planAfterSCCRelease, {"status", "session", "tile_plan"}],
      planAfterSCCRelease],
    "ColumnAfterSCCRelease" -> If[AssociationQ[columnAfterSCCRelease],
      KeyTake[columnAfterSCCRelease, {"status", "session", "local"}],
      columnAfterSCCRelease],
    "Transport" -> If[AssociationQ[transport],
      KeyTake[transport, {"status", "session", "transport_state",
        "matches", "tiles", "final_local"}], transport],
    "TransportStats" -> transportStats,
    "TransportRelease" -> transportRelease,
    "AfterTransport" -> afterTransport,
    "Run" -> If[AssociationQ[run],
      KeyTake[run, {"Type", "Lower", "Upper", "NativeSummary"}], run],
    "Exported" -> If[AssociationQ[exported],
      <|"status" -> Lookup[exported, "status", None],
        "min" -> Lookup[Lookup[exported, "value", <||>], "min", None],
        "max" -> Lookup[Lookup[exported, "value", <||>], "max", None]|>,
      exported], "Decoded" -> decoded,
    "EpsilonZero" -> epsilonZero, "Expected" -> expected,
    "AfterRun" -> afterRun, "Released" -> released, "After" -> after,
    "Checks" -> {preparedQ, opaqueBasisQ, planOwnsSCCQ,
      strongOwnershipQ, transportQ, nativeExecutionQ, cleanupQ, closedQ}|>]];
  Exit[1]];
