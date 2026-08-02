(* Focused prepare-only native regular chart-owner regression. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

If[!TrueQ[DiffExp2`CppBackend`BackendAvailableQ[]],
  Print["SKIP: compiled DiffExp2 backend is not available"];
  Exit[If[Environment["DE2_REQUIRE_CPP"] === "1", 1, 0]]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];
SetAttributes[catchDE2, HoldFirst];
catchDE2[expression_] := Quiet[Catch[expression, "DiffExp2Error"]];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 1, "DivisionOrder" -> 2,
  "Variables" -> {}, "Verbosity" -> 0}];

x = Global`x;
system = DiffExp2`LoadSystem[<|
  "Matrix" -> {{0, 1}, {1, 0}}, "Variable" -> x|>];
plan = DiffExp2`Transport`SegmentLine[system, {0, 1/4}];
chartSystem = DiffExp2`Solve`PrepareChart[system, First[plan["Charts"]]];
request = <|"EpsWindow" -> <|"Min" -> -1, "CompleteMax" -> 1|>,
  "TOrder" -> 10|>;
zero = DiffExp2`EpsSeries`ESNew[-1, {0, 0, 0}];
one = DiffExp2`EpsSeries`ESNew[-1, {0, 1, 0}];

sessionStats[session_] := Lookup[
  DiffExp2`CppBackend`PersistentSessionInformation[], session, <||>];
ownerChartStats[stats_, chart_] := SelectFirst[
  Lookup[stats, "chart_stats", {}],
  Lookup[#, "chart", None] === chart &, <||>];

localSolveCalls = 0;
owner = catchDE2[Block[{
    DiffExp2`CppBackend`RunPersistentLocalSolve =
      Function[{ignoredRequest, ignoredMetadata, ignoredLocalMetadata},
        localSolveCalls++;
        Failure["UnexpectedLocalSolve", <||>]]},
  DiffExp2`Solve`PrepareNativeRegularBasisOwner[chartSystem, request]]];
ownerStats = If[AssociationQ[owner], sessionStats[owner["Session"]], <||>];
ownerChart = If[AssociationQ[owner],
  ownerChartStats[ownerStats, owner["NativeChart"]], <||>];

assert["prepare_only_regular_owner_publishes_no_witness_local",
  AssociationQ[owner] &&
    Lookup[owner, "Type", None] === "DiffExp2NativeRegularBasisOwner" &&
    localSolveCalls === 0 &&
    FreeQ[Keys[owner], "Witness" | "Local"] &&
    StringQ[Lookup[owner, "ChartIdentity", None]] &&
    Lookup[ownerStats, "charts", -1] === 1 &&
    Lookup[ownerStats, "locals", -1] === 0 &&
    Lookup[ownerStats, "local_solves", -1] === 0 &&
    Lookup[ownerStats, "pending_local_solves", -1] === 0 &&
    Lookup[ownerChart, "local_solves", -1] === 0];

(* The backend chart seam must not regain a dependency on per-run fields.
   Reconstruct the same certified capture, remove every dynamic field, and
   prove that static-only preparation touches the existing owner without a
   solve. *)
runKeys = {"nmax", "p", "has_initial", "adaptive_probe", "a_target",
  "b_target", "a_shift_min", "a_shifts", "schedule", "initial",
  "initial_validity", "source", "return_u"};
physical = DiffExp2`Solve`Private`regularPhysicalChartSystem[chartSystem];
capture = catchDE2[
  DiffExp2`Solve`Private`prepareNativeLocalFamilyRun[physical, request,
    <|"a" -> 0, "b" -> 0, "p" -> 0|>, {{one, zero}}]];
staticOwner = If[AssociationQ[capture],
  DiffExp2`CppBackend`PreparePersistentChart[
    KeyDrop[capture["Request"], runKeys],
    capture["PersistentMetadata"]], capture];
staticStats = If[AssociationQ[owner], sessionStats[owner["Session"]], <||>];
assert["static_chart_preparation_is_independent_of_dynamic_run_fields",
  AssociationQ[staticOwner] && AssociationQ[owner] &&
    staticOwner["Session"] === owner["Session"] &&
    staticOwner["Chart"] === owner["NativeChart"] &&
    staticOwner["ChartIdentity"] === owner["ChartIdentity"] &&
    Lookup[staticStats, "charts", -1] === 1 &&
    Lookup[staticStats, "locals", -1] === 0 &&
    Lookup[staticStats, "local_solves", -1] === 0];

value = If[FailureQ[owner], owner, catchDE2[
  DiffExp2`Solve`SolveNativeValueRegular[
    chartSystem, request, {one, zero}]]];
valueStats = If[AssociationQ[owner], sessionStats[owner["Session"]], <||>];
valueChart = If[AssociationQ[owner],
  ownerChartStats[valueStats, owner["NativeChart"]], <||>];

assert["later_native_value_run_reuses_prepare_only_owner_chart",
  AssociationQ[value] && AssociationQ[owner] &&
    value["Session"] === owner["Session"] &&
    value["NativeChart"] === owner["NativeChart"] &&
    Lookup[valueStats, "charts", -1] === 1 &&
    Lookup[valueStats, "locals", -1] === 1 &&
    Lookup[valueStats, "local_solves", -1] === 1 &&
    Lookup[valueChart, "local_solves", -1] === 1];

valueRelease = If[AssociationQ[value],
  DiffExp2`CppBackend`ReleasePersistentLocal[value], value];
basis = If[FailureQ[owner], owner, catchDE2[
  DiffExp2`Solve`SolveNativeRegularBasis[
    chartSystem, request, 1]]];
basisStats = If[AssociationQ[owner], sessionStats[owner["Session"]], <||>];
basisChart = If[AssociationQ[owner],
  ownerChartStats[basisStats, owner["NativeChart"]], <||>];

assert["fallback_monolithic_basis_reuses_prepare_only_owner_chart",
  AssociationQ[basis] && AssociationQ[owner] &&
    Lookup[basis, "Type", None] === "DiffExp2NativeRegularBasis" &&
    basis["Session"] === owner["Session"] &&
    basis["NativeChart"] === owner["NativeChart"] &&
    Lookup[basis["Columns"], "NativeChart", {}] ===
      ConstantArray[owner["NativeChart"], 2] &&
    Lookup[valueRelease, "status", "error"] === "ok" &&
    Lookup[basisStats, "charts", -1] === 1 &&
    Lookup[basisStats, "locals", -1] === 2 &&
    Lookup[basisStats, "local_solves", -1] === 3 &&
    Lookup[basisChart, "local_solves", -1] === 3];

If[AssociationQ[basis],
  Scan[DiffExp2`CppBackend`ReleasePersistentLocal, basis["Columns"]]];

wideFrameCalls = 0;
lightOwner = catchDE2[Block[{
    DiffExp2`Solve`Private`prepareNativeLocalFamilyRun =
      Function[Null, wideFrameCalls++;
        Failure["UnexpectedWideFrame", <||>]]},
  DiffExp2`Solve`PrepareNativeRegularBasisOwner[
    chartSystem, request, owner]]];
lightStats = If[AssociationQ[lightOwner],
  sessionStats[lightOwner["Session"]], <||>];
assert["anchor_session_owner_is_frame_independent_and_frame_free",
  AssociationQ[lightOwner] &&
    lightOwner["Session"] === owner["Session"] &&
    Lookup[lightOwner, "OwnerKind", None] ===
      "FrameIndependentRegularPhysicalEquation" &&
    StringStartsQ[Lookup[lightOwner, "NativeEquationOwner", ""],
      "eq:"] &&
    Lookup[lightOwner["ValueSolver"], "schema", None] ===
      "diffexp2-native-ordinary-physical-value-solver-v1" &&
    Lookup[lightOwner["PhysicalPreparation"], "Schema", None] ===
      "DiffExp2RegularOwnerPhysicalPreparation/v1" &&
    lightOwner["PhysicalPreparation", "EquationIdentity"] ===
      lightOwner["EquationIdentity"] &&
    lightOwner["PhysicalPreparation", "PhysicalPayloadIdentity"] ===
      lightOwner["PhysicalPayloadIdentity"] &&
    Sort[Keys[lightOwner["ValueSolver"]]] === Sort[{
      "schema", "taylor_complete_max", "metadata",
      "relative_accuracy_max_exact"}] &&
    wideFrameCalls === 0 &&
    Lookup[lightStats, "regular_equation_owners", -1] === 1 &&
    Lookup[lightStats, "charts", -1] === 1 &&
    Lookup[lightStats, "locals", -1] === 0];

secondClearCalls = 0;
lightBasis = If[FailureQ[lightOwner], lightOwner, catchDE2[Block[{
    DiffExp2`Solve`Private`clearedSymbolic =
      Function[Null, secondClearCalls++;
        Failure["UnexpectedSecondChartClear", <||>]]},
  DiffExp2`Solve`SolveNativeRegularBasis[
    chartSystem, request, 1, True, lightOwner]]]];
lightBasisStats = If[AssociationQ[lightOwner],
  sessionStats[lightOwner["Session"]], <||>];
assert["framed_fallback_is_published_under_frame_independent_owner",
  AssociationQ[lightBasis] && AssociationQ[lightOwner] &&
    secondClearCalls === 0 &&
    lightBasis["Session"] === lightOwner["Session"] &&
    lightBasis["NativeChart"] === lightOwner["NativeEquationOwner"] &&
    Lookup[lightBasis["Columns"], "NativeChart", {}] ===
      ConstantArray[lightOwner["NativeEquationOwner"], 2] &&
    Lookup[lightBasisStats, "regular_equation_owners", -1] === 1 &&
    Lookup[lightBasisStats, "locals", -1] === 2];
If[AssociationQ[lightBasis],
  Scan[DiffExp2`CppBackend`ReleasePersistentLocal,
    lightBasis["Columns"]]];
lightPreparedRelease = If[AssociationQ[lightBasis],
  DiffExp2`CppBackend`ReleasePersistentPreparedToken[
    lightBasis["PreparedToken"]], lightBasis];
afterLightPreparedRelease = If[AssociationQ[lightOwner],
  sessionStats[lightOwner["Session"]], <||>];
assert["consumed_regular_basis_token_releases_transient_recurrence_chart",
  AssociationQ[lightBasis] &&
    StringQ[Lookup[lightBasis, "PreparedToken", None]] &&
    lightPreparedRelease === Null &&
    Lookup[afterLightPreparedRelease, "charts", -1] === 0 &&
    Lookup[afterLightPreparedRelease, "locals", -1] === 0 &&
    Lookup[afterLightPreparedRelease, "regular_equation_owners", -1] === 1];
tamperedOwner = If[AssociationQ[lightOwner], Join[lightOwner, <|
  "PhysicalPreparation" -> Join[lightOwner["PhysicalPreparation"], <|
    "Request" -> Join[request, <|"TOrder" -> request["TOrder"] + 1|>]|>]|>],
  lightOwner];
tamperedOwnerBasis = If[AssociationQ[tamperedOwner], catchDE2[
  DiffExp2`Solve`SolveNativeRegularBasis[
    chartSystem, request, 1, True, tamperedOwner]], tamperedOwner];
assert["tampered retained physical preparation is rejected before solving",
  FailureQ[tamperedOwnerBasis] && AssociationQ[lightOwner] &&
    Lookup[sessionStats[lightOwner["Session"]], "locals", -1] === 0];
lightOwnerRelease = If[AssociationQ[lightOwner],
  DiffExp2`CppBackend`ReleasePersistentRegularEquationOwner[lightOwner],
  lightOwner];
afterLightRelease = If[AssociationQ[owner],
  sessionStats[owner["Session"]], <||>];
assert["frame_independent_owner_public_token_releases_cleanly",
  AssociationQ[lightOwnerRelease] &&
    Lookup[lightOwnerRelease, "status", "error"] === "ok" &&
    Lookup[afterLightRelease, "regular_equation_owners", -1] === 0 &&
    Lookup[afterLightRelease, "locals", -1] === 0];

DiffExp2`Solve`ClearSolveCaches[];
DiffExp2`CppBackend`ClearPersistentSessions[];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1],
  Print["All tests PASSED!"]; Exit[0]];
