(* Focused Wolfram orchestration checks for the retained regular-basis batch. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

(* Exercise the public CppBackend wrapper without opening a native session:
   the private preparation seam is replaced by an exact opaque chart token,
   leaving the protocol validation and payload construction under test. *)
staticKeys = DiffExp2`CppBackend`Private`$persistentRequiredStaticKeys;
runKeys = DiffExp2`CppBackend`Private`$persistentRunKeys;
baseRequest = AssociationThread[Join[staticKeys, runKeys],
  Range[Length[staticKeys] + Length[runKeys]]];
requests = {Append[baseRequest, "nmax" -> 3],
  Append[baseRequest, "nmax" -> 4]};
localMetadata = {<|"checkpoint_identity" -> "column-1"|>,
  <|"checkpoint_identity" -> "column-2"|>};
capturedPayload = None;
wrapperResult = Block[{
    DiffExp2`CppBackend`Private`preparePersistentRequest =
      Function[{request, metadata}, <|"Session" -> "s:batch",
        "Chart" -> "c:batch",
        "Run" -> KeyTake[request, runKeys]|>],
    DiffExp2`CppBackend`Private`persistentRegularEquationOwnerHandles =
      Function[owner, <|"Session" -> "s:batch",
        "EquationOwner" -> "eq:batch"|>],
    DiffExp2`CppBackend`RunRequest = Function[payload,
      capturedPayload = payload; <|"status" -> "ok"|>]},
  DiffExp2`CppBackend`RunPersistentLocalSolves[requests, <||>,
    localMetadata, 2, <|"Session" -> "s:batch",
      "NativeEquationOwner" -> "eq:batch"|>]];

assert["backend_wrapper_emits_one_owner_bound_ordered_batch",
  Lookup[wrapperResult, "status", "error"] === "ok" &&
    AssociationQ[capturedPayload] &&
    capturedPayload["op"] === "local.solve_batch" &&
    capturedPayload["session"] === "s:batch" &&
    capturedPayload["chart"] === "c:batch" &&
    capturedPayload["equation_owner"] === "eq:batch" &&
    capturedPayload["threads"] === 2 &&
    capturedPayload["runs"] ===
      (KeyTake[#, runKeys] & /@ requests) &&
    capturedPayload["metadata"] === localMetadata];

badWrapper = DiffExp2`CppBackend`RunPersistentLocalSolves[
  requests, <||>, Take[localMetadata, 1], 2];
assert["backend_wrapper_rejects_nonparallel_metadata_before_native_call",
  FailureQ[badWrapper]];

(* Isolate SolveNativeRegularBasis orchestration.  The expensive shared and
   lightweight dynamic preparation seams are counted independently; final
   native summaries are represented by opaque ordered fixture handles. *)
d = 3;
cs = <|"SystemSize" -> d,
  "IndicialData" -> <|"Regular" -> True|>,
  "Center" -> 0, "ChartMap" -> <||>, "Radius" -> 2,
  "Prescriptions" -> {}|>;
req = <|"EpsWindow" -> <|"Min" -> -1, "CompleteMax" -> 1|>,
  "TOrder" -> 5|>;
equationOwner = <|"Session" -> "s:basis",
  "NativeEquationOwner" -> "eq:basis",
  "EquationIdentity" -> "de2-equation-basis"|>;
sharedCalls = 0; dynamicCalls = 0; batchCalls = 0;
recordedInitials = {}; recordedThreads = None; recordedOwner = None;

basis = Block[{
    DiffExp2`Solve`Private`regularPhysicalChartSystem = Function[input, input],
    DiffExp2`Solve`Private`prepareNativeLocalFamilyShared =
      Function[{physical, request, tag, prototype}, sharedCalls++;
        <|"Shared" -> True|>],
    DiffExp2`Solve`Private`prepareNativeLocalFamilyDynamic =
      Function[{shared, init, retain}, dynamicCalls++;
        AppendTo[recordedInitials, init];
        <|"Dimension" -> d, "RequestedMin" -> -1,
          "RequestedMax" -> 1,
          "Tag" -> <|"a" -> 0, "b" -> 0, "p" -> 0|>,
          "Request" -> <|"column" -> dynamicCalls|>,
          "PersistentMetadata" -> <|"PreparedToken" -> "operator"|>,
          "LocalMetadata" ->
            <|"checkpoint_identity" -> ToString[dynamicCalls]|>,
          "CheckpointIdentity" -> ToString[dynamicCalls]|>],
    DiffExp2`CppBackend`RunPersistentLocalSolves,
    DiffExp2`Solve`Private`nativeLocalFamilyFinalize =
      Function[{physical, request, spec, response, session, chart, operator},
        <|"Type" -> "DiffExp2NativeLocalFamily",
          "Session" -> response["session"], "Local" -> response["local"],
          "NativeChart" -> response["chart"]|>]},
  DiffExp2`CppBackend`RunPersistentLocalSolves[
      nativeRequests_, metadata_, nativeMetadata_, threads_, owner_] :=
    (batchCalls++; recordedThreads = threads; recordedOwner = owner;
      <|"status" -> "ok", "session" -> "s:basis",
        "chart" -> "c:framed", "results" ->
          Table[<|"status" -> "ok", "session" -> "s:basis",
            "local" -> ("l:" <> ToString[index]),
            "chart" -> "eq:basis", "column" -> index|>,
            {index, d}],
        "attempted" -> d, "succeeded" -> d, "failed" -> 0,
        "requested_threads" -> threads, "worker_threads" -> threads,
        "thread_limit" -> 64, "atomic_retention" -> True,
        "json_coefficients" -> 0, "elapsed_ms" -> 1.|>);
  DiffExp2`Solve`SolveNativeRegularBasis[
    cs, req, 3, True, equationOwner]];

unitMatrix = Table[
  DiffExp2`EpsSeries`ESCoefficient[
    recordedInitials[[basisIndex, 1, component]], 0],
  {basisIndex, d}, {component, d}];
assert["regular_basis_hoists_immutable_preparation_once",
  AssociationQ[basis] && sharedCalls === 1 && dynamicCalls === d &&
    batchCalls === 1 && recordedThreads === 3 &&
    recordedOwner === equationOwner && unitMatrix === IdentityMatrix[d]];
assert["regular_basis_reports_native_parallel_batch_and_owner_handles",
  basis["NativeChart"] === "eq:basis" &&
    basis["PreparedToken"] === "operator" &&
    Lookup[basis["Columns"], "BasisIndex", {}] === Range[d] &&
    basis["NativeSummary", "worker_threads"] === 3 &&
    TrueQ[basis["NativeSummary", "atomic_retention"]] &&
    basis["NativeSummary", "execution_capability"] ===
      "retained-regular-monolithic-unit-basis-v2" &&
    FreeQ[basis, "assembled" | "coefficients" | "u" | "validity"]];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Exit[1], Exit[0]];
