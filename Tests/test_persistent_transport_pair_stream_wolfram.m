(* Focused wire/memory-shape coverage for the Wolfram paired tile stream. *)

repo = DirectoryName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{repo, "DiffExp2.m"}]];

passed = 0; failed = 0;
assert[label_String, condition_, detail_:None] := If[TrueQ[condition],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label,
    If[detail === None, "", ": " <> ToString[detail, InputForm]]]];

lowerState = <|"session" -> "mock-session",
  "transport_state" -> "transport:lower",
  "checkpoint_identity" -> "lower-checkpoint",
  "provenance_identity" -> "lower-provenance", "arm" -> "lower",
  "tiles" -> 3,
  "epsilon" -> <|"required_complete_max" -> 2|>|>;
upperState = <|"session" -> "mock-session",
  "transport_state" -> "transport:upper",
  "checkpoint_identity" -> "upper-checkpoint",
  "provenance_identity" -> "upper-provenance", "arm" -> "upper",
  "tiles" -> 2,
  "epsilon" -> <|"required_complete_max" -> 2|>|>;
epsilon = <|"Min" -> -1, "Max" -> 2,
  "RequiredCompleteMax" -> 1|>;
observable = <|"Identity" -> "streamed-observable",
  "CheckpointIdentity" -> "streamed-observable-checkpoint",
  "CoefficientVector" -> {1}, "Epsilon" -> epsilon,
  "TailPolicy" -> "stored", "MinimumEpsilonShift" -> -4|>;
row[id_String] := <|
  "schema" -> "diffexp2-prepared-rational-local-row-v1",
  "columns" -> 1, "exact_identity" -> id, "entries" -> {}|>;
recipe[id_String] := <|"ID" -> id,
  "Shape" -> <|"EpsWindow" -> <|"Min" -> -1,
      "CompleteMax" -> 7|>,
    "TWindow" -> <|"CompleteMax" -> 4|>, "Dimension" -> 1|>|>;
lowerSourceEpsilon = {
  <|"Min" -> -5, "CompleteMax" -> 4|>,
  <|"Min" -> 0, "CompleteMax" -> 2|>,
  <|"Min" -> 0, "CompleteMax" -> 2|>};
upperSourceEpsilon = {
  <|"Min" -> 0, "CompleteMax" -> 2|>,
  <|"Min" -> 0, "CompleteMax" -> 2|>};

requests = {}; cacheClears = 0; preparedWindows = {};
result = Block[{
    DiffExp2`NativeTransport`Private`nativePrepareArmRecipeRow =
      Function[{recipe, ignoredCoefficients, ignoredVar, ignoredDomain},
        AppendTo[preparedWindows, recipe["Shape", "EpsWindow"]];
        row[recipe["ID"]]],
    DiffExp2`NativeTransport`Private`nativeDropRationalMultiplierPreparationCache =
      Function[Null, cacheClears++; Null],
    DiffExp2`CppBackend`RunRequest = Function[request,
      AppendTo[requests, request];
      Switch[request["op"],
        "transport.contract_pair_stream_begin",
          <|"status" -> "ok", "session" -> "mock-session",
            "stream" -> "pair-stream:mock-session:1",
            "stream_checkpoint_identity" -> "stream-checkpoint"|>,
        "transport.contract_pair_stream_add_tile",
          <|"status" -> "ok"|>,
        "transport.contract_pair_stream_finish",
          <|"status" -> "ok", "lines" -> {
            <|"session" -> "mock-session", "line" -> "line:1",
              "checkpoint_identity" ->
                "streamed-observable-checkpoint",
              "provenance_identity" -> "line-provenance",
              "request_index" -> 0,
              "observable_identity" -> "streamed-observable"|>}|>,
        "transport.contract_pair_stream_abort",
          <|"status" -> "ok", "aborted" -> True|>]]},
  DiffExp2`NativeTransport`Private`nativeContractStoredPairObservableStreamed[
    lowerState, upperState, observable,
    recipe /@ {"lower-0", "lower-1", "lower-2"},
    recipe /@ {"upper-0", "upper-1"},
    Global`x, "acb", "mock-root", lowerSourceEpsilon,
    upperSourceEpsilon]];

ops = Lookup[requests, "op"];
adds = Select[requests,
  Lookup[#, "op", None] ===
    "transport.contract_pair_stream_add_tile" &];
begin = First[requests];
assert["stream_begin_is_metadata_only_and_finish_is_atomic",
  Lookup[result, "status", "error"] === "ok" &&
    ops === Join[{"transport.contract_pair_stream_begin"},
      ConstantArray["transport.contract_pair_stream_add_tile", 5],
      {"transport.contract_pair_stream_finish"}] &&
    Sort[Keys[begin["observable"]]] ===
      Sort[{"identity", "checkpoint_identity", "epsilon",
        "tail_policy"}] &&
    FreeQ[begin, "integrand_rows" | "row" | "rows"] &&
    Length[Lookup[result, "lines", {}]] === 1,
  {ops, Lookup[begin, "observable", None]}];
assert["stream_adds_exactly_one_row_in_strict_lower_upper_order",
  Lookup[adds, "side"] ===
    {"lower", "lower", "lower", "upper", "upper"} &&
    Lookup[adds, "tile"] === {0, 1, 2, 0, 1} &&
    Lookup[Lookup[adds, "row"], "exact_identity"] ===
      {"lower-0", "lower-1", "lower-2", "upper-0", "upper-1"} &&
    First[preparedWindows] ===
      <|"Min" -> -1, "CompleteMax" -> 8|> &&
    Rest[preparedWindows] === ConstantArray[
      <|"Min" -> -1, "CompleteMax" -> 7|>, 4] &&
    AllTrue[adds, Sort[Keys[#]] === Sort[{"schema", "op", "session",
      "stream", "stream_checkpoint_identity", "side", "tile", "row"}] &] &&
    cacheClears >= 6,
  {KeyTake[#, {"side", "tile", "row"}] & /@ adds, cacheClears}];

failureRequests = {};
failure = Block[{
    DiffExp2`NativeTransport`Private`nativePrepareArmRecipeRow =
      Function[{recipe, ignoredCoefficients, ignoredVar, ignoredDomain},
        row[recipe["ID"]]],
    DiffExp2`NativeTransport`Private`nativeDropRationalMultiplierPreparationCache =
      Function[Null, Null],
    DiffExp2`CppBackend`RunRequest = Function[request,
      AppendTo[failureRequests, request];
      Switch[request["op"],
        "transport.contract_pair_stream_begin",
          <|"status" -> "ok", "session" -> "mock-session",
            "stream" -> "pair-stream:mock-session:2",
            "stream_checkpoint_identity" -> "stream-checkpoint-2"|>,
        "transport.contract_pair_stream_add_tile",
          If[Count[Lookup[failureRequests, "op"],
                "transport.contract_pair_stream_add_tile"] === 2,
            <|"status" -> "error", "detail" -> "injected"|>,
            <|"status" -> "ok"|>],
        "transport.contract_pair_stream_abort",
          <|"status" -> "ok", "aborted" -> True|>]]},
  DiffExp2`NativeTransport`Private`nativeContractStoredPairObservableStreamed[
    lowerState, upperState, observable,
    recipe /@ {"lower-0", "lower-1", "lower-2"},
    recipe /@ {"upper-0", "upper-1"},
    Global`x, "acb", "mock-failure-root", lowerSourceEpsilon,
    upperSourceEpsilon]];
failureOps = Lookup[failureRequests, "op"];
assert["stream_failure_aborts_once_and_never_finishes",
  Lookup[failure, "status", "ok"] === "error" &&
    Count[failureOps, "transport.contract_pair_stream_abort"] === 1 &&
    FreeQ[failureOps, "transport.contract_pair_stream_finish"],
  failureOps];

DiffExp2`LoadConfiguration[{
  "RecurrenceBackend" -> "Cpp", "WorkingPrecision" -> 60,
  "ChopPrecision" -> 30, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 2, "DivisionOrder" -> 2,
  "Variables" -> {}, "Verbosity" -> 0}];
assert["metadata_epsilon_shift_ignores_structural_zero_components",
  DiffExp2`NativeTransport`Private`nativeIntegrandMinimumShift[
    {0, DiffExp2`Config`CanonicalEps[]^2}, Global`x, 2] === 2 &&
  DiffExp2`NativeTransport`Private`nativeIntegrandMinimumShift[
    {0, 0}, Global`x, 2] === 0];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Print["Some tests FAILED."]; Exit[1]];
