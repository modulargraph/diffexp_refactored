(* Replay one diagnostic native checkpoint saved immediately before the final
   singular transport hop.  This is intentionally opt-in: it exists to make
   terminal-frame development deterministic without remarching every earlier
   chart.

   Required environment:
     DE2_TERMINAL_REPLAY_MANIFEST=/absolute/path/prehop.de2cp.manifest.wl

   DE2_CPP_LIBRARY must select the library build under test. *)

replayManifestPath = Environment["DE2_TERMINAL_REPLAY_MANIFEST"];
If[!StringQ[replayManifestPath] ||
    StringLength[StringTrim[replayManifestPath]] == 0,
  Print["DE2_TERMINAL_REPLAY_MANIFEST is required"];
  Exit[2]];
replayManifestPath = ExpandFileName[replayManifestPath];
If[!FileExistsQ[replayManifestPath],
  Print["terminal replay manifest does not exist: ",
    replayManifestPath];
  Exit[2]];

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

manifest = Get[replayManifestPath];
postHopReplay = TrueQ[Lookup[manifest, "PostHop", False]];
requiredManifestKeys = Join[{
  "CheckpointPath", "CheckpointIdentity", "Plan", "Arm",
  "Anchor", "TileSources", "Epsilon", "CheckpointRoot", "Refinement",
  "PreparedObservables", "RowRecipes", "Variable", "Domain",
  "ObservableCheckpointRoot"},
  If[postHopReplay, {}, {"Index", "Basis", "Incoming"}]];
If[!AssociationQ[manifest] ||
    !AllTrue[requiredManifestKeys, KeyExistsQ[manifest, #] &],
  Print["terminal replay manifest is incomplete: ",
    InputForm[If[AssociationQ[manifest],
      Complement[requiredManifestKeys, Keys[manifest]], manifest]]];
  Exit[2]];
hopEpsilon = manifest["Epsilon"];
(* Historical diagnostic manifests saved the three-field per-hop epsilon
   contract.  State publication uses the four-field batch contract; its
   match requirement is necessarily the already-maximized public requirement
   for this batch. *)
publishEpsilon = If[
  KeyExistsQ[hopEpsilon, "match_required_complete_max"],
  hopEpsilon,
  Append[hopEpsilon, "match_required_complete_max" ->
    hopEpsilon["required_complete_max"]]];
Print["TERMINAL REPLAY MANIFEST LOADED"];

replayCheckpointIdentity[handle_Association] :=
  Lookup[handle, "checkpoint_identity",
    Lookup[handle, "CheckpointIdentity", None]];
replayPlanHandle[handle_Association, session_String] := Module[
  {plan = Lookup[handle, "tile_plan",
      Lookup[handle, "TilePlan", None]],
   checkpoint = replayCheckpointIdentity[handle]},
  If[!StringQ[plan] || !StringQ[checkpoint] ||
      StringLength[checkpoint] == 0,
    Return[Failure["TerminalReplay", <|
      "Detail" -> "saved tile-plan handle is incomplete"|>], Module]];
  <|"session" -> session, "tile_plan" -> plan,
    "checkpoint_identity" -> checkpoint|>];
replayLocalHandle[handle_Association, session_String] := Module[
  {local = Lookup[handle, "local", Lookup[handle, "Local", None]],
   checkpoint = replayCheckpointIdentity[handle]},
  If[!StringQ[local] || !StringQ[checkpoint] ||
      StringLength[checkpoint] == 0,
    Return[Failure["TerminalReplay", <|
      "Detail" -> "saved local handle is incomplete"|>], Module]];
  <|"session" -> session, "local" -> local,
    "checkpoint_identity" -> checkpoint|>];

restored = DiffExp2`CppBackend`RestorePersistentCheckpoint[
  manifest["CheckpointPath"], manifest["CheckpointIdentity"]];
If[FailureQ[restored] || !AssociationQ[restored] ||
    Lookup[restored, "status", "error"] =!= "ok" ||
    !StringQ[Lookup[restored, "session", None]],
  Print["terminal replay checkpoint restore failed: ",
    InputForm[restored]];
  Exit[3]];
restoredSession = restored["session"];
Print["TERMINAL REPLAY CHECKPOINT RESTORED session=", restoredSession];
(* Opaque object tokens are stable across checkpoint restore, but every
   restored object belongs to the newly allocated session.  Reconstruct the
   tiny handle envelopes instead of recursively rewriting the saved
   associations: Association keys are not Rule expressions, and the saved
   handles can contain very large diagnostic summaries. *)
planHandle = replayPlanHandle[manifest["Plan"], restoredSession];
anchorHandle = replayLocalHandle[manifest["Anchor"], restoredSession];
tileSourceHandles =
  replayLocalHandle[#, restoredSession] & /@ manifest["TileSources"];
basisHandles =
  If[postHopReplay, {},
    replayLocalHandle[#, restoredSession] & /@ manifest["Basis"]];
incomingHandle = If[postHopReplay, None,
  replayLocalHandle[manifest["Incoming"], restoredSession]];
handleFailures = Cases[
  {planHandle, anchorHandle, tileSourceHandles, basisHandles,
    incomingHandle}, _Failure, Infinity];
If[handleFailures =!= {},
  Print["terminal replay handle rebinding failed: ",
    InputForm[First[handleFailures]]];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
  Exit[3]];
manifest = KeyDrop[
  manifest, {"Plan", "Anchor", "TileSources", "Basis", "Incoming"}];

If[postHopReplay,
  Print["TERMINAL REPLAY POSTHOP STATE LOADED"];
  tileSources = tileSourceHandles,
  hop = DiffExp2`CppBackend`ConsumePersistentTransportHop[
    planHandle, manifest["Arm"], manifest["Index"],
    basisHandles, incomingHandle, hopEpsilon,
    manifest["CheckpointRoot"], manifest["Refinement"]];
  If[FailureQ[hop] || !AssociationQ[hop] ||
      Lookup[hop, "status", "error"] =!= "ok" ||
      !AssociationQ[Lookup[hop, "next_local", None]],
    Print["terminal replay hop failed: ", InputForm[hop]];
    Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
    Exit[4]];
  Print["TERMINAL REPLAY FINAL HOP CONSUMED"];
  matchDiagnostics = Lookup[hop["match_reference"],
    "diagnostic_native_match_summary", <||>];
  If[AssociationQ[matchDiagnostics] && matchDiagnostics =!= <||>,
    compactMatchDiagnostics = KeyTake[matchDiagnostics, {
      "basis_point_exact", "incoming_point_exact",
      "physical_match_point_exact", "epsilon", "refinement",
      "weight_windows", "transformed_weight_windows",
      "normal_frame_attempt", "materialization_association"}];
    compactMatchDiagnostics = Join[compactMatchDiagnostics, <|
      "matching_taylor_widths" -> Lookup[
        Lookup[matchDiagnostics, "basis", {}],
        "matching_taylor_width", {}],
      "exact_lattice" -> KeyTake[
        Lookup[matchDiagnostics, "exact_lattice", <||>], {
          "schema", "witness_schema", "canonical_witness_bytes",
          "transformation_terms", "transformation_min_power",
          "initial_column_shifts", "normalized_determinant_valuation",
          "initial_leading_rank", "final_leading_rank",
          "saturation_actions"}],
      "residual" -> KeyTake[
        Lookup[matchDiagnostics, "residual", <||>], {
          "verdict", "complete_window", "required_complete_max",
          "complete_through_required", "coefficient_diagnostics",
          "coefficient_verdicts", "detail", "scope",
          "frame_identity"}]|>];
    Print["TERMINAL REPLAY MATCH DIAGNOSTICS ", InputForm[
      compactMatchDiagnostics]]];

  nextLocal = hop["next_local"];
  If[!KeyExistsQ[nextLocal, "session"],
    nextLocal = Append[nextLocal, "session" -> restoredSession]];
  tileSources = Append[tileSourceHandles, nextLocal]];
published =
  DiffExp2`CppBackend`PublishPersistentConsumedTransportState[
    planHandle, anchorHandle, manifest["Arm"], tileSources,
    publishEpsilon, manifest["CheckpointRoot"],
    manifest["Refinement"]];
state = If[AssociationQ[published],
  Lookup[Lookup[published, "states", <||>], manifest["Arm"], None],
  None];
If[FailureQ[published] || !AssociationQ[published] ||
    Lookup[published, "status", "error"] =!= "ok" ||
    !AssociationQ[state],
  Print["terminal replay state publication failed: ",
    InputForm[published]];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
  Exit[5]];
Print["TERMINAL REPLAY STATE PUBLISHED"];
stateStats =
  DiffExp2`CppBackend`PersistentTransportArmStatistics[state];
If[Environment["DE2_DIAGNOSTIC_TERMINAL_STATE"] === "1",
  If[KeyExistsQ[stateStats, "terminal_diagnostic"],
    terminalDiagnostic = stateStats["terminal_diagnostic"];
    compactTerminalDiagnostic = <|
      "Schema" -> Lookup[terminalDiagnostic, "schema", None],
      "Match" -> KeyTake[
        Lookup[terminalDiagnostic, "match", <||>], {
          "basis_point_exact", "incoming_point_exact",
          "physical_match_point_exact", "epsilon",
          "materialization_association", "residual"}],
      "Columns" -> Lookup[terminalDiagnostic, "columns", {}]|>;
    Print["TERMINAL REPLAY STATE DIAGNOSTICS ",
      InputForm[compactTerminalDiagnostic]],
    Print["TERMINAL REPLAY STATE DIAGNOSTICS UNAVAILABLE ",
      InputForm[KeyTake[stateStats, {
        "transport_state", "arm", "tiles", "terminal_factorized_match",
        "epsilon", "final_local"}]]]]];

observables = Select[manifest["PreparedObservables"],
  Lookup[#, "Operation", None] === "integrate" &];
If[Length[observables] =!= 1 ||
    Length[manifest["RowRecipes"]] =!=
      Lookup[state, "tiles", -1],
  Print["terminal replay currently requires exactly one integrate ",
    "observable and one row recipe per tile"];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
  Exit[6]];
observable = First[observables];
consumerContracts =
  DiffExp2`NativeTransport`Private`nativeTransportTileConsumerEpsilonContracts[
    state];
rows = MapThread[
  Function[{recipe, contract},
    DiffExp2`NativeTransport`Private`nativePrepareArmRecipeRow[
      DiffExp2`NativeTransport`Private`nativeConsumerRowRecipe[
        recipe, contract, observable, 1],
      observable["CoefficientVector"], manifest["Variable"],
      manifest["Domain"]]],
  {manifest["RowRecipes"], consumerContracts}];

contractObservable = Join[<|
    "Identity" -> observable["Identity"] <> ":terminal-replay:" <>
      manifest["Arm"],
    "CheckpointIdentity" ->
      observable["CheckpointIdentity"] <> ":terminal-replay:" <>
        manifest["Arm"],
    "IntegrandRows" -> rows,
    "Epsilon" -> observable["Epsilon"],
    "TailPolicy" -> observable["TailPolicy"]|>,
  If[KeyExistsQ[observable, "DivergentCancellation"],
    <|"DivergentCancellation" ->
      observable["DivergentCancellation"]|>, <||>]];
contracted =
  DiffExp2`CppBackend`ContractPersistentTransportObservables[
    state, {contractObservable},
    manifest["ObservableCheckpointRoot"] <> ":terminal-replay:" <>
      manifest["Arm"]];
line = If[AssociationQ[contracted] &&
    Length[Lookup[contracted, "lines", {}]] == 1,
  First[contracted["lines"]], None];
If[FailureQ[contracted] || !AssociationQ[contracted] ||
    Lookup[contracted, "status", "error"] =!= "ok" ||
    !AssociationQ[line],
  Print["terminal replay contraction failed: ",
    InputForm[contracted]];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
  Exit[7]];

exported = DiffExp2`CppBackend`ExportPersistentLineIntegral[
  line, line["checkpoint_identity"], 50];
If[FailureQ[exported] || !AssociationQ[exported] ||
    Lookup[exported, "status", "error"] =!= "ok",
  Print["terminal replay export failed: ", InputForm[exported]];
  Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
  Exit[8]];

Print["TERMINAL REPLAY ARM ", manifest["Arm"]];
Print["TERMINAL REPLAY VALUE ", InputForm[exported["value"]]];
Print["TERMINAL REPLAY LINE ", InputForm[
  KeyTake[First[contracted["lines"]], {
    "request_index", "observable_identity", "line",
    "checkpoint_identity", "scope", "epsilon", "tiles",
    "elapsed_ms"}]]];
Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
Exit[0];
