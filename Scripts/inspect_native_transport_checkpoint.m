(* Restore a completed Feynman-trick native transport ladder checkpoint and
   print compact terminal-mode diagnostics for its retained arm states.

   Required environment:
     DE2_NATIVE_LADDER_CHECKPOINT=/absolute/path/example_levelN_native_transport.mx
     DE2_CPP_LIBRARY=/absolute/path/diffexp2_librarylink
*)

ladderPath = Environment["DE2_NATIVE_LADDER_CHECKPOINT"];
If[!StringQ[ladderPath] || StringLength[StringTrim[ladderPath]] == 0,
  Print["DE2_NATIVE_LADDER_CHECKPOINT is required"];
  Exit[2]];
ladderPath = ExpandFileName[ladderPath];
If[!FileExistsQ[ladderPath],
  Print["native ladder checkpoint does not exist: ", ladderPath];
  Exit[2]];

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
SetEnvironment["FT_RUNNER_DEFINITIONS_ONLY" -> "1"];
Get[FileNameJoin[{repoRoot, "Scripts", "run_ft_stepwise2.m"}]];

Clear[Global`$FT2LadderCheckpoint];
Get[ladderPath];
ladder = Global`$FT2LadderCheckpoint;
If[!AssociationQ[ladder] ||
    Lookup[ladder, "Kind", None] =!= "NativeTransport" ||
    !AssociationQ[Lookup[ladder, "NativeTransportCheckpoint", None]],
  Print["native ladder checkpoint payload is malformed"];
  Exit[2]];
stateManifest = ladder["NativeTransportCheckpoint", "State"];
If[!AssociationQ[stateManifest] ||
    !AssociationQ[Lookup[stateManifest, "StateHandles", None]],
  Print["native ladder checkpoint has no completed state manifest"];
  Exit[2]];

restored = DiffExp2`CppBackend`RestorePersistentCheckpoint[
  stateManifest["Path"], stateManifest["CheckpointIdentity"]];
If[FailureQ[restored] || !AssociationQ[restored] ||
    Lookup[restored, "status", "error"] =!= "ok",
  Print["native checkpoint restore failed: ", InputForm[restored]];
  Exit[3]];
session = restored["session"];
Print["NATIVE CHECKPOINT RESTORED session=", session];

compactTerminal[stats_Association] := Module[
  {terminal, match, residual, finalLocal},
  terminal = Lookup[stats, "terminal_diagnostic", None];
  finalLocal = KeyTake[
    Lookup[stats, "final_local", <||>], {
      "local", "chart", "checkpoint_identity"}];
  If[!AssociationQ[terminal],
    Return[<|"transport_state" -> Lookup[stats, "transport_state", None],
      "arm" -> Lookup[stats, "arm", None],
      "tiles" -> Lookup[stats, "tiles", None],
      "terminal_factorized_match" ->
        Lookup[stats, "terminal_factorized_match", None],
      "epsilon" -> Lookup[stats, "epsilon", None],
      "final_local" -> finalLocal|>, Module]];
  match = Lookup[terminal, "match", <||>];
  residual = KeyTake[Lookup[match, "residual", <||>], {
    "verdict", "complete_window", "required_complete_max",
    "complete_through_required", "coefficient_diagnostics",
    "coefficient_verdicts", "detail", "scope"}];
  weightSummary[weight_Association] := Module[
    {coefficients, logarithmicScales},
    coefficients = Lookup[weight, "coefficients", {}];
    logarithmicScales = Map[
      Function[coefficient, Module[{upper},
        upper = Lookup[coefficient, "absolute_upper", 0];
        <|"power" -> Lookup[coefficient, "power", Missing["power"]],
          "log10_absolute_upper" ->
            If[TrueQ[upper == 0], "-Infinity",
              ToString[
                NumberForm[N[Log[10, upper], 8],
                  {10, 3}, NumberPadding -> {"", "0"}],
                InputForm]]|>]],
      coefficients];
    <|"min" -> Lookup[weight, "min", None],
      "max" -> Lookup[weight, "max", None],
      "log10_absolute_upper_by_power" -> logarithmicScales|>];
  compactColumn[column_Association] :=
    <|"column" -> Lookup[column, "column", None],
      "local" -> Lookup[column, "local", None],
      "epsilon" -> Lookup[column, "epsilon", None],
      "taylor_complete_max" ->
        Lookup[column, "taylor_complete_max", None],
      "sectors" -> Lookup[column, "sectors", {}],
      "physical_weight" ->
        weightSummary[Lookup[column, "physical_weight", <||>]],
      "transformed_weight" ->
        weightSummary[Lookup[column, "transformed_weight", <||>]]|>;
  compactMode[mode_Association] := Join[
    KeyDrop[mode, {
      "absolute_upper_approx", "coefficient_scales"}],
    <|"coefficient_scales" ->
      Map[KeyDrop[#, {"absolute_upper_approx"}] &,
        Lookup[mode, "coefficient_scales", {}]]|>];
  <|"State" -> <|
       "transport_state" -> Lookup[stats, "transport_state", None],
       "arm" -> Lookup[stats, "arm", None],
       "tiles" -> Lookup[stats, "tiles", None],
       "terminal_factorized_match" ->
         Lookup[stats, "terminal_factorized_match", None],
       "epsilon" -> Lookup[stats, "epsilon", None],
       "final_local" -> finalLocal|>,
    "Match" -> KeyTake[match, {
       "basis_point_exact", "incoming_point_exact",
       "physical_match_point_exact", "epsilon",
       "matching_frame_identity", "residual_frame_identity",
       "exact_lattice", "refinement", "normal_frame_attempt",
       "materialization_association"}],
    "Residual" -> residual,
    "OwnerNormalFrameDiagnostic" ->
      Lookup[terminal, "owner_normal_frame_diagnostic", None],
    "Columns" ->
      Map[compactColumn, Lookup[terminal, "columns", {}]],
    "PhysicalEndpointModes" ->
      Map[compactMode,
        Lookup[terminal, "physical_endpoint_modes", {}]]|>];

terminalReplays = <||>;
restoredStateHandles = <||>;
Scan[Function[side, Module[{saved, handle, stats},
    saved = stateManifest["StateHandles", side];
    handle = <|"session" -> session,
      "transport_state" -> saved["Handle"],
      "checkpoint_identity" -> saved["CheckpointIdentity"]|>;
    restoredStateHandles[side] = handle;
    stats = DiffExp2`CppBackend`PersistentTransportArmStatistics[handle];
    If[FailureQ[stats] || !AssociationQ[stats] ||
        Lookup[stats, "status", "error"] =!= "ok",
      Print["NATIVE CHECKPOINT STATE FAIL side=", side, " ",
        InputForm[stats]],
      Print["NATIVE CHECKPOINT STATE side=", side, " ",
        InputForm[
          If[Environment["DE2_REBUILD_TERMINAL_MATCH"] === "1",
            KeyTake[compactTerminal[stats], {
              "State", "Match", "Residual",
              "OwnerNormalFrameDiagnostic"}],
            compactTerminal[stats]]]];
      If[AssociationQ[Lookup[stats, "terminal_diagnostic", None]],
        terminalReplays[side] =
          stats["terminal_diagnostic", "replay"]]]]], {"lower", "upper"}];

If[Environment["DE2_REBUILD_TERMINAL_MATCH"] === "1",
  replayState = Lookup[restoredStateHandles, "upper", None];
  If[!AssociationQ[replayState],
    Print["NATIVE TERMINAL REBUILD state unavailable"],
    rebuildRequest = <|
      "schema" -> 2,
      "op" -> "transport.diagnostic_rebuild_terminal_match",
      "session" -> session,
      "transport_state" -> replayState["transport_state"]|>;
    scanPoints = Environment["DE2_TERMINAL_MATCH_SCAN_POINTS"];
    If[StringQ[scanPoints] &&
        StringLength[StringTrim[scanPoints]] > 0,
      rebuildRequest["receiving_local_points"] =
        StringTrim /@ StringSplit[scanPoints, ","]];
    scanPrecision = Environment["DE2_TERMINAL_MATCH_PRECISION_BITS"];
    If[StringQ[scanPrecision] &&
        StringMatchQ[StringTrim[scanPrecision],
          DigitCharacter ..],
      rebuildRequest["precision_bits"] =
        FromDigits[StringTrim[scanPrecision]]];
    rebuilt = DiffExp2`CppBackend`RunRequest[rebuildRequest];
    If[AssociationQ[rebuilt],
      Print["NATIVE TERMINAL REBUILD ",
        InputForm[If[ListQ[Lookup[rebuilt, "scan", None]],
          Map[
            Function[item,
              If[AssociationQ[Lookup[item, "match", None]],
                Join[KeyTake[item, {
                    "receiving_local", "incoming_local",
                    "physical", "status"}],
                  <|"match" -> KeyTake[item["match"], {
                      "matching_frame_identity",
                      "normal_frame_attempt",
                      "materialization_association",
                      "residual"}]|>],
                item]],
            rebuilt["scan"]],
          If[AssociationQ[Lookup[rebuilt, "rebuilt", None]],
          KeyTake[
            rebuilt["rebuilt"], {
              "basis_point_exact", "incoming_point_exact",
              "physical_match_point_exact", "epsilon",
              "normal_frame_attempt", "materialization_association",
              "residual"}],
          rebuilt]]]],
      Print["NATIVE TERMINAL REBUILD ", InputForm[rebuilt]]]]];

Quiet[DiffExp2`CppBackend`ClosePersistentSession[restored]];
Exit[0];
