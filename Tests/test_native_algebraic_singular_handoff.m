(* An algebraic chart anywhere in an arm activates certified physical
   geometry.  Its terminal regular-to-singular handoff must be regenerated
   after scale bridging instead of retaining SegmentLine's pre-bridge point. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 80, "ExpansionOrder" -> 10,
  "EpsilonOrder" -> 0, "RecurrenceBackend" -> "Cpp",
  "Variables" -> {}}]];

roc = DiffExp2`Config`CFG["RadiusOfConvergence"];
algebraicCenter = (3 - Sqrt[3])/6;
algebraicScale = Sqrt[3]/2;
algebraicRadius = RootReduce[algebraicScale/2];
staleTerminalPoint = 99/100;

anchor = <|"Center" -> 0, "Scale" -> 1, "Radius" -> 1,
  "LocalRadius" -> 1, "MatchRadius" -> roc,
  "Singular" -> False, "Prescriptions" -> {}|>;
algebraic = <|"Center" -> algebraicCenter, "Scale" -> algebraicScale,
  "Radius" -> algebraicRadius, "LocalRadius" -> 1/2,
  "MatchRadius" -> RootReduce[algebraicScale roc],
  "Singular" -> False, "Prescriptions" -> {},
  "IncomingMatchPoint" -> algebraicCenter/2|>;
clearance = <|"Center" -> 3/4, "Scale" -> 1/4,
  "Radius" -> 1/4, "LocalRadius" -> 1,
  "MatchRadius" -> roc/4, "Singular" -> False,
  "Prescriptions" -> {}, "IncomingMatchPoint" -> 2/3|>;
terminal = <|"Center" -> 1, "Scale" -> 1, "Radius" -> 2,
  "LocalRadius" -> 2, "MatchRadius" -> roc,
  "Singular" -> True, "Prescriptions" -> {},
  "IncomingMatchPoint" -> staleTerminalPoint|>;
plan = <|"From" -> 0, "To" -> 1, "Direction" -> 1,
  "Charts" -> {anchor, algebraic, clearance, terminal}|>;

bridged = catchDE2[
  DiffExp2`NativeTransport`Private`bridgeNativeRegularPlanScales[plan]];
expected = If[AssociationQ[bridged],
  DiffExp2`Transport`Private`singularMatchPoint[
    bridged["Charts"][[3, "Center"]], bridged["Charts"][[3, "Radius"]],
    bridged["Charts"][[3, "MatchRadius"]],
    bridged["Charts"][[4, "Center"]], bridged["Charts"][[4, "Radius"]],
    bridged["Charts"][[4, "MatchRadius"]], 1,
    DiffExp2`Config`CFG["DivisionOrder"]],
  $Failed];
actual = If[AssociationQ[bridged],
  Lookup[bridged["Charts"][[4]],
    "NativeCertifiedIncomingMatchPoint", None], None];
certified = If[AssociationQ[bridged], catchDE2[
  DiffExp2`NativeTransport`Private`nativeCertifiedArmGeometry[bridged]],
  bridged];

exactSameQ = AssociationQ[bridged] &&
  TrueQ[PossibleZeroQ[RootReduce[actual - expected]]];
ok = AssociationQ[bridged] && AssociationQ[certified] &&
  exactSameQ &&
  !TrueQ[PossibleZeroQ[RootReduce[actual - staleTerminalPoint]]] &&
  TrueQ[0 < actual < 1] &&
  TrueQ[Abs[actual - bridged["Charts"][[3, "Center"]]] <
    bridged["Charts"][[3, "Radius"]]] &&
  TrueQ[Abs[actual - 1] <=
    DiffExp2`Transport`Private`matchOffset[
      bridged["Charts"][[4, "MatchRadius"]],
      bridged["Charts"][[4, "Radius"]],
      DiffExp2`Config`CFG["DivisionOrder"]]] &&
  certified["matches"][[3, "physical", "exact"]] ===
    ToString[RootReduce[actual], InputForm];

Print[If[TrueQ[ok], "PASS", "FAIL"],
  ": algebraic certified geometry rebalances terminal singular handoff"];
If[!TrueQ[ok], Print[InputForm[{
  "Bridged" -> bridged, "Expected" -> expected, "Actual" -> actual,
  "Certified" -> certified}]]];
Exit[If[TrueQ[ok], 0, 1]];
