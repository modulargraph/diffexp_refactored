(* DiffExp2/Config.m — the single validated option store.
   Spec: Docs/specs/Config.md (binding); decisions: Docs/specs/DECISIONS-M0.md.
   May call ONLY DiffExp2`Tolerances`. *)

BeginPackage["DiffExp2`Config`", {"DiffExp2`Tolerances`"}];

LoadConfiguration::usage = "LoadConfiguration[rules] resets every key to its schema default, then applies rules. Returns CurrentConfiguration[].";
UpdateConfiguration::usage = "UpdateConfiguration[rules] validates and merges rules all-or-nothing; recomputes and installs the tolerance state.";
CurrentConfiguration::usage = "CurrentConfiguration[] gives all known keys and current (resolved) values, String-keyed.";
CFG::usage = "CFG[key] is THE validated configuration read. key may be a String or Symbol (canonicalized by SymbolName). Unknown key or unloaded configuration is a loud error, never Missing.";
ConfigSchema::usage = "ConfigSchema[] gives a read-only copy of the configuration schema.";
CanonicalKey::usage = "CanonicalKey[key] gives the canonical String key for a String or Symbol key.";
PinnedVariable::usage = "PinnedVariable[sym] returns the symbol with the same name in the pinned context Global`.";
ConfiguredQ::usage = "ConfiguredQ[] gives True after the first successful LoadConfiguration/UpdateConfiguration.";
PrintInfo::usage = "PrintInfo[level, args...] prints args when CFG[\"Verbosity\"] >= level. The only verbosity-gated print helper in DiffExp2.";
PrintWarning::usage = "PrintWarning[args...] prints a warning; never gated below visibility of Verbosity level 1.";
EpsSymbols::usage = "EpsSymbols[] gives the accepted regulator symbols {Global`eps, Global`\\[Epsilon]}.";
CanonicalEps::usage = "CanonicalEps[] gives the canonical regulator symbol Global`eps.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id, Join[<|"Module" -> "Config"|>, payload]];

PinnedVariable[sym_Symbol] := Symbol["Global`" <> SymbolName[sym]];
CanonicalKey[key_String] := key;
CanonicalKey[key_Symbol] := SymbolName[key];
EpsSymbols[] := {Global`eps, Global`\[Epsilon]};
CanonicalEps[] := Global`eps;

pinAll[expr_] := expr /. s_Symbol /; Context[s] =!= "System`" :> PinnedVariable[s];

(* DeltaPrescriptions: parse, validate, sign-aware canonicalization (spec 3.3, DEC-16). *)
parseOnePrescription[{poly_, sign_}] := {pinAll[poly], sign};
parseOnePrescription[expr_] := Module[{d = Global`\[Delta], c, poly},
  If[!PolynomialQ[expr, d] || Exponent[expr, d] =!= 1,
    err["E7", <|"Element" -> expr, "Message" -> "delta must appear exactly linearly with coefficient ±I"|>]];
  c = Coefficient[expr, d, 1];
  poly = expr /. d -> 0;
  If[!FreeQ[c, d] || PossibleZeroQ[c],
    err["E7", <|"Element" -> expr, "Message" -> "delta must appear exactly linearly with coefficient ±I"|>]];
  {pinAll[poly], Simplify[c/I]}];
canonOrient[{poly_, sign_}] := Module[{vars, lead, nc},
  vars = Sort[DeleteDuplicates[Cases[poly, s_Symbol /; Context[s] === "Global`", {0, Infinity}]]];
  lead = First[MonomialList[poly, vars]];
  nc = lead /. s_Symbol /; Context[s] === "Global`" :> 1;
  If[TrueQ[nc < 0], {Expand[-poly], -sign}, {Expand[poly], sign}]];
validatePrescription[{poly_, sign_}] := Module[{fl, vars},
  If[!MemberQ[{1, -1}, sign],
    err["E7", <|"Element" -> {poly, sign}, "Message" -> "prescription sign must be +1 or -1"|>]];
  If[!FreeQ[poly, Global`\[Delta]],
    err["E7", <|"Element" -> {poly, sign},
      "Message" -> "reserved delta may appear only in the poly ± I delta expression form"|>]];
  vars = DeleteDuplicates[Cases[poly,
    s_Symbol /; Context[s] === "Global`", {0, Infinity}]];
  If[TrueQ[PossibleZeroQ[poly]] || vars === {} || !PolynomialQ[poly, vars],
    err["E7", <|"Element" -> {poly, sign},
      "Message" -> "prescription factor must be a nonzero polynomial in at least one pinned variable"|>]];
  fl = Select[FactorList[poly], !NumericQ[First[#]] &];
  If[Length[fl] > 1 || (Length[fl] === 1 && fl[[1, 2]] > 1),
    err["E6", <|"Polynomial" -> poly, "Message" -> "Physical singularities should be irreducible polynomials!"|>]];
  {poly, sign}];
parsePrescriptions[list_List] := Module[{norm, grouped},
  norm = canonOrient[validatePrescription[parseOnePrescription[#]]] & /@ list;
  grouped = GatherBy[norm, First];
  Map[(If[Length[DeleteDuplicates[#[[All, 2]]]] > 1,
    err["E7", <|"Elements" -> #, "Message" -> "same canonical factor with OPPOSITE i-delta sides (sign-aware dedup, DEC-16)"|>]];
    First[#]) &, grouped]];
parsePrescriptions[x_] := err["E3", <|"Key" -> "DeltaPrescriptions", "Value" -> x, "Expected" -> "a List"|>];

symbolListQ[l_] := ListQ[l] && AllTrue[l, MatchQ[#, _Symbol] &];
exactPosQ[v_] := (IntegerQ[v] || Head[v] === Rational) && v > 0;

$schema = <|
  "WorkingPrecision" -> <|"Type" -> (IntegerQ[#] && # >= 20 &), "Default" -> 500, "Normalize" -> None|>,
  "ChopPrecision" -> <|"Type" -> (# === Automatic || (IntegerQ[#] && # > 0) &), "Default" -> Automatic, "Normalize" -> None|>,
  "LinearSolveChopPrecision" -> <|"Type" -> (# === Automatic || (IntegerQ[#] && # > 0) &), "Default" -> Automatic, "Normalize" -> None|>,
  "AccuracyGoal" -> <|"Type" -> (# === "?" || (IntegerQ[#] && # > 0) &), "Default" -> "?", "Normalize" -> None|>,
  "AccuracyGoalValidate" -> <|"Type" -> BooleanQ, "Default" -> False,
    "Normalize" -> (Replace[#, {"Before" -> True, "After" -> True, None -> False}] &)|>,
  "DeltaPrescriptions" -> <|"Type" -> ListQ, "Default" -> {}, "Normalize" -> parsePrescriptions|>,
  "DivisionOrder" -> <|"Type" -> (IntegerQ[#] && # >= 2 &), "Default" -> 3, "Normalize" -> None|>,
  (* Compatibility key retained for old configurations and checkpoint
     metadata.  The active classic GetCPL/GetCPR planner couples placement
     and matching through DivisionOrder, so StepDivisionOrder no longer
     changes SegmentLine geometry. *)
  "StepDivisionOrder" -> <|"Type" -> (# === Automatic ||
      ((IntegerQ[#] || Head[#] === Rational) && 1 <= # <= 16) &),
    "Default" -> 3, "Normalize" -> None|>,
  "EpsilonOrder" -> <|"Type" -> (IntegerQ[#] && # >= 0 &), "Default" -> 4, "Normalize" -> None|>,
  "EstimateError" -> <|"Type" -> (MemberQ[{False, True, "Fast"}, #] &), "Default" -> "Fast", "Normalize" -> None|>,
  "ExpansionOrder" -> <|"Type" -> (IntegerQ[#] && # >= DiffExp2`Tolerances`$MinExpansionOrder &), "Default" -> 50, "Normalize" -> None|>,
  "RecurrenceBackend" -> <|"Type" -> (MemberQ[{"Wolfram", "Cpp"}, #] &),
    "Default" -> "Wolfram", "Normalize" -> None|>,
  "LineParameter" -> <|"Type" -> (MatchQ[#, _Symbol] &), "Default" -> Global`x, "Normalize" -> PinnedVariable|>,
  "MatrixDirectory" -> <|"Type" -> StringQ, "Default" -> "", "Normalize" -> None|>,
  "RadiusOfConvergence" -> <|"Type" -> exactPosQ, "Default" -> 1, "Normalize" -> None|>,
  "RationalizationTolerance" -> <|"Type" -> (# === Automatic || (Head[#] === Rational && 0 < # < 1) &), "Default" -> Automatic, "Normalize" -> None|>,
  "SegmentationStrategy" -> <|"Type" -> (# === "Predivision" &), "Default" -> "Predivision", "Normalize" -> None|>,
  "UsePade" -> <|"Type" -> BooleanQ, "Default" -> False, "Normalize" -> None|>,
  "Variables" -> <|"Type" -> symbolListQ, "Default" -> {}, "Normalize" -> (Map[PinnedVariable, #] &)|>,
  "Verbosity" -> <|"Type" -> (IntegerQ[#] && # >= 0 &), "Default" -> 1, "Normalize" -> None|>,
  "VerbosityDebug" -> <|"Type" -> (IntegerQ[#] && # >= 0 &), "Default" -> 0, "Normalize" -> None|>,
  "SaveExpansionsCompress" -> <|"Type" -> BooleanQ, "Default" -> False, "Normalize" -> None|>,
  "SaveExpansionsCompressDirectory" -> <|"Type" -> (StringQ[#] || # === None &), "Default" -> None, "Normalize" -> None|>,
  "SaveExpansionsOrder" -> <|"Type" -> (IntegerQ[#] || # === None &), "Default" -> None, "Normalize" -> None|>,
  "AbortOnAnalyticContinuationFail" -> <|"Type" -> BooleanQ, "Default" -> True, "Normalize" -> None|>
|>;

$droppedReasons = <|
  "IntegrationStrategy" -> "the strategy stack is deleted (RewritePlan I2); DiffExp2 has ONE solver",
  "UseRationalRecurrence" -> "the denominator-cleared recursion is the only path in Solve.m",
  "InvWronskSolver" -> "the Wronskian machinery is deleted",
  "HomogeneousSolve" -> "symbolic-eps Frobenius has one homogeneous path",
  "KeepMatrixExpansions" -> "cache policy is internal to DiffExp2",
  "Parallel" -> "no consumer existed in the old core",
  "IgnoreIndicialCheck" -> "the I1 indicial contract violation is a hard error and must not be suppressible",
  "CrosscheckLevel" -> "DiffExp2 invariants are always-on and not configurable",
  "CrosscheckFlags" -> "DiffExp2 invariants are always-on and not configurable",
  "LogFile" -> "no in-repo consumer; session logging is the shell's job",
  "UseMobius" -> "Mobius chart maps are dropped (DEC-18): every DiffExp2 chart is affine; RoC rescaling survives as RadiusOfConvergence"
|>;

$store = None;
$userDeltaPrescriptions = {};  (* v1: equals the effective list; kept for the v1.1 sqrt auto-prescription seam *)

ConfiguredQ[] := $store =!= None;
ConfigSchema[] := $schema;

resolveRead[name_String] := Module[{raw = $store[name]},
  Which[
    name === "ChopPrecision" && raw === Automatic, DiffExp2`Tolerances`ChopDigits[$store["WorkingPrecision"]],
    name === "LinearSolveChopPrecision" && raw === Automatic, resolveRead["ChopPrecision"],
    name === "RationalizationTolerance" && raw === Automatic, DiffExp2`Tolerances`SnapTol[$store["WorkingPrecision"]],
    True, raw]];

CFG[key : (_String | _Symbol)] := Module[{name = CanonicalKey[key]},
  If[$store === None,
    err["E11", <|"Key" -> name, "Message" -> "no configuration loaded"|>]];
  If[KeyExistsQ[$droppedReasons, name],
    err["E9", <|"Key" -> name, "Message" -> name <> " was removed in DiffExp2: " <> $droppedReasons[name]|>]];
  If[!KeyExistsQ[$schema, name],
    err["E1", <|"Key" -> name, "ValidKeys" -> Keys[$schema],
      "Message" -> "unknown configuration key" <>
        If[MatchQ[key, _Symbol], " (symbol context: " <> Context[key] <>
          " — a non-Global context here usually means a package-path symbol-resolution trap)", ""]|>]];
  resolveRead[name]];
CFG[x___] := err["E1", <|"Arguments" -> {x}, "Message" -> "CFG expects one String or Symbol key"|>];

CurrentConfiguration[] := If[$store === None,
  err["E11", <|"Message" -> "no configuration loaded"|>],
  AssociationMap[resolveRead, Keys[$schema]]];

doUpdate[base_Association, rules_] := Module[
  {assoc, names, candidate, wp, cd, md, snap, lpName, varNames},
  assoc = Association[Flatten[{rules}]];
  names = CanonicalKey /@ Keys[assoc];
  assoc = AssociationThread[names, Values[assoc]];
  Do[
    If[KeyExistsQ[$droppedReasons, n],
      err["E9", <|"Key" -> n, "Message" -> n <> " was removed in DiffExp2: " <> $droppedReasons[n] <> ". No keys were updated."|>]];
    If[!KeyExistsQ[$schema, n],
      err["E2", <|"Key" -> n, "ValidKeys" -> Keys[$schema], "Message" -> "unknown configuration key; no keys were updated"|>]],
    {n, names}];
  (* normalize, then validate types on normalized values *)
  assoc = Association[KeyValueMap[
    Function[{k, v}, k -> If[$schema[k, "Normalize"] === None, v, $schema[k, "Normalize"][v]]], assoc]];
  Do[
    If[n === "SegmentationStrategy" && assoc[n] === "Dynamic",
      err["E10", <|"Key" -> n, "Message" -> "Dynamic segmentation is not ported in DiffExp2 v1 (RewritePlan 3.2); use Predivision"|>]];
    If[!TrueQ[$schema[n, "Type"][assoc[n]]],
      err["E3", <|"Key" -> n, "Value" -> assoc[n], "Message" -> "invalid value (after normalization) for " <> n|>]],
    {n, names}];
  candidate = Join[base, assoc];
  (* LinearSolveChopPrecision auto-sync: explicit ChopPrecision without explicit LSCP in the SAME update resets LSCP *)
  If[MemberQ[names, "ChopPrecision"] && !MemberQ[names, "LinearSolveChopPrecision"],
    candidate["LinearSolveChopPrecision"] = Automatic];
  (* cross-field invariants on resolved values *)
  wp = candidate["WorkingPrecision"];
  cd = candidate["ChopPrecision"] /. Automatic -> DiffExp2`Tolerances`ChopDigits[wp];
  md = candidate["LinearSolveChopPrecision"] /. Automatic -> cd;
  If[cd >= wp,
    err["E4", <|"ChopPrecision" -> cd, "WorkingPrecision" -> wp,
      "Message" -> "The value of ChopPrecision should be smaller than the value of WorkingPrecision."|>]];
  If[md > cd,
    err["E4", <|"LinearSolveChopPrecision" -> md, "ChopPrecision" -> cd,
      "Message" -> "LinearSolveChopPrecision must not exceed ChopPrecision"|>]];
  lpName = SymbolName[candidate["LineParameter"]];
  varNames = SymbolName /@ candidate["Variables"];
  If[MemberQ[varNames, lpName],
    err["E5", <|"LineParameter" -> lpName, "Variables" -> varNames,
      "Message" -> "The symbol for the line parameter can't be equal to one of the kinematic invariants or masses."|>]];
  If[TrueQ[candidate["AccuracyGoalValidate"]] && !IntegerQ[candidate["AccuracyGoal"]],
    err["E8", <|"AccuracyGoal" -> candidate["AccuracyGoal"],
      "Message" -> "AccuracyGoalValidate -> True requires a numeric AccuracyGoal"|>]];
  snap = candidate["RationalizationTolerance"] /. Automatic -> DiffExp2`Tolerances`SnapTol[wp];
  (* commit: install tolerances, then the store (validate-then-commit, F-h) *)
  DiffExp2`Tolerances`InstallToleranceState[<|
    "WorkingPrecision" -> wp, "ChopDigits" -> cd, "MatchDigits" -> md,
    "ChopFloor" -> DiffExp2`Tolerances`ChopFloor[cd],
    "MatchTol" -> DiffExp2`Tolerances`MatchTol[md],
    "SnapTol" -> snap,
    "RankTol" -> DiffExp2`Tolerances`RankTol[wp],
    "LaurentLeadTol" -> DiffExp2`Tolerances`LaurentLeadTol[cd],
    "ResidTol" -> DiffExp2`Tolerances`ResidTol[wp]|>];
  $store = candidate;
  If[MemberQ[names, "DeltaPrescriptions"],
    $userDeltaPrescriptions = candidate["DeltaPrescriptions"]];
  CurrentConfiguration[]];

UpdateConfiguration[rules : ({___Rule} | _Association)] :=
  doUpdate[If[$store === None, AssociationMap[$schema[#, "Default"] &, Keys[$schema]], $store], Normal[rules]];
UpdateConfiguration[rules__Rule] := UpdateConfiguration[{rules}];
UpdateConfiguration[x___] := err["E2", <|"Arguments" -> {x}, "Message" -> "UpdateConfiguration expects rules or an Association"|>];

LoadConfiguration[rules : ({___Rule} | _Association)] :=
  doUpdate[AssociationMap[$schema[#, "Default"] &, Keys[$schema]], Normal[rules]];
LoadConfiguration[rules__Rule] := LoadConfiguration[{rules}];
LoadConfiguration[] := LoadConfiguration[{}];
LoadConfiguration[x___] := err["E2", <|"Arguments" -> {x}, "Message" -> "LoadConfiguration expects rules or an Association"|>];

PrintInfo[level_Integer, args__] :=
  If[TrueQ[ConfiguredQ[]] && resolveRead["Verbosity"] >= level, Print[args]];
PrintWarning[args__] :=
  If[!TrueQ[ConfiguredQ[]] || resolveRead["Verbosity"] >= 1, Print["DiffExp2 warning: ", args]];

End[];
EndPackage[];
