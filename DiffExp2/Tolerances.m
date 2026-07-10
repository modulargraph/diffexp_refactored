(* DiffExp2/Tolerances.m — the single home for every numeric threshold.
   Spec: Docs/specs/Tolerances.md (binding); decisions: Docs/specs/DECISIONS-M0.md.
   Bottom module: depends on NOTHING. *)

BeginPackage["DiffExp2`Tolerances`"];

ChopDigits::usage = "ChopDigits[wp] gives the default chop-digit count Floor[wp/2].";
ChopFloor::usage = "ChopFloor[chopDigits] gives the absolute series noise floor 10^-chopDigits.";
MatchTol::usage = "MatchTol[matchDigits] gives the absolute matching/pivot zero test 10^-matchDigits.";
SnapTol::usage = "SnapTol[wp] gives the snap tolerance 10^(-Floor[wp/2]) for COMPUTED values only.";
InputSnapTol::usage = "InputSnapTol[v] gives the input-scaled snap tolerance 10^(-Floor[3*Min[Precision[v],Accuracy[v]]/4]) for an INEXACT external input v.";
RankTol::usage = "RankTol[wp] gives the relative rank/nullspace threshold 10^(-Floor[wp/4]).";
GeomGuardTol::usage = "GeomGuardTol[wp] gives the radius-vs-pole comparison guard 10^(-Floor[wp/2]).";
ChopReserve::usage = "ChopReserve[wp, chopDigits] gives the digit reserve wp - chopDigits.";
LaurentLeadTol::usage = "LaurentLeadTol[chopDigits] gives the RELATIVE leading-coefficient zero tolerance Max[10^(-Floor[chopDigits/2]), 10^-24]; the 10^-24 floor is a hard constant (DEC-2).";
ResidTol::usage = "ResidTol[wp] gives the relative ODE-residual spot-check threshold 10^(-Floor[wp/10]).";
NumericMagnitude::usage = "NumericMagnitude[c, digits] gives a stable Euclidean magnitude of the numeric scalar c without allowing an uncertain zero complex component to poison the resolved component.";
NumericMagnitudeBounds::usage = "NumericMagnitudeBounds[c, digits] gives conservative componentwise-uncertainty {lower, upper} bounds on the modulus of numeric scalar c.";
NumericallyZeroQ::usage = "NumericallyZeroQ[c, scale, tol, context, bandDecades, binaryFloor] is THE library-wide uncertainty-aware zero classification: True | False | loud ambiguity error; Laurent callers explicitly pass binaryFloor=True for the 10^-24 floor.";
DE2Error::usage = "DE2Error[id, payload] prints a one-line summary and throws Failure[\"DiffExp2\", payload] on tag \"DiffExp2Error\". The library-wide error primitive (DEC-1).";
InstallToleranceState::usage = "InstallToleranceState[assoc] atomically installs the nine-key tolerance record (called by Config.m only).";
Tol::usage = "Tol[name] is the validated read of the installed tolerance state. Unknown name or no installed state is a loud error, never a default.";
ToleranceStateInstalledQ::usage = "ToleranceStateInstalledQ[] gives True if a tolerance state is installed.";
EvalErrorSeriesDecrease::usage = "EvalErrorSeriesDecrease[couplingDepth] gives the error-probe order reduction Ceiling[0.7*couplingDepth] + 2.";
$SafetyDigits::usage = "$SafetyDigits = 2: digit-budget safety margin (old ISafetyDigits).";
$InputPrecisionFactor::usage = "$InputPrecisionFactor = 2: inputs are raised to this multiple of WorkingPrecision.";
$MaxExtraPrecisionValue::usage = "$MaxExtraPrecisionValue = 1000: value for $MaxExtraPrecision in evaluation blocks.";
$MinExpansionOrder::usage = "$MinExpansionOrder = 10: floor for the ExpansionOrder validator.";
$AmbiguityBandDecades::usage = "$AmbiguityBandDecades = 4: half-width in decades of the NumericallyZeroQ loud-error band (exempt at the 10^-24 floor).";
$NearSingularityGuardDecades::usage = "$NearSingularityGuardDecades = 6: outer edge of the loud near-singularity guard zone for inexact targets.";

Begin["`Private`"];

$SafetyDigits = 2;
$InputPrecisionFactor = 2;
$MaxExtraPrecisionValue = 1000;
$MinExpansionOrder = 10;
$AmbiguityBandDecades = 4;
$NearSingularityGuardDecades = 6;

DE2Error[id_String, payload_Association] := Module[
  (* "Message" is a reserved Failure property; store free text under "Detail" *)
  {p = Join[<|"ID" -> id|>, KeyMap[Replace[#, "Message" -> "Detail"] &, payload]]},
  Print["DiffExp2 error ", id, ": ",
    StringRiffle[KeyValueMap[ToString[#1] <> "=" <> ToString[#2, InputForm] &, p], "; "]];
  Throw[Failure["DiffExp2", p], "DiffExp2Error"]
];

badArg[fn_String, x___] := DE2Error["E6",
  <|"Module" -> "Tolerances", "Function" -> fn, "Arguments" -> {x},
    "Message" -> "argument must be a positive machine Integer"|>];

ChopDigits[wp_Integer?Positive] := Floor[wp/2];
ChopDigits[x___] := badArg["ChopDigits", x];
ChopFloor[cd_Integer?Positive] := 10^-cd;
ChopFloor[x___] := badArg["ChopFloor", x];
MatchTol[md_Integer?Positive] := 10^-md;
MatchTol[x___] := badArg["MatchTol", x];
SnapTol[wp_Integer?Positive] := 10^(-Floor[wp/2]);
SnapTol[x___] := badArg["SnapTol", x];
RankTol[wp_Integer?Positive] := 10^(-Floor[wp/4]);
RankTol[x___] := badArg["RankTol", x];
GeomGuardTol[wp_Integer?Positive] := 10^(-Floor[wp/2]);
GeomGuardTol[x___] := badArg["GeomGuardTol", x];
ChopReserve[wp_Integer?Positive, cd_Integer?Positive] := wp - cd;
ChopReserve[x___] := badArg["ChopReserve", x];
LaurentLeadTol[cd_Integer?Positive] := Max[10^(-Floor[cd/2]), 10^-24];
LaurentLeadTol[x___] := badArg["LaurentLeadTol", x];
ResidTol[wp_Integer?Positive] := 10^(-Floor[wp/10]);
ResidTol[x___] := badArg["ResidTol", x];
EvalErrorSeriesDecrease[cdep_Integer?Positive] := Ceiling[7*cdep/10] + 2;
EvalErrorSeriesDecrease[x___] := badArg["EvalErrorSeriesDecrease", x];

InputSnapTol[v_?InexactNumberQ] := 10^(-Floor[3*Min[Precision[v], Accuracy[v]]/4]);
InputSnapTol[x___] := DE2Error["E6",
  <|"Module" -> "Tolerances", "Function" -> "InputSnapTol", "Arguments" -> {x},
    "Message" -> "argument must be an inexact number; exact inputs use exact membership only"|>];

(* Stable TRUE complex modulus.  Arbitrary-precision Complex arithmetic can
   attach poor accuracy to a centered-zero component; Abs then lets that
   irrelevant uncertainty poison Sqrt[re^2+im^2] and can collapse a clearly
   nonzero other component to 0``negative-accuracy.  Inspect the stored
   component centers independently, then use a scaled Euclidean norm without
   dropping any nonzero center.  Scaling avoids overflow/underflow and, unlike the
   tempting infinity norm, preserves threshold semantics for values such as
   (8 + 8 I) 10^-25 at the 10^-24 Laurent floor.  Precision uncertainty is a
   separate proof obligation at residual call sites; this function computes
   the stable central magnitude only. *)
NumericMagnitude[c_?NumericQ, digits_Integer?Positive] := Module[
  {z = N[c, digits], parts, mags, centers, hiAt, hi, lo},
  parts = If[Head[z] === Complex, {Re[z], Im[z]}, {z}];
  mags = Abs /@ parts;
  (* Rationalize[...,0] is used only to inspect/order the stored central
     values.  Re-numericizing those centers prevents their original
     uncertainty metadata from contaminating this CENTRAL-value helper; the
     modulus itself remains a scaled floating computation, so this hot helper
     does not build and square giant exact rationals. *)
  centers = Rationalize[mags, 0];
  mags = N[centers, digits];
  Which[
    Max[centers] === 0, 0,
    Length[mags] === 1, First[mags],
    True,
      hiAt = First[Ordering[centers, -1]];
      hi = mags[[hiAt]];
      lo = mags[[3 - hiAt]];
      If[centers[[3 - hiAt]] === 0, hi,
        hi*Sqrt[1 + (lo/hi)^2]]]];
NumericMagnitude[c_, digits_] := DE2Error["E6",
  <|"Module" -> "Tolerances", "Function" -> "NumericMagnitude",
    "Arguments" -> {c, digits},
    "Message" -> "expected a numeric scalar and a positive machine Integer digit count"|>];

(* Exact stored-center magnitude and componentwise uncertainty enclosure.
   A centered-zero imaginary component contributes zero to the LOWER bound,
   but its uncertainty cannot erase a separately resolved real component.
   Arbitrary-precision powers of ten preserve the Accuracy-sized uncertainty
   without machine underflow above 308 digits. *)
centralMagnitudeExact[c_?NumericQ, digits_Integer?Positive] := Module[
  {q = If[FreeQ[c, _?InexactNumberQ], c,
      Rationalize[N[c, digits], 0]]},
  Sqrt[Re[q]^2 + Im[q]^2]];

NumericMagnitudeBounds[c_?NumericQ, digits_Integer?Positive] := Module[
  {z = N[c, digits], parts, centers, exactInput, uncertainty, lo, hi},
  exactInput = FreeQ[c, _?InexactNumberQ];
  If[exactInput, z = c];
  parts = If[Head[z] === Complex, {Re[z], Im[z]}, {z}];
  centers = Abs[If[exactInput, parts, Rationalize[parts, 0]]];
  uncertainty = If[exactInput, ConstantArray[0, Length[parts]],
    Map[Module[{acc = Accuracy[#]},
      If[acc === Infinity, 0,
        If[NumericQ[acc],
          N[10^(-SetPrecision[acc, 30]), 30], Infinity]]] &, parts]];
  lo = MapThread[Max[0, #1 - #2] &, {centers, uncertainty}];
  hi = centers + uncertainty;
  {Sqrt[Total[lo^2]], Sqrt[Total[hi^2]]}];
NumericMagnitudeBounds[c_, digits_] := DE2Error["E6",
  <|"Module" -> "Tolerances", "Function" -> "NumericMagnitudeBounds",
    "Arguments" -> {c, digits},
    "Message" -> "expected a numeric scalar and a positive machine Integer digit count"|>];

NumericallyZeroQ[c_, scale_, tol_, context_String,
    bandDecades_:$AmbiguityBandDecades, binaryFloor_:False] :=
  Which[
    FreeQ[c, _?InexactNumberQ] && TrueQ[PossibleZeroQ[c]], True,
    !NumericQ[c], False,
    !NumericQ[scale], False,
    !MemberQ[{True, False}, binaryFloor], DE2Error["E6",
      <|"Module" -> "Tolerances", "Function" -> "NumericallyZeroQ",
        "BinaryFloor" -> binaryFloor, "Context" -> context,
        "Message" -> "binaryFloor must be True or False"|>],
    True,
    Module[{digits = Tol["ChopDigits"], bounds, lo, hi, scaleCenter},
      bounds = NumericMagnitudeBounds[c, digits];
      {lo, hi} = bounds;
      scaleCenter = centralMagnitudeExact[scale, digits];
      If[TrueQ[scaleCenter === 0], Return[False, Module]];
      Which[
        TrueQ[binaryFloor], TrueQ[hi <= tol*scaleCenter],
        TrueQ[hi < tol*scaleCenter/10^bandDecades], True,
        TrueQ[lo > tol*scaleCenter*10^bandDecades], False,
        True, DE2Error["E5",
          <|"Module" -> "Tolerances", "Coefficient" -> N[c, 6], "Scale" -> N[scale, 6],
            "MagnitudeLowerBound" -> N[lo, 6],
            "MagnitudeUpperBound" -> N[hi, 6],
            "Tolerance" -> tol, "BandDecades" -> bandDecades,
            "BinaryFloor" -> binaryFloor, "Context" -> context,
            "Message" -> "cannot classify coefficient against scale within the ambiguity band"|>]
      ]
    ]
  ];

$tolKeys = {"WorkingPrecision", "ChopDigits", "MatchDigits", "ChopFloor",
  "MatchTol", "SnapTol", "RankTol", "LaurentLeadTol", "ResidTol"};
$intKeys = {"WorkingPrecision", "ChopDigits", "MatchDigits"};
$tolState = None;

InstallToleranceState[assoc_Association] := Module[
  {missing = Complement[$tolKeys, Keys[assoc]],
   extra = Complement[Keys[assoc], $tolKeys], badTypes},
  If[missing =!= {} || extra =!= {},
    DE2Error["E3", <|"Module" -> "Tolerances", "MissingKeys" -> missing,
      "ExtraKeys" -> extra, "Message" -> "tolerance record must contain exactly the nine schema keys"|>]];
  badTypes = Select[$tolKeys, !If[MemberQ[$intKeys, #],
    IntegerQ[assoc[#]] && assoc[#] > 0,
    Head[assoc[#]] === Rational && 0 < assoc[#] < 1] &];
  If[badTypes =!= {},
    DE2Error["E3", <|"Module" -> "Tolerances", "BadTypeKeys" -> badTypes,
      "Values" -> (assoc /@ badTypes),
      "Message" -> "Integer expected for WorkingPrecision/ChopDigits/MatchDigits; exact positive Rational < 1 for tolerance values"|>]];
  If[assoc["ChopDigits"] >= assoc["WorkingPrecision"],
    DE2Error["E4", <|"Module" -> "Tolerances", "ChopDigits" -> assoc["ChopDigits"],
      "WorkingPrecision" -> assoc["WorkingPrecision"],
      "Message" -> "The value of ChopPrecision should be smaller than the value of WorkingPrecision"|>]];
  If[assoc["MatchDigits"] > assoc["ChopDigits"],
    DE2Error["E3", <|"Module" -> "Tolerances", "MatchDigits" -> assoc["MatchDigits"],
      "ChopDigits" -> assoc["ChopDigits"],
      "Message" -> "MatchDigits must not exceed ChopDigits (the auto-sync invariant)"|>]];
  $tolState = assoc;
  Null
];
InstallToleranceState[x___] := DE2Error["E3",
  <|"Module" -> "Tolerances", "Arguments" -> {x},
    "Message" -> "InstallToleranceState expects one Association"|>];

Tol[name_String] := Which[
  $tolState === None,
  DE2Error["E1", <|"Module" -> "Tolerances", "Key" -> name,
    "Message" -> "no tolerance state installed; load a configuration first"|>],
  !MemberQ[$tolKeys, name],
  DE2Error["E2", <|"Module" -> "Tolerances", "Key" -> name, "ValidNames" -> $tolKeys,
    "Message" -> "unknown tolerance name"|>],
  True, $tolState[name]];
Tol[x___] := DE2Error["E2",
  <|"Module" -> "Tolerances", "Arguments" -> {x}, "Message" -> "Tol expects one String name"|>];

ToleranceStateInstalledQ[] := $tolState =!= None;

End[];
EndPackage[];
