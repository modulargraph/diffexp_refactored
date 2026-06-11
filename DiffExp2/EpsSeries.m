(* DiffExp2/EpsSeries.m — truncated eps-Laurent arrays with honest windows.
   Spec: Docs/specs/EpsSeries.md (binding); decisions: Docs/specs/DECISIONS-M0.md.
   May call only Tolerances (predicates, DE2Error) and Config (pinned eps symbols). *)

BeginPackage["DiffExp2`EpsSeries`", {"DiffExp2`Tolerances`", "DiffExp2`Config`"}];

ESNew::usage = "ESNew[kmin, coeffs] builds the EpsSeries Σ coeffs[[k-kmin+1]] eps^k with CompleteMax = kmin + Length[coeffs] - 1.";
ESZero::usage = "ESZero[kmax] gives the canonical zero series known complete through kmax.";
ESQ::usage = "ESQ[expr] gives True if expr is a structurally valid EpsSeries.";
ESWindow::usage = "ESWindow[s] gives <|\"Min\" -> kmin, \"CompleteMax\" -> kmax|> verbatim in the RewritePlan 3.1 shape.";
ESMinPower::usage = "ESMinPower[s] gives the window Min.";
ESCompleteMax::usage = "ESCompleteMax[s] gives the window CompleteMax.";
ESCoefficient::usage = "ESCoefficient[s, k] gives the eps^k coefficient: stored in window, exact 0 below Min (the zero guarantee), LOUD ERROR above CompleteMax.";
ESCoefficientList::usage = "ESCoefficientList[s, k1, k2] gives coefficients of eps^k1..eps^k2; loud on incomplete range or on silently dropping nonzero content below k1.";
ESLeading::usage = "ESLeading[s] gives {k, coefficient} of the lowest non-negligible order, or None.";
ESAdd::usage = "ESAdd[a, b, ...] gives the sum with window Min[Min]/Min[CompleteMax] — never a zero-padded union.";
ESScale::usage = "ESScale[c, s] multiplies every coefficient by the eps-free scalar c; window unchanged.";
ESShift::usage = "ESShift[s, j] multiplies by eps^j exactly: both window edges move by j.";
ESTimes::usage = "ESTimes[a, b] gives the Cauchy product with the computed honest window.";
ESInvert::usage = "ESInvert[d] gives 1/d with window [-L, dCM - 2L], L the leading index; an eps^0-vanishing denominator is a window SHIFT, not an error.";
ESDivide::usage = "ESDivide[a, b] gives a/b = ESTimes[a, ESInvert[b]].";
ESTruncate::usage = "ESTruncate[s, newCM] lowers CompleteMax; extending or emptying the window is a loud error.";
ESTrim::usage = "ESTrim[s] advances Min past negligible leading coefficients; CompleteMax never changes; all-negligible gives ESZero[CompleteMax].";
ESCoeffZeroQ::usage = "ESCoeffZeroQ[c, scale, tol] is the EpsSeries wrapper over Tolerances`NumericallyZeroQ (tol Automatic -> LaurentLeadTol).";
ESSameQ::usage = "ESSameQ[a, b] gives equality at matchTol on the shared complete range; window agreement is not required.";
ESFromExpression::usage = "ESFromExpression[expr, epsSym, kmax] gives the exact Laurent expansion complete through kmax; loud on failure, fractional powers, or residual eps. API boundary only.";
ESToExpression::usage = "ESToExpression[s, epsSym] gives the plain expression Σ c_k epsSym^k (window metadata is lost). API boundary only.";
ESMap::usage = "ESMap[f, s] applies f to every stored coefficient; window unchanged; result revalidated.";
$ESErrorContext::usage = "$ESErrorContext is the caller-set context string (chart/sector/order) embedded in every EpsSeries error.";

Begin["`Private`"];

$ESErrorContext = "(no context)";

esError[id_String, payload_Association] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "EpsSeries", "Context" -> $ESErrorContext|>, payload]];

(* ---- structure ---- *)

ESQ[s_] := AssociationQ[s] && Sort[Keys[s]] === {"Coeffs", "EpsWindow"} &&
  AssociationQ[s["EpsWindow"]] &&
  Sort[Keys[s["EpsWindow"]]] === {"CompleteMax", "Min"} &&
  IntegerQ[s["EpsWindow", "Min"]] && IntegerQ[s["EpsWindow", "CompleteMax"]] &&
  ListQ[s["Coeffs"]] &&
  Length[s["Coeffs"]] === s["EpsWindow", "CompleteMax"] - s["EpsWindow", "Min"] + 1 &&
  Length[s["Coeffs"]] >= 1;

(* single choke point (I-1, I-2); no validation of coefficient content here *)
mkSeries[kmin_Integer, kmax_Integer, coeffs_List] := (
  If[kmin > kmax || Length[coeffs] =!= kmax - kmin + 1,
    esError["ERR-WINDOW-EMPTY", <|"Min" -> kmin, "CompleteMax" -> kmax,
      "CoeffLength" -> Length[coeffs], "Detail" -> "internal window assertion failed"|>]];
  <|"EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>, "Coeffs" -> coeffs|>);

$badHeads = {$Failed, Indeterminate, ComplexInfinity, Overflow[]};
badCoeffQ[c_] := MemberQ[$badHeads, c] || MatchQ[c, DirectedInfinity[___]];
epsFreeQ[c_] := AllTrue[DiffExp2`Config`EpsSymbols[], FreeQ[c, #] &];

validateCoeffs[fn_String, coeffs_List] := Module[{bad},
  bad = SelectFirst[Transpose[{Range[Length[coeffs]], coeffs}],
    badCoeffQ[#[[2]]] || !epsFreeQ[#[[2]]] &, None];
  If[bad =!= None,
    esError["ERR-BAD-CONSTRUCT", <|"Function" -> fn, "Index" -> bad[[1]],
      "Coefficient" -> StringTake[ToString[bad[[2]], InputForm], UpTo[200]],
      "Detail" -> "coefficient is invalid or eps-dependent"|>]];
  coeffs];

ESNew[kmin_Integer, coeffs_List /; Length[coeffs] >= 1] :=
  mkSeries[kmin, kmin + Length[coeffs] - 1, validateCoeffs["ESNew", coeffs]];
ESNew[x___] := esError["ERR-BAD-CONSTRUCT",
  <|"Function" -> "ESNew", "Arguments" -> StringTake[ToString[{x}, InputForm], UpTo[200]],
    "Detail" -> "expected ESNew[kmin_Integer, coeffs_List (non-empty)]"|>];

ESZero[kmax_Integer] := mkSeries[kmax, kmax, {0}];
ESZero[x___] := esError["ERR-BAD-CONSTRUCT",
  <|"Function" -> "ESZero", "Arguments" -> {x}, "Detail" -> "expected ESZero[kmax_Integer]"|>];

(* ---- accessors ---- *)

ESWindow[s_?ESQ] := s["EpsWindow"];
ESMinPower[s_?ESQ] := s["EpsWindow", "Min"];
ESCompleteMax[s_?ESQ] := s["EpsWindow", "CompleteMax"];

ESCoefficient[s_?ESQ, k_Integer] := Which[
  k < ESMinPower[s], 0,
  k <= ESCompleteMax[s], s["Coeffs"][[k - ESMinPower[s] + 1]],
  True, esError["ERR-WINDOW-READ", <|"Requested" -> k,
    "Min" -> ESMinPower[s], "CompleteMax" -> ESCompleteMax[s],
    "Detail" -> "coefficient requested above the complete window"|>]];

(* series scale: max numeric coefficient magnitude (symbolic excluded) *)
numMag[c_] := If[NumericQ[c], Abs[N[c, 10]], Nothing];
seriesScale[s_?ESQ] := Max[0, Sequence @@ (numMag /@ s["Coeffs"])];

ESCoeffZeroQ[c_, scale_, tol_:Automatic] :=
  DiffExp2`Tolerances`NumericallyZeroQ[c, scale,
    tol /. Automatic -> DiffExp2`Tolerances`Tol["LaurentLeadTol"], $ESErrorContext];

ESLeading[s_?ESQ] := Module[{scale = seriesScale[s], kmin = ESMinPower[s]},
  Do[If[!ESCoeffZeroQ[s["Coeffs"][[i]], scale], Return[{kmin + i - 1, s["Coeffs"][[i]]}, Module]],
    {i, Length[s["Coeffs"]]}];
  None];

ESCoefficientList[s_?ESQ, k1_Integer, k2_Integer] := Module[
  {kmin = ESMinPower[s], kmax = ESCompleteMax[s], scale},
  If[k1 > k2,
    (* the ONE sanctioned empty slice: the polar part of a pole-free series *)
    If[k2 === -1 && k1 === kmin && kmin >= 0, Return[{}, Module]];
    esError["ERR-RANGE", <|"k1" -> k1, "k2" -> k2, "Min" -> kmin, "CompleteMax" -> kmax,
      "Detail" -> "empty coefficient range requested"|>]];
  If[k2 > kmax,
    esError["ERR-RANGE", <|"k1" -> k1, "k2" -> k2, "Min" -> kmin, "CompleteMax" -> kmax,
      "Detail" -> "range extends above the complete window"|>]];
  If[k1 > kmin,
    scale = seriesScale[s];
    Do[If[!ESCoeffZeroQ[ESCoefficient[s, k], scale],
      esError["ERR-DROP-BELOW", <|"DroppedOrder" -> k, "k1" -> k1,
        "Magnitude" -> If[NumericQ[ESCoefficient[s, k]], N[Abs[ESCoefficient[s, k]], 6], "Symbolic"],
        "Scale" -> N[scale, 6], "Min" -> kmin, "CompleteMax" -> kmax,
        "Detail" -> "slicing would silently drop non-negligible content below k1"|>]],
      {k, kmin, k1 - 1}]];
  Table[ESCoefficient[s, k], {k, k1, k2}]];

(* ---- arithmetic ---- *)

norm[c_] := If[NumericQ[c], c, Together[Expand[c]]];

ESAdd[a_?ESQ, b_?ESQ] := Module[
  {kmin = Min[ESMinPower[a], ESMinPower[b]],
   kmax = Min[ESCompleteMax[a], ESCompleteMax[b]]},
  mkSeries[kmin, kmax,
    Table[norm[ESCoefficient[a, k] + ESCoefficient[b, k]], {k, kmin, kmax}]]];
ESAdd[a_?ESQ, b_?ESQ, rest__] := ESAdd[ESAdd[a, b], rest];
ESAdd[a_?ESQ] := a;

ESScale[c_, s_?ESQ] := (
  validateCoeffs["ESScale", {c}];
  mkSeries[ESMinPower[s], ESCompleteMax[s], norm[c*#] & /@ s["Coeffs"]]);

ESShift[s_?ESQ, j_Integer] :=
  mkSeries[ESMinPower[s] + j, ESCompleteMax[s] + j, s["Coeffs"]];

ESTimes[a_?ESQ, b_?ESQ] := Module[
  {aMin = ESMinPower[a], bMin = ESMinPower[b],
   kmin, kmax},
  kmin = aMin + bMin;
  kmax = Min[ESCompleteMax[a] + bMin, ESCompleteMax[b] + aMin];
  mkSeries[kmin, kmax,
    Table[norm[Sum[ESCoefficient[a, i]*ESCoefficient[b, k - i],
      {i, Max[aMin, k - ESCompleteMax[b]], Min[ESCompleteMax[a], k - bMin]}]],
      {k, kmin, kmax}]]];

ESInvert[d_?ESQ] := Module[{lead, L, dL, relCM, e},
  lead = ESLeading[d];
  If[lead === None,
    esError["ERR-DIV-ZERO", <|"Min" -> ESMinPower[d], "CompleteMax" -> ESCompleteMax[d],
      "Scale" -> N[seriesScale[d], 6],
      "LaurentLeadTol" -> DiffExp2`Tolerances`Tol["LaurentLeadTol"],
      "Detail" -> "denominator is entirely negligible"|>]];
  {L, dL} = lead;
  relCM = ESCompleteMax[d] - 2 L;  (* result window: [-L, dCM - 2L] *)
  e = Table[0, {relCM + L + 1}];   (* relative orders 0 .. relCM + L *)
  e[[1]] = norm[1/dL];
  Do[e[[m + 1]] = norm[-(1/dL) Sum[ESCoefficient[d, L + j]*e[[m - j + 1]], {j, 1, m}]],
    {m, 1, relCM + L}];
  mkSeries[-L, relCM, e]];

ESDivide[a_?ESQ, b_?ESQ] := ESTimes[a, ESInvert[b]];

(* Worked shift example (normative comment per spec 3.3): dividing 1 + O(eps^6)
   (window [0,5]) by 3 eps + O(eps^7) (window [1,6]) gives L = 1, window
   [-1, Min[5-1, 0+6-2]] = [-1, 4], coefficients {1/3, 0, 0, 0, 0, 0}:
   exact 1/3, BOTH edges shifted down by one. *)

ESTruncate[s_?ESQ, newCM_Integer] := Which[
  newCM > ESCompleteMax[s],
  esError["ERR-TRUNCATE-EXTEND", <|"Requested" -> newCM,
    "Min" -> ESMinPower[s], "CompleteMax" -> ESCompleteMax[s],
    "Detail" -> "truncation cannot EXTEND the window (no padding exists)"|>],
  newCM < ESMinPower[s],
  esError["ERR-TRUNCATE-EXTEND", <|"Requested" -> newCM,
    "Min" -> ESMinPower[s], "CompleteMax" -> ESCompleteMax[s],
    "Detail" -> "truncation below Min would empty the window"|>],
  True,
  mkSeries[ESMinPower[s], newCM, Take[s["Coeffs"], newCM - ESMinPower[s] + 1]]];

ESTrim[s_?ESQ] := Module[{scale = seriesScale[s], kmin = ESMinPower[s], first = None},
  Do[If[!ESCoeffZeroQ[s["Coeffs"][[i]], scale], first = i; Break[]],
    {i, Length[s["Coeffs"]]}];
  If[first === None, ESZero[ESCompleteMax[s]],
    mkSeries[kmin + first - 1, ESCompleteMax[s], Drop[s["Coeffs"], first - 1]]]];

(* ---- comparison ---- *)

ESSameQ[a_?ESQ, b_?ESQ] := Module[
  {k1 = Min[ESMinPower[a], ESMinPower[b]],
   k2 = Min[ESCompleteMax[a], ESCompleteMax[b]], scale, mtol, diff},
  If[k2 < k1, Return[True]];  (* shared range inside both zero certificates *)
  scale = Max[seriesScale[a], seriesScale[b]];
  mtol = DiffExp2`Tolerances`Tol["MatchTol"];
  (* "Never errors; returns False honestly" (spec): an ambiguity-band abort
     inside the classification means not-provably-equal -> False *)
  AllTrue[Range[k1, k2], Function[k,
    diff = ESCoefficient[a, k] - ESCoefficient[b, k];
    If[NumericQ[diff],
      TrueQ[Catch[ESCoeffZeroQ[diff, scale, mtol], "DiffExp2Error"]],
      TrueQ[PossibleZeroQ[diff]]]]]];

(* ---- conversion (API boundary only) ---- *)

ESFromExpression[expr_, epsSym_Symbol, kmax_Integer] := Module[
  {ser, nmin, nmax, den, raw, coeffs, cm},
  If[PossibleZeroQ[expr], Return[ESZero[kmax]]];
  If[FreeQ[expr, epsSym],
    Return[mkSeries[0, kmax, validateCoeffs["ESFromExpression",
      Join[{Together[expr]}, Table[0, {kmax}]]]]]];
  ser = Quiet[Check[Series[expr, {epsSym, 0, kmax}], $Failed]];
  If[ser === $Failed || Head[ser] =!= SeriesData,
    esError["ERR-EXPAND-FAIL", <|"Expression" -> StringTake[ToString[expr, InputForm], UpTo[200]],
      "EpsSymbol" -> epsSym, "RequestedMax" -> kmax,
      "Got" -> Head[ser], "Detail" -> "Series did not produce a SeriesData expansion"|>]];
  {nmin, nmax, den} = {ser[[4]], ser[[5]], ser[[6]]};
  If[den =!= 1,
    esError["ERR-EXPAND-FAIL", <|"Expression" -> StringTake[ToString[expr, InputForm], UpTo[200]],
      "EpsSymbol" -> epsSym, "LeadingPower" -> nmin/den,
      "Detail" -> "fractional eps-power expansion (old code Floor-ed this; forbidden F7)"|>]];
  raw = PadRight[ser[[3]], nmax - nmin, 0];
  cm = Min[nmax - 1, kmax];
  If[cm < nmin, Return[ESZero[cm]]];
  coeffs = Together /@ Take[raw, cm - nmin + 1];
  If[!AllTrue[coeffs, FreeQ[#, epsSym] &],
    esError["ERR-EXPAND-FAIL", <|"Expression" -> StringTake[ToString[expr, InputForm], UpTo[200]],
      "EpsSymbol" -> epsSym, "Detail" -> "residual eps inside an extracted coefficient"|>]];
  mkSeries[nmin, cm, validateCoeffs["ESFromExpression", coeffs]]];

ESToExpression[s_?ESQ, epsSym_Symbol] :=
  Sum[ESCoefficient[s, k]*epsSym^k, {k, ESMinPower[s], ESCompleteMax[s]}];

ESMap[f_, s_?ESQ] :=
  mkSeries[ESMinPower[s], ESCompleteMax[s], validateCoeffs["ESMap", f /@ s["Coeffs"]]];

End[];
EndPackage[];
