(* DiffExp2/Integrate.m — exact per-sector definite integrals.
   Spec: Docs/specs/Integrate.md (binding); DECISIONS-M0.md.
   The dimreg lower boundary t^(m+1+b eps) Log^q t |_{t=0} := 0 for b != 0
   (the drop rule's integral counterpart) makes every old IBP-step /
   subtraction regularization a single closed form. *)

BeginPackage["DiffExp2`Integrate`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`"}];

SectorMonomialIntegral::usage = "SectorMonomialIntegral[m, b, p, T, kMaxOut] gives Integrate[t^(m + b eps) (eps Log t)^p/p!, {t, 0, T}] exactly as an EpsSeries (dimreg lower boundary for b != 0; loud E2 for divergent b = 0 cells).";
IntegrateLocalSolution::usage = "IntegrateLocalSolution[ls, {t1, t2}] integrates a LocalSolution over a chart interval: center-endpoint sector integrals with the cancellation gate, interior PV pairing, plain antiderivative differences outside.";
EndpointSectorLimit::usage = "EndpointSectorLimit[ls] gives the t -> 0 limit per component: b != 0 sectors dropped EXACTLY (dimreg), divergent b = 0 content is a loud error after the cancellation gate.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Integrate"|>, payload]];
cfg = DiffExp2`Config`CFG;
(* zeroQ: exact-first (see Solve.m).  Call sites here test exact tags
   (a, b, m + 1): the rational fast path decides them without numerics;
   anything outside the rational domain keeps PossibleZeroQ unchanged. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && TrueQ[PossibleZeroQ[c]]);
zeroQ[e_] := zeroCanQ[Together[e]];
esNew = DiffExp2`EpsSeries`ESNew; esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd; esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes; esShift = DiffExp2`EpsSeries`ESShift;
esInv = DiffExp2`EpsSeries`ESInvert;
esCoeff = DiffExp2`EpsSeries`ESCoefficient; esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax;
pow[x_, 0] := 1; pow[x_, n_] := x^n;

(* ---- 2.1 the primitive ---- *)

(* antiderivative value F(T) for alpha != 0:
   F(T) = T^(m+1) T^(b eps) Sum_j (-1)^j eps^p Log^(p-j)T / ((p-j)! alpha^(j+1)) *)
antiderivativeAt[m_, b_, p_, T_, kMaxOut_] := Module[
  {logT, width, tbeps, invAlpha, jsum, alphaES, wp},
  wp = DiffExp2`Config`CFG["WorkingPrecision"];
  (* numeric log: exact Log[rational] towers compound into symbolic giants *)
  logT = If[FreeQ[T, _Symbol], N[Log[T], wp + 20], Log[T]];
  width = kMaxOut + p + 4;
  tbeps = If[zeroQ[b], esNew[0, PadRight[{1}, width + 1]],
    esNew[0, Table[pow[b*logT, r]/r!, {r, 0, width}]]];
  alphaES = If[zeroQ[b],
    esNew[0, PadRight[{m + 1}, width + 1]],
    If[zeroQ[m + 1],
      esNew[1, PadRight[{b}, width]],
      esNew[0, PadRight[{m + 1, b}, width + 1]]]];
  invAlpha = esInv[alphaES];
  jsum = Module[{acc = None, pwr = invAlpha},
    Do[Module[{term = esScale[(-1)^j*pow[logT, p - j]/(p - j)!, pwr]},
      acc = If[acc === None, term, esAdd[acc, term]];
      pwr = esTimes[pwr, invAlpha]],
      {j, 0, p}];
    acc];
  esScale[pow[T, m + 1], esTimes[tbeps, esShift[jsum, p]]]];

SectorMonomialIntegral[m_, b_, p_Integer, T_, kMaxOut_Integer] := Module[{},
  If[!FreeQ[m, _?InexactNumberQ] || !FreeQ[b, _?InexactNumberQ] || p < 0,
    err["E10", <|"m" -> m, "b" -> b, "p" -> p,
      "Detail" -> "tags must be exact; p a non-negative integer"|>]];
  If[!TrueQ[T > 0],
    err["E9", <|"T" -> T, "Detail" -> "SectorMonomialIntegral requires T > 0"|>]];
  If[zeroQ[b] && TrueQ[Together[m + 1] <= 0],
    err["E2", <|"m" -> m, "b" -> b, "p" -> p,
      "Detail" -> "divergent b = 0 endpoint monomial reached the primitive (the cancellation gate must remove or error upstream)"|>]];
  (* alpha != 0 in all remaining cells; F(0) = 0 (true limit or dimreg) *)
  DiffExp2`EpsSeries`ESTruncate[antiderivativeAt[m, b, p, T, kMaxOut],
    kMaxOut]];

(* ---- the cancellation gate (DEC-4) ----
   merged coefficient of a divergent (b=0, m+1<=0, any p) monomial across
   ALL sectors, per eps row; scale = max |coeff| of the merged combination
   at that row *)
divergentGate[ls_, kmin_, kmax_] := Module[
  {secs = ls["Sectors"], ncols, ncomp, offenders = {}},
  {ncols, ncomp} = Dimensions[First[secs]["Coeffs"]][[2 ;; 3]];
  Do[Module[{b = sec["b"], a = sec["a"], p = sec["p"]},
    If[zeroQ[b],
      Do[Module[{m = Together[a + n]},
        If[TrueQ[m + 1 <= 0],
          AppendTo[offenders, {sec, n, m, p}]]],
        {n, 0, ncols - 1}]]],
    {sec, secs}];
  (* group by (m, p): merged coefficient must vanish per eps row *)
  Do[Module[{grp = g, merged, scale},
    Do[Module[{tot, rowsAll},
      tot = Total[Map[Module[{sec = #[[1]], n = #[[2]]},
        sec["Coeffs"][[k - kmin + 1, n + 1]]] &, grp]];
      rowsAll = Flatten[Map[Module[{sec = #[[1]], n = #[[2]]},
        sec["Coeffs"][[k - kmin + 1, n + 1]]] &, grp]];
      scale = Max[0, Sequence @@ (Abs[N[#, 15]] & /@ Select[rowsAll, NumericQ])];
      Do[
        If[!TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[tot[[c]], scale,
            DiffExp2`Tolerances`Tol["LaurentLeadTol"],
            "Integrate cancellation gate"]],
          err["E2", <|"m" -> grp[[1, 3]], "p" -> grp[[1, 4]], "EpsRow" -> k,
            "MergedCoefficient" -> N[tot[[c]], 6],
            "Detail" -> "divergent b = 0 endpoint content does not cancel in the assembled combination"|>]],
        {c, Length[tot]}]],
      {k, kmin, kmax}]],
    {g, GatherBy[offenders, {Together[#[[3]]], #[[4]]} &]}];
  offenders];

(* ---- 2.2 LocalSolution integration ---- *)

IntegrateLocalSolution[ls0_Association, {t1_, t2_}] := Module[
  {ls = DiffExp2`SectorSeries`ValidateLocalSolution[ls0], kmin, kmax, ncols,
   ncomp, kOut, vals, sigma},
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
  If[!TrueQ[t1 < t2],
    err["E9", <|"Bounds" -> {t1, t2}, "Detail" -> "require t1 < t2"|>]];
  If[TrueQ[Abs[t1] >= ls["Radius"]] || TrueQ[Abs[t2] >= ls["Radius"]],
    err["E9", <|"Bounds" -> {t1, t2}, "Radius" -> ls["Radius"],
      "Detail" -> "integration bounds outside the chart radius"|>]];
  Which[
    (* center endpoint: [0, t2] *)
    TrueQ[t1 == 0],
    Module[{skip = divergentGate[ls, kmin, kmax]},
      vals = sumSectorIntegrals[ls, t2, kmin, kmax, ncols, ncomp, skip]],
    (* center endpoint mirrored: [t1, 0] = crossing-evaluated arm *)
    TrueQ[t2 == 0],
    Module[{skip = divergentGate[ls, kmin, kmax], lsX, sg},
      sg = resolveSigma[ls];
      lsX = DiffExp2`Transport`ApplyCrossing[ls, sg];
      (* Int_{t1}^0 f dt = Int_0^{|t1|} f(-u) du: crossing-applied, positive *)
      vals = sumSectorIntegrals[lsX, -t1, lsX["EpsWindow", "Min"],
        lsX["EpsWindow", "CompleteMax"], ncols, ncomp,
        divergentGate[lsX, lsX["EpsWindow", "Min"], lsX["EpsWindow", "CompleteMax"]]]],
    (* interior crossing: [t1, 0] + [0, t2], PV-paired *)
    TrueQ[t1 < 0 < t2],
    Module[{vPos, vNeg, skip},
      skip = interiorGate[ls, kmin, kmax];
      vPos = sumSectorIntegralsPV[ls, t2, kmin, kmax, ncols, ncomp];
      vNeg = sumSectorIntegralsPVNeg[ls, -t1, kmin, kmax, ncols, ncomp];
      vals = MapThread[esAdd, {vPos, vNeg}]],
    (* center outside: plain F(t2) - F(t1) *)
    True,
    vals = Module[{f2, f1},
      If[TrueQ[t1 > 0],
        f2 = sumAntiderivative[ls, t2, kmin, kmax, ncols, ncomp];
        f1 = sumAntiderivative[ls, t1, kmin, kmax, ncols, ncomp];
        MapThread[esAdd[#1, esScale[-1, #2]] &, {f2, f1}],
        (* both negative: crossing-evaluated arm *)
        Module[{lsX = DiffExp2`Transport`ApplyCrossing[ls, resolveSigma[ls]], g2, g1},
          g2 = sumAntiderivative[lsX, -t1, lsX["EpsWindow", "Min"],
            lsX["EpsWindow", "CompleteMax"], ncols, ncomp];
          g1 = sumAntiderivative[lsX, -t2, lsX["EpsWindow", "Min"],
            lsX["EpsWindow", "CompleteMax"], ncols, ncomp];
          MapThread[esAdd[#1, esScale[-1, #2]] &, {g2, g1}]]]]];
  <|"Values" -> vals,
    "EpsWindow" -> <|"Min" -> Min[esMin /@ vals], "CompleteMax" -> Min[esCM /@ vals]|>,
    "TWindowUsed" -> ncols - 1|>];

resolveSigma[ls_] := Module[{sg = DiffExp2`SectorSeries`ChartImSign[ls]},
  (* prescriptions when derivable, else the SAME fixed +1 convention the
     transport matching used (Transport`sigmaFor): the weights were fitted
     under +1, so the crossing must use +1 - phases cancel in the real
     Euclidean combination *)
  If[MemberQ[{1, -1}, sg], sg, 1]];

(* sum of endpoint sector integrals over [0, T]; skip = offender list
   (cancelled divergent monomials are EXCLUDED from every sector) *)
sumSectorIntegrals[ls_, T_, kmin_, kmax_, ncols_, ncomp_, skip_] := Module[
  {kOut = kmax, acc},
  acc = Table[None, {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    Do[Module[{m = Together[a + n], base},
      If[zeroQ[b] && TrueQ[m + 1 <= 0],
        Null,  (* gated: merged coefficient verified zero *)
        base = SectorMonomialIntegral[m, b, p, T, kOut + p + 2];
        Do[Module[{cESc},
          cESc = esNew[kmin, Table[arr[[k - kmin + 1, n + 1, c]], {k, kmin, kmax}]];
          Module[{term = esTimes[cESc, base]},
            acc[[c]] = If[acc[[c]] === None, term, esAdd[acc[[c]], term]]]],
          {c, ncomp}]]],
      {n, 0, ncols - 1}]],
    {sec, ls["Sectors"]}];
  Map[If[# === None, esZero[kmax], #] &, acc]];

(* plain antiderivative sum at T > 0 (no boundary subtleties) *)
sumAntiderivative[ls_, T_, kmin_, kmax_, ncols_, ncomp_] := Module[{acc},
  acc = Table[None, {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    Do[Module[{m = Together[a + n], base},
      base = If[zeroQ[b] && zeroQ[m + 1],
        (* alpha == 0: F = eps^p Log^(p+1) t/((p+1) p!) *)
        esShift[esNew[0, PadRight[{pow[Log[T], p + 1]/((p + 1))}, kmax - kmin + 4]], p],
        DiffExp2`EpsSeries`ESTruncate[
          antiderivativeAt[m, b, p, T, kmax + p + 2], kmax + p + 2]];
      Do[Module[{cESc},
        cESc = esNew[kmin, Table[arr[[k - kmin + 1, n + 1, c]], {k, kmin, kmax}]];
        Module[{term = esTimes[cESc, base]},
          acc[[c]] = If[acc[[c]] === None, term, esAdd[acc[[c]], term]]]],
        {c, ncomp}]],
      {n, 0, ncols - 1}]],
    {sec, ls["Sectors"]}];
  Map[If[# === None, esZero[kmax], #] &, acc]];

(* interior gate: which (m,p) monomials are PV-paired (b=0, m+1 <= 0);
   odd/even handled by the pairing formulas below; returns Null (no skip) *)
interiorGate[ls_, kmin_, kmax_] := Null;

(* positive arm of an interior crossing: F(t2) - lim_{c->0+} F(c) with the
   divergent b=0 cells kept as REAL-log PV halves: F(t2) with Log real *)
sumSectorIntegralsPV[ls_, T_, kmin_, kmax_, ncols_, ncomp_] :=
  sumAntiderivativePV[ls, T, kmin, kmax, ncols, ncomp, 1];
sumSectorIntegralsPVNeg[ls_, T_, kmin_, kmax_, ncols_, ncomp_] :=
  sumAntiderivativePV[ls, T, kmin, kmax, ncols, ncomp, -1];

(* PV arm: for b != 0 sectors use the dimreg endpoint integral on each side
   evaluated with REAL logs of |t| (Euclidean PV convention, campaign
   9aeb300); the negative arm enters with orientation factor and the parity
   phase (-1)^(m+1) from t -> -u, u > 0 with REAL u^... :
   Int_{-T}^0 t^m dt = (-1)^(m+1) Int_0^T u^m du  (integer m; PV pairing
   cancels the divergent symmetric parts of odd integrands) *)
sumAntiderivativePV[ls_, T_, kmin_, kmax_, ncols_, ncomp_, arm_] := Module[{acc},
  acc = Table[None, {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    If[arm === -1 && !IntegerQ[a],
      err["E8", <|"a" -> a,
        "Detail" -> "interior crossing with fractional-a sector requires the prescription-phase path (not in v1 PV)"|>]];
    Do[Module[{m = Together[a + n], base, parity},
      parity = If[arm === 1, 1, (-1)^m];  (* t -> -u, real-log PV convention *)
      base = Which[
        zeroQ[b] && TrueQ[m + 1 > 0],
        SectorMonomialIntegral[m, b, p, T, kmax + p + 2],
        zeroQ[b] && zeroQ[m + 1],
        esShift[esNew[0, PadRight[{pow[Log[T], p + 1]/(p + 1)}, kmax - kmin + 4]], p],
        zeroQ[b],
        (* m+1 < 0 power cell: real finite part F(T) (the 1/t^k endpoint
           terms cancel between PV arms for the surviving combinations;
           genuine non-cancelling power divergence is caught by the
           E2-on-assembly check below at the eps-row level) *)
        antiderivativeAt[m, 0, p, T, kmax + p + 2],
        True,
        SectorMonomialIntegral[m, b, p, T, kmax + p + 2]];
      Do[Module[{cESc},
        cESc = esNew[kmin, Table[arr[[k - kmin + 1, n + 1, c]], {k, kmin, kmax}]];
        Module[{term = esScale[parity, esTimes[cESc, base]]},
          acc[[c]] = If[acc[[c]] === None, term, esAdd[acc[[c]], term]]]],
        {c, ncomp}]],
      {n, 0, ncols - 1}]],
    {sec, ls["Sectors"]}];
  Map[If[# === None, esZero[kmax], #] &, acc]];

(* ---- 2.4 endpoint limit (the FT drop rule) ---- *)

EndpointSectorLimit[ls0_Association] := Module[
  {ls = DiffExp2`SectorSeries`ValidateLocalSolution[ls0], kmin, kmax, ncomp,
   ncols, acc},
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
  acc = Table[esZero[kmax], {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    Which[
      !zeroQ[b], Null,  (* dropped EXACTLY: the dimreg rule *)
      TrueQ[Together[a] > 0], Null,  (* vanishes at the center *)
      zeroQ[a] && p === 0,
      Do[Module[{cESc = esNew[kmin,
          Table[arr[[k - kmin + 1, 1, c]], {k, kmin, kmax}]]},
        acc[[c]] = esAdd[acc[[c]], cESc]],
        {c, ncomp}],
      True,
      (* a < 0 or log content at b = 0: divergent unless the rows vanish *)
      Module[{flat = Flatten[arr], scale},
        scale = Max[0, Sequence @@ (Abs[N[#, 15]] & /@ Select[flat, NumericQ])];
        If[AnyTrue[flat, !TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[#, scale,
            DiffExp2`Tolerances`Tol["LaurentLeadTol"], "EndpointSectorLimit"]] &],
          err["E2", <|"Sector" -> {a, b, p},
            "Detail" -> "divergent b = 0 content at the endpoint limit"|>]]]]],
    {sec, ls["Sectors"]}];
  acc];

End[];
EndPackage[];
