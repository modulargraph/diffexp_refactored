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
numMag = DiffExp2`Tolerances`NumericMagnitude;
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
   F(T) = T^(m+1) T^(b eps) Sum_j (-1)^j eps^p Log^(p-j)T / ((p-j)! alpha^(j+1))
   MEMOIZED on the exact argument tuple (+ WP, which sets the numeric log):
   a tile's base towers are identical across components, masters and cvec
   entries — LineIntegral re-integrates the same tile once per nonzero
   coefficient, and the per-n towers dominated the tile profile
   (Docs/SpeedIdeas.md R1).  Pure function of exact inputs: same values,
   bounded cache (flush at cap, the $shCache pattern). *)
$adCache = <||>; $adCacheMax = 4096;
antiderivativeAt[m_, b_, p_, T_, kMaxOut_] := Module[
  {key = {m, b, p, T, kMaxOut, DiffExp2`Config`CFG["WorkingPrecision"]}},
  If[KeyExistsQ[$adCache, key], $adCache[key],
    If[Length[$adCache] >= $adCacheMax, $adCache = <||>];
    $adCache[key] = antiderivativeAtCore[m, b, p, T, kMaxOut]]];

antiderivativeAtCore[m_, b_, p_, T_, kMaxOut_] := Module[{logT, wp},
  wp = DiffExp2`Config`CFG["WorkingPrecision"];
  (* numeric log: exact Log[rational] towers compound into symbolic giants *)
  logT = If[FreeQ[T, _Symbol],
    N[Log[T], DiffExp2`Tolerances`$InputPrecisionFactor*wp], Log[T]];
  antiderivativeAtLog[m, b, p, T, logT, 1, kMaxOut]];

(* Antiderivative at a point of modulus T with an already branch-resolved
   logarithm logT.  phase0 is the eps-independent part of t^(m+1+b eps).
   Keeping this sibling of antiderivativeAtCore lets the interior-crossing
   path evaluate the negative endpoint with the full i-delta phase, without
   sending branch data through a symbolic Series[]. *)
antiderivativeAtLog[m_, b_, p_, T_, logT_, phase0_, kMaxOut_] := Module[
  {width, tbeps, invAlpha, jsum, alphaES},
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
  esScale[phase0*pow[T, m + 1], esTimes[tbeps, esShift[jsum, p]]]];

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
      scale = Max[1, Sequence @@ (numMag[#, 15] & /@
        Select[rowsAll, NumericQ])];
      Do[
        If[!TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[tot[[c]], scale,
            DiffExp2`Tolerances`Tol["LaurentLeadTol"],
            "Integrate cancellation gate",
            DiffExp2`Tolerances`$AmbiguityBandDecades,
            SameQ[DiffExp2`Tolerances`Tol["LaurentLeadTol"], 10^-24]]],
          err["E2", <|"m" -> grp[[1, 3]], "p" -> grp[[1, 4]], "EpsRow" -> k,
            "MergedCoefficient" -> N[tot[[c]], 6], "Scale" -> scale,
            "Detail" -> "divergent b = 0 endpoint content does not cancel in the assembled combination"|>]],
        {c, Length[tot]}]],
      {k, kmin, kmax}]],
    {g, GatherBy[offenders, {Together[#[[3]]], #[[4]]} &]}];
  offenders];

requireEndpointIntegralOrder[ls_, ncols_Integer] := Module[{missing = {}},
  Do[If[zeroQ[sec["b"]] && TrueQ[sec["a"] + 1 <= 0],
    Module[{needed = Floor[-sec["a"] - 1]},
      If[IntegerQ[needed] && needed >= ncols,
        AppendTo[missing, <|"Sector" -> {sec["a"], sec["b"], sec["p"]},
          "NeededTOrder" -> needed|>]]]],
    {sec, ls["Sectors"]}];
  If[missing =!= {},
    err["E10", <|"Sectors" -> missing, "AvailableTOrder" -> ncols - 1,
      "Detail" -> "endpoint integral is not certified: unseen Taylor cells may still be nonintegrable"|>]]];

(* ---- 2.2 LocalSolution integration ---- *)

IntegrateLocalSolution[ls0_Association, {t1_, t2_}] := Module[
  {ls = DiffExp2`SectorSeries`ValidateLocalSolution[ls0], kmin, kmax, ncols,
   ncomp, vals},
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
      requireEndpointIntegralOrder[ls, ncols];
      vals = sumSectorIntegrals[ls, t2, kmin, kmax, ncols, ncomp, skip]],
    (* center endpoint mirrored: [t1, 0] = crossing-evaluated arm *)
    TrueQ[t2 == 0],
    Module[{skip = divergentGate[ls, kmin, kmax], lsX, sg},
      requireEndpointIntegralOrder[ls, ncols];
      sg = resolveSigma[ls];
      lsX = DiffExp2`Transport`ApplyCrossing[ls, sg];
      (* Int_{t1}^0 f dt = Int_0^{|t1|} f(-u) du: crossing-applied, positive *)
      vals = sumSectorIntegrals[lsX, -t1, lsX["EpsWindow", "Min"],
        lsX["EpsWindow", "CompleteMax"], ncols, ncomp,
        divergentGate[lsX, lsX["EpsWindow", "Min"], lsX["EpsWindow", "CompleteMax"]]]],
    (* interior crossing: [t1, 0] + [0, t2], PV-paired *)
    TrueQ[t1 < 0 < t2],
    Module[{skip},
      skip = interiorGate[ls, kmin, kmax];
      vals = sumInteriorPaired[ls, -t1, t2, kmin, kmax, ncols, ncomp]],
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
  {span = kmax - kmin, acc},
  acc = Table[None, {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    Do[Module[{m = Together[a + n], base},
      If[zeroQ[b] && TrueQ[m + 1 <= 0],
        Null,  (* gated: merged coefficient verified zero *)
        (* The primitive is multiplied by a coefficient window of width
           span+1.  Its private depth is relative to that width, not to the
           absolute (possibly very negative) CompleteMax. *)
        base = SectorMonomialIntegral[m, b, p, T, span + p + 2];
        Do[Module[{cESc},
          cESc = esNew[kmin, Table[arr[[k - kmin + 1, n + 1, c]], {k, kmin, kmax}]];
          Module[{term = esTimes[cESc, base]},
            acc[[c]] = If[acc[[c]] === None, term, esAdd[acc[[c]], term]]]],
          {c, ncomp}]]],
      {n, 0, ncols - 1}]],
    {sec, ls["Sectors"]}];
  Map[If[# === None, esZero[kmax], #] &, acc]];

(* plain antiderivative sum at T > 0 (no boundary subtleties) *)
sumAntiderivative[ls_, T_, kmin_, kmax_, ncols_, ncomp_] := Module[
  {acc, span = kmax - kmin},
  acc = Table[None, {ncomp}];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    Do[Module[{m = Together[a + n], base},
      base = If[zeroQ[b] && zeroQ[m + 1],
        (* alpha == 0: F = eps^p Log^(p+1) t/((p+1) p!) *)
        esShift[esNew[0, PadRight[{pow[Log[T], p + 1]/((p + 1)*p!)}, kmax - kmin + 4]], p],
        DiffExp2`EpsSeries`ESTruncate[
          antiderivativeAt[m, b, p, T, span + p + 2], span + p + 2]];
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

(* Strict sign lookup for a branch-sensitive interior crossing.  Unlike the
   regular transport-side convention, an i-delta integral has no meaningful
   default side: the two choices are complex conjugates. *)
interiorSigma[ls_] := Module[{sg = DiffExp2`SectorSeries`ChartImSign[ls]},
  If[!MemberQ[{1, -1}, sg],
    err["E3", <|"Center" -> ls["Center"],
      "Prescriptions" -> ls["Prescriptions"],
      "Detail" -> "b != 0 interior crossing requires a derivable i-delta prescription sign"|>]];
  sg];

(* The pole cell m=-1 must be paired BEFORE epsilon-window arithmetic.  With
     LB = Log B, LN = Log A + i Pi sigma, beta = b eps,
     H(beta) = (Exp[beta LB] - Exp[beta LN])/beta,
   the normalized p-log monomial is eps^p H^(p)(beta)/p!.  Its manifestly
   regular coefficients are
     eps^(p+q) b^q (LB^(p+q+1)-LN^(p+q+1)) /
       (p! q! (p+q+1)), q>=0.
   In particular p=0 has no spurious eps^-1 term. *)
pairedPoleIntegral[b_, p_, A_, B_, sigma_, width_Integer] := Module[
  {wp = DiffExp2`Config`CFG["WorkingPrecision"], logB, logNeg},
  logB = If[FreeQ[B, _Symbol],
    N[Log[B], DiffExp2`Tolerances`$InputPrecisionFactor*wp], Log[B]];
  logNeg = If[FreeQ[A, _Symbol],
      N[Log[A], DiffExp2`Tolerances`$InputPrecisionFactor*wp], Log[A]] +
    sigma*I*Pi;
  esNew[p, Table[
    pow[b, q]*(pow[logB, p + q + 1] - pow[logNeg, p + q + 1])/
      (p!*q!*(p + q + 1)),
    {q, 0, width}]]];

(* Branch-resolved paired form away from the pole cell.  Since m+1 != 0,
   both endpoint antiderivatives are epsilon-regular and their direct
   difference cannot manufacture a Laurent pole. *)
pairedRegularIntegral[m_, b_, p_, A_, B_, sigma_, kMaxOut_] := Module[
  {wp = DiffExp2`Config`CFG["WorkingPrecision"], logNeg, pos, neg},
  logNeg = If[FreeQ[A, _Symbol],
      N[Log[A], DiffExp2`Tolerances`$InputPrecisionFactor*wp], Log[A]] +
    sigma*I*Pi;
  pos = antiderivativeAt[m, b, p, B, kMaxOut];
  neg = antiderivativeAtLog[m, b, p, A, logNeg,
    Exp[sigma*I*Pi*(m + 1)], kMaxOut];
  esAdd[pos, esScale[-1, neg]]];

(* Single-owner interior pairing.  b=0 retains the real-log PV/Hadamard
   finite-part convention exactly.  b!=0 uses the chart prescription on the
   negative arm, with the m=-1 pole cell expanded in its combined regular
   form above rather than as two independently shifted Laurent series. *)
sumInteriorPaired[ls_, A_, B_, kmin_, kmax_, ncols_, ncomp_] := Module[
  {acc = Table[None, {ncomp}], sigma = None, needSigma, span = kmax - kmin},
  needSigma = AnyTrue[ls["Sectors"], Function[sec,
    !zeroQ[sec["b"]] && !AllTrue[Flatten[sec["Coeffs"]], # === 0 &]]];
  If[needSigma, sigma = interiorSigma[ls]];
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    If[zeroQ[b] && !IntegerQ[a] &&
        !AllTrue[Flatten[arr], # === 0 &],
      err["E8", <|"a" -> a,
        "Detail" -> "real-log PV for a fractional-a b=0 sector is undefined"|>]];
    Do[Module[{m = Together[a + n], base},
      (* An exact-zero coefficient slab is structurally inactive: it neither
         requires a branch sign nor constrains the result window. *)
      If[AllTrue[Flatten[arr[[All, n + 1, All]]], # === 0 &], Continue[]];
      base = Which[
        zeroQ[b] && zeroQ[m + 1],
          esShift[esNew[0, PadRight[{
            (pow[Log[B], p + 1] - pow[Log[A], p + 1])/((p + 1)*p!)},
            span + 4]], p],
        zeroQ[b],
          (* Int_{-A}^B with real Log|t|: the negative-arm substitution
             contributes (-1)^m.  For m<-1 this is the defined Hadamard
             finite part; for m>-1 it is the ordinary integral. *)
          esAdd[antiderivativeAt[m, 0, p, B, span + p + 2],
            esScale[(-1)^m, antiderivativeAt[m, 0, p, A, span + p + 2]]],
        zeroQ[m + 1],
          pairedPoleIntegral[b, p, A, B, sigma, span + 3],
        True,
          (* The base-series depth is relative to the coefficient-window
             width, not to its absolute (possibly negative) upper order. *)
          pairedRegularIntegral[m, b, p, A, B, sigma, p + span + 3]];
      Do[Module[{cESc, term},
        cESc = esNew[kmin,
          Table[arr[[k - kmin + 1, n + 1, c]], {k, kmin, kmax}]];
        term = esTimes[cESc, base];
        acc[[c]] = If[acc[[c]] === None, term, esAdd[acc[[c]], term]]],
        {c, ncomp}]],
      {n, 0, ncols - 1}]],
    {sec, ls["Sectors"]}];
  Map[If[# === None, esZero[kmax], #] &, acc]];

(* ---- 2.4 endpoint limit (the FT drop rule) ---- *)

EndpointSectorLimit[ls0_Association] := Module[
  {ls, kmin, kmax, ncomp, ncols, finite = {}, divergent = {},
   incomplete = {}, finiteRows},
  (* Work on the original sectors and merge by ABSOLUTE monomial below.
     CanonicalizeLocalSolution uses a fixed Taylor width, so shifting a
     higher-a tower into a much lower-a tower can legitimately place its
     finite t^0 cell beyond that width.  Endpoint classification must not
     discard such a cell merely to normalize the representation. *)
  ls = DiffExp2`SectorSeries`ValidateLocalSolution[ls0];
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
  (* Classify each coefficient by its ABSOLUTE monomial power m=a+n.
     Positive powers vanish even when their sector's leading a is negative;
     this is what permits a finite t^0 coefficient to be buried at n=-a. *)
  Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"]},
    If[zeroQ[b],
      (* A negative-a tower needs every coefficient through n=floor(-a)
         to certify its limit.  A finite zero prefix is not an exact-zero
         certificate for the unknown tail.  Record an incomplete tower now;
         report it only after any already-visible divergence has had the
         chance to raise the more specific E2 below. *)
      If[TrueQ[a < 0],
        Module[{needed = Floor[-a]},
          If[IntegerQ[needed] && needed >= ncols,
            AppendTo[incomplete, <|"Sector" -> {a, b, p},
              "NeededTOrder" -> needed|>]]]];
      Do[Module[{m = Together[a + n], cell = {sec, n}},
        Which[
          TrueQ[m > 0], Null,
          zeroQ[m] && p === 0, AppendTo[finite, cell],
          TrueQ[m < 0] || (zeroQ[m] && p > 0),
            AppendTo[divergent, {m, p, sec, n}],
          True,
            err["E10", <|"Sector" -> {a, b, p}, "AbsolutePower" -> m,
              "Detail" -> "could not classify exact endpoint monomial power"|>]]],
        {n, 0, ncols - 1}],
    Null  (* b != 0: dropped exactly by the dimensional-regulator rule *)]],
    {sec, ls["Sectors"]}];
  (* Divergent monomials may cancel only against the SAME absolute power
     and Log depth.  Check their merged coefficient per epsilon row and
     component; unrelated positive-power and finite cells never enter this
     gate. *)
  Do[Module[{grp = g},
    Do[
      Do[Module[{terms, total, scale},
        terms = Map[#[[3]]["Coeffs"][[k - kmin + 1, #[[4]] + 1, c]] &, grp];
        total = Total[terms];
        scale = Max[1, Sequence @@ (numMag[#, 15] & /@
          Select[terms, NumericQ])];
        If[!TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[total, scale,
              DiffExp2`Tolerances`Tol["LaurentLeadTol"],
              "EndpointSectorLimit cancellation gate",
              DiffExp2`Tolerances`$AmbiguityBandDecades,
              SameQ[DiffExp2`Tolerances`Tol["LaurentLeadTol"], 10^-24]]],
          err["E2", <|"AbsolutePower" -> grp[[1, 1]],
            "LogPower" -> grp[[1, 2]], "EpsRow" -> k, "Component" -> c,
            "MergedCoefficient" -> total, "Scale" -> scale,
            "SectorTags" -> ({#[[3]]["a"], #[[3]]["b"], #[[3]]["p"]} & /@ grp),
            "Detail" -> "divergent b = 0 content at the endpoint limit does not cancel"|>]]],
        {c, ncomp}],
      {k, kmin, kmax}]],
    {g, GatherBy[divergent, {Together[#[[1]]], #[[2]]} &]}];
  If[incomplete =!= {},
    err["E10", <|"Sectors" -> incomplete,
      "AvailableTOrder" -> ncols - 1,
      "Detail" -> "endpoint limit is not certified: a negative-power tower needs Taylor coefficients beyond the stored T window"|>]];
  If[finite === {}, Return[Table[esZero[kmax], {ncomp}], Module]];
  finiteRows = Total[Map[#[[1]]["Coeffs"][[All, #[[2]] + 1]] &, finite]];
  Table[esNew[kmin, finiteRows[[All, c]]], {c, ncomp}]];

End[];
EndPackage[];
