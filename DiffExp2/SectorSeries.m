(* DiffExp2/SectorSeries.m — the LocalSolution object and its closed algebra.
   Spec: Docs/specs/SectorSeries.md (binding); decisions: DECISIONS-M0.md.
   Tags are exact data, never inferred; numeric smallness never changes a
   tag, a window, or a sector count. *)

BeginPackage["DiffExp2`SectorSeries`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`"}];

ValidateLocalSolution::usage = "ValidateLocalSolution[ls] checks every structural invariant and returns ls unchanged; loud error otherwise.";
CanonicalizeLocalSolution::usage = "CanonicalizeLocalSolution[ls] merges identical and integer-spaced same-(b,p) sectors, drops syntactically zero sectors, sorts by (a,b,p).";
ChartImSign::usage = "ChartImSign[ls] derives the chart Im-sign (+1|-1|None) from the Prescriptions list; conflicting odd-multiplicity entries are a loud error. THE one sign derivation in DiffExp2.";
EvaluateLocalSolution::usage = "EvaluateLocalSolution[ls, tval, opts] evaluates at the chart point tval. Options: \"UsePade\", \"TOrderReduction\", \"ImSign\". Returns <|\"Value\" -> EpsSeries, \"PadeFallbacks\" -> ..., \"TailEstimates\" -> ...|>.";
MultiplyRational::usage = "MultiplyRational[ls, c] multiplies by a rational c(t, eps): center poles shift a, far poles fold into Taylor parts, interior poles are a loud error, eps-denominators shift the windows.";
ReexpandLocalSolution::usage = "ReexpandLocalSolution[ls, newCenter, targetOrder, opts] re-expands around a regular point inside the chart, producing a single-(0,0,0)-sector LocalSolution with the explicit truncation contract.";
DifferentiateLocalSolution::usage = "DifferentiateLocalSolution[ls] gives the exact chart-coordinate derivative via tag algebra.";
SectorDecomposition::usage = "SectorDecomposition[ls] gives the exact sector list + windows (the named replacement for the old DecomposeSingularity fit).";
CombineLocalSolutions::usage = "CombineLocalSolutions[weights, lss] gives the exact linear combination with EpsSeries weights.";
ParseTaggedPower::usage = "ParseTaggedPower[expr, var, epsSym] parses c*var^(a+b*eps)*Log[var]^p into exact tags, or returns $FailedParse.";
$FailedParse::usage = "$FailedParse is the inert marker returned by ParseTaggedPower for non-tagged-power input.";
$ForcePadeFail::usage = "$ForcePadeFail (test hook): force the Pade failure path.";

Begin["`Private`"];

T = Hold;  (* no-op placeholder to avoid bare-context syntax mistakes *)
err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "SectorSeries"|>, payload]];
cfg = DiffExp2`Config`CFG;
numMag = DiffExp2`Tolerances`NumericMagnitude;
numMagBounds = DiffExp2`Tolerances`NumericMagnitudeBounds;
$ForcePadeFail = False;

exactQ[e_] := FreeQ[e, _?InexactNumberQ];
pow[x_, 0] := 1;   (* 0^0 = 1 convention for the empty product *)
pow[x_, n_] := x^n;
(* zeroQ: exact-first (see Solve.m).  Together-canonical rational functions
   over Q(i) in non-numeric symbols: === 0 decides; inexact/radical/
   numeric-constant forms keep the PossibleZeroQ fallback unchanged. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && TrueQ[PossibleZeroQ[c]]);
zeroQ[e_] := zeroCanQ[Together[e]];
structuralZeroQ[e_] := FreeQ[e, _?InexactNumberQ] && zeroQ[e];
badHeadQ[c_] := !FreeQ[c, _SeriesData | _SeriesCoefficient] ||
  MemberQ[{$Failed, Indeterminate, ComplexInfinity}, c];

tagOf[sec_] := {sec["a"], sec["b"], sec["p"]};
lsMin[ls_] := ls["EpsWindow", "Min"];
lsCM[ls_] := ls["EpsWindow", "CompleteMax"];
chartName[ls_] := ToString[ls["Center"], InputForm];

(* ---- 2.1 validation ---- *)

ValidateLocalSolution[ls_Association] := Module[{dims, ncomp},
  If[!AllTrue[{"Center", "Sectors", "EpsWindow", "TWindow", "Radius",
      "ErrorEstimate", "Prescriptions"}, KeyExistsQ[ls, #] &],
    err["dims", <|"Chart" -> chartName[ls], "Detail" -> "missing LocalSolution keys",
      "Keys" -> Keys[ls]|>]];
  If[!IntegerQ[lsMin[ls]] || !IntegerQ[lsCM[ls]] || lsMin[ls] > lsCM[ls],
    err["window", <|"Chart" -> chartName[ls], "EpsWindow" -> ls["EpsWindow"],
      "Detail" -> "bad eps window"|>]];
  dims = Dimensions[#["Coeffs"]] & /@ ls["Sectors"];
  If[Length[DeleteDuplicates[dims]] =!= 1 || Length[First[dims]] =!= 3 ||
      First[First[dims]] =!= lsCM[ls] - lsMin[ls] + 1,
    err["dims", <|"Chart" -> chartName[ls], "Dimensions" -> dims,
      "Expected" -> {lsCM[ls] - lsMin[ls] + 1, "ncols", "Ncomp"},
      "Detail" -> "sector Coeffs dimensions mismatch"|>]];
  Do[
    If[!(exactQ[sec["a"]] && exactQ[sec["b"]] && IntegerQ[sec["p"]] && sec["p"] >= 0),
      err["badtag", <|"Chart" -> chartName[ls], "Sector" -> tagOf[sec],
        "Detail" -> "sector tags must be exact; p a non-negative integer"|>]];
    Module[{bad = Position[sec["Coeffs"], _?badHeadQ, {3}, 1]},
      If[bad =!= {},
        err["badcoeff", <|"Chart" -> chartName[ls], "Sector" -> tagOf[sec],
          "Position" -> First[bad],
          "Detail" -> "compound/invalid head in Coeffs"|>]]],
    {sec, ls["Sectors"]}];
  ls];

(* ---- 2.2 canonicalization ---- *)

shiftCols[arr_, sh_Integer] := Module[{d = Dimensions[arr]},
  (* shift t-columns right by sh (a decreased by sh): row content moves up in n *)
  Map[Join[ConstantArray[0, sh], #][[;; d[[2]]]] &, arr, {1}]];

shiftFitsQ[arr_, sh_Integer] := Module[{n = Dimensions[arr][[2]], tail},
  If[sh <= 0, Return[True, Module]];
  tail = If[sh >= n, arr, arr[[All, n - sh + 1 ;; n, All]]];
  AllTrue[Flatten[tail], # === 0 &]];

CanonicalizeLocalSolution[ls_Association] := Module[{secs, grouped},
  secs = ls["Sectors"];
  (* merge identical tags *)
  grouped = GatherBy[secs, {Together[#["a"] ], Together[#["b"]], #["p"]} &];
  secs = Map[If[Length[#] === 1, First[#],
    Append[First[#], "Coeffs" -> Total[#[[All, "Coeffs"]]]]] &, grouped];
  (* merge same-(b,p) integer-spaced a into minimal a *)
  grouped = GatherBy[secs, {Together[#["b"]], #["p"]} &];
  secs = Flatten[Map[Module[{g = #},
    (* greedy integer-spacing clusters *)
    Module[{rem = SortBy[g, #["a"] &], out = {}, cur, mates, candidates},
      While[rem =!= {},
        cur = First[rem]; rem = Rest[rem];
        candidates = Select[rem,
          IntegerQ[RootReduce[Together[#["a"] - cur["a"]]]] &];
        (* A fixed-width slab cannot absorb a shifted nonzero tail without
           erasing known data.  Keep such a tower separate; downstream
           operations already understand multiple integer-spaced sectors. *)
        mates = Select[candidates, shiftFitsQ[#["Coeffs"],
          RootReduce[Together[#["a"] - cur["a"]]]] &];
        rem = Complement[rem, mates];
        If[mates === {}, AppendTo[out, cur],
          Module[{base = cur["a"], acc = cur["Coeffs"]},
            Do[acc = acc + shiftCols[m["Coeffs"],
                RootReduce[Together[m["a"] - base]]], {m, mates}];
            AppendTo[out, Append[cur, "Coeffs" -> acc]]]]];
      out]] &, grouped], 1];
  (* Drop exact-zero sectors, except an under-expanded negative b=0 tower:
     its finite zero prefix does not certify the endpoint-relevant tail.
     Retaining that tag lets Integrate raise E10 consistently instead of
     making the answer depend on whether canonicalization ran first. *)
  Module[{ncols = Dimensions[First[secs]["Coeffs"]][[2]], kept},
    kept = Select[secs, Function[s, Module[{allZero, needed, preserve},
      allZero = AllTrue[Flatten[s["Coeffs"]], Function[c, c === 0]];
      needed = If[TrueQ[s["a"] < 0], Floor[-s["a"]], -1];
      preserve = allZero && zeroQ[s["b"]] && TrueQ[s["a"] < 0] &&
        (!IntegerQ[needed] || needed >= ncols);
      !allZero || preserve]]];
    If[kept === {}, kept = {First[secs]}];  (* never drop the last sector *)
    secs = kept];
  secs = SortBy[secs, tagOf];
  Append[ls, "Sectors" -> secs]];

(* ---- 2.3 ChartImSign ---- *)

ChartImSign[ls_Association] := Module[{odd, sigmas},
  odd = Select[ls["Prescriptions"], OddQ[#["Multiplicity"]] &];
  If[odd === {}, Return[None]];
  sigmas = DeleteDuplicates[#["Sign"]*#["LeadingCoeffSign"] & /@ odd];
  If[Length[sigmas] > 1,
    err["branchconflict", <|"Chart" -> chartName[ls],
      "Factors" -> (#["Factor"] & /@ odd), "Sigmas" -> sigmas,
      "Detail" -> "conflicting odd-multiplicity prescription signs"|>]];
  First[sigmas]];

(* ---- 2.4 evaluation ---- *)

Options[EvaluateLocalSolution] = {"UsePade" -> Automatic,
  "TOrderReduction" -> 0, "ImSign" -> Automatic};

(* private: diagonal Pade evaluation of one coefficient vector at point *)
padeEvaluate[coeffs_List, point_, ctx_] := Module[{u, mm, poly, pa},
  mm = Floor[Length[coeffs]/2];   (* [m/m] from 2m+1 <= ncoeffs (DEC-12) *)
  poly = Sum[coeffs[[n + 1]]*u^n, {n, 0, Length[coeffs] - 1}];
  pa = If[TrueQ[$ForcePadeFail], PadeApproximant[Unevaluated[Sequence[]]],
    PadeApproximant[poly, {u, 0, {mm, mm}}]];
  If[!FreeQ[pa, _PadeApproximant],
    DiffExp2`Config`PrintWarning["SectorSeries::padefail at ", ctx,
      " - falling back to the direct partial sum (recorded)"];
    {Total[coeffs*Table[pow[point, nn], {nn, 0, Length[coeffs] - 1}]], True},
    {pa /. u -> point, False}]];

EvaluateLocalSolution[ls0_Association, tval_, OptionsPattern[]] := Module[
  {ls = ValidateLocalSolution[ls0], usePade, tred, sigma, wp, secs, ncols,
   ncomp, kmin, kmax, Lv, fallbacks = {}, needsBranch},
  wp = cfg["WorkingPrecision"];
  usePade = OptionValue["UsePade"] /. Automatic :> cfg["UsePade"];
  tred = OptionValue["TOrderReduction"];
  If[!(NumericQ[tval] && Im[tval] == 0),
    err["radius", <|"Chart" -> chartName[ls], "Point" -> tval,
      "Detail" -> "evaluation point must be real numeric (v1)"|>]];
  secs = ls["Sectors"]; kmin = lsMin[ls]; kmax = lsCM[ls];
  ncols = Dimensions[First[secs]["Coeffs"]][[2]] - tred;
  ncomp = Dimensions[First[secs]["Coeffs"]][[3]];
  If[ncols < 1, err["window", <|"Chart" -> chartName[ls],
    "Detail" -> "TOrderReduction removes every column"|>]];
  needsBranch = AnyTrue[secs, !IntegerQ[#["a"]] || !zeroQ[#["b"]] || #["p"] > 0 &];
  Which[
    TrueQ[tval == 0] && needsBranch,
    err["originlimit", <|"Chart" -> chartName[ls],
      "Tags" -> tagOf /@ Select[secs, !IntegerQ[#["a"]] || !zeroQ[#["b"]] || #["p"] > 0 &],
      "Detail" -> "evaluation at the singular center: use Integrate`EndpointSectorLimit via API EndpointLimit"|>],
    TrueQ[Abs[tval] >= ls["Radius"]],
    err["radius", <|"Chart" -> chartName[ls], "Point" -> tval,
      "Radius" -> ls["Radius"], "Detail" -> "point outside the chart radius"|>]];
  sigma = OptionValue["ImSign"] /. Automatic :> ChartImSign[ls];
  If[TrueQ[tval < 0] && needsBranch && !MemberQ[{1, -1}, sigma],
    err["branchmissing", <|"Chart" -> chartName[ls], "Point" -> tval,
      "Detail" -> "negative point on a multivalued chart with no derivable Im-sign"|>]];
  Lv = If[TrueQ[tval < 0], Log[-tval] + I*Pi*sigma, If[TrueQ[tval == 0], 0, Log[tval]]];
  (* value window: each sector is complete to (kmax + p); honest min/max *)
  Module[{valueCM, valueMin2, value, tails},
    valueCM = Min @@ (kmax + #["p"] & /@ secs);
    valueMin2 = valueCM;
    (* first pass: per-sector first nonzero row -> value min *)
    Do[Module[{firstRow},
      firstRow = SelectFirst[Range[kmin, kmax],
        !AllTrue[Flatten[sec["Coeffs"][[# - kmin + 1]]], structuralZeroQ] &,
        kmax];
      valueMin2 = Min[valueMin2, firstRow + sec["p"]]],
      {sec, secs}];
    value = ConstantArray[0, {valueCM - valueMin2 + 1, ncomp}];
    tails = ConstantArray[0, valueCM - valueMin2 + 1];
    Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"],
        ta, alphas},
      ta = If[TrueQ[tval < 0], (-tval)^a*Exp[I*Pi*sigma*a],
        If[TrueQ[tval == 0], If[a === 0, 1, 0], tval^a]];
      Module[{tpows = Table[pow[tval, nn], {nn, 0, ncols - 1}]},
      alphas = If[TrueQ[usePade],
        Table[Module[{slab = Transpose[arr[[k - kmin + 1, 1 ;; ncols]]]},
          Table[Module[{cs = slab[[c2]], v, fb},
            If[!FreeQ[cs, s_Symbol /; Context[s] === "Global`"],
              err["indetpade", <|"Chart" -> chartName[ls], "Sector" -> tagOf[sec],
                "Detail" -> "Pade requested on symbolic coefficients"|>]];
            {v, fb} = padeEvaluate[cs, tval, {chartName[ls], tagOf[sec], k, c2}];
            If[fb, AppendTo[fallbacks,
              <|"Sector" -> tagOf[sec], "EpsRow" -> k, "Component" -> c2|>]];
            v], {c2, ncomp}]],
          {k, kmin, kmax}],
        (* fast path: ONE contraction per sector — the (k x n x c) slab
           against the t-power vector along n; element (k, c2) equals the
           former per-entry arr[[k, 1;;ncols, c2]] . tpows verbatim *)
        Transpose[arr[[All, 1 ;; ncols]], {1, 3, 2}] . tpows]];
      Do[Module[{Kc},
        Kc = Sum[If[kmin <= K - p - j <= kmax,
          alphas[[K - p - j - kmin + 1]]*pow[b*Lv, j]/j!, 0],
          {j, 0, K - p - kmin}];
        If[valueMin2 <= K <= valueCM,
          value[[K - valueMin2 + 1]] += ta*pow[Lv, p]/p!*Kc]],
        {K, kmin + p, valueCM}];
      (* tail estimate: geometric at the actual point *)
      If[ls["Radius"] =!= Infinity && TrueQ[Abs[tval] > 0],
        Do[Module[{nums, topc, q},
          nums = Select[Flatten[arr[[k - kmin + 1, ncols]]], NumericQ];
          topc = If[nums === {}, 0,
            Max[Last[numMagBounds[#, 10]] & /@ nums]];
          If[topc > 0 && kmin <= k + p && k + p <= valueCM && valueMin2 <= k + p,
            q = Abs[N[tval/ls["Radius"], 10]];
            tails[[k + p - valueMin2 + 1]] +=
              topc*Abs[N[tval, 10]]^(ncols - 1)*q/(1 - q)]],
          {k, kmin, kmax}]]],
      {sec, secs}];
    <|"Value" -> DiffExp2`EpsSeries`ESNew[valueMin2, value],
      "PadeFallbacks" -> fallbacks,
      "TailEstimates" -> tails|>]];

(* ---- 2.5 MultiplyRational ---- *)

epsValPoly[p_, eps_] := If[zeroQ[p], Infinity, Exponent[p, eps, Min]];
tValuation[e_, t_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    Exponent[Numerator[c], t, Min] - Exponent[Denominator[c], t, Min]]];

MultiplyRational[ls0_Association, c_, var_Symbol] := Module[
  {ls = ValidateLocalSolution[ls0], eps = DiffExp2`Config`CanonicalEps[],
   cT, num, den, nv, dv, jmin, jcount, cj, d0, troots, wp, M, Q, secs,
   kmin, kmax, ncols, ncomp, newSecs},
  wp = cfg["WorkingPrecision"];
  cT = Together[c /. var -> var];
  If[!PolynomialQ[Numerator[cT], {var, eps}] || !PolynomialQ[Denominator[cT], {var, eps}],
    err["nonrational", <|"Chart" -> chartName[ls], "Expression" -> cT,
      "Detail" -> "multiplier must be rational in (chart variable, eps)"|>]];
  num = Numerator[cT]; den = Denominator[cT];
  If[!exactQ[den],
    err["inexactdenominator", <|"Chart" -> chartName[ls],
      "Detail" -> "denominator must be exact for pole classification"|>]];
  kmin = lsMin[ls]; kmax = lsCM[ls];
  ncols = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2]];
  ncomp = Dimensions[First[ls["Sectors"]]["Coeffs"]][[3]];
  nv = epsValPoly[num, eps]; dv = epsValPoly[den, eps];
  jmin = nv - dv;
  jcount = kmax - kmin + 1;
  (* eps-Laurent coefficients c_j(t), j = jmin .. jmin + jcount - 1 *)
  Module[{Nc, Dc},
    Nc = Table[Coefficient[num, eps, nv + i], {i, 0, jcount - 1}];
    Dc = Table[Coefficient[den, eps, dv + i], {i, 0, jcount - 1}];
    cj = Table[0, {jcount}];
    cj[[1]] = Together[Nc[[1]]/Dc[[1]]];
    Do[cj[[m + 1]] = Together[
        (Nc[[m + 1]] - Sum[Dc[[l + 1]]*cj[[m - l + 1]], {l, 1, m}])/Dc[[1]]],
      {m, 1, jcount - 1}];
    d0 = Dc[[1]]];
  (* classify t-poles of the leading eps denominator *)
  Module[{dred = Cancel[d0/var^Exponent[d0, var, Min]]},
    troots = If[FreeQ[dred, var], {},
      DeleteDuplicates[var /. Solve[dred == 0, var]]];
    Do[Module[{exact = exactQ[r0] && exactQ[ls["Radius"]], absr, rad, diff, interior},
      If[exact,
        (* exact comparison decides; |t_i| == Radius counts as FAR (other chart) *)
        interior = TrueQ[RootReduce[Abs[r0]] < ls["Radius"]] &&
          !TrueQ[PossibleZeroQ[RootReduce[Abs[r0]] - ls["Radius"]]],
        absr = numMag[r0, wp]; rad = N[ls["Radius"], wp];
        diff = absr - rad;
        If[Abs[diff] < DiffExp2`Tolerances`GeomGuardTol[wp]*Max[absr, rad],
          err["geomambiguous", <|"Chart" -> chartName[ls], "Pole" -> r0,
            "Radius" -> ls["Radius"],
            "Detail" -> "pole modulus vs radius numerically ambiguous"|>]];
        interior = diff < 0];
      If[interior && !TrueQ[PossibleZeroQ[r0]],
        err["interiorpole", <|"Chart" -> chartName[ls], "Pole" -> r0,
          "Radius" -> ls["Radius"],
          "Detail" -> "interior pole: segmentation must place a chart here (no radius-shrink fallback)"|>]]],
      {r0, troots}]];
  (* center-pole order and analytic Taylor parts *)
  M = Max[0, Max[Map[Max[0, -tValuation[#, var]] &, cj]]];
  (* one-pass exact Taylor via the division recursion (SeriesCoefficient
     per order is prohibitively slow on big rationals), then numericize at
     2x WP (structure already decided on exact data) *)
  Q = Table[Module[{qj = Cancel[Together[cj[[j]]*var^M]], num1, den1, nc, dc, csr},
      If[zeroCanQ[qj], ConstantArray[0, ncols],
        num1 = Numerator[qj]; den1 = Denominator[qj];
        nc = Table[Coefficient[num1, var, n], {n, 0, ncols - 1}];
        dc = Table[Coefficient[den1, var, n], {n, 0, ncols - 1}];
        csr = ConstantArray[0, ncols];
        csr[[1]] = Together[nc[[1]]/dc[[1]]];
        Do[csr[[m + 1]] = Together[
            (nc[[m + 1]] - Sum[dc[[l + 1]]*csr[[m - l + 1]], {l, 1, m}])/dc[[1]]],
          {m, 1, ncols - 1}];
        (* ByteCount-gated: small exact coefficients stay exact (I-6);
           giants numericize (R6 performance) *)
        Map[If[# === 0 || ByteCount[#] <= 500, #,
          N[#, DiffExp2`Tolerances`$InputPrecisionFactor*wp]] &, csr]]],
    {j, 1, jcount}];
  (* convolve into each sector; a -> a - M; windows shift by jmin *)
  newSecs = Map[Module[{arr = #["Coeffs"], out},
    (* vectorized t-convolution: per (kp, jrow) one ListConvolve per
       component along n (the per-element Sum was the last interpreted
       hot loop - banana: 5 masters x 12 tiles x 1.5M ops each) *)
    out = Table[Module[{acc = ConstantArray[0, {ncols, ncomp}]},
      Do[Module[{krow = kp - (jrow - 1) - jmin},
        If[kmin <= krow <= kmax && !AllTrue[Q[[jrow]], # === 0 &],
          Module[{slabT = Transpose[arr[[krow - kmin + 1]]]},
            acc += Transpose[Table[
              Take[ListConvolve[Q[[jrow]], slabT[[cc]], {1, -1}, 0], ncols],
              {cc, ncomp}]]]]],
        {jrow, 1, jcount}];
      acc], {kp, kmin + jmin, kmax + jmin}];
    <|"a" -> Together[#["a"] - M], "b" -> #["b"], "p" -> #["p"],
      "Coeffs" -> If[FreeQ[out, _Symbol], out, Map[Together, out, {3}]]|>] &,
    ls["Sectors"]];
  CanonicalizeLocalSolution[Join[ls, <|"Sectors" -> newSecs,
    "EpsWindow" -> <|"Min" -> kmin + jmin, "CompleteMax" -> kmax + jmin|>|>]]];

(* ---- 2.6 re-expansion ---- *)

Options[ReexpandLocalSolution] = {"ImSign" -> Automatic};
ReexpandLocalSolution[ls0_Association, Dc_, NN_Integer, OptionsPattern[]] := Module[
  {ls = ValidateLocalSolution[ls0], eps, sigma, kmin, kmax, ncols, ncomp,
   singularQ, rho, LogD, w, out, tails, divOrd, newCenter},
  eps = DiffExp2`Config`CanonicalEps[];
  kmin = lsMin[ls]; kmax = lsCM[ls];
  ncols = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2]];
  ncomp = Dimensions[First[ls["Sectors"]]["Coeffs"]][[3]];
  singularQ = AnyTrue[ls["Sectors"],
    !IntegerQ[#["a"]] || #["a"] < 0 || !zeroQ[#["b"]] || #["p"] > 0 &];
  If[TrueQ[Dc == 0] && singularQ,
    err["resingular", <|"Chart" -> chartName[ls],
      "Detail" -> "re-expansion around the singular center is forbidden; tags come only from fresh Frobenius"|>]];
  If[TrueQ[Abs[Dc] >= ls["Radius"]],
    err["radius", <|"Chart" -> chartName[ls], "Point" -> Dc,
      "Radius" -> ls["Radius"], "Detail" -> "new center outside the chart"|>]];
  sigma = OptionValue["ImSign"] /. Automatic :> ChartImSign[ls];
  If[TrueQ[Dc < 0] && singularQ && !MemberQ[{1, -1}, sigma],
    err["branchmissing", <|"Chart" -> chartName[ls], "Point" -> Dc,
      "Detail" -> "negative re-expansion center on a multivalued chart"|>]];
  LogD = If[TrueQ[Dc < 0], Log[-Dc] + I*Pi*sigma, Log[Dc]];
  rho = Min[ls["Radius"] - Abs[Dc], If[singularQ, Abs[Dc], Infinity]];
  (* assemble per sector via factored (eps-order, w-power) tables:
     value = Da * e^{b eps LogD} * (1+w)^(a+b eps) * (eps(LogD+Log[1+w]))^p/p!
                * Sum_k eps^k Sum_n c[k,n] D^n (1+w)^n,   w = (t-D)/D *)
  Module[{log1pw},
    log1pw = Table[If[m == 0, 0, (-1)^(m + 1)/m], {m, 0, NN}];
    out = ConstantArray[0, {kmax - kmin + 1, NN + 1, ncomp}];
    Do[Module[{a = sec["a"], b = sec["b"], p = sec["p"], arr = sec["Coeffs"],
        Da, P3, G, binTab, Bj2m, j2max},
      Da = If[TrueQ[Dc < 0], (-Dc)^a*Exp[I*Pi*sigma*a], Dc^a];
      (* P3[[m+1]]: w-coeffs of (LogD + Log[1+w])^p / p! *)
      P3 = Module[{acc = ConstantArray[0, NN + 1]},
        acc[[1]] = 1;
        Do[acc = Table[Sum[acc[[m1 + 1]]*
            If[m - m1 == 0, LogD, log1pw[[m - m1 + 1]]], {m1, 0, m}],
          {m, 0, NN}], {p}];
        acc/p!];
      (* G[[k-kmin+1, m+1, cc]]: w-coeffs of the inner t-sum at eps row k *)
      G = Table[Sum[arr[[k - kmin + 1, n + 1, cc]]*Dc^n*Binomial[n, m],
          {n, 0, ncols - 1}],
        {k, kmin, kmax}, {m, 0, NN}, {cc, ncomp}];
      (* Bj2m[[j2+1, m+1]]: eps^j2 coefficient of Binomial[a + b eps, m] *)
      j2max = kmax - kmin;
      binTab = Table[Expand[Binomial[a + b*eps, m]], {m, 0, NN}];
      Bj2m = Table[Coefficient[binTab[[m + 1]], eps, j2],
        {j2, 0, j2max}, {m, 0, NN}];
      Do[
        Module[{contrib},
          contrib = Sum[
            If[kmin <= K - p - j1 - j2 <= kmax,
              pow[b*LogD, j1]/j1!*Sum[
                Bj2m[[j2 + 1, m1 + 1]]*P3[[m2 + 1]]*
                  G[[K - p - j1 - j2 - kmin + 1, m - m1 - m2 + 1, cc]],
                {m1, 0, m}, {m2, 0, m - m1}],
              0],
            {j1, 0, Max[0, K - p - kmin]}, {j2, 0, Min[j2max, Max[0, K - p - kmin]]}];
          out[[K - kmin + 1, m + 1, cc]] += Da*contrib],
        {K, kmin, kmax}, {m, 0, NN}, {cc, ncomp}]],
      {sec, ls["Sectors"]}]];
  (* columns are w-powers (w = t'/D); convert to t'-powers: divide by D^m *)
  out = Table[Together[out[[i, m + 1, cc]]/Dc^m],
    {i, Length[out]}, {m, 0, NN}, {cc, ncomp}];
  divOrd = cfg["DivisionOrder"];
  tails = Table[Module[{topc},
    topc = Max[0, Sequence @@ (Last[numMagBounds[#, 10]] & /@
      Select[Flatten[out[[K - kmin + 1, NN + 1]]], NumericQ])];
    topc*N[(rho/divOrd), 10]^NN/(divOrd - 1)], {K, kmin, kmax}];
  <|"Center" -> ls["Center"] + Dc, "ChartMap" -> ls["ChartMap"],
    "Radius" -> rho,
    "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0, "Coeffs" -> out|>},
    "EpsWindow" -> ls["EpsWindow"],
    "TWindow" -> <|"CompleteMax" -> NN|>,
    "ErrorEstimate" -> If[ListQ[ls["ErrorEstimate"]],
      ls["ErrorEstimate"] + tails, tails],
    "Prescriptions" -> {},
    "TailEstimates" -> tails|>];

(* ---- 2.7 derivative ---- *)

DifferentiateLocalSolution[ls0_Association] := Module[
  {ls = ValidateLocalSolution[ls0], kmin, kmax, ncols, ncomp, newSecs},
  kmin = lsMin[ls]; kmax = lsCM[ls];
  {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
  newSecs = Flatten[Map[Module[{a = #["a"], b = #["b"], p = #["p"], arr = #["Coeffs"], s1, s2},
    s1 = Table[
      (a + n)*arr[[k - kmin + 1, n + 1, cc]] +
        b*If[k - 1 >= kmin, arr[[k - kmin, n + 1, cc]], 0],
      {k, kmin, kmax}, {n, 0, ncols - 1}, {cc, ncomp}];
    s2 = If[p >= 1,
      {<|"a" -> a - 1, "b" -> b, "p" -> p - 1,
        "Coeffs" -> Table[If[k - 1 >= kmin, arr[[k - kmin, n + 1, cc]], 0],
          {k, kmin, kmax}, {n, 0, ncols - 1}, {cc, ncomp}]|>}, {}];
    Join[{<|"a" -> a - 1, "b" -> b, "p" -> p, "Coeffs" -> Map[Together, s1, {3}]|>}, s2]] &,
    ls["Sectors"]], 1];
  CanonicalizeLocalSolution[Join[ls, <|"Sectors" -> newSecs|>]]];

(* ---- 2.8 SectorDecomposition ---- *)

SectorDecomposition[ls0_Association] := Module[{ls = CanonicalizeLocalSolution[
    ValidateLocalSolution[ls0]]},
  <|"Sectors" -> ls["Sectors"], "EpsWindow" -> ls["EpsWindow"],
    "TWindow" -> ls["TWindow"], "Center" -> ls["Center"],
    "ChartMap" -> ls["ChartMap"], "Radius" -> ls["Radius"]|>];

(* ---- 2.10 linear combination with EpsSeries weights ---- *)

CombineLocalSolutions[weights_List, lss_List] := Module[
  {n = Length[lss], base, active, kmin, kmax, tmax, ncols, ncomp,
   secs = {}},
  If[Length[weights] =!= n || n == 0,
    err["dims", <|"Detail" -> "weights/solutions length mismatch"|>]];
  base = ValidateLocalSolution[First[lss]];
  Do[Module[{l = ValidateLocalSolution[lss[[i]]]},
    If[!(l["Center"] === base["Center"] && l["Radius"] === base["Radius"] &&
        l["ChartMap"] === base["ChartMap"] &&
        l["Prescriptions"] === base["Prescriptions"]),
      err["dims", <|"Detail" -> "CombineLocalSolutions requires identical charts",
        "Centers" -> {base["Center"], l["Center"]}|>]];
    If[Dimensions[First[l["Sectors"]]["Coeffs"]][[3]] =!=
        Dimensions[First[base["Sectors"]]["Coeffs"]][[3]],
      err["dims", <|"Detail" -> "component dimensions differ in linear combination"|>]]],
    {i, n}];
  If[AnyTrue[Select[weights, !DiffExp2`EpsSeries`ESQ[#] &],
      Function[w, !AllTrue[DiffExp2`Config`EpsSymbols[], FreeQ[w, #] &]]],
    err["badcoeff", <|"Detail" ->
      "plain linear-combination weights must be epsilon-independent; use EpsSeries for epsilon-dependent weights"|>]];
  tmax = Min[#["TWindow", "CompleteMax"] & /@ lss];
  ncols = tmax + 1;
  ncomp = Dimensions[First[base["Sectors"]]["Coeffs"]][[3]];
  (* A plain scalar is an exact eps-independent multiplier, not the finite
     series ESNew[0,{c}].  Treating it as the latter capped every product at
     lsMin and collapsed multi-sector particular windows.  A plain exact
     zero is inactive; an EpsSeries zero remains active because its missing
     coefficients above CompleteMax are unknown and constrain the sum. *)
  active = Select[Range[n], Function[i,
    If[DiffExp2`EpsSeries`ESQ[weights[[i]]],
      True,
      Together[weights[[i]]] =!= 0]]];
  If[active === {},
    kmax = Min[lsCM /@ lss];
    Return[Join[base, <|
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> ConstantArray[0, {1, ncols, ncomp}]|>},
      "EpsWindow" -> <|"Min" -> kmax, "CompleteMax" -> kmax|>,
      "TWindow" -> <|"CompleteMax" -> tmax|>,
      "ErrorEstimate" -> {0}|>], Module]];
  (* output window: exact scalars preserve their operand window; genuine
     EpsSeries weights use the honest Cauchy-product intersection. *)
  kmin = Min @@ Map[lsMin[lss[[#]]] +
      If[DiffExp2`EpsSeries`ESQ[weights[[#]]],
        DiffExp2`EpsSeries`ESMinPower[weights[[#]]], 0] &, active];
  kmax = Min @@ Map[If[DiffExp2`EpsSeries`ESQ[weights[[#]]],
      Min[lsCM[lss[[#]]] + DiffExp2`EpsSeries`ESMinPower[weights[[#]]],
        DiffExp2`EpsSeries`ESCompleteMax[weights[[#]]] + lsMin[lss[[#]]]],
      lsCM[lss[[#]]]] &, active];
  If[kmax < kmin,
    err["window", <|"Min" -> kmin, "CompleteMax" -> kmax,
      "Detail" -> "linear combination has an empty eps window"|>]];
  Do[Module[{l = lss[[i]], w = weights[[i]], lkmin, lkmax, wmin, wmax, wc},
    lkmin = lsMin[l]; lkmax = lsCM[l];
    If[DiffExp2`EpsSeries`ESQ[w],
      wmin = DiffExp2`EpsSeries`ESMinPower[w];
      wmax = DiffExp2`EpsSeries`ESCompleteMax[w];
      wc = Table[DiffExp2`EpsSeries`ESCoefficient[w, j], {j, wmin, wmax}],
      wmin = 0; wmax = 0; wc = {w}];
    Do[Module[{arr = sec["Coeffs"][[All, 1 ;; ncols]], out},
      (* slab convolution: out[kp] += wc[j]*arr[k], kp = k + j *)
      out = ConstantArray[0, {kmax - kmin + 1, ncols, ncomp}];
      Do[Module[{kp = k + j},
        If[kmin <= kp <= kmax && wc[[j - wmin + 1]] =!= 0,
          out[[kp - kmin + 1]] += wc[[j - wmin + 1]]*arr[[k - lkmin + 1]]]],
        {k, lkmin, lkmax}, {j, wmin, wmax}];
      AppendTo[secs, <|"a" -> sec["a"], "b" -> sec["b"], "p" -> sec["p"],
        "Coeffs" -> If[FreeQ[out, _Symbol], out, Map[Together, out, {3}]]|>]],
      {sec, l["Sectors"]}]],
    {i, active}];
  CanonicalizeLocalSolution[Join[base, <|"Sectors" -> secs,
    "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> kmax|>,
    "TWindow" -> <|"CompleteMax" -> tmax|>,
    "ErrorEstimate" -> Module[{outErr = ConstantArray[0, kmax - kmin + 1]},
      Do[Module[{l, w, lkmin, lkmax, wmin, wmax, wc, ee},
        l = lss[[i]]; w = weights[[i]];
        lkmin = lsMin[l]; lkmax = lsCM[l];
        If[DiffExp2`EpsSeries`ESQ[w],
          wmin = DiffExp2`EpsSeries`ESMinPower[w];
          wmax = DiffExp2`EpsSeries`ESCompleteMax[w];
          wc = Table[DiffExp2`EpsSeries`ESCoefficient[w, j], {j, wmin, wmax}],
          wmin = 0; wmax = 0; wc = {w}];
        ee = If[ListQ[l["ErrorEstimate"]],
          PadRight[l["ErrorEstimate"], lkmax - lkmin + 1],
          ConstantArray[0, lkmax - lkmin + 1]];
        Do[Module[{kp = k + j},
          If[kmin <= kp <= kmax,
            outErr[[kp - kmin + 1]] +=
              If[NumericQ[wc[[j - wmin + 1]]],
                Last[numMagBounds[wc[[j - wmin + 1]], 20]],
                Abs[wc[[j - wmin + 1]]]]*ee[[k - lkmin + 1]]]],
          {k, lkmin, lkmax}, {j, wmin, wmax}]],
        {i, active}];
      outErr]|>]]];

(* ---- 2.11 boundary tagged-power parser ---- *)

ParseTaggedPower[expr0_, var_Symbol, epsSym_Symbol] := Module[
  {expr, lp = 0, rest, expo, coef, a, b, logPoly, logCoef},
  expr = expr0 /. Log[k_*var] :> Log[k] + Log[var];
  rest = expr;
  (* Peel a SINGLE Log[var]^p monomial.  After Log[k var] is expanded,
     mixed log depths (for example Log[2 var]) are a sum of sectors and
     must be handled by the caller's polynomial path, not collapsed here. *)
  logPoly = Expand[rest /. Log[var] -> \[FormalL]];
  lp = Exponent[logPoly, \[FormalL]];
  If[!IntegerQ[lp] || lp < 0, Return[$FailedParse]];
  If[lp > 0,
    logCoef = Coefficient[logPoly, \[FormalL], lp];
    If[!zeroQ[Expand[logPoly - logCoef*\[FormalL]^lp]], Return[$FailedParse]];
    rest = logCoef];
  If[!FreeQ[rest, Log[var]], Return[$FailedParse]];
  (* peel var^expo *)
  Which[
    MatchQ[rest, var^e_ * c_. /; FreeQ[c, var]],
    {expo, coef} = {First[Cases[rest, var^e_ :> e, {0, 2}]],
      rest /. var^_ -> 1},
    MatchQ[rest, var * c_. /; FreeQ[c, var]],
    {expo, coef} = {1, rest/var},
    FreeQ[rest, var],
    {expo, coef} = {0, rest},
    True, Return[$FailedParse]];
  If[!FreeQ[coef, var], Return[$FailedParse]];
  b = Together[D[expo, epsSym]];
  If[!zeroQ[D[b, epsSym]],
    err["badtag", <|"Exponent" -> expo,
      "Detail" -> "non-affine eps exponent in boundary condition"|>]];
  a = Together[expo /. epsSym -> 0];
  <|"a" -> a, "b" -> b, "p" -> lp, "Coefficient" -> coef|>];

End[];
EndPackage[];
