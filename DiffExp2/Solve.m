(* DiffExp2/Solve.m — THE one local solver: symbolic-eps Frobenius with
   EpsSeries coefficient arithmetic.  Spec: Docs/specs/Solve.md (binding);
   decisions: DECISIONS-M0.md.
   True resonance -> explicit log-chain ladder (CASE R); pseudo-resonance ->
   exact Laurent shift (CASE P); regular charts run the SAME recursion
   (CASE T only, V = I).
   CASE-P quotients keep their honest Laurent windows internally, then the
   polar coefficient is cancelled by the already-built target-family
   column (the joint construction of Solve.md 3.7).  The delivered
   multi-sector column is therefore regular as a VALUE while its sector
   representation honestly retains the compensating Laurent pieces. *)

BeginPackage["DiffExp2`Solve`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`"}];

PrepareChart::usage = "PrepareChart[sys, chart] applies the chart map, runs ChartIndicial, and assembles the ChartSystem (theta matrix, gauge, V/VInv spectral frame, families).";
SolveHomogeneous::usage = "SolveHomogeneous[chartSystem, req] gives the FundamentalSystem: one LocalSolution column per indicial sector spec.";
SolveParticular::usage = "SolveParticular[chartSystem, source, req] gives THE particular solution (canonical kernel choice) for a sector-native theta-form source.";
SolveChart::usage = "SolveChart[chartSystem, req, source] gives <|\"Basis\", \"Particular\", \"CouplingDepth\"|>.";
SolveValueRegular::usage = "SolveValueRegular[chartSystem, req, vals] propagates an incoming VALUE vector (one EpsSeries per component: the solution value AT THE CHART CENTER t = 0) through a REGULAR chart with ONE d-dimensional recursion (init = vals); no basis, no matching. The delivered eps-window is capped by the incoming window. Loud error on non-regular charts. (Value-transport prototype; see Docs/PerfGapAnalysis.md lever 1.)";
ClearSolveCaches::usage = "ClearSolveCaches[] empties the PrepareChart and SolveHomogeneous memo caches. Called by API`LoadSystem; the SolveHomogeneous cache additionally self-flushes whenever the chart's SystemHash changes and is entry-capped.";
ODEResidualCheck::usage = "ODEResidualCheck[chartSystem, sol, source, probe] checks the theta-form ODE residual at an interior probe point; loud error above ResidTol.";

Begin["`Private`"];

err[id_, cs_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Solve",
    "Chart" -> ToString[Lookup[cs, "Center", "?"], InputForm]|>, payload]];
cfg = DiffExp2`Config`CFG;
numMag = DiffExp2`Tolerances`NumericMagnitude;
(* zeroQ: exact-first.  Together-canonical rational functions over the
   Gaussian rationals in non-numeric symbols are a zero-DECISION domain:
   === 0 decides.  Only forms outside it (inexact numbers, radicals/Root,
   numeric symbol constants like Pi) fall through to PossibleZeroQ, whose
   numeric ztest meprec-storms and can mis-answer on giant exact rationals.
   zeroCanQ takes input the caller already Together/Cancel-canonicalized. *)
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
esTimes = DiffExp2`EpsSeries`ESTimes; esDiv = DiffExp2`EpsSeries`ESDivide;
esShift = DiffExp2`EpsSeries`ESShift; esCoeff = DiffExp2`EpsSeries`ESCoefficient;
esMin = DiffExp2`EpsSeries`ESMinPower; esCM = DiffExp2`EpsSeries`ESCompleteMax;
esFrom = DiffExp2`EpsSeries`ESFromExpression; esTrim = DiffExp2`EpsSeries`ESTrim;

(* exact t-Laurent coefficient of a rational function at t = 0 (local copy
   of the Indicial recursion; entries exact) *)
(* contract: p is 0 or a polynomial part of a CANCELED fraction, so the
   syntactic test in zeroCanQ is complete *)
polyMinDeg[p_, t_] := If[zeroCanQ[p], Infinity, Exponent[p, t, Min]];
tVal[e_, t_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    polyMinDeg[Numerator[c], t] - polyMinDeg[Denominator[c], t]]];
tLaurent[e_, t_, k_Integer] := Module[{c, v, num, den, nc, dc, ord, csr},
  c = Cancel[Together[e]];
  If[zeroCanQ[c] || k < (v = tVal[c, t]), Return[0]];
  num = Numerator[c]; den = Denominator[c];
  num = Cancel[num/t^polyMinDeg[num, t]]; den = Cancel[den/t^polyMinDeg[den, t]];
  ord = k - v;
  nc = Table[Coefficient[num, t, j], {j, 0, ord}];
  dc = Table[Coefficient[den, t, j], {j, 0, ord}];
  csr = Table[0, {ord + 1}];
  csr[[1]] = Cancel[Together[nc[[1]]/dc[[1]]]];
  Do[csr[[m + 1]] = Cancel[Together[
      (nc[[m + 1]] - Sum[dc[[j + 1]]*csr[[m - j + 1]], {j, 1, m}])/dc[[1]]]],
    {m, 1, ord}];
  csr[[ord + 1]]];

(* ---- PrepareChart ---- *)

$pcCache = <||>;
familyCollisionDepth[roots_List, collisions_List] := Module[{perRoot},
  perRoot = Table[Module[{events, layers},
    (* Homogeneous recursion starts at n=1; same-a n=0 symmetrization is
       diagnostic degeneracy, not an executed division layer. *)
    events = Select[collisions,
      #["Type"] === "LaurentShift" && #["LowerIdx"] === i && #["n"] >= 1 &];
    layers = GatherBy[events, #["n"] &];
    Total[Map[Max[roots[[#[[All, "UpperIdx"]]]][[All, "BlockSize"]]] &, layers]]],
    {i, Length[roots]}];
  If[perRoot === {}, 0, Max[perRoot]]];

PrepareChart[sys_Association, chart_Association] := Module[
  {pcKey = {Hash[sys["Matrix"]], chart["Center"], Lookup[chart, "Scale", 1],
    Lookup[chart, "Radius", None], Lookup[chart, "LocalRadius", None],
    Lookup[chart, "Prescriptions", {}]}},
  If[KeyExistsQ[$pcCache, pcKey], Return[$pcCache[pcKey]]];
  $pcCache[pcKey] = prepareChartCore[sys, chart]];

prepareChartCore[sys_Association, chart_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t, x0, beta, A, Achart, idata, d,
   cols, V, VInv, fams, detV, colIdx},
  t = chart["ChartVar"]; x0 = chart["Center"]; beta = Lookup[chart, "Scale", 1];
  A = sys["Matrix"];
  Achart = Map[Cancel[Together[#]] &,
    beta*(A /. sys["Variable"] -> x0 + beta*t), {2}];
  idata = DiffExp2`Indicial`ChartIndicial[Achart, t, eps,
    <|"Name" -> Lookup[chart, "Name", "chart@" <> ToString[x0, InputForm]],
      "Center" -> x0, "Variable" -> sys["Variable"]|>];
  d = idata["Dimension"];
  cols = {}; fams = {}; colIdx = 0;
  Do[Module[{members, roots = {}, colRange0 = colIdx + 1, collisions = {}},
    members = Reverse[SortBy[fam["Members"], {#["a"], #["b"]} &]];
    Do[
      Do[Module[{chain = mem["Chains"][[ci]]},
        AppendTo[roots, <|"a" -> mem["a"], "b" -> mem["b"],
          "BlockSize" -> Length[chain]|>];
        Do[AppendTo[cols, chain[[q]]]; colIdx++, {q, Length[chain]}]],
        {ci, Length[mem["Chains"]]}],
      {mem, members}];
    Do[Module[{ri = roots[[i]], rj = roots[[j]], nn},
      If[i =!= j,
        nn = RootReduce[Together[rj["a"] - ri["a"]]];
        If[IntegerQ[nn] && nn >= 0 &&
            !(zeroQ[ri["a"] - rj["a"]] && zeroQ[ri["b"] - rj["b"]]),
          AppendTo[collisions, <|"n" -> nn, "LowerIdx" -> i, "UpperIdx" -> j,
            "DeltaB" -> RootReduce[Together[ri["b"] - rj["b"]]],
            "Type" -> If[zeroQ[ri["b"] - rj["b"]], "Log", "LaurentShift"]|>]]]],
      {i, Length[roots]}, {j, Length[roots]}];
    AppendTo[fams, <|"Roots" -> roots,
      "ColumnRange" -> Range[colRange0, colIdx],
      "Collisions" -> collisions,
      "CollisionDepth" -> familyCollisionDepth[roots, collisions]|>]],
    {fam, idata["Families"]}];
  V = Transpose[cols];
  detV = Together[Det[V]];
  If[zeroCanQ[detV], err["E2", chart, <|"Detail" -> "spectral frame V is singular"|>]];
  VInv = Map[Cancel[Together[#]] &, Inverse[V], {2}];
  <|"ChartVar" -> t, "Center" -> x0,
    "ChartMap" -> <|"Center" -> x0, "Scale" -> beta|>,
    (* LocalSolution evaluation is in t, not the physical line variable. *)
    "Radius" -> Lookup[chart, "LocalRadius", chart["Radius"]],
    "SystemHash" -> Hash[sys["Matrix"]],   (* solve-cache flush tag *)
    "Prescriptions" -> Lookup[chart, "Prescriptions", {}],
    "SystemSize" -> d,
    "ThetaMatrix" -> idata["Reduction"]["ThetaMatrix"],
    (* original theta matrix t*A(chart): the residual check runs on the
       GAUGED-BACK solution f = T.g, which solves the ORIGINAL system -
       checking it against the reduced B would fail every gauge chart *)
    "ThetaOriginal" -> Map[Cancel[Together[t*#]] &, Achart, {2}],
    "Gauge" -> idata["Reduction"]["Gauge"],
    "GaugeInverse" -> idata["Reduction"]["GaugeInverse"],
    "Residue" -> idata["Residue"],
    "V" -> V, "VInv" -> VInv, "Families" -> fams,
    "IndicialData" -> idata|>];

(* block list in column order: one entry per Jordan block *)
blockList[cs_] := Flatten[Map[Function[fam,
  Module[{off = First[fam["ColumnRange"]] - 1, bl = {}},
    Do[AppendTo[bl, <|"a" -> r["a"], "b" -> r["b"], "q" -> r["BlockSize"],
        "Cols" -> Range[off + 1, off + r["BlockSize"]]|>];
      off += r["BlockSize"],
      {r, fam["Roots"]}];
    bl]], cs["Families"]], 1];

(* log ceiling for tag (a, b): own p0 + same-b block sizes at higher
   integer-spaced a (the generalized MaxLogPowers formula) *)
logCeiling[cs_, aT_, bT_, p0_, includeEqual_:False] := p0 + Total[Map[
  Module[{dd = RootReduce[Together[#["a"] - aT]]},
    If[zeroQ[#["b"] - bT] && IntegerQ[dd] &&
        (TrueQ[dd > 0] || (includeEqual && dd === 0)), #["q"], 0]] &,
  blockList[cs]]];

(* CASE-P divisions at one Taylor order act on disjoint spectral blocks in
   parallel; only distinct increasing n layers compose.  A size-q Jordan
   target costs q epsilon orders at that layer. *)
pseudoDepthForTag[cs_, aT_, bT_, nStart_Integer] := Module[{events, groups},
  events = Select[blockList[cs], Function[blk, Module[{nn =
      RootReduce[Together[blk["a"] - aT]]},
    IntegerQ[nn] && nn >= nStart && !zeroQ[bT - blk["b"]]]]];
  groups = GatherBy[events,
    RootReduce[Together[#["a"] - aT]] &];
  Total[Map[Max[#[[All, "q"]]] &, groups]]];

(* ---- the recursion ---- *)

(* ===================== PACKED-FRAME RECURSION CORE =====================
   Coefficient vectors are PLAIN LISTS over a fixed eps-frame
   [fb, fb+W-1] (index i <-> eps^(fb+i-1)); products via ListConvolve.
   This bypasses the per-op EpsSeries object layer (Association churn +
   per-coefficient Together) that cost ~18s/chart on 3-master systems.
   Honest windows: the frame is sized by the deterministic work-window
   formula; a per-(n,log,component) integer shadow follows completeness
   through the actual dependency DAG.  EpsSeries objects appear only at
   the boundaries (sources in, LocalSolutions out). *)

(* fast exact eps-expansion of a RATIONAL expr into a frame list;
   ByteCount-gated numericization at 2x WP *)
ratEpsList[expr_, eps_, fb_, W_] := Module[
  {c, num, den, vn, vd, v, rel, out, top},
  out = ConstantArray[0, W];
  c = Cancel[Together[expr]];
  If[zeroCanQ[c], Return[out]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min];
  v = vn - vd;
  If[v < fb, err["E4", <|"Center" -> "(frame)"|>,
    <|"Valuation" -> v, "FrameBase" -> fb,
      "Detail" -> "eps-valuation below the work frame (buffer formula violated)"|>]];
  num = num/eps^vn // Cancel; den = den/eps^vd // Cancel;
  top = fb + W - 1;
  If[v > top, err["E4", <|"Center" -> "(frame)"|>,
    <|"Valuation" -> v, "FrameTop" -> top,
      "Detail" -> "eps-valuation above the work frame (buffer formula violated)"|>]];
  rel = top - v;
  Module[{ncl = Table[Coefficient[num, eps, j], {j, 0, rel}],
    dcl = Table[Coefficient[den, eps, j], {j, 0, rel}], csr},
    csr = ConstantArray[0, rel + 1];
    csr[[1]] = Together[ncl[[1]]/dcl[[1]]];
    Do[csr[[m + 1]] = Together[
        (ncl[[m + 1]] - Sum[dcl[[j + 1]]*csr[[m - j + 1]], {j, 1, m}])/dcl[[1]]],
      {m, 1, rel}];
    csr = Map[preparedEpsCoefficient, csr];
    Do[out[[v - fb + 1 + m]] = csr[[m + 1]], {m, 0, rel}]];
  out];

(* Adaptive singular homogeneous solves use a deliberately narrow lower
   frame first.  A genuine lower underflow throws this private retry signal;
   the same recurrence is then rerun on a wider rectangle.  Outside that
   narrowly scoped Block, the established public E4 contract is unchanged. *)
$adaptiveLowerFrameProbe = False;
$disableAdaptiveLowerFrames = False;   (* private parity/benchmark seam *)
$disableRationalDenominatorFusion = False;  (* private parity/benchmark seam *)
lowerFrameUnderflow[cs_, payload_Association] :=
  If[TrueQ[$adaptiveLowerFrameProbe],
    Throw[Failure["AdaptiveLowerFrame", Join[
      <|"ID" -> "AdaptiveLowerFrame"|>, payload]],
      "DiffExp2AdaptiveLowerFrame"],
    err["E4", cs, payload]];

(* One (d x W) epsilon-frame shift.  This is shared by the ordinary sparse
   polynomial path and the grouped rational path below, so both retain the
   same strict lower-frame contract. *)
shiftFrameBlock[M_, s_Integer, fb_Integer, W_Integer, cs_] := Which[
  s === 0, M,
  s > 0 && s < W,
    Join[ConstantArray[0, {Length[M], s}], M[[All, 1 ;; W - s]], 2],
  s >= W, ConstantArray[0, Dimensions[M]],
  -s < W,
    If[AnyTrue[Flatten[M[[All, 1 ;; -s]]], # =!= 0 &],
      lowerFrameUnderflow[cs, <|"FrameBase" -> fb, "Shift" -> s,
        "Detail" -> "epsilon shift would discard nonzero lower-frame content"|>]];
    Join[M[[All, -s + 1 ;; W]], ConstantArray[0, {Length[M], -s}], 2],
  True,
    If[AnyTrue[Flatten[M], # =!= 0 &],
      lowerFrameUnderflow[cs, <|"FrameBase" -> fb, "Shift" -> s,
        "Detail" -> "epsilon shift lies wholly below the work frame"|>]];
    ConstantArray[0, Dimensions[M]]];

matrixShiftProduct[A_, M_, s_Integer, fb_Integer, W_Integer, cs_] :=
  If[TrueQ[$adaptiveLowerFrameProbe] && s < 0,
    shiftFrameBlock[A . M, s, fb, W, cs],
    A . shiftFrameBlock[M, s, fb, W, cs]];

exactNonRationalNumberQ[e_] := NumericQ[e] && Precision[e] === Infinity &&
  !IntegerQ[e] && Head[e] =!= Rational;

$numericizeAllPreparedNumbers = False;
preparedNumericizationRequiredQ[e_] :=
  exactNonRationalNumberQ[e] ||
  (TrueQ[$numericizeAllPreparedNumbers] && e =!= 0 && NumericQ[e] &&
    Precision[e] === Infinity);

groupedEpsExactSafeQ[e_] := e === 0 ||
  (ByteCount[e] <= 500 && !preparedNumericizationRequiredQ[e]);

preparedEpsCoefficient[e_] := Module[{wp2},
  If[groupedEpsExactSafeQ[e], Return[e]];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*
    DiffExp2`Config`CFG["WorkingPrecision"];
  N[e, wp2]];

(* Canonical epsilon-rational entry.  eps^Valuation is kept outside P/Q;
   Q is normalized to Q(0)=1 so denominator equality is an exact grouping
   key and coefficient division is causal on the finite frame. *)
epsRationalData[expr_, eps_, fb_Integer, top_Integer, cs_] := Module[
  {c, num, den, vn, vd, v, num0, den0, q0, p, q, qd},
  c = Cancel[Together[expr]];
  If[zeroCanQ[c], Return[<|"Zero" -> True, "Valuation" -> Infinity|>]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min];
  v = vn - vd;
  If[v < fb,
    err["E4", cs, <|"Valuation" -> v, "FrameBase" -> fb,
      "Detail" -> "eps-valuation below the work frame (buffer formula violated)"|>]];
  If[v > top,
    err["E4", cs, <|"Valuation" -> v, "FrameTop" -> top,
      "Detail" -> "eps-valuation above the work frame (buffer formula violated)"|>]];
  num0 = Cancel[num/eps^vn]; den0 = Cancel[den/eps^vd];
  q0 = Coefficient[den0, eps, 0];
  If[zeroCanQ[q0],
    err["E5", cs, <|"Denominator" -> den0,
      "Detail" -> "epsilon-unit denominator has zero constant coefficient"|>]];
  p = Expand[Cancel[Together[num0/q0]]];
  q = Expand[Cancel[Together[den0/q0]]];
  If[!PolynomialQ[p, eps] || !PolynomialQ[q, eps],
    err["E5", cs, <|"Expression" -> c,
      "Detail" -> "prepared epsilon-rational entry is not polynomial-over-polynomial"|>]];
  qd = Exponent[q, eps];
  <|"Zero" -> False, "Valuation" -> v, "P" -> p, "Q" -> q,
    "DenominatorDegree" -> qd|>];

(* Q=1 entries need no second Together/series-division pass after
   epsRationalData has already canonicalized them. *)
polynomialEpsListFromData[data_Association, eps_, fb_Integer, W_Integer] :=
 Module[{out = ConstantArray[0, W], v = data["Valuation"], pCoeffs, kmax},
  pCoeffs = CoefficientList[data["P"], eps];
  kmax = Min[Length[pCoeffs] - 1, fb + W - 1 - v];
  Do[out[[v - fb + 1 + k]] = preparedEpsCoefficient[pCoeffs[[k + 1]]],
    {k, 0, kmax}];
  out];

(* Let Mbar be the Laurent series of one rational entry truncated exactly to
   [fb, top], as ratEpsList stores it.  Instead of storing all O(W) entries of
   Mbar, store B = Q Mbar.  B consists of the finite numerator plus at most
   deg(Q) boundary coefficients above top.  Causal division by Q after the
   matrix application therefore reproduces Mbar.U, including the old upper
   cutoff when U has negative epsilon powers. *)
rationalBoundaryData[rec_List, qCoeffs_List, eps_, top_Integer] := Module[
  {v = rec[[4]], p = rec[[5]], qd = Length[qCoeffs] - 1, rel,
   pCoeffs, csr, coeffAt, low, boundary, n0, reason},
  rel = top - v;
  pCoeffs = CoefficientList[p, eps];
  csr = ConstantArray[0, rel + 1];
  Do[
    n0 = If[k + 1 <= Length[pCoeffs], pCoeffs[[k + 1]], 0];
    csr[[k + 1]] = Together[(n0 - Sum[
      qCoeffs[[i + 1]]*csr[[k - i + 1]], {i, 1, Min[k, qd]}])/
        qCoeffs[[1]]],
    {k, 0, rel}];
  coeffAt[k_Integer] := If[0 <= k <= rel, csr[[k + 1]], 0];
  low = Table[{v + k, pCoeffs[[k + 1]]},
    {k, 0, Min[Length[pCoeffs] - 1, rel]}];
  boundary = Table[{top + h, Together[Sum[
      qCoeffs[[i + 1]]*coeffAt[rel + h - i], {i, h, qd}]]},
    {h, 1, qd}];
  reason = Which[
    !AllTrue[pCoeffs, groupedEpsExactSafeQ], "NumeratorByteCount",
    !AllTrue[csr, groupedEpsExactSafeQ], "ImpulseByteCount",
    !AllTrue[boundary[[All, 2]], groupedEpsExactSafeQ], "BoundaryByteCount",
    True, None];
  <|"Pairs" -> Select[Join[low, boundary], #[[2]] =!= 0 &],
    "Safe" -> (reason === None), "RouteReason" -> reason|>];

buildRationalMatrixGroup[records_List, eps_, top_Integer, d_Integer] := Module[
  {q = records[[1, 1]], qCoeffs, boundaryData, bad, decorated,
   shifts, mats, rules},
  qCoeffs = CoefficientList[q, eps];
  If[!groupedEpsExactSafeQ[q] || !AllTrue[qCoeffs, groupedEpsExactSafeQ],
    Return[<|"GroupedSafe" -> False, "Records" -> records,
      "RouteReason" -> "DenominatorByteCount"|>]];
  boundaryData = rationalBoundaryData[#, qCoeffs, eps, top] & /@ records;
  bad = SelectFirst[boundaryData, !TrueQ[#["Safe"]] &, None];
  If[bad =!= None,
    Return[<|"GroupedSafe" -> False, "Records" -> records,
      "RouteReason" -> bad["RouteReason"]|>]];
  decorated = MapThread[{#1[[2]], #1[[3]], #2["Pairs"]} &,
    {records, boundaryData}];
  shifts = Union @@ Map[#[[3, All, 1]] &, decorated];
  mats = Table[
    rules = Map[Function[rr, With[
        {z = SelectFirst[rr[[3]], #[[1]] === s &, None]},
        If[z === None, Nothing, {rr[[1]], rr[[2]]} -> z[[2]]]]],
      decorated];
    {s, SparseArray[rules, {d, d}]},
    {s, shifts}];
  <|"GroupedSafe" -> True, "Denominator" -> q,
    "DenominatorCoefficients" -> qCoeffs,
    "NumeratorSp" -> mats, "EntryCount" -> Length[records]|>];

(* Hybrid preparation: epsilon-Laurent polynomials stay on the established
   sparse-shift path.  Only entries with a genuinely nonconstant epsilon
   denominator are grouped, avoiding a regression on charts such as banana
   x=0 where Nhat is already finite. *)
prepareNhatHybrid[exprs_List, eps_, fb_Integer, W_Integer, cs_] := Module[
  {top = fb + W - 1, nj = Length[exprs], d = Length[First[exprs]],
   polySp = {}, ratGroups = {}, valuations, stats = {}, analyzed,
   polyL, idxs, records, builtGroups, groups, unsafeGroups,
   legacyRecords, legacyIndices},
  valuations = ConstantArray[Infinity, {nj, d, d}];
  Do[
    analyzed = Table[
      epsRationalData[exprs[[j, r, c]], eps, fb, top, cs],
      {r, d}, {c, d}];
    valuations[[j]] = Map[#["Valuation"] &, analyzed, {2}];
    records = Flatten[Table[
      If[TrueQ[analyzed[[r, c, "Zero"]]] ||
          analyzed[[r, c, "DenominatorDegree"]] === 0,
        Nothing,
        {analyzed[[r, c, "Q"]], r, c,
          analyzed[[r, c, "Valuation"]], analyzed[[r, c, "P"]]}],
      {r, d}, {c, d}], 1];
    builtGroups = If[records === {}, {},
      Map[buildRationalMatrixGroup[#, eps, top, d] &,
        GatherBy[records, First]]];
    groups = Select[builtGroups, TrueQ[#["GroupedSafe"]] &];
    unsafeGroups = Select[builtGroups, !TrueQ[#["GroupedSafe"]] &];
    legacyRecords = If[unsafeGroups === {}, {},
      Flatten[unsafeGroups[[All, "Records"]], 1]];
    legacyIndices = Association[Table[
      ((rec[[2]] - 1)*d + rec[[3]]) -> True, {rec, legacyRecords}]];
    polyL = Table[Which[
      TrueQ[analyzed[[r, c, "Zero"]]], ConstantArray[0, W],
      analyzed[[r, c, "DenominatorDegree"]] === 0,
        polynomialEpsListFromData[analyzed[[r, c]], eps, fb, W],
      KeyExistsQ[legacyIndices, (r - 1)*d + c],
        ratEpsList[exprs[[j, r, c]], eps, fb, W],
      True, ConstantArray[0, W]],
      {r, d}, {c, d}];
    idxs = Select[Range[W], Module[{i = #},
      AnyTrue[Flatten[polyL, 1], #[[i]] =!= 0 &]] &];
    AppendTo[polySp,
      Table[{i + fb - 1, Map[#[[i]] &, polyL, {2}]}, {i, idxs}]];
    AppendTo[ratGroups, groups];
    AppendTo[stats, <|
      "PolynomialShifts" -> Length[idxs],
      "LegacySparseShifts" -> Length[idxs],
      "RationalEntries" -> Length[records],
      "GroupedRationalEntries" -> Total[Lookup[groups, "EntryCount", 0]],
      "RationalGroups" -> Length[groups],
      "RationalNumeratorShifts" -> Total[
        Length[#["NumeratorSp"]] & /@ groups],
      "LegacyRationalEntries" -> Length[legacyRecords],
      "LegacyRationalGroups" -> Length[unsafeGroups],
      "LegacyRouteReasons" -> If[unsafeGroups === {}, <||>,
        Counts[unsafeGroups[[All, "RouteReason"]]]],
      "LegacyRationalShiftUpperBound" -> If[records === {}, 0,
        Length[Union @@ (Range[#[[4]], top] & /@ records)]],
      "FrameWidth" -> W|>],
    {j, nj}];
  <|"PolynomialSp" -> polySp, "RationalGroups" -> ratGroups,
    "Valuations" -> valuations, "Stats" -> stats|>];

(* B.U for one grouped rational matrix Mbar, where B = Q Mbar includes the
   finite-top boundary correction.  Keep this operation per (lag, group): a
   lower-frame witness must be raised before any cancellation with another
   lag is allowed. *)
rationalMatrixGroupNumerator[group_Association, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {rhs = ConstantArray[0, Dimensions[U]]},
  (* During an adaptive probe, apply the coefficient matrix before the
     lower-frame shift so banana's rank-one/square-zero annihilation is
     visible before the underflow decision.  Terminal/ordinary runs retain
     the established operation order bit-for-bit. *)
  Do[rhs += matrixShiftProduct[sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, group["NumeratorSp"]}];
  rhs];

(* Causal finite-frame division by the normalized epsilon unit Q.  Q(0) is
   nonzero, so this map is linear and lower triangular. *)
divideRationalMatrixRHS[rhs_List, q_List, W_Integer] := Module[
  {qd = Length[q] - 1, y = ConstantArray[0, Dimensions[rhs]]},
  y[[All, 1]] = rhs[[All, 1]]/q[[1]];
  Do[y[[All, i]] = (rhs[[All, i]] - Sum[
      q[[h + 1]]*y[[All, i - h]], {h, 1, Min[i - 1, qd]}])/q[[1]],
    {i, 2, W}];
  y];

(* Unfused single-lag form retained as the equivalence/benchmark seam. *)
applyRationalMatrixGroups[groups_List, U_List, fb_Integer, W_Integer, cs_] :=
 Module[{acc = ConstantArray[0, Dimensions[U]], rhs},
  Do[
    rhs = rationalMatrixGroupNumerator[group, U, fb, W, cs];
    acc += divideRationalMatrixRHS[
      rhs, group["DenominatorCoefficients"], W],
    {group, groups}];
  acc];

(* White-box benchmark/equivalence seam used by the focused tests. *)
applyPreparedNhat[polySp_List, groups_List, U_List,
    fb_Integer, W_Integer, cs_] := Module[
  {acc = ConstantArray[0, Dimensions[U]]},
  Do[acc += matrixShiftProduct[sp[[2]], U, sp[[1]], fb, W, cs],
    {sp, polySp}];
  If[groups =!= {},
    acc += applyRationalMatrixGroups[groups, U, fb, W, cs]];
  acc];

(* framed product: a, b based at fb -> product re-framed at fb.  A nonzero
   product below the physical frame is a budgeting failure, never a
   sanctioned truncation.  Exact zero and one-term frames occur pervasively
   in the spectral V/VInv and gauge transforms.  Their product is exactly a
   zero or an epsilon shift-and-scale, so do not pay for a dense W-by-W
   ListConvolve at every Taylor row/basis column.  Support classification is
   structural (=!= 0) only; genuinely multi-term frames retain the original
   convolution path unchanged. *)
frConv[a_, b_, fb_, W_] := Module[{sa, sb, ia, ib, lead, shiftMono},
  sa = Select[Range[W], a[[#]] =!= 0 &];
  sb = Select[Range[W], b[[#]] =!= 0 &];
  If[sa === {} || sb === {}, Return[ConstantArray[0, W], Module]];
  ia = First[sa]; ib = First[sb];
  lead = 2 fb + ia + ib - 2;
  If[lead < fb,
    lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
      <|"FrameBase" -> fb, "ProductLead" -> lead,
        "Detail" -> "framed convolution would discard nonzero lower-epsilon content"|>]];
  shiftMono[mono_, pos_, other_] := Module[
    {out = ConstantArray[0, W], shift = fb + pos - 1, dstLo, dstHi},
    dstLo = Max[1, 1 + shift];
    dstHi = Min[W, W + shift];
    If[dstLo <= dstHi,
      out[[dstLo ;; dstHi]] =
        mono*other[[dstLo - shift ;; dstHi - shift]]];
    out];
  If[Length[sa] === 1, Return[shiftMono[a[[ia]], ia, b], Module]];
  If[Length[sb] === 1, Return[shiftMono[b[[ib]], ib, a], Module]];
  Take[ListConvolve[a, b, {1, -1}, 0], {1 - fb, W - fb}]];

(* framed inverse: leading nonzero at frame position pos0 (eps-power
   fb+pos0-1); the inverse leads at the negated power (must lie in-frame) *)
frInv[dlist_, fb_, W_] := Module[{pos0, drel, m, e, out, lead, startIdx},
  pos0 = SelectFirst[Range[W], dlist[[#]] =!= 0 &, None];
  If[pos0 === None,
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "framed inverse of an identically zero series"|>]];
  lead = fb + pos0 - 1;
  startIdx = -lead - fb + 1;   (* frame index of eps^(-lead) *)
  If[startIdx < 1,
    err["E4", <|"Center" -> "(frame)"|>,
      <|"Lead" -> lead, "FrameBase" -> fb,
        "Detail" -> "inverse leads below the work frame"|>]];
  If[startIdx > W,
    err["E4", <|"Center" -> "(frame)"|>,
      <|"Lead" -> lead, "FrameTop" -> fb + W - 1,
        "Detail" -> "inverse leads above the work frame"|>]];
  drel = dlist[[pos0 ;;]];
  m = W - startIdx;            (* relative orders we can store *)
  e = ConstantArray[0, m + 1];
  e[[1]] = If[NumericQ[drel[[1]]], 1/drel[[1]], Together[1/drel[[1]]]];
  Do[e[[k + 1]] = -e[[1]]*Sum[
      If[j + 1 <= Length[drel], drel[[j + 1]], 0]*e[[k - j + 1]], {j, 1, k}],
    {k, 1, m}];
  out = ConstantArray[0, W];
  out[[startIdx ;; startIdx + m]] = e;
  out];

(* x eps and / eps as frame shifts *)
frShiftUp[a_] := Prepend[Most[a], 0];
frDivEps[a_] := Module[{},
  If[First[a] =!= 0,
    lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
      <|"Detail" -> "division by epsilon would discard the lowest framed coefficient"|>]];
  Append[Rest[a], 0]];

epsValuation[e_, eps_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    Exponent[Numerator[c], eps, Min] - Exponent[Denominator[c], eps, Min]]];

(* Symbolic finite-width data.  D is defined only up to multiplication by
   a nonzero function of eps.  Removing its polynomial content in t
   before the d0(eps=0) check: a global 1/eps coefficient is a legitimate
   Laurent multiplier, while a moving factor such as (t-eps) still fails
   E3 because its primitive d0 vanishes at eps=0. *)
clearedSymbolic[cs_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], den,
   denCoeffs, denContent, num, d0, dD, dN, dExpr, NhatExpr},
  den = Together[PolynomialLCM @@ (Denominator[Together[#]] & /@ Flatten[cs["ThetaMatrix"]])];
  denCoeffs = Select[CoefficientList[den, t], !zeroCanQ[#] &];
  If[denCoeffs === {},
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator is identically zero"|>]];
  denContent = If[Length[denCoeffs] === 1, First[denCoeffs],
    Fold[PolynomialGCD, First[denCoeffs], Rest[denCoeffs]]];
  den = Cancel[Together[den/denContent]];
  num = Map[Cancel[Together[#*den]] &, cs["ThetaMatrix"], {2}];
  d0 = Together[den /. t -> 0];
  If[zeroQ[d0 /. eps -> 0],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  dExpr = Table[Cancel[Together[Coefficient[den, t, j]]], {j, 0, dD}];
  NhatExpr = Table[Module[{Nj},
    Nj = Map[Coefficient[#, t, j] &, num, {2}];
    Map[Cancel[Together[#]] &, cs["VInv"] . Nj . cs["V"], {2}]],
    {j, 0, dN}];
  <|"dExpr" -> dExpr, "NhatExpr" -> NhatExpr,
    "dD" -> dD, "dN" -> dN|>];

(* A negative-valued j-th Taylor multiplier can be used again after every
   j Taylor steps.  The old one-hit epsPoleDepth misses this composition.
   This longest-path bound is deterministic and tighter than nmax times the
   worst pole when the pole occurs only at j>1. *)
recurrencePoleDepth[data_Association, nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], dD = data["dD"],
   dN = data["dN"], maxJ, stepCost, depth},
  If[nmax <= 0, Return[0]];
  maxJ = Min[nmax, Max[dD, dN]];
  If[maxJ <= 0, Return[0]];
  stepCost = Table[Module[{vals = {}},
      If[j <= dD, AppendTo[vals, epsValuation[data["dExpr"][[j + 1]], eps]]];
      If[j <= dN, vals = Join[vals,
        Flatten[Map[epsValuation[#, eps] &, data["NhatExpr"][[j + 1]], {2}]]]];
      vals = Select[vals, IntegerQ];
      If[vals === {}, 0, Max[0, -Min[vals]]]],
    {j, 1, maxJ}];
  depth = ConstantArray[0, nmax + 1];
  Do[depth[[n + 1]] = Max[Join[{0}, Table[
      depth[[n - j + 1]] + stepCost[[j]], {j, 1, Min[n, maxJ]}]]],
    {n, 1, nmax}];
  depth[[-1]]];

(* Deepest pole of ONE prepared Taylor multiplier, without composing it
   along the recurrence.  prepareCleared materializes every dExpr/NhatExpr
   lag before runRecursion starts, so the adaptive first rectangle must
   contain even a deep high-j operator (for example t^5/eps^10 when j=1 is
   pole-free).  The longer recurrencePoleDepth remains the terminal bound. *)
recurrenceSingleUsePoleDepth[data_Association] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], vals},
  vals = Join[
    Map[epsValuation[#, eps] &, data["dExpr"]],
    Flatten[Map[epsValuation[#, eps] &, data["NhatExpr"], {3}]]];
  vals = Select[vals, IntegerQ];
  If[vals === {}, 0, Max[0, -Min[vals]]]];

(* Exact epsilon valuation of the t-Laurent coefficients that applyGauge
   will use.  Evaluating T at a generic numeric t is not sufficient: a
   mixed denominator such as eps+t can have finite generic valuation while
   its expansion at the chart origin contains increasingly deep eps poles.
   Keeping the coefficient expressions also lets applyGauge and the work
   budget use precisely the same structural valuation. *)
gaugeCoefficientData[cs_, nmax_Integer] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"],
   T = cs["Gauge"], gv, coeffs, vals},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[<|"TMin" -> 0, "EpsValuation" -> 0, "Coefficients" -> None|>]];
  gv = Min[tVal[#, t] & /@ Flatten[T]];
  (* applyGauge convolves output orders np=0..nmax with
     m=gv..gv+nmax.  Budget and expand exactly that same coefficient range,
     including the allowed (though uncommon) gv>0 case. *)
  coeffs = Table[Map[tLaurent[#, t, m] &, T, {2}],
    {m, gv, gv + nmax}];
  vals = Select[Map[epsValuation[#, eps] &, Flatten[coeffs]], IntegerQ];
  <|"TMin" -> gv,
    "EpsValuation" -> If[vals === {}, 0, Min[0, Min[vals]]],
    "Coefficients" -> coeffs|>];

matrixEpsPoleDepth[m_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], vVals, vMin},
  vVals = Select[Map[epsValuation[#, eps] &, Flatten[m]], IntegerQ];
  vMin = If[vVals === {}, 0, Min[0, Min[vVals]]];
  -vMin];

spectralTransformPoleDepth[cs_] := matrixEpsPoleDepth[cs["V"]];
inverseSpectralTransformPoleDepth[cs_] := matrixEpsPoleDepth[cs["VInv"]];

(* assembleSolution applies V and then T as two finite Cauchy products.
   Their honest-window losses therefore add, even if entries of the formal
   product T.V would cancel.  This is private work headroom only: the
   delivered window is still capped to the caller's request. *)
finalTransformPoleDepth[cs_, nmax_Integer] :=
  spectralTransformPoleDepth[cs] -
    gaugeCoefficientData[cs, nmax]["EpsValuation"];

prepareCleared[cs_, fb_, W_, symbolic_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], data, dD, dN, dL, nhatPrep,
   d0Expr, d0InvScalar, ratGroups, ratDenominators, ratDenominatorKeys,
   numericGeometryQ},
  data = If[symbolic === Automatic, clearedSymbolic[cs], symbolic];
  dD = data["dD"]; dN = data["dN"];
  (* Exact algebraic centers/scales are valuable to the planner, but keeping
     their transformed scalar operators exact through dozens of recurrence
     steps causes severe rational/algebraic expression swell.  On a regular
     chart all structural decisions are already complete, so ground every
     exact numeric operator coefficient once at 2x WP.  Symbolic analytic
     regulators are not NumberQ and remain exact. *)
  numericGeometryQ = TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    AnyTrue[{cs["Center"], cs["ChartMap", "Scale"]},
      exactNonRationalNumberQ];
  Block[{$numericizeAllPreparedNumbers = numericGeometryQ},
    dL = Map[ratEpsList[#, eps, fb, W] &, data["dExpr"]];
    nhatPrep = prepareNhatHybrid[data["NhatExpr"], eps, fb, W, cs]];
  (* Cross-lag fusion keys.  epsRationalData already normalized every Q to
     Q(0)=1, so SameQ of the exact denominator expression is the canonical
     equivalence relation.  Store the small integer index once rather than
     rediscovering it in every Taylor recurrence. *)
  ratGroups = nhatPrep["RationalGroups"];
  ratDenominators = DeleteDuplicatesBy[
    Flatten[ratGroups, 1], #["Denominator"] &];
  ratDenominators = KeyTake[#,
      {"Denominator", "DenominatorCoefficients"}] & /@ ratDenominators;
  ratDenominatorKeys = Lookup[ratDenominators, "Denominator", {}];
  ratGroups = Map[Function[groups, Map[Function[group,
      Append[group, "DenominatorIndex" -> SelectFirst[
        Range[Length[ratDenominatorKeys]],
        ratDenominatorKeys[[#]] === group["Denominator"] &]]], groups]],
    ratGroups];
  (* Sparse-offset forms.  Epsilon-Laurent-polynomial Nhat entries use the
     established shift/matrix representation.  Nonconstant denominators use
     finite grouped numerators plus a causal coefficient recurrence. *)
  (* d0 is an epsilon unit by the cleared-system contract.  In the common
     (and banana) case it is actually epsilon-independent.  Cache the
     reciprocal of the already-expanded coefficient so the recursion keeps
     ratEpsList's exact/2x-WP numericization choice and never constructs a
     dense framed inverse for a scalar.  A genuinely epsilon-dependent d0
     retains the general frInv/frConv path. *)
  d0Expr = Cancel[Together[data["dExpr"][[1]]]];
  d0InvScalar = If[FreeQ[d0Expr, eps],
    Module[{i0 = 1 - fb, c0},
      If[i0 < 1 || i0 > W,
        err["E4", cs, <|"FrameBase" -> fb, "FrameTop" -> fb + W - 1,
          "Detail" -> "epsilon-independent d0 lies outside the work frame"|>]];
      c0 = dL[[1, i0]];
      If[zeroCanQ[c0],
        err["E5", cs, <|"Detail" -> "cleared d0 is zero"|>]];
      If[NumericQ[c0], 1/c0, Together[1/c0]]],
    None];
  Module[{dSp},
    dSp = Table[Module[{lst = dL[[j + 1]]},
      Table[{i + fb - 1, lst[[i]]}, {i, Select[Range[W], lst[[#]] =!= 0 &]}]],
      {j, 0, dD}];
    <|"dL" -> dL, "dD" -> dD, "dN" -> dN,
      "d0InvScalar" -> d0InvScalar,
      "dSp" -> dSp,
      "NhatSp" -> nhatPrep["PolynomialSp"],
      "NhatRationalGroups" -> ratGroups,
      "NhatRationalDenominators" -> ratDenominators,
      "NhatValuations" -> nhatPrep["Valuations"],
      "NhatStats" -> nhatPrep["Stats"]|>]];

(* EpsSeries -> frame list *)
esToFrame[s_, fb_, W_] := Module[{out = ConstantArray[0, W], m1, m2},
  m1 = esMin[s]; m2 = esCM[s];
  Do[If[1 <= k - fb + 1 <= W, out[[k - fb + 1]] = esCoeff[s, k]], {k, m1, m2}];
  out];

(* Completeness shadow helpers.  Infinity means a structurally exact zero,
   which never constrains a sum.  Finite integers are the highest eps power
   known complete for one framed scalar. *)
validShift[Infinity, _Integer, _Integer] := Infinity;
validShift[k_Integer, s_Integer, top_Integer] := Min[top, k + s];
validMin[x_List] := If[x === {}, Infinity, Min @@ x];
frameValuation[x_List, fb_Integer] := Module[{p},
  p = SelectFirst[Range[Length[x]], x[[#]] =!= 0 &, None];
  If[p === None, Infinity, fb + p - 1]];

(* The old validity loop visited every nonzero epsilon coefficient of every
   Nhat entry.  validShift is monotone in the shift, so its minimum is exactly
   the entry's epsilon valuation.  This O(d^2) update is coefficient/window
   identical and avoids restoring an O(W) loop beside the grouped product. *)
updateNhatValidity[accValid_List, valuations_?MatrixQ,
    inputValid_List, frameTop_Integer] := Table[
  Min[accValid[[r]], validMin[Table[
    If[IntegerQ[valuations[[r, c]]],
      validShift[inputValid[[c]], valuations[[r, c]], frameTop],
      Infinity],
    {c, Length[inputValid]}]]],
  {r, Length[accValid]}];

(* Exact O(W) division by the affine spectral offset dA + dB eps.
   CASE T (dA != 0) is a coefficient recurrence; CASE P is an exact
   Laurent shift.  No tolerance-based structural decisions occur here. *)
affineDivFrame[rhs_List, dA_, dB_, fb_Integer, W_Integer] := Module[
  {top = fb + W - 1, invA, invB, y, i},
  If[Length[rhs] =!= W,
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "affine frame divisor received the wrong width"|>]];
  Which[
    !zeroCanQ[dA],
      If[!(fb <= 0 <= top),
        err["E4", <|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "FrameTop" -> top,
            "Detail" -> "CASE-T affine inverse exponent zero is outside the work frame"|>]];
      invA = If[NumericQ[dA], 1/dA, Together[1/dA]];
      If[zeroCanQ[dB],
        invA*rhs,
        y = ConstantArray[0, W];
        y[[1]] = invA*rhs[[1]];
        Do[y[[i]] = invA*(rhs[[i]] - dB*y[[i - 1]]), {i, 2, W}];
        y],
    !zeroCanQ[dB],
      If[!(fb <= -1 && top >= 1),
        err["E4", <|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "FrameTop" -> top,
            "Detail" -> "CASE-P epsilon inverse is outside the work frame"|>]];
      invB = If[NumericQ[dB], 1/dB, Together[1/dB]];
      invB*frDivEps[rhs],
    True,
      err["E5", <|"Center" -> "(frame)"|>,
        <|"Detail" -> "zero affine divisor belongs to CASE R"|>]]];

applyInvD0Frame[v_List, invD0_, d0InvScalar_, fb_Integer, W_Integer] :=
  If[d0InvScalar === None, frConv[v, invD0, fb, W], d0InvScalar*v];

(* Solve d0(eps) ((dA+dB eps) I - N_q) y = rhs.  N_q has a unit
   superdiagonal, so the affine/Jordan solve is bottom-up and d0^-1 is
   applied ONCE per row (not once per Jordan step).  This replaces the old
   dense powers of a truncated framed inverse with q linear recurrences. *)
blockSolveTPFrame[rhs_List, dA_, dB_, invD0_, d0InvScalar_, q_Integer,
    fb_Integer, W_Integer] := Module[{top = fb + W - 1, z, r, j, v},
  If[q < 1 || Length[rhs] =!= q || !AllTrue[rhs, ListQ[#] && Length[#] === W &],
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "malformed affine/Jordan block frame"|>]];
  If[zeroCanQ[dA] && zeroCanQ[dB],
    err["E5", <|"Center" -> "(frame)"|>,
      <|"Detail" -> "zero affine divisor belongs to CASE R"|>]];
  If[zeroCanQ[dA],
    (* Preserve the strict legacy underflow contract before bottom-up sums
       can conceal an individually out-of-frame rhs_j/eps^j term.  The old
       implementation also built inv^q even for an exact-zero RHS. *)
    If[fb > -q || top < 1,
      lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
        <|"FrameBase" -> fb, "FrameTop" -> top, "JordanSize" -> q,
          "Detail" -> "CASE-P Jordan inverse exceeds the lower work frame"|>]];
    Do[
      v = frameValuation[rhs[[j]], fb];
      If[IntegerQ[v] && v - j < fb,
        lowerFrameUnderflow[<|"Center" -> "(frame)"|>,
          <|"FrameBase" -> fb, "RHSComponent" -> j,
            "RHSValuation" -> v, "RequiredShift" -> -j,
            "Detail" -> "CASE-P Jordan solve would discard nonzero lower-frame content"|>]],
      {j, q}]];
  z = Table[ConstantArray[0, W], {q}];
  z[[q]] = affineDivFrame[rhs[[q]], dA, dB, fb, W];
  Do[z[[r]] = affineDivFrame[rhs[[r]] + z[[r + 1]], dA, dB, fb, W],
    {r, q - 1, 1, -1}];
  applyInvD0Frame[#, invD0, d0InvScalar, fb, W] & /@ z];

(* the recursion on frames.  UValid mirrors U without an eps axis:
   UValid[[n+1,l+1,r]] is the honest complete top of that scalar, or
   Infinity for a structural exact zero.  This shadow follows the actual
   dependency DAG; pseudo/log divisions at parallel blocks or log levels
   therefore take a maximum rather than being globally counted. *)
runRecursionFramedSrc[args___] := runRecursion[args];
runRecursion[cs_, prep_, aT_, bT_, P_, nmax_, srcHat_, fb_, W_, init_] := Module[
  {d = cs["SystemSize"], dL = prep["dL"],
   dD = prep["dD"], dN = prep["dN"], invD0, d0InvScalar, blocks, U,
   hits = {}, n0,
   UValid, topValid, frameTop = fb + W - 1, lamL, zeroV, srcData},
  d0InvScalar = Lookup[prep, "d0InvScalar", None];
  invD0 = If[d0InvScalar === None, frInv[dL[[1]], fb, W], None];
  blocks = blockList[cs];
  zeroV = ConstantArray[0, W];
  lamL[m_] := Module[{l = zeroV}, l[[-fb + 1]] = Together[aT + m];
    If[-fb + 2 <= W, l[[-fb + 2]] = Together[bT]]; l];
  srcData[n_, l_] := If[srcHat === None, None,
    Module[{sv = srcHat[n, l]},
      Which[
        sv === None, None,
        AssociationQ[sv] && KeyExistsQ[sv, "Frames"] && KeyExistsQ[sv, "Validity"],
          {sv["Frames"], sv["Validity"]},
        (* already frame lists (vectors of plain lists of length W) *)
        VectorQ[First[sv], AtomQ] || (ListQ[First[sv]] && Length[First[sv]] === W &&
          !AssociationQ[First[sv]]), {sv, ConstantArray[frameTop, d]},
        True, {Map[esToFrame[#, fb, W] &, sv], esCM /@ sv}]]];
  U = Table[Table[Table[zeroV, {d}], {l, 0, P + 1}], {n, 0, nmax}];
  UValid = ConstantArray[Infinity, {nmax + 1, P + 2, d}];
  n0 = If[init === None, 0, 1];
  If[init =!= None,
    Do[
      U[[1, l + 1]] = Map[esToFrame[#, fb, W] &, init[[l + 1]]];
      (* Seed from the input contracts, including finite zero series. *)
      UValid[[1, l + 1]] = esCM /@ init[[l + 1]],
      {l, 0, Min[P, Length[init] - 1]}]];
  (* shift a (d x W) block by s eps-slots (s > 0: content moves UP) *)
  Module[{shB, dSp = prep["dSp"], NhatSp = prep["NhatSp"],
      NhatRationalGroups = prep["NhatRationalGroups"],
      NhatRationalDenominators = prep["NhatRationalDenominators"],
      NhatValuations = prep["NhatValuations"], fuseRational},
  fuseRational = NhatRationalDenominators =!= {} &&
    !TrueQ[$disableRationalDenominatorFusion];
  shB[M_, s_] := shiftFrameBlock[M, s, fb, W, cs];
  Do[Module[{RData, R, RValid, rhsFull, rhsValid},
    RData = Table[Module[{acc = ConstantArray[0, {d, W}],
        accValid = ConstantArray[Infinity, d], rationalRHS, rationalUsed,
        denominatorIndex},
      rationalRHS = If[fuseRational,
        Table[ConstantArray[0, {d, W}], {Length[NhatRationalDenominators]}], {}];
      rationalUsed = If[fuseRational,
        ConstantArray[False, Length[NhatRationalDenominators]], {}];
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]]},
          Do[
            acc += matrixShiftProduct[
              sp[[2]], Ublk, sp[[1]], fb, W, cs],
            {sp, NhatSp[[j + 1]]}];
          If[NhatRationalGroups[[j + 1]] =!= {},
            If[fuseRational,
              Do[
                denominatorIndex = group["DenominatorIndex"];
                (* Accumulate only after the per-contribution numerator
                   application has performed its own lower-frame check. *)
                rationalRHS[[denominatorIndex]] +=
                  rationalMatrixGroupNumerator[group, Ublk, fb, W, cs];
                rationalUsed[[denominatorIndex]] = True,
                {group, NhatRationalGroups[[j + 1]]}],
              acc += applyRationalMatrixGroups[
                NhatRationalGroups[[j + 1]], Ublk, fb, W, cs]]];
          accValid = updateNhatValidity[accValid,
            NhatValuations[[j + 1]],
            UValid[[n - j + 1, l + 1]], frameTop]]],
        {j, 1, Min[n, dN]}];
      If[fuseRational,
        Do[If[TrueQ[rationalUsed[[denominatorIndex]]],
          acc += divideRationalMatrixRHS[
            rationalRHS[[denominatorIndex]],
            NhatRationalDenominators[[denominatorIndex,
              "DenominatorCoefficients"]], W]],
          {denominatorIndex, Length[NhatRationalDenominators]}]];
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]], term, a0, termValid},
          (* lam*U + eps*U_above: lam = (aT + n - j) + bT*eps *)
          a0 = Together[aT + n - j];
          term = a0*Ublk + shB[Together[bT]*Ublk, 1] +
            shB[U[[n - j + 1, l + 2]], 1];
          termValid = Table[validMin[{
              If[zeroCanQ[a0], Infinity, UValid[[n - j + 1, l + 1, r]]],
              If[zeroCanQ[bT], Infinity,
                validShift[UValid[[n - j + 1, l + 1, r]], 1, frameTop]],
              validShift[UValid[[n - j + 1, l + 2, r]], 1, frameTop]}],
            {r, d}];
          Do[
            acc -= sp[[2]]*shB[term, sp[[1]]];
            Do[accValid[[r]] = Min[accValid[[r]],
                validShift[termValid[[r]], sp[[1]], frameTop]], {r, d}],
            {sp, dSp[[j + 1]]}]]],
        {j, 1, Min[n, dD]}];
      If[srcHat =!= None,
        Do[If[n - j >= 0,
          Module[{sd = srcData[n - j, l], sv, svValid},
            If[sd =!= None,
              sv = sd[[1]]; svValid = sd[[2]];
              Do[
                acc += sp[[2]]*shB[sv, sp[[1]]];
                Do[accValid[[r]] = Min[accValid[[r]],
                    validShift[svValid[[r]], sp[[1]], frameTop]], {r, d}],
                {sp, dSp[[j + 1]]}]]]],
          {j, 0, Min[n, dD]}]];
      {acc, accValid}], {l, 0, P}];
    R = RData[[All, 1]];
    RValid = RData[[All, 2]];
    (* solve top-down in l *)
    Do[Module[{},
      rhsFull = R[[l + 1]];
      rhsValid = RValid[[l + 1]];
      Do[
        rhsFull -= sp[[2]]*shB[shB[U[[n + 1, l + 2]], 1], sp[[1]]];
        Do[rhsValid[[r]] = Min[rhsValid[[r]],
            validShift[UValid[[n + 1, l + 2, r]], 1 + sp[[1]], frameTop]],
          {r, d}],
        {sp, dSp[[1]]}];
      rhsFull = Table[rhsFull[[r]], {r, d}];
      Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"],
          dA, dB, solvedValid},
        dA = Together[aT + n - aI]; dB = Together[bT - bI];
        Which[
          !zeroCanQ[dA],
          U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], dA, dB,
              invD0, d0InvScalar, q, fb, W]]];
          solvedValid = Table[validMin[rhsValid[[colsB[[r ;; q]]]]], {r, q}];
          UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
            Thread[colsB -> solvedValid]],
          !zeroCanQ[dB],
          (U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], dA, dB,
              invD0, d0InvScalar, q, fb, W]]];
          solvedValid = Table[validMin[Table[
              validShift[rhsValid[[colsB[[r + m]]]], -(m + 1), frameTop],
              {m, 0, q - r}]], {r, q}];
          UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
            Thread[colsB -> solvedValid]];
          If[l == 0, AppendTo[hits, <|"n" -> n, "Cols" -> colsB,
            "DeltaB" -> dB, "FrameBase" -> fb,
            "GammaFrames" -> U[[n + 1, l + 1, colsB]],
            "GammaValidity" -> UValid[[n + 1, l + 1, colsB]]|>]]),
          True, Null]],
        {blk, blocks}]],
      {l, P, 0, -1}];
    (* CASE R ladder *)
    Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"]},
      If[zeroQ[aT + n - aI] && zeroQ[bT - bI],
        Module[{Rt, RtValid, assigned = Association[], assignedValid = Association[]},
          Rt = Table[Table[
            applyInvD0Frame[R[[l + 1, colsB[[r]]]], invD0, d0InvScalar,
              fb, W], {r, q}], {l, 0, P}];
          RtValid = RValid[[All, colsB]];
          Do[If[l + 1 <= P,
            assigned[{l + 1, q}] = frDivEps[Rt[[l + 1, q]]];
            assignedValid[{l + 1, q}] =
              validShift[RtValid[[l + 1, q]], -1, frameTop]],
            {l, 0, P - 1}];
          Do[
            Do[Module[{upr = Lookup[assigned, Key[{l + 1, r}], zeroV], uprValid},
              uprValid = Lookup[assignedValid, Key[{l + 1, r}], Infinity];
              If[!KeyExistsQ[assigned, {l, r + 1}],
                assigned[{l, r + 1}] = frShiftUp[upr] - Rt[[l + 1, r]];
                assignedValid[{l, r + 1}] = Min[
                  validShift[uprValid, 1, frameTop], RtValid[[l + 1, r]]]]],
              {r, 1, q - 1}],
            {l, P, 0, -1}];
          Do[If[KeyExistsQ[assigned, {l, r}],
            U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
              colsB[[r]] -> assigned[{l, r}]];
            UValid[[n + 1, l + 1]] = ReplacePart[UValid[[n + 1, l + 1]],
              colsB[[r]] -> assignedValid[{l, r}]]],
            {l, 0, P}, {r, 1, q}]]]],
      {blk, blocks}]],
    {n, n0, nmax}]];
  topValid = Module[{f = Select[Flatten[UValid], IntegerQ]},
    If[f === {}, frameTop, Min @@ f]];
  <|"U" -> U, "Validity" -> UValid, "Hits" -> hits, "P" -> P,
    "FrameBase" -> fb, "TopValid" -> topValid|>];

(* ---- assembly ---- *)

(* frame U -> original-frame LocalSolution with sectors (a, b, l) *)
assembleSolution[cs_, aT_, bT_, rec_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], d = cs["SystemSize"], U = rec["U"],
   UValid = rec["Validity"], P = rec["P"], fb = rec["FrameBase"], W,
   frameTop, kminO, kmaxO, VexpL, VVal, WL, WValid, secs, ls, firstNZ},
  W = Length[U[[1, 1, 1]]];
  frameTop = fb + W - 1;
  VexpL = Map[ratEpsList[Together[#], eps, fb, W] &, cs["V"], {2}];
  VVal = Map[frameValuation[#, fb] &, VexpL, {2}];
  WL = Table[
    Table[Sum[frConv[VexpL[[r, c]], U[[n + 1, l + 1, c]], fb, W], {c, d}], {r, d}],
    {n, 0, nmax}, {l, 0, P}];
  WValid = Table[
    Table[validMin[Table[
      If[VVal[[r, c]] === Infinity, Infinity,
        validShift[UValid[[n + 1, l + 1, c]], VVal[[r, c]], frameTop]],
      {c, d}]], {r, d}],
    {n, 0, nmax}, {l, 0, P}];
  kmaxO = Module[{f = Select[Flatten[WValid], IntegerQ]},
    If[f === {}, frameTop, Min @@ f]];
  (* content min: an all-zero result has no first coefficient and is
     represented canonically at its honest complete top. *)
  firstNZ = SelectFirst[Range[W],
    Function[i, AnyTrue[Flatten[WL, 2], Function[fr, fr[[i]] =!= 0]]], None];
  kminO = If[firstNZ === None, kmaxO, fb + firstNZ - 1];
  If[firstNZ =!= None && kmaxO < kminO,
    err["E6", cs, <|"FrameBase" -> fb, "FrameTop" -> frameTop,
      "ContentMin" -> kminO, "CompleteMax" -> kmaxO,
      "Tag" -> {aT, bT},
      "Detail" -> "recursion completeness is exhausted below the first stored coefficient"|>]];
  secs = If[firstNZ === None,
    Table[<|"a" -> aT, "b" -> bT, "p" -> l,
      "Coeffs" -> ConstantArray[0, {1, nmax + 1, d}]|>, {l, 0, P}],
    Table[<|"a" -> aT, "b" -> bT, "p" -> l,
      "Coeffs" -> Table[Table[Table[WL[[n + 1, l + 1, r, k - fb + 1]], {r, d}],
        {n, 0, nmax}], {k, kminO, kmaxO}]|>, {l, 0, P}]];
  ls = <|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
    "Radius" -> cs["Radius"], "Sectors" -> secs,
    "EpsWindow" -> <|"Min" -> kminO, "CompleteMax" -> kmaxO|>,
    "TWindow" -> <|"CompleteMax" -> nmax|>,
    "ErrorEstimate" -> ConstantArray[0, kmaxO - kminO + 1],
    "Prescriptions" -> cs["Prescriptions"]|>;
  ls = applyGauge[cs, ls, nmax];
  DiffExp2`SectorSeries`CanonicalizeLocalSolution[ls]];

(* gauge composition f = T.g: T exact rational, poles only at t = 0.
   Frame-list implementation: the k-axis (eps) stays a plain index, the
   t-convolution runs per (np, r) with ListConvolve on eps-frame lists.
   A negative eps-lead in T pays the window top (same rule as the V
   multiply in assembleSolution). *)
applyGauge[cs_, ls_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], T = cs["Gauge"], d,
   gv, kmin, kmax, newSecs, vT, fbT, WG, TexpL, gaugeData, Tcoeffs},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[ls]];
  d = Length[T];
  gaugeData = gaugeCoefficientData[cs, nmax];
  gv = gaugeData["TMin"];
  Tcoeffs = gaugeData["Coefficients"];
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  vT = gaugeData["EpsValuation"];
  fbT = vT; WG = (kmax - kmin) + (-vT) + 1;
  (* TexpL[m - gv + 1][r][c] = eps-frame list (base fbT, width WG) of the
     t^m Laurent coefficient of T_rc *)
  TexpL = Map[
    If[zeroCanQ[#], ConstantArray[0, WG], ratEpsList[#, eps, fbT, WG]] &,
    Tcoeffs, {3}];
  newSecs = Map[Module[{arr = #["Coeffs"], aS = #["a"], outF, kminN, kmaxN,
      ncolsS = Dimensions[#["Coeffs"]][[2]], topValidG, firstNZ},
    (* arr columns as eps-frame lists on [kmin, kmin + WG - 1] (pad top) *)
    outF = Table[ConstantArray[0, WG], {nmax + 1}, {d}];
    Do[Module[{n = np - (m - gv)},
      If[0 <= n <= nmax && n < ncolsS,
        Do[Module[{arrCol = PadRight[arr[[All, n + 1, c]], WG], tl},
          Do[
            tl = TexpL[[m - gv + 1, r, c]];
            If[!AllTrue[tl, # === 0 &],
              outF[[np + 1, r]] += Take[
                ListConvolve[tl, arrCol, {1, -1}, 0], WG]],
            {r, d}]],
          {c, d}]]],
      {np, 0, nmax}, {m, gv, gv + nmax}];
    topValidG = kmax + vT;  (* pay |vT| at the top (vT <= 0) *)
    firstNZ = SelectFirst[Range[WG],
      AnyTrue[Flatten[outF[[All, All, #]]], Function[z, z =!= 0]] &, None];
    kminN = If[firstNZ === None, topValidG, kmin + fbT + firstNZ - 1];
    If[firstNZ =!= None && topValidG < kminN,
      err["E6", cs, <|"InputWindow" -> {kmin, kmax},
        "GaugeValuation" -> vT, "ContentMin" -> kminN,
        "CompleteMax" -> topValidG, "Tag" -> {aS, #["b"], #["p"]},
        "Detail" -> "gauge multiplication exhausted completeness below the first stored coefficient"|>]];
    kmaxN = topValidG;
    <|"a" -> Together[aS + gv], "b" -> #["b"], "p" -> #["p"],
      "Coeffs" -> Table[Table[Table[
        Module[{idx = k - (kmin + fbT) + 1},
          If[1 <= idx <= WG, outF[[n + 1, r, idx]], 0]], {r, d}],
        {n, 0, nmax}], {k, kminN, kmaxN}],
      "KMin" -> kminN, "KMax" -> kmaxN|>] &, ls["Sectors"]];
  Module[{kminA = Min[#["KMin"] & /@ newSecs], kmaxA = Min[#["KMax"] & /@ newSecs]},
    newSecs = Map[Module[{sh = #["KMin"] - kminA},
      <|"a" -> #["a"], "b" -> #["b"], "p" -> #["p"],
        "Coeffs" -> Table[
          If[kminA + i - 1 < #["KMin"],
            Map[0 &, #["Coeffs"][[1]], {2}],
            #["Coeffs"][[i - sh]]],
          {i, 1, kmaxA - kminA + 1}]|>] &, newSecs];
    Join[ls, <|"Sectors" -> newSecs,
      "EpsWindow" -> <|"Min" -> kminA, "CompleteMax" -> kmaxA|>,
      "ErrorEstimate" -> ConstantArray[0, kmaxA - kminA + 1]|>]]];

(* Joint CASE-P compensation (Solve.md 3.7).  At a collision the raw
   quotient coefficient gamma along target block r multiplies the target
   root's leading chain vector.  Its negative-epsilon polar part is not a
   new solution: it must be cancelled by -gammaPolar times the already
   built target-family column.  Zeros above eps^-1 in that polynomial are
   structural and may be certified through the product's honest top.

   Target columns are kept in their wider work windows until every source
   column has been compensated.  A depth-one collision consumes one upper
   order in the target product; parallel collisions consume it only once.
   The deterministic CollisionDepth work budget guarantees the requested
   top survives. *)
compensatePseudoColumn[cs_, raw_Association, hits_List, workColumns_List,
    reqMax_Integer] := Module[{out = raw, records = {}},
  Do[Module[{fb = hit["FrameBase"], gammas = hit["GammaFrames"], cols = hit["Cols"]},
    Do[Module[{gamma = gammas[[r]], negIdx, targetIndex, target, wMin,
        targetMin, targetMax, desiredMax, wMax, weight},
      negIdx = Select[Range[Length[gamma]],
        fb + # - 1 < 0 && gamma[[#]] =!= 0 &];
      If[negIdx =!= {},
        targetIndex = cols[[r]];
        If[targetIndex > Length[workColumns],
          err["E5", cs, <|"TargetColumn" -> targetIndex,
            "BuiltColumns" -> Length[workColumns], "Hit" -> KeyDrop[hit,
              {"GammaFrames", "GammaValidity"}],
            "Detail" -> "pseudo-resonant compensation target is not yet built (family ordering violated)"|>]];
        target = workColumns[[targetIndex]];
        wMin = fb + First[negIdx] - 1;
        targetMin = target["EpsWindow", "Min"];
        targetMax = target["EpsWindow", "CompleteMax"];
        desiredMax = Min[out["EpsWindow", "CompleteMax"], targetMax + wMin];
        If[desiredMax < reqMax,
          err["E4", cs, <|"RequestedCompleteMax" -> reqMax,
            "AvailableCompleteMax" -> desiredMax, "WeightMin" -> wMin,
            "TargetWindow" -> target["EpsWindow"],
            "Detail" -> "pseudo-resonant compensation exhausted the deterministic work window"|>]];
        wMax = desiredMax - targetMin;
        weight = esNew[wMin, Table[
          If[k < 0 && 1 <= k - fb + 1 <= Length[gamma],
            -gamma[[k - fb + 1]], 0], {k, wMin, wMax}]];
        out = DiffExp2`SectorSeries`CombineLocalSolutions[
          {1, weight}, {out, target}];
        AppendTo[records, <|"n" -> hit["n"], "TargetColumn" -> targetIndex,
          "DeltaB" -> hit["DeltaB"], "WeightWindow" ->
            <|"Min" -> wMin, "CompleteMax" -> wMax|>|>]]],
      {r, Length[cols]}]],
    {hit, hits}];
  {out, records}];

pseudoHitPolarOrders[hit_Association] := If[KeyExistsQ[hit, "PolarOrders"],
  hit["PolarOrders"],
  Module[{fb = hit["FrameBase"], gammas = hit["GammaFrames"],
      cols = hit["Cols"]},
    DeleteDuplicates[Flatten[Table[
      Map[{cols[[r]], fb + # - 1} &,
        Select[Range[Length[gammas[[r]]]],
          fb + # - 1 < 0 && gammas[[r, #]] =!= 0 &]],
      {r, Length[cols]}], 1]]]];

publicPseudoHit[hit_Association] := Append[
  KeyDrop[hit, {"GammaFrames", "GammaValidity", "FrameBase"}],
  "PolarOrders" -> pseudoHitPolarOrders[hit]];

(* Spectral-coordinate probe retained as a diagnostic/test seam.  It is NOT
   the CASE-P certificate: when VInv is epsilon-singular, a finite physical
   solution can have legitimate polar spectral coordinates. *)
spectralProbeValue[cs_, ls_, t0_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], value, fcomps,
   H, hTop, out},
  value = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
    "UsePade" -> False]["Value"];
  fcomps = Table[esNew[esMin[value],
      Table[esCoeff[value, k][[c]], {k, esMin[value], esCM[value]}]],
    {c, cs["SystemSize"]}];
  H = Map[Cancel[Together[#]] &,
    cs["VInv"] . (cs["GaugeInverse"] /. t -> t0), {2}];
  hTop = Max[0, -1 - Min[esMin /@ fcomps]];
  out = Table[Module[{acc = None},
      Do[If[!zeroCanQ[H[[r, c]]],
        Module[{term = esTimes[esFrom[H[[r, c]], eps, hTop], fcomps[[c]]]},
          acc = If[acc === None, term, esAdd[acc, term]]]],
        {c, cs["SystemSize"]}];
      If[acc === None, esZero[esCM[value]], acc]],
    {r, cs["SystemSize"]}];
  out];

(* Gauge-reduced physical value g = GaugeInverse.f.  Unlike the individual
   coordinates of u = VInv.g, this frame does not acquire artificial Laurent
   poles from an epsilon-singular spectral inverse. *)
reducedProbeValue[cs_, ls_, t0_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], value, fcomps,
   H, hTop, out},
  value = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
    "UsePade" -> False]["Value"];
  fcomps = Table[esNew[esMin[value],
      Table[esCoeff[value, k][[c]], {k, esMin[value], esCM[value]}]],
    {c, cs["SystemSize"]}];
  H = Map[Cancel[Together[#]] &,
    cs["GaugeInverse"] /. t -> t0, {2}];
  hTop = Max[0, -1 - Min[esMin /@ fcomps]];
  out = Table[Module[{acc = None},
      Do[If[!zeroCanQ[H[[r, c]]],
        Module[{term = esTimes[esFrom[H[[r, c]], eps, hTop], fcomps[[c]]]},
          acc = If[acc === None, term, esAdd[acc, term]]]],
        {c, cs["SystemSize"]}];
      If[acc === None, esZero[esCM[value]], acc]],
    {r, cs["SystemSize"]}];
  out];

(* Map every registered spectral collision order through the exact V column.
   The certificate checks the corresponding polar orders in the reduced
   physical components.  Requiring the VInv coordinates themselves to be
   pole-free is wrong when VInv has epsilon poles: the banana L2 x=0 frame is
   the minimal counterexample. *)
reducedPolarOrders[cs_, hits_List] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], V = cs["V"], pairs},
  pairs = DeleteDuplicates[Flatten[pseudoHitPolarOrders /@ hits, 1]];
  DeleteDuplicates[Flatten[Table[Module[
      {target = pair[[1]], k = pair[[2]], v},
      Flatten[Table[
        v = epsValuation[V[[c, target]], eps];
        If[IntegerQ[v] && k + v < 0,
          Table[{c, ko, target}, {ko, k + v, -1}], {}],
        {c, cs["SystemSize"]}], 1]],
    {pair, pairs}], 1]]];

certifyPseudoCompensation[cs_, ls_, hits_List, label_] := Module[
  {polarPairs, probes, nmax, minA, tailOrder, fracRaw, frac, wp2,
   ltol = DiffExp2`Tolerances`Tol["LaurentLeadTol"]},
  If[hits === {}, Return[True]];
  polarPairs = reducedPolarOrders[cs, hits];
  If[polarPairs === {}, Return[True]];
  (* A compensated target whose integer a is shifted relative to the raw
     sector has one unmatched top Taylor tail at finite TOrder.  Probe deeply
     enough inside the disk that this known truncation tail lies below the
     structural floor; the two probes remain deterministic and nonzero. *)
  nmax = ls["TWindow", "CompleteMax"];
  minA = Min[Floor[#["a"]] & /@ ls["Sectors"]];
  tailOrder = Max[1, nmax + minA];
  fracRaw = N[(ltol/100)^(1/tailOrder), 40];
  frac = Min[1/5, Rationalize[fracRaw, Abs[fracRaw]/50]];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"];
  probes = N[{cs["Radius"]*frac/2, 2*cs["Radius"]*frac/3}, wp2];
  Do[Module[{gv = reducedProbeValue[cs, ls, t0]},
    Do[Module[{r = pair[[1]], k = pair[[2]], target = pair[[3]], allVals,
        scale, z, mag, uncertainty, acc},
      If[k > esCM[gv[[r]]],
        (* The source may not certify this higher polar order.  Its honest
           CompleteMax is retained, so skip rather than infer a value. *)
        Continue[]];
      allVals = Table[If[k < esMin[gv[[c]]], 0, esCoeff[gv[[c]], k]],
        {c, Length[gv]}];
      scale = Max[1, Sequence @@
        (numMag[#, 20] & /@ Select[Flatten[allVals], NumericQ])];
      z = If[k < esMin[gv[[r]]], 0, esCoeff[gv[[r]], k]];
      uncertainty = If[InexactNumberQ[z],
        acc = Accuracy[z];
        If[NumericQ[acc] && acc =!= Infinity, 10^-Floor[acc], 0], 0];
      Which[
        NumericQ[z],
          mag = numMag[z, 20];
          If[TrueQ[mag + uncertainty > ltol*scale],
            err["E5", cs, <|"Certificate" -> label, "Probe" -> t0,
              "EpsOrder" -> k, "TargetColumn" -> target,
              "ReducedComponent" -> r,
              "Residual" -> z, "ResidualUncertainty" -> uncertainty,
              "Scale" -> scale, "Tolerance" -> ltol,
              "Detail" -> "joint pseudo-resonant polar value did not cancel"|>]],
        !zeroQ[z],
          err["E5", cs, <|"Certificate" -> label, "Probe" -> t0,
            "EpsOrder" -> k, "TargetColumn" -> target,
            "ReducedComponent" -> r, "Residual" -> z,
            "Detail" -> "joint pseudo-resonant symbolic polar value did not cancel"|>],
        True, Null]],
      {pair, polarPairs}]],
    {t0, probes}];
  True];

(* ---- public functions ---- *)

(* SolveHomogeneous memo: the fundamental basis depends ONLY on the
   ChartSystem and the request window — solveHomogeneousCore reads cs
   (blocks/families, the prepared cleared system, V/gauge) and req (TOrder,
   EpsWindow); boundary values never enter (they meet the basis later in
   MatchWeights).  The lo/hi endpoint transports of one level share the
   anchor chart (and any coinciding centers), so the second direction's
   solve is a pure replay: memoize on (Hash[cs], request window, WP).
   PrepareChart already dedups cs itself ($pcCache), so shared charts hash
   identically.  WorkingPrecision is in the key because ratEpsList
   numericizes large coefficients at 2x WP.  Memory policy: the cache holds
   ONE system's charts (flushed when SystemHash changes — each FT level is
   a fresh system), is entry-capped as a runaway guard, and is cleared by
   ClearSolveCaches[] from API`LoadSystem.  A memo hit skips the
   ODEResidualCheck rerun: the identical result already passed it when
   computed. *)
$shCache = <||>; $shSysTag = None; $shCacheMax = 64;
SolveHomogeneous[cs_Association, req_Association] := Module[
  {tag = Lookup[cs, "SystemHash", None], key},
  key = {Hash[cs], req["TOrder"], req["EpsWindow", "Min"],
    req["EpsWindow", "CompleteMax"], cfg["WorkingPrecision"],
    TrueQ[$disableAdaptiveLowerFrames],
    TrueQ[$disableRationalDenominatorFusion]};
  If[tag =!= $shSysTag, $shCache = <||>; $shSysTag = tag];
  If[KeyExistsQ[$shCache, key], Return[$shCache[key]]];
  If[Length[$shCache] >= $shCacheMax, $shCache = <||>];
  $shCache[key] = solveHomogeneousCore[cs, req]];

ClearSolveCaches[] := ($pcCache = <||>; $shCache = <||>; $shSysTag = None;);

solveHomogeneousCore[cs_Association, req_Association] := Module[
  {d = cs["SystemSize"], blocks = blockList[cs], nmax, reqMin, reqMax,
   columns = {}, workColumns = {}, specs = {}, hitsAll = {}, compAll = {},
   fams = cs["Families"], colCursor = 0, wideTop, Pmax, cdMax,
   symbolic, poleDepth, singleUseDepth, spectralDepth, transformDepth,
   startFb, terminalFb, adaptiveQ, prepCache = <||>, prepFor,
   adaptiveDiags = {}, certs = {}},
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"]; reqMax = req["EpsWindow", "CompleteMax"];
  Pmax = Max[0, Max[Table[logCeiling[cs, b["a"], b["b"], b["q"] - 1], {b, blocks}]]];
  cdMax = Max[0, Max[#["CollisionDepth"] & /@ fams]];
  symbolic = clearedSymbolic[cs];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  singleUseDepth = recurrenceSingleUsePoleDepth[symbolic];
  spectralDepth = spectralTransformPoleDepth[cs];
  transformDepth = finalTransformPoleDepth[cs, nmax];
  wideTop = reqMax + Pmax + cdMax + 2 - Min[0, reqMin - Pmax - 2];
  (* Keep the established scratch halo, enlarging it only when the exact
     longest recurrence path composes more negative-eps shifts. *)
  wideTop = Max[wideTop, reqMax + Pmax + cdMax + poleDepth];
  wideTop += transformDepth;
  terminalFb = Min[Min[reqMin, 0] - Pmax - cdMax - 2,
      reqMin - Pmax - cdMax - poleDepth] - spectralDepth;
  startFb = Min[Min[reqMin, 0] - Pmax - cdMax - 2,
      reqMin - Pmax - cdMax - singleUseDepth] - spectralDepth;
  (* The upper frame remains the proven scalar longest-path budget, so
     CompleteMax/UValid are unchanged.  Only singular homogeneous lower
     support is adaptive.  The terminal retry is the former rectangle and
     uses the same strict recurrence -- never an alternate solver. *)
  adaptiveQ = !TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
    startFb > terminalFb && !TrueQ[$disableAdaptiveLowerFrames];
  If[!adaptiveQ, startFb = terminalFb];
  prepFor[ff_Integer] := If[KeyExistsQ[prepCache, ff], prepCache[ff],
    AssociateTo[prepCache, ff -> prepareCleared[
      cs, ff, wideTop - ff + 1, symbolic]]; prepCache[ff]];
  Do[Module[{fIdx = fi, fam = fams[[fi]]},
    Do[Module[{blk, root, q},
      (* identify the block for this root via the running cursor *)
      root = fam["Roots"][[ri]]; q = root["BlockSize"];
      Do[Module[{P, init, rec, ls, comp, fbRun = startFb, WRun,
          attempts = 0, underflow, used, nextUsed},
        P = logCeiling[cs, root["a"], root["b"], qpos];
        (* init ladder: U[0, l] = e_{qpos+1-l} eps^{-l}, l = 0..qpos *)
        init = Table[
          Module[{vv = Table[esZero[wideTop], {d}]},
            If[l <= qpos,
              vv[[colCursor + (qpos + 1 - l)]] =
                esShift[esNew[0, PadRight[{1}, wideTop + 1 + l]], -l]];
            vv],
          {l, 0, P}];
        While[True,
          WRun = wideTop - fbRun + 1;
          attempts++;
          underflow = Catch[
            Block[{$adaptiveLowerFrameProbe = adaptiveQ && fbRun > terminalFb},
              runRecursion[cs, prepFor[fbRun], root["a"], root["b"],
                P, nmax, None, fbRun, WRun, init]],
            "DiffExp2AdaptiveLowerFrame"];
          If[!FailureQ[underflow], rec = underflow; Break[]];
          (* Monotone geometric widening.  Since terminalFb is exactly the
             previous scalar bound, termination and the old error behavior
             are preserved even for scalar/idempotent repeated poles. *)
          used = startFb - fbRun;
          nextUsed = If[used === 0, Max[1, singleUseDepth], 2 used];
          fbRun = Max[terminalFb, startFb - nextUsed]];
        AppendTo[adaptiveDiags, <|
          "Tag" -> {root["a"], root["b"], qpos},
          "Adaptive" -> adaptiveQ, "Attempts" -> attempts,
          "FrameBase" -> fbRun, "FrameTop" -> wideTop,
          "FrameWidth" -> WRun,
          "TerminalFrameBase" -> terminalFb,
          "TerminalFrameWidth" -> wideTop - terminalFb + 1|>];
        ls = assembleSolution[cs, root["a"], root["b"], rec, nmax];
        comp = compensatePseudoColumn[cs, ls, rec["Hits"], workColumns, reqMax];
        ls = comp[[1]];
        AppendTo[certs, certifyPseudoCompensation[cs, ls, rec["Hits"],
          <|"Kind" -> "Homogeneous", "Tag" ->
            {root["a"], root["b"], qpos}|>]];
        compAll = Join[compAll, Map[Append[#, "SourceColumn" ->
          (Length[workColumns] + 1)] &, comp[[2]]]];
        hitsAll = Join[hitsAll, publicPseudoHit /@ rec["Hits"]];
        (* Work-frame headroom is scratch space for Cauchy products and
           collision divisions, not a delivered accuracy claim. *)
        If[ls["EpsWindow", "CompleteMax"] < reqMax,
          err["E6", cs, <|"RequestedCompleteMax" -> reqMax,
            "DeliveredWindow" -> ls["EpsWindow"],
            "Tag" -> {root["a"], root["b"], qpos},
            "Detail" -> "homogeneous work budget did not reach the requested epsilon order"|>]];
        AppendTo[workColumns, ls];
        ls = capWindow[cs, ls, reqMax];
        AppendTo[columns, ls];
        AppendTo[specs, <|"a" -> root["a"], "b" -> root["b"], "p" -> qpos,
          "Family" -> fIdx, "Root" -> ri, "ChainPos" -> qpos|>]],
        {qpos, 0, q - 1}];
      colCursor += q],
      {ri, Length[fam["Roots"]]}]],
    {fi, Length[fams]}];
  Module[{fs = <|"Columns" -> columns, "Specs" -> specs,
      "Diagnostics" -> <|"PseudoCollisionsHit" -> hitsAll,
        "PseudoCompensations" -> compAll,
        "PseudoCollisionsCompensated" -> And @@ certs,
        "AdaptiveLowerFrames" -> adaptiveDiags,
        "CouplingDepth" -> 0|>|>},
    ODEResidualCheck[cs, fs];
    fs]];

SolveParticular[cs_Association, source_Association, req_Association] := Module[
  {d = cs["SystemSize"], nmax, reqMin, reqMax, parts = {}, wideTop, prep,
   eps = DiffExp2`Config`CanonicalEps[], hitsAll = {}, compAll = {}, certs = {},
   homTargets = None, symbolic, poleDepth, spectralDepth,
   inverseSpectralDepth, sourceFrameDepth},
  nmax = Min[req["TOrder"], source["TWindow", "CompleteMax"]];
  reqMin = req["EpsWindow", "Min"]; reqMax = req["EpsWindow", "CompleteMax"];
  If[source["Sectors"] === {},
    Return[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"],
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> Table[0, {reqMax - reqMin + 1}, {nmax + 1}, {d}]|>},
      "EpsWindow" -> <|"Min" -> reqMin, "CompleteMax" -> reqMax|>,
      "TWindow" -> <|"CompleteMax" -> nmax|>,
      "ErrorEstimate" -> ConstantArray[0, reqMax - reqMin + 1],
      "Prescriptions" -> cs["Prescriptions"]|>]];
  If[cs["Gauge"] =!= IdentityMatrix[cs["SystemSize"]],
    err["E8", cs, <|"Detail" ->
      "nonzero particular sources on a rank-reduced chart require the GaugeInverse source transform, which is not implemented"|>]];
  symbolic = clearedSymbolic[cs];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  spectralDepth = spectralTransformPoleDepth[cs];
  inverseSpectralDepth = inverseSpectralTransformPoleDepth[cs];
  sourceFrameDepth = spectralDepth + inverseSpectralDepth;
  Do[Module[{aS = sec["a"], bS = sec["b"], pS = sec["p"], arr = sec["Coeffs"],
      P, srcMin, srcMax, VInvExp, VInvVal, bHat, bHatValid, srcFn, rec, ls,
      wideTop2, prep2, pseudoDepth, comp, desiredMax},
    P = logCeiling[cs, aS, bS, pS, True];  (* sources: Z>=0 incl. the same-a hit *)
    pseudoDepth = pseudoDepthForTag[cs, aS, bS, 0];
    srcMin = source["EpsWindow", "Min"]; srcMax = source["EpsWindow", "CompleteMax"];
    wideTop2 = srcMax + P + pseudoDepth + 2 - Min[0, srcMin - P - 2];
    wideTop2 = Max[wideTop2, srcMax + P + pseudoDepth + poleDepth];
    wideTop2 += sourceFrameDepth;
    Module[{fb2 = Min[Min[srcMin, 0] - P - pseudoDepth - 2 - 2,
          srcMin - P - pseudoDepth - poleDepth], Wd2},
      (* The source first crosses VInv, then the solved coefficient crosses
         V on assembly.  Both matrices live in the shared frame, so their
         lower pole depths and upper scratch losses add. *)
      fb2 -= sourceFrameDepth;
      Wd2 = wideTop2 - fb2 + 1;
      prep2 = prepareCleared[cs, fb2, Wd2, symbolic];
      VInvExp = Map[ratEpsList[Together[#], eps, fb2, Wd2] &, cs["VInv"], {2}];
      VInvVal = Map[frameValuation[#, fb2] &, VInvExp, {2}];
      (* J-frame source as FRAME LISTS: bHat[n+1][r] *)
      bHat = Table[Module[{vc = Table[Module[{fl = ConstantArray[0, Wd2]},
          Do[If[1 <= k - fb2 + 1 <= Wd2,
            fl[[k - fb2 + 1]] = arr[[k - srcMin + 1, n + 1, c]]],
            {k, srcMin, srcMax}]; fl], {c, d}]},
        Table[Sum[frConv[VInvExp[[r, c]], vc[[c]], fb2, Wd2], {c, d}], {r, d}]],
        {n, 0, Min[nmax, Length[arr[[1]]] - 1]}];
      (* The source is finite data: zero frame slots above srcMax are
         unknown, not an extension of completeness.  Transform its actual
         top through the exact VInv entry valuations. *)
      bHatValid = Table[Table[validMin[Table[
          If[VInvVal[[r, c]] === Infinity, Infinity,
            validShift[srcMax, VInvVal[[r, c]], fb2 + Wd2 - 1]],
          {c, d}]], {r, d}], {Length[bHat]}];
      (* hand frame lists straight through: wrap as a function returning
         pre-framed vectors via a marker *)
      (* fresh formal names: bare Private n/l carry values from package
         evaluation and corrupt Function parameter lists *)
      srcFn = With[{bb = bHat, bv = bHatValid, pp = pS},
        Function[{srcN, srcL},
          If[srcL === pp && srcN + 1 <= Length[bb],
            <|"Frames" -> bb[[srcN + 1]],
              "Validity" -> bv[[srcN + 1]]|>, None]]];
      rec = runRecursionFramedSrc[cs, prep2, aS, bS, P, nmax, srcFn, fb2, Wd2, None]];
    hitsAll = Join[hitsAll, publicPseudoHit /@ rec["Hits"]];
    ls = assembleSolution[cs, aS, bS, rec, nmax];
    desiredMax = Min[ls["EpsWindow", "CompleteMax"], reqMax];
    If[rec["Hits"] =!= {},
      If[homTargets === None,
        Module[{hreq = Join[req, <|"EpsWindow" ->
            Join[req["EpsWindow"], <|"CompleteMax" ->
              Max[reqMax, source["EpsWindow", "CompleteMax"]]|>]|>]},
          homTargets = SolveHomogeneous[cs, hreq]["Columns"]]];
      comp = compensatePseudoColumn[cs, ls, rec["Hits"], homTargets, desiredMax];
      ls = comp[[1]];
      compAll = Join[compAll, Map[Append[#, "SourceSector" ->
        {aS, bS, pS}] &, comp[[2]]]];
      AppendTo[certs, certifyPseudoCompensation[cs, ls, rec["Hits"],
        <|"Kind" -> "Particular", "Tag" -> {aS, bS, pS}|>]],
      AppendTo[certs, True]];
    AppendTo[parts, ls]],
    {sec, source["Sectors"]}];
  Module[{ls = If[Length[parts] === 1, First[parts],
      DiffExp2`SectorSeries`CombineLocalSolutions[
        ConstantArray[1, Length[parts]], parts]]},
    ls = capWindow[cs, ls, reqMax];
    ls = Join[ls, <|"Diagnostics" -> <|
      "PseudoCollisionsHit" -> hitsAll,
      "PseudoCompensations" -> compAll,
      "PseudoCollisionsCompensated" -> And @@ certs|>|>];
    ODEResidualCheck[cs, ls, source];
    ls]];

SolveChart[cs_Association, req_Association, source_:None] := Module[{basis, part},
  basis = SolveHomogeneous[cs, req];
  part = If[source === None, None, SolveParticular[cs, source, req]];
  <|"Basis" -> basis, "Particular" -> part, "CouplingDepth" -> 0|>];

(* ---- value-vector propagation (regular charts; prototype) ----
   The incoming VALUE at the chart center is the t^0 Cauchy datum of the
   transported solution: ONE runRecursion with init = vals replaces the
   d-column basis + MatchWeights + CombineLocalSolutions of the basis
   path (the old engine's Currbcs chaining, DiffExp/Transport.m).
   Regular charts only: the single (0,0,0) family makes every n >= 1
   step CASE T (dA = n != 0).  Negative-valued Taylor multipliers can
   nevertheless lower the honest CompleteMax and are included in the work
   frame.  Invariants kept:
   exact (0,0,0) sector tag, honest window (capped at the incoming
   CompleteMax), no matching solve to go ambiguous, and the always-on
   ODE residual check runs on the propagated solution itself. *)
SolveValueRegular[cs_Association, req_Association, vals_List] := Module[
  {d = cs["SystemSize"], nmax, vMin, vCM, fb, wideTop, Wd, prep, rec, ls,
   symbolic, poleDepth},
  If[!TrueQ[Lookup[cs["IndicialData"], "Regular", False]],
    err["E8", cs, <|"Detail" ->
      "SolveValueRegular requires a regular chart (pole order 0); singular charts keep the basis+matching path"|>]];
  If[Length[vals] =!= d || !AllTrue[vals, DiffExp2`EpsSeries`ESQ],
    err["E8", cs, <|"Components" -> Length[vals], "Dimension" -> d,
      "Detail" -> "value vector must be d EpsSeries components"|>]];
  nmax = req["TOrder"];
  vMin = Min[esMin /@ vals]; vCM = Min[esCM /@ vals];
  If[AllTrue[vals, Function[v,
      AllTrue[Table[esCoeff[v, k], {k, esMin[v], esCM[v]}], # === 0 &]]],
    (* exactly-zero incoming value: the propagated solution is zero *)
    Return[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
      "Radius" -> cs["Radius"],
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> Table[0, {vCM - vMin + 1}, {nmax + 1}, {d}]|>},
      "EpsWindow" -> <|"Min" -> vMin, "CompleteMax" -> vCM|>,
      "TWindow" -> <|"CompleteMax" -> nmax|>,
      "ErrorEstimate" -> ConstantArray[0, vCM - vMin + 1],
      "Prescriptions" -> cs["Prescriptions"]|>]];
  (* frame: the SolveHomogeneous work-window shape with Pmax = cdMax = 0
     (regular chart), based on the INCOMING window instead of req *)
  symbolic = clearedSymbolic[cs];
  poleDepth = recurrencePoleDepth[symbolic, nmax];
  fb = Min[vMin, 0] - 2;
  wideTop = vCM + 2 - Min[0, vMin - 2];
  fb = Min[fb, vMin - poleDepth];
  wideTop = Max[wideTop, vCM + poleDepth];
  Wd = wideTop - fb + 1;
  prep = prepareCleared[cs, fb, Wd, symbolic];
  rec = runRecursion[cs, prep, 0, 0, 0, nmax, None, fb, Wd, {vals}];
  ls = assembleSolution[cs, 0, 0, rec, nmax];
  ls = capWindow[cs, ls, vCM];
  ODEResidualCheck[cs, ls];
  ls];

(* honesty cap: nothing above the incoming CompleteMax is complete (the
   frame's buffer slots above it absorb eps-shifts but are not delivered) *)
capWindow[cs_, ls_, capCM_] := Module[
  {kmin = ls["EpsWindow", "Min"], kmax = ls["EpsWindow", "CompleteMax"],
   keep, ncols, ncomp},
  If[kmax <= capCM, Return[ls]];
  If[capCM < kmin,
    {ncols, ncomp} = Dimensions[First[ls["Sectors"]]["Coeffs"]][[2 ;; 3]];
    Return[Join[ls, <|
      "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
        "Coeffs" -> ConstantArray[0, {1, ncols, ncomp}]|>},
      "EpsWindow" -> <|"Min" -> capCM, "CompleteMax" -> capCM|>,
      "ErrorEstimate" -> {0}|>], Module]];
  keep = capCM - kmin + 1;
  Join[ls, <|
    "Sectors" -> Map[Append[#, "Coeffs" -> Take[#["Coeffs"], keep]] &,
      ls["Sectors"]],
    "EpsWindow" -> <|"Min" -> kmin, "CompleteMax" -> capCM|>,
    "ErrorEstimate" -> If[ListQ[ls["ErrorEstimate"]],
      Take[PadRight[ls["ErrorEstimate"], kmax - kmin + 1], keep],
      ls["ErrorEstimate"]]|>]];

(* ---- residual check ---- *)

(* Residual-check numeric handoff (the Transport`numHandoff policy): the
   evaluated f / theta-f VALUES feeding the d^2 esTimes grid numericize at
   2x WP.  Exact evaluation outputs otherwise compound into exact-rational
   giants across the grid (measured: the check was ~3.1 s of a 13 s chart
   at d = 7, Docs/SpeedIdeas.md §2).  Values only — the residual compare
   below is numeric by construction and uses componentwise uncertainty
   bounds; windows untouched. *)
numV[s_] := esNew[esMin[s],
  N[Table[esCoeff[s, k], {k, esMin[s], esCM[s]}],
    DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]];

ODEResidualCheck[cs_Association, sol_Association, source_:None, probe_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], t0, sols, maxRel = 0,
   Bsub, bMinVal, bt0Cache = <||>, bt0For, win, rtol},
  rtol = DiffExp2`Tolerances`Tol["ResidTol"];
  sols = If[KeyExistsQ[sol, "Columns"], sol["Columns"], {sol}];
  (* truncation-aware probe: the residual of a degree-nmax truncation is
     O((t0/R)^(nmax+1)); place t0 so that tail sits below rtol/100 *)
  t0 = If[probe === Automatic,
    Module[{nmaxS, frac, fracRaw, raw},
      nmaxS = Min[#["TWindow", "CompleteMax"] & /@ sols];
      fracRaw = N[(rtol/100)^(1/(nmaxS + 1)), 30];
      (* A fixed absolute Rationalize tolerance rounded the small-order
         probe to exactly zero (for example nmax=0..2 at ResidTol=10^-10),
         sending singular sectors through the forbidden origin limit. *)
      frac = Min[1/4, Rationalize[fracRaw, Abs[fracRaw]/50]];
      raw = cs["Radius"]*frac*(2 + Mod[Hash[cs["Center"]], 5])/10;
      (* a SIMPLE RATIONAL probe: algebraic radii otherwise drag every
         downstream evaluation into exact algebraic arithmetic *)
      Rationalize[N[raw, 20], N[raw, 20]/50]],
    probe];
  (* B(t0) is column-independent: hoist the d^2 exact substitutions out of
     the per-column loop (they dominated the check at d = 7: ~7x recompute)
     and memoize the eps-frame expansion per (k1, k2) window — columns
     share windows in the common case.  Values are identical per column;
     pure cost hoist, no honesty change. *)
  Bsub = Map[Together[# /. t -> t0] &,
    Lookup[cs, "ThetaOriginal", cs["ThetaMatrix"]], {2}];
  bMinVal = Module[{v = Select[
      epsValuation[#, eps] & /@ Flatten[Bsub], IntegerQ]},
    If[v === {}, 0, Min[0, Min[v]]]];
  bt0For[k1_, k2_] := Module[
    {opTop = Max[0, k2 - k1], key},
    key = {bMinVal, opTop};
    If[!KeyExistsQ[bt0Cache, key],
      bt0Cache[key] = Map[Module[
          {v = epsValuation[#, eps], fl},
        If[v === Infinity || v > opTop, esZero[opTop],
          fl = ratEpsList[#, eps, bMinVal, opTop - bMinVal + 1];
          (* Leading frame zeros are padding, not part of the operator's
             mathematical Laurent support. *)
          esTrim[esNew[bMinVal, fl]]]] &, Bsub, {2}]];
    bt0Cache[key]];
  Do[Module[{f, df, lhs, rhs, srcv, k1, k2, scale, Bt0},
    If[Environment["DEBUG_RESID"] === "1",
      Print["RESID col window=", ls["EpsWindow"], " tags=",
        {#["a"], #["b"], #["p"]} & /@ ls["Sectors"]]];
    f = numV[DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0,
      "UsePade" -> False]["Value"]];
    df = numV[DiffExp2`SectorSeries`EvaluateLocalSolution[
      DiffExp2`SectorSeries`DifferentiateLocalSolution[ls], t0,
      "UsePade" -> False]["Value"]];
    (* Coefficients below an EpsSeries Min are certified zero.  Residuals
       must therefore start at the UNION of supports, not their intersection;
       otherwise a low-order B.f term can be hidden by a zero derivative. *)
    k1 = Min[esMin[f], esMin[df]];
    k2 = Min[esCM[f], esCM[df]];
    Bt0 = bt0For[esMin[f], k2];
    (* theta f = t f' *)
    Module[{thetaF, Bf, resid},
      thetaF = esScale[t0, df];
      Bf = Module[{fv = f, comps},
        comps = Length[esCoeff[fv, esMin[fv]]];
        Table[Module[{s = None},
          Do[Module[{cESr},
            If[!AllTrue[Bt0[[r, c]]["Coeffs"], # === 0 &],
              cESr = esNew[esMin[fv],
                Table[esCoeff[fv, k][[c]], {k, esMin[fv], esCM[fv]}]];
              Module[{term = esTimes[Bt0[[r, c]], cESr]},
                s = If[s === None, term, esAdd[s, term]]]]],
            {c, comps}];
          If[s === None, esZero[k2], s]], {r, comps}]];
      srcv = If[source =!= None && AssociationQ[source] && KeyExistsQ[source, "Sectors"],
        Module[{sLS = Join[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
            "Radius" -> cs["Radius"], "ErrorEstimate" -> 0,
            "Prescriptions" -> cs["Prescriptions"],
            "EpsWindow" -> source["EpsWindow"], "TWindow" -> source["TWindow"]|>,
            <|"Sectors" -> source["Sectors"]|>]},
          DiffExp2`SectorSeries`EvaluateLocalSolution[sLS, t0, "UsePade" -> False]["Value"]],
        None];
      Do[Module[{tf = esCoeff[thetaF, k], bfk, sv, residualEntries,
          badEntries, rv, sc},
        bfk = Table[If[esMin[Bf[[r]]] <= k <= esCM[Bf[[r]]], esCoeff[Bf[[r]], k], 0],
          {r, Length[Bf]}];
        sv = If[srcv =!= None && esMin[srcv] <= k <= esCM[srcv],
          esCoeff[srcv, k], ConstantArray[0, Length[Bf]]];
        residualEntries = Flatten[{tf - bfk - sv}];
        badEntries = Select[residualEntries,
          !NumericQ[#] &&
            !(FreeQ[#, _?InexactNumberQ] && TrueQ[PossibleZeroQ[#]]) &];
        If[badEntries =!= {},
          err["E7", cs, <|"Probe" -> t0, "EpsOrder" -> k,
            "Residual" -> First[badEntries],
            "Detail" -> "ODE residual contains uncertified nonnumeric content"|>]];
        rv = Select[residualEntries, NumericQ];
        sc = Max[1, Sequence @@ (numMag[#, 20] & /@
          Select[Flatten[{tf, bfk, sv}], NumericQ])];
        Module[{rk = Max[0, Sequence @@
              (Last[DiffExp2`Tolerances`NumericMagnitudeBounds[#, 20]] & /@ rv)]/sc},
          If[Environment["DEBUG_RESID"] === "1" && rk > 10^-12,
            Print["RESID k=", k, " rel=", N[rk, 4]]];
          maxRel = Max[maxRel, rk]]],
        {k, Min[k1, esMin[thetaF], Min[esMin /@ Bf],
            If[srcv === None, k1, esMin[srcv]]],
          Min[k2, esCM[thetaF], Min[esCM /@ Bf],
            If[srcv === None, k2, esCM[srcv]]]}]]],
    {ls, sols}];
  If[maxRel > rtol,
    err["E7", cs, <|"Probe" -> t0, "MaxRelativeResidual" -> N[maxRel, 6],
      "ResidTol" -> rtol, "Detail" -> "ODE residual check failed"|>]];
  maxRel];

End[];
EndPackage[];
