(* DiffExp2/Solve.m — THE one local solver: symbolic-eps Frobenius with
   EpsSeries coefficient arithmetic.  Spec: Docs/specs/Solve.md (binding);
   decisions: DECISIONS-M0.md.
   True resonance -> explicit log-chain ladder (CASE R); pseudo-resonance ->
   exact Laurent shift (CASE P); regular charts run the SAME recursion
   (CASE T only, V = I).
   V1 DEVIATION (recorded; spec 3.7): pseudo-resonant columns are delivered
   with honest LAURENT eps-windows (hits in Diagnostics) instead of the
   compensated eps-regular combination; Transport's eps-graded matching
   consumes Laurent columns by design.  Compensation lands with the M5
   evidence if matching needs it. *)

BeginPackage["DiffExp2`Solve`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`"}];

PrepareChart::usage = "PrepareChart[sys, chart] applies the chart map, runs ChartIndicial, and assembles the ChartSystem (theta matrix, gauge, V/VInv spectral frame, families).";
SolveHomogeneous::usage = "SolveHomogeneous[chartSystem, req] gives the FundamentalSystem: one LocalSolution column per indicial sector spec.";
SolveParticular::usage = "SolveParticular[chartSystem, source, req] gives THE particular solution (canonical kernel choice) for a sector-native theta-form source.";
SolveChart::usage = "SolveChart[chartSystem, req, source] gives <|\"Basis\", \"Particular\", \"CouplingDepth\"|>.";
ODEResidualCheck::usage = "ODEResidualCheck[chartSystem, sol, source, probe] checks the theta-form ODE residual at an interior probe point; loud error above ResidTol.";

Begin["`Private`"];

err[id_, cs_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Solve",
    "Chart" -> ToString[Lookup[cs, "Center", "?"], InputForm]|>, payload]];
cfg = DiffExp2`Config`CFG;
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
PrepareChart[sys_Association, chart_Association] := Module[
  {pcKey = {Hash[sys["Matrix"]], chart["Center"], Lookup[chart, "Scale", 1],
    Lookup[chart, "Radius", None]}},
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
      "CollisionDepth" -> Length[Select[collisions, #["Type"] === "LaurentShift" &]]|>]],
    {fam, idata["Families"]}];
  V = Transpose[cols];
  detV = Together[Det[V]];
  If[zeroCanQ[detV], err["E2", chart, <|"Detail" -> "spectral frame V is singular"|>]];
  VInv = Map[Cancel[Together[#]] &, Inverse[V], {2}];
  <|"ChartVar" -> t, "Center" -> x0,
    "ChartMap" -> <|"Center" -> x0, "Scale" -> beta|>,
    "Radius" -> chart["Radius"],
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

(* ---- the recursion ---- *)

(* ===================== PACKED-FRAME RECURSION CORE =====================
   Coefficient vectors are PLAIN LISTS over a fixed eps-frame
   [fb, fb+W-1] (index i <-> eps^(fb+i-1)); products via ListConvolve.
   This bypasses the per-op EpsSeries object layer (Association churn +
   per-coefficient Together) that cost ~18s/chart on 3-master systems.
   Honest windows: the frame is sized by the deterministic work-window
   formula; eps-divisions decrement a topValid watermark; the delivered
   EpsWindow is [content-min, topValid].  EpsSeries objects appear only
   at the boundaries (sources in, LocalSolutions out). *)

(* fast exact eps-expansion of a RATIONAL expr into a frame list;
   ByteCount-gated numericization at 2x WP *)
ratEpsList[expr_, eps_, fb_, W_] := Module[
  {c, num, den, vn, vd, v, nc, dc, rel, out, wp2, top},
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
  rel = top - v;
  Module[{ncl = Table[Coefficient[num, eps, j], {j, 0, rel}],
    dcl = Table[Coefficient[den, eps, j], {j, 0, rel}], csr},
    csr = ConstantArray[0, rel + 1];
    csr[[1]] = Together[ncl[[1]]/dcl[[1]]];
    Do[csr[[m + 1]] = Together[
        (ncl[[m + 1]] - Sum[dcl[[j + 1]]*csr[[m - j + 1]], {j, 1, m}])/dcl[[1]]],
      {m, 1, rel}];
    wp2 = DiffExp2`Config`CFG["WorkingPrecision"] + 20;
    csr = Map[If[# === 0 || ByteCount[#] <= 500, #, N[#, wp2]] &, csr];
    Do[out[[v - fb + 1 + m]] = csr[[m + 1]], {m, 0, rel}]];
  out];

(* framed product: a, b based at fb -> product re-framed at fb *)
frConv[a_, b_, fb_, W_] :=
  Take[ListConvolve[a, b, {1, -1}, 0], {1 - fb, W - fb}];

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
frDivEps[a_] := Append[Rest[a], 0];

prepareCleared[cs_, fb_, W_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], d = cs["SystemSize"],
   den, num, d0, dD, dN, dL, NhatL},
  den = Together[PolynomialLCM @@ (Denominator[Together[#]] & /@ Flatten[cs["ThetaMatrix"]])];
  num = Map[Cancel[Together[#*den]] &, cs["ThetaMatrix"], {2}];
  d0 = Together[den /. t -> 0];
  If[zeroQ[d0 /. eps -> 0],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  dL = Table[ratEpsList[Coefficient[den, t, j], eps, fb, W], {j, 0, dD}];
  NhatL = Table[Module[{Nj, Nhat},
    Nj = Map[Coefficient[#, t, j] &, num, {2}];
    Nhat = Map[Cancel[Together[#]] &, cs["VInv"] . Nj . cs["V"], {2}];
    Map[ratEpsList[#, eps, fb, W] &, Nhat, {2}]], {j, 0, dN}];
  (* sparse-offset forms: den coefficients and Nhat are eps-POLYNOMIAL
     (cleared system) -> 1-3 nonzero eps-orders each.  Each multiplier
     becomes a short list of {eps-shift, scalar} / {eps-shift, d x d
     matrix}; the recursion then runs on whole (d x W) blocks with
     matrix Dot + column shifts instead of per-entry convolutions. *)
  Module[{dSp, NhatSp},
    dSp = Table[Module[{lst = dL[[j + 1]]},
      Table[{i + fb - 1, lst[[i]]}, {i, Select[Range[W], lst[[#]] =!= 0 &]}]],
      {j, 0, dD}];
    NhatSp = Table[Module[{idxs},
      idxs = Select[Range[W], Module[{i = #},
        AnyTrue[Flatten[NhatL[[j + 1]], 1], #[[i]] =!= 0 &]] &];
      Table[{i + fb - 1, Map[#[[i]] &, NhatL[[j + 1]], {2}]}, {i, idxs}]],
      {j, 0, dN}];
    <|"dL" -> dL, "NhatL" -> NhatL, "dD" -> dD, "dN" -> dN,
      "dSp" -> dSp, "NhatSp" -> NhatSp|>]];

(* EpsSeries -> frame list *)
esToFrame[s_, fb_, W_] := Module[{out = ConstantArray[0, W], m1, m2},
  m1 = esMin[s]; m2 = esCM[s];
  Do[If[1 <= k - fb + 1 <= W, out[[k - fb + 1]] = esCoeff[s, k]], {k, m1, m2}];
  out];

(* CASE T/P block solve on frames *)
blockSolveTPFrame[rhs_, deltaList_, invD0_, q_, fb_, W_] := Module[
  {inv = frInv[deltaList, fb, W]},
  Table[Module[{acc = ConstantArray[0, W], pw = inv},
    Do[
      If[r + m <= q, acc = acc + frConv[pw, rhs[[r + m]], fb, W]];
      If[m < q - 1, pw = frConv[pw, inv, fb, W]],
      {m, 0, q - 1}];
    frConv[acc, invD0, fb, W]], {r, q}]];

(* the recursion on frames.  Returns <|"U" -> U (frame lists), "Hits",
   "P", "FrameBase", "TopValid"|>.  runRecursionFramedSrc = same entry
   point; srcF auto-detects pre-framed source vectors. *)
runRecursionFramedSrc[args___] := runRecursion[args];
runRecursion[cs_, prep_, aT_, bT_, P_, nmax_, srcHat_, fb_, W_, init_] := Module[
  {d = cs["SystemSize"], dL = prep["dL"], NhatL = prep["NhatL"],
   dD = prep["dD"], dN = prep["dN"], invD0, blocks, U, hits = {}, n0,
   topValid, lamL, zeroV, srcF},
  invD0 = frInv[dL[[1]], fb, W];
  blocks = blockList[cs];
  (* negative eps-leads in the sparse multipliers (Nhat = VInv.N.V is
     rational in eps; the frame can carry eps-poles) erode the top |s|
     slots of every product - a STABLE depth, paid once *)
  topValid = fb + W - 1 - Max[0, Max[Map[Module[{offs = #[[All, 1]]},
    If[offs === {}, 0, -Min[offs]]] &,
    Join[prep["dSp"], prep["NhatSp"]]]]];
  zeroV = ConstantArray[0, W];
  lamL[m_] := Module[{l = zeroV}, l[[-fb + 1]] = Together[aT + m];
    If[-fb + 2 <= W, l[[-fb + 2]] = Together[bT]]; l];
  srcF[n_, l_] := If[srcHat === None, None,
    Module[{sv = srcHat[n, l]},
      Which[
        sv === None, None,
        (* already frame lists (vectors of plain lists of length W) *)
        VectorQ[First[sv], AtomQ] || (ListQ[First[sv]] && Length[First[sv]] === W &&
          !AssociationQ[First[sv]]), sv,
        True, Map[esToFrame[#, fb, W] &, sv]]]];
  U = Table[Table[Table[zeroV, {d}], {l, 0, P + 1}], {n, 0, nmax}];
  n0 = If[init === None, 0, 1];
  If[init =!= None,
    Do[U[[1, l + 1]] = Map[esToFrame[#, fb, W] &, init[[l + 1]]],
      {l, 0, Min[P, Length[init] - 1]}]];
  (* shift a (d x W) block by s eps-slots (s > 0: content moves UP) *)
  Module[{shB, dSp = prep["dSp"], NhatSp = prep["NhatSp"]},
  shB[M_, s_] := Which[s === 0, M,
    s > 0, Join[ConstantArray[0, {Length[M], s}], M[[All, 1 ;; W - s]], 2],
    True, Join[M[[All, -s + 1 ;; W]], ConstantArray[0, {Length[M], -s}], 2]];
  Do[Module[{R, rhsFull},
    R = Table[Module[{acc = ConstantArray[0, {d, W}]},
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]]},
          Do[acc += sp[[2]] . shB[Ublk, sp[[1]]], {sp, NhatSp[[j + 1]]}]]],
        {j, 1, Min[n, dN]}];
      Do[If[n - j >= 0,
        Module[{Ublk = U[[n - j + 1, l + 1]], term},
          (* lam*U + eps*U_above: lam = (aT + n - j) + bT*eps *)
          term = Together[aT + n - j]*Ublk + shB[Together[bT]*Ublk, 1] +
            shB[U[[n - j + 1, l + 2]], 1];
          Do[acc -= sp[[2]]*shB[term, sp[[1]]], {sp, dSp[[j + 1]]}]]],
        {j, 1, Min[n, dD]}];
      If[srcHat =!= None,
        Do[If[n - j >= 0,
          Module[{sv = srcF[n - j, l]},
            If[sv =!= None,
              Do[acc += sp[[2]]*shB[sv, sp[[1]]], {sp, dSp[[j + 1]]}]]]],
          {j, 0, Min[n, dD]}]];
      acc], {l, 0, P}];
    (* solve top-down in l *)
    Do[Module[{},
      rhsFull = R[[l + 1]];
      Do[rhsFull -= sp[[2]]*shB[shB[U[[n + 1, l + 2]], 1], sp[[1]]],
        {sp, dSp[[1]]}];
      rhsFull = Table[rhsFull[[r]], {r, d}];
      Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"],
          dA, dB, deltaList},
        dA = Together[aT + n - aI]; dB = Together[bT - bI];
        Which[
          !zeroCanQ[dA],
          deltaList = Module[{l2 = zeroV}, l2[[-fb + 1]] = dA;
            If[-fb + 2 <= W, l2[[-fb + 2]] = dB]; l2];
          U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], deltaList,
              invD0, q, fb, W]]],
          !zeroCanQ[dB],
          (deltaList = Module[{l2 = zeroV},
            If[-fb + 2 <= W, l2[[-fb + 2]] = dB]; l2];
          U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTPFrame[rhsFull[[colsB]], deltaList,
              invD0, q, fb, W]]];
          topValid -= q;
          If[l == 0, AppendTo[hits, <|"n" -> n, "Cols" -> colsB, "DeltaB" -> dB|>]]),
          True, Null]],
        {blk, blocks}]],
      {l, P, 0, -1}];
    (* CASE R ladder *)
    Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"]},
      If[zeroQ[aT + n - aI] && zeroQ[bT - bI],
        Module[{Rt, assigned = Association[]},
          Rt = Table[Table[
            frConv[R[[l + 1, colsB[[r]]]], invD0, fb, W], {r, q}], {l, 0, P}];
          Do[If[l + 1 <= P,
            assigned[{l + 1, q}] = frDivEps[Rt[[l + 1, q]]];
            topValid -= 1],
            {l, 0, P - 1}];
          Do[
            Do[Module[{upr = Lookup[assigned, Key[{l + 1, r}], zeroV]},
              If[!KeyExistsQ[assigned, {l, r + 1}],
                assigned[{l, r + 1}] = frShiftUp[upr] - Rt[[l + 1, r]]]],
              {r, 1, q - 1}],
            {l, P, 0, -1}];
          Do[If[KeyExistsQ[assigned, {l, r}],
            U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
              colsB[[r]] -> assigned[{l, r}]]],
            {l, 0, P}, {r, 1, q}]]]],
      {blk, blocks}]],
    {n, n0, nmax}]];
  <|"U" -> U, "Hits" -> hits, "P" -> P, "FrameBase" -> fb, "TopValid" -> topValid|>];

(* ---- assembly ---- *)

(* frame U -> original-frame LocalSolution with sectors (a, b, l) *)
assembleSolution[cs_, aT_, bT_, rec_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], d = cs["SystemSize"], U = rec["U"],
   P = rec["P"], fb = rec["FrameBase"], W, kminO, kmaxO, VexpL, WL, secs, ls,
   nzQ, rec2TopShift = 0},
  W = Length[U[[1, 1, 1]]];
  VexpL = Map[ratEpsList[Together[#], eps, fb, W] &, cs["V"], {2}];
  (* a NEGATIVE-lead V entry (rational-in-eps chain vectors) shifts content
     down AND makes the product's top |lead| slots incomplete: the delivered
     window pays for it *)
  rec2TopShift = Max[0, Max[Flatten[Map[
    Function[entry, Module[{pos = SelectFirst[Range[Length[entry]],
        entry[[#]] =!= 0 &, None]},
      If[pos === None, 0, Max[0, -(fb + pos - 1)]]]],
    VexpL, {2}]]]];
  WL = Table[
    Table[Sum[frConv[VexpL[[r, c]], U[[n + 1, l + 1, c]], fb, W], {c, d}], {r, d}],
    {n, 0, nmax}, {l, 0, P}];
  (* content min: lowest frame index with a nonzero entry anywhere *)
  nzQ[i_] := AnyTrue[Flatten[WL, 2], #[[i]] =!= 0 &];
  kminO = Module[{i = 1}, While[i < W && !nzQ[i], i++]; fb + i - 1];
  kmaxO = Min[rec["TopValid"], fb + W - 1] - rec2TopShift;
  If[kmaxO < kminO, kmaxO = kminO];
  secs = Table[<|"a" -> aT, "b" -> bT, "p" -> l,
    "Coeffs" -> Table[Table[Table[WL[[n + 1, l + 1, r, k - fb + 1]], {r, d}],
      {n, 0, nmax}], {k, kminO, kmaxO}]|>,
    {l, 0, P}];
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
   gv, kmin, kmax, newSecs, vT, fbT, WG, TexpL},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[ls]];
  d = Length[T];
  gv = Min[tVal[#, t] & /@ Flatten[T]];
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  (* eps-valuation of T entries (rational in eps) *)
  vT = Min[0, Min[Map[Module[{c2 = Cancel[Together[#]]},
    If[zeroCanQ[c2], 0,
      Exponent[Numerator[c2] /. t -> 1/2, eps, Min] -
        Exponent[Denominator[c2] /. t -> 1/2, eps, Min]]] &, Flatten[T]]]];
  fbT = vT; WG = (kmax - kmin) + (-vT) + 1;
  (* TexpL[m - gv + 1][r][c] = eps-frame list (base fbT, width WG) of the
     t^m Laurent coefficient of T_rc *)
  TexpL = Table[Map[Module[{lc = tLaurent[#, t, m]},
      (* tLaurent output is canonical (0 or a Cancel[Together[...]]) *)
      If[zeroCanQ[lc], ConstantArray[0, WG], ratEpsList[lc, eps, fbT, WG]]] &,
      T, {2}],
    {m, gv, nmax}];
  newSecs = Map[Module[{arr = #["Coeffs"], aS = #["a"], outF, kminN, kmaxN,
      ncolsS = Dimensions[#["Coeffs"]][[2]], topValidG},
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
    kminN = Module[{i = 1, found = False, fbG = kmin + fbT},
      While[i <= WG && !found,
        If[AnyTrue[Flatten[outF[[All, All, i]]], # =!= 0 &], found = True, i++]];
      kmin + fbT + i - 1];
    kmaxN = Max[topValidG, kminN];
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

(* ---- public functions ---- *)

SolveHomogeneous[cs_Association, req_Association] := Module[
  {d = cs["SystemSize"], blocks = blockList[cs], nmax, reqMin, reqMax,
   columns = {}, specs = {}, hitsAll = {}, fams = cs["Families"], colCursor = 0,
   wideTop, prep, Pmax, cdMax, fb, Wd},
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"]; reqMax = req["EpsWindow", "CompleteMax"];
  Pmax = Max[0, Max[Table[logCeiling[cs, b["a"], b["b"], b["q"] - 1], {b, blocks}]]];
  cdMax = Max[0, Max[#["CollisionDepth"] & /@ fams]];
  wideTop = reqMax + Pmax + cdMax + 2 - Min[0, reqMin - Pmax - 2];
  fb = Min[reqMin, 0] - Pmax - cdMax - 2;   (* frame base; eps-poles of the
    cleared system are caught by ratEpsList's below-frame assert *)
  Wd = wideTop - fb + 1;
  prep = prepareCleared[cs, fb, Wd];
  Do[Module[{fIdx = fi, fam = fams[[fi]]},
    Do[Module[{blk, root, q},
      (* identify the block for this root via the running cursor *)
      root = fam["Roots"][[ri]]; q = root["BlockSize"];
      Do[Module[{P, init, rec, ls},
        P = logCeiling[cs, root["a"], root["b"], qpos];
        (* init ladder: U[0, l] = e_{qpos+1-l} eps^{-l}, l = 0..qpos *)
        init = Table[
          Module[{vv = Table[esZero[wideTop], {d}]},
            If[l <= qpos,
              vv[[colCursor + (qpos + 1 - l)]] =
                esShift[esNew[0, PadRight[{1}, wideTop + 1 + l]], -l]];
            vv],
          {l, 0, P}];
        rec = runRecursion[cs, prep, root["a"], root["b"], P, nmax, None, fb, Wd, init];
        hitsAll = Join[hitsAll, rec["Hits"]];
        ls = assembleSolution[cs, root["a"], root["b"], rec, nmax];
        AppendTo[columns, ls];
        AppendTo[specs, <|"a" -> root["a"], "b" -> root["b"], "p" -> qpos,
          "Family" -> fIdx, "Root" -> ri, "ChainPos" -> qpos|>]],
        {qpos, 0, q - 1}];
      colCursor += q],
      {ri, Length[fam["Roots"]]}]],
    {fi, Length[fams]}];
  Module[{fs = <|"Columns" -> columns, "Specs" -> specs,
      "Diagnostics" -> <|"PseudoCollisionsHit" -> hitsAll,
        "CouplingDepth" -> 0|>|>},
    ODEResidualCheck[cs, fs];
    fs]];

SolveParticular[cs_Association, source_Association, req_Association] := Module[
  {d = cs["SystemSize"], nmax, reqMin, reqMax, parts = {}, wideTop, prep,
   eps = DiffExp2`Config`CanonicalEps[], hitsAll = {}},
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
  Do[Module[{aS = sec["a"], bS = sec["b"], pS = sec["p"], arr = sec["Coeffs"],
      P, srcMin, srcMax, VInvExp, bHat, srcFn, rec, ls, wideTop2, prep2},
    P = logCeiling[cs, aS, bS, pS, True];  (* sources: Z>=0 incl. the same-a hit *)
    srcMin = source["EpsWindow", "Min"]; srcMax = source["EpsWindow", "CompleteMax"];
    wideTop2 = srcMax + P + 2 - Min[0, srcMin - P - 2];
    Module[{fb2 = Min[srcMin, 0] - P - 2 - 2, Wd2},
      Wd2 = wideTop2 - fb2 + 1;
      prep2 = prepareCleared[cs, fb2, Wd2];
      VInvExp = Map[ratEpsList[Together[#], eps, fb2, Wd2] &, cs["VInv"], {2}];
      (* J-frame source as FRAME LISTS: bHat[n+1][r] *)
      bHat = Table[Module[{vc = Table[Module[{fl = ConstantArray[0, Wd2]},
          Do[If[1 <= k - fb2 + 1 <= Wd2,
            fl[[k - fb2 + 1]] = arr[[k - srcMin + 1, n + 1, c]]],
            {k, srcMin, srcMax}]; fl], {c, d}]},
        Table[Sum[frConv[VInvExp[[r, c]], vc[[c]], fb2, Wd2], {c, d}], {r, d}]],
        {n, 0, Min[nmax, Length[arr[[1]]] - 1]}];
      (* hand frame lists straight through: wrap as a function returning
         pre-framed vectors via a marker *)
      (* fresh formal names: bare Private n/l carry values from package
         evaluation and corrupt Function parameter lists *)
      srcFn = With[{bb = bHat, pp = pS},
        Function[{srcN, srcL},
          If[srcL === pp && srcN + 1 <= Length[bb], bb[[srcN + 1]], None]]];
      rec = runRecursionFramedSrc[cs, prep2, aS, bS, P, nmax, srcFn, fb2, Wd2, None]];
    hitsAll = Join[hitsAll, rec["Hits"]];
    AppendTo[parts, assembleSolution[cs, aS, bS, rec, nmax]]],
    {sec, source["Sectors"]}];
  Module[{ls = If[Length[parts] === 1, First[parts],
      DiffExp2`SectorSeries`CombineLocalSolutions[
        ConstantArray[1, Length[parts]], parts]]},
    ODEResidualCheck[cs, ls, source];
    ls]];

SolveChart[cs_Association, req_Association, source_:None] := Module[{basis, part},
  basis = SolveHomogeneous[cs, req];
  part = If[source === None, None, SolveParticular[cs, source, req]];
  <|"Basis" -> basis, "Particular" -> part, "CouplingDepth" -> 0|>];

(* ---- residual check ---- *)

ODEResidualCheck[cs_Association, sol_Association, source_:None, probe_:Automatic] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], t0, sols, maxRel = 0,
   Bt0, win, rtol},
  rtol = DiffExp2`Tolerances`Tol["ResidTol"];
  sols = If[KeyExistsQ[sol, "Columns"], sol["Columns"], {sol}];
  (* truncation-aware probe: the residual of a degree-nmax truncation is
     O((t0/R)^(nmax+1)); place t0 so that tail sits below rtol/100 *)
  t0 = If[probe === Automatic,
    Module[{nmaxS, frac, raw},
      nmaxS = Min[#["TWindow", "CompleteMax"] & /@ sols];
      frac = Min[1/4, Rationalize[N[(rtol/100)^(1/(nmaxS + 1))], 10^-3]];
      raw = cs["Radius"]*frac*(2 + Mod[Hash[cs["Center"]], 5])/10;
      (* a SIMPLE RATIONAL probe: algebraic radii otherwise drag every
         downstream evaluation into exact algebraic arithmetic *)
      Rationalize[N[raw, 20], N[raw, 20]/50]],
    probe];
  Do[Module[{f, df, lhs, rhs, srcv, k1, k2, scale},
    If[Environment["DEBUG_RESID"] === "1",
      Print["RESID col window=", ls["EpsWindow"], " tags=",
        {#["a"], #["b"], #["p"]} & /@ ls["Sectors"]]];
    f = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0, "UsePade" -> False];
    df = DiffExp2`SectorSeries`EvaluateLocalSolution[
      DiffExp2`SectorSeries`DifferentiateLocalSolution[ls], t0, "UsePade" -> False];
    k1 = Max[esMin[f["Value"]], esMin[df["Value"]]];
    k2 = Min[esCM[f["Value"]], esCM[df["Value"]]];
    Bt0 = Map[Together[# /. t -> t0] &,
      Lookup[cs, "ThetaOriginal", cs["ThetaMatrix"]], {2}];
    Bt0 = Map[Module[{fl = ratEpsList[#, eps, Min[k1, 0], k2 - Min[k1, 0] + 1]},
      DiffExp2`EpsSeries`ESNew[Min[k1, 0], fl]] &, Bt0, {2}];
    (* theta f = t f' *)
    Module[{thetaF, Bf, resid},
      thetaF = esScale[t0, df["Value"]];
      Bf = Module[{fv = f["Value"], comps},
        comps = Length[esCoeff[fv, esMin[fv]]];
        Table[Module[{s = None},
          Do[Module[{cESr},
            cESr = esNew[esMin[fv], Table[esCoeff[fv, k][[c]], {k, esMin[fv], esCM[fv]}]];
            Module[{term = esTimes[Bt0[[r, c]], cESr]},
              s = If[s === None, term, esAdd[s, term]]]],
            {c, comps}];
          s], {r, comps}]];
      srcv = If[source =!= None && AssociationQ[source] && KeyExistsQ[source, "Sectors"],
        Module[{sLS = Join[<|"Center" -> cs["Center"], "ChartMap" -> cs["ChartMap"],
            "Radius" -> cs["Radius"], "ErrorEstimate" -> 0,
            "Prescriptions" -> cs["Prescriptions"],
            "EpsWindow" -> source["EpsWindow"], "TWindow" -> source["TWindow"]|>,
            <|"Sectors" -> source["Sectors"]|>]},
          DiffExp2`SectorSeries`EvaluateLocalSolution[sLS, t0, "UsePade" -> False]["Value"]],
        None];
      Do[Module[{tf = esCoeff[thetaF, k], bfk, sv, rv, sc},
        bfk = Table[If[esMin[Bf[[r]]] <= k <= esCM[Bf[[r]]], esCoeff[Bf[[r]], k], 0],
          {r, Length[Bf]}];
        sv = If[srcv =!= None && esMin[srcv] <= k <= esCM[srcv],
          esCoeff[srcv, k], ConstantArray[0, Length[Bf]]];
        rv = Select[Flatten[{tf - bfk - sv}], NumericQ];
        sc = Max[1, Sequence @@ (Abs[N[#, 20]] & /@ Select[Flatten[{tf, bfk, sv}], NumericQ])];
        Module[{rk = Max[0, Sequence @@ (Abs[N[#, 20]] & /@ rv)]/sc},
          If[Environment["DEBUG_RESID"] === "1" && rk > 10^-12,
            Print["RESID k=", k, " rel=", N[rk, 4]]];
          maxRel = Max[maxRel, rk]]],
        {k, Max[k1, esMin[thetaF]], Min[k2, esCM[thetaF],
          Min[esCM /@ Bf]]}]]],
    {ls, sols}];
  If[maxRel > rtol,
    err["E7", cs, <|"Probe" -> t0, "MaxRelativeResidual" -> N[maxRel, 6],
      "ResidTol" -> rtol, "Detail" -> "ODE residual check failed"|>]];
  maxRel];

End[];
EndPackage[];
