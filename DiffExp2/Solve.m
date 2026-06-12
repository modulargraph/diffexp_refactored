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
zeroQ[e_] := TrueQ[PossibleZeroQ[Together[e]]];

esNew = DiffExp2`EpsSeries`ESNew; esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd; esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes; esDiv = DiffExp2`EpsSeries`ESDivide;
esShift = DiffExp2`EpsSeries`ESShift; esCoeff = DiffExp2`EpsSeries`ESCoefficient;
esMin = DiffExp2`EpsSeries`ESMinPower; esCM = DiffExp2`EpsSeries`ESCompleteMax;
esFrom = DiffExp2`EpsSeries`ESFromExpression; esTrim = DiffExp2`EpsSeries`ESTrim;

(* exact t-Laurent coefficient of a rational function at t = 0 (local copy
   of the Indicial recursion; entries exact) *)
polyMinDeg[p_, t_] := If[zeroQ[p], Infinity, Exponent[p, t, Min]];
tVal[e_, t_] := Module[{c = Cancel[Together[e]]},
  If[zeroQ[c], Infinity,
    polyMinDeg[Numerator[c], t] - polyMinDeg[Denominator[c], t]]];
tLaurent[e_, t_, k_Integer] := Module[{c, v, num, den, nc, dc, ord, csr},
  c = Cancel[Together[e]];
  If[zeroQ[c] || k < (v = tVal[c, t]), Return[0]];
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

PrepareChart[sys_Association, chart_Association] := Module[
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
  If[zeroQ[detV], err["E2", chart, <|"Detail" -> "spectral frame V is singular"|>]];
  VInv = Map[Cancel[Together[#]] &, Inverse[V], {2}];
  <|"ChartVar" -> t, "Center" -> x0,
    "ChartMap" -> <|"Center" -> x0, "Scale" -> beta|>,
    "Radius" -> chart["Radius"],
    "Prescriptions" -> Lookup[chart, "Prescriptions", {}],
    "SystemSize" -> d,
    "ThetaMatrix" -> idata["Reduction"]["ThetaMatrix"],
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

(* prepared eps-expansions of the cleared system, shared per solve call *)
prepareCleared[cs_, wideTop_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], d = cs["SystemSize"],
   den, num, d0, dD, dN, dES, NhatES},
  den = Together[PolynomialLCM @@ (Denominator[Together[#]] & /@ Flatten[cs["ThetaMatrix"]])];
  num = Map[Cancel[Together[#*den]] &, cs["ThetaMatrix"], {2}];
  d0 = Together[den /. t -> 0];
  If[zeroQ[d0 /. eps -> 0],
    err["E3", cs, <|"Denominator" -> den,
      "Detail" -> "cleared denominator degenerates onto the chart origin at eps = 0"|>]];
  dD = Exponent[den, t];
  dN = Max[0, Max[Exponent[#, t] & /@ Flatten[num]]];
  (* R6: numericize the recursion coefficient expansions at 2x WP -
     structure (tags, windows, case selection) is decided on EXACT data
     upstream; exact-rational coefficient growth grinds multi-master
     systems.  Exact zeros stay exact. *)
  Module[{wp2 = 2*DiffExp2`Config`CFG["WorkingPrecision"], numz},
    numz[s_] := DiffExp2`EpsSeries`ESMap[
      If[# === 0 || ByteCount[#] <= 500, #, N[#, wp2]] &, s];
    dES = Table[numz[esFrom[Coefficient[den, t, j], eps, wideTop]], {j, 0, dD}];
    NhatES = Table[Module[{Nj, Nhat},
      Nj = Map[Coefficient[#, t, j] &, num, {2}];
      Nhat = Map[Cancel[Together[#]] &, cs["VInv"] . Nj . cs["V"], {2}];
      Map[If[zeroQ[#], esZero[wideTop], numz[esFrom[#, eps, wideTop]]] &, Nhat, {2}]],
      {j, 0, dN}]];
  <|"dES" -> dES, "NhatES" -> NhatES, "dD" -> dD, "dN" -> dN|>];

(* matrix(of ES).vector(of ES) *)
mDotV[mES_, v_] := Table[
  Module[{s = None},
    Do[Module[{term = esTimes[mES[[r, c]], v[[c]]]},
      s = If[s === None, term, esAdd[s, term]]],
      {c, Length[v]}];
    s], {r, Length[mES]}];

(* CASE T/P block solve: uhat_i = (1/d0) Sum_{m=0}^{q-1} N^m rhs / delta^(m+1);
   (N^m rhs)_r = rhs_{r+m} *)
blockSolveTP[rhs_, deltaES_, d0ES_, q_] := Module[{inv = DiffExp2`EpsSeries`ESInvert[deltaES]},
  Table[Module[{acc = None, pw = inv},
    Do[
      If[r + m <= q,
        acc = If[acc === None, esTimes[pw, rhs[[r + m]]],
          esAdd[acc, esTimes[pw, rhs[[r + m]]]]]];
      pw = esTimes[pw, inv],
      {m, 0, q - 1}];
    esDiv[acc, d0ES]], {r, q}]];

(* runRecursion: theta g = B g + s in the J-frame for tag lambda = a + b eps.
   srcHat: None | function srcHat[n, l] -> d-vector of EpsSeries (theta-form,
   J-frame).  init: None | initial U[0][l] table.  Returns
   <|"U", "Hits", "P"|> with U[[n+1, l+1]] a d-vector of EpsSeries. *)
runRecursion[cs_, prep_, aT_, bT_, P_, nmax_, srcHat_, wideTop_, init_] := Module[
  {d = cs["SystemSize"], dES = prep["dES"], NhatES = prep["NhatES"],
   dD = prep["dD"], dN = prep["dN"], d0ES, blocks, U, hits = {}, lamPoly, n0},
  d0ES = dES[[1]];
  blocks = blockList[cs];
  lamPoly[m_] := esNew[0, PadRight[{Together[aT + m], Together[bT]}, wideTop + 1]];
  U = Table[Table[Table[esZero[wideTop], {d}], {l, 0, P + 1}], {n, 0, nmax}];
  n0 = If[init === None, 0, 1];
  If[init =!= None,
    Do[U[[1, l + 1]] = init[[l + 1]], {l, 0, Min[P, Length[init] - 1]}]];
  Do[Module[{R, rhsFull},
    (* R[l] from HISTORY (t-orders < n) + source *)
    R = Table[Module[{acc = Table[esZero[wideTop], {d}]},
      Do[If[n - j >= 0,
        acc = vAddLocal[acc, mDotV[NhatES[[j + 1]], U[[n - j + 1, l + 1]]]]],
        {j, 1, Min[n, dN]}];
      Do[If[n - j >= 0,
        acc = vAddLocal[acc, Table[esScale[-1, esAdd[
          esTimes[dES[[j + 1]], esTimes[lamPoly[n - j], U[[n - j + 1, l + 1, r]]]],
          esTimes[dES[[j + 1]], esShift[U[[n - j + 1, l + 2, r]], 1]]]], {r, d}]]],
        {j, 1, Min[n, dD]}];
      If[srcHat =!= None,
        Do[If[n - j >= 0,
          Module[{sv = srcHat[n - j, l]},
            If[sv =!= None,
              acc = vAddLocal[acc, Table[esTimes[dES[[j + 1]], sv[[r]]], {r, d}]]]]],
          {j, 0, Min[n, dD]}]];
      acc], {l, 0, P}];
    (* solve top-down in l: L_n U[n,l] = R[l] - d0 eps U[n,l+1] *)
    Do[Module[{},
      rhsFull = Table[esAdd[R[[l + 1, r]],
        esScale[-1, esTimes[d0ES, esShift[U[[n + 1, l + 2, r]], 1]]]], {r, d}];
      Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"],
          dA, dB},
        dA = Together[aT + n - aI]; dB = Together[bT - bI];
        Which[
          !zeroQ[dA],
          U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTP[rhsFull[[colsB]],
              esTrim[esNew[0, PadRight[{dA, dB}, wideTop + 1]]], d0ES, q]]],
          !zeroQ[dB],
          (U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
            Thread[colsB -> blockSolveTP[rhsFull[[colsB]],
              esNew[1, PadRight[{dB}, wideTop]], d0ES, q]]];
           If[l == 0, AppendTo[hits,
             <|"n" -> n, "Cols" -> colsB, "DeltaB" -> dB|>]]),
          True, Null]],  (* CASE R after the level loop *)
        {blk, blocks}]],
      {l, P, 0, -1}];
    (* CASE R ladder *)
    Do[Module[{aI = blk["a"], bI = blk["b"], q = blk["q"], colsB = blk["Cols"]},
      If[zeroQ[aT + n - aI] && zeroQ[bT - bI],
        Module[{Rt, assigned = Association[]},
          Rt = Table[Table[esDiv[R[[l + 1, colsB[[r]]]], d0ES], {r, q}], {l, 0, P}];
          (* row q, level l: U[n,l+1]_q = Rt[l]_q / eps  (the log bump) *)
          Do[If[l + 1 <= P,
            assigned[{l + 1, q}] =
              esDiv[Rt[[l + 1, q]], esNew[1, PadRight[{1}, wideTop]]]],
            {l, 0, P - 1}];
          (* rows r<q top-down: U[n,l]_{r+1} = eps U[n,l+1]_r - Rt[l]_r *)
          Do[
            Do[Module[{upr = Lookup[assigned, Key[{l + 1, r}], esZero[wideTop]]},
              If[!KeyExistsQ[assigned, {l, r + 1}],
                assigned[{l, r + 1}] =
                  esAdd[esShift[upr, 1], esScale[-1, Rt[[l + 1, r]]]]]],
              {r, 1, q - 1}],
            {l, P, 0, -1}];
          Do[If[KeyExistsQ[assigned, {l, r}],
            U[[n + 1, l + 1]] = ReplacePart[U[[n + 1, l + 1]],
              colsB[[r]] -> assigned[{l, r}]]],
            {l, 0, P}, {r, 1, q}]]]],
      {blk, blocks}]],
    {n, n0, nmax}];
  <|"U" -> U, "Hits" -> hits, "P" -> P|>];

vAddLocal[u_, v_] := MapThread[esAdd, {u, v}];

(* ---- assembly ---- *)

(* J-frame U -> original-frame LocalSolution with sectors (a, b, l) *)
assembleSolution[cs_, aT_, bT_, rec_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], d = cs["SystemSize"], U = rec["U"],
   P = rec["P"], kminO, kmaxO, Vexp, W, secs, gv, ls},
  (* original-frame vectors W[n+1][l+1] = V . U *)
  Module[{top = Max[esCM /@ Flatten[U, 2]],
    wp2 = 2*DiffExp2`Config`CFG["WorkingPrecision"]},
    Vexp = Map[If[zeroQ[#], esZero[top],
      DiffExp2`EpsSeries`ESMap[If[# === 0 || ByteCount[#] <= 500, #, N[#, wp2]] &,
        esFrom[Together[#], eps, top]]] &, cs["V"], {2}]];
  W = Table[mDotV[Vexp, U[[n + 1, l + 1]]], {n, 0, nmax}, {l, 0, P}];
  kminO = Min[esMin /@ Flatten[W, 2]];
  kmaxO = Min[esCM /@ Flatten[W, 2]];
  secs = Table[<|"a" -> aT, "b" -> bT, "p" -> l,
    "Coeffs" -> Table[Table[Table[esCoeff[W[[n + 1, l + 1, r]], k], {r, d}],
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

(* gauge composition f = T.g: T exact rational, poles only at t = 0 *)
applyGauge[cs_, ls_, nmax_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], t = cs["ChartVar"], T = cs["Gauge"], d,
   gv, Texp, kmin, kmax, newSecs, top},
  If[T === IdentityMatrix[cs["SystemSize"]] || T === IdentityMatrix[Length[T]],
    Return[ls]];
  d = Length[T];
  gv = Min[tVal[#, t] & /@ Flatten[T]];
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  top = kmax;
  (* Texp[m][r][c] = EpsSeries of the t^m Laurent coefficient of T_{rc} *)
  Texp = Table[Map[Module[{lc = tLaurent[#, t, m]},
      If[zeroQ[lc], esZero[top], esFrom[Together[lc], eps, top]]] &, T, {2}],
    {m, gv, nmax}];
  newSecs = Map[Module[{arr = #["Coeffs"], aS = #["a"], out, kminS, kmaxS, esArr},
    (* per (n, comp): EpsSeries from the rows *)
    esArr = Table[Table[esNew[kmin, Table[arr[[k - kmin + 1, n + 1, c]],
        {k, kmin, kmax}]], {c, d}], {n, 0, nmax}];
    out = Table[Table[Module[{s = None},
      Do[Module[{n = np - (m - gv)},  (* m + n = np + gv *)
        If[0 <= n <= nmax,
          Do[Module[{term = esTimes[Texp[[m - gv + 1, r, c]], esArr[[n + 1, c]]]},
            s = If[s === None, term, esAdd[s, term]]],
            {c, d}]]],
        {m, gv, gv + np}];
      If[s === None, esZero[kmax], s]], {r, d}], {np, 0, nmax}];
    Module[{kminN = Min[esMin /@ Flatten[out]], kmaxN = Min[esCM /@ Flatten[out]]},
      <|"a" -> Together[aS + gv], "b" -> #["b"], "p" -> #["p"],
        "Coeffs" -> Table[Table[Table[esCoeff[out[[n + 1, r]], k], {r, d}],
          {n, 0, nmax}], {k, kminN, kmaxN}],
        "KMin" -> kminN, "KMax" -> kmaxN|>]] &, ls["Sectors"]];
  Module[{kminA = Min[#["KMin"] & /@ newSecs], kmaxA = Min[#["KMax"] & /@ newSecs]},
    newSecs = Map[Module[{sh = #["KMin"] - kminA, padTop = kmaxA - #["KMin"] + 1},
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
   wideTop, prep, Pmax, cdMax},
  nmax = req["TOrder"];
  reqMin = req["EpsWindow", "Min"]; reqMax = req["EpsWindow", "CompleteMax"];
  Pmax = Max[0, Max[Table[logCeiling[cs, b["a"], b["b"], b["q"] - 1], {b, blocks}]]];
  cdMax = Max[0, Max[#["CollisionDepth"] & /@ fams]];
  wideTop = reqMax + Pmax + cdMax + 4 - Min[0, reqMin - Pmax - 2];
  prep = prepareCleared[cs, wideTop];
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
        rec = runRecursion[cs, prep, root["a"], root["b"], P, nmax, None, wideTop, init];
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
    wideTop2 = srcMax + P + 4 - Min[0, srcMin - P - 2];
    prep2 = prepareCleared[cs, wideTop2];
    VInvExp = Map[If[zeroQ[#], esZero[wideTop2],
      esFrom[Together[#], eps, wideTop2]] &, cs["VInv"], {2}];
    (* J-frame source: bHat[n+1] = VInv . (source eps-vectors at t-order n) *)
    bHat = Table[mDotV[VInvExp,
      Table[esNew[srcMin, Table[arr[[k - srcMin + 1, n + 1, c]],
        {k, srcMin, srcMax}]], {c, d}]],
      {n, 0, Min[nmax, Length[arr[[1]]] - 1]}];
    srcFn = Function[{n, l},
      If[l === pS && n + 1 <= Length[bHat], bHat[[n + 1]], None]];
    rec = runRecursion[cs, prep2, aS, bS, P, nmax, srcFn, wideTop2, None];
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
    Module[{nmaxS, frac},
      nmaxS = Min[#["TWindow", "CompleteMax"] & /@ sols];
      frac = Min[1/4, Rationalize[N[(rtol/100)^(1/(nmaxS + 1))], 10^-3]];
      cs["Radius"]*frac*(2 + Mod[Hash[cs["Center"]], 5])/10],
    probe];
  Do[Module[{f, df, lhs, rhs, srcv, k1, k2, scale},
    f = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, t0, "UsePade" -> False];
    df = DiffExp2`SectorSeries`EvaluateLocalSolution[
      DiffExp2`SectorSeries`DifferentiateLocalSolution[ls], t0, "UsePade" -> False];
    k1 = Max[esMin[f["Value"]], esMin[df["Value"]]];
    k2 = Min[esCM[f["Value"]], esCM[df["Value"]]];
    Bt0 = Map[Together[# /. t -> t0] &, cs["ThetaMatrix"], {2}];
    Bt0 = Map[If[zeroQ[#], esZero[k2], esFrom[#, eps, k2]] &, Bt0, {2}];
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
        maxRel = Max[maxRel, Max[0, Sequence @@ (Abs[N[#, 20]] & /@ rv)]/sc]],
        {k, Max[k1, esMin[thetaF]], Min[k2, esCM[thetaF],
          Min[esCM /@ Bf]]}]]],
    {ls, sols}];
  If[maxRel > rtol,
    err["E7", cs, <|"Probe" -> t0, "MaxRelativeResidual" -> N[maxRel, 6],
      "ResidTol" -> rtol, "Detail" -> "ODE residual check failed"|>]];
  maxRel];

End[];
EndPackage[];
