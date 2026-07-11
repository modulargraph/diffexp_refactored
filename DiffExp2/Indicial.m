(* DiffExp2/Indicial.m — exact sector-tag derivation at a chart center.
   Exact algebra only: no numerics, no tolerances, no fallbacks. *)

BeginPackage["DiffExp2`Indicial`", {"DiffExp2`Tolerances`", "DiffExp2`Config`"}];

ChartIndicial::usage = "ChartIndicial[A, t, eps, chartRef] classifies the chart-local system f' = A.f exactly: pole data, rank reduction, an epsilon-affine spectrum, Jordan chains, and resonance families. Returns IndicialData.";
MatrixPoleData::usage = "MatrixPoleData[A, t] gives <|\"PoleOrder\" -> r, \"Coefficients\" -> <|-r -> A_-r, ..., -1 -> A_-1|>|> exactly.";
FuchsianReduce::usage = "FuchsianReduce[A, t, eps, chartRef] reduces a pole order r >= 2 system to Fuchsian form by an exact Moser/shearing gauge transform. Returns ReductionData.";
EpsDegenerateFamilies::usage = "EpsDegenerateFamilies[indicialData] selects the families whose eps -> 0 eigenvectors collide (EpsZeroDegeneracy > 0), for Transport's RecombineBasis.";

Begin["`Private`"];

err[id_, chartRef_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Indicial",
    "Chart" -> Lookup[chartRef, "Name", "(unnamed)"],
    "Center" -> Lookup[chartRef, "Center", None],
    "Variable" -> Lookup[chartRef, "Variable", None]|>, payload]];

(* exact zero test on rational/algebraic scalars; no banned tokens.
   Together-canonical rational functions over the Gaussian rationals in
   non-numeric symbols are a zero-DECISION domain: === 0 decides.  Only
   forms outside it (radicals/Root, numeric symbol constants) need the
   exact RootReduce test; both branches stay exact.  zeroCanQ takes input
   the caller has already Together/Cancel-canonicalized. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && RootReduce[c] === 0);
zeroQ[e_] := zeroCanQ[Together[e]];
matZeroQ[m_] := AllTrue[Flatten[m], zeroQ];

(* ---- exact t-valuation and Laurent heads of rational entries ---- *)

(* min t-degree of a polynomial; entry must be polynomial in t.
   contract: p is 0 or a polynomial part of a CANCELED fraction, so the
   syntactic test in zeroCanQ is complete *)
polyMinDeg[p_, t_] := If[zeroCanQ[p], Infinity, Exponent[p, t, Min]];
tValuation[e_, t_] := Module[{c = Cancel[Together[e]]},
  If[zeroCanQ[c], Infinity,
    polyMinDeg[Numerator[c], t] - polyMinDeg[Denominator[c], t]]];
matMinValuation[m_, t_] := Min[tValuation[#, t] & /@ Flatten[m]];

(* exact Laurent coefficient of a rational function at t = 0, order k:
   e = t^v num1/den1 with num1(0), den1(0) != 0; coefficients by the
   exact division recursion (never Series on possibly-inexact objects) *)
laurentCoeff[e_, t_, k_Integer] := Module[{c, v, num, den, nc, dc, ord, cs},
  c = Cancel[Together[e]];
  If[zeroCanQ[c], Return[0]];
  v = tValuation[c, t];
  If[k < v, Return[0]];
  num = Numerator[c]; den = Denominator[c];
  num = Cancel[num/t^polyMinDeg[num, t]];
  den = Cancel[den/t^polyMinDeg[den, t]];
  ord = k - v;
  nc = Table[Coefficient[num, t, j], {j, 0, ord}];
  dc = Table[Coefficient[den, t, j], {j, 0, ord}];
  cs = Table[0, {ord + 1}];
  cs[[1]] = Cancel[Together[nc[[1]]/dc[[1]]]];
  Do[cs[[m + 1]] = Cancel[Together[
      (nc[[m + 1]] - Sum[dc[[j + 1]]*cs[[m - j + 1]], {j, 1, m}])/dc[[1]]]],
    {m, 1, ord}];
  cs[[ord + 1]]];

(* canonicalize once: tValuation and laurentCoeff each re-Cancel their input *)
leadingCoeff[e_, t_] := Module[{c = Cancel[Together[e]]},
  laurentCoeff[c, t, tValuation[c, t]]];

MatrixPoleData[A_?MatrixQ, t_Symbol] := Module[{r, heads},
  r = Max[0, Max[-matMinValuation[A, t], 0]];
  (* decrement past identically-zero candidate heads *)
  While[r > 0 && matZeroQ[Map[laurentCoeff[#, t, -r] &, A, {2}]], r--];
  heads = If[r == 0, <||>,
    Association[Table[(-k) -> Map[laurentCoeff[#, t, -k] &, A, {2}], {k, r, 1, -1}]]];
  <|"PoleOrder" -> r, "Coefficients" -> heads|>];

(* ---- FuchsianReduce: exact Moser/shearing rank reduction ---- *)

connectionMatrix[M_, T_, t_] :=
  Map[Cancel[Together[#]] &, Inverse[T] . (M . T - t*D[T, t]), {2}];

FuchsianReduce[A_?MatrixQ, t_Symbol, eps_Symbol, chartRef_Association] := Module[
  {d = Length[A], M, v, L, T, B, steps = 0, maxSteps = 200, trimmed = False,
   j, u, vals, m, cvec, pivot, Emat, uNew, U, changed, pass, TInv, residue,
   maxDeg, r},
  M = Map[Cancel[Together[#]] &, t*A, {2}];
  v = matMinValuation[M, t];
  r = 1 - v;
  (* Step 1: Moser pre-check — leading theta coefficient must be nilpotent *)
  L = Map[laurentCoeff[#, t, v] &, M, {2}];
  If[!matZeroQ[MatrixPower[L, d]],
    err["E3", chartRef, <|"PoleOrder" -> r, "LeadingCoefficient" -> L,
      "NonNilpotencyWitness" -> Map[If[zeroQ[#], 0, 1] &, MatrixPower[L, d], {2}],
      "Detail" -> "irregular singular point: leading Laurent matrix is not nilpotent (Moser condition)"|>]];
  (* Step 2: lattice saturation *)
  T = IdentityMatrix[d];
  While[True,
    B = connectionMatrix[M, T, t];
    If[matMinValuation[B, t] >= 0, Break[]];
    If[steps >= maxSteps,
      err["E4", chartRef, <|"MaxSteps" -> maxSteps,
        "FinalMinValuation" -> matMinValuation[B, t],
        "GaugeMaxDegree" -> Max[Exponent[#, t] & /@ Flatten[Together /@ Flatten[T]]],
        "PoleOrder" -> r,
        "Detail" -> "rank reduction did not terminate"|>]];
    (* column with most negative min valuation; ties -> smallest j *)
    j = First[Ordering[Table[matMinValuation[{B[[All, c]]}, t], {c, d}], 1]];
    u = B[[All, j]];
    vals = tValuation[#, t] & /@ u;
    m = Min[vals];
    cvec = Table[If[vals[[i]] === m, leadingCoeff[u[[i]], t], 0], {i, d}];
    pivot = SelectFirst[Range[d], !zeroQ[cvec[[#]]] &, None];
    If[pivot === None,
      err["E6", chartRef, <|"Stage" -> "Pivot", "Column" -> j,
        "Detail" -> "no nonzero leading coefficient in adjoined lattice vector"|>]];
    Emat = IdentityMatrix[d];
    Emat[[pivot, pivot]] = 1/cvec[[pivot]];
    Do[If[i =!= pivot, Emat[[i, pivot]] = -cvec[[i]]/cvec[[pivot]]], {i, d}];
    uNew = Cancel[Together[#]] & /@ (Emat . u);
    U = Transpose[ReplacePart[Transpose[IdentityMatrix[d]], pivot -> uNew]];
    T = Map[Cancel[Together[#]] &, (T . Inverse[Emat]) . U, {2}];
    steps++];
  (* Step 3: trim pass — ONLY columns that actually carry a pole (N-1b) *)
  Do[
    changed = False;
    Do[
      If[matMinValuation[{T[[All, p]]}, t] < 0,
        Module[{Tp = T, Bp},
          Tp[[All, p]] = t*Tp[[All, p]];
          Bp = connectionMatrix[M, Map[Cancel[Together[#]] &, Tp, {2}], t];
          If[matMinValuation[Bp, t] >= 0,
            T = Map[Cancel[Together[#]] &, Tp, {2}]; trimmed = True; changed = True]]],
      {p, d}];
    If[!changed, Break[]],
    {pass, 50}];
  (* Step 4: degree guard *)
  maxDeg = Max[Flatten[Map[
    {Exponent[Numerator[Together[#]], t], Exponent[Denominator[Together[#]], t]} &,
    Flatten[T]]]];
  If[maxDeg > 4 d + 8,
    err["E5", chartRef, <|"MaxDegree" -> maxDeg, "Bound" -> 4 d + 8,
      "Dimension" -> d, "Steps" -> steps,
      "Detail" -> "degenerate gauge lattice (degree guard)"|>]];
  (* Step 5: finalize *)
  TInv = Map[Cancel[Together[#]] &, Inverse[T], {2}];
  If[!matZeroQ[T . TInv - IdentityMatrix[d]],
    err["E6", chartRef, <|"Stage" -> "GaugeInverse",
      "Detail" -> "T.TInv != Identity"|>]];
  B = connectionMatrix[M, T, t];
  If[matMinValuation[B, t] < 0,
    err["E6", chartRef, <|"Stage" -> "GaugeInverse",
      "Detail" -> "final theta matrix not holomorphic"|>]];
  residue = Map[Cancel[Together[#]] &, B /. t -> 0, {2}];
  <|"PoleOrder" -> r, "Gauge" -> T, "GaugeInverse" -> TInv,
    "ThetaMatrix" -> B, "Residue" -> residue, "Steps" -> steps,
    "Trimmed" -> trimmed,
    "GaugeValuation" -> matMinValuation[T, t]|>];

(* ---- exact epsilon-affine spectrum ---- *)

AffineSpectrum[R_?MatrixQ, eps_Symbol, chartRef_Association] := Module[
  {lam, chi, coeffs, lcm, chiC, facs, recs = {}, roots, b, a, merged, lc, closure},
  chi = CharacteristicPolynomial[R, lam];
  coeffs = CoefficientList[Together[chi], lam];
  lcm = PolynomialLCM @@ (Denominator[Together[#]] & /@ coeffs);
  chiC = Together[chi*lcm];
  facs = FactorList[chiC, Extension -> Automatic];
  Do[Module[{fct = f[[1]], e = f[[2]], degL},
    degL = Exponent[fct, lam];
    Which[
      degL === 0, Null,
      degL === 1,
      AppendTo[recs, {Together[-Coefficient[fct, lam, 0]/Coefficient[fct, lam, 1]], e}],
      True,
      roots = lam /. Solve[fct == 0, lam];
      If[!ListQ[roots] || Length[roots] =!= degL,
        err["E2", chartRef, <|"CharPoly" -> chiC, "Factor" -> fct,
          "Detail" -> "no exact root extraction for irreducible factor of degree >= 2"|>]];
      Do[AppendTo[recs, {r0, e}], {r0, roots}]]],
    {f, facs}];
  (* affine verification per root *)
  recs = Map[Module[{r0 = #[[1]], e = #[[2]], bb, aa},
    bb = Together[D[r0, eps]];
    If[!zeroQ[D[bb, eps]],
      err["E2", chartRef, <|"CharPoly" -> chiC, "Root" -> r0,
        "Detail" -> "indicial eigenvalue is not affine in epsilon; nonlinear epsilon exponents are unsupported"|>]];
    aa = Together[r0 /. eps -> 0];
    If[!zeroQ[r0 - (aa + bb*eps)],
      err["E2", chartRef, <|"CharPoly" -> chiC, "Root" -> r0,
        "Detail" -> "affine reconstruction failed"|>]];
    {RootReduce[aa], RootReduce[bb], e}] &, recs];
  (* merge identical (a, b) *)
  merged = Map[<|"a" -> #[[1, 1]], "b" -> #[[1, 2]],
      "Multiplicity" -> Total[#[[All, 3]]]|> &,
    GatherBy[recs, {#[[1]], #[[2]]} &]];
  (* closure assert I-3 *)
  lc = Last[CoefficientList[chiC, lam]];
  closure = Together[lc*Product[(lam - m["a"] - m["b"]*eps)^m["Multiplicity"],
      {m, merged}] - chiC];
  If[!zeroQ[closure] || Total[merged[[All, "Multiplicity"]]] =!= Length[R],
    err["E6", chartRef, <|"Stage" -> "CharPolyClosure",
      "Detail" -> "product of affine factors does not reproduce the characteristic polynomial"|>]];
  SortBy[merged, {#["a"], #["b"]} &]];

rnk[{}] := 0;
rnk[m_List] := MatrixRank[m];

(* ---- JordanChains: exact confluence detection ---- *)

JordanChains[R_?MatrixQ, spectrum_List, eps_Symbol, chartRef_Association] := Module[
  {d = Length[R]},
  Map[Module[{lam = #["a"] + #["b"]*eps, m = #["Multiplicity"], Mk, nulls, nullDims,
      kmax, blocks, chains = {}, claimed, lvl, want, basis, cand, tops, vrec = #},
    Mk[j_] := Mk[j] = Map[Cancel[Together[#]] &, MatrixPower[R - lam*IdentityMatrix[d], j], {2}];
    nullDims = {0};
    kmax = 0;
    While[Last[nullDims] < m,
      kmax++;
      nulls = NullSpace[Mk[kmax]];
      AppendTo[nullDims, Length[nulls]];
      If[kmax > m || Length[nullDims] >= 2 &&
          nullDims[[-1]] < nullDims[[-2]],
        err["E6", chartRef, <|"Stage" -> "Nullity", "Eigenvalue" -> {vrec["a"], vrec["b"]},
          "NullitySequence" -> Rest[nullDims],
          "Detail" -> "inconsistent nested nullity sequence"|>]]];
    (* blocks of size >= j: nullDims[[j+1]] - nullDims[[j]] *)
    blocks = Reverse[Sort[Flatten[Table[
      ConstantArray[j, (nullDims[[j + 1]] - nullDims[[j]]) -
        If[j + 2 <= Length[nullDims], nullDims[[j + 2]] - nullDims[[j + 1]], 0]],
      {j, 1, kmax}]]]];
    If[Total[blocks] =!= m,
      err["E6", chartRef, <|"Stage" -> "Nullity", "Eigenvalue" -> {vrec["a"], vrec["b"]},
        "BlockSizes" -> blocks, "Detail" -> "block sizes do not sum to multiplicity"|>]];
    (* chains: pick tops level by level, longest first *)
    claimed = {};
    Do[
      want = Count[blocks, lvl];
      If[want > 0,
        basis = NullSpace[Mk[lvl]];
        (* independence vs Null^(lvl-1) AND the level-lvl members of longer chains *)
        cand = Join[If[lvl > 1, NullSpace[Mk[lvl - 1]], {}],
          Map[#[[1, lvl]] &, Select[claimed, Length[#[[1]]] >= lvl &]]];
        tops = {};
        Do[
          If[Length[tops] < want &&
              rnk[Join[cand, tops, {bv}]] > rnk[Join[cand, tops]],
            AppendTo[tops, bv]],
          {bv, basis}];
        If[Length[tops] < want,
          err["E6", chartRef, <|"Stage" -> "Chain", "Eigenvalue" -> {vrec["a"], vrec["b"]},
            "Detail" -> "could not select independent chain tops"|>]];
        Do[Module[{top, chain, fnz},
          fnz = SelectFirst[top0, !zeroQ[#] &, 1];
          top = Map[Cancel[Together[#/fnz]] &, top0];
          chain = NestList[Map[Cancel[Together[#]] &, Mk[1] . #] &, top, lvl - 1];
          chain = Reverse[chain];  (* eigenvector FIRST *)
          If[!matZeroQ[{Mk[1] . First[chain]}] ||
             !AllTrue[Range[2, lvl], matZeroQ[{Mk[1] . chain[[#]] - chain[[# - 1]]}] &],
            err["E6", chartRef, <|"Stage" -> "Chain",
              "Eigenvalue" -> {vrec["a"], vrec["b"]},
              "Detail" -> "chain identities failed"|>]];
          AppendTo[chains, chain];
          AppendTo[claimed, {chain, top}]],
          {top0, tops}]],
      {lvl, kmax, 1, -1}];
    Join[vrec, <|"BlockSizes" -> blocks, "Chains" -> chains|>]] &, spectrum]];

(* ---- PartitionResonanceFamilies ---- *)

PartitionResonanceFamilies[spectrum_List, eps_Symbol, chartRef_Association] := Module[
  {n = Length[spectrum], adj, comps, fams},
  adj = Table[IntegerQ[RootReduce[spectrum[[i]]["a"] - spectrum[[j]]["a"]]], {i, n}, {j, n}];
  comps = ConnectedComponents[AdjacencyGraph[Boole[adj]]];
  fams = Map[Module[{members = SortBy[spectrum[[#]], {#["a"], #["b"]} &], bs, sectors,
      collisions = {}, class, joint, logmax, r0},
    bs = DeleteDuplicates[RootReduce /@ members[[All, "b"]]];
    (* sector specs: per member, per block, per chain position *)
    sectors = Flatten[Map[Function[mem,
      Flatten[Map[Function[bsz, Table[
        <|"a" -> mem["a"], "b" -> mem["b"],
          "p" -> q + Total[Map[If[zeroQ[#["b"] - mem["b"]] &&
              IntegerQ[RootReduce[#["a"] - mem["a"]]] &&
              TrueQ[RootReduce[#["a"] - mem["a"]] > 0],
            #["Multiplicity"], 0] &, members]]|>,
        {q, 0, bsz - 1}]], mem["BlockSizes"]], 1]], members], 1];
    sectors = SortBy[sectors, {#["a"], #["b"], #["p"]} &];
    (* collisions *)
    Do[Module[{ai = members[[i]]["a"], bi = members[[i]]["b"],
        aj = members[[j]]["a"], bj = members[[j]]["b"], nn},
      If[i =!= j,
        nn = RootReduce[aj - ai];
        If[IntegerQ[nn] && nn >= 0 && !(zeroQ[ai - aj] && zeroQ[bi - bj]) &&
            (nn > 0 || i < j),
          AppendTo[collisions, <|"n" -> nn, "LowerIdx" -> i, "UpperIdx" -> j,
            "DeltaB" -> RootReduce[bi - bj],
            "Type" -> If[zeroQ[bi - bj], "Log", "LaurentShift"]|>]]]],
      {i, Length[members]}, {j, Length[members]}];
    joint = Length[bs] >= 2;
    class = Which[
      joint, "Pseudo",
      Length[members] >= 2, "TrueResonant",
      Max[members[[1]]["BlockSizes"]] >= 2, "Confluent",
      True, "Single"];
    logmax = Max[sectors[[All, "p"]]];
    (* EpsZeroDegeneracy (2.6 step 7): chain tops, eps-valuation normalized, at eps = 0 *)
    r0 = If[!joint, 0,
      Module[{tops, cols},
        tops = Flatten[Map[Function[mem, Last /@ mem["Chains"]], members], 1];
        cols = Map[Function[v, Module[{ev},
          ev = Min[Map[Module[{c = Cancel[Together[#]]},
            If[zeroCanQ[c], Infinity,
              polyMinDeg[Numerator[c], eps] - polyMinDeg[Denominator[c], eps]]] &, v]];
          Map[Cancel[Together[#/eps^ev]] /. eps -> 0 &, v]]], tops];
        Length[cols] - MatrixRank[cols]]];
    <|"Sectors" -> sectors, "Members" -> members, "Class" -> class,
      "JointSolve" -> joint, "LogMax" -> logmax,
      "EpsZeroDegeneracy" -> r0, "Collisions" -> collisions|>] &, comps];
  SortBy[fams, {First[SortBy[#["Members"][[All, "a"]], Identity]],
    First[SortBy[#["Members"][[All, "b"]], Identity]]} &]];

(* ---- orchestrator ---- *)

validateExact[A_, chartRef_] := Module[{pos},
  pos = Position[A, _?InexactNumberQ, {0, Infinity}, 1];
  If[pos =!= {},
    err["E1", chartRef, <|"Position" -> First[pos],
      "Entry" -> Extract[A, Take[First[pos], Min[2, Length[First[pos]]]]],
      "Detail" -> "inexact entry; exact eps-rational input required (the d<var>_full.m exact export, not eps-truncated slices)"|>]];
  pos = Position[A, _SeriesData, {0, Infinity}, 1];
  If[pos =!= {},
    err["E1", chartRef, <|"Position" -> First[pos],
      "Detail" -> "SeriesData entry: epsilon-truncated slice exports cannot certify the exact indicial spectrum; use the exact full export d<var>_full.m (ExportGeneralMatrix)"|>]];];

ChartIndicial[A_?MatrixQ, t_Symbol, eps_Symbol, chartRef_Association] := Module[
  {d = Length[A], pole, red, R, spec, fams, idMat, regular, Anorm},
  If[d == 0 || Length[A] =!= Length[First[A]] || t === eps ||
      !KeyExistsQ[chartRef, "Name"],
    err["E7", chartRef, <|"Shape" -> Dimensions[A], "t" -> t, "eps" -> eps,
      "Detail" -> "bad shape or chart reference"|>]];
  validateExact[A, chartRef];
  (* canonicalize entries ONCE: MatrixPoleData re-walks every entry per
     Laurent order (valuation + laurentCoeff), each walk re-Cancels *)
  Anorm = Map[Cancel[Together[#]] &, A, {2}];
  pole = MatrixPoleData[Anorm, t];
  idMat = IdentityMatrix[d];
  regular = pole["PoleOrder"] == 0;
  Which[
    regular,
    red = <|"PoleOrder" -> 0, "Gauge" -> idMat, "GaugeInverse" -> idMat,
      "ThetaMatrix" -> Map[Cancel[Together[#]] &, t*Anorm, {2}],
      "Residue" -> ConstantArray[0, {d, d}], "Steps" -> 0,
      "Trimmed" -> False, "GaugeValuation" -> 0|>,
    pole["PoleOrder"] == 1,
    red = <|"PoleOrder" -> 1, "Gauge" -> idMat, "GaugeInverse" -> idMat,
      "ThetaMatrix" -> Map[Cancel[Together[#]] &, t*Anorm, {2}],
      "Residue" -> pole["Coefficients"][-1], "Steps" -> 0,
      "Trimmed" -> False, "GaugeValuation" -> 0|>,
    True,
    red = FuchsianReduce[Anorm, t, eps, chartRef]];
  R = red["Residue"];
  If[regular,
    (* one (0,0,0) family with d-dimensional coefficient space *)
    spec = {<|"a" -> 0, "b" -> 0, "Multiplicity" -> d,
      "BlockSizes" -> ConstantArray[1, d],
      "Chains" -> Table[{idMat[[i]]}, {i, d}]|>};
    fams = {<|"Sectors" -> ConstantArray[<|"a" -> 0, "b" -> 0, "p" -> 0|>, d],
      "Members" -> spec, "Class" -> "Single", "JointSolve" -> False,
      "LogMax" -> 0, "EpsZeroDegeneracy" -> 0, "Collisions" -> {}|>},
    spec = AffineSpectrum[R, eps, chartRef];
    spec = JordanChains[R, spec, eps, chartRef];
    fams = PartitionResonanceFamilies[spec, eps, chartRef]];
  (* I-2 invariants *)
  If[Total[spec[[All, "Multiplicity"]]] =!= d ||
     Total[Length[#["Sectors"]] & /@ fams] =!= d ||
     !AllTrue[spec, Total[#["BlockSizes"]] === #["Multiplicity"] &],
    err["E6", chartRef, <|"Stage" -> "Nullity",
      "Detail" -> "multiplicity/sector counts do not sum to the dimension"|>]];
  <|"Chart" -> chartRef, "Dimension" -> d, "PoleData" -> pole,
    "Reduction" -> red, "Residue" -> R, "Spectrum" -> spec,
    "Families" -> fams, "Regular" -> regular|>];

EpsDegenerateFamilies[data_Association] :=
  Map[<|"FamilyIndex" -> #, "EpsZeroDegeneracy" -> data["Families"][[#]]["EpsZeroDegeneracy"]|> &,
    Select[Range[Length[data["Families"]]],
      data["Families"][[#]]["EpsZeroDegeneracy"] > 0 &]];

End[];
EndPackage[];
