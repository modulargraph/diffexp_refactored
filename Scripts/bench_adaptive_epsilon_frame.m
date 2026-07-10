(* Exact lower-Laurent support prototype for singular recurrence frames.

   This is deliberately OUTSIDE DiffExp2/Solve.m.  It measures the support
   that the complete coupled recurrence actually reaches at the banana L1
   endpoints, without allocating the scalar longest-path frame used by the
   production solver.  Coefficients are exact; a coefficient is discarded
   only after exact Cancel/Together simplification.  Each Taylor/log row is
   completed through all spectral blocks (including CASE R) before its
   support is recorded.

   This exact-rational E shadow is a proof prototype, not the production
   representation: it intentionally keeps full rational functions so that
   omitted-tail influence is certified.  The n=20 scaling run was stopped
   after 150 seconds; larger orders require explicit opt-in and are design
   experiments, not speed claims.

   Environment:
     AEF_ORDERS=2              report checkpoints (safe default)
     AEF_ALLOW_SLOW=1          required when Max[AEF_ORDERS] > 6
     AEF_GUARD_TOP=16          retained exact upper guard
     AEF_DELIVERED_WIDTH=1     K in the scalar proxy W=K+4 n+5

   Output is one compact JSON object per line, prefixed by AEF. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envOr[name_, def_] := Module[{v = Environment[name]}, If[StringQ[v], v, def]];
orders = Sort@DeleteDuplicates@Map[ToExpression,
  StringSplit[envOr["AEF_ORDERS", "2"], ","]];
If[orders === {} || !AllTrue[orders, IntegerQ[#] && # >= 1 &],
  Print["AEF invalid AEF_ORDERS"]; Quit[2]];
maxOrder = Max[orders];
If[maxOrder > 6 && envOr["AEF_ALLOW_SLOW", "0"] =!= "1",
  Print["AEF exact rational shadow above order 6 requires AEF_ALLOW_SLOW=1"];
  Quit[2]];
guardTop = ToExpression[envOr["AEF_GUARD_TOP", "16"]];
deliveredWidth = ToExpression[envOr["AEF_DELIVERED_WIDTH", "1"]];
If[!IntegerQ[guardTop] || guardTop < 4 ||
    !IntegerQ[deliveredWidth] || deliveredWidth < 1,
  Print["AEF invalid guard/delivered width"]; Quit[2]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[e_] := Quiet[Catch[e, "DiffExp2Error"]];
catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> Max[10, maxOrder],
  "EpsilonOrder" -> 0, "DivisionOrder" -> 3,
  "StepDivisionOrder" -> 3, "Variables" -> {}}]];

eps = DiffExp2`Config`CanonicalEps[];
t = Global`t;
fixture = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench", "banana_L1.m"}]];
var = fixture["Variable"];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fixture["Matrix"], "Variable" -> var|>]];
If[FailureQ[sys], Print["AEF system load failed: ", sys]; Quit[1]];

exactZeroQ[e_] := TrueQ[DiffExp2`Solve`Private`zeroQ[
  Cancel[Together[e]]]];
simp[e_] := Cancel[Together[e]];

(* A scalar exact Laurent frame.  Finite=True certifies exact zero above
   the stored support; otherwise coefficients above guardTop are unknown.
   Empty+Finite is therefore an identically-zero formal series. *)
zFrame[] := <|"C" -> <||>, "Finite" -> True, "E" -> 0|>;
monoFrame[k_Integer, c_] := If[exactZeroQ[c], zFrame[],
  <|"C" -> <|k -> simp[c]|>, "Finite" -> True,
    "E" -> simp[c eps^k]|>];
frameZeroAllQ[f_] := exactZeroQ[f["E"]];
frameMin[f_] := epsVal[f["E"]];

frameFromTerms[terms_List, finite_, expr_] := Module[{keys, vals, out, ex},
  ex = simp[expr];
  If[terms === {}, Return[<|"C" -> <||>, "Finite" -> TrueQ[finite],
      "E" -> ex|>]];
  keys = Union @@ (Keys[#] & /@ terms);
  vals = Table[simp[Total[Lookup[#, k, 0] & /@ terms]], {k, keys}];
  out = Association@MapThread[Rule, {keys, vals}];
  out = KeySelect[out, !exactZeroQ[out[#]] &];
  <|"C" -> out, "Finite" -> TrueQ[finite], "E" -> ex|>];

fAdd[fs_List] := frameFromTerms[Lookup[fs, "C"],
  And @@ Lookup[fs, "Finite"], Total[Lookup[fs, "E"]]];
fNeg[f_] := <|"C" -> Map[-# &, f["C"]], "Finite" -> f["Finite"],
  "E" -> -f["E"]|>;
fScaleShift[f_, c_, s_Integer, cap_:guardTop] :=
  If[exactZeroQ[c] || frameZeroAllQ[f], zFrame[],
    <|"C" -> Association@KeyValueMap[(#1 + s) -> simp[c #2] &,
        KeySelect[f["C"], # + s <= cap &]], "Finite" -> f["Finite"],
      "E" -> simp[c eps^s f["E"]]|>];

epsVal[e_] := Module[{c = simp[e], num, den},
  If[exactZeroQ[c], Return[Infinity]];
  num = Numerator[c]; den = Denominator[c];
  Exponent[num, eps, Min] - Exponent[den, eps, Min]];

laurentPolynomialQ[e_] := Module[{c = simp[e], v, den0},
  If[exactZeroQ[c], Return[True]];
  v = epsVal[c];
  den0 = simp[Denominator[c]/eps^Exponent[Denominator[c], eps, Min]];
  FreeQ[den0, eps]];

(* Exact rational Laurent coefficients, copied in small standalone form
   from ratEpsList but without its numericization gate. *)
exactLaurent[e_, hi_Integer] := Module[
  {c = simp[e], num, den, vn, vd, v, n0, d0, rel, nc, dc, cs},
  If[exactZeroQ[c], Return[<||>]];
  num = Numerator[c]; den = Denominator[c];
  vn = Exponent[num, eps, Min]; vd = Exponent[den, eps, Min]; v = vn - vd;
  If[v > hi, Return[<||>]];
  n0 = simp[num/eps^vn]; d0 = simp[den/eps^vd]; rel = hi - v;
  nc = Table[Coefficient[n0, eps, j], {j, 0, rel}];
  dc = Table[Coefficient[d0, eps, j], {j, 0, rel}];
  cs = ConstantArray[0, rel + 1];
  cs[[1]] = simp[nc[[1]]/dc[[1]]];
  Do[cs[[m + 1]] = simp[(nc[[m + 1]] -
      Sum[dc[[j + 1]] cs[[m - j + 1]], {j, 1, m}])/dc[[1]]],
    {m, 1, rel}];
  Association@Table[v + j -> cs[[j + 1]], {j, 0, rel}]];

scalarOp[e_, hi_Integer] := <|
  "C" -> exactLaurent[e, hi], "Finite" -> laurentPolynomialQ[e],
  "Min" -> epsVal[e], "E" -> simp[e]|>;

matrixOp[m_, hi_Integer] := Module[{d = Length[m], entries, keys, mats},
  entries = Map[exactLaurent[#, hi] &, m, {2}];
  keys = Union[Flatten[Map[Keys, entries, {2}]]];
  mats = Association@Table[k -> Table[Lookup[entries[[r, c]], k, 0],
      {r, d}, {c, d}], {k, keys}];
  <|"C" -> mats,
    "Finite" -> And @@ Flatten[Map[laurentPolynomialQ, m, {2}]],
    "E" -> m,
    "Min" -> Module[{v = Select[Flatten[Map[epsVal, m, {2}]], IntegerQ]},
      If[v === {}, Infinity, Min[v]]]|>];

fScalarOp[op_, f_, cap_Integer] := Module[{terms},
  If[frameZeroAllQ[f] || exactZeroQ[op["E"]], Return[zFrame[]]];
  terms = Flatten@Table[
    If[k + s <= cap, <|k + s -> simp[a b]|>, <||>],
    {s, Keys[op["C"]]}, {k, Keys[f["C"]]},
    {a, {op["C"][s]}}, {b, {f["C"][k]}}];
  frameFromTerms[terms, op["Finite"] && f["Finite"], op["E"] f["E"]]];

fMatrixOp[op_, vf_List, cap_Integer] := Module[
  {d = Length[vf], terms, rowFinite},
  Table[
    terms = Flatten@Table[Module[{a = op["C"][s][[r, c]]},
        If[exactZeroQ[a], {},
          KeyValueMap[If[#1 + s <= cap, <|#1 + s -> simp[a #2]|>, <||>] &,
            vf[[c, "C"]]]]],
      {s, Keys[op["C"]]}, {c, d}];
    rowFinite = op["Finite"] && And @@ Lookup[vf, "Finite"];
    frameFromTerms[Flatten[terms], rowFinite,
      Sum[op["E"][[r, c]] vf[[c, "E"]], {c, d}]],
    {r, d}]];

vfAdd[vfs_List] := Map[fAdd, Transpose[vfs]];
vfNeg[vf_List] := fNeg /@ vf;
vfScaleShift[vf_List, c_, s_Integer] := fScaleShift[#, c, s] & /@ vf;

(* Divide a scalar frame by dA+dB eps.  CASE T is an exact nonnegative
   coefficient recurrence; CASE P is the exact Laurent shift. *)
affineDiv[f_, dA_, dB_, cap_Integer] := Which[
  frameZeroAllQ[f], zFrame[],
  !exactZeroQ[dA], fScalarOp[scalarOp[1/(dA + dB eps), cap - frameMin[f] + 2], f, cap],
  !exactZeroQ[dB], fScaleShift[f, 1/dB, -1],
  True, zFrame[]];

componentFrame[vf_List, r_Integer] := vf[[r]];
replaceComponents[vf_List, cols_List, vals_List] := ReplacePart[vf,
  Thread[cols -> vals]];

blockSolve[rhs_List, dA_, dB_, d0Inv_, q_Integer, cap_Integer] := Module[
  {z = Table[zFrame[], {q}], r},
  z[[q]] = affineDiv[rhs[[q]], dA, dB, cap];
  Do[z[[r]] = affineDiv[fAdd[{rhs[[r]], z[[r + 1]]}], dA, dB, cap],
    {r, q - 1, 1, -1}];
  fScalarOp[d0Inv, #, cap] & /@ z];

(* The exact recurrence shadow.  C is a sparse ragged Laurent jet used to
   model the prospective payload.  E is the exact rational influence
   shadow and is propagated through every complete coupled transition;
   therefore support read from E remains valid even when content initially
   above guardTop crosses the payload boundary after many negative shifts. *)
runColumnShadow[cs_, data_, root_, qpos_Integer, colCursor_Integer,
    nmax_Integer, cap_Integer, checkpoints_List] := Module[
  {d, blocks, P, dD = data["dD"], dN = data["dN"],
   opHi, dOps, nOps, d0Inv, negDepth, maxQ, G, U, R, n, l, j,
   rhs, term, blk, cols, dA, dB, q, solved, assigned, Rt,
   minsByN = <||>, componentDepth, columnDepth = 0, rowMins,
   shadowMismatches = {}, initIdx, zeroVF, addTo, getAssigned, fr,
   payloadMin},
  d = cs["SystemSize"];
  componentDepth = ConstantArray[0, d];
  blocks = DiffExp2`Solve`Private`blockList[cs];
  P = DiffExp2`Solve`Private`logCeiling[cs, root["a"], root["b"], qpos];
  maxQ = Max[blocks[[All, "q"]]];
  negDepth = Max[0, -Min[0, Sequence @@ Select[
      Join[Map[epsVal, Rest[data["dExpr"]]],
        Flatten[Map[epsVal, Rest[data["NhatExpr"]], {3}]]], IntegerQ]]];
  G = negDepth + maxQ + 1;
  opHi = cap + G + 2;
  dOps = scalarOp[#, opHi] & /@ data["dExpr"];
  nOps = matrixOp[#, opHi] & /@ data["NhatExpr"];
  d0Inv = scalarOp[1/data["dExpr"][[1]], opHi];
  zeroVF[] := Table[zFrame[], {d}];
  U = Table[Table[zeroVF[], {P + 2}], {nmax + 1}];
  Do[If[l <= qpos,
      initIdx = colCursor + qpos + 1 - l;
      U[[1, l + 1, initIdx]] = monoFrame[-l, 1]],
    {l, 0, P}];
  addTo[a_, b_] := vfAdd[{a, b}];
  Do[
    R = Table[zeroVF[], {P + 1}];
    Do[
      Do[If[n - j >= 0,
        R[[l + 1]] = addTo[R[[l + 1]],
          fMatrixOp[nOps[[j + 1]], U[[n - j + 1, l + 1]], cap]]],
        {j, 1, Min[n, dN]}];
      Do[If[n - j >= 0,
        term = vfAdd[{
          vfScaleShift[U[[n - j + 1, l + 1]], root["a"] + n - j, 0],
          vfScaleShift[U[[n - j + 1, l + 1]], root["b"], 1],
          vfScaleShift[U[[n - j + 1, l + 2]], 1, 1]}];
        R[[l + 1]] = addTo[R[[l + 1]], vfNeg[
          fScalarOp[dOps[[j + 1]], #, cap] & /@ term]]],
        {j, 1, Min[n, dD]}],
      {l, 0, P}];
    (* Complete TP blocks, top-down in log degree. *)
    Do[
      rhs = addTo[R[[l + 1]], vfNeg[
        fScalarOp[dOps[[1]], #, cap] & /@
          vfScaleShift[U[[n + 1, l + 2]], 1, 1]]];
      Do[
        cols = blk["Cols"]; q = blk["q"];
        dA = simp[root["a"] + n - blk["a"]];
        dB = simp[root["b"] - blk["b"]];
        If[!exactZeroQ[dA] || !exactZeroQ[dB],
          solved = blockSolve[componentFrame[rhs, #] & /@ cols,
            dA, dB, d0Inv, q, cap];
          U[[n + 1, l + 1]] = replaceComponents[
            U[[n + 1, l + 1]], cols, solved]],
        {blk, blocks}],
      {l, P, 0, -1}];
    (* Complete CASE-R assignments before inspecting support. *)
    Do[
      If[exactZeroQ[root["a"] + n - blk["a"]] &&
          exactZeroQ[root["b"] - blk["b"]],
        cols = blk["Cols"]; q = blk["q"];
        Rt = Table[fScalarOp[d0Inv, R[[l + 1, cols[[r]]]], cap],
          {l, 0, P}, {r, q}];
        assigned = <||>;
        Do[If[l + 1 <= P,
          AssociateTo[assigned, {l + 1, q} -> fScaleShift[Rt[[l + 1, q]], 1, -1]]],
          {l, 0, P - 1}];
        getAssigned[key_] := Lookup[assigned, Key[key], zFrame[]];
        Do[Do[If[!KeyExistsQ[assigned, {l, r + 1}],
            AssociateTo[assigned, {l, r + 1} -> fAdd[{
              fScaleShift[getAssigned[{l + 1, r}], 1, 1],
              fNeg[Rt[[l + 1, r]]]}]]],
          {r, 1, q - 1}], {l, P, 0, -1}];
        Do[If[KeyExistsQ[assigned, {l, r}],
            U[[n + 1, l + 1, cols[[r]]]] = assigned[{l, r}]],
          {l, 0, P}, {r, 1, q}]],
      {blk, blocks}];
    rowMins = Table[frameMin[U[[n + 1, l + 1, r]]],
      {l, 0, P}, {r, d}];
    Do[If[IntegerQ[rowMins[[l + 1, r]]],
        componentDepth[[r]] = Max[componentDepth[[r]],
          -rowMins[[l + 1, r]]]], {l, 0, P}, {r, d}];
    columnDepth = Max[columnDepth, Sequence @@ componentDepth];
    If[MemberQ[checkpoints, n], minsByN[n] = <|
      "Depth" -> columnDepth,
      "ComponentDepths" -> componentDepth|>];
    (* Cross-check the sparse payload wherever the exact influence shadow's
       leading term lies inside the retained jet.  E itself is the horizon-
       independent certificate; this check catches implementation drift
       between the prospective ragged payload and its exact shadow. *)
    Do[fr = U[[n + 1, l + 1, r]];
      payloadMin = If[fr["C"] === <||>, Infinity, Min[Keys[fr["C"]]]];
      If[(IntegerQ[rowMins[[l + 1, r]]] &&
            rowMins[[l + 1, r]] <= cap &&
            payloadMin =!= rowMins[[l + 1, r]]) ||
          (rowMins[[l + 1, r]] === Infinity && payloadMin =!= Infinity),
        AppendTo[shadowMismatches,
          {n, l, r, rowMins[[l + 1, r]], payloadMin}]],
      {l, 0, P}, {r, d}],
    {n, 1, nmax}];
  <|"ByN" -> minsByN, "GuardDepth" -> G,
    "ShadowPayloadMismatches" -> shadowMismatches,
    "LogCeiling" -> P|>];

endpointChart[center_] := Module[{sys2, plan, chart},
  sys2 = Join[sys, <|"ExtraSingularFactors" ->
    Select[fixture["ExtraSingularFactors"], !FreeQ[#, var] &]|>];
  plan = catchDE2[DiffExp2`Transport`SegmentLine[sys2, {11/23, center}]];
  If[FailureQ[plan], Return[plan]];
  chart = SelectFirst[Reverse[plan["Charts"]],
    TrueQ[# ["Singular"]] && exactZeroQ[# ["Center"] - center] &, None];
  If[chart === None,
    Return[Failure["EndpointChart", <|"Endpoint" -> center,
      "Detail" -> "SegmentLine did not produce the singular endpoint chart"|>]]];
  chart];

columnSpecs[cs_] := Module[{out = {}, cursor = 0},
  Do[Do[AppendTo[out, <|"Root" -> root, "QPos" -> qpos,
        "Cursor" -> cursor, "Label" -> Length[out] + 1|>],
      {qpos, 0, root["BlockSize"] - 1}]; cursor += root["BlockSize"],
    {fam, cs["Families"]}, {root, fam["Roots"]}];
  out];

emit[a_] := Print["AEF ", ExportString[a, "RawJSON", "Compact" -> True]];

runEndpoint[center_] := Module[
  {chart, cs, data, specs, runs, elapsed, g, scalarW, adaptiveW,
   depths, compDepths, mismatchCount, endpointDepth},
  chart = endpointChart[center];
  If[FailureQ[chart], emit[<|"Endpoint" -> center,
      "Failure" -> ToString[chart]|>]; Return[$Failed]];
  cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart]];
  If[FailureQ[cs], emit[<|"Endpoint" -> center, "Failure" -> ToString[cs]|>];
    Return[$Failed]];
  data = DiffExp2`Solve`Private`clearedSymbolic[cs];
  specs = columnSpecs[cs];
  elapsed = First@AbsoluteTiming[runs = Map[
      runColumnShadow[cs, data, #["Root"], #["QPos"], #["Cursor"],
        maxOrder, guardTop, orders] &, specs]];
  g = Max[runs[[All, "GuardDepth"]]];
  mismatchCount = Total[Length /@ runs[[All, "ShadowPayloadMismatches"]]];
  endpointDepth = Max[Flatten@Table[
    Lookup[Lookup[runs[[All, "ByN"]], n, <||>], "Depth", 0],
    {n, orders}]];
  Do[
    depths = Lookup[Lookup[runs[[All, "ByN"]], n, <||>], "Depth", 0];
    compDepths = Map[Lookup[Lookup[#["ByN"], n, <||>],
        "ComponentDepths", ConstantArray[0, cs["SystemSize"]]] &, runs];
    scalarW = deliveredWidth + 4 n + 5;
    adaptiveW = deliveredWidth + 2 Max[depths] + 5;
    emit[<|"Kind" -> "BananaEndpoint", "Endpoint" -> center,
      "TaylorOrder" -> n, "Columns" -> Length[specs],
      "SystemSize" -> cs["SystemSize"], "PerColumnDepths" -> depths,
      "PerColumnComponentDepths" -> compDepths,
      "ObservedDepth" -> Max[depths], "OneStepGuard" -> g,
      "GuardTop" -> guardTop,
      "ExactTailInfluenceCertified" -> True,
      "ShadowPayloadMismatches" -> mismatchCount,
      "ChartScale" -> N[chart["Scale"], 12],
      "ChartLocalRadius" -> N[chart["LocalRadius"], 12],
      "ScalarWidth" -> scalarW, "AdaptiveWidthProxy" -> adaptiveW,
      "WidthRatio" -> N[scalarW/adaptiveW, 8],
      "QuadraticCostRatio" -> N[(scalarW/adaptiveW)^2, 8],
      "TotalShadowSeconds" -> N[elapsed, 8]|>],
    {n, orders}];
  If[mismatchCount =!= 0,
    emit[<|"Kind" -> "AssertionFailure", "Endpoint" -> center,
      "ShadowPayloadMismatches" -> mismatchCount|>];
    Return[$Failed]];
  <|"CS" -> cs, "Runs" -> runs, "Certified" -> True,
    "MaxObservedDepth" -> endpointDepth|>];

(* Exact controls.  These are intentionally simple enough that the support
   result is a proof, not a numerical comparison. *)
controlRun[name_, A_, seed_, multiplier_, nmax_] := Module[
  {v = AssociationThread[Range[Length[seed]], monoFrame[0, #] & /@ seed],
   mins = {}, next, n, op, depth},
  op = <|"C" -> <|-multiplier -> A|>, "Finite" -> True,
    "Min" -> -multiplier, "E" -> eps^-multiplier A|>;
  Do[next = fMatrixOp[op, Values[v], guardTop];
    v = AssociationThread[Range[Length[seed]], next];
    AppendTo[mins, frameMin /@ Values[v]], {n, nmax}];
  depth = Max[0, Sequence @@ Select[-Flatten[mins], IntegerQ]];
  emit[<|"Kind" -> "Control", "Name" -> name, "TaylorOrder" -> nmax,
    "Depth" -> depth, "FinalMins" ->
      Map[If[IntegerQ[#], #, ToString[#, InputForm]] &, Last[mins]]|>];
  depth];

delayedTailControl[] := Module[{f, steps = guardTop + 2, n, depth},
  f = <|"C" -> exactLaurent[eps^(guardTop + 1)/(1 - eps), guardTop],
    "Finite" -> False, "E" -> eps^(guardTop + 1)/(1 - eps)|>;
  Do[f = fScaleShift[f, 1, -1], {n, steps}];
  depth = Max[0, -frameMin[f]];
  emit[<|"Kind" -> "Control", "Name" -> "DelayedRationalTailCrossesGuard",
    "Transitions" -> steps, "Depth" -> depth,
    "SparsePayloadEmpty" -> (f["C"] === <||>),
    "ExactShadowMin" -> frameMin[f]|>];
  depth];

cScalar = controlRun["ScalarRepeatedPole", {{1}}, {1}, 1, maxOrder];
cIdempotent = controlRun["RankOneIdempotent", {{1, 0}, {0, 0}},
  {1, 0}, 1, maxOrder];
cNilpotent = controlRun["SquareZeroNilpotent", {{0, 1}, {0, 0}},
  {0, 1}, 1, maxOrder];
cRegulator = controlRun["AnalyticRegulatorIsNonzero", {{Global`delta}},
  {1}, 1, maxOrder];
cDelayedTail = delayedTailControl[];
If[!TrueQ[cScalar === maxOrder && cIdempotent === maxOrder &&
    cNilpotent === 1 && cRegulator === maxOrder && cDelayedTail === 1],
  emit[<|"Kind" -> "AssertionFailure", "ControlDepths" ->
    {cScalar, cIdempotent, cNilpotent, cRegulator, cDelayedTail},
    "Expected" -> {maxOrder, maxOrder, 1, maxOrder, 1}|>];
  Quit[1]];

res0 = runEndpoint[0];
res1 = runEndpoint[1];
If[res0 === $Failed || res1 === $Failed, Quit[1]];
If[!TrueQ[res0["Certified"] && res1["Certified"] &&
    res0["MaxObservedDepth"] === 2 && res1["MaxObservedDepth"] === 2],
  emit[<|"Kind" -> "AssertionFailure",
    "EndpointDepths" -> {res0["MaxObservedDepth"],
      res1["MaxObservedDepth"]}, "Expected" -> {2, 2}|>];
  Quit[1]];
Quit[0];
