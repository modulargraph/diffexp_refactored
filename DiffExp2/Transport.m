(* DiffExp2/Transport.m — segmentation, marching, eps-graded matching,
   crossing.  Spec: Docs/specs/Transport.md (binding); DECISIONS-M0.md.
   V1 notes (recorded): MatchWeights implements the Laurent-graded solve via
   per-column eps-normalization (consumes Solve's honest Laurent columns
   directly); RecombineBasis runs on demand when the eps^0 system is
   singular; the error probe is the full-vs-reduced evaluation at the
   outgoing match point. *)

BeginPackage["DiffExp2`Transport`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`", "DiffExp2`Solve`"}];

FindSingularities::usage = "FindSingularities[sys] gives <|\"All\", \"Real\", \"Factors\"|>: exact deduplicated roots of the matrix singular factors plus per-call ExtraSingularFactors.";
ChartRadius::usage = "ChartRadius[center, allSingularities] gives the exact/numeric complex-plane distance to the nearest OTHER singularity.";
SegmentLine::usage = "SegmentLine[sys, {from, to}] gives the SegmentPlan: charts, radii, match points, digit budget.";
TransportLine::usage = "TransportLine[sys, boundary, plan] runs the marching loop and returns the TransportResult.";
MatchWeights::usage = "MatchWeights[basisValues, incoming, label] solves the eps-graded (Laurent) weight system with loud residual asserts.";
ApplyCrossing::usage = "ApplyCrossing[ls, sigma] applies the crossing operator (phase times unipotent log-chain mixing) so the far side evaluates at positive chart coordinate.";
SegmentErrorProbe::usage = "SegmentErrorProbe[ls, tOut, couplingDepth] gives the full-vs-reduced evaluation error per eps order.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Transport"|>, payload]];
cfg = DiffExp2`Config`CFG;
(* zeroQ: exact-first (see Solve.m).  Together-canonical rational functions
   over Q(i) in non-numeric symbols: === 0 decides; inexact/radical/
   numeric-constant forms (RootReduce'd radii etc.) keep the PossibleZeroQ
   fallback unchanged. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && TrueQ[PossibleZeroQ[c]]);
zeroQ[e_] := zeroCanQ[Together[e]];
esQ = DiffExp2`EpsSeries`ESQ;
esNew = DiffExp2`EpsSeries`ESNew; esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd; esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes; esShift = DiffExp2`EpsSeries`ESShift;
esCoeff = DiffExp2`EpsSeries`ESCoefficient; esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax;

(* ---- 2.1 singularities ---- *)

FindSingularities[sys_Association] := Module[
  {var = sys["Variable"], facs, extra, all, roots},
  facs = Lookup[sys, "SingularFactors", {}];
  extra = Lookup[sys, "ExtraSingularFactors", {}];
  facs = DeleteDuplicates[Join[facs, extra]];
  roots = Map[# -> DeleteDuplicates[var /. Solve[# == 0, var]] &, facs];
  all = DeleteDuplicates[Flatten[Last /@ roots],
    TrueQ[PossibleZeroQ[RootReduce[#1 - #2]]] &];
  <|"All" -> all,
    "Real" -> Select[all, zeroQ[Im[RootReduce[#]]] &],
    "Factors" -> Association[roots]|>];

ChartRadius[center_, all_List] := Module[
  {others = Select[all, !TrueQ[PossibleZeroQ[RootReduce[# - center]]] &], diffs},
  If[others === {}, Infinity,
    (* EXACT distances when inputs are exact: inexact radii leak into match
       points and destroy precision through symbolic-log numericization *)
    diffs = Map[If[FreeQ[{#, center}, _?InexactNumberQ],
      RootReduce[Abs[# - center]], Abs[N[# - center, 40]]] &, others];
    Min[diffs]]];

(* ---- 2.4 GetCPL/GetCPR geometry ---- *)

nextCenter[xb_, sings_, dir_, k_, allSings_] := Module[
  {caps, zcap, s, xnew, rad, rb},
  (* nearest real on-path singularity ahead *)
  caps = Select[sings, TrueQ[dir*(N[#, 40] - N[xb, 40]) > 0] &];
  zcap = If[caps === {}, dir*Infinity,
    First[SortBy[caps, Abs[N[# - xb, 40]] &]]];
  If[zcap === dir*Infinity,
    (* unbounded: fixed stride of k * |xb|-scale or 1 *)
    s = Max[1, Abs[xb]]/2; xnew = xb + dir*s/k,
    s = k*Abs[N[zcap - xb, 40]]/(1 + k); xnew = xb + dir*s/k];
  (* CHAIN INVARIANT: the next chart's match point (at r_new/k from its
     center) must lie INSIDE this chart's disk.  The step geometry above
     only sees singularities AHEAD - a singularity just BEHIND the path
     (banana: x = 1/2 next to the anchor) caps r_b and breaks coverage.
     step <= r_b/2 guarantees |step| + r_new/k <= (1/2 + (3/2)/k) r_b < r_b
     for k >= 4 (radius is 1-Lipschitz in the center). *)
  rb = ChartRadius[xb, allSings];
  If[rb =!= Infinity && TrueQ[s/k > rb/2],
    s = k*rb/2; xnew = xb + dir*s/k];
  (* complex-distance cap: re-solve placement against the capped radius *)
  rad = ChartRadius[xnew, allSings];
  If[rad < s,
    s = k*rad/(1 + k); xnew = xb + dir*s/k];
  (* exact-ify the center COARSELY: the center is a placement CHOICE; a
     simple nearby rational keeps Indicial's exact algebra on small
     fractions (20-digit-denominator centers ground PossibleZeroQ into
     meprec storms).  Geometry self-corrects: the true radius and match
     point are recomputed from the actual center downstream. *)
  xnew = Module[{cand = Rationalize[N[xnew, 20], Abs[N[s, 20]]/(8 k)]},
    (* never land ON a singularity *)
    If[AnyTrue[allSings, TrueQ[PossibleZeroQ[RootReduce[# - cand]]] &],
      Rationalize[N[xnew, 20], Abs[N[s, 20]]/(64 k)], cand]];
  xnew];

(* ---- 2.5 digit budget ---- *)

DigitBudget[ag_, nseg_Integer] := Module[{wp = cfg["WorkingPrecision"], dn, cd},
  cd = cfg["ChopPrecision"];
  dn = If[ag === "?", cd,
    ag + Ceiling[Log10[Max[nseg, 1]]] + DiffExp2`Tolerances`$SafetyDigits];
  If[dn + DiffExp2`Tolerances`ChopReserve[wp, cd] > wp,
    err["E3", <|"DigitsNeeded" -> dn, "ChopReserve" -> wp - cd,
      "WorkingPrecision" -> wp,
      "Detail" -> "digit budget exceeds WorkingPrecision; raise WP or lower AccuracyGoal"|>]];
  dn];

(* ---- 2.3 segmentation ---- *)

SegmentLine[sys_Association, {from_, to_}] := Module[
  {sings, real, dir, k = cfg["DivisionOrder"], charts = {}, cur, interior,
   endpointSingular, all, guard = 0},
  sings = FindSingularities[sys];
  all = sings["All"]; real = sings["Real"];
  dir = Sign[to - from];
  If[dir === 0, err["E1", <|"From" -> from, "To" -> to, "Detail" -> "empty line"|>]];
  endpointSingular = AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - to]]] &];
  If[AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - from]]] &],
    err["E1", <|"From" -> from,
      "Detail" -> "transport FROM a singular point requires a singular boundary object (not supported in v1 marching start)"|>]];
  interior = Sort[Select[real, TrueQ[dir*(N[#, 40] - N[from, 40]) > 0] &&
      TrueQ[dir*(N[to, 40] - N[#, 40]) > 0] &], dir*N[#1 - #2, 40] < 0 &];
  cur = from;
  Module[{targets = Join[interior, {to}], reached = False},
    Do[Module[{target = targets[[ti]], targetSingular},
      targetSingular = ti < Length[targets] || endpointSingular;
      While[True,
        guard++;
        If[guard > 500, err["E1", <|"Detail" -> "segmentation runaway"|>]];
        Module[{radCur, nxt, radTarget, matchTarget},
          radTarget = Min[ChartRadius[target, all], 2*Abs[to - from]];
          (* the point the LAST chart must reach: the target chart's
             incoming match point (radius/k before the target) *)
          matchTarget = If[targetSingular, target - dir*radTarget/k, target];
          nxt = nextCenter[cur, real, dir, k, all];
          radCur = Min[ChartRadius[nxt, all], 2*Abs[to - from]];
          If[TrueQ[Abs[N[matchTarget - nxt, 40]] <= radCur/k] ||
             TrueQ[dir*(N[nxt, 40] - N[matchTarget, 40]) >= 0],
            (* nxt covers the target's match point: append it as the last
               regular chart; a singular target additionally gets its own
               chart (matched from nxt) *)
            AppendTo[charts, <|"Center" -> nxt, "Singular" -> False|>];
            If[targetSingular,
              AppendTo[charts, <|"Center" -> target, "Singular" -> True|>]];
            Break[],
            AppendTo[charts, <|"Center" -> nxt, "Singular" -> False|>];
            cur = nxt]]];
      cur = target],
      {ti, Length[Join[interior, {to}]]}]];
  (* attach radii, match points, names; radii capped at line scale
     (a validity bound: capping is conservative; uncapped Infinity poisons
     the match-point arithmetic on singularity-free systems) *)
  charts = MapIndexed[Module[{c = #1, rad},
    rad = Min[ChartRadius[c["Center"], all], 2*Abs[to - from]];
    Join[c, <|"Radius" -> rad, "ChartVar" -> Global`t,
      "Name" -> "seg" <> ToString[First[#2]] <> "@" <>
        ToString[N[c["Center"], 6]],
      "Prescriptions" -> {}|>]] &, charts];
  <|"From" -> from, "To" -> to, "Direction" -> dir,
    "Charts" -> charts, "EndpointIsSingular" -> endpointSingular,
    "DigitsNeeded" -> DigitBudget[cfg["AccuracyGoal"], Length[charts]],
    "Singularities" -> sings|>];

(* ---- 2.9 crossing ---- *)

ApplyCrossing[ls_Association, sigma_] := Module[
  {eps = DiffExp2`Config`CanonicalEps[], kmin, kmax, nk, out},
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  nk = kmax - kmin + 1;
  (* per sector (a,b,p): phase e^(sigma I Pi (a + b eps)) times
     M_{p -> p-j} = (sigma I Pi eps)^j / j!.  The phase's b-part is an
     eps-exponential: expand to the window width.  Slab-vectorized (the
     CombineLocalSolutions pattern): the per-(k,n,c) Sum becomes weighted
     SHIFTED-SLAB adds over the k-axis — one whole-array op per prefactor
     eps-order.  The crossing phase is e^(sigma I Pi (a + n + b eps)) PER
     t-COLUMN: the integer part depends on n -> the (-1)^n per-column
     factor, pre-applied to the source slabs (exact, jj-independent). *)
  out = Flatten[Map[Module[{a = #["a"], b = #["b"], p = #["p"], arr = #["Coeffs"],
      phaseA, dims, arrS},
    phaseA = Exp[sigma*I*Pi*a];
    dims = Dimensions[arr];
    arrS = Module[{signs = Table[(-1)^n, {n, 0, dims[[2]] - 1}]},
      Map[signs*# &, arr]];
    Table[Module[{prefES, newArr},
      (* contribution to target (a, b, p - j):
         phaseA * e^(sigma I Pi b eps) * (sigma I Pi eps)^j / j! *)
      prefES = DiffExp2`EpsSeries`ESFromExpression[
        phaseA*Exp[sigma*I*Pi*b*eps]*(sigma*I*Pi*eps)^j/j!, eps, kmax - kmin + 2];
      newArr = ConstantArray[0, {nk, dims[[2]], dims[[3]]}];
      (* newArr[k] += pref[jj]*arrS[k - jj] <=> target slabs [1+jj, nk]
         from source slabs [1, nk - jj] (the same index set as the old
         per-element Sum, in the same ascending-jj addition order) *)
      Do[Module[{pc = esCoeff[prefES, jj], lo, hi},
        If[pc =!= 0,
          lo = Max[1, 1 - jj]; hi = Min[nk, nk - jj];
          If[lo <= hi,
            newArr[[lo + jj ;; hi + jj]] += pc*arrS[[lo ;; hi]]]]],
        {jj, esMin[prefES], Min[esCM[prefES], kmax - kmin]}];
      <|"a" -> a, "b" -> b, "p" -> p - j, "Coeffs" -> newArr|>],
      {j, 0, p}]] &, ls["Sectors"]], 1];
  DiffExp2`SectorSeries`CanonicalizeLocalSolution[
    Join[ls, <|"Sectors" -> out|>]]];

(* ---- 2.8 RecombineBasis: remove eps=0 degeneracy of the fundamental
   matrix (the log x = (x^(2eps)-1)/(2eps) class).  Tag-driven off
   Indicial`EpsDegenerateFamilies (DEC-7): per degenerate family, ordered
   by b: B_1 = S_1, B_m = (S_m - S_1)/((b_m - b_1) eps); recursively on
   {B_2, ...} up to the family's recorded degeneracy depth. *)

recombineDegenerate[cs_, basis_List, specs_List] := Module[
  {degs, newBasis = basis, width},
  degs = DiffExp2`Indicial`EpsDegenerateFamilies[cs["IndicialData"]];
  If[Environment["DEBUG_RECOMBINE"] === "1",
    Print["RECOMBINE chart ", cs["Center"], " degs=", degs,
      " specs=", {#["a"], #["b"], #["Family"], #["ChainPos"]} & /@ specs]];
  If[degs === {}, Return[newBasis]];
  width = 4 + Max[0, Max[Map[#["EpsWindow", "CompleteMax"] -
    #["EpsWindow", "Min"] &, basis]]];
  Do[Module[{fi = rec["FamilyIndex"], r0 = rec["EpsZeroDegeneracy"], cols, bs},
    cols = Select[Range[Length[specs]], specs[[#]]["Family"] === fi &];
    (* value-level degeneracy pairs columns with the SAME exact a and
       DISTINCT b (different-a columns are different functions); per
       same-a group, the recursive (S_m - S_1)/((b_m - b_1) eps) ladder *)
    Module[{groups},
      groups = GatherBy[Select[cols, specs[[#]]["ChainPos"] === 0 &],
        Together[specs[[#]]["a"]] &];
      Do[Module[{bs = SortBy[grp, N[specs[[#]]["b"], 20] &], active, pass = 0},
        active = bs;
        While[Length[active] >= 2 && pass < Max[r0, 1],
          Module[{base = First[active], rest = Rest[active]},
            Do[Module[{m = rest[[ri]], db, wPlus, wMinus},
              db = Together[specs[[m]]["b"] - specs[[base]]["b"]];
              (* exact tag difference, already Together'd *)
              If[!zeroCanQ[db],
                wPlus = DiffExp2`EpsSeries`ESNew[-1, PadRight[{1/db}, width]];
                wMinus = DiffExp2`EpsSeries`ESNew[-1, PadRight[{-1/db}, width]];
                newBasis[[m]] = DiffExp2`SectorSeries`CombineLocalSolutions[
                  {wPlus, wMinus}, {newBasis[[m]], newBasis[[base]]}]]],
              {ri, Length[rest]}];
            active = rest; pass++]]],
        {grp, Select[groups, Length[#] >= 2 &]}]]],
    {rec, degs}];
  newBasis];

(* ---- 2.7 matching ---- *)

(* endpoint-matching log-branch convention: prescriptions when derivable,
   else the fixed +1 convention (weights absorb it; Euclidean FT results
   stay real - verified against the oracle) *)
sigmaFor[ls_] := Module[{s = DiffExp2`SectorSeries`ChartImSign[ls]},
  If[MemberQ[{1, -1}, s], s, 1]];

(* ===================== MATCHWEIGHTS FRAME KERNELS =====================
   Gaussian elimination over the eps-Laurent field on PLAIN WINDOWED
   COEFFICIENT LISTS {min, coeffs} (index i <-> eps^(min+i-1); the honest
   complete window is [min, min + Length[coeffs] - 1]) — the runRecursion
   packed-list style applied to the matching solve.  This bypasses the
   per-op EpsSeries object layer (Association churn, per-coefficient
   window reads) that made the nb^2..nb^3 row operations the per-chart
   hot spot at d = 7.  EpsSeries objects appear only at the boundaries
   (matrix/rhs in, weights out).
   WINDOW HONESTY: every kernel implements EXACTLY the EpsSeries window
   rule it replaces (ESAdd min-rules; the ESTimes Cauchy window; the
   ESInvert/ESDivide structure shift [-L, CM - 2L] for division by a
   pivot leading at relative index L; ESTrim's value-classified leading
   advance with CompleteMax untouched), so the delivered weight windows
   equal the former per-object path's windows order for order. *)

(* seriesScale mirror: max |numeric coefficient|, symbolic excluded *)
mwScale[c_List] := Max[0, Sequence @@ (If[NumericQ[#], Abs[N[#, 10]], Nothing] & /@ c)];

(* EpsSeries norm mirror: numerics pass through; symbolics Together *)
mwNorm[c_List] := If[FreeQ[c, _Symbol], c,
  Map[If[NumericQ[#], #, Together[Expand[#]]] &, c]];

(* ESTrim mirror: advance min past negligible leading coefficients under
   the SAME relative gate (ESCoeffZeroQ -> Tol["LaurentLeadTol"], scale
   over the stored window, the ambiguity band stays loud); the top never
   changes; all-negligible gives the canonical zero {top, {0}}.  This is
   the trim-before-normalize guard: numerically-tiny leading storage
   never becomes a pivot. *)
mwTrim[{min_, c_}] := Module[{scale = mwScale[c], n = Length[c], i = 1},
  While[i <= n && TrueQ[DiffExp2`EpsSeries`ESCoeffZeroQ[c[[i]], scale]], i++];
  If[i > n, {min + n - 1, {0}}, {min + i - 1, Drop[c, i - 1]}]];

(* trimmed-entry zero test (the ESLeading === None equivalent) *)
mwZeroQ[{_, c_}] := c === {0};

mwNeg[{m_, c_}] := {m, -c};

(* ESAdd mirror: min = Min[mins], top = Min[tops] — never zero-padded *)
mwAdd[{am_, ac_}, {bm_, bc_}] := Module[
  {min = Min[am, bm], top = Min[am + Length[ac], bm + Length[bc]] - 1, out},
  out = ConstantArray[0, top - min + 1];
  Module[{lo = am - min + 1, hi = Min[am + Length[ac] - 1, top] - min + 1},
    If[lo <= hi, out[[lo ;; hi]] += Take[ac, hi - lo + 1]]];
  Module[{lo = bm - min + 1, hi = Min[bm + Length[bc] - 1, top] - min + 1},
    If[lo <= hi, out[[lo ;; hi]] += Take[bc, hi - lo + 1]]];
  {min, mwNorm[out]}];

(* ESTimes mirror via ListConvolve: mins add; honest length =
   Min[operand lengths] (the Cauchy-product window rule) *)
mwMul[{am_, ac_}, {bm_, bc_}] :=
  {am + bm, mwNorm[Take[ListConvolve[ac, bc, {1, -1}, 0],
    Min[Length[ac], Length[bc]]]]};

(* ESDivide mirror as the direct quotient recursion (b q = a through the
   honest order — coefficient-identical to ESTimes[a, ESInvert[b]]):
   min = aMin - bMin, length = Min[operand lengths].  Division by a pivot
   leading at index L > 0 shifts content down by L AND costs completeness
   at the top, exactly ESInvert's [-L, CM - 2L] window rule.  The caller
   guarantees a trimmed denominator with a non-negligible (hence nonzero)
   leading entry. *)
mwDiv[{am_, ac_}, {bm_, bc_}] := Module[
  {len = Min[Length[ac], Length[bc]], b0 = First[bc], sym, q},
  sym = !FreeQ[{ac, bc}, _Symbol];
  q = ConstantArray[0, len];
  q[[1]] = If[sym, Together[ac[[1]]/b0], ac[[1]]/b0];
  Do[q[[m + 1]] = Module[
      {v = (ac[[m + 1]] - Sum[bc[[j + 1]]*q[[m - j + 1]],
          {j, 1, Min[m, Length[bc] - 1]}])/b0},
      If[sym, Together[v], v]],
    {m, 1, len - 1}];
  {am - bm, q}];

mwFromES[s_] := {esMin[s], s["Coeffs"]};
mwToES[{m_, c_}] := esNew[m, c];

MatchWeights[Fmat_List, vIn_List, label_String] := Module[
  {nb = Length[Fmat], FF, vv, perm, w},
  (* Gaussian elimination over the eps-Laurent field on frame lists:
     pivots by minimal leading eps-order (then largest leading magnitude),
     quotient-recursion row operations.  Honest windows propagate through
     the mirrored EpsSeries rules; a column with no usable pivot =
     genuinely singular system (E5). *)
  FF = Map[mwTrim[mwFromES[#]] &, Fmat, {2}];
  vv = mwFromES /@ vIn;
  perm = Range[nb];   (* row order after pivoting *)
  Do[Module[{cands, pivRow},
    cands = Select[Range[col, nb], !mwZeroQ[FF[[perm[[#]], col]]] &];
    If[cands === {},
      err["E5", <|"Chart" -> label, "Column" -> col,
        "Detail" -> "matching system singular over the eps-Laurent field"|>]];
    pivRow = First[SortBy[cands, Module[{e = FF[[perm[[#]], col]]},
      {e[[1]], -Abs[N[e[[2, 1]], 20]]}] &]];
    If[pivRow =!= col, perm[[{col, pivRow}]] = perm[[{pivRow, col}]]];
    Do[Module[{entry = FF[[perm[[r]], col]], factor},
      If[!mwZeroQ[entry],
        factor = mwTrim[mwDiv[entry, FF[[perm[[col]], col]]]];
        Do[FF[[perm[[r]], c2]] = mwTrim[mwAdd[FF[[perm[[r]], c2]],
          mwNeg[mwMul[factor, FF[[perm[[col]], c2]]]]]],
          {c2, col, nb}];
        vv[[perm[[r]]]] = mwTrim[mwAdd[vv[[perm[[r]]]],
          mwNeg[mwMul[factor, vv[[perm[[col]]]]]]]]]],
      {r, col + 1, nb}]],
    {col, nb}];
  (* back substitution *)
  w = Table[None, {nb}];
  Do[Module[{rhs = vv[[perm[[col]]]]},
    Do[rhs = mwTrim[mwAdd[rhs,
      mwNeg[mwMul[FF[[perm[[col]], c2]], w[[c2]]]]]],
      {c2, col + 1, nb}];
    w[[col]] = mwTrim[mwDiv[rhs, FF[[perm[[col]], col]]]]],
    {col, nb, 1, -1}];
  mwToES /@ w];

(* ---- 2.10 probe ---- *)

SegmentErrorProbe[ls_Association, tOut_, couplingDepth_Integer] := Module[
  {dec = DiffExp2`Tolerances`EvalErrorSeriesDecrease[Max[couplingDepth, 1]],
   full, red},
  full = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut, "UsePade" -> False,
    "ImSign" -> sigmaFor[ls]];
  red = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut, "ImSign" -> sigmaFor[ls],
    "UsePade" -> False, "TOrderReduction" -> dec];
  Table[Module[{kf = esCoeff[full["Value"], k],
      kr = If[esMin[red["Value"]] <= k <= esCM[red["Value"]],
        esCoeff[red["Value"], k], 0*esCoeff[full["Value"], k]]},
    Max[0, Sequence @@ (Abs[N[#, 20]] & /@ Select[Flatten[{kf - kr}], NumericQ])]],
    {k, esMin[full["Value"]], esCM[full["Value"]]}]];

(* ---- 2.6 marching ---- *)

TransportLine[sys_Association, boundary_, plan_Association] := Module[
  {charts = plan["Charts"], dir = plan["Direction"], current, errAcc,
   req, expOrd = cfg["ExpansionOrder"], epsOrd = cfg["EpsilonOrder"],
   lastSingular = False, lastLS = None, lastChart = None, kept = {}},
  req = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> epsOrd|>,
    "TOrder" -> expOrd|>;
  (* boundary: LocalSolution (anchored at plan From) or plain values matrix *)
  current = If[AssociationQ[boundary] && KeyExistsQ[boundary, "Sectors"], boundary,
    Module[{vals = boundary, d, kmax},
      d = Length[vals]; kmax = Length[First[vals]] - 1;
      <|"Center" -> plan["From"], "ChartMap" -> <|"Center" -> plan["From"], "Scale" -> 1|>,
        "Radius" -> Min[ChartRadius[plan["From"], plan["Singularities"]["All"]],
          2*Abs[plan["To"] - plan["From"]]],
        "Sectors" -> {<|"a" -> 0, "b" -> 0, "p" -> 0,
          "Coeffs" -> Table[Table[Table[vals[[c, k + 1]], {c, d}], {n, 0, 0}],
            {k, 0, kmax}]|>},
        "EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> kmax|>,
        "TWindow" -> <|"CompleteMax" -> 0|>,
        "ErrorEstimate" -> ConstantArray[0, kmax + 1],
        "Prescriptions" -> {}|>]];
  errAcc = None;
  Do[Module[{chart = charts[[ci]], cs, sol, matchPt, tIn, vvals, F, w, ls,
      basis, probeErrs, valueMode, couplingDepth = 0},
    If[Environment["DEBUG_CHART"] === "1",
      Print["CHART ", chart["Name"], " prep start t=", SessionTime[]]];
    cs = DiffExp2`Solve`PrepareChart[sys, chart];
    (* PROTOTYPE (env-gated, default off): value-vector propagation for
       REGULAR interior charts (Docs/PerfGapAnalysis.md lever 1).  The
       incoming object is evaluated AT THIS CHART'S CENTER — inside the
       previous chart's disk by the predivision geometry (step/radius <=
       1/(1+k)); EvaluateLocalSolution's radius assert is the loud
       backstop — and that value is the t^0 Cauchy datum of ONE
       d-dimensional recursion (Solve`SolveValueRegular), replacing the
       d-column basis + MatchWeights + CombineLocalSolutions.  Singular
       charts and the first chart (anchor-only incoming data) keep the
       basis+matching path unchanged. *)
    valueMode = Environment["DE2_VALUE_TRANSPORT"] === "1" && ci > 1 &&
      !TrueQ[chart["Singular"]] &&
      TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
      (* conservative geometry pre-check: the center must sit WELL inside
         the previous object's disk (margin 9/10); otherwise fall back to
         the basis path (a performance choice, not an ambiguity) *)
      TrueQ[Abs[N[chart["Center"] - current["Center"], 30]] <
        9/10*N[current["Radius"], 30]];
    (* match point: the boundary anchor for the FIRST chart (the incoming
       object is only valid at its anchor); thereafter at radius/k of THIS
       chart on the incoming side *)
    matchPt = If[ci === 1, plan["From"],
      Module[{raw = chart["Center"] - dir*chart["Radius"]/cfg["DivisionOrder"]},
        (* simple rational match point: algebraic radii otherwise force
           exact algebraic arithmetic through every evaluation *)
        Rationalize[N[raw, 20],
          N[chart["Radius"], 20]/(100*cfg["DivisionOrder"])]]];
    (* incoming value in the PREVIOUS object's chart coordinate: at the
       match point (basis path) or at this chart's center (value mode) *)
    tIn = If[valueMode, chart["Center"], matchPt] - current["Center"];
    Module[{prevEval, sigma},
      (* crossing: if the previous chart was singular and matchPt lies on its
         far side (sign of tIn relative to approach), apply the operator *)
      If[lastSingular && TrueQ[N[tIn, 30] dir > 0],
        sigma = DiffExp2`SectorSeries`ChartImSign[current];
        If[!MemberQ[{1, -1}, sigma],
          If[AnyTrue[current["Sectors"],
              !IntegerQ[#["a"]] || !zeroQ[#["b"]] || #["p"] > 0 &],
            err["E8", <|"Chart" -> chart["Name"],
              "Detail" -> "crossing a multivalued singular chart without a derivable Im-sign (missing DeltaPrescriptions)"|>],
            sigma = 1]];
        current = ApplyCrossing[current, sigma];
        tIn = -tIn  (* far side evaluates at positive u *)];
      prevEval = DiffExp2`SectorSeries`EvaluateLocalSolution[current, tIn,
        "UsePade" -> False, "ImSign" -> sigmaFor[current]];
      vvals = Module[{vv = prevEval["Value"], d2 = cs["SystemSize"]},
        Table[esNew[esMin[vv], Table[esCoeff[vv, k][[c]],
          {k, esMin[vv], esCM[vv]}]], {c, d2}]]];
    If[valueMode,
      If[Environment["DEBUG_CHART"] === "1",
        Print["CHART value-solve start t=", SessionTime[]]];
      ls = DiffExp2`Solve`SolveValueRegular[cs, req, vvals];
      If[Environment["DEBUG_CHART"] === "1",
        Print["CHART value-solve done t=", SessionTime[]]];
      ,
      If[Environment["DEBUG_CHART"] === "1",
        Print["CHART solve start t=", SessionTime[]]];
      sol = DiffExp2`Solve`SolveChart[cs, req];
      If[Environment["DEBUG_CHART"] === "1",
        Print["CHART solve done t=", SessionTime[]]];
      basis = recombineDegenerate[cs, sol["Basis"]["Columns"],
        sol["Basis"]["Specs"]];
      couplingDepth = sol["CouplingDepth"];
      (* basis values at the same point, in THIS chart's coordinate *)
      Module[{tLoc = matchPt - chart["Center"], Feval},
        Feval = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#,
          tLoc, "UsePade" -> False, "ImSign" -> sigmaFor[#]]["Value"] &, basis];
        F = Table[esNew[esMin[Feval[[i]]],
          Table[esCoeff[Feval[[i]], k][[c]], {k, esMin[Feval[[i]]], esCM[Feval[[i]]]}]],
          {c, cs["SystemSize"]}, {i, Length[basis]}]];
      w = MatchWeights[F, vvals, chart["Name"]];
      ls = DiffExp2`SectorSeries`CombineLocalSolutions[w, basis]];
    (* probe on the APPROACH side (positive chart coordinate after any
       crossing handling; singular charts are one-sided here) *)
    probeErrs = SegmentErrorProbe[ls,
      (matchPt - chart["Center"])/2, couplingDepth];
    errAcc = If[errAcc === None, probeErrs,
      Module[{l1 = Length[errAcc], l2 = Length[probeErrs]},
        PadRight[errAcc, Max[l1, l2]] + PadRight[probeErrs, Max[l1, l2]]]];
    If[Max[errAcc] > 1,
      err["E10", <|"Chart" -> chart["Name"], "Errors" -> N[errAcc, 4],
        "Detail" -> "accumulated error estimate exceeds 1; transport aborted"|>]];
    current = ls;
    AppendTo[kept, <|"Chart" -> chart, "LocalSolution" -> ls|>];
    lastSingular = TrueQ[chart["Singular"]];
    lastChart = chart],
    {ci, Length[charts]}];
  <|"Final" -> current,
    "Charts" -> kept,
    "EndpointIsSingular" -> plan["EndpointIsSingular"],
    "ErrorEstimate" -> errAcc,
    "Value" -> If[plan["EndpointIsSingular"], None,
      DiffExp2`SectorSeries`EvaluateLocalSolution[current,
        plan["To"] - current["Center"], "UsePade" -> False,
        "ImSign" -> sigmaFor[current]]["Value"]]|>];

End[];
EndPackage[];
