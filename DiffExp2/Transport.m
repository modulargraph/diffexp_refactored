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
zeroQ[e_] := TrueQ[PossibleZeroQ[Together[e]]];
esQ = DiffExp2`EpsSeries`ESQ;
esNew = DiffExp2`EpsSeries`ESNew; esZero = DiffExp2`EpsSeries`ESZero;
esAdd = DiffExp2`EpsSeries`ESAdd; esScale = DiffExp2`EpsSeries`ESScale;
esTimes = DiffExp2`EpsSeries`ESTimes; esShift = DiffExp2`EpsSeries`ESShift;
esCoeff = DiffExp2`EpsSeries`ESCoefficient; esMin = DiffExp2`EpsSeries`ESMinPower;
esCM = DiffExp2`EpsSeries`ESCompleteMax; esTrim = DiffExp2`EpsSeries`ESTrim;

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
  {caps, zcap, s, xnew, rad},
  (* nearest real on-path singularity ahead *)
  caps = Select[sings, TrueQ[dir*(N[#, 40] - N[xb, 40]) > 0] &];
  zcap = If[caps === {}, dir*Infinity,
    First[SortBy[caps, Abs[N[# - xb, 40]] &]]];
  If[zcap === dir*Infinity,
    (* unbounded: fixed stride of k * |xb|-scale or 1 *)
    s = Max[1, Abs[xb]]/2; xnew = xb + dir*s/k,
    s = k*Abs[N[zcap - xb, 40]]/(1 + k); xnew = xb + dir*s/k];
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
  {eps = DiffExp2`Config`CanonicalEps[], secs, kmin, kmax, out, grouped},
  kmin = ls["EpsWindow", "Min"]; kmax = ls["EpsWindow", "CompleteMax"];
  (* per sector (a,b,p): phase e^(sigma I Pi (a + b eps)) times
     M_{p -> p-j} = (sigma I Pi eps)^j / j!.  The phase's b-part is an
     eps-exponential: expand to the window width. *)
  out = Flatten[Map[Module[{a = #["a"], b = #["b"], p = #["p"], arr = #["Coeffs"],
      phaseA, width = kmax - kmin},
    phaseA = Exp[sigma*I*Pi*a];
    Table[Module[{shiftTot, prefES, newArr},
      (* contribution to target (a, b, p - j):
         phaseA * e^(sigma I Pi b eps) * (sigma I Pi eps)^j / j! *)
      prefES = DiffExp2`EpsSeries`ESFromExpression[
        phaseA*Exp[sigma*I*Pi*b*eps]*(sigma*I*Pi*eps)^j/j!, eps, kmax - kmin + 2];
      newArr = Module[{d = Dimensions[arr]},
        Table[
          Module[{s},
            s = Sum[Module[{kk = k - jj},
              If[kmin <= kk <= kmax && esMin[prefES] <= jj <= esCM[prefES],
                esCoeff[prefES, jj]*arr[[kk - kmin + 1, n + 1, c]], 0]],
              {jj, esMin[prefES], Min[esCM[prefES], k - kmin]}];
            s],
          {k, kmin, kmax}, {n, 0, d[[2]] - 1}, {c, d[[3]]}]];
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
              If[!TrueQ[PossibleZeroQ[db]],
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

MatchWeights[Fmat_List, vIn_List, label_String] := Module[
  {nb = Length[Fmat], FF, vv, perm, mtol, w},
  mtol = DiffExp2`Tolerances`Tol["MatchTol"];
  (* Gaussian elimination over the eps-Laurent field: EpsSeries entries,
     pivots by minimal leading eps-order (then largest leading magnitude),
     ESDivide row operations.  Honest windows propagate; a column with no
     usable pivot = genuinely singular system (E5). *)
  FF = Map[esTrim, Fmat, {2}];
  vv = vIn;
  perm = Range[nb];   (* row order after pivoting *)
  Do[Module[{cands, pivRow, pivLead},
    cands = Select[Range[col, nb], Module[{ld},
      ld = DiffExp2`EpsSeries`ESLeading[FF[[perm[[#]], col]]];
      ld =!= None] &];
    If[cands === {},
      err["E5", <|"Chart" -> label, "Column" -> col,
        "Detail" -> "matching system singular over the eps-Laurent field"|>]];
    pivRow = First[SortBy[cands, Module[{ld =
        DiffExp2`EpsSeries`ESLeading[FF[[perm[[#]], col]]]},
      {ld[[1]], -Abs[N[ld[[2]], 20]]}] &]];
    If[pivRow =!= col, perm[[{col, pivRow}]] = perm[[{pivRow, col}]]];
    Do[Module[{entry = FF[[perm[[r]], col]], factor},
      If[DiffExp2`EpsSeries`ESLeading[entry] =!= None,
        factor = esTrim[DiffExp2`EpsSeries`ESDivide[entry, FF[[perm[[col]], col]]]];
        Do[FF[[perm[[r]], c2]] = esTrim[esAdd[FF[[perm[[r]], c2]],
          esScale[-1, esTimes[factor, FF[[perm[[col]], c2]]]]]],
          {c2, col, nb}];
        vv[[perm[[r]]]] = esTrim[esAdd[vv[[perm[[r]]]],
          esScale[-1, esTimes[factor, vv[[perm[[col]]]]]]]]]],
      {r, col + 1, nb}]],
    {col, nb}];
  (* back substitution *)
  w = Table[None, {nb}];
  Do[Module[{rhs = vv[[perm[[col]]]]},
    Do[rhs = esTrim[esAdd[rhs,
      esScale[-1, esTimes[FF[[perm[[col]], c2]], w[[c2]]]]]],
      {c2, col + 1, nb}];
    w[[col]] = esTrim[DiffExp2`EpsSeries`ESDivide[rhs, FF[[perm[[col]], col]]]]],
    {col, nb, 1, -1}];
  w];

(* ---- 2.10 probe ---- *)

SegmentErrorProbe[ls_Association, tOut_, couplingDepth_Integer] := Module[
  {dec = DiffExp2`Tolerances`EvalErrorSeriesDecrease[Max[couplingDepth, 1]],
   full, red},
  full = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut, "UsePade" -> False];
  red = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut,
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
      basis, probeErrs},
    cs = DiffExp2`Solve`PrepareChart[sys, chart];
    sol = DiffExp2`Solve`SolveChart[cs, req];
    basis = recombineDegenerate[cs, sol["Basis"]["Columns"],
      sol["Basis"]["Specs"]];
    (* match point: the boundary anchor for the FIRST chart (the incoming
       object is only valid at its anchor); thereafter at radius/k of THIS
       chart on the incoming side *)
    matchPt = If[ci === 1, plan["From"],
      chart["Center"] - dir*chart["Radius"]/cfg["DivisionOrder"]];
    (* incoming value at matchPt in the PREVIOUS object's chart coordinate *)
    tIn = matchPt - current["Center"];
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
        "UsePade" -> False];
      vvals = Module[{vv = prevEval["Value"], d2 = cs["SystemSize"]},
        Table[esNew[esMin[vv], Table[esCoeff[vv, k][[c]],
          {k, esMin[vv], esCM[vv]}]], {c, d2}]]];
    (* basis values at the same point, in THIS chart's coordinate *)
    Module[{tLoc = matchPt - chart["Center"], Feval},
      Feval = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#,
        tLoc, "UsePade" -> False]["Value"] &, basis];
      F = Table[esNew[esMin[Feval[[i]]],
        Table[esCoeff[Feval[[i]], k][[c]], {k, esMin[Feval[[i]]], esCM[Feval[[i]]]}]],
        {c, cs["SystemSize"]}, {i, Length[basis]}]];
    w = MatchWeights[F, vvals, chart["Name"]];
    ls = DiffExp2`SectorSeries`CombineLocalSolutions[w, basis];
    (* probe on the APPROACH side (positive chart coordinate after any
       crossing handling; singular charts are one-sided here) *)
    probeErrs = SegmentErrorProbe[ls,
      (matchPt - chart["Center"])/2, sol["CouplingDepth"]];
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
        plan["To"] - current["Center"], "UsePade" -> False]["Value"]]|>];

End[];
EndPackage[];
