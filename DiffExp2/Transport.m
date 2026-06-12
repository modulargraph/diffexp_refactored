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
ValidatePlan::usage = "ValidatePlan[plan] statically audits the chart chain: every incoming match point (the shared chartMatchPoint formula) must lie strictly inside the producing chart's disk, singular handoffs must be geometrically possible. Loud E8 on violation; returns the plan.";
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

(* stepDivisor: the chart PLACEMENT stride divisor k_eff (<= DivisionOrder
   k when Automatic).  Two distinct roles were coupled in one k:
     (a) the MATCH-POINT ratio: match points sit at radius/k — truncation
         tail there is (1/k)^(ExpansionOrder+1) ((1/4)^41 ~ 1e-25 at the FT
         defaults k = 4, ExpansionOrder = 40: ~11 decades of headroom below
         the 1e-14 design line);
     (b) the STRIDE: marching toward a dominating singularity, radius s and
         step g = s/k_eff satisfy s = (distance) - g, so the remaining
         distance shrinks by r = k_eff/(1+k_eff) per chart and the chart
         count grows like log(d0/dstop)/log((1+k_eff)/k_eff) — k = 4 gives
         r = 4/5 and ~10.3 charts per decade of approach; LARGER divisors
         mean MORE charts.
   Decoupling: placement uses k_eff <= k (bigger steps), match points keep
   radius/k.  Evaluation-ratio inventory under the nextCenter invariants
   (G1)+(G2) below, each ratio in units of the OWNING chart's radius:
     - incoming match point in the new chart:            1/k
     - incoming match point in the PREVIOUS chart:      <= ~0.41
       ((G2) clamps h = g - R_new/k <= (9/10) clamp prevRad directly,
       clamp = Max[1/(2 k_eff), 1/k]; the h < 0 side is <= prevRad/(k-1))
     - LineIntegral tile edges (midpoints between centers, API.m):
       backward edge g/(2 R_new) <= 1/(2 k_eff) when (G1) holds
       (R_new >= k_eff g); in the back-binding case R_new >= g keeps it
       <= 1/2 — the SAME class the old geometry produced next to singular
       charts; forward edge g/(2 prevRad) <= ~0.41 by (G2) + the
       1-Lipschitz bound R_new <= prevRad + g
   so the controlled worst ratio is ~Max[1/(2 k_eff), 1/k, 0.41-class] and
   the per-evaluation tail is ratio^(ExpansionOrder+1).  Automatic picks
   the largest stride whose 1/(2 k_eff) ratio keeps that tail below 1e-14,
     k_eff = Min[k, Max[5/4, ceil16(10^(14/(eo+1))/2)]]
   (ceil16 = round UP to sixteenths; the 5/4 floor is where the 0.41-class
   ratios stop improving).  At eo = 40 this gives 5/4: controlled ratios
   <= 0.41, tails <= 0.41^41 ~ 1.3e-16, and the approach ratio drops from
   4/5 to 5/9 — ln(9/5)/ln(5/4) ~ 2.6x fewer marching charts per decade.
   At eo <= ~14 the formula exceeds k and the Min restores the classic
   coupled geometry exactly.  TransportLine's SegmentErrorProbe probes at
   the design ratio Max[1/(2 k_eff), 1/k], so the accumulated error
   accounting sees these wider handoffs honestly. *)
stepDivisor[k_] := Module[{raw = cfg["StepDivisionOrder"], eo},
  eo = cfg["ExpansionOrder"];
  If[raw === Automatic,
    Min[k, Max[5/4, Ceiling[N[8*10^(14/(eo + 1)), 20]]/16]],
    raw]];

(* nextCenter: place the next chart center marching dir-ward from xb.
   k = match-point divisor (DivisionOrder), keff = stepDivisor stride
   divisor, prevRad = the previous chart's (line-capped) radius,
   matchTarget = the current target's incoming match point, lineCap = the
   2|to - from| radius cap used throughout SegmentLine. *)
nextCenter[xb_, sings_, dir_, k_, keff_, allSings_, prevRad_, matchTarget_,
    lineCap_] := Module[
  {caps, zcap, s, g, xnew, rad, radC, clamp, h, tgtGap, it},
  (* nearest real on-path singularity ahead *)
  caps = Select[sings, TrueQ[dir*(N[#, 40] - N[xb, 40]) > 0] &];
  zcap = If[caps === {}, dir*Infinity,
    First[SortBy[caps, Abs[N[# - xb, 40]] &]]];
  If[zcap === dir*Infinity,
    (* unbounded: fixed stride scale of |xb| or 1 *)
    s = Max[1, Abs[xb]]/2,
    (* dominated march: radius s touches the cap and xb sits at s/keff of
       the new chart; remaining distance shrinks by keff/(1+keff) *)
    s = keff*Abs[N[zcap - xb, 40]]/(1 + keff)];
  g = s/keff;
  (* never step past the target's incoming match point: the cover test and
     the endpoint handoff assume |matchTarget - center| <= radius-scale *)
  tgtGap = dir*(N[matchTarget, 40] - N[xb, 40]);
  If[TrueQ[tgtGap > 0] && TrueQ[g > tgtGap], g = tgtGap; s = keff*g];
  xnew = xb + dir*g;
  (* (G1) complex-distance cap, ITERATED toward its fixed point rad >= s
     (i.e. R_new >= keff*g; the old one-shot re-solve could leave
     rad << s).  The contraction g <- rad/(1+keff) converges only when the
     binding singularity lies AHEAD (rad grows as g shrinks).  When it lies
     BEHIND xb (receding from a just-crossed singular chart) rad ~ u_b + g
     SHRINKS with g and the iteration would contract g to 0 — there
     R_new >= keff*g is unsatisfiable for keff > 1, so STOP at the first
     non-improving contraction and accept the back-binding placement:
     R_new >= g keeps the backward tile edge <= 1/2 (the pre-existing class
     next to singular charts) and (G2) below enforces the honest handoff
     bound, which also caps the forward edge in the 0.41 class. *)
  rad = ChartRadius[xnew, allSings];
  it = 0;
  While[TrueQ[N[rad, 40] < N[s, 40]] && it < 10,
    (* iterate NUMERICALLY: intermediate placements are choices, and exact
       algebraic intermediate centers would drag ChartRadius into
       RootReduce algebra; only the final center is (coarsely) exact *)
    Module[{s2 = keff*N[rad, 40]/(1 + keff), g2, x2, rad2},
      g2 = s2/keff; x2 = N[xb, 40] + dir*g2;
      rad2 = ChartRadius[x2, allSings];
      If[TrueQ[N[rad2, 40] < N[rad, 40]], Break[]];
      s = s2; g = g2; xnew = x2; rad = rad2];
    it++];
  (* (G2) prev-disk clamp: the new chart's incoming match point
     xnew - dir*Radius/k is where the PREVIOUS solution gets evaluated; it
     must sit well inside the previous disk.  Enforce directly:
       h := g - Radius/k <= (9/10)*clamp*prevRad,
       clamp = Max[1/(2 keff), 1/k]
     (the 9/10 absorbs the coarse center rationalization below; the h < 0
     side needs no clamp: h < 0 forces g < prevRad/(k-1), so
     |h| <= R_new/k <= (prevRad + g)/k stays in the 1/(k-1) class).  h is
     increasing in g with slope in [1 - 1/k, 1 + 1/k] (ChartRadius is
     1-Lipschitz), so the unit-slope Newton step contracts the excess by a
     factor <= 1/k; shrinking g preserves (G1) because s = keff*g falls at
     least as fast as rad (keff >= 1).  Uses the line-CAPPED radius radC:
     that is the radius attached downstream, hence the actual match-point
     offset.  With keff = k the dominated march has h = 0 exactly and the
     clamp is inert (classic behavior). *)
  clamp = Max[1/(2*keff), 1/k];
  radC = Min[rad, lineCap];
  h = N[g - radC/k, 40];
  it = 0;
  While[TrueQ[h > (9/10)*clamp*N[prevRad, 40]] && it < 8,
    g = g - (h - (9/10)*clamp*N[prevRad, 40]);
    s = keff*g; xnew = xb + dir*g;
    rad = ChartRadius[xnew, allSings]; radC = Min[rad, lineCap];
    h = N[g - radC/k, 40]; it++];
  (* exact-ify the center COARSELY: the center is a placement CHOICE; a
     simple nearby rational keeps Indicial's exact algebra on small
     fractions (20-digit-denominator centers ground PossibleZeroQ into
     meprec storms).  Geometry self-corrects: the true radius and match
     point are recomputed from the actual center downstream.  Tolerance
     g/8 equals the former s/(8k) at keff = k. *)
  xnew = Module[{cand = Rationalize[N[xnew, 20], Abs[N[g, 20]]/8]},
    (* never land ON a singularity *)
    If[AnyTrue[allSings, TrueQ[PossibleZeroQ[RootReduce[# - cand]]] &],
      Rationalize[N[xnew, 20], Abs[N[g, 20]]/64], cand]];
  xnew];

(* singularMatchPoint: the incoming match point of a SINGULAR chart
   approached from the producing chart (prevCenter, prevRad) — the
   FixWithin clip (spec 2.3; old Transport.m:1124-1136).  The naive point
   z - dir*radTarget/k is anchored to the singular chart ALONE and can
   land far outside the producing disk (banana level 1: anchor radius
   1/46, singular radius ~0.45 — the naive point sits ~4 anchor radii
   BEHIND the anchor).  Required instead: a point of the intersection of
     (a) |m - prevCenter| <= margin*prevRad  (margin 9/10: strictly inside
         the producing disk, headroom absorbing the coarse rationalization
         below), and
     (b) dir*(z - m) > 0 and |m - z| < radTarget  (approach side of the
         singular chart, inside its disk).
   Since prevRad <= gap := dir*(z - prevCenter) (z is a singularity, so
   the producing disk never contains it), the intersection is nonempty iff
   gap < den := margin*prevRad + radTarget, and the BALANCED point
     m* = prevCenter + dir*gap*margin*prevRad/den
   equalizes the two normalized evaluation ratios
     |m - prevCenter|/(margin*prevRad) = |m - z|/radTarget = gap/den < 1,
   which (max of two opposed V-shaped functions of m) MINIMIZES
   max[ratioPrev/margin, ratioSing] over the intersection — the producing
   chart's truncation tail and the singular chart's matching/read tail are
   jointly as small as the geometry allows.  Regimes: tiny producing disk
   far inside a wide singular disk (banana: gap ~ prevRad ~ 0.022,
   radTarget ~ 0.45) gives BOTH ratios ~ gap/radTarget ~ 0.05; tiny
   singular disk reproduces the classic near-z handoff.  SegmentLine's
   cover test admits the singular chart only once |m - prevCenter| <=
   prevRad/k, which caps the executed ratios at 1/k (producing side) and
   ~(9/8)/(margin*k) (singular side) — the design classes.  Returns None
   when the intersection is empty: SegmentLine reads "not coverable from
   here, keep marching"; TransportLine/ValidatePlan read E8.  The point is
   a coarse simple rational like every match point: Rationalize at 1/8 of
   the least distance from m* to a constraint boundary, so (a) and (b)
   survive rationalization with 7/8 of their slack. *)
singularMatchPoint[prevCenter_, prevRad_, z_, radTarget_, dir_] := Module[
  {margin = 9/10, gap, den, mstar, bnd},
  gap = dir*(N[z, 40] - N[prevCenter, 40]);
  den = margin*N[prevRad, 40] + N[radTarget, 40];
  If[!TrueQ[0 < gap < den], Return[None]];
  mstar = N[prevCenter, 40] + dir*gap*margin*N[prevRad, 40]/den;
  bnd = Min[Min[margin*N[prevRad, 40], N[radTarget, 40]]*(den - gap),
    gap*N[radTarget, 40]]/den;
  Rationalize[N[mstar, 20], N[bnd, 20]/8]];

(* chartMatchPoint: THE incoming-match-point formula — the single source
   shared by TransportLine (evaluation), SegmentLine (cover target), and
   ValidatePlan (static audit), so plan and execution cannot drift.
   First chart (prev === None): the boundary anchor `from`, matched at
   t = 0.  Regular chart: center - dir*Radius/k, coarsely rationalized
   (algebraic radii otherwise force exact algebra through every
   evaluation).  Singular chart: the FixWithin clip against the producing
   chart (None when the handoff is impossible -> E8 at the consumers). *)
chartMatchPoint[None, chart_Association, from_, dir_, k_] := from;
chartMatchPoint[prev_Association, chart_Association, from_, dir_, k_] :=
  If[TrueQ[chart["Singular"]],
    singularMatchPoint[prev["Center"], prev["Radius"], chart["Center"],
      chart["Radius"], dir],
    Module[{raw = chart["Center"] - dir*chart["Radius"]/k},
      Rationalize[N[raw, 20], N[chart["Radius"], 20]/(100*k)]]];

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
  {sings, real, dir, k = cfg["DivisionOrder"], keff, charts, cur, interior,
   endpointSingular, all, guard = 0, lineCap, prevRad},
  sings = FindSingularities[sys];
  all = sings["All"]; real = sings["Real"];
  dir = Sign[to - from];
  If[dir === 0, err["E1", <|"From" -> from, "To" -> to, "Detail" -> "empty line"|>]];
  endpointSingular = AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - to]]] &];
  If[AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - from]]] &],
    err["E1", <|"From" -> from,
      "Detail" -> "transport FROM a singular point requires a singular boundary object (not supported in v1 marching start)"|>]];
  keff = stepDivisor[k];
  lineCap = 2*Abs[to - from];
  interior = Sort[Select[real, TrueQ[dir*(N[#, 40] - N[from, 40]) > 0] &&
      TrueQ[dir*(N[to, 40] - N[#, 40]) > 0] &], dir*N[#1 - #2, 40] < 0 &];
  (* the ANCHOR CHART sits exactly at `from` (regular by the check above):
     the boundary is matched at t = 0 — exact and perfectly conditioned —
     and the chart is plan-independent, so the lo/hi endpoint transports
     from one anchor share its PrepareChart AND SolveHomogeneous cache
     entries (one anchor solve per level instead of two near-copies). *)
  charts = {<|"Center" -> from, "Singular" -> False|>};
  cur = from;
  prevRad = Min[ChartRadius[from, all], lineCap];
  Module[{targets = Join[interior, {to}]},
    Do[Module[{target = targets[[ti]], targetSingular, radTarget, matchTarget},
      targetSingular = ti < Length[targets] || endpointSingular;
      radTarget = Min[ChartRadius[target, all], lineCap];
      While[True,
        guard++;
        If[guard > 500, err["E1", <|"Detail" -> "segmentation runaway"|>]];
        (* the point the LAST chart must reach: the target chart's incoming
           match point.  Singular target: the SHARED FixWithin clip
           (singularMatchPoint — the same formula TransportLine evaluates
           and ValidatePlan audits); it depends on the producing chart
           (cur, prevRad), hence recomputed per placement; None = the
           validity intersection is empty from here, not coverable yet.
           Regular target: the endpoint itself. *)
        matchTarget = If[targetSingular,
          singularMatchPoint[cur, prevRad, target, radTarget, dir],
          target];
        (* cover check on the LAST placed chart: matchTarget within
           prevRad/k keeps the handoff evaluation in the design-ratio
           class in BOTH charts (clipped point: <= 1/k producing-side,
           <= ~(9/8)/(margin k) singular-side) and no further regular
           chart is needed.  Only a REGULAR last chart may cover:
           consecutive singular targets can never cover each other
           (radius = distance to the nearest OTHER singularity forbids
           it), and a singular-chart-to-endpoint shortcut would put the
           final evaluation on an uncrossed branch — keep one regular
           chart in between. *)
        If[matchTarget =!= None && !TrueQ[Last[charts]["Singular"]] &&
           (TrueQ[Abs[N[matchTarget - cur, 40]] <= N[prevRad, 40]/k] ||
            TrueQ[dir*(N[cur, 40] - N[matchTarget, 40]) >= 0]),
          If[targetSingular,
            AppendTo[charts, <|"Center" -> target, "Singular" -> True|>]];
          Break[]];
        Module[{nxt},
          nxt = nextCenter[cur, real, dir, k, keff, all, prevRad,
            If[matchTarget === None, target - dir*radTarget/k, matchTarget],
            lineCap];
          AppendTo[charts, <|"Center" -> nxt, "Singular" -> False|>];
          cur = nxt;
          prevRad = Min[ChartRadius[nxt, all], lineCap]]];
      cur = target;
      If[targetSingular, prevRad = radTarget]],
      {ti, Length[targets]}]];
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

(* ValidatePlan: static pre-pass over the chart chain (TransportLine runs
   it before any solve).  Every chart's incoming match point — computed by
   the SAME chartMatchPoint formula the marching loop evaluates — must lie
   strictly inside the producing chart's disk (EvaluateLocalSolution's
   validity bound is |t| < Radius); a singular chart's point must in
   addition sit on the approach side strictly inside its own disk; the
   first chart's point is the boundary anchor plan["From"], inside the
   first disk.  A violation is a planner bug: E8 with the full chain
   geometry, raised before the expensive solves instead of deep inside
   them. *)
ValidatePlan[plan_Association] := Module[
  {charts = plan["Charts"], dir = plan["Direction"],
   k = cfg["DivisionOrder"], chain},
  chain = Map[<|"Center" -> #["Center"], "Radius" -> #["Radius"],
    "Singular" -> #["Singular"]|> &, charts];
  Do[Module[{chart = charts[[ci]], prev, mp, bad},
    prev = If[ci === 1, None, charts[[ci - 1]]];
    mp = chartMatchPoint[prev, chart, plan["From"], dir, k];
    If[mp === None,
      err["E8", <|"Chart" -> chart["Name"],
        "Detail" -> "singular-chart handoff impossible: the producing disk (margin 9/10) does not intersect the singular chart's approach interval",
        "PrevCenter" -> prev["Center"], "PrevRadius" -> prev["Radius"],
        "SingularCenter" -> chart["Center"],
        "SingularRadius" -> chart["Radius"],
        "AttemptedPoint" -> chart["Center"] - dir*chart["Radius"]/k,
        "Chain" -> chain|>]];
    bad = Which[
      ci === 1,
      !TrueQ[Abs[N[mp - chart["Center"], 40]] < N[chart["Radius"], 40]],
      True,
      !TrueQ[Abs[N[mp - prev["Center"], 40]] < N[prev["Radius"], 40]] ||
      (TrueQ[chart["Singular"]] &&
        !(TrueQ[dir*N[chart["Center"] - mp, 40] > 0] &&
          TrueQ[Abs[N[mp - chart["Center"], 40]] < N[chart["Radius"], 40]]))];
    If[bad,
      err["E8", <|"Chart" -> chart["Name"], "MatchPoint" -> mp,
        "PrevCenter" -> If[prev === None, None, prev["Center"]],
        "PrevRadius" -> If[prev === None, None, prev["Radius"]],
        "Center" -> chart["Center"], "Radius" -> chart["Radius"],
        "Chain" -> chain,
        "Detail" -> "incoming match point outside the producing chart's disk (planner bug)"|>]]],
    {ci, Length[charts]}];
  plan];

TransportLine[sys_Association, boundary_, plan_Association] := Module[
  {charts = plan["Charts"], dir = plan["Direction"], current, errAcc,
   req, expOrd = cfg["ExpansionOrder"], epsOrd = cfg["EpsilonOrder"],
   lastSingular = False, lastLS = None, lastChart = None, kept = {}},
  ValidatePlan[plan];
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
       previous chart's disk (step/radius <= 1/(1+k_eff) on dominated
       marches; EvaluateLocalSolution's radius assert is the loud
       backstop) — and that value is the t^0 Cauchy datum of ONE
       d-dimensional recursion (Solve`SolveValueRegular), replacing the
       d-column basis + MatchWeights + CombineLocalSolutions.  Singular
       charts and the first chart (anchor-only incoming data) keep the
       basis+matching path unchanged. *)
    valueMode = Environment["DE2_VALUE_TRANSPORT"] === "1" && ci > 1 &&
      !TrueQ[chart["Singular"]] &&
      TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
      (* conservative geometry pre-check: the center must sit WELL inside
         the previous object's disk; otherwise fall back to the basis path
         (a performance choice, not an ambiguity).  Margin: 9/10 under the
         classic coupled geometry (stepDivisor = DivisionOrder); under the
         wide-stride geometry the center handoff must ALSO keep the
         ratio^(ExpansionOrder+1) truncation tail below the 1e-14 design
         line — receding-leg gaps (up to ~0.8 of the previous radius)
         would otherwise pass the 9/10 check with ~1e-4 tails. *)
      Module[{margin = If[stepDivisor[cfg["DivisionOrder"]] < cfg["DivisionOrder"],
          Min[9/10, N[10^(-14/(cfg["ExpansionOrder"] + 1)), 30]], 9/10]},
        TrueQ[Abs[N[chart["Center"] - current["Center"], 30]] <
          margin*N[current["Radius"], 30]]];
    (* match point: the boundary anchor for the FIRST chart (the incoming
       object is only valid at its anchor); thereafter the shared
       chartMatchPoint formula — radius/k of THIS chart on the incoming
       side, or the FixWithin clip when THIS chart is singular
       (ValidatePlan certified every point against the producing disk) *)
    matchPt = chartMatchPoint[If[ci === 1, None, charts[[ci - 1]]], chart,
      plan["From"], dir, cfg["DivisionOrder"]];
    If[matchPt === None,
      err["E8", <|"Chart" -> chart["Name"],
        "Detail" -> "singular-chart handoff impossible (unreachable past ValidatePlan)"|>]];
    (* incoming value in the PREVIOUS object's chart coordinate: at the
       match point (basis path) or at this chart's center (value mode) *)
    tIn = If[valueMode, chart["Center"], matchPt] - current["Center"];
    Module[{prevEval, sigma},
      (* crossing: if the previous chart was singular and matchPt lies on its
         far side (sign of tIn relative to approach), apply the operator *)
      If[lastSingular && TrueQ[N[tIn, 30] dir > 0],
        sigma = DiffExp2`SectorSeries`ChartImSign[current];
        If[!MemberQ[{1, -1}, sigma],
          (* magnitude-aware gate: a multivalued sector blocks the crossing
             only if it carries MATERIAL content (box L2's apparent chart at
             7/11 has a t^(-1+eps) sector of physically zero weight; the
             syntactic tag test alone would E8 every such crossing).
             Ambiguous magnitudes stay material -> still loud. *)
          Module[{secs = current["Sectors"], scale, materialQ},
            scale = Max[1*^-300, Sequence @@ (Abs[N[#, 20]] & /@
              Select[Flatten[#["Coeffs"] & /@ secs], NumericQ])];
            materialQ[sec_] := Module[{mx},
              mx = Max[0, Sequence @@ (Abs[N[#, 20]] & /@
                Select[Flatten[sec["Coeffs"]], NumericQ])];
              !TrueQ[mx <= DiffExp2`Tolerances`Tol["LaurentLeadTol"]*scale]];
            If[AnyTrue[Select[secs, materialQ],
                !IntegerQ[#["a"]] || !zeroQ[#["b"]] || #["p"] > 0 &],
              err["E8", <|"Chart" -> chart["Name"],
                "Detail" -> "crossing a multivalued singular chart without a derivable Im-sign (missing DeltaPrescriptions)"|>],
              sigma = 1]]];
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
    (* probe on the INCOMING side (sign of matchPt - center keeps singular
       charts one-sided; the anchor chart, whose matchPt IS its center,
       probes the outgoing side) at the DESIGN evaluation ratio
       rho = Max[1/(2 k_step), 1/k]: downstream consumers — the next
       chart's match point under the stepDivisor geometry, the
       LineIntegral tile edges — evaluate THIS solution at |t| up to
       ~rho*Radius, so the honest truncation probe sits there.  The old
       half-incoming-match-point probe (Radius/(2k)) underestimates those
       tails by many decades once strides are wider than radius/k, and is
       identically 0 on the anchor chart. *)
    probeErrs = Module[{rho = Max[1/(2*stepDivisor[cfg["DivisionOrder"]]),
        1/cfg["DivisionOrder"]], sgn, raw},
      sgn = Sign[N[matchPt - chart["Center"], 30]];
      If[!MemberQ[{-1, 1}, sgn], sgn = dir];
      raw = rho*N[chart["Radius"], 20];
      SegmentErrorProbe[ls, sgn*Rationalize[raw, raw/50], couplingDepth]];
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
