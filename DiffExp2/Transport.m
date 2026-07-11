(* DiffExp2/Transport.m — segmentation, marching, eps-graded matching,
   crossing.  Spec: Docs/specs/Transport.md (binding); DECISIONS-M0.md.
   V1 notes (recorded): MatchWeights implements the honest Laurent-graded
   solve; the marching path applies tag recombination plus determinant-series
   epsilon-lattice saturation before matching.  The error probe is the
   full-vs-reduced evaluation at the outgoing match point. *)

BeginPackage["DiffExp2`Transport`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`", "DiffExp2`Solve`"}];

FindSingularities::usage = "FindSingularities[sys] gives <|\"All\", \"Real\", \"Factors\"|>: exact deduplicated roots of the matrix singular factors plus per-call ExtraSingularFactors.";
ChartRadius::usage = "ChartRadius[center, allSingularities] gives the exact/numeric complex-plane distance to the nearest OTHER singularity.";
SegmentLine::usage = "SegmentLine[sys, {from, to}] gives the SegmentPlan: charts, radii, match points, digit budget.";
TransportLine::usage = "TransportLine[sys, boundary, plan] runs the marching loop and returns the TransportResult.";
ValidatePlan::usage = "ValidatePlan[plan] statically audits the chart chain: every incoming match point (the shared chartMatchPoint formula) must lie inside both adjacent physical disks and their half-radius error-probe envelopes; singular handoffs must approach from the correct side. Loud E8 on violation; returns the plan.";
MatchWeights::usage = "MatchWeights[basisValues, incoming, label] solves the eps-graded (Laurent) weight system with loud residual asserts.";
ApplyCrossing::usage = "ApplyCrossing[ls, sigma] applies the crossing operator (phase times unipotent log-chain mixing) so the far side evaluates at positive chart coordinate.";
SegmentErrorProbe::usage = "SegmentErrorProbe[ls, tOut, couplingDepth] gives the full-vs-reduced evaluation error per eps order.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "Transport"|>, payload]];
cfg = DiffExp2`Config`CFG;
SetAttributes[plannerProfile, HoldRest];
plannerProfile[label_String, expr_] := If[
  Environment["DE2_PLANNER_PROFILE"] === "1",
  Module[{time, value},
    time = First@AbsoluteTiming[value = expr];
    Print["DE2PLAN ", label, " ", N[time, 8]];
    value],
  expr];
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

(* Private, default-off parity seam.  The constant right-frame
   normalization is already the production path for regular charts.  Its
   singular extension is kept gated while the real banana match and tagged
   resonant/log frames are exercised: it changes only the basis coordinates,
   never a tolerance, epsilon valuation, or transport point. *)
$enableSingularMatchPrecondition = False;
mwUseConstantMatchPreconditionQ[cs_Association] :=
  TrueQ[Lookup[cs["IndicialData"], "Regular", False]] ||
    TrueQ[$enableSingularMatchPrecondition];

(* Opt-in replay seam for expensive matches.  The path is read at the last
   point before any constant preconditioning or Laurent solve, after the
   actual LocalSolution basis has passed epsilon-lattice saturation.  A
   same-directory temporary plus overwrite rename keeps an interrupted write
   from destroying the previous replay.  With no environment path this is a
   pure no-op and creates neither globals nor files. *)
$matchFixtureSchema = "DiffExp2.MatchFixture/v1";
mwMaybeDumpMatchFixture[basis_List, Fmat_List, vIn_List, tLoc_,
    requiredTop_Integer, chart_Association, matchPoint_] := Module[
  {file = Quiet[Environment["DE2_MATCH_FIXTURE_FILE"]], dir, tmp, payload,
   wrote, renamed},
  If[!StringQ[file] || StringLength[StringTrim[file]] == 0,
    Return[Null, Module]];
  file = ExpandFileName[file];
  dir = DirectoryName[file];
  If[!DirectoryQ[dir],
    Quiet[Check[CreateDirectory[dir, CreateIntermediateDirectories -> True],
      err["E5", <|"Chart" -> Lookup[chart, "Name", "(unknown)"],
        "FixtureFile" -> file,
        "Detail" -> "could not create match-fixture directory"|>]]]];
  payload = <|
    "Schema" -> $matchFixtureSchema,
    "Label" -> Lookup[chart, "Name", "(unknown)"],
    "Chart" -> KeyTake[chart, {"Name", "Center", "Scale", "Radius",
      "LocalRadius", "MatchRadius", "Singular", "IncomingMatchPoint"}],
    "MatchPoint" -> matchPoint, "LocalCoordinate" -> tLoc,
    "RequiredTop" -> requiredTop,
    "Config" -> AssociationMap[Function[key, cfg[key]],
      {"WorkingPrecision", "ExpansionOrder", "EpsilonOrder",
       "DivisionOrder", "StepDivisionOrder", "RadiusOfConvergence"}],
    "Basis" -> basis, "F" -> Fmat, "V" -> vIn|>;
  Global`$DE2MatchFixture = payload;
  tmp = file <> ".tmp-" <> ToString[$ProcessID] <> "-" <>
    StringReplace[CreateUUID[], "-" -> ""] <> ".mx";
  wrote = Quiet[Check[
    DumpSave[tmp, Global`$DE2MatchFixture]; FileExistsQ[tmp], False]];
  Clear[Global`$DE2MatchFixture];
  If[!TrueQ[wrote],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    err["E5", <|"Chart" -> payload["Label"], "FixtureFile" -> file,
      "Detail" -> "could not write atomic match fixture"|>]];
  renamed = Quiet[Check[
    RenameFile[tmp, file, OverwriteTarget -> True]; True, False]];
  If[!TrueQ[renamed],
    If[FileExistsQ[tmp], Quiet[DeleteFile[tmp]]];
    err["E5", <|"Chart" -> payload["Label"], "FixtureFile" -> file,
      "Detail" -> "could not install atomic match fixture"|>]];
  file];

(* Stable true modulus shared by all coefficient/residual decisions. *)
numMag = DiffExp2`Tolerances`NumericMagnitude;
numMagBounds = DiffExp2`Tolerances`NumericMagnitudeBounds;
(* Exact point identity is used only after a numeric candidate prefilter.
   This is the singularity-dedup contract from Transport.md 2.1: a numeric
   key may nominate a pair, but only RootReduce may merge it.  The structural
   and rational exits cover nearly every ordinary center without entering an
   algebraic number field. *)
exactSamePointQ[a_, b_] := SameQ[a, b] || Module[{delta = Together[a - b]},
  delta === 0 || (!ratExprQ[delta] &&
    TrueQ[PossibleZeroQ[RootReduce[delta]]])];

(* At 70 probe digits, the 10^-50 relative window is wide enough to catch
   separately represented equal algebraic numbers, yet merely makes very
   close distinct roots exact-comparison CANDIDATES.  It can never merge
   them.  Keeping this threshold independent of WorkingPrecision also keeps
   pure planner cost independent of a recurrence precision such as WP1000. *)
numericSameCandidateQ[a_, b_] := Module[{na, nb, scale},
  If[SameQ[a, b], Return[True, Module]];
  {na, nb} = Quiet[Check[N[{a, b}, 70], {$Failed, $Failed}]];
  If[!And[NumericQ[na], NumericQ[nb]], Return[True, Module]];
  scale = Max[1, Abs[na], Abs[nb]];
  TrueQ[Abs[na - nb] <= 10^-50 scale]];

exactDeduplicatePoints[points_List] := Module[{out = {}},
  Do[If[!AnyTrue[Select[out, numericSameCandidateQ[p, #] &],
        exactSamePointQ[p, #] &], AppendTo[out, p]], {p, points}];
  out];

exactPointMemberQ[p_, points_List] := MemberQ[points, p] ||
  AnyTrue[Select[points, numericSameCandidateQ[p, #] &],
    exactSamePointQ[p, #] &];

(* Order exact real planner points without imposing a fixed absolute
   resolution.  Separated points take the cheap endpoint-first numeric path;
   only an ambiguous pair enters RootReduce/sign isolation.  In particular,
   1 and 1+10^-80 must not become equal merely because a 40-digit rendering
   rounds both endpoints to the same decimal. *)
pointOrderSign[a_, b_, digits_Integer:40] := Module[
  {delta = Together[a - b], na, nb, nd, scale, reduced, sgn},
  If[delta === 0, Return[0, Module]];
  If[ratExprQ[delta] && NumericQ[delta], Return[Sign[delta], Module]];
  {na, nb} = Quiet[Check[N[{a, b}, digits + 20], {$Failed, $Failed}]];
  If[And[NumericQ[na], NumericQ[nb]],
    nd = na - nb;
    scale = Max[1, Abs[na], Abs[nb]];
    If[TrueQ[Abs[nd] > 10^-digits scale],
      sgn = Quiet[Sign[Re[nd]]];
      If[MemberQ[{-1, 1}, sgn], Return[sgn, Module]]]];
  reduced = Quiet[RootReduce[delta]];
  If[reduced === 0, Return[0, Module]];
  sgn = Quiet[Sign[reduced]];
  If[MemberQ[{-1, 1}, sgn], Return[sgn, Module]];
  nd = Quiet[Check[N[reduced, digits], $Failed]];
  sgn = If[NumericQ[nd], Quiet[Sign[Re[nd]]], $Failed];
  If[MemberQ[{-1, 1}, sgn], sgn,
    err["E8", <|"PointA" -> a, "PointB" -> b,
      "ExactDifference" -> reduced,
      "Detail" -> "could not order exact real planner points"|>]]];

(* Numerically compare planner points by evaluating the endpoints first.
   The previous unconditional RootReduce[a-b] made an O(n^2) spacing scan
   build a fresh degree-8 number field for every pair.  Rational differences
   remain exact (pinning the 10^-100 regression); algebraic expressions are
   evaluated separately, which avoids catastrophic symbolic cancellation.
   Exact zero decisions never use this helper. *)
numericDistance[a_, b_, digits_Integer:40] := Module[
  {delta = Together[a - b], value, reduced},
  If[delta === 0, Return[0, Module]];
  If[ratExprQ[delta] && NumericQ[delta],
    Return[Abs[N[delta, digits]], Module]];
  value = Quiet[Check[Block[{$MaxExtraPrecision = Max[1000,
        2*DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]},
      Abs[N[a, digits + 20] - N[b, digits + 20]]], $Failed]];
  If[value === $Failed || !NumericQ[value],
    err["E8", <|"PointA" -> a, "PointB" -> b,
      "ExactDifference" -> delta,
      "Detail" -> "could not resolve an exact algebraic planner distance"|>]];
  (* Rare close-cancellation fallback.  Endpoint-first evaluation is fast for
     separated points, but two distinct algebraic numbers can agree through
     all probe digits.  Low output precision nominates that pair for the
     expensive exact reduction; it does not make a zero decision. *)
  If[FreeQ[{a, b}, _?InexactNumberQ] &&
      (TrueQ[PossibleZeroQ[value]] || !NumericQ[Precision[value]] ||
        TrueQ[Precision[value] < digits]),
    reduced = Quiet[RootReduce[delta]];
    If[reduced === 0, Return[0, Module]];
    value = Quiet[Check[Block[{$MaxExtraPrecision = Max[1000,
          2*DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]},
        N[Abs[reduced], digits]], $Failed]];
    If[value === $Failed || !NumericQ[value],
      err["E8", <|"PointA" -> a, "PointB" -> b,
        "ExactDifference" -> reduced,
        "Detail" -> "could not resolve a close exact algebraic planner distance"|>]]];
  N[value, digits]];

(* Pipeline numeric handoff: evaluated VALUES entering the matching seam
   numericize at 2x WP — the pinned $InputPrecisionFactor policy and the
   runner's boundary/level-handoff policy (M5c).  The extra headroom is
   essential for ill-conditioned honest Laurent bases: matching may lose
   many digits while still meeting MatchTol at the configured precision.
   applied inside the marching loop.  EvaluateLocalSolution itself stays
   exact-in/exact-out (a module contract, pinned by test_sectorseries);
   but exact evaluation outputs fed into MatchWeights/Combine grind
   exact-giant field arithmetic (measured d = 7 hot spot: 430/560 exact
   coefficients, MatchWeights 3.0 s -> see Docs/SpeedIdeas.md §2).
   Values only; tags, windows and structural decisions untouched;
   symbols pass through N unchanged; idempotent on bignums. *)
numHandoff[x_] := N[x,
  DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]];

(* A value solve cannot recover significance already lost while evaluating
   the preceding truncated solution at the next chart center.  In
   particular N[x, 2 WP] does not add digits to an inexact x.  Require the
   incoming Cauchy data's tracked uncertainty to sit below the same relative
   scale used by ODEResidualCheck, with the repository safety margin.  When
   this preflight fails TransportLine chooses the ordinary match-point basis
   path before solving; this is a conditioning choice between equivalent
   transports, not a fallback from a recurrence failure. *)
valueHandoffAccurateQ[vals_List] := Module[
  {rtol = DiffExp2`Tolerances`Tol["ResidTol"],
   guard = DiffExp2`Tolerances`$SafetyDigits, nums},
  nums = Cases[Flatten[Table[
      Table[esCoeff[v, k], {k, esMin[v], esCM[v]}], {v, vals}]],
    _?InexactNumberQ, Infinity];
  AllTrue[nums, Function[z, Module[{acc = Accuracy[z], uncertainty, scale},
    If[acc === Infinity, Return[True, Module]];
    If[!NumericQ[acc], Return[False, Module]];
    (* Exact decade: cannot machine-underflow for high-precision inputs. *)
    uncertainty = 10^-Floor[acc];
    scale = Max[1, numMag[z, 20]];
    TrueQ[uncertainty <= rtol*scale/10^guard]]]]];

(* A center handoff evaluates a truncated Taylor solution farther from its
   origin than the ordinary match point.  Precision metadata cannot see
   that truncation error, so admit it only when the geometric tail proxy is
   two decades below the structural Laurent floor. *)
valueCenterMargin[expansionOrder_Integer] := Min[9/10,
  N[(DiffExp2`Tolerances`Tol["LaurentLeadTol"]/100)^
    (1/(expansionOrder + 1)), 30]];

(* The classic +1/k,-1/k handoff keeps basis matching well conditioned, but
   its adjacent CENTERS can be much farther apart than a truncated solution
   may safely be evaluated for value-vector propagation.  In value mode only,
   refine regular receivers until both contracts hold:

     (1) the receiver center is strictly inside the producing chart's
         truncation-tail margin; and
     (2) a balanced ordinary match point lies in the 1/k class of BOTH disks,
         so the certified basis path remains available if the value datum
         later fails its significance check.

   Singular receivers are never changed: their canonical FixWithin point,
   sector decomposition, and crossing semantics remain owned by the original
   planner.  Midpoints of exact centers are exact (RootReduce); the inexact
   branch is rationalized before it can become a chart center. *)
valueRefineRegularChain[charts_List, all_List, dir_, k_Integer,
    lineCap_] := Module[
  {margin = valueCenterMargin[cfg["ExpansionOrder"]], safety = 99/100,
   radius, exactCenter, fallbackPoint, hopSafeQ, refinePair, out},
  radius[c_] := Min[ChartRadius[c, all], lineCap];
  exactCenter[z_] := If[FreeQ[z, _?InexactNumberQ], RootReduce[z],
    Rationalize[N[z, 40], Max[10^-40, Abs[N[z, 40]]*10^-35]]];
  fallbackPoint[left_, right_] := Module[
    {cl = left["Center"], cr = right["Center"], rl, rr},
    rl = radius[cl]; rr = radius[cr];
    RootReduce[(cl*rr + cr*rl)/(rl + rr)]];
  hopSafeQ[left_, right_] := Module[
    {cl = left["Center"], cr = right["Center"], rl, rr, gap},
    rl = radius[cl]; rr = radius[cr];
    gap = numericDistance[cr, cl, 40];
    TrueQ[dir*pointOrderSign[cr, cl] > 0] &&
      TrueQ[gap < safety*margin*N[rl, 40]] &&
      (* This is stronger than the half-disk condition audited by
         ValidatePlan and pins the fallback to the established design class. *)
      TrueQ[gap <= N[(rl + rr)/k, 40]]];
  refinePair[left_, right_, depth_Integer] := Module[
    {mid, bridge, tail, p},
    If[hopSafeQ[left, right],
      p = fallbackPoint[left, right];
      Return[{Join[right, <|"IncomingMatchPoint" -> p,
        "SymmetricMatch" -> False|>]}, Module]];
    If[depth >= 64,
      err["E8", <|"LeftCenter" -> left["Center"],
        "RightCenter" -> right["Center"], "Margin" -> margin,
        "Detail" -> "value-aware regular chart refinement exceeded 64 bisections"|>]];
    mid = exactCenter[(left["Center"] + right["Center"])/2];
    If[TrueQ[PossibleZeroQ[RootReduce[mid - left["Center"]]]] ||
        TrueQ[PossibleZeroQ[RootReduce[mid - right["Center"]]]],
      err["E8", <|"LeftCenter" -> left["Center"],
        "RightCenter" -> right["Center"], "Midpoint" -> mid,
        "Detail" -> "value-aware refinement could not construct a distinct exact midpoint"|>]];
    bridge = <|"Center" -> mid, "Singular" -> False|>;
    tail = refinePair[left, bridge, depth + 1];
    Join[tail, refinePair[Last[tail], right, depth + 1]]];
  out = {First[charts]};
  Do[
    If[TrueQ[right["Singular"]],
      AppendTo[out, right],
      out = Join[out, refinePair[Last[out], right, 0]]],
    {right, Rest[charts]}];
  out];

(* EvaluateLocalSolution's origin rule is intentionally strict for tagged
   powers.  Admit the first anchor to value mode only when its incoming object
   is manifestly center-evaluable: nonnegative integer Taylor powers, no
   epsilon-dependent exponent, and no logarithm.  Plain boundary matrices are
   wrapped as the single (0,0,0) sector and therefore qualify. *)
centerValueDatumQ[ls_Association] := KeyExistsQ[ls, "Sectors"] &&
  AllTrue[ls["Sectors"], Function[sec,
    IntegerQ[sec["a"]] && TrueQ[sec["a"] >= 0] &&
      zeroQ[sec["b"]] && sec["p"] === 0]];
centerValueDatumQ[_] := False;

(* ---- 2.1 singularities ---- *)

projectComplexRoots[all_List, real_List] := Module[{data, projected},
  (* Port the old DiffExp ghost-waypoint construction faithfully, but keep
     the true complex roots as the sole convergence-radius alphabet.  For a
     conjugate pair z=re+/-i h this contributes the real waypoints
     re-h,re,re+h unless another singularity's real projection already lies
     in the corresponding open strip.  These are REGULAR chart centers, not
     fake poles; they stabilize the real march around a nearby complex pair
     and carry no indicial/branch semantics of their own. *)
  data = Map[Function[root, Module[{z = root, re, im, h},
      (* Keep exact algebraic projections as unevaluated field expressions.
         RootReduce here is unnecessary: these are numericized/rationalized
         before becoming regular chart centers.  For conjugate Root pairs,
         Wolfram canonical arithmetic makes the resulting re,h expressions
         structurally identical, so the cheap DeleteDuplicates below removes
         them before the exact-candidate fallback. *)
      re = (z + Conjugate[z])/2;
      im = (z - Conjugate[z])/(2 I);
      h = If[exactPointMemberQ[z, real], 0,
        If[TrueQ[Re[N[im, 70]] < 0], -im, im]];
      (* Independent numerical evaluations of conjugate Root objects can
         leave a harmless 10^-p imaginary center in an exactly real re.
         Strip that numerical noise only in the occupancy key; the exact re
         expression above is retained in the returned waypoint. *)
      {z, re, h, Re[N[re, 70]], Abs[N[h, 70]]}]], all];
  data = DeleteDuplicatesBy[data, {#[[2]], #[[3]]} &];
  projected = Flatten[Map[Function[row, Module[{re = row[[2]], h = row[[3]],
        reN = row[[4]], hN = row[[5]], leftOccupied, rightOccupied},
      If[TrueQ[PossibleZeroQ[h]], {re},
        leftOccupied = AnyTrue[data,
          Function[other, TrueQ[reN - hN < other[[4]] < reN]]];
        rightOccupied = AnyTrue[data,
          Function[other, TrueQ[reN < other[[4]] < reN + hN]]];
        Join[If[leftOccupied, {}, {re - h}], {re},
          If[rightOccupied, {}, {re + h}]]]]], data]];
  Sort[exactDeduplicatePoints[DeleteDuplicates[projected]],
    pointOrderSign[#1, #2, 70] < 0 &]];

FindSingularities[sys_Association] := Module[
  {var = sys["Variable"], facs, extra, all, roots, real},
  facs = Lookup[sys, "SingularFactors", {}];
  extra = Lookup[sys, "ExtraSingularFactors", {}];
  facs = DeleteDuplicates[Join[facs, extra]];
  roots = plannerProfile["FindRoots",
    Map[# -> DeleteDuplicates[var /. Solve[# == 0, var]] &, facs]];
  all = plannerProfile["DeduplicateRoots",
    exactDeduplicatePoints[Flatten[Last /@ roots]]];
  real = plannerProfile["ClassifyRealRoots",
    Select[all, zeroQ[Im[RootReduce[#]]] &]];
  <|"All" -> all,
    "Real" -> real,
    "Projected" -> plannerProfile["ProjectComplexRoots",
      projectComplexRoots[all, real]],
    "Factors" -> Association[roots]|>];

chartPrescriptions[center_, var_Symbol] := Module[
  {entries = cfg["DeltaPrescriptions"], t = Global`t, out = {}},
  Do[Module[{poly = entry[[1]], sign = entry[[2]], shifted, mult, lead,
      leadSign},
    shifted = Together[poly /. s_Symbol /;
        SymbolName[s] === SymbolName[var] :> center + t];
    If[zeroQ[shifted /. t -> 0],
      If[!PolynomialQ[shifted, t],
        err["E8", <|"Center" -> center, "Factor" -> poly,
          "Detail" -> "prescription did not become a polynomial in the chart coordinate"|>]];
      mult = Exponent[shifted, t, Min];
      lead = RootReduce[Coefficient[shifted, t, mult]];
      leadSign = Sign[lead];
      If[!MemberQ[{1, -1}, leadSign],
        err["E8", <|"Center" -> center, "Factor" -> poly,
          "LeadingCoefficient" -> lead,
          "Detail" -> "could not derive a real nonzero leading-coefficient sign for DeltaPrescriptions"|>]];
      AppendTo[out, <|"Factor" -> poly, "Sign" -> sign,
        "Multiplicity" -> mult, "LeadingCoeffSign" -> leadSign|>]]],
    {entry, entries}];
  out];

ChartRadius[center_, all_List] := Module[
  {others, diffs, digits = cfg["WorkingPrecision"]},
  (* Exact membership is candidate-filtered: far algebraic roots never need
     a number-field reduction, while a separately represented equal point is
     still excluded only after exact RootReduce confirmation. *)
  others = Select[all, !exactPointMemberQ[#, {center}] &];
  If[others === {}, Infinity,
    (* Transport.md 2.2 permits a >=WP numerical radius when the exact
       algebraic modulus is expensive.  Rational/radical distances stay
       exact (pinning the classic +/-1/k formulas); a Root-valued distance is
       evaluated from its exact endpoints at WP instead of constructing the
       often enormous number field of Abs[root-center].  Radius is geometry
       metadata only: regular chart Scale continues to come from the small
       exact/rationalized projection alphabet. *)
    diffs = Map[If[FreeQ[{#, center}, _Root | _?InexactNumberQ],
      RootReduce[Abs[# - center]], numericDistance[#, center, digits]] &,
      others];
    Min[diffs]]];

simpleProjectionWaypoints[projected_List, real_List, lineCap_] := Module[
  {pts, numericProjected = N[projected, 70]},
  pts = Map[Function[p,
    If[exactPointMemberQ[p, real], p,
      Module[{seps, tol, cand},
        (* projected is already exactly deduplicated.  This scan chooses only
           a rationalization tolerance, so one cached numeric vector is the
           faithful old DiffExp operation and avoids O(n^2) RootReduce. *)
        seps = Select[Abs[N[p, 70] - #] & /@ numericProjected, # > 0 &];
        tol = Min[N[lineCap, 40], If[seps === {}, N[lineCap, 40], Min[seps]]]/64;
        cand = Rationalize[N[p, 30], tol];
        If[exactPointMemberQ[cand, real],
          Rationalize[N[p, 30], tol/64], cand]]]], projected];
  Sort[exactDeduplicatePoints[DeleteDuplicates[pts]],
    pointOrderSign[#1, #2, 70] < 0 &]];

projectionRadius[center_, projected_List, lineCap_] := Module[{others},
  others = Select[projected,
    !TrueQ[PossibleZeroQ[RootReduce[# - center]]] &];
  If[others === {}, lineCap,
    Min[lineCap, Min[RootReduce[Abs[# - center]] & /@ others]]]];

(* Matching and the error probe share one hard geometric envelope.  Keep the
   old projected-geometry 1/k point whenever it lies within half of the true
   complex convergence radius; k=2 plus rationalized algebraic projections
   can differ by a tiny amount, in which case the true half-radius wins. *)
matchOffset[matchRadius_, trueRadius_, k_Integer] :=
  Block[{$MaxExtraPrecision = Max[1000,
      2*DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]},
  Quiet[Module[{projected = RootReduce[matchRadius/k],
      cap = RootReduce[trueRadius/2], delta, numericDelta},
    delta = RootReduce[cap - projected];
    Which[
      delta === 0, projected,
      TrueQ[Quiet[Sign[delta]] === 1], projected,
      TrueQ[Quiet[Sign[delta]] === -1], cap,
      True,
      numericDelta = Quiet[Check[Block[{
          $MaxExtraPrecision = Max[1000,
            2*DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]]},
          N[delta, DiffExp2`Tolerances`$InputPrecisionFactor*
            cfg["WorkingPrecision"]]], $Failed]];
      Which[TrueQ[numericDelta > 0], projected,
        TrueQ[numericDelta < 0], cap,
        True, err["E8", <|"MatchRadius" -> matchRadius,
          "TrueRadius" -> trueRadius, "DivisionOrder" -> k,
          "Difference" -> delta,
          "Detail" -> "could not decide projected 1/k versus true half-radius match bound"|>]]]]]];

(* ---- 2.4 GetCPL/GetCPR geometry ---- *)

(* Exact port of the old non-Mobius GetCPL/GetCPR construction.  Let xb be
   the current chart center and r its projected-real geometry radius.  The
   outgoing point is q=xb+dir r/k.  The next center is chosen so q is also
   exactly centerNext-dir rNext/k.  Thus every regular-to-regular match is
   evaluated at +1/k in the producing segment and -1/k in the receiving
   segment, instead of accumulating one-sided, ill-conditioned matches. *)
classicNextCenter[xb_, projected_List, dir_, k_Integer, lineCap_,
    matchTarget_] := Module[
  {r, q, left, right, steps, g, xnew},
  r = projectionRadius[xb, projected, lineCap];
  q = xb + dir*r/k;
  left = Select[projected, pointOrderSign[#, q] < 0 &];
  right = Select[projected, pointOrderSign[#, q] > 0 &];
  left = If[left === {}, -Infinity, Last[left]];
  right = If[right === {}, Infinity, First[right]];
  steps = If[dir > 0,
    DeleteCases[{If[right === Infinity, Infinity, (right - q)/(k + 1)],
      If[left === -Infinity, Infinity, (q - left)/(k - 1)]}, Infinity],
    DeleteCases[{If[left === -Infinity, Infinity, (q - left)/(k + 1)],
      If[right === Infinity, Infinity, (right - q)/(k - 1)]}, Infinity]];
  g = If[steps === {}, lineCap/k, Min[steps]];
  xnew = q + dir*g;
  If[TrueQ[dir*pointOrderSign[matchTarget, xb] > 0] &&
      TrueQ[dir*pointOrderSign[xnew, matchTarget] >= 0],
    xnew = matchTarget];
  (* All non-real projections were simplified to nearby rationals above;
     retain exact real singularities verbatim and keep ordinary centers in
     the same small exact field as the old predivision formulas. *)
  If[FreeQ[xnew, _?InexactNumberQ], RootReduce[xnew],
    Rationalize[N[xnew, 30], Max[10^-30, Abs[N[g, 30]]/256]]]];

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
  {margin = 9/10, gap, den, offset, bnd, offsetRat, exactGap, exactDen},
  If[dir*pointOrderSign[z, prevCenter] <= 0, Return[None, Module]];
  If[AllTrue[{prevCenter, prevRad, z, radTarget},
      ratExprQ[#] && NumericQ[#] &],
    exactGap = dir*(z - prevCenter);
    exactDen = margin*prevRad + radTarget;
    If[TrueQ[0 < exactGap < exactDen],
      Return[Together[prevCenter +
        (z - prevCenter)*margin*prevRad/exactDen], Module],
      Return[None, Module]]];
  gap = numericDistance[z, prevCenter, 40];
  den = margin*N[prevRad, 40] + N[radTarget, 40];
  If[!TrueQ[0 < gap < den], Return[None]];
  offset = gap*margin*N[prevRad, 40]/den;
  bnd = Min[Min[margin*N[prevRad, 40], N[radTarget, 40]]*(den - gap),
    gap*N[radTarget, 40]]/den;
  (* Rationalize the small RELATIVE offset before adding it to an exact
     center.  Rationalizing an absolute number such as 1+10^-80 at fixed
     precision erases the displacement and can make the planner loop. *)
  offsetRat = Rationalize[N[offset, 20], N[bnd, 20]/8];
  If[FreeQ[prevCenter, _?InexactNumberQ],
    RootReduce[prevCenter + dir*offsetRat],
    N[prevCenter, 40] + dir*N[offset, 40]]];

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
  If[KeyExistsQ[chart, "IncomingMatchPoint"],
    chart["IncomingMatchPoint"],
    If[TrueQ[chart["Singular"]],
    singularMatchPoint[prev["Center"], prev["Radius"], chart["Center"],
      chart["Radius"], dir],
    chart["Center"] - dir*matchOffset[
      Lookup[chart, "MatchRadius", chart["Radius"]], chart["Radius"], k]]];

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
  {sings, real, projected, dir, k = cfg["DivisionOrder"], charts,
   cur, interior, endpointSingular, all, guard = 0, lineCap, prevRad,
   prevMatchRad, var = sys["Variable"]},
  sings = plannerProfile["FindSingularities", FindSingularities[sys]];
  all = sings["All"]; real = sings["Real"];
  dir = Sign[to - from];
  If[dir === 0, err["E1", <|"From" -> from, "To" -> to, "Detail" -> "empty line"|>]];
  endpointSingular = AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - to]]] &];
  If[AnyTrue[real, TrueQ[PossibleZeroQ[RootReduce[# - from]]] &],
    err["E1", <|"From" -> from,
      "Detail" -> "transport FROM a singular point requires a singular boundary object (not supported in v1 marching start)"|>]];
  lineCap = 2*Abs[to - from];
  projected = plannerProfile["SimpleProjectionWaypoints",
    simpleProjectionWaypoints[sings["Projected"], real, lineCap]];
  interior = Sort[Select[projected,
      dir*pointOrderSign[#, from] > 0 &&
      dir*pointOrderSign[to, #] > 0 &],
    dir*pointOrderSign[#1, #2] < 0 &];
  (* the ANCHOR CHART sits exactly at `from` (regular by the check above):
     the boundary is matched at t = 0 — exact and perfectly conditioned —
     and the chart is plan-independent, so the lo/hi endpoint transports
     from one anchor share its PrepareChart AND SolveHomogeneous cache
     entries (one anchor solve per level instead of two near-copies). *)
  charts = {<|"Center" -> from, "Singular" -> False|>};
  cur = from;
  prevRad = Min[ChartRadius[from, all], lineCap];
  prevMatchRad = projectionRadius[from, projected, lineCap];
  Module[{targets = Join[Map[Function[q, {q, AnyTrue[real,
          Function[r, TrueQ[PossibleZeroQ[RootReduce[r - q]]]]]}], interior],
        {{to, endpointSingular}}]},
    Do[Module[{target = targets[[ti, 1]],
        targetSingular = targets[[ti, 2]], forcedWaypoint,
        radTarget, matchTarget},
      forcedWaypoint = ti < Length[targets] ||
        (!targetSingular && AnyTrue[projected,
          TrueQ[PossibleZeroQ[RootReduce[# - target]]] &]);
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
          If[forcedWaypoint,
            target - dir*matchOffset[
              projectionRadius[target, projected, lineCap], radTarget, k],
            target]];
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
           (TrueQ[numericDistance[matchTarget, cur, 40] <=
                N[matchOffset[prevMatchRad, prevRad, k], 40]] ||
            TrueQ[dir*pointOrderSign[cur, matchTarget] >= 0]) &&
           (!targetSingular ||
             TrueQ[numericDistance[matchTarget, target, 40] <= N[radTarget, 40]/2]),
          If[targetSingular,
            AppendTo[charts, <|"Center" -> target, "Singular" -> True,
              "IncomingMatchPoint" -> matchTarget|>],
            If[forcedWaypoint &&
                !TrueQ[PossibleZeroQ[RootReduce[Last[charts]["Center"] - target]]],
              AppendTo[charts, <|"Center" -> target, "Singular" -> False,
                "IncomingMatchPoint" -> matchTarget|>]]];
          Break[]];
        Module[{nxt, incoming, producerIncoming, receiverIncoming,
            nxtMatchRad, nxtTrueRad, symmetric, receiverSafeQ,
            clipRefinements = 0},
          nxt = classicNextCenter[cur, projected, dir, k, lineCap,
            If[matchTarget === None, target - dir*radTarget/k, matchTarget]];
          nxtMatchRad = projectionRadius[nxt, projected, lineCap];
          nxtTrueRad = Min[ChartRadius[nxt, all], lineCap];
          producerIncoming = cur + dir*matchOffset[prevMatchRad, prevRad, k];
          receiverIncoming = nxt - dir*matchOffset[
            nxtMatchRad, nxtTrueRad, k];
          receiverSafeQ[] :=
            TrueQ[numericDistance[receiverIncoming, cur, 40] <=
              N[matchOffset[prevMatchRad, prevRad, k], 40]] &&
            TrueQ[numericDistance[receiverIncoming, cur, 40] < N[prevRad, 40]] &&
            TrueQ[numericDistance[receiverIncoming, nxt, 40] <= N[nxtTrueRad, 40]/2];
          (* classicNextCenter may clip its natural symmetric center to a
             nearby target (notably the mandatory regular chart immediately
             after a singular crossing).  The producer's old +R/k point then
             belongs to the UNCLIPPED chart and can sit near/outside the new
             receiver disk.  Recompute from the actual receiver's -R/k side.
             If that point is not yet in the producer's design disk, insert a
             midpoint chart and repeat; never store a stale pre-clip point. *)
          If[TrueQ[PossibleZeroQ[Together[
                producerIncoming - receiverIncoming]]] && receiverSafeQ[],
            incoming = producerIncoming,
            While[!receiverSafeQ[] && clipRefinements < 64,
              nxt = RootReduce[(cur + nxt)/2];
              nxtMatchRad = projectionRadius[nxt, projected, lineCap];
              nxtTrueRad = Min[ChartRadius[nxt, all], lineCap];
              receiverIncoming = nxt - dir*matchOffset[
                nxtMatchRad, nxtTrueRad, k];
              clipRefinements++];
            If[!receiverSafeQ[],
              err["E8", <|"CurrentCenter" -> cur,
                "CurrentRadius" -> prevRad,
                "CurrentMatchRadius" -> prevMatchRad,
                "CandidateCenter" -> nxt,
                "CandidateRadius" -> nxtTrueRad,
                "CandidateMatchRadius" -> nxtMatchRad,
                "CandidateIncomingPoint" -> receiverIncoming,
                "Detail" -> "could not place a clipped receiver match point inside both adjacent design disks"|>]];
            incoming = receiverIncoming];
          symmetric = TrueQ[PossibleZeroQ[Together[
              (incoming - cur)/prevMatchRad - dir/k]]] &&
            TrueQ[PossibleZeroQ[Together[
              (incoming - nxt)/nxtMatchRad + dir/k]]];
          AppendTo[charts, <|"Center" -> nxt, "Singular" -> False,
            "IncomingMatchPoint" -> incoming,
            "SymmetricMatch" -> symmetric|>];
          cur = nxt;
          prevRad = Min[ChartRadius[nxt, all], lineCap];
          prevMatchRad = nxtMatchRad]];
      cur = target;
      prevRad = radTarget;
      prevMatchRad = projectionRadius[target, projected, lineCap]],
      {ti, Length[targets]}]];
  (* Value propagation needs Cauchy data at the next chart CENTER.  The
     classic chart chain was designed for +/-1/k matching points instead, so
     refine it only under the prototype flag.  Flag-off plans remain exactly
     the pre-existing plans, including centers and match points. *)
  If[Environment["DE2_VALUE_TRANSPORT"] === "1",
    charts = valueRefineRegularChain[charts, all, dir, k, lineCap]];
  (* attach radii, match points, names; radii capped at line scale
     (a validity bound: capping is conservative; uncapped Infinity poisons
     the match-point arithmetic on singularity-free systems) *)
  charts = MapIndexed[Module[{c = #1, rad, matchRad, scale, localRad,
      roc = cfg["RadiusOfConvergence"]},
    rad = Min[ChartRadius[c["Center"], all], 2*Abs[to - from]];
    matchRad = projectionRadius[c["Center"], projected,
      2*Abs[to - from]];
    (* Faithful affine part of old GetLineRescaled: projected geometry sets
       the physical scale, while the true complex distance remains the
       validity bound.  In local t, x = center + scale t, so high Taylor
       coefficients stay O(1) and the true convergence radius is rad/scale.
       RadiusOfConvergence is a coordinate normalization, not a change to
       the physical chart coverage or match points. *)
    scale = RootReduce[matchRad/roc];
    localRad = RootReduce[rad/scale];
    Join[c, <|"Radius" -> rad, "ChartVar" -> Global`t,
      "MatchRadius" -> matchRad, "Scale" -> scale,
      "LocalRadius" -> localRad,
      "Name" -> "seg" <> ToString[First[#2]] <> "@" <>
        ToString[N[c["Center"], 6]],
      "Prescriptions" -> chartPrescriptions[c["Center"], var]|>]] &, charts];
  <|"From" -> from, "To" -> to, "Direction" -> dir,
    "Charts" -> charts, "SegmentCount" -> Length[charts],
    "EndpointIsSingular" -> endpointSingular,
    "DigitsNeeded" -> DigitBudget[cfg["AccuracyGoal"], Length[charts]],
    "Singularities" -> Append[sings, "ProjectionWaypoints" -> projected]|>];

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
   {B_2, ...} up to the family's recorded degeneracy depth.

   SAFETY GATE: an uncompensated CASE-P quotient must not be divided by eps
   again.  Solve certifies its joint polar cancellation with
   Diagnostics["PseudoCollisionsCompensated"] -> True.  Once certified,
   distinct same-a eps-degeneracy is a separate obligation and still needs
   this divided-difference recombination.  Banana has exactly this mixed
   case: an n=1 compensated hit alongside a same-a pair whose unrecombined
   columns are catastrophically ill-conditioned. *)

recombineDegenerate[cs_, basis_List, specs_List, diagnostics_:<||>] := Module[
  {degs, newBasis = basis, width,
   hits = If[AssociationQ[diagnostics] &&
       KeyExistsQ[diagnostics, "PseudoCollisionsHit"],
     diagnostics["PseudoCollisionsHit"], Missing["NotAvailable"]],
   compensated = AssociationQ[diagnostics] &&
     TrueQ[Lookup[diagnostics, "PseudoCollisionsCompensated", False]]},
  (* Missing/malformed diagnostics are not evidence that CASE-P was
     discharged.  A nonempty hit list is safe only with Solve's explicit
     joint-compensation certificate. *)
  If[!ListQ[hits] || (hits =!= {} && !compensated),
    Return[newBasis]];
  degs = DiffExp2`Indicial`EpsDegenerateFamilies[cs["IndicialData"]];
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
mwScale[c_List] := Max[0,
  Sequence @@ (If[NumericQ[#], numMag[#, 10], Nothing] & /@ c)];

(* A frame whose stored central values are ALL zero has no relative scale.
   NumericallyZeroQ intentionally refuses scale == 0, since an arbitrary
   caller has supplied no meaningful comparison.  Matching does have an
   absolute contract in this one case: use unit scale so a well-resolved
   centered-zero frame can canonicalize, while an underresolved zero remains
   nonzero/loud under the existing uncertainty gate.  Do not Max[1,scale]
   generally: a resolved 10^-200 coefficient is rank data at its own scale. *)
mwDecisionScale[c_List] := Module[{scale = mwScale[c]},
  If[scale === 0, 1, scale]];

(* EpsSeries norm mirror: numerics pass through; symbolics Together *)
mwNorm[c_List] := If[FreeQ[c, _Symbol], c,
  Map[If[NumericQ[#], #, Together[Expand[#]]] &, c]];

(* Matching-rank trim: advance min past coefficients certified zero at the
   matching pivot tolerance RankTol, not at LaurentLeadTol.  The latter is
   the structural/source-leading floor; using its fixed 10^-24 threshold in
   every Gaussian row operation discarded resolved rank information (at
   WP1000 RankTol is 10^-250) and left coherent endpoint residuals.  The
   binding Transport contract assigns matching pivot/lead decisions to the
   library-wide ternary NumericallyZeroQ RankTol gate.  The top never
   changes; all-zero gives canonical {top,{0}}. *)
mwTrim[{min_, c_}, context_String:"matching Laurent elimination"] := Module[
  {scale = mwDecisionScale[c], n = Length[c], i = 1,
   rtol = DiffExp2`Tolerances`Tol["RankTol"]},
  While[i <= n && TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[
      c[[i]], scale, rtol, context,
      DiffExp2`Tolerances`$AmbiguityBandDecades, False]], i++];
  If[i > n, {min + n - 1, {0}}, {min + i - 1, Drop[c, i - 1]}]];

(* Input normalization has two distinct proof domains.  Negative epsilon
   orders are compensated polar content and may be removed at the calibrated
   Laurent floor.  At and above eps^0 they are rank data, so the much tighter
   RankTol gate owns the decision (a resolved 1e-40 coefficient at WP500 is
   not structural zero).  Keeping these roles separate prevents polar noise
   from entering the Laurent solve without erasing a small F[0] pivot. *)
mwInputTrim[{min_, c_}, context_String] := Module[
  {scale = mwDecisionScale[c], n = Length[c], i = 1,
   ltol = DiffExp2`Tolerances`Tol["LaurentLeadTol"], entry},
  While[i <= n && min + i - 1 < 0 &&
      TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[c[[i]], scale, ltol,
        context <> ": compensated polar input",
        DiffExp2`Tolerances`$AmbiguityBandDecades,
        SameQ[ltol, 10^-24]]], i++];
  entry = If[i > n, {min + n - 1, {0}},
    {min + i - 1, Drop[c, i - 1]}];
  mwTrim[entry, context <> ": rank input"]];

(* Row arithmetic must discard only a certified centered zero.  Scaling a
   Schur series by the maximum over its entire future window is invalid when
   coefficients grow by many decades per epsilon order: it erased genuine
   early pivots and collapsed the banana endpoint from 13 coefficients to 3.
   Initial structural/rank classification is completed by mwInputTrim; after
   that, only an exact zero or an inexact centered zero whose entire
   uncertainty ball lies below the matching residual contract may advance
   the formal valuation.
   In particular, an underresolved 0``acc is not evidence of cancellation. *)
mwCenteredZeroQ[c_, context_String, scale_:1] := Module[
  {digits = DiffExp2`Tolerances`Tol["ChopDigits"], bounds, lower, upper,
   center, mtol},
  Which[
    NumericQ[c],
      If[!InexactNumberQ[c], Return[TrueQ[PossibleZeroQ[c]], Module]];
      bounds = numMagBounds[c, digits];
      {lower, upper} = bounds;
      center = numMag[c, digits];
      mtol = Max[DiffExp2`Tolerances`Tol["MatchTol"],
        DiffExp2`Tolerances`Tol["LaurentLeadTol"]];
      Which[
        (* Input trimming has already classified structural smallness.  A
           resolved nonzero Schur coefficient remains formal data no matter
           how small it is; only its uncertainty ball relative to zero is
           relevant here. *)
        !SameQ[center, 0] && TrueQ[lower > 0], False,
        SameQ[center, 0] && TrueQ[upper <= mtol*scale], True,
        True,
        err["E5", <|"Context" -> context, "Coefficient" -> c,
          "UncertaintyLowerBound" -> lower,
          "UncertaintyUpperBound" -> upper, "Scale" -> scale,
          "EffectiveTolerance" -> mtol,
          "Detail" -> "coefficient uncertainty overlaps zero in matching elimination"|>]],
    True, zeroQ[c]]];

mwCancellationTrim[{min_, c_}, context_String:"matching row cancellation",
    sourceFrames_List:{}] := Module[{n = Length[c], i = 1, scaleAt},
    scaleAt[k_] := If[sourceFrames === {}, 1,
      Max[1, Sequence @@ (numMag[#, 20] & /@
        Select[(mwCoeff[#, k] & /@ sourceFrames), NumericQ])]];
    While[i <= n && mwCenteredZeroQ[c[[i]], context,
        scaleAt[min + i - 1]], i++];
    If[i > n, {min + n - 1, {0}},
      {min + i - 1, Drop[c, i - 1]}]];

(* trimmed-entry zero test (the ESLeading === None equivalent) *)
mwZeroQ[{_, c_}] := c === {0};

mwNeg[{m_, c_}] := {m, -c};
mwScaleBy[a_, {m_, c_}] := {m, mwNorm[a*c]};

mwCoeff[{m_, c_}, k_Integer] := If[m <= k <= m + Length[c] - 1,
  c[[k - m + 1]], 0];
mwTop[{m_, c_}] := m + Length[c] - 1;

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

(* A rank decision at ONE epsilon order may use a scale from that same
   coefficient matrix.  It must never look at later epsilon coefficients:
   their growth is unrelated to the formal valuation. *)
mwRankZeroQ[c_, scale_, context_String] := Which[
  FreeQ[c, _?InexactNumberQ] && TrueQ[PossibleZeroQ[c]], True,
  NumericQ[c], TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[
    c, scale, DiffExp2`Tolerances`Tol["RankTol"], context,
    DiffExp2`Tolerances`$AmbiguityBandDecades, False]],
  True, zeroQ[c]];

(* Full-pivot elimination of a single coefficient matrix.  When deficient,
   return one certified right-null relation in ORIGINAL column order.  This
   is deliberately local rather than NullSpace[..., Tolerance -> ...]: every
   zero/nonzero decision goes through the package ambiguity gate, and tiny
   coordinates are snapped only after the same certified RankTol test. *)
mwNullRelation[g0_List, label_String] := Module[
  {n = Length[g0], original = g0, A, colScales, colPerm, rank = 0,
   pos, active, nums, scale, cands, piv, factor, y, x, target,
   residual, terms, rowScale},
  (* Constant column magnitudes are units of the epsilon-adic coefficient
     field, not determinant valuations.  Normalize them before rank
     elimination so an ordinary but ill-scaled fundamental matrix (for
     example columns differing by 10^70) is not mistaken for an epsilon
     degeneracy.  Map the null relation back afterward. *)
  colScales = Table[Max[10^-300,
      Sequence @@ (numMag[#,
        DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]] & /@
          Select[g0[[All, c]], NumericQ])],
    {c, n}];
  A = Table[g0[[r, c]]/colScales[[c]], {r, n}, {c, n}];
  If[!MatrixQ[A] || Dimensions[A] =!= {n, n},
    err["E5", <|"Chart" -> label,
      "Detail" -> "matching saturation received a nonsquare leading matrix"|>]];
  colPerm = Range[n];
  For[pos = 1, pos <= n, pos++,
    active = Flatten[Table[{r, c}, {r, pos, n}, {c, pos, n}], 1];
    nums = Select[(A[[#[[1]], #[[2]]]] & /@ active), NumericQ];
    scale = Max[1, Sequence @@ (numMag[#, 20] & /@ nums)];
    cands = Select[active, !mwRankZeroQ[
        A[[#[[1]], #[[2]]]], scale,
        label <> ": saturation rank pivot"] &];
    If[cands === {}, Break[]];
    piv = First[SortBy[cands, If[NumericQ[A[[#[[1]], #[[2]]]]],
        -numMag[A[[#[[1]], #[[2]]]], 20], 0] &]];
    If[piv[[1]] =!= pos, A[[{pos, piv[[1]]}]] = A[[{piv[[1]], pos}]]];
    If[piv[[2]] =!= pos,
      A[[All, {pos, piv[[2]]}]] = A[[All, {piv[[2]], pos}]];
      colPerm[[{pos, piv[[2]]}]] = colPerm[[{piv[[2]], pos}]]];
    Do[
      factor = A[[r, pos]]/A[[pos, pos]];
      A[[r, pos]] = 0;
      Do[A[[r, c]] = If[FreeQ[{A[[r, c]], factor, A[[pos, c]]}, _Symbol],
          A[[r, c]] - factor*A[[pos, c]],
          Together[Expand[A[[r, c]] - factor*A[[pos, c]]]]],
        {c, pos + 1, n}],
      {r, pos + 1, n}];
    rank++];
  If[rank === n, Return[<|"Rank" -> n, "Vector" -> None|>, Module]];
  y = ConstantArray[0, n];
  y[[rank + 1]] = 1;
  Do[y[[r]] = If[FreeQ[A, _Symbol],
      -Sum[A[[r, c]]*y[[c]], {c, r + 1, n}]/A[[r, r]],
      Together[-Sum[A[[r, c]]*y[[c]], {c, r + 1, n}]/A[[r, r]]]],
    {r, rank, 1, -1}];
  x = ConstantArray[0, n];
  Do[x[[colPerm[[c]]]] = y[[c]], {c, n}];
  x = MapThread[#1/#2 &, {x, colScales}];
  (* Do not threshold null-vector coordinates in isolation.  A component
     of size 10^-70 can be essential beside a column of size 10^70; only
     the certified matrix-vector residual decides whether the relation is
     valid. *)
  target = First[Ordering[Map[If[NumericQ[#], numMag[#, 20], 1] &, x], -1]];
  If[zeroQ[x[[target]]],
    err["E5", <|"Chart" -> label, "Rank" -> rank,
      "Detail" -> "matching saturation produced a zero null vector"|>]];
  x = mwNorm[x/x[[target]]];
  Do[
    terms = Table[original[[r, c]]*x[[c]], {c, n}];
    residual = If[FreeQ[terms, _Symbol], Total[terms], Together[Total[terms]]];
    rowScale = Max[1, Sequence @@ (numMag[#, 20] & /@
      Select[terms, NumericQ])];
    If[!mwRankZeroQ[residual, rowScale,
        label <> ": saturation null-vector residual"],
      err["E5", <|"Chart" -> label, "Rank" -> rank, "Row" -> r,
        "Residual" -> residual, "Scale" -> rowScale,
        "Detail" -> "candidate saturation relation is not in the certified nullspace"|>]],
    {r, n}];
  <|"Rank" -> rank, "Vector" -> x, "Target" -> target|>];

(* Certify that one column combination has no eps^0 value coefficient, then
   divide its honest finite window by epsilon.  The scale is formed only
   from terms at eps^0; later large coefficients cannot change the result. *)
mwSaturationDivide[rowFrames_List, a_List, label_String] := Module[
  {active, terms, combo, c0, term0, scale},
  active = Select[Range[Length[a]], !zeroQ[a[[#]]] &];
  terms = mwScaleBy[a[[#]], rowFrames[[#]]] & /@ active;
  If[terms === {},
    err["E5", <|"Chart" -> label,
      "Detail" -> "empty saturation column combination"|>]];
  combo = Fold[mwAdd, First[terms], Rest[terms]];
  If[combo[[1]] < 0,
    err["E5", <|"Chart" -> label, "Window" -> {combo[[1]], mwTop[combo]},
      "Detail" -> "saturation combination retained a material epsilon pole"|>]];
  If[combo[[1]] === 0,
    c0 = First[combo[[2]]];
    term0 = Table[a[[j]]*mwCoeff[rowFrames[[j]], 0], {j, active}];
    scale = Max[1, Sequence @@ (numMag[#, 20] & /@
      Select[term0, NumericQ])];
    If[!mwRankZeroQ[c0, scale, label <> ": saturation eps^0 divisibility"],
      err["E5", <|"Chart" -> label, "Coefficient" -> c0, "Scale" -> scale,
        "Detail" -> "null-column combination is not divisible by epsilon"|>]];
    If[Length[combo[[2]]] <= 1,
      err["E5", <|"Chart" -> label,
        "Detail" -> "insufficient complete epsilon window for saturation division"|>]];
    combo = {1, Rest[combo[[2]]]}];
  combo = mwCancellationTrim[combo, label <> ": saturation quotient"];
  {combo[[1]] - 1, combo[[2]]}];

(* Determinant of a nonnegative epsilon-frame matrix in polynomial time.
   Newton's identities compute det from Tr[M^k] using only ring addition,
   multiplication, and division by the exact integers 1..n.  This replaces
   the former Leibniz sum over n! permutations (which is already impossible
   for the 16-master kite endpoint).

   After column-valuation normalization every entry starts at eps^0 or
   above.  If L is the shortest stored nonzero frame length, coefficients
   eps^0..eps^(L-1) of EVERY determinant term are complete: each product's
   honest Cauchy length is at least L and its valuation is nonnegative.
   Working on that common rectangle is conservative and keeps the same
   finite-window contract without enumerating perfect matchings. *)
mwDetFrame[M_List, label_String] := Module[
  {n = Length[M], nonzero, len, zero, one, A, power, traces = {}, elementary,
   polyMul, polyAdd, matMul, trace, p, s, det},
  nonzero = Select[Flatten[M, 1], !mwZeroQ[#] &];
  If[nonzero === {},
    err["E5", <|"Chart" -> label,
      "Detail" -> "matching determinant is identically zero in the complete window"|>]];
  If[AnyTrue[nonzero, First[#] < 0 &],
    err["E5", <|"Chart" -> label,
      "Detail" -> "determinant frame received a negative epsilon valuation after normalization"|>]];
  len = Min[Length[#[[2]]] & /@ nonzero];
  zero = ConstantArray[0, len];
  one = ReplacePart[zero, 1 -> 1];
  polyAdd[a_List, b_List] := mwNorm[a + b];
  polyMul[a_List, b_List] := mwNorm[
    Take[ListConvolve[a, b, {1, -1}, 0], len]];
  matMul[X_List, Y_List] := Table[Fold[polyAdd, zero,
      Table[polyMul[X[[r, j]], Y[[j, c]]], {j, n}]],
    {r, n}, {c, n}];
  trace[X_List] := Fold[polyAdd, zero, Table[X[[i, i]], {i, n}]];
  A = Map[Table[mwCoeff[#, k], {k, 0, len - 1}] &, M, {2}];
  power = A;
  elementary = {one};
  Do[
    p = trace[power];
    AppendTo[traces, p];
    s = Fold[polyAdd, zero, Table[
      mwNorm[(-1)^(i - 1)*polyMul[elementary[[k - i + 1]], traces[[i]]]],
      {i, 1, k}]];
    AppendTo[elementary, mwNorm[s/k]];
    If[k < n, power = matMul[power, A]],
    {k, 1, n}];
  det = Last[elementary];
  mwCancellationTrim[{0, det}, label <> ": determinant cancellation"]];

(* Construct a sequence of elementary epsilon-adic column operations which
   turns the match-point value matrix into a regular full-rank frame.
   Initial shifts normalize each column valuation.  Every later action
   replaces column j by (Sum_i a_i column_i)/eps, lowering the determinant
   valuation by one while consuming exactly one honest top coefficient. *)
mwSaturationPlan[Fmat_List, label_String] := Module[
  {n = Length[Fmat], G, nonzero, q, shifts, actions = {}, g0, rel, a,
   target, detFrame, detValuation = None, maxSteps, steps = 0},
  G = Map[mwInputTrim[mwFromES[#], label <> ": saturation input"] &,
    Fmat, {2}];
  q = Table[
    nonzero = Select[G[[All, j]], !mwZeroQ[#] &];
    If[nonzero === {},
      err["E5", <|"Chart" -> label, "Column" -> j,
        "Detail" -> "zero matching column before epsilon saturation"|>]];
    Min[First /@ nonzero],
    {j, n}];
  shifts = -q;
  Do[G[[All, j]] = Map[{#[[1]] + shifts[[j]], #[[2]]} &, G[[All, j]]],
    {j, n}];
  maxSteps = n*(1 + Max[0, Min[mwTop /@ Flatten[G]]]);
  While[True,
    g0 = Table[mwCoeff[G[[r, c]], 0], {r, n}, {c, n}];
    rel = Catch[mwNullRelation[g0,
      label <> "#sat" <> ToString[steps]], "DiffExp2Error"];
    If[!FailureQ[rel] && rel["Rank"] === n,
      If[IntegerQ[detValuation] && steps =!= detValuation,
        err["E5", <|"Chart" -> label, "Steps" -> steps,
          "DeterminantValuation" -> detValuation,
          "Detail" -> "saturation reached full rank before consuming the certified determinant valuation"|>]];
      Break[]];
    If[detValuation === None,
      detFrame = mwInputTrim[mwDetFrame[G, label],
        label <> ": determinant valuation"];
      If[mwZeroQ[detFrame],
        err["E5", <|"Chart" -> label,
          "Detail" -> "matching determinant is unresolved in the complete epsilon window"|>]];
      detValuation = detFrame[[1]];
      If[detValuation < 0,
        err["E5", <|"Chart" -> label,
          "DeterminantValuation" -> detValuation,
          "Detail" -> "normalized matching determinant retained an epsilon pole"|>]];
      (* Valuation zero means a constant but possibly ill-conditioned unit.
         It must not trigger epsilon division. *)
      If[detValuation === 0, Break[]]];
    If[FailureQ[rel], Throw[rel, "DiffExp2Error"]];
    If[steps >= detValuation,
      err["E5", <|"Chart" -> label, "Steps" -> steps,
        "DeterminantValuation" -> detValuation,
        "Detail" -> "certified saturation steps did not produce a full-rank leading matrix"|>]];
    a = rel["Vector"];
    target = rel["Target"];
    G[[All, target]] = Table[mwSaturationDivide[G[[r]], a,
      label <> "#sat" <> ToString[steps + 1] <> "/row" <> ToString[r]],
      {r, n}];
    AppendTo[actions, <|"Target" -> target, "Vector" -> a|>];
    steps++;
    If[steps > maxSteps,
      err["E5", <|"Chart" -> label, "Steps" -> steps,
        "Detail" -> "epsilon saturation did not reach a full-rank leading matrix"|>]]];
  <|"Matrix" -> Map[mwToES, G, {2}], "InitialShifts" -> shifts,
    "Actions" -> actions, "Steps" -> steps, "Label" -> label|>];

mwShiftLocalSolutionEps[ls_Association, shift_Integer] := Join[ls, <|
  "EpsWindow" -> <|"Min" -> ls["EpsWindow", "Min"] + shift,
    "CompleteMax" -> ls["EpsWindow", "CompleteMax"] + shift|>|>];

mwApplySaturationPlan[basis_List, plan_Association] := Module[
  {out = MapThread[mwShiftLocalSolutionEps, {basis, plan["InitialShifts"]}],
   a, target, active, combo},
  Do[
    a = action["Vector"]; target = action["Target"];
    active = Select[Range[Length[a]], !zeroQ[a[[#]]] &];
    combo = DiffExp2`SectorSeries`CombineLocalSolutions[a[[active]], out[[active]]];
    out[[target]] = mwShiftLocalSolutionEps[combo, -1],
    {action, plan["Actions"]}];
  out];

(* Once epsilon-lattice saturation has made the leading value matrix a unit,
   normalize the solution frame at THIS match point by the epsilon-independent
   right action P = F_eps0^-1.  This is ordinary fundamental-matrix frame
   normalization: whole LocalSolutions (including all fractional/log sectors
   and their prescriptions) are mixed together, so no tag, branch, valuation,
   or epsilon-completeness rule changes.  Keeping P local to the march avoids
   making a sheet-specific frame part of the global chart cache. *)
mwLeadingValueMatrix[Fmat_List, label_String] := Module[{n = Length[Fmat]},
  If[n == 0 || !MatrixQ[Fmat] || Dimensions[Fmat] =!= {n, n},
    err["E5", <|"Chart" -> label,
      "Detail" -> "match preconditioner received a nonsquare value matrix"|>]];
  If[AnyTrue[Flatten[Fmat], esMin[#] < 0 &],
    err["E5", <|"Chart" -> label,
      "Detail" -> "match preconditioner received a value matrix with epsilon poles after saturation"|>]];
  Table[If[esMin[Fmat[[r, c]]] <= 0 <= esCM[Fmat[[r, c]]],
      esCoeff[Fmat[[r, c]], 0], 0], {r, n}, {c, n}]];

mwContractZeroQ[res_, contributions_List, label_String] := Module[
  {tol = Max[DiffExp2`Tolerances`Tol["MatchTol"],
      DiffExp2`Tolerances`Tol["LaurentLeadTol"]], scale},
  scale = Max[1, Sequence @@ (numMag[#, 20] & /@
      Select[Flatten[contributions], NumericQ])];
  Which[
    FreeQ[res, _?InexactNumberQ] && TrueQ[PossibleZeroQ[res]], True,
    NumericQ[res], TrueQ[DiffExp2`Tolerances`NumericallyZeroQ[
      res, scale, tol, label, DiffExp2`Tolerances`$AmbiguityBandDecades,
      True]],
    True, zeroQ[res]]];

mwIdentityMatrixQ[M_List, label_String] := Module[{n = Length[M]},
  MatrixQ[M] && Dimensions[M] === {n, n} &&
    And @@ Flatten[Table[Module[{target = KroneckerDelta[r, c], res},
      res = M[[r, c]] - target;
      If[!FreeQ[res, _Symbol], res = Together[Expand[res]]];
      mwContractZeroQ[res, {M[[r, c]], target},
        label <> "/" <> ToString[r] <> "," <> ToString[c]]],
      {r, n}, {c, n}]]];

mwProductIdentityAssert[left_List, right_List, label_String] := Module[
  {n = Length[left]},
  Do[Module[{terms = Table[left[[r, j]]*right[[j, c]], {j, n}],
      target = KroneckerDelta[r, c], res},
    res = Total[terms] - target;
    If[!FreeQ[res, _Symbol], res = Together[Expand[res]]];
    If[!mwContractZeroQ[res, Join[terms, {target}],
        label <> "/" <> ToString[r] <> "," <> ToString[c]],
      err["E5", <|"Chart" -> label, "Row" -> r, "Column" -> c,
        "Residual" -> res,
        "Detail" -> "F_eps0.P failed the identity contract"|>]]],
    {r, n}, {c, n}];
  Null];

mwConstantMatchPrecondition[basis_List, Fmat_List, basisValues_,
    requiredTop_Integer, label_String] := Module[
  {fEps0 = mwLeadingValueMatrix[Fmat, label], n = Length[Fmat], P, out,
   Fprime, fPrimeEps0, verify, sharedTop, sharedTopPrime},
  (* An already normalized frame (notably the anchor chart) needs no dense
     recombination.  The uncertainty-aware contract makes this skip safe for
     both exact and arbitrary-precision numeric matrices. *)
  If[mwIdentityMatrixQ[fEps0, label <> "#already-normalized"],
    Return[<|"Basis" -> basis, "Matrix" -> Fmat,
      "LeadingInverse" -> IdentityMatrix[n], "Applied" -> False|>, Module]];
  P = Quiet[Check[LinearSolve[fEps0, IdentityMatrix[n]], $Failed]];
  If[P === $Failed || !MatrixQ[P] || Dimensions[P] =!= {n, n} ||
      !FreeQ[P, Indeterminate | ComplexInfinity | DirectedInfinity[_]],
    err["E5", <|"Chart" -> label,
      "Detail" -> "could not construct a finite constant match preconditioner"|>]];
  (* Do not threshold small entries of P: a tiny coefficient can be essential
     next to a large column.  Certify the full contribution-aware product. *)
  mwProductIdentityAssert[fEps0, P, label <> "#raw-preconditioner"];
  sharedTop = Min[# ["EpsWindow", "CompleteMax"] & /@ basis];
  out = Table[DiffExp2`SectorSeries`CombineLocalSolutions[P[[All, j]], basis],
    {j, n}];
  sharedTopPrime = Min[# ["EpsWindow", "CompleteMax"] & /@ out];
  If[sharedTopPrime < sharedTop || sharedTopPrime < requiredTop,
    err["E5", <|"Chart" -> label, "RequiredTop" -> requiredTop,
      "SharedCompleteMaxBefore" -> sharedTop,
      "SharedCompleteMaxAfter" -> sharedTopPrime,
      "Detail" -> "constant match preconditioning reduced the honest epsilon window"|>]];
  (* Re-evaluation of the transformed LocalSolutions is authoritative: it
     checks the actual tagged objects, not just matrix algebra on old values. *)
  Fprime = basisValues[out];
  fPrimeEps0 = mwLeadingValueMatrix[Fprime, label <> "#transformed"];
  If[!mwIdentityMatrixQ[fPrimeEps0, label <> "#transformed-identity"],
    err["E5", <|"Chart" -> label,
      "Detail" -> "re-evaluated preconditioned basis is not identity at epsilon order zero"|>]];
  verify = mwSaturationPlan[Fprime, label <> "#precondition-verify"];
  If[verify["Steps"] =!= 0 ||
      AnyTrue[verify["InitialShifts"], # =!= 0 &],
    err["E5", <|"Chart" -> label,
      "VerificationShifts" -> verify["InitialShifts"],
      "VerificationSteps" -> verify["Steps"],
      "Detail" -> "constant match preconditioner did not preserve the regular epsilon lattice"|>]];
  <|"Basis" -> out, "Matrix" -> Fprime, "LeadingInverse" -> P,
    "Applied" -> True|>];

(* Residual vector for iterative refinement.  It is formed against the
   ORIGINAL, untrimmed matching matrix: the primary elimination may classify
   individually sub-LaurentLeadTol leading coefficients as structural zero,
   but several such discarded terms can add coherently above the tolerance.
   Solving that small residual on its own rescales it to O(1), so the same
   relative structural gate retains the correction without changing either
   the pivot convention or the public tolerance. *)
matchingResidualVector[Fmat_List, vIn_List, weights_List] := Module[
  {nb = Length[Fmat]},
  Table[Module[{terms, lhs},
    terms = Table[esTimes[Fmat[[comp, col]], weights[[col]]], {col, nb}];
    lhs = Fold[esAdd, First[terms], Rest[terms]];
    esAdd[vIn[[comp]], esScale[-1, lhs]]],
    {comp, nb}]];

mwDebugSummary[x_] := Module[{nums, accs, mags},
  nums = Cases[x, z_?InexactNumberQ :> z, Infinity];
  If[nums === {}, Return[<|"Count" -> 0|>, Module]];
  accs = Accuracy /@ nums;
  mags = numMag[#, 20] & /@ nums;
  <|"Count" -> Length[nums], "AccuracyMinMax" -> MinMax[accs],
    "MagnitudeMinMax" -> MinMax[mags]|>];

(* Non-throwing form of the residual proof.  MatchWeights uses it to decide
   whether a bounded refinement is needed; matchingResidualAssert remains the
   single loud public proof obligation after the final correction. *)
matchingResidualFailure[Fmat_List, vIn_List, weights_List,
    label_String] := Module[
  {nb = Length[Fmat], mtol = DiffExp2`Tolerances`Tol["MatchTol"],
   ltol = DiffExp2`Tolerances`Tol["LaurentLeadTol"], effectiveTol},
  effectiveTol = Max[mtol, ltol];
  Catch[
    Do[Module[{terms, lhs, rhs = vIn[[comp]], kmin, kmax},
      terms = Table[esTimes[Fmat[[comp, col]], weights[[col]]], {col, nb}];
      lhs = Fold[esAdd, First[terms], Rest[terms]];
      kmin = Min[esMin[lhs], esMin[rhs]];
      kmax = Min[esCM[lhs], esCM[rhs]];
      Do[Module[{lv = esCoeff[lhs, k], rv = esCoeff[rhs, k], residual,
          termVals, nums, scale, mag, uncertainty, bad},
        residual = If[FreeQ[{lv, rv}, _Symbol], lv - rv,
          Together[Expand[lv - rv]]];
        termVals = esCoeff[#, k] & /@ terms;
        nums = Select[Flatten[{lv, rv, termVals}], NumericQ];
        scale = Max[1,
          Sequence @@ (numMag[#, 20] & /@ nums)];
        uncertainty = If[InexactNumberQ[residual],
          Module[{acc = Accuracy[residual]},
            (* Accuracy is commonly returned as a machine real even for a
               many-hundred-digit number.  Evaluating 10^-acc would therefore
               underflow at about 308 digits and silently turn the uncertainty
               allowance into zero.  The exact decade below is conservative
               (at most a factor ten larger) and cannot underflow. *)
            If[NumericQ[acc] && acc =!= Infinity,
              10^-Floor[acc], 0]], 0];
        bad = Which[
          NumericQ[residual],
            mag = numMag[residual, 20];
            TrueQ[mag + uncertainty > effectiveTol*scale],
          TrueQ[PossibleZeroQ[residual]], False,
          (* Exact symbolic nonzero content is a proof.  An inexact symbolic
             residue has no parameter domain or coefficient-wise significance
             certificate here, so accepting it would be an unproved match. *)
          FreeQ[residual, _?InexactNumberQ], !zeroQ[residual],
          True, True];
        If[TrueQ[bad],
          Throw[<|"Chart" -> label, "EpsOrder" -> k,
            "Component" -> comp, "Residual" -> residual,
            "ResidualUncertainty" -> uncertainty, "Scale" -> scale,
            "MatchTol" -> mtol, "LaurentLeadTol" -> ltol,
            "EffectiveTolerance" -> effectiveTol,
            "Detail" -> "matching residual F.w-v exceeds tolerance"|>,
            "MatchingResidualFailure"]]],
        {k, kmin, kmax}]],
      {comp, nb}];
    None,
    "MatchingResidualFailure"]];

(* Always-on proof obligation for the Laurent matching solve.  Compare the
   reconstructed F.w with v on the full shared complete window, including
   certified-zero orders below either operand's Min.  Gaussian elimination
   deliberately uses the LaurentLeadTol structural-zero gate while trimming
   rows; therefore the strongest residual contract for the ORIGINAL,
   untrimmed F is Max[MatchTol, LaurentLeadTol].  Demanding MatchTol below
   that floor is internally inconsistent: content the solve was required to
   discard would immediately fail its checker.  Inexact residuals include a
   conservative 10^-Accuracy uncertainty
   allowance; low precision never turns a resolved violation into a pass.
   Exact/symbolic nonzero residuals remain algebraically checked. *)
matchingResidualAssert[Fmat_List, vIn_List, weights_List, label_String,
    checkedFailure_:Automatic] := Module[
  (* MatchWeights passes the result of the immediately preceding check so a
     successful seam does not reconstruct F.w twice.  Direct callers retain
     the four-argument always-compute contract. *)
  {failure = If[checkedFailure === Automatic,
      matchingResidualFailure[Fmat, vIn, weights, label], checkedFailure]},
  If[AssociationQ[failure], err["E6", failure]];
  Null];

(* One Laurent-field elimination.  It intentionally has no residual policy:
   MatchWeights owns the original-system proof and may call this kernel on a
   residual right-hand side during refinement. *)
mwSolve[Fmat_List, vIn_List, label_String] := Module[
  {nb = Length[Fmat], FF, vv, scaleFrames, rowPerm, colPerm, wPerm, w, wES},
  (* Gaussian elimination over the eps-Laurent field on frame lists:
     FULL row+column pivoting by minimal leading eps-order (then largest
     leading magnitude), followed by quotient-recursion row operations.
     Column pivoting is essential for compensated Laurent bases: a fixed
     column order can manufacture a near-zero final Schur pivot and lose all
     significance even though a well-conditioned pivot exists elsewhere in
     the remaining submatrix.  Unknowns are mapped back to their original
     basis-column order after back substitution.  Honest windows propagate
     through the mirrored EpsSeries rules; no usable submatrix pivot means a
     genuinely singular system (E5). *)
  FF = Map[mwInputTrim[mwFromES[#], label <> ": input"] &, Fmat, {2}];
  vv = mwFromES /@ vIn;
  (* Same-order input scale retained throughout elimination.  A zero can
     have modest absolute Accuracy after canceling 10^N-sized operands yet
     still be certified to thousands of relative digits.  Keeping the
     original coefficient frames lets cancellation trimming recover that
     scale without ever consulting later epsilon orders. *)
  scaleFrames = Join[Flatten[FF], vv];
  rowPerm = Range[nb];
  colPerm = Range[nb];
  Do[Module[{cands, piv},
    cands = Select[Tuples[{Range[col, nb], Range[col, nb]}],
      !mwZeroQ[FF[[rowPerm[[#[[1]]]], colPerm[[#[[2]]]]]]] &];
    If[cands === {},
      err["E5", <|"Chart" -> label, "Column" -> col,
        "Detail" -> "matching system singular over the eps-Laurent field"|>]];
    piv = First[SortBy[cands, Module[
      {e = FF[[rowPerm[[#[[1]]]], colPerm[[#[[2]]]]]]},
      {e[[1]], -numMag[e[[2, 1]], 20]}] &]];
    If[piv[[1]] =!= col,
      rowPerm[[{col, piv[[1]]}]] = rowPerm[[{piv[[1]], col}]]];
    If[piv[[2]] =!= col,
      colPerm[[{col, piv[[2]]}]] = colPerm[[{piv[[2]], col}]]];
    Do[Module[{entry = FF[[rowPerm[[r]], colPerm[[col]]]], factor},
      If[!mwZeroQ[entry],
        factor = mwCancellationTrim[mwDiv[entry,
          FF[[rowPerm[[col]], colPerm[[col]]]]],
          label <> ": elimination factor", scaleFrames];
        Do[Module[{old = FF[[rowPerm[[r]], colPerm[[c2]]]], term},
          term = mwMul[factor, FF[[rowPerm[[col]], colPerm[[c2]]]]];
          term = mwNeg[term];
          FF[[rowPerm[[r]], colPerm[[c2]]]] = mwCancellationTrim[
            mwAdd[old, term], label <> ": row elimination",
            Join[{old, term}, scaleFrames]]],
          {c2, col, nb}];
        Module[{old = vv[[rowPerm[[r]]]], term},
          term = mwMul[factor, vv[[rowPerm[[col]]]]];
          term = mwNeg[term];
          vv[[rowPerm[[r]]]] = mwCancellationTrim[mwAdd[old, term],
            label <> ": rhs elimination",
            Join[{old, term}, scaleFrames]]]]],
      {r, col + 1, nb}]],
    {col, nb}];
  (* back substitution *)
  wPerm = Table[None, {nb}];
  Do[Module[{rhs = vv[[rowPerm[[col]]]]},
    Do[Module[{old = rhs, term},
      term = mwMul[FF[[rowPerm[[col]], colPerm[[c2]]]], wPerm[[c2]]];
      term = mwNeg[term];
      rhs = mwCancellationTrim[mwAdd[old, term],
        label <> ": back substitution", Join[{old, term}, scaleFrames]]],
      {c2, col + 1, nb}];
    wPerm[[col]] = mwCancellationTrim[mwDiv[rhs,
      FF[[rowPerm[[col]], colPerm[[col]]]]],
      label <> ": pivot division", scaleFrames]],
    {col, nb, 1, -1}];
  w = Table[None, {nb}];
  Do[w[[colPerm[[col]]]] = wPerm[[col]], {col, nb}];
  wES = mwToES /@ w;
  wES];

MatchWeights[Fmat_List, vIn_List, label_String] := Module[
  {w = mwSolve[Fmat, vIn, label], failure, residual, correction,
   refinementSteps = 0, maxRefinementSteps = 2},
  failure = matchingResidualFailure[Fmat, vIn, w, label];
  If[Environment["DEBUG_MATCHACC"] === "1",
    Print["MATCHDBG ", label, " inputF=", mwDebugSummary[Fmat],
      " inputV=", mwDebugSummary[vIn], " initialW=", mwDebugSummary[w],
      " initialFailure=", failure]];
  While[AssociationQ[failure] && refinementSteps < maxRefinementSteps,
    residual = matchingResidualVector[Fmat, vIn, w];
    If[Environment["DEBUG_MATCHACC"] === "1",
      Print["MATCHDBG ", label, " refine=", refinementSteps + 1,
        " residual=", mwDebugSummary[residual]]];
    correction = mwSolve[Fmat, residual,
      label <> "#refine" <> ToString[refinementSteps + 1]];
    (* ESAdd takes the intersection of complete upper windows.  Refinement
       must never manufacture a pass by shrinking the original weight
       window past a failing order; if the correction cannot cover every
       current complete order, leave w unchanged and let the final proof
       raise the original failure. *)
    If[AnyTrue[MapThread[esCM[#1] < esCM[#2] &,
        {correction, w}], TrueQ], Break[]];
    w = MapThread[esAdd, {w, correction}];
    refinementSteps++;
    failure = matchingResidualFailure[Fmat, vIn, w, label]];
  matchingResidualAssert[Fmat, vIn, w, label, failure];
  w];

(* ---- 2.10 probe ---- *)

SegmentErrorProbe[ls_Association, tOut_, couplingDepth_Integer] := Module[
  {dec = DiffExp2`Tolerances`EvalErrorSeriesDecrease[Max[couplingDepth, 1]],
   full, red},
  full = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut, "UsePade" -> False,
    "ImSign" -> sigmaFor[ls], "ComputeTailEstimates" -> False];
  red = DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tOut, "ImSign" -> sigmaFor[ls],
    "UsePade" -> False, "TOrderReduction" -> dec,
    "ComputeTailEstimates" -> False];
  Table[Module[{kf = esCoeff[full["Value"], k],
      kr = If[esMin[red["Value"]] <= k <= esCM[red["Value"]],
        esCoeff[red["Value"], k], 0*esCoeff[full["Value"], k]]},
    Max[0, Sequence @@ (Last[numMagBounds[#, 20]] & /@
      Select[Flatten[{kf - kr}], NumericQ])]],
    {k, esMin[full["Value"]], esCM[full["Value"]]}]];

(* ---- 2.6 marching ---- *)

(* ValidatePlan: static pre-pass over the chart chain (TransportLine runs
   it before any solve).  Every chart's incoming match point — computed by
   the SAME chartMatchPoint formula the marching loop evaluates — must lie
   strictly inside both adjacent charts' physical disks and no farther than
   one half-radius out (the SegmentErrorProbe certification envelope)
   (EvaluateLocalSolution's validity bound is |t| < Radius after affine
   conversion); a singular chart's point must in addition sit on the
   approach side; the
   first chart's point is the boundary anchor plan["From"], inside the
   first disk.  A violation is a planner bug: E8 with the full chain
   geometry, raised before the expensive solves instead of deep inside
   them. *)
ValidatePlan[plan_Association] := Module[
  {charts = plan["Charts"], dir = plan["Direction"],
   k = cfg["DivisionOrder"], chain},
  If[KeyExistsQ[plan, "SegmentCount"] &&
      plan["SegmentCount"] =!= Length[charts],
    err["E8", <|"SegmentCount" -> plan["SegmentCount"],
      "ChartCount" -> Length[charts],
      "Detail" -> "plan SegmentCount does not equal the stored chart count"|>]];
  chain = Map[<|"Center" -> #["Center"], "Radius" -> #["Radius"],
    "Singular" -> #["Singular"]|> &, charts];
  Do[Module[{chart = charts[[ci]], prev, mp, expectedSingular, bad},
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
    If[ci > 1 && TrueQ[chart["Singular"]],
      expectedSingular = singularMatchPoint[prev["Center"], prev["Radius"],
        chart["Center"], chart["Radius"], dir];
      If[expectedSingular === None ||
          !TrueQ[PossibleZeroQ[Together[mp - expectedSingular]]],
        err["E8", <|"Chart" -> chart["Name"], "MatchPoint" -> mp,
          "ExpectedMatchPoint" -> expectedSingular,
          "PrevCenter" -> prev["Center"], "PrevRadius" -> prev["Radius"],
          "Center" -> chart["Center"], "Radius" -> chart["Radius"],
          "Chain" -> chain,
          "Detail" -> "stored singular handoff differs from the canonical balanced approach point"|>]]];
    bad = Which[
      ci === 1,
      !TrueQ[numericDistance[mp, chart["Center"], 40] < N[chart["Radius"], 40]] ||
      !TrueQ[numericDistance[mp, chart["Center"], 40] <= N[chart["Radius"], 40]/2],
      True,
      !TrueQ[numericDistance[mp, prev["Center"], 40] < N[prev["Radius"], 40]] ||
      !TrueQ[numericDistance[mp, chart["Center"], 40] < N[chart["Radius"], 40]] ||
      (* SegmentErrorProbe certifies every handed-off LocalSolution at one
         half of its true physical convergence radius.  The actual shared
         point must be no farther out in EITHER adjacent disk; otherwise the
         accumulated error estimate would be nonconservative even though
         analytic convergence alone still holds. *)
      !TrueQ[numericDistance[mp, prev["Center"], 40] <= N[prev["Radius"], 40]/2] ||
      !TrueQ[numericDistance[mp, chart["Center"], 40] <= N[chart["Radius"], 40]/2] ||
      (TrueQ[chart["Singular"]] &&
        !TrueQ[dir*N[chart["Center"] - mp, 40] > 0])];
    If[bad,
      err["E8", <|"Chart" -> chart["Name"], "MatchPoint" -> mp,
        "PrevCenter" -> If[prev === None, None, prev["Center"]],
        "PrevRadius" -> If[prev === None, None, prev["Radius"]],
        "Center" -> chart["Center"], "Radius" -> chart["Radius"],
        "Chain" -> chain,
        "Detail" -> "incoming match point outside a producing/receiving disk, beyond the half-radius error-probe envelope, or on the wrong singular approach side (planner bug)"|>]]],
    {ci, Length[charts]}];
  If[!TrueQ[Lookup[plan, "EndpointIsSingular", False]],
    Module[{last = Last[charts], endpoint = plan["To"]},
      If[!TrueQ[numericDistance[endpoint, last["Center"], 40] <
            N[last["Radius"], 40]] ||
          !TrueQ[numericDistance[endpoint, last["Center"], 40] <=
            N[last["Radius"], 40]/2],
        err["E8", <|"Chart" -> last["Name"], "Endpoint" -> endpoint,
          "Center" -> last["Center"], "Radius" -> last["Radius"],
          "Chain" -> chain,
          "Detail" -> "regular final endpoint lies beyond the half-radius error-probe envelope"|>]]]];
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
    Do[Module[{chart = charts[[ci]], cs, sol, matchPt, tIn, centerTIn,
      matchTIn, vvals, F, w, ls, basis, satPlan, satVerify,
      probeErrs, valueMode,
      couplingDepth = 0},
    If[Environment["DEBUG_CHART"] === "1",
      Print["CHART ", chart["Name"], " prep start t=", SessionTime[]]];
    cs = DiffExp2`Solve`PrepareChart[sys, chart];
    (* PROTOTYPE (env-gated, default off): value-vector propagation for
       REGULAR interior charts (Docs/PerfGapAnalysis.md lever 1).  The
       incoming object is evaluated AT THIS CHART'S CENTER — inside the
       previous chart's disk (the classic +1/k,-1/k construction keeps
       ordinary handoffs safely interior; EvaluateLocalSolution's radius assert is the loud
       backstop) — and that value is the t^0 Cauchy datum of ONE
       d-dimensional recursion (Solve`SolveValueRegular), replacing the
       d-column basis + MatchWeights + CombineLocalSolutions.  Singular
       charts keep the basis+matching path unchanged.  The first regular
       chart may use the value path when the incoming object is manifestly an
       evaluable Cauchy datum at that same center. *)
    valueMode = Environment["DE2_VALUE_TRANSPORT"] === "1" &&
      !TrueQ[chart["Singular"]] &&
      TrueQ[Lookup[cs["IndicialData"], "Regular", False]] &&
      Module[{margin = valueCenterMargin[cfg["ExpansionOrder"]],
          currentScale = current["ChartMap", "Scale"], sameCenter},
        sameCenter = TrueQ[PossibleZeroQ[RootReduce[Together[
          chart["Center"] - current["Center"]]]]];
        If[ci === 1,
          (* A plain boundary is already the exact Cauchy value at the anchor;
             no fundamental matrix or matching solve is needed there. *)
          sameCenter && centerValueDatumQ[current],
          (* conservative geometry pre-check: the center must sit WELL inside
             the previous object's disk and ratio^(ExpansionOrder+1) must stay
             two decades below LaurentLeadTol.  Value-aware SegmentLine plans
             refine regular hops to satisfy this strictly; retain the runtime
             check as a loud planner/foreign-plan backstop. *)
          TrueQ[Abs[N[(chart["Center"] - current["Center"])/currentScale, 30]] <
            margin*N[current["Radius"], 30]]]];
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
       match point (basis path) or at this chart's center (value mode).
       Across a singular chart both candidate points must be on the same
       sheet; otherwise the conservative basis path owns the handoff. *)
    centerTIn = Together[(chart["Center"] - current["Center"])/
      current["ChartMap", "Scale"]];
    matchTIn = Together[(matchPt - current["Center"])/
      current["ChartMap", "Scale"]];
    If[valueMode && lastSingular &&
        !SameQ[TrueQ[N[centerTIn, 30] dir > 0],
          TrueQ[N[matchTIn, 30] dir > 0]],
      valueMode = False];
    tIn = If[valueMode, centerTIn, matchTIn];
    Module[{sigma, crossed = False, valuesAt},
      (* A far-side transition always owns the prescription/material-content
         proof, in both directions.  ApplyCrossing itself is specifically the
         t<0 -> u=-t representation change: rightward far-side t is already
         positive and must not receive a second phase; leftward t is negative
         and is reflected exactly once. *)
      If[lastSingular && TrueQ[N[tIn, 30] dir > 0],
        sigma = DiffExp2`SectorSeries`ChartImSign[current];
        If[!MemberQ[{1, -1}, sigma],
          (* magnitude-aware gate: a multivalued sector blocks the crossing
             only if it carries MATERIAL content (box L2's apparent chart at
             7/11 has a t^(-1+eps) sector of physically zero weight; the
             syntactic tag test alone would E8 every such crossing).
             Ambiguous magnitudes stay material -> still loud. *)
          Module[{secs = current["Sectors"], scale, materialQ,
              coeffZeroCertifiedQ, ltol, floorQ, zeroLimit},
            scale = Max[1*^-300, Sequence @@ (numMag[#, 20] & /@
              Select[Flatten[#["Coeffs"] & /@ secs], NumericQ])];
            ltol = DiffExp2`Tolerances`Tol["LaurentLeadTol"];
            floorQ = SameQ[ltol, 10^-24];
            zeroLimit = ltol*scale/If[floorQ, 1,
              10^DiffExp2`Tolerances`$AmbiguityBandDecades];
            coeffZeroCertifiedQ[z_] := Which[
              FreeQ[z, _?InexactNumberQ] && TrueQ[PossibleZeroQ[z]], True,
              NumericQ[z], With[{upper = Last[numMagBounds[z, 20]]},
                If[floorQ, TrueQ[upper <= zeroLimit],
                  TrueQ[upper < zeroLimit]]],
              True, False];
            materialQ[sec_] := !AllTrue[Flatten[sec["Coeffs"]],
              coeffZeroCertifiedQ];
            If[AnyTrue[Select[secs, materialQ],
                !IntegerQ[#["a"]] || !zeroQ[#["b"]] || #["p"] > 0 &],
              err["E8", <|"Chart" -> chart["Name"],
                "Detail" -> "crossing a multivalued singular chart without a derivable Im-sign (missing DeltaPrescriptions)"|>],
              sigma = 1]]];
        If[TrueQ[N[tIn, 30] < 0],
          current = ApplyCrossing[current, sigma];
          crossed = True;
          tIn = -tIn]];  (* a reflected far side evaluates at positive u *)
      valuesAt[tt_] := Module[{ev, vv, d2 = cs["SystemSize"]},
        ev = DiffExp2`SectorSeries`EvaluateLocalSolution[current, tt,
          "UsePade" -> False, "ImSign" -> sigmaFor[current],
          "ComputeTailEstimates" -> False];
        vv = ev["Value"];
        Table[esNew[esMin[vv], numHandoff[Table[esCoeff[vv, k][[c]],
          {k, esMin[vv], esCM[vv]}]]], {c, d2}]];
      vvals = valuesAt[tIn];
      (* Repeated center-to-center value handoffs can exhaust a finite
         boundary's guard digits even though every local recurrence is
         algebraically sound.  Re-evaluate at the much nearer standard
         match point and use the basis path before that happens. *)
      If[valueMode && !valueHandoffAccurateQ[vvals],
        valueMode = False;
        If[Environment["DEBUG_CHART"] === "1",
          Print["CHART value handoff significance fallback t=", SessionTime[]]];
        tIn = matchTIn;
        If[crossed, tIn = -tIn];
        vvals = valuesAt[tIn]]];
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
        sol["Basis"]["Specs"], sol["Basis"]["Diagnostics"]];
      couplingDepth = sol["CouplingDepth"];
      (* basis values at the same point, in THIS chart's coordinate *)
      Module[{tLoc = Together[(matchPt - chart["Center"])/chart["Scale"]],
          basisValues, pre},
        basisValues[bb_List] := Module[{Feval},
          Feval = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#,
            tLoc, "UsePade" -> False, "ImSign" -> sigmaFor[#],
            "ComputeTailEstimates" -> False]["Value"] &, bb];
          Table[esNew[esMin[Feval[[i]]],
            numHandoff[Table[esCoeff[Feval[[i]], k][[c]],
              {k, esMin[Feval[[i]]], esCM[Feval[[i]]]}]]],
            {c, cs["SystemSize"]}, {i, Length[bb]}]];
        F = basisValues[basis];
        satPlan = mwSaturationPlan[F, chart["Name"]];
        If[satPlan["Steps"] > 0 ||
            AnyTrue[satPlan["InitialShifts"], # =!= 0 &],
          basis = mwApplySaturationPlan[basis, satPlan];
          (* Re-evaluate the ACTUAL transformed LocalSolutions.  The frame
             algebra constructs the certified column operations, but this
             second evaluation is the authoritative check that applying
             them to the sector objects produced the same regular lattice. *)
          F = basisValues[basis];
          satVerify = mwSaturationPlan[F, chart["Name"] <> "#verify"];
          If[satVerify["Steps"] =!= 0 ||
              AnyTrue[satVerify["InitialShifts"], # =!= 0 &],
            err["E5", <|"Chart" -> chart["Name"],
              "InitialShifts" -> satPlan["InitialShifts"],
              "SaturationSteps" -> satPlan["Steps"],
              "VerificationShifts" -> satVerify["InitialShifts"],
              "VerificationSteps" -> satVerify["Steps"],
              "Detail" -> "transformed LocalSolution basis did not verify as a regular epsilon lattice"|>]]];
        mwMaybeDumpMatchFixture[basis, F, vvals, tLoc,
          req["EpsWindow", "CompleteMax"], chart, matchPt];
        (* Ordinary charts use constant right-frame normalization in
           production.  The private seam extends the same certified GL(d),
           epsilon-independent coordinate change to singular/resonant
           charts.  CombineLocalSolutions intersects unequal operand tops,
           so mwConstantMatchPrecondition proves that the common honest top
           still covers the requested window, re-evaluates the actual tagged
           objects, and re-runs epsilon-lattice saturation before matching. *)
        If[mwUseConstantMatchPreconditionQ[cs],
          pre = mwConstantMatchPrecondition[basis, F, basisValues,
            req["EpsWindow", "CompleteMax"], chart["Name"]];
          basis = pre["Basis"];
          F = pre["Matrix"];
          If[Environment["DEBUG_MATCHACC"] === "1",
            Print["MATCHDBG ", chart["Name"],
              " preconditionedF=", mwDebugSummary[F]]]]];
      w = MatchWeights[F, vvals, chart["Name"]];
      If[AnyTrue[w, esMin[#] < 0 &],
        err["E5", <|"Chart" -> chart["Name"],
          "WeightWindows" -> ({esMin[#], esCM[#]} & /@ w),
          "Detail" -> "epsilon-saturated matching returned Laurent weights"|>]];
      ls = DiffExp2`SectorSeries`CombineLocalSolutions[w, basis]];
    (* Probe BOTH sides of a regular chart and the incoming side of a
       singular chart (sign of matchPt - center keeps it one-sided) just
       outside half radius.  Under classic
       GetCPL/GetCPR geometry the +/-1/k joins lie no farther out, while
       LineIntegral midpoint tile edges can reach the 1/2 class near a
       change of limiting projection.  Probing at 51/100 is therefore a
       conservative cover of the shared R/2 envelope and is nonzero on the
       anchor chart. *)
    probeErrs = Module[{rho = 51/100, sgn, raw, points, probes,
        twoSidedQ},
      sgn = Sign[N[matchPt - chart["Center"], 30]];
      If[!MemberQ[{-1, 1}, sgn], sgn = dir];
      (* Probe just OUTSIDE the planner/API half-radius envelope.  Keeping
         this numeric at the 2x-WP handoff precision avoids an algebraic
         evaluation grind without a Rationalize step that could move the
         probe inward past an accepted R/2 handoff/tile edge. *)
      raw = rho*N[ls["Radius"],
        DiffExp2`Tolerances`$InputPrecisionFactor*cfg["WorkingPrecision"]];
      (* Regular charts and INTERIOR singular charts can both be consumed on
         either side: LineIntegral owns a two-arm tile around an interior
         singularity and applies the prescription-aware pairing directly to
         this kept LocalSolution.  Only a singular endpoint is genuinely
         one-sided. *)
      twoSidedQ = TrueQ[Lookup[cs["IndicialData"], "Regular", False]] ||
        (TrueQ[chart["Singular"]] &&
          !(TrueQ[plan["EndpointIsSingular"]] && ci === Length[charts]));
      points = If[twoSidedQ,
        {-raw, raw}, {sgn*raw}];
      probes = SegmentErrorProbe[ls, #, couplingDepth] & /@ points;
      If[Length[probes] === 1, First[probes], MapThread[Max, probes]]];
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
    "SegmentCount" -> Length[kept],
    "EndpointIsSingular" -> plan["EndpointIsSingular"],
    "ErrorEstimate" -> errAcc,
    "Value" -> If[plan["EndpointIsSingular"], None,
      DiffExp2`SectorSeries`EvaluateLocalSolution[current,
        Together[(plan["To"] - current["Center"])/
          current["ChartMap", "Scale"]], "UsePade" -> False,
        "ImSign" -> sigmaFor[current],
        "ComputeTailEstimates" -> False]["Value"]]|>];

End[];
EndPackage[];
