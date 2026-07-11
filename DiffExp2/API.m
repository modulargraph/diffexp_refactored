(* DiffExp2/API.m — typed implementation entry points shared by the public
   umbrella and Feynman-trick runtime. *)

BeginPackage["DiffExp2`API`",
  {"DiffExp2`Tolerances`", "DiffExp2`Config`", "DiffExp2`EpsSeries`",
   "DiffExp2`SectorSeries`", "DiffExp2`Indicial`", "DiffExp2`Solve`",
   "DiffExp2`Transport`", "DiffExp2`Integrate`"}];

LoadSystem::usage = "LoadSystem[spec] loads an exact eps-rational system: <|\"Matrix\" -> m, \"Variable\" -> x|> (closed form), or <|\"FullMatrixFile\" -> path, \"Variable\" -> x|> (the d<var>_full.m exact export). Returns the unchanged exact Matrix with SingularFactors (the epsilon-zero planner alphabet) and SingularFactorsExact (the full epsilon-dependent denominator factors).";
TransportEndpoint::usage = "TransportEndpoint[sys, bvals, from, to, opts] transports plain boundary values from a regular anchor to `to` (singular endpoints return the LocalSolution). Options: \"ExtraSingularFactors\".";
LineIntegral::usage = "LineIntegral[sys, bvals, from, {lo, hi}, c, opts] gives Integrate[c(x, eps) . f(x), {x, lo, hi}] for the transported solution vector f: chart tiling, per-chart rational multiply, exact sector integrals. c is a coefficient VECTOR (one rational function per component).";
EndpointLimitValues::usage = "EndpointLimitValues[transportResult, cvec] gives the dimreg endpoint limit of c . f at a singular endpoint (drop rule + divergence gate). cvec entries are epsilon-free endpoint scalars; substitute the line variable before calling.";

Begin["`Private`"];

err[id_, payload_] := DiffExp2`Tolerances`DE2Error[id,
  Join[<|"Module" -> "API"|>, payload]];
cfg = DiffExp2`Config`CFG;
(* zeroQ: exact-first (see Solve.m).  Call sites test user coefficient
   functions: rational ones get the exact decision; anything outside the
   rational domain keeps PossibleZeroQ unchanged. *)
ratExprQ[c_] := Switch[c,
  _Integer | _Rational, True,
  _Complex, !InexactNumberQ[c],
  _Symbol, !NumericQ[c],
  _Plus | _Times, AllTrue[List @@ c, ratExprQ],
  Power[_, _Integer], ratExprQ[First[c]],
  _, False];
zeroCanQ[c_] := c === 0 || (!ratExprQ[c] && TrueQ[PossibleZeroQ[c]]);
zeroQ[e_] := zeroCanQ[Together[e]];
esAdd = DiffExp2`EpsSeries`ESAdd;
esZero = DiffExp2`EpsSeries`ESZero;

(* Geometry remains exact through planning/tiling.  At the integration seam,
   however, exact quadratic-algebraic local coordinates can grow into huge
   a+b Sqrt[n] quotients; asking N[Log[T],2 WP] to resolve their cancellation
   with the kernel's small default $MaxExtraPrecision is both slow and can
   lose most guard digits.  Root-reduce and ground once at 2x WP. *)
numericLocalCoordinate[z_] := Module[{wp, wp2},
  If[!NumericQ[z] || Precision[z] =!= Infinity, Return[z]];
  wp = cfg["WorkingPrecision"];
  wp2 = DiffExp2`Tolerances`$InputPrecisionFactor*wp;
  Block[{$MaxExtraPrecision = Max[$MaxExtraPrecision,
      DiffExp2`Tolerances`$MaxExtraPrecisionValue, 2 wp2]},
    Chop[SetPrecision[N[RootReduce[z], wp2], wp2],
      DiffExp2`Tolerances`Tol["ChopFloor"]]]];

(* ---- LoadSystem ---- *)

singularFactorsOf[A_, var_] := Module[{dens},
  dens = DeleteDuplicates[Flatten[Map[
    FactorList[Denominator[Together[#]]][[All, 1]] &, Flatten[A]]]];
  Select[dens, !FreeQ[#, var] &]];

LoadSystem[spec_Association] := Module[
  {A, var, exactFactors, plannerFactors},
  var = spec["Variable"];
  A = Which[
    KeyExistsQ[spec, "Matrix"], spec["Matrix"],
    KeyExistsQ[spec, "FullMatrixFile"],
    Module[{raw = Get[spec["FullMatrixFile"]]},
      If[!MatrixQ[raw], err["E6", <|"File" -> spec["FullMatrixFile"],
        "Detail" -> "full matrix file did not contain a matrix"|>]];
      raw],
    True, err["E6", <|"Keys" -> Keys[spec],
      "Detail" -> "LoadSystem needs \"Matrix\" or \"FullMatrixFile\""|>]];
  If[!FreeQ[A, _?InexactNumberQ],
    err["E6", <|"Detail" -> "inexact entries: the exact eps-rational matrix is required (d<var>_full.m, not eps-truncated slices)"|>]];
  If[!FreeQ[A, Power[b_, e_] /; !FreeQ[b, var] && !IntegerQ[e]],
    err["E6", <|"Detail" -> "irrational x-dependence (sqrt matrices) is out of scope v1"|>]];
  (* loading a system is the natural memory boundary for the solver memo
     caches (PrepareChart, SolveHomogeneous): entries for previous systems
     are dead weight from here on *)
  DiffExp2`Solve`ClearSolveCaches[];
  exactFactors = singularFactorsOf[A, var];
  plannerFactors = DiffExp2`Transport`EpsilonZeroSingularFactors[
    exactFactors, var];
  <|"Matrix" -> A, "Variable" -> var,
    "SingularFactors" -> plannerFactors,
    "SingularFactorsExact" -> exactFactors|>];

(* ---- TransportEndpoint ---- *)

Options[TransportEndpoint] = {"ExtraSingularFactors" -> {}};
TransportEndpoint[sys_Association, bvals_, from_, to_, OptionsPattern[]] := Module[
  {sys2, plan},
  sys2 = Join[sys, <|"ExtraSingularFactors" ->
    Select[OptionValue["ExtraSingularFactors"], !FreeQ[#, sys["Variable"]] &]|>];
  plan = DiffExp2`Transport`SegmentLine[sys2, {from, to}];
  DiffExp2`Transport`TransportLine[sys2, bvals, plan]];

(* ---- LineIntegral ----
   Integrate[c(x,eps).f(x), {x, lo, hi}]: transport from the anchor across
   [lo, hi] keeping charts; tile the interval by certified half-disk overlap
   ownership boundaries; per tile: MultiplyRational by
   each c-component IN CHART COORDINATES (x = center + scale t), including
   the affine dx = scale dt Jacobian, assemble the
   scalar LocalSolution, then integrate ONCE.  Combine-before-integrate is
   required when endpoint/PV divergences cancel between master components. *)

halfRadiusTiles[kept_List, lo_, hi_] := Module[
  {ordered = SortBy[kept, N[#["Chart", "Center"], 30] &], cs, rs, bps},
  cs = #["Chart", "Center"] & /@ ordered;
  rs = #["Chart", "Radius"] & /@ ordered;
  bps = Join[{lo}, Table[Module[{mid, lA, lB, raw, gap, cand},
      (* The transport error certificate probes each kept solution at R/2.
         Ownership boundaries must stay in the overlap of those certified
         half-disks; the former 0.9 R clamp was convergent but could make the
         reported error estimate nonconservative. *)
      lA = cs[[i]] + rs[[i]]/2;
      lB = cs[[i + 1]] - rs[[i + 1]]/2;
      If[TrueQ[N[lB, 40] > N[lA, 40]],
        err["E9", <|"Charts" -> {cs[[i]], cs[[i + 1]]},
          "Radii" -> {rs[[i]], rs[[i + 1]]},
          "HalfDiskGap" -> N[lB - lA, 20],
          "Detail" -> "adjacent kept charts have no certified half-radius tiling overlap"|>]];
      (* Preserve global monotonicity across unequal radii and duplicate
         anchor charts: start from the geometric center midpoint, then clip
         it into the certified half-disk overlap.  The overlap midpoint can
         move backward at the next pair even when centers are sorted. *)
      mid = (cs[[i]] + cs[[i + 1]])/2;
      raw = RootReduce[Max[Min[mid, lA], lB]];
      gap = N[lA - lB, 40];
      cand = If[TrueQ[gap > 0],
        Rationalize[N[raw, 30], gap/16], raw];
      If[TrueQ[N[lB, 40] <= N[cand, 40] <= N[lA, 40]], cand, raw]],
    {i, Length[cs] - 1}], {hi}];
  Do[If[!TrueQ[N[bps[[i]], 40] <= N[bps[[i + 1]], 40]],
    err["E9", <|"Breakpoints" -> {bps[[i]], bps[[i + 1]]},
      "Detail" -> "non-monotone tile breakpoints (tiling inversion)"|>]],
    {i, Length[bps] - 1}];
  Do[If[!TrueQ[Max[Abs[N[bps[[i]] - cs[[i]], 40]],
          Abs[N[bps[[i + 1]] - cs[[i]], 40]]] <= N[rs[[i]], 40]/2],
    err["E9", <|"Chart" -> cs[[i]], "Radius" -> rs[[i]],
      "TileBounds" -> {bps[[i]], bps[[i + 1]]},
      "Detail" -> "tile extends beyond the chart's certified half-radius envelope"|>]],
    {i, Length[cs]}];
  Table[{ordered[[i]], bps[[i]], bps[[i + 1]]}, {i, Length[cs]}]];

Options[LineIntegral] = {"ExtraSingularFactors" -> {},
  "PrecomputedCharts" -> None};
LineIntegral[sys_Association, bvals_, from_, {lo_, hi_}, cvec_List,
    OptionsPattern[]] := Module[
  {sys2, var = sys["Variable"], res, tiles, total = None,
   planLo, planHi, keptAll},
  If[Length[cvec] =!= Length[sys["Matrix"]],
    err["E9", <|"Detail" -> "line-integral coefficient vector has the wrong dimension",
      "Coefficients" -> Length[cvec], "SystemSize" -> Length[sys["Matrix"]]|>]];
  sys2 = Join[sys, <|"ExtraSingularFactors" ->
    Select[OptionValue["ExtraSingularFactors"], !FreeQ[#, var] &]|>];
  keptAll = OptionValue["PrecomputedCharts"];
  If[keptAll === None,
  (* transport BOTH ways from the anchor to cover [lo, hi] *)
  keptAll = {};
  If[TrueQ[lo < from],
    planLo = DiffExp2`Transport`SegmentLine[sys2, {from, lo}];
    res = DiffExp2`Transport`TransportLine[sys2, bvals, planLo];
    keptAll = Join[keptAll, res["Charts"]]];
  If[TrueQ[hi > from],
    planHi = DiffExp2`Transport`SegmentLine[sys2, {from, hi}];
    res = DiffExp2`Transport`TransportLine[sys2, bvals, planHi];
    keptAll = Join[keptAll, res["Charts"]]]];
  If[keptAll === {},
    err["E9", <|"Detail" -> "empty integration range or anchor outside"|>]];
  If[Environment["DEBUG_LI"] === "1",
    Print["    LI: transports done t=", SessionTime[], " ncharts=", Length[keptAll]]];
  tiles = halfRadiusTiles[keptAll, lo, hi];
  Do[Module[{entry = tile[[1]], a = tile[[2]], b2 = tile[[3]], ls, center,
      scale, t1, t2},
    If[TrueQ[b2 > a],
      ls = entry["LocalSolution"]; center = entry["Chart"]["Center"];
      scale = ls["ChartMap", "Scale"];
      t1 = numericLocalCoordinate[Together[(a - center)/scale]];
      t2 = numericLocalCoordinate[Together[(b2 - center)/scale]];
      (* Project each selected component to a scalar LocalSolution, multiply,
         then combine BEFORE integration so the endpoint/PV cancellation
         gate sees the assembled scalar integrand. *)
      Module[{pieces = {}},
        Do[Module[{cc = cvec[[ci]], lsP, lsM},
          If[!zeroQ[cc],
            lsP = Join[ls, <|"Sectors" -> Map[
              Join[#, <|"Coeffs" -> #["Coeffs"][[All, All, {ci}]]|>] &,
              ls["Sectors"]]|>];
            lsM = DiffExp2`SectorSeries`MultiplyRational[lsP,
              Together[scale*(cc /. var -> center + scale*Global`t)],
              Global`t];
            AppendTo[pieces, lsM];
            If[Environment["DEBUG_LI"] === "1",
              Print["      tile mul done t=", SessionTime[]]]]],
          {ci, Length[cvec]}];
        If[pieces =!= {},
          Module[{scalarLS, ri2, value},
            scalarLS = If[Length[pieces] === 1, First[pieces],
              DiffExp2`SectorSeries`CombineLocalSolutions[
                ConstantArray[1, Length[pieces]], pieces]];
            ri2 = DiffExp2`Integrate`IntegrateLocalSolution[scalarLS, {t1, t2}];
            If[Environment["DEBUG_LI"] === "1",
              Print["      tile int done t=", SessionTime[]]];
            value = ri2["Values"][[1]];
            total = If[total === None, value, esAdd[total, value]]]]]]],
    {tile, tiles}];
  If[total === None, esZero[cfg["EpsilonOrder"]], total]];

(* ---- EndpointLimitValues ---- *)

EndpointLimitValues[tres_Association, cvec_List] := Module[
  {ls = tres["Final"], pieces = {}, weights = {}, scalar},
  If[!TrueQ[tres["EndpointIsSingular"]],
    err["E9", <|"Detail" -> "EndpointLimitValues requires a singular-endpoint transport result"|>]];
  If[Length[cvec] =!= Dimensions[First[ls["Sectors"]]["Coeffs"]][[3]],
    err["E9", <|"Detail" -> "endpoint coefficient vector has the wrong dimension",
      "Coefficients" -> Length[cvec]|>]];
  (* Combine eps-free endpoint coefficients before applying the drop and
     divergence rules.  A limit of a scalar observable may exist even when
     its individual master terms do not. *)
  Do[If[!zeroQ[cvec[[ci]]],
    AppendTo[pieces, Join[ls, <|"Sectors" -> Map[
      Join[#, <|"Coeffs" -> #["Coeffs"][[All, All, {ci}]]|>] &,
      ls["Sectors"]]|>]];
    AppendTo[weights, cvec[[ci]]]],
    {ci, Length[cvec]}];
  If[pieces === {}, Return[esZero[ls["EpsWindow", "CompleteMax"]], Module]];
  scalar = DiffExp2`SectorSeries`CombineLocalSolutions[weights, pieces];
  DiffExp2`Integrate`EndpointSectorLimit[scalar][[1]]];

End[];
EndPackage[];
