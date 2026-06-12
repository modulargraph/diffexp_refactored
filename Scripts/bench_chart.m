(* Scripts/bench_chart.m — per-chart pipeline micro-benchmark for DiffExp2.
   ======================================================================
   README
   ------
   Purpose: time the per-chart transport pipeline PHASE BY PHASE on the
   real FT level-1 systems WITHOUT the FT layer (no FeynmanTrick, no
   FIRE), so optimization work can iterate in seconds instead of hours.

   Fixtures (committed; regenerate only if the FT layer changes):
     Tests/refs/bench/sunrise_L1.m   d = 3 system in xx1
     Tests/refs/bench/banana_L1.m    d = 7 system in xx1
   generated ONCE by Scripts/gen_bench_fixtures.m (the only script that
   needs FeynmanTrick + FIRE).

   Settings are FT-realistic (run_ft_stepwise2.m at level 1):
     WorkingPrecision 120, ExpansionOrder 40, EpsilonOrder 5,
     DivisionOrder 4, StepDivisionOrder 4 (the runner's pinned stride).

   Phases (AbsoluteTiming, BENCH_REPS reps, MIN reported; the bench chart
   is the first interior REGULAR chart of the [11/23 -> 0] plan):
     Plan              SegmentLine for [11/23 -> 0] (+ chart count)
     PrepareChartCold  exact rational shift + ChartIndicial + V/VInv
                       (solve caches cleared before every rep)
     PrepareChartMemo  the $pcCache hit
     SolveChartCold    full d-column fundamental basis (includes the
                       always-on ODEResidualCheck = d+1 evaluations)
     SolveChartMemo    the $shCache hit (the lo/hi-shared-anchor replay)
     EvalBasis         d basis-column evaluations at the match point
     MatchWeights      eps-Laurent Gaussian elimination (d x d)
     Combine           CombineLocalSolutions of weights x basis
     EvalCombined      1 evaluation of the combined LocalSolution
     Validate          ValidateLocalSolution alone (the per-read scan
                       repeated by every Evaluate/Multiply call)
     ErrorProbe        SegmentErrorProbe (2 reduced/full evaluations)
     MultiplyRational  1 tile multiply by an IBP-like rational c(x)
     IntegrateTileCold IntegrateLocalSolution over a radius/8 tile with
                       the antiderivative memo flushed (first-touch cost)
     IntegrateTile     same, warm memo (the steady state across the
                       cvec components/masters that share a tile)
     TileSliced1Comp   warm multiply+integrate on a 1-component
                       projection (the LineIntegral per-component-waste
                       A/B: the runner keeps Values[[ci]] only)
     PrepareSingCold / SolveSingCold (reps = 1): the same for the plan's
                       singular endpoint chart (BENCH_SINGULAR=0 skips)

   SolveChart decomposition (uses DiffExp2`Solve`Private` internals — a
   benchmark-only seam; revisit if Solve.m refactors):
     PrepCleared       prepareCleared: PolynomialLCM + exact eps-expansion
                       of every cleared (t-degree, entry) pair
     Recursion1Col     ONE d-dimensional framed recursion column
                       (basis cost ~ d x this)
     Assemble1Col      assembleSolution for that column (V multiply + gauge)
     ResidualCheck     ODEResidualCheck on the full basis (d+1 evaluations)

   Derived rows:
     PerChartBasis  = PrepareChartCold + SolveChartCold + EvalBasis
                      + MatchWeights + Combine + EvalCombined + ErrorProbe
                      (the marching cost of ONE regular interior chart)
     LineEstimate   = Charts x PerChartBasis (upper bound; the anchor
                      chart's solve is shared lo/hi)

   Output: one compact JSON row per (fixture, phase), prefixed "BENCH ".
   Full default run target: <= ~3 min on one kernel.

   Env knobs:
     BENCH_FIXTURES=sunrise_L1,banana_L1   which fixtures
     BENCH_WP=120 BENCH_EXP_ORDER=40 BENCH_EPS_ORDER=5 BENCH_REPS=3
     BENCH_SINGULAR=1   also time the singular endpoint chart (1 rep)
     BENCH_TRANSPORT=1  additionally run ONE full TransportLine over the
                        plan (1 rep; adds minutes for banana — opt in)

   Run (under the shared kernel lock):
     env WolframKernel='/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel' \
       wolframscript -file Scripts/bench_chart.m
   ====================================================================== *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

envOr[name_, default_] := Module[{v = Environment[name]},
  If[StringQ[v] && StringLength[StringTrim[v]] > 0, v, default]];

wp = ToExpression[envOr["BENCH_WP", "120"]];
expOrd = ToExpression[envOr["BENCH_EXP_ORDER", "40"]];
epsOrd = ToExpression[envOr["BENCH_EPS_ORDER", "5"]];
reps = ToExpression[envOr["BENCH_REPS", "3"]];
doSing = envOr["BENCH_SINGULAR", "1"] === "1";
doTransport = envOr["BENCH_TRANSPORT", "0"] === "1";
fixNames = StringTrim /@ StringSplit[
  envOr["BENCH_FIXTURES", "sunrise_L1,banana_L1"], ","];

esNew = DiffExp2`EpsSeries`ESNew; esAdd = DiffExp2`EpsSeries`ESAdd;
esMin = DiffExp2`EpsSeries`ESMinPower; esCM = DiffExp2`EpsSeries`ESCompleteMax;
esCoeff = DiffExp2`EpsSeries`ESCoefficient;
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

row[fix_, phase_, secs_, extra_:<||>] := Print["BENCH ", ExportString[
  Join[<|"Fixture" -> fix, "Phase" -> phase,
    "Seconds" -> Round[N[secs], 1.*^-6]|>, extra],
  "RawJSON", "Compact" -> True]];

(* min-of-reps timing; clearer[] runs untimed before each rep *)
timeIt[body_, n_, clearer_:None] := Module[{ts},
  ts = Table[If[clearer =!= None, clearer[]]; First[AbsoluteTiming[body]], {n}];
  Min[ts]];
SetAttributes[timeIt, HoldFirst];

simpleRat[x_] := Rationalize[N[x, 20], Abs[N[x, 20]]/100];

benchFixture[fname_String] := Module[
  {path, fix, var, sys, sys2, plan, tPlan, charts, chart, dir, rad, tLoc,
   cs, req, sol, basis, d, Feval, F, vvals, w, ls, lsM, cc, ccT, sings,
   perChart = 0.0, ph},
  (* NB: never localize `t` here — Module would capture the literal
     Global`t passed to MultiplyRational (chart-variable shadowing) *)
  path = FileNameJoin[{repoRoot, "Tests", "refs", "bench", fname <> ".m"}];
  If[!FileExistsQ[path], Print["MISSING FIXTURE ", path,
    " — run Scripts/gen_bench_fixtures.m first"]; Quit[1]];
  fix = Get[path];
  var = fix["Variable"];
  Print["== fixture ", fname, " d=", Length[fix["Matrix"]],
    " var=", var, " t=", SessionTime[]];
  catch2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> wp, "ExpansionOrder" -> expOrd,
    "EpsilonOrder" -> epsOrd, "DivisionOrder" -> 4,
    "StepDivisionOrder" -> 4, "Variables" -> {}}]];
  sys = catch2[DiffExp2`API`LoadSystem[
    <|"Matrix" -> fix["Matrix"], "Variable" -> var|>]];
  If[FailureQ[sys], Print["LOAD FAIL ", sys]; Quit[1]];
  sys2 = Join[sys, <|"ExtraSingularFactors" ->
    Select[fix["ExtraSingularFactors"], !FreeQ[#, var] &]|>];

  (* ---- Plan: SegmentLine [11/23 -> 0] ---- *)
  tPlan = timeIt[plan = catch2[DiffExp2`Transport`SegmentLine[sys2, {11/23, 0}]], reps];
  If[FailureQ[plan], Print["PLAN FAIL ", plan]; Quit[1]];
  charts = plan["Charts"];
  row[fname, "Plan", tPlan, <|"Charts" -> Length[charts],
    "SingularCharts" -> Count[charts, c_ /; TrueQ[c["Singular"]]],
    "d" -> Length[fix["Matrix"]], "WP" -> wp, "ExpansionOrder" -> expOrd,
    "EpsilonOrder" -> epsOrd|>];

  (* bench chart: first interior REGULAR chart (the marching steady state) *)
  chart = SelectFirst[Rest[charts], !TrueQ[#["Singular"]] &, First[charts]];
  dir = plan["Direction"]; rad = chart["Radius"];
  tLoc = simpleRat[-dir*rad/4];   (* the incoming match-point offset *)
  Print["  bench chart ", chart["Name"], " radius=", N[rad, 6],
    " tLoc=", N[tLoc, 6]];

  (* ---- PrepareChart cold / memo ---- *)
  ph = timeIt[cs = DiffExp2`Solve`PrepareChart[sys2, chart], reps,
    DiffExp2`Solve`ClearSolveCaches];
  row[fname, "PrepareChartCold", ph]; perChart += ph;
  ph = timeIt[DiffExp2`Solve`PrepareChart[sys2, chart], reps];
  row[fname, "PrepareChartMemo", ph];

  (* ---- SolveChart cold / memo ---- *)
  d = cs["SystemSize"];
  req = <|"EpsWindow" -> <|"Min" -> 0, "CompleteMax" -> epsOrd|>,
    "TOrder" -> expOrd|>;
  ph = timeIt[sol = DiffExp2`Solve`SolveChart[cs, req], reps,
    Function[{}, DiffExp2`Solve`ClearSolveCaches[];
      cs = DiffExp2`Solve`PrepareChart[sys2, chart]]];
  row[fname, "SolveChartCold", ph]; perChart += ph;
  ph = timeIt[DiffExp2`Solve`SolveChart[cs, req], reps];
  row[fname, "SolveChartMemo", ph];
  basis = sol["Basis"]["Columns"];
  Print["  basis EpsWindow ", basis[[1]]["EpsWindow"]];
  If[Length[basis] =!= d, Print["BASIS SIZE MISMATCH ", Length[basis]]; Quit[1]];

  (* ---- SolveChart decomposition (Solve`Private seam; the bench chart is
     regular: Pmax = cdMax = 0, so the homogeneous frame is fb = -2,
     wideTop = epsOrd + 4 — the solveHomogeneousCore formula) ---- *)
  Module[{fb = -2, wideTop = epsOrd + 4, Wd, prep, init, rec1},
    Wd = wideTop - fb + 1;
    ph = timeIt[prep = DiffExp2`Solve`Private`prepareCleared[cs, fb, Wd], reps];
    row[fname, "PrepCleared", ph];
    init = {Module[{vv = Table[DiffExp2`EpsSeries`ESZero[wideTop], {d}]},
      vv[[1]] = esNew[0, PadRight[{1}, wideTop + 1]]; vv]};
    ph = timeIt[rec1 = DiffExp2`Solve`Private`runRecursion[cs, prep, 0, 0, 0,
      expOrd, None, fb, Wd, init], reps];
    row[fname, "Recursion1Col", ph, <|"FrameW" -> Wd|>];
    ph = timeIt[DiffExp2`Solve`Private`assembleSolution[cs, 0, 0, rec1, expOrd],
      reps];
    row[fname, "Assemble1Col", ph];
    ph = timeIt[DiffExp2`Solve`ODEResidualCheck[cs, sol["Basis"]], reps];
    row[fname, "ResidualCheck", ph]];

  (* ---- EvalBasis: d columns at the match point ---- *)
  ph = timeIt[
    Feval = Map[DiffExp2`SectorSeries`EvaluateLocalSolution[#, tLoc,
      "UsePade" -> False, "ImSign" -> 1]["Value"] &, basis], reps];
  row[fname, "EvalBasis", ph, <|"Columns" -> Length[basis]|>]; perChart += ph;

  (* component-major EpsSeries matrix + realistic incoming values
     (vvals = F . ones, so MatchWeights solves a genuine system) *)
  F = Table[esNew[esMin[Feval[[i]]],
      Table[esCoeff[Feval[[i]], k][[c]], {k, esMin[Feval[[i]]], esCM[Feval[[i]]]}]],
    {c, d}, {i, Length[basis]}];
  vvals = Table[Fold[esAdd, F[[c, 1]], Rest[F[[c]]]], {c, d}];

  (* ---- MatchWeights ---- *)
  ph = timeIt[w = catch2[DiffExp2`Transport`MatchWeights[F, vvals, "bench"]], reps];
  If[FailureQ[w], Print["MATCH FAIL ", w]; Quit[1]];
  row[fname, "MatchWeights", ph, Module[
    {fc = Flatten[Map[Function[s, Table[esCoeff[s, k],
        {k, esMin[s], esCM[s]}]], Flatten[{F, vvals}]]]},
    (* content diagnosis: exact-rational coefficients with large ByteCount
       force exact-giant field arithmetic inside the elimination *)
    <|"CoeffMaxBytes" -> Max[ByteCount /@ fc],
      "CoeffExactNonzero" -> Count[fc, x_ /; x =!= 0 && FreeQ[x, _?InexactNumberQ]],
      "CoeffTotal" -> Length[fc]|>]]; perChart += ph;

  (* ---- Combine ---- *)
  ph = timeIt[ls = DiffExp2`SectorSeries`CombineLocalSolutions[w, basis], reps];
  row[fname, "Combine", ph]; perChart += ph;

  (* ---- EvalCombined ---- *)
  ph = timeIt[DiffExp2`SectorSeries`EvaluateLocalSolution[ls, tLoc,
    "UsePade" -> False, "ImSign" -> 1], reps];
  row[fname, "EvalCombined", ph]; perChart += ph;

  (* ---- Validate: the per-read ValidateLocalSolution scan that EVERY
     Evaluate/Multiply call repeats (lever: validate-at-construction) ---- *)
  ph = timeIt[DiffExp2`SectorSeries`ValidateLocalSolution[ls], reps];
  row[fname, "Validate", ph];

  (* ---- ErrorProbe (the always-on per-chart accounting) ---- *)
  ph = timeIt[DiffExp2`Transport`SegmentErrorProbe[ls, tLoc, 0], reps];
  row[fname, "ErrorProbe", ph]; perChart += ph;

  (* ---- tile phase: MultiplyRational + IntegrateTile ----
     c(x) IBP-like: pole at a system singular factor whose roots all sit
     well OUTSIDE the bench chart's disk (>= 1.2 radius), else polynomial *)
  sings = DiffExp2`Transport`FindSingularities[sys2];
  cc = Module[{fac = SelectFirst[Keys[sings["Factors"]],
      Function[f, AllTrue[sings["Factors"][f],
        TrueQ[Abs[N[# - chart["Center"], 30]] > 1.2*N[rad, 30]] &]], None]},
    If[fac === None, 1 + var^2, (1 + var^2)/fac]];
  ccT = Together[cc /. var -> chart["Center"] + Global`t];
  Print["  tile multiplier c = ", InputForm[cc]];
  ph = timeIt[lsM = catch2[DiffExp2`SectorSeries`MultiplyRational[ls, ccT, Global`t]],
    reps];
  If[FailureQ[lsM], Print["MUL FAIL ", lsM]; Quit[1]];
  row[fname, "MultiplyRational", ph];
  Module[{t1 = simpleRat[-rad/8], t2 = simpleRat[rad/8], phCold},
    (* the antiderivative memo makes repeat tile integrals warm — the
       production steady state across cvec components/masters on a tile.
       Report cold (flushed) and warm (min of reps) separately. *)
    phCold = timeIt[catch2[DiffExp2`Integrate`IntegrateLocalSolution[lsM, {t1, t2}]],
      1, Function[{}, DiffExp2`Integrate`Private`$adCache = <||>]];
    ph = timeIt[catch2[DiffExp2`Integrate`IntegrateLocalSolution[lsM, {t1, t2}]],
      reps];
    row[fname, "IntegrateTileCold", phCold];
    row[fname, "IntegrateTile", ph];
    (* sliced variant: LineIntegral consumes ONE component per cvec entry
       (ri2["Values"][[ci]]) yet multiplies/integrates the FULL d-vector.
       Time the 1-component projection -> multiply -> integrate chain (the
       candidate API.m fix); identical arithmetic on the kept row. *)
    Module[{lsP, lsMP},
      lsP = Join[ls, <|"Sectors" -> Map[
        Append[#, "Coeffs" -> #["Coeffs"][[All, All, {1}]]] &,
        ls["Sectors"]]|>];
      ph = timeIt[
        lsMP = catch2[DiffExp2`SectorSeries`MultiplyRational[lsP, ccT, Global`t]];
        catch2[DiffExp2`Integrate`IntegrateLocalSolution[lsMP, {t1, t2}]],
        reps];
      row[fname, "TileSliced1Comp", ph]]];

  (* ---- derived ---- *)
  row[fname, "PerChartBasis", perChart];
  row[fname, "LineEstimate", perChart*Length[charts],
    <|"Charts" -> Length[charts]|>];

  (* ---- singular endpoint chart (1 rep) ---- *)
  If[doSing,
    Module[{sc = SelectFirst[Reverse[charts], TrueQ[#["Singular"]] &, None],
      cs2, ts1, ts2},
      If[sc =!= None,
        DiffExp2`Solve`ClearSolveCaches[];
        ts1 = First[AbsoluteTiming[
          cs2 = catch2[DiffExp2`Solve`PrepareChart[sys2, sc]]]];
        If[FailureQ[cs2], Print["SING PREP FAIL ", cs2],
          row[fname, "PrepareSingCold", ts1];
          ts2 = First[AbsoluteTiming[
            catch2[DiffExp2`Solve`SolveChart[cs2, req]]]];
          row[fname, "SolveSingCold", ts2]];
        (* restore caches for any later phase *)
        DiffExp2`Solve`ClearSolveCaches[]]]];

  (* ---- optional: ONE full TransportLine over the plan ---- *)
  If[doTransport,
    Module[{bvals, tr, tt},
      bvals = Table[N[(c + 1)/(c + k + 2), wp + 20], {c, d}, {k, 0, epsOrd}];
      tt = First[AbsoluteTiming[
        tr = catch2[DiffExp2`Transport`TransportLine[sys2, bvals, plan]]]];
      If[FailureQ[tr], Print["TRANSPORT FAIL ", tr],
        row[fname, "TransportLineFull", tt,
          <|"Charts" -> Length[charts],
            "MaxErr" -> Round[N[Max[tr["ErrorEstimate"]]], 0.0001]|>]]]];
  ];

Do[benchFixture[f], {f, fixNames}];
Print["BENCH DONE t=", SessionTime[]];
Quit[0];
