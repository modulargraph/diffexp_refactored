(* DiffExp2 FT cutover runner: the M5 ladder.
   Mirrors Scripts/run_ft_stepwise.m but the transport/integration chain is
   DiffExp2 (sector-native, no fitting): per level the in-memory exact
   DiffMatrix is loaded directly (no slice export round-trip), boundary
   values are transported from the anchor 11/23 to both endpoints, and the
   boundary cases run through LineIntegral / EndpointLimitValues / direct
   convolution.  Output: STEPWISE/FINAL rows compatible with the old
   comparator. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
$Path = DeleteDuplicates[Prepend[$Path, repoRoot]];

Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

envOrDefault[name_, default_] := Module[{value = Environment[name]},
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, value, default]];

wp = ToExpression[envOrDefault["FT_WORKING_PRECISION", "100"]];
epsOrder = ToExpression[envOrDefault["FT_EPS_ORDER", "0"]];
expansionOrder = ToExpression[envOrDefault["FT_EXPANSION_ORDER", "40"]];
boundaryExtraOrder = ToExpression[envOrDefault["FT_BOUNDARY_EXTRA_ORDER", "4"]];
anchor = 11/23;

FeynmanTrick`SetFTOption["Threads", 1];
FeynmanTrick`SetFTOption["FThreads", 1];
FeynmanTrick`SetFTOption["Verbosity", 0];
FeynmanTrick`SetFTOption["WorkingPrecision", wp];
FeynmanTrick`SetFTOption["ReductionCache", True];
FeynmanTrick`SetFTOption["FIRETimeoutSeconds", 1800];

eps = Global`eps;
esC = DiffExp2`EpsSeries`ESCoefficient;
esMn = DiffExp2`EpsSeries`ESMinPower;
esCMx = DiffExp2`EpsSeries`ESCompleteMax;
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

cleanNumber[value_] := Module[{n = N[value, 50], re, im},
  If[!NumericQ[n], Return[ToString[InputForm[value]]]];  (* JSON-safe *)
  re = Re[n]; im = Im[n];
  If[Abs[im] < 10^-30, re, <|"Re" -> re, "Im" -> im|>]];

printRows[example_, level_, masters_, rawES_List, prefactors_] := Module[{},
  Do[Module[{r = rawES[[i]], rowMin},
    rowMin = esMn[r];
    Print["STEPWISE ", ExportString[<|
      "Example" -> example, "Level" -> level, "Master" -> masters[[i]],
      "EpsPrefactor" -> prefactors[[i]], "RawMinPower" -> rowMin,
      "Coefficients" -> Table[{p, cleanNumber[esC[r, p]]},
        {p, rowMin, Min[0, esCMx[r]]}]|> /. x_Rational :> N[x, 50],
      "RawJSON", "Compact" -> True]]],
    {i, Length[masters]}]];

(* combined endpoint limit: lim Sum_j c_j(x) f_j(x) at the chart center *)
limitCombined[tres_, cvec_, var_] := Module[{ls = tres["Final"], out = None},
  Do[Module[{cc = cvec[[j]], lsM, lim},
    If[!PossibleZeroQ[Together[cc]],
      lsM = DiffExp2`SectorSeries`MultiplyRational[ls,
        Together[cc /. var -> ls["Center"] + Global`t], Global`t];
      lim = DiffExp2`Integrate`EndpointSectorLimit[lsM][[j]];
      out = If[out === None, lim, DiffExp2`EpsSeries`ESAdd[out, lim]]]],
    {j, Length[cvec]}];
  If[out === None, DiffExp2`EpsSeries`ESZero[0], out]];

runExample[name_String] := Module[
  {topology, ftData, outputDir, nLevels, boundaryOrder, deepBoundary,
   currentBCs, currentPrefactors, finalRaw = None},
  Print["EXAMPLE ", name];
  FeynmanTrick`SetFTOption["DimensionExpression", FTExampleDimension[name]];
  topology = FTExampleTopology[name, "step"];
  If[topology === $Failed, Return[$Failed]];
  ftData = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    topology, FTExampleSequence[name], {}];
  outputDir = FileNameJoin[{$TemporaryDirectory,
    "FT2_" <> name <> "_" <> ToString[$ProcessID]}];
  If[DirectoryQ[outputDir], DeleteDirectory[outputDir, DeleteContents -> True]];
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True];
  ftData = FeynmanTrick`FeynmanTrickIteration`RunFullIteration[ftData, outputDir];
  If[ftData === $Failed, Return[$Failed]];
  nLevels = ftData["NumLevels"];
  boundaryOrder = epsOrder + nLevels + boundaryExtraOrder;
  deepBoundary = FeynmanTrick`BoundaryConditions`DeepestLevelBoundary[
    ftData, boundaryOrder];
  If[!AssociationQ[deepBoundary], Return[$Failed]];
  (* coefficients are the only numerics (tags stay exact): numericize the
     deep boundary at 2x WP so the chain runs at arbitrary precision instead
     of exact-symbolic (Log/Gamma giants grind the Laurent-field algebra) *)
  currentBCs = N[deepBoundary["BoundaryValues"], 2*wp];
  currentPrefactors = deepBoundary["EpsPrefactors"];
  printRows[name, nLevels, ftData["Levels"][nLevels]["Masters"],
    Table[DiffExp2`EpsSeries`ESShift[
      DiffExp2`EpsSeries`ESNew[0, currentBCs[[i]]], -currentPrefactors[[i]]],
      {i, Length[currentBCs]}],
    currentPrefactors];

  Module[{abortRes},
  abortRes = Catch[
  Do[Module[
    {levelData = ftData["Levels"][level], levelBelow = ftData["Levels"][level - 1],
     var, A, sys, mastersBelow, mastersHere, requests, neededVecs, reductions,
     extraFacs, rawES, trims, rawMin, shift, kmaxAvail, ftEps, dimVar,
     dimExpr, normalizeFT},
    var = levelData["FeynmanParameter"];
    (* normalize FT-layer symbols at the seam: dimension d -> 2-2eps form,
       FT epsilon symbol -> the DiffExp2 canonical Global`eps *)
    ftEps = FeynmanTrick`Private`$FTConfig["EpsilonSymbol"];
    dimVar = FeynmanTrick`Private`$FTConfig["DimensionVariable"];
    dimExpr = FeynmanTrick`Private`DimensionExpression[];
    normalizeFT[e_] := ((e /. dimVar -> dimExpr /. Global`d -> dimExpr) /.
      ftEps -> Global`eps);
    A = normalizeFT[levelData["DiffMatrix"]];
    Print["LEVEL ", level, " var=", var, " d=", Length[A]];
    sys = catch2[DiffExp2`API`LoadSystem[<|"Matrix" -> A, "Variable" -> var|>]];
    If[FailureQ[sys], Print["LOAD FAIL ", sys]; Return[$Failed, Module]];
    mastersHere = levelData["Masters"];
    mastersBelow = levelBelow["Masters"];
    requests = FeynmanTrick`DiffExpIntegration`Private`BoundaryRequestRecords[
      mastersBelow, levelData["CombinedPositions"]];
    neededVecs = DeleteDuplicates[#["NeededVec"] & /@ requests];
    reductions = FeynmanTrick`FIREInterface`ReduceIntegrals[
      levelData["Topology"], neededVecs];
    If[reductions === $Failed, Print["FIRE FAIL"]; Return[$Failed, Module]];
    extraFacs = normalizeFT[
      FeynmanTrick`DiffExpIntegration`CollectLevelIBPSingularFactors[
        ftData, level - 1]];
    (* configure DiffExp2 for this level *)
    catch2[DiffExp2`Config`LoadConfiguration[{
      "WorkingPrecision" -> wp, "ExpansionOrder" -> expansionOrder,
      "EpsilonOrder" -> Max[esCMxLevel = epsOrder + level + boundaryExtraOrder, 1],
      "DivisionOrder" -> 4, "Variables" -> {}}]];
    (* per lower master: dispatch the boundary case *)
    rawES = Table[Module[
      {req = requests[[mi]], expr2, cvecBase, cvec, case, res2, vi, vj, gammaFac},
      case = req["Case"]; vi = req["Vi"]; vj = req["Vj"];
      Print["  master ", mi, " case=", case, " vi=", vi, " vj=", vj,
        " t=", SessionTime[]];
      If[!KeyExistsQ[reductions, req["NeededVec"]],
        Print["MISSING REDUCTION ", req["NeededVec"]]; Return[$Failed, Module]];
      expr2 = normalizeFT[reductions[req["NeededVec"]]];
      cvecBase = Table[
        Together[Coefficient[expr2, Global`G[1, mastersHere[[j]]]]]/
          eps^currentPrefactors[[j]],
        {j, Length[mastersHere]}];
      Switch[case,
        "integrate",
        Module[{w},
          gammaFac = Gamma[vi + vj]/(Gamma[vi]*Gamma[vj]);
          cvec = Table[Together[
            var^(vi - 1)*(1 - var)^(vj - 1)*cvecBase[[j]]], {j, Length[mastersHere]}];
          w = catch2[DiffExp2`API`LineIntegral[sys, currentBCs, anchor, {0, 1},
            cvec, "ExtraSingularFactors" -> extraFacs]];
          If[FailureQ[w], Print["INTEGRATE FAIL master ", mi, ": ", w];
            Return[$Failed, Module]];
          DiffExp2`EpsSeries`ESScale[gammaFac, w]],
        "limitUpper",
        Module[{tr},
          tr = catch2[DiffExp2`API`TransportEndpoint[sys, currentBCs, anchor, 1,
            "ExtraSingularFactors" -> extraFacs]];
          If[FailureQ[tr], Print["LIMIT1 FAIL ", tr]; Return[$Failed, Module]];
          If[TrueQ[tr["EndpointIsSingular"]],
            limitCombined[tr, cvecBase, var],
            Module[{vv = tr["Value"], out = None},
              Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
                Module[{cES = DiffExp2`EpsSeries`ESFromExpression[
                    Together[cvecBase[[j]] /. var -> 1], eps, esCMx[vv]],
                    comp},
                  comp = DiffExp2`EpsSeries`ESNew[esMn[vv],
                    Table[esC[vv, k][[j]], {k, esMn[vv], esCMx[vv]}]];
                  Module[{term = DiffExp2`EpsSeries`ESTimes[cES, comp]},
                    out = If[out === None, term,
                      DiffExp2`EpsSeries`ESAdd[out, term]]]]],
                {j, Length[mastersHere]}];
              out]]],
        "limitLower",
        Module[{tr},
          tr = catch2[DiffExp2`API`TransportEndpoint[sys, currentBCs, anchor, 0,
            "ExtraSingularFactors" -> extraFacs]];
          If[FailureQ[tr], Print["LIMIT0 FAIL ", tr]; Return[$Failed, Module]];
          If[TrueQ[tr["EndpointIsSingular"]],
            limitCombined[tr, cvecBase, var],
            Module[{vv = tr["Value"], out = None},
              Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
                Module[{cES = DiffExp2`EpsSeries`ESFromExpression[
                    Together[cvecBase[[j]] /. var -> 0], eps, esCMx[vv]],
                    comp},
                  comp = DiffExp2`EpsSeries`ESNew[esMn[vv],
                    Table[esC[vv, k][[j]], {k, esMn[vv], esCMx[vv]}]];
                  Module[{term = DiffExp2`EpsSeries`ESTimes[cES, comp]},
                    out = If[out === None, term,
                      DiffExp2`EpsSeries`ESAdd[out, term]]]]],
                {j, Length[mastersHere]}];
              out]]],
        "direct",
        Module[{out = None, kmax = epsOrder + level + boundaryExtraOrder},
          Do[If[!PossibleZeroQ[Together[cvecBase[[j]]]],
            Module[{cES, bES, term},
              cES = DiffExp2`EpsSeries`ESFromExpression[
                Together[cvecBase[[j]]], eps, kmax];
              bES = DiffExp2`EpsSeries`ESNew[0, currentBCs[[j]]];
              term = DiffExp2`EpsSeries`ESTimes[cES, bES];
              out = If[out === None, term, DiffExp2`EpsSeries`ESAdd[out, term]]]],
            {j, Length[mastersHere]}];
          If[out === None, DiffExp2`EpsSeries`ESZero[0], out]]]],
      {mi, Length[mastersBelow]}];
    If[MemberQ[rawES, $Failed], Throw[$Failed, "FT2Abort"]];
    rawES = DiffExp2`EpsSeries`ESTrim /@ rawES;
    printRows[name, level - 1, mastersBelow, rawES,
      ConstantArray[0, Length[mastersBelow]]];
    (* shift to finite for the next level's transport *)
    rawMin = Min[esMn /@ rawES];
    shift = Max[0, -rawMin];
    currentBCs = Table[Module[{r = rawES[[i]], top},
      top = Min[esCMx[r], epsOrder + (level - 1) + boundaryExtraOrder - shift];
      Table[esC[r, k], {k, -shift, top}]],
      {i, Length[rawES]}];
    currentPrefactors = ConstantArray[shift, Length[rawES]];
    finalRaw = rawES],
    {level, nLevels, 1, -1}], "FT2Abort"];
  If[abortRes === $Failed, Return[$Failed]]];
  If[finalRaw === None || MemberQ[finalRaw, $Failed], Return[$Failed]];

  Print["FINAL ", ExportString[<|
    "Example" -> name,
    "Finite" -> cleanNumber[esC[finalRaw[[1]], 0]],
    "RawMinPower" -> esMn[finalRaw[[1]]]|> /. x_Rational :> N[x, 50],
    "RawJSON", "Compact" -> True]];
  True];

requested = StringTrim /@ StringSplit[envOrDefault["FT_EXAMPLES", "bubble"], ","];
Do[
  If[runExample[name] === $Failed, Print["FAILED ", name]; Quit[1]],
  {name, requested}];
Quit[0];
