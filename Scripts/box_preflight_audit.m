(* Box preflight, phase B (DiffExp2 only, no FIRE): read the vendored
   Tests/refs/bench/box_L<n>.m fixtures and statically audit everything the
   stepwise run will exercise, BEFORE any expensive transport:
     - SegmentLine plans for [anchor -> 0] and [anchor -> 1] at FT settings
       (WP 120, ExpansionOrder 40, DivisionOrder 4, classic stride 4):
       chart centers/radii/singular flags + the incoming match-point chain;
     - ValidatePlan on every plan (geometry E8s are the whole point);
     - FindSingularities layout vs the anchor and [0, 1];
     - apparent-singularity candidates: matrix denominator roots NOT among
       the IBP singular-factor roots;
     - per singular chart: PrepareChart -> IndicialData (pole order,
       FuchsianReduce trigger, spectrum, family classes, degeneracy).
   Output: AUDIT-prefixed JSON rows + readable prints.  Runs in ~a minute;
   no boundary values are touched. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

anchor = 11/23;
benchDir = FileNameJoin[{repoRoot, "Tests", "refs", "bench"}];
catch2[expr_] := Catch[expr, "DiffExp2Error"];
SetAttributes[catch2, HoldFirst];

(* recursive: containers walked, every leaf VALUE made JSON-encodable as a
   whole (Map[..., {-1}] only reaches atoms and leaves Plus/Abs heads in
   place -> Export::jsonstrictencoding) *)
jsonSafe[l_List] := jsonSafe /@ l;
jsonSafe[a_Association] := Map[jsonSafe, a];
jsonSafe[s_String] := s;
jsonSafe[True] = True; jsonSafe[False] = False; jsonSafe[None] = Null;
jsonSafe[e_] := Which[
  IntegerQ[e], e,
  Head[e] === Rational, N[e, 25],
  NumericQ[e] && Head[e] =!= Complex && FreeQ[e, _Symbol], N[e, 25],
  True, ToString[InputForm[e]]];

emit[tag_String, payload_Association] := Print["AUDIT ", ExportString[
  jsonSafe[payload], "RawJSON", "Compact" -> True]];

chartMP = DiffExp2`Transport`Private`chartMatchPoint;

auditLevel[level_Integer] := Module[
  {fix, var, A, extraFacs, sys, sys2, sings, realIn01, matrixFacs, ibpFacs,
   matrixRoots, ibpRoots, apparent, plans},
  fix = Get[FileNameJoin[{benchDir, "box_L" <> ToString[level] <> ".m"}]];
  var = fix["Variable"]; A = fix["Matrix"];
  extraFacs = fix["ExtraSingularFactors"];
  Print["\n================ LEVEL ", level, " var=", var, " d=", Length[A],
    " ================"];
  catch2[DiffExp2`Config`LoadConfiguration[{
    "WorkingPrecision" -> fix["FTSettings"]["WorkingPrecision"],
    "ExpansionOrder" -> fix["FTSettings"]["ExpansionOrder"],
    "EpsilonOrder" -> fix["FTSettings"]["EpsilonOrder"],
    "DivisionOrder" -> fix["FTSettings"]["DivisionOrder"],
    "StepDivisionOrder" -> fix["FTSettings"]["StepDivisionOrder"],
    "Variables" -> {}, "Verbosity" -> 0}]];
  sys = catch2[DiffExp2`API`LoadSystem[<|"Matrix" -> A, "Variable" -> var|>]];
  If[FailureQ[sys], emit["loadfail", <|"Level" -> level, "Fail" -> sys|>];
    Return[$Failed, Module]];
  sys2 = Join[sys, <|"ExtraSingularFactors" ->
    Select[extraFacs, !FreeQ[#, var] &]|>];

  (* ---- singularity layout ---- *)
  sings = catch2[DiffExp2`Transport`FindSingularities[sys2]];
  If[FailureQ[sings], emit["singfail", <|"Level" -> level, "Fail" -> sings|>];
    Return[$Failed, Module]];
  matrixFacs = sys["SingularFactors"];
  ibpFacs = Select[extraFacs, !FreeQ[#, var] &];
  matrixRoots = DeleteDuplicates[Flatten[Map[
    var /. Solve[# == 0, var] &, matrixFacs]],
    TrueQ[PossibleZeroQ[RootReduce[#1 - #2]]] &];
  ibpRoots = DeleteDuplicates[Flatten[Map[
    var /. Solve[# == 0, var] &, ibpFacs]],
    TrueQ[PossibleZeroQ[RootReduce[#1 - #2]]] &];
  apparent = Select[matrixRoots, Function[r,
    !AnyTrue[ibpRoots, TrueQ[PossibleZeroQ[RootReduce[# - r]]] &]]];
  realIn01 = Select[sings["Real"], TrueQ[0 <= N[#, 30] <= 1] &];
  emit["singularities", <|"Level" -> level,
    "MatrixFactors" -> matrixFacs, "IBPFactors" -> ibpFacs,
    "AllRoots" -> sings["All"], "RealRoots" -> sings["Real"],
    "RealRootsIn01" -> realIn01,
    "MatrixOnlyRoots_apparentCandidates" -> apparent,
    "DistancesToAnchor" -> Map[N[Abs[# - anchor], 20] &, sings["Real"]]|>];

  (* ---- plans + ValidatePlan + match-point chains ---- *)
  plans = Association[];
  Do[Module[{plan, vp, chain, prev = None},
    plan = catch2[DiffExp2`Transport`SegmentLine[sys2, {anchor, target}]];
    If[FailureQ[plan],
      emit["planfail", <|"Level" -> level, "Target" -> target,
        "Fail" -> plan|>],
      vp = catch2[DiffExp2`Transport`ValidatePlan[plan]];
      chain = Table[Module[{c = plan["Charts"][[ci]], mp},
        prev = If[ci === 1, None, plan["Charts"][[ci - 1]]];
        mp = chartMP[prev, c, plan["From"], plan["Direction"],
          DiffExp2`Config`CFG["DivisionOrder"]];
        If[!FreeQ[mp, chartMP],
          emit["mpfail", <|"Level" -> level, "Detail" ->
            "chartMatchPoint did not evaluate (context trap)"|>]];
        <|"Index" -> ci, "Name" -> c["Name"], "Singular" -> c["Singular"],
          "Center" -> c["Center"], "CenterN" -> N[c["Center"], 20],
          "Radius" -> c["Radius"], "RadiusN" -> N[c["Radius"], 20],
          "MatchPoint" -> mp,
          "MatchPointN" -> If[mp === None, None, N[mp, 20]]|>],
        {ci, Length[plan["Charts"]]}];
      emit["plan", <|"Level" -> level, "Target" -> target,
        "NCharts" -> Length[plan["Charts"]],
        "NSingularCharts" -> Count[plan["Charts"], c_ /; TrueQ[c["Singular"]]],
        "EndpointIsSingular" -> plan["EndpointIsSingular"],
        "DigitsNeeded" -> plan["DigitsNeeded"],
        "ValidatePlan" -> If[FailureQ[vp], vp, "OK"],
        "Chain" -> chain|>];
      plans[target] = plan]],
    {target, {0, 1}}];

  (* ---- indicial audit at every singular chart of both plans ---- *)
  Module[{singCharts},
    singCharts = DeleteDuplicatesBy[
      Select[Flatten[Table[Lookup[plans, t, <|"Charts" -> {}|>]["Charts"],
        {t, Keys[plans]}], 1], TrueQ[#["Singular"]] &],
      #["Center"] &];
    Do[Module[{cs, idata, pole, famSummary},
      cs = catch2[DiffExp2`Solve`PrepareChart[sys2, chart]];
      If[FailureQ[cs],
        emit["indicialfail", <|"Level" -> level,
          "Center" -> chart["Center"], "Fail" -> cs|>],
        idata = cs["IndicialData"];
        pole = idata["PoleData"];
        famSummary = Map[<|"Class" -> #["Class"],
          "JointSolve" -> #["JointSolve"], "LogMax" -> #["LogMax"],
          "EpsZeroDegeneracy" -> #["EpsZeroDegeneracy"],
          "NCollisions" -> Length[#["Collisions"]],
          "Sectors" -> Map[{#["a"], #["b"], #["p"]} &, #["Sectors"]]|> &,
          idata["Families"]];
        emit["indicial", <|"Level" -> level, "Center" -> chart["Center"],
          "PoleOrder" -> pole["PoleOrder"],
          "FuchsianReduceTriggered" -> (pole["PoleOrder"] >= 2),
          "GaugeSteps" -> idata["Reduction"]["Steps"],
          "Regular" -> idata["Regular"],
          "Spectrum" -> Map[{#["a"], #["b"], #["Multiplicity"]} &,
            idata["Spectrum"]],
          "Families" -> famSummary|>]]],
      {chart, singCharts}]];

  (* ---- boundary cases (from the fixture, for the report) ---- *)
  emit["cases", <|"Level" -> level,
    "Requests" -> Map[<|"MasterVec" -> #["MasterVec"], "Case" -> #["Case"],
      "Vi" -> #["Vi"], "Vj" -> #["Vj"], "NeededVec" -> #["NeededVec"]|> &,
      fix["Requests"]]|>];
  True];

Do[auditLevel[level], {level, 3, 1, -1}];
Print["\nAUDIT DONE t=", SessionTime[]];
Quit[0];
