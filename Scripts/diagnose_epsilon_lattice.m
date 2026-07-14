repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Get[FileNameJoin[{repoRoot, "DiffExp2", "DiffExp2.m"}]];

SetAttributes[catchDE2, HoldFirst];
catchDE2[expr_] := Quiet[Catch[expr, "DiffExp2Error"]];

catchDE2[DiffExp2`Config`LoadConfiguration[{
  "WorkingPrecision" -> 100, "ExpansionOrder" -> 50,
  "EpsilonOrder" -> 4, "DivisionOrder" -> 4,
  "StepDivisionOrder" -> 4, "RecurrenceBackend" -> "Cpp",
  "Variables" -> {}}]];

fix = Get[FileNameJoin[{repoRoot, "Tests", "refs", "bench",
    "banana_L1.m"}]];
sys = catchDE2[DiffExp2`API`LoadSystem[
  <|"Matrix" -> fix["Matrix"], "Variable" -> fix["Variable"]|>]];
plan = catchDE2[DiffExp2`Transport`SegmentLine[sys, {11/23, 0}]];
chart = SelectFirst[Reverse[plan["Charts"]], TrueQ[#["Singular"]] &];
cs = catchDE2[DiffExp2`Solve`PrepareChart[sys, chart]];
blocks = DiffExp2`Solve`Private`sccBlockChartSystem[cs, #] & /@
  Range[Length[cs["IntegrationSequence", "Components"]]];

Print["components=", InputForm[cs["IntegrationSequence", "Components"]]];
Print["SEED_HALOS=", InputForm[
  DiffExp2`Solve`Private`sccSeedWorkHalos[cs, blocks, 56]]];
Do[
  data = DiffExp2`Solve`Private`clearedSymbolic[blocks[[block]]];
  certificate = DiffExp2`Solve`Private`epsilonRegularPrincipalCertificate[
    blocks[[block]], 56];
  framePlan = DiffExp2`Solve`Private`homogeneousFramePlan[blocks[[block]],
    <|"EpsWindow" -> <|"Min" -> -1, "CompleteMax" -> 4|>,
      "TOrder" -> 56|>];
  nvals = Map[DiffExp2`Solve`Private`epsValuation[#,
      DiffExp2`Config`CanonicalEps[]] &, data["NhatExpr"], {3}];
  dvals = Map[DiffExp2`Solve`Private`epsValuation[#,
      DiffExp2`Config`CanonicalEps[]] &, data["dExpr"]];
  Print["block=", block, " dim=", blocks[[block, "SystemSize"]],
    " poleDepth56=",
    DiffExp2`Solve`Private`recurrencePoleDepth[data, 56],
    " singleUse=",
    DiffExp2`Solve`Private`recurrenceSingleUsePoleDepth[data],
    " epsilonRegular=", InputForm[KeyDrop[certificate, "Symbolic"]],
    " framePlan=", InputForm[KeyTake[framePlan, {"PoleDepth",
      "InitialSourceDepth", "FrameTop", "TerminalFrameBase",
      "EpsilonRegularPrincipal"}]],
    " dvals=", InputForm[dvals],
    " nvals=", InputForm[nvals]];
  If[DiffExp2`Solve`Private`recurrencePoleDepth[data, 56] > 0,
    Print["OFFENDING_RESIDUE=",
      InputForm[Normal[blocks[[block, "Residue"]]]]];
    Print["OFFENDING_V=", InputForm[Normal[blocks[[block, "V"]]]]];
    Print["OFFENDING_V_INVERSE=",
      InputForm[Normal[blocks[[block, "VInv"]]]]];
    Print["OFFENDING_DET_V=",
      InputForm[Cancel[Together[Det[Normal[blocks[[block, "V"]]]]]]]];
    Print["OFFENDING_FAMILIES=",
      InputForm[blocks[[block, "Families"]]]];
    Print["OFFENDING_INDICIAL_FAMILIES=",
      InputForm[blocks[[block, "IndicialData", "Families"]]]];
    Print["OFFENDING_GAUGE=", InputForm[Normal[blocks[[block, "Gauge"]]]]];
    Print["OFFENDING_GAUGE_INVERSE=",
      InputForm[Normal[blocks[[block, "GaugeInverse"]]]]];
    Print["OFFENDING_THETA=",
      InputForm[Normal[blocks[[block, "ThetaMatrix"]]]]];
    Print["OFFENDING_D_EXPR=", InputForm[data["dExpr"]]];
    Print["OFFENDING_NHAT_EXPR=", InputForm[data["NhatExpr"]]]],
  {block, Length[blocks]}];

If[Environment["DE2_DIAGNOSE_COLUMN_PLANS"] === "1",
  prepared = catchDE2[DiffExp2`Solve`PrepareNativeSCCComposite[cs,
    <|"EpsWindow" -> <|"Min" -> -1, "CompleteMax" -> 4|>,
      "TOrder" -> 10|>]];
  records = Values[DiffExp2`Solve`Private`$nativeSCCCompositeCache];
  Print["PREPARED=", InputForm[prepared], " RECORD_COUNT=", Length[records]];
  If[records =!= {}, Print["COLUMN_PLANS=", InputForm[
    Last[records]["Execution", "ColumnPlans"]]]]];
