(* Structural pins for the extended Euclidean FeynmanTrick examples. *)

repoRoot = ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[repoRoot];
Quiet[Get["FeynmanTrick/FeynmanTrick.m"], {General::shdw, Symbol::shdw}];
Get[FileNameJoin[{repoRoot, "Scripts", "FTExamples.m"}]];

passed = 0; failed = 0;
assert[label_, cond_] := If[TrueQ[cond],
  passed++; Print["  PASS: ", label],
  failed++; Print["  FAIL: ", label]];

kite = FTExampleSpec["kite"];
banana4 = FTExampleSpec["banana4"];
banana4Unequal = FTExampleSpec["banana4_unequal"];
bananaUnequal = FTExampleSpec["banana_unequal"];
pentagonMassive = FTExampleSpec["pentagon_massive"];
bananaUnequalFT = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  FTExampleTopology["banana_unequal", "spec"],
  FTExampleSequence["banana_unequal"], {}];
banana4UnequalFT = FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
  FTExampleTopology["banana4_unequal", "spec"],
  FTExampleSequence["banana4_unequal"], {}];
anchorVectorFT = Block[{
    FeynmanTrick`Private`$FTConfig = Join[
      FeynmanTrick`Private`$FTConfig,
      <|"FixedParameterValues" -> {1/5, 3/10, 2/5}|>]},
  Fold[
    FeynmanTrick`FeynmanTrickIteration`BuildLevel[#1, #2] &,
    FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
      FTExampleTopology["banana_unequal", "anchor-vector-spec"],
      FTExampleSequence["banana_unequal"], {}],
    Range[3]]];
invalidAnchorVectorFT = Block[{
    FeynmanTrick`Private`$FTConfig = Join[
      FeynmanTrick`Private`$FTConfig,
      <|"FixedParameterValues" -> {1/5, 3/10}|>]},
  FeynmanTrick`FeynmanTrickIteration`DefineFTIteration[
    FTExampleTopology["banana_unequal", "invalid-anchor-vector-spec"],
    FTExampleSequence["banana_unequal"], {}]];

assert["ft_examples_extended_names",
  ContainsAll[FTExampleNames[],
    {"kite", "banana4", "banana4_unequal", "banana_unequal",
      "pentagon_massive"}]];
assert["ft_banana_unequal_euclidean_definition",
  Length[bananaUnequal["LoopMomenta"]] === 3 &&
  Length[bananaUnequal["ExternalMomenta"]] === 1 &&
  Length[bananaUnequal["Propagators"]] === 4 &&
  bananaUnequal["Replacements"] === {Global`p^2 -> -1} &&
  bananaUnequal["Dimension"] === 2 - 2 FeynmanTrick`FTeps &&
  ((# /. Thread[Join[bananaUnequal["LoopMomenta"],
      bananaUnequal["ExternalMomenta"]] -> 0]) & /@
      bananaUnequal["Propagators"]) === {2, 3/2, 4/3, 1}];
assert["ft_banana_unequal_selects_scalar_top_integral",
  bananaUnequalFT["Levels"][0]["Masters"] === {{1, 1, 1, 1}}];
assert["ft_kite_fully_massive_euclidean_definition",
  Length[kite["LoopMomenta"]] === 2 &&
  Length[kite["ExternalMomenta"]] === 1 &&
  Length[kite["Propagators"]] === 5 &&
  kite["Replacements"] === {Global`p^2 -> -1} &&
  kite["Dimension"] === 2 - 2 FeynmanTrick`FTeps &&
  AllTrue[kite["Propagators"],
    (# /. Thread[Join[kite["LoopMomenta"], kite["ExternalMomenta"]] -> 0]) ===
      1 &]];
assert["ft_banana4_four_loop_five_line_definition",
  Length[banana4["LoopMomenta"]] === 4 &&
  Length[banana4["ExternalMomenta"]] === 1 &&
  Length[banana4["Propagators"]] === 5 &&
  banana4["Replacements"] === {Global`p^2 -> -1} &&
  banana4["Dimension"] === 2 - 2 FeynmanTrick`FTeps];
assert["ft_banana4_unequal_euclidean_definition",
  Length[banana4Unequal["LoopMomenta"]] === 4 &&
  Length[banana4Unequal["ExternalMomenta"]] === 1 &&
  Length[banana4Unequal["Propagators"]] === 5 &&
  banana4Unequal["Replacements"] === {Global`p^2 -> -1} &&
  banana4Unequal["Dimension"] === 2 - 2 FeynmanTrick`FTeps &&
  ((# /. Thread[Join[banana4Unequal["LoopMomenta"],
      banana4Unequal["ExternalMomenta"]] -> 0]) & /@
      banana4Unequal["Propagators"]) === {2, 3/2, 4/3, 5/4, 1}];
assert["ft_banana4_unequal_selects_scalar_top_integral",
  banana4UnequalFT["Levels"][0]["Masters"] === {{1, 1, 1, 1, 1}}];
assert["ft_level_specific_anchors_fix_each_preceding_parameter",
  anchorVectorFT["FixedParamValues"] === {1/5, 3/10, 2/5} &&
    anchorVectorFT["Levels"][1]["FixedParams"] === {} &&
    anchorVectorFT["Levels"][2]["FixedParams"] ===
      {Global`xx1 -> 1/5} &&
    anchorVectorFT["Levels"][3]["FixedParams"] ===
      {Global`xx1 -> 1/5, Global`xx2 -> 3/10} &&
    invalidAnchorVectorFT === $Failed];
assert["ft_massive_pentagon_strictly_euclidean_definition",
  Length[pentagonMassive["LoopMomenta"]] === 1 &&
  Length[pentagonMassive["ExternalMomenta"]] === 4 &&
  Length[pentagonMassive["Propagators"]] === 5 &&
  pentagonMassive["Dimension"] === 4 - 2 FeynmanTrick`FTeps &&
  ((# /. Thread[Join[pentagonMassive["LoopMomenta"],
      pentagonMassive["ExternalMomenta"]] -> 0]) & /@
      pentagonMassive["Propagators"]) === {1, 3/2, 4/3, 5/4, 6/5} &&
  Sort[Eigenvalues[-Table[
      pentagonMassive["ExternalMomenta"][[i]]*
        pentagonMassive["ExternalMomenta"][[j]] /.
          pentagonMassive["Replacements"],
      {i, 4}, {j, 4}]]] === {1/4, 5/4, 5/4, 5/4} &&
  (Expand[(Total[pentagonMassive["ExternalMomenta"]])^2] /.
      pentagonMassive["Replacements"]) === -1];
assert["ft_extended_left_to_right_sequences",
  FTExampleSequence["kite"] === {{1, 2}, {1, 3}, {1, 4}, {1, 5}} &&
  FTExampleSequence["banana4"] === {{1, 2}, {1, 3}, {1, 4}, {1, 5}} &&
  FTExampleSequence["banana4_unequal"] ===
    {{1, 2}, {1, 3}, {1, 4}, {1, 5}} &&
  FTExampleSequence["pentagon_massive"] ===
    {{1, 2}, {1, 3}, {1, 4}, {1, 5}} &&
  FTExampleSequence["banana_unequal"] === {{1, 2}, {1, 3}, {1, 4}}];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
