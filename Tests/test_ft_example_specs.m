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

assert["ft_examples_extended_names",
  ContainsAll[FTExampleNames[], {"kite", "banana4"}]];
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
assert["ft_extended_left_to_right_sequences",
  FTExampleSequence["kite"] === {{1, 2}, {1, 3}, {1, 4}, {1, 5}} &&
  FTExampleSequence["banana4"] === {{1, 2}, {1, 3}, {1, 4}, {1, 5}}];

Print["Results: ", passed, " / ", passed + failed, " tests passed"];
If[failed > 0, Quit[1], Quit[0]];
