(* FeynmanTrick release loader. *)
Module[{root = DirectoryName[$InputFileName]},
  If[!ValueQ[FeynmanTrick`$FeynmanTrickVersion],
    Block[{$ContextPath},
      Quiet[Get[FileNameJoin[{root, "FeynmanTrick", "FeynmanTrick.m"}]],
        {General::shdw, Symbol::shdw}]]];
  $ContextPath = DeleteDuplicates[Prepend[$ContextPath, "FeynmanTrick`"]];
];
