(* FeynmanTrick Package Root Loader *)
(* Loads FeynmanTrick/FeynmanTrick.m from the same directory *)
Block[{$ContextPath},
  Get[FileNameJoin[{DirectoryName[$InputFileName], "FeynmanTrick", "FeynmanTrick.m"}]];
];
