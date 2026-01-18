(* DiffExp Package Loader *)
(* This file loads the DiffExp package from the DiffExp subdirectory *)

Block[{$ContextPath},
  $Path = Prepend[$Path, FileNameJoin[{DirectoryName[$InputFileName], "DiffExp"}]];
  Get["DiffExp.m"];
];
