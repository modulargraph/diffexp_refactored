(* DiffExp2 release loader.

   The implementation loader remains available at DiffExp2/DiffExp2.m for
   low-level development tests.  Users should Get this root file: it installs
   the umbrella DiffExp2` API and a validated default configuration. *)

Module[{root = DirectoryName[$InputFileName]},
  If[!ValueQ[DiffExp2`$DiffExp2Version],
    Block[{$ContextPath},
      Get[FileNameJoin[{root, "DiffExp2", "DiffExp2.m"}]];
      Quiet[Get[FileNameJoin[{root, "DiffExp2", "PublicAPI.m"}]],
        {General::shdw, Symbol::shdw}]]];
  If[!TrueQ[DiffExp2`Config`ConfiguredQ[]],
    DiffExp2`LoadConfiguration[]];
  $ContextPath = DeleteDuplicates[Prepend[$ContextPath, "DiffExp2`"]];
];
