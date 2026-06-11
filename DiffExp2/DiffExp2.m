(* DiffExp2 loader: Get the modules in dependency order.
   Later milestones append SectorSeries, Indicial, Solve, Transport,
   Integrate, API. *)
Module[{dir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{dir, "Tolerances.m"}]];
  Get[FileNameJoin[{dir, "Config.m"}]];
  Get[FileNameJoin[{dir, "EpsSeries.m"}]];
  Get[FileNameJoin[{dir, "SectorSeries.m"}]];
  Get[FileNameJoin[{dir, "Indicial.m"}]];
  Get[FileNameJoin[{dir, "Solve.m"}]];
];
