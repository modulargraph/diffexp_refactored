(* DiffExp2 implementation loader: Get the modules in dependency order.
   The release/user loader is ../DiffExp2.m, which adds PublicAPI.m and
   installs the validated defaults. *)
Module[{dir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{dir, "Tolerances.m"}]];
  Get[FileNameJoin[{dir, "Config.m"}]];
  Get[FileNameJoin[{dir, "EpsSeries.m"}]];
  Get[FileNameJoin[{dir, "SectorSeries.m"}]];
  Get[FileNameJoin[{dir, "Indicial.m"}]];
  Get[FileNameJoin[{dir, "CppBackend.m"}]];
  Get[FileNameJoin[{dir, "Solve.m"}]];
  Get[FileNameJoin[{dir, "Transport.m"}]];
  Get[FileNameJoin[{dir, "NativeTransport.m"}]];
  Get[FileNameJoin[{dir, "Integrate.m"}]];
  Get[FileNameJoin[{dir, "API.m"}]];
];
