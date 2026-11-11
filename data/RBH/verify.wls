(* Command-line entry point. Returns nonzero for failed or incomplete runs. *)
report = Check[Get[FileNameJoin[{DirectoryName[$InputFileName], "run_all.wls"}]], $Failed];
Exit[If[AssociationQ[report] && TrueQ[report["AllPassed"]], 0, 1]];
