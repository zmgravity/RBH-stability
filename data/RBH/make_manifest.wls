(* Rebuild the release manifest after a successful notebook execution.
   Run this last: the manifest covers all delivered files except itself. *)
Module[{root=DirectoryName[$InputFileName],report,notebook,files,relative,hashes,result,sourceHashes,
  releaseFiles,manifestPath},
  report=Import[FileNameJoin[{root,"outputs","verification_report.json"}],"RawJSON"];
  notebook=Import[FileNameJoin[{root,"outputs","notebook_check.json"}],"RawJSON"];
  If[!AssociationQ[report]||!AssociationQ[notebook]||
    !TrueQ[report["AllPassed"]]||!TrueQ[notebook["AllPassed"]],
    Print["Successful verification and notebook reports are required."];Exit[1]];
  sourceHashes=Lookup[report,"SourceSHA256",<||>];
  If[Length[sourceHashes]==0 || !And@@KeyValueMap[
    Function[{name,digest},With[{path=FileNameJoin[{root,name}]},
      FileExistsQ[path] && IntegerString[FileHash[path,"SHA256"],16,64]===digest]],sourceHashes],
    Print["Sources changed after verification. Rerun check_notebook.wls first."];Exit[1]];
  (* A checkout can contain a .git directory or a worktree's .git file. *)
  releaseFiles[dir_] := Flatten[Map[If[DirectoryQ[#],releaseFiles[#],#]&,
    Select[FileNames["*",dir],FileNameTake[#]=!=".git"&]]];
  manifestPath=FileNameJoin[{root,"MANIFEST.json"}];
  files=Sort[Select[releaseFiles[root],FileType[#]===File && #=!=manifestPath&]];
  relative[path_]:=StringReplace[
    FileNameJoin[Drop[FileNameSplit[path],Length[FileNameSplit[root]]]],"\\"->"/"];
  hashes=Association@Table[relative[file]->IntegerString[FileHash[file,"SHA256"],16,64],{file,files}];
  result=Export[manifestPath,
    <|"version"->report["PackageVersion"],"date"->report["ManuscriptReference"]["package_date"],
      "mathematica_checks"->report["CheckCount"],"notebook_cells"->notebook["InputCells"],
      "files"->hashes|>,"RawJSON"];
  Print["Manifest written for ",Length[files]," files."];
  Exit[If[StringQ[result],0,1]]
]
