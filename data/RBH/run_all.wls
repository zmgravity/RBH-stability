(* Evaluate with Get in Mathematica, or use verify.wls at the command line.
   This file never exits the kernel. Outputs are confined to outputs/.
   Every module must return an Association, register checks, and emit no
   messages. A failure in any module makes the complete run fail. *)
Module[{root = DirectoryName[$InputFileName], out, files, records = {},
  derived = <||>, result, file, elapsed, messages, before, checks, passed,
  report, hashes, sourceFiles, exportResult, exportOK, summary, manifest},
  out = FileNameJoin[{root, "outputs"}];
  If[!DirectoryQ[out], CreateDirectory[out]];
  Get[FileNameJoin[{root, "00_Common.wl"}]];
  TransportNoGo`ResetChecks[];
  files = {"01_QuadraticIdentity.wl", "02_TensorPrincipal.wl",
    "03_FoliationTurningPoint.wl", "04_HorndeskiDisformal.wl",
    "05_ExteriorExamples.wl", "06_NumericsFigures.wl", "07_ManuscriptChecks.wl"};
  Block[{$GeometricNoGoOutputDirectory = out},
    Do[
      Print["Running ", file];
      before = Length[TransportNoGo`CheckResults[]];
      Block[{$MessageList = {}},
        {elapsed, result} = AbsoluteTiming[
          CheckAbort[Get[FileNameJoin[{root, file}]], $Aborted]];
        messages = ToString[#, InputForm] & /@ $MessageList];
      AssociateTo[derived, file -> result];
      AppendTo[records, <|"File" -> file, "Seconds" -> elapsed,
        "Checks" -> (Length[TransportNoGo`CheckResults[]] - before),
        "ReturnedAssociation" -> AssociationQ[result], "Messages" -> messages|>];
      Print["  ", Last[records]["Checks"], " checks; ",
        If[AssociationQ[result] && messages === {}, "completed", "FAILED"],
        " (", ToString[NumberForm[elapsed, {8, 2}], OutputForm], " s)"], {file, files}]];
  checks = TransportNoGo`CheckResults[];
  passed = Length[checks] > 0 && And @@ Lookup[checks, "Passed"]
    && DuplicateFreeQ[Lookup[checks, "ID"]]
    && AllTrue[records, TrueQ[#["ReturnedAssociation"]] &&
      #["Checks"] > 0 && #["Messages"] === {} &];
  sourceFiles = Sort[Join[FileNames["*.wl", root], FileNames["*.wls", root]]];
  hashes = Association@Table[FileNameTake[file] ->
    IntegerString[FileHash[file, "SHA256"], 16, 64], {file, sourceFiles}];
  manifest = If[FileExistsQ[FileNameJoin[{root, "manuscript_reference.json"}]],
    Import[FileNameJoin[{root, "manuscript_reference.json"}], "RawJSON"], <||>];
  report = <|"Package" -> "Geometric obstruction: analytical reproduction",
    "PackageVersion" -> "1.1.3", "Kernel" -> $Version, "SystemID" -> $SystemID,
    "UTC" -> DateString[Now, {"ISODate", "T", "Time", "Z"}, TimeZone -> 0],
    "AllPassed" -> passed, "CheckCount" -> Length[checks],
    "UniqueCheckIDs" -> DuplicateFreeQ[Lookup[checks, "ID"]],
    "FailedCheckIDs" -> Lookup[Select[checks, !TrueQ[#["Passed"]] &], "ID", {}],
    "Modules" -> records, "SourceSHA256" -> hashes,
    "ManuscriptReference" -> manifest, "Checks" -> checks,
    "Scope" -> "Local analytical derivations and exterior benchmarks; global endpoint arguments and general functional-inverse hypotheses are specified in the paper."|>;
  summary = "Geometric obstruction - Mathematica reproducibility package\n" <>
    "Kernel: " <> $Version <> "\nAll passed: " <> ToString[passed] <>
    "\nChecks: " <> ToString[Length[checks]] <> "\n\n" <>
    StringRiffle[(#["ID"] <> " : " <> If[TrueQ[#["Passed"]], "PASS", "FAIL"]) & /@ checks, "\n"];
  exportOK = Check[
    Put[derived, FileNameJoin[{out, "derived_expressions.wl"}]];
    exportResult = Export[FileNameJoin[{out, "verification_summary.txt"}], summary, "Text"];
    StringQ[exportResult], False];
  AssociateTo[report, {"AllPassed" -> (passed && TrueQ[exportOK]),
    "OutputWritesSucceeded" -> TrueQ[exportOK]}];
  exportResult = Check[Export[FileNameJoin[{out, "verification_report.json"}], report, "RawJSON"], $Failed];
  If[!StringQ[exportResult], AssociateTo[report, "AllPassed" -> False]];
  Print["All passed: ", report["AllPassed"], "; checks: ", Length[checks]];
  Print["Report: ", FileNameJoin[{out, "verification_report.json"}]];
  report
]
