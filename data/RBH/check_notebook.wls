(* Execute actual saved notebook input boxes in order in a fresh kernel.
   Only NotebookDirectory[] is replaced, because no front end is open. *)
Module[{root = DirectoryName[$InputFileName], nb, boxes, held, value,
  messages, log = {}, n = 0, result, ok, written},
  nb = Get[FileNameJoin[{root, "Start_Here.nb"}]];
  boxes = Cases[nb, Cell[BoxData[b_], "Input", ___] :> b, Infinity];
  Do[
    n++;
    Block[{$MessageList = {}},
      held = ToExpression[box, StandardForm, HoldComplete];
      held = held /. HoldPattern[NotebookDirectory[]] -> root;
      value = CheckAbort[ReleaseHold[held], $Aborted];
      messages = ToString[#, InputForm] & /@ $MessageList];
    AppendTo[log, <|"Cell" -> n, "Messages" -> messages,
      "Aborted" -> (value === $Aborted)|>], {box, boxes}];
  ok = Length[boxes] > 0 && TrueQ[value] &&
    AllTrue[log, #["Messages"] === {} && !TrueQ[#["Aborted"]] &];
  result = <|"Kernel" -> $Version, "InputCells" -> Length[boxes],
    "AllPassed" -> ok, "FinalCellReturnedTrue" -> TrueQ[value], "Cells" -> log,
    "Method" -> "Saved input boxes evaluated sequentially; NotebookDirectory[] replaced by package path."|>;
  written = Export[FileNameJoin[{root, "outputs", "notebook_check.json"}], result, "RawJSON"];
  Print["Notebook cells: ", Length[boxes], "; all passed: ", ok];
  Exit[If[ok && StringQ[written], 0, 1]]
]
