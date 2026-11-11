(* Shared exact checks. No external packages are required. *)
BeginPackage["TransportNoGo`"];
ResetChecks::usage = "ResetChecks[] starts a new verification run.";
AssertZero::usage = "AssertZero[id, expression] checks an exact algebraic residual.";
AssertTrue::usage = "AssertTrue[id, statement] records whether statement is exactly True.";
AssertNear::usage = "AssertNear[id, residual, tolerance] records a numerical residual and its absolute tolerance.";
CheckResults::usage = "CheckResults[] returns all checks, including failed residuals.";
Begin["`Private`"];
checks = {};
ResetChecks[] := (checks = {});
AssertZero[id_String, expression_] := Module[{residual, ok, reduce},
  reduce[e_] := If[ListQ[e], reduce /@ e, Cancel[Together[e]]];
  residual = reduce[expression];
  ok = And @@ (TrueQ[# == 0] & /@ Flatten[{residual}]);
  AppendTo[checks, <|"ID" -> id, "Passed" -> ok,
    "Residual" -> ToString[residual, InputForm]|>];
  If[!ok, Print["FAILED: ", id, " (residual retained in report)"]];
  ok
];
AssertTrue[id_String, statement_] := Module[{ok = TrueQ[statement]},
  AppendTo[checks, <|"ID" -> id, "Passed" -> ok,
    "Residual" -> ToString[statement, InputForm]|>];
  If[!ok, Print["FAILED: ", id, "  ", statement]];
  ok
];
CheckResults[] := checks;
AssertNear[id_String, residual_, tolerance_] := Module[{ok},
  ok = NumericQ[residual] && NumericQ[tolerance] && TrueQ[tolerance > 0]
    && TrueQ[Abs[residual] <= tolerance];
  AppendTo[checks, <|"ID" -> id, "Passed" -> ok,
    "Residual" -> ToString[residual, InputForm],
    "Tolerance" -> ToString[tolerance, InputForm]|>];
  If[!ok, Print["FAILED: ", id, " (numerical residual retained in report)"]];
  ok
];
End[];
EndPackage[];
