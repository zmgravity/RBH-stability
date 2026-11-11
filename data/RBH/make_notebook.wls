(* Build a portable, annotated notebook without launching the front end. *)
Module[{root = DirectoryName[$InputFileName], input, cells, nb},
  SetAttributes[input, HoldAllComplete];
  input[expr_] := Cell[BoxData[ToBoxes[Defer[expr]]], "Input"];
  cells = {
    Cell["A Geometric No-Go Theorem for Regular Black Holes in Scalar-Tensor Theories", "Title"],
    Cell["Analytical reproduction and numerical consistency checks | Version 1.1.3", "Subtitle"],
    Cell["Open this notebook from the extracted package and evaluate its input cells in order. The calculation uses the adjacent Wolfram Language sources. It writes the verification report, exact expressions, numerical tables and a rendering of Letter Figure 2 to outputs/.", "Text"],
    Cell["1. Conventions and complete run", "Section"],
    Cell["Signature (-,+,+,+), X = grad(phi)^2, chi = -X/2. The scalar is real and inherits the spherical and Killing symmetries. K_T = F2 - X A1. The interior squared speed is F2/K_T. The exterior radial squared speed is K_T/F2.", "Text"],
    input[base = NotebookDirectory[]; report = Get[FileNameJoin[{base, "run_all.wls"}]];],
    input[KeyTake[report, {"Kernel", "AllPassed", "CheckCount", "FailedCheckIDs"}]],
    input[Dataset[report["Modules"]]],
    input[derived = Get[FileNameJoin[{base, "outputs", "derived_expressions.wl"}]];],
    Cell["2. General quadratic action: Letter Eq. (12), SM S3", "Section"],
    Cell["The source constructs the connection, curvature and five Hessian contractions directly from the spherical metric. It varies each action term under a volume-preserving shape change and independently checks the Euler-Lagrange combination. Couplings are arbitrary functions along the background; their radial derivatives are retained. The shape variation fixes phi and X, so this covers arbitrary original dependence on (phi,X). Both metric-sign patches are checked.", "Text"],
    input[KeyTake[derived["01_QuadraticIdentity.wl"], {"RicciScalar", "HessianEigenvalues", "OffShellIdentity", "OnShellIdentity"}]],
    Cell["3. Quadratic tensor principal symbol: SM S3", "Section"],
    Cell["The tensor calculation starts from the standard Einstein Fierz-Pauli principal action and constructs the scalar-Hessian perturbation. It checks the gauge null vector and physical characteristic factor in timelike and spacelike scalar-gradient frames. The Einstein principal action is the stated standard input.", "Text"],
    input[derived["02_TensorPrincipal.wl"]],
    Cell["4. Foliation dynamics and the turning point: Letter Eqs. (6)-(10), SM S1", "Section"],
    Cell["The foliation module derives the shear invariants and shape variation, then evaluates the undivided equation at zero shear. The odd shift constraint removes the longitudinal gauge mode. The resulting kinetic and gradient coefficients are the same ones that govern the shape acceleration.", "Text"],
    input[derived["03_FoliationTurningPoint.wl"]],
    Cell["5. Horndeski, field transformations and mixing: SM S2, S4-S6", "Section"],
    Cell["The module checks the Horndeski shape equation, the interior turning-point speed, the clock and shape transformations, and the explicit generalized-disformal inverse. It also checks the cosmological Jacobian and the wormhole mixing limit in S2. MANUSCRIPT_MAP.md identifies the analytical and published inputs.", "Text"],
    input[derived["04_HorndeskiDisformal.wl"]],
    Cell["6. Extremal exterior and examples: Letter Eqs. (15)-(17), SM S7", "Section"],
    Cell["The exact calculation derives the extremal coefficient relation and the Hayward exclusion threshold, then treats the arctangent profile. Numerical tests compare coefficients from the supporting action against independently computed horizon geometry, at two working precisions. These examples use their real-scalar exterior domains.", "Text"],
    input[derived["05_ExteriorExamples.wl"]],
    input[derived["06_NumericsFigures.wl"]],
    Cell["7. Current manuscript: additional exact checks", "Section"],
    Cell["K-essence Einstein components, the explicit U-DHOST example, ordinary-disformal inverse and Jacobian, full speed transport, finite-order extremal coefficient zeros, and the arctangent action dictionary. See MANUSCRIPT_MAP.md for regularity and seed-action assumptions.", "Text"],
    input[derived["07_ManuscriptChecks.wl"]],
    Cell["8. Inspect checks and files", "Section"],
    input[Dataset[report["Checks"]]],
    input[FileNames["*", FileNameJoin[{base, "outputs"}]]],
    Cell["The verification report records the kernel version, execution time, checks, numerical residuals and source hashes. Numerical agreement supports the explicit calculations; the global endpoint proof and the hypotheses for general field transformations remain mathematical arguments in the Letter and supplement.", "Text"],
    input[TrueQ[report["AllPassed"]]]
  };
  nb = Notebook[cells, WindowTitle -> "Geometric obstruction - reproducibility",
    StyleDefinitions -> "Default.nb", Saveable -> True];
  Put[nb, FileNameJoin[{root, "Start_Here.nb"}]];
  Print["Notebook written; input cells: ", Length[Cases[cells, Cell[_, "Input", ___]]]];
]
