# Mathematica reproducibility package

**A Geometric No-Go Theorem for Regular Black Holes in Scalar-Tensor Theories**  

## Execution

Extract the entire folder and open **Start_Here.nb** in Mathematica.
Evaluate its input cells in order. A fresh kernel is recommended.

From another notebook in the same directory:

```wl
Get[FileNameJoin[{NotebookDirectory[], "run_all.wls"}]]
```

From a terminal in the extracted folder:

```text
wolframscript -file verify.wls
```

If WolframScript does not automatically locate the local kernel, pass its
installed path with `-local`. For example, the validated Windows command was:

```powershell
& 'C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe' -local 'C:\Program Files\Wolfram Research\Wolfram\14.3\WolframKernel.exe' -file .\verify.wls
```

The wrapper exits with status 0 only when all checks pass, each module
returns an Association and registers checks, no module emits messages,
check IDs are unique, and the verification report is written.
`run_all.wls` never exits the kernel. Generated files replace the contents
with the same names in `outputs/`.

## Validated release

The delivered reports record **324 grouped checks** and **14 notebook input
cells**, executed using Mathematica 14.3.0 on Windows. Check groups may contain
whole matrices or operator lists; they are not counts of independent theorems.
Actual residuals, tolerances, kernel version, and source hashes are recorded
in the reports. The exported plot was visually inspected.

## Calculations

| Source | Calculation |
| --- | --- |
| `01_QuadraticIdentity.wl` | Constructs the metric connection, curvature and scalar Hessian; varies all quadratic HOST operators and checks the off-shell shape identity in both causal patches. |
| `02_TensorPrincipal.wl` | Derives the scalar-Hessian correction to the odd principal action, its gauge null vector and physical characteristic. The standard Einstein Fierz-Pauli action is an identified input. |
| `03_FoliationTurningPoint.wl` | Derives the arbitrary-lapse shape equation, undivided turning relation and odd shift-constraint reduction. |
| `04_HorndeskiDisformal.wl` | Varies the full Horndeski action, converts the odd characteristic to the interior frame and checks the proper-time turning identity. Also checks smooth scalar stalls, mixing factors, shape and clock transport, inverse maps and direct extensions. |
| `05_ExteriorExamples.wl` | Derives the extremal limit, Hayward exclusion threshold and arctangent identities. |
| `06_NumericsFigures.wl` | Solves original extremality equations at 40 and 70 digits, checks convergence and exterior bands, and exports the Figure 2 data and plot. |
| `07_ManuscriptChecks.wl` | Derives the k-essence Einstein and stress-tensor components. Tests the explicit U-DHOST example, ordinary-disformal Jacobian and speed transport, finite-order extremal zeros and the arctangent coefficient from the published action. |

`MANUSCRIPT_MAP.md` maps the printed equation numbers to computations and states
which conclusions use analytical reasoning or published inputs.
`SOURCES.txt` identifies these inputs. `00_Common.wl` records exact residuals
and numerical tolerances; unresolved expressions fail rather than pass.

## Numerical results and figures

Original metric horizon equations are solved at 40 and 70 digits. Metric
derivatives hold the mass and couplings fixed before substituting the horizon.
The tensor-coefficient calculation and metric-curvature calculation provide
two routes to the extremal radial speed within Mathematica:

```text
coefficient route: c_r,T^2 = K_T/F2
geometric route:   c_r,T^2 = -2/[r_e^2 B''(r_e)]
```

| Quantity | Value |
| --- | --- |
| Hayward upper mass ratio | `2 (2/5)^(2/3) = 1.08576704663796...` |
| Hayward tensor-zero radius ratio | `(5/2)^(1/3) = 1.35720880829745...` |
| Hayward extremal radial speed squared | `-1` |
| Arctangent p=3 extremal radial speed squared | `-0.67771376185156...` |
| Arctangent p=4 extremal radial speed squared | `-0.48277866773445...` |

Figure 2 data are regenerated from the analytic coefficient and horizon
functions. Panel (a) shows the finite ratio F2/KT from r/r_e = 1 through
7/4. Panel (b) shades the excluded exterior band and the region inside the
outer horizon. Both panels use a white background and the exact boundary
labels `(5/2)^(1/3)` and `2(2/5)^(2/3)`. Letter Figure 1 is a
schematic of the turning-point argument. The supplement's comparison of bounces
and wormholes uses the algebraic checks in module 04.

## Outputs

| File in `outputs/` | Contents |
| --- | --- |
| `verification_report.json` | Pass/fail records, residuals, module times, source hashes and manuscript metadata. |
| `verification_summary.txt` | Compact check registry. |
| `derived_expressions.wl` | Symbolic results from all seven modules. |
| `numerical_results.json` | Numerical results, working precisions and tolerances. |
| `extremal_hayward.csv`, `extremal_arctangent.csv` | Extremal roots, speeds and residuals. |
| `hayward_mass_scan.csv` | Horizon roots and exterior-band tests. |
| `fig2_horizon.csv`, `fig2_coefficient.csv`, `fig2_exclusion_band.csv` | Figure 2 data. |
| `hayward_exterior.pdf`, `hayward_exterior.png` | Regenerated quantitative figure. |
| `notebook_check.json` | Sequential execution of the saved notebook input boxes in a fresh kernel. |

CSV values use decimal/scientific notation. High-precision benchmark values
are decimal strings retaining up to 30 digits; the working precision and
tolerance are separate metadata.

## Scope and assumptions

The signature is (-,+,+,+), X = grad(phi)^2, and chi = -X/2. The scalar is
real and inherits the static and spherical symmetries. The quadratic shape
identity uses arbitrary couplings, before imposing degeneracy relations.
Their radial derivatives are retained.

The interior squared speed is F2/K_T. The exterior radial squared speed is K_T/F2.
The explicit Hayward and arctangent examples apply on their open, real-scalar
static exteriors. At F2=0 the exterior principal form degenerates, so the
figure displays the finite reciprocal F2/K_T.

The foliation proof uses a regular scalar clock. Scalar stalls are checked
in the covariant quadratic and Horndeski formulas with finite smooth coupling
jets and a nonzero limiting tensor coefficient. Degenerate coefficients are
handled with the undivided relation; 0/0 is not assigned a speed.

The compact-interval maximum proof, continuity to a near-extremal family,
general parity selection and smooth functional invertibility are analytical
arguments. The code verifies their local formulas and explicit map examples.
For mimetic additions, the verified scope is a seed with an off-shell shape
identity preserved by the added constraint; an on-shell definition of TPST
alone does not establish closure under this addition.
