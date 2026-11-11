# Manuscript equations and Mathematica calculations

The file hashes are in
`manuscript_reference.json`. Module numbers refer to the numbered `.wl` files.

| Printed location | Calculation and scope |
| --- | --- |
| Letter (1)–(3), Fig. 1. SM S1 A, (S1)–(S2) | Modules 01 and 03 use the spherical metric and the scale ratio. With positive interior ratio, continuous zero endpoint values and a regular proper-time coordinate at the extremum, the existence of a positive interior maximum follows analytically. Figure 1 illustrates this argument. |
| Letter (4)–(5). SM S1 A, (S3)–(S4) | Module 07 constructs the Einstein tensor, the covariant and contravariant scalar gradients, and the mixed k-essence stress tensor. Equality of the three spatial pressures gives the shape equation. The mixed temporal component is `P - 2 X P_X`. |
| Letter (6)–(10). SM S1 B–C, (S5)–(S11) | Module 03 derives the shape equation with an arbitrary lapse and reduces the odd shift constraint at zero shear. It evaluates the undivided turning relation at a stationary point. If the kinetic coefficient vanishes, the gradient coefficient vanishes there as well. The scalar foliation requires a regular timelike scalar gradient. |
| SM S2, (S12) and the mixing comparison | The cosmological monotonicity relation is an input from SM Ref. [4]. Module 04 checks the pure-disformal kinetic Jacobian and its numerator using the cosmological relation in SM Ref. [6], Eqs. (23)–(24). It checks the Horndeski limit of the separate static wormhole mixing factor in SM Ref. [5]. The even-sector stability arguments are given in those references. |
| Letter (11)–(12). SM S3 A, (S13)–(S20) | Module 01 constructs and varies all five quadratic operators and the curvature term in both causal regions. For spherical slices it finds `delta R = 4 epsilon/J^2 - 4 B(n-j) epsilon'` and the stated off-shell radial identity with arbitrary coupling profiles. |
| Letter (9), (12). SM S3 B, (S21) | Module 02 computes the scalar-Hessian correction and the gauge-reduced characteristic in both causal regions using the Einstein Fierz–Pauli principal action as input. The interior squared speed is `F2/KT`. The exterior radial squared speed is `KT/F2`. |
| Letter after (12). SM S3 C, (S22) | Module 07 verifies the explicit U-DHOST coefficient choice, its unitary degeneracy condition, the homogeneous cancellation and the tensor coefficients. For regular covariant couplings the quadratic scalar-stall limit gives `KT=F2`, with unit squared speed if `F2 != 0`. Module 04 checks the corresponding smooth Horndeski limit. |
| SM S4, (S23)–(S27) | Module 04 varies the full G2–G5 action in both causal patches and checks an independent combination of published background equations. With the normalization of (S27), the quadratic dictionary is `(KT,F2) = (H/2,F/2)`. |
| SM S4, (S28)–(S29) | Module 04 checks the G−H factorization and equality at the turning point. Using the published exterior ratios G/F and G/H, it converts the characteristic to the interior proper-time frame. The interior longitudinal and angular squared speeds are F/G and F/H, respectively. They coincide at F/H when G=H and the coefficients are nonzero. The code also converts the radial current to the undivided proper-time relation (S29). |
| Letter (13)–(14). SM S5 A, (S30)–(S31) | Module 04 checks the clock/shape transformation and congruence algebra. Module 07 checks the orthonormal speed rescaling. The general argument uses solution correspondence, parity closure, and invariance of inertia under an invertible congruence. |
| SM S5 B, (S32) | Module 07 verifies the inverse metric, kinetic invariant, Jacobian for C(X), D(X), and the algebraic odd map. |
| SM S5 C–D, (S33)–(S38) | Smooth functional invertibility, positive finite C and E, and odd-sector parity selection give the conditional result (S33). For the explicit map (S34), module 04 checks (S35)–(S36), including the inverse, determinant and composition, in both gradient signatures. It also checks the Weyl inverse Jacobian and action in (S37)–(S38), using the four-dimensional Weyl conformal weight. |
| SM S6, Table I and (S39)–(S40) | Module 04 checks the spatial-gradient example with homogeneous and linear-odd gradient zeros, and the vanishing shape variation of the added mimetic constraint. The constraint contains no odd principal derivatives. The extension applies to seeds whose turning relation follows from the preserved shape equation. |
| Letter (15). SM S7 A, (S41)–(S45) | Module 05 expands the quadratic HOST radial identity at a smooth double-zero outer horizon. Module 07 derives the finite-order leading coefficient and negative exterior radial squared speed for arbitrary nonnegative integer order and checks Taylor limits at orders 0 through 3. Persistence of an exterior interval in a jointly continuous family is an analytical argument. |
| Letter (16)–(17), Fig. 2. SM S7 B, (S46)–(S48) | Module 05 obtains the Hayward tensor coefficients, mass curve and exact exclusion threshold. Module 06 solves the original horizon equations, tests the exterior band and exports the two panels of Letter Figure 2 with their numerical data. Panel (a) shows F2/KT and panel (b) shows the excluded exterior interval as mass varies. |
| SM S7 C, (S49)–(S51) | Module 05 checks the arctangent background identity, extremality conditions and ingredients of the sign argument for p ≥ 3. Module 07 obtains KT from the published action coefficients. Module 06 independently solves the original metric equations at 40 and 70 digits for p=3 and p=4. |

## Reading the reports

`AssertZero` simplifies each scalar or matrix residual and requires every
entry to vanish. `AssertTrue` requires the result `True`. `AssertNear` records
a numerical residual and its tolerance. The run also checks that each module
returns an association, registers checks, and finishes without messages or
aborts. Each check has a unique ID in `outputs/verification_report.json`.
