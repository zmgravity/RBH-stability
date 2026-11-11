(* Letter eq:quadidentity (Eq. 12), Supplemental Material S3.
   Signature (-,+,+,+); X = grad(phi)^2, with no factor -1/2.
   Couplings c_i[x] are arbitrary functions ALONG the background. Their
   derivatives are retained. Since only N and J vary, both phi and X stay
   fixed: this includes arbitrary original dependence on (phi,X).
   No degeneracy conditions or scalar equation are used in this file. *)

Module[{x, tt, th, az, nn, jj, bb, zz, ph, sig, kap, eps, eta,
  c0, cq, cf, ca1, ca2, ca3, ca4, ca5, coupling, coords, metric, inverse,
  gamma, ricci, hessian, mixed, grad, kinetic, box, contraction, basis,
  expectedRicci, expectedEigen, sphericalRicci, sphericalBasis, angular,
  results = <||>, sign, curvature, geometry, measure, s, b, xx,
  operators, couplings, labels, varRules, vary, euler, lag, change,
  expected, current, source, kt, u, v, allL, offshell, nlog, jlog,
  metricSigned, assumption, k, ii, aa, mu, nu, rho},

  coords = {tt, x, th, az};
  expectedRicci = -2 bb[x] nn''[x]/nn[x] - bb'[x] nn'[x]/nn[x]
    + 2 kap/jj[x]^2 - 2 bb[x] jj'[x]^2/jj[x]^2
    - 4 (bb[x] jj''[x] + (bb[x] nn'[x]/nn[x] + bb'[x]/2) jj'[x])/jj[x];
  expectedEigen = {bb[x] nn'[x] ph'[x]/nn[x],
    bb[x] ph''[x] + bb'[x] ph'[x]/2,
    bb[x] jj'[x] ph'[x]/jj[x], bb[x] jj'[x] ph'[x]/jj[x]};

  (* Curvature and all five contractions are constructed from the metric.
     sig is nonzero; the connections are valid for either causal sign. *)
  Do[
    angular = Switch[curvature, 1, Sin[th], 0, 1, -1, Sinh[th]];
    metric = DiagonalMatrix[{-sig nn[x]^2, 1/bb[x], jj[x]^2,
      jj[x]^2 angular^2}];
    inverse = Inverse[metric];
    gamma = Table[Cancel[Sum[inverse[[aa, rho]]/2
      (D[metric[[rho, nu]], coords[[mu]]] + D[metric[[rho, mu]], coords[[nu]]]
       - D[metric[[mu, nu]], coords[[rho]]]), {rho, 4}]],
      {aa, 4}, {mu, 4}, {nu, 4}];
    ricci = FullSimplify[Sum[inverse[[mu, mu]] (
      Sum[D[gamma[[aa, mu, mu]], coords[[aa]]]
        - D[gamma[[aa, mu, aa]], coords[[mu]]], {aa, 4}]
      + Sum[gamma[[aa, aa, rho]] gamma[[rho, mu, mu]]
        - gamma[[aa, mu, rho]] gamma[[rho, mu, aa]], {aa, 4}, {rho, 4}]), {mu, 4}]];
    AssertZero["geometry/Ricci/kappa=" <> ToString[curvature],
      ricci - (expectedRicci /. kap -> curvature)];
    hessian = Table[D[ph[x], coords[[mu]], coords[[nu]]]
      - gamma[[2, mu, nu]] ph'[x], {mu, 4}, {nu, 4}];
    mixed = Simplify[inverse.hessian];
    AssertZero["geometry/Hessian/kappa=" <> ToString[curvature],
      mixed - DiagonalMatrix[expectedEigen]];
    grad = inverse.{0, ph'[x], 0, 0};
    kinetic = Cancel[grad.metric.grad];
    box = Tr[mixed];
    contraction = Cancel[grad.hessian.grad];
    basis = Cancel /@ {Tr[mixed.mixed], box^2, box contraction,
      grad.hessian.inverse.hessian.grad, contraction^2};
    AssertZero["geometry/operators/kappa=" <> ToString[curvature],
      basis - {Total[expectedEigen^2], Total[expectedEigen]^2,
        kinetic expectedEigen[[2]] Total[expectedEigen],
        kinetic expectedEigen[[2]]^2, kinetic^2 expectedEigen[[2]]^2}];
    If[curvature == 1, sphericalRicci = ricci; sphericalBasis = basis],
    {curvature, {1, 0, -1}}];

  couplings = {c0[x], cq[x], cf[x], ca1[x], ca2[x], ca3[x], ca4[x], ca5[x]};
  labels = {"F0", "Q", "F2", "A1", "A2", "A3", "A4", "A5"};

  (* A = sign N^2 and B = sign Z^2 avoid square-root branch ambiguity. *)
  Do[
    b = sign zz[x]^2; xx = b ph'[x]^2;
    s = nn[x]/zz[x]; measure = s jj[x]^2;
    assumption = nn[x] > 0 && zz[x] > 0 && jj[x] > 0 && 0 < th < Pi;
    metricSigned = DiagonalMatrix[{-sign nn[x]^2, 1/b, jj[x]^2, jj[x]^2 Sin[th]^2}];
    AssertZero["measure/sign=" <> ToString[sign],
      FullSimplify[Sqrt[-Det[metricSigned]]/Sin[th] - measure,
        Assumptions -> assumption]];
    geometry = Join[{1, Total[expectedEigen], sphericalRicci + 2 (kap - 1)/jj[x]^2},
      sphericalBasis] /. bb -> Function[{x}, sign zz[x]^2];
    operators = Cancel /@ geometry;
    varRules = Flatten[Table[{
      Derivative[k][nn][x] -> D[Exp[2 eta eps[x]] nn[x], {x, k}],
      Derivative[k][jj][x] -> D[Exp[-eta eps[x]] jj[x], {x, k}]}, {k, 0, 2}]];
    vary[expression_] := Cancel[Together[D[expression /. varRules, eta] /. eta -> 0]];
    euler[expression_, field_] := D[expression, field[x]]
      - D[D[expression, field'[x]], x] + D[D[expression, field''[x]], {x, 2}];
    AssertZero["variation/volume/sign=" <> ToString[sign], vary[measure]];
    AssertZero["variation/X/sign=" <> ToString[sign], vary[xx]];
    AssertZero["variation/arbitrary-phi-X-coupling/sign=" <> ToString[sign],
      vary[coupling[ph[x], xx]]];
    nlog = nn'[x]/nn[x]; jlog = jj'[x]/jj[x];

    Do[
      lag = measure couplings[[ii]] operators[[ii]];
      change = vary[lag];
      {expected, current, source} = Which[
        ii == 3, {4 kap s cf[x] eps[x] - 4 measure b cf[x] (nlog-jlog) eps'[x],
          2 measure b cf[x] (nlog-jlog), 4 kap s cf[x]},
        ii == 4, {4 measure b xx ca1[x] (nlog-jlog) eps'[x],
          -2 measure b xx ca1[x] (nlog-jlog), 0},
        True, {0, 0, 0}];
      AssertZero["variation/" <> labels[[ii]] <> "/sign=" <> ToString[sign], change - expected];
      AssertZero["Euler/" <> labels[[ii]] <> "/sign=" <> ToString[sign],
        2 nn[x] euler[lag, nn] - jj[x] euler[lag, jj] - 2 D[current, x] - source],
      {ii, Length[operators]}];

    kt = cf[x] - xx ca1[x]; u = sign nn[x]^2/jj[x]^2;
    v = jj[x]^4 kt D[u, x]/s;
    AssertZero["geometric-current/sign=" <> ToString[sign],
      v - 2 measure b kt (nlog-jlog)];
    allL = measure (couplings.operators);
    offshell = 2 nn[x] euler[allL, nn] - jj[x] euler[allL, jj];
    AssertZero["full-off-shell-identity/sign=" <> ToString[sign],
      offshell - 2 D[v, x] - 4 kap s cf[x]];
    AssociateTo[results, ToString[sign] -> <|"ReducedAction" -> allL,
      "Current" -> v, "Source" -> -2 kap s cf[x]|>],
    {sign, {1, -1}}];
  <|"RicciScalar" -> (expectedRicci /. kap -> 1),
    "HessianEigenvalues" -> expectedEigen, "SignedPatches" -> results,
    "OffShellIdentity" -> "2 N E_N - J E_J = 2 V' + 4 kappa s F2",
    "OnShellIdentity" -> "V' = -2 kappa s F2",
    "DegeneracyConditionsUsed" -> None|>
]
