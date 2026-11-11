(* Letter eq:oddchar and eq:quadidentity; Supplemental Material S3.
   Freeze the local background in an orthonormal frame. Keep h_{a3}=u_a,
   a=0,1,2, and partial_3=0. This is the odd polarization, NOT the full
   scalar/even perturbation problem. The Einstein term starts from its
   standard Fierz-Pauli principal action, up to a common normalization. *)

Module[{eta = {-1,1,1,1}, du, uu, kk, p, g, aone, dh, f01, f02, f12,
  eh, maxwell, fourierRules, normal, label, xx, kt, deltaH, deltaL1,
  action, expected, symbol, matrix, cone, adj, results = <||>,
  aa, bb, dd, mu, nu},
  dh = Table[Which[mu <= 3 && nu == 4 && dd <= 3, du[mu,dd],
    nu <= 3 && mu == 4 && dd <= 3, du[nu,dd], True, 0],
    {dd,4}, {mu,4}, {nu,4}];
  f01 = du[2,1]-du[1,2]; f02 = du[3,1]-du[1,3]; f12 = du[3,2]-du[2,3];
  eh = -g/4 Sum[eta[[dd]] eta[[mu]] eta[[nu]] dh[[dd,mu,nu]]^2,
      {dd,4},{mu,4},{nu,4}]
    + g/2 Sum[eta[[nu]] Sum[eta[[mu]] dh[[mu,mu,nu]], {mu,4}]^2, {nu,4}];
  maxwell = g/2 (f01^2+f02^2-f12^2);
  (* Common factors of i and the overall Fourier sign do not change zeros.
     This also tests equality modulo a total derivative for the quadratic action. *)
  fourierRules = Flatten[Table[du[aa,dd] -> kk[dd] uu[aa], {aa,3},{dd,3}]];
  AssertZero["tensor/Einstein-Maxwell-modulo-boundary", (eh-maxwell) /. fourierRules];

  Do[
    label = If[normal==1, "timelike", "spacelike"];
    xx = eta[[normal]] p^2; kt = g-xx aone;
    (* delta(phi_mn) = -delta(Gamma^rho_mn) phi_rho, with phi_normal=p. *)
    deltaH = Table[-eta[[normal]] p/2 (dh[[mu,nu,normal]]
      + dh[[nu,mu,normal]] - dh[[normal,mu,nu]]), {mu,4},{nu,4}];
    AssertZero["tensor/"<>label<>"/trace", Sum[eta[[mu]] deltaH[[mu,mu]], {mu,4}]];
    AssertZero["tensor/"<>label<>"/gradient-contraction", eta[[normal]] p deltaH[[normal]]];
    Do[AssertZero["tensor/"<>label<>"/compact-Hessian-"<>ToString[aa],
      deltaH[[aa,4]] + eta[[normal]] p/2
        (dh[[aa,normal,4]]-dh[[normal,aa,4]])], {aa,3}];
    (* L2-L5 have no derivative-quadratic term in this odd polarization:
       the trace and gradient contractions just tested vanish. Variations
       of couplings and inverse metrics carry fewer derivatives. *)
    deltaL1 = Sum[eta[[mu]] eta[[nu]] deltaH[[mu,nu]]^2, {mu,4},{nu,4}];
    action = Expand[maxwell+aone deltaL1];
    expected = If[normal==1, (kt (f01^2+f02^2)-g f12^2)/2,
      (kt f01^2+g f02^2-kt f12^2)/2];
    AssertZero["tensor/"<>label<>"/parent-action", action-expected];
    symbol = action /. fourierRules;
    matrix = Table[D[symbol,uu[aa],uu[bb]], {aa,3},{bb,3}];
    cone = If[normal==1, -kt kk[1]^2+g (kk[2]^2+kk[3]^2),
      -g kk[1]^2+kt kk[2]^2+g kk[3]^2];
    AssertZero["tensor/"<>label<>"/gauge-null", matrix.Array[kk,3]];
    adj = Table[(-1)^(aa+bb) Det[Drop[matrix,{bb},{aa}]], {aa,3},{bb,3}];
    AssertZero["tensor/"<>label<>"/characteristic-adjugate",
      adj+kt cone Outer[Times,Array[kk,3],Array[kk,3]]];
    AssociateTo[results, label -> <|"X"->xx,"KT"->kt,"ParentAction"->action,
      "GaugeSingularSymbol"->matrix,"PhysicalCharacteristic"->cone|>],
    {normal,{1,2}}];
  results
]
