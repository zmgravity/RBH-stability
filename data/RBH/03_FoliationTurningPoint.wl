(* Letter Eqs. (6)-(10), SM S1 B-C, Eqs. (S5)-(S11).
   Derive the homogeneous shape equation with arbitrary lapse and arbitrary
   smooth L(eta,N,K,S,T,R,A). At zero shear, derive the odd two-derivative
   action from delta K_ij and the spatial Einstein curvature density.
   The latter is built from the perturbed three-metric; equality modulo a
   spatial boundary term is tested by its Euler derivatives.
   A regular scalar foliation, a,b,N>0, and frozen principal coefficients
   are used. The shift reduction requires nonzero k and KT; degenerate
   limits are checked separately in the undivided action/equation. *)

Module[{eta, a0, b0, al, be, lapse, lag, a, b, vol, q, kmat, trace,
  shear, sinv, tinv, rinv, args, reduced, ls, lt, lr, momentum, source,
  expectedMomentum, shapeEOM, turnRules, turnEOM, acceleration,
  qs, rad, kt, ft, acc, undivided, zeroK, flat,
  x, y, z, u, v, shift, eps, coords, h1, inverse1, conn1, conn2,
  ricci1, ricci2, scalar1, scalar2, curl, spatialTarget, spatialBoundary,
  euler2, f01, f02, f12, deltaK, deltaSigma, s2, t2,
  kbg, rbg, eta0, n0, generalized, quadraticCoefficient,
  du1, du2, d01, d02, d12, parentAction,
  k1, k2, kn, ulr, uli, utr, uti, vlr, vli, vtr, vti, u0r, u0i,
  wave, khat, ehat, avec, davec, u0, electric, magnetic, abs2,
  reality, fourierAction, splitAction, constraintRules, constrained,
  physical, omega, amplitude, characteristic, kineticHessian,
  i, j, m, ell, results},

  (* Homogeneous geometry, retaining the lapse until after variation. *)
  a = a0 Exp[al[eta] + 2 be[eta]];
  b = b0 Exp[al[eta] - be[eta]];
  vol = Simplify[a b^2]; q = Simplify[a/b];
  kmat = DiagonalMatrix[{D[a,eta]/a, D[b,eta]/b, D[b,eta]/b}]/lapse[eta];
  trace = Simplify[Tr[kmat]];
  shear = Simplify[kmat - IdentityMatrix[3] trace/3];
  sinv = Simplify[Tr[shear.shear]];
  tinv = Simplify[Tr[shear.shear.shear]];
  rinv = 2/b^2;
  AssertZero["foliation/volume-independent-of-shape", D[vol,be[eta]]];
  AssertZero["foliation/K", trace - 3 al'[eta]/lapse[eta]];
  AssertZero["foliation/S", sinv - 6 (be'[eta]/lapse[eta])^2];
  AssertZero["foliation/T", tinv - 6 (be'[eta]/lapse[eta])^3];
  AssertZero["foliation/curvature-shape-variation", D[rinv,be[eta]] - 2 rinv];
  AssertZero["foliation/log-q-proper-time", D[q,eta]/(lapse[eta] q) - 3 be'[eta]/lapse[eta]];

  args = {eta,lapse[eta],trace,sinv,tinv,rinv,0};
  reduced = lapse[eta] vol Apply[lag,args];
  ls = Apply[Derivative[0,0,0,1,0,0,0][lag],args];
  lt = Apply[Derivative[0,0,0,0,1,0,0][lag],args];
  lr = Apply[Derivative[0,0,0,0,0,1,0][lag],args];
  momentum = D[reduced,be'[eta]];
  source = D[reduced,be[eta]];
  expectedMomentum = vol (12 ls be'[eta]/lapse[eta] + 18 lt (be'[eta]/lapse[eta])^2);
  AssertZero["foliation/shape-momentum", Simplify[momentum-expectedMomentum]];
  AssertZero["foliation/shape-source", Simplify[source-2 lapse[eta] vol rinv lr]];
  shapeEOM = D[momentum,eta] - source;
  turnRules = {be'[eta]->0};
  turnEOM = Simplify[(shapeEOM/(lapse[eta] vol)) /. turnRules];
  AssertZero["foliation/zero-shear-equation",
    Simplify[turnEOM - ((12 ls be''[eta]/lapse[eta]^2-2 rinv lr) /. turnRules)]];
  acceleration = Simplify[(D[D[q,eta]/lapse[eta],eta]/(lapse[eta] q)) /. turnRules];
  AssertZero["foliation/q-acceleration-at-turn", acceleration-3 be''[eta]/lapse[eta]^2];
  AssertZero["foliation/undivided-turning-identity",
    Simplify[turnEOM/4 - ((ls acceleration-lr/b^2) /. turnRules)]];

  (* Keep the equation undivided when one coefficient vanishes. *)
  undivided = kt acc - ft/rad^2;
  zeroK = First[Solve[(undivided /. kt->0)==0,ft]];
  flat = First[Solve[(undivided /. acc->0)==0,ft]];
  AssertZero["foliation/zero-kinetic-forces-zero-gradient", ft /. zeroK];
  AssertZero["foliation/flat-turn-forces-zero-gradient", ft /. flat];
  AssertTrue["foliation/Lorentzian-cone-excludes-maximum",
    Reduce[rad>0 && kt!=0 && ft/kt>0 && undivided==0 && acc<=0,
      {rad,kt,ft,acc},Reals]===False];

  (* Local spatial metric h_ij=delta_ij+eps h1_ij. The odd direction is z,
     and partial_z=0. Compute both orders of the connection explicitly. *)
  coords = {x,y,z};
  h1 = {{0,0,u[x,y]},{0,0,v[x,y]},{u[x,y],v[x,y],0}};
  inverse1 = -h1;
  conn1 = Table[(D[h1[[ell,j]],coords[[i]]] + D[h1[[ell,i]],coords[[j]]]
    - D[h1[[i,j]],coords[[ell]]])/2,{ell,3},{i,3},{j,3}];
  conn2 = Table[Sum[inverse1[[ell,m]] conn1[[m,i,j]],{m,3}],{ell,3},{i,3},{j,3}];
  ricci1 = Table[Sum[D[conn1[[m,i,j]],coords[[m]]]-D[conn1[[m,i,m]],coords[[j]]],{m,3}],{i,3},{j,3}];
  ricci2 = Table[Sum[D[conn2[[m,i,j]],coords[[m]]]-D[conn2[[m,i,m]],coords[[j]]],{m,3}]
    + Sum[conn1[[m,m,ell]] conn1[[ell,i,j]]-conn1[[m,j,ell]] conn1[[ell,i,m]],{m,3},{ell,3}],{i,3},{j,3}];
  scalar1 = Simplify[Tr[ricci1]];
  scalar2 = Simplify[Tr[ricci2] + Sum[inverse1[[i,j]] ricci1[[i,j]],{i,3},{j,3}]];
  AssertZero["foliation/odd-linear-spatial-curvature", scalar1];
  AssertZero["foliation/odd-linear-volume", Tr[h1]];
  (* R0=0 in this frozen local calculation and delta R=0; sqrt(h) has no
     effect on the quadratic two-derivative term. Euler derivatives test
     equality of the derived density to -curl^2/2 modulo a boundary. *)
  curl = D[v[x,y],x]-D[u[x,y],y];
  spatialTarget = -curl^2/2;
  spatialBoundary = Expand[scalar2-spatialTarget];
  euler2[expr_,field_] := D[expr,field[x,y]]
    - D[D[expr,D[field[x,y],x]],x] - D[D[expr,D[field[x,y],y]],y]
    + D[D[expr,D[field[x,y],{x,2}]],{x,2}]
    + D[D[expr,D[field[x,y],x,y]],x,y]
    + D[D[expr,D[field[x,y],{y,2}]],{y,2}];
  AssertZero["foliation/spatial-curvature-boundary-u", Simplify[euler2[spatialBoundary,u]]];
  AssertZero["foliation/spatial-curvature-boundary-v", Simplify[euler2[spatialBoundary,v]]];

  (* Principal part of delta K: shift subtraction and dot h/2. *)
  f01=du1-D[shift[x,y],x]; f02=du2-D[shift[x,y],y]; f12=curl;
  deltaK={{0,0,f01},{0,0,f02},{f01,f02,0}}/2;
  deltaSigma=deltaK-IdentityMatrix[3] Tr[deltaK]/3;
  s2=Expand[Tr[deltaSigma.deltaSigma]];
  t2=Coefficient[Tr[MatrixPower[eps deltaSigma,3]],eps,2];
  AssertZero["foliation/odd-linear-K", Tr[deltaK]];
  AssertZero["foliation/odd-quadratic-S", s2-(f01^2+f02^2)/2];
  AssertZero["foliation/odd-quadratic-T", t2];
  (* At this principal order delta N=delta A=delta K=delta S=delta T=
     delta R=0. The arbitrary smooth L therefore contributes only its
     first derivatives L_S and L_R multiplying the quadratic invariants. *)
  generalized=lag[eta0,n0,kbg,eps^2 s2,0,rbg+eps^2 spatialTarget,0];
  quadraticCoefficient=Simplify[(D[generalized,{eps,2}]/2) /. eps->0];
  parentAction=quadraticCoefficient /. {
    Derivative[0,0,0,1,0,0,0][lag][eta0,n0,kbg,0,0,rbg,0]->kt,
    Derivative[0,0,0,0,0,1,0][lag][eta0,n0,kbg,0,0,rbg,0]->ft};
  AssertZero["foliation/general-smooth-L-odd-action",
    parentAction-(kt (f01^2+f02^2)-ft f12^2)/2];

  (* Complex Fourier amplitudes with real and imaginary parts independent.
     This checks the actual |...|^2 action, including the factor i. *)
  wave={k1,k2}; kn=Sqrt[k1^2+k2^2]; khat=wave/kn; ehat={-k2,k1}/kn;
  avec=khat (ulr+I uli)+ehat (utr+I uti);
  davec=khat (vlr+I vli)+ehat (vtr+I vti); u0=u0r+I u0i;
  electric=davec-I wave u0; magnetic=I (k1 avec[[2]]-k2 avec[[1]]);
  reality=Element[{k1,k2,ulr,uli,utr,uti,vlr,vli,vtr,vti,u0r,u0i,kt,ft},Reals]
    && k1^2+k2^2>0;
  abs2[expression_] := FullSimplify[ComplexExpand[Conjugate[expression] expression],Assumptions->reality];
  fourierAction=FullSimplify[kt Total[abs2 /@ electric]/2-ft abs2[magnetic]/2,Assumptions->reality];
  splitAction=kt (abs2[vlr+I vli-I kn u0]+vtr^2+vti^2)/2
    -ft kn^2 (utr^2+uti^2)/2;
  AssertZero["foliation/Fourier-longitudinal-transverse", FullSimplify[fourierAction-splitAction,Assumptions->reality]];
  AssertZero["foliation/shift-equation-real", Simplify[D[splitAction,u0r]-kt kn (kn u0r-vli)]];
  AssertZero["foliation/shift-equation-imaginary", Simplify[D[splitAction,u0i]-kt kn (kn u0i+vlr)]];
  constraintRules={u0r->vli/kn,u0i->-vlr/kn};
  constrained=Simplify[splitAction /. constraintRules];
  physical=kt (vtr^2+vti^2)/2-ft kn^2 (utr^2+uti^2)/2;
  AssertZero["foliation/physical-odd-action", constrained-physical];
  characteristic=kt omega^2-ft kn^2;
  AssertZero["foliation/physical-characteristic",
    D[physical,vtr,vtr] omega^2 + D[physical,utr,utr]-characteristic];
  kineticHessian=Table[D[physical,i,j],{i,{vtr,vti}},{j,{vtr,vti}}];
  AssertZero["foliation/KT-zero-kinetic-degeneracy", kineticHessian /. kt->0];
  AssertZero["foliation/FT-zero-spatial-degeneracy", D[physical,utr,utr] /. ft->0];

  <|"HomogeneousInvariants"-><|"K"->trace,"S"->sinv,"T"->tinv,"R"->rinv|>,
    "ShapeMomentum"->momentum,"ShapeSource"->source,
    "ZeroShearEquation"->turnEOM,"UndividedTurningRelation"->(undivided==0),
    "SpatialCurvatureQuadraticDensity"->scalar2,
    "SpatialCurvatureBoundaryTerm"->spatialBoundary,
    "OddParentAction"->parentAction,"FourierReducedAction"->physical,
    "PhysicalCharacteristic"->characteristic,
    "Assumptions"->{"Regular timelike scalar foliation; a,b,N positive.",
      "Principal coefficients frozen locally at zero shear.",
      "Shift elimination uses KT!=0 and k!=0; degenerate cases checked undivided.",
      "The compact-interval maximum argument is the analytical proof in the Letter."}|>
]
