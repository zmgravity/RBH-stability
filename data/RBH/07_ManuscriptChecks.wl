(* Additional local derivations for the November 2026 manuscript.
   The compact-interval maximum argument, continuity arguments and general
   functional-inverse/parity hypotheses remain analytical arguments.
   Smooth finite covariant couplings are required at a scalar stall.
   The U-DHOST coefficient conditions are inputs from arXiv:2110.03194v2,
   Eqs. (4)-(6); the tests below evaluate them for the displayed model. *)
Module[{t,z,theta,az,a,b,ph,coords,metric,inv,connection,ricci,scalarR,
  einstein,mu,nu,rho,lam,q,vol,gradient,stress,pressure,px,shape,
  x,cf,ca1,mpl,xi,ax,unitaryResidual,udRules,lx,udCorrection,
  eta,normal,p,phi,upperPhi,c,d,gbar,invbar,kin,h,dgi,deltaX,oddMap,
  conf,disf,barX,jac,clock,space,kt,ft,w,k,omega2,barCone,speed,
  re,a2,a3,b2,b3,kp,k1,fp,f1,order,aa,bb,ss,rr,current,
  currentFactor,scaledDerivative,leading,rule,finiteSpeed,concrete,
  denominator,numerator,arctanA1,arctanF,arctanK,kinetic,n},

  (* S1 A, Eqs. (S3)-(S4): evaluate the Einstein tensor and scalar stress
     in the proper-time Kantowski--Sachs metric. *)
  coords={t,z,theta,az};
  metric=DiagonalMatrix[{-1,a[t]^2,b[t]^2,b[t]^2 Sin[theta]^2}];
  inv=Inverse[metric];
  connection=Table[Cancel[Sum[inv[[lam,rho]]/2
    (D[metric[[rho,nu]],coords[[mu]]]+D[metric[[rho,mu]],coords[[nu]]]
      -D[metric[[mu,nu]],coords[[rho]]]),{rho,4}]],{lam,4},{mu,4},{nu,4}];
  ricci=Table[FullSimplify[Sum[D[connection[[lam,mu,nu]],coords[[lam]]]
      -D[connection[[lam,mu,lam]],coords[[nu]]],{lam,4}]
    +Sum[connection[[lam,lam,rho]]connection[[rho,mu,nu]]
      -connection[[lam,nu,rho]]connection[[rho,mu,lam]],{lam,4},{rho,4}]],
    {mu,4},{nu,4}];
  scalarR=FullSimplify[Tr[inv.ricci]];
  einstein=FullSimplify[inv.ricci-IdentityMatrix[4] scalarR/2];
  AssertZero["kessence/Einstein-longitudinal-from-metric",
    einstein[[2,2]]+2 b''[t]/b[t]+b'[t]^2/b[t]^2+1/b[t]^2];
  AssertZero["kessence/Einstein-angular-from-metric",
    einstein[[3,3]]+a''[t]/a[t]+b''[t]/b[t]+a'[t] b'[t]/(a[t] b[t])];
  AssertZero["kessence/Einstein-angular-equality",einstein[[3,3]]-einstein[[4,4]]];
  gradient=D[ph[t],#]& /@ coords;
  stress=inv.(pressure metric-2 px Outer[Times,gradient,gradient]);
  AssertZero["kessence/proper-time-gradient-index-sign",
    {gradient, inv.gradient} - {{ph'[t],0,0,0}, {-ph'[t],0,0,0}}];
  AssertZero["kessence/mixed-temporal-stress",
    stress[[1,1]] - (pressure-2 px (gradient.inv.gradient))];
  AssertZero["kessence/isotropic-spatial-stress",
    {stress[[2,2]],stress[[3,3]],stress[[4,4]]}-{pressure,pressure,pressure}];
  q=a[t]/b[t]; vol=a[t] b[t]^2;
  shape=D[vol D[q,t]/q,t]/vol-1/b[t]^2;
  AssertZero["kessence/Einstein-difference-shape-current",
    einstein[[3,3]]-einstein[[2,2]]+shape];
  AssertZero["kessence/stationary-acceleration",
    (shape-D[q,{t,2}]/q+1/b[t]^2)/.a'[t]->a[t] b'[t]/b[t]];

  (* The covariant quadratic coefficients at X=0 need finite A1 and F2.
     This calculation supplies no continuation of a singular unitary clock. *)
  kinetic=cf-x ca1;
  AssertZero["quadratic/smooth-scalar-stall-KT",(kinetic/.x->0)-cf];
  AssertZero["quadratic/smooth-nondegenerate-scalar-stall-speed",
    (cf/kinetic/.x->0)-1];

  (* Explicit U-DHOST example, including its unitary degeneracy equation. *)
  udRules={cf->mpl^2/2,ax[1]->0,ax[2]->0,ax[3]->0,ax[4]->-xi x,ax[5]->xi};
  unitaryResidual=4 (x (ax[1]+3 ax[2])+2 cf)
    (ax[1]+ax[2]+x (ax[3]+ax[4])+x^2 ax[5])-3 x (2 ax[2]+x ax[3])^2;
  AssertZero["UDHOST/example-unitary-degeneracy",unitaryResidual/.udRules];
  udCorrection=(-xi x) (x lx^2)+xi (x^2 lx^2);
  AssertZero["UDHOST/example-homogeneous-deformation",udCorrection];
  AssertZero["UDHOST/example-tensor-coefficients",
    ({cf-x ax[1],cf}/.udRules)-{mpl^2/2,mpl^2/2}];
  AssertTrue["UDHOST/example-noncovariant-degeneracy-for-nonzero-xi",
    FullSimplify[(-xi x!=0),Assumptions->Element[{xi,x},Reals]&&xi!=0&&x<0]];

  (* S5 B: ordinary disformal inverse and full C(X), D(X) Jacobian. *)
  eta=DiagonalMatrix[{-1,1,1,1}];
  h={{0,0,0,ax[1]},{0,0,0,ax[2]},{0,0,0,ax[3]},
    {ax[1],ax[2],ax[3],0}};
  Do[
    phi=p UnitVector[4,normal]; upperPhi=eta.phi; kin=phi.eta.phi;
    gbar=c eta+d Outer[Times,phi,phi];
    invbar=eta/c-d Outer[Times,upperPhi,upperPhi]/(c (c+d kin));
    AssertZero["ordinary-map/inverse/normal="<>ToString[normal],gbar.invbar-IdentityMatrix[4]];
    AssertZero["ordinary-map/kinetic/normal="<>ToString[normal],phi.invbar.phi-kin/(c+d kin)];
    dgi=-eta.h.eta; deltaX=phi.dgi.phi;
    oddMap=c h+(conf'[x] eta+disf'[x] Outer[Times,phi,phi]) deltaX;
    AssertZero["ordinary-map/odd-algebraic/normal="<>ToString[normal],oddMap-c h],
    {normal,{1,2}}];
  barX=x/(conf[x]+disf[x] x); jac=D[barX,x];
  AssertZero["ordinary-map/general-kinetic-Jacobian",
    jac-(conf[x]-x conf'[x]-x^2 disf'[x])/(conf[x]+disf[x] x)^2];
  barCone=-kt w^2/clock+ft k^2/space;
  speed=(omega2/k^2)/.First[Solve[(barCone/.w^2->omega2)==0,omega2]];
  AssertZero["transport/orthonormal-speed-rescaling",speed-(clock/space) ft/kt];
  AssertZero["transport/full-turning-identity",
    clock (ft/kt)/(space b[t]^2)-speed/b[t]^2];

  (* S7 A: arbitrary finite order, with subleading metric/coefficient jets.
     Factor the current as rho^(order+1) times a smooth function. This
     proves the leading coefficient for every integer order >=0 without
     differentiating a truncated O(rho^n) expression as if it were exact. *)
  rr=re+rho; aa=a2+a3 rho; bb=b2+b3 rho; ss=Sqrt[aa/bb];
  current=rr^4 kp rho^order (1+k1 rho)/ss D[rho^2 aa/rr^2,rho];
  currentFactor=rr^4 kp (1+k1 rho)/ss (2 aa/rr^2+rho D[aa/rr^2,rho]);
  AssertZero["exterior/finite-order-current-factorization",
    FullSimplify[current-rho^(order+1) currentFactor,Assumptions->rho>0]];
  scaledDerivative=(order+1) currentFactor+rho D[currentFactor,rho];
  AssertZero["exterior/finite-order-current-derivative",
    FullSimplify[D[current,rho]/rho^order-scaledDerivative,Assumptions->rho>0]];
  leading=FullSimplify[(scaledDerivative+2 ss (fp+f1 rho))/.rho->0,
    Assumptions->a2>0&&b2>0&&re>0]/(2 Sqrt[a2/b2]);
  AssertZero["exterior/finite-order-leading-identity",
    FullSimplify[leading-(fp+(order+1) re^2 b2 kp),Assumptions->a2>0&&b2>0&&re>0]];
  rule=First[Solve[leading==0,fp]]; finiteSpeed=FullSimplify[kp/fp/.rule];
  AssertZero["exterior/finite-order-speed",finiteSpeed+1/((order+1) re^2 b2)];
  AssertTrue["exterior/finite-order-speed-negative",
    FullSimplify[finiteSpeed<0,Assumptions->re>0&&b2>0&&order>=0]];
  Do[
    concrete=FullSimplify[Limit[(D[current,rho]+2 ss rho^order (fp+f1 rho))/rho^order/.order->n,
      rho->0,Direction->"FromAbove"],Assumptions->re>0&&a2>0&&b2>0];
    AssertZero["exterior/finite-order-Taylor/n="<>ToString[n],
      FullSimplify[concrete-2 Sqrt[a2/b2] (fp+(n+1) re^2 b2 kp),
        Assumptions->re>0&&a2>0&&b2>0]],{n,{0,1,2,3}}];

  (* S7 C: derive K_T from the published arctangent A1 and F2, rather
     than treating its quoted tensor coefficient as a separate input. *)
  denominator=p Sin[Pi x]-3 Pi x;
  numerator=p^2 Sin[2 Pi x]-8 p Sin[Pi x]+6 Pi x;
  arctanF=numerator/denominator^2;
  arctanA1=2 p Sin[Pi x] (p Cos[Pi x]-3)/(x denominator^2);
  arctanK=Together[arctanF-x arctanA1];
  AssertZero["arctangent/KT-from-published-A1",TrigExpand[arctanK+2/denominator]];

  <|"KessenceEinsteinTensor"->einstein,"KessenceMixedStress"->stress,"KessenceShapeResidual"->shape,
    "OrdinaryDisformalJacobian"->jac,"TransportedSpeedSquared"->speed,
    "FiniteOrderLeadingEquation"->(leading==0),"FiniteOrderCoefficientRule"->rule,
    "FiniteOrderExtremalSpeedSquared"->finiteSpeed,
    "ArctangentA1"->arctanA1,
    "Scope"->{"Smooth finite couplings and a regular nondegenerate tensor sector at scalar stalls.",
      "The foliation derivation itself requires a regular scalar clock.",
      "General mimetic closure requires an off-shell seed shape identity; the added constraint alone does not prove closure for every on-shell-defined TPST action.",
      "Continuity of a near-extremal solution family and existence of an interior maximum are analytical arguments."}|>
]
