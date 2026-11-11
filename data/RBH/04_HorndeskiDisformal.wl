(* Letter Eqs. (12)-(14), SM S4-S6 and the comparison in S2.
   Horndeski convention chi=-X/2. The full Horndeski action is reduced
   from the metric and varied term by term (G2, G3, G4, G5) first.
   A second route uses standard background equations and
   odd coefficients below are literature inputs (2204.13837v2, Eqs. 5--8
   and its odd coefficients). Their off-shell combination, coordinate
   transformation and turning-point limit are computed independently.
   The map identities are exact matrix calculations for arbitrary local
   non-null scalar gradient and arbitrary projected gradient. The
   extension to parity-even differential invariants uses the parity and
   smooth-functional-inverse hypotheses stated in the manuscript.
   No external Mathematica package is required. *)

Module[{r,b,bp,p,p2,ta,da,g2,g3x,g3p,g4,g4x,g4xx,g4p,g4px,g4pp,
  g5x,g5xx,g5p,g5px,g5pp,a1,a2,a3,a4,a5,a6,a7,a8,a10,
  ett,eang,calf,calh,calg,chip,chain,hprime,sp,wp,xx,j,j1,j2,
  br,bpr,pr,qr,tr,fr,hr,gr,hx,t,aa,bb,clock,conf,q,barq,clockD,
  eta,normal,phi,u,dx,kappa,z,kin,lam,mu,gbar,invbar,invExpected,
  upperPhi,upperU,xbar,projector,ubar,zbar,back,composed,direct,
  h,dgi,deltaX,deltaP,oddMap,mat,principal,symbol,vec,detM,
  w,om,omega,curv,grad2,lap,boxLog,weylJac,rawEH,reducedEH,
  nu,eps,nn,jj,lapse,scalar,mult,kinetic,constraint,shapeConstraint,
  gradR0,gradR1,gradR2,gradRSeries,invisible,cc,gamma,barX,
  tensorMix,mix,mb,mwh,ft4,results=<||>,maps=<||>,i,k,assumptions,
  actionDerivation, properTimeCurrent, radialRate, hTensor, fTensor},

  (* Direct action route. Compute the Ricci tensor from the connection,
     form each Horndeski invariant and take the Euler shape combination.
     A=sign N^2 and B=sign Z^2 with positive N,Z,J avoid square-root
     branch assumptions. No scalar equation or division by phi', J'
     enters this derivation. G3 uses the convention -G3 Box(phi). *)
  actionDerivation=Module[{x,nf,jf,zf,pf,g2f,g3f,g4f,g5f,sign,bc,s,measure,
    chi,lambda,box,ricci,scalarR,einstein,lags,hs,fs,euler,idx,off,current,
    u,c4,c5,tt,theta,az,coords,metric,inverse,connection,ricciDirect,
    ri,a,b,c,patches=<||>,terms,labels={"G2","G3","G4","G5"}},
    coords={tt,x,theta,az};
    euler[lag_,field_]:=D[lag,field[x]]-D[D[lag,field'[x]],x]
      +D[D[lag,field''[x]],{x,2}];
    Do[
      bc=sign zf[x]^2; s=nf[x]/zf[x]; measure=s jf[x]^2;
      chi=-bc pf'[x]^2/2;
      metric=DiagonalMatrix[{-sign nf[x]^2,1/bc,jf[x]^2,jf[x]^2 Sin[theta]^2}];
      inverse=Inverse[metric];
      connection=Table[Cancel[Sum[inverse[[a,c]]/2
        (D[metric[[c,b]],coords[[ri]]]+D[metric[[c,ri]],coords[[b]]]
          -D[metric[[ri,b]],coords[[c]]]),{c,4}]],{a,4},{ri,4},{b,4}];
      ricciDirect=Table[FullSimplify[inverse[[ri,ri]](
        Sum[D[connection[[a,ri,ri]],coords[[a]]]
          -D[connection[[a,ri,a]],coords[[ri]]],{a,4}]
        +Sum[connection[[a,a,b]]connection[[b,ri,ri]]
          -connection[[a,ri,b]]connection[[b,ri,a]],{a,4},{b,4}])],{ri,4}];
      lambda={bc nf'[x]pf'[x]/nf[x],bc pf''[x]+D[bc,x]pf'[x]/2,
        bc jf'[x]pf'[x]/jf[x],bc jf'[x]pf'[x]/jf[x]};
      box=Total[lambda];
      ricci={-bc nf''[x]/nf[x]-(D[bc,x]/2+2 bc jf'[x]/jf[x])nf'[x]/nf[x],
        -bc nf''[x]/nf[x]-D[bc,x]nf'[x]/(2 nf[x])-2 bc jf''[x]/jf[x]
          -D[bc,x]jf'[x]/jf[x],
        1/jf[x]^2-bc jf'[x]^2/jf[x]^2-bc jf''[x]/jf[x]
          -(bc nf'[x]/nf[x]+D[bc,x]/2)jf'[x]/jf[x],
        1/jf[x]^2-bc jf'[x]^2/jf[x]^2-bc jf''[x]/jf[x]
          -(bc nf'[x]/nf[x]+D[bc,x]/2)jf'[x]/jf[x]};
      AssertZero["Horndeski/action-Ricci-from-metric/sign="<>ToString[sign],ricciDirect-ricci];
      scalarR=Total[ricciDirect]; einstein=ricciDirect-scalarR/2;
      c4=g4f[pf[x],chi]; c5=g5f[pf[x],chi];
      lags={g2f[pf[x],chi],-g3f[pf[x],chi]box,
        c4 scalarR+Derivative[0,1][g4f][pf[x],chi](box^2-Total[lambda^2]),
        c5 einstein.lambda-Derivative[0,1][g5f][pf[x],chi]
          (box^3-3 box Total[lambda^2]+2 Total[lambda^3])/6};
      hs={0,0,2 c4+2 bc pf'[x]^2 Derivative[0,1][g4f][pf[x],chi],
        -bc pf'[x]^2 Derivative[1,0][g5f][pf[x],chi]
          -bc^2 jf'[x]/jf[x]pf'[x]^3 Derivative[0,1][g5f][pf[x],chi]};
      fs={0,0,2 c4,bc pf'[x]^2 Derivative[1,0][g5f][pf[x],chi]
        -bc D[bc,x]pf'[x]^3 Derivative[0,1][g5f][pf[x],chi]/2
        -bc^2 pf'[x]^2 pf''[x]Derivative[0,1][g5f][pf[x],chi]};
      u=sign nf[x]^2/jf[x]^2; terms=<||>;
      Do[
        off=2 nf[x]euler[measure lags[[idx]],nf]-jf[x]euler[measure lags[[idx]],jf];
        current=jf[x]^4 hs[[idx]] D[u,x]/s;
        AssertZero["Horndeski/action-off-shell/"<>labels[[idx]]<>"/sign="<>ToString[sign],
          off-D[current,x]-2 s fs[[idx]]];
        AssociateTo[terms,labels[[idx]]-><|"H"->hs[[idx]],"F"->fs[[idx]],
          "Current"->current|>],{idx,4}];
      AssociateTo[patches,ToString[sign]-><|"Terms"->terms,
        "RicciEigenvalues"->ricciDirect,"ScalarHessianEigenvalues"->lambda|>],
      {sign,{1,-1}}];
    <|"SignedPatches"->patches,
      "OffShellIdentity"->"2 N E_N-J E_J = W' + 2 s F; W=J^4 H U'/s",
      "Couplings"->"Arbitrary G2(phi,chi), G3(phi,chi), G4(phi,chi), G5(phi,chi)"|>
  ];
  AssociateTo[results,"HorndeskiActionDerivation"->actionDerivation];

  (* Independent jets: ta=A'/A, da=A''/A, p=phi', p2=phi'', bp=B'.
     Each Gi partial derivative is an independent symbol. The chain
     rule retains both phi and chi dependence. *)
  a1=-b^2 (g3x-2 g4px) p^2-2 g4p b;
  a2=2 b^3 (2 g4xx-g5px) p^3-4 b^2 (g4x-g5p) p;
  a3=-b^4 g5xx p^4+b^2 g5x (3 b-1) p^2;
  a4=b^2 (2 g4xx-g5px) p^4+b (3 g5p-4 g4x) p^2-2 g4;
  a5=-(g5xx b^3 p^5-b g5x (5 b-1) p^3)/2;
  a6=b (g3p-2 g4pp) p^2+g2;
  a7=-2 b^2 (2 g4px-g5pp) p^3-4 g4p b p;
  a8=g5px b^3 p^4-b (2 g4x b-g5p b-g5p) p^2-2 g4 (b-1);
  a10=g5px b^3 p^4/2-b^2 (2 g4x-g5p) p^2/2-g4 b;
  ett=(a1+a2/r+a3/r^2) p2+(p a1/(2 b)+a4/r+a5/r^2) bp
    +a6+a7/r+a8/r^2;
  eang=((a2+((2 b-1) p a3+2 b a5)/(b p r)) ta/4+a1+a2/(2 r)) p2
    +(2 b a4-p a2+(2 b a5-p a3)/r) (da-ta^2/2)/4
    +(a4+(2 b (2 b+1) a5-p a3)/(2 b^2 r)) ta bp/4
    +(a7/4+a10/r) ta+(p a1/b+a4/r) bp/2+a6+a7/(2 r);
  calf=2 g4+b p^2 g5p-b bp p^3 g5x/2-b^2 p^2 p2 g5x;
  calh=2 g4+b p^2 (2 g4x-g5p)-b^2 p^3 g5x/r;
  calg=2 g4+b p^2 (2 g4x-g5p)-b^2 ta p^3 g5x/2;
  chip=-bp p^2/2-b p p2;
  chain={r->1,b->bp,p->p2,g4->g4p p+g4x chip,
    g4x->g4px p+g4xx chip,g5p->g5pp p+g5px chip,
    g5x->g5px p+g5xx chip};
  hprime=Total[(D[calh,First[#]] Last[#])& /@ chain];
  sp=ta/2-bp/(2 b);
  (* W=r H(r A'-2 A)/s; wp denotes W'/s. *)
  wp=b ((calh+r hprime-r calh sp) (r ta-2)+r calh (r da-ta));
  AssertZero["Horndeski/full-off-shell-background-identity",
    wp+2 calf-2 r^2 (ett-eang)];
  AssertZero["Horndeski/G-minus-H-factorization",
    calg-calh+b^2 p^3 g5x (ta-2/r)/2];
  AssertZero["Horndeski/turning-G-equals-H",(calg-calh)/.ta->2/r];
  (* Scalar stalls are covered when the covariant coupling jets are finite.
     No continuation of a singular scalar-clock representation is assumed. *)
  AssertZero["Horndeski/smooth-scalar-stall-coefficients",
    ({calf,calg,calh}/.p->0)-{2 g4,2 g4,2 g4}];
  AssertZero["Horndeski/stationary-background-equation",
    (2 r^2 (ett-eang)-b r^2 calh (da-2/r^2)-2 calf)/.ta->2/r];
  AssertZero["Horndeski/G2-G3-shape-cancellation",
    Table[D[ett-eang,xx],{xx,{g2,g3p,g3x}}]];

  (* Transform from areal radius r to J(x). Apparent J' denominators
     cancel in all three coefficients, so the final formula is also
     finite at a smooth J'=0 point. *)
  br=b j1^2; bpr=bp j1+2 b j2; pr=p/j1;
  qr=p2/j1^2-p j2/j1^3; tr=ta/j1;
  fr=2 g4+br pr^2 g5p-br bpr pr^3 g5x/2-br^2 pr^2 qr g5x;
  hr=2 g4+br pr^2 (2 g4x-g5p)-br^2 pr^3 g5x/j;
  gr=2 g4+br pr^2 (2 g4x-g5p)-br^2 tr pr^3 g5x/2;
  hx=2 g4+b p^2 (2 g4x-g5p)-b^2 p^3 j1 g5x/j;
  AssertZero["Horndeski/general-radial-F",fr-calf];
  AssertZero["Horndeski/general-radial-H",hr-hx];
  AssertZero["Horndeski/general-radial-G",gr-calg];
  AssertZero["Horndeski/general-radial-turn",(calg-hx)/.ta->2 j1/j];
  AssertTrue["Horndeski/no-Jprime-pole",
    And@@(FreeQ[Denominator[Cancel[#]],j1]& /@ {fr,hr,gr})];
  (* Check the derivative expression for U=A/J^2 directly. *)
  AssertZero["Horndeski/U-stationary-logarithmic-condition",
    D[aa[t]/bb[t]^2,t]-aa[t]/bb[t]^2 (aa'[t]/aa[t]-2 bb'[t]/bb[t])];
  AssociateTo[results,"Horndeski"-><|"F"->calf,"H"->calh,"G"->calg,
    "OffShellResidual"->Cancel[wp+2 calf-2 r^2 (ett-eang)],
    "BackgroundEquationsSource"->"2204.13837v2, Eqs. (5)-(8)",
    "Identity"->"W'/s+2 F=2 r^2(Ett-Eangular); W=J^4 H U'/s"|>];

  (* The static-exterior speed ratios G/F and G/H are inputs from
     arXiv:2204.13837, Eq. (2.16). Write their characteristic in one
     coordinate basis and evaluate it in the two orthonormal frames.
     kAngular2 is the squared physical angular wave number in either frame.
     The ratios and the overall normalization here assume nonzero F,G,H. *)
  AssociateTo[results, "HorndeskiPropagation" -> Module[
    {fTensorCoeff, gTensorCoeff, hTensorCoeff, metricA, metricB,
     ktCoord2, kxCoord2, kAngular2, kLong2, omega2, coordinateCone,
     exterior, interior, speeds, exteriorSpeeds, interiorSpeeds, coefficients},
    coordinateCone = -fTensorCoeff hTensorCoeff ktCoord2/metricA
      + gTensorCoeff hTensorCoeff metricB kxCoord2 + fTensorCoeff gTensorCoeff kAngular2;
    exterior = coordinateCone /. {ktCoord2 -> metricA omega2, kxCoord2 -> kLong2/metricB};
    interior = coordinateCone /. {ktCoord2 -> -metricA kLong2, kxCoord2 -> -omega2/metricB};
    speeds[cone_] := -{Coefficient[cone, kLong2], Coefficient[cone, kAngular2]}/Coefficient[cone, omega2];
    exteriorSpeeds = speeds[exterior];
    interiorSpeeds = speeds[interior];
    AssertZero["Horndeski/exterior-radial-and-angular-speeds",
      exteriorSpeeds - {gTensorCoeff/fTensorCoeff, gTensorCoeff/hTensorCoeff}];
    AssertZero["Horndeski/interior-longitudinal-and-angular-speeds",
      interiorSpeeds - {fTensorCoeff/gTensorCoeff, fTensorCoeff/hTensorCoeff}];
    coefficients = {fTensorCoeff -> calf, gTensorCoeff -> calg, hTensorCoeff -> hx};
    AssertZero["Horndeski/interior-isotropic-turning-symbol",
      (interior/hTensorCoeff + hTensorCoeff omega2 - fTensorCoeff (kLong2 + kAngular2))
        /. coefficients /. ta -> 2 j1/j];
    <|"ExteriorRadialAndAngularSpeedsSquared" -> (exteriorSpeeds /. coefficients),
      "InteriorLongitudinalAndAngularSpeedsSquared" -> (interiorSpeeds /. coefficients),
      "InteriorTurningSpeedSquared" -> calf/hx|>
  ]];

  (* Eq. (S27) inside the horizons: xdot=radialRate>0, s=a/xdot,
     U=-q^2 and d/dx=(1/xdot)d/dt, where t here denotes proper time.
     Convert the on-shell current equation before setting qdot=0. *)
  q=aa[t]/bb[t];
  properTimeCurrent=bb[t]^4 hTensor[t] (D[-q^2,t]/radialRate[t])/(aa[t]/radialRate[t]);
  AssertZero["Horndeski/proper-time-turning-identity",
    (-(D[properTimeCurrent,t]+2 aa[t] fTensor[t])/(2 aa[t] bb[t]^2)
      -hTensor[t] D[q,{t,2}]/q+fTensor[t]/bb[t]^2) /. aa'[t]->aa[t] bb'[t]/bb[t]];

  (* Shape and clock transformation. C=conf^2, E=clock^2 use the
     positive square-root branch; clock>0 is the time orientation. *)
  barq=(conf[t] aa[t])/(conf[t] bb[t]);
  clockD[expr_]:=D[expr,t]/clock[t];
  AssertZero["transport/shape-invariance",barq-q];
  AssertZero["transport/clock-second-derivative",
    clockD[clockD[barq]]-D[q,{t,2}]/clock[t]^2
      +clock'[t] D[q,t]/clock[t]^3];
  AssertZero["transport/stationary-acceleration",
    (clockD[clockD[barq]]-D[q,{t,2}]/clock[t]^2)/.
      Derivative[1][aa][t]->aa[t] bb'[t]/bb[t]];
  AssertTrue["transport/positive-clock-preserves-acceleration-sign",
    FullSimplify[Sign[nu/cc^2]==Sign[nu],Assumptions->Element[nu,Reals]&&cc>0]];

  (* Congruence identities for a generic symmetric principal matrix and
     a generic invertible upper-triangular field map. Sylvester's law
     supplies inertia invariance for every real invertible map. *)
  mat={{cc[1],cc[4],cc[5]},{0,cc[2],cc[6]},{0,0,cc[3]}};
  principal={{symbol[1],symbol[4],symbol[5]},
    {symbol[4],symbol[2],symbol[6]},{symbol[5],symbol[6],symbol[3]}};
  vec=Array[nu,3]; detM=Det[mat];
  AssertZero["transport/principal-quadratic-form",
    vec.(Transpose[mat].principal.mat).vec-(mat.vec).principal.(mat.vec)];
  AssertZero["transport/principal-determinant",
    Det[Transpose[mat].principal.mat]-detM^2 Det[principal]];
  AssertZero["transport/principal-inverse-congruence",
    Transpose[Inverse[mat]].(Transpose[mat].principal.mat).Inverse[mat]-principal];

  (* Explicit T_lambda, using both causal types of non-null phi. No
     restriction is placed on the norm or direction of the projected u.
     dx=u+kappa phi is a generic covector, with u perpendicular to phi. *)
  eta=DiagonalMatrix[{-1,1,1,1}];
  Do[
    phi=p UnitVector[4,normal];
    u=Table[If[i==normal,0,nu[i]],{i,4}];
    dx=u+kappa phi;
    kin=phi.eta.phi; z=u.eta.u;
    upperPhi=eta.phi; upperU=eta.u;
    AssertZero["map/orthogonality/normal="<>ToString[normal],phi.eta.u];
    gbar=eta+lam^2 z Outer[Times,phi,phi]
      +lam (Outer[Times,phi,u]+Outer[Times,u,phi]);
    invExpected=eta-lam (Outer[Times,upperPhi,upperU]+Outer[Times,upperU,upperPhi])
      +lam^2 kin Outer[Times,upperU,upperU];
    (* A multiplication proof verifies the inverse without asking
       Inverse to expand a large generic rational matrix. *)
    AssertZero["map/inverse-metric/normal="<>ToString[normal],gbar.invExpected-IdentityMatrix[4]];
    AssertZero["map/determinant/normal="<>ToString[normal],Det[gbar]/Det[eta]-1];
    invbar=invExpected;
    xbar=Cancel[phi.invbar.phi];
    AssertZero["map/X-invariant/normal="<>ToString[normal],xbar-kin];
    projector=IdentityMatrix[4]-Outer[Times,phi,invbar.phi]/xbar;
    ubar=Cancel /@ (projector.dx);
    AssertZero["map/projected-gradient/normal="<>ToString[normal],ubar-u-lam z phi];
    zbar=Cancel[ubar.invbar.ubar];
    AssertZero["map/Z-invariant/normal="<>ToString[normal],zbar-z];
    back=gbar+lam^2 zbar Outer[Times,phi,phi]
      -lam (Outer[Times,phi,ubar]+Outer[Times,ubar,phi]);
    AssertZero["map/functional-inverse/normal="<>ToString[normal],back-eta];
    composed=gbar+mu^2 zbar Outer[Times,phi,phi]
      +mu (Outer[Times,phi,ubar]+Outer[Times,ubar,phi]);
    direct=eta+(lam+mu)^2 z Outer[Times,phi,phi]
      +(lam+mu) (Outer[Times,phi,u]+Outer[Times,u,phi]);
    AssertZero["map/composition-group/normal="<>ToString[normal],composed-direct];
    AssociateTo[maps,ToString[normal]-><|"Metric"->gbar,"InverseMetric"->invbar,
      "X"->kin,"ProjectedNorm"->z,"ProjectedGradientAfterMap"->ubar|>],
    {normal,{1,2}}];
  AssociateTo[results,"ExplicitMap"->maps];

  (* On a radial background dx is parallel to phi. Compute odd inverse
     metric and projector variations to verify that the T_lambda map
     acts identically on the odd sector. *)
  h={{0,0,0,nu[1]},{0,0,0,nu[2]},{0,0,0,nu[3]},
    {nu[1],nu[2],nu[3],0}};
  Do[
    phi=p UnitVector[4,normal]; dx=kappa phi; kin=phi.eta.phi;
    dgi=-eta.h.eta; deltaX=phi.dgi.phi;
    projector=IdentityMatrix[4]-Outer[Times,phi,eta.phi]/kin;
    deltaP=-Outer[Times,phi,dgi.phi]/kin
      +deltaX Outer[Times,phi,eta.phi]/kin^2;
    AssertZero["odd-map/radial-projection/normal="<>ToString[normal],projector.dx];
    AssertZero["odd-map/delta-X/normal="<>ToString[normal],deltaX];
    AssertZero["odd-map/delta-projected-gradient/normal="<>ToString[normal],deltaP.dx];
    oddMap=h+lam (Outer[Times,phi,deltaP.dx]+Outer[Times,deltaP.dx,phi]);
    AssertZero["odd-map/identity/normal="<>ToString[normal],oddMap-h],
    {normal,{1,2}}];

  (* Weyl conformal map. The four-dimensional conformal laws are inputs:
     Wbar=W/Omega^2 and Rbar=Omega^-1[R-3 Box log Omega
       -(3/2)(grad log Omega)^2]. Everything below follows by calculus. *)
  weylJac=D[w/om[w]^2,w];
  AssertZero["Weyl/kinetic-inverse-Jacobian",
    weylJac-(om[w]-2 w om'[w])/om[w]^3];
  AssertTrue["Weyl/linear-Omega-admissible-domain",
    FullSimplify[1+lam w>0 && 1-lam w!=0,
      Assumptions->Element[{lam,w},Reals]&&-1<lam w<1]];
  boxLog=lap/omega-grad2/omega^2;
  rawEH=Expand[omega (curv-3 boxLog-3 grad2/(2 omega^2))];
  reducedEH=omega curv+3 grad2/(2 omega);
  AssertZero["Weyl/EH-pullback-up-to-divergence",rawEH-reducedEH+3 lap];
  AssociateTo[results,"Weyl"-><|"Jacobian"->weylJac,
    "EinsteinHilbertDensityModuloBoundary"->reducedEH,
    "OmittedDivergence"->-3 lap|>];

  (* Mimetic term: the exact volume-preserving rescaling leaves X
     unchanged because B and phi are fixed. *)
  kinetic=b p^2;
  constraint=nn jj^2 mult (kinetic-scalar);
  shapeConstraint=constraint/.{nn->Exp[2 eps] nn,jj->Exp[-eps] jj};
  AssertZero["extension/mimetic-shape-invariance",FullSimplify[shapeConstraint-constraint]];
  (* If grad R at order zero and one vanish, a squared-gradient operator
     begins at fourth perturbative order, including metric variations. *)
  gradRSeries=eps^2 Array[gradR2,3]/2;
  invisible=(1+eps cc) gradRSeries.(IdentityMatrix[3]+eps mat).gradRSeries;
  AssertZero["extension/curvature-gradient-no-odd-quadratic",
    D[invisible,{eps,2}]/.eps->0];

  (* Pure-disformal kinetic map and the mixing-zero Jacobian. The
     model-specific Gamma_X relation is an input of the cited
     cosmological construction, not a universal identity. *)
  barX=xx/(1+gamma[xx] xx);
  AssertZero["mixing/pure-disformal-kinetic-Jacobian",
    D[barX,xx]-(1-xx^2 gamma'[xx])/(1+gamma[xx] xx)^2];
  (* arXiv:2204.05889, Eq. (23): Gamma_X=-(Mb-GT)/(X^2 GT).
     Check the Eq. (24) identification from this cited input. *)
  AssertZero["mixing/model-specific-Jacobian-numerator",
    (1-xx^2 gamma'[xx]/.gamma'[xx]->-(mb-tensorMix)/(xx^2 tensorMix))-mb/tensorMix];
  mwh=tensorMix-2 ft4 b^2 p^4;
  AssertZero["mixing/wormhole-Horndeski-limit",(mwh/.ft4->0)-tensorMix];
  AssociateTo[results,"PureDisformal"-><|"Xbar"->barX,
    "KineticJacobian"->D[barX,xx],
    "MixingIdentificationSource"->"2204.05889, pure-disformal cosmological construction"|>];
  AssociateTo[results,"Scope"->"Exact local algebra. General parity selection, Sylvester inertia invariance, and the existence of the forced maximum are mathematical arguments in the manuscript."];
  results
]
