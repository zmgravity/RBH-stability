(* Letter, extremal exterior test; Supplemental Material S7.
   Inputs: the BCCH Hayward couplings and the Baake et al. arctangent
   tensor coefficients cited in SOURCES.txt. Differentiate metrics at fixed
   mass and fixed couplings before imposing the horizon equations.
   Both examples are used only on their real radial-scalar exteriors. *)
Module[{rho, re, hh2, hh3, ff2, ff3, ke, k1, fe, f1, hh, ff, ss,
  curvatureResidual, leading, extremalRule, xx, r, mu, alpha, beta, gamma,
  h0, fCoupling, a1, a3, kt, profile, metric, current, massOfRadius,
  radiusE, massE, radius0, threshold, y, scale, coefficientRules,
  z, p, sigma, xz, dx, dxx, dz, dp, np, fc, kc, metricZ,
  trigRules, arctanCurrent, arctanResidual, logSlope, secondAtHorizon,
  trigAssumptions, rootRule, de, ne, metricAtExtremality, speed,
  horizonMassDimensionless, normalizedMetric, radiusRatio, massRatio,
  scalarSlope, xPositive, fNegative, t, sincNumerator, rootFunction},

  (* Independent smooth Taylor coefficients retain the potentially relevant
     first derivatives. The leading identity loses all of them. *)
  hh = rho^2 (hh2 + hh3 rho);
  ff = rho^2 (ff2 + ff3 rho);
  ss = Sqrt[(hh2 + hh3 rho)/(ff2 + ff3 rho)];
  curvatureResidual = D[(re+rho)^4 (ke+k1 rho)/ss D[hh/(re+rho)^2,rho],rho]
    + 2 ss (fe+f1 rho);
  leading = FullSimplify[(curvatureResidual /. rho->0)/(2 Sqrt[hh2/ff2]),
    Assumptions->re>0 && hh2>0 && ff2>0];
  AssertZero["exterior/extremal-leading-identity",leading-(fe+re^2 ff2 ke)];
  extremalRule = First[Solve[leading==0,fe]];
  AssertZero["exterior/extremal-speed",(ke/fe /. extremalRule)+1/(re^2 ff2)];
  AssertTrue["exterior/extremal-negative-speed",
    FullSimplify[(ke/fe /. extremalRule)<0,Assumptions->re>0 && ff2>0 && ke!=0]];
  AssertZero["exterior/zero-X-smooth-couplings",
    (leading /. ke->fe)-fe (1+re^2 ff2)];

  (* Hayward: derive K_T from the action dictionary, not from its target ratio. *)
  fCoupling = h0 alpha (5 alpha-6 beta xx)/xx^2;
  a1 = 6 h0 alpha (alpha-beta xx)/xx^3;
  kt = Cancel[fCoupling-xx a1];
  a3 = -4 D[fCoupling,xx]/xx+2 a1/xx-8 kt/(3 xx^2);
  AssertZero["Hayward/KT-from-couplings",kt+h0 alpha^2/xx^2];
  AssertZero["Hayward/outside-GLPV",xx a3+4 D[fCoupling,xx]-2 a1+8 kt/(3 xx)];
  profile = alpha r^3/(beta r^3+gamma);
  metric = 1-mu alpha r^2/(beta r^3+gamma);
  current = Cancel[r^4 (kt /. xx->profile) D[metric/r^2,r]];
  AssertZero["Hayward/background-identity",D[current,r]+2 (fCoupling /. xx->profile)];
  coefficientRules = {xx->alpha y/beta,h0->-scale/beta^2};
  AssertZero["Hayward/normalized-KT",(kt /. coefficientRules)/scale-1/y^2];
  AssertZero["Hayward/normalized-F2",(fCoupling /. coefficientRules)/scale-(6 y-5)/y^2];
  AssertZero["Hayward/speed-ratio",kt/fCoupling-1/(6 beta xx/alpha-5)];
  massOfRadius = mu /. First[Solve[metric==0,mu]];
  AssertZero["Hayward/horizon-mass-derivative",D[massOfRadius,r]-(beta-2 gamma/r^3)/alpha];
  radiusE=(2 gamma/beta)^(1/3); radius0=(5 gamma/beta)^(1/3);
  massE=FullSimplify[massOfRadius /. r->radiusE,Assumptions->alpha>0 && beta>0 && gamma>0];
  AssertZero["Hayward/extremal-mass",FullSimplify[massE-3 beta^(2/3) gamma^(1/3)/(2^(2/3) alpha),
    Assumptions->alpha>0 && beta>0 && gamma>0]];
  AssertZero["Hayward/extremal-horizon",FullSimplify[metric /. {r->radiusE,mu->massE},
    Assumptions->alpha>0 && beta>0 && gamma>0]];
  AssertZero["Hayward/extremal-stationarity",FullSimplify[D[metric,r] /. {r->radiusE,mu->massE},
    Assumptions->alpha>0 && beta>0 && gamma>0]];
  AssertZero["Hayward/extremal-curvature",FullSimplify[(r^2 D[metric,{r,2}]/2) /.
    {r->radiusE,mu->massE},Assumptions->alpha>0 && beta>0 && gamma>0]-1];
  AssertZero["Hayward/extremal-coefficient-speed",FullSimplify[(kt/fCoupling /. xx->profile) /.
    r->radiusE,Assumptions->alpha>0 && beta>0 && gamma>0]+1];
  AssertZero["Hayward/F2-zero",FullSimplify[(fCoupling /. xx->profile) /. r->radius0,
    Assumptions->alpha>0 && beta>0 && gamma>0]];
  threshold=FullSimplify[(massOfRadius /. r->radius0)/massE,
    Assumptions->alpha>0 && beta>0 && gamma>0];
  AssertZero["Hayward/mass-window-endpoint",threshold-2 (2/5)^(2/3)];
  normalizedMetric=FullSimplify[metric /. {r->radiusE radiusRatio,mu->massE massRatio},
    Assumptions->alpha>0 && beta>0 && gamma>0];
  horizonMassDimensionless=massRatio /. First[Solve[normalizedMetric==0,massRatio]];
  AssertZero["Hayward/normalized-horizon-curve",horizonMassDimensionless-(2 radiusRatio^3+1)/(3 radiusRatio^2)];
  AssertTrue["Hayward/outer-branch-monotone",
    FullSimplify[D[horizonMassDimensionless,radiusRatio]>0,Assumptions->radiusRatio>1]];
  AssertTrue["Hayward/exterior-sign-interval",
    FullSimplify[(6 y-5)/y^2<0 && 1/y^2>0,Assumptions->2/3<y<5/6]];
  AssertTrue["examples/trapped-real-scalar-exclusion",
    Reduce[xPositive>0 && fNegative<0 && scalarSlope^2==xPositive/fNegative,
      {xPositive,fNegative,scalarSlope},Reals]===False];

  (* z = Pi r^p/(2 sigma^(p-1)); total radial derivative keeps sigma,p fixed.
     This representation avoids branch-ambiguous inverse substitutions. *)
  xz=2 ArcTan[z]/Pi;
  dz[expr_] := D[expr,r]+p z/r D[expr,z];
  dx=dz[xz]; dxx=dz[dx];
  trigRules={Sin[2 ArcTan[z]]->2 z/(1+z^2),
    Sin[4 ArcTan[z]]->4 z (1-z^2)/(1+z^2)^2};
  dp=p Sin[Pi xx]-3 Pi xx;
  np=p^2 Sin[2 Pi xx]-8 p Sin[Pi xx]+6 Pi xx;
  fc=np/dp^2; kc=-2/dp;
  metricZ=1-2 mu xz/r;
  arctanCurrent=Together[r^4 ((kc /. xx->xz) /. trigRules) dz[metricZ/r^2]];
  arctanResidual=Together[dz[arctanCurrent]+2 ((fc /. xx->xz) /. trigRules)];
  AssertZero["arctangent/background-identity",arctanResidual];
  AssertZero["arctangent/logarithmic-profile-derivative",
    Together[r dx/xz-p ((Sin[Pi xx]/(Pi xx) /. xx->xz) /. trigRules)]];
  (* f=0 fixes mu=r/(2X). Its derivative is then (1-r X'/X)/r. *)
  AssertZero["arctangent/horizon-extremality-condition",
    Together[(r dz[metricZ] /. mu->r/(2 xz))-(1-r dx/xz)]];
  secondAtHorizon=Together[(r^2 dz[dz[metricZ]]/2) /. mu->r/(2 xz)];
  trigAssumptions=0<xx<1 && p>=3;
  logSlope=FullSimplify[(r dx/xz) /. z->Tan[Pi xx/2],Assumptions->trigAssumptions];
  AssertZero["arctangent/trigonometric-extremality",logSlope-p Sin[Pi xx]/(Pi xx)];
  rootRule=Sin[Pi xx]->Pi xx/p;
  de=dp /. rootRule;
  ne=Expand[TrigExpand[np]] /. rootRule;
  AssertZero["arctangent/extremal-D",de+2 Pi xx];
  AssertZero["arctangent/extremal-N",ne-2 Pi xx (p Cos[Pi xx]-1)];
  metricAtExtremality=FullSimplify[secondAtHorizon /. z->Tan[Pi xx/2],
    Assumptions->trigAssumptions];
  metricAtExtremality=FullSimplify[TrigExpand[metricAtExtremality] /. rootRule,
    Assumptions->trigAssumptions];
  (* A second form based on logarithmic differentiation has only sin(Pi X).
     This equality is off shell in X, before imposing the extremality root. *)
  AssertZero["arctangent/metric-second-derivative-profile",
    Together[secondAtHorizon-(-r^2 dxx/(2 xz)+r dx/xz-1)]];
  AssertZero["arctangent/extremal-curvature",FullSimplify[metricAtExtremality-(1-p Cos[Pi xx])/2,
    Assumptions->trigAssumptions]];
  speed=Cancel[(-2/de)/(ne/de^2)];
  AssertZero["arctangent/extremal-speed",speed-2/(p Cos[Pi xx]-1)];
  AssertZero["arctangent/geometric-coefficient-agreement",speed+1/((1-p Cos[Pi xx])/2)];
  (* Calculus bridge for every real p>=3: g(0)=0 and g'>0 imply
     sinc'(t)=-g(t)/t^2<0 on (0,Pi). Thus sinc(t)=1/p has one
     root. The endpoint signs below place it in (Pi/2,Pi), where
     2/(p Cos(t)-1)<0. The code verifies the ingredients; the stated
     monotonicity/continuity inference is ordinary real calculus. *)
  sincNumerator=Sin[t]-t Cos[t];
  rootFunction=p Sin[t]-t;
  AssertZero["arctangent/sign-proof/g-at-zero",sincNumerator /. t->0];
  AssertZero["arctangent/sign-proof/g-derivative",D[sincNumerator,t]-t Sin[t]];
  AssertTrue["arctangent/sign-proof/g-derivative-positive",
    FullSimplify[D[sincNumerator,t]>0,Assumptions->0<t<Pi]];
  AssertZero["arctangent/sign-proof/sinc-derivative",D[Sin[t]/t,t]+sincNumerator/t^2];
  AssertTrue["arctangent/sign-proof/root-endpoint-signs",
    FullSimplify[(rootFunction /. t->Pi/2)>0 && (rootFunction /. t->Pi)<0,
      Assumptions->p>=3]];
  AssertTrue["arctangent/sign-proof/negative-speed",
    FullSimplify[2/(p Cos[t]-1)<0,Assumptions->p>=3 && Pi/2<t<Pi]];

  <|"ExtremalLeadingIdentity"->leading,"ExtremalCoefficientRule"->extremalRule,
    "Hayward"-><|"F2"->fCoupling,"A1"->a1,"KT"->kt,"X"->profile,
      "Metric"->metric,"Current"->current,"ExtremalRadius"->radiusE,
      "ExtremalMass"->massE,"TensorZeroRadius"->radius0,
      "MassRatioThreshold"->threshold,"NormalizedHorizonCurve"->horizonMassDimensionless|>,
    "Arctangent"-><|"ProfileInZ"->xz,"F2"->fc,"KT"->kc,
      "ExtremalityEquation"->HoldForm[p Sin[Pi xx]==Pi xx],
      "ExtremalRadiusSquaredF2"->(1-p Cos[Pi xx])/2,"ExtremalSpeedSquared"->speed,
      "SignArgument"->"On (0,Pi), g=sin(t)-t cos(t) has g(0)=0 and g'=t sin(t)>0, so sinc strictly decreases. The unique root p sin(t)=t lies in (Pi/2,Pi) for p>=3, giving negative extremal speed squared."|>,
    "Domain"->"Real symmetry-inheriting scalar on the open static exterior. Published tensor coefficients are inputs; the code checks their identities and consequences."|>
]
