(* Numerical checks for SM S7 and the two panels of Letter Figure 2.
   FindRoot solves the original metric horizon conditions at two working
   precisions. Coefficient ratios and radial metric derivatives are evaluated
   independently. No field evolution or real trapped scalar is assumed.
   Outputs are confined to $GeometricNoGoOutputDirectory (runner setting),
   or to the adjacent outputs/ directory when evaluated directly. *)
Module[{root=DirectoryName[$InputFileName],out,wp,tol,r,mu,xx,p,sigma,
  parameters,alpha,beta,gamma,ell,profile,metric,derivative,second,solution,
  rn,mn,xn,coeffSpeed,geometrySpeed,originalResidual,radiusExact,massExact,
  profileRoot,radiusFromProfile,denominator,numerator,fCoupling,kCoupling,
  precisionRuns=<||>,arctanRecords={},haywardRecords={},scanRecords={},
  key,record,id,fmt,near,profiles,values,low,high,m,rho,roots,rplus,
  boundary=(5/2)^(1/3),threshold=2 (2/5)^(2/3),sample,coef,band,
  horizonRows,coefficientRadii,coefficientRows,bandRows,negativeRows,positiveRows,horizonPlot,coefficientPlot,
  radiusBoundaryLabel,massBoundaryLabel,
  combinedFigure,red,ink,gray,bandFill,insideFill,paperStyle,tableNames,exports={},exportChecked,
  numericReport,plotMass,plotRadius,plotCoefficient,files,mm,rr,hh},
  out=If[StringQ[$GeometricNoGoOutputDirectory],$GeometricNoGoOutputDirectory,
    FileNameJoin[{root,"outputs"}]];
  If[!DirectoryQ[out],CreateDirectory[out,CreateIntermediateDirectories->True]];
  fmt[value_] := ToString[NumberForm[N[value,30],30,
    NumberFormat->(If[#3=="",#1,#1<>"e"<>#3]&)],OutputForm];
  near[name_,residual_,tolerance_] := AssertNear[name,Max[Abs[Flatten[{residual}]]],tolerance];
  exportChecked[name_,expr_,format_] := Module[{path=FileNameJoin[{out,name}],result},
    result=Export[path,expr,format];
    AssertTrue["exports/"<>name,StringQ[result] && FileExistsQ[path] && FileByteCount[path]>0];
    AppendTo[exports,name];result];

  Do[
    tol=10^(-Floor[wp/2]);
    (* Couplings remain exact until passed to the high-precision solver. *)
    Do[
      {alpha,beta,gamma}=parameters;
      ell=(gamma/beta)^(1/3);
      profile=alpha r^3/(beta r^3+gamma);
      metric=1-mu alpha r^2/(beta r^3+gamma);
      derivative=D[metric,r]; second=D[metric,{r,2}];
      solution=FindRoot[{metric==0,derivative==0},
        {{r,N[6 ell/5,wp]},{mu,N[2 beta ell/alpha,wp]}},
        WorkingPrecision->wp,AccuracyGoal->wp-10,PrecisionGoal->wp-10,MaxIterations->100];
      rn=r /. solution; mn=mu /. solution; xn=profile /. r->rn;
      (* h0=-1 is a harmless nonzero common coefficient normalization. *)
      fCoupling=-alpha (5 alpha-6 beta xx)/xx^2;
      kCoupling=Cancel[fCoupling-xx (-6 alpha (alpha-beta xx)/xx^3)];
      coeffSpeed=(kCoupling/fCoupling) /. xx->xn;
      geometrySpeed=-2/(rn^2 (second /. solution));
      originalResidual=Max[Abs[N[{metric,ell derivative} /. solution,wp]]];
      radiusExact=(2 gamma/beta)^(1/3);
      massExact=3 beta^(2/3) gamma^(1/3)/(2^(2/3) alpha);
      id="numerics/Hayward/"<>ToString[parameters,InputForm]<>"/wp="<>ToString[wp];
      near[id<>"/original-horizon-equations",originalResidual,tol];
      near[id<>"/radius",rn/radiusExact-1,tol];
      near[id<>"/mass",mn/massExact-1,tol];
      near[id<>"/coefficient-versus-metric",coeffSpeed-geometrySpeed,tol];
      near[id<>"/extremal-minus-one",coeffSpeed+1,tol];
      record=<|"WorkingPrecision"->wp,"alpha"->fmt[alpha],"beta"->fmt[beta],"gamma"->fmt[gamma],
        "ExtremalRadius"->fmt[rn],"ExtremalMass"->fmt[mn],
        "CoefficientSpeedSquared"->fmt[coeffSpeed],"MetricSpeedSquared"->fmt[geometrySpeed],
        "OriginalEquationResidual"->fmt[originalResidual],"AgreementResidual"->fmt[Abs[coeffSpeed-geometrySpeed]],
        "Tolerance"->fmt[tol]|>;
      AppendTo[haywardRecords,record];
      key="Hayward-"<>ToString[parameters,InputForm]<>"-"<>ToString[wp];
      AssociateTo[precisionRuns,key->{rn,mn,coeffSpeed}],
      {parameters,{{1,1,1},{2,3,5},{7/10,4,1/5}}}];

    Do[
      sigma=1;
      profile=2/Pi ArcTan[Pi r^p/(2 sigma^(p-1))];
      metric=1-2 mu profile/r;
      derivative=D[metric,r]; second=D[metric,{r,2}];
      solution=FindRoot[{metric==0,derivative==0},
        {{r,N[6/5,wp]},{mu,N[4/5,wp]}},
        WorkingPrecision->wp,AccuracyGoal->wp-10,PrecisionGoal->wp-10,MaxIterations->100];
      rn=r /. solution; mn=mu /. solution; xn=profile /. r->rn;
      (* Second root solve uses the analytically reduced scalar equation,
         and never provides initial values to the original metric solver. *)
      profileRoot=xx /. FindRoot[p Sin[Pi xx]==Pi xx,{xx,N[3/4,wp]},
        WorkingPrecision->wp,AccuracyGoal->wp-10,PrecisionGoal->wp-10];
      radiusFromProfile=(2 sigma^(p-1) Tan[Pi profileRoot/2]/Pi)^(1/p);
      denominator=p Sin[Pi xn]-3 Pi xn;
      numerator=p^2 Sin[2 Pi xn]-8 p Sin[Pi xn]+6 Pi xn;
      coeffSpeed=(-2/denominator)/(numerator/denominator^2);
      geometrySpeed=-2/(rn^2 (second /. solution));
      originalResidual=Max[Abs[N[{metric,derivative} /. solution,wp]]];
      id="numerics/arctangent/p="<>ToString[p]<>"/wp="<>ToString[wp];
      near[id<>"/original-horizon-equations",originalResidual,tol];
      near[id<>"/independent-profile-root",xn-profileRoot,tol];
      near[id<>"/independent-radius",rn-radiusFromProfile,tol];
      near[id<>"/coefficient-versus-metric",coeffSpeed-geometrySpeed,tol];
      AssertTrue[id<>"/negative-speed-and-ordinary-extremum",coeffSpeed<0 && (second /. solution)>0 && 0<xn<1];
      record=<|"p"->p,"sigma"->sigma,"WorkingPrecision"->wp,
        "ExtremalX"->fmt[xn],"ExtremalRadius"->fmt[rn],"ExtremalMass"->fmt[mn],
        "CoefficientSpeedSquared"->fmt[coeffSpeed],"MetricSpeedSquared"->fmt[geometrySpeed],
        "OriginalEquationResidual"->fmt[originalResidual],"AgreementResidual"->fmt[Abs[coeffSpeed-geometrySpeed]],
        "Tolerance"->fmt[tol]|>;
      AppendTo[arctanRecords,record];
      key="arctangent-"<>ToString[p]<>"-"<>ToString[wp];
      AssociateTo[precisionRuns,key->{rn,mn,coeffSpeed}],
      {p,{3,4}}],
    {wp,{40,70}}];
  (* Convergence is measured on independently executed 40/70 digit runs. *)
  Do[
    low=precisionRuns["Hayward-"<>ToString[parameters,InputForm]<>"-40"];
    high=precisionRuns["Hayward-"<>ToString[parameters,InputForm]<>"-70"];
    near["numerics/precision-convergence/Hayward/"<>ToString[parameters,InputForm],low-high,10^-20],
    {parameters,{{1,1,1},{2,3,5},{7/10,4,1/5}}}];
  Do[
    near["numerics/precision-convergence/arctangent/p="<>ToString[p],
      precisionRuns["arctangent-"<>ToString[p]<>"-40"]-precisionRuns["arctangent-"<>ToString[p]<>"-70"],10^-20],
    {p,{3,4}}];

  (* Direct cubic roots supply an independent check of the parametric plot. *)
  Do[
    roots=Sort[rho /. NSolve[2 rho^3-3 m rho^2+1==0 && rho>0,rho,Reals,WorkingPrecision->70]];
    AssertTrue["numerics/Hayward/scan/two-horizons/m="<>ToString[m,InputForm],Length[roots]==2];
    rplus=Last[roots];
    near["numerics/Hayward/scan/metric-residual/m="<>ToString[m,InputForm],
      1-3 m rplus^2/(2 rplus^3+1),10^-35];
    near["numerics/Hayward/scan/parametric-agreement/m="<>ToString[m,InputForm],
      (2 rplus^3+1)/(3 rplus^2)-m,10^-35];
    band=TrueQ[rplus<boundary];
    AssertTrue["numerics/Hayward/scan/exclusion-equivalence/m="<>ToString[m,InputForm],
      band===TrueQ[m<threshold]];
    sample=If[band,(rplus+boundary)/2,101 rplus/100];
    coef=(2 sample^3-5)/(2 sample^3+1);
    AssertTrue["numerics/Hayward/scan/exterior-sign/m="<>ToString[m,InputForm],
      sample>rplus && If[band,coef<0,coef>0]];
    AppendTo[scanRecords,<|"MassRatio"->fmt[m],"InnerRadiusRatio"->fmt[First[roots]],
      "OuterRadiusRatio"->fmt[rplus],"TensorZeroRadiusRatio"->fmt[boundary],
      "ExcludedExteriorExists"->band,"ExteriorSampleRadiusRatio"->fmt[sample],
      "F2OverKTAtSample"->fmt[coef]|>],
    {m,{1001/1000,21/20,27/25,109/100,23/20}}];

  (* Sample the same dimensionless analytic functions as the Letter Fig. 2. *)
  horizonRows=Table[N[{(2 rr^3+1)/(3 rr^2),rr},24],{rr,1,9/5,1/1500}];
  coefficientRadii=Append[Range[1,7/4,7/4000],7/4];
  coefficientRows=Table[N[{rr,(2 rr^3-5)/(2 rr^3+1)},24],{rr,coefficientRadii}];
  bandRows=Table[N[{(2 rr^3+1)/(3 rr^2),rr,boundary},24],{rr,Subdivide[1,boundary,400]}];
  exportChecked["fig2_horizon.csv",Prepend[horizonRows,{"mu_over_mu_e","r_plus_over_r_e"}],"CSV"];
  exportChecked["fig2_coefficient.csv",Prepend[coefficientRows,{"r_over_r_e","F2_over_KT"}],"CSV"];
  exportChecked["fig2_exclusion_band.csv",Prepend[bandRows,{"mu_over_mu_e","r_plus_over_r_e","r0_over_re"}],"CSV"];
  exportChecked["extremal_hayward.csv",Prepend[Values /@ haywardRecords,Keys[First[haywardRecords]]],"CSV"];
  exportChecked["extremal_arctangent.csv",Prepend[Values /@ arctanRecords,Keys[First[arctanRecords]]],"CSV"];
  exportChecked["hayward_mass_scan.csv",Prepend[Values /@ scanRecords,Keys[First[scanRecords]]],"CSV"];

  red=RGBColor[.702,.267,.212]; ink=RGBColor[.141,.216,.275]; gray=RGBColor[.486,.529,.565];
  bandFill=RGBColor[.961,.867,.843]; insideFill=RGBColor[.906,.929,.941];
  paperStyle={Frame->{{True,False},{True,False}}, Axes->False, Background->White,
    FrameStyle->gray, FrameTicksStyle->Directive[ink,12], LabelStyle->Directive[ink,12],
    BaseStyle->{FontFamily->"Arial",FontSize->12,FontColor->ink}, PlotRangePadding->None,
    ImagePadding->{{72,22},{43,28}}, ImageSize->520};
  (* Keep the analytic boundary labels exact and unevaluated for display. *)
  radiusBoundaryLabel=Superscript[Row[{"(",5,"/",2,")"}],"1/3"];
  massBoundaryLabel=Row[{2,Superscript[Row[{"(",2,"/",5,")"}],"2/3"]}];
  horizonPlot=ListLinePlot[{horizonRows,{{1,boundary},{9/8,boundary}}},
    PlotStyle->{Directive[ink,AbsoluteThickness[2]],Directive[red,Dashed,AbsoluteThickness[1]]},
    FrameLabel->{Row[{"\[Mu] / ",Subscript["\[Mu]","e"],"  (fixed couplings)"}],Row[{"r / ",Subscript["r","e"]}]},
    FrameTicks->{{{1,1.2,{boundary,radiusBoundaryLabel},1.6},None},
      {{1,1.04,{threshold,massBoundaryLabel},1.12},None}},
    PlotRange->{{1,9/8},{.985,1.61}},AspectRatio->.56,
    PlotLabel->Style["(b)  Exterior exclusion band",15,Bold,Black],
    Prolog->{EdgeForm[None],insideFill,
      Polygon[Join[{{1,.985}},horizonRows,{{horizonRows[[-1,1]],.985}}]],
      bandFill,Polygon[Join[bandRows[[All,{1,2}]],Reverse[bandRows[[All,{1,3}]]]]]},
    Epilog->{{gray,Dotted,Line[N[{{threshold,.985},{threshold,boundary}}]]},
      red,PointSize[.014],Point[N[{threshold,boundary}]],
      Line[N[{{1.066,1.51},{threshold,boundary}}]],
      Text[Style[Row[{Subscript["r","+"]," = ",Subscript["r","0"]}],12,red],{1.062,1.55}],
      Text[Style["Excluded\nexterior",12,red,TextAlignment->Center],{1.022,1.29}],
      Text[Style["Inside outer horizon",12,ink],{1.067,1.095}],
      Text[Rotate[Style[Row[{Subscript["r","+"],"(\[Mu])"}],12,ink],17 Degree],{1.104,1.445}]},
    Evaluate[Sequence@@paperStyle]];
  negativeRows=Append[Select[coefficientRows,#[[1]]<boundary&],N[{boundary,0},24]];
  positiveRows=Prepend[Select[coefficientRows,#[[1]]>boundary&],N[{boundary,0},24]];
  coefficientPlot=ListLinePlot[{negativeRows,positiveRows},
    PlotStyle->{Directive[red,AbsoluteThickness[2]],Directive[ink,AbsoluteThickness[2]]},
    Filling->{1->{0,bandFill}},
    FrameLabel->{Row[{"r / ",Subscript["r","e"]}],Row[{Subscript["F",2]," / ",Subscript["K","T"]}]},
    FrameTicks->{{{-1,0,0.5},None},{{1,{boundary,radiusBoundaryLabel},1.7},None}},
    PlotRange->{{1,7/4},{-1.12,.57}},AspectRatio->.42,
    PlotLabel->Style["(a)  Finite tensor ratio",15,Bold,Black],
    Epilog->{gray,Line[{{1,0},{7/4,0}}],
      {red,Dotted,Line[N[{{boundary,-1.12},{boundary,.57}}]]},
      red,PointSize[.014],Point[N[{boundary,0}]],
      Text[Style[Subscript["r","0"],12,red],{boundary-.035,.40}],
      Text[Style[Row[{Subscript["F",2]," < 0 < ",Subscript["K","T"]}],12,red],{1.135,-.17}]},
    Evaluate[Sequence@@paperStyle]];
  combinedFigure=GraphicsColumn[{coefficientPlot,horizonPlot},Background->White,ImageSize->600,Spacings->6];
  exportChecked["hayward_exterior.pdf",combinedFigure,"PDF"];
  exportChecked["hayward_exterior.png",combinedFigure,"PNG"];
  numericReport=<|"WorkingPrecisions"->{40,70},"NumericalToleranceRule"->"10^(-Floor[WorkingPrecision/2])",
    "HighPrecisionValues"->"Portable decimal strings retain up to 30 digits; working precision and test tolerance are recorded separately.",
    "Hayward"->haywardRecords,"Arctangent"->arctanRecords,"MassScan"->scanRecords,
    "MassRatioThreshold"->fmt[threshold],"TensorZeroRadiusRatio"->fmt[boundary],
    "Domain"->"Real radial-scalar exterior. Curves sample analytic functions; no numerical spacetime evolution.",
    "Files"->exports|>;
  exportChecked["numerical_results.json",numericReport,"RawJSON"];
  <|"HaywardExtremalChecks"->haywardRecords,"ArctangentExtremalChecks"->arctanRecords,
    "MassScan"->scanRecords,"Exports"->exports,"OutputDirectory"->"outputs/",
    "Figure"->combinedFigure,"MassRatioThreshold"->threshold|>
]
