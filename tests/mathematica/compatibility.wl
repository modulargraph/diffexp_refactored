Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "DiffExp.m"}]];
failures={}; timings=<||>;
check[condition_,name_] := If[!TrueQ[condition],AppendTo[failures,name]];
close[a_,b_,tolerance_:10^-20] := NumberQ[a] && TrueQ[Abs[a-b]<tolerance];
measured[name_,expression_] := Null;
SetAttributes[measured,HoldRest];
measured[name_,expression_] := Module[{seconds,value},{seconds,value}=AbsoluteTiming[TimeConstrained[expression,60,$Aborted]];timings[name]=seconds;value];
folder=CreateDirectory[];
Put[{{1/(1+t)}},FileNameJoin[{folder,"dt_0.m"}]];
Put[{{1/(1+t)}},FileNameJoin[{folder,"dt_1.m"}]];
config=LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->2,WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
check[AssociationQ[config] && config[Variables]==={t},"configuration matrix variables"];
UpdateConfiguration["ExpansionOrder"->44];check[CurrentConfiguration[][ExpansionOrder]===44,"string option update"];
b=PrepareBoundaryConditions[{{1,0,0}},{t->0}];
r=measured["regular",TransportTo[b,{t->1}]];
check[ListQ[r] && Length[r]===3,"ordinary return shape"];
If[ListQ[r],check[And@@MapThread[close,{r[[2,1]],{2,2 Log[2],Log[2]^2}}],"ordinary epsilon coefficients"]];
r2=measured["chained",TransportTo[r,{t->3}]];
If[ListQ[r2],check[close[r2[[2,1,1]],4],"chained boundary"],check[False,"chained result"]];
saved=measured["saved",TransportTo[b,{t->x},1,True]];
functions=ToPiecewise[saved,True];
check[MatrixQ[functions] && close[functions[[1,1]][1/2],3/2],"saved physical-point replay"];
check[close[functions[[1,1]][x]/.x->1/3,4/3],"plot substitution preserves saved metadata"];
check[AssociationQ[DiffExpLastTimings[]],"timing accessor"];
check[close[functions[[1,1]][0.25],5/4],"machine-precision plotting coordinate"];
check[DiffExp`Private`$lastRequest["boundary"][[1,1]]=!="1","plot reuses a nearby physical sample"];
check[close[functions[[1,1]][N[1/7,10000]],8/7],"plot precision bounded by integration precision"];
plot=measured["plot",Plot[Evaluate[functions[[1,1]][x]],{x,0.1,0.9},PlotPoints->3,MaxRecursion->0]];
check[Head[plot]===Graphics && !FreeQ[plot,_Line],"original Plot function workflow"];

check[FailureQ[PrepareBoundaryConditions[{1/eps},{t->0}]],"negative epsilon window rejected"];
check[FailureQ[IntegrateSystem[{t->x}]],"unsupported symbolic general solution explicit"];
check[FailureQ[Catch[DiffExp`Private`decimalValue["1;Print[42]"],"DiffExpCompatibility"]],"decimal response cannot execute code"];
check[close[DiffExp`Private`decimalValue["-1.234e-17"],-1234/10^20,10^-35],"decimal precision parser"];
DeleteFile[FileNames["*.m",folder]];
Put[{{s}},FileNameJoin[{folder,"dt_0.m"}]];
LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->0,WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
r=measured["nonlinear_path",TransportTo[PrepareBoundaryConditions[{1},{t->0,s->0}],{t->x,s->x^2}]];
If[ListQ[r],check[close[r[[2,1,1]],Exp[1/3],10^-16],"explicit nonlinear path preserved"],check[False,"nonlinear path result"]];
DeleteFile[FileNames["*.m",folder]];
Put[{{Log[1+t]}},FileNameJoin[{folder,"d_1.m"}]];
LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->2,WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
r=measured["canonical",TransportTo[PrepareBoundaryConditions[{1},{t->0}],{t->1}]];
If[ListQ[r],check[close[r[[2,1,3]],Log[2]^2/2],"canonical dlog weights"],check[False,"canonical result"]];
(* Explicit algebraic basis data survives process dispatch, chaining and saves. *)
DeleteFile[FileNames["*.m",folder]];
Put[{{1/(2 Sqrt[t](2+Sqrt[t]))}},FileNameJoin[{folder,"dt_0.m"}]];
LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->0,WorkingPrecision->90,ExpansionOrder->40,
 AccuracyGoal->30,BasisPrefactors->{2+Sqrt[t]},Verbosity->0}];
loop=1-8 x(1-x)+8 I x(1-x)(1-2 x);
b=PrepareBoundaryConditions[{1},{t->1}];
saved=measured["algebraic_basis_loop",TransportTo[b,{t->loop},1,True]];
check[ListQ[saved] && close[saved[[1,2,1,1]],1],"principal endpoint basis after root winding"];
If[ListQ[saved],r=TransportTo[saved[[1]],{t->4}];
 check[ListQ[r] && close[r[[2,1,1]],4/3],"principal basis chained boundary"]];
functions=ToPiecewise[saved];UpdateConfiguration[BasisPrefactors->{}];
check[MatrixQ[functions] && close[functions[[1,1]][1/4],(2+Sqrt[loop/.x->1/4])/3],"saved basis configuration"];
DeleteFile[FileNames["*.m",folder]];
Put[{{0,1/t},{1/t,0}},FileNameJoin[{folder,"dt_0.m"}]];
LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->0,WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
b=PrepareBoundaryConditions[{x+O[x]^2,"?"},{t->x}];
r=measured["partial_asymptotic",TransportTo[b,{t->1/2}]];
If[ListQ[r],check[close[r[[2,1,1]],1/2]&&close[r[[2,2,1]],1/2],"partial asymptotic coupled solve"],check[False,"partial asymptotic result"]];

(* Full shuffle-based multiple-polylogarithm workflow, with a bounded example. *)
Get[FileNameJoin[{DirectoryName[$InputFileName],"mpl_workflow.wl"}]];
mpl=measured["mpl_shuffle",GEvaluate[1,0,1/2]];
check[close[mpl,Log[1/2]Log[1-1/2]+PolyLog[2,1/2],10^-16],"full MPL shuffle workflow"];

(* Prepare the original partial gamma/log boundary without integrating a large system. *)
DeleteFile[FileNames["*.m",folder]];Put[ConstantArray[0,{4,4}],FileNameJoin[{folder,"dt_0.m"}]];
LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->4,WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
partial=measured["gamma_partial_preparation",PrepareBoundaryConditions[{"?","?",
 eps(1+3 eps)(1+4 eps)(-4 Exp[3 EulerGamma eps] Gamma[eps]^3/t+
 6 Exp[3 EulerGamma eps](-1/t)^(1+eps) eps Gamma[-eps]^2 Gamma[eps]^3/Gamma[-2 eps]+
 8 Exp[3 EulerGamma eps](-1/t)^(1+2 eps) eps Gamma[-eps]^3 Gamma[eps] Gamma[2 eps]/Gamma[-3 eps]+
 3 Exp[3 EulerGamma eps](-1/t)^(1+3 eps) eps Gamma[-eps]^4 Gamma[3 eps]/Gamma[-4 eps]),
 Exp[3 EulerGamma eps]eps^3 Gamma[eps]^3},{t->-1/x}]];
check[ListQ[partial] && Dimensions[partial[[2]]]==={4,5},"gamma partial boundary shape"];
If[ListQ[partial],asym=Catch[DiffExp`Private`asymptoticData[partial[[2]],x],"DiffExpCompatibility"];check[AssociationQ[asym]&&Length[asym["constraints"]]>0,"gamma power-log normalization"]];

(* Kinematic x must not collide with the native path parameter. *)
DeleteFile[FileNames["*.m",folder]];Put[{{1/(1-x)}},FileNameJoin[{folder,"dx_0.m"}]];
LoadConfiguration[{Variables->{x},LineParameter->t,MatrixDirectory->folder,EpsilonOrder->0,AccuracyGoal->25,
 WorkingPrecision->90,ChopPrecision->60,ExpansionOrder->40,Verbosity->0}];
b=PrepareBoundaryConditions[{1},{x->0}];r=measured["kinematic_x",TransportTo[b,{x->1/2}]];
check[ListQ[r]&&close[r[[2,1,1]],2],"kinematic x with separate line parameter"];

family=DiffExpFamilyTemplate["bubble"];
check[AssociationQ[family] && family["schema"]==="DiffExp.FeynmanFamily/v1","editable FT template"];
If[AssociationQ[family],family["name"]="custom-wrapper-test";
 ft=measured["generic_ft",DiffExpFeynmanTrick[family,{"--epsilon-order","0","--leaf-digits","12"}]];
 check[AssociationQ[ft],"generic FT association dispatch"];
 check[NumberQ[DiffExpLastTimings[]["total_seconds"]],"generic FT timing accessor"];
 familyFile=FileNameJoin[{folder,"family configuration.json"}];Export[familyFile,family,"RawJSON"];
 ftFile=measured["generic_ft_file",DiffExpFeynmanTrick[familyFile,{"--epsilon-order","0","--leaf-digits","12"}]];
 check[AssociationQ[ftFile],"generic FT JSON file dispatch"]];
DeleteDirectory[folder,DeleteContents->True];
Print["Compatibility timings: ",InputForm[timings]];
If[failures==={},Print["Native Mathematica compatibility: all checks passed"];Exit[0],Print["Failures: ",failures];Exit[1]];
