root=ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName],"..",".."}]];
Get[FileNameJoin[{root,"DiffExp.m"}]];
assert[c_,m_]:=If[!TrueQ[c],If[AssociationQ[DiffExp`Private`$lastRequest],Export[FileNameJoin[{root,"build-reference","banana-wrapper-request.json"}],DiffExp`Private`$lastRequest,"RawJSON"]];Print["FAIL: ",m];Quit[1]];
started=AbsoluteTime[];
base=FileNameJoin[{root,"examples","original"}];
configuration={MatrixDirectory->FileNameJoin[{base,"Data","Banana","EqualMass"}],EpsilonOrder->4,WorkingPrecision->150,ExpansionOrder->60,DivisionOrder->4,DeltaPrescriptions->{t-16+I δ},UseMobius->True,UsePade->True,Verbosity->1};
loaded=LoadConfiguration[configuration];assert[!FailureQ[loaded],loaded];
{preparation,bc}=AbsoluteTiming[PrepareBoundaryConditions[{"?","?",
 eps(1+3 eps)(1+4 eps)(-4 Exp[3 EulerGamma eps] Gamma[eps]^3/t+
 6 Exp[3 EulerGamma eps](-1/t)^(1+eps) eps Gamma[-eps]^2 Gamma[eps]^3/Gamma[-2 eps]+
 8 Exp[3 EulerGamma eps](-1/t)^(1+2 eps) eps Gamma[-eps]^3 Gamma[eps] Gamma[2 eps]/Gamma[-3 eps]+
 3 Exp[3 EulerGamma eps](-1/t)^(1+3 eps) eps Gamma[-eps]^4 Gamma[3 eps]/Gamma[-4 eps]),
 Exp[3 EulerGamma eps]eps^3 Gamma[eps]^3},{t->-1/x}]];
assert[!FailureQ[bc],bc];
Print["Original asymptotic input prepared in ",preparation," s"];
minus=TransportTo[bc,{t->-1}];assert[!FailureQ[minus],minus];
reference=Take[#,5]&/@Get[FileNameJoin[{base,"Reference","BananaBoundaryAtMinusOneEps7.m"}]];
error=Max[Abs[minus[[2]]-reference]];Print["Equal banana at -1 discrepancy: ",error];assert[error<10^-20,error];
Put[minus,FileNameJoin[{root,"build-reference","banana-wrapper-minus-one.m"}]];
report=<|"schema"->"DiffExp.OriginalWrapperValidation/v1","family"->"banana-equal-asymptotic","status"->"passed","max_discrepancy"->ToString[error,InputForm],"boundary_preparation_seconds"->preparation,"timings"->Join[DiffExpLastTimings[],<|"test_total_seconds"->AbsoluteTime[]-started|>]|>;
Export[FileNameJoin[{root,"build-reference","banana-wrapper-asymptotic.json"}],report,"RawJSON"];
Print[ExportString[report,"RawJSON"]];Quit[0];
