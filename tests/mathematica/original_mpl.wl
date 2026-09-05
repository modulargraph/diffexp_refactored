root=ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName],"..",".."}]];
Get[FileNameJoin[{root,"DiffExp.m"}]];
$MPLWorkingPrecision=250;$MPLExpansionOrder=75;
Get[FileNameJoin[{DirectoryName[$InputFileName],"mpl_workflow.wl"}]];
folder=CreateDirectory[];report={};
assert[c_,m_]:=If[!TrueQ[c],Print["FAIL: ",m];If[AssociationQ[DiffExp`Private`$lastRequest],Export[FileNameJoin[{root,"build-reference","mpl-wrapper-request.json"}],DiffExp`Private`$lastRequest,"RawJSON"]];Quit[1]];
cases={{{1,0,1,4},-6.7782180257804207212554826775005988168291802221955692129682`55+I 0.9250147943833369547396749852220309435917997631163983727603`55},
 {{1,-10,0,4},-0.0191508840720296721365611597236750922866172732200324064383`55-I 0.3066358899483403657463434439014286049874538907865239005438`55},
 {{10,-10+I,-1/2,-50,1},-0.0000098802442781507281548895360764863423574760704710022738`50-I 0.0000009352314628872620198852585457725779647560856726857223`50},
 {Append[Range[20],21],0.00000000000513066731719907533179589918813462949766948546641803107466616580097`60}};
Do[If[Length[test[[1]]]>10,$MPLExpansionOrder=100];{seconds,value}=AbsoluteTiming[GEvaluate@@test[[1]]];assert[NumberQ[value],value];error=Abs[value-test[[2]]];assert[error<10^-30,error];row=<|"word_and_endpoint"->ToString[test[[1]],InputForm],"seconds"->seconds,"discrepancy"->ToString[error,InputForm],"status"->"passed"|>;AppendTo[report,row];Print[row],{test,cases}];
Export[FileNameJoin[{root,"build-reference","mpl-wrapper-results.json"}],report,"RawJSON"];
DeleteDirectory[folder,DeleteContents->True];Quit[0];
