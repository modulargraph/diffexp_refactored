(* Original notebook computations, with local ancillary data and the process wrapper. *)
root=ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName],"..",".."}]];
Get[FileNameJoin[{root,"DiffExp.m"}]];
family=Replace[Environment["DIFFEXP_ORIGINAL_FAMILY"],Except[_String]->"1loop"];
highPrecision=Environment["DIFFEXP_HIGH_PRECISION"]==="1";
work=FileNameJoin[{root,"build-reference","mathematica-original",family<>If[highPrecision,"-high",""]}];
If[!DirectoryQ[work],CreateDirectory[work,CreateIntermediateDirectories->True]];
assert[condition_,message_]:=If[!TrueQ[condition],If[AssociationQ[DiffExp`Private`$lastRequest],Export[FileNameJoin[{work,"request.json"}],DiffExp`Private`$lastRequest,"RawJSON"]];Print["FAIL: ",message];Quit[1]];
(* Optional split execution lets a long native run release the Wolfram license.
   Replay validates the exact request before accepting the external response. *)
requestOnly=Environment["DIFFEXP_REQUEST_ONLY"]==="1";
responseFile=Environment["DIFFEXP_RESPONSE_FILE"];
If[requestOnly,
 Clear[DiffExp`Private`nativeTransport];
 DiffExp`Private`nativeTransport[request_Association]:=(Export[FileNameJoin[{work,"request.json"}],request,"RawJSON"];Print["Native request prepared."];Quit[0])];
If[StringQ[responseFile] && responseFile=!="",
 Clear[DiffExp`Private`nativeTransport];
 DiffExp`Private`nativeTransport[request_Association]:=Module[{expected,data},
  expected=Import[FileNameJoin[{work,"request.json"}],"RawJSON"];assert[request===expected,"External response request does not match the wrapper input."];
  data=Import[responseFile,"RawJSON"];assert[AssociationQ[data]&&Lookup[data,"schema",""]==="DiffExp.TransportResult/v1","Invalid external transport response."];
  DiffExp`Private`$lastRequest=request;DiffExp`Private`$lastTimings=Join[Lookup[data,"timings",<||>],<|"external_native_process"->True|>];data]];
t0=AbsoluteTime[];
If[MemberQ[{"1loop","zmz","mzz","zzz"},family],
 data=FileNameJoin[{root,"examples","original","Data","PlanarOneMass"}];
 roots={tr5->Sqrt[(-s12 s15+s12 s23+p1s s34+s15 s45-s34 s45-s23 s34)^2-4 s23 s34 s45(p1s-s12-s15+s34)],
 sqrtG3->Sqrt[p1s^2+(s23-s45)^2-2 p1s(s23+s45)],sqrtG3nc->Sqrt[(s12+s15)^2-4 p1s s34]};
 alphabet=Get[FileNameJoin[{data,"alphabet.m"}]]/.roots;
 tensor=Get[FileNameJoin[{data,family,"diffEq-"<>family<>".m"}]];
 Put[Table[Log[alphabet[[i]]],{i,58}].Normal[tensor],FileNameJoin[{work,"d_1.m"}]];
 samples=Get[FileNameJoin[{data,family,"numIntegrals-"<>family<>".m"}]];
 ph1=MapAt[Transpose,samples[[1]],2];ph6=MapAt[Transpose,samples[[If[highPrecision,2,6]]],2];
 config={AccuracyGoal->15,ExpansionOrder->40,DeltaPrescriptions->{p1s+I δ,s12+I δ,s15+I δ,s23+I δ,s34+I δ,s45+I δ},MatrixDirectory->work,WorkingPrecision->150,ChopPrecision->100};
 If[highPrecision,config=Join[config,{AccuracyGoal->128,DivisionOrder->4,ExpansionOrder->230,WorkingPrecision->280,ChopPrecision->230}]];
 loaded=LoadConfiguration[config];assert[!FailureQ[loaded],loaded];
 boundary=PrepareBoundaryConditions[ph1[[2]],ph1[[1]]];assert[!FailureQ[boundary],boundary];
 result=TransportTo[boundary,ph6[[1]]];assert[!FailureQ[result],result];
 discrepancy=Max[Abs[result[[2]]-ph6[[2]]]];assert[discrepancy<If[highPrecision,10^-126,10^-14],discrepancy];,
 If[family==="henn",
 data=FileNameJoin[{root,"examples","original","Data","HennNonplanar"}];
 deltaGram=v1^2(v2-v5)^2+(v2 v3+v4(-v3+v5))^2+2 v1(-v2^2 v3+v4(v3-v5)v5+v2(v3 v4+(v3+v4)v5));
 sr=Sqrt[deltaGram];
 alphabet={v1,v2,v3,v4,v5,v3+v4,v4+v5,v1+v5,v1+v2,v2+v3,v1-v4,v2-v5,-v1+v3,-v2+v4,-v3+v5,v1+v2-v4,v2+v3-v5,-v1+v3+v4,-v2+v4+v5,v1-v3+v5,-v1-v2+v3+v4,-v2-v3+v4+v5,v1-v3-v4+v5,v1+v2-v4-v5,-v1+v2+v3-v5};
 bases={v1 v2-v2 v3+v3 v4-v1 v5-v4 v5,-v1 v2+v2 v3-v3 v4-v1 v5+v4 v5,-v1 v2-v2 v3+v3 v4+v1 v5-v4 v5,v1 v2-v2 v3-v3 v4-v1 v5+v4 v5,-v1 v2+v2 v3-v3 v4+v1 v5-v4 v5};
 alphabet=Join[alphabet,(#-sr)/(#+sr)&/@bases,{sr}];
 matrix=Get[FileNameJoin[{data,"XB_Atilde.txt"}]]/.Table[W[i]->alphabet[[i]],{i,31}];Put[matrix,FileNameJoin[{work,"d_1.m"}]];
 config={AccuracyGoal->30,ExpansionOrder->80,DeltaPrescriptions->{v1+I δ,v2+I δ,v3+I δ,v4+I δ,v5+I δ},MatrixDirectory->work,WorkingPrecision->150,ChopPrecision->100,UseMobius->True,UsePade->True};
 loaded=LoadConfiguration[config];assert[!FailureQ[loaded],loaded];
 boundary=PrepareBoundaryConditions[Get[FileNameJoin[{data,"XB_Boundary_values_X0.txt"}]],{v1->3,v2->-1,v3->1,v4->1,v5->-1}];
 result=TransportTo[boundary,{v1->4,v2->-113/47,v3->281/149,v4->349/257,v5->-863/541}];assert[!FailureQ[result],result];
 reference=Get[FileNameJoin[{data,"XB_Boundary_values_X1.txt"}]];discrepancy=Max[Abs[result[[2]]-reference]];assert[discrepancy<10^-28,discrepancy];,
 Print["Unknown original family"];Quit[1]]];
Export[FileNameJoin[{work,"request.json"}],DiffExp`Private`$lastRequest,"RawJSON"];
report=<|"schema"->"DiffExp.OriginalWrapperValidation/v1","family"->(family<>If[highPrecision,"-high",""]),"status"->"passed","max_discrepancy"->ToString[discrepancy,InputForm],"max_estimated_error"->ToString[Max[result[[3]]],InputForm],"timings"->Join[DiffExpLastTimings[],<|"test_total_seconds"->AbsoluteTime[]-t0|>]|>;
Export[FileNameJoin[{work,"result.json"}],report,"RawJSON"];
Print[ExportString[report,"RawJSON"]];Quit[0];
