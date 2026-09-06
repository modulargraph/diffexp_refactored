root=DirectoryName[$InputFileName];SetDirectory[root];
If[Environment["DIFFEXP_HOME"]===$Failed,Print["Set DIFFEXP_HOME to the package checkout"];Exit[1]];
args=Take[$CommandLine,-2]; mode=If[Length[args]>0,args[[1]],"original"]; fam=If[Length[args]>1,args[[2]],"RL1"];
Get[If[mode==="original",root<>"DiffExp-original.m",FileNameJoin[{Replace[Environment["DIFFEXP_HOME"],$Failed->"REQUIRE_DIFFEXP_HOME"],"DiffExp.m"}]]];
folder=root<>"prepared-"<>fam;If[!DirectoryQ[folder],CreateDirectory[folder]];
{prep,ignored}=AbsoluteTiming[
 m=Get[root<>"ancillary/"<>fam<>"/MIs_DEs/dAtildexyz.m"]; a=Get[root<>"ancillary/GlobalDefinitions/alphabetxyz.m"];
 Do[Put[m/.dlog[v_]:>D[Log[v/.a],var],folder<>"/d"<>ToString[var]<>"_1.m"],{var,{x,y,z}}];
 bc=Get[root<>"ancillary/"<>fam<>"/MIs_DEs/boundariesxyz.m"];
 config=LoadConfiguration[{Variables->{x,y,z},DeltaPrescriptions->{},MatrixDirectory->folder,UseMobius->False,UsePade->False,EpsilonOrder->6,AccuracyGoal->40,ExpansionOrder->50,WorkingPrecision->115,ChopPrecision->100,LineParameter->t}];
 initial=PrepareBoundaryConditions[bc,<|x->3/2,y->1/5,z->1/2|>];];
Print["PREPARATION ",prep];
{elapsed,result}=AbsoluteTiming[TimeConstrained[TransportTo[initial,<|x->5/3,y->1/11,z->1/7|>],570,$Aborted]];
If[result===$Aborted||FailureQ[result],Print[result];Exit[1]];
Put[result,root<>fam<>"-"<>mode<>"-result.m"];
str[q_]:=ToString[q,InputForm];
Export[root<>fam<>"-"<>mode<>"-result.json",<|"seconds"->elapsed,"preparation_seconds"->prep,"values"->Map[str,result[[2]],{2}]|>,"RawJSON"];
Print["COMPLETED ",mode," ",fam," ",elapsed];
