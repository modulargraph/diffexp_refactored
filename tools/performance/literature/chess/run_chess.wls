$HistoryLength=0; SetDirectory[Environment["CHESS_BENCHMARK_DIR"]];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->1];
nodes=ToExpression[Last[$CommandLine]];
{prep,unused}=AbsoluteTiming[Get["Chess.wl"];dLettersLine=Get["data/DP/dLettersLine.txt"];Atilde=SparseArray[Get["data/DP/Atilde.txt"]];CHESSAtildeLinearData[];bc=N[Get["boundary.wl"],60];ref=Get["reference.wl"]];
{elapsed,result}=AbsoluteTiming[TimeConstrained[SpectralPropagate[nAt,bc,{0,1},"Nodes"->nodes,"Precision"->60,"WorkingPrecisionA"->80,"ParallelEvaluation"->False,"ParallelKernels"->1,"RegularizedEndpoints"->{},"SolveStrategy"->"SequentialEpsilon"],590,$Aborted]];
If[result===$Aborted||result===$Failed,Export["chess-n"<>ToString[nodes]<>".json",<|"status"->ToString[result],"preparation_seconds"->prep,"seconds"->elapsed|>,"RawJSON"];Quit[2]];
final=CHESSFinalState[result];
Export["chess-n"<>ToString[nodes]<>".json",<|"status"->"ok","nodes"->nodes,"precision"->60,"working_precision_A"->80,"preparation_seconds"->prep,"seconds"->elapsed,"max_absolute_error"->ToString[Max[Abs[Flatten[final-ref]]],InputForm],"max_scaled_error"->ToString[Max[Flatten[Abs[final-ref]/(1+Abs[ref])]],InputForm]|>,"RawJSON"];
Put[final,"chess-n"<>ToString[nodes]<>"-values.wl"];Print["nodes ",nodes," seconds ",elapsed];Quit[];
