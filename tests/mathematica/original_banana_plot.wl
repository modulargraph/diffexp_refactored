(* Saved original-banana functions: reference, nearby replay and plotting. *)
root=ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName],"..",".."}]];
Get[FileNameJoin[{root,"DiffExp.m"}]];
assert[c_,m_]:=If[!TrueQ[c],Print["FAIL: ",m];Quit[1]];
saved=Get[FileNameJoin[{root,"build-reference","banana-wrapper-equal-saved.m"}]];
assert[ListQ[saved],"Run original_banana_routes.wl first."];
f=ToPiecewise[saved,True];assert[MatrixQ[f],f];
{referenceSeconds,at20}=AbsoluteTiming[Table[f[[i,k]][20],{i,4},{k,5}]];
reference=Get[FileNameJoin[{root,"examples","original","Reference","BananaEqualMassAt20.m"}]];
referenceError=Max[Abs[at20-reference]];assert[referenceError<10^-9,referenceError];
{nearbySeconds,nearby}=AbsoluteTiming[Table[f[[i,k]][20.1],{i,4},{k,5}]];
meta=saved[[2,1]];LoadConfiguration[meta["configuration"]];
{directSeconds,direct}=AbsoluteTiming[TransportTo[meta["boundary"],meta["line"],20.1]];
assert[!FailureQ[direct],direct];nearbyError=Max[Abs[nearby-direct[[2]]]];assert[nearbyError<10^-20,nearbyError];
{plotSeconds,plot}=AbsoluteTiming[ReImPlot[Evaluate[f[[3,4]][x]],{x,20,20.2},WorkingPrecision->100,PlotPoints->3,MaxRecursion->0]];
assert[Head[plot]===Graphics && !FreeQ[plot,_Line],"Original ReImPlot workflow did not produce curves."];
report=<|"schema"->"DiffExp.OriginalWrapperValidation/v1","family"->"banana-plot","status"->"passed","reference_discrepancy"->ToString[referenceError,InputForm],"nearby_replay_discrepancy"->ToString[nearbyError,InputForm],"reference_evaluation_seconds"->referenceSeconds,"nearby_evaluation_seconds"->nearbySeconds,"direct_evaluation_seconds"->directSeconds,"plot_seconds"->plotSeconds|>;
Export[FileNameJoin[{root,"build-reference","banana-wrapper-plot.json"}],report,"RawJSON"];Print[ExportString[report,"RawJSON"]];Quit[0];
