(* Original DiffExp notebook shuffle and matrix-input workflow. *)
If[!ValueQ[$MPLWorkingPrecision],$MPLWorkingPrecision=90];
If[!ValueQ[$MPLExpansionOrder],$MPLExpansionOrder=40];
ListShuffleProduct[u:{a_,xx___},v:{bb_,yy___},c___] := Join[ListShuffleProduct[{xx},v,c,a],ListShuffleProduct[u,{yy},c,bb]];
ListShuffleProduct[{xx___},{yy___},c___] := {{c,xx,yy}};
ShuffleG[G[a__,z_],G[bb__,z_]] := Total[G[Sequence@@#,z]&/@ListShuffleProduct[{a},{bb}]];
ClearAll[GFB];
GFB[a___,a1_,bb___,0,z_] /; !a1===0 && bb===0 := G[a,a1,bb,0,z]+G[a,a1,z]G[bb,0,z]-ShuffleG[G[a,a1,z],G[bb,0,z]]/.G->GFB;
GFB[a___,a1_,z_] /; !a1===0 := G[a,a1,z];
GFB[a___,z_] /; a===0 := Log[z]^Length[{a}]/Length[{a}]!;
GEvaluateFiniteBasepoint[a__,endpoint_] /; Im[endpoint]===0 := Module[{matrix,l=Length[{a}],boundary,result},
 matrix=Append[PadLeft[#,l+1,0]&/@DiagonalMatrix[1/(t-#)&/@{a}],ConstantArray[0,l+1]];
 DeleteFile[FileNames["*.m",folder]];Put[matrix,FileNameJoin[{folder,"dt_0.m"}]];
 LoadConfiguration[{MatrixDirectory->folder,EpsilonOrder->0,WorkingPrecision->$MPLWorkingPrecision,ChopPrecision->60,ExpansionOrder->$MPLExpansionOrder,Verbosity->0,
  DeltaPrescriptions->Table[t-ind+If[ind>0,-I,I] \[Delta],{ind,Select[{a},Im[#]===0&]}]}];
 boundary=PrepareBoundaryConditions[Append[ConstantArray[0+O[x]^(1/2),l],1+O[x]^(1/2)],{t->endpoint x}];
 result=TransportTo[boundary,{t->endpoint}];If[ListQ[result],result[[2,1,1]],result]];
GEvaluate[a__,endpoint_] /; Im[endpoint]===0 := Expand[GFB[a,endpoint]/.G[a1__,bb_]/;SameQ[a1]:>(Log[1-bb/First[{a1}]]^Length[{a1}]/Length[{a1}]!)/.G[args__]:>GEvaluateFiniteBasepoint[args]];
