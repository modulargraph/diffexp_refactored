root=DirectoryName[$InputFileName]; SetDirectory[root];
str[v_]:=ToString[v,InputForm];
a=Get["ancillary/GlobalDefinitions/alphabetxyz.m"]/.x->u;
Do[
 m=Get["ancillary/"<>fam<>"/MIs_DEs/dAtildexyz.m"];
 bc=Get["ancillary/"<>fam<>"/MIs_DEs/boundariesxyz.m"];
 entries=Flatten[Table[With[{c=Coefficient[m[[i,j]],dlog[W[k]]]},If[c===0,Nothing,<|"row"->i-1,"column"->j-1,"epsilon"->1,"variable"->"dlog","expression"->str[W[k]/.a],"coefficient"->str[c]|>]],{i,Length[m]},{j,Length[m]},{k,15}]];
 req=<|"schema"->"DiffExp.Transport/v1","dimension"->Length[m],"paths"-><|"u"->"3/2+x/6","y"->"1/5-6*x/55","z"->"1/2-5*x/14"|>,"entries"->entries,"boundary"->Table[str[Coefficient[bc[[i]],eps,k]],{i,Length[m]},{k,0,6}],"epsilon_order"->6,"taylor_order"->80,"working_bits"->384,"accuracy_goal"->40,"division_order"->4|>;
 Export[fam<>"-request.json",req,"RawJSON"];
 Print[fam," dimension ",Length[m]," entries ",Length[entries]," boundary precision ",Min[Cases[bc,v_Real:>Precision[v],Infinity]]];
 ,{fam,{"RL1","RL2"}}];
