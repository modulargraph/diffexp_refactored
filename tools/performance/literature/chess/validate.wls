$HistoryLength=0;SetDirectory[Environment["CHESS_BENCHMARK_DIR"]];
ref=Get["reference.wl"];
num[s_]:=ToExpression[StringReplace[s,{"e"->"*^","E"->"*^"}]];
metrics[v_]:=<|"components"->Dimensions[v],"coefficient_count"->Length[Flatten[v]],"max_absolute_error"->ToString[Max[Flatten[Abs[v-ref]]],InputForm],"max_normalized_error"->ToString[Max[Flatten[Abs[v-ref]/Map[Max[1,Abs[#]]&,ref,{2}]]],InputForm],"max_scaled_error"->ToString[Max[Flatten[Abs[v-ref]/(1+Abs[ref])]],InputForm]|>;
Do[r=Import[file,"RawJSON"];v=Map[num[#["real_midpoint"]]+I*num[#["imaginary_midpoint"]]&,r["values"],{2}];Export[StringDrop[file,-7]<>"-validation.json",metrics[v],"RawJSON"],{file,Join[FileNames["spectral*-r*.stdout"],FileNames["native*.stdout"]]}];
Do[v=Get[file];Export[StringDrop[file,-10]<>"-validation.json",metrics[v],"RawJSON"],{file,FileNames["chess*-r*-values.wl"]}];
Quit[];
