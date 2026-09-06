SetDirectory[DirectoryName[$InputFileName]];
start=AbsoluteTime[];
g=EulerGamma;
bubble=Table[NIntegrate[Switch[n,0,1,1,-g-Log[1+x(1-x)],2,((g+Log[1+x(1-x)])^2/2+Pi^2/12)]/(1+x(1-x)),{x,0,1},WorkingPrecision->65,AccuracyGoal->45,PrecisionGoal->45],{n,0,2}];
sunrise=4 NIntegrate[r BesselJ[0,r] BesselK[0,r]^3,{r,0,1,Infinity},WorkingPrecision->65,AccuracyGoal->42,PrecisionGoal->42];
ls=Log[1];lt=Log[1/3];
box=N[{12,-6(ls+lt+2g),3(2ls lt+2g(ls+lt)+2g^2-4Pi^2/3)},55];
report=<|"bubble"->(ToString[#,InputForm]&/@bubble),"sunrise"->{ToString[sunrise,InputForm]},"box"->(ToString[#,InputForm]&/@box),"seconds"->AbsoluteTime[]-start,"method"->"Independent Feynman-parameter quadrature for bubble, coordinate-space Bessel quadrature for sunrise, analytic massless-box Laurent coefficients; same positive scalar normalization, no EulerGamma exponential."|>;
Export["references.json",report,"RawJSON"];Print[ExportString[report,"RawJSON"]];Quit[];
