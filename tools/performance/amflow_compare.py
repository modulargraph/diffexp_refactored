#!/usr/bin/env python3
"""Three bounded, matched FT benchmarks; see docs/amflow-comparison.md.
Requires AMFLOW_HOME, WOLFRAM_KERNEL, FIRE7_BINARY. Builds/installs nothing.
AMFLOW_BENCH_ROOT defaults to /tmp/diffexp-amflow-comparison.
DIFFEXP_HOME defaults to this repository. Run examples sequentially.
"""
import argparse, json, os, signal, subprocess, threading, time
from pathlib import Path
LITERED_SETUP=os.environ.get("AMFLOW_LITERED_SETUP", "")
ROOT=Path(os.environ.get("AMFLOW_BENCH_ROOT", "/tmp/diffexp-amflow-comparison")).resolve()
AMFLOW=Path(os.environ["AMFLOW_HOME"]).resolve()
DIFFEXP=Path(os.environ.get("DIFFEXP_HOME",str(Path(__file__).resolve().parents[2]))).resolve()
KERNEL=os.environ["WOLFRAM_KERNEL"]
FIRE=os.environ["FIRE7_BINARY"]

def prepare(case):
    c=json.loads((DIFFEXP/'examples/feynman'/f'{case}.json').read_text())
    c['epsilon_order']=2 if case=='bubble' else 0
    c['name']='benchmark_'+case
    return c

def wl_driver(config, folder, mode, label, digits):
    # AMFlow requires a complete scalar-product family, including zero-power ISPs.
    # Its physical denominators are q^2-m^2; DiffExp uses the opposite sign.
    return f'''SetDirectory[{json.dumps(str(folder))}];
start=AbsoluteTime[];
Get[{json.dumps(str(AMFLOW/'AMFlow.m'))}];
SetReductionOptions["IBPReducer"->"FIRE+LiteRed"];
{'FIRE`Private`$LiteRedPath='+json.dumps(LITERED_SETUP)+';' if LITERED_SETUP else ''}
config=Import["family.json","RawJSON"];
loops=Table[Symbol["k"<>ToString[i]],{{i,config["loops"]}}]; legs=Table[Symbol["p"<>ToString[i]],{{i,Length[config["external_gram"]]}}];
read[s_String]:=ToExpression[s,InputForm];
props=Map[(Dot[read/@#["loop_coefficients"],loops]+Dot[read/@#["external_coefficients"],legs])^2-read[#["mass_squared"]]&,config["propagators"]];
indices=ConstantArray[1,Length[props]];
If[config["loops"]===2 && Length[legs]===1,
 props=Join[props,Map[(#+First[legs])^2&,loops]];indices=Join[indices,{{0,0}}]];
AMFlowInfo["Family"]=bench;
AMFlowInfo["Loop"]=loops; AMFlowInfo["Leg"]=legs;
AMFlowInfo["Conservation"]={{}};
AMFlowInfo["Replacement"]=Flatten[Table[legs[[i]] legs[[j]]->read[config["external_gram"][[i,j]]],{{i,Length[legs]}},{{j,i,Length[legs]}}]];
AMFlowInfo["Propagator"]=props;AMFlowInfo["Numeric"]={{}};
AMFlowInfo["NThread"]=1;
SetAMFOptions["RecursionMode"->{json.dumps(mode)},"DESolver"->"CPP","D0"->config["dimension_at_epsilon_zero"],"UseCache"->True,"CacheName"->"cache","ChopPre"->{digits}];
SeedRandom[1729];
benchmarkTarget=j[bench,Sequence@@indices];
order=config["epsilon_order"]+2 config["loops"];
numconfig=GenerateNumericalConfig[{digits},order];
{{seconds,sol}}=AbsoluteTiming[CheckAbort[SolveIntegrals[{{benchmarkTarget}},{digits},order],$Aborted]];
Put[sol,{json.dumps(label+".raw.wl")}];
If[sol===$Aborted || !ListQ[sol] || !MemberQ[First/@sol,benchmarkTarget],Print["BENCHMARK FAILED"];Quit[1]];
value=benchmarkTarget/.sol;
coefs=Table[(-1)^Total[indices] Coefficient[value,eps,n],{{n,-2 config["loops"],config["epsilon_order"]}}];
If[!VectorQ[coefs,NumberQ],Print["NONNUMERICAL COEFFICIENTS: ",InputForm[coefs]];Quit[1]];
report=<|"engine"->"AMFlow 2.0","recursion_mode"->{json.dumps(mode)},"case"->config["name"],"status"->"completed","goal_digits"->{digits},"solve_seconds"->seconds,"kernel_session_seconds"->AbsoluteTime[]-start,
 "epsilon_low"->(-2 config["loops"]),"epsilon_high"->config["epsilon_order"],"normalized_coefficients"->(ToString[#,InputForm]&/@coefs),"normalization"->"Multiply AMFlow by (-1)^(sum of propagator powers); no EulerGamma prefactor on either side.",
 "numerical_configuration"-><|"epsilon_samples"->Length[numconfig[[1]]],"working_decimal_digits"->numconfig[[2]],"taylor_order"->numconfig[[3]],"extra_order"->numconfig[[4]]|>,"configured_threads"->1,"litered_setup_override"->{json.dumps(LITERED_SETUP)}|>;
Export[{json.dumps(label+'.result.json')},report,"RawJSON"];Put[sol,{json.dumps(label+'.result.wl')}];Print["BENCHMARK RESULT: ",ExportString[report,"RawJSON"]];Quit[0];
'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('engine',choices=['amflow','diffexp']);ap.add_argument('case',choices=['bubble','sunrise','box']);ap.add_argument('--mode',default='FT',choices=['FT','AMF']);ap.add_argument('--label',default='cold');ap.add_argument('--group',default='comparison');ap.add_argument('--digits',type=int,default=20);ap.add_argument('--timeout',type=int,default=600);ap.add_argument('--endpoint-order',type=int);ap.add_argument('--ordinary-order',type=int);ap.add_argument('--working-bits',type=int)
    a=ap.parse_args();assert 0<a.timeout<=600
    folder=ROOT/'runs'/a.group/(a.engine+'-'+a.mode.lower())/a.case;folder.mkdir(parents=True,exist_ok=True)
    config=prepare(a.case)
    for key in ('endpoint_order','ordinary_order','working_bits'):
        if getattr(a,key) is not None:config['numerical'][key]=getattr(a,key)
    if (folder/(a.label+'.process.json')).exists():raise SystemExit('Refusing to overwrite a completed observation; choose a new label/group.')
    if a.label=='cold' and (folder/'cache').exists():raise SystemExit('Cold run requires an unused preparation cache; choose a new group.')
    (folder/'family.json').write_text(json.dumps(config,indent=2)+'\n')
    if a.engine=='amflow':
        script=folder/'run.wl';script.write_text(wl_driver(config,folder,a.mode,a.label,a.digits));cmd=[KERNEL,'-noprompt','-script',str(script)]
    else: cmd=[str(DIFFEXP/'build/diffexp'),'ft',str(folder/'family.json'),'--fire',FIRE,'--cache',str(folder/'cache'),'--json','--total-seconds','550','--level-seconds','300','--fire-seconds','120']
    env=dict(os.environ);env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1';env['VECLIB_MAXIMUM_THREADS']='1'
    start=time.monotonic();code=None
    with (folder/(a.label+'.stdout.log')).open('w') as out,(folder/(a.label+'.stderr.log')).open('w') as err:
        proc=subprocess.Popen(cmd,cwd=folder,env=env,stdout=out,stderr=err,start_new_session=True)
        timed_out=threading.Event()
        def stop():
            if proc.poll() is None:
                timed_out.set()
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
        watchdog=threading.Timer(a.timeout,stop);watchdog.daemon=True;watchdog.start()
        code=proc.wait();watchdog.cancel()
        if timed_out.is_set():code='timeout'
    result={'engine':a.engine,'mode':a.mode,'case':a.case,'label':a.label,'exit':code,'wall_seconds':time.monotonic()-start,'time_limit_seconds':a.timeout,'command':cmd,'root':str(folder)}
    (folder/(a.label+'.process.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if a.engine=='diffexp' and code==0:
        data=json.loads((folder/(a.label+'.stdout.log')).read_text());(folder/(a.label+'.result.json')).write_text(json.dumps(data,indent=2)+'\n')
    raise SystemExit(0 if code==0 and (folder/(a.label+'.result.json')).exists() else 1)
if __name__=='__main__':main()
