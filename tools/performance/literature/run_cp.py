#!/usr/bin/env python3
"""Run Canko–Pozzoli RL1/RL2 through the native CLI and Mathematica wrapper.
Download/extract the authors' ancillary archive first; see docs/literature-benchmarks.md.
"""
import argparse,pathlib,subprocess,os,shutil,time,json,signal
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--kernel',required=True)
p.add_argument('--package',type=pathlib.Path,required=True)
p.add_argument('--original',type=pathlib.Path,required=True)
a=p.parse_args();source=pathlib.Path(__file__).resolve().parent
if a.output.exists():p.error('output must be a fresh directory')
a.output.mkdir(parents=True);out=a.output.resolve()
# Copy only equation/boundary files; analytic-solution archives are not required.
for f in ['RL1','RL2']:
 shutil.copytree(a.data/f/'MIs_DEs',out/'ancillary'/f/'MIs_DEs')
shutil.copytree(a.data/'GlobalDefinitions',out/'ancillary'/'GlobalDefinitions')
for name in ['prepare.wls','compare_cp.wls','analyze_cp.py','epsilon_sampling.py']:
 shutil.copy(source/name,out/name)
shutil.copy(a.original,out/'DiffExp-original.m')
env=dict(os.environ,DIFFEXP_HOME=str(a.package.resolve()),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
observations=[]
def run(label,command,stdout_name=None):
 start=time.monotonic()
 with open(out/(stdout_name or label+'.stdout'),'w') as stdout,open(out/(label+'.stderr'),'w') as stderr:
  child=subprocess.Popen(command,cwd=out,env=env,stdout=stdout,stderr=stderr,start_new_session=True)
  try: child.wait(timeout=600);status=child.returncode
  except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait();status='timeout'
 record={'label':label,'status':status,'wall_seconds':time.monotonic()-start,'command':command};observations.append(record)
 (out/'processes.json').write_text(json.dumps(observations,indent=2));print(record,flush=True)
 if status!=0:raise RuntimeError(label+' failed; inspect preserved logs')
run('prepare',[a.kernel,'-noprompt','-script',str(out/'prepare.wls')])
for family in ['RL1','RL2']:
 for mode in ['original','wrapper']:
  run(family+'-'+mode,[a.kernel,'-noprompt','-script',str(out/'compare_cp.wls'),mode,family])
 run(family+'-native',[str(a.package.resolve()/'build/diffexp'),'transport',str(out/(family+'-request.json'))],family+'-native.json')
run('analysis',['python3',str(out/'analyze_cp.py')])
