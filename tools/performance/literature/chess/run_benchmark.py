#!/usr/bin/env python3
"""Reproduce the standalone spectral experiment, with at most 600s per solve."""
import argparse,json,os,pathlib,shutil,statistics,subprocess,time
from fetch import fetch
HERE=pathlib.Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('work_dir',type=pathlib.Path)
 p.add_argument('--kernel',default=os.environ.get('WOLFRAM_KERNEL','/Applications/Wolfram.app/Contents/MacOS/WolframKernel'))
 p.add_argument('--repo',type=pathlib.Path,default=HERE.parents[3]);p.add_argument('--prefix',default='/opt/homebrew')
 p.add_argument('--cxx',default=os.environ.get('CXX','c++'));p.add_argument('--repeats',type=int,default=3)
 p.add_argument('--include-native',action='store_true');p.add_argument('--skip-fetch',action='store_true')
 a=p.parse_args();root=a.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
 if not a.skip_fetch:fetch(root)
 env=dict(os.environ,CHESS_BENCHMARK_DIR=str(root),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 def run(name,cmd):
  start=time.monotonic()
  with (root/(name+'.stdout')).open('w') as out,(root/(name+'.stderr')).open('w') as err:
   try:status=subprocess.run([str(x) for x in cmd],cwd=root,env=env,stdout=out,stderr=err,timeout=600).returncode
   except subprocess.TimeoutExpired:status='timeout600'
  (root/(name+'-process.json')).write_text(json.dumps({'command':[str(x) for x in cmd],'status':status,'wall_seconds':time.monotonic()-start},indent=2))
  if status!=0:raise RuntimeError(f'{name} failed: {status}; inspect {root/(name+".stderr")}')
 run('conversion',[a.kernel,'-script',HERE/'convert.wls'])
 run('build',[a.cxx,'-O3','-DNDEBUG','-std=c++20','-I'+str(a.repo.resolve()/'include'),'-I'+a.prefix+'/include',HERE/'spectral.cpp','-L'+a.prefix+'/lib','-lflint','-lgmp','-lmpfr','-lboost_json','-lboost_container','-o',root/'spectral'])
 req=json.loads((root/'dp-shared-request.json').read_text())
 for entry in req['entries']:entry['variable']='form'
 for n,goal in [(16,20),(32,40)]:
  req.update(taylor_order=n,accuracy_goal=goal)
  (root/f'dp-form-n{n}-g{goal}.json').write_text(json.dumps(req,separators=(',',':')))
 records=[]
 for rep in range(1,a.repeats+1):
  for n,goal in [(16,20),(32,40)]:
   run(f'spectral{n}-r{rep}',[root/'spectral',root/f'dp-form-n{n}-g{goal}.json',n])
   run(f'chess{n}-r{rep}',[a.kernel,'-script',HERE/'run_chess.wls',n])
   shutil.copy2(root/f'chess-n{n}.json',root/f'chess{n}-r{rep}.json')
   shutil.copy2(root/f'chess-n{n}-values.wl',root/f'chess{n}-r{rep}-values.wl')
 if a.include_native:
  for n,goal in [(16,20),(32,40)]:run(f'native{n}',[a.repo.resolve()/'build/diffexp','transport',root/f'dp-form-n{n}-g{goal}.json'])
 run('validation',[a.kernel,'-script',HERE/'validate.wls'])
 for n in [16,32]:
  for method in ['chess','spectral']:
   samples=[]
   for rep in range(1,a.repeats+1):
    name=f'{method}{n}-r{rep}';result=json.loads((root/(name+('.json' if method=='chess' else '.stdout'))).read_text())
    times={'preparation_seconds':result['preparation_seconds'],'numerical_seconds':result['seconds'],'total_seconds':result['preparation_seconds']+result['seconds']} if method=='chess' else result['timings']
    samples.append({'timings':times,'validation':json.loads((root/(name+'-validation.json')).read_text()),'process':json.loads((root/(name+'-process.json')).read_text())})
   records.append({'method':method,'nodes':n,'median_seconds':{key:statistics.median(s['timings'][key] for s in samples) for key in ['preparation_seconds','numerical_seconds','total_seconds']},'samples':samples})
 (root/'summary.json').write_text(json.dumps({'source_commit':'4a701fc1332f29f6237d14427336e60615b966e3','conversion':json.loads((root/'conversion-timing.json').read_text()),'records':records},indent=2))
 print(root/'summary.json')
if __name__=='__main__':main()
