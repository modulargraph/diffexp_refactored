#!/usr/bin/env python3
"""Sequential baseline/shadow tests, each compile/run bounded to 60 seconds."""
from pathlib import Path
import json,subprocess,time
root=Path(__file__).resolve().parents[2]
subprocess.run(['python3','tools/performance/make_guarded_shadow.py'],cwd=root,check=True,timeout=60)
results=[];start=time.monotonic()
for name in ['conditioned_adjoint','centered_adjoint','polynomial_transport','adjoint_transport','adjoint_scaling']:
 for mode in ['baseline','shadow']:
  if time.monotonic()-start>540:raise RuntimeError('overall budget exhausted')
  binary=f'/tmp/diffexp-{name}-{mode}'
  cmd=['/usr/bin/c++','-std=c++20','-O3','-DNDEBUG']
  cmd+=['-Itools/performance/'+mode]
  cmd+=['-Iinclude','-I/opt/homebrew/include',f'tests/test_{name}.cpp','-L/opt/homebrew/lib','-lflint','-lmpfr','-lgmp','-lboost_json','-o',binary]
  t=time.monotonic();compiled=subprocess.run(cmd,cwd=root,capture_output=True,text=True,timeout=60,check=True);compile_seconds=time.monotonic()-t
  t=time.monotonic();ran=subprocess.run([binary],cwd=root,capture_output=True,text=True,timeout=60);elapsed=time.monotonic()-t
  result=dict(test=name,mode=mode,compile_seconds=compile_seconds,run_seconds=elapsed,status='pass' if ran.returncode==0 else 'fail',stdout=ran.stdout.strip(),stderr=ran.stderr.strip());results.append(result)
  print(json.dumps(result),flush=True)
  (root/'docs/validation/performance-guarded-dot.json').write_text(json.dumps(dict(status='pass' if all(r['status']=='pass' for r in results) else 'fail',elapsed_seconds=time.monotonic()-start,tests=results),indent=2)+'\n')
  if ran.returncode:raise RuntimeError('test failed')
