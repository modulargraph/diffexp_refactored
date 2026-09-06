"""Controlled transport-only epsilon sampling: same RL1 matrix and boundary polynomial.
This is not AMFlow's implementation and not a full FT benchmark.
"""
import json, pathlib, re, time, subprocess, signal, os, copy
import mpmath as mp
mp.mp.dps=120
root=pathlib.Path(__file__).resolve().parent; out=root/'epsilon-study';out.mkdir(exist_ok=True)
base=json.loads((root/'RL1-request.json').read_text()); high=2
base.update(epsilon_order=high,taylor_order=64,working_bits=384,accuracy_goal=40)
base['boundary']=[row[:high+1] for row in base['boundary']]
def val(s):
 s=re.sub(r'`{1,2}[+-]?[0-9.]+','',s).replace('*^','e').replace(' ','')
 terms=re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?(?:\*I)?',s)
 if ''.join(terms)!=s: raise ValueError(s)
 return sum((mp.mpc(0,mp.mpf(t[:-2])) if t.endswith('*I') else mp.mpc(mp.mpf(t))) for t in terms)
bc=[[val(s) for s in row] for row in base['boundary']]
def run(label,req):
 p=out/(label+'.request.json');p.write_text(json.dumps(req));st=time.monotonic()
 with (out/(label+'.json')).open('w') as stdout,(out/(label+'.log')).open('w') as stderr:
  proc=subprocess.Popen([str(pathlib.Path(os.environ['DIFFEXP_HOME'])/'build/diffexp'),'transport',str(p)],stdout=stdout,stderr=stderr,start_new_session=True)
  try:proc.wait(timeout=600)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();raise RuntimeError('timeout')
 if proc.returncode:raise RuntimeError(label+' failed')
 result=json.loads((out/(label+'.json')).read_text());return result,time.monotonic()-st
summary={'description':__doc__,'dimension':base['dimension'],'epsilon_high':high,'working_bits':384,'taylor_order':64,'goal_digits':40,'runs':[]}
for method in ['series','auto']:
 req=dict(base,recurrence=method);r,wall=run('direct-'+method,req)
 if method=='series': reference=[[mp.mpc(c['real_midpoint'],c['imaginary_midpoint']) for c in row] for row in r['values']]
 summary['runs'].append({'method':'direct-'+method,'wall_seconds':wall,'timings':r['timings']})
 print(summary['runs'][-1],flush=True)
for count in [12,16]:
 radius=mp.mpf('0.0001');samples=[];wall=0;phases=[]
 for j in range(count):
  z=radius*mp.exp(2j*mp.pi*j/count)
  # Exact rational approximations to the complex sampling nodes.
  den=10**110; zr=int(mp.nint(z.real*den));zi=int(mp.nint(z.imag*den));z=mp.mpc(mp.mpf(zr)/den,mp.mpf(zi)/den)
  zs=f'({zr}/{den}+I*({zi}/{den}))'
  req=copy.deepcopy(base);req['epsilon_order']=0;req['recurrence']='series'
  for e in req['entries']:e['epsilon']=0;e['coefficient']='('+e['coefficient']+')*'+zs
  req['boundary']=[[mp.nstr(v.real,110).replace('e','*^')+'+I*('+mp.nstr(v.imag,110).replace('e','*^')+')'] for v in [sum(c*z**k for k,c in enumerate(row)) for row in bc]]
  req['boundary_errors']=[['1*^-88'] for _ in bc]
  r,w=run(f'sample-{count}-{j}',req);wall+=w;phases.append(r['timings']);samples.append([mp.mpc(row[0]['real_midpoint'],row[0]['imaginary_midpoint']) for row in r['values']])
 result=[[sum(samples[j][i]*mp.exp(-2j*mp.pi*j*k/count) for j in range(count))/count/radius**k for k in range(high+1)] for i in range(len(bc))]
 error=max(abs(a-b)/max(1,abs(b)) for row,ref in zip(result,reference) for a,b in zip(row,ref))
 observation={'method':'cauchy-sampling','samples':count,'radius':str(radius),'wall_seconds':wall,'timings_sum':{k:sum(p[k] for p in phases) for k in phases[0]},'maximum_scaled_difference':str(error),'passes_40_digits':bool(error<mp.mpf('1e-40'))}
 summary['runs'].append(observation);print(observation,flush=True)
 (out/'summary.json').write_text(json.dumps(summary,indent=2))
(out/'summary.json').write_text(json.dumps(summary,indent=2))
