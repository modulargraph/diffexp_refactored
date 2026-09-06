import subprocess,time,json,pathlib,statistics,decimal,re,hashlib
import argparse
a=argparse.ArgumentParser();a.add_argument('--fire',required=True);a.add_argument('--output',type=pathlib.Path,required=True);a.add_argument('--repeats',type=int,default=3);args=a.parse_args()
if not 1<=args.repeats<=20:a.error('--repeats must be 1..20')
root=pathlib.Path(__file__).resolve().parents[3];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);exe=root/'build/diffexp';fire=args.fire
report={'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'scope':'sequential cold caches, identical complete FT requests, no numerical checkpoint cache; 61-bit native reconstruction versus symbolic FIRE7','cases':[]}
def midpoint(x):
 x=x.strip('[] ')
 if x.startswith('+/-'):return decimal.Decimal(0)
 return decimal.Decimal(x.split('+/-')[0].strip())
def compare(a,b):
 worst=decimal.Decimal(0);decimal.getcontext().prec=100
 assert len(a['coefficients'])==len(b['coefficients'])
 for r in range(len(a['coefficients'])):
  for e in range(min(a['epsilon_low'],b['epsilon_low']),a['epsilon_high']+1):
   for part in ['real','imaginary']:
    x=midpoint(a['coefficients'][r][e-a['epsilon_low']][part]) if e>=a['epsilon_low'] else 0
    y=midpoint(b['coefficients'][r][e-b['epsilon_low']][part]) if e>=b['epsilon_low'] else 0
    worst=max(worst,abs(x-y))
 return str(worst)
for family in ['sunrise','box','box_triangle']:
 repeats=args.repeats
 case={'family':family,'runs':[]};values={}
 for rep in range(repeats):
  values={}
  for provider in (['ibp-solver','fire'] if rep%2==0 else ['fire','ibp-solver']):
   label=f'{family}-{provider}-{rep}';cmd=[str(exe),'ft',str(root/f'examples/feynman/{family}.json'),'--ibp-provider',provider,'--cache',str(out/(label+'-cache')),'--total-seconds','90','--level-seconds','70','--no-numerical-cache','--json']
   if provider=='fire':cmd+=['--fire',fire]
   start=time.perf_counter()
   with (out/(label+'.json')).open('w') as stdout,(out/(label+'.log')).open('w') as stderr:
    try:p=subprocess.run(cmd,stdout=stdout,stderr=stderr,timeout=150);code=p.returncode
    except subprocess.TimeoutExpired:code=124
   wall=time.perf_counter()-start
   try:v=json.loads((out/(label+'.json')).read_text())
   except Exception:v={}
   r={'provider':provider,'repeat':rep,'process_seconds':wall,'returncode':code,'status':v.get('status','timeout'),'timings':v.get('timings'),'ibp_statistics':v.get('ibp_statistics'),'message':v.get('message')};case['runs'].append(r)
   print(family,provider,rep,json.dumps(r),flush=True)
   if code==0:values[provider]=v
  if len(values)==2:case['maximum_coefficient_midpoint_difference']=str(max(decimal.Decimal(case.get('maximum_coefficient_midpoint_difference','0')),decimal.Decimal(compare(values['ibp-solver'],values['fire']))))
 report['cases'].append(case);(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')

for case in report['cases']:
 if any(run['returncode']!=0 for run in case['runs']):raise SystemExit('A complete FT run failed; see summary.json')
 if decimal.Decimal(case['maximum_coefficient_midpoint_difference'])>decimal.Decimal('1e-30'):raise SystemExit('Coefficient agreement failed; see summary.json')
