#!/usr/bin/env python3
"""Bounded, sequential first-chart experiment; never runs full transport."""
import argparse, gzip, hashlib, json, os, pathlib, subprocess
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repo', type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[3])
p.add_argument('--prefix', type=pathlib.Path, default=pathlib.Path('/opt/homebrew'))
p.add_argument('--orders', type=int, nargs='+', default=[115,230,460])
a=p.parse_args()
if len(a.orders)>4 or any(n<8 or n>460 for n in a.orders):
    p.error('Choose at most four orders between 8 and 460')
src=pathlib.Path(__file__).resolve().parent
out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
def run(cmd,name,timeout=180,env=None):
    with (out/(name+'.json')).open('w') as stdout, (out/(name+'.log')).open('w') as stderr:
        subprocess.run([str(x) for x in cmd],stdout=stdout,stderr=stderr,check=True,timeout=timeout,env=env)
data=gzip.decompress((src/'zzz-high-request.json.gz').read_bytes())
assert hashlib.sha256(data).hexdigest()=='f59f627573206c1d8d3362221317eecb3a72c2c647b8709e8d93bb7249d0479d'
(out/'request.json').write_bytes(data)
for name in ['test','probe','benchmark']:
    run([os.environ.get('CXX','c++'),'-std=c++20','-O3','-DNDEBUG','-I'+str(a.repo/'include'),'-I'+str(a.prefix/'include'),src/(name+'.cpp'),'-L'+str(a.prefix/'lib'),'-lflint','-lmpfr','-lgmp','-lboost_json','-o',out/name],'build-'+name,120)
run([out/'test'],'tests',60)
run([out/'probe',out/'request.json'],'decomposition')
for n in a.orders:
    run([out/'benchmark',out/'request.json',out/'decomposition.json',n],'order'+str(n))
print('Results and full coefficient arrays saved in',out)

if 230 in a.orders and 460 in a.orders:
    from decimal import Decimal, localcontext
    low=json.loads((out/'order230.json').read_text())
    high=json.loads((out/'order460.json').read_text())
    def distance(left,right):
        with localcontext() as ctx:
            ctx.prec=310
            worst=Decimal(0)
            for lr,rr in zip(left,right):
                for l,r in zip(lr,rr):
                    keys=('real_midpoint','imaginary_midpoint')
                    delta=sum(abs(Decimal(l[k])-Decimal(r[k])) for k in keys)
                    scale=max(Decimal(1),sum(abs(Decimal(r[k])) for k in keys))
                    worst=max(worst,delta/scale)
            return str(worst)
    comparisons={
        'candidate230_vs_candidate460':distance(low['candidate_values'],high['candidate_values']),
        'baseline230_vs_baseline460':distance(low['baseline_values'],high['baseline_values']),
        'candidate230_vs_baseline460':distance(low['candidate_values'],high['baseline_values']),
        'candidate460_vs_baseline460':distance(high['candidate_values'],high['baseline_values'])}
    (out/'convergence.json').write_text(json.dumps(comparisons,indent=2)+'\n')
