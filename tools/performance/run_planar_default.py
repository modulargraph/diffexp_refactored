#!/usr/bin/env python3
"""Bounded full-path comparison on the original 13-master one-loop family."""
import argparse, gzip, hashlib, json, pathlib, subprocess
from decimal import Decimal, localcontext
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--executable',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--orders',type=int,nargs='+',default=[40,200],choices=[40,200])
a=p.parse_args()
if len(a.orders)>2 or len(set(a.orders))!=len(a.orders): p.error("Choose each order at most once")
out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
source=pathlib.Path(__file__).resolve().parent/'planar-one-loop-request.json.gz'
raw=gzip.decompress(source.read_bytes());request=json.loads(raw)
summary={'family':'original planar 1loop','request_sha256':hashlib.sha256(raw).hexdigest(),'comparisons':[]}
for order in a.orders:
    observations={}
    for method in ['series','auto']:
        name=f'order{order}-{method}'
        r=dict(request,recurrence=method,taylor_order=order,working_bits=499,accuracy_goal=15 if order==40 else 40)
        request_file=out/(name+'-request.json');request_file.write_text(json.dumps(r))
        with (out/(name+'.json')).open('w') as stdout,(out/(name+'.log')).open('w') as stderr:
            subprocess.run([str(a.executable.resolve()),'transport',str(request_file)],stdout=stdout,stderr=stderr,check=True,timeout=180)
        observations[method]=json.loads((out/(name+'.json')).read_text())
    with localcontext() as context:
        context.prec=180;worst=Decimal(0)
        for row,reference in zip(observations['auto']['values'],observations['series']['values']):
            for value,ref in zip(row,reference):
                keys=('real_midpoint','imaginary_midpoint')
                delta=sum(abs(Decimal(value[k])-Decimal(ref[k])) for k in keys)
                scale=max(Decimal(1),sum(abs(Decimal(ref[k])) for k in keys))
                worst=max(worst,delta/scale)
        if worst >= Decimal(10)**(-r['accuracy_goal']): raise RuntimeError('endpoint comparison failed')
    summary['comparisons'].append({'order':order,'maximum_normalized_midpoint_difference':str(worst),
        'observations':{method:{key:res[key] for key in ['timings','charts','recurrence','omitted_tails_certified']} for method,res in observations.items()}})
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(order,summary['comparisons'][-1],flush=True)
