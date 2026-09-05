#!/usr/bin/env python3
"""Run finite chart experiments sequentially; every child has a 60-second limit."""
import argparse,json,platform,statistics,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary');p.add_argument('--output',default='docs/validation/performance-charts.json');a=p.parse_args()
start=time.monotonic();records=[]
methods=['baseline','inplace_horner','raw_scalar','dot','precise_dot']
for trial in range(3):
    # Rotate invocation order to reduce systematic thermal/order bias.
    ordered=methods[trial:]+methods[:trial]
    for dense in (0,1):
        for method in ordered:
            if time.monotonic()-start>540:raise RuntimeError('overall benchmark budget exhausted')
            r=subprocess.run([a.binary,method,str(dense),'5'],text=True,capture_output=True,timeout=60,check=True)
            record=json.loads(r.stdout);record['trial']=trial+1;records.append(record)
            print(f"trial={trial+1} dense={dense} method={method} warm_cpu={record['warm_cpu_seconds']:.6f}",flush=True)
summary=[]
for dense in (False,True):
    baseline=statistics.median(r['warm_cpu_seconds'] for r in records if r['dense']==dense and r['method']=='baseline')
    for method in methods:
        rows=[r for r in records if r['dense']==dense and r['method']==method]
        cpu=statistics.median(r['warm_cpu_seconds'] for r in rows)
        summary.append(dict(dense=dense,method=method,warm_median_cpu_seconds=cpu,warm_median_wall_seconds=statistics.median(r['warm_wall_seconds'] for r in rows),first_median_wall_seconds=statistics.median(r['first_wall_seconds'] for r in rows),speedup_cpu=baseline/cpu))
result=dict(status='pass',recorded_at=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),machine=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),platform=platform.platform(),compiler=subprocess.check_output(['/usr/bin/c++','--version'],text=True).splitlines()[0],flags='-std=c++20 -O3 -DNDEBUG',flint_version='3.4.0',precision_bits=384,taylor_order=80,epsilon_high=5,dimension=41,entries=124,observables=4,chart_workers=1,first_definition='First chart after input preparation in a fresh process; no chart warm-up; not an OS-cache-cold claim.',warm_definition='Mean of five chart calls after one initial chart; summaries are medians over three fresh processes.',scope='Retained Taylor polynomial only; no omitted-tail certification or full-example speed claim.',elapsed_seconds=time.monotonic()-start,summary=summary,measurements=records)
Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
