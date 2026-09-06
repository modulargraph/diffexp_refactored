#!/usr/bin/env python3
"""Validate the three benchmark coefficient windows against independent references.
Usage: analyze_amflow.py BENCH_ROOT OUTPUT_JSON
Run references.wl first; use final/{amflow-ft,diffexp-ft}, corrected_v2/amflow-ft/sunrise.
"""
import decimal, json, os, re, sys
from pathlib import Path
D=decimal.Decimal;decimal.getcontext().prec=75
root=Path(sys.argv[1]).resolve();output=Path(sys.argv[2])
package=Path(os.environ.get('DIFFEXP_HOME', str(Path(__file__).resolve().parents[2]))).resolve()
fire=Path(os.environ['FIRE7_BINARY']).resolve().parent.parent if 'FIRE7_BINARY' in os.environ else None
def wl(s):
    s=re.sub(r'`{1,2}[0-9.]*','',s).replace('*^','e').strip()
    return D(s) # These Euclidean cases have real AMFlow coefficients after Chop.
def ball(s):
    if s=='0':return D(0)
    s=s.strip('[] ')
    if s.startswith('+/-'):return D(0)
    return D(s.split('+/-')[0].strip())
def read(p):return json.loads(p.read_text())
def sanitize(x):
    if isinstance(x,str):
        x=x.replace('\\/','/').replace(str(root),'<benchmark-root>').replace(str(package),'<diffexp>')
        return x.replace(str(fire),'<fire>') if fire else x
    if isinstance(x,list):return [sanitize(i) for i in x]
    if isinstance(x,dict):return {k:sanitize(v) for k,v in x.items()}
    return x
refs=read(root/'references.json')
report={'schema':'DiffExp.AMFlowComparison/v1','target_digits':20,'acceptance':'max |computed-reference|/max(1,|reference|) <= 1e-20; imaginary native midpoints checked too; general FT tails are not certified','cases':[]}
for case,orders in [('bubble',[0,1,2]),('sunrise',[0]),('box',[-2,-1,0])]:
    entry={'case':case,'epsilon_orders':orders,'configuration':read(root/'runs/final/amflow-ft'/case/'family.json'),'observations':[]}
    reference=[wl(s) for s in refs[case]]
    for engine in ['diffexp','amflow']:
        group='corrected_v2' if case=='sunrise' and engine=='amflow' else 'final'
        folder=root/'runs'/group/(engine+'-ft')/case
        for label in ['cold','warm']:
            source_label=label
            if label=='warm' and not (folder/(label+'.result.json')).exists() and (folder/'warm_retry1.result.json').exists():source_label='warm_retry1'
            data=read(folder/(source_label+'.result.json'));process=read(folder/(source_label+'.process.json'))
            if process['exit']!=0:raise SystemExit(f'{case} {engine} {label} did not finish')
            low=data['epsilon_low']
            values=[(ball(data['coefficients'][0][o-low]['real']) if engine=='diffexp' else wl(data['normalized_coefficients'][o-low])) for o in orders]
            imag=[ball(data['coefficients'][0][o-low]['imaginary']) for o in orders] if engine=='diffexp' else [D(0)]*len(orders)
            errors=[max(abs(a-b),abs(c))/max(D(1),abs(b)) for a,b,c in zip(values,reference,imag)]
            maximum=max(errors)
            obs={'engine':engine,'label':label,'max_scaled_error':str(maximum),'accuracy_passed':maximum<=D('1e-20'),'process':process,'result':data}
            if engine=='amflow':
                log=(folder/(source_label+'.stdout.log')).read_text()
                times=re.findall(r'BlackBoxFTSingle: FT systems solved in ([0-9]+)s',log)
                if times:obs['numerical_phase_seconds_rounded_up']=int(times[-1])
            entry['observations'].append(obs)
            if maximum>D('1e-20'):raise SystemExit(f'Accuracy failed: {case} {engine} {label}: {maximum}')
    # Check all native negative coefficients vanish for the finite D=2 sunrise.
    if case=='sunrise':
        for label in ['cold','warm']:
            d=read(root/'runs/final/diffexp-ft/sunrise'/(label+'.result.json'))
            polemax=max(abs(ball(c[part])) for c in d['coefficients'][0][:-1] for part in ['real','imaginary'])
            if polemax>D('1e-20'):raise SystemExit('Native sunrise spurious pole')
    report['cases'].append(entry)
report['independent_references']=refs
report['stock_sunrise_failure']={'status':'failed_accuracy_validation','explanation':'Bundled LiteRed 1.83 marks fewer-propagator-than-loop sectors zero before its actual scaleless test. A merged FT quadratic is a counterexample. Corrected measurements use an isolated two-line dependency change.','process':read(root/'runs/final/amflow-ft/sunrise/cold.process.json'),'result':read(root/'runs/final/amflow-ft/sunrise/cold.result.json'),'diagnostic':read(root/'litered-diagnostic.json')}
report['excluded_failed_warm_startups']=[read(f) for f in root.glob('runs/corrected_v2/amflow-ft/sunrise/warm.process.json') if read(f)['exit']!=0]
report['litered_patch']=read(root/'litered-ft-setup/ft-compatibility-patch.json')
report['all_accepted_observations_passed']=True
output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(sanitize(report),indent=2)+'\n')
for c in report['cases']:
 print(c['case'])
 for o in c['observations']:print(' ',o['engine'],o['label'],round(o['process']['wall_seconds'],6),'s; error',o['max_scaled_error'])
