#!/usr/bin/env python3
"""Sequential generic FT examples with independent checks and process-group caps."""
import argparse, hashlib, json, os, pathlib, signal, subprocess, time
REPO=pathlib.Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--executable',type=pathlib.Path,default=REPO/'build/diffexp')
p.add_argument('--checker',type=pathlib.Path,default=REPO/'build/diffexp_check_ft_laptop')
p.add_argument('--fire',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
p.add_argument('--cases',nargs='+',choices=['pentagon_massive','box_triangle','banana_unequal','double_box_planar'],default=['pentagon_massive','box_triangle','banana_unequal','double_box_planar'])
p.add_argument('--timeout',type=int,default=600)
a=p.parse_args()
if not 1<=a.timeout<=600 or len(a.cases)!=len(set(a.cases)):p.error('Timeout must be 1..600 seconds; select each case once')
root=a.output.resolve()
if root.exists() and any(root.iterdir()):p.error('Use an empty output directory so preparation caches start cold')
root.mkdir(parents=True,exist_ok=True)
for key in ['executable','checker','fire']:
    path=getattr(a,key).resolve()
    if not path.is_file():p.error(f'{key} executable does not exist: {path}')
    setattr(a,key,path)
summary={'schema':'DiffExp.LaptopFT/v1','persistent_numerical_cache':False,'cold_definition':'Empty exact preparation cache; no OS cache flush','cases':[]}
def run(cmd,folder,label,limit):
    start=time.monotonic()
    with (folder/(label+'.json')).open('w') as stdout,(folder/(label+'.log')).open('w') as stderr:
        proc=subprocess.Popen([str(x) for x in cmd],stdout=stdout,stderr=stderr,start_new_session=True)
        try:
            code=proc.wait(timeout=limit);status='completed' if code==0 else 'failed'
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=3)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            code=proc.returncode;status='timeout'
    record={'status':status,'exit_code':code,'wall_seconds':time.monotonic()-start,'timeout_seconds':limit}
    (folder/(label+'-run.json')).write_text(json.dumps(record,indent=2)+'\n')
    return record
for case in a.cases:
    folder=root/case;folder.mkdir()
    family=json.loads((REPO/'examples/feynman'/f'{case}.json').read_text())
    family['name']='laptop_'+case
    config=folder/'family.json';config.write_text(json.dumps(family,indent=2)+'\n')
    record={'case':case,'configuration_sha256':hashlib.sha256(config.read_bytes()).hexdigest()}
    record['run']=run([a.executable,'ft',config,'--fire',a.fire,'--cache',folder/'cache','--total-seconds','480','--fire-threads','1','--fire-simplifier-threads','1','--no-numerical-cache','--json'],folder,'cold',a.timeout)
    if record['run']['status']=='completed':
        record['result']=json.loads((folder/'cold.json').read_text())
        record['reference_run']=run([a.checker,case,config,folder/'cold.json'],folder,'reference',120)
        if record['reference_run']['status']=='completed':record['validation']=json.loads((folder/'reference.json').read_text())
    summary['cases'].append(record)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(case,record['run'],record.get('validation',{}).get('status','not validated'),flush=True)
if any(case.get('validation',{}).get('status')!='pass' for case in summary['cases']):raise SystemExit(1)
