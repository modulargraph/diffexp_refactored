import argparse,gzip,hashlib,json,os,pathlib,subprocess,time
p=argparse.ArgumentParser();p.add_argument('--request',type=pathlib.Path,required=True);p.add_argument('--executable',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--repeats',type=int,default=3);p.add_argument('--control',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
q=json.loads(a.request.read_text());records=[]
for repeat in range(a.repeats):
 for goal in [20,40]:
  for mode in (['auto','taylor'] if a.control and repeat==0 else ['auto']):
   q.update(working_bits=200,accuracy_goal=goal,recurrence=mode)
   payload=json.dumps(q);started=time.monotonic();child=subprocess.Popen([str(a.executable),'transport','-'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
   try:out,err=child.communicate(payload,timeout=600)
   except subprocess.TimeoutExpired:
    import signal
    os.killpg(child.pid,signal.SIGKILL);out,err=child.communicate();raise RuntimeError('600-second benchmark limit exceeded')
   elapsed=time.monotonic()-started;name=f'{mode}-{goal}-r{repeat+1}'
   with gzip.open(a.output/(name+'.json.gz'),'wt') as f:f.write(out)
   (a.output/(name+'.stderr')).write_text(err)
   if child.returncode:raise RuntimeError(err[-2000:])
   r=json.loads(out);record={'name':name,'goal':goal,'repeat':repeat+1,'mode':mode,'working_bits':200,'taylor_order':q['taylor_order'],'wall_seconds':elapsed,'timings':r['timings'],'spectral':r['spectral'],'recurrence':r['recurrence']};records.append(record);print(json.dumps(record),flush=True)
   report={'request_sha256':hashlib.sha256(a.request.read_bytes()).hexdigest(),'executable_sha256':hashlib.sha256(a.executable.read_bytes()).hexdigest(),'records':records,'scope':'Native total includes exact compilation and adaptive refinement; wall adds process startup, JSON input/output. Input conversion, reference validation and compilation excluded. Sequential runs; no concurrent compiler.'};(a.output/'summary.json').write_text(json.dumps(report,indent=2))
