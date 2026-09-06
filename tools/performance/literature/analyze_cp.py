import pathlib,json,re,mpmath as mp
mp.mp.dps=115
root=pathlib.Path(__file__).resolve().parent

def value(s):
 s=re.sub(r'`{1,2}[+-]?[0-9.]+','',s).replace('*^','e').replace(' ','')
 terms=re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?(?:\*I)?',s)
 if ''.join(terms)!=s:raise ValueError(s)
 return sum((mp.mpc(0,mp.mpf(t[:-2])) if t.endswith('*I') else mp.mpc(mp.mpf(t))) for t in terms)
report=[]
for f in ['RL1','RL2']:
 native=json.loads((root/(f+'-native.json')).read_text());original=json.loads((root/(f+'-original-result.json')).read_text())
 references=[[value(s) for s in row] for row in original['values']]
 if len(native['values'])!=len(references) or any(len(a)!=len(b) for a,b in zip(native['values'],references)):raise ValueError('endpoint shape mismatch')
 vals=[[mp.mpc(c['real_midpoint'],c['imaginary_midpoint']) for c in row] for row in native['values']]
 worst=max(abs(a-b)/max(1,abs(b)) for row,ref in zip(vals,references) for a,b in zip(row,ref))
 case={'family':f,'masters':len(vals),'epsilon_order':6,'original_transport_seconds':original['seconds'],'native_timings':native['timings'],'recurrence':native['recurrence'],'maximum_scaled_cross_method_difference':str(worst),'passes_40_digits':bool(worst<mp.mpf('1e-40'))}
 path=root/(f+'-wrapper-result.json')
 if path.exists():
  wrapper=json.loads(path.read_text());wv=[[value(s) for s in row] for row in wrapper['values']];err=max(abs(a-b)/max(1,abs(b)) for row,ref in zip(wv,references) for a,b in zip(row,ref))
  case.update(wrapper_seconds=wrapper['seconds'],wrapper_difference=str(err),wrapper_passes_40_digits=bool(err<mp.mpf('1e-40')))
 report.append(case)
print(json.dumps(report,indent=2));(root/'cp-summary.json').write_text(json.dumps(report,indent=2))

if not all(c["passes_40_digits"] and c.get("wrapper_passes_40_digits",True) for c in report):raise SystemExit("40-digit comparison failed")
