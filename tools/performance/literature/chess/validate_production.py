import argparse,gzip,json,pathlib,re,mpmath as mp
mp.mp.dps=115

def number(text):
 text=re.sub(r'`{1,2}[+-]?[0-9.]+','',text).replace('*^','e')
 terms=re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?(?:\*I)?',text)
 if ''.join(terms)!=text:raise ValueError('Unsupported numeric reference: '+text[:80])
 return sum((mp.mpc(0,mp.mpf(t[:-2])) if t.endswith('*I') else mp.mpc(mp.mpf(t))) for t in terms)

def main():
 p=argparse.ArgumentParser();p.add_argument('--reference',type=pathlib.Path,required=True);p.add_argument('--results',type=pathlib.Path,required=True);a=p.parse_args()
 raw=re.sub(r'\s+','',a.reference.read_text().replace('\\\n',''))
 rows=re.findall(r'\{([^{}]*)\}',raw);ref=[[number(x) for x in row.split(',')] for row in rows]
 if len(ref)!=267 or any(len(row)!=7 for row in ref):raise ValueError('Expected all267x7 reference coefficients')
 checks=[]
 for f in sorted(list(a.results.glob('auto-*.json.gz'))+list(a.results.glob('taylor-*.json.gz'))):
  result=json.loads(gzip.decompress(f.read_bytes()));vals=[[mp.mpc(c['real_midpoint'],c['imaginary_midpoint']) for c in row] for row in result['values']]
  if len(vals)!=len(ref) or any(len(v)!=len(r) for v,r in zip(vals,ref)):raise ValueError('Endpoint dimension mismatch')
  absolute=max(abs(v-r) for row,rr in zip(vals,ref) for v,r in zip(row,rr));normalized=max(abs(v-r)/max(1,abs(r)) for row,rr in zip(vals,ref) for v,r in zip(row,rr));goal=int(f.name.split('-')[1]);checks.append({'file':f.name,'coefficient_count':1869,'max_absolute_error':str(absolute),'max_normalized_error':str(normalized),'requested_digits':goal,'passed':bool(normalized<mp.mpf(10)**-goal)})
 (a.results/'validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
 if not checks or not all(c['passed'] for c in checks):raise SystemExit('Reference check failed')
if __name__=='__main__':main()
