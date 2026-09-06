#!/usr/bin/env python3
"""Audit exact FT homogeneous stages without computing endpoint or integral values."""
import argparse,json,pathlib,subprocess,tempfile,os
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('caches',nargs='+',help='FAMILY=CACHE_DIRECTORY or FAMILY=CLOSURE_FILE')
p.add_argument('--output',required=True)
p.add_argument('--kernel',default='/Applications/Wolfram.app/Contents/MacOS/WolframKernel')
a=p.parse_args(); records=[]
for spec in a.caches:
 name,location=spec.split('=',1);root=pathlib.Path(location)
 for path in ([root] if root.is_file() else root.glob('*/*.json')):
  try:d=json.loads(path.read_text())
  except (ValueError,OSError):continue
  if d.get('payload',{}).get('schema')!='DiffExp3.ExactLevelClosure/v1':continue
  records.append(dict(family=name,path=str(path.resolve()),physical_count=d['identity']['family']['physical_count'],matrix=d['payload']['matrix'],basis=d['payload']['ordered_basis']))
if not records:raise SystemExit('No exact closure records found; scalar leaves do not have these matrices.')
with tempfile.TemporaryDirectory(prefix='ft-structure-') as td:
 source=pathlib.Path(td)/'closures.json';source.write_text(json.dumps(records))
 env=dict(os.environ,DIFFEXP_FT_AUDIT_INPUT=str(source),DIFFEXP_FT_AUDIT_OUTPUT=str(pathlib.Path(a.output).resolve()))
 subprocess.run([a.kernel,'-script',str(pathlib.Path(__file__).with_suffix('.wls')),str(source),str(pathlib.Path(a.output).resolve())],check=True,timeout=120,env=env)
 if not pathlib.Path(a.output).is_file():raise SystemExit("Wolfram diagnostic produced no result")
