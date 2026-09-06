#!/usr/bin/env python3
"""Run a product-recurrence experiment with a sixty-second process cap."""
import argparse,json,subprocess
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('binary',type=Path);ap.add_argument('--cached',action='store_true');ap.add_argument('--output',type=Path);a=ap.parse_args()
cmd=[str(a.binary.resolve()),'--run']
if a.cached:
    f=Path(__file__).resolve().parent/'fixtures';cmd += [str(f/'closure.json'),str(f/'rows.json')]
try:r=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
except subprocess.TimeoutExpired:raise SystemExit('Sixty-second cap exceeded; no result accepted.')
print(r.stdout,end='');print(r.stderr,end='')
if r.returncode:raise SystemExit(r.returncode)
data=json.loads(r.stdout)
if data.get('validation')!='passed':raise SystemExit('Validation did not pass.')
if a.output:a.output.write_text(json.dumps(data,indent=2)+'\n')
