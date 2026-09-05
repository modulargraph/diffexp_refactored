#!/usr/bin/env python3
"""Bounded accepted/rejected chart guard probes; binary uses shadow headers."""
from pathlib import Path
import json,subprocess,time
root=Path(__file__).resolve().parents[2];records=[];start=time.monotonic()
for trial in range(3):
 for kind in range(3):
  for method in ['baseline','guarded_dot']:
   r=subprocess.run(['/tmp/diffexp-performance-guarded',method,str(kind),'3'],cwd=root,text=True,capture_output=True,timeout=60,check=True)
   row=json.loads(r.stdout);row.update(trial=trial+1,input_kind=kind,guard_probe=r.stderr.strip());records.append(row);print(json.dumps(row),flush=True)
p=root/'docs/validation/performance-guarded-dot.json';data=json.loads(p.read_text());data.update(chart_probe_seconds=time.monotonic()-start,chart_probes=records);p.write_text(json.dumps(data,indent=2)+'\n')
