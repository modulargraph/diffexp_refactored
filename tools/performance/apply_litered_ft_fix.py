#!/usr/bin/env python3
"""Copy a LiteRed Setup tree and remove an invalid prefilter for merged FT quadratics.
The genuine Feynman-parametric/IBP scaleless tests remain unchanged.
Never modifies the installed source; refuses an existing destination.
"""
import argparse, hashlib, json, shutil
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('destination',type=Path);a=ap.parse_args()
source=a.source.resolve();destination=a.destination.resolve()
old='st=(Plus@@#<nloops)&/@sectors;'
sym_old='nm/:jSymmetries[Sequence@@ujs]=({m2rule[ujs,#]}&/@us)'
sym_new='nm/:jSymmetries[Sequence@@ujs]=If[Total[Rest[ujs]]<Length[LMs[nm]],{},({m2rule[ujs,#]}&/@us)]'
new='st=(Plus@@#==0)&/@sectors;'
f=source/'RNL/LiteRed183.m';raw=f.read_bytes();text=raw.decode()
if text.count(old)!=1 or text.count(sym_old)!=1:raise SystemExit('Expected exactly one known LiteRed 1.83 prefilter; inspect this version manually.')
if destination.exists():raise SystemExit('Destination must not exist.')
shutil.copytree(source,destination)
p=destination/'RNL/LiteRed183.m';p.write_text(text.replace(old,new).replace(sym_old,sym_new))
manifest={'source_sha256':hashlib.sha256(raw).hexdigest(),'patched_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'removed':old,'replacement':new,'optional_symmetries_removed':sym_old,'optional_symmetries_replacement':sym_new,'purpose':'Skip optional internal symmetries for low-denominator merged sectors: LiteRed otherwise emits free continuous-transformation parameters that FIRE cannot parse. Full IBP identities remain. Only the empty sector is unconditionally zero. Mixed FT quadratics must reach the unchanged actual scaleless test.'}
(destination/'ft-compatibility-patch.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest))
