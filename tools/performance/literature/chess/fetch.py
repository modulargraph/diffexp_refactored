#!/usr/bin/env python3
"""Fetch only the pinned public DP inputs; never clone the full repository."""
import argparse,hashlib,json,pathlib,urllib.request
SHA='4a701fc1332f29f6237d14427336e60615b966e3'
FILES=['Chess.wl','README.md','data/DP/Atilde.txt','data/DP/dLettersLine.txt','data/DP/boundary_pr1_weight6.txt','data/DP/reference_pr2_weight6.txt']
def fetch(root):
 root.mkdir(parents=True,exist_ok=True);hashes={}
 for name in FILES:
  target=root/name;target.parent.mkdir(parents=True,exist_ok=True)
  urllib.request.urlretrieve(f'https://raw.githubusercontent.com/Alice-Shimada/CHESS/{SHA}/{name}',target)
  hashes[name]=hashlib.sha256(target.read_bytes()).hexdigest()
 (root/'source-manifest.json').write_text(json.dumps({'repository':'https://github.com/Alice-Shimada/CHESS','commit':SHA,'sha256':hashes},indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('work_dir',type=pathlib.Path);args=parser.parse_args();fetch(args.work_dir.resolve())
