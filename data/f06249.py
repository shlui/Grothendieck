#!/usr/bin/env python3
"""Portable complete replay of new R3 gates and independently expanded prices."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import argparse, hashlib, json, os, subprocess, sys, time
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'replay_outputs_expanded'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def check_manifest():
 d=json.loads((ROOT/'MANIFEST.json').read_text())
 for rel,h in d['sha256'].items():assert digest(ROOT/rel)==h,rel
 return d

def run(tag,args):
 env=dict(os.environ);env['PYTHONINTMAXSTRDIGITS']='0';env['PYTHONPATH']=os.pathsep.join(str(ROOT/p)for p in ['tonight/audit','tonight/theory','tonight/operator/projection/certificate'])
 with (OUT/(tag+'.log')).open('w')as log:
  p=subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError('Replay failed: '+tag)
 return tag

def price_sources(specs):
 sys.path.insert(0,str(ROOT/'tonight/audit/global_b7171'))
 import replay_R3_corners as R
 from replay_middlebank_R3_corners import fresh_tail
 result=[]
 for spec in specs:
  source=ROOT/'tonight/operator/projection/targeted_certificate'/(spec['ledger']+'_source_box.json');d=R.load(source)
  assert d['status']=='PASS_DIRECTED_ROBUST_RESTRICTED_R3_SOURCE_BOXES'
  for box in d['boxes']:
   for rel,h in box['sha256'].items():assert digest(ROOT/rel)==h,rel
   profile=source.parent/box['profile'];p=R.exact(R.load(profile)['profile'])
   row=dict(name=box['name'],profile=profile,p=p,kind='rank',reserve0=F(box['reserve0']),reserve6=F(box['reserve6']))
   coefficient=R.coefficients(row,R.B);R.spectral(row,R.B)
   assert coefficient['r_quadratic']==R.B>0 and F(box['tail_start'])==spec['R']
   tail=fresh_tail(row,F(spec['R']));assert all(F(t)>=0 for t in tail['shifted_coefficients'])
   points=set(product(*[tuple(map(F,box[k]))for k in ('u','r','kappa','shape')]))
   largest=F(0)
   for saved in box['vertices']:
    point=tuple(F(saved[k])for k in ('u','r','kappa','shape'));assert point in points;points.remove(point)
    u,r,k,z=point;actual=R.price(row,coefficient,k,z,u,r);assert actual.below(F(saved['upper']))
    largest=max(largest,F(actual.hi.text(2,100)))
   assert not points and largest<F(box['cap'])<F(451,800)
   result.append(dict(name=spec['label'],corners=16,cap=box['cap'],independent_upper=str(largest)))
 return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--all',action='store_true');ap.add_argument('--workers',type=int,default=4);a=ap.parse_args();OUT.mkdir(exist_ok=True);start=time.time();check_manifest()
 specs=json.loads((ROOT/'EXPANDED_SPECS.json').read_text());jobs=[]
 for spec in specs:
  for part in ('body','exterior','tail'):
   args=[ROOT/'tonight/audit'/spec['checker'],'--profile',spec['label'],'--branch',spec['branch'],'--part',part,'--tail-R',spec['R'],'--data-dir',ROOT/'tonight/operator/projection/targeted_certificate','--output-dir',OUT]
   jobs.append((spec['label']+'_'+part,args))
 with ThreadPoolExecutor(a.workers)as pool:
  for tag in pool.map(lambda job:run(*job),jobs):print('PASS',tag,flush=True)
 prices=price_sources(specs)
 if a.all:run('earlier_targeted_complete',[ROOT/'verify_targeted.py','--workers',a.workers])
 check_manifest();report=dict(status='PASS_PORTABLE_EXPANDED_RESTRICTED_R3_REPLAY',fresh_partitions=len(specs),exact_tails=len(specs),fresh_source_corners=16*len(specs),seconds=time.time()-start,prices=prices,scope='Accepted whole-q gates and only the listed restricted source boxes. Complete global amplitude/source coverage is separate.')
 (OUT/'result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
