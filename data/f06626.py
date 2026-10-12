"""Restricted source boxes, bound to independently replayed whole-q gates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib,gzip,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'certificate'))
from certify_local_point import freeze
from price_frozen_profile import price
from mpfr_interval import UP

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 cases=[dict(name='interval_lowerface_u7175_k053_08',u=['143/200','18/25'],r=['44/125','177/500'],kappa=['53/100','4/5'],reserve0='1/10000',reserve6='1/100000000',cap='281493/500000',branch='lowerface',audit='lowerface_u7175_R3_recheck')]
 rows=[]
 for c in cases:
  path=HERE/f"profile_{c['name']}.json";p=freeze(path);assert p['eps']==Q(7171,10000)
  assert p['G']-p['d0']-p['zeta']**2>=0 and p['major'][3]+Q(c['reserve6'])+p['zeta']**2>=0
  files=[path,Path(__file__)];replays=[];regions=[]
  for part in ('body','exterior','tail'):
   audit=ROOT/'tonight/audit'/c['audit']/f"independent_{c['name']}_{c['branch']}_{part}.json"
   d=json.loads(audit.read_text());assert d['status'].startswith('PASS_INDEPENDENT_');files.append(audit)
   if part!='tail':
    tree=HERE/f"{part}_{c['name']}.json.gz";files.append(tree)
    assert d['tree_sha256']==sha(tree) and d['profile_sha256']==sha(path)
    t=json.load(gzip.open(tree,'rt'));assert t['pending']==0
    assert Q(t['reserve0'])==Q(c['reserve0']) and Q(t['reserve'])==Q(c['reserve6'])
    regions.append(dict(q=[t['qa'],t['qb']],accepted=d['accepted']))
   else:
    tail=HERE/f"tail_{c['name']}.json";files.append(tail);t=json.loads(tail.read_text())
    assert d['R']==t['R']=='32' and d['difference_polynomial']==t['difference_polynomial'] and d['shifted_coefficients']==t['shifted_coefficients']
   replays.append(str(audit.relative_to(ROOT)))
  vertices=[]
  for u,r,k,z in product(*[[Q(t)for t in c[s]]for s in ('u','r','kappa')],(Q(2,5),Q(1))):
   v=price(p,u,r,k,z,Q(c['reserve6']),Q(c['reserve0']));assert v.below(Q(c['cap']))
   vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=v.hi.text(UP,70)))
  maximum=max(Q(t['upper'])for t in vertices);assert maximum<Q(c['cap'])<Q(451,800)
  rows.append(dict(**c,profile=path.name,kind='rank',operator_b=str(p['eps']),shape=['2/5','1'],r_quadratic=str(p['eps']),tail_start='32',regions=regions,independent_q_replays=replays,maximum_directed_upper=str(maximum),vertices=vertices,sha256={str(f.relative_to(ROOT)):sha(f)for f in files}))
 out=dict(status='PASS_DIRECTED_ROBUST_RESTRICTED_R3_SOURCE_BOXES',scope='Only these continuous source boxes. Complete amplitude/source union is separate.',convexity='Positive rank-one r curvature; u quadratic positive. Cumulant price separately convex in kappa and shared shape, including positive covariance term.',boxes=rows)
 (HERE/'lowerface_source_box.json').write_text(json.dumps(out,indent=2));print(out['status']);[print(t['name'],float(Q(t['maximum_directed_upper'])))for t in rows]
if __name__=='__main__':main()
