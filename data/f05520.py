"""Price frozen accepted R3 gates on exact leftover native source boxes.

Local q comparisons and spectral gates are fixed accepted inputs. This script
only performs rational/directed source transport and continuous box pricing.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];CORE=HERE.parent/'certificate'
sys.path.insert(0,str(CORE))
from certify_local_point import freeze as rankfreeze,validate
from price_frozen_profile import price as rankprice
from general_R3_cover import freeze as fullfreeze
from general_R3_price import price as fullprice
from price_finer_support_cells import bind
from mpfr_interval import UP

B=Q(7171,10000);TARGET=Q(40,71);FALLBACK=Q(451,800)

def profilebank():
 bank=[];files=[]
 for label in ['68','70']:
  p,paths,regions,t=bind(label);assert B<=p['eps'];p=validate({**p,'eps':B});files+=paths
  bank.append(dict(name='rank'+label,p=p,r0=Q(1,50000),r6=Q(1,5000000),g=B,kind='rank'))
 for label in ['68','70']:
  path=HERE/f'general_profile_u{label}_high_exact.json';p=fullfreeze(path);p=fullfreeze(path) # original exact six gates
  assert B<=p['eps'];p={**p,'eps':B};g=B-p['E00']+p['E02']**2/p['E22'];files.append(path)
  for prefix in ['body','exterior']:
   f=HERE/f'general_u{label}_high_{prefix}.json.gz';files.append(f)
   with gzip.open(f,'rt')as h:d=json.load(h)
   assert d['status']=='PASS_GENERAL_PHYSICAL_R3_COMPACT_COVER'and d['pending']==0
  bank.append(dict(name='full'+label,p=p,r0=Q(1,1000),r6=Q(1,10**6),g=g,kind='full'))
 for name,path,r0,r6 in [
  ('critical755',CORE/'frozen_R3_profile_exact.json',Q(1,50000),Q(1,1000000)),
  ('diffuse7478',HERE.parent/'diffuse_certificate/profile_interval_u7478_k0_01.json',Q(1,50000),Q(1,5000000)),
  ('diffuse7505',HERE.parent/'diffuse_certificate/profile_interval_b7171_newworst.json',Q(1,50000),Q(1,5000000))]:
  p=rankfreeze(path);assert B<=p['eps'];p=validate({**p,'eps':B});files.append(path)
  bank.append(dict(name=name,p=p,r0=r0,r6=r6,g=B,kind='rank'))
 return bank,files

def boxprice(profile,box):
 ends=[[Q(x)for x in box[k]]for k in ['u_interval','r_interval','kappa_interval']]+[[Q(2,5),Q(1)]]
 reserve=max(Q(0),-profile['g'])*(ends[1][1]-ends[1][0])**2/4
 values=[]
 for u,r,k,z in product(*ends):
  p=profile['p']
  v=rankprice(p,u,r,k,z,profile['r6'],profile['r0'])if profile['kind']=='rank'else fullprice(p,u,r,k,z,profile['r0'],profile['r6'])
  values.append(Q(v.hi.text(UP,60)))
 return max(values)+reserve,reserve

def main():
 source=ROOT/'tonight/recovery/profile_cover/native_upper_b7171_target_union.json';d=json.loads(source.read_text());assert Q(d['operator_b'])==B
 bank,files=profilebank();rows=[]
 for index,box in enumerate(d['unclosed_boxes']):
  candidates=[]
  for p in bank:
   upper,reserve=boxprice(p,box);candidates.append(dict(profile=p['name'],upper=str(upper),curvature_reserve=str(reserve)))
  winner=min(candidates,key=lambda t:Q(t['upper']));upper=Q(winner['upper'])
  rows.append(dict(native_box_index=index,**box,winner=winner['profile'],upper=winner['upper'],target_pass=upper<TARGET,fallback_pass=upper<FALLBACK,prices=candidates))
 files += [source,Path(__file__),HERE/'general_R3_cover.py',HERE/'general_R3_price.py',*[(CORE/x)for x in ['certify_local_point.py','price_frozen_profile.py','mpfr_interval.py']]]
 out=dict(status='PASS_DIRECTED_FROZEN_GATE_NATIVE_REMAINDER_PRICING',scope='R3 price acceptance only on listed continuous source boxes; input native complementary gates and full exact domain union require separate audit. Failed prices do not prove scalar norm counterexamples.',operator_b=str(B),source_native_union=str(source.relative_to(ROOT)),target=str(TARGET),fallback=str(FALLBACK),boxes=rows,target_accepted=sum(r['target_pass']for r in rows),fallback_accepted=sum(r['fallback_pass']for r in rows),sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
 (HERE/'native_upper_remainder_R3_all_accepted_prices.json').write_text(json.dumps(out,indent=2));print(out['status'],len(rows),out['target_accepted'],out['fallback_accepted'])
 from collections import defaultdict
 summary=defaultdict(lambda:[0,0,0,0])
 for row in rows:
  key=','.join(row['u_interval']);summary[key][0]+=1;summary[key][1]+=row['target_pass'];summary[key][2]+=row['fallback_pass'];summary[key][3]=max(summary[key][3],float(Q(row['upper'])))
 for k,v in summary.items():print(k,v)
if __name__=='__main__':main()
