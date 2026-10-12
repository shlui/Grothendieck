"""Directed continuous source boxes for two independently accepted q gates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];CORE=HERE.parent/'certificate';GENERAL=HERE.parent/'lower_amp_certificate'
sys.path.insert(0,str(CORE));sys.path.insert(0,str(GENERAL))
from certify_local_point import freeze as rankfreeze
from general_R3_cover import freeze as fullfreeze
from price_frozen_profile import price as rankprice
from general_R3_price import price as fullprice
from mpfr_interval import UP

def main():
 rows=[]
 cases=[dict(name='highk_global_worst_box',kind='full',profile='general_profile_highk_u7575_r259_k2_4.json',u=['151/200','19/25'],r=['129/500','5199/20000'],kappa=['11/10','4'],reserve0='1/200',reserve6='1/100000',cap='279/500',trees=['general_highk_body.json.gz','general_highk_exterior16_loose.json.gz'],audit='tonight/theory/targeted_highk_whole_q_audit.json',tail='general_highk_tail16_loose.json'),dict(name='lowr_center_box',kind='rank',profile='profile_interval_lowr_u76_r16_k0_2.json',u=['37/50','19/25'],r=['1/10','21/100'],kappa=['0','2'],reserve0='1/1000',reserve6='1/1000000',cap='559/1000',trees=['body_interval_lowr_u76_r16_k0_2.json.gz','exterior_interval_lowr_u76_r16_k0_2.json.gz'],audit='tonight/audit/targeted_lowr_R3_recheck/ACCEPTANCE.json',tail='tail_interval_lowr_u76_r16_k0_2.json')]
 for case in cases:
  path=HERE/case['profile'];p=(fullfreeze if case['kind']=='full'else rankfreeze)(path);r0,r6=Q(case['reserve0']),Q(case['reserve6']);g=p['eps']-p['E00']+p['E02']**2/p['E22']if case['kind']=='full'else p['eps'];assert g>=0
  files=[path,HERE/case['tail'],ROOT/case['audit'],Path(__file__)];regions=[]
  for tree in case['trees']:
   f=HERE/tree;files.append(f);d=json.load(gzip.open(f,'rt'));assert d['pending']==0 and d['status'].startswith('PASS_')
   local0=Q(d.get('reserve0'));local6=Q(d.get('reserve6',d.get('reserve')));assert local0<=r0 and local6<=r6
   regions.append(dict(q=[d['qa'],d['qb']],leaves=d['accepted'],recorded_reserve0=str(local0),recorded_reserve6=str(local6),global_majorant_domination='Nonnegative coefficientwise increment to q^0 and q^6.'))
  vertices=[]
  for u,r,k,z in product(*[[Q(x)for x in case[t]]for t in ['u','r','kappa']], [Q(2,5),Q(1)]):
   v=fullprice(p,u,r,k,z,r0,r6)if case['kind']=='full'else rankprice(p,u,r,k,z,r6,r0);assert v.below(Q(case['cap']))
   vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=v.hi.text(UP,60)))
  maximum=max(Q(v['upper'])for v in vertices);assert maximum<Q(40,71)
  rows.append(dict(**case,operator_b=str(p['eps']),shape=['2/5','1'],r_quadratic=str(g),maximum_directed_upper=str(maximum),margin_to_40_71=str(Q(40,71)-maximum),regions=regions,vertices=vertices,sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files}))
 out=dict(status='PASS_DIRECTED_CONTINUOUS_TARGETED_R3_SOURCE_BOXES',scope='Only these source boxes. Each fixed comparison is independently accepted on the complete q line; global all-source norm additionally requires domain union with native and other R3 gates.',boxes=rows)
 (HERE/'targeted_source_boxes.json').write_text(json.dumps(out,indent=2));print(out['status']);[print(x['name'],x['maximum_directed_upper'],x['cap'])for x in rows]
if __name__=='__main__':main()
