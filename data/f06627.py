"""Directed restricted source box for the independently accepted middle gate."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,sys,hashlib,gzip
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'certificate'))
from certify_local_point import freeze
from price_frozen_profile import price
from mpfr_interval import UP

profile=HERE/'profile_interval_middle_u7475_r267_k08_1.json'
p=freeze(profile);r0=Q(1,10000);r6=Q(1,1000000)
audit=ROOT/'tonight/audit/targeted_middle_R3_recheck/ACCEPTANCE.json'
assert json.loads(audit.read_text())['status'].startswith('PASS_')
files=[profile,audit,Path(__file__)]
regions=[]
for part in ('body','exterior'):
 f=HERE/f'{part}_interval_middle_u7475_r267_k08_1.json.gz';files.append(f)
 d=json.load(gzip.open(f,'rt'));assert d['pending']==0 and d['status'].startswith('PASS_')
 assert Q(d['reserve0'])==r0 and Q(d.get('reserve6',d.get('reserve')))==r6
 regions.append(dict(q=[d['qa'],d['qb']],accepted=d['accepted']))
tail=HERE/'tail_interval_middle_u7475_r267_k08_1.json';files.append(tail)
vertices=[]
for u,r,k,z in product((Q(149,200),Q(3,4)),(Q(133,500),Q(67,250)),(Q(4,5),Q(1)),(Q(2,5),Q(1))):
 v=price(p,u,r,k,z,r6,r0);assert v.below(Q(56302,100000))
 vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=v.hi.text(UP,70)))
maximum=max(Q(v['upper'])for v in vertices)
out=dict(status='PASS_DIRECTED_CONTINUOUS_MIDDLE_R3_SOURCE_BOX',scope='Only this source box; complete source coverage is separate.',operator_b=str(p['eps']),profile=str(profile.relative_to(ROOT)),u=['149/200','3/4'],r=['133/500','67/250'],kappa=['4/5','1'],shape=['2/5','1'],reserve0=str(r0),reserve6=str(r6),cap='28151/50000',maximum_directed_upper=str(maximum),regions=regions,tail_start='16',vertices=vertices,sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
(HERE/'middle_source_box.json').write_text(json.dumps(out,indent=2))
print(out['status'],float(maximum))
