"""Assemble independently replayable all-q profile coverage and outer cell.

This does not declare a global operator norm: only the specified outer cell is
certified.  The global uncovered source region is explicitly retained.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,gzip,json,hashlib
from itertools import product
from certify_local_point import freeze,ROOT
from certify_local_region import tail
from price_frozen_profile import price
from mpfr_interval import UP

def read(path):
    if str(path).endswith('.gz'):
        with gzip.open(path,'rt') as f:return json.load(f)
    return json.loads(Path(path).read_text())

def exactprofile(p):return {k:[str(t) for t in v] if isinstance(v,list) else str(v) for k,v in p.items()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--profile',default=str(Path(__file__).with_name('frozen_R3_profile_exact.json')));ap.add_argument('--body',default=str(ROOT/'results'/'directed_bernstein_region_0_8.json.gz'));ap.add_argument('--exterior',default=str(ROOT/'results'/'directed_region_8_64.json.gz'));ap.add_argument('--output',default=str(ROOT/'results'/'R3_complete_local_and_outer_cell.json'));args=ap.parse_args()
    p=freeze(args.profile);body=read(args.body);ex=read(args.exterior)
    assert body['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM' and body['pending']==0 and body['qa']=='0' and body['qb']=='8'
    assert ex['status']=='PASS_COMPACT_LOCAL_SUPREMUM' and ex['pending']==0 and ex['qa']=='8' and ex['qb']=='64'
    reserve6=Q(body['reserve']);reserve0=Q(body['reserve0'])
    pbody={k:list(v) if isinstance(v,list) else v for k,v in p.items()};pbody['major'][0]+=reserve0;pbody['major'][3]+=reserve6
    assert exactprofile(pbody)==body['profile']
    # The old exterior majorant is stronger: raising only its constant/sextic
    # coefficient preserves its directed gate pointwise on all real q.
    for k,v in ex['profile'].items():
        if k!='major':assert v==exactprofile(pbody)[k]
    em=[Q(t) for t in ex['profile']['major']]
    assert em[1:3]==pbody['major'][1:3] and em[0]<=pbody['major'][0] and em[3]<=pbody['major'][3]
    t=tail(pbody,Q(64));assert t['status']=='PASS'
    cap=Q(28167,50000) # .56334, comfortably below 40/71.
    rows=[];largest=Q(-1)
    for u,r,k,z in product([Q(7549,10000),Q(7551,10000)],[Q(1287,5000),Q(2577,10000)],[Q(51,100),Q(11,20)],[Q(2,5),Q(1)]):
        v=price(p,u,r,k,z,reserve6,reserve0);assert v.below(cap)
        U=Q(v.hi.text(UP,50));largest=max(largest,U);rows.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=str(U)))
    assert cap<Q(40,71)
    files=[Path(args.profile),Path(args.body),Path(args.exterior),Path(__file__),Path(__file__).with_name('certify_local_point.py'),Path(__file__).with_name('certify_local_bernstein_region.py'),Path(__file__).with_name('rational_phi_bounds.py'),Path(__file__).with_name('bernstein_integer_bounds.py'),Path(__file__).with_name('gaussian_max_certificate.py'),Path(__file__).with_name('mpfr_interval.py'),Path(__file__).with_name('price_frozen_profile.py')]
    out=dict(status='PASS_UNIFORM_R3_OUTER_CELL',scope='All Gaussian dimensions, all physical signs whose polarization/source parameters fall in the specified cell, all shared tensor shapes and every real scalar q. This is not a global operator norm or a new Grothendieck lower bound.',operator='P1-(18/25)P3',cell=dict(u=['7549/10000','7551/10000'],r=['1287/5000','2577/10000'],kappa=['51/100','11/20'],shape=['2/5','1']),rounded_exact_upper=str(cap),maximum_directed_upper=str(largest),margin_to_40_71=str(Q(40,71)-cap),reserve0=str(reserve0),reserve6=str(reserve6),all_q_regions=[dict(q=['0','8'],method='physical moment-center quadratic completion and exact rational Bernstein curvature',nodes=body['nodes'],leaves=body['accepted']),dict(q=['8','64'],method='full dual-box product corner interpolation',nodes=ex['nodes'],leaves=ex['accepted']),dict(q=['64','infinity'],method='exact shifted sextic polynomial',tail=t)],remaining_gate='Cover every admissible outer u,r,kappa source parameter; combine with exterior/amplitude certificates. The present cell alone gives no global norm cap.',profile=exactprofile(p),rows=rows,sha256={str(f.resolve()):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    Path(args.output).write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['profile','rows','sha256','all_q_regions']},indent=2))

if __name__=='__main__':main()
