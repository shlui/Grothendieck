"""Replay directed vertex prices of the accepted profile's affine polygons.

The completed whole-q profile is reused by exact y(r) transport.  These are
source cells, and do not cover all admissible r at any amplitude.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse, hashlib, json
from certify_local_point import freeze, ROOT
from price_frozen_profile import price
from mpfr_interval import UP

POLYGONS = {
    'target': (['.7538', '.7562'], ['-.00025', '.00025'], ['.3', '.8'], Q(28168,50000)),
    'fallback': (['.75', '.76'], ['-.001', '.001'], ['.3', '.8'], Q(451,800)),
    'fallback_wide_kappa': (['.7505', '.7595'], ['-.001', '.001'], ['.15', '1.1'], Q(451,800)),
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--profile', default=str(Path(__file__).with_name('frozen_R3_profile_exact.json')))
    ap.add_argument('--output', default=str(ROOT/'results'/'moving_polygons_replay.json'))
    args=ap.parse_args(); p=freeze(args.profile)
    # Exact rank-one even source matrix, nonnegative scalar curvature costs,
    # and positive tensor envelope coefficients give vertex reduction.
    assert p['major'][3]+Q(1,10**6)+p['zeta']**2 >= 0
    assert p['G']-p['d0']-p['zeta']**2 >= 0
    records=[]
    for name,(us,ts,ks,cap) in POLYGONS.items():
        rows=[]
        for u,t,k,z in product(*[[Q(x) for x in us],[Q(x) for x in ts],[Q(x) for x in ks],[Q(2,5),Q(1)]]):
            r=Q(5151,20000)-3*(u-Q(151,200))+t
            val=price(p,u,r,k,z,Q(1,10**6),Q(1,50000))
            assert val.below(cap)
            rows.append(dict(u=str(u),t=str(t),r=str(r),kappa=str(k),shape=str(z),upper=val.hi.text(UP,50)))
        upper=max(Q(row['upper']) for row in rows)
        records.append(dict(name=name,u=us,t=ts,kappa=ks,shape=['2/5','1'],
                            parameterization='r=5151/20000-3*(u-151/200)+t',
                            rounded_exact_cap=str(cap),maximum_directed_upper=str(upper),
                            margin=str(cap-upper),rows=rows))
    out=dict(status='PASS_DIRECTED_AFFINE_POLYGON_PRICES',
             scope='Uses the separately completed whole-q local profile. Each source polygon only; no global operator norm or all-r amplitude cover.',
             reserve0='1/50000',reserve6='1/1000000',polygons=records,
             profile_sha256=hashlib.sha256(Path(args.profile).read_bytes()).hexdigest())
    Path(args.output).write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='polygons'},indent=2))
    for row in records: print(row['name'],row['rounded_exact_cap'],row['maximum_directed_upper'])

if __name__=='__main__':main()
