"""All-q R3 coverage of the physical r>=.2 domain at small cumulants.

Only u in [.682,.69], kappa in [0,1] is covered. Native cells must separately
cover r<=.2; larger cumulants and the rest of the amplitudes remain separate.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;CORE=HERE.parent/'certificate'
sys.path.insert(0,str(CORE))
from certify_local_point import freeze
from certify_local_bernstein_region import tail
from price_frozen_profile import price
from mpfr_interval import I,UP

def bind(label):
    path=HERE/f'profile_u.{label}.json';p=freeze(path)
    inflated={k:list(v)if isinstance(v,list)else v for k,v in p.items()};inflated['major'][0]+=Q(1,50000);inflated['major'][3]+=Q(1,5000000)
    expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in inflated.items()}
    paths=[];regions=[]
    for a,b,prefix in [('0','8','body'),('8','64','exterior')]:
        f=HERE/f'{prefix}_u{label}_finer.json.gz'
        with gzip.open(f,'rt')as fh:d=json.load(fh)
        assert d['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'and d['pending']==0
        assert d['qa']==a and d['qb']==b and Q(d['reserve0'])==Q(1,50000)and Q(d['reserve'])==Q(1,5000000)
        assert d['profile']==expected;paths.append(f);regions.append(dict(q=[a,b],leaves=d['accepted'],events=d['nodes']))
    t=tail(inflated,Q(64));assert t['status']=='PASS'
    return p,[path,*paths],regions,t

def main():
    l1path=HERE.parents[1]/'cubic_L1_refinement/kappa_shape_caps.json'
    l1=json.loads(l1path.read_text());assert l1['status']=='PASS_CERTIFIED_KAPPA_SHAPE_CUBIC_L1_CAPS'
    row=next(v for v in l1['table']if v['name']=='coarse0_1');assert row['kappa_interval']==['0','1']
    rho=Q(830644,1000000);assert Q(row['rho'])<=rho
    a=Q(9,5);B=Q(1649921,1000000);c=(I(2)/I.pi()).sqrt();rr=I(rho);aa=I(a);x=1/rr
    Ja=aa*c*(-x*x/(2*aa*aa)).exp()+x*(x/(I(2).sqrt()*aa)).erf()
    bridge=(1-rr*rr)*aa*c+rr*rr*Ja;assert bridge.below(B)
    # Exact dimension-free mass/source bridge au+r<=E max(a|Z|,|Q|)
    # and the complete monotonicity quadratic majorant for J_a give this gate.
    cells=[];files=[l1path,Path(__file__)]
    for label,cap in [('70',Q(5631,10000)),('68',Q(28166,50000))]:
        p,paths,regions,t=bind(label);files.extend(paths)
        vertices=[]
        for u in [Q(341,500),Q(69,100)]:
            rs=[Q(1,5),Q(199,500)]if label=='70'else[Q(199,500),B-a*u]
            assert rs[0]<=rs[1]
            for r,k,z in product(rs,[Q(0),Q(1)],[Q(2,5),Q(1)]):
                val=price(p,u,r,k,z,Q(1,5000000),Q(1,50000));assert val.below(cap)
                vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=val.hi.text(UP,50)))
        largest=max(Q(v['upper'])for v in vertices)
        cells.append(dict(profile=label,u=['341/500','69/100'],r='[1/5,199/500]'if label=='70'else'[199/500,1649921/1000000-(9/5)u]',kappa=['0','1'],shape=['2/5','1'],rounded_exact_cap=str(cap),maximum_upper=str(largest),margin_to_40_71=str(Q(40,71)-largest),regions=regions,tail=t,vertices=vertices))
    files.extend(CORE/name for name in ['certify_local_bernstein_region.py','certify_local_point.py','rational_phi_bounds.py','bernstein_integer_bounds.py','gaussian_max_certificate.py','mpfr_interval.py','price_frozen_profile.py'])
    out=dict(status='PASS_ALL_Q_R3_SMALL_CUMULANT_HIGH_SOURCE_DOMAIN',
             scope='All physical signs and Gaussian dimensions with u in [.682,.69], kappa in [0,1], and r>=.2. The valid source bridge bounds every remaining r above; this does not cover r<.2 or kappa>1, or prove a global operator norm.',
             reserve0='1/50000',reserve6='1/5000000',rho=str(rho),support_gate=dict(a=str(a),upper=str(B),directed_price_upper=bridge.hi.text(UP,50),proof_input='Complete J_a quadratic Hermite interpolation bridge, using E Q^2=1 and E|Q|<=rho.'),cells=cells,
             sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
    (HERE/'lower_source_domain_finer.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['cells','sha256']},indent=2))
    for v in cells:print(v['profile'],v['rounded_exact_cap'],v['maximum_upper'])

if __name__=='__main__':main()
