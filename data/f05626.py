"""Completed all-q profile and exact continuous high-source masks near CLT.

The mask covers the specified upper r domains, not every r at those amplitudes.
The smaller-b spectral transport is checked exactly; no local rerun is needed.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;CORE=HERE.parent/'certificate'
sys.path.insert(0,str(CORE))
from certify_local_point import freeze,validate
from certify_local_bernstein_region import tail
from price_frozen_profile import price
from mpfr_interval import I,UP

def main():
    profile=HERE/'profile_interval_b7171_newworst.json';p=freeze(profile)
    inflated={k:list(v)if isinstance(v,list)else v for k,v in p.items()};inflated['major'][0]+=Q(1,50000);inflated['major'][3]+=Q(1,5000000)
    expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in inflated.items()};files=[profile,Path(__file__)]
    regions=[]
    for a,b,prefix in [('0','8','body'),('8','64','exterior')]:
        f=HERE/f'{prefix}_interval_b7171_newworst.json.gz'
        with gzip.open(f,'rt')as fh:d=json.load(fh)
        assert d['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'and d['pending']==0 and d['qa']==a and d['qb']==b and d['profile']==expected
        files.append(f);regions.append(dict(q=[a,b],leaves=d['accepted'],events=d['nodes']))
    T=tail(inflated,Q(64));assert T['status']=='PASS';(HERE/'tail_interval_b7171_newworst.json').write_text(json.dumps(T,indent=2))
    l1path=HERE.parents[1]/'cubic_L1_refinement/kappa_shape_caps.json';l1=json.loads(l1path.read_text());row=next(x for x in l1['table']if x['name']=='coarse0_1');rho=Q(830644,10**6);assert Q(row['rho'])<=rho;files.append(l1path)
    c=(I(2)/I.pi()).sqrt();a=I(3);rr=I(rho);x=1/rr;J=a*c*(-x*x/(2*a*a)).exp()+x*(x/(I(2).sqrt()*a)).erf();support=(1-rr*rr)*a*c+rr*rr*J
    assert (support-3*I(Q(743,1000))).below(Q(3,10))
    masks=[];witness_u=Q('0.750515528144');witness_r=Q('0.270685360072')
    for b,cap in [(Q(7171,10000),Q(40,71)-Q(1,1000000))]:
        assert 0<b<=p['eps'];pb={**p,'eps':b};validate(pb)
        # Local values and all Stein/tail comparisons do not use epsilon.
        # Six PSD gates monotonically persist when epsilon decreases.
        constant=max(Q(price(pb,Q(0),Q(0),k,z,Q(1,5000000),Q(1,50000)).hi.text(UP,60))for k,z in product([Q(0),Q(1,10)],[Q(2,5),Q(1)]))
        vy=p['y']/p['p2']-p['p0']*p['r'];gamma=p['gamma']-p['p0']*vy
        assert gamma>b*Q(3,10)
        def F(u,r):return u*u-2*p['ell']*u+b*r*r-2*gamma*r+constant
        endpoints=[]
        for u in [Q(j,1000)for j in range(743,756)]:
            # Exact integer bisection proposes the lower r boundary, then the
            # saved rational point is verified against the directed constant.
            lo,hi=0,300000000000
            assert F(u,Q(hi,10**12))<cap
            while hi-lo>1:
                mid=(lo+hi)//2
                if F(u,Q(mid,10**12))<cap:hi=mid
                else:lo=mid
            r=Q(hi,10**12)+Q(1,10**12);assert F(u,r)<cap
            endpoints.append((u,r))
        cells=[]
        for (ua,ra),(ub,rb)in zip(endpoints[:-1],endpoints[1:]):
            vertices=[(ua,ra),(ub,rb),(ua,Q(3,10)),(ub,Q(3,10))]
            upper=max(F(u,r)for u,r in vertices);assert upper<cap
            cells.append(dict(u=[str(ua),str(ub)],lower_r_vertices=[str(ra),str(rb)],upper_r='3/10',cap=str(cap),maximum_upper=str(upper)))
        actual=max(Q(price(pb,witness_u,witness_r,k,z,Q(1,5000000),Q(1,50000)).hi.text(UP,50))for k,z in product([Q(0),Q(1,10)],[Q(2,5),Q(1)]))
        masks.append(dict(operator_b=str(b),cap=str(cap),kappa=['0','1/10'],shape=['2/5','1'],source_constant_upper=str(constant),effective_gamma=str(gamma),cells=cells,example_rational_point=dict(u=str(witness_u),r=str(witness_r),upper=str(actual)),free_b_price_change='(b-7171/10000)*r^2+6*beta_raw^2*(1/b-10000/7171)'))
    files.extend(CORE/name for name in ['certify_local_bernstein_region.py','certify_local_point.py','rational_phi_bounds.py','bernstein_integer_bounds.py','gaussian_max_certificate.py','mpfr_interval.py','price_frozen_profile.py'])
    out=dict(status='PASS_ALL_Q_CONTINUOUS_CRITICAL_HIGH_SOURCE_MASKS',scope='These explicit high-r source polygons only. Sources below their lower boundary, kappa>.1 and other amplitudes require complementary gates. No global norm or improved Grothendieck endpoint.',reserve0='1/50000',reserve6='1/5000000',regions=regions,tail=T,support_upper=support.hi.text(UP,50),support_r_bound='All physical r<3/10 for u>=743/1000 and kappa<=1, by accepted cubic L1+J_a bridge.',masks=masks,sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
    (HERE/'critical_source_upper_masks_b7171_newworst.json').write_text(json.dumps(out,indent=2));print(out['status'])
    for m in masks:print(m['operator_b'],m['cap'],m['example_rational_point'])

if __name__=='__main__':main()
