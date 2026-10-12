"""Accepted full-E local profile transported over continuous source bands."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import gzip,hashlib,json
from general_R3_cover import freeze
from general_R3_price import price,tail
from mpfr_interval import UP

HERE=Path(__file__).resolve().parent

def main():
    profile=HERE/'general_profile_u68_high_exact.json';p=freeze(profile)
    inflated={k:list(v)if isinstance(v,list)else v for k,v in p.items()};inflated['major'][0]+=Q(1,1000);inflated['major'][3]+=Q(1,10**6)
    expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in inflated.items()};regions=[];files=[profile,Path(__file__),HERE/'general_R3_cover.py',HERE/'general_R3_price.py']
    for a,b,prefix in [('0','8','body'),('8','64','exterior')]:
        f=HERE/f'general_u68_high_{prefix}.json.gz'
        with gzip.open(f,'rt')as fh:d=json.load(fh)
        assert d['status']=='PASS_GENERAL_PHYSICAL_R3_COMPACT_COVER'and d['pending']==0
        assert d['qa']==a and d['qb']==b and d['profile']==expected
        regions.append(dict(q=[a,b],events=d['nodes'],accepted=d['accepted'],excluded=d['excluded']));files.append(f)
    T=tail(p,Q(64),Q(1,1000),Q(1,10**6));assert T['status']=='PASS'
    g=p['eps']-p['E00']+p['E02']**2/p['E22'];slope=Q(9,5);polys=[]
    cases=[('target_lower68',['.682','.69'],['-.005','.02'],['0','1'],Q(539,1000))]
    for name,us,ts,ks,cap in cases:
        ue,te,ke=([Q(t)for t in xs]for xs in [us,ts,ks]);vertices=[]
        # Product corner interpolation of the sole negative quadratic g*r^2
        # contributes (-g)(slope^2*width_u^2+width_t^2)/4. No cross covariance.
        reserve=max(Q(0),-g)*(slope*slope*(ue[1]-ue[0])**2+(te[1]-te[0])**2)/4
        for u,t,k,z in product(ue,te,ke,[Q(2,5),Q(1)]):
            r=Q(851,2000)-slope*(u-Q(17,25))+t
            val=price(p,u,r,k,z,Q(1,1000),Q(1,10**6));upper=Q(val.hi.text(UP,50))
            assert upper+reserve<cap
            vertices.append(dict(u=str(u),t=str(t),r=str(r),kappa=str(k),shape=str(z),upper=str(upper)))
        largest=max(Q(v['upper'])for v in vertices)+reserve
        polys.append(dict(name=name,u=us,t=ts,kappa=ks,shape=['2/5','1'],parameterization='r=851/2000-(9/5)*(u-17/25)+t',negative_r_quadratic=str(g),curvature_reserve=str(reserve),rounded_cap=str(cap),maximum_directed_upper=str(largest),margin=str(cap-largest),vertices=vertices))
    out=dict(status='PASS_FULL_E_ALL_Q_AND_CONTINUOUS_SOURCE_BANDS',scope='These explicit source bands only; no unrestricted amplitude or all-kappa cover and no global norm cap.',profile=profile.name,reserve0='1/1000',reserve6='1/1000000',regions=regions,tail=T,polygons=polys,sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
    (HERE/'general68_source_polygons.json').write_text(json.dumps(out,indent=2));print(out['status'])
    for case in polys:print(case['name'],case['rounded_cap'],float(Q(case['maximum_directed_upper'])))

if __name__=='__main__':main()
