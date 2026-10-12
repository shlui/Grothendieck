"""Bind completed all-q gates and direct continuous outer source-cell prices."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;CORE=HERE.parent/'certificate'
sys.path.insert(0,str(CORE))
from certify_local_point import freeze
from certify_local_bernstein_region import tail
from price_frozen_profile import price
from mpfr_interval import UP

def load(path):
    with gzip.open(path,'rt')as f:return json.load(f)

def main():
    rows=[];files=[]
    for label,us,rs,ks,cap in [
        ('68',['.682','.69'],['.37','.43'],['0','.53'],Q(563,1000)),
        ('70',['.682','.69'],['.2','.398'],['0','1'],Q(40,71))]:
        profile=HERE/f'profile_u.{label}.json';p=freeze(profile)
        bodypath=HERE/f'body_u{label}_tight.json.gz';expath=HERE/f'exterior_u{label}_tight.json.gz'
        body=load(bodypath);ex=load(expath)
        for region,a,b in [(body,'0','8'),(ex,'8','64')]:
            assert region['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'and region['pending']==0
            assert region['qa']==a and region['qb']==b and Q(region['reserve0'])==Q(1,20000) and Q(region['reserve'])==Q(1,10**6)
            expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in p.items()};expected['major']=list(expected['major']);expected['major'][0]=str(Q(expected['major'][0])+Q(1,20000));expected['major'][3]=str(Q(expected['major'][3])+Q(1,10**6));assert expected==region['profile']
        inflated={k:list(v)if isinstance(v,list)else v for k,v in p.items()};inflated['major'][0]+=Q(1,20000);inflated['major'][3]+=Q(1,10**6)
        exacttail=tail(inflated,Q(64));assert exacttail['status']=='PASS'
        vertices=[]
        for u,r,k,z in product(*[[Q(t)for t in us],[Q(t)for t in rs],[Q(t)for t in ks],[Q(2,5),Q(1)]]):
            value=price(p,u,r,k,z,Q(1,10**6),Q(1,20000));assert value.below(cap)
            vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=value.hi.text(UP,50)))
        maximum=max(Q(v['upper'])for v in vertices)
        rows.append(dict(profile=label,status='PASS_ALL_Q_SOURCE_CELL',u=us,r=rs,kappa=ks,shape=['2/5','1'],rounded_cap=str(cap),maximum_upper=str(maximum),margin_to_40_71=str(Q(40,71)-maximum),body_leaves=body['accepted'],exterior_leaves=ex['accepted'],tail=exacttail,vertices=vertices))
        files.extend([profile,bodypath,expath])
    files.extend(CORE/name for name in ['certify_local_bernstein_region.py','certify_local_point.py','rational_phi_bounds.py','bernstein_integer_bounds.py','gaussian_max_certificate.py','mpfr_interval.py','price_frozen_profile.py'])
    out=dict(status='PASS_LOWER_AMPLITUDE_SOURCE_CELLS',scope='All Gaussian dimensions and physical signs with parameters in these source cells; not unrestricted amplitude strips or a global norm cap.',reserve0='1/20000',reserve6='1/1000000',cells=rows,sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
    (HERE/'lower_source_cells_tight.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['cells','sha256']},indent=2))
    for row in rows:print(row['profile'],row['rounded_cap'],row['maximum_upper'])

if __name__=='__main__':main()
