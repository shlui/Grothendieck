"""Independent expanded-price replay of smaller-b native remainder boxes."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import gzip,hashlib,json,sys
ROOT=Path(__file__).resolve().parent
BASE=ROOT/'tonight/operator/projection/lower_amp_certificate'
sys.path.insert(0,str(ROOT/'tonight/audit/lower_amp_R3'));sys.path.insert(0,str(ROOT/'tonight/audit'))
from independent_source_price import exact_profile,coefficients,price,bind_local_metadata

def main():
 p=exact_profile(BASE/'profile_u.70.json');assert p['eps']==F(18,25)
 metadata=bind_local_metadata('70',p,'finer',F(1,50000),F(1,5000000));rows=[]
 for name in ['b7171_native_remainder_box.json','b7171_native_remainder_extended_box.json']:
  d=json.loads((BASE/name).read_text());b=F(d['operator_b']);assert 0<b<=p['eps']
  for path,h in d['sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
  pb={**p,'eps':b};assert pb['p2']**2>=b*(1+pb['p0']**2)/2
  for a,c,e in [(pb['R11'],pb['R13'],pb['R33']),(pb['R11']+1,pb['R13'],pb['R33']-b/6),(pb['R11']-b,pb['R13'],pb['R33'])]:assert a>=0 and e>=0 and a*e>=c*c
  c=coefficients(pb,F(1,50000),F(1,5000000));vertices=set(product(*[[F(t)for t in d[k]]for k in ['u','r','kappa','shape']]))
  maximum=F(-10**6)
  for row in d['vertices']:
   key=tuple(F(row[k])for k in ['u','r','kappa','shape']);assert key in vertices;vertices.remove(key)
   v=price(pb,c,*key);assert v.below(F(row['upper'])) and v.below(F(d['rounded_exact_cap']));maximum=max(maximum,F(v.hi.text(2)))
  assert not vertices and maximum<=F(d['maximum_upper'])
  rows.append({'ledger':name,'maximum_independent_upper':str(maximum),'cap':d['rounded_exact_cap']})
 print(json.dumps({'status':'PASS_INDEPENDENT_B7171_RESTRICTED_BOX_PRICES','boxes':rows},indent=2))
if __name__=='__main__':main()
