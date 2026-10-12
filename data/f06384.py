"""Exact smaller-b transport of accepted profile70 over native remainder.

This is one restricted source box, not a global norm certificate.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from price_finer_support_cells import bind,CORE
from certify_local_point import validate
from price_frozen_profile import price
from mpfr_interval import UP

def main():
 p,paths,regions,tail=bind('70');b=Q(7171,10000);baseline=p['eps'];assert baseline==Q(18,25) and b<=baseline
 pb=validate({**p,'eps':b});cap=Q(40,71);vertices=[]
 for u,r,k,z in product([Q(137,200),Q(1381,2000)],[Q(47,125),Q(39,100)],[Q(0),Q(1)],[Q(2,5),Q(1)]):
  v=price(pb,u,r,k,z,Q(1,5000000),Q(1,50000));assert v.below(cap)
  vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),upper=v.hi.text(UP,60)))
 maximum=max(Q(v['upper'])for v in vertices);rounded=Q(5633,10000);assert maximum<rounded
 files=[Path(__file__),*paths,*[CORE/name for name in ['certify_local_bernstein_region.py','certify_local_point.py','rational_phi_bounds.py','bernstein_integer_bounds.py','gaussian_max_certificate.py','mpfr_interval.py','price_frozen_profile.py']]]
 out=dict(status='PASS_ALL_Q_B7171_NATIVE_REMAINDER_BOX',scope='Restricted source box only; global closure additionally requires exact accepted native complementary masks.',operator_b=str(b),baseline_b=str(baseline),u=['137/200','1381/2000'],r=['47/125','39/100'],kappa=['0','1'],shape=['2/5','1'],reserve0='1/50000',reserve6='1/5000000',regions=regions,tail=tail,rounded_exact_cap=str(rounded),maximum_upper=str(maximum),margin_to_target=str(cap-maximum),free_b_transport='(b-baseline_b)r^2+6 beta_raw^2(1/b-1/baseline_b); all six PSD gates monotone.',vertices=vertices,sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files})
 (HERE/'b7171_native_remainder_extended_box.json').write_text(json.dumps(out,indent=2));print(out['status'],out['maximum_upper'])
if __name__=='__main__':main()
