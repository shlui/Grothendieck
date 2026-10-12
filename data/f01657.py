"""Coefficient-independent ceiling for pure cubic reciprocal-norm sequences.
Two exact rational reflected sign witnesses bracket a zero cubic pairing.
Continuity and monotonicity supply an exact intervening witness, not a norm cap.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,time
import certify_cubic720_norm_witness as C
from mpfr_interval import I
HERE=Path(__file__).resolve().parent
start=time.time();C.B=F(40812698,10**8);C.EPS=F(0)
records=[]
for a in [F(264042705,10**8),F(264042707,10**8)]:
 C.A=a;records.append(C.compute(check_fixed_witness=False))
def interval(r,k):return I(r[k]['lower'],r[k]['upper'])
lo,hi=records
assert interval(lo,'P3_pairing').below(0)
assert interval(hi,'P3_pairing').above(0)
constant=interval(lo,'P1_pairing');ceiling=1/constant
assert constant.above('0.56323352359')
assert ceiling.below('1.775462501')
result={'status':'PASS_DIRECTED_COEFFICIENT_INDEPENDENT_CUBIC_FAMILY_CEILING','parameters':{'a_low':lo['parameters']['a'],'a_high':hi['parameters']['a'],'c':lo['parameters']['c']},'P3_pairing_low':lo['P3_pairing'],'P3_pairing_high':hi['P3_pairing'],'universal_scalar_norm_lower':constant.record(),'universal_reciprocal_norm_ceiling':ceiling.record(),'ceiling_above_1_775':(ceiling-I(F(71,40))).record(),'precision_argument':'At fixed positive c, the cubic pairing is continuous in a and changes sign between these rational endpoints. There is an exact a_star with cubic pairing zero. The first Hermite moment is increasing in a, because its derivative is sqrt(2/pi)*E[Z^2 exp(-q(Z)^2/2)]>0. Therefore this exact threshold has first-chaos pair at least the lower endpoint first-chaos pair, uniformly for every real cubic coefficient b.','scope':'For every real b, sup_d ||P1-bP3||_(infinity-to-1) is at least this constant. Its reciprocal only limits this particular pure-cubic reciprocal-norm construction; it is neither a universal Grothendieck upper bound nor an achieved lower bound. Target1.775 remains possible.','endpoint_records':records,'seconds':time.time()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'cubic_family_ceiling.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','P3_pairing_low','P3_pairing_high','universal_scalar_norm_lower','universal_reciprocal_norm_ceiling','seconds']},indent=2))
