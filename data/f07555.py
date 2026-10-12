"""Exact coefficient transfer of the accepted complete outer amplitudes.

The physical partitions and scalar/local primitives are unchanged. All native
outer profiles have normalized longitudinal square beta^2<=1/16; the sole R3
replacement has raw longitudinal square 6*beta^2<=1/16. The exact transport
changes a pairing price by (b-b0)*r^2+beta_normalized^2*(1/b-1/b0).
Dropping its nonpositive first term gives one explicit outer cap for every
0<b<=b0. This does not certify a middle amplitude interval by itself.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[2]
B0=F(7171,10000)
PACKAGES=[
 'tonight/regional_certificates/'+x for x in
 ('u014_042','u030_0485','u040_054','u050_0638','u060_0662','u0662_0682','u077_c','u079_0797885')]
PACKAGES+=['tonight/native_frames/'+x for x in ('repriced0685_0700','repriced0700_0720')]
R3='tonight/operator/projection/lower_amp_certificate/profile_u.70.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    proofs={}
    for rel in ('tonight/theory/same_b_native_R3_glue_audit.json',
                'tonight/theory/global_native_R3_envelope_algebra_audit.json'):
        p=ROOT/rel;d=read(p);assert d['status'].startswith('PASS_') and F(d['operator_b'])==B0
        proofs[rel]=sha(p)
    records=[]
    for rel in PACKAGES:
        p=ROOT/rel/'data/local_summary.json';d=read(p)
        assert d['status']=='PASS'
        source=ROOT/rel/'data/source_profile.json'
        assert sha(source)==d['source_sha256']
        beta=F(d['profile']['beta']);assert beta*beta<=F(1,16)
        assert F(d['profile']['epsilon'])>=B0
        records.append({'profile':str(p.relative_to(ROOT)),'sha256':sha(p),
                        'source_sha256':sha(source),'normalization':'normalized_native',
                        'beta_square':str(beta*beta)})
    p=ROOT/R3;d=read(p)['profile'];beta=F(d['beta'])
    assert F(d['eps'])>=B0 and 6*beta*beta<=F(1,16)
    records.append({'profile':R3,'sha256':sha(p),'normalization':'raw_R3',
                    'beta_square':str(6*beta*beta)})
    b=F(713,1000);old_cap=F(56338,100000)
    shift=(1/b-1/B0)/16;cap=old_cap+shift
    assert cap==F(144154052987,255646150000)
    out={'status':'PASS_EXACT_B713_OUTER_TRANSPORT_FROM_ACCEPTED_B7171_BRANCHES',
         'operator_b':str(b),'baseline_b':str(B0),
         'scope':'All amplitudes [0,.6905] and [.775,sqrt(2/pi)]; not the middle interval.',
         'general_formula':'56338/100000+(1/b-1/(7171/10000))/16 for every0<b<=7171/10000',
         'baseline_uniform_cap':str(old_cap),'uniform_beta_square_cap':'1/16',
         'transport_shift':str(shift),'outer_cap':str(cap),
         'accepted_coverage_proofs':proofs,'profile_budget_checks':records,
         'analytic_prerequisite':'Exact native/R3 smaller-b transport; elementary low-u bound decreases with b.',
         'claim':'Both unchanged complete physical source covers transfer to one common coefficient b=.713.'}
    Path(__file__).with_name('smaller_b_outer_transport_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k!='profile_budget_checks'},indent=2))
if __name__=='__main__':main()
