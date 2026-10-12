#!/usr/bin/env python3
"""Independent read-only exact arithmetic for the universal upper corollary."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse,hashlib,json
R=Path(__file__).resolve().parents[2];H=R/'new_round_2026-10-02/audit'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def bind(p):return {'path':str(p.relative_to(R)),'sha256':sha(p)}
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--lower-acceptance',type=Path);a=ap.parse_args()
review_path=H/'ELEMENTARY_KRIVINE_UPPER_INDEPENDENT_AUDIT.json';review=json.loads(review_path.read_text())
assert review['status']=='PASS_INDEPENDENT_ANALYTIC_AND_EXACT_KRIVINE_UPPER_REVIEW'
for record in review['bindings'].values():assert sha(R/record['path'])==record['sha256']
x=F(141,160);s=F(1);term=F(1)
for k in range(1,13):term*=x/k;s+=term
term*=x/13;U=s+term/(1-x/14)
assert U<F(239,99) and F(140,99)**2<2
P=list(map(F,[4,0,-4,0,5,-4,1]));product=[F(0)]*9
for i,c in enumerate(P):product[i]+=c;product[i+2]+=c
product[0]-=4
assert product==list(map(F,[0,0,0,0,1,-4,6,-4,1]))
assert sum((c/(i+1)for i,c in enumerate(P)),F(0))==F(22,7)
upper=F(22,7)/(2*x);assert upper==F(1760,987)<F(357,200)
source_path=R/review['bindings']['ELEMENTARY_KRIVINE_UPPER_EXACT_CERTIFICATE.json']['path'];source=json.loads(source_path.read_text())
assert F(source['exponential_upper'])==U and F(source['universal_constant_strict_upper'])==upper
out={'status':'PASS_INDEPENDENT_READ_ONLY_ELEMENTARY_KRIVINE_UPPER_REPLAY','universal_strict_upper':str(upper),
 'strict_below_rounding_endpoint':'357/200','exponential_upper':str(U),'exact_rational_arithmetic':True,
 'three_significant_digit_lower_required':'71/40','three_significant_digit_rounding_proved':False,
 'bindings':{'accepted_analytic_review':bind(review_path)},'entrypoint_sha256':sha(Path(__file__)),
 'scope':'Unconditional universal upper bound. Rounding1.78 additionally needs an independently accepted complete constructive lower bound at least71/40.'}
if a.lower_acceptance:
 p=a.lower_acceptance.resolve();d=json.loads(p.read_text());assert d['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
 lower=F(d['constructive_G_R_lower_bound']);cap=F(d['norm_upper_bound']);assert lower==1/cap and lower>=F(71,40)
 out['bindings']['complete_constructive_lower']=bind(p);out['constructive_lower']=str(lower)
 out['three_significant_digit_rounding_proved']=True;out['rounding_result']='1.78'
 out['certified_rounding_interval']=['71/40','1760/987 (strict)']
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
