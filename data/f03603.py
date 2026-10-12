#!/usr/bin/env python3
"""Certify saved inverse caps via independent FORWARD monotone source gate."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from mpfr_interval import I,BITS
ROOT=Path(__file__).resolve().parents[3]
a=argparse.ArgumentParser();a.add_argument('ledger',type=Path);a.add_argument('--output',type=Path,required=True);arg=a.parse_args()
src=arg.ledger.resolve();d=json.loads(src.read_text());assert d['physical_u_clip']
nu=(I(2)/I.pi()).sqrt();rho=F(847,1000);rows=[]
for b in d['physical_u_bounds']:
 r,u=F(b['r_lower']),F(b['upper']);assert 0<=r<rho
 if nu.below(u):method='physical u<=nu endpoint'
 else:
  assert r>0
  z=(I(u)/nu-1+rho*rho)/(rho*rho);assert z.above(0)and z.below(1)
  forward=I(rho)*(-z.log()).sqrt().erf();assert forward.below(r)
  method='R_rho(u_cap)<=r_lo via forward log/sqrt/erf'
 rows.append({'r_lower':str(r),'u_upper':str(u),'method':method})
out={'status':'PASS_INDEPENDENT_FORWARD_CERTIFICATION_OF_PHYSICAL_INVERSE_CAPS','input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'precision_bits':BITS,'bounds_checked':len(rows),'source_rho':str(rho),'inverse_function_used':False,'producer_code_imported':False,'rows':rows,'analytic_review':'R_rho is decreasing. Every physical r>=r_lo obeys r_lo<=R_rho(u), hence the independent forward comparison R_rho(u_cap)<=r_lo implies u<=u_cap. Endpoint u_cap>=nu follows directly. These caps are bounding rectangles; continuous source coverage is checked separately.'}
arg.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='rows'},indent=2))
