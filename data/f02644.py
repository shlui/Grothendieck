#!/usr/bin/env python3
"""Export source-price endpoints without discarding the shared Malliavin shape."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ctx.prec=512
A=lambda x:arb(F(x).numerator)/F(x).denominator
src=HERE/'J_source_k053_fresh_directed.json'
if not src.exists():src=HERE/'J_source_fresh_directed.json'
data=json.loads(src.read_text());out={'status':'PASS_FRESH_J_SOURCE_SHAPE_PRICE_EXPORT','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_path':str(src),'formula':'A(k,V)=c*sqrt(1+a^2)+q0+q2+(3+k)*q4+[15+15*k+(194/sqrt(10))*k^(3/2)+(633*7/45)*k^(3/2)]*q6+[eta/2-633*sqrt(k)*q6]*V','shape_convention':'V=k*(7/45+zeta/15), 0<=zeta<=1','rows':[]}
c=(2/arb.pi()).sqrt();K=194/arb(10).sqrt()
for row in data['rows']:
 assert row['status']=='PASS_DIRECTED_GENERIC_J_WHOLE_LINE'
 k=F(row['kappa_discovery']);a=F(row['a']);u=F(row['u']);q=list(map(F,row['majorant']));eta=F(row['eta'])
 baseline=c*(1+A(a*a)).sqrt()+A(q[0]+q[1]+q[2]*(3+k))+A(q[3])*(15+15*A(k)+K*A(k)*A(k).sqrt())+A(eta*k*F(7,90))
 slope=A(eta*k/F(30))-A(F(633,15)*q[3])*A(k)*A(k).sqrt()
 prices=[]
 for z in [0,1]:
  p=baseline+z*slope;prices.append({'zeta':z,'V':str(k*(F(7,45)+F(z,15))),'A_upper':str(F(str(p.upper().fmpq()))),'A_display':float(p),'r_upper_at_recorded_u':str(F(str(p.upper().fmpq()))-a*u),'r_display_at_recorded_u':float(p)-float(a*u)})
 out['rows'].append({'a':str(a),'kappa':str(k),'u':str(u),'eta':row['eta'],'majorant':row['majorant'],'phi_coefficients':row['phi_coefficients'],'degree':row['degree'],'scale':row['scale'],'zeta_slope_display':float(slope),'endpoint_prices':prices})
 print(json.dumps({key:out['rows'][-1][key]for key in ['a','kappa','u','zeta_slope_display','endpoint_prices']}),flush=True)
(HERE/'Ja_shape_prices.json').write_text(json.dumps(out,indent=2)+'\n')
