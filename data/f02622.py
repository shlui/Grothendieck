#!/usr/bin/env python3
"""Exact all-halfline cubic-head potential majorant plus directed moment prices.

Checks a dimension-free analytic derivative inequality, not sampled points.
Source geometry and Taylor integration are proved in CUBIC_HEAD_SOURCE_GATE.md.
No Stein scalar fit or numerical Gaussian quadrature enters acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx
HERE=Path(__file__).resolve().parent
ctx.prec=512
mctx=fmpq_mpoly_ctx.get(('A','B','w'));A,B,w=mctx.gens();P=mctx.constant(1)
D1=A+B*w;D3=A+3*B*w;r=0;s=1;derivatives=[]
for j in range(1,4):
 P=(fmpq(1,2)*P-P.derivative(2))*D1*D3+r*B*P*D3+3*s*B*P*D1
 r+=2;s+=2
 assert all(c>=0 for c in P.coeffs()) and any(c>0 for c in P.coeffs())
 derivatives.append({'order':j,'sign':(-1)**j,'denominator_powers':[r,s],'positive_numerator':str(P),'coefficient_count':len(P.coeffs()),'coefficient_minimum':str(min(P.coeffs()))})
A0=F(97,40);B0=F(1,6);a0=A0+3*B0
q=[a0-B0,1/(2*A0),-(A0+6*B0)/(24*A0**4),(A0*A0+20*A0*B0+120*B0*B0)/(240*A0**7)]
def aa(x):
 x=F(x);return arb(x.numerator)/x.denominator
c=(2/arb.pi()).sqrt();K=194/arb(10).sqrt();d=arb(6).sqrt()*aa(B0)
k=F(53,100);u=F(151,200);v=F(-1,4);prices=[]
for zeta in [0,1]:
 V=k*(F(7,45)+F(zeta,15))
 sixth=15+15*aa(k)+(K-aa(F(633*zeta,15)))*aa(k)*aa(k).sqrt()
 price=c*(aa(q[0]+q[1]+q[2]*(3+k))+aa(q[3])*sixth)
 source=price-aa(a0*u)-d*aa(v)
 high=F(str(price.upper().fmpq()));sr=F(str(source.upper().fmpq()))
 assert high< (F(235915722,10**8) if zeta==0 else F(235725855,10**8))
 assert sr<(F(252844293,10**9) if zeta==0 else F(250945617,10**9))
 prices.append({'zeta':zeta,'V':str(V),'A_upper':str(high),'A_display':float(price),'source_r_upper_at_u755_vminus025':str(sr),'source_r_display':float(source),'source_r_simple_upper':str(F(252844293,10**9) if zeta==0 else F(250945617,10**9))})
out={'status':'PASS_EXACT_CUBIC_HEAD_ALL_HALF_LINE_MAJORANT','scope':'Analytic physicalsource gate, notglobalGrothendiecknorm endpoint. Momentpriceusesincludedfresh633theorem.','derivative_check':derivatives,'parameters':{'A':str(A0),'B':str(B0),'a':str(a0),'d':'sqrt(6)/6'},'polynomial_majorant_coefficients_without_c':list(map(str,q)),'moment_price_formula':'c*[q0+q2+(3+k)*q4+q6*(15+15*k+(194/sqrt(10))*k^(3/2)-633*sqrt(k)*(V-7*k/45))]','shape_convention':'V=k*(7/45+zeta/15),0<=zeta<=1','example':{'kappa':str(k),'u':str(u),'v':str(v),'endpoint_prices':prices},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'cubic_head_polynomial_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'q_without_c':out['polynomial_majorant_coefficients_without_c'],'examples':[{'zeta':p['zeta'],'price':p['A_display'],'r':p['source_r_display']}for p in prices]}),flush=True)
