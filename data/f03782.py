#!/usr/bin/env python3
"""Exact global cubic-Hermite-contact upper bound for J(sqrt(s)).

The derivative sign is proved by exact multivariate polynomial coefficients;
the interpolation remainder has sign s^2(s-t)^2*J_fourth<=0 everywhere.
All sextic x-coefficients are rounded upward to exact rationals, and all
source moment prices are directed. No scalar samples or optimizer premise.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from flint import arb,ctx,fmpq_mpoly_ctx

HERE=Path(__file__).resolve().parent
ctx.prec=512
A0=F(26412631,10**7)
B0=F(4084047,10**7)
CONTACTS=[F(3,5),F(7,10),F(13,10),F(2),F(3),F(189,50),F(17,4)]


def aa(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def upper(v,den=10**20):
    r=F(str(v.upper().fmpq()))
    return F(-((-r.numerator*den)//r.denominator),den)


def record(v):
    return {'lower':str(F(str(v.lower().fmpq()))),
            'upper':str(F(str(v.upper().fmpq()))),'display':float(v)}


def derivative_proof():
    ring=fmpq_mpoly_ctx.get(('A','B','w','tau'))
    A,B,w,tau=ring.gens()
    D1=A+B*w
    D3=A+3*B*w
    P=ring.constant(1)
    # K_tau(s)=exp(-tau*w(s))/(A+B*w(s)), s=w(A+Bw)^2.
    r,s=1,0
    rows=[]
    for j in range(1,4):
        P=(tau*P-P.derivative(2))*D1*D3+r*B*P*D3+3*s*B*P*D1
        r+=2
        s+=2
        assert all(c>=0 for c in P.coeffs())
        assert any(c>0 for c in P.coeffs())
        rows.append({'order':j,'sign':(-1)**j,'positive_numerator':str(P),
                     'denominator_powers':[r,s],
                     'coefficient_count':len(P.coeffs()),
                     'minimum_coefficient':str(min(P.coeffs()))})
    return rows


def row(w):
    t=w*(A0+B0*w)**2
    x=aa(w).sqrt()*aa(A0+B0*w)
    nu=(2/arb.pi()).sqrt()
    erf=(aa(w)/2).sqrt().erf()
    J0=nu*aa(A0+2*B0)
    Jp0=nu/aa(2*A0)
    Jt=x*erf+nu*(-aa(w)/2).exp()*aa(A0+B0*(w+2))
    Jpt=erf/(2*x)
    delta=Jt-J0-Jp0*aa(t)
    derivative_delta=Jpt-Jp0
    q6=(aa(t)*derivative_delta-2*delta)/aa(t**3)
    q4=(3*delta-aa(t)*derivative_delta)/aa(t**2)
    assert q6>0
    coefficients=[J0,Jp0,q4,q6]
    rationals=[upper(v) for v in coefficients]
    assert all(aa(q)>=v for q,v in zip(rationals,coefficients))
    assert rationals[3]>0
    return {'w':str(w),'contact_s':str(t),
            'contact_J':record(Jt),'contact_Jprime_s':record(Jpt),
            'polynomial_coefficients_q0_q2_q4_q6':list(map(str,rationals)),
            'unrounded_coefficients':list(map(record,coefficients))}


def price(q,k,shape):
    k=F(k)
    q0,q2,q4,q6=map(F,q)
    sixth=15+15*aa(k)+(194/arb(10).sqrt()-aa(F(633)*F(shape)/15))*aa(k)*aa(k).sqrt()
    return aa(q0+q2+(3+k)*q4)+aa(q6)*sixth


def main():
    proof=derivative_proof()
    rows=[row(w) for w in CONTACTS]
    for r in rows:
        q=r['polynomial_coefficients_q0_q2_q4_q6']
        r['prices']=[{'kappa':str(k),'shape':shape,
                      'upper':str(upper(price(q,k,shape),10**16)),
                      'price':record(price(q,k,shape))}
                    for k in [F(0),F(1,10),F(53,100),F(1),F(2),F(3),F(4)]
                    for shape in [0,1]]
        r['uniform_price_on_kappa_2_4_all_shapes']=str(max(
            F(p['upper']) for p in r['prices'] if F(p['kappa']) in [F(2),F(4)]))
    result={'status':'PASS_EXACT_ALL_HALFLINE_DIRECT_CUBIC_HEAD_CONTACT_MAJORANT',
            'bits':ctx.prec,'parameters':{'A':str(A0),'B':str(B0),
                  'a':str(A0+3*B0),'d':'sqrt(6)*'+str(B0)},
            'derivative_positivity':proof,'rows':rows,
            'global_majorant':'J(x)<=q0+q2*x^2+q4*x^4+q6*x^6 for every real x.',
            'source_gate':'a*u+d*v+r <= price(kappa,shape), '
                          'v=E[H3(Z)h], r=||P3^W E_Z k||2.',
            'moment_price_formula':'q0+q2+(3+kappa)*q4+'
               '[15+15*kappa+(194/sqrt(10)-633*shape/15)*kappa^(3/2)]*q6',
            'continuous_price_proof':'q6>=0; 194/sqrt(10)-633*shape/15>0 '
                 'on shape[0,1], hence convex in kappa and affine in shape. '
                 'A fixed row is priced at four corners of each rectangle.',
            'moment_dependency':'reduced_search/energy/633/PROOF.md',
            'scope':'Exact additional physical source constraint; no global norm endpoint.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert 194/arb(10).sqrt()-arb(633)/15>0
    (HERE/'direct_contact_majorants.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'parameters':result['parameters'],
                     'rows':[{'w':r['w'],'contact_s':r['contact_s'],
                         'coefficients':r['polynomial_coefficients_q0_q2_q4_q6'],
                         'kappa2_4_uniform':r['uniform_price_on_kappa_2_4_all_shapes']}
                         for r in rows]},indent=2))


if __name__=='__main__':
    main()
