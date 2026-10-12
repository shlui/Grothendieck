#!/usr/bin/env python3
"""Exact global off-origin Hermite-contact upper bound for a cubic-head source.

No sampled scalar comparison or Stein fit enters this proof. Exact polynomial
positivity proves F'''<=0; the quadratic Hermite interpolation error has the
sign of s(s-t)^2 F''' throughout s>=0. Directed Arb rounds all resulting
coefficients upward to rationals and prices included fresh633 moment bounds.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from flint import arb,ctx,fmpq,fmpq_mpoly_ctx

HERE=Path(__file__).resolve().parent
ctx.prec=512
A0=F(26412631,10**7)
B0=F(4084047,10**7)
CONTACTS=[F(1,4),F(5,8),F(1)]


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
    ring=fmpq_mpoly_ctx.get(('A','B','w'))
    A,B,w=ring.gens()
    D1=A+B*w
    D3=A+3*B*w
    P=ring.constant(1)
    r,s=0,1
    rows=[]
    for j in range(1,4):
        P=(fmpq(1,2)*P-P.derivative(2))*D1*D3+r*B*P*D3+3*s*B*P*D1
        r+=2
        s+=2
        assert all(c>=0 for c in P.coeffs())
        assert any(c>0 for c in P.coeffs())
        rows.append({'order':j,'sign':(-1)**j,'numerator':str(P),
                     'denominator_powers':[r,s],
                     'coefficient_count':len(P.coeffs()),
                     'minimum_coefficient':str(min(P.coeffs()))})
    return rows


def row(w):
    D1=A0+B0*w
    D3=A0+3*B0*w
    t=w*D1*D1
    exp=(-aa(w)/2).exp()
    f=exp/aa(D3)
    fp=-exp*aa(D3/2+3*B0)/aa(D1*D3**3)
    f0=aa(1/A0)
    U2=(f0-f+aa(t)*fp)/aa(t*t)
    U1=fp-2*aa(t)*U2
    assert U2>0
    nu=(2/arb.pi()).sqrt()
    real_coeff=[nu*aa(A0+2*B0),nu*f0/2,nu*U1/12,nu*U2/30]
    rational_coeff=[upper(v) for v in real_coeff]
    assert all(aa(q)>=v for q,v in zip(rational_coeff,real_coeff))
    assert rational_coeff[3]>0
    return {'w':str(w),'contact_s':str(t),
            'F0':record(f0),'F_contact':record(f),'Fprime_contact':record(fp),
            'U1':record(U1),'U2':record(U2),
            'polynomial_coefficients_q0_q2_q4_q6':list(map(str,rational_coeff)),
            'unrounded_coefficients':list(map(record,real_coeff))}


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
        r['prices']=[{'kappa':str(k),'shape':shape,'upper':str(upper(price(q,k,shape),10**16)),
                      'price':record(price(q,k,shape))}
                     for k in [F(0),F(1,10),F(53,100),F(1),F(2),F(4)] for shape in [0,1]]
    result={'status':'PASS_EXACT_ALL_HALFLINE_CUBIC_HEAD_HERMITE_CONTACT_MAJORANT',
            'bits':ctx.prec,'parameters':{'A':str(A0),'B':str(B0),'a':str(A0+3*B0),
                                          'd':'sqrt(6)*'+str(B0)},
            'derivative_positivity':proof,'rows':rows,
            'global_majorant':'J(x)<=q0+q2*x^2+q4*x^4+q6*x^6 for every real x.',
            'source_gate':'a*u+d*v+r <= price(kappa,shape), v=E[H3(Z)h], '
                          'r=||P3^W E_Z k||2.',
            'moment_price_formula':'q0+q2+(3+kappa)*q4+'
                '[15+15*kappa+(194/sqrt(10)-633*shape/15)*kappa^(3/2)]*q6',
            'moment_dependency':'reduced_search/energy/633/PROOF.md',
            'scope':'Exact additional physical source constraint; no global norm endpoint.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'contact_majorants.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'parameters':result['parameters'],
                      'rows':[{'w':r['w'],'coefficients':r['polynomial_coefficients_q0_q2_q4_q6'],
                               'prices':[{'kappa':p['kappa'],'shape':p['shape'],'upper':p['upper']} for p in r['prices']]}
                              for r in rows]},indent=2))


if __name__=='__main__':
    main()
