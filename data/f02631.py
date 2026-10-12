#!/usr/bin/env python3
"""Exact physical exclusion of the κ[2,4], u>=.755, r>=.258 bottleneck.

Uses the already proved global direct-contact family with B=0. No Gaussian
law assumption: prices normalized third-chaos Q using the fresh633 theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

import certify_direct_contact_majorant as C

HERE=Path(__file__).resolve().parent
A=F(147,50)
B=F(0)
W=F(11,2)
ULO=F(151,200)
RLO=F(129,500)
KLO,KHI=F(2),F(4)


def main():
    C.A0=A
    C.B0=B
    contact=C.row(W)
    q=contact['polynomial_coefficients_q0_q2_q4_q6']
    corners=[]
    for k in [KLO,KHI]:
        for shape in [0,1]:
            price=C.price(q,k,shape)
            bound=C.upper(price,10**20)
            corners.append({'kappa':str(k),'shape':shape,
                            'price':C.record(price),'price_upper':str(bound),
                            'source_upper_at_ulo':str(bound-A*ULO)})
    price_upper=max(F(c['price_upper']) for c in corners)
    source_upper=price_upper-A*ULO
    clean_source=F(257485167,10**9)
    assert source_upper<clean_source<RLO
    result={'status':'PASS_EXACT_GLOBAL_LINEAR_CONTACT_PHYSICAL_SOURCE_EXCLUSION',
            'bits':C.ctx.prec,'parameters':{'A':str(A),'B':str(B),'w':str(W)},
            'contact_row':contact,'moment_price_corners':corners,
            'uniform_price_upper':str(price_upper),
            'u_lower':str(ULO),'kappa_interval':[str(KLO),str(KHI)],
            'shape_interval':['0','1'],'source_r_uniform_upper':str(source_upper),
            'clean_source_r_upper':str(clean_source),
            'excluded_r_lower':str(RLO),
            'exclusion_margin':str(RLO-source_upper),
            'proof':'For the fixed globally valid B=0 row, A*u+r<=E J_A(Q)'
                    '<=P(kappa,shape). Its q6>=0 price is convex in kappa and '
                    'affine in shape, hence the four corner upper prices '
                    'bound all kappa[2,4],shape[0,1]. Since A>0, '
                    'u>=151/200 gives r<=uniform_price_upper-A*151/200'
                    '<257485167/10^9<129/500.',
            'scope':'A physical source-domain exclusion, independent of operator b. '
                    'Not a complete norm cap or Grothendieck endpoint.',
            'analytic_family_source_sha256':hashlib.sha256(
                (HERE/'certify_direct_contact_majorant.py').read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'linear_contact_exclusion.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','parameters',
                        'source_r_uniform_upper','clean_source_r_upper',
                        'exclusion_margin']},indent=2))


if __name__=='__main__':
    main()
