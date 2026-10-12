"""Directed arithmetic for the smaller-grid explicit cubic matrix sequence.

An all-dimension scalar norm cap is a separate proof obligation. No reported
arithmetic value alone certifies a Grothendieck lower bound.
"""
from argparse import ArgumentParser
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
from flint import arb,ctx

def qa(x):
    x=F(x);return arb(x.numerator)/x.denominator
def floor_decimal(x,places=15):
    q=F(str(x.lower().fmpq()));scale=10**places
    return F((q.numerator*scale)//q.denominator,scale)
def ceil_decimal(x,places=60):
    q=F(str(x.upper().fmpq()));scale=10**places
    return F(-((-q.numerator*scale)//q.denominator),scale)
def dec(q,places=15):
    q=F(q);scale=10**places;k=q.numerator*scale//q.denominator
    return f'{k//scale}.{k%scale:0{places}d}'

def row(n,power,b,B):
    assert n>=2
    d=n**power;h=F(1,d);eta=F(1,d*d)
    # Floating nomination only; the accepted integer is verified below.
    radius=max(1,math.ceil(math.sqrt(10*math.log(d))))
    logd=qa(d).log()
    while radius>1 and (qa((radius-1)**2)-10*logd).lower()>0:radius-=1
    while not (qa(radius*radius)-10*logd).lower()>0:radius+=1
    assert (qa((radius-1)**2)-10*logd).upper()<0
    phi=(-qa(radius*radius)/2).exp()/(2*arb.pi()).sqrt()
    prob=qa(F(2*d,radius))*phi
    assert prob.upper()<1
    alpha=qa(1-b/2)+qa(b)/2*qa(F(d+4,d+2)).sqrt()
    loss=qa(d*h*h*(F(1-b/2,d-1)+F(b,2*(d+2))+F(b,d-2)))+4*alpha*prob.sqrt()+4*qa(b)*prob
    # Classical chi-mean inequality a_d >= 1-1/(2d), as in the baseline.
    radial=qa((1-F(1,2*d))*(1-F(b,2*(d+2))))
    numerator=radial-loss
    ratio=(numerator-qa(eta))/(qa(B)+qa(eta))
    M=2*radius*d+2
    return dict(n=n,dimension=d,mesh=str(h),integer_cutoff=radius,
        cutoff_inequality='R^2 >= 10 log(d)',cutoff_difference=str(qa(radius*radius)-10*logd),
        intervals_per_coordinate=M,seed_rows=f'{M}^{d}/2',rationalization_error=str(eta),
        tail_probability_upper=str(ceil_decimal(prob)),
        quadratic_loss_upper=str(ceil_decimal(loss)),
        radial_lower=str(floor_decimal(radial,60)),numerator_lower=str(floor_decimal(numerator,60)),
        finite_ratio_lower=str(floor_decimal(ratio)),finite_ratio_decimal=dec(floor_decimal(ratio)))

def main():
    p=ArgumentParser(description=__doc__);p.add_argument('--b',required=True,type=F)
    p.add_argument('--cap',required=True,type=F);p.add_argument('--dimension-power',type=int,default=2)
    p.add_argument('--norm-proof',type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--n',nargs='+',type=int,default=[8,16,32,64,128,1024]);a=p.parse_args()
    assert 0<=a.b<1 and a.cap>0 and a.dimension_power in (2,4)
    ctx.prec=768
    proof=None
    if a.norm_proof:
        proof=dict(path=str(a.norm_proof.resolve()),sha256=hashlib.sha256(a.norm_proof.read_bytes()).hexdigest())
    data=dict(status='PASS_DIRECTED_FINITE_FORMULA_ARITHMETIC',operator_b=str(a.b),required_global_norm_cap=str(a.cap),
        limit_of_certified_lower_formula=str(1/a.cap),dimension_power=a.dimension_power,
        proof_binding=proof,precision_bits=ctx.prec,
        scope='Explicit smaller-grid sequence arithmetic; requires an independently accepted complete scalar norm proof. Actual ratio limit is 1/N(b); displayed lower formulas tend to 1/B.',
        rows=[row(n,a.dimension_power,a.b,a.cap) for n in a.n])
    a.output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:data[k]for k in ('status','operator_b','required_global_norm_cap','limit_of_certified_lower_formula')},indent=2))
    print(json.dumps([{k:r[k]for k in ('n','dimension','integer_cutoff','finite_ratio_decimal')} for r in data['rows']],indent=2))
if __name__=='__main__':main()
