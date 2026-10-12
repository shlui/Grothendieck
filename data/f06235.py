"""Directed point Gaussian integration of absolute normalized Hermite cubics.

Parameters are exact decimal rationals or Arb values. Root isolation is
certified by Arb; ambiguous/multiple roots raise, never pass silently.
"""
from flint import arb,acb,acb_poly,ctx
from fractions import Fraction
ctx.prec=256

def exact(x):
    if isinstance(x,arb):return x
    z=Fraction(str(x))
    return arb(z.numerator)/z.denominator

def abs_cubic(c):
    c=[exact(z) for z in c]
    power=[c[0]-c[2]/arb(2).sqrt(),c[1]-3*c[3]/arb(6).sqrt(),
           c[2]/arb(2).sqrt(),c[3]/arb(6).sqrt()]
    while len(power)>1 and power[-1].is_zero():power.pop()
    if len(power)==1:return abs(power[0])
    if not (power[-1]>0 or power[-1]<0):raise ArithmeticError('ambiguous degree')
    roots=acb_poly(power).roots(tol=arb('1e-60'),maxprec=512)
    real=[]
    for i,z in enumerate(roots):
        if not z.imag.contains(0):continue
        # A real-coefficient polynomial has conjugation-stable roots.
        # Its conjugate must lie in exactly this isolated root box.
        conj=z.conjugate()
        hits=[j for j,w in enumerate(roots) if w.overlaps(conj)]
        if hits!=[i]:raise ArithmeticError('ambiguous real-root classification')
        real.append(z.real)
    real.sort(key=lambda r:float(r.mid()))
    for a,b in zip(real[:-1],real[1:]):
        if not a<b:raise ArithmeticError('unordered root isolators')
    lead=1 if power[-1]>0 else -1
    sign=lead*((-1)**(len(power)-1))
    endsign=sign*((-1)**len(real))
    value=endsign*c[0]
    inv_sqrt_2pi=1/(2*arb.pi()).sqrt()
    for r in real:
        density=(-r*r/2).exp()*inv_sqrt_2pi
        cdf=(1+(r/arb(2).sqrt()).erf())/2
        primitive=c[0]*cdf-density*(c[1]+c[2]*r/arb(2).sqrt()+c[3]*(r*r-1)/arb(6).sqrt())
        value+=2*sign*primitive
        sign=-sign
    if not value>=0:raise ArithmeticError('uncertain nonnegative absolute integral')
    return value

def local_point(profile,t,x):
    t=exact(t);x=[exact(z) for z in x]
    M=[[exact(z) for z in row] for row in profile['M']]
    ell,beta,delta,nu,eps=[exact(profile[k]) for k in ['ell','beta','delta','nu','epsilon']]
    v=[sum(M[i][j]*x[j] for j in range(4)) for i in range(4)]
    v[0]+=t*delta;v[2]-=t*nu
    plus=[v[0],v[1]+ell,v[2],v[3]+beta]
    minus=[v[0],v[1]-ell,v[2],v[3]-beta]
    quadratic=sum(x[i]*M[i][j]*x[j] for i in range(4) for j in range(4))
    return abs_cubic(plus)+abs_cubic(minus)-quadratic+beta*beta/eps

if __name__=='__main__':
    import json,sys,numpy as np
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parent));import discover_middle_profiles as D
    p=json.load(open(sys.argv[1]))['profile'];rng=np.random.default_rng(73)
    worst=0
    for _ in range(100):
        t=float(rng.uniform(0,2.5));x=rng.uniform(-1,1,4)
        bound=local_point(p,t,x);nominal=D.N.local(p,t,x)[0]
        worst=max(worst,abs(float(bound.mid())-nominal))
        if abs(float(bound.mid())-nominal)>1e-10:raise AssertionError((bound,nominal))
    print(json.dumps({'status':'DIRECTED_POINT_BACKEND_CONSISTENCY_PASSED',
                     'points':100,'precision_bits':ctx.prec,'max_float_discrepancy':worst}))
