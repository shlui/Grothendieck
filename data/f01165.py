"""Exact rational Stein-function derivatives and Bernstein interval bounds.

All floating fit coefficients have already been frozen to Fraction values.
The Chebyshev rational expression is expanded once; acceptance bounds use
exact Bernstein coefficients on each rational interval, not sampled values.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from mpfr_interval import I

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def add(p,r):return trim([(p[j] if j<len(p) else Q(0))+(r[j] if j<len(r) else Q(0)) for j in range(max(len(p),len(r)))])
def scale(p,c):return trim([c*t for t in p])
def mul(p,r):
    out=[Q(0)]*(len(p)+len(r)-1)
    for j,a in enumerate(p):
        for k,b in enumerate(r):out[j+k]+=a*b
    return trim(out)
def derivative(p):return trim([j*p[j] for j in range(1,len(p))] or [Q(0)])

@lru_cache(maxsize=16)
def numerators(coeff):
    n=len(coeff);D=[Q(1),Q(0),Q(1)];T=[Q(-1),Q(0),Q(1)];DD=mul(D,D)
    Ps=[[Q(1)]]
    if n>1:Ps.append(T)
    for j in range(2,n):Ps.append(add(scale(mul(T,Ps[-1]),2),scale(mul(DD,Ps[-2]),-1)))
    powers=[[Q(1)]]
    for j in range(n):powers.append(mul(powers[-1],D))
    numerator=[Q(0)]
    for j,c in enumerate(coeff):numerator=add(numerator,scale(mul(Ps[j],powers[n-j-1]),c))
    numerator=mul([Q(0),Q(1)],numerator)
    out=[numerator]
    for k in range(3):out.append(add(mul(D,derivative(out[-1])),scale(mul([Q(0),Q(2)],out[-1]),-(n+k))))
    return tuple(tuple(t) for t in out),n

def bernstein_bounds_reference(p,a,b):
    """The polynomial image is in the convex hull of these exact coefficients."""
    m=len(p)-1;w=b-a
    aa=[Q(1)];ww=[Q(1)]
    for j in range(m):aa.append(aa[-1]*a);ww.append(ww[-1]*w)
    power=[ww[i]*sum((p[j]*comb(j,i)*aa[j-i] for j in range(i,m+1)),Q(0)) for i in range(m+1)]
    Bs=[sum((power[i]*Q(comb(k,i),comb(m,i)) for i in range(k+1)),Q(0)) for k in range(m+1)]
    return min(Bs),max(Bs)

from bernstein_integer_bounds import bernstein_bounds

def phi_bounds(p,a,b):
    assert 0<=a<=b
    nums,n=numerators(tuple(p['phi']));out=[]
    for k,num in enumerate(nums):
        lo,hi=bernstein_bounds(num,a,b);dl=(1+a*a)**(n+k);dh=(1+b*b)**(n+k)
        interval=I(I(lo).lo,I(hi).hi)/I(I(dl).lo,I(dh).hi);out.append(interval)
    return tuple(out)

if __name__=='__main__':
    from certify_local_point import freeze
    from certify_local_region import phi_jet
    from pathlib import Path
    from mpfr_interval import UP,DOWN
    import json,time
    p=freeze(Path(__file__).resolve().parents[1]/'results'/'shape_residual_n16_clip1_res1_u0.755_r0.25755_k0.53_pub0_pe0.15.json')
    for a,b in [(Q(0),Q(1,10)),(Q(1,5),Q(21,100)),(Q(1),Q(11,10)),(Q(3),Q(31,10)),(Q(8),Q(9))]:
        start=time.time();vals=phi_bounds(p,a,b);jet=phi_jet(p,a,b)
        print(json.dumps(dict(a=str(a),b=str(b),seconds=time.time()-start,bernstein=[[v.lo.text(DOWN,12),v.hi.text(UP,12)]for v in vals],jet=[[v.lo.text(DOWN,12),v.hi.text(UP,12)]for v in jet])),flush=True)
