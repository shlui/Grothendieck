"""Directed whole-line Gaussian maximum of rational cubics.
Numerical roots propose rational surrogates only. Polynomial discrepancies
are paid using a global L2 bound; no numerical root is assumed exact.
"""
from fractions import Fraction as Q
from functools import lru_cache
import numpy as np
from mpfr_interval import I, UP, DOWN

ZERO=Q(0); ROOT_DEN=10**12
C=(I(2)*I.pi()).sqrt().inv(); SQRT2=I(2).sqrt()

def trim(p):
    p=list(p)
    while len(p)>1 and not p[-1]:p.pop()
    return p

def ev(p,x):
    s=ZERO
    for a in reversed(p):s=s*x+a
    return s

def add(p,q):
    return trim([(p[i] if i<len(p) else ZERO)+(q[i] if i<len(q) else ZERO) for i in range(max(len(p),len(q)))])

def mul(p,q):
    r=[ZERO]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return trim(r)

def neg(p):return [-a for a in p]

def rat(x):return Q(int(round(float(x)*ROOT_DEN)),ROOT_DEN)

def coefficient_error(p,q):
    # ||z^j||_2 = 1,1,sqrt(3),sqrt(15), bounded by 1,1,2,4.
    d=add(p,neg(q));return sum(abs(x)*[1,1,2,4][j] for j,x in enumerate(d))

@lru_cache(maxsize=128)
def surrogate_cached(p):return surrogate(list(p))

def surrogate(p):
    p=trim(p);d=len(p)-1
    if d==0:return tuple(p),(),ZERO
    if d==1:
        return tuple(p),(-p[0]/p[1],),ZERO
    roots=np.roots([float(x) for x in reversed(p)])
    if d==2:
        disc=p[1]**2-4*p[0]*p[2]
        if disc<0:return tuple(p),(),ZERO
        if disc==0:return tuple(p),(-p[1]/(2*p[2]),),ZERO
        rr=sorted(rat(x.real) for x in roots)
    elif d==3:
        realroot=min(roots,key=lambda x:abs(x.imag)).real
        r=rat(realroot)
        # Exact division by z-r, with exact residual p(r).
        c=p[3];b=p[2]+c*r;a=p[1]+b*r
        if b*b-4*c*a<0:
            qq=list(p);err=ev(p,r);qq[0]-=err
            return tuple(qq),(r,),abs(err)
        rr=sorted(rat(x.real) for x in roots)
    else:raise ValueError('Cubics only')
    pp=[p[-1]]
    for r in rr:pp=mul(pp,[-r,Q(1)])
    return tuple(pp),tuple(sorted(set(rr))),coefficient_error(p,pp)

def primitives(x):
    if x is None:return (I(0),I(0),I(0),I(0))
    if x=='inf':return (I(1),I(0),I(0),I(0))
    z=I(x);phi=C*(-z.square()/2).exp()
    Phi=(I(1)+(z/SQRT2).erf())/2
    return (Phi,-phi,-z*phi,-(z.square()-1)*phi)

def integrate_power(p,diff):
    p=list(p)+[ZERO]*(4-len(p));h=[p[0]+p[2],p[1]+3*p[3],p[2],p[3]]
    return sum((I(x)*v for x,v in zip(h,diff) if x),I(0))

def max_integral(polys, detailed=False):
    p=[trim(x) for x in polys];n=len(p);sur={};roots=set();err=ZERO
    for i in range(n):
        for j in range(i+1,n):
            d=tuple(add(p[j],neg(p[i])))
            qq,rr,ee=surrogate_cached(d)
            sur[i,j]=qq;sur[j,i]=tuple(neg(qq));roots.update(rr);err+=ee
    rr=sorted(roots);bounds=[None]+rr+['inf']
    prim=[primitives(x) for x in bounds]
    lower=I(0);upper=I(0)
    for k in range(len(bounds)-1):
        a,b=bounds[k:k+2]
        if a is None and b=='inf':x=ZERO
        elif a is None:x=b-1
        elif b=='inf':x=a+1
        else:x=(a+b)/2
        i=max(range(n),key=lambda j:ev(p[j],x))
        qp=list(p[i])
        for j in range(n):
            if j!=i and ev(sur[i,j],x)>0:qp=add(qp,sur[i,j])
        diff=[v-u for u,v in zip(prim[k],prim[k+1])]
        lower+=integrate_power(p[i],diff)
        upper+=integrate_power(qp,diff)
    upper+=I(err)
    if detailed:return {'lower':lower.lo.text(DOWN,40),'upper':upper.hi.text(UP,40),'surrogate_L2_allowance':str(err),'intervals':len(bounds)-1}
    return upper

if __name__=='__main__':
    import time
    from scipy.integrate import quad
    tests=[[[Q(0),Q(1)],[Q(0),Q(-1)]],[[Q(-1),Q(0),Q(1)],[Q(1),Q(0),Q(-1)]],[[Q(-1,125),Q(-1,3),Q(0),Q(1,9)],[Q(1,5),Q(1,4),Q(1,7),Q(1,8)],[Q(0)]]]
    for pp in tests:
        t=time.time();z=max_integral(pp,True);print(z,'seconds',time.time()-t)
        f=lambda x:max(float(ev(q,Q(str(x)))) for q in pp)*np.exp(-x*x/2)/np.sqrt(2*np.pi)
        num=quad(f,-12,12,epsabs=1e-11,points=list(np.arange(-5,6)),limit=1000)[0];print('NUM',num)
