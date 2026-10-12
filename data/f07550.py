"""Independent exact algebra checks for the frozen R3 local certificate.

These checks audit formulas and directed derivative evaluations.  They do not
constitute a replay of the complete local cover or an outer-parameter cover.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
import hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'operator/projection/certificate'
sys.path.insert(0, str(CERT))
import certify_local_region as producer
from mpfr_interval import DOWN, UP


def add(a,b):
    return [ (a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0)
             for j in range(max(len(a),len(b))) ]


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c


def deriv(a):
    return [j*a[j] for j in range(1,len(a))] or [0]


def ev(a,q):
    v=F(0)
    for c in a[::-1]:v=v*q+c
    return v


def quotient_derivatives(num,den,q,order=3):
    ns=[];ds=[]
    for _ in range(order+1):
        ns.append(ev(num,q));ds.append(ev(den,q))
        num=deriv(num);den=deriv(den)
    fs=[]
    for j in range(order+1):
        fs.append((ns[j]-sum(comb(j,k)*ds[k]*fs[j-k]
                            for k in range(1,j+1)))/ds[0])
    return fs


def feature_polynomials(order):
    # T_j((q^2-1)/(q^2+1)) = A_j(q)/(q^2+1)^j.
    A=[[1],[-1,0,1]]
    for j in range(2,order):
        A.append(add(mul([-2,0,2],A[-1]),
                     [-x for x in mul([1,0,2,0,1],A[-2])]))
    D=[1]
    out=[]
    for j in range(order):
        D=mul(D,[1,0,1])
        out.append(([0]+A[j],D[:]))
    return out


NAMES=('q','theta','eta','inveta','res','rp','rpp',
       'phi','fp','fpp','fppp','k','V','gamma')
class Poly:
    """Small independent multivariate Fraction polynomial ring."""
    def __init__(self,x=0):
        self.c=dict(x.c) if isinstance(x,Poly) else {
            (0,)*len(NAMES):F(x)} if x else {}
    @staticmethod
    def var(name):
        x=Poly();m=[0]*len(NAMES);m[NAMES.index(name)]=1
        x.c[tuple(m)]=F(1);return x
    def __add__(self,x):
        x=Poly(x);z=Poly(self)
        for m,v in x.c.items():
            z.c[m]=z.c.get(m,F(0))+v
            if not z.c[m]:del z.c[m]
        return z
    __radd__=__add__
    def __neg__(self):
        z=Poly();z.c={m:-v for m,v in self.c.items()};return z
    def __sub__(self,x):return self+-Poly(x)
    def __rsub__(self,x):return Poly(x)+-self
    def __mul__(self,x):
        x=Poly(x);z=Poly()
        for m,v in self.c.items():
            for n,w in x.c.items():
                key=tuple(a+b for a,b in zip(m,n))
                z.c[key]=z.c.get(key,F(0))+v*w
        z.c={m:v for m,v in z.c.items() if v};return z
    __rmul__=__mul__
    def __pow__(self,n):
        z=Poly(1)
        for _ in range(n):z=z*self
        return z
    def __truediv__(self,x):return self*F(1,x)
    def derivative(self,rules):
        out=Poly()
        for m,v in self.c.items():
            for name,dx in rules.items():
                j=NAMES.index(name)
                if m[j]:
                    n=list(m);n[j]-=1;term=Poly();term.c[tuple(n)]=v*m[j]
                    out+=term*dx
        return out


def main():
    q,theta,eta,inveta,residual,rp,rpp,phi,fp,fpp,fppp,k,V,gamma=(
        Poly.var(n) for n in NAMES)
    dq={'q':1,'res':rp,'rp':rpp,'phi':fp,'fp':fpp,'fpp':fppp}
    low=-1-theta*(q*q-1)
    inside=(residual**2*inveta/2).derivative(dq).derivative(dq)
    assert not (inside-(rp**2+residual*rpp)*inveta).c
    outside=(low*residual-eta*low**2/2).derivative(dq).derivative(dq)
    wanted=(-2*theta*residual-4*theta*q*rp+low*rpp-
            eta*(4*theta**2*q*q-2*theta*low))
    assert not (outside-wanted).c
    stein=-q*phi+(1+theta*(q*q-1))*fp
    wanted=-2*fp-q*fpp+2*theta*fp+4*theta*q*fpp+(1+theta*(q*q-1))*fppp
    assert not (stein.derivative(dq).derivative(dq)-wanted).c

    order=16
    coeff=[F((-1)**j*(j+1),17*(j+2)) for j in range(order)]
    polynomials=feature_polynomials(order)
    points=[F(0),F(1,11),F(1,3),F(1),F(7,3),F(8),F(64)]
    checks=0
    for x in points:
        exact=[F(0)]*4
        for c,(num,den) in zip(coeff,polynomials):
            for j,val in enumerate(quotient_derivatives(num,den,x)):
                exact[j]+=c*val
        jet=producer.phi_jet({'phi':coeff},x,x)
        for j,val in enumerate(exact):
            lo=F(jet[j].lo.text(DOWN,90))
            hi=F(jet[j].hi.text(UP,90))
            assert lo<=val<=hi, (x,j,lo,val,hi)
            checks+=1

    cov=k/3-theta*(k+2)
    variance=V-2*theta*k/3+theta**2*(k+2)
    # E[U_theta phi'] = E[U_theta residual] + gamma*Cov.
    assert not ((gamma*cov).derivative({'gamma':1})-cov).c
    assert not (variance.derivative({'theta':1})+2*cov).c
    files=['certify_local_region.py','certify_local_point.py',
           'gaussian_max_certificate.py','mpfr_interval.py']
    out={'status':'PASS_FORMULA_AUDIT_ONLY',
         'scope':'Exact derivative and frozen-theta algebra; not cover replay.',
         'symbolic_second_derivatives':3,
         'exact_chebyshev_rational_derivative_enclosures':checks,
         'chebyshev_order':order,
         'tested_rational_q':[str(x) for x in points],
         'frozen_theta_nuisance_sign':'plus gamma2*Cov',
         'source_hashes':{f:hashlib.sha256((CERT/f).read_bytes()).hexdigest()
                          for f in files}}
    (Path(__file__).parent/'r3_local_symbolic_audit.json').write_text(
        json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
