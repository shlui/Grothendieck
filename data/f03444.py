#!/usr/bin/env python3
"""Independent exact replay of the 633 comparison.

No producer modules, expanded polynomials, or Bernstein functions are
imported. The comparison is evaluated in factored form at a 5x5 rational
grid. A separately inverted degree-four Bernstein collocation matrix
recovers every coefficient by exact interpolation.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json


@dataclass(frozen=True)
class E:
    p: F=F(0)
    q: F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'p',F(self.p))
        object.__setattr__(self,'q',F(self.q))
    @staticmethod
    def as_e(x): return x if isinstance(x,E) else E(F(x))
    def __add__(self,x):
        x=self.as_e(x);return E(self.p+x.p,self.q+x.q)
    __radd__=__add__
    def __neg__(self): return E(-self.p,-self.q)
    def __sub__(self,x): return self+-self.as_e(x)
    def __rsub__(self,x): return self.as_e(x)+-self
    def __mul__(self,x):
        x=self.as_e(x)
        return E(self.p*x.p+10*self.q*x.q,self.p*x.q+self.q*x.p)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self.as_e(x);d=x.p*x.p-10*x.q*x.q
        if not d: raise ZeroDivisionError
        return self*E(x.p/d,-x.q/d)
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        v=E(1)
        for _ in range(n): v=v*self
        return v
    def nonnegative(self):
        if self.p==0: return self.q>=0
        if self.q==0: return self.p>=0
        if self.p>0 and self.q>0: return True
        if self.p<0 and self.q<0: return False
        if self.p>0: return self.p*self.p>=10*self.q*self.q
        return 10*self.q*self.q>=self.p*self.p
    def record(self): return {'rational':str(self.p),'sqrt10':str(self.q)}


def inverse(a):
    n=len(a)
    z=[list(row)+[F(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if z[i][j])
        z[j],z[pivot]=z[pivot],z[j]
        w=z[j][j];z[j]=[x/w for x in z[j]]
        for i in range(n):
            if i!=j:
                w=z[i][j]
                z[i]=[u-w*v for u,v in zip(z[i],z[j])]
    return [row[n:] for row in z]


NODES=[F(i,4) for i in range(5)]
COLLOCATION=[[F(comb(4,j))*t**j*(1-t)**(4-j) for j in range(5)] for t in NODES]
INV=inverse(COLLOCATION)
for i in range(5):
    for j in range(5):
        assert sum(INV[i][k]*COLLOCATION[k][j] for k in range(5))==int(i==j)

PAIRS=[(F(17,33),F(-1)),(F(62,105),F(-3,4)),
       (F(28,39),F(-1,2)),(F(3,4),F(-1))]
K=E(0,F(97,5))


def value(x,z,mode,parameters):
    """Factored RHS-minus-upper-bound; no expanded producer polynomial."""
    y,a=1-x,1-z
    y0,l0,tp,tq,tr,tg=parameters
    d=1+F(3,2)*y*y
    if mode=='corner': L=E(0,F(2+3*y,10))
    else: L=E((1+F(3,2)*y0*y)/l0)
    rhs=L*((K+F(422,15))*d-F(211,3))
    root_arg=(1+2*a**4)/3
    root_bound=(root_arg+tr*tr)/(2*tr)
    cross=F(135,2)*y*root_bound
    if mode in ['j0','j1','j2','j3','corner']:
        p,q=PAIRS[3 if mode=='corner' else int(mode[1])]
        cp=(F(55,2)-15*q)/(p*(p-q))
        cq=(F(55,2)-15*p)/(q*(q-p))
        zeta=F(135,2)-cp-cq
        beta=20+cp*p**3+cq*q**3
        assert cp>=0 and cq>=0 and zeta>=0
        Xp=p*p+(1+2*p)*y*y
        Xq=q*q+(1+2*q)*y*y
        assert Xp>=0 and Xq>=0
        first=Xp*(Xp+tp*tp)/(2*tp)
        second=Xq/512 if mode=='corner' else Xq*(Xq+tq*tq)/(2*tq)
        major=cp*first+cq*second+zeta*y**3-beta*a**3+cross
    elif mode in ['m0','m1']:
        if mode=='m0': mixed=y*a*(a+tg*tg)/(2*tg)
        else:
            gram_root=(tg*(1-y*y)+(1-a*a)/tg)/2
            mixed=y*y*a**3+y*a*a*gram_root
        major=F(135,2)*y**3+cross+45*a*y*y-20*a**3+F(165,2)*mixed
    else: raise ValueError(mode)
    return rhs-major


def recover_coefficients(rect,mode,parameters):
    xl,xh,zl,zh=rect
    grid=[[value(xl+(xh-xl)*u,zl+(zh-zl)*v,mode,parameters)
           for v in NODES] for u in NODES]
    return [[sum((grid[k][l]*(INV[i][k]*INV[j][l])
                  for k in range(5) for l in range(5)),E())
             for j in range(5)] for i in range(5)]


def bernstein_value(coeffs,u,v):
    return sum((coeffs[i][j]*F(comb(4,i)*comb(4,j))*
                u**i*(1-u)**(4-i)*v**j*(1-v)**(4-j)
                for i in range(5) for j in range(5)),E())


def main():
    path=Path(__file__).with_name('regenerated_cover.json')
    raw=path.read_bytes();certificate=json.loads(raw)
    leaves=certificate['leaves']
    boxes=[];records=[];counts={};zeros=0;minimum_positive=None
    for number,leaf in enumerate(leaves):
        rect=tuple(F(s) for s in leaf['rectangle'])
        parameters=tuple(F(s) for s in leaf['parameters'])
        assert len(rect)==4 and len(parameters)==6
        xl,xh,zl,zh=rect
        y0,l0,tp,tq,tr,tg=parameters
        assert 0<=xl<xh<=1 and 0<=zl<zh<=1
        assert 0<=y0<=1 and l0>0 and l0*l0>=1+F(3,2)*y0*y0
        assert min(tp,tq,tr,tg)>0
        mode=leaf['mode']
        if mode=='corner':
            assert xh<=F(1,524288)
            assert tp==F(7,4) and tr==1
            # For every point in the box: 1-y²<=2x<=1/262144.
            assert 2*xh<=F(1,262144)
        bs=recover_coefficients(rect,mode,parameters)
        assert all(b.nonnegative() for row in bs for b in row),(number,mode)
        for row in bs:
            for b in row:
                if b==E(): zeros+=1
                elif minimum_positive is None or (minimum_positive-b).nonnegative():
                    minimum_positive=b
        # Extra exact evaluation checks independently exercise interpolation.
        for u,v in [(F(1,7),F(2,7)),(F(3,7),F(5,7))]:
            assert bernstein_value(bs,u,v)==value(xl+(xh-xl)*u,zl+(zh-zl)*v,mode,parameters)
        boxes.append(rect);counts[mode]=counts.get(mode,0)+1
        records.append({'leaf':number,'mode':mode,
                        'bernstein_coefficients':[[b.record() for b in row] for row in bs]})
    area=sum((xh-xl)*(zh-zl) for xl,xh,zl,zh in boxes)
    assert area==1
    for j,(xl,xh,zl,zh) in enumerate(boxes):
        for ul,uh,vl,vh in boxes[j+1:]:
            assert min(xh,uh)<=max(xl,ul) or min(zh,vh)<=max(zl,vl)
    out={'status':'PASS_INDEPENDENT_EXACT_COLLOCATION_REPLAY',
         'source_sha256':hashlib.sha256(raw).hexdigest(),
         'method':'factored exact formulas; inverse Bernstein collocation; no producer imports',
         'field':'Q(sqrt(10))','rectangles':len(boxes),'coefficients':25*len(boxes),
         'zero_coefficients':zeros,'minimum_positive_coefficient':minimum_positive.record(),
         'exact_area':str(area),'mode_counts':counts,'records':records}
    target=Path(__file__).with_name('independent_replay.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))


if __name__=='__main__': main()
