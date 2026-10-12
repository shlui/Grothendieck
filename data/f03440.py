#!/usr/bin/env python3
"""Independent exact replay of the freshly generated scalar633 cover.

Imports no producer code. Constructs each polynomial directly in that box's
local coordinates, from the TeX J/M formulas, then converts to Bernstein.
All arithmetic is Fraction in Q(sqrt(10)); there is no floating acceptance.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from math import comb
from pathlib import Path
import json,hashlib

HERE=Path(__file__).resolve().parent


@dataclass(frozen=True)
class Quadratic:
    rational:F=F(0)
    radical:F=F(0)
    def __add__(self,b):
        b=asq(b)
        return Quadratic(self.rational+b.rational,self.radical+b.radical)
    __radd__=__add__
    def __neg__(self):return Quadratic(-self.rational,-self.radical)
    def __sub__(self,b):return self+-asq(b)
    def __rsub__(self,b):return asq(b)+-self
    def __mul__(self,b):
        b=asq(b)
        return Quadratic(self.rational*b.rational+10*self.radical*b.radical,
                         self.rational*b.radical+self.radical*b.rational)
    __rmul__=__mul__
    def nonnegative(self):
        p,r=self.rational,self.radical
        if p>=0 and r>=0:return True
        if p<0 and r<=0:return False
        if p>=0:return p*p>=10*r*r
        return 10*r*r>=p*p
    def encoding(self):return [str(self.rational),str(self.radical)]


def asq(x):return x if isinstance(x,Quadratic) else Quadratic(F(x))


class Polynomial:
    def __init__(self,data):
        if isinstance(data,dict):self.data={k:asq(v) for k,v in data.items() if asq(v)!=Quadratic()}
        else:self.data={(0,0):asq(data)} if asq(data)!=Quadratic() else {}
    def __add__(self,b):
        b=asp(b); out=self.data.copy()
        for k,v in b.data.items():out[k]=out.get(k,Quadratic())+v
        return Polynomial(out)
    __radd__=__add__
    def __neg__(self):return Polynomial({k:-v for k,v in self.data.items()})
    def __sub__(self,b):return self+-asp(b)
    def __rsub__(self,b):return asp(b)+-self
    def __mul__(self,b):
        b=asp(b);out={}
        for (i,j),v in self.data.items():
            for (k,l),w in b.data.items():
                ij=(i+k,j+l);out[ij]=out.get(ij,Quadratic())+v*w
        return Polynomial(out)
    __rmul__=__mul__
    def __pow__(self,n):
        out=Polynomial(1)
        for _ in range(n):out=out*self
        return out


def asp(p):return p if isinstance(p,Polynomial) else Polynomial(p)


PAIRS=[(F(17,33),F(-1)),(F(62,105),F(-3,4)),
       (F(28,39),F(-1,2)),(F(3,4),F(-1))]


def cube_coefficients(p,q):
    cp=(F(55,2)-15*q)/(p*(p-q))
    cq=(F(55,2)-15*p)/(q*(q-p))
    zeta=F(135,2)-cp-cq
    beta=20+cp*p**3+cq*q**3
    assert cp>0 and cq>0 and zeta>=0
    assert 3*(cp*p+cq*q)==45
    assert 3*(cp*p*p+cq*q*q)==F(165,2)
    assert cp+cq+zeta==F(135,2)
    assert cp*p**3+cq*q**3-beta==-20
    return cp,cq,zeta,beta


def tangent(P,t):
    assert t>0
    return (P+t*t)*F(1,2)/t


# Polynomial division is needed only by rational scalars.
Polynomial.__truediv__=lambda self,t:self*F(1,F(t))


def reconstruct_local(leaf):
    xl,xh,zl,zh=map(F,leaf['rectangle'])
    y0,l0,tp,tq,tr,tg=map(F,leaf['parameters'])
    mode=leaf['mode']
    assert 0<=xl<xh<=1 and 0<=zl<zh<=1
    assert 0<=y0<=1 and l0>0 and l0*l0>=1+F(3,2)*y0*y0
    assert min(tp,tq,tr,tg)>0
    # Build directly in each rectangle's normalized Bernstein coordinates.
    y=Polynomial({(0,0):1-xl,(1,0):-(xh-xl)})
    a=Polynomial({(0,0):1-zl,(0,1):-(zh-zl)})
    d=1+F(3,2)*y**2
    bracket=d*Quadratic(F(422,15),F(97,5))-F(211,3)
    if mode=='corner':
        assert xh<=F(1,524288) and tp==F(7,4) and tr==1
        lower=(2+3*y)*Quadratic(F(0),F(1,10))
    else:
        lower=(1+F(3,2)*y0*y)/l0
    common_root=tangent((1+2*a**4)/3,tr)
    if mode=='corner' or mode.startswith('j'):
        idx=3 if mode=='corner' else int(mode[1:])
        assert 0<=idx<=3
        p,q=PAIRS[idx];cp,cq,zeta,beta=cube_coefficients(p,q)
        xp=p*p+(1+2*p)*y**2
        xq=q*q+(1+2*q)*y**2
        upper_p=xp*tangent(xp,tp)
        upper_q=xq/F(512) if mode=='corner' else xq*tangent(xq,tq)
        major=cp*upper_p+cq*upper_q+zeta*y**3-beta*a**3+F(135,2)*y*common_root
    else:
        common=F(135,2)*y**3+F(135,2)*y*common_root+45*a*y**2-20*a**3
        if mode=='m0':
            mixed=y*a*tangent(a,tg)
        else:
            assert mode=='m1'
            mixed=y**2*a**3+y*a**2*(tg*(1-y**2)+(1-a**2)/tg)/2
        major=common+F(165,2)*mixed
    result=lower*bracket-major
    assert all(i<=4 and j<=4 for i,j in result.data)
    return result


def bernstein_coefficients(local):
    result=[]
    for i in range(5):
        for j in range(5):
            c=Quadratic()
            for (k,l),v in local.data.items():
                if k<=i and l<=j:
                    c+=v*F(comb(i,k)*comb(j,l),comb(4,k)*comb(4,l))
            result.append(c)
    return result


if __name__=='__main__':
    for pair in PAIRS:cube_coefficients(*pair)
    assert cube_coefficients(F(3,4),F(-1))==(F(680,21),F(65,7),F(155,6),F(195,8))
    assert all(cube_coefficients(*pair)[2]==0 for pair in PAIRS[:3])
    # Exact positivity of the common RHS bracket.
    assert Quadratic(-F(211,5),F(97,5)).nonnegative()
    source=HERE/'regenerated_cover.json';payload=source.read_bytes();cert=json.loads(payload)
    rectangles=[];area=F(0);counts={};zeros=0;leaves=[]
    for index,leaf in enumerate(cert['leaves']):
        coefficients=bernstein_coefficients(reconstruct_local(leaf))
        assert all(c.nonnegative() for c in coefficients),('BERNSTEIN_FAIL',index,leaf)
        zeros+=sum(c==Quadratic() for c in coefficients)
        rect=tuple(map(F,leaf['rectangle']));rectangles.append(rect)
        area+=(rect[1]-rect[0])*(rect[3]-rect[2])
        counts[leaf['mode']]=counts.get(leaf['mode'],0)+1
        encoded=[c.encoding() for c in coefficients]
        leaves.append({'index':index,'mode':leaf['mode'],
                       'coefficient_sha256':hashlib.sha256(json.dumps(encoded).encode()).hexdigest()})
    assert area==1
    for index,(xl,xh,zl,zh) in enumerate(rectangles):
        for pl,ph,ql,qh in rectangles[:index]:
            assert min(xh,ph)<=max(xl,pl) or min(zh,qh)<=max(zl,ql)
    result={'status':'PASS_INDEPENDENT_FRACTION_QSQRT10_REPLAY',
            'source_sha256':hashlib.sha256(payload).hexdigest(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'rectangle_count':len(rectangles),'coefficient_count':25*len(rectangles),
            'zero_coefficient_count':zeros,'mode_counts':counts,'exact_area':str(area),
            'scope':'Complete scalar633 comparison from source J/M formulas; not a Wick-identity replay',
            'leaves':leaves}
    (HERE/'audit_cover_independent_fraction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='leaves'},indent=2))
