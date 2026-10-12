#!/usr/bin/env python3
"""Generate and replay a rational Bernstein proof of the scalar 633 bound.

Floating arithmetic chooses candidates only. Acceptance, subdivision,
coverage and every polynomial coefficient test use Fraction arithmetic
in Q(sqrt(10)). No recovered rectangle data are imported.
"""
from fractions import Fraction as F
from math import comb, isqrt, sqrt
from pathlib import Path
import json, sys, time

OUT = Path(__file__).with_name('regenerated_cover.json')

def q(v): return (F(v), F(0))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def neg(a): return (-a[0],-a[1])
def mul(a,b): return (a[0]*b[0]+10*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,r): return (a[0]*r,a[1]*r)
def nonnegative(a):
    p,r=a
    if r==0: return p>=0
    if p>=0 and r>=0: return True
    if p<=0 and r<=0: return False
    if p>=0: return p*p>=10*r*r
    return 10*r*r>=p*p
def scalarpoly(v): return {(0,0):q(v)}
def padd(*ps):
    out={}
    for p in ps:
        for e,c in p.items(): out[e]=add(out.get(e,q(0)),c)
    return {e:c for e,c in out.items() if c!=q(0)}
def pscale(p,r): return {e:scale(c,r) for e,c in p.items() if r}
def pmul(p,r):
    out={}
    for (i,j),c in p.items():
        for (k,l),v in r.items():
            e=(i+k,j+l);out[e]=add(out.get(e,q(0)),mul(c,v))
    return {e:c for e,c in out.items() if c!=q(0)}
def ppow(p,n):
    r=scalarpoly(1)
    for _ in range(n):r=pmul(r,p)
    return r

Y={(0,0):q(1),(1,0):q(-1)}
A={(0,0):q(1),(0,1):q(-1)}
Y2=ppow(Y,2); A2=ppow(A,2); A3=ppow(A,3); A4=ppow(A,4)
ONE=scalarpoly(1)
D=padd(ONE,pscale(Y2,F(3,2)))
ROOTARG=pscale(padd(ONE,pscale(A4,2)),F(1,3))
PAIRS=[(F(17,33),F(-1)),(F(62,105),F(-3,4)),
       (F(28,39),F(-1,2)),(F(3,4),F(-1))]
K=(F(0),F(97,5)) # 194/sqrt(10)
BRACKET=padd({e:mul(c,add(K,q(F(422,15)))) for e,c in D.items()},scalarpoly(F(-211,3)))

def rational_sqrt_upper(v, bits=36):
    den=1<<bits
    k=isqrt((v.numerator*den*den)//v.denominator)
    if F(k*k,den*den)<v:k+=1
    return F(k,den)
def tangent(x,t):
    assert t>0
    return pscale(padd(x,scalarpoly(t*t)),F(1,2)/t)
def cube_upper(x,t):return pmul(x,tangent(x,t))
def coefficients(p,r):
    cp=(F(55,2)-15*r)/(p*(p-r))
    cq=(F(55,2)-15*p)/(r*(r-p))
    return cp,cq,F(135,2)-cp-cq,20+cp*p**3+cq*r**3
def make_polynomial(mode,pars):
    # pars=(y0,l0,tp,tq,tr,tg), all rational; corner mode overrides L.
    y0,l0,tp,tq,tr,tg=pars
    if mode=='corner':
        L={e:mul(c,(F(0),F(1,10))) for e,c in padd(scalarpoly(2),pscale(Y,3)).items()}
    else:L=pscale(padd(ONE,pscale(Y,F(3,2)*y0)),1/l0)
    rhs=pmul(L,BRACKET)
    if mode.startswith('j') or mode=='corner':
        k=3 if mode=='corner' else int(mode[1:])
        p,r=PAIRS[k];cp,cq,zeta,beta=coefficients(p,r)
        xp=padd(scalarpoly(p*p),pscale(Y2,1+2*p))
        xq=padd(scalarpoly(r*r),pscale(Y2,1+2*r))
        uq=pscale(xq,F(1,512)) if mode=='corner' else cube_upper(xq,tq)
        major=padd(pscale(cube_upper(xp,tp),cp),pscale(uq,cq),
                   pscale(ppow(Y,3),zeta),pscale(A3,-beta),
                   pscale(pmul(Y,tangent(ROOTARG,tr)),F(135,2)))
    else:
        major=padd(pscale(ppow(Y,3),F(135,2)),
                   pscale(pmul(Y,tangent(ROOTARG,tr)),F(135,2)),
                   pscale(pmul(A,Y2),45),pscale(A3,-20))
        if mode=='m0':
            mixed=pmul(pmul(A,Y),tangent(A,tg))
        elif mode=='m1':
            x=padd(ONE,pscale(Y2,-1));z=padd(ONE,pscale(A2,-1))
            rootmajor=pscale(padd(pscale(x,tg),pscale(z,1/tg)),F(1,2))
            mixed=padd(pmul(Y2,A3),pmul(pmul(Y,A2),rootmajor))
        else:raise ValueError(mode)
        major=padd(major,pscale(mixed,F(165,2)))
    diff=padd(rhs,pscale(major,-1))
    assert all(i<=4 and j<=4 for i,j in diff)
    return diff

def bernstein(p,rect):
    xl,xh,zl,zh=rect;dx=xh-xl;dz=zh-zl
    local={}
    for (i,j),v in p.items():
        for k in range(i+1):
            for l in range(j+1):
                r=F(comb(i,k)*comb(j,l))*xl**(i-k)*zl**(j-l)*dx**k*dz**l
                e=(k,l);local[e]=add(local.get(e,q(0)),scale(v,r))
    bs=[]
    for i in range(5):
        for j in range(5):
            v=q(0)
            for (k,l),c in local.items():
                if k<=i and l<=j:
                    v=add(v,scale(c,F(comb(i,k),comb(4,k))*F(comb(j,l),comb(4,l))))
            bs.append(v)
    return bs
def fr(v):return f'{v.numerator}/{v.denominator}'
def encode(rect,mode,pars):
    return {'rectangle':[fr(v) for v in rect],'mode':mode,'parameters':[fr(v) for v in pars]}
def params(rect,mode):
    xl,xh,zl,zh=rect;y0=1-(xl+xh)/2;a0=1-(zl+zh)/2
    l0=rational_sqrt_upper(1+F(3,2)*y0*y0)
    tr=rational_sqrt_upper((1+2*a0**4)/3)
    if mode=='corner':return y0,l0,F(7,4),F(1),F(1),F(1)
    if mode.startswith('j'):
        p,r=PAIRS[int(mode[1:])]
        tp=rational_sqrt_upper(p*p+(1+2*p)*y0*y0)
        tq=rational_sqrt_upper(r*r+(1+2*r)*y0*y0)
        return y0,l0,tp,tq,tr,F(1)
    if mode=='m0':tg=rational_sqrt_upper(a0)
    else:
        # Minimum is attained at sqrt((1-a^2)/(1-y^2)).
        ratio=(1-a0*a0)/(1-y0*y0)
        tg=rational_sqrt_upper(ratio)
    return y0,l0,F(1),F(1),tr,tg
def value(mode,y,a):
    d=1+1.5*y*y
    rhs=sqrt(d)*((194/sqrt(10)+422/15)*d-211/3)
    if mode.startswith('j'):
        pair=PAIRS[int(mode[1:])]
        p,r=map(float,pair);cp,cq,zeta,beta=map(float,coefficients(*pair))
        bound=cp*max(0,p*p+(1+2*p)*y*y)**1.5+cq*max(0,r*r+(1+2*r)*y*y)**1.5+zeta*y**3-beta*a**3+67.5*y*sqrt((1+2*a**4)/3)
    else:
        mixed=a**1.5*y if mode=='m0' else y*y*a**3+y*a*a*sqrt(max(0,(1-y*y)*(1-a*a)))
        bound=67.5*y**3+67.5*y*sqrt((1+2*a**4)/3)+45*a*y*y-20*a**3+82.5*mixed
    return rhs-bound
def accept(rect):
    xl,xh,zl,zh=rect;y=float(1-(xl+xh)/2);a=float(1-(zl+zh)/2)
    modes=sorted(['j0','j1','j2','j3','m0','m1'],key=lambda m:value(m,y,a),reverse=True)
    if xl==0 and zl==0 and xh<=F(1,524288):modes=['corner']+modes
    for mode in modes:
        pars=params(rect,mode);poly=make_polynomial(mode,pars)
        bs=bernstein(poly,rect)
        if all(nonnegative(b) for b in bs):return encode(rect,mode,pars)
    return None
def replay(cert):
    leaves=cert['leaves'];rects=[];area=F(0);counts={};minimum=None
    for leaf in leaves:
        rect=tuple(map(F,leaf['rectangle']));pars=tuple(map(F,leaf['parameters']));mode=leaf['mode']
        xl,xh,zl,zh=rect
        assert 0<=xl<xh<=1 and 0<=zl<zh<=1
        y0,l0,tp,tq,tr,tg=pars
        assert l0>0 and l0*l0>=1+F(3,2)*y0*y0 and 0<=y0<=1
        assert tp>0 and tq>0 and tr>0 and tg>0
        if mode=='corner':assert xh<=F(1,524288) and tp==F(7,4) and tr==1
        bs=bernstein(make_polynomial(mode,pars),rect)
        assert all(nonnegative(v) for v in bs),(mode,rect)
        rects.append(rect);area+=(xh-xl)*(zh-zl);counts[mode]=counts.get(mode,0)+1
    # Every box stays in the unit square. Pairwise disjoint interiors
    # plus exact total area1 certifies coverage up to boundary points;
    # finitely many closed rectangles then cover the whole square.
    assert area==1
    for i,r in enumerate(rects):
        for s in rects[:i]:
            assert min(r[1],s[1])<=max(r[0],s[0]) or min(r[3],s[3])<=max(r[2],s[2])
    return {'status':'PASS','leaf_count':len(leaves),'bernstein_coefficients':25*len(leaves),'mode_counts':counts,'exact_area':str(area)}
def discover():
    # Isolate the equality corner, where ordinary tangent errors cannot
    # dominate a margin that vanishes exactly at (1,1).
    eps=F(1,524288)
    corner=(F(0),eps,F(0),eps)
    queue=[(corner,0),((eps,F(1),F(0),F(1)),0),((F(0),eps,eps,F(1)),0)]
    leaves=[];tries=0;start=time.time();maxdepth=0
    while queue:
        rect,depth=queue.pop();tries+=1;maxdepth=max(maxdepth,depth)
        leaf=accept(rect)
        if leaf:leaves.append(leaf)
        else:
            if depth>=32:raise RuntimeError(f'Failed depth32 at {rect}')
            xl,xh,zl,zh=rect
            if xh-xl>=zh-zl:
                mid=(xl+xh)/2;children=[(xl,mid,zl,zh),(mid,xh,zl,zh)]
            else:
                mid=(zl+zh)/2;children=[(xl,xh,zl,mid),(xl,xh,mid,zh)]
            queue.extend((r,depth+1) for r in children)
        if tries%100==0:print(json.dumps({'tries':tries,'leaves':len(leaves),'pending':len(queue),'maxdepth':maxdepth,'seconds':round(time.time()-start,1)}),flush=True)
    cert={'theorem':'scalar633comparison','field':'Q(sqrt(10))','degrees':[4,4],'generator':'discover_cover.py','leaves':leaves}
    cert['replay']=replay(cert);cert['discovery']={'tries':tries,'maxdepth':maxdepth,'seconds':time.time()-start}
    OUT.write_text(json.dumps(cert,indent=2)+'\n');print(json.dumps(cert['replay']),flush=True)

if __name__=='__main__':
    if '--replay' in sys.argv:print(json.dumps(replay(json.loads(OUT.read_text()))))
    else:discover()
