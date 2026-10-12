#!/usr/bin/env python3
"""Independent direct-Chebyshev replay of the frozen cubic L1 Stein rows.

No code is imported from the discovery or archived verifier. Uses 768-bit Arb,
sixth-order Taylor remainder on [0,80], and independent Fraction polynomial
reconstruction of the tail with the wider rational Mills multiplier c<4/5.
Analytic Stein and sharp633 tensor premises are external to this replay.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
from functools import lru_cache
import hashlib
import json
import time
from flint import arb, arb_series, ctx

HERE = Path(__file__).resolve().parent
ctx.prec = 768
ctx.cap = 8
SPECS = [(0,F(0),F(39,100)), (2,F(39,100),F(39,25)),
         (3,F(39,25),F(4)), (4,F(4),F(34,5)),
         ('k8',F(34,5),F(48,5)), (5,F(48,5),F(76,5)),
         (6,F(76,5),F(28)), (7,F(28),F(90))]


def aa(q):
    q = F(q)
    return arb(q.numerator)/q.denominator


def trim(p):
    p = list(p)
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def add(p,q):
    z = [F(0)]*max(len(p),len(q))
    for i,c in enumerate(p): z[i] += c
    for i,c in enumerate(q): z[i] += c
    return trim(z)


def scale(p,c): return trim([x*c for x in p])


def mul(p,q):
    z=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): z[i+j] += a*b
    return trim(z)


def power(p,n):
    z=[F(1)]
    for _ in range(n): z=mul(z,p)
    return z


def derivative(p): return trim([F(j)*p[j] for j in range(1,len(p))] or [F(0)])


def translate(p,R):
    return [sum((p[j]*comb(j,i)*R**(j-i) for j in range(i,len(p))),F(0))
            for i in range(len(p))]


@lru_cache(None)
def mills(q):
    x=aa(q)
    return (arb.pi()/2).sqrt()*(x*x/2).exp()*(x/arb(2).sqrt()).erfc()


def direct_gap_jet(row, location, R):
    """Degree-six gap jet from normalized Chebyshev recurrence directly."""
    X=arb_series([location,1])
    L=aa(row['scale'])
    z=(X-L)/(X+L)
    vc=list(map(aa,row['phi_coefficients']))
    T=[arb_series([1]),z]
    for _ in range(2,row['degree']+1): T.append(2*z*T[-1]-T[-2])
    phi=sum((vc[j]*X**(j+1)/(j+1) for j in range(4)),arb_series([0]))
    phi += X/(X+L)*sum((vc[j+4]*T[j] for j in range(len(T))),arb_series([0]))
    dp=phi.derivative()
    rr=[R,location*R-1]
    for j in range(1,6): rr.append((location*rr[-1]+rr[-2])/(j+1))
    G=1-X*arb_series(rr)
    fhp=-(2/arb.pi()).sqrt()*G
    q=sum((aa(c)*X**(2*j) for j,c in enumerate(row['majorant'])),arb_series([0]))
    return q-dp+X*phi-(fhp-dp)**2/(2*aa(row['eta']))


def compact(row,R=F(80)):
    stack=[(F(0),R,0)]
    leaf=0; maximum_depth=0; min_gap=None
    while stack:
        l,r,d=stack.pop()
        mid=(l+r)/2; rad=(r-l)/2
        center=direct_gap_jet(row,aa(mid),mills(mid))
        whole=direct_gap_jet(row,aa(l).union(aa(r)),mills(r).union(mills(l)))
        gap=center[0]-sum((abs(center[j])*aa(rad**j) for j in range(1,6)),arb(0))
        gap -= abs(whole[6])*aa(rad**6)
        if gap>0:
            leaf+=1;maximum_depth=max(maximum_depth,d)
            low=F(str(gap.lower().fmpq()))
            min_gap=low if min_gap is None else min(min_gap,low)
        else:
            assert d<40,('COMPACT_FAILURE',str(l),str(r),str(gap))
            stack += [(mid,r,d+1),(l,mid,d+1)]
    return {'domain':['0',str(R)],'Taylor_order':6,'leaves':leaf,
            'maximum_depth':maximum_depth,'minimum_gap_rational':str(min_gap)}


def exact_tail(row,R=F(80)):
    """Independent Fraction reconstruction, no flint rational polynomials."""
    n=row['degree']; L=F(row['scale']); v=list(map(F,row['phi_coefficients']))
    d=[L,F(1)]; z=[-L,F(1)]; X=[F(0),F(1)]
    A=[[F(1)],z]
    for _ in range(2,n+1):
        A.append(add(scale(mul(z,A[-1]),F(2)),scale(mul(power(d,2),A[-2]),F(-1))))
    phi_polynomial=[F(0)]+[v[j]/(j+1) for j in range(4)]
    N=mul(phi_polynomial,power(d,n+1))
    for j in range(n+1): N=add(N,scale(mul(mul(X,A[j]),power(d,n-j)),v[j+4]))
    B=add(mul(derivative(N),d),scale(N,F(-n-1)))
    D=power(d,n+2)
    Q=[F(0)]*7
    for j,c in enumerate(row['majorant']): Q[2*j]=F(c)
    base=scale(mul(mul(power(X,4),D),
                   add(add(mul(Q,D),scale(B,F(-1))),mul(mul(X,N),d))),
               2*F(row['eta']))
    checks=[]
    for endpoint in (0,1):
        residual=add(scale(mul(power(X,2),B),F(-1)),scale(D,-endpoint*F(4,5)))
        P=add(base,scale(power(residual,2),F(-1)))
        shifted=translate(P,R)
        assert shifted[0]>0 and all(c>=0 for c in shifted),('TAIL_FAILURE',endpoint)
        checks.append({'endpoint':endpoint,'degree':len(P)-1,
                       'minimum_translated_coefficient':str(min(shifted)),
                       'polynomial_sha256':hashlib.sha256('|'.join(map(str,P)).encode()).hexdigest()})
    return {'domain':[str(R),'infinity'],'c_upper_rational':'4/5','checks':checks}


def prices(row,lo,hi):
    q=list(map(F,row['majorant'])); eta=F(row['eta'])
    assert q[3]>0 and eta>0
    K=194/arb(10).sqrt()
    assert K>aa(F(211,5))
    vals=[]
    for k in (lo,hi):
        for coef,hv in ((K,F(7,90)),(K-aa(F(211,5)),F(1,9))):
            value=(2/arb.pi()).sqrt()+aa(q[0]+q[1]+q[2]*(3+k))
            value += aa(q[3])*(15+15*aa(k)+coef*aa(k)*aa(k).sqrt())
            value += aa(eta*k*hv)
            vals.append(F(str(value.upper().fmpq())))
    maximum=max(vals)
    assert maximum < F(849,1000),('PRICE_FAILURE',float(maximum))
    return {'kappa_interval':[str(lo),str(hi)],'corner_upper_rational':str(maximum),
            'corner_upper_display':float(maximum)}


def arctan(q,n): return sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(n+1)),F(0))


if __name__=='__main__':
    # A separate exact Mills multiplier c=sqrt(2/pi)<4/5.
    pi_lower=16*arctan(F(1,5),41)-4*arctan(F(1,239),40)
    assert pi_lower>F(25,8)
    results=[]
    start=time.time()
    for i,lo,hi in SPECS:
        path=HERE/f'stein_L1_row_{i}.json'
        content=path.read_bytes(); row=json.loads(content)
        assert row['odd_weights']==['1','0','0'] and F(row['sixth_offset'])==0
        assert len(row['phi_coefficients'])==row['degree']+5
        one={'row':i,'row_sha256':hashlib.sha256(content).hexdigest(),
             'compact':compact(row),'tail':exact_tail(row),**prices(row,lo,hi)}
        results.append(one)
        print(json.dumps({'row':i,'compact_leaves':one['compact']['leaves'],
                          'price':one['corner_upper_display']}),flush=True)
    assert SPECS[0][1]==0 and SPECS[-1][2]==90
    assert all(SPECS[j][2]==SPECS[j+1][1] for j in range(len(SPECS)-1))
    maximum=max(F(r['corner_upper_rational']) for r in results)
    report={'status':'PASS_INDEPENDENT_COMPLETE_HALF_LINE_CUBIC_L1',
            'conditional_analytic_inputs':['fixed normalized third-chaos Stein identity',
                                          'sharp633 coupled sixth-cumulant theorem and variance range',
                                          '0<=kappa4<=90'],
            'precision_bits':ctx.prec,'target':'849/1000',
            'maximum_price_upper_rational':str(maximum),
            'maximum_price_display':float(maximum),'rows':results,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'seconds':time.time()-start}
    (HERE/'audit_stein_L1_independent.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','target','maximum_price_display','seconds']}))
