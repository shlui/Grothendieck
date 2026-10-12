"""Exact independent audit of rational Stein derivatives and Bernstein bounds."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json,sys,time

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'operator/projection/certificate'
sys.path.insert(0,str(CERT))
from certify_local_point import freeze
from rational_phi_bounds import numerators,bernstein_bounds as production_bound
from bernstein_integer_bounds import bernstein_bounds as fast


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b):out[j+k]+=x*y
    return out
def add(a,b):
    return [(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0)
            for j in range(max(len(a),len(b)))]
def derivative(a):return [j*a[j] for j in range(1,len(a))] or [0]
def evaluate(a,q):
    out=F(0)
    for c in a[::-1]:out=out*q+c
    return out


def features(order):
    A=[[1],[-1,0,1]]
    for j in range(2,order):
        A.append(add(mul([-2,0,2],A[-1]),
                     [-x for x in mul([1,0,2,0,1],A[-2])]))
    D=[1];out=[]
    for j in range(order):
        D=mul(D,[1,0,1]);out.append(([0]+A[j],D[:]))
    return out


def quotient_values(N,D,q,order=3):
    ns=[];ds=[]
    for _ in range(order+1):
        ns.append(evaluate(N,q));ds.append(evaluate(D,q))
        N=derivative(N);D=derivative(D)
    out=[]
    for j in range(order+1):
        out.append((ns[j]-sum(comb(j,k)*ds[k]*out[j-k]
                             for k in range(1,j+1)))/ds[0])
    return out


def reference(poly,a,b):
    """Independent direct affine power transform and basis conversion."""
    degree=len(poly)-1;w=b-a
    powers=[sum(poly[j]*comb(j,i)*a**(j-i)*w**i
                for j in range(i,degree+1)) for i in range(degree+1)]
    beta=[sum(powers[i]*F(comb(k,i),comb(degree,i))
              for i in range(k+1)) for k in range(degree+1)]
    return min(beta),max(beta)


def main():
    profile=ROOT/'operator/projection/results/shape_residual_n16_clip1_res1_u0.755_r0.25755_k0.53_pub0_pe0.15.json'
    p=freeze(profile);ns,n=numerators(tuple(p['phi']))
    fp=features(n);identities=0
    # After multiplying by (1+q^2)^(n+k), degree <=2n-1+k.
    # 2n+3 distinct zeros therefore prove all four rational identities.
    points=[F(j,n+1) for j in range(-(n+1),n+2)]
    assert len(points)>=2*n+3
    for q in points:
        expected=[F(0)]*4
        for c,(N,D) in zip(p['phi'],fp):
            for k,v in enumerate(quotient_values(N,D,q)):
                expected[k]+=c*v
        for k,poly in enumerate(ns):
            assert len(poly)-1<=2*n-1+k
            actual=evaluate(poly,q)/(1+q*q)**(n+k)
            assert actual==expected[k],(q,k)
            identities+=1

    intervals=[(F(0),F(1,10)),(F(1,5),F(21,100)),
               (F(1),F(11,10)),(F(3),F(31,10)),(F(8),F(9)),
               (F(-3),F(-1)),(F(7,53),F(77,529))]
    checks=0;slow=quick=0
    for a,b in intervals:
        for poly in ns:
            t=time.perf_counter();old=reference(poly,a,b);slow+=time.perf_counter()-t
            t=time.perf_counter();new=fast(poly,a,b);quick+=time.perf_counter()-t
            assert old==new,(a,b);checks+=1
            assert production_bound(poly,a,b)==new
    # Deterministic unrelated rational polynomials cover endpoints, signs,
    # cancellation, exact degeneracy, and degrees up to 40.
    extra=0
    for degree in [0,1,2,3,7,16,34,40]:
        poly=[F((-1)**j*(j*j+3),j+2) for j in range(degree+1)]
        for a,b in intervals+[(F(2,7),F(2,7))]:
            assert reference(poly,a,b)==fast(poly,a,b)
            assert production_bound(poly,a,b)==fast(poly,a,b)
            extra+=1
    out={'status':'PASS_EXACT_RATIONAL_IDENTITIES_AND_BERNSTEIN_RANGES',
         'derivative_orders':4,'chebyshev_order':n,
         'distinct_identity_points':len(points),
         'exact_identity_checks':identities,
         'active_range_checks':checks,'additional_range_checks':extra,
         'reference_seconds':slow,'integer_seconds':quick,
         'measured_range_speedup':slow/quick,
         'proof_scope':'All rational derivative identities by degree/zero count; exact Bernstein recurrences and identical ranges.',
         'source_hashes':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in [CERT/'rational_phi_bounds.py',
                                       Path(__file__).parent/'bernstein_integer_bounds.py',profile]}}
    path=Path(__file__).parent/'rational_phi_bounds_audit.json'
    path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
