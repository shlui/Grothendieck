"""Directed local supremum for one rational q; no uniform operator certificate.

The even quadratic form is frozen to an exact rational rank-one matrix.  This
reduces four dual moment coordinates to three, with no numerical maximization
used for acceptance.  Product corner interpolation supplies the box bound.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import argparse, json, math, time, hashlib
from mpfr_interval import I, UP, DOWN, BITS, selftest
from gaussian_max_certificate import max_integral

ROOT=Path(__file__).resolve().parents[1]

def rational(x,places=11):
    return Q(f'{x:.{places}f}')

def validate(p):
    assert p['eps']>0 and p['kappa']>=0
    assert p['eta']>0 and p['G']>0
    assert p['p2']**2>=p['eps']*(1+p['p0']**2)/2
    assert p['p0']**2+2*p['p2']**2<1
    for a,b,c in [(p['R11'],p['R13'],p['R33']),
                  (p['R11']+1,p['R13'],p['R33']-p['eps']/6),
                  (p['R11']-p['eps'],p['R13'],p['R33'])]:
        assert a>=0 and c>=0 and a*c>=b*b
    return p

def freeze(path):
    d=json.loads(Path(path).read_text())
    if d.get('kind')=='EXACT_RATIONAL_R3_PROFILE':
        return validate({k:[Q(t) for t in v] if isinstance(v,list) else Q(v) for k,v in d['profile'].items()})
    x=d['history'][-1]['x']; n=d['order']
    assert not d.get('phi_quintic',False) and not d.get('mixed_mu',0)
    p0=rational(d['p'][0]); p2=rational(d['p'][1]/math.sqrt(2))+Q(1,10**10)
    eps=Q(str(.72 if '_e0.72_' in Path(path).name or '_e' not in Path(path).name else .7171))
    # Exact rank-one even matrix obeys both spectral and projection gates.
    assert p2*p2 >= eps*(1+p0*p0)/2
    r11=rational(x[5])+Q(1,10**8); r13=rational(x[6]/math.sqrt(6)); r33=rational(x[7]/6)+Q(1,10**8)
    for a,b,c in [(r11,r13,r33),(r11+1,r13,r33-eps/6),(r11-eps,r13,r33)]:
        assert a>=0 and c>=0 and a*c>=b*b
    ia=14+n
    ans=dict(ell=rational(x[0]),beta=rational(x[1]/math.sqrt(6)),
                p0=p0,p2=p2,R11=r11,R13=r13,R33=r33,
                gamma=rational(x[8]),y=rational(x[9]/math.sqrt(2)),r=Q(str(d['r'])),
                kappa=Q(str(d['kappa'])),eps=eps,
                major=[rational(t) for t in x[10:14]],
                phi=[rational(t) for t in x[14:14+n]],
                nuisance0=rational(x[ia]),nuisance2=rational(x[ia+1]),eta=rational(x[ia+2]),
                zeta=rational(x[ia+3]),c=rational(x[ia+4]),d0=rational(x[ia+5]),G=rational(x[ia+6]))
    return validate(ans)

def stein_features(q,n):
    z=q*q; t=(z-1)/(z+1); dt=4*q/(1+z)**2
    base=q/(1+z); db=(1-z)/(1+z)**2
    T=[Q(1)];D=[Q(0)]
    if n>1:T.append(t);D.append(Q(1))
    for j in range(2,n):
        T.append(2*t*T[-1]-T[-2]);D.append(2*T[-2]+2*t*D[-1]-D[-2])
    return [base*v for v in T],[db*v+base*w*dt for v,w in zip(T,D)]

def majorant(p,q):
    ph,dp=stein_features(q,len(p['phi'])); f=sum(a*b for a,b in zip(p['phi'],ph));fp=sum(a*b for a,b in zip(p['phi'],dp))
    hh=q*q-1;theta=p['kappa']/(3*(p['kappa']+2));low=-1-theta*hh
    res=fp-p['nuisance0']-p['nuisance2']*hh;eta=p['eta']
    pen=res*res/(2*eta) if res>=eta*low else low*res-eta*low*low/2
    return sum(a*q**(2*j) for j,a in enumerate(p['major']))+q*f-(1+theta*hh)*fp-pen

def context(p,q):
    S=q**3-3*q;g=p['c']*q+p['d0']*S;rho=1+p['zeta']**2/p['G']
    coeff=p['zeta']*(-S+g/(2*p['G']))
    b0=q*p['gamma']+coeff*p['p0'];b2=q*(p['y']-p['p0']*p['p2']*p['r'])+coeff*p['p2']
    return dict(rho=rho,b0=b0,b2=b2,constant=g*g/(4*p['G']))

def value(p,C,m):
    s,a,b=m
    h0=C['b0']+C['rho']*p['p0']*s;h2=C['b2']+C['rho']*p['p2']*s
    h1=p['R11']*a+p['R13']*b;h3=p['R13']*a+p['R33']*b
    B=[h0-h2,h1-3*h3,h2,h3]
    L=[Q(0),p['ell']-3*p['beta'],Q(0),p['beta']]
    polys=[L,[-x for x in L],B,[-x for x in B]]
    quad=C['rho']*s*s+p['R11']*a*a+2*p['R13']*a*b+p['R33']*b*b
    return max_integral(polys)*2-I(quad)+I(C['constant'])

def upper(p,C,box,cache):
    vals=[]
    for m in product(*[(a,b) for a,b in box]):
        if m not in cache:cache[m]=value(p,C,m)
        vals.append(cache[m])
    U=max(v.hi for v in vals)
    diag=[C['rho'],p['R11'],p['R33']]
    allowance=sum(d*(b-a)**2/4 for d,(a,b) in zip(diag,box))
    return I(U,U)+I(allowance)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('profile');ap.add_argument('--q',default='1');ap.add_argument('--inflate',default='0.00001');ap.add_argument('--max-nodes',type=int,default=100000);ap.add_argument('--output');args=ap.parse_args()
    p=freeze(args.profile);q=Q(args.q);C=context(p,q);cap=majorant(p,q)+Q(args.inflate)
    # |p0*k0+p2*k2| <= sqrt(p0²+2p2²)<1; |k1|<=1, |k3|<=sqrt(6)<5/2.
    assert p['p0']**2+2*p['p2']**2<1
    stack=[[(Q(-1),Q(1)),(Q(-1),Q(1)),(Q(-5,2),Q(5,2))]];cache={};events=[];count=accepted=0;start=time.time();largest=Q(-10**9)
    while stack and count<args.max_nodes:
        box=stack.pop();count+=1;v=upper(p,C,box,cache)
        if v.below(cap):
            accepted+=1;u=Q(v.hi.text(UP,30));largest=max(largest,u);events.append('C')
        else:
            diag=[C['rho'],p['R11'],p['R33']];axis=max(range(3),key=lambda j:diag[j]*(box[j][1]-box[j][0])**2)
            a,b=box[axis];mid=(a+b)/2;l=list(box);r=list(box);l[axis]=(a,mid);r[axis]=(mid,b);stack.extend([r,l]);events.append(str(axis))
        if count%2000==0:print(json.dumps(dict(nodes=count,accepted=accepted,stack=len(stack),corners=len(cache),seconds=time.time()-start)),flush=True)
    out=dict(status='PASS_FIXED_Q_LOCAL_SUPREMUM' if not stack else 'INCOMPLETE_NO_CERTIFICATE',scope='One rational q only. No q continuum, source-parameter cover or operator norm certificate.',q=str(q),cap=str(cap),inflate=args.inflate,bits=BITS,nodes=count,accepted=accepted,corners=len(cache),pending=len(stack),largest_accepted=str(largest),seconds=time.time()-start,interval_selftest=selftest(),profile={k:[str(t) for t in v] if isinstance(v,list) else str(v) for k,v in p.items()},events=events,source_sha256=hashlib.sha256(Path(args.profile).read_bytes()).hexdigest())
    output=Path(args.output) if args.output else ROOT/'results'/f'directed_local_point_q{args.q.replace("/","_")}.json';output.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['events','profile']},indent=2),flush=True)

if __name__=='__main__':main()
