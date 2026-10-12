"""Directed prices and exact tail for the general four-moment R3 profile."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,math
from general_R3_cover import freeze
from mpfr_interval import I,UP

def price(p,u,r,k,z,reserve0,reserve6):
    theta=p['kappa']/(3*(p['kappa']+2));q0,q2,q4,q6=p['major'];q0+=reserve0;q6+=reserve6
    A=q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    assert A>=0 and B>=0
    c6=(I(194)/I(10).sqrt()-I(Q(633,9)*(z-Q(2,5))))*I(k)*I(k).sqrt()
    D=I(k*k)+I(Q(9,4))*I(k*z)*I(k*z).sqrt()
    thetaVar=k*(1+z)/9-2*theta*k/3+theta**2*(k+2)
    cov=k/3-theta*(k+2)
    y=p['y']+p['E02']*(r-p['r'])
    rational=(u*u-2*p['ell']*u+6*p['beta']**2/p['eps']+
              (p['eps']-p['E00'])*r*r-2*p['gamma']*r+y*y/p['E22']+
              q0+q2+(3+k)*q4+q6*(15+15*k)+p['eta']*thetaVar/2+
              p['nuisance2']*cov+p['zeta']**2*(6+9*k)-p['c']*k)
    return I(rational)+I(A)*c6+I(B)*D

def tail(p,R,reserve0,reserve6):
    qmajor=list(p['major']);qmajor[0]+=reserve0;qmajor[3]+=reserve6
    Ap=sum(abs(c)for c in p['phi']);Ad=sum(j*j*abs(c)for j,c in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2));dp=Ap/R**2+4*Ad/R**4
    Dtheta=theta*Ap+Ap/R**2+4*theta*Ad/R**2+4*Ad/R**4
    A=abs(p['nuisance0']-p['nuisance2'])+dp;B=abs(p['nuisance2']);eta=p['eta']
    cc=p['c']-3*p['d0'];d=p['d0']
    l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0']
    l2=p['y']-p['E02']*p['r']+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    a3=abs(p['zeta']*(-1+d/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    trace=p['E00']+2*p['E22']+p['R11']+6*p['R33']+p['zeta']**2/p['G']*(p['p0']**2+2*p['p2']**2)
    constant=2*(abs(p['ell'])+Q(5,2)*abs(p['beta']))+trace
    pol=[Q(0)]*7
    for j,c in enumerate(qmajor):pol[2*j]=c
    pol[0]-=constant+Ap+Dtheta+A*A/(2*eta)
    pol[1]-=2*(abs(l1)+abs(l2));pol[2]-=A*B/eta+cc*cc/(4*p['G'])
    pol[3]-=2*a3;pol[4]-=B*B/(2*eta)+abs(cc*d)/(2*p['G']);pol[6]-=d*d/(4*p['G'])
    shifted=[sum(pol[j]*math.comb(j,k)*R**(j-k)for j in range(k,7))for k in range(7)]
    return dict(status='PASS'if all(x>=0 for x in shifted)else'FAIL',R=str(R),difference_polynomial=[str(x)for x in pol],shifted_coefficients=[str(x)for x in shifted])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('profile');ap.add_argument('--u',default='.7');ap.add_argument('--r',default='.389');ap.add_argument('--kappa',default='.53');ap.add_argument('--reserve0',default='.001');ap.add_argument('--reserve6',default='.000001');ap.add_argument('--tail',default='64');ap.add_argument('--output');args=ap.parse_args()
    p=freeze(args.profile);rows=[price(p,Q(args.u),Q(args.r),Q(args.kappa),z,Q(args.reserve0),Q(args.reserve6)).hi.text(UP,50)for z in [Q(2,5),Q(1)]]
    out=dict(status='DIRECTED_PRICE_LOCAL_GATE_SEPARATE',maximum_upper=str(max(Q(x)for x in rows)),tail=tail(p,Q(args.tail),Q(args.reserve0),Q(args.reserve6)),r_quadratic=str(p['eps']-p['E00']+p['E02']**2/p['E22']))
    if args.output:Path(args.output).write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
