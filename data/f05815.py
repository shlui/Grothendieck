"""Directed outer-price rectangle, CONDITIONAL on the whole-q local gate.

The scalar centering theta is frozen at the fitting cumulant.  The nuisance
covariance correction is kept explicitly at every other cumulant.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse,json
from mpfr_interval import I,UP,DOWN
from certify_local_point import freeze,ROOT

def price(p,u,r,k,z,reserve,reserve0=None):
    theta=p['kappa']/(3*(p['kappa']+2));q0,q2,q4,q6=p['major'];q0+=reserve if reserve0 is None else reserve0;q6+=reserve
    A=q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    # These signs make the reduced price convex in both kappa and shape.
    assert A>=0 and B>=0
    K=I(194)/I(10).sqrt()-I(Q(633,9)*(z-Q(2,5)));assert I(0).lo<=K.lo
    kp=I(k)*I(k).sqrt();c6=K*kp;D=I(k*k)+I(Q(9,4))*I(k*z)*I(k*z).sqrt()
    Vc=k*(1+z)/9-2*theta*k/3+theta*theta*(k+2)
    covariance=k/3-theta*(k+2)
    # Transport y with r so y(r)-E02*r, and therefore the ENTIRE local
    # potential, remain exactly frozen.  No local Lipschitz reserve is needed.
    y=p['y']+p['p0']*p['p2']*(r-p['r'])
    val=I(u*u-2*p['ell']*u+6*p['beta']**2/p['eps']+(p['eps']-p['p0']**2)*r*r-2*p['gamma']*r+y*y/p['p2']**2)
    val+=I(q0+q2+(3+k)*q4+q6*(15+15*k))+I(q6)*c6
    val+=I(p['eta']*Vc/2+p['nuisance2']*covariance+p['zeta']**2*(6+9*k)-p['c']*k)+I(p['zeta']**2)*c6+I(B)*D
    return val

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('profile');ap.add_argument('--u',nargs=2,default=['.755','.755']);ap.add_argument('--r',nargs=2,default=['.25755','.25755']);ap.add_argument('--kappa',nargs=2,default=['.53','.53']);ap.add_argument('--reserve',default='.000001');ap.add_argument('--reserve0');ap.add_argument('--output');args=ap.parse_args();p=freeze(args.profile);rows=[]
    for u,r,k,z in product(*[[Q(t) for t in args.u],[Q(t) for t in args.r],[Q(t) for t in args.kappa],[Q(2,5),Q(1)]]):
        v=price(p,u,r,k,z,Q(args.reserve),Q(args.reserve0) if args.reserve0 else None);rows.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),lower=v.lo.text(DOWN,45),upper=v.hi.text(UP,45)))
    largest=max(rows,key=lambda t:Q(t['upper']));U=Q(largest['upper']);out=dict(status='CONDITIONAL_ON_COMPLETE_WHOLE_Q_LOCAL_GATE',scope='Directed outer rectangle, all shapes. The local gate must separately pass the compact q cover and exact q tail for this frozen profile and reserve.',rectangle=dict(u=args.u,r=args.r,kappa=args.kappa),reserve=args.reserve,reserve0=args.reserve0 or args.reserve,centering_theta=str(p['kappa']/(3*(p['kappa']+2))),maximum_upper=str(U),maximum_upper_decimal=largest['upper'],margin_to_40_71=str(Q(40,71)-U),rows=rows)
    out['r_transport']='y(r)=y0+E02*(r-r0); local potential stays exactly frozen.'
    path=Path(args.output or ROOT/'results'/'frozen_outer_price.json');path.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
