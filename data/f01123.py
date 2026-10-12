"""Numerically nominate a local witness; certify it with a genuine integral lower.

The Gaussian upper expression's lower endpoint is never used as a lower
witness. Only integration of a chosen polynomial on rational intervals is.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,sys
import numpy as np
from scipy.optimize import minimize
from certify_local_point import freeze,context,majorant
from gaussian_max_certificate import max_integral
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from operator_gaussian_moments import moments


def main(a):
    p=freeze(a.profile);p['major'][0]+=F(a.reserve0);p['major'][3]+=F(a.reserve6)
    f={k:float(v) for k,v in p.items() if not isinstance(v,list)}
    sq2,sq6=math.sqrt(2),math.sqrt(6)
    rho=1+f['zeta']**2/f['G']
    best=(-math.inf,None,None)
    for q in np.arange(a.q_min,a.q_max+a.q_step/2,a.q_step):
        C=context(p,F(str(round(float(q),10))))
        c0,c2,cconst=map(float,[C['b0'],C['b2'],C['constant']])
        def fg(m):
            s,x,y=m
            coeff=np.array([c0+rho*f['p0']*s,f['R11']*x+f['R13']*y,
                sq2*(c2+rho*f['p2']*s),sq6*(f['R13']*x+f['R33']*y)])
            row=moments(f['ell'],f['beta']*sq6,coeff);k=row[2:6]
            quad=rho*s*s+f['R11']*x*x+2*f['R13']*x*y+f['R33']*y*y
            value=2*(f['ell']*row[0]+f['beta']*sq6*row[1]+coeff@k)-quad+cconst
            grad=2*np.array([rho*(f['p0']*k[0]+sq2*f['p2']*k[2]-s),
                f['R11']*(k[1]-x)+f['R13']*(sq6*k[3]-y),
                f['R13']*(k[1]-x)+f['R33']*(sq6*k[3]-y)])
            return -value,-grad
        seeds=[np.zeros(3)]+[t*np.eye(3)[j] for j in range(3) for t in [-.75,.75]]
        for seed in seeds:
            sol=minimize(fg,seed,jac=True,method='L-BFGS-B',bounds=[(-1,1),(-1,1),(-2.5,2.5)],options={'maxiter':150,'ftol':1e-14,'gtol':1e-10})
            gap=-fg(sol.x)[0]-float(majorant(p,F(str(round(float(q),10)))))
            if gap>best[0]:best=(gap,q,sol.x)
    gap,q,m=best
    q=F(str(round(float(q),10)))
    m=[F(f'{x:.10f}') for x in m]
    if m[1]<0:m[1],m[2]=-m[1],-m[2]
    C=context(p,q);s,x,y=m
    h0=C['b0']+C['rho']*p['p0']*s;h2=C['b2']+C['rho']*p['p2']*s
    h1=p['R11']*x+p['R13']*y;h3=p['R13']*x+p['R33']*y
    B=[h0-h2,h1-3*h3,h2,h3];L=[F(0),p['ell']-3*p['beta'],F(0),p['beta']]
    detail=max_integral([L,[-v for v in L],B,[-v for v in B]],detailed=True)
    quad=C['rho']*s*s+p['R11']*x*x+2*p['R13']*x*y+p['R33']*y*y
    lower=2*F(detail['lower'])-quad+C['constant']-majorant(p,q)
    physical_measure=s*s/(p['p0']**2+2*p['p2']**2)+x*x+y*y/6
    report=dict(status='TRUE_DIRECTED_LOCAL_COUNTEREXAMPLE' if lower>0 and physical_measure<=1 else 'FINITE_NUMERICAL_SCREEN_NO_COUNTEREXAMPLE',
        profile=str(Path(a.profile).resolve()),reserve0=a.reserve0,reserve6=a.reserve6,
        q=str(q),moment=list(map(str,m)),numerical_gap=gap,true_gap_lower=str(lower),true_gap_lower_display=float(lower),
        physical_measure=str(physical_measure),physical_measure_display=float(physical_measure),
        integral_details=detail,scope='Finite numerical nomination is discovery only. A positive genuine lower-integral gap at a rational physical point certifies that the indicated local gate fails. Absence of a positive lower gap proves no continuum statement.')
    Path(a.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','q','numerical_gap','true_gap_lower_display','physical_measure_display')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('profile');p.add_argument('--reserve0',required=True);p.add_argument('--reserve6',default='1/100000000')
    p.add_argument('--q-min',type=float,default=0);p.add_argument('--q-max',type=float,default=8);p.add_argument('--q-step',type=float,default=.05)
    p.add_argument('--output',required=True);main(p.parse_args())
