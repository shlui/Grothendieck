#!/usr/bin/env python3
"""Independent exact six-gate and full-E tail audit (no producer imports).

This does not replay the compact four-moment cover or establish global
source/amplitude coverage. All arithmetic below is Python Fraction.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json

ROOT=Path(__file__).resolve().parents[2]
DIRECTORY=ROOT/'tonight/operator/projection/lower_amp_certificate'

def load(path):
    raw=json.loads(path.read_text())
    assert raw['kind']=='EXACT_RATIONAL_GENERAL_R3_PROFILE'
    return {k:[F(x) for x in v] if isinstance(v,list) else F(v)
            for k,v in raw['profile'].items()}

def gates(p):
    b,p0,p2=p['eps'],p['p0'],p['p2']
    assert b>0 and p['eta']>0 and p['G']>0 and p['kappa']>=0
    triples=[('E',p['E00'],p['E02'],p['E22']),
      ('E-longitudinal-one',p['E00']+1,p['E02'],p['E22']-b/2),
      ('E-ppT',p['E00']-p0*p0,p['E02']-p0*p2,p['E22']-p2*p2),
      ('R',p['R11'],p['R13'],p['R33']),
      ('R-transverse-zero',p['R11']+1,p['R13'],p['R33']-b/6),
      ('R-transverse-two',p['R11']-b,p['R13'],p['R33'])]
    result=[]
    for name,a,c,d in triples:
        det=a*d-c*c
        assert a>=0 and d>=0 and det>=0,(name,a,c,d,det)
        result.append({'gate':name,'first_diagonal':str(a),'second_diagonal':str(d),'determinant':str(det)})
    assert p['E22']>0
    return result

def tail(p,R,reserve0,reserve6):
    assert R>=1
    theta=p['kappa']/(3*(p['kappa']+2))
    Ap=sum(map(abs,p['phi']));Ad=sum(j*j*abs(c) for j,c in enumerate(p['phi']))
    dp=Ap/R**2+4*Ad/R**4
    Dtheta=theta*Ap+Ap/R**2+4*theta*Ad/R**2+4*Ad/R**4
    A=abs(p['nuisance0']-p['nuisance2'])+dp;B=abs(p['nuisance2'])
    cc=p['c']-3*p['d0'];dd=p['d0'];h=p['zeta']**2/p['G']
    l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0']
    l2=p['y']-p['E02']*p['r']+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    C3=abs(p['zeta']*(-1+dd/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    trace=p['E00']+2*p['E22']+p['R11']+6*p['R33']+h*(p['p0']**2+2*p['p2']**2)
    pol=[F(0)]*7
    for j,x in enumerate(p['major']):pol[2*j]=x
    pol[0]+=reserve0;pol[6]+=reserve6
    pol[0]-=2*abs(p['ell'])+5*abs(p['beta'])+trace+Ap+Dtheta+A*A/(2*p['eta'])
    pol[1]-=2*(abs(l1)+abs(l2))
    pol[2]-=A*B/p['eta']+cc*cc/(4*p['G'])
    pol[3]-=2*C3
    pol[4]-=B*B/(2*p['eta'])+abs(cc*dd)/(2*p['G'])
    pol[6]-=dd*dd/(4*p['G'])
    shifted=[sum(pol[j]*comb(j,k)*R**(j-k) for j in range(k,7)) for k in range(7)]
    assert all(x>=0 for x in shifted)
    return {'status':'PASS_EXACT_FULL_E_TAIL','R':str(R),
      'reserve0':str(reserve0),'reserve6':str(reserve6),'trace':str(trace),
      'difference_polynomial':list(map(str,pol)),'shifted_coefficients':list(map(str,shifted))}

def main():
    results=[]
    for label in ('70','68'):
        file=DIRECTORY/f'general_profile_u{label}_high_exact.json'
        p=load(file)
        t=tail(p,F(64),F(1,1000),F(1,10**6))
        rcurvature=p['eps']-p['E00']+p['E02']**2/p['E22']
        results.append({'profile':str(file.relative_to(ROOT)),
          'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),
          'gates':gates(p),'r_quadratic_coefficient':str(rcurvature),
          'requires_concave_r_chord_reserve':rcurvature<0,'tail':t})
        if label=='70':
            prior=json.loads((DIRECTORY/'general_u70_high_price_and_tail.json').read_text())
            assert prior['tail']['difference_polynomial']==t['difference_polynomial']
            assert prior['tail']['shifted_coefficients']==t['shifted_coefficients']
            assert prior['r_quadratic']==str(rcurvature)
    out={'status':'PASS_INDEPENDENT_SIX_GATES_AND_FULL_E_TAIL',
      'scope':'Analytic gates and q>=64 tail only; compact cover and global outer coverage are separate.',
      'profiles':results}
    Path(__file__).with_name('general_r3_parameters_tail_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'profiles':[
      {'profile':r['profile'],'tail':r['tail']['status'],
       'requires_concave_r_chord_reserve':r['requires_concave_r_chord_reserve']}
      for r in results]},indent=2))

if __name__=='__main__':main()
