#!/usr/bin/env python3
"""Independent 48-vertex full-E source-band price and curvature audit."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json
from mpfr_interval import I,UP
from audit_general_r3_parameters_tail import ROOT,DIRECTORY,load,gates

def price(p,u,r,k,z,reserve0,reserve6):
    theta=p['kappa']/(3*(p['kappa']+2))
    q0,q2,q4,q6=p['major'];q0+=reserve0;q6+=reserve6
    A=q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    assert A>=0 and B>=0
    theta_var=k*(1+z)/9-2*theta*k/3+theta*theta*(k+2)
    theta_cov=k/3-theta*(k+2)
    y=p['y']+p['E02']*(r-p['r'])
    rational=(u*u-2*p['ell']*u+6*p['beta']**2/p['eps']+
      (p['eps']-p['E00'])*r*r-2*p['gamma']*r+y*y/p['E22']+
      q0+q2+(3+k)*q4+(15+15*k)*q6+p['eta']*theta_var/2+
      p['nuisance2']*theta_cov+p['zeta']**2*(6+9*k)-p['c']*k)
    c6=(194/I(10).sqrt()-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt()
    D=I(k*k)+F(9,4)*I(k*z)*I(k*z).sqrt()
    return I(rational)+A*c6+B*D

def main():
    profile=DIRECTORY/'general_profile_u70_high_exact.json';p=load(profile);gates(p)
    ledgerpath=DIRECTORY/'general_source_polygons.json';ledger=json.loads(ledgerpath.read_text())
    audit=json.loads(Path(__file__).with_name('general_r3_u70_whole_q_audit.json').read_text())
    assert audit['status']=='PASS_INDEPENDENT_GENERAL_R3_WHOLE_Q_PROFILE'
    assert audit['profile_sha256']==hashlib.sha256(profile.read_bytes()).hexdigest()
    for name,h in audit['certificate_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    q0,q6=F(ledger['reserve0']),F(ledger['reserve6'])
    g=p['eps']-p['E00']+p['E02']**2/p['E22'];slope=F(9,5)
    A=p['major'][3]+q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    assert A>=0 and B>=0
    assert (194/I(10).sqrt()-F(633,9)*(1-F(2,5))).above(0)
    rows=[]
    for case in ledger['polygons']:
        ue,te,ke,ze=([F(x) for x in case[name]] for name in ('u','t','kappa','shape'))
        assert ze==[F(2,5),F(1)]
        reserve=max(F(0),-g)*(slope*slope*(ue[1]-ue[0])**2+(te[1]-te[0])**2)/4
        assert reserve==F(case['curvature_reserve']) and g==F(case['negative_r_quadratic'])
        maximum=F(-10**10);stored=[];vertices=[]
        expected=set(product(ue,te,ke,ze));assert len(case['vertices'])==16
        for v in case['vertices']:
            u,t,k,z=(F(v[name]) for name in ('u','t','kappa','shape'))
            assert (u,t,k,z) in expected;expected.remove((u,t,k,z))
            r=F(389,1000)-slope*(u-F(7,10))+t;assert r==F(v['r'])
            upper=F(price(p,u,r,k,z,q0,q6).hi.text(UP,65))
            assert upper<=F(v['upper'])
            maximum=max(maximum,upper+reserve);stored.append(F(v['upper']))
            vertices.append({'u':str(u),'t':str(t),'r':str(r),'kappa':str(k),
              'shape':str(z),'independent_directed_upper':str(upper)})
        assert not expected
        assert max(stored)+reserve==F(case['maximum_directed_upper'])
        cap=F(case['rounded_cap']);assert maximum<cap
        rows.append({'name':case['name'],'u':list(map(str,ue)),
          't':list(map(str,te)),'kappa':list(map(str,ke)),
          'negative_r_quadratic':str(g),'exact_curvature_reserve':str(reserve),
          'maximum_independent_upper':str(maximum),'cap':str(cap),
          'margin':str(cap-maximum),'vertices':vertices})
    out={'status':'PASS_INDEPENDENT_FULL_E_CONTINUOUS_SOURCE_BANDS',
      'scope':'Listed constrained-r bands only; no all-r amplitude coverage or global norm cap.',
      'profile_sha256':audit['profile_sha256'],
      'ledger_sha256':hashlib.sha256(ledgerpath.read_bytes()).hexdigest(),
      'whole_q_audit':'general_r3_u70_whole_q_audit.json',
      'vertex_checks':48,'bands':rows,
      'proof':'Product corner interpolation exactly preserves the mixed u*t term. Its r variance is slope^2 Var(u)+Var(t); negative g*r^2 costs at most (-g)*(slope^2 width_u^2+width_t^2)/4. The omitted positive u^2 variance only improves this upper bound. Source price is separately convex in kappa and shape.'}
    Path(__file__).with_name('general_r3_source_bands_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'bands':[
      {'name':r['name'],'cap':r['cap'],'upper':float(F(r['maximum_independent_upper']))}
      for r in rows]},indent=2))

if __name__=='__main__':main()
