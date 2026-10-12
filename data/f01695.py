#!/usr/bin/env python3
"""Independent endpoint assembly, without importing a producer verifier."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from mpfr_interval import I, UP
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
source=BASE/'tonight/operator/high_u_literature/source_profile.json'
summary=BASE/'tonight/recovery/profile_cover/high_u_literature_summary.json'
data=json.loads(source.read_text());local=json.loads(summary.read_text());p=local['profile'];raw=data['profile']
assert local['status']=='PASS' and len(local['nodes'])==51
assert local['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
assert all(x['status']=='PASS' and x['pending_boxes']==0 for x in local['nodes'])
M=[[F(z) for z in row] for row in p['M']]
for i in range(4):
 for j in range(4):assert M[i][j]==F(str(raw['M'][i][j]))+(F(1,10000) if i==j else 0)
for new,old in zip(p['nodal_majorant'],raw['nodal_majorant']):assert F(new)==F(str(old))+F(1,1000)
assert F(p['tail_slope'])==F(str(raw['tail_slope']))+F(1,10**6)
for key in ('ell','beta','delta','nu','epsilon'):assert F(p[key])==F(str(raw[key]))
ell,delta,nu,b=[F(p[k]) for k in ('ell','delta','nu','epsilon')];assert b==F(18,25)
gates=[]
for name,a,o,c in [('E',M[0][0],M[0][2],M[2][2]),('R',M[1][1],M[1][3],M[3][3]),('E-D13',M[0][0]+1,M[0][2],M[2][2]-b),('R-D13',M[1][1]+1,M[1][3],M[3][3]-b),('R-D2',M[1][1]-b,M[1][3],M[3][3])]:
 assert a>0 and c>0 and a*c-o*o>0;gates.append({'name':name,'determinant':str(a*c-o*o)})
assert F(p['tail_slope'])**2>=4*(delta**2+nu**2)
t=[F(str(x)) for x in data['t']];h=[F(str(x)) for x in raw['nodal_majorant']]
s=[(c-a)/(y-x) for x,y,a,c in zip(t[:-1],t[1:],h[:-1],h[1:])]
assert s[0]>=0 and all(a<=c for a,c in zip(s[:-1],s[1:])) and F(str(raw['tail_slope']))>=s[-1]
rho=F(849,1000);C=(I(2)/I.pi()).sqrt();lo,hi=[F(str(x)) for x in raw['amplitude_interval']]
grid=[F(str(x)) for x in data['r_grid']];assert grid[0]==0 and all(a<c for a,c in zip(grid[:-1],grid[1:]))
z=(I(lo)/C-(1-rho*rho))/(rho*rho);assert z.above(0) and z.below(1)
source_cap=I(rho)*(-z.log()).sqrt().erf();assert source_cap.below(grid[-1])
g=(nu-M[0][2])**2/M[2][2]-M[0][0]+b-2*delta
reserve=max(-g,F(0))*max((y-x)**2 for x,y in zip(grid[:-1],grid[1:]))/4
tailcharge=F(1,10**6)*rho*grid[-1]
ledger=[];best=F(p['nodal_majorant'][0]);active='zero_source';groups={r:[] for r in grid[1:]}
for index,rr in enumerate(data['rows']):
 cert=json.loads((HERE/'high_u_scalars'/f'continuous_scalar_row_{index:02d}.json').read_text());fr=cert['frozen_row']
 assert cert['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE' and cert['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
 assert F(fr['added_q0'])==F(1,5000) and F(fr['added_q6'])==F(1,10**9)
 r=F(fr['r']);assert r==F(str(rr['r'])) and r in groups
 interval=tuple(F(str(x)) for x in rr['kappa_interval']);assert tuple(map(F,fr['kappa_interval']))==interval;groups[r].append(interval)
 assert F(cert['tail']['domain'][0])>=max(t[-1]/r,F(12))
 slopes=[r*s[0]]+[r*(s[j]-s[j-1]) for j in range(1,len(s))]+[r*(F(str(raw['tail_slope']))-s[-1])]
 mean=I(h[0])
 for knot,weight in zip(t,slopes):
  v=I(knot/r);hinge=C*(-v*v/2).exp()-v*(1-(v/I(2).sqrt()).erf());mean+=I(weight)*hinge
 q=list(map(F,fr['majorant']));assert all(q[j]==F(str(rr['majorant'][j]))+(F(1,5000) if j==0 else F(1,10**9) if j==3 else 0) for j in range(4))
 mu,eta=F(fr['lambda_L1']),F(fr['eta']);assert mu==F(str(rr['lambda_L1'])) and eta==F(str(rr['eta'])) and mu>=0 and eta>0 and q[3]>0
 assert (I(194)/I(10).sqrt()-F(633,15)).above(0)
 for k in interval:
  for shape in (0,1):
   sixth=15+15*I(k)+(I(194)/I(10).sqrt()-F(633*shape,15))*I(k)*I(k).sqrt()
   price=mean+q[0]+q[1]+q[2]*(3+k)+I(q[3])*sixth+I(eta*k*(F(7,90)+F(shape,30)))+I(mu)*(I(rho)-C)
   value=price+g*r*r+F(1,1000)+tailcharge
   upper=F(value.hi.text(UP,180));ledger.append({'row':index,'r':str(r),'kappa':str(k),'shape':shape,'independent_price_upper':str(upper)})
   if upper>best:best=upper;active={'row':index,'kappa':str(k),'shape':shape}
for intervals in groups.values():
 intervals.sort();assert intervals[0][0]==0 and intervals[-1][1]==90 and all(a[1]==c[0] for a,c in zip(intervals[:-1],intervals[1:]))
outer=max(lo*lo-2*ell*lo,hi*hi-2*ell*hi);total=best+reserve+outer;assert total<F(9,16)<F(40,71)
out={'status':'PASS_INDEPENDENT_HIGH_U_ENDPOINT_ASSEMBLY','total_upper':str(total),'total_display':float(total),'rational_cap':'9/16','active_corner':active,'exact_gram_cost':str(g),'source_chord_reserve':str(reserve),'source_cap':source_cap.record(),'endpoint_ledger':ledger,'PSD_gates':gates,'local_producer_boxes':local['total_boxes'],'scope':'Regional operator bound only. Selected physical trees are independently replayed in separate artifact.'}
(HERE/'high_u_independent_endpoint_ledger.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','total_display','rational_cap','active_corner')},indent=2))
