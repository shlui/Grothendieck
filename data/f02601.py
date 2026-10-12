#!/usr/bin/env python3
"""Independent exact derivative and MPFR coefficient/source-price replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import re
from mpfr_interval import I,UP,BITS
BASE=Path(__file__).resolve().parents[2];src=BASE/'tonight/moments/cubic_head_contact/direct_contact_majorants.json'
d=json.loads(src.read_text());assert d['status']=='PASS_EXACT_ALL_HALFLINE_DIRECT_CUBIC_HEAD_CONTACT_MAJORANT'
def add(*polys):
    d={}
    for p in polys:
        for m,c in p.items():d[m]=d.get(m,F(0))+c
    return {m:c for m,c in d.items()if c}
def scale(p,c):return {m:v*c for m,v in p.items()if v*c}
def mul(p,q):
    out={}
    for m,c in p.items():
        for n,e in q.items():
            v=tuple(a+b for a,b in zip(m,n));out[v]=out.get(v,F(0))+c*e
    return {m:c for m,c in out.items()if c}
def derw(p):
    return { (m[0],m[1],m[2]-1,m[3]):c*m[2]for m,c in p.items()if m[2] }
def recorded_poly(text):
    out={};names=['A','B','w','tau']
    for term in text.split(' + '):
        coeff=F(1);powers=[0]*4
        for factor in term.split('*'):
            parts=factor.split('^')
            if parts[0]in names:powers[names.index(parts[0])]+=int(parts[1])if len(parts)>1 else 1
            else:coeff*=F(factor)
        m=tuple(powers);assert m not in out;out[m]=coeff
    return out
A={(1,0,0,0):F(1)};B={(0,1,0,0):F(1)};w={(0,0,1,0):F(1)};tau={(0,0,0,1):F(1)}
D1=add(A,mul(B,w));D3=add(A,scale(mul(B,w),3));P={(0,0,0,0):F(1)};proof=[]
for j in range(1,4):
    oldj=j-1
    P=add(mul(add(mul(tau,P),scale(derw(P),-1)),mul(D1,D3)),scale(mul(mul(B,P),D3),2*oldj+1),scale(mul(mul(B,P),D1),6*oldj))
    cs=list(P.values());assert all(c>0 for c in cs)
    recorded=d['derivative_positivity'][j-1]
    assert P==recorded_poly(recorded['positive_numerator'])
    assert len(cs)==recorded['coefficient_count'] and recorded['denominator_powers']==[2*j+1,2*j]
    proof.append({'order':j,'coefficient_count':len(cs),'minimum_coefficient':str(min(cs))})
a,b=F(d['parameters']['A']),F(d['parameters']['B']);assert a==F(26412631,10**7)and b==F(4084047,10**7)
assert F(d['parameters']['a'])==a+3*b
nu=(I(2)/I.pi()).sqrt();K=I(194)/I(10).sqrt();assert(K-F(633,15)).above(0)
rows=[];prices=0
for row in d['rows']:
    W=F(row['w']);s=W*(a+b*W)**2;assert s==F(row['contact_s'])
    x=I(W).sqrt()*(a+b*W);erf=(I(W)/2).sqrt().erf()
    J= x*erf+nu*(-I(W)/2).exp()*(a+b*(W+2));Jp=erf/(2*x)
    p0=nu*(a+2*b);p2=nu/(2*a);delta=J-p0-p2*s;de=Jp-p2
    p4=(3*delta-s*de)/(s*s);p6=(s*de-2*delta)/(s**3);assert p6.above(0)
    exact=[p0,p2,p4,p6];q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']))
    for Q,v in zip(q,exact):assert v.below(Q)
    assert q[3]>0
    maxima=[]
    for price in row['prices']:
        k=F(price['kappa']);z=F(price['shape']);assert 0<=k<=90 and 0<=z<=1
        value=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(K-F(633,15)*z)*I(k)*I(k).sqrt())
        assert value.below(F(price['upper']));prices+=1
        if k in(F(2),F(4)):maxima.append(F(price['upper']))
    assert max(maxima)==F(row['uniform_price_on_kappa_2_4_all_shapes'])
    rows.append({'contact_w':str(W),'coefficients':list(map(str,q)),'uniform_kappa2_4_price':row['uniform_price_on_kappa_2_4_all_shapes']})
out={'status':'PASS_INDEPENDENT_ALL_LINE_CUBIC_HEAD_CONTACT_GATE','input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'precision_bits':BITS,'positive_derivative_numerators':proof,'directed_coefficients_checked':4*len(rows),'directed_moment_prices_checked':prices,'rows':rows,'producer_code_imported':False,'analytic_review':['Gprime=nu/2 integral_0^1 exp(-wt^2/2)/(A+Bw) dt.','Exact recurrence proves Ktau derivatives signs through order3 and hence G fourth derivative nonpositive globally.','Repeated-node Hermite remainder applies outside contact interval; all-line majorant follows.','All coefficients already contain nu; componentwise rational upward rounding preserves comparison.','Physical source inequality uses diamond and independence of transverse Q and longitudinal Z.','Fixed nonnegative lambda credit has PLUS beta*d/b and retains local gates; no pointwise-min convexity inference permitted.'],'scope':'Accepted additional physical source gate and fixed-lambda credit. Does not alone prove a global norm bound.'}
path=Path(__file__).parent/'cubic_head_contact_independent.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ['rows','analytic_review']},indent=2))
