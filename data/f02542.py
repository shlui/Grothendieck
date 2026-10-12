#!/usr/bin/env python3
"""Independent Arb Taylor integration of the cubic .7171 norm lower witness.

Uses no producer code or MPFR adapter. The signed ordinary-third compact
integral and all Gaussian tails are recomputed from exact rational data.
"""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_series,ctx
import hashlib,json,time

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
SOURCE=BASE/'reduced_search/operator/cubic7171/cubic7171_norm_witness.json'
A=F(26412631,10**7);B=F(4084047,10**7);EPS=F(7171,10000)
ctx.prec=512;ORDER=12;CELLS=400;L=F(5);ctx.cap=ORDER+1

def aa(x):
    x=F(x);return arb(x.numerator)/x.denominator
def jets(pos):
    z=arb_series([pos,1]);q=aa(A)*z+aa(B)*z**3
    e=(-(z*z+q*q)/2).exp();qprime=aa(A)+3*aa(B)*z*z
    return(e,qprime*e,(z*z-1)*qprime*e)
def interval(v):
    return{'lower':str(F(str(v.lower().fmpq()))),'upper':str(F(str(v.upper().fmpq())))}

def main():
    start=time.time();h=L/(2*CELLS);sums=[arb(0)]*3
    for k in range(CELLS):
        mid=F(2*k+1)*h;center=jets(aa(mid));whole=jets(aa(mid-h).union(aa(mid+h)))
        for j in range(3):
            sums[j]+=sum((center[j][i]*aa(2*h**(i+1)/F(i+1))for i in range(0,ORDER,2)),arb(0))
            sums[j]+=whole[j][ORDER]*aa(2*h**(ORDER+1)/F(ORDER+1))
    # q(z)^2>=A^2 z^2. Recurrence for Gaussian tail moments at c=1+A^2.
    c=aa(1+A*A);decay=(-c*aa(L*L)/2).exp();J0=decay/(c*aa(L))
    J2=aa(L)*decay/c+J0/c;J4=aa(L**3)*decay/c+3*J2/c
    factors=[2/arb.pi(),2/arb.pi(),2/(arb.pi()*arb(6).sqrt())]
    tails=[factors[0]*J0,factors[1]*(aa(A)*J0+3*aa(B)*J2),factors[2]*(aa(A)*(J2+J0)+3*aa(B)*(J4+J2))]
    moments=[factors[j]*sums[j]+arb(0).union(tails[j].upper())for j in range(3)]
    my,m1,m3=moments;P1=m1*m1;P3=m3*m3-my*my;value=P1-aa(EPS)*P3;inverse=1/value
    assert value>aa(F(5632335,10**7)) and value<aa(F(40,71))
    assert inverse>aa(F(177546,10**5)) and inverse<aa(F(177547,10**5))
    claimed=json.loads(SOURCE.read_text());assert claimed['status']=='PASS_DIRECTED_CLT_SCALAR_NORM_LOWER_WITNESS'
    assert {k:F(v)for k,v in claimed['parameters'].items()}=={'a':A,'c':B,'epsilon':EPS}
    for key,new in zip(['diffuse_cubic','first','ordinary_third'],moments):
        old=claimed['moments'][key];assert new.overlaps(aa(F(old['lower'])).union(aa(F(old['upper']))))
    for key,new in [('scalar_norm_lower_witness',value),('reciprocal_norm_ceiling',inverse)]:
        old=claimed[key];assert new.overlaps(aa(F(old['lower'])).union(aa(F(old['upper']))))
    result={'status':'PASS_INDEPENDENT_ARB_CUBIC7171_NORM_LOWER_WITNESS','producer_modules_imported':False,'bits':ctx.prec,'cells':CELLS,'Taylor_order':ORDER,'domain':['0',str(L)],'parameters':{'a':str(A),'c':str(B),'b':str(EPS)},'moments':{k:interval(v)for k,v in zip(['diffuse_cubic','first','ordinary_third'],moments)},'tails':list(map(interval,tails)),'P1_pairing':interval(P1),'P3_pairing':interval(P3),'scalar_norm_lower_witness':interval(value),'reciprocal_norm_ceiling':interval(inverse),'display':{'witness':float(value),'ceiling':float(inverse)},'seconds':time.time()-start,'producer_json_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'scope':'Limiting scalar norm lower test and reciprocal-norm method ceiling only. Does not prove any Grothendieck lower endpoint and does not exclude target1.775.'}
    (HERE/'cubic7171_witness_independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['status','display','bits','cells','Taylor_order','seconds']},indent=2))

if __name__=='__main__':main()
