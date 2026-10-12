"""Discovery-only reflected bivariate cubic signs with diffuse third chaos.

The new transverse Gaussian coordinate is an original degree-one variable;
the independent Y channel is a CLT of original degree-three Hermites.
All six first/cubic polynomial coefficients and the diffuse cubic energy
are evaluated, preserving transverse parity in the reflected sign.
These floating computations do not certify an operator norm or endpoint.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm, erf
from scipy.optimize import minimize

def engine(nodes):
    x,w=roots_hermitenorm(nodes); w=w/np.sqrt(2*np.pi)
    z=x[:,None]+np.zeros((1,nodes)); x=x[None,:]+np.zeros((nodes,1))
    weights=(w[:,None]*w[None,:]).ravel()
    h2z=(z*z-1)/np.sqrt(2);h2x=(x*x-1)/np.sqrt(2)
    h3z=(z*z*z-3*z)/np.sqrt(6);h3x=(x*x*x-3*x)/np.sqrt(6)
    B=np.array([z,x,h3z,h3x,h2z*x,z*h2x]).reshape(6,-1)
    def evaluate(theta,b,grad=True):
        q=theta@B
        density=np.exp(-q*q/2)/np.sqrt(2*np.pi)
        m=B@(weights*erf(q/np.sqrt(2)))
        my=2*weights@density
        signs=np.array([1,-1,-b,b,b,-b])
        value=np.dot(signs,m*m)+b*my*my
        if not grad:return float(value),m.tolist(),float(my)
        dm=2*(B*(weights*density))@B.T
        dmy=-2*B@(weights*q*density)
        g=2*(signs*m)@dm+2*b*my*dmy
        return -float(value),-g
    return evaluate

if __name__=='__main__':
    rng=np.random.default_rng(714)
    e=engine(112);check=engine(224)
    records=[]
    base=np.array([2.555+3*.393,0,np.sqrt(6)*.393,0,0,0])
    for b in [.72,.7171]:
        for j in range(25):
            start=base.copy() if j==0 else base+rng.normal(0,[.8,.5,.4,.5,.4,.4])
            opt=minimize(lambda p:e(p,b),start,jac=True,method='L-BFGS-B',
                         bounds=[(.1,12),(-6,6),(-5,5),(-5,5),(-5,5),(-5,5)],
                         options={'maxiter':200,'ftol':2e-14,'gtol':1e-9})
            value,m,my=check(opt.x,b,False)
            row={'b':b,'start':j,'parameters':opt.x.tolist(),
                 'discovery_value':-float(opt.fun),'check_value':value,
                 'moments':m,'diffuse_cubic':my,'success':bool(opt.success)}
            records.append(row)
            print(json.dumps(row),flush=True)
            Path('tonight/adversarial_bivariate_signs.json').write_text(json.dumps({
              'status':'DISCOVERY_ONLY','nodes':[112,224],'records':records},indent=2)+'\n')
