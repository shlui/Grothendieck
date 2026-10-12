"""Floating discovery-only general cubic thresholds for two reflected Y signs."""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm,erf
from scipy.optimize import minimize

def engine(nodes):
    x,w=roots_hermitenorm(nodes);w=w/np.sqrt(2*np.pi)
    z=x[:,None]+np.zeros((1,nodes));x=x[None,:]+np.zeros((nodes,1))
    weights=(w[:,None]*w[None,:]).ravel()
    basis=np.array([z,x,(z**3-3*z)/np.sqrt(6),(x**3-3*x)/np.sqrt(6),
      (z*z-1)*x/np.sqrt(2),z*(x*x-1)/np.sqrt(2)]).reshape(6,-1)
    def moments(t):
        q=t@basis;phi=np.exp(-q*q/2)/np.sqrt(2*np.pi)
        m=basis@(weights*erf(q/np.sqrt(2)));my=2*weights@phi
        jac=2*(basis*(weights*phi))@basis.T
        gy=-2*basis@(weights*q*phi)
        return m,my,jac,gy
    def evaluate(t,b,gradient=True):
        f,yf,jf,df=moments(t[:6]);g,yg,jg,dg=moments(t[6:])
        weight=np.array([1,1,-b,-b,-b,-b])
        value=np.dot(weight*f,g)+b*yf*yg
        if not gradient:return float(value),f.tolist(),g.tolist(),float(yf),float(yg)
        grad=np.r_[(weight*g)@jf+b*yg*df,(weight*f)@jg+b*yf*dg]
        return -float(value),-grad
    return evaluate

if __name__=='__main__':
    rng=np.random.default_rng(9031);e=engine(96);check=engine(256)
    base=np.array([3.73075,0,.961717,0,0,0]*2);rows=[]
    for b in [.72,.7171]:
        for i in range(30):
            initial=base.copy() if i==0 else base+rng.normal(0,np.tile([.5,.6,.4,.6,.4,.4],2))
            r=minimize(lambda t:e(t,b),initial,jac=True,method='L-BFGS-B',
             bounds=[(.1,12),(-6,6),(-5,5),(-5,5),(-5,5),(-5,5)]*2,
             options={'maxiter':400,'ftol':1e-14,'gtol':1e-9})
            value,f,g,yf,yg=check(r.x,b,False)
            row={'b':b,'start':i,'parameters':r.x.tolist(),'discovery':-r.fun,
                 'check':value,'f':f,'g':g,'Yf':yf,'Yg':yg,'success':bool(r.success)}
            rows.append(row)
            Path('tonight/adversarial_asymmetric_bivariate.json').write_text(json.dumps({
              'status':'DISCOVERY_ONLY','nodes':[96,256],'rows':rows},indent=2)+'\n')
        print(json.dumps({'b':b,'best':max(x['check'] for x in rows if x['b']==b)}),flush=True)
