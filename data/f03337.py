"""Independent floating streamed GH and randomized Sobol candidate checks.

This does not import the discovery engine and does not constitute a rigorous
integration error certificate.  The complete sign-class scope is in README.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.special import erf, eval_hermitenorm, ndtri, roots_hermitenorm
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
INDEX = [(1,0,0), (0,1,0), (0,0,1), (3,0,0), (0,3,0), (0,0,3),
         (2,1,0), (2,0,1), (1,2,0), (0,2,1), (1,0,2), (0,1,2), (1,1,1)]


def basis(coords):
    hs = [[eval_hermitenorm(j, x)/math.sqrt(math.factorial(j)) for j in range(4)]
          for x in coords]
    return np.array([hs[0][a]*hs[1][b]*hs[2][c] for a,b,c in INDEX])


def weighted_moments(B, w, parameters, symmetric):
    fcoef = np.array(parameters[:13])
    gcoef = fcoef if symmetric else np.array(parameters[13:])
    fq, gq = fcoef @ B, gcoef @ B
    f = B @ (w*erf(fq/math.sqrt(2)))
    g = B @ (w*erf(gq/math.sqrt(2)))
    yf = math.sqrt(2/math.pi)*(w @ np.exp(-fq*fq/2))
    yg = math.sqrt(2/math.pi)*(w @ np.exp(-gq*gq/2))
    return np.r_[f, g, yf, yg]


def objective(moments, b):
    f, g, yf, yg = moments[:13], moments[13:26], moments[26], moments[27]
    return float(f[:3]@g[:3]-b*(f[3:]@g[3:])+b*yf*yg)


def gh(nodes, row):
    x, w = roots_hermitenorm(nodes)
    w /= math.sqrt(2*math.pi)
    x2, x3 = np.meshgrid(x, x, indexing='ij')
    x2, x3 = x2.ravel(), x3.ravel()
    weights23 = (w[:, None]*w[None, :]).ravel()
    total = np.zeros(28)
    for i in range(nodes):
        B = basis([np.full_like(x2,x[i]), x2, x3])
        total += weighted_moments(B, w[i]*weights23, row['parameters'], row['symmetric'])
    return {'method': 'streamed_GH', 'nodes': nodes, 'value': objective(total,row['b']),
            'moments': total.tolist()}


def sobol(power, seed, row):
    points = qmc.Sobol(3, scramble=True, seed=seed).random_base2(power)
    coords = ndtri(points).T
    total = np.zeros(28)
    block = 1 << 14
    for start in range(0, len(points), block):
        B = basis(coords[:, start:start+block])
        total += weighted_moments(B, np.full(B.shape[1],1/len(points)),
                                  row['parameters'], row['symmetric'])
    return {'method': 'randomized_Sobol', 'power': power, 'seed': seed,
            'value': objective(total,row['b'])}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nodes', type=int, nargs='+', default=[192,256,320])
    p.add_argument('--power', type=int, default=18)
    p.add_argument('--replicates', type=int, default=8)
    p.add_argument('--input', default='results.json')
    p.add_argument('--output', default='verification.json')
    args = p.parse_args()
    data = json.loads((HERE/args.input).read_text())
    selected=[]
    for b in [.72,.7171]:
        rows=[r for r in data['refined'] if r['b']==b]
        selected.append(max(rows,key=lambda r:r['checks'][-1]['value']))
    for b in [.72,.7171]:
        rows=[r for r in data['refined'] if r['b']==b]
        maximum=max(rows,key=lambda r:r['discovery_value'])
        if maximum not in selected:
            selected.append(maximum)
    output={'status':'DISCOVERY_ONLY', 'notes':'Separate recurrence basis; no rigorous quadrature errors.',
            'arguments':vars(args), 'rows':[]}
    for row in selected:
        entry={'b':row['b'],'coarse_start':row['coarse_start'], 'symmetric':row['symmetric'],
               'parameters':row['parameters'],'checks':[]}
        output['rows'].append(entry)
        for n in args.nodes:
            tic=time.monotonic(); result=gh(n,row); result['seconds']=time.monotonic()-tic
            entry['checks'].append(result)
            (HERE/args.output).write_text(json.dumps(output,indent=2)+'\n')
            print(json.dumps({'b':row['b'],'start':row['coarse_start'],
                              'nodes':n,'value':result['value'],'seconds':result['seconds']}),flush=True)
        entry['sobol']=[sobol(args.power,71021+i,row) for i in range(args.replicates)]
        values=[r['value'] for r in entry['sobol']]
        entry['sobol_summary']={'mean':float(np.mean(values)),
                                'sample_std':float(np.std(values,ddof=1))}
        (HERE/args.output).write_text(json.dumps(output,indent=2)+'\n')
    output['finished']=True
    (HERE/args.output).write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':
    main()
