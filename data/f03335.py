"""Discovery-only 3D asymmetric cubic thresholds plus one degree-three CLT channel.

For S=P1-bP3, use f=sgn(qf(X)+Y), g=sgn(qg(X)-Y).
X consists of three independent original Gaussians; Y is an independent
Gaussian limit of normalized diffuse original-degree-three chaoses.  Y is
integrated analytically.  All quadrature and optimization here are floating.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import erf, roots_hermitenorm

HERE = Path(__file__).resolve().parent
MULTIINDICES = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
                (3, 0, 0), (0, 3, 0), (0, 0, 3),
                (2, 1, 0), (2, 0, 1), (1, 2, 0),
                (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]
C = math.sqrt(2 / math.pi)
TARGET = 40 / 71


class Engine:
    def __init__(self, nodes):
        self.nodes = int(nodes)
        x, w = roots_hermitenorm(self.nodes)
        w /= math.sqrt(2 * math.pi)
        coords = np.meshgrid(x, x, x, indexing='ij')
        xx = [z.ravel() for z in coords]
        self.weights = (w[:, None, None] * w[None, :, None] * w[None, None, :]).ravel()
        h = [[np.ones_like(z), z, (z*z - 1)/math.sqrt(2),
              (z*z*z - 3*z)/math.sqrt(6)] for z in xx]
        self.basis = np.ascontiguousarray([h[0][a] * h[1][b] * h[2][c]
                                          for a, b, c in MULTIINDICES])

    def moments(self, coeff):
        q = coeff @ self.basis
        phi = np.exp(-q*q/2) / math.sqrt(2*math.pi)
        m = self.basis @ (self.weights * erf(q/math.sqrt(2)))
        y = 2 * (self.weights @ phi)
        return m, float(y), q, phi

    def evaluate(self, coeff, b, symmetric=False, gradient=True):
        f, yf, qf, phif = self.moments(coeff[:13])
        g, yg, qg, phig = ((f, yf, qf, phif) if symmetric else
                           self.moments(coeff[13:]))
        weights = np.r_[np.ones(3), -b*np.ones(10)]
        value = float((weights*f) @ g + b*yf*yg)
        if not gradient:
            return {'value': value, 'first_pairing': float(f[:3] @ g[:3]),
                    'cubic_pairing': float(f[3:] @ g[3:]),
                    'diffuse_pairing': float(yf*yg), 'f_moments': f.tolist(),
                    'g_moments': g.tolist(), 'Yf': yf, 'minus_Yg': yg}
        # Direct matvec gradient avoids forming the 13x13 moment Jacobian.
        gradf = 2 * self.basis @ (self.weights * phif *
                                 ((weights*g) @ self.basis - b*yg*qf))
        if symmetric:
            grad = 2 * gradf
        else:
            gradg = 2 * self.basis @ (self.weights * phig *
                                     ((weights*f) @ self.basis - b*yf*qg))
            grad = np.r_[gradf, gradg]
        return -value, -grad


def rotated_seed(v, linear=3.73075, cubic=.961717):
    v = np.asarray(v, dtype=float)
    v /= np.linalg.norm(v)
    out = np.zeros(13)
    out[:3] = linear*v
    out[3:6] = cubic*v**3
    for j, alpha in enumerate(MULTIINDICES[6:12], 6):
        i2 = alpha.index(2)
        i1 = alpha.index(1)
        out[j] = cubic*math.sqrt(3)*v[i2]**2*v[i1]
    out[12] = cubic*math.sqrt(6)*np.prod(v)
    return out


def optimize(engine, coeff, b, symmetric, maxiter):
    tic = time.monotonic()
    result = minimize(lambda t: engine.evaluate(t, b, symmetric), coeff,
                      jac=True, method='L-BFGS-B', bounds=[(-12., 12.)]*len(coeff),
                      options={'maxiter': maxiter, 'ftol': 2e-14, 'gtol': 2e-9,
                               'maxls': 40})
    return {'b': b, 'symmetric': symmetric, 'nodes': engine.nodes,
            'parameters': result.x.tolist(), 'discovery_value': -float(result.fun),
            'success': bool(result.success), 'message': str(result.message),
            'iterations': int(result.nit), 'seconds': time.monotonic()-tic}


def selftest():
    e = Engine(12)
    gram = (e.basis*e.weights) @ e.basis.T
    rng = np.random.default_rng(8172)
    t = np.r_[rotated_seed([1, 0, 0]), rotated_seed([1, 0, 0])] + rng.normal(0, .08, 26)
    value, gradient = e.evaluate(t, .72)
    fd = np.array([(e.evaluate(t+np.eye(26)[j]*1e-6, .72)[0] -
                    e.evaluate(t-np.eye(26)[j]*1e-6, .72)[0])/2e-6 for j in range(26)])
    return {'gram_max_error': float(np.max(abs(gram-np.eye(13)))),
            'gradient_max_error': float(np.max(abs(fd-gradient))),
            'value': -value}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nodes', type=int, default=32)
    p.add_argument('--starts', type=int, default=24)
    p.add_argument('--refine-nodes', type=int, default=64)
    p.add_argument('--refine-top', type=int, default=3)
    p.add_argument('--maxiter', type=int, default=400)
    p.add_argument('--checks', type=int, nargs='+', default=[80, 112, 144])
    p.add_argument('--b', type=float, nargs='+', default=[.72, .7171])
    p.add_argument('--output', default='results.json')
    p.add_argument('--diverse', action='store_true',
                   help='After seeded starts, explore independent and cubic-antipodal thresholds.')
    args = p.parse_args()
    output = HERE / args.output
    rng = np.random.default_rng(8301904)
    data = {'status': 'DISCOVERY_ONLY', 'target': TARGET,
            'basis_multiindices': MULTIINDICES, 'selftest': selftest(),
            'arguments': vars(args), 'coarse': [], 'refined': []}
    assert data['selftest']['gram_max_error'] < 1e-12
    assert data['selftest']['gradient_max_error'] < 1e-7
    def save():
        output.write_text(json.dumps(data, indent=2)+'\n')
    e = Engine(args.nodes)
    for b in args.b:
        for j in range(args.starts):
            symmetric = j % 3 == 0
            direction = np.array([1., 0., 0.]) if j < 3 else rng.normal(size=3)
            seed = rotated_seed(direction)
            scale = [0., .12, .3, .7][j % 4]
            if symmetric:
                initial = seed + rng.normal(0, scale, 13)
            else:
                initial = np.r_[seed, seed] + rng.normal(0, scale, 26)
            if args.diverse and j >= 12:
                # These starts also enter sectors far from the known reflected
                # rank-one extremizer.  No quadrant or sign bound is imposed.
                amplitude = [.3, .8, 1.5, 3.][j % 4]
                first = rng.normal(0, amplitude, 13)
                first[:3] *= [0., 1., 2.][j % 3]
                if symmetric:
                    initial = first
                else:
                    second = rng.normal(0, amplitude, 13)
                    if j % 5 == 0:
                        second = -first + rng.normal(0, .1, 13)
                    elif j % 5 == 1:
                        second[:3] = first[:3]
                        second[3:] = -first[3:]
                    initial = np.r_[first, second]
            row = optimize(e, initial, b, symmetric, args.maxiter)
            row['start'] = j
            data['coarse'].append(row)
            save()
        print(json.dumps({'stage': 'coarse', 'b': b,
                          'best': max(r['discovery_value'] for r in data['coarse'] if r['b']==b)}), flush=True)
    del e
    refine = Engine(args.refine_nodes)
    for b in args.b:
        rows = sorted((r for r in data['coarse'] if r['b']==b),
                      key=lambda r: r['discovery_value'], reverse=True)[:args.refine_top]
        # Always include the original one-dimensional reflected threshold.
        rows.append({'parameters': rotated_seed([1., 0., 0.]).tolist(), 'symmetric': True, 'start': -1})
        for row in rows:
            rr = optimize(refine, np.array(row['parameters']), b,
                          row['symmetric'], args.maxiter)
            rr['coarse_start'] = row['start']
            rr['checks'] = []
            data['refined'].append(rr)
            save()
        print(json.dumps({'stage': 'refined', 'b': b,
                          'best': max(r['discovery_value'] for r in data['refined'] if r['b']==b)}), flush=True)
    del refine
    for n in args.checks:
        check = Engine(n)
        for row in data['refined']:
            checked = check.evaluate(np.array(row['parameters']), row['b'], row['symmetric'], False)
            checked['nodes'] = n
            row['checks'].append(checked)
        save()
        print(json.dumps({'stage': 'check', 'nodes': n,
                          'best_by_b': {str(b): max(r['checks'][-1]['value'] for r in data['refined'] if r['b']==b)
                                        for b in args.b}}), flush=True)
        del check
    data['finished'] = True
    save()


if __name__ == '__main__':
    main()
