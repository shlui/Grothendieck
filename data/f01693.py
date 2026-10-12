#!/usr/bin/env python3
"""Exact no-optimizer assembly of the freshly checked .30--.31 region."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
from mpfr_interval import I, selftest

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SOURCE = BASE/'tonight/recovery/profile_cover/continuous_e7200_u3000_3100_r6_d10.json'
LOCAL = BASE/'tonight/recovery/profile_cover/local_physical_d256_inflation1_125_summary.json'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--local-summary',type=Path,default=LOCAL)
    parser.add_argument('--scalar-dir',type=Path,default=HERE)
    args=parser.parse_args()
    data = json.loads(args.source.read_text())
    local = json.loads(args.local_summary.read_text())
    assert local['status'] == 'PASS' and len(local['nodes']) == 51
    assert local['source_sha256'] == hashlib.sha256(args.source.read_bytes()).hexdigest()
    assert sorted(n['node_index'] for n in local['nodes']) == list(range(51))
    assert all(n['status'] == 'PASS' and n['pending_boxes'] == 0 for n in local['nodes'])
    p = local['profile']
    m = [[F(x) for x in row] for row in p['M']]
    original = data['profile']
    for i in range(4):
        for j in range(4):
            assert m[i][j] == F(str(original['M'][i][j]))+(F(1,1000) if i == j else 0)
    ell,beta,delta,nu,b = [F(p[k]) for k in ('ell','beta','delta','nu','epsilon')]
    assert b == F(18,25)
    for k,v in zip(('ell','beta','delta','nu','epsilon'), (ell,beta,delta,nu,b)):
        assert v == F(str(original[k]))
    for new,old in zip(p['nodal_majorant'], original['nodal_majorant']):
        assert F(new) == F(str(old))+F(1,125)
    assert F(p['tail_slope']) == F(str(original['tail_slope']))+F(1,10**6)
    assert F(p['tail_slope'])**2 >= 4*(delta*delta+nu*nu)
    gates = []
    for name,a,off,c in (('E',m[0][0],m[0][2],m[2][2]),
                        ('R',m[1][1],m[1][3],m[3][3]),
                        ('E-minus-D13',m[0][0]+1,m[0][2],m[2][2]-b),
                        ('R-minus-D13',m[1][1]+1,m[1][3],m[3][3]-b),
                        ('R-minus-D2',m[1][1]-b,m[1][3],m[3][3])):
        determinant = a*c-off*off
        assert a > 0 and c > 0 and determinant > 0
        gates.append({'name':name,'first_diagonal':str(a),'second_diagonal':str(c),
                      'off_diagonal':str(off),'determinant':str(determinant)})
    gram = (nu-m[0][2])**2/m[2][2]-m[0][0]+b-2*delta
    coarse_curvature=F(str(original['gram_negative_bound']))
    assert gram >= -coarse_curvature
    grid = list(map(lambda x:F(str(x)), data['r_grid']))
    assert grid[0] == 0 and all(a < c for a,c in zip(grid[:-1],grid[1:]))
    rho = F(849,1000)
    lo,hi = map(lambda x:F(str(x)), original['amplitude_interval'])
    allowed={(F(3,10),F(31,100)):F(541105,10**6),
             (F(21,100),F(3,10)):F(546399,10**6)}
    assert (lo,hi) in allowed
    c = (I(2)/I.pi()).sqrt()
    z = (I(lo)/c-(1-rho*rho))/(rho*rho)
    if z.above(0):
        assert z.below(1)
        source_cap = I(rho)*(-z.log()).sqrt().erf()
        assert source_cap.below(grid[-1])
    else:
        assert I(lo).below(c*(1-rho*rho))
        source_cap=I(rho)
        assert grid[-1]>=rho
    max_width_squared=max((c-a)**2 for a,c in zip(grid[:-1],grid[1:]))
    coarse_reserve = coarse_curvature*max_width_squared/4
    reserve = max(-gram,F(0))*max_width_squared/4
    tail_charge = F(1,10**6)*grid[-1]*rho
    upper = F(str(original['nodal_majorant'][0]))
    ledger = []
    active = {'kind':'zero_source'}
    for index in range(len(data['rows'])):
        path = args.scalar_dir/f'continuous_scalar_row_{index:02d}.json'
        row = json.loads(path.read_text())
        assert row['status'] == 'PASS_FRESH_COMPLETE_SCALAR_HALF_LINE'
        assert row['source_sha256'] == hashlib.sha256(args.source.read_bytes()).hexdigest()
        raw = data['rows'][index]
        r = F(row['frozen_row']['r'])
        assert r == F(str(raw['r'])) == grid[index//3+1]
        assert list(map(F,row['frozen_row']['kappa_interval'])) == list(map(F,raw['kappa_interval']))
        assert F(row['frozen_row']['added_q0']) == F(1,500)
        assert F(row['frozen_row']['added_q6']) == F(1,20_000_000)
        entries = []
        for corner in row['prices']:
            value = F(corner['price_upper'])+gram*r*r
            entries.append({'kappa':corner['kappa'],'shape':corner['shape'],
                            'scalar_price_upper':corner['price_upper'],
                            'with_gram_upper':str(value)})
            if value > upper:
                upper = value
                active = {'row':index,'r':str(r),'kappa':corner['kappa'],'shape':corner['shape']}
        ledger.append({'row':index,'r':str(r),
                       'kappa_interval':row['frozen_row']['kappa_interval'],
                       'corners':entries,'file_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    amplitude = [{'u':str(u),'correction':str(u*u-2*ell*u)} for u in (lo,hi)]
    amplitude_max = max(F(x['correction']) for x in amplitude)
    cap = upper+reserve+F(1,125)+tail_charge+amplitude_max
    coarse_cap = upper+coarse_reserve+F(1,125)+tail_charge+amplitude_max
    simple_cap = allowed[lo,hi]
    assert cap < simple_cap < F(40,71)
    coarse_simple_cap=F(556109,10**6) if lo==F(3,10) else F(547351,10**6)
    assert coarse_cap < coarse_simple_cap
    result = {'status':'PASS_EXACT_REGIONAL_ASSEMBLY', 'engine':selftest(),
              'interval':[str(lo),str(hi)],'operator':'P1-(18/25)P3',
              'profile':p,'PSD_gates':gates,'gram_cost':str(gram),
              'r_grid':list(map(str,grid)), 'source_cap_interval':source_cap.record(),
              'zero_source_upper':str(F(str(original['nodal_majorant'][0]))),
              'endpoint_ledger':ledger,'maximum_endpoint_with_gram':str(upper),
              'active_corner':active,'quadratic_chord_reserve':str(reserve),
              'chord_curvature_bound':str(max(-gram,F(0))),
              'coarse_curvature_bound':str(coarse_curvature),
              'coarse_reserve':str(coarse_reserve),
              'coarse_assembled_cap':str(coarse_cap),
              'coarse_strict_rational_cap':str(coarse_simple_cap),
              'uniform_height_inflation':'1/125','tail_slope_charge':str(tail_charge),
              'amplitude_endpoints':amplitude,'amplitude_correction':str(amplitude_max),
              'exact_assembled_cap':str(cap),'cap_display':float(cap),
              'strict_rational_cap':str(simple_cap),
              'local_cover_nodes':51,'local_boxes':local['total_boxes'],
              'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
              'local_summary_sha256':hashlib.sha256(args.local_summary.read_bytes()).hexdigest(),
              'scope':'Regional bilinear operator cap only; no global Grothendieck endpoint.'}
    filename=f'regional_u{int(lo*10000):04d}_{int(hi*10000):04d}_endpoint_ledger.json'
    (HERE/filename).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'cap_display':result['cap_display'],
                      'rational_cap':str(simple_cap),'active_corner':active},indent=2))


if __name__ == '__main__':
    main()
