#!/usr/bin/env python3
"""Independent numeric premises of a source-polygon envelope.

Uses the accepted all-line contact coefficients, not the producer's pricing
functions. A single fixed polynomial is priced at the four cumulant/shape
corners before any half-plane is used. All inverse caps are forward checked.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, sys, gzip
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mpfr_interval import I, BITS
ROOT = Path(__file__).resolve().parents[3]
ap = argparse.ArgumentParser()
ap.add_argument('ledger', type=Path)
ap.add_argument('--output', type=Path, required=True)
args = ap.parse_args()
src = args.ledger.resolve()
d = json.loads(gzip.decompress(src.read_bytes()) if src.suffix=='.gz' else src.read_text())
assert d['source_polygon_pricing'] and 0<F(d['operator_b'])<=F(7171,10000)
bind = d['linear_source_bank_binding']
bankpath = ROOT / bind['path']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(bankpath) == bind['sha256']
bank = json.loads(bankpath.read_text())
accepted = json.loads((ROOT/'tonight/audit/linear_contact_bank_independent.json').read_text())
assert accepted['status'] == 'PASS_INDEPENDENT_COMPLETE_LINEAR_CONTACT_SOURCE_BANK'
assert accepted['input_sha256'] == sha(bankpath)
K = I(194)/I(10).sqrt()
assert (K-F(633,15)).above(0)
cache = {}
checks = 0
def check(row_index, klo, khi, cap, a):
    global checks
    row = bank['rows'][row_index]
    assert F(row['a']) == a > 0
    qs = list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']))
    assert qs[3] > 0 and 0 <= klo <= khi <= 90
    key=(row_index,klo,khi)
    if key not in cache:
        values=[]
        for k in (klo,khi):
            for s in (F(0),F(1)):
                v = I(qs[0]+qs[1]+(3+k)*qs[2]) + I(qs[3])*(15+15*k+(K-F(633,15)*s)*I(k)*I(k).sqrt())
                values.append(v)
                checks += 1
        cache[key]=values
    assert all(v.below(cap) for v in cache[key]), (row_index,klo,khi,cap)
bands = []
plane_count=0
for group in d['source_polygon_planes']:
    klo,khi=map(F,group['kappa_interval'])
    bands.append((klo,khi))
    heads=set()
    for p in group['planes']:
        a,cap=F(p['a']),F(p['cap'])
        assert a not in heads
        heads.add(a)
        check(p['source_row'],klo,khi,cap,a)
        plane_count+=1
    assert len(heads)==10
assert bands[0][0]==0 and bands[-1][1]==90
assert all(left[1]==right[0] for left,right in zip(bands,bands[1:]))
cell_plane_rows=0
for cell in d['cells']+d['physically_empty_cells']:
    klo,khi=map(F,cell['kappa_interval'])
    if 'source_planes' in cell:
        p=cell['source_planes']
        for side in ('u','r'):
            check(p[side+'_row'],klo,khi,F(p[side+'_price_upper']),F(p[side+'_a']))
            cell_plane_rows+=1
nu=(I(2)/I.pi()).sqrt()
rho=F(847,1000)
inverse_rows=[]
for bound in d['physical_u_bounds']:
    r,u=F(bound['r_lower']),F(bound['upper'])
    assert 0<=r<rho
    if nu.below(u): method='u<=sqrt(2/pi)'
    else:
        assert r>0
        z=(I(u)/nu-1+rho*rho)/(rho*rho)
        assert z.above(0) and z.below(1)
        assert (I(rho)*(-z.log()).sqrt().erf()).below(r)
        method='forward R_rho(u_cap)<=r_lower'
    inverse_rows.append(dict(r_lower=str(r),u_upper=str(u),method=method))
out=dict(status='PASS_INDEPENDENT_UNIFORM_POLYGON_PLANES_AND_INVERSE_CAPS',
         input_sha256=sha(src),operator_b=d['operator_b'],precision_bits=BITS,bank_sha256=sha(bankpath),
         kappa_bands=len(bands),uniform_polygon_planes=plane_count,
         recorded_cell_plane_rows=cell_plane_rows,distinct_row_bands=len(cache),
         fresh_directed_corner_prices=checks,inverse_caps=len(inverse_rows),
         producer_code_imported=False,
         analytic_review='Each fixed q6-positive all-line row has a price convex in kappa and affine in shape; four endpoint prices uniformly bound E J_A(Q) over the whole band. Thus a*u+r<=cap is a valid physical plane before clipping. Inverse caps are independently checked via the monotone forward source function.',
         scope='Numeric premises only. Complete polygon/source coverage, all selected native/R3 prices, and final global union require separate exact audits.')
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
