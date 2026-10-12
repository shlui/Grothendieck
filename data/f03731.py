#!/usr/bin/env python3
"""Independent directed MPFR replay of the complete ten-row L1 cover.

Reads exact rational coefficients and partitions, imports no producer code,
and checks newly certified scalar rows against independent whole-line audits.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
from functools import lru_cache
import hashlib
import json
import os
import sys

os.environ.setdefault('GR_MPFR_BITS', '384')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from mpfr_interval import I, selftest

SOURCE = ROOT/'tonight/operator/cubic_L1_refinement/refit847/complete_price_cover.json'
OLD_SOURCE = SOURCE.parent.parent/'kappa_shape_caps.json'
ENERGY = ROOT/'reduced_search/energy'
NAMES = ('0', '2', '3', '4', 'k8', '5', '6', '7', '68', '96')
BIN_INTERVALS = (('0','1'),('1','10'),('10','90'),('0','3/10'),
                 ('3/10','11/10'),('11/10','4'),('4','10'),
                 ('10','20'),('20','40'),('40','60'),('60','90'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upper_ceiling(interval):
    x = F(interval.hi.text(2))
    scale = 10**12
    return F(-((-x.numerator*scale)//x.denominator), scale)


def main():
    engine = selftest()
    initial_sha = sha(SOURCE)
    source = json.loads(SOURCE.read_text())
    assert source['status'] == 'PASS_COMPLETE_CUBIC_L1_REFIT_PRICE_COVER'
    assert F(source['target']) == F(847,1000)
    assert source['source_cover_sha256'] == sha(OLD_SOURCE)
    prior = json.loads(OLD_SOURCE.read_text())
    accepted_path = ENERGY/'cubic_L1_certificate.json'
    accepted = json.loads(accepted_path.read_text())
    assert accepted['status'] == 'PASS_DIRECTED_COMPLETE_CUBIC_L1_BOUND'
    old_hashes = {p['row_path']:p['row_sha256'] for p in accepted['rows']}
    assert len(old_hashes) == 8
    assert set(source['rows']) == set(NAMES)
    rows, row_audit = {}, {}
    for name in NAMES:
        record = source['rows'][name]
        path = (SOURCE.parent/record['path']).resolve()
        assert sha(path) == record['sha256']
        row = json.loads(path.read_text())
        assert row['status'] == 'PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        assert tuple(map(F,row['odd_weights'])) == (F(1),F(0),F(0))
        assert F(row['sixth_offset']) == 0
        q, eta = tuple(map(F,row['majorant'])), F(row['eta'])
        assert len(q) == 4 and q[3] > 0 and eta > 0
        assert tuple(map(F,record['majorant'])) == q
        assert F(record['eta']) == eta
        proof = 'previously accepted global row'
        independent_path = None
        if name in ('68','96'):
            independent_path = HERE/f'independent_row_{name}.json'
            independent = json.loads(independent_path.read_text())
            assert independent['status'] == 'PASS_FRESH_INDEPENDENT_WHOLE_LINE_CUBIC_L1_REFIT'
            assert independent['source_sha256'] == record['sha256']
            assert independent['compact']['domain'] == ['0','80']
            assert independent['tail']['domain'] == ['80','infinity']
            assert tuple(map(F,independent['majorant'])) == q
            assert F(independent['eta']) == eta
            proof = 'fresh independent 768-bit Taylor6 + exact Fraction tail'
        else:
            assert path.parent == ENERGY.resolve()
            assert old_hashes[path.name] == record['sha256']
            old = prior['rows'][name]
            assert tuple(map(F,old['majorant'])) == q
            assert F(old['eta']) == eta
            assert old['sha256'] == record['sha256']
        rows[name] = q,eta
        row_audit[name] = {
            'path':str(path.relative_to(ROOT)), 'sha256':record['sha256'],
            'majorant':list(map(str,q)), 'eta':str(eta), 'scalar_proof':proof,
            'independent_proof_sha256':sha(independent_path) if independent_path else None}

    c = (I(2)/I.pi()).sqrt()
    K = I(194)/I(10).sqrt()
    assert (K-I(F(633,15))).above(0)

    @lru_cache(None)
    def k32(k):
        return I(k)*I(k).sqrt()

    @lru_cache(None)
    def price(name,k,s):
        q0,q2,q4,q6 = rows[name][0]
        eta = rows[name][1]
        rational = (q0+q2+3*q4+15*q6
                    +k*(q4+15*q6+F(7,90)*eta+s*eta/F(30)))
        return c+I(rational)+I(q6)*(K-I(F(633,15)*s))*k32(k)

    @lru_cache(None)
    def rounded(name,k,s):
        return upper_ceiling(price(name,k,s))

    def corners(name,lo,hi):
        return tuple(rounded(name,k,s) for k in (lo,hi) for s in (F(0),F(1)))

    leaves = source['leaves']
    assert len(leaves) == len(prior['leaves']) == 2780
    previous, results = F(0), []
    for idx,(leaf,old_leaf) in enumerate(zip(leaves,prior['leaves'])):
        lo,hi = map(F,leaf['kappa_interval'])
        assert (lo,hi) == tuple(map(F,old_leaf['kappa_interval']))
        assert lo == previous and lo < hi <= 90
        previous = hi
        name = leaf['row']
        selected = corners(name,lo,hi)
        assert selected == tuple(map(F,leaf['corners']))
        assert max(selected) == F(leaf['cap'])
        checks = []
        for j,(k,s) in enumerate((k,s) for k in (lo,hi) for s in (F(0),F(1))):
            actual = price(name,k,s)
            assert actual.below(selected[j])
            checks.append({'kappa':str(k),'shape':str(s),
                           'price_interval':actual.record(),
                           'saved_ceiling':str(selected[j])})
        candidates = {n:max(corners(n,lo,hi)) for n in NAMES}
        optimal = min((cap,n) for n,cap in candidates.items())
        assert (F(leaf['cap']),name) == optimal
        results.append({'index':idx,'kappa_interval':[str(lo),str(hi)],
                        'row':name,'cap':str(optimal[0]),
                        'selected_corner_checks':checks,
                        'all_ten_candidate_caps':{n:str(v) for n,v in candidates.items()}})
    assert previous == 90

    table = []
    assert len(source['table']) == len(BIN_INTERVALS)
    for item,(a,b) in zip(source['table'],BIN_INTERVALS):
        lo,hi = F(a),F(b)
        assert tuple(map(F,item['kappa_interval'])) == (lo,hi)
        pieces = [(i,p) for i,p in enumerate(leaves)
                  if F(p['kappa_interval'][0]) >= lo
                  and F(p['kappa_interval'][1]) <= hi]
        assert pieces and len(pieces) == item['pieces']
        assert F(pieces[0][1]['kappa_interval'][0]) == lo
        assert F(pieces[-1][1]['kappa_interval'][1]) == hi
        maximum,active = max((F(p['cap']),i) for i,p in pieces)
        assert maximum == F(item['cap'])
        table.append({'kappa_interval':[str(lo),str(hi)],'cap':str(maximum),
                      'active_leaf':active,'pieces':len(pieces)})

    maximum,active = max((F(p['cap']),i) for i,p in enumerate(leaves))
    assert maximum == F(source['maximum_price'])
    assert maximum < F(847,1000)
    assert sha(SOURCE) == initial_sha
    points = []
    for k in map(F,('0','53/100','4','271/40','383/40','10','90')):
        cap,name = min((max(corners(n,k,k)),n) for n in NAMES)
        points.append({'kappa':str(k),'cap':str(cap),'row':name,
                       'shape0':price(name,k,F(0)).record(),
                       'shape1':price(name,k,F(1)).record()})
    output = {
        'status':'PASS_INDEPENDENT_COMPLETE_CUBIC_L1_REFIT847_PRICE_COVER',
        'engine':engine,'source_path':str(SOURCE.relative_to(ROOT)),
        'source_sha256':initial_sha,'old_source_sha256':sha(OLD_SOURCE),
        'accepted_energy_certificate_sha256':sha(accepted_path),
        'rows':row_audit,'leaf_count':len(leaves),
        'selected_corner_count':4*len(leaves),
        'all_ten_candidate_corner_checks':10*4*len(leaves),
        'unique_mpfr_price_evaluations':price.cache_info().currsize,
        'row_leaf_counts':dict(Counter(p['row'] for p in leaves)),
        'coverage':{'interval':['0','90'],'exact_rational':True,'no_gaps':True},
        'maximum_price':str(maximum),'active_leaf':active,
        'active_kappa_interval':leaves[active]['kappa_interval'],
        'target':'847/1000','squared_energy_bound':'717409/1000000',
        'analytic_validity': {
            'formula':'c+q0+q2+3q4+15q6+k(q4+15q6+7eta/90+s eta/30)+q6(K-633s/15)k^(3/2)',
            'convexity':'q6>0 and K-633/15>0 imply convexity in kappa, including 0 by continuity.',
            'shape':'Affine in shape; endpoint checks cover all shape0..1.',
            'new_scalar_rows':'Fresh complete half-line proofs with matched SHA; old pricing intervals do not restrict global scalar inequalities.',
            'premises':'Standard normalized third-chaos Stein identity, sharp633 coupled moments, and 0<=kappa<=90.'},
        'table':table,'points':points,'leaves':results,
        'scope':'Independent global cubic L1 cap .847, with two freshly independently certified scalar rows. This is not a global Grothendieck operator norm cap.'}
    (HERE/'independent_complete_price_cover.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in (
        'status','engine','leaf_count','selected_corner_count',
        'all_ten_candidate_corner_checks','unique_mpfr_price_evaluations',
        'row_leaf_counts','maximum_price','active_leaf','active_kappa_interval',
        'target','table','points')},indent=2),flush=True)


if __name__ == '__main__':
    main()
