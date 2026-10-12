#!/usr/bin/env python3
"""Producer-free directed sixteen-vertex prices for accepted rank q gates."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse, hashlib, json, sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE/'global_b7171'))
import replay_R3_corners as R

def main(certificate,acceptance,output):
    d=R.load(certificate); gate=R.load(acceptance)
    assert gate['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    path=ROOT/gate['profile_path'];raw=R.load(path)
    assert R.sha(path)==gate['profile_sha256']==d['profile_binding']['sha256']
    row=dict(kind='rank',p=R.exact(raw['profile']),profile=path,
             reserve0=F(gate['reserves'][0]),reserve6=F(gate['reserves'][1]))
    assert (row['reserve0'],row['reserve6'])==(F(d['reserve0']),F(d['reserve6']))
    for item in gate['reports'].values():
        assert R.sha(ROOT/item['path'])==item['sha256']
    bounds={k:tuple(map(F,v))for k,v in d['source_rectangle'].items()}
    assert set(bounds)=={'u','r','kappa','shape'}
    assert all(len(v)==2 and v[0]<=v[1] for v in bounds.values())
    results=[]
    for old in d['results']:
        b=F(old['operator_b']);c=R.coefficients(row,b);spectral=R.spectral(row,b)
        assert c['A']>0 and c['B']>0
        saved={tuple(F(x[k])for k in ('u','r','kappa','shape')):F(x['upper']) for x in old['rows']}
        assert len(saved)==16
        corners=[]
        for u,r,k,z in product(*(bounds[x]for x in ('u','r','kappa','shape'))):
            actual=R.price(row,c,k,z,u,r)
            assert actual.below(saved[u,r,k,z]),('SAVED_DIRECTED_UPPER',b,u,r,k,z,actual.record())
            assert actual.below(F(old['clean_rational_upper']))
            corners.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),price=actual.record()))
        results.append(dict(operator_b=str(b),clean_rational_upper=old['clean_rational_upper'],
                            spectral=spectral,corners=corners))
    out=dict(status='PASS_INDEPENDENT_COMPLETE_RANK_REGIONAL_SOURCE_BOX',
             source=R.file_record(certificate),accepted_whole_q_gate=R.file_record(acceptance),
             profile=R.file_record(path),source_rectangle=d['source_rectangle'],
             engine=R.selftest(),results=results,directed_prices=16*len(results),
             convexity='Amplitude and radial quadratic coefficients are positive. Expanded cumulant/shape coefficients A and B are positive, so each price is separately convex in kappa and shape. All sixteen vertices enclose the whole box; target coefficient spectral conditions are checked exactly.',
             scope='Only these source boxes, at the displayed common operator coefficient. Global coverage requires a separate exhaustive ledger and exact assembly.')
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ('status','source_rectangle','directed_prices')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('certificate',type=Path)
    a.add_argument('acceptance',type=Path);a.add_argument('--output',type=Path,required=True)
    x=a.parse_args();main(x.certificate.resolve(),x.acceptance.resolve(),x.output.resolve())
