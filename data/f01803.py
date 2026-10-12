#!/usr/bin/env python3
"""Compose complete SAME-b global acceptance with a producer-free R3 delta."""
from pathlib import Path
from fractions import Fraction as F
import argparse, gzip, hashlib, json, os
os.environ['GR_MPFR_BITS']='768'
import replay_R3_corners as R

ROOT=Path(__file__).resolve().parents[3]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def rec(p): return {'path':str(p.resolve().relative_to(ROOT)), 'sha256':sha(p)}

def main(x):
    old=read(x.parent_global); d=read(x.envelope); g=read(x.geometry)
    assert old['status'] in ('PASS_INDEPENDENT_COMPLETE_NATIVE_KAPPA_REFINED_GLOBAL_SCALAR_CERTIFICATE',
                             'PASS_INDEPENDENT_COMPLETE_MONOTONE_REFINED_GLOBAL_SCALAR_CERTIFICATE')
    assert g['status']=='PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE'
    digest=sha(x.envelope); b=F(d['operator_b'])
    assert b==F(old['operator_b'])==F(g['operator_b'])
    assert digest==g['source_sha256']
    parent=g['parent']; pp=ROOT/parent['path']; pg=ROOT/parent['accepted_geometry_audit']
    assert sha(pp)==parent['sha256']==old['envelope']['sha256']
    assert sha(pg)==parent['audit_sha256']
    assert parent['audit_sha256'] in {v['sha256'] for v in old['accepted_geometry_proof_chain']}
    assert g['identical_priced_cells']==old['full_upper_priced_cells']==len(d['cells'])
    assert g['identical_empty_cells']==old['proved_physically_empty_full_upper_cells']==len(d['physically_empty_cells'])
    assert g['amplitude_interval']==d['amplitude_interval']==read(pp)['amplitude_interval']
    numeric=g['numeric_cache_audit']
    assert numeric['status']=='PASS_FRESH_PRODUCER_FREE_NEW_MPFR768_CORNERS'
    assert g['improved_cells']==sum(g['new_winner_counts'].values())
    whole=[]; fresh_corners=0
    old_source=read(pp)
    names=set(g['new_profiles'])
    assert [v for v in d['cached_corner_constants'] if v['profile'] not in names]==old_source['cached_corner_constants']
    for name, bind in g['new_profiles'].items():
        profile=ROOT/bind['path']; gate=ROOT/bind['report_path']; a=read(gate)
        assert sha(profile)==bind['sha256']==a['profile_sha256']
        assert sha(gate)==bind['report_sha256'] and a['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        assert list(map(F,a['reserves']))==list(map(F,bind['reserves']))
        row={'kind':'rank','p':R.exact(read(profile)['profile']),
             'reserve0':F(bind['reserves'][0]),'reserve6':F(bind['reserves'][1])}
        coefficients=R.coefficients(row,b); R.spectral(row,b)
        assert coefficients['A']>0 and coefficients['B']>0
        for key in ('r_quadratic','r_linear','ell'):
            assert F(d['profile_coefficients'][name][key])==coefficients[key]
        for table in d['cached_corner_constants']:
            if table['profile']!=name: continue
            ka,kz=map(F,table['kappa_interval']); assert 0<=ka<kz<=90
            for k in (ka,kz):
                for shape in (F(2,5),F(1)):
                    assert R.price(row,coefficients,k,shape).below(F(table['upper']))
                    fresh_corners+=1
        for part, r in a['reports'].items():
            rp=ROOT/r['path']; assert sha(rp)==r['sha256'] and read(rp)['status'].startswith('PASS_INDEPENDENT_')
        whole.append({**rec(gate),'profile':name})
    assert whole and numeric['fresh_corners']==fresh_corners>=4*len(whole)
    cap=F(d['strict_rational_cap']); maximum=F(d['maximum_cell_cap'])
    assert maximum==F(g['maximum_cell_cap'])<cap==F(g['strict_rational_cap'])<=F(old['norm_upper_bound'])
    assert F(old['lower_amplitude_upper_bound'])<cap
    out={**old,
         'status':'PASS_INDEPENDENT_COMPLETE_MONOTONE_REFINED_GLOBAL_SCALAR_CERTIFICATE',
         'norm_upper_bound':str(cap),'constructive_G_R_lower_bound':str(1/cap),
         'constructive_G_R_lower_bound_display':float(1/cap),'envelope':rec(x.envelope),
         'full_upper_amplitude_upper_bound':str(maximum),
         'complete_whole_q_profiles':len(d['profile_coefficients']),
         'monotone_bindings':{'parent_global':rec(x.parent_global),'geometry':rec(x.geometry),
                              'new_complete_whole_q_gates':whole},
         'independent_new_corner_prices':fresh_corners,
         'independent_new_polygon_prices':g['improved_polygon_prices'],
         'accepted_geometry_proof_chain':[rec(x.geometry)]+old['accepted_geometry_proof_chain'],
         'analytic_review':'The complete same-coefficient parent proof is retained on the identical exhaustive physical partition. Every changed price is independently recomputed from newly accepted complete whole-line R3 gates, with fresh directed moment corners and exact polygon edge/knot extrema; unchanged parent prices remain valid. The complete lower-amplitude parent bound remains below the final cap. This finite delta uses no inherited historical Grothendieck endpoint premise.'}
    x.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','operator_b','norm_upper_bound','constructive_G_R_lower_bound_display','complete_whole_q_profiles')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser()
    for key in ('parent-global','envelope','geometry','output'):
        a.add_argument('--'+key,required=True,type=Path)
    main(a.parse_args())
