#!/usr/bin/env python3
"""Combine a complete same-b global proof and an exact native-only refinement."""
from pathlib import Path
from fractions import Fraction as F
import argparse, gzip, hashlib, json

ROOT=Path(__file__).resolve().parents[3]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def rec(p): return {'path':str(p.resolve().relative_to(ROOT)), 'sha256':sha(p)}

def main(x):
    old=read(x.parent_global); d=read(x.envelope); g=read(x.geometry)
    n=read(x.native); p=read(x.plane); baseline=read(x.baseline)
    assert old['status']=='PASS_INDEPENDENT_COMPLETE_FULLUPPER_GLOBAL_SCALAR_CERTIFICATE'
    assert g['status']=='PASS_PRODUCER_FREE_NATIVE_CUMULANT_REFINEMENT'
    assert all(v['status'].startswith('PASS') for v in (n,p,baseline))
    digest=sha(x.envelope); b=F(d['operator_b']); parent=d['native_kappa_refinement']
    assert b==F(old['operator_b'])==F(g['operator_b'])==F(n['operator_b'])==F(p['operator_b'])
    assert digest==g['source_sha256']==n['input_sha256']==p['input_sha256']
    assert old['envelope']['sha256']==parent['parent_sha256']==g['parent']['sha256']
    assert sha(ROOT/parent['parent_path'])==parent['parent_sha256']
    gp=ROOT/g['parent']['accepted_audit']
    assert sha(gp)==g['parent']['audit_sha256']
    assert g['parent']['audit_sha256'] in {v['sha256'] for v in old['accepted_geometry_proof_chain']}
    assert baseline['status']=='PASS_INDEPENDENT_ALL_FINAL68_KAPPA_NATIVE_SELECTED_ENDPOINT_PRICES'
    assert baseline['precision_bits']>=768 and baseline['selected_endpoint_prices_recomputed']==13936
    assert n['baseline_numeric_audit']['sha256']==sha(x.baseline)
    assert g['new_native_numeric_audit']['sha256']==sha(x.baseline)
    assert n['accepted_baseline_directed_prices']==13936 and n['kappa_bands']==67
    assert n['target_spectral_gates_checked'] and len(n['frames'])==7
    assert p['new_kappa_plane_bands']==g['new_plane_bands']==14
    assert p['accepted_parent_source_sha256']==parent['parent_sha256']
    assert p['fresh_directed_corner_prices']>=g['new_plane_directed_corner_checks']==560
    assert g['refined_parents']+g['unchanged_price_cells']==old['full_upper_priced_cells']
    assert g['priced_children']+g['exact_empty_children']==g['exact_child_intervals']
    assert g['total_priced_cells']==g['unchanged_price_cells']+g['priced_children']==len(d['cells'])
    assert g['total_empty_cells']==old['proved_physically_empty_full_upper_cells']+g['exact_empty_children']==len(d['physically_empty_cells'])
    assert d['amplitude_interval']==g['amplitude_interval']
    assert list(map(F,d['amplitude_interval']))==[F(1381,2000),F(159577,200000)]
    cap=F(d['strict_rational_cap']); maximum=F(d['maximum_cell_cap'])
    assert maximum==F(g['maximum_cell_cap'])<cap==F(g['strict_rational_cap'])<F(old['norm_upper_bound'])
    assert F(old['lower_amplitude_upper_bound'])<cap
    out={**old,
         'status':'PASS_INDEPENDENT_COMPLETE_NATIVE_KAPPA_REFINED_GLOBAL_SCALAR_CERTIFICATE',
         'norm_upper_bound':str(cap),'constructive_G_R_lower_bound':str(1/cap),
         'constructive_G_R_lower_bound_display':float(1/cap),'envelope':rec(x.envelope),
         'full_upper_amplitude_upper_bound':str(maximum),'cumulant_points':68,
         'full_upper_original_cells':g['total_priced_cells']+g['total_empty_cells'],
         'full_upper_priced_cells':g['total_priced_cells'],
         'proved_physically_empty_full_upper_cells':g['total_empty_cells'],
         'new_native_endpoint_prices':13936,'new_uniform_source_planes':140,
         'refinement_bindings':{k:rec(getattr(x,k)) for k in ('parent_global','geometry','native','plane','baseline')},
         'accepted_geometry_proof_chain':[rec(x.geometry)]+old['accepted_geometry_proof_chain'],
         'analytic_review':'The complete SAME-b parent norm proof is retained. Exactly selected cumulant cells are partitioned into children, all inherited bounds remain valid, and each new native bound uses fresh directed endpoints and valid uniform analytic source planes. The independent exact geometry audit checks every child, source polygon, radial knot and edge stationary point. The complete lower-amplitude branch is unchanged and lies below the new cap. No new R3 integration is required or assumed.'}
    x.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','operator_b','norm_upper_bound','constructive_G_R_lower_bound_display','full_upper_original_cells')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser()
    for key in ('parent-global','envelope','geometry','native','plane','baseline','output'):
        a.add_argument('--'+key,required=True,type=Path)
    main(a.parse_args())
