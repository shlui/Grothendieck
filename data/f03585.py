#!/usr/bin/env python3
"""Combine independently accepted full-upper and same-b lower proofs."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,os,sys
os.environ['GR_MPFR_BITS']='768'
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
from mpfr_interval import I

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def record(p):return {'path':str(p.resolve().relative_to(ROOT)),'sha256':sha(p)}

def complete_geometry(path,chain):
    """Expand accepted monotone/refinement proof chains without losing scope."""
    raw=read(path);assert raw['status'].startswith('PASS_PRODUCER_FREE_')
    chain.append(record(path))
    if raw['status']=='PASS_PRODUCER_FREE_EXACT_POLYGON_GLOBAL_MIDDLE_COVER':return raw
    assert raw['status']in ('PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE','PASS_PRODUCER_FREE_ADAPTIVE_POLYGON_REFINEMENT')
    p=raw['parent'];audit=ROOT/(p.get('accepted_geometry_audit')or p['accepted_audit'])
    assert sha(audit)==p['audit_sha256'] and sha(ROOT/p['path'])==p['sha256']
    old=complete_geometry(audit,chain)
    assert old['source_sha256']==p['sha256'] and old['operator_b']==raw['operator_b']
    assert old['amplitude_interval']==raw['amplitude_interval']
    out={**old,**{k:raw[k]for k in ('source_sha256','operator_b','amplitude_interval','maximum_cell_cap','strict_rational_cap')}}
    if raw['status']=='PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE':
        assert raw['identical_priced_cells']==old['accepted_price_cells'] and raw['identical_empty_cells']==old['exact_physically_empty_cells']
    else:
        assert raw['refined_parents']+raw['unchanged_price_cells']==old['accepted_price_cells']
        assert raw['priced_children']+raw['exact_empty_children']==raw['exact_child_intervals']
        out['accepted_price_cells']=raw['total_priced_cells'];out['exact_physically_empty_cells']=raw['total_empty_cells']
        assert out['accepted_price_cells']==raw['unchanged_price_cells']+raw['priced_children']
        assert out['exact_physically_empty_cells']==old['exact_physically_empty_cells']+raw['exact_empty_children']
        out['original_cells']=out['accepted_price_cells']+out['exact_physically_empty_cells']
    return out

def main(x):
    source=read(x.envelope);digest=sha(x.envelope);b=F(source['operator_b'])
    assert 0<b<=F(7171,10000) and source['status'].startswith('PASS')
    cap=x.cap or F(source['strict_rational_cap']);assert 0<cap<1
    files={k:getattr(x,k)for k in ('native_baseline','native','cache','plane','contact','geometry','lower')}
    proofs={k:read(p)for k,p in files.items()};assert all(v['status'].startswith('PASS')for v in proofs.values())
    baseline,native,cache,plane,contact,geometry,lower=[proofs[k]for k in files]
    geometry_chain=[];geometry=complete_geometry(x.geometry,geometry_chain)
    assert baseline['status']=='PASS_INDEPENDENT_ALL_FINAL57_KAPPA_NATIVE_SELECTED_ENDPOINT_PRICES'
    assert baseline['precision_bits']>=768 and baseline['selected_endpoint_prices_recomputed']==11648
    assert native['source']['sha256']==cache['source']['sha256']==digest
    assert plane['input_sha256']==contact['input_sha256']==geometry['source_sha256']==digest
    assert all(F(v['operator_b'])==b for v in (native,cache,plane,geometry,lower))
    assert native['baseline_numeric_audit']['sha256']==sha(x.native_baseline)
    assert native['accepted_baseline_directed_prices']==baseline['selected_endpoint_prices_recomputed']
    assert native['target_spectral_gates_checked'] is True
    nb=len(source['kappa_grid'])-1;np=len(source['profile_coefficients'])
    assert nb==native['kappa_bands']==cache['kappa_bands']==plane['kappa_bands']==geometry['kappa_bands']
    assert cache['profiles']==np and cache['corner_tables']==np*nb and cache['directed_corner_prices']==4*np*nb
    assert plane['uniform_polygon_planes']==10*nb and plane['fresh_directed_corner_prices']==40*nb
    assert plane['inverse_caps']==geometry['directed_inverse_caps']==len(source['physical_u_bounds'])
    assert contact['selected_contact_cells']+contact.get('inactive_contact_metadata_cells',0)==sum('contact_credit'in v for v in source['cells'])
    if contact.get('inactive_contact_metadata_cells',0):
        premise=contact['accepted_parent_contact_premise'];pp=ROOT/premise['path']
        assert sha(pp)==premise['sha256'] and read(pp)['input_sha256']==premise['input_sha256']
        ancestor_sources={read(ROOT/v['path'])['source_sha256']for v in geometry_chain}
        assert premise['input_sha256']in ancestor_sources
    assert geometry['accepted_price_cells']==len(source['cells'])
    assert geometry['exact_physically_empty_cells']==len(source['physically_empty_cells'])
    assert geometry['original_cells']==len(source['cells'])+len(source['physically_empty_cells'])
    assert geometry['amplitude_interval']==source['amplitude_interval']
    a,z=map(F,source['amplitude_interval']);assert a==F(1381,2000)
    assert (I(2)/I.pi()).sqrt().below(z)
    assert lower['status']=='PASS_INDEPENDENT_COMPLETE_SMALLER_B_LOWER_TARGET_BRANCH'
    assert list(map(F,lower['complete_amplitude_coverage']))==[F(0),a]
    assert lower['source_rho']=='847/1000' and lower['native_cells']==39929 and lower['R3_cells']==41
    upper=F(geometry['maximum_cell_cap']);low=F(lower['clean_rational_upper'])
    assert upper==F(source['maximum_cell_cap'])<cap and low<cap
    assert F(source['strict_rational_cap'])>upper
    assert sha(x.envelope)==digest
    out={'status':'PASS_INDEPENDENT_COMPLETE_FULLUPPER_GLOBAL_SCALAR_CERTIFICATE',
         'operator_b':str(b),'norm_upper_bound':str(cap),
         'constructive_G_R_lower_bound':str(1/cap),'constructive_G_R_lower_bound_display':float(1/cap),
         'envelope':record(x.envelope),'full_amplitude_coverage':['0','sqrt(2/pi)'],
         'lower_amplitude_upper_bound':str(low),'full_upper_amplitude_upper_bound':str(upper),
         'source_rho':'847/1000','cumulant_points':nb+1,'complete_whole_q_profiles':np,
         'full_upper_original_cells':geometry['original_cells'],
         'full_upper_priced_cells':geometry['accepted_price_cells'],
         'proved_physically_empty_full_upper_cells':geometry['exact_physically_empty_cells'],
         'full_upper_complete_groups':geometry['complete_groups'],
         'baseline_native_endpoint_prices':baseline['selected_endpoint_prices_recomputed'],
         'independent_R3_corner_prices':cache['directed_corner_prices'],
         'independent_plane_corner_prices':plane['fresh_directed_corner_prices'],
         'uniform_physical_source_planes':plane['uniform_polygon_planes'],
         'forward_verified_inverse_source_caps':plane['inverse_caps'],
         'fixed_contact_credit_cells':contact['selected_contact_cells'],
         'bindings':{k:{**record(p),'status':proofs[k]['status']}for k,p in files.items()},
         'accepted_geometry_proof_chain':geometry_chain,
         'historical_1773_or_1780498_premises_used':False,
         'analytic_review':'The low-amplitude proof and exhaustive full-upper proof use the SAME operator coefficient. Their intervals meet at .6905 and jointly contain the entire physical amplitude range [0,sqrt(2/pi)]. Each proof covers all physical source, cumulant and tensor-shape values, including exact positive-head exclusions. Every numeric premise and complete geometric assembly is separately accepted and hash-bound. Absolute norm follows by sign symmetry; odd-sign reduction was proved in the construction theorem.',
         'matrix_sequence_limit':'Actual limit is 1/N(b), at least the certified reciprocal cap; equality to the cap and convergence to unrestricted G_R are not asserted for this fixed family.',
         'scope':'Complete fresh all-dimension scalar norm and explicit matrix-sequence lower bound. This is not a three-significant-digit enclosure of G_R.'}
    x.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('bindings','analytic_review')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser()
    for name in ('envelope','native','cache','plane','contact','geometry','lower','output'):
        a.add_argument('--'+name,required=True,type=Path)
    a.add_argument('--native-baseline',type=Path,default=Path(__file__).with_name('final57_kappa_native_endpoint_replay.json'))
    a.add_argument('--cap',type=F);main(a.parse_args())
