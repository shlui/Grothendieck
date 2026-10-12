#!/usr/bin/env python3
"""Compose an accepted all-source baseline with a uniform b=.711 transport."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[3]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(p):return {'path':str(p.resolve().relative_to(ROOT)),'sha256':sha(p)}

def main(x):
    old=read(x.parent_global);g=read(x.transport);s=read(x.spectral);low=read(x.lower)
    assert old['status']=='PASS_INDEPENDENT_COMPLETE_MONOTONE_REFINED_GLOBAL_SCALAR_CERTIFICATE'
    assert F(old['operator_b'])==F(713,1000)
    assert g['status']=='PASS_PRODUCER_FREE_COMPLETE_ASSIGNED_B711_CELL_TRANSPORT'
    assert s['status']=='PASS_INDEPENDENT_ALL_ASSIGNED_PROFILE_BETA_SPECTRAL_BINDINGS'
    assert low['status']=='PASS_INDEPENDENT_COMPLETE_SMALLER_B_LOWER_TARGET_BRANCH'
    b=F(711,1000);assert all(F(v['operator_b'])==b for v in (g,s,low))
    assert F(g['baseline_b'])==F(s['baseline_b'])==F(old['operator_b'])
    assert s['parent_global']==rec(x.parent_global)
    assert g['source_scan']==s['scan']
    assert g['parent']['sha256']==s['envelope']['sha256']==old['envelope']['sha256']
    assert sha(ROOT/g['parent']['path'])==g['parent']['sha256']
    assert sha(ROOT/g['parent']['accepted_audit'])==g['parent']['audit_sha256']
    assert g['parent']['audit_sha256'] in {v['sha256']for v in old['accepted_geometry_proof_chain']}
    assert g['cells']==old['full_upper_priced_cells']==272579
    assert g['preserved_empty_cells']==old['proved_physically_empty_full_upper_cells']==8737
    assert g['ordinary_cells']+g['contact_cells']==g['cells']
    assert g['exact_profile_bindings']==s['selected_profiles']==43
    assert g['directed_precision_bits']>=768
    assert len(g['exact_R3_target_spectra'])==29 and all(v=='PASS_EXACT_SIX_PSD_GATES_AT_711'for v in g['exact_R3_target_spectra'].values())
    assert list(map(F,low['complete_amplitude_coverage']))==[F(0),F(1381,2000)]
    assert low['source_rho']=='847/1000' and low['native_cells']==39929 and low['R3_cells']==41
    maximum=F(g['maximum_independent_upper']);cap=F(g['strict_rational_cap'])
    assert maximum==F(563914263671,10**12)<cap==F(112782853,200000000)
    lower_cap=F(low['clean_rational_upper']);assert lower_cap<cap
    out={**old,'status':'PASS_INDEPENDENT_COMPLETE_B711_GLOBAL_SCALAR_CERTIFICATE',
         'operator_b':str(b),'norm_upper_bound':str(cap),
         'constructive_G_R_lower_bound':str(1/cap),'constructive_G_R_lower_bound_display':float(1/cap),
         'lower_amplitude_upper_bound':str(lower_cap),'full_upper_amplitude_upper_bound':str(maximum),
         'reference_envelope_operator_b':'713/1000',
         'bindings':{key:rec(getattr(x,key))for key in ('parent_global','transport','spectral','lower')},
         'analytic_review':'The accepted .713 physical partition and complete local gates are retained. A SINGLE new coefficient .711 is applied to all upper cells. Every ordinary completion shifts by beta_norm^2*(1/b-1/b0)+(b-b0)r^2; fixed contact completion replaces beta_norm by beta_norm+lambda*d/2. For the negative radial difference, each full physical cell is uniformly bounded at its nonnegative radial lower endpoint. All 43 exact beta/profile bindings and target spectral gates are independently accepted, and all 272579 shifted cell caps are independently recomputed. The 8737 physical empty cells do not depend on the operator coefficient. The separately accepted .711 lower branch covers [0,.6905], completing the same-coefficient global norm bound.',
         'scope':'Complete fresh all-dimension scalar norm and explicit matrix-sequence lower bound at .711. The reference envelope is the .713 proof source; every selected cell is explicitly transported to .711. No coefficient choices are mixed across cells and no historical1.773 premise is used.'}
    x.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ('status','operator_b','norm_upper_bound','constructive_G_R_lower_bound_display')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser()
    for key in ('parent-global','transport','spectral','lower','output'):
        a.add_argument('--'+key,required=True,type=Path)
    main(a.parse_args())
