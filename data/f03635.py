#!/usr/bin/env python3
"""Producer-free exact beta, coefficient and spectral bindings for b transport."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,sys,os
os.environ['GR_MPFR_BITS']='768'
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
import replay_R3_corners as R

def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def rec(p):return {'path':str(p.resolve().relative_to(ROOT)),'sha256':N.sha(p)}

def main(x):
    scan=read(x.scan);g=read(x.parent_global);src=ROOT/scan['parent_path'];d=read(src)
    assert g['status']=='PASS_INDEPENDENT_COMPLETE_MONOTONE_REFINED_GLOBAL_SCALAR_CERTIFICATE'
    assert N.sha(src)==scan['parent_sha256']==g['envelope']['sha256']
    b0=F(d['operator_b']);b=x.b;assert 0<b<b0==F(g['operator_b'])
    catalog={};gate_bindings=[];walk=g
    while walk['status']=='PASS_INDEPENDENT_COMPLETE_MONOTONE_REFINED_GLOBAL_SCALAR_CERTIFICATE':
        for binding in walk['monotone_bindings']['new_complete_whole_q_gates']:
            gatepath=ROOT/binding['path'];assert N.sha(gatepath)==binding['sha256']
            gate=read(gatepath);assert gate['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
            path=ROOT/gate['profile_path'];assert N.sha(path)==gate['profile_sha256']
            catalog[binding['profile']]={'profile':rec(path),'kind':'rank','target_reserves':gate['reserves']}
            gate_bindings.append(rec(gatepath))
        p=walk['monotone_bindings']['parent_global'];pp=ROOT/p['path'];assert N.sha(pp)==p['sha256'];walk=read(pp)
    cachepath=ROOT/g['bindings']['cache']['path'];assert N.sha(cachepath)==g['bindings']['cache']['sha256']
    cache=read(cachepath);assert cache['status'].startswith('PASS_INDEPENDENT_')
    catalog.update(cache['complete_whole_q_bindings'])
    for name,binding in list(catalog.items()):
        if 'profile' not in binding:
            prior=binding['accepted_prior_ninebank_audit'];assert N.sha(ROOT/prior['path'])==prior['sha256']
            frozen=binding['frozen_gate_binding'];path=ROOT/frozen['profile']['path']
            assert N.sha(path)==frozen['profile']['sha256']
            kind='full' if 'E00' in read(path)['profile'] else 'rank'
            catalog[name]={'profile':frozen['profile'],'kind':kind,'target_reserves':frozen['global_reserves']}
        else:
            binding=dict(binding)
            if 'target_reserves' not in binding:binding['target_reserves']=binding['reserves']
            if 'kind' not in binding:binding['kind']='full' if 'E00' in read(ROOT/binding['profile']['path'])['profile'] else 'rank'
            if 'accepted_whole_q_gate' in binding:
                accepted=binding['accepted_whole_q_gate'];assert N.sha(ROOT/accepted['path'])==accepted['sha256']
            catalog[name]=binding
    assert set(catalog)==set(d['profile_coefficients'])
    reported=scan['profile_beta_bindings'];rows=[]
    for name,binding in catalog.items():
        path=ROOT/binding['profile']['path'];assert N.sha(path)==binding['profile']['sha256']
        p=R.exact(read(path)['profile']);r0,r6=map(F,binding['target_reserves'])
        row={'kind':binding['kind'],'p':p,'reserve0':r0,'reserve6':r6}
        old=R.coefficients(row,b0);new=R.coefficients(row,b);spectral=R.spectral(row,b)
        for key in ('r_quadratic','r_linear','ell'):
            assert old[key]==F(d['profile_coefficients'][name][key])
        assert new['r_quadratic']==old['r_quadratic']+b-b0
        assert new['r_linear']==old['r_linear'] and new['ell']==old['ell']
        beta2=6*p['beta']**2;key='R3:'+name;rb=reported[key]
        assert F(rb['beta_norm_squared'])==beta2
        assert N.sha(ROOT/rb['profile_path'])==rb['profile_sha256']==N.sha(path)
        assert len(rb['all_matching_profiles'])==1
        match=rb['all_matching_profiles'][0]
        assert N.sha(ROOT/match['path'])==match['sha256']==N.sha(path)
        rows.append({'name':key,'profile':rec(path),'normalized_beta_squared':str(beta2),'spectral':spectral})
    native=d['native_endpoint_inputs']+d['native_kappa_refinement']['new_native_bindings']
    for binding in native:
        package=N.package(binding['package']);p=package['profile'];M=package['M'];beta=F(p['beta'])
        assert b<=F(p['epsilon'])
        for a,o,c in ((M[0][0]+1,M[0][2],M[2][2]-b),(M[1][1]+1,M[1][3],M[3][3]-b),(M[1][1]-b,M[1][3],M[3][3])):N.psd(a,o,c)
        key='native:'+binding['mask'];rb=reported[key];path=ROOT/binding['package']/'data/local_summary.json'
        assert F(rb['beta_norm_squared'])==beta**2
        assert N.sha(path)==rb['profile_sha256']==N.sha(ROOT/rb['profile_path'])
        rows.append({'name':key,'profile':rec(path),'normalized_beta_squared':str(beta**2),'target_PSD_checked':True})
    assert len(rows)==len(reported)==43 and {v['name']for v in rows}==set(reported)
    out={'status':'PASS_INDEPENDENT_ALL_ASSIGNED_PROFILE_BETA_SPECTRAL_BINDINGS',
         'operator_b':str(b),'baseline_b':str(b0),'parent_global':rec(x.parent_global),
         'scan':rec(x.scan),'envelope':rec(src),'bindings':rows,'selected_profiles':43,
         'accepted_R3_cache':rec(cachepath),'accepted_new_whole_q_gates':gate_bindings,
         'producer_modules_imported':False,
         'scope':'Exact profile normalization and target-coefficient spectral premises. Complete all-cell transport and same-b lower-branch assembly remain separate.'}
    x.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(rows))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--scan',type=Path,required=True);a.add_argument('--parent-global',type=Path,required=True)
    a.add_argument('--b',type=F,required=True);a.add_argument('--output',type=Path,required=True);main(a.parse_args())
