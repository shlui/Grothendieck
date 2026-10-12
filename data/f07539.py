"""Producer-free complete assigned-completion transport: exact Fraction + MPFR768."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import json,sys
from audit_fullupper_polygon_envelope import ROOT,read,sha,ceil12
from audit_general_r3_parameters_tail import gates
from mpfr_interval import I,UP

parentpath=ROOT/'tonight/recovery/profile_cover/global_upper_b713_native_kappa68_twentyninebank_envelope.json.gz'
scanpath=ROOT/'tonight/recovery/profile_cover/assigned_twenty-nine_coefficient_scan.json'
parentaudit=ROOT/'tonight/theory/fullupper_b713_native_kappa68_twentyninebank_monotone_audit.json'
d=read(parentpath);scan=read(scanpath);pa=read(parentaudit)
assert pa['status']=='PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE'and pa['source_sha256']==sha(parentpath)
assert scan['parent_sha256']==sha(parentpath)and Path(scan['parent_path']).resolve()==parentpath.resolve()
b0=F(d['operator_b']);b=F(711,1000);assert b0==F(713,1000)
delta=b-b0;reciprocal=1/b-1/b0;assert delta<0 and reciprocal>0
bindings={};spectra={}
for name,item in scan['profile_beta_bindings'].items():
    path=ROOT/item['profile_path'];assert sha(path)==item['profile_sha256'];raw=read(path)['profile'];beta=F(raw['beta'])
    if name.startswith('native:'):
        normsq=beta**2;assert F(raw['epsilon'])>=b0>b
    else:
        assert name.startswith('R3:');normsq=6*beta**2
        p={k:list(map(F,v))if isinstance(v,list)else F(v)for k,v in raw.items()};assert p['eps']>=b0
        p['eps']=b
        if 'E00'not in p:p.update(E00=p['p0']**2,E02=p['p0']*p['p2'],E22=p['p2']**2)
        gates(p);spectra[name]='PASS_EXACT_SIX_PSD_GATES_AT_711'
    assert normsq==F(item['beta_norm_squared']);bindings[name]=(beta,normsq,path,sha(path))
assert len(bindings)==43
contactpath=ROOT/d['contact_gate_binding']['path'];assert sha(contactpath)==d['contact_gate_binding']['sha256'];contact=read(contactpath)
B=F(contact['parameters']['B']);assert B>0;d2=6*B**2
maximum=F(0);worst=None;ordinary=0;contacts=0;inherited=0;names=Counter();shifts={}
for index,row in enumerate(d['cells']):
    name=row['winner']
    while name.startswith('inherited:'):name=name[len('inherited:'):];inherited+=1
    iscontact=name.startswith(('native-contact:','R3-contact:'))
    name=name.replace('native-contact:','native:',1).replace('R3-contact:','R3:',1)
    assert name in bindings;beta,normsq,_,_=bindings[name]
    lo,hi=map(F,row['r_interval']);hi=min(hi,F(row['physical_r_upper']));assert 0<=lo<=hi
    shift=normsq*reciprocal+delta*lo*lo
    if iscontact:
        lam=F(row['contact_credit']['lambda_value']);assert lam>=0
        if name.startswith('native:'):
            key=(beta,lam)
            if key not in shifts:
                cross=I(reciprocal)*(I(lam*beta*B)*I(6).sqrt()+I(lam*lam*d2/4))
                shifts[key]=F(cross.hi.text(UP,240))
            extra=shifts[key]
        else:extra=reciprocal*(lam*6*beta*B+lam*lam*d2/4)
        shift+=extra;contacts+=1
    else:ordinary+=1
    cap=ceil12(F(row['cap'])+shift);names[name]+=1
    if cap>maximum:maximum=cap;worst={'parent_index':index,'u_interval':row['u_interval'],'r_interval':row['r_interval'],'kappa_interval':row['kappa_interval'],'winner':row['winner'],'independent_upper':str(cap)}
entry=[v for v in scan['results']if F(v['operator_b'])==b];assert len(entry)==1;entry=entry[0]
assert entry['failed_spectral_profiles']==[]and entry['triangle_fallback_cells']==0
assert maximum==F(entry['maximum'])<F(entry['strict_upper_cap'])==F(112782853,200000000)
assert len(d['cells'])==ordinary+contacts==272579
out={'status':'PASS_PRODUCER_FREE_COMPLETE_ASSIGNED_B711_CELL_TRANSPORT','operator_b':str(b),'baseline_b':str(b0),'parent':{'path':str(parentpath.relative_to(ROOT)),'sha256':sha(parentpath),'accepted_audit':str(parentaudit.relative_to(ROOT)),'audit_sha256':sha(parentaudit)},'source_scan':{'path':str(scanpath.relative_to(ROOT)),'sha256':sha(scanpath)},'exact_profile_bindings':len(bindings),'directed_precision_bits':768,'ordinary_cells':ordinary,'contact_cells':contacts,'inherited_prefixes_resolved':inherited,'cells':len(d['cells']),'preserved_empty_cells':len(d['physically_empty_cells']),'fresh_native_contact_cross_values':len(shifts),'exact_R3_target_spectra':spectra,'maximum_independent_upper':str(maximum),'strict_rational_cap':entry['strict_upper_cap'],'worst_cell':worst,'proof':'Every fixed physical cell uses the same operator .711. The frozen ordinary completion changes by delta_b*r^2+beta_norm^2*delta_reciprocal. Fixed contact completion adds its exact effective-beta cross/square terms. Since delta_b<0, r_lower gives a uniform radial upper. Accepted parent caps bound all physical points, including inherited children. All original physical exclusions and partitions are unchanged.'}
outpath=ROOT/'tonight/theory/assigned_b711_shifted_cells_audit.json';outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
