#!/usr/bin/env python3
"""Producer-free coupled native/R3/contact replay, fixed lambda per cell."""
from bisect import bisect_right
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import sys
import time

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
import replay_native_R3_unions as N
from mpfr_interval import I,selftest
import replay_R3_corners as R

SOURCE=ROOT/'tonight/operator/contact_credit_joint/new_full_linear_pairing_with_snapshots.json'
OUTPUT=HERE/'joint_contact_full_718_independent.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def record(path):return dict(path=str(path.relative_to(ROOT)),sha256=sha(path))
def key(cell):return tuple(tuple(map(F,cell[field]))for field in('u_interval','r_interval','kappa_interval'))
def bound(binding):
    path=Path(binding['path'])
    if not path.is_absolute():path=ROOT/path
    assert sha(path)==binding['sha256'],('BOUND_HASH',path)
    return path,read(path)


def native_curve(frame,binding,cell):
    ulo,uhi=map(F,cell['u_interval']);rlo,rhi=map(F,cell['r_interval'])
    uhi=min(uhi,F(cell['physical_u_upper']));rhi=min(rhi,F(cell['physical_r_upper']))
    klo,khi=map(F,cell['kappa_interval']);rg=frame['rg']
    assert rg[0]<=rlo<rhi<=rg[-1]and ulo<=uhi
    ep=frame['endpoints'][klo,khi]
    splits=[rlo]+[r for r in rg if rlo<r<rhi]+[rhi]
    curves=[]
    for left,right in zip(splits,splits[1:]):
        j=min(bisect_right(rg,left)-1,len(rg)-2)
        assert rg[j]<=left<right<=rg[j+1]
        for corner in range(4):
            linear=(ep[j+1][corner]-ep[j][corner])/(rg[j+1]-rg[j])
            constant=ep[j][corner]-linear*rg[j]
            curves.append(dict(corner=corner,constant=str(constant),linear=str(linear),
                quadratic=str(frame['g']),r_interval=list(map(str,(left,right)))))
    p=N.package(binding['package'])['profile']
    return dict(profile='native:'+binding['mask'],kind='native',ell=str(frame['ell']),
        beta=str(F(p['beta'])),u_interval=list(map(str,(ulo,uhi))),
        kappa_interval=list(map(str,(klo,khi))),curves=curves)


def contact_row(source,credit):
    index=credit['contact_row'];scale=F(credit['contact_scale']);assert scale>0
    if credit['contact_family']=='cubic':
        row=source['cubic']['rows'][index]
        q=[F(v)*scale**(1-2*j)for j,v in enumerate(row['polynomial_coefficients_q0_q2_q4_q6'])]
        a=F(source['cubic']['parameters']['a'])*scale
        B=F(source['cubic']['parameters']['B'])*scale
    else:
        assert credit['contact_family']=='linear'
        row=source['linear']['rows'][index];assert F(row['scale'])==scale
        q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']))
        a=F(row['a']);B=F(0)
    assert F(credit['contact_w'])==F(row['w'])and a>0 and B>=0 and q[3]>=0
    expected=dict(a=str(a),B=str(B),d_squared=str(6*B*B),q0_q2_q4_q6=list(map(str,q)))
    assert credit['contact_exact_parameters']==expected
    return a,B,q


def main():
    start=time.time();digest=sha(SOURCE);data=read(SOURCE)
    assert data['status']=='PASS_EXACT_FIXED_MULTIPLIER_PAIRED_CORNER_REPRICING_REQUIRES_INDEPENDENT_ASSEMBLY_AUDIT'
    assert not data['worst_only'];b=F(data['operator_b']);assert b==F(7171,10000)
    basepath,base=bound(data['envelope_binding'])
    assert sha(basepath)=='efda4aa16d57edc8c67eb94c4121a376b026ee8a52d141de8dace0460711f561'
    assert base['status']=='PASS_CONTINUOUS_ALL_SOURCE_NATIVE_R3_ENVELOPE_PRICING'
    assert F(base['operator_b'])==b
    mainpath,mainledger=bound(data['source_assembly_binding'])
    assert len(mainledger['cells'])==len(data['cells'])==len(base['cells'])==39351
    for fresh,saved in zip(data['cells'],mainledger['cells']):
        stripped=dict(fresh)
        if stripped.get('joint_credit_improved'):
            stripped['joint_contact_credit']=dict(stripped['joint_contact_credit'])
            for field in('source_curve_snapshot','contact_exact_parameters'):
                if field not in saved['joint_contact_credit']:stripped['joint_contact_credit'].pop(field)
        assert stripped==saved,('SNAPSHOT_CHANGED_MAIN_CELL',key(fresh))
    sources={};proofs={}
    for name,field,proof in(
        ('cubic','contact_binding','tonight/audit/cubic_head_contact_independent.json'),
        ('linear','linear_bank_binding','tonight/audit/linear_contact_bank_independent.json')):
        path,sources[name]=bound(data[field]);proofpath=ROOT/proof;accepted=read(proofpath)
        assert accepted['input_sha256']==sha(path)and accepted['status'].startswith('PASS_INDEPENDENT')
        proofs[name]=record(proofpath)
    assert sources['cubic']['status']=='PASS_EXACT_ALL_HALFLINE_DIRECT_CUBIC_HEAD_CONTACT_MAJORANT'
    assert sources['linear']['status']=='PASS_EXACT_CONTACT_PHYSICAL_SOURCE_BANK_AND_FROZEN_BOX_EXCLUSIONS'
    # Reconstruct every used native frame from primitive packages and fresh
    # directed prices, rather than trusting saved source-curve snapshots.
    wanted=Counter(cell['joint_contact_credit']['source_curve_snapshot']['profile']
        for cell in data['cells']if cell.get('joint_credit_improved'))
    frames={};native_bindings={}
    for binding in data['native_bindings']:
        assert sha(ROOT/binding['mask'])==binding['mask_sha256']
        package=N.package(binding['package'])
        assert package['source_sha256']==binding['source_sha256']
        assert package['replay_sha256']==binding['replay_sha256']
        if 'native:'+binding['mask']in wanted:
            frames[binding['mask']]=N.frame(binding,b)
            native_bindings[binding['mask']]=binding
            print(json.dumps(dict(stage='FRESH_NATIVE_FRAME_PASS',frame=binding['mask'],point_prices=N.POINTS)),flush=True)
    # Bind every used R3 profile to its independent whole-q proof. The
    # actual cumulant is free in the accepted scalar gate; theta is frozen.
    accepted_bindings=[]
    for item in data['accepted_R3_binding']:
        path=Path(item['path']);path=path if path.is_absolute()else ROOT/path
        assert sha(path)==item['sha256'];accepted_bindings.append(record(path))
    bank=R.bank();r3={};r3_bindings={}
    for name,row in bank.items():
        if 'R3:'+name not in wanted:continue
        r3_bindings[name]=R.binding(row);coeff=R.coefficients(row,b)
        assert base['profile_coefficients'][name]=={k:str(coeff[k])for k in('r_quadratic','r_linear','ell')}
        assert coeff['A']>0 and coeff['B']>0
        r3[name]=(row,coeff)
    assert set(wanted)=={'native:'+name for name in frames}|{'R3:'+name for name in r3}
    original={key(cell):cell for cell in base['cells']};assert len(original)==39351
    seen=set();results=[];families=Counter();unchanged=0;curvecount=0;freshcontact=0
    contact_cache={};r3_cache={}
    mutable={'cap','winner','target_pass','fallback_pass'}
    additions={'joint_credit_improved','joint_credit_best_new_cap','joint_contact_credit','previous_cap'}
    for number,cell in enumerate(data['cells']):
        cellkey=key(cell);assert cellkey not in seen and cellkey in original;seen.add(cellkey)
        old=original[cellkey]
        assert {k:v for k,v in cell.items()if k not in mutable|additions}=={k:v for k,v in old.items()if k not in mutable}
        if not cell.get('joint_credit_improved'):
            assert {k:v for k,v in cell.items()if k not in additions}==old
            unchanged+=1;continue
        credit=cell['joint_contact_credit'];families[credit['contact_family']]+=1
        assert F(cell['previous_cap'])==F(old['cap'])and F(cell['cap'])<F(old['cap'])
        assert credit['corner_order']==['klo,s0','klo,s1','khi,s0','khi,s1']
        snap=credit['source_curve_snapshot'];klo,khi=map(F,cell['kappa_interval'])
        ulo,uhi=map(F,cell['u_interval']);rlo,rhi=map(F,cell['r_interval'])
        uhi=min(uhi,F(cell['physical_u_upper']));rhi=min(rhi,F(cell['physical_r_upper']))
        assert ulo<=uhi and rlo<rhi
        name=snap['profile'];assert cell['winner']=='paired-contact:'+name
        if snap['kind']=='native':
            mask=name[len('native:'):]
            rebuilt=native_curve(frames[mask],native_bindings[mask],cell)
            assert rebuilt==snap,('NATIVE_SNAPSHOT',number)
        else:
            assert snap['kind']=='R3';profile=name[len('R3:'):];row,c=r3[profile]
            assert snap['ell']==str(c['ell'])and snap['beta']==str(row['p']['beta'])
            assert tuple(map(F,snap['u_interval']))==(ulo,uhi)
            assert tuple(map(F,snap['kappa_interval']))==(klo,khi)
            assert len(snap['curves'])==4
            for corner,(k,s)in enumerate(product((klo,khi),(F(0),F(1)))):
                curve=snap['curves'][corner];assert curve['corner']==corner
                assert F(curve['linear'])==c['r_linear']and F(curve['quadratic'])==c['r_quadratic']
                assert tuple(map(F,curve['r_interval']))==(rlo,rhi)
                cachekey=profile,k,s
                if cachekey not in r3_cache:r3_cache[cachekey]=R.price(row,c,k,F(2,5)+F(3,5)*s)
                assert r3_cache[cachekey].below(F(curve['constant'])),('R3_SNAPSHOT_PRICE',number,corner)
        assert tuple(map(F,snap['u_interval']))==(ulo,uhi)
        assert tuple(map(F,snap['kappa_interval']))==(klo,khi)
        a,B,q=contact_row(sources,credit);lam=F(credit['lambda_value']);assert lam>0
        beta=F(snap['beta']);saved_beta=F(credit['beta_d_over_b_upper'])
        if B==0:assert saved_beta==0
        elif snap['kind']=='R3':assert 6*beta*B/b<=saved_beta,('BETA_COMPLETION',number)
        else:
            beta_actual=I(beta*B/b)*I(6).sqrt()
            assert beta_actual.below(saved_beta),('BETA_COMPLETION',number)
        if credit['contact_family']=='linear':assert B==saved_beta==0
        saved_contact=list(map(F,credit['contact_corner_prices']));assert len(saved_contact)==4
        for corner,(k,s)in enumerate(product((klo,khi),(F(0),F(1)))):
            cachekey=credit['contact_family'],credit['contact_row'],F(credit['contact_scale']),k,s
            if cachekey not in contact_cache:
                contact_cache[cachekey]=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(N.K-F(633,15)*s)*I(k)*I(k).sqrt())
                freshcontact+=1
            assert contact_cache[cachekey].below(saved_contact[corner]),('CONTACT_CORNER',number,corner)
        ell=F(snap['ell']);outer=max(u*u-2*ell*u-lam*a*u for u in(ulo,uhi))
        lam2=lam*lam*6*B*B/(4*b);values=[]
        for curve in snap['curves']:
            corner=curve['corner'];left,right=map(F,curve['r_interval'])
            assert rlo<=left<right<=rhi and 0<=corner<4
            constant=F(curve['constant'])+lam*(saved_contact[corner]+saved_beta)+lam2
            linear=F(curve['linear'])-lam;g=F(curve['quadratic'])
            values.append(N.qmax(constant,linear,g,left,right)+outer);curvecount+=1
        maximum=max(values);cap=N.decimal_ceiling(maximum,12)
        assert maximum==F(credit['exact_unrounded_cap']),('EXACT_CELL_MAX',number)
        assert cap==F(cell['cap']),('ROUNDED_CAP',number)
        assert bool(cell['target_pass'])==(cap<F(40,71))
        assert bool(cell['fallback_pass'])==(cap<F(451,800))
        results.append(dict(cell_index=number,source_profile=name,contact_family=credit['contact_family'],
            row=credit['contact_row'],scale=credit['contact_scale'],fixed_lambda=str(lam),
            exact_unrounded_cap=str(maximum),rounded_cap=str(cap),
            previous_cap=cell['previous_cap'],snapshot_sha256=hashlib.sha256(json.dumps(snap,sort_keys=True).encode()).hexdigest()))
    assert seen==set(original)and len(results)==data['snapshot_count']==data['improved_cells']==718
    assert families==Counter(cubic=341,linear=377)and unchanged==38633
    worst=max(data['cells'],key=lambda x:F(x['cap']))
    assert worst==data['worst_cell']and F(worst['cap'])==F(data['maximum_cell_cap'])==F('0.564828452412')
    assert sha(SOURCE)==digest
    output=dict(status='PASS_INDEPENDENT_FULL_718_FIXED_LAMBDA_PAIRED_CONTACT_ASSEMBLY',
        source=record(SOURCE),base_envelope=record(basepath),source_assembly=record(mainpath),engine=selftest(),
        producer_modules_imported=False,operator_b=str(b),baseline_cells_retained=39351,
        unchanged_baseline_cells=unchanged,improved_cells=len(results),contact_families=dict(families),
        unique_fresh_contact_corners=freshcontact,unique_fresh_R3_corners=len(r3_cache),
        fresh_native_endpoint_prices=N.POINTS,exact_quadratic_curve_maxima=curvecount,
        maximum_cell_cap=str(F(worst['cap'])),worst_cell=worst,
        native_package_bindings=list(native_bindings.values()),accepted_contact_proofs=proofs,
        accepted_R3_ledger=accepted_bindings,used_R3_whole_q_bindings=r3_bindings,
        source_choice_counts=dict(wanted),selected_cell_replays=results,
        coupling='Each selected nonnegative rational multiplier, contact row and positive scale are fixed across its complete original cell. Native/R3 and contact prices use the same kappa/shape corner before maximization.',
        completion='The positive nuisance completion terms are lambda*beta*d/b+lambda^2*d^2/(4*b), with native beta normalization d=sqrt(6)*B and rank/full raw beta term 6*beta*B/b. Linear B0 has exactly zero d and beta correction.',
        convexity='Native effective sixth coefficients and contact sixth coefficients are nonnegative. R3 has positive kappa quadratic coefficient and nonnegative kappa^(3/2) coefficient on all shapes. The fixed contact augmentation preserves separate kappa convexity and shape convexity (native shape affine). The u quadratic is convex; exact r maximization includes any concave interior vertex. No convexity of the minimum of frames is used.',
        scope='Assembly and baseline retention are accepted conditional on the exact bound base envelope efda4aa1… and its complete physical-domain/primitive acceptance. This report does not independently certify new base coverage or claim a final global value for G_R.')
    output['seconds']=time.time()-start
    OUTPUT.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k]for k in('status','baseline_cells_retained','unchanged_baseline_cells','improved_cells','contact_families','fresh_native_endpoint_prices','unique_fresh_contact_corners','unique_fresh_R3_corners','maximum_cell_cap','seconds')},indent=2))


if __name__=='__main__':main()
