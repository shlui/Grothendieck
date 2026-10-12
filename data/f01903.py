#!/usr/bin/env python3
"""Independent two-new-rank cache pricing and accepted whole-q bindings."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import argparse
import json
import os
import sys

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest
import replay_R3_corners as R


def fresh_tail(row,radius):
    p=row['p'];E00,E02,E22=R.even(row);q=p['major'][:]
    q[0]+=row['reserve0'];q[3]+=row['reserve6']
    ap=sum(abs(v)for v in p['phi']);ad=sum(j*j*abs(v)for j,v in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2))
    derivative=ap/radius**2+4*ad/radius**4
    stein=theta*ap+ap/radius**2+4*theta*ad/radius**2+4*ad/radius**4
    e0=abs(p['nuisance0']-p['nuisance2'])+derivative;e2=abs(p['nuisance2'])
    cc=p['c']-3*p['d0'];dd=p['d0'];h=p['zeta']**2/p['G'];v=p['y']-E02*p['r']
    l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0']
    l2=v+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    a3=abs(p['zeta']*(-1+dd/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    trace=E00+2*E22+p['R11']+6*p['R33']+h*(p['p0']**2+2*p['p2']**2)
    constant=2*(abs(p['ell'])+F(5,2)*abs(p['beta']))+trace
    poly=[F(0)]*7
    for j,value in enumerate(q):poly[2*j]=value
    poly[0]-=constant+ap+stein+e0*e0/(2*p['eta'])
    poly[1]-=2*(abs(l1)+abs(l2));poly[2]-=e0*e2/p['eta']+cc*cc/(4*p['G'])
    poly[3]-=2*a3;poly[4]-=e2*e2/(2*p['eta'])+abs(cc*dd)/(2*p['G'])
    poly[6]-=dd*dd/(4*p['G'])
    shifted=[sum((poly[j]*comb(j,k)*radius**(j-k)for j in range(k,7)),F(0))for k in range(7)]
    assert all(x>=0 for x in shifted)
    return dict(R=str(radius),difference_polynomial=list(map(str,poly)),shifted_coefficients=list(map(str,shifted)))


def main(source_path,output_path):
    digest=R.sha(source_path);source=R.load(source_path);assert F(source['operator_b'])==R.B
    priorpath=HERE/'ninebank_R3_corner_replay.json';prior=R.load(priorpath)
    assert prior['status']=='PASS_INDEPENDENT_NINEBANK_504_R3_CORNER_PRICES'
    prior_source=R.load(ROOT/prior['source']['path']);assert R.sha(ROOT/prior['source']['path'])==prior['source']['sha256']
    oldtables={(x['profile'],tuple(map(F,x['kappa_interval']))):x for x in prior_source['cached_corner_constants']}
    for name,c in prior['expanded_coefficients'].items():
        assert source['profile_coefficients'][name]=={k:c[k]for k in('r_quadratic','r_linear','ell')}
    for item in source['accepted_R3_binding']:
        path=Path(item['path']);path=path if path.is_absolute()else ROOT/path
        assert R.sha(path)==item['sha256']
    rows={};bindings={};target=ROOT/'tonight/operator/projection/targeted_certificate'
    oldrows=R.bank()
    for name in('lowr76','highk7575'):
        accepted=prior['complete_whole_q_bindings'][name]
        path=ROOT/accepted['profile']['path'];raw=R.load(path)
        oldrows[name]=dict(name=name,profile=path,p=R.exact(raw['profile']),
            kind='full'if raw['kind']=='EXACT_RATIONAL_GENERAL_R3_PROFILE'else'rank',
            reserve0=F(accepted['global_reserves'][0]),reserve6=F(accepted['global_reserves'][1]))
    for name,row in oldrows.items():
        accepted=prior['complete_whole_q_bindings'][name]
        assert R.sha(row['profile'])==accepted['profile']['sha256']
        row['coefficients']=R.coefficients(row,R.B)
        assert {k:str(v)for k,v in row['coefficients'].items()}==prior['expanded_coefficients'][name]
        R.spectral(row,R.B)
    for name,label,folder in(
        ('middle7475','interval_middle_u7475_r267_k08_1','targeted_middle_R3_recheck'),
        ('middlehigh7475','interval_middlehigh_u7475_r275_k11_2','targeted_middlehigh_R3_recheck')):
        path=target/f'profile_{label}.json';raw=R.load(path)
        assert raw['kind']=='EXACT_RATIONAL_R3_PROFILE'
        row=dict(name=name,kind='rank',profile=path,p=R.exact(raw['profile']),reserve0=F(1,10000),reserve6=F(1,1000000))
        ledgerpath=ROOT/f'tonight/audit/{folder}/ACCEPTANCE.json';ledger=R.load(ledgerpath)
        assert ledger['status'].startswith('PASS_INDEPENDENT_COMPLETE_TARGETED_')
        assert ledger['profile_sha256']==R.sha(path)
        assert list(map(F,ledger['reserves']))==[row['reserve0'],row['reserve6']]
        adapter=ROOT/'tonight/audit/replay_targeted_rank_q_trees.py'
        assert R.sha(adapter)==ledger['adapter_sha256']
        regions=[];previous=F(0)
        for part in('body','exterior'):
            reportpath=ROOT/ledger['reports'][part]['path'];report=R.load(reportpath)
            assert R.sha(reportpath)==ledger['reports'][part]['sha256']
            assert report['status']=='PASS_INDEPENDENT_PHYSICAL_Q_TREE'
            assert report['profile_sha256']==R.sha(path)
            treepath=target/f'{part}_{label}.json.gz';tree=R.load(treepath)
            assert R.sha(treepath)==report['tree_sha256']
            assert tree['status'].startswith('PASS')and tree['pending']==0
            assert tree['accepted']==report['accepted']and F(tree['qa'])==previous
            previous=F(tree['qb']);p=R.exact(tree['profile'])
            for key in row['p']:
                if key!='major':assert p[key]==row['p'][key]
            assert p['major']==[v+(row['reserve0']if j==0 else row['reserve6']if j==3 else F(0))for j,v in enumerate(row['p']['major'])]
            assert report['cubic_integral_backend_sha256']==R.sha(ROOT/'tonight/recovery/profile_cover/directed_cubic_integral.py')
            regions.append(dict(tree=R.file_record(treepath),accepted_proof=R.file_record(reportpath),q=[tree['qa'],tree['qb']]))
        assert previous==16
        tailpath=ROOT/ledger['reports']['tail']['path'];tail=R.load(tailpath)
        assert R.sha(tailpath)==ledger['reports']['tail']['sha256']
        assert tail['status']=='PASS_INDEPENDENT_EXACT_R3_TAIL'
        fresh=fresh_tail(row,F(16))
        for key in('R','difference_polynomial','shifted_coefficients'):assert fresh[key]==tail[key]
        coeff=R.coefficients(row,R.B);row['coefficients']=coeff
        assert source['profile_coefficients'][name]=={k:str(coeff[k])for k in('r_quadratic','r_linear','ell')}
        original=R.coefficients(row,row['p']['eps'])
        assert coeff['constant']-original['constant']==6*row['p']['beta']**2*(1/R.B-1/row['p']['eps'])
        bindings[name]=dict(profile=R.file_record(path),whole_q_acceptance=R.file_record(ledgerpath),regions=regions,
            tail_proof=R.file_record(tailpath),fresh_exact_tail=fresh,spectral_gates=R.spectral(row,R.B),
            expanded_coefficients={k:str(v)for k,v in coeff.items()},
            frozen_theta=str(coeff['theta']),actual_kappa_scope='Every actual kappa>=0; the local scalar gate has frozen theta and no actual cumulant.')
        rows[name]=row
    tables=[];keys=set();bands={name:set()for name in rows};oldcount=0
    for item in source['cached_corner_constants']:
        name=item['profile'];a,b=map(F,item['kappa_interval']);saved=F(item['upper'])
        assert(name,a,b)not in keys;keys.add((name,a,b))
        if name not in rows:
            assert item==oldtables[name,(a,b)];oldcount+=1;row=oldrows[name]
        else:
            bands[name].add((a,b));row=rows[name]
        corners=[]
        for k,z in product((a,b),(F(2,5),F(1))):
            actual=R.price(row,row['coefficients'],k,z);assert actual.below(saved)
            corners.append(dict(kappa=str(k),shape=str(z),price=actual.record(),
                plus_covariance_correction=str(row['p']['nuisance2']*(k/3-row['coefficients']['theta']*(k+2)))))
        tables.append(dict(profile=name,kappa_interval=list(map(str,(a,b))),saved_upper=str(saved),corners=corners))
    assert oldcount==126 and len(tables)==154 and len(keys)==154
    for name,intervals in bands.items():
        previous=F(0);assert len(intervals)==14
        for a,b in sorted(intervals):assert a==previous;previous=b
        assert previous==90
    assert R.sha(source_path)==digest
    out=dict(status='PASS_INDEPENDENT_ELEVENBANK_616_R3_CORNER_PRICES',source=R.file_record(source_path),
        engine=selftest(),operator_b=str(R.B),profiles=11,accepted_prior_corner_prices=504,
        fresh_new_corner_prices=112,fresh_all_corner_prices=616,all_cached_tables=154,prior_ninebank_audit=R.file_record(priorpath),
        old_cached_constants_exactly_unchanged=True,complete_new_whole_q_bindings=bindings,new_tables=tables,
        scope='R3 source pricing and complete whole-q bindings only. Exact-polygon physical coverage and the full global envelope need their separate independent audit. Frozen original and strengthened portable bundles unchanged.')
    output_path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in('status','source','profiles','fresh_all_corner_prices','fresh_new_corner_prices','all_cached_tables')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'tonight/recovery/profile_cover/global_middle_b7171_elevenbank_exactpolygon_envelope.json')
    parser.add_argument('--output',type=Path,default=HERE/'elevenbank_R3_corner_replay.json');args=parser.parse_args()
    main(args.source.resolve(),args.output.resolve())
