#!/usr/bin/env python3
"""Isolated next-grid producer-free all-profile R3 cache replay, including exact tail32."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import gzip
import json
import os
import sys

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest
import replay_R3_corners as R
from replay_middlebank_R3_corners import fresh_tail


def check_file_records(value):
    if isinstance(value,dict):
        if 'path'in value and 'sha256'in value:
            path=Path(value['path']);path=path if path.is_absolute()else ROOT/path
            assert R.sha(path)==value['sha256'],('FROZEN_PREMISE_SHA',path)
        for child in value.values():check_file_records(child)
    elif isinstance(value,list):
        for child in value:check_file_records(child)


def add_rank(name,folder):
    ledgerpath=ROOT/f'tonight/audit/{folder}/ACCEPTANCE.json';ledger=R.load(ledgerpath)
    assert ledger['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    path=ROOT/ledger['profile_path'];raw=R.load(path)
    assert raw['kind']=='EXACT_RATIONAL_R3_PROFILE'and R.sha(path)==ledger['profile_sha256']
    reserves=list(map(F,ledger['reserves']));assert len(reserves)==2 and all(v>=0 for v in reserves)
    row=dict(name=name,kind='rank',profile=path,p=R.exact(raw['profile']),reserve0=reserves[0],reserve6=reserves[1])
    # Bind the exact independent adapter used to produce the accepted proof.
    adapters=[p for p in(HERE.parent/'replay_targeted_rank_q_trees.py',
        HERE.parent/'replay_R3_robust_q_trees.py',HERE.parent/'replay_R3_robust_current_q_trees.py',
        HERE.parent/'replay_R3_highk_q_trees.py',HERE.parent/'replay_R3_lowerface_q_trees.py',
        HERE.parent/'replay_R3_nextface_q_trees.py',HERE.parent/'replay_R3_q50_q_trees.py',
        HERE.parent/'replay_R3_highu_q_trees.py',HERE.parent/'replay_R3_curve_q_trees.py')if p.exists()]
    matching=[p for p in adapters if R.sha(p)==ledger['adapter_sha256']];assert len(matching)==1
    if 'adapter_path'in ledger:assert ROOT/ledger['adapter_path']==matching[0]
    aliases=ledger.get('byte_identical_input_aliases',{})
    for rel,binding in aliases.items():
        assert R.sha(ROOT/rel)==R.sha(ROOT/binding['original'])==binding['sha256']
    regions=[];previous=F(0)
    for part in('body','exterior'):
        reportpath=ROOT/ledger['reports'][part]['path'];report=R.load(reportpath)
        assert R.sha(reportpath)==ledger['reports'][part]['sha256']
        assert report['status']=='PASS_INDEPENDENT_PHYSICAL_Q_TREE'and report['profile_sha256']==R.sha(path)
        filename=f'{part}_{report["profile"]}.json.gz';treepath=path.parent/filename
        if not treepath.exists():
            choices=[ROOT/rel for rel in aliases if Path(rel).name==filename]
            assert len(choices)==1,('MISSING_EXACT_TREE_ALIAS',name,filename)
            treepath=choices[0]
        tree=R.load(treepath)
        assert R.sha(treepath)==report['tree_sha256']
        assert tree['status'].startswith('PASS')and tree['pending']==0 and tree['accepted']==report['accepted']
        assert F(tree['qa'])==previous;previous=F(tree['qb']);p=R.exact(tree['profile'])
        for key in row['p']:
            if key!='major':assert p[key]==row['p'][key]
        assert p['major']==[v+(reserves[0]if j==0 else reserves[1]if j==3 else F(0))for j,v in enumerate(row['p']['major'])]
        assert report['cubic_integral_backend_sha256']==R.sha(ROOT/'tonight/recovery/profile_cover/directed_cubic_integral.py')
        regions.append(dict(tree=R.file_record(treepath),proof=R.file_record(reportpath),q=[tree['qa'],tree['qb']]))
    assert previous in(F(16),F(20),F(32))
    tailpath=ROOT/ledger['reports']['tail']['path'];tail=R.load(tailpath)
    assert R.sha(tailpath)==ledger['reports']['tail']['sha256']and tail['status']=='PASS_INDEPENDENT_EXACT_R3_TAIL'
    fresh=fresh_tail(row,previous)
    for key in('R','difference_polynomial','shifted_coefficients'):assert fresh[key]==tail[key]
    assert ledger['q_cover']==[['0','8'],['8',str(previous)],[str(previous),'infinity']]
    binding=dict(profile=R.file_record(path),accepted_whole_q_gate=R.file_record(ledgerpath),
        accepted_adapter=R.file_record(matching[0]),q_regions=regions,
        accepted_tail=R.file_record(tailpath),fresh_exact_tail=fresh,reserves=list(map(str,reserves)),
        byte_identical_input_aliases=aliases,
        spectral_gates=R.spectral(row,R.B))
    return row,binding


def expected_native_bands(source):
    records=[];bandsets=[]
    for item in source.get('native_endpoint_inputs',[]):
        path=ROOT/item['mask'];data=R.load(path)
        assert R.sha(path)==item['mask_sha256']
        assert data['status']=='PASS_DIRECTED_NATIVE_CELL_PRICING_BASE_LEAF_CERTIFICATE_REQUIRED'
        assert 0<R.B<=F(data['operator_b'])
        bands=sorted(tuple(map(F,x['kappa_interval']))for x in data['endpoint_prices'])
        assert len(bands)==len(set(bands))and bands[0][0]==0 and bands[-1][1]==90
        assert all(a[1]==b[0]for a,b in zip(bands,bands[1:]))
        records.append(dict(mask=R.file_record(path),native_table_operator_b=data['operator_b'],target_operator_b=str(R.B),kappa_bands=[list(map(str,x))for x in bands]))
        bandsets.append(set(bands))
    assert records,'Final cache requires frozen native endpoint ledger bindings'
    # A finer R3 partition may subdivide already accepted native bands. Its
    # explicit points are required; each native band endpoint must remain.
    points=source.get('kappa_points')or source.get('kappa_grid')
    if points is not None:
        points=list(map(F,points));assert points[0]==0 and points[-1]==90 and all(a<b for a,b in zip(points,points[1:]))
        for bands in bandsets:assert {x for band in bands for x in band}<=set(points)
        expected=set(zip(points,points[1:]))
        return expected,records,'explicit source grid refining every native band'
    assert all(bands==bandsets[0]for bands in bandsets),('NATIVE_BAND_PARTITIONS_DIFFER',list(map(len,bandsets)))
    return bandsets[0],records,'identical frozen native endpoint ledger bands'


def main(source_path,output_path):
    digest=R.sha(source_path);source=R.load(source_path)
    assert source['status'].startswith('PASS_')
    R.B=F(source['operator_b']);assert 0<R.B<=F(7171,10000)
    expected,gridbindings,gridreason=expected_native_bands(source)
    rows=R.bank();bindings={name:R.binding(row)for name,row in rows.items()}
    priorpath=HERE/'ninebank_R3_corner_replay.json';prior=R.load(priorpath)
    assert prior['status']=='PASS_INDEPENDENT_NINEBANK_504_R3_CORNER_PRICES'
    for name in('lowr76','highk7575'):
        accepted=prior['complete_whole_q_bindings'][name];check_file_records(accepted)
        path=ROOT/accepted['profile']['path'];raw=R.load(path)
        rows[name]=dict(name=name,profile=path,p=R.exact(raw['profile']),
            kind='full'if raw['kind']=='EXACT_RATIONAL_GENERAL_R3_PROFILE'else'rank',
            reserve0=F(accepted['global_reserves'][0]),reserve6=F(accepted['global_reserves'][1]))
        bindings[name]=dict(accepted_prior_ninebank_audit=R.file_record(priorpath),frozen_gate_binding=accepted,
            target_spectral_gates=R.spectral(rows[name],R.B))
    new_ranks=(('middle7475','targeted_middle_R3_recheck'),('middlehigh7475','targeted_middlehigh_R3_recheck'),
        ('interval_robust_middlehigh_box','robust_middlehigh_R3_recheck'),
        ('interval_robust_currentworst_box','robust_currentworst_R3_recheck'),
        ('interval_robust_broadcurrent_box','robust_broadcurrent_R3_recheck'),
        ('robust_broadcurrent','robust_broadcurrent_R3_recheck'),
        ('moderatek','robust_moderatek_R3_recheck'),
        ('highk7625','highk_lowr_R3_recheck'),
        ('secondary7475','secondary_highu_R3_recheck'),
        ('primary7425','primary_u7425_R3_recheck'),
        ('lowerface7175','lowerface_u7175_R3_recheck'),
        ('nextface7375','nextface_u7375_R3_recheck'),
        ('lowerdiffuse72125','lowerdiffuse_u72125_R3_recheck'),
        ('moderatefaceq50','moderateface_q50_R3_recheck'),
        ('moderateface7375q50','moderateface_q50_R3_recheck'),
        ('finalface7475','finalface_u7475_R3_recheck'),
        ('highu7575','highu_highk_R3_recheck'),
        ('highu7775','highu_highk_R3_recheck'),
        ('lastface7575','lastface_u7575_R3_recheck'),
        ('refinedface75125','refinedface_u75125_R3_recheck'),
        ('curve12_20','curve12_20_R3_recheck'))
    required=set(source['profile_coefficients'])
    for name,folder in new_ranks:
        if name in required:rows[name],bindings[name]=add_rank(name,folder)
    assert required<=set(rows),('UNKNOWN_PROFILE',required-set(rows))
    rows={name:row for name,row in rows.items()if name in required}
    bindings={name:value for name,value in bindings.items()if name in required}
    for item in source['accepted_R3_binding']:check_file_records(item)
    coefficients={}
    for name,row in rows.items():
        c=R.coefficients(row,R.B);row['coefficients']=c
        assert source['profile_coefficients'][name]=={k:str(c[k])for k in('r_quadratic','r_linear','ell')}
        assert c['A']>0 and c['B']>0
        E00,E02,E22=R.even(row);p=row['p']
        for r in(F(0),F(1,3),F(1)):
            y=p['y']+E02*(r-p['r'])
            assert y-E02*r==c['yc']
            assert (R.B-E00)*r*r-2*p['gamma']*r+y*y/E22==c['r_quadratic']*r*r+c['r_linear']*r+c['yc']**2/E22
        original=R.coefficients(row,p['eps'])
        assert c['constant']-original['constant']==6*p['beta']**2*(1/R.B-1/p['eps'])
        coefficients[name]={k:str(v)for k,v in c.items()}
    seen=set();bands={name:set()for name in rows};tables=[]
    for item in source['cached_corner_constants']:
        name=item['profile'];a,b=map(F,item['kappa_interval']);saved=F(item['upper'])
        assert name in rows and(name,a,b)not in seen and(a,b)in expected
        seen.add((name,a,b));bands[name].add((a,b));row=rows[name];corners=[]
        for k,z in product((a,b),(F(2,5),F(1))):
            actual=R.price(row,row['coefficients'],k,z);assert actual.below(saved),(name,a,b,k,z,actual.record())
            corners.append(dict(kappa=str(k),shape=str(z),price=actual.record(),
                plus_covariance=str(row['p']['nuisance2']*(k/3-row['coefficients']['theta']*(k+2)))))
        tables.append(dict(profile=name,kappa_interval=list(map(str,(a,b))),saved_upper=str(saved),corners=corners))
    assert all(intervals==expected for intervals in bands.values())
    assert len(tables)==len(rows)*len(expected)and R.sha(source_path)==digest
    out=dict(status='PASS_INDEPENDENT_GENERIC_COMPLETE_R3_KAPPA_CACHE',source=R.file_record(source_path),engine=selftest(),
        operator_b=str(R.B),profiles=len(rows),kappa_bands=len(expected),corner_tables=len(tables),
        directed_corner_prices=4*len(tables),expected_kappa_bands=[list(map(str,x))for x in sorted(expected)],
        grid_derivation=gridreason,native_grid_bindings=gridbindings,complete_whole_q_bindings=bindings,
        expanded_coefficients=coefficients,corners=tables,
        convexity='Separate cumulant/shape convexity is established by A>0, B>0 and 194/sqrt(10)-633/15>0. Frozen theta, positive nuisance covariance and exact y(r) transport are reconstructed from each profile. Native/subdivided grid equality and complete profile-by-band coverage are checked exactly.',
        scope='All R3 cache prices and accepted whole-q gates only. Exact physical polygon/native/plane coverage and global bound assembly require their separate independent acceptance. Existing portable bundles unchanged.')
    content=json.dumps(out,indent=2)+'\n'
    if output_path.suffix=='.gz':
        with gzip.open(output_path,'wt')as stream:stream.write(content)
    else:output_path.write_text(content)
    print(json.dumps({k:out[k]for k in('status','source','profiles','kappa_bands','corner_tables','directed_corner_prices')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path)
    parser.add_argument('--output',type=Path,default=HERE/'next57_R3_corner_replay.json');args=parser.parse_args()
    main(args.source.resolve(),args.output.resolve())
