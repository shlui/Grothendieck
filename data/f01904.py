#!/usr/bin/env python3
"""Fresh all-nine bank R3 corner prices for the separate upgraded envelope."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import argparse
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('independent_R3_expansion',HERE/'replay_R3_corners.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)


def main(source_path,output_path):
    digest=R.sha(source_path);source=R.load(source_path);rows=R.bank()
    assert F(source['operator_b'])==R.B
    targeted=R.load(ROOT/'tonight/operator/projection/targeted_certificate/targeted_source_boxes.json')
    outer=R.load(HERE/'targeted_32_outer_vertices.json')
    assert outer['status']=='PASS_INDEPENDENT_TARGETED_32_OUTER_SOURCE_VERTICES'
    assert outer['source']['sha256']==R.sha(ROOT/'tonight/operator/projection/targeted_certificate/targeted_source_boxes.json')
    bindings={name:R.binding(row)for name,row in rows.items()}
    for name,box_name in (('lowr76','lowr_center_box'),('highk7575','highk_global_worst_box')):
        box=next(x for x in targeted['boxes']if x['name']==box_name)
        path=ROOT/'tonight/operator/projection/targeted_certificate'/box['profile']
        row=dict(name=name,p=R.exact(R.load(path)['profile']),profile=path,kind=box['kind'],
                 reserve0=F(box['reserve0']),reserve6=F(box['reserve6']))
        accepted=next(x for x in outer['boxes']if x['name']==box_name)
        assert accepted['profile']['sha256']==R.sha(path)
        assert accepted['whole_q_acceptance']['sha256']==R.sha(ROOT/box['audit'])
        for rel,hh in box['sha256'].items():assert R.sha(ROOT/rel)==hh
        audit=R.load(ROOT/box['audit']);assert audit['status'].startswith('PASS_INDEPENDENT')
        if row['kind']=='full':
            assert F(audit['global_reserve0'])==row['reserve0']
            assert F(audit['global_reserve6'])==row['reserve6']
            for rel,hh in audit['frozen_source_hashes'].items():assert R.sha(ROOT/rel)==hh
        else:
            assert audit['profile_sha256']==R.sha(path)
            assert list(map(F,audit['reserves']))==[row['reserve0'],row['reserve6']]
            for report in audit['reports'].values():assert R.sha(ROOT/report['path'])==report['sha256']
        regions=[]
        previous=F(0)
        for filename in box['trees']:
            tree_path=path.parent/filename;tree=R.load(tree_path);p=R.exact(tree['profile'])
            assert tree['status'].startswith('PASS')and tree['pending']==0
            assert F(tree['qa'])==previous;previous=F(tree['qb'])
            assert set(p)==set(row['p'])
            for key in p:
                if key!='major':assert p[key]==row['p'][key]
            for index,reserve in ((0,row['reserve0']),(3,row['reserve6'])):
                assert 0<=p['major'][index]-row['p']['major'][index]<=reserve
            assert p['major'][1:3]==row['p']['major'][1:3]
            regions.append(R.file_record(tree_path))
        assert previous==(F(16)if row['kind']=='full'else F(64))
        bindings[name]=dict(profile=R.file_record(path),accepted_whole_q_audit=R.file_record(ROOT/box['audit']),
            independent_outer_audit=R.file_record(HERE/'targeted_32_outer_vertices.json'),
            q_regions=regions,global_reserves=list(map(str,(row['reserve0'],row['reserve6']))),
            spectral_gates=R.spectral(row,R.B),tail_start=str(previous),
            scope='Accepted fixed-theta whole-q comparison, valid for every actual kappa>=0; the original fitting range constrains only discovery.')
        rows[name]=row
    for item in source['accepted_R3_binding']:
        path=Path(item['path'])
        if not path.is_absolute():path=ROOT/path
        assert R.sha(path)==item['sha256']
    coefficients={}
    for name,row in rows.items():
        c=R.coefficients(row,R.B);row['coefficients']=c
        assert source['profile_coefficients'][name]=={k:str(c[k])for k in ('r_quadratic','r_linear','ell')}
        coefficients[name]={k:str(v)for k,v in c.items()}
    keys=set();intervals={name:set()for name in rows};tables=[]
    for item in source['cached_corner_constants']:
        name=item['profile'];a,b=map(F,item['kappa_interval']);saved=F(item['upper']);row=rows[name]
        assert (name,a,b)not in keys and 0<=a<b<=90
        keys.add((name,a,b));intervals[name].add((a,b));values=[]
        for k,z in product((a,b),(F(2,5),F(1))):
            actual=R.price(row,row['coefficients'],k,z)
            assert actual.below(saved),(name,a,b,k,z)
            values.append(dict(kappa=str(k),shape=str(z),price=actual.record(),
                independent_upper_80_significant_digits=str(F(actual.hi.text(2,80))),
                plus_covariance_correction=str(row['p']['nuisance2']*(k/3-row['coefficients']['theta']*(k+2)))))
        fresh=max(F(v['independent_upper_80_significant_digits'])for v in values)
        assert fresh<=saved
        tables.append(dict(profile=name,kappa_interval=list(map(str,(a,b))),saved_upper=str(saved),
                           fresh_upper=str(fresh),corners=values))
    assert len(keys)==126
    for name,bands in intervals.items():
        previous=F(0);assert len(bands)==14
        for a,b in sorted(bands):assert a==previous;previous=b
        assert previous==90
    assert R.sha(source_path)==digest
    out=dict(status='PASS_INDEPENDENT_NINEBANK_504_R3_CORNER_PRICES',source=R.file_record(source_path),
        operator_b=str(R.B),engine=R.selftest(),profiles=9,corner_tables=126,directed_corner_prices=504,
        new_targeted_tables=28,new_targeted_corner_prices=112,complete_whole_q_bindings=bindings,
        expanded_coefficients=coefficients,corners=tables,
        scope='All nine frozen R3 bank prices and whole-q bindings accepted. Physical/contact/linear exclusion and the complete new envelope algebra require their separate independent audit. Original portable global bundle is unchanged.')
    output_path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ('status','profiles','corner_tables','directed_corner_prices','new_targeted_corner_prices','source')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'tonight/recovery/profile_cover/global_middle_b7171_critical_ninebank_contact_envelope.json')
    parser.add_argument('--output',type=Path,default=HERE/'ninebank_R3_corner_replay.json');args=parser.parse_args();main(args.source.resolve(),args.output.resolve())
