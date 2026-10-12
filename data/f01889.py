#!/usr/bin/env python3
"""Expanded independent prices for accepted new rank-profile source boxes."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import json
import replay_generic_R3_cache_57 as M

R=M.R;HERE=M.HERE;ROOT=M.ROOT
FOLDERS={
    'interval_primary_u7425_k4_10':'primary_u7425_R3_recheck',
    'interval_lowerface_u7175_k053_08':'lowerface_u7175_R3_recheck',
    'interval_nextface_u7375_k4_10':'nextface_u7375_R3_recheck',
    'interval_lowerdiffuse_u72125_k0_01':'lowerdiffuse_u72125_R3_recheck',
}

def main(sources,output):
    boxes=[];inputs=[]
    for source in sources:
        digest=R.sha(source);data=R.load(source)
        assert data['status']=='PASS_DIRECTED_ROBUST_RESTRICTED_R3_SOURCE_BOXES'
        inputs.append(R.file_record(source))
        for box in data['boxes']:
            name=box['name'];folder=FOLDERS[name];assert box['audit']==folder
            for rel,h in box['sha256'].items():assert R.sha(ROOT/rel)==h
            R.B=F(box['operator_b']);assert R.B==F(7171,10000)
            row,binding=M.add_rank(name,folder)
            assert row['profile']==source.parent/box['profile']
            assert row['reserve0']==F(box['reserve0'])and row['reserve6']==F(box['reserve6'])
            assert F(binding['fresh_exact_tail']['R'])==F(box['tail_start'])==32
            assert [x['q']for x in binding['q_regions']]==[x['q']for x in box['regions']]
            for region,saved in zip(binding['q_regions'],box['regions']):
                assert R.load(ROOT/region['proof']['path'])['accepted']==saved['accepted']
            expected_reports=[x['proof']['path']for x in binding['q_regions']]+[binding['accepted_tail']['path']]
            assert expected_reports==box['independent_q_replays']
            c=R.coefficients(row,R.B)
            assert c['r_quadratic']==F(box['r_quadratic'])==R.B>0
            assert c['A']>0 and c['B']>0 and F(194)**2>10*F(633,15)**2
            p=row['p'];E00,E02,E22=R.even(row)
            for r in map(F,box['r']):
                y=p['y']+E02*(r-p['r'])
                assert y-E02*r==c['yc']
                assert (R.B-E00)*r*r-2*p['gamma']*r+y*y/E22==c['r_quadratic']*r*r+c['r_linear']*r+c['yc']**2/E22
            bounds=[tuple(map(F,box[k]))for k in('u','r','kappa','shape')]
            assert all(a<b for a,b in bounds)and bounds[3]==(F(2,5),F(1))
            required=set(product(*bounds));vertices=[];maximum=F(0)
            assert len(box['vertices'])==16
            for saved in box['vertices']:
                point=tuple(F(saved[k])for k in('u','r','kappa','shape'))
                assert point in required;required.remove(point)
                u,r,k,z=point;actual=R.price(row,c,k,z,u,r)
                assert actual.below(F(saved['upper']))
                upper=F(actual.hi.text(2,235));maximum=max(maximum,upper)
                vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),
                    price=actual.record(),saved_upper=saved['upper'],
                    independent_upper=str(upper),positive_covariance_term=str(p['nuisance2']*(k/3-c['theta']*(k+2)))))
            assert not required
            assert maximum<=F(box['maximum_directed_upper'])<F(box['cap'])<F(40,71)
            boxes.append(dict(name=name,source=R.file_record(source),operator_b=str(R.B),
                domain={k:box[k]for k in('u','r','kappa','shape')},
                strict_cap=box['cap'],maximum_independent_upper=str(maximum),
                accepted_whole_q_binding=binding,expanded_coefficients={k:str(v)for k,v in c.items()},
                spectral_gates=R.spectral(row,R.B),vertices=vertices))
        assert R.sha(source)==digest
    result=dict(status='PASS_INDEPENDENT_ADDITIONAL_RESTRICTED_R3_SOURCE_BOXES',
        sources=inputs,engine=M.selftest(),directed_vertex_prices=16*len(boxes),boxes=boxes,
        producer_imports=False,
        convexity='u curvature is one and r curvature is b>0. Each frozen-theta price is separately convex in kappa and shape by A>0, B>0 and 194/sqrt(10)>633/15. Each continuous four-dimensional source box is covered by its sixteen corners.',
        scope='Only these explicitly stated source boxes and accepted whole-q gates. No global amplitude coverage, new global cap, or change to the frozen preferred fourteen-profile package is claimed.')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],directed_vertex_prices=result['directed_vertex_prices'],
        boxes=[dict(name=x['name'],strict_cap=x['strict_cap'],maximum=float(F(x['maximum_independent_upper'])))for x in boxes]),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('sources',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,default=HERE/'additional_64_outer_vertices.json');args=parser.parse_args()
    main([p.resolve()for p in args.sources],args.output.resolve())
