#!/usr/bin/env python3
"""Fresh expanded outer prices for two independently accepted robust R3 gates."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import os
import sys
import argparse

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest
import replay_R3_corners as R
from replay_middlebank_R3_corners import fresh_tail

SOURCE=ROOT/'tonight/operator/projection/targeted_certificate/robust_source_boxes.json'


def main(source_path=SOURCE,output_path=HERE/'robust_32_outer_vertices.json',threshold=F(40,71)):
    global SOURCE
    SOURCE=source_path
    digest=R.sha(SOURCE);data=R.load(SOURCE)
    assert data['status']=='PASS_DIRECTED_ROBUST_RESTRICTED_R3_SOURCE_BOXES'
    boxes=[];total=0
    for box in data['boxes']:
        for rel,h in box['sha256'].items():assert R.sha(ROOT/rel)==h
        path=SOURCE.parent/box['profile'];p=R.exact(R.load(path)['profile'])
        row=dict(name=box['name'],profile=path,p=p,kind='rank',reserve0=F(box['reserve0']),reserve6=F(box['reserve6']))
        assert F(box['operator_b'])==R.B and F(box['tail_start'])==20
        reports=[];previous=F(0)
        for part,reportrel in zip(('body','exterior'),box['independent_q_replays'][:2]):
            reportpath=ROOT/reportrel;report=R.load(reportpath)
            assert report['status']=='PASS_INDEPENDENT_PHYSICAL_Q_TREE'
            assert report['profile_sha256']==R.sha(path)
            treepath=SOURCE.parent/f'{part}_{box["name"]}.json.gz';tree=R.load(treepath)
            assert R.sha(treepath)==report['tree_sha256']
            assert tree['status'].startswith('PASS')and tree['pending']==0 and tree['accepted']==report['accepted']
            assert F(tree['qa'])==previous;previous=F(tree['qb']);inflated=R.exact(tree['profile'])
            for k in p:
                if k!='major':assert inflated[k]==p[k]
            assert inflated['major']==[v+(row['reserve0']if j==0 else row['reserve6']if j==3 else F(0))for j,v in enumerate(p['major'])]
            assert report['cubic_integral_backend_sha256']==R.sha(ROOT/'tonight/recovery/profile_cover/directed_cubic_integral.py')
            reports.append(dict(tree=R.file_record(treepath),accepted_proof=R.file_record(reportpath),accepted=report['accepted']))
        assert previous==20
        tailpath=ROOT/box['independent_q_replays'][2];tail=R.load(tailpath)
        assert tail['status']=='PASS_INDEPENDENT_EXACT_R3_TAIL'
        fresh=fresh_tail(row,F(20))
        for k in('R','difference_polynomial','shifted_coefficients'):assert fresh[k]==tail[k]
        c=R.coefficients(row,R.B);assert c['r_quadratic']==R.B==F(box['r_quadratic'])>0
        gates=R.spectral(row,R.B)
        ends=[tuple(map(F,box[k]))for k in('u','r','kappa','shape')]
        assert ends[3]==(F(2,5),F(1))
        todo=set(product(*ends));vertices=[];largest=F(0)
        for item in box['vertices']:
            point=tuple(F(item[k])for k in('u','r','kappa','shape'));assert point in todo;todo.remove(point)
            u,r,k,z=point;actual=R.price(row,c,k,z,u,r);saved=F(item['upper'])
            assert actual.below(saved)
            independent=F(actual.hi.text(2,235));largest=max(largest,independent)
            vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),price=actual.record(),saved_upper=str(saved),
                plus_covariance=str(p['nuisance2']*(k/3-c['theta']*(k+2)))))
        assert not todo and len(vertices)==16
        cap=F(box['cap']);assert largest<=F(box['maximum_directed_upper'])<cap<threshold
        total+=16
        boxes.append(dict(name=box['name'],profile=R.file_record(path),reserves=[str(row['reserve0']),str(row['reserve6'])],
            scope={k:box[k]for k in('u','r','kappa','shape')},cap=str(cap),maximum_independent_upper=str(largest),
            frozen_theta=str(c['theta']),expanded_coefficients={k:str(v)for k,v in c.items()},spectral_gates=gates,
            accepted_q_regions=reports,accepted_tail_proof=R.file_record(tailpath),fresh_exact_tail=fresh,vertices=vertices))
    assert total==16*len(data['boxes']) and R.sha(SOURCE)==digest
    out=dict(status=f'PASS_INDEPENDENT_ROBUST_{total}_OUTER_SOURCE_VERTICES',source=R.file_record(SOURCE),engine=selftest(),
        strict_threshold=str(threshold),
        operator_b=str(R.B),directed_outer_prices=total,boxes=boxes,
        convexity='The r coefficient is positive, the u coefficient is one, and the kappa/shape price is separately convex with A>0 and B>0. The continuous box is bounded by its sixteen corners. Frozen-theta covariance enters with a positive sign.',
        scope='Only the two stated continuous source boxes and their accepted whole-q prerequisites. No complete amplitude coverage or final global norm bound is claimed. Existing frozen bundles unchanged.')
    output_path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],directed_outer_prices=total,boxes=[dict(name=x['name'],cap=x['cap'])for x in boxes]),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--output',type=Path,default=HERE/'robust_32_outer_vertices.json')
    parser.add_argument('--threshold',default='40/71');args=parser.parse_args()
    main(args.source.resolve(),args.output.resolve(),F(args.threshold))
