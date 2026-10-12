#!/usr/bin/env python3
"""Separate targeted upgrade audit; never changes frozen original envelope."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'tonight/operator/projection/targeted_certificate'
SOURCE=BASE/'targeted_source_boxes.json'
spec=importlib.util.spec_from_file_location('independent_R3_expansion',HERE/'replay_R3_corners.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)


def main():
    digest=R.sha(SOURCE);source=R.load(SOURCE)
    assert source['status']=='PASS_DIRECTED_CONTINUOUS_TARGETED_R3_SOURCE_BOXES'
    results=[]
    for box in source['boxes']:
        for rel,digest0 in box['sha256'].items():assert R.sha(ROOT/rel)==digest0
        path=BASE/box['profile'];row=dict(p=R.exact(R.load(path)['profile']),kind=box['kind'],
                                       reserve0=F(box['reserve0']),reserve6=F(box['reserve6']))
        b=F(box['operator_b']);assert b==F(7171,10000)
        coefficients=R.coefficients(row,b);spectral=R.spectral(row,b)
        assert coefficients['r_quadratic']==F(box['r_quadratic'])>=0
        audit=R.load(ROOT/box['audit']);assert audit['status'].startswith('PASS_INDEPENDENT')
        if box['kind']=='full':
            assert F(audit['operator_b'])==b
            assert F(audit['global_reserve0'])==row['reserve0']
            assert F(audit['global_reserve6'])==row['reserve6']
            for rel,hh in audit['frozen_source_hashes'].items():assert R.sha(ROOT/rel)==hh
        else:
            assert audit['profile_sha256']==R.sha(path)
            assert list(map(F,audit['reserves']))==[row['reserve0'],row['reserve6']]
            for report in audit['reports'].values():assert R.sha(ROOT/report['path'])==report['sha256']
        bounds=[tuple(map(F,box[key]))for key in ('u','r','kappa','shape')]
        assert bounds[3]==(F(2,5),F(1))
        vertices=[];assert len(box['vertices'])==16
        for saved,(u,r,k,z)in zip(box['vertices'],product(*bounds)):
            assert tuple(F(saved[key])for key in ('u','r','kappa','shape'))==(u,r,k,z)
            value=R.price(row,coefficients,k,z,u,r)
            assert value.below(F(saved['upper']))and value.below(F(box['cap']))
            rounded=F(value.hi.text(2,80))
            vertices.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),price=value.record(),
                independent_upper_80_digits=str(rounded),saved_upper=saved['upper'],
                plus_covariance_term=str(row['p']['nuisance2']*(k/3-coefficients['theta']*(k+2)))))
        maximum=max(F(v['saved_upper'])for v in vertices)
        assert maximum==F(box['maximum_directed_upper'])<F(box['cap'])<F(40,71)
        results.append(dict(name=box['name'],profile=R.file_record(path),whole_q_acceptance=R.file_record(ROOT/box['audit']),
            b=str(b),u=box['u'],r=box['r'],kappa=box['kappa'],shape=box['shape'],
            reserves=list(map(str,(row['reserve0'],row['reserve6']))),
            frozen_theta=str(coefficients['theta']),exact_spectral_gates=spectral,
            expanded_coefficients={key:str(v)for key,v in coefficients.items()},
            cap=box['cap'],maximum_saved_upper=str(maximum),vertices=vertices))
    assert len(results)==2 and R.sha(SOURCE)==digest
    out=dict(status='PASS_INDEPENDENT_TARGETED_32_OUTER_SOURCE_VERTICES',
             source=R.file_record(SOURCE),engine=R.selftest(),boxes=results,
             total_vertices=32,producer_imports=False,
             convexity='u and r convex (r curvature nonnegative); kappa/shape separately convex by A>0,B>0,H(shape)>0. Thus16 endpoints cover each continuous4Dbox.',
             scope='Only these two restricted source boxes, using their separately accepted whole-q gates. The frozen original global cap is unchanged; stronger global claims require a new exhaustive envelope audit.')
    target=HERE/'targeted_32_outer_vertices.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'total_vertices':32,'boxes':[{'name':x['name'],'cap':x['cap'],'maximum_saved_upper':x['maximum_saved_upper']}for x in results]},indent=2))


if __name__=='__main__':main()
