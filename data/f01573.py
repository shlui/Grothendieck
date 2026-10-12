#!/usr/bin/env python3
"""Producer-free rank additions; supports exact <=40/71 completion.

The original frozen checker remains unchanged. This copy adds only the
zero-failure endpoint case, retaining every exact geometry computation.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import audit_strong_shape_geometry as H


def main(path,out,parent_report):
    d=H.read(path);assert d['status']=='COMPLETE_MONOTONE_RANK_CHILD_UPGRADE_PENDING_INDEPENDENT_AUDIT'
    pp=Path(d['parent']['path']);assert H.sha(pp)==d['parent']['sha256'];parent=H.read(pp)
    pr=H.read(parent_report);assert pr['status']=='PASS_PRODUCER_FREE_COMPLETE_STRONG_SHAPE_U_REFINEMENT_GEOMETRY'
    assert pr['source']['sha256']==H.sha(pp)
    assert pr['entrypoint_sha256']==H.sha(Path(__file__).with_name('audit_u_refinement_geometry.py'))
    bp=Path(d['base_global']['path']);assert H.sha(bp)==d['base_global']['sha256'];base=H.read(bp)
    assert d['base_global']==parent['parent'];b=F(d['operator_b']);assert b==F(parent['operator_b'])==F(base['operator_b'])==F(711,1000)
    assert set(base['profiles']).isdisjoint(d['profiles'])
    coefficients={n:H.parse_profile(r,b)for n,r in {**base['profiles'],**d['profiles']}.items()}
    for n in d['profiles']:assert coefficients[n][1]==b
    constants={(r['profile'],tuple(map(F,r['kappa_interval']))):F(r['upper'])for r in base['cache']}
    freshkeys=set()
    for r in d['cache']:
        n=r['profile'];assert n in d['profiles'];ka,kb=map(F,r['kappa_interval']);za,zb=map(F,r['shape_interval']);assert[za,zb]==[F(2,5),F(1)]
        points=H.physical_moment_points(ka,kb,za,zb)
        assert points and points==[tuple(map(F,p))for p in r['evaluation_points']]
        prices=list(map(F,r['corner_uppers']));assert len(prices)==len(points)and F(r['upper'])==max(prices)
        key=(n,(ka,kb));assert key not in constants;constants[key]=max(prices);freshkeys.add(key)
    changes={(r['parent_index'],r['child_index']):r for r in d['changes']}
    assert len(changes)==len(d['changes'])==d['changed_children']
    threshold=F(40,71);maximum=threshold;worst=None;seen=set();regions=vertices=failed=nonempty=0
    for p in parent['refinements']:
        i=p['parent_index'];kband=tuple(map(F,p['kappa_interval']))
        for j,c in enumerate(p['children']):
            key=(i,j)
            if c['status'].startswith('EXACT_EMPTY'):assert key not in changes;continue
            nonempty+=1;old=F(c['upper']);cap=old
            if key in changes:
                seen.add(key);change=changes[key];assert old>threshold and F(change['old_upper'])==old
                poly=[tuple(map(F,x))for x in c['physical_polygon']];assert poly
                offer=change['new_offer'];assert offer['profile']=='POINTWISE_RANK_ENVELOPE';proof=offer['dominance_proof']
                candidates=[];signatures=set()
                for r in proof['candidates']:
                    n=r['name'];ell,g,l=coefficients[n];C=constants[(n,kband)]
                    assert g==b and F(r['u_linear'])==-2*ell and F(r['r_linear'])==l and F(r['constant'])==C
                    signature=(-2*ell,l,C);assert signature not in signatures;signatures.add(signature);candidates.append((n,-2*ell,l,C))
                assert candidates;own=[];overall=None;att=set()
                for n,au,ar,C in candidates:
                    part=list(poly)
                    for other,bu,br,D in candidates:
                        if other!=n:part=H.clip_general(part,au-bu,ar-br,D-C)
                    if not part:continue
                    values=[u*u+b*r*r+au*u+ar*r+C for u,r in part];v=max(values);own.append((n,part,values,v))
                    if overall is None or v>overall:overall=v;att=set()
                    if v==overall:att.update(p for p,vv in zip(part,values)if vv==v)
                assert own and len(own)==len(proof['regions'])
                for(n,part,values,v),r in zip(own,proof['regions']):
                    assert n==r['profile']and part==[tuple(map(F,x))for x in r['vertices']]
                    assert values==list(map(F,r['vertex_values']))and F(r['upper'])==v;regions+=1;vertices+=len(part)
                assert tuple(map(F,proof['attaining_point']))in att
                cap=H.ceil_dec(overall);assert cap==F(offer['upper'])<old
            if cap>threshold:
                failed+=1
                if cap>maximum:maximum=cap;worst=key
    assert seen==set(changes)and failed==len(d['failed_children'])
    assert maximum==F(d['maximum'])
    if failed:
        assert maximum<F(d['strict_upper_cap'])
        assert worst==(d['worst_child']['parent_index'],d['worst_child']['child_index'])and maximum==F(d['worst_child']['cap'])
        cap_relation='STRICT_UPPER_BOUND'
    else:
        assert maximum==threshold==F(d['strict_upper_cap'])
        assert worst is None and d['worst_child'] is None
        cap_relation='LEQ_EXACT_TARGET_40_OVER_71'
    result=dict(status='PASS_PRODUCER_FREE_COMPLETE_MONOTONE_RANK_CHILD_GEOMETRY',source=dict(path=str(Path(path).resolve()),sha256=H.sha(path)),parent=dict(path=str(pp),sha256=H.sha(pp),accepted_geometry_report=str(Path(parent_report).resolve()),accepted_geometry_report_sha256=H.sha(parent_report)),operator_b=str(b),original_parent_cells=parent['cells'],closed_u_children=sum(len(r['children'])for r in parent['refinements']),nonempty_u_children=nonempty,changed_children=len(changes),fresh_profiles=len(d['profiles']),fresh_uniform_moment_tables=len(freshkeys),rank_dominance_regions=regions,rank_dominance_vertices=vertices,maximum=str(maximum),strict_upper_cap=d['strict_upper_cap'],failed_children=failed,cap_relation=cap_relation,producer_imports=False,arithmetic='fractions.Fraction',entrypoint_sha256=H.sha(__file__),scope='Exact unchanged-parent partition and rank-offer assembly only. All new whole-q and directed price inputs still require independent acceptance before norm promotion.')
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--parent-report',required=True);a=p.parse_args();main(a.input,a.output,a.parent_report)
