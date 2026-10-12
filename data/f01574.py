#!/usr/bin/env python3
"""Producer-free shared moment refinement, with exact <=40/71 case.

The original49 checker remains frozen. This copy accepts the separately
audited monotone v2 parent and an exact unrounded target endpoint.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import audit_strong_shape_geometry as H


def main(path,out,parent_report):
    d=H.read(path);assert d['status']=='COMPLETE_PARENT_LINKED_KAPPA_SHAPE_RANK_REFINEMENT_PENDING_INDEPENDENT_AUDIT'
    pp=Path(d['parent']['path']);assert H.sha(pp)==d['parent']['sha256'];parent=H.read(pp)
    pr=H.read(parent_report);assert pr['status']=='PASS_PRODUCER_FREE_COMPLETE_MONOTONE_RANK_CHILD_GEOMETRY'
    allowed_parent_entries={H.sha(Path(__file__).with_name(n))for n in ('audit_monotone_rank_children.py','audit_monotone_rank_children_v2.py')}
    assert pr['source']['sha256']==H.sha(pp)and pr['entrypoint_sha256']in allowed_parent_entries
    up=Path(d['u_parent']['path']);assert H.sha(up)==d['u_parent']['sha256'];u=H.read(up)
    assert d['u_parent']==parent['parent'];bp=Path(d['base_global']['path']);assert H.sha(bp)==d['base_global']['sha256'];base=H.read(bp)
    assert d['base_global']==u['parent']==parent['base_global'];b=F(d['operator_b']);assert b==F(parent['operator_b'])==F(u['operator_b'])==F(base['operator_b'])==F(711,1000)
    coefficients={n:H.parse_profile(r,b)for n,r in {**base['profiles'],**parent['profiles']}.items()}
    constants={}
    for r in d['cache']:
        ka,kb=map(F,r['kappa_interval']);za,zb=map(F,r['shape_interval']);key=(r['profile'],ka,kb,za,zb)
        assert key not in constants and key[0]in coefficients
        points=H.physical_moment_points(ka,kb,za,zb)
        assert points and points==[tuple(map(F,p))for p in r['evaluation_points']]
        prices=list(map(F,r['corner_uppers']));assert len(prices)==len(points)and F(r['upper'])==max(prices);constants[key]=max(prices)
    bankp=Path(d['source_plane_bank']['path']);assert H.sha(bankp)==d['source_plane_bank']['sha256'];bank=H.read(bankp)
    planes={}
    for r in d['plane_cache']:
        ka,kb=map(F,r['kappa_interval']);za,zb=map(F,r['shape_interval']);key=(ka,kb,za,zb)
        assert key not in planes and r['physical_constraint']=='4*kappa*zeta-kappa<=54'
        points=H.physical_moment_points(ka,kb,za,zb)
        assert points and points==[tuple(map(F,p))for p in r['evaluation_points']];rows=[]
        for p in r['planes']:
            raw=bank['rows'][p['source_row']];assert F(raw['a'])==F(p['a'])>0 and F(raw['polynomial_coefficients_q0_q2_q4_q6'][3])>=0
            prices=list(map(F,p['corner_uppers']));assert len(prices)==len(points)and F(p['cap'])==max(prices);rows.append((F(p['a']),F(p['cap'])))
        assert rows;planes[key]=rows
    upgrades={(r['parent_index'],r['child_index']):F(r['new_offer']['upper'])for r in parent['changes']}
    refinements={(r['parent_index'],r['child_index']):r for r in d['refinements']}
    assert len(refinements)==len(d['refinements'])==d['selected_u_children']
    kp,zp=d['kappa_parts'],d['shape_parts'];assert kp>0 and zp>0
    threshold=F(40,71);maximum=threshold;worst=None;empty=changed=failed=children=rankoffers=regions=vertices=0;seen=set()
    for p in u['refinements']:
        i=p['parent_index'];ka,kb=map(F,p['kappa_interval']);assert ka<kb
        knots=[ka+(kb-ka)*j/kp for j in range(kp+1)];shapes=[F(2,5)+F(3,5)*j/zp for j in range(zp+1)]
        expected=[(a,z,c,e)for a,z in zip(knots,knots[1:])for c,e in zip(shapes,shapes[1:])]
        for j,oldchild in enumerate(p['children']):
            key=(i,j)
            if oldchild['status'].startswith('EXACT_EMPTY'):assert key not in refinements;continue
            old=upgrades.get(key,F(oldchild['upper']));assert(key in refinements)==(old>threshold)
            if key not in refinements:continue
            seen.add(key);r=refinements[key];assert F(r['old_upper'])==old
            assert r['u_interval']==oldchild['u_interval']and r['r_interval']==p['r_interval']and r['kappa_interval']==p['kappa_interval']
            assert list(map(F,r['shape_interval']))==[F(2,5),F(1)]and len(r['children'])==len(expected)
            original=[tuple(map(F,x))for x in oldchild['physical_polygon']];assert original;cap=F(0)
            for child,(a,z,c,e)in zip(r['children'],expected):
                children+=1;box=(a,z,c,e)
                assert list(map(F,child['kappa_interval']))==[a,z]and list(map(F,child['shape_interval']))==[c,e]
                points=H.physical_moment_points(a,z,c,e)
                if not points:
                    assert child['status']=='EXACT_EMPTY_PHYSICAL_MOMENT_SUBCELL'and child['physical_polygon']==[];empty+=1;continue
                poly=list(original)
                for aa,cc in planes[box]:poly=H.clip(poly,aa,cc)
                assert poly==[tuple(map(F,x))for x in child['physical_polygon']]
                if not poly:assert child['status']=='EXACT_EMPTY_SOURCE_SUBCELL';empty+=1;continue
                assert child['status']=='EXACT_DIRECTED_RANK_MOMENT_SUBCELL_PENDING_AUDIT'and F(child['old_upper'])==old
                value=old;offer=child['new_offer']
                if offer:
                    assert offer['profile']=='POINTWISE_RANK_ENVELOPE';proof=offer['dominance_proof'];candidates=[];sig=set()
                    for q in proof['candidates']:
                        n=q['name'];ell,g,l=coefficients[n];C=constants[(n,*box)]
                        assert g==b and F(q['u_linear'])==-2*ell and F(q['r_linear'])==l and F(q['constant'])==C
                        ss=(-2*ell,l,C);assert ss not in sig;sig.add(ss);candidates.append((n,-2*ell,l,C))
                    assert candidates;own=[];overall=None;att=set()
                    for n,au,ar,C in candidates:
                        part=list(poly)
                        for other,bu,br,D in candidates:
                            if other!=n:part=H.clip_general(part,au-bu,ar-br,D-C)
                        if not part:continue
                        vv=[u*u+b*r*r+au*u+ar*r+C for u,r in part];v=max(vv);own.append((n,part,vv,v))
                        if overall is None or v>overall:overall=v;att=set()
                        if v==overall:att.update(q for q,w in zip(part,vv)if w==v)
                    assert own and len(own)==len(proof['regions'])
                    for(n,part,vv,v),q in zip(own,proof['regions']):
                        assert n==q['profile']and part==[tuple(map(F,x))for x in q['vertices']]
                        assert vv==list(map(F,q['vertex_values']))and v==F(q['upper']);regions+=1;vertices+=len(part)
                    assert tuple(map(F,proof['attaining_point']))in att;value=H.ceil_dec(overall)
                    if overall<=threshold and value>threshold:value=threshold
                    assert value==F(offer['upper'])<old;rankoffers+=1
                assert value==F(child['upper']);cap=max(cap,value)
                if value>threshold:
                    failed+=1
                    if value>maximum:maximum=value;worst=(i,j,box)
            assert cap==F(r['new_parent_upper'])and cap<=old
            if cap<old:changed+=1
    assert seen==set(refinements)and children==len(refinements)*kp*zp
    assert empty==d['empty_moment_children']and changed==d['improved_u_children']and failed==d['failed_moment_children']
    assert maximum==F(d['maximum'])
    if failed:
        assert maximum<F(d['strict_upper_cap'])
        w=d['worst_child'];assert worst==(w['parent_index'],w['child_index'],(*map(F,w['kappa_interval']),*map(F,w['shape_interval'])))
        cap_relation='STRICT_UPPER_BOUND'
    else:
        assert maximum==threshold==F(d['strict_upper_cap'])and worst is None and d['worst_child']is None
        cap_relation='LEQ_EXACT_TARGET_40_OVER_71'
    result=dict(status='PASS_PRODUCER_FREE_COMPLETE_RANK_CHILD_MOMENT_REFINEMENT_GEOMETRY',source=dict(path=str(Path(path).resolve()),sha256=H.sha(path)),parent=dict(path=str(pp),sha256=H.sha(pp),accepted_geometry_report=str(Path(parent_report).resolve()),accepted_geometry_report_sha256=H.sha(parent_report)),operator_b=str(b),original_physical_parent_cells=u['cells'],closed_u_children=sum(len(p['children'])for p in u['refinements']),selected_u_children=len(refinements),moment_children=children,exact_empty_moment_children=empty,improved_u_children=changed,rank_offers=rankoffers,rank_dominance_regions=regions,rank_dominance_vertices=vertices,fresh_uniform_moment_tables=len(constants),fresh_uniform_plane_boxes=len(planes),maximum=str(maximum),strict_upper_cap=d['strict_upper_cap'],failed_moment_children=failed,cap_relation=cap_relation,producer_imports=False,arithmetic='fractions.Fraction',entrypoint_sha256=H.sha(__file__),scope='Exact complete parent-linked shared cumulant/shape partition, source clipping and rank-price assembly. Fresh moment/source numeric caches remain independent premises.')
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--parent-report',required=True);a=p.parse_args();main(a.input,a.output,a.parent_report)
