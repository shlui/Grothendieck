#!/usr/bin/env python3
"""Producer-free closed cumulant/shared-shape descendants of an affine cover.

Every local and directed numerical premise has a separate independent audit.
This checker reconstructs their geometry and the exact all-parent assembly.
No producer module is imported; earlier accepted checkers remain immutable.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import audit_strong_shape_geometry as H
from audit_affine_hybrid_geometry import replay


def audit(path, output, parent_report):
    d=H.read(path)
    assert d['status']=='COMPLETE_PARENT_LINKED_CLOSED_MOMENT_AFFINE_REFINEMENT_PENDING_INDEPENDENT_AUDIT'
    pp=Path(d['parent']['path']);assert H.sha(pp)==d['parent']['sha256'];parent=H.read(pp)
    pr=H.read(parent_report)
    assert pr['status']=='PASS_PRODUCER_FREE_COMPLETE_AFFINE_HYBRID_DELTA_GEOMETRY'
    assert pr['source']['sha256']==H.sha(pp)
    assert pr['entrypoint_sha256']==H.sha(Path(__file__).with_name('audit_affine_hybrid_delta_v2.py'))
    mp=Path(d['moment_ancestor']['path'])
    assert d['moment_ancestor']==parent['parent'] and H.sha(mp)==d['moment_ancestor']['sha256']
    moment=H.read(mp)
    up=Path(moment['parent']['path']);assert H.sha(up)==moment['parent']['sha256'];upper=H.read(up)
    up2=Path(upper['parent']['path']);assert H.sha(up2)==upper['parent']['sha256'];uparent=H.read(up2)
    bp=Path(uparent['parent']['path']);assert H.sha(bp)==uparent['parent']['sha256'];base=H.read(bp)
    b=F(d['operator_b']);assert b==F(711,1000)==F(parent['operator_b'])==F(moment['operator_b'])==F(base['operator_b'])
    metadata={**base['profiles'],**upper['profiles']}
    for binding in parent['extra_bindings']:
        p=Path(binding['path']);assert H.sha(p)==binding['sha256'];manifest=H.read(p)
        for record in manifest['profiles']:
            name=record['name'];assert name not in metadata
            assert H.sha(record['acceptance_path'])==record['acceptance_sha256']
            assert F(record['reserve0'])>=0 and F(record['reserve6'])>=0
            for gate in record['gate_parts']:assert H.sha(gate['path'])==gate['sha256']
            metadata[name]=record
    coefficients={n:H.parse_profile(p,b)for n,p in metadata.items()}
    constants={}
    for row in d['cache']:
        ka,kb=map(F,row['kappa_interval']);za,zb=map(F,row['shape_interval'])
        key=(row['profile'],ka,kb,za,zb)
        assert key not in constants and key[0]in coefficients
        points=H.physical_moment_points(ka,kb,za,zb)
        assert points and points==[tuple(map(F,p))for p in row['evaluation_points']]
        values=list(map(F,row['corner_uppers']))
        assert len(values)==len(points) and F(row['upper'])==max(values)
        constants[key]=max(values)
    bankp=Path(d['source_plane_bank']['path'])
    assert H.sha(bankp)==d['source_plane_bank']['sha256'];bank=H.read(bankp)
    planes={}
    for row in d['plane_cache']:
        ka,kb=map(F,row['kappa_interval']);za,zb=map(F,row['shape_interval']);key=(ka,kb,za,zb)
        assert key not in planes and row['physical_constraint']=='4*kappa*zeta-kappa<=54'
        points=H.physical_moment_points(*key)
        assert points and points==[tuple(map(F,p))for p in row['evaluation_points']]
        rows=[]
        for plane in row['planes']:
            raw=bank['rows'][plane['source_row']]
            assert F(raw['a'])==F(plane['a'])>0
            assert F(raw['polynomial_coefficients_q0_q2_q4_q6'][3])>=0
            values=list(map(F,plane['corner_uppers']))
            assert len(values)==len(points) and F(plane['cap'])==max(values)
            rows.append((F(plane['a']),F(plane['cap'])))
        assert rows;planes[key]=rows
    assert d['native_binding']==parent['native_binding']
    np=Path(d['native_binding']['path']);assert H.sha(np)==d['native_binding']['sha256'];native=H.read(np)
    assert F(native['operator_b'])==b and native['physical_constraint']=='kappa*zeta<=36'
    groups={}
    for i,row in enumerate(native['tables']):
        grid=list(map(F,row['r_grid']))
        assert len(grid)==len(row['endpoint_uniform_caps'])>=2
        assert grid[0]==0 and all(x<y for x,y in zip(grid,grid[1:]))
        groups.setdefault(row['frame'],[]).append((i,row))
    native_cache={}
    def native_tables(ka,kb):
        if (ka,kb)in native_cache:return native_cache[ka,kb]
        tables=[]
        for name,records in groups.items():
            selected=sorted([(i,r)for i,r in records
                if F(r['kappa_interval'][0])<kb and ka<F(r['kappa_interval'][1])],
                key=lambda p:tuple(map(F,p[1]['kappa_interval'])))
            assert selected;cursor=ka
            for i,r in selected:
                assert F(r['kappa_interval'][0])<=cursor
                cursor=max(cursor,F(r['kappa_interval'][1]))
            assert cursor>=kb;first=selected[0][1]
            signature=(first['r_grid'],first['ell'],first['r_quadratic'],
                       first['source_package'],first['source_sha256'],first['native_mask_sha256'])
            assert all((r['r_grid'],r['ell'],r['r_quadratic'],r['source_package'],
                        r['source_sha256'],r['native_mask_sha256'])==signature for i,r in selected)
            grid=list(map(F,first['r_grid']))
            caps=[max(F(r['endpoint_uniform_caps'][j])for i,r in selected)for j in range(len(grid))]
            tables.append(dict(name=name,table_indices=[i for i,r in selected],grid=grid,caps=caps,
                               ell=F(first['ell']),g=F(first['r_quadratic'])))
        native_cache[ka,kb]=tables;return tables
    parts=d['parts'];assert type(parts)is int and parts>0
    threshold=F(40,71)
    ancestor_groups={(r['parent_index'],r['child_index']):r for r in moment['refinements']}
    failed={(r['parent_index'],r['child_index'],r['moment_child_index']):r for r in parent['failed_children']}
    assert len(failed)==parent['failed_moment_children']==pr['failed_moment_children']
    refinements={(r['parent_index'],r['child_index'],r['moment_child_index']):r for r in d['refinements']}
    assert len(refinements)==len(d['refinements'])==d['selected_moment_children']
    assert set(refinements)<=set(failed)
    replacements={};child_count=empty=improved=0
    totals=dict(radial_pieces=0,native_piece_bindings=0,affine_dominance_regions=0,affine_dominance_vertices=0)
    needed_boxes=set();needed_constants=set()
    for key,record in refinements.items():
        pi,ci,mi=key;group=ancestor_groups[pi,ci];kid=group['children'][mi]
        assert not kid['status'].startswith('EXACT_EMPTY')
        old=F(failed[key]['cap']);assert old>threshold and F(record['old_upper'])==old
        assert record['u_interval']==group['u_interval'] and record['r_interval']==group['r_interval']
        assert record['kappa_interval']==kid['kappa_interval'] and record['shape_interval']==kid['shape_interval']
        assert record['physical_polygon']==kid['physical_polygon']
        original=[tuple(map(F,p))for p in kid['physical_polygon']];assert original
        ka,kb=map(F,kid['kappa_interval']);za,zb=map(F,kid['shape_interval']);assert ka<kb and za<zb
        ks=[ka+(kb-ka)*i/parts for i in range(parts+1)]
        zs=[za+(zb-za)*i/parts for i in range(parts+1)]
        expected=[(a,c,z,e)for a,c in zip(ks,ks[1:])for z,e in zip(zs,zs[1:])]
        assert len(record['children'])==len(expected)
        cap=F(0)
        for child,box in zip(record['children'],expected):
            child_count+=1;a,c,z,e=box
            assert list(map(F,child['kappa_interval']))==[a,c]
            assert list(map(F,child['shape_interval']))==[z,e]
            assert F(child['old_upper'])==old
            points=H.physical_moment_points(*box)
            if not points:
                assert child['status']=='EXACT_EMPTY_NECESSARY_MOMENT_CHILD' and child['physical_polygon']==[]
                empty+=1;continue
            needed_boxes.add(box);poly=list(original)
            for aa,cc in planes[box]:poly=H.clip(poly,aa,cc)
            assert poly==[tuple(map(F,p))for p in child['physical_polygon']]
            if not poly:
                assert child['status']=='EXACT_EMPTY_SOURCE_CHILD';empty+=1;continue
            assert child['status']=='EXACT_DIRECTED_AFFINE_MOMENT_CHILD_PENDING_AUDIT'
            rows=[]
            for name,(ell,g,l)in coefficients.items():
                ck=(name,*box);needed_constants.add(ck)
                rows.append(dict(name=name,u_linear=-2*ell,r_linear=l,r_quadratic=g,constant=constants[ck]))
            exact,stats=replay(child['physical_polygon'],child['hybrid_proof'],rows,native_tables(a,c))
            assert exact==F(child['raw_affine_upper'])
            value=min(old,threshold if exact<=threshold else H.ceil_dec(exact))
            assert value==F(child['upper']);cap=max(cap,value)
            for name,count in stats.items():totals[name]+=count
        assert cap==F(record['new_parent_upper']) and cap<=old
        replacements[key]=cap
        if cap<old:improved+=1
    assert child_count==len(refinements)*parts*parts
    assert set(planes)==needed_boxes and set(constants)==needed_constants
    remaining=[];maximum=threshold
    for key,summary in failed.items():
        cap=replacements.get(key,F(summary['cap']))
        if cap>threshold:
            remaining.append({**summary,'cap':str(cap)})
            maximum=max(maximum,cap)
    remaining.sort(key=lambda q:F(q['cap']),reverse=True)
    assert improved==d['improved_moment_children'] and len(remaining)==d['failed_moment_children']
    assert maximum==F(d['maximum'])
    assert len(remaining)==len(d['failed_children'])
    for own,saved in zip(remaining,d['failed_children']):
        assert all(saved[k]==v for k,v in own.items()if k!='cap_display')
    if remaining:
        strict=F(d['strict_upper_cap'])
        expected_cap=min(F(parent['strict_upper_cap']),H.ceil_dec(maximum,9)+F(1,10**9))
        assert maximum<strict==expected_cap and d['upper_cap_relation']=='STRICT_LT'
        assert all(d['worst_child'][k]==v for k,v in remaining[0].items()if k!='cap_display')
        relation='STRICT_UPPER_BOUND'
    else:
        assert maximum==threshold==F(d['strict_upper_cap']) and d['upper_cap_relation']=='LEQ' and d['worst_child']is None
        relation='LEQ_EXACT_TARGET_40_OVER_71'
    result=dict(status='PASS_PRODUCER_FREE_COMPLETE_CLOSED_MOMENT_AFFINE_REFINEMENT_GEOMETRY',
        source=dict(path=str(Path(path).resolve()),sha256=H.sha(path)),
        parent=dict(path=str(pp),sha256=H.sha(pp),accepted_geometry_report=str(Path(parent_report).resolve()),
                    accepted_geometry_report_sha256=H.sha(parent_report)),operator_b=str(b),
        original_physical_parent_cells=pr['original_physical_parent_cells'],closed_u_children=pr['closed_u_children'],
        inherited_closed_moment_children=pr['closed_moment_children'],selected_moment_children=len(refinements),
        added_closed_moment_children=child_count,exact_empty_children=empty,improved_moment_children=improved,
        failed_moment_children=len(remaining),profile_count=len(coefficients),fresh_uniform_moment_tables=len(constants),
        fresh_uniform_plane_boxes=len(planes),maximum=str(maximum),strict_upper_cap=d['strict_upper_cap'],
        cap_relation=relation,producer_imports=False,arithmetic='fractions.Fraction',
        entrypoint_sha256=H.sha(__file__),affine_replay_sha256=H.sha(Path(__file__).with_name('audit_affine_hybrid_geometry.py')),
        scope='Complete closed cumulant/shared-shape descendant partition and affine source geometry; all inherited caps retained. Directed prices/local gates remain separately audited premises.',**totals)
    Path(output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--parent-report',required=True)
    a=p.parse_args();audit(a.input,a.output,a.parent_report)
