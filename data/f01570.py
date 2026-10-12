#!/usr/bin/env python3
"""Producer-free affine delta replay with an exact target-cap branch.

The original73/75 checker remains unchanged. This separate copy accepts
the legitimate fixed rational40/71 upper when the unrounded piece value
already proves that endpoint, matching the frozen producer v2.
"""
import argparse, json
from fractions import Fraction as F
from pathlib import Path
import audit_strong_shape_geometry as H
from audit_affine_hybrid_geometry import replay


def audit(path, output, parent_report):
    d = H.read(path)
    assert d['status'] == 'COMPLETE_PARENT_LINKED_AFFINE_HYBRID_PENDING_INDEPENDENT_AUDIT'
    pp = Path(d['parent']['path'])
    assert H.sha(pp) == d['parent']['sha256']
    parent = H.read(pp)
    report = H.read(parent_report)
    assert report['status'] == 'PASS_PRODUCER_FREE_COMPLETE_RANK_CHILD_MOMENT_REFINEMENT_GEOMETRY'
    assert report['source']['sha256'] == H.sha(pp)
    assert report['entrypoint_sha256'] == H.sha(Path(__file__).with_name('audit_rank_child_moment_refinement_v2.py'))
    upper_path = Path(d['upper_parent']['path'])
    assert d['upper_parent'] == parent['parent']
    assert H.sha(upper_path) == d['upper_parent']['sha256']
    upper = H.read(upper_path)
    upath = Path(upper['parent']['path'])
    assert H.sha(upath) == upper['parent']['sha256']
    u = H.read(upath)
    bp = Path(u['parent']['path'])
    assert H.sha(bp) == u['parent']['sha256']
    base = H.read(bp)
    b = F(d['operator_b'])
    assert b == F(711,1000) == F(parent['operator_b']) == F(upper['operator_b']) == F(base['operator_b'])
    metadata = {**base['profiles'], **upper['profiles']}
    for binding in d['extra_bindings']:
        p = Path(binding['path'])
        assert H.sha(p) == binding['sha256']
        manifest = H.read(p)
        for record in manifest['profiles']:
            name = record['name']
            assert name not in metadata
            assert H.sha(record['acceptance_path']) == record['acceptance_sha256']
            assert F(record['reserve0']) >= 0 and F(record['reserve6']) >= 0
            for gate in record['gate_parts']:
                assert H.sha(gate['path']) == gate['sha256']
            metadata[name] = record
    coefficients = {n:H.parse_profile(p,b) for n,p in metadata.items()}
    constants = {}
    for cache in (parent['cache'], d['cache']):
        for row in cache:
            ka,kb = map(F,row['kappa_interval'])
            za,zb = map(F,row['shape_interval'])
            key = (row['profile'],ka,kb,za,zb)
            assert key not in constants and key[0] in coefficients
            points = H.physical_moment_points(ka,kb,za,zb)
            assert points and points == [tuple(map(F,p)) for p in row['evaluation_points']]
            values = list(map(F,row['corner_uppers']))
            assert len(points) == len(values) and F(row['upper']) == max(values)
            constants[key] = max(values)
    np = Path(d['native_binding']['path'])
    assert H.sha(np) == d['native_binding']['sha256']
    native = H.read(np)
    assert F(native['operator_b']) == b
    assert native['physical_constraint'] == 'kappa*zeta<=36'
    groups = {}
    for i,row in enumerate(native['tables']):
        grid = list(map(F,row['r_grid']))
        assert len(grid) == len(row['endpoint_uniform_caps']) >= 2
        assert grid[0] == 0 and all(x<y for x,y in zip(grid,grid[1:]))
        groups.setdefault(row['frame'],[]).append((i,row))
    native_cache = {}
    def native_tables(ka,kb):
        if (ka,kb) in native_cache:
            return native_cache[ka,kb]
        tables=[]
        for name,records in groups.items():
            selected=sorted([(i,r) for i,r in records
                if F(r['kappa_interval'][0])<kb and ka<F(r['kappa_interval'][1])],
                key=lambda p:tuple(map(F,p[1]['kappa_interval'])))
            assert selected
            cursor=ka
            for i,r in selected:
                assert F(r['kappa_interval'][0])<=cursor
                cursor=max(cursor,F(r['kappa_interval'][1]))
            assert cursor>=kb
            first=selected[0][1]
            signature=(first['r_grid'],first['ell'],first['r_quadratic'],
                first['source_package'],first['source_sha256'],first['native_mask_sha256'])
            assert all((r['r_grid'],r['ell'],r['r_quadratic'],r['source_package'],
                r['source_sha256'],r['native_mask_sha256'])==signature for i,r in selected)
            grid=list(map(F,first['r_grid']))
            caps=[max(F(r['endpoint_uniform_caps'][j])for i,r in selected)for j in range(len(grid))]
            tables.append(dict(name=name,table_indices=[i for i,r in selected],
                grid=grid,caps=caps,ell=F(first['ell']),g=F(first['r_quadratic'])))
        native_cache[ka,kb]=tables
        return tables
    changes={(r['parent_index'],r['child_index'],r['moment_child_index']):r for r in d['changes']}
    assert len(changes)==len(d['changes'])==d['changed_moment_children']
    selected=failed=changed=0
    maximum=F(40,71);worst=None;seen=set();summaries=[]
    totals=dict(radial_pieces=0,native_piece_bindings=0,
                affine_dominance_regions=0,affine_dominance_vertices=0)
    for group in parent['refinements']:
        for mi,kid in enumerate(group['children']):
            key=(group['parent_index'],group['child_index'],mi)
            if kid['status'].startswith('EXACT_EMPTY'):
                assert key not in changes
                continue
            old=F(kid['upper']);value=old
            if old<=F(40,71):
                assert key not in changes
                continue
            selected+=1
            if key in changes:
                seen.add(key);record=changes[key]
                assert F(record['old_upper'])==old and record['physical_polygon']==kid['physical_polygon']
                ka,kb=map(F,kid['kappa_interval']);za,zb=map(F,kid['shape_interval'])
                rows=[]
                for name,(ell,g,l)in coefficients.items():
                    C=constants[name,ka,kb,za,zb]
                    rows.append(dict(name=name,u_linear=-2*ell,r_linear=l,r_quadratic=g,constant=C))
                exact,stats=replay(kid['physical_polygon'],record['hybrid_proof'],rows,native_tables(ka,kb))
                value=min(old,F(40,71) if exact<=F(40,71) else H.ceil_dec(exact))
                assert value==F(record['new_upper'])<old
                changed+=1
                for name,count in stats.items():totals[name]+=count
            if value>F(40,71):
                failed+=1
                summary=dict(parent_index=key[0],child_index=key[1],moment_child_index=key[2],
                    u_interval=group['u_interval'],r_interval=group['r_interval'],
                    kappa_interval=kid['kappa_interval'],shape_interval=kid['shape_interval'],cap=str(value))
                summaries.append(summary)
                if value>maximum:maximum,worst=value,key
    assert seen==set(changes)
    assert (selected,changed,failed)==(d['selected_moment_children'],d['changed_moment_children'],d['failed_moment_children'])
    assert maximum==F(d['maximum'])
    if failed:
        assert maximum<F(d['strict_upper_cap'])
        raw=d['worst_child'];assert worst==(raw['parent_index'],raw['child_index'],raw['moment_child_index'])
        assert F(raw['cap'])==maximum
        cap_relation='STRICT_UPPER_BOUND'
    else:
        assert maximum==F(40,71)==F(d['strict_upper_cap']) and d['worst_child']is None
        cap_relation='LEQ_EXACT_TARGET_40_OVER_71'
    for saved,own in zip(d['failed_children'],sorted(summaries,key=lambda r:F(r['cap']),reverse=True)):
        assert all(saved[k]==v for k,v in own.items())
    assert len(d['failed_children'])==len(summaries)
    result=dict(status='PASS_PRODUCER_FREE_COMPLETE_AFFINE_HYBRID_DELTA_GEOMETRY',
        source=dict(path=str(Path(path).resolve()),sha256=H.sha(path)),
        parent=dict(path=str(pp),sha256=H.sha(pp),accepted_geometry_report=str(Path(parent_report).resolve()),
                    accepted_geometry_report_sha256=H.sha(parent_report)),operator_b=str(b),
        original_physical_parent_cells=report['original_physical_parent_cells'],
        closed_u_children=report['closed_u_children'],closed_moment_children=report['moment_children'],
        selected_moment_children=selected,changed_moment_children=changed,failed_moment_children=failed,
        profile_count=len(coefficients),fresh_uniform_moment_tables=len(d['cache']),
        maximum=str(maximum),strict_upper_cap=d['strict_upper_cap'],cap_relation=cap_relation,
        producer_imports=False,arithmetic='fractions.Fraction',
        entrypoint_sha256=H.sha(__file__),affine_replay_sha256=H.sha(Path(__file__).with_name('audit_affine_hybrid_geometry.py')),
        scope='Complete unchanged source/moment partition, old-cap retention, native radial knots and affine quadratic-price assembly. Fresh profile gates and directed moment values remain separate independent premises.',**totals)
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--parent-report',required=True)
    a=p.parse_args();audit(a.input,a.output,a.parent_report)
