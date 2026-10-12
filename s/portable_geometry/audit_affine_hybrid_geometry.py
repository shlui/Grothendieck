#!/usr/bin/env python3
"""New independent exact affine-envelope replay for the recovered proof.
The archived helper was absent. This is a transparent reconstruction from
its mathematical specification, not a reconstruction of its byte identity.
All arithmetic uses Fraction; all selected native chords and dominance
regions are recomputed, including closed ties and degenerate sets.
"""
from fractions import Fraction as F
import audit_strong_shape_geometry as H


def points(raw):
    return [tuple(map(F,p))for p in raw]


def row_coefficients(row):
    return tuple(F(row[k])for k in ('u_linear','r_linear','r_quadratic','constant'))


def replay(original_raw, proof, rows, native_tables):
    original=points(original_raw)
    assert original and proof['pieces']
    source_rows={r['name']:row_coefficients(r)for r in rows}
    assert len(source_rows)==len(rows)
    natives={('NATIVE:'+r['name']):r for r in native_tables}
    assert len(natives)==len(native_tables)
    minimum=min(r for u,r in original);maximum_r=max(r for u,r in original)
    intervals=[];bounds=[]
    stats=dict(radial_pieces=0,native_piece_bindings=0,
               affine_dominance_regions=0,affine_dominance_vertices=0)
    for piece in proof['pieces']:
        ra,rb=map(F,piece['r_interval']);assert ra<=rb
        assert rb>=minimum and ra<=maximum_r
        poly=H.clip_general(original,F(0),F(-1),-ra)
        poly=H.clip_general(poly,F(0),F(1),rb)
        assert poly and poly==points(piece['polygon'])
        intervals.append((ra,rb));stats['radial_pieces']+=1
        actual_min=min(r for u,r in poly);actual_max=max(r for u,r in poly)
        native_rows={}
        names=set()
        for rec in piece['native_bindings']:
            name=rec['name'];assert name not in names and name in natives;names.add(name)
            table=natives[name];grid=list(map(F,table['grid']));caps=list(map(F,table['caps']))
            assert len(grid)==len(caps)>=2 and all(a<b for a,b in zip(grid,grid[1:]))
            index=rec['piece_index'];assert type(index)is int and 0<=index<len(grid)-1
            a,b=grid[index:index+2];ca,cb=caps[index:index+2]
            assert a<=actual_min<=actual_max<=b
            assert list(map(F,rec['source_interval']))==[a,b]
            assert list(map(F,rec['endpoint_caps']))==[ca,cb]
            assert rec['table_indices']==table['table_indices']
            # A native chord is used only on its actual saved source interval;
            # no interior native knot may be crossed by this polygon piece.
            assert not any(actual_min<k<actual_max for k in grid)
            slope=(cb-ca)/(b-a)
            native_rows[name]=(-2*F(table['ell']),slope,F(table['g']),ca-slope*a)
            stats['native_piece_bindings']+=1
        dp=piece['dominance_proof'];candidates=[];candidate_names=set()
        for rec in dp['candidates']:
            name=rec['name'];assert name not in candidate_names;candidate_names.add(name)
            if name in source_rows:
                au,ar,g,C=source_rows[name]
            else:
                assert name in native_rows;au,ar,g,C=native_rows[name]
            under=rec['underlying'];assert under['name']==name
            assert row_coefficients(under)==(au,ar,g,C)
            ulo,uhi,rlo,rhi=map(F,rec['rectangle']);assert ulo<=uhi and rlo<=rhi
            assert all(ulo<=u<=uhi and rlo<=r<=rhi for u,r in poly)
            Uu=ulo+uhi+au
            if g>=0:
                Ur=ar+g*(rlo+rhi);Uc=C-ulo*uhi-g*rlo*rhi
                assert rec['r_majorant']=='SECANT'
            else:
                center=(rlo+rhi)/2;Ur=ar+2*g*center;Uc=C-ulo*uhi-g*center*center
                assert rec['r_majorant']=='TANGENT'
            assert tuple(F(rec[k])for k in ('u_linear','r_linear','constant'))==(Uu,Ur,Uc)
            candidates.append((name,Uu,Ur,Uc))
        assert candidates
        regions=[];overall=None;attaining=set()
        for name,au,ar,C in candidates:
            domain=list(poly)
            for other,bu,br,D in candidates:
                if other!=name:domain=H.clip_general(domain,au-bu,ar-br,D-C)
            if not domain:continue
            # Every half-plane is closed. Identical affine bounds retain their
            # full common domain, and zero-area domains retain their points.
            values=[au*u+ar*r+C for u,r in domain];upper=max(values)
            regions.append((name,domain,values,upper))
            if overall is None or upper>overall:overall=upper;attaining=set()
            if upper==overall:attaining.update(p for p,v in zip(domain,values)if v==upper)
        assert regions and len(regions)==len(dp['regions'])
        for (name,domain,values,upper),rec in zip(regions,dp['regions']):
            assert rec['profile']==name and points(rec['vertices'])==domain
            assert list(map(F,rec['vertex_values']))==values and F(rec['upper'])==upper
            stats['affine_dominance_regions']+=1
            stats['affine_dominance_vertices']+=len(domain)
        assert tuple(map(F,dp['attaining_point']))in attaining
        assert F(piece['upper'])==overall
        bounds.append(overall)
    # Projection of a convex polygon, segment, or point onto the r-axis is
    # this closed interval. Covering it by the checked clipped pieces proves
    # that every original parameter point belongs to a checked piece.
    cursor=minimum
    for a,b in sorted(intervals):
        assert a<=cursor
        cursor=max(cursor,b)
    assert cursor>=maximum_r
    return max(bounds),stats
