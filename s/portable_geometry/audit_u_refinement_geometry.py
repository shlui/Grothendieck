#!/usr/bin/env python3
"""Producer-free exact u-subdivision of an accepted frozen price mosaic."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import audit_strong_shape_geometry as H


def main(path,out,parent_report):
    d=H.read(path);assert d['status']=='COMPLETE_PARENT_LINKED_U_REFINEMENT_PENDING_INDEPENDENT_AUDIT'
    pp=Path(d['parent']['path']);assert H.sha(pp)==d['parent']['sha256'];parent=H.read(pp)
    pr=H.read(parent_report);assert pr['status']=='PASS_PRODUCER_FREE_COMPLETE_STRONG_SHAPE_NATIVE_CURVED_RANK_GEOMETRY'
    assert pr['source']['sha256']==H.sha(pp)and pr['entrypoint_sha256']==H.sha(H.__file__)
    sp=Path(parent['parent']['path']);assert H.sha(sp)==parent['parent']['sha256'];source=H.read(sp)
    scanp=Path(parent['protected_transport_scan']['path']);assert H.sha(scanp)==parent['protected_transport_scan']['sha256'];scan=H.read(scanp)
    root=Path(str(scanp).split('/tonight/recovery/')[0]);b,b0=F(d['operator_b']),F(source['operator_b']);assert b==F(parent['operator_b'])==F(711,1000)<b0
    coeff={n:H.parse_profile(p,b)for n,p in parent['profiles'].items()}
    constants={(r['profile'],tuple(map(F,r['kappa_interval']))):F(r['upper'])for r in parent['cache']}
    planes={tuple(map(F,r['kappa_interval'])):[(F(p['a']),F(p['cap']))for p in r['planes']]for r in parent['source_plane_cache']}
    oldplanes=[(tuple(map(F,r['kappa_interval'])),[(F(p['a']),F(p['cap']))for p in r['planes']])for r in source['source_polygon_planes']]
    cp=root/source['contact_gate_binding']['path'];assert H.sha(cp)==source['contact_gate_binding']['sha256'];CB=F(H.read(cp)['parameters']['B'])
    betas={}
    for n,r in scan['profile_beta_bindings'].items():
        p=root/r['profile_path'];assert H.sha(p)==r['profile_sha256'];raw=F(H.read(p)['profile']['beta']);sq=F(r['beta_norm_squared'])
        assert sq==raw*raw*(6 if n.startswith('R3:')else 1);betas[n]=(sq,raw)
    upgrades={r['parent_index']:F(r['new_offer']['cap'])for r in parent['changes']}
    np=Path(parent['native_curved_binding']['path']);assert H.sha(np)==parent['native_curved_binding']['sha256'];native=H.read(np)
    groups={}
    for j,r in enumerate(native['tables']):
        ka,kb=map(F,r['kappa_interval']);grid=tuple(map(F,r['r_grid']));signature=(grid,F(r['ell']),F(r['r_quadratic']),r['source_package'],r['source_sha256'],r['native_mask_sha256'])
        groups.setdefault(r['frame'],[]).append((j,(ka,kb),signature,list(map(F,r['endpoint_uniform_caps']))))
    refinements={r['parent_index']:r for r in d['refinements']};assert len(refinements)==len(d['refinements'])==d['selected_parents']
    threshold=F(40,71);step=F(d['u_step']);assert step>0
    ranks=rankregions=rankvertices=natives=pieces=children=empties=failed=0;maximum=F(0);worst=None
    def offerprice(offer,poly,kband):
        nonlocal ranks,rankregions,rankvertices,natives,pieces
        n=offer['profile']
        if n=='POINTWISE_RANK_ENVELOPE':
            proof=offer['dominance_proof'];candidates=[];seen=set()
            for p in proof['candidates']:
                name=p['name'];ell,g,l=coeff[name];C=constants[(name,kband)]
                assert g==b and F(p['u_linear'])==-2*ell and F(p['r_linear'])==l and F(p['constant'])==C
                sig=(-2*ell,l,C);assert sig not in seen;seen.add(sig);candidates.append((name,-2*ell,l,C))
            assert candidates;regions=[];overall=None;att=set()
            for name,au,ar,C in candidates:
                reg=list(poly)
                for other,bu,br,D in candidates:
                    if other!=name:reg=H.clip_general(reg,au-bu,ar-br,D-C)
                if not reg:continue
                vals=[u*u+b*r*r+au*u+ar*r+C for u,r in reg];v=max(vals);regions.append((name,reg,vals,v))
                if overall is None or v>overall:overall=v;att=set()
                if v==overall:att.update(p for p,vv in zip(reg,vals)if vv==v)
            assert regions and len(regions)==len(proof['regions'])
            for (name,reg,vals,v),saved in zip(regions,proof['regions']):
                assert saved['profile']==name and reg==[tuple(map(F,p))for p in saved['vertices']]
                assert vals==list(map(F,saved['vertex_values']))and v==F(saved['upper']);rankregions+=1;rankvertices+=len(reg)
            assert tuple(map(F,proof['attaining_point']))in att;ranks+=1;value=H.ceil_dec(overall)
        elif n.startswith('CURVED_NATIVE:'):
            name=n[len('CURVED_NATIVE:'):];proof=offer['native_proof'];ka,kb=kband
            tables=sorted([r for r in groups[name]if r[1][0]<kb and ka<r[1][1]],key=lambda x:x[1])
            assert tables and [r[0]for r in tables]==proof['table_indices'];cursor=ka
            for j,(lo,hi),signature,caps in tables:assert lo<=cursor;cursor=max(cursor,hi)
            assert cursor>=kb;signature=tables[0][2];assert all(r[2]==signature for r in tables)
            grid,ell,g,_,_,_=signature;caps=[max(r[3][j]for r in tables)for j in range(len(grid))]
            assert F(proof['u_linear'])==-2*ell and F(proof['r_quadratic'])==g
            rlo=min(r for u,r in poly);rhi=max(r for u,r in poly);assert 0<=rlo<=rhi<=grid[-1]
            own=[];overall=None;att=set()
            for j,(ra,rb)in enumerate(zip(grid,grid[1:])):
                if rb<rlo or rhi<ra:continue
                part=H.clip_general(poly,0,-1,-ra);part=H.clip_general(part,0,1,rb)
                if not part:continue
                slope=(caps[j+1]-caps[j])/(rb-ra);C=caps[j]-slope*ra;v,points=H.quadratic_max(part,ell,g,slope,C)
                own.append((j,ra,rb,part,slope,C,v,points))
                if overall is None or v>overall:overall=v;att=set()
                if v==overall:att.update(points)
            assert own and len(own)==len(proof['pieces'])
            for (j,ra,rb,part,slope,C,v,points),saved in zip(own,proof['pieces']):
                assert j==saved['source_piece']and list(map(F,saved['source_interval']))==[ra,rb]
                assert part==[tuple(map(F,p))for p in saved['polygon']]
                assert F(saved['r_linear'])==slope and F(saved['constant'])==C and F(saved['upper'])==v
                assert tuple(map(F,saved['attaining_point']))in points;pieces+=1
            assert tuple(map(F,offer['attaining_point']))in att;natives+=1;value=H.ceil_dec(overall)
        else:
            ell,g,l=coeff[n];v,att=H.quadratic_max(poly,ell,g,l,constants[(n,kband)]);value=H.ceil_dec(v)
            assert tuple(map(F,offer['attaining_point']))in att
        assert value==F(offer['cap']);return value
    for i,row in enumerate(source['cells']):
        name=H.canonical(row['winner']);sq,raw=betas[name];rlo=F(row['r_interval'][0]);shift=sq*(1/b-1/b0)+(b-b0)*rlo*rlo
        if 'contact:'in row['winner']:
            assert name.startswith('R3:');lam=F(row['contact_credit']['lambda_value']);assert lam>=0
            shift+=(1/b-1/b0)*(lam*6*raw*CB+lam*lam*6*CB*CB/4)
        old=min(H.ceil_dec(F(row['cap'])+shift),upgrades.get(i,F(100)));cap=old
        assert(i in refinements)==(old>threshold)
        if i in refinements:
            rr=refinements[i];assert F(rr['old_offer'])==old
            for k in ['u_interval','r_interval','kappa_interval']:assert rr[k]==row[k]
            ka,kb=map(F,row['kappa_interval']);kband=(ka,kb)
            ulo,uhi=map(F,row['u_interval']);uhi=min(uhi,F(row.get('physical_u_upper',uhi)))
            rlo,rhi=map(F,row['r_interval']);rhi=min(rhi,F(row.get('physical_r_upper',rhi)))
            if'physical_polygon'in row:poly=[tuple(map(F,p))for p in row['physical_polygon']]
            else:
                poly=[(ulo,rlo),(uhi,rlo),(uhi,rhi),(ulo,rhi)];matches=[p for (a,z),p in oldplanes if a<=ka<=kb<=z];assert matches
                for a,c in matches[0]:poly=H.clip(poly,a,c)
            for a,c in planes[kband]:poly=H.clip(poly,a,c)
            assert poly;knots=[ulo]
            while knots[-1]<uhi:knots.append(min(uhi,knots[-1]+step))
            assert len(knots)>=2 and len(rr['children'])==len(knots)-1;cap=F(0)
            for child,(ua,ub)in zip(rr['children'],zip(knots,knots[1:])):
                children+=1;assert list(map(F,child['u_interval']))==[ua,ub]
                part=H.clip_general(poly,-1,0,-ua);part=H.clip_general(part,1,0,ub)
                assert part==[tuple(map(F,p))for p in child['physical_polygon']]
                if not part:assert child['status']=='EXACT_EMPTY_SOURCE_CHILD';empties+=1;continue
                assert child['status']=='EXACT_CHILD_PRICE_PENDING_AUDIT'and F(child['old_offer'])==old
                value=old
                if child['new_offer']:value=offerprice(child['new_offer'],part,kband);assert value<old
                assert value==F(child['upper']);cap=max(cap,value)
                if value>threshold:failed+=1
            assert cap==F(rr['new_parent_upper'])and cap<=old
        if cap>maximum:maximum=cap;worst=i
    assert len(source['cells'])==d['cells']==272579 and len(d['failed_children'])==failed
    assert maximum==F(d['maximum'])and maximum<F(d['strict_upper_cap'])and worst==d['worst_cell']['parent_index']
    result=dict(status='PASS_PRODUCER_FREE_COMPLETE_STRONG_SHAPE_U_REFINEMENT_GEOMETRY',source=dict(path=str(Path(path).resolve()),sha256=H.sha(path)),parent=dict(path=str(pp),sha256=H.sha(pp),accepted_geometry_report=str(Path(parent_report).resolve()),accepted_geometry_report_sha256=H.sha(parent_report)),operator_b=str(b),original_parent_cells=len(source['cells']),selected_parents=len(refinements),closed_u_children=children,empty_children=empties,rank_envelopes=ranks,rank_dominance_regions=rankregions,rank_dominance_vertices=rankvertices,native_offers=natives,native_source_pieces=pieces,maximum=str(maximum),strict_upper_cap=d['strict_upper_cap'],producer_imports=False,arithmetic='fractions.Fraction',entrypoint_sha256=H.sha(__file__),scope='Complete exact source geometry and reused accepted-parent price assembly; parent numeric and whole-q binding acceptance remain separate prerequisites.')
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--parent-report',required=True);a=p.parse_args();main(a.input,a.output,a.parent_report)
