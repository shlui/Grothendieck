#!/usr/bin/env python3
"""Producer-free exact geometry/min audit of a fresh frozen price envelope.

This imports neither the producer nor its polygon/price helpers. Directed
moment cache values and underlying whole-q gates need separate acceptance.
"""
import argparse, gzip, hashlib, json
from fractions import Fraction as F
from pathlib import Path


def read(p):
    with gzip.open(p,'rt') if str(p).endswith('.gz') else open(p) as f:
        return json.load(f)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ceil_dec(x, n=15):
    scale=10**n
    return F(-((-x.numerator*scale)//x.denominator),scale)


def canonical(s):
    while s.startswith('inherited:'):
        s=s[len('inherited:'):]
    return s.replace('native-contact:','native:',1).replace('R3-contact:','R3:',1)


def clip(polygon,a,cap):
    return clip_general(polygon,a,F(1),cap)


def clip_general(polygon,au,ar,cap):
    if not polygon:
        return []
    if au==0 and ar==0:
        return list(polygon) if cap>=0 else []
    output=[]
    prev=polygon[-1]
    fp=au*prev[0]+ar*prev[1]-cap
    for point in polygon:
        fv=au*point[0]+ar*point[1]-cap
        if (fp<=0)!=(fv<=0):
            t=fp/(fp-fv)
            output.append((prev[0]+t*(point[0]-prev[0]),
                           prev[1]+t*(point[1]-prev[1])))
        if fv<=0:
            output.append(point)
        prev,fp=point,fv
    dedup=[]
    for point in output:
        if not dedup or point!=dedup[-1]:
            dedup.append(point)
    if len(dedup)>1 and dedup[0]==dedup[-1]:
        dedup.pop()
    return dedup


def physical_moment_points(ka,kb,za,zb):
    """Independent lower-edge/horizontal-edge/curve-edge endpoint rule."""
    assert 0<=ka<=kb<=90 and F(2,5)<=za<=zb<=1
    c,d=F(1,4),F(27,2)
    km=min(kb,d/(za-c))
    if ka>km:return []
    candidates={(ka,za),(km,za)}
    kt=d/(zb-c)
    if ka<=min(km,kt):
        candidates.update([(ka,zb),(min(km,kt),zb)])
    kc=max(ka,kt)
    if kc<=km:
        candidates.update([(kc,c+d/kc),(km,c+d/km)])
    return sorted(candidates)


def quadratic_max(poly,ell,g,l,C):
    """Positive u-square rules out interior maxima; retain all edge vertices."""
    assert poly
    def value(u,r):
        return u*u-2*ell*u+g*r*r+l*r+C
    offers=[(value(*point),point) for point in poly]
    for a,b in zip(poly,poly[1:]+poly[:1]):
        du,dr=b[0]-a[0],b[1]-a[1]
        qa=du*du+g*dr*dr
        qb=(2*a[0]-2*ell)*du+(2*g*a[1]+l)*dr
        if qa<0:
            t=-qb/(2*qa)
            if 0<t<1:
                point=(a[0]+t*du,a[1]+t*dr)
                offers.append((value(*point),point))
    return max(v for v,p in offers), {p for v,p in offers if v==max(x[0] for x in offers)}


def psd(a,c,e):
    assert a>=0 and e>=0 and a*e-c*c>=0


def parse_profile(rec,b):
    path=Path(rec['profile_path'])
    assert sha(path)==rec['profile_sha256']
    p=read(path)['profile']
    p={key:[F(x) for x in v] if isinstance(v,list) else F(v)
       for key,v in p.items()}
    E=[p[x] for x in ['E00','E02','E22']] if 'E00' in p else [p['p0']**2,p['p0']*p['p2'],p['p2']**2]
    a,c,e=E
    assert p['eta']>0 and p['G']>0 and 0<b<=p['eps']
    # Extra accepted q0/q6 allowances are nonnegative; these lower guards
    # already prove the curved-price coefficients A,D nonnegative.
    assert p['major'][3]+p['zeta']**2>=0
    assert p['G']-p['d0']-p['zeta']**2>=0
    for triple in [(a,c,e),(a+1,c,e-b/2),
                   (a-p['p0']**2,c-p['p0']*p['p2'],e-p['p2']**2),
                   (p['R11'],p['R13'],p['R33']),
                   (p['R11']+1,p['R13'],p['R33']-b/6),
                   (p['R11']-b,p['R13'],p['R33'])]:
        psd(*triple)
    return p['ell'],b-a+c*c/e,-2*p['gamma']+2*(p['y']-c*p['r'])*c/e


def audit(path,output):
    data=read(path)
    assert data['status']=='COMPLETE_NEW_PRICING_PENDING_INDEPENDENT_AUDIT'
    basepath=Path(data['parent']['path']);assert sha(basepath)==data['parent']['sha256']
    source=read(basepath)
    scanpath=Path(data['protected_transport_scan']['path']);assert sha(scanpath)==data['protected_transport_scan']['sha256']
    scan=read(scanpath);assert scan['parent_sha256']==sha(basepath)
    b,b0=F(data['operator_b']),F(source['operator_b']);assert b==F(711,1000) and b<b0
    # The accepted portable base root is identified by the recorded scan path.
    root=Path(str(scanpath).split('/tonight/recovery/')[0])
    contactpath=root/source['contact_gate_binding']['path']
    assert sha(contactpath)==source['contact_gate_binding']['sha256']
    CB=F(read(contactpath)['parameters']['B'])
    coeff={name:parse_profile(rec,b) for name,rec in data['profiles'].items()}
    assert data['physical_kappa_shape_constraint'] and data['pointwise_rank_envelope']
    assert data['physical_constraint']=='4*kappa*zeta-kappa<=54'
    constants={}
    for row in data['cache']:
        assert list(map(F,row['shape_interval']))==[F(2,5),F(1)]
        key=(row['profile'],tuple(row['kappa_interval']))
        assert key not in constants
        assert row['physical_kappa_shape_constraint']
        ka,kb=map(F,row['kappa_interval']);za,zb=map(F,row['shape_interval'])
        points=physical_moment_points(ka,kb,za,zb)
        assert points==[tuple(map(F,p))for p in row['evaluation_points']]
        corners=list(map(F,row['corner_uppers']));assert len(corners)==len(points)
        C=F(row['upper']);assert C==max(corners)
        constants[key]=C
    changes={x['parent_index']:x for x in data['changes']}
    assert len(changes)==data['changed_cells']==len(data['changes'])
    plane_rows=[(tuple(map(F,r['kappa_interval'])),[(F(p['a']),F(p['cap']))for p in r['planes']])
                for r in source['source_polygon_planes']]
    bankpath=root/'tonight/moments/linear_contact_bank/source_bank.json'
    assert str(bankpath)==str(root/source['linear_source_bank_binding']['path'])
    assert sha(bankpath)==source['linear_source_bank_binding']['sha256']
    bank=read(bankpath)
    fresh_planes={}
    for row in data['source_plane_cache']:
        assert row['physical_constraint']=='4*kappa*zeta-kappa<=54'
        assert list(map(F,row['shape_interval']))==[F(2,5),F(1)]
        ka,kb=map(F,row['kappa_interval'])
        points=physical_moment_points(ka,kb,F(2,5),F(1))
        assert points==[tuple(map(F,p))for p in row['evaluation_points']]
        key=tuple(row['kappa_interval']);assert key not in fresh_planes
        used=[]
        for plane in row['planes']:
            qrow=bank['rows'][plane['source_row']]
            assert F(qrow['a'])==F(plane['a'])>0
            assert F(qrow['polynomial_coefficients_q0_q2_q4_q6'][3])>=0
            prices=list(map(F,plane['corner_uppers']))
            assert len(prices)==len(points) and F(plane['cap'])==max(prices)
            used.append((F(plane['a']),F(plane['cap'])))
        assert used
        fresh_planes[key]=used
    def polygon(row):
        if 'physical_polygon' in row:
            return [tuple(map(F,point))for point in row['physical_polygon']]
        ulo,uhi=map(F,row['u_interval']);uhi=min(uhi,F(row.get('physical_u_upper',uhi)))
        rlo,rhi=map(F,row['r_interval']);rhi=min(rhi,F(row.get('physical_r_upper',rhi)))
        poly=[(ulo,rlo),(uhi,rlo),(uhi,rhi),(ulo,rhi)]
        ka,kb=map(F,row['kappa_interval'])
        matches=[p for (a,z),p in plane_rows if a<=ka<=kb<=z]
        assert matches
        for a,c in matches[0]:
            poly=clip(poly,a,c)
        return poly
    betas={}
    for name,rec in scan['profile_beta_bindings'].items():
        profilepath=root/rec['profile_path'];assert sha(profilepath)==rec['profile_sha256']
        rawbeta=F(read(profilepath)['profile']['beta'])
        beta2=F(rec['beta_norm_squared'])
        assert beta2==rawbeta*rawbeta*(6 if name.startswith('R3:') else 1)
        betas[name]=(beta2,rawbeta)
    maximum=F(0);worst=None;contacts=0;rank_offers=0;rank_regions=0;rank_vertices=0;fresh_empty=0
    native_binding=data['native_curved_binding'];nativepath=Path(native_binding['path'])
    assert sha(nativepath)==native_binding['sha256']
    native_data=read(nativepath);assert F(native_data['operator_b'])==b
    assert native_data['physical_constraint']=='kappa*zeta<=36'
    native_tables=native_data['tables'];native_groups={}
    for idx,row in enumerate(native_tables):
        ka,kb=map(F,row['kappa_interval']);assert 0<=ka<kb<=90
        grid=list(map(F,row['r_grid']));assert grid[0]==0 and all(x<y for x,y in zip(grid,grid[1:]))
        assert len(grid)==len(row['endpoint_uniform_caps'])
        signature=(tuple(grid),F(row['ell']),F(row['r_quadratic']),
                   row['source_package'],row['source_sha256'],row['native_mask_sha256'])
        native_groups.setdefault(row['frame'],[]).append((idx,(ka,kb),signature,list(map(F,row['endpoint_uniform_caps']))))
    native_offers=0;native_pieces=0;native_vertices=0
    threshold=F(data['threshold']);oldfail=newfail=0
    for i,row in enumerate(source['cells']):
        name=canonical(row['winner']);beta2,beta=betas[name]
        rlo=F(row['r_interval'][0]);shift=beta2*(1/b-1/b0)+(b-b0)*rlo*rlo
        if 'contact:' in row['winner']:
            # Protected final29 has no selected native-contact square-root term.
            assert name.startswith('R3:')
            lam=F(row['contact_credit']['lambda_value']);assert lam>=0
            shift+=(1/b-1/b0)*(lam*6*beta*CB+lam*lam*6*CB*CB/4)
            contacts+=1
        old=ceil_dec(F(row['cap'])+shift)
        cap=old
        if old>threshold:oldfail+=1
        if i in changes:
            change=changes[i]
            for k in ['u_interval','r_interval','kappa_interval']:
                assert change[k]==row[k]
            assert change['old_winner']==row['winner']
            assert F(change['old_b711_offer'])==old and old>threshold
            poly=polygon(row);assert poly
            for a,c in fresh_planes[tuple(row['kappa_interval'])]:
                poly=clip(poly,a,c)
            assert poly==[tuple(map(F,point))for point in change['physical_polygon']]
            offer=change['new_offer'];n=offer['profile']
            if n=='EXACT_EMPTY_PHYSICAL_SOURCE_CELL':
                assert not poly and F(offer['cap'])==0 and offer['attaining_point']==[]
                fresh_empty+=1
            elif n.startswith('CURVED_NATIVE:'):
                assert poly;name=n[len('CURVED_NATIVE:'):];proof=offer['native_proof']
                ka,kb=map(F,row['kappa_interval']);assert ka<kb
                selected=sorted([t for t in native_groups[name] if t[1][0]<kb and ka<t[1][1]],key=lambda t:t[1])
                assert selected and [t[0]for t in selected]==proof['table_indices']
                cursor=ka
                for idx,band,signature,caps in selected:
                    assert band[0]<=cursor
                    cursor=max(cursor,band[1])
                assert cursor>=kb
                signature=selected[0][2]
                assert all(t[2]==signature for t in selected)
                grid,ell,g,_,_,_=signature
                caps=[max(t[3][j] for t in selected)for j in range(len(grid))]
                assert -2*ell==F(proof['u_linear']) and g==F(proof['r_quadratic'])
                rlo=min(r for u,r in poly);rhi=max(r for u,r in poly)
                assert 0<=rlo<=rhi<=grid[-1]
                pieces=[];overall=None;points=set()
                for j,(ra,rb) in enumerate(zip(grid,grid[1:])):
                    if rb<rlo or rhi<ra:continue
                    piece=clip_general(poly,0,-1,-ra)
                    piece=clip_general(piece,0,1,rb)
                    if not piece:continue
                    slope=(caps[j+1]-caps[j])/(rb-ra)
                    C=caps[j]-slope*ra
                    exact,attaining=quadratic_max(piece,ell,g,slope,C)
                    pieces.append((j,ra,rb,piece,slope,C,exact,attaining))
                    if overall is None or exact>overall:overall=exact;points=set()
                    if exact==overall:points.update(attaining)
                assert pieces and len(pieces)==len(proof['pieces'])
                for own,saved in zip(pieces,proof['pieces']):
                    j,ra,rb,piece,slope,C,exact,attaining=own
                    assert j==saved['source_piece']
                    assert [ra,rb]==list(map(F,saved['source_interval']))
                    assert piece==[tuple(map(F,p))for p in saved['polygon']]
                    assert slope==F(saved['r_linear']) and C==F(saved['constant'])
                    assert exact==F(saved['upper'])
                    assert tuple(map(F,saved['attaining_point'])) in attaining
                    native_pieces+=1;native_vertices+=len(piece)
                assert ceil_dec(overall)==F(offer['cap'])<old
                assert tuple(map(F,offer['attaining_point'])) in points
                native_offers+=1
            elif n=='POINTWISE_RANK_ENVELOPE':
                assert poly;proof=offer['dominance_proof'];ranks=proof['candidates']
                candidates=[];seen=set()
                for c in ranks:
                    name=c['name'];ell,g,l=coeff[name]
                    C=constants[(name,tuple(row['kappa_interval']))]
                    assert g==b and F(c['u_linear'])==-2*ell
                    assert F(c['r_linear'])==l and F(c['constant'])==C
                    signature=(-2*ell,l,C);assert signature not in seen
                    seen.add(signature);candidates.append((name,-2*ell,l,C))
                assert candidates
                regions=[];overall=None;points=set()
                for name,au,ar,C in candidates:
                    region=list(poly)
                    for other,bu,br,D in candidates:
                        if other==name:continue
                        region=clip_general(region,au-bu,ar-br,D-C)
                    if not region:continue
                    values=[u*u+b*r*r+au*u+ar*r+C for u,r in region]
                    exact=max(values)
                    regions.append((name,region,values,exact))
                    if overall is None or exact>overall:overall=exact;points=set()
                    if exact==overall:points.update(p for p,v in zip(region,values)if v==exact)
                assert regions and len(regions)==len(proof['regions'])
                for own,saved in zip(regions,proof['regions']):
                    name,region,values,exact=own
                    assert saved['profile']==name
                    assert region==[tuple(map(F,p))for p in saved['vertices']]
                    assert values==list(map(F,saved['vertex_values']))
                    assert exact==F(saved['upper'])
                    rank_regions+=1;rank_vertices+=len(region)
                assert ceil_dec(overall)==F(offer['cap'])<old
                assert tuple(map(F,offer['attaining_point'])) in points
                assert offer['attaining_point']==proof['attaining_point']
                rank_offers+=1
            else:
                assert poly
                ell,g,l=coeff[n];C=constants[(n,tuple(row['kappa_interval']))]
                exact,points=quadratic_max(poly,ell,g,l,C)
                assert ceil_dec(exact)==F(offer['cap'])<old
                assert tuple(map(F,offer['attaining_point'])) in points
            cap=F(offer['cap'])
        if cap>threshold:newfail+=1
        if cap>maximum:maximum=cap;worst=i
    assert len(source['cells'])==data['cells']==272579
    assert oldfail==data['old_failed_cells'] and newfail==data['new_failed_cells']
    assert maximum==F(data['maximum'])
    assert F(data['strict_upper_cap'])>maximum
    assert data['worst_cell']['parent_index']==worst
    assert F(data['worst_cell']['cap'])==maximum
    result={'status':'PASS_PRODUCER_FREE_COMPLETE_STRONG_SHAPE_NATIVE_CURVED_RANK_GEOMETRY',
       'scope':'Exact price/geometry assembly; fresh directed cache and whole-q input acceptance remain separate prerequisites.',
       'source':{'path':str(Path(path).resolve()),'sha256':sha(path)},
       'protected_parent':{'path':str(basepath),'sha256':sha(basepath)},
       'operator_b':str(b),'original_physical_parent_cells':len(source['cells']),
       'physical_moment_constraint':'4*kappa*zeta-kappa<=54',
       'new_selected_polygon_offers':len(changes),'selected_contact_transports':contacts,
       'six_PSD_profile_checks':len(coeff),'directed_constant_records':len(constants),
       'uniform_source_plane_records':len(fresh_planes),
       'exact_rank_envelope_offers':rank_offers,'exact_rank_dominance_regions':rank_regions,
       'exact_rank_dominance_vertices':rank_vertices,'fresh_exact_empty_cells':fresh_empty,
       'native_curved_input':{'path':str(nativepath),'sha256':sha(nativepath)},
       'curved_native_selected_offers':native_offers,'curved_native_source_pieces':native_pieces,
       'curved_native_piece_vertices':native_vertices,
       'source_plane_bank':{'path':str(bankpath),'sha256':sha(bankpath)},
       'maximum':str(maximum),'strict_upper_cap':data['strict_upper_cap'],
       'worst_parent_index':worst,'old_failed_cells':oldfail,'new_failed_cells':newfail,
       'geometry_and_native_knot_basis':'All protected final29 parent cells remain. Physical-shape plane clips and closed rank dominance subpolygons retain every physical point; each new native offer splits every inherited radial knot and first maximizes endpoint caps over a complete closed cumulant-band cover.',
       'producer_imports':False,'arithmetic':'Python fractions.Fraction',
       'entrypoint_sha256':sha(__file__)}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output')
    a=p.parse_args();audit(a.input,a.output)
