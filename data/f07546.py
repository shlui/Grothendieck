"""Producer-free monotone delta audit; parent cover remains the geometric proof."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,sys
from audit_fullupper_polygon_envelope import ROOT,read,sha,key,iv,ceil12,polygon_max,qmax
from audit_general_r3_parameters_tail import gates
from audit_general_r3_source_bands import price
from mpfr_interval import UP

def main(parent_path,path,parent_audit,cache_audit,output):
    parent=read(parent_path);d=read(path);pa=read(parent_audit);ca=read(cache_audit) if cache_audit else None
    assert pa['status'].startswith('PASS_PRODUCER_FREE_')
    assert (pa['source_sha256'] if 'source_sha256' in pa else pa['source']['sha256'])==sha(parent_path)
    if ca:
        assert ca['status'].startswith('PASS_INDEPENDENT_')
        assert ca['source']=={'path':str(path.relative_to(ROOT)),'sha256':sha(path)}
    up=d['monotone_upgrade'];assert Path(up['parent_path']).resolve()==parent_path.resolve()
    assert up['parent_sha256']==sha(parent_path)
    b=F(d['operator_b']);assert d['operator_b']==parent['operator_b']
    new=set(up['new_profiles']);assert new==set(d['profile_coefficients'])-set(parent['profile_coefficients'])
    assert len(new)==len(up['new_profiles'])
    for k in ('native_inputs','native_endpoint_inputs','physical_u_bounds','physically_empty_cells','source_polygon_planes','linear_source_bank_binding','linear_exclusion_binding','contact_gate_binding','amplitude_interval'):
        assert d[k]==parent[k],k
    assert d['accepted_R3_binding'][:len(parent['accepted_R3_binding'])]==parent['accepted_R3_binding']
    assert d['cached_corner_constants'][:len(parent['cached_corner_constants'])]==parent['cached_corner_constants']
    for n in parent['profile_coefficients']:assert d['profile_coefficients'][n]==parent['profile_coefficients'][n]
    reports=[];profiles={}
    for rec in d['accepted_R3_binding'][len(parent['accepted_R3_binding']):]:
        p=ROOT/rec['path'];assert sha(p)==rec['sha256'];v=read(p)
        if 'profile' in v:profiles[sha(p)]=(p,v['profile'])
        else:reports.append((p,v))
    assert len(profiles)==len(reports)==len(new)
    accepted={};new_parameters={}
    for reportpath,v in reports:
        assert v['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        pth,raw=profiles[v['profile_sha256']]
        assert Path(v['profile_path']).resolve()==pth.resolve()
        p={k:list(map(F,x)) if isinstance(x,list) else F(x) for k,x in raw.items()}
        assert b<=p['eps'];p['eps']=b
        p.update(E00=p['p0']**2,E02=p['p0']*p['p2'],E22=p['p2']**2);gates(p)
        r0,r6=map(F,v['reserves']);assert r0>=0 and r6>=0
        A=p['major'][3]+r6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
        assert A>=0 and B>=0
        yc=p['y']-p['E02']*p['r'];coef=(b,-2*p['gamma']+2*yc*p['E02']/p['E22'],p['ell'])
        matching=[n for n in new if tuple(F(d['profile_coefficients'][n][k])for k in('r_quadratic','r_linear','ell'))==coef]
        assert len(matching)==1;new_parameters[matching[0]]=(p,r0,r6);accepted[matching[0]]={'path':str(pth.relative_to(ROOT)),'sha256':sha(pth),'reserves':list(map(str,(r0,r6))), 'report_path':str(reportpath.relative_to(ROOT)),'report_sha256':sha(reportpath)}
    cache={}
    for t in d['cached_corner_constants']:
        k=t['profile'],iv(t['kappa_interval']);assert k not in cache;cache[k]=F(t['upper'])
    checked_corners=0
    if ca:
        tables=ca.get('new_tables',ca.get('corners'));assert len(tables)==len(cache)
        for t in tables:assert cache[t['profile'],iv(t['kappa_interval'])]==F(t['saved_upper']) and len(t['corners'])==4
    else:
        for (name,k),bound in cache.items():
            if name not in new:continue
            p,r0,r6=new_parameters[name]
            for kap in k:
                for z in (F(2,5),F(1)):
                    upper=F(price(p,F(0),F(0),kap,z,r0,r6).hi.text(UP,240))
                    assert upper<=bound,(name,k,z,upper,bound)
                    checked_corners+=1
        assert checked_corners==4*sum(name in new for name,k in cache)
    assert len(parent['cells'])==len(d['cells'])
    mutable={'cap','winner','target_pass','fallback_pass','previous_cap','previous_winner','contact_credit'}
    improved=0;maximum=F(0);newwinner={};polygons=0
    for old,row in zip(parent['cells'],d['cells']):
        assert key(old)==key(row)
        assert {k:v for k,v in old.items()if k not in mutable}=={k:v for k,v in row.items()if k not in mutable}
        cap=F(row['cap']);oldcap=F(old['cap']);assert cap<=oldcap
        if cap==oldcap:assert row==old
        else:
            improved+=1;assert row['previous_cap']==old['cap'] and row['previous_winner']==old['winner']
            assert 'contact_credit' not in row
            kind,name=row['winner'].split(':',1);assert kind=='R3' and name in new
            ulo,uhi=iv(row['u_interval']);rlo,rhi=iv(row['r_interval']);k=iv(row['kappa_interval'])
            uhi=min(uhi,F(row['physical_u_upper']));rhi=min(rhi,F(row['physical_r_upper']))
            assert ulo<=uhi and rlo<=rhi
            a,l,e=(F(d['profile_coefficients'][name][z])for z in('r_quadratic','r_linear','ell'))
            if row['polygon_pricing']:
                points=[tuple(map(F,p))for p in row['physical_polygon']]
                value=ceil12(polygon_max(points,-2*e,a,l,cache[name,k])[0]);polygons+=1
            else:value=ceil12(cache[name,k]+max(u*u-2*e*u for u in (ulo,uhi))+qmax(a,l,F(0),rlo,rhi))
            assert value==cap,(key(row),value,cap)
            newwinner[name]=newwinner.get(name,0)+1
        assert row['target_pass']==(cap<F(40,71)) and row['fallback_pass']==(cap<F(451,800))
        maximum=max(maximum,cap)
    assert improved==up['improved_cells'] and maximum==F(d['maximum_cell_cap'])
    assert maximum<F(d['strict_rational_cap'])
    assert d['target_cells']==sum(c['target_pass']for c in d['cells'])
    assert d['fallback_cells']==sum(c['fallback_pass']for c in d['cells'])
    out={'status':'PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE','source':{'path':str(path.relative_to(ROOT)),'sha256':sha(path)},'source_sha256':sha(path),'parent':{'path':str(parent_path.relative_to(ROOT)),'sha256':sha(parent_path),'accepted_geometry_audit':str(parent_audit.relative_to(ROOT)),'audit_sha256':sha(parent_audit)},'numeric_cache_audit':({'path':str(cache_audit.relative_to(ROOT)),'sha256':sha(cache_audit)} if cache_audit else {'status':'PASS_FRESH_PRODUCER_FREE_NEW_MPFR768_CORNERS','fresh_corners':checked_corners,'unchanged_parent_cache_entries':len(parent['cached_corner_constants'])}),'operator_b':str(b),'amplitude_interval':d['amplitude_interval'],'identical_priced_cells':len(d['cells']),'identical_empty_cells':len(d['physically_empty_cells']),'improved_cells':improved,'new_winner_counts':newwinner,'improved_polygon_prices':polygons,'new_profiles':accepted,'maximum_cell_cap':str(maximum),'strict_rational_cap':d['strict_rational_cap'],'proof':'The accepted parent provides the identical physical partition. Every changed price is an independently accepted rank-profile polynomial maximum on that same polygon or containing rectangle. All unchanged old prices remain valid.'}
    output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent');p.add_argument('path');p.add_argument('--parent-audit',required=True);p.add_argument('--cache-audit');p.add_argument('--output',required=True);a=p.parse_args()
    main(ROOT/a.parent,ROOT/a.path,ROOT/a.parent_audit,ROOT/a.cache_audit if a.cache_audit else None,ROOT/a.output)
