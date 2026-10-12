"""Producer-free exact cumulant-child partition, plane and native price audit."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
from audit_fullupper_polygon_envelope import ROOT,read,sha,key,iv,verify_polygon,native_polygon_price,moment_price,frac_upper

def main(parent_path,path,parent_audit,native_audit,output):
    old=read(parent_path);d=read(path);a=read(parent_audit);na=read(native_audit)
    assert a['status'].startswith('PASS_PRODUCER_FREE_')and a['source_sha256']==sha(parent_path)
    assert na['status'].startswith('PASS_INDEPENDENT_ALL_')and na['precision_bits']==768
    accepted={t['mask']:t for t in na['frames']}
    chain=d['native_kappa_refinement'];assert Path(chain['parent_path']).resolve()==parent_path.resolve()and chain['parent_sha256']==sha(parent_path)
    assert F(chain['prior_strict_cap'])==F(old['strict_rational_cap'])
    jobs=ROOT/chain['native_jobs'];assert sha(jobs)==chain['native_jobs_sha256'];j=read(jobs);points=list(map(F,j['kappa_points']))
    assert points[0]==0 and points[-1]==90 and all(v<w for v,w in zip(points,points[1:]))
    assert len(points)==68
    for k in ('operator_b','amplitude_interval','native_inputs','native_endpoint_inputs','physical_u_bounds','source_polygon_planes','linear_source_bank_binding','linear_exclusion_binding','contact_gate_binding','accepted_R3_binding','cached_corner_constants','profile_coefficients'):
        assert d[k]==old[k],k
    assert d['physically_empty_cells'][:len(old['physically_empty_cells'])]==old['physically_empty_cells']
    b=F(d['operator_b']);assert b<=F(na['operator_b'])==F(7171,10000)
    frames={};oldpackages={t['package']:t for t in old['native_endpoint_inputs']}
    for t in chain['new_native_bindings']:
        mask=ROOT/t['mask'];pkg=ROOT/t['package'];assert sha(mask)==t['mask_sha256']and t['package']in oldpackages
        v=accepted[t['mask']];assert v['mask_sha256']==t['mask_sha256']and v['package']==t['package']
        sourcehash=sha(pkg/'data/source_profile.json');assert v['source_sha256']==sourcehash==oldpackages[t['package']]['source_sha256']
        assert sha(pkg/'replay.json')==oldpackages[t['package']]['replay_sha256']
        m=read(mask);p=read(pkg/'data/local_summary.json')['profile'];src=read(pkg/'data/source_profile.json');baseline=F(m['operator_b']);assert baseline==F(na['operator_b'])and b<=F(p['epsilon'])
        assert list(map(F,m['kappa_grid']))==points
        shift=F(p['beta'])**2*(1/b-1/baseline)
        frames[t['mask']]={'rgrid':[F(str(v))for v in src['r_grid']],'endpoints':{iv(z['kappa_interval']):[[F(v)+shift for v in row]for row in z['prices']]for z in m['endpoint_prices']},'g':F(m['exact_gram_cost'])+b-baseline,'ell':F(p['ell'])}
    assert len(frames)==7 and set(frames)=={x['mask']for x in j['jobs']}
    bankpath=ROOT/d['linear_source_bank_binding']['path'];assert sha(bankpath)==d['linear_source_bank_binding']['sha256'];bank=read(bankpath)
    planes={};planechecks=0
    for band in chain['new_plane_cache']:
        k=iv(band['kappa_interval']);assert k not in planes;vals=[]
        for x in band['planes']:
            row=bank['rows'][x['source_row']];a,cap=F(x['a']),F(x['cap']);assert a==F(row['a'])>0
            q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']));assert q[3]>=0
            for kap in k:
                for s in (F(0),F(1)):
                    assert frac_upper(moment_price(q,kap,s))<=cap;planechecks+=1
            vals.append((a,cap))
        assert len(vals)==10 and len({v[0]for v in vals})==10;planes[k]=vals
    refs={};expected=set()
    for x in chain['refinements']:
        n=x['parent_index'];assert isinstance(n,int)and 0<=n<len(old['cells'])and n not in refs;p=old['cells'][n];k=iv(p['kappa_interval']);knots=list(map(F,x['kappa_knots']))
        assert knots==[k[0]]+[v for v in points if k[0]<v<k[1]]+[k[1]]and len(knots)>2
        assert F(p['cap'])>F(chain['threshold']);refs[n]=p
        for lo,hi in zip(knots,knots[1:]):expected.add((n,lo,hi));assert (lo,hi)in planes
    untouched={key(p):p for n,p in enumerate(old['cells'])if n not in refs};assert len(untouched)==len(old['cells'])-len(refs)
    seen=set();seenuntouched=set();maximum=F(0);pricedchildren=0;emptychildren=0;knotsplit=0;stationary=0;newnative=0;inherited=0
    def child(row,empty=False):
        n=row['kappa_refined_parent_index'];assert n in refs;p=refs[n];k=iv(row['kappa_interval']);tag=(n,*k);assert tag in expected and tag not in seen;seen.add(tag)
        assert row['u_interval']==p['u_interval']and row['r_interval']==p['r_interval']and list(map(F,row['new_plane_band']))==list(k)
        ul,uh=iv(p['u_interval']);rl,rh=iv(p['r_interval']);uh=min(uh,F(p['physical_u_upper']));rh=min(rh,F(p['physical_r_upper']))
        if empty:assert row['polygon_empty'];polygon=verify_polygon([],(ul,uh),(rl,rh),planes[k])
        else:
            assert row['physical_u_upper']==p['physical_u_upper']and row['physical_r_upper']==p['physical_r_upper']and row['polygon_pricing']is True
            polygon=verify_polygon([tuple(map(F,z))for z in row['physical_polygon']],(ul,uh),(rl,rh),planes[k]);assert polygon
        return p,(ul,uh),(rl,rh),k,polygon
    for row in d['physically_empty_cells'][len(old['physically_empty_cells']):]:child(row,True);emptychildren+=1
    for row in d['cells']:
        cap=F(row['cap'])
        if 'kappa_refined_parent_index'not in row:
            k=key(row);assert k in untouched and k not in seenuntouched and row==untouched[k];seenuntouched.add(k)
        else:
            p,ui,ri,k,poly=child(row);pricedchildren+=1;assert row['inherited_kappa_cap']==p['cap']and cap<=F(p['cap'])
            if row['winner'].startswith('inherited:'):assert row['winner']=='inherited:'+p['winner']and cap==F(p['cap']);inherited+=1
            else:
                kind,name=row['winner'].split(':',1);assert kind=='native'and name in frames and 'contact_credit'not in row
                val,n,s=native_polygon_price(frames[name],ui,ri,k,poly,planes[k]);assert val==cap,(key(row),val,cap)
                knotsplit+=n;stationary+=s;newnative+=1
        assert row['target_pass']==(cap<F(40,71))and row['fallback_pass']==(cap<F(451,800));maximum=max(maximum,cap)
    assert seen==expected and seenuntouched==set(untouched)
    assert maximum==F(d['maximum_cell_cap'])and maximum<F(d['strict_rational_cap'])<=F(old['strict_rational_cap'])
    assert d['target_cells']==sum(x['target_pass']for x in d['cells'])and d['fallback_cells']==sum(x['fallback_pass']for x in d['cells'])
    out={'status':'PASS_PRODUCER_FREE_NATIVE_CUMULANT_REFINEMENT','source':{'path':str(path.relative_to(ROOT)),'sha256':sha(path)},'source_sha256':sha(path),'parent':{'path':str(parent_path.relative_to(ROOT)),'sha256':sha(parent_path),'accepted_audit':str(parent_audit.relative_to(ROOT)),'audit_sha256':sha(parent_audit)},'operator_b':str(b),'amplitude_interval':d['amplitude_interval'],'new_native_numeric_audit':{'path':str(native_audit.relative_to(ROOT)),'sha256':sha(native_audit)},'cumulant_grid_points':len(points),'refined_parents':len(refs),'exact_child_intervals':len(expected),'unchanged_price_cells':len(untouched),'priced_children':pricedchildren,'exact_empty_children':emptychildren,'total_priced_cells':len(d['cells']),'total_empty_cells':len(d['physically_empty_cells']),'new_plane_bands':len(planes),'new_uniform_planes':10*len(planes),'new_plane_directed_corner_checks':planechecks,'new_native_selected_children':newnative,'inherited_selected_children':inherited,'native_knot_pieces':knotsplit,'edge_stationary_points':stationary,'maximum_cell_cap':str(maximum),'strict_rational_cap':d['strict_rational_cap'],'proof':'The complete accepted parent cover is retained or exactly partitioned in cumulant. All new planes are independently replayed valid fixed analytic rows. Every selected old parent cap remains physically valid; every new native cap is independently rederived on its complete new source polygon, splitting every radial knot and including edge stationary points.'}
    output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent');p.add_argument('path');p.add_argument('--parent-audit',required=True);p.add_argument('--native-audit',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(ROOT/a.parent,ROOT/a.path,ROOT/a.parent_audit,ROOT/a.native_audit,ROOT/a.output)
