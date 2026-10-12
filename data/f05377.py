"""Producer-free parent-linked refinement audit of a complete accepted cover."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter,defaultdict
import argparse,json
from audit_fullupper_polygon_envelope import ROOT,read,sha,key,iv,ceil12,polygon_max,verify_polygon,native_polygon_price

def main(parent_path,path,parent_audit,output):
    old=read(parent_path);d=read(path);a=read(parent_audit)
    assert a['status'].startswith('PASS_PRODUCER_FREE_')
    h=a['source_sha256'] if 'source_sha256' in a else a['source']['sha256'];assert h==sha(parent_path)
    chain=d['adaptive_amplitude_refinement'];assert Path(chain['parent_path']).resolve()==parent_path.resolve()
    assert chain['parent_sha256']==sha(parent_path) and F(chain['u_step'])>0
    assert F(chain['prior_strict_cap'])==F(old['strict_rational_cap'])
    for k in ('operator_b','amplitude_interval','native_inputs','native_endpoint_inputs','physical_u_bounds','source_polygon_planes','linear_source_bank_binding','linear_exclusion_binding','contact_gate_binding','accepted_R3_binding','cached_corner_constants','profile_coefficients'):
        assert d[k]==old[k],k
    assert d['physically_empty_cells'][:len(old['physically_empty_cells'])]==old['physically_empty_cells']
    refinements={};intervals={};expected=set()
    for x in chain['refinements']:
        n=x['parent_index'];assert isinstance(n,int)and 0<=n<len(old['cells'])and n not in refinements
        p=old['cells'][n];u=iv(p['u_interval']);knots=list(map(F,x['u_knots']))
        assert knots[0]==u[0]and knots[-1]==u[1]and all(z<w for z,w in zip(knots,knots[1:]))
        assert F(p['cap'])>F(chain['threshold'])and u[1]-u[0]>F(chain['u_step'])
        for z,w in zip(knots,knots[1:]):assert w-z<=F(chain['u_step']);expected.add((n,z,w))
        refinements[n]=p;intervals[n]=knots
    untouched={key(p):p for n,p in enumerate(old['cells'])if n not in refinements};assert len(untouched)==len(old['cells'])-len(refinements)
    b=F(d['operator_b']);frames={}
    for t in d['native_endpoint_inputs']:
        mp=ROOT/t['mask'];pk=ROOT/t['package'];assert sha(mp)==t['mask_sha256']and sha(pk/'data/source_profile.json')==t['source_sha256']and sha(pk/'replay.json')==t['replay_sha256']
        m=read(mp);p=read(pk/'data/local_summary.json')['profile'];src=read(pk/'data/source_profile.json');baseline=F(m['operator_b']);assert b<=baseline
        shift=F(p['beta'])**2*(1/b-1/baseline)
        end={iv(z['kappa_interval']):[[F(v)+shift for v in row]for row in z['prices']]for z in m['endpoint_prices']}
        frames[t['mask']]={'rgrid':list(map(lambda z:F(str(z)),src['r_grid'])),'endpoints':end,'g':F(m['exact_gram_cost'])+b-baseline,'ell':F(p['ell'])}
    cache={(t['profile'],iv(t['kappa_interval'])):F(t['upper'])for t in d['cached_corner_constants']}
    planes={iv(x['kappa_interval']):[(F(z['a']),F(z['cap']))for z in x['planes']]for x in d['source_polygon_planes']}
    seen=set();seen_untouched=set();maximum=F(0);winner_counts=Counter();priced_children=0;empty_children=0;knots=0;stationary=0
    def child_geometry(row,empty=False):
        n=row['refined_parent_index'];assert n in refinements;p=refinements[n]
        lo,hi=iv(row['u_interval']);tag=(n,lo,hi);assert tag in expected and tag not in seen;seen.add(tag)
        assert row['r_interval']==p['r_interval']and row['kappa_interval']==p['kappa_interval']
        assert lo>=F(p['u_interval'][0])and hi<=F(p['u_interval'][1])
        uh=min(hi,F(p['physical_u_upper']));rl,rh=iv(p['r_interval']);rh=min(rh,F(p['physical_r_upper']));kk=iv(p['kappa_interval'])
        if uh<lo:
            assert empty and row.get('inherited_parent_exclusion')and F(row['physical_u_upper'])==uh
            return p,(lo,uh),(rl,rh),kk,None
        if empty:
            assert row.get('inherited_parent_exclusion')and row.get('polygon_empty');verify_polygon([],(lo,uh),(rl,rh),planes[kk]);polygon=[]
        else:
            assert F(row['physical_u_upper'])==uh and F(row['physical_r_upper'])==rh and row['polygon_pricing'] is True
            polygon=verify_polygon([tuple(map(F,z))for z in row['physical_polygon']],(lo,uh),(rl,rh),planes[kk]);assert polygon
        return p,(lo,uh),(rl,rh),kk,polygon
    for row in d['physically_empty_cells'][len(old['physically_empty_cells']):]:child_geometry(row,True);empty_children+=1
    for row in d['cells']:
        cap=F(row['cap'])
        if 'refined_parent_index'not in row:
            k=key(row);assert k in untouched and k not in seen_untouched and row==untouched[k];seen_untouched.add(k)
        else:
            p,ui,ri,kk,polygon=child_geometry(row);priced_children+=1
            assert row['inherited_previous_cap']==p['cap']and row['inherited_previous_winner']==p['winner']and cap<=F(p['cap'])
            w=row['winner']
            if w.startswith('inherited:'):assert w=='inherited:'+p['winner']and cap==F(p['cap'])
            else:
                kind,name=w.split(':',1);assert 'contact_credit'not in row
                if kind=='native':
                    val,n,s=native_polygon_price(frames[name],ui,ri,kk,polygon,planes[kk]);knots+=n;stationary+=s
                else:
                    assert kind=='R3';cf=d['profile_coefficients'][name];g,lin,ell=(F(cf[z])for z in('r_quadratic','r_linear','ell'))
                    v,s=polygon_max(polygon,-2*ell,g,lin,cache[name,kk]);val=ceil12(v);stationary+=s
                assert val==cap,(key(row),w,val,cap)
        assert row['target_pass']==(cap<F(40,71))and row['fallback_pass']==(cap<F(451,800))
        maximum=max(maximum,cap);winner_counts[row['winner']]+=1
    assert seen==expected and seen_untouched==set(untouched)
    assert maximum==F(d['maximum_cell_cap'])and maximum<F(d['strict_rational_cap'])<=F(old['strict_rational_cap'])
    assert d['target_cells']==sum(c['target_pass']for c in d['cells'])and d['fallback_cells']==sum(c['fallback_pass']for c in d['cells'])
    out={'status':'PASS_PRODUCER_FREE_ADAPTIVE_POLYGON_REFINEMENT','source':{'path':str(path.relative_to(ROOT)),'sha256':sha(path)},'source_sha256':sha(path),'parent':{'path':str(parent_path.relative_to(ROOT)),'sha256':sha(parent_path),'accepted_audit':str(parent_audit.relative_to(ROOT)),'audit_sha256':sha(parent_audit)},'operator_b':str(b),'amplitude_interval':d['amplitude_interval'],'refined_parents':len(refinements),'unchanged_price_cells':len(untouched),'exact_child_intervals':len(expected),'priced_children':priced_children,'exact_empty_children':empty_children,'total_priced_cells':len(d['cells']),'total_empty_cells':len(d['physically_empty_cells']),'native_knot_pieces':knots,'edge_stationary_points':stationary,'maximum_cell_cap':str(maximum),'strict_rational_cap':d['strict_rational_cap'],'proof':'Every accepted parent cell is retained unchanged or exactly partitioned into its recorded child amplitude intervals. Original physical constraints are inherited; each new clipped polygon is independently verified. Selected old parent bounds or recomputed exact cached polynomial maxima bound every nonempty child.'}
    output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent');p.add_argument('path');p.add_argument('--parent-audit',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(ROOT/a.parent,ROOT/a.path,ROOT/a.parent_audit,ROOT/a.output)
