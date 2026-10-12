#!/usr/bin/env python3
"""Independent simplex-vertex enumeration and exact polygon edge maxima."""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import json
import os
import sys

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest
import replay_native_R3_unions as N
import replay_R3_corners as R
from replay_joint_contact_full import bound,native_curve


def polygon(constraints):
    points=set()
    for (a,b,c),(d,e,f)in combinations(constraints,2):
        det=a*e-b*d
        if not det:continue
        p=((c*e-b*f)/det,(a*f-c*d)/det)
        if all(x*p[0]+y*p[1]<=z for x,y,z in constraints):points.add(p)
    # Exact monotone convex hull, including a segment or point degeneration.
    if len(points)<3:return sorted(points)
    def cross(a,b,c):return(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    lower=[];upper=[]
    for p in sorted(points):
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    for p in sorted(points,reverse=True):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return lower[:-1]+upper[:-1]


def value(q,u,r):return q[0]*u*u+q[1]*u*r+q[2]*r*r+q[3]*u+q[4]*r+q[5]
def maximum(q,vertices):
    assert q[0]==1
    values=[value(q,*p)for p in vertices]
    for p,z in zip(vertices,vertices[1:]+vertices[:1]):
        mid=tuple((a+b)/2 for a,b in zip(p,z))
        y0=value(q,*p);yh=value(q,*mid);y1=value(q,*z)
        quadratic=2*(y1-2*yh+y0);linear=y1-y0-quadratic
        if quadratic<0:
            t=-linear/(2*quadratic)
            if 0<t<1:values.append(value(q,*(a+t*(b-a)for a,b in zip(p,z))))
    return max(values)


def main():
    prior=R.load(HERE/'joint_contact_full_718_independent.json')
    assert prior['status']=='PASS_INDEPENDENT_FULL_718_FIXED_LAMBDA_PAIRED_CONTACT_ASSEMBLY'
    linear=R.load(ROOT/'tonight/moments/linear_contact_bank/source_bank.json')
    linearproof=ROOT/'tonight/audit/linear_contact_bank_independent.json'
    assert R.load(linearproof)['input_sha256']==R.sha(ROOT/'tonight/moments/linear_contact_bank/source_bank.json')
    cases=[]
    for filename in('shared_polygon_actual_worst.json','shared_polygon_active_B0.json','shared_polygon_active_native_B0.json'):
        path=ROOT/'tonight/operator/contact_credit_joint'/filename;digest=R.sha(path);data=R.load(path)
        assert data['status']=='PASS_EXACT_SHARED_CORNER_SIMPLEX_POLYGON_PRICING_REQUIRES_INDEPENDENT_PROOF_AND_ASSEMBLY_AUDIT'
        envelopepath,envelope=bound(data['envelope_binding']);cell=data['source_cell']
        if cell.get('domain_override'):
            assert cell['physical_u_upper']==cell['u_interval'][1]and cell['physical_r_upper']==cell['r_interval'][1]
        else:assert cell==envelope['worst_cell']
        _,bank=bound(data['linear_bank_binding']);assert bank==linear
        b=F(data['operator_b']);assert b==R.B
        row=linear['rows'][data['linear_row']];a=F(row['a']);q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']))
        assert q[3]>=0 and data['contact_polynomial']==list(map(str,q))
        cert=data['certificate'];assert F(cert['a'])==a and cert['row']==data['linear_row']
        snapshot=data['source_curve_snapshot'];name=snapshot['name']
        assert name==data['profile']==cert['profile']
        ulo,uhi=map(F,cell['u_interval']);rlo,rhi=map(F,cell['r_interval']);klo,khi=map(F,cell['kappa_interval'])
        uhi=min(uhi,F(cell['physical_u_upper']));rhi=min(rhi,F(cell['physical_r_upper']))
        assert tuple(map(F,(snapshot['ulo'],snapshot['uhi'],snapshot['klo'],snapshot['khi'])))==(ulo,uhi,klo,khi)
        if snapshot['kind']=='native':
            mask=name[len('native:'):];binding=next(x for x in data['native_bindings']if x['mask']==mask)
            frame=N.frame(binding,b);rebuilt=native_curve(frame,binding,cell)
            assert F(snapshot['ell'])==F(rebuilt['ell'])and F(snapshot['beta'])==F(rebuilt['beta'])
            expected=[[v['corner'],v['constant'],v['linear'],v['quadratic'],*v['r_interval']]for v in rebuilt['curves']]
            assert snapshot['curves']==expected
            sourceproof=dict(kind='native',frame=binding)
        else:
            assert snapshot['kind']=='R3'and name=='R3:rank70'
            rr=R.bank()['rank70'];sourceproof=R.binding(rr);c=R.coefficients(rr,b)
            assert F(snapshot['ell'])==c['ell']and F(snapshot['beta'])==rr['p']['beta']
            assert len(snapshot['curves'])==4
            for corner,(k,z)in enumerate(product((klo,khi),(F(2,5),F(1)))):
                cc,Ic,Sc,g,lo,hi=snapshot['curves'][corner]
                assert cc==corner and F(Sc)==c['r_linear']and F(g)==c['r_quadratic']
                assert (F(lo),F(hi))==(rlo,rhi)
                assert R.price(rr,c,k,z).below(F(Ic))
        prices=list(map(F,cert['prices']));assert len(prices)==4
        for i,(k,s)in enumerate(product((klo,khi),(F(0),F(1)))):
            actual=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(N.K-F(633,15)*s)*I(k)*I(k).sqrt())
            assert actual.below(prices[i])
        pieces={}
        for corner,Ic,Sc,g,lo,hi in snapshot['curves']:
            piece=F(lo),F(hi),F(g);assert rlo<=piece[0]<piece[1]<=rhi
            pieces.setdefault(piece,{})[corner]=F(Ic),F(Sc)
        branches=[];ell=F(snapshot['ell'])
        for (lo,hi,g),C in pieces.items():
            assert set(C)==set(range(4))
            rectangle=[(F(-1),F(0),-ulo),(F(1),F(0),uhi),(F(0),F(-1),-lo),(F(0),F(1),hi)]
            # Enumerate all simplex vertices: four one-corner vertices and
            # each two-corner intersection with the active source plane.
            for i in range(4):
                poly=polygon(rectangle+[(a,F(1),prices[i])])
                if not poly:continue
                quad=[F(1),F(0),g,-2*ell,C[i][1],C[i][0]]
                branches.append(dict(type='corner',corners=[i],source_piece=[lo,hi],quadratic=quad,polygon=poly,cap=maximum(quad,poly)))
            for i,j in combinations(range(4),2):
                if prices[i]==prices[j]:continue
                if prices[i]>prices[j]:i,j=j,i
                D=prices[j]-prices[i];dI=C[j][0]-C[i][0];dS=C[j][1]-C[i][1]
                quad=[F(1),a*dS/D,g+dS/D,-2*ell+a*dI/D,C[i][1]+dI/D-prices[i]*dS/D,C[i][0]-prices[i]*dI/D]
                poly=polygon(rectangle+[(-a,F(-1),-prices[i]),(a,F(1),prices[j])])
                if not poly:continue
                branches.append(dict(type='two_corner_source_equality',corners=[i,j],source_piece=[lo,hi],quadratic=quad,polygon=poly,cap=maximum(quad,poly)))
        saved=cert['branches'];assert len(saved)==len(branches)
        for own,record in zip(branches,saved):
            for field in('type','corners'):assert own[field]==record[field]
            for field in('source_piece','quadratic'):assert own[field]==list(map(F,record[field]))
            assert set(own['polygon'])==set(tuple(map(F,p))for p in record['polygon'])
            assert own['cap']==F(record['cap'])
        cap=N.decimal_ceiling(max(x['cap']for x in branches),12)
        assert cap==F(data['cap'])and F(data['previous_cap'])==F(cell['cap'])
        assert R.sha(path)==digest
        cases.append(dict(source=R.file_record(path),source_envelope=R.file_record(envelopepath),
            profile=name,contact_row=data['linear_row'],branches=len(branches),
            nonzero_ur_branches=sum(q['quadratic'][1]!=0 for q in branches),
            cap=str(cap),baseline_cap=data['previous_cap'],actual_clipped_domain=dict(u=list(map(str,(ulo,uhi))),r=list(map(str,(rlo,rhi))),kappa=list(map(str,(klo,khi)))),
            accepted_source_primitives=sourceproof,
            scope='Only this exact clipped rectangle; complete physical-domain cover is conditional on the source envelope.') )
    out=dict(status='PASS_INDEPENDENT_SHARED_CORNER_SIMPLEX_THEOREM_AND_THREE_SELECTED_POLYGON_ASSEMBLIES',
        engine=selftest(),producer_modules_imported=False,accepted_linear_contact=R.file_record(linearproof),
        fresh_native_endpoint_prices=N.POINTS,cases=cases,
        theorem='Sequential separate convexity yields the same product moment-corner weights for profile and pure-linear contact prices. Relaxing to all nonnegative simplex weights constrained by sum(w_i*P_i)>=a*u+r only enlarges the objective. Vertices are single corners or two-corner intersections with the active source plane; equal-P pairs add no new extremum. Each mixed objective has u^2 coefficient one, so a nondegenerate polygon has no interior maximum. Exact rational hull vertices and every concave edge stationary point therefore certify its maximum. Point/segment degeneracies are included.',
        mathematical_cautions='The claimed source rectangle uses the saved physical_u_upper and physical_r_upper clips. Product-weight relaxation alone neither proves those clips nor the global source cover. The source row and profile remain fixed throughout each selected cell, and no convexity of their pointwise minimum is invoked.',
        scope='The simplex theorem and these three conditional selected source cells only. No new complete global norm bound or Grothendieck approximation follows without exhaustive base-domain acceptance.')
    (HERE/'shared_corner_polygon_independent.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],cases=[{k:c[k]for k in('profile','contact_row','branches','nonzero_ur_branches','cap')}for c in cases]),indent=2))


if __name__=='__main__':main()
