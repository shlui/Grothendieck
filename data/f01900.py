#!/usr/bin/env python3
"""Independent native/contact moment coupling with one fixed cell multiplier."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from bisect import bisect_right
import hashlib
import json
import os
import sys

os.environ['GR_MPFR_BITS']='768'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
import replay_native_R3_unions as N
from mpfr_interval import I,selftest

SOURCE=ROOT/'tonight/operator/contact_credit_joint/worst_cell_old_frames.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def record(path):return dict(path=str(path.relative_to(ROOT)),sha256=sha(path))


def qmaximum(c,l,g,a,b):
    choices=[(a,c+l*a+g*a*a),(b,c+l*b+g*b*b)]
    if g<0:
        x=-l/(2*g)
        if a<x<b:choices.append((x,c+l*x+g*x*x))
    x,v=max(choices,key=lambda t:t[1])
    return v,x


def main():
    digest=sha(SOURCE);source=read(SOURCE)
    assert source['status']=='PASS_EXACT_FIXED_MULTIPLIER_PAIRED_CORNER_REPRICING_REQUIRES_INDEPENDENT_ASSEMBLY_AUDIT'
    assert source['worst_only']and len(source['cells'])==1
    cell=source['cells'][0];assert cell==source['worst_cell']
    b=F(source['operator_b']);assert b==F(7171,10000)
    binding=next(x for x in source['native_bindings']if x['mask']==cell['frame'])
    frame=N.frame(binding,b);package=N.package(binding['package'])
    assert sha(ROOT/binding['mask'])==binding['mask_sha256']
    assert package['source_sha256']==binding['source_sha256']
    primary=package['profile'];beta=F(primary['beta']);assert beta>0
    credit=cell['joint_contact_credit'];lam=F(credit['lambda_value'])
    assert lam==F(20885397,10**9)>=0
    contact_path=ROOT/source['contact_binding']['path'];assert sha(contact_path)==source['contact_binding']['sha256']
    contact=read(contact_path);assert contact['status']=='PASS_EXACT_ALL_HALFLINE_DIRECT_CUBIC_HEAD_CONTACT_MAJORANT'
    accepted_path=ROOT/'tonight/audit/cubic_head_contact_independent.json';accepted=read(accepted_path)
    assert accepted['status']=='PASS_INDEPENDENT_ALL_LINE_CUBIC_HEAD_CONTACT_GATE'
    assert accepted['input_sha256']==sha(contact_path)
    scale=F(credit['contact_scale']);assert scale>0
    row=contact['rows'][credit['contact_row']];assert F(row['w'])==F(credit['contact_w'])==2
    q=[F(v)*scale**(1-2*j)for j,v in enumerate(row['polynomial_coefficients_q0_q2_q4_q6'])]
    assert q[3]>0
    A=F(contact['parameters']['A'])*scale;B=F(contact['parameters']['B'])*scale
    a=A+3*B;d2=6*B*B
    assert F(contact['parameters']['a'])*scale==a
    beta_term=I(beta*B/b)*I(6).sqrt()
    saved_beta=F(credit['beta_d_over_b_upper']);assert beta_term.below(saved_beta)
    ulo,uhi=map(F,cell['u_interval']);rlo,rhi=map(F,cell['r_interval']);klo,khi=map(F,cell['kappa_interval'])
    assert (ulo,uhi,rlo,rhi,klo,khi)==tuple(map(F,('.745','.75','.266','.268','.8','1')))
    assert F(cell['physical_u_upper'])==uhi and F(cell['physical_r_upper'])==rhi
    order=list(product((klo,khi),(F(0),F(1))))
    assert credit['corner_order']==['klo,s0','klo,s1','khi,s0','khi,s1']
    saved_contact=list(map(F,credit['contact_corner_prices']));assert len(saved_contact)==4
    contact_values=[]
    for saved,(k,s)in zip(saved_contact,order):
        value=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(N.K-F(633,15)*s)*I(k)*I(k).sqrt())
        assert value.below(saved)
        contact_values.append(dict(kappa=str(k),shape=str(s),price=value.record(),saved_upper=str(saved)))
    rg=frame['rg'];ep=frame['endpoints'][klo,khi]
    splits=[rlo]+[r for r in rg if rlo<r<rhi]+[rhi]
    outer_values=[(u,u*u-2*frame['ell']*u-lam*a*u)for u in (ulo,uhi)]
    outer_u,outer=max(outer_values,key=lambda t:t[1])
    lambda_square=lam*lam*d2/(4*b)
    endpoint_checks=[];curves=[];combined=[]
    mask=read(ROOT/binding['mask'])
    band=next(x for x in mask['endpoint_prices']if tuple(map(F,x['kappa_interval']))==(klo,khi))
    for left,right in zip(splits,splits[1:]):
        j=min(bisect_right(rg,left)-1,len(rg)-2)
        assert rg[j]<=left<right<=rg[j+1]
        for endpoint in (j,j+1):
            choice=band['selected_rows'][endpoint];assert choice!='r=0'
            name=choice['source_package'];local=N.package(name)
            scalar=N.scalar(name,choice['native_row']);assert scalar['r']==rg[endpoint]
            qq=scalar['q'][:]
            if choice['l1_row']!='constant':
                pure,eta=N.pure(choice['l1_row']);qq=[v+scalar['mu']*w for v,w in zip(qq,pure)]
            assert qq[3]>=0
            for corner,(k,s)in enumerate(order):
                actual=N.price(name,choice['native_row'],k,s,F(band['rho']),choice['l1_row'])
                actual+=local['inflation']+local['tail']*F(849,1000)*rg[endpoint]
                actual+=F(mask['local_constant_transport'])
                assert actual.below(ep[endpoint][corner])
                endpoint_checks.append(dict(source_r=str(rg[endpoint]),corner=corner,kappa=str(k),shape=str(s),
                    source_package=name,native_row=choice['native_row'],pure_L1_row=choice['l1_row'],
                    effective_sixth_coefficient=str(qq[3]),native_price=actual.record(),saved_upper=str(ep[endpoint][corner])))
        for corner,(k,s)in enumerate(order):
            slope=(ep[j+1][corner]-ep[j][corner])/(rg[j+1]-rg[j])
            intercept=ep[j][corner]-slope*rg[j]
            constant=intercept+lam*(saved_contact[corner]+saved_beta)+lambda_square
            linear=slope-lam
            value,active_r=qmaximum(constant,linear,frame['g'],left,right)
            value+=outer
            curves.append(dict(corner=corner,kappa=str(k),shape=str(s),r_interval=list(map(str,(left,right))),
                native_constant=str(intercept),native_linear=str(slope),quadratic=str(frame['g']),
                fixed_lambda=str(lam),contact_corner_upper=str(saved_contact[corner]),
                augmented_constant=str(constant),augmented_linear=str(linear),
                maximizing_r=str(active_r),price=str(value)))
            combined.append(value)
    maximum=max(combined);cap=N.decimal_ceiling(maximum,12)
    assert maximum==F(credit['exact_unrounded_cap'])
    assert cap==F(cell['cap'])==F(source['maximum_cell_cap'])==F('0.565027324915')
    assert cap<F(cell['previous_cap'])
    assert sha(SOURCE)==digest
    output=dict(status='PASS_INDEPENDENT_SELECTED_CELL_PAIRED_NATIVE_CONTACT_CREDIT',
        source=record(SOURCE),engine=selftest(),producer_modules_imported=False,
        native_frame=record(ROOT/binding['mask']),native_package=binding['package'],
        accepted_contact_gate=record(accepted_path),contact_input=record(contact_path),
        operator_b=str(b),cell=dict(u=cell['u_interval'],r=cell['r_interval'],kappa=cell['kappa_interval'],shape=['0','1']),
        fixed_contact_row=credit['contact_row'],fixed_contact_w=credit['contact_w'],fixed_scale=str(scale),fixed_lambda=str(lam),
        native_beta=str(beta),contact_parameters=dict(A=str(A),B=str(B),a=str(a),d_squared=str(d2),q=list(map(str,q))),
        beta_d_over_b=dict(price=beta_term.record(),saved_upper=str(saved_beta)),
        contact_corner_prices=contact_values,selected_native_endpoint_checks=endpoint_checks,
        native_gram_coefficient=str(frame['g']),paired_source_curves=curves,
        amplitude_maximizer=str(outer_u),exact_unrounded_cap=str(maximum),rounded_cap=str(cap),
        improvement=str(F(cell['previous_cap'])-cap),
        coupling='One contact row/scale and one nonnegative rational lambda are fixed across the entire cell. Each native moment corner is added to the contact price at the SAME kappa and shape before maximizing; independent component maxima are not substituted.',
        completion='Adding lambda*(Aq-a*u-d*v-r) before -b*v^2-2beta*v completion gives PLUSlambda*beta*d/b +lambda^2*d^2/(4b), together with -lambda*a*u-lambda*r.',
        convexity='Every selected native endpoint has nonnegative effective sixth coefficient. Adding lambda times the positive contact sixth coefficient preserves convexity in kappa and affine dependence on shape. Convexity in u gives endpoint maximization. Exact quadratic maximization in r includes an interior vertex when negative curvature requires it. No joint convexity of (kappa,shape), or convexity of the pointwise minimum of frames, is assumed.',
        scope='This selected continuous source cell only, conditional on the accepted native and contact scalar primitives. A complete joint-credit envelope needs separate exhaustive coverage and every selected credit audit. Frozen original bundle unchanged.')
    (HERE/'joint_contact_selected_cell_independent.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k]for k in ('status','fixed_lambda','rounded_cap','improvement','native_gram_coefficient')},indent=2))


if __name__=='__main__':main()
