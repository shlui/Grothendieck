#!/usr/bin/env python3
"""Fresh whole-line scalar acceptance for the continuous .30--.31 profile.

No discovery code is imported. Exact decimal inputs are frozen as Fractions.
Compact pieces use outward Arb Taylor bounds. Tails use exact polynomial
translation. The clipped Young target is uncentered, with U >= -1.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import json
import time
from flint import arb, arb_series, fmpq, fmpq_poly, ctx

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1]/'tonight/recovery/profile_cover/continuous_e7200_u3000_3100_r6_d10.json'
HINGES = tuple(map(F, ('0.5','1','1.5','2','2.25','2.5','2.75','3',
                        '3.25','3.4','3.6','3.75','4','4.5','5','6','7','8','10','12')))
Q0_SLACK, Q6_SLACK = F(1,500), F(1,20_000_000)


def aa(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def qq(x):
    x = F(x)
    return fmpq(x.numerator, x.denominator)


def pp(coefficients):
    return fmpq_poly([qq(x) for x in coefficients])


def horner(coefficients, x):
    value = arb_series([0])
    for coefficient in reversed(coefficients):
        value = value*x+coefficient
    return value


@lru_cache(None)
def mills(x):
    x = aa(x)
    return (arb.pi()/2).sqrt()*(x*x/2).exp()*(x/arb(2).sqrt()).erfc()


@lru_cache(None)
def hinge_constants(b):
    bb = aa(b)
    g = 1-bb*mills(b)
    h = g*(-bb*bb/2).exp()
    c = (2/arb.pi()).sqrt()
    return c*h, h


def rational_bounds(value, digits=60):
    scale = 10**digits
    lo = F(str(value.lower().fmpq()))
    hi = F(str(value.upper().fmpq()))
    return F(lo.numerator*scale//lo.denominator, scale), F(-((-hi.numerator*scale)//hi.denominator), scale)


def base_polynomials(values, degree=10, length=F(2)):
    x = pp([0,1])
    den, z = x+qq(length), x-qq(length)
    cheb = [pp([1]), z]
    for _ in range(2,degree+1):
        cheb.append(2*z*cheb[-1]-den*den*cheb[-2])
    numerator = sum((qq(values[j]/(j+1))*x**(j+1) for j in range(4)), pp([0]))*den**(degree+1)
    for j, coefficient in enumerate(values[4:degree+5]):
        numerator += qq(coefficient)*x*cheb[j]*den**(degree-j)
    deriv_numerator = numerator.derivative()*den-(degree+1)*numerator
    return numerator, deriv_numerator, den**(degree+2)


def freeze(data, index):
    raw = data['rows'][index]
    p = data['profile']
    r, mu, eta = F(str(raw['r'])), F(str(raw['lambda_L1'])), F(str(raw['eta']))
    assert r > 0 and mu >= 0 and eta > 0
    t = tuple(F(str(x)) for x in data['t'])
    heights = tuple(F(str(x)) for x in p['nodal_majorant'])
    tau = F(str(p['tail_slope']))
    slopes = [(heights[i+1]-heights[i])/(t[i+1]-t[i]) for i in range(len(t)-1)]
    weights = [r*slopes[0]] + [r*(slopes[i]-slopes[i-1]) for i in range(1,len(slopes))] + [r*(tau-slopes[-1])]
    target = tuple(zip((x/r for x in t),weights))
    coefficients = tuple(F(str(x)) for x in raw['phi_coefficients'])
    assert len(coefficients) == 35
    regression = tuple(zip(HINGES,coefficients[15:]))
    q = [F(str(x)) for x in raw['majorant']]
    q[0] += Q0_SLACK
    q[3] += Q6_SLACK
    assert q[3] > 0
    numerator, deriv_numerator, denominator = base_polynomials(coefficients)
    row = dict(index=index,r=r,mu=mu,eta=eta,t=t,heights=heights,tau=tau,
               target=target,regression=regression,q=q,coefficients=coefficients,
               N=numerator,B=deriv_numerator,D=denominator,
               N_arb=[aa(str(x)) for x in numerator.coeffs()],
               B_arb=[aa(str(x)) for x in deriv_numerator.coeffs()],
               kappa_interval=list(map(lambda x:F(str(x)),raw['kappa_interval'])))
    row['target_mean'] = aa(heights[0])+sum((aa(w)*hinge_constants(b)[0] for b,w in target),arb(0))
    row['target_u'] = sum((aa(w)*(aa(b)+hinge_constants(b)[0]) for b,w in target),arb(0))-aa(mu)*(2/arb.pi()).sqrt()
    row['phi_u'] = sum((aa(w)*(aa(b)+hinge_constants(b)[0]) for b,w in regression),arb(0))
    row['phi_v'] = sum((w for _,w in regression),F(0))
    return row


def segment_constants(row,left):
    phi_u = sum((aa(w)*(hinge_constants(b)[0]+(aa(b) if b<=left else arb(0)))
                 for b,w in row['regression']),arb(0))
    target_u = sum((aa(w)*(hinge_constants(b)[0]+(aa(b) if b<=left else arb(0)))
                   for b,w in row['target']),arb(0))
    exponent = sum((aa(w)*hinge_constants(b)[1] for b,w in row['regression'] if b>left),arb(0))
    exponent -= sum((aa(w)*hinge_constants(b)[1] for b,w in row['target'] if b>left),arb(0))
    residual_g = phi_u-target_u+aa(row['mu'])*(2/arb.pi()).sqrt()
    phi_v = sum((w for b,w in row['regression'] if b<=left),F(0))
    return phi_u,residual_g,exponent,phi_v


def jets(row,pos,r_mills,constants):
    x = arb_series([pos,1])
    den = x+2
    phi = horner(row['N_arb'],x)/den**11
    derivative = horner(row['B_arb'],x)/den**12
    rr = [r_mills,pos*r_mills-1]
    for j in range(1,4):
        rr.append((pos*rr[-1]+rr[-2])/(j+1))
    g = 1-x*arb_series(rr)
    phi_u,rg,exponent,phi_v = constants
    residual = rg*g-derivative+x*(x*x/2).exp()*exponent
    regression = derivative-x*phi+aa(phi_v)*x-phi_u
    majorant = sum((aa(q)*x**(2*j) for j,q in enumerate(row['q'])),arb_series([0]))
    return majorant,regression,residual


def compact(row,cutoff):
    breaks = sorted({F(0),cutoff} | {b for b,_ in row['target'] if 0<b<cutoff}
                    | {b for b,_ in row['regression'] if 0<b<cutoff})
    leaves,minimum,linear,quadratic = [],None,0,0
    eta = aa(row['eta'])
    for left,right in zip(breaks[:-1],breaks[1:]):
        constants = segment_constants(row,left)
        stack = [(left,right,0)]
        while stack:
            lo,hi,depth = stack.pop()
            mid,rad = (lo+hi)/2,(hi-lo)/2
            ctx.cap = 5
            try:
                whole = jets(row,aa(lo).union(aa(hi)),mills(lo).union(mills(hi)),constants)
                center = jets(row,aa(mid),mills(mid),constants)
                clipped = whole[2][0] >= eta
                def gap(parts):
                    majorant,regression,residual = parts
                    penalty = residual-eta/2 if clipped else residual*residual/(2*eta)
                    return majorant-regression-penalty
                cc,ww = gap(center),gap(whole)
                lower = cc[0]-sum((abs(cc[j])*aa(rad**j) for j in range(1,4)),arb(0))-abs(ww[4])*aa(rad**4)
            except (ValueError,ZeroDivisionError):
                lower,clipped = arb(-1),False
            if lower > 0:
                bound = F(str(lower.lower().fmpq()))
                minimum = bound if minimum is None else min(minimum,bound)
                leaves.append([str(lo),str(hi),depth,'linear' if clipped else 'quadratic'])
                if clipped:linear+=1
                else:quadratic+=1
            else:
                assert depth<45,('COMPACT_FAIL',row['index'],str(lo),str(hi),str(lower))
                stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
    return dict(domain=['0',str(cutoff)],leaves=leaves,leaf_count=len(leaves),
                minimum_directed_gap=str(minimum),linear_leaves=linear,quadratic_leaves=quadratic)


def positive_translate(poly,cutoff):
    shifted = poly(pp([cutoff,1]))
    coefficients = shifted.coeffs()
    return bool(coefficients and coefficients[0]>0 and all(x>=0 for x in coefficients)),coefficients


def tail(row):
    x = pp([0,1]);D,N,B = row['D'],row['N'],row['B'];eta=row['eta']
    # The algebra below is the all-active hinge formula. Its unbounded
    # domain must begin after every target and regression breakpoint.
    last_hinge = max(b for b,_ in row['target']+row['regression'])
    lower_u,upper_u = rational_bounds(row['phi_u'])
    lower_c,upper_c = rational_bounds(row['phi_u']-row['target_u'])
    lower_c,upper_c = min(lower_c,F(0)),max(upper_c,F(0))
    q = pp([row['q'][0],0,row['q'][1],0,row['q'][2],0,row['q'][3]])
    v = row['phi_v']
    selector = -x*x*B+qq(lower_c)*D-qq(eta)*x*x*D
    linear_gap = x*x*((q-qq(v)*x+qq(lower_u+eta/2))*D+x*N*(x+2))-qq(upper_c)*D
    base = q*D-B+x*N*(x+2)-qq(v)*x*D+qq(lower_u)*D
    quadratic_gaps = [qq(2*eta)*x**4*D*base-(-x*x*B+qq(c)*D)**2
                      for c in (lower_c,upper_c)]
    for cutoff in map(F,(64,128,256,512,1024,2048,4096)):
        if cutoff < last_hinge:
            continue
        choose_linear,sc = positive_translate(selector,cutoff)
        candidates = [linear_gap] if choose_linear else quadratic_gaps
        proofs = [positive_translate(poly,cutoff) for poly in candidates]
        if all(ok for ok,_ in proofs):
            return cutoff,dict(domain=[str(cutoff),'infinity'],branch='linear' if choose_linear else 'quadratic_endpoints',
                last_hinge=str(last_hinge),all_hinges_active=True,
                phi_u_interval=[str(lower_u),str(upper_u)],
                residual_g_interval=[str(lower_c),str(upper_c)],
                selector_minimum=str(min(sc)) if choose_linear else None,
                checks=[dict(degree=poly.degree(),translated_coefficients=len(cc),minimum_coefficient=str(min(cc)),
                             polynomial_sha256=hashlib.sha256('|'.join(map(str,poly.coeffs())).encode()).hexdigest())
                        for poly,(_,cc) in zip(candidates,proofs)])
    raise AssertionError(('TAIL_FAIL',row['index']))


def prices(row):
    c=(2/arb.pi()).sqrt();big=194/arb(10).sqrt();rho=F(849,1000)
    records=[]
    for kappa in row['kappa_interval']:
        for shape in (0,1):
            sixth=15+15*aa(kappa)+(big-aa(F(633*shape,15)))*aa(kappa)*aa(kappa).sqrt()
            q=row['q']
            price=row['target_mean']+aa(q[0]+q[1]+q[2]*(3+kappa))+aa(q[3])*sixth
            price+=aa(row['eta']*kappa*(F(7,90)+F(shape,30)))
            price+=aa(row['mu'])*(aa(rho)-c)
            records.append(dict(kappa=str(kappa),shape=shape,price_upper=str(F(str(price.upper().fmpq()))),
                                price_display=float(price)))
    return records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--rows',type=int,nargs='*')
    parser.add_argument('--bits',type=int,default=512)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    args=parser.parse_args();ctx.prec=args.bits;ctx.cap=5
    args.output_dir.mkdir(parents=True,exist_ok=True)
    data=json.loads(args.source.read_text())
    indices=args.rows if args.rows is not None else list(range(len(data['rows'])))
    summary=[];start=time.time()
    for index in indices:
        row_start=time.time();row=freeze(data,index);cutoff,tail_proof=tail(row)
        compact_proof=compact(row,cutoff)
        frozen=dict(r=str(row['r']),lambda_L1=str(row['mu']),eta=str(row['eta']),
                    kappa_interval=list(map(str,row['kappa_interval'])),majorant=list(map(str,row['q'])),
                    phi_coefficients=list(map(str,row['coefficients'])),degree=10,scale='2',
                    regression_hinges=list(map(str,HINGES)),added_q0=str(Q0_SLACK),added_q6=str(Q6_SLACK))
        result=dict(status='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE',row_index=index,precision_bits=args.bits,
                    frozen_row=frozen,compact=compact_proof,tail=tail_proof,prices=prices(row),
                    seconds=time.time()-row_start,
                    source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    scope='Whole-line scalar majorant for the exact frozen original profile; local inflation/tail change and operator assembly are separate.')
        out=args.output_dir/f'continuous_scalar_row_{index:02d}.json';out.write_text(json.dumps(result,indent=2)+'\n')
        item=dict(row=index,status=result['status'],leaves=compact_proof['leaf_count'],tail_cutoff=str(cutoff),
                  largest_corner=max(p['price_display'] for p in result['prices']),seconds=result['seconds'],file=out.name)
        summary.append(item);print(json.dumps(item),flush=True)
        (args.output_dir/'continuous_scalar_summary.json').write_text(json.dumps(dict(status='PASS' if len(summary)==len(data['rows']) else 'PARTIAL',
            source=str(args.source),rows=summary,seconds=time.time()-start),indent=2)+'\n')


if __name__=='__main__':
    main()
