#!/usr/bin/env python3
"""Independent R3 corner prices and complete frozen-gate bindings.

No producer module is imported. Exact profile, reserve, transport, covariance,
PSD and tail algebra is reconstructed here; prices use our MPFR adapter.
Previously fresh compact-q replays are bound by profile and tree hashes.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import gzip
import hashlib
import json
import os
import sys

os.environ['GR_MPFR_BITS']='384'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest

B=F(7171,10000)
LOWER=ROOT/'tonight/operator/projection/lower_amp_certificate'
DIFFUSE=ROOT/'tonight/operator/projection/diffuse_certificate'
CORE=ROOT/'tonight/operator/projection/certificate'
RESULTS=ROOT/'tonight/operator/projection/results'
SOURCE=ROOT/'tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):
    if str(path).endswith('.gz'):
        with gzip.open(path,'rt')as stream:return json.load(stream)
    return json.loads(path.read_text())
def relative(path):return str(path.relative_to(ROOT))
def file_record(path):return dict(path=relative(path),sha256=sha(path))
def exact(raw):
    return {k:[F(x)for x in v]if isinstance(v,list)else F(v)for k,v in raw.items()}


def bank():
    rows=[]
    for label in ('68','70'):
        rows.append(dict(name='rank'+label,profile=LOWER/f'profile_u.{label}.json',
            kind='rank',reserve0=F(1,50000),reserve6=F(1,5000000),
            body=LOWER/f'body_u{label}_finer.json.gz',
            exterior=LOWER/f'exterior_u{label}_finer.json.gz',
            audit_root=HERE.parent/'lower_amp_R3',audit_label=label))
        rows.append(dict(name='full'+label,profile=LOWER/f'general_profile_u{label}_high_exact.json',
            kind='full',reserve0=F(1,1000),reserve6=F(1,1000000),
            body=LOWER/f'general_u{label}_high_body.json.gz',
            exterior=LOWER/f'general_u{label}_high_exterior.json.gz',
            whole_audit=ROOT/f'tonight/theory/general_r3_u{label}_whole_q_audit.json'))
    rows.append(dict(name='critical755',profile=CORE/'frozen_R3_profile_exact.json',
        kind='rank',reserve0=F(1,50000),reserve6=F(1,1000000),
        body=RESULTS/'directed_bernstein_region_0_8.json.gz',
        exterior=RESULTS/'directed_region_8_64.json.gz',tail=RESULTS/'directed_local_tail.json'))
    for name,label,audit_folder in (
        ('diffuse7478','interval_u7478_k0_01','diffuse_R3'),
        ('diffuse7505','interval_b7171_newworst','diffuse_b7171_R3')):
        rows.append(dict(name=name,profile=DIFFUSE/f'profile_{label}.json',kind='rank',
            reserve0=F(1,50000),reserve6=F(1,5000000),
            body=DIFFUSE/f'body_{label}.json.gz',exterior=DIFFUSE/f'exterior_{label}.json.gz',
            tail=DIFFUSE/f'tail_{label}.json',audit_root=HERE.parent/audit_folder,audit_label=label))
    for row in rows:
        raw=load(row['profile'])
        assert raw['kind']==('EXACT_RATIONAL_GENERAL_R3_PROFILE'if row['kind']=='full'else'EXACT_RATIONAL_R3_PROFILE')
        row['p']=exact(raw['profile'])
    return {row['name']:row for row in rows}


def psd(a,b,c):
    det=a*c-b*b
    assert a>=0 and c>=0 and det>=0
    return dict(diagonal=list(map(str,(a,c))),off_diagonal=str(b),determinant=str(det))


def spectral(row,b):
    p=row['p'];assert 0<b<=p['eps'] and p['eta']>0 and p['G']>0
    assert p['kappa']>=0 and p['p0']**2+2*p['p2']**2<1
    E00,E02,E22=even(row)
    matrices=[psd(E00,E02,E22),psd(E00+1,E02,E22-b/2),
              psd(E00-p['p0']**2,E02-p['p0']*p['p2'],E22-p['p2']**2),
              psd(p['R11'],p['R13'],p['R33']),
              psd(p['R11']+1,p['R13'],p['R33']-b/6),
              psd(p['R11']-b,p['R13'],p['R33'])]
    assert E22>0
    return dict(b=str(b),six_PSD_matrices=matrices)


def even(row):
    p=row['p']
    if row['kind']=='full':return p['E00'],p['E02'],p['E22']
    return p['p0']**2,p['p0']*p['p2'],p['p2']**2


def coefficients(row,b):
    p=row['p'];E00,E02,E22=even(row)
    theta=p['kappa']/(3*(p['kappa']+2))
    q0,q2,q4,q6=p['major'];q0+=row['reserve0'];q6+=row['reserve6']
    A=q6+p['zeta']**2;D=p['G']-p['d0']-p['zeta']**2
    assert A>0 and D>0
    yc=p['y']-E02*p['r']
    constant=(6*p['beta']**2/b+yc**2/E22+q0+q2+3*q4+15*q6
              +p['eta']*theta**2-2*p['nuisance2']*theta+6*p['zeta']**2)
    return dict(theta=theta,A=A,B=D,constant=constant,q4=q4,q6=q6,yc=yc,
                r_quadratic=b-E00+E02**2/E22,
                r_linear=-2*p['gamma']+2*yc*E02/E22,ell=p['ell'],b=b)


def price(row,c,k,z,u=F(0),r=F(0)):
    p=row['p'];theta=c['theta']
    linear=(c['q4']+15*c['q6']+p['eta']/2*((1+z)/9-2*theta/3+theta**2)
            +p['nuisance2']*(F(1,3)-theta)+9*p['zeta']**2-p['c'])
    H=I(c['A'])*(I(194)/I(10).sqrt()-I(F(633,9)*(z-F(2,5))))
    H+=I(F(9,4)*c['B'])*I(z)*I(z).sqrt()
    assert H.above(0)
    rational=(c['constant']+u*u-2*c['ell']*u+c['r_quadratic']*r*r+
              c['r_linear']*r+linear*k+c['B']*k*k)
    return I(rational)+I(k)*I(k).sqrt()*H


def tail(row,reserve0=None,reserve6=None):
    # Reconstruct exact half-line majorization; coefficients after q=64+s
    # being nonnegative prove the difference is nonnegative for every s>=0.
    p=row['p'];R=F(64);E00,E02,E22=even(row)
    q=list(p['major']);q[0]+=row['reserve0']if reserve0 is None else reserve0
    q[3]+=row['reserve6']if reserve6 is None else reserve6
    ap=sum(abs(v)for v in p['phi']);ad=sum(j*j*abs(v)for j,v in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2))
    derivative=ap/R**2+4*ad/R**4
    stein=theta*ap+ap/R**2+4*theta*ad/R**2+4*ad/R**4
    e0=abs(p['nuisance0']-p['nuisance2'])+derivative;e2=abs(p['nuisance2'])
    cc=p['c']-3*p['d0'];dd=p['d0'];h=p['zeta']**2/p['G']
    v=p['y']-E02*p['r']
    l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0']
    l2=v+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    a3=abs(p['zeta']*(-1+dd/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    trace=E00+2*E22+p['R11']+6*p['R33']+h*(p['p0']**2+2*p['p2']**2)
    constant=2*(abs(p['ell'])+F(5,2)*abs(p['beta']))+trace
    poly=[F(0)]*7
    for j,value in enumerate(q):poly[2*j]=value
    poly[0]-=constant+ap+stein+e0*e0/(2*p['eta'])
    poly[1]-=2*(abs(l1)+abs(l2));poly[2]-=e0*e2/p['eta']+cc*cc/(4*p['G'])
    poly[3]-=2*a3;poly[4]-=e2*e2/(2*p['eta'])+abs(cc*dd)/(2*p['G'])
    poly[6]-=dd*dd/(4*p['G'])
    shifted=[sum((poly[j]*comb(j,k)*R**(j-k)for j in range(k,7)),F(0))for k in range(7)]
    assert all(x>=0 for x in shifted)
    return dict(status='PASS_FRESH_EXACT_TAIL',R=str(R),difference_polynomial=list(map(str,poly)),
                shifted_coefficients=list(map(str,shifted)))


def match_tail(actual,record):
    assert F(actual['R'])==F(record['R'])==64
    for key in ('difference_polynomial','shifted_coefficients'):
        assert list(map(F,actual[key]))==list(map(F,record[key]))


def check_proof_hashes(report):
    for path,digest in report.get('shared_source_hashes',{}).items():assert sha(ROOT/path)==digest
    for path,digest in report.get('source_hashes',{}).items():assert sha(CORE/path)==digest


def binding(row):
    p=row['p'];name=row['name'];proofs=[];regions=[]
    if row['kind']=='full':
        audit=load(row['whole_audit']);assert audit['status']=='PASS_INDEPENDENT_GENERAL_R3_WHOLE_Q_PROFILE'
        assert audit['profile_sha256']==sha(row['profile']);check_proof_hashes(audit)
        proofs.append(file_record(row['whole_audit']))
    for part,qa,qb in (('body',F(0),F(8)),('exterior',F(8),F(64))):
        path=row[part];record=load(path);rp=exact(record['profile'])
        assert record['status'].startswith('PASS')and record['pending']==0
        assert tuple(map(F,(record['qa'],record['qb'])))==(qa,qb)
        assert set(rp)==set(p)
        for key in p:
            if key!='major':assert rp[key]==p[key],(name,part,key)
        dq0=rp['major'][0]-p['major'][0];dq6=rp['major'][3]-p['major'][3]
        assert rp['major'][1:3]==p['major'][1:3]
        assert 0<=dq0<=row['reserve0']and 0<=dq6<=row['reserve6']
        if name=='critical755':
            proof=ROOT/('tonight/theory/r3_bernstein_center_cover_audit.json'if part=='body'
                        else'tonight/theory/r3_local_partition_audit.json')
            report=load(proof);assert report['status']in(
                'PASS_FRESH_PHYSICAL_CENTER_PARTITION_AND_WHOLE_Q_GLUE',
                'PASS_RECORDED_COMPACT_PARTITION_AND_RATIONAL_TAIL')
            assert report['certificate_sha256']==sha(path);check_proof_hashes(report)
        elif row['kind']=='full':
            assert audit['certificate_sha256'][relative(path)]==sha(path)
            proof=row['whole_audit'];report=audit
        else:
            proof=row['audit_root']/f"independent_{row['audit_label']}_finer_{part}.json"
            report=load(proof);assert report['status']=='PASS_INDEPENDENT_PHYSICAL_Q_TREE'
            assert report['tree_sha256']==sha(path)and report['profile_sha256']==sha(row['profile'])
            assert report['accepted']==record['accepted']
            backend=ROOT/'tonight/recovery/profile_cover/directed_cubic_integral.py'
            assert report['cubic_integral_backend_sha256']==sha(backend)
        if file_record(proof)not in proofs:proofs.append(file_record(proof))
        regions.append({**file_record(path),'q':list(map(str,(qa,qb))),
                        'gate_reserve0':str(dq0),'gate_reserve6':str(dq6),
                        'target_minus_gate_reserve0':str(row['reserve0']-dq0),
                        'target_minus_gate_reserve6':str(row['reserve6']-dq6),
                        'accepted_leaves':record['accepted']})
    actual=tail(row)
    if name=='critical755':
        old=tail(row,F(1,1000000),F(1,1000000));match_tail(old,load(row['tail']))
        delta=row['reserve0']-F(1,1000000)
        assert F(actual['difference_polynomial'][0])-F(old['difference_polynomial'][0])==delta
        assert actual['difference_polynomial'][1:]==old['difference_polynomial'][1:]
        tail_binding=file_record(row['tail']);tail_binding['positive_constant_domination']=str(delta)
    elif row['kind']=='full':
        match_tail(actual,audit['tail']);tail_binding=file_record(row['whole_audit'])
    else:
        proof=row['audit_root']/f"independent_{row['audit_label']}_finer_tail.json"
        report=load(proof);assert report['status']=='PASS_INDEPENDENT_EXACT_R3_TAIL'
        match_tail(actual,report);proofs.append(file_record(proof))
        if name.startswith('rank'):
            ledger=LOWER/'lower_source_domain_finer.json'
            record=next(x['tail']for x in load(ledger)['cells']if x['profile']==name[4:])
            tail_binding={**file_record(ledger),'tail_profile':name[4:]}
        else:record=load(row['tail']);tail_binding=file_record(row['tail'])
        match_tail(actual,record)
    return dict(profile=file_record(row['profile']),kind=row['kind'],
        target_reserves=list(map(str,(row['reserve0'],row['reserve6']))),
        q_cover=regions,tail_binding=tail_binding,fresh_reconstructed_tail=actual,
        accepted_independent_proofs=proofs,
        original_spectral_gates=spectral(row,p['eps']),transported_spectral_gates=spectral(row,B),
        scope='The scalar local gate has frozen theta and no actual kappa; its accepted whole-q inequality may be priced at every kappa>=0. Actual moment inputs occur only in the separately convex outer price.')


def main():
    initial_sha=sha(SOURCE);source=load(SOURCE);rows=bank()
    assert F(source['operator_b'])==B
    checked_ledger=[]
    for item in source['accepted_R3_binding']:
        path=Path(item['path'])
        if not path.is_absolute():path=ROOT/path
        assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))
    bindings={name:binding(row)for name,row in rows.items()}
    coefficients_record={}
    for name,row in rows.items():
        c=coefficients(row,B);row['coefficients']=c
        expected={k:str(c[k])for k in ('r_quadratic','r_linear','ell')}
        assert source['profile_coefficients'][name]==expected
        assert c['theta']==(F(1,123)if name.startswith('diffuse')else F(53,759))
        for r in (F(0),F(1,3),F(1)):
            p=row['p'];E00,E02,E22=even(row);y=p['y']+E02*(r-p['r'])
            assert y-E02*r==c['yc']
            assert ((B-E00)*r*r-2*p['gamma']*r+y*y/E22==
                    c['r_quadratic']*r*r+c['r_linear']*r+c['yc']**2/E22)
        baseline=coefficients(row,row['p']['eps'])
        assert c['constant']-baseline['constant']==6*row['p']['beta']**2*(1/B-1/row['p']['eps'])
        coefficients_record[name]={k:str(v)for k,v in c.items()}
    corners=[];keys=set();intervals={name:set()for name in rows}
    for item in source['cached_corner_constants']:
        name=item['profile'];row=rows[name];a,b=map(F,item['kappa_interval']);saved=F(item['upper'])
        assert 0<=a<b<=90 and (name,a,b)not in keys
        keys.add((name,a,b));intervals[name].add((a,b));values=[]
        for k,z in product((a,b),(F(2,5),F(1))):
            actual=price(row,row['coefficients'],k,z)
            assert actual.below(saved),(name,a,b,k,z,saved,actual.record())
            rounded=F(actual.hi.text(2,80));assert rounded<=saved
            covariance=k/3-row['coefficients']['theta']*(k+2)
            values.append(dict(kappa=str(k),shape=str(z),price=actual.record(),
                independent_upper_80_significant_digits=str(rounded),
                plus_covariance_correction=str(row['p']['nuisance2']*covariance)))
        fresh=max(F(v['independent_upper_80_significant_digits'])for v in values)
        corners.append(dict(profile=name,kappa_interval=list(map(str,(a,b))),saved_upper=str(saved),
                            fresh_upper=str(fresh),saved_minus_fresh=str(saved-fresh),corners=values))
    assert len(keys)==98 and len(rows)==7
    for name,leaves in intervals.items():
        previous=F(0);assert len(leaves)==14
        for a,b in sorted(leaves):assert a==previous;previous=b
        assert previous==90
    assert all(v==intervals['rank68']for v in intervals.values())
    used=Counter(c['winner']for c in source['cells']if c['winner'].startswith('R3:'))
    assert all(name[3:]in rows for name in used)
    assert sha(SOURCE)==initial_sha
    output=dict(status='PASS_INDEPENDENT_384BIT_ALL_R3_CORNER_PRICES_AND_WHOLE_Q_BINDINGS',
        source_path=relative(SOURCE),source_sha256=initial_sha,operator_b=str(B),
        engine=selftest(),independent_script=file_record(Path(__file__)),
        independent_MPFR_backend=file_record(HERE.parent/'mpfr_interval.py'),
        profiles=7,corner_tables=98,directed_corner_prices=392,
        kappa_interval=['0','90'],shape_interval=['2/5','1'],
        checked_producer_binding_ledger=checked_ledger,complete_independent_bindings=bindings,
        expanded_coefficients=coefficients_record,corner_tables_replayed=corners,
        R3_winner_counts=dict(used),all_used_R3_profiles_whole_q_bound=True,
        covariance_sign='The nuisance2 covariance correction is ADDED: nuisance2*(kappa/3-theta*(kappa+2)).',
        convexity='A>0,B>0 and K(z)>0. The kappa second derivative is 2B+3H(z)/(4sqrt(kappa))>=0; the shape second derivative is 27B*kappa^(3/2)/(16sqrt(z))>=0. Continuity covers kappa=0.',
        transport='y(r)=y0+E02*(r-r0) fixes y-E02*r exactly. Lowering b only relaxes the spectral gates; local q comparison is unchanged. Outer change is (b-b0)*r^2+6beta^2*(1/b-1/b0).',
        no_source_credit_used=True,
        scope='Accepts every saved R3 corner upper bound and binds all seven gates on the whole scalar q line. Finite envelope min/quadratic/source-domain union and lower/high amplitude branches require the separate independent audits; no standalone global norm claim.')
    (HERE/'R3_corner_replay.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k]for k in ('status','profiles','corner_tables','directed_corner_prices','R3_winner_counts','source_sha256')},indent=2))


if __name__=='__main__':main()
