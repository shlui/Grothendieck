#!/usr/bin/env python3
"""Bind the accepted global cap to the two independently checked finite tables."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    here=Path(__file__).resolve().parent
    norm=here.parent/'audit'/'FINAL_GOAL_AFFINE_MOMENT_GLOBAL_ACCEPTANCE.json'
    arb=here/'FINAL_GOAL_LOG_CUTOFF_FINITE_ARB768.json'
    mpfr=here/'FINAL_GOAL_LOG_CUTOFF_FINITE_MPFR_AUDIT.json'
    rational=here/'FINAL_GOAL_RATIONAL_FINITE_BOUND.json'
    proof=json.loads(norm.read_text());a=json.loads(arb.read_text())
    m=json.loads(mpfr.read_text());r=json.loads(rational.read_text())
    assert proof['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    assert proof['target_1775_reached'] and F(proof['operator_b'])==F(711,1000)
    assert F(proof['norm_upper_bound'])==F(40,71)
    assert F(proof['constructive_G_R_lower_bound'])==F(71,40)
    assert F(a['operator_b'])==F(r['operator_b'])==F(m['operator_b'])==F(711,1000)
    assert F(a['required_global_norm_cap'])==F(r['scalar_cap_input'])==F(m['norm_cap'])==F(40,71)
    assert a['proof_binding']['sha256']==r['norm_proof']['sha256']==sha(norm)
    assert m['status']=='PASS_INDEPENDENT_MPFR768_LOG_CUTOFF_ARITHMETIC'
    assert m['source']['sha256']==sha(arb) and a['dimension_power']==4
    assert len(a['rows'])==len(m['rows'])==6
    for x,y in zip(a['rows'],m['rows']):
        assert x['n']==y['n'] and x['dimension']==x['n']**4==y['d']
        assert x['finite_ratio_decimal']==y['root_floor']
        assert 0<F(x['finite_ratio_decimal'])<F(71,40)
    assert r['status']=='PASS_EXACT_RATIONAL_FINITE_ARITHMETIC'
    assert F(r['endpoint_lower'])==F(71,40)
    result={'status':'PASS_COMPLETE_GOAL_NORM_AND_FINITE_ARITHMETIC_BINDING',
        'operator_b':'711/1000','norm_upper_bound':'40/71','G_R_lower_bound':'71/40',
        'directed_power4_rows':[{k:x[k]for k in('n','dimension','integer_cutoff','finite_ratio_decimal')}for x in a['rows']],
        'exact_rational_rows':[{k:x[k]for k in('dimension','ratio_lower_floor12')}for x in r['finite_bounds']],
        'bindings':{k:{'path':str(p),'sha256':sha(p)}for k,p in
                    [('accepted_global_norm',norm),('Arb768_table',arb),('independent_MPFR768_audit',mpfr),('exact_Fraction_table',rational)]},
        'scope':'Certified finite lower formulas tend to71/40. The actual reduced matrix ratio tends to1/N(.711), at least71/40; no rate or equality for the actual denominator is asserted.',
        'entrypoint_sha256':sha(__file__)}
    output=here/'FINAL_GOAL_FINITE_FORMULA_ACCEPTANCE.json'
    assert not output.exists()
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'output':str(output),'sha256':sha(output)},indent=2))


if __name__=='__main__':
    main()
