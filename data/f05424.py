#!/usr/bin/env python3
"""Portable optimizer-free replay of the accepted strengthened global cap."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'replay_outputs'


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(work,tag,args,bits=768):
    env=dict(os.environ);env['GR_MPFR_BITS']=str(bits)
    env['PYTHONPATH']=os.pathsep.join(str(work/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
    with(OUT/(tag+'.log')).open('w')as stream:
        result=subprocess.run([sys.executable,*map(str,args)],cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'{tag} failed; inspect replay_outputs/{tag}.log')
    print(json.dumps({'stage':'PASS','job':tag}),flush=True)


def refresh_metadata_binding(path,key):
    data=read(path)
    for value in data[key].values():value['sha256']=sha(path.parents[3]/value['path'])
    path.write_text(json.dumps(data,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args();start=time.time()
    OUT.mkdir(exist_ok=True)
    manifest=read(ROOT/'MANIFEST.json')
    assert all(sha(ROOT/p)==h for p,h in manifest['sha256'].items())
    provenance=read(ROOT/'PROVENANCE.json');strong=provenance['strengthened']
    with gzip.open(ROOT/strong['archive'],'rb')as stream:original=stream.read()
    assert hashlib.sha256(original).hexdigest()==strong['original_envelope_sha256']
    envelope='tonight/recovery/profile_cover/global_middle_b7171_critical_ninebank_contact_envelope.json'
    portable=read(ROOT/envelope);original_record=json.loads(original)
    assert {k:v for k,v in original_record.items()if k!='accepted_R3_binding'}=={
           k:v for k,v in portable.items()if k!='accepted_R3_binding'}
    # The copied original entrypoint replays every old branch, native price,
    # seven-gate binding and exact partition. It leaves the evidence immutable.
    command=[ROOT/'verify_original_components.py']
    if args.full:command+=['--full','--workers',str(args.workers)]
    run(ROOT,'original_prerequisites_and_prices',command)
    with tempfile.TemporaryDirectory(prefix='grothendieck_strengthened_replay_')as temporary:
        work=Path(temporary)/'proof';shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('replay_outputs','__pycache__'))
        if args.full:
            for part in('body','exterior','tail'):
                run(work,'targeted_lowr_'+part,[work/'tonight/audit/replay_R3_generous_q_trees.py',
                    '--profile','interval_lowr_u76_r16_k0_2','--branch','generous','--part',part,
                    '--data-dir',work/'tonight/operator/projection/targeted_certificate',
                    '--output-dir',work/'tonight/audit/targeted_lowr_R3_recheck'],384)
            acceptance=work/'tonight/audit/targeted_lowr_R3_recheck/ACCEPTANCE.json'
            data=read(acceptance)
            for value in data['reports'].values():value['sha256']=sha(work/value['path'])
            acceptance.write_text(json.dumps(data,indent=2)+'\n')
            run(work,'targeted_highk_whole_q',[work/'tonight/theory/replay_targeted_highk.py','--workers',str(args.workers)],384)
            ledger=work/'tonight/operator/projection/targeted_certificate/targeted_source_boxes.json';data=read(ledger)
            for box in data['boxes']:
                for rel in box['sha256']:box['sha256'][rel]=sha(work/rel)
            ledger.write_text(json.dumps(data,indent=2)+'\n')
            run(work,'targeted_outer_32',[work/'tonight/audit/global_b7171/replay_targeted_outer_boxes.py'],384)
            data=read(work/envelope)
            for item in data['accepted_R3_binding']:
                path=Path(item['path']);path=path if path.is_absolute()else work/path
                item['sha256']=sha(path)
            (work/envelope).write_text(json.dumps(data,indent=2)+'\n')
        run(work,'critical_native_672_endpoints',[work/'tonight/audit/global_b7171/replay_critical_native_endpoints.py'])
        run(work,'nine_R3_504_corners',[work/'tonight/audit/global_b7171/replay_ninebank_R3_corners.py'],384)
        run(work,'contact_all_line_coefficients',[work/'tonight/audit/verify_cubic_head_contact.py'])
        run(work,'linear_physical_source_exclusion',[work/'tonight/audit/verify_linear_contact_exclusion.py'])
        audit=work/'tonight/audit/global_b7171'
        run(work,'inverse_238_forward_bounds',[audit/'check_physical_inverse_bounds.py',work/envelope,
            '--output',audit/'critical_ninebank_inverse_bounds.json'])
        run(work,'selected_226_fixed_lambda_credits',[audit/'check_selected_contact_credits.py',work/envelope,
            '--output',audit/'critical_ninebank_contact_prices.json'])
        run(work,'exact_clipped_middle_cover',[work/'tonight/theory/audit_clipped_contact_envelope.py',envelope])
        run(work,'strengthened_global_assembly',[audit/'assemble_strengthened_acceptance.py'])
        result=read(audit/'STRENGTHENED_GLOBAL_B7171_ACCEPTANCE.json')
        assert result['status']=='PASS_INDEPENDENT_STRENGTHENED_COMPLETE_GLOBAL_B7171_SCALAR_BOUND'
        assert F(result['norm_upper_bound'])==F(56513883,10**8)
        assert F(result['constructive_G_R_lower_bound'])==F(10**8,56513883)
        for name,binding in result['bindings'].items():shutil.copy2(work/binding['path'],OUT/('strengthened_'+name+'.json'))
        summary={k:result[k]for k in('norm_upper_bound','constructive_G_R_lower_bound','middle_price_cells',
            'proved_physically_empty_middle_cells','independent_native_endpoint_prices','independent_R3_corner_prices',
            'independent_physical_inverse_caps','selected_contact_credit_cells')}
        summary.update(status='PASS_PORTABLE_STRENGTHENED_GLOBAL_B7171',
            original_strengthened_envelope_sha256=strong['original_envelope_sha256'],
            math_content_equal=True,immutable_evidence_manifest_pass=True,
            mode='complete_scalar_local_primitives_and_prices'if args.full else'accepted_fresh_primitive_reports_plus_fresh_prices_and_complete_cover',
            scalar_local_primitives_freshly_replayed_this_run=args.full,seconds=time.time()-start,
            scope='Accepted strengthened nine-profile/contact result only. Later joint-credit/middle-profile research upgrades are not global premises here.')
        (OUT/'result.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
