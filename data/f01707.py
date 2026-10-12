#!/usr/bin/env python3
"""Fresh exact half-line certification with explicit per-row reserves.

Runs the audited common verifier without discovery code; failed reserves do
not constitute acceptance. Existing correctly bound row certificates are kept.
"""
import argparse,hashlib,json,time
from pathlib import Path
from fractions import Fraction as F
import certify_continuous_profile_scalars as S

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--q0',type=F,default=F(1,20000))
    ap.add_argument('--q6',type=F,default=F(1,100000000));ap.add_argument('--bits',type=int,default=512)
    args=ap.parse_args();assert args.q0>=0 and args.q6>0
    S.ctx.prec=args.bits;S.ctx.cap=5;args.output_dir.mkdir(parents=True,exist_ok=True)
    data=json.loads(args.source.read_text());sha=hashlib.sha256(args.source.read_bytes()).hexdigest()
    attempt_path=args.output_dir/'scalar_repair_attempts.json'
    attempts=json.loads(attempt_path.read_text())if attempt_path.exists()else[]
    reserves=sorted(set([args.q6,F(3,100000000),F(1,10000000),F(3,10000000),F(1,1000000),F(2,1000000),F(3,1000000),F(1,100000)]))
    for index in range(len(data['rows'])):
        out=args.output_dir/f'continuous_scalar_row_{index:02d}.json'
        if out.exists():
            prior=json.loads(out.read_text());assert prior['source_sha256']==sha and prior['row_index']==index and prior['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE'
            continue
        for q6 in reserves:
            S.Q0_SLACK=args.q0;S.Q6_SLACK=q6;start=time.time()
            try:
                row=S.freeze(data,index);cutoff,tail=S.tail(row);compact=S.compact(row,cutoff)
            except AssertionError as error:
                attempts.append({'row':index,'q0':str(args.q0),'q6':str(q6),'status':'FAIL_NOT_ACCEPTED','error':str(error)})
                print(json.dumps(attempts[-1]),flush=True)
                continue
            frozen=dict(r=str(row['r']),lambda_L1=str(row['mu']),eta=str(row['eta']),kappa_interval=list(map(str,row['kappa_interval'])),majorant=list(map(str,row['q'])),phi_coefficients=list(map(str,row['coefficients'])),degree=10,scale='2',regression_hinges=list(map(str,S.HINGES)),added_q0=str(args.q0),added_q6=str(q6))
            result=dict(status='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE',row_index=index,precision_bits=args.bits,frozen_row=frozen,compact=compact,tail=tail,prices=S.prices(row),seconds=time.time()-start,source_sha256=sha,script_sha256=hashlib.sha256(Path(S.__file__).read_bytes()).hexdigest(),scope='Whole-line scalar inequality for a frozen exact profile. Native physical gates and source coverage are separate.')
            out.write_text(json.dumps(result,indent=2)+'\n')
            attempts.append({'row':index,'q0':str(args.q0),'q6':str(q6),'status':'PASS','leaves':compact['leaf_count'],'tail_cutoff':str(cutoff)})
            print(json.dumps(attempts[-1]),flush=True);break
        else:raise AssertionError(('NO_ACCEPTED_RESERVE',index))
        (args.output_dir/'scalar_repair_attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
    print(json.dumps({'status':'COMPLETE_FRESH_ROW_CERTIFICATES','rows':len(data['rows']),'source_sha256':sha}),flush=True)

if __name__=='__main__':main()
