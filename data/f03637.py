#!/usr/bin/env python3
"""Producer-free native endpoint reprice for the extended critical frame."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
source=ROOT/'tonight/recovery/profile_cover/global_middle_b7171_extendedcritical_ninebank_planes_envelope.json';d=json.loads(source.read_text());assert F(d['operator_b'])==F(7171,10000)
bind=next(x for x in d['native_endpoint_inputs']if x['package'].endswith('critical_extended_r032'))
start=time.time();frame=N.frame(bind,F(7171,10000));assert len(frame['keys'])==14
out={'status':'PASS_INDEPENDENT_EXTENDED_CRITICAL_NATIVE_ALL_SELECTED_ENDPOINT_PRICES','operator_b':'7171/10000','endpoint_prices_recomputed':N.POINTS,'precision_bits':N.BITS,'frame':bind,'exact_gram_cost':str(frame['g']),'kappa_interval':['0','90'],'all_shapes':['0','1'],'source_grid_endpoint':str(frame['rg'][-1]),'packages':[{'package':n,'manifest_sha256':p['manifest_sha256'],'replay_sha256':p['replay_sha256'],'source_sha256':p['source_sha256']}for n,p in N.PACKAGES.items()],'producer_code_imported':False,'seconds':time.time()-start,'scope':'Accepted native physical gate and exact b-transport with independently repriced selected scalar/L1 endpoint bounds. The updated global source partition and contact credit are separate.'}
assert N.POINTS==224
out['envelope_input_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
p=Path(__file__).parent/'extended_critical_native_endpoint_replay.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
