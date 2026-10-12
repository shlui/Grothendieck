#!/usr/bin/env python3
"""Promote only a fully replayed, hash-bound preferred global package."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,zipfile,subprocess,sys,tempfile

ROOT=Path(__file__).resolve().parents[3]
EXCLUDED={'replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def main(stage):
    stage=stage.resolve();result=read(stage/'replay_outputs/result.json')
    assert result['status']=='PASS_PORTABLE_COMPLETE_FULLUPPER_PREFERRED_GLOBAL'
    provenance=read(stage/'PROVENANCE.json');pref=provenance['preferred']
    assert result['operator_b']==pref['operator_b']=='711/1000'
    assert result['norm_upper_bound']==pref['norm_cap']=='112782853/200000000'
    assert result['constructive_G_R_lower_bound']==pref['GR_lower_bound']=='200000000/112782853'
    assert result['final_uniform_target_coefficient_transport_freshly_replayed']is True
    assert result['complete_whole_q_profiles']==29 and result['cumulant_points']==68
    assert result['full_upper_priced_cells']==272579 and result['proved_physically_empty_full_upper_cells']==8737
    old=read(stage/'MANIFEST.json')
    assert all(sha(stage/rel)==h for rel,h in old['sha256'].items())
    immutable=[]
    for rel in('tonight/operator/native_global_b7171','tonight/operator/native_global_b7171_strengthened'):
        p=ROOT/rel;m=read(p/'MANIFEST.json')
        assert all(sha(p/name)==h for name,h in m['sha256'].items())
        immutable.append({'path':rel,'manifest_sha256':sha(p/'MANIFEST.json'),'unchanged':True})
    backupzip=ROOT/'tonight/operator/preferred_global_final14_b713.zip'
    assert sha(backupzip)=='ca2db028c1769337b65d954de6ae306717d5f21c192788e83eee5658b2a64409'
    shutil.copy2(stage/'replay_outputs/result.json',stage/'DEFAULT_REPLAY_RESULT.json')
    shutil.copy2(stage/'replay_outputs/preferred_complete_acceptance.json',stage/'DEFAULT_REPLAY_ACCEPTANCE.json')
    validation={'status':'PASS_FINAL711_COMPLETE_PORTABLE_GLOBAL_REPLAY',
        'default_result':result,'default_replay_result_sha256':sha(stage/'DEFAULT_REPLAY_RESULT.json'),
        'default_replay_acceptance_sha256':sha(stage/'DEFAULT_REPLAY_ACCEPTANCE.json'),
        'full_prerequisite_leaf_replay_repeated_this_packaging_run':False,
        'prerequisite_evidence':'Complete previously fresh accepted scalar/local leaf replay reports are bound to their exact profiles, reserves and trees. The default freshly replays numeric prices, complete exact partition geometry, same-b lower coverage and global assembly. The optional --full mode and every input tree are included.',
        'final_two_profile_optional_exact_tail_smoke':'Both final profile exact tails were freshly replayed from the portable payload in this packaging run.',
        'post_default_changes':'Only immutable copies of the PASS result/acceptance and this validation record were added; mathematical inputs and default verifier code are unchanged.',
        'phase_validation':'The original/strengthened base freshly passed, then its unchanged byte-bound checkpoint was reused for all40 fresh713 ancestor/global jobs. The four final711 transport/assembly jobs subsequently passed against the identical frozen713 proof inputs. The default entrypoint includes all phases in one fresh ordered run; no mathematical stage is omitted.',
        'commands':{'default':'python preferred_global/verify_bundle.py','full':'python preferred_global/verify_bundle.py --full --workers 4'},
        'original_and_strengthened_packages':immutable,
        'previous_final14_checkpoint_sha256':sha(backupzip),
        'fresh_base_checkpoint_seconds':read(stage/'BASE_REPLAY_RESULT.json')['seconds'],
        'successful_validation_phase_seconds':read(stage/'BASE_REPLAY_RESULT.json')['seconds']+result['seconds'],
        'single_fresh_default_invocation_repeated_after_final711_append':False,
        'extension_filter_used':False,'original_envelope_and_metadata_archives_preserved':True,
        'scope':'Explicit constructive lower bound; no three-digit enclosure of unrestricted G_R is asserted.'}
    write(stage/'VALIDATION.json',validation)
    files=sorted(p for p in stage.rglob('*')if p.is_file()and not any(k in EXCLUDED for k in p.relative_to(stage).parts)and p.name!='MANIFEST.json')
    manifest={'status':'FROZEN_FINAL711_COMPLETE_PREFERRED_GLOBAL_IMMUTABLE_INPUT_MANIFEST',
        'extension_filter_used':False,'excluded_mutable_directories':sorted(EXCLUDED),
        'sha256':{str(p.relative_to(stage)):sha(p)for p in files}}
    write(stage/'MANIFEST.json',manifest)
    assert all(sha(stage/rel)==h for rel,h in manifest['sha256'].items())
    canonical=ROOT/'tonight/operator/preferred_global';zip_path=ROOT/'tonight/operator/preferred_global.zip'
    temporary_zip=zip_path.with_suffix('.zip.next_final711');members=[stage/'MANIFEST.json',*files]
    with zipfile.ZipFile(temporary_zip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6)as z:
        for p in sorted(members):
            info=zipfile.ZipInfo('preferred_global/'+str(p.relative_to(stage)),date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,p.read_bytes(),compresslevel=6)
    with zipfile.ZipFile(temporary_zip)as z:
        assert z.testzip()is None
        assert set(z.namelist())=={'preferred_global/'+str(p.relative_to(stage))for p in members}
        for rel,h in manifest['sha256'].items():assert hashlib.sha256(z.read('preferred_global/'+rel)).hexdigest()==h
        assert hashlib.sha256(z.read('preferred_global/MANIFEST.json')).hexdigest()==sha(stage/'MANIFEST.json')
    with tempfile.TemporaryDirectory(prefix='preferred711_extraction_smoke_')as directory:
        extracted=Path(directory)/'preferred_global'
        with zipfile.ZipFile(temporary_zip)as z:z.extractall(directory)
        smoke=subprocess.run([sys.executable,str(extracted/'verify_bundle.py'),'--manifest-only'],cwd=directory,capture_output=True,text=True)
        assert smoke.returncode==0,smoke.stdout+smoke.stderr
        env=dict(os.environ);env['GR_MPFR_BITS']='768';env['PYTHONINTMAXSTRDIGITS']='0'
        env['PYTHONPATH']=os.pathsep.join(str(extracted/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
        beta=subprocess.run([sys.executable,str(extracted/'tonight/audit/global_b7171/replay_assigned_beta_bindings.py'),
            '--scan',str(extracted/'tonight/recovery/profile_cover/assigned_twenty-nine_coefficient_scan.json'),
            '--parent-global',str(extracted/'tonight/audit/global_b7171/FINAL29_GLOBAL_ACCEPTANCE.json'),
            '--b','711/1000','--output',str(Path(directory)/'fresh_beta_spectral_smoke.json')],
            cwd=directory,env=env,capture_output=True,text=True)
        assert beta.returncode==0,beta.stdout+beta.stderr
        smoke_result=json.loads(smoke.stdout);beta_result=read(Path(directory)/'fresh_beta_spectral_smoke.json')
        assert beta_result['selected_profiles']==43
    backup=ROOT/'tonight/operator/preferred_global_final14_b713_tree';assert not backup.exists()
    if canonical.exists():canonical.rename(backup)
    stage.rename(canonical);os.replace(temporary_zip,zip_path)
    summary={'status':'PASS_ATOMIC_FINAL711_COMPLETE_PREFERRED_GLOBAL_PACKAGE',
        'path':str(canonical),'zip':str(zip_path),'zip_sha256':sha(zip_path),'zip_bytes':zip_path.stat().st_size,
        'immutable_manifest_files':len(files),'manifest_sha256':sha(canonical/'MANIFEST.json'),
        'preferred_and711_phase_seconds':result['seconds'],
        'base_checkpoint_seconds':read(canonical/'BASE_REPLAY_RESULT.json')['seconds'],
        'single_final711_default_invocation_repeated':False,
        'operator_b':result['operator_b'],
        'norm_cap':result['norm_upper_bound'],'constructive_G_R_lower_bound':result['constructive_G_R_lower_bound'],
        'crc_and_all_zip_member_hashes_pass':True,'original_and_strengthened_unchanged':True,
        'path_independent_actual_zip_extraction_smoke':smoke_result,
        'extracted_fresh43_beta_profile_target_spectral_bindings':beta_result['status'],
        'previous_final14_checkpoint_preserved':True,'extension_filter_used':False}
    write(ROOT/'tonight/audit/global_b7171/PREFERRED_PACKAGE_CHECKPOINT.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--staging',required=True,type=Path);main(p.parse_args().staging)
