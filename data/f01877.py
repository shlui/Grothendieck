#!/usr/bin/env python3
"""Generate an ordered optimizer-free replay of the accepted geometry chain."""
from pathlib import Path
import argparse
import gzip
import json

ROOT=Path(__file__).resolve().parents[3]
AUDIT='tonight/audit/global_b7171/'
THEORY='tonight/theory/'

def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def path(value):
    p=Path(value);return p if p.is_absolute()else ROOT/p
def rel(p):return str(path(p).resolve().relative_to(ROOT))
def source_of(report):
    v=report['source'];return v['path']if isinstance(v,dict)else v

def main(x):
    accepted=read(x.acceptance);assert accepted['status']=='PASS_INDEPENDENT_COMPLETE_FULLUPPER_GLOBAL_SCALAR_CERTIFICATE'
    nodes=[];seen=set()
    def visit(report_path):
        rp=path(report_path);r=read(rp);assert r['status'].startswith('PASS_PRODUCER_FREE_')
        key=rel(rp)
        if key in seen:return
        if 'parent'in r:
            parent=r['parent'];visit(parent.get('accepted_geometry_audit',parent.get('accepted_audit')))
        seen.add(key);nodes.append((key,r))
    visit(accepted['bindings']['geometry']['path'])
    jobs=[]
    def job(name,args):jobs.append(dict(name=name,args=args,bits=768))
    for name,script in [('complete_cubic_L1_moment_prices','tonight/audit/refit847/replay_complete_refit_cover.py'),
        ('all_line_cubic_contact_coefficients','tonight/audit/verify_cubic_head_contact.py'),
        ('all_line_linear_source_coefficients','tonight/audit/verify_linear_contact_bank.py'),
        ('all_final57_native_endpoint_prices',AUDIT+'replay_final57_native_endpoints.py')]:job(name,[script])
    final_source=accepted['envelope']['path'];cache_paths={};plane_paths={};native_paths={};outputs={}
    for number,(report_path,r)in enumerate(nodes):
        source=rel(source_of(r));suffix=f'portable_chain_{number:02d}'
        if source==final_source:
            cache=accepted['bindings']['cache']['path'];plane=accepted['bindings']['plane']['path'];native=accepted['bindings']['native']['path']
        else:
            cache=AUDIT+suffix+'_R3_cache.json';plane=AUDIT+suffix+'_plane_prices.json';native=AUDIT+suffix+'_native_transport.json'
        cache_paths[source]=cache;plane_paths[source]=plane;native_paths[source]=native
        job(suffix+'_native_transport',[AUDIT+'check_final57_native_transport.py',source,'--output',native])
        job(suffix+'_all_R3_corners',[AUDIT+'replay_generic_R3_cache_57.py',source,'--output',cache])
        job(suffix+'_all_source_planes_and_inverse',[AUDIT+'check_fullupper57_plane_prices.py',source,'--output',plane])
        outputs[source]=report_path
        if r['status']=='PASS_PRODUCER_FREE_EXACT_POLYGON_GLOBAL_MIDDLE_COVER':
            job(suffix+'_complete_geometry',[THEORY+'audit_fullupper_polygon_envelope.py',source,
                '--cache-audit',cache,'--plane-audit',plane,'--native-audit',accepted['bindings']['native_baseline']['path'],
                '--output',report_path])
        elif r['status']=='PASS_PRODUCER_FREE_MONOTONE_R3_UPGRADE':
            p=r['parent'];parent_source=rel(p['path']);parent_report=p.get('accepted_geometry_audit',p.get('accepted_audit'))
            job(suffix+'_monotone_geometry',[THEORY+'audit_monotone_R3_upgrade_fresh.py',parent_source,source,
                '--parent-audit',parent_report,'--cache-audit',cache,'--output',report_path])
        elif r['status']=='PASS_PRODUCER_FREE_ADAPTIVE_POLYGON_REFINEMENT':
            p=r['parent'];job(suffix+'_adaptive_child_partition',[THEORY+'audit_adaptive_polygon_refinement.py',rel(p['path']),source,
                '--parent-audit',p['accepted_audit'],'--output',report_path])
        else:raise AssertionError(('UNSUPPORTED_GEOMETRY_TYPE',r['status']))
    contact_path=accepted['bindings']['contact']['path'];contact=read(path(contact_path))
    # Active final credits and inherited parent credits must both be freshly
    # checked. Inactive metadata have no mathematical role in the child bound.
    parent=contact.get('accepted_parent_contact_premise')
    if parent:
        pp=path(parent['path']);pr=read(pp);parent_source=None
        for _,g in nodes:
            if g['source_sha256']==pr['input_sha256']:parent_source=rel(source_of(g));break
        assert parent_source is not None,'Parent contact proof must have a node in the accepted exact geometry chain'
        job('parent_contact_credits',[AUDIT+'check_polygon_contact_credits.py',parent_source,
            '--r3-report',cache_paths[parent_source],'--output',rel(pp)])
        contact_checker=AUDIT+'check_adaptive_contact_credits.py'
    else:contact_checker=AUDIT+'check_polygon_contact_credits.py'
    job('final_active_and_inherited_contact_credits',[contact_checker,final_source,
        '--r3-report',accepted['bindings']['cache']['path'],'--output',contact_path])
    job('complete_chosen_b_lower_branch',[rel(x.lower_checker),'--b',accepted['operator_b'],
        '--output',accepted['bindings']['lower']['path']])
    assembly=[AUDIT+'assemble_fullupper_global_acceptance.py','--envelope',final_source]
    for key in('native_baseline','native','cache','plane','contact','geometry','lower'):
        assembly+=['--'+key.replace('_','-'),accepted['bindings'][key]['path']]
    assembly+=['--cap',accepted['norm_upper_bound'],'--output',rel(x.acceptance)]
    config=dict(status='FROZEN_COMPLETE_FULLUPPER_PARENT_CHAIN_REPLAY_CONFIG',
        accepted_cache_report=accepted['bindings']['cache']['path'],geometry_nodes=len(nodes),
        jobs=jobs,assembly_args=assembly,
        scope='Fresh native endpoints, every chain R3 cache/physical plane/inverse table, exact complete base and every monotone/adaptive child partition, active/inherited contact prices, same-b lower branch, and global assembly. Complete scalar/local leaf replay is optional and separately explicit.')
    x.output.write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps(dict(status=config['status'],geometry_nodes=len(nodes),jobs=len(jobs),output=str(x.output)),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--acceptance',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--lower-checker',type=Path,default=Path(AUDIT+'replay_lower_target_transport.py'))
    main(ap.parse_args())
