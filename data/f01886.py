from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
p=HERE/'NEXT_FINAL68_REPLAY_CONFIG.json';d=json.loads(p.read_text())
a=json.loads((HERE/'FULLUPPER68_NATIVE_REFINED_GLOBAL_ACCEPTANCE.json').read_text())
A='tonight/audit/global_b7171/';T='tonight/theory/'
r=a['refinement_bindings'];source=a['envelope']['path']
d['jobs'].append(dict(name='complete_parent27_global_assembly',args=d['assembly_args'],bits=768))
jobs=[
    ('all_final68_native_endpoint_prices',[A+'replay_final68_native_endpoints.py']),
    ('final68_target_b_native_transport',[A+'check_final68_native_transport.py',source,'--output',r['native']['path']]),
    ('final68_all_extra_physical_planes',[A+'check_native68_overlay_planes.py',source,'--output',r['plane']['path']]),
    ('exact_final68_child_partition_geometry',[T+'audit_native_kappa_refinement.py',
        'tonight/recovery/profile_cover/global_upper_b713_finekappa57_twentysevenbank_refined_envelope.json.gz',
        source,'--parent-audit',a['bindings']['geometry']['path'],'--native-audit',r['baseline']['path'],
        '--output',r['geometry']['path']]),
]
for name,args in jobs:d['jobs'].append(dict(name=name,args=args,bits=768))
args=[A+'assemble_native_kappa_refinement_acceptance.py','--envelope',source]
for key in('parent_global','geometry','native','plane','baseline'):args+=['--'+key.replace('_','-'),r[key]['path']]
args+=['--output',A+'FULLUPPER68_NATIVE_REFINED_GLOBAL_ACCEPTANCE.json']
d['assembly_args']=args;d['geometry_nodes']+=1
d['scope']+=' The final native-only cumulant overlay also freshly checks all13936 new native endpoints, exact target-b transport, every new source plane, all child intervals and exact polygon/knot/edge prices.'
p.write_text(json.dumps(d,indent=2)+'\n')
print('Final68 config jobs',len(d['jobs']))
