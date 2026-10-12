"""Manual exact adversarial fixtures for the NEW affine replay helper.

The fixture polygons and expected regions are written explicitly, without
using the helper's clipping or coefficient routines to generate expectations.
"""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).parent/'portable_geometry'))
from audit_affine_hybrid_geometry import replay

def row(name,au=0,ar=0,g=0,C=0):
 return dict(name=name,u_linear=str(au),r_linear=str(ar),r_quadratic=str(g),constant=str(C))
def cand(under,rect,au,ar,C,kind='SECANT'):
 return dict(name=under['name'],underlying=under,rectangle=list(map(str,rect)),u_linear=str(au),r_linear=str(ar),constant=str(C),r_majorant=kind)
def region(name,poly,values,upper):
 return dict(profile=name,vertices=[[str(x),str(y)]for x,y in poly],vertex_values=list(map(str,values)),upper=str(upper))
def piece(poly,ra,rb,candidates,regions,upper,point,bindings=None):
 return dict(polygon=[[str(x),str(y)]for x,y in poly],r_interval=[str(ra),str(rb)],upper=str(upper),native_bindings=bindings or [],dominance_proof=dict(candidates=candidates,regions=regions,attaining_point=list(map(str,point))))
def raw(poly):return [[str(x),str(y)]for x,y in poly]
cases=[]
square=[(0,0),(1,0),(1,1),(0,1)]
a=row('a')
p=piece(square,0,1,[cand(a,[0,1,0,1],1,0,0)],[region('a',square,[0,1,1,0],1)],1,(1,0))
cases.append(('square_secant',square,dict(pieces=[p]),[a],[],F(1)))
seg=[(0,0),(1,0)];a=row('a',-1);b=row('b',-1)
p=piece(seg,0,0,[cand(a,[0,1,0,0],0,0,0),cand(b,[0,1,0,0],0,0,0)],[region('a',seg,[0,0],0),region('b',seg,[0,0],0)],0,(0,0))
cases.append(('identical_closed_segment',seg,dict(pieces=[p]),[a,b],[],F(0)))
point=[(F(1,2),F(1,2))];a=row('a',g=-2)
p=piece(point,F(1,2),F(1,2),[cand(a,[F(1,2)]*4,1,-2,F(1,4),'TANGENT')],[region('a',point,[F(-1,4)],F(-1,4))],F(-1,4),point[0])
cases.append(('point_negative_tangent',point,dict(pieces=[p]),[a],[],F(-1,4)))
tri=[(0,0),(1,0),(0,1)];a=row('a');b=row('b',-2,C=1)
left=[(0,0),(F(1,2),0),(F(1,2),F(1,2)),(0,1)]
right=[(F(1,2),0),(1,0),(F(1,2),F(1,2))]
p=piece(tri,0,1,[cand(a,[0,1,0,1],1,0,0),cand(b,[0,1,0,1],-1,0,1)],[region('a',left,[0,F(1,2),F(1,2),0],F(1,2)),region('b',right,[F(1,2),0,F(1,2)],F(1,2))],F(1,2),(F(1,2),0))
cases.append(('interior_dominance_tie',tri,dict(pieces=[p]),[a,b],[],F(1,2)))
native=dict(name='fixture',table_indices=[0],grid=[0,F(1,2),1],caps=[0,F(1,2),0],ell=0,g=-1)
nname='NATIVE:fixture';n1=row(nname,ar=1,g=-1);n2=row(nname,ar=-1,g=-1,C=1)
poly1=[(0,F(1,2)),(0,0),(1,0),(1,F(1,2))]
poly2=[(0,F(1,2)),(1,F(1,2)),(1,1),(0,1)]
bind1=dict(name=nname,table_indices=[0],piece_index=0,source_interval=['0','1/2'],endpoint_caps=['0','1/2'])
bind2=dict(name=nname,table_indices=[0],piece_index=1,source_interval=['1/2','1'],endpoint_caps=['1/2','0'])
p1=piece(poly1,0,F(1,2),[cand(n1,[0,1,0,F(1,2)],1,F(1,2),F(1,16),'TANGENT')],[region(nname,poly1,[F(5,16),F(1,16),F(17,16),F(21,16)],F(21,16))],F(21,16),(1,F(1,2)),[bind1])
p2=piece(poly2,F(1,2),1,[cand(n2,[0,1,F(1,2),1],1,F(-5,2),F(25,16),'TANGENT')],[region(nname,poly2,[F(5,16),F(21,16),F(1,16),F(-15,16)],F(21,16))],F(21,16),(1,F(1,2)),[bind2])
cases.append(('native_knot_closed_cover',square,dict(pieces=[p1,p2]),[],[native],F(21,16)))
accepted=[]
for name,poly,proof,rows,natives,expected in cases:
 got,stats=replay(raw(poly),proof,rows,natives)
 assert got==expected
 accepted.append(dict(name=name,exact_upper=str(got),stats=stats))
def reject(name,poly,proof,rows,natives):
 try:replay(raw(poly),proof,rows,natives)
 except AssertionError:return name
 raise AssertionError('tamper accepted: '+name)
rejections=[]
q=deepcopy(cases[-1][2]);q['pieces'].pop();rejections.append(reject('missing_radial_half',square,q,[],[native]))
q=deepcopy(cases[-1][2]);q['pieces'][0]['native_bindings'][0]['endpoint_caps'][1]='0';rejections.append(reject('wrong_native_endpoint_cap',square,q,[],[native]))
q=deepcopy(cases[3][2]);q['pieces'][0]['dominance_proof']['regions'].pop();rejections.append(reject('missing_closed_dominance_region',tri,q,cases[3][3],[]))
q=deepcopy(cases[2][2]);q['pieces'][0]['dominance_proof']['candidates'][0]['constant']='0';rejections.append(reject('wrong_negative_tangent_constant',point,q,cases[2][3],[]))
q=deepcopy(cases[0][2]);q['pieces'][0]['dominance_proof']['candidates'][0]['rectangle'][1]='1/2';rejections.append(reject('majorant_rectangle_omits_vertices',square,q,cases[0][3],[]))
out=dict(status='PASS_MANUAL_EXACT_AFFINE_GEOMETRY_FIXTURES',accepted=accepted,rejected_tamper_cases=rejections,scope='Independent manually specified fixtures; supplements complete frozen-data replay, does not replace it.')
print(json.dumps(out,indent=2))
(Path(__file__).parent/'AFFINE_HELPER_MANUAL_FIXTURES.json').write_text(json.dumps(out,indent=2)+'\n')
