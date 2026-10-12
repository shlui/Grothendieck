"""Producer-free exact polygon coverage, knot splitting, and price audit.

No polygon clipping/maximization or envelope producer is imported. Vertex
feasibility and facet identities prove completeness of each recorded polygon.
Edge extrema are found using slope coordinates, independently of the producer's
parametric edge formula. Native chord intervals are split at every radial knot.
Numerical local gates/endpoints are independently replayed in bound ledgers.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse, bisect, gzip, hashlib, json, os, sys

os.environ.setdefault("GR_MPFR_BITS", "768")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tonight/operator/projection/certificate"))
from mpfr_interval import I, UP
from audit_general_r3_parameters_tail import gates
from audit_global_native_R3_envelope import PROFILES, RESERVES

PROFILES.update({
 "lowr76": "tonight/operator/projection/targeted_certificate/profile_interval_lowr_u76_r16_k0_2.json",
 "highk7575": "tonight/operator/projection/targeted_certificate/general_profile_highk_u7575_r259_k2_4.json"})
RESERVES.update({"lowr76": (F(1,1000), F(1,10**6)), "highk7575": (F(1,200), F(1,100000))})
PROFILES.update({
 "middle7475": "tonight/operator/projection/targeted_certificate/profile_interval_middle_u7475_r267_k08_1.json",
 "middlehigh7475": "tonight/operator/projection/targeted_certificate/profile_interval_middlehigh_u7475_r275_k11_2.json"})
RESERVES.update({"middle7475": (F(1,10000),F(1,10**6)),
                 "middlehigh7475": (F(1,10000),F(1,10**6))})
PROFILES.update({
 "interval_robust_middlehigh_box":"tonight/operator/projection/targeted_certificate/profile_interval_robust_middlehigh_box.json",
 "interval_robust_currentworst_box":"tonight/operator/projection/targeted_certificate/profile_interval_robust_currentworst_box.json",
 "robust_broadcurrent":"tonight/operator/projection/targeted_certificate/profile_interval_robust_broadcurrent_box.json"})
RESERVES.update({"interval_robust_middlehigh_box":(F(1,200000),F(1,10**8)),
 "interval_robust_currentworst_box":(F(1,50000),F(1,10**7)),
 "robust_broadcurrent":(F(1,50000),F(1,10**7))})

def read(p):
    with (gzip.open(p, "rt") if p.suffix == ".gz" else p.open()) as h:
        return json.load(h)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def frac_upper(x): return F(x.hi.text(UP, 240))
def iv(x):
    a,b=map(F,x);assert a<b;return a,b
def key(c): return tuple(iv(c[k]) for k in ("u_interval","r_interval","kappa_interval"))
def ceil12(x): return F(-((-x.numerator*10**12)//x.denominator),10**12)
def qmax(a,b,c,lo,hi):
    xs=[lo,hi]
    if a<0 and lo<-b/(2*a)<hi:xs.append(-b/(2*a))
    return max(a*x*x+b*x+c for x in xs)

def hull(points):
    """Canonical convex hull; point/segment cases are retained."""
    pts=sorted(set(points))
    if len(pts)<=1:return pts
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[]
    for p in pts:
        while len(lower)>1 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    upper=[]
    for p in reversed(pts):
        while len(upper)>1 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return lower[:-1]+upper[:-1]

def verify_polygon(points,ui,ri,planes):
    """Prove recorded hull equals rectangle/half-plane intersection.

    For full polygons all vertices are feasible and each hull facet is a
    boundary of an original constraint. Therefore both set inclusions hold.
    Since a>0, lower-left is feasible iff the intersection is nonempty;
    equality to a positive-slope plane forces a singleton there. Degenerate
    axis-aligned rectangles are checked separately by exact endpoints.
    """
    ulo,uhi=ui;rlo,rhi=ri
    assert ulo<=uhi and rlo<=rhi
    assert all(a>0 for a,c in planes)
    p=hull(points)
    if any(a*ulo+rlo>c for a,c in planes):assert not p;return p
    if any(a*ulo+rlo==c for a,c in planes):assert p==[(ulo,rlo)];return p
    if ulo==uhi:
        top=min([rhi]+[c-a*ulo for a,c in planes]);want=hull([(ulo,rlo),(ulo,top)])
        assert p==want;return p
    if rlo==rhi:
        right=min([uhi]+[(c-rlo)/a for a,c in planes]);want=hull([(ulo,rlo),(right,rlo)])
        assert p==want;return p
    assert len(p)>=3
    constraints=[(-F(1),F(0),-ulo),(F(1),F(0),uhi),
                 (F(0),-F(1),-rlo),(F(0),F(1),rhi)]+[(a,F(1),c)for a,c in planes]
    assert all(A*u+B*r<=C for u,r in p for A,B,C in constraints)
    for start,end in zip(p,p[1:]+p[:1]):
        assert any(A*start[0]+B*start[1]==C==A*end[0]+B*end[1]for A,B,C in constraints)
    return p

def strip_polygon(points,lo,hi):
    """Independent horizontal strip clipping by edge intersections."""
    out=[p for p in points if lo<=p[1]<=hi]
    for p,q in zip(points,points[1:]+points[:1]):
        if p[1]==q[1]:continue
        for r in (lo,hi):
            if min(p[1],q[1])<=r<=max(p[1],q[1]):
                u=p[0]+(q[0]-p[0])*(r-p[1])/(q[1]-p[1]);out.append((u,r))
    return hull(out)

def polygon_max(points,Lu,g,Lr,C=F(0)):
    """Exact boundary maximum using r=m*u+c or vertical-edge coordinates."""
    assert points
    def price(u,r):return u*u+Lu*u+g*r*r+Lr*r+C
    values=[price(*p)for p in points];interior=0
    for p,q in zip(points,points[1:]+points[:1]):
        if p==q:continue
        if p[0]==q[0]:
            if g<0:
                r=-Lr/(2*g)
                if min(p[1],q[1])<r<max(p[1],q[1]):values.append(price(p[0],r));interior+=1
        else:
            m=(q[1]-p[1])/(q[0]-p[0]);c=p[1]-m*p[0]
            A=F(1)+g*m*m;B=Lu+2*g*m*c+Lr*m
            if A<0:
                u=-B/(2*A)
                if min(p[0],q[0])<u<max(p[0],q[0]):values.append(price(u,m*u+c));interior+=1
    return max(values),interior

def native_polygon_price(p,ui,ri,ki,polygon,planes,lam=F(0),contact_constant=F(0),beta_term=F(0),a=F(0),d2=F(0),b=F(1)):
    lo,hi=ri;rg=p['rgrid'];e=p['endpoints'][ki]
    assert F(0)<=lo<=hi<=rg[-1]
    cuts=sorted(set([lo,hi]+[r for r in rg if lo<r<hi]))
    pieces=list(zip(cuts,cuts[1:])) if lo<hi else [(lo,hi)]
    shift=lam*(contact_constant+beta_term)+lam*lam*d2/(4*b)
    values=[];stationary=0
    for x,y in pieces:
        sub=strip_polygon(polygon,x,y)
        verify_polygon(sub,ui,(x,y),planes)
        if not sub:continue
        j=min(bisect.bisect_right(rg,x)-1,len(rg)-2)
        assert rg[j]<=x<=y<=rg[j+1]
        for corner in range(4):
            slope=(e[j+1][corner]-e[j][corner])/(rg[j+1]-rg[j])
            C=e[j][corner]-slope*rg[j]+shift
            value,n=polygon_max(sub,-2*p['ell']-lam*a,p['g'],slope-lam,C)
            values.append(value);stationary+=n
    assert values
    return ceil12(max(values)),len(pieces),stationary

def moment_price(q,k,s):
    sixth=I(15+15*k)+(I(194)/I(10).sqrt()-I(F(633,15)*s))*I(k)*I(k).sqrt()
    return I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*sixth

def native_price(p,ui,ri,ki,lam=F(0),contact_constant=F(0),beta_term=F(0),a=F(0),d2=F(0),b=F(1)):
    lo,hi=ri;rg=p["rgrid"];e=p["endpoints"][ki]
    assert F(0)<=lo<=hi<=rg[-1]
    cuts=sorted(set([lo,hi]+[r for r in rg if lo<r<hi]))
    pieces=list(zip(cuts,cuts[1:])) if lo<hi else [(lo,hi)]
    ell=p["ell"]+lam*a/2
    outer=max(u*u-2*ell*u for u in ui)
    shift=lam*(contact_constant+beta_term)+lam*lam*d2/(4*b)
    vals=[]
    for x,y in pieces:
        j=min(bisect.bisect_right(rg,x)-1,len(rg)-2)
        assert rg[j]<=x<=y<=rg[j+1]
        for corner in range(4):
            slope=(e[j+1][corner]-e[j][corner])/(rg[j+1]-rg[j])
            constant=e[j][corner]-slope*rg[j]
            vals.append(qmax(p["g"],slope-lam,constant+shift,x,y)+outer)
    return ceil12(max(vals))

def main(path,cache_path,plane_path,native_path,output):
    d=read(path);b=F(d["operator_b"]);assert F(0)<b<=F(7171,10000)
    assert d["status"]=="PASS_CONTINUOUS_ALL_SOURCE_NATIVE_R3_ENVELOPE_PRICING"
    assert not d.get('pending_gate_profiles',[])
    lo,hi=iv(d["amplitude_interval"])
    names=set(d['profile_coefficients']);assert names<=set(PROFILES)
    # Check independent cached-price status and exact identity of input file.
    cache_audit=read(cache_path)
    assert cache_audit["status"].startswith("PASS_INDEPENDENT_")
    assert cache_audit["source"]=={"path":str(path.relative_to(ROOT)),"sha256":sha(path)}
    plane_audit=read(plane_path)
    assert plane_audit['status']=='PASS_INDEPENDENT_UNIFORM_POLYGON_PLANES_AND_INVERSE_CAPS'
    assert plane_audit['input_sha256']==sha(path) and plane_audit['precision_bits']==768
    native_audit=read(native_path)
    assert native_audit['status'].startswith('PASS_INDEPENDENT_ALL_')
    native_baseline=F(native_audit['operator_b'])
    assert b<=native_baseline==F(7171,10000) and native_audit['precision_bits']==768
    native_frames={v['mask']:v for v in native_audit['frames']}
    coeff={};raws={};convexity={}
    binding={str(Path(v["path"]).resolve()):v["sha256"] for v in d["accepted_R3_binding"]}
    for p,h in binding.items():assert sha(Path(p))==h
    for name in sorted(names):
        rel=PROFILES[name]
        pth=ROOT/rel;assert sha(pth)==binding[str(pth)]
        raw=read(pth)["profile"]
        p={k:list(map(F,v)) if isinstance(v,list) else F(v) for k,v in raw.items()}
        assert b<=p["eps"];p["eps"]=b
        if "E00" not in p:p.update(E00=p["p0"]**2,E02=p["p0"]*p["p2"],E22=p["p2"]**2)
        gates(p);raws[name]=p
        yc=p["y"]-p["E02"]*p["r"]
        aa=b-p["E00"]+p["E02"]**2/p["E22"]
        bb=-2*p["gamma"]+2*yc*p["E02"]/p["E22"]
        coeff[name]=(aa,bb,p["ell"])
        stored=d["profile_coefficients"][name]
        assert (aa,bb,p["ell"])==tuple(F(stored[k]) for k in ("r_quadratic","r_linear","ell"))
        A=p["major"][3]+RESERVES[name][1]+p["zeta"]**2
        B=p["G"]-p["d0"]-p["zeta"]**2;assert A>=0 and B>=0
        convexity[name]="PASS_SEPARATELY_CONVEX_KAPPA_AND_SHAPE"
    assert set(raws)==names
    cache={}
    for t in d["cached_corner_constants"]:
        k=(t["profile"],iv(t["kappa_interval"]));assert k not in cache;cache[k]=F(t["upper"])
    replay_tables=cache_audit.get('new_tables',cache_audit.get('corners'))
    assert len(replay_tables)==len(cache)
    for t in replay_tables:
        k=(t["profile"],iv(t["kappa_interval"]))
        assert cache[k]==F(t["saved_upper"])
        assert len(t['corners'])==4
        # The independent report checks MPFR upper endpoints before formatting.
        # Its shorter decimal display may round above the saved rational by a
        # display ulp; it is not the numerical certificate's exact endpoint.
    frames={}
    for t in d["native_endpoint_inputs"]:
        mp=ROOT/t["mask"];pk=ROOT/t["package"]
        accepted=native_frames[t['mask']]
        assert all(accepted[k]==t[k] for k in ('mask','package','mask_sha256','source_sha256'))
        assert sha(mp)==t["mask_sha256"] and sha(pk/"data/source_profile.json")==t["source_sha256"]
        assert sha(pk/"replay.json")==t["replay_sha256"]
        m=read(mp);local=read(pk/"data/local_summary.json");p=local["profile"]
        assert local["status"]=="PASS" and local["source_sha256"]==t["source_sha256"]
        assert read(pk/"replay.json")["status"] in {"PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES","PASS_COMPLETE_REGIONAL_BOUND"}
        baseline=F(m['operator_b']);assert baseline==native_baseline
        assert F(p['epsilon'])>=baseline
        src=read(pk/"data/source_profile.json");rg=[F(str(x)) for x in src["r_grid"]]
        shift=F(p['beta'])**2*(1/b-1/baseline)
        endpoints={iv(x["kappa_interval"]):[[F(v)+shift for v in row] for row in x["prices"]] for x in m["endpoint_prices"]}
        M=[[F(v) for v in row] for row in p["M"]]
        g=(F(p["nu"])-M[0][2])**2/M[2][2]-M[0][0]+b-2*F(p["delta"])
        assert g==F(m["exact_gram_cost"])+b-baseline
        frames[t["mask"]]={"rgrid":rg,"endpoints":endpoints,"g":g,"ell":F(p["ell"]),"beta":F(p["beta"])}
    contact=read(ROOT/d["contact_gate_binding"]["path"])
    assert sha(ROOT/d["contact_gate_binding"]["path"])==d["contact_gate_binding"]["sha256"]
    a,headB=F(contact["parameters"]["a"]),F(contact["parameters"]["B"]);d2=6*headB**2
    assert headB>0
    bank=read(ROOT/d['linear_source_bank_binding']['path'])
    assert sha(ROOT/d['linear_source_bank_binding']['path'])==d['linear_source_bank_binding']['sha256']==plane_audit['bank_sha256']
    plane_rows={}
    def valid_row(row_index,kappa,cap,head):
        row=bank['rows'][row_index];assert F(row['a'])==head>0
        q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']));assert q[3]>=0
        k=(row_index,kappa,cap)
        if k not in plane_rows:
            assert all(frac_upper(moment_price(q,x,s))<=cap for x in kappa for s in (F(0),F(1)))
            plane_rows[k]=True
    planes={}
    for entry in d['source_polygon_planes']:
        kk=iv(entry['kappa_interval']);assert kk not in planes
        ps=[]
        for row in entry['planes']:
            A,C=F(row['a']),F(row['cap']);valid_row(row['source_row'],kk,C,A);ps.append((A,C))
        assert len(ps)==10 and len({A for A,C in ps})==10
        planes[kk]=ps
    ex=read(ROOT/d["linear_exclusion_binding"]["path"])
    assert sha(ROOT/d["linear_exclusion_binding"]["path"])==d["linear_exclusion_binding"]["sha256"]
    assert ex["status"]=="PASS_EXACT_GLOBAL_LINEAR_CONTACT_PHYSICAL_SOURCE_EXCLUSION"
    exa=F(ex["parameters"]["A"]);exc=max(F(v["price_upper"]) for v in ex["moment_price_corners"])
    assert list(map(F,ex["kappa_interval"]))==[F(2),F(4)] and list(map(F,ex["shape_interval"]))==[F(0),F(1)]
    bounds={F(v["r_lower"]):F(v["upper"]) for v in d["physical_u_bounds"]}
    rho=F(847,1000);c=(I(2)/I.pi()).sqrt()
    for r,u in bounds.items():
        if r==0:assert frac_upper(c)<=u
        else:
            v=(I(u)/c-I(1-rho*rho))/I(rho*rho)
            R=I(rho)*(-v.log()).sqrt().erf()
            assert frac_upper(R)<=r
    expected={};summaries={}
    for t in d["native_inputs"]:
        p=ROOT/t["path"];assert sha(p)==t["sha256"];src=read(p)
        assert F(src["operator_b"])==native_baseline and F(src["source_rho"])==rho
        for cell in src["cells"]:
            k=key(cell)
            if lo<=k[0][0]<k[0][1]<=hi:assert k not in expected;expected[k]=cell
        for s in src["summary"]:
            u=iv(s["u_interval"]);k=(u,iv(s["kappa_interval"]))
            if lo<=u[0]<u[1]<=hi:assert k not in summaries;summaries[k]=s
    ub=sorted({k[0] for k in expected});kb=sorted({k[2] for k in expected})
    assert ub[0][0]==lo and ub[-1][1]==hi and all(x[1]==y[0] for x,y in zip(ub,ub[1:]))
    assert kb[0][0]==0 and kb[-1][1]==90 and all(x[1]==y[0] for x,y in zip(kb,kb[1:]))
    assert set(summaries)=={(u,k) for u in ub for k in kb}
    assert set(cache)=={(name,k) for name in names for k in kb}
    assert set(planes)==set(kb)
    groups=defaultdict(list)
    for u,r,k in expected:groups[u,k].append(r)
    source_outward_checks={}
    for g,rs in groups.items():
        rs.sort();assert rs[0][0]==0 and all(x[1]==y[0] for x,y in zip(rs,rs[1:]));assert rs[-1][1]==F(summaries[g]["source_upper"])
        ulo=g[0][0];upper=rs[-1][1]
        if ulo not in source_outward_checks:
            v=(I(ulo)/c-I(1-rho*rho))/I(rho*rho)
            R=I(rho)*(-v.log()).sqrt().erf()
            assert frac_upper(R)<=upper
            source_outward_checks[ulo]=upper
        else:assert source_outward_checks[ulo]==upper
    seen=set();winner_counts=Counter();maximum=F(0);contact_checks={};knots=0;edge_stationary=0;polygon_count=0
    def clipped(k,record):
        u,r,kk=k;uh=min(u[1],bounds[r[0]]);rh=r[1]
        if F(2)<=kk[0]<=kk[1]<=F(4):uh=min(uh,(exc-r[0])/exa);rh=min(rh,exc-exa*u[0])
        sc=record['source_planes']
        pA,pC=F(sc['u_a']),F(sc['u_price_upper']);qA,qC=F(sc['r_a']),F(sc['r_price_upper'])
        valid_row(sc['u_row'],kk,pC,pA);valid_row(sc['r_row'],kk,qC,qA)
        uh=min(uh,(pC-r[0])/pA);rh=min(rh,qC-qA*u[0])
        return (u[0],uh),(r[0],rh)
    for empty in d["physically_empty_cells"]:
        k=key(empty);assert k in expected and k not in seen;seen.add(k)
        u,r=clipped(k,empty)
        if 'physical_u_upper' in empty:assert F(empty['physical_u_upper'])==u[1]
        if 'physical_r_upper' in empty:assert F(empty['physical_r_upper'])==r[1]
        if empty.get('polygon_empty',False):
            assert u[0]<=u[1] and r[0]<=r[1]
            verify_polygon([],u,r,planes[k[2]])
        else:assert u[1]<u[0] or r[1]<r[0]
    for cell in d["cells"]:
        k=key(cell);assert k in expected and k not in seen;seen.add(k)
        u,r=clipped(k,cell);assert u[0]<=u[1] and r[0]<=r[1]
        assert F(cell["physical_u_upper"])==u[1] and F(cell["physical_r_upper"])==r[1]
        kk=k[2];kind,name=cell["winner"].split(":",1);is_contact=kind.endswith("-contact")
        use_polygon=cell.get('polygon_pricing','physical_polygon' in cell)
        if use_polygon:
            polygon=verify_polygon([tuple(map(F,p))for p in cell['physical_polygon']],u,r,planes[kk]);assert polygon
            polygon_count+=1
        lam=Aq=beta=F(0)
        if is_contact:
            cc=cell["contact_credit"];lam=F(cc["lambda_value"]);assert lam>=0
            index=cc["contact_row"];row=contact["rows"][index];assert row["w"]==cc["contact_w"]
            q=list(map(F,row["polynomial_coefficients_q0_q2_q4_q6"]));assert q[3]>=0
            Aq=F(cc["gate_price_upper"]);ck=(index,kk,Aq)
            if ck not in contact_checks:
                assert all(frac_upper(moment_price(q,x,s))<=Aq for x in kk for s in (F(0),F(1)));contact_checks[ck]=True
            beta=F(cc["beta_d_over_b_upper"])
            if kind.startswith("native"):
                truebeta=I(frames[name]["beta"])*I(6).sqrt()*I(headB)/I(b)
                assert frac_upper(truebeta)<=beta
            else:
                assert 6*raws[name]["beta"]*headB/b<=beta
            assert F(cc["S_upper"])==Aq-a*u[0]-r[0]+beta
        if kind.startswith("native"):
            if use_polygon:
                value,n,s=native_polygon_price(frames[name],u,r,kk,polygon,planes[kk],lam,Aq,beta,a,d2,b)
                knots+=n;edge_stationary+=s
            else:value=native_price(frames[name],u,r,kk,lam,Aq,beta,a,d2,b)
        else:
            aa,bb,ell=coeff[name];ell+=lam*a/2
            shift=lam*(Aq+beta)+lam*lam*d2/(4*b)
            if use_polygon:
                v,s=polygon_max(polygon,-2*ell,aa,bb-lam,cache[name,kk]+shift);value=ceil12(v);edge_stationary+=s
            else:value=ceil12(cache[name,kk]+shift+max(x*x-2*ell*x for x in u)+qmax(aa,bb-lam,F(0),*r))
        assert value==F(cell["cap"]),(kind,name,k,value,cell["cap"])
        assert cell["target_pass"]==(value<F(40,71)) and cell["fallback_pass"]==(value<F(451,800))
        maximum=max(maximum,value);winner_counts[kind+":"+name]+=1
    assert seen==set(expected)
    assert maximum==F(d["maximum_cell_cap"])<F(d["strict_rational_cap"])
    out={"status":"PASS_PRODUCER_FREE_EXACT_POLYGON_GLOBAL_MIDDLE_COVER",
         "scope":"Complete physical middle interval only. A global norm additionally requires accepted outer branches at this same coefficient.",
         "source":str(path.relative_to(ROOT)),"source_sha256":sha(path),"operator_b":str(b),
         "amplitude_interval":d["amplitude_interval"],"u_bands":len(ub),"kappa_bands":len(kb),
         "complete_groups":len(groups),"original_cells":len(expected),"accepted_price_cells":len(d["cells"]),
         "exact_physically_empty_cells":len(d["physically_empty_cells"]),"directed_inverse_caps":len(bounds),
         "recorded_complete_polygons":polygon_count,"native_knot_pieces":knots,"checked_edge_stationary_points":edge_stationary,
         "uniform_source_planes":sum(map(len,planes.values())),"independently_checked_plane_row_band_groups":len(plane_rows),
         "independent_outward_source_grid_endpoint_checks":len(source_outward_checks),
         "contact_corner_groups":len(contact_checks),"winners":dict(winner_counts),
         "maximum_cell_cap":str(maximum),"strict_rational_cap":d["strict_rational_cap"],
         "cache_replay_status":cache_audit["status"],"convexity":convexity}
    out['R3_profiles']=len(names)
    out['native_endpoint_baseline_b']=str(native_baseline)
    out['native_transport_formula']='Each normalized endpoint gains beta^2*(1/b-1/b0); its source quadratic gains b-b0.'
    out['cached_cumulant_tables']=len(cache)
    out['whole_line_binding_files']=len(binding)
    out['numeric_dependency_replays']={str(p.relative_to(ROOT)):sha(p) for p in (cache_path,plane_path,native_path)}
    output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("path")
    parser.add_argument('--cache-audit',required=True);parser.add_argument('--plane-audit',required=True)
    parser.add_argument('--native-audit',required=True)
    parser.add_argument('--output',default='tonight/theory/preferred_polygon_envelope_audit.json')
    args=parser.parse_args();main(ROOT/args.path,ROOT/args.cache_audit,ROOT/args.plane_audit,ROOT/args.native_audit,ROOT/args.output)
