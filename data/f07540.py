"""Producer-free exact physical cover and selected contact-price audit.

The native endpoints and nine R3 corner caches have independent replays.
This reconstructs every chosen bound, proves every clipped or discarded
cell retains all physical points, and checks exhaustive original partitions.
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

def main(path):
    d=read(path);b=F(d["operator_b"]);assert b==F(7171,10000)
    assert d["status"]=="PASS_CONTINUOUS_ALL_SOURCE_NATIVE_R3_ENVELOPE_PRICING"
    lo,hi=iv(d["amplitude_interval"])
    # Check independent cached-price status and exact identity of input file.
    cache_audit=read(ROOT/"tonight/audit/global_b7171/ninebank_R3_corner_replay.json")
    assert cache_audit["status"]=="PASS_INDEPENDENT_NINEBANK_504_R3_CORNER_PRICES"
    assert cache_audit["source"]=={"path":str(path.relative_to(ROOT)),"sha256":sha(path)}
    coeff={};raws={};convexity={}
    binding={str(Path(v["path"]).resolve()):v["sha256"] for v in d["accepted_R3_binding"]}
    for p,h in binding.items():assert sha(Path(p))==h
    for name,rel in PROFILES.items():
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
    assert set(d["profile_coefficients"])==set(PROFILES)
    cache={}
    for t in d["cached_corner_constants"]:
        k=(t["profile"],iv(t["kappa_interval"]));assert k not in cache;cache[k]=F(t["upper"])
    assert len(cache_audit["corners"])==len(cache)
    for t in cache_audit["corners"]:
        k=(t["profile"],iv(t["kappa_interval"]))
        assert cache[k]==F(t["saved_upper"]) and F(t["fresh_upper"])<=cache[k]
    frames={}
    for t in d["native_endpoint_inputs"]:
        mp=ROOT/t["mask"];pk=ROOT/t["package"]
        assert sha(mp)==t["mask_sha256"] and sha(pk/"data/source_profile.json")==t["source_sha256"]
        assert sha(pk/"replay.json")==t["replay_sha256"]
        m=read(mp);local=read(pk/"data/local_summary.json");p=local["profile"]
        assert local["status"]=="PASS" and local["source_sha256"]==t["source_sha256"]
        assert read(pk/"replay.json")["status"] in {"PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES","PASS_COMPLETE_REGIONAL_BOUND"}
        assert F(m["operator_b"])==b
        src=read(pk/"data/source_profile.json");rg=[F(str(x)) for x in src["r_grid"]]
        endpoints={iv(x["kappa_interval"]):[list(map(F,row)) for row in x["prices"]] for x in m["endpoint_prices"]}
        M=[[F(v) for v in row] for row in p["M"]]
        g=(F(p["nu"])-M[0][2])**2/M[2][2]-M[0][0]+b-2*F(p["delta"])
        assert g==F(m["exact_gram_cost"])
        frames[t["mask"]]={"rgrid":rg,"endpoints":endpoints,"g":g,"ell":F(p["ell"]),"beta":F(p["beta"])}
    contact=read(ROOT/d["contact_gate_binding"]["path"])
    assert sha(ROOT/d["contact_gate_binding"]["path"])==d["contact_gate_binding"]["sha256"]
    a,headB=F(contact["parameters"]["a"]),F(contact["parameters"]["B"]);d2=6*headB**2
    assert headB>0
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
        assert F(src["operator_b"])==b and F(src["source_rho"])==rho
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
    assert set(cache)=={(name,k) for name in PROFILES for k in kb}
    groups=defaultdict(list)
    for u,r,k in expected:groups[u,k].append(r)
    for g,rs in groups.items():
        rs.sort();assert rs[0][0]==0 and all(x[1]==y[0] for x,y in zip(rs,rs[1:]));assert rs[-1][1]==F(summaries[g]["source_upper"])
    seen=set();winner_counts=Counter();maximum=F(0);contact_checks={}
    def clipped(k):
        u,r,kk=k;uh=min(u[1],bounds[r[0]]);rh=r[1]
        if F(2)<=kk[0]<=kk[1]<=F(4):uh=min(uh,(exc-r[0])/exa);rh=min(rh,exc-exa*u[0])
        return (u[0],uh),(r[0],rh)
    for empty in d["physically_empty_cells"]:
        k=key(empty);assert k in expected and k not in seen;seen.add(k)
        u,r=clipped(k);assert u[1]<u[0] or r[1]<r[0]
    for cell in d["cells"]:
        k=key(cell);assert k in expected and k not in seen;seen.add(k)
        u,r=clipped(k);assert u[0]<=u[1] and r[0]<=r[1]
        assert F(cell["physical_u_upper"])==u[1] and F(cell["physical_r_upper"])==r[1]
        kk=k[2];kind,name=cell["winner"].split(":",1);is_contact=kind.endswith("-contact")
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
            value=native_price(frames[name],u,r,kk,lam,Aq,beta,a,d2,b)
        else:
            aa,bb,ell=coeff[name];ell+=lam*a/2
            shift=lam*(Aq+beta)+lam*lam*d2/(4*b)
            value=ceil12(cache[name,kk]+shift+max(x*x-2*ell*x for x in u)+qmax(aa,bb-lam,F(0),*r))
        assert value==F(cell["cap"]),(kind,name,k,value,cell["cap"])
        assert cell["target_pass"]==(value<F(40,71)) and cell["fallback_pass"]==(value<F(451,800))
        maximum=max(maximum,value);winner_counts[kind+":"+name]+=1
    assert seen==set(expected)
    assert maximum==F(d["maximum_cell_cap"])<F(d["strict_rational_cap"])
    out={"status":"PASS_EXACT_CLIPPED_CONTACT_GLOBAL_MIDDLE_COVER",
         "scope":"Complete physical middle interval; lower and upper same-b branches complete the global norm.",
         "source":str(path.relative_to(ROOT)),"source_sha256":sha(path),"operator_b":str(b),
         "amplitude_interval":d["amplitude_interval"],"u_bands":len(ub),"kappa_bands":len(kb),
         "complete_groups":len(groups),"original_cells":len(expected),"accepted_price_cells":len(d["cells"]),
         "exact_physically_empty_cells":len(d["physically_empty_cells"]),"directed_inverse_caps":len(bounds),
         "contact_corner_groups":len(contact_checks),"winners":dict(winner_counts),
         "maximum_cell_cap":str(maximum),"strict_rational_cap":d["strict_rational_cap"],
         "cache_replay_status":cache_audit["status"],"convexity":convexity}
    output=Path(__file__).with_name("clipped_contact_envelope_audit.json");output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("path",nargs="?",default="tonight/recovery/profile_cover/global_middle_b7171_critical_ninebank_contact_envelope.json")
    args=parser.parse_args();main(ROOT/args.path)
