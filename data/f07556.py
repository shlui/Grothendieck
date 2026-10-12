"""Independent expanded raw-Hermite price replay for the accepted middle box.

No producer price code is imported.  The rank-one profile is represented
as E=ppT in the already independently derived full-even price identity.
This verifies the stated source box only, not a global operator norm.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json, os, sys

os.environ.setdefault("GR_MPFR_BITS", "768")
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tonight/operator/projection/certificate"))
from mpfr_interval import I, UP
from audit_general_r3_parameters_tail import gates
from audit_general_r3_source_bands import price

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ledgerpath=ROOT/"tonight/operator/projection/targeted_certificate/middle_source_box.json"
    acceptedpath=ROOT/"tonight/audit/targeted_middle_R3_recheck/ACCEPTANCE.json"
    d=json.loads(ledgerpath.read_text());accepted=json.loads(acceptedpath.read_text())
    assert d["status"]=="PASS_DIRECTED_CONTINUOUS_MIDDLE_R3_SOURCE_BOX"
    assert accepted["status"]=="PASS_INDEPENDENT_COMPLETE_TARGETED_MIDDLE_R3_SCALAR_GATE"
    profile=ROOT/d["profile"];assert sha(profile)==accepted["profile_sha256"]
    assert accepted["profile_path"]==d["profile"]
    p={k:list(map(F,v)) if isinstance(v,list) else F(v)
       for k,v in json.loads(profile.read_text())["profile"].items()}
    assert p["eps"]==F(d["operator_b"])==F(7171,10000)
    p.update(E00=p["p0"]**2,E02=p["p0"]*p["p2"],E22=p["p2"]**2)
    gate_rows=gates(p)
    reserves=tuple(map(F,(d["reserve0"],d["reserve6"])))
    assert reserves==tuple(map(F,accepted["reserves"]))==(F(1,10000),F(1,1000000))
    assert accepted["q_cover"]==[["0","8"],["8","16"],["16","infinity"]]
    for report in accepted["reports"].values():
        assert sha(ROOT/report["path"])==report["sha256"]
        assert json.loads((ROOT/report["path"]).read_text())["status"]==report["status"]
    ends=[tuple(map(F,d[k])) for k in ("u","r","kappa","shape")]
    assert ends==[(F(149,200),F(3,4)),(F(133,500),F(67,250)),
                  (F(4,5),F(1)),(F(2,5),F(1))]
    A=p["major"][3]+reserves[1]+p["zeta"]**2
    B=p["G"]-p["d0"]-p["zeta"]**2
    assert A>=0 and B>=0
    assert (I(194)/I(10).sqrt()-I(F(633,9)*F(3,5))).above(0)
    ar=p["eps"]-p["E00"]+p["E02"]**2/p["E22"]
    assert ar==p["eps"]>0
    todo=set(product(*ends));rows=[];largest=F(0)
    for v in d["vertices"]:
        point=tuple(F(v[k]) for k in ("u","r","kappa","shape"))
        assert point in todo;todo.remove(point)
        upper=F(price(p,*point,*reserves).hi.text(UP,240))
        assert upper<=F(v["upper"])
        rows.append({"u":str(point[0]),"r":str(point[1]),"kappa":str(point[2]),"shape":str(point[3]),
                     "independent_upper":str(upper),"saved_upper":v["upper"]})
        largest=max(largest,upper)
    assert not todo and len(rows)==16
    cap=F(d["cap"]);assert largest<=F(d["maximum_directed_upper"])<cap<F(40,71)
    out={"status":"PASS_INDEPENDENT_TARGETED_MIDDLE_CONTINUOUS_SOURCE_BOX",
         "scope":"Only u[.745,.75],r[.266,.268],kappa[.8,1],all shapes; no complete amplitude or global norm claim.",
         "operator_b":str(p["eps"]),"profile":d["profile"],"profile_sha256":sha(profile),
         "ledger_sha256":sha(ledgerpath),"whole_q_acceptance_sha256":sha(acceptedpath),
         "reserve0":str(reserves[0]),"reserve6":str(reserves[1]),"corners":16,
         "maximum_independent_upper":str(largest),"cap":str(cap),"six_gates":gate_rows,
         "convexity":"Positive r quadratic and separately convex u,kappa,shape; complete box maximum is bounded by the sixteen corners.",
         "vertices":rows}
    output=Path(__file__).with_name("targeted_middle_source_box_audit.json")
    output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ("status","profile_sha256","corners","cap")},indent=2))

if __name__=="__main__":main()
