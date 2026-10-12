#!/usr/bin/env python3
"""Independent targeted full-E whole-q audit with different group reserves.

The accepted producer-free 4D walker and directed Gaussian/Stein primitives
are reused.  Body and exterior are replayed at their OWN frozen majorants,
then dominated by one global sextic majorant.  The q>=16 tail is freshly
reconstructed with exact Fractions.  No producer-cover module is imported.
"""
from fractions import Fraction as F
from pathlib import Path
from multiprocessing import Pool
import argparse
import hashlib
import json
import time

import replay_general_r3_cover as R
from audit_general_r3_parameters_tail import load, gates, tail

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tonight/operator/projection/targeted_certificate"
GLOBAL_R0, GLOBAL_R6 = F(1, 200), F(1, 100000)


def audit_group(path, raw, workers):
    start = time.time()
    d, p, leaves, excluded = R.reconstruct(path)
    assert all(p[k] == raw[k] for k in p if k != "major")
    reserve0, reserve6 = F(d["reserve0"]), F(d["reserve6"])
    assert p["major"] == [raw["major"][j] + (reserve0 if j == 0 else reserve6 if j == 3 else 0)
                           for j in range(4)]
    assert 0 <= reserve0 <= GLOBAL_R0 and 0 <= reserve6 <= GLOBAL_R6
    N = R.build_matrix(p)
    keys, intervals = set(), set()
    for bx in leaves:
        center = tuple((a+b)/2 for a, b in bx[1:])
        intervals.add(bx[0])
        keys.update((q,)+center for q in bx[0])
    qnodes = {key[0] for key in keys}
    majors = {q: R.majorant(p, q) for q in qnodes}
    curvature = {ab: max(F(0), -R.second_lower(p, *ab)) for ab in intervals}
    print(json.dumps({"stage": "PARTITION_AND_CURVATURE_PASS", "path": path.name,
                      "leaves": len(leaves), "excluded": excluded, "endpoints": len(keys)}), flush=True)
    endpoints = {}
    with Pool(workers, initializer=R.init_worker, initargs=(p, N, majors)) as pool:
        for count, (key, value) in enumerate(pool.imap_unordered(R.directed_endpoint, keys, chunksize=100), 1):
            endpoints[key] = value
            if count % 4000 == 0:
                print(json.dumps({"stage": "DIRECTED_ENDPOINT_REPLAY", "path": path.name,
                                  "endpoints": count, "seconds": time.time()-start}), flush=True)
    C3 = abs(p["zeta"]*(-1+p["d0"]/(2*p["G"])))*(abs(p["p0"])+abs(p["p2"]))
    largest = F(-10**20)
    for bx in leaves:
        a, b = bx[0]
        center = tuple((x+y)/2 for x, y in bx[1:])
        width = [(y-x)/2 for x, y in bx[1:]]
        radius = sum(N[j][j]*width[j]**2 for j in range(4))
        radius += 2*abs(N[0][2])*width[0]*width[2] + 2*abs(N[1][3])*width[1]*width[3]
        qreserve = (F(3, 2)*b*C3 + curvature[a, b]/8)*(b-a)**2
        residual = max(endpoints[(a,)+center], endpoints[(b,)+center]) + radius + qreserve
        assert residual < 0, (path, bx, residual)
        largest = max(largest, residual)
    return {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "q_interval": [d["qa"], d["qb"]], "reserve0": str(reserve0), "reserve6": str(reserve6),
            "events": len(d["events"]), "leaves": len(leaves), "excluded": excluded,
            "directed_endpoints": len(keys), "q_nodes": len(qnodes), "largest_leaf_residual": str(largest),
            "global_majorant_dominance": "PASS", "seconds": time.time()-start}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    start = time.time()
    profile = DIRECTORY / "general_profile_highk_u7575_r259_k2_4.json"
    paths = [DIRECTORY / "general_highk_body.json.gz", DIRECTORY / "general_highk_exterior16_loose.json.gz"]
    sources = [profile, *paths, Path(__file__), Path(R.__file__),
               ROOT / "tonight/theory/audit_general_r3_parameters_tail.py"]
    sources += [R.CERT/name for name in ("mpfr_interval.py", "gaussian_max_certificate.py", "certify_local_point.py",
                                        "certify_local_bernstein_region.py", "rational_phi_bounds.py", "bernstein_integer_bounds.py")]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    raw = load(profile)
    assert raw["eps"] == F(7171, 10000)
    gate_rows = gates(raw)
    groups = [audit_group(path, raw, args.workers) for path in paths]
    assert groups[0]["q_interval"] == ["0", "8"] and groups[1]["q_interval"] == ["8", "16"]
    t = tail(raw, F(16), GLOBAL_R0, GLOBAL_R6)
    previous = json.loads((DIRECTORY / "general_highk_tail16_loose.json").read_text())
    pt = previous.get("tail", previous)
    assert pt["difference_polynomial"] == t["difference_polynomial"]
    assert pt["shifted_coefficients"] == t["shifted_coefficients"]
    assert all(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h for p, h in hashes.items())
    out = {"status": "PASS_INDEPENDENT_TARGETED_FULL_E_WHOLE_Q_PROFILE",
           "scope": "One frozen full-E operator profile valid on the entire q line; outer source/u/kappa coverage is separate.",
           "operator_b": str(raw["eps"]), "global_reserve0": str(GLOBAL_R0), "global_reserve6": str(GLOBAL_R6),
           "groups": groups, "events": sum(g["events"] for g in groups),
           "leaves": sum(g["leaves"] for g in groups), "excluded": sum(g["excluded"] for g in groups),
           "six_gates": gate_rows, "tail": t, "frozen_source_hashes": hashes, "seconds": time.time()-start}
    output = Path(__file__).with_name("targeted_highk_whole_q_audit.json")
    output.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({k: out[k] for k in ("status", "events", "leaves", "excluded", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
