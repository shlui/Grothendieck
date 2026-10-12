"""Exact independent envelope algebra and complete physical-source audit.

Local Gaussian/Stein inequalities and the directed cached corner constants
are audited separately.  No producer module is imported.  This script
reconstructs the transported quadratic coefficients from exact profiles,
checks every option/minimum, and proves exhaustive continuous coverage.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
B = F(7171, 10000)
RHO = F(847, 1000)
PROFILES = {
    "rank68": "tonight/operator/projection/lower_amp_certificate/profile_u.68.json",
    "rank70": "tonight/operator/projection/lower_amp_certificate/profile_u.70.json",
    "full68": "tonight/operator/projection/lower_amp_certificate/general_profile_u68_high_exact.json",
    "full70": "tonight/operator/projection/lower_amp_certificate/general_profile_u70_high_exact.json",
    "critical755": "tonight/operator/projection/certificate/frozen_R3_profile_exact.json",
    "diffuse7478": "tonight/operator/projection/diffuse_certificate/profile_interval_u7478_k0_01.json",
    "diffuse7505": "tonight/operator/projection/diffuse_certificate/profile_interval_b7171_newworst.json",
}
RESERVES = {name: (F(1, 1000), F(1, 10**6)) if name.startswith("full") else
            (F(1, 50000), F(1, 10**6) if name == "critical755" else F(1, 5000000))
            for name in PROFILES}


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def iv(x):
    a, b = map(F, x)
    assert a < b
    return a, b


def quadratic_max(a, c, lo, hi):
    xs = [lo, hi]
    if a < 0 and lo < -c/(2*a) < hi:
        xs.append(-c/(2*a))
    return max(a*x*x + c*x for x in xs)


def ceil12(x):
    return F(-((-x.numerator*10**12)//x.denominator), 10**12)


def main(path):
    d = read(path)
    assert d["status"] == "PASS_CONTINUOUS_ALL_SOURCE_NATIVE_R3_ENVELOPE_PRICING"
    assert F(d["operator_b"]) == B
    lo, hi = iv(d["amplitude_interval"])
    proof_binding = {str(Path(t["path"]).resolve()): t["sha256"] for t in d["accepted_R3_binding"]}
    for p, h in proof_binding.items():
        assert sha(Path(p)) == h
    coeff = {}
    convexity = {}
    for name, rel in PROFILES.items():
        pp = ROOT / rel
        assert proof_binding[str(pp)] == sha(pp)
        p = read(pp)["profile"]
        q = {k: F(v) for k, v in p.items() if not isinstance(v, list)}
        assert 0 < B <= q["eps"]
        E00, E02, E22 = ((q["E00"], q["E02"], q["E22"]) if "E00" in q else
                        (q["p0"]**2, q["p0"]*q["p2"], q["p2"]**2))
        assert E22 > 0
        yc = q["y"] - E02*q["r"]
        ar = B - E00 + E02**2/E22
        br = -2*q["gamma"] + 2*yc*E02/E22
        coeff[name] = (ar, br, q["ell"])
        stored = d["profile_coefficients"][name]
        assert ar == F(stored["r_quadratic"])
        assert br == F(stored["r_linear"])
        assert q["ell"] == F(stored["ell"])
        reserve0, reserve6 = RESERVES[name]
        A = F(p["major"][3]) + reserve6 + q["zeta"]**2
        C = q["G"] - q["d0"] - q["zeta"]**2
        assert A >= 0 and C >= 0
        # Price is separately convex in kappa and z: its nonlinear terms
        # are A*K(z)*kappa^(3/2)+C*kappa^2+(9C/4)*(kappa*z)^(3/2).
        # K(z)=194/sqrt(10)-(633/9)(z-2/5)>0 on z[2/5,1].
        assert F(194, 4) - F(633, 9)*F(3, 5) > 0
        convexity[name] = {"sixth_coefficient": str(A), "residual_coefficient": str(C),
                           "reserve0": str(reserve0), "reserve6": str(reserve6)}
    assert set(d["profile_coefficients"]) == set(PROFILES)
    cache = {}
    for t in d["cached_corner_constants"]:
        key = (t["profile"], iv(t["kappa_interval"]))
        assert key not in cache
        cache[key] = F(t["upper"])
    expected = {}
    summaries = {}
    native_records = []
    for binding in d["native_inputs"]:
        p = ROOT / binding["path"]
        assert sha(p) == binding["sha256"]
        n = read(p)
        assert F(n["operator_b"]) == B and F(n["source_rho"]) == RHO
        assert n["status"] == "PASS_EXACT_CONTINUOUS_NATIVE_FRAME_UNION_PRICING"
        for c in n["cells"]:
            ui = iv(c["u_interval"])
            if lo <= ui[0] < ui[1] <= hi:
                key = (ui, iv(c["r_interval"]), iv(c["kappa_interval"]))
                assert key not in expected
                expected[key] = c
        for s in n["summary"]:
            ui = iv(s["u_interval"])
            if lo <= ui[0] < ui[1] <= hi:
                key = (ui, iv(s["kappa_interval"]))
                assert key not in summaries
                summaries[key] = s
        native_records.append({"path": binding["path"], "sha256": binding["sha256"]})
    groups = defaultdict(list)
    ubands = sorted({ui for ui, _, _ in expected})
    kbands = sorted({ki for _, _, ki in expected})
    assert ubands[0][0] == lo and ubands[-1][1] == hi
    assert all(a[1] == b[0] for a, b in zip(ubands, ubands[1:]))
    assert kbands[0][0] == 0 and kbands[-1][1] == 90
    assert all(a[1] == b[0] for a, b in zip(kbands, kbands[1:]))
    assert set(summaries) == {(ui, ki) for ui in ubands for ki in kbands}
    assert set(cache) == {(name, ki) for name in PROFILES for ki in kbands}
    winners = Counter()
    maximum = F(0)
    seen = set()
    for cell in d["cells"]:
        ui, ri, ki = [iv(cell[k]) for k in ("u_interval", "r_interval", "kappa_interval")]
        key = (ui, ri, ki)
        assert key in expected and key not in seen
        seen.add(key)
        native = expected[key]
        assert cell["native_cap"] == native["cap"]
        options = []
        if native["cap"] is not None:
            options.append((F(native["cap"]), "native:" + native["frame"]))
        for name, (a, c, ell) in coeff.items():
            outer = max(x*x - 2*ell*x for x in ui)
            source = quadratic_max(a, c, *ri)
            cap = ceil12(cache[name, ki] + outer + source)
            options.append((cap, "R3:" + name))
        value, winner = min(options)
        assert value == F(cell["cap"]) and winner == cell["winner"]
        assert cell["target_pass"] == (value < F(40, 71))
        assert cell["fallback_pass"] == (value < F(451, 800))
        winners[winner] += 1
        maximum = max(maximum, value)
        groups[ui, ki].append(ri)
    assert seen == set(expected) and set(groups) == set(summaries)
    source_by_u = {}
    for (ui, ki), ranges in groups.items():
        ranges.sort()
        assert ranges[0][0] == 0
        assert all(a[1] == b[0] for a, b in zip(ranges, ranges[1:]))
        source = F(summaries[ui, ki]["source_upper"])
        assert ranges[-1][1] == source
        source_by_u.setdefault(ui[0], source)
        assert source_by_u[ui[0]] == source
    assert maximum == F(d["maximum_cell_cap"]) < F(d["strict_rational_cap"])
    assert F(d["resulting_reciprocal_lower_bound"]) == 1/F(d["strict_rational_cap"])
    out = {"status": "PASS_EXACT_ENVELOPE_ALGEBRA_AND_EXHAUSTIVE_PHYSICAL_SOURCE_COVERAGE",
           "scope": "Continuous regional middle bound; cached directed R3 corner constants and local/scalar primitive gates require their separate independent replays. A global norm also requires lower and upper amplitude branches.",
           "operator_b": str(B), "amplitude_interval": d["amplitude_interval"],
           "u_bands": len(ubands), "kappa_bands": len(kbands), "complete_groups": len(groups),
           "cell_count": len(seen), "maximum_cell_cap": str(maximum), "strict_rational_cap": d["strict_rational_cap"],
           "winners": dict(winners), "convexity_and_reserves": convexity,
           "source_upper_at_lower_u": {str(k): str(v) for k, v in sorted(source_by_u.items())},
           "native_inputs": native_records, "envelope_path": str(path.relative_to(ROOT)), "envelope_sha256": sha(path)}
    target = Path(__file__).with_name("updated_global_native_R3_envelope_algebra_audit.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in ["status", "cell_count", "complete_groups", "maximum_cell_cap", "strict_rational_cap"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--extra-profile", action="append", default=[], help="NAME=PROFILE_PATH,RESERVE0,RESERVE6")
    parser.add_argument("path", nargs="?", default="tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json")
    args = parser.parse_args()
    for spec in args.extra_profile:
        name, value = spec.split("=", 1)
        path, r0, r6 = value.rsplit(",", 2)
        assert name not in PROFILES
        PROFILES[name] = path
        RESERVES[name] = (F(r0), F(r6))
    main(ROOT / args.path)
