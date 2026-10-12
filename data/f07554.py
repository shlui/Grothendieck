"""Independent exact coverage audit, not a re-run of scalar integrations.

The local/scalar leaf inequalities and directed endpoint prices have separate
independent replays.  This script checks their frozen hash bindings, complete
u/kappa/source partitions, same-operator transport, and R3 remainder inclusion.
Run from the project root.  Only Python's standard library is used.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "tonight/recovery/profile_cover"
B = F(7171, 10000)
TARGET = F(40, 71)
RHO = F(847, 1000)
OUT = Path(__file__).with_name("same_b_native_R3_glue_audit.json")


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval(x):
    lo, hi = map(F, x)
    assert lo < hi
    return lo, hi


def check_binding(record):
    path = ROOT / record["path"]
    assert sha(path) == record["sha256"], str(path)
    return read(path)


def check_mask_bindings(native):
    result = []
    for binding in native["bindings"]:
        mask_path = ROOT / binding["mask"]
        package = ROOT / binding["package"]
        assert sha(mask_path) == binding["mask_sha256"]
        mask = read(mask_path)
        assert mask["source_sha256"] == binding["source_sha256"]
        assert sha(package / "data/source_profile.json") == binding["source_sha256"]
        assert F(mask["operator_b"]) == B
        assert F(mask["original_operator_b"]) == F(18, 25)
        assert F(mask["analytic_L1_rho"]) == RHO
        # The old .849 correction is an outward tail payment, not the
        # physical source restriction employed by the union.
        assert F(mask["source_bridge_rho"]) == F(849, 1000)
        local = read(package / "data/local_summary.json")
        p = local["profile"]
        beta = F(p["beta"])
        assert F(mask["local_constant_transport"]) == beta**2 * (1/B - F(25, 18))
        M = [[F(x) for x in row] for row in p["M"]]
        g = (F(p["nu"]) - M[0][2])**2 / M[2][2] - M[0][0] + B - 2*F(p["delta"])
        assert F(mask["exact_gram_cost"]) == g
        assert F(mask["ell"]) == F(p["ell"])
        assert local["status"] == "PASS"
        assert local["source_sha256"] == binding["source_sha256"]
        replay = read(package / "replay.json")
        assert replay["status"] in {"PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES", "PASS_COMPLETE_REGIONAL_BOUND"}
        endpoint_bands = [interval(e["kappa_interval"]) for e in mask["endpoint_prices"]]
        assert endpoint_bands[0][0] == 0 and endpoint_bands[-1][1] == 90
        assert all(a[1] == b[0] for a, b in zip(endpoint_bands, endpoint_bands[1:]))
        for alt in mask["alternative_scalar_bindings"]:
            ap = ROOT / alt["package"]
            assert sha(ap / "data/source_profile.json") == alt["source_sha256"]
            assert read(ap / "replay.json")["status"] == "PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES"
        result.append({"mask": binding["mask"], "package": binding["package"], "same_b_transport": "PASS"})
    return result


def check_partition(native, ulo, uhi, r3=None):
    assert native["status"] == "PASS_EXACT_CONTINUOUS_NATIVE_FRAME_UNION_PRICING"
    assert F(native["operator_b"]) == B
    assert F(native["threshold"]) == TARGET
    assert F(native["source_rho"]) == RHO
    ugrid = list(map(F, native["u_grid"]))
    assert ulo in ugrid and uhi in ugrid
    assert all(x < y for x, y in zip(ugrid, ugrid[1:]))
    ubands = [(a, b) for a, b in zip(ugrid, ugrid[1:]) if ulo <= a and b <= uhi]
    assert ubands[0][0] == ulo and ubands[-1][1] == uhi
    summaries = {}
    for s in native["summary"]:
        ui, ki = interval(s["u_interval"]), interval(s["kappa_interval"])
        if ui in ubands:
            assert (ui, ki) not in summaries
            summaries[ui, ki] = s
    kbands = sorted({ki for _, ki in summaries})
    assert kbands[0][0] == 0 and kbands[-1][1] == 90
    assert all(a[1] == b[0] for a, b in zip(kbands, kbands[1:]))
    assert set(summaries) == {(ui, ki) for ui in ubands for ki in kbands}
    groups = defaultdict(list)
    for cell in native["cells"]:
        ui, ki = interval(cell["u_interval"]), interval(cell["kappa_interval"])
        if ui in ubands:
            assert (ui, ki) in summaries
            groups[ui, ki].append(cell)
    assert set(groups) == set(summaries)
    n_native = n_r3 = 0
    caps = []
    source_upper_by_u = {}
    r3_box = {key + "_interval": interval(r3[key]) for key in ("u", "r", "kappa")} if r3 else None
    for key, s in summaries.items():
        cells = sorted(groups[key], key=lambda c: F(c["r_interval"][0]))
        endpoints = [interval(c["r_interval"]) for c in cells]
        assert endpoints[0][0] == 0
        assert all(a[1] == b[0] for a, b in zip(endpoints, endpoints[1:]))
        outward_source = F(s["source_upper"])
        assert endpoints[-1][1] == outward_source
        assert 0 < outward_source <= RHO
        # The directed source price is separately re-evaluated by the
        # independent endpoint-pricing audit.  Here we require its exact
        # equality across all cumulant bands at each lower u endpoint.
        ui = key[0]
        source_upper_by_u.setdefault(ui[0], outward_source)
        assert source_upper_by_u[ui[0]] == outward_source
        for c in cells:
            if c["accepted"]:
                assert c["frame"] in {x["mask"] for x in native["bindings"]}
                cap = F(c["cap"])
                assert cap < TARGET
                n_native += 1
            else:
                assert r3 is not None
                for field in r3_box:
                    a, b = interval(c[field]); lo, hi = r3_box[field]
                    assert lo <= a < b <= hi, (field, c)
                cap = F(r3["rounded_exact_cap"])
                assert cap < TARGET
                n_r3 += 1
            caps.append(cap)
    assert all(a >= b for a, b in zip(source_upper_by_u.values(), list(source_upper_by_u.values())[1:]))
    return {"amplitude_interval": list(map(str, (ulo, uhi))), "u_bands": len(ubands),
            "kappa_bands": len(kbands), "complete_product_groups": len(groups),
            "native_cells": n_native, "r3_cells": n_r3, "maximum_cell_cap": str(max(caps)),
            "source_upper_at_lower_u": {str(k): str(v) for k, v in source_upper_by_u.items()}}


def early_transfers():
    d = read(BASE / "native_b7171_regional_transfer_assembly.json")
    assert d["status"] == "PASS_EXACT_ACCEPTED_INPUT_NATIVE_OPERATOR_TRANSFER"
    assert F(d["operator_b"]) == B
    out = []
    for r in d["rows"]:
        if not r["passed_target"]:
            continue
        assert F(r["operator_b"]) == B and F(r["target"]) == TARGET
        assert F(r["directed_total_upper"]) < TARGET
        p = ROOT / r["package"]
        binding = r["bindings"]
        assert sha(p / "data/source_profile.json") == binding["source_sha256"]
        assert sha(p / "data/local_summary.json") == binding["local_summary_sha256"]
        assert sha(p / "replay.json") == binding["replay_sha256"]
        for n, h in binding["scalars"].items():
            assert sha(p / "data/scalars" / n) == h
        out.append({"amplitude_interval": r["amplitude_interval"], "package": r["package"],
                    "directed_total_display": r["directed_total_display"]})
    assert len(out) == 4
    intervals = sorted(interval(r["amplitude_interval"]) for r in out)
    right = F(7, 50)
    for lo, hi in intervals:
        assert lo <= right
        right = max(right, hi)
    assert right == F(319, 500)
    return out


def main():
    union = read(BASE / "native_R3_b7171_u682_6905_union.json")
    assert union["status"] == "PASS_COMPLETE_SAME_OPERATOR_NATIVE_R3_AMPLITUDE_UNION"
    assert F(union["operator_b"]) == B and F(union["target"]) == TARGET
    assert F(union["physical_source_rho"]) == RHO
    native = check_binding(union["native_binding"])
    r3 = check_binding(union["R3_binding"])
    assert r3["status"] == "PASS_ALL_Q_B7171_NATIVE_REMAINDER_BOX"
    assert F(r3["operator_b"]) == B
    assert list(map(F, r3["shape"])) == [F(2, 5), F(1)]
    assert F(r3["maximum_upper"]) < F(r3["rounded_exact_cap"]) < TARGET
    middle = check_partition(native, F(341, 500), F(1381, 2000), r3)
    assert middle["native_cells"] == union["native_cells"]
    assert middle["r3_cells"] == union["r3_cells"]
    assert F(middle["maximum_cell_cap"]) == F(union["maximum_rounded_cell_cap"])
    base = read(BASE / "native_lower_base_b7171_target_union.json")
    lower = check_partition(base, F(3, 5), F(341, 500))
    bindings = check_mask_bindings(native) + check_mask_bindings(base)
    early = early_transfers()
    # The elementary energy branch is conservative at b0=.72 and hence
    # valid at b=.7171 as well.  Its upper endpoint overlaps the first
    # transferred frame; the .60 base union overlaps the fourth.
    lowcap = F(541, 2500)**2 + F(18, 25)*RHO**2
    assert lowcap < TARGET and F(541, 2500) > F(7, 50)
    assert F(319, 500) > F(3, 5)
    out = {"status": "PASS_EXACT_ANALYTIC_COVERAGE_AND_FROZEN_BINDINGS",
           "scope": "Complete upper-pairing amplitude cover only, conditional on the separately independently replayed scalar/local primitives and endpoint prices; not a global norm bound.",
           "operator_b": str(B), "target": str(TARGET),
           "complete_amplitude_interval": ["0", "1381/2000"], "physical_L1_rho": str(RHO),
           "native_shape": ["0", "1"], "R3_shape": ["2/5", "1"],
           "exact_shape_map": "z=2/5+(3/5)s", "low_branch_cap": str(lowcap),
           "early_transfers": early, "lower_union": lower, "native_R3_union": middle,
           "mask_binding_checks": bindings,
           "frozen_inputs": {str(p.relative_to(ROOT)): sha(p) for p in [
               BASE / "native_b7171_regional_transfer_assembly.json", BASE / "native_lower_base_b7171_target_union.json",
               BASE / "native_R3_b7171_u682_6905_union.json", ROOT / union["native_binding"]["path"], ROOT / union["R3_binding"]["path"]]}}
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "lower_cells": lower["native_cells"],
                      "middle_cells": middle["native_cells"] + middle["r3_cells"], "output": str(OUT)}))


if __name__ == "__main__":
    main()
