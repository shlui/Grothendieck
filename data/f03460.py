#!/usr/bin/env python3
"""Replay only the lower amplitude R3 remainder used by the global theorem."""
from pathlib import Path
import importlib.util
import json

ROOT=Path(__file__).resolve().parent
path=ROOT/'tonight/audit/diffuse_b7171_R3/independent_source_pricing.py'
spec=importlib.util.spec_from_file_location('independent_lower_R3_price',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
result=module.audit_remainder_box()
print(json.dumps({k:result[k]for k in ('status','operator_b','maximum_upper','scope')},indent=2))
