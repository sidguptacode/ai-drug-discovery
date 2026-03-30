#!/usr/bin/env python3
"""
Validate Step 4b adjudication artifacts (schema + CSV/JSON consistency).

Does NOT choose labels — substantive adjudication is performed by the controller agent
following .cursor/skills/st-lr-bioinformatics-pipeline/adjudication_artifacts.md.

Usage:
  python3 scripts/validate_step4b_artifacts.py --out-dir outputs/<run_id>
  python3 scripts/validate_step4b_artifacts.py --config runs/<id>/config_merged.yml

Optional:
  --forbid-regex PATTERN   Repeatable; any adjudication_label matching any pattern fails.
  --require-approval       Require OUT_DIR/step4_adjudication_approved.flag to exist (non-empty optional).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore


TIER_OK = frozenset(
    {
        "tier1_strong_consensus",
        "tier2_partial_consensus",
        "tier3_conflicted",
    }
)


def load_out_dir_from_config(path: str) -> str | None:
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if yaml is not None:
        cfg = yaml.safe_load(text)
    else:
        cfg = {}
        for raw in text.splitlines():
            line = raw.split("#", 1)[0].strip()
            if line.startswith("out_dir:"):
                _, _, rest = line.partition(":")
                cfg["out_dir"] = rest.strip().strip('"').strip("'")
                break
    return cfg.get("out_dir") if isinstance(cfg, dict) else None


def parse_bool_cell(x: str) -> bool:
    v = str(x).strip().upper()
    return v in ("TRUE", "T", "1", "YES", "Y")


def main() -> int:
    ap = argparse.ArgumentParser(description="Validate step4 adjudication CSV/JSON.")
    ap.add_argument("--out-dir", help="Pipeline OUT_DIR (outputs/<run_id>)")
    ap.add_argument("--config", help="YAML with out_dir (merged or base config)")
    ap.add_argument(
        "--forbid-regex",
        action="append",
        default=[],
        help="Repeatable; adjudication_label must not match (case-insensitive)",
    )
    ap.add_argument(
        "--require-approval",
        action="store_true",
        help="Require step4_adjudication_approved.flag in OUT_DIR",
    )
    args = ap.parse_args()

    out_dir = args.out_dir
    if not out_dir and args.config:
        out_dir = load_out_dir_from_config(args.config)
    if not out_dir:
        print("validate_step4b: need --out-dir or --config with out_dir", file=sys.stderr)
        return 2

    labels_path = os.path.join(out_dir, "step4_adjudication_labels.csv")
    report_path = os.path.join(out_dir, "step4_adjudication_report.json")

    errs: list[str] = []

    if not os.path.isfile(labels_path):
        errs.append("missing: step4_adjudication_labels.csv")
    if not os.path.isfile(report_path):
        errs.append("missing: step4_adjudication_report.json")

    if args.require_approval:
        flag = os.path.join(out_dir, "step4_adjudication_approved.flag")
        if not os.path.isfile(flag):
            errs.append("missing: step4_adjudication_approved.flag (--require-approval)")

    if errs:
        for e in errs:
            print("validate_step4b: %s" % e, file=sys.stderr)
        return 1

    req_csv = (
        "cluster",
        "original_label",
        "adjudication_label",
        "override_applied",
        "tier",
        "confidence",
    )

    with open(labels_path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            errs.append("CSV: no header row")
            rows = []
        else:
            miss = [c for c in req_csv if c not in reader.fieldnames]
            if miss:
                errs.append("CSV missing columns: %s" % ", ".join(miss))
            rows = list(reader)

    with open(report_path, encoding="utf-8") as f:
        try:
            rep: dict[str, Any] = json.load(f)
        except json.JSONDecodeError as e:
            errs.append("JSON parse error: %s" % e)
            rep = {}

    if isinstance(rep, dict):
        clusters = rep.get("clusters")
        if not isinstance(clusters, list):
            errs.append("JSON missing or invalid key: clusters (array required)")
            clusters = []
        if isinstance(clusters, list):
            by_c: dict[str, dict] = {}
            for i, obj in enumerate(clusters):
                if not isinstance(obj, dict):
                    errs.append("JSON clusters[%d] not an object" % i)
                    continue
                cl = str(obj.get("cluster", "")).strip()
                if not cl:
                    errs.append("JSON clusters[%d] missing cluster" % i)
                    continue
                by_c[cl] = obj
                for k in ("override_applied", "tier", "chosen_label"):
                    if k not in obj:
                        errs.append("JSON cluster %s missing %s" % (cl, k))
                if obj.get("tier") not in TIER_OK:
                    errs.append("JSON cluster %s invalid tier: %r" % (cl, obj.get("tier")))
    else:
        errs.append("JSON root must be an object")
        by_c = {}

    csv_by: dict[str, dict[str, str]] = {}
    for r in rows:
        cl = str(r.get("cluster", "")).strip()
        if not cl:
            errs.append("CSV row with empty cluster")
            continue
        csv_by[cl] = r

    if isinstance(rep, dict) and isinstance(rep.get("clusters"), list):
        for cl, r in csv_by.items():
            jo = by_c.get(cl)
            if jo is None:
                errs.append("cluster %s in CSV but not in JSON clusters[]" % cl)
                continue
            oa_csv = parse_bool_cell(r.get("override_applied", ""))
            oa_j = bool(jo.get("override_applied"))
            if oa_csv != oa_j:
                errs.append("cluster %s override_applied mismatch CSV=%s JSON=%s" % (cl, oa_csv, oa_j))
            adj = (r.get("adjudication_label") or "").strip()
            ch = str(jo.get("chosen_label") or "").strip()
            if adj != ch:
                errs.append("cluster %s adjudication_label != chosen_label" % cl)
            if oa_csv and not adj:
                errs.append("cluster %s override_applied true but empty adjudication_label" % cl)
            tier = (r.get("tier") or "").strip()
            if tier not in TIER_OK:
                errs.append("cluster %s invalid tier in CSV: %r" % (cl, tier))

        for cl, jo in by_c.items():
            if cl not in csv_by:
                errs.append("cluster %s in JSON but not in CSV" % cl)

    compiled: list[re.Pattern[str]] = []
    for pat in args.forbid_regex:
        try:
            compiled.append(re.compile(pat, re.I))
        except re.error as e:
            errs.append("invalid --forbid-regex %r: %s" % (pat, e))

    for cl, r in csv_by.items():
        lab = (r.get("adjudication_label") or "").strip()
        for cre in compiled:
            if cre.search(lab):
                errs.append("cluster %s adjudication_label matches forbidden regex %r" % (cl, cre.pattern))

    if errs:
        print("validate_step4b: FAILED (%d issue(s)):" % len(errs), file=sys.stderr)
        for e in errs:
            print("  - %s" % e, file=sys.stderr)
        return 1

    print("validate_step4b: OK — %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
