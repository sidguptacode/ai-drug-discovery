#!/usr/bin/env python3
"""
Validate Step 4b adjudication artifacts (schema + CSV/JSON consistency).

Step 5 requires these files. Does NOT choose labels — substantive adjudication is performed
by the reflecting agent following adjudication_artifacts.md in the ST-LR skill.

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

# Grounded explanations: rationale must cite evidence; CSV label_evidence must not be trivial.
MIN_RATIONALE_LEN = 40
MIN_TIER_RATIONALE_LEN = 15
MIN_LABEL_EVIDENCE_LEN = 25


def load_expected_clusters(out_dir: str) -> list[str] | None:
    """Cluster ids from step4_collated_candidates.json or step4_enrichr_marker_source.csv."""
    collated = os.path.join(out_dir, "step4_collated_candidates.json")
    if os.path.isfile(collated):
        try:
            with open(collated, encoding="utf-8") as f:
                doc = json.load(f)
        except json.JSONDecodeError:
            return None
        cl = doc.get("clusters")
        if not isinstance(cl, list):
            return None
        out: list[str] = []
        for obj in cl:
            if isinstance(obj, dict):
                c = str(obj.get("cluster", "")).strip()
                if c:
                    out.append(c)
        return out or None
    src = os.path.join(out_dir, "step4_enrichr_marker_source.csv")
    if not os.path.isfile(src):
        return None
    try:
        with open(src, encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None or "cluster" not in reader.fieldnames:
                return None
            return sorted({str(r.get("cluster", "")).strip() for r in reader if str(r.get("cluster", "")).strip()})
    except OSError:
        return None


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
        "label_evidence",
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
                for k in ("override_applied", "tier", "chosen_label", "rationale", "tier_rationale"):
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

    expected = load_expected_clusters(out_dir)
    if expected is not None:
        for cl in expected:
            if cl not in csv_by:
                errs.append("cluster %s from step 4 collation/marker source missing in CSV" % cl)

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
            tier = (r.get("tier") or "").strip()
            if tier not in TIER_OK:
                errs.append("cluster %s invalid tier in CSV: %r" % (cl, tier))

        for cl, jo in by_c.items():
            if cl not in csv_by:
                errs.append("cluster %s in JSON but not in CSV" % cl)

    for cl, r in csv_by.items():
        if not (r.get("adjudication_label") or "").strip():
            errs.append("cluster %s empty adjudication_label (required for step 5)" % cl)
        lev = (r.get("label_evidence") or "").strip()
        if len(lev) < MIN_LABEL_EVIDENCE_LEN:
            errs.append(
                "cluster %s label_evidence too short or empty (need >=%d chars citing markers/terms)"
                % (cl, MIN_LABEL_EVIDENCE_LEN)
            )

    if isinstance(rep, dict) and isinstance(rep.get("clusters"), list):
        for cl, jo in by_c.items():
            rat = str(jo.get("rationale") or "").strip()
            if len(rat) < MIN_RATIONALE_LEN:
                errs.append(
                    "cluster %s rationale too short or empty (need >=%d chars; cite markers/EnrichR)"
                    % (cl, MIN_RATIONALE_LEN)
                )
            tr = str(jo.get("tier_rationale") or "").strip()
            if len(tr) < MIN_TIER_RATIONALE_LEN:
                errs.append(
                    "cluster %s tier_rationale too short or empty (need >=%d chars)" % (cl, MIN_TIER_RATIONALE_LEN)
                )
            kmg = jo.get("key_marker_genes")
            cand = jo.get("candidates")
            n_kmg = len(kmg) if isinstance(kmg, list) else 0
            n_cand = len(cand) if isinstance(cand, list) else 0
            if n_kmg == 0 and n_cand == 0:
                errs.append(
                    "cluster %s needs key_marker_genes (>=1 gene) and/or candidates (>=1 enrichr row)"
                    % cl
                )
            elif n_kmg > 0:
                ok_gene = False
                for x in kmg:
                    if isinstance(x, str) and x.strip():
                        ok_gene = True
                        break
                if not ok_gene:
                    errs.append("cluster %s key_marker_genes must list at least one non-empty gene symbol" % cl)

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
