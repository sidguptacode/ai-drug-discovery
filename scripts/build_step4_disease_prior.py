#!/usr/bin/env python3
"""
Build step4_disease_marker_prior.json for step 4 marker filtering.

Combines:
  - Open Targets disease match + therapeutic areas + target–disease associations
  - Optional agent-provided tissues + canonical cell types (YAML or JSON — no LLM)
  - Optional agent-curated extra gene symbols merged after OT associations

Tissues and cell types MUST be supplied by the controller agent (merged config or
disease_prior_context_json). The script does not infer them.

Cap: default 100 genes total (annotation.disease_prior_max_genes).

Run from run_step_4.sh or:
  STEP4_PRIOR_FORCE=1 python3 scripts/build_step4_disease_prior.py --config config.yml
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore

GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DEFAULT_OUT_NAME = "step4_disease_marker_prior.json"
PAGE_SIZE = 100
SCHEMA_VERSION = 3


def graphql(query: str) -> dict:
    body = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(
        GRAPHQL_URL,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode("utf-8"))


def search_disease(query_string: str, n: int = 5) -> list[dict]:
    q = query_string.replace("\\", "\\\\").replace('"', '\\"')
    gql = (
        "query { search(queryString: \"%s\", entityNames: [\"disease\"], "
        "page: {index: 0, size: %d}) { hits { id name } } }" % (q, n)
    )
    data = graphql(gql)
    err = data.get("errors")
    if err:
        raise RuntimeError("Open Targets search error: %s" % err)
    hits = (data.get("data") or {}).get("search") or {}
    return list(hits.get("hits") or [])


def fetch_disease_context(efo_id: str) -> dict[str, Any]:
    gql = (
        "query { disease(efoId: \"%s\") { id name description "
        "therapeuticAreas { id name } } }" % efo_id.replace('"', '\\"')
    )
    data = graphql(gql)
    err = data.get("errors")
    if err:
        raise RuntimeError("Open Targets disease query error: %s" % err)
    dis = (data.get("data") or {}).get("disease")
    if not dis:
        return {}
    return dis


def fetch_associated_symbols(efo_id: str, max_genes: int) -> tuple[list[str], int]:
    symbols: list[str] = []
    seen: set[str] = set()
    index = 0
    total_count = 0
    while len(symbols) < max_genes:
        gql = (
            "query { disease(efoId: \"%s\") { associatedTargets(page: {index: %d, size: %d}) "
            "{ count rows { target { approvedSymbol } score } } } }"
            % (efo_id, index, PAGE_SIZE)
        )
        data = graphql(gql)
        err = data.get("errors")
        if err:
            raise RuntimeError("Open Targets disease query error: %s" % err)
        dis = (data.get("data") or {}).get("disease")
        if not dis:
            break
        at = dis.get("associatedTargets") or {}
        total_count = int(at.get("count") or 0)
        rows = at.get("rows") or []
        if not rows:
            break
        for row in rows:
            sym = (row.get("target") or {}).get("approvedSymbol")
            if not sym:
                continue
            u = sym.upper()
            if u not in seen:
                seen.add(u)
                symbols.append(sym)
                if len(symbols) >= max_genes:
                    break
        index += 1
        if len(rows) < PAGE_SIZE:
            break
    return symbols, total_count


def load_config_simple(path: str) -> dict:
    keys = ("out_dir", "disease", "species", "dataset_name")
    out: dict = {}
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line or line.startswith("---"):
                continue
            if ":" not in line:
                continue
            key, _, rest = line.partition(":")
            key = key.strip()
            if key not in keys:
                continue
            val = rest.strip()
            if not val or val.startswith("["):
                continue
            if (val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'")):
                val = val[1:-1]
            out[key] = val
    return out


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if yaml is not None:
        return yaml.safe_load(text)
    low = path.lower()
    if low.endswith((".yml", ".yaml")):
        print(
            "build_step4_disease_prior: PyYAML not installed; nested keys under "
            "`annotation:` are ignored. Use conda env ai-drug-discovery or install PyYAML.",
            file=sys.stderr,
        )
    return load_config_simple(path)


def _annotation_block(cfg: dict) -> dict:
    ann = cfg.get("annotation")
    if isinstance(ann, dict):
        return ann
    return {}


def _as_str_list(x: Any) -> list[str]:
    if x is None:
        return []
    if isinstance(x, str):
        return [s.strip() for s in x.split(",") if s.strip()]
    if isinstance(x, (list, tuple)):
        return [str(i).strip() for i in x if str(i).strip()]
    return []


def _normalize_gene_symbols(raw: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for g in raw:
        g = str(g).strip()
        if not g:
            continue
        if not re.match(r"^[A-Za-z0-9][A-Za-z0-9\-]*$", g):
            continue
        u = g.upper()
        if u in seen:
            continue
        seen.add(u)
        out.append(g)
    return out


def resolve_context_path(rel_or_abs: str, out_dir: str, config_path: str) -> str | None:
    if not rel_or_abs or not str(rel_or_abs).strip():
        return None
    p = str(rel_or_abs).strip()
    if os.path.isabs(p) and os.path.isfile(p):
        return p
    config_dir = os.path.dirname(os.path.abspath(config_path))
    for base in (out_dir, config_dir, os.getcwd()):
        cand = os.path.join(base, p)
        if os.path.isfile(cand):
            return cand
    return None


def load_agent_disease_context(
    ann: dict,
    out_dir: str,
    config_path: str,
) -> tuple[list[str], list[str], list[str], str | None]:
    """
    Tissues, cell types, extra genes from:
      - annotation.disease_prior_context_json (optional JSON file)
      - annotation.disease_prior_tissues / disease_prior_cell_types / disease_prior_extra_genes
    JSON may contain: tissues, cell_types, extra_genes (aliases: prior_tissues, prior_cell_types).
    """
    tissues: list[str] = []
    ctypes: list[str] = []
    extra_genes: list[str] = []
    note: str | None = None

    jpath = ann.get("disease_prior_context_json")
    if jpath:
        resolved = resolve_context_path(str(jpath), out_dir, config_path)
        if resolved and os.path.isfile(resolved):
            try:
                with open(resolved, encoding="utf-8") as f:
                    doc = json.load(f)
                if isinstance(doc, dict):
                    tissues.extend(_as_str_list(doc.get("tissues") or doc.get("prior_tissues")))
                    ctypes.extend(_as_str_list(doc.get("cell_types") or doc.get("prior_cell_types")))
                    extra_genes.extend(_normalize_gene_symbols(_as_str_list(doc.get("extra_genes"))))
            except (OSError, json.JSONDecodeError) as e:
                note = "Could not read disease_prior_context_json: %s" % e
        else:
            note = "disease_prior_context_json path not found: %s" % jpath

    tissues.extend(_as_str_list(ann.get("disease_prior_tissues")))
    ctypes.extend(_as_str_list(ann.get("disease_prior_cell_types")))
    extra_genes.extend(_normalize_gene_symbols(_as_str_list(ann.get("disease_prior_extra_genes"))))

    # de-dup preserve order
    def dedupe_str(xs: list[str]) -> list[str]:
        seen: set[str] = set()
        out: list[str] = []
        for a in xs:
            k = a.strip().lower()
            if not k or k in seen:
                continue
            seen.add(k)
            out.append(a.strip())
        return out

    return dedupe_str(tissues), dedupe_str(ctypes), dedupe_str(extra_genes), note


def merge_gene_lists(ot_genes: list[str], agent_genes: list[str], cap: int) -> list[str]:
    """OT genes first (association rank), then agent extra genes not already present."""
    seen: set[str] = set()
    ordered: list[str] = []
    for g in ot_genes:
        u = g.upper()
        if u in seen:
            continue
        seen.add(u)
        ordered.append(g)
        if len(ordered) >= cap:
            return ordered
    for g in agent_genes:
        u = g.upper()
        if u in seen:
            continue
        seen.add(u)
        ordered.append(g)
        if len(ordered) >= cap:
            break
    return ordered


def main() -> int:
    ap = argparse.ArgumentParser(description="Build step4 disease marker prior JSON.")
    ap.add_argument("--config", required=True, help="Pipeline YAML (merged or base config)")
    ap.add_argument("--force", action="store_true", help="Overwrite existing prior file")
    ap.add_argument(
        "--max-genes",
        type=int,
        default=None,
        help="Override max genes cap (default: annotation.disease_prior_max_genes or 100)",
    )
    args = ap.parse_args()
    if os.environ.get("STEP4_PRIOR_FORCE", "").strip() in ("1", "true", "yes"):
        args.force = True

    cfg = load_config(args.config)
    ann = _annotation_block(cfg)
    cap = args.max_genes
    if cap is None:
        cap = ann.get("disease_prior_max_genes")
    if cap is None:
        cap = 100
    cap = int(cap)
    if cap < 1:
        cap = 100

    out_dir = cfg.get("out_dir")
    if not out_dir:
        print("build_step4_disease_prior: config missing out_dir; skip.", file=sys.stderr)
        return 0
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, DEFAULT_OUT_NAME)
    if os.path.isfile(out_path) and not args.force:
        print("build_step4_disease_prior: %s exists; skip (use --force to rebuild)." % out_path)
        return 0

    disease = (cfg.get("disease") or "").strip()
    species = (cfg.get("species") or "").strip()
    if not disease:
        doc = {
            "schema_version": SCHEMA_VERSION,
            "disease_query": "",
            "disease_matched_id": None,
            "disease_matched_name": None,
            "species": species or None,
            "disease_prior_tissues": [],
            "disease_prior_cell_types": [],
            "open_targets_therapeutic_areas": [],
            "source": None,
            "source_url": "https://platform.opentargets.org",
            "genes": [],
            "association_target_count": 0,
            "built_at": datetime.now(timezone.utc).isoformat(),
            "notes": "No disease string in config; add cfg disease: \"...\" or hand-edit genes in this file.",
        }
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=2)
            f.write("\n")
        print("build_step4_disease_prior: wrote empty prior (no disease in config): %s" % out_path)
        return 0

    try:
        hits = search_disease(disease, n=5)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, RuntimeError) as e:
        print("build_step4_disease_prior: Open Targets search failed: %s" % e, file=sys.stderr)
        return 1
    if not hits:
        print("build_step4_disease_prior: no disease match for %r; not writing prior." % disease, file=sys.stderr)
        return 1

    top = hits[0]
    efo_id = top["id"]
    efo_name = top.get("name") or ""

    dis_ctx: dict[str, Any] = {}
    try:
        dis_ctx = fetch_disease_context(efo_id)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, RuntimeError) as e:
        print("build_step4_disease_prior: disease context fetch failed: %s" % e, file=sys.stderr)

    description = (dis_ctx.get("description") or "").strip()
    ta_list = dis_ctx.get("therapeuticAreas") or []
    therapeutic_areas = []
    ta_names: list[str] = []
    for ta in ta_list:
        if isinstance(ta, dict):
            therapeutic_areas.append({"id": ta.get("id"), "name": ta.get("name")})
            if ta.get("name"):
                ta_names.append(str(ta["name"]))

    agent_tissues, agent_cell_types, agent_extra_genes, ctx_note = load_agent_disease_context(
        ann, out_dir, args.config
    )

    reserved_agent = min(40, max(0, cap // 2))
    ot_budget = max(1, cap - reserved_agent)

    try:
        ot_genes, assoc_count = fetch_associated_symbols(efo_id, max_genes=ot_budget)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, RuntimeError) as e:
        print("build_step4_disease_prior: associations fetch failed: %s" % e, file=sys.stderr)
        ot_genes, assoc_count = [], 0

    ot_genes = ot_genes[:ot_budget]

    final_genes = merge_gene_lists(ot_genes, agent_extra_genes, cap)
    ot_set = {g.upper() for g in ot_genes}
    merge_info = {
        "open_targets_genes_requested": ot_budget,
        "open_targets_genes_in_final": sum(1 for g in final_genes if g.upper() in ot_set),
        "agent_extra_genes_in_final": sum(1 for g in final_genes if g.upper() not in ot_set),
        "cap": cap,
    }

    notes_parts = [
        "Prior genes: Open Targets associations first, then optional agent-curated symbols "
        "(annotation.disease_prior_extra_genes or disease_prior_context_json). "
        "Tissues and cell types are recorded from the agent (YAML/JSON), not inferred here. "
        "step_4.R keeps marker rows whose gene symbol matches (case-insensitive)."
    ]
    if not agent_tissues and not agent_cell_types:
        notes_parts.append(
            "No disease_prior_tissues / disease_prior_cell_types (or context JSON) set — "
            "set these in merged config or add step4_disease_prior_context.json under OUT_DIR."
        )
    if ctx_note:
        notes_parts.append(ctx_note)

    doc = {
        "schema_version": SCHEMA_VERSION,
        "disease_query": disease,
        "disease_matched_id": efo_id,
        "disease_matched_name": efo_name,
        "disease_description": description or None,
        "species": species or None,
        "disease_prior_tissues": agent_tissues,
        "disease_prior_cell_types": agent_cell_types,
        "open_targets_therapeutic_areas": therapeutic_areas,
        "source": "Open Targets Platform (target–disease associations) + agent context (tissues/cell types/extra genes)",
        "source_url": "https://platform.opentargets.org/disease/%s/associations" % efo_id,
        "genes": final_genes,
        "association_target_count": assoc_count,
        "genes_capped_at": cap,
        "gene_merge": merge_info,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "notes": " ".join(notes_parts),
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(
        "build_step4_disease_prior: wrote %d genes (OT %d + agent extras) -> %s"
        % (len(final_genes), len(ot_genes), out_path)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
