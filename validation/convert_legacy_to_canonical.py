#!/usr/bin/env python3
"""Audit the migration ledger; refuse conversion until reviewed adapters exist.

The first pass's generic column copy silently skipped missing branches,
selected the first alias and confused event rows with object rows. This entry
point validates the review contract without creating a misleading ROOT file.
Exact layout-specific identity, explode/join and validity adapters remain TODO.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

DISPOSITIONS = {"STORED", "DERIVED", "POST-TREE", "ARCHIVE", "REVIEW"}


def audit_ledger(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    problems, seen = [], set()
    if not rows:
        raise ValueError("empty field ledger")
    for line, row in enumerate(rows, 2):
        identity = (row["legacy_tree"], row["legacy_branch"])
        if identity in seen:
            problems.append({"line": line, "problem": "duplicate legacy field identity"})
        seen.add(identity)
        disposition = row["disposition"]
        required = ["evidence", "review_status"]
        if disposition == "STORED": required += ["canonical_table", "canonical_field", "validity_requirements"]
        if disposition == "DERIVED": required += ["source_fields", "formula", "validity_requirements"]
        if disposition == "POST-TREE": required += ["downstream_owner"]
        if disposition == "ARCHIVE": required += ["archive_reason"]
        if disposition not in DISPOSITIONS or disposition == "REVIEW":
            problems.append({"line": line, "problem": "unresolved disposition"})
        elif row.get("review_status") != "accepted" or any(not row.get(k, "").strip() for k in required):
            problems.append({"line": line, "problem": "missing accepted evidence/definition"})
    return {"schema": "MigrationLedgerAuditV1", "ledger": str(path.resolve()),
            "ledger_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "rows": len(rows),
            "dispositions": dict(Counter(r["disposition"] for r in rows)), "problems": problems,
            "ledger_complete": not problems, "converter_implemented": False,
            "missing": "reviewed row/identity/validity adapters and complete target-schema writer"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--layout")
    args = parser.parse_args(argv)
    report = audit_ledger(args.ledger)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with args.report.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k != "problems"}, indent=2))
    if not args.audit_only:
        parser.error("conversion blocked: " + report["missing"] + "; no ROOT output created")
    return 0 if report["ledger_complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
