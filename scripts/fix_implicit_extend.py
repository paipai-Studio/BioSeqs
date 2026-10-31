#!/usr/bin/env python3
"""Explicitly promote derived trait methods as regular methods.

Each `} derive(Eq, Debug, ...)` generates an impl whose methods are implicitly
promoted as dot-methods.  MoonBit deprecates this (`implicit_impl_as_method`)
and asks for an explicit `pub extend X with Trait::{methods}` (or `extend ...`
for private types).  This script inserts those declarations right after each
`} derive(...)` line, keeping the current dot-method behaviour.

Usage:
  python3 scripts/fix_implicit_extend.py --dry-run   # preview only
  python3 scripts/fix_implicit_extend.py             # write the changes
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TRAITS = {
    "Eq": "Eq::{equal, not_equal}",
    "Debug": "@moonbitlang/core/debug.Debug::{to_repr}",
    "Hash": "Hash::{hash, hash_combine}",
    "Default": "Default::{default}",
}

TYPE_RE = re.compile(r"^(pub(?:\([^)]*\))? )?(struct|enum)\s+([A-Za-z0-9_]+)")
DERIVE_RE = re.compile(r"} derive\(([^)]*)\)")


def process(path: str, dry_run: bool) -> int:
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().split("\n")

    out = []
    changes = 0
    current_name = None
    current_pub = False
    for line in lines:
        out.append(line)
        m = TYPE_RE.match(line)
        if m:
            current_pub = m.group(1) is not None
            current_name = m.group(3)
            continue
        dm = DERIVE_RE.search(line)
        if dm and current_name:
            traits = [t.strip() for t in dm.group(1).split(",") if t.strip() in TRAITS]
            for t in traits:
                vis = "pub " if current_pub else ""
                out.append(f"{vis}extend {current_name} with {TRAITS[t]}")
                changes += 1

    if changes and not dry_run:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(out))
    return changes


def main() -> int:
    dry_run = "--dry-run" in sys.argv
    total = 0
    changed_files = 0
    for name in sorted(os.listdir(os.path.join(REPO, "src"))):
        if not name.endswith(".mbt"):
            continue
        p = os.path.join(REPO, "src", name)
        n = process(p, dry_run)
        if n:
            changed_files += 1
            total += n
    print(f"{'[DRY RUN] ' if dry_run else ''}{total} extend 声明将插入，涉及 {changed_files} 个文件")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())