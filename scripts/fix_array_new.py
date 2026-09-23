#!/usr/bin/env python3
"""Migrate deprecated `Array::new(...)` to the `Array(capacity=...)` constructor.

`Array::new` was deprecated in favour of the builtin `Array(capacity=...)`
constructor.  This script rewrites:
  Array::new()            -> Array(capacity=0)
  Array::new(capacity=X)  -> Array(capacity=X)   (single- or multi-line)

A word-boundary guard avoids touching substring matches such as
`BitArray::new`, `DelayedArray::new`, `SuffixArray::new`, `LCPArray::new`.

Usage:
  python3 scripts/fix_array_new.py --dry-run   # preview only
  python3 scripts/fix_array_new.py             # write the changes
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOTS = ["src", "test/moonbit", "examples", "cmd"]

EMPTY = re.compile(r"(?<![A-Za-z0-9_])Array::new\(\)")
CAP = re.compile(r"(?<![A-Za-z0-9_])Array::new\(\s*capacity=")


def walk_mbt(root: str):
    base = os.path.join(REPO, root)
    if not os.path.isdir(base):
        return []
    out = []
    for dirpath, _dirs, files in os.walk(base):
        for f in files:
            if f.endswith(".mbt"):
                out.append(os.path.join(dirpath, f))
    return out


def process(path: str, dry_run: bool) -> int:
    with open(path, encoding="utf-8") as fh:
        s = fh.read()
    new = EMPTY.sub("Array(capacity=0)", s)
    new = CAP.sub("Array(capacity=", new)
    n = s.count("Array::new(") - new.count("Array::new(")  # rough delta
    if new != s and not dry_run:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(new)
    return 1 if new != s else 0


def main() -> int:
    dry_run = "--dry-run" in sys.argv
    changed = 0
    total = 0
    for root in ROOTS:
        for p in walk_mbt(root):
            n = 0
            with open(p, encoding="utf-8") as fh:
                s = fh.read()
            new = EMPTY.sub("Array(capacity=0)", s)
            new = CAP.sub("Array(capacity=", new)
            if new != s:
                total += 1
                if not dry_run:
                    with open(p, "w", encoding="utf-8") as fh:
                        fh.write(new)
    print(f"{'[DRY RUN] ' if dry_run else ''}{total} 个文件涉及 Array::new 迁移")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())