#!/usr/bin/env bash
# compare.sh — Run the MoonBit reimplementation and BioPython's Bio.Seq on the
# same battery of test cases and report any differences.
#
# Usage:
#   ./tests/compare.sh
#
# Exit status:
#   0  all cases match
#   1  one or more cases differ (or a prerequisite failed)
#
# Requirements:
#   - moon (MoonBit toolchain) on PATH
#   - python3 with biopython installed (pip install biopython)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
MB_JSON="$TMPDIR/moonbit.json"
PY_JSON="$TMPDIR/biopython.json"

echo "==> Building MoonBit project"
moon build

echo "==> Running MoonBit companion program"
moon run cmd/main > "$MB_JSON"

echo "==> Running BioPython reference script"
python3 tests/python_reference.py > "$PY_JSON"

echo "==> Comparing outputs"
python3 - "$MB_JSON" "$PY_JSON" <<'PYEOF'
import json, sys

mb_path, py_path = sys.argv[1], sys.argv[2]
with open(mb_path) as f:
    mb = json.load(f)
with open(py_path) as f:
    py = json.load(f)

print(f"MoonBit cases : {len(mb)}")
print(f"BioPython cases: {len(py)}")

if len(mb) != len(py):
    print(f"MISMATCH: different number of cases ({len(mb)} vs {len(py)})")
    sys.exit(1)

mismatches = 0
for i, (m, p) in enumerate(zip(mb, py)):
    if m["name"] != p["name"]:
        print(f"[{i}] NAME mismatch: moonbit={m['name']!r} biopython={p['name']!r}")
        mismatches += 1
    elif m["result"] != p["result"]:
        print(f"[{i}] RESULT mismatch ({m['name']}): moonbit={m['result']!r} biopython={p['result']!r}")
        mismatches += 1

if mismatches == 0:
    print(f"\nPASS: all {len(mb)} cases match BioPython exactly.")
    sys.exit(0)
else:
    print(f"\nFAIL: {mismatches} mismatch(es) out of {len(mb)} cases.")
    sys.exit(1)
PYEOF
