#!/usr/bin/env bash
# mdanalysis_compare.sh — Run the MoonBit MD-analysis implementation and
# MDAnalysis on the same battery of test cases and report any differences.
#
# Usage:
#   ./test/python/mdanalysis_compare.sh
#
# Exit status:
#   0  all cases match (within tolerance)
#   1  one or more cases differ (or a prerequisite failed)
#   2  MDAnalysis is not installed
#
# Requirements:
#   - moon (MoonBit toolchain) on PATH or set MOON_BIN environment variable
#   - python3 with numpy + MDAnalysis installed (pip install MDAnalysis)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

MOON="${MOON_BIN:-moon}"

if ! python3 -c "import MDAnalysis" 2>/dev/null; then
  echo "MDAnalysis is not installed (pip install MDAnalysis) — skipping."
  exit 2
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
MB_JSON="$TMPDIR/moonbit.json"
PY_JSON="$TMPDIR/mdanalysis.json"

echo "==> Running MoonBit companion program"
"$MOON" run cmd/md_compare/main.mbt > "$MB_JSON"

echo "==> Running MDAnalysis reference script"
python3 test/python/python_mdanalysis_reference.py > "$PY_JSON"

echo "==> Comparing outputs"
python3 - "$MB_JSON" "$PY_JSON" <<'PYEOF'
import json, sys

mb_path, py_path = sys.argv[1], sys.argv[2]
with open(mb_path) as f:
    mb = {c["name"]: c["result"] for c in json.load(f)}
with open(py_path) as f:
    py = {c["name"]: c["result"] for c in json.load(f)}

if set(mb) != set(py):
    print(f"MISMATCH: different case sets: moonbit={sorted(mb)} mdanalysis={sorted(py)}")
    sys.exit(1)

TOL = 1e-4
mismatches = 0
for name in sorted(py):
    a = [float(v) for v in mb[name].split(",")]
    b = [float(v) for v in py[name].split(",")]
    if len(a) != len(b) or any(abs(x - y) > TOL for x, y in zip(a, b)):
        print(f"[{name}] RESULT mismatch: moonbit={mb[name]} mdanalysis={py[name]}")
        mismatches += 1

if mismatches == 0:
    print(f"\nPASS: all {len(py)} cases match MDAnalysis within {TOL}.")
    sys.exit(0)
print(f"\nFAIL: {mismatches} mismatch(es) out of {len(py)} cases.")
sys.exit(1)
PYEOF
