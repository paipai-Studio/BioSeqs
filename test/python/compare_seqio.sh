#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

MOON="${MOON_BIN:-moon}"

echo "=== Running MoonBit SeqIO tests ==="
"$MOON" run cmd/seqio_main/main.mbt 2>/dev/null > /tmp/moonbit_seqio.txt

echo "=== Running BioPython reference ==="
python3 test/python/python_seqio_reference.py > /tmp/biopython_seqio.txt 2>&1

echo "=== Comparing outputs ==="
if diff -u /tmp/moonbit_seqio.txt /tmp/biopython_seqio.txt; then
    echo ""
    echo "ALL TESTS PASSED: MoonBit SeqIO matches BioPython exactly"
    exit 0
else
    echo ""
    echo "TESTS FAILED: outputs differ"
    exit 1
fi
