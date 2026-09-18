#!/usr/bin/env python3
"""Sync or verify the README.md module / example / test-file / test-case counts.

Usage:
  python3 scripts/check_readme_counts.py          # verify only (exit 1 on drift)
  python3 scripts/check_readme_counts.py --fix    # update README.md in place
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def count_files(root: str, ext: str) -> int:
    n = 0
    for dirpath, _dirs, files in os.walk(os.path.join(REPO, root)):
        for f in files:
            if f.endswith(ext):
                n += 1
    return n


def count_examples() -> int:
    p = os.path.join(REPO, "examples")
    return sum(
        1 for name in os.listdir(p)
        if os.path.isdir(os.path.join(p, name))
    )


def count_test_cases() -> int:
    d = os.path.join(REPO, "test", "moonbit")
    n = 0
    for f in os.listdir(d):
        if not f.endswith(".mbt"):
            continue
        with open(os.path.join(d, f), encoding="utf-8") as fh:
            for line in fh:
                if line.startswith('test "'):
                    n += 1
    return n


def actual_counts() -> dict:
    return {
        "modules": count_files("src", ".mbt"),
        "examples": count_examples(),
        "test_files": count_files("test/moonbit", ".mbt"),
        "test_cases": count_test_cases(),
    }


# README count locations, as (regex, replacement-template) pairs.
PATTERNS = [
    (re.compile(r"源代码（\d+ 个 \.mbt 模块）"), "源代码（{modules} 个 .mbt 模块）"),
    (re.compile(r"示例程序（\d+ 个演示 demo）"), "示例程序（{examples} 个演示 demo）"),
    (re.compile(r"单元测试（\d+ 个测试文件，\d+ 用例）"),
     "单元测试（{test_files} 个测试文件，{test_cases} 用例）"),
    (re.compile(r"\d+ 个测试全部通过"), "{test_cases} 个测试全部通过"),
    (re.compile(r"\d+ 个测试用例"), "{test_cases} 个测试用例"),
]


def main() -> int:
    fix = "--fix" in sys.argv[1:]
    counts = actual_counts()
    path = os.path.join(REPO, "README.md")
    text = open(path, encoding="utf-8").read()
    new_text = text
    for pattern, tmpl in PATTERNS:
        if not pattern.search(new_text):
            print("!! pattern not found in README:", pattern.pattern)
            return 2
        new_text = pattern.sub(tmpl.format(**counts), new_text)

    if new_text == text:
        print("README counts OK:", counts)
        return 0
    if fix:
        open(path, "w", encoding="utf-8").write(new_text)
        print("README counts updated:", counts)
        return 0
    print("README counts DRIFT. Actual:", counts)
    print("Run with --fix to update.")
    return 1


if __name__ == "__main__":
    sys.exit(main())