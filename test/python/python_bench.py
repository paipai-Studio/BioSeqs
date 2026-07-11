#!/usr/bin/env python3
"""Performance benchmark for BioPython's Bio.Seq.Seq.

Runs the same operations as cmd/bench/main.mbt so the two implementations
can be compared. Uses time.perf_counter for wall-clock timing.
"""
import time
import warnings

from Bio.Seq import Seq

warnings.filterwarnings("ignore")


def gen_dna(n):
    """Deterministic pseudo-random DNA sequence of length n."""
    bases = "ACGT"
    state = 12345
    buf = []
    for _ in range(n):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        buf.append(bases[state % 4])
    return "".join(buf)


def bench(label, iters, f):
    start = time.perf_counter()
    for _ in range(iters):
        f()
    elapsed = time.perf_counter() - start
    per_iter = elapsed / iters * 1e6  # us
    print(f"{label}: {elapsed*1e6:.1f} us total, {per_iter:.1f} us/iter ({iters} iters)")
    return elapsed


def main():
    dna = gen_dna(100_000)
    rna = dna.replace("T", "U")
    s = Seq(dna)
    sr = Seq(rna)
    iters_small = 20
    iters_find = 200

    print("=== BioPython Seq benchmark ===")
    print(f"sequence length: {len(dna)} chars")
    print(f"iterations (ops): {iters_small}")
    print("")

    # complement
    bench("complement", iters_small, lambda: s.complement())
    bench("complement_rna", iters_small, lambda: s.complement_rna())
    bench("reverse_complement", iters_small, lambda: s.reverse_complement())
    bench("reverse", iters_small, lambda: s[::-1])

    # transcribe
    bench("transcribe", iters_small, lambda: s.transcribe())
    bench("back_transcribe", iters_small, lambda: sr.back_transcribe())

    # translate
    bench("translate", iters_small, lambda: s.translate())

    # upper / lower
    bench("upper", iters_small, lambda: s.upper())
    bench("lower", iters_small, lambda: s.lower())

    # count / find
    bench("count AAA", iters_find, lambda: s.count("AAA"))
    bench("count_overlap AAA", iters_find, lambda: s.count_overlap("AAA"))
    bench("find ACGT", iters_find, lambda: s.find("ACGT"))
    bench("rfind ACGT", iters_find, lambda: s.rfind("ACGT"))

    # startswith / endswith / contains
    bench("startswith ATG", iters_find, lambda: s.startswith("ATG"))
    bench("endswith TAA", iters_find, lambda: s.endswith("TAA"))
    bench("contains GGG", iters_find, lambda: "GGG" in s)

    # slice
    bench("slice 0..1000", iters_find, lambda: s[0:1000])

    # replace
    bench("replace A->N", iters_small, lambda: s.replace("A", "N"))

    print("")
    print("done.")


if __name__ == "__main__":
    main()
