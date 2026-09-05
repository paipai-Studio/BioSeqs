///|
/// pqsfinder demo: G-quadruplex (PQS) detection on DNA sequences, ported
/// from Bioconductor pqsfinder. Shows the default scoring system, strand
/// handling and option tuning.
///
/// Run: moon run examples/pqs_demo
fn main {
  println("=== pqsfinder Demo (G-quadruplex detection) ===\n")

  // 1. The classic pqsfinder vignette sequence
  println("1. Basic search on the vignette sequence")
  let hits = @bio.pqs_finder("GGGAGGGTGGGTGGG") catch { _ => [] }
  println("   GGGAGGGTGGGTGGG -> \{hits.length()} PQS found")
  for p in hits {
    println(
      "   [\{p.start},\{p.end}) strand=\{p.strand} score=\{p.score} nt=\{p.nt} nb=\{p.nb} nm=\{p.nm} loops=\{p.ll1},\{p.ll2},\{p.ll3}",
    )
  }

  // 2. A longer sequence with several quadruplexes
  println("\n2. Search on a longer sequence")
  let seq = "TTGGGCCTGGGCGGGAGGGTCCTGGGTTGGGGTCCCAAGGGT"
  let hits2 = @bio.pqs_finder(seq) catch { _ => [] }
  println("   \{seq}")
  println("   \{hits2.length()} PQS found (non-overlapping, sorted by start):")
  for p in hits2 {
    println(
      "   [\{p.start},\{p.end}) score=\{p.score} runs=\{p.rl1},\{p.rl2},\{p.rl3} loops=\{p.ll1},\{p.ll2},\{p.ll3}",
    )
  }

  // 3. Both strands
  println("\n3. Strand-specific searches")
  // G-rich part forms a PQS on the plus strand; its reverse complement
  // (C-rich part) forms one on the minus strand
  let two = "GGGAGGGTGGGTGGGAAAAAACCCACCCACCCACTCCC"
  let plus = @bio.pqs_finder(two, strand="+") catch { _ => [] }
  let minus = @bio.pqs_finder(two, strand="-") catch { _ => [] }
  println("   \{two}")
  println("   plus strand: \{plus.length()} PQS")
  for p in plus {
    println("   + [\{p.start},\{p.end}) score=\{p.score}")
  }
  println("   minus strand (reverse complement): \{minus.length()} PQS")
  for p in minus {
    println("   - [\{p.start},\{p.end}) score=\{p.score}")
  }

  // 4. Custom scoring: require two tetrads minimum, stricter defects
  println("\n4. Custom options: longer runs, higher min_score")
  let strict = @bio.pqs_finder(
    seq,
    options=@bio.PqsOptions::of(run_min_len=3, min_score=60),
  ) catch { _ => [] }
  println("   run_min_len=3, min_score=60 -> \{strict.length()} PQS")
  for p in strict {
    println("   [\{p.start},\{p.end}) score=\{p.score} runs=\{p.rl1},\{p.rl2},\{p.rl3}")
  }

  println("\nDone.")
}
