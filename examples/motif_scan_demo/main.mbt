///|
/// FIMO-like motif scanning demo.
///
/// Demonstrates building a PSSM from a Position Weight Matrix, scanning
/// biological sequences for motif occurrences, and reporting matches with
/// log-likelihood scores and p-values.
fn main {
  println("=== FIMO Motif Scanning Demo ===\n")

  // 1. Create a sample PWM (6-position motif, consensus ATGCAA).
  println("1. Create Position Weight Matrix:")
  let pwm = @bio.motif_scan_sample_pwm()
  println("   Motif ID: \{pwm.motif_id}")
  println("   Motif name: \{pwm.motif_name}")
  println("   Motif width: \{pwm.matrix[0].length()}")
  println("   Strand: \{pwm.strand}")
  println("   Background: [\{pwm.background[0]}, \{pwm.background[1]}, \{pwm.background[2]}, \{pwm.background[3]}]")

  // 2. Build a PSSM (log-likelihood scores) from the PWM.
  println("\n2. Build PSSM (log-likelihood scores):")
  let pssm = @bio.build_pssm(pwm)
  println("   Min possible score: \{pssm.min_score}")
  println("   Max possible score: \{pssm.max_score}")

  // 3. Score the consensus sequence and a non-matching sequence.
  println("\n3. Score sequences against the PSSM:")
  let perfect_score = @bio.compute_score(pssm, "ATGCAA")
  println("   Score of consensus 'ATGCAA': \{perfect_score} (max score)")
  let weak_score = @bio.compute_score(pssm, "CCCCCC")
  println("   Score of 'CCCCCC': \{weak_score}")

  // 4. Compute p-values for those scores.
  println("\n4. Compute p-values:")
  let pval_perfect = @bio.compute_pvalue(pssm, perfect_score)
  println("   p-value for perfect match: \{pval_perfect}")
  let pval_weak = @bio.compute_pvalue(pssm, weak_score)
  println("   p-value for weak match: \{pval_weak}")

  // 5. Scan sample sequences for motif occurrences (p-value threshold 0.001).
  println("\n5. Scan sample sequences (p-value threshold = 0.001):")
  let sequences = @bio.motif_scan_sample_sequences()
  for i in 0..<sequences.length() {
    let pair = sequences[i]
    println("   \{pair.0}: \{pair.1} (length \{pair.1.length()})")
  }
  let matches = @bio.scan_sequences(pssm, sequences, 0.001)
  println("\n   Found \{matches.length()} match(es):")
  for i in 0..<matches.length() {
    let m = matches[i]
    println(
      "   - \{m.sequence_id} | \{m.strand} | pos \{m.start}-\{m.end} | score=\{m.score} | p=\{m.p_value} | \{m.matched_sequence}",
    )
  }

  // 6. Demonstrate building a PWM from a count matrix.
  println("\n6. Build PWM from a count matrix (with pseudocount 0.5):")
  // Counts for a 4-position motif. Rows: A, C, G, T. Columns: positions.
  let counts = [
    [18, 1, 1, 2],
    [2, 17, 0, 1],
    [1, 1, 18, 1],
    [1, 3, 3, 18],
  ]
  let pwm2 = @bio.motif_scan_pwm_from_counts("M0002", "FromCounts", counts, 0.5)
  println("   Motif ID: \{pwm2.motif_id}")
  println("   Motif name: \{pwm2.motif_name}")
  println("   Width: \{pwm2.matrix[0].length()}")
  let pssm2 = @bio.build_pssm(pwm2)
  println("   PSSM min score: \{pssm2.min_score}")
  println("   PSSM max score: \{pssm2.max_score}")

  // 7. Demonstrate reverse complement.
  println("\n7. Reverse complement:")
  let seq = "ATGCAA"
  let rc = @bio.reverse_complement_sequence(seq)
  println("   \{seq} -> \{rc}")

  // 8. Print a full summary of the scan results.
  println("\n8. Full summary:")
  let summary = @bio.motif_scan_summary(matches)
  println(summary)

  // 9. Filter matches by a stricter p-value threshold.
  println("9. Filter matches with stricter threshold (p <= 0.0001):")
  let strict = @bio.filter_matches(matches, 0.0001)
  println("   \{strict.length()} match(es) survive stricter filtering")

  println("\n=== Demo Complete ===")
}
