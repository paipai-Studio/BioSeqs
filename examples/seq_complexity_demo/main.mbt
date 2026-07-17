///|
/// Sequence Complexity Demo - Sequence composition and complexity analysis.
///
/// Demonstrates Shannon entropy, linguistic complexity, GC/AT content,
/// GC/AT skew, DUST score, sequence signature, nucleotide frequencies,
/// Chaos Game Representation, and sequence similarity.
///
/// Run: moon run examples/seq_complexity_demo/main.mbt
fn main {
  println("=== Sequence Complexity Demo ===\n")

  // Test sequences
  let seq_random = "ATCGATCGATCGATCGATCG"
  let seq_repeat = "ATATATATATATATATATAT"
  let seq_gc = "GCGCGCGCGCGCGCGCGCGC"
  let seq_low_complexity = "AAAAAAAAAAAAAAAAAAAA"
  println("Test sequences:")
  println("  Random-like: " + seq_random)
  println("  Repetitive : " + seq_repeat)
  println("  GC-rich    : " + seq_gc)
  println("  Low complex: " + seq_low_complexity)

  // === 1. Shannon entropy ===
  println("\n1. Shannon Entropy (bits, higher = more complex):")
  println("==================================================")
  let entropy_random = @src.shannon_entropy_bits(seq_random)
  let entropy_repeat = @src.shannon_entropy_bits(seq_repeat)
  let entropy_gc = @src.shannon_entropy_bits(seq_gc)
  let entropy_low = @src.shannon_entropy_bits(seq_low_complexity)
  println("  Random-like: " + entropy_random.to_string())
  println("  Repetitive : " + entropy_repeat.to_string())
  println("  GC-rich    : " + entropy_gc.to_string())
  println("  Low complex: " + entropy_low.to_string())

  // === 2. Linguistic complexity ===
  println("\n2. Linguistic Complexity (k=2, higher = more complex):")
  println("=======================================================")
  let lc_random = @src.linguistic_complexity(seq_random, k=2)
  let lc_repeat = @src.linguistic_complexity(seq_repeat, k=2)
  let lc_gc = @src.linguistic_complexity(seq_gc, k=2)
  let lc_low = @src.linguistic_complexity(seq_low_complexity, k=2)
  println("  Random-like: " + lc_random.to_string())
  println("  Repetitive : " + lc_repeat.to_string())
  println("  GC-rich    : " + lc_gc.to_string())
  println("  Low complex: " + lc_low.to_string())

  // === 3. GC/AT content ===
  println("\n3. GC/AT Content (percentage):")
  println("===============================")
  let sequences = [seq_random, seq_repeat, seq_gc, seq_low_complexity]
  let labels = ["Random-like", "Repetitive ", "GC-rich    ", "Low complex"]
  let mut si = 0
  while si < sequences.length() {
    let gc = @src.gc_content_percent(sequences[si])
    let at = @src.at_content(sequences[si])
    println(
      "  " +
      labels[si] +
      ": GC=" +
      gc.to_string() +
      "%, AT=" +
      at.to_string() +
      "%",
    )
    si = si + 1
  }

  // === 4. GC/AT skew ===
  println("\n4. GC/AT Skew:")
  println("===============")
  let seq_g_skew = "GGGGGGGGCCCCCCCC"
  let seq_c_skew = "CCCCCCCCGGGGGGGG"
  let seq_a_skew = "AAAAAAAATTTTTTTT"
  let seq_t_skew = "TTTTTTTTAAAAAAAA"
  println(
    "  " + seq_g_skew + " GC skew: " + @src.gc_skew(seq_g_skew).to_string(),
  )
  println(
    "  " + seq_c_skew + " GC skew: " + @src.gc_skew(seq_c_skew).to_string(),
  )
  println(
    "  " + seq_a_skew + " AT skew: " + @src.at_skew(seq_a_skew).to_string(),
  )
  println(
    "  " + seq_t_skew + " AT skew: " + @src.at_skew(seq_t_skew).to_string(),
  )

  // === 5. DUST score ===
  println("\n5. DUST Score (higher = lower complexity):")
  println("===========================================")
  let dust_random = @src.dust_score(seq_random)
  let dust_repeat = @src.dust_score(seq_repeat)
  let dust_low = @src.dust_score(seq_low_complexity)
  println("  Random-like: " + dust_random.to_string())
  println("  Repetitive : " + dust_repeat.to_string())
  println("  Low complex: " + dust_low.to_string())

  // === 6. Sequence signature ===
  println("\n6. Sequence Signature (k-mer frequency, k=2):")
  println("==============================================")
  let sig = @src.sequence_signature("ATCGATCG", k=2)
  println("  Unique 2-mers: " + sig.size().to_string())
  let sig_kmers = ["AT", "TC", "CG", "GA"]
  let mut ski = 0
  while ski < sig_kmers.length() {
    let kmer = sig_kmers[ski]
    match sig.get(kmer) {
      Some(freq) => println("  " + kmer + ": " + freq.to_string())
      None => println("  " + kmer + ": not found")
    }
    ski = ski + 1
  }

  // === 7. Nucleotide frequencies ===
  println("\n7. Nucleotide Frequencies:")
  println("==========================")
  let (a, t, g, c) = @src.nucleotide_frequencies(seq_random)
  println("  Sequence: " + seq_random)
  println("  A: " + a.to_string())
  println("  T: " + t.to_string())
  println("  G: " + g.to_string())
  println("  C: " + c.to_string())
  let (a2, t2, g2, c2) = @src.nucleotide_frequencies(seq_gc)
  println("  Sequence: " + seq_gc)
  println(
    "  A: " +
    a2.to_string() +
    ", T: " +
    t2.to_string() +
    ", G: " +
    g2.to_string() +
    ", C: " +
    c2.to_string(),
  )

  // === 8. Chaos Game Representation ===
  println("\n8. Chaos Game Representation (CGR coordinates):")
  println("===============================================")
  let cgr_seqs = ["AAAA", "TTTT", "GGGG", "CCCC", "ATCG"]
  let mut ci = 0
  while ci < cgr_seqs.length() {
    let (x, y) = @src.chaos_game_representation(cgr_seqs[ci])
    println(
      "  " +
      cgr_seqs[ci] +
      " -> (x=" +
      x.to_string() +
      ", y=" +
      y.to_string() +
      ")",
    )
    ci = ci + 1
  }

  // === 9. Sequence similarity ===
  println("\n9. Sequence Similarity (cosine similarity of k-mer vectors, k=2):")
  println("=================================================================")
  let sim_identical = @src.sequence_similarity(seq_random, seq_random, k=2)
  let sim_similar = @src.sequence_similarity("ATCGATCG", "ATCGATCGATCG", k=2)
  let sim_different = @src.sequence_similarity("ATCGATCG", "GCTAGCTA", k=2)
  println(
    "  'ATCGATCG' vs 'ATCGATCG' (identical): " + sim_identical.to_string(),
  )
  println(
    "  'ATCGATCG' vs 'ATCGATCGATCG' (similar): " + sim_similar.to_string(),
  )
  println(
    "  'ATCGATCG' vs 'GCTAGCTA' (different): " + sim_different.to_string(),
  )

  println("\n=== Sequence Complexity Demo Complete ===")
}
