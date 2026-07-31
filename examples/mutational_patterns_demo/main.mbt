///|
/// MutationalPatterns module example - Somatic mutation pattern analysis.
///
/// Run: moon run examples/mutational_patterns_demo
fn main {
  println("=== MutationalPatterns Module Demo ===")
  println("Somatic Mutation Pattern Analysis (96-channel)")
  println("")

  // ---- Section 1: Mutation categories ----
  println("=== 1. Mutation Categories (Pyrimidine-based) ===")
  let types = @bio.all_mutation_types()
  println("Total categories: " + types.length().to_string())
  let mut ti = 0
  while ti < types.length() {
    println("  " + types[ti].to_string())
    ti = ti + 1
  }
  println("")

  // ---- Section 2: Trinucleotide contexts ----
  println("=== 2. Trinucleotide Contexts ===")
  let contexts = @bio.trinucleotide_contexts()
  println("Total contexts: " + contexts.length().to_string())
  let mut ci = 0
  while ci < contexts.length() {
    println("  " + contexts[ci])
    ci = ci + 1
  }
  println("")

  // ---- Section 3: Channel labels (96 channels) ----
  println("=== 3. 96-Channel Labels ===")
  let labels = @bio.get_channel_labels()
  println("Total channels: " + labels.length().to_string())
  // Show first 6 and last 6 labels
  println("First 6: " + labels[0] + ", " + labels[1] + ", " + labels[2] + ", " + labels[3] + ", " + labels[4] + ", " + labels[5])
  println("Last 6:  " + labels[90] + ", " + labels[91] + ", " + labels[92] + ", " + labels[93] + ", " + labels[94] + ", " + labels[95])
  println("")

  // ---- Section 4: Mutation type detection ----
  println("=== 4. Mutation Type Detection ===")
  // C>T transition
  let mtype_ct = @bio.get_mutation_type("C", "T")
  println("  C>T -> " + mtype_ct.to_string())
  // G>A should be reverse-complemented to C>T (pyrimidine-based)
  let mtype_ga = @bio.get_mutation_type("G", "A")
  println("  G>A -> " + mtype_ga.to_string() + " (pyrimidine-equivalent)")
  // C>A transversion
  let mtype_ca = @bio.get_mutation_type("C", "A")
  println("  C>A -> " + mtype_ca.to_string())
  println("")

  // ---- Section 5: Context normalization ----
  println("=== 5. Trinucleotide Context Normalization ===")
  // ACT (middle is C, pyrimidine) -> stays as ACT
  let norm1 = @bio.normalize_context("ACT")
  println("  ACT -> " + norm1 + " (pyrimidine middle, unchanged)")
  // AGT (middle is G, purine) -> reverse complement -> ACT
  let norm2 = @bio.normalize_context("AGT")
  println("  AGT -> " + norm2 + " (purine middle, reverse-complemented)")
  println("")

  // ---- Section 6: Build mutation matrix ----
  println("=== 6. Build 96-channel Mutation Matrix ===")
  let mutations = @bio.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @bio.build_mutation_matrix(mutations, sample_ids)
  println("Samples: " + matrix.sample_ids.length().to_string())
  println("Channels: " + matrix.channel_labels.length().to_string())
  let mut si = 0
  while si < matrix.sample_ids.length() {
    println("  " + matrix.sample_ids[si] + ": " + matrix.total_mutations[si].to_string() + " total mutations")
    si = si + 1
  }
  println("")

  // ---- Section 7: Mutation spectrum summary ----
  println("=== 7. Mutation Spectrum Summary ===")
  let summary = matrix.spectrum_summary()
  println(summary)
  println("")

  // ---- Section 8: Sample-sample cosine similarity ----
  println("=== 8. Cosine Similarity Between Samples ===")
  let sim_same = matrix.sample_cosine_similarity(0, 0)
  let sim_diff = matrix.sample_cosine_similarity(0, 1)
  println("  S1 vs S1 (identical): " + sim_same.to_string())
  println("  S1 vs S2 (different): " + sim_diff.to_string())
  println("")

  // ---- Section 9: Relative frequencies ----
  println("=== 9. Relative Mutation Frequencies (per sample) ===")
  let freqs = matrix.to_relative_frequencies()
  let mut sum0 = 0.0
  let mut fi = 0
  while fi < 96 {
    sum0 = sum0 + freqs[fi][0]
    fi = fi + 1
  }
  println("  Sum of S1 frequencies (should be ~1.0): " + sum0.to_string())
  println("")

  // ---- Section 10: COSMIC signatures and signature fitting ----
  println("=== 10. Known COSMIC Signatures ===")
  let sigs = @bio.known_cosmic_signatures()
  let mut sj = 0
  while sj < sigs.length() {
    let sig = sigs[sj]
    let mut total = 0.0
    let mut pi = 0
    while pi < sig.probabilities.length() {
      total = total + sig.probabilities[pi]
      pi = pi + 1
    }
    println("  " + sig.id + ": " + sig.probabilities.length().to_string() + " channels, total prob = " + total.to_string())
    sj = sj + 1
  }
  println("")

  // ---- Section 11: Fit signatures to a sample ----
  println("=== 11. Fit Signatures to Sample S1 ===")
  // Use S1 frequencies as sample profile
  let s1_profile : Array[Double] = Array::new()
  let mut pp = 0
  while pp < 96 {
    s1_profile.push(freqs[pp][0])
    pp = pp + 1
  }
  let contributions = @bio.fit_signatures(s1_profile, sigs)
  println("Signature contributions:")
  let mut ck = 0
  while ck < contributions.length() {
    println("  " + sigs[ck].id + ": " + contributions[ck].to_string())
    ck = ck + 1
  }
  let mut total_contrib = 0.0
  let mut tk = 0
  while tk < contributions.length() {
    total_contrib = total_contrib + contributions[tk]
    tk = tk + 1
  }
  println("  Total (should be ~1.0): " + total_contrib.to_string())
  println("")

  // ---- Section 12: Channel index lookup ----
  println("=== 12. Channel Index Lookup ===")
  let idx1 = @bio.get_channel_index("C", "T", "ACT")
  let idx2 = @bio.get_channel_index("C", "A", "CCA")
  println("  C>T at ACT context -> channel " + idx1.to_string())
  println("  C>A at CCA context -> channel " + idx2.to_string())
  println("")

  println("=== Demo complete ===")
}
