///|
/// dada2 demo - Amplicon Sequence Variant inference pipeline.
/// Demonstrates: dereplicate → learnErrors → dada → mergePairs → removeChimeras → makeTable → assignTaxonomy.

fn main {
  println("=== dada2 Demo: Amplicon Sequence Variant Inference ===")
  println("")

  // Step 1: Dereplication
  println("--- Step 1: Dereplication ---")
  let seqs = [
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGA",
    "ACGTACGTACGTACGA",
    "ACGTACGTACGTACGT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAC",
    "ACGTACGTACGTACGT",
  ]
  let derep = @src.dereplicate(seqs, [])
  println("Input reads: " + derep.total_reads.to_string())
  println("Unique sequences: " + derep.uniques.length().to_string())
  let mut i = 0
  while i < derep.uniques.length() {
    let e = derep.uniques[i]
    println("  " + e.sequence + " (abundance=" + e.abundance.to_string() + ")")
    i = i + 1
  }
  println("")

  // Step 2: Learn Error Model
  println("--- Step 2: Learn Error Model ---")
  let model = @src.learn_errors(derep, 40)
  println("Max quality: " + model.max_qual.to_string())
  println("Reads processed: " + model.n_reads.to_string())
  println("Error rate (A->C at Q20): " + model.err_mat[@src.transition_index(0, 1)][20].to_string())
  println("Match rate (A->A at Q20): " + model.err_mat[@src.transition_index(0, 0)][20].to_string())
  println("")

  // Step 3: DADA Inference
  println("--- Step 3: DADA Inference ---")
  let dada = @src.dada_inference(derep, model)
  println("Input reads: " + dada.n_input_reads.to_string())
  println("ASV clusters: " + dada.n_clusters.to_string())
  let mut ci = 0
  while ci < dada.partitions.length() {
    let p = dada.partitions[ci]
    println("  ASV " + (ci + 1).to_string() + ": " + p.center + " (abundance=" + p.abundance.to_string() + ", members=" + p.members.length().to_string() + ")")
    ci = ci + 1
  }
  println("")

  // Step 4: Merge Paired-End Reads
  println("--- Step 4: Merge Paired-End Reads ---")
  let fwd_reads = ["ACGTACGTACGTACGTAAAA", "ACGTACGTACGTACGTAAAA"]
  let rev_reads = ["TTTTACGTACGTACGTACGT", "TTTTACGTACGTACGTACGT"]
  let merged = @src.merge_pairs(fwd_reads, rev_reads, min_overlap=12)
  let mut mi = 0
  while mi < merged.length() {
    let m = merged[mi]
    println("  Read " + (mi + 1).to_string() + ": " + (if m.accepted { "accepted" } else { "rejected" }))
    if m.accepted {
      println("    Merged: " + m.sequence + " (overlap=" + m.overlap.to_string() + ", mismatch=" + m.mismatch.to_string() + ")")
    }
    mi = mi + 1
  }
  println("")

  // Step 5: Remove Chimeras
  println("--- Step 5: Remove Chimeras ---")
  let parent_a = "AAAAAAATTTTTTT"
  let parent_b = "GGGGGGGCCCCCCC"
  let chimera = "AAAAAAACCCCCCC"
  let asvs = [parent_a, parent_b, chimera]
  let abundances = [100, 100, 10]
  let (clean_asvs, clean_abund, _flags) = @src.remove_chimeras(asvs, abundances)
  println("Input ASVs: " + asvs.length().to_string())
  println("Chimeras removed: " + (asvs.length() - clean_asvs.length()).to_string())
  println("Clean ASVs: " + clean_asvs.length().to_string())
  let mut ai = 0
  while ai < clean_asvs.length() {
    println("  " + clean_asvs[ai] + " (abundance=" + clean_abund[ai].to_string() + ")")
    ai = ai + 1
  }
  println("")

  // Step 6: Make Sequence Table
  println("--- Step 6: Make Sequence Table ---")
  let derep2 = @src.create_example_derep()
  let model2 = @src.learn_errors(derep2, 40)
  let dada2 = @src.dada_inference(derep2, model2)
  let table = @src.make_sequence_table([("sample1", dada2), ("sample2", dada)])
  println("ASVs: " + table.n_asvs.to_string())
  println("Samples: " + table.n_samples.to_string())
  let sample_totals = @src.seq_table_sample_totals(table)
  let mut si = 0
  while si < table.n_samples {
    println("  " + table.sample_names[si] + " total: " + sample_totals[si].to_string())
    si = si + 1
  }
  println("")

  // Step 7: Assign Taxonomy
  println("--- Step 7: Assign Taxonomy ---")
  let ref_db = @src.create_example_ref_db()
  let tax_results = @src.assign_taxonomy(clean_asvs, ref_db)
  let mut ti = 0
  while ti < tax_results.length() {
    let r = tax_results[ti]
    println("  ASV: " + r.asv)
    println("    Taxonomy: " + r.taxonomy.join(" > "))
    println("    Confidence: " + r.confidence.to_string())
    ti = ti + 1
  }
  println("")

  // Step 8: Quality Profile
  println("--- Step 8: Quality Profile ---")
  let quals = [
    [30.0, 31.0, 32.0, 33.0, 34.0],
    [28.0, 29.0, 30.0, 31.0, 32.0],
    [32.0, 33.0, 34.0, 35.0, 36.0],
  ]
  let profile = @src.quality_profile(quals)
  println("Positions: " + profile.positions.to_string())
  let mut qi = 0
  while qi < profile.positions {
    println("  Pos " + qi.to_string() + ": mean=" + profile.mean[qi].to_string() + ", min=" + profile.min[qi].to_string() + ", max=" + profile.max[qi].to_string())
    qi = qi + 1
  }
  println("")

  println("=== dada2 Demo Complete ===")
}
