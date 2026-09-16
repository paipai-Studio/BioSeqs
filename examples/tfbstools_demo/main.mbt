///|
/// Demo for tfbstools module: TF binding site analysis.

fn main {
  println("=== TFBSTools Demo: TF Binding Site Analysis ===")
  println("")

  // Example 1: PFM normalization
  println("--- Step 1: PFM Normalization ---")
  let counts = [
    [8.0, 2.0, 0.0, 0.0],
    [0.0, 10.0, 0.0, 0.0],
    [1.0, 1.0, 8.0, 0.0],
    [0.0, 0.0, 2.0, 8.0],
  ]
  let pfm = @src.tfbs_pfm_normalize(counts, pseudocount=0.8)
  println("PFM: " + pfm.length().to_string() + " positions")
  println("")

  // Example 2: PFM to PWM
  println("--- Step 2: PFM to PWM ---")
  let pwm = @src.tfbs_pfm_to_pwm(pfm)
  println("PWM: " + pwm.length().to_string() + " positions")
  for i in 0..<pwm.length() {
    println("  pos " + (i+1).to_string() + ": A=" + pwm[i][0].to_string() + " C=" + pwm[i][1].to_string() +
            " G=" + pwm[i][2].to_string() + " T=" + pwm[i][3].to_string())
  }
  println("")

  // Example 3: Consensus sequence
  println("--- Step 3: Consensus Sequence ---")
  let cons = @src.tfbs_consensus(pfm)
  println("Consensus: " + cons)
  println("")

  // Example 4: Information content
  println("--- Step 4: Information Content ---")
  let (per_pos, total) = @src.tfbs_info_content(pfm)
  println("Total IC: " + total.to_string() + " bits")
  println("Per-position IC: " + per_pos.to_string())
  println("")

  // Example 5: Sequence scanning
  println("--- Step 5: Sequence Scanning ---")
  let seq = "ACGTACGTACGTACGT"
  let hits = @src.tfbs_scan_sequence(seq, pwm, threshold=0.0)
  println("Sequence: " + seq)
  println("Hits found: " + hits.length().to_string())
  for i in 0..<hits.length() {
    let h = hits[i]
    println("  " + h.sequence_id + ":" + h.start.to_string() + "-" + h.end.to_string() +
            " " + h.strand + " score=" + h.score.to_string() + " " + h.sequence)
  }
  println("")

  // Example 6: Batch scanning
  println("--- Step 6: Batch Scanning ---")
  let seqs = [("seq1", "ACGTACGT"), ("seq2", "TGCAACGT")]
  let batch_hits = @src.tfbs_scan_batch(seqs, pwm, threshold=0.0)
  println("Batch hits: " + batch_hits.length().to_string())
  println("")

  // Example 7: Motif enrichment
  println("--- Step 7: Motif Enrichment ---")
  let r = @src.tfbs_enrich_hypergeometric(
    "MA001", "TF1",
    20, 100, 5, 100,
  )
  println("TF: " + r.tf_name + " (" + r.matrix_id + ")")
  println("  Targets: " + r.n_targets_hit.to_string() + "/" + r.total_targets.to_string())
  println("  Background: " + r.n_background_hit.to_string() + "/" + r.total_background.to_string())
  println("  Fold enrichment: " + r.fold_enrichment.to_string())
  println("  p-value: " + r.p_value.to_string())
  println("")

  // Example 8: Multiple enrichment with BH correction
  println("--- Step 8: BH Correction ---")
  let r1 = @src.tfbs_enrich_hypergeometric("M1", "TF1", 10, 100, 5, 100)
  let r2 = @src.tfbs_enrich_hypergeometric("M2", "TF2", 20, 100, 10, 100)
  let r3 = @src.tfbs_enrich_hypergeometric("M3", "TF3", 5, 100, 15, 100)
  let results = @src.tfbs_adjust_pvalues([r1, r2, r3])
  println("Adjusted results:")
  for i in 0..<results.length() {
    let res = results[i]
    println("  " + res.tf_name + ": p=" + res.p_value.to_string() + " adj_p=" + res.adj_p.to_string())
  }
  println("")

  // Example 9: Example JASPAR matrices
  println("--- Step 9: Example JASPAR Matrices ---")
  let mats = @src.tfbs_example_jaspar()
  println("Example matrices: " + mats.length().to_string())
  println("")

  println("=== TFBSTools Demo Complete ===")
}
