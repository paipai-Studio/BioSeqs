///|
/// motif_matrix demo - Unified motif matrix abstraction (PFM/PPM/PWM/PSSM).
/// Inspired by Biopython Bio.motifs.matrix.

fn main {
  println("=== motif_matrix Demo: PFM/PPM/PWM/PSSM ===")
  println("")

  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGTT", "ACGTT", "ACGTA", "TCGTT", "ACGCT"]

  // Step 1: PFM
  println("--- Step 1: Position Frequency Matrix ---")
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  println(@src.motif_pfm_to_string(pfm))
  println("")

  // Step 2: PPM
  println("--- Step 2: Position Probability Matrix ---")
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=1.0)
  println("Consensus: " + @src.motif_ppm_consensus(ppm))
  let ic = @src.motif_ppm_ic_per_position(ppm)
  println("IC per position (bits): " + ic.to_string())
  println("Total IC: " + @src.motif_ppm_total_ic(ppm).to_string() + " bits")
  println("")

  // Step 3: PWM
  println("--- Step 3: Position Weight Matrix ---")
  let pwm = @src.motif_pfm_to_pwm(pfm, pseudocount=1.0)
  println(@src.motif_pwm_to_string(pwm))
  let score = @src.motif_pwm_score(pwm, "ACGTT")
  println("Score for ACGTT: " + score.to_string())
  println("")

  // Step 4: PSSM
  println("--- Step 4: Position Specific Scoring Matrix ---")
  let pssm = @src.motif_pwm_to_pssm(pwm, scale=100.0)
  let iscore = @src.motif_pssm_score(pssm, "ACGTT")
  println("Integer score for ACGTT: " + iscore.to_string())
  println("")

  // Step 5: Scan
  println("--- Step 5: Scanning ---")
  let hits = @src.motif_pwm_scan(pwm, "TTACGTTTACGTTT", 0.0)
  println("Hits in TTACGTTTACGTTT: " + hits.length().to_string())
  let mut i = 0
  while i < hits.length() {
    println("  pos=" + hits[i].0.to_string() + " score=" + hits[i].1.to_string())
    i = i + 1
  }
  println("")

  // Step 6: Reverse complement
  println("--- Step 6: Reverse Complement PWM ---")
  let rc = @src.motif_pwm_reverse_complement(pwm)
  println("RC PWM rows: " + rc.values.length().to_string())
  println("")

  println("=== motif_matrix Demo Complete ===")
}