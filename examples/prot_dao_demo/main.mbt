///|
/// prot_dao module example - Protein disorder prediction (IUPred).
///
/// Run: moon run examples/prot_dao_demo
fn main {
  println("=== Bio.SeqUtils.ProtDao Module Demo ===")
  println("Protein Disorder Prediction (IUPred-like Algorithm)")
  println("")

  // Amino acid scores
  println("=== Amino Acid Disorder Scores ===")
  let disordered_aa : Array[UInt16] = ['E', 'P', 'K', 'S']
  let ordered_aa : Array[UInt16] = ['I', 'W', 'Y', 'C']
  let mut i = 0
  while i < disordered_aa.length() {
    let aa = disordered_aa[i]
    let score = @bio.prot_dao_disorder_score(aa)
    let energy = @bio.prot_dao_energy_score(aa)
    println("  " + aa.to_string() + ": disorder_score=" + score.to_string() + ", energy_score=" + energy.to_string())
    i = i + 1
  }
  println("  ---")
  let mut j = 0
  while j < ordered_aa.length() {
    let aa = ordered_aa[j]
    let score = @bio.prot_dao_disorder_score(aa)
    let energy = @bio.prot_dao_energy_score(aa)
    println("  " + aa.to_string() + ": disorder_score=" + score.to_string() + ", energy_score=" + energy.to_string())
    j = j + 1
  }
  println("")

  // Predict disorder for a protein
  println("=== Disorder Prediction ===")
  let sequence = @bio.prot_dao_sample_sequence()
  println("Protein sequence:")
  println("  " + sequence)
  println("  Length: " + sequence.length().to_string())
  println("")

  let result = @bio.prot_dao_predict(sequence)
  println(result.summary())
  println("")

  // Disordered residues
  println("=== Disordered Residues ===")
  println("  Total residues: " + result.sequence.length().to_string())
  println("  Disordered: " + result.get_n_disordered().to_string())
  println("  Long disordered: " + result.get_n_disordered_long().to_string())
  println("  Fraction disordered: " + result.get_fraction_disordered().to_string())
  println("")

  // Disordered sequence view
  println("=== Disordered Sequence View ===")
  let ds = result.disordered_sequence()
  println("  Disordered residues (marked with *, ordered with -):")
  println("  " + ds)
  println("")

  // Disordered regions
  println("=== Disordered Regions ===")
  let regions = result.get_regions()
  println("  Number of regions: " + result.get_n_regions().to_string())
  let mut k = 0
  while k < regions.length() {
    let reg = regions[k]
    println("  Region " + (k + 1).to_string() + ": positions " + reg.start.to_string() + "-" + reg.end_.to_string() + ", length=" + reg.length.to_string() + ", avg_score=" + reg.avg_score.to_string() + ", type=" + reg.region_type)
    k = k + 1
  }
  println("")

  // Longest region
  println("=== Longest Disordered Region ===")
  let longest = result.get_longest_region()
  println("  Start: " + longest.start.to_string())
  println("  End: " + longest.end_.to_string())
  println("  Length: " + longest.length.to_string())
  println("  Avg score: " + longest.avg_score.to_string())
  println("  Type: " + longest.region_type)
  println("")

  // Threshold comparison
  println("=== Threshold Comparison ===")
  let low_thresh = @bio.prot_dao_predict(sequence, threshold_disordered=0.3)
  let high_thresh = @bio.prot_dao_predict(sequence, threshold_disordered=0.7)
  println("  Low threshold (0.3): " + low_thresh.get_n_disordered().to_string() + " disordered residues")
  println("  High threshold (0.7): " + high_thresh.get_n_disordered().to_string() + " disordered residues")
  println("")

  // Scores
  println("=== Disorder Scores (first 20) ===")
  let scores = result.get_scores()
  let mut s = ""
  let mut si = 0
  while si < 20 && si < scores.length() {
    s = s + scores[si].to_string() + " "
    si = si + 1
  }
  println("  " + s)
  println("")

  // ASCII visualization
  println("=== ASCII Visualization ===")
  println(result.to_ascii())
}
