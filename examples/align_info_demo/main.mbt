///|
/// Alignment statistics demo inspired by Biopython Bio.Align.AlignInfo.
///
/// Run: moon run examples/align_info_demo/main.mbt
fn main {
  println("=== Bio.Align.AlignInfo Demo ===")

  let seqs = ["ATCGATCGATCG", "ATGGATCGATCG", "ATCGATCGATCA", "ATCGATCGATCG"]

  println("\n1. Input sequences:")
  let mut i = 0
  while i < seqs.length() {
    println("   Seq\{i + 1}: \{seqs[i]}")
    i = i + 1
  }

  println("\n2. Consensus sequence:")
  let consensus = @bio.alignment_consensus(seqs)
  println("   Consensus: \{consensus}")

  println("\n3. Majority sequence:")
  let majority = @bio.majority_sequence(seqs)
  println("   Majority: \{majority}")

  println("\n4. Strict consensus:")
  let strict = @bio.strict_consensus(seqs)
  println("   Strict: \{strict}")

  println("\n5. Conservation profile:")
  let profile = @bio.conservation_profile(seqs)
  let mut j = 0
  while j < profile.length() {
    println("   Position \{j}: \{profile[j]}")
    j = j + 1
  }

  println("\n6. Entropy profile:")
  let entropy = @bio.entropy_profile(seqs)
  let mut k = 0
  while k < entropy.length() {
    println("   Position \{k}: \{entropy[k]}")
    k = k + 1
  }

  println("\n7. Pairwise identity:")
  let identity = @bio.pairwise_identity(seqs)
  for pair in identity {
    println("   \{pair.0} vs \{pair.1}: \{pair.2}%")
  }

  println("\n8. Conserved positions:")
  let conserved = @bio.find_conserved_positions(seqs, min_conservation=0.8)
  for pos in conserved {
    println("   Position \{pos.0}: \{pos.1} (\{pos.2})")
  }

  println("\n9. Variable positions:")
  let variable = @bio.find_variable_positions(seqs, max_conservation=0.5)
  for pos in variable {
    println("   Position \{pos.0}: entropy=\{pos.1}")
  }
}
