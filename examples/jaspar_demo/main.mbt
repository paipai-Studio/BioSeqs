///|
/// JASPAR PFM Format Parsing Demo.
///
/// Demonstrates parsing JASPAR Position Frequency Matrices (PFM),
/// computing Position Weight Matrices (PWM), and scanning sequences.
///
/// Run: moon run examples/jaspar_demo/main.mbt
fn main {
  println("=== JASPAR PFM Format Demo ===\n")

  // Example JASPAR PFM content
  let jaspar_content = ">MA0001.1\tFOXA1\n" +
    "    A  [    10    15    20     5    10 ]\n" +
    "    C  [    20    10     5    15    20 ]\n" +
    "    G  [    30    25    30    35    25 ]\n" +
    "    T  [    40    50    45    45    45 ]\n" +
    "    >MA0002.1\tCTCF\n" +
    "    A  [    20    30    10    15    20 ]\n" +
    "    C  [    10    20    30    25    15 ]\n" +
    "    G  [    40    35    25    30    35 ]\n" +
    "    T  [    30    15    35    30    30 ]\n"

  // Parse JASPAR
  println("1. Parse JASPAR PFM content:")
  let matrices = @bio.parse_jaspar(jaspar_content)
  println("   Found \{matrices.length()} matrices")

  for matrix in matrices {
    println("\n2. Matrix: \{matrix.name} (\{matrix.id})")
    println("   Width: \{matrix.width}")
    println("   Total counts: \{matrix.total_counts}")

    // Consensus sequence
    let consensus = matrix.consensus()
    println("   Consensus: \{consensus}")

    // Information content
    let ic = matrix.total_information()
    println("   Total information content: \{ic} bits")

    // PWM
    println("   PWM (first 3 columns):")
    let pwm = matrix.to_pwm(pseudocount=0.01)
    println("     A: \{pwm[0][0]}, \{pwm[0][1]}, \{pwm[0][2]}")
    println("     C: \{pwm[1][0]}, \{pwm[1][1]}, \{pwm[1][2]}")
    println("     G: \{pwm[2][0]}, \{pwm[2][1]}, \{pwm[2][2]}")
    println("     T: \{pwm[3][0]}, \{pwm[3][1]}, \{pwm[3][2]}")

    // Scan sequence
    println("   Sequence scanning:")
    let test_seq = "ACGTACGTACGTACGT"
    let hits = matrix.scan_sequence(test_seq)
    println("     Sequence: \{test_seq}")
    println("     Hits found: \{hits.length()}")
  }

  println("\n=== Demo Complete ===")
}
