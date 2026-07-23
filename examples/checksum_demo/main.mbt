///|
/// CheckSum demo.
///
/// Run: moon run examples/checksum_demo/main.mbt
fn main {
  println("=== CheckSum Demo ===\n")

  let seq1 = "ACDEFGHIKLMNPQRSTVWY"
  let seq2 = "ACGT"
  let seq3 = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

  println("1. GCG checksum:")
  let result1 = @bio.checksum_gcg(seq1)
  println("   Seq: \{seq1}")
  println("   Checksum: \{result1.checksum}")
  println("   Method: \{result1.method}")

  println("\n2. SEGUID checksum:")
  let result2 = @bio.checksum_seguid(seq2)
  println("   Seq: \{seq2}")
  println("   Checksum: \{result2.checksum}")
  println("   Method: \{result2.method}")

  let result3 = @bio.checksum_seguid(seq3)
  println("\n   Seq (long): \{seq3.length()} bp")
  println("   Checksum: \{result3.checksum}")
  println("   Method: \{result3.method}")

  println("\n3. Verify checksum:")
  let verified_gcg = @bio.verify_checksum(seq1, result1.checksum, "gcg")
  println("   GCG verification: \{verified_gcg}")
  let verified_seguid = @bio.verify_checksum(seq2, result2.checksum, "seguid")
  println("   SEGUID verification: \{verified_seguid}")

  println("\n4. Example sequence:")
  let example_seq = @bio.create_example_sequence()
  println("   Example: \{example_seq}")

  println("\n=== Demo Complete ===")
}