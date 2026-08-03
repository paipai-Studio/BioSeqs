///|
/// NibIO module demo - UCSC nib 2-bit binary sequence format reader/writer.
///
/// Run: moon run examples/nib_io_demo
fn main {
  println("=== NibIO Module Demo ===")
  println("UCSC nib 2-bit binary sequence format")
  println("Reference: BioPython Bio.SeqIO.NibIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Base encoding table (T=0, C=1, A=2, G=3)
  // ------------------------------------------------------------------
  println("=== 1. Base encoding (T=0, C=1, A=2, G=3) ===")
  let bases = ['T', 'C', 'A', 'G']
  for i in 0..<bases.length() {
    let code = @bio.nib_base_to_code(bases[i])
    println("  \{bases[i]} -> \{code}")
  }
  println("Reverse (code -> base):")
  for c in 0..<4 {
    let base = @bio.nib_code_to_base(c)
    println("  \{c} -> \{base}")
  }
  let bad = @bio.nib_base_to_code('N')
  println("  N -> \{bad} (invalid, returns -1)")
  println("")

  // ------------------------------------------------------------------
  // 2. Create a NibSequence from a DNA string
  // ------------------------------------------------------------------
  println("=== 2. Create a NibSequence ===")
  let dna = "ATGCATGCATGC"
  let nib = @bio.NibSequence::new("demo1", dna)
  println("Sequence: \{dna}")
  println(nib.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Display hex encoding of the packed data
  // ------------------------------------------------------------------
  println("=== 3. Hex encoding ===")
  let hex = @bio.nib_to_hex(nib)
  println("Packed hex: \{hex}")
  println("Packed bytes: \{nib.packed_data().length()}")
  let pb = nib.packed_data()
  let pack_sb = StringBuilder::new()
  pack_sb.write_string("  bytes: ")
  for i in 0..<pb.length() {
    if i > 0 {
      pack_sb.write_string(", ")
    }
    pack_sb.write_string(pb[i].to_string())
  }
  println(pack_sb.to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Decode back to a DNA string
  // ------------------------------------------------------------------
  println("=== 4. Decode back ===")
  let decoded = @bio.nib_to_string(nib)
  println("Decoded: \{decoded}")
  println("Round-trip OK: \{decoded == dna}")
  println("")

  // ------------------------------------------------------------------
  // 5. Extract a subsequence
  // ------------------------------------------------------------------
  println("=== 5. Subsequence extraction ===")
  let sub = @bio.nib_subsequence(nib, 2, 6)
  println("Subsequence (start=2, length=6): \{sub}")
  println("Single base at index 5: \{@bio.nib_get_base(nib, 5)}")
  println("")

  // ------------------------------------------------------------------
  // 6. Reverse complement
  // ------------------------------------------------------------------
  println("=== 6. Reverse complement ===")
  let rc = @bio.nib_reverse_complement(nib)
  println("Original:  \{dna}")
  println("Rev-comp:  \{@bio.nib_to_string(rc)}")
  println("Rev-comp id: \{rc.seq_id()}")
  println("")

  // ------------------------------------------------------------------
  // 7. GC content
  // ------------------------------------------------------------------
  println("=== 7. GC content ===")
  let gc = @bio.nib_gc_content(nib)
  println("GC fraction: \{gc}")
  println("GC percentage: \{gc * 100.0}%")
  println("")

  // ------------------------------------------------------------------
  // 8. Compression ratio
  // ------------------------------------------------------------------
  println("=== 8. Compression ratio ===")
  let ratio = @bio.nib_compression_ratio(nib)
  println("Original size: \{@bio.nib_size(nib)} bytes (1 byte/base)")
  println("Packed size: \{@bio.nib_compressed_size(nib)} bytes (4 bases/byte)")
  println("Compression ratio: \{ratio}")
  println("")

  // ------------------------------------------------------------------
  // 9. Hex round-trip
  // ------------------------------------------------------------------
  println("=== 9. Hex round-trip ===")
  let nib2 = @bio.nib_from_hex("from_hex", hex, dna.length())
  println("From hex: \{@bio.nib_to_string(nib2)}")
  println("Hex round-trip OK: \{@bio.nib_to_string(nib2) == dna}")
  println("")

  // ------------------------------------------------------------------
  // 10. Lowercase handling
  // ------------------------------------------------------------------
  println("=== 10. Lowercase handling ===")
  let lower = @bio.NibSequence::new("lower", "atgcatgcatgc")
  println("Lowercase input: atgcatgcatgc")
  println("is_uppercase flag: \{lower.is_uppercase()}")
  println("Decoded: \{@bio.nib_to_string(lower)}")
  println("Hex (same as uppercase): \{@bio.nib_to_hex(lower)}")
  println("")

  // ------------------------------------------------------------------
  // 11. Load sample data
  // ------------------------------------------------------------------
  println("=== 11. Sample data ===")
  let sample = @bio.nib_sample_data()
  println(sample.to_string())
  println("Decoded: \{@bio.nib_to_string(sample)}")
  println("Hex: \{@bio.nib_to_hex(sample)}")
  println("GC content: \{@bio.nib_gc_content(sample) * 100.0}%")
  println("")

  // ------------------------------------------------------------------
  // 12. Summary of odd-length encoding
  // ------------------------------------------------------------------
  println("=== 12. Odd-length encoding ===")
  let odd = @bio.NibSequence::new("odd", "ATGCA")
  println("Sequence 'ATGCA' (5 bases, not a multiple of 4)")
  println("Packed bytes: \{@bio.nib_compressed_size(odd)} (ceil(5/4) = 2)")
  println("Hex: \{@bio.nib_to_hex(odd)}")
  println("Decoded: \{@bio.nib_to_string(odd)}")
  println("Compression ratio: \{@bio.nib_compression_ratio(odd)}")
  println("")

  println("=== Demo complete ===")
}
