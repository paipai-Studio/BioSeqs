///|
/// TwoBitIO demo - UCSC .2bit binary genome format reader/writer.
///
/// Run: moon run examples/twobit_demo
fn main {
  println("=== TwoBitIO Demo ===")
  println("UCSC 2bit binary sequence format")
  println("Reference: UCSC Genome Browser 2bit format specification")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample 2bit file
  // ------------------------------------------------------------------
  println("=== 1. Load sample 2bit file ===")
  let tbf = @bio.twobit_sample_file()
  println("Loaded sample 2bit file")
  println("")

  // ------------------------------------------------------------------
  // 2. Print file info (n_sequences, sequence names, dna sizes)
  // ------------------------------------------------------------------
  println("=== 2. File info ===")
  println("n_sequences: \{tbf.n_sequences()}")
  let names = tbf.sequence_names()
  println("sequence names count: \{names.length()}")
  for n in names {
    match tbf.get_sequence(n) {
      Some(s) => println("  \{n}: dna_size=\{s.dna_size()}")
      None => println("  \{n}: not found")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Extract sequences as strings and print
  // ------------------------------------------------------------------
  println("=== 3. Extract sequences as strings ===")
  for n in names {
    match @bio.twobit_get_sequence_string(tbf, n) {
      Some(s) => println("  \{n}: \{s}")
      None => println("  \{n}: not found")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Show N blocks for seq2
  // ------------------------------------------------------------------
  println("=== 4. N blocks for seq2 ===")
  match tbf.get_sequence("seq2") {
    Some(seq2) => {
      let n_blocks = seq2.n_blocks()
      println("seq2 has \{n_blocks.length()} N-block(s)")
      for b in n_blocks {
        println("  start=\{b.start()}, size=\{b.size()}, end=\{b.end()}")
      }
    }
    None => println("seq2 not found")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Demonstrate pack/unpack round-trip
  // ------------------------------------------------------------------
  println("=== 5. Pack/unpack round-trip ===")
  let dna = "ATGCNNNNATGC"
  let (packed, n_blocks, mask_blocks) = @bio.twobit_pack_sequence(dna)
  let unpacked = @bio.twobit_unpack_sequence(
    packed,
    dna.length(),
    n_blocks,
    mask_blocks,
  )
  println("Original:  \{dna}")
  println("Unpacked:  \{unpacked}")
  println("Match: \{dna == unpacked}")
  println("")

  // ------------------------------------------------------------------
  // 6. Serialize to bytes and back, verify
  // ------------------------------------------------------------------
  println("=== 6. Serialize to bytes and back ===")
  let bytes = @bio.twobit_to_bytes(tbf)
  println("Serialized to \{bytes.length()} bytes")
  let restored = @bio.twobit_from_bytes(bytes)
  println("Restored n_sequences: \{restored.n_sequences()}")
  for n in restored.sequence_names() {
    match @bio.twobit_get_sequence_string(restored, n) {
      Some(s) => println("  \{n}: \{s}")
      None => println("  \{n}: not found")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize to hex string and back, verify
  // ------------------------------------------------------------------
  println("=== 7. Serialize to hex string and back ===")
  let hex = @bio.twobit_write(tbf)
  println("Hex string length: \{hex.length()}")
  let hex_restored = @bio.twobit_read(hex)
  println("Restored n_sequences: \{hex_restored.n_sequences()}")
  for n in hex_restored.sequence_names() {
    match @bio.twobit_get_sequence_string(hex_restored, n) {
      Some(s) => println("  \{n}: \{s}")
      None => println("  \{n}: not found")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Convert to SeqRecord array and show
  // ------------------------------------------------------------------
  println("=== 8. Convert to SeqRecord array ===")
  let recs = @bio.twobit_to_seq_records(tbf)
  println("Converted to \{recs.length()} SeqRecord(s)")
  for rec in recs {
    println("  \{rec.repr()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 9. Create from a DNA string and verify
  // ------------------------------------------------------------------
  println("=== 9. Create from a DNA string ===")
  let new_dna = "GGGATCCCATGC"
  let new_rec = @bio.SeqRecord::new(
    @bio.Seq::new(new_dna),
    id="custom",
    name="custom",
    description="Created from a DNA string",
  )
  let new_file = @bio.twobit_from_seq_records([new_rec])
  match @bio.twobit_get_sequence_string(new_file, "custom") {
    Some(s) => {
      println("Original:  \{new_dna}")
      println("Extracted: \{s}")
      println("Match: \{new_dna == s}")
    }
    None => println("custom sequence not found")
  }
  println("")

  println("=== Demo complete ===")
}
