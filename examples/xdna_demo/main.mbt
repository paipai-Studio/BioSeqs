///|
/// XDNA module demo - Geneious XDNA binary sequence format parser and writer.
///
/// Run: moon run examples/xdna_demo
fn main {
  println("=== XDNA Module Demo ===")
  println("Geneious XDNA binary sequence format (2-bit packed DNA)")
  println("Reference: BioPython Bio.SeqIO.XdnaIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample XDNA file
  // ------------------------------------------------------------------
  println("=== 1. Load sample XDNA file ===")
  let file = @bio.xdna_sample_file()
  println("Loaded sample XDNA file")
  println("")

  // ------------------------------------------------------------------
  // 2. Print file info (n_records, version)
  // ------------------------------------------------------------------
  println("=== 2. File info ===")
  println("n_records: \{file.n_records()}")
  println("version: \{file.version()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show all records (name, sequence, type, checksum)
  // ------------------------------------------------------------------
  println("=== 3. Records ===")
  let records = file.records()
  for i in 0..<records.length() {
    let rec = records[i]
    println("Record \{i}:")
    println("  name: \{rec.name()}")
    println("  sequence: \{rec.sequence()}")
    println("  type: \{rec.seq_type().to_string()} (\{rec.seq_type().description()})")
    println("  checksum: \{rec.checksum()}")
    println("  annotations: '\{rec.annotations()}'")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Demonstrate base encoding (A->0, C->1, G->2, T->3)
  // ------------------------------------------------------------------
  println("=== 4. Base encoding (A->0, C->1, G->2, T->3) ===")
  let bases = ['A', 'C', 'G', 'T']
  for i in 0..<bases.length() {
    let code = @bio.xdna_base_to_code(bases[i])
    println("  \{bases[i]} -> \{code}")
  }
  // Round-trip: code back to base.
  println("Reverse (code -> base):")
  for c in 0..<4 {
    let base = @bio.xdna_code_to_base(c)
    println("  \{c} -> \{base}")
  }
  // Invalid base handling.
  let bad_code = @bio.xdna_base_to_code('N')
  println("  N -> \{bad_code} (invalid, returns -1)")
  println("")

  // ------------------------------------------------------------------
  // 5. Pack/unpack round-trip for a DNA sequence
  // ------------------------------------------------------------------
  println("=== 5. Pack/unpack round-trip ===")
  let dna = "ATGCATGC"
  let packed = @bio.xdna_pack_sequence(dna)
  println("Sequence: \{dna}")
  println("Packed bytes (\{packed.length()} byte(s)):")
  let pack_sb = StringBuilder::new()
  pack_sb.write_string("  ")
  for i in 0..<packed.length() {
    if i > 0 {
      pack_sb.write_string(", ")
    }
    pack_sb.write_string(packed[i].to_string())
  }
  println(pack_sb.to_string())
  let unpacked = @bio.xdna_unpack_sequence(packed, dna.length())
  println("Unpacked: \{unpacked}")
  println("Round-trip OK: \{unpacked == dna}")
  println("")

  // ------------------------------------------------------------------
  // 6. Serialize to bytes and back, verify
  // ------------------------------------------------------------------
  println("=== 6. Serialize to bytes and back ===")
  let bytes = @bio.xdna_to_bytes(file)
  println("Serialized to \{bytes.length()} byte(s)")
  // Show first 20 bytes as a preview.
  let preview_len = if bytes.length() < 20 { bytes.length() } else { 20 }
  let byte_sb = StringBuilder::new()
  byte_sb.write_string("  first ")
  byte_sb.write_string(preview_len.to_string())
  byte_sb.write_string(" byte(s): ")
  for i in 0..<preview_len {
    if i > 0 {
      byte_sb.write_string(", ")
    }
    byte_sb.write_string(bytes[i].to_string())
  }
  println(byte_sb.to_string())
  let file2 = @bio.xdna_from_bytes(bytes)
  println("Deserialized: n_records=\{file2.n_records()}, version=\{file2.version()}")
  let recs2 = file2.records()
  let mut byte_ok = true
  if recs2.length() != records.length() {
    byte_ok = false
  } else {
    for i in 0..<recs2.length() {
      if recs2[i].name() != records[i].name() ||
        recs2[i].sequence() != records[i].sequence() {
        byte_ok = false
      }
    }
  }
  println("Bytes round-trip verify: \{byte_ok}")
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize to hex string and back, verify
  // ------------------------------------------------------------------
  println("=== 7. Serialize to hex string and back ===")
  let hex = @bio.xdna_write(file)
  println("Hex string length: \{hex.length()} chars")
  // Show first 40 hex chars as a preview.
  let hex_preview_len = if hex.length() < 40 { hex.length() } else { 40 }
  println("  first \{hex_preview_len} chars: \{hex[0:hex_preview_len]}")
  let file3 = @bio.xdna_read(hex)
  println("Deserialized from hex: n_records=\{file3.n_records()}")
  let recs3 = file3.records()
  let mut hex_ok = true
  if recs3.length() != records.length() {
    hex_ok = false
  } else {
    for i in 0..<recs3.length() {
      if recs3[i].name() != records[i].name() ||
        recs3[i].sequence() != records[i].sequence() {
        hex_ok = false
      }
    }
  }
  println("Hex round-trip verify: \{hex_ok}")
  println("")

  // ------------------------------------------------------------------
  // 8. Convert to SeqRecord array and show
  // ------------------------------------------------------------------
  println("=== 8. Convert to SeqRecord array ===")
  let seq_records = @bio.xdna_to_seq_records(file)
  println("Converted \{seq_records.length()} XdnaRecord(s) to SeqRecord(s)")
  for i in 0..<seq_records.length() {
    let rec = seq_records[i]
    println("  Record \{i}: id=\{rec.id}, name=\{rec.name}, length=\{rec.seq.to_string().length()}")
  }
  println("")

  // Bonus: build an XdnaFile from SeqRecords using XdnaSeqType::from_int.
  println("=== Bonus: Build XdnaFile from SeqRecords ===")
  let new_records = [
    @bio.SeqRecord::new(@bio.Seq::new("ACGTACGT"), id="built1"),
    @bio.SeqRecord::new(@bio.Seq::new("TTTTGGGG"), id="built2"),
  ]
  let built_file = @bio.xdna_from_seq_records(
    new_records,
    @bio.XdnaSeqType::from_int(0),
  )
  println("Built file: n_records=\{built_file.n_records()}")
  let built_recs = built_file.records()
  for i in 0..<built_recs.length() {
    println("  Record \{i}: name=\{built_recs[i].name()}, sequence=\{built_recs[i].sequence()}, type=\{built_recs[i].seq_type().to_string()}")
  }
  println("")

  println("=== Demo complete ===")
}
