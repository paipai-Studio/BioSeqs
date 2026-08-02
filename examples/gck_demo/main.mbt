///|
/// GckIO module demo - Gene Construction Kit (GCK) binary format parser and writer.
///
/// Run: moon run examples/gck_demo --target wasm-gc
fn main {
  println("=== GckIO Module Demo ===")
  println("Gene Construction Kit (GCK) binary sequence format")
  println("Reference: BioPython Bio.SeqIO.GckIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample GCK file
  // ------------------------------------------------------------------
  println("=== 1. Load sample GCK file ===")
  let file = @bio.gck_sample()
  println("Loaded sample GCK file")
  println("")

  // ------------------------------------------------------------------
  // 2. Print file info
  // ------------------------------------------------------------------
  println("=== 2. File info ===")
  println("Version: \{file.version}")
  println("Type: \{@bio.GckSeqType::to_string(file.seq_type)}")
  println("Circular: \{file.is_circular}")
  println("Sequence Length: \{@bio.gck_seq_length(file)} bp")
  println("Features: \{file.features.length()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show features
  // ------------------------------------------------------------------
  println("=== 3. Features ===")
  for i in 0..<file.features.length() {
    let f = file.features[i]
    let dir_str = if f.direction == 1 {
      "forward"
    } else if f.direction == 2 {
      "reverse"
    } else {
      "none"
    }
    println("Feature \{i}:")
    println("  name: \{f.name}")
    println("  type: \{f.type_}")
    println("  direction: \{dir_str}")
    println("  segments: \{f.segments.length()}")
    for j in 0..<f.segments.length() {
      let (start, end) = f.segments[j]
      println("    segment \{j}: \{start}..\{end}")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Count features by type
  // ------------------------------------------------------------------
  println("=== 4. Count features by type ===")
  let cds_count = @bio.gck_count_features_by_type(file, "CDS")
  let promoter_count = @bio.gck_count_features_by_type(file, "promoter")
  let misc_count = @bio.gck_count_features_by_type(file, "misc")
  println("CDS: \{cds_count}")
  println("promoter: \{promoter_count}")
  println("misc: \{misc_count}")
  println("")

  // ------------------------------------------------------------------
  // 5. Sequence type enum conversion
  // ------------------------------------------------------------------
  println("=== 5. Sequence type conversion ===")
  let type_codes = [0, 1, 2, 3]
  for i in 0..<type_codes.length() {
    let t = @bio.GckSeqType::from_int(type_codes[i])
    let code = @bio.GckSeqType::to_int(t)
    let name = @bio.GckSeqType::to_string(t)
    let back = @bio.GckSeqType::from_int(code)
    println("  \{name} (code=\{code}) -> round-trip OK: \{back == t}")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Serialize to bytes and back
  // ------------------------------------------------------------------
  println("=== 6. Serialize to bytes and back ===")
  let bytes = @bio.gck_to_bytes(file)
  println("Serialized to \{bytes.length()} byte(s)")
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
  // Verify magic bytes
  println("Magic: \{bytes[0].to_char()} \{bytes[1].to_char()} \{bytes[2].to_char()} (should be G C K)")
  let file2 = @bio.gck_from_bytes(bytes)
  println("Deserialized: sequence length=\{@bio.gck_seq_length(file2)}")
  println("Bytes round-trip verify: \{file2.sequence == file.sequence}")
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize to hex string and back
  // ------------------------------------------------------------------
  println("=== 7. Serialize to hex string and back ===")
  let hex = @bio.gck_write(file)
  println("Hex string length: \{hex.length()} chars")
  let hex_preview_len = if hex.length() < 40 { hex.length() } else { 40 }
  println("  first \{hex_preview_len} chars: \{hex[0:hex_preview_len]}")
  let file3 = @bio.gck_parse(hex)
  println("Deserialized from hex: sequence length=\{@bio.gck_seq_length(file3)}")
  println("Hex round-trip verify: \{file3.sequence == file.sequence}")
  println("")

  // ------------------------------------------------------------------
  // 8. Build a custom GCK file
  // ------------------------------------------------------------------
  println("=== 8. Build a custom GCK file ===")
  let custom = @bio.GckFile::new()
  custom.set_sequence("ATGAAATAAGCATGCATGCATGC")
  custom.set_circular(false)
  custom.set_version(1)
  custom.set_seq_type(@bio.GckSeqType::from_int(0))
  custom.add_feature(
    @bio.GckFeature::new(
      name="geneA",
      type_="CDS",
      direction=1,
      segments=[(0, 8)],
    ),
  )
  custom.add_feature(
    @bio.GckFeature::new(
      name="promoter1",
      type_="promoter",
      direction=1,
      segments=[(10, 22)],
    ),
  )
  println("Custom file built:")
  println("  sequence: \{custom.sequence}")
  println("  circular: \{custom.is_circular}")
  println("  features: \{custom.features.length()}")
  let custom_hex = @bio.gck_write(custom)
  let custom_parsed = @bio.gck_parse(custom_hex)
  println("Round-trip verify: \{custom_parsed.sequence == custom.sequence}")
  println("  parsed features: \{custom_parsed.features.length()}")
  println("  parsed circular: \{custom_parsed.is_circular}")
  println("")

  // ------------------------------------------------------------------
  // 9. Convert to SeqRecord
  // ------------------------------------------------------------------
  println("=== 9. Convert to SeqRecord ===")
  let rec = @bio.gck_to_seqrecord(file)
  println("SeqRecord id: \{rec.id}")
  println("SeqRecord description: \{rec.description}")
  println("SeqRecord sequence length: \{rec.seq.data.length()}")
  // Convert back
  let back = @bio.gck_from_seqrecord(rec)
  println("Back-converted GckFile: sequence length=\{@bio.gck_seq_length(back)}")
  println("")

  // ------------------------------------------------------------------
  // 10. Summary
  // ------------------------------------------------------------------
  println("=== 10. Summary ===")
  println(@bio.gck_summary(file))
  println("=== Demo complete ===")
}
