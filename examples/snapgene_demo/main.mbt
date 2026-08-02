///|
/// SnapGeneIO module demo - SnapGene .dna binary format parser and writer.
///
/// Run: moon run examples/snapgene_demo --target wasm-gc
fn main {
  println("=== SnapGeneIO Module Demo ===")
  println("SnapGene .dna binary sequence format (packet-based)")
  println("Reference: BioPython Bio.SeqIO.SnapGeneIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample SnapGene file
  // ------------------------------------------------------------------
  println("=== 1. Load sample SnapGene file ===")
  let file = @bio.snapgene_sample()
  println("Loaded sample SnapGene file")
  println("")

  // ------------------------------------------------------------------
  // 2. Print file info
  // ------------------------------------------------------------------
  println("=== 2. File info ===")
  println("Type: \{@bio.SnapgeneSeqType::to_string(file.seq_type)}")
  println("Circular: \{file.is_circular}")
  println("Export Version: \{file.export_version}")
  println("Comments: \{file.comments}")
  println("Sequence Length: \{@bio.snapgene_seq_length(file)} bp")
  println("Features: \{file.features.length()}")
  println("Primers: \{file.primers.length()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show features
  // ------------------------------------------------------------------
  println("=== 3. Features ===")
  for i in 0..<file.features.length() {
    let f = file.features[i]
    let dir_str = if f.direction == 1 {
      "forward"
    } else if f.direction == -1 {
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
  // 4. Show primers
  // ------------------------------------------------------------------
  println("=== 4. Primers ===")
  for i in 0..<file.primers.length() {
    let p = file.primers[i]
    let bind_str = if p.bind_to == 1 { "top strand" } else { "bottom strand" }
    println("Primer \{i}:")
    println("  name: \{p.name}")
    println("  sequence: \{p.sequence}")
    println("  description: \{p.description}")
    println("  binds to: \{bind_str}")
    println("  binds at: \{p.binds_at_position}")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Sequence type enum conversion
  // ------------------------------------------------------------------
  println("=== 5. Sequence type conversion ===")
  let type_names = ["DNA", "RNA", "protein", "unknown"]
  for i in 0..<type_names.length() {
    let t = @bio.SnapgeneSeqType::from_string(type_names[i])
    let name = @bio.SnapgeneSeqType::to_string(t)
    let back = @bio.SnapgeneSeqType::from_string(name)
    println("  \{name} -> round-trip OK: \{back == t}")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Serialize to bytes and back
  // ------------------------------------------------------------------
  println("=== 6. Serialize to bytes and back ===")
  let bytes = @bio.snapgene_to_bytes(file)
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
  let file2 = @bio.snapgene_from_bytes(bytes)
  println("Deserialized: sequence length=\{@bio.snapgene_seq_length(file2)}")
  println("Bytes round-trip verify: \{file2.sequence == file.sequence}")
  println("")

  // ------------------------------------------------------------------
  // 7. Serialize to hex string and back
  // ------------------------------------------------------------------
  println("=== 7. Serialize to hex string and back ===")
  let hex = @bio.snapgene_write(file)
  println("Hex string length: \{hex.length()} chars")
  let hex_preview_len = if hex.length() < 40 { hex.length() } else { 40 }
  println("  first \{hex_preview_len} chars: \{hex[0:hex_preview_len]}")
  let file3 = @bio.snapgene_parse(hex)
  println("Deserialized from hex: sequence length=\{@bio.snapgene_seq_length(file3)}")
  println("Hex round-trip verify: \{file3.sequence == file.sequence}")
  println("")

  // ------------------------------------------------------------------
  // 8. Build a custom SnapGene file
  // ------------------------------------------------------------------
  println("=== 8. Build a custom SnapGene file ===")
  let custom = @bio.SnapgeneFile::new()
  custom.set_sequence("ATGAAATAA")
  custom.set_circular(true)
  custom.set_comments("Custom demo plasmid")
  custom.set_export_version("2.0")
  custom.set_seq_type(@bio.SnapgeneSeqType::from_string("DNA"))
  let q : Map[String, String] = Map([], capacity=0)
  q["label"] = "demo"
  custom.add_feature(
    @bio.SnapgeneFeature::new(
      name="gene1",
      type_="CDS",
      direction=1,
      segments=[(0, 8)],
      qualifiers=q,
    ),
  )
  custom.add_primer(
    @bio.SnapgenePrimer::new(
      name="M13-F",
      sequence="GTAAAACGACGGCCAG",
      description="M13 forward primer",
      bind_to=1,
      binds_at_position=0,
    ),
  )
  println("Custom file built:")
  println("  sequence: \{custom.sequence}")
  println("  circular: \{custom.is_circular}")
  println("  features: \{custom.features.length()}")
  println("  primers: \{custom.primers.length()}")
  let custom_hex = @bio.snapgene_write(custom)
  let custom_parsed = @bio.snapgene_parse(custom_hex)
  println("Round-trip verify: \{custom_parsed.sequence == custom.sequence}")
  println("  parsed features: \{custom_parsed.features.length()}")
  println("  parsed primers: \{custom_parsed.primers.length()}")
  println("")

  // ------------------------------------------------------------------
  // 9. Convert to SeqRecord
  // ------------------------------------------------------------------
  println("=== 9. Convert to SeqRecord ===")
  let rec = @bio.snapgene_to_seqrecord(file)
  println("SeqRecord id: \{rec.id}")
  println("SeqRecord description: \{rec.description}")
  println("SeqRecord sequence length: \{rec.seq.data.length()}")
  println("")

  // ------------------------------------------------------------------
  // 10. Summary
  // ------------------------------------------------------------------
  println("=== 10. Summary ===")
  println(@bio.snapgene_summary(file))
  println("=== Demo complete ===")
}
