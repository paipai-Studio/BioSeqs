///|
/// IgIO module demo - IntelliGenetics (MASE) sequence format parser and writer.
///
/// Run: moon run examples/ig_demo
fn main {
  println("=== IgIO Module Demo ===")
  println("IntelliGenetics (MASE) sequence format parser and writer")
  println("Reference: Biopython Bio.SeqIO.IgIO module")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample IntelliGenetics records
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Records ===")
  let text = @bio.ig_sample_text()
  println("Sample IG file text:")
  println(text)
  let records = @bio.ig_parse(text)
  println("Parsed \{records.length()} record(s)")
  println("")

  // ------------------------------------------------------------------
  // 2. Display each record
  // ------------------------------------------------------------------
  println("=== 2. Record Details ===")
  for i in 0..<records.length() {
    let r = records[i]
    println("Record \{i + 1}:")
    println("  Title: \{r.title()}")
    println("  Comment: \{r.comment()}")
    println("  Sequence length: \{r.length()} nt")
    let seq = r.sequence()
    let preview_len = if seq.length() > 40 { 40 } else { seq.length() }
    println("  Sequence preview: \{seq.substring(start=0, end=preview_len)}...")
    println("")
  }

  // ------------------------------------------------------------------
  // 3. Parse a custom IG string
  // ------------------------------------------------------------------
  println("=== 3. Parse Custom IG String ===")
  let custom = ";;\n;; Custom IntelliGenetics file\n;HIV-1 group M\n;subtype B\nB_HXB2R\nATGGGTGCTGTGATGAAGCAGAGACGGGTAAGAGCTC\n1\n"
  println("Input:")
  println(custom)
  let custom_records = @bio.ig_parse(custom)
  println("Parsed \{custom_records.length()} record(s)")
  if custom_records.length() > 0 {
    let r = custom_records[0]
    println("  Title: \{r.title()}")
    println("  Comment: \{r.comment()}")
    println("  Sequence: \{r.sequence()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Serialize a record to IG format
  // ------------------------------------------------------------------
  println("=== 4. Serialize Record to IG Format ===")
  let r = @bio.ig_sample_record()
  println("Original record: \{r.to_string()}")
  let ig_text = @bio.ig_record_to_string(r)
  println("Serialized IG format:")
  println(ig_text)
  println("")

  // ------------------------------------------------------------------
  // 5. Serialize multiple records with file header
  // ------------------------------------------------------------------
  println("=== 5. Serialize Multiple Records ===")
  let all_records = @bio.ig_sample_records()
  let header = @bio.ig_sample_file_header()
  let full_text = @bio.ig_records_to_string(all_records, Some(header))
  println("Full IG file output:")
  println(full_text)
  println("")

  // ------------------------------------------------------------------
  // 6. Round-trip: serialize -> parse -> compare
  // ------------------------------------------------------------------
  println("=== 6. Round-Trip Verification ===")
  let original = @bio.ig_sample_record_2()
  let serialized = @bio.ig_record_to_string(original)
  let parsed = @bio.ig_parse(serialized)
  if parsed.length() > 0 {
    let round_tripped = parsed[0]
    let title_ok = round_tripped.title() == original.title()
    let seq_ok = round_tripped.sequence() == original.sequence()
    let comment_ok = round_tripped.comment() == original.comment()
    println("Title match: \{title_ok}")
    println("Sequence match: \{seq_ok}")
    println("Comment match: \{comment_ok}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Convert to SeqRecord
  // ------------------------------------------------------------------
  println("=== 7. SeqRecord Interop ===")
  let ig_rec = @bio.ig_sample_record()
  let seq_rec = @bio.ig_to_seq_record(ig_rec)
  println("Converted to SeqRecord:")
  println("  id: \{seq_rec.id}")
  println("  name: \{seq_rec.name}")
  println("  seq: \{seq_rec.seq.to_string()}")
  println("")

  // Convert back
  let back = @bio.ig_from_seq_record(seq_rec)
  println("Converted back to IgRecord:")
  println("  title: \{back.title()}")
  println("  sequence: \{back.sequence()}")
  println("")

  // ------------------------------------------------------------------
  // 8. Custom line width for sequence wrapping
  // ------------------------------------------------------------------
  println("=== 8. Custom Line Width ===")
  let long_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
  let long_rec = @bio.IgRecord::new("LONG_SEQ", "Long sequence demo", long_seq)
  println("Sequence length: \{long_seq.length()} nt")
  println("Wrapped at 20 chars:")
  let wrapped = @bio.ig_record_to_string_width(long_rec, 20)
  println(wrapped)
  println("")

  println("=== Demo complete ===")
}
