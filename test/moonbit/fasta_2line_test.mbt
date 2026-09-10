// Whitebox tests for the fasta-2line format (spec 5.2.1).

///|
test "fasta_2line_parse_two_records" {
  let content = ">rec1 First record\nMVKVYAPASS\n>rec2\nACGT"
  let records = @src.parse_fasta_2line(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "rec1")
  assert_eq(records[0].description, "rec1 First record")
  assert_eq(records[0].seq.to_string(), "MVKVYAPASS")
  assert_eq(records[1].id, "rec2")
  assert_eq(records[1].description, "rec2")
  assert_eq(records[1].seq.to_string(), "ACGT")
}

///|
test "fasta_2line_write_single_line_sequence" {
  // A long sequence must be written on exactly one line (no wrapping)
  let long_seq = "M".repeat(200)
  let rec = @src.SeqRecord::new(
    @src.Seq::new(long_seq),
    id="long1",
    description="long1 Long sequence",
  )
  let out = @src.write_fasta_2line([rec])
  let lines = out.split("\n").to_array()
  // header, one sequence line, and the trailing newline's empty remainder
  assert_eq(lines.length(), 3)
  assert_eq(lines[0], ">long1 Long sequence")
  assert_eq(lines[1], long_seq)
  assert_eq(lines[2], "")
}

///|
test "fasta_2line_roundtrip_semantic_equivalence" {
  let original = ">r1 Alpha\nMVKV\n>r2 Beta\nACGTACGT"
  let records = @src.parse_fasta_2line(original)
  let written = @src.write_fasta_2line(records)
  let reparsed = @src.parse_fasta_2line(written)
  assert_eq(reparsed.length(), records.length())
  for i = 0; i < records.length(); i = i + 1 {
    assert_eq(reparsed[i].id, records[i].id)
    assert_eq(reparsed[i].seq.to_string(), records[i].seq.to_string())
  }
  assert_eq(written, ">r1 Alpha\nMVKV\n>r2 Beta\nACGTACGT\n")
}

///|
test "fasta_2line_unexpected_wrap_handled_as_fasta" {
  // Documented trade-off (spec 5.2.1.1b): an unexpectedly wrapped
  // sequence is accepted and concatenated like standard multi-line FASTA.
  let wrapped = ">w1 Wrapped\nAAA\nTTT"
  let records = @src.parse_fasta_2line(wrapped)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AAATTT")
}