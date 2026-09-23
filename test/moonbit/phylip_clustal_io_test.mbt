///|
/// Round-trip tests for the PHYLIP and ClustalW alignment I/O modules
/// (`phylip_io` / `clustal_io`), operating on `MultipleSeqAlignment`.

///|
fn make_msa() -> @src.MultipleSeqAlignment {
  let rec = fn(id : String, seq : String) -> @src.SeqRecord {
    @src.SeqRecord::new(@src.Seq::new(seq), id=id, name=id)
  }
  @src.MultipleSeqAlignment::new([
    rec("seq1", "ACGTACGTACGT"),
    rec("seq2", "ACGTTCGTACGT"),
    rec("seq3", "ACGTACGTTGGT"),
  ]) catch {
    AlignIOError(m) => abort("MSA 构造失败: " + m)
  }
}

///|
fn parse_phylip_ok(text : String) -> @src.MultipleSeqAlignment {
  @src.parse_phylip(text) catch {
    AlignIOError(m) => abort("parse_phylip 失败: " + m)
  }
}

///|
test "phylip_io_interlaced_roundtrip" {
  let msa = make_msa()
  let text = @src.write_phylip(msa)
  let back = parse_phylip_ok(text)
  assert_eq(back.num_records(), 3)
  assert_eq(back.get_alignment_length(), 12)
  let r = back.iter()
  assert_eq(r[0].id, "seq1")
  assert_eq(r[0].seq.to_string(), "ACGTACGTACGT")
  assert_eq(r[2].id, "seq3")
  assert_eq(r[2].seq.to_string(), "ACGTACGTTGGT")
}

///|
test "phylip_io_sequential_roundtrip" {
  let msa = make_msa()
  let text = @src.write_phylip_sequential(msa)
  let back = @src.parse_phylip_sequential(text) catch {
    AlignIOError(m) => abort("parse_phylip_sequential 失败: " + m)
  }
  assert_eq(back.num_records(), 3)
  assert_eq(back.get_alignment_length(), 12)
  let r = back.iter()
  assert_eq(r[1].id, "seq2")
  assert_eq(r[1].seq.to_string(), "ACGTTCGTACGT")
}

///|
test "clustal_io_roundtrip" {
  let msa = make_msa()
  let text = @src.write_clustal(msa)
  let back = @src.parse_clustal(text) catch {
    AlignIOError(m) => abort("parse_clustal 失败: " + m)
  }
  assert_eq(back.num_records(), 3)
  assert_eq(back.get_alignment_length(), 12)
  let r = back.iter()
  assert_eq(r[0].id, "seq1")
  assert_eq(r[2].seq.to_string(), "ACGTACGTTGGT")
}