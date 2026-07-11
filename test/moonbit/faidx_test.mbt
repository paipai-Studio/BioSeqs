///|
test "faidx_build_index" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let idx = @src.build_fai(fasta)
  assert_eq(idx.num_sequences(), 2)
  assert_eq(idx.get_names(), ["chr1", "chr2"])
  assert_eq(idx.contains("chr1"), true)
  assert_eq(idx.contains("chr3"), false)
}

///|
test "faidx_parse_fai" {
  let fai = "chr1\t12\t6\t6\t7\nchr2\t8\t20\t8\t9\n"
  try {
    let parsed = @src.parse_fai(fai)
    assert_eq(parsed.num_sequences(), 2)
    let rec1 = parsed.get_record("chr1")
    match rec1 {
      Some(r) => {
        assert_eq(r.length, 12)
        assert_eq(r.offset, 6)
        assert_eq(r.line_len, 6)
        assert_eq(r.line_blen, 7)
      }
      None => assert_eq(true, false)
    }
  } catch {
    _ => assert_eq(true, false)
  }
}

///|
test "faidx_write_fai" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let idx = @src.build_fai(fasta)
  let fai_str = @src.write_fai(idx)
  try {
    let parsed = @src.parse_fai(fai_str)
    assert_eq(parsed.num_sequences(), 2)
    assert_eq(parsed.get_names(), ["chr1", "chr2"])
  } catch {
    _ => assert_eq(true, false)
  }
}

///|
test "faidx_fetch_full_sequence" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let fa = @src.Fasta::from_content(fasta)
  let seq = fa.get_seq("chr1")
  match seq {
    Some(s) => assert_eq(s.to_string(), "ACGTACGTACGT")
    None => assert_eq(true, false)
  }
  let seq2 = fa.get_seq("chr2")
  match seq2 {
    Some(s) => assert_eq(s.to_string(), "TTTTAAAA")
    None => assert_eq(true, false)
  }
}

///|
test "faidx_fetch_subsequence" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let fa = @src.Fasta::from_content(fasta)
  try {
    let sub = fa.fetch("chr1", 0, 4)
    match sub {
      Some(seq) => assert_eq(seq.to_string(), "ACGT")
      None => assert_eq(true, false)
    }
    let sub2 = fa.fetch("chr1", 4, 8)
    match sub2 {
      Some(seq) => assert_eq(seq.to_string(), "ACGT")
      None => assert_eq(true, false)
    }
  } catch {
    _ => assert_eq(true, false)
  }
}

///|
test "faidx_fetch_across_lines" {
  let fasta = ">chr1\nACGTAC\nGTACGT\n>chr2\nTTTTAAAA\n"
  let fa = @src.Fasta::from_content(fasta)
  try {
    let sub = fa.fetch("chr1", 3, 9)
    match sub {
      Some(seq) => assert_eq(seq.to_string(), "TACGTA")
      None => assert_eq(true, false)
    }
  } catch {
    _ => assert_eq(true, false)
  }
}

///|
test "faidx_get_length" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let fa = @src.Fasta::from_content(fasta)
  assert_eq(fa.get_length("chr1"), Some(12))
  assert_eq(fa.get_length("chr2"), Some(8))
  assert_eq(fa.get_length("chr3"), None)
}

///|
test "faidx_contains" {
  let fasta = build_fasta("chr1", "ACGTACGTACGT", "chr2", "TTTTAAAA")
  let fa = @src.Fasta::from_content(fasta)
  assert_eq(fa.contains("chr1"), true)
  assert_eq(fa.contains("chr2"), true)
  assert_eq(fa.contains("chr3"), false)
}

///|
test "faidx_multiline_fasta" {
  let fasta = ">chr1\nACGT\nACGT\nACGT\n>chr2\nTT\nTTAA\nAA\n"
  let fa = @src.Fasta::from_content(fasta)
  let seq = fa.get_seq("chr1")
  match seq {
    Some(s) => {
      assert_eq(s.length(), 12)
      assert_eq(s.to_string(), "ACGTACGTACGT")
    }
    None => assert_eq(true, false)
  }
  let seq2 = fa.get_seq("chr2")
  match seq2 {
    Some(s) => {
      assert_eq(s.length(), 8)
      assert_eq(s.to_string(), "TTTTAAAA")
    }
    None => assert_eq(true, false)
  }
}

///|
fn build_fasta(
  name1 : String,
  seq1 : String,
  name2 : String,
  seq2 : String,
) -> String {
  ">" + name1 + "\n" + seq1 + "\n>" + name2 + "\n" + seq2 + "\n"
}
