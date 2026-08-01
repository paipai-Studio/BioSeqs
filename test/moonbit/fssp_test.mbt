///|
/// Tests for Bio.FSSP structure alignment database parser.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "fssp_header_creation" {
  let h = @src.FsspHeader::new(
    "1dfa_A", "30-jul-1998", "crystal structure",
    "17beta-hsd", "human", "holm",
    25, 3, 2.0,
  )
  assert_eq(h.pdbid(), "1dfa_A")
  assert_eq(h.pdbid_date(), "30-jul-1998")
  assert_eq(h.title(), "crystal structure")
  assert_eq(h.compound(), "17beta-hsd")
  assert_eq(h.source(), "human")
  assert_eq(h.author(), "holm")
  assert_eq(h.seq_length(), 25)
  assert_eq(h.n_alignments(), 3)
  assert_eq(h.threshold(), "2")
}

test "fssp_alignment_creation" {
  let a = @src.FsspAlignment::new(
    "1csa_A", "1csa_A", 1, "MNIFVHEKDLFRTIVS",
    12.1, 2.0, 16, 68.8,
  )
  assert_eq(a.pdbid(), "1csa_A")
  assert_eq(a.alignment_id(), "1csa_A")
  assert_eq(a.start_position(), 1)
  assert_eq(a.aligned_sequence(), "MNIFVHEKDLFRTIVS")
  assert_eq(a.z_score(), 12.1)
  assert_eq(a.rmsd(), 2.0)
  assert_eq(a.lali(), 16)
  assert_eq(a.pid(), 68.8)
}

test "fssp_data_creation" {
  let h = @src.FsspHeader::new(
    "1dfa_A", "30-jul-1998", "", "", "", "",
    25, 2, 2.0,
  )
  let a1 = @src.FsspAlignment::new(
    "1dfa_A", "1dfa_A", 1, "ACGT", 15.0, 1.5, 4, 100.0,
  )
  let a2 = @src.FsspAlignment::new(
    "1csa_A", "1csa_A", 1, "ACGT", 10.0, 2.5, 4, 75.0,
  )
  let data = @src.FsspData::new(h, [a1, a2])
  assert_eq(data.n_alignments(), 2)
  assert_true(data.reference().unwrap().pdbid() == "1dfa_A")
}

test "fssp_filter_by_zscore" {
  let h = @src.FsspHeader::new(
    "1dfa_A", "", "", "", "", "", 10, 3, 2.0,
  )
  let a1 = @src.FsspAlignment::new("a1", "a1", 1, "AAA", 15.0, 1.0, 3, 100.0)
  let a2 = @src.FsspAlignment::new("a2", "a2", 1, "AAA", 10.0, 2.0, 3, 80.0)
  let a3 = @src.FsspAlignment::new("a3", "a3", 1, "AAA", 5.0, 3.0, 3, 60.0)
  let data = @src.FsspData::new(h, [a1, a2, a3])
  let filtered = data.filter_by_zscore(10.0)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0].pdbid(), "a1")
  assert_eq(filtered[1].pdbid(), "a2")
}

test "fssp_empty_alignments" {
  let h = @src.FsspHeader::new(
    "1dfa_A", "", "", "", "", "", 0, 0, 0.0,
  )
  let data = @src.FsspData::new(h, [])
  assert_eq(data.n_alignments(), 0)
  assert_true(data.reference() is None)
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "fssp_sample_text_has_header" {
  let text = @src.fssp_sample_text()
  assert_true(text.contains("HEADER"))
  assert_true(text.contains("1dfa_A"))
}

test "fssp_sample_text_has_alignments" {
  let text = @src.fssp_sample_text()
  assert_true(text.contains("## ALIGNMENTS"))
  assert_true(text.contains("1csa_A"))
  assert_true(text.contains("1hsd_A"))
}

test "fssp_sample_text_has_threshold" {
  let text = @src.fssp_sample_text()
  assert_true(text.contains("THRESHOLD"))
  assert_true(text.contains("SEQLENGTH"))
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

test "fssp_parse_header" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let h = data.header()
  assert_eq(h.pdbid(), "1dfa_A")
  assert_eq(h.seq_length(), 25)
  assert_eq(h.n_alignments(), 3)
  assert_eq(h.threshold(), "2")
}

test "fssp_parse_header_title" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let h = data.header()
  assert_true(h.title().contains("17beta-hydroxysteroid"))
  assert_eq(h.compound(), "17beta-hsd")
  assert_eq(h.source(), "human")
  assert_eq(h.author(), "holm")
}

test "fssp_parse_alignments_count" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  assert_eq(data.n_alignments(), 3)
}

test "fssp_parse_alignment_details" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let alignments = data.alignments()
  // First alignment is the reference
  let ref_aln = alignments[0]
  assert_eq(ref_aln.pdbid(), "1dfa_A")
  assert_eq(ref_aln.start_position(), 1)
  assert_eq(ref_aln.aligned_sequence(), "MNIFVHEKDLFRTIVS")
  assert_eq(ref_aln.z_score(), 15.3)
  assert_eq(ref_aln.rmsd(), 1.5)
  assert_eq(ref_aln.lali(), 16)
  assert_eq(ref_aln.pid(), 75.0)
}

test "fssp_parse_second_alignment" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let a2 = data.alignments()[1]
  assert_eq(a2.pdbid(), "1csa_A")
  assert_eq(a2.z_score(), 12.1)
  assert_eq(a2.pid(), 68.8)
}

test "fssp_parse_reference" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let ref = data.reference().unwrap()
  assert_eq(ref.pdbid(), "1dfa_A")
}

test "fssp_parse_empty_input" {
  let data = @src.fssp_parse("")
  assert_eq(data.n_alignments(), 0)
}

test "fssp_parse_header_only" {
  let text = "HEADER \\_1abc_A 1 01-jan-2000\nSEQLENGTH 10\n"
  let data = @src.fssp_parse(text)
  let h = data.header()
  assert_eq(h.pdbid(), "1abc_A")
  assert_eq(h.seq_length(), 10)
  assert_eq(data.n_alignments(), 0)
}

test "fssp_parse_filter_after_parse" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let high_z = data.filter_by_zscore(10.0)
  assert_eq(high_z.length(), 2)
  let low_z = data.filter_by_zscore(20.0)
  assert_eq(low_z.length(), 0)
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

test "fssp_header_to_string" {
  let h = @src.FsspHeader::new(
    "1dfa_A", "30-jul-1998", "", "", "", "", 25, 3, 2.0,
  )
  let s = h.to_string()
  assert_true(s.contains("1dfa_A"))
  assert_true(s.contains("25"))
}

test "fssp_alignment_to_string" {
  let a = @src.FsspAlignment::new(
    "1csa_A", "1csa_A", 1, "ACGT", 12.1, 2.0, 4, 75.0,
  )
  let s = a.to_string()
  assert_true(s.contains("1csa_A"))
  assert_true(s.contains("ACGT"))
  assert_true(s.contains("12.1"))
}

test "fssp_data_to_string" {
  let text = @src.fssp_sample_text()
  let data = @src.fssp_parse(text)
  let s = data.to_string()
  assert_true(s.contains("FSSP"))
  assert_true(s.contains("1dfa_A"))
  assert_true(s.contains("Alignments"))
}
