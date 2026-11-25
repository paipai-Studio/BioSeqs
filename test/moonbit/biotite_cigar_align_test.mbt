///|
/// Tests for biotite_cigar_align.mbt — CIGAR parsing and lengths.
test "biotite_cigar_lengths" {
  let cigar = @src.BiotiteCigar::from_string("2M1I3M1D")
  assert_eq(cigar.reference_length(), 6)
  assert_eq(cigar.read_length(), 6)
  assert_eq(cigar.alignment_length(), 7)
  assert_eq(cigar.to_string(), "2M1I3M1D")
}

///|
test "biotite_cigar_from_alignment" {
  let cigar = @src.biotite_alignment_to_cigar("ACGT", "AC-T")
  assert_eq(cigar.reference_length(), 4)
  assert_eq(cigar.read_length(), 3)
  assert_eq(cigar.alignment_length(), 4)
  let (ref_aln, read_aln) = cigar.to_alignment("ACGT", "ACT")
  assert_eq(ref_aln, "ACGT")
  assert_eq(read_aln, "AC-T")
}

///|
test "biotite_cigar_soft_and_hard_clips" {
  let cigar = @src.BiotiteCigar::from_string("5H10M2S")
  assert_eq(cigar.reference_length(), 10)
  assert_eq(cigar.read_length(), 12)
}

///|
test "biotite_cigar_empty" {
  let cigar = @src.BiotiteCigar::from_string("")
  assert_eq(cigar.reference_length(), 0)
  assert_eq(cigar.read_length(), 0)
  assert_eq(cigar.to_string(), "")
}
