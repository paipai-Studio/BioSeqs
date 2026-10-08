test "banded_nw_identical" {
  // Identical sequences — band of 0 should suffice
  let r = @src.banded_needleman_wunsch("ACGTACGT", "ACGTACGT", 2, 1.0, -1.0, -2.0)
  assert_eq(r.aligned_seq1, "ACGTACGT")
  assert_eq(r.aligned_seq2, "ACGTACGT")
  assert_eq(r.score, 8.0)
}

test "banded_nw_close_homolog" {
  // One mismatch in the middle — small band is enough
  let r = @src.banded_needleman_wunsch("ACGTACGT", "ACGAACGT", 3, 1.0, -1.0, -2.0)
  assert_eq(r.aligned_seq1, "ACGTACGT")
  assert_eq(r.aligned_seq2, "ACGAACGT")
  // 3 matches + 1 mismatch + 4 matches = 7
  assert_eq(r.score, 6.0)
}

test "banded_nw_with_gaps" {
  // Small insertion — band of 3 covers it
  let r = @src.banded_needleman_wunsch("ACGT", "ACAGT", 3, 2.0, -1.0, -2.0)
  // Should align as ACG-T / ACAGT
  assert_eq(r.aligned_seq1.length(), r.aligned_seq2.length())
  assert_eq(r.score, 6.0) // 4 matches * 2 - 1 gap * 2 = 6
}

test "banded_nw_widens_band" {
  // Length difference > k → band auto-widens
  let r = @src.banded_needleman_wunsch("ACGT", "ACGTACGT", 1, 1.0, -1.0, -2.0)
  // Should still produce a valid alignment
  assert_true(r.aligned_seq1.length() == r.aligned_seq2.length())
  assert_true(r.band_used >= 4) // |8-4| = 4, so band >= 4+1 = 5
}

test "banded_nw_matches_full" {
  // Banded result should match full NW when band is wide enough
  let seq1 = "THEQUICKBROWNFOX"
  let seq2 = "THEQUIKCBROWNFOX"
  let banded = @src.banded_needleman_wunsch(seq1, seq2, 5, 1.0, -1.0, -2.0)
  let full = @src.needleman_wunsch(seq1, seq2, @src.NWScoring::new(1.0, -1.0, -2.0))
  assert_eq(banded.score, full.score)
}

test "banded_sw_identical" {
  let r = @src.banded_smith_waterman("ACGTACGT", "ACGTACGT", 2, 2.0, -1.0, -2.0)
  assert_eq(r.aligned_seq1, "ACGTACGT")
  assert_eq(r.aligned_seq2, "ACGTACGT")
  assert_eq(r.score, 16.0)
}

test "banded_sw_partial" {
  // Local alignment of a substring
  let r = @src.banded_smith_waterman("XXXXACGTYYYY", "ACGT", 5, 2.0, -1.0, -2.0)
  assert_eq(r.aligned_seq1, "ACGT")
  assert_eq(r.aligned_seq2, "ACGT")
  assert_eq(r.score, 8.0)
}

test "banded_nw_protein" {
  // Protein-like sequences with a few substitutions
  let seq1 = "MVLSPADKTNVKAAW"
  let seq2 = "MVLSPADKSNVKAAW"
  let r = @src.banded_needleman_wunsch(seq1, seq2, 3, 1.0, -1.0, -2.0)
  // One substitution (T→S), 14 matches + 1 mismatch
  assert_eq(r.score, 13.0)
}
