///|
test "bwa_suffix_array basic" {
  let sa = @src.bwa_suffix_array("banana")
  assert_eq(sa.length(), 6)
  assert_eq(sa[0], 5)
  assert_eq(sa[1], 3)
  assert_eq(sa[2], 1)
  assert_eq(sa[3], 0)
  assert_eq(sa[4], 4)
  assert_eq(sa[5], 2)
}

///|
test "bwa_suffix_array single char" {
  let sa = @src.bwa_suffix_array("a")
  assert_eq(sa.length(), 1)
  assert_eq(sa[0], 0)
}

///|
test "bwa_bwt basic" {
  let text = "banana"
  let sa = @src.bwa_suffix_array(text)
  let bwt = @src.bwa_bwt(text, sa)
  assert_eq(bwt.length(), text.length())
}

///|
test "bwa_bwt_from_text" {
  let bwt = @src.bwa_bwt_from_text("banana")
  assert_eq(bwt.length(), 6)
}

///|
test "fm_index_build basic" {
  let index = @src.fm_index_build("banana")
  assert_true(index.bwt.length() > 0)
  assert_true(index.sa.length() > 0)
}

///|
test "fm_backward_search exact match" {
  let index = @src.fm_index_build("banana")
  let (lo, hi) = @src.fm_backward_search(index, "ana")
  assert_eq(hi - lo, 2)
}

///|
test "fm_backward_search no match" {
  let index = @src.fm_index_build("banana")
  let (lo, hi) = @src.fm_backward_search(index, "xyz")
  assert_eq(hi - lo, 0)
}

///|
test "fm_backward_search single char" {
  let index = @src.fm_index_build("banana")
  let (lo, hi) = @src.fm_backward_search(index, "a")
  assert_eq(hi - lo, 3)
  let (lo2, hi2) = @src.fm_backward_search(index, "n")
  assert_eq(hi2 - lo2, 2)
}

///|
test "fm_locate finds positions" {
  let index = @src.fm_index_build("banana")
  let positions = @src.fm_locate(index, "ana")
  assert_eq(positions.length(), 2)
  let mut has_1 = false
  let mut has_3 = false
  for p in positions {
    if p == 1 { has_1 = true }
    if p == 3 { has_3 = true }
  }
  assert_true(has_1)
  assert_true(has_3)
}

///|
test "fm_locate whole string" {
  let index = @src.fm_index_build("banana")
  let positions = @src.fm_locate(index, "banana")
  assert_eq(positions.length(), 1)
  assert_eq(positions[0], 0)
}

///|
test "bwa_index_build" {
  let index = @src.bwa_index_build("ref1", "ACGTACGTACGTACGT")
  assert_true(index.fm.bwt.length() > 0)
  assert_eq(index.ref_name, "ref1")
}

///|
test "bwa_align exact match" {
  let index = @src.bwa_index_build("ref1", "ACGTACGTACGTACGT")
  let alignments = @src.bwa_align(index, "ACGTACGT")
  assert_true(alignments.length() > 0)
  let aln = alignments[0]
  assert_eq(aln.ref_name, "ref1")
  assert_true(aln.pos >= 0)
}

///|
test "bwa_align with mismatch" {
  let index = @src.bwa_index_build("ref1", "ACGTACGTACGTACGT")
  let alignments = @src.bwa_align(index, "ACGTACCTACGT")
  assert_true(alignments.length() > 0)
}

///|
test "bwa_to_sam produces valid SAM" {
  let index = @src.bwa_index_build("ref1", "ACGTACGTACGTACGT")
  let alignments = @src.bwa_align(index, "ACGTACGT")
  assert_true(alignments.length() > 0)
  let sam = @src.bwa_to_sam(alignments[0], "read1", "ACGTACGT")
  assert_true(sam.length() > 0)
  assert_true(sam.contains("\t"))
}

///|
test "bwa_align_batch multiple reads" {
  let index = @src.bwa_index_build("ref1", "ACGTACGTACGTACGT")
  let reads = [("r1", "ACGTACGT"), ("r2", "ACGTACGT"), ("r3", "TTACGTAC")]
  let sam = @src.bwa_align_batch(index, reads)
  assert_true(sam.length() > 0)
}

///|
test "fm_backward_search DNA sequence" {
  let index = @src.fm_index_build("ACGTACGTACGT")
  let (lo, hi) = @src.fm_backward_search(index, "ACGT")
  assert_eq(hi - lo, 3)
  let (lo2, hi2) = @src.fm_backward_search(index, "CGTA")
  assert_eq(hi2 - lo2, 2)
}
