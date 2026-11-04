///|
/// Tests for Burrows-Wheeler Transform and FM-index.

// ---------------------------------------------------------------------------
// BWT tests
// ---------------------------------------------------------------------------

///|
test "bwt_simple" {
  // "banana$" → BWT is "annb$aa"  (classic example)
  let result = @src.bwt("banana")
  // Expected BWT: last column of sorted rotations
  // Rotations sorted:
  //   $banana   → last char 'a'
  //   a$banan   → last char 'n'
  //   ana$ban   → last char 'n'
  //   anana$b   → last char 'b'
  //   banana$   → last char '$'
  //   na$bana   → last char 'a'
  //   nana$ba   → last char 'a'
  // So BWT = "annb$aa"
  assert_eq(result, "annb$aa")
}

///|
test "bwt_single_char" {
  let result = @src.bwt("A")
  // "$A" → sorted: $A, A$ → BWT = "A$"
  assert_eq(result, "A$")
}

///|
test "bwt_empty" {
  let result = @src.bwt("")
  assert_eq(result, "$")
}

///|
test "bwt_all_same" {
  let result = @src.bwt("AAA")
  // "$AAA" < "A$AA" < "AA$A" < "AAA$"
  // BWT = last chars: A, A, A, $ = "AAA$"
  assert_eq(result, "AAA$")
}

///|
test "bwt_inverse_simple" {
  // BWT("banana") = "annb$aa", inverse should give "banana"
  let bwt_str = @src.bwt("banana")
  let original = @src.bwt_inverse(bwt_str)
  assert_eq(original, "banana")
}

///|
test "bwt_inverse_single" {
  let bwt_str = @src.bwt("A")
  let original = @src.bwt_inverse(bwt_str)
  assert_eq(original, "A")
}

///|
test "bwt_inverse_all_same" {
  let bwt_str = @src.bwt("GGG")
  let original = @src.bwt_inverse(bwt_str)
  assert_eq(original, "GGG")
}

///|
test "bwt_inverse_dna" {
  let seq = "ACGTACGT"
  let bwt_str = @src.bwt(seq)
  let original = @src.bwt_inverse(bwt_str)
  assert_eq(original, seq)
}

///|
test "bwt_inverse_empty" {
  let bwt_str = @src.bwt("")
  let original = @src.bwt_inverse(bwt_str)
  assert_eq(original, "")
}

///|
test "bwt_length_preserved" {
  // BWT string length = original length + 1 (for sentinel)
  let s = "ACGTACGTACGT"
  let bwt_str = @src.bwt(s)
  assert_eq(bwt_str.length(), s.length() + 1)
}

// ---------------------------------------------------------------------------
// Suffix Array tests
// ---------------------------------------------------------------------------

///|
test "suffix_array_basic" {
  let sa = @src.suffix_array("ab")
  // "ab$" → sorted suffixes: "$ab" (2), "ab$" (0), "b$a" (1)
  assert_eq(sa.length(), 3)
  assert_eq(sa[0], 2) // $ab
  assert_eq(sa[1], 0) // ab$
  assert_eq(sa[2], 1) // b$a
}

///|
test "suffix_array_length" {
  let sa = @src.suffix_array("ACGT")
  assert_eq(sa.length(), 5) // 4 + sentinel
}

// ---------------------------------------------------------------------------
// FM-index tests
// ---------------------------------------------------------------------------

///|
test "fm_index_bwt" {
  let fm = @src.FMIndex::new("banana")
  assert_eq(fm.bwt(), "annb$aa")
}

///|
test "fm_index_len" {
  let fm = @src.FMIndex::new("ACGT")
  assert_eq(fm.len(), 5) // 4 + sentinel
}

///|
test "fm_count_exact_match" {
  let fm = @src.FMIndex::new("banana")
  // "ana" appears in positions 1, 3 → count = 2
  let count = fm.count("ana")
  assert_eq(count, 2)
}

///|
test "fm_count_single_char" {
  let fm = @src.FMIndex::new("banana")
  // 'a' appears 3 times
  assert_eq(fm.count("a"), 3)
  // 'n' appears 2 times
  assert_eq(fm.count("n"), 2)
  // 'b' appears 1 time
  assert_eq(fm.count("b"), 1)
}

///|
test "fm_count_no_match" {
  let fm = @src.FMIndex::new("banana")
  assert_eq(fm.count("xyz"), 0)
  assert_eq(fm.count("bnn"), 0)
}

///|
test "fm_count_full_string" {
  let fm = @src.FMIndex::new("banana")
  assert_eq(fm.count("banana"), 1)
}

///|
test "fm_locate_basic" {
  let fm = @src.FMIndex::new("banana")
  let positions = fm.locate("ana")
  // "ana" occurs at positions 1 and 3 (0-based) in "banana"
  assert_eq(positions.length(), 2)
  assert_eq(positions[0], 1)
  assert_eq(positions[1], 3)
}

///|
test "fm_locate_single" {
  let fm = @src.FMIndex::new("banana")
  let positions = fm.locate("banana")
  assert_eq(positions.length(), 1)
  assert_eq(positions[0], 0)
}

///|
test "fm_locate_no_match" {
  let fm = @src.FMIndex::new("banana")
  let positions = fm.locate("xyz")
  assert_eq(positions.length(), 0)
}

///|
test "fm_count_dna_sequence" {
  let fm = @src.FMIndex::new("ACGTACGTACGT")
  // "ACGT" appears 3 times in "ACGTACGTACGT"
  assert_eq(fm.count("ACGT"), 3)
  // "CGT" appears 3 times
  assert_eq(fm.count("CGT"), 3)
  // "TT" appears 0 times
  assert_eq(fm.count("TT"), 0)
}

///|
test "fm_locate_dna_sequence" {
  let fm = @src.FMIndex::new("ACGTACGTACGT")
  let positions = fm.locate("ACGT")
  assert_eq(positions.length(), 3)
  assert_eq(positions[0], 0)
  assert_eq(positions[1], 4)
  assert_eq(positions[2], 8)
}

///|
test "fm_lf_mapping" {
  // LF(i) maps last column position i to first column position
  // For "banana" BWT = "annb$aa":
  // First column (sorted): $, a, a, a, b, n, n
  // LF mapping should be consistent: applying LF repeatedly should cycle
  let fm = @src.FMIndex::new("banana")
  // Start at any position, apply LF n times, should come back (kind of)
  // Actually LF is a permutation, applying n times covers all positions
  let start = 2
  let lf1 = fm.lf(start)
  assert_true(lf1 >= 0)
  assert_true(lf1 < fm.len())
}

///|
test "fm_count_empty_pattern" {
  let fm = @src.FMIndex::new("ACGT")
  // Empty pattern matches all suffixes (the whole text)
  // Actually backward search with empty pattern returns [0, n-1] = n matches
  assert_eq(fm.count(""), 5) // includes sentinel
}

///|
test "fm_locate_single_char" {
  let fm = @src.FMIndex::new("ACGTACGT")
  let positions = fm.locate("C")
  assert_eq(positions.length(), 2)
}

///|
test "fm_count_overlapping" {
  let fm = @src.FMIndex::new("AAAA")
  // "AA" appears 3 times (overlapping: positions 0,1,2)
  assert_eq(fm.count("AA"), 3)
  // "AAAA" appears 1 time
  assert_eq(fm.count("AAAA"), 1)
}

///|
test "fm_locate_overlapping" {
  let fm = @src.FMIndex::new("AAAA")
  let positions = fm.locate("AA")
  assert_eq(positions.length(), 3)
  assert_eq(positions[0], 0)
  assert_eq(positions[1], 1)
  assert_eq(positions[2], 2)
}
