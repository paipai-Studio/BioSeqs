///|
/// Tests for the pyfaidx Sequence slice model, cross-validated against
/// pyfaidx 0.9.0.4 `Sequence.__getitem__` and Python `slice.indices`.

///|
fn idx_is(start : Int?, stop : Int?, step : Int?, ea : Int, eb : Int, ec : Int) -> Bool {
  try {
    let (a, b, c) = @src.pyfaidx_slice_indices(start, stop, step, 6)
    a == ea && b == eb && c == ec
  } catch {
    PyfaidxError(_) => false
  }
}

///|
fn seq_is(
  s : @src.PyfaidxSequence,
  seq : String,
  start : Int?,
  end : Int?,
  comp : Bool,
) -> Bool {
  s.seq == seq && s.start == start && s.end == end && s.comp == comp
}

///|
test "pyfaidx_slice_indices_matches_python" {
  assert_true(idx_is(None, None, None, 0, 6, 1))
  assert_true(idx_is(None, None, Some(-1), 5, -1, -1))
  assert_true(idx_is(None, None, Some(2), 0, 6, 2))
  assert_true(idx_is(None, None, Some(-3), 5, -1, -3))
  assert_true(idx_is(Some(1), Some(-1), Some(1), 1, 5, 1))
  assert_true(idx_is(Some(-100), Some(100), Some(1), 0, 6, 1))
  assert_true(idx_is(None, Some(0), Some(-1), 5, 0, -1))
  assert_true(idx_is(Some(2), None, Some(-1), 2, -1, -1))
  assert_true(idx_is(Some(-2), None, Some(1), 4, 6, 1))
  assert_true(idx_is(Some(0), Some(0), Some(1), 0, 0, 1))
}

///|
test "pyfaidx_sequence_slice_1based" {
  let x = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(1), Some(6), false)
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, None, Some(3), None), "ATC", Some(1), Some(3), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, Some(3), None, None), "GTA", Some(4), Some(6), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, Some(1), Some(-1), None), "TCGT", Some(2), Some(5), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, None, None, Some(-1)), "ATGCTA", Some(6), Some(1), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, None, None, Some(2)), "ACT", None, None, false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(x, None, None, Some(-3)), "AC", None, None, false))
}

///|
test "pyfaidx_sequence_slice_0based" {
  let y = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(0), Some(6), false)
  assert_true(seq_is(@src.pyfaidx_seq_slice(y, None, Some(3), None), "ATC", Some(0), Some(3), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(y, Some(3), None, None), "GTA", Some(3), Some(6), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(y, Some(1), Some(-1), None), "TCGT", Some(1), Some(5), false))
  assert_true(seq_is(@src.pyfaidx_seq_slice(y, None, None, Some(-1)), "ATGCTA", Some(6), Some(0), false))
}

///|
test "pyfaidx_sequence_getitem_int" {
  let x = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(1), Some(6), false)
  assert_true(seq_is(@src.pyfaidx_seq_getitem(x, 0), "A", Some(1), Some(1), false))
  assert_true(seq_is(@src.pyfaidx_seq_getitem(x, -2), "T", Some(5), Some(5), false))
}

///|
test "pyfaidx_sequence_neg_complement_reverse" {
  let x = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(1), Some(6), false)
  // reverse complement
  assert_true(seq_is(@src.pyfaidx_seq_neg(x), "TACGAT", Some(6), Some(1), true))
  // complement
  assert_true(seq_is(@src.pyfaidx_seq_complement(x), "TAGCAT", Some(1), Some(6), true))
  // reverse
  assert_true(seq_is(@src.pyfaidx_seq_reverse(x), "ATGCTA", Some(6), Some(1), false))
  // fancy name
  assert_eq(@src.pyfaidx_seq_fancy_name(@src.pyfaidx_seq_neg(x)), "chr1:6-1 (complement)")
}

///|
test "pyfaidx_sequence_orientation" {
  let x = @src.pyfaidx_seq_new("chr1", "ATCGTA", Some(1), Some(6), false)
  assert_eq(@src.pyfaidx_seq_orientation(x), Some(1))
  assert_eq(@src.pyfaidx_seq_orientation(@src.pyfaidx_seq_neg(x)), Some(-1))
  assert_eq(@src.pyfaidx_seq_orientation(@src.pyfaidx_seq_reverse(x)), None)
  assert_eq(@src.pyfaidx_seq_orientation(@src.pyfaidx_seq_complement(x)), None)
}