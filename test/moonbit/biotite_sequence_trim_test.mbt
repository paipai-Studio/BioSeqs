///|
/// Tests for biotite_sequence_trim.mbt.
test "trim single strand left end" {
  let (start, stop) = @src.biotite_trim_single_strand("ACGTACGTAA", "ACGTAC")
  assert_eq(start, 6)
  assert_eq(stop, 10)
}

///|
test "trim single strand right end" {
  let (start, stop) = @src.biotite_trim_single_strand("TTACGTAC", "ACGTAC")
  assert_eq(start, 0)
  assert_eq(stop, 2)
}

///|
test "trim single strand both ends prefer left on tie" {
  // both ends match "ACGT" fully; the left end wins on equal score
  let (start, stop) = @src.biotite_trim_single_strand("ACGTACGT", "ACGT")
  assert_eq(start, 4)
  assert_eq(stop, 8)
}

///|
test "trim single strand no match keeps sequence" {
  let (start, stop) = @src.biotite_trim_single_strand("ACGTACGT", "GGGG")
  assert_eq(start, 0)
  assert_eq(stop, 8)
}

///|
test "trim single strand mismatch tolerance" {
  // right end: the 8-base overlap carries one mismatch (allowed with
  // max_mismatch=1) and outscores any shorter overlap on this sequence
  let (start, stop) = @src.biotite_trim_single_strand(
    "TTACGTACGT",
    "ACGTTCGT",
    max_mismatch=1,
  )
  assert_eq(start, 0)
  assert_eq(stop, 2)
}

///|
test "trim single strand zero mismatch blocks mismatched overlap" {
  // the o=8 overlap has one mismatch -> rejected, o=4 remains valid
  let (start, stop) = @src.biotite_trim_single_strand(
    "ACGTACGT",
    "ACGTTCGT",
    max_mismatch=0,
  )
  assert_eq(start, 0)
  assert_eq(stop, 4)
}

///|
test "trim single strand min overlap" {
  // primer shares only 2 bases with the left end; min_overlap=3 forbids
  // the tiny overlap
  let (start, stop) = @src.biotite_trim_single_strand(
    "ACAAAAAAAA",
    "AC",
    min_overlap=3,
  )
  assert_eq(start, 0)
  assert_eq(stop, 10)
}

///|
test "trim double strand independent primers" {
  let (t1, t2) = @src.biotite_trim_double_strand(
    "TTTTACGTACGTAAAA", "TTTT", "AAAA",
  )
  // primer_1 matches the left end (4 T), primer_2 the right end (4 A)
  assert_eq(t1.0, 4)
  assert_eq(t1.1, 16)
  assert_eq(t2.0, 0)
  assert_eq(t2.1, 12)
}

///|
test "trim distance weight penalizes mismatches" {
  // heavy weight: the mismatched 6-base overlap scores 5 - 1.5 = 3.5
  // and loses against the perfect 4-base overlap (score 4.0)
  let (start, stop) = @src.biotite_trim_single_strand(
    "ACGTACGT",
    "ACGTTC",
    max_mismatch=1,
    distance_weight=1.5,
  )
  assert_eq(start, 0)
  assert_eq(stop, 4)
}

///|
test "trim light weight keeps longer overlap" {
  // light weight: the mismatched 6-base overlap scores 5 - 0.5 = 4.5
  // and beats the perfect 4-base overlap (score 4.0)
  let (start, stop) = @src.biotite_trim_single_strand(
    "ACGTACGT",
    "ACGTTC",
    max_mismatch=1,
    distance_weight=0.5,
  )
  assert_eq(start, 6)
  assert_eq(stop, 8)
}
