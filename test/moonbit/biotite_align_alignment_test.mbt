///|
/// Tests for biotite_align_alignment.mbt.
test "biotite_alignment_from_strings" {
  let alignment = @src.biotite_alignment_from_strings(["ACGT", "AC-T"], [
    "a", "b",
  ])
  assert_eq(@src.biotite_alignment_num_sequences(alignment), 2)
  assert_eq(@src.biotite_alignment_length(alignment), 4)
  assert_true(@moonbit.assert_num_close(alignment.score, 0.75, 1.0e-12))
  assert_true(@src.biotite_alignment_is_gapped(alignment))
}

///|
test "biotite_alignment_getters" {
  let alignment = @src.biotite_alignment_from_strings(["ACGT", "ACGT"], [
    "a", "b",
  ])
  match @src.biotite_alignment_get_sequence(alignment, 0) {
    Some(seq) => assert_eq(seq, "ACGT")
    None => fail("expected sequence")
  }
  match @src.biotite_alignment_get_name(alignment, 1) {
    Some(name) => assert_eq(name, "b")
    None => fail("expected name")
  }
  assert_true(@src.biotite_alignment_get_sequence(alignment, 9) is None)
  assert_true(!@src.biotite_alignment_is_gapped(alignment))
}

///|
test "biotite_alignment_identity" {
  let identical = @src.biotite_alignment_from_strings(["ACGT", "ACGT"], [
    "a", "b",
  ])
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_alignment_identity(identical),
      1.0,
      1.0e-12,
    ),
  )
  let half = @src.biotite_alignment_from_strings(["ACGT", "AGGT"], ["a", "b"])
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_alignment_identity(half),
      0.75,
      1.0e-12,
    ),
  )
}

///|
test "biotite_alignment_empty" {
  let alignment = @src.biotite_alignment_from_strings([], [])
  assert_eq(@src.biotite_alignment_length(alignment), 0)
  assert_eq(@src.biotite_alignment_identity(alignment), 1.0)
}
