///|
/// Tests for biotite_seq_profile.mbt — position-specific scoring matrices.
test "biotite_seq_profile_build" {
  let profile = @src.BiotiteSequenceProfile::from_alignment(["ACGT", "ACGT"], [
    "A", "C", "G", "T",
  ])
  assert_eq(profile.length(), 4)
  assert_eq(profile.consensus_sequence(), "ACGT")
  assert_eq(profile.position_scores(0).length(), 4)
  assert_eq(profile.n_sequences, 2)
}

///|
test "biotite_seq_profile_score" {
  let profile = @src.BiotiteSequenceProfile::from_alignment(["ACGT", "ACGT"], [
    "A", "C", "G", "T",
  ])
  assert_true(profile.score("ACGT") > profile.score("TTTT"))
}

///|
test "biotite_seq_profile_gap_ignored" {
  let profile = @src.BiotiteSequenceProfile::from_alignment(["ACGT", "----"], [
    "A", "C", "G", "T",
  ])
  assert_true(profile.is_ignored(0, [0.6]))
  assert_true(!profile.is_ignored(0, [0.1]))
  assert_true(profile.is_ignored(9, [0.0]))
}

///|
test "biotite_seq_profile_empty" {
  let profile = @src.BiotiteSequenceProfile::from_alignment([], [
    "A", "C", "G", "T",
  ])
  assert_eq(profile.length(), 0)
  assert_eq(profile.n_sequences, 0)
}
