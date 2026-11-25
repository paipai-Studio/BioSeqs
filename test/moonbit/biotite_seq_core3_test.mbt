///|
/// Tests for biotite_seq_core3.mbt — alphabets, profiles, annotations.
test "biotite3_alphabet_mapper" {
  let mapper = @src.biotite3_alphabet_mapper(["A", "C", "G", "T"], [
    "A", "C", "G", "U",
  ])
  let mapped = @src.biotite3_alphabet_mapper_map(mapper, [0, 1, 2, 3])
  assert_eq(mapped.length(), 4)
  assert_eq(mapped[0], 0)
  assert_eq(mapped[3], -1)
  let symbols = @src.biotite3_alphabet_mapper_map_symbols(mapper, ["A", "T"])
  assert_eq(symbols[0], "A")
  assert_eq(symbols[1], "")
}

///|
test "biotite3_common_alphabet_preserves_order" {
  let common = @src.biotite3_common_alphabet([["A", "C"], ["C", "G"]])
  assert_eq(common.length(), 3)
  assert_eq(common[0], "A")
  assert_eq(common[1], "C")
  assert_eq(common[2], "G")
}

///|
test "biotite3_sequence_profile" {
  let profile = @src.biotite3_sequence_profile(["A", "C"], [[0, 0], [0, 1]])
  assert_eq(profile.frequencies.length(), 2)
  let consensus = @src.biotite3_profile_consensus(profile)
  assert_eq(consensus[0], 0)
  let entropy = @src.biotite3_profile_entropy(profile)
  assert_true(@moonbit.assert_num_close(entropy[0], 0.0, 1.0e-12))
  assert_true(entropy[1] > 0.0)
}

///|
test "biotite3_positional_sequence" {
  let seq = @src.biotite3_pure_positional_sequence(3)
  let (symbol, position) = @src.biotite3_positional_sequence_get(seq, 1)
  assert_eq(symbol, "1")
  assert_true(@moonbit.assert_num_close(position, 1.0, 1.0e-12))
  let sliced = @src.biotite3_positional_sequence_slice(seq, 1, 3)
  assert_eq(sliced.symbols.length(), 2)
  assert_eq(sliced.symbols[0], "1")
}

///|
test "biotite3_location_helpers" {
  let outer = @src.biotite3_location(
    1,
    10,
    1,
    @src.biotite3_location_defect_none(),
  )
  let inner = @src.biotite3_location(
    3,
    5,
    1,
    @src.biotite3_location_defect_none(),
  )
  let disjoint = @src.biotite3_location(20, 30, 1, 0)
  assert_true(@src.biotite3_location_overlap(outer, inner))
  assert_true(!@src.biotite3_location_overlap(outer, disjoint))
  assert_true(@src.biotite3_location_contains(outer, inner))
  assert_eq(@src.biotite3_location_to_string(outer), "1..10(+)")
}

///|
test "biotite3_annotation_merge_split" {
  let annot = @src.biotite3_annotation_add(
    @src.biotite3_annotation(),
    @src.biotite3_feature("gene", @src.biotite3_location(1, 5, 1, 0), []),
  )
  let merged = @src.biotite3_annotation_merge(annot, "gene")
  assert_eq(merged.features.length(), 1)
  let (left, right) = @src.biotite3_annotation_split(merged, 3)
  assert_eq(left.features.length(), 1)
  assert_eq(right.features.length(), 1)
}

///|
test "biotite3_annotated_sequence" {
  let annot = @src.biotite3_annotation_add(
    @src.biotite3_annotation(),
    @src.biotite3_feature("gene", @src.biotite3_location(0, 4, 1, 0), []),
  )
  let seq = @src.biotite3_annotated_sequence("AACCGGTT", annot)
  assert_eq(
    @src.biotite3_annotated_sequence_get_features(seq, "gene").length(),
    1,
  )
  let sub = @src.biotite3_annotated_sequence_subsequence(seq, 0, 4)
  assert_eq(sub.sequence, "AACC")
}
