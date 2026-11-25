///|
/// Tests for biotite_sequence_bundles.mbt.
test "biotite_bundles_sequence_add_get" {
  let bundles = @src.biotite_bundles_new()
  assert_eq(@src.biotite_bundles_length(bundles), 0)
  let one = @src.biotite_bundles_add_sequence(bundles, "seq1", "ACGT")
  assert_eq(@src.biotite_bundles_length(one), 1)
  match @src.biotite_bundles_get_sequence(one, "seq1") {
    Some(seq) => assert_eq(seq, "ACGT")
    None => fail("expected seq1")
  }
  assert_true(@src.biotite_bundles_get_sequence(one, "missing") is None)
}

///|
test "biotite_bundles_sequence_replace" {
  let bundles = @src.biotite_bundles_add_sequence(
    @src.biotite_bundles_new(),
    "seq1",
    "ACGT",
  )
  let replaced = @src.biotite_bundles_add_sequence(bundles, "seq1", "TTTT")
  assert_eq(@src.biotite_bundles_length(replaced), 1)
  match @src.biotite_bundles_get_sequence(replaced, "seq1") {
    Some(seq) => assert_eq(seq, "TTTT")
    None => fail("expected seq1")
  }
}

///|
test "biotite_bundles_structure_add_get" {
  let arr = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      element="C",
    ),
  ])
  let bundles = @src.biotite_bundles_add_structure(
    @src.biotite_bundles_new(),
    "model",
    arr,
  )
  assert_eq(@src.biotite_bundles_structure_length(bundles), 1)
  match @src.biotite_bundles_get_structure(bundles, "model") {
    Some(structure) => assert_eq(structure.len(), 1)
    None => fail("expected model")
  }
}
