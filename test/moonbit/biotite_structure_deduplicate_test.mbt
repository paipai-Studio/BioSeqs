///|
/// Tests for biotite_structure_deduplicate.mbt.
fn bsd_atom(name : String, res_id : Int32, x : Double) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname="LIG",
    chainid='A',
    resseq=res_id,
    element="C",
  )
}

///|
test "biotite_deduplicate_structure" {
  let arr = @src.AtomArray::from_atoms([
    bsd_atom("C1", 1, 0.0),
    bsd_atom("C1", 1, 0.0),
    bsd_atom("C2", 1, 1.0),
  ])
  let dedup = @src.biotite_deduplicate_structure(arr)
  assert_eq(dedup.len(), 2)
}

///|
test "biotite_filter_identical" {
  let a = @src.AtomArray::from_atoms([
    bsd_atom("C1", 1, 0.0),
    bsd_atom("C2", 1, 1.0),
  ])
  let b = @src.AtomArray::from_atoms([bsd_atom("C1", 1, 0.0)])
  let mask = @src.biotite_filter_identical(a, b)
  assert_eq(mask.length(), 2)
  assert_true(mask[0])
  assert_true(!mask[1])
}

///|
test "biotite_deduplicate_structures" {
  let a = @src.AtomArray::from_atoms([bsd_atom("C1", 1, 0.0)])
  let b = @src.AtomArray::from_atoms([bsd_atom("C1", 1, 0.0)])
  let c = @src.AtomArray::from_atoms([bsd_atom("C1", 1, 2.0)])
  let result = @src.biotite_deduplicate_structures([a, b, c])
  assert_eq(result.length(), 2)
}

///|
test "biotite_filter_indices" {
  let indices = @src.biotite_filter_indices([true, false, true])
  assert_eq(indices.length(), 2)
  assert_eq(indices[0], 0)
  assert_eq(indices[1], 2)
}
