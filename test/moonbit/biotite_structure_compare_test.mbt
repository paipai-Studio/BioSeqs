///|
/// Tests for biotite_structure_compare.mbt.
fn bsc_atom(name : String, x : Double) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
  )
}

///|
test "biotite_compare_identical" {
  let a = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0), bsc_atom("CB", 1.0)])
  let b = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0), bsc_atom("CB", 1.0)])
  let result = @src.biotite_compare(a, b)
  assert_eq(result.n_identical, 2)
  assert_eq(result.n_different, 0)
  assert_true(@moonbit.assert_num_close(result.rmsd, 0.0, 1.0e-12))
}

///|
test "biotite_compare_shifted" {
  let a = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0), bsc_atom("CB", 1.0)])
  let b = @src.AtomArray::from_atoms([bsc_atom("CA", 0.5), bsc_atom("CB", 1.5)])
  let result = @src.biotite_compare(a, b)
  assert_eq(result.n_identical, 0)
  assert_eq(result.n_different, 2)
  assert_true(@moonbit.assert_num_close(result.rmsd, 0.5, 1.0e-12))
}

///|
test "biotite_compare_length_mismatch" {
  let a = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0)])
  let b = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0), bsc_atom("CB", 1.0)])
  let result = @src.biotite_compare(a, b)
  assert_eq(result.n_identical, 1)
  assert_eq(result.n_different, 1)
  assert_eq(result.identical.length(), 2)
}

///|
test "biotite_get_identical_atoms" {
  let a = @src.AtomArray::from_atoms([bsc_atom("CA", 0.0), bsc_atom("CB", 1.0)])
  let b = @src.AtomArray::from_atoms([bsc_atom("CB", 1.0), bsc_atom("CA", 0.0)])
  let pairs = @src.biotite_get_identical_atoms(a, b)
  assert_eq(pairs.length(), 2)
}
