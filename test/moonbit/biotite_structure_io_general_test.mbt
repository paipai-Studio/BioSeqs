///|
/// Tests for biotite_structure_io_general.mbt.
fn bigs_atom(name : String, element : String, x : Double) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element~,
  )
}

///|
test "biotite_guess_format" {
  assert_eq(@src.biotite_guess_format("a.pdb"), "pdb")
  assert_eq(@src.biotite_guess_format("a.xyz"), "xyz")
  assert_eq(@src.biotite_guess_format("a.mmcif"), "mmcif")
  assert_eq(@src.biotite_guess_format("a.gro"), "gro")
  assert_eq(@src.biotite_guess_format("a.txt"), "")
}

///|
test "biotite_load_save_pdb" {
  let arr = @src.AtomArray::from_atoms([
    bigs_atom("CA", "C", 0.0),
    bigs_atom("CB", "C", 1.0),
  ])
  let pdb = @src.biotite_save_structure(arr, "pdb")
  assert_true(pdb.has_prefix("ATOM"))
  let loaded = @src.biotite_load_structure(pdb, "pdb")
  assert_eq(loaded.len(), 2)
}

///|
test "biotite_load_save_xyz" {
  let arr = @src.AtomArray::from_atoms([bigs_atom("C1", "C", 0.0)])
  let xyz = @src.biotite_save_structure(arr, "xyz")
  let loaded = @src.biotite_load_structure(xyz, "xyz")
  assert_eq(loaded.len(), 1)
}

///|
test "biotite_load_save_mmcif" {
  let arr = @src.AtomArray::from_atoms([
    bigs_atom("CA", "C", 0.0),
    bigs_atom("CB", "C", 1.0),
  ])
  let cif = @src.biotite_save_structure(arr, "mmcif")
  let loaded = @src.biotite_load_structure(cif, "mmcif")
  assert_eq(loaded.len(), 2)
}

///|
test "biotite_load_unknown_format" {
  assert_eq(@src.biotite_load_structure("", "unknown").len(), 0)
  assert_eq(@src.biotite_save_structure(@src.AtomArray::new(), "unknown"), "")
}
