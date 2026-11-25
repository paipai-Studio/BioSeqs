///|
/// Tests for biotite_structure_annotation.mbt.
fn bsa_atom(
  name : String,
  res_name : String,
  res_id : Int32,
  x : Double,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname=res_name,
    chainid='A',
    resseq=res_id,
    element="C",
  )
}

///|
test "biotite_annotate_amino_acids" {
  let arr = @src.AtomArray::from_atoms([
    bsa_atom("CA", "ALA", 1, 0.0),
    bsa_atom("CA", "GLY", 2, 1.0),
    bsa_atom("O", "HOH", 3, 2.0),
  ])
  let flags = @src.biotite_annotate_amino_acids(arr)
  assert_eq(flags.length(), 3)
  assert_true(flags[0])
  assert_true(flags[1])
  assert_true(!flags[2])
}

///|
test "biotite_annotate_residues" {
  let arr = @src.AtomArray::from_atoms([
    bsa_atom("CA", "ALA", 1, 0.0),
    bsa_atom("CA", "GLY", 2, 1.0),
  ])
  let residues = @src.biotite_annotate_residues(arr)
  assert_eq(residues.length(), 2)
  assert_eq(residues[0].0, "A")
  assert_eq(residues[0].1, 1)
  assert_eq(residues[0].2, "ALA")
}

///|
test "biotite_annotate_sse" {
  let arr = @src.AtomArray::from_atoms([
    bsa_atom("CA", "ALA", 1, 0.0),
    bsa_atom("CA", "GLY", 2, 1.0),
  ])
  assert_eq(@src.biotite_annotate_sse(arr).length(), 2)
  assert_eq(@src.biotite_annotate_sse_by_residue(arr).length(), 2)
}
