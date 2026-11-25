///|
/// Tests for biotite_structure_filter.mbt — filters and AtomArray adapters.
fn bf_info(
  atom_name : String,
  res_name : String,
  chain_id : String,
  res_id : Int32,
  hetero : Bool,
) -> @src.BiotiteAtomInfo {
  { atom_name, res_name, chain_id, res_id, hetero, }
}

///|
test "biotite_filter_residue_classes" {
  let atoms : Array[@src.BiotiteAtomInfo] = [
    bf_info("CA", "ALA", "A", 1, false),
    bf_info("O", "HOH", "A", 2, true),
    bf_info("C1'", "DA", "B", 1, false),
  ]
  assert_eq(@src.biotite_filter_amino_acids(atoms).length(), 1)
  assert_eq(@src.biotite_filter_backbone(atoms).length(), 1)
  assert_eq(@src.biotite_filter_nucleotides(atoms).length(), 1)
  assert_eq(@src.biotite_filter_hetero(atoms).length(), 1)
}

///|
test "biotite_filter_chain_and_residue" {
  let atoms : Array[@src.BiotiteAtomInfo] = [
    bf_info("CA", "ALA", "A", 1, false),
    bf_info("CA", "GLY", "B", 1, false),
  ]
  assert_eq(@src.biotite_filter_chain(atoms, "B").length(), 1)
  assert_eq(@src.biotite_filter_residue(atoms, "ALA").length(), 1)
}

///|
test "biotite_segment_residue_chain_starts" {
  let atoms : Array[@src.BiotiteAtomInfo] = [
    bf_info("CA", "ALA", "A", 1, false),
    bf_info("CB", "ALA", "A", 1, false),
    bf_info("CA", "GLY", "A", 2, false),
    bf_info("CA", "GLY", "B", 1, false),
  ]
  assert_eq(@src.biotite_get_chain_starts(atoms).length(), 2)
  assert_eq(@src.biotite_get_residue_starts(atoms).length(), 3)
  assert_eq(@src.biotite_get_segment_starts(atoms).length(), 3)
}

///|
test "biotite_filter_atom_array_adapters" {
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
    ),
    @src.Atom::new(
      name="O",
      coord=@src.Vector3::new(1.0, 0.0, 0.0),
      resname="HOH",
      chainid='A',
      resseq=2,
    ),
    @src.Atom::new(
      name="C1'",
      coord=@src.Vector3::new(2.0, 0.0, 0.0),
      resname="DA",
      chainid='B',
      resseq=1,
    ),
  ]
  let array = @src.AtomArray::from_atoms(atoms)
  assert_eq(@src.biotite_filter_amino_acids_array(array).len(), 1)
  assert_eq(@src.biotite_filter_nucleotides_array(array).len(), 1)
  assert_eq(@src.biotite_filter_residue_array(array, "HOH").len(), 1)
  assert_eq(@src.biotite_get_chain_starts_array(array).length(), 2)
  assert_eq(@src.biotite_get_residue_starts_array(array).length(), 3)
  let infos = @src.biotite_atom_info_from_array(array)
  assert_eq(infos.length(), 3)
  assert_eq(infos[1].res_name, "HOH")
}
