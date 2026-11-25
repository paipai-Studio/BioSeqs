///|
/// Tests for biotite_interface_mdtraj.mbt.
fn bim_template() -> @src.AtomArray {
  @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    @src.Atom::new(
      name="CB",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
  ])
}

///|
test "biotite_mdtraj_topology" {
  let text = @src.biotite_mdtraj_topology(bim_template())
  assert_true(text.contains("ALA 1 CA C"))
  assert_true(text.contains("ALA 1 CB C"))
}

///|
test "biotite_mdtraj_coords_roundtrip" {
  let template = bim_template()
  let coords = @src.biotite_mdtraj_coords(template)
  assert_eq(coords.length(), 2)
  assert_eq(coords[0].length(), 3)
}

///|
test "biotite_mdtraj_from_coords" {
  let template = bim_template()
  let rebuilt = @src.biotite_mdtraj_from_coords(template, [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ])
  assert_eq(rebuilt.len(), 2)
  assert_true(@moonbit.assert_num_close(rebuilt.get_coord(0).x, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rebuilt.get_coord(1).z, 6.0, 1.0e-12))
  assert_eq(rebuilt.get_name(0), "CA")
}

///|
test "biotite_mdtraj_from_coords_short_input" {
  let rebuilt = @src.biotite_mdtraj_from_coords(bim_template(), [[1.0, 2.0]])
  assert_eq(rebuilt.len(), 1)
  assert_true(@moonbit.assert_num_close(rebuilt.get_coord(0).z, 0.0, 1.0e-12))
}
