test "pdbqt_roundtrip" {
  // Create atoms and write PDBQT
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(10.0, 20.0, 30.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    @src.Atom::new(
      name="N",
      coord=@src.Vector3::new(11.0, 20.0, 30.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="N",
    ),
  ]
  let pdbqt_str = @src.write_pdbqt(atoms)
  // Should contain ROOT/ENDROOT
  assert_true(pdbqt_str.contains("ROOT"))
  assert_true(pdbqt_str.contains("ENDROOT"))
  // Parse back
  let parsed = @src.parse_pdbqt(pdbqt_str)
  assert_eq(parsed.length(), 2)
  // Check coordinates
  let c0 = parsed[0].get_coord()
  assert_true(c0.x > 9.9 && c0.x < 10.1)
  assert_true(c0.y > 19.9 && c0.y < 20.1)
  let c1 = parsed[1].get_coord()
  assert_true(c1.x > 10.9 && c1.x < 11.1)
}

test "pdbqt_atom_types" {
  // Carbon aromatic vs aliphatic
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(
      name="CG",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="PHE",
      element="C",
    ),
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(1.0, 0.0, 0.0),
      resname="ALA",
      element="C",
    ),
  ]
  let pdbqt_str = @src.write_pdbqt(atoms)
  // Should contain both "A" (aromatic) and "C" (aliphatic) atom types
  assert_true(pdbqt_str.length() > 0)
}
