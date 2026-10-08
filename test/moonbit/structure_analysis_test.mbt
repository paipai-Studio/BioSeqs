test "detect_molecules_water" {
  // Two water molecules far apart
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(name="O", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="HOH", element="O"),
    @src.Atom::new(name="H1", coord=@src.Vector3::new(0.9, 0.0, 0.0), resname="HOH", element="H"),
    @src.Atom::new(name="H2", coord=@src.Vector3::new(0.0, 0.9, 0.0), resname="HOH", element="H"),
    @src.Atom::new(name="O", coord=@src.Vector3::new(20.0, 20.0, 20.0), resname="HOH", element="O"),
    @src.Atom::new(name="H1", coord=@src.Vector3::new(20.9, 20.0, 20.0), resname="HOH", element="H"),
    @src.Atom::new(name="H2", coord=@src.Vector3::new(20.0, 20.9, 20.0), resname="HOH", element="H"),
  ]
  let mols = @src.detect_molecules(atoms)
  assert_eq(mols.length(), 2)
  // Largest molecule first
  assert_true(mols[0].n_atoms >= 3)
}

test "detect_molecules_single_bonded" {
  // Two bonded carbons
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="LIG", element="C"),
    @src.Atom::new(name="C2", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="LIG", element="C"),
  ]
  let mols = @src.detect_molecules(atoms)
  assert_eq(mols.length(), 1)
  assert_eq(mols[0].n_atoms, 2)
}

test "rdf_basic" {
  // Simple cubic lattice of atoms
  let atoms : Array[@src.Atom] = Array(capacity=27)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      for k = 0; k < 3; k = k + 1 {
        atoms.push(@src.Atom::new(
          name="C",
          coord=@src.Vector3::new(i.to_double() * 3.0, j.to_double() * 3.0, k.to_double() * 3.0),
          resname="C",
          element="C",
        ))
      }
    }
  }
  let result = @src.compute_rdf(atoms, r_max=10.0, n_bins=50)
  assert_eq(result.n_bins, 50)
  assert_eq(result.r.length(), 50)
  assert_eq(result.g_r.length(), 50)
  // First bin should be 0 (no atoms at r≈0.1)
  assert_true(result.g_r[0] == 0.0)
  // Some bin should have non-zero g(r) (at r≈3.0 Å)
  let mut found_nonzero = false
  for b = 0; b < 50; b = b + 1 {
    if result.g_r[b] > 0.0 {
      found_nonzero = true
    }
  }
  assert_true(found_nonzero)
}
