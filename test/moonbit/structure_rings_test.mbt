test "detect_rings_benzene" {
  // Benzene ring: 6 carbons in a hexagonal arrangement
  let r = 1.39  // C-C bond length in benzene
  let atoms : Array[@src.Atom] = [
    @src.Atom::new(name="C1", coord=@src.Vector3::new(r, 0.0, 0.0), resname="BEN", element="C"),
    @src.Atom::new(name="C2", coord=@src.Vector3::new(r / 2.0, r * 0.866, 0.0), resname="BEN", element="C"),
    @src.Atom::new(name="C3", coord=@src.Vector3::new(-r / 2.0, r * 0.866, 0.0), resname="BEN", element="C"),
    @src.Atom::new(name="C4", coord=@src.Vector3::new(-r, 0.0, 0.0), resname="BEN", element="C"),
    @src.Atom::new(name="C5", coord=@src.Vector3::new(-r / 2.0, -r * 0.866, 0.0), resname="BEN", element="C"),
    @src.Atom::new(name="C6", coord=@src.Vector3::new(r / 2.0, -r * 0.866, 0.0), resname="BEN", element="C"),
  ]
  let rings = @src.detect_rings(atoms)
  // Should find at least one 6-membered ring
  let mut found_six = false
  for ring in rings {
    if ring.size == 6 {
      found_six = true
    }
  }
  assert_true(found_six)
}
