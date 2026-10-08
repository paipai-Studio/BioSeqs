test "detect_basepairs_wc" {
  // Create a simple A-U Watson-Crick base pair
  // Place N1 of U and N9 of A at reasonable distance
  let atoms : Array[@src.Atom] = [
    // Adenine representative atom (N9)
    @src.Atom::new(name="N9", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="A", chainid='A', resseq=1, element="N"),
    @src.Atom::new(name="N1", coord=@src.Vector3::new(1.0, 0.0, 0.0), resname="A", chainid='A', resseq=1, element="N"),
    @src.Atom::new(name="N6", coord=@src.Vector3::new(2.0, 0.0, 0.0), resname="A", chainid='A', resseq=1, element="N"),
    // Uracil representative atom (N1)
    @src.Atom::new(name="N1", coord=@src.Vector3::new(4.0, 0.0, 0.0), resname="U", chainid='B', resseq=2, element="N"),
    @src.Atom::new(name="O4", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="U", chainid='B', resseq=2, element="O"),
    @src.Atom::new(name="O2", coord=@src.Vector3::new(3.0, 0.0, 0.0), resname="U", chainid='B', resseq=2, element="O"),
  ]
  let bps = @src.detect_basepairs(atoms)
  // Should detect at least one base pair
  assert_true(bps.length() >= 1)
  // Check pair type
  let mut has_wc = false
  for bp in bps {
    if bp.pair_type == "WC" {
      has_wc = true
    }
  }
  assert_true(has_wc)
}
