///|
/// Tests for biotite_structure_rings_ext.mbt.
fn bsr_atom(index : Int32, x : Double, y : Double) -> @src.Atom {
  @src.Atom::new(
    name="C" + index.to_string(),
    coord=@src.Vector3::new(x, y, 0.0),
    resname="BEN",
    element="C",
  )
}

///|
fn benzene() -> Array[@src.Atom] {
  let r = 1.39
  [
    bsr_atom(1, r, 0.0),
    bsr_atom(2, r / 2.0, r * 0.866),
    bsr_atom(3, -r / 2.0, r * 0.866),
    bsr_atom(4, -r, 0.0),
    bsr_atom(5, -r / 2.0, -r * 0.866),
    bsr_atom(6, r / 2.0, -r * 0.866),
  ]
}

///|
test "biotite_ring_orders_benzene" {
  let atoms = benzene()
  let orders = @src.biotite_get_ring_orders(atoms)
  assert_eq(orders.length(), 6)
  for i in 0..<6 {
    assert_eq(orders[i], 6)
  }
}

///|
test "biotite_ring_histogram_benzene" {
  let histogram = @src.biotite_ring_histogram(benzene())
  assert_eq(histogram.length(), 8)
  assert_eq(histogram[6], 1)
  assert_eq(histogram[5], 0)
}

///|
test "biotite_aromatic_ring_count" {
  assert_eq(@src.biotite_count_aromatic_rings(benzene()), 1)
}

///|
test "biotite_ring_atoms_and_acyclic" {
  assert_eq(@src.biotite_ring_atoms(benzene()).length(), 6)
  let chain : Array[@src.Atom] = [
    bsr_atom(1, 0.0, 0.0),
    bsr_atom(2, 1.5, 0.0),
    bsr_atom(3, 3.0, 0.0),
  ]
  assert_eq(@src.biotite_get_ring_orders(chain).length(), 3)
  assert_eq(@src.biotite_ring_atoms(chain).length(), 0)
}
