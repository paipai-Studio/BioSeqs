///|
test "neighbor_search_new" {
  let atoms = [
    @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    @src.Atom::new(name="CB", coord=@src.Vector3::new(1.0, 0.0, 0.0)),
    @src.Atom::new(name="CC", coord=@src.Vector3::new(2.0, 0.0, 0.0)),
  ]

  let ns = @src.NeighborSearch::new(atoms)
  assert_eq(ns.size(), 3)
}

///|
test "neighbor_search_search_radius" {
  let atoms = [
    @src.Atom::new(name="A1", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    @src.Atom::new(name="A2", coord=@src.Vector3::new(1.0, 0.0, 0.0)),
    @src.Atom::new(name="A3", coord=@src.Vector3::new(3.0, 0.0, 0.0)),
    @src.Atom::new(name="A4", coord=@src.Vector3::new(0.5, 0.5, 0.5)),
  ]

  let ns = @src.NeighborSearch::new(atoms)
  let center = @src.Vector3::new(0.0, 0.0, 0.0)
  let result = ns.search(center, 1.5)

  assert_true(result.length() >= 2)
}

///|
test "neighbor_search_nearest" {
  let atoms = [
    @src.Atom::new(name="A1", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    @src.Atom::new(name="A2", coord=@src.Vector3::new(2.0, 0.0, 0.0)),
    @src.Atom::new(name="A3", coord=@src.Vector3::new(5.0, 0.0, 0.0)),
  ]

  let ns = @src.NeighborSearch::new(atoms)
  let center = @src.Vector3::new(0.5, 0.0, 0.0)
  let nearest = ns.nearest(center)

  assert_true(nearest.is_some())
}

///|
test "neighbor_search_empty" {
  let atoms : Array[@src.Atom] = []
  let ns = @src.NeighborSearch::new(atoms)

  assert_eq(ns.size(), 0)

  let center = @src.Vector3::new(0.0, 0.0, 0.0)
  let result = ns.search(center, 10.0)
  assert_eq(result.length(), 0)

  let nearest = ns.nearest(center)
  assert_true(nearest.is_none())
}

///|
test "neighbor_search_large_radius" {
  let atoms = [
    @src.Atom::new(name="A1", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    @src.Atom::new(name="A2", coord=@src.Vector3::new(10.0, 10.0, 10.0)),
    @src.Atom::new(name="A3", coord=@src.Vector3::new(-5.0, 3.0, 2.0)),
  ]

  let ns = @src.NeighborSearch::new(atoms)
  let center = @src.Vector3::new(0.0, 0.0, 0.0)
  let result = ns.search(center, 100.0)

  assert_eq(result.length(), 3)
}

///|
test "neighbor_search_search_pairs" {
  let atoms = [
    @src.Atom::new(name="A1", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    @src.Atom::new(name="A2", coord=@src.Vector3::new(1.0, 0.0, 0.0)),
    @src.Atom::new(name="A3", coord=@src.Vector3::new(10.0, 0.0, 0.0)),
    @src.Atom::new(name="A4", coord=@src.Vector3::new(0.5, 0.5, 0.0)),
  ]

  let ns = @src.NeighborSearch::new(atoms)
  let pairs = ns.search_all(2.0)

  assert_true(pairs.length() >= 1)
}
