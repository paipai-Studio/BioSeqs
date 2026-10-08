test "celllist_neighbors" {
  // Place atoms in a grid
  let coords : Array[@src.Vector3] = Array(capacity=8)
  for i = 0; i < 2; i = i + 1 {
    for j = 0; j < 2; j = j + 1 {
      for k = 0; k < 2; k = k + 1 {
        coords.push(@src.Vector3::new(i.to_double() * 3.0, j.to_double() * 3.0, k.to_double() * 3.0))
      }
    }
  }
  let cl = @src.CellList::new(coords, cell_size=5.0)
  // Find neighbors of atom 0 within 3.5 Å
  let neighbors = cl.get_neighbors(coords, 0, radius=3.5)
  // Atom 0 at (0,0,0): neighbors at (3,0,0), (0,3,0), (0,0,3) — all distance 3.0
  assert_eq(neighbors.length(), 3)
}

test "celllist_pairs" {
  let coords : Array[@src.Vector3] = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.5, 0.0, 0.0),
    @src.Vector3::new(10.0, 10.0, 10.0),
  ]
  let cl = @src.CellList::new(coords, cell_size=3.0)
  let pairs = cl.get_all_pairs(coords, radius=2.0)
  // Only (0,1) should be within 2.0 Å
  assert_eq(pairs.length(), 1)
  assert_eq(pairs[0].0, 0)
  assert_eq(pairs[0].1, 1)
}

test "celllist_coordination" {
  // Central atom with 3 neighbors at distance 1.0
  let coords : Array[@src.Vector3] = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 1.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 1.0),
  ]
  let cl = @src.CellList::new(coords, cell_size=2.0)
  let cn = cl.coordination_numbers(coords, radius=1.1)
  // Atom 0 has 3 neighbors at distance 1.0
  assert_eq(cn[0], 3)
  // Atoms 1-3 each have 1 neighbor (atom 0) at distance 1.0
  // (distance between atom 1 and atom 2 is sqrt(2) ≈ 1.41 > 1.1)
  assert_eq(cn[1], 1)
  assert_eq(cn[2], 1)
  assert_eq(cn[3], 1)
}
