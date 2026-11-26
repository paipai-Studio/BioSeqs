///|
/// Tests for biotite_structure_celllist.mbt (official `CellList` API:
/// get_atoms / get_atoms_in_cells / create_adjacency_matrix with the
/// semantics of the official Rust implementation, verified against
/// the official test_celllist.py expectations).
fn bcl_vec(x : Double, y : Double, z : Double) -> @src.Vector3 {
  @src.Vector3::new(x, y, z)
}

///|
/// Five atoms along the z axis (official test_get_atoms molecule).
fn bcl_z_axis() -> Array[@src.Vector3] {
  [
    bcl_vec(0.0, 0.0, 0.0),
    bcl_vec(0.0, 0.0, 1.0),
    bcl_vec(0.0, 0.0, 2.0),
    bcl_vec(0.0, 0.0, 3.0),
    bcl_vec(0.0, 0.0, 4.0),
  ]
}

///|
test "celllist get atoms official expectations" {
  // Official test_get_atoms: results must be independent of the cell
  // size (parametrized 0.5, 1, 2, 5, 10 upstream).
  let sizes = [0.5, 1.0, 2.0, 5.0, 10.0]
  for cs in sizes {
    let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=cs)
    let r1 = cl.get_atoms(bcl_vec(0.0, 0.0, 0.1), 1.0)
    assert_eq(r1.length(), 2)
    assert_eq(r1[0], 0)
    assert_eq(r1[1], 1)
    let r2 = cl.get_atoms(bcl_vec(0.0, 0.0, 1.1), 1.0)
    assert_eq(r2.length(), 2)
    assert_eq(r2[0], 1)
    assert_eq(r2[1], 2)
    let r3 = cl.get_atoms(bcl_vec(0.0, 0.0, 1.1), 2.0)
    assert_eq(r3.length(), 4)
    assert_eq(r3[0], 0)
    assert_eq(r3[1], 1)
    assert_eq(r3[2], 2)
    assert_eq(r3[3], 3)
  }
}

///|
test "celllist get atoms pairs official order" {
  // Official MAPPING order over multiple positions (radius 2):
  // queries 0/1/2 find atoms [0,1,2], [0,1,2,3] and [3,4] in cell
  // traversal order.
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=1.0)
  let pos = [
    bcl_vec(0.0, 0.0, 0.1),
    bcl_vec(0.0, 0.0, 1.1),
    bcl_vec(0.0, 0.0, 4.1),
  ]
  let pairs = cl.get_atoms_pairs(pos, radius=Some(2.0))
  // (query, atom) tuples in official order.
  let expected : Array[(Int32, Int32)] = [
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 3),
    (2, 4),
  ]
  assert_eq(pairs.length(), expected.length())
  for k = 0; k < expected.length(); k = k + 1 {
    assert_eq(pairs[k].0, expected[k].0)
    assert_eq(pairs[k].1, expected[k].1)
  }
}

///|
test "celllist per coordinate radii" {
  // Official single-vs-multiple radii: rad = [1.0, 2.0, 3.0]
  // yields flat atom indices [0,1, 0,1,2,3, 2,3,4].
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=1.0)
  let pos = [
    bcl_vec(0.0, 0.0, 0.1),
    bcl_vec(0.0, 0.0, 1.1),
    bcl_vec(0.0, 0.0, 4.1),
  ]
  let pairs = cl.get_atoms_pairs(pos, radii=Some([1.0, 2.0, 3.0]))
  let expected : Array[(Int32, Int32)] = [
    (0, 0),
    (0, 1),
    (1, 0),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (2, 4),
  ]
  assert_eq(pairs.length(), expected.length())
  for k = 0; k < expected.length(); k = k + 1 {
    assert_eq(pairs[k].0, expected[k].0)
    assert_eq(pairs[k].1, expected[k].1)
  }
}

///|
test "celllist mask format" {
  // Official MASK result: mask[q][j] = dist(q, atom j) <= radius.
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=2.0)
  let mask = cl.get_atoms_mask([bcl_vec(0.0, 0.0, 0.1)], 1.0)
  assert_eq(mask.length(), 1)
  assert_eq(mask[0].length(), 5)
  assert_true(mask[0][0])
  assert_true(mask[0][1])
  assert_false(mask[0][2])
  assert_false(mask[0][3])
  assert_false(mask[0][4])
}

///|
test "celllist adjacency matrix matches brute force" {
  // Official test_adjacency_matrix semantics:
  // m[i,j] = dist(i,j) <= threshold (boundary inclusive).
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=2.0)
  let m = cl.create_adjacency_matrix(1.0)
  assert_eq(m.length(), 5)
  // Self pairs: distance 0 <= 1.
  for i = 0; i < 5; i = i + 1 {
    assert_true(m[i][i])
  }
  // Adjacent atoms at distance exactly 1.0 are included (<=).
  assert_true(m[0][1])
  assert_true(m[1][0])
  // Distance 2.0 exceeds the threshold.
  assert_false(m[0][2])
  // Distance 3.0 and 4.0 excluded.
  assert_false(m[0][3])
  assert_false(m[0][4])
  // Threshold 2.0 includes the distance-2.0 pairs.
  let m2 = cl.create_adjacency_matrix(2.0)
  assert_true(m2[0][2])
  assert_true(m2[2][0])
  assert_false(m2[0][3])
  // Threshold 10 includes everything.
  let m10 = cl.create_adjacency_matrix(10.0)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 5; j = j + 1 {
      assert_true(m10[i][j])
    }
  }
  // Negative threshold: all-False matrix (official raises ValueError).
  let mneg = cl.create_adjacency_matrix(-1.0)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 5; j = j + 1 {
      assert_false(mneg[i][j])
    }
  }
}

///|
test "celllist outside location yields empty" {
  // Official test_outside_location: a query far outside the grid
  // finds no atoms.
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=5.0)
  assert_eq(cl.get_atoms(bcl_vec(1000.0, 1000.0, 1000.0), 5.0).length(), 0)
  assert_eq(cl.get_atoms(bcl_vec(-1000.0, -1000.0, -1000.0), 5.0).length(), 0)
}

///|
test "celllist get atoms in cells" {
  // Official get_atoms_in_cells: no distance filter, cell_radius
  // counts cells (default 1: the cell plus its 26 neighbors).
  let cl = @src.biotite_celllist_new(bcl_z_axis(), cell_size=1.0)
  // Query at z=1.1 sits in cell z=1; radius 0 -> only that cell,
  // which holds just atom 1 (z=1).
  let r0 = cl.get_atoms_in_cells(bcl_vec(0.0, 0.0, 1.1), cell_radius=0)
  assert_eq(r0.length(), 1)
  assert_eq(r0[0], 1)
  // Default radius 1 -> cells z=0..2 -> atoms 0..2 in traversal order.
  let r1 = cl.get_atoms_in_cells(bcl_vec(0.0, 0.0, 1.1))
  assert_eq(r1.length(), 3)
  assert_eq(r1[0], 0)
  assert_eq(r1[2], 2)
}

///|
test "celllist selection stores masked atoms only" {
  // Official test_selection: masked atoms are not stored, but the
  // returned indices still refer to the original array.
  let sel = [true, false, true, false, true]
  let cl = @src.biotite_celllist_new(
    bcl_z_axis(),
    cell_size=2.0,
    selection=Some(sel),
  )
  let r = cl.get_atoms(bcl_vec(0.0, 0.0, 0.1), 2.0)
  // Atoms 0 (d=0.1) and 2 (d=1.9) qualify by distance; atom 1 is
  // masked out.
  assert_eq(r.length(), 2)
  assert_eq(r[0], 0)
  assert_eq(r[1], 2)
  // Adjacency matrix: rows/columns of unselected atoms are all False.
  let m = cl.create_adjacency_matrix(10.0)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 5; j = j + 1 {
      let expected = sel[i] && sel[j]
      assert_eq(m[i][j], expected)
    }
  }
}

///|
test "celllist empty input yields empty results" {
  // Official test_empty_coordinates: empty queries give empty
  // outputs. Construction without any coordinate raises upstream
  // ("Coordinates must not be empty"); here it yields an empty grid.
  let cl = @src.biotite_celllist_new([], cell_size=2.0)
  assert_eq(cl.get_atoms(bcl_vec(0.0, 0.0, 0.0), 1.0).length(), 0)
  assert_eq(cl.get_atoms_in_cells(bcl_vec(0.0, 0.0, 0.0)).length(), 0)
  assert_eq(cl.get_atoms_pairs([], radius=Some(1.0)).length(), 0)
  let mask = cl.get_atoms_mask([], 1.0)
  assert_eq(mask.length(), 0)
}

///|
test "celllist boundary inclusive radius" {
  // distance == radius is included (official adjacency semantics:
  // m[i,j] = dist <= threshold; get_atoms uses the same filter).
  let cl = @src.biotite_celllist_new(
    [bcl_vec(0.0, 0.0, 0.0), bcl_vec(2.0, 0.0, 0.0)],
    cell_size=2.0,
  )
  let r = cl.get_atoms(bcl_vec(0.0, 0.0, 0.0), 2.0)
  assert_eq(r.length(), 2)
  assert_eq(r[0], 0)
  assert_eq(r[1], 1)
  let r_small = cl.get_atoms(bcl_vec(0.0, 0.0, 0.0), 1.999)
  assert_eq(r_small.length(), 1)
}
