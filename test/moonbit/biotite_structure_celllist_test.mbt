///|
/// Tests for biotite_structure_celllist.mbt.
test "celllist adjacent atoms" {
  let cl = @src.biotite_celllist_new(
    [
      @src.Vector3::new(0.0, 0.0, 0.0),
      @src.Vector3::new(1.0, 0.0, 0.0),
      @src.Vector3::new(10.0, 0.0, 0.0),
    ],
    cell_size=2.0,
  )
  let adj = cl.adjacent_atoms(0, 1.5)
  assert_eq(adj.length(), 1)
  assert_eq(adj[0], 1)
  // atom 2 is far away
  let adj2 = cl.adjacent_atoms(0, 5.0)
  assert_eq(adj2.length(), 1)
  // far atom not included even at bigger radius
  let adj3 = cl.adjacent_atoms(0, 8.0)
  assert_eq(adj3.length(), 1)
}

///|
test "celllist nearest neighbors" {
  let cl = @src.biotite_celllist_new(
    [
      @src.Vector3::new(0.0, 0.0, 0.0),
      @src.Vector3::new(1.0, 0.0, 0.0),
      @src.Vector3::new(3.0, 0.0, 0.0),
      @src.Vector3::new(7.0, 0.0, 0.0),
    ],
    cell_size=2.0,
  )
  let (idx, dists) = cl.nearest_neighbors(0, number=2)
  assert_eq(idx.length(), 2)
  assert_eq(idx[0], 1)
  assert_eq(idx[1], 2)
  assert_true(@moonbit.assert_num_close(dists[0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(dists[1], 3.0, 1.0e-12))
  // with distance cap
  let (idx2, dists2) = cl.nearest_neighbors(0, number=3, distance=Some(3.5))
  assert_eq(idx2.length(), 2)
  assert_true(@moonbit.assert_num_close(dists2[1], 3.0, 1.0e-12))
  // out-of-range index
  let (idx3, _) = cl.nearest_neighbors(-1, number=1)
  assert_eq(idx3.length(), 0)
}

///|
test "celllist conflict" {
  let cl = @src.biotite_celllist_new(
    [
      @src.Vector3::new(0.0, 0.0, 0.0),
      @src.Vector3::new(0.5, 0.0, 0.0),
      @src.Vector3::new(5.0, 0.0, 0.0),
    ],
    cell_size=2.0,
  )
  let conflicts = cl.conflict(1.0)
  assert_eq(conflicts.length(), 1)
  let (i, j) = conflicts[0]
  assert_eq(i, 0)
  assert_eq(j, 1)
  // larger threshold still only that pair (5.0 apart)
  assert_eq(cl.conflict(2.0).length(), 1)
}

///|
test "celllist create_by_boxes" {
  let cl = @src.biotite_celllist_new(
    [
      @src.Vector3::new(0.0, 0.0, 0.0),
      @src.Vector3::new(1.0, 0.0, 0.0),
      @src.Vector3::new(5.0, 5.0, 5.0),
    ],
    cell_size=2.0,
  )
  let boxes = cl.create_by_boxes([
    (@src.Vector3::new(-1.0, -1.0, -1.0), @src.Vector3::new(1.5, 1.5, 1.5)),
    (@src.Vector3::new(4.0, 4.0, 4.0), @src.Vector3::new(6.0, 6.0, 6.0)),
  ])
  assert_eq(boxes.length(), 2)
  assert_eq(boxes[0].length(), 2)
  assert_eq(boxes[1].length(), 1)
  assert_eq(boxes[1][0], 2)
}

///|
test "celllist empty positions" {
  let cl = @src.biotite_celllist_new([], cell_size=2.0)
  assert_eq(cl.adjacent_atoms(0, 1.0).length(), 0)
  assert_eq(cl.conflict(1.0).length(), 0)
  assert_eq(cl.create_by_boxes([]).length(), 0)
}
