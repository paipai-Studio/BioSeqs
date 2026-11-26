///|
/// Tests for biotite_structure_density.mbt (official histogramdd
/// semantics: even delta grid, per-bin counts/weights, optional
/// probability-density normalization).
fn bsd_elem_atom(
  element : String,
  x : Double,
  y : Double,
  z : Double,
) -> @src.Atom {
  @src.Atom::new(
    name="A",
    coord=@src.Vector3::new(x, y, z),
    resname="LIG",
    chainid='A',
    resseq=1,
    element~,
  )
}

///|
test "density single atom histogram" {
  let arr = @src.AtomArray::from_atoms([bsd_elem_atom("C", 0.0, 0.0, 0.0)])
  let d = @src.biotite_density(arr.coords, delta=1.0)
  let (nx, ny, nz) = d.shape()
  // one atom at the origin: single bin in each axis
  assert_eq(nx, 1)
  assert_eq(ny, 1)
  assert_eq(nz, 1)
  assert_true(@moonbit.assert_num_close(d.total(), 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.max(), 1.0, 1.0e-12))
  // edges start at the coordinate minimum
  assert_true(@moonbit.assert_num_close(d.get_edges_x()[0], 0.0, 1.0e-12))
}

///|
test "density separated atoms count separately" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 3.0, 0.0, 0.0),
  ])
  let d = @src.biotite_density(arr.coords, delta=1.0)
  let (nx, ny, nz) = d.shape()
  // span 3.0 with delta 1.0 -> 3 bins along x
  assert_eq(nx, 3)
  assert_eq(ny, 1)
  assert_eq(nz, 1)
  assert_true(@moonbit.assert_num_close(d.total(), 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.max(), 1.0, 1.0e-12))
}

///|
test "density weights accumulate" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 0.1, 0.0, 0.0),
  ])
  // both atoms share the first bin
  let d = @src.biotite_density(arr.coords, delta=1.0)
  assert_true(@moonbit.assert_num_close(d.max(), 2.0, 1.0e-12))
  // weighted: one heavy, one light atom in the same bin
  let dw = @src.biotite_density(arr.coords, delta=1.0, weights=Some([5.0, 1.0]))
  assert_true(@moonbit.assert_num_close(dw.max(), 6.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(dw.total(), 6.0, 1.0e-12))
}

///|
test "density probability density normalization" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 3.0, 0.0, 0.0),
  ])
  let d = @src.biotite_density(arr.coords, delta=1.0, density=true)
  // two atoms, two occupied bins of volume 1 -> density 0.5 per bin
  assert_true(@moonbit.assert_num_close(d.total(), 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.max(), 0.5, 1.0e-12))
}

///|
test "density z axis binning" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 0.0, 0.0, 2.5),
  ])
  let d = @src.biotite_density(arr.coords, delta=1.0)
  let (nx, ny, nz) = d.shape()
  assert_eq(nx, 1)
  assert_eq(ny, 1)
  // span 2.5 -> 3 bins along z (0-1, 1-2, 2-3)
  assert_eq(nz, 3)
  assert_true(@moonbit.assert_num_close(d.total(), 2.0, 1.0e-12))
  // z edges: 0, 1, 2, 3
  assert_eq(d.get_edges_z().length(), 4)
  assert_true(@moonbit.assert_num_close(d.get_edges_z()[3], 3.0, 1.0e-12))
}

///|
test "density edges cover coordinate range" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", -2.0, 0.0, 0.0),
    bsd_elem_atom("C", 2.0, 0.0, 0.0),
  ])
  let d = @src.biotite_density(arr.coords, delta=1.0)
  // span 4.0 -> 4 bins along x
  let (nx, _ny, _nz) = d.shape()
  assert_eq(nx, 4)
  // edges from -2.0 stepping by 1.0
  assert_true(@moonbit.assert_num_close(d.get_edges_x()[0], -2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.get_edges_x()[4], 2.0, 1.0e-12))
}

///|
test "density empty coordinates" {
  let d = @src.biotite_density([], delta=1.0)
  let (nx, ny, nz) = d.shape()
  assert_eq(nx, 1)
  assert_eq(ny, 1)
  assert_eq(nz, 1)
  assert_true(@moonbit.assert_num_close(d.total(), 0.0, 1.0e-12))
}

///|
test "density non-integer span rounds up" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 2.5, 0.0, 0.0),
  ])
  let d = @src.biotite_density(arr.coords, delta=1.0)
  let (nx, _ny, _nz) = d.shape()
  // 2.5 span needs ceil(2.5) = 2..3 bins; atom at 2.5 clamps into the
  // last bin
  assert_true(nx == 3)
  assert_true(@moonbit.assert_num_close(d.total(), 2.0, 1.0e-12))
}
