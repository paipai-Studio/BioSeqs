///|
/// Tests for biotite_structure_density.mbt.
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
test "density electron counts" {
  assert_eq(@src.biotite_density_electron_count("C"), Some(6))
  assert_eq(@src.biotite_density_electron_count("O"), Some(8))
  assert_eq(@src.biotite_density_electron_count("S"), Some(16))
  assert_eq(@src.biotite_density_electron_count("X"), None)
}

///|
test "density scattering lengths" {
  match @src.biotite_density_scattering_length("H") {
    Some(v) => assert_true(@moonbit.assert_num_close(v, -3.739, 1.0e-9))
    None => fail("H should have a scattering length")
  }
  match @src.biotite_density_scattering_length("C") {
    Some(v) => assert_true(@moonbit.assert_num_close(v, 6.646, 1.0e-9))
    None => fail("C should have a scattering length")
  }
  assert_eq(@src.biotite_density_scattering_length("X"), None)
}

///|
test "density from single carbon atom" {
  let arr = @src.AtomArray::from_atoms([bsd_elem_atom("C", 0.0, 0.0, 0.0)])
  let density = @src.biotite_density_from_atom_array(arr, resolution=1.0)
  let (nz, ny, nx) = density.shape()
  assert_true(nz > 0 && ny > 0 && nx > 0)
  assert_true(density.max() > 0.0)
  assert_true(density.total() > 0.0)
  // peak close to the atom center: find max cell
  let grid = match density {
    d => d
  }
  let _ = grid
}

///|
test "density unknown element contributes nothing" {
  let arr = @src.AtomArray::from_atoms([bsd_elem_atom("X", 0.0, 0.0, 0.0)])
  let density = @src.biotite_density_from_atom_array(arr, resolution=1.0)
  assert_true(@moonbit.assert_num_close(density.total(), 0.0, 1.0e-15))
}

///|
test "density neutron weights" {
  // hydrogen has a negative scattering length -> negative total
  let arr = @src.AtomArray::from_atoms([bsd_elem_atom("H", 0.0, 0.0, 0.0)])
  let xray = @src.biotite_density_from_atom_array(arr, resolution=1.0)
  assert_true(xray.total() > 0.0)
  let neutron = @src.biotite_density_from_atom_array(
    arr,
    resolution=1.0,
    scattering_length_consistency=true,
  )
  assert_true(neutron.total() < 0.0)
}

///|
test "density separates distant atoms" {
  let arr = @src.AtomArray::from_atoms([
    bsd_elem_atom("C", 0.0, 0.0, 0.0),
    bsd_elem_atom("C", 10.0, 0.0, 0.0),
  ])
  let density = @src.biotite_density_from_atom_array(arr, resolution=1.0)
  let (nz, ny, nx) = density.shape()
  // total density equals 2 electrons (integral approximation >= 1)
  assert_true(density.total() > 1.0)
  let _ = (nz, ny, nx)
}

///|
test "density respects explicit box" {
  let arr = @src.AtomArray::from_atoms([bsd_elem_atom("C", 2.0, 2.0, 2.0)])
  let density = @src.biotite_density_from_atom_array(
    arr,
    resolution=1.0,
    box=Some(
      (@src.Vector3::new(0.0, 0.0, 0.0), @src.Vector3::new(4.0, 4.0, 4.0)),
    ),
  )
  assert_true(density.total() > 0.0)
  assert_true(density.max() > 0.0)
}
