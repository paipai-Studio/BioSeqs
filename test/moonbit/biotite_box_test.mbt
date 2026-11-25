///|
/// Tests for biotite_box.mbt — unit cell boxes and periodic boundaries.
test "biotite_box_unit_vectors_and_volume" {
  let box = @src.BiotiteBox::from_vectors([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [
    0.0, 0.0, 1.0,
  ])
  assert_true(@moonbit.assert_num_close(box.volume(), 1.0, 1.0e-9))
  assert_true(box.is_orthogonal())
}

///|
test "biotite_box_from_unitcell_roundtrip" {
  let box = @src.BiotiteBox::from_unitcell(2.0, 3.0, 4.0, 90.0, 90.0, 90.0)
  assert_true(@moonbit.assert_num_close(box.volume(), 24.0, 1.0e-6))
  let (a, b, c, alpha, beta, gamma) = box.to_unitcell()
  assert_true(@moonbit.assert_num_close(a, 2.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(b, 3.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(c, 4.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(alpha, 90.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(beta, 90.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(gamma, 90.0, 1.0e-6))
}

///|
test "biotite_box_fractional_coordinates" {
  let box = @src.BiotiteBox::from_vectors([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [
    0.0, 0.0, 10.0,
  ])
  let frac = box.coord_to_fraction([2.5, 5.0, 10.0])
  assert_true(@moonbit.assert_num_close(frac[0], 0.25, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frac[1], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frac[2], 1.0, 1.0e-9))
  let cart = box.fraction_to_coord(frac)
  assert_true(@moonbit.assert_num_close(cart[0], 2.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(cart[1], 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(cart[2], 10.0, 1.0e-9))
}
