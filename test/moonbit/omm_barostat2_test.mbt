///|
test "omm membrane barostat neutral rng keeps box and coords" {
  let coords = [(1.0, 2.0, 3.0)]
  let (new_coords, new_box) = @src.omm_membrane_barostat(
    coords,
    (10.0, 20.0, 30.0),
    (1.0, 1.0),
    300.0,
    1,
    0.001,
    rng=(0.5, 0.5),
  )
  assert_eq(new_coords.length(), 1)
  assert_true(@moonbit.assert_num_close(new_coords[0].0, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_coords[0].1, 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_coords[0].2, 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.0, 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.1, 20.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.2, 30.0, 1.0e-12))
}

///|
test "omm membrane barostat scales xy jointly and z independently" {
  // scale_xy = 1 + (0.6-0.5)*0.001 = 1.0001, scale_z = 1 + (0.4-0.5)*0.001 = 0.9999
  let coords = [(1.0, 2.0, 3.0)]
  let (new_coords, new_box) = @src.omm_membrane_barostat(
    coords,
    (10.0, 20.0, 30.0),
    (1.0, 1.0),
    300.0,
    1,
    0.001,
    rng=(0.6, 0.4),
  )
  assert_true(@moonbit.assert_num_close(new_coords[0].0, 1.0001, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_coords[0].1, 2.0002, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_coords[0].2, 2.9997, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.0, 10.001, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.1, 20.002, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.2, 29.997, 1.0e-9))
}

///|
test "omm anisotropic barostat scales each axis independently" {
  // scales: x=1.0001, y=1.0, z=0.9999
  let coords = [(1.0, 1.0, 1.0)]
  let (new_coords, new_box) = @src.omm_anisotropic_barostat(
    coords,
    (10.0, 10.0, 10.0),
    (1.0, 1.0, 1.0),
    300.0,
    1,
    0.001,
    rng=(0.6, 0.5, 0.4),
  )
  assert_true(@moonbit.assert_num_close(new_coords[0].0, 1.0001, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_coords[0].1, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_coords[0].2, 0.9999, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.0, 10.001, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.1, 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.2, 9.999, 1.0e-9))
}

///|
test "omm membrane barostat virial zero move accepted with unchanged box" {
  let (_accept_unused, new_box) = @src.omm_membrane_barostat_virial(
    [(0.0, 0.0, 0.0)],
    (10.0, 10.0, 10.0),
    (0.0, 0.0, 0.0),
    (1.0, 1.0),
    0.0,
    300.0,
    1,
    rng=0.5,
  )
  assert_true(@moonbit.assert_num_close(new_box.0, 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.1, 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(new_box.2, 10.0, 1.0e-12))
}

///|
test "omm membrane barostat virial rejects unfavorable compression" {
  // delta_lnV = -0.005, dE = -P*dV*2 = +10 > 0 -> Metropolis rejects at rng=0
  let (accept, new_box) = @src.omm_membrane_barostat_virial(
    [(0.0, 0.0, 0.0)],
    (10.0, 10.0, 10.0),
    (0.0, 0.0, 0.0),
    (1.0, 1.0),
    0.0,
    300.0,
    1,
    rng=0.0,
  )
  assert_false(accept)
  assert_true(
    @moonbit.assert_num_close(new_box.0, 10.0 * @math.exp(-0.0025), 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(new_box.1, 10.0 * @math.exp(-0.0025), 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(new_box.2, 10.0 * @math.exp(-0.005), 1.0e-9),
  )
}

///|
test "omm anisotropic barostat virial rejects x compression" {
  // axis=0, delta=-0.005, dE = -(P - V/L^3)*delta*L^3 = +5 -> reject
  let (accept, new_box) = @src.omm_anisotropic_barostat_virial(
    (10.0, 10.0, 10.0),
    (0.0, 0.0, 0.0),
    (1.0, 1.0, 1.0),
    0.0,
    300.0,
    1,
    axis=0,
    rng=0.0,
  )
  assert_false(accept)
  assert_true(
    @moonbit.assert_num_close(new_box.0, 10.0 * @math.exp(-0.005), 1.0e-9),
  )
  assert_true(@moonbit.assert_num_close(new_box.1, 10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.2, 10.0, 1.0e-9))
}

///|
test "omm anisotropic barostat virial accepts z expansion" {
  // axis=2, delta=+0.004, dE = -(1-0)*0.004*1000 = -4 -> arg>0 -> accept
  let (accept, new_box) = @src.omm_anisotropic_barostat_virial(
    (10.0, 10.0, 10.0),
    (0.0, 0.0, 0.0),
    (1.0, 1.0, 1.0),
    0.0,
    300.0,
    1,
    axis=2,
    rng=0.9,
  )
  assert_true(accept)
  assert_true(@moonbit.assert_num_close(new_box.0, 10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(new_box.1, 10.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(new_box.2, 10.0 * @math.exp(0.004), 1.0e-9),
  )
}
