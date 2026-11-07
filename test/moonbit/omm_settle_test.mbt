///|
test "omm_settle_ideal_water_unchanged" {
  let d_oh = 0.9572
  let d_hh = 1.5139
  let half = @math.asin(d_hh / (2.0 * d_oh))
  let hx = d_oh * @math.sin(half)
  let hz = d_oh * @math.cos(half)
  let coords = [(0.0, 0.0, 0.0), (hx, 0.0, hz), (-hx, 0.0, hz)]
  let zero = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  let masses = [15.9994, 1.008, 1.008]
  let waters = [(0, 1, 2)]
  let (nc, nv) = @src.omm_settle(coords, zero, masses, waters, 0.001)
  for i = 0; i < 3; i = i + 1 {
    let (x0, y0, z0) = coords[i]
    let (x1, y1, z1) = nc[i]
    assert_true(@moonbit.assert_num_close(x1, x0, 1.0e-9))
    assert_true(@moonbit.assert_num_close(y1, y0, 1.0e-9))
    assert_true(@moonbit.assert_num_close(z1, z0, 1.0e-9))
    let (vx, vy, vz) = nv[i]
    assert_true(@moonbit.assert_num_close(vx, 0.0, 1.0e-9))
    assert_true(@moonbit.assert_num_close(vy, 0.0, 1.0e-9))
    assert_true(@moonbit.assert_num_close(vz, 0.0, 1.0e-9))
  }
}

///|
test "omm_settle_restores_stretched_water" {
  // O-H2 stretched to 1.5 while O-H1 is ideal and non-collinear
  let coords = [(0.0, 0.0, 0.0), (0.9572, 0.0, 0.0), (0.0, 1.5, 0.0)]
  let zero = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  let masses = [15.9994, 1.008, 1.008]
  let waters = [(0, 1, 2)]
  let (nc, _nv) = @src.omm_settle(coords, zero, masses, waters, 0.001)
  let (ox, oy, oz) = nc[0]
  let (h1x, h1y, h1z) = nc[1]
  let (h2x, h2y, h2z) = nc[2]
  // O stays at origin with zero velocities
  assert_true(@moonbit.assert_num_close(ox, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(oy, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(oz, 0.0, 1.0e-9))
  // O-H1 constraint restored exactly
  let d_oh = @math.pow(
    (h1x - ox) * (h1x - ox) + (h1y - oy) * (h1y - oy) + (h1z - oz) * (h1z - oz),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(d_oh, 0.9572, 1.0e-9))
  // H1-H2 constraint restored exactly
  let d_hh = @math.pow(
    (h2x - h1x) * (h2x - h1x) +
    (h2y - h1y) * (h2y - h1y) +
    (h2z - h1z) * (h2z - h1z),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(d_hh, 1.5139, 1.0e-9))
}

///|
test "omm_settle_skips_out_of_range_molecule" {
  let coords = [(0.0, 0.0, 0.0), (0.9572, 0.0, 0.0)]
  let zero = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  let masses = [16.0, 1.0]
  // indices beyond the coordinate array must be ignored
  let waters = [(0, 1, 99)]
  let (nc, _nv) = @src.omm_settle(coords, zero, masses, waters, 0.001)
  let (x1, _, _) = nc[1]
  assert_true(@moonbit.assert_num_close(x1, 0.9572, 1.0e-12))
}

///|
test "omm_ccma_converges_to_target_distance" {
  let coords = [(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)]
  let zero = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  let masses = [1.0, 1.0]
  let constraints = [(0, 1, 1.0)]
  let (nc, nv) = @src.omm_ccma_constraints(
    coords, zero, masses, constraints, 0.001,
  )
  let dx = nc[1].0 - nc[0].0
  assert_true(@moonbit.assert_num_close(dx, 1.0, 1.0e-9))
  // equal masses: center of mass preserved
  assert_true(@moonbit.assert_num_close(nc[0].0 + nc[1].0, 1.2, 1.0e-9))
  // velocities consistent with displacement
  let (vx0, _, _) = nv[0]
  assert_true(@moonbit.assert_num_close(vx0, 100.0, 1.0e-6))
  let (vx1, _, _) = nv[1]
  assert_true(@moonbit.assert_num_close(vx1, -100.0, 1.0e-6))
}

///|
test "omm_ccma_no_constraints_passthrough" {
  let coords = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
  let vels = [(0.1, 0.2, 0.3), (-0.1, 0.0, 0.5)]
  let masses = [1.0, 2.0]
  let (nc, nv) = @src.omm_ccma_constraints(coords, vels, masses, [], 0.01)
  // without constraints the positions are a pure Euler drift
  let (x, y, z) = nc[1]
  assert_true(@moonbit.assert_num_close(x, 1.0 - 0.1 * 0.01, 1.0e-12))
  assert_true(@moonbit.assert_num_close(y, 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(z, 3.0 + 0.5 * 0.01, 1.0e-12))
  // velocities are reconstructed as (new - old) / dt
  let (vx, _, _) = nv[0]
  assert_true(@moonbit.assert_num_close(vx, 0.1, 1.0e-12))
  let (_, _, vz1) = nv[1]
  assert_true(@moonbit.assert_num_close(vz1, 0.5, 1.0e-12))
}
