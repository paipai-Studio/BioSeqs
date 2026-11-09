///|
/// Tests for mda_hole.mbt (pore radius profile along an axis).

///|
test "mda_hole_profile_no_atoms_open" {
  // No atoms: perpendicular distance defaults to ~1e10 -> huge radius.
  let profile = @src.mda_hole_profile(
    [],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 4.0],
    n_samples=5,
  )
  assert_eq(profile.length(), 5)
  // Positions are evenly spaced along the 4-long axis.
  assert_true(@moonbit.assert_num_close(profile[0].position(), 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(profile[1].position(), 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(profile[4].position(), 4.0, 1.0e-9))
  assert_true(profile[0].radius() > 1000.0)
}

///|
test "mda_hole_profile_single_wall_atom_constant_perp" {
  // Atom at (3,0,0); axis is the z-line through the origin. The
  // perpendicular distance is 3 at every sample -> radius = 3 - 1 = 2.
  let coords = [[3.0, 0.0, 0.0]]
  let profile = @src.mda_hole_profile(
    coords,
    [0.0, 0.0, -2.0],
    [0.0, 0.0, 2.0],
    n_samples=5,
  )
  assert_eq(profile.length(), 5)
  for p in profile {
    assert_true(@moonbit.assert_num_close(p.radius(), 2.0, 1.0e-9))
  }
  assert_true(@moonbit.assert_num_close(profile[2].position(), 2.0, 1.0e-9))
}

///|
test "mda_hole_profile_radius_clamped_to_zero" {
  // Atom only 0.5 away from the axis: raw radius = 0.5 - 1.0 < 0 -> 0.
  let coords = [[0.5, 0.0, 0.0]]
  let profile = @src.mda_hole_profile(
    coords,
    [0.0, 0.0, -2.0],
    [0.0, 0.0, 2.0],
    n_samples=3,
  )
  assert_eq(profile.length(), 3)
  for p in profile {
    assert_true(@moonbit.assert_num_close(p.radius(), 0.0, 1.0e-12))
  }
}

///|
test "mda_hole_profile_degenerate_axis" {
  // Zero-length axis -> empty profile.
  let profile = @src.mda_hole_profile([[3.0, 0.0, 0.0]], [1.0, 1.0, 1.0], [
    1.0, 1.0, 1.0,
  ])
  assert_eq(profile.length(), 0)
}

///|
test "mda_hole_min_radius_and_bottleneck" {
  // Atoms at perpendicular distances 1 and 2 from the axis; probe 0.5
  // -> every sample has radius 1 - 0.5 = 0.5. The bottleneck is the
  // first sample (position 0) since all radii tie.
  let coords = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
  let profile = @src.mda_hole_profile(
    coords,
    [0.0, 0.0, -2.0],
    [0.0, 0.0, 2.0],
    n_samples=5,
    probe_radius=0.5,
  )
  assert_true(
    @moonbit.assert_num_close(@src.mda_hole_min_radius(profile), 0.5, 1.0e-9),
  )
  let (pos, r) = @src.mda_hole_bottleneck(profile)
  assert_true(@moonbit.assert_num_close(pos, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r, 0.5, 1.0e-9))
}

///|
test "mda_hole_min_radius_empty" {
  assert_true(
    @moonbit.assert_num_close(@src.mda_hole_min_radius([]), 0.0, 1.0e-12),
  )
  let (pos, r) = @src.mda_hole_bottleneck([])
  assert_true(@moonbit.assert_num_close(pos, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r, 0.0, 1.0e-12))
}

///|
test "mda_hole_conductance_hagen_poiseuille" {
  // G = pi * r^4 / (8 * eta * L) = pi / 0.008 for r=1, L=1, eta=0.001.
  let g = @src.mda_hole_conductance(1.0, 1.0)
  assert_true(@moonbit.assert_num_close(g, 3.141592653589793 / 0.008, 1.0e-9))
  // Invalid inputs -> 0.
  assert_true(
    @moonbit.assert_num_close(@src.mda_hole_conductance(0.0, 1.0), 0.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.mda_hole_conductance(1.0, 0.0), 0.0, 1.0e-12),
  )
}

///|
test "mda_hole_classify_thresholds" {
  assert_eq(@src.mda_hole_classify(0.5), "Closed")
  assert_eq(@src.mda_hole_classify(1.0), "Narrow")
  assert_eq(@src.mda_hole_classify(1.5), "Narrow")
  assert_eq(@src.mda_hole_classify(2.0), "Open")
  assert_eq(@src.mda_hole_classify(3.9), "Open")
  assert_eq(@src.mda_hole_classify(4.0), "Wide")
  assert_eq(@src.mda_hole_classify(10.0), "Wide")
}
