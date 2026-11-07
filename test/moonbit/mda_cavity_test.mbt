///|
/// Tests for mda_cavity.mbt (grid-based cavity detection).

///|
test "mda_detect_cavities_empty_input" {
  let cavities = @src.mda_detect_cavities([], [])
  assert_eq(cavities.length(), 0)
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_total_cavity_volume(cavities),
      0.0,
      1.0e-9,
    ),
  )
  match @src.mda_largest_cavity(cavities) {
    Some(_) => fail("expected no cavity for empty input")
    None => ()
  }
}

///|
test "mda_detect_cavities_single_atom_no_cavity" {
  // One atom at the origin, radius 1, no probe. Only the center grid point
  // is occupied; every other point lies on the boundary of the 3x3x3 grid
  // and is reachable from outside, so no enclosed cavity exists.
  let coords = [[0.0, 0.0, 0.0]]
  let radii = [1.0]
  let cavities = @src.mda_detect_cavities(coords, radii, probe_radius=0.0)
  assert_eq(cavities.length(), 0)
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_total_cavity_volume(cavities),
      0.0,
      1.0e-9,
    ),
  )
}

///|
test "mda_detect_cavities_enclosed_void" {
  // 26 atoms of radius 0.45 on the {-1,0,1}^3 lattice (center removed) with
  // probe radius 0.35 -> combined exclusion radius 0.8. The bounding box is
  // [-1.8, 1.8]^3 giving a 4x4x4 grid. Interior grid points at {-0.8,0.2}^3
  // are all occupied (nearest atom 0.2*sqrt(3) ~= 0.35 away) except the one
  // at (0.2,0.2,0.2), whose nearest atoms (0,0,±1)/(0,±1,0)/(±1,0,0) are
  // sqrt(0.72) ~= 0.85 away — it stays free and is fully surrounded by
  // occupied grid points, forming exactly one enclosed cavity.
  let coords : Array[Array[Double]] = Array(capacity=26)
  let radii : Array[Double] = Array(capacity=26)
  for x = -1; x <= 1; x = x + 1 {
    for y = -1; y <= 1; y = y + 1 {
      for z = -1; z <= 1; z = z + 1 {
        if x == 0 && y == 0 && z == 0 {
          continue
        }
        coords.push([x.to_double(), y.to_double(), z.to_double()])
        radii.push(0.45)
      }
    }
  }
  let cavities = @src.mda_detect_cavities(coords, radii, probe_radius=0.35)
  assert_eq(cavities.length(), 1)
  let c = cavities[0]
  assert_eq(c.id, 0)
  assert_eq(c.n_grid_points, 1)
  assert_true(@moonbit.assert_num_close(c.volume, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(c.center[0], 0.2, 1.0e-9))
  assert_true(@moonbit.assert_num_close(c.center[1], 0.2, 1.0e-9))
  assert_true(@moonbit.assert_num_close(c.center[2], 0.2, 1.0e-9))
  assert_eq(c.grid_indices.length(), 1)
  // Aggregated helpers.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_total_cavity_volume(cavities),
      1.0,
      1.0e-9,
    ),
  )
  match @src.mda_largest_cavity(cavities) {
    Some(l) => assert_eq(l.id, 0)
    None => fail("expected a largest cavity")
  }
}

///|
test "mda_detect_cavities_grid_spacing_zero_total_volume" {
  // A single tiny atom with a large probe occupies no grid point when the
  // combined radius is smaller than the grid spacing offset; the bounding
  // box is 1x1x1 grid points, all exterior -> no cavities.
  let coords = [[5.0, 5.0, 5.0]]
  let radii = [0.1]
  let cavities = @src.mda_detect_cavities(
    coords,
    radii,
    grid_spacing=2.0,
    probe_radius=0.0,
  )
  assert_eq(cavities.length(), 0)
}
