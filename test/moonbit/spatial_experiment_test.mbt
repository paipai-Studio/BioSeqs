///|
/// Tests for SpatialExperiment module.

test "SpatialExperiment creation" {
  let se = @src.create_example_spatial_experiment()
  assert_eq(@src.se_num_rows(se), 3)
  assert_eq(@src.se_num_cols(se), 6)
}

test "SpatialExperiment spatial range" {
  let se = @src.create_example_spatial_experiment()
  let (min_x, max_x, min_y, max_y) = @src.se_get_spatial_range(se)
  assert_eq(min_x, 100.0)
  assert_eq(max_x, 350.0)
  assert_eq(min_y, 150.0)
  assert_eq(max_y, 250.0)
}

test "SpatialExperiment filter spots" {
  let se = @src.create_example_spatial_experiment()
  let filtered = @src.se_filter_spots_by_range(se, 150.0, 300.0, 150.0, 250.0)
  assert_eq(@src.se_num_cols(filtered), 4)
}

test "SpatialCoord new_2d" {
  let coord = @src.SpatialCoord::new_2d(100.0, 200.0)
  assert_eq(coord.x, 100.0)
  assert_eq(coord.y, 200.0)
  assert_eq(coord.z, 0.0)
}