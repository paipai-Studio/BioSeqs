///|
/// SpatialExperiment demo.

pub fn main() {
  let se = create_example_spatial_experiment()
  println("SpatialExperiment created with")
  println("- rows: " + se_num_rows(se).to_string())
  println("- cols: " + se_num_cols(se).to_string())
  println("- assays: " + se_num_assays(se).to_string())
  
  let (min_x, max_x, min_y, max_y) = se_get_spatial_range(se)
  println("\nSpatial range:")
  println("- X: [" + min_x.to_string() + ", " + max_x.to_string() + "]")
  println("- Y: [" + min_y.to_string() + ", " + max_y.to_string() + "]")
  
  let filtered = se_filter_spots_by_range(se, 150.0, 300.0, 150.0, 250.0)
  println("\nFiltered spots: " + se_num_cols(filtered).to_string())
}