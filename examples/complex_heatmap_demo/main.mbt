///|
fn main {
  let data = @bio.ch_create_example_data()
  let options = @bio.ch_create_example_options()

  let result = @bio.ch_draw_heatmap(data, options)
  println(result.text)

  let row_clusters = @bio.ch_cluster_rows(data.matrix, "complete")
  println(
    "Row clustering completed: " + row_clusters.length().to_string() + " rows",
  )

  let scaled_data = @bio.ch_scale_matrix(data.matrix, "row", true)
  println("Scaled matrix has " + scaled_data.length().to_string() + " rows")
}
