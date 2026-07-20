///|
fn main {
  let motifs = @bio.chromvar_create_example_motifs()
  let data = @bio.chromvar_create_example_data()

  let bias = @bio.BiasCorrection::new()
  let deviations = @bio.chromvar_compute_deviations(data, motifs, bias)

  println(
    "Deviation results: " +
    deviations.motif_names.length().to_string() +
    " motifs",
  )
  println("Number of cells: " + deviations.cell_names.length().to_string())

  let variability = @bio.chromvar_variability(deviations.deviations)
  println(
    "Variability analysis completed: " +
    variability.length().to_string() +
    " results",
  )

  let clusters = @bio.chromvar_cluster_cells(deviations.deviations, "kmeans")
  println(
    "Cell clustering completed: " +
    clusters.length().to_string() +
    " cells assigned",
  )
}
