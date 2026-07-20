///|
fn main {
  let data = @bio.scdf_create_example_data()

  println(
    "Single cell data: " +
    data.cell_names.length().to_string() +
    " cells, " +
    data.gene_names.length().to_string() +
    " genes",
  )

  let scores = @bio.scdf_compute_doublet_score(data, 10)
  println(
    "Doublet scores computed for " + scores.length().to_string() + " cells",
  )

  let doublets = @bio.scdf_detect_doublets(scores, 0.3)
  println("Detected doublets: " + doublets.length().to_string())

  let summary = @bio.scdf_doublet_summary(scores)
  println(summary)

  let filtered = @bio.scdf_filter_doublets(data, scores)
  println(
    "After filtering: " + filtered.cell_names.length().to_string() + " cells",
  )

  let pca = @bio.scdf_compute_pca(data, 5)
  println("PCA computed: " + pca.length().to_string() + " cells, 5 components")
}
